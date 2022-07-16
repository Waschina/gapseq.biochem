#' Get reaction equations
#'
#' @description Nicely formatted reaction equations with compound names or IDs.
#'
#' @param db Database of class \link(gapseqDB)
#' @param rxn Reaction ID(s)
#' @param highlight Compound IDs whose names/IDs should be highlighted in red.
#' (only works if 'format.style' is TRUE)
#' @param use.ids Use compound IDs instead of compound names?
#' @param format.styles Should the output be formatted (highlights and grey
#' background for external metabolites)?
#'
#' @import crayon
#'
#' @export
getRxnEquation <- function(db, rxn, highlight = NULL, use.ids = FALSE,
                           format.style = FALSE) {
  if(length(rxn) == 0)
    return(character(0L))

  # sanity checks
  rxn <- sanityCheckRxns(db, rxn)

  res <- rep(NA_character_, length(rxn))
  names(res) <- rxn

  for(rxni in rxn) {
    # get reversibility
    revi <- db@rxn$main[rxni, reversibility]
    revi <- ifelse(revi == "=","<=>",ifelse(revi == ">","––>","<––"))

    # get stoichiometries
    dttmp <- db@rxn$stoich[rxni][order(stoichiometry)]
    cpd_str <- db@cpd$main[dttmp$compound, name]
    cpd_ids   <- dttmp$compound
    cpd_comp  <- ifelse(dttmp$compartment == 0, "c",
                        ifelse(dttmp$compartment == 1, "e", "p"))
    if(use.ids)
      cpd_str <- cpd_ids

    cpd_str   <- paste0("(",bold(abs(dttmp$stoichiometry)),") ",
                        cpd_str,
                        "[",cpd_comp,"]")

    # highlighting anything?
    ind_hl <- which(cpd_ids %in% highlight)
    if(length(ind_hl) > 0)
      cpd_str[ind_hl] <- red(cpd_str[ind_hl])

    # any external metabolites
    ind_ex <- which(dttmp$compartment != 0)
    if(length(ind_ex) > 0)
      cpd_str[ind_ex] <- bgWhite(cpd_str[ind_ex])

    # combine equation
    res[rxni] <- paste(paste(cpd_str[dttmp$stoichiometry < 0], collapse = " + "),
                       bold(revi),
                       paste(cpd_str[dttmp$stoichiometry > 0], collapse = " + "))


  }

  if(format.style == FALSE)
    res <- strip_style(res)

  return(res)
}

#' Get all reactions with a specific compound
#'
#' @description Find all reaction where a specific compound participates grouped
#' by consuming, producing, consuming+producing, transport.
#'
#' @param db Database of class \link(gapseqDB)
#' @param cpd Compound ID
#'
#' @returns A list with four character vectors (consuming, producing,
#' consuming+producing, transport)
#'
#' @import crayon
#'
#' @export
cpdReactions <- function(db, cpd) {
  if(length(cpd) > 1)
    stop("More than one compound ID provided, which is not supported (yet)")

  if(!all(cpd %in% db@cpd$main$id))
    stop(paste0("Compound '", cpd,"' not part of DB."))

  dttmp <- db@rxn$stoich[compound == cpd]

  # transports
  ids.trans <- dttmp[, .N, by = id][N > 1, id]

  # producing
  ids.prod   <- dttmp[(stoichiometry > 0 & reversibility == ">") |
                        (stoichiometry < 0 & reversibility == "<"), id]
  ids.prod <- ids.prod[!(ids.prod %in% ids.trans)]

  # consuming
  ids.cons   <- dttmp[(stoichiometry > 0 & reversibility == "<") |
                        (stoichiometry < 0 & reversibility == ">"), id]
  ids.cons <- ids.cons[!(ids.cons %in% ids.trans)]

  # consuming + producing
  ids.prco   <- dttmp[reversibility == "=", id]
  ids.prco <- ids.prco[!(ids.prco %in% ids.trans)]

  return(list(consuming = ids.cons,
              producing = ids.prod,
              `consuming+producing` = ids.prco,
              transporters = ids.trans))
}

#' Calculate charge balance
#'
#' @description Calculates the charge balance for reaction(s)
#'
#' @param db Database of class \link(gapseqDB)
#' @param rxn Reaction ID(s)
#'
#' @export
getChargeBalance <- function(db, rxn) {
  if(length(rxn) == 0)
    return(character(0L))

  # sanity checks
  rxn <- sanityCheckRxns(db, rxn)

  out <- rep(NA_real_, length(rxn))
  names(out) <- rxn
  for(rxni in rxn) {
    dttmp <- db@rxn$stoich[rxni]
    chargetmp <- db@cpd$main[dttmp$compound, charge]
    out[rxni] <- sum(dttmp$stoichiometry * chargetmp)
  }

  return(out)
}

#' Calculate mass balance
#'
#' @description Calculates the mass balance for reaction(s)
#'
#' @param db Database of class \link(gapseqDB)
#' @param rxn Reaction ID(s)
#'
#' @importFrom CHNOSZ makeup
#'
#' @export
getMassBalance <- function(db, rxn) {
  if(length(rxn) == 0)
    return(character(0L))

  # sanity checks
  rxn <- sanityCheckRxns(db, rxn)

  out <- rep(NA_character_, length(rxn))
  names(out) <- rxn
  for(rxni in rxn) {
    dttmp <- db@rxn$stoich[rxni]
    mftmp <- db@cpd$main[dttmp$compound, formula]

    if(any(is.na(mftmp))) {
      out[rxni] <- "unknown"
    } else {
      suppressWarnings(makeup.tmp <- makeup(mftmp, multiplier = dttmp$stoichiometry,
                                            count.zero = TRUE))

      tmp_mbal <- Reduce(`+`, makeup.tmp)

      # wild guess elements
      wge <- c("R","X","Y","Z")

      if(any(names(tmp_mbal) %in% wge))
        warning("At least one compound formula contains 'R','X','Y','Z'. Mass balance might be faulty.")

      tmp_mbal <- tmp_mbal[!(names(tmp_mbal) %in% wge)]

      if(all(tmp_mbal == 0)) {
        out[rxni] <- "Ok"
      } else {
        tmp_mbal <- tmp_mbal[tmp_mbal != 0]
        out[rxni] <- paste(paste(names(tmp_mbal), tmp_mbal, sep = ":"),
                           collapse = ", ")
      }
    }
  }

  return(out)
}

#' Find reaction(s)
#'
#' @description Find a reaction by single or combined criteria.
#'
#' @param db Database of class \link(gapseqDB)
#' @param cpd Compound IDs to look for.
#' @param cpd.link "AND" or "OR". "AND" requires that all IDs in 'cpd.link'
#' occur in a reaction. "OR" requires only one compound to be present.
#' @param ec EC number to look for.
#' @param grep.str Regular expression applied to the reaction equation.
#'
#' @export
findReactions <- function(db,
                          cpd = NULL, cpd.link = "AND",
                          ec = NULL,
                          grep.str = NULL,
                          global.link = "AND") {
  highlights <- cpd

  #––––––––––––––#
  # compound ids #
  #––––––––––––––#
  cpd.hits <- NA
  if(!is.null(cpd)) {
    dttmp <- db@rxn$stoich[compound %in% cpd][!duplicated(paste0(id, compound))]
    if(cpd.link == "AND") {
      cpd.hits <- dttmp[,.N, by = id][N == length(cpd), id]
    }
    if(cpd.link == "OR") {
      cpd.hits <- dttmp[, unique(id)]
    }
  }

  #––––––––––––––#
  # compound ids #
  #––––––––––––––#
  grep.hits <- NA
  if(!is.null(grep.str)) {

  }

  comb.hits <- cpd.hits
  comb.hits <- sort(comb.hits)

  if(length(comb.hits) > 0) {
    rxn.hits <- colorizeRxnID(db, comb.hits)
    hits.equation <- getRxnEquation(db, comb.hits,
                                    highlight = highlights,
                                    format.style = TRUE)
    cat(paste0(rxn.hits,": ", hits.equation), sep = "\n")
  }

}

sanityCheckRxns <- function(db, rxn) {
  tmp <- rxn[rxn %in% db@rxn$main$id]
  if(length(tmp) == 0)
    stop("Reaction ID(s) not in database")
  if(length(tmp) < length(rxn)) {
    tmp.abs <- rxn[!(rxn %in% tmp)]
    warning(paste0("Some reactions IDs are not in database:\n  ",
                   paste0(tmp.abs, collapse = "\n  ")))
  }
  return(tmp)
}
