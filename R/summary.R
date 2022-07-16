#' Summary for a compound/reaction from gapseq biochemistry database
#'
#' @description Prints a summary for a reaction or compound to console prompt
#'
#' @param object Database of class \link(gapseqDB)
#' @param id ID of a compound or reaction
#'
#' @details Compound IDs should have the prefix "cpd" and reactions "rxn"
#'
#' @import stringr
#' @import crayon
#'
#' @export
setMethod(f          = "summary",
          signature  = signature(object = "gapseqDB"),
          definition = function(object, id) {

            db <- object

            mode <- str_match(id, "^cpd|^rxn")
            if(mode != "rxn" & mode != "cpd")
              stop("Invalid ID. ID should have the prefix 'rxn' or 'cpd'.")

            if(mode == "cpd")
              summaryCPD(object, id)

            if(mode == "rxn")
              summaryRXN(object, id)

          }
)

#' @import crayon
summaryCPD <- function(db, id) {
  #––––––––––––#
  # Metabolite #
  #------------#
  id <- sanityCheckCpds(db, id)
  cat(bold(id),"–",db@cpd$main[id, name],"\n\n")

  cat("Mol. formula: \t",db@cpd$main[id,formula],"\n")
  cat("Charge: \t",db@cpd$main[id,charge],"\n")
  cat("Abbreviation: \t",db@cpd$main[id,abbreviation],"\n")
  cat("\n")

  cpdrxn <- cpdReactions(db, id)

  # consuming reactions
  c.rxns <- getRxnEquation(db, cpdrxn$consuming, highlight = id,
                           format.style = TRUE)
  c.rxns.ids <- names(c.rxns)
  c.rxns.ids <- colorizeRxnID(db, c.rxns.ids)

  if(length(c.rxns) > 0) {
    c.rxns <- paste0(c.rxns.ids, ": ", c.rxns)
    cat(bold("Reactions consuming " %+% id %+% ":\n"))
    cat(c.rxns, sep = "\n")
    cat("\n")
  }

  # producing reactions
  c.rxns <- getRxnEquation(db, cpdrxn$producing, highlight = id,
                           format.style = TRUE)
  c.rxns.ids <- names(c.rxns)
  c.rxns.ids <- colorizeRxnID(db, c.rxns.ids)

  if(length(c.rxns) > 0) {
    c.rxns <- paste0(c.rxns.ids, ": ", c.rxns)
    cat(bold("Reactions producing " %+% id %+% ":\n"))
    cat(c.rxns, sep = "\n")
    cat("\n")
  }

  # consuming+producing reactions
  c.rxns <- getRxnEquation(db, cpdrxn$`consuming+producing`, highlight = id,
                           format.style = TRUE)
  c.rxns.ids <- names(c.rxns)
  c.rxns.ids <- colorizeRxnID(db, c.rxns.ids)

  if(length(c.rxns) > 0) {
    c.rxns <- paste0(c.rxns.ids, ": ", c.rxns)
    cat(bold("Reactions consuming+producing " %+% id %+% ":\n"))
    cat(c.rxns, sep = "\n")
    cat("\n")
  }


  # transporters
  c.rxns <- getRxnEquation(db, cpdrxn$`transporter`, highlight = id,
                           format.style = TRUE)
  c.rxns.ids <- names(c.rxns)
  c.rxns.ids <- colorizeRxnID(db, c.rxns.ids)

  if(length(c.rxns) > 0) {
    c.rxns <- paste0(c.rxns.ids, ": ", c.rxns)
    cat(bold("Transporters of " %+% id %+% ":\n"))
    cat(c.rxns, sep = "\n")
    cat("\n")
  }
}

#' @import crayon
summaryRXN <- function(db, id) {

  #––––––––––––#
  # Reaction   #
  #------------#
  id <- sanityCheckRxns(db, id)
  cat(bold(colorizeRxnID(db, id)),"–",db@rxn$main[id, name],"\n\n")

  cat("Abbreviation:\t", db@rxn$main[id, abbreviation],"\n")
  cat("gapseq status:\t", db@rxn$main[id, gapseq.status],"\n")
  cat("EC:\t\t", paste(db@rxn$ec[id, EC], collapse = "; "),"\n")

  # chemical balances
  cb <- getChargeBalance(db, id)
  cb <- ifelse(cb == 0, green(cb), red(bold(cb)))
  cat("Charge balance:\t", cb,"\n")
  mb <- getMassBalance(db, id)
  mb <- ifelse(mb == "Ok", green(mb), ifelse(mb == "unknown",yellow(mb),
                                             red(bold(mb))))
  cat("Mass balance:\t", mb,"\n")
  cat("\n")

  cat("Equation:\n",getRxnEquation(db, id, format.style = TRUE),"\n\n")
  cat("Equation (IDs):\n",getRxnEquation(db, id, use.ids = TRUE,
                                         format.style = TRUE),"\n\n")

  # Participating compound infos: id, sum.formula, charge
  dttmp <- copy(db@rxn$stoich[id])
  dttmp$charge <- db@cpd$main[dttmp$compound, charge]
  dttmp$metsf <- db@cpd$main[dttmp$compound, formula]
  dttmp[, cpd.str := paste0(compound," (",ifelse(charge>0,"+",""),charge," * ",abs(stoichiometry)," = ",charge*abs(stoichiometry),")\t", metsf)]

  cat("LHS:\n")
  cat(paste(dttmp[stoichiometry < 0, cpd.str], collapse = "\n"),"\n\n")

  cat("RHS:\n")
  cat(paste(dttmp[stoichiometry > 0, cpd.str], collapse = "\n"),"\n\n")
}
