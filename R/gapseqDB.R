#' S4 class for gapseq's biochemistry database
#'
#' @slot rxn A list of tables (\link(data.table)) with reaction information
#' @slot cpd A list of tables (\link(data.table)) with compound information
#'
#' @export
setClass("gapseqDB",

         slots = c(
           rxn = "list",
           cpd = "list"
         ),

         prototype = list(
           rxn = list(),
           cpd = list()
         )
)

#' Initialize a gapseq biochemistry database object
#'
#' @description Biochemistry data is retrieved either from gapseq's main github
#' repository or from a local directory
#'
#' @param src Location of gagpseq (see details)
#'
#' @details
#' If src is NULL, the latest database files are retrieved from
#' the latest commit on the master branch of gapseq from the github repository
#' (https://github.com/jotech/gapseq).
#'
#' If 'src' is a git commit SHA fingerprint, the database files are retrieved
#' from the master branch from github for the specified commit.
#'
#' If 'src' is a path to a local gapseq installation directory, database
#' files are directly retrieved from there.
#'
#' @import stringr
#'
#' @export
initDB <- function(src = NULL) {

  # Case A: Retrieve db from latest commit on master branch from github
  if(is.null(src)) {
    pfx <- "https://raw.githubusercontent.com/jotech/gapseq/master/"
  }

  # Case B: Retrieve db from specific commit on master branch from github
  if(!is.null(src) && !dir.exists(src)) {
    pfx <- paste0("https://raw.githubusercontent.com/jotech/gapseq/",src,"/")
  }

  # Case C: Retrieve db from local gapseq installation directory
  if(!is.null(src) && dir.exists(src)) {
    pfx <- src
    if(!grepl("/$",pfx))
      pfx <- paste0(pfx,"/")
  }

  #–––––––#
  # links #
  #–––––––#

  # mets
  link.mets <- paste0(pfx,"dat/seed_metabolites_edited.tsv")

  # reactions
  link.rxns <- paste0(pfx,"dat/seed_reactions_corrected.tsv")
  link.ec   <- paste0(pfx,"dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv")

  #–––––––––––––#
  # Metabolites #
  #–––––––––––––#
  dl.mets <- fread(link.mets)

  mets <- dl.mets[, .(id, abbreviation, name, formula, mass, source, charge,
                      smiles, InChI, InChIKey)]
  mets[formula == "null", formula := NA_character_]
  setkey(mets, "id")

  xrefdbs <- c("MNX_ID", "hmdbID","reactomeID","chebiID","keggID","biggID",
               "biocycID")

  mets.xref <- dl.mets[, c("id", xrefdbs), with = F]
  mets.xref <- melt(mets.xref, id.vars = "id", measure.vars = xrefdbs,
                    variable.name = "resource", value.name = "resourceID")
  mets.xref[resourceID == "", resourceID := NA_character_]
  setkey(mets.xref, id)
  mets.xref <- mets.xref[, list(resourceID = unlist(strsplit(resourceID, ";"))),
                         by=.(id, resource)]

  #–––––––––––––#
  # Reactions   #
  #–––––––––––––#
  dl.rxns <- fread(link.rxns)

  rxns <- dl.rxns[, .(id, abbreviation, name, stoichiometry, reversibility,
                      notes, is_copy_of, gapseq.status)]
  setkey(rxns, "id")
  rxns.stoich <- rxns[, list(stoichiometry = unlist(strsplit(stoichiometry, ";"))),
                      by=.(id)]
  rxns.stoich[, stoichiometry := str_match(stoichiometry, ".+\\:cpd[0-9]{5}\\:[0-9]")[,1]]
  rxns.stoich[, compound := str_match(stoichiometry, "cpd[0-9]{5}")]
  rxns.stoich[, compartment := as.numeric(str_match(stoichiometry, "[0-9]$"))]
  rxns.stoich[, stoichiometry := as.numeric(str_match(stoichiometry, "^[0-9]+|^-[0-9]+"))]
  rxns.stoich <- merge(rxns.stoich, rxns[,.(id, reversibility)], by = "id")
  setkey(rxns.stoich, "id")

  dl.ecs <- fread(link.ec)
  dl.ecs <- dl.ecs[, list(id = unlist(strsplit(`MS ID`, "\\|"))),
                   by=.(`External ID`, Source)]
  rxns.ec <- dl.ecs[, .(id, EC = `External ID`, source = Source)]
  setkey(rxns.ec, "id")

  #––––––––––––––––––––––––#
  # construct S4 DB object #
  #––––––––––––––––––––––––#
  out <- new("gapseqDB",
             rxn = list(main   = rxns,
                        ec     = rxns.ec,
                        stoich = rxns.stoich),
             cpd = list(main   = mets,
                        xref   = mets.xref))


  return(out)
}

#' @import crayon
colorizeRxnID <- function(db, id) {
  gs.corrected <- make_style("palegreen3", bg = TRUE)
  gs.approved <- make_style("palegreen1", bg = TRUE)
  gs.removed <- make_style("pink", bg = TRUE)
  gs.notassessed <- make_style("khaki", bg = TRUE)

  # corrected
  ind <- which(id %in% db@rxn$main[gapseq.status == "corrected", id])
  if(length(ind)>0)
    id[ind] <- gs.corrected(id[ind])

  # corrected
  ind <- which(id %in% db@rxn$main[gapseq.status == "approved", id])
  if(length(ind)>0)
    id[ind] <- gs.approved(id[ind])

  # removed
  ind <- which(id %in% db@rxn$main[gapseq.status == "removed", id])
  if(length(ind)>0)
    id[ind] <- gs.removed(id[ind])

  # removed
  ind <- which(id %in% db@rxn$main[gapseq.status == "not.assessed", id])
  if(length(ind)>0)
    id[ind] <- gs.notassessed(id[ind])

  return(id)
}
