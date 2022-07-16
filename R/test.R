# library(stringr)
#
#
# # mets
# url.mets <- "https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites_edited.tsv"
#
# # reactions
# url.rxns <- "https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions_corrected.tsv"
# url.ec   <- "https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"
#
#
# #–––––––––––––#
# # Metabolites #
# #–––––––––––––#
# dl.mets <- fread(url.mets)
#
# mets <- dl.mets[, .(id, abbreviation, name, formula, mass, source, charge,
#                     smiles, InChI, InChIKey)]
#
# xrefdbs <- c("MNX_ID", "hmdbID","reactomeID","chebiID","keggID","biggID",
#              "biocycID")
#
# mets.xref <- dl.mets[, c("id", xrefdbs), with = F]
# mets.xref <- melt(mets.xref, id.vars = "id", measure.vars = xrefdbs,
#                   variable.name = "resource", value.name = "resourceID")
# mets.xref[resourceID == "", resourceID := NA_character_]
# setkey(mets.xref, id)
# mets.xref <- mets.xref[, list(resourceID = unlist(strsplit(resourceID, ";"))),
#                        by=.(id, resource)]
# #mets.xref[id == "cpd19507"]
#
# #–––––––––––––#
# # Reactions   #
# #–––––––––––––#
# dl.rxns <- fread(url.rxns)
#
# rxns <- dl.rxns[, .(id, abbreviation, name, stoichiometry, reversibility,
#                     notes, is_copy_of, gapseq.status)]
# setkey(rxns, "id")
# rxns.stoich <- rxns[, list(stoichiometry = unlist(strsplit(stoichiometry, ";"))),
#                     by=.(id)]
# rxns.stoich[, stoichiometry := str_match(stoichiometry, ".+\\:cpd[0-9]{5}\\:[0-9]")[,1]]
# rxns.stoich[, compound := str_match(stoichiometry, "cpd[0-9]{5}")]
# rxns.stoich[, compartment := str_match(stoichiometry, "[0-9]$")]
# rxns.stoich[, stoichiometry := str_match(stoichiometry, "^[0-9]+|^-[0-9]+")]
#
# dl.ecs <- fread(url.ec)
#
# dl.ecs <- dl.ecs[, list(id = unlist(strsplit(`MS ID`, "\\|"))),
#                  by=.(`External ID`, Source)]
