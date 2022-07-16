sanityCheckCpds <- function(db, cpd) {
  tmp <- cpd[cpd %in% db@cpd$main$id]
  if(length(tmp) == 0)
    stop("Compound ID(s) not in database")
  if(length(tmp) < length(cpd)) {
    tmp.abs <- cpd[!(cpd %in% tmp)]
    warning(paste0("Some compound IDs are not in database:\n  ",
                   paste0(tmp.abs, collapse = "\n  ")))
  }
  return(tmp)
}
