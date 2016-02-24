##' Filter GOs list
##'
##' Filter GOs list and append GO names when appropriate
##'
##' @param go.env GO or gene set list
##' @param min.size Minimum size for number of genes in a gene set (default: 5)
##' @param max.size Maximum size for number of genes in a gene set (default: 5000)
##' @param annot Whether to append GO annotations for easier interpretation (default: FALSE)
##'
##' @return a filtered GO list
##'
##' @examples
##' \donttest{
##' # 10 sample GOs
##' library(org.Hs.eg.db)
##' go.env <- mget(ls(org.Hs.egGO2ALLEGS)[1:10], org.Hs.egGO2ALLEGS)
##' # Filter this list and append names for easier interpretation
##' go.env <- clean.gos(go.env)
##' }
##'
##' @export
clean.gos <- function(go.env, min.size = 5, max.size = 5000, annot = FALSE) {
  go.env <- as.list(go.env)
  size <- unlist(lapply(go.env, length))
  go.env <- go.env[size > min.size & size < max.size]
  # If we have GO.db installed, then add the term to each GO code.
  if (annot && "GO.db" %in% installed.packages()[,1]) {
    desc <- select(
      GO.db,
      keys = names(go.env),
      columns = c("TERM"),
      multiVals = 'CharacterList'
    )
    stopifnot(all(names(go.env) == desc$GOID))
    names(go.env) <- paste(names(go.env), desc$TERM)
  }
  return(go.env)
}

##' Filter counts matrix
##'
##' Filter counts matrix based on gene and cell requirements
##'
##' @param counts read count matrix. The rows correspond to genes, columns correspond to individual cells
##' @param min.lib.size Minimum number of genes detected in a cell. Cells with fewer genes will be removed (default: 1.8e3)
##' @param min.reads Minimum number of reads per gene. Genes with fewer reads will be removed (default: 10)
##' @param min.detected Minimum number of cells a gene must be seen in. Genes not seen in a sufficient number of cells will be removed (default: 5)
##'
##' @return a filtered read count matrix
##'
##' @examples
##' data(pollen)
##' dim(pollen)
##' cd <- clean.counts(pollen)
##' dim(cd)
##'
##' @export
clean.counts <- function(counts, min.lib.size = 1.8e3, min.reads = 10, min.detected = 5) {
    # filter out low-gene cells
    counts <- counts[, colSums(counts>0)>min.lib.size]
    # remove genes that don't have many reads
    counts <- counts[rowSums(counts)>min.reads, ]
    # remove genes that are not seen in a sufficient number of cells
    counts <- counts[rowSums(counts>0)>min.detected, ]
    return(counts)
}

sn <- function(x) {
    names(x) <- x
    return(x)
}

# convert R color to a web hex representation
col2hex <- function(col) {
    unlist(lapply(col, function(c) {
        c <- col2rgb(c)
        sprintf("#%02X%02X%02X", c[1], c[2], c[3])
    }))
}

##' wrapper around different mclapply mechanisms
##'
##' Abstracts out mclapply implementation, and defaults to lapply when only one core is requested (helps with debugging)
##' @param ... parameters to pass to lapply, mclapply, bplapply, etc.
##' @param n.cores number of cores. If 1 core is requested, will default to lapply
papply <- function(...,n.cores=detectCores()) {
  if(n.cores>1) {
    # bplapply implementation
    if(is.element("parallel", installed.packages()[,1])) {
      mclapply(...,mc.cores=n.cores)
    } else {
      # last resort
      bplapply(... , BPPARAM = MulticoreParam(workers = n.cores))
    }
  } else { # fall back on lapply
    lapply(...);
  }
}


