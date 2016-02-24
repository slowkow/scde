##' Single-cell Differential Expression (with Pathway And Gene set Overdispersion Analysis)
##'
##' The scde package implements a set of statistical methods for analyzing single-cell RNA-seq data.
##' scde fits individual error models for single-cell RNA-seq measurements. These models can then be used for
##' assessment of differential expression between groups of cells, as well as other types of analysis.
##' The scde package also contains the pagoda framework which applies pathway and gene set overdispersion analysis
##' to identify and characterize putative cell subpopulations based on transcriptional signatures.
##' See vignette("diffexp") for a brief tutorial on differential expression analysis.
##' See vignette("pagoda") for a brief tutorial on pathway and gene set overdispersion analysis to identify and characterize cell subpopulations.
##' More extensive tutorials are available at \url{http://pklab.med.harvard.edu/scde/index.html}.
##'  (test)
##' @name scde
##' @docType package
##' @author Peter Kharchenko \email{Peter_Kharchenko@@hms.harvard.edu}
##' @author Jean Fan \email{jeanfan@@fas.harvard.edu}
NULL

# clean up stale web server reference
.onAttach <- function(...) {

    if(exists("___scde.server", envir = globalenv())) {
        old.server <- get("___scde.server", envir = globalenv())
        n.apps <- length(old.server$appList)-1
        # TODO fix server rescue...
        packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. removing.")
        # remove
        rm("___scde.server", envir = globalenv())
        return(TRUE)

        if(n.apps > 0) {
            require(Rook)
            require(rjson)
            packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. restarting.")
            rm("___scde.server", envir = globalenv()) # remove old instance (apparently saved Rook servers can't just be restarted ... we'll make a new one and re-add all of the apps

            tryCatch( {
                server <- get.scde.server(ip = old.server$listenAddr, port = old.server$listenPort) # launch a new server
                if(!is.null(server)) {
                    lapply(old.server$appList[-1], function(sa) {
                        server$add(app = sa$app, name = sa$name)
                    })
                }
            }, error = function(e) message(e))

        } else {
            packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. removing.")
            # remove
            rm("___scde.server", envir = globalenv())
        }
    }
}

.onUnload <- function(libpath) {
    library.dynam.unload("scde", libpath, verbose = TRUE)
}

