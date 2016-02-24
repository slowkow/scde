################################# SCDE Methods

##' Fit single-cell error/regression models
##'
##' Fit error models given a set of single-cell data (counts) and an optional grouping factor (groups). The cells (within each group) are first cross-compared to determine a subset of genes showing consistent expression. The set of genes is then used to fit a mixture model (Poisson-NB mixture, with expression-dependent concomitant).
##'
##' Note: the default implementation has been changed to use linear-scale fit with expression-dependent NB size (overdispersion) fit. This represents an interative improvement on the originally published model. Use linear.fit=F to revert back to the original fitting procedure.
##'
##' @param counts read count matrix. The rows correspond to genes (should be named), columns correspond to individual cells. The matrix should contain integer counts
##' @param groups an optional factor describing grouping of different cells. If provided, the cross-fits and the expected expression magnitudes will be determined separately within each group. The factor should have the same length as ncol(counts).
##' @param min.nonfailed minimal number of non-failed observations required for a gene to be used in the final model fitting
##' @param threshold.segmentation use a fast threshold-based segmentation during cross-fit (default: TRUE)
##' @param min.count.threshold the number of reads to use to guess which genes may have "failed" to be detected in a given measurement during cross-cell comparison (default: 4)
##' @param zero.count.threshold threshold to guess the initial value (failed/non-failed) during error model fitting procedure (defaults to the min.count.threshold value)
##' @param zero.lambda the rate of the Poisson (failure) component (default: 0.1)
##' @param save.crossfit.plots whether png files showing cross-fit segmentations should be written out (default: FALSE)
##' @param save.model.plots whether pdf files showing model fits should be written out (default = TRUE)
##' @param n.cores number of cores to use
##' @param min.size.entries minimum number of genes to use when determining expected expression magnitude during model fitting
##' @param max.pairs maximum number of cross-fit comparisons that should be performed per group (default: 5000)
##' @param min.pairs.per.cell minimum number of pairs that each cell should be cross-compared with
##' @param verbose 1 for increased output
##' @param linear.fit Boolean of whether to use a linear fit in the regression (default: TRUE).
##' @param local.theta.fit Boolean of whether to fit the overdispersion parameter theta, ie. the negative binomial size parameter, based on local regression (default: set to be equal to the linear.fit parameter)
##' @param theta.fit.range Range of valid values for the overdispersion parameter theta, ie. the negative binomial size parameter (default: c(1e-2, 1e2))
##'
##' @return a model matrix, with rows corresponding to different cells, and columns representing different parameters of the determined models
##'
##' @useDynLib scde
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' \donttest{
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' }
##'
##' @export
scde.error.models <- function(counts, groups = NULL, min.nonfailed = 3, threshold.segmentation = TRUE, min.count.threshold = 4, zero.count.threshold = min.count.threshold, zero.lambda = 0.1, save.crossfit.plots = FALSE, save.model.plots = TRUE, n.cores = 12, min.size.entries = 2e3, max.pairs = 5000, min.pairs.per.cell = 10, verbose = 0, linear.fit = TRUE, local.theta.fit = linear.fit, theta.fit.range = c(1e-2, 1e2)) {
    # default same group
    if(is.null(groups)) {
        groups <- as.factor(rep("cell", ncol(counts)))
    }
    # check for integer counts
    if(any(!unlist(lapply(counts,is.integer)))) {
      stop("Some of the supplied counts are not integer values (or stored as non-integer types). Aborting!\nThe method is designed to work on read counts - do not pass normalized read counts (e.g. FPKM values). If matrix contains read counts, but they are stored as numeric values, use counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) to recast.");
    }

    # crossfit
    if(verbose) {
        cat("cross-fitting cells.\n")
    }
    cfm <- calculate.crossfit.models(counts, groups, n.cores = n.cores, threshold.segmentation = threshold.segmentation, min.count.threshold = min.count.threshold, zero.lambda = zero.lambda, max.pairs = max.pairs, save.plots = save.crossfit.plots, min.pairs.per.cell = min.pairs.per.cell, verbose = verbose)
    # error model for each cell
    if(verbose) {
        cat("building individual error models.\n")
    }
    ifm <- calculate.individual.models(counts, groups, cfm, min.nonfailed = min.nonfailed, zero.count.threshold = zero.count.threshold, n.cores = n.cores, save.plots = save.model.plots, linear.fit = linear.fit, return.compressed.models = TRUE, verbose = verbose, min.size.entries = min.size.entries, local.theta.fit = local.theta.fit, theta.fit.range = theta.fit.range)
    rm(cfm)
    gc()
    return(ifm)
}


##' Estimate prior distribution for gene expression magnitudes
##'
##' Use existing count data to determine a prior distribution of genes in the dataset
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts count matrix
##' @param length.out number of points (resolution) of the expression magnitude grid (default: 400). Note: larger numbers will linearly increase memory/CPU demands.
##' @param show.plot show the estimate posterior
##' @param pseudo.count pseudo-count value to use (default 1)
##' @param bw smoothing bandwidth to use in estimating the prior (default: 0.1)
##' @param max.quantile determine the maximum expression magnitude based on a quantile (default : 0.999)
##' @param max.value alternatively, specify the exact maximum expression magnitude value
##'
##' @return a structure describing expression magnitude grid ($x, on log10 scale) and prior ($y)
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##'
##' @export
scde.expression.prior <- function(models, counts, length.out = 400, show.plot = FALSE, pseudo.count = 1, bw = 0.1, max.quantile = 1-1e-3, max.value = NULL) {
    fpkm <- scde.expression.magnitude(models, counts)
    fail <- scde.failure.probability(models, counts = counts)
    fpkm <- log10(exp(as.matrix(fpkm))+1)
    wts <- as.numeric(as.matrix(1-fail[, colnames(fpkm)]))
    wts <- wts/sum(wts)

    # fit density on a mirror image
    if(is.null(max.value)) {
        x <- as.numeric(fpkm)
        max.value <- as.numeric(quantile(x[x<Inf], p = max.quantile))
    }
    md <- density(c(-1*as.numeric(fpkm), as.numeric(fpkm)), bw = bw, weights = c(wts/2, wts/2), n = 2*length.out+1, from = -1*max.value, to = max.value)

    gep <- data.frame(x = md$x[-seq_len(length.out)], y = md$y[-seq_len(length.out)])
    gep$y[is.na(gep$y)] <- 0
    gep$y <- gep$y+pseudo.count/nrow(fpkm) # pseudo-count
    gep$y <- gep$y/sum(gep$y)
    if(show.plot) {
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
        plot(gep$x, gep$y, col = 4, panel.first = abline(h = 0, lty = 2), type = 'l', xlab = "log10( signal+1 )", ylab = "probability density", main = "signal prior")
    }
    gep$lp <- log(gep$y)

    # grid weighting (for normalization)
    gep$grid.weight <- diff(10^c(gep$x[1], gep$x+c(diff(gep$x)/2, 0))-1)

    return(gep)
    plot(x)
}


##' Test for expression differences between two sets of cells
##'
##' Use the individual cell error models to test for differential expression between two groups of cells.
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts read count matrix
##' @param prior gene expression prior as determined by \code{\link{scde.expression.prior}}
##' @param groups a factor determining the two groups of cells being compared. The factor entries should correspond to the rows of the model matrix. The factor should have two levels. NAs are allowed (cells will be omitted from comparison).
##' @param batch a factor (corresponding to rows of the model matrix) specifying batch assignment of each cell, to perform batch correction
##' @param n.randomizations number of bootstrap randomizations to be performed
##' @param n.cores number of cores to utilize
##' @param batch.models (optional) separate models for the batch data (if generated using batch-specific group argument). Normally the same models are used.
##' @param return.posteriors whether joint posterior matrices should be returned
##' @param verbose integer verbose level (1 for verbose)
##'
##' @return \subsection{default}{
##' a data frame with the following fields:
##' \itemize{
##' \item{lb, mle, ub} {lower bound, maximum likelihood estimate, and upper bound of the 95% confidence interval for the expression fold change on log2 scale.}
##' \item{ce} { conservative estimate of expression-fold change (equals to the min(abs(c(lb, ub))), or 0 if the CI crosses the 0}
##' \item{Z} { uncorrected Z-score of expression difference}
##' \item{cZ} {expression difference Z-score corrected for multiple hypothesis testing using Holm procedure}
##' }
##'  If batch correction has been performed (\code{batch} has been supplied), analogous data frames are returned in slots \code{$batch.adjusted} for batch-corrected results, and \code{$batch.effect} for the differences explained by batch effects alone.
##' }}
##' \subsection{return.posteriors = TRUE}{
##' A list is returned, with the default results data frame given in the \code{$results} slot.
##' \code{difference.posterior} returns a matrix of estimated expression difference posteriors (rows - genes, columns correspond to different magnitudes of fold-change - log2 values are given in the column names)
##' \code{joint.posteriors} a list of two joint posterior matrices (rows - genes, columns correspond to the expression levels, given by prior$x grid)
##' }
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' \donttest{
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # make sure groups corresponds to the models (o.ifm)
##' groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
##' names(groups) <- row.names(o.ifm)
##' ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = n.cores, verbose = 1)
##' }
##'
##' @export
scde.expression.difference <- function(models, counts, prior, groups = NULL, batch = NULL, n.randomizations = 150, n.cores = 10, batch.models = models, return.posteriors = FALSE, verbose = 0) {
    if(!all(rownames(models) %in% colnames(counts))) {
        stop("ERROR: provided count data does not cover all of the cells specified in the model matrix")
    }

    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[, ci])

    if(is.null(groups)) { # recover groups from models
        groups <- as.factor(attr(models, "groups"))
        if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
        names(groups) <- rownames(models)
    }
    if(length(levels(groups)) != 2) {
        stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
    }

    correct.batch <- FALSE
    if(!is.null(batch)) {
        if(length(levels(batch)) > 1) {
            correct.batch <- TRUE
        } else {
            if(verbose) {
                cat("WARNING: only one batch level detected. Nothing to correct for.")
            }
        }
    }

    # batch control
    if(correct.batch) {
        batch <- as.factor(batch)
        # check batch-group interactions
        bgti <- table(groups, batch)
        bgti.ft <- fisher.test(bgti)
        if(verbose) {
            cat("controlling for batch effects. interaction:\n")
            print(bgti)
        }
        #if(any(bgti == 0)) {
        #  cat("ERROR: cannot control for batch effect, as some batches are found only in one group:\n")
        #  print(bgti)
        #}
        if(bgti.ft$p.value < 1e-3) {
            cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
            print(bgti.ft)
        }

        # calculate batch posterior
        if(verbose) {
            cat("calculating batch posteriors\n")
        }
        batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
            scde.posteriors(models = batch.models, counts = counts, prior = prior, batch = batch, composition = table(batch[ii]), n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = FALSE)
        })
        if(verbose) {
            cat("calculating batch differences\n")
        }
        batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], batch.jpl[[2]], prior, n.cores = n.cores)
        batch.bdiffp.rep <- quick.distribution.summary(batch.bdiffp)
    } else {
        if(verbose) {
            cat("comparing groups:\n")
            print(table(as.character(groups)))
        }
    }


    # fit joint posteriors for each group
    jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
        scde.posteriors(models = models[ii, , drop = FALSE], counts = counts[, ii, drop = FALSE], prior = prior, n.cores = n.cores, n.randomizations = n.randomizations)
    })
    if(verbose) {
        cat("calculating difference posterior\n")
    }
    # calculate difference posterior
    bdiffp <- calculate.ratio.posterior(jpl[[1]], jpl[[2]], prior, n.cores = n.cores)

    if(verbose) {
        cat("summarizing differences\n")
    }
    bdiffp.rep <- quick.distribution.summary(bdiffp)

    if(correct.batch) {
        if(verbose) {
            cat("adjusting for batch effects\n")
        }
        # adjust for batch effects
        a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, prior = data.frame(x = as.numeric(colnames(bdiffp)), y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE, n.cores = n.cores)
        a.bdiffp.rep <- quick.distribution.summary(a.bdiffp)

        # return with batch correction info
        if(return.posteriors) {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep, difference.posterior = bdiffp, batch.adjusted.difference.posterior = a.bdiffp, joint.posteriors = jpl))
        } else {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep))
        }
    } else {
        # no batch correction return
        if(return.posteriors) {
            return(list(results = bdiffp.rep, difference.posterior = bdiffp, joint.posteriors = jpl))
        } else {
            return(bdiffp.rep)
        }
    }
}


##' View differential expression results in a browser
##'
##' Launches a browser app that shows the differential expression results, allowing to sort, filter, etc.
##' The arguments generally correspond to the \code{scde.expression.difference()} call, except that the results of that call are also passed here. Requires \code{Rook} and \code{rjson} packages to be installed.
##'
##' @param results result object returned by \code{scde.expression.difference()}. Note to browse group posterior levels, use \code{return.posteriors = TRUE} in the \code{scde.expression.difference()} call.
##' @param models model matrix
##' @param counts count matrix
##' @param prior prior
##' @param groups group information
##' @param batch batch information
##' @param geneLookupURL The URL that will be used to construct links to view more information on gene names. By default (if can't guess the organism) the links will forward to ENSEMBL site search, using \code{geneLookupURL = "http://useast.ensembl.org/Multi/Search/Results?q = {0}"}. The "{0}" in the end will be substituted with the gene name. For instance, to link to GeneCards, use \code{"http://www.genecards.org/cgi-bin/carddisp.pl?gene = {0}"}.
##' @param server optional previously returned instance of the server, if want to reuse it.
##' @param name app name (needs to be altered only if adding more than one app to the server using \code{server} parameter)
##' @param port Interactive browser port
##'
##' @return server instance, on which $stop() function can be called to kill the process.
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' \donttest{
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # make sure groups corresponds to the models (o.ifm)
##' groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
##' names(groups) <- row.names(o.ifm)
##' ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 10, verbose = 1)
##' scde.browse.diffexp(ediff, o.ifm, cd, o.prior, groups = groups, geneLookupURL="http://www.informatics.jax.org/searchtool/Search.do?query={0}")  # creates browser
##' }
##'
##' @export
scde.browse.diffexp <- function(results, models, counts, prior, groups = NULL, batch = NULL, geneLookupURL = NULL, server = NULL, name = "scde", port = NULL) {
    #require(Rook)
    #require(rjson)
    if(is.null(server)) { server <- get.scde.server(port) }
    sa <- ViewDiff$new(results, models, counts, prior, groups = groups, batch = batch, geneLookupURL = geneLookupURL)
    server$add(app = sa, name = name)
    browseURL(paste(server$full_url(name), "index.html", sep = "/"))
    return(server)
}


##' View PAGODA application
##'
##' Installs a given pagoda app (or any other rook app) into a server, optionally
##' making a call to show it in the browser.
##'
##' @param app pagoda app (output of make.pagoda.app()) or another rook app
##' @param name URL path name for this app
##' @param browse whether a call should be made for browser to show the app
##' @param port optional port on which the server should be initiated
##' @param ip IP on which the server should listen (typically localhost)
##' @param server an (optional) Rook server instance (defaults to ___scde.server)
##'
##' @examples
##' \donttest{
##' app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols=col.cols, cell.clustering=hc, title="NPCs")
##' # show app in the browser (port 1468)
##' show.app(app, "pollen", browse = TRUE, port=1468)
##' }
##'
##' @return Rook server instance
##'
##' @export
show.app <- function(app, name, browse = TRUE, port = NULL, ip = '127.0.0.1', server = NULL) {
    # replace special characters
    name <- gsub("[^[:alnum:.]]", "_", name)
    
    if (tools:::httpdPort() !=0 && tools:::httpdPort() != port) {
        cat("ERROR: port is already being used. The PAGODA app is currently incompatible with RStudio. Please try running the interactive app in the R console.")
    }
    if(is.null(server)) { server <- get.scde.server(port) }
    server$add(app = app, name = name)
    if(browse) {
        browseURL(paste(server$full_url(name), "index.html", sep = "/"))
    }
    return(server)
}
# get SCDE server from saved session
get.scde.server <- function(port = NULL, ip = '127.0.0.1') {
    if(exists("___scde.server", envir = globalenv())) {
        server <- get("___scde.server", envir = globalenv())
    } else {
        require(Rook)
        server <- Rhttpd$new()
        assign("___scde.server", server, envir = globalenv())
        server$start(listen = ip, port = port)
    }
    return(server)
}


# calculate individual and joint posterior information
# models - all or a subset of models belonging to a particular group
#
##' Calculate joint expression magnitude posteriors across a set of cells
##'
##' Calculates expression magnitude posteriors for the individual cells, and then uses bootstrap resampling to calculate a joint expression posterior for all the specified cells. Alternatively during batch-effect correction procedure, the joint posterior can be calculated for a random composition of cells of different groups (see \code{batch} and \code{composition} parameters).
##'
##' @param models models models determined by \code{\link{scde.error.models}}
##' @param counts read count matrix
##' @param prior gene expression prior as determined by \code{\link{scde.expression.prior}}
##' @param n.randomizations number of bootstrap iterations to perform
##' @param batch a factor describing which batch group each cell (i.e. each row of \code{models} matrix) belongs to
##' @param composition a vector describing the batch composition of a group to be sampled
##' @param return.individual.posteriors whether expression posteriors of each cell should be returned
##' @param return.individual.posterior.modes whether modes of expression posteriors of each cell should be returned
##' @param ensemble.posterior Boolean of whether to calculate the ensemble posterior (sum of individual posteriors) instead of a joint (product) posterior. (default: FALSE)
##' @param n.cores number of cores to utilize
##'
##' @return \subsection{default}{ a posterior probability matrix, with rows corresponding to genes, and columns to expression levels (as defined by \code{prior$x})
##' }
##' \subsection{return.individual.posterior.modes}{ a list is returned, with the \code{$jp} slot giving the joint posterior matrix, as described above. The \code{$modes} slot gives a matrix of individual expression posterior mode values on log scale (rows - genes, columns -cells)}
##' \subsection{return.individual.posteriors}{ a list is returned, with the \code{$post} slot giving a list of individual posterior matrices, in a form analogous to the joint posterior matrix, but reported on log scale }
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate joint posteriors
##' jp <- scde.posteriors(o.ifm, cd, o.prior, n.cores = 1)
##'
##' @export
scde.posteriors <- function(models, counts, prior, n.randomizations = 100, batch = NULL, composition = NULL, return.individual.posteriors = FALSE, return.individual.posterior.modes = FALSE, ensemble.posterior = FALSE, n.cores = 20) {
    if(!all(rownames(models) %in% colnames(counts))) { stop("ERROR: provided count data does not cover all of the cells specified in the model matrix") }
    if(!is.null(batch)) { # calculating batch-sampled posteriors instead of evenly sampled ones
        if(is.null(composition)) { stop("ERROR: group composition must be provided if the batch argument is passed") }
        batchil <- tapply(c(1:nrow(models))-1, batch, I)
    }
    # order counts according to the cells
    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[, ci, drop = FALSE])
    marginals <- 10^prior$x - 1
    marginals[marginals<0] <- 0
    marginals <- log(marginals)

    min.slope <- 1e-10
    if(any(models$corr.a<min.slope)) {
        cat("WARNING: the following cells have negatively-correlated or 0-slope fits: ", paste(rownames(models)[models$corr.a<min.slope], collapse = " "), ". Setting slopes to 1e-10.\n")
        models$corr.a[models$corr.a<min.slope] <- min.slope
    }

    postflag <- 0
    if(return.individual.posteriors) {
        postflag <- 2
        if(return.individual.posterior.modes) {
            postflag <- 3
        }
    } else if(return.individual.posterior.modes) {
        postflag <- 1
    }

    ensembleflag <- ifelse(ensemble.posterior, 1, 0)

    localthetaflag <- "corr.ltheta.b" %in% colnames(models)
    squarelogitconc <- "conc.a2" %in% colnames(models)

    # prepare matrix models
    mn <- c("conc.b", "conc.a", "fail.r", "corr.b", "corr.a", "corr.theta", "corr.ltheta.b", "corr.ltheta.t", "corr.ltheta.m", "corr.ltheta.s", "corr.ltheta.r", "conc.a2")
    mc <- match(c(mn), colnames(models))
    mm <- matrix(NA, nrow(models), length(mn))
    mm[, which(!is.na(mc))] <- as.matrix(models[, mc[!is.na(mc)], drop = FALSE])

    chunk <- function(x, n) split(x, sort(rank(x) %% n.cores))
    if(n.cores > 1 && nrow(counts) > n.cores) { # split by genes
        xl <- papply(chunk(seq_len(nrow(counts)), n.cores), function(ii) {
            ucl <- lapply(seq_len(ncol(counts)), function(i) as.vector(unique(counts[ii, i, drop = FALSE])))
            uci <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) match(counts[ii, i, drop = FALSE], ucl[[i]])-1))
            #x <- logBootPosterior(models, ucl, uci, marginals, n.randomizations, 1, postflag)
            if(!is.null(batch)) {
                x <- .Call("logBootBatchPosterior", mm, ucl, uci, marginals, batchil, composition, n.randomizations, ii[1], postflag, localthetaflag, squarelogitconc, PACKAGE = "scde")
            } else {
                x <- .Call("logBootPosterior", mm, ucl, uci, marginals, n.randomizations, ii[1], postflag, localthetaflag, squarelogitconc, ensembleflag, PACKAGE = "scde")
            }
        }, n.cores = n.cores)
        if(postflag == 0) {
            x <- do.call(rbind, xl)
        } else if(postflag == 1) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), modes = do.call(rbind, lapply(xl, function(d) d$modes)))
        } else if(postflag == 2) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), post = lapply(seq_along(xl[[1]]$post), function(pi) { do.call(rbind, lapply(xl, function(d) d$post[[pi]])) }))
        } else if(postflag == 3) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), modes = do.call(rbind, lapply(xl, function(d) d$modes)), post = lapply(seq_along(xl[[1]]$post), function(pi) { do.call(rbind, lapply(xl, function(d) d$post[[pi]])) }))
        }
        rm(xl)
        gc()
    } else {
        # unique count lists with matching indices
        ucl <- lapply(seq_len(ncol(counts)), function(i) as.vector(unique(counts[, i, drop = FALSE])))
        uci <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) match(counts[, i, drop = FALSE], ucl[[i]])-1))
        #x <- logBootPosterior(models, ucl, uci, marginals, n.randomizations, 1, postflag)
        if(!is.null(batch)) {
            x <- .Call("logBootBatchPosterior", mm, ucl, uci, marginals, batchil, composition, n.randomizations, 1, postflag, localthetaflag, squarelogitconc, PACKAGE = "scde")
        } else {
            x <- .Call("logBootPosterior", mm, ucl, uci, marginals, n.randomizations, 1, postflag, localthetaflag, squarelogitconc, ensembleflag, PACKAGE = "scde")
        }
    }
    if(postflag == 0) {
        rownames(x) <- rownames(counts)
        colnames(x) <- as.character(exp(marginals))
    } else if(postflag == 1) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        rownames(x$modes) <- rownames(counts)
        colnames(x$modes) <- rownames(models)
    } else if(postflag == 2) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        names(x$post) <- rownames(models)
        x$post <- lapply(x$post, function(d) {
            rownames(d) <- rownames(counts)
            colnames(d) <- as.character(exp(marginals))
            return(d)
        })
    } else if(postflag == 3) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        rownames(x$modes) <- rownames(counts)
        colnames(x$modes) <- rownames(models)
        names(x$post) <- rownames(models)
        x$post <- lapply(x$post, function(d) {
            rownames(d) <- rownames(counts)
            colnames(d) <- as.character(exp(marginals))
            return(d)
        })
    }
    return(x)
}


# get estimates of expression magnitude for a given set of models
# models - entire model matrix, or a subset of cells (i.e. select rows) of the model matrix for which the estimates should be obtained
# counts - count data that covers the desired set of genes (rows) and all specified cells (columns)
# return - a matrix of log(FPM) estimates with genes as rows and cells  as columns (in the model matrix order).
##' Return scaled expression magnitude estimates
##'
##' Return point estimates of expression magnitudes of each gene across a set of cells, based on the regression slopes determined during the model fitting procedure.
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts count matrix
##'
##' @return a matrix of expression magnitudes on a log scale (rows - genes, columns - cells)
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' # get expression magnitude estimates
##' lfpm <- scde.expression.magnitude(o.ifm, cd)
##'
##' @export
scde.expression.magnitude <- function(models, counts) {
    if(!all(rownames(models) %in% colnames(counts))) { stop("ERROR: provided count data does not cover all of the cells specified in the model matrix") }
    t((t(log(counts[, rownames(models), drop = FALSE]))-models$corr.b)/models$corr.a)
}


# calculate drop-out probability given either count data or magnitudes (log(FPM))
# magnitudes can either be a per-cell matrix or a single vector of values which will be evaluated for each cell
# returns a probability of a drop out event for every gene (rows) for every cell (columns)
##' Calculate drop-out probabilities given a set of counts or expression magnitudes
##'
##' Returns estimated drop-out probability for each cell (row of \code{models} matrix), given either an expression magnitude
##' @param models models determined by \code{\link{scde.error.models}}
##' @param magnitudes a vector (\code{length(counts) == nrows(models)}) or a matrix (columns correspond to cells) of expression magnitudes, given on a log scale
##' @param counts a vector (\code{length(counts) == nrows(models)}) or a matrix (columns correspond to cells) of read counts from which the expression magnitude should be estimated
##'
##' @return a vector or a matrix of drop-out probabilities
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate probability of observing a drop out at a given set of magnitudes in different cells
##' mags <- c(1.0, 1.5, 2.0)
##' p <- scde.failure.probability(o.ifm, magnitudes = mags)
##' # calculate probability of observing the dropout at a magnitude corresponding to the
##' # number of reads actually observed in each cell
##' self.p <- scde.failure.probability(o.ifm, counts = cd)
##'
##' @export
scde.failure.probability <- function(models, magnitudes = NULL, counts = NULL) {
    if(is.null(magnitudes)) {
        if(!is.null(counts)) {
            magnitudes <- scde.expression.magnitude(models, counts)
        } else {
            stop("ERROR: either magnitudes or counts should be provided")
        }
    }
    if(is.matrix(magnitudes)) { # a different vector for every cell
        if(!all(rownames(models) %in% colnames(magnitudes))) { stop("ERROR: provided magnitude data does not cover all of the cells specified in the model matrix") }
        if("conc.a2" %in% names(models)) {
            x <- t(1/(exp(t(magnitudes)*models$conc.a +t(magnitudes^2)*models$conc.a2 + models$conc.b)+1))
        } else {
            x <- t(1/(exp(t(magnitudes)*models$conc.a + models$conc.b)+1))
        }
    } else { # a common vector of magnitudes for all cells
        if("conc.a2" %in% names(models)) {
            x <- t(1/(exp((models$conc.a %*% t(magnitudes)) + (models$conc.a2 %*% t(magnitudes^2)) + models$conc.b)+1))
        } else {
            x <- t(1/(exp((models$conc.a %*% t(magnitudes)) + models$conc.b)+1))
        }
    }
    x[is.nan(x)] <- 0
    colnames(x) <- rownames(models)
    x
}


##' Test differential expression and plot posteriors for a particular gene
##'
##' The function performs differential expression test and optionally plots posteriors for a specified gene.
##'
##' @param gene name of the gene to be tested
##' @param models models
##' @param counts read count matrix (must contain the row corresponding to the specified gene)
##' @param prior expression magnitude prior
##' @param groups a two-level factor specifying between which cells (rows of the models matrix) the comparison should be made
##' @param batch optional multi-level factor assigning the cells (rows of the model matrix) to different batches that should be controlled for (e.g. two or more biological replicates). The expression difference estimate will then take into account the likely difference between the two groups that is explained solely by their difference in batch composition. Not all batch configuration may be corrected this way.
##' @param batch.models optional set of models for batch comparison (typically the same as models, but can be more extensive, or recalculated within each batch)
##' @param n.randomizations number of bootstrap/sampling iterations that should be performed
##' @param show.plots whether the plots should be shown
##' @param return.details whether the posterior should be returned
##' @param verbose set to T for some status output
##' @param ratio.range optionally specifies the range of the log2 expression ratio plot
##' @param show.individual.posteriors whether the individual cell expression posteriors should be plotted
##' @param n.cores number of cores to use (default = 1)
##'
##' @return by default returns MLE of log2 expression difference, 95% CI (upper, lower bound), and a Z-score testing for expression difference. If return.details = TRUE, a list is returned containing the above structure, as well as the expression fold difference posterior itself.
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)
##'
##' @export
scde.test.gene.expression.difference <- function(gene, models, counts, prior, groups = NULL, batch = NULL, batch.models = models, n.randomizations = 1e3, show.plots = TRUE, return.details = FALSE, verbose = FALSE, ratio.range = NULL, show.individual.posteriors = TRUE, n.cores = 1) {
    if(!gene %in% rownames(counts)) {
        stop("ERROR: specified gene (", gene, ") is not found in the count data")
    }

    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[gene, ci, drop = FALSE])


    if(is.null(groups)) { # recover groups from models
        groups <- as.factor(attr(models, "groups"))
        if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
        names(groups) <- rownames(models)
    }
    if(length(levels(groups)) != 2) {
        stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
    }

    if(verbose) {
        cat("comparing gene ", gene, " between groups:\n")
        print(table(as.character(groups)))
    }

    # calculate joint posteriors
    jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
        scde.posteriors(models = models[ii, , drop = FALSE], counts = counts[, ii, drop = FALSE], prior = prior, n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = TRUE)
    })

    bdiffp <- calculate.ratio.posterior(jpl[[1]]$jp, jpl[[2]]$jp, prior, n.cores = n.cores)

    bdiffp.rep <- quick.distribution.summary(bdiffp)

    nam1 <- levels(groups)[1]
    nam2 <- levels(groups)[2]

    # batch control
    correct.batch <- !is.null(batch) && length(levels(batch)) > 1
    if(correct.batch) {
        batch <- as.factor(batch)
        # check batch-group interactions
        bgti <- table(groups, batch)
        bgti.ft <- fisher.test(bgti)
        if(verbose) {
            cat("controlling for batch effects. interaction:\n")
        }
        if(any(bgti == 0)) {
            cat("ERROR: cannot control for batch effect, as some batches are found only in one group:\n")
            print(bgti)
        }
        if(bgti.ft$p.value<1e-3) {
            cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
            print(bgti)
            print(bgti.ft)
        }
        # calculate batch posterior
        batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
            scde.posteriors(models = batch.models, counts = counts, prior = prior, batch = batch, composition = table(batch[ii]), n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = FALSE)
        })
        batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], batch.jpl[[2]], prior, n.cores = n.cores)
        a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, prior = data.frame(x = as.numeric(colnames(bdiffp)), y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE)
        a.bdiffp.rep <- quick.distribution.summary(a.bdiffp)
    }


    if(show.plots) {
        # show each posterior
        layout(matrix(c(1:3), 3, 1, byrow = TRUE), heights = c(2, 1, 2), widths = c(1), FALSE)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        #par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)

        pp <- exp(do.call(rbind, lapply(jpl[[1]]$post, as.numeric)))
        cols <- rainbow(nrow(pp), s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, pp)), xlab = "expression level", ylab = "individual posterior", main = nam1)
        if(show.individual.posteriors) {
            lapply(seq_len(nrow(pp)), function(i) lines(prior$x, pp[i, ], col = rgb(1, 0.5, 0, alpha = 0.25)))
        }
        #legend(x = ifelse(which.max(na.omit(pjpc)) > length(pjpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = rownames(pp), lty = rep(1, nrow(pp)))
        if(correct.batch) {
            par(new = TRUE)
            plot(prior$x, batch.jpl[[1]][1, ], axes = FALSE, ylab = "", xlab = "", type = 'l', col = 8, lty = 1, lwd = 2)
        }
        pjpc <- jpl[[1]]$jp
        par(new = TRUE)
        jpr <- range(c(0, na.omit(pjpc)))
        plot(prior$x, pjpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)


        # ratio plot
        if(is.null(ratio.range)) { ratio.range <- range(as.numeric(colnames(bdiffp))/log10(2)) }

        par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        rv <- as.numeric(colnames(bdiffp))/log10(2)
        rp <- as.numeric(bdiffp[1, ])
        plot(rv, rp, xlab = "log2 expression ratio", ylab = "ratio posterior", type = 'l', lwd = ifelse(correct.batch, 1, 2), main = "", axes = FALSE, xlim = ratio.range, ylim = c(0, max(bdiffp)))
        axis(1, pretty(ratio.range, 5), col = 1)
        abline(v = 0, lty = 2, col = 8)
        if(correct.batch) { # with batch correction
            # show batch difference
            par(new = TRUE)
            plot(as.numeric(colnames(batch.bdiffp))/log10(2), as.numeric(batch.bdiffp[1, ]), xlab = "", ylab = "", type = 'l', lwd = 1, main = "", axes = FALSE, xlim = ratio.range, col = 8, ylim = c(0, max(batch.bdiffp)))
            # fill out the a.bdiffp confidence interval
            par(new = TRUE)
            rv <- as.numeric(colnames(a.bdiffp))/log10(2)
            rp <- as.numeric(a.bdiffp[1, ])
            plot(rv, rp, xlab = "", ylab = "", type = 'l', lwd = 2, main = "", axes = FALSE, xlim = ratio.range, col = 2, ylim = c(0, max(rp)))
            axis(2, pretty(c(0, max(a.bdiffp)), 2), col = 1)
            r.lb <- which.min(abs(rv-a.bdiffp.rep$lb))
            r.ub <- which.min(abs(rv-a.bdiffp.rep$ub))
            polygon(c(rv[r.lb], rv[r.lb:r.ub], rv[r.ub]), y = c(-10, rp[r.lb:r.ub], -10), col = rgb(1, 0, 0, alpha = 0.2), border = NA)
            abline(v = a.bdiffp.rep$mle, col = 2, lty = 2)
            abline(v = c(rv[r.ub], rv[r.lb]), col = 2, lty = 3)

            legend(x = ifelse(a.bdiffp.rep$mle > 0, "topleft", "topright"), legend = c(paste("MLE: ", round(a.bdiffp.rep$mle, 2), sep = ""), paste("95% CI: ", round(a.bdiffp.rep$lb, 2), " : ", round(a.bdiffp.rep$ub, 2), sep = ""), paste("Z = ", round(a.bdiffp.rep$Z, 2), sep = ""), paste("cZ = ", round(a.bdiffp.rep$cZ, 2), sep = "")), bty = "n")

        } else {  # without batch correction
            # fill out the bdiffp confidence interval
            axis(2, pretty(c(0, max(bdiffp)), 2), col = 1)

            r.lb <- which.min(abs(rv-bdiffp.rep$lb))
            r.ub <- which.min(abs(rv-bdiffp.rep$ub))
            polygon(c(rv[r.lb], rv[r.lb:r.ub], rv[r.ub]), y = c(-10, rp[r.lb:r.ub], -10), col = rgb(1, 0, 0, alpha = 0.2), border = NA)
            abline(v = bdiffp.rep$mle, col = 2, lty = 2)
            abline(v = c(rv[r.ub], rv[r.lb]), col = 2, lty = 3)

            legend(x = ifelse(bdiffp.rep$mle > 0, "topleft", "topright"), legend = c(paste("MLE: ", round(bdiffp.rep$mle, 2), sep = ""), paste("95% CI: ", round(bdiffp.rep$lb, 2), " : ", round(bdiffp.rep$ub, 2), sep = ""), paste("Z = ", round(bdiffp.rep$Z, 2), sep = ""), paste("aZ = ", round(bdiffp.rep$cZ, 2), sep = "")), bty = "n")
        }

        # distal plot
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        #par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        dp <- exp(do.call(rbind, lapply(jpl[[2]]$post, as.numeric)))
        cols <- rainbow(nrow(dp), s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, dp)), xlab = "expression level", ylab = "individual posterior", main = nam2)
        if(show.individual.posteriors) {
            lapply(seq_len(nrow(dp)), function(i) lines(prior$x, dp[i, ], col = rgb(0, 0.5, 1, alpha = 0.25)))
        }
        if(correct.batch) {
            par(new = TRUE)
            plot(prior$x, batch.jpl[[2]][1, ], axes = FALSE, ylab = "", xlab = "", type = 'l', col = 8, lty = 1, lwd = 2)
        }
        djpc <- jpl[[2]]$jp
        #legend(x = ifelse(which.max(na.omit(djpc)) > length(djpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = rownames(dp), lty = rep(1, nrow(dp)))
        par(new = TRUE)
        jpr <- range(c(0, na.omit(djpc)))
        plot(prior$x, djpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)
    }

    if(return.details) {
        if(correct.batch) { # with batch correction
            return(list(results = a.bdiffp.rep, difference.posterior = a.bdiffp, results.nobatchcorrection = bdiffp.rep))
        } else {
            return(list(results = bdiffp.rep, difference.posterior = bdiffp, posteriors = jpl))
        }
    } else {
        if(correct.batch) { # with batch correction
            return(a.bdiffp.rep)
        } else {
            return(bdiffp.rep)
        }
    }
}


# fit models to external (bulk) reference
##' Fit scde models relative to provided set of expression magnitudes
##'
##' If group-average expression magnitudes are available (e.g. from bulk measurement), this method can be used
##' to fit individual cell error models relative to that reference
##'
##' @param counts count matrix
##' @param reference a vector of expression magnitudes (read counts) corresponding to the rows of the count matrix
##' @param min.fpm minimum reference fpm of genes that will be used to fit the models (defaults to 1). Note: fpm is calculated from the reference count vector as reference/sum(reference)*1e6
##' @param n.cores number of cores to use
##' @param zero.count.threshold read count to use as an initial guess for the zero threshold
##' @param nrep number independent of mixture fit iterations to try (default = 1)
##' @param save.plots whether to write out a pdf file showing the model fits
##' @param plot.filename model fit pdf filename
##' @param verbose verbose level
##'
##' @return matrix of scde models
##'
##' @examples
##' data(es.mef.small)
##' cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
##' \donttest{
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate joint posteriors across all cells
##' jp <- scde.posteriors(models = o.ifm, cd, o.prior, n.cores = 10, return.individual.posterior.modes = TRUE, n.randomizations = 100)
##' # use expected expression magnitude for each gene
##' av.mag <- as.numeric(jp$jp %*% as.numeric(colnames(jp$jp)))
##' # translate into counts
##' av.mag.counts <- as.integer(round(av.mag))
##' # now, fit alternative models using av.mag as a reference (normally this would correspond to bulk RNA expression magnitude)
##' ref.models <- scde.fit.models.to.reference(cd, av.mag.counts, n.cores = 1)
##' }
##'
##' @export
scde.fit.models.to.reference <- function(counts, reference, n.cores = 10, zero.count.threshold = 1, nrep = 1, save.plots = FALSE, plot.filename = "reference.model.fits.pdf", verbose = 0, min.fpm = 1) {
    return.compressed.models <- TRUE
    verbose <- 1
    ids <- colnames(counts)
    ml <- papply(seq_along(ids), function(i) {
        df <- data.frame(count = counts[, ids[i]], fpm = reference/sum(reference)*1e6)
        df <- df[df$fpm > min.fpm, ]
        m1 <- fit.nb2.mixture.model(df, nrep = nrep, verbose = verbose, zero.count.threshold = zero.count.threshold)
        if(return.compressed.models) {
            v <- get.compressed.v1.model(m1)
            cl <- clusters(m1)
            rm(m1)
            gc()
            return(list(model = v, clusters = cl))
        } else {
            return(m1)
        }
    }, n.cores = n.cores)
    names(ml) <- ids

    # check if there were errors in the multithreaded portion
    lapply(seq_along(ml), function(i) {
        if(class(ml[[i]]) == "try-error") {
            message("ERROR encountered in building a model for cell ", ids[i], ":")
            message(ml[[i]])
            tryCatch(stop(paste("ERROR encountered in building a model for cell ", ids[i])), error = function(e) stop(e))
        }
    })

    if(save.plots) {
        # model fits
        #CairoPNG(file = paste(group, "model.fits.png", sep = "."), width = 1024, height = 300*length(ids))
        pdf(file = plot.filename, width = 13, height = 4)
        #l <- layout(matrix(seq(1, 4*length(ids)), nrow = length(ids), byrow = TRUE), rep(c(1, 1, 1, 0.5), length(ids)), rep(1, 4*length(ids)), FALSE)
        l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, 0.5), 1), rep(1, 4), FALSE)
        par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
        invisible(lapply(seq_along(ids), function(i) {
            df <- data.frame(count = counts[, ids[i]], fpm = reference/sum(reference)*1e6)
            df <- df[df$fpm > min.fpm, ]
            plot.nb2.mixture.fit(ml[[i]], df, en = ids[i], do.par = FALSE, compressed.models = return.compressed.models)
        }))
        dev.off()
    }

    if(return.compressed.models) {
        # make a joint model matrix
        jmm <- data.frame(do.call(rbind, lapply(ml, function(m) m$model)))
        rownames(jmm) <- names(ml)
        jmm
        return(jmm)
    } else {
        return(ml)
    }
}


##' Determine principal components of a matrix using per-observation/per-variable weights
##'
##' Implements a weighted PCA
##'
##' @param mat matrix of variables (columns) and observations (rows)
##' @param matw  corresponding weights
##' @param npcs number of principal components to extract
##' @param nstarts number of random starts to use
##' @param smooth smoothing span
##' @param em.tol desired EM algorithm tolerance
##' @param em.maxiter maximum number of EM iterations
##' @param seed random seed
##' @param center whether mat should be centered (weighted centering)
##' @param n.shuffles optional number of per-observation randomizations that should be performed in addition to the main calculations to determine the lambda1 (PC1 eigenvalue) magnitude under such randomizations (returned in $randvar)
##'
##' @return a list containing eigenvector matrix ($rotation), projections ($scores), variance (weighted) explained by each component ($var), total (weighted) variance of the dataset ($totalvar)
##'
##' @examples
##' set.seed(0)
##' mat <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random matrix
##' base.pca <- bwpca(mat)  # non-weighted pca, equal weights set automatically
##' matw <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random weight matrix
##' matw <- abs(matw)/max(matw)
##' base.pca.weighted <- bwpca(mat, matw)  # weighted pca
##'
##' @export
bwpca <- function(mat, matw = NULL, npcs = 2, nstarts = 1, smooth = 0, em.tol = 1e-6, em.maxiter = 25, seed = 1, center = TRUE, n.shuffles = 0) {
    if(smooth<4) { smooth <- 0 }
    if(any(is.nan(matw))) {
      stop("bwpca: weight matrix contains NaN values")
    }
    if(any(is.nan(mat))) {
      stop("bwpca: value matrix contains NaN values")
    }
    if(is.null(matw)) {
        matw <- matrix(1, nrow(mat), ncol(mat))
        nstarts <- 1
    }
    if(center) { mat <- t(t(mat)-colSums(mat*matw)/colSums(matw)) }

    res <- .Call("baileyWPCA", mat, matw, npcs, nstarts, smooth, em.tol, em.maxiter, seed, n.shuffles, PACKAGE = "scde")
    #res <- bailey.wpca(mat, matw, npcs, nstarts, smooth, em.tol, em.maxiter, seed)
    rownames(res$rotation) <- colnames(mat)
    rownames(res$scores) <- rownames(mat)
    colnames(res$rotation) <- paste("PC", seq(1:ncol(res$rotation)), sep = "")
    res$sd <- t(sqrt(res$var))
    res
}


##' Winsorize matrix
##'
##' Sets the ncol(mat)*trim top outliers in each row to the next lowest value same for the lowest outliers
##'
##' @param mat matrix
##' @param trim fraction of outliers (on each side) that should be Winsorized, or (if the value is  >= 1) the number of outliers to be trimmed on each side
##'
##' @return Winsorized matrix
##'
##' @examples
##' set.seed(0)
##' mat <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random matrix
##' mat[1,1] <- 1000  # make outlier
##' range(mat)  # look at range of values
##' win.mat <- winsorize.matrix(mat, 0.1)
##' range(win.mat)  # note outliers removed
##'
##' @export
winsorize.matrix <- function(mat, trim) {
    if(trim  >  0.5) { trim <- trim/ncol(mat)  }
    wm <- .Call("winsorizeMatrix", mat, trim, PACKAGE = "scde")
    rownames(wm) <- rownames(mat)
    colnames(wm) <- colnames(mat)
    return(wm)
}

