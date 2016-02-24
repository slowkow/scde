##' Build error models for heterogeneous cell populations, based on K-nearest neighbor cells.
##'
##' Builds cell-specific error models assuming that there are multiple subpopulations present
##' among the measured cells. The models for each cell are based on average expression estimates
##' obtained from K closest cells within a given group (if groups = NULL, then within the entire
##' set of measured cells). The method implements fitting of both the original log-fit models
##' (when linear.fit = FALSE), or newer linear-fit models (linear.fit = TRUE, default) with locally
##' fit overdispersion coefficient (local.theta.fit = TRUE, default).
##'
##' @param counts count matrix (integer matrix, rows- genes, columns- cells)
##' @param groups optional groups partitioning known subpopulations
##' @param cor.method correlation measure to be used in determining k nearest cells
##' @param k number of nearest neighbor cells to use during fitting. If k is set sufficiently high, all of the cells within a given group will be used.
##' @param min.nonfailed minimum number of non-failed measurements (within the k nearest neighbor cells) required for a gene to be taken into account during error fitting procedure
##' @param min.size.entries minimum number of genes to use for model fitting
##' @param min.count.threshold minimum number of reads required for a measurement to be considered non-failed
##' @param save.model.plots whether model plots should be saved (file names are (group).models.pdf, or cell.models.pdf if no group was supplied)
##' @param max.model.plots maximum number of models to save plots for (saves time when there are too many cells)
##' @param n.cores number of cores to use through the calculations
##' @param min.fpm optional parameter to restrict model fitting to genes with group-average expression magnitude above a given value
##' @param verbose level of verbosity
##' @param fpm.estimate.trim trim fraction to be used in estimating group-average gene expression magnitude for model fitting (0.5 would be median, 0 would turn off trimming)
##' @param linear.fit whether newer linear model fit with zero intercept should be used (T), or the log-fit model published originally (F)
##' @param local.theta.fit whether local theta fitting should be used (only available for the linear fit models)
##' @param theta.fit.range allowed range of the theta values
##' @param alpha.weight.power 1/theta weight power used in fitting theta dependency on the expression magnitude
##'
##' @return a data frame with parameters of the fit error models (rows- cells, columns- fitted parameters)
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' }
##'
##' @export
knn.error.models <- function(counts, groups = NULL, k = round(ncol(counts)/2), min.nonfailed = 5, min.count.threshold = 1, save.model.plots = TRUE, max.model.plots = 50, n.cores = parallel::detectCores(), min.size.entries = 2e3, min.fpm = 0, cor.method = "pearson", verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE, local.theta.fit = linear.fit, theta.fit.range = c(1e-2, 1e2), alpha.weight.power = 1/2) {
    threshold.prior = 1-1e-6

    # check for integer counts
    if(any(!unlist(lapply(counts,is.integer)))) {
      stop("Some of the supplied counts are not integer values (or stored as non-integer types). Aborting!\nThe method is designed to work on read counts - do not pass normalized read counts (e.g. FPKM values). If matrix contains read counts, but they are stored as numeric values, use counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) to recast.");
    }

    # TODO:
    #  - implement check for k >= n.cells (to avoid correlation calculations)
    #  - implement error reporting/handling for failed cell fits

    if(is.null(groups)) {
        groups <- as.factor(rep("cell", ncol(counts)))
    }
    names(groups) <- colnames(counts)

    if(k >  ncol(counts)-1) {
        message("the value of k (", k, ") is too large, setting to ", (ncol(counts)-1))
        k <- ncol(counts)-1
    }

    ls <- estimate.library.sizes(counts, NULL, groups, min.size.entries, verbose = verbose, return.details = TRUE, vil = counts >= min.count.threshold)
    ca <- counts
    ca[ca<min.count.threshold] <- NA # a version of counts with all "drop-out" components set to NA
    mll <- tapply(colnames(counts), groups, function(ids) {
        # use Spearman rank correlation on pairwise complete observations to establish distance relationships between cells
        group <- as.character(groups[ids[1]])

        if(verbose > 0) {
            cat(group, ": calculating cell-cell similarities ...")
        }

        #if(n.cores > 1) { allowWGCNAThreads(n.cores) } else { disableWGCNAThreads() }
        #celld <- WGCNA::cor(log10(matrix(as.numeric(as.matrix(ca)), nrow = nrow(ca), ncol = ncol(ca))+1), method = cor.method, use = "p", nThreads = n.cores)
        if(is.element("WGCNA", installed.packages()[, 1])) {
            celld <- WGCNA::cor(sqrt(matrix(as.numeric(as.matrix(ca[, ids])), nrow = nrow(ca), ncol = length(ids))), method = cor.method, use = "p", nThreads = n.cores)
        } else {
            celld <- stats::cor(sqrt(matrix(as.numeric(as.matrix(ca[, ids])), nrow = nrow(ca), ncol = length(ids))), method = cor.method, use = "p")
        }
        rownames(celld) <- colnames(celld) <- ids

        if(verbose > 0) {
            cat(" done\n")
        }

        # TODO: correct for batch effect in cell-cell similarity matrix
        if(FALSE) {
            # number batches 10^(seq(0, n)) compute matrix of id sums, NA the diagonal,
            bid <- 10^(as.integer(batch)-1)
            bm <- matrix(bid, byrow = TRUE, nrow = length(bid), ncol = length(bid))+bid
            diag(bm) <- NA

            # use tapply to calculate means shifts per combination reconstruct shift vector, matrix, subtract
            # select the upper triangle, tapply to it to correct celld vector directly
        }

        if(verbose)  message(paste("fitting", group, "models:"))

        ml <- papply(seq_along(ids), function(i) { try({
            if(verbose)  message(paste(group, '.', i, " : ", ids[i], sep = ""))
            # determine k closest cells
            oc <- ids[-i][order(celld[ids[i], -i, drop = FALSE], decreasing = TRUE)[1:min(k, length(ids)-1)]]
            #set.seed(i)   oc <- sample(ids[-i], k)
            # determine a subset of genes that show up sufficiently often
            #fpm <- rowMeans(t(t(counts[, oc, drop = FALSE])/(ls$ls[oc])))
            fpm <- apply(t(ca[, oc, drop = FALSE])/(ls$ls[oc]), 2, mean, trim = fpm.estimate.trim, na.rm = TRUE)
            # rank genes by the number of non-zero occurrences, take top genes
            vi <- which(rowSums(counts[, oc] > min.count.threshold)  >=  min(ncol(oc)-1, min.nonfailed) & fpm > min.fpm)
            if(length(vi)<40)  message("WARNING: only ", length(vi), " valid genes were found to fit ", ids[i], " model")
            df <- data.frame(count = counts[vi, ids[i]], fpm = fpm[vi])

            # determine failed-component posteriors for each gene
            #fp <- ifelse(df$count <=  min.count.threshold, threshold.prior, 1-threshold.prior)
            fp <- ifelse(df$count <=  min.count.threshold & df$fpm  >=  median(df$fpm[df$count <=  min.count.threshold]), threshold.prior, 1-threshold.prior)
            cp <- cbind(fp, 1-fp)

            if(linear.fit) {
                # use a linear fit (nb2gth)
                m1 <- fit.nb2gth.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = min.count.threshold, full.theta.range = theta.fit.range, theta.fit.range = theta.fit.range, use.constant.theta.fit = !local.theta.fit, alpha.weight.power = alpha.weight.power)

            }  else {
                # mixture fit (the originally published method)
                m1 <- fit.nb2.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = min.count.threshold)
            }
            v <- get.compressed.v1.model(m1)
            cl <- clusters(m1)
            m1<-list(model = v, clusters = cl)
            #plot.nb2.mixture.fit(m1, df, en = ids[i], do.par = FALSE, compressed.models = TRUE)
            return(m1)
            #})
        })}, n.cores = n.cores)
        vic <- which(unlist(lapply(seq_along(ml), function(i) {
            if(class(ml[[i]]) == "try-error") {
                message("ERROR encountered in building a model for cell ", ids[i], " - skipping the cell. Error:")
                message(ml[[i]])
                #tryCatch(stop(paste("ERROR encountered in building a model for cell ", ids[i])), error = function(e) stop(e))
                return(FALSE);
            }
            return(TRUE);
        })))
        ml <- ml[vic]; names(ml) <- ids[vic];

        if(length(vic)<length(ids)) {
          message("ERROR fitting of ", (length(ids)-length(vic)), " out of ", length(ids), " cells resulted in errors reporting remaining ", length(vic), " cells")
        }
        if(save.model.plots) {
                # model fits
                if(verbose)  message("plotting ", group, " model fits... ")
                tryCatch( {
                    pdf(file = paste(group, "model.fits.pdf", sep = "."), width = ifelse(local.theta.fit, 13, 15), height = 4)
                    l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, ifelse(local.theta.fit, 1, 0.5)), 1), rep(1, 4), FALSE)
                    par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
                    invisible(lapply(vic[1:min(max.model.plots, length(vic))], function(i) {
                        oc <- ids[-i][order(celld[ids[i], -i, drop = FALSE], decreasing = TRUE)[1:min(k, length(ids)-1)]]
                        #set.seed(i) oc <- sample(ids[-i], k)
                        # determine a subset of genes that show up sufficiently often
                        #fpm <- rowMeans(t(t(counts[, oc, drop = FALSE])/(ls$ls[oc])))
                        fpm <- apply(t(ca[, oc, drop = FALSE])/(ls$ls[oc]), 2, mean, trim = fpm.estimate.trim, na.rm = TRUE)
                        vi <- which(rowSums(counts[, oc] > min.count.threshold)  >=  min(ncol(oc)-1, min.nonfailed) & fpm > min.fpm)
                        df <- data.frame(count = counts[vi, ids[i]], fpm = fpm[vi])
                        plot.nb2.mixture.fit(ml[[ids[i]]], df, en = ids[i], do.par = FALSE, compressed.models = TRUE)
                    }))
                    dev.off()
                }, error = function(e) {
                    message("ERROR encountered during model fit plot outputs:")
                    message(e)
                    dev.off()
                })
        }

        return(ml)
    })


    # make a joint model matrix
    jmm <- data.frame(do.call(rbind, lapply(mll, function(tl) do.call(rbind, lapply(tl, function(m) m$model)))))
    rownames(jmm) <- unlist(lapply(mll, names))
    # reorder in the original cell order
    attr(jmm, "groups") <- rep(names(mll), unlist(lapply(mll, length)))
    return(jmm)
}


##' Normalize gene expression variance relative to transcriptome-wide expectations
##'
##' Normalizes gene expression magnitudes to ensure that the variance follows chi-squared statistics
##' with respect to its ratio to the transcriptome-wide expectation as determined by local regression
##' on expression magnitude (and optionally gene length). Corrects for batch effects.
##'
##' @param models model matrix (select a subset of rows to normalize variance within a subset of cells)
##' @param counts read count matrix
##' @param batch measurement batch (optional)
##' @param trim trim value for Winsorization (optional, can be set to 1-3 to reduce the impact of outliers, can be as large as 5 or 10 for datasets with several thousand cells)
##' @param prior expression magnitude prior
##' @param fit.genes a vector of gene names which should be used to establish the variance fit (default is NULL: use all genes). This can be used to specify, for instance, a set spike-in control transcripts such as ERCC.
##' @param plot whether to plot the results
##' @param minimize.underdispersion whether underdispersion should be minimized (can increase sensitivity in datasets with high complexity of population, however cannot be effectively used in datasets where multiple batches are present)
##' @param n.cores number of cores to use
##' @param n.randomizations number of bootstrap sampling rounds to use in estimating average expression magnitude for each gene within the given set of cells
##' @param weight.k k value to use in the final weight matrix
##' @param verbose verbosity level
##' @param weight.df.power power factor to use in determining effective number of degrees of freedom (can be increased for datasets exhibiting particularly high levels of noise at low expression magnitudes)
##' @param smooth.df degrees of freedom to be used in calculating smoothed local regression between coefficient of variation and expression magnitude (and gene length, if provided). Leave at -1 for automated guess.
##' @param max.adj.var maximum value allowed for the estimated adjusted variance (capping of adjusted variance is recommended when scoring pathway overdispersion relative to randomly sampled gene sets)
##' @param theta.range valid theta range (should be the same as was set in knn.error.models() call
##' @param gene.length optional vector of gene lengths (corresponding to the rows of counts matrix)
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' }
##'
##' @return a list containing the following fields:
##' \itemize{
##' \item{mat} {adjusted expression magnitude values}
##' \item{matw} { weight matrix corresponding to the expression matrix}
##' \item{arv} { a vector giving adjusted variance values for each gene}
##' \item{avmodes} {a vector estimated average expression magnitudes for each gene}
##' \item{modes} {a list of batch-specific average expression magnitudes for each gene}
##' \item{prior} {estimated (or supplied) expression magnitude prior}
##' \item{edf} { estimated effective degrees of freedom}
##' \item{fit.genes} { fit.genes parameter }
##' }
##'
##' @export
pagoda.varnorm <- function(models, counts, batch = NULL, trim = 0, prior = NULL, fit.genes=NULL, plot = TRUE, minimize.underdispersion = FALSE, n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9, verbose = 0, weight.df.power = 1, smooth.df = -1, max.adj.var = 10, theta.range = c(1e-2, 1e2), gene.length = NULL) {

    cd <- counts

    min.edf <- 1
    weight.k.internal <- 1
    use.mean.fpm <- FALSE
    use.expected.value <- TRUE
    cv.fit <- TRUE
    edf.damping <- 1

    # load NB extensions
    data(scde.edff, envir = environment())

    # subset cd to the cells occurring in the models
    if(verbose) { cat("checking counts ... ") }
    if(!all(rownames(models) %in% colnames(cd))) {
        stop(paste("supplied count matrix (cd) is missing data for the following cells:[", paste(rownames(models)[!rownames(models) %in% colnames(cd)], collapse = ", "), "]", sep = ""))
    }
    if(!length(rownames(models)) == length(colnames(cd)) || !all(rownames(models) == colnames(cd))) {
        cd <- cd[, match(rownames(models), colnames(cd))]
    }
    if(verbose) { cat("done\n") }

    # trim counts according to the extreme fpm values
    if(trim > 0) {
        if(verbose) { cat("Winsorizing count matrix ... ") }
        fpm <- t((t(log(cd))-models$corr.b)/models$corr.a)
        #tfpm <- log(winsorize.matrix(exp(fpm), trim = trim))
        tfpm <- winsorize.matrix(fpm, trim)
        rn <- rownames(cd)
        cn <- colnames(cd)
        cd <- round(exp(t(t(tfpm)*models$corr.a+models$corr.b)))
        cd[cd<0] <- 0
        rownames(cd) <- rn
        colnames(cd) <- cn
        rm(fpm, tfpm)
        cd <- cd[rowSums(cd) > 0, ] # omit genes without any data after Winsorization
        if(verbose) { cat("done\n") }
    }

    # check/fix batch vector
    if(verbose) { cat("checking batch ... ") }
    if(!is.null(batch)) {
        if(!is.factor(batch)) {
            batch <- as.factor(batch)
        }
        if(is.null(names(batch))) {
            if(length(batch) != nrow(models)) {
                stop("invalid batch vector supplied: length differs from nrow(models)!")
            }
            names(batch) <- rownames(models)
        } else {
            if(!all(rownames(models) %in% names(batch))) {
                stop(paste("invalid batch vector supplied: the following cell(s) are not present: [", paste(rownames(models)[!rownames(models) %in% names(batch)], collapse = ", "), "]", sep = ""))
            }
            batch <- batch[rownames(models)]
        }

        bt <- table(batch)
        min.batch.level <- 2
        if(any(bt<min.batch.level)) {
            if(verbose) { cat("omitting small batch levels [", paste(names(bt)[bt<min.batch.level], collapse = " "), "] ... ") }
            batch[batch %in% names(bt)[bt<min.batch.level]] <- names(bt)[which.max(bt)]
        }
    }
    if(verbose) { cat("ok\n") }

    # recalculate modes as needed
    if(verbose) { cat("calculating modes ... ") }
    if(is.null(prior)) {
        if(verbose) { cat("prior ") }
        prior <- scde.expression.prior(models = models, counts = cd, length.out = 400, show.plot = FALSE)
    }
    # dataset-wide mode
    if(use.mean.fpm) { # use mean fpm across cells
        avmodes <- modes <- rowMeans(exp(scde.expression.magnitude(models, cd)))
    } else { # use joint posterior mode/expected value
        jp <- scde.posteriors(models = models, cd, prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
        if(use.expected.value) {
            avmodes <- modes <- (jp$jp %*% as.numeric(colnames(jp$jp)))[, 1]
        } else { # use mode
            avmodes <- modes <- (as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
        }
    }
    if(verbose) { cat(". ") }

    # batch-specific modes, if necessary
    if(!is.null(batch) && length(levels(batch)) > 1) {
        # calculate mode for each batch
        if(verbose) { cat("batch: [ ") }
        modes <- tapply(seq_len(nrow(models)), batch, function(ii) {
            if(verbose) { cat(as.character(batch[ii[1]]), " ") }
            if(use.mean.fpm) { # use mean fpm across cells
                modes <- rowMeans(exp(scde.expression.magnitude(models[ii, ], cd[, ii])))
            } else { # use joint posterior mode
                jp <- scde.posteriors(models = models[ii, ], cd[, ii], prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
                if(use.expected.value) {
                    modes <- (jp$jp %*% as.numeric(colnames(jp$jp)))[, 1]
                } else { # use mode
                    modes <- (as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
                }
            }
        })
        # set dataset-wide mode
        #if(use.mean.fpm) { # use mean fpm across cells
        #  avmodes <- colMeans(do.call(rbind, modes)*as.vector(unlist(tapply(1:length(batch), batch, length))))*length(levels(batch))/length(batch)
        #jp <- scde.posteriors(models = models, cd, prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
        if(verbose) { cat("] ") }
    }
    if(verbose) { cat("done\n") }

    # check/calculate weights
    if(verbose) { cat("calculating weight matrix ... ") }

    # calculate default weighting scheme
    if(verbose) { cat("calculating ... ") }

    # dataset-wide version of matw (disregarding batch)
    sfp <- do.call(cbind, lapply(seq_len(ncol(cd)), function(i) ppois(cd[, i]-1, exp(models[i, "fail.r"]), lower.tail = FALSE)))
    mfp <- scde.failure.probability(models = models, magnitudes = log(avmodes))
    ofpT <- do.call(cbind, lapply(seq_len(ncol(cd)), function(i) { # for each cell
        lfpm <- log(avmodes)
        mu <- models$corr.b[i] + models$corr.a[i]*lfpm
        thetas <- get.corr.theta(models[i, ], lfpm, theta.range)
        pnbinom(1, size = thetas, mu = exp(mu), lower.tail = TRUE)
    }))
    matw <- 1-weight.k.internal*mfp*sfp # only mode failure probability
    # mode failure or NB failure
    #tmfp <- 1-(1-mfp)*(1-ofpT)
    #matw <- 1-weight.k.internal*tmfp*sfp


    # calculate batch-specific version of the weight matrix if needed
    if(!is.null(batch) && length(levels(batch)) > 1) { # with batch correction
        # save the dataset-wide one as avmatw
        # calculate mode for each batch
        if(verbose) { cat("batch: [ ") }
        bmatw <- do.call(cbind, tapply(seq_len(nrow(models)), batch, function(ii) {
            if(verbose) { cat(as.character(batch[ii[1]]), " ") }
            # set self-fail probability to p(count|background)
            # total mode failure (including overdispersion dropouts)
            #sfp <- do.call(cbind, lapply(ii, function(i) dpois(cd[, i], exp(models[i, "fail.r"]), log = FALSE)))
            sfp <- do.call(cbind, lapply(ii, function(i) ppois(cd[, i]-1, exp(models[i, "fail.r"]), lower.tail = FALSE)))

            mfp <- scde.failure.probability(models = models[ii, ], magnitudes = log(modes[[batch[ii[1]]]]))
            ofpT <- do.call(cbind, lapply(ii, function(i) { # for each cell
                lfpm <- log(modes[[batch[i]]])
                mu <- models$corr.b[i] + models$corr.a[i]*lfpm
                thetas <- get.corr.theta(models[i, ], lfpm, theta.range)
                pnbinom(1, size = thetas, mu = exp(mu), lower.tail = TRUE)
            }))

            x <- 1-weight.k.internal*mfp*sfp # only mode failure probability
            # mode failure or NB failure
            #tmfp <- 1-(1-mfp)*(1-ofpT)
            #x <- 1-weight.k.internal*tmfp*sfp
        }))
        # reorder
        bmatw <- bmatw[, rownames(models)]
        if(verbose) { cat("] ") }
    }
    if(verbose) { cat("done\n") }

    # calculate effective degrees of freedom
    # total effective degrees of freedom per gene
    if(verbose) { cat("calculating effective degrees of freedom ..") }
    ids <- 1:ncol(cd)
    names(ids) <- colnames(cd)
    # dataset-wide version
    edf.mat <- do.call(cbind, papply(ids, function(i) {
        v <- models[i, ]
        lfpm <- log(avmodes)
        mu <- exp(lfpm*v$corr.a + v$corr.b)
        # adjust very low mu levels except for those that have 0 counts (to avoid inf values)

        thetas <- get.corr.theta(v, lfpm, theta.range)
        edf <- exp(predict(scde.edff, data.frame(lt = log(thetas))))
        edf[thetas > 1e3] <- 1
        edf
    }, n.cores = n.cores))
    if(edf.damping != 1) {
        edf.mat <- ((edf.mat/ncol(edf.mat))^edf.damping) * ncol(edf.mat)
    }

    # incorporate weight into edf
    #edf.mat <- ((matw^weight.df.power)*edf.mat)
    edf.mat <- (matw*edf.mat)^weight.df.power
    #edf <- rowSums(matw*edf.mat)+1.5 # summarize eDF per gene
    edf <- rowSums(edf.mat)+1 # summarize eDF per gene
    if(verbose) { cat(".") }

    # batch-specific version if necessary
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        bedf.mat <- do.call(cbind, papply(ids, function(i) {
            v <- models[i, ]
            lfpm <- log(modes[[batch[i]]])
            mu <- exp(lfpm*v$corr.a + v$corr.b)
            # adjust very low mu levels except for those that have 0 counts (to avoid inf values)

            thetas <- get.corr.theta(v, lfpm, theta.range)
            edf <- exp(predict(scde.edff, data.frame(lt = log(thetas))))
            edf[thetas > 1e3] <- 1
            return(edf)
        }, n.cores = n.cores))
        if(edf.damping != 1) { bedf.mat <-  ((bedf.mat/ncol(bedf.mat))^edf.damping) * ncol(edf.mat) }

        # incorporate weight into edf
        #bedf.mat <- ((bmatw^weight.df.power)*bedf.mat)
        bedf.mat <- (bmatw*bedf.mat)^weight.df.power
        bedf <- rowSums(bedf.mat)+1 # summarize eDF per gene
        if(verbose) { cat(".") }
    }

    if(verbose) { cat(" done\n") }

    if(verbose) { cat("calculating normalized expression values ... ") }
    # evaluate negative binomial deviations and effective degrees of freedom
    ids <- 1:ncol(cd)
    names(ids) <- colnames(cd)
    mat <- do.call(cbind, papply(ids, function(i) {
        v <- models[i, ]
        lfpm <- log(avmodes)
        mu <- exp(lfpm*v$corr.a + v$corr.b)
        # adjust very low mu levels except for those that have 0 counts (to avoid inf values)
        thetas <- get.corr.theta(v, lfpm, theta.range)

        #matw[, i]*edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
        #x <- (cd[, i]-mu)^2/(mu+mu^2/thetas)
        #edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
        # considering Poisson-nb mixture
        fail.lambda <- exp(as.numeric(v["fail.r"]))
        #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*(mu+mu^2/thetas) + (1-matw[, i])*((mu-fail.lambda)^2 + fail.lambda))
        edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas +  fail.lambda)

        #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*mu+(mu^2)*((1-matw[, i])+matw[, i]/thetas))
    }, n.cores = n.cores))
    rownames(mat) <- rownames(cd)
    if(verbose) { cat(".") }
    # batch-specific version of mat
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        bmat <- do.call(cbind, papply(ids, function(i) {
            v <- models[i, ]
            lfpm <- log(modes[[batch[i]]])
            mu <- exp(lfpm*v$corr.a + v$corr.b)
            # adjust very low mu levels except for those that have 0 counts (to avoid inf values)
            thetas <- get.corr.theta(v, lfpm, theta.range)

            #matw[, i]*edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
            #x <- (cd[, i]-mu)^2/(mu+mu^2/thetas)
            #edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
            #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*mu+(mu^2)*((1-matw[, i])+matw[, i]/thetas))
            fail.lambda <- exp(as.numeric(v["fail.r"]))
            #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*(mu+mu^2/thetas) + (1-matw[, i])*((mu-fail.lambda)^2 + fail.lambda))
            edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas +  fail.lambda)
        }, n.cores = n.cores))
        rownames(bmat) <- rownames(cd)

        if(verbose) { cat(".") }
    }
    if(verbose) { cat(" done\n") }

    # do a model fit on the weighted standard deviation (as a function of the batch-average expression mode)
    wvar <- rowSums(mat)/rowSums(edf.mat)

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        # estimate the ratio of the batch-specific variance to the total dataset variance
        bwvar <- rowSums(bmat)/rowSums(bedf.mat)
        bwvar.ratio <- bwvar/wvar
        wvar <- bwvar # replace wvar now that we have the ratio of
        matw <- bmatw # replace matw with the batch-specific one that will be used from here on
        # ALTERNATIVE: could adjust wvar for the bwvar.ratio here, before fitting expression dependency
      }
    fvi <- vi <- rowSums(matw) > 0 & is.finite(wvar) & wvar > 0
    if(!is.null(fit.genes)) { fvi <- fvi & rownames(mat) %in% fit.genes }
    if(!any(fvi)) { stop("unable to find a set of valid genes to establish the variance fit") }

    # s = mgcv:::s
    s = mgcv::s
    if(cv.fit) {
        #x <- gam(as.formula("cv2 ~ s(lev)"), data = df[vi, ], weights = rowSums(matw[vi, ]))
        if(is.null(gene.length)) {
            df <- data.frame(lev = log10(avmodes), cv2 = log10(wvar/avmodes^2))
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        } else {
            df <- data.frame(lev = log10(avmodes), cv2 = log10(wvar/avmodes^2), len = gene.length[rownames(cd)])
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df) + s(len, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        }
        #x <- lm(cv2~lev, data = df[vi, ], weights = rowSums(matw[vi, ]))

        zval.m <- 10^(df$cv2[vi]-predict(x, newdata = df[vi, ]))

        if(plot) {
            par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
            #smoothScatter(df$lev[vi], log(wvar[vi]), nbin = 256, xlab = "expression magnitude (log10)", ylab = "wvar (log)") abline(h = 0, lty = 2, col = 2)
            #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], log(wvar[paste("g", diff.exp.gene.ids, sep = "")]), col = 2)

            smoothScatter(df$lev[vi], df$cv2[vi], nbin = 256, xlab = "expression magnitude (log10)", ylab = "cv^2 (log10)")
            lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 2, pch = ".", cex = 1)
            if(!is.null(fit.genes)) { # show genes used for the fit
              points(df$lev[fvi],df$cv2[fvi],pch=".",col="green",cex=1)
            }

            #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], df[paste("g", diff.exp.gene.ids, sep = ""), "cv2"], col = 2)
        }

        # optional : re-weight to minimize the underdispersed points
        if(minimize.underdispersion) {
            pv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = FALSE, lower.tail = TRUE)
            pv[edf[vi]<= min.edf] <- 0
            pv <- p.adjust(pv)
            #x <- gam(as.formula("cv2 ~ s(lev)"), data = df[vi, ], weights = (pmin(10, -log(pv))+1)*rowSums(matw[vi, ]))
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df), data = df[fvi, ], weights = (pmin(10, -log(pv))+1)*rowSums(matw[fvi, ]))
            zval.m <- 10^(df$cv2[vi]-predict(x,newdata=df[vi,]))
            if(plot) {
              lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 4, pch = ".", cex = 1)
            }
        }
    } else {
        df <- data.frame(lev = log10(avmodes), sd = sqrt(wvar))
        #x <- gam(as.formula("sd ~ s(lev)"), data = df[vi, ], weights = rowSums(matw[vi, ]))
        x <- mgcv::gam(sd ~ s(lev, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        zval.m <- (as.numeric((df$sd[vi])/pmax(min.sd, predict(x,newdata=df[vi,]))))^2

        if(plot) {
            par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
            smoothScatter(df$lev[vi], df$sd[vi], nbin = 256, xlab = "expression magnitude", ylab = "weighted sdiv")
            lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 2, pch = ".", cex = 1)
            if(!is.null(fit.genes)) { # show genes used for the fit
              points(df$lev[fvi],df$sd[fvi],pch=".",col="green",cex=1)
            }
        }

        # optional : re-weight to minimize the underdispersed points
        if(minimize.underdispersion) {
            pv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = FALSE, lower.tail = TRUE)
            pv[edf[vi]<= min.edf] <- 0
            pv <- p.adjust(pv)
            #x <- gam(as.formula("sd ~ s(lev)"), data = df[vi, ], weights = (pmin(20, -log(pv))+1)*rowSums(matw[vi, ]))
            x <- mgcv::gam(sd ~ s(lev, k = smooth.df), data = df[fvi, ], weights = (pmin(20, -log(pv))+1)*rowSums(matw[fvi, ]))
            zval.m <- (as.numeric((df$sd[vi])/pmax(x$fitted.values, min.sd)))^2
            if(plot) {
              lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 4, pch = ".", cex = 1)
            }
        }
    }

    # adjust for inter-batch variance
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        #zval.m <- zval.m*pmin(bwvar.ratio[vi], 1) # don't increase zval.m even if batch-specific specific variance is higher than the dataset-wide variance
        zval.m <- zval.m*pmin(bwvar.ratio[vi], 1/bwvar.ratio[vi]) # penalize for strong deviation in either direction
    }

    # calculate adjusted variance
    qv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = TRUE, lower.tail = FALSE)
    qv[edf[vi]<= min.edf] <- 0
    qv[abs(qv)<1e-10] <- 0
    arv <- rep(NA, length(vi))
    arv[vi] <- qchisq(qv, ncol(matw)-1, lower.tail = FALSE, log.p = TRUE)/ncol(matw)
    arv <- pmin(max.adj.var, arv)
    names(arv) <- rownames(cd)
    if(plot) {
        smoothScatter(df$lev[vi], arv[vi], xlab = "expression magnitude (log10)", ylab = "adjusted variance (log10)", nbin = 256)
        abline(h = 1, lty = 2, col = 8)
        abline(h = max.adj.var, lty = 3, col = 2)
        if(!is.null(fit.genes)) {
          points(df$lev[fvi],arv[fvi],pch=".",col="green",cex=1)
        }
        #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], arv[paste("g", diff.exp.gene.ids, sep = "")], col = 2)
        #points(df$lev[vi], arv[vi], col = 2, pch = ".", cex = 2)
    }

    # Wilcox score upper bound
    wsu <- function(k, n, z = qnorm(0.975)) {
        p <- k/n
        pmin(1, (2*n*p+z^2+(z*sqrt(z^2-1/n+4*n*p*(1-p)-(4*p-2)) +1))/(2*(n+z^2)))
    }

    # use milder weight matrix
    #matw <- 1-0.9*((1-matw)^2) # milder weighting for the the PCA (1-0.9*sp*mf)
    matw <- 1-weight.k*(1-matw) # milder weighting for the the PCA (1-0.9*sp*mf)
    matw <- matw/rowSums(matw)
    mat <- log10(exp(scde.expression.magnitude(models, cd))+1)

    # estimate observed variance (for scaling) before batch adjustments
    #varm <- sqrt(arv/pmax(weightedMatVar(mat, matw, batch = batch), 1e-5)) varm[varm<1e-5] <- 1e-5 mat <- mat*varm
    ov <- weightedMatVar(mat, matw)
    vr <- arv/ov
    vr[ov <=  0] <- 0

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        # adjust proportion of zeros
        # determine lowest upper bound of non-zero measurement probability among the batch (for each gene)
        nbub <- apply(do.call(cbind, tapply(seq_len(ncol(mat)), batch, function(ii) {
            wsu(rowSums(mat[, ii] > 0), length(ii), z = qnorm(1-1e-2))
        })), 1, min)

        # decrease the batch weights for each gene to match the total
        # expectation of the non-zero measurements
        nbo <- do.call(cbind, tapply(seq_len(ncol(mat)), batch, function(ii) {
            matw[, ii]*pmin(1, ceiling(nbub*length(ii))/rowSums(mat[, ii] > 0))
        }))
        nbo <- nbo[, colnames(matw)]
        matw <- nbo

        ## # center 0 and non-0 observations between batches separately
        ## amat <- mat amat[amat == 0] <- NA
        ## amat.av <- rowMeans(amat, na.rm = TRUE) # dataset means
        ## # adjust each batch by the mean of its non-0 measurements
        ## amat <- do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   amat[, ii]-rowMeans(amat[, ii], na.rm = TRUE)
        ## }))
        ## amat <- amat[, colnames(mat)] # fix the ordering
        ## # shift up each gene by the dataset mean
        ## amat <- amat+amat.av
        ## amat[is.na(amat)] <- 0
        ## mat <- amat

        amat <- mat
        nr <- ncol(matw)/rowSums(matw)
        amat.av <- rowMeans(amat*matw)*nr # dataset means
        amat <- do.call(cbind, tapply(seq_len(ncol(amat)), batch, function(ii) {
            amat[, ii]-(rowMeans(amat[, ii]*matw[, ii]*nr, na.rm = TRUE))
        }))
        amat <- amat[, colnames(matw)]
        mat <- amat+amat.av

        # alternative: actually zero-out entries in mat
        ## nbub <- rowMin(do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   wsu(rowSums(amat[, ii] > 0), length(ii), z = qnorm(1-1e-2))
        ## })))
        ## set.seed(0)

        ## # decrease the batch weights for each gene to match the total
        ## # expectation of the non-zero measurements
        ## matm <- do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   # number of entries to zero-out per gene
        ##   nze <- rowSums(amat[, ii] > 0) - ceiling(nbub*length(ii))
        ##   # construct mat multiplier submatrix
        ##   sa <- rep(1, length(ii))
        ##   smatm <- do.call(rbind, lapply(1:length(nze), function(ri) {
        ##     if(nze[ri]<1) { return(sa) }
        ##     vi <- which(mat[ri, ii] > 0)
        ##     a <- sa a[vi[sample.int(length(vi), nze[ri])]] <- 0
        ##     a
        ##   }))
        ##   colnames(smatm) <- colnames(mat[, ii])
        ##   rownames(smatm) <- rownames(mat)
        ##   smatm
        ## }))
        ## matm <- matm[, colnames(mat)]
        ## mat <- mat*matm
        ## matw <- matw*matm
    }

    # center (no batch)
    mat <- weightedMatCenter(mat, matw)
    mat <- mat*sqrt(vr)

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        return(list(mat = mat, matw = matw, arv = arv, modes = modes, avmodes = avmodes, prior = prior, edf = edf, batch = batch, trim = trim, bwvar.ratio = bwvar.ratio))
    } else {
        return(list(mat = mat, matw = matw, arv = arv, modes = modes, avmodes = avmodes, prior = prior, edf = edf, batch = batch, trim = trim))
    }
}


##' Control for a particular aspect of expression heterogeneity in a given population
##'
##' Similar to subtracting n-th principal component, the current procedure determines
##' (weighted) projection of the expression matrix onto a specified aspect (some pattern
##' across cells, for instance sequencing depth, or PC corresponding to an undesired process
##' such as ribosomal pathway variation) and subtracts it from the data so that it is controlled
##' for in the subsequent weighted PCA analysis.
##'
##' @param varinfo normalized variance info (from pagoda.varnorm())
##' @param aspect a vector giving a cell-to-cell variation pattern that should be controlled for (length should be corresponding to ncol(varinfo$mat))
##' @param center whether the matrix should be re-centered following pattern subtraction
##'
##' @return a modified varinfo object with adjusted expression matrix (varinfo$mat)
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' # create go environment
##' library(org.Hs.eg.db)
##' # translate gene names to ids
##' ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
##' rids <- names(ids); names(rids) <- ids
##' go.env <- lapply(mget(ls(org.Hs.egGO2ALLEGS), org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
##' # clean GOs
##' go.env <- clean.gos(go.env)
##' # convert to an environment
##' go.env <- list2env(go.env)
##' # subtract the pattern
##' cc.pattern <- pagoda.show.pathways(ls(go.env)[1:2], varinfo, go.env, show.cell.dendrogram = TRUE, showRowLabels = TRUE)  # Look at pattern from 2 GO annotations
##' varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)
##' }
##'
##' @export
pagoda.subtract.aspect <- function(varinfo, aspect, center = TRUE) {
    if(length(aspect) != ncol(varinfo$mat)) { stop("aspect should be a numeric vector of the same length as the number of cells (i.e. ncol(varinfo$mat))") }
    v <- aspect
    v <- v-mean(v)
    v <- v/sqrt(sum(v^2))
    nr <- ((varinfo$mat * varinfo$matw) %*% v)/(varinfo$matw %*% v^2)
    mat.c <- varinfo$mat - t(v %*% t(nr))
    if(center) {
        mat.c <- weightedMatCenter(mat.c, varinfo$matw) # this commonly re-introduces some background dependency because of the matw
    }
    varinfo$mat <- mat.c
    varinfo
}


##' Run weighted PCA analysis on pre-annotated gene sets
##'
##' For each valid gene set (having appropriate number of genes) in the provided environment (setenv),
##' the method will run weighted PCA analysis, along with analogous analyses of random gene sets of the
##' same size, or shuffled expression magnitudes for the same gene set.
##'
##' @param varinfo adjusted variance info from pagoda.varinfo() (or pagoda.subtract.aspect())
##' @param setenv environment listing gene sets (contains variables with names corresponding to gene set name, and values being vectors of gene names within each gene set)
##' @param n.components number of principal components to determine for each gene set
##' @param n.cores number of cores to use
##' @param min.pathway.size minimum number of observed genes that should be contained in a valid gene set
##' @param max.pathway.size maximum number of observed genes in a valid gene set
##' @param n.randomizations number of random gene sets (of the same size) to be evaluated in parallel with each gene set (can be kept at 5 or 10, but should be increased to 50-100 if the significance of pathway overdispersion will be determined relative to random gene set models)
##' @param n.internal.shuffles number of internal (independent row shuffles) randomizations of expression data that should be evaluated for each gene set (needed only if one is interested in gene set coherence P values, disabled by default; set to 10-30 to estimate)
##' @param n.starts number of random starts for the EM method in each evaluation
##' @param center whether the expression matrix should be recentered
##' @param batch.center whether batch-specific centering should be used
##' @param proper.gene.names alternative vector of gene names (replacing rownames(varinfo$mat)) to be used in cases when the provided setenv uses different gene names
##' @param verbose verbosity level
##'
##' @return a list of weighted PCA info for each valid gene set
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' # create go environment
##' library(org.Hs.eg.db)
##' # translate gene names to ids
##' ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
##' rids <- names(ids); names(rids) <- ids
##' go.env <- lapply(mget(ls(org.Hs.egGO2ALLEGS), org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
##' # clean GOs
##' go.env <- clean.gos(go.env)
##' # convert to an environment
##' go.env <- list2env(go.env)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' }
##'
##' @export
pagoda.pathway.wPCA <- function(varinfo, setenv, n.components = 2, n.cores = detectCores(), min.pathway.size = 10, max.pathway.size = 1e3, n.randomizations = 10, n.internal.shuffles = 0, n.starts = 10, center = TRUE, batch.center = TRUE, proper.gene.names = NULL, verbose = 0) {
    mat <- varinfo$mat
    matw <- varinfo$matw
    gsl <- NULL
    return.gsl <- FALSE
    smooth <- 0
    if(batch.center) { batch <- varinfo$batch } else { batch <- NULL }

    if(is.null(proper.gene.names)) { proper.gene.names <- rownames(mat) }


    if(center) {
        mat <- weightedMatCenter(mat, matw, batch = batch)
    }

    vi <- apply(mat, 1, function(x) sum(abs(diff(x))) > 0)
    vi[is.na(vi)] <- FALSE
    mat <- mat[vi, , drop = FALSE] # remove constant rows
    matw <- matw[vi, , drop = FALSE]
    proper.gene.names <- proper.gene.names[vi]

    if(is.null(gsl)) {
        gsl <- ls(envir = setenv)
        gsl.ng <- unlist(lapply(sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names)))
        gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
        names(gsl) <- gsl
    }
    if(verbose) {
        message("processing ", length(gsl), " valid pathways")
    }
    if(return.gsl) return(gsl)


    # transpose mat to save a bit of calculations
    mat <- t(mat)
    matw <- t(matw)

    mcm.pc <- papply(gsl, function(x) {
        lab <- proper.gene.names %in% get(x, envir = setenv)
        if(sum(lab)<1) { return(NULL) }

        #smooth <- round(sum(lab)*smooth.fraction)
        #smooth <- max(sum(lab), smooth)

        #xp <- pca(d, nPcs = n.components, center = TRUE, scale = "none")
        #xp <- epca(mat[, lab], ncomp = n.components, center = FALSE, nstarts = n.starts)
        xp <- bwpca(mat[, lab, drop = FALSE], matw[, lab, drop = FALSE], npcs = n.components, center = FALSE, nstarts = n.starts, smooth = smooth, n.shuffles = n.internal.shuffles)

        # get standard deviations for the random samples
        ngenes <- sum(lab)
        z <- do.call(rbind, lapply(seq_len(n.randomizations), function(i) {
            si <- sample(1:ncol(mat), ngenes)
            #epca(mat[, si], ncomp = 1, center = FALSE, nstarts = n.starts)$sd
            xp <- bwpca(mat[, si, drop = FALSE], matw[, si, drop = FALSE], npcs = 1, center = FALSE, nstarts = n.starts, smooth = smooth)$sd
        }))

        # flip orientations to roughly correspond with the means
        cs <- unlist(lapply(seq_len(ncol(xp$scores)), function(i) sign(cor(xp$scores[, i], colMeans(t(mat[, lab, drop = FALSE])*abs(xp$rotation[, i]))))))

        xp$scores <- t(t(xp$scores)*cs)
        xp$rotation <- t(t(xp$rotation)*cs)

        # local normalization of each component relative to sampled PC1 sd
        avar <- pmax(0, (xp$sd^2-mean(z[, 1]^2))/sd(z[, 1]^2))
        xv <- t(xp$scores)
        xv <- xv/apply(xv, 1, sd)*sqrt(avar)
        return(list(xv = xv, xp = xp, z = z, sd = xp$sd, n = ngenes))
    }, n.cores = n.cores)
}


##' Estimate effective number of cells based on lambda1 of random gene sets
##'
##' Examines the dependency between the amount of variance explained by the first principal component
##' of a gene set and the number of genes in a gene set to determine the effective number of cells
##' for the Tracy-Widom distribution
##'
##' @param pwpca result of the pagoda.pathway.wPCA() call with n.randomizations > 1
##' @param start optional starting value for the optimization (if the NLS breaks, trying high starting values usually fixed the local gradient problem)
##'
##' @return effective number of cells
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' pagoda.effective.cells(pwpca)
##' }
##'
##' @export
pagoda.effective.cells <- function(pwpca, start = NULL) {
    n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
    var <- unlist(lapply(pwpca, function(x) x$z[, 1]))^2
    if(is.null(start)) { start <- nrow(pwpca[[1]]$xp$scores)*10 }

    n.cells <- nrow(pwpca[[1]]$xp$scores)
    of <- function(p, v, sp) {
        sn <- p[1]
        vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
        residuals <- (v-vfit)^2
        return(sum(residuals))
    }
    x <- nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(n.cells))
    return((x$par)^2+1/2)
}


##' Determine de-novo gene clusters and associated overdispersion info
##'
##' Determine de-novo gene clusters, their weighted PCA lambda1 values, and random matrix expectation.
##'
##' @param varinfo varinfo adjusted variance info from pagoda.varinfo() (or pagoda.subtract.aspect())
##' @param trim additional Winsorization trim value to be used in determining clusters (to remove clusters that group outliers occurring in a given cell). Use higher values (5-15) if the resulting clusters group outlier patterns
##' @param n.clusters number of clusters to be determined (recommended range is 100-200)
##' @param cor.method correlation method ("pearson", "spearman") to be used as a distance measure for clustering
##' @param n.samples number of randomly generated matrix samples to test the background distribution of lambda1 on
##' @param n.starts number of wPCA EM algorithm starts at each iteration
##' @param n.internal.shuffles number of internal shuffles to perform (only if interested in set coherence, which is quite high for clusters by definition, disabled by default; set to 10-30 shuffles to estimate)
##' @param n.cores number of cores to use
##' @param verbose verbosity level
##' @param plot whether a plot showing distribution of random lambda1 values should be shown (along with the extreme value distribution fit)
##' @param show.random whether the empirical random gene set values should be shown in addition to the Tracy-Widom analytical approximation
##' @param n.components number of PC to calculate (can be increased if the number of clusters is small and some contain strong secondary patterns - rarely the case)
##' @param method clustering method to be used in determining gene clusters
##' @param secondary.correlation whether clustering should be performed on the correlation of the correlation matrix instead
##' @param n.cells number of cells to use for the randomly generated cluster lambda1 model
##' @param old.results optionally, pass old results just to plot the model without recalculating the stats
##'
##' @return a list containing the following fields:
##' \itemize{
##' \item{clusters} {a list of genes in each cluster values}
##' \item{xf} { extreme value distribution fit for the standardized lambda1 of a randomly generated pattern}
##' \item{tci} { index of a top cluster in each random iteration}
##' \item{cl.goc} {weighted PCA info for each real gene cluster}
##' \item{varm} {standardized lambda1 values for each randomly generated matrix cluster}
##' \item{clvlm} {a linear model describing dependency of the cluster lambda1 on a Tracy-Widom lambda1 expectation}
##' }
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' clpca <- pagoda.gene.clusters(varinfo, trim=7.1/ncol(varinfo$mat), n.clusters=150, n.cores=10, plot=FALSE)
##' }
##'
##' @export
pagoda.gene.clusters <- function(varinfo, trim = 3.1/ncol(varinfo$mat), n.clusters = 150, n.samples = 60, cor.method = "p", n.internal.shuffles = 0, n.starts = 10, n.cores = detectCores(), verbose = 0, plot = FALSE, show.random = FALSE, n.components = 1, method = "ward.D", secondary.correlation = FALSE, n.cells = ncol(varinfo$mat), old.results = NULL) {

    smooth <- 0
    mat <- varinfo$mat
    matw <- varinfo$matw
    batch = varinfo$batch

    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    if(!is.null(batch)) {
        # center mat by batch
        mat <- weightedMatCenter(mat, matw, batch)
    }


    if(!is.null(old.results)) {
        if(verbose) { cat ("reusing old results for the observed clusters\n")}
        gcls <- old.results$clusters
        cl.goc <- old.results$cl.goc
    } else {
        if(verbose) { cat ("determining gene clusters ...")}
        # actual clusters
        vi<-which(abs(apply(mat, 1, function(x) sum(abs(diff(x))))) > 0)
        if(is.element("WGCNA", installed.packages()[, 1])) {
            gd <- as.dist(1-WGCNA::cor(t(mat)[, vi], method = cor.method, nThreads = n.cores))
        } else {
            gd <- as.dist(1-cor(t(mat)[, vi], method = cor.method))
        }

        if(secondary.correlation) {
            if(is.element("WGCNA", installed.packages()[, 1])) {
                gd <- as.dist(1-WGCNA::cor(as.matrix(gd), method = "p", nThreads = n.cores))
            } else {
                gd <- as.dist(1-cor(as.matrix(gd), method = "p"))
            }
        }

        if(is.element("fastcluster", installed.packages()[, 1])) {
            gcl <- fastcluster::hclust(gd, method = method)
        } else {
            gcl <- stats::hclust(gd, method = method)
        }
        gcll <- cutree(gcl, n.clusters)
        gcls <- tapply(rownames(mat)[vi], as.factor(gcll), I)
        names(gcls) <- paste("geneCluster", names(gcls), sep = ".")

        rm(gd, gcl)
        gc()

        # determine PC1 for the actual clusters
        if(verbose) { cat (" cluster PCA ...")}
        il <- tapply(vi, factor(gcll, levels = c(1:length(gcls))), I)
        cl.goc <- papply(il, function(ii) {
            xp <- bwpca(t(mat[ii, , drop = FALSE]), t(matw[ii, , drop = FALSE]), npcs = n.components, center = FALSE, nstarts = n.starts, smooth = smooth, n.shuffles = n.internal.shuffles)

            cs <- unlist(lapply(seq_len(ncol(xp$scores)), function(i) sign(cor(xp$scores[, i], colMeans(mat[ii, , drop = FALSE]*abs(xp$rotation[, i]))))))

            xp$scores <- t(t(xp$scores)*cs)
            xp$rotation <- t(t(xp$rotation)*cs)

            return(list(xp = xp, sd = xp$sd, n = length(ii)))
        }, n.cores = n.cores)
        names(cl.goc) <- paste("geneCluster", names(cl.goc), sep = ".")

        if(verbose) { cat ("done\n")}
    }

    # sampled variation
    if(!is.null(old.results) && !is.null(old.results$varm)) {
        if(verbose) { cat ("reusing old results for the sampled clusters\n")}
        varm <- old.results$varm } else {
            if(verbose) { cat ("generating", n.samples, "randomized samples ")}
            varm <- do.call(rbind, papply(seq_len(n.samples), function(i) { # each sampling iteration
                set.seed(i)
                # generate random normal matrix
                # TODO: use n.cells instead of ncol(matw)
                m <- matrix(rnorm(nrow(mat)*n.cells), nrow = nrow(mat), ncol = n.cells)
                #m <- weightedMatCenter(m, matw, batch = batch)

                if(show.random) {
                    full.m <- t(m) # save untrimmed version of m for random gene set controls
                }

                if(trim > 0) {
                    m <- winsorize.matrix(m, trim = trim)
                }

                vi<-which(abs(apply(m, 1, function(x) sum(diff(abs(x))))) > 0)
                if(is.element("WGCNA", installed.packages()[, 1])) {
                    gd <- as.dist(1-WGCNA::cor(t(m[vi, ]), method = cor.method, nThreads = 1))
                } else {
                    gd <- as.dist(1-cor(t(m[vi, ]), method = cor.method))
                }
                if(secondary.correlation) {
                    if(is.element("WGCNA", installed.packages()[, 1])) {
                        gd <- as.dist(1-WGCNA::cor(as.matrix(gd), method = "p", nThreads = 1))
                    } else {
                        gd <- as.dist(1-cor(as.matrix(gd), method = "p"))
                    }
                }

                if(is.element("fastcluster", installed.packages()[, 1])) {
                    gcl <- fastcluster::hclust(gd, method = method)
                } else {
                    gcl <- stats::hclust(gd, method = method)
                }
                gcll <- cutree(gcl, n.clusters)
                rm(gd, gcl)
                gc()

                # transpose to save time
                m <- t(m) # matw <- t(matw)

                sdv <- tapply(vi, gcll, function(ii) {
                    #as.numeric(bwpca(m[, ii], matw[, ii], npcs = 1, center = FALSE, nstarts = n.starts, smooth = smooth)$sd)^2
                    pcaMethods::sDev(pcaMethods::pca(m[, ii], nPcs = 1, center = FALSE))^2
                })

                pathsizes <- unlist(tapply(vi, gcll, length))
                names(pathsizes) <- pathsizes

                if(show.random) {
                    rsdv <- unlist(lapply(names(pathsizes), function(s) {
                        vi <- sample(1:ncol(full.m), as.integer(s))
                        pcaMethods::sDev(pcaMethods::pca(full.m[, vi], nPcs = 1, center = FALSE))^2
                    }))
                    if(verbose) { cat (".")}
                    return(data.frame(n = as.integer(pathsizes), var = unlist(sdv), round = i, rvar = rsdv))
                }

                if(verbose) { cat (".")}
                data.frame(n = as.integer(pathsizes), var = unlist(sdv), round = i)

            }, n.cores = n.cores))
            if(verbose) { cat ("done\n")}
        }

    # score relative to Tracey-Widom distribution
    #require(RMTstat)
    x <- RMTstat::WishartMaxPar(n.cells, varm$n)
    varm$pm <- x$centering-(1.2065335745820)*x$scaling # predicted mean of a random set
    varm$pv <- (1.607781034581)*x$scaling # predicted variance of a random set
    #clvlm <- lm(var~pm, data = varm)
    clvlm <- lm(var~0+pm+n, data = varm)
    varm$varst <- (varm$var-predict(clvlm))/sqrt(varm$pv)
    #varm$varst <- as.numeric(varm$var - (cbind(1, varm$pm) %*% coef(clvlm)))/sqrt(varm$pv)
    #varm$varst <- as.numeric(varm$var - (varm$pm* coef(clvlm)[2]))/sqrt(varm$pv)

    #varm$varst <- (varm$var-varm$pm)/sqrt(varm$pv)
    tci <- tapply(seq_len(nrow(varm)), as.factor(varm$round), function(ii) ii[which.max(varm$varst[ii])])

    #xf <- fevd(varm$varst[tci], type = "Gumbel") # fit on top clusters
    xf <- extRemes::fevd(varm$varst, type = "Gumbel") # fit on all clusters

    if(plot) {
        require(extRemes)
        par(mfrow = c(1, 2), mar = c(3.5, 3.5, 3.5, 1.0), mgp = c(2, 0.65, 0), cex = 0.9)
        smoothScatter(varm$n, varm$var, main = "simulations", xlab = "cluster size", ylab = "PC1 variance")
        if(show.random) {
            points(varm$n, varm$rvar, pch = ".", col = "red")
        }
        #pv <- predict(rsm, newdata = data.frame(n = sort(varm$n)), se.fit = TRUE)
        on <- order(varm$n, decreasing = TRUE)
        lines(varm$n[on], predict(clvlm)[on], col = 4, lty = 3)
        lines(varm$n[on], varm$pm[on], col = 2)
        lines(varm$n[on], (varm$pm+1.96*sqrt(varm$pv))[on], col = 2, lty = 2)
        lines(varm$n[on], (varm$pm-1.96*sqrt(varm$pv))[on], col = 2, lty = 2)
        legend(x = "bottomright", pch = c(1, 19, 19), col = c(1, 4, 2), legend = c("top clusters", "clusters", "random"), bty = "n")

        points(varm$n[tci], varm$var[tci], col = 1)
        extRemes::plot.fevd(xf, type = "density", main = "Gumbel fit")
        abline(v = 0, lty = 3, col = 4)
    }

    #pevd(9, loc = xf$results$par[1], scale = xf$results$par[2], lower.tail = FALSE)
    #xf$results$par

    return(list(clusters = gcls, xf = xf, tci = tci, cl.goc = cl.goc, varm = varm, clvlm = clvlm, trim = trim))
}


##' Score statistical significance of gene set and cluster overdispersion
##'
##' Evaluates statistical significance of the gene set and cluster lambda1 values, returning
##' either a text table of Z scores, etc, a structure containing normalized values of significant
##' aspects, or a set of genes underlying the significant aspects.
##'
##' @param pwpca output of pagoda.pathway.wPCA()
##' @param clpca output of pagoda.gene.clusters() (optional)
##' @param n.cells effective number of cells (if not provided, will be determined using pagoda.effective.cells())
##' @param z.score Z score to be used as a cutoff for statistically significant patterns (defaults to 0.05 P-value
##' @param return.table whether a text table showing
##' @param return.genes whether a set of genes driving significant aspects should be returned
##' @param plot whether to plot the cv/n vs. dataset size scatter showing significance models
##' @param adjust.scores whether the normalization of the aspect patterns should be based on the adjusted Z scores - qnorm(0.05/2, lower.tail = FALSE)
##' @param score.alpha significance level of the confidence interval for determining upper/lower bounds
##' @param use.oe.scale whether the variance of the returned aspect patterns should be normalized using observed/expected value instead of the default chi-squared derived variance corresponding to overdispersion Z score
##' @param effective.cells.start starting value for the pagoda.effective.cells() call
##'
##' @return if return.table = FALSE and return.genes = FALSE (default) returns a list structure containing the following items:
##' \itemize{
##' \item{xv} {a matrix of normalized aspect patterns (rows- significant aspects, columns- cells}
##' \item{xvw} { corresponding weight matrix }
##' \item{gw} { set of genes driving the significant aspects }
##' \item{df} { text table with the significance testing results }
##' }
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' }
##'
##' @export
pagoda.top.aspects <- function(pwpca, clpca = NULL, n.cells = NULL, z.score = qnorm(0.05/2, lower.tail = FALSE), return.table = FALSE, return.genes = FALSE, plot = FALSE, adjust.scores = TRUE, score.alpha = 0.05, use.oe.scale = FALSE, effective.cells.start = NULL) {
    basevar = 1

    if(is.null(n.cells)) {
        n.cells <- pagoda.effective.cells(pwpca, start = effective.cells.start)
    }


    vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
        vars <- as.numeric((pwpca[[i]]$sd)^2)
        shz <- NA
        if(!is.null(pwpca[[i]]$xp$randvar)) { shz <- (vars - mean(pwpca[[i]]$xp$randvar))/sd(pwpca[[i]]$xp$randvar) }
        cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$scores)), shz = shz)
    })))

    # fix p-to-q mistake in qWishartSpike
    qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
        params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
        qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
    }

    # add right tail approximation to ptw, which gives up quite early
    pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
        params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
        q.tw <- (q - params$centering)/(params$scaling)
        p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
        p[p == -Inf] <- pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
        p
    }


    #bi <- which.max(unlist(lapply(pwpca, function(x) x$n)))
    #vshift <- mean(pwpca[[bi]]$z[, 1]^2)/pwpca[[bi]]$n
    #ev <- ifelse(spike > 0, qWishartSpikeFixed(0.5, spike, n.cells, pwpca[[bi]]$n, var = basevar, lower.tail = FALSE), RMTstat::qWishartMax(0.5, n.cells, pwpca[[bi]]$n, var = basevar, lower.tail = FALSE))/pwpca[[bi]]$n
    #cat("vshift = ", vshift)
    vshift <- 0
    ev <- 0

    vdf$var <- vdf$var-(vshift-ev)*vdf$n

    #vdf$var[vdf$npc == 1] <- vdf$var[vdf$npc == 1]-(vshift-ev)*vdf$n[vdf$npc == 1]
    vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    vdf$z <- qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    vdf$cz <- qnorm(bh.adjust(pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
    vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)

    if(!is.null(clpca)) {
        clpca$xf <- extRemes::fevd(varst, data = clpca$varm, type = "Gumbel")
        #clpca$xf <- fevd(clpca$varm$varst[clpca$tci], type = "Gumbel")
        clpca$xf$results$par <- c(clpca$xf$results$par, c(shape = 0))
        #plot(xf)

        clvdf <- data.frame(do.call(rbind, lapply(seq_along(clpca$cl.goc), function(i)  {
            vars <- as.numeric((clpca$cl.goc[[i]]$sd)^2)
            shz <- NA
            if(!is.null(clpca$cl.goc[[i]]$xp$randvar)) {
                shz <- (vars - mean(clpca$cl.goc[[i]]$xp$randvar))/sd(clpca$cl.goc[[i]]$xp$randvar)
            }
            cbind(i = i, var = vars, n = clpca$cl.goc[[i]]$n, npc = seq(1:ncol(clpca$cl.goc[[i]]$xp$scores)), shz = shz)
        })))

        clvdf$var <- clvdf$var-(vshift-ev)*clvdf$n

        x <- RMTstat::WishartMaxPar(n.cells, clvdf$n)
        clvdf$pm <- x$centering-(1.2065335745820)*x$scaling # predicted mean of a random set
        clvdf$pv <- (1.607781034581)*x$scaling # predicted variance of a random set
        pvar <- predict(clpca$clvlm, newdata = clvdf)
        clvdf$varst <- (clvdf$var-pvar)/sqrt(clvdf$pv)
        clvdf$exp <- clpca$xf$results$par[1]*sqrt(clvdf$pv)+pvar
        #clvdf$varst <- (clvdf$var-clvdf$pm)/sqrt(clvdf$pv)
        #clvdf$exp <- clpca$xf$results$par[1]*sqrt(clvdf$pv)+clvdf$pm

        lp <- pgev.upper.log(clvdf$varst, clpca$xf$results$par[1], clpca$xf$results$par[2], rep(clpca$xf$results$par[3], nrow(clvdf)))
        clvdf$z <- qnorm(lp, lower.tail = FALSE, log.p = TRUE)
        clvdf$cz <- qnorm(bh.adjust(pnorm(as.numeric(clvdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)

        # CI relative to the background
        clvdf$ub <- extRemes::qevd(score.alpha/2, loc = clpca$xf$results$par[1], scale = clpca$xf$results$par[2], shape = clpca$xf$results$par[3], lower.tail = FALSE)*sqrt(clvdf$pv) + pvar
        clvdf$ub.stringent <- extRemes::qevd(score.alpha/2/nrow(clvdf), loc = clpca$xf$results$par[1], scale = clpca$xf$results$par[2], shape = clpca$xf$results$par[3], lower.tail = FALSE)*sqrt(clvdf$pv) + pvar

    }

    if(plot) {
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
        un <- sort(unique(vdf$n))
        on <- order(vdf$n, decreasing = FALSE)
        pccol <- colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
        plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = pccol[vdf$npc])
        lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
        lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)

        if(!is.null(clpca)) {
            pccol <- colorRampPalette(c("darkgreen", "lightgreen"), space = "Lab")(max(clvdf$npc))
            points(clvdf$n, clvdf$var/clvdf$n, col = pccol[clvdf$npc], pch = 1)

            #clvm <- clpca$xf$results$par[1]*sqrt(pmax(1e-3, predict(vm, data.frame(n = un)))) + predict(mm, data.frame(n = un))
            on <- order(clvdf$n, decreasing = FALSE)

            lines(clvdf$n[on], (clvdf$exp/clvdf$n)[on], col = "darkgreen")
            lines(clvdf$n[on], (clvdf$ub.stringent/clvdf$n)[on], col = "darkgreen", lty = 2)
        }
        #mi<-which.max(vdf$n) sv<- (vdf$var/vdf$n)[mi] - (vdf$exp/vdf$n)[mi]
        #lines(vdf$n[on], (vdf$exp/vdf$n)[on]+sv, col = 2, lty = 3)
        #lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on]+sv, col = 2, lty = 2)
    }


    if(!is.null(clpca)) { # merge in cluster stats based on their own model

        # merge pwpca, psd and pm
        # all processing from here is common
        clvdf$i <- clvdf$i+length(pwpca) # shift cluster ids
        pwpca <- c(pwpca, clpca$cl.goc)
        vdf <- rbind(vdf, clvdf[, c("i", "var", "n", "npc", "exp", "cz", "z", "ub", "ub.stringent", "shz")])
    }

    vdf$adj.shz <- qnorm(bh.adjust(pnorm(as.numeric(vdf$shz), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
    #vdf$oe <- vdf$var/vdf$exp
    rs <- (vshift-ev)*vdf$n
    #rs <- ifelse(vdf$npc == 1, (vshift-ev)*vdf$n, 0)
    vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
    #vdf$oe[vdf$oe<0] <- 0
    #vdf$oec <- (vdf$var-vdf$ub.stringent+vdf$exp)/vdf$exp
    #vdf$oec <- (vdf$var-vdf$ub+vdf$exp)/vdf$exp
    #vdf$oec <- (vdf$var-vdf$ub+vdf$exp+rs)/(vdf$exp+rs)
    vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)
    #vdf$oec[vdf$oec<0] <- 0
    #vdf$z[vdf$z<0] <- 0



    df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, sh.z = vdf$shz, adj.sh.z = vdf$adj.shz, stringsAsFactors = FALSE)
    if(adjust.scores) {
        vdf$valid <- vdf$cz  >=  z.score
    } else {
        vdf$valid <- vdf$z  >=  z.score
    }

    if(return.table) {
        df <- df[vdf$valid, ]
        df <- df[order(df$score, decreasing = TRUE), ]
        return(df)
    }

    # determine genes driving significant pathways
    # return genes within top 2/3rds of PC loading
    gl <- lapply(which(vdf$valid), function(i) { s <- abs(pwpca[[vdf[i, "i"]]]$xp$rotation[, vdf[i, "npc"]] )
    s[s >= max(s)/3] })
    gw <- tapply(abs(unlist(gl)), as.factor(unlist(lapply(gl, names))), max)
    if(return.genes) {
        return(gw)
    }
    # return combined data structure

    # weight
    xvw <- do.call(rbind, lapply(pwpca, function(x) {
        xm <- t(x$xp$scoreweights)
    }))
    vi <- vdf$valid
    xvw <- xvw[vi, ]/rowSums(xvw[vi, ])

    # return scaled patterns
    xmv <- do.call(rbind, lapply(pwpca, function(x) {
        xm <- t(x$xp$scores)
    }))

    if(use.oe.scale) {
        xmv <- (xmv[vi, ] -rowMeans(xmv[vi, ]))* (as.numeric(vdf$oe[vi])/sqrt(apply(xmv[vi, ], 1, var)))
    } else {
        # chi-squared
        xmv <- (xmv[vi, ]-rowMeans(xmv[vi, ])) * sqrt((qchisq(pnorm(vdf$z[vi], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv[vi, ], 1, var))
    }
    rownames(xmv) <- paste("#PC", vdf$npc[vi], "# ", names(pwpca)[vdf$i[vi]], sep = "")

    return(list(xv = xmv, xvw = xvw, gw = gw, df = df))

}


##' Collapse aspects driven by the same combinations of genes
##'
##' Examines PC loading vectors underlying the identified aspects and clusters aspects based
##' on a product of loading and score correlation (raised to corr.power). Clusters of aspects
##' driven by the same genes are determined based on the distance.threshold and collapsed.
##'
##' @param tam output of pagoda.top.aspects()
##' @param pwpca output of pagoda.pathway.wPCA()
##' @param clpca output of pagoda.gene.clusters() (optional)
##' @param plot whether to plot the resulting clustering
##' @param cluster.method one of the standard clustering methods to be used (fastcluster::hclust is used if available or stats::hclust)
##' @param distance.threshold similarity threshold for grouping interdependent aspects
##' @param corr.power power to which the product of loading and score correlation is raised
##' @param abs Boolean of whether to use absolute correlation
##' @param n.cores number of cores to use during processing
##' @param ... additional arguments are passed to the pagoda.view.aspects() method during plotting
##'
##' @return a list structure analogous to that returned by pagoda.top.aspects(), but with addition of a $cnam element containing a list of aspects summarized by each row of the new (reduced) $xv and $xvw
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)
##' }
##'
##' @export
pagoda.reduce.loading.redundancy <- function(tam, pwpca, clpca = NULL, plot = FALSE, cluster.method = "complete", distance.threshold = 0.01, corr.power = 4, n.cores = 1, abs = TRUE, ...) {
    pclc <- pathway.pc.correlation.distance(c(pwpca, clpca$cl.goc), tam$xv, target.ndf = 100, n.cores = n.cores)
    cda <- cor(t(tam$xv))
    if(abs) {
        cda <- abs(cda)
    } else {
        cda[cda<0] <- 0
    }
    cda <- as.dist(1-cda)
    cc <- (1-sqrt((1-pclc)*(1-cda)))^corr.power

    if(is.element("fastcluster", installed.packages()[, 1])) {
        y <- fastcluster::hclust(cc, method = cluster.method)
    } else {
        y <- stats::hclust(cc, method = cluster.method)
    }
    ct <- cutree(y, h = distance.threshold)
    ctf <- factor(ct, levels = sort(unique(ct)))
    xvl <- collapse.aspect.clusters(tam$xv, tam$xvw, ct, pick.top = FALSE, scale = TRUE)

    if(plot) {
        sc <- sample(colors(), length(levels(ctf)), replace = TRUE)
        view.aspects(tam$xv, row.clustering = y, row.cols = sc[as.integer(ctf)], ...)
    }

    # collapsed names
    if(!is.null(tam$cnam)) { # already has collapsed names
        cnam <- tapply(rownames(tam$xv), ctf, function(xn) unlist(tam$cnam[xn]))
    } else {
        cnam <- tapply(rownames(tam$xv), ctf, I)
    }
    names(cnam) <- rownames(xvl$d)
    tam$xv <- xvl$d
    tam$xvw <- xvl$w
    tam$cnam <- cnam
    return(tam)
}


##' Collapse aspects driven by similar patterns (i.e. separate the same sets of cells)
##'
##' Examines PC loading vectors underlying the identified aspects and clusters aspects based on score correlation. Clusters of aspects driven by the same patterns are determined based on the distance.threshold.
##'
##' @param tamr output of pagoda.reduce.loading.redundancy()
##' @param distance.threshold similarity threshold for grouping interdependent aspects
##' @param cluster.method one of the standard clustering methods to be used (fastcluster::hclust is used if available or stats::hclust)
##' @param distance distance matrix
##' @param weighted.correlation Boolean of whether to use a weighted correlation in determining the similarity of patterns
##' @param plot Boolean of whether to show plot
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param trim Winsorization trim to use prior to determining the top aspects
##' @param abs Boolean of whether to use absolute correlation
##' @param ... additional arguments are passed to the pagoda.view.aspects() method during plotting
##'
##' @return a list structure analogous to that returned by pagoda.top.aspects(), but with addition of a $cnam element containing a list of aspects summarized by each row of the new (reduced) $xv and $xvw
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)
##' tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
##' }
##'
##' @export
pagoda.reduce.redundancy <- function(tamr, distance.threshold = 0.2, cluster.method = "complete", distance = NULL, weighted.correlation = TRUE, plot = FALSE, top = Inf, trim = 0, abs = FALSE, ...) {
    if(is.null(distance)) {
        if(weighted.correlation) {
            distance <- .Call("matWCorr", t(tamr$xv), t(tamr$xvw), PACKAGE = "scde")
            rownames(distance) <- colnames(distance) <- rownames(tamr$xv)
            if(abs) {
                distance <- stats::as.dist(1-abs(distance), upper = TRUE)
            } else {
                distance <- stats::as.dist(1-distance, upper = TRUE)
            }
        } else {
            if(abs) {
                distance <- stats::as.dist(1-abs(cor(t(tamr$xv))))
            } else {
                distance <- stats::as.dist(1-cor(t(tamr$xv)))
            }
        }
    }
    if(is.element("fastcluster", installed.packages()[, 1])) {
        y <- fastcluster::hclust(distance, method = cluster.method)
    } else {
        y <- stats::hclust(distance, method = cluster.method)
    }

    ct <- cutree(y, h = distance.threshold)
    ctf <- factor(ct, levels = sort(unique(ct)))
    xvl <- collapse.aspect.clusters(tamr$xv, tamr$xvw, ct, pick.top = FALSE, scale = TRUE)

    if(plot) {
        sc <- sample(colors(), length(levels(ctf)), replace = TRUE)
        view.aspects(tamr$xv, row.clustering = y, row.cols = sc[as.integer(ctf)], ...)
    }

    # collapsed names
    if(!is.null(tamr$cnam)) { # already has collapsed names
        cnam <- tapply(rownames(tamr$xv), ctf, function(xn) unlist(tamr$cnam[xn]))
    } else {
        cnam <- tapply(rownames(tamr$xv), ctf, I)
    }
    names(cnam) <- rownames(xvl$d)

    if(trim > 0) { xvl$d <- winsorize.matrix(xvl$d, trim) } # trim prior to determining the top sets

    rcmvar <- apply(xvl$d, 1, var)
    vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]

    tamr2 <- tamr
    tamr2$xv <- xvl$d[vi, ]
    tamr2$xvw <- xvl$w[vi, ]
    tamr2$cnam <- cnam[vi]
    return(tamr2)
}


##' Determine optimal cell clustering based on the genes driving the significant aspects
##'
##' Determines cell clustering (hclust result) based on a weighted correlation of genes
##' underlying the top aspects of transcriptional heterogeneity. Branch orientation is optimized
##' if 'cba' package is installed.
##'
##' @param tam result of pagoda.top.aspects() call
##' @param varinfo result of pagoda.varnorm() call
##' @param method clustering method ('ward.D' by default)
##' @param verbose 0 or 1 depending on level of desired verbosity
##' @param include.aspects whether the aspect patterns themselves should be included alongside with the individual genes in calculating cell distance
##' @param return.details Boolean of whether to return just the hclust result or a list containing the hclust result plus the distance matrix and gene values
##'
##' @return hclust result
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' hc <- pagoda.cluster.cells(tam, varinfo)
##' plot(hc)
##' }
##'
##' @export
pagoda.cluster.cells <- function(tam, varinfo, method = "ward.D", include.aspects = FALSE, verbose = 0, return.details = FALSE) {
    # gene clustering
    gw <- tam$gw
    gw <- gw[(rowSums(varinfo$matw)*varinfo$arv)[names(gw)] > 1]

    gw <- gw/gw
    mi <- match(names(gw), rownames(varinfo$mat))
    wgm <- varinfo$mat[mi, ]
    wgm <- wgm*as.numeric(gw)
    wgwm <- varinfo$matw[mi, ]

    if(include.aspects) {
        if(verbose) { message("clustering cells based on ", nrow(wgm), " genes and ", nrow(tam$xv), " aspect patterns")}
        wgm <- rbind(wgm, tam$xv)
        wgwm <- rbind(wgwm, tam$xvw)
    } else {
        if(verbose) { message("clustering cells based on ", nrow(wgm), " genes")}
    }

    snam <- sample(colnames(wgm))

    dm <- .Call("matWCorr", wgm, wgwm, PACKAGE = "scde")
    dm <- 1-dm
    rownames(dm) <- colnames(dm) <- colnames(wgm)
    wcord <- stats::as.dist(dm, upper = TRUE)
    hc <- hclust(wcord, method = method)

    if(is.element("cba", installed.packages()[, 1])) {
        co <- cba::order.optimal(wcord, hc$merge)
        hc$merge <- co$merge
        hc$order <- co$order
    }
    if(return.details) {
        return(list(clustering = hc, distance = wcord, genes = gw))
    } else {
        return(hc)
    }
}


##' View PAGODA output
##'
##' Create static image of PAGODA output visualizing cell hierarchy and top aspects of transcriptional heterogeneity
##'
##' @param tamr Combined pathways that show similar expression patterns. Output of \code{\link{pagoda.reduce.redundancy}}
##' @param row.clustering Dendrogram of combined pathways clustering
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param ... additional arguments are passed to the \code{\link{view.aspects}} method during plotting
##'
##' @return PAGODA heatmap
##'
##' @examples
##' data(pollen)
##' cd <- clean.counts(pollen)
##' \donttest{
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' pagoda.view.aspects(tam)
##' }
##'
##' @export
pagoda.view.aspects <- function(tamr, row.clustering = hclust(dist(tamr$xv)), top = Inf, ...) {
    if(is.finite(top)) {
        rcmvar <- apply(tamr$xv, 1, var)
        vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]
        tamr$xv <- tamr$xv[vi, ]
        tamr$xvw <- tamr$xvw[vi, ]
        tamr$cnam <- tamr$cnam[vi]
    }

    view.aspects(tamr$xv, row.clustering = row.clustering, ... )
}


##' View heatmap
##'
##' Internal function to visualize aspects of transcriptional heterogeneity as a heatmap. Used by \code{\link{pagoda.view.aspects}}.
##'
##' @param mat Numeric matrix
##' @param row.clustering Row dendrogram
##' @param cell.clustering Column dendrogram
##' @param zlim Range of the normalized gene expression levels, inputted as a list: c(lower_bound, upper_bound). Values outside this range will be Winsorized. Useful for increasing the contrast of the heatmap visualizations. Default, set to the 5th and 95th percentiles.
##' @param row.cols  Matrix of row colors.
##' @param col.cols  Matrix of column colors. Useful for visualizing cell annotations such as batch labels.
##' @param cols Heatmap colors
##' @param show.row.var.colors Boolean of whether to show row variance as a color track
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param ... additional arguments for heatmap plotting
##'
##' @return A heatmap
##'
view.aspects <- function(mat, row.clustering = NA, cell.clustering = NA, zlim = c(-1, 1)*quantile(mat, p = 0.95), row.cols = NULL, col.cols = NULL, cols = colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(1024), show.row.var.colors = TRUE, top = Inf, ...) {
    #row.cols, col.cols are matrices for now
    rcmvar <- apply(mat, 1, var)
    mat[mat<zlim[1]] <- zlim[1]
    mat[mat > zlim[2]] <- zlim[2]
    if(class(row.clustering) == "hclust") { row.clustering <- as.dendrogram(row.clustering) }
    if(class(cell.clustering) == "hclust") { cell.clustering <- as.dendrogram(cell.clustering) }
    if(show.row.var.colors) {
        if(is.null(row.cols)) {
            icols <- colorRampPalette(c("white", "black"), space = "Lab")(1024)[1023*(rcmvar/max(rcmvar))+1]
            row.cols <- cbind(var = icols)
        }
    }
    my.heatmap2(mat, Rowv = row.clustering, Colv = cell.clustering, zlim = zlim, RowSideColors = row.cols, ColSideColors = col.cols, col = cols, ...)
}


##' Make the PAGODA app
##'
##' Create an interactive user interface to explore output of PAGODA.
##'
##' @param tamr Combined pathways that show similar expression patterns. Output of \code{\link{pagoda.reduce.redundancy}}
##' @param tam Combined pathways that are driven by the same gene sets. Output of \code{\link{pagoda.reduce.loading.redundancy}}
##' @param varinfo Variance information. Output of \code{\link{pagoda.varnorm}}
##' @param env Gene sets as an environment variable.
##' @param pwpca Weighted PC magnitudes for each gene set provided in the \code{env}. Output of \code{\link{pagoda.pathway.wPCA}}
##' @param clpca Weighted PC magnitudes for de novo gene sets identified by clustering on expression. Output of \code{\link{pagoda.gene.clusters}}
##' @param col.cols  Matrix of column colors. Useful for visualizing cell annotations such as batch labels. Default NULL.
##' @param cell.clustering Dendrogram of cell clustering. Output of \code{\link{pagoda.cluster.cells} } . Default   NULL.
##' @param row.clustering Dendrogram of combined pathways clustering. Default NULL.
##' @param title Title text to be used in the browser label for the app. Default, set as 'pathway clustering'
##' @param zlim Range of the normalized gene expression levels, inputted as a list: c(lower_bound, upper_bound). Values outside this range will be Winsorized. Useful for increasing the contrast of the heatmap visualizations. Default, set to the 5th and 95th percentiles.
##'
##' @return PAGODA app
##'
##' @export
make.pagoda.app <- function(tamr, tam, varinfo, env, pwpca, clpca = NULL, col.cols = NULL, cell.clustering = NULL, row.clustering = NULL, title = "pathway clustering", zlim = c(-1, 1)*quantile(tamr$xv, p = 0.95),inchlib=T) {
    # rcm - xv
    
    # matvar
    if(is.null(cell.clustering)) {
        cell.clustering <- pagoda.cluster.cells(tam, varinfo)
    }
    if(is.null(row.clustering) || is.null(row.clustering$order)) {
      row.clustering <- hclust(dist(tamr$xv))
    } else if(class(row.clustering)!="hclust") {
      # make a fake clustering to match the provided order
      or <- row.clustering$order;
      row.clustering <- hclust(dist(tamr3$xv),method='single')
      names(or) <- as.character(-1*row.clustering$order)
      nmm <- -1*or[as.character(row.clustering$merge)]
      nmm[is.na(nmm)] <- as.character(row.clustering$merge)[is.na(nmm)]
      row.clustering$merge <- matrix(as.integer(nmm),ncol=ncol(row.clustering$merge))
      row.clustering$order <- as.integer(or);
    }

    #fct - which tam row in which tamr$xv cluster.. remap tamr$cnams
    cn <- tamr$cnam
    fct <- rep(1:length(cn), lapply(cn, length))
    names(fct) <- unlist(cn)
    fct <- fct[rownames(tam$xv)]
    rcm <- tamr$xv
    #rownames(rcm) <- as.character(1:nrow(rcm))
    fres <- list(hvc = cell.clustering, tvc = row.clustering, rcm = rcm, zlim2 = zlim, matvar = apply(tam$xv, 1, sd), ct = fct, matrcmcor = rep(1, nrow(tam$xv)), cols = colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(1024), colcol = col.cols)

    # gene df
    gene.df <- data.frame(var = varinfo$arv*rowSums(varinfo$matw))
    gene.df$gene <- rownames(varinfo$mat)
    gene.df <- gene.df[order(gene.df$var, decreasing = TRUE), ]

    # prepare pathway df
    df <- tamr$df
    if(exists("myGOTERM", envir = globalenv())) {
        df$desc <- mget(df$name, get("myGOTERM", envir = globalenv()), ifnotfound = "")
    } else {
        df$desc <- ""
    }
    min.z <- -9
    df$z[df$z<min.z] <- min.z
    df$adj.z[df$adj.z<min.z] <- min.z
    df$sh.z[df$sh.z<min.z] <- min.z
    df$adj.sh.z[df$adj.sh.z<min.z] <- min.z
    df <- data.frame(id = paste("#PC", df$npc, "# ", df$name, sep = ""), npc = df$npc, n = df$n, score = df$score, Z = df$z, aZ = df$adj.z, sh.Z = df$sh.z, sh.aZ = df$adj.sh.z, name = paste(df$name, df$desc))

    df <- df[order(df$score, decreasing = TRUE), ]

    # merge go.env
    if(!is.null(clpca)) {
        set.env <- list2env(c(as.list(env), clpca$clusters))
    } else {
        set.env <- env
    }
  if(inchlib) {
    sa <- ViewPagodaApp$new(fres, df, gene.df, varinfo$mat, varinfo$matw, set.env, name = title, trim = 0, batch = varinfo$batch)
  } else {
    sa <- ViewPagodaAppOld$new(fres, df, gene.df, varinfo$mat, varinfo$matw, set.env, name = title, trim = 0, batch = varinfo$batch)
  }
}

