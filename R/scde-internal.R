
one.sided.test.id <- function(id, nam1, nam2, ifm, dm, prior, difference.prior = 0.5, bootstrap = TRUE, n.samples = 1e3, show.plots = TRUE, return.posterior = FALSE, return.both = FALSE) {
    gr <- 10^prior$x - 1
    gr[gr<0] <- 0
    lpp <- get.rep.set.general.model.logposteriors(ifm[[nam1]], dm[rep(id, length(gr)), names(ifm[[nam1]])], data.frame(fpm = gr), grid.weight = prior$grid.weight)
    ldp <- get.rep.set.general.model.logposteriors(ifm[[nam2]], dm[rep(id, length(gr)), names(ifm[[nam2]])], data.frame(fpm = gr), grid.weight = prior$grid.weight)

    if(bootstrap) {
        pjp <- do.call(cbind, lapply(seq_along(n.samples), function(i) {
            pjp <- rowSums(lpp[, sample(1:ncol(lpp), replace = TRUE)])
            pjp <- exp(pjp-max(pjp))
            pjp <- pjp/sum(pjp)
            return(pjp)
        }))
        pjp <- rowSums(pjp)
        pjp <- log(pjp/sum(pjp))

        djp <- do.call(cbind, lapply(seq_along(n.samples), function(i) {
            djp <- rowSums(ldp[, sample(1:ncol(ldp), replace = TRUE)])
            djp <- exp(djp-max(djp))
            djp <- djp/sum(djp)
            return(djp)
        }))
        djp <- rowSums(djp)
        djp <- log(djp/sum(djp))
    } else {
        pjp <- rowSums(lpp)
        djp <- rowSums(ldp)
    }

    dpy <- exp(prior$lp+djp)
    mpgr <- sum(exp(prior$lp+pjp+log(c(0, cumsum(dpy)[-length(dpy)])))) # m1
    mpls <- sum(exp(prior$lp+pjp+log(sum(dpy)-cumsum(dpy)))) # m0
    mpls/mpgr

    pjpc <- exp(prior$lp+pjp)
    pjpc <- pjpc/sum(pjpc)
    djpc <- exp(prior$lp+djp)
    djpc <- djpc/sum(djpc)

    if(show.plots || return.posterior || return.both) {
        # calculate log-fold-change posterior
        n <- length(pjpc)
        rp <- c(unlist(lapply(n:2, function(i) sum(pjpc[1:(n-i+1)]*djpc[i:n]))), unlist(lapply(seq_along(n), function(i) sum(pjpc[i:n]*djpc[1:(n-i+1)]))))
        rv <- seq(prior$x[1]-prior$x[length(prior$x)], prior$x[length(prior$x)]-prior$x[1], length = length(prior$x)*2-1)
        fcp <- data.frame(v = rv, p = rp)
    }

    if(show.plots) {
        # show each posterior
        layout(matrix(c(1:3), 3, 1, byrow = TRUE), heights = c(2, 1, 2), widths = c(1), FALSE)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        jpr <- range(c(0, pjpc), na.rm = TRUE)
        pp <- exp(lpp)
        cols <- rainbow(dim(pp)[2], s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, pp)), xlab = "expression level", ylab = "individual posterior", main = nam1)
        lapply(seq_len(ncol(pp)), function(i) lines(prior$x, pp[, i], col = cols[i]))
        legend(x = ifelse(which.max(na.omit(pjpc)) > length(pjpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = colnames(pp), lty = rep(1, dim(pp)[2]))
        par(new = TRUE)
        plot(prior$x, pjpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)

        # ratio plot
        par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        plot(fcp$v, fcp$p, xlab = "log10 expression ratio", ylab = "ratio posterior", type = 'l', lwd = 2, main = "")
        r.mle <- fcp$v[which.max(fcp$p)]
        r.lb <- max(which(cumsum(fcp$p)<0.025))
        r.ub <- min(which(cumsum(fcp$p) > (1-0.025)))
        polygon(c(fcp$v[r.lb], fcp$v[r.lb:r.ub], fcp$v[r.ub]), y = c(-10, fcp$p[r.lb:r.ub], -10), col = "grey90")
        abline(v = r.mle, col = 2, lty = 2)
        abline(v = c(fcp$v[r.ub], fcp$v[r.lb]), col = 2, lty = 3)
        box()
        legend(x = ifelse(r.mle > 0, "topleft", "topright"), legend = c(paste("MLE = ", round(10^(r.mle), 1), " (", round(r.mle, 2), " in log10)", sep = ""), paste("95% CI: ", round(10^(fcp$v[r.lb]), 1), " - ", round(10^(fcp$v[r.ub]), 1), sep = ""), paste(" log10: ", round(fcp$v[r.lb], 2), " - ", round(fcp$v[r.ub], 2), sep = "")), bty = "n")

        # distal plot
        dp <- exp(ldp)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        jpr <- range(c(0, djpc), na.rm = TRUE)
        cols <- rainbow(dim(dp)[2], s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, dp)), xlab = "expression level", ylab = "individual posterior", main = nam2)
        lapply(seq_len(ncol(dp)), function(i) lines(prior$x, dp[, i], col = cols[i]))
        legend(x = ifelse(which.max(na.omit(djpc)) > length(djpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = colnames(dp), lty = rep(1, dim(dp)[2]))

        par(new = TRUE)
        plot(prior$x, djpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)
    }

    lbf <- mpls/mpgr
    lbf <- (difference.prior*lbf)/(difference.prior*lbf+1-difference.prior)
    #return(c(equal = qnorm(ebf, lower.tail = TRUE), less = qnorm(lbf, lower.tail = TRUE)))
    if(return.both) {
        return(list(z = qnorm(lbf, lower.tail = TRUE), post = fcp))
    } else if(return.posterior) {
        return(fcp)
    } else {
        return(qnorm(lbf, lower.tail = TRUE))
    }
}

# counts - data frame with fragment counts (rows - fragments columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# min.count.threshold - the number of reads used to make an initial guess for the failed component
# threshold.segmentation - use min.count.threshold to perform very quick identification of the drop-outs
# threshold.prior - prior that should be associated with threshold segmentation
calculate.crossfit.models <- function(counts, groups, min.count.threshold = 4, nrep = 1, verbose = 0, min.prior = 1e-5, n.cores = 12, save.plots = TRUE, zero.lambda = 0.1, old.cfm = NULL, threshold.segmentation = FALSE, threshold.prior = 1-1e-6, max.pairs = 1000, min.pairs.per.cell = 10) {
    names(groups) <- colnames(counts)
    # enumerate cross-fit pairs within each group
    cl <- do.call(cbind, tapply(colnames(counts), groups, function(ids) {
        cl <- combn(ids, 2)
        min.pairs.per.cell <- min(length(ids)*(length(ids)-1)/2, min.pairs.per.cell)
        if(verbose) {
            cat("number of pairs: ", ncol(cl), "\n")
        }
        if(ncol(cl) > max.pairs) {
            if(verbose) {
                cat("reducing to a random sample of ", max.pairs, " pairs\n")
            }

            # make sure there's at least min.pairs.per.cell pairs for each cell
            cl <- cl[, unique(c(sample(1:ncol(cl), max.pairs),
                                unlist(lapply(ids, function(id) sample(which(colSums(cl == id) > 0), min.pairs.per.cell)))))]
        }
        cl
    }))

    orl <- c()
    if(!is.null(old.cfm)) {
        # check which pairs have already been fitted in compared in old.cfm
        pn1 <- unlist(apply(cl, 2, function(ii) paste(ii, collapse = ".vs.")))
        pn2 <- unlist(apply(cl, 2, function(ii) paste(rev(ii), collapse = ".vs."))) ### %%% use rev() to revert element order
        vi <- (pn1 %in% names(old.cfm)) | (pn2 %in% names(old.cfm))
        cl <- cl[, !vi, drop = FALSE]
        orl <- old.cfm[names(old.cfm) %in% c(pn1, pn2)]
    }
    if(verbose) {
        cat("total number of pairs: ", ncol(cl), "\n")
    }

    if(dim(cl)[2] > 0) {
        if(verbose)  message(paste("cross-fitting", ncol(cl), "pairs:"))
        rl <- papply(seq_len(ncol(cl)), function(cii) {
            ii <- cl[, cii]
            df <- data.frame(c1 = counts[, ii[1]], c2 = counts[, ii[2]])
            vi <- which(rowSums(df) > 0, )
            if(!threshold.segmentation) {
                if(verbose) {
                    message("fitting pair [", paste(ii, collapse = " "), "]")
                }
                mo1 <- FLXMRglmCf(c1~1, family = "poisson", components = c(1), mu = log(zero.lambda))
                mo2 <- FLXMRnb2glmC(c1~1+I(log(c2+1)), components = c(2))
                mo3 <- FLXMRnb2glmC(c2~1+I(log(c1+1)), components = c(2))
                mo4 <- FLXMRglmCf(c2~1, family = "poisson", components = c(3), mu = log(zero.lambda))
                m1 <- mc.stepFlexmix(c1~1, data = df[vi, ], k = 3, model = list(mo1, mo2, mo3, mo4), control = list(verbose = verbose, minprior = min.prior), concomitant = FLXPmultinom(~I((log(c1+1)+log(c2+1))/2)+1), cluster = cbind(df$c1[vi]<= min.count.threshold, df$c1[vi] > min.count.threshold & df$c2[vi] > min.count.threshold, df$c2[vi]<= min.count.threshold), nrep = nrep)

                # reduce return size
                m1@posterior <- lapply(m1@posterior, function(m) {
                    rownames(m) <- NULL
                    return(m)
                })
                #rownames(m1@concomitant@x) <- NULL
                m1@concomitant@x <- matrix()
                m1@model <- lapply(m1@model, function(mod) {
                    mod@x <- matrix()
                    mod@y <- matrix()
                    #rownames(mod@x) <- NULL
                    #rownames(mod@y) <- NULL
                    return(mod)
                })

                #parent.env(environment(m1@components[[1]][[1]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[1]][[2]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[2]][[1]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[2]][[2]]@logLik)) <- globalenv()

                names(vi) <- NULL
                pm <- posterior(m1)[, c(1, 3)]
                rownames(pm) <- NULL
                cl <- clusters(m1)
                names(cl) <- NULL
                gc()
            } else {
                # use min.count.threshold to quickly segment the points
                cl <- rep(2, length(vi))
                cl[df[vi, 1]<min.count.threshold] <- 1
                cl[df[vi, 2]<min.count.threshold] <- 3
                cl[df[vi, 1]<min.count.threshold & df[vi, 2]<min.count.threshold] <- 0
                names(cl) <- NULL
                pm <- cbind(ifelse(cl == 1, threshold.prior, 1-threshold.prior), ifelse(cl == 3, threshold.prior, 1-threshold.prior))
                rownames(pm) <- NULL
            }
            rli <- list(ii = ii, clusters = cl, posterior = pm, vi = vi)
            #message("return object size for pair [", paste(ii, collapse = " "), "] is ", round(object.size(rli)/(1024^2), 3), " MB")
            return(rli)
        }, n.cores = round(n.cores/nrep))
        #, mc.preschedule = threshold.segmentation) # mclapply function has preschedule
        names(rl) <- apply(cl, 2, paste, collapse = ".vs.")
        # clean up invalid entries
        rl <- rl[!unlist(lapply(rl, is.null))]
        rl <- rl[unlist(lapply(rl, is.list))]
        #names(rl) <- unlist(lapply(rl, function(d) paste(d$ii, collapse = ".vs.")))
    } else {
        rl <- c()
    }

    if(!is.null(old.cfm)) rl <- c(rl, orl)

    if(save.plots) {
        #require(Cairo) require(RColorBrewer)
        tapply(colnames(counts), groups, function(ids) {
            cl <- combn(ids, 2)
            group <- as.character(groups[ids[1]])
            # log-scale hist
            t.pairs.panel.hist <- function(x, i = NULL, ...) {
                usr <- par("usr")
                on.exit(par(usr))
                par(usr = c(usr[1:2], 0, 1.5) )
                vi <- x > 0
                h <- hist(x, plot = FALSE)
                breaks <- h$breaks
                nB <- length(breaks)
                y <- log10(h$counts)
                y <- y/max(y)
                rect(breaks[-nB], 0, breaks[-1], y, col = "gray60", ...)
            }
            t.pairs.smoothScatter.spearman <- function(x, y, i = NULL, j = NULL, cex = 0.8, ...) {
                vi <- x > 0 | y > 0
                smoothScatter(x[vi], y[vi], add = TRUE, useRaster = TRUE, ...)
                legend(x = "bottomright", legend = paste("sr = ", round(cor(x[vi], y[vi], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
            }
            # component assignment scatter
            t.panel.component.scatter <- function(x, y, i, j, cex = 0.8, ...) {
                if(!is.null(rl[[paste(ids[i], "vs", ids[j], sep = ".")]])) {
                    m1 <- rl[[paste(ids[i], "vs", ids[j], sep = ".")]]
                    # forward plot
                    vi <- which(x > 0 | y > 0)
                    ci <- vi[m1$clusters == 1]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
                    }

                    ci <- vi[m1$clusters == 3]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Greens")[-(1:3)])), cex = 2)
                    }
                    ci <- vi[m1$clusters == 2]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
                    }
                    legend(x = "topleft", pch = c(19), col = "blue", legend = paste("sr = ", round(cor(x[ci], y[ci], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
                    legend(x = "bottomright", pch = c(rep(19, 3)), col = c("red", "blue", "green"), legend = paste(round(unlist(tapply(m1$clusters, factor(m1$clusters, levels = c(1, 2, 3)), length))*100/length(vi), 1), "%", sep = ""), bty = "n", cex = cex)

                } else if(!is.null(rl[[paste(ids[i], "vs", ids[j], sep = ".")]])) {
                    m1 <- rl[[paste(ids[j], "vs", ids[i], sep = ".")]]
                    # reverse plot
                    vi <- which(x > 0 | y > 0)
                    ci <- vi[m1$clusters == 3]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
                    }

                    ci <- vi[m1$clusters == 1]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Greens")[-(1:3)])), cex = 2)
                    }
                    ci <- vi[m1$clusters == 2]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
                    }
                    legend(x = "topleft", pch = c(19), col = "blue", legend = paste("sr = ", round(cor(x[ci], y[ci], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
                    legend(x = "bottomright", pch = c(rep(19, 3)), col = c("red", "blue", "green"), legend = paste(round(unlist(tapply(m1$clusters, factor(m1$clusters, levels = c(3, 2, 1)), length))*100/length(vi), 1), "%", sep = ""), bty = "n", cex = cex)
                } else {
                    #message(paste("ERROR: unable to find model for i = ", i, "j = ", j))
                    message(paste("INFO: cross-fit plots: skipping model for i = ", i, "j = ", j, " (increase max.pairs parameter if needed"))
                }
            }
            #pdf(file = paste(group, "crossfits.pdf", sep = "."), width = 3*length(ids), height = 3*length(ids))
            CairoPNG(filename = paste(group, "crossfits.png", sep = "."), width = 250*length(ids), height = 250*length(ids))
            pairs.extended(log10(counts[, ids]+1), lower.panel = t.pairs.smoothScatter.spearman, upper.panel = t.panel.component.scatter, diag.panel = t.pairs.panel.hist, cex = 1.5)
            dev.off()
        })
    }

    return(rl)
}

# estimates library sizes based on the correlated components
# min.size.entries - minimal number of entries (genes) used to determine scaling factors for individual experiments
# counts - data frame with fragment counts (rows - fragments columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# cfm - cross-fit models (return of calculate.crossfit.models())
# vil - optional binary matrix (corresponding to counts) with 0s marking likely drop-out observations
# return value - library size vector in millions of reads
estimate.library.sizes <- function(counts, cfm, groups, min.size.entries = min(nrow(counts), 2e3), verbose = 0, return.details = FALSE, vil = NULL, ...) {
    #require(edgeR)
    names(groups) <- colnames(counts)
    # determine the set fragments that were not attributed to failure in any cross-comparison
    if(is.null(vil)) {
        #x <- lapply(cfm, function(d) { ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters != 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters != 3)]) names(ll) <- d$ii return(ll) })
        x <- lapply(cfm, function(d) { ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters > 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters %% 3  != 0)])
        names(ll) <- d$ii
        return(ll)
        })
        vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = colnames(counts)[!is.na(groups)]), function(l) {
            x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0
            x[is.na(x)] <- FALSE
            return(x)
        }))
    }

    # order entries by the number of non-failed experiments,
    # select entries for library size estimation
    ni <- cbind(1:nrow(counts), rowSums(vil))
    ni <- ni[order(ni[, 2], decreasing = TRUE), ]
    if(nrow(ni)<min.size.entries) {
        stop("The number of valid genes (", nrow(ni), ") is lower then the specified min.size.entries (", min.size.entries, "). Please either increase min.size.entries or lower min.nonfailed parameter to increase the number of valid genes")
    }
    if(ni[min.size.entries, 2]<ncol(vil)) {
        # if the min.size.entries -th gene has failures, take only min.size.entries genes
        gis <- ni[1:min.size.entries, 1]
    } else {
        # otherwise take all genes that have not failed in any experiment
        gis <- ni[ni[, 2] == ncol(vil), 1]
    }

    if(verbose)  message(paste("adjusting library size based on", length(gis), "entries"))
    f <- calcNormFactors(as.matrix(counts[gis, !is.na(groups)]), ...)
    f <- f/exp(mean(log(f)))
    ls <- colSums(counts[gis, !is.na(groups)])*f/1e6
    if(return.details) { return(list(ls = ls, vil = vil)) } else { return(ls) }
}

# an alternative prior estimation procedure that weights down contributions by failure probability
# and uses pre-scaled fpkm guesses for magnitude estimates
estimate.signal.prior <- function(fpkm, fail, length.out = 400, show.plot = FALSE, pseudo.count = 1, bw = 0.1, max.quantile = 0.999, max.value = NULL) {
    fpkm <- log10(exp(as.matrix(fpkm))+1)
    wts <- as.numeric(as.matrix(1-fail[, colnames(fpkm)]))
    wts <- wts/sum(wts)
    # fit density on a mirror image
    if(is.null(max.value)) {
        x <- as.numeric(fpkm)
        max.value <- as.numeric(quantile(x[x<Inf], p = max.quantile))
    }
    md <- density(c(-1*as.numeric(fpkm), as.numeric(fpkm)), bw = bw, weights = c(wts/2, wts/2), n = 2*length.out+1, from = -1*max.value, to = max.value)

    gep <- data.frame(x = md$x[-c(1:length.out)], y = md$y[-c(1:length.out)])
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

# counts - data frame with gene counts (rows - genes columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# cfm - cross-fit models (return of calculate.crossfit.models())
# min.nonfailed - minimal number of non-failed observations required for a gene to be used in the final model fitting
#  A minimum of either the specified value or number of experiments -1 will be used.
calculate.individual.models <- function(counts, groups, cfm, nrep = 1, verbose = 0, n.cores = 12, min.nonfailed = 2, min.size.entries = 2e3, zero.count.threshold = 10, save.plots = TRUE, linear.fit = TRUE, return.compressed.models = FALSE,  local.theta.fit = FALSE, theta.fit.range = c(1e-2, 1e2), ...) {
    names(groups) <- colnames(counts)
    # determine library size discarding non-zero entries
    ls <- estimate.library.sizes(counts, cfm, groups, min.size.entries, verbose = verbose, return.details = TRUE)

    # fit three-component models to unique pairs within each group
    mll <- tapply(colnames(counts), groups, function(ids) {
        cl <- combn(ids, 2)
        group <- as.character(groups[ids[1]])

        # incorporate cross-fit pairs from cfm
        pn1 <- unlist(apply(cl, 2, function(ii) paste(ii, collapse = ".vs.")))
        pn2 <- unlist(apply(cl, 2, function(ii) paste(rev(ii), collapse = ".vs."))) ### %%% use rev() to revert element order
        vi <- (pn1 %in% names(cfm)) | (pn2 %in% names(cfm)) # check both reverse and forward pairing
        #if(!all(vi)) stop("unable to find cross-fit models for the following pairs : ", paste(pn1[!vi]))
        if(!all(vi)) {
            if(verbose > 0) {
                if(verbose > 1) {
                    cat(paste("WARNING: unable to find cross-fit models for the following pairs : ", paste(pn1[!vi], collapse = " ")), "\n")
                } else {
                    cat("WARNING: unable to find cross-fit models for ", sum(!vi), " out of ", length(vi), " pairs. Using a subset.\n")
                }
            }
            # use a subset
            if(sum(vi) > 3) {
                pn1 <- pn1[vi]
                pn2 <- pn2[vi]
                vi <- vi[vi]
            } else {
                stop("less than 3 valid cross-fit pairs are available! giving up.")
            }
        }

        #rl <- cfm[vi]
        vi.names<-names(cfm)[names(cfm) %in% c(pn1, pn2)] ### a similar selection was done like this in calculate.crossfit.models() function
        rl <- cfm[vi.names]  ### with this sub-selection we select only sample pairs within the current group (e.g. pairs of ES)

        # determine the set genes that were not attributed to failure in any cross-comparison
        x <- lapply(rl, function(d) {
            ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters > 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters %% 3  != 0)])
            names(ll) <- d$ii
            return(ll)
        })
        vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = ids), function(l) {
            x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0
            x[is.na(x)] <- FALSE
            return(x)
        }))

        #x <- lapply(rl, function(d) { ll <- list((d$failures == 1), (d$failures == 2)) names(ll) <- d$ii return(ll) })
        #vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = ids), function(l) { x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0 x[is.na(x)] <- FALSE return(x) }))

        t.ls <- ls$ls[ids]
        adjust <- NULL
        if(!is.null(ls$adjustments)) { ls$adjustments[[groups[ids[1]]]] }
        # fit two-NB2 mixture for each experiment
        if(verbose) { message(paste("fitting", group, "models:")) }
        gc()

        # pair cell name matrix
        nm <- do.call(rbind, lapply(rl, function(x) x$ii))

        ml <- papply(seq_along(ids), function(i) { try({
            if(verbose)  message(paste(i, ":", ids[i]))
            # determine genes with sufficient number of non-failed observations in other experiments
            vi <- which(rowSums(vil[, -i, drop = FALSE]) >= min(length(ids)-1, min.nonfailed))
            fpm <- rowMeans(t(t(counts[vi, ids[-i], drop = FALSE])/(t.ls[-i])))
            if(!is.null(adjust)) { fpm <- adjust(fpm)  } # adjust for between-group systematic differences
            df <- data.frame(count = counts[vi, ids[i]], fpm = fpm)

            # reconstruct failure prior for the cell by averaging across
            # cross-cell comparisons where the cell did participate
            cp <- exp(rowMeans(log(cbind(
                do.call(cbind, lapply(rl[which(nm[, 1] == ids[i])], function(d) {
                    ivi <- rep(NA, nrow(counts))
                    ivi[d$vi] <- 1:length(d$vi)
                    d$posterior[ivi[vi], 1]
                })),
                do.call(cbind, lapply(rl[which(nm[, 2] == ids[i])], function(d) {
                    ivi <- rep(NA, nrow(counts))
                    ivi[d$vi] <- 1:length(d$vi)
                    d$posterior[ivi[vi], 2]
                }))
            )), na.rm = TRUE))
            cp <- cbind(cp, 1-cp)

            nai <- which(is.na(cp[, 1]))
            cp[nai, 1] <- 1-(1e-10)
            cp[nai, 2] <- (1e-10)
            if(linear.fit) {
                m1 <- fit.nb2gth.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = zero.count.threshold, full.theta.range = theta.fit.range, theta.fit.range = theta.fit.range, use.constant.theta.fit = !local.theta.fit, ...)
            } else {
                m1 <- fit.nb2.mixture.model(df, prior = cp, nrep = nrep, verbose = verbose, zero.count.threshold = zero.count.threshold, ...)
            }

            if(return.compressed.models) {
                v <- get.compressed.v1.model(m1)
                cl <- clusters(m1)
                rm(m1)
                gc()
                return(list(model = v, clusters = cl))
            }

            # otherwise try to reduce the size of a full model
            # reduce return size
            #m1@posterior <- lapply(m1@posterior, function(m) { rownames(m) <- NULL return(m)})
            m1@posterior <- NULL
            #rownames(m1@concomitant@x) <- NULL
            m1@concomitant@x <- matrix()
            m1@model <- lapply(m1@model, function(mod) {
                mod@x <- matrix()
                mod@y <- matrix()
                #rownames(mod@x) <- NULL
                #rownames(mod@y) <- NULL
                return(mod)
            })

            # make a clean copy of the internal environment
            t.cleanenv <- function(comp) {
                el <- list2env(as.list(environment(comp@logLik), all.names = TRUE), parent = globalenv())
                ep <- list2env(as.list(environment(comp@predict), all.names = TRUE), parent = globalenv())
                pf <- get("predict", envir = el)
                environment(pf) <- ep
                assign("predict", pf, envir = el)
                pf <- get("predict", envir = ep)
                environment(pf) <- ep
                assign("predict", pf, envir = ep)

                pf <- get("logLik", envir = el)
                environment(pf) <- el
                assign("logLik", pf, envir = el)
                pf <- get("logLik", envir = ep)
                environment(pf) <- el
                assign("logLik", pf, envir = ep)

                environment(comp@logLik) <- el
                environment(comp@predict) <- ep
                comp
            }
            m1@components <- lapply(m1@components, function(cl) lapply(cl, t.cleanenv))

            # clean up the formula environment (was causing multithreading problems)
            rm(list = ls(env = attr(m1@concomitant@formula, ".Environment")), envir = attr(m1@concomitant@formula, ".Environment"))
            gc()
            #rm(list = ls(env = attr(m1@formula, ".Environment")), envir = attr(m1@formula, ".Environment"))
            return(m1)
        })}, n.cores = n.cores) # end cell iteration

        # check if there were errors in the multithreaded portion
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

        if(save.plots && length(ml)>0) {
            # model fits
            #CairoPNG(filename = paste(group, "model.fits.png", sep = "."), width = 1024, height = 300*length(ids))
            pdf(file = paste(group, "model.fits.pdf", sep = "."), width = ifelse(linear.fit, 15, 13), height = 4)
            #l <- layout(matrix(seq(1, 4*length(ids)), nrow = length(ids), byrow = TRUE), rep(c(1, 1, 1, 0.5), length(ids)), rep(1, 4*length(ids)), FALSE)
            l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, ifelse(linear.fit, 1, 0.5)), 1), rep(1, 4), FALSE)
            par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
            invisible(lapply(seq_along(vic), function(j) {
                i <- vic[j];
                vi <- which(rowSums(vil[, -i, drop = FALSE]) >= min(length(ids)-1, min.nonfailed))
                df <- data.frame(count = counts[vi, ids[i]], fpm = rowMeans(t(t(counts[vi, ids[-i], drop = FALSE])/(t.ls[-i]))))
                plot.nb2.mixture.fit(ml[[j]], df, en = ids[i], do.par = FALSE, compressed.models = return.compressed.models)
            }))
            dev.off()
        }

        return(ml)

    }) # end group iteration

    if(return.compressed.models) {
        # make a joint model matrix
        jmm <- data.frame(do.call(rbind, lapply(mll, function(tl) do.call(rbind, lapply(tl, function(m) m$model)))))
        rownames(jmm) <- unlist(lapply(mll, names))
        # reorder in the original cell order
        attr(jmm, "groups") <- rep(names(mll), unlist(lapply(mll, length)))
        return(jmm)
    } else {
        return(mll)
    }
}


#######
## V1 optimized methods
#######

# gets an array summary of gam model structure (assumes a flat ifm list)
get.compressed.v1.models <- function(ifml) {
    data.frame(do.call(rbind, lapply(ifml, get.compressed.v1.model)))
}

# get a vector representation of a given model
get.compressed.v1.model <- function(m1) {
    if(class(m1@model[[2]]) == "FLXMRnb2gthC") { # linear fit model
        v <- c(m1@concomitant@coef[c(1:2), 2], get("coef", environment(m1@components[[1]][[1]]@predict)))
        names(v) <- c("conc.b", "conc.a", "fail.r")
        vth <- m1@components[[2]][[2]]@parameters$coef
        # translate mu regression from linear to log model
        v <- c(v, c("corr.b" = log(as.numeric(vth["corr.a"])), "corr.a" = 1), vth[-match("corr.a", names(vth))], "conc.a2" = m1@concomitant@coef[3, 2])
    } else { # original publication model
        v <- c(m1@concomitant@coef[, 2], get("coef", environment(m1@components[[1]][[1]]@predict)), m1@components[[2]][[2]]@parameters$coef, get("theta", environment(m1@components[[2]][[2]]@predict)))
        names(v) <- c("conc.b", "conc.a", "fail.r", "corr.b", "corr.a", "corr.theta")
    }
    v
}

# calculates posterior matrices (log scale) for a set of ifm models
calculate.posterior.matrices <- function(dat, ifm, prior, n.cores = 32, inner.cores = 4, outer.cores = round(n.cores/inner.cores)) {
    marginals <- data.frame(fpm = 10^prior$x - 1)
    marginals$fpm[marginals$fpm<0] <- 0
    lapply(ifm, function(group.ifm) {
        papply(sn(names(group.ifm)), function(nam) {
            df <- get.exp.logposterior.matrix(group.ifm[[nam]], dat[, nam], marginals, n.cores = inner.cores, grid.weight = prior$grid.weight)
            rownames(df) <- rownames(dat)
            colnames(df) <- as.character(prior$x)
            return(df)
        }, n.cores = n.cores)
    })
}

sample.posterior <- function(dat, ifm, prior, n.samples = 1, n.cores = 32) {
    marginals <- data.frame(fpm = 10^prior$x - 1)
    lapply(ifm, function(group.ifm) {
        papply(sn(names(group.ifm)), function(nam) {
            get.exp.sample(group.ifm[[nam]], dat[, nam], marginals, prior.x = prior$x, n = n.samples)
        }, n.cores = n.cores)
    })
}

# calculate joint posterior matrix for a given group of experiments
# lmatl - list of posterior matrices (log scale) for individual experiments
calculate.joint.posterior.matrix <- function(lmatl, n.samples = 100, bootstrap = TRUE, n.cores = 15) {
    if(bootstrap) {
        jpl <- papply(seq_len(n.cores), function(i) jpmatLogBoot(Matl = lmatl, Nboot = ceiling(n.samples/n.cores), Seed = i), n.cores = n.cores)
        jpl <- Reduce("+", jpl)
        jpl <- jpl/rowSums(jpl)
    } else {
        jpl <- Reduce("+", lmatl)
        jpl <- exp(jpl-log.row.sums(jpl))
    }
    rownames(jpl) <- rownames(lmatl[[1]])
    colnames(jpl) <- colnames(lmatl[[1]])
    jpl
}

# calculate joint posterior of a group defined by a composition vector
# lmatll - list of posterior matrix lists (as obtained from calculate.posterior.matrices)
# composition - a named vector, indicating the number of samples that should be drawn from each element of lmatll to compose a group
calculate.batch.joint.posterior.matrix <- function(lmatll, composition, n.samples = 100, n.cores = 15) {
    # reorder composition vector to match lmatll names
    jpl <- papply(seq_len(n.cores), function(i) jpmatLogBatchBoot(lmatll, composition[names(lmatll)], ceiling(n.samples/n.cores), i), n.cores = n.cores)
    jpl <- Reduce("+", jpl)
    jpl <- jpl/rowSums(jpl)
    #jpl <- jpmatLogBatchBoot(lmatll, composition[names(lmatll)], n.samples, n.cores)
    rownames(jpl) <- rownames(lmatll[[1]][[1]])
    colnames(jpl) <- colnames(lmatll[[1]][[1]])
    jpl
}

# calculates the likelihood of expression difference based on
# two posterior matrices (not adjusted for prior)
calculate.ratio.posterior <- function(pmat1, pmat2, prior, n.cores = 15, skip.prior.adjustment = FALSE) {
    n <- length(prior$x)
    if(!skip.prior.adjustment) {
        pmat1 <- t(t(pmat1)*prior$y)
        pmat2 <- t(t(pmat2)*prior$y)
    }

    chunk <- function(x, n) split(x, sort(rank(x) %% n))
    if(n.cores > 1) {
        x <- do.call(rbind, papply(chunk(1:nrow(pmat1), n.cores*5), function(ii) matSlideMult(pmat1[ii, , drop = FALSE], pmat2[ii, , drop = FALSE]), n.cores = n.cores))
    } else {
        x <- matSlideMult(pmat1, pmat2)
    }
    x <- x/rowSums(x)

    rv <- seq(prior$x[1]-prior$x[length(prior$x)], prior$x[length(prior$x)]-prior$x[1], length = length(prior$x)*2-1)
    colnames(x) <- as.character(rv)
    rownames(x) <- rownames(pmat1)
    return(x)
}

# quick utility function to get the difference Z score from the ratio posterior
get.ratio.posterior.Z.score <- function(rpost, min.p = 1e-15) {
    rpost <- rpost+min.p
    rpost <- rpost/rowSums(rpost)
    zi <- which.min(abs(as.numeric(colnames(rpost))))
    gs <- rowSums(rpost[, 1:(zi-1), drop = FALSE])
    zl <- pmin(0, qnorm(gs, lower.tail = FALSE))
    zg <- pmax(0, qnorm(gs+rpost[, zi, drop = FALSE], lower.tail = FALSE))
    z <- ifelse(abs(zl) > abs(zg), zl, zg)
}

# calculate a joint posterior matrix with bootstrap
jpmatLogBoot <- function(Matl, Nboot, Seed) {
    .Call("jpmatLogBoot", Matl, Nboot, Seed, PACKAGE = "scde")
}

# similar to the above, but compiles joint by sampling a pre-set
# number of different types (defined by Comp factor)
jpmatLogBatchBoot <- function(Matll, Comp, Nboot, Seed) {
    .Call("jpmatLogBatchBoot", Matll, Comp, Nboot, Seed, PACKAGE = "scde")
}

matSlideMult <- function(Mat1, Mat2) {
    .Call("matSlideMult", Mat1, Mat2, PACKAGE = "scde")
}

calculate.failure.p <- function(dat, ifm, n.cores = 32) {
    lapply(ifm, function(group.ifm) {
        lapply(sn(names(group.ifm)), function(nam) {
            get.concomitant.prob(group.ifm[[nam]], counts = dat[, nam])
        })
    })
}

# calculate failure probabilities across all cells for a given set
# of levels (lfpm - log(fpm) vector for all genes
calculate.failure.lfpm.p <- function(lfpm, ifm, n.cores = 32) {
    lapply(ifm, function(group.ifm) {
        lapply(sn(names(group.ifm)), function(nam) {
            get.concomitant.prob(group.ifm[[nam]], lfpm = lfpm)
        })
    })
}

# get expected fpm from counts
get.fpm.estimates <- function(m1, counts) {
    if(class(m1@components[[2]][[2]]) == "FLXcomponentE") {
        # gam do inverse interpolation
        b1 <- get("b1", envir = environment(m1@components[[2]][[2]]@predict))
        z <- approx(x = b1$fitted.values, y = b1$model$x, xout = counts, rule = 1:2)$y
        z[is.na(z)] <- -Inf
        z
    } else {
        # linear model
        par <- m1@components[[2]][[2]]@parameters
        if(!is.null(par[["linear"]])) {
            log((counts-par$coef[[1]])/par$coef[[2]])
        } else {
            (log(counts)-par$coef[[1]])/par$coef[[2]]
        }
    }
}

# rdf : count/fpm data frame
fit.nb2.mixture.model <- function(rdf, zero.count.threshold = 10, prior = cbind(rdf$count<= zero.count.threshold, rdf$count > zero.count.threshold), nrep = 3, iter = 50, verbose = 0, background.rate = 0.1, ...) {
    #mo1 <- FLXMRnb2glmC(count~1, components = c(1), theta.range = c(0.5, Inf))
    #mo1 <- FLXMRglmCf(count~1, components = c(1), family = "poisson", mu = 0.01)
    #mo1 <- FLXMRglmC(count~1, components = c(1), family = "poisson")
    mo1 <- FLXMRglmCf(count~1, family = "poisson", components = c(1), mu = log(background.rate))
    mo2 <- FLXMRnb2glmC(count~1+I(log(fpm)), components = c(2), theta.range = c(0.5, Inf))

    m1 <- mc.stepFlexmix(count~1, data = rdf, k = 2, model = list(mo1, mo2), control = list(verbose = verbose, minprior = 0, iter = iter), concomitant = FLXPmultinom(~I(log(fpm))+1), cluster = prior, nrep = nrep, ...)

    # check if the theta was underfit
    if(get("theta", envir = environment(m1@components[[2]][[2]]@logLik)) == 0.5) {
        # refit theta
        sci <- clusters(m1) == 2
        fit <- glm.nb.fit(m1@model[[2]]@x[sci, , drop = FALSE], m1@model[[2]]@y[sci], weights = rep(1, sum(sci)), offset = c(), init.theta = 0.5)
        assign("theta", value = fit$theta, envir = environment(m1@components[[2]][[2]]@logLik))
        m1@components[[2]][[2]]@parameters$coef <- fit$coefficients
        assign("coef", value = fit$coefficients, envir = environment(m1@components[[2]][[2]]@logLik))
        message("WARNING: theta was underfit, new theta = ", fit$theta)
    }

    return(m1)
}

fit.nb2gth.mixture.model <- function(rdf, zero.count.threshold = 10, prior = as.integer(rdf$count >= zero.count.threshold | rdf$fpm<median(rdf$fpm[rdf$count<zero.count.threshold]))+1, nrep = 0, verbose = 0 , full.theta.range = c(1e-2, 1e2), theta.fit.range = full.theta.range, theta.sp = 1e-2, use.constant.theta.fit = FALSE, alpha.weight.power = 1/2, iter = 50) {
    #mo1 <- FLXMRglmC(count~1, components = c(1), family = "poisson")
    #matrix(cbind(ifelse(rdf$count<= zero.count.threshold, 0.95, 0.05), ifelse(rdf$count > zero.count.threshold, 0.95, 0.05)))
    mo1 <- FLXMRglmCf(count~1, family = "poisson", components = c(1), mu = log(0.1))
    mo2 <- FLXMRnb2gthC(count~0+fpm, components = c(2), full.theta.range = full.theta.range, theta.fit.range = theta.fit.range, theta.fit.sp = theta.sp, constant.theta = use.constant.theta.fit, alpha.weight.power = alpha.weight.power)
    m1 <- mc.stepFlexmix(count~1, data = rdf, k = 2, model = list(mo1, mo2), control = list(verbose = verbose, minprior = 0, iter = iter), concomitant = FLXPmultinom(~I(log(fpm))+I(log(fpm)^2)+1), cluster = prior, nrep = nrep)
    return(m1)
}

# rdf : count/fpm data frame
# en : experiment name for plotting
# n.zero.windows - number of windows to visualize failure model fit
# m1 - fitted model
plot.nb2.mixture.fit <- function(m1, rdf, en, do.par = TRUE, n.zero.windows = 50, compressed.models = FALSE, bandwidth = 0.05) {
    #require(Cairo) require(RColorBrewer)
    if(do.par) {
        CairoPNG(filename = paste(en, "model.fit.png", sep = "."), width = 800, height = 300)
        l <- layout(matrix(c(1:4), 1, 4, byrow = TRUE), c(1, 1, 1, 0.5), rep(1, 4), FALSE)
        par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
    }
    smoothScatter(log10(rdf$fpm+1), log10(rdf$count+1), xlab = "expected FPM", ylab = "observed counts", main = paste(en, "scatter", sep = " : "), bandwidth = bandwidth)

    plot(c(), c(), xlim = range(log10(rdf$fpm+1)), ylim = range(log10(rdf$count+1)), xlab = "expected FPM", ylab = "observed counts", main = paste(en, "components", sep = " : "))
    if(compressed.models) {
        vpi <- m1$clusters == 1
    } else {
        vpi <- clusters(m1) == 1
    }
    if(sum(vpi) > 2){
        points(log10(rdf$fpm[vpi]+1), log10(rdf$count[vpi]+1), pch = ".", col = densCols(log10(rdf$fpm[vpi]+1), log10(rdf$count[vpi]+1), colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
    }
    if(sum(!vpi) > 2){
        points(log10(rdf$fpm[!vpi]+1), log10(rdf$count[!vpi]+1), pch = ".", col = densCols(log10(rdf$fpm[!vpi]+1), log10(rdf$count[!vpi]+1), colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
        # show fit
        fpmo <- order(rdf$fpm[!vpi], decreasing = FALSE)
        if(compressed.models) {
            #rf <- scde.failure.probability(data.frame(t(m1$model)), magnitudes = log(rdf$fpm))
            lines(log10(rdf$fpm[!vpi]+1)[fpmo], log10(exp(m1$model[["corr.a"]]*log(rdf$fpm[!vpi])[fpmo]+m1$model[["corr.b"]])+1))
            if("corr.ltheta.b" %in% names(m1$model)) {
                # show 95% CI for the non-constant theta fit
                xval <- range(log(rdf$fpm[!vpi]))
                xval <- seq(xval[1], xval[2], length.out = 100)
                thetas <- get.corr.theta(m1$model, xval)
                #thetas <- exp(m1$model[["corr.ltheta.i"]]+m1$model[["corr.ltheta.lfpm"]]*xval)
                #thetas <- (1+exp((m1$model["corr.ltheta.lfpm.m"] - xval)/m1$model["corr.ltheta.lfpm.s"]))/m1$model["corr.ltheta.a"]
                alpha <- 0.05
                yval <- exp(m1$model[["corr.a"]]*xval + m1$model[["corr.b"]])
                lines(log10(exp(xval)+1), log10(qnbinom(alpha/2, size = thetas, mu = yval)+1), col = 1, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(1-alpha/2, size = thetas, mu = yval)+1), col = 1, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(alpha/2, size = m1$model[["corr.theta"]], mu = yval)+1), col = 8, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(1-alpha/2, size = m1$model[["corr.theta"]], mu = yval)+1), col = 8, lty = 2)
            }
        } else {
            lines(log10(rdf$fpm[!vpi]+1)[fpmo], log10(m1@components[[2]][[2]]@predict(cbind(1, log(rdf$fpm[!vpi])))+1)[fpmo], col = 4)
        }
    }
    legend(x = "topleft", col = c("red", "blue"), pch = 19, legend = c("failure component", "correlated component"), bty = "n", cex = 0.9)

    # zero fit
    if(n.zero.windows > nrow(rdf)) { n.zero.windows <- nrow(rdf) }
    bw <- floor(nrow(rdf)/n.zero.windows)
    if(compressed.models) {
        rdf$cluster <- m1$clusters
    } else {
        rdf$cluster <- clusters(m1)
    }
    rdf <- rdf[order(rdf$fpm, decreasing = FALSE), ]
    fdf <- data.frame(y = rowMeans(matrix(log10(rdf$fpm[1:(n.zero.windows*bw)]+1), ncol = bw, byrow = TRUE)), zf = rowMeans(matrix(as.integer(rdf$cluster[1:(n.zero.windows*bw)] == 1), ncol = bw, byrow = TRUE)))
    plot(zf~y, fdf, ylim = c(0, 1), xlim = range(na.omit(log10(rdf$fpm+1))), xlab = "expected FPM", ylab = "fraction of failures", main = "failure model", pch = 16, cex = 0.5)
    ol <- order(rdf$fpm, decreasing = TRUE)
    if(compressed.models) {
        fp <- scde.failure.probability(data.frame(t(m1$model)), magnitudes = log(rdf$fpm))
        lines(log10(rdf$fpm[ol]+1), fp[ol], col = 2)
    } else {
        mt <- terms(m1@concomitant@formula, data = rdf)
        mf <- model.frame(delete.response(mt), data = rdf, na.action = NULL)
        cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
        cm0 <- cm0/rowSums(cm0)
        lines(log10(rdf$fpm[ol]+1), cm0[ol, 1], col = 2)
    }


    # show thetas
    #tl <- c(fail = get("theta", envir = environment(m1@components[[1]][[1]]@logLik)), corr = get("theta", envir = environment(m1@components[[2]][[2]]@logLik)))
    if(compressed.models) {
        if("corr.ltheta.b" %in% names(m1$model)) {
            p <- exp(m1$model[["corr.a"]]*log(rdf$fpm[!vpi])+m1$model[["corr.b"]])
            alpha <- ((rdf$count[!vpi]/p-1)^2 - 1/p)
            trng <- log(range(c(m1$model[["corr.theta"]], thetas))) + 0.5*c(-1, 1)
            # restrict the alpha to the confines of the estimated theta values
            alpha[alpha > exp(-trng[1])] <- exp(-trng[1])
            alpha[alpha<exp(-trng[2])] <- exp(-trng[2])

            smoothScatter(log10(rdf$fpm[!vpi]+1), -log10(alpha), ylim = trng*log10(exp(1)), xlab = "FPM", ylab = "log10(theta)", main = "overdispersion", bandwidth = bandwidth)
            xval <- range(log(rdf$fpm[!vpi]))
            xval <- seq(xval[1], xval[2], length.out = 100)
            #thetas <- exp(m1$model[["corr.ltheta.i"]]+m1$model[["corr.ltheta.lfpm"]]*xval)
            thetas <- get.corr.theta(m1$model, xval)
            #plot(log10(exp(xval)+1), log(thetas), ylim = trng, type = 'l', xlab = "FPM", ylab = "log(theta)", main = "overdispersion")
            lines(log10(exp(xval)+1), log10(thetas))
            abline(h = log10(m1$model[["corr.theta"]]), col = 1, lty = 2)
        } else {
            tl <- c(fail = c(), corr = m1$model[["corr.theta"]])
            barplot(tl, beside = TRUE, las = 2, col = c("dodgerblue1", "indianred1"), ylab = "magnitude", main = "theta")
        }

    } else {
        tl <- c(fail = c(0), corr = get("theta", envir = environment(m1@components[[2]][[2]]@logLik)))
        barplot(tl, beside = TRUE, las = 2, col = c("dodgerblue1", "indianred1"), ylab = "magnitude", main = "theta")
    }
    box()
    if(do.par) {   dev.off() }
}

## from nb2.crossmodels.r
mc.stepFlexmix <- function(..., nrep = 5, n.cores = nrep, return.all = FALSE) {
    if(nrep < 2) {
        return(flexmix(...))
    } else {
        ml <- papply(seq_len(nrep), function(m) {
            x = try(flexmix(...))
        }, n.cores = n.cores)
        if(return.all) { return(ml) }
        ml <- ml[unlist(lapply(ml, function(x) !is(x, "try-error")))]
        logLiks <- unlist(lapply(ml, logLik))
        ml[[which.max(logLiks)]]
    }
}

# df: count matrix
# xr: expression level for each row in the matrix
# ml: fitted model list for a replicate
get.rep.set.posteriors <- function(xr, df, ml, rescale = TRUE) {
    pl <- do.call(cbind, lapply(seq_along(ml), function(i) {
        edf <- data.frame(y = df[, i], xr = xr)
        m1 <- ml[[i]]
        x <- FLXgetModelmatrix(m1@model[[1]], edf, m1@model[[1]]@formula)
        #cx <- FLXgetModelmatrix(m1@concomitant, edf, m1@concomitant@formula)
        cm <- (1/(1+exp(-1*(x@x %*% m1@concomitant@coef[, -1]))))
        p1 <- (1-cm)*exp(FLXdeterminePostunscaled(x, m1@components[[1]]))
        p2 <- cm*exp(FLXdeterminePostunscaled(x, m1@components[[2]]))
        tpr <- as.numeric(p1+p2)
        if(rescale) {tpr <- tpr/sum(tpr) }
        return(tpr)
    }))
    colnames(pl) <- names(ml)

    return(pl)
}

# evaluates likelihood for a list of models and a set of
# corresponding counts
# ml - model list
# counts - observed count matrix corresponding to the models
# marginals - marginal info, to which model-specific count will be appended
get.rep.set.general.model.posteriors <- function(ml, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, min.p = 0) {
    pl <- do.call(cbind, lapply(seq_along(ml), function(i) {
        marginals$count <- counts[, i]
        rowSums(get.component.model.lik(ml[[i]], marginals))+min.p
    }))
    if(rescale) {
        #pl <- pl*grid.weight+min.p
        pl <- t(t(pl)/colSums(pl))
    }
    colnames(pl) <- names(ml)
    return(pl)
}

get.rep.set.general.model.logposteriors <- function(ml, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE) {
    pl <- do.call(rbind, lapply(seq_along(ml), function(i) {
        marginals$count <- counts[, i]
        log.row.sums(get.component.model.loglik(ml[[i]], marginals))
    }))
    if(rescale) {
        pl <- pl-log.row.sums(pl)
    }
    rownames(pl) <- names(ml)
    return(t(pl))
}

# evaluate likelihood on a mixed model with a binomial concomitant
# returns posterior probability for each component: rowSums(return) gives
# total likelihood. (note it's not on a log scale!)
get.component.model.lik <- function(m1, newdata) {
    # core models
    cp <- exp(do.call("+", lapply(seq_along(m1@model), function(i) {
        y <- posterior(m1@model[[i]], newdata, lapply(m1@components, "[[", i))
    })))
    # concomitant

    # no groups!
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
    cm0 <- cm0/rowSums(cm0)
    cm0[!is.finite(cm0)] <- 1
    return(cp*cm0)
}

# same as above, but keeping log resolution
get.component.model.loglik <- function(m1, newdata) {
    # core models
    cp <- do.call("+", lapply(seq_along(m1@model), function(i) {
        y <- posterior(m1@model[[i]], newdata, lapply(m1@components, "[[", i))
    }))
    cp[!is.finite(cp)] <- sign(cp[!is.finite(cp)])*.Machine$double.xmax
    # concomitant
    # no groups!
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- model.matrix(mt, data = mf) %*% m1@concomitant@coef
    cm0[is.nan(cm0)] <- 1
    cm0[!is.finite(cm0)] <- sign(cm0[!is.finite(cm0)])*.Machine$double.xmax

    cm0 <- cm0-log.row.sums(cm0)
    return(cp+cm0)
}

# returns a matrix of posterior values, with rows corresponding to genes, and
# columns to marginal values (prior fpkm grid)
# m1 - model
# counts - vector of per-gene counts for a given experiment
# marginals - fpm data frame
get.exp.posterior.matrix <- function(m1, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, n.cores = 32, min.p = 0) {
    uc <- unique(counts)
    #message(paste("get.exp.posterior.matrix() :", round((1-length(uc)/length(counts))*100, 3), "% savings"))
    cat(".")
    df <- do.call(rbind, papply(uc, function(x) {
        rowSums(get.component.model.lik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))+min.p
    }, n.cores = n.cores))
    if(rescale) {
        #df <- t(t(df)*grid.weight)+min.p
        df <- df/rowSums(df)
    }
    df <- df[match(counts, uc), , drop = FALSE]
    rownames(df) <- names(counts)
    df
}

get.exp.logposterior.matrix <- function(m1, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, n.cores = 32) {
    uc <- unique(counts)
    #message(paste("get.exp.logposterior.matrix() :", round((1-length(uc)/length(counts))*100, 3), "% savings"))
    cat(".")
    df <- do.call(rbind, papply(uc, function(x) {
        log.row.sums(get.component.model.loglik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))
    }, n.cores = n.cores))
    if(rescale) {
        df <- df-log.row.sums(df)
    }
    df <- df[match(counts, uc), , drop = FALSE]
    rownames(df) <- names(counts)
    df
}

# similar to get.exp.posterior.matrix(), but returns inverse ecdf list
# note that x must be supplied
get.exp.posterior.samples <- function(pmatl, prior, n.samples = 1, n.cores = 32) {
    sl <- papply(seq_along(pmatl), function(i) t(apply(pmatl[[i]], 1, function(d) approxfun(cumsum(d), prior$x, rule = 2)(runif(n.samples)))), n.cores = n.cores)
    names(sl) <- names(pmatl)
    sl
}
# similar to get.exp.posterior.matrix(), but returns inverse ecdf list
# note that x must be supplied
get.exp.sample <- function(m1, counts, marginals, prior.x, n, rescale = TRUE) {
    do.call(rbind, papply(counts, function(x) {
        tpr <- log.row.sums(get.component.model.loglik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))
        if(rescale)  {
            tpr <- exp(tpr-max(tpr))
            tpr <- tpr/sum(tpr)
        }
        return(approxfun(cumsum(tpr), prior.x, rule = 2)(runif(n)))
    }, n.cores=1))
}

# gets a probability of failed detection for a given observation
# optional vector of fpm values (log) can be supplied to evaluate mixing probability
# at a point other than MLE fpm
get.concomitant.prob <- function(m1, counts = NULL, lfpm = NULL) {
    if(is.null(lfpm)) {
        lfpm <- get.fpm.estimates(m1, counts)
    }
    newdata <- data.frame(fpm = exp(lfpm))
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
    cm0[is.nan(cm0)] <- 1
    cm0 <- cm0/rowSums(cm0)
    return(as.numeric(cm0[, 1]))
}

# copied from flexmix
log.row.sums <- function(m) {
    M <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))] # "random" doesn't work!
    M + log(rowSums(exp(m - M)))
}

# variation of negative.binomial family that keeps theta value accessible
negbin.th <- function (theta = stop("'theta' must be specified"), link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt"))
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family available links are \"identity\", \"log\" and \"sqrt\"")
    }
    .Theta <- theta
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/.Theta
    validmu <- function(mu) all(mu  >  0)
    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1,
                                                             y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
            lgamma(.Theta + y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y ==  0)/6
    })
    simfun <- function(object, nsim) {
        ftd <- fitted(object)
        val <- rnegbin(nsim * length(ftd), ftd, .Theta)
    }
    environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
    famname <- paste("Negative Binomial(", format(round(theta,
                                                        4)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                   aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                   validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = theta),
              class = "family")
}


# .fit version of the glm.nb
glm.nb.fit <- function(x, y, weights = rep(1, nobs), control = list(trace = 0, maxit = 20), offset = rep(0, nobs), etastart = NULL, start = NULL, mustart = NULL, init.theta = NULL, link = "log", method = "glm.fit", intercept = TRUE, theta.range = c(0, 1e5), ...) {
    method <- "custom.glm.fit"
    #require(MASS)
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
                                                            y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                                     log(mu + (y ==  0)) - (th + y) * log(th + mu)))
    # link <- substitute(link)

    Call <- match.call()
    control <- do.call("glm.control", control)
    n <- length(y)
    # family for the initial guess
    fam0 <- if (missing(init.theta) | is.null(init.theta))
        do.call("poisson", list(link = link))
    else
        do.call("negative.binomial", list(theta = init.theta, link = link))

    # fit function
    if (!missing(method)) {
        #message(paste("glm.nb.fit: method = ", method))
        if (!exists(method, mode = "function"))
            stop("unimplemented method: ", sQuote(method))
        glm.fitter <- get(method)
    }
    else {
        #message("glm.nb.fit: using default glm.fit")
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
        #glm.fitter <- custom.glm.fit
    }

    if (control$trace  >  1) {
        message("Initial fit:")
    }
    fit <- glm.fitter(x = x, y = y, weights = weights, start = start, etastart = etastart, mustart = mustart, offset = offset, family = fam0, control = control, intercept = intercept)
    class(fit) <- c("glm", "lm")

    mu <- fit$fitted.values
    th <- as.vector(theta.ml(y, mu, sum(weights), weights, limit = control$maxit, trace = control$trace  >  2))
    if(!is.null(theta.range)) {
        if(th<theta.range[1]) {
            if (control$trace  >  1)
                message("adjusting theta from ", signif(th), " to ", signif(theta.range[1]), " to fit the specified range")
            th <- theta.range[1]
        } else if(th > theta.range[2]) {
            if (control$trace  >  1)
                message("adjusting theta from ", signif(th), " to ", signif(theta.range[2]), " to fit the specified range")
            th <- theta.range[2]
        }
    }
    if (control$trace  >  1)
        message("Initial value for theta:", signif(th))
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, y, weights)
    Lm0 <- Lm + 2 * d1
    while ((iter <- iter + 1) <=  control$maxit && (abs(Lm0 - Lm)/d1 + abs(del)/d2)  >  control$epsilon) {
        eta <- g(mu)
        fit <- glm.fitter(x = x, y = y, weights = weights, etastart = eta, offset = offset, family = fam, control = list(maxit = control$maxit*10, epsilon = control$epsilon, trace = control$trace  >  1), intercept = intercept)
        t0 <- th
        th <- theta.ml(y, mu, sum(weights), weights, limit = control$maxit, trace = control$trace  >  2)
        if(!is.null(theta.range)) {
            if(th<theta.range[1]) {
                if (control$trace  >  1)
                    message("adjusting theta from ", signif(th), " to ", signif(theta.range[1]), " to fit the specified range")
                th <- theta.range[1]
            } else if(th > theta.range[2]) {
                if (control$trace  >  1)
                    message("adjusting theta from ", signif(th), " to ", signif(theta.range[2]), " to fit the specified range")
                th <- theta.range[2]
            }
        }
        fam <- do.call("negative.binomial", list(theta = th, link = link))
        mu <- fit$fitted.values
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, y, weights)
        if (control$trace) {
            Ls <- loglik(n, th, y, y, weights)
            Dev <- 2 * (Ls - Lm)
            message("Theta(", iter, ")  = ", signif(th), ", 2(Ls - Lm)  = ",  signif(Dev))
        }
    }
    if (!is.null(attr(th, "warn")))
        fit$th.warn <- attr(th, "warn")
    if (iter  >  control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && intercept) {
        null.deviance <- if ("(Intercept)" %in% colnames(x))
            glm.fitter(x[, "(Intercept)", drop = FALSE], y, weights = weights, offset = offset, family = fam, control = list(maxit = control$maxit*10, epsilon = control$epsilon, trace = control$trace  >   1), intercept = TRUE)$deviance
        else
            fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin.th", "glm", "lm")
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    fit$x <- x
    fit$y <- y
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}

custom.glm.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                            mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                            control = list(), intercept = TRUE, alpha = 0)
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars ==  0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                 call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0L)
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start)  !=  nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                              nvars, paste(deparse(xnames), collapse = ", ")),
                     domain = NA)
        else {
            coefold <- start
            offset + as.vector(if (NCOL(x) ==  1)
                x * start
                else x %*% start)
        }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                 call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights  >  0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu)))
                stop("NAs in V(mu)")
            if (any(varmu ==  0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights  >  0) & (mu.eta.val  !=  0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ",
                        iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            #z <- (eta - offset)[good] + (family$linkfun(y[good]) - family$linkfun(mu[good]))/family$linkfun(mu.eta.val[good])

            # attempting to be robust here, trowing out fraction with highest abs(z)
            if(alpha > 0) {
                qv <- quantile(abs(z), probs = c(alpha/2, 1.0-alpha/2))
                gvi <- which(good)[which(abs(z)<qv[1] | abs(z) > qv[2])]
                good[gvi] <- FALSE
                if (all(!good)) {
                    conv <- FALSE
                    warning("no observations informative at iteration ",
                            iter)
                    break
                }
                z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            }

            w <- mu.eta.val[good]*sqrt(weights[good]/variance(mu)[good])

            ngoodobs <- as.integer(nobs - sum(!good))
            fit <- .Fortran("dqrls", qr = x[good, ] * w, n = ngoodobs,
                            p = nvars, y = w * z, ny = 1L, tol = min(1e-07,
                                                                     control$epsilon/1000), coefficients = double(nvars),
                            residuals = double(ngoodobs), effects = double(ngoodobs),
                            rank = integer(1L), pivot = 1L:nvars, qraux = double(nvars),
                            work = double(2 * nvars)) # , PACKAGE = "base"
            #browser()
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                                 iter), domain = NA)
                break
            }
            if (nobs < fit$rank)
                stop(gettextf("X matrix has rank %d, but only %d observations",
                              fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance  = ", dev, "Iterations -", iter,
                    "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                    stop("no valid set of coefficients has been found: please supply starting values",
                         call. = FALSE)
                warning("step size truncated due to divergence",
                        call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                    if (ii  >  control$maxit)
                        stop("inner loop 1 cannot correct step size",
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            # require deviance to go down
            if ((!is.null(coefold)) & (dev - devold)/(0.1 + abs(dev))  >  3*control$epsilon) {
                warning("step size truncated due to increasing divergence", call. = FALSE)
                ii <- 1
                while ((dev - devold)/(0.1 + abs(dev))  >  3*control$epsilon) {
                    if (ii  >  control$maxit)   {
                        warning("inner loop 1 cannot correct step size", call. = FALSE)
                        break
                    }
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                    stop("no valid set of coefficients has been found: please supply starting values",
                         call. = FALSE)
                warning("step size truncated: out of bounds",
                        call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                    if (ii  >  control$maxit)
                        stop("inner loop 2 cannot correct step size",
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }

            devold <- dev
            coef <- coefold <- start
        }
        if (!conv)
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary)
            warning("glm.fit: algorithm stopped at boundary value",
                    call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family ==  "binomial") {
            if (any(mu  >  1 - eps) || any(mu < eps))
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                        call. = FALSE)
        }
        if (family$family ==  "poisson") {
            if (any(mu < eps))
                warning("glm.fit: fitted rates numerically 0 occurred",
                        call. = FALSE)
        }
        if (fit$rank < nvars)
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat)  >  col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
                                                                    sum(good) - fit$rank))
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights ==  0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, length(y), mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
         effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
         rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
                                                       "qraux", "pivot", "tol")], class = "qr"), family = family,
         linear.predictors = eta, deviance = dev, aic = aic.model,
         null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
         df.residual = resdf, df.null = nulldf, y = y, converged = conv,
         boundary = boundary)
}

# copied from limma
weighted.median.scde <- function (x, w, na.rm = FALSE)
    #       Weighted median
    #       Gordon Smyth
    #       30 June 2005
{
    if (missing(w))
        w <- rep.int(1, length(x))
    else {
        if(length(w)  !=  length(x)) stop("'x' and 'w' must have the same length")
        if(any(is.na(w))) stop("NA weights not allowed")
        if(any(w<0)) stop("Negative weights not allowed")
    }
    if(is.integer(w))
        w <- as.numeric(w)
    if(na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    if(all(w == 0)) {
        warning("All weights are zero")
        return(NA)
    }
    o <- order(x)
    x <- x[o]
    w <- w[o]
    p <- cumsum(w)/sum(w)
    n <- sum(p<0.5)
    if(p[n+1]  >  0.5)
        x[n+1]
    else
        (x[n+1]+x[n+2])/2
}

# given a set of pdfs (columns), calculate summary statistics (mle, 95% CI, Z-score deviations from 0)
quick.distribution.summary <- function(s.bdiffp) {
    diffv <- as.numeric(colnames(s.bdiffp))
    dq <- t(apply(s.bdiffp, 1, function(p) {
        mle <- which.max(p)
        p <- cumsum(p)
        return(diffv[c(lb = max(c(1, which(p<0.025))), mle, min(c(length(p), which(p > (1-0.025)))))])
    }))/log10(2)
    colnames(dq) <- c("lb", "mle", "ub")
    cq <- rep(0, nrow(dq))
    cq[dq[, 1] > 0] <- dq[dq[, 1] > 0, 1]
    cq[dq[, 3]<0] <- dq[dq[, 3]<0, 3]
    z <- get.ratio.posterior.Z.score(s.bdiffp)
    za <- sign(z)*qnorm(p.adjust(pnorm(abs(z), lower.tail = FALSE), method = "BH"), lower.tail = FALSE)
    data.frame(dq, "ce" = as.numeric(cq), "Z" = as.numeric(z), "cZ" = as.numeric(za))
}

