# performs weighted centering of mat rows (mat - rowSums(mat*weights)/rowSums(weights))
# possibly accounting for batch effects (i.e. centering each batch separately
weightedMatCenter <- function(mat, matw, batch = NULL) {
    if(is.null(batch)) {
        return(mat-rowSums(mat*matw)/rowSums(matw))
    } else {
        cmat <- mat
        invisible(tapply(seq_len(ncol(mat)), as.factor(batch), function(ii) {
            cmat[, ii] <<- cmat[, ii, drop = FALSE] - rowSums(cmat[, ii, drop = FALSE]*matw[, ii, drop = FALSE])/rowSums(matw[, ii, drop = FALSE])
        }))
        return(cmat)
    }
}

# per-experiment/per-gene weighted variance estimate
# weight matrix should have the same dimensions as the data matrix
weightedMatVar <- function(mat, matw, batch = NULL, center = TRUE, min.weight = 0, normalize.weights = TRUE) {
    # normalize weights
    #matw <- matw/rowSums(matw)
    #matw <- matw/rowSums(matw)*ncol(matw)
    if(center) {
        mat <- weightedMatCenter(mat, matw, batch)
    }

    #weightedMatVar.Rcpp(mat, matw)
    #return(rowSums(mat*mat*matw) / (1-rowSums(matw*matw)))
    #return(rowSums(mat*mat*matw))

    #return(rowSums(mat*mat*matw) * rowSums(matw) /pmax(rowSums(matw)^2 - rowSums(matw*matw), rep(min.weight, nrow(matw))))

    v<- rowSums(mat*mat*matw)
    if(normalize.weights) { v <- v/rowSums(matw) }
    v
}

# GEV t() function
gev.t <- function(x, loc, scale, shape = rep(0, length(loc)), log = FALSE) {
    if(log) {
        pmin(0, ifelse(shape == 0, -(x-loc)/scale, (-1/shape)*log(pmax(0, 1+shape*(x-loc)/scale))))
    } else {
        pmin(1, ifelse(shape == 0, exp(-(x-loc)/scale), ((pmax(0, 1+shape*(x-loc)/scale))^(-1/shape))))
    }
}
# returns upper tail of GEV in log scale
pgev.upper.log <- function(x, loc, scale, shape = rep(0, length(loc))) {
    tv <- gev.t(x, loc, scale, shape, log = TRUE)
    tv[tv >  -5 & tv<0] <- log(-expm1(-exp(tv[tv >  -5 & tv<0])))
    tv
}

# BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
    nai <- which(!is.na(x))
    ox <- x
    x<-x[nai]
    id <- order(x, decreasing = FALSE)
    if(log) {
        q <- x[id] + log(length(x)/seq_along(x))
    } else {
        q <- x[id]*length(x)/seq_along(x)
    }
    a <- rev(cummin(rev(q)))[order(id)]
    ox[nai]<-a
    ox
}

pathway.pc.correlation.distance <- function(pcc, xv, n.cores = 1, target.ndf = NULL) {
    # all relevant gene names
    rotn <- unique(unlist(lapply(pcc[gsub("^#PC\\d+# ", "", rownames(xv))], function(d) rownames(d$xp$rotation))))
    # prepare an ordered (in terms of genes) and centered version of each component
    pl <- papply(rownames(xv), function(nam) {
        pnam <- gsub("^#PC\\d+# ", "", nam)
        pn <- as.integer(gsub("^#PC(\\d+)# .*", "\\1", nam))
        rt <- pcc[[pnam]]$xp$rotation[, pn]
        # order names/values according to increasing name match index
        mi <- match(names(rt), rotn)
        mo <- order(mi, decreasing = FALSE)
        rt <- as.numeric(rt)-mean(rt)
        return(list(i = mi[mo], v = rt[mo]))
    }, n.cores = n.cores)

    x <- .Call("plSemicompleteCor2", pl, PACKAGE = "scde")

    if(!is.null(target.ndf)) {
        r <- x$r[upper.tri(x$r)]
        n <- x$n[upper.tri(x$n)]
        suppressWarnings(tv <- r*sqrt((n-2)/(1-r^2)))
        z <- pt(tv, df = n-2, lower.tail = FALSE, log.p = TRUE)
        nr <- qt(z, df = target.ndf-2, lower.tail = FALSE, log.p = TRUE)
        nr <- nr/sqrt(target.ndf-2+nr^2)
        nr[is.nan(nr)] <- r[is.nan(nr)]

        cr <- x$r
        cr[upper.tri(cr)] <- nr
        cr[lower.tri(cr)] <- t(cr)[lower.tri(cr)]
    } else {
        cr <- x$r
    }

    rownames(cr) <- colnames(cr) <- rownames(xv)
    d <- stats::as.dist(1-abs(cr))
    d[d<0] <- 0
    d

}

collapse.aspect.clusters <- function(d, dw, ct, scale = TRUE, pick.top = FALSE) {
    xvm <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        if(length(ii) == 1) return(d[ii, ])
        if(pick.top) {
            return(d[ii[which.max(apply(d[ii, ], 1, var))], ])
        }
        xp <- pcaMethods::pca(t(d[ii, ]), nPcs = 1, center = TRUE, scale = "none")
        xv <- pcaMethods::scores(xp)[, 1]
        if(sum(abs(diff(xv))) > 0 && cor(xv, colMeans(d[ii, ]*abs(pcaMethods::loadings(xp)[, 1])))<0) { xv <- -1*xv }
        #set scale at top pathway?
        if(sum(abs(diff(xv))) > 0) {
            if(scale) {
                xv <- xv*sqrt(max(apply(d[ii, ], 1, var)))/sqrt(var(xv))
            }
            if(sum(abs(xv)) == 0) { xv <- abs(rnorm(length(xv), sd = 1e-6)) }
        } else {
            xv <- abs(rnorm(length(xv), sd = 1e-6))
        }
        #xv <- xv/sqrt(length(ii))
        xv
    }))
    rownames(xvm) <- unlist(tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        if(length(ii) == 1) return(rownames(d)[ii])
        return(rownames(d)[ii[which.max(apply(d[ii, ], 1, var))]])
    }))

    xvmw <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        w <- colSums(dw[ii, , drop = FALSE]*apply(d[ii, , drop = FALSE], 1, sd))
        w <- w/sum(w)
    }))

    return(list(d = xvm, w = xvmw))
}
