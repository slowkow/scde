
#######
## from nb1glm.R
#######

# nb2 glm implementation
setClass("FLXMRnb2gam", contains = "FLXMRglm", package = "flexmix")

setClass("FLXcomponentE",
         representation(refitTheta = "function",
                        theta.fit = "ANY"),
         contains = "FLXcomponent", package = "flexmix")

# nb2 glm implementation
setClass("FLXMRnb2glm", contains = "FLXMRglm", package = "flexmix")

FLXMRnb2glm <- function(formula = . ~ .,  offset = NULL, init.theta = NULL, theta.range = c(0, 1e3), ...) {
    #require(MASS)
    family <- "negative.binomial"
    glmrefit <- function(x, y, w) {
        #message("FLXRnb2glm:refit:nb2")
        fit <- c(glm.nb.fit(x, y, weights = w, offset = offset, init.theta = init.theta, theta.range = theta.range),
                 list(call = sys.call(), offset = offset, control = eval(formals(glm.fit)$control), method = "glm.fit")
        )
        fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
        fit$df.residual <- sum(w) - fit$rank
        fit$x <- x
        fit
    }

    z <- new("FLXMRnb2glm", weighted = TRUE, formula = formula,
             name = "FLXMRnb2glm", offset = offset,
             family = family, refit = glmrefit)
    z@preproc.y <- function(x) {
        if (ncol(x)  >  1)
            stop(paste("for the", family, "family y must be univariate"))
        x
    }


    z@defineComponent <- expression({
        predict <- function(x, ...) {
            dotarg = list(...)
            #message("FLXRnb2glm:predict:nb2")
            if("offset" %in% names(dotarg)) offset <- dotarg$offset
            p <- x%*%coef
            if (!is.null(offset)) p <- p + offset
            negative.binomial(theta)$linkinv(p)
        }
        logLik <- function(x, y, ...) {
            r <- dnbinom(y, size = theta, mu = predict(x, ...), log = TRUE)
            #message(paste("FLXRnb2glm:loglik:nb2", theta))
            return(r)
        }

        new("FLXcomponent",
            parameters = list(coef = coef),
            logLik = logLik, predict = predict,
            df = df)
    })

    z@fit <- function(x, y, w){
        #message("FLXRnb2glm:fit:nb2")
        w[y<= 1] <- w[y<= 1]/1e6 # focus the fit on non-failed genes
        fit <- glm.nb.fit(x, y, weights = w, offset = offset, init.theta = init.theta, theta.range = theta.range)
        # an ugly hack to restrict to non-negative slopes
        cf <- coef(fit)
        if(cf[2]<0) { cf <- c(mean(y*w)/sum(w), 0) }
        with(list(coef = cf, df = ncol(x), theta = fit$theta, offset = offset), eval(z@defineComponent))
    }

    return(z)
}

# component-specific version of the nb2glm
# nb2 glm implementation
setClass("FLXMRnb2glmC", representation(vci = "ANY"), contains = "FLXMRnb2glm", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRnb2glmC <- function(... , components = NULL) {
    #require(MASS)
    z <- new("FLXMRnb2glmC", FLXMRnb2glm(...), vci = components)
    z
}

# get values of theta for a given set of models and expression (log-scale) magnitudes
get.corr.theta <- function(model, lfpm, theta.range = NULL) {
    if("corr.ltheta.b" %in% names(model)) {
        #th <- exp(-1*(model[["corr.ltheta.a"]]/(1+exp((model[["corr.ltheta.lfpm.m"]] - lfpm)/model[["corr.ltheta.lfpm.s"]])) + log(model[["corr.ltheta.b"]])))
        th <- exp(-1*(model[["corr.ltheta.b"]]+(model[["corr.ltheta.t"]]-model[["corr.ltheta.b"]])/(1+10^((model[["corr.ltheta.m"]]-lfpm)*model[["corr.ltheta.s"]]))^model[["corr.ltheta.r"]]))
    } else {
        if(length(lfpm) > 1) {
            th <- rep(model[["corr.theta"]], length(lfpm))
        } else {
            th <- model[["corr.theta"]]
        }
    }
    if(!is.null(theta.range)) {
        th[th<theta.range[1]] <- theta.range[1]
        th[th > theta.range[2]] <- theta.range[2]
        th[is.nan(th)] <- theta.range[1]
    }
    th
}

# nb2 implementation with a simple trimmed-mean/median slope, and a gam theta fit
setClass("FLXMRnb2gth", contains = "FLXMRglm", package = "flexmix")

FLXMRnb2gth <- function(formula = . ~ .,  offset = NULL, full.theta.range = c(1e-3, 1e3), theta.fit.range = full.theta.range*c(1e-1, 1e1), theta.fit.sp = c(-1), constant.theta = FALSE, slope.mean.trim = 0.4, alpha.weight.power = 1/2, ...) {
    if(slope.mean.trim<0) { slope.mean.trim <- 0 }
    if(slope.mean.trim > 0.5) { slope.mean.trim <- 0.5 }

    family <- "negative.binomial"
    glmrefit <- function(x, y, w) {
        message("ERROR: FLXRnb2gth:glmrefit: NOT IMPLEMENTED")
        return(NULL)
    }

    z <- new("FLXMRnb2gth", weighted = TRUE, formula = formula,
             name = "FLXMRnb2gth", offset = offset,
             family = family, refit = glmrefit)
    z@preproc.y <- function(x) {
        if (ncol(x)  >  1)
            stop(paste("for the", family, "family y must be univariate"))
        x
    }

    z@defineComponent <- expression({
        predict <- function(x, ...) {
            dotarg = list(...)
            #message("FLXRnb2gth:predict:nb2")
            coef["corr.a"]*x
        }
        logLik <- function(x, y, ...) {
            dotarg = list(...)
            #message("FLXRnb2gth:logLik")
            if(constant.theta) {
                th <- coef["corr.theta"]
            } else {
                #th <- exp(coef["corr.ltheta.i"] + coef["corr.ltheta.lfpm"]*log(x))
                #th <- exp(-1*(coef["corr.ltheta.a"]/(1+exp((coef["corr.ltheta.lfpm.m"] - log(x))/coef["corr.ltheta.lfpm.s"])) + log(coef["corr.ltheta.b"])))
                th <- get.corr.theta(coef, log(x))

            }
            # restrict theta to the pre-defined range
            th[th  >  full.theta.range[2]] <- full.theta.range[2]
            th[th < full.theta.range[1]] <- full.theta.range[1]
            # evaluate NB
            r <- dnbinom(y, size = th, mu = coef["corr.a"]*x, log = TRUE)
        }

        new("FLXcomponent",
            parameters = list(coef = coef, linear = TRUE),
            logLik = logLik, predict = predict, df = df)
    })

    z@fit <- function(x, y, w){
        # message("FLXRnb2gth:fit")

        # estimate slope using weighted trimmed mean
        #w[y == 0] <- w[y == 0]/1e6
        #r <- y/x
        #ro <- order(r)
        ## cumulative weight sum along the ratio order (to figure out where to trim)
        #cs <- cumsum(w[ro])/sum(w)
        #lb <- min(which(cs > slope.mean.trim))
        #ub <- max(which(cs<(1-slope.mean.trim)))
        #ro <- ro[lb:ub]
        ## slope fit
        #a <- weighted.mean(r[ro], w[ro])

        a <- as.numeric(coef(glm(y~0+x, family = poisson(link = "identity"), start = weighted.mean(y/x, w), weights = w))[1])

        # predicted values
        p <- a*x

        #te <- p^2/((y-p)^2 - p) # theta point estimates
        #te[te<theta.fit.range[1]] <- theta.fit.range[1] te[te > theta.fit.range[2]] <- theta.fit.range[2]
        alpha <- ((y/p-1)^2 - 1/p)
        alpha[alpha<1/theta.fit.range[2]] <- 1/theta.fit.range[2]
        alpha[alpha > 1/theta.fit.range[1]] <- 1/theta.fit.range[1]

        #theta <- MASS::theta.ml(y, p, sum(w), w)
        theta <- MASS::theta.md(y, p, sum(w)-1, w)
        theta <- pmin(pmax(theta.fit.range[1], theta), theta.fit.range[2])
        if(constant.theta) {
            v <- c("corr.a" = a, "corr.theta" = theta)
        } else {
            # fit theta linear model in the log space
            #theta.l <- glm(log(te)~log(x), weights = w)
            #ac <- tryCatch( {

            mw <- w*(as.numeric(alpha)^(alpha.weight.power))
            lx <- log(x)
            lx.rng <- range(lx)
            mid.s <- (sum(lx.rng))/2
            low <- log(x) < mid.s
            lalpha <- log(alpha)
            bottom.s <- quantile(lalpha[low], 0.025, na.rm = TRUE)
            top.s <- quantile(lalpha[!low], 0.975, na.rm = TRUE)

            wsr <- function(p, x, y, w = rep(1, length(y))) {
                # borrowing nplr approach here
                bottom <- p[1]
                top <- p[2]
                xmid <- p[3]
                scal <- p[4]
                r <- p[5]
                yfit <- bottom+(top-bottom)/(1+10^((xmid-x)*scal))^r
                #exp(-1*(model[["corr.ltheta.b"]]+(model[["corr.ltheta.t"]]-model[["corr.ltheta.b"]])/(1+10^((model[["corr.ltheta.m"]]-lfpm)*model[["corr.ltheta.s"]]))^model[["corr.ltheta.r"]]))
                residuals <- (y - yfit)^2
                return(sum(w*residuals))
            }

            #po <- nlm(f = wsr, p = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = w*(as.numeric(alpha)^(1/2)))
            #ac <- po$estimate

            #po <- nlminb(objective = wsr, start = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = mw)
            po <- nlminb(objective = wsr, start = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = mw, lower = c(-100, -10, -100, -100, 0.1), upper = c(10, 100, 100, 0, 20))
            ac <- po$par

            #smoothScatter(log(x), log(alpha))
            #p <- ac points(log(x), p[1]+(p[2]-p[1])/(1+10^((p[3]-log(x))*p[4]))^p[5], col = 2, pch = ".")
            #browser()
            #}, error = function(e) {
            #  message("encountered error trying to fit logistic model with guessed parameters.")
            #  # fit with fewer parameters?
            #})

            v <- c(a, theta, ac)
            names(v) <- c("corr.a", "corr.theta", "corr.ltheta.b", "corr.ltheta.t", "corr.ltheta.m", "corr.ltheta.s", "corr.ltheta.r")
        }

        with(list(coef = v, full.theta.range = full.theta.range, df = ncol(x), offset = offset), eval(z@defineComponent))
    }

    return(z)
}

# component-specific version of the nb2gth
# nb2 gam implementation
setClass("FLXMRnb2gthC", representation(vci = "ANY"), contains = "FLXMRnb2gth", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRnb2gthC <- function(... , components = NULL) {
    #require(mgcv)
    z <- new("FLXMRnb2gthC", FLXMRnb2gth(...), vci = components)
    z
}

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRnb2glmC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
})

setMethod("FLXmstep", signature(model = "FLXMRnb2glmC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRnb2glmC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRnb2glm"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRnb2glmC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# same for gth
setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRnb2gthC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2gthC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2gthC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
})

setMethod("FLXmstep", signature(model = "FLXMRnb2gthC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(q1 = list(coefficients = c(1)), coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRnb2gthC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRnb2gth"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRnb2gthC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# component-specific version of the nb2glm
# nb2 glm implementation
setClass("FLXMRglmC", representation(vci = "ANY"), contains = "FLXMRglm", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRglmC <- function(... , components = NULL) {
    #require(MASS)
    z <- new("FLXMRglmC", FLXMRglm(...), vci = components)
    z
}

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRglmC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
    #message("FLXMRnb2glmC:FLXdeterminePostunscaled : ")
    #message(m)
    #browser()
    m
})

setMethod("FLXmstep", signature(model = "FLXMRglmC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRglmC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRglm"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRglmC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# mu-fixed version
setClass("FLXMRglmCf", representation(mu = "numeric"), contains = "FLXMRglmC", package = "flexmix")

FLXMRglmCf <- function(... , family = c("binomial", "poisson"), mu = 0) {
    #require(MASS)
    family <- match.arg(family)
    z <- new("FLXMRglmCf", FLXMRglmC(..., family = family), mu = mu)
    z
}

setMethod("FLXmstep", signature(model = "FLXMRglmCf"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- c(model@mu, rep(0, ncol(model@x)-1))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        list(comp.1)
    }), recursive = FALSE)
})

# a magnitude-weighted version of the FLXPmultinom (to down-weight low-fpkm points during concomitant fit)
# alternatively: some kind of non-decreasing function could be used
setClass("FLXPmultinomW", contains = "FLXPmultinom")

FLXPmultinomW <- function(formula = ~1) {
    z <- new("FLXPmultinom", name = "FLXPmultinom", formula = formula)
    multinom.fit <- function(x, y, w, ...) {
        r <- ncol(x)
        p <- ncol(y)
        if (p < 2) stop("Multinom requires at least two components.")
        mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
        #if(missing(w)) w <- rep(1, nrow(y))
        #w <- round(exp(x[, 2]))
        nnet::nnet.default(x, y, w, mask = mask, size = 0,
                           skip = TRUE, softmax = TRUE, censored = FALSE,
                           rang = 0, trace = FALSE, ...)
    }
    z@fit <- function(x, y, w, ...) multinom.fit(x, y, w, ...)$fitted.values
    z@refit <- function(x, y, w, ...) {
        if (missing(w) || is.null(w)) w <- rep(1, nrow(y))
        #w <- round(exp(x[, 2]))
        fit <- nnet::multinom(y ~ 0 + x, weights = w, data = list(y = y, x = x), Hess = TRUE, trace = FALSE)
        fit$coefnames <- colnames(x)
        fit$vcoefnames <- fit$coefnames[seq_along(fit$coefnames)]
        dimnames(fit$Hessian) <- lapply(dim(fit$Hessian) / ncol(x), function(i) paste(rep(seq_len(i) + 1, each = ncol(x)), colnames(x), sep = ":"))
        fit
    }
    z
}

