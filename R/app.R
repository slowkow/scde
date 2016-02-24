
# rook class for browsing differential expression results

ViewDiff <- setRefClass(
    'ViewDiff',
    fields = c('gt', 'models', 'counts', 'prior', 'groups', 'batch', 'geneLookupURL'),
    methods = list(

        initialize = function(results, models, counts, prior, groups = NULL, batch = NULL, geneLookupURL = NULL) {
            if(!is.null(results$results)) {
                gt <<- results$results
            } else {
                gt <<- results
            }
            # add raw names if this wasn't a batch-corrected sample
            if("mle" %in% colnames(gt)) {
                colnames(gt) <<- paste(colnames(gt), "raw", sep = "_")
            }
            if(!is.null(results$batch.adjusted)) {
                df <- results$batch.adjusted
                colnames(df) <- paste(colnames(df), "cor", sep = "_")
                gt <<- cbind(gt, df)
            }
            if(!is.null(results$batch.effect)) {
                df <- results$batch.effect
                colnames(df) <- paste(colnames(df), "bat", sep = "_")
                gt <<- cbind(gt, df)
            }
            colnames(gt) <<- tolower(colnames(gt))

            # append expression levels to the results
            if(!is.null(results$joint.posteriors)) {
                gt$lev1 <<- log10(as.numeric(colnames(results$joint.posteriors[[1]]))[max.col(results$joint.posteriors[[1]])]+1)
                gt$lev2 <<- log10(as.numeric(colnames(results$joint.posteriors[[2]]))[max.col(results$joint.posteriors[[2]])]+1)
            }
            gt$gene <<- rownames(gt)
            gt <<- data.frame(gt)

            # guess gene lookup for common cases
            if(is.null(geneLookupURL)) {
                # human
                if( any(grepl("ENSG\\d+", gt$gene[1])) || any(c("CCLU1", "C22orf45") %in% gt$gene)) {
                    geneLookupURL <<- "http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g = {0}"
                } else if( any(grepl("ENSMUSG\\d+", gt$gene[1]))) {
                    geneLookupURL <<- "http://useast.ensembl.org/Mus_musculus/Gene/Summary?g = {0}"
                } else if( any(c("Foxp2", "Sept1", "Lrrc34") %in% gt$gene)) {
                    # mouse MGI
                    geneLookupURL <<- "http://www.informatics.jax.org/searchtool/Search.do?query = {0}"
                } else if( any(grepl("FBgn\\d+", gt$gene[1])) || any(c("CG3680", "CG8290") %in% gt$gene)) {
                    # flybase
                    geneLookupURL <<- "http://flybase.org/cgi-bin/uniq.html?db = fbgn&GeneSearch = {0}&context = {1}&species = Dmel&cs = yes&caller = genejump"
                } else {
                    # default, forward to ensemble search
                    geneLookupURL <<- "http://useast.ensembl.org/Multi/Search/Results?q = {0}site = ensembl"
                }
            } else {
                geneLookupURL <<- geneLookupURL
            }



            gt <<- gt[gt$z_raw != "NA", ]
            gt <<- gt[!is.na(gt$z_raw), ]


            models <<- models
            counts <<- counts
            prior <<- prior
            if(is.null(groups)) { # recover groups from models
                groups <<- as.factor(attr(models, "groups"))
                if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
                names(groups) <<- rownames(models)
            } else {
                groups <<- groups
            }
            if(length(levels(groups)) != 2) {
                stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
            }

            batch <<- batch
            callSuper()
        },
        call = function(env){
            path <- env[['PATH_INFO']]
            req <- Request$new(env)
            res <- Response$new()

            switch(path,
                   # INDEX
                   '/index.html' = {
                       body <- paste('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd" >
                                     <html >
                                     <head >
                                     <meta http-equiv = "Content-Type" content = "text/html charset = iso-8859-1" >
                                     <title > SCDE: ', paste(levels(groups), collapse = " vs. "), '</title >
                                     <!-- ExtJS -- >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/resources/ext-theme-neptune/ext-theme-neptune-all.css" / >

                                     <!-- Shared -- >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/ext-4.2.1.883/examples/shared/example.css" / >

                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/additional.css" / >
                                     <!-- GC -- >

                                     <style type = "text/css" >
                                     .x-panel-framed {
                                     padding: 0
                                     }
                                     </style >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/ext-4.2.1.883/ext-all.js" > </script >

                                     <script type = "text/javascript" > var geneLookupURL = "', geneLookupURL, '"</script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/viewembed.js" > </script >

                                     </head >
                                     <body style = "margin-top:0padding-top:10px" >
                                     <div id = "example-grid" > </div >
                                     </body >
                                     </html >
                                     ', sep = "")
                       res$header('"Content-Type": "text/html"')
                       res$write(body)
                   },
                   # GENE TABLE
                   '/genetable.json' = {
                       lgt <- gt
                       if(!is.null(req$params()$filter)) {
                           fl <- rjson::fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               if(is.null(lgt$z_cor)) { lgt <- lgt[order(abs(lgt$z_raw), decreasing = TRUE), ] } else { lgt <- lgt[order(abs(lgt$z_cor), decreasing = TRUE), ] }
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- rjson::toJSON(list(totalCount = trows, genes = ol))
                       res$header('"Content-Type": "application/json"')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   # POSTERIOR PLOT
                   '/posterior.png' = {
                       gene <- ifelse(is.null(req$params()$gene), sample(gt$gene), req$params()$gene)
                       bootstrap <- ifelse(is.null(req$params()$bootstrap), TRUE, req$params()$bootstrap == "T")
                       show.individual.posteriors <- ifelse(is.null(req$params()$show.individual.posteriors), TRUE, req$params()$show.individual.posteriors == "true")

                       t <- tempfile()
                       #require(Cairo)
                       CairoPNG(filename = TRUE, width = 350, height = 560)
                       scde.test.gene.expression.difference(gene = gene, models = models, counts = counts, groups = groups, prior = prior, batch = batch, ratio.range = c(-10, 10), show.individual.posteriors = show.individual.posteriors, verbose = FALSE)
                       dev.off()
                       res$header('Content-type', 'image/png')
                       res$body <- t
                       names(res$body) <- 'file'
                   },
                   # GENE EXPRESSION LEVELS
                   '/elevels.html' = {
                       geneName <- ifelse(is.null(req$params()$geneName), gt$gene[[1]], req$params()$geneName)
                       gc <- counts[rownames(counts) == geneName, , drop = FALSE]
                       fpm <- exp(scde.expression.magnitude(models, counts = gc))
                       df <- rbind(FPM = gc, level = fpm)
                       df <- round(df, 2)
                       # order columns according to groups
                       df <- df[, unlist(tapply(seq_len(ncol(df)), groups, I))]
                       cell.col <- rep(c("#E9A994", "#66CCFF"), as.integer(table(groups)))

                       render.row <- function(nam, val, col) {
                           paste("<tr > ", "<th > ", nam, "</th > ", paste("<td bgcolor = ", col, " > ", val, "</td > ", sep = "", collapse = " "), "</tr > ", sep = "")
                       }

                       sh <- paste("<tr > ", paste("<th > ", c(" ", colnames(df)), "</th > ", sep = "", collapse = " "), "</tr > ")
                       #sb <- paste(render.row("cells", colnames(df), cell.col), render.row("FPKM", df[1, ], cell.col), render.row("mode", df[2, ], cell.col), collapse = "\n")
                       sb <- paste(render.row("counts", df[1, ], cell.col), render.row("FPM", df[2, ], cell.col), collapse = "\n")
                       res$header('"Content-Type": "text/html"')
                       res$write(paste("<table id = \"elevels\" > ", sh, sb, "</table > "))
                   },
{
    res$write('default')
}
                       )
            res$finish()
        }
            )
    )

t.view.pathways <- function(pathways, mat, matw, env, proper.names = rownames(mat), colcols = NULL, zlim = NULL, labRow = NA, vhc = NULL, cexCol = 1, cexRow = 1, n.pc = 1, nstarts = 50, row.order = NULL, show.Colv = TRUE, plot = TRUE, trim = 1.1/ncol(mat), bwpca = TRUE, ...) {
    # retrieve gis
    lab <- which(proper.names %in% na.omit(unlist(mget(pathways, envir = env, ifnotfound = NA))))

    if(length(lab) == 0) {
        # try genes
        lab <- which(proper.names %in% pathways)
    }
    if(length(lab) == 0)
        return(NULL)
    #t.quick.show.mat(mat[lab, ], normalize.rows = TRUE)
    #table(rownames(mat) %in% mget(pathways, envir = env))

    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    d <- mat[lab, , drop = FALSE]
    dw <- matw[lab, , drop = FALSE]
    if(length(lab) > 2) {
        xp <- bwpca(t(d), t(dw), npcs = n.pc, center = FALSE, nstarts = 3)
    } else {
        xp <- list()
    }
    hc <- NULL;
    d <- d-rowMeans(d)
    dd <- as.dist(1-abs(cor(t(as.matrix(d)))))
    dd[is.na(dd)] <- 1
    if(is.null(row.order)) {
        if(length(lab) > 2) {
            if(is.element("fastcluster", installed.packages()[, 1])) {
                hc <- fastcluster::hclust(dd, method = "ward.D")
            } else {
                hc <- stats::hclust(dd, method = "ward.D")
            }
            row.order <- hc$order
        } else {
            row.order <- c(seq_along(lab))
            # make up a fake hc
            if(length(lab)>1) {
              hc<-list();
              attributes(hc)<-list(members=length(lab),height=1);
              class(hc)<-"dendrogram";
              hc[[1]] <- list();
              attributes(hc[[1]]) <- list(members=1,height=0,label=rownames(mat)[lab[1]],leaf=T)
              hc[[2]] <- list();
              attributes(hc[[2]]) <- list(members=1,height=0,label=rownames(mat)[lab[2]],leaf=T)
            } else { 
              hc <- list(); attributes(hc) <- list(members=1,height=0,label=rownames(mat)[lab[1]],leaf=T); class(hc) <- "dendrogram";
            }
        }
    }

    if(is.null(vhc)) {
        vd <- as.dist(1-cor(as.matrix(d)))
        vd[is.na(vd)] <- 1
        if(is.element("fastcluster", installed.packages()[, 1])) {
            vhc <- fastcluster::hclust(vd, method = "ward.D")
        } else {
            vhc <- stats::hclust(vd, method = "ward.D")
        }

    }

    if(is.null(zlim)) { zlim <- quantile(d, p = c(0.01, 0.99)) }
    vmap <- d
    vmap[vmap<zlim[1]] <- zlim[1]
    vmap[vmap > zlim[2]] <- zlim[2]
    rownames(vmap) <- rownames(d)

    aval <- colSums(d*dw)/colSums(dw)
    if(!is.null(xp$scores)) {
        oc <- xp$scores[, n.pc]
        if(cor(oc, aval, method = "spearman")<0) {
            oc <- -1*oc
            xp$scores[, n.pc] <- -1*xp$scores[, n.pc]
            xp$rotation[, n.pc] <- -1*xp$rotation[, n.pc]
        }
        xp$oc <- oc
        z <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc/max(abs(oc))*49)+50]
        ld <- xp$rotation[row.order, , drop = FALSE]
        ld <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(ld/max(abs(ld))*49)+50]
    } else {
        ld <- z <- NULL
    }

    aval <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(100)[round(aval/max(abs(aval))*49)+50]



    # oc2 <- xp$scores[, 2]
    # oc2 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc2/max(abs(oc2))*49)+50]
    # oc3 <- xp$scores[, 3]
    # oc3 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc3/max(abs(oc3))*49)+50]


    #z <- do.call(rbind, list(aval, oc))
    #z <- rbind(oc3, oc2, oc)

    col <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(256)
    if(plot) {
        if(!is.null(colcols)) { z <- rbind(colcols, z) }
        if(show.Colv) { Colv <- as.dendrogram(vhc) } else { Colv <- NA }
        my.heatmap2(vmap[row.order, , drop = FALSE], Rowv = NA, Colv = Colv, zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z, labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
    }
    xp$vhc <- vhc
    xp$lab <- lab
    xp$hc <- as.dendrogram(hc)
    xp$row.order <- row.order
    xp$col <- col
    xp$oc.col <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(256)
    xp$vmap <- vmap
    xp$zlim <- zlim
    return(invisible(xp))
}

##' View pathway or gene weighted PCA
##'
##' Takes in a list of pathways (or a list of genes), runs weighted PCA, optionally showing the result.
##' @param pathways character vector of pathway or gene names
##' @param varinfo output of pagoda.varnorm()
##' @param goenv environment mapping pathways to genes
##' @param n.genes number of genes to show
##' @param two.sided whether the set of shown genes should be split among highest and lowest loading (T) or if genes with highest absolute loading (F) should be shown
##' @param n.pc optional integer vector giving the number of principal component to show for each listed pathway
##' @param colcols optional column color matrix
##' @param zlim optional z color limit
##' @param showRowLabels controls whether row labels are shown in the plot
##' @param cexCol column label size (cex)
##' @param cexRow row label size (cex)
##' @param nstarts number of random starts for the wPCA
##' @param cell.clustering cell clustering
##' @param show.cell.dendrogram whether cell dendrogram should be shown
##' @param plot whether the plot should be shown
##' @param box whether to draw a box around the plotted matrix
##' @param trim optional Winsorization trim that should be applied
##' @param return.details whether the function should return the matrix as well as full PCA info instead of just PC1 vector
##' @param ... additional arguments are passed to the \code{c.view.pathways}
##' @return cell scores along the first principal component of shown genes (returned as invisible)
##' @export
pagoda.show.pathways <- function(pathways, varinfo, goenv = NULL, n.genes = 20, two.sided = FALSE, n.pc = rep(1, length(pathways)), colcols = NULL, zlim = NULL, showRowLabels = FALSE, cexCol = 1, cexRow = 1, nstarts = 10, cell.clustering = NULL, show.cell.dendrogram = TRUE, plot = TRUE, box = TRUE, trim = 0, return.details = FALSE , ...) {
    labRow <- NA
    if(showRowLabels) { labRow <- NULL }
    x <- c.view.pathways(pathways, varinfo$mat, varinfo$matw, goenv, batch = varinfo$batch, n.genes = n.genes, two.sided = two.sided, n.pc = n.pc, colcols = colcols, zlim = zlim, labRow = labRow, cexCol = cexCol, cexRow = cexRow, trim = trim, show.Colv = show.cell.dendrogram, plot = plot, vhc = cell.clustering, labCol = NA, box = TRUE, ...)
    if(return.details) {
        invisible(x)
    } else {
        invisible(x$scores[, 1])
    }
}

# takes in a list of pathways with a list of corresponding PC numbers
# recalculates PCs for each individual pathway, weighting gene loading in each pathway and then by total
# pathway variance over the number of genes (rough approximation)
c.view.pathways <- function(pathways, mat, matw, goenv = NULL, batch = NULL, n.genes = 20, two.sided = TRUE, n.pc = rep(1, length(pathways)), colcols = NULL, zlim = NULL, labRow = NA, vhc = NULL, cexCol = 1, cexRow = 1, nstarts = 50, row.order = NULL, show.Colv = TRUE, plot = TRUE, trim = 1.1/ncol(mat), showPC = TRUE,  ...) {
    # are these genes or pathways being passed?
    if(!is.null(goenv)) {
        x <- pathways %in% ls(goenv)
    } else {
        x <- rep(FALSE, length(pathways))
    }
    if(sum(x) > 0) { # some pathways matched
      if(!all(x)) {
        message("WARNING: partial match to pathway names. The following entries did not match: ", paste(pathways[!x], collapse = " "))
        }
        # look up genes for each pathway
        pathways <- pathways[x]
        p.genes <- mget(pathways, goenv, ifnotfound = NA)
    } else { # try as genes
        x <- pathways %in% rownames(mat)
        if(sum(x) > 0) {
            if(!all(x)) {
                message("WARNING: partial match to gene names. The following entries did not match: ", paste(pathways[!x], collapse = " "))
            }
            p.genes <- list("genes" = pathways[x])
            pathways <- c("genes");
        } else { # neither genes nor pathways are passed
            stop("ERROR: provided names do not match either gene nor pathway names (if the pathway environment was provided)")
        }
    }
    gvi <- rownames(mat) %in% unlist(p.genes)
    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    # recalculate wPCA for each pathway
    ppca <- pagoda.pathway.wPCA(varinfo = list(mat = mat[gvi, , drop = FALSE], matw = matw[gvi, , drop = FALSE], batch = batch), setenv = list2env(p.genes), n.cores = 1, n.randomizations = 0, n.starts = 2, n.components = max(n.pc), verbose = FALSE, min.pathway.size = 0, max.pathway.size = Inf, n.internal.shuffles = 0)

    if(length(ppca) > 1) { # if more than one pathway was supplied, combine genes using appropriate loadings and use consensus PCA (1st PC) as a pattern
        # score top loading genes for each desired PC, scaling by the sd/sqrt(n) (so that ^2 is = var/n)
        scaled.gene.loadings <- unlist(lapply(seq_along(pathways), function(i) {
            gl <- ppca[[pathways[i]]]$xp$rotation[, n.pc[i], drop = TRUE]*as.numeric(ppca[[pathways[i]]]$xp$sd)[n.pc[i]]/sqrt(ppca[[pathways[i]]]$n)
            names(gl) <- rownames(ppca[[pathways[i]]]$xp$rotation)
            gl
        }))


        if(two.sided) {
            # positive
            reduced.gene.loadings <- sort(tapply(scaled.gene.loadings, as.factor(names(scaled.gene.loadings)), max), decreasing = TRUE)
            selected.genes.pos <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), round(n.genes/2))]

            # negative
            reduced.gene.loadings <- sort(tapply(scaled.gene.loadings, as.factor(names(scaled.gene.loadings)), min), decreasing = FALSE)
            selected.genes.neg <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), round(n.genes/2))]
            selected.genes <- c(selected.genes.pos, selected.genes.neg)
            selected.genes <- selected.genes[match(unique(names(selected.genes)), names(selected.genes))]

        } else {
            reduced.gene.loadings <- sort(tapply(abs(scaled.gene.loadings), as.factor(names(scaled.gene.loadings)), max), decreasing = TRUE)
            selected.genes <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), n.genes)]
        }

        # consensus pattern
        #lab <- match(names(selected.genes), rownames(mat))
        lab <- names(selected.genes);

        if(length(lab) == 0)
            return(NULL)
        if(length(lab)<3) { return(NULL) }
        if(trim > 0) {
            rn <- rownames(mat)
            cn <- colnames(mat)
            mat <- winsorize.matrix(mat, trim = trim)
            rownames(mat) <- rn
            colnames(mat) <- cn
        }
        d <- mat[lab, , drop = FALSE]
        dw <- matw[lab, , drop = FALSE]

        #d <- d*abs(as.numeric(selected.genes))
        xp <- bwpca(t(d), t(dw), npcs = 1, center = FALSE, nstarts = 3)

        consensus.npc = 1 # use first PC as a pattern
    } else { # only one pathway was provided
        xp <- ppca[[1]]$xp
        lab <- rownames(xp$rotation)[order(abs(xp$rotation[, n.pc[1]]), decreasing = TRUE)]
        if(length(lab) > n.genes) {
            lab <- lab[1:n.genes]
            xp$rotation <- xp$rotation[lab, , drop = FALSE]
        }

        d <- mat[lab, , drop = FALSE]
        dw <- matw[lab, , drop = FALSE]
        consensus.npc = n.pc[1] # use specified PC as a pattern
    }

    d <- d-rowMeans(d)
    dd <- as.dist(1-abs(cor(t(as.matrix(d)))))
    dd[is.na(dd)] <- 1
    if(is.null(row.order)) {
        if(length(lab) > 2) {
            if(is.element("fastcluster", installed.packages()[, 1])) {
                hc <- fastcluster::hclust(dd, method = "ward.D")
            } else {
                hc <- stats::hclust(dd, method = "ward.D")
            }
            row.order <- hc$order
        } else {
            row.order <- c(seq_along(lab))
            if(length(lab)>1) {
              hc<-list();
              attributes(hc)<-list(members=length(lab),height=1);
              class(hc)<-"dendrogram";
              hc[[1]] <- list();
              attributes(hc[[1]]) <- list(members=1,height=0,label=lab[1],leaf=T)
              hc[[2]] <- list();
              attributes(hc[[2]]) <- list(members=1,height=0,label=lab[2],leaf=T)
            } else { 
              hc <- list(); attributes(hc) <- list(members=1,height=0,label=lab[1],leaf=T); class(hc) <- "dendrogram";
            }
        }
    }

    if(is.null(vhc)) {
        vd <- as.dist(1-cor(as.matrix(d)))
        vd[is.na(vd)] <- 1
        if(is.element("fastcluster", installed.packages()[, 1])) {
            vhc <- fastcluster::hclust(vd, method = "ward.D")
        } else {
            vhc <- stats::hclust(vd, method = "ward.D")
        }
    }

    if(is.null(zlim)) { zlim <- quantile(d, p = c(0.01, 0.99)) }
    vmap <- d
    vmap[vmap<zlim[1]] <- zlim[1]
    vmap[vmap > zlim[2]] <- zlim[2]
    rownames(vmap) <- rownames(d)

    aval <- colSums(d*dw*as.numeric(abs(xp$rotation[, consensus.npc])))/colSums(dw)
    oc <- xp$scores[, consensus.npc]
    if(cor(oc, aval, method = "p")<0) {
        oc <- -1*oc
        xp$scores[, consensus.npc] <- -1*xp$scores[, consensus.npc]
        xp$rotation[, consensus.npc] <- -1*xp$rotation[, consensus.npc]
    }

    aval <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(100)[round(aval/max(abs(aval))*49)+50]
    z <- rbind(colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc/max(abs(oc))*49)+50])

    ld <- xp$rotation[lab[row.order], consensus.npc]
    ld <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(ld/max(abs(ld))*49)+50]

    # oc2 <- xp$scores[, 2]
    # oc2 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc2/max(abs(oc2))*49)+50]
    # oc3 <- xp$scores[, 3]
    # oc3 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc3/max(abs(oc3))*49)+50]

    #z <- do.call(rbind, list(aval, oc))
    #z <- rbind(oc3, oc2, oc)
    if((!showPC) || length(lab)<= 1) {
        z <- NULL
    }
    col <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(256)

    if(!is.null(colcols)) {
      if(is.null(z)) {
        z <- colcols;
      } else {
        z <- rbind(colcols, z)
      }
    }

    if(plot) {
        if(show.Colv) {
            my.heatmap2(vmap[row.order, , drop = FALSE], Rowv = NA, Colv = as.dendrogram(vhc), zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z, labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
        } else {
            my.heatmap2(vmap[row.order, vhc$order, drop = FALSE], Rowv = NA, Colv = NA, zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z[,vhc$order], labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
        }

    }
    xp$vhc <- vhc
    xp$lab <- lab
    xp$hc <- as.dendrogram(hc)
    xp$row.order <- row.order
    xp$oc <- oc
    xp$col <- col
    xp$oc.col <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(256)
    xp$vmap <- vmap
    xp$zlim <- zlim
    xp$consensus.pc <- consensus.npc
    return(invisible(xp))
}

# returns enriched categories for a given gene list as compared with a given universe
# returns a list with over and under fields containing list of over and underrepresented terms
calculate.go.enrichment <- function(genelist, universe, pvalue.cutoff = 1e-3, mingenes = 3, env = go.env, subset = NULL, list.genes = FALSE, over.only = FALSE) {
    genelist <- unique(genelist)
    all.genes <- unique(ls(env))
    # determine sizes
    universe <- unique(c(universe, genelist))
    universe <- universe[universe != ""]
    genelist <- genelist[genelist != ""]
    ns <- length(intersect(genelist, all.genes))
    us <- length(intersect(universe, all.genes))
    #pv <- lapply(go.map, function(gl) { nwb <- length(intersect(universe, gl[[1]])) if(nwb<mingenes) { return(0.5)} else { p <- phyper(length(intersect(genelist, gl[[1]])), nwb, us-nwb, ns) return(ifelse(p > 0.5, 1.0-p, p)) }})

    # compile count vectors
    stab <- table(unlist(mget(as.character(genelist), env, ifnotfound = NA), recursive = TRUE))
    utab <- table(unlist(mget(as.character(universe), env, ifnotfound = NA), recursive = TRUE))
    if(!is.null(subset)) {
        stab <- stab[names(stab) %in% subset]
        utab <- utab[names(utab) %in% subset]
    }

    tabmap <- match(rownames(stab), rownames(utab))

    cv <- data.frame(cbind(utab, rep(0, length(utab))))
    names(cv) <- c("u", "s")
    cv$s[match(rownames(stab), rownames(utab))] <- as.vector(stab)
    cv <- na.omit(cv)
    cv <- cv[cv$u > mingenes, ]

    if(over.only) {
        lpr <- phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE)
    } else {
        pv <- phyper(cv$s, cv$u, us-cv$u, ns, lower.tail = FALSE)
        lpr <- ifelse(pv<0.5, phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE), phyper(cv$s+1, cv$u, us-cv$u, ns, lower.tail = TRUE, log.p = TRUE))
    }
    lpr <- phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE)
    lpra <- bh.adjust(lpr, log = TRUE)
    z <- qnorm(lpr, lower.tail = FALSE, log.p = TRUE)
    za <- qnorm(lpra, lower.tail = FALSE, log.p = TRUE)
    # correct for multiple hypothesis
    mg <- length(which(cv$u > mingenes))
    if(over.only) {
        if(pvalue.cutoff<1) {
            ovi <- which(lpra<= log(pvalue.cutoff))
            uvi <- c()
        } else {
            ovi <- which((lpr+mg)<= log(pvalue.cutoff))
            uvi <- c()
        }
    } else {
        if(pvalue.cutoff<1) {
            ovi <- which(pv<0.5 & lpra<= log(pvalue.cutoff))
            uvi <- which(pv > 0.5 & lpra<= log(pvalue.cutoff))
        } else {
            ovi <- which(pv<0.5 & (lpr+mg)<= log(pvalue.cutoff))
            uvi <- which(pv > 0.5 & (lpr+mg)<= log(pvalue.cutoff))
        }
    }
    ovi <- ovi[order(lpr[ovi])]
    uvi <- uvi[order(lpr[uvi])]

    #return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], p = pr[ovi]*mg), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], p = pr[uvi]*mg)))
    if(list.genes) {
        x <- mget(as.character(genelist), env, ifnotfound = NA)
        df <- data.frame(id = rep(names(x), unlist(lapply(x, function(d) length(na.omit(d))))), go = na.omit(unlist(x)), stringsAsFactors = FALSE)
        ggl <- tapply(df$id, as.factor(df$go), I)
        ovg <- as.character(unlist(lapply(ggl[rownames(cv)[ovi]], paste, collapse = " ")))
        uvg <- as.character(unlist(lapply(ggl[rownames(cv)[uvi]], paste, collapse = " ")))
        return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], Za = za, fe = cv$s[ovi]/(ns*cv$u[ovi]/us), genes = ovg), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], Za = za, fe = cv$s[uvi]/(ns*cv$u[uvi]/us), genes = uvg)))
    } else {
        return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], p.raw = exp(lpr[ovi]), fdr = exp(lpra)[ovi], Z = z[ovi], Za = za[ovi], fe = cv$s[ovi]/(ns*cv$u[ovi]/us), fer = cv$s[ovi]/(length(genelist)*cv$u[ovi]/length(universe))), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], p.raw = exp(lpr[uvi]), fdr = exp(lpra)[uvi], Z = z[uvi], Za = za[uvi], fe = cv$s[uvi]/(ns*cv$u[uvi]/us))))
    }
}

##' A Reference Class to represent the PAGODA application
##'
##' This ROOK application class enables communication with the client-side ExtJS framework and Inchlib HTML5 canvas libraries to create the graphical user interface for PAGODA
##' Refer to the code in \code{\link{make.pagoda.app}} for usage example
##'
##' @field results Output of the pathway clustering and redundancy reduction
##' @field genes List of genes to display in the Detailed clustering panel
##' @field pathways
##' @field mat Matrix of posterior mode count estimates
##' @field matw Matrix of weights associated with each estimate in \code{mat}
##' @field goenv Gene set list as an environment
##' @field renv Global environment
##' @field name Name of the application page; for display as the page title
##' @field trim Trim quantity used for Winsorization for visualization
##' @field batch Any batch or other known confounders to be included in the visualization as a column color track
##'
ViewPagodaAppOld <- setRefClass(
    'ViewPagodaAppOld',
    fields = c('results', 'genes', 'pathways', 'mat', 'matw', 'goenv', 'renv', 'name', 'trim', 'batch'),
    methods = list(

        initialize = function(results, pathways, genes, mat, matw, goenv, batch = NULL, name = "pathway overdispersion", trim = 1.1/ncol(mat)) {
            results <<- results
            #results$tvc$order <<- rev(results$tvc$order);
            results$tvc$labels <<- as.character(1:nrow(results$rcm));
            rownames(results$rcm) <<- as.character(1:nrow(results$rcm));
            genes <<- genes
            genes$svar <<- genes$var/max(genes$var)
            genes <<- genes
            mat <<- mat
            matw <<- matw
            batch <<- batch
            goenv <<- goenv
            pathways <<- pathways
            name <<- name
            trim <<- trim
            # reverse lookup environment
            renvt <- new.env(parent = globalenv())
            xn <- ls(envir = goenv)
            xl <- mget(xn, envir = goenv)
            gel <- tapply(rep(xn, unlist(lapply(xl, length))), unlist(xl), I)
            gel <- gel[nchar(names(gel)) > 0]
            x <- lapply(names(gel), function(n) assign(n, gel[[n]], envir = renvt))
            renv <<- renvt
            rm(xn, xl, x, gel, renvt)
            gc()
            callSuper()
        },
        getgenecldata = function(genes = NULL, gcl = NULL, ltrim = 0) { # helper function to get the heatmap data for a given set of genes
            if(is.null(gcl)) {
                gcl <- t.view.pathways(genes, mat = mat, matw = matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim)
            }

            matrix <- gcl$vmap[rev(gcl$row.order), results$hvc$order, drop = FALSE]
            matrix <- list(data = as.numeric(t(matrix)),
                           dim = dim(matrix),
                           rows = rownames(matrix),
                           cols = colnames(matrix),
                           colors = gcl$col,
                           domain = seq.int(gcl$zlim[1], gcl$zlim[2], length.out = length(gcl$col))
            )

            ol <- list(matrix = matrix)
            if(nrow(gcl$vmap) > 2) {
                rcmvar <- matrix(gcl$rotation[rev(gcl$row.order), , drop = FALSE], ncol = 1)
                rowcols <- list(data = as.numeric(t(rcmvar)),
                                dim = dim(rcmvar),
                                colors = gcl$oc.col,
                                domain = seq.int(-1*max(abs(rcmvar)), max(abs(rcmvar)), length.out = length(gcl$oc.col))
                )

                colcols <- matrix(gcl$oc[results$hvc$order], nrow = 1)
                colcols <- list(data = as.numeric(t(colcols)),
                                dim = dim(colcols),
                                colors = gcl$oc.col,
                                domain = seq.int(-1*max(abs(colcols)), max(abs(colcols)), length.out = length(gcl$oc.col))
                )
                ol <- c(ol, list(rowcols = rowcols, colcols = colcols))
            }
            ol
        },
        call = function(env){
            path <- env[['PATH_INFO']]
            req <- Request$new(env)
            res <- Response$new()
            switch(path,
                   # INDEX
                   '/index.html' = {
                       body <- paste('<!DOCTYPE html >
                                     <meta charset = "utf-8" >
                                     <html >
                                     <head >
                                     <title > ', name, '</title >
                                     <meta http-equiv = "Content-Type" content = "text/html charset = iso-8859-1" >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/resources/ext-theme-neptune/ext-theme-neptune-all.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/examples/shared/example.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/pathcl.css" / >
                                     <head profile = "http://www.w3.org/2005/10/profile" >
                                     <link rel = "icon" type = "image/png" href = "http://pklab.med.harvard.edu/sde/pagoda.png" >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/extjs/ext-all.js" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/jquery-1.11.1.min.js" > </script >
                                     <script src = "http://d3js.org/d3.v3.min.js" charset = "utf-8" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/pathcl.js" > </script >
                                     </head >
                                     <body > </body >
                                     </html >
                                     ', sep = "")
                       res$header('"Content-Type": "text/html"')
                       res$write(body)
                   },
                   '/pathcl.json' = { # report pathway clustering heatmap data
                       # column dendrogram
                       t <- paste(tempfile(), "svg", sep = ".")
                       svg(file = t, width = 1, height = 1) # will be rescaled later
                       par(mar = rep(0, 4), mgp = c(2, 0.65, 0), cex = 1, oma = rep(0, 4))
                       #plot(results$hvc, main = "", sub = "", xlab = "", ylab = "", axes = FALSE, labels = FALSE, xaxs = "i", yaxs = "i", hang = 0.02)
                       plot(as.dendrogram(results$hvc), axes = FALSE, yaxs = "i", xaxs = "i", xlab = "", ylab = "", sub = "", main = "", leaflab = "none")
                       dev.off()
                       x <- readLines(t)
                       treeg <- paste(x[-c(1, 2, length(x))], collapse = "")

                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       matrix <- list(data = as.numeric(t(matrix)),
                                      dim = dim(matrix),
                                      rows = rownames(matrix),
                                      cols = colnames(matrix),
                                      colors = results$cols,
                                      domain = seq.int(results$zlim2[1], results$zlim2[2], length.out = length(results$cols)),
                                      range = range(matrix)
                       )


                       icols <- colorRampPalette(c("white", "black"), space = "Lab")(256)
                       rcmvar <- matrix(apply(results$rcm[rev(results$tvc$order), , drop = FALSE], 1, var), ncol = 1)
                       rowcols <- list(data = as.numeric(t(rcmvar)),
                                       # TODO: add annotation
                                       dim = dim(rcmvar),
                                       colors = icols,
                                       domain = seq.int(0, max(rcmvar), length.out = length(icols))
                       )
                       colcols <- list(data = unlist(lapply(as.character(t(results$colcol[nrow(results$colcol):1, results$hvc$order, drop = FALSE])), col2hex)),
                                       dim = dim(results$colcol)
                       )
                       ol <- list(matrix = matrix, rowcols = rowcols, colcols = colcols, coldend = treeg, trim = trim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genecl.json' = { # report heatmap data for a selected set of genes
                       selgenes <- fromJSON(req$POST()$genes)
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       ol <- getgenecldata(selgenes, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/pathwaygenes.json' = { # report heatmap data for a selected set of pathways
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       pws <- fromJSON(req$POST()$genes)
                       n.pcs <- as.integer(gsub("^#PC(\\d+)# .*", "\\1", pws))
                       n.pcs[is.na(n.pcs)]<-1
                       x <- c.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, goenv = goenv, n.pc = n.pcs, n.genes = ngenes, two.sided = twosided, vhc = results$hvc, plot = FALSE, trim = ltrim, batch = batch)
                       #x <- t.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim, n.pc = 1)
                       ##rsc <- as.vector(rowSums(matw[rownames(x$rotation), ]))*x$rotation[, 1]
                       #rsc <- x$rotation[, 1]
                       #if(twosided) {
                       #  extgenes <- unique(c(names(sort(rsc))[1:min(length(rsc), round(ngenes/2))], names(rev(sort(rsc)))[1:min(length(rsc), round(ngenes/2))]))
                       # } else {
                       #   extgenes <- names(sort(abs(rsc), decreasing = TRUE))[1:min(length(rsc), ngenes)]
                       #}
                       #ol <- getgenecldata(extgenes, ltrim = ltrim)
                       ol <- getgenecldata(genes = NULL, gcl = x, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/patterngenes.json' = { # report heatmap of genes most closely matching a given pattern
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       pat <- fromJSON(req$POST()$pattern)
                       # reorder the pattern back according to column clustering
                       pat[results$hvc$order] <- pat
                       patc <- .Call("matCorr", as.matrix(t(mat)), as.matrix(pat, ncol = 1) , PACKAGE = "scde")
                       if(twosided) { patc <- abs(patc) }
                       mgenes <- rownames(mat)[order(as.numeric(patc), decreasing = TRUE)[1:ngenes]]
                       ol <- getgenecldata(mgenes, ltrim = ltrim)
                       ol$pattern <- pat
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/clinfo.json' = {
                       pathcl <- ifelse(is.null(req$params()$pathcl), 1, as.integer(req$params()$pathcl))
                       ii <- which(results$ct == pathcl)
                       tpi <- order(results$matvar[ii], decreasing = TRUE)
                       #tpi <- tpi[seq(1, min(length(tpi), 15))]
                       npc <- gsub("^#PC(\\d+)#.*", "\\1", names(ii[tpi]))
                       nams <- gsub("^#PC\\d+# ", "", names(ii[tpi]))
                       if(exists("myGOTERM", envir = globalenv())) {
                           tpn <- paste(nams, mget(nams, get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           tpn <- nams;
                       }

                       lgt <- data.frame(do.call(rbind, lapply(seq_along(tpn), function(i) c(id = names(ii[tpi[i]]), name = tpn[i], npc = npc[i], od = as.numeric(results$matvar[ii[tpi[i]]])/max(results$matvar), sign = as.numeric(results$matrcmcor[ii[tpi[i]]]), initsel = as.integer(results$matvar[ii[tpi[i]]] >= results$matvar[ii[tpi[1]]]*0.8)))))

                       # process additional filters
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 100, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ]
                       lgt$od <- format(lgt$od, nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))

                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genes.json' = {
                       lgt <- genes
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/pathways.json' = {
                       lgt <- pathways
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/testenr.json' = { # run an enrichment test
                       selgenes <- fromJSON(req$POST()$genes)
                       lgt <- calculate.go.enrichment(selgenes, rownames(mat), pvalue.cutoff = 0.99, env = renv, over.only = TRUE)$over
                       if(exists("myGOTERM", envir = globalenv())) {
                           lgt$nam <- paste(lgt$t, mget(as.character(lgt$t), get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           lgt$name <- lgt$t
                       }
                       lgt <- data.frame(id = paste("#PC1#", lgt$t), name = lgt$nam, o = lgt$o, u = lgt$u, Z = lgt$Z, Za = lgt$Za, fe = lgt$fe, stringsAsFactors = FALSE)

                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/celltable.txt' = {
                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       body <- paste(capture.output(write.table(round(matrix, 1), sep = "\t")), collapse = "\n")
                       res$header('Content-Type', 'text/plain')
                       #res$header('"Content-disposition": attachment')
                       res$write(body)
                   },
                   {
                     res$header('Location', 'index.html')
                     res$write('Redirecting to <a href = "index.html" > index.html</a >  for interactive browsing.')
                   }
                   )
            res$finish()
        }
    )
)


ViewPagodaApp <- setRefClass(
    'ViewPagodaApp',
    fields = c('results', 'genes', 'pathways', 'mat', 'matw', 'goenv', 'renv', 'name', 'trim', 'batch'),
    methods = list(

        initialize = function(results, pathways, genes, mat, matw, goenv, batch = NULL, name = "pathway overdispersion", trim = 1.1/ncol(mat)) {
            results <<- results
            #results$tvc$order <<- rev(results$tvc$order);
            results$tvc$labels <<- as.character(1:nrow(results$rcm));
            rownames(results$rcm) <<- as.character(1:nrow(results$rcm));
            genes <<- genes
            genes$svar <<- genes$var/max(genes$var)
            genes <<- genes
            mat <<- mat
            matw <<- matw
            batch <<- batch
            goenv <<- goenv
            pathways <<- pathways
            name <<- name
            trim <<- trim
            # reverse lookup environment
            renvt <- new.env(parent = globalenv())
            xn <- ls(envir = goenv)
            xl <- mget(xn, envir = goenv)
            gel <- tapply(rep(xn, unlist(lapply(xl, length))), unlist(xl), I)
            gel <- gel[nchar(names(gel)) > 0]
            x <- lapply(names(gel), function(n) assign(n, gel[[n]], envir = renvt))
            renv <<- renvt
            rm(xn, xl, x, gel, renvt)
            gc()
            callSuper()
        },
        getgenecldata = function(genes = NULL, gcl = NULL, ltrim = 0) { # helper function to get the heatmap data for a given set of genes
            if(is.null(gcl)) {
                gcl <- t.view.pathways(genes, mat = mat, matw = matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim)
            }

            matrix <- gcl$vmap[rev(gcl$row.order), results$hvc$order, drop = FALSE]
            matrix <- list(data = as.numeric(t(matrix)),
                           dim = dim(matrix),
                           rows = rownames(matrix),
                           cols = colnames(matrix),
                           colors = gcl$col,
                           zlim = as.numeric(gcl$zlim)
                           )
                                      
            ol <- list(matrix = matrix)
            if(nrow(gcl$vmap) > 2) {
                rcmvar <- matrix(gcl$rotation[rev(gcl$row.order), , drop = FALSE], ncol = 1)
                rowcols <- list(data = as.numeric(t(rcmvar)),
                                dim = dim(rcmvar),
                                colors = gcl$oc.col,
                                zlim = c(-1,1)*max(abs(rcmvar))
                )

                colcols <- matrix(gcl$oc[results$hvc$order], nrow = 1)
                colcols <- list(data = as.numeric(t(colcols)),
                                dim = dim(colcols),
                                colors = gcl$oc.col,
                                zlim = c(-1,1)*max(abs(colcols))
                )
                ol <- c(ol, list(rowcols = rowcols, colcols = colcols))
            }
            ol
        },
        call = function(env){
            path <- env[['PATH_INFO']]
            req <- Request$new(env)
            res <- Response$new()
            switch(path,
                   # INDEX
                   '/index.html' = {
                       body <- paste('<!DOCTYPE html >
                                     <meta charset = "utf-8" >
                                     <html >
                                     <head >
                                     <title > ', name, '</title >
                                     <meta http-equiv = "Content-Type" content = "text/html charset = iso-8859-1" >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/resources/ext-theme-neptune/ext-theme-neptune-all.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/examples/shared/example.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/pathcl.css" / >
                                     <head profile = "http://www.w3.org/2005/10/profile" >
                                     <link rel = "icon" type = "image/png" href = "http://pklab.med.harvard.edu/sde/pagoda.png" >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/extjs/ext-all.js" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/jquery-1.11.1.min.js" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/pathcl_canvas.js" > </script >
                                     </head >
                                     <body > </body >
                                     </html >
                                     ', sep = "")
                       res$header('"Content-Type": "text/html"')
                       res$write(body)
                   },
                   '/pathcl.json' = { # report pathway clustering heatmap data
                       # column dendrogram
                     treeg <- list(merge=as.vector(t(results$hvc$merge)),height=results$hvc$height,order=match(1:length(results$hvc$order),results$hvc$order))

                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       matrix <- list(data = as.numeric(t(matrix)),
                                      dim = dim(matrix),
                                      rows = rownames(matrix),
                                      cols = colnames(matrix),
                                      colors = results$cols,
                                      zlim = as.numeric(results$zlim2)
                       )


                       icols <- colorRampPalette(c("white", "black"), space = "Lab")(256)
                       rcmvar <- matrix(apply(results$rcm[rev(results$tvc$order), , drop = FALSE], 1, var), ncol = 1)
                       rowcols <- list(data = as.numeric(t(rcmvar)),
                                       dim = dim(rcmvar),
                                       colors = icols,
                                       zlim = c(0, max(rcmvar))
                       )
                       colcols <- list(data = unlist(lapply(as.character(t(results$colcol[nrow(results$colcol):1, results$hvc$order, drop = FALSE])), col2hex)),
                                       dim = dim(results$colcol),
                                       rows=rev(rownames(results$colcol))
                       )
                       ol <- list(matrix = matrix, rowcols = rowcols, colcols = colcols, coldend = treeg, trim = trim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genecl.json' = { # report heatmap data for a selected set of genes
                       selgenes <- fromJSON(req$POST()$genes)
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       ol <- getgenecldata(selgenes, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/pathwaygenes.json' = { # report heatmap data for a selected set of pathways
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       pws <- fromJSON(req$POST()$genes)
                       n.pcs <- as.integer(gsub("^#PC(\\d+)# .*", "\\1", pws))
                       n.pcs[is.na(n.pcs)]<-1
                       x <- c.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, goenv = goenv, n.pc = n.pcs, n.genes = ngenes, two.sided = twosided, vhc = results$hvc, plot = FALSE, trim = ltrim, batch = batch)
                       #x <- t.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim, n.pc = 1)
                       ##rsc <- as.vector(rowSums(matw[rownames(x$rotation), ]))*x$rotation[, 1]
                       #rsc <- x$rotation[, 1]
                       #if(twosided) {
                       #  extgenes <- unique(c(names(sort(rsc))[1:min(length(rsc), round(ngenes/2))], names(rev(sort(rsc)))[1:min(length(rsc), round(ngenes/2))]))
                       # } else {
                       #   extgenes <- names(sort(abs(rsc), decreasing = TRUE))[1:min(length(rsc), ngenes)]
                       #}
                       #ol <- getgenecldata(extgenes, ltrim = ltrim)
                       ol <- getgenecldata(genes = NULL, gcl = x, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/patterngenes.json' = { # report heatmap of genes most closely matching a given pattern
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 0/ncol(mat), as.numeric(req$params()$trim))
                       pat <- fromJSON(req$POST()$pattern)
                       # reorder the pattern back according to column clustering
                       pat[results$hvc$order] <- pat
                       patc <- .Call("matCorr", as.matrix(t(mat)), as.matrix(pat, ncol = 1) , PACKAGE = "scde")
                       if(twosided) { patc <- abs(patc) }
                       mgenes <- rownames(mat)[order(as.numeric(patc), decreasing = TRUE)[1:ngenes]]
                       ol <- getgenecldata(mgenes, ltrim = ltrim)
                       ol$pattern <- pat
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/clinfo.json' = {
                       pathcl <- ifelse(is.null(req$params()$pathcl), 1, as.integer(req$params()$pathcl))
                       ii <- which(results$ct == pathcl)
                       tpi <- order(results$matvar[ii], decreasing = TRUE)
                       #tpi <- tpi[seq(1, min(length(tpi), 15))]
                       npc <- gsub("^#PC(\\d+)#.*", "\\1", names(ii[tpi]))
                       nams <- gsub("^#PC\\d+# ", "", names(ii[tpi]))
                       if(exists("myGOTERM", envir = globalenv())) {
                           tpn <- paste(nams, mget(nams, get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           tpn <- nams;
                       }

                       lgt <- data.frame(do.call(rbind, lapply(seq_along(tpn), function(i) c(id = names(ii[tpi[i]]), name = tpn[i], npc = npc[i], od = as.numeric(results$matvar[ii[tpi[i]]])/max(results$matvar), sign = as.numeric(results$matrcmcor[ii[tpi[i]]]), initsel = as.integer(results$matvar[ii[tpi[i]]] >= results$matvar[ii[tpi[1]]]*0.8)))))

                       # process additional filters
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 100, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ]
                       lgt$od <- format(lgt$od, nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))

                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genes.json' = {
                       lgt <- genes
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/pathways.json' = {
                       lgt <- pathways
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/testenr.json' = { # run an enrichment test
                       selgenes <- fromJSON(req$POST()$genes)
                       lgt <- calculate.go.enrichment(selgenes, rownames(mat), pvalue.cutoff = 0.99, env = renv, over.only = TRUE)$over
                       if(exists("myGOTERM", envir = globalenv())) {
                           lgt$nam <- paste(lgt$t, mget(as.character(lgt$t), get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           lgt$name <- lgt$t
                       }
                       lgt <- data.frame(id = paste("#PC1#", lgt$t), name = lgt$nam, o = lgt$o, u = lgt$u, Z = lgt$Z, Za = lgt$Za, fe = lgt$fe, stringsAsFactors = FALSE)

                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/celltable.txt' = {
                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       body <- paste(capture.output(write.table(round(matrix, 1), sep = "\t")), collapse = "\n")
                       res$header('Content-Type', 'text/plain')
                       #res$header('"Content-disposition": attachment')
                       res$write(body)
                   },
                   {
                     res$header('Location', 'index.html')
                     res$write('Redirecting to <a href = "index.html" > index.html</a >  for interactive browsing.')
                   }
                   )
            res$finish()
        }
    )
)

