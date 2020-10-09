#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
limmaDE <- function(formula=NULL, data=NULL, pheno=NULL, ...){
                    call <- match.call()
                    formula <- formula(call)
                    terms <- terms(formula)
                    design <- model.matrix(formula, data=pheno)
                    fit <- limma_fit(object=data, design=design)
                    asgn <- fit$assign[fit$qr$pivot[1L:fit$rank]]
                    uasgn <- unique(asgn)
                    nterms <- length(uasgn)
                    rdf <- fit$df.residual
                    nmeffect <- c('(intercept)', attr(terms, 'term.labels'))
                    resid <- as.matrix(fit$residuals)
                    effects <- as.matrix(fit$effects)
                    ss.list <- vector('list', nrow(effects))
                    df.list <- vector('list', nrow(effects))
                    pval.list <- vector('list', nrow(effects))
                    names(ss.list) <- names(df.list) <- names(pval.list) <- rownames(effects)
                    for(i in 1L:nrow(effects)){
                        roweff <- effects[i, ]
                        rowres <- resid[i, ]
                        if(!is.null(roweff)) roweff <- as.matrix(roweff)[seq_along(asgn), , drop=FALSE]
                        if(!is.null(rowres)) rowres <- as.matrix(rowres)
                        df <- ss <- ms <- numeric()
                        for(j in seq(nterms)){
                            if(is.null(roweff)){
                               df <- c(df, NA)
                               ss <- c(ss, NA)
                               }
                            else{aj <- (asgn==uasgn[j])
                            df <- c(df, sum(aj))
                            ss <- c(ss, sum(roweff[aj, 1]^2))
                                }
                            }
                        if(rdf[i]>0L){
                           df <- c(df, rdf[i])
                           ss <- c(ss, sum(rowres[ , 1]^2))
                           }
                        nt <- length(df)
                        ms <- ifelse(df>0L, ss/df, NA)
                        if(rdf[i]>0L){
                           Fstat <- ms/ms[nt]
                           Fpval <- pf(Fstat, df, rdf[i], lower.tail=FALSE)
                           Fstat[nt] <- Fpval[nt] <- NA
                           }
                        ss.list[[i]] <- ss
                        pval.list[[i]] <- Fpval
                        }
                    SS <- as.data.frame(ss.list)
                    DF <- df
                    Fpval <- as.data.frame(pval.list)
                    rownames(SS) <- rownames(Fpval) <- names(DF) <- c(nmeffect, 'residuals')
                    SS <- t(SS)
                    Fpval <- t(Fpval)
                    Fpval <- Fpval[ , -nt]
                    fit$sum.squares <- SS
                    fit$df <- DF
                    fit$F.pvalue <- Fpval
                    fit$effects <- fit$residuals <- NULL
                    fit
                    }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
limma_fit <- function(object, design=NULL, ndups=1, spacing=1, block=NULL, correlation, weights=NULL, method='ls', ...){
                      y <- extEAWP(object)
                      if(is.null(design)) design <- y$design
                      if(is.null(design)) design <- matrix(1, ncol(y$exprs), 1)
                      else{design <- as.matrix(design)
                      if(mode(design)!='numeric') stop('design must be a numeric matrix')
                      if(nrow(design)!=ncol(y$exprs)) stop("row dimension of design doesn't match column dimension of data object")
                      }
                      ne <- nonEst(design)
                      if(!is.null(ne)) cat('coefficients not estimable:', paste(ne, collapse=' '), '\n')
                      if(missing(ndups) && !is.null(y$printer$ndups)) ndups <- y$printer$ndups
                      if(missing(spacing) && !is.null(y$printer$spacing)) spacing <- y$printer$spacing
                      if(missing(weights) && !is.null(y$weights)) weights <- y$weights
                      method <- match.arg(method, c('ls', 'robust'))
                      if(ndups>1){
                         if(!is.null(y$probes)) y$probes <- uniquegenelist(y$probes, ndups=ndups, spacing=spacing)
                         if(!is.null(y$Amean)) y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), ndups=ndups, spacing=spacing), na.rm=TRUE)
                         }
                      if(method=='robust') fit <- limma_mrlm(y$exprs, design=design, ndups=ndups, spacing=spacing, weights=weights, ...)
                      else if(ndups<2 && is.null(block)) fit <- limma_lm(y$exprs, design=design, ndups=ndups, spacing=spacing, weights=weights)
                      else{if(missing(correlation)) stop('the correlation must be set, see duplicateCorrelation')
                           fit <- limma_gls(y$exprs, design=design, ndups=ndups, spacing=spacing, block=block, correlation=correlation, weights=weights, ...)
                           }
                      if(NCOL(fit$coefficients)>1){
                         n <- rowSums(is.na(fit$coefficients))
                         n <- sum(n>0 & n<NCOL(fit$coefficients))
                         if(n>0) warning('partial NA coefficients for ', n, ' probe(s)', call.=FALSE)
                         }
                      fit$genes <- y$probes
                      fit$Amean <- y$Amean
                      fit$method <- method
                      fit$design <- design
                      new('MArrayLM', fit)
                      }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
limma_lm <- function(M, design=NULL, ndups=1, spacing=1, weights=NULL){
                     M <- as.matrix(M)
                     narrays <- ncol(M)
                     if(is.null(design)) design <- matrix(1, narrays, 1)
                     else design <- as.matrix(design)
                     nbeta <- ncol(design)
                     coef.names <- colnames(design)
                     if(is.null(coef.names)) coef.names <- paste('x', 1:nbeta, sep='')
                     if(!is.null(weights)){
                        weights <- asMatrixWeights(weights, dim(M))
                        weights[weights<=0] <- NA
                        M[!is.finite(weights)] <- NA
                        }
                     if(ndups>1){
                        M <- unwrapdups(M, ndups=ndups, spacing=spacing)
                        design <- design %x% rep_len(1, ndups)
                        if(!is.null(weights)) weights <- unwrapdups(weights, ndups=ndups, spacing=spacing)
                        }
                     ngenes <- nrow(M)
                     stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames=list(rownames(M), coef.names))
                     NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, 'arrayweights')))
                     if(NoProbeWts){
                        if(is.null(weights)) fit <- lm.fit(design, t(M))
                        else{fit <- lm.wfit(design, t(M), weights[1, ])
                        fit$weights <- NULL
                        }
                        if(fit$df.residual>0){
                           if(is.matrix(fit$effects)) fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank+1):narrays, , drop=FALSE]^2))
                           else fit$sigma <- sqrt(mean(fit$effects[(fit$rank+1):narrays]^2))
                           }
                        else fit$sigma <- rep_len(NA_real_, ngenes)
                        fit$effects <- t(fit$effects)
                        fit$residuals <- t(fit$residuals)
                        fit$fitted.values <- NULL
                        fit$coefficients <- t(fit$coefficients)
                        fit$cov.coefficients <- chol2inv(fit$qr$qr, size=fit$qr$rank)
                        est <- fit$qr$pivot[1:fit$qr$rank]
                        dimnames(fit$cov.coefficients) <- list(coef.names[est], coef.names[est])
                        stdev.unscaled[ , est] <- matrix(sqrt(diag(fit$cov.coefficients)), ngenes, fit$qr$rank, byrow=TRUE)
                        fit$stdev.unscaled <- stdev.unscaled
                        fit$df.residual <- rep_len(fit$df.residual, length.out=ngenes)
                        dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
                        fit$pivot <- fit$qr$pivot
                        return(fit)
                        }
                     beta <- stdev.unscaled
                     sigma <- rep_len(NA_real_, ngenes)
                     df.residual <- rep_len(0, ngenes)
                     for(i in 1:ngenes){
                         y <- as.vector(M[i, ])
                         obs <- is.finite(y)
                         if(sum(obs)>0){
                            X <- design[obs, , drop=FALSE]
                            y <- y[obs]
                            if(is.null(weights)) out <- lm.fit(X, y)
                            else{w <- as.vector(weights[i, obs])
                                 out <- lm.wfit(X, y, w)
                                 }
                            est <- !is.na(out$coefficients)
                            beta[i, ] <- out$coefficients
                            stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, size=out$rank)))
                            df.residual[i] <- out$df.residual
                            if(df.residual[i]>0) sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
                            }
                         }
                     QR <- qr(design)
                     cov.coef <- chol2inv(QR$qr, size=QR$rank)
                     est <- QR$pivot[1:QR$rank]
                     dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
                     list(coefficients=beta, stdev.unscaled=stdev.unscaled, sigma=sigma,
                          df.residual=df.residual, cov.coefficients=cov.coef, pivot=QR$pivot, rank=QR$rank)
                     }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
limma_mrlm <- function(M, design=NULL, ndups=1, spacing=1, weights=NULL, ...){
                       if(!requireNamespace('MASS', quietly=TRUE)) stop("MASS package required but is not installed (or can't be loaded)")
                       M <- as.matrix(M)
                       narrays <- ncol(M)
                       samples.names <- colnames(M)
                       if(is.null(design)) design <- matrix(1, narrays, 1)
                       design <- as.matrix(design)
                       coef.names <- colnames(design)
                       nbeta <- ncol(design)
                       if(!is.null(weights)){
                          weights <- asMatrixWeights(weights, dim(M))
                          weights[weights<=0] <- NA
                          M[!is.finite(weights)] <- NA
                          }
                       if(ndups>1){
                          M <- unwrapdups(M, ndups=ndups, spacing=spacing)
                          design <- design %x% rep_len(1, ndups)
                          if(!is.null(weights)) weights <- unwrapdups(weights, ndups=ndups, spacing=spacing)
                          }
                       ngenes <- nrow(M)
                       stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames=list(rownames(M), coef.names))
                       eff <- res <- matrix(NA, ngenes, narrays, dimnames=list(rownames(M), samples.names))
                       sigma <- rep_len(NA_real_, ngenes)
                       df.residual <- rep_len(0, ngenes)
                       for(i in 1:ngenes){
                           y <- as.vector(M[i, ])
                           obs <- is.finite(y)
                           X <- design[obs, , drop=FALSE]
                           y <- y[obs]
                           if(is.null(weights)) w <- rep_len(1, length(y))
                           else w <- as.vector(weights[i, obs])
                           if(length(y)>nbeta){
                              out <- MASS::rlm(x=X, y=y, weights=w, ...)
                              beta[i, ] <- coef(out)
                              stdev.unscaled[i, ] <- sqrt(diag(chol2inv(out$qr$qr)))
                              eff[i, ] <- out$effects
                              res[i, ] <- out$residuals
                              df.residual[i] <- length(y)-out$rank
                              if(df.residual[i]>0) sigma[i] <- out$s
                              }
                           }
                       QR <- qr(design)
                       cov.coef <- chol2inv(QR$qr, size=QR$rank)
                       est <- QR$pivot[1:QR$rank]
                       dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
                       list(coefficients=beta, stdev.unscaled=stdev.unscaled, sigma=sigma,
                            df.residual=df.residual, cov.coefficients=cov.coef,
                            pivot=QR$pivot, rank=QR$rank, effects=eff, residuals=res)
                       }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
limma_gls <- function(M, design=NULL, ndups=2, spacing=1, block=NULL, correlation=NULL, weights=NULL, ...){
                      M <- as.matrix(M)
                      ngenes <- nrow(M)
                      narrays <- ncol(M)
                      if(is.null(design)) design <- matrix(1, narrays, 1)
                      design <- as.matrix(design)
                      if(nrow(design)!=narrays) stop('number of rows of design matrix does not match number of arrays')
                      nbeta <- ncol(design)
                      coef.names <- colnames(design)
                      if(is.null(correlation)) correlation <- duplicateCorrelation(M, design=design,
                                                                                   ndups=ndups,
                                                                                   spacing=spacing,
                                                                                   block=block,
                                                                                   weights=weights, ...)$consensus.correlation
                      if(abs(correlation)>=1) stop('correlation is 1 or -1, so the model is degenerate')
                      if(!is.null(weights)){
                         weights[is.na(weights)] <- 0
                         weights <- asMatrixWeights(weights, dim(M))
                         M[weights<1e-15] <- NA
                         weights[weights<1e-15] <- NA
                         }
                      if(is.null(block)){
                         if(ndups<2){
                            warning('no duplicates (ndups<2)')
                            ndups <- 1
                            correlation <- 0
                            }
                        cormatrix <- diag(rep_len(correlation, narrays), nrow=narrays, ncol=narrays) %x% array(1, c(ndups, ndups))
                        if(is.null(spacing)) spacing <- 1
                        M <- unwrapdups(M, ndups=ndups, spacing=spacing)
                        if(!is.null(weights)) weights <- unwrapdups(weights, ndups=ndups, spacing=spacing)
                        design <- design %x% rep_len(1, ndups)
                        colnames(design) <- coef.names
                        ngenes <- nrow(M)
                        narrays <- ncol(M)
                        }
                      else{ndups <- spacing <- 1
                      block <- as.vector(block)
                      if(length(block)!=narrays) stop('length of block does not match number of arrays')
                      ub <- unique(block)
                      nblocks <- length(ub)
                      Z <- matrix(block, narrays, nblocks)==matrix(ub, narrays, nblocks, byrow=TRUE)
                      cormatrix <- Z %*% (correlation*t(Z))
                      }
                      diag(cormatrix) <- 1
                      stdev.unscaled <- matrix(NA, ngenes, nbeta, dimnames=list(rownames(M), coef.names))
                      NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, 'arrayweights')))
                      if(NoProbeWts){
                         V <- cormatrix
                         if(!is.null(weights)){
                            wrs <- 1/sqrt(weights[1, ])
                            V <- wrs*t(wrs*t(V))
                            }
                         cholV <- chol(V)
                         y <- backsolve(cholV, t(M), transpose=TRUE)
                         dimnames(y) <- rev(dimnames(M))
                         X <- backsolve(cholV, design, transpose=TRUE)
                         dimnames(X) <- dimnames(design)
                         fit <- lm.fit(X, y)
                         if(fit$df.residual>0){
                            if(is.matrix(fit$effects)) fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank), , drop=FALSE]^2))
                            else fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
                            }
                         else fit$sigma <- rep_len(NA_real_, ngenes)
                         fit$effects <- t(fit$effects)
                         fit$residuals <- t(fit$residuals)
                         fit$fitted.values <- NULL
                         fit$coefficients <- t(fit$coefficients)
                         fit$cov.coefficients <- chol2inv(fit$qr$qr, size=fit$qr$rank)
                         est <- fit$qr$pivot[1:fit$qr$rank]
                         dimnames(fit$cov.coefficients) <- list(coef.names[est], coef.names[est])
                         stdev.unscaled[ , est] <- matrix(sqrt(diag(fit$cov.coefficients)), ngenes, fit$qr$rank, byrow=TRUE)
                         fit$stdev.unscaled <- stdev.unscaled
                         fit$df.residual <- rep_len(fit$df.residual, length.out=ngenes)
                         dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
                         fit$pivot <- fit$qr$pivot
                         fit$ndups <- ndups
                         fit$spacing <- spacing
                         fit$block <- block
                         fit$correlation <- correlation
                         return(fit)
                         }
                      beta <- stdev.unscaled
                      sigma <- rep_len(NA_real_, ngenes)
                      df.residual <- rep_len(0, ngenes)
                      for(i in 1:ngenes){
                          y <- drop(M[i, ])
                          o <- is.finite(y)
                          y <- y[o]
                          n <- length(y)
                          if(n>0){
                             X <- design[o, , drop=FALSE]
                             V <- cormatrix[o, o]
                             if(!is.null(weights)){
                                wrs <- 1/sqrt(drop(weights[i, o]))
                                V <- wrs*t(wrs*t(V))
                                }
                             cholV <- chol(V)
                             y <- backsolve(cholV, y, transpose=TRUE)
                             if(all(X==0)){
                                df.residual[i] <- n
                                sigma[i] <- sqrt(array(1/n, c(1, n)) %*% y^2)
                                }
                             else{X <- backsolve(cholV, X, transpose=TRUE)
                                  out <- lm.fit(X, y)
                                  est <- !is.na(out$coefficients)
                                  beta[i, ] <- out$coefficients
                                  stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, size=out$rank)))
                                  df.residual[i] <- out$df.residual
                                  if(df.residual[i]>0) sigma[i] <- sqrt(array(1/out$df.residual, c(1, n)) %*% out$residuals^2)
                                  }
                             }
                          }
                      cholV <- chol(cormatrix)
                      QR <- qr(backsolve(cholV, design, transpose=TRUE))
                      cov.coef <- chol2inv(QR$qr, size=QR$rank)
                      est <- QR$pivot[1:QR$rank]
                      dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
                      list(coefficients=beta, stdev.unscaled=stdev.unscaled, sigma=sigma,
                           df.residual=df.residual, ndups=ndups, spacing=spacing, block=block,
                           correlation=correlation, cov.coefficients=cov.coef, pivot=QR$pivot, rank=QR$rank)
                      }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
extEAWP <- function(object){
                    y <- list()
                    if(is(object, 'list')){
                       if(is(object, 'EList')){
                          y$exprs <- as.matrix(object$E)
                          y$Amean <- rowMeans(y$exprs, na.rm=TRUE)
                          }
                      else{if(is(object, 'EListRaw')) stop('EListRaw object: please normalize first')
                           if(is(object, 'RGList')) stop('RGList object: please normalize first')
                           y$printer <- object$printer
                           if(is.null(object$M)) stop("data object isn't of a recognized data class")
                           y$exprs <- as.matrix(object$M)
                           if(!is.null(object$A)) y$Amean <- rowMeans(as.matrix(object$A), na.rm=TRUE)
                           }
                      y$weights <- object$weights
                      y$probes <- object$genes
                      y$design <- object$design
                      }
                    else{if(is(object, 'ExpressionSet')){
                            if(!requireNamespace('Biobase', quietly=TRUE)) stop('Biobase package required but is not available')
                            y$exprs <- Biobase::exprs(object)
                            if(length(object@featureData@data)) y$probes <- object@featureData@data
                            y$Amean <- rowMeans(y$exprs, na.rm=TRUE)
                            if('weights' %in% Biobase::assayDataElementNames(object)) y$weights <- Biobase::assayDataElement(object, 'weights')
                            }
                         else{if(is(object, 'PLMset')){
                                 y$exprs <- object@chip.coefs
                                 if(length(y$exprs)==0) stop('chip.coefs has length zero')
                                 if(length(object@se.chip.coefs)) y$weights <- 1/pmax(object@se.chip.coefs, 1e-05)^2
                                 y$Amean <- rowMeans(y$exprs, na.rm=TRUE)
                                 }
                              else{if(is(object, 'marrayNorm')){
                                      y$exprs <- object@maM
                                      if(length(object@maW)) y$weights <- object@maW
                                      if(length(object@maGnames@maInfo)){
                                         y$probes <- object@maGnames@maInfo
                                         attr(y$probes, 'Notes') <- object@maGnames@maNotes
                                         }
                                      if(length(object@maA)) y$Amean <- rowMeans(object@maA, na.rm=TRUE)
                                      }
                                   else{if(is(object, 'eSet')){
                                           if(!requireNamespace('Biobase', quietly=TRUE)) stop('Biobase package required but is not available')
                                           y$exprs <- Biobase::assayDataElement(object, 'exprs')
                                           if(length(object@featureData@data)) y$probes <- object@featureData@data
                                           y$Amean <- rowMeans(y$exprs, na.rm=TRUE)
                                           if('weights' %in% Biobase::assayDataElementNames(object)) y$weights <- Biobase::assayDataElement(object, 'weights')
                                           }
                                        else{if(is.vector(object)) y$exprs <- matrix(object, nrow=1)
                                             else y$exprs <- as.matrix(object)
                                             y$Amean <- rowMeans(y$exprs, na.rm=TRUE)
                                             }
                                        }
                                   }
                              }
                         }
                    if(mode(y$exprs)!='numeric') stop("data object doesn't contain numeric expression values")
                    if(is.null(rownames(y$exprs)) && !is.null(row.names(y$probes))) rownames(y$exprs) <- row.names(y$probes)
                    y
                    }

#' function to retrieve the contrasts stored in a dataset and fit an ANOVA on
#'   those
#'
#' this function is a wrapper and an adapter around \code{'\link[stats]{aov}'},
#'   it retrieves the contrasts for each independent variable specified in an
#'   aov() model and stored in the original dataset, and fits an ANOVA model on
#'   all contrasts for each independent variable, in order to obtained the
#'   computed SS and MS for further calculations, otherwise not available from the
#'   original ANOVA model as displayed by summary() or summary.lm(), the model
#'   fitted for all contrasts of one variable is the same as displayed by
#'   summary.lm() on the original ANOVA model
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the fitted ANOVA models for all the contrasts stored for each variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{fit_aov}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' fit_aov(goggles.model)
#' fit_aov(simple.model)
#' fit_aov(depression.model)
#'
#' @importFrom dplyr arrange mutate
#' @importFrom stats aov
#'
#' @export
nonEst <- function(x){
                   x <- as.matrix(x)
                   p <- ncol(x)
                   QR <- qr(x)
                   if(QR$rank<p){
                      n <- colnames(x)
                      if(is.null(n)) n <- as.character(1:p)
                      notest <- n[QR$pivot[(QR$rank+1):p]]
                      blank <- notest==''
                      if(any(blank)) notest[blank] <- as.character(((QR$rank+1):p)[blank])
                      return(notest)
                      }
                   else{return(NULL)}
                   }
