#' function to check if a set of contrasts for a variable are orthogonal and consistent
#'
#' this function is a wrapper around the functions
#' \code{'\link[seqinr]{dist.alignment}'}, \code{'\link[ape]{dist.dna}'},
#' \code{'\link[ape]{nj}'}, \code{'\link[ape]{bionj}'},
#' \code{'\link[ape]{fastme.bal}'}, \code{'\link[ape]{fastme.ols}'},
#' \code{'\link[phangorn]{pml}'} and \code{'\link[phangorn]{optim.pml}'}. it takes a sequences alignment in
#' format 'alignment' of 'DNAbin' matrix and perform all transformations and steps
#' to calculate a phylogenetic distance matrix based on similarity or identity
#' in the case of proteins or based in evolutionary models in the case of DNA or
#' RNA, to perform a likelihood-based phylogenetic clustering and to optimise the phylogeny by a
#' maximum likelihood algorithm
#'
#' @param variable
#'
#' @return the function returns an object of class 'pml' of the 'phangorn'
#'   package. Advanced and elaborated plots can be drawn in later steps based on
#'   the tree data of the pml class object
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{omega_factorial}'} \code{'\link{es}'}
#'
#' @aliases \alias{check_contrasts}
#'
#' @examples
#' data(fastaRNA)
#' data(phylipProt)
#' mytree <- max_likelihood(fastaRNA, type=RNA, clustering=fastme.bal,
#'                         pml.model=GTR, clean=FALSE)
#' \dontrun{mytree <- max_likelihood(phylipProt, type=protein, pml.model=Blosum62,
#'                         outgroup=YP_0010399)}
#' \dontrun{plot.phylo(mytree, type='u')}
#'
#' @export
check_contrasts <- function(variable){
                            dim <- dim(attr(variable, 'contrasts'))
                            rowsfactors <- vector()
                            print('row factors to check consistency:', quote=FALSE)
                            for(i in 1:dim[1]){
                                ortho <- 1
                                consis <- 1
                                for(j in 1:dim[2]){
                                    weight <- attr(variable, 'contrasts')[i, j]
                                    if(weight==0){ortho <- weight*ortho}
                                    else{ortho <- weight*ortho
                                         consis <- weight*consis}
                                    }
                                print(paste(rownames(attr(variable, 'contrasts'))[i], consis, sep='    '), quote=FALSE)
                                rowsfactors[i] <- ortho
                                }
                            sum <- sum(rowsfactors)
                            if(sum==0){
                               print(' ', quote=FALSE)
                               print('the constrasts are orthogonal', quote=FALSE)}
                            else{print(' ', quote=FALSE)
                                 print('caution: the contrasts are NOT orthogonal', quote=FALSE)}
                            }

#' function to calculate the effect size of a comparison
#'
#'
#' @param anova an object of class anova data.frame with the results of an ANOVA
#'   type III returned by Anova(model, type='III'), the model is created by
#'   aov()
#'
#' @return the function prints the omega squared for two factors and one
#'   interaction
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{es}'}
#'
#' @aliases \alias{omega_factorial}
#'
#' @examples
#' data(anova)
#' omega_factorial(anova)
#'
#' @export
omega_factorial <- function(anova){
                            a <- anova$Df[[2]]
                            b <- anova$Df[[3]]
                            ab <- anova$Df[[4]]
                            r <- anova$Df[[5]]
                            n <- sum(anova$Df)/((a+1)*(b+1))
                            SSa <- anova[[1]][2]
                            SSb <- anova[[1]][3]
                            SSab <- anova[[1]][4]
                            SSr <- anova[[1]][5]
                            # calculate mean squares by dividing the sums of squares by the degrees of freedom
                            MSa <- SSa/a
                            MSb <- SSb/b
                            MSab <- SSab/ab
                            MSr <- SSr/r
                            # calculate the variance estimates
                            varA <- (a*(MSa-MSr))/(n*(a+1)*(b+1))
                            varB <- (b*(MSb-MSr))/(n*(a+1)*(b+1))
                            varAB <- (a*b*(MSab-MSr))/(n*(a+1)*(b+1))
                            varTotal <- varA+varB+varAB+MSr
                            # calculate omega^2 by dividing each variance estimate by the total variance estimate
                            print(paste('omega-squared A : ', varA/varTotal), quote=FALSE)
                            print(paste('omega-squared B : ', varB/varTotal), quote=FALSE)
                            print(paste('omega-squared AB: ', varAB/varTotal), quote=FALSE)
                            }

#' function to calculate the effect size of a comparison
#'
#' this function is a wrapper and an adapter around the functions
#'   \code{'\link[pastecs]{stat.desc}'} and \code{'\link[compute.es]{mes}'} it
#'   takes a data.frame with variables and calculate the effect size of a
#'   comparison between two chunks of variation, they could be two levels of one
#'   factor, the combined effect of more than one level of a factor vs another
#'   level or between two levels of a factor constrained to one level of another
#'   factor, what is called simple effects analysis, all this possible comparisons
#'   between two chunks can be analysed with the same function by using orthogonal
#'   contrasts coded inside the data.frame as columns of dummy variables with the
#'   weights representing the comparison worth to be analysed closer
#'
#' @param data a character string without '' specifying a data.frame object with
#'   the data, each variable must be in only one column, one dummy variable with
#'   weights (in one column) for each contrast to be analysed must be provided
#' @param dep a character string without '' specifying the name of the dependent
#'   variable, this must be at the same time a column name in the data object
#' @param contrast a character string without '' specifying the name of the
#'   contrast to be analysed, this must be at the same time the name of a column
#'   for a dummy variable with weights specifying which samples should be
#'   compared
#' @param dig numeric an integer specifying the number of decimal digits to be
#'   printed out and also invisibly returned by the mes() function
#'
#' @return the function returns invisibly a data.frame with all coefficients
#'   from the calculation of effect size, in addition it prints out a summary
#'   with the most important coefficients
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{es}
#'
#' @examples
#' data(gogglesDataES)
#' data(depressionDataES)
#' es(gogglesDataES, attractiveness, alcEffect1)
#' es(gogglesDataES, attractiveness, alcEffect2)
#' es(gogglesDataES, attractiveness, gender_none)
#' es(gogglesDataES, attractiveness, gender_twoPints)
#' es(gogglesDataES, attractiveness, gender_fourPints)
#' es(depressionDataES, diff, all_vs_NoTreatment, dig=4)
#' es(depressionDataES, diff, treatment_vs_placebo)
#' es(depressionDataES, diff, old_vs_newDrug, dig=2)
#' es(depressionDataES, diff, old_vs_oldDrug)
#'
#' @importFrom dplyr filter count
#' @importFrom pastecs stat.desc
#' @importFrom compute.es mes
#'
#' @export
es <- function(data, dep, contrast, dig=3){
               filter <- substitute(contrast)
               compare <- data %>% filter(!eval(filter)==0)
               count <- compare %>% count(eval(filter))
               stats <- by(eval(substitute(compare$dep)),
                           eval(substitute(compare$contrast)),
                           stat.desc, basic=FALSE)
               m.1 <- stats[[1]][[2]]
               m.2 <- stats[[2]][[2]]
               sd.1 <- stats[[1]][[6]]
               sd.2 <- stats[[2]][[6]]
               n.1 <- count$n[[1]]
               n.2 <- count$n[[2]]
               effect_size <- mes(m.1, m.2, sd.1, sd.2, n.1, n.2, dig=dig, verbose=FALSE)
               out1 <- c(d=effect_size$d, var.d=effect_size$var.d, p.value=effect_size$pval.d,
                         U3=effect_size$U3.d, CLES=effect_size$cl.d, Cliff=effect_size$cliffs.d)
               out2 <- c(r=effect_size$r, var.r=effect_size$var.r, p.value=effect_size$pval.r,
                         fisher.z=effect_size$fisher.z, var.z=effect_size$var.z)
               out3 <- c(OR=effect_size$OR, p.value=effect_size$pval.or, logOR=effect_size$lOR)
               out4 <- round(c(NNT=effect_size$NNT, n.chunk1=n.1, n.chunk2=n.2), digits=1)
               out <- list('MeanDifference'=out1, 'Correlation'=out2, 'OddsRatio'=out3, 'others'=out4)
               cat('effect size for contrast', filter, '\n', '\n')
               print(out)
               invisible(effect_size)
               }

#' function to calculate correlation coefficients, confidence intervals and
#'   p-values for one or more pairs of variables
#'
#' this function is an extension of cor.test() that calculates correlation
#'   coefficients, p-values and confidence intervals for a pair of variables or
#'   for all pairs of variables if more than two are provided in a data.frame. It
#'   can calculate Pearson, Kendall or Spearman correlation coefficients and
#'   p-values for a two-tailed or one-tailed alternative hypothesis. For Kendall
#'   or Spearman correlation, confidence intervals are calculated by bootstrap
#'   samling. In case of NAs the complete row can be taken out or only for the
#'   affected comparisons
#'
#' @param data a character string without '' specifying a data.frame object with
#'   the data, each variable must be in only one column, or a subset of a
#'   data.frame with two or more variables or two vectors in a list
#' @param method a character string specifying the correlation method to be
#'   used, options are pearson, spearman or kendall, names can be abbreviated
#' @param alternative a character string indicating the alternative hypothesis
#'   and must be one of 'two.sided', 'greater' or 'less', 'greater' corresponds
#'   to a positive association, 'less' to a negative association, names can be
#'   abbreviated
#' @param conf.level numeric an decimal specifying the confidence level for the
#'   returned confidence interval
#' @param use an character string giving a method for computing covariances in
#'   the presence of missing values, options are (an abbreviation of) one of the
#'   strings 'complete.obs' or 'pairwise.complete.obs' (default)
#' @param boots a single positive integer indicating the number of bootstrap
#'   replicates, default=2000
#' @param verbose a logical indicating whether a summary of the output should be
#'   printed to the console
#' @param exact a logical indicating whether an exact p-value should be computed
#'   (used for Kendall's tau and Spearman's rho), see 'details' for the meaning
#'   of NULL (the default)
#' @param continuity logical: if true, a continuity correction is used for
#'   Kendall's tau and Spearman's rho when not computed exactly
#'
#' @return if the input are two vectors or two columns of a data.frame and the
#'   method is pearson the function returns invisibly and print to the console
#'   the output of cor.test(), if the method is kendall or spearman the function
#'   returns invisibly a list and prints to the console (if verbose=TRUE) the
#'   output of cor.test() and the confidence interval from boot.ci(), if the
#'   input is a data.frame with more than two variables the function returns
#'   invisibly a list with the call, relevant input parameters and the matrices
#'   with correlation coefficients, p-values and lower and upper limits of
#'   confidence intervals for all combinations of variables, if verbose=TRUE the
#'   function prints to the console the matrices with correlation coefficients,
#'   p-values and lower and upper limits of confidence intervals for all
#'   combinations of variables
#'
#' @details the three methods each estimate the association between paired
#'   samples and compute a test of the value being zero. They use different
#'   measures of association, all in the range [-1, 1] with 0 indicating no
#'   association, these are sometimes referred to as tests of no correlation,
#'   but that term is often confined to the default method. If method is
#'   'pearson', the test statistic is based on Pearson's product moment
#'   correlation coefficient cor(x, y) and follows a t distribution with
#'   length(x)-2 degrees of freedom if the samples follow independent normal
#'   distributions. If there are at least 4 complete pairs of observation, an
#'   asymptotic confidence interval is given based on Fisher's Z transform. If
#'   method is 'kendall' or 'spearman', Kendall's tau or Spearman's rho
#'   statistic is used to estimate a rank-based measure of association, these
#'   tests may be used if the data do not necessarily come from a bivariate
#'   normal distribution. For Kendall's test, by default (if exact is NULL), an
#'   exact p-value is computed if there are less than 50 paired samples
#'   containing finite values and there are no ties, otherwise, the test
#'   statistic is the estimate scaled to zero mean and unit variance, and is
#'   approximately normally distributed. For Spearman's test, p-values are
#'   computed using algorithm AS 89 for n<1290 and exact=TRUE, otherwise via the
#'   asymptotic t approximation. Note that these are ‘exact’ for n<10, and use
#'   an Edgeworth series approximation for larger sample sizes (the cutoff has
#'   been changed from the original paper)
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{corr.test}
#'
#' @examples
#' data(examData)
#' data(depressionDataES)
#' corr.test(examData[ , c(2, 4)], 'pearson', 'greater', 0.99)
#' corr.test(list(examData$Revise, examData$Anxiety), 'kendall', 'two.sided', 0.99, boots=1000)
#' corr.test(examData[ , 2:4], 'spearman', 'less', 0.95, 'complete.obs', exact=FALSE, verbose=FALSE)
#'
#' @importFrom boot boot boot.ci
#'
#' @export
corr.test <- function(data, method='pearson', alternative='two.sided',
                      conf.level=0.95, use='pairwise.complete.obs',
                      boots=2000, verbose=TRUE, exact=NULL, continuity=FALSE){
                      method <- match.arg(method, choices=c('pearson', 'spearman', 'kendall'))
                      alternative <- match.arg(alternative, choices=c('two.sided', 'greater', 'less'))
                      use <- match.arg(use, choices=c('complete.obs', 'pairwise.complete.obs'))
                      Call <- match.call()
                      if(is.null(dim(data)) || dim(data)[[2]]==2){
                         data <- data.frame(data[[1]], data[[2]])
                         corr <- cor.test(data[[1]], data[[2]],
                                          alternative=alternative,
                                          method=method,
                                          exact=exact,
                                          conf.level=conf.level,
                                          continuity=continuity)
                         if(method!='pearson'){
                            bootR <- function(data, i) cor(data[[1]][i], data[[2]][i],
                                                           use=use, method=method)
                            boot.out <- boot(data, bootR, R=boots)
                            ci <- boot.ci(boot.out, conf=conf.level, type='bca')
                            out <- list(correlation=corr, confidence.interval=ci)
                            if(verbose==TRUE){
                               print(corr)
                               print(ci)
                               }
                            invisible(out)
                            }
                         else{return(corr)}
                         }
                      else{if(use=='complete.obs'){data <- na.omit(data)
                           }
                           col <- dim(data)[[2]]
                           s.est <- vector('list', col)
                           p.val <- vector('list', col)
                           lower.ci <- vector('list', col)
                           upper.ci <- vector('list', col)
                           for(i in 1:col){
                               estimate <- c()
                               p.value <- c()
                               lower <- c()
                               upper <- c()
                               for(j in 1:col){
                                   corr <- cor.test(data[[i]], data[[j]],
                                                    alternative=alternative,
                                                    method=method,
                                                    exact=exact,
                                                    conf.level=conf.level,
                                                    continuity=continuity)
                                   estimate[j] <- corr$estimate
                                   p.value[j] <- corr$p.value
                                   if(method=='pearson'){
                                      lower[j] <- corr$conf.int[1]
                                      upper[j] <- corr$conf.int[2]
                                      }
                                   else{if(i!=j){
                                           bootR <- function(data, k) cor(data[[i]][k],
                                                                          data[[j]][k],
                                                                          use=use,
                                                                          method=method)
                                           boot.out <- boot(data, bootR, R=boots)
                                           ci <- boot.ci(boot.out,
                                                         conf=conf.level,
                                                         type='bca')
                                           lower[j] <- ci$bca[4]
                                           upper[j] <- ci$bca[5]
                                           }
                                        else{lower[j] <- 1
                                             upper[j] <- 1
                                             }
                                        }
                                   }
                               s.est[[i]] <- estimate
                               p.val[[i]] <- p.value
                               lower.ci[[i]] <- lower
                               upper.ci[[i]] <- upper
                               }
                           names(s.est) <- colnames(data)
                           names(p.val) <- colnames(data)
                           names(lower.ci) <- colnames(data)
                           names(upper.ci) <- colnames(data)
                           s.est <- as.data.frame(s.est)
                           p.val <- as.data.frame(p.val)
                           lower.ci <- as.data.frame(lower.ci)
                           upper.ci <- as.data.frame(upper.ci)
                           rownames(s.est) <- colnames(data)
                           rownames(p.val) <- colnames(data)
                           rownames(lower.ci) <- colnames(data)
                           rownames(upper.ci) <- colnames(data)
                           s.est <- as.matrix(s.est)
                           p.val <- as.matrix(p.val)
                           lower.ci <- as.matrix(lower.ci)
                           upper.ci <- as.matrix(upper.ci)
                           if(verbose==TRUE){
                              cat('\n', corr$method, '\n', '\n', 'sample estimates', '\n')
                              print(s.est)
                              cat('\n', 'p-values', '\n')
                              print(p.val)
                              cat('\n', 'lower limit of', conf.level*100,
                                  '% confidence interval', '\n')
                              print(lower.ci)
                              cat('\n', 'upper limit of', conf.level*100,
                                  '% confidence interval', '\n')
                              print(upper.ci)
                              if(alternative=='two.sided'){
                                 cat('\n',
                                     'alternative hypothesis: true correlation is not equal to 0')}
                              else{cat('\n', 'alternative hypothesis is: true correlation is',
                                       corr$alternative, 'than 0')}
                              }
                           out <- list(call=Call, method=corr$method, alternative=corr$alternative,
                                       use=use, conf.level=conf.level, estimates=s.est,
                                       p.values=p.val, lowerCI=lower.ci, upperCI=upper.ci)
                           invisible(out)
                           }
                      }

#' function to calculate the omega squared and f squared effect size of independent
#'   variables in an ANOVA model
#'
#' this function uses the function fit_aov() to retrieve the contrasts for each
#'   independent variable specified in an aov() model and stored in the original
#'   dataset, and to fit an ANOVA model on all contrasts for each independent
#'   variable, in order to obtained the computed SS and MS for further calculations,
#'   otherwise not available from the original ANOVA model as displayed by summary()
#'   or summary.lm(), the model fitted for all contrasts of one variable is the same
#'   as displayed by summary.lm() on the original ANOVA model
#'   this function calculates the omega squared and f squared coefficients from the SS,
#'   MS and df calculated by aov() for each independent variable and for each contrast
#'   specified for each independent variable
#'
#' @param anova an object of class anova data.frame with the ANOVA model created
#'   by aov(), in the original dataset each variable must be in only one
#'   column, and all contrasts for each variable should be loaded as contrasts
#'   for that variable
#'
#' @return the function returns the omega squared and f squared effect sizes for
#'   each independent variable specified in the main ANOVA model and for each
#'   contrast in the ANOVA model fitted for each independent variable
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
#'
#' @aliases \alias{omega}
#'
#' @examples
#' data(gogglesData)
#' data(depressionData)
#' goggles.model <- aov(attractiveness ~ gender + alcohol + gender:alcohol, data=gogglesData)
#' simple.model <- aov(attractiveness ~ simple, data=gogglesData)
#' depression.model <- aov(diff ~ treat, data=depressionData)
#' omega(goggles.model)
#' omega(simple.model)
#' omega(depression.model)
#'
#' @export
omega <- function(anova){
                  contr <- fit_aov(anova)
                  ncontr <- length(contr)
                  colect <- vector('list', length=ncontr+1)
                  var.nam <- c('main_ANOVA', names(contr))
                  out <- vector('list', length=ncontr+1)
                  colect[[1]] <- anova
                  for(i in 1:ncontr){
                      colect[[i+1]] <- contr[[i]]
                      names(colect[[i+1]]) <- names(contr[[i]])
                      }
                  nloop <- length(colect)
                  for(j in 1:nloop){
                      loop <- colect[[j]]
                      if(length(loop)==0) next
                      model <- summary(loop)
                      SSm <- model[[1]][[2]]
                      MSr <- model[[1]][[3]][[length(SSm)]]
                      SSt <- sum(SSm)
                      df <- model[[1]][[1]]
                      omega.square <- c()
                      f.square <- c()
                      for(k in 1:(length(SSm)-1)){
                          i.SSm <- SSm[k]
                          i.df <- df[k]
                          num <- i.SSm-(i.df*MSr)
                          den <- SSt+MSr
                          omega <- num/den
                          omega.square[k] <- omega
                          f <- omega/(1-omega)
                          f.square[k] <- f
                          }
                      names <- rownames(model[[1]])
                      names <- names[1:length(SSm)-1]
                      names(omega.square) <- names
                      res <- as.matrix(as.data.frame(omega.square))
                      res <- cbind(res, f.square)
                      out[[j]] <- res
                      }
                  names(out) <- var.nam
                  out
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
#' @importFrom dplyr arrange
#' @importFrom stats aov
#'
#' @export
fit_aov <- function(anova){
                    aovcall <- anova$call
                    formula <- formula(aovcall)
                    y <- deparse(formula[[2]])
                    model <- anova[length(anova)]
                    nvar <- dim(model$model)[2]
                    Ntotal <- dim(model$model)[1]
                    fits <- vector('list', nvar)
                    names(fits) <- colnames(model$model)
                    for(i in 1:nvar){
                        var <- attr(model$model[[i]], 'contrasts')
                        if(is.null(var)) next
                        ncol <- dim(var)[2]
                        contrasts <- vector('list', ncol)
                        names(contrasts) <- colnames(var)
                        n <- Ntotal/length(levels(model$model[[i]]))
                        for(j in 1:ncol){
                            contrast <- var[ , j]
                            weights <- rep(contrast, each=n)
                            contrasts[[j]] <- weights
                            }
                        contrasts <- as.data.frame(contrasts)
                        dataset <- arrange(model$model, model$model[[i]])
                        dep <- dataset[y]
                        data <- cbind(dep, contrasts)
                        indep <- colnames(contrasts)
                        ind.for <- indep[1]
                        if(length(indep)>=2){
                           for(k in 2:length(indep)){
                               ind.for <- paste0(ind.for, ' + ', indep[k])
                               }
                           }
                        aovcall[[1L]] <- quote(stats::aov)
                        aovcall$formula <- as.formula(paste0(y, ' ~ ', ind.for))
                        aovcall$data <- quote(data)
                        fit <- eval(substitute(aovcall))
                        fits[[i]] <- fit
                        out <- fits[-1]
                        }
                    out
                    }
