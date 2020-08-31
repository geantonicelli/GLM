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
#' @param anova an object of class anova data.frame with the results of an ANOVA type III returned by Anova(model, type='III'), the model is created by aov()
#'
#' @return the function prints the omega squared for two factors and one interaction
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
#' @importFrom boot boot boot.ci
#'
#' @export
corr.test <- function(data, method='pearson', alternative='two.sided',
                      conf.level=0.95, use='pairwise.complete.obs',
                      boots=1000, verbose=TRUE){
                      if(is.null(dim(data)) || dim(data)[[2]]==2){
                         data <- data.frame(data[[1]], data[[2]])
                         corr <- cor.test(data[[1]], data[[2]],
                                          alternative=alternative,
                                          method=method,
                                          conf.level=conf.level)
                         if(method!='pearson'){
                            bootR <- function(data, i) cor(data[[1]][i], data[[2]][i], use=use, method=method)
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
                                                    conf.level=conf.level)
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
                           out <- list(method=corr$method, alternative=corr$alternative,
                                       use=use, conf.level=conf.level,
                                       estimates=s.est, p.values=p.val,
                                       lowerCI=lower.ci, upperCI=upper.ci)
                           invisible(out)
                           }
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
#' @importFrom boot boot boot.ci
#'
#' @export
omega <- function(anova){
  ee <<- environment() # trick to see what's inside of the function's execution environment
  print(environment())
  # print(ls.str()) # alternatively
  # ------------------------------------------------------------------------------------------
  model <- summary(anova)
  SSm <- model[[1]][[2]]
  MSr <- model[[1]][[3]][[length(SSm)]]
  SSt <- sum(SSm)
  df <- model[[1]][[1]]
  omega.square <- c()
  f.square <- c()
  for(i in 1:length(SSm)-1){
    i.SSm <- SSm[i]
    i.df <- df[i]
    num <- i.SSm-(i.df*MSr)
    den <- SSt+MSr
    omega <- num/den
    omega.square[i] <- omega
    f <- omega/(1-omega)
    f.square[i] <- f
  }
  names <- rownames(model[[1]])
  names <- names[1:length(SSm)-1]
  names(omega.square) <- names
  out <- as.matrix(as.data.frame(omega.square))
  out <- cbind(out, f.square)
  return(out)
}
