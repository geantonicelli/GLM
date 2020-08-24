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
#' @param anova
#'
#' @return the function returns an object of class 'pml' of the 'phangorn'
#'   package. Advanced and elaborated plots can be drawn in later steps based on
#'   the tree data of the pml class object
#'
#' @author gerardo esteban antonicelli
#'
#' @seealso \code{'\link{check_contrasts}'} \code{'\link{es}'}
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

#' function to check if a set of contrasts for a variable are orthogonal and
#' consistent
#'
#' this function is a wrapper around the functions
#' \code{'\link[seqinr]{dist.alignment}'}, \code{'\link[ape]{dist.dna}'},
#' \code{'\link[ape]{nj}'}, \code{'\link[ape]{bionj}'},
#' \code{'\link[ape]{fastme.bal}'}, \code{'\link[ape]{fastme.ols}'},
#' \code{'\link[phangorn]{pml}'} and \code{'\link[phangorn]{optim.pml}'}. it
#' takes a sequences alignment in format 'alignment' of 'DNAbin' matrix and
#' perform all transformations and steps to calculate a phylogenetic distance
#' matrix based on similarity or identity in the case of proteins or based in
#' evolutionary models in the case of DNA or RNA, to perform a likelihood-based
#' phylogenetic clustering and to optimise the phylogeny by a maximum likelihood
#' algorithm
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
#' mytree <- max_likelihood(fastaRNA, type=RNA, clustering=fastme.bal,
#'                         pml.model=GTR, clean=FALSE)
#' \dontrun{mytree <- max_likelihood(phylipProt, type=protein, pml.model=Blosum62,
#'                         outgroup=YP_0010399)}
#' \dontrun{plot.phylo(mytree, type='u')}
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
