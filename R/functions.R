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
#' @seealso \code{'\link{retrieve_seqs}'} \code{'\link{print_alignment}'}
#'   \code{'\link{clean_alignment}'} \code{'\link{load_alignment}'}
#'   \code{'\link{make_tree}'} \code{'\link{max_parsimony}'}
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
#' @importFrom seqinr dist.alignment as.matrix.alignment
#' @importFrom ape as.alignment as.DNAbin dist.dna nj bionj fastme.bal
#'   fastme.ols root
#' @importFrom phangorn as.phyDat pml optim.pml
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
