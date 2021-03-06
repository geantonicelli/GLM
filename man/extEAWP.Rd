% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limma.R
\name{extEAWP}
\alias{extEAWP}
\alias{\alias{extEAWP}}
\title{extract basic data from expression data objects}
\usage{
extEAWP(object)
}
\arguments{
\item{object}{any matrix-like object containing log-expression values, it can
be an object of class MAList, EList, marrayNorm, PLMset, vsn, or any class
inheriting from ExpressionSet, or any object that can be coerced to a
numeric matrix}
}
\value{
a list with components * exprs numeric matrix of log-ratios,
log-intensities or log-expression values * weights numeric matrix of
weights * probes data.frame of probe-annotation * Amean numeric vector of
average log-expression for each probe

exprs is the only required component, the other components will be NULL if
not found in the input object
}
\description{
given an expression data object of any known class, get the expression
values, weights, probe annotation and A-values that are needed for linear
modelling. This function is called by the linear modelling functions in the
limma package or by limmaDE() in the GLM package
}
\details{
rows correspond to probes and columns to RNA samples in the case of
two-color microarray data objects (MAList or marrayNorm), Amean is the
vector of row means of the matrix of A-values. For other data objects,
Amean is the vector of row means of the matrix of expression values. from
April 2013, the rownames of the output exprs matrix are required to be
unique. If object has no row names, then the output rownames of exprs are
1:nrow(object). If object has row names but with duplicated names, then the
rownames of exprs are set to 1:nrow(object) and the original row names are
preserved in the ID column of probes. object should be a normalized data
object. getEAWP will return an error if object is a non-normalized data
object such as RGList or EListRaw, because these do not contain
log-expression values
}
\seealso{
\code{'\link{check_contrasts}'} \code{'\link{omega_factorial}'}
}
\author{
gerardo esteban antonicelli
}
