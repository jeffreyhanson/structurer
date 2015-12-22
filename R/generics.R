#' @include structurer-internal.R misc.R
NULL
 
#' Number of loci
#'
#' This function returns the number of loci in a \code{Structure} object.
#'
#' @param x \code{StructureData}, \code{StructureAnalysis}.
#' @return \code{integer}.
#' @export
n.loci <- function(x) UseMethod('n.loci')

#' Negative Log Likelihood
#'
#' This function returns the negative loglikelihood of a Structure analysis.
#'
#' @param x \code{StructureReplicate}, \code{StructureResults}, \code{StructureAnalysis} or \code{StructureCollection} object.
#' @return \code{numeric}.
#' @export
loglik <- function(x) UseMethod('loglik')

#' Number of samples
#'
#' This function returns the number of samples in a \code{Structure} object.
#' 
#' @param x \code{StructureData}, \code{StructureAnalysis}.
#' @return \code{integer}.
#' @export
n.samples <- function(x) UseMethod('n.samples')

#' Sample membership
#'
#' This function returns the membership of samples in a \code{Structure} object.
#' 
#' @param x \code{StructureReplicate}, \code{StructureResults}, \code{StructureAnalysis}.
#' @param threshold \code{numeric} threshold to assign samples to populations. If \code{NULL} then probabilities for each population are returned. Defaults to \code{NULL}.
#' @return \code{integer}.
#' @export
sample.membership <- function(x, threshold) UseMethod('sample.membership')

#' Number of populations
#'
#' This function returns the number of populations sin a \code{Structure} object.
#' 
#' @param x \code{StructureReplicate}, \code{StructureResults}, or \code{StructureAnalysis}.
#' @return \code{integer}.
#' @export
n.pop <- function(x) UseMethod('n.pop')

#' Sample names
#'
#' This function returns the sample names in a \code{Structure} object.
#' 
#' @param x \code{StructureData}, \code{StructureAnalysis}.
#' @param value \code{character} vector of new names;
#' @return \code{character}.
#' @export
sample.names <- function(x) UseMethod('sample.names')

#' @rdname sample.names
#' @export
`sample.names<-` <- function(x, value) UseMethod('sample.names<-')

#' Loci names
#'
#' This function returns the loci names in a \code{Structure} object.
#' 
#' @param x \code{StructureData}, \code{StructureAnalysis}.
#' @param value \code{character} new loci names.
#' @return \code{character}.
#' @export
loci.names <- function(x) UseMethod('loci.names')

#' @rdname loci.names
#' @export
`loci.names<-` <- function(x, value) UseMethod('loci.names<-')

#' Subset samples
#'
#' This function returns a subset of the sampels in a \code{Structure} object.
#' 
#' @param x \code{StructureData}, \code{StructureAnalysis}.
#' @param samples \code{character}, \code{numeric}, or code{logical} indicating the samples to return.
#' @return \code{StructureData}.
#' @export
sample.subset <- function(x, samples) UseMethod('sample.subset')

#' Subset loci
#'
#' This function returns a subset of loci in a \code{Structure} object.
#' 
#' @param x \code{StructureData}, or \code{StructureAnalysis}.
#' @param loci \code{character}, \code{numeric}, or code{logical} indicating the loci to return.
#' @return \code{StructureData}.
#' @export
loci.subset <- function(x, loci) UseMethod('loci.subset')

#' Print objects
#'
#' This function prints objects.
#'
#' @param x the object to print.
#' @param header \code{logical} should object header be shown?
#' @param ... not used.
#' @name print
NULL

#' Show objects
#'
#' This function shows objects.
#'
#' @param object the object to show.
#' @name show
NULL

#' Plot negative log-likelihood values
#'
#' This function plots the negative log-likelihood values in a \code{StructureCollection} object.
#' 
#' @param x \code{StructureCollection} object.
#' @param main \code{character} plot title.
#' @seealso \code{\link{StructureCollection}}.
loglik.plot <- function(x, main) UseMethod('loglik.plot')


#' Plot delta-K values
#'
#' This function plots the delta-K values in a \code{StructureCollection} object.
#' 
#' @param x \code{StructureCollection} object.
#' @param main \code{character} plot title.
#' @seealso \code{\link{StructureCollection}}.
delta.k.plot <- function(x, main) UseMethod('delta.k.plot')
