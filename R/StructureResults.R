#' @include structurer-internal.R misc.R generics.R StructureReplicate.R
NULL

#' StructureResults: An S4 class to results from Structure
#'
#' This class stores results from the Structure program.
#'
#' @slot summary \code{matrix} object containing average membership probabilities across replicates.
#' @slot replicates \code{list} of \code{StructureReplicate} objects.
#' @seealso \code{\link{StructureResults}}, \code{\link{StructureReplicate}}.
#' @export
setClass(
	"StructureResults",
	representation(
		summary='matrix',
		replicates='list'
	),
	validity=function(object) {
		# check that all replicates have the same properties
		expect_equal(length(unique(sapply(object@replicates, n.samples))), 1)
		expect_equal(length(unique(sapply(object@replicates, n.pop))), 1)
		return(TRUE)
	}
)

#' Create StructureResults object
#'
#' This function creates a new \code{StructureResults} object.
#'
#' @param summary \code{data.frame} object containing overall results from Structure replicates.
#' @param replicates \code{list} of \code{StructureReplicate} objects.
#' @seealso \code{\link{StructureReplicate-class}}.
#' @return \code{\link{StructureReplicate}}.
#' @export
StructureResults<-function(summary=NULL, replicates) {
	# compute summary results
	if (is.null(summary)) {
		if (length(replicates)>1) {
			dat <- lapply(replicates, slot, name='matrix')
			summary <- sapply(seq_len(ncol(dat[[1]])), function(i) {
				rowMeans(sapply(dat, function(x) {x[,i]}))
			})
		} else {
			summary <- replicates[[1]]@matrix
		}
	}
	# return new object
	x<-new("StructureResults", summary=summary, replicates=replicates)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.pop
#' @method n.pop StructureResults
#' @export
n.pop.StructureResults <- function(x) {
	return(n.pop(x@replicates[[1]]))
}

#' @rdname n.samples
#' @method n.samples StructureResults
#' @export
n.samples.StructureResults <- function(x) {
	return(n.samples(x@replicates[[1]]))
}

#' @rdname sample.names
#' @method sample.names StructureResults
#' @export
sample.names.StructureResults <- function(x) {
	return(x@replicates[[1]]@sample.names)
}

#' @rdname sample.membership
#' @method sample.membership StructureResults
#' @export
sample.membership.StructureResults <- function(x, threshold=NULL) {
	if (is.null(threshold))
		return(x@summary)
	pops <- apply(x@summary, 1, which.max)
	pops[which(apply(x@summary, 1, function(i) {all(i < threshold)}))] <- NA
	return(pops)
}

#' @rdname loglik
#' @method loglik StructureResults
#' @export
loglik.StructureResults <- function(x) {
	return(mean(sapply(x@replicates, loglik.StructureReplicate)))
}

#' @method print StructureResults
#' @rdname print
#' @export
print.StructureResults=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureResults object.\n")
	cat('  K:',n.pop(x),'\n')
	cat('  loglik:',loglik.StructureResults(x),'\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureResults',
	function(object)
		print.StructureResults(object)
)

 
 

 
