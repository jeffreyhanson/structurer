#' @include structurer-internal.R misc.R generics.R
NULL
 
#' StructureReplicate: An S4 class to results from Structure
#'
#' This class stores results from the Structure program.
#'
#' @slot loglik \code{numeric} mean negative loglikelihood of model.
#' @slot var_loglik \code{numeric} variance in negative loglikelihood of model.
#' @slot alpha \code{numeric} mean value of alpha.
#' @slot matrix \code{matrix} population membership probabilities. Each row is an individual. Each column is for a different population.
#' @slot sample.names \code{character} name of samples.
#' @slot log \code{character} log file.
#' @seealso \code{\link{StructureReplicate}}.
#' @export
setClass(
	"StructureReplicate",
	representation(
		loglik='numeric',
		var_loglik='numeric',
		alpha='numeric',
		matrix='matrix',
		sample.names='character',
		log='character'
	),
	validity=function(object) {
		# check slots are finite
		sapply(c('loglik', 'var_loglik'),
			function(x) {
				expect_true(is.finite(slot(object, x)))
			}
		)
		# check dimensions of matrix match log file
		pars <- gsub(' ', '', gsub(',', '', strsplit(grep('NUMINDS', object@log, value=TRUE),'\t')[[1]], fixed=TRUE))
		n.inds <- as.numeric(gsub('NUMINDS=', '', grep('NUMINDS', pars, fixed=TRUE, value=TRUE), fixed=TRUE))
		n.pop <- as.numeric(gsub('MAXPOPS=', '', grep('MAXPOPS', pars, fixed=TRUE, value=TRUE), fixed=TRUE))
		expect_equal(ncol(object@matrix), n.pop)
		expect_equal(nrow(object@matrix), n.inds)
		expect_equal(nrow(object@matrix), length(object@sample.names))
		# check alpha
		if (n.pop > 1) {
			expect_true(is.finite(slot(object, 'alpha')))
		} else {
			expect_true(is.na(slot(object, 'alpha')))
		}
		return(TRUE)
	}
)

#' Create StructureReplicate object
#'
#' This function creates a new \code{StructureReplicate} object.
#'
#' @param loglik \code{numeric} mean negative loglikelihood of model.
#' @param var_loglik \code{numeric} variance in negative loglikelihood of model.
#' @param alpha \code{numeric} mean value of alpha.
#' @param matrix \code{matrix} population membership probabilities. Each row is an individual. Each column is for a different population.
#' @param sample.names \code{character} name of samples.
#' @param log \code{character} log file.
#' @seealso \code{\link{StructureReplicate-class}}.
#' @return \code{\link{StructureReplicate}}.
#' @export
StructureReplicate<-function(loglik, var_loglik, alpha, matrix, sample.names, log) {
	x<-new("StructureReplicate", loglik=loglik, var_loglik=var_loglik, alpha=alpha, matrix=matrix, sample.names=sample.names, log=log)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.pop
#' @method n.pop StructureReplicate
#' @export
n.pop.StructureReplicate <- function(x) {
	return(ncol(x@matrix))
}

#' @rdname n.samples
#' @method n.samples StructureReplicate
#' @export
n.samples.StructureReplicate <- function(x) {
	return(nrow(x@matrix))
}

#' @rdname sample.names
#' @method sample.names StructureReplicate
#' @export
sample.names.StructureReplicate <- function(x) {
	return(x@sample.names)
}

#' @rdname sample.membership
#' @method sample.membership StructureReplicate
sample.membership.StructureReplicate <- function(x) {
	return(apply(x@matrix, 1, which.max))
}

#' @rdname loglik
#' @method loglik StructureReplicate
loglik.StructureReplicate <- function(x) {
	return(x@loglik)
}

#' Read Structure run
#'
#' This function reads the results of a single run of the Structure program.
#'
#' @param file \code{character} file path of output file.
#' @seealso \code{\link{StructureReplicate-class}}.
#' @return \code{\link{StructureReplicate}}.
#' @export
read.StructureReplicate<-function(file) {
	# load file
	logfile <- suppressWarnings(readLines(file))
	# load matrix
	pars <- gsub(' ', '', gsub(',', '', strsplit(grep('NUMINDS', logfile, value=TRUE),'\t')[[1]], fixed=TRUE))
	n.inds <- as.numeric(gsub('NUMINDS=', '', grep('NUMINDS', pars, fixed=TRUE, value=TRUE), fixed=TRUE))
	start.line <- grep('Inferred ancestry of individuals', logfile, fixed=TRUE)+1
	mat <- read.table(file, skip=start.line, nrows=n.inds)
	# load alpha
	alpha <-as.numeric(gsub('Mean value of alpha         = ', '', grep('Mean value of alpha', logfile, fixed=TRUE, value=TRUE), fixed=TRUE)) 
	if (length(alpha)==0)
		alpha <- NA_real_
	# return object
	StructureReplicate(
		loglik=as.numeric(gsub('Mean value of ln likelihood = ', '', grep('Mean value of ln likelihood', logfile, fixed=TRUE, value=TRUE), fixed=TRUE)),
		var_loglik=as.numeric(gsub('Variance of ln likelihood   = ', '', grep('Variance of ln likelihood', logfile, fixed=TRUE, value=TRUE), fixed=TRUE)),
		alpha=alpha,
		matrix=as.matrix(mat[,c(-1, -2, -3, -4),drop=FALSE]),
		sample.names=as.character(mat[[2]]),
		log=logfile
	)
}

#' @method print StructureReplicate
#' @rdname print
#' @export
print.StructureReplicate=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureReplicate object.\n")
	cat('  K:',n.pop(x),'\n')
	cat('  loglik:',loglik(x),'\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureReplicate',
	function(object)
		print.StructureReplicate(object)
)

