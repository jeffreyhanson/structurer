#' @include structurer-internal.R misc.R generics.R
NULL
 
#' ClumppOpts: An S4 class to represent Structure parameters
#'
#' This class stores input parameters for the Structure program.
#'
#' @slot M \code{character} name of search method. Valid arguments are 'FullSearch', 'Greedy', or 'LargeKGreedy'. Defaults to 'Greedy'.
#' @slot W \code{logical} weight populations by number of individuals? Defaults to \code{TRUE}. 
#' @slot S \code{logical} if \code{TRUE} the \code{G} matrix similarity statistic ised used. Else the \code{G\\'} statistic is used. Defaults to code{FALSE}.
#' @slot REPEATS \code{numeric} Number of random input orders tested. Defaults to 1000.
#' @seealso \code{\link{ClumppOpts}}.
#' @export
setClass(
	"ClumppOpts",
	representation(
		M="character",
		W="logical",
		S="logical",
		REPEATS="numeric"
	),
	prototype=list(
		M='Greedy',
		W=TRUE,
		S=FALSE,
		REPEATS=1000
	),
	validity=function(object) {
		# M
		match.arg(object@M, c('FullSearch','Greedy','LargeKGreedy'))
		# check that parameters are not NA
		sapply(
			c('W','S', 'REPEATS'),
			function(x) {
				expect_true(!is.na(slot(object, x)))
				return(invisible())
		})
		# REPEATS
		expect_true(object@REPEATS > 0)
		return(TRUE)
	}
)

#' Create ClumppOpts object
#'
#' This function creates a new \code{ClumppOpts} object.
#'
#' @param M \code{character} name of search method. Valid arguments are 'FullSearch', 'Greedy', or 'LargeKGreedy'. Defaults to 'Greedy'.
#' @param W \code{logical} weight populations by number of individuals? Defaults to \code{TRUE}. 
#' @param S \code{logical} if \code{TRUE} the \code{G} matrix similarity statistic ised used. Else the \code{G'} statistic is used. Defaults to code{FALSE}.
#' @param REPEATS \code{numeric} Number of random input orders tested. Defaults to 1000.
#' @seealso \code{\link{ClumppOpts-class}}.
#' @examples 
#' ClumppOpts(M='Greedy', W=TRUE, S=FALSE, REPEATS=1000)
#' @export
ClumppOpts<-function(M='Greedy', W=TRUE, S=FALSE, REPEATS=1000) {
	x<-new("ClumppOpts", M=M, W=W, S=S, REPEATS=REPEATS)
	validObject(x, test=FALSE)
	return(x)
}

#' @method print ClumppOpts
#' @rdname print
#' @export
print.ClumppOpts=function(x, ..., header=TRUE) {
	if (header)
		cat("ClumppOpts object.\n")
	cat('  M:',x@M,'\n')
	cat('  W:',x@W,'\n')
	cat('  S:',x@S,'\n')
	cat('  REPEATS:',x@REPEATS,'\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'ClumppOpts',
	function(object)
		print.ClumppOpts(object)
)

#' Write Clumpp parameters to file
#'
#' This function writes a \code{ClumppOpts} object to file.
#' @param x \code{ClumppOpts} object to file.
#' @param file \code{character} file path to save file.
#' @seealso \code{ClumppOpts}.
write.ClumppOpts <- function(x, file) {
	cat(paste0("
# --------------- Main parameters ---------------------------------------------

# file path parameters
DATATYPE 1
POPFILE popfile.txt
OUTFILE outfile.txt
MISCFILE miscfile.txt

# data-dependent parameters - these are overwritten using commandline arguments
K 1
C 1

# opts dependent parameters-  these are hard coded into the file
R 1
M ",switch(x@M, 'FullSearch'='1', 'Greedy'='2', 'LargeKGreedy'='3'),"
W ",c('2','1')[as.integer(x@W)+1],"
S ",c('2','1')[as.integer(x@S)+1],"

# - Additional options for the Greedy and LargeKGreedy algorithm (M = 2 or 3) -

GREEDY_OPTION 2
REPEATS ",x@REPEATS,"

# --------------- Optional outputs --------------------------------------------

PRINT_PERMUTED_DATA 1
PERMUTED_DATAFILE perm_data.txt
PRINT_EVERY_PERM 0
PRINT_RANDOM_INPUTORDER 0

# --------------- Advanced options --------------------------------------------

OVERRIDE_WARNINGS 0
ORDER_BY_RUN 1

"), file=file)
}

