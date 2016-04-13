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
#' @param replicates \code{list} of \code{StructureReplicate} objects.
#' @seealso \code{\link{StructureReplicate-class}}.
#' @return \code{\link{StructureReplicate}}.
#' @export
StructureResults<-function(replicates, opts=ClumpOpts(), dir=tempdir()) {
	### clumpp analysis to combine replicates
	# sink replicates to file
	if (length(replicates)>1 & n.pop(replicates[[1]])>1) {
		# init
		clumpp.path <- switch(
			Sys.info()['sysname'],
			'Linux'=system.file('bin', 'CLUMPP_linux', package='structurer'),
			'Darwin'=system.file('bin', 'CLUMPP_mac', package='structurer'),
			'Windows'=system.file('bin', 'CLUMPP_win.exe', package='structurer')
		)
		# update permissions
		if (!grepl(basename(clumpp.path), 'win'))
			system(paste0('chmod 700 ',clumpp.path))
		# sink replicates to file
		write.ClumppOpts(opts, file.path(dir, 'paramfile'))
		write.ClumppReplicates(replicates, file.path(dir,'popfile.txt'))
		# run clumpp
		o<-system(paste0(clumpp.path,' ',file.path(dir, 'paramfile'),' -p ',file.path(dir, 'popfile.txt'),' -o ',file.path(dir, 'outfile.txt'),' -k ',n.pop(replicates[[1]]),' -c ',n.samples(replicates[[1]]),' -r ',length(replicates),' -j ',file.path(dir, 'miscfile.txt'), ' > ',file.path(dir, 'clumpp.log'),' 2>&1'), intern=TRUE)
		# delete extra files
		if (file.exists('perm_data.txt')) unlink('perm_data.txt')
		if (file.exists(file.path(dir,'perm_data.txt'))) unlink(file.path(dir,'perm_data.txt'))
		# load in replicates
		summary.dat <- try(read.ClumppReplicates(file.path(dir,'outfile.txt')))
		if (inherits(summary.dat, 'try-error'))
		stop(paste(readLines(file.path(dir, 'clumpp.log')), collapse='\n'))
	} else {
		summary.dat <- replicates[[1]]@matrix
	}
	### return new object
	x<-new("StructureResults", summary=summary.dat, replicates=replicates)
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

 
 
