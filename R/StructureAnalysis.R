#' @include structurer-internal.R misc.R generics.R StructureOpts.R StructureData.R StructureResults.R
NULL

#' StructureAnalysis: An S4 class to represent inputs and outputs from Structure
#'
#' This class stores input data and associated output results from the Structure program.
#'
#' @slot opts \code{StructureOpts} object with parameters used to run Structure .
#' @slot data \code{StructureData} object with input data used for analysis.
#' @slot results \code{StructureResults} object with results from analysis.
#' @seealso \code{\link{StructureAnalysis}}.
#' @export
setClass(
	"StructureAnalysis",
	representation(
		opts='StructureOpts',
		data='StructureData',
		results='StructureResults'
	),
	validity=function(object) {
		## opts
		# checks are internal
		## data
		# checks are internal
		## results
		# checks are internal
		## cross-object checks
		expect_equal(n.samples(object@data), n.samples(object@results))
		return(TRUE)
	}
)

#' Create StructureAnalysis object
#'
#' This function creates a new \code{StructureAnalysis} object.
#'
#' @param opts \code{StructureOpts} object with parameters used to run Structure .
#' @param data \code{StructureData} object with input data used for analysis.
#' @param results \code{StructureResults} object with results from analysis.
#' @seealso \code{\link{StructureAnalysis-class}}, \code{\link{StructureData}}, \code{\link{StructureData}}, \code{\link{StructureResults}}.
#' @export
StructureAnalysis<-function(opts, data, results) {
	x<-new("StructureAnalysis", opts=opts, data=data, results=results)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.loci
#' @method n.loci StructureAnalysis
#' @export
n.loci.StructureAnalysis <- function(x) {
	return(n.loci(x@data))
}

#' @rdname n.pop
#' @method n.pop StructureAnalysis
#' @export
n.pop.StructureAnalysis <- function(x) {
	return(n.pop(x@opts))
}

#' @rdname n.samples
#' @method n.samples StructureAnalysis
#' @export
n.samples.StructureAnalysis <- function(x) {
	return(n.samples(x@data))
}

#' @rdname sample.names
#' @method sample.names StructureAnalysis
#' @export
sample.names.StructureAnalysis <- function(x) {
	return(sample.names(x@data))
}

#' @rdname sample.names
#' @method sample.names<- StructureAnalysis
#' @export
`sample.names<-.StructureAnalysis` <- function(x, value) {
	return(
		StructureAnalysis(
			opts=x@opt,
			data=`sample.names<-`(x@data, value),
			results=`sample.names<-`(x@results, value)
		)
	)
}

#' @rdname loci.names
#' @method loci.names StructureAnalysis
#' @export
loci.names.StructureAnalysis <- function(x) {
	return(loci.names(x@data))
}

#' @rdname loci.names
#' @method loci.names<- StructureAnalysis
#' @export
`loci.names<-.StructureAnalysis` <- function(x, value) {
	return(
		StructureAnalysis(
			opts=x@opts,
			`loci.names<-`(x@data, value),
			results=x@results
		)
	)
}

#' @rdname sample.membership
#' @method sample.membership StructureAnalysis
#' @export
sample.membership.StructureAnalysis <- function(x, threshold=NULL) {
	return(sample.membership(x@results, threshold))
}


#' Run Structure anaylsis for a single MAXPOPS parameter.
#'
#' This function analyses data using Structure. It generates replicate runs for a single MAXPOPS parameter.
#'
#' @param x \code{StructureData} object.
#' @inheritParams StructureOpts
#' @param dir \code{character} with directory to use for analysis.
#' @param clean \code{logical} should input and output files be deleted after analysis is finished?
#' @seealso \code{StructureData}, \code{StructureOpts}.
#' @examples
#' # run Structure using low number of iterations
#' dat <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
#' x <- run.single.Structure(dat, NUMRUNS=1, MAXPOPS=2, BURNIN=10,
#'	NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10)
#' @export
run.single.Structure<-function(x, NUMRUNS=2, MAXPOPS=2, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, SEED=NA_real_, dir=tempdir(), clean=TRUE) {
	## initialization
	# argument checks
	opts <- StructureOpts(NUMRUNS=NUMRUNS, MAXPOPS=MAXPOPS, BURNIN=BURNIN, NUMREPS=NUMREPS, NOADMIX=NOADMIX, ADMBURNIN=ADMBURNIN, SEED=SEED)
	expect_is(x, 'StructureData')
	# identify structure path
	structure.path <- switch(
		Sys.info()['sysname'],
		'Linux'=system.file('bin', 'structure_linux', package='structurer'),
		'Darwin'=system.file('bin', 'structure_mac', package='structurer'),
		'Windows'=system.file('bin', 'structure_win.exe', package='structurer')
	)
	# update permissions
	if (!grepl(basename(structure.path), 'win'))
		system(paste0('chmod 700 ',structure.path))
	# sink opts and data to file
	write.StructureOpts(opts,dir)
	write.StructureData(x,file.path(dir, 'data.txt'))
	# run BayeScan
	return(
		StructureAnalysis(
			results=StructureResults(
				replicates=lapply(
					seq_len(opts@NUMRUNS),
					function(i) {
						system(paste0(structure.path, ' ', '-m ',file.path(dir, 'mainparams.txt'),' -e ',file.path(dir, 'extraparams.txt'),' -K ',2,' -L ',n.loci(x),' -N ',n.samples(x),' -i ',file.path(dir, 'data.txt'),' -o ',file.path(dir, 'output.txt')))
						return(read.StructureReplicate(file.path(dir, 'output.txt_f')))
					}
				)
			),
			data=x,
			opts=opts
		)
	)
}

#' @method print StructureAnalysis
#' @rdname print
#' @export
print.StructureAnalysis=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureAnalysis object.\n\n")
	cat('Options','\n')
	print(x@opts, header=FALSE)
	cat('Data','\n')
	print(x@data, header=FALSE)
	cat('Results','\n')
	print(x@results, header=FALSE)
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureAnalysis',
	function(object)
		print.StructureAnalysis(object)
)

