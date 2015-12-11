#' @include structurer-internal.R misc.R generics.R StructureCollection.R
NULL

#' StructureCollection: An S4 class to represent a collection of inputs and outputs from Structure
#'
#' This class stores a collection of input data and associated output results from the Structure program.
#'
#' @slot summary \code{data.frame}  summary of analyses.
#' @slot analyses \code{list} of \code{StructureCollection} objects run using different .
#' @slot best \code{numeric} index of \code{StructureCollection} with best negative log-likelihood.
#' @seealso \code{\link{StructureCollection}}.
#' @export
setClass(
	"StructureCollection",
	representation(
		summary='data.frame',
		analyses='list',
		best='numeric'
	),
	validity=function(object) {
		# analysis
		o1 <<- object
		sapply(object@analyses, function(x) {
			expect_is(x, 'StructureAnalysis')
		})
		expect_equal(length(unique(sapply(object@analyses, n.samples))),1)
		# delta.k
		expect_equal(names(object@summary), c('k','mean.loglik','sd.loglik','delta.k','n.replicates'))
		expect_equal(nrow(object@summary),length(object@analyses))
		# best
		expect_true(is.finite(object@best))
		expect_true(object@best >= 1)
		expect_true(object@best <= length(object@analyses))
		return(TRUE)
	}
)

#' Create StructureCollection object
#'
#' This function creates a new \code{StructureCollection} object.
#'
#' @param summary \code{data.frame}  summary of analyses.
#' @param analyses \code{list} of \code{StructureCollection} objects run using different .
#' @param best \code{numeric} index of \code{StructureCollection} with best negative log-likelihood.
#' @seealso \code{\link{StructureCollection-class}}, \code{\link{StructureData}}, \code{\link{StructureData}}, \code{\link{StructureResults}}.
#' @details The delta-k values are calculated using the method used by Structure Harvestor \url{http://taylor0.biology.ucla.edu/structureHarvester/faq.html.}
#' @export
StructureCollection<-function(analyses, summary=NULL, best=NULL) {
	if (is.null(summary)) {
		## k is calued using method used in structureHarvestor, see 
		# extract ll
		ll <- sapply(analyses, function(x) {sapply(x@results@replicates, loglik.StructureReplicate)})
		ll.mean <- colMeans(ll)
		ll.sd <- apply(ll, 2, sd)
		# delta k calculations
		prevK <- ll.mean[1:(length(ll.mean)-2)]
		nextK <- ll.mean[3:(length(ll.mean))]
		thisK <- ll.mean[2:(length(ll.mean)-1)]
		thisKsd <- ll.sd[2:(length(ll.mean)-1)]
		delta.k <- abs(nextK - (2*thisK ) + prevK) / thisKsd
		# generate data.frame
		summary <- data.frame(
			k=sapply(analyses, function(x) {x@opts@MAXPOPS}),
			mean.loglik=ll.mean,
			sd.loglik=ll.sd,
			delta.k=c(NA, delta.k,NA),
			n.replicates=sapply(analyses, function(x) {length(x@results@replicates)})
		)
	}
	if (is.null(best)) {
		best <- which.max(summary$delta.k)
	}
	x<-new("StructureCollection", summary=summary, best=best, analyses=analyses)
	validObject(x, test=FALSE)
	return(x)
}

#' Run Structure anaylsis
#'
#' This function analyses data using Structure.
#'
#' @param x \code{StructureData} object.
#' @inheritParams StructureOpts
#' @param dir \code{character} with directory to use for analysis.
#' @param clean \code{logical} should input and output files be deleted after analysis is finished?
#' @seealso \code{StructureData}, \code{StructureOpts}.
#' @examples
#' # run Structure using low number of iterations
#' dat <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
#' x <- run.Structure(dat, NUMRUNS=2, MAXPOPS=1:3, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10)
#' print(x)
#' @export
run.Structure<-function(x, NUMRUNS=2, MAXPOPS=1:10, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, SEED=NA_real_, dir=tempdir(), clean=TRUE) {
	# run analysis
	return(
		StructureCollection(
			analyses=lapply(MAXPOPS, function(n) {
				run.single.Structure(x, NUMRUNS=NUMRUNS, MAXPOPS=n, BURNIN=BURNIN, NUMREPS=NUMREPS, NOADMIX=NOADMIX, ADMBURNIN=ADMBURNIN, SEED=SEED)
			}),
			best=NULL,
			summary=NULL
		)
	)
}

#' @rdname n.loci
#' @method n.loci StructureCollection
#' @export
n.loci.StructureCollection <- function(x) {
	return(n.loci(x@analyses[[1]]@data))
}

#' @rdname n.pop
#' @method n.pop StructureCollection
#' @export
n.pop.StructureCollection <- function(x) {
	return(n.pop(x@analyses[[x@best]]@opts))
}

#' @rdname n.samples
#' @method n.samples StructureCollection
#' @export
n.samples.StructureCollection <- function(x) {
	return(n.samples(x@analyses[[1]]@data))
}

#' @rdname sample.names
#' @method sample.names StructureCollection
#' @export
sample.names.StructureCollection <- function(x) {
	return(sample.names(x@analyses[[1]]@data))
}

#' @rdname sample.membership
#' @method sample.membership StructureCollection
#' @export
sample.membership.StructureCollection <- function(x) {
	return(x@analyses[[x@best]]@results@summary)
}

#' @rdname sample.names
#' @method sample.names<- StructureCollection
#' @export
`sample.names<-.StructureCollection` <- function(x, value) {
	return(
		StructureCollection(
			summary=x@summary,
			best=x@best,
			analyses=lapply(x@anaylses, `sample.names<-.StructureAnalysis`, value=value)
		)
	)
}

#' @rdname loci.names
#' @method loci.names StructureCollection
#' @export
loci.names.StructureCollection <- function(x) {
	return(loci.names(x@analyses[[1]]@data))
}

#' @rdname loci.names
#' @method loci.names<- StructureCollection
#' @export
`loci.names<-.StructureCollection` <- function(x, value) {
	return(
		StructureCollection(
			summary=x@summary,
			best=x@best,
			analyses=lapply(x@anaylses, `loci.names<-.StructureAnalysis`, value=value)
		)
	)
}

#' @method print StructureCollection
#' @rdname print
#' @export
print.StructureCollection=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureCollection object.\n\n")
	cat('Best K:',x@analyses[[x@best]]@opts@MAXPOPS,'\n')
	cat('Options','\n')
	print(x@analyses[[x@best]]@opts, header=FALSE)
	cat('Data','\n')
	print(x@analyses[[x@best]]@data, header=FALSE)
	cat('Results','\n')
	print(x@analyses[[x@best]]@results, header=FALSE)
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureCollection',
	function(object)
		print.StructureCollection(object)
)

#' @rdname loglik.plot
#' @method loglik.plot StructureCollection
#' @export
loglik.plot.StructureCollection <- function(x, main='') {
	## make plot
	dat <- x@summary
	dat$lower <- dat$mean.loglik - dat$sd.loglik
	dat$upper <- dat$mean.loglik + dat$sd.loglik
	ggplot(data=x@summary) +
		geom_point(aes_string(y='mean.loglik', x='k')) +
		geom_errorbar(aes_string(ymin='lower', ymax='upper')) +
		theme_classic() +
		xlab('Number of populations (k)') +
		ylab('Negative log-likelihood') +
		ggtitle(main)
}

#' @rdname delta.plot
#' @method delta.plot StructureCollection
#' @export
delta.plot.StructureCollection <- function(x, main='Delta-K plot') {
	## make plot
	ggplot(data=x@summary) +
		geom_point(aes_string(y='delta.k', x='k')) +
		geom_line(aes_string(y='delta.k', x='k')) +
		theme_classic() +
		xlab('Number of populations') +
		ylab('Relative support (delta-K)') +
		ggtitle(main)
}

