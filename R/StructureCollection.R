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
		sapply(object@analyses, function(x) {
			expect_is(x, 'StructureAnalysis')
		})
		expect_equal(length(unique(sapply(object@analyses, n.samples))),1)
		# delta.k
		expect_equal(names(object@summary), c('k','mean.lnprob','sd.lnprob','delta.k','n.replicates'))
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
#' @param analyses \code{list} of \code{StructureCollection} objects run using different .
#' @seealso \code{\link{StructureCollection-class}}, \code{\link{StructureData}}, \code{\link{StructureData}}, \code{\link{StructureResults}}.
#' @details The delta-k values are calculated using the method used by Structure Harvestor \url{http://taylor0.biology.ucla.edu/structureHarvester/faq.html.}
#' @export
StructureCollection<-function(analyses) {
	## k is calued using method used in structureHarvestor
	# calculate summary
	# extract ln
	ln <- sapply(analyses, function(x) {sapply(x@results@replicates, lnprob.StructureReplicate)})
	ln.mean <- colMeans(ln)
	ln.sd <- apply(ln, 2, sd)
	# delta k calculations
	thisK <- ln.mean[2:(length(ln.mean)-1)]
	this.K.sd <- ln.sd[2:(length(ln.sd)-1)]
	prevK <- ln.mean[1:(length(ln.mean)-2)]
	nextK <- ln.mean[3:(length(ln.mean))]
	delta.k <- abs(nextK - (2 * thisK) + prevK) / this.K.sd
	# generate data.frame
	summary.data <- data.frame(
		k=sapply(analyses, function(x) {x@opts@MAXPOPS}),
		mean.lnprob=ln.mean,
		sd.lnprob=ln.sd,
		delta.k=c(NA, delta.k,NA),
		n.replicates=sapply(analyses, function(x) {length(x@results@replicates)})
	)
	# determine best
	best <- summary.data$k[which.max(summary.data$delta.k)]
	# create object
	x<-new("StructureCollection", summary=summary.data, best=best, analyses=analyses)
	validObject(x, test=FALSE)
	return(x)
}


#' @rdname run.Structure
#' @method run.Structure StructureData
#' @export
run.Structure.StructureData<-function(x, NUMRUNS=2, MAXPOPS=1:10, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, FREQSCORR=TRUE, UPDATEFREQ=max(floor(BURNIN+NUMREPS)/1000,1), M='Greedy', W=TRUE, S=FALSE, REPEATS=1000, dir=tempdir(), clean=TRUE, verbose=FALSE, threads=1, SEED=sample.int(1e5,NUMRUNS*length(MAXPOPS))) {
	return(run.Structure.list(
		list(x), NUMRUNS=NUMRUNS, MAXPOPS=MAXPOPS, BURNIN=BURNIN, NUMREPS=NUMREPS, NOADMIX=NOADMIX, ADMBURNIN=ADMBURNIN, FREQSCORR=FREQSCORR,
		M=M, W=W, S=S, REPEATS=REPEATS, dir=dir, clean=clean, verbose=verbose, threads=threads, SEED=SEED
	)[[1]])
}

#' @rdname run.Structure
#' @method run.Structure list
#' @export
run.Structure.list<-function(x, NUMRUNS=2, MAXPOPS=1:10, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, FREQSCORR=TRUE, UPDATEFREQ=max(floor(BURNIN+NUMREPS)/1000,1),
	M='Greedy', W=TRUE, S=FALSE, REPEATS=1000, dir=file.path(tempdir(), seq_along(x)), clean=TRUE, verbose=FALSE, threads=1, SEED=sample.int(1e5,NUMRUNS*length(MAXPOPS)*length(x))) {
	## init
	# init tests
	test_that('argument to MAXPOPS have at least 3 elements', expect_true(length(MAXPOPS) >= 3))
	sapply(x, expect_is, 'StructureData')
	test_that('arguments to x and dir must have equal lengths', expect_equal(length(dir), length(x)))
	test_that('arguments dir must be unique', expect_false(any(duplicated(dir))))
	# create table with parameters for all runs
	run.DF <- expand.grid(run=seq_len(NUMRUNS), k=MAXPOPS, species=seq_along(x))
	run.DF$SEED <- SEED
	run.DF$spp.dir <- dir[run.DF$species]
	run.DF$k.dir <- file.path(run.DF$spp.dir, run.DF$k)
	run.DF$i <- seq_len(nrow(run.DF))
	# create opts objects
	StructureOpts.LST<- dlply(run.DF, c('species'), function(z1) {
		z3 <- dlply(z1, c('species', 'k'), function(z2) {
			StructureOpts(NUMRUNS=NUMRUNS, MAXPOPS=z2$k[1], BURNIN=BURNIN, NUMREPS=NUMREPS, NOADMIX=NOADMIX, ADMBURNIN=ADMBURNIN,
			FREQSCORR=FREQSCORR, UPDATEFREQ=UPDATEFREQ, SEED=z2$SEED)
		})
		attributes(z3) <- NULL
		return(z3)
	})
	attributes(StructureOpts.LST) <- NULL
	clummp.opts <- ClumppOpts(M=M, W=W, S=S, REPEATS=REPEATS)
	## prelim
	# identify structure path
	structure.path <- system.file('bin', 'structure', package='structurer')
	# update permissions
	if (!grepl(basename(structure.path), 'win'))
		system(paste0('chmod 700 ',structure.path))
	# create dirs
	sapply(run.DF$k.dir, dir.create, showWarnings=FALSE, recursive=TRUE)
	# save data to each directory
	ret <- dlply(run.DF, c('species', 'k'), function(z) {
		write.StructureOpts(StructureOpts.LST[[z$species[1]]][[z$k[1]]],z$k.dir[1])
		write.StructureData(x[[z$species[1]]],file.path(z$k.dir[1], 'data.txt'))
	})
	## main
	# initialize cluster
	is.parallel.run <- (is.numeric(threads) && (threads>1)) | (is.character(threads) && (length(threads)>1))
	if (is.parallel.run) {
		clust <- makeCluster(threads, 'PSOCK')
		clusterEvalQ(clust, library(structurer))
		clusterExport(clust, c('run.DF', 'structure.path'), envir=environment())
		registerDoParallel(clust)
	}
	# generate replicates
	StructureReplicates.LST <- llply(
		seq_len(nrow(run.DF)),
		function(i) {
			o<-system(paste0('"', structure.path, '" ', '-m ',file.path(run.DF$k.dir[i], 'mainparams.txt'),' -e ',file.path(run.DF$k.dir[i], 'extraparams.txt'),' -K ',run.DF$k[i],' -L ',n.loci(x[[run.DF$species[i]]]),' -N ',n.samples(x[[run.DF$species[i]]]),' -i ',file.path(run.DF$k.dir[i], 'data.txt'),' -o ', paste0(run.DF$k.dir[i], '/output_run_',run.DF$run[i],'.txt'), ' -D ',run.DF$SEED[i], ' > ', run.DF$k.dir[i], '/structure_run_',run.DF$run[i],'_log.txt 2>&1'), intern=TRUE)
			# delete extra files created by structure
			if (file.exists('seed.txt')) unlink('seed.txt')
			if (file.exists(file.path(run.DF$k.dir[i],'seed.txt'))) unlink(file.path(run.DF$k.dir[i],'seed.txt'))
			# read results
			return(read.StructureReplicate(paste0(run.DF$k.dir[i], '/output_run_',run.DF$run[i],'.txt_f'), paste0(run.DF$k.dir[i], '/structure_run_',run.DF$run[i],'_log.txt')))
		},
		.parallel=is.parallel.run
	)
	# kill cluster
	if (is.parallel.run)
		clust <- stopCluster(clust)
	# combine into StrucutreCollection objects
	StructureCollection.LST <- dlply(run.DF, c('species'), function(z1) {
		curr.analyses <- dlply(z1, c('species', 'k'), function(z2) {
			StructureAnalysis(
				results=StructureResults(replicates=StructureReplicates.LST[z2$i], clummp.opts, dir=z1$k.dir[1]),
				opts=StructureOpts.LST[[z2$species[1]]][[z2$k[1]]],
				data=x[[z2$species[1]]]
			)
		})
		attributes(curr.analyses) <- NULL
		return(StructureCollection(curr.analyses))
	})
	attributes(StructureCollection.LST) <- NULL
	# return results
	return(StructureCollection.LST)
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
	return(n.pop(x@analyses[[1]]@opts))
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
sample.membership.StructureCollection <- function(x, threshold=NULL, k=x@best, ...) {
	pos <- which(sapply(x@analyses, function(z) {z@opts@MAXPOPS==k}))
	if (length(pos)==0)
		stop('Specified k is not in the StructureCollection object.')
	return(sample.membership(x@analyses[[pos]], threshold))
}

#' @rdname sample.names
#' @method sample.names<- StructureCollection
#' @export
`sample.names<-.StructureCollection` <- function(x, value) {
	return(
		new("StructureCollection",
			summary=x@summary, best=x@best,
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
		new("StructureCollection",
			summary=x@summary, best=x@best,
			analyses=lapply(x@analyses, `loci.names<-.StructureAnalysis`, value=value)
		)
	)
}

#' @method print StructureCollection
#' @rdname print
#' @export
print.StructureCollection=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureCollection object.\n\n")
	print(x@analyses[[which(sapply(x@analyses, function(z) {z@opts@MAXPOPS==x@best}))]], header=FALSE)
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureCollection',
	function(object)
		print.StructureCollection(object)
)

#' @rdname lnprob.plot
#' @method lnprob.plot StructureCollection
#' @export
lnprob.plot.StructureCollection <- function(x, main='', ...) {
	# extract data
	dat <- lapply(x@analyses, function(y) {
		lapply(y@results@replicates, function(z) {
			data.frame(k=as.character(ncol(z@matrix)), lnprob=z@lnprob)
		})
	})
	dat <- do.call(rbind, unlist(dat, recursive=FALSE))
	dat.summary <- data.frame(k=as.numeric(as.character(unique(dat$k))), lnprob=tapply(dat$lnprob, dat$k, mean))
	# make plot
	ggplot() +
		geom_violin(aes_string(x='k', y='lnprob'), data=dat) +
		geom_line(aes_string(x='k', y='lnprob'), data=dat.summary) +
		geom_point(aes_string(x='k', y='lnprob'), data=dat.summary) +
		theme_classic() +
		xlab('Number of populations (k)') +
		ylab('Estimated Ln prob of data') +
		theme(axis.line.x=element_line(), axis.line.y=element_line()) +
		ggtitle(main)
}

#' @rdname delta.k.plot
#' @method delta.k.plot StructureCollection
#' @export
delta.k.plot.StructureCollection <- function(x, main='Delta-K plot') {
	## make plot
	ggplot(data=x@summary) +
		geom_point(aes_string(y='delta.k', x='k')) +
		geom_line(aes_string(y='delta.k', x='k')) +
		theme_classic() +
		xlab('Number of populations') +
		ylab('Relative support (delta-K)') +
		theme(axis.line.x=element_line(), axis.line.y=element_line()) +
		ggtitle(main)
}

#' @rdname traceplot
#' @method traceplot StructureCollection
#' @export
traceplot.StructureCollection <- function(x, k=x@best, ...) {
	pos <- which(sapply(x@analyses, function(z) {z@opts@MAXPOPS==k}))
	if (length(pos)==0)
		stop('Specified k is not in the StructureCollection object.')
	return(traceplot.StructureAnalysis(x@analyses[[pos]]))
}

#' @rdname gelman.diag
#' @method gelman.diag StructureCollection
#' @export
gelman.diag.StructureCollection <- function(x, k=x@best, ...) {
	pos <- which(sapply(x@analyses, function(z) {z@opts@MAXPOPS==k}))
	if (length(pos)==0)
		stop('Specified k is not in the StructureCollection object.')
	return(gelman.diag.StructureAnalysis(x@analyses[[pos]]))
}

