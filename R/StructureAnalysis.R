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
sample.membership.StructureAnalysis <- function(x, threshold=NULL, ...) {
	return(sample.membership(x@results, threshold))
}

#' @rdname logLik
#' @method logLik StructureAnalysis
#' @export
logLik.StructureAnalysis <- function(object, ...) {
	return(logLik(object@results))
}

#' @rdname lnprob
#' @method lnprob StructureAnalysis
#' @export
lnprob.StructureAnalysis <- function(x, ...) {
	return(lnprob(x@results))
}


#' Run Structure anaylsis for a single MAXPOPS parameter.
#'
#' This function analyses data using Structure. It generates replicate runs for a single MAXPOPS parameter.
#'
#' @param x \code{StructureData} object.
#' @inheritParams StructureOpts
#' @inheritParams ClumppOpts
#' @inheritParams run.Structure
#' @param threads \code{numeric} number of threads to use for processing. Defaults to 1.
#'
#' @seealso \code{StructureData}, \code{StructureOpts}.
#' @examples
#' # run Structure using low number of iterations
#' dat <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
#' x <- run.single.Structure(dat, NUMRUNS=1, MAXPOPS=2, BURNIN=10,
#'	NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10)
#' @export
run.single.Structure<-function(x, NUMRUNS=2, MAXPOPS=2, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, FREQSCORR=TRUE,
	UPDATEFREQ=max(floor(BURNIN+NUMREPS)/1000,1), M='Greedy', W=TRUE, S=FALSE, REPEATS=1000, dir=tempdir(), clean=TRUE, verbose=FALSE, threads=1, SEED=sample.int(1e5,NUMRUNS))
{
	## initialization
	# argument checks
	opts <- StructureOpts(NUMRUNS=NUMRUNS, MAXPOPS=MAXPOPS, BURNIN=BURNIN, NUMREPS=NUMREPS, NOADMIX=NOADMIX, ADMBURNIN=ADMBURNIN, FREQSCORR=FREQSCORR, SEED=SEED, UPDATEFREQ=UPDATEFREQ)
	opts2 <- ClumppOpts(M=M, W=W, S=S, REPEATS=REPEATS)
	expect_is(x, 'StructureData')
	# identify structure path
	structure.path <- system.file('bin', 'structure', package='structurer')
	# update permissions
	if (!grepl(basename(structure.path), 'win'))
		system(paste0('chmod 700 ',structure.path))
	# sink opts and data to file
	write.StructureOpts(opts,dir)
	write.StructureData(x,file.path(dir, 'data.txt'))
	# setup cluster
	is.parallel.run <- (is.numeric(threads) && (threads>1)) | (is.character(threads) && (length(threads)>1))
	if (is.parallel.run) {
		clust <- makeCluster(threads,type='PSOCK')
		clusterEvalQ(clust, {library(structurer)})
		clusterExport(clust, c('structure.path','dir','MAXPOPS','x', 'verbose'), envir=environment())
		registerDoParallel(clust)
	}
	# run Structure
	ret <- StructureAnalysis(
			results=StructureResults(
				replicates=llply(
					seq_len(opts@NUMRUNS),
					function(i) {
						if (verbose) cat('\tstarting structure replicate ',i,'\n')
						o<-system(paste0('"', structure.path, '" ', '-m ',file.path(dir, 'mainparams.txt'),' -e ',file.path(dir, 'extraparams.txt'),' -K ',MAXPOPS,' -L ',n.loci(x),' -N ',n.samples(x),' -i ',file.path(dir, 'data.txt'),' -o ', paste0(dir, '/output_run_',i,'.txt'), ' -D ',opts@SEED[i], ' > ', dir, '/structure_run_',i,'_log.txt 2>&1'), intern=TRUE)
						# delete extra files created by structure
						if (file.exists('seed.txt')) unlink('seed.txt')
						if (file.exists(file.path(dir,'seed.txt'))) unlink(file.path(dir,'seed.txt'))
						# read results
						return(read.StructureReplicate(paste0(dir, '/output_run_',i,'.txt_f'), paste0(dir, '/structure_run_',i,'_log.txt')))
					},
					.parallel=is.parallel.run
				),
				opts2,
				dir=dir
			),
			data=x,
			opts=opts
		)
	# kill cluster
	if (is.parallel.run) {
		clust <- stopCluster(clust)
	}
	# if clean then delete files
	if (clean) unlink(dir)
	# return results
	return(ret)
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
	cat('  Gelman-Ruben R:', gelman.diag(x)$psrf[[1]], '\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'StructureAnalysis',
	function(object)
		print.StructureAnalysis(object)
)

#' @rdname traceplot
#' @method traceplot StructureAnalysis
#' @export
traceplot.StructureAnalysis <- function(x, parameter='Ln.Like', ...) {
	# create fake vars to avoid CRAN notes
	iteration <- NULL
	chain <- NULL
	param <- NULL
	# check that parameter is valid
	expect_length(parameter,1)
	expect_true(
		parameter %in% names(x@results@replicates[[1]]@mcmc)[-1],
		info=paste0('argument to parameter must be one of ',paste(names(x@results@replicates[[1]]@mcmc)[-1], collapse=','))
	)
	# extract logliks
	ll <- data.frame(iteration=x@results@replicates[[1]]@mcmc$Rep)
	ll <- cbind(ll, data.frame(sapply(x@results@replicates, function(y) {y@mcmc[[parameter]]})))
	ll <- ll[rowSums(ll[,-1])<0 || is.na(ll[,-1]),]
	if (ncol(ll)==2) names(ll)[2] <- 'X1'
	ll <- gather(ll, chain, param, -iteration)
	ll$chain <- as.factor(as.numeric(gsub('X', '', ll$chain, fixed=TRUE)))
	curr.burnin <- ifelse(sum(grepl('Admixture Burnin complete', x@results@replicates[[1]]@log, fixed=TRUE))>0, x@opts@ADMBURNIN, 0)+x@opts@BURNIN
	# make plot
	ggplot(data=ll, aes(x=iteration, y=param, color=chain)) +
		coord_cartesian(xlim=range(ll$iteration), ylim=range(ll$param, na.rm=TRUE)) +
		geom_rect(xmin=-10, xmax=curr.burnin,
			ymin=min(ll$param,na.rm=TRUE)*1.5, ymax=max(ll$param,na.rm=TRUE)*0.5,
			color='grey80', fill='grey80') +
		geom_line() + xlab('Iteration') + ylab(parameter) +
		theme_classic() + theme(axis.line.x=element_line(), axis.line.y=element_line())
}

#' @rdname gelman.diag
#' @method gelman.diag StructureAnalysis
#' @export
gelman.diag.StructureAnalysis <- function(x, ...) {
	# if only one replicate then cannot calculate diagnostics
	if (length(x@results@replicates)>1) {
		# extract -logliks
		ll <- sapply(x@results@replicates, function(y) {y@mcmc$Ln.Like})
		ll <- ll[rowSums(ll[-1,])<0,]
		ll2 <- list()
		for (i in seq_len(ncol(ll))) {
			curr.burnin <- ifelse(sum(grepl('Admixture Burnin complete', x@results@replicates[[i]]@log, fixed=TRUE))>0, x@opts@ADMBURNIN, 0)+x@opts@BURNIN
			ll2[[i]] <- mcmc(ll[,i], start=curr.burnin+1, thin=x@opts@UPDATEFREQ)
		}
		# return object
		return(coda::gelman.diag(mcmc.list(ll2), autoburnin=FALSE))
	}
	# if only one then return NA object
	warning('Gelman-Rubin statistics cannot be calculated for only a single Structure run')
	return(structure(list(psrf = structure(c(NA_real_, NA_real_), .Dim = 1:2, .Dimnames = list(NULL, c("Point est.", "Upper C.I."))),
		mpsrf = NULL), .Names = c("psrf", "mpsrf"), class = "gelman.diag"))	
}
