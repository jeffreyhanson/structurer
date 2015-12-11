#' @include structurer-internal.R misc.R generics.R
NULL
 
#' StructureData: An S4 class to represent input data for Structure
#'
#' This class stores input data for the Structure program.
#'
#' @slot matrix \code{matrix} with binary loci data. Each row is for a different sample; each column is for a different loci.
#' @slot loci.names \code{character} labels used to idenitfy loci.
#' @slot sample.names \code{character} labels used to identify labels.
#' @seealso \code{\link{StructureData}}.
#' @export
setClass(
	"StructureData",
	representation(
		matrix='matrix',
		loci.names='character',
		sample.names='character'
	),
	validity=function(object) {
		# matrix
		expect_true(all(object@matrix[] %in% c(0, 1, NA)))
		# check that dimensions > 0
		expect_true(all(dim(object@matrix)>0))
		# loci.names
		expect_is(object@loci.names,'character')
		expect_true(all(!is.na(object@loci.names)))
		# sample.names
		expect_is(object@loci.names,'character')
		expect_true(all(!is.na(object@sample.names)))
		# cross-object checks
		expect_equal(length(object@loci.names),ncol(object@matrix))
		expect_equal(length(object@sample.names),nrow(object@matrix))
		return(TRUE)
	}
)

#' Create StructureData object
#'
#' This function creates a new \code{StructureData} object.
#'
#' @param matrix \code{matrix} with binary loci data. Each row is for a different sample; each column is for a different loci.
#' @param loci.names \code{character} labels used to idenitfy loci.
#' @param sample.names \code{character} labels used to identify labels.
#' @seealso \code{\link{StructureData-class}}.
#' @export
StructureData<-function(matrix, loci.names, sample.names) {
	x<-new("StructureData", matrix=matrix, loci.names=loci.names, sample.names=sample.names)
	validObject(x, test=FALSE)
	return(x)
}
 

#' @method print StructureData
#' @rdname print
#' @export
print.StructureData=function(x, ..., header=TRUE) {
	if (header)
		cat("StructureData object.\n")
	cat('  samples:',nrow(x@matrix),'\n')
	cat('  loci:',ncol(x@matrix),'\n')
}

#' @rdname show
#' @method show StructureData
#' @export
setMethod(
	'show',
	'StructureData',
	function(object)
		print.StructureData(object)
)

#' Write data for Structure
#'
#' This function writes data for Structure.
#'
#' @param x \code{StructureData} object.
#' @param file \code{character} file path to write data.
#' @seealso \code{StructureData}.
#' @examples
#' x <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
#' write.StructureData(x, tempfile(fileext='txt'))
#' @export
write.StructureData <- function(x,file) {
	# set options
	oldSAF <- options()$stringsAsFactors
	options(stringsAsFactors=FALSE)
	# replace NAs with -9999
	mat <- x@matrix
	mat[which(is.na(mat))] <- -9999
	mat <- matrix(as.character(mat), ncol=ncol(mat))
	# save loci header
	write.table(
		as.data.frame(matrix(c('', as.character(x@loci.names), '', c(rep('0', n.loci(x)))), byrow=TRUE, nrow=2)),
		file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE
	)
	# write file
	write.table(
		cbind(data.frame(name=x@sample.names), as.data.frame(mat)),
		file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE
	)
	# reset options
	options(stringsAsFactors=oldSAF)
}


#' @rdname n.loci
#' @method n.loci StructureData
#' @export
n.loci.StructureData <- function(x) {
	return(ncol(x@matrix))
}

#' @rdname n.samples
#' @method n.samples StructureData
#' @export
n.samples.StructureData <- function(x) {
	return(nrow(x@matrix))
}

#' @rdname sample.names
#' @method sample.names StructureData
#' @export
sample.names.StructureData <- function(x) {
	return(x@sample.names)
}

#' @rdname sample.names
#' @method sample.names<- StructureData
#' @export
`sample.names<-.StructureData` <- function(x, value) {
	x@sample.names <- value
	return(x)
}

#' @rdname loci.names
#' @method loci.names StructureData
#' @export
loci.names.StructureData <- function(x) {
	return(x@loci.names)
}

#' @rdname loci.names
#' @method loci.names<- StructureData
#' @export
`loci.names<-.StructureData` <- function(x, value) {
	x@loci.names <- value
	return(x)
}


#' @rdname sample.subset
#' @method sample.subset StructureData
#' @export
sample.subset.StructureData <- function(x, samples) {
	if (is.character(samples))
		samples <- which(samples %in% x@sample.names)
	if (is.logical(samples))
		samples <- which(samples)
	return(
		StructureData(
			matrix=x@matrix[samples,,drop=FALSE],
			sample.names=x@sample.names[samples],
			loci.names=x@loci.names
		)
	)
}

#' @rdname loci.subset
#' @method loci.subset StructureData
#' @export
loci.subset.StructureData <- function(x, loci) {
	if (is.character(loci))
		loci <- which(loci %in% x@loci.names)
	if (is.logical(loci))
		loci <- which(loci)
	return(
		StructureData(
			matrix=x@matrix[,loci,drop=FALSE],
			sample.names=x@sample.names,
			loci.names=x@loci.names[loci]
		)
	)
}


#' Read FSTAT data for Structure
#'
#' This function reads FSTAT data for Structure.
#'
#' @param x \code{character} file path with data
#' @seealso \code{StructureData}.
#' @return \code{StructureData}.
#' @examples
#' x <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
#' @export
read.StructureData <- function(x) {
	# read first line
	meta <- as.numeric(strsplit(readLines(x, n=1), '\t')[[1]])
	# read in locinames
	locinames <- scan(x, skip=1, n=meta[2], what='character', quiet=TRUE)
	# read in samplenames
	dat <- fread(x, skip=meta[2]+1, data.table=FALSE)
	na.col <- apply(dat, 2, function(x) {all(is.na(x))})
	if (any(na.col))
		dat <- dat[,-which(na.col),drop=FALSE]
	# substitute 0, 1, NA
	mat <- as.matrix(dat[,-1,drop=FALSE])
	mat[which(mat[]==0)] <- NA
	mat[which(mat[]==1)] <- 0
	mat[which(mat[]==2)] <- 1
	# return object
	return(
		StructureData(matrix=mat, loci.names=locinames, sample.names=as.character(dat[[1]]))
	)
}



