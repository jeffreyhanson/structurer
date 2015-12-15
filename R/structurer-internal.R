#' @include dependencies.R
NULL

#' Write replicate data to file for clumpp
#'
#' This function writes replicate Structure runs to file for processing in clumpp.
#' @param x \code{list} with \code{StructureReplicate} objects.
#' @param file \code{character} file path to save data.
write.ClumppReplicates<-function(x, file) {
	lapply(
		seq_along(x),
		function(i) {
			write.table(
				do.call(
					cbind,
					list(
						data.frame(id=paste0(seq_len(nrow(x[[i]]@matrix)), ':')),
						as.data.frame(x[[i]]@matrix),
						data.frame(n=rep(1, nrow(x[[i]]@matrix)))
					)
				),
				file,
				sep=' ',
				row.names=FALSE,
				col.names=FALSE,
				append=(i!=1),
				quote=FALSE
			)
			cat('\n', file=file, append=TRUE)
		}
	)
	return(invisible())
}

#' Read replicate data to file for clumpp
#'
#' This function reads results from clumpp.
#' @param file \code{character} file path to save data.
#' @return \code{list} of \code{StructureReplicate} objects.
read.ClumppReplicates<-function(file) {
	x<-fread(file, data.table=FALSE)
	return(as.matrix(x[,c(-1, -ncol(x)),drop=FALSE]))
}


