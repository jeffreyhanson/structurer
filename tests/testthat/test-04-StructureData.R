test_that('StructureData', {
	# simualate data
	locinames <- letters[1:5]
	samplenames <- LETTERS[1:10]
	matrix <- matrix(sample(0:1, size=50, replace=TRUE), ncol=5, nrow=10)
	matrix[sample(1:50, size=10, replace=FALSE)] <- NA
	# tests implicit
	x <- StructureData(
		matrix=matrix,
		loci.names=locinames,
		sample.names=samplenames
	)
	# test methods
	expect_equal(n.samples(x), nrow(x@matrix))
	expect_equal(n.loci(x), ncol(x@matrix))
	expect_equal(loci.names(x), x@loci.names)
	expect_equal(sample.names(x), x@sample.names)
	print(x)
	x
})

test_that('read.StructureData', {
	in.ssd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
})

test_that('write.StructureData', {
	path <- tempfile(fileext='.txt')
	ssd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	write.StructureData(ssd, path)
})

test_that('sample.subset.StructureData', {
	# load data
	bsd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	# subset data
	sample.names(bsd) <- as.character(seq_len(n.samples(bsd)))
	bsd2 <- sample.subset(bsd, as.character(1:10))
	# tests
	expect_equal(n.samples(bsd2), 10)
	expect_equal(sample.names(bsd2), as.character(1:10))
})

test_that('loci.subset.StructureData', {
	# load data
	bsd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	# subset data
	bsd2 <- loci.subset(bsd, 1:10)
	# tests
	expect_equal(n.loci(bsd2), 10)
	expect_equal(length(bsd2@loci.names), 10)
})
