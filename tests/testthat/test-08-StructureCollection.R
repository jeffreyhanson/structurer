context('StructureCollection')

test_that('run.Structure', {
	# make Structure object
	sd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sc <- run.Structure(sd, NUMRUNS=2, MAXPOPS=1:5, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, M='LargeKGreedy',S=TRUE,REPEATS=100)
	# methods
	print(sc)
	sc
	n.samples(sc)
	n.pop(sc)
	sample.membership(sc)
	sample.membership(sc, 0.75)
	summary(sc)
	loglik.plot(sc)
	delta.k.plot(sc)
	traceplot(sc, K=1)
	gelman.diag(sc, K=1)
})

test_that('run.Structure (parallel)', {
	# make Structure object
	sd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sc <- run.Structure(sd, NUMRUNS=2, MAXPOPS=1:5, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, M='LargeKGreedy',S=TRUE,REPEATS=100, threads=2)
})

test_that('run.Structure (parallel, list input)', {
	# make Structure object
	sd1 <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sd2 <- StructureData(matrix(sample(0:1, size=100, replace=TRUE), ncol=5), as.character(1:5), as.character(1:20))
	sc <- run.Structure(list(sd1,sd2), NUMRUNS=2, MAXPOPS=1:5, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, M='LargeKGreedy',S=TRUE,REPEATS=100, threads=2)
	# check that objects are correctly parsed
	expect_equal(n.loci(sd1), n.loci(sc[[1]]))
	expect_equal(n.loci(sd2), n.loci(sc[[2]]))
})
