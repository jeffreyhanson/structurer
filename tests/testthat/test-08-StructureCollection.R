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
