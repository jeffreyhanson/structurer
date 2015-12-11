test_that('run.Structure', {
	# make Structure object
	sd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sc <- run.Structure(sd, NUMRUNS=2, MAXPOPS=1:5, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, SEED=NA_real_)
	# methods
	print(sc)
	sc
	n.samples(sc)
	n.pop(sc)
	sample.membership(sc)
	summary(sc)
	loglik.plot(sc)
	delta.plot(sc)
})

 
