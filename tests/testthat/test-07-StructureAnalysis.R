test_that('run.single.Structure', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=2, MAXPOPS=2, BURNIN=10000, NUMREPS=10000, NOADMIX=FALSE, ADMBURNIN=10, SEED=NA_real_, clean=FALSE)
	# methods
	print(sa)
	sa
	n.samples(sa)
	sample.membership(sa)
	n.pop(sa)
})

