context('StructureAnalysis')

test_that('run.single.Structure (default parameters)', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=10, MAXPOPS=2, BURNIN=20, NUMREPS=30, NOADMIX=FALSE, ADMBURNIN=10, clean=FALSE)
	# methods
	print(sa)
	sa
	n.samples(sa)
	sample.membership(sa)
	n.pop(sa)
	loglik(sa)
	lnprob(sa)
	traceplot(sa)
	gelman.diag(sa)
})


test_that('run.single.Structure (NOADMIX=FALSE)', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=10, MAXPOPS=2, BURNIN=20, NUMREPS=30, NOADMIX=TRUE, ADMBURNIN=10, clean=FALSE)
})

test_that('run.single.Structure (FREQSCORR=FALSE)', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=10, MAXPOPS=2, BURNIN=20, NUMREPS=30, NOADMIX=FALSE, FREQSCORR=FALSE, ADMBURNIN=10, clean=FALSE)
})

test_that('run.single.Structure (NOADMIX=TRUE; FREQSCORR=FALSE)', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=10, MAXPOPS=2, BURNIN=20, NUMREPS=30, NOADMIX=TRUE, FREQSCORR=FALSE, ADMBURNIN=10, clean=FALSE)
})

test_that('run.single.Structure (parallel)', {
	# make Structure object
	sda <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sa <- run.single.Structure(sda, NUMRUNS=10, MAXPOPS=2, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, clean=FALSE, threads=2)
})

