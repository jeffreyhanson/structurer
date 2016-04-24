test_that('StructureOpts', {
	# tests implicit
	x <- StructureOpts(MAXPOPS=2, BURNIN=10000, NUMREPS=20000, NOADMIX=FALSE, ADMBURNIN=500, SEED=1:2)
	# methods
	print(x)
	x
})

