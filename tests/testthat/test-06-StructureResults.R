test_that('StructureResults', {
	# create data
	dir<-tempdir()
	so <- StructureOpts(MAXPOPS=2, BURNIN=10, NUMREPS=10, NOADMIX=FALSE, ADMBURNIN=10, SEED=1:2)
	write.StructureOpts(so,dir)
	sd <- read.StructureData(system.file('extdata', 'example_fstat_aflp.dat', package='structurer'))
	sample.names(sd) <- as.character(seq_len(n.samples(sd)))
	write.StructureData(sd,file.path(dir, 'data.txt'))
	# identify bayescan path
	structure.path <- system.file('bin', 'structure', package='structurer')
	# run BayeScan
	system(paste0(structure.path, ' ', ' -m ',file.path(dir, 'mainparams.txt'),' -e ',file.path(dir, 'extraparams.txt'),' -K ',2,' -L ',n.loci(sd),' -N ',n.samples(sd),' -i ',file.path(dir, 'data.txt'),' -o ',file.path(dir, 'output.txt')))
	# try reading results back into R
	results <- StructureResults(replicates=list(read.StructureReplicate(file.path(dir, 'output.txt_f'))))
	# methods
	print(results)
	results
	n.pop(results)
	n.samples(results)
	loglik(results)
	sample.membership(results)
	sample.names(results)
})
