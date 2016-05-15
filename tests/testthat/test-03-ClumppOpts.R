context('ClumppOpts')

test_that('ClumppOpts', {
	# tests implicit
	x <- ClumppOpts(M='Greedy', W=TRUE, S=FALSE, REPEATS=1000)
	# methods
	print(x)
	x
})

