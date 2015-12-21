#### Parameters
optimization = -O3

all: structure
	cp structure ../inst/bin/structure

clean:
	rm *.o structure

##### structure makefile
structure: structure.o params.o datain.o output.o ran.o mymath.o
	gcc -o structure structure.o params.o datain.o output.o ran.o mymath.o $(optimization) -lm

STRAT: STRAT.o params.o datain.o ran.o mymath.o
	gcc -o STRAT STRAT.o params.o datain.o ran.o mymath.o  -lm $(optimization) 

STRAT.o: STRAT.c
	gcc -c STRAT.c $(optimization)

structure.o: structure.c
	gcc -c structure.c $(optimization) 

output.o: output.c
	gcc -c output.c $(optimization) 

datain.o: datain.c
	gcc -c datain.c $(optimization) 

params.o: params.c
	gcc -c params.c $(optimization) 

ran.o: ran.c
	gcc -c ran.c $(optimization) 

mymath.o: mymath.c
	gcc -c mymath.c $(optimization) 

	
