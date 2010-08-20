all: gnu run

gnu:
	gcc -O3 main.c -o gscale -lm -DGNU
nopt:
	icc -O3 main.c -o scale_slow
opt:
	icc -O3 -fast main.c -o scale
prof:
	icc -O3 -opt-report -fast main.c -o scale

run:
	./gscale 100
