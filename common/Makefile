all: bamfilter.out divbam.out

bamfilter.out : bamfilter.c Makefile
	gcc -O3 -o bamfilter.out bamfilter.c ../../ext/htslib-1.9-lowdep/libhts.a -I ../../ext/htslib-1.9-lowdep/ -pthread -lm -lz -lbz2 -llzma
divbam.out : divbam.c Makefile
	gcc -O3 -o divbam.out divbam.c ../../ext/htslib-1.9-lowdep/libhts.a -I ../../ext/htslib-1.9-lowdep/ -pthread -lm -lz -lbz2 -llzma

.PHONY clean:
	rm bamfilter.out divbam.out
	
