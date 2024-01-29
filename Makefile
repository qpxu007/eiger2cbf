# This Makefile was contributed by Harry Powell (MRC-LMB)
# QX https://stackoverflow.com/questions/19901934/libpthread-so-0-error-adding-symbols-dso-missing-from-command-line

CBFLIB=/usr/lib/x86_64-linux-gnu/
CBFINC=/usr/include/cbflib/
HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/
CC=/usr/bin/gcc -O3

all:	
	${CC} -std=c99 -o eiger2cbf-g  -g  \
	-I${CBFINC} -I/usr/include/hdf5/serial/ -Wl,--copy-dt-needed-entries \
	-L${CBFLIB} -Ilz4 \
	eiger2cbf-g.c \
	lz4/lz4.c lz4/h5zlz4.c \
	bitshuffle/bshuf_h5filter.c \
	bitshuffle/bshuf_h5plugin.c \
	bitshuffle/bitshuffle.c \
	${HDF5LIB}/libhdf5_hl.a \
	${HDF5LIB}/libhdf5.a \
	-lcbf -lm -lpthread -lz -ldl


	${CC} -std=c99 -o eiger2cbf-omp  -fopenmp -g  \
	-I${CBFINC} -I/usr/include/hdf5/serial/ -Wl,--copy-dt-needed-entries \
	-L${CBFLIB} -Ilz4 \
	eiger2cbf-omp.c \
	lz4/lz4.c lz4/h5zlz4.c \
	bitshuffle/bshuf_h5filter.c \
	bitshuffle/bshuf_h5plugin.c \
	bitshuffle/bitshuffle.c \
	${HDF5LIB}/libhdf5_hl.a \
	${HDF5LIB}/libhdf5.a \
	-lcbf -lm -lpthread -lz -ldl 


clean: 
	rm -f *.o minicbf
