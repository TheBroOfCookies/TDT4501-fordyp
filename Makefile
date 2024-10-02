CC=mpicc
OUTPUT=#-DSAVE_RECEIVERS
CFLAGS_gcc=-ffast-math -fopenmp
#CFLAGS_icc=-fast -no-inline-max-size -no-inline-max-total-size -qopenmp -qopt-report=5
CFLAGS_icc=-qopenmp
CFLAGS+=-std=c99 -O3 ${OUTPUT} ${CFLAGS_${CC}} -DSAVE_RECEIVERS=ON
LDLIBS+=-lm
TARGETS=main_mpi

all: dirs ${TARGETS}
dirs:
	mkdir -p receivers
clean:
	-rm -fr ${TARGETS} receivers
