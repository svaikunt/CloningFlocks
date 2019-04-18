GCC=g++
LDFLAGS= -lm -lgsl -lgslcblas -lfftw3 -I/usr/local/include/ -L /usr/local/lib/
OFLAGS=-o cloning.out
SOURCE=Cloning_run.cpp Langevin_dynamics.cpp

all:
	$(GCC) $(INC) $(LIB) $(SOURCE) $(LDFLAGS) $(OFLAGS)
