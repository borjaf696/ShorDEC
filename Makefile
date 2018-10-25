ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LDFLAGS = -lz -lemon -lboost_system -lboost_filesystem -fopenmp
export OPENMPFLAG = -fopenmp

all: 
	make release -C src
test:
	make tests -C src
clean:
	make clean -C src
