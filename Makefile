ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LDFLAGS = -lz -lboost_system -lboost_filesystem -fopenmp -lboost_regex
export OPENMPFLAG = -fopenmp

all: 
	make release -C Src
test:
	make tests -C Src
clean:
	make clean -C Src
