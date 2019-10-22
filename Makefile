ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

GATB=$(ROOT_DIR)/lib/gatb-core/gatb-core/
GATB_LIB=$(GATB)/build/lib
export GATBFLAGS = -I$(GATB)/src -I$(GATB)/build/include -I$(GATB)/thirdparty
export LDFLAGS = -lz -lboost_system -lboost_filesystem -fopenmp -lboost_regex -lboost_program_options
export OPENMPFLAG = -fopenmp

all: 
	make release -C Src
test:
	make tests -C Src
clean:
	make clean -C Src
