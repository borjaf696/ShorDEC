ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LDFLAGS = -lz -lemon

all: 
	make release -C src
clean:
	make clean -C src
