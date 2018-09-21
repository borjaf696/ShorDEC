ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LDFLAGS = -lz -lemon -lboost_system -lboost_filesystem

all: 
	make release -C src
clean:
	make clean -C src
