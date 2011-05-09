ROOT_DIR = $(shell pwd)
export ROOT_DIR

.PHONY : install clean

all :
	cd m3dc1_lib ; make
	cd trace_lib ; make
	cd convert ; make

install :
	cd trace_lib ; make install
	cd convert ; make install

clean : 
	cd m3dc1_lib ; make clean
	cd trace_lib ; make clean
	cd convert ; make clean
