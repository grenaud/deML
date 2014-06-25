SHELL := /bin/bash



all: 	src/deML

src/deML: SimpleJSON libgab bamtools
	make -C src

SimpleJSON:
	if [ ! -d SimpleJSON ]; then git clone --recursive https://github.com/MJPA/SimpleJSON.git; fi	
	make -C SimpleJSON

libgab: bamtools
	if [ ! -d libgab ]; then git clone --recursive https://github.com/grenaud/libgab.git; fi
	make -C libgab

bamtools:
	if [ ! -d bamtools ]; then git clone --recursive https://github.com/pezmaster31/bamtools.git; fi
	cd bamtools/ && mkdir build/  && cd build/ && cmake .. && make && cd ../..

clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
