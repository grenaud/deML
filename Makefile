SHELL := /bin/bash



all: 	src/deML

src/deML: SimpleJSON/obj/JSONValue.o libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

SimpleJSON/obj/JSONValue.o:
	git clone --recursive https://github.com/MJPA/SimpleJSON.git
	make -C SimpleJSON

libgab/utils.o: bamtools/lib/libbamtools.so
	git clone --recursive https://github.com/grenaud/libgab.git
	make -C libgab

bamtools/lib/libbamtools.so:
	git clone --recursive https://github.com/pezmaster31/bamtools.git
	cd bamtools/ && mkdir build/  && cd build/ && cmake .. && make && cd ../..

clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
