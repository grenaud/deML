SHELL := /bin/bash



all: 	src/deML

src/deML: SimpleJSON/obj/JSONValue.o libgab/utils.o bamtools/build/src/api/libbamtools.a
	make -C src

SimpleJSON/obj/JSONValue.h:
	rm -rf SimpleJSON/
	git clone --recursive https://github.com/MJPA/SimpleJSON.git


SimpleJSON/obj/JSONValue.o: SimpleJSON/obj/JSONValue.h
	make -C SimpleJSON

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git


libgab/utils.o: bamtools/build/src/api/libbamtools.a  libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git


bamtools/build/src/api/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..



clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
