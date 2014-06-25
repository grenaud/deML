SHELL := /bin/bash



all: 	src/deML

src/deML: SimpleJSON/obj/JSONValue.o libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

SimpleJSON/obj/JSONValue.h:
	git clone --recursive https://github.com/MJPA/SimpleJSON.git


SimpleJSON/obj/JSONValue.o: SimpleJSON/obj/JSONValue.h
	make -C SimpleJSON

libgab/utils.h:
	git clone --recursive https://github.com/grenaud/libgab.git


libgab/utils.o: bamtools/lib/libbamtools.so  libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	git clone --recursive https://github.com/pezmaster31/bamtools.git


bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir build/  && cd build/ && cmake .. && make && cd ../..



clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
