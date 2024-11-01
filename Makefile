SHELL := /bin/bash
CXXFLAGS += -std=c++17

all: 	src/deML

src/deML: SimpleJSON/obj/JSONValue.o libgab/libgab.a bamtools/build/src/api/libbamtools.a
	make -C src

SimpleJSON/obj/JSONValue.h:
	rm -rf SimpleJSON/
	git clone --recursive https://github.com/MJPA/SimpleJSON.git


SimpleJSON/obj/JSONValue.o: SimpleJSON/obj/JSONValue.h
	make -C SimpleJSON

libgab/libgab.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git


libgab/libgab.a: bamtools/build/src/api/libbamtools.a  libgab/libgab.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git


bamtools/build/src/api/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ && git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35 && mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi && make && cd ../..

clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
