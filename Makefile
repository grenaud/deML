SHELL := /bin/bash



all: 	src/deML

src/deML: SimpleJSON/obj/JSONValue.o libgab/libgab.a libgab/gzstream/libgzstream.a htslib/libhts.a
	make -C src

SimpleJSON/obj/JSONValue.h:
	rm -rf SimpleJSON/
	git clone --recursive https://github.com/MJPA/SimpleJSON.git


SimpleJSON/obj/JSONValue.o: SimpleJSON/obj/JSONValue.h
	make -C SimpleJSON

libgab/libgab.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

# libgab.a is the bamtools-free core (libgab.o, FastQObj.o, FastQParser.o,
# PutProgramInHeaderHTS.o); building just this target avoids needing bamtools.
libgab/libgab.a: libgab/libgab.h
	make -C libgab libgab.a

# gzstream is a separate archive under libgab that deML links directly; the
# libgab.a target does not build it, so build it explicitly here.
libgab/gzstream/libgzstream.a: libgab/libgab.h
	make -C libgab gzstream/libgzstream.a

# htslib is vendored as a git submodule; clone (recursively, for htscodecs) and
# build the static libhts.a. Only zlib is needed for BAM/BGZF, so the optional
# bz2/lzma/libcurl backends are disabled to keep the link dependencies minimal.
htslib/Makefile:
	git submodule update --init --recursive htslib || \
	  git clone --recursive https://github.com/samtools/htslib.git htslib

htslib/libhts.a: htslib/Makefile
	cd htslib && autoreconf -i && ./configure --disable-bz2 --disable-lzma --disable-libcurl && $(MAKE) libhts.a



clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C src clean


.PHONY: all
