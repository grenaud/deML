
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)
SIMPLEJSON= $(realpath ../SimpleJSON/)



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${SIMPLEJSON}/src/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:    assignRG

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

assignRG: assignRG.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)





clean:
	rm -f assignRG *.o

