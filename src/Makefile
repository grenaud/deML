
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)
SIMPLEJSON= $(realpath ../SimpleJSON/)



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/src/ -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${SIMPLEJSON}/src/
LDLIBS+= ${BAMTOOLS}/build/src/api/libbamtools.a ${LIBGAB}/utils.o ${LIBGAB}/gzstream/libgzstream.a  -lm -lz

CXX=g++ #-g -pg

all:    deML

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

deML: deML.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)





clean:
	rm -f deML *.o

