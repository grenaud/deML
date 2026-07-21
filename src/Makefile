
HTSLIB= $(realpath ../htslib/)
LIBGAB= $(realpath ../libgab/)
SIMPLEJSON= $(realpath ../SimpleJSON/)



CXXFLAGS+= -O3 -Wall -std=c++11 -pthread -I${HTSLIB} -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${SIMPLEJSON}/src/
LDLIBS+= ${HTSLIB}/libhts.a ${LIBGAB}/libgab.a ${LIBGAB}/gzstream/libgzstream.a  -lm -lz -ldl
LDFLAGS+= -pthread

# For the dynamic build, -pthread (above) links libpthread correctly. A fully
# static binary, however, drops the pthread objects unless they are forced in
# with --whole-archive; without this std::thread's pthread_create resolves to a
# weak stub and crashes at runtime whenever -t>1 is used.
PTHREAD_DYNAMIC= -lpthread
PTHREAD_STATIC = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

CXX=g++ #-g -pg

all:    deML

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^

deML: deML.o PrefixTree.o RGAssign.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PTHREAD_DYNAMIC)

static: deML.o PrefixTree.o RGAssign.o ${SIMPLEJSON}/obj/JSON.o ${SIMPLEJSON}/obj/JSONValue.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} -static $(LDFLAGS) -o deML $^ $(LDLIBS) $(PTHREAD_STATIC)



clean:
	rm -f deML *.o

