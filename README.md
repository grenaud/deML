==========================================================
  deML: Maximum likelihood demultiplexing for NGS data
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
----------------------

deML is a program for maximum likelihood demultiplexing


Downloading:
----------------------

Go to https://github.com/grenaud/deML and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/deML.git"


Installation:
----------------------

1) Build Bamtools first:

    cd bamtools/   
    mkdir build/   
    cd build/
    cmake ..
    make 
    cd ../..

2) Build the submodules and main code by typing :

    make



running the program:
----------------------

To launch the program simply type:

    src/deML


Test data:
----------------------

For raw BAM files:

    src/deML -i testData/index.txt -o testData/demultiplexed.bam testData/todemultiplex.bam

For FASTQ files

    src/deML -i testData/index.txt -f testData/todemultiplex.fq1.gz  -r testData/todemultiplex.fq2.gz -if1 testData/todemultiplex.i1.gz  -if2 testData/todemultiplex.i2.gz  -o testData/demultiplex


-o testData/demultiplexed.bam testData/todemultiplex.bam
 
Format for input
----------------------

You can either specify the actual sequences used for multiplexing:

#Index1	Index2	Name
AATTCAA	CATCCGG	RG1
CGCGCAG	TCATGGT	RG2
AAGGTCT	AGAACCG	RG3
ACTGGAC	TGGAATA	RG4
AGCAGGT	CAGGAGG	RG5
GTACCGG	AATACCT	RG6
GGTCAAG	CGAATGC	RG7
AATGATG	TTCGCAA	RG8
AGTCAGA	AATTCAA	RG9

or you can also specify the raw indices:

#Index1	Index2	Name
341	33	RG1
342	34	RG2
343	35	RG3
344	36	RG4
345	37	RG5
346	38	RG6
347	39	RG7
348	40	RG8
349	41	RG9

   
