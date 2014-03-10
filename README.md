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

    src/assignRG



Test data:
----------------------

For raw BAM files:

    src/assignRG -i testData/index.txt -o testData/demultiplexed.bam testData/todemultiplex.bam

For FASTQ files
    



   
