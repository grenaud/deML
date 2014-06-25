==========================================================
  deML: Maximum likelihood demultiplexing for NGS data
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
----------------------

deML is a program for maximum likelihood demultiplexing
of next-generation sequencing data. 


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






 
Format for input
----------------------

The two main inputs that are required are the 

1) Input sequences
2) Indices used in the experiment

----  Input sequences  ----

deML needs your sequences along with the sequenced index in either two formats:

a) BAM where the first index is specified as the XI tag and the quality as YI. 
   XJ and YJ contain the sequence and quality for the second index if it's there

b) fastq file where the forward reads are in one file, reverse in another file.
   The index must be in fastq format. If a second index is present, it must be in its own
   file. If you lost the quality information for the indices (they are present in the defline
   but no quality), you can still use deML but it will not provide optimal results. 
   Use a baseline quality for the bases of the index.


----  Indices used in the experiment  ----

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

However, these numbers must match those in the webform/config.json file. You can alwys modify that file if you use different indices.


Test data:
----------------------

For raw BAM files:

    src/deML -i testData/index.txt -o testData/demultiplexed.bam testData/todemultiplex.bam

For FASTQ files

    src/deML -i testData/index.txt -f testData/todemultiplex.fq1.gz  -r testData/todemultiplex.fq2.gz -if1 testData/todemultiplex.i1.gz  -if2 testData/todemultiplex.i2.gz   -o testData/demultiplexed.bam testData/todemultiplex.bam


Explanation for the scores
----------------------

deML works by computing the likelihood of stemming from potential samples and assigns a read to the most likely sample. To measure the confidence in the assignment, deML reports 3 different scores (Z_0, Z_1 and Z_2). If you have BAM as input, the output BAM file will include 3 different flags. 


Z_0  The likelihood for the top read group on a PHRED scale. Considering the sequence of this read group as the template to sequence and using the principle of independence between bases, this likelihood is computed as:

      Z_0= -10 * log_10{ PROD_{1...length indices} P(base read i|base template i)  }

The P(base read i|base template i) is given as (1-e_i) if both bases match or e_i otherwise where e_i is the predicted error rate for base i.

Since it is on a likelihood of assignment on a PHRED scale, the lower the Z_0 score, the higher the confidence. Around 0, the probability in assignment to that read group is around 1.



Z_1 Is the likelihood of misassignment on a PHRED scale. Let M be the number of potential read groups are present and let t be the read group with the highestZ_0. If Z_0_i is the Z_0 score for sample i, the Z_1 score is computed by:

                      |    SUM_i={1..M} (10^(Z_0_i)) - 10^(Z_0_t)   |
  Z_1= -10*log_10     |    --------------------------------------   |
                      |          SUM_i={1..M} (10^(Z_0_i))          |

The higher the Z_1 score, the lower the probability of misassignment and the higher the confidence. This score is not reported if only a single read group is found.



Z_2 Is a log odds ratio of mispairing on a logarithmic scale. It is only computed when two indices are used (P7 and P5). It is computed as:


                  |   SUM_{i!=j} ( P7_i) (P5_j)  )   |
Z_2  = 10* log_10 |   ------------------------------ |
                  |   SUM_{i==j} ( P7_i) (P5_j)  )   |


An approximation is made to speed up computations.  The higher the Z_2 score, the higher the odds ratio of mispairing and the lower the confidence. This score is not reported for runs with single indices and where the log ratio is less than 0 thus making the mispairing scenario unlikely.  

