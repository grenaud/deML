# deML: Demultiplexer for Illumina Sequences (Multithreaded Version)

## Project Overview
This project is an extension of the original deML tool, which is designed to demultiplex Illumina sequencing reads based on their index sequences. The main enhancement in this version is the implementation of multithreading to improve performance on multi-core systems.

## Key Features
1. **Multithreading Support**: The tool now utilizes multiple threads to process reads in parallel, potentially significantly improving performance on multi-core systems.
2. **User-Defined Thread Count**: The number of threads to be used is set by the user. (Default number of threads is, 'All available hardware threads - 2')
3. **Support for BAM and FASTQ**: The tool can process both BAM and FASTQ input files.
4. **Performance Measurement**: Includes built-in timing and memory usage measurements to compare single-threaded and multi-threaded performance.
5. **Thread-Safe Operations**: Implemented mutex locks to ensure thread-safe access to shared resources.
6. **Preserved Read Order**: Maintains the exact same order of reads as the input file using a rank-based system.

## Major Changes
1. **Worker Thread Implementation**: Created a worker thread function that handles both single-end and paired-end reads.
2. **Thread Pool**: Implemented a simple thread pool using POSIX threads (pthread).
3. **Work Queue**: Added a queue to manage work items for the threads.
4. **Mutex Protection**: Added mutex locks for shared resources (unknownSeq, wrongSeq, conflictSeq) and for BAM writing operations.
5. **Performance Logging**: Implemented timing and memory usage measurements for single-thread and multi-thread executions.
6. **Ordered Output**: Added rank tracking and priority queue system to ensure reads are written in the same order as input.

## UsageY
To compile the project:
```bash
make
```

To Run the Tool:

```
./deML [options] <input_file>
```

Options:
```
-t, --threads <num>    # Set the number of threads to use (default: Total available threads - 2)
```
