
GO : Out of Core Graph Partitioning for large irregular graphs
=============================

GO is a single level graph partitioner that can successfully partition large graphs on a single machine in a memory constrained manner.
Instead of maintining multiple copies of the graph, GO performs just two passes over the entire input graph - partition creation pass that 
creates balanced partitions and partition refinement pass that reduces edgecuts. Both passes function in a memory constrained manner via 
disk-based processing. GO successfully partitions large graphs for which Mt-Metis runs out of memory. GO algorithm details and results
can be found in the [paper](https://www.cs.ucr.edu/~gupta/research/Publications/Comp/NAS2021.pdf).



Compiling 
-----------------------------
Following compilers/libraries are required:

* C++14 Compiler.
* [Protocol Buffers](https://protobuf.dev/)

   apt-get install libprotobuf-dev protobuf-compiler libgoogle-perftools-dev


To compile, go to the programs directory and run:

    make

It will create a 'go' executable.



The 'go' Executable
-----------------------------

To partition a graph with go into 32 parts:

    go <graphfolderpath> <fileType> <nvertices> <hDegree> <nthreads> <nparts> <batchsize> <outputprefix>
    go  graph adjlist 32767 10000 32 32 0 test   

*filetype* - adjlist or edgelist.  
*nvertices* - number of vertices in the graph.  
*graphfolderpath* - path for the input graph files.  
*hDegree* - enter the maximum size for the adjacency list for the irregular graph. For regular graphs, enter zero.  
*nparts* - number of graph parititions to be produced. Should be 2 or greater.  
*batchsize* - this parameter gives you ability to limit the memory size used by algorith by specifying a smaller batchsize.
              Enter zero for default option. 
*outputprefix* - path where the partitioned output files will be stored. Example - C:\partitionDir\


Script to split the input graph into multiple files is located under programs/scripts directory.  



Reference
-----------------------------

To cite GO, use the following BibTeX entry:

    @INPROCEEDINGS{9605433,
     author={Kaur, Gurneet and Gupta, Rajiv},
     booktitle={2021 IEEE International Conference on Networking, Architecture and Storage (NAS)},
     title={GO: Out-Of-Core Partitioning of Large Irregular Graphs},
     year={2021},
     volume={},
     number={},
     pages={1-10},
     keywords={Runtime;Conferences;Memory management;Partitioning algorithms;irregular graphs;out-of-core processing;multilevel graph partitioning},
     doi={10.1109/NAS51552.2021.9605433}}




<!-- Including GO API
-------------------------------

The file [mtmetis.h](@ref mtmetis.h) is the header that should be included
by external programs wishing link to mt-Metis. There are two high level
functions, mtmetis_partkway() for partitioning, and mtmetis_nd() for generating
orderings. At this time mt-Metis is highly experimental, and its
API is subject to change.
-->
