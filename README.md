# {#mainpage}

GO : Out of Core Graph Partitioning for large irregular graphs
=============================

GO is a single level graph partitioner that can successfully partition large graphs on a single machine in a memory constrained manner.
Instead of maintining multiple copies of the graph, GO performs just two passes over the entire input graph - partition creation pass that 
creates balanced partitions and partition refinement pass that reduces edgecuts. Both passes function in a memory constrained manner via 
disk-based processing. GO successfully partitions large graphs for which Mt-Metis runs out of memory.



The 'go' Executable
-----------------------------

To partition a graph with mt-Metis into 16 parts:

    go <graphfolderpath> <fileType> <nvertices> <hDegree> <nthreads> <nparts> <batchsize> <outputprefix>
    go  graph adjlist 32767 10000 32 32 0 test   

<filetype> - adjlist or edgelist
<graphfolderpath> - path for the input graph files.
<hDegree> - enter the maximum size for the adjacency list for the irregular graph. For regular graphs, enter zero.
<nparts> - number of graph parititions to be produced. Should be 2 or greater.
<batchsize> - this parameter gives you ability to limit the memory size used by algorith by specifying a smaller batchsize. 
              Enter zero for default option.

Script to split the input graph into multiple files is located under programs/scripts directory.  



<!-- Including GO API
-------------------------------

The file [mtmetis.h](@ref mtmetis.h) is the header that should be included
by external programs wishing link to mt-Metis. There are two high level
functions, mtmetis_partkway() for partitioning, and mtmetis_nd() for generating
orderings. At this time mt-Metis is highly experimental, and its
API is subject to change.
-->
