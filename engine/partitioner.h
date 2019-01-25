#ifndef __PARTITIONER_H__
#define __PARTITIONER_H__
#include "infinimem/fileIO.h"
#include <utility> //GK
#include <map>
#include <set>
#include <vector>
#include <stack>
#include <queue>
#include "../programs/graph.h"

/*#ifdef USE_STRING_HASH
#define hashKey(str) stringHash(str)
#endif
*/
#ifdef USE_NUMERICAL_HASH
#define hashKey(number) (number)
#endif
#define HTLENGTH		((1<<11)-1)
#define UNMATCHED  -1
#define COARSEN_FRACTION	0.85	/* Node reduction between succesive coarsening levels */
//-*-*-*-*-
/*template <typename KeyType, typename ValueType>
using InMemoryContainer = std::map<KeyType, std::vector<ValueType> >;

template <typename KeyType, typename ValueType>
using InMemoryContainerIterator = typename InMemoryContainer<KeyType, ValueType>::iterator; 

template <typename KeyType, typename ValueType>
using InMemoryContainerConstIterator = typename InMemoryContainer<KeyType, ValueType>::const_iterator;
//-*-*-*-*-
*/

typedef std::map<unsigned, std::vector<unsigned> > InMemoryContainer;  
typedef InMemoryContainer::iterator InMemoryContainerIterator; 
typedef InMemoryContainer::const_iterator InMemoryConstIterator;
typedef double real_t;

std::vector<double> writeBuf_times;
std::vector<double> flushResidues_times;
std::vector<double> infinimem_read_times;
std::vector<double> infinimem_write_times;
std::vector<uint64_t> localCombinedPairs; 

void* combine(const unsigned& key, std::vector<unsigned>& to, const std::vector<unsigned>& from);


class Partitioner
{
  public:
    void initg(unsigned nCoarseners, unsigned nRefiners, unsigned bSize, unsigned kItems, unsigned nParts);
//    void coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const  unsigned int* mapSupRowstoRows);
    void writeBuf(const unsigned tid, const IdType to, const std::vector<unsigned>& from); 
    graph_t* initsubgraph(const unsigned tid, const unsigned buffer);
    void performWrite(const unsigned tid, const unsigned buffer, const IdType to, const std::vector<unsigned>& from);


   // Coarsening related 
    graph_t* coarsen(const unsigned tid, graph_t *origG, unsigned CoarsenTo);
    IdType Match_RM(graph_t *graph);
    void CreateCoarseGraph(graph_t *graph, IdType cnvtxs, IdType *match);

   // Init Partition
    void AllocateKWayPartitionMemory(graph_t *graph);
    void InitKWayPartitioning(graph_t *graph);
    void PartGraphRecursive(IdType*, IdType*, IdType*, IdType*, IdType*, 
 			    IdType*, IdType*, unsigned int*, IdType*, IdType*);
    IdType MlevelRecursiveBisection(graph_t *graph, IdType nparts, 
          IdType *part, IdType fpart);
    IdType MultilevelBisect(graph_t *graph);
    void SplitGraphPart(graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph);

   // Graph related
    graph_t *CreateGraph(void);
    void InitGraph(graph_t *graph);
    graph_t* SetupCoarseGraph(graph_t *graph, IdType cnvtxs, IdType dovsize);
    graph_t* SetupGraph(IdType nvtxs, IdType ncon, IdType *xadj,
             IdType *adjncy, IdType *vwgt, IdType *vsize, IdType *adjwgt);
    void SetupGraph_label(graph_t *graph);
    void FreeGraph(graph_t **r_graph);

    void releaseMapStructures();
    void shutdown();
    void flushBResidues(const unsigned tid);
    unsigned long long merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end);

    bool getWrittenToDisk() { return writtenToDisk; }

  private:
   
    unsigned nvertices; 
    unsigned nRows;
    unsigned nCols;
    unsigned batchSize;  //GK
    unsigned kBItems;  //GK
    unsigned nparts;  //GK
    bool firstInit;
    //IdType* cTotalKeys; //GK
    IdType* nItems; //GK
    IdType* nEdges; //GK
    InMemoryContainer* outBufMap;  

    IdType* totalKeysInFile;
    std::vector<pthread_mutex_t> locks;
    FileIO<RecordType> *io;  //GK
    bool writtenToDisk;
};

#endif // __PARTITIONER_H__

