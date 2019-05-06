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

//typedef std::map<unsigned, std::vector<unsigned> > InMemoryContainer;  
//typedef std::vector<Edge> InMemoryContainer;  
typedef std::map<unsigned, std::vector<unsigned> > InMemoryContainer;  
typedef InMemoryContainer::iterator InMemoryContainerIterator; 
typedef InMemoryContainer::const_iterator InMemoryConstIterator;
typedef std::map<unsigned, std::vector<unsigned> > LookUpTable;
typedef double real_t;

std::vector<double> writeBuf_times;
std::vector<double> flushResidues_times;
std::vector<double> infinimem_read_times;
std::vector<double> infinimem_write_times;
std::vector<uint64_t> localCombinedPairs; 
std::vector<unsigned> gWhere;
// *bndind, *bndptr;

void* combine(const unsigned& key, std::vector<unsigned>& to, const std::vector<unsigned>& from);

class InMemoryReductionState {
  public:
  std::vector<InMemoryConstIterator > begins;
  std::vector<InMemoryConstIterator > ends;

  InMemoryReductionState(unsigned size) : begins(size), ends(size) { }
};

class Partitioner
{
  public:
    void initg(unsigned nVertices, unsigned nMemParts, unsigned bSize, unsigned kItems, unsigned nReduceParts);
//    void coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const  unsigned int* mapSupRowstoRows);
    void writeInit(const unsigned tid);
    unsigned writeBuf(const unsigned tid, const unsigned to, const unsigned from, unsigned k, unsigned* adjncy, std::vector<unsigned>& where, std::vector<unsigned>& whereDst, bool firstRep);
//    graph_t* initsubgraph(const unsigned tid, const unsigned buffer);
    unsigned performWrite(const unsigned tid, const unsigned buffer, const unsigned to, const unsigned from, unsigned k, unsigned* adjncy, bool firstRep);
 
    void writeToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, const InMemoryContainer& inMemMap);

    void bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end);

    void setNum(const unsigned tid, std::vector<unsigned>& where, unsigned num);
    void gCopy(const unsigned tid, std::vector<unsigned>& gWhere, std::vector<unsigned>& where, const std::vector<unsigned> whereDst);
    void sCopy(const unsigned tid, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, std::vector<unsigned>& partitionBndind, std::vector<unsigned>& partitionBndPtr);
    void BNDInsert(const unsigned tid, unsigned n, std::vector<unsigned>& bndind, std::vector<unsigned> bndptr, unsigned i, bool first);
    void BNDDelete(const unsigned tid, unsigned n, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, unsigned i);

    void readInit(const unsigned tid);
    bool read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted);
    bool read(const unsigned tid);    

    bool getNextMinKey(InMemoryReductionState* state, InMemoryContainer* record);
    InMemoryReductionState initiateInMemoryReduce(unsigned tid);

    void ComputeBECut(const unsigned tid, const unsigned buffer, std::vector<unsigned>& where, std::vector<unsigned>& bndind, std::vector<unsigned> bndptr);
    void releaseInMemStructures();
    void shutdown();
    void flushBResidues(const unsigned tid);
    unsigned long long merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end);

    bool isPresent(unsigned k){
         return (k>=0);
    }

    bool getWrittenToDisk() { return writtenToDisk; }

    InMemoryContainer* readBufMap;
    LookUpTable* lookUpTable;
    std::set<unsigned>* fetchBatchIds;

    std::vector<unsigned>* partitionBndInd;
    std::vector<unsigned>* partitionBndPtr;
    std::vector<unsigned long long>* readNextInBatch;
    std::vector<bool>* batchesCompleted;
    std::vector<unsigned>* keysPerBatch;
  private:
   
    unsigned nVtces;   
    unsigned nRows;
    unsigned nCols;
    unsigned batchSize;  
    unsigned kBItems;  
//    unsigned nparts;  
    bool firstInit;
    IdType* totalPECuts; 
    IdType* nItems; 
    IdType* nCuts; 
    IdType* nEdges; 
    InMemoryContainer* outBufMap;  
//    std::vector<IdType>* xadj;  
//    std::vector<IdType>* adjncy;  

    IdType* totalKeysInFile;
    std::vector<pthread_mutex_t> locks;
    FileIO<RecordType> *io;  //GK
    bool writtenToDisk;
};

#endif // __PARTITIONER_H__

