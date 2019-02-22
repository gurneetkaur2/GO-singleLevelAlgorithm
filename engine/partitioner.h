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
typedef std::map<unsigned, std::vector<unsigned> > LookUpTable;
typedef double real_t;

std::vector<double> writeBuf_times;
std::vector<double> flushResidues_times;
std::vector<double> infinimem_read_times;
std::vector<double> infinimem_write_times;
std::vector<uint64_t> localCombinedPairs; 

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
    void initg(unsigned nMemParts, unsigned bSize, unsigned kItems, unsigned nReduceParts);
//    void coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const  unsigned int* mapSupRowstoRows);
    void writeBuf(const unsigned tid, const IdType to, const std::vector<unsigned>& from); 
//    graph_t* initsubgraph(const unsigned tid, const unsigned buffer);
    void performWrite(const unsigned tid, const unsigned buffer, const IdType to, const std::vector<unsigned>& from);
 
    void writeToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, const InMemoryContainer& inMemMap);

    void bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end);

    void readInit(const unsigned tid);
    bool read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted);
    bool read(const unsigned tid);    

    bool getNextMinKey(InMemoryReductionState* state, InMemoryContainer* record);
    InMemoryReductionState initiateInMemoryReduce(unsigned tid);

    void releaseInMemStructures();
    void shutdown();
    void flushBResidues(const unsigned tid);
    unsigned long long merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end);

    bool getWrittenToDisk() { return writtenToDisk; }

    InMemoryContainer* readBufMap;
    LookUpTable* lookUpTable;
    std::set<unsigned>* fetchBatchIds;

    std::vector<unsigned long long>* readNextInBatch;
    std::vector<bool>* batchesCompleted;
    std::vector<unsigned>* keysPerBatch;
  private:
   
    unsigned nRows;
    unsigned nCols;
    unsigned batchSize;  //GK
    unsigned kBItems;  //GK
//    unsigned nparts;  //GK
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

