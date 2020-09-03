#ifndef __PARTITIONER_H__
#define __PARTITIONER_H__
#include "infinimem/fileIO.h"
#include <utility> //GK
#include <map>
#include <set>
#include <vector>
#include <stack>
#include <queue>
//#include "../programs/graph.h"

/*#ifdef USE_STRING_HASH
#define hashKey(str) stringHash(str)
#endif
*/
#ifdef USE_NUMERICAL_HASH
#define hashKey(number) (number)
#endif
/*#define HTLENGTH		((1<<11)-1)
#define UNMATCHED  -1
#define COARSEN_FRACTION	0.85 */	/* Node reduction between succesive coarsening levels */
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
typedef std::map<unsigned, unsigned > InMemTable;
typedef double real_t;
    std::vector<unsigned> gWhere;

std::vector<double> writeBuf_times;
std::vector<double> flushResidues_times;
std::vector<double> infinimem_read_times;
std::vector<double> infinimem_write_times;
std::vector<uint64_t> localCombinedPairs; 
std::vector<double> infinimem_cread_times;
std::vector<double> infinimem_cwrite_times;
// *bndind, *bndptr;

void* combine(const unsigned& key, std::vector<unsigned>& to, const std::vector<unsigned>& from);


class Partitioner
{

  public:
    void initg(unsigned nVertices, unsigned hDegree, unsigned nMemParts, unsigned bSize, unsigned kItems, unsigned nReduceParts, unsigned nrefiners);
//    void coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const  unsigned int* mapSupRowstoRows);
    //void writeInit(const unsigned tid);
    void writeInit();
    void writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hiDegree);
//    graph_t* initsubgraph(const unsigned tid, const unsigned buffer);
    void performWrite(const unsigned tid, const unsigned buffer, const unsigned to, const unsigned from);
 
    void writeToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, const InMemoryContainer& inMemMap);

    void bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end);

    void cWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end);
    void setNum(const unsigned tid, std::vector<unsigned>& where, unsigned num);
    void gCopy(const unsigned tid);
    void gCopy(const unsigned tid, std::vector<unsigned>& gWhere);
    void sCopy(const unsigned tid, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, std::vector<unsigned>& partitionBndind, std::vector<unsigned>& partitionBndPtr);

    void readInit(const unsigned tid);
    void readClear(const unsigned tid);
    bool read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted);
    bool read(const unsigned tid);    
    bool readInMem(const unsigned tid);    
    void readMemMap(const unsigned tid);
    void ctotalEdgeCuts(const unsigned tid);

//    bool refine(const unsigned tid, const IdType& totalCombined);
 
  //  bool getNextMinKey(InMemoryReductionState* state, InMemoryContainer* record);
    void initiateInMemoryRefine(unsigned tid);
    void addDVals(const unsigned tid);
 
    void ComputeBECut(const unsigned tid, const std::vector<unsigned>& where, LookUpTable& bndind, const InMemoryContainer& inMemMap);
    void ComputeBECut(const unsigned tid, const InMemoryContainer& inMemMap);
    void cWrite(const unsigned tid, unsigned noItems, InMemoryConstIterator end);
    unsigned countTotalPECut(const unsigned tid);
    unsigned maxPECut(const unsigned tid);
    void bRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, const bool ret);
    void bRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gWhere, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin, const bool ret);
    void inMemRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, const bool ret);
    void inMemRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gWhere, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin, const bool ret);


    void refineInit(const unsigned tid);
    bool refine(const unsigned tid);
    void cread(const unsigned tid);

    void computeDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax, const unsigned long long k);
    void updateDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax, unsigned src, unsigned dst);
    unsigned computeGain(const unsigned tid, const unsigned hipart, const unsigned whereMax, const unsigned long long k, const InMemoryContainer& inMemMap);
    unsigned computeGain(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin, const unsigned long long k, const InMemoryContainer& inMemMap);

 //   void changeWhere(const unsigned tid, const unsigned hipart, const unsigned chVtx );
   void writePartInfo(const unsigned tid, const unsigned hipart, const unsigned whereMax );
    bool checkPIDStarted(const unsigned tid, const unsigned hipart, const unsigned whereMax);

    void changeWhere(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gwhere, const unsigned maxVtx, const unsigned minVtx);
    void printParts(const unsigned tid, std::string outputPrefix);
    void clearMemorystructures(const unsigned tid);
    void releaseInMemStructures();
    void releaseReadPartStructures();
    void shutdown();
    void flushBResidues(const unsigned tid);
    unsigned long long merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end);

    bool isPresent(unsigned k){
         return (k>=0);
    }

 //   thread_local std::ofstream ofile;
 //   thread_local double stime;

    bool getWrittenToDisk() { return writtenToDisk; }
    bool getBndSet() { return bndSet; }
    bool getPartRefine() { return partRefine; }
    IdType getTotalPECuts(const unsigned tid) { return totalPECuts[tid]; }
    IdType setTotalPECuts(const unsigned tid) { totalPECuts[tid] = 0; }
    unsigned getTotalCuts(const unsigned tid) { if (tid == 0) return totalCuts; else return 0;  }
    void setTotalCuts(const unsigned tid) { if (tid == 0) totalCuts = 0; }
   
    pthread_barrier_t barRefinepart;
    pthread_barrier_t barEdgeCuts;
    pthread_barrier_t barCompute;

    InMemoryContainer* readBufMap;
    InMemoryContainer* refineMap;
    IdType* totalCombined;
    IdType* readNext; 
    //std::map<unsigned, unsigned>* bndindmap;
    LookUpTable* bndIndMap;
    LookUpTable* lookUpTable;
    InMemTable* dTable;
 //   InMemTable gainTable;

    std::set<unsigned>* fetchBatchIds;
    std::set<unsigned>* fetchPIds;
    std::map<unsigned, unsigned> hIds;

//    std::vector<unsigned>* partitionBndInd;
//    std::vector<unsigned>* partitionBndPtr;
    std::vector<unsigned long long>* readNextInBatch;
    std::vector<bool>* batchesCompleted;
    std::vector<bool>* pIdsCompleted;
    std::vector<bool> pIdCompleted;
    std::vector<unsigned>* keysPerBatch;
    std::map<unsigned, unsigned>* markMax;
    std::map<unsigned, unsigned>* markMin;
    std::map<unsigned, unsigned> pIdStarted;
  private:
   
    unsigned nVtces;   
    unsigned nRows;
    unsigned nCols;
    unsigned hiDegree;
    unsigned batchSize;  
    unsigned kBItems;  
    unsigned kAItems;  
    unsigned totalCuts; 
    //unsigned nRefiners; 
    unsigned nparts;  

    bool firstInit;
    IdType* totalPECuts; 
    IdType* totalPCuts; 
    IdType* nItems; 
    IdType* nEdges; 
    std::vector<unsigned>* where;
    InMemoryContainer* outBufMap;  
//    std::vector<IdType>* xadj;  
//    std::vector<IdType>* adjncy;  

    IdType* totalKeysInFile;
//    IdType* totalKeysRead;
    std::vector<pthread_mutex_t> locks;
    FileIO<RecordType> *io;  //GK
    FileIO<RecordType> *cio;  //GK
    bool writtenToDisk;
    bool partRefine;
    bool bndSet;
};

#endif // __PARTITIONER_H__

