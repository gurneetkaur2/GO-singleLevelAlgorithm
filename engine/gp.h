#include "util.h"
#include "infinimem/fileIO.h"
#include "partitioner.h"


//RecordType* edgeLists = NULL;

void* doCoarsen(void* arg);

//void* doRefine(void* arg);

//void* doInMemoryRefine(void* arg);
//EdgeList* edgeLists = NULL;

class GraphParts
{
  public:
    virtual void* createMParts(const unsigned tid, const std::string& input, const std::string& type, const unsigned lineId, const unsigned hiDegree) = 0; 
    virtual void* refine(const unsigned tid, const unsigned& rank, const std::vector<unsigned>& nbrs) = 0; 
    virtual void* beforeRefine(const unsigned tid) { };
    virtual void* afterRefine(const unsigned tid, const unsigned nVertices) { };
         
    // System provided default; overridable by user
    virtual void run();
    void init(const std::string input, const std::string type, const unsigned nvertices, const unsigned hdegree, const unsigned nthreads, const unsigned nparts, const unsigned bSize, const unsigned kItems);

    void setRefiners(const unsigned refiners);
    void writeInit();
    void writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hidegree);
    bool read(const unsigned tid);
    void printParts(const unsigned tid, std::string outputPrefix);
    void readInit(const unsigned buffer);
    bool refine(const unsigned tid);
    void refineInit(const unsigned buffer);
    void subtractRefineTimes(const unsigned tid, const double stime);
    void ComputeBECut(const unsigned tid, const InMemoryContainer& inMemMap);
    void cWrite(const unsigned tid, unsigned noItems, InMemoryConstIterator end);
    unsigned countTotalPECut(const unsigned tid);

    // Variables. Ideally, make these private and provide getters/setters.
    unsigned nVertices;
    //unsigned nEdges;
    unsigned hDegree;
    unsigned nThreads;
    unsigned batchSize;  //Number of items in a batch
    unsigned kBItems;  //Top-k items to be fetched from in memory coarsen
    unsigned nParts;
    unsigned nrefiners;
    IdType numLines;
    std::vector<unsigned> end_read;
    std::vector<double> mparts_times;
    std::vector<double> refine_times;
    std::vector<double> aftr_refine_times;

    std::vector<std::string> fileList;

    pthread_barrier_t barMParts;
    pthread_barrier_t barRead;
    pthread_barrier_t barRefine;
    pthread_barrier_t barWriteInfo;
    pthread_barrier_t barAfterRefine;
    pthread_barrier_t barClear;

    friend void* doMParts(void* arg);
    friend void* doRefine(void* arg);
    friend void* doInMemoryRefine(void* arg);

  private:
    // Variables
    std::string inputFileName;
    std::string inType;
    size_t bytesPerFile;
    size_t linesPerThread;
    Partitioner partitioner;

    //Methods
    void partitionInputForParallelReads();
};

