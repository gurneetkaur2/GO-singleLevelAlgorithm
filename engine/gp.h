#include "util.h"
#include "infinimem/fileIO.h"
#include "partitioner.h"


//RecordType* edgeLists = NULL;

void* doCoarsen(void* arg);

//void* doRefine(void* arg);

//void* doInMemoryRefine(void* arg);

class GraphParts
{
  public:
    virtual void* MParts(const unsigned tid, const std::string& input) { };
    virtual void* refine(const unsigned tid, const unsigned& rank, const std::vector<unsigned>& nbrs) = 0; 
    virtual void* beforeRefine(const unsigned tid) { };
    virtual void* afterRefine(const unsigned tid) { };
         
    // System provided default; overridable by user
    virtual void run();
    void init(const std::string input, const unsigned inmemparts, const unsigned reducers, const unsigned bSize, const unsigned kItems, const unsigned nParts);
    void writeBuf(const unsigned tid, const unsigned to, const std::vector<unsigned>& from);
    bool read(const unsigned tid);
    void readInit(const unsigned buffer);
    void subtractRefineTimes(const unsigned tid, const double stime);

    // Variables. Ideally, make these private and provide getters/setters.
 //   unsigned nVertices;
    unsigned nInMemParts;
    unsigned nReducers;
    unsigned batchSize;  //Number of items in a batch
    unsigned kBItems;  //Top-k items to be fetched from in memory coarsen
    unsigned nparts;
    std::vector<double> mparts_times;
    std::vector<double> refine_times;

    std::vector<std::string> fileList;

    pthread_barrier_t barMParts;
    pthread_barrier_t barRefine;

    friend void* doMParts(void* arg);
    friend void* doRefine(void* arg);
    friend void* doInMemoryRefine(void* arg);

  private:
    // Variables
    std::string inputFileName;
    size_t bytesPerFile;
    Partitioner partitioner;

    //Methods
    void partitionInputForParallelReads();
};

