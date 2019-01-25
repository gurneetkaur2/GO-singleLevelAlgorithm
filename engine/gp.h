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
    virtual void* coarsenG(const unsigned tid, const std::string& input) { };
//    virtual void* refine(const unsigned tid, const unsigned& rank, const std::vector<unsigned>& nbrs) = 0; 
     
    // System provided default; overridable by user
    virtual void run();
    void init(const std::string input, const unsigned coarseners, const unsigned refiners, const unsigned bSize, const unsigned kItems, const unsigned nParts);
    void writeBuf(const unsigned tid, const unsigned to, const std::vector<unsigned>& from);

    // Variables. Ideally, make these private and provide getters/setters.
 //   unsigned nVertices;
    unsigned nCoarseners;
    unsigned nRefiners;
    unsigned batchSize;  //Number of items in a batch
    unsigned kBItems;  //Top-k items to be fetched from in memory coarsen
    unsigned nparts;
    graph_t graph;
    std::vector<double> coarsen_times;
//    std::vector<double> refine_times;

    std::vector<std::string> fileList;

    pthread_barrier_t barCoarsen;
//    pthread_barrier_t barRefine;

    friend void* doCoarsen(void* arg);
 //   friend void* doRefine(void* arg);
 //   friend void* doInMemoryRefine(void* arg);

  private:
    // Variables
    std::string inputFileName;
    size_t bytesPerFile;
    Partitioner partitioner;

    //Methods
    void partitionInputForParallelReads();
};

