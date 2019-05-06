#ifndef __HM_H_
#define __HM_H_
#include "fileIO.h"
#include "util.h"
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cassert>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <cmath>

#include <unistd.h> // sysconf(PAGE_SIZE)
#include <errno.h>
#include <stdlib.h> // getenv(BIG_MEMORY_LIMIT_MB)

#include <boost/lexical_cast.hpp>
#include<signal.h>

//------------------------------------------------------------------
#ifdef USE_METIS
#define METIS_BASE_PATH "../../../tools/partitioner/inputs/parts_"
#define GRAPH ".newGraph"
#define SRC ".src"
#define BOUNDARIES ".boundaries"

#else
#define GRAPH_BASE_PATH "../../../graphs/csr/"
#endif

extern const char* GRAPH_NAME; // define in main.cpp

//------------------------------------------------------------------
#define MAX_THREADS 128
#define TFLAG_IDLE 0
#define TFLAG_WORK 1
#define TFLAG_CLOSE 2
#define LFLAG (thread_flags[local_id])

const size_t mc_max_item_size = sizeof(char)*1024*128; // 128M.

template<typename VertexType, typename EdgeListType>
class HugeMem {
 public:
  static bool Init(int *argc, char ***argv, IdType wload);
  static void Shutdown();
  static fileIoReturn storeVertex(IdType id, VertexType* vertex);
  static fileIoReturn storeVertices(IdType start, IdType count, VertexType* vertices);
  static fileIoReturn storeEdgeList(IdType id, EdgeListType* edgeList);
  static fileIoReturn storeEdgeLists(IdType start, IdType count, EdgeListType* edgeLists);

  static fileIoReturn fetchVertex(IdType id, VertexType& v);
  static VertexType* fetchVertex(IdType id);
  static fileIoReturn fetchVertices(IdType id, unsigned count, VertexType* vertices);
  static fileIoReturn fetchEdgeList(IdType id, EdgeListType& el);
  static EdgeListType* fetchEdgeList(IdType id);
  static fileIoReturn fetchEdgeLists(IdType startId, IdType count, EdgeListType* edgeList);
  static void fetchNeighbors(std::vector<IdType>& nbrIds, VertexType** neighbors);
  static void fetchNeighbor(IdType nbrIds, VertexType* neighbor);
  static void updateVertex(IdType id);
  static void updateEdgeList(IdType id);

  static void workStart() { gettimeofday(&workStartTime, NULL); }
  static void workEnd() { gettimeofday(&workEndTime, NULL); }

  //--------------------
  static void SIGINT_handler(int signo) {
    if (signo == SIGINT) {
      if(signalled) //Ctrl-C twice instantly kills program
        exit(-1);
      signalled = true;
    }
  }

  //--------------------
  static unsigned long MemoryUsageMB() {
    typedef struct  {
      unsigned long size,res,share,text,lib,data,dt;
    } statm_t;

    statm_t result;
    const char* statm_path = "/proc/self/statm";

    FILE *f = fopen(statm_path,"r");
    if(!f) {
      perror(statm_path);
      abort();
    }

    if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld", &result.size, &result.res, &result.share, &result.text, &result.lib, &result.data, &result.dt)) {
      perror(statm_path);
      abort();
    }
    fclose(f);
    return (((result.size*pagesize)/1024UL)/1024UL);
  }

  inline static const std::string& getPreProcessDirFor(const char* g) {
    static const std::string dir(std::string("/tmp/hm/preprocess")
                          + std::string("/")
                          + std::string(g)
                          + std::string("/")
                          + boost::lexical_cast<std::string>(getCompThreads())
                          + std::string("/"));
    return dir;
  }

#ifdef USE_METIS
  inline static std::pair<IdType, IdType> get_bounds() { return std::make_pair(local_loop_begin, local_loop_end); }
  inline static const std::vector<IdType>& get_bounds_vector() { return bounds; }
  inline static const char* getMetisPathFor(std::string TYPE) {
    return (std::string(METIS_BASE_PATH) + boost::lexical_cast<std::string>(getCompThreads()) + "/" + std::string(GRAPH_NAME) + std::string(TYPE)).c_str();
  }
#else
  inline static const char* getVerticesPath(std::string PATH) {
    return (std::string(GRAPH_BASE_PATH) + std::string(GRAPH_NAME) + std::string(".numEdges")).c_str();
  }
  inline static const char* getGraphPath(std::string PATH) {
    return (std::string(GRAPH_BASE_PATH) + std::string(GRAPH_NAME)).c_str();
  }
#endif
#if defined(USE_METIS) || defined(GRAPHCHI)
  inline static const std::vector<IdType>& get_bounds_vector() { return bounds; }
#endif
#ifdef GRAPHCHI
  inline static void addInterval(IdType b) { bounds.push_back(b); }
#endif
  static void set_workload(unsigned long wload) { workload = wload; }
  inline static int get_tid() { return local_id; }
  inline static bool getDone() { return don; }
  inline static void done() { don = true; }
  inline static void notDone(EdgeListType *edgeList) {
    don = false;
#ifdef USE_ACTIVATIONS
    activateNeighborsForNextIteration(edgeList);
#endif // USE_ACTIVATIONS
  }
  inline static void reInitializeActivations();
#ifdef USE_ACTIVATIONS
  inline static IdType getActiveCount() {
    IdType ac = 0;
    for(int i=0; i<compthreads; i++)
      ac += nactive[i];
    return ac;
  }
#endif // USE_ACTIVATIONS

  static void waitForThread(int id);
  static void waitForAllThreads();
  static void startAllThreads(unsigned long lb, unsigned long le, void (*func)());
  static void closeAllThreads();
  static void startWorkTimes();

  static void setLoopIndexRange();
  static bool getLoopIndex(IdType *index);
  inline static int getThreadId() { return local_id; }
  inline static int getCompThreads() { return compthreads; }
  inline static int getWorkloadPerThread() { return workloadPerThread; }
  inline static IdType getWorkload() { return workload; }
  inline static int getBatchSize() { return vertexListSize; }
  inline static unsigned int getIterations() { return nIterations; }

 private:
  HugeMem(int *argc, char ***argv);
  HugeMem(HugeMem const&) {};
  HugeMem& operator=(HugeMem const&);
  ~HugeMem();

  static HugeMem *dsm;
  static time_t t;
  static timeval initTime;
  static timeval shutdownTime;
  static timeval workStartTime;
  static timeval workEndTime;
  static IdType workload;
#if defined(USE_METIS) || defined(GRAPHCHI)
  static std::vector<IdType> bounds;
#endif

  static bool signalled;
  static unsigned long pagesize;
  static unsigned int memoryLimitMB;

  static __thread VertexType* vertexList;
  static __thread EdgeListType* edgeLists;
#ifdef USE_THREAD_CACHE
  static __thread std::map<IdType, VertexType>* vertexMap;
#endif
  static __thread IdType vertexOffset;
  static __thread IdType edgeListsOffset;
  static __thread bool vertexListDirty;
  static __thread bool edgeListDirty;
  static IdType vertexListSize;
  static IdType edgeListsSize;
  static unsigned int vertexMapSize;
  static unsigned int nIterations;

  static void *Worker(void *argPtr);

  static int uberIteration;

  static int compthreads; // number of computing threads
  static int threads;   // number of threads

  static pthread_barrier_t bar; // barrier used to synchronize computation threads
  static __thread int local_id; // local thread ID
  static void (*thread_func)(); // work function for computation threads
  static pthread_t thread_handles[MAX_THREADS];
  static int thread_ids[MAX_THREADS];
  static int thread_flags[MAX_THREADS];
  static IdType loop_begin;
  static IdType loop_end;
  static __thread IdType local_loop_begin;   // beginning of local loop exec
  static __thread IdType local_loop_index;     // local loop index
  static __thread IdType local_loop_end;   // end of local loop exec

  static std::vector<double> fetch_times;
  static std::vector<unsigned> fetch_counts;
  static FileIO<VertexType> *vertexIO;
  static FileIO<EdgeListType> *edgeListIO;

  static bool don;

#ifdef USE_THREAD_CACHE
  static void insertToVertexMap(IdType id, VertexType* vertex);
#endif

  static IdType elementsPerThread();
  static IdType workloadPerThread;

#ifdef PART_MEMORY
  static VertexType* memoryVertices;
#endif

#ifdef USE_ACTIVATIONS
  static bool **activated;
  static IdType *nactive;
  static void activateNeighborsForNextIteration(EdgeListType *edgeList);
#endif // USE_ACTIVATIONS
};

#endif // __HM_H_


