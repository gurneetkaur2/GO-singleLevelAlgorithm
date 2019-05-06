#ifndef _HM_HPP_
#define _HM_HPP_

#include "hm.h"
#include "twoPhaseFileIO.hpp"
#include "onePhaseFileIO.hpp"

#include <fstream>
#include <cstring> // for memset()

#ifdef PART_MEMORY
template<typename VertexType, typename EdgeListType>
VertexType* HugeMem<VertexType, EdgeListType>::memoryVertices = NULL;
#endif

template<typename VertexType, typename EdgeListType>
bool HugeMem<VertexType, EdgeListType>::signalled = false;

template<typename VertexType, typename EdgeListType>
unsigned long HugeMem<VertexType, EdgeListType>::pagesize = 0;

// Define the declared static variables.
template<typename VertexType, typename EdgeListType>
HugeMem<VertexType, EdgeListType> *HugeMem<VertexType, EdgeListType>::dsm = NULL;

template<typename VertexType, typename EdgeListType>
time_t HugeMem<VertexType, EdgeListType>::t;

template<typename VertexType, typename EdgeListType>
timeval HugeMem<VertexType, EdgeListType>::initTime;

template<typename VertexType, typename EdgeListType>
timeval HugeMem<VertexType, EdgeListType>::shutdownTime;

template<typename VertexType, typename EdgeListType>
timeval HugeMem<VertexType, EdgeListType>::workStartTime;

template<typename VertexType, typename EdgeListType>
timeval HugeMem<VertexType, EdgeListType>::workEndTime;

template<typename VertexType, typename EdgeListType>
int HugeMem<VertexType, EdgeListType>::compthreads = 0;

template<typename VertexType, typename EdgeListType>
int HugeMem<VertexType, EdgeListType>::threads = 0;

template<typename VertexType, typename EdgeListType>
pthread_barrier_t HugeMem<VertexType, EdgeListType>::bar;

template<typename VertexType, typename EdgeListType>
__thread int HugeMem<VertexType, EdgeListType>::local_id = -1;

template<typename VertexType, typename EdgeListType>
void (*HugeMem<VertexType, EdgeListType>::thread_func)() = NULL;

template<typename VertexType, typename EdgeListType>
pthread_t HugeMem<VertexType, EdgeListType>::thread_handles[MAX_THREADS];

template<typename VertexType, typename EdgeListType>
int HugeMem<VertexType, EdgeListType>::thread_ids[MAX_THREADS];

template<typename VertexType, typename EdgeListType>
int HugeMem<VertexType, EdgeListType>::thread_flags[MAX_THREADS];

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::loop_begin = 0;

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::loop_end = 0;

template<typename VertexType, typename EdgeListType>
__thread IdType HugeMem<VertexType, EdgeListType>::local_loop_begin = 0;

template<typename VertexType, typename EdgeListType>
__thread IdType HugeMem<VertexType, EdgeListType>::local_loop_index = 0;

template<typename VertexType, typename EdgeListType>
__thread IdType HugeMem<VertexType, EdgeListType>::local_loop_end = 0;

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::workload = 0;

#if defined(USE_METIS) || defined(GRAPHCHI)
template<typename VertexType, typename EdgeListType>
std::vector<IdType>  HugeMem<VertexType, EdgeListType>::bounds;
#endif

template<typename VertexType, typename EdgeListType>
int HugeMem<VertexType, EdgeListType>::uberIteration = -1;

template<typename VertexType, typename EdgeListType>
std::vector<double> HugeMem<VertexType, EdgeListType>::fetch_times;

template<typename VertexType, typename EdgeListType>
std::vector<unsigned> HugeMem<VertexType, EdgeListType>::fetch_counts;

template<typename VertexType, typename EdgeListType>
FileIO<VertexType>* HugeMem<VertexType, EdgeListType>::vertexIO = NULL;

template<typename VertexType, typename EdgeListType>
FileIO<EdgeListType>* HugeMem<VertexType, EdgeListType>::edgeListIO = NULL;

template<typename VertexType, typename EdgeListType>
bool HugeMem<VertexType, EdgeListType>::don = false;

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::vertexListSize = 0;

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::edgeListsSize = 0;

template<typename VertexType, typename EdgeListType>
__thread VertexType* HugeMem<VertexType, EdgeListType>::vertexList = NULL;

template<typename VertexType, typename EdgeListType>
__thread EdgeListType* HugeMem<VertexType, EdgeListType>::edgeLists = NULL;

template<typename VertexType, typename EdgeListType>
__thread IdType HugeMem<VertexType, EdgeListType>::vertexOffset = MAX_ID_VALUE;

template<typename VertexType, typename EdgeListType>
__thread IdType HugeMem<VertexType, EdgeListType>::edgeListsOffset = MAX_ID_VALUE;

template<typename VertexType, typename EdgeListType>
__thread bool HugeMem<VertexType, EdgeListType>::vertexListDirty = false;

template<typename VertexType, typename EdgeListType>
__thread bool HugeMem<VertexType, EdgeListType>::edgeListDirty = false;

template<typename VertexType, typename EdgeListType>
unsigned int HugeMem<VertexType, EdgeListType>::vertexMapSize = 0;

template<typename VertexType, typename EdgeListType>
unsigned int HugeMem<VertexType, EdgeListType>::nIterations = UINT_MAX;

#ifdef USE_THREAD_CACHE
template<typename VertexType, typename EdgeListType>
__thread std::map<IdType, VertexType>* HugeMem<VertexType, EdgeListType>::vertexMap = false;
#endif

template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::workloadPerThread = 0;

#ifdef USE_ACTIVATIONS
template<typename VertexType, typename EdgeListType>
bool **HugeMem<VertexType, EdgeListType>::activated = NULL;

template<typename VertexType, typename EdgeListType>
IdType *HugeMem<VertexType, EdgeListType>::nactive = NULL;
#endif // USE_ACTIVATIONS

//------------------------------------------------------------------
  static void remove_arg(int *argc, char **argv, int i) {
    if (i<*argc)
      (*argc) --;

    for(;i<*argc;i++) {
      argv[i] = argv[i+1];
    }
  }

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
bool HugeMem<VertexType, EdgeListType>::Init(int *argc, char ***argv, IdType wload) {
  gettimeofday(&initTime, NULL);

  workload = wload;
  errno = 0; // global system error variable.
  pagesize = sysconf(_SC_PAGESIZE);
  assert( errno == 0 );

  assert(signal(SIGINT, SIGINT_handler) != SIG_ERR);

  eng.seed((unsigned int)time(NULL));

  if(dsm) return true;
  dsm = new HugeMem<VertexType, EdgeListType>(argc, argv);
  return dsm ? true : false;
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::Shutdown() {
  gettimeofday(&shutdownTime, NULL);
  fprintf(stderr, "Data-WallTime: %.3lf\n", timevalToDouble(workStartTime)-timevalToDouble(initTime));
  fprintf(stderr, "Work-WallTime: %.3lf\n", timevalToDouble(workEndTime)-timevalToDouble(workStartTime));
  fprintf(stderr, "Run-WallTime: %.3lf\n", timevalToDouble(shutdownTime)-timevalToDouble(initTime));
  delete dsm;
  dsm = NULL;
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
HugeMem<VertexType, EdgeListType>::HugeMem(int *argc, char ***argv) {
  // get the number of cores on the current machine
  threads = sysconf( _SC_NPROCESSORS_ONLN );
  assert(threads <= MAX_THREADS);
  compthreads = threads;

  // parse the argv to get the thread number
  for(int i=1; i<*argc; i++) {
    if ( strncmp((*argv)[i], "--prodsm-compthreads=", 21) == 0) {
      sscanf((*argv)[i]+21, "%d", &compthreads);
      fprintf(stderr, "set computation threads = %d\n", compthreads);
      threads = compthreads;
      assert(threads <= MAX_THREADS);
      remove_arg(argc, *argv,i);
      i--;
    }
    else if  ( strncmp((*argv)[i], "--prodsm-vertexlist=", 20) == 0) {
      sscanf((*argv)[i] + 20, "%lu", &vertexListSize);
      fprintf(stderr, "set vertexListSize to: %lu\n", vertexListSize);
      remove_arg(argc, *argv,i);
      i--;
    }
    else if  ( strncmp((*argv)[i], "--prodsm-edgelist=", 18) == 0) {
      sscanf((*argv)[i] + 18, "%lu", &edgeListsSize);
      fprintf(stderr, "set edgeListsSize to: %lu\n", edgeListsSize);
      remove_arg(argc, *argv,i);
      i--;
    }
    else if  ( strncmp((*argv)[i], "--prodsm-vertexmap=", 19) == 0) {
      sscanf((*argv)[i] + 19, "%u", &vertexMapSize);
      fprintf(stderr, "set vertexMapSize to: %u\n", vertexMapSize);
      remove_arg(argc, *argv,i);
      i--;
    }
    else if  ( strncmp((*argv)[i], "--prodsm-iters=", 15) == 0) {
      sscanf((*argv)[i] + 15, "%u", &nIterations);
      fprintf(stderr, "set nIterations to: %u\n", nIterations);
      remove_arg(argc, *argv,i);
      i--;
    }
  }

#ifdef USE_ACTIVATIONS
  fprintf(stderr, "USE_ACTIVATIONS: true\n");
#endif

  assert(vertexListSize != 0);
  assert(edgeListsSize != 0);
  assert(vertexMapSize != 0);

  for(int i=0; i<compthreads; ++i) {
    fetch_counts.push_back(0);
    fetch_times.push_back(0.0);
  }

  // create barrier for computation threads
  pthread_barrier_init(&bar, NULL, compthreads);

  // initialize the current thread, which is also used as a computation thread
  thread_ids[0] = 0;
  thread_flags[0] = TFLAG_WORK;
  local_id = 0;
  vertexList = new VertexType[vertexListSize];
  edgeLists = new EdgeListType[edgeListsSize];
#ifdef USE_THREAD_CACHE
  vertexMap = new std::map<IdType, VertexType>();;
#endif

#ifndef GRAPHCHI
  // create threads
  for(int i=1; i<compthreads; i++) {
    thread_ids[i] = i;
    thread_flags[i] = TFLAG_IDLE;
    pthread_create(&(thread_handles[i]), NULL, Worker, (void *)&thread_ids[i]);
  }
#endif // GRAPCHI

#ifdef USE_METIS
  std::ifstream boundaries(getMetisPathFor(BOUNDARIES));
  assert(boundaries.is_open());
  IdType b;
  while( boundaries >> b) {
    //efprintf(stderr, "F: %s, B: %lu\n", getMetisPathFor(BOUNDARIES), b);
    bounds.push_back(b);
  }
#endif

  workloadPerThread = elementsPerThread();

#if defined(USE_METIS) || defined(GRAPHCHI)
  vertexIO = new OnePhaseFileIO<VertexType>("/tmp/hm/data/vertices/", compthreads, workloadPerThread, bounds);
  edgeListIO = new TwoPhaseFileIO<EdgeListType>("/tmp/hm/data/edgeLists/", compthreads, workloadPerThread, bounds);
#else 
  vertexIO = new OnePhaseFileIO<VertexType>("/tmp/hm/data/vertices/", compthreads, workloadPerThread);
  edgeListIO = new TwoPhaseFileIO<EdgeListType>("/tmp/hm/data/edgeLists/", compthreads, workloadPerThread);
#endif

#ifdef PART_MEMORY
  memoryVertices = new VertexType[PART_MEMORY];
#endif

#ifdef USE_ACTIVATIONS
  // vertex activations
  activated = new bool*[2]; // one for current iteration, other for next iteration
  activated[0] = new bool[workload];
  activated[1] = new bool[workload];
  
  nactive = new IdType[compthreads];
  for(int i=0; i<compthreads; i++) nactive[i] = 0;
  reInitializeActivations();
#endif // USE_ACTIVATIONS

  time(&t);
  fprintf(stderr, "HugeMem init-ed: %s\n", ctime(&t));
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
HugeMem<VertexType, EdgeListType>::~HugeMem() {
  unsigned fc = 0;
  double ft = 0.0;
  for(int i=0; i<compthreads; ++i) {
    ft += fetch_times[i];
    fc += fetch_counts[i];
  }

  fprintf(stderr, "Total fetches: %u, Total fetch times: %.3lf, Average time/fetch: %.3lf\n", fc, ft, ft/fc);

  delete vertexIO;
  delete edgeListIO;

#ifdef PART_MEMORY
  delete [] memoryVertices;
#endif

  closeAllThreads();
  pthread_barrier_destroy(&bar);
}

#ifdef USE_METIS
//------------------------------------------------------------------
// return the length of the longest interval from bounds
template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::elementsPerThread() {
  IdType max = bounds[0]-1; // size of 0-th interval
  for(unsigned i=1; i<bounds.size(); i++)
    if(max < bounds[i]-bounds[i-1])
      max = bounds[i]-bounds[i-1];
  if(max < workload-bounds[bounds.size()]) // last interval
    max = workload-bounds[bounds.size()];
  return max;
}
#else
//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
IdType HugeMem<VertexType, EdgeListType>::elementsPerThread() {
  IdType rl = workload / compthreads;
  return rl + (workload % compthreads);
}
#endif

//------------------------------------------------------------------
#ifdef USE_THREAD_CACHE
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::insertToVertexMap(IdType id, VertexType* vertex) {
  if(vertexMap->size() == vertexMapSize)
     vertexMap->erase(vertexMap->end()); // TODO: need a better heuristic
   vertexMap->insert(std::pair<IdType, VertexType>(id, *vertex));
}
#endif

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::storeVertex(IdType id, VertexType* vertex) {
  return vertexIO->file_set(id, *vertex);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::storeVertices(IdType start, IdType count, VertexType* vertices) {
  //efprintf(stderr, "TID: %d, start: %llu, count: %llu\n", get_tid(), start, count);
  return vertexIO->file_set_batch(start, count, vertices);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::storeEdgeList(IdType id, EdgeListType* edgeList) {
  return edgeListIO->file_set(id, *edgeList);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::storeEdgeLists(IdType start, IdType count, EdgeListType* edgeLists) {
  return edgeListIO->file_set_batch(start, count, edgeLists);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::fetchVertex(IdType id, VertexType& v) {
  return vertexIO->file_get(id, v);
}

#ifdef USE_METIS
#define BATCH_SIZE(start, listSize) min(bounds[local_id]-start, listSize)
#else
#define BATCH_SIZE(start, listSize) min(min((local_id + 1) * workloadPerThread, workload) - start, listSize)
#endif

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
VertexType* HugeMem<VertexType, EdgeListType>::fetchVertex(IdType id) {
  if(id >= vertexOffset && id < vertexOffset + vertexListSize)
    return &(vertexList[id-vertexOffset]);

  if(vertexListDirty)
    assert(vertexIO->file_set_batch(vertexOffset, BATCH_SIZE(vertexOffset, vertexListSize), vertexList) == FILEIO_SUCCESS);

  vertexListDirty = false;
  vertexOffset = id;
  assert(vertexIO->file_get_batch(vertexOffset, BATCH_SIZE(vertexOffset, vertexListSize), vertexList) == FILEIO_SUCCESS);

  //efprintf(stderr, "Vertex %llu is %llu\n", id, vertexList[0].distance);
  return &vertexList[0];
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::fetchVertices(IdType start, unsigned count, VertexType* vertices) {
  return vertexIO->file_get_batch(start, count, vertices);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::fetchEdgeList(IdType id, EdgeListType& el) {
  return edgeListIO->file_get(id, el);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
fileIoReturn HugeMem<VertexType, EdgeListType>::fetchEdgeLists(IdType startId, IdType count, EdgeListType* edgeList) {
  return edgeListIO->file_get_batch(startId, count, edgeList);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
EdgeListType* HugeMem<VertexType, EdgeListType>::fetchEdgeList(IdType id) {
  if(id >= edgeListsOffset && id < edgeListsOffset + edgeListsSize)
    return &(edgeLists[id-edgeListsOffset]);

  // if(edgeListDirty)
  //   assert(edgeListIO->file_set_batch(edgeListsOffset, BATCH_SIZE(edgeListsOffset, vertexListSize), edgeLists) == FILEIO_SUCCESS);

  edgeListDirty = false;
  edgeListsOffset = id;
  assert(edgeListIO->file_get_batch(edgeListsOffset, BATCH_SIZE(edgeListsOffset, edgeListsSize), edgeLists) == FILEIO_SUCCESS);
  return &edgeLists[0];
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::fetchNeighbors(std::vector<IdType>& nbrIds, VertexType** neighbors) {
  unsigned size = nbrIds.size();
  for(unsigned i=0; i<size; ++i) {
    if(nbrIds[i] >= vertexOffset && nbrIds[i] < vertexOffset + vertexListSize)
      *(neighbors + i) = &vertexList[nbrIds[i]-vertexOffset];
#ifdef USE_THREAD_CACHE
    else {
      typename std::map<IdType, VertexType>::iterator it = vertexMap->find(nbrIds[i]);
      if(it != vertexMap->end())
        *(neighbors + i) = it->second;
      else {
#endif
        assert(vertexIO->file_get(nbrIds[i], 1, *(neighbors + i)));
#ifdef USE_THREAD_CACHE
        HugeMem<VertexType, EdgeListType>::insertToVertexMap(nbrIds[i], neighbors + i);
      }
    }
#endif
  }
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::fetchNeighbor(IdType nbrId, VertexType* neighbor) {
#ifdef PART_MEMORY
  if(nbrId >= 0 && nbrId < PART_MEMORY) {
    *neighbor = memoryVertices[nbrId];
    return;
  }
  //efprintf(stderr, "nbrId: %lu not in memory\n", nbrId);
#endif

  if(nbrId >= vertexOffset && nbrId < vertexOffset + vertexListSize) {
    *neighbor = vertexList[nbrId-vertexOffset];
    return;
  }

#ifdef USE_THREAD_CACHE
  typename std::map<IdType, VertexType>::iterator it = vertexMap->find(nbrId);
  if(it != vertexMap->end()) {
    *neighbor = it->second;
    return;
  }
#endif
  assert(vertexIO->file_get(nbrId, *neighbor) == 0);
#ifdef USE_THREAD_CACHE
  HugeMem<VertexType, EdgeListType>::insertToVertexMap(nbrId, neighbor);
#endif
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::updateVertex(IdType id) {
  vertexListDirty = true;
#ifdef PART_MEMORY
  memoryVertices[id] = vertexList[id];
#endif

#ifdef USE_THREAD_CACHE
  typename std::map<IdType, VertexType>::iterator it = vertexMap->find(id);
  if(it != vertexMap->end())
    it->second = vertexList[id];
#endif
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::updateEdgeList(IdType id) {
  edgeListDirty = true;
}

//------------------------------------------------------------------
// computation threads
template<typename VertexType, typename EdgeListType>
void *HugeMem<VertexType, EdgeListType>::Worker(void *argPtr) {
  local_id = *(int *)argPtr;
  vertexList = new VertexType[vertexListSize];
  edgeLists = new EdgeListType[edgeListsSize];
#ifdef USE_THREAD_CACHE
  vertexMap = new std::map<IdType, VertexType>();
#endif

  while (1) {
    pthread_barrier_wait(&bar); 

    if ( LFLAG == TFLAG_CLOSE ) {
      break;
    }
    else if ( LFLAG == TFLAG_WORK ) {
      setLoopIndexRange();			// distribute workload
      thread_func();
      LFLAG = TFLAG_IDLE;
      pthread_barrier_wait(&bar);		// join all computation threads after parallel computations
    }
    else 
      assert(0);
  }

  delete[] vertexList;
  delete[] edgeLists;
#ifdef USE_THREAD_CACHE
  delete vertexMap;
#endif
  return NULL;
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::startWorkTimes() {
  gettimeofday(&workStartTime, NULL);
  vertexIO->clearStats();
  edgeListIO->clearStats();
  reInitializeActivations();
#ifdef PART_MEMORY
  assert(vertexIO->file_get_batch(0, PART_MEMORY, memoryVertices) == FILEIO_SUCCESS);
#endif
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::waitForThread(int id) {
  while( thread_flags[id] == TFLAG_WORK )
    time(0);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::waitForAllThreads() {
  assert( local_id == 0 );
  for (int i=1; i<compthreads; i++)
    waitForThread(i);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::startAllThreads(unsigned long lb, unsigned long le, void (*func)()) {
  loop_begin = lb;
  loop_end = le;
  thread_func = func;
  LFLAG = TFLAG_WORK;

  for(int i=1; i<compthreads; i++) {
    thread_flags[i] = TFLAG_WORK;
  }

#ifdef USE_ACTIVATIONS
  std::memset(activated[1], false, workload * sizeof(bool)); // will be marked by udpateVertex() in current iteration
  for(int i=0; i<compthreads; i++) nactive[i] = 0;
#endif // USE_ACTIVATIONS

  pthread_barrier_wait(&bar);		// wake up all computation threads 

  setLoopIndexRange();			// distribute workload
  thread_func();
  LFLAG = TFLAG_IDLE;

  pthread_barrier_wait(&bar);		// join all computation threads after parallel computations

#ifdef USE_ACTIVATIONS
  bool *tmp = activated[0];
  activated[0] = activated[1];
  activated[1] = tmp;
#endif // USE_ACTIVATIONS

  gettimeofday(&workEndTime, NULL);
  fprintf(stderr, "Work elapsed time: %.3lf\n", timevalToDouble(workEndTime)-timevalToDouble(workStartTime));
  if(signalled) Shutdown(), exit(0);
}

//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::closeAllThreads() {
#ifndef GRAPHCHI
  waitForAllThreads();
  for(int i=1; i<compthreads; i++)
    thread_flags[i] = TFLAG_CLOSE;

  pthread_barrier_wait(&bar); 		// wake up all computation threads and let them terminate
#endif
}

#ifdef USE_METIS
//------------------------------------------------------------------
//     [0, bounds[0]), [bounds[0], bounds[1]), ...
// === [0, bounds[0]-1], [bounds[0], bounds[1]-1], ...
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::setLoopIndexRange() {
  local_loop_begin = local_loop_index = (local_id == 0) ? 0 : bounds[local_id-1];
  local_loop_end = (local_id == compthreads-1) ? workload : bounds[local_id]-1;
}
#else
//------------------------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::setLoopIndexRange() {
  local_loop_begin = local_loop_index = workloadPerThread * local_id;
  local_loop_end = local_loop_begin + workloadPerThread;
  if (local_id == compthreads - 1)	// if last thread
    local_loop_end = workload;
}
#endif

//------------------------------------------------------------------
// false if loop ends
template<typename VertexType, typename EdgeListType>
bool HugeMem<VertexType, EdgeListType>::getLoopIndex(IdType* index) {
#ifdef USE_ACTIVATIONS
  while(local_loop_index < local_loop_end && !activated[0][local_loop_index]) local_loop_index++;
#endif // USE_ACTIVATIONS
  if (local_loop_index < local_loop_end) {
    *index = local_loop_index;
    ++local_loop_index;
    return true;
  }
  return false;
}

#ifdef USE_ACTIVATIONS
//--------------------------------------------------
template<typename VertexType, typename EdgeListType>
void HugeMem<VertexType, EdgeListType>::activateNeighborsForNextIteration(EdgeListType *edgeList) {
  for(int i=0; i<edgeList->nbrs_size(); i++)
    activated[1][edgeList->nbrs(i)] = true;
  nactive[local_id] += edgeList->nbrs_size();
}
#endif // USE_ACTIVATIONS

//--------------------------------------------------
template<typename VertexType, typename EdgeListType>
inline void HugeMem<VertexType, EdgeListType>::reInitializeActivations() {
#ifdef USE_ACTIVATIONS
    std::memset(activated[0], true, workload * sizeof(bool));  // compute everything for 0-th iteration
    std::memset(activated[1], false, workload * sizeof(bool));  // will be marked in udpateVertex()
#endif // USE_ACTIVATIONS
  }

#endif // _HM_HPP_
