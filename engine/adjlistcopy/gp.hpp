#include "gp.h"
#include "partitioner.hpp"
#include "../programs/graph.h"
#include <utility> // for std::pair
#include <fstream>
#include <iostream> //GK
#include <pthread.h>
#include <atomic>
#include <vector>
#include <string>
#include <ctime>
#include <math.h>

//--------------------------------------------
// Helper NON-member Functions
// Facade for BSP style parallel execution for Coarseners and then Refiners
void parallelExecute(void *(*func)(void *), void* arg, unsigned threads)
{
  pthread_t thread_handles[threads];
  std::pair<unsigned, void*> args[threads];

  for (unsigned i = 0; i < threads; i++)
  {
    args[i] = std::make_pair(i, arg);
    pthread_create(&thread_handles[i], NULL, func, &args[i]);
  }

  for (unsigned i = 0; i < threads; i++)
  {
    pthread_join(thread_handles[i], NULL);
  }

}

//---------------------------------------------
// Coarsening driver
void* doMParts(void* arg)
{
  double time_mparts = 0.0;
  unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
  GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
  Partitioner& partitioner = mr->partitioner;
  //std::cout << "DoCoarsen tid, *mr:" << tid << "\n";  //GK

    mr->writeInit(tid);

   std::ifstream infile(mr->inputFileName.c_str()); 
   assert(infile.is_open());
   fprintf(stderr,"\n Input file: %s\n",mr->inputFileName.c_str());
   infile.seekg(tid*mr->bytesPerFile);

   std::string line;
   fprintf(stderr, "Creating memory partitions nVertices: %d, nEdges: %d\n", mr->nVertices, mr->nEdges);
   unsigned k, i;
   std::vector<unsigned> where (mr->nVertices, -1);
   std::vector<unsigned> whereDst (mr->nVertices, -1);    
   unsigned* xadj;
   unsigned* adncy;
   bool firstRep = false;
   xadj = (unsigned *) malloc((mr->nVertices + 1) * sizeof(unsigned));
   adncy = (unsigned *) malloc(((mr->nEdges * 2) + 1) * sizeof(unsigned));

   xadj[0]=0, k=0, i=0;
   while(std::getline(infile, line)) {
        time_mparts -= getTimer();
        std::stringstream inputStream(line);
        unsigned to, from;
        inputStream >> to;
        firstRep = true;

        while(inputStream >> from){
             fprintf(stderr,"\nTO: %zu FROM; %zu ", to, from);
             k = mr->writeBuf(tid, to, from, k, adncy, where, whereDst, firstRep);
             firstRep = false;
        }
        xadj[i+1] = k;
        fprintf(stderr,"xadj: %d, k: %d, i: %d\n", xadj[i+1], k, i);
        i++;
        
          time_mparts += getTimer();
    }
  
//  fprintf(stderr, "Written to disk: %s \n", partitioner.getWrittenToDisk() );
//  fprintf(stderr, "thread %u waiting for others to finish work\n", tid);
  time_mparts -= getTimer();
//  pthread_barrier_wait(&(mr->barCoarsen));
  if(partitioner.getWrittenToDisk())
    mr->partitioner.flushBResidues(tid);

  time_mparts += getTimer();

  mr->mparts_times[tid] = time_mparts;
 
     fprintf(stderr,"\nSuccessfully Uploaded the Graph\n");
  free(xadj); xadj=NULL;
  free(adncy); adncy=NULL;
  return NULL;
}


//---------------------------------------------
// Refiner driver
void* doRefine(void* arg)
{

/*  double time_refine = -getTimer();

  unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
  GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
  Partitioner& partitioner = mr->partitioner;

  mr->beforeRefine(tid);

//  fprintf(stderr, "\nDoRefine: Initializing read parameter\n");
  mr->readInit(tid);

  while(true) {
//  fprintf(stderr, "\nDoRefine: Calling Read\n");
    bool execLoop = mr->read(tid);
    if(execLoop == false) {
      for(InMemoryConstIterator it = partitioner.readBufMap[tid].begin(); it != partitioner.readBufMap[tid].end(); ++it)
        mr->refine(tid, it->first, it->second);
      break;
    }

    unsigned counter = 0; 
    InMemoryContainerIterator it;
    for (it = partitioner.readBufMap[tid].begin(); it != partitioner.readBufMap[tid].end(); ++it) {
      if (counter >= mr->kBItems)
        break;

      mr->refine(tid, it->first, it->second);

      const unsigned rank = it->first;
      auto pos = partitioner.lookUpTable[tid].find(rank);
      assert(pos != partitioner.lookUpTable[tid].end());
      const std::vector<unsigned>& lookVal = pos->second;
      for(unsigned val=0; val<lookVal.size(); val++) {
        partitioner.fetchBatchIds[tid].insert(lookVal[val]);
        partitioner.keysPerBatch[tid][lookVal[val]] += 1; 
      }
      partitioner.lookUpTable[tid].erase(rank);
      counter++;
    }

    partitioner.readBufMap[tid].erase(partitioner.readBufMap[tid].begin(), it);
  }

  mr->afterRefine(tid);

  time_refine += getTimer();
  mr->refine_times[tid] += time_refine;
*/
  return NULL;

}

//---------------------------------------------
// In Memory Refiner driver
//---------------------------------------------
void* doInMemoryRefine(void* arg) {

  /*double time_refine = -getTimer();

  unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
  GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
  Partitioner& partitioner = mr->partitioner;

  mr->beforeRefine(tid);
  InMemoryReductionState state = partitioner.initiateInMemoryReduce(tid); 
  InMemoryContainer record;
  while(partitioner.getNextMinKey(&state, &record)) {
    mr->refine(tid, record.begin()->first, record.begin()->second);
    record.clear();
  }
  mr->afterRefine(tid);
  time_refine += getTimer();
  mr->refine_times[tid] += time_refine;
*/
  return NULL;
}

//=============================================
// Member Functions
//+++++++++++++++++++++++++++++++++++++++++++++
void GraphParts::run()
{
  fprintf(stderr, "initializing\n");

  fprintf(stderr, "Init Graph partitioners in-Memory Buffers\n");
  double run_time = -getTimer();
  partitioner.initg(nVertices, nThreads, batchSize, kBItems, nParts); // GK 


  fprintf(stderr, "Partitioning input for Parallel Reads\n");
  partitionInputForParallelReads();

  mparts_times.resize(nThreads, 0.0);
  refine_times.resize(nParts, 0.0);



  fprintf(stderr, "Reading Graph from file\n");
  parallelExecute(doMParts, this, nThreads);

//  fprintf(stderr, "Running Coarseners\n");
//  parallelExecute(doCoarsen, this, nCoarseners);

  if(!partitioner.getWrittenToDisk()) {
    fprintf(stderr, "Running InMemoryRefiners\n");
    parallelExecute(doInMemoryRefine, this, nParts);
    partitioner.releaseInMemStructures();
  } else {
    partitioner.releaseInMemStructures();
    fprintf(stderr, "Running Refiners\n");
//    parallelExecute(doRefine, this, nRefiners);
  }

  fprintf(stderr, "Graph partitioned. Shutting down.\n");

  fprintf(stderr, "--------------------------------------\n");

  //parts.shutdown();

  std::cout << "------- Final Time ---------" << std::endl;
  std::cout << " Total time : " << run_time << " (msec)" << std::endl;

  run_time += getTimer();

  auto mparts_time = max_element(std::begin(mparts_times), std::end(mparts_times)); 
  std::cout << " Partitioning time : " << *mparts_time/1000 << " (sec)" << std::endl;

  auto refine_time = max_element(std::begin(refine_times), std::end(refine_times));
  std::cout << " Refine time : " << *refine_time << " msec" << std::endl;

  std::cout << std::endl;
}

//--------------------------------------------
void GraphParts::partitionInputForParallelReads() {
  // Get size of input file
  std::ifstream in(inputFileName.c_str(), std::ifstream::ate | std::ifstream::binary); assert(in.is_open());
  size_t fileSizeInBytes = in.tellg();
  efprintf(stderr, "fileSizeInBytes: %zu\n", fileSizeInBytes);
  in.close();

  bytesPerFile = fileSizeInBytes/nThreads;
}

//--------------------------------------------
void GraphParts::init(const std::string input, const unsigned nvertices, const unsigned nedges, const unsigned nthreads, const unsigned nparts, const unsigned bSize, const unsigned kItems) {

  inputFileName = input;
  std::cout << "Input file name: " << inputFileName << std::endl;
  
  //nInMemParts and nParts are same
  nVertices = nvertices;
  nEdges = nedges;
  nThreads = nthreads;
  nParts = nparts;
  batchSize = bSize;
  kBItems = kItems;
 
  
  //setRefiners(std::min(nCoarseners, nRefiners));
//TODO:need to check if I need this --  nParts = std::min(nThreads, nParts);

  std::cout << "nVertices: " << nVertices << std::endl;
  std::cout << "nEdges: " << nEdges << std::endl;
  std::cout << "nThreads: " << nThreads << std::endl;
  std::cout << "nParts: " << nParts << std::endl;
  std::cout << "batchSize: " << batchSize << std::endl;
  std::cout << "topk: " << kBItems << std::endl;

//  pthread_barrier_init(&barMParts, NULL, nInMemParts);
//  pthread_barrier_init(&barReduce, NULL, nRefiners);
}

//--------------------------------------------  GK
//void GraphParts::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const unsigned int* mapSupRowstoRows)
unsigned GraphParts::writeBuf(const unsigned tid, const unsigned to, const unsigned from, unsigned k, unsigned* adjncy, std::vector<unsigned>& where, std::vector<unsigned>& whereDst, bool firstRep){
   k = partitioner.writeBuf(tid, to, from, k, adjncy, where, whereDst, firstRep);
   return k;
}

/*void GraphParts::coarsen(const unsigned tid, graph_t *graph, graph_t *cgraph, unsigned CoarsenTo)
{
  partitioner.coarsen(tid, graph, cgraph, CoarsenTo);
//  partitioner.coarsen(tid, cgraph, CoarsenTo, numEdgesSupRowstoRows, mapSupRowstoRows);
}
*/
bool GraphParts::read(const unsigned tid) {
  return partitioner.read(tid);
}

//--------------------------------------------
void GraphParts::writeInit(const unsigned tid) {
  return partitioner.writeInit(tid);
}

//--------------------------------------------
void GraphParts::readInit(const unsigned tid) {
  return partitioner.readInit(tid);
}

//--------------------------------------------
void GraphParts::subtractRefineTimes(const unsigned tid, double stime) {
  refine_times[tid] -= stime;
}
