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
void* doCoarsen(void* arg)
{
  double time_coarsen = 0.0;
  unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
  GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
  Partitioner& partitioner = mr->partitioner;
  //std::cout << "DoCoarsen tid, *mr:" << tid << "\n";  //GK

    std::ifstream infile(mr->inputFileName.c_str()); 
    assert(infile.is_open());
    fprintf(stderr,"\n Input file: %s\n",mr->inputFileName.c_str());
    infile.seekg(tid*mr->bytesPerFile);

    std::string line;
    unsigned nvts=0;
 //   edgeLists = new RecordType[mr->graph.v];

    fprintf(stderr,"\n Input line: %s\n",line.c_str());
    fprintf(stderr, "Reading edge lists\n");
    while(std::getline(infile, line)) {
         time_coarsen -= getTimer();
          mr->coarsenG(tid, line);
          time_coarsen += getTimer();
          ++nvts;
    }
  
//  fprintf(stderr, "Written to disk: %s \n", partitioner.getWrittenToDisk() );
//  fprintf(stderr, "thread %u waiting for others to finish work\n", tid);
  time_coarsen -= getTimer();
//  pthread_barrier_wait(&(mr->barCoarsen));
  if(partitioner.getWrittenToDisk())
    mr->partitioner.flushBResidues(tid);

  time_coarsen += getTimer();

  mr->coarsen_times[tid] = time_coarsen;
 
     fprintf(stderr,"\nSuccessfully Read Graph from file\n");
  return NULL;
}


//---------------------------------------------
// Refiner driver
/*void* doRefine(void* arg)
{


   return NULL;
}

//---------------------------------------------
// In Memory Refiner driver
//---------------------------------------------
void* doInMemoryRefine(void* arg) {
  
   
   return NULL;
}
*/
//=============================================
// Member Functions
//+++++++++++++++++++++++++++++++++++++++++++++
void GraphParts::run()
{
  fprintf(stderr, "initializing");

  fprintf(stderr, "Init Graph partitioners in-Memory Buffers\n");
  double init_time = -getTimer();
  partitioner.initg(nCoarseners, nRefiners, batchSize, kBItems, nparts); // GK 

  fprintf(stderr, "Partitioning input for Parallel Reads\n");
  partitionInputForParallelReads();

  coarsen_times.resize(nCoarseners, 0.0);
  init_time += getTimer();

  fprintf(stderr, "Reading Graph from file\n");
  parallelExecute(doCoarsen, this, nCoarseners);

//  fprintf(stderr, "Running Coarseners\n");
//  parallelExecute(doCoarsen, this, nCoarseners);

  fprintf(stderr, "Graph partitioned. Shutting down.\n");

  fprintf(stderr, "--------------------------------------\n");

  //parts.shutdown();

  std::cout << "------- Final Time ---------" << std::endl;
  std::cout << " Init time : " << init_time << " (msec)" << std::endl;

  auto coarsen_time = max_element(std::begin(coarsen_times), std::end(coarsen_times)); 
  std::cout << " Coarsen time : " << *coarsen_time << " (msec)" << std::endl;

  std::cout << std::endl;
}

//--------------------------------------------
void GraphParts::partitionInputForParallelReads() {
  // Get size of input file
  std::ifstream in(inputFileName.c_str(), std::ifstream::ate | std::ifstream::binary); assert(in.is_open());
  size_t fileSizeInBytes = in.tellg();
  efprintf(stderr, "fileSizeInBytes: %zu\n", fileSizeInBytes);
  in.close();

  bytesPerFile = fileSizeInBytes/nCoarseners;
}

//--------------------------------------------
void GraphParts::init(const std::string input, const unsigned coarseners, const unsigned refiners, const unsigned bSize, const unsigned kItems, const unsigned nParts) {

  inputFileName = input;
  std::cout << "Input file name: " << inputFileName << std::endl;
  
  nCoarseners = coarseners;
  nRefiners = refiners;
  batchSize = bSize;
  kBItems = kItems;
 // nVertices = v;
  nparts = nParts;
 
  
  //setRefiners(std::min(nCoarseners, nRefiners));

  std::cout << "nCoarseners: " << nCoarseners << std::endl;
  std::cout << "nRefiners: " << nRefiners << std::endl;
  std::cout << "batchSize: " << batchSize << std::endl;
  std::cout << "topk: " << kBItems << std::endl;

//  pthread_barrier_init(&barCoarsen, NULL, nCoarseners);
//  pthread_barrier_init(&barRefine, NULL, nRefiners);
}

//--------------------------------------------  GK
//void GraphParts::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const unsigned int* mapSupRowstoRows)
void GraphParts::writeBuf(const unsigned tid, const unsigned to, const std::vector<unsigned>& from){
   partitioner.writeBuf(tid, to, from);
}

/*void GraphParts::coarsen(const unsigned tid, graph_t *graph, graph_t *cgraph, unsigned CoarsenTo)
{
  partitioner.coarsen(tid, graph, cgraph, CoarsenTo);
//  partitioner.coarsen(tid, cgraph, CoarsenTo, numEdgesSupRowstoRows, mapSupRowstoRows);
}
*/
