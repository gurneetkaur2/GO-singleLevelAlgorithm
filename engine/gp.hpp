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

#include <sys/mman.h>

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
	infile.seekg(std::ios::beg); //need this for Tid = 0
        unsigned lineId = tid*mr->linesPerThread + tid;//0, 4, 8
	//infile.seekg(tid*lineId, infile.beg);
       // mr->end_read[tid] = (tid+1)*mr->linesPerThread + 1; //4, 7, 9
        mr->end_read[tid] = (tid+1)*mr->linesPerThread + tid; //4, 7, 9
	fprintf(stderr, "\nCreating memory partitions nVertices: %d, nEdges: %d\n", mr->nVertices, mr->nEdges);
        
	std::string line;
 
       if(tid > 0){
        //   infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	   infile.seekg(std::ios::beg);
    	    for(int i=0; i < lineId; ++i){  //remove -1 if starting from 0
        	infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
        }
	while(std::getline(infile, line, '\n')){
                //fprintf(stderr,"\nTID: %d, lineID %d THREADCT %d\n", tid, lineId, threadCt);
        //        fprintf(stderr,"\n ********** TID: %d, lineID %d, end_read: %d \n", tid, lineId+1, mr->end_read[tid]);
           if(lineId <= mr->nVertices && lineId <= mr->end_read[tid]){
		time_mparts -= getTimer();
		mr->createMParts(tid, line, ++lineId, mr->hDegree);     
	        time_mparts += getTimer();
      //         fprintf(stderr,"\nTID: %d THreadCT %d ", tid, threadCt);
	  }
          else
            break;
        }

	//  fprintf(stderr, "Written to disk: %s \n", partitioner.getWrittenToDisk() );
	  fprintf(stderr, "thread %u waiting for others to finish work\n", tid);
	//copy the local partition to global 

	time_mparts -= getTimer();
	pthread_barrier_wait(&(mr->barMParts));

	//  mr->partitioner.gCopy(tid);

	if(partitioner.getWrittenToDisk())
		mr->partitioner.flushBResidues(tid);


	if (tid == 0) {
		mr->partitioner.gCopy(tid);
	}
	time_mparts += getTimer();

	mr->mparts_times[tid] = time_mparts;
	// fprintf(stderr,"\nAfter flushign residues");

	return NULL;
}


//---------------------------------------------
// Refine driver
void* doRefine(void* arg)
{

	double time_refine = -getTimer();

	unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
	GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
	Partitioner& partitioner = mr->partitioner;

	//  fprintf(stderr, "\nDoRefine: Initializing read parameter\n");
	mr->readInit(tid);
       // for(unsigned i = 0; i < mr->nParts; i++){
	for(auto it = partitioner.fetchPIds[tid].begin(); it != partitioner.fetchPIds[tid].end(); ++it) {
	    unsigned hipart = tid;
	    unsigned whereMax = *it; 
	    if(whereMax == tid){
 	       partitioner.pIdsCompleted[tid][*it] = true;
               continue;
	     }
     //   partitioner.fetchPIds.insert(std::pair<std::pair<unsigned, unsigned>, bool>(std::make_pair(hipart, whereMax), 0);
        bool ret = partitioner.checkPIDStarted(tid, hipart, whereMax);
	time_refine += getTimer();

        fprintf(stderr,"\n-- COMBINE TID: %d, hipart: %d whereMax: %d\n", tid, hipart, whereMax);
	time_refine = -getTimer();
       // fprintf(stderr,"\n PIDS: %d ", partitioner.pIdsCompleted[hipart][whereMax]);
            partitioner.bRefine(tid, hipart, whereMax, ret);
 	    partitioner.pIdsCompleted[hipart][whereMax] = true;

           // fprintf(stderr,"\n WhereMax: %d, hipart: %d, PIDS: %d ", whereMax, hipart, partitioner.pIdsCompleted[tid][*it]);
  }
	pthread_barrier_wait(&(mr->barRefine));
	  mr->afterRefine(tid, mr->nVertices);
	time_refine += getTimer();
	mr->refine_times[tid] += time_refine;
       
  	fprintf(stderr, "thread %u finished Combining\n", tid);
	return NULL;

}

//---------------------------------------------
// In Memory Refiner driver
//---------------------------------------------
void* doInMemoryRefine(void* arg) {

	double time_refine = -getTimer();

	  unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
	  GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
	  Partitioner& partitioner = mr->partitioner;

	  partitioner.initiateInMemoryRefine(tid); 
       partitioner.ComputeBECut(tid, partitioner.readBufMap[tid]);
	time_refine += getTimer();
	pthread_barrier_wait(&(mr->barRefine));
       if(tid == 0)
       fprintf(stderr,"\n Total EdgeCuts Before: %d\n", partitioner.countTotalPECut(tid));
	time_refine = -getTimer();
	  mr->refineInit(tid);
 
	for(auto it = partitioner.fetchPIds[tid].begin(); it != partitioner.fetchPIds[tid].end(); ++it) {
	    unsigned hipart = tid;
	    unsigned whereMax = *it; 
	    if(whereMax == tid){
 	       partitioner.pIdsCompleted[tid][*it] = true;
               continue;
	     }
     //   partitioner.fetchPIds.insert(std::pair<std::pair<unsigned, unsigned>, bool>(std::make_pair(hipart, whereMax), 0);
        bool ret = partitioner.checkPIDStarted(tid, hipart, whereMax);
	time_refine += getTimer();

        fprintf(stderr,"\n-- InMemory TID: %d, hipart: %d whereMax: %d\n", tid, hipart, whereMax);
	time_refine = -getTimer();
       // fprintf(stderr,"\n PIDS: %d ", partitioner.pIdsCompleted[hipart][whereMax]);
            partitioner.inMemRefine(tid, hipart, whereMax, ret);
 	    partitioner.pIdsCompleted[hipart][whereMax] = true;

  }
        
	pthread_barrier_wait(&(mr->barRefine));
	mr->afterRefine(tid, mr->nVertices);
	time_refine += getTimer();
	mr->refine_times[tid] += time_refine;
   //    partitioner.ComputeBECut(tid, partitioner.readBufMap[tid]);
	pthread_barrier_wait(&(mr->barRefine));
       if(tid == 0)
        fprintf(stderr,"\n Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
  	fprintf(stderr, "thread %u finished InMemory Refining\n", tid);
	 
	return NULL;
}

//=============================================
// Member Functions
//+++++++++++++++++++++++++++++++++++++++++++++
void GraphParts::run()
{
	fprintf(stderr, "initializing\n");
	fprintf(stderr, "Init Graph partitioners in-Memory Buffers\n");

        double init_time = -getTimer();
	partitioner.initg(nVertices, hDegree, nThreads, batchSize, kBItems, nParts); // GK 


	fprintf(stderr, "Partitioning input for Parallel Reads\n");
	partitionInputForParallelReads();

	end_read.resize(nThreads, 0.0);
	mparts_times.resize(nThreads, 0.0);
	refine_times.resize(nParts, 0.0);
        writeBuf_times.resize(nThreads, 0.0);
        flushResidues_times.resize(nThreads, 0.0);
        infinimem_read_times.resize(nParts, 0.0);
        infinimem_write_times.resize(nThreads, 0.0);
        localCombinedPairs.resize(nThreads, uint64_t(0));
 
 	init_time += getTimer();

	fprintf(stderr, "Reading Graph from file\n");
	parallelExecute(doMParts, this, nThreads);
	fprintf(stderr,"\nSuccessfully Uploaded the Graph\n");

	//  fprintf(stderr, "Running Coarseners\n");
	//  parallelExecute(doCoarsen, this, nCoarseners);
	//      fprintf(stderr,"\nGP.HPP before Running refiners");

	 if(!partitioner.getWrittenToDisk()) {
	    fprintf(stderr, "Running InMemoryRefiners\n");
	    parallelExecute(doInMemoryRefine, this, nParts);
	    partitioner.releaseInMemStructures();
	    } else {
	      partitioner.releaseInMemStructures();
	fprintf(stderr, "\nRunning Combiners\n");
	parallelExecute(doRefine, this, nParts);

	partitioner.releaseReadPartStructures();
//	fprintf(stderr, "\nRunning Refiners\n");
//	parallelExecute(doRefine, this, nParts);
	 }

	fprintf(stderr, "Graph partitioned. Shutting down.\n");

	fprintf(stderr, "--------------------------------------\n");

	partitioner.shutdown();

	std::cout << "------- Final Time ---------" << std::endl;
	std::cout << " Init time : " << init_time << " (msec)" << std::endl;


	auto mparts_time = max_element(std::begin(mparts_times), std::end(mparts_times)); 
	std::cout << " Initial Partitioning time : " << *mparts_time/1000 << " (sec)" << std::endl;

	auto refine_time = max_element(std::begin(refine_times), std::end(refine_times));
	std::cout << " Refine time : " << *refine_time << " msec" << std::endl;

        std::cout<< " Partition+Refine time: " << ((*mparts_time) + (*refine_time)) << " msec" <<std::endl;
        std::cout<< " Init+Partition+Refine time: " << ((init_time) + (*mparts_time) + (*refine_time)) << " msec" <<std::endl;

//	auto refine_time = max_element(std::begin(refine_times), std::end(refine_times));
//	std::cout << " Refine time : " << *refine_time << " msec" << std::endl;

	uint64_t combinedPairs = std::accumulate(std::begin(localCombinedPairs), std::end(localCombinedPairs), uint64_t(0));
  std::cout << " Total combined pairs : " << combinedPairs << std::endl;

        auto writeBuf_time = max_element(std::begin(writeBuf_times), std::end(writeBuf_times)); 
        std::cout << " writeBuf time : " << *writeBuf_time << " (msec)" << std::endl;

        auto flushResidues_time = max_element(std::begin(flushResidues_times), std::end(flushResidues_times));
        std::cout << " flushResidues time : " << *flushResidues_time << " (msec)" << std::endl;

        auto infinimem_read_time = max_element(std::begin(infinimem_read_times), std::end(infinimem_read_times));
  	std::cout << " InfiniMem Read time: " << *infinimem_read_time << " (msec)" << std::endl;

	auto infinimem_write_time = max_element(std::begin(infinimem_write_times), std::end(infinimem_write_times));
  	std::cout << " InfiniMem Write time: " << *infinimem_write_time << " (msec)" << std::endl;

	std::cout << std::endl;
}

//--------------------------------------------
void GraphParts::partitionInputForParallelReads() {
	// Get size of input file
	std::ifstream in(inputFileName.c_str(), std::ifstream::ate | std::ifstream::binary); assert(in.is_open());
	size_t fileSizeInBytes = in.tellg();
	fprintf(stderr, "fileSizeInBytes: %zu\n", fileSizeInBytes);
	//fprintf(stderr, "\nnumlines: %u\n", numLines);
	in.close();
	fprintf(stderr,"\n NumLines: %zu", numLines);
	linesPerThread = numLines/nThreads;
	fprintf(stderr,"\nLinesPerThread: %zu \n", linesPerThread);
	bytesPerFile = fileSizeInBytes/nThreads + 1; //fileSizeInBytes/nThreads + 1;
}

//--------------------------------------------
void GraphParts::init(const std::string input, const unsigned nvertices, const unsigned nedges, const unsigned hdegree, const unsigned nthreads, const unsigned nparts, const unsigned bSize, const unsigned kItems) {

	inputFileName = input;
	std::cout << "Input file name: " << inputFileName << std::endl;
	numLines = getNumLines(inputFileName);
	//  fprintf(stderr,"\nINIT NumLines: %zu\n", numLines);

	//nInMemParts and nParts are same
	nVertices = nvertices;
	nEdges = nedges;
 	hDegree = hdegree;
	nThreads = nthreads;
	nParts = nparts;
	batchSize = bSize;
	kBItems = kItems;


	//  setRefiners(std::min(nThreads, nRefiners));
	//TODO:need to check if I need this --  nParts = std::min(nThreads, nParts);

	std::cout << "nVertices: " << nVertices << std::endl;
	std::cout << "nEdges: " << nEdges << std::endl;
	std::cout << "hDegree: " << hDegree << std::endl;
	std::cout << "nThreads: " << nThreads << std::endl;
	std::cout << "nParts: " << nParts << std::endl;
	std::cout << "batchSize: " << batchSize << std::endl;
	std::cout << "topk: " << kBItems << std::endl;

	pthread_barrier_init(&barMParts, NULL, nThreads);
	pthread_barrier_init(&barRead, NULL, nThreads);
	pthread_barrier_init(&barRefine, NULL, nThreads);
}

//--------------------------------------------  GK
void GraphParts::writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hidegree = 0){
	partitioner.writeBuf(tid, to, from, hidegree);
}

//--------------------------------------------
bool GraphParts::read(const unsigned tid) {
	return partitioner.read(tid);
}

//--------------------------------------------
void GraphParts::printParts(const unsigned tid, std::string outputPrefix) {
	return partitioner.printParts(tid, outputPrefix);
}

//--------------------------------------------
void GraphParts::ComputeBECut(const unsigned tid, const InMemoryContainer& inMemMap) {
	return partitioner.ComputeBECut(tid, inMemMap);
}

//--------------------------------------------
void GraphParts::cWrite(const unsigned tid, unsigned noItems, InMemoryConstIterator end) {
	return partitioner.cWrite(tid, noItems, end);
}

//--------------------------------------------
unsigned GraphParts::countTotalPECut(const unsigned tid) {
	return partitioner.countTotalPECut(tid);
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
void GraphParts::refineInit(const unsigned tid) {
	return partitioner.refineInit(tid);
}

//--------------------------------------------
void GraphParts::subtractRefineTimes(const unsigned tid, double stime) {
	refine_times[tid] -= stime;
}
