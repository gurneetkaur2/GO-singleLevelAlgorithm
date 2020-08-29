#include "gp.h"
#include "partitioner.hpp"
//#include "../programs/graph.h"
#include <utility> // for std::pair
#include <fstream>
#include <iostream> //GK
#include <numeric>
#include <pthread.h>
#include <atomic>
#include <vector>
#include <string>
#include <ctime>
#include <math.h>

#include <sys/mman.h>

  static int commentLines = 0; 
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
//	fprintf(stderr, "\n DoMparts tid %d  ", tid);  //GK
   
	//  mr->writeInit(rand() % mr->nParts);
//	fprintf(stderr, "\n After writeinit tid %d  ", tid);  //GK

	std::ifstream infile(mr->inputFileName.c_str()); 
	assert(infile.is_open());
	  //      time_mparts += getTimer();
//	fprintf(stderr,"\n Input file: %s\n",mr->inputFileName.c_str());
	//	time_mparts -= getTimer();
	infile.seekg(std::ios::beg); //need this for Tid = 0
        unsigned lineId = tid*mr->linesPerThread + tid;//0, 4, 8
	//infile.seekg(tid*lineId, infile.beg);
       // mr->end_read[tid] = (tid+1)*mr->linesPerThread + 1; //4, 7, 9
        mr->end_read[tid] = (tid+1)*mr->linesPerThread + tid; //4, 7, 9
	  //      time_mparts += getTimer();
	//fprintf(stderr, "\nCreating memory partitions nVertices: %d, Partitions: %d\n", mr->nVertices, mr->nParts);
		//time_mparts -= getTimer();
        
	std::string line;
       if(tid > 0){
        //   infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	   infile.seekg(std::ios::beg);
    	    for(int i=0; i < lineId; ++i){  //remove -1 if starting from 0
        	infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
        }
   //     fprintf(stderr,"\n Type : %s \n", mr->inType.c_str());
	while(std::getline(infile, line, '\n')){
                //fprintf(stderr,"\nTID: %d, lineID %d THREADCT %d\n", tid, lineId, threadCt);
      //          fprintf(stderr,"\n ********** TID: %d, lineID %d, end_read: %d Type: %s \n", tid, lineId+1, mr->end_read[tid], mr->inType.c_str());
           if( lineId < mr->nVertices && lineId <= mr->end_read[tid]){

		time_mparts -= getTimer();
		mr->createMParts(tid, line, mr->inType, ++lineId, mr->hDegree);     
	        time_mparts += getTimer();
      //         fprintf(stderr,"\nTID: %d returned ", tid);
	  }
    /*      else if(mr->inType == "edgelist" && lineId <= mr->end_read[tid]){

		time_mparts -= getTimer();
		mr->createMParts(tid, line, mr->inType, ++lineId, mr->hDegree);     
	        time_mparts += getTimer();
      //         fprintf(stderr,"\nTID: %d THreadCT %d ", tid, threadCt);
	  }*/
          else
            break;
        }

	//  fprintf(stderr, "Written to disk: %s \n", partitioner.getWrittenToDisk() );
	  fprintf(stderr, "thread %u waiting for others to finish work lineId %d \n", tid, lineId);
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
        int k = 0;
       // for(unsigned i = 0; i < mr->nParts; i++){
//	for(auto it = partitioner.fetchPIds[tid].begin(); it != partitioner.fetchPIds[tid].end(); ++it) {
	for(auto it = partitioner.fetchPIds[tid].begin(); it != partitioner.fetchPIds[tid].end(); ++it) {
	    unsigned hipart = tid;
	    unsigned whereMax = *it; 
	    if(whereMax == tid){
 	       partitioner.pIdsCompleted[tid][*it] = true;
               continue;
	     }
       else if (hipart < mr->nParts && whereMax >= mr->nParts){
          partitioner.pIdsCompleted[tid][*it] = true;
          continue;
       }
       else if ( hipart >= mr->nParts && whereMax < mr->nParts) {
         partitioner.pIdsCompleted[tid][*it] = true;
         continue;
       }
        bool ret = partitioner.checkPIDStarted(tid, hipart, whereMax);
	//time_refine += getTimer();

    //    fprintf(stderr,"\n-- COMBINE TID: %d, hipart: %d whereMax: %d, ret: %d \n", tid, hipart, whereMax, ret);
//	time_refine = -getTimer();
       // fprintf(stderr,"\n PIDS: %d ", partitioner.pIdsCompleted[hipart][whereMax]);
            partitioner.bRefine(tid, hipart, whereMax, ret);
            if(ret == true) {
        // fprintf(stderr,"\nTID %d going to write part ", tid);
               partitioner.writePartInfo(tid, hipart, whereMax);
          // pthread_mutex_unlock(&locks[tid]);
            }
	pthread_barrier_wait(&(mr->barWriteInfo));
            partitioner.clearMemorystructures(tid);
 	    partitioner.pIdsCompleted[hipart][whereMax] = true;
   
       //TODO:: may be wait before starting next iter?
//	pthread_barrier_wait(&(mr->barRefine));
	//time_refine += getTimer();
    //    fprintf(stderr,"\nTID %d finished iter %d  using disk", tid,++k);
	//time_refine = -getTimer();
      // if(tid ==0) partitioner.pIdStarted.clear();
  }
	time_refine += getTimer();
	fprintf(stderr, "\nthread %u waiting for others to finish Refine\n", tid);
	time_refine = -getTimer();
	pthread_barrier_wait(&(mr->barAfterRefine));
//  if(tid == 0)
	double time_aftr_refine = -getTimer();
	  mr->afterRefine(tid, mr->nVertices);
	time_aftr_refine += getTimer();
	mr->aftr_refine_times[tid] += time_aftr_refine;

	time_refine += getTimer();
	mr->refine_times[tid] += time_refine;

//	pthread_barrier_wait(&(mr->barRefine));
       /*if(tid == 0){
        partitioner.setTotalCuts(tid);
        fprintf(stderr,"\n\n Total EdgeCuts BEFORE: %d\n", partitioner.countTotalPECut(tid));
	}*/
	pthread_barrier_wait(&(mr->barClear));
    if(tid == 0)
		mr->partitioner.gCopy(tid);
        partitioner.readClear(tid); 
	partitioner.refineInit(tid);
        partitioner.cread(tid);
	pthread_barrier_wait(&(mr->barRefine));
       if(tid == 0){
        partitioner.setTotalCuts(tid);
//	time_refine += getTimer();
       fprintf(stderr,"\n Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
//	time_refine = -getTimer();
	}
       
//	time_refine += getTimer();
  	//fprintf(stderr, "\nthread %u finished Refining\n", tid);
	//time_refine = -getTimer();
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
	pthread_barrier_wait(&(mr->barRefine));
        if(tid == 0)  //Need this condition when used in multi thread envt.
	    partitioner.releaseInMemStructures();
	//time_refine += getTimer();
//	  mr->refineInit(tid);
        partitioner.setTotalCuts(tid);
	//time_refine = -getTimer();
	  mr->refineInit(tid);
          int k = 0; 
	for(auto it = partitioner.fetchPIds[tid].begin(); it != partitioner.fetchPIds[tid].end(); ++it) {
	    unsigned hipart = tid;
	    unsigned whereMax = *it; 
	    if(whereMax == tid){
 	       partitioner.pIdsCompleted[tid][*it] = true;
               continue;
	     }
       else if (hipart < mr->nParts && whereMax >= mr->nParts){
          partitioner.pIdsCompleted[tid][*it] = true;
          continue;
       }
       else if ( hipart >= mr->nParts && whereMax < mr->nParts) {
         partitioner.pIdsCompleted[tid][*it] = true;
         continue;
       }
        bool ret = partitioner.checkPIDStarted(tid, hipart, whereMax);
	//time_refine += getTimer();

    //    fprintf(stderr,"\n-- InMemory TID: %d, hipart: %d whereMax: %d, ret: %d \n", tid, hipart, whereMax, ret);
	//time_refine = -getTimer();
       // fprintf(stderr,"\n PIDS: %d ", partitioner.pIdsCompleted[hipart][whereMax]);
            partitioner.inMemRefine(tid, hipart, whereMax, ret);
        // fprintf(stderr,"\nTID %d back after inMem Refine ret: %d ", tid, ret);
          //  pthread_barrier_wait(&(mr->barRefine));
            if(ret == true) {
//         fprintf(stderr,"\nTID %d going to write part ", tid);
               partitioner.writePartInfo(tid, hipart, whereMax);
          // pthread_mutex_unlock(&locks[tid]);
            }
	pthread_barrier_wait(&(mr->barWriteInfo));
            partitioner.clearMemorystructures(tid);
 	    partitioner.pIdsCompleted[hipart][whereMax] = true;
//	time_refine += getTimer();
  //      fprintf(stderr,"\nTID %d finished iter %d  in Memory", tid,++k);
	//time_refine = -getTimer();
       //if(tid ==0) partitioner.pIdStarted.clear();

  }
        
	time_refine += getTimer();
	fprintf(stderr, "\nthread %u waiting for others to finish InMemory Refine\n", tid);
	time_refine = -getTimer();
	pthread_barrier_wait(&(mr->barAfterRefine));
	double time_aftr_refine = -getTimer();
//  if(tid == 0)
  	mr->afterRefine(tid, mr->nVertices);

	time_aftr_refine += getTimer();
	mr->aftr_refine_times[tid] += time_aftr_refine;

	time_refine += getTimer();
	mr->refine_times[tid] += time_refine;
	pthread_barrier_wait(&(mr->barClear));
    if(tid == 0)
		   mr->partitioner.gCopy(tid);
       //   partitioner.setTotalCuts(tid);
	pthread_barrier_wait(&(mr->barRefine));
    if(tid == 0)
        fprintf(stderr,"\n\n Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
          mr->refineInit(tid);
    // fprintf(stderr,"\nTID %d RefineMap size: %d Total Cuts: %d", tid, partitioner.refineMap[tid].size(), partitioner.countTotalPECut(tid));
      //  partitioner.ComputeBECut(tid, partitioner.refineMap[tid]);
          partitioner.cread(tid);
	  pthread_barrier_wait(&(mr->barRefine));
       if(tid == 0){
         partitioner.setTotalCuts(tid);
        // time_refine += getTimer();

        fprintf(stderr,"\n\n Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
	//   time_refine = -getTimer();
  	}

	//time_refine += getTimer();
//	fprintf(stderr, "thread %u finished InMemory Refining\n", tid);
//	time_refine = -getTimer();
	 
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
	partitioner.initg(nVertices, hDegree, nThreads, batchSize, kBItems, nParts, nrefiners); // GK 


	fprintf(stderr, "Partitioning input for Parallel Reads\n");
	partitionInputForParallelReads();

	end_read.resize(nThreads, 0);
	mparts_times.resize(nThreads, 0.0);
// if(nParts < 10){
	refine_times.resize(nrefiners, 0.0);
	aftr_refine_times.resize(nrefiners, 0.0);
        writeBuf_times.resize(nThreads, 0.0);
        flushResidues_times.resize(nThreads, 0.0);
        infinimem_read_times.resize(nrefiners, 0.0);
        infinimem_write_times.resize(nThreads, 0.0);
        infinimem_cread_times.resize(nrefiners, 0.0);
        infinimem_cwrite_times.resize(nrefiners, 0.0);
/*    }
    else{
        refine_times.resize(nParts, 0.0);
        writeBuf_times.resize(nThreads, 0.0);
        flushResidues_times.resize(nThreads, 0.0);
        infinimem_read_times.resize(nParts, 0.0);
        infinimem_write_times.resize(nThreads, 0.0);
        infinimem_cread_times.resize(nParts, 0.0);
        infinimem_cwrite_times.resize(nParts, 0.0);
   }*/
        localCombinedPairs.resize(nThreads, uint64_t(0));

 	init_time += getTimer();

	fprintf(stderr, "Reading Graph from file\n");
  partitioner.writeInit();
	parallelExecute(doMParts, this, nThreads);
	fprintf(stderr,"\nSuccessfully Uploaded the Graph\n");

	//  fprintf(stderr, "Running Coarseners\n");
	//  parallelExecute(doCoarsen, this, nCoarseners);
	//      fprintf(stderr,"\nGP.HPP before Running refiners");

	 if(!partitioner.getWrittenToDisk()) {
	    fprintf(stderr, "Running InMemoryRefiners\n");
	    parallelExecute(doInMemoryRefine, this, nrefiners);
//	    partitioner.releaseInMemStructures();
	    } else {
	      partitioner.releaseInMemStructures();
	fprintf(stderr, "\nRunning Combiners\n");
	parallelExecute(doRefine, this, nrefiners);

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

	auto aftr_refine_time = max_element(std::begin(aftr_refine_times), std::end(aftr_refine_times));
	std::cout << " Writing to output time : " << *aftr_refine_time << " msec" << std::endl;

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
       
        auto infinimem_cread_time = max_element(std::begin(infinimem_cread_times), std::end(infinimem_cread_times));
  	std::cout << " InfiniMem Sequential Read time: " << *infinimem_cread_time << " (msec)" << std::endl;

	auto infinimem_cwrite_time = max_element(std::begin(infinimem_cwrite_times), std::end(infinimem_cwrite_times));
  	std::cout << " InfiniMem Sequential Write time: " << *infinimem_cwrite_time << " (msec)" << std::endl;

	std::cout << std::endl;
}

//--------------------------------------------
void GraphParts::partitionInputForParallelReads() {
	// Get size of input file
/*	std::ifstream in(inputFileName.c_str(), std::ifstream::ate | std::ifstream::binary); assert(in.is_open());
	size_t fileSizeInBytes = in.tellg();
	fprintf(stderr, "fileSizeInBytes: %zu\n", fileSizeInBytes);
	//fprintf(stderr, "\nnumlines: %u\n", numLines);
	in.close();
*/	fprintf(stderr,"\n NumLines: %zu", numLines);
	linesPerThread = numLines/nThreads;
	fprintf(stderr,"\nLinesPerThread: %zu \n", linesPerThread);
	//bytesPerFile = fileSizeInBytes/nThreads + 1; //fileSizeInBytes/nThreads + 1;
}

//--------------------------------------------
void GraphParts::init(const std::string input, const std::string type, const unsigned nvertices, const unsigned hdegree, const unsigned nthreads, const unsigned nparts, const unsigned mSize, const unsigned kItems) {

	inputFileName = input;
	std::cout << "Input file name: " << inputFileName << std::endl;
//	numLines = getNumLines(inputFileName);
	//  fprintf(stderr,"\nINIT NumLines: %zu\n", numLines);

	//nInMemParts and nParts are same
	nVertices = nvertices;
//	nEdges = nedges;
 	hDegree = hdegree;
        inType = type;
	nThreads = nthreads;

  unsigned wload = mSize / (nThreads * nparts);
	batchSize = wload + mSize % (nThreads * nparts) ;
	//batchSize = mSize; // / (nThreads * nparts);
	fprintf(stderr, "batch size: %zu\n", batchSize);
//  if(nparts <= 16 && batchSize >=5000){
 //     nrefiners = nparts * 2;
 /*   if(nVertices < 5000000){
    if(nparts <= 2 && bSize < 16000)
      nrefiners = nparts * 2;
    else if (nparts >= 4 && bSize < 16000)
      nrefiners = nparts * 2;
    else
        nrefiners = nparts;
    }
    else{
    //if(nparts <= 2 && bSize < 16000)
    //  nrefiners = nparts * 8;
     if (nparts <= 4 )
      nrefiners = nparts * 2;
    else if (nparts <= 16 )
      nrefiners = nparts * 2;
    }*/
/*  }
   else if(nparts <= 8 && batchSize < 5000){ // for smaller graphs upto 8 parts
    nrefiners = nparts * 2;
  }
  else
  */
  nrefiners = nparts;
	kBItems = kItems;
 
  nParts = nparts;
   
 // if(inType == "adj")
     numLines = nVertices;
               
//  if(inType == "edge")
//     numLines = nEdges;
                                    
                                          
        if (inType == "adjlist" && numLines != nVertices){
		fprintf(stderr, "\nNo. of Vertices %d not correct numlines %d \n", nVertices, numLines);
                assert(false);
	}

/*        if (inType == "edge" && numLines != nEdges){
		fprintf(stderr, "\nNo. of Edges %d not correct\n", numLines);
                assert(false);
	}
*/
	//setRefiners(std::min(nThreads, nParts));
	//TODO:need to check if I need this --  nParts = std::min(nThreads, nParts);

	std::cout << "nVertices: " << nVertices << std::endl;
	//std::cout << "nEdges: " << nEdges << std::endl;
	std::cout << "hDegree: " << hDegree << std::endl;
	std::cout << "nThreads: " << nThreads << std::endl;
	std::cout << "nParts: " << nParts << std::endl;
	std::cout << "batchSize: " << batchSize << std::endl;
	std::cout << "topk: " << kBItems << std::endl;

	pthread_barrier_init(&barMParts, NULL, nThreads);
	pthread_barrier_init(&barRead, NULL, nrefiners);
	pthread_barrier_init(&barRefine, NULL, nrefiners);
	pthread_barrier_init(&barWriteInfo, NULL, nrefiners);
	pthread_barrier_init(&barClear, NULL, nrefiners);
	pthread_barrier_init(&barAfterRefine, NULL, nrefiners);
}

//--------------------------------------------  GK
void GraphParts::writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hidegree = 0){
	partitioner.writeBuf(tid, to, from, hidegree);
}

//--------------------------------------------
void GraphParts::setRefiners(const unsigned nparts) {
       nParts = nparts;             
    
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
void GraphParts::writeInit() {
	return partitioner.writeInit();
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
