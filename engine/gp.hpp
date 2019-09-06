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
        unsigned threadCt = 0; // mr->linesPerThread;
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
                 threadCt++;
      //         fprintf(stderr,"\nTID: %d THreadCT %d ", tid, threadCt);
	  }
          else
            break;
        }

	//  fprintf(stderr, "Written to disk: %s \n", partitioner.getWrittenToDisk() );
	//  fprintf(stderr, "thread %u waiting for others to finish work\n", tid);
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
// Combine driver
void* doCombine(void* arg)
{

	double time_combine = -getTimer();

	unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
	GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
	Partitioner& partitioner = mr->partitioner;

	//  fprintf(stderr, "\nDoRefine: Initializing read parameter\n");
	mr->readInit(tid);

	while(true) {
		//  fprintf(stderr, "\nDoRefine: Calling Read\n");
		bool execLoop = mr->read(tid);
		if(execLoop == false) {
			InMemoryContainerIterator fit;
			for(fit = partitioner.readBufMap[tid].begin(); fit != partitioner.readBufMap[tid].end(); ++fit){
				//        mr->ComputeBECut(tid);
				//        mr->refine(tid, it->first, it->second);
		//		fprintf(stderr,"\ntid: %d, Key: %d", tid, fit->first);
				//partitioner.totalCombined[tid]++;
			}
	//		partitioner.totalCombined[tid] += partitioner.readBufMap[tid].size();
		//	fprintf(stderr,"\n-- tid: %d, readMap Size: %d\n", tid, partitioner.readBufMap[tid].size());
			//Write combined records to a new partition    
			mr->cWrite(tid, partitioner.readBufMap[tid].size(), fit);
			partitioner.readBufMap[tid].erase(partitioner.readBufMap[tid].begin(), fit);
			break;
		}

		unsigned counter = 0; 
		InMemoryContainerIterator it;
		for (it = partitioner.readBufMap[tid].begin(); it != partitioner.readBufMap[tid].end(); ++it) {
			if (counter >= mr->kBItems)
				break;

			//      mr->refine(tid, it->first, it->second);

			const unsigned rank = it->first;
		//	fprintf(stderr,"\n TID: %d, Key: %d", tid, rank);
			auto pos = partitioner.lookUpTable[tid].find(rank);
			assert(pos != partitioner.lookUpTable[tid].end());
			const std::vector<unsigned>& lookVal = pos->second; // find the batch of the vertex
			for(unsigned val=0; val<lookVal.size(); val++) {
				partitioner.fetchBatchIds[tid].insert(lookVal[val]);
				partitioner.keysPerBatch[tid][lookVal[val]] += 1; 
			}
			partitioner.lookUpTable[tid].erase(rank);
			counter++;
		}

	//	partitioner.totalCombined[tid] += mr->kBItems;
	//	fprintf(stderr,"\n-- tid: %d, readMap Size: %d\n", tid, partitioner.readBufMap[tid].size());
		//Write combined records to a new partition    
		mr->cWrite(tid, mr->kBItems, it);

		partitioner.readBufMap[tid].erase(partitioner.readBufMap[tid].begin(), it);
	}

	time_combine += getTimer();
	mr->combine_times[tid] += time_combine;

	return NULL;

}

//---------------------------------------------
// Refiner driver
void* doRefine(void* arg)
{

	double time_refine = -getTimer();

	unsigned tid = static_cast<unsigned>(static_cast<std::pair<unsigned, void*>*>(arg)->first);
	GraphParts *mr = static_cast<GraphParts *>(static_cast<std::pair<unsigned, void*>*>(arg)->second);
	Partitioner& partitioner = mr->partitioner;
	mr->refineInit(tid);
	partitioner.cread(tid);
	pthread_barrier_wait(&(mr->barRefine));
	// Count the total edge cuts and also check the partition with max edgecuts
	if(tid == 0){

		for (unsigned i = 0; i < mr->nThreads; i++){
			partitioner.fetchPIds.insert(i);  //Partition ids  - 0,1,2 .. 
		}
        unsigned tCuts = partitioner.countTotalPECut(tid);
        fprintf(stderr,"\nBefore refining, Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));

		while(partitioner.fetchPIds.size() > 1){
			unsigned hipart = -1;
                        fprintf(stderr,"\nFinding Next partition to refine");
			hipart = partitioner.maxPECut(tid);
			fprintf(stderr,"\n----Refining partition %d with max cuts\n", hipart);
                        partitioner.gainTable.clear();
			// Keep refining the partition until it reaches min edgecuts and all the bnd vtx list is exhausted
                        unsigned newtCuts = partitioner.countTotalPECut(tid);
//                        if(partitioner.bndIndMap[hipart].size() > 0 && newtCuts >=tCuts){
				partitioner.refinePart(tid, hipart, newtCuts);
//			}
			fprintf(stderr,"\nFinished refining the partition %d", hipart);
  //                      }
			// remove the refined partition
			auto id = partitioner.fetchPIds.find(hipart);
			partitioner.fetchPIds.erase(id);
			fprintf(stderr,"\nRemoving %d from fetchPIds size %d\n", *id, partitioner.fetchPIds.size());
		}

//        partitioner.printParts(tid);
		mr->partitioner.gCopy(tid);
        fprintf(stderr,"\nFinished IO refining, Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
	}
 //       partitioner.cread(tid);
	pthread_barrier_wait(&(mr->barRefine));
         if(tid == 0)
	  mr->afterRefine(tid, mr->nVertices);
	time_refine += getTimer();
	mr->refine_times[tid] += time_refine;

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

	  mr->beforeRefine(tid);
	  mr->refineInit(tid);
	  partitioner.initiateInMemoryRefine(tid); 
         partitioner.addDVals(tid);
 
 	 fprintf(stderr,"\nIn Memory REFINE- tid: %d, Computing edgecuts with Map Size %d", tid, partitioner.refineMap[tid].size());
         partitioner.ComputeBECut(tid);
	pthread_barrier_wait(&(mr->barRefine));
	// Count the total edge cuts and also check the partition with max edgecuts
	if(tid == 0){
	 partitioner.releaseInMemStructures();

		for (unsigned i = 0; i < mr->nThreads; i++){
			partitioner.fetchPIds.insert(i);  //Partition ids  - 0,1,2 .. 
		}

        unsigned tCuts = partitioner.countTotalPECut(tid);
        fprintf(stderr,"\nBefore refining, Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
		while(partitioner.fetchPIds.size() > 1){
			unsigned hipart = -1;
                        fprintf(stderr,"\nFinding Next partition to refine");
			hipart = partitioner.maxPECut(tid);
			fprintf(stderr,"\n----Refining partition %d with max cuts\n", hipart);
                        partitioner.gainTable.clear();

                        unsigned newtCuts = partitioner.countTotalPECut(tid);
				partitioner.refinePart(tid, hipart, newtCuts);
			fprintf(stderr,"\nFinished refining the partition %d", hipart);
			// remove the refined partition
			auto id = partitioner.fetchPIds.find(hipart);
			partitioner.fetchPIds.erase(id);
			fprintf(stderr,"\nRemoving %d from fetchPIds size %d\n", *id, partitioner.fetchPIds.size());
		}

  //      partitioner.printParts(tid);
		mr->partitioner.gCopy(tid);
        fprintf(stderr,"\nFinished In-mem refining, Total EdgeCuts: %d\n", partitioner.countTotalPECut(tid));
	}
        
//         partitioner.ComputeBECut(tid);
	pthread_barrier_wait(&(mr->barRefine));
         if(tid == 0)
	  mr->afterRefine(tid, mr->nVertices);
	  time_refine += getTimer();
	  mr->refine_times[tid] += time_refine;
	 
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
	partitioner.initg(nVertices, hDegree, nThreads, batchSize, kBItems, nParts); // GK 


	fprintf(stderr, "Partitioning input for Parallel Reads\n");
	partitionInputForParallelReads();

	end_read.resize(nThreads, 0.0);
	mparts_times.resize(nThreads, 0.0);
	combine_times.resize(nParts, 0.0);
	refine_times.resize(nParts, 0.0);



	fprintf(stderr, "Reading Graph from file\n");
	parallelExecute(doMParts, this, nThreads);
	fprintf(stderr,"\nSuccessfully Uploaded the Graph\n");

	//  fprintf(stderr, "Running Coarseners\n");
	//  parallelExecute(doCoarsen, this, nCoarseners);
	//      fprintf(stderr,"\nGP.HPP before Running refiners");

	 if(!partitioner.getWrittenToDisk()) {
	    fprintf(stderr, "Running InMemoryRefiners\n");
	    parallelExecute(doInMemoryRefine, this, nParts);
	  //  partitioner.releaseInMemStructures();
	    fprintf(stderr, "\nSuccess !!!!!!!!!!!!\n");
	    } else {
	      partitioner.releaseInMemStructures();
	fprintf(stderr, "\nRunning Combiners\n");
	parallelExecute(doCombine, this, nParts);

	partitioner.releaseReadPartStructures();
	fprintf(stderr, "\nRunning Refiners\n");
	parallelExecute(doRefine, this, nParts);
	 }

	fprintf(stderr, "Graph partitioned. Shutting down.\n");

	fprintf(stderr, "--------------------------------------\n");

	partitioner.shutdown();

	std::cout << "------- Final Time ---------" << std::endl;
	std::cout << " Total time : " << run_time << " (msec)" << std::endl;

	run_time += getTimer();

	auto mparts_time = max_element(std::begin(mparts_times), std::end(mparts_times)); 
	std::cout << " Partitioning time : " << *mparts_time/1000 << " (sec)" << std::endl;

	auto combine_time = max_element(std::begin(combine_times), std::end(combine_times));
	std::cout << " Combine time : " << *combine_time << " msec" << std::endl;

	auto refine_time = max_element(std::begin(refine_times), std::end(refine_times));
	std::cout << " Refine time : " << *refine_time << " msec" << std::endl;

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
//void GraphParts::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const unsigned int* mapSupRowstoRows)
void GraphParts::writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hidegree = 0){
	partitioner.writeBuf(tid, to, from, hidegree);
}

/*void GraphParts::coarsen(const unsigned tid, graph_t *graph, graph_t *cgraph, unsigned CoarsenTo)
  {
  partitioner.coarsen(tid, graph, cgraph, CoarsenTo);
//  partitioner.coarsen(tid, cgraph, CoarsenTo, numEdgesSupRowstoRows, mapSupRowstoRows);
}
 */
//--------------------------------------------
bool GraphParts::read(const unsigned tid) {
	return partitioner.read(tid);
}

//--------------------------------------------
void GraphParts::printParts(const unsigned tid, std::string outputPrefix) {
	return partitioner.printParts(tid, outputPrefix);
}

//--------------------------------------------
bool GraphParts::refine(const unsigned tid) {
	return partitioner.refine(tid);
}

//--------------------------------------------
void GraphParts::ComputeBECut(const unsigned tid) {
	return partitioner.ComputeBECut(tid);
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
