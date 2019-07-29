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
	infile.seekg(tid*mr->bytesPerFile);

	std::string line;

	//   std::vector<unsigned> where (mr->nVertices+1, -1);
	//   std::vector<unsigned>* where;
	//   std::vector<unsigned> whereDst (mr->nVertices+1, -1);    
	fprintf(stderr, "Creating memory partitions nVertices: %d, nEdges: %d\n", mr->nVertices, mr->nEdges);

	//std::getline(infile, line);
	//  unsigned long long bytesRead = 0; 
	// unsigned long long linesRead = 0; 

	//while(std::getline(infile, line) && bytesRead <= mr->bytesPerFile) {
	while(std::getline(infile, line)) {
		time_mparts -= getTimer();
		mr->createMParts(tid, line);     
		/*        std::stringstream inputStream(line);
			  unsigned to, from;
			  inputStream >> to;
		//             mr->writeBuf(tid, to, to);

		while(inputStream >> from){
		fprintf(stderr,"\nTID: %d picked TO: %zu FROM: %zu \n",tid, to, from);
		mr->writeBuf(tid, to, from);
		//       mr->writeBuf(tid, from, to);
		}
		// TODO         mr->writeBuf(tid, to, to);
		//          bytesRead = line.length();        
		 */          time_mparts += getTimer();
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
				fprintf(stderr,"\ntid: %d, Key: %d", tid, fit->first);
				//partitioner.totalCombined[tid]++;
			}
			partitioner.totalCombined[tid] += partitioner.readBufMap[tid].size();
			fprintf(stderr,"\n-- tid: %d, totalCombined: %d, Size: %d\n", tid, partitioner.totalCombined[tid], partitioner.readBufMap[tid].size());
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
			fprintf(stderr,"\n TID: %d, Key: %d", tid, rank);
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

		partitioner.totalCombined[tid] += mr->kBItems;
		fprintf(stderr,"\n-- tid: %d, totalCombined: %d\n", tid, partitioner.totalCombined[tid]);
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

//	mr->beforeRefine(tid);

	//  fprintf(stderr, "\nDoRefine: Initializing read parameter\n");
	mr->refineInit(tid);

	partitioner.cread(tid);
	/*while(true) {
	//  fprintf(stderr, "\nDoRefine: Calling Read\n");
	bool execLoop = mr->refine(tid);
	fprintf(stderr, "\nExecloop is %d", execLoop);

	if(execLoop == false) {
	for(InMemoryConstIterator it = partitioner.refineMap[tid].begin(); it != partitioner.refineMap[tid].end(); ++it){
	fprintf(stderr,"\nREFINE- tid: %d, Key: %d", tid, it->first);
	//        mr->refine(tid, it->first, it->second);
	mr->ComputeBECut(tid);
	partitioner.refineMap[tid].erase(partitioner.refineMap[tid].begin(), it);

	}
	break;
	}

	// Read the kitems from infinimem and compute the edgecuts for that part
	mr->ComputeBECut(tid);

	partitioner.refineMap[tid].erase(partitioner.refineMap[tid].begin(), partitioner.refineMap[tid].end());
	}
	 */
	pthread_barrier_wait(&(mr->barRead));
	// Count the total edge cuts and also check the partition with max edgecuts
	if(tid == 0){

		for (unsigned i = 0; i < mr->nThreads; i++){
			partitioner.fetchPIds.insert(i);  //Partition ids  - 0,1,2 .. 
		}

		while(partitioner.fetchPIds.size() > 1){
//			unsigned tCuts = 0;
//			unsigned chVtx = -1;
			unsigned hipart = -1;
//			tCuts = mr->countTotalPECut(tid);
//			fprintf(stderr,"\n----TID: %d, Total Edge cuts : %d\n", tid, tCuts);
                        fprintf(stderr,"\nFinding Next partition to refine");
			hipart = partitioner.maxPECut(tid);
			fprintf(stderr,"\n----Refining partition %d with max cuts\n", hipart);

//			unsigned pCut = partitioner.getTotalPECuts(hipart);
			// Keep refining the partition until it reaches min edgecuts and all the bnd vtx list is exhausted
//			while(partitioner.bndIndMap[hipart].size() > 0){
                        if(partitioner.bndIndMap[hipart].size() > 0){
				partitioner.refinePart(tid, hipart);
					// TODO:: changed where but how to move the vertex to another on disk? or may be not?
//				}
//			}
			fprintf(stderr,"\nFinished refining the partition %d", hipart);
                        }
			// remove the refined partition
			auto id = partitioner.fetchPIds.find(hipart);
			partitioner.fetchPIds.erase(id);
			fprintf(stderr,"\nRemoving %d from fetchPIds size %d\n", *id, partitioner.fetchPIds.size());
			// Calculate the next partition with max edgecuts
	//		hipart = partitioner.maxPECut(tid);
	//		fprintf(stderr,"\n----Next Partition with max cuts : %d\n", hipart);
			//  if(hipart < 0) break;
		}

	}
	//Check if moving the vertex with max bnd value reduce the num of cuts in the parition and how much does it impact in the current partition also check the total edgecuts
	  mr->afterRefine(tid, mr->nVertices);
        fprintf(stderr,"\nFinished refining ......\n");
  //      partitioner.printParts(tid);
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
 
 	 fprintf(stderr,"\nIn Memory REFINE- tid: %d, Computing edgecuts with Map Size %d", tid, partitioner.refineMap[tid].size());
         partitioner.ComputeBECut(tid);
	 partitioner.releaseInMemStructures();
	pthread_barrier_wait(&(mr->barRefine));
	// Count the total edge cuts and also check the partition with max edgecuts
	if(tid == 0){

		for (unsigned i = 0; i < mr->nThreads; i++){
			partitioner.fetchPIds.insert(i);  //Partition ids  - 0,1,2 .. 
		}

		while(partitioner.fetchPIds.size() > 1){
			unsigned hipart = -1;
                        fprintf(stderr,"\nFinding Next partition to refine");
			hipart = partitioner.maxPECut(tid);
			fprintf(stderr,"\n----Refining partition %d with max cuts\n", hipart);

                        if(partitioner.bndIndMap[hipart].size() > 0){
				partitioner.refinePart(tid, hipart);
			fprintf(stderr,"\nFinished refining the partition %d", hipart);
                        }
			// remove the refined partition
			auto id = partitioner.fetchPIds.find(hipart);
			partitioner.fetchPIds.erase(id);
			fprintf(stderr,"\nRemoving %d from fetchPIds size %d\n", *id, partitioner.fetchPIds.size());
		}

	}
	//Check if moving the vertex with max bnd value reduce the num of cuts in the parition and how much does it impact in the current partition also check the total edgecuts
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
	partitioner.initg(nVertices, nThreads, batchSize, kBItems, nParts); // GK 


	fprintf(stderr, "Partitioning input for Parallel Reads\n");
	partitionInputForParallelReads();

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
	    parallelExecute(doInMemoryRefine, this, nThreads);
	  //  partitioner.releaseInMemStructures();
	    } else {
	      partitioner.releaseInMemStructures();
	fprintf(stderr, "\nRunning Combiners\n");
	parallelExecute(doCombine, this, nThreads);

	partitioner.releaseReadPartStructures();
	fprintf(stderr, "\nRunning Refiners\n");
	parallelExecute(doRefine, this, nThreads);
	 }

	fprintf(stderr, "Graph partitioned. Shutting down.\n");

	fprintf(stderr, "--------------------------------------\n");

	//parts.shutdown();

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
	//  fprintf(stderr, "\nnumlines: %u\n", numLines);
	in.close();
	//fprintf(stderr,"\n NumLines: %zu", numLines);
	//linesPerFile = numLines/nThreads + 1;
	// fprintf(stderr,"\nLinesPerfile: %zu", linesPerFile);
	bytesPerFile = fileSizeInBytes/nThreads + 1; //fileSizeInBytes/nThreads + 1;
}

//--------------------------------------------
void GraphParts::init(const std::string input, const unsigned nvertices, const unsigned nedges, const unsigned nthreads, const unsigned nparts, const unsigned bSize, const unsigned kItems) {

	inputFileName = input;
	std::cout << "Input file name: " << inputFileName << std::endl;
	numLines = getNumLines(inputFileName);
	//  fprintf(stderr,"\nINIT NumLines: %zu\n", numLines);

	//nInMemParts and nParts are same
	nVertices = nvertices;
	nEdges = nedges;
	nThreads = nthreads;
	nParts = nparts;
	batchSize = bSize;
	kBItems = kItems;


	//  setRefiners(std::min(nThreads, nRefiners));
	//TODO:need to check if I need this --  nParts = std::min(nThreads, nParts);

	std::cout << "nVertices: " << nVertices << std::endl;
	std::cout << "nEdges: " << nEdges << std::endl;
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
void GraphParts::writeBuf(const unsigned tid, const unsigned to, const unsigned from){
	partitioner.writeBuf(tid, to, from);
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
