#include "partitioner.h"
//#include "coarsen.h"
//#include "kway.h"


/*#ifdef USE_ONE_PHASE_IO
#include "infinimem/onePhaseFileIO.hpp"
#else
#include "infinimem/twoPhaseFileIO.hpp"
#endif*/
#include "infinimem/twoPhaseFileIO.hpp"
#include "infinimem/onePhaseFileIO.hpp"

#include <vector>
#include <set>
#include <ctype.h> // for tolower()
#include <algorithm> // for std::sort()
#include <stack> // for lookuptable
#include <fstream>  // GK
#include<iostream> //GK
#include<iterator> //GK
#include <utility>//GK

//pthread_mutex_t lock_buffer = PTHREAD_MUTEX_INITIALIZER;
//------------------------------------------------- GK
// Initialize in-memory Buffers
void Partitioner::initg(unsigned nVertices, unsigned hDegree, unsigned nThreads, unsigned bSize, unsigned kItems, unsigned nParts)
{
  //nBuffers = pow(buffers, 2); // number of buffers is square of number of threads
  //nBuffers = nMappers * nReducers;

  nVtces = nVertices;
  nRows = nThreads;
  hiDegree = hDegree;
  nCols = nParts; //TODO: should be two by default
  writtenToDisk = false;
  partRefine = false;
  batchSize = bSize;
  kBItems = kItems;
  cbartid = 0;
//  nparts = nParts;
 // nvertices = nVertices;

  totalPECuts = new IdType[nCols];
  totalKeysInFile = new IdType[nCols];
  totalCombined = new IdType[nCols];
  nItems = new IdType[nRows * nCols];
  //nCuts = new IdType[nCols];
//  nEdges = new IdType[nRows * nCols];
    where = new std::vector<unsigned>[nCols];
//  partitionBndInd = new std::vector<unsigned>[nCols];
//  partitionBndPtr = new std::vector<unsigned>[nCols];
  
  outBufMap = new InMemoryContainer[nRows * nCols];
  readBufMap = new InMemoryContainer[nCols];
  bndIndMap = new LookUpTable[nCols];
  lookUpTable = new LookUpTable[nCols];
//  gainTable = new InMemTable[nCols];
  dTable = new InMemTable[nCols];
  fetchBatchIds = new std::set<unsigned>[nCols];
  fetchPIds = new std::set<unsigned>[nCols];
  markMax = new std::map<unsigned, unsigned>[nCols];
  markMin = new std::map<unsigned, unsigned>[nCols];
  readNextInBatch = new std::vector<unsigned long long>[nCols];
  batchesCompleted = new std::vector<bool>[nCols];
  pIdsCompleted = new std::vector<bool>[nCols];
  keysPerBatch = new std::vector<unsigned>[nCols];
  
  pthread_barrier_init(&barRefinepart, NULL, nParts);
  for (unsigned i=0; i<nRows * nCols; ++i){ 
    nItems[i] = 0;
    }

  for(unsigned i=0; i<nCols; ++i) {  
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    locks.push_back(mutex);
    
    totalKeysInFile[i] = 0;
  }

  // setup FileIO
//#ifdef USE_ONE_PHASE_IO
  //io = new OnePhaseFileIO<RecordType>("/tmp/gkaur007/mrdata/", nCols, 0/*UNUSED*/);
//#else
  io = new TwoPhaseFileIO<RecordType>("/tmp/gkaur007/mrdata/", nCols, 0/*UNUSED*/);
  cio = new TwoPhaseFileIO<RecordType>("/tmp/gkaur007/combdata/", nCols, 0/*UNUSED*/);
//#endif
}

//-------------------------------------------------
void Partitioner::releaseInMemStructures()
{
  for (unsigned i = 0; i < nCols; i++)
    pthread_mutex_destroy(&locks[i]);

   for (unsigned i = 0; i < nRows * nCols; i++)
    outBufMap[i].clear();

  delete[] nItems;
//  delete[] nCuts;
//  delete[] nEdges;
  delete[] outBufMap;
//  delete[] bndind;
//  delete[] bndptr;
}
//-------------------------------------------------
void Partitioner::shutdown()
{

  //gainTable->clear();
  dTable->clear();
  delete cio;
   for (unsigned i = 0; i < nCols; i++){
    bndIndMap[i].clear();
  }
  delete[] totalPECuts;
  delete[] totalCombined;
  delete[] bndIndMap;
  delete[] dTable;
}

//-------------------------------------------------
void Partitioner::releaseReadPartStructures()
{
  delete io;
  fetchBatchIds->clear();
  fetchPIds->clear();
  markMax->clear();
  markMin->clear();
  batchesCompleted->clear();
  pIdsCompleted->clear();
  lookUpTable->clear();
  keysPerBatch->clear();
   for (unsigned i = 0; i < nCols; i++){
       readBufMap[i].clear();
       readNextInBatch[i].clear();
  }


  delete[] readBufMap;
  delete[] readNextInBatch;
  delete[] fetchBatchIds;
  delete[] fetchPIds;
  delete[] markMax;
  delete[] markMin;
  delete[] batchesCompleted;
  delete[] pIdsCompleted;
  delete[] keysPerBatch;
  delete[] lookUpTable;
  delete[] totalKeysInFile;
}

//--------------------------------------------------
void Partitioner::writeInit(const unsigned tid) {
//  unsigned j=0;
  for (unsigned i = 0; i<=nVtces; ++i) {
//       partitionBndInd[tid].push_back(NULL);
//       partitionBndPtr[tid].push_back(-1);
         where[tid].push_back(-1); 
         if(tid == 0){
           gWhere.push_back(-1);
         } 
 //TODO : not sure if I need to initialize the below loop
//    pthread_mutex_lock(&locks[tid]);
//       gWhere[i] = 1;
//    pthread_mutex_unlock(&locks[tid]);
  }

}


//-------------------------------------------------
//void Partitioner::writeBuf(const unsigned tid, const unsigned to, const unsigned from, std::vector<unsigned>& where){
void Partitioner::writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hIdSize = 0){
    double timeWBF = -getTimer();

    unsigned bufferId = hashKey(to) % nCols;
    unsigned buffer = tid * nCols + bufferId; 
    if(hIdSize != 0){
      hIds.emplace(to, hIdSize); 
      // Put the records into different buffer randomly
      unsigned bufferId = hashKey(from) % nCols; 
      unsigned buffer = tid * nCols + bufferId; 
// less chances of adjlist vertices being on boundary if they are hashed to the partition which has similar vertices as key; hence less edgecuts 
   }

    if(where[tid].at(to) == -1){
      where[tid].at(to) = bufferId;
  }
//  unsigned buffer = tid * nCols + bufferId;  
    unsigned whereFrom = hashKey(from) % nCols; 
    if(where[tid].at(from) == -1){
      where[tid].at(from) = whereFrom;
  }

  //fprintf(stderr,"\nwhere[%d]: %d, where[%d]: %d, bufferID: %d\n\n", to, where[to], from, where[from], buffer);

  if (outBufMap[buffer].size() >= batchSize) {

    infinimem_write_times[tid] -= getTimer();     
    pthread_mutex_lock(&locks[bufferId]);
//    fprintf(stderr,"\nTID %d outbufmap buffer %d full with %d records noItems:%d\n", tid, buffer, outBufMap[buffer].size(), nItems[buffer]);
   
    writeToInfinimem(bufferId, totalKeysInFile[bufferId], outBufMap[buffer].size(), outBufMap[buffer]);
    totalKeysInFile[bufferId] += nItems[buffer];
    pthread_mutex_unlock(&locks[bufferId]);
    infinimem_write_times[tid] += getTimer();

    outBufMap[buffer].clear();
    nItems[buffer] = 0;
    writtenToDisk = true;
  }
 
  performWrite(tid, buffer, to, from);

  timeWBF += getTimer();
  writeBuf_times[tid] += timeWBF;
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

fprintf(stderr,"\nFlushing buffer %d residues\n", tid);

if(tid >= nCols) 
    return;

double frTime = -getTimer();
assert(writtenToDisk);

  if(nRows == 1) {
    infinimem_write_times[tid] -= getTimer();
    writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
    infinimem_write_times[tid] += getTimer();
    outBufMap[tid].clear();
    totalKeysInFile[tid] += nItems[tid];
    nItems[tid] = 0;
  }
  else {
    unsigned b1 = 0, b2 = 0;
    InMemoryContainerIterator b2Iter, b2End;
    unsigned long long b2Merged = 0;
    bool findB1 = true, findB2 = true;
    unsigned i = 0;
    while(i < nRows - 1) {
      if(findB1) {
        b1 = (i * nCols) + tid;
        findB1 = false;
        ++i;
      }
      if(findB2) {
        b2 = (i * nCols) + tid;
        b2Iter = outBufMap[b2].begin();
        b2End = outBufMap[b2].end();
        b2Merged = 0;
        findB2 = false;
        ++i;
      }

      b2Merged += merge(outBufMap[b1], b1, tid, b2Iter, b2End);

      if(outBufMap[b1].size() == 0) {
        findB1 = true;
      }
      if(b2Iter == b2End) {
        findB2 = true;
      }
    }

    if(i == nRows - 1) {
      if(findB1 && findB2) {
        infinimem_write_times[tid] = -getTimer();
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
        infinimem_write_times[tid] += getTimer();
        outBufMap[i].clear();
        totalKeysInFile[tid] += nItems[i];
        nItems[i] = 0;
      }
      else if(findB1 || findB2) {
        if(findB1) {
          b1 = (i * nCols) + tid;
          findB1 = false;
          ++i;
        } else {
          b2 = (i * nCols) + tid;
          b2Iter = outBufMap[b2].begin();
          b2End = outBufMap[b2].end();
          b2Merged = 0;
          findB2 = false;
          ++i;
        }  

        b2Merged += merge(outBufMap[b1], b1, tid, b2Iter, b2End);
        if(outBufMap[b1].size() == 0) {
          if(b2Iter != b2End) {
            infinimem_write_times[tid] = -getTimer();
            bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
	    infinimem_write_times[tid] += getTimer();
            outBufMap[b2].clear();
            totalKeysInFile[tid] += (nItems[b2] - b2Merged);
            nItems[b2] = 0;
          }
        } 
        if(b2Iter == b2End) {
          if(outBufMap[b1].size() != 0) {
     	    infinimem_write_times[tid] = -getTimer();
            writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
	    infinimem_write_times[tid] += getTimer();
            outBufMap[b1].clear();
            totalKeysInFile[tid] += nItems[b1];
            nItems[b1] = 0;
          }
        }
      } else
          assert(false);
    } else if(i == nRows) {
      if(outBufMap[b1].size() > 0) {
	infinimem_write_times[tid] = -getTimer();
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
	infinimem_write_times[tid] += getTimer();
        outBufMap[b1].clear();
        totalKeysInFile[tid] += nItems[b1];
        nItems[b1] = 0;
      } else {  
        infinimem_write_times[tid] = -getTimer();
        bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
	infinimem_write_times[tid] += getTimer();
        outBufMap[b2].clear();
        totalKeysInFile[tid] += (nItems[b2] - b2Merged);
        nItems[b2] = 0;
      }
    } else
      assert(false);
  }

  frTime += getTimer();
  flushResidues_times[tid] += frTime;
}

//--------------------------------------------------
unsigned long long Partitioner::merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end) {
  unsigned long long ct = 0;
  while(begin != end) {
    if(toMap.size() >= batchSize) {
      infinimem_write_times[tid] -= getTimer();
      writeToInfinimem(tid, totalKeysInFile[tid], toMap.size(), toMap);
      infinimem_write_times[tid] += getTimer();
      toMap.clear();
      totalKeysInFile[tid] += nItems[whichMap];
      nItems[whichMap] = 0; 
      return ct;
    }
    InMemoryContainerIterator it = toMap.find(begin->first);
    if(it != toMap.end()){
      combine(it->first, it->second, begin->second);
      ++localCombinedPairs[tid];
    }
    else {
      toMap.emplace(begin->first, begin->second);
      ++nItems[whichMap];
    }
    ++ct;
    ++begin;
  }
  return ct;
}

//--------------------------------------------------
void Partitioner::performWrite(const unsigned tid, const unsigned buffer, const unsigned to, const unsigned from) {
  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it_to = outBufMap[buffer].find(to); 
  IdType edge;
//  InMemoryContainerIterator it_dst = std::find(outBufMap[buffer].begin(), outBufMap[buffer].end(), e); 

  if(it_to != outBufMap[buffer].end()){
     combine(to, it_to->second, nbrs);
     ++localCombinedPairs[tid];
   //  fprintf(stderr,"\nTID %d Combining for key: %d, value: %d", tid, to, from);
//     outBufMap[buffer][to].push_back(from);
//     nEdges[buffer]++;
      }
  else {
    outBufMap[buffer].emplace(to, nbrs); 
    nItems[buffer]++;
  //  nEdges[buffer]++;
  }
}

//--------------------------------------------------
void Partitioner::writeToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, const InMemoryContainer& inMemMap) {
  RecordType* records = new RecordType[noItems]; 
  unsigned ct = 0;
//  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {
//      fprintf(stderr,"\nBuffer: %d full with startkey: %d, noItmes:%d \t", buffer, startKey, inMemMap.size());
    
  /*    for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit)
      fprintf(stderr,"%d\t", *vit); 

  }
 */ for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {
      records[ct].set_rank(it->first);
//      fprintf(stderr,"\n \n WTI- TID: %d,  Key: %d\t, Values: ", buffer, it->first); 
      for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
//      fprintf(stderr,"added %d\t", *vit); 
          records[ct].add_nbrs(*vit);
//      fprintf(stderr,"ct %d\t", ct); 
      }
   ++ct;

  }
    if (ct != noItems){
  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {              //fprintf(stderr,"\n \n WTI- TID: %d,  Key: %d\n", buffer, it->first);
     } 
//     fprintf(stderr,"\nBuffer: %d CT: %d noItems: %d InMemMap Size: %d\n", buffer, ct, noItems, inMemMap.size());
}
  assert(ct == noItems);
  //    fprintf(stderr,"\nWRITING TO INFINIMEM buffer: %d, startKey: %d, totalKeys: %d \n", buffer, startKey, totalKeysInFile[buffer]); 
  io->file_set_batch(buffer, startKey, noItems, records);
    //  fprintf(stderr,"\nBUFFER %d WRITTEN %d records\n", buffer, noItems); 

  delete[] records;
}

//--------------------------------------------------
void Partitioner::bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end) {
  RecordType* records = new RecordType[noItems]; 
  unsigned ct = 0;

   //  fprintf(stderr,"\n BWTI- TID: %d, startKey: %d ", buffer, startKey); 
  for (InMemoryConstIterator it = begin; it != end; ++it) {
     records[ct].set_rank(it->first);
  //    fprintf(stderr,"\n BWTI- TID: %d, startKey: %d, Key: %d\t, Values: ", buffer, startKey, it->first); 

    for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
      records[ct].add_nbrs(*vit);
    //  fprintf(stderr,"%d\t", *vit); 
      }
      ++ct;
      totalCombined[buffer]++;
  }

  assert(ct == noItems);
  io->file_set_batch(buffer, startKey, noItems, records);
 delete[] records;
}

//--------------------------------------------------
void Partitioner::readInit(const unsigned tid) {
  unsigned j=0;
  for (unsigned long long i = 0; i <= totalKeysInFile[tid]; i+= batchSize) {
    readNextInBatch[tid].push_back(i); //start position to read
    fetchBatchIds[tid].insert(j++);  //batch ids inside partition - 0,1,2 .. 
    batchesCompleted[tid].push_back(false); 
    keysPerBatch[tid].push_back(kBItems); //Keys read in each batch
  }
    totalCombined[tid] = 0;
//  std::vector<unsigned> bndptr (nVtces+1,-1); 
    refineInit(tid);
}


//--------------------------------------------------
void Partitioner::refineInit(const unsigned tid) {
  for (unsigned i = 0; i < nCols; i++) {
    fetchPIds[tid].insert(i);  //partition ids - 0,1,2 .. 
        pIdsCompleted[tid].push_back(false);
  }
  bndIndMap[tid].clear();
  totalPECuts[tid] = 0; 
}

//--------------------------------------------------
bool Partitioner::read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted) {

//  fprintf(stderr,"\nInside Partitioner::read \n");
  infinimem_read_times[tid] -= getTimer();

  RecordType* records = new RecordType[kBItems];
  unsigned batch = 0;
  for(auto it = fetchBatchIds.begin(); it != fetchBatchIds.end(); ++it) {
    batch = *it ;
    unsigned long long batchBoundary = std::min(static_cast<unsigned long long>((batch + 1) * batchSize), static_cast<unsigned long long>(totalKeysInFile[tid]));

    if (readNextInBatch[batch] >= batchBoundary) {
      batchesCompleted[batch] = true;
      continue;
    }

    keysPerBatch[batch] = std::min(keysPerBatch[batch], static_cast<unsigned>(batchBoundary - readNextInBatch[batch]));

    if (keysPerBatch[batch] > 0 && readNextInBatch[batch] < batchBoundary)
      io->file_get_batch(tid, readNextInBatch[batch], keysPerBatch[batch], records); 

    for (unsigned i = 0; i < keysPerBatch[batch]; i++) {
      lookUpTable[records[i].rank()].push_back(batch);
      readBufMap[records[i].rank()];
//      fprintf(stderr,"\nREAD- TID: %d, Key: %d\t Values: ", tid, records[i].rank()); 

      for (unsigned k = 0; k < records[i].nbrs_size(); k++){
        readBufMap[records[i].rank()].push_back(records[i].nbrs(k));

 //     fprintf(stderr,"%d\t", records[i].nbrs(k)); 
     }
//     dTable[tid][records[i].rank()] = records[i].nbrs_size();
    }

    readNextInBatch[batch] += keysPerBatch[batch];
    keysPerBatch[batch] = 0;

    if (readNextInBatch[batch] >= batchBoundary)
      batchesCompleted[batch] = true;
  }

  fetchBatchIds.clear();

  bool ret = false;
  for (unsigned readBatch = 0; readBatch < batchesCompleted.size(); readBatch++)
    if (batchesCompleted[readBatch] == false) {
      ret = true;
      break;
    }

  infinimem_read_times[tid] += getTimer();

  delete[] records;
  return ret;
}

//--------------------------------------------------
bool Partitioner::read(const unsigned tid) {
  return read(tid, readBufMap[tid], keysPerBatch[tid], lookUpTable[tid], fetchBatchIds[tid], readNextInBatch[tid], batchesCompleted[tid]);
}

//--------------------------------------------------
void Partitioner::initiateInMemoryRefine(unsigned tid) {
// fprintf(stderr,"\nInitiating In-Memory refine TID %d\n", tid);
   for(unsigned i=0; i<nRows; ++i) {
     readBufMap[tid].insert(outBufMap[tid + nCols * i].begin(), outBufMap[tid + nCols * i].end());
   }
}

//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid, const InMemoryContainer& inMemMap) {
     ComputeBECut(tid, gWhere, bndIndMap[tid], inMemMap);
}
//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid, const std::vector<unsigned>& where, LookUpTable& bndind, const InMemoryContainer& inMemMap) {

     unsigned src;
    // bool first = 1;
     std::vector<unsigned> bndvert;     

  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {
      src = it->first;
      for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
       //   if (std::find(bndvert.begin(), bndvert.end(), *vit) == bndvert.end()){
              bndvert.push_back(*vit);
         // }
       }
         IdType ebnd=0;
     for(unsigned i=0; i<bndvert.size(); i++){
      //   fprintf(stderr,"ComputeCut- i: %d\n", i);
         IdType dst = bndvert[i];
        //   fprintf(stderr,"\nTID: %d, src :%d, dst: %d", tid, src, dst); 
        //   fprintf(stderr,"\n"); 
  //               fprintf(stderr,"where[%d]: %d, where[%d]: %d\n", src, where[src], dst ,where[dst]);
                 if( where[dst] != -1 && where[src] != where[dst] ) {
               //  if( where[src] != where[dst] ) {
           //         fprintf(stderr,"Edge CUT tid: %d\n", tid);
                    ebnd++;
    		    totalPECuts[tid]++;
                   // dst = adjncy[j];
          //   fprintf(stderr,"\nBoundary Vertices tid: %d, dst: %d\n", tid, dst);
          // store the boundary vertices along with their actual location
              bndind[dst].push_back(where[dst]); // TODO: should this be where[dst]?
  //   fprintf(stderr,"\nSize bndind[%d]: %d in tid %d \n", dst, bndind[dst].size(), tid);
              //      }
                 }
             // }
          }
 /*         if(xadj[src] == xadj[src+1]){
  */	bndvert.clear();
     }
     fprintf(stderr,"\ntPECuts[%d]: %d\n", tid, totalPECuts[tid]/2);

}

//--------------------------------------------------
void Partitioner::setNum(const unsigned tid, std::vector<unsigned>& where, unsigned num){
     for(unsigned i=0; i<nVtces; ++i){
         where[i] = num;
     }
}

//--------------------------------------------------
void Partitioner::gCopy(const unsigned tid){
     gCopy(tid, gWhere);
 
 }
//--------------------------------------------------
void Partitioner::gCopy(const unsigned tid, std::vector<unsigned>& gWhere){
      bool first = 1;
      for(unsigned i=0; i<nCols; ++i){
         for(unsigned j=0; j<=nVtces; ++j){
             if(first){  // All Values of first thread will be copied
//         fprintf(stderr,"\nwhere[%d][%d]: %d,", i, j, where[i][j]);
                gWhere[j] = where[i][j];
             }
             else {
                if(where[i][j] != -1){
                   gWhere[j] = where[i][j];
               }
            }
  //       fprintf(stderr,"\nGWHERE[%d]: %d", j, gWhere[j]);
         }
         first = 0;
     }
}

//--------------------------------------------------
unsigned Partitioner::countTotalPECut(const unsigned tid) {
      totalCuts = 0;
      for(unsigned i=0; i<nCols; i++){
      //for(auto i=fetchPIds.begin(); i != fetchPIds.end(); ++i){
          totalCuts += totalPECuts[i];
      }
   return (totalCuts/2);
}

//--------------------------------------------------
unsigned Partitioner::maxPECut(const unsigned tid) {
/*      unsigned low = 0, hipart = -1;
//    fprintf(stderr,"\nMaxPECut fetchpids size %d", fetchPIds.size());
      for(unsigned i=0; i<nCols; i++){
//      for(auto i=fetchPIds.begin(); i != fetchPIds.end(); ++i){
          if(low < totalPECuts[*i]){
             low = totalPECuts[*i];
             hipart = *i;
          }
      }
   fprintf(stderr,"\n Partition: %d has max cuts: %d", hipart, (totalPECuts[hipart]/2));
   return hipart;
*/ }

//--------------------------------------------------
void Partitioner::updateDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax, unsigned src, unsigned dst){
//Go through adjacency (boundary) vertices of the masked vertices and update their DVals
  auto it_map = readBufMap[hipart].find(src);
 fprintf(stderr,"\nUpdateDVAL src: %d  size %d ", src, it_map->second.size());
  for ( auto vit = it_map->second.begin(); vit != it_map->second.end(); ++vit ) {
       unsigned adjvtx = *vit;
     //    fprintf(stderr,"\nCHECK VIT %d ", *vit);
       if(gWhere.at(*vit) == whereMax){ // the vertex should belong to the other partition being refined with
      // fprintf(stderr,"\nVIT: %d ", adjvtx);
   /*       auto it_dt = dTable[whereMax].find(adjvtx);
	  if(it_dt == dTable[whereMax].end()){
             continue;
          }
     */  // these vertices are connected to src .. need to find their connect with dst
          auto it_dst = readBufMap[whereMax].find(dst); //findSrcMap(tid, readBufMap[whereMax], dst);    //readBufMap[whereMax].find(dst);
          unsigned conDst = 0;
         // bool stat = findSrcVec(tid, 
         // if (std::find(readBufMap[whereMax][it_dst].begin(), readBufMap[whereMax][it_dst].end(), *vit) == readBufMap[whereMax][it_dst].end()){
          if (std::find(it_dst->second.begin(), it_dst->second.end(), *vit) == it_dst->second.end()){
        
	      conDst = 0;
         }        
         else
              conDst = 1;
		
             unsigned connect = conDst - 1;
        	unsigned currval = dTable[whereMax].at(*vit);    
		int dVal = currval + 2 * connect;
           // fprintf(stderr,"\nSRC: %d DST: %d vtx: %d CONNECT: %d dVal: %d ", src, dst, *vit, connect, dVal);
                dTable[whereMax].at(adjvtx) = dVal;
             
          // fprintf(stderr,"\nUpdated DVAL for vertex %d is %d   ", *vit, dTable[whereMax].at(*vit));
          }
  }

}

//--------------------------------------------------
void Partitioner::computeDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax, const unsigned long long k) {
//Go through each key in the selected partition to be refined and update the DVals
  auto begin = std::next(dTable[hipart].begin(), k);
 fprintf(stderr,"\nTID %d ComputeDVAL BEGIN: %d , k: %d ", tid, *begin, k);
  for (auto it = begin; it != dTable[hipart].end(); ) {
      unsigned src = it->first;

      bool stat = findSrc(tid, markMax[hipart], src);
      if (stat == true ){              
         ++it;  
         continue;
	}
      else{
     //  fprintf(stderr,"\nSRC: %d ", src);
//        fprintf(stderr,"\n-- tid: %d, k: %d hipart: %d whereMax: %d", tid, k, hipart, whereMax);
           //   auto it_bnd = bndIndMap[whereMax].find(src); 
           //     if(it_bnd != bndIndMap[whereMax].end()){
                if(bndIndMap[whereMax][src].size() > 0){
                   unsigned inDeg = dTable[hipart][src] - bndIndMap[whereMax][src].size();                   int dVal = bndIndMap[whereMax][src].size() - inDeg;
   
                      dTable[hipart].at(src) = dVal;
  //   fprintf(stderr,"\nDVAL for vertex %d is %d in part %d \n", src, dTable[hipart].at(it_bnd->first), hipart);
			++it;
 		}
               else {
               // not a boundary vertex
              //  fprintf(stderr,"\nDVAL src %d is not a bndry vertex ", src);
                dTable[hipart].erase(it++);
               }
     
   }
   }
}

//--------------------------------------------------
bool Partitioner::findSrc(const unsigned tid, std::map<unsigned, unsigned>& mark, unsigned src){
     bool stat = false;
     for(std::map<unsigned, unsigned>::const_iterator it = mark.begin(); it!=mark.end(); ++it){
     //  fprintf(stderr,"\nFinding src: %d ", src);
        if(src == it->first || src == it->second)
           stat = true;
           break;
    }

  return stat;
}

//--------------------------------------------------
bool Partitioner::findSrcVec(const unsigned tid, const std::vector<unsigned>& val, unsigned src){
     bool stat = false;
     for(unsigned i = 0; i < val.size(); ++i){
     //  fprintf(stderr,"\nFinding src: %d ", src);
        if(src == val[i]){
           stat = true;
           break;
        }
    }

  return stat;
}
//--------------------------------------------------
unsigned findSrcMap(const unsigned tid, InMemoryContainer& inMemMap, unsigned src){
   unsigned ind = -1;
   unsigned counter = 0; 
   unsigned mid = inMemMap.size()/2;
   std::vector<unsigned> val;
//   InMemoryContainerIterator begin = inMemMap.begin();
//   InMemoryContainerIterator end = inMemMap.end();
   if(src < std::next(inMemMap.begin(),mid)->first){
      auto begin = inMemMap.begin();
      auto end = std::next(inMemMap.begin(), mid);
    }
   else{
   	auto begin = std::next(inMemMap.begin(), mid);
   	auto end = inMemMap.end();
//      fit = std::next(inMemMap.begin(), ind);
     for (auto fit = begin ; fit != end; ++fit ) {
	 if(src == fit->first){
            ind = counter; //val = fit->second;
         fprintf(stderr, "\n Index in Map %d ", ind);
           break;           
        }
         ++counter;
     }

  }
 return ind;
}
//--------------------------------------------------
unsigned Partitioner::computeGain(const unsigned tid, const unsigned hipart, const unsigned whereMax, const unsigned long long k, const InMemoryContainer& inMemMap){
  int ct = -1;
  int  maxG = 0;
  int maxvtx = -1, minvtx = -1;
  unsigned vtx_ind = -1;
 fprintf(stderr,"\nComputeGAIN tid: %d ", tid, k );
//Go through each key in the selected partition to be refined and update the gain 
  for (auto it = dTable[hipart].begin(); it != dTable[hipart].end(); ++it) {
      unsigned src = it->first; //maxvtx;
      std::map<unsigned, unsigned>::const_iterator it_max = markMax[hipart].find(src);
      if(it_max != markMax[hipart].end()){ 
/*      bool stat = findSrc(tid, markMax[hipart], src);
      if (stat == true ){              
       */  continue;
	}
      else{
         if(dTable[whereMax].size() > 0 ){
         auto begin = std::next(dTable[whereMax].begin(), k);
         for (auto it_hi = begin; it_hi != dTable[whereMax].end(); ++it_hi) {
            unsigned dst = it_hi->first;
            bool connect = 0;
	    unsigned dsrc = dTable[hipart][src];
            unsigned ddst = dTable[whereMax][dst];
          if(dsrc >= 0 || ddst >= 0){// || (dsrc >=0 && ddst > 0) || (dsrc > 0 && ddst >=0)){  
           std::map<unsigned, unsigned>::const_iterator it_min = markMin[hipart].find(dst);
           if(it_min != markMin[hipart].end()){
             continue;
	    }
           /*bool stat2 = findSrc(tid, markMax[hipart], dst);//MarkMax has both src and dst
           if (stat2 == true ){              
              continue;
	   }*/
	   else{
                ct++;
                
               auto it_map = inMemMap.find(src);
             //const std::vector<unsigned>& val = findSrcMap(tid, inMemMap, src);
             
              if (std::find(it_map->second.begin(), it_map->second.end(), dst) == it_map->second.end()){
              // std::vector<unsigned> val = it_map->second;
              // bool statDst = findSrcVec(tid, val, dst);
              //  if(statDst != true){
                  // if(gWhere[src] != gWhere[dst]){
			connect = 0;
                 }      
                 else
                       connect = 1;
		
    //     fprintf(stderr,"\ntid: %d SRC: %d DST: %d CONNECT: %d ", tid, src, dst, connect);
                  if(!connect)      
             	    gainTable[ct] = dsrc + ddst;
           	  else
              	    gainTable[ct] = dsrc + ddst - 2;

                 int currGain = gainTable[ct];
             if(currGain > maxG){
//         fprintf(stderr,"\nGainTable: %d MaxG: %d ", gainTable[ct], maxG);
               maxG = currGain;
               maxvtx = src;
               minvtx = dst;
               vtx_ind = ct;
             }
   //   fprintf(stderr,"\nGAIN for vertex %d with dst %d is %d\n", src, dst, gainTable[ct]);
              }
            }
	 }
        }
       }
    }
             
    if(maxvtx != -1 && minvtx != -1){
  //  fprintf(stderr,"\nMASKING %d and %d with max gain %d \n", maxvtx, minvtx, maxG);
    markMax[hipart][maxvtx] = minvtx;  
    markMin[hipart][minvtx] = maxG;  
    if(!getWrittenToDisk()){
      updateDVals(tid, hipart, whereMax, maxvtx, minvtx); 
      updateDVals(tid, whereMax, hipart, minvtx, maxvtx);
    }
    return maxG;
   }
return -1;
}

//--------------------------------------------------
void Partitioner::bRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, const int ret) {
  //   fprintf(stderr,"\n Barrier threads %d ", cbartid); 
//     pthread_barrier_init(&barRefinepart, NULL, cbartid);
     bRefine(tid, hipart, whereMax, gWhere, markMax[hipart], markMin[hipart], ret);

}

//--------------------------------------------------
void Partitioner::bRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gWhere, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin, const int ret) {

  int maxGain = -1;
  unsigned long long k = 0;
  while(true) {
                //  fprintf(stderr, "\nDoRefine: Calling Read\n");
     bool execLoop = read(tid);
     if(execLoop == false) {
       

     for (InMemoryContainerIterator fit = readBufMap[tid].begin(); fit != readBufMap[tid].end(); ++fit) {
      //   fprintf(stderr,"\nVTx %d in memory ", fit->first);
     		dTable[tid][fit->first] = fit->second.size();
	}
//     fprintf(stderr,"\n-- tid: %d, readMap HipartSize: %d, k: %d hipart: %d whereMax: %d\n", tid, readBufMap[tid].size(), k, hipart, whereMax);
     maxGain = computeInfo(tid, hipart, whereMax, ret, k, maxGain);
     if (maxGain == -1){
        fprintf(stderr,"\n OR Partition %d is refined with whereMax partition %d ***** ", hipart, whereMax);
      }
     if(ret != 2)
       readBufMap[tid].erase(readBufMap[tid].begin(), readBufMap[tid].end());
     break;
     }

     unsigned counter = 0;
     InMemoryContainerIterator it;
     for (it = readBufMap[tid].begin(); it != readBufMap[tid].end(); ++it) {
          if (counter >= kBItems)
              break;

//         fprintf(stderr,"\nVTx %d in memory  size %d", it->first, it->second.size());
     		dTable[tid][it->first] = it->second.size();
          const unsigned rank = it->first;
     dTable[tid][rank] = it->second.size();
                //      fprintf(stderr,"\n TID: %d, Key: %d", tid, rank);
          auto pos = lookUpTable[tid].find(rank);
          assert(pos != lookUpTable[tid].end());
          const std::vector<unsigned>& lookVal = pos->second; // find the batch of the vertex
          for(unsigned val=0; val<lookVal.size(); val++) {
              fetchBatchIds[tid].insert(lookVal[val]);
              keysPerBatch[tid][lookVal[val]] += 1;
          }
          lookUpTable[tid].erase(rank);
          counter++;
     }

     maxGain = computeInfo(tid, hipart, whereMax, ret, k, maxGain);
     if (maxGain == -1){
        fprintf(stderr,"\n Partition %d is refined with whereMax partition %d ***** ", hipart, whereMax);
          break;
      }
     //Write combined records to a new partition    
     //cWrite(tid, kBItems, it);

     if(ret != 2)
     readBufMap[tid].erase(readBufMap[tid].begin(), it);
     }

   // pthread_mutex_lock(&locks[tid]);
     if(ret == 1) {
        writePartInfo(tid, hipart, whereMax);
   // pthread_mutex_unlock(&locks[tid]);
     }
}

//--------------------------------------------------
int Partitioner::computeInfo(const unsigned tid, const unsigned hipart, const unsigned whereMax, const int ret, unsigned long long k, int maxGain){
     if(ret != 2)
        ComputeBECut(tid, readBufMap[tid]);
     pthread_barrier_wait(&(barRefinepart));
     if(ret != 2)
     computeDVals(tid, hipart, whereMax, k);
     pthread_barrier_wait(&(barRefinepart));
    //fprintf(stderr, "\ntid: %d, dTableSize: %d  ", tid, dTable[hipart].size());
    
//     fprintf(stderr,"\n**********************TID %d , ret: %d  ******** ", tid, ret);
     if(ret == 1) {
        do{
              gainTable.clear();
              maxGain = computeGain(tid, hipart, whereMax, k, readBufMap[hipart]);
	} while(maxGain > 0);
     }
     if(ret != 2)
        k += dTable[hipart].size();

return maxGain;
}
//--------------------------------------------------
void Partitioner::inMemRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, const int ret) {
    // pthread_barrier_init(&barRefinepart, NULL, cbartid);
     inMemRefine(tid, hipart, whereMax, gWhere, markMax[hipart], markMin[hipart], ret);

}

//--------------------------------------------------
void Partitioner::inMemRefine(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gWhere, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin, const int ret) {

  int maxGain = -1;
  unsigned long long k = 0;
  while(true){
   if(readBufMap[tid].size() <= kBItems){
     addDVals(tid);
     maxGain = computeInfo(tid, hipart, whereMax, ret, k, maxGain);
     if (maxGain == -1){
        fprintf(stderr,"\n Last Partition %d is refined with whereMax partition %d ***** ", hipart, whereMax);
      }
      break;
   }

   InMemoryContainerIterator it;
   unsigned counter = 0;
   for (it = readBufMap[tid].begin(); it != readBufMap[tid].end(); ++it) {
       if (counter >= kBItems)
          break;

        dTable[tid][it->first] = it->second.size();
        counter++;
     } 
   maxGain = computeInfo(tid, hipart, whereMax, ret, k, maxGain);
     if (maxGain == -1){
        fprintf(stderr,"\n Partition %d is refined with whereMax partition %d ***** ", hipart, whereMax);
          break;
      }
    
     if(ret != 2)
     readBufMap[tid].erase(readBufMap[tid].begin(), it);
    }
     if(ret == 1) {
        writePartInfo(tid, hipart, whereMax);
   // pthread_mutex_unlock(&locks[tid]);
     }
     pthread_barrier_wait(&(barRefinepart));
     //   ComputeBECut(tid, readBufMap[tid]);
     //   ComputeBECut(tid, readBufMap[tid]);
}
//--------------------------------------------------
void Partitioner::writePartInfo(const unsigned tid, const unsigned hipart, const unsigned whereMax) {
     for(std::map<unsigned, unsigned>::const_iterator it = markMax[hipart].begin(); it!=markMax[hipart].end(); ++it){
         unsigned vtx1 = it->first;
         unsigned vtx2 = it->second;
         gWhere.at(vtx1) = whereMax;
         gWhere.at(vtx2) = hipart;
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx1, whereMax, hipart);
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx2, hipart, whereMax);
          changeWhere(tid, hipart, whereMax, gWhere, vtx1, vtx2);
        }
        bndIndMap[hipart].clear();
        //bndIndMap[whereMax].clear();
        totalPECuts[hipart] = 0;
       // totalPECuts[whereMax] = 0;
        markMax[hipart].clear();
        markMin[hipart].clear();
        dTable[hipart].clear();
       // dTable[whereMax].clear();
       // if(getWrittenToDisk())
	 //  refineMap->clear();
}

//--------------------------------------------------
void Partitioner::cWrite(const unsigned tid, unsigned noItems, InMemoryConstIterator end) {

  unsigned buffer = tid % nCols;
 // fprintf(stderr, "\nThread %d cWrite to partition %d\n", tid, buffer); 

  infinimem_cwrite_times[tid] -= getTimer();
    pthread_mutex_lock(&locks[buffer]);
  cWriteToInfinimem(buffer, totalCombined[tid], noItems, readBufMap[tid].begin(), end);
    pthread_mutex_unlock(&locks[buffer]);
    infinimem_cwrite_times[tid] += getTimer();
  //  fprintf(stderr,"\n*********TID: %d Total Combined : %d \n\n", tid, totalCombined[buffer]); 
}


//--------------------------------------------------
void Partitioner::cWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end) {
  RecordType* records = new RecordType[noItems]; 
  unsigned ct = 0;

  for (InMemoryConstIterator it = begin; it != end; ++it) {
     records[ct].set_rank(it->first);
  //    fprintf(stderr,"\n BWTI- TID: %d, startKey: %d, Key: %d\t, Values: ", buffer, startKey, it->first); 

    for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
      records[ct].add_nbrs(*vit);
    //  fprintf(stderr,"%d\t", *vit); 
      }
      ++ct;
      totalCombined[buffer]++;
  }

  assert(ct == noItems);
  cio->file_set_batch(buffer, startKey, noItems, records);
 delete[] records;
}

//--------------------------------------------------
void Partitioner::changeWhere(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gwhere, const unsigned maxVtx, const unsigned minVtx ) {
     fprintf(stderr,"\nPerm Changing WHERE for %d and %d\n", maxVtx, minVtx);
     where[hipart].at(maxVtx) = gwhere[maxVtx];
     where[whereMax].at(maxVtx) = gwhere[maxVtx];
     where[whereMax].at(minVtx) = gwhere[minVtx];
     where[hipart].at(minVtx) = gwhere[minVtx];
}

  static thread_local std::ofstream ofile;
 // static thread_local double stime;
//--------------------------------------------------
void Partitioner::printParts(const unsigned tid, std::string fileName) {
//       std::cout<<std::endl<<"Partition "<< tid <<std::endl;
//  std::string outputPrefix = "testing";
  ofile.open(fileName);
  assert(ofile.is_open());
 // stime = 0.0;
  for(unsigned i = 0; i <= nVtces; ++i){
     if(gWhere[i] != -1 && gWhere[i] == tid){
     //if(where[tid][i] != -1){// && where[tid][i] == tid){
//       std::cout<<"\t"<<i << "\t" << gWhere[i]<< std::endl;
   //    stime -= getTimer();
       ofile<<i << "\t" << gWhere[i]<< std::endl;
     //  stime += getTimer();
     }
  }
  ofile.close();
}

//--------------------------------------------------
void Partitioner::addDVals(const unsigned tid) {

   for (auto it = readBufMap[tid].begin(); it != readBufMap[tid].end(); ++it) {

        dTable[tid][it->first] = it->second.size();
   }
//   fprintf(stderr,"\nDtable tid %d size %d " , tid, dTable[tid].size());
}

//--------------------------------------------------
int Partitioner::checkPIDStarted(const unsigned tid, const unsigned hipart, const unsigned whereMax) {
 int ret = -1;
 pthread_mutex_lock(&locks[tid]);
 auto it_hi = pIdStarted.find(hipart); 
 auto it_wh = pIdStarted.find(whereMax); 
 if (it_hi != pIdStarted.end() || it_wh != pIdStarted.end()){   //that means key is present
   unsigned key1 = it_hi->first;
   unsigned val1 = it_hi->second;
   unsigned key2 = it_wh->first;
   unsigned val2 = it_wh->second;
    if(key1 == hipart && key2 == whereMax){ //both are present as keys with diff tids
     ret = 2;
    }
   else if (key1 == whereMax && key2 == hipart){
     ret = 2;
    }
    else if ((key1 == hipart && val1 == whereMax) || (key2 == whereMax && val2 == hipart)){
     fprintf(stderr,"\n TID %d partitioner %d is already present to refine with other thread  \n", tid, hipart);
     ret = 0;
    }
   else
      ret = 0;  //key present with a diff value
  }
 else{
     fprintf(stderr,"\n Emplacing %d %d ", hipart, whereMax);
     pIdStarted.emplace(hipart, whereMax);
     ret = 1;     // this will compute gain
 }
 
 pthread_mutex_unlock(&locks[tid]);
 return ret;
}
