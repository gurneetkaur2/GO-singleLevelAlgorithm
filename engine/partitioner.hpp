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
void Partitioner::initg(unsigned nVertices, unsigned nThreads, unsigned bSize, unsigned kItems, unsigned nParts)
{
  //nBuffers = pow(buffers, 2); // number of buffers is square of number of threads
  //nBuffers = nMappers * nReducers;

  nVtces = nVertices;
  nRows = nThreads;
  nCols = nParts; //TODO: should be two by default
  writtenToDisk = false;
  batchSize = bSize;
  kBItems = kItems;
//  nparts = nParts;
 // nvertices = nVertices;

  totalPECuts = new IdType[nCols];
  totalKeysInFile = new IdType[nCols];
  totalCombined = new IdType[nCols];
  readNext = new IdType[nCols];
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nCols];
  //nCuts = new IdType[nCols];
//  nEdges = new IdType[nRows * nCols];
    where = new std::vector<unsigned>[nCols];
//  partitionBndInd = new std::vector<unsigned>[nCols];
//  partitionBndPtr = new std::vector<unsigned>[nCols];
  
  outBufMap = new InMemoryContainer[nCols];
  readBufMap = new InMemoryContainer[nCols];
  refineMap = new InMemoryContainer[nCols];
  bndIndMap = new LookUpTable[nCols];
  lookUpTable = new LookUpTable[nCols];
  fetchBatchIds = new std::set<unsigned>[nCols];
  readNextInBatch = new std::vector<unsigned long long>[nCols];
  batchesCompleted = new std::vector<bool>[nCols];
  keysPerBatch = new std::vector<unsigned>[nCols];
  
  
  for (unsigned i=0; i<nRows * nCols; ++i){ 
    nItems[i] = 0;
    //nCuts[i] = 0;
    //nEdges[i] = 0;
    }

  //where = (IdType *) malloc((nVtces + 1) * sizeof(IdType));    
  for(unsigned i=0; i<nCols; ++i) {  
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    locks.push_back(mutex);
    
//    totalPECuts[i] = 0;
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

   for (unsigned i = 0; i < nCols; i++)
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
   for (unsigned i = 0; i < nCols; i++){
//        partitionBndInd[i].clear();
//        partitionBndPtr[i].clear();
  }

  delete cio;
   for (unsigned i = 0; i < nCols; i++){
    refineMap[i].clear();
    bndIndMap[i].clear();
  }
  delete[] refineMap;
  delete[] readNext;
  delete[] totalPECuts;
  delete[] totalCombined;
  delete[] bndIndMap;
//  delete[] partitionBndInd;
//  delete[] partitionBndPtr;
  //delete[] nReadKeys;
}

//-------------------------------------------------
void Partitioner::releaseReadPartStructures()
{
  delete io;
  fetchBatchIds->clear();
  batchesCompleted->clear();
  lookUpTable->clear();
  keysPerBatch->clear();
   for (unsigned i = 0; i < nCols; i++){
       readBufMap[i].clear();
       readNextInBatch[i].clear();
  }


  delete[] readBufMap;
  delete[] readNextInBatch;
  delete[] fetchBatchIds;
  delete[] batchesCompleted;
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
void Partitioner::writeBuf(const unsigned tid, const unsigned to, const unsigned from){

  unsigned buffer = hashKey(to) % nCols; 
//  unsigned buffer = tid * nCols + bufferId;  
//  std::vector<unsigned> where (nVtces,0);  //TODO: Ideally it should be batchsize here
//  std::vector<unsigned> whereDst (nVtces);
//   fprintf(stderr,"\nTID: %d assigning TO: %d Where: %d", tid, to, buffer); 
    if(where[tid].at(to) == -1){
      where[tid].at(to) = buffer;
  }

  // if the vertex assigned a partition is not added already then add it
//  if (std::find(whereDst.begin(), whereDst.end(), to) == whereDst.end()){
//     fprintf(stderr,"\nAdding to: %d in whereDst\n", to); 
//      whereDst.push_back(to);  //record all the vertices which have been assigned a partition avoiding duplicates
//   }

  //fprintf(stderr,"\nwhere[%d]: %d, where[%d]: %d, bufferID: %d\n\n", to, where[to], from, where[from], buffer);

  if (outBufMap[buffer].size() >= batchSize) {
    fprintf(stderr,"\noutbufmap buffer %d full with %d records noItems:%d\n", buffer, outBufMap[buffer].size(), nItems[buffer]);
     
  //  ComputeBECut(tid, buffer, where, bndind, bndptr, outBufMap[buffer]);
//    sCopy(tid, bndind, bndptr, partitionBndInd[tid], partitionBndPtr[tid]);
    pthread_mutex_lock(&locks[buffer]);
   
   // gCopy(tid, gWhere[buffer], where, whereDst);
 
 //   totalPECuts[buffer] += nCuts[buffer];
    writeToInfinimem(buffer, totalKeysInFile[buffer], outBufMap[buffer].size(), outBufMap[buffer]);
    totalKeysInFile[buffer] += nItems[buffer];
    pthread_mutex_unlock(&locks[buffer]);

//    setNum(tid, nVtces, &where, 1);
//     fprintf(stderr,"\ntotalPECuts[%d]: %d\n\n", buffer, totalPECuts[buffer]);
    outBufMap[buffer].clear();
    nItems[buffer] = 0;
 //   nCuts[buffer] = 0;
    //nEdges[buffer] = 0;
    writtenToDisk = true;
  }

  performWrite(tid, buffer, to, from);
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
void Partitioner::performWrite(const unsigned tid, const unsigned buffer, const unsigned to, const unsigned from) {
  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it_to = outBufMap[buffer].find(to); 
  IdType edge;
//  InMemoryContainerIterator it_dst = std::find(outBufMap[buffer].begin(), outBufMap[buffer].end(), e); 

  if(it_to != outBufMap[buffer].end()){
     combine(to, it_to->second, nbrs);
     fprintf(stderr,"\nCombining for key: %d, value: %d", to, from);
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
//TODO: Race condition when boht the threads write to same buffer
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
     fprintf(stderr,"\nCT: %d noItems: %d InMemMap Size: %d\n", ct, noItems, inMemMap.size());
}
  assert(ct == noItems);
//      fprintf(stderr,"\nWRITING TO INFINIMEM start: %d\n", startKey); 
  io->file_set_batch(buffer, startKey, noItems, records);
      fprintf(stderr,"\nWRITTEN %d records\n", noItems); 

  delete[] records;
}

//--------------------------------------------------
void Partitioner::bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end) {
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
/*  RecordType* record = new RecordType[noItems]; 
      fprintf(stderr,"\nBetter READ- tid: %d, StartKey: %d ", buffer, startKey); 
  cio->file_get_batch(buffer, startKey, noItems, record);
    for (unsigned i = 0; i < noItems; i++) {
      fprintf(stderr,"\nBetter READ- TID: %d, Key: %d\t Values: ", buffer, record[i].rank()); 

      for (unsigned k = 0; k < record[i].nbrs_size(); k++){

      fprintf(stderr,"%d\t", record[i].nbrs(k)); 
     }
    }
  delete[] record;
*/
}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

  fprintf(stderr,"\nFlushing buffer residues tid: %d\n", tid);
     
  if(tid >= nCols) 
    return;
  
//  if(nCols == 1) {
    writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
    outBufMap[tid].clear();
    totalKeysInFile[tid] += nItems[tid];
    nItems[tid] = 0;
//  }
   
/*  else {
    unsigned b1 = 0, b2 = 0;
    InMemoryContainerIterator b2Iter, b2End;
    unsigned long long b2Merged = 0;
    bool findB1 = true, findB2 = true;
    unsigned i = 0;
    if(i == nCols - 1) {
      fprintf(stderr,"\nCheck1 TID: %d, i: %d", tid, i);
      if(findB1 && findB2) {
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
        outBufMap[i].clear();
        totalKeysInFile[tid] += nItems[i];
        nItems[i] = 0;
      }
    }
    while(i < nCols - 1) {
      fprintf(stderr,"\nCheck2 TID: %d, i: %d", tid, i);
      if(findB1) {
        b1 = i % nCols;
        findB1 = false;
        ++i;
      }
      if(findB2) {
        b2 = i % nCols;
        b2Iter = outBufMap[b2].begin();
        b2End = outBufMap[b2].end();
        b2Merged = 0;
        findB2 = false;
        ++i;
      }
      fprintf(stderr,"\nC1 Merging Buffer: %d, buffer: %d", b1, b2);
      b2Merged += merge(outBufMap[b1], b1, tid, b2Iter, b2End);

      if(outBufMap[b1].size() == 0) {
        findB1 = true;
      }
      if(b2Iter == b2End) {
        findB2 = true;
      }
    }

    if(i == nCols - 1) {
      fprintf(stderr,"\nCheck3 TID: %d, i: %d", tid, i);
      if(findB1 && findB2) {
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
        outBufMap[i].clear();
        totalKeysInFile[tid] += nItems[i];
        nItems[i] = 0;
      }
      else if(findB1 || findB2) {
        if(findB1) {
          b1 = i % nCols;
          findB1 = false;
          ++i;
        } else {
          b2 = i % nCols;
          b2Iter = outBufMap[b2].begin();
          b2End = outBufMap[b2].end();
          b2Merged = 0;
          findB2 = false;
          ++i;
        }  

      fprintf(stderr,"\nC2 Merging Buffer: %d, buffer: %d", b1, b2);
        b2Merged += merge(outBufMap[b1], b1, tid, b2Iter, b2End);
        if(outBufMap[b1].size() == 0) {
          if(b2Iter != b2End) {
            bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
            outBufMap[b2].clear();
            totalKeysInFile[tid] += (nItems[b2] - b2Merged);
            nItems[b2] = 0;
          }
        } 
        if(b2Iter == b2End) {
          if(outBufMap[b1].size() != 0) {
            writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
            outBufMap[b1].clear();
            totalKeysInFile[tid] += nItems[b1];
            nItems[b1] = 0;
          }
        }
      } else
          assert(false);
    } else if(i == nCols) {
      if(outBufMap[b1].size() > 0) {
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
        outBufMap[b1].clear();
        totalKeysInFile[tid] += nItems[b1];
        nItems[b1] = 0;
      } else {  
        bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
        outBufMap[b2].clear();
        totalKeysInFile[tid] += (nItems[b2] - b2Merged);
        nItems[b2] = 0;
      }
    } else
      assert(false);
  }*/
}

/*//--------------------------------------------------
unsigned long long Partitioner::merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end) {
  unsigned long long ct = 0;
  fprintf(stderr,"\nMerging to Buffer: %d, tid: %d", whichMap, tid);
  while(begin != end) {
    if(toMap.size() >= batchSize) {
      writeToInfinimem(tid, totalKeysInFile[tid], toMap.size(), toMap);
      toMap.clear();
      totalKeysInFile[tid] += nItems[whichMap];
      nItems[whichMap] = 0; 
      return ct;
    }
//  fprintf(stderr,"\nMerging- adding records");
    InMemoryContainerIterator it = toMap.find(begin->first);
//  fprintf(stderr,"\nMerging- Check 1");
    if(it != toMap.end()){
//  fprintf(stderr,"\nMerging- Check 2");
      combine(it->first, it->second, begin->second);
      }
    else {
//  fprintf(stderr,"\nMerging- Check 3");
      toMap.emplace(begin->first, begin->second);
      ++nItems[whichMap];
//      ++nEdges[whichMap];
    }
//  fprintf(stderr,"\nMerging- Check 4");
    ++ct;
    ++begin;
  }
  return ct;
}
*/
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
}


//--------------------------------------------------
bool Partitioner::read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted) {

//  fprintf(stderr,"\nInside Partitioner::read \n");
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

     // fprintf(stderr,"\nREAD- TID: %d, startKey: %d ", tid, readNextInBatch[batch]); 
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

  delete[] records;
  return ret;
}

//--------------------------------------------------
bool Partitioner::read(const unsigned tid) {
  return read(tid, readBufMap[tid], keysPerBatch[tid], lookUpTable[tid], fetchBatchIds[tid], readNextInBatch[tid], batchesCompleted[tid]);
}

//--------------------------------------------------
void Partitioner::initiateInMemoryRefine(unsigned tid) {

     refineMap[tid].insert(outBufMap[tid].begin(), outBufMap[tid].end());

}

//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid) {
     ComputeBECut(tid, gWhere, bndIndMap[tid], refineMap[tid]);
}
//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid, const std::vector<unsigned>& where, LookUpTable& bndind, const InMemoryContainer& inMemMap) {

//TODO bndind and bndptr needs to be locked when used by more than one thread
// do not need first flag now
//  pthread_mutex_lock(&locks[tid]);
//  gCopy(tid);
//  pthread_mutex_unlock(&locks[tid]);
    fprintf(stderr,"\nInside ComputeBECUT -------"); 
     IdType i, j, nbnd=0;
     unsigned pid = tid % nCols;
//     IdType *bndind, *bndptr;
     unsigned src;
     bool first = 1;
     std::vector<unsigned> bndvert;     

  //  fprintf(stderr,"\nBefore looping ComputeBECUT -------"); 
  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {
      src = it->first;
      for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
       //   if (std::find(bndvert.begin(), bndvert.end(), *vit) == bndvert.end()){
              bndvert.push_back(*vit);
         // }
       }
   //       const std::vector<unsigned>& nbrs = i->second;
        //  unsigned max = *max_element(std::begin(nbrs), std::end(nbrs));
  // }
    //     fprintf(stderr,"ComputeCut- bndvert size: %d\n", bndvert.size());
     for(unsigned i=0; i<bndvert.size(); i++){
      //   fprintf(stderr,"ComputeCut- i: %d\n", i);
         IdType dst = bndvert[i];
        //   fprintf(stderr,"\nTID: %d, src :%d, dst: %d", tid, src, dst); 
        //   fprintf(stderr,"\n"); 
 //       for (j=xadj[src]; j<xadj[src+1]; j++){

 //               IdType dst = adjncy[j];
//            fprintf(stderr,"PECUT tid: %d, src: %d, j: %d, xadj[src+1]: %d, adjncy[j]: %d\n", tid, src, j,xadj[src+1], dst);
          //       fprintf(stderr,"where[%d]: %d, where[%d]: %d\n", src, where[src], dst ,where[dst]);
         //    if(!(adjncy[j] > max)){
                 if( where[dst] != -1 && where[src] != where[dst] ) {
               //  if( where[src] != where[dst] ) {
            //        fprintf(stderr,"Edge CUT tid: %d\n", tid);
                    nbnd++;
    		    totalPECuts[tid]++;
                   // dst = adjncy[j];
          //   fprintf(stderr,"\nBoundary Vertices tid: %d, dst: %d\n", tid, dst);
          // store the boundary vertices along with their actual location
              bndind[dst].push_back(where[dst]); // TODO: should this be where[dst]?
     fprintf(stderr,"\nSize bndind[%d]: %d\n", dst, bndind[dst].size());
              //      first = 0;
              //      }
                 }
             // }
          }
 /*         if(xadj[src] == xadj[src+1]){
          // store the boundary vertices
             fprintf(stderr,"\nBoundary source: %d, cut: %d\n", src, cut);
             BNDInsert(cut, bndind, bndptr, src);
          }
  */	bndvert.clear();
     }
//     fprintf(stderr,"\nNbnd: %d\n", nbnd);
     fprintf(stderr,"\ntPECuts[%d]: %d, MapSize: %d\n", tid, totalPECuts[tid], refineMap[tid].size());

}

//--------------------------------------------------
void Partitioner::setNum(const unsigned tid, std::vector<unsigned>& where, unsigned num){
     for(unsigned i=0; i<nVtces; ++i){
  fprintf(stderr,"\nBreaking 3\n"); 
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
             if(first){
//         fprintf(stderr,"\nwhere[%d][%d]: %d,", i, j, where[i][j]);
                gWhere[j] = where[i][j];
             }
             else {
                if(where[i][j] != -1){
                   gWhere[j] = where[i][j];
               }
            }
   //      fprintf(stderr,"\nGWHERE[%d]: %d", j, gWhere[j]);
         }
         first = 0;
     }
}

//--------------------------------------------------
void Partitioner::sCopy(const unsigned tid, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, std::vector<unsigned>& partitionBndInd, std::vector<unsigned>& partitionBndPtr){
      for(unsigned i=0; i<nVtces; ++i){
         partitionBndInd.at(i) = bndind[i];
         partitionBndPtr.at(i) = bndptr[i];
     }
}

//--------------------------------------------------
unsigned Partitioner::countTotalPECut(const unsigned tid) {
      totalCuts = 0;
      //for(unsigned i=0; i<nCols; i++){
      for(auto i=fetchPIds.begin(); i != fetchPIds.end(); ++i){
          totalCuts += totalPECuts[*i];
      }
   return (totalCuts/2);
}

//--------------------------------------------------
unsigned Partitioner::maxPECut(const unsigned tid) {
      unsigned low = 0, hipart = -1;
    fprintf(stderr,"\nMaxPECut fetchpids size %d", fetchPIds.size());
      //for(unsigned i=0; i<nCols; i++){
      for(auto i=fetchPIds.begin(); i != fetchPIds.end(); ++i){
          if(low < totalPECuts[*i]){
             low = totalPECuts[*i];
             hipart = *i;
          }
      }
   fprintf(stderr,"\n Partition: %d has max cuts: %d", hipart, totalPECuts[hipart]);
   return hipart;
}

//--------------------------------------------------
void Partitioner::refinePart(const unsigned tid, const unsigned hipart, unsigned tCuts) {
     refinePart(tid, hipart, tCuts, gWhere);

}

//--------------------------------------------------
void Partitioner::refinePart(const unsigned tid, const unsigned hipart, unsigned tCuts, std::vector<unsigned>& gWhere) {
//Check if moving the vertex with max bnd value reduce the num of cuts in the parition and how much does it impact in the current partition also check the total edgecuts
  unsigned invalidmoves = 0;
  // in the partition with highest edge cuts, find the bnd vertex with max occurrence
  unsigned maxvtx = -1;
//  unsigned minvtx = -1;
//  unsigned tCuts = countTotalPECut(tid);
  fprintf(stderr,"\n----TID: %d, Total Edge cuts : %d\n", tid, tCuts);
  unsigned pCut = getTotalPECuts(hipart);  //total partition edge cuts
  std::vector<unsigned> boundpart {nCols};  // to keep track of partitions
  
  for (InMemoryConstIterator it = bndIndMap[hipart].begin(); it != bndIndMap[hipart].end(); ++it) {
      	boundpart.push_back(it->first);
  }

 // Loop through all the boundary vertices
 while(boundpart.size() > 0){
    maxvtx = maxBound(tid, bndIndMap[hipart]);
      std::map<unsigned, unsigned>::iterator it_max = markMax.find(maxvtx);
 //   unsigned extmax = bndIndMap[hipart][maxvtx].size();
    fprintf(stderr," in partition %d \n", hipart);
      if(maxvtx != -1 && it_max == markMax.end()){
         unsigned whereMax = gWhere[maxvtx];  //where is the vertex with max boundary positions in this selected partition
//    unsigned intmax = refineMap[whereMax][maxvtx].size();
//    unsigned dmax = extmax - intmax;
         unsigned epCut = getTotalPECuts(whereMax); // get total edgecuts of the partition where that vertex is
         fprintf(stderr,"\nREFINEPART Max vertex %d is in partition %d with MaxECuts %d ", maxvtx, whereMax, epCut);
     unsigned oldpreduction = -2;// set it to random negative value
   // Check the edgecuts in the partition where this Maxvertex is currently
   // Get the vtx with min edgecuts in this partition
        unsigned minvtx = minBound(tid, bndIndMap[whereMax], hipart);
    fprintf(stderr," in partition %d \n", whereMax);
      std::map<unsigned, unsigned>::iterator it_min = markMin.find(minvtx);
      if(minvtx != -1 && it_min == markMin.end() && whereMax != hipart){
  //      unsigned extmin = bndIndMap[whereMax][minvtx].size();
  //  unsigned intmin = refineMap[hipart][minvtx].size();
  //  unsigned dmin = extmin - intmin;
// move the above found vertex from the partition where it exists to the current partition(hipart) with max edgecuts and remove it from bndindMap
     gWhere[maxvtx] = hipart;
     gWhere[minvtx] = whereMax;
     fprintf(stderr,"\nREFINEPART AFTER  where[maxvtx]: %d ", gWhere[maxvtx]);
   fprintf(stderr,"\nREFINEPART AFTER  where[minvtx]: %d ", gWhere[minvtx]);
  //Iterate through the partition and recompute edge cuts after moving the vertex 
     refineInit(hipart); refineInit(whereMax); bndIndMap[whereMax].clear(); bndIndMap[hipart].clear();
     cread(hipart); cread(whereMax);
   // Calculate the total and partition edgecuts after iterating through combined records
     unsigned oldTCuts = tCuts;
     tCuts = countTotalPECut(tid); 
     fprintf(stderr,"\nRefining partition %d next iteration and found max cuts: %d\n", hipart, tCuts);
     unsigned newPCut = getTotalPECuts(hipart); //new total partition cuts
     unsigned newEPCut = getTotalPECuts(whereMax);  //new total other partition cuts
     unsigned diff = oldTCuts - tCuts;   //difference in total number of edgecuts
     unsigned preduction = (pCut - newPCut) + (epCut - newEPCut) - 2;
    // unsigned preduction = extmax + (epCut - newEPCut) - 2;
     if(preduction > oldpreduction){  //if edgecuts decresed for current partition
        markMax[maxvtx] = preduction;// mark the max and min vtces being swapped
        markMin[minvtx] = preduction;
       if(oldTCuts < tCuts && invalidmoves < 50){   //total edgecuts should alse decrease
          fprintf(stderr,"\nMoving %d to %d from %d ", maxvtx, hipart, whereMax);
          fprintf(stderr,"\nMoving %d to %d from %d ", minvtx, whereMax, hipart);
          changeWhere(tid, hipart, whereMax, gWhere, maxvtx, minvtx);
          setTotalPECuts(hipart);
          pCut = newPCut;
          epCut = newEPCut;
          oldpreduction = preduction;
          invalidmoves = 0;  //reset
      }
    }
     else{
          fprintf(stderr,"\nUn-Moving %d to %d from %d ", maxvtx, whereMax, hipart);
          fprintf(stderr,"\nUn-Moving %d to %d from %d ", minvtx, hipart, whereMax);
         gWhere[maxvtx] = whereMax;
         gWhere[minvtx] = hipart;
         invalidmoves++;   //if the total edgecut did not decrease from previous iteration
         if(invalidmoves >= 50)
		break;
     }
   }
     else{
          fprintf(stderr,"\nNo Swap vertex from other partition ***** ");
          bndIndMap[hipart].erase(maxvtx);
        }
   }
   else{
        fprintf(stderr,"\nPartition %d is refined ***** ", hipart);
        bndIndMap[hipart].erase(maxvtx);
        boundpart.clear();
      }
      // Delete the boundary vertices processed and update again
      deletebndvert(tid, hipart, markMax);
     fprintf(stderr,"\nBndInd has %d elements after iteration", bndIndMap[hipart].size());
     if(boundpart.size() > 0){
      boundpart.erase(std::remove(boundpart.begin(), boundpart.end(), maxvtx), boundpart.end());
      fprintf(stderr,"\nremoved %d from Boundpart size: %d", maxvtx, boundpart.size());
    }
 }
}

//--------------------------------------------------
void Partitioner::deletebndvert(const unsigned tid, const unsigned hipart, std::map<unsigned, unsigned>& markMax ) {
  for(auto it = markMax.begin(); it!= markMax.end(); ++it){
     bndIndMap[hipart].erase(it->first);
     fprintf(stderr,"\nDeleting %d from bndIndMap size %d", it->first, bndIndMap[hipart].size());
  }

}

//--------------------------------------------------
unsigned Partitioner::maxBound(const unsigned tid, LookUpTable& bndind ) {
//TODO: Add a vertex which stores all the vertices with same number of edge cuts and return that vertex
unsigned currentMax = 0;
unsigned vtx_max = -1;
for(auto it = bndind.cbegin(); it != bndind.cend(); ++it ) {
    if (it ->second.size() > currentMax) {
        vtx_max = it->first;
        currentMax = it->second.size();
    }
}
fprintf(stderr, "\nVertex %d occurs %d times as boundary vertex", vtx_max, currentMax);

return vtx_max;
}

//--------------------------------------------------
unsigned Partitioner::minBound(const unsigned tid, LookUpTable& bndind, const unsigned hipart ) {
//TODO: Add a vertex which stores all the vertices with same number of edge cuts and return that vertex
unsigned currentMin = 2147483647;  //set it to some highest value initially
unsigned vtx_min = -1;
for(auto it = bndind.cbegin(); it != bndind.cend(); ++it ) {
    //const std::vector<unsigned>& lookwhere = it->second;
   //If the the vtx_min belongs to the partition being refined (hipart)
 //   if (std::find(lookwhere.begin(), lookwhere.end(), hipart) != lookwhere.end()){
     // then consider it for swap
      if (it ->second.size() < currentMin) {
          vtx_min = it->first;
          currentMin = it->second.size();
      }
  // }
}
fprintf(stderr, "\nMin Vertex %d occurs %d times as boundary vertex\n", vtx_min, currentMin);

return vtx_min;
}

//--------------------------------------------------
void Partitioner::cWrite(const unsigned tid, unsigned noItems, InMemoryConstIterator end) {

  unsigned buffer = tid % nCols;
 // fprintf(stderr, "\nThread %d cWrite to partition %d\n", tid, buffer); 
  
    pthread_mutex_lock(&locks[buffer]);
  bWriteToInfinimem(buffer, totalCombined[tid], noItems, readBufMap[tid].begin(), end);
    pthread_mutex_unlock(&locks[buffer]);
     
}


//--------------------------------------------------
/*bool Partitioner::refine(const unsigned tid) {
  return refine(tid, refineMap[tid]);
}
*/
//--------------------------------------------------
bool Partitioner::refine(const unsigned tid) {

  unsigned partition = tid;

 
    unsigned partbound = min(totalCombined[tid]-readNext[tid], kBItems); 
  RecordType* parts = new RecordType[partbound];
  fprintf(stderr,"\nREFINE tid: %d, totalCombined: %d, readNext[tid]: %d, keys to read: %d \n", tid, totalCombined[tid], readNext[tid], partbound);
 //for(unsigned ckey = readNextInBatch[partition]; ckey < totalCombined[partition]; ckey += batchSize){
    if (partbound > 0 && readNext[tid] < totalCombined[tid])
      cio->file_get_batch(tid, readNext[tid], partbound, parts); 

  //  fprintf(stderr,"\nREFINE- tid: %d  BATCHSIZE: %d, Map size: %d, readNext: %d \n", tid, batchSize, refineMap[tid].size(), readNext[tid]);
    for (unsigned i = 0; i < partbound; i++) {
      refineMap[tid][parts[i].rank()];
//      fprintf(stderr,"\nREFINE- TID: %d, Key: %d\t Values: ", tid, parts[i].rank()); 

      for (unsigned k = 0; k < parts[i].nbrs_size(); k++){
        refineMap[tid][parts[i].rank()].push_back(parts[i].nbrs(k));

  //    fprintf(stderr,"%d\t", parts[i].nbrs(k)); 
     }
    }
    
//    fprintf(stderr,"\nREFINE - tid: %d, RefineMap size: %d", tid, refineMap[tid].size());
    readNext[tid] += partbound;
//    fprintf(stderr,"\nREFINE update-  tid: %d Total Combined: %d, readNext: %d \n", tid, totalCombined[tid], readNext[tid]);

  bool ret = false;
    if (readNext[tid] < totalCombined[tid]){
        ret = true;
  //      fprintf(stderr, "\nREFINE - still reading - tid: %d, readNext: %d, totalCombined: %d, ret: %d \n", tid, readNext[tid], totalCombined[tid], ret);
    }

  delete[] parts;
  return ret;
}

//--------------------------------------------------
void Partitioner::refineInit(const unsigned tid) {
  
        readNext[tid] = 0;
        totalPECuts[tid] = 0;

 // fprintf(stderr,"\nREFINEInit tid: %d, readNextInBatch[tid]: %d \n", tid, readNextInBatch[tid]);
//  std::map<unsigned, unsigned> bndind;
}

//--------------------------------------------------
void Partitioner::cread(const unsigned tid) {
 
while(true) {
//  fprintf(stderr, "\nDoRefine: Calling Read\n");
    bool execLoop;
    if(!getWrittenToDisk()){
      execLoop = 0;
    } 
   else{
       execLoop = refine(tid);
    }
    fprintf(stderr, "\nExecloop is %d", execLoop);

    if(execLoop == false) {
    //  for(InMemoryConstIterator it = refineMap[tid].begin(); it != refineMap[tid].end(); ++it){
    //     fprintf(stderr,"\nCREAD- tid: %d, Computing edgecuts with Map Size %d", tid, refineMap[tid].size());
         ComputeBECut(tid);
       if(getWrittenToDisk()){
         refineMap[tid].clear(); //erase(refineMap[tid].begin(), it) erase when writing to disk
      }
       break;
    }

   // Read the kitems from infinimem and compute the edgecuts for that part
    ComputeBECut(tid);

    refineMap[tid].erase(refineMap[tid].begin(), refineMap[tid].end());
  }

}

//--------------------------------------------------
//void Partitioner::changeWhere(const unsigned tid, const unsigned hipart, const unsigned chVtx ) {

//     changeWhere(tid, hipart, gWhere, chVtx);
//}
  
//--------------------------------------------------
void Partitioner::changeWhere(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::vector<unsigned>& gwhere, const unsigned maxVtx, const unsigned minVtx ) {
     where[hipart].at(maxVtx) = gwhere[maxVtx];
     where[whereMax].at(minVtx) = gwhere[minVtx];
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
       std::cout<<"\t"<<i << "\t" << gWhere[i]<< std::endl;
   //    stime -= getTimer();
       ofile<<i << "\t" << gWhere[i]<< std::endl;
     //  stime += getTimer();
     }
  }
  ofile.close();
}
