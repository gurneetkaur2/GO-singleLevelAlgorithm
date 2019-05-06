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
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nCols];
  nCuts = new IdType[nCols];
//  nEdges = new IdType[nRows * nCols];
 // gWhere = new std::vector<unsigned long long>[nCols];
  partitionBndInd = new std::vector<unsigned>[nCols];
  partitionBndPtr = new std::vector<unsigned>[nCols];
  
  outBufMap = new InMemoryContainer[nCols];
  readBufMap = new InMemoryContainer[nCols];
  lookUpTable = new LookUpTable[nCols];
  fetchBatchIds = new std::set<unsigned>[nCols];
  readNextInBatch = new std::vector<unsigned long long>[nCols];
  batchesCompleted = new std::vector<bool>[nCols];
  keysPerBatch = new std::vector<unsigned>[nCols];
  
  
  for (unsigned i=0; i<nRows * nCols; ++i){ 
    nItems[i] = 0;
    nCuts[i] = 0;
    //nEdges[i] = 0;
    }

  //where = (IdType *) malloc((nVtces + 1) * sizeof(IdType));    
  for(unsigned i=0; i<nCols; ++i) {  
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    locks.push_back(mutex);
    
    totalPECuts[i] = 0;
    totalKeysInFile[i] = 0;
  }

  // setup FileIO
//#ifdef USE_ONE_PHASE_IO
  //io = new OnePhaseFileIO<RecordType>("/tmp/gkaur007/mrdata/", nCols, 0/*UNUSED*/);
//#else
  io = new TwoPhaseFileIO<RecordType>("/tmp/gkaur007/mrdata/", nCols, 0/*UNUSED*/);
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
  delete[] nCuts;
//  delete[] nEdges;
  delete[] outBufMap;
//  delete[] bndind;
//  delete[] bndptr;
}
//-------------------------------------------------
void Partitioner::shutdown()
{
   for (unsigned i = 0; i < nCols; i++){
        partitionBndInd[i].clear();
        partitionBndPtr[i].clear();
  }

  delete io;
  delete[] totalPECuts;
  delete[] totalKeysInFile;
  delete[] partitionBndInd;
  delete[] partitionBndPtr;
  //delete[] nReadKeys;
}

//--------------------------------------------------
void Partitioner::writeInit(const unsigned tid) {
  unsigned j=0;
  for (unsigned i = 0; i<nVtces; ++i) {
       partitionBndInd[tid].push_back(NULL);
       partitionBndPtr[tid].push_back(-1);
  //TODO : not sure if I need to initialize the below loop
//    pthread_mutex_lock(&locks[tid]);
//       gWhere[i] = 1;
//    pthread_mutex_unlock(&locks[tid]);
  }

}


//-------------------------------------------------
unsigned Partitioner::writeBuf(const unsigned tid, const unsigned to, const unsigned from, unsigned k, unsigned* adjncy, std::vector<unsigned>& where, std::vector<unsigned>& whereDst, bool firstRep){

  unsigned xadjVal = k;
  unsigned buffer = tid % nCols; 
//  unsigned buffer = tid * nCols + bufferId;  
//  std::vector<unsigned> where (nVtces,0);  //TODO: Ideally it should be batchsize here
//  std::vector<unsigned> whereDst (nVtces);
  std::vector<unsigned> bndind (nVtces,0);
  std::vector<unsigned> bndptr (nVtces,-1); 

   where.at(to) = buffer;

  // if the vertex assigned a partition is not added already then add it
  if (std::find(whereDst.begin(), whereDst.end(), to) == whereDst.end()){
     fprintf(stderr,"\nAdding to: %d in whereDst\n", to); 
      whereDst.push_back(to);  //record all the vertices which have been assigned a part
   }

  fprintf(stderr,"\nwhere[%d]: %d, bufferID: %d\n\n", to, where[to], buffer);

  if (outBufMap[buffer].size() >= batchSize) {
    fprintf(stderr,"outbufmap buffer full with %d records\n", outBufMap[buffer].size());
     
//    ComputeBECut(tid, buffer, where, bndind, bndptr);
//    sCopy(tid, bndind, bndptr, partitionBndInd[tid], partitionBndPtr[tid]);
    pthread_mutex_lock(&locks[buffer]);
   
   // gCopy(tid, gWhere, where, whereDst);
 
    totalPECuts[buffer] += nCuts[buffer];
    writeToInfinimem(buffer, totalKeysInFile[buffer], outBufMap[buffer].size(), outBufMap[buffer]);
    totalKeysInFile[buffer] += nItems[buffer];
    pthread_mutex_unlock(&locks[buffer]);

//    setNum(tid, nVtces, &where, 1);
     fprintf(stderr,"\ntotalPECuts[%d]: %d\n\n", buffer, totalPECuts[buffer]);
    outBufMap[buffer].clear();
    nItems[buffer] = 0;
    nCuts[buffer] = 0;
    //nEdges[buffer] = 0;
    writtenToDisk = true;
  }

  xadjVal = performWrite(tid, buffer, to, from, k, adjncy, firstRep);
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
  return xadjVal;
}

//--------------------------------------------------
unsigned Partitioner::performWrite(const unsigned tid, const unsigned buffer, const unsigned to, const unsigned from, unsigned k, unsigned* adjncy, bool firstRep) {
  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it_to = outBufMap[buffer].find(to); 
  IdType edge;
//  InMemoryContainerIterator it_dst = std::find(outBufMap[buffer].begin(), outBufMap[buffer].end(), e); 

  if(it_to != outBufMap[buffer].end()){
     combine(to, it_to->second, nbrs);
     if(firstRep){
        edge = to;
        adjncy[k] = edge - 1;
    fprintf(stderr,"\nsrc: %d, adjncy: %d, edge-1: %d\n", to, adjncy[k], edge-1);
        k++;
     }
     edge = from;
     adjncy[k] = edge - 1;
     fprintf(stderr,"\nsrc: %d, adjncy: %d, edge-1: %d\n", from, adjncy[k], edge-1);
     k++;
//     outBufMap[buffer][to].push_back(from);
//     nEdges[buffer]++;
      }
  else {
    outBufMap[buffer].emplace(to, nbrs); 
    if(firstRep){
      edge = to;
      adjncy[k] = edge - 1;
    fprintf(stderr,"\nsrc: %d, adjncy: %d, edge-1: %d\n", to, adjncy[k], edge-1);
      k++;
     }
    edge = from;
    adjncy[k] = edge - 1;
    fprintf(stderr,"\nsrc: %d, adjncy: %d, edge-1: %d\n", from, adjncy[k], edge-1);
    k++;
    nItems[buffer]++;
  //  nEdges[buffer]++;
  }
 return k;
}

//--------------------------------------------------
void Partitioner::writeToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, const InMemoryContainer& inMemMap) {
  RecordType* records = new RecordType[noItems]; 
  unsigned ct = 0;

  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {
      records[ct].set_rank(it->first);
      for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
          records[ct].add_nbrs(*vit);
      }
   ++ct;

  }
    if (ct != noItems){
     fprintf(stderr,"\nCT: %d noItems: %d InMemMap Size: %d\n", ct, noItems, inMemMap.size());
}
  assert(ct == noItems);
  io->file_set_batch(buffer, startKey, noItems, records);

  delete[] records;
}

//--------------------------------------------------
void Partitioner::bWriteToInfinimem(const unsigned buffer, const IdType startKey, unsigned noItems, InMemoryConstIterator begin, InMemoryConstIterator end) {
/*  RecordType* records = new RecordType[noItems]; 
  unsigned ct = 0;

  for (InMemoryConstIterator it = begin; it != end; ++it) {
     records[ct].set_rank(it->first);

    for (std::vector<unsigned>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit)
      records[ct].add_nbrs(*vit);

      ++ct;
  }

  assert(ct == noItems);
  io->file_set_batch(buffer, startKey, noItems, records);

  delete[] records;
*/}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

  //fprintf(stderr,"\nFlushing buffer residues\n");
    /* 
  if(tid >= nCols) 
    return;

  if(nRows == 1) {
    ComputeBECut(tid, tid, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
    totalPECuts[tid] += nCuts[tid];
    writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
    outBufMap[tid].clear();
    totalKeysInFile[tid] += nItems[tid];
    nItems[tid] = 0;
    nCuts[tid] = 0;
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
        ComputeBECut(tid, i, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
        totalPECuts[tid] += nCuts[i];
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
        outBufMap[i].clear();
        totalKeysInFile[tid] += nItems[i];
        nItems[i] = 0;
        nCuts[i] = 0;
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
            ComputeBECut(tid, b2, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
            totalPECuts[tid] += (nCuts[b2] - b2Merged);
            bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
            outBufMap[b2].clear();
            totalKeysInFile[tid] += (nItems[b2] - b2Merged);
            nCuts[b2] = 0;
            nItems[b2] = 0;
          }
        } 
        if(b2Iter == b2End) {
          if(outBufMap[b1].size() != 0) {
            ComputeBECut(tid, b1, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
            totalPECuts[tid] += nCuts[b1];
            writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
            outBufMap[b1].clear();
            totalKeysInFile[tid] += nItems[b1];
            nItems[b1] = 0;
            nCuts[b1] = 0;
          }
        }
      } else
          assert(false);
    } else if(i == nRows) {
      if(outBufMap[b1].size() > 0) {
        ComputeBECut(tid, b1, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
        totalPECuts[tid] += nCuts[b1];
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
        outBufMap[b1].clear();
        totalKeysInFile[tid] += nItems[b1];
        nItems[b1] = 0;
        nCuts[b1] = 0;
      } else {  
        ComputeBECut(tid, b2, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
        totalPECuts[tid] += (nCuts[b2] - b2Merged);
        bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
        outBufMap[b2].clear();
        totalKeysInFile[tid] += (nItems[b2] - b2Merged);
        nItems[b2] = 0;
        nCuts[b2] = 0;
      }
    } else
      assert(false);
  }*/
}

//--------------------------------------------------
unsigned long long Partitioner::merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end) {
/*  unsigned long long ct = 0;
  while(begin != end) {
    if(toMap.size() >= batchSize) {
      ComputeBECut(tid, whichMap, gWhere, partitionBndInd[tid], partitionBndPtr[tid]);
      totalPECuts[tid] += nCuts[whichMap];
      writeToInfinimem(tid, totalKeysInFile[tid], toMap.size(), toMap);
      toMap.clear();
      totalKeysInFile[tid] += nItems[whichMap];
      nItems[whichMap] = 0; 
      nCuts[whichMap] = 0;
      return ct;
    }
    InMemoryContainerIterator it = toMap.find(begin->first);
    if(it != toMap.end())
      combine(it->first, it->second, begin->second);
    else {
      toMap.emplace(begin->first, begin->second);
      ++nItems[whichMap];
      ++nEdges[whichMap];
    }
    ++ct;
    ++begin;
  }
  return ct;
*/}

//--------------------------------------------------
void Partitioner::readInit(const unsigned tid) {
  unsigned j=0;
  for (unsigned long long i = 0; i <= totalKeysInFile[tid]; i+= batchSize) {
    readNextInBatch[tid].push_back(i); 
    fetchBatchIds[tid].insert(j++); 
    batchesCompleted[tid].push_back(false); 
    keysPerBatch[tid].push_back(kBItems); 
  }
}


//--------------------------------------------------
bool Partitioner::read(const unsigned tid, InMemoryContainer& readBufMap, std::vector<unsigned>& keysPerBatch, LookUpTable& lookUpTable, std::set<unsigned>& fetchBatchIds, std::vector<unsigned long long>& readNextInBatch, std::vector<bool>& batchesCompleted) {
/*
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

    if (keysPerBatch[batch] > 0 && readNextInBatch[batch] < batchBoundary)
      io->file_get_batch(tid, readNextInBatch[batch], keysPerBatch[batch], records); 

    for (unsigned i = 0; i < keysPerBatch[batch]; i++) {
      lookUpTable[records[i].rank()].push_back(batch);
      readBufMap[records[i].rank()];

      for (unsigned k = 0; k < records[i].nbrs_size(); k++)
        readBufMap[records[i].rank()].push_back(records[i].nbrs(k));
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
*/}

//--------------------------------------------------
bool Partitioner::read(const unsigned tid) {
  return read(tid, readBufMap[tid], keysPerBatch[tid], lookUpTable[tid], fetchBatchIds[tid], readNextInBatch[tid], batchesCompleted[tid]);
}

//--------------------------------------------------
InMemoryReductionState Partitioner::initiateInMemoryReduce(unsigned tid) {
  InMemoryReductionState state(nRows); 
/*  for(unsigned i=0; i<nRows; ++i) {
    state.begins[i] = outBufMap[tid + nCols * i].begin();
    state.ends[i] = outBufMap[tid + nCols * i].end();
  }
*/
  return state;
}

//--------------------------------------------------
bool Partitioner::getNextMinKey(InMemoryReductionState* state, InMemoryContainer* record) {
/*  std::vector<unsigned> minIds;
  unsigned minKey;
  bool found = false;

  for(unsigned i=0; i<nRows; ++i) {
    if(state->begins[i] == state->ends[i])
      continue;

    if(!found) {
      minKey = state->begins[i]->first;
      minIds.push_back(i);
      found = true;
    } else {
      if(state->begins[i]->first < minKey) {
        minKey = state->begins[i]->first;
        minIds.clear();
        minIds.push_back(i);
      } else if(state->begins[i]->first == minKey) {
        minIds.push_back(i);
      }
    }
  }

  if(!found)
    return false;

  std::vector<unsigned>& values = (*record)[minKey];
  for(std::vector<unsigned>::iterator it = minIds.begin(); it != minIds.end(); ++it) {
    for(std::vector<unsigned>::const_iterator vit = state->begins[*it]->second.begin(); vit != state->begins[*it]->second.end(); ++vit) 
      values.push_back(*vit);

    ++state->begins[*it];
  }
*/
  return true;
}


//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid, const unsigned buffer, std::vector<unsigned>& where, std::vector<unsigned>& bndind, std::vector<unsigned> bndptr) {

//TODO bndind and bndptr needs to be locked when used by more than one thread
// do not need first flag now
/*     IdType i, j, nbnd=0;
//     IdType *bndind, *bndptr;

     bool first = 1;
     std::vector<unsigned> bndvert;     

     for(IdType i=0; i<outBufMap[buffer].size(); i++) {
         IdType src = outBufMap[buffer][i];

         if (std::find(bndvert.begin(), bndvert.end(), src) == bndvert.end()){
             bndvert.push_back(src);
         }
   //       const std::vector<unsigned>& nbrs = i->second;
        //  unsigned max = *max_element(std::begin(nbrs), std::end(nbrs));
       }
     for(unsigned i=0; i<bndvert.size(); i++){
           IdType src = bndvert[i];
           fprintf(stderr,"\nTID: %d, src :%d, buffer: %d", tid, src, buffer); 
           fprintf(stderr,"\n"); 
        for (j=xadj[src]; j<xadj[src+1]; j++){

                IdType dst = adjncy[j];
            fprintf(stderr,"PECUT tid: %d, src: %d, j: %d, xadj[src+1]: %d, adjncy[j]: %d\n", tid, src, j,xadj[src+1], dst);
                 fprintf(stderr,"where[%d]: %d, where[%d]: %d\n", src, where[src], dst ,where[dst]);
         //    if(!(adjncy[j] > max)){
                 fprintf(stderr,"\nChecking if the where is same\n");
                 if( where[dst] != -1 && where[src] != where[dst] ) {
                    fprintf(stderr,"Edge CUT tid: %d\n", tid);
                    nbnd++;
    		    nCuts[buffer]++;
                   // dst = adjncy[j];
             fprintf(stderr,"\nBoundary Vertices tid: %d, dst: %d, nbnd: %d, bndptr: %d, first: %d\n", tid, dst, nbnd, bndptr[dst], first);
          // store the boundary vertices
                    BNDInsert(tid, nbnd, bndind, bndptr, dst, first);
                    first = 0;
                 }
             // }
          }
*/ /*         if(xadj[src] == xadj[src+1]){
          // store the boundary vertices
             fprintf(stderr,"\nBoundary source: %d, cut: %d\n", src, cut);
             BNDInsert(cut, bndind, bndptr, src);
          }
   *///  }
  	
//     fprintf(stderr,"\nNbnd: %d\n", nbnd);
//     fprintf(stderr,"\nnCuts[%d]: %d, MapSize: %d\n", buffer, nCuts[buffer], outBufMap[buffer].size());

}


//--------------------------------------------------
void Partitioner::BNDInsert(const unsigned tid, unsigned n, std::vector<unsigned>& bndind, std::vector<unsigned> bndptr, unsigned i, bool first){
     do {
//     if(first){ 
        assert(bndptr[i] == -1); 
//     }
     bndind[n] = i; 
     fprintf(stderr,"\nbndind[%d]: %d\n",n, i);
     bndptr[i] = (n)++;
   } while(0); 
}

//--------------------------------------------------
void Partitioner::BNDDelete(const unsigned tid, unsigned n, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, unsigned i){
     do {
     assert(bndptr[i] != -1); 
     bndind[bndptr[i]] = bndind[--(n)]; 
     bndptr[bndind[n]] = bndptr[i]; 
     bndptr[i] = -1;  
   } while(0); 
}

//--------------------------------------------------
void Partitioner::setNum(const unsigned tid, std::vector<unsigned>& where, unsigned num){
     for(unsigned i=0; i<nVtces; ++i){
  fprintf(stderr,"\nBreaking 3\n"); 
         where[i] = num;
     }
}

//--------------------------------------------------
void Partitioner::gCopy(const unsigned tid, std::vector<unsigned>& gWhere, std::vector<unsigned>& where, const std::vector<unsigned> whereDst){
     for(std::vector<unsigned>::const_iterator vit=whereDst.begin(); vit != whereDst.end(); ++vit){
         gWhere[*vit] = where[*vit];
     }
}

//--------------------------------------------------
void Partitioner::sCopy(const unsigned tid, std::vector<unsigned>& bndind, std::vector<unsigned>& bndptr, std::vector<unsigned>& partitionBndInd, std::vector<unsigned>& partitionBndPtr){
      for(unsigned i=0; i<nVtces; ++i){
         partitionBndInd.at(i) = bndind[i];
         partitionBndPtr.at(i) = bndptr[i];
     }
}
