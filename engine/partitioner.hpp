#include "partitioner.h"
//#include "coarsen.h"
//#include "kway.h"


/*#ifdef USE_ONE_PHASE_IO
#include "infinimem/onePhaseFileIO.hpp"
#else
#include "infinimem/twoPhaseFileIO.hpp"
#endif*/
#include "infinimem/twoPhaseFileIO.hpp"

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
void Partitioner::initg(unsigned nMemParts, unsigned bSize, unsigned kItems, unsigned nReduceParts)
{
  //nBuffers = pow(buffers, 2); // number of buffers is square of number of threads
  //nBuffers = nMappers * nReducers;

  nRows = nMemParts;
  nCols = nReduceParts;
  writtenToDisk = false;
  batchSize = bSize;
  kBItems = kItems;
//  nparts = nParts;
 // nvertices = nVertices;

  //cTotalKeys = new IdType[nBuffers];
  totalKeysInFile = new IdType[nCols];
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nRows * nCols];
  nEdges = new IdType[nRows * nCols];
 
  outBufMap = new InMemoryContainer[nRows * nCols];
  readBufMap = new InMemoryContainer[nCols];
  lookUpTable = new LookUpTable[nCols];
  fetchBatchIds = new std::set<unsigned>[nCols];
  readNextInBatch = new std::vector<unsigned long long>[nCols];
  batchesCompleted = new std::vector<bool>[nCols];
  keysPerBatch = new std::vector<unsigned>[nCols];
  
  for (unsigned i=0; i<nRows * nCols; ++i){ 
    nItems[i] = 0;
    nEdges[i] = 0;
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
//#endif
}

//-------------------------------------------------
void Partitioner::releaseInMemStructures()
{
  for (unsigned i = 0; i < nCols; i++)
    pthread_mutex_destroy(&locks[i]);

   for (unsigned i = 0; i < nRows * nCols; i++)
    outBufMap[i].clear();

  //delete[] cTotalKeys;
  delete[] nItems;
  delete[] nEdges;
  delete[] outBufMap;
  
}
//-------------------------------------------------
void Partitioner::shutdown()
{
  delete io;
  delete[] totalKeysInFile;
  //delete[] nReadKeys;
}

//-------------------------------------------------
void Partitioner::writeBuf(const unsigned tid, const IdType to, const std::vector<unsigned>& from) {

//     fprintf(stderr,"\nInside WriteBuf\n");
  unsigned bufferId = hashKey(to) % nCols; 
  unsigned buffer = tid * nCols + bufferId;  

  if (outBufMap[buffer].size() >= batchSize) {
 //    fprintf(stderr,"outbufmap buffer full with %d records\n", outBufMap[buffer].size());
     
    pthread_mutex_lock(&locks[bufferId]);
      writeToInfinimem(bufferId, totalKeysInFile[bufferId], outBufMap[buffer].size(), outBufMap[buffer]);
    totalKeysInFile[bufferId] += nItems[buffer];
    pthread_mutex_unlock(&locks[bufferId]);

    outBufMap[buffer].clear();
    nItems[buffer] = 0;
    nEdges[buffer] = 0;
    writtenToDisk = true;
  }

  performWrite(tid, buffer, to, from);
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
void Partitioner::performWrite(const unsigned tid, const unsigned buffer, const IdType to, const std::vector<unsigned>& from) {
  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it_to = outBufMap[buffer].find(to); 

  if(it_to != outBufMap[buffer].end()){
     combine(to, it_to->second, from);
//      outBufMap[buffer][to].push_back(from);
      nEdges[buffer]++;
      }
  else {
    outBufMap[buffer].emplace(to, from); 
    nItems[buffer]++;
    nEdges[buffer]++;
  }
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
  RecordType* records = new RecordType[noItems]; 
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
}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

  //fprintf(stderr,"\nFlushing buffer residues\n");
     
  if(tid >= nCols) 
    return;

  if(nRows == 1) {
    writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
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
        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
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
    } else if(i == nRows) {
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
  }
}

//--------------------------------------------------
unsigned long long Partitioner::merge(InMemoryContainer& toMap, unsigned whichMap, unsigned tid, InMemoryContainerIterator& begin, InMemoryConstIterator end) {
  unsigned long long ct = 0;
  while(begin != end) {
    if(toMap.size() >= batchSize) {
      writeToInfinimem(tid, totalKeysInFile[tid], toMap.size(), toMap);
      toMap.clear();
      totalKeysInFile[tid] += nItems[whichMap];
      nItems[whichMap] = 0; 
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
}

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
}

//--------------------------------------------------
bool Partitioner::read(const unsigned tid) {
  return read(tid, readBufMap[tid], keysPerBatch[tid], lookUpTable[tid], fetchBatchIds[tid], readNextInBatch[tid], batchesCompleted[tid]);
}

//--------------------------------------------------
InMemoryReductionState Partitioner::initiateInMemoryReduce(unsigned tid) {
  InMemoryReductionState state(nRows); 
  for(unsigned i=0; i<nRows; ++i) {
    state.begins[i] = outBufMap[tid + nCols * i].begin();
    state.ends[i] = outBufMap[tid + nCols * i].end();
  }

  return state;
}

//--------------------------------------------------
bool Partitioner::getNextMinKey(InMemoryReductionState* state, InMemoryContainer* record) {
  std::vector<unsigned> minIds;
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

  return true;
}
