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
//  nparts = nParts;
 // nvertices = nVertices;

  totalPECuts = new IdType[nCols];
  totalKeysInFile = new IdType[nCols];
  totalCombined = new IdType[nCols];
  readNext = new IdType[nCols];
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nRows * nCols];
  //nCuts = new IdType[nCols];
//  nEdges = new IdType[nRows * nCols];
    where = new std::vector<unsigned>[nCols];
//  partitionBndInd = new std::vector<unsigned>[nCols];
//  partitionBndPtr = new std::vector<unsigned>[nCols];
  
  outBufMap = new InMemoryContainer[nRows * nCols];
  readBufMap = new InMemoryContainer[nCols];
  refineMap = new InMemoryContainer[nCols];
  bndIndMap = new LookUpTable[nCols];
  lookUpTable = new LookUpTable[nCols];
//  gainTable = new InMemTable[nCols];
  dTable = new InMemTable[nCols];
  fetchBatchIds = new std::set<unsigned>[nCols];
  readNextInBatch = new std::vector<unsigned long long>[nCols];
  batchesCompleted = new std::vector<bool>[nCols];
  keysPerBatch = new std::vector<unsigned>[nCols];
  
  
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
    refineMap[i].clear();
    bndIndMap[i].clear();
  }
  delete[] refineMap;
  delete[] readNext;
  delete[] totalPECuts;
  delete[] totalCombined;
  delete[] bndIndMap;
//  delete[] gainTable;
  delete[] dTable;
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
void Partitioner::writeBuf(const unsigned tid, const unsigned to, const unsigned from, const unsigned hIdSize = 0){
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
     
    pthread_mutex_lock(&locks[bufferId]);
//    fprintf(stderr,"\nTID %d outbufmap buffer %d full with %d records noItems:%d\n", tid, buffer, outBufMap[buffer].size(), nItems[buffer]);
   
    writeToInfinimem(bufferId, totalKeysInFile[bufferId], outBufMap[buffer].size(), outBufMap[buffer]);
    totalKeysInFile[bufferId] += nItems[buffer];
    pthread_mutex_unlock(&locks[bufferId]);

    outBufMap[buffer].clear();
    nItems[buffer] = 0;
    writtenToDisk = true;
  }

  performWrite(tid, buffer, to, from);
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

fprintf(stderr,"\nFlushing buffer %d residues\n", tid);
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
  for (InMemoryConstIterator it = inMemMap.begin(); it != inMemMap.end(); ++it) {              fprintf(stderr,"\n \n WTI- TID: %d,  Key: %d\n", buffer, it->first);
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
/*void Partitioner::flushBResidues(const unsigned tid) {

  fprintf(stderr,"\nFlushing buffer residues tid: %d\n", tid);
     
  if(tid >= nCols) 
    return;
  
//  if(nCols == 1) {
    writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
    outBufMap[tid].clear();
    totalKeysInFile[tid] += nItems[tid];
    nItems[tid] = 0;
//  }
   
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
// fprintf(stderr,"\nInitiating In-Memory refine TID %d\n", tid);
   for(unsigned i=0; i<nRows; ++i) {
     refineMap[tid].insert(outBufMap[tid + nCols * i].begin(), outBufMap[tid + nCols * i].end());
   }
}

//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid) {
     ComputeBECut(tid, gWhere, bndIndMap[tid], refineMap[tid]);
}
//--------------------------------------------------
void Partitioner::ComputeBECut(const unsigned tid, const std::vector<unsigned>& where, LookUpTable& bndind, const InMemoryContainer& inMemMap) {

// do not need first flag now
//  pthread_mutex_lock(&locks[tid]);
//  gCopy(tid);
//  pthread_mutex_unlock(&locks[tid]);
   // fprintf(stderr,"\nInside ComputeBECUT TID %d -------\n", tid); 
   //  unsigned pid = tid % nCols;
//     IdType *bndind, *bndptr;
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
              bndind[dst].push_back(src); // TODO: should this be where[dst]?
//     fprintf(stderr,"\nSize bndind[%d]: %d\n", dst, bndind[dst].size());
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
   fprintf(stderr,"\n Partition: %d has max cuts: %d", hipart, (totalPECuts[hipart]/2));
   return hipart;
}

//--------------------------------------------------
void Partitioner::refinePart(const unsigned tid, const unsigned hipart, unsigned tCuts) {
     refinePart(tid, hipart, tCuts, gWhere, markMax, markMin);

}

//--------------------------------------------------
void Partitioner::updateDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax, unsigned src, unsigned dst){
//Go through adjacency (boundary) vertices of the masked vertices and update their DVals
  //fprintf(stderr,"\nUPDVAL SRC: %d, DST: %d\n", src, dst);
  auto it_map = refineMap[hipart].find(src);
  for (auto vit = it_map->second.begin(); vit != it_map->second.end(); ++vit) {
    //   unsigned adjvtx = *vit;
       if(gWhere.at(*vit) == whereMax){ // the vertex should belong to the other partition being refined with
         auto it_bnd = bndIndMap[hipart].find(*vit); 
         if(it_bnd != bndIndMap[hipart].end()){
       // these vertices are connected to src .. need to find their connect with dst
          auto it_dst = refineMap[whereMax].find(dst);
          unsigned conDst = 0;
          if (std::find(it_dst->second.begin(), it_dst->second.end(), *vit) == it_dst->second.end()){
        
	      conDst = 0;
         }        
         else
              conDst = 1;
		
             unsigned connect = conDst - 1;
            fprintf(stderr,"\nSRC: %d DST: %d vtx: %d CONNECT: %d ", src, dst, *vit, connect);
        	unsigned currval = dTable[whereMax].at(*vit);    
             dTable[whereMax].at(*vit) = currval + 2 * connect;
 
     //      fprintf(stderr,"\nUpdated DVAL for vertex %d is %d   ", *vit, dTable[whereMax].at(*vit));
          }
          else {
          // not a boundary vertex
       //     fprintf(stderr,"\nDVAL %d is not a bndry vertex in part %d ", *vit, hipart);
            continue;
          }
      }
  }

}

//--------------------------------------------------
void Partitioner::computeDVals(const unsigned tid, const unsigned hipart, const unsigned whereMax) {
//Go through each key in the selected partition to be refined and update the DVals
//  fprintf(stderr,"\nDVAL HIPART: %d, whereMax: %d\n", hipart, whereMax);
  for (auto it = dTable[hipart].begin(); it != dTable[hipart].end(); ++it) {
      unsigned src = it->first;
//      std::map<unsigned, unsigned>::const_iterator it_max = markMax.find(src);
//      std::map<unsigned, unsigned>::const_iterator it_min = markMin.find(it_max->second);
//              if(it_max != markMax.end() || it_min != markMin.end()){
//                 continue;
//		}
//	      else{
              auto it_bnd = bndIndMap[whereMax].find(src); 
                if(it_bnd != bndIndMap[whereMax].end()){
                   unsigned inDeg = dTable[hipart][src] - bndIndMap[whereMax][src].size();
                   dTable[hipart].at(it_bnd->first) = bndIndMap[whereMax][src].size() - inDeg;
    //  fprintf(stderr,"\nDVAL for vertex %d is %d   ", src, dTable[hipart].at(it_bnd->first));
                }
               else {
               // not a boundary vertex
//                fprintf(stderr,"\nDVAL src %d is not a bndry vertex ", src);
                dTable[hipart].erase(src);
               }
     //   }
   }

}

//--------------------------------------------------
unsigned Partitioner::computeGain(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin){
  int ct = -1;
  int  maxG = 0;
  int maxvtx = -1, minvtx = -1;
  unsigned vtx_ind = -1;
//fprintf(stderr,"\nCOMPUTING Gain hipart:%d, whereMax: %d\n", hipart, whereMax);
//Go through each key in the selected partition to be refined and update the gain 
  for (auto it = dTable[hipart].begin(); it != dTable[hipart].end(); ++it) {
      unsigned src = it->first; //maxvtx;
  //    auto it_map = refineMap[whereMax].find(src);
  //    if(it_map != refineMap[whereMax].end()){
      std::map<unsigned, unsigned>::const_iterator it_max = markMax.find(src);
      if(it_max != markMax.end()){ 
         continue;
	}
      else{
         for (auto it_hi = dTable[whereMax].begin(); it_hi != dTable[whereMax].end(); ++it_hi) {
            unsigned dst = it_hi->first;
            bool connect = 0;
      std::map<unsigned, unsigned>::const_iterator it_min = markMin.find(dst);
      if(it_min != markMin.end()){
         continue;
	}
	else{
             ct++;
                
          auto it_map = refineMap[hipart].find(src);
          if (std::find(it_map->second.begin(), it_map->second.end(), dst) == it_map->second.end()){
                  // if(gWhere[src] != gWhere[dst]){
			connect = 0;
                 }      
                 else
                       connect = 1;
		
    //     fprintf(stderr,"\nSRC: %d DST: %d CONNECT: %d ", src, dst, connect);
		unsigned dsrc = dTable[hipart][src];
                unsigned ddst = dTable[whereMax][dst];
                if(!connect)      
             	  gainTable[ct] = dsrc + ddst;
           	else
              	  gainTable[ct] = dsrc + ddst - 2;

               int currGain = gainTable[ct];
     //    fprintf(stderr,"\nGainTable: %d MaxG: %d ", gainTable[ct], maxG);
             if(currGain > maxG){
               maxG = currGain;
   //            fprintf(stderr,"\n MAX GAIN till now %d  ", maxG);
               maxvtx = src;
               minvtx = dst;
               vtx_ind = ct;
             }
 //     fprintf(stderr,"\nGAIN for vertex %d with dst %d is %d\n", src, dst, gainTable[ct]);
              }
          }
       }
    }
       //   else
//	    fprintf(stderr,"\nATTENTION - %d vertex not present in-memory\n", src); 
             
//  }
    if(maxvtx != -1 && minvtx != -1){
    fprintf(stderr,"\nMASKING %d and %d with max gain %d \n", maxvtx, minvtx, maxG);
    markMax[maxvtx] = minvtx;  
    markMin[minvtx] = maxG;  
    return maxvtx;
   }

return -1;
}

//--------------------------------------------------
void Partitioner::refinePart(const unsigned tid, const unsigned hipart, unsigned tCuts, std::vector<unsigned>& gWhere, std::map<unsigned, unsigned>& markMax, std::map<unsigned, unsigned>& markMin) {
//Check if moving the vertex with max bnd value reduce the num of cuts in the parition and how much does it impact in the current partition also check the total edgecuts
  unsigned invalidmoves = 0;
  partRefine = false; 
  // in the partition with highest edge cuts, find the bnd vertex with max occurrence
  unsigned maxvtx = -1, whereMax, vtx1;
//  unsigned minvtx = -1;
//  unsigned tCuts = countTotalPECut(tid);
  fprintf(stderr,"\n----TID: %d, Total Edge cuts : %d\n", tid, tCuts);
  unsigned pCut = getTotalPECuts(hipart);  //total partition edge cuts
  std::vector<unsigned> boundpart {nCols};  // to keep track of partitions
  
  for (InMemoryConstIterator it = bndIndMap[hipart].begin(); it != bndIndMap[hipart].end(); ++it) {
      	boundpart.push_back(it->first);
  }

    maxvtx = maxBound(tid, bndIndMap[hipart]); // start with vertex with max boundaries
  //    std::map<unsigned, unsigned>::iterator it_max = markMax.find(maxvtx);
    fprintf(stderr," in partition %d \n", hipart);
      if(maxvtx != -1 ){
         whereMax = gWhere[maxvtx];  //where is the vertex with max boundary positions in this selected partition
         fprintf(stderr,"\nREFINEPART Max vertex %d belongs to partition %d ", maxvtx, whereMax);
   if(getWrittenToDisk()){ 
   partRefine = true; 
   refineInit(hipart); cread(hipart); 
   refineInit(whereMax); cread(whereMax);
  }
       // compute dVals for boundary vertices in both the selected partitions
       computeDVals(tid, hipart, whereMax); computeDVals(tid, whereMax, hipart);

 // Loop through all the boundary vertices in the selected partition (hipart)
 while(boundpart.size() > 0){
       unsigned vtx_ind = computeGain(tid, hipart, whereMax, markMax, markMin);
   // Check the edgecuts in the partition where this Maxvertex is currently
   // Get the vtx with min edgecuts in this partition
   if(vtx_ind != -1){
      vtx1 = vtx_ind;
         bndIndMap[hipart].erase(vtx1);
        fprintf(stderr,"\nREMOVED vtx %d from part %d bndIndMap size %d ", vtx1, hipart, bndIndMap[hipart].size());
          bndIndMap[whereMax].erase(markMax.at(vtx1));
        fprintf(stderr,"\nREMOVED vtx %d from part %d bndIndMap size %d ", markMax.at(vtx1), whereMax, bndIndMap[whereMax].size());
          gainTable.clear();
     
   fprintf(stderr,"\nDtable hipart %d size %d " , hipart, dTable[hipart].size());
   fprintf(stderr,"\nDtable whereMax %d size %d " , whereMax, dTable[whereMax].size());
     updateDVals(tid, hipart, whereMax, vtx1, markMax.at(vtx1)); 
     updateDVals(tid, whereMax, hipart, markMax.at(vtx1), vtx1);
     fprintf(stderr,"\nBndInd %d has %d elements after iteration", hipart, bndIndMap[hipart].size());
     if(boundpart.size() > 0){
      boundpart.erase(std::remove(boundpart.begin(), boundpart.end(), vtx1), boundpart.end());
     // fprintf(stderr,"\nremoved %d from Boundpart size: %d", maxvtx, boundpart.size());
    }
 }
  else {
       fprintf(stderr,"\n NO positive gain found !!!!");
        fprintf(stderr,"\n OR Partition %d is refined ***** ", hipart);
     for(std::map<unsigned, unsigned>::const_iterator it = markMax.begin(); it!=markMax.end(); ++it){
         vtx1 = it->first;
         unsigned vtx2 = it->second;
         gWhere.at(vtx1) = whereMax;
         gWhere.at(vtx2) = hipart;
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx1, whereMax, hipart);
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx2, hipart, whereMax);
          changeWhere(tid, hipart, whereMax, gWhere, vtx1, vtx2);
     fprintf(stderr,"\nREFINEPART AFTER  where[%d]: %d ", vtx1, gWhere.at(vtx1));
   fprintf(stderr,"\nREFINEPART AFTER  where[%d]: %d ", vtx2, gWhere.at(vtx2));
        }
        bndIndMap[hipart].clear();
        boundpart.clear();
        markMax.clear();
        markMin.clear();
        dTable->clear();
        if(getWrittenToDisk())
	   refineMap->clear();
      }
  }
      // Delete the boundary vertices processed and update again
    //  deletebndvert(tid, hipart, whereMax, markMax);
 }
   else{
     for(std::map<unsigned, unsigned>::const_iterator it = markMax.begin(); it!=markMax.end(); ++it){
         vtx1 = it->first;
         unsigned vtx2 = it->second;
         gWhere.at(vtx1) = whereMax;
         gWhere.at(vtx2) = hipart;
     fprintf(stderr,"\nREFINEPART AFTER  where[%d]: %d ", vtx1, gWhere[vtx1]);
   fprintf(stderr,"\nREFINEPART AFTER  where[%d]: %d ", vtx2, gWhere[vtx2]);
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx1, whereMax, hipart);
          fprintf(stderr,"\nMoving %d to %d from %d ", vtx2, hipart, whereMax);
          changeWhere(tid, hipart, whereMax, gWhere, vtx1, vtx2);
        }
        fprintf(stderr,"\nPartition %d is refined ***** ", hipart);
        bndIndMap[hipart].clear();
        boundpart.clear();
        markMax.clear();
        markMin.clear();
        dTable->clear();
        if(getWrittenToDisk())
	   refineMap->clear();
      }
}
//--------------------------------------------------
void Partitioner::deletebndvert(const unsigned tid, const unsigned hipart, const unsigned whereMax, std::map<unsigned, unsigned>& markMax ) {
  for(auto it = markMax.begin(); it!= markMax.end(); ++it){
     bndIndMap[hipart].erase(it->first);
     bndIndMap[whereMax].erase(it->second);
     fprintf(stderr,"\nDeleting %d, bndIndMap Hipart size %d", it->first, bndIndMap[hipart].size());
     fprintf(stderr,"\nDeleting %d, bndIndMap whereMax size %d", it->second, bndIndMap[whereMax].size());
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
unsigned Partitioner::findMaxGain(const unsigned tid){
  unsigned maxG = -1;
  unsigned ct = -1;
//  unsigned vtx = -1
  for(auto it = gainTable.begin(); it != gainTable.end(); ++it){
//	maxG = it->second;
     fprintf(stderr,"\nit->second: %d MaxG: %d ", it->second, maxG);
	if(it->second > maxG){
	   ct = it->first;
           maxG = it->second;
	}
  }

//unsigned vtx_ind = -1;
/*  if(ct != -1){
    unsigned dTSize = gainTable.size()/dTable[hipart].size();
    vtx_ind = ct > dTSize ? dTSize % ct : vtx_ind;

  fprintf(stderr, "\nVertex %d has max Gain of %d ", dTable[hipart][vtx_ind], maxG);

  return dTable[hipart][vtx_ind];
 }
*/
 
  fprintf(stderr, "\nVertex index %d has max Gain of %d ", ct, maxG);
return ct;
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
 // fprintf(stderr,"\nREFINE tid: %d, totalCombined: %d, readNext[tid]: %d, keys to read: %d \n", tid, totalCombined[tid], readNext[tid], partbound);
 //for(unsigned ckey = readNextInBatch[partition]; ckey < totalCombined[partition]; ckey += batchSize){
    if (partbound > 0 && readNext[tid] < totalCombined[tid])
      cio->file_get_batch(tid, readNext[tid], partbound, parts); 

  //  fprintf(stderr,"\nREFINE- tid: %d  BATCHSIZE: %d, Map size: %d, readNext: %d \n", tid, batchSize, refineMap[tid].size(), readNext[tid]);
    for (unsigned i = 0; i < partbound; i++) {
      refineMap[tid][parts[i].rank()];
//      fprintf(stderr,"\nREFINE- TID: %d, Key: %d\t Values: ", tid, parts[i].rank()); 

      for (unsigned k = 0; k < parts[i].nbrs_size(); k++){
        refineMap[tid][parts[i].rank()].push_back(parts[i].nbrs(k));
         if(getPartRefine()){
           dTable[tid][parts[i].rank()] = parts[i].nbrs_size();
//           fprintf(stderr,"\nREFINE - DTABLE key: %d value: %d\n", parts[i].rank(), parts[i].nbrs_size());
         }
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
    // for in-memory reads
      execLoop = 0;
    } 
   else{
       execLoop = refine(tid);
    }
  //  fprintf(stderr, "\nExecloop is %d, TID %d", execLoop, tid);

    if(execLoop == false) {
    //  for(InMemoryConstIterator it = refineMap[tid].begin(); it != refineMap[tid].end(); ++it){
    //     fprintf(stderr,"\nCREAD- tid: %d, Computing edgecuts with Map Size %d", tid, refineMap[tid].size());
       if(getWrittenToDisk() && !getPartRefine()){
         ComputeBECut(tid);
         refineMap[tid].clear(); //erase(refineMap[tid].begin(), it) erase when writing to disk
      }
       break;
    }

   // Read the kitems from infinimem and compute the edgecuts for that part

  if(!getPartRefine())
    ComputeBECut(tid);
    refineMap[tid].erase(refineMap[tid].begin(), refineMap[tid].end());
  }
//fprintf(stderr,"\nCREAD - RefineMap %d size: %d\n", tid, refineMap[tid].size());

}

//--------------------------------------------------
//void Partitioner::changeWhere(const unsigned tid, const unsigned hipart, const unsigned chVtx ) {

//     changeWhere(tid, hipart, gWhere, chVtx);
//}
  
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
     if(gWhere[i] != -1 ){ //&& gWhere[i] == tid){
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
   for (auto it = refineMap[tid].begin(); it != refineMap[tid].end(); ++it) {
          dTable[tid][it->first] = it->second.size();
   }
   fprintf(stderr,"\nDtable tid %d size %d " , tid, dTable[tid].size());
}
