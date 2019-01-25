#include "partitioner.h"
#include "coarsen.h"
#include "kway.h"


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
void Partitioner::initg(unsigned nCoarseners, unsigned nRefiners, unsigned bSize, unsigned kItems, unsigned nParts)
{
  //nBuffers = pow(buffers, 2); // number of buffers is square of number of threads
  //nBuffers = nMappers * nReducers;

  nRows = nCoarseners;
  nCols = nRefiners;
  writtenToDisk = false;
  batchSize = bSize;
  kBItems = kItems;
  nparts = nParts;
 // nvertices = nVertices;

  //cTotalKeys = new IdType[nBuffers];
  totalKeysInFile = new IdType[nCols];
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nRows * nCols];
  nEdges = new IdType[nRows * nCols];
 
  outBufMap = new InMemoryContainer[nRows * nCols];
  
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
void Partitioner::releaseMapStructures()
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

     fprintf(stderr,"\nInside WriteBuf\n");
  unsigned bufferId = hashKey(to) % nCols; 
  unsigned buffer = tid * nCols + bufferId;  

  if (outBufMap[buffer].size() >= batchSize) {
     graph_t *cgraph, *localg; 
     unsigned CoarsenTo;
     fprintf(stderr,"outbufmap buffer full with %d records\n", outBufMap[buffer].size());
     
     localg = (graph_t *)malloc((outBufMap[buffer].size()) * sizeof(IdType));
     localg = initsubgraph(tid, buffer);

    CoarsenTo = max((localg->v)/(20*gk_log2(nparts)), 30*(nparts));
     
  fprintf(stderr,"Calling coarsen CoarsenTo; %d vertices: %lu\n", CoarsenTo, localg->v);
     cgraph = coarsen(tid, localg, 2);

   // pthread_mutex_lock(&locks[bufferId]);
   // totalKeysInFile[bufferId] += nItems[buffer];
   // pthread_mutex_unlock(&locks[bufferId]);

    outBufMap[buffer].clear();
    nItems[buffer] = 0;
    nEdges[buffer] = 0;
    writtenToDisk = true;
  }

  performWrite(tid, buffer, to, from);
 // fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
graph_t* Partitioner::initsubgraph(const unsigned tid, const unsigned buffer){
  IdType *xadj, *adjncy;
// *vwgt, *adjwgt;
  IdType edge, i ,k;
     fprintf(stderr,"\nInside initsubgraph\n");
     graph_t *localg;
     localg = CreateGraph();

     localg->v = outBufMap[buffer].size();
     localg->e = nEdges[buffer];

  fprintf(stderr,"localg.v: %d, localg.e: %d\n", localg->v, localg->e);
     xadj = (IdType *) malloc((localg->v + 1) * sizeof(IdType));
     adjncy = (IdType *) malloc(((localg->e * 2) + 1) * sizeof(IdType));

     xadj[0]=0, i=0, k=0;
    
 
     for (InMemoryContainerIterator it = outBufMap[buffer].begin(); it != outBufMap[buffer].end(); ++it) {
           fprintf(stderr,"\n First : %d\n", it->first);
	   //Add the "to" in the adjncy as well
           adjncy[k++] = it->first - 1;  
           std::vector<unsigned> vtemp = it->second;
           for (unsigned v=0; v<vtemp.size(); v++){
		edge = vtemp[v];
                adjncy[k] = edge - 1;
		k++;       
           }
  fprintf(stderr,"i: %d, k: %d\n", i, k);
           xadj[i+1] = k;
           i++;
    }
     
     localg->xadj = xadj;
     localg->adjncy = adjncy;
     localg->adjwgt = (IdType *) malloc(localg->e * sizeof(IdType));

  for(int i = 0; i < localg->e; i++)
                localg->adjwgt[i] = 1;


   localg->vwgt = (IdType *) malloc(localg->v * sizeof(IdType));

   for(int i = 0; i < localg->v; i++)
                localg->vwgt[i] = 1;

    localg->vsize = (IdType *)malloc((localg->v) * sizeof(IdType)); 
      for(int i = 0; i < localg->v; i++)
                localg->vsize[i] = 1;

//   localg->vsize     = NULL;
   localg->cmap      = NULL;
   localg->coarser   = NULL;
   localg->finer     = NULL;
   localg->ncon      = 1;

   /* by default these are set to true, but the can be explicitly changed afterwards */
  localg->free_xadj   = 1;
  localg->free_vwgt   = 1;
  localg->free_vsize  = 1;
  localg->free_adjncy = 1;
  localg->free_adjwgt = 1;

   localg->label     = NULL;
   localg->where     = NULL;
   localg->pwgts     = NULL;
   localg->bndptr    = NULL;
   localg->bndind    = NULL;

   SetupGraph_label(localg);

   return localg;
}

//--------------------------------------------------
void Partitioner::performWrite(const unsigned tid, const unsigned buffer, const IdType to, const std::vector<unsigned>& from) {
//  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it = outBufMap[buffer].find(to); 

  if(it != outBufMap[buffer].end()){
     combine(to, it->second, from);
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
void Partitioner::flushBResidues(const unsigned tid) {

  fprintf(stderr,"\nFlushing buffer residues\n");
     graph_t *cgraph, *localg; 
     unsigned CoarsenTo;
     
     localg = (graph_t *)malloc((outBufMap[tid].size()) * sizeof(IdType));

    CoarsenTo = max((localg->v)/(20*gk_log2(nparts)), 30*(nparts));
     
  if(tid >= nCols) 
    return;

  if(nRows == 1) {
      localg = initsubgraph(tid, tid);
  fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
     cgraph = coarsen(tid, localg, 2);

  //  writeToInfinimem(tid, totalKeysInFile[tid], static_cast<unsigned>(outBufMap[tid].size()), outBufMap[tid]);
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
          localg = initsubgraph(tid, i);
  fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
     cgraph = coarsen(tid, localg, 2);

//        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[i].size(), outBufMap[i]);
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
               localg = initsubgraph(tid, b2);
  		fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
     		cgraph = coarsen(tid, localg, 2);

 //           bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
            outBufMap[b2].clear();
            totalKeysInFile[tid] += (nItems[b2] - b2Merged);
            nItems[b2] = 0;
          }
        } 
        if(b2Iter == b2End) {
          if(outBufMap[b1].size() != 0) {
              localg = initsubgraph(tid, b1);
  		fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
       		cgraph = coarsen(tid, localg, 2);

//            writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
            outBufMap[b1].clear();
            totalKeysInFile[tid] += nItems[b1];
            nItems[b1] = 0;
          }
        }
      } else
          assert(false);
    } else if(i == nRows) {
      if(outBufMap[b1].size() > 0) {
          localg = initsubgraph(tid, b1);
  	  fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
          cgraph = coarsen(tid, localg, 2);

//        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
        outBufMap[b1].clear();
        totalKeysInFile[tid] += nItems[b1];
        nItems[b1] = 0;
      } else {  
           localg = initsubgraph(tid, b2);
  	   fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
     	   cgraph = coarsen(tid, localg, 2);

//        bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
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
     graph_t *cgraph, *localg; 
     unsigned CoarsenTo;
     
     localg = (graph_t *)malloc((toMap.size()) * sizeof(IdType));

    CoarsenTo = max((localg->v)/(20*gk_log2(nparts)), 30*(nparts));

      localg = initsubgraph(tid, tid);   // TODO::buffer maynot be correct
     fprintf(stderr,"Calling coarsen CoarsenTo; %d\n", CoarsenTo);
      cgraph = coarsen(tid, localg, 2);

 //     writeToInfinimem(tid, totalKeysInFile[tid], toMap.size(), toMap);
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
