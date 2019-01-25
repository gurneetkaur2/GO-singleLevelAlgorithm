#include "partitioner.h"

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
void Partitioner::initg(unsigned nVertices, unsigned nCoarseners, unsigned nRefiners, unsigned bSize, unsigned kItems, unsigned nParts)
{
  //nBuffers = pow(buffers, 2); // number of buffers is square of number of threads
  //nBuffers = nMappers * nReducers;

  nRows = nCoarseners;
  nCols = nRefiners;
  writtenToDisk = false;
  batchSize = bSize;
  kBItems = kItems;
  nparts = nParts;
  nvertices = nVertices;
  //cTotalKeys = new IdType[nBuffers];
  totalKeysInFile = new IdType[nCols];
  //nReadKeys = new IdType[nBuffers];
  nItems = new IdType[nRows * nCols];

  
  for (unsigned i=0; i<nRows * nCols; ++i) 
    nItems[i] = 0;
  
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

  //delete[] cTotalKeys;
  delete[] nItems;
}
//-------------------------------------------------
void Partitioner::shutdown()
{
  delete io;
  delete[] totalKeysInFile;
  //delete[] nReadKeys;
}

//------------------------------------------------- GK
// Write the Map output to in-memory Buffer
// Gets Hashing function from Application Programmer
//void Partitioner::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const unsigned int* mapSupRowstoRows){

void Partitioner::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo) {
     unsigned int *numEdgesSupRowsToRows = NULL;
     unsigned int *mapSupRowstoRows = (unsigned int *)malloc(100 * sizeof(unsigned int));

     std::vector<unsigned int *> mappingVectors;
     std::vector<unsigned int> vectorSizes;

     fprintf(stderr,"Inside Coarsen\n");
      graph_t localG;
      localG.v = nvertices;

      unsigned int *perm;
      while(localG.v > CoarsenTo) {
	perm = (unsigned int *) malloc(localG.v * sizeof(unsigned int));
     fprintf(stderr,"Calling MATCH_RM\n");
        Match_RM(localG, perm, cgraph);
        vectorSizes.push_back( localG.v);
        mappingVectors.push_back(perm);
        localG = cgraph;
      }

      //Find the final mapping from the coarsest graph to the fine graph 
      findFinalMapping(localG.v, vectorSizes, mappingVectors, numEdgesSupRowsToRows, mapSupRowstoRows);

	fprintf(stderr, " number of coarsest vertices: %d \n", cgraph.v);

	for(int i = 0; i < cgraph.v; i++) {
		fprintf(stderr, "Coarsest vertex: %d constains vertices: \n", i);
		for(int j = numEdgesSupRowsToRows[i]; j < numEdgesSupRowsToRows[i+1]; j++)
			fprintf(stderr,"%d \t", mapSupRowstoRows[j]);
		fprintf(stderr,"\n");
	}


      //free memory
        free( numEdgesSupRowsToRows );
      free( mapSupRowstoRows );

       for(int i = 0; i < mappingVectors.size(); i++){
           if( mappingVectors[i] != NULL )
               free( mappingVectors[i] );
      }
}

//------------------------------------------------- GK
/* Given intermediate mappinf of vertices find the final mapping from the coarsest vertices to the finest vertices
 * Input: numSupRows - number of coarsest vertices
 * 	  vectorSizes - vector containing number of vertices in intermediate coarsened graphs
 *	  mappingVectors - vector containing the mapping of fine vertices to coarse vertices 
 * Output: numEdgesSupRowsToRows - points to mapSupRowstoRows to indicate where each coarsest row starts
 * 	   mapSupRowstoRows - contains the vertices of the finest graph
 */
//------------------------------------------------- GK

void Partitioner::findFinalMapping(unsigned int numSupRows, std::vector<unsigned int> & vectorSizes, std::vector<unsigned int *> & mappingVectors, unsigned int* & numEdgesSupRowsToRows, unsigned int* & mapSupRowstoRows) {

 std::vector< std::vector<unsigned int> > finalAdjLists;
 std::vector< std::vector<unsigned int> > tempVectAdjLists;

	for(int i = vectorSizes.size() -1; i >= 0; i--) {
	
		if( i == vectorSizes.size() -1 ) {
			finalAdjLists.resize( numSupRows );

			for(int j = 0; j < vectorSizes[i]; j++) 
				finalAdjLists[ mappingVectors[i][j] ].push_back( j );
		}
		else {
			tempVectAdjLists.resize(vectorSizes[i+1]);
			for(int j = 0; j < vectorSizes[i+1]; j++) 
				tempVectAdjLists[ j ].resize(0);
			

			for(int j = 0; j < vectorSizes[i]; j++) 
				tempVectAdjLists[ mappingVectors[i][j] ].push_back( j );
		
			for(int j = 0; j < numSupRows; j++) {
				std::vector<unsigned int> curAdjList;
				for(int k = 0; k < finalAdjLists[j].size(); k++) {
					int finerRowIndex = finalAdjLists[j][k];
					for(int l = 0; l < tempVectAdjLists[finerRowIndex].size(); l++) 
						curAdjList.push_back( tempVectAdjLists[finerRowIndex][l] );
				}
				finalAdjLists[j] = curAdjList;
			}
		}

	}


	numEdgesSupRowsToRows = (unsigned int *)malloc( (numSupRows + 1) * sizeof(unsigned int));

	numEdgesSupRowsToRows[0] = 0;
	int mapIndex = 0;
	for(int i = 0; i < numSupRows; i++) {
		numEdgesSupRowsToRows[i+1] = numEdgesSupRowsToRows[i] + finalAdjLists[i].size();
		for(int j = 0; j < finalAdjLists[i].size(); j++) {
			mapSupRowstoRows[mapIndex] = finalAdjLists[i][j];
			mapIndex++;
		}
	}	

}


/*************************************************************************/
/*! This function finds a matching by randomly selecting one of the 
    unmatched adjacent vertices. 
 */
/**************************************************************************/
void Partitioner::Match_RM(graph_t g, unsigned int *perm, graph_t coarsenedG) {

  idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx, last_unmatched;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
  idx_t *match, *cmap, *perm;
  size_t nunmatched=0;


//  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->MatchTmr));

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;

  maxvwgt  = ctrl->maxvwgt;

  match = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
  perm  = iwspacemalloc(ctrl, nvtxs);

  irandArrayPermute(nvtxs, perm, nvtxs/8, 1);

  for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) {
          last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }
        else {
          /* Find a random matching, subject to maxvwgt constraints */
          if (ncon == 1) {
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
                maxidx = k;
                break;
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && 3*vwgt[i] < maxvwgt[0]) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
          else {
            /* multi-constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  ivecaxpylez(ncon, 1, vwgt+i*ncon, vwgt+k*ncon, maxvwgt)) {
                maxidx = k;
                break;
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && ivecaxpylez(ncon, 2, vwgt+i*ncon, vwgt+i*ncon, maxvwgt)) {
              nunmatched++;
              maxidx = UNMATCHED;
            }
          }
        }
      }

      if (maxidx != UNMATCHED) {
        cmap[i]  = cmap[maxidx] = cnvtxs++;
        match[i] = maxidx;
        match[maxidx] = i;
      }
    }
  }

  //printf("nunmatched: %zu\n", nunmatched);

  /* see if a 2-hop matching is required/allowed */
  if (!ctrl->no2hop && nunmatched > UNMATCHEDFOR2HOP*nvtxs) 
    cnvtxs = Match_2Hop(ctrl, graph, perm, match, cnvtxs, nunmatched);


  /* match the final unmatched vertices with themselves and reorder the vertices 
     of the coarse graph for memory-friendly contraction */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED) {
      match[i] = i;
      cmap[i]  = cnvtxs++;
    }
    else {
      if (i <= match[i]) 
        cmap[i] = cmap[match[i]] = cnvtxs++;
    }
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match);

  WCOREPOP;

  return cnvtxs;
}
