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
void Partitioner::initg(graph_t graph, unsigned nCoarseners, unsigned nRefiners, unsigned bSize, unsigned kItems, unsigned nParts)
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

 //Initialize graph
  graph_t g = graph;


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
void Partitioner::writeBuf(const unsigned tid, const IdType to, const IdType from) {

  unsigned bufferId = hashKey(to) % nCols; 
  unsigned buffer = tid * nCols + bufferId;  

  if (outBufMap[buffer].size() >= batchSize) {
    fprintf(stderr,"outbufmap buffer full with %d records\n", outBufMap[buffer].size());
     initsubgraph(tid, buffer);
   // pthread_mutex_lock(&locks[bufferId]);
   // totalKeysInFile[bufferId] += nItems[buffer];
   // pthread_mutex_unlock(&locks[bufferId]);

    outBufMap[buffer].clear();
    nItems[buffer] = 0;
   // nEdges[buffer] = 0;
    writtenToDisk = true;
  }

  performWrite(tid, buffer, to, from);
  fprintf(stderr,"total edges in buffer: %d\n", nEdges[buffer]);
}

//--------------------------------------------------
void Partitioner::initsubgraph(const unsigned tid, const unsigned buffer){
  IdType *xadj, *adjncy;
// *vwgt, *adjwgt;
  IdType edge, i ,k;
     graph_t localg, cgraph;
     localg.v = outBufMap[buffer].size();
     localg.e = nEdges[buffer];

  fprintf(stderr,"localg.v: %d, localg.e: %d\n", localg.v, localg.e);
     xadj = (IdType *) malloc((localg.v + 1) * sizeof(IdType));
     adjncy = (IdType *) malloc(((localg.e * 2) + 1) * sizeof(IdType));

     xadj[0]=0, i=0, k=0;
  fprintf(stderr,"Allocated memory to xadj\n");
    
 
     for (InMemoryContainerIterator it = outBufMap[buffer].begin(); it != outBufMap[buffer].end(); ++it) {
  fprintf(stderr,"Inside for loop\n");
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
     
     localg.num_edges = xadj;
     localg.adj = adjncy;
     localg.edgeWeight = (IdType *) malloc(localg.e * sizeof(IdType));

  for(int i = 0; i < localg.e; i++)
                localg.edgeWeight[i] = 1;


   localg.vertexWeight = (IdType *) malloc(localg.v * sizeof(IdType));;

   for(int i = 0; i < localg.v; i++)
                localg.vertexWeight[i] = 1;



  fprintf(stderr,"Calling coarsen\n");
     coarsen(tid, localg, cgraph, 2);

// free (localg.num_edges); free(localg.adj); free(localg.edgeWeight); free(localg.vertexWeight);
// free (cgraph.num_edges); free(cgraph.adj); free(cgraph.edgeWeight); free(cgraph.vertexWeight);
}

//--------------------------------------------------
void Partitioner::flushBResidues(const unsigned tid) {

  fprintf(stderr,"\nFlushing buffer residues\n");
  if(tid >= nCols) 
    return;

  if(nRows == 1) {
      initsubgraph(tid, tid);
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
          initsubgraph(tid, i);
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
               initsubgraph(tid, b2);
 //           bWriteToInfinimem(tid, totalKeysInFile[tid], outBufMap[b2].size() - b2Merged, b2Iter, b2End);
            outBufMap[b2].clear();
            totalKeysInFile[tid] += (nItems[b2] - b2Merged);
            nItems[b2] = 0;
          }
        } 
        if(b2Iter == b2End) {
          if(outBufMap[b1].size() != 0) {
              initsubgraph(tid, b1);
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
          initsubgraph(tid, b1);
//        writeToInfinimem(tid, totalKeysInFile[tid], outBufMap[b1].size(), outBufMap[b1]);
        outBufMap[b1].clear();
        totalKeysInFile[tid] += nItems[b1];
        nItems[b1] = 0;
      } else {  
           initsubgraph(tid, b2);
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
      initsubgraph(tid, tid);   // TODO::buffer maynot be correct
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
//--------------------------------------------------
void Partitioner::performWrite(const unsigned tid, const unsigned buffer, const IdType to, const IdType from) {
  std::vector<unsigned> nbrs{from};
  InMemoryContainerIterator it = outBufMap[buffer].find(to); 

  if(it != outBufMap[buffer].end()){
     combine(to, it->second, nbrs);
//      outBufMap[buffer][to].push_back(from);
      nEdges[buffer]++;
      }
  else {
    outBufMap[buffer].emplace(to, nbrs); 
    nItems[buffer]++;
    nEdges[buffer]++;
  }
}

//------------------------------------------------- GK
// Write the Map output to in-memory Buffer
// Gets Hashing function from Application Programmer
//void Partitioner::coarsen(const unsigned tid, const graph_t cgraph, const unsigned CoarsenTo, const unsigned int* numEdgesSupRowsToRows, const unsigned int* mapSupRowstoRows){

void Partitioner::coarsen(const unsigned tid, const graph_t g, const graph_t& cgraph, const unsigned CoarsenTo) {
     IdType *numEdgesSupRowsToRows = NULL;
     IdType *mapSupRowstoRows = (IdType *)malloc(100 * sizeof(IdType));

     std::vector<IdType *> mappingVectors;
     std::vector<IdType> vectorSizes;
     

     fprintf(stderr,"\nInside Coarsen\n");
      graph_t localG = g;

      IdType *perm;
      while(localG.v > CoarsenTo) {
	perm = (IdType *) malloc(localG.v * sizeof(IdType));
     fprintf(stderr,"Calling MATCH_RM LocalG.v: %d, CoarsenTo: %d, perm: %d\n", localG.v, CoarsenTo, perm);
        Match_RM(localG, perm, cgraph);

        vectorSizes.push_back( localG.v);
        mappingVectors.push_back(perm);
	fprintf(stderr,"Updating the vectors Size vectorSizes: %d, mappingVectors: %d \n", vectorSizes.size(), mappingVectors.size());
        localG = cgraph;
      }

      //Find the final mapping from the coarsest graph to the fine graph 
    //  findFinalMapping(localG.v, vectorSizes, mappingVectors, numEdgesSupRowsToRows, mapSupRowstoRows);

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

void Partitioner::findFinalMapping(IdType numSupRows, std::vector<IdType> & vectorSizes, std::vector<IdType *> & mappingVectors, IdType* & numEdgesSupRowsToRows, IdType* & mapSupRowstoRows) {

 std::vector< std::vector<IdType> > finalAdjLists;
 std::vector< std::vector<IdType> > tempVectAdjLists;

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
				std::vector<IdType> curAdjList;
				for(int k = 0; k < finalAdjLists[j].size(); k++) {
					int finerRowIndex = finalAdjLists[j][k];
					for(int l = 0; l < tempVectAdjLists[finerRowIndex].size(); l++) 
						curAdjList.push_back( tempVectAdjLists[finerRowIndex][l] );
				}
				finalAdjLists[j] = curAdjList;
			}
		}

	}


	numEdgesSupRowsToRows = (IdType *)malloc( (numSupRows + 1) * sizeof(IdType));

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

//------------------------------------------------- GK
void Partitioner::Match_RM(graph_t g, IdType *perm, graph_t coarsenedG) {

     fprintf(stderr,"\nInside MATCH_RM\n");
	//Matched edges   
        std::set<Edge> matchEdges; 

	//Initially vertices are not matched 
        bool isMatched[g.v];

	srand ( unsigned ( std::time(0) ) );
	std::vector<IdType> myvector;

     fprintf(stderr,"Inside MATCH_RM- before initializing g.v: %d\n", g.v);
	//Initialize these values
	for (int i = 0; i < g.v; ++i) {
                isMatched[i] = false; 
		myvector.push_back(i);
        }
	fprintf(stderr,"MATCH_RM - after initialization  myvector: %d \n", myvector.size());
 
	//Randomly shuffles the vertices
	std::random_shuffle ( myvector.begin(), myvector.end() );

     fprintf(stderr,"Inside MATCH_RM - before visiting a vertex\n");
	//Visit a vertex at random and randomly match it with one of its neighbors
	for (int i = 0; i < g.v; ++i) {
                int currVtx = myvector[i];
     fprintf(stderr,"Inside MATCH_RM- currVtx: %d, myvector[i]: %d, isMatched[currVtx]: %d\n", currVtx, myvector[i], isMatched[currVtx]);

                if( !isMatched[currVtx] ) {
			//Make a vector of unmatched neighbors	
			std::vector<IdType> neighbors;
			
     fprintf(stderr,"\nInside MATCH_RM- making vector of unmatched neighbors g.num_edges[currVtx]: %d\n", g.num_edges[currVtx]);
			for(int j = g.num_edges[currVtx]; j < g.num_edges[currVtx+1]; j++) {
//				fprintf(stderr,"currVtx: %d, g.adj[%d]: %d, isMatched[g.adj[j]]: %d \n", currVtx, j, g.adj[j], isMatched[g.adj[j]]);
				if( (g.adj[j] != currVtx) && !isMatched[ g.adj[j] ] ) {
					neighbors.push_back( g.adj[j] );
  //   				fprintf(stderr,"Inside MATCH_RM- neighbors[%d]: %d\n", i, neighbors[i]);
				}
			}

			//Randomly select one of its unmatched neighbors	
     fprintf(stderr,"Inside MATCH_RM- randomly select one of the unmatched neighbors\n");
			int numVtxs = neighbors.size();
			if( numVtxs == 1) {	
				isMatched[ currVtx ] = true;
                                isMatched [ neighbors[0] ] = true;
                                matchEdges.insert( Edge(currVtx, neighbors[0]) );
			}
			else if( numVtxs > 1) {
				int randIndex = rand() % numVtxs;	
				isMatched[ currVtx ] = true;
                                isMatched [  neighbors[ randIndex] ] = true;         
                                matchEdges.insert( Edge(currVtx, neighbors[ randIndex]) ); 
			}	
                }
        }

        //Count the vertices in the coarsened graph 
        //Store a mapping from coarsend vertex to the vertices in the finer graph        
        int numConargenedVertices = 0;
        std::vector< std::vector<IdType> > superRowsToRows;


	/*	
        //Match edges
        cout<<"Matched edges\n";
        for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
                    unsigned int u = it->either();
                    unsigned int v = it->other(u);
                    cout<<"Edge("<<u<<","<<v<<") "<<"\n";
	}*/	

     fprintf(stderr,"Inside MATCH_RM- counting vertices in the coarsened graph\n");
	for (std::set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
		    IdType u = it->either();
		    IdType n = it->other(u);
	
                    perm[u] = numConargenedVertices;
                    perm[n] = numConargenedVertices;

                    std::vector<IdType> supRowContains(2,0);
		    supRowContains[0] = u;
		    supRowContains[1] = n;
		    superRowsToRows.push_back(supRowContains);
                    numConargenedVertices++;
        }
	
	/*	
	cout<<"Unmatched vertices\n";	
        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] )
			cout<<i<<"\n"; 
	} */   


        for(int i = 0; i < g.v; ++i) {
               if( !isMatched[i] ) {
			perm[i] = numConargenedVertices;
			std::vector<IdType> supRowContains(1,i);
	                superRowsToRows.push_back(supRowContains);
			numConargenedVertices++;
               }
        }

        coarsenedG.v = numConargenedVertices;

	//Calculate the vertex weights of the coarsened graph
	coarsenedG.vertexWeight = (IdType *)malloc( numConargenedVertices * sizeof(IdType));
     fprintf(stderr,"Inside MATCH_RM- calcultaing vertex weight\n");
	for(int i = 0; i < coarsenedG.v; i++) { 
		IdType currVtxWeight = 0;
		for(int j = 0; j < superRowsToRows[i].size(); j++) {
			 currVtxWeight = currVtxWeight + g.vertexWeight[ superRowsToRows[i][j] ];
		}
		coarsenedG.vertexWeight[ i ] = currVtxWeight;
	}

        coarsenedG.num_edges = (IdType *)malloc( (coarsenedG.v+1) * sizeof(IdType) );

	//vectors to compute adjacency of coarsened graph and edge weights
	std::vector< IdType > adjVector;
	std::vector< IdType > weightVector;
        
	coarsenedG.num_edges[0] = 0;
	
     fprintf(stderr,"Inside MATCH_RM- computing adjacency\n");
	for(int i = 0; i < coarsenedG.v; i++) {

		//Maps neighbor to weight
		std::map<IdType, IdType> adjNeighborWeightMap;

		for(int j = 0; j < superRowsToRows[i].size(); j++) {
                      IdType vtxInCurSuperRow = superRowsToRows[i][j];
		      for(int j = g.num_edges[vtxInCurSuperRow]; j < g.num_edges[vtxInCurSuperRow+1]; j ++) {
				IdType adjVtx = g.adj[j];
				IdType connectedSupRow = perm[ adjVtx ];
			
				std::map<IdType, IdType>::iterator it = adjNeighborWeightMap.find( connectedSupRow );
				
				if( it != adjNeighborWeightMap.end() ) {	
					int totalWeight = g.edgeWeight[j] +  it->second;
					it->second = totalWeight;
				}
				else {
					adjNeighborWeightMap.insert( std::pair<IdType, IdType>(connectedSupRow, g.edgeWeight[j]) );
				}
                      }
                }
		
		//Sort the neighbors according to weight in increasing order
		std::vector< std::pair<IdType, IdType> > adjWeightVector( adjNeighborWeightMap.begin(), adjNeighborWeightMap.end());
		std::sort(adjWeightVector.begin(), adjWeightVector.end(), sortBySecond<IdType, IdType>()); 

                //Put the neighbors in the graph_t struct
                coarsenedG.num_edges[ i + 1 ] =  coarsenedG.num_edges[i] + adjWeightVector.size();
     fprintf(stderr,"Inside MATCH_RM- putting neighbors in the graph adjWeightVector.size(): %d\n", adjWeightVector.size());

                for(int j = 0; j < adjWeightVector.size(); j++) {
			adjVector.push_back( adjWeightVector[j].first );
			weightVector.push_back( adjWeightVector[j].second ); 
		}
	}

	coarsenedG.e = adjVector.size();
	coarsenedG.orgM = coarsenedG.e;

	assert( coarsenedG.e == coarsenedG.num_edges[coarsenedG.v] );

        coarsenedG.adj = (IdType *)malloc( coarsenedG.e  * sizeof(IdType) );
        coarsenedG.edgeWeight = (IdType *)malloc( coarsenedG.e * sizeof(IdType) );

	for(int j = 0; j < coarsenedG.e; j++) {
		coarsenedG.adj [j] = adjVector[j];
		coarsenedG.edgeWeight [j] = weightVector[j];
	}
}

