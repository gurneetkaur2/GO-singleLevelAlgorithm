//#define MAX_WORD_SIZE 32
//#define MAX_REAL_WORD_SIZE 32

#define USE_NUMERICAL_HASH
#include "edgeList.pb.h"
#include "recordtype.pb.h"
#include "adjacencyList.pb.h"

#include "../engine/gp.hpp"
#include "graph.h"
//#include "../engine/infinimem/hm.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "pthread.h"
#include <ctime>
#include <cstdlib>
#include <google/malloc_extension.h>

using namespace std;

//Meta data for partition processing
/*typedef struct __partitionInfo {
  IdType lbEdgeCount;
  IdType ubEdgeCount;
  IdType edgeCount;

  IdType lbIndex;
  IdType ubIndex;
  IdType indexCount;
} PartitionInfo;
PartitionInfo *pi = NULL;
*/
//-------------------------------------------------
// Let the game begin 
//-------------------------------------------------

static std::string outputPrefix = "";
class KParts : public GraphParts
{
  static thread_local std::ofstream ofile;
  static thread_local double stime;
  public:
    
  //  void* readEdgeLists(const unsigned tid, graph_t origG, RecordType* edgeLists, const std::string& input)
    void* createMParts(const unsigned tid, const std::string& input, const unsigned lineId, const unsigned hiDegree) 
    {
	std::stringstream inputStream(input);
        std::vector<unsigned> from;
        unsigned hid = 0; 
        unsigned token;
//       std::vector<IdType> from;
 
        while(inputStream >> token){
	    from.push_back(token);
         }
          for(unsigned i = 0; i < from.size(); ++i){
    //        fprintf(stderr,"\nVID: %d FROM: %zu size: %zu", lineId, from[i], from.size());
            if(from.size() < hiDegree){
              writeBuf(tid, lineId, from[i], hid);
             }
            else{
              hid = from.size();
              writeBuf(tid, lineId, from[i], hid);
             }
        }
      return NULL;
  }

  void* beforeRefine(const unsigned tid) {
    //ofile.open(outputPrefix + std::to_string(tid));
    //stime = 0.0;
    return NULL;
  }
    void* refine(const unsigned tid, const unsigned& rank, const std::vector<unsigned>& nbrs) {
//	uint64_t total = std::accumulate(nbrs.begin(), nbrs.end(), 0);
//        stime -= getTimer();
  //  	ofile << rank << " " << total << std::endl;
 //   	ofile << rank << std::endl;
 //   	stime += getTimer();     
          return NULL;
    }

   void* afterRefine(const unsigned tid, const unsigned nVertices) {
  //  ofile.open(outputPrefix + std::to_string(tid));
    stime = 0.0;
    std::string fileName = outputPrefix + std::to_string(tid);
    printParts(tid, fileName.c_str());
//    ofile.close();
    this->subtractRefineTimes(tid, stime);
    return NULL;
  }
};

thread_local std::ofstream KParts::ofile;
thread_local double KParts::stime;

void* combine(const unsigned& key, std::vector<unsigned>& to, const std::vector<unsigned>& from) {
//  assert(to.size() == 1);
//  assert(from.size() == 1);
  to.insert(std::end(to), std::begin(from), std::end(from));
 // to[0] += from[0];
  return NULL;
}
/*
//-------------------------------------------------
void readEdgeLists(EdgeList* edgeLists, string fileName, IdType nVertices, IdType nEdges) {
  std::ifstream infile(fileName);
  assert(infile.is_open());
  string line;
  
  IdType k, i;
  xadj = (IdType *) malloc((nVertices + 1) * sizeof(IdType));
  adjncy = (IdType *) malloc(((nEdges * 2) + 1) * sizeof(IdType));

  xadj[0]=0, k=0, i=0;
  while (getline(infile,line)){
        std::istringstream iss(line);
        IdType to, from, edge;
  
        iss >> to;
        edge = to;
        adjncy[k] = edge - 1;
        fprintf(stderr,"\nEdgeLists TO: %d", to);
        fprintf(stderr,"\ndst: %d, adjncy: %d, edge-1: %d\n", to, adjncy[k], edge-1);
        k++;
        while (iss >> from) {
              edge = from;
              edgeLists[to].add_nbrs(from);
              totalEdges++;
              fprintf(stderr,"\tFROM: %d", from);
              adjncy[k] = edge - 1;
              fprintf(stderr,"\nsrc: %d, adjncy: %d, edge-1: %d\n", from, adjncy[k], edge-1);
              k++;
             // from.push_back(value);
        }
        xadj[i+1] = k;
        fprintf(stderr,"xadj: %d, k: %d, i: %d\n", xadj[i+1], k, i);
        i++;
  }
}
*/
//-------------------------------------------------
int main(int argc, char** argv)
{
  KParts kp;
  if (argc != 9)
  {
std::cout << "Usage: " << argv[0] << " <fileName> <nvertices> <nedges> <hDegree> <nparts> <batchsize> <kitems> <outputprefix>" << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  IdType nvertices = atoi(argv[2]);
  IdType nedges = atoi(argv[3]);
  unsigned nthreads = atoi(argv[5]);
  unsigned batchSize = atoi(argv[6]);
  unsigned kitems = atoi(argv[7]);
  unsigned nparts = atoi(argv[5]);
  unsigned hDegree = atoi(argv[4]);
  outputPrefix = argv[8];
  unsigned edgesPerMPart = (nedges/nparts) + 1;

  fprintf(stderr, "total partitions: %zu\n", nparts);
  fprintf(stderr, "total edges: %zu\n", nedges);
  fprintf(stderr, "Edges per partition: %zu\n", edgesPerMPart);

  
  assert(batchSize > 0);

  kp.init(fileName, nvertices, nedges, hDegree, nthreads, nparts, batchSize, kitems);
  fprintf(stderr,"\nCreating partitions ..");
  double runTime = -getTimer();
  kp.run(); 
  runTime += getTimer();

  std::cout << endl << "Main::Run time : " << runTime << " (msec)" << std::endl;

  return 0;
}

