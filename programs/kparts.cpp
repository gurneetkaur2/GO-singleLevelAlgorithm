//#define MAX_WORD_SIZE 32
//#define MAX_REAL_WORD_SIZE 32

#include "edgeList.pb.h"
//#include "graph.h"
#define USE_NUMERICAL_HASH

#include "../engine/gp.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "pthread.h"
#include <ctime>
#include <cstdlib>

//RecordType* edgeLists = NULL;
static std::string outputPrefix = "";
//-------------------------------------------------
// Let the game begin 
//-------------------------------------------------

class KParts : public GraphParts
{
  static thread_local std::ofstream ofile;
  static thread_local double stime;
  public:
    
  //  void* readEdgeLists(const unsigned tid, graph_t origG, RecordType* edgeLists, const std::string& input)
    void* MParts(const unsigned tid, const std::string& input)
    {
      if(input[0] == '%' || input[0] == '#')
      return NULL;

      std::stringstream inputStream(input);
      IdType to;
      std::vector<unsigned> from;
      std::vector<unsigned> one;
      one.push_back(1);
     // while(inputStream >> to >> from) {
      while(inputStream >> to) {
//      while(inputStream.peek() != '\n' && inputStream.peek() != EOF){
//           inputStream >> to;
//           int i = 0;
//           IdType token;
//           while(inputStream.peek() != '\n' && inputStream.peek() != EOF){
//		inputStream >> token;
//		from.push_back(token);
//           }
//           fprintf(stderr,"\nto: %d\t", to);
//            edgeLists[to].add_nbrs(from);
	  writeBuf(tid, to, one);
//           for(int i=0; i<from.size(); i++){
//   		fprintf(stderr,"from: %d\t", from[i]); 
//                writeBuf(tid, from[i], one);
//           }
      }
      return NULL;
    
    }

  void* beforeRefine(const unsigned tid) {
    ofile.open(outputPrefix + std::to_string(tid));
    stime = 0.0;
    return NULL;
  }
    void* refine(const unsigned tid, const unsigned& rank, const std::vector<unsigned>& nbrs) {
//	uint64_t total = std::accumulate(nbrs.begin(), nbrs.end(), 0);
        stime -= getTimer();
  //  	ofile << rank << " " << total << std::endl;
    	ofile << rank << std::endl;
    	stime += getTimer();     
          return NULL;
    }

   void* afterRefine(const unsigned tid) {
    ofile.close();
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

//-------------------------------------------------
int main(int argc, char** argv)
{
  KParts kp;
  if (argc != 9)
  {
std::cout << "Usage: " << argv[0] << " <fileName> <nvertices> <nedges> <nparts> <nrefiners> <batchsize> <kitems> <outputprefix>" << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  IdType nvertices = atoi(argv[2]);
  IdType nedges = atoi(argv[3]);
  int ninmemparts = atoi(argv[4]);
  int nreducers = atoi(argv[5]);
  int batchSize = atoi(argv[6]);
  int kitems = atoi(argv[7]);
  int nparts = atoi(argv[4]);
  outputPrefix = argv[8];

  assert(batchSize > 0);
  
  kp.init(fileName, ninmemparts, nreducers, batchSize, kitems, nparts);

  double runTime = -getTimer();
  kp.run(); 
  runTime += getTimer();
  
  std::cout << "Main::Run time : " << runTime << " (msec)" << std::endl;

  return 0;
}

