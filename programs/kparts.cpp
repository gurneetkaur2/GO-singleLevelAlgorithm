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
//-------------------------------------------------
// Let the game begin 
//-------------------------------------------------

class KParts : public GraphParts
{

  public:
    
  //  void* readEdgeLists(const unsigned tid, graph_t origG, RecordType* edgeLists, const std::string& input)
    void* coarsenG(const unsigned tid, const std::string& input)
    {
      if(input[0] == '%' || input[0] == '#')
      return NULL;

      std::stringstream inputStream(input);
      IdType to;
      std::vector<unsigned> from;

     // while(inputStream >> to >> from) {
      while(inputStream) {
           inputStream >> to;
           int i = 0;
           IdType token;
           while(inputStream.peek() != '\n' && inputStream.peek() != EOF){
		inputStream >> token;
		from.push_back(token);
           }
           fprintf(stderr,"\nto: %d\t", to);
           for(int i=0; i<from.size(); i++){
   		fprintf(stderr,"from: %d\t", from[i]); 
           }
//            edgeLists[to].add_nbrs(from);
	  writeBuf(tid, to, from);
      }
      return NULL;
    
    }

/*    void* refine(const unsigned tid, const unsigned& key, const std::vector<unsigned>& values) {
          unsigned int number;
          return NULL;
    }*/
};

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
std::cout << "Usage: " << argv[0] << " <fileName> <nvertices> <nedges> <ncoarseners> <nrefiners> <batchsize> <kitems> <nparts>" << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  IdType nvertices = atoi(argv[2]);
  IdType nedges = atoi(argv[3]);
  int ncoarseners = atoi(argv[4]);
  int nrefiners = atoi(argv[5]);
  int batchSize = atoi(argv[6]);
  int kitems = atoi(argv[7]);
  int nparts = atoi(argv[8]);

  assert(batchSize > 0);
  
  kp.init(fileName, ncoarseners, nrefiners, batchSize, kitems, nparts);

  double runTime = -getTimer();
  kp.run(); 
  runTime += getTimer();
  
  std::cout << "Main::Run time : " << runTime << " (msec)" << std::endl;

  return 0;
}

