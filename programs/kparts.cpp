//#define MAX_WORD_SIZE 32
//#define MAX_REAL_WORD_SIZE 32

#define USE_NUMERICAL_HASH
//#include "edgeList.pb.h"
#include "recordtype.pb.h"

#include "../engine/gp.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "pthread.h"
#include <ctime>
#include <cstdlib>
#include <google/malloc_extension.h>

using namespace std;

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
    void* createMParts(const unsigned tid, const std::string& input, const std::string& type, const unsigned lineId, const unsigned hiDegree) 
    {
	std::stringstream inputStream(input);
        std::vector<unsigned> from;
        unsigned hid = 0; 
        unsigned token, to;
//       std::vector<IdType> from;

        if (type == "edge"){
        	inputStream >> to;
        }

        while(inputStream >> token){
	    from.push_back(token);
         }
        if (type == "edge"){
          for(unsigned i = 0; i < from.size(); ++i){
    //        fprintf(stderr,"\nVID: %d FROM: %zu size: %zu", lineId, from[i], from.size());
            if(from.size() < hiDegree){
              writeBuf(tid, to, from[i], hid);
             }
            else{
              hid = from.size();
              writeBuf(tid, to, from[i], hid);
             }
           }
        }

 //   if (type == "adj") {
      else {
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


//-------------------------------------------------
int main(int argc, char** argv)
{
  KParts kp;
  if (argc != 10)
  {
    std::cout << "Usage: " << argv[0] << " <fileName> <fileType> <nvertices> <nedges> <hDegree> <nparts> <batchsize> <kitems> <outputprefix>" << std::endl;
    return 0;
  }

  std::string fileName = "";
  std::string fileType = "";
  fileName = argv[1];
  fileType = argv[2]; 
  IdType nvertices = atoi(argv[3]);
  IdType nedges = atoi(argv[4]);
  unsigned nthreads = atoi(argv[6]);
  unsigned batchSize = atoi(argv[7]);
  unsigned kitems = atoi(argv[8]);
  unsigned nparts = atoi(argv[6]);
  unsigned hDegree = atoi(argv[5]);
  outputPrefix = argv[9];
  //  unsigned edgesPerMPart = (nedges/nparts) + 1;

  if(fileType != "edge" && fileType != "adj" ){
     fprintf(stderr, "\nFile Type %s not accepted, please select edge or adj \n", fileType.c_str());
     return 0;
  }

  fprintf(stderr, "total partitions: %zu\n", nparts);
  fprintf(stderr, "total edges: %zu\n", nedges);
  //fprintf(stderr, "Edges per partition: %zu\n", edgesPerMPart);


  assert(batchSize > 0);

  kp.init(fileName, fileType, nvertices, nedges, hDegree, nthreads, nparts, batchSize, kitems);
  fprintf(stderr,"\nCreating partitions ..");
  double runTime = -getTimer();
  kp.run(); 
  runTime += getTimer();

  std::cout << endl << "Main::Run time : " << runTime << " (msec)" << std::endl;

  return 0;
}

