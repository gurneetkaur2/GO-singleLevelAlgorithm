
#define USE_NUMERICAL_HASH
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

  void* createMParts(const unsigned tid, const std::string& input, const std::string& type, const unsigned fileId, const unsigned hiDegree) 
  {
    std::stringstream inputStream(input);
    std::vector<unsigned> from;
    unsigned hid = 0; 
    unsigned token, to;
    std::map<unsigned, unsigned> keys;
    unsigned long long num_vertices = 0; 

    if (type == "edgelist"){
      inputStream >> to;
      inputStream >> token;
      std::map<unsigned, unsigned>::iterator it_to = keys.find(to);

      // This map is just to calculate the hidegree
      if(it_to != keys.end()){
        keys[to] = it_to->second + 1;

        if((it_to->second + 1) > hiDegree)
          hid = it_to->second;
      }
      else {
        keys[to] = 1;
      }

      num_vertices = max(num_vertices, max(to, token));
      writeBuf(tid, to, token, hid);
    }

    // Adjacency list format
    else {
      inputStream >> to;
      while(inputStream >> token){
        from.push_back(token);
      }
      for(unsigned i = 0; i < from.size(); ++i){
        if(from.size() < hiDegree){
          writeBuf(tid, to, from[i], hid);
        }
        else{
          hid = from.size();
          writeBuf(tid, to, from[i], hid);
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
    return NULL;
  }


  //-------------------------------------------------
  int main(int argc, char** argv)
  {
    KParts kp;
    if (argc != 9)
    {
      std::cout << "Usage: " << argv[0] << " <folderpath> <fileType> <nvertices> <hDegree> <nthreads> <nparts> <batchsize> <outputprefix>" << std::endl;
      return 0;
    }

    std::string fileType = "";
    std::string folderpath = argv[1];
    fileType = argv[2]; 
    IdType nvertices = atoi(argv[3]);
    unsigned nthreads = atoi(argv[5]);
    unsigned memSize = atoi(argv[7]);
    unsigned kitems = 20;
    unsigned nparts = atoi(argv[6]);
    unsigned hDegree = atoi(argv[4]);
    outputPrefix = argv[8];

    if(fileType != "edgelist" && fileType != "adjlist" ){
      fprintf(stderr, "\nFile Type %s not accepted, please select edgelist or adjlist format \n", fileType.c_str());
      return 0;
    }
    
    if(nparts < 2){
      fprintf(stderr, "\n Number of graph partitions cannot be less than 2 \n")
      return 0;
      }

    fprintf(stderr, "total partitions: %zu\n", nparts);
    
    if(hDegree < 1000)
       hDegree = 1000;

    if(memsize <= 0)
       memsize = nVertices;

    //assert(memSize > 0);

    kp.init(folderpath, fileType, nvertices, hDegree, nthreads, nparts, memSize, kitems);
    fprintf(stderr,"\nCreating partitions ..");
    double runTime = -getTimer();
    kp.run(); 
    runTime += getTimer();

    std::cout << endl << "Main::Run time : " << runTime << " (msec)" << std::endl;

    return 0;
  }

