#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

//#include "metis.h"

//#define BASE_PATH "../inputs/newgraphs_24/"
#define BASE_PATH "../newgraphs_24/"
//#define PARTS_PATH "parts_"

#define SOURCE 0

using namespace std;
typedef int64_t idx_t;
//typedef int32_t idx_t;
//typedef float real_t;

typedef struct edgeList
{
  set<unsigned long long> edges;
} EdgeList;

int main(int argc, char* argv[])
{
  assert(argc == 2);
  string graphName = argv[1];
  //  string transformedGraphName = graphName + ".newGraph";
  string adjGraphName = "adj_" + graphName;
  //string boundariesName = graphName + ".boundaries";
  // string commName = graphName + ".comm";
  string sourceName = graphName + ".src";

  unsigned long long num_edges = 0, num_vertices = 0;
  {
    unsigned long long to, from, edge;
    string graphPath = BASE_PATH + graphName;
    std::ifstream infile(graphPath.c_str());
    //infile >> to >> from;
    while (infile >> to >> from ) {
      num_vertices = max(num_vertices, max(to, from));
      ++num_edges;
    }
  }
 // cout <<"\nV: " << num_vertices-1 << " " << num_edges << "\n";  //<< " 110 1\n"; 
  ++num_vertices; 
  num_edges *= 2;

  EdgeList* edgeLists = new EdgeList[num_vertices];
 // EdgeList* ogEdgeLists = new EdgeList[num_vertices];

  string graphPath = BASE_PATH + graphName;

  std::ifstream infile(graphPath.c_str());
  unsigned long long to, from;
 //  unsigned edge;

   // infile >> to >> from;
  while (infile >> to >> from ) {
    edgeLists[to].edges.insert(from);
    edgeLists[from].edges.insert(to);
   // ogEdgeLists[to].edges.insert(from);
  }

  unsigned long long numEdges = 0;
  for(unsigned long long i=0; i<num_vertices; ++i)
  {
    numEdges += edgeLists[i].edges.size();
    // vsize[i] = 1;
    // vwgt[i] = edgeLists[i].edges.size();  
  }

  string adjPath = BASE_PATH + adjGraphName;
  ofstream graphFile;
  cout <<"\nWriting to AdjPath: " <<adjPath <<"\n";
  graphFile.open(adjPath.c_str());

  graphFile << num_vertices-1 << " " << numEdges/2 << "\n";  //<< " 110 1\n"; 
  //cout <<"\nV: " << num_vertices-1 << " " << numEdges/2 << "\n";  //<< " 110 1\n"; 

 // unsigned long long adj_i = 0;
  for(unsigned long long i=1; i<num_vertices; ++i)
  {
   // if(i == 0)
     // graphFile << "1 1";
    set<unsigned long long>::iterator it = edgeLists[i].edges.begin();
    while(it != edgeLists[i].edges.end())
    {
//      cout << "\n i: " << i <<" it: " <<(*it);
      graphFile << (*it)  << " ";
      ++it;
    }
    graphFile << "\n";
  }
  graphFile.close();

  delete[] edgeLists;
 // delete[] ogEdgeLists;
}
