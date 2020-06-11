#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

//#include "./metis/include/metis.h"

#define BASE_PATH "../inputs/"
#define PARTS_PATH "parts_"

#define SOURCE 0

using namespace std;
typedef int64_t idx_t;

// Convert the graph and its partitions to a transformed graph

typedef struct edgeList
{
  set<unsigned long long> edges;
} EdgeList;

int main(int argc, char* argv[])
{
  if(argc !=4){
    cout << "\nUsage: <input adj_graph> <nparts> <num_vertices>\n";
  assert(argc == 4);
  }
  string graphName = argv[1];
  string transformedGraphName = graphName + ".newGraph";
 // string adjGraphName = "adj_" + graphName;
  string boundariesName = graphName + ".boundaries";
//  string commName = graphName + ".comm";
  string sourceName = graphName + ".src";

  idx_t nparts = atoi(argv[2]);
  unsigned long long  num_vertices = atoi(argv[3]);

  string partsDir = "";
  partsDir += PARTS_PATH;
  partsDir += argv[2];
  partsDir += "/";

  mkdir((BASE_PATH + partsDir).c_str(), 0777); 

 /* unsigned long long num_edges = 0, num_vertices = 0;
  {
    unsigned long long to, from, edge;
    string graphPath = BASE_PATH + graphName;
    std::ifstream infile(graphPath.c_str());
    while (infile >> to >> from >> edge) {
      num_vertices = max(num_vertices, max(to, from));
      ++num_edges;
    }
  }
  ++num_vertices;
  num_edges *= 2;

  idx_t nvtxs = num_vertices;
  idx_t ncon = 1; // TODO

  idx_t* xadj = new idx_t[num_vertices + 1];
  idx_t* adjncy = new idx_t[num_edges]; 
  idx_t* vsize = new idx_t[num_vertices];
  idx_t* vwgt = new idx_t[num_vertices];
  //real_t ubvec[] = {1.5}; // TODO
  real_t* ubvec = NULL;
*/
 // EdgeList* edgeLists = new EdgeList[num_vertices];
  ++num_vertices;
  EdgeList* ogEdgeLists = new EdgeList[num_vertices];

  string graphPath = BASE_PATH + graphName;

  std::ifstream infile(graphPath.c_str());
 // std::vector<unsigned> from;
  unsigned long long from;
  unsigned long long to = 0;
 // unsigned edge;
  std::string line;
//  ogEdgeLists[to].edges.insert(to); // for vertex 0
  while (std::getline(infile, line)){
          to++;
     //     cout <<"\n TO: " << to ;
         std::istringstream iss(line);
         while (iss >> from){
//  cout << "\t from: " << from << endl;
        ogEdgeLists[to].edges.insert(from);
          } 
  }
/*  unsigned long long numEdges = 0;
  for(unsigned long long i=0; i<num_vertices; ++i)
  {
    numEdges += edgeLists[i].edges.size();  
    vsize[i] = 1;
    vwgt[i] = edgeLists[i].edges.size();  
  }

  string adjPath = BASE_PATH + adjGraphName;
  ofstream graphFile;
  graphFile.open(adjPath.c_str());

  graphFile << num_vertices << " " << numEdges/2 << " 110 1\n"; 

  unsigned long long adj_i = 0;
  for(unsigned long long i=0; i<num_vertices; ++i)
  {
    graphFile << "1 1";
    xadj[i] = adj_i;
    set<unsigned long long>::iterator it = edgeLists[i].edges.begin();
    while(it != edgeLists[i].edges.end())
    {
      graphFile << " " << (*it) + 1;
      adjncy[adj_i++] = *it;
      ++it;    
    }
    graphFile << "\n";
  }
  xadj[num_vertices] = adj_i;

  graphFile.close();
*/
  /*
     idx_t options[METIS_NOPTIONS];

     options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
     options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
     options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
     options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM;
     options[METIS_OPTION_NUMBERING] = 0;
     options[METIS_OPTION_NITER] = 10; // TODO
     options[METIS_OPTION_SEED] = 29;
     options[METIS_OPTION_CONTIG] = 0;
     options[METIS_OPTION_MINCONN] = 0;
     options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
     */

/*  idx_t* options = NULL;

  idx_t objval;
  idx_t* part = new idx_t[nvtxs];
  part[0] = 256;

  if(nparts <= 8)
    METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, vsize, NULL, &nparts, NULL, ubvec, options, &objval, part);
  else
    METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, vsize, NULL, &nparts, NULL, ubvec, options, &objval, part);

  string commPath = BASE_PATH + partsDir + commName;
  ofstream commFile;
  commFile.open(commPath.c_str());
  commFile << "Communication cost: " << objval << endl; // << "Partitions:" << endl;
  //for(unsigned long long i=0; i<nvtxs; ++i)
  //cout << part[i] << endl;
*/

  // new vertex IDs of vertices in the same partition
  unsigned long long* vertexToNewVertex = new unsigned long long[num_vertices];
  // stores the old vertex IDs 
  unsigned long long* newVertexToVertex = new unsigned long long[num_vertices];
  unsigned long long new_id = 1;

  unsigned long long* part_boundaries = new unsigned long long[nparts];
  idx_t* part = new idx_t[num_vertices];
  //unsigned v;

  string partsPath = BASE_PATH + partsDir + graphName + to_string(nparts);
  cout <<"\nPartPath: " << partsPath  << "\tfileName " << graphName;
  std::ifstream partFile(partsPath.c_str());
  // partInFile.open(partsPath.c_str());
   //assert(partInFile.is_open());
   unsigned long long p = 1;
     unsigned long long v;

     // READ part file and match *******
  //   part[p] = 0;
    // ++p;
   while (partFile >> v){
    // if(p >= num_vertices) break;
    //  cout <<"\n part " << v << " p " << p << endl;
       part[p] = v;
      ++p;
   }
   partFile.close();
     /*
    while (p <= num_vertices) {
       partInFile >> v;
      cout <<"\n part " << v << " p " << p << endl;
         part[p] = v;
      p++;
    }*/
  for(idx_t src=0; src < nparts; ++src)
  {
  //  cout << "\n src: " << src << endl;
    for(unsigned long long i=1; i<num_vertices; ++i){
      if(part[i] == src)   // match the parts read with the src (number of parts)
      {
        vertexToNewVertex[i] = new_id;  //mapping new id at index i
  //    cout <<"\n i: " << i << " part[i]: " <<part[i] << " new_id: " << new_id << endl;
  //    cout <<  "\nvtoNewV: " << vertexToNewVertex[i] << endl;
        newVertexToVertex[new_id] = i;  //stores actual vertex
  //    cout <<"\n ** newVtoV: " << newVertexToVertex[new_id] << endl;
        ++new_id;
      }
     // ++i;
    }
    part_boundaries[src] = new_id;
  }

  for(unsigned long long i=1; i<num_vertices; ++i)
  {
    set<unsigned long long>::iterator it = ogEdgeLists[i].edges.begin();
    set<unsigned long long> new_set;
    while(it != ogEdgeLists[i].edges.end())
    {
      new_set.insert(vertexToNewVertex[*it]);
//      cout <<"\ni: " << i << "\t *it: " << *it << endl;
      cout <<"\n swapping " << *it << " with " <<  vertexToNewVertex[*it] << endl;
      ++it;
    }
    ogEdgeLists[i].edges.clear();
    ogEdgeLists[i].edges.swap(new_set);
  }

  string transPath = BASE_PATH + partsDir + transformedGraphName;
  cout <<"\ntransPath: " << transPath  << "\n";
  ofstream transformedGraphFile;
  transformedGraphFile.open(transPath.c_str());
  for(unsigned long long i=1; i<num_vertices; ++i)
  {
 // cout <<"\nGoing to for loop " << "\n";
    set<unsigned long long>::iterator it = ogEdgeLists[newVertexToVertex[i]].edges.begin();
    while(it != ogEdgeLists[newVertexToVertex[i]].edges.end())
    {
      transformedGraphFile << i << " " << *it << endl;
   //   cout <<"\n" << i << " " << *it << endl;
      ++it; 
    }
  }

  transformedGraphFile.close();

  string srcPath = BASE_PATH + sourceName;
  ofstream srcFile;
  srcFile.open(srcPath.c_str());
  srcFile << vertexToNewVertex[SOURCE] << endl;
  srcFile.close();

  string boundsPath = BASE_PATH + boundariesName;
  ofstream boundariesFile;
  boundariesFile.open(boundsPath.c_str());
  for(int i=0; i<nparts; ++i)
    boundariesFile << part_boundaries[i] << endl;
  boundariesFile.close();

  delete[] ogEdgeLists;
  delete[] part;
  delete[] vertexToNewVertex;
  delete[] newVertexToVertex;
  delete[] part_boundaries;
}
