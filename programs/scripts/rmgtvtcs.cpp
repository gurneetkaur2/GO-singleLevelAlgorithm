#include<iostream>
#include<fstream>
#include <sstream>
#include<map>
#include<vector>
#include<cassert>

using namespace std;
// Remove the max vertex from each adjacency list in a file
int main(int argc, char** argv)
{
  //   ifstream ifile("/home/gkaur007/work/datasets/links-anon.txt");
  //  ofstream ofile("/home/gkaur007/work/datasets/link-anon.txt_done");
  if (argc != 3)
  {
    cout << "Usage: " << argv[0] << " <fileName> <maxVtx> " << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  unsigned maxvtx = atoi(argv[2]);
  
  std::string outFile = fileName.c_str();
  outFile += "_cnvt";

  ifstream ifile(fileName.c_str());
  assert(ifile.is_open());
  ofstream ofile(outFile.c_str());
  assert(ofile.is_open());

  std::map<unsigned long long, std::vector<unsigned long long> > adjL;
  unsigned long long num_vertices = 0;
  unsigned long long num_edges = 0;
  std::string line;
  ifile.seekg(std::ios::beg);
  unsigned long long i = 0;
 
  std::getline(ifile, line);
  while(std::getline(ifile, line)){
    unsigned long long from;
    std::stringstream iss(line);
    ++i;
    while(iss >> from) {
      std::vector<unsigned> vals(from);
  //      if(to < maxvtx){
  //        cout << "\n to: " << i << " from: " << from << endl;
          if(from <= maxvtx){
         // cout << "\n From: " << from << endl;
           adjL[i].push_back(from);
           num_edges++;
          }
    //    }
     }
  }
  
  
  //for(std::map<unsigned long long, std::vector<unsigned long long> >::const_iterator it = adjL.begin(); it != adjL.end(); ++it){
    // num_edges += it->second.size();

 // }
  ofile << maxvtx << " " << num_edges/2 << "\n";
  //unsigned long long i = 0;
  unsigned long long numEdges = 0;
  for(std::map<unsigned long long, std::vector<unsigned long long> >::const_iterator it = adjL.begin(); it != adjL.end(); ++it){
    //  i++;
     numEdges += it->second.size();
     for (std::vector<unsigned long long>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
            ofile << " " << *vit << " ";
    // std::cout <<"Edge " << *vit << "\tsize "<< it->second.size() << "\n";
    //  sum++;
     }
    ofile << "\n";
    

  }
  //  sum = keys.size();
  fprintf(stderr,"\nFile Name: %s ",fileName.c_str());
  fprintf(stderr,"\n Total Vertices %llu,  num_edges %llu, Edges: %llu ", i, num_edges, numEdges/2);
  ofile.close();
  ifile.close();
  return 0;
}

