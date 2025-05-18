#include<iostream>
#include<fstream>
#include <sstream>
#include<map>
#include<vector>
#include<cassert>

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    cout << "Usage: " << argv[0] << " <fileName> <graph-type> " << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  std::string type = "";
  type = argv[2];
  if(type != "edgelist"){
    cout<<"\nOnly 'edgelist' to adjacency list format conversion supported\n";
    return 0;
  }

  
  std::string outFile = fileName.c_str();
  outFile += "_done";

  ifstream ifile(fileName.c_str());
  assert(ifile.is_open());
  ofstream ofile(outFile.c_str());
  assert(ofile.is_open());

  std::map<unsigned long long, std::vector<unsigned long long> > adjL;
  unsigned long long num_vertices = 0;
  unsigned long long num_edges = 0;
  std::string line;
  ifile.seekg(std::ios::beg);

  while(std::getline(ifile, line)){
    unsigned long long from, to;
    std::stringstream iss(line);
    if (type == "edgelist" ){
    while(iss >> to >> from) {
      std::vector<unsigned> vals(from);
       num_vertices = max(num_vertices, max(to, from));
        adjL[to].push_back(from);
        num_edges++;
     }
    }
  }
  
  ofile << adjL.size() << " " << num_edges/2 << "\n";
  unsigned long long i = 0;
  unsigned long long numEdges = 0;
  for(std::map<unsigned long long, std::vector<unsigned long long> >::const_iterator it = adjL.begin(); it != adjL.end(); ++it){
      i++;
     numEdges += it->second.size();
     for (std::vector<unsigned long long>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
         if(*vit <= adjL.size())
            ofile << " " << *vit << " ";
     }
    ofile << "\n";
    
  fprintf(stderr,"\nFile Name: %s , type: %s ",fileName.c_str(), type.c_str());
  fprintf(stderr,"\n Total Vertices %llu, i: %llu,  num_edges %llu, Edges: %llu ", num_vertices, i, num_edges, numEdges/2);
  ofile.close();
  ifile.close();
  return 0;
}

