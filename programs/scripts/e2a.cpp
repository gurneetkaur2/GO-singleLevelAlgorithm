#include<iostream>
#include<fstream>
#include <sstream>
#include<map>
#include<vector>
#include<cassert>

using namespace std;

int main(int argc, char** argv)
{
  //   ifstream ifile("/home/gkaur007/work/datasets/links-anon.txt");
  //  ofstream ofile("/home/gkaur007/work/datasets/link-anon.txt_done");
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

//  std::map<int, int> keys;
  std::map<unsigned long long, std::vector<unsigned long long> > adjL;
  unsigned long long num_vertices = 0;
  unsigned long long num_edges = 0;
  //unsigned long long sum = 0;
  std::string line;
  ifile.seekg(std::ios::beg);

  while(std::getline(ifile, line)){
    unsigned long long from, to;
    std::stringstream iss(line);
    if (type == "edgelist" ){
    while(iss >> to >> from) {
      std::vector<unsigned> vals(from);
     // std::map<unsigned long long, std::vector<unsigned long long> >::iterator it_to = adjL.find(to);
       num_vertices = max(num_vertices, max(to, from));
    //  if(it_to != adjL.end()){
        //to.insert(std::end(to), std::begin(vals), std::end(vals));
     //std::cout <<"Adding " << from <<"\n";
        adjL[to].push_back(from);
     // }
     // else {
     //   adjL.emplace(to, from);
    // std::cout <<"EMPLACE " << from <<"\n";
        num_edges++;
        //  fprintf(stderr, "\nWriting %d to file ", to);
  //      ofile << to << endl;
     // }

     }
    }
  }
  
  
  //for(std::map<unsigned long long, std::vector<unsigned long long> >::const_iterator it = adjL.begin(); it != adjL.end(); ++it){
    // num_edges += it->second.size();

 // }
  ofile << adjL.size() << " " << num_edges/2 << "\n";
  unsigned long long i = 0;
  unsigned long long numEdges = 0;
  for(std::map<unsigned long long, std::vector<unsigned long long> >::const_iterator it = adjL.begin(); it != adjL.end(); ++it){
    // std::cout <<"first " << it->first <<"\t i: " << i <<"\n";
    // ofile << it->first <<"  ";
   // if(it->first != i){
   //   ofile << " " << 1 << "\n";
      i++;
    // }
     numEdges += it->second.size();
     for (std::vector<unsigned long long>::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
         //if(*vit <= adjL.size())
            ofile << " " << *vit << " ";
    // std::cout <<"Edge " << *vit << "\tsize "<< it->second.size() << "\n";
    //  sum++;
     }
    ofile << "\n";
    

 /*   if(i > num_vertices)
      break;
    ++i;
 */ }
  //  sum = keys.size();
  fprintf(stderr,"\nFile Name: %s , type: %s ",fileName.c_str(), type.c_str());
  fprintf(stderr,"\n Total Vertices %llu, i: %llu,  num_edges %llu, Edges: %llu ", num_vertices, i, num_edges, numEdges/2);
  ofile.close();
  ifile.close();
  return 0;
}

