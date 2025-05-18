#include<iostream>
#include<fstream>
#include <sstream>
#include<map>
#include<vector>
#include<cassert>

using namespace std;

int main(int argc, char** argv)
{
  int sum = 0;
  if (argc != 3)
  {
    cout << "Usage: " << argv[0] << " <fileName> <graph-type> " << std::endl;
    return 0;
  }

  std::string fileName = "";
  fileName = argv[1];
  std::string type = "";
  type = argv[2];

  std::string outFile = fileName.c_str();
  outFile += "_done";

  ifstream ifile(fileName.c_str());
  assert(ifile.is_open());
  ofstream ofile(outFile.c_str());
  assert(ofile.is_open());

  std::map<int, int> keys;
  std::string line;
  ifile.seekg(std::ios::beg);

  while(std::getline(ifile, line)){
    unsigned from, to;
    std::stringstream iss(line);
    if (type == "adjlist" ){
    while(iss >> to) {
      std::map<int, int>::iterator it_to = keys.find(to);

      if(it_to != keys.end()){
        keys[to] = it_to->second + 1;
      }
      else {
        keys[to] = 1;
        sum++;
      }

     }
    }
    if (type == "edgelist" ){
    while(iss >> to >> from) {
      std::map<int, int>::iterator it_to = keys.find(to);

      if(it_to != keys.end()){
        keys[to] = it_to->second + 1;
      }
      else {
        keys[to] = 1;
        sum++;
      }

     }
    }
  }
  fprintf(stderr,"\nFile Name: %s , Keys %d , type: %s ",fileName.c_str(), keys.size(), type.c_str());
  fprintf(stderr,"\n Total Vertices %d ", sum);
  ofile.close();
  ifile.close();
  return 0;
}

