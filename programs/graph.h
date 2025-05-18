#ifndef __GRAPH_H_
#define __GRAPH_H_

#include <climits>
#include "edgeList.pb.h"
#define MAX_VERTEX_VALUE (ULLONG_MAX)

//------------------------------------------------------------------
// Warning: for OnePhaseFileIO, do not use STL structures with variable sizes
//Structure to represent a graph

struct Edge {
  unsigned long long src;
  unsigned long long dst;

} ;

   bool operator==(const Edge &lhs, const Edge &rhs) {
        return lhs.dst == rhs.dst && lhs.src == rhs.src;
    }


#endif // __GRAPH_H_
