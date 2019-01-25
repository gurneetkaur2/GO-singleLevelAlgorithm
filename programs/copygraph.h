#ifndef __GRAPH_H_
#define __GRAPH_H_

#include <climits>

#define MAX_VERTEX_VALUE (ULLONG_MAX)

//------------------------------------------------------------------
// Warning: for OnePhaseFileIO, do not use STL structures with variable sizes
//Structure to represent a graph
class graph_t {

	public:
 
        IdType v;
        IdType e;
        IdType orgM;

        IdType *num_edges;
        IdType *adj;
        IdType *edgeWeight;
        IdType *vertexWeight;

 };

//Structure to represent an edge
struct Edge {
   IdType u;
   IdType v;

   Edge(IdType inU, IdType inV) : u(inU), v(inV) {}

   IdType either() const { return u; }
   IdType other(IdType inU) const {
        if (inU == u )
                return v;
        else
                return u;
   }

   bool operator < (const Edge right) const {
        if( u < right.u )
                return true;
        else if( v < right.v)
                return true;
        return false;
   }
};
/*typedef struct vertex {
  long double rank; // Distance of this node from source vertex.
  unsigned numNeighbors;
} Vertex;

//------------------------------------------------------------------
typedef GraphPart<Vertex, EdgeList> RunTime;
*/
#endif // __GRAPH_H_
