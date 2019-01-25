#ifndef __GRAPH_H_
#define __GRAPH_H_

#include <climits>

#define MAX_VERTEX_VALUE (ULLONG_MAX)

//------------------------------------------------------------------
// Warning: for OnePhaseFileIO, do not use STL structures with variable sizes
//Structure to represent a graph

typedef struct graph_t {
  IdType v, e;	/* The # of vertices and edges in the graph */
  IdType ncon;		/* The # of constrains */ 
  IdType *xadj;		/* Pointers to the locally stored vertices */
  IdType *vwgt;		/* Vertex weights */
  IdType *vsize;		/* Vertex sizes for min-volume formulation */
  IdType *adjncy;        /* Array that stores the adjacency lists of nvtxs */
  IdType *adjwgt;        /* Array that stores the weights of the adjacency lists */
//  IdType *tvwgt;         /* The sum of the vertex weights in the graph */
//  real_t *invtvwgt;     /* The inverse of the sum of the vertex weights in the graph */


  IdType *cmap;
  /* These are to keep track control if the corresponding fields correspond to
     application or library memory */
  int free_xadj, free_vwgt, free_vsize, free_adjncy, free_adjwgt;

  IdType *label;

//  idx_t *cmap;

  /* Partition parameters */
//  idx_t mincut, minvol;
  IdType *where, *pwgts;
//  idx_t nbnd;
  IdType *bndptr, *bndind;

  /* Bisection refinement parameters */
//  idx_t *id, *ed;

  /* K-way refinement parameters */
//  ckrinfo_t *ckrinfo;   /*!< The per-vertex cut-based refinement info */
//  vkrinfo_t *vkrinfo;   /*!< The per-vertex volume-based refinement info */

  /* Node refinement information */
//  nrinfo_t *nrinfo;

  struct graph_t *coarser, *finer;
} graph_t;

#endif // __GRAPH_H_
