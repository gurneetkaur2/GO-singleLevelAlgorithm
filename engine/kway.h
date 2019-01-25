/*************************************************************************/
/*! This function allocates memory for the k-way cut-based refinement */
/*************************************************************************/
void Partitioner::AllocateKWayPartitionMemory(graph_t *graph)
{

  graph->pwgts  = (IdType *) malloc((nparts*graph->ncon) * sizeof(IdType));
  graph->where  = (IdType *) malloc((graph->v) * sizeof(IdType));
  graph->bndptr = (IdType *) malloc((graph->v) * sizeof(IdType));
  graph->bndind = (IdType *) malloc((graph->v) * sizeof(IdType));

}

/*************************************************************************/
/*! This function computes the initial k-way partitioning 
*/
/*************************************************************************/
void Partitioner::InitKWayPartitioning(graph_t *graph) {
	IdType curobj=0;
//	double *ubvec=NULL;
//	int status;

//	ubvec = (double *) malloc((graph->ncon) * sizeof(double));

//TODO: not sure if needed
//	for (i=0; i<graph->ncon; i++) 
//          ubvec[i] = (real_t)pow(ctrl->ubfactors[i], 1.0/log(ctrl->nparts));

//	status = PartGraphRecursive(&graph->v, &graph->ncon, 
	 PartGraphRecursive(&graph->v, &graph->ncon, 
                   graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                   graph->adjwgt, &nparts, &curobj, graph->where);

  /*      if(status != 1){
		fprintf(stderr,"\nFailed during initial partitioning\n");
//		break;
	}*/

//  	free((void **)&ubvec);
}

/*************************************************************************/

void Partitioner::PartGraphRecursive(IdType* nvtxs, IdType* ncon, IdType* xadj, IdType* adjncy, IdType* vwgt, IdType* vsize, IdType* adjwgt, unsigned int* nparts, IdType* objval, IdType* part){

     graph_t *graph;

     graph = SetupGraph(*nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

     *objval = MlevelRecursiveBisection(graph, *nparts, part, 0);

}

/*************************************************************************/
/*! This function is the top-level driver of the recursive bisection 
    routine. */
/*************************************************************************/
IdType Partitioner::MlevelRecursiveBisection(graph_t *graph, IdType nparts,                   IdType *part, IdType fpart) {

    IdType i, j, nvtxs, ncon, objval;
    IdType *label, *where;
    graph_t *lgraph, *rgraph;

    if ((nvtxs = graph->v) == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n"
           "\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

   ncon = graph->ncon;

   /* perform the bisection */
  objval = MultilevelBisect(graph);

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

//  if (nparts > 2) 
//    SplitGraphPart(graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  FreeGraph(&graph);

  /* Do the recursive call */
  if (nparts > 3) {
    objval += MlevelRecursiveBisection(lgraph, (nparts>>1), part, fpart);
    objval += MlevelRecursiveBisection(rgraph, nparts-(nparts>>1), part, fpart+(nparts>>1));
  }
  else if (nparts == 3) {
    FreeGraph(&lgraph);
    objval += MlevelRecursiveBisection(rgraph, nparts-(nparts>>1), part, 
               fpart+(nparts>>1));
  }


  return objval;

}

/*************************************************************************/
/*! This function performs a multilevel bisection */
/*************************************************************************/
IdType Partitioner::MultilevelBisect(graph_t *graph)
{
  IdType i, niparts, bestobj=0, curobj=0, *bestwhere=NULL;
  graph_t *cgraph;
  real_t bestbal=0.0, curbal=0.0;

  return bestobj;
}

/*************************************************************************/
/*! This function splits a graph into two based on its bisection */
/*************************************************************************/
void SplitGraphPart(graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph)
{



}
/*************************************************************************/
/*! This function sets up the graph from the user input */
/*************************************************************************/
graph_t* Partitioner::SetupGraph(IdType nvtxs, IdType ncon, IdType *xadj, 
             IdType *adjncy, IdType *vwgt, IdType *vsize, IdType *adjwgt) 
{
  IdType i, j, k, sum;
  double *nvwgt;
  graph_t *graph;

  /* allocate the graph and fill in the fields */
  graph = CreateGraph();

  graph->v  = nvtxs;
  graph->e = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj      = xadj;
  graph->free_xadj = 0;

  graph->adjncy      = adjncy;
  graph->free_adjncy = 0;


  /* setup the vertex weights */
  if (vwgt) {
    graph->vwgt      = vwgt;
    graph->free_vwgt = 0;
  }
  else {
    vwgt = graph->vwgt = (IdType *)malloc((ncon*nvtxs) * sizeof(IdType));
    for(int i = 0; i < ncon*nvtxs; i++)
                graph->vwgt[i] = 1;
	
  }


//  if (ctrl->objtype == METIS_OBJTYPE_VOL) { 
    /* Setup the vsize */
    if (vsize) {
      graph->vsize      = vsize;
      graph->free_vsize = 0;
    }
    else {
      vsize = graph->vsize = (IdType *)malloc((nvtxs) * sizeof(IdType)); 
      for(int i = 0; i < nvtxs; i++)
                graph->vsize[i] = 1;
    }

    /* Allocate memory for edge weights and initialize them to the sum of the vsize */
/*    adjwgt = graph->adjwgt = (IdType *) malloc(graph->e * sizeof(IdType));
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        adjwgt[j] = 1+vsize[i]+vsize[adjncy[j]];
    }*/
//  }
//  else { /* For edgecut minimization */
    /* setup the edge weights */
    if (adjwgt) {
      graph->adjwgt      = adjwgt;
      graph->free_adjwgt = 0;
    }
    else {
      adjwgt = graph->adjwgt = (IdType *) malloc(graph->e * sizeof(IdType));
      for(int i = 0; i < graph->e; i++)
                graph->adjwgt[i] = 1;	
    }
//  }
   SetupGraph_label(graph);


  return graph;
}

/*************************************************************************/
/*! This function deallocates any memory stored in a graph */
/*************************************************************************/
void Partitioner::FreeGraph(graph_t **r_graph) 
{
  graph_t *graph;

  graph = *r_graph;

  /* free graph structure */
  if (graph->free_xadj)
    free((void **)&graph->xadj);
  if (graph->free_vwgt)
    free((void **)&graph->vwgt);
  if (graph->free_vsize)
    free((void **)&graph->vsize);
  if (graph->free_adjncy)
    free((void **)&graph->adjncy);
  if (graph->free_adjwgt)
    free((void **)&graph->adjwgt);
    
  /* free partition/refinement structure */
  free((void **)&graph->where);
  free((void **)&graph->pwgts);
  free((void **)&graph->bndptr); 
  free((void **)&graph->bndind);

  free((void **)&graph->label);
  free((void **)&graph->cmap);
  free((void **)&graph);

  *r_graph = NULL;
}


/*************************************************************************/
/*! Set's up the label info */
/*************************************************************************/
void Partitioner::SetupGraph_label(graph_t *graph)
{
  IdType i;

  if (graph->label == NULL)
    graph->label = (IdType *)malloc(graph->v * sizeof(IdType));

  for (i=0; i<graph->v; i++)
    graph->label[i] = i;
}

