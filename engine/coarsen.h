
/*************************************************************************/
/*! This function takes a graph and creates a sequence of coarser graphs.
    It implements the coarsening phase of the multilevel paradigm. 
 */
/*************************************************************************/
graph_t* Partitioner::coarsen(const unsigned tid, graph_t *graph, unsigned CoarsenTo) {

     fprintf(stderr,"\nInside Coarsen\n");
// IdType *maxvwgt;
 IdType i, eqewgts, level=0;

     fprintf(stderr,"\nSetting adjwgt with E: %d\n", graph->e);
 for (eqewgts=1, i=1; i<graph->e; i++) {
     if (graph->adjwgt[0] != graph->adjwgt[i]) {
       eqewgts = 0;
       break;
     }
   }

 /* set the maximum allowed coarsest vertex weight */
//  for (i=0; i<graph->ncon; i++)
//    maxvwgt[i] = 1.5*graph->tvwgt[i]/CoarsenTo;

 /* allocate memory for cmap, if it has not already been done due to
       multiple cuts */
    if (graph->cmap == NULL)
      graph->cmap = (IdType *) malloc((graph->v) * sizeof(IdType));

 do{
     fprintf(stderr,"\nBefore MATCH_RM \n");

    Match_RM(graph);
   
    graph = graph->coarser;
    eqewgts = 0;
    level++;
     fprintf(stderr,"\nAfter MATCH_RM V: %d, E: %d, Cfraction: %d\n", graph->v, graph->e, COARSEN_FRACTION*graph->finer->v);
   
   } while(graph->v > CoarsenTo && 
           graph->v < COARSEN_FRACTION*graph->finer->v && 
           graph->e > graph->v/2);

 fprintf(stderr,"\nCoarsening done\n");
 return graph;
}



/*************************************************************************/
/*! This function finds a matching by randomly selecting one of the 
    unmatched adjacent vertices. 
 */
/**************************************************************************/
IdType Partitioner::Match_RM(graph_t *graph) {
  IdType i, pi, ii, j, k, nvtxs, cnvtxs, maxidx, last_unmatched;
  IdType *xadj, *vwgt, *adjncy, *adjwgt;
  IdType *match, *cmap, *perm;
  IdType nunmatched=0;

     fprintf(stderr,"\nInside Match_RM\n");
  nvtxs  = graph->v;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;
  
  perm = (IdType *) malloc(graph->v * sizeof(IdType));
  match = (IdType *) malloc(graph->v * sizeof(IdType));

     fprintf(stderr,"\nInitialize Match Vertices V: %d\n", graph->v);
  for (IdType z = 0; z < graph->v; ++z) {
       match[z] = UNMATCHED; 
       perm[z] = z;
  }
  
     fprintf(stderr,"\nBefore Shuffling: %d\n", perm[0]);
//  std::random_shuffle(std::begin(perm), std::end(perm));
  std::random_shuffle(&perm[0], &perm[nvtxs-1]);
  
  for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
 
    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      //if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) {
          last_unmatched = max(pi, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }

      if (maxidx != UNMATCHED) {
        cmap[i]  = cmap[maxidx] = cnvtxs++;
        match[i] = maxidx;
        match[maxidx] = i;
      }
    }
  }

 /* match the final unmatched vertices with themselves and reorder the vertices 
     of the coarse graph for memory-friendly contraction */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED) {
      match[i] = i;
      cmap[i]  = cnvtxs++;
    }
    else {
      if (i <= match[i]) 
        cmap[i] = cmap[match[i]] = cnvtxs++;
    }
  }
 
 CreateCoarseGraph(graph, cnvtxs, match);
 
 return cnvtxs;
}

/*************************************************************************/
/*! This function creates the coarser graph. It uses a simple hash-table 
    for identifying the adjacent vertices that get collapsed to the same
    node. The hash-table can have conflicts, which are handled via a
    linear scan. 
 */
/*************************************************************************/
void Partitioner::CreateCoarseGraph(graph_t *graph, IdType cnvtxs, IdType *match)
{
  IdType j, jj, k, kk, l, m, istart, iend, nvtxs, ncon, nedges, cnedges, v, u, mask, dovsize;
  IdType *xadj, *vwgt, *vsize, *adjncy, *adjwgt;
  IdType *cmap, *htable;
  IdType *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt;
  graph_t* cgraph;

     fprintf(stderr,"\nInside CreateCoarseGraph\n");
 //TODO: dovsize = (ctrl->objtype == METIS_OBJTYPE_VOL ? 1 : 0); May not need it
  dovsize = 0;

  mask = HTLENGTH;
  nvtxs = graph->v;
  xadj  = graph->xadj;
  
  vwgt    = graph->vwgt;
  vsize   = graph->vsize;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  cmap    = graph->cmap;

  ncon    = graph->ncon;
  /* Initialize the coarser graph */
  cgraph   = SetupCoarseGraph(graph, cnvtxs, dovsize);
  cxadj    = cgraph->xadj;
  cvwgt    = cgraph->vwgt;
  cvsize   = cgraph->vsize;
  cadjncy  = cgraph->adjncy;
  cadjwgt  = cgraph->adjwgt;

 // htable = iset(min(cnvtxs+1, mask+1), -1, iwspacemalloc(ctrl, mask+1)); 
  htable = (IdType *) malloc((mask+1) * sizeof(IdType));
  
  for(int i=0; i<cnvtxs+1; i++){
      htable[i] = -1;   
  }
  
  cxadj[0] = cnvtxs = cnedges = 0;
  for (v=0; v<nvtxs; v++) {
    if ((u = match[v]) < v)
      continue;

    assert(cmap[v] == cnvtxs);
    assert(cmap[match[v]] == cnvtxs);
   
    cvwgt[cnvtxs] = vwgt[v];
  
    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    nedges = 0;

    istart = xadj[v];
    iend   = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k  = cmap[adjncy[j]];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges]   = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

   if (v != u) { 
      cvwgt[cnvtxs] += vwgt[u];

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      istart = xadj[u];
      iend   = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k  = cmap[adjncy[j]];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk]      = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges]   = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
     
     /* This 2nd check is needed for non-adjacent matchings */
      if (jj >= 0 && jj < nedges && cadjncy[jj] == cnvtxs) { 
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }
 
    /* Zero out the htable */
    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  
    htable[cnvtxs&mask] = -1;

    cnedges         += nedges;
    cxadj[++cnvtxs]  = cnedges;
    cadjncy         += nedges;
    cadjwgt         += nedges;
  }

  cgraph->e = cnedges;

  //no for looping -  ncon will always be 1 in my case
/* for (j=0; j<ncon; j++) {
    IdType temp = *cgraph->vwgt;
    temp +=j;
    cgraph->tvwgt[j]    = isum(cgraph->v, temp, ncon);
//    cgraph->invtvwgt[j] = 1.0/(cgraph->tvwgt[j] > 0 ? cgraph->tvwgt[j] : 1);
  }*/

}

/*************************************************************************/
/*! This function creates and initializes a graph_t data structure */
/*************************************************************************/
graph_t* Partitioner::CreateGraph(void)
{
  graph_t *graph;

 //TODO:: something to always allocate some memory may be one byte
  graph = (graph_t *)malloc(sizeof(graph_t));

  InitGraph(graph);

  return graph;
}

/*************************************************************************/
/*! This function initializes a graph_t data structure */
/*************************************************************************/
void Partitioner::InitGraph(graph_t *graph) 
{
  memset((void *)graph, 0, sizeof(graph_t));

  /* graph size constants */
  graph->v     = -1;
  graph->e    = -1;

/*  graph->mincut    = -1;
  graph->minvol    = -1;
  graph->nbnd      = -1;
*/
  /* memory for the graph structure */
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->vsize     = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->cmap      = NULL;
  graph->label     = NULL;
//  graph->tvwgt     = NULL;
//  graph->invtvwgt  = NULL;

  /* by default these are set to true, but the can be explicitly changed afterwards */
  graph->free_xadj   = 1;
  graph->free_vwgt   = 1;
  graph->free_vsize  = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;

  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->bndptr    = NULL;
  graph->bndind    = NULL;

  /* linked-list structure */
  graph->coarser   = NULL;
  graph->finer     = NULL;
}

/*************************************************************************/
/*! Setup the various arrays for the coarse graph 
 */
/*************************************************************************/
graph_t* Partitioner::SetupCoarseGraph(graph_t *graph, IdType cnvtxs, IdType dovsize)
{
  graph_t *cgraph;

  cgraph = CreateGraph();

  cgraph->v = cnvtxs;

  cgraph->finer  = graph;
  graph->coarser = cgraph;


  /* Allocate memory for the coarser graph */
  cgraph->xadj     = (IdType *) malloc((cnvtxs+1) * sizeof(IdType));
  cgraph->adjncy   = (IdType *) malloc((graph->e) * sizeof(IdType));
  cgraph->adjwgt   = (IdType *) malloc((graph->e) * sizeof(IdType));
  cgraph->vwgt     = (IdType *) malloc((cnvtxs) * sizeof(IdType));
 // cgraph->tvwgt     = (IdType *) malloc((1) * sizeof(IdType));

  if (dovsize)
    cgraph->vsize = (IdType *) malloc((cnvtxs) * sizeof(IdType));

  return cgraph;
}

