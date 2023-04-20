//------------------------------------------------------------------------------
// LAGraph/src/benchmark/tc_demo.c: benchmark for LAGr_TriangleCount
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Usage:  test_tc < matrixmarketfile.mtx
//         test_tc matrixmarketfile.mtx
//         test_tc matrixmarketfile.grb

//  Known triangle counts:
//      kron:       106873365648
//      urand:      5378
//      twitter:    34824916864
//      web:        84907041475
//      road:       438804

#include "LAGraph_demo.h"

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&A) ;                 \
}

char t [256] ;

char *method_name (int method, int sorting)
{
  char *s ;
  switch (method)
  {
    case LAGr_MaskedSpGEMM_AutoMethod:  s = "default (No Mask)           " ; break ;
    case LAGr_MaskedSpGEMM_NoMask:      s = "No Mask: Use no mask" ; break ;
    case LAGr_MaskedSpGEMM_AllZeroMask: s = "All-Zero Mask: Use an all-zero mask" ; break ;
    case LAGr_MaskedSpGEMM_UseMask:     s = "Use Mask: Use the provided mask, same as All-Zero Mask" ; break ;
    default: abort ( ) ;
  }

//  if (sorting == LAGr_MaskedSpGEMM_Descending) sprintf (t, "%s sort: descending degree", s) ;
//  else if (sorting == LAGr_MaskedSpGEMM_Ascending) sprintf (t, "%s ascending degree", s) ;
//  else if (sorting == LAGr_MaskedSpGEMM_AutoSort) sprintf (t, "%s auto-sort", s) ;
//  else sprintf (t, "%s sort: none", s) ;

  if (sorting == 0) sprintf (t, "%s sort: none", s);
  else sprintf (t, "%s sort: no way", s);
  return (t) ;
}


void print_method (FILE *f, int method, int sorting)
{
  fprintf (f, "%s\n", method_name (method, sorting)) ;
}

/// Added by Zhen Peng on 12/28/2022
//------------------------------------------------------------------------------
// readmask: read a GAP problem from a file
//------------------------------------------------------------------------------
#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&A2) ;                \
    GrB_free (&M) ;                 \
    if (f != NULL) fclose (f) ;     \
    f = NULL ;                      \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    LAGraph_Delete (G, NULL) ;      \
    GrB_free (src_nodes) ;          \
}

static int readmask          // returns 0 if successful, -1 if failure
    (
        // output
        LAGraph_Graph *G,           // graph from the file
        GrB_Matrix *src_nodes,      // source nodes
        // inputs
        bool make_symmetric,        // if true, always return G as undirected
        bool remove_self_edges,     // if true, remove self edges
        bool structural,            // if true, return G->A as bool (all true)
        GrB_Type pref,              // if non-NULL, typecast G->A to this type
        bool ensure_positive,       // if true, ensure all entries are > 0
        int argc,                   // input to main test program
        char **argv                 // input to main test program
    )
{

  //--------------------------------------------------------------------------
  // check inputs
  //--------------------------------------------------------------------------

  char msg [LAGRAPH_MSG_LEN] ;
  msg [0] = '\0' ;
  GrB_Matrix A = NULL, A2 = NULL, M = NULL ;
  GrB_Type atype = NULL ;
  FILE *f = NULL ;
  if (G == NULL) CATCH (GrB_NULL_POINTER) ;
  (*G) = NULL ;
  if (src_nodes != NULL) (*src_nodes) = NULL ;
  GrB_Type src_type = NULL;

  //--------------------------------------------------------------------------
  // read in a matrix from a file
  //--------------------------------------------------------------------------

  double t_read = LAGraph_WallClockTime ( ) ;

  if (argc > 2)
  {
    // Usage:
    //      ./test_whatever matrixfile.mtx [sources.mtx]
    //      ./test_whatever matrixfile.grb [sources.mtx]

    // read in the file in Matrix Market format from the input file
    char *filename = argv [2] ;
//    char *filename = argv [1] ;
    printf ("matrix: %s\n", filename) ;

    // find the filename extension
    size_t len = strlen (filename) ;
    char *ext = NULL ;
    for (int k = len-1 ; k >= 0 ; k--)
    {
      if (filename [k] == '.')
      {
        ext = filename + k ;
        printf ("[%s]\n", ext) ;
        break ;
      }
    }

    bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

    if (is_binary)
    {
      printf ("Reading binary file: %s\n", filename) ;
      f = fopen (filename, "r") ;
      if (f == NULL)
      {
        printf ("Binary file not found: [%s]\n", filename) ;
        exit (1) ;
      }
      if (binread (&A, f) < 0) CATCH (-1001) ;    // file I/O error
      fclose (f) ;
      f = NULL ;
    }
    else
    {
      printf ("Reading matrix market file: %s\n", filename) ;
      f = fopen (filename, "r") ;
      if (f == NULL)
      {
        printf ("Matrix market file not found: [%s]\n", filename) ;
        exit (1) ;
      }
      int result = LAGraph_MMRead (&A, f, msg) ;
      if (result != GrB_SUCCESS)
      {
        printf ("LAGraph_MMRead failed to read matrix: %s\n",
                filename) ;
        printf ("result: %d msg: %s\n", result, msg) ;
      }
      LAGRAPH_TRY (result) ;
      fclose (f) ;
      f = NULL ;
    }

    // read in source nodes in Matrix Market format from the input file
    if (argc > 2 && src_nodes != NULL)
    {
      // do not read in the file if the name starts with "-"
      filename = argv [2] ;
      if (filename [0] != '-')
      {
        printf ("sources: %s\n", filename) ;
        f = fopen (filename, "r") ;
        if (f == NULL)
        {
          printf ("Source node file not found: [%s]\n", filename) ;
          exit (1) ;
        }
        int result = LAGraph_MMRead (src_nodes, f, msg) ;
        if (result != GrB_SUCCESS)
        {
          printf ("LAGraph_MMRead failed to read source nodes"
                  " from: %s\n", filename) ;
          printf ("result: %d msg: %s\n", result, msg) ;
        }
        LAGRAPH_TRY (result) ;
        fclose (f) ;
        f = NULL ;
      }
    }
  }
  else
  {

    // Usage:  ./test_whatever < matrixfile.mtx
    printf ("matrix: from stdin\n") ;

    // read in the file in Matrix Market format from stdin
    int result = LAGraph_MMRead (&A, stdin, msg) ;
    if (result != GrB_SUCCESS)
    {
      printf ("LAGraph_MMRead failed to read: stdin\n") ;
      printf ("result: %d msg: %s\n", result, msg) ;
    }
    LAGRAPH_TRY (result) ;
  }

  //--------------------------------------------------------------------------
  // get the size of the problem.
  //--------------------------------------------------------------------------

  GrB_Index nrows, ncols ;
  GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
  GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
  GrB_Index n = nrows ;
  if (nrows != ncols) CATCH (GrB_DIMENSION_MISMATCH) ;    // A must be square

  //--------------------------------------------------------------------------
  // typecast, if requested
  //--------------------------------------------------------------------------

  GRB_TRY (GxB_Matrix_type (&atype, A)) ;

  if (structural)
  {
    // convert to boolean, with all entries true
    atype = GrB_BOOL ;
    LAGRAPH_TRY (LAGraph_Matrix_Structure (&A2, A, msg)) ;
  }
  else if (pref != NULL && atype != pref)
  {
    // convert to the requested type
    GRB_TRY (GrB_Matrix_new (&A2, pref, n, n)) ;
    atype = pref ;

    GrB_UnaryOp op = NULL ;
    if      (pref == GrB_BOOL  ) op = GrB_IDENTITY_BOOL ;
    else if (pref == GrB_INT8  ) op = GrB_IDENTITY_INT8 ;
    else if (pref == GrB_INT16 ) op = GrB_IDENTITY_INT16 ;
    else if (pref == GrB_INT32 ) op = GrB_IDENTITY_INT32 ;
    else if (pref == GrB_INT64 ) op = GrB_IDENTITY_INT64 ;
    else if (pref == GrB_UINT8 ) op = GrB_IDENTITY_UINT8 ;
    else if (pref == GrB_UINT16) op = GrB_IDENTITY_UINT16 ;
    else if (pref == GrB_UINT32) op = GrB_IDENTITY_UINT32 ;
    else if (pref == GrB_UINT64) op = GrB_IDENTITY_UINT64 ;
    else if (pref == GrB_FP32  ) op = GrB_IDENTITY_FP32 ;
    else if (pref == GrB_FP64  ) op = GrB_IDENTITY_FP64 ;
#if 0
      else if (pref == GxB_FC32  ) op = GxB_IDENTITY_FC32 ;
        else if (pref == GxB_FC64  ) op = GxB_IDENTITY_FC64 ;
#endif
    else CATCH (GrB_NOT_IMPLEMENTED) ;    // unsupported type

    GRB_TRY (GrB_apply (A2, NULL, NULL, op, A, NULL)) ;
  }

  if (A2 != NULL)
  {
    GrB_free (&A) ;
    A = A2 ;
    A2 = NULL ;
    GRB_TRY (GrB_wait (A, GrB_MATERIALIZE)) ;
  }

  //--------------------------------------------------------------------------
  // construct the initial graph
  //--------------------------------------------------------------------------

  bool A_is_symmetric =
      (n == 134217726 ||  // HACK for kron
       n == 134217728) ;  // HACK for urand

  LAGraph_Kind G_kind = A_is_symmetric ?  LAGraph_ADJACENCY_UNDIRECTED :
                        LAGraph_ADJACENCY_DIRECTED ;
  LAGRAPH_TRY (LAGraph_New (G, &A, G_kind, msg)) ;
  // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

  //--------------------------------------------------------------------------
  // remove self-edges, if requested
  //--------------------------------------------------------------------------

  if (remove_self_edges)
  {
    LAGRAPH_TRY (LAGraph_DeleteSelfEdges (*G, msg)) ;
  }
  // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

  //--------------------------------------------------------------------------
  // ensure all entries are > 0, if requested
  //--------------------------------------------------------------------------

  if (!structural && ensure_positive)
  {
    // drop explicit zeros (FUTURE: make this a utility function)
    GrB_IndexUnaryOp idxop = NULL ;
    if      (atype == GrB_BOOL  ) idxop = GrB_VALUENE_BOOL ;
    else if (atype == GrB_INT8  ) idxop = GrB_VALUENE_INT8 ;
    else if (atype == GrB_INT16 ) idxop = GrB_VALUENE_INT16 ;
    else if (atype == GrB_INT32 ) idxop = GrB_VALUENE_INT32 ;
    else if (atype == GrB_INT64 ) idxop = GrB_VALUENE_INT64 ;
    else if (atype == GrB_UINT8 ) idxop = GrB_VALUENE_UINT8 ;
    else if (atype == GrB_UINT16) idxop = GrB_VALUENE_UINT16 ;
    else if (atype == GrB_UINT32) idxop = GrB_VALUENE_UINT32 ;
    else if (atype == GrB_UINT64) idxop = GrB_VALUENE_UINT64 ;
    else if (atype == GrB_FP32  ) idxop = GrB_VALUENE_FP32 ;
    else if (atype == GrB_FP64  ) idxop = GrB_VALUENE_FP64 ;
#if 0
    else if (atype == GxB_FC32  ) idxop = GxB_VALUENE_FC32 ;
        else if (atype == GxB_FC64  ) idxop = GxB_VALUENE_FC64 ;
#endif
    if (idxop != NULL)
    {
      GRB_TRY (GrB_select ((*G)->A, NULL, NULL, idxop, (*G)->A, 0, NULL));
    }

    // A = abs (A)
    GrB_UnaryOp op = NULL ;
    if      (atype == GrB_INT8  ) op = GrB_ABS_INT8 ;
    else if (atype == GrB_INT16 ) op = GrB_ABS_INT16 ;
    else if (atype == GrB_INT32 ) op = GrB_ABS_INT32 ;
    else if (atype == GrB_INT64 ) op = GrB_ABS_INT64 ;
    else if (atype == GrB_FP32  ) op = GrB_ABS_FP32 ;
    else if (atype == GrB_FP64  ) op = GrB_ABS_FP64 ;
#if 0
    else if (atype == GxB_FC32  ) op = GxB_ABS_FC32 ;
        else if (atype == GxB_FC64  ) op = GxB_ABS_FC64 ;
#endif
    if (op != NULL)
    {
      GRB_TRY (GrB_apply ((*G)->A, NULL, NULL, op, (*G)->A, NULL)) ;
    }
  }

  //--------------------------------------------------------------------------
  // determine the graph properies
  //--------------------------------------------------------------------------

  // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

  if (!A_is_symmetric)
  {
    // compute G->AT and determine if A has a symmetric structure
    LAGRAPH_TRY (LAGraph_Cached_IsSymmetricStructure (*G, msg)) ;
    if (((*G)->is_symmetric_structure == LAGraph_TRUE) && structural)
    {
      // if G->A has a symmetric structure, declare the graph undirected
      // and free G->AT since it isn't needed.
      (*G)->kind = LAGraph_ADJACENCY_UNDIRECTED ;
      GRB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
    }
    else if (make_symmetric)
    {
      // make sure G->A is symmetric
      bool sym ;
      LAGRAPH_TRY (LAGraph_Matrix_IsEqual (&sym, (*G)->A, (*G)->AT, msg));
      if (!sym)
      {
        printf ("forcing G-> to be symmetric (via A = A+A')\n") ;
        GrB_BinaryOp op = NULL ;
        GrB_Type type ;
        if      (atype == GrB_BOOL  ) op = GrB_LOR ;
        else if (atype == GrB_INT8  ) op = GrB_PLUS_INT8 ;
        else if (atype == GrB_INT16 ) op = GrB_PLUS_INT16 ;
        else if (atype == GrB_INT32 ) op = GrB_PLUS_INT32 ;
        else if (atype == GrB_INT64 ) op = GrB_PLUS_INT64 ;
        else if (atype == GrB_UINT8 ) op = GrB_PLUS_UINT8 ;
        else if (atype == GrB_UINT16) op = GrB_PLUS_UINT16 ;
        else if (atype == GrB_UINT32) op = GrB_PLUS_UINT32 ;
        else if (atype == GrB_UINT64) op = GrB_PLUS_UINT64 ;
        else if (atype == GrB_FP32  ) op = GrB_PLUS_FP32 ;
        else if (atype == GrB_FP64  ) op = GrB_PLUS_FP64 ;
#if 0
          else if (type == GxB_FC32  ) op = GxB_PLUS_FC32 ;
                else if (type == GxB_FC64  ) op = GxB_PLUS_FC64 ;
#endif
        else CATCH (GrB_NOT_IMPLEMENTED) ;    // unknown type
        GRB_TRY (GrB_eWiseAdd ((*G)->A, NULL, NULL, op,
                               (*G)->A, (*G)->AT, NULL)) ;
      }
      // G->AT is not required
      GRB_TRY (GrB_Matrix_free (&((*G)->AT))) ;
      (*G)->kind = LAGraph_ADJACENCY_UNDIRECTED ;
      (*G)->is_symmetric_structure = LAGraph_TRUE ;
    }
  }
  // LAGRAPH_TRY (LAGraph_Graph_Print (*G, 2, stdout, msg)) ;

  //--------------------------------------------------------------------------
  // generate 64 random source nodes, if requested but not provided on input
  //--------------------------------------------------------------------------

#define NSOURCES 64

  if (src_nodes != NULL && (*src_nodes == NULL))
  {
    src_type = GrB_UINT64;
    GRB_TRY (GrB_Matrix_new (src_nodes, src_type, NSOURCES, 1)) ;
    srand (1) ;
    for (int k = 0 ; k < NSOURCES ; k++)
    {
      uint64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
      // src_nodes [k] = i
      GRB_TRY (GrB_Matrix_setElement (*src_nodes, i, k, 0)) ;
    }
  }

  if (src_nodes != NULL)
  {
    GRB_TRY (GrB_wait (*src_nodes, GrB_MATERIALIZE)) ;
  }

  //--------------------------------------------------------------------------
  // free workspace, print a summary of the graph, and return result
  //--------------------------------------------------------------------------

  t_read = LAGraph_WallClockTime ( ) - t_read ;
  printf ("read time: %g\n", t_read) ;

  LG_FREE_WORK ;
  // LAGRAPH_TRY (LAGraph_Graph_Print (*G, LAGraph_SHORT, stdout, msg)) ;
  return (GrB_SUCCESS) ;
}

#undef  LG_FREE_WORK
#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

int main (int argc, char **argv)
{

  //--------------------------------------------------------------------------
  // initialize LAGraph and GraphBLAS
  //--------------------------------------------------------------------------

  char msg [LAGRAPH_MSG_LEN] ;

  GrB_Matrix A = NULL ;
  LAGraph_Graph G = NULL ;

  // start GraphBLAS and LAGraph
  bool burble = false ;
  demo_init (burble) ;

//  int ntrials = 3 ;
  int ntrials = 10 ;
  // ntrials = 1 ;        // HACK
  printf ("# of trials: %d\n", ntrials) ;

  int nt = NTHREAD_LIST ;
  int Nthreads [20] = { 0, THREAD_LIST } ;

  /// Changed by Zhen Peng on 12/28/2022
  int nthreads_max, nthreads_outer, nthreads_inner ;
  LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
  nthreads_max = nthreads_outer * nthreads_inner ;
//  int nthreads_max = 1;
  if (Nthreads [1] == 0)
  {
    // create thread list automatically
    Nthreads [1] = nthreads_max ;
    for (int t = 2 ; t <= nt ; t++)
    {
      Nthreads [t] = Nthreads [t-1] / 2 ;
      if (Nthreads [t] == 0) nt = t-1 ;
    }
  }
  printf ("threads to test: ") ;
  for (int t = 1 ; t <= nt ; t++)
  {
    int nthreads = Nthreads [t] ;
    if (nthreads > nthreads_max) continue ;
    printf (" %d", nthreads) ;
  }
  printf ("\n") ;

  //--------------------------------------------------------------------------
  // read in the graph
  //--------------------------------------------------------------------------

  char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
//  LAGRAPH_TRY (readproblem (&G, NULL,
//                            true, true, true, NULL, false, argc, argv)) ;
  LAGRAPH_TRY (readproblem (&G, NULL,
                            false, false, false, GrB_FP64, false, argc, argv)) ;
  LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

  // determine the cached out degree property
  LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

  GrB_Index n, nvals ;
  GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
  GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

//  LAGraph_Graph M;
////  LAGRAPH_TRY (readmask (&M, NULL,
////                            true, true, true, NULL, false, argc, argv)) ;
//  LAGRAPH_TRY (readmask (&M, NULL,
//                            false, false, false, GrB_FP64, false, argc, argv)) ;
//  LAGRAPH_TRY (LAGraph_Graph_Print (M, LAGraph_SHORT, stdout, msg)) ;

  //--------------------------------------------------------------------------
  // triangle counting
  //--------------------------------------------------------------------------

  GrB_Index ntriangles, ntsimple = 0 ;

#if 0
  // check # of triangles
    double tsimple = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LG_check_tri (&ntsimple, G, NULL)) ;
    tsimple = LAGraph_WallClockTime ( ) - tsimple ;
    printf ("# of triangles: %" PRId64 " slow time: %g sec\n",
        ntsimple, tsimple) ;
#endif

//  // warmup for more accurate timing, and also print # of triangles
//  double ttot = LAGraph_WallClockTime ( ) ;
//  printf ("\nwarmup method: ") ;
////  int presort = LAGr_TriangleCount_AutoSort ; // = 0 (auto selection)
//  int presort = 0;
//  print_method (stdout, 0, presort) ;
//
//  // warmup method:
//  // LAGr_TriangleCount_Sandia_ULT: sum (sum ((U * L') .* U))
//  LAGr_MaskedSpGEMM_Method method = LAGr_MaskedSpGEMM_NoMask ;
////  LAGRAPH_TRY (LAGr_MaskedSpGEMM (&ntriangles, G, &method, &presort, msg)) ;
//  LAGRAPH_TRY (LAGr_MaskedSpGEMM (&ntriangles, G, M, &method, msg)) ;
//  printf ("# of triangles: %" PRIu64 "\n", ntriangles) ;
//  print_method (stdout, 0, presort) ;
//  ttot = LAGraph_WallClockTime ( ) - ttot ;
//  printf ("nthreads: %3d time: %12.6f rate: %6.2f (NoMask, one trial)\n",
//          nthreads_max, ttot, 1e-6 * nvals / ttot) ;

  int presort = 0;
  ntriangles = -1;

#if 0
  if (ntriangles != ntsimple)
    {
        printf ("wrong # triangles: %g %g\n", (double) ntriangles,
            (double) ntsimple) ;
        abort ( ) ;
    }
#endif

  double t_best = INFINITY ;
  int method_best = -1 ;
  int nthreads_best = -1 ;
//  int sorting_best = 0 ;

  // kron: input graph: nodes: 134217726 edges: 4223264644
  // fails on methods 3 and 4.

  // just try methods 5 and 6
  // for (int method = 5 ; method <= 6 ; method++)

  // try all methods 3 to 5
//    for (int method = 3 ; method <= 5 ; method++)
//    for (int method = 1 ; method <= 6 ; method++)
//    for (int method = 1; method <= 7 ; method++)
//  for (int method = 7; method >= 1 ; method--)
//  for (int method = 1; method <= 3; ++method)

/// No Mask
  { /// No Mask
    LAGr_MaskedSpGEMM_Method method = LAGr_MaskedSpGEMM_NoMask;
    // for (int sorting = -1 ; sorting <= 2 ; sorting++)

//    int sorting = LAGr_TriangleCount_AutoSort ; // just use auto-sort
    {
//      printf ("\nMethod: ") ; print_method (stdout, method, sorting) ;
//            if (n == 134217726 && method < 5)
//            {
//                printf ("kron fails on method %d; skipped\n", method) ;
//                continue ;
//            }
//            if (n != 134217728 && method < 5)
//            {
//                printf ("all but urand is slow with method %d: skipped\n",
//                        method) ;
//                continue ;
//            }

      for (int t = 1 ; t <= nt ; t++)
      {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;
        GrB_Index nt2 ;
        double ttot = 0, ttrial [100] ;
//        LAGr_TriangleCount_Presort p ;
        LAGr_MaskedSpGEMM_Method m ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
//          {
//            printf("trial: %d\n", trial); //test
//          }
          double tt = LAGraph_WallClockTime ( ) ;
          m = method ;
//          p = sorting ;
          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, G, &m, msg));
//          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, M, &m, msg));
//          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, &m, &p, msg));
          ttrial [trial] = LAGraph_WallClockTime ( ) - tt ;
          ttot += ttrial [trial] ;
          printf ("trial %2d: %12.6f sec rate %6.2f  # triangles: "
                  "%g\n", trial, ttrial [trial],
                  1e-6 * nvals / ttrial [trial], (double) nt2) ;
        }
        ttot = ttot / ntrials ;
        printf ("nthreads: %3d time: %12.6f rate: %6.2f", nthreads,
                ttot, 1e-6 * nvals / ttot) ;
        printf ("   # of triangles: %" PRId64 " presort: %d\n",
                ntriangles, (int) 0) ;
        if (nt2 != ntriangles)
        {
          printf ("Test failure!\n") ;
          abort ( ) ;
        }
        fprintf (stderr, "\nMethod used: ") ;
        print_method (stderr, m, presort) ;
        fprintf (stderr, "Avg: MaskedSpGEEMM method%d %3d: %10.3f sec: %s\n",
                 method, nthreads, ttot, matrix_name) ;
//        fprintf (stderr, "Avg: TC method%d.%d %3d: %10.3f sec: %s\n",
//                 method, sorting, nthreads, ttot, matrix_name) ;

        if (ttot < t_best)
        {
          t_best = ttot ;
          method_best = method ;
          nthreads_best = nthreads ;
//          sorting_best = sorting ;
        }
      }
    }
  }
/// End No Mask

///// Use Mask
//  { /// Use Mask
//    LAGr_MaskedSpGEMM_Method method = LAGr_MaskedSpGEMM_UseMask;
//    // for (int sorting = -1 ; sorting <= 2 ; sorting++)
//
////    int sorting = LAGr_TriangleCount_AutoSort ; // just use auto-sort
//    {
////      printf ("\nMethod: ") ; print_method (stdout, method, sorting) ;
////            if (n == 134217726 && method < 5)
////            {
////                printf ("kron fails on method %d; skipped\n", method) ;
////                continue ;
////            }
////            if (n != 134217728 && method < 5)
////            {
////                printf ("all but urand is slow with method %d: skipped\n",
////                        method) ;
////                continue ;
////            }
//
//      for (int t = 1 ; t <= nt ; t++)
//      {
//        int nthreads = Nthreads [t] ;
//        if (nthreads > nthreads_max) continue ;
//        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;
//        GrB_Index nt2 ;
//        double ttot = 0, ttrial [100] ;
////        LAGr_TriangleCount_Presort p ;
//        LAGr_MaskedSpGEMM_Method m ;
//        for (int trial = 0 ; trial < ntrials ; trial++)
//        {
////          {
////            printf("trial: %d\n", trial); //test
////          }
//          double tt = LAGraph_WallClockTime ( ) ;
//          m = method ;
////          p = sorting ;
//          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, G, &m, msg));
////          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, M, &m, msg));
////          LAGRAPH_TRY(LAGr_MaskedSpGEMM (&nt2, G, &m, &p, msg));
//          ttrial [trial] = LAGraph_WallClockTime ( ) - tt ;
//          ttot += ttrial [trial] ;
//          printf ("trial %2d: %12.6f sec rate %6.2f  # triangles: "
//                  "%g\n", trial, ttrial [trial],
//                  1e-6 * nvals / ttrial [trial], (double) nt2) ;
//        }
//        ttot = ttot / ntrials ;
//        printf ("nthreads: %3d time: %12.6f rate: %6.2f", nthreads,
//                ttot, 1e-6 * nvals / ttot) ;
//        printf ("   # of triangles: %" PRId64 " presort: %d\n",
//                ntriangles, (int) 0) ;
//        if (nt2 != ntriangles)
//        {
//          printf ("Test failure!\n") ;
//          abort ( ) ;
//        }
//        fprintf (stderr, "\nMethod used: ") ;
//        print_method (stderr, m, presort) ;
//        fprintf (stderr, "Avg: MaskedSpGEEMM method%d %3d: %10.3f sec: %s\n",
//                 method, nthreads, ttot, matrix_name) ;
////        fprintf (stderr, "Avg: TC method%d.%d %3d: %10.3f sec: %s\n",
////                 method, sorting, nthreads, ttot, matrix_name) ;
//
//        if (ttot < t_best)
//        {
//          t_best = ttot ;
//          method_best = method ;
//          nthreads_best = nthreads ;
////          sorting_best = sorting ;
//        }
//      }
//    }
//  }
///// End Use Mask

  printf ("\nBest method: ") ;
  print_method (stdout, method_best, presort) ;
  printf ("nthreads: %3d time: %12.6f rate: %6.2f\n",
          nthreads_best, t_best, 1e-6 * nvals / t_best) ;
  LG_FREE_ALL ;
  LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
  return (GrB_SUCCESS) ;
}

