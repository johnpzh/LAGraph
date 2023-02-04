//------------------------------------------------------------------------------
// LAGr_TriangleCount: Triangle counting using various methods
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Count the number of triangles in a graph,

// This is an Advanced algorithm (G->nself_edges, G->out_degree,
// G->is_symmetric_structure are required).

// Given a symmetric graph A with no-self edges, LAGr_TriangleCount counts the
// number of triangles in the graph.  A triangle is a clique of size three,
// that is, 3 nodes that are all pairwise connected.

// One of 6 methods are used, defined below where L and U are the strictly
// lower and strictly upper triangular parts of the symmetrix matrix A,
// respectively.  Each method computes the same result, ntri:

//  0:  default:    use the default method (currently method Sandia_LUT)
//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia_LL:  ntri = sum (sum ((L * L) .* L))
//  4:  Sandia_UU:  ntri = sum (sum ((U * U) .* U))
//  5:  Sandia_LUT: ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  Sandia_ULT: ntri = sum (sum ((U * L') .* U)).  Note that U=L'.

// A is a square symmetric matrix, of any type.  Its values are ignored.
// Results are undefined for methods 1 and 2 if self-edges exist in A.  Results
// are undefined for all methods if A is unsymmetric.

// The Sandia_* methods all tend to be faster than the Burkhardt or Cohen
// methods.  For the largest graphs, Sandia_LUT tends to be fastest, except for
// the GAP-urand matrix, where the saxpy-based Sandia_LL method (L*L.*L) is
// fastest.  For many small graphs, the saxpy-based Sandia_LL and Sandia_UU
// methods are often faster that the dot-product-based methods.

// Reference for the Burkhardt method:  Burkhardt, Paul. "Graphing Trillions of
// Triangles." Information Visualization 16, no. 3 (July 2017): 157–66.
// https://doi.org/10.1177/1473871616666393.

// Reference for the Cohen method:  J. Cohen, "Graph twiddling in a mapreduce
// world," Computing in Science & Engineering, vol. 11, no. 4, pp. 29–41, 2009.
// https://doi.org/10.1109/MCSE.2009.120

// Reference for the "Sandia_*" methods: Wolf, Deveci, Berry, Hammond,
// Rajamanickam, "Fast linear algebra- based triangle counting with
// KokkosKernels", IEEE HPEC'17, https://dx.doi.org/10.1109/HPEC.2017.8091043

#define LG_FREE_ALL             \
{                               \
    GrB_free (L) ;              \
    GrB_free (U) ;              \
}

#include "LG_internal.h"

//------------------------------------------------------------------------------
// Added by Zhen Peng 12/5/2022
// get_matrix_trace: return the matrix's trace.
//------------------------------------------------------------------------------
static int64_t get_matrix_trace
    (
        GrB_Matrix *L,
        GrB_Matrix *U,
        GrB_Matrix A,
        char *msg
    )
{
  GrB_Index num_rows;
  GRB_TRY (GrB_Matrix_nrows (&num_rows, A)) ;
  int64_t trace = 0;
  for (int64_t r_i = 0; r_i < num_rows; ++r_i) {
    int64_t val;
    GRB_TRY (GrB_Matrix_extractElement (&val, A, r_i, r_i)) ;
    trace += val;
  }

  return trace;
}

//------------------------------------------------------------------------------
// tricount_prep: construct L and U for LAGr_TriangleCount
//------------------------------------------------------------------------------

//static int tricount_prep
//    (
//        GrB_Matrix *L,      // if present, compute L = tril (A,-1)
//        GrB_Matrix *U,      // if present, compute U = triu (A, 1)
//        GrB_Matrix A,       // input matrix
//        char *msg
//    )
//{
//  GrB_Index n ;
//  GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
//
//  if (L != NULL)
//  {
//    // L = tril (A,-1)
//    GRB_TRY (GrB_Matrix_new (L, GrB_BOOL, n, n)) ;
//    GRB_TRY (GrB_select (*L, NULL, NULL, GrB_TRIL, A, (int64_t) (-1),
//                         NULL)) ;
//    GRB_TRY (GrB_Matrix_wait (*L, GrB_MATERIALIZE)) ;
//  }
//
//  if (U != NULL)
//  {
//    // U = triu (A,1)
//    GRB_TRY (GrB_Matrix_new (U, GrB_BOOL, n, n)) ;
//    GRB_TRY (GrB_select (*U, NULL, NULL, GrB_TRIU, A, (int64_t) 1, NULL)) ;
//    GRB_TRY (GrB_Matrix_wait (*U, GrB_MATERIALIZE)) ;
//  }
//  return (GrB_SUCCESS) ;
//}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    GrB_free (&C) ;                         \
    GrB_free (&L) ;                         \
    GrB_free (&T) ;                         \
    GrB_free (&U) ;                         \
    LAGraph_Free ((void **) &P, NULL) ;     \
}

int LAGr_MaskedSpGEMM
    (
        // output:
        uint64_t *ntriangles,
        // input:
        const LAGraph_Graph G,
        const LAGraph_Graph M, // Mask
        LAGr_MaskedSpGEMM_Method *p_method,
//        LAGr_TriangleCount_Presort *p_presort,
        char *msg
    )
{

  //--------------------------------------------------------------------------
  // check inputs
  //--------------------------------------------------------------------------

  LG_CLEAR_MSG ;
  GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL ;
  int64_t *P = NULL ;

  // get the method
  LAGr_MaskedSpGEMM_Method method ;
  method = (p_method == NULL) ? LAGr_MaskedSpGEMM_AutoMethod : (*p_method) ;
  LG_ASSERT_MSG (
      method == LAGr_MaskedSpGEMM_AutoMethod ||  // 0: use auto method
      method == LAGr_MaskedSpGEMM_NoMask  ||  // 1: sum (sum ((A^2) .* A))/6
      method == LAGr_MaskedSpGEMM_AllZeroMask      ||  // 2: sum (sum ((L * U) .*A))/2
      method == LAGr_MaskedSpGEMM_UseMask,  // 3: sum (sum ((L * L) .* L))
      GrB_INVALID_VALUE, "method is invalid") ;

//  // get the presort
//  LAGr_TriangleCount_Presort presort ;
//  presort = (p_presort == NULL) ? LAGr_TriangleCount_AutoSort : (*p_presort) ;
//  LG_ASSERT_MSG (
//      presort == LAGr_TriangleCount_NoSort     ||
//      presort == LAGr_TriangleCount_Ascending  ||
//      presort == LAGr_TriangleCount_Descending ||
//      presort == LAGr_TriangleCount_AutoSort,
//      GrB_INVALID_VALUE, "presort is invalid") ;

  LG_TRY (LAGraph_CheckGraph (G, msg)) ;
  LG_ASSERT (ntriangles != NULL, GrB_NULL_POINTER) ;
//  LG_ASSERT (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED) ;
//
//  LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
//                  (G->kind == LAGraph_ADJACENCY_DIRECTED &&
//                   G->is_symmetric_structure == LAGraph_TRUE)),
//                 LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
//                 "G->A must be known to be symmetric") ;

  if (method == LAGr_MaskedSpGEMM_AutoMethod)
  {
    // AutoMethod: use default, Sandia_LUT: sum (sum ((L * U') .* L))
    method = LAGr_MaskedSpGEMM_UseMask ;
  }

//  // only the Sandia_* methods can benefit from the presort
//  bool method_can_use_presort =
//      method == LAGr_TriangleCount_Sandia_LL || // sum (sum ((L * L) .* L))
//      method == LAGr_TriangleCount_Sandia_UU || // sum (sum ((U * U) .* U))
//      method == LAGr_TriangleCount_Sandia_LUT || // sum (sum ((L * U') .* L))
//      method == LAGr_TriangleCount_Sandia_ULT ; // sum (sum ((U * L') .* U))
//
  GrB_Matrix A = G->A ;
  GrB_Vector Degree = G->out_degree ;
  GrB_Matrix Mask = M->A;

//  bool auto_sort = (presort == LAGr_TriangleCount_AutoSort) ;
//  if (auto_sort && method_can_use_presort)
//  {
//    LG_ASSERT_MSG (Degree != NULL,
//                   LAGRAPH_NOT_CACHED, "G->out_degree is required") ;
//  }

  //--------------------------------------------------------------------------
  // initializations
  //--------------------------------------------------------------------------

  GrB_Index n ;
  GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
  GRB_TRY (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
#if LAGRAPH_SUITESPARSE
  GrB_Semiring semiring = GxB_PLUS_TIMES_FP64 ;
//  GrB_Semiring semiring = GxB_PLUS_PAIR_INT64 ;
#else
  GrB_Semiring semiring = LAGraph_plus_one_int64 ;
#endif
  GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;

//  //--------------------------------------------------------------------------
//  // heuristic sort rule
//  //--------------------------------------------------------------------------
//
//  if (!method_can_use_presort)
//  {
//    // no sorting for the Burkhardt and Cohen methods: presort parameter
//    // is ignored.
//    presort = LAGr_TriangleCount_NoSort ;
//  }
//  else if (auto_sort)
//  {
//    // auto selection of sorting method for Sandia_* methods
//    presort = LAGr_TriangleCount_NoSort ; // default is not to sort
//
//    if (method_can_use_presort)
//    {
//      // This rule is very similar to Scott Beamer's rule in the GAP TC
//      // benchmark, except that it is extended to handle the ascending
//      // sort needed by methods 3 and 5.  It also uses a stricter rule,
//      // since the performance of triangle counting in SuiteSparse:
//      // GraphBLAS is less sensitive to the sorting as compared to the
//      // GAP algorithm.  This is because the dot products in SuiteSparse:
//      // GraphBLAS use binary search if one vector is very sparse
//      // compared to the other.  As a result, SuiteSparse:GraphBLAS needs
//      // the sort for fewer matrices, as compared to the GAP algorithm.
//
//      // With this rule, the GAP-kron and GAP-twitter matrices are
//      // sorted, and the others remain unsorted.  With the rule in the
//      // GAP tc.cc benchmark, GAP-kron and GAP-twitter are sorted, and so
//      // is GAP-web, but GAP-web is not sorted here.
//
//#define NSAMPLES 1000
//      GrB_Index nvals ;
//      GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
//      if (n > NSAMPLES && ((double) nvals / ((double) n)) >= 10)
//      {
//        // estimate the mean and median degrees
//        double mean, median ;
//        LG_TRY (LAGr_SampleDegree (&mean, &median,
//                                   G, true, NSAMPLES, n, msg)) ;
//        // sort if the average degree is very high vs the median
//        if (mean > 4 * median)
//        {
//          switch (method)
//          {
//            case LAGr_TriangleCount_Sandia_LL:
//              // 3:sum (sum ((L * L) .* L))
//              presort = LAGr_TriangleCount_Ascending  ;
//              break ;
//            case LAGr_TriangleCount_Sandia_UU:
//              // 4: sum (sum ((U * U) .* U))
//              presort = LAGr_TriangleCount_Descending ;
//              break ;
//            default:
//            case LAGr_TriangleCount_Sandia_LUT:
//              // 5: sum (sum ((L * U') .* L))
//              presort = LAGr_TriangleCount_Ascending  ;
//              break ;
//            case LAGr_TriangleCount_Sandia_ULT:
//              // 6: sum (sum ((U * L') .* U))
//              presort = LAGr_TriangleCount_Descending ;
//              break ;
//          }
//        }
//      }
//    }
//  }

  //--------------------------------------------------------------------------
  // sort the input matrix, if requested
  //--------------------------------------------------------------------------

//  if (presort != LAGr_TriangleCount_NoSort)
//  {
//    // P = permutation that sorts the rows by their degree
//    LG_TRY (LAGr_SortByDegree (&P, G, true,
//                               presort == LAGr_TriangleCount_Ascending, msg)) ;
//
//    // T = A (P,P) and typecast to boolean
//    GRB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
//    GRB_TRY (GrB_extract (T, NULL, NULL, A, (GrB_Index *) P, n,
//                          (GrB_Index *) P, n, NULL)) ;
//    A = T ;
//
//    // free workspace
//    LG_TRY (LAGraph_Free ((void **) &P, NULL)) ;
//  }

  //--------------------------------------------------------------------------
  // count triangles
  //--------------------------------------------------------------------------

  int64_t ntri = -1;

  /// Added by Zhen Peng on 1/9/2023
  /// Turning off the sorting
  /**
   *  GrB_Info GrB_Descriptor_set // set a parameter in a descriptor
      (
        GrB_Descriptor desc, // descriptor to modify:: the descriptor object to be used during the op, i.e., gemm.
        GrB_Desc_Field field, // parameter to change:: GxB_SORT
        GrB_Desc_Value val // value to change it to:: 0
      ) ;
   */
  GrB_Descriptor_set(GrB_DESC_S, GxB_SORT, 0);

  switch (method)
  {
//    case LAGr_MaskedSpGEMM_NoMask:  // 7: sum (sum ((A^2) .* A)) / 6, no masking, elementwise multiplication
//
//    GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_TIMES_INT64, C, A, NULL)) ;
////            GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_BAND_INT64, C, A, NULL)) ; // Wrong answers! Not sure why.
////            ntri = get_matrix_trace(&L, &U, C, msg);
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      ntri /= 6 ;
//      break ;

    case LAGr_MaskedSpGEMM_NoMask:  // 1: Use no mask

      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, GrB_DESC_S)) ;
      GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_TIMES_FP64, C, Mask, NULL)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      ntri /= 6 ;
      break ;

    case LAGr_MaskedSpGEMM_AllZeroMask: // 2: Read a mask from files

//    LG_TRY (tricount_prep (&L, &U, A, msg)) ;
//      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, GrB_DESC_S)) ;
      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      ntri /= 2 ;
      break ;

    case LAGr_MaskedSpGEMM_UseMask: // 3: Currently uses the input matrix A as the mask

      // using the masked saxpy3 method
//      LG_TRY (tricount_prep (&L, NULL, A, msg)) ;
      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
      break ;

//    case LAGr_TriangleCount_Sandia_UU: // 4: sum (sum ((U * U) .* U))
//
//      // using the masked saxpy3 method
//    LG_TRY (tricount_prep (NULL, &U, A, msg)) ;
//      GRB_TRY (GrB_mxm (C, U, NULL, semiring, U, U, GrB_DESC_S)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      break ;
//
//    default:
//    case LAGr_TriangleCount_Sandia_LUT: // 5: sum (sum ((L * U') .* L))
//
//      // This tends to be the fastest method for most large matrices, but
//      // the Sandia_ULT method is also very fast.
//
//      // using the masked dot product
//    LG_TRY (tricount_prep (&L, &U, A, msg)) ;
//      GRB_TRY (GrB_mxm (C, L, NULL, semiring, L, U, GrB_DESC_ST1)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      break ;
//
//    case LAGr_TriangleCount_Sandia_ULT: // 6: sum (sum ((U * L') .* U))
//
//      // using the masked dot product
//    LG_TRY (tricount_prep (&L, &U, A, msg)) ;
//      GRB_TRY (GrB_mxm (C, U, NULL, semiring, U, L, GrB_DESC_ST1)) ;
//      GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
//      break ;
  }

  //--------------------------------------------------------------------------
  // return result
  //--------------------------------------------------------------------------

  LG_FREE_ALL ;
  if (p_method != NULL) (*p_method) = method ;
//  if (p_presort != NULL) (*p_presort) = presort ;
  (*ntriangles) = (uint64_t) ntri ;
  return (GrB_SUCCESS) ;
}


int LAGr_MaskedSpGEMM_print_matrix
    (
//        // output:
//        GrB_Matrix *C_output,
        // input:
        const LAGraph_Graph G,
        const LAGraph_Graph M, // Mask
        LAGr_MaskedSpGEMM_Method *p_method,
//    LAGr_TriangleCount_Presort *presort,
        char *msg
)
{

  //--------------------------------------------------------------------------
  // check inputs
  //--------------------------------------------------------------------------

  LG_CLEAR_MSG ;
//  GrB_Matrix L = NULL, U = NULL, T = NULL ;
  GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL ;
  int64_t *P = NULL ;

  // get the method
  LAGr_MaskedSpGEMM_Method method ;
  method = (p_method == NULL) ? LAGr_MaskedSpGEMM_AutoMethod : (*p_method) ;
  LG_ASSERT_MSG (
      method == LAGr_MaskedSpGEMM_AutoMethod ||  // 0: use auto method
      method == LAGr_MaskedSpGEMM_NoMask  ||  // 1: sum (sum ((A^2) .* A))/6
      method == LAGr_MaskedSpGEMM_AllZeroMask      ||  // 2: sum (sum ((L * U) .*A))/2
      method == LAGr_MaskedSpGEMM_UseMask,  // 3: sum (sum ((L * L) .* L))
      GrB_INVALID_VALUE, "method is invalid") ;

  LG_TRY (LAGraph_CheckGraph (G, msg)) ;
//  LG_ASSERT (ntriangles != NULL, GrB_NULL_POINTER) ;
//  LG_ASSERT (C_output != NULL, GrB_NULL_POINTER) ;
//  LG_ASSERT (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED) ;
//

  if (method == LAGr_MaskedSpGEMM_AutoMethod)
  {
    // AutoMethod: use default, Sandia_LUT: sum (sum ((L * U') .* L))
    method = LAGr_MaskedSpGEMM_NoMask ;
  }

//
  GrB_Matrix A = G->A ;
  GrB_Vector Degree = G->out_degree ;
  GrB_Matrix Mask = M->A;

  //--------------------------------------------------------------------------
  // initializations
  //--------------------------------------------------------------------------

  GrB_Index n ;
  GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
  GRB_TRY (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
//  GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
//  GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
#if LAGRAPH_SUITESPARSE
  GrB_Semiring semiring = GxB_PLUS_TIMES_FP64 ;
//  GrB_Semiring semiring = GxB_PLUS_PAIR_FP64 ;
//  GrB_Semiring semiring = GxB_PLUS_PAIR_INT64 ;
#else
  GrB_Semiring semiring = LAGraph_plus_one_int64 ;
#endif
  GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;

  //--------------------------------------------------------------------------
  // count triangles
  //--------------------------------------------------------------------------

  int64_t ntri = -1;

  /// Added by Zhen Peng on 1/9/2023
  /// Turning off the sorting
  /**
   *  GrB_Info GrB_Descriptor_set // set a parameter in a descriptor
      (
        GrB_Descriptor desc, // descriptor to modify:: the descriptor object to be used during the op, i.e., gemm.
        GrB_Desc_Field field, // parameter to change:: GxB_SORT
        GrB_Desc_Value val // value to change it to:: 0
      ) ;
   */
  GrB_Descriptor_set(GrB_DESC_S, GxB_SORT, 0);

  switch (method)
  {

    case LAGr_MaskedSpGEMM_NoMask:  // 1: Use no mask

      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, GrB_DESC_S)) ;
      GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_TIMES_FP64, C, Mask, NULL)) ;

//      printf("#### Output C ####\n");
//      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
      {
        char filename[] = "output.LAGr_MaskedSpGEMM_NoMask.txt";
        FILE *fout = fopen(filename, "w");
        if (!fout) {
          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
          exit(EXIT_FAILURE);
        }
        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
        fclose(fout);
      }
      break ;

    case LAGr_MaskedSpGEMM_AllZeroMask: // 2: Read a mask from files

      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;

//      printf("#### Output C ####\n");
//      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
      {
        char filename[] = "output.LAGr_MaskedSpGEMM_AllZeroMask.txt";
        FILE *fout = fopen(filename, "w");
        if (!fout) {
          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
          exit(EXIT_FAILURE);
        }
        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
        fclose(fout);
      }
      break ;

    case LAGr_MaskedSpGEMM_UseMask: // 3: Currently uses the input matrix A as the mask

      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;

//      printf("#### Output C ####\n");
//      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_SHORT, stdout, msg));
      {
        char filename[] = "output.LAGr_MaskedSpGEMM_UseMask.txt";
        FILE *fout = fopen(filename, "w");
        if (!fout) {
          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
          exit(EXIT_FAILURE);
        }
        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
        fclose(fout);
      }
      break ;

  }

//  printf("#### Output C ####\n");
//  LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
//  LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE, stdout, msg));
  //--------------------------------------------------------------------------
  // return result
  //--------------------------------------------------------------------------

  LG_FREE_ALL ;
  if (p_method != NULL) (*p_method) = method ;
//  if (p_presort != NULL) (*p_presort) = presort ;
//  (*ntriangles) = (uint64_t) ntri ;
  return (GrB_SUCCESS) ;
}
