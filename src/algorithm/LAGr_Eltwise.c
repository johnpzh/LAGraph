//------------------------------------------------------------------------------
// LAGr_Eltwise: Elementwise Multiplication
//------------------------------------------------------------------------------

#define LG_FREE_ALL             \
{                               \
    GrB_free (L) ;              \
    GrB_free (U) ;              \
}

#include "LG_internal.h"


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


int LAGr_Eltwise
(
        // output:
        uint64_t *ntriangles,
        // input:
        const LAGraph_Graph G1,
        const LAGraph_Graph G2,
//        const LAGraph_Graph M, // Mask
        LAGr_Eltwise_Method *p_method,
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
  LAGr_Eltwise_Method method ;
  method = (p_method == NULL) ? LAGr_Eltwise_NoMask : (*p_method) ;
  LG_ASSERT_MSG (
      method == LAGr_Eltwise_NoMask,  // 0: use no mask
      GrB_INVALID_VALUE, "method is invalid") ;

  LG_TRY (LAGraph_CheckGraph (G1, msg)) ;
  LG_TRY (LAGraph_CheckGraph (G2, msg)) ;
  LG_ASSERT (ntriangles != NULL, GrB_NULL_POINTER) ;

  GrB_Matrix A = G1->A ;
  GrB_Vector Degree_A = G1->out_degree ;
  GrB_Matrix B = G2->A ;
  GrB_Vector Degree_B = G2->out_degree ;
//  GrB_Matrix Mask = M->A;

  //--------------------------------------------------------------------------
  // initializations
  //--------------------------------------------------------------------------

  GrB_Index n ;
  GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
  GRB_TRY (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
//#if LAGRAPH_SUITESPARSE
  GrB_Semiring semiring = GxB_PLUS_TIMES_FP64 ;
//  GrB_Semiring semiring = GxB_PLUS_PAIR_INT64 ;
//#else
//  GrB_Semiring semiring = LAGraph_plus_one_int64 ;
//#endif
//  GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;

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
  /// Changed by Zhen Peng on 4/12/2023
  /// Turn on the sorting because COMET uses sorting in SpGEMM.
//  GrB_Descriptor_set(GrB_DESC_S, GxB_SORT, 1);
//  GrB_Descriptor_set(GrB_DESC_S, GxB_SORT, 0);
//  GxB_Desc_set(GrB_DESC_S, GxB_SORT, 0);
//  GxB_Desc_set(GrB_DESC_S, GxB_SORT, 1);
//  GxB_set(GrB_DESC_S, GxB_SORT, 1);
  GrB_Descriptor DESC_FOR_SORT;
  GrB_Descriptor_new(&DESC_FOR_SORT);
  GxB_set(DESC_FOR_SORT, GxB_SORT, 17);  // sort the multiplication
//  GxB_set(DESC_FOR_SORT, GxB_SORT, 0);  // lazy sort


  switch (method)
  {
    case LAGr_Eltwise_NoMask:  // 0: Use no mask

//      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, DESC_FOR_SORT)) ;
      GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_TIMES_FP64, A, B, NULL)) ;

//      printf("#### Output C ####\n");
//      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
//      {
//        char filename[] = "output.LAGr_MaskedSpGEMM_NoMask.txt";
//        printf("[ Writing to %s ... ]\n", filename);
//        FILE *fout = fopen(filename, "w");
//        if (!fout) {
//          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
//          exit(EXIT_FAILURE);
//        }
//        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
//        fclose(fout);
//      }

      break ;

    case LAGr_Eltwise_UseMask:
      break ;
//    case LAGr_MaskedSpGEMM_UseMask: // 3: Currently uses the input matrix A as the mask
//
//      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, DESC_FOR_SORT)) ;
//      break ;

  }

  //--------------------------------------------------------------------------
  // return result
  //--------------------------------------------------------------------------

//  LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));

  LG_FREE_ALL ;
  if (p_method != NULL) (*p_method) = method ;
//  if (p_presort != NULL) (*p_presort) = presort ;
  (*ntriangles) = (uint64_t) ntri ;
  /// Added by Zhen Peng on 4/27/2023
  GrB_Descriptor_free(&DESC_FOR_SORT);
  return (GrB_SUCCESS) ;
}

//
//int LAGr_MaskedSpGEMM_print_matrix
//    (
////        // output:
////        GrB_Matrix *C_output,
//        // input:
//        const LAGraph_Graph G,
//        const LAGraph_Graph M, // Mask
//        LAGr_MaskedSpGEMM_Method *p_method,
////    LAGr_TriangleCount_Presort *presort,
//        char *msg
//)
//{
//
//  //--------------------------------------------------------------------------
//  // check inputs
//  //--------------------------------------------------------------------------
//
//  LG_CLEAR_MSG ;
////  GrB_Matrix L = NULL, U = NULL, T = NULL ;
//  GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL ;
//  int64_t *P = NULL ;
//
//  // get the method
//  LAGr_MaskedSpGEMM_Method method ;
//  method = (p_method == NULL) ? LAGr_MaskedSpGEMM_AutoMethod : (*p_method) ;
//  LG_ASSERT_MSG (
//      method == LAGr_MaskedSpGEMM_AutoMethod ||  // 0: use auto method
//      method == LAGr_MaskedSpGEMM_NoMask  ||  // 1: sum (sum ((A^2) .* A))/6
//      method == LAGr_MaskedSpGEMM_AllZeroMask      ||  // 2: sum (sum ((L * U) .*A))/2
//      method == LAGr_MaskedSpGEMM_UseMask,  // 3: sum (sum ((L * L) .* L))
//      GrB_INVALID_VALUE, "method is invalid") ;
//
//  LG_TRY (LAGraph_CheckGraph (G, msg)) ;
////  LG_ASSERT (ntriangles != NULL, GrB_NULL_POINTER) ;
////  LG_ASSERT (C_output != NULL, GrB_NULL_POINTER) ;
////  LG_ASSERT (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED) ;
////
//
//  if (method == LAGr_MaskedSpGEMM_AutoMethod)
//  {
//    // AutoMethod: use default, Sandia_LUT: sum (sum ((L * U') .* L))
//    method = LAGr_MaskedSpGEMM_NoMask ;
//  }
//
////
//  GrB_Matrix A = G->A ;
//  GrB_Vector Degree = G->out_degree ;
//  GrB_Matrix Mask = M->A;
//
//  //--------------------------------------------------------------------------
//  // initializations
//  //--------------------------------------------------------------------------
//
//  GrB_Index n ;
//  GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
//  GRB_TRY (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
////  GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
////  GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
//#if LAGRAPH_SUITESPARSE
//  GrB_Semiring semiring = GxB_PLUS_TIMES_FP64 ;
////  GrB_Semiring semiring = GxB_PLUS_PAIR_FP64 ;
////  GrB_Semiring semiring = GxB_PLUS_PAIR_INT64 ;
//#else
//  GrB_Semiring semiring = LAGraph_plus_one_int64 ;
//#endif
//  GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;
//
//  //--------------------------------------------------------------------------
//  // count triangles
//  //--------------------------------------------------------------------------
//
//  int64_t ntri = -1;
//
//  /// Added by Zhen Peng on 1/9/2023
//  /// Turning off the sorting
//  /**
//   *  GrB_Info GrB_Descriptor_set // set a parameter in a descriptor
//      (
//        GrB_Descriptor desc, // descriptor to modify:: the descriptor object to be used during the op, i.e., gemm.
//        GrB_Desc_Field field, // parameter to change:: GxB_SORT
//        GrB_Desc_Value val // value to change it to:: 0
//      ) ;
//   */
////  GrB_Descriptor_set(GrB_DESC_S, GxB_SORT, 0);
////  GxB_set(GrB_DESC_S, GxB_SORT, 27);
//  GxB_Desc_set(GrB_DESC_S, GxB_SORT, 27);
//  {// test
//    int sort_set = -1;
////    GxB_get(GrB_DESC_S, GxB_SORT, &sort_set);
//    GxB_Desc_get(GrB_DESC_S, GxB_SORT, &sort_set);
//    fprintf(stderr, "[! before GrB_DESC_S GxB_SORT: %d ]\n", sort_set);
//  }
////  {//test
////    int other_set = 17;
////    GxB_get(GrB_DESC_S, GrB_STRUCTURE, &other_set);
////    fprintf(stderr, " [! before GrB_DESC_S GrB_STRUCTURE: %d ]\n", other_set);
////    GxB_set(GrB_DESC_S, GrB_STRUCTURE, 27);
////    GxB_get(GrB_DESC_S, GrB_STRUCTURE, &other_set);
////    fprintf(stderr, " [! chaged GrB_DESC_S GrB_STRUCTURE: %d ]\n", other_set);
////  }
////  GrB_Descriptor DESC_FOR_SORT = GrB_DESC_S;
////  GrB_Descriptor DESC_FOR_SORT = (GrB_Descriptor) malloc(sizeof(struct GB_Descriptor_opaque));
////  GrB_Descriptor DESC_FOR_SORT = GrB_Descriptor_new(GrB_DESC_S);
////  GrB_Descriptor DESC_FOR_SORT = GrB_Descriptor_new(GrB_DESC_T1);
//  GrB_Descriptor DESC_FOR_SORT;
//  GrB_Descriptor_new(&DESC_FOR_SORT);
////  *DESC_FOR_SORT = *GrB_DESC_S;
////  DESC_FOR_SORT->magic = 1;
////  GxB_Desc_set(DESC_FOR_SORT, GxB_SORT, 17);
//  GxB_set(DESC_FOR_SORT, GxB_SORT, 17);
////  GxB_set(DESC_FOR_SORT, GxB_SORT, 0);
//  {// test
//    int sort_set = -1;
//    GxB_get(DESC_FOR_SORT, GxB_SORT, &sort_set);
//    fprintf(stderr, "[! before DESC_FOR_SORT GxB_SORT: %d ]\n", sort_set);
//  }
//
//  switch (method)
//  {
//
//    case LAGr_MaskedSpGEMM_NoMask:  // 1: Use no mask
//
////      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, A, DESC_FOR_SORT)) ;
////      GRB_TRY (GrB_eWiseMult(C, NULL, NULL, GrB_TIMES_FP64, C, Mask, NULL)) ;
//
////      printf("#### Output C ####\n");
////      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
//      {
//        char filename[] = "output.LAGr_MaskedSpGEMM_NoMask.txt";
//        printf("[ Writing to %s ... ]\n", filename);
//        FILE *fout = fopen(filename, "w");
//        if (!fout) {
//          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
//          exit(EXIT_FAILURE);
//        }
//        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
//        fclose(fout);
//      }
//      break ;
//
//    case LAGr_MaskedSpGEMM_AllZeroMask: // 2: Read a mask from files
//
////      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, DESC_FOR_SORT)) ;
//
////      printf("#### Output C ####\n");
////      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
//      {
//        char filename[] = "output.LAGr_MaskedSpGEMM_AllZeroMask.txt";
//        printf("[ Writing to %s ... ]\n", filename);
//        FILE *fout = fopen(filename, "w");
//        if (!fout) {
//          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
//          exit(EXIT_FAILURE);
//        }
//        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
//        fclose(fout);
//      }
//      break ;
//
//    case LAGr_MaskedSpGEMM_UseMask: // 3: Currently uses the input matrix A as the mask
//
////      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, GrB_DESC_S)) ;
//      GRB_TRY (GrB_mxm (C, Mask, NULL, semiring, A, A, DESC_FOR_SORT)) ;
//
////      printf("#### Output C ####\n");
////      LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_SHORT, stdout, msg));
//      {
//        char filename[] = "output.LAGr_MaskedSpGEMM_UseMask.txt";
//        printf("[ Writing to %s ... ]\n", filename);
//        FILE *fout = fopen(filename, "w");
//        if (!fout) {
//          fprintf(stderr, "Error %s:%d: cannot create file %s .\n", __FILE__, __LINE__, filename);
//          exit(EXIT_FAILURE);
//        }
//        LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, fout, msg));
//        fclose(fout);
//      }
//      break ;
//
//  }
//  {// test
//    int sort_set = -1;
//    GxB_get(GrB_DESC_S, GxB_SORT, &sort_set);
//    fprintf(stderr, "[! after GrB_DESC_S GxB_SORT: %d ]\n", sort_set);
//  }
//  {// test
//    int sort_set = -1;
//    GxB_get(DESC_FOR_SORT, GxB_SORT, &sort_set);
//    fprintf(stderr, "[! after DESC_FOR_SORT GxB_SORT: %d ]\n", sort_set);
//  }
//
//
////  printf("#### Output C ####\n");
////  LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE_VERBOSE, stdout, msg));
////  LAGRAPH_TRY (LAGraph_Matrix_Print(C, LAGraph_COMPLETE, stdout, msg));
//  //--------------------------------------------------------------------------
//  // return result
//  //--------------------------------------------------------------------------
//
//  LG_FREE_ALL ;
//
//  /// Add by Zhen Peng on 4/27/2023
////  free(DESC_FOR_SORT);
//  GrB_Descriptor_free(&DESC_FOR_SORT);
//  /// End add
//  if (p_method != NULL) (*p_method) = method ;
////  if (p_presort != NULL) (*p_presort) = presort ;
////  (*ntriangles) = (uint64_t) ntri ;
//  return (GrB_SUCCESS) ;
//}
