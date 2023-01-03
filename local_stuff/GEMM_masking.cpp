// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_mxm.c
    //--------------------------------------------------------------------------
    // T = A*B, A'*B, A*B', or A'*B', also using the mask if present
    //--------------------------------------------------------------------------
    GB_OK (GB_AxB_meta (T, C, C_replace, C->is_csc, MT, &M_transposed, M,
        Mask_comp, Mask_struct, accum, A, B, semiring, A_transpose,
        B_transpose, flipxy, &mask_applied, &done_in_place, AxB_method,
        do_sort, Context)) ;

    //--------------------------------------------------------------------------
    // C<M> = accum (C,T): accumulate the results into C via the mask
    //--------------------------------------------------------------------------
        // C<M> = accum (C,T)
        // GB_accum_mask also conforms C to its desired hypersparsity.
        info = GB_accum_mask (C, M, (M_transposed) ? MT : NULL, accum, &T,
            C_replace, Mask_comp, Mask_struct, Context) ;


// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_meta.c
        GB_AxB_meta_adotb_control (&tentative_axb_method, C_in, M_in,
            Mask_comp, B_in, A_in, accum, semiring_in, flipxy, can_do_in_place,
            allow_scale, A_in_is_diagonal, AxB_method, Context) ;
        GB_AxB_meta_adotb_control (&tentative_axb_method, C_in, M_in,
            Mask_comp, A_in, B_in, accum, semiring_in, flipxy, can_do_in_place,
            allow_scale, B_in_is_diagonal, AxB_method, Context) ;

        GB_OK (GB_transpose_cast (MT, GrB_BOOL, C_is_csc, M_in, Mask_struct,
            Context)) ;

        //----------------------------------------------------------------------
        // select the method for C<M>=A'*B
        //----------------------------------------------------------------------

        GB_AxB_meta_adotb_control (&axb_method, C_in, M,
            Mask_comp, A, B, accum, semiring, flipxy, can_do_in_place,
            allow_scale, B_is_diagonal, AxB_method, Context) ;

                // C<M>=A'*B via dot, or C_in<M>+=A'*B if in-place
                GB_OK (GB_AxB_dot (C, (can_do_in_place) ? C_in : NULL,
                    M, Mask_comp, Mask_struct, accum, A, B, semiring, flipxy,
                    mask_applied, done_in_place, Context)) ;

                // C = A'*B via saxpy: Gustavson + Hash method
                GBURBLE ("C%s=A'*B, saxpy (transposed %s) ", M_str, A_str) ;
                GB_OK (GB_AxB_saxpy (C, can_do_in_place ? C_in : NULL, M,
                    Mask_comp, Mask_struct, accum, AT, B, semiring, flipxy,
                    mask_applied, done_in_place, AxB_method, do_sort,
                    Context)) ;

                // C<M>=A*B' via dot product, or C_in<M>+=A*B' if in-place
                GB_OK (GB_AxB_dot (C, (can_do_in_place) ? C_in : NULL,
                    M, Mask_comp, Mask_struct, accum, AT, BT, semiring, flipxy,
                    mask_applied, done_in_place, Context)) ;

                // C = A*B' via saxpy: Gustavson + Hash method
                GB_OK (GB_AxB_saxpy (C, can_do_in_place ? C_in : NULL, M,
                    Mask_comp, Mask_struct, accum, A, BT, semiring, flipxy,
                    mask_applied, done_in_place, AxB_method, do_sort,
                    Context)) ;

                GB_AxB_saxpy_sparsity (&ignore, &saxpy_method, M, Mask_comp,
                    A, B, Context) ;

        switch (axb_method)
        {
            case GB_USE_COLSCALE : 
                // C = A*D, column scale
                GBURBLE ("C%s=A*B, colscale ", M_str) ;
                GB_OK (GB_AxB_colscale (C, A, B, semiring, flipxy, Context)) ;
                break ;

            case GB_USE_ROWSCALE : 
                // C = D*B, row scale
                GBURBLE ("C%s=A*B, rowscale ", M_str) ;
                GB_OK (GB_AxB_rowscale (C, A, B, semiring, flipxy, Context)) ;
                break ;

            case GB_USE_DOT : 
                // C<M>=A*B via dot product, or C_in<M>+=A*B if in-place.
                GBURBLE ("C%s=A*B', dot_product (transposed %s) ",
                    M_str, A_str) ;
                GB_CLEAR_STATIC_HEADER (AT, &AT_header) ;
                GB_OK (GB_transpose_cast (AT, atype_cast, true, A, A_is_pattern,
                    Context)) ;
                GB_OK (GB_AxB_dot (C, (can_do_in_place) ? C_in : NULL,
                    M, Mask_comp, Mask_struct, accum, AT, B, semiring, flipxy,
                    mask_applied, done_in_place, Context)) ;
                break ;

            default : 
                // C = A*B via saxpy: Gustavson + Hash method
                GBURBLE ("C%s=A*B, saxpy ", M_str) ;
                GB_OK (GB_AxB_saxpy (C, can_do_in_place ? C_in : NULL, M,
                    Mask_comp, Mask_struct, accum, A, B, semiring, flipxy,
                    mask_applied, done_in_place, AxB_method, do_sort,
                    Context)) ;
                break ;
        }


// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_meta_adotb_control.c
//------------------------------------------------------------------------------
// GB_AxB_meta_adotb_control: determine method for computing C=A'*B
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------
void GB_AxB_meta_adotb_control
(
    // output:
    int *axb_method,
    // input:
    const GrB_Matrix C_in,
    const GrB_Matrix M,
    bool Mask_comp,
    const GrB_Matrix A,
    const GrB_Matrix B,
    const GrB_BinaryOp accum,
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    bool flipxy,
    bool can_do_in_place,
    bool allow_scale,
    bool B_is_diagonal,
    GrB_Desc_Value AxB_method,
    GB_Context Context
)
{

    // C=A'*B is being computed: use the dot product without computing A'
    // or use the saxpy (Gustavson) method

    // use saxpy by default, unless selecting other methods below
    (*axb_method) = GB_USE_SAXPY ;

    // If the mask is present, only entries for which M(i,j)=1 are
    // computed, which makes this method very efficient when the mask is
    // very sparse (triangle counting, for example).  Each entry C(i,j) for
    // which M(i,j)=1 is computed via a dot product, C(i,j) =
    // A(:,i)'*B(:,j).  If the mask is not present, the dot-product method
    // is very slow in general, and thus the saxpy method is usually used
    // instead.

    if (allow_scale && M == NULL
        && !GB_IS_BITMAP (A)     // TODO: A'*D colscale with A bitmap
        && B_is_diagonal)
    { 
        // C = A'*D, col scale
        (*axb_method) = GB_USE_COLSCALE ;
    }
    else if (allow_scale && M == NULL
        && !GB_IS_BITMAP (B)     // TODO: D*B rowscale with B bitmap
        && GB_is_diagonal (A, Context))
    { 
        // C = D*B, row scale
        (*axb_method) = GB_USE_ROWSCALE ;
    }
    else if (AxB_method == GxB_DEFAULT)
    {
        // auto selection for A'*B
        bool C_out_iso = false ;    // ignored unless C can be done in-place
        if (can_do_in_place && C_in != NULL)
        { 
            // check if C will be iso on output (for dot4 control only).
            // Ignored if dot4 C_in is not present or C cannot be
            // computed in-place.
            C_out_iso = GB_iso_AxB (NULL, A, B, A->vlen, semiring, flipxy,
                false) ;
        }
        if (GB_AxB_dot4_control (C_out_iso, can_do_in_place ? C_in : NULL,
            M, Mask_comp, accum, semiring))
        { 
            // C+=A'*B can be done with dot4
            (*axb_method) = GB_USE_DOT ;
        }
        else if (GB_AxB_dot3_control (M, Mask_comp))
        { 
            // C<M>=A'*B uses the masked dot product method (dot3)
            (*axb_method) = GB_USE_DOT ;
        }
        else if (GB_AxB_dot2_control (A, B, Context))
        { 
            // C=A'*B or C<!M>=A'B* can efficiently use the dot2 method
            (*axb_method) = GB_USE_DOT ;
        }
    }
    else if (AxB_method == GxB_AxB_DOT)
    { 
        // user selection for A'*B
        (*axb_method) = GB_USE_DOT ;
    }
}


// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_dot.c
//------------------------------------------------------------------------------
// GB_AxB_dot: C<M>=A'*B using dot products
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// Parallel matrix-matrix multiply, A'*B, with optional mask M.  This
// method is used by GrB_mxm, GrB_vxm, and GrB_mxv.  For both of the latter two
// methods, B on input will be an nrows-by-1 column vxector.

// This function, and the matrices C, M, A, and B are all CSR/CSC agnostic.
// For this discussion, suppose they are CSC, with vlen = # of rows, and vdim =
// # of columns.

// C=A'*B, C<M>=A'*B or C<!M>=A'*B is being computed.  A has not been
// transposed yet (and will not be).  A and B must have the same vector length,
// vlen (as if both A and B are CSC matrices with the same number of rows, for
// example).  GB_AxB_dot2 and GB_AxB_dot3 operate on A' without forming it.
// GB_AxB_dot2 computes C=A'*B and C<!M>=A'*B, and it takes Omega(m*n) time,
// if C is m-by-n.  It is thus only suitable for cases when A and B are large,
// and C is small.  GB_AxB_dot3 computes C<M>=A'*B, and it only needs to
// examine entries in M, taking Omega(nnz(M)) time.  It can thus be used for
// very large matrices C.  GB_AxB_dot4 computes C+=A'*B when C is dense.

// The output matrix C has not been allocated.  It is an uninitialzed static
// header on input.  The mask M is optional.

// If the result is computed in-place, then the C parameter is ignored, and the
// result is computed in C_in instead.  This case requires the accum operator
// to match the monoid of the semiring.

// The semiring defines C=A*B.  flipxy modifies how the semiring multiply
// operator is applied.  If false, then fmult(aik,bkj) is computed.  If true,
// then the operands are swapped, and fmult(bkj,aij) is done instead.

// Context: the GB_Context containing the # of threads to use, a string of the
// user-callable function that is calling this function (GrB_mxm, GrB_mxv, or
// GxB_vxm) and detailed error reports.


    //--------------------------------------------------------------------------
    // in-place C+=A'*B.  mask is not present (and not applied)
    //--------------------------------------------------------------------------

    if (GB_AxB_dot4_control (C_iso, C_in, M, Mask_comp, accum, semiring))
    { 
        // C_in must be as-if-full on input.  M must be NULL and not
        // complemented.  the C iso case is not handled (where C is iso on
        // output), but C_in might be iso on input.

        #ifdef GB_DEBUGIFY_DEFN
        GB_debugify_mxm (C_iso, GxB_FULL, ztype, M, Mask_struct,
            Mask_comp, semiring, flipxy, A, B) ;
        #endif

        (*mask_applied) = false ;    // no mask to apply
        info = GB_AxB_dot4 (C_in, A, B, semiring, flipxy, done_in_place,
            Context) ;
        if (info != GrB_NO_VALUE)
        { 
            // return if dot4 has handled this case, otherwise fall through
            // to dot2 or dot3 below.
            return (info) ;
        }
    }

    //--------------------------------------------------------------------------
    // C<M>=A'*B: general case
    //--------------------------------------------------------------------------

    if (GB_AxB_dot3_control (M, Mask_comp))
    { 

        // use dot3 if M is present and not complemented, and either sparse or
        // hypersparse
        GBURBLE ("(%sdot3) ", iso_kind) ;
        (*mask_applied) = true ;    // mask is always applied
        (*done_in_place) = false ;
        GrB_Info info ;

        // construct the hyper hashes for A and B
        GB_OK (GB_hyper_hash_build (A, Context)) ;
        GB_OK (GB_hyper_hash_build (B, Context)) ;

        GBURBLE ("(%s%s%s%s = %s'*%s) ",
            GB_sparsity_char_matrix (M),    // C has the same sparsity as M
            Mask_struct ? "{" : "<",
            GB_sparsity_char_matrix (M),
            Mask_struct ? "}" : ">",
            GB_sparsity_char_matrix (A),
            GB_sparsity_char_matrix (B)) ;

        #ifdef GB_DEBUGIFY_DEFN
        GB_debugify_mxm (C_iso, GB_sparsity (M), ztype, M,
            Mask_struct, Mask_comp, semiring, flipxy, A, B) ;
        #endif

        #if defined ( GBCUDA )
        if (!C_iso &&   // fixme for CUDA, remove and create C iso on output
            GB_AxB_dot3_cuda_branch (M, Mask_struct, A, B, semiring,
            flipxy, Context))
        {
            info = (GB_AxB_dot3_cuda (C, M, Mask_struct, A, B, semiring,
                flipxy, Context)) ;
        }
        else
        #endif
        { 
            // use the CPU
            info = (GB_AxB_dot3 (C, C_iso, cscalar, M, Mask_struct, A, B,
                semiring, flipxy, Context)) ;
        }
        return (info) ;
    }

    //--------------------------------------------------------------------------
    // general case: C<M>=A'*B, C<!M>=A'*B, or C=A'*B, not in-place
    //--------------------------------------------------------------------------

    GBURBLE ("(%sdot2) ", iso_kind) ;
    (*mask_applied) = (M != NULL) ; // mask applied if present
    (*done_in_place) = false ;      // TODO: allow dot2 to work in-place
    return (GB_AxB_dot2 (C, C_iso, cscalar, M, Mask_comp, Mask_struct,
        false, A, B, semiring, flipxy, Context)) ;


// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_dot3.c

GB_PUBLIC
GrB_Info GB_AxB_dot3                // C<M> = A'*B using dot product method
(
    GrB_Matrix C,                   // output matrix, static header
    const bool C_iso,               // true if C is iso
    const GB_void *cscalar,         // iso value of C
    const GrB_Matrix M,             // mask matrix
    const bool Mask_struct,         // if true, use the only structure of M
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    GB_Context Context
)
{

    //--------------------------------------------------------------------------
    // phase1: estimate the work to compute each entry in C
    //--------------------------------------------------------------------------

    // The work to compute C(i,j) is held in Cwork [p], if C(i,j) appears in
    // as the pth entry in C.

    #define GB_DOT3
    #define GB_DOT3_PHASE1

    if (M_is_sparse && Mask_struct)
    {
        // special case: M is sparse and structural
        #define GB_MASK_SPARSE_AND_STRUCTURAL
        #include "GB_meta16_factory.c"
        #undef GB_MASK_SPARSE_AND_STRUCTURAL
        // TODO: skip phase1 if A and B are both bitmap/full.
    }
    else
    {
        // general case: M sparse/hyper, structural/valued
        #include "GB_meta16_factory.c"
    }

    #undef GB_DOT3
    #undef GB_DOT3_PHASE1

    //--------------------------------------------------------------------------
    // free the current tasks and construct the tasks for the second phase
    //--------------------------------------------------------------------------

    GB_FREE_WORK (&TaskList, TaskList_size) ;
    GB_OK (GB_AxB_dot3_slice (&TaskList, &TaskList_size, &ntasks, &nthreads,
        C, Context)) ;

    GBURBLE ("nthreads %d ntasks %d ", nthreads, ntasks) ;

    //--------------------------------------------------------------------------
    // C<M> = A'*B, via masked dot product method and built-in semiring
    //--------------------------------------------------------------------------

    if (C_iso)
    { 

        //----------------------------------------------------------------------
        // C is iso; compute the pattern of C<M>=A'*B with the any_pair semiring
        //----------------------------------------------------------------------

        memcpy (C->x, cscalar, ctype->size) ;
        GB_OK (GB (_Adot3B__any_pair_iso) (C, M, Mask_struct, A, B,
            TaskList, ntasks, nthreads)) ;

    }
    else
    {

        //----------------------------------------------------------------------
        // C is non-iso
        //----------------------------------------------------------------------

        bool done = false ;

        #ifndef GBCUDA_DEV

            //------------------------------------------------------------------
            // define the worker for the switch factory
            //------------------------------------------------------------------

            #define GB_Adot3B(add,mult,xname) \
                GB (_Adot3B_ ## add ## mult ## xname)

            #define GB_AxB_WORKER(add,mult,xname)                           \
            {                                                               \
                info = GB_Adot3B (add,mult,xname) (C, M, Mask_struct, A, B, \
                    TaskList, ntasks, nthreads) ;                           \
                done = (info != GrB_NO_VALUE) ;                             \
            }                                                               \
            break ;

            //------------------------------------------------------------------
            // launch the switch factory
            //------------------------------------------------------------------

            GB_Opcode mult_binop_code, add_binop_code ;
            GB_Type_code xcode, ycode, zcode ;
            if (GB_AxB_semiring_builtin (A, A_is_pattern, B, B_is_pattern,
                semiring, flipxy, &mult_binop_code, &add_binop_code, &xcode,
                &ycode, &zcode))
            { 
                #include "GB_AxB_factory.c"
            }

        #endif

        //----------------------------------------------------------------------
        // C<M> = A'*B, via masked dot product method and typecasting
        //----------------------------------------------------------------------

        if (!done)
        { 
            #define GB_DOT3_GENERIC
            GB_BURBLE_MATRIX (C, "(generic C<M>=A'*B) ") ;
            #include "GB_AxB_dot_generic.c"
        }
    }


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_meta16_factory.c
            //------------------------------------------------------------------
            // both A and B are sparse
            //------------------------------------------------------------------

            #define GB_A_IS_SPARSE 1
            #define GB_A_IS_HYPER  0
            #define GB_A_IS_BITMAP 0
            #define GB_A_IS_FULL   0
            #define GB_B_IS_SPARSE 1
            #define GB_B_IS_HYPER  0
            #define GB_B_IS_BITMAP 0
            #define GB_B_IS_FULL   0
            #include "GB_meta16_methods.c"

// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_meta16_methods.c
    // dot product methods
    #if defined ( GB_DOT4 )
    #include "GB_AxB_dot4_template.c"
    #elif defined ( GB_DOT3_PHASE1 )
    #include "GB_AxB_dot3_phase1_template.c"
    #elif defined ( GB_DOT3_PHASE2 )
    #include "GB_AxB_dot3_template.c"
    #elif defined ( GB_DOT2 )
    #include "GB_AxB_dot2_template.c"


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_dot3_template.c
//------------------------------------------------------------------------------
// GB_AxB_dot3_template: C<M>=A'*B via dot products, where C is sparse/hyper
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// C and M are both sparse or hyper, and C->h is a copy of M->h.
// M is present, and not complemented.  It may be valued or structural.

{

    int tid ;
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) \
        reduction(+:nzombies)
    for (tid = 0 ; tid < ntasks ; tid++)
    {

        //----------------------------------------------------------------------
        // get the task descriptor
        //----------------------------------------------------------------------

        int64_t kfirst = TaskList [tid].kfirst ;
        int64_t klast  = TaskList [tid].klast ;
        int64_t pC_first = TaskList [tid].pC ;
        int64_t pC_last  = TaskList [tid].pC_end ;
        int64_t task_nzombies = 0 ;     // # of zombies found by this task

        //----------------------------------------------------------------------
        // compute all vectors in this task
        //----------------------------------------------------------------------

        for (int64_t k = kfirst ; k <= klast ; k++)
        {

            //------------------------------------------------------------------
            // get C(:,k) and M(:k)
            //------------------------------------------------------------------

            #if defined ( GB_MASK_SPARSE_AND_STRUCTURAL )
            // M and C are sparse
            const int64_t j = k ;
            #else
            // M and C are either both sparse or both hypersparse
            const int64_t j = GBH (Ch, k) ;
            #endif

            int64_t pC_start = Cp [k] ;
            int64_t pC_end   = Cp [k+1] ;
            if (k == kfirst)
            { 
                // First vector for task; may only be partially owned.
                pC_start = pC_first ;
                pC_end   = GB_IMIN (pC_end, pC_last) ;
            }
            else if (k == klast)
            { 
                // Last vector for task; may only be partially owned.
                pC_end   = pC_last ;
            }
            else
            { 
                // task completely owns this vector C(:,k).
            }

            //------------------------------------------------------------------
            // get B(:,j)
            //------------------------------------------------------------------

            #if GB_B_IS_HYPER
                // B is hyper: find B(:,j) using the B->Y hyper hash
                int64_t pB_start, pB_end ;
                GB_hyper_hash_lookup (Bp, B_Yp, B_Yi, B_Yx, B_hash_bits,
                    j, &pB_start, &pB_end) ;
            #elif GB_B_IS_SPARSE
                // B is sparse
                const int64_t pB_start = Bp [j] ;
                const int64_t pB_end = Bp [j+1] ;
            #else
                // B is bitmap or full
                const int64_t pB_start = j * vlen ;
            #endif

            #if (GB_B_IS_SPARSE || GB_B_IS_HYPER)
                const int64_t bjnz = pB_end - pB_start ;
                if (bjnz == 0)
                {
                    // no work to do if B(:,j) is empty, except for zombies
                    task_nzombies += (pC_end - pC_start) ;
                    for (int64_t pC = pC_start ; pC < pC_end ; pC++)
                    { 
                        // C(i,j) is a zombie
                        int64_t i = Mi [pC] ;
                        Ci [pC] = GB_FLIP (i) ;
                    }
                    continue ;
                }
                #if (GB_A_IS_SPARSE || GB_A_IS_HYPER)
                    // Both A and B are sparse; get first and last in B(:,j)
                    const int64_t ib_first = Bi [pB_start] ;
                    const int64_t ib_last  = Bi [pB_end-1] ;
                #endif
            #endif

            //------------------------------------------------------------------
            // C(:,j)<M(:,j)> = A(:,i)'*B(:,j)
            //------------------------------------------------------------------

            for (int64_t pC = pC_start ; pC < pC_end ; pC++)
            {

                //--------------------------------------------------------------
                // get C(i,j) and M(i,j)
                //--------------------------------------------------------------

                bool cij_exists = false ;
                GB_CIJ_DECLARE (cij) ;

                // get the value of M(i,j)
                int64_t i = Mi [pC] ;
                #if !defined ( GB_MASK_SPARSE_AND_STRUCTURAL )
                // if M is structural, no need to check its values
                if (GB_mcast (Mx, pC, msize))
                #endif
                { 

                    //----------------------------------------------------------
                    // the mask allows C(i,j) to be computed
                    //----------------------------------------------------------

                    #if GB_A_IS_HYPER
                    // A is hyper: find A(:,i) using the A->Y hyper hash
                    int64_t pA, pA_end ;
                    GB_hyper_hash_lookup (Ap, A_Yp, A_Yi, A_Yx, A_hash_bits,
                        i, &pA, &pA_end) ;
                    const int64_t ainz = pA_end - pA ;
                    if (ainz > 0)
                    #elif GB_A_IS_SPARSE
                    // A is sparse
                    int64_t pA = Ap [i] ;
                    const int64_t pA_end = Ap [i+1] ;
                    const int64_t ainz = pA_end - pA ;
                    if (ainz > 0)
                    #else
                    // A is bitmap or full
                    const int64_t pA = i * vlen ;
                    #endif
                    { 
                        // C(i,j) = A(:,i)'*B(:,j)
                        #include "GB_AxB_dot_cij.c"
                    }
                }

                if (!GB_CIJ_EXISTS)
                { 
                    // C(i,j) is a zombie
                    task_nzombies++ ;
                    Ci [pC] = GB_FLIP (i) ;
                }
            }
        }
        nzombies += task_nzombies ;
    }
}

#undef GB_A_IS_SPARSE
#undef GB_A_IS_HYPER
#undef GB_A_IS_BITMAP
#undef GB_A_IS_FULL
#undef GB_B_IS_SPARSE
#undef GB_B_IS_HYPER
#undef GB_B_IS_BITMAP
#undef GB_B_IS_FULL





// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_saxpy.c
//------------------------------------------------------------------------------
// GB_AxB_saxpy: compute C=A*B, C<M>=A*B, or C<!M>=A*B
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

#include "GB_mxm.h"
#include "GB_bitmap_AxB_saxpy.h"
#include "GB_stringify.h"

// TODO: allow bitmap multiply to work in-place as well

GrB_Info GB_AxB_saxpy               // C = A*B using Gustavson/Hash/Bitmap
(
    GrB_Matrix C,                   // output, static header
    GrB_Matrix C_in,                // original input matrix
    const GrB_Matrix M,             // optional mask matrix
    const bool Mask_comp,           // if true, use !M
    const bool Mask_struct,         // if true, use the only structure of M
    const GrB_BinaryOp accum,
    const GrB_Matrix A,             // input matrix A
    const GrB_Matrix B,             // input matrix B
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    bool *mask_applied,             // if true, then mask was applied
    bool *done_in_place,            // if true, C was computed in-place 
    const GrB_Desc_Value AxB_method,
    const int do_sort,              // if nonzero, try to sort in saxpy3
    GB_Context Context
)
{

        //----------------------------------------------------------------------
        // saxpy3: general-purpose Gustavson/Hash method, C is sparse/hyper
        //----------------------------------------------------------------------

        // C is sparse or hypersparse

        // This method allocates its own workspace, which is very small if the
        // Hash method is used.  The workspace for Gustavson's method is
        // larger, but saxpy3 selects that method only if the total work is
        // high enough so that the time to initialize the space.  C is sparse
        // or hypersparse.

        #ifdef GB_DEBUGIFY_DEFN
        GB_debugify_mxm (C_iso, C_sparsity, ztype, M,
            Mask_struct, Mask_comp, semiring, flipxy, A, B) ;
        #endif

        ASSERT (C_sparsity == GxB_HYPERSPARSE || C_sparsity == GxB_SPARSE) ;
        info = GB_AxB_saxpy3 (C, C_iso, cscalar, C_sparsity, M, Mask_comp,
            Mask_struct, A, B, semiring, flipxy, mask_applied, AxB_method,
            do_sort, Context) ;

// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_saxpy3.c

        //----------------------------------------------------------------------
        // generic saxpy3 method
        //----------------------------------------------------------------------

        if (!done)
        { 
            info = GB_AxB_saxpy_generic (C, M, Mask_comp, Mask_struct,
                M_in_place, A, A_is_pattern, B, B_is_pattern, semiring,
                flipxy, GB_SAXPY_METHOD_3,
                SaxpyTasks, ntasks, nfine, nthreads, do_sort,
                Context) ;
        }

// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_saxpy_generic.c

        //----------------------------------------------------------------------
        // generic semirings with standard multiply operators
        //----------------------------------------------------------------------

        GB_BURBLE_MATRIX (C, "(generic C=A*B) ") ;

        if (opcode == GB_FIRST_binop_code)
        {
            // t = A(i,k)
            // fmult is not used and can be NULL.  This is required for
            // GB_reduce_to_vector for user-defined types.
            ASSERT (!flipxy) ;
            ASSERT (B_is_pattern) ;
            if (saxpy_method == GB_SAXPY_METHOD_3)
            { 
                // C is sparse or hypersparse
                info = GB_AxB_saxpy3_generic_first 
                     (C, M, Mask_comp, Mask_struct, M_in_place,
                      A, A_is_pattern, B, B_is_pattern, semiring,
                      SaxpyTasks, ntasks, nfine, nthreads, do_sort,
                      Context) ;
            }
            else
            { 
                // C is bitmap or full
                info = GB_bitmap_AxB_saxpy_generic_first 
                     (C, M, Mask_comp, Mask_struct, M_in_place,
                      A, A_is_pattern, B, B_is_pattern, semiring,
                      NULL, 0, 0, 0, 0,
                      Context) ;
            }
        }

// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_saxpy_generic_method.c
GrB_Info GB_AXB_SAXPY_GENERIC_METHOD
(
    GrB_Matrix C,                   // any sparsity
    const GrB_Matrix M,
    bool Mask_comp,
    const bool Mask_struct,
    const bool M_in_place,          // ignored if C is bitmap
    const GrB_Matrix A,
    bool A_is_pattern,
    const GrB_Matrix B,
    bool B_is_pattern,
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    // for saxpy3 only:
    GB_saxpy3task_struct *restrict SaxpyTasks, // NULL if C is bitmap
    int ntasks,
    int nfine,
    int nthreads,
    const int do_sort,              // if true, sort in saxpy3
    GB_Context Context
)
{
        #if OP_IS_FIRST
        { 
            // t = A(i,k)
            ASSERT (B_is_pattern) ;
            #undef  GB_MULT
            #define GB_MULT(t, aik, bkj, i, k, j) memcpy (t, aik, csize)
            #include "GB_AxB_saxpy_generic_template.c"
        }


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_saxpy_generic_template.c
{
    #if C_IS_SPARSE_OR_HYPERSPARSE            
    {
        // C is sparse or hypersparse
        ASSERT (GB_IS_SPARSE (C) || GB_IS_HYPERSPARSE (C)) ;
        if (M == NULL)
        { 
            // C = A*B, no mask
            #define GB_NO_MASK 1
            #define GB_MASK_COMP 0
            #include "GB_AxB_saxpy3_template.c"
        }
        else if (!Mask_comp)
        { 
            // C<M> = A*B
            #define GB_NO_MASK 0
            #define GB_MASK_COMP 0
            #include "GB_AxB_saxpy3_template.c"
        }
        else
        { 
            // C<!M> = A*B
            #define GB_NO_MASK 0
            #define GB_MASK_COMP 1
            #include "GB_AxB_saxpy3_template.c"
        }
    }
    #else
    {
        // C is bitmap or full
        #include "GB_bitmap_AxB_saxpy_template.c"
    }
    #endif
}

// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_saxpy3_template.c

                    #include "GB_AxB_saxpy3_fineHash_M_phase2.c"

    //==========================================================================
    // phase3/phase4: count nnz(C(:,j)) for fine tasks, cumsum of Cp
    //==========================================================================

    GB_AxB_saxpy3_cumsum (C, SaxpyTasks, nfine, chunk, nthreads, Context) ;



// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_saxpy3_fineHash_M_phase2.c
//------------------------------------------------------------------------------
// GB_AxB_saxpy3_fineHash_M_phase2: C<M>=A*B, fine Hash, phase2
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

{
    //--------------------------------------------------------------------------
    // phase2: fine hash task, C(:,j)<M(:,j)>=A*B(:,j)
    //--------------------------------------------------------------------------

    // M is sparse and scattered into Hf

    // Given Hf [hash] split into (h,f)

    // h == 0  , f == 0 : unlocked, unoccupied. C(i,j) ignored
    // h == i+1, f == 1 : unlocked, occupied by M(i,j)=1.
    //                    C(i,j) has not been seen.
    //                    Hx is not initialized.
    // h == i+1, f == 2 : unlocked, occupied by C(i,j), M(i,j)=1
    //                    Hx is initialized.
    // h == ..., f == 3 : locked.

    // 0 -> 0 : to ignore, if M(i,j)=0
    // 1 -> 3 : to lock, if i seen for first time
    // 2 -> 3 : to lock, if i seen already
    // 3 -> 2 : to unlock; now i has been seen

    GB_GET_M_j_RANGE (16) ;     // get first and last in M(:,j)
    for ( ; pB < pB_end ; pB++)     // scan B(:,j)
    { 
        GB_GET_B_kj_INDEX ;         // get index k of B(k,j)
        GB_GET_A_k ;                // get A(:,k)
        if (aknz == 0) continue ;
        GB_GET_B_kj ;               // bkj = B(k,j)
        #define GB_IKJ                                                        \
        {                                                                     \
            GB_MULT_A_ik_B_kj ;      /* t = A(i,k) * B(k,j) */                \
            int64_t i1 = i + 1 ;     /* i1 = one-based index */               \
            int64_t i_unlocked = (i1 << 2) + 2 ;  /* (i+1,2) */               \
            for (GB_HASH (i))        /* find i in hash table */               \
            {                                                                 \
                int64_t hf ;                                                  \
                GB_ATOMIC_READ                                                \
                hf = Hf [hash] ;        /* grab the entry */                  \
                if (GB_HAS_ATOMIC && (hf == i_unlocked))                      \
                {                                                             \
                    /* Hx [hash] += t */                                      \
                    GB_ATOMIC_UPDATE_HX (hash, t) ;                           \
                    break ;     /* C(i,j) has been updated */                 \
                }                                                             \
                if (hf == 0) break ; /* M(i,j)=0; ignore Cij */               \
                if ((hf >> 2) == i1) /* if true, i found */                   \
                {                                                             \
                    do /* lock the entry */                                   \
                    {                                                         \
                        /* do this atomically: */                             \
                        /* { hf = Hf [hash] ; Hf [hash] |= 3 ; }*/            \
                        GB_ATOMIC_CAPTURE_INT64_OR (hf, Hf [hash], 3) ;       \
                    } while ((hf & 3) == 3) ; /* own: f=1,2 */                \
                    if ((hf & 3) == 1) /* f == 1 */                           \
                    {                                                         \
                        /* C(i,j) is a new entry in C(:,j) */                 \
                        GB_ATOMIC_WRITE_HX (hash, t) ; /* Hx [hash] = t */    \
                    }                                                         \
                    else /* f == 2 */                                         \
                    {                                                         \
                        /* C(i,j) already appears in C(:,j) */                \
                        GB_ATOMIC_UPDATE_HX (hash, t) ; /* Hx [hash] += t */  \
                    }                                                         \
                    GB_ATOMIC_WRITE                                           \
                    Hf [hash] = i_unlocked ; /* unlock entry */               \
                    break ;                                                   \
                }                                                             \
            }                                                                 \
        }
        GB_SCAN_M_j_OR_A_k (A_ok_for_binary_search) ;
        #undef GB_IKJ
    }
}





// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_dot_generic.c
//------------------------------------------------------------------------------
// GB_AxB_dot_generic: generic template for all dot-product methods
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2022, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// This template serves the dot2 and dot3 methods, but not dot4.  The
// #including file defines GB_DOT2_GENERIC or GB_DOT3_GENERIC.



////////////////////////////////////////////////////////////////
//// for GB_AxB_dot4.c, not using mask
////////////////////////////////////////////////////////////////



// /Users/peng599/pppp/CLion/GraphBLAS/Source/GB_AxB_dot4.c
    //--------------------------------------------------------------------------
    // define the worker for the switch factory
    //--------------------------------------------------------------------------

    info = GrB_NO_VALUE ;

    #define GB_Adot4B(add,mult,xname) GB (_Adot4B_ ## add ## mult ## xname)
    #define GB_AxB_WORKER(add,mult,xname)                           \
    {                                                               \
        info = GB_Adot4B (add,mult,xname) (C, A, A_slice, naslice,  \
            B, B_slice, nbslice, nthreads, Context) ;               \
    }                                                               \
    break ;

    //--------------------------------------------------------------------------
    // launch the switch factory
    //--------------------------------------------------------------------------

    // disabled the ANY monoid
    #define GB_NO_ANY_MONOID
    #include "GB_AxB_factory.c"



// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_factory.c
// This switch factory is not used to call the ANY_PAIR iso semiring.

ASSERT (mult_binop_code != GB_ANY_binop_code) ;

{
    //--------------------------------------------------------------------------
    // launch the switch factory
    //--------------------------------------------------------------------------

    switch (mult_binop_code)
    {

        //----------------------------------------------------------------------
        case GB_FIRST_binop_code   :    // z = x
        //----------------------------------------------------------------------

            // 61 semirings with FIRST:
            // 50: (min,max,plus,times,any) for 10 non-boolean real
            // 5: (or,and,xor,eq,any) for boolean
            // 6: (plus,times,any) for 2 complex
            #define GB_MNAME _first
            #define GB_COMPLEX
            #include "GB_AxB_type_factory.c"
            break ;

        //----------------------------------------------------------------------
        case GB_SECOND_binop_code  :    // z = y
        //----------------------------------------------------------------------

            // 61 semirings with SECOND:
            // 50: (min,max,plus,times,any) for 10 real non-boolean
            // 5: (or,and,xor,eq,any) for boolean
            // 6: (plus,times,any) for 2 complex
            #define GB_MNAME _second
            #define GB_COMPLEX
            #include "GB_AxB_type_factory.c"
            break ;

        //----------------------------------------------------------------------
        case GB_MIN_binop_code     :    // z = min(x,y)
        //----------------------------------------------------------------------

            // 50 semirings: (min,max,plus,times,any) for 10 real non-boolean
            // MIN == TIMES == AND for boolean
            #define GB_NO_BOOLEAN
            #define GB_MNAME _min
            #include "GB_AxB_type_factory.c"
            break ;



// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_type_factory.c
if (xcode != GB_BOOL_code)
{
    switch (add_binop_code)
    {

        // disable the MIN, MAX, ANY, and TIMES monoids for some multops
        #ifndef GB_NO_MIN_MAX_ANY_TIMES_MONOIDS

        case GB_MIN_binop_code:

            switch (xcode)
            {
                // 10 real, non-boolean types
                case GB_INT8_code   : GB_AxB_WORKER (_min, GB_MNAME, _int8  )
                case GB_INT16_code  : GB_AxB_WORKER (_min, GB_MNAME, _int16 )
                case GB_INT32_code  : GB_AxB_WORKER (_min, GB_MNAME, _int32 )
                case GB_INT64_code  : GB_AxB_WORKER (_min, GB_MNAME, _int64 )
                case GB_UINT8_code  : GB_AxB_WORKER (_min, GB_MNAME, _uint8 )
                case GB_UINT16_code : GB_AxB_WORKER (_min, GB_MNAME, _uint16)
                case GB_UINT32_code : GB_AxB_WORKER (_min, GB_MNAME, _uint32)
                case GB_UINT64_code : GB_AxB_WORKER (_min, GB_MNAME, _uint64)
                case GB_FP32_code   : GB_AxB_WORKER (_min, GB_MNAME, _fp32  )
                case GB_FP64_code   : GB_AxB_WORKER (_min, GB_MNAME, _fp64  )
                default: ;
            }
            break ;

        case GB_MAX_binop_code:

            switch (xcode)
            {
                // 10 real, non-boolean types
                case GB_INT8_code   : GB_AxB_WORKER (_max, GB_MNAME, _int8  )
                case GB_INT16_code  : GB_AxB_WORKER (_max, GB_MNAME, _int16 )
                case GB_INT32_code  : GB_AxB_WORKER (_max, GB_MNAME, _int32 )
                case GB_INT64_code  : GB_AxB_WORKER (_max, GB_MNAME, _int64 )
                case GB_UINT8_code  : GB_AxB_WORKER (_max, GB_MNAME, _uint8 )
                case GB_UINT16_code : GB_AxB_WORKER (_max, GB_MNAME, _uint16)
                case GB_UINT32_code : GB_AxB_WORKER (_max, GB_MNAME, _uint32)
                case GB_UINT64_code : GB_AxB_WORKER (_max, GB_MNAME, _uint64)
                case GB_FP32_code   : GB_AxB_WORKER (_max, GB_MNAME, _fp32  )
                case GB_FP64_code   : GB_AxB_WORKER (_max, GB_MNAME, _fp64  )
                default: ;
            }
            break ;

        case GB_TIMES_binop_code:

            switch (xcode)
            {
                // 10 real, non-boolean types, plus 2 complex
                case GB_INT8_code   : GB_AxB_WORKER (_times, GB_MNAME, _int8  )
                case GB_INT16_code  : GB_AxB_WORKER (_times, GB_MNAME, _int16 )
                case GB_INT32_code  : GB_AxB_WORKER (_times, GB_MNAME, _int32 )
                case GB_INT64_code  : GB_AxB_WORKER (_times, GB_MNAME, _int64 )
                case GB_UINT8_code  : GB_AxB_WORKER (_times, GB_MNAME, _uint8 )
                case GB_UINT16_code : GB_AxB_WORKER (_times, GB_MNAME, _uint16)
                case GB_UINT32_code : GB_AxB_WORKER (_times, GB_MNAME, _uint32)
                case GB_UINT64_code : GB_AxB_WORKER (_times, GB_MNAME, _uint64)
                case GB_FP32_code   : GB_AxB_WORKER (_times, GB_MNAME, _fp32  )
                case GB_FP64_code   : GB_AxB_WORKER (_times, GB_MNAME, _fp64  )
                #if defined ( GB_COMPLEX ) && !defined (GB_NO_NONATOMIC_MONOID)
                // the TIMES monoid is non-atomic for complex types
                case GB_FC32_code   : GB_AxB_WORKER (_times, GB_MNAME, _fc32  )
                case GB_FC64_code   : GB_AxB_WORKER (_times, GB_MNAME, _fc64  )
                #endif
                default: ;
            }
            break ;

        #ifndef GB_NO_ANY_MONOID
        case GB_ANY_binop_code:

            switch (xcode)
            {
                // 10 real, non-boolean types, plus 2 complex
                case GB_INT8_code   : GB_AxB_WORKER (_any, GB_MNAME, _int8  )
                case GB_INT16_code  : GB_AxB_WORKER (_any, GB_MNAME, _int16 )
                case GB_INT32_code  : GB_AxB_WORKER (_any, GB_MNAME, _int32 )
                case GB_INT64_code  : GB_AxB_WORKER (_any, GB_MNAME, _int64 )
                case GB_UINT8_code  : GB_AxB_WORKER (_any, GB_MNAME, _uint8 )
                case GB_UINT16_code : GB_AxB_WORKER (_any, GB_MNAME, _uint16)
                case GB_UINT32_code : GB_AxB_WORKER (_any, GB_MNAME, _uint32)
                case GB_UINT64_code : GB_AxB_WORKER (_any, GB_MNAME, _uint64)
                case GB_FP32_code   : GB_AxB_WORKER (_any, GB_MNAME, _fp32  )
                case GB_FP64_code   : GB_AxB_WORKER (_any, GB_MNAME, _fp64  )
                #if defined ( GB_COMPLEX ) && !defined (GB_NO_NONATOMIC_MONOID)
                // the ANY monoid is non-atomic for complex types
                case GB_FC32_code   : GB_AxB_WORKER (_any, GB_MNAME, _fc32  )
                case GB_FC64_code   : GB_AxB_WORKER (_any, GB_MNAME, _fc64  )
                #endif
                default: ;
            }
            break ;
        #endif
        #endif

        case GB_PLUS_binop_code:

            switch (xcode)
            {
                // 10 real, non-boolean types, plus 2 complex
                case GB_INT8_code   : GB_AxB_WORKER (_plus, GB_MNAME, _int8  )
                case GB_INT16_code  : GB_AxB_WORKER (_plus, GB_MNAME, _int16 )
                case GB_INT32_code  : GB_AxB_WORKER (_plus, GB_MNAME, _int32 )
                case GB_INT64_code  : GB_AxB_WORKER (_plus, GB_MNAME, _int64 )
                case GB_UINT8_code  : GB_AxB_WORKER (_plus, GB_MNAME, _uint8 )
                case GB_UINT16_code : GB_AxB_WORKER (_plus, GB_MNAME, _uint16)
                case GB_UINT32_code : GB_AxB_WORKER (_plus, GB_MNAME, _uint32)
                case GB_UINT64_code : GB_AxB_WORKER (_plus, GB_MNAME, _uint64)
                case GB_FP32_code   : GB_AxB_WORKER (_plus, GB_MNAME, _fp32  )
                case GB_FP64_code   : GB_AxB_WORKER (_plus, GB_MNAME, _fp64  )
                #if defined ( GB_COMPLEX )
                // only the PLUS monoid is atomic for complex types
                case GB_FC32_code   : GB_AxB_WORKER (_plus, GB_MNAME, _fc32  )
                case GB_FC64_code   : GB_AxB_WORKER (_plus, GB_MNAME, _fc64  )
                #endif
                default: ;
            }
            break ;

        default: ;
    }
}


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Generated2/GB_AxB__plus_first_int32.c
    GrB_Info GB (_Adot4B__plus_first_int32)
    (
        GrB_Matrix C,
        const GrB_Matrix A, int64_t *restrict A_slice, int naslice,
        const GrB_Matrix B, int64_t *restrict B_slice, int nbslice,
        const int nthreads,
        GB_Context Context
    )
    { 
        #if GB_DISABLE
        return (GrB_NO_VALUE) ;
        #else
        #include "GB_AxB_dot4_meta.c"
        return (GrB_SUCCESS) ;
        #endif
    }


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_dot4_meta.c

    //--------------------------------------------------------------------------
    // C += A'*B
    //--------------------------------------------------------------------------

    #include "GB_meta16_factory.c"


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_meta16_factory.c

#define GB_META16

{
    if (A_is_sparse)
    {

        if (B_is_sparse)
        { 

            //------------------------------------------------------------------
            // both A and B are sparse
            //------------------------------------------------------------------

            #define GB_A_IS_SPARSE 1
            #define GB_A_IS_HYPER  0
            #define GB_A_IS_BITMAP 0
            #define GB_A_IS_FULL   0
            #define GB_B_IS_SPARSE 1
            #define GB_B_IS_HYPER  0
            #define GB_B_IS_BITMAP 0
            #define GB_B_IS_FULL   0
            #include "GB_meta16_methods.c"

        }
        else if (B_is_hyper)
        { 

            //------------------------------------------------------------------
            // A is sparse and B is hyper
            //------------------------------------------------------------------

            #define GB_A_IS_SPARSE 1
            #define GB_A_IS_HYPER  0
            #define GB_A_IS_BITMAP 0
            #define GB_A_IS_FULL   0
            #define GB_B_IS_SPARSE 0
            #define GB_B_IS_HYPER  1
            #define GB_B_IS_BITMAP 0
            #define GB_B_IS_FULL   0
            #include "GB_meta16_methods.c"

        }
        else if (B_is_bitmap)
        { 

            //------------------------------------------------------------------
            // A is sparse and B is bitmap
            //------------------------------------------------------------------

            #define GB_A_IS_SPARSE 1
            #define GB_A_IS_HYPER  0
            #define GB_A_IS_BITMAP 0
            #define GB_A_IS_FULL   0
            #define GB_B_IS_SPARSE 0
            #define GB_B_IS_HYPER  0
            #define GB_B_IS_BITMAP 1
            #define GB_B_IS_FULL   0
            #include "GB_meta16_methods.c"

        }
        else
        { 

            //------------------------------------------------------------------
            // A is sparse and B is full
            //------------------------------------------------------------------

            #define GB_A_IS_SPARSE 1
            #define GB_A_IS_HYPER  0
            #define GB_A_IS_BITMAP 0
            #define GB_A_IS_FULL   0
            #define GB_B_IS_SPARSE 0
            #define GB_B_IS_HYPER  0
            #define GB_B_IS_BITMAP 0
            #define GB_B_IS_FULL   1
            #include "GB_meta16_methods.c"

        }


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_meta16_methods.c

{

    // declare macros that depend on the sparsity of A and B
    #include "GB_meta16_definitions.h"

    // dot product methods
    #if defined ( GB_DOT4 )
    #include "GB_AxB_dot4_template.c"
    #elif defined ( GB_DOT3_PHASE1 )
    #include "GB_AxB_dot3_phase1_template.c"
    #elif defined ( GB_DOT3_PHASE2 )
    #include "GB_AxB_dot3_template.c"
    #elif defined ( GB_DOT2 )
    #include "GB_AxB_dot2_template.c"

    #else
    #error "method undefined"
    #endif

    // undefine the macros that define the A and B sparsity
    #undef GB_A_IS_SPARSE
    #undef GB_A_IS_HYPER
    #undef GB_A_IS_BITMAP
    #undef GB_A_IS_FULL
    #undef GB_B_IS_SPARSE
    #undef GB_B_IS_HYPER
    #undef GB_B_IS_BITMAP
    #undef GB_B_IS_FULL
}


// /Users/peng599/pppp/CLion/GraphBLAS/Source/Template/GB_AxB_dot2_template.c

#if ( !GB_A_IS_HYPER && !GB_B_IS_HYPER )
{

    //--------------------------------------------------------------------------
    // C=A'*B, C<M>=A'*B, or C<!M>=A'*B where C is bitmap
    //--------------------------------------------------------------------------

    int tid ;
    #if GB_C_IS_FULL
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
    #else
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) \
        reduction(+:cnvals)
    #endif
    for (tid = 0 ; tid < ntasks ; tid++)
    {

        //----------------------------------------------------------------------
        // get the task descriptor
        //----------------------------------------------------------------------

        const int a_tid = tid / nbslice ;
        const int b_tid = tid % nbslice ;
        const int64_t kA_start = A_slice [a_tid] ;
        const int64_t kA_end   = A_slice [a_tid+1] ;
        const int64_t kB_start = B_slice [b_tid] ;
        const int64_t kB_end   = B_slice [b_tid+1] ;
        #if (!GB_C_IS_FULL)
        int64_t task_cnvals = 0 ;
        #endif

        //----------------------------------------------------------------------
        // C=A'*B, C<M>=A'*B, or C<!M>=A'*B via dot products
        //----------------------------------------------------------------------

        for (int64_t j = kB_start ; j < kB_end ; j++)
        {

            //------------------------------------------------------------------
            // get C(:,j)
            //------------------------------------------------------------------

            const int64_t pC_start = j * cvlen ;

            //------------------------------------------------------------------
            // get B(:,j)
            //------------------------------------------------------------------

            #if GB_B_IS_SPARSE
                // B is sparse (never hypersparse)
                const int64_t pB_start = Bp [j] ;
                const int64_t pB_end = Bp [j+1] ;
                const int64_t bjnz = pB_end - pB_start ;
                if (bjnz == 0)
                { 
                    // no work to do if B(:,j) is empty, except to clear Cb
                    memset (&Cb [pC_start + kA_start], 0, kA_end - kA_start) ;
                    continue ;
                }
                #if GB_A_IS_SPARSE
                    // Both A and B are sparse; get first and last in B(:,j)
                    const int64_t ib_first = Bi [pB_start] ;
                    const int64_t ib_last  = Bi [pB_end-1] ;
                #endif
            #else
                // B is bitmap or full
                const int64_t pB_start = j * vlen ;
            #endif

            //------------------------------------------------------------------
            // C(:,j)<#M(:,j)> = A'*B(:,j), or C(:,j) = A'*B(:,j) if no mask
            //------------------------------------------------------------------

            for (int64_t i = kA_start ; i < kA_end ; i++)
            {

                //--------------------------------------------------------------
                // get C(i,j), M(i,j), and clear the C(i,j) bitmap
                //--------------------------------------------------------------

                int64_t pC = pC_start + i ;     // C is bitmap

                #if defined ( GB_ANY_SPECIALIZED )
                // M is bitmap and structural; Mask_comp true
                Cb [pC] = 0 ;
                if (!Mb [pC])
                #elif defined ( GB_MASK_IS_PRESENT )
                bool mij ;
                if (M_is_bitmap)
                { 
                    // M is bitmap
                    mij = Mb [pC] && GB_mcast (Mx, pC, msize) ;
                }
                else if (M_is_full)
                { 
                    // M is full
                    mij = GB_mcast (Mx, pC, msize) ;
                }
                else // M is sparse or hyper
                { 
                    // M has been scattered into the C bitmap
                    mij = (Cb [pC] > 1) ;
                }
                Cb [pC] = 0 ;
                if (mij ^ Mask_comp)
                #elif GB_C_IS_FULL
                // C is full; nothing to do
                #else
                // M is not present
                Cb [pC] = 0 ;
                #endif
                { 

                    //----------------------------------------------------------
                    // the mask allows C(i,j) to be computed
                    //----------------------------------------------------------

                    #if GB_A_IS_SPARSE
                        // A is sparse
                        int64_t pA = Ap [i] ;
                        const int64_t pA_end = Ap [i+1] ;
                        const int64_t ainz = pA_end - pA ;
                        #if (!GB_C_IS_FULL)
                        if (ainz > 0)       // skip this test if C is full
                        #endif
                    #else
                        // A is bitmap or full
                        #ifdef GB_A_NOT_TRANSPOSED
                        // A(i,:) starts at position i
                        const int64_t pA = i ;
                        #else
                        // A(:,i) starts at position i * vlen
                        const int64_t pA = i * vlen ;
                        #endif
                    #endif
                    { 
                        // C(i,j) = A(:,i)'*B(:,j) or A(i,:)*B(:,j)
                        bool cij_exists = false ;
                        GB_CIJ_DECLARE (cij) ;
                        #include "GB_AxB_dot_cij.c"
                    }
                }
            }
        }
        #if (!GB_C_IS_FULL)
        cnvals += task_cnvals ;
        #endif
    }
}
#endif