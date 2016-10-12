/*! \file itsolver_pnp_stokes.c
 *  \brief Iterative solvers for pnp+stokes system (main file)
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"




/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_bdcsr_krylov_pnp_stokes (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 *
 */
INT fasp_solver_bdcsr_krylov_pnp_stokes (block_dCSRmat *A,
                                      dvector *b,
                                      dvector *x,
                                      itsolver_param *itparam,
                                      itsolver_param *itparam_pnp,
                                      AMG_param *amgparam_pnp,
                                      itsolver_ns_param *itparam_stokes,
                                      AMG_ns_param *amgparam_stokes)
{
    const SHORT prtlvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
    
    //const SHORT max_levels = amgparam->max_levels;
    INT m, n, nnz, i;
    
    // local variables
    dCSRmat A_pnp_csr;
    dCSRmat A_stokes_csr;
    
    //AMG_data **mgl = NULL;
    
#if WITH_UMFPACK
    void **LU_diag = (void **)fasp_mem_calloc(2, sizeof(void *));
#endif
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    /* setup preconditioner */
    fasp_gettime(&setup_start);
    
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
        
#if WITH_UMFPACK
        // Need to sort the diagonal blocks for UMFPACK format
        
        // pnp block
        A_pnp_csr = fasp_dcsr_create(A->blocks[0]->row, A->blocks[0]->col, A->blocks[0]->nnz);
        fasp_dcsr_transz(A->blocks[0], NULL, &A_pnp_csr);
        
        printf("Factorization for pnp diagonal block: \n");
        LU_diag[0] = fasp_umfpack_factorize(&A_pnp_csr, prtlvl);
        
        // stokes block
        A_stokes_csr = fasp_dcsr_create(A->blocks[3]->row, A->blocks[3]->col, A->blocks[3]->nnz);
        fasp_dcsr_transz(A->blocks[3], NULL, &A_stokes_csr);
        
        printf("Factorization for pnp diagonal block: \n");
        LU_diag[1] = fasp_umfpack_factorize(&A_stokes_csr, prtlvl);
        
#endif
        
    }
    
    /* diagonal blocks are solved by other method */

    
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    precond_pnp_stokes_data precdata;
    precdata.Abcsr = A;
    precdata.A_pnp_csr = &A_pnp_csr;
    precdata.A_stokes_csr = &A_stokes_csr;

    precdata.r = fasp_dvec_create(b->row);
    
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
#if WITH_UMFPACK
        precdata.LU_diag = LU_diag;
#endif
    }
    /* diagonal blocks are solved by AMG */

    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    precond prec; prec.data = &precdata;
    
    switch (precond_type)
    {
        case 21:
            prec.fct = fasp_precond_pnp_stokes_diag;
            break;
         
        case 22:
            prec.fct = fasp_precond_pnp_stokes_lower;
            break;
            
        case 23:
            prec.fct = fasp_precond_pnp_stokes_upper;
            break;
            
        default:
            fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
            break;
    }
    
    if ( prtlvl >= PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }
    
    
    // solver part
    fasp_gettime(&solver_start);
    
    status=fasp_solver_bdcsr_itsolver(A,b,x, &prec,itparam);
    //status=fasp_solver_bdcsr_itsolver(A,b,x, NULL,itparam);

    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( prtlvl >= PRINT_MIN )
        print_cputime("Krylov method totally", solver_duration);
    
    // clean
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
#if WITH_UMFPACK
        for (i=0; i<2; i++) fasp_umfpack_free_numeric(LU_diag[i]);
        
        fasp_dcsr_free(&A_pnp_csr);
        fasp_dcsr_free(&A_stokes_csr);
#endif
    }
    /* diagonal blocks are solved by AMG */
    
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}








/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
