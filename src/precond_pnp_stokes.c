/*! \file precond_pnp_stokes.c
 *  \brief Preconditioners for pnp+stokes problems
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/


/**
 * \fn void fasp_precond_pnp_stokes_diag (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 */
void fasp_precond_pnp_stokes_diag (REAL *r,
                                REAL *z,
                                void *data)
{
    
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_pnp_csr->row;
    const INT N1 = A_stokes_csr->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_array_cp(N, r, tempr->val);
    fasp_array_set(N, z, 0.0);
    
    // prepare
#if  WITH_UMFPACK
    void **LU_diag = precdata->LU_diag;
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
#endif
    
    // Preconditioning pnp block
#if  WITH_UMFPACK
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
#endif
    
    // Preconditioning A11 block
#if  WITH_UMFPACK
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_stokes_csr, &r1, &z1, LU_diag[1], 0);
#endif
    
    // restore r
    fasp_array_cp(N, tempr->val, r);
    
}

/**
 * \fn void fasp_precond_pnp_stokes_lower (REAL *r, REAL *z, void *data)
 * \brief block lower triangular preconditioning (3x3 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 */
void fasp_precond_pnp_stokes_lower (REAL *r,
                                 REAL *z,
                                 void *data)
{
    
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_pnp_csr->row;
    const INT N1 = A_stokes_csr->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_array_cp(N, r, tempr->val);
    fasp_array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
    
    // Preconditioning pnp block
#if  WITH_UMFPACK
    void **LU_diag = precdata->LU_diag;
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
#endif
    
    // r1 = r1 - A3*z0
    fasp_blas_dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
    
    // Preconditioning stokes block
#if  WITH_UMFPACK
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_stokes_csr, &r1, &z1, LU_diag[1], 0);
#endif
    
    // restore r
    fasp_array_cp(N, tempr->val, r);
    
}

/**
 * \fn void fasp_precond_pnp_stokes_upper (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (3x3 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 */
void fasp_precond_pnp_stokes_upper (REAL *r,
                                 REAL *z,
                                 void *data)
{
    
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_pnp_csr->row;
    const INT N1 = A_pnp_csr->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_array_cp(N, r, tempr->val);
    fasp_array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
    
    // Preconditioning stokes block
#if  WITH_UMFPACK
    void **LU_diag = precdata->LU_diag;
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_stokes_csr, &r1, &z1, LU_diag[1], 0);
#endif
    
    // r1 = r1 - A5*z2
    fasp_blas_dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    
    // Preconditioning pnp block
#if  WITH_UMFPACK
    /* use UMFPACK direct solver */
    fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
#endif
    
    // restore r
    fasp_array_cp(N, tempr->val, r);
    
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
