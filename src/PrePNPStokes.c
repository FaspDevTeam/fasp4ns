/*! \file  PrePNPStokes.c
 *
 *  \brief Preconditioners for PNP+Stokes problems
 *
 *  \note  This file contains Level-4 (Pre) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--Present by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
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
 *
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
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
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
    fasp_darray_cp(N, tempr->val, r);
}

/**
 * \fn void fasp_precond_pnp_stokes_lower (REAL *r, REAL *z, void *data)
 *
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
    dBLCmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_pnp_csr->row;
    const INT N1 = A_stokes_csr->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
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
    fasp_darray_cp(N, tempr->val, r);
}

/**
 * \fn void fasp_precond_pnp_stokes_upper (REAL *r, REAL *z, void *data)
 *
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
    dBLCmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_pnp_csr->row;
    const INT N1 = A_pnp_csr->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
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
    fasp_darray_cp(N, tempr->val, r);
}

/**
 * \fn void fasp_precond_pnp_stokes_diag_inexact (REAL *r, REAL *z, void *data)
 *
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 */
void fasp_precond_pnp_stokes_diag_inexact (REAL *r,
                                           REAL *z,
                                           void *data)
{
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dBSRmat *A_pnp_bsr = precdata->A_pnp_bsr;
    dBLCmat *A_stokes_bcsr = precdata->A_stokes_bcsr;
    dvector *tempr = &(precdata->r);
    
    void **LU_diag = precdata->LU_diag;
    precond_data_bsr *precdata_pnp = precdata->precdata_pnp;
    precond_ns_data  *precdata_stokes = precdata->precdata_stokes;
    
    const INT N0 = A_pnp_bsr->ROW*A_pnp_bsr->nb;
    const INT N1 = A_stokes_bcsr->blocks[0]->row + A_stokes_bcsr->blocks[2]->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
    
    // Preconditioning pnp block
    //precond prec_pnp;
    
    //prec_pnp.data = precdata_pnp;
    //prec_pnp.fct = precdata->pnp_fct;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_ilu;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_ilu;
    
    //prec_pnp.data = precdata->diag_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_diag;
    
    //fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
    fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, NULL, 1e-3,1e-8, 50, 50, 1, 0);
    //fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    //fasp_solver_dcsr_pvgmres(A->blocks[0], &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    
    // Preconditioning A11 block
    precond prec_stokes;
    prec_stokes.data = precdata_stokes;
    prec_stokes.fct = precdata->stokes_fct;
    
    fasp_solver_dblc_pvfgmres(A_stokes_bcsr, &r1, &z1, &prec_stokes, 1e-3,1e-8, 100, 100, 1, 0);
    
    // restore r
    fasp_darray_cp(N, tempr->val, r);
}

/**
 * \fn void fasp_precond_pnp_stokes_lower_inexact (REAL *r, REAL *z, void *data)
 *
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
void fasp_precond_pnp_stokes_lower_inexact (REAL *r,
                                            REAL *z,
                                            void *data)
{
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    dBLCmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dBSRmat *A_pnp_bsr = precdata->A_pnp_bsr;
    //dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    dBLCmat *A_stokes_bcsr = precdata->A_stokes_bcsr;
    
    dvector *tempr = &(precdata->r);
    
    void **LU_diag = precdata->LU_diag;
    precond_data_bsr *precdata_pnp= precdata->precdata_pnp;
    precond_ns_data  *precdata_stokes = precdata->precdata_stokes;
    
    const INT N0 = A_pnp_bsr->ROW*A_pnp_bsr->nb;
    const INT N1 = A_stokes_bcsr->blocks[0]->row + A_stokes_bcsr->blocks[2]->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
    
    // Preconditioning pnp block
    //precond prec_pnp;
    
    //prec_pnp.data = precdata_pnp;
    //prec_pnp.fct = precdata->pnp_fct;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_ilu;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_ilu;
    
    //prec_pnp.data = precdata->diag_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_diag;
    
    //fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
    fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, NULL, 1e-3,1e-8,50, 50, 1, 0);
    //fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    //fasp_solver_dcsr_pvgmres(A->blocks[0], &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    
    // r1 = r1 - A3*z0
    fasp_blas_dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
    
    // Preconditioning stokes block
    precond prec_stokes;
    prec_stokes.data = precdata_stokes;
    prec_stokes.fct = precdata->stokes_fct;
    
    fasp_solver_dblc_pvfgmres(A_stokes_bcsr, &r1, &z1, &prec_stokes, 1e-3, 1e-8,100, 100, 1, 0);
    
    // restore r
    fasp_darray_cp(N, tempr->val, r);
}

/**
 * \fn void fasp_precond_pnp_stokes_upper_inexact (REAL *r, REAL *z, void *data)
 *
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
void fasp_precond_pnp_stokes_upper_inexact (REAL *r,
                                            REAL *z,
                                            void *data)
{
    precond_pnp_stokes_data *precdata=(precond_pnp_stokes_data *)data;
    dBLCmat *A = precdata->Abcsr;
    dCSRmat *A_pnp_csr = precdata->A_pnp_csr;
    dBSRmat *A_pnp_bsr = precdata->A_pnp_bsr;
    //dCSRmat *A_stokes_csr = precdata->A_stokes_csr;
    dBLCmat *A_stokes_bcsr = precdata->A_stokes_bcsr;
    
    dvector *tempr = &(precdata->r);
    
    void **LU_diag = precdata->LU_diag;
    precond_data_bsr *precdata_pnp= precdata->precdata_pnp;
    precond_ns_data  *precdata_stokes = precdata->precdata_stokes;
    
    const INT N0 = A_pnp_bsr->ROW*A_pnp_bsr->nb;
    const INT N1 = A_stokes_bcsr->blocks[0]->row + A_stokes_bcsr->blocks[2]->row;
    const INT N = N0 + N1;
    
    // back up r, setup z;
    fasp_darray_cp(N, r, tempr->val);
    fasp_darray_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, z0, z1;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    
    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);
    
    // Preconditioning stokes block
    precond prec_stokes;
    prec_stokes.data = precdata_stokes;
    prec_stokes.fct = precdata->stokes_fct;
    
    fasp_solver_dblc_pvfgmres(A_stokes_bcsr, &r1, &z1, &prec_stokes, 1e-3, 1e-8,100, 100, 1, 0);
    
    // r1 = r1 - A5*z2
    fasp_blas_dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    
    // Preconditioning pnp block
    //precond prec_pnp;
    
    //prec_pnp.data = precdata_pnp;
    //prec_pnp.fct = precdata->pnp_fct;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_ilu;
    
    //prec_pnp.data = precdata->ILU_pnp;
    //prec_pnp.fct = fasp_precond_ilu;
    
    //prec_pnp.data = precdata->diag_pnp;
    //prec_pnp.fct = fasp_precond_dbsr_diag;
    
    //fasp_umfpack_solve(A_pnp_csr, &r0, &z0, LU_diag[0], 0);
    //fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, NULL, 1e-3, 50, 50, 1, 1);
    //fasp_solver_dbsr_pvgmres(A_pnp_bsr, &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    fasp_solver_dcsr_pvgmres(A->blocks[0], &r0, &z0, NULL, 1e-3, 1e-8,50, 50, 1, 0);
    //fasp_solver_dcsr_pvgmres(A->blocks[0], &r0, &z0, &prec_pnp, 1e-3, 100, 100, 1, 1);
    
    // restore r
    fasp_darray_cp(N, tempr->val, r);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
