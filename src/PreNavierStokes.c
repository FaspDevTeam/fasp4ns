/*! \file  PreNavierStokes.c
 *
 *  \brief Preconditioners for (Navier-)Stokes problems
 *
 *  \note  This file contains Level-4 (Pre) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--Present by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  TODO: Clean up these functions! --Chensong
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

#define USE_EXACT_SOLVE OFF

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_ns_bdiag (REAL *r, REAL *z, void *data)
 *
 * \brief Block diagonal preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu, Lu Wang
 * \date   10/20/2013
 *
 * Modified by Lu Wang on 02/12/2014
 * Modified by Xiaozhe Hu on 02/21/2014
 * Modified by Xiaozhe Hu on 05/27/2014
 */
void fasp_precond_ns_bdiag (REAL *r,
                            REAL *z,
                            void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    
    const INT col = predata->col, colA = predata->colA, colB = predata->colB;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // setup z;
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
    // prepare AMG preconditioner
    AMG_data  *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    const REAL  tolv = itparam_v->tol;
    const REAL  abstolv = itparam_v->abstol;
    const INT   maxitv = itparam_v->maxit, restartv = itparam_v->restart;
    const SHORT prtlvlv = itparam_v->print_level;

#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v, amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A,&rv,&zv,&pc_v,tolv,abstolv,maxitv,restartv,1,prtlvlv);
    
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif

    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    const REAL  tolp = itparam_p->tol;
    const REAL  abstolp = itparam_p->abstol;
    const INT   maxitp = itparam_p->maxit, restartp = itparam_p->restart;
    const SHORT prtlvlp = itparam_p->print_level;

#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;

        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_s,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 2){
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A,&rs,&zs,&pc_p,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 4) {
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_ilu,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }
    
#else
    
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    if ( prtlvlv > 0 ) printf(COLOR_GREEN "\n");
}

/**
 * \fn void fasp_precond_ns_low_btri (double *r, double *z, void *data)
 *
 * \brief Upper block diagonal preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu, Lu Wang
 * \date   10/20/2013
 *
 * Modified by Lu Wang on 02/11/2014
 * Modified by Xiaozhe Hu on 02/21/2014
 * Modified by Xiaozhe Hu on 05/27/2014
 */
void fasp_precond_ns_low_btri (REAL *r,
                               REAL *z,
                               void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    const REAL  tolv = itparam_v->tol;
    const REAL  abstolv = itparam_v->abstol;
    const INT   maxitv = itparam_v->maxit, restartv = itparam_v->restart;
    const SHORT prtlvlv = itparam_v->print_level;

    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A,&rv,&zv,&pc_v,tolv,abstolv,maxitv,restartv,1,prtlvlv);

#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    const REAL  tolp = itparam_p->tol;
    const REAL  abstolp = itparam_p->abstol;
    const INT   maxitp = itparam_p->maxit, restartp = itparam_p->restart;
    const SHORT prtlvlp = itparam_p->print_level;

#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_s,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 2){
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A,&rs,&zs,&pc_p,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 4) {
        ILU_data *LU_p = predata->ILU_p;
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_ilu,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }
    
#else

    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    fasp_darray_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_up_btri (double *r, double *z, void *data)
 *
 * \brief Upper block diagonal preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date   10/20/2013
 *
 * Modified by Lu Wang on 02/11/2014
 * Modified by Xiaozhe Hu on 02/21/2014
 * Modified by Xiaozhe Hu on 05/27/2014
 */
void fasp_precond_ns_up_btri (REAL *r,
                              REAL *z,
                              void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    const REAL  tolp = itparam_p->tol;
    const REAL  abstolp = itparam_p->abstol;
    const INT   maxitp = itparam_p->maxit, restartp = itparam_p->restart;
    const SHORT prtlvlp = itparam_p->print_level;

#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_s,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 2) {
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A,&rs,&zs,&pc_p,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }

    else if (itparam_p->precond_type == 4) {
        ILU_data *LU_p = predata->ILU_p;
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S,&rs,&zs,&pc_ilu,tolp,abstolp,maxitp,restartp,1,prtlvlp);
    }
    
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
    // Solve velocity
    //-------------------
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    const REAL  tolv = itparam_v->tol;
    const REAL  abstolv = itparam_v->abstol;
    const INT   maxitv = itparam_v->maxit, restartv = itparam_v->restart;
    const SHORT prtlvlv = itparam_v->print_level;

#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v,tolv,abstolv,maxitv,restartv,1,prtlvlv);
    
#else

    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    // restore r
    fasp_darray_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_blu (double *r, double *z, void *data)
 *
 * \brief Block LU preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   04/15/2016
 */
void fasp_precond_ns_blu (REAL *r,
                          REAL *z,
                          void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol, itparam_v->abstol, itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, &zs, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    fasp_darray_set(colA, zv.val, 0.0);
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol,itparam_v->abstol, itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    // restore r
    fasp_darray_cp(col, tempr, r);
    
}

/**
 * \fn void fasp_precond_ns_simple (double *r, double *z, void *data)
 *
 * \brief SIMPLE block diagonal preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   04/22/2016
 */
void fasp_precond_ns_simple (REAL *r,
                             REAL *z,
                             void *data)
{
    
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    INT i;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol,itparam_v->abstol, itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, &zs, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute zu = zu - D^{-1}B^T zs
    //-------------------
    fasp_blas_dcsr_mxv(predata->Bt, zs.val, rv.val); // rv = B^T zs
    
    for (i=0;i<colA;i++)
    {
        if (predata->diag_A->val[i] > SMALLREAL) rv.val[i]   = rv.val[i]/predata->diag_A->val[i]; // rv = D^{-1}rv
    }
    
    fasp_blas_darray_axpy (colA, -1.0, rv.val, zv.val); // zu = zu - rv
    
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    //-------------------
    // restore r
    //-------------------
    fasp_darray_cp(col, tempr, r);
    
}

/**
 * \fn void fasp_precond_ns_simpler (double *r, double *z, void *data)
 *
 * \brief SIMPLER block preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   04/29/2016
 */
void fasp_precond_ns_simpler (REAL *r,
                              REAL *z,
                              void *data)
{
    
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    INT i;
    
    dvector *deltaS = predata->sp;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Compute rs = rs - B D^{-1} rv
    //-------------------
    for (i=0;i<colA;i++)
    {
        if (predata->diag_A->val[i] > SMALLREAL) zv.val[i]   = rv.val[i]/predata->diag_A->val[i]; // zv = D^{-1}rv
    }
    
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val); // rs = rs - B zv
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, &zs, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol, itparam_v->abstol,itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    
    //-------------------
    // restore r
    //-------------------
    fasp_darray_cp(col, tempr, r);
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, deltaS, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, deltaS, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, deltaS, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, deltaS, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute zs = zs + deltaS
    //-------------------
    fasp_blas_darray_axpy(colB, -1.0, deltaS->val, zs.val);
    
    //-------------------
    // Compute zu = zu - D^{-1}B^T deltaS
    //-------------------
    fasp_blas_dcsr_mxv(predata->Bt, deltaS->val, rv.val); // rv = B^T deltaS
    
    for (i=0;i<colA;i++)
    {
        if (predata->diag_A->val[i] > SMALLREAL) rv.val[i]   = rv.val[i]/predata->diag_A->val[i]; // rv = D^{-1}rv
    }
    
    fasp_blas_darray_axpy (colA, -1.0, rv.val, zv.val); // zu = zu - rv
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    //-------------------
    // restore r
    //-------------------
    fasp_darray_cp(col, tempr, r);
    
    // free
    //fasp_dvec_free (&deltaS);
    
}

/**
 * \fn void fasp_precond_ns_uzawa (double *r, double *z, void *data)
 *
 * \brief Uzawa preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/03/2016
 */
void fasp_precond_ns_uzawa (REAL *r,
                            REAL *z,
                            void *data)
{
    
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol,itparam_v->abstol, itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute B zv - rs
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------
    // Compute zs = omega*(-1)*(rs)
    //-------------------
    REAL omega = -1.0;
    fasp_blas_darray_axpy(colB, omega, rs.val,zs.val);
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    // restore r
    fasp_darray_cp(col, tempr, r);
    
}

/**
 * \fn void fasp_precond_ns_projection (double *r, double *z, void *data)
 *
 * \brief Projection preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/08/2016
 */
void fasp_precond_ns_projection (REAL *r,
                                 REAL *z,
                                 void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol,itparam_v->abstol,itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    fasp_darray_cp(colB, rs.val, predata->sp->val);
    
    //-------------------------
    // Solve Schur complement -B*B^T
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, &zs, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute zv = zv - B^T * zs
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, zv.val);
    
    //-------------------
    // Compute zs = sp
    //-------------------
    fasp_darray_cp(colB, predata->sp->val, zs.val);
    
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    // restore r
    fasp_darray_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_DGS (double *r, double *z, void *data)
 *
 * \brief DGS preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/13/2016
 */
void fasp_precond_ns_DGS (REAL *r,
                          REAL *z,
                          void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol, itparam_v->abstol,itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement -BB^T
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, predata->sp, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, predata->sp, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, predata->sp, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, predata->sp, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute zv = zv - B^T * sp
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, predata->sp->val, zv.val);
    
    //-------------------
    // Compute zs = zs + BB^T * sp
    //-------------------
    fasp_blas_dcsr_aAxpy(1.0, predata->S, predata->sp->val, zs.val);
    
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    // restore r
    fasp_darray_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_LSCDGS (double *r, double *z, void *data)
 *
 * \brief LSCDGS preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/13/2016
 */
void fasp_precond_ns_LSCDGS (REAL *r,
                             REAL *z,
                             void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl_v = predata->mgl_data_v;
    AMG_param *amgparam_v = predata->param_v;
    ITS_param *itparam_v = predata->ITS_param_v;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    //-------------------
    // Solve velocity
    //-------------------
#if USE_EXACT_SOLVE==OFF
    
    precond_data pcdata_v;
    fasp_param_amg_to_prec(&pcdata_v,amgparam_v);
    pcdata_v.max_levels = mgl_v[0].num_levels;
    pcdata_v.mgl_data = predata->mgl_data_v;
    precond pc_v; pc_v.data = &pcdata_v;
    pc_v.fct = fasp_precond_amg;
    
    if (itparam_v->print_level > 0) printf(COLOR_RESET "\n");
    
    fasp_solver_dcsr_pvfgmres(&mgl_v[0].A, &rv, &zv, &pc_v, itparam_v->tol, itparam_v->abstol,itparam_v->maxit, itparam_v->restart, 1, itparam_v->print_level);
#else
    
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl_v[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    ITS_param *itparam_p = predata->ITS_param_p;
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, predata->sp, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, predata->sp, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, predata->sp, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
    /*
     fasp_dcoo_write("Ap.dat",  predata->S);
     fasp_dvec_write("rp.dat", &rs);
     getchar();
     */
    
#else
    
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, predata->sp, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    // change the sign of the solution
    fasp_blas_darray_ax(predata->sp->row, -1.0, predata->sp->val);
    
    //-------------------
    // Compute zv = zv + B^T * sp
    //-------------------
    fasp_blas_dcsr_aAxpy(1.0, predata->Bt, predata->sp->val, zv.val);
    
    //-------------------
    // Compute rs = BAB^t * sp
    //-------------------
    fasp_blas_dcsr_mxv(predata->BABt, predata->sp->val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    
#if USE_EXACT_SOLVE==OFF
    
    if (itparam_p->precond_type == 1) {
        precond pc_s;
        pc_s.data = predata->diag_S;
        pc_s.fct  = fasp_precond_diag;
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 2){
        // prepare  AMG preconditioner for S
        AMG_data *mgl_p = predata->mgl_data_p;
        AMG_param *amgparam_p = predata->param_p;
        
        precond_data pcdata_p;
        fasp_param_amg_to_prec(&pcdata_p,amgparam_p);
        pcdata_p.max_levels = mgl_p[0].num_levels;
        pcdata_p.mgl_data = predata->mgl_data_p;
        precond pc_p; pc_p.data = &pcdata_p;
        pc_p.fct = fasp_precond_amg;
        
        fasp_solver_dcsr_pvfgmres(&mgl_p[0].A, &rs, &zs, &pc_p, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
    }
    else if (itparam_p->precond_type == 4) {
        
        ILU_data *LU_p = predata->ILU_p;
        
        precond pc_ilu;
        pc_ilu.data = LU_p;
        pc_ilu.fct  = fasp_precond_ilu;
        
        fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_ilu, itparam_p->tol,itparam_p->abstol,itparam_p->maxit, itparam_p->restart, 1, itparam_p->print_level);
        
    }
    
#else
    
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
    
#endif
    
    //-------------------
    // Compute zs = -zs
    //-------------------
    //fasp_blas_darray_ax(colB,-1.0,zs.val);
    
    
    if (itparam_v->print_level > 0) printf(COLOR_GREEN "\n");
    // restore r
    fasp_darray_cp(col, tempr, r);
}

#if 0
/**
 * \fn void fasp_precond_ns_sym_btri (double *r, double *z, void *data)
 *
 * \brief Upper block diagonal preconditioner for the NS equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 *
 * Xiaozhe Hu modified on 02/21/2014
 * Xiaozhe Hu modified on 04/15/2016
 */
void fasp_precond_ns_sym_btri (REAL *r,
                               REAL *z,
                               void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    //const int maxit = predata->maxit;
    //double *diagptr=predata->diag_S->val;
    
    // local variables
    double	*tempr = predata->w;
    
    // prepare	AMG preconditioner
    AMG_data *mgl = predata->mgl_data;
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation      = predata->relaxation;
    amgparam.coarse_scaling  = predata->coarse_scaling;
    amgparam.ILU_levels      = predata->mgl_data->ILU_levels;
    
    dvector rv; rv.row = colA; rv.val = r;
    dvector zv; zv.row = colA; zv.val = z;
    dvector rs; rs.row = colB; rs.val = r+colA;
    dvector zs; zs.row = colB; zs.val = z+colA;
    dvector z1v; z1v.row = colA; z1v.val = predata->w1;
    dvector z1s; z1s.row = colB; z1s.val = (predata->w1)+colA;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    //-------------------
    // Solve velocity
    //-------------------
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,&amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = predata->mgl_data;
    precond pc; pc.data = &pcdata;
    pc.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    precond pc_s;
    pc_s.data = predata->diag_S;
    pc_s.fct  = fasp_precond_diag;
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);
    
    //-------------------
    // Compute residule
    //-------------------
    //fasp_blas_dcsr_aAxpy(1.0, predata->C, zs.val, rs.val);
    
    //-------------------------
    // Solve Schur complement
    //-------------------------
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);
    
    //-------------------
    // Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
    // Solve velocity
    //-------------------
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    
    //fasp_blas_dvec_axpy(1.0,&z1v,&zv);
    //fasp_blas_dvec_axpy(1.0,&z1s,&zs);
    
    // restore r
    fasp_darray_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_lsc (REAL *r, REAL *z, void *data)
 *
 * \brief LSC preconditioner for ns equation
 *
 * \param r     pointer to residual
 * \param z     pointer to preconditioned residual
 * \param data  pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 */
void fasp_precond_ns_lsc (REAL *r,
                          REAL *z,
                          void *data)
{
    precond_ns_data *predata=(precond_ns_data *)data;
    const int col = predata->col, colA = predata->colA, colB = predata->colB;
    //const int maxit = predata->maxit;
    //double *diagptr=predata->diag_S->val;
    // local variables
    double	*tempr = predata->w;
    int i, status;
    //double tmp = 0;
    // prepare	 AMG preconditioner for A
    AMG_data *mgl = predata->mgl_data;
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.smooth_order = 1;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation      = predata->relaxation;
    amgparam.coarse_scaling  = predata->coarse_scaling;
    amgparam.tentative_smooth = 1.0;
    amgparam.ILU_levels      = predata->mgl_data->ILU_levels;
    
    // prepare iterative parameter
    ITS_param  itparam;fasp_param_solver_init (&itparam);
    itparam.print_level = 0;
    itparam.itsolver_type = SOLVER_VFGMRES;
    itparam.precond_type   = 3;
    itparam.maxit          = 100;
    itparam.restart        = 100;
    itparam.tol            = predata->amg_tol;
    
    // prepare	 AMG preconditioner for S
    //AMG_data *mgl_s = predata->mgl_data_p;
    AMG_param amgparam_s; fasp_param_amg_init(&amgparam_s);
    amgparam_s.cycle_type = predata->cycle_type;
    amgparam_s.smoother   = predata->smoother;
    amgparam_s.presmooth_iter  = predata->presmooth_iter;
    amgparam_s.postsmooth_iter = predata->postsmooth_iter;
    amgparam_s.relaxation      = predata->relaxation;
    amgparam_s.coarse_scaling  = predata->coarse_scaling;
    amgparam_s.ILU_levels      = predata->mgl_data_p->ILU_levels;
    
    // back up r, setup z;
    fasp_darray_cp(col, r, tempr);
    fasp_darray_set(col, z, 0.0);
    
    // Solve M
    //mgl_s->b.row=colB;
    //for (i=0;i<colB;++i) {
    //    mgl_s->b.val[i]=r[colA+i];
    //}
    //tmp = fasp_blas_dvec_norm2(&(mgl_s->b));
    // printf("norm of b=%e\n",tmp);
    //mgl_s->x.row=colB; fasp_dvec_set(colB,&mgl_s->x,0.0);
    //for (i=0;i<10;i++) fasp_solver_mgcycle(mgl_s,&amgparam_s);
    //for (i=0;i<colB;++i) {
    //    z[colA+i]=  mgl_s->x.val[i];
    //}
    
    dvector *diag_S = predata->diag_S;
    
    //fasp_blas_dcsr_mxv (predata->P,mgl_s->x.val,mgl_s->b.val);
    
    //for (i=0;i<10;i++) fasp_solver_mgcycle(mgl_s,&amgparam_s);
    
    
    for (i=0;i<colB;++i) {
        z[colA+i]=r[colA+i]*diag_S->val[i];
    }
    
    //for (i=0;i<colB;++i) {
    //	if (ABS(diagptr[i])>SMALLREAL) z[colA+i]=r[colA+i]/diagptr[i];
    //}
    
    //tmp  = fasp_blas_darray_norm2(col,z);
    //printf("norm of z2=%e\n",tmp);
    
    // Solve A by AMG
    mgl->b.row=colA; fasp_darray_cp(colA,r,mgl->b.val); // residual is an input
    mgl->x.row=colA; fasp_dvec_set(colA,&mgl->x,0.0);
    
    fasp_blas_dcsr_aAxpy(-1.0,predata->Bt,(z+colA),mgl->b.val);
    
    dCSRmat     *ptrA=&mgl[0].A;
    dvector     *b=&mgl[0].b, *x=&mgl[0].x; //*rr=&mgl[0].w;
    //fasp_dvec_set(colA,rr,0.0);
    //fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,rr->val);
    //double normb = fasp_blas_dvec_norm2(b);
    /*
     //for (i=0;i<100;++i)
     int iter = 0;
     while (++iter<200)
     {
     fasp_solver_mgcycle(mgl,&amgparam);
     fasp_dvec_cp(b,rr);
     fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,rr->val);
     double absres  = fasp_blas_dvec_norm2(rr); // residual ||r||
     
     if (absres<1e-1) break;
     if (absres/normb<1e-1) break;
     //printf("%d: res:%e\n",iter,absres);
     }
     */
    
    status = fasp_solver_dcsr_krylov_amg(ptrA, b, x, &itparam, &amgparam);
    //if (iter > 150) fasp_dvec_write("rhs.dat",b);
    for (i=0;i<colA;++i) {
        z[i]=mgl->x.val[i];
    }	
    
    // restore r
    fasp_darray_cp(col, tempr, r);
}
#endif 

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
