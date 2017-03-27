/*! \file  SolNavierStokes.c
 *
 *  \brief Iterative solvers for Navier-Stokes matrices (main file)
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         PreNavierStokes.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  // TODO: Fix Doxygen. --Chensong
 *  // TODO: Shorten function names. --Chensong
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

static inline void get_schur_diagA (dCSRmat *, dCSRmat *, dCSRmat *,
                                    dCSRmat *, dCSRmat *);
static inline void get_schur_massP (dCSRmat *, dCSRmat *, dCSRmat *,
                                    dCSRmat *, REAL, dCSRmat *);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_ns_solver_itsolver (dCSRmat *A, dvector *b, dvector *x,
 *                                    precond *prec, itsolver_ns_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods
 *
 * \param *A        pointer to the block dCSRmat matrix
 * \param *b        pointer to the dvector of right hand side
 * \param *x        pointer to the dvector of dofs
 * \param *prec     pointer to the preconditioner data
 * \param *itparam  pointer to parameters for iterative solvers
 * \return          the number of iterations
 *
 * \author Lu Wang
 * \date 03/02/2012
 *
 * Modified by Xiaozhe on 10/20/2013
 */
SHORT fasp_ns_solver_itsolver (dBLCmat *A,
                               dvector *b,
                               dvector *x,
                               precond *prec,
                               itsolver_ns_param *itparam)
{
    const SHORT PrtLvl     = itparam->print_level;
    const SHORT SolverType = itparam->itsolver_type;
    const SHORT StopType   = itparam->stop_type;
    const SHORT MaxIt      = itparam->maxit;
    const SHORT restart    = itparam->restart;
    const REAL  tol        = itparam->tol;
    
    clock_t solver_start=clock();
    INT iter;
    switch (SolverType) {
            
        case SOLVER_BiCGstab:
            if (PrtLvl>0) printf("BiCGstab method (Block CSR format) ...\n");
            iter=fasp_solver_dblc_pbcgs(A, b, x, prec, tol, MaxIt, StopType, PrtLvl);
            break;
            
        case SOLVER_MinRes:
            if (PrtLvl>0) printf("Calling MinRes solver (Block CSR format) ...\n");
            iter=fasp_solver_dblc_pminres(A, b, x, prec, tol, MaxIt, StopType, PrtLvl);
            break;
            
        case SOLVER_GMRES:
            if (PrtLvl>0) printf("Calling GMRES solver (Block CSR format) ...\n");
            iter=fasp_solver_dblc_pvgmres(A, b, x, prec, tol, MaxIt, restart, StopType, PrtLvl);
            break;
            
        case SOLVER_FGMRES:
            if (PrtLvl>0) printf("Calling FGMRES solver (Block CSR format) ...\n");
            iter=fasp_solver_dblc_pvfgmres(A, b, x, prec, tol, MaxIt, restart, StopType, PrtLvl);
            break;
            
        case SOLVER_GCR:
            if (PrtLvl>0) printf("Calling GCR solver (Block CSR format) ...\n");
            iter=fasp_solver_dblc_pgcr(A, b, x, prec, tol, MaxIt, restart, StopType, PrtLvl);
            break;
            
        default:
            printf("Error: wrong itertive solver type %d!\n", SolverType);
            iter = ERROR_SOLVER_TYPE;
    }
    
    if ((PrtLvl>1) && (iter >= 0)) {
        clock_t solver_end=clock();
        REAL solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
        printf("Iterative solver costs %f seconds.\n", solver_duration);
    }
    
    return iter;
}

/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes (dBLCmat *Mat, dvector *b, dvector *x,
 *                                                itsolver_ns_param *itparam,
 *                                                AMG_ns_param *amgparam,
 *                                                ILU_param *iluparam,
 *                                                SWZ_param *schparam)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param *A:	       pointer to the dBLCmat matrix
 * \param *b:	       pointer to the dvector of right hand side
 * \param *x:	       pointer to the dvector of dofs
 * \param *itparam:  pointer to parameters for iterative solvers
 * \param *precdata: pionter to preconditioner data for ns
 * \return           number of iterations
 *
 * \author Lu Wang
 * \date 03/02/2012
 *
 * Modified by Xiaozhe on 05/31/2016
 */
SHORT fasp_solver_dblc_krylov_navier_stokes (dBLCmat *Mat,
                                             dvector *b,
                                             dvector *x,
                                             itsolver_ns_param *itparam,
                                             AMG_ns_param *amgnsparam,
                                             ILU_param *iluparam,
                                             SWZ_param *schparam)
{
    
    // parameters
    const SHORT PrtLvl  = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    const INT schwarz_mmsize = schparam->SWZ_mmsize;
    const INT schwarz_maxlvl = schparam->SWZ_maxlvl;
    const INT schwarz_type   = schparam->SWZ_type;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    // Navier-Stokes 2 by 2 matrix
    dCSRmat *A  = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B  = Mat->blocks[2];
    dCSRmat *C  = Mat->blocks[3];
    const INT n = A->row, m = B->row, nnzA = A->nnz;
    
    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S;
    SWZ_data schwarz_data;
    dvector diag_A;
    dvector diag_S;
    dCSRmat BABt;
    
    // local variable
    clock_t solver_start, solver_end, setup_start, setup_end;
    REAL solver_duration, setup_duration;
    SHORT status=FASP_SUCCESS;
    
    //------ setup phase ------//
    setup_start = clock();
    
    //-----------------------//
    // setup AMG for velocity
    //-----------------------//
    
    AMG_data *mgl_v=fasp_amg_data_create(amgnsparam->param_v.max_levels);
    mgl_v[0].A=fasp_dcsr_create(n,n,nnzA);
    if (precond_type > 10) {
        
        dCSRmat BtB;
        fasp_blas_dcsr_mxm(Bt, B, &BtB);
        
        REAL gamma = 10;
        fasp_blas_dcsr_add (A, 1.0, &BtB, gamma, &mgl_v[0].A);
        
        fasp_dcsr_free(&BtB);
        
    }
    else {
        fasp_dcsr_cp(A,&mgl_v[0].A);
    }
    mgl_v[0].b=fasp_dvec_create(n); mgl_v[0].x=fasp_dvec_create(n);
    
    // setup AMG
    switch (amgnsparam->param_v.AMG_type) {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_v, &amgnsparam->param_v);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_v, &amgnsparam->param_v);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_v, &amgnsparam->param_v);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d!\n",amgnsparam->param_v.AMG_type);
            exit(ERROR_INPUT_PAR);
    }
    
    // get diagonal of A
    fasp_dcsr_getdiag(n, &mgl_v[0].A, &diag_A);
    
    //-------------------------//
    // setup Schur complement S
    //-------------------------//
    
    if (precond_type == 8 || precond_type == 9 || precond_type == 18 || precond_type == 19){
        fasp_blas_dcsr_mxm(B, Bt, &S);
        
        // change the sign of the BB^T
        fasp_blas_dcsr_axm(&S, -1.0);
        
        // make it non-singular
        INT i,k,j,ibegin,iend;
        
        for (i=0;i<S.row;++i) {
            ibegin=S.IA[i]; iend=S.IA[i+1];
            for (k=ibegin;k<iend;++k) {
                j=S.JA[k];
                if ((j-i)==0) {
                    S.val[k] = S.val[k] + 1e-8; break;
                } // end if
            } // end for k
        } // end for i
        
    }
    else if (precond_type == 10 || precond_type == 20){
        fasp_blas_dcsr_mxm(B, Bt, &S);
        fasp_blas_dcsr_rap(B, A, Bt, &BABt);
        
        // change the sign of the BB^T
        fasp_blas_dcsr_axm(&S, -1.0);
        
        // make it non-singular
        INT i,k,j,ibegin,iend;
        
        for (i=0;i<S.row;++i) {
            ibegin=S.IA[i]; iend=S.IA[i+1];
            for (k=ibegin;k<iend;++k) {
                j=S.JA[k];
                if ((j-i)==0) {
                    S.val[k] = S.val[k] + 1e-8; break;
                } // end if
            } // end for k
        } // end for i
        
    }
    else {
        get_schur_diagA(B,Bt,A,C,&S);
    }
    
    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);
    
    AMG_data *mgl_p;
    ILU_data LU_p;
    
    if (itparam->precond_type_p == 1){
        fasp_dcsr_getdiag(0,&S,&diag_S);
    }
    else if (itparam->precond_type_p == 2) {
        // Setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p=fasp_amg_data_create(amgnsparam->param_p.max_levels);
        mgl_p[0].A=fasp_dcsr_create(m,m,nnzS); fasp_dcsr_cp(As,&mgl_p[0].A);
        mgl_p[0].b=fasp_dvec_create(m); mgl_p[0].x=fasp_dvec_create(m);
        // setup AMG
        switch (amgnsparam->param_p.AMG_type) {
            case CLASSIC_AMG:
                fasp_amg_setup_rs(mgl_p, &amgnsparam->param_p);
                break;
            case SA_AMG:
                fasp_amg_setup_sa(mgl_p, &amgnsparam->param_p);
                break;
            case UA_AMG:
                fasp_amg_setup_ua(mgl_p, &amgnsparam->param_p);
                break;
            default:
                printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                       amgnsparam->param_p.AMG_type);
                exit(ERROR_INPUT_PAR);
        }
    }
    else if (itparam->precond_type_p == 4) {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }
    
    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itparam->print_level_v;
    ITS_param_v.itsolver_type = itparam->itsolver_type_v;
    ITS_param_v.restart = itparam->pre_restart_v;
    ITS_param_v.tol = itparam->pre_tol_v;
    ITS_param_v.maxit = itparam->pre_maxit_v;
    ITS_param_v.precond_type = itparam->precond_type_v;
    
    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itparam->print_level_p;
    ITS_param_p.itsolver_type = itparam->itsolver_type_p;
    ITS_param_p.restart = itparam->pre_restart_p;
    ITS_param_p.tol = itparam->pre_tol_p;
    ITS_param_p.maxit = itparam->pre_maxit_p;
    ITS_param_p.precond_type = itparam->precond_type_p;
    
    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;
    
    precdata.colA = n;
    precdata.colB = m;
    precdata.col  = n+m;
    precdata.M  = M;
    precdata.B  = B;
    precdata.Bt = Bt;
    precdata.C = C;
    precdata.BABt  = &BABt;
    
    precdata.param_v = &amgnsparam->param_v;
    precdata.param_p = &amgnsparam->param_p;
    precdata.ITS_param_v = &ITS_param_v;
    precdata.ITS_param_p = &ITS_param_p;
    precdata.mgl_data_v       = mgl_v;
    precdata.mgl_data_p       = mgl_p;
    precdata.ILU_p            = &LU_p;
    
    precdata.max_levels     = mgl_v[0].num_levels;
    precdata.print_level    = amgnsparam->param_v.print_level;
    precdata.maxit          = amgnsparam->param_v.maxit;
    precdata.amg_tol        = amgnsparam->param_v.tol;
    precdata.cycle_type     = amgnsparam->param_v.cycle_type;
    precdata.smoother       = amgnsparam->param_v.smoother;
    precdata.presmooth_iter = amgnsparam->param_v.presmooth_iter;
    precdata.postsmooth_iter= amgnsparam->param_v.postsmooth_iter;
    precdata.relaxation     = amgnsparam->param_v.relaxation;
    precdata.coarse_scaling = amgnsparam->param_v.coarse_scaling;
    
    precdata.diag_A = &diag_A;
    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;
    
    precdata.w = (REAL *)fasp_mem_calloc(precdata.col,sizeof(double));
    
    switch (precond_type) {
        case 1:
            prec.fct = fasp_precond_ns_bdiag;
            break;
        case 2:
            prec.fct = fasp_precond_ns_low_btri;
            break;
        case 3:
            prec.fct = fasp_precond_ns_up_btri;
            break;
        case 4:
            prec.fct = fasp_precond_ns_blu;
            break;
        case 5:
            prec.fct = fasp_precond_ns_simple;
            break;
        case 6:
            prec.fct = fasp_precond_ns_simpler;
            break;
        case 7:
            prec.fct = fasp_precond_ns_uzawa;
            break;
        case 8:
            prec.fct = fasp_precond_ns_projection;
            break;
        case 9:
            prec.fct = fasp_precond_ns_DGS;
            break;
        case 10:
            prec.fct = fasp_precond_ns_LSCDGS;
            break;
        case 11:
            prec.fct = fasp_precond_ns_bdiag;
            break;
        case 12:
            prec.fct = fasp_precond_ns_low_btri;
            break;
        case 13:
            prec.fct = fasp_precond_ns_up_btri;
            break;
        case 14:
            prec.fct = fasp_precond_ns_blu;
            break;
        case 15:
            prec.fct = fasp_precond_ns_simple;
            break;
        case 16:
            prec.fct = fasp_precond_ns_simpler;
            break;
        case 17:
            prec.fct = fasp_precond_ns_uzawa;
            break;
        case 18:
            prec.fct = fasp_precond_ns_projection;
            break;
        case 19:
            prec.fct = fasp_precond_ns_DGS;
            break;
        case 20:
            prec.fct = fasp_precond_ns_LSCDGS;
            break;
        default:
            printf("### ERROR: Unknown preconditioner type!\n");
            exit(ERROR_SOLVER_PRECTYPE);
    }
    
    setup_end = clock();
    
    if (PrtLvl>0) {
        setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
        printf("Setup costs %f.\n", setup_duration);
    }
    
    //------ solver phase ------//
    solver_start=clock();
    status=fasp_ns_solver_itsolver(Mat,b,x,&prec,itparam);
    //status=fasp_ns_solver_itsolver(Mat,b,x,NULL,itparam);
    solver_end=clock();
    
    if (PrtLvl>0) {
        solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
        printf(COLOR_RESET);
        printf("Solver costs %f seconds.\n", solver_duration);
        printf("Total costs %f seconds.\n", setup_duration + solver_duration);
    }
    
    //FINISHED:
    // clean up memory
    if (mgl_v) fasp_amg_data_free(mgl_v,&amgnsparam->param_v);
    if (itparam->precond_type_p == 1){fasp_dvec_free(&diag_S);}
    if (itparam->precond_type_p == 2) fasp_amg_data_free(mgl_p,&amgnsparam->param_p);
    
    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dcsr_free(&S);
    if (precond_type == 10 || precond_type == 20) fasp_dcsr_free(&BABt);
    fasp_dvec_free(&diag_A);
    
    return status;
}


/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes_with_pressure_mass (dBLCmat *Mat, dvector *b, dvector *x, itsolver_ns_param *itparam, AMG_ns_param *amgparam, ILU_param *iluparam, SWZ_param *schparam, dCSRmat *Mp)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param *A:	       pointer to the dBLCmat matrix
 * \param *b:	       pointer to the dvector of right hand side
 * \param *x:	       pointer to the dvector of dofs
 * \param *itparam:    pointer to parameters for iterative solvers
 * \param *precdata:   pionter to preconditioner data for ns
 * \param *Mp:         pointer to dCSRmat of the pressure mass matrix
 * \return           number of iterations
 *
 * \author Xiaozhe Hu
 * \date 017/07/2014
 *
 * \note In general, this is for purely Stokes problem, NS problem with div-div
 *       stablization -- Xiaozhe
 */
SHORT fasp_solver_dblc_krylov_navier_stokes_with_pressure_mass (dBLCmat *Mat,
                                                                dvector *b,
                                                                dvector *x,
                                                                itsolver_ns_param *itparam,
                                                                AMG_ns_param *amgnsparam,
                                                                ILU_param *iluparam,
                                                                SWZ_param *schparam,
                                                                dCSRmat *Mp)
{
    // parameters
    const SHORT PrtLvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    const INT schwarz_mmsize = schparam->SWZ_mmsize;
    const INT schwarz_maxlvl = schparam->SWZ_maxlvl;
    const INT schwarz_type   = schparam->SWZ_type;
    
    // Navier-Stokes 4 by 4 matrix
    dCSRmat *A  = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B  = Mat->blocks[2];
    dCSRmat *C  = Mat->blocks[3];
    const INT n = A->row, m = B->row, nnzA = A->nnz;
    
    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S,P;
    SWZ_data schwarz_data;
    dvector diag_S;
    
    // local variable
    clock_t solver_start, solver_end, setup_start, setup_end;
    REAL    solver_duration, setup_duration;
    SHORT   status = FASP_SUCCESS;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    //------ setup phase ------//
    setup_start = clock();
    
    //-----------------------//
    // setup AMG for velocity
    //-----------------------//
    
    AMG_data *mgl_v=fasp_amg_data_create(amgnsparam->param_v.max_levels);
    mgl_v[0].A=fasp_dcsr_create(n,n,nnzA); fasp_dcsr_cp(A,&mgl_v[0].A);
    mgl_v[0].b=fasp_dvec_create(n); mgl_v[0].x=fasp_dvec_create(n);
    
    // setup AMG
    switch (amgnsparam->param_v.AMG_type) {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_v, &amgnsparam->param_v);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_v, &amgnsparam->param_v);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_v, &amgnsparam->param_v);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d!\n",amgnsparam->param_v.AMG_type);
            exit(ERROR_INPUT_PAR);
    }
    
    //-------------------------//
    // setup Schur complement S using pressure mass
    //-------------------------//
    
    fasp_dcsr_alloc(Mp->row, Mp->col, Mp->nnz, &S);
    fasp_dcsr_cp(Mp, &S);
    
    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);
    
    AMG_data *mgl_p;
    ILU_data LU_p;
    
    if ( itparam->precond_type_p == 1 ) {
        fasp_dcsr_getdiag(0,&S,&diag_S);
    }
    else if ( itparam->precond_type_p == 2 ) {
        // setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p=fasp_amg_data_create(amgnsparam->param_p.max_levels);
        mgl_p[0].A=fasp_dcsr_create(m,m,nnzS); fasp_dcsr_cp(As,&mgl_p[0].A);
        mgl_p[0].b=fasp_dvec_create(m); mgl_p[0].x=fasp_dvec_create(m);
        // setup AMG
        switch ( amgnsparam->param_p.AMG_type ) {
            case CLASSIC_AMG:
                fasp_amg_setup_rs(mgl_p, &amgnsparam->param_p);
                break;
            case SA_AMG:
                fasp_amg_setup_sa(mgl_p, &amgnsparam->param_p);
                break;
            case UA_AMG:
                fasp_amg_setup_ua(mgl_p, &amgnsparam->param_p);
                break;
            default:
                printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                       amgnsparam->param_p.AMG_type);
                exit(ERROR_INPUT_PAR);
        }
    }
    else if ( itparam->precond_type_p == 4 ) {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }
    
    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itparam->print_level_v;
    ITS_param_v.itsolver_type = itparam->itsolver_type_v;
    ITS_param_v.restart = itparam->pre_restart_v;
    ITS_param_v.tol = itparam->pre_tol_v;
    ITS_param_v.maxit = itparam->pre_maxit_v;
    ITS_param_v.precond_type = itparam->precond_type_v;
    
    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itparam->print_level_p;
    ITS_param_p.itsolver_type = itparam->itsolver_type_p;
    ITS_param_p.restart = itparam->pre_restart_p;
    ITS_param_p.tol = itparam->pre_tol_p;
    ITS_param_p.maxit = itparam->pre_maxit_p;
    ITS_param_p.precond_type = itparam->precond_type_p;
    
    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;
    
    precdata.colA = n;
    precdata.colB = m;
    precdata.col  = n+m;
    precdata.M    = M;
    precdata.B    = B;
    precdata.Bt   = Bt;
    precdata.C    = C;
    
    precdata.param_v        = &amgnsparam->param_v;
    precdata.param_p        = &amgnsparam->param_p;
    precdata.ITS_param_v    = &ITS_param_v;
    precdata.ITS_param_p    = &ITS_param_p;
    precdata.mgl_data_v     = mgl_v;
    precdata.mgl_data_p     = mgl_p;
    precdata.ILU_p          = &LU_p;
    
    precdata.max_levels     = mgl_v[0].num_levels;
    precdata.print_level    = amgnsparam->param_v.print_level;
    precdata.maxit          = amgnsparam->param_v.maxit;
    precdata.amg_tol        = amgnsparam->param_v.tol;
    precdata.cycle_type     = amgnsparam->param_v.cycle_type;
    precdata.smoother       = amgnsparam->param_v.smoother;
    precdata.presmooth_iter = amgnsparam->param_v.presmooth_iter;
    precdata.postsmooth_iter= amgnsparam->param_v.postsmooth_iter;
    precdata.relaxation     = amgnsparam->param_v.relaxation;
    precdata.coarse_scaling = amgnsparam->param_v.coarse_scaling;
    
    
    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;
    
    precdata.w = (REAL *)fasp_mem_calloc(precdata.col,sizeof(double));
    
    switch (precond_type) {
        case 1:
            prec.fct = fasp_precond_ns_bdiag; break;
        case 2:
            prec.fct = fasp_precond_ns_low_btri; break;
        case 3:
            prec.fct = fasp_precond_ns_up_btri; break;
        case 4:
            prec.fct = fasp_precond_ns_blu; break;
        default:
            printf("### ERROR: Unknown preconditioner type!\n");
            exit(ERROR_SOLVER_PRECTYPE);
    }
    
    setup_end = clock();
    
    if (PrtLvl>0) {
        setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
        printf("Setup costs %f.\n", setup_duration);
    }
    
    //------ solver phase ------//
    solver_start=clock();
    status=fasp_ns_solver_itsolver(Mat,b,x,&prec,itparam);
    solver_end=clock();
    
    if (PrtLvl>0) {
        solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
        printf(COLOR_RESET);
        printf("Solver costs %f seconds.\n", solver_duration);
        printf("Total costs %f seconds.\n", setup_duration + solver_duration);
    }
    
    // clean up memory
    if (mgl_v) fasp_amg_data_free(mgl_v,&amgnsparam->param_v);
    if (itparam->precond_type_p == 1) fasp_dvec_free(&diag_S);
    if (itparam->precond_type_p == 2) fasp_amg_data_free(mgl_p,&amgnsparam->param_p);
    
    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dcsr_free(&S);
    
    return status;
}

/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes_schur_complement_with_pressure_mass (dBLCmat *Mat, dvector *b, dvector *x, itsolver_ns_param *itparam, AMG_ns_param *amgparam, ILU_param *iluparam, SWZ_param *schparam, dCSRmat *Mp)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param *A:	       pointer to the dBLCmat matrix
 * \param *b:	       pointer to the dvector of right hand side
 * \param *x:	       pointer to the dvector of dofs
 * \param *itparam:    pointer to parameters for iterative solvers
 * \param *precdata:   pionter to preconditioner data for ns
 * \param *Mp:         pointer to dCSRmat of the pressure mass matrix
 * \return           number of iterations
 *
 * \author Xiaozhe Hu
 * \date 017/07/2014
 *
 * \note In general, this is for NS problems without div-div stablization and pressure stablization (pressure block is zero), moreover, pressure mass matrix is provided -- Xiaozhe
 *
 */
SHORT fasp_solver_dblc_krylov_navier_stokes_schur_complement_with_pressure_mass (dBLCmat *Mat,
                                                                                 dvector *b,
                                                                                 dvector *x,
                                                                                 itsolver_ns_param *itparam,
                                                                                 AMG_ns_param *amgnsparam,
                                                                                 ILU_param *iluparam,
                                                                                 SWZ_param *schparam,
                                                                                 dCSRmat *Mp)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    // parameters
    const SHORT PrtLvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    const INT schwarz_mmsize = schparam->SWZ_mmsize;
    const INT schwarz_maxlvl = schparam->SWZ_maxlvl;
    const INT schwarz_type   = schparam->SWZ_type;
    
    // Navier-Stokes 4 by 4 matrix
    dCSRmat *A  = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B  = Mat->blocks[2];
    dCSRmat *C  = Mat->blocks[3];
    const INT n = A->row, m = B->row, nnzA = A->nnz;
    
    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S,P;
    SWZ_data schwarz_data;
    dvector diag_S;
    
    // local variable
    clock_t solver_start, solver_end, setup_start, setup_end;
    REAL solver_duration, setup_duration;
    SHORT status=FASP_SUCCESS;
    
    //------ setup phase ------//
    setup_start = clock();
    
    //-----------------------//
    // setup AMG for velocity
    //-----------------------//
    
    AMG_data *mgl_v=fasp_amg_data_create(amgnsparam->param_v.max_levels);
    
    mgl_v[0].A=fasp_dcsr_create(n,n,nnzA); fasp_dcsr_cp(A,&mgl_v[0].A);
    //get_schur_massP(Bt, B, Mp, A, 1e1, &mgl_v[0].A);
    
    mgl_v[0].b=fasp_dvec_create(n); mgl_v[0].x=fasp_dvec_create(n);
    
    // setup AMG
    switch (amgnsparam->param_v.AMG_type) {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_v, &amgnsparam->param_v);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_v, &amgnsparam->param_v);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_v, &amgnsparam->param_v);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d!\n",amgnsparam->param_v.AMG_type);
            exit(ERROR_INPUT_PAR);
    }
    //-------------------------//
    // setup Schur complement S using pressure mass
    //-------------------------//
    
    //fasp_dcsr_alloc(Mp->row, Mp->col, Mp->nnz, &S);
    //fasp_dcsr_cp(Mp, &S);
    
    get_schur_massP(B, Bt, &mgl_v[0].A, Mp, 1e5, &S);  // 1e5 is a parameter can be tuned, not sure how to tune this one now.... -- Xiaozhe
    
    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);
    
    AMG_data *mgl_p;
    ILU_data LU_p;
    
    if (itparam->precond_type_p == 1){
        fasp_dcsr_getdiag(0,&S,&diag_S);
    }
    else if (itparam->precond_type_p == 2) {
        // Setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p=fasp_amg_data_create(amgnsparam->param_p.max_levels);
        mgl_p[0].A=fasp_dcsr_create(m,m,nnzS); fasp_dcsr_cp(As,&mgl_p[0].A);
        mgl_p[0].b=fasp_dvec_create(m); mgl_p[0].x=fasp_dvec_create(m);
        // setup AMG
        switch (amgnsparam->param_p.AMG_type) {
            case CLASSIC_AMG:
                fasp_amg_setup_rs(mgl_p, &amgnsparam->param_p);
                break;
            case SA_AMG:
                fasp_amg_setup_sa(mgl_p, &amgnsparam->param_p);
                break;
            case UA_AMG:
                fasp_amg_setup_ua(mgl_p, &amgnsparam->param_p);
                break;
            default:
                printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                       amgnsparam->param_p.AMG_type);
                exit(ERROR_INPUT_PAR);
        }
    }
    else if (itparam->precond_type_p == 4) {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }
    
    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itparam->print_level_v;
    ITS_param_v.itsolver_type = itparam->itsolver_type_v;
    ITS_param_v.restart = itparam->pre_restart_v;
    ITS_param_v.tol = itparam->pre_tol_v;
    ITS_param_v.maxit = itparam->pre_maxit_v;
    ITS_param_v.precond_type = itparam->precond_type_v;
    
    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itparam->print_level_p;
    ITS_param_p.itsolver_type = itparam->itsolver_type_p;
    ITS_param_p.restart = itparam->pre_restart_p;
    ITS_param_p.tol = itparam->pre_tol_p;
    ITS_param_p.maxit = itparam->pre_maxit_p;
    ITS_param_p.precond_type = itparam->precond_type_p;
    
    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;
    
    precdata.colA = n;
    precdata.colB = m;
    precdata.col  = n+m;
    precdata.M    = M;
    precdata.B    = B;
    precdata.Bt   = Bt;
    precdata.C    = C;
    
    precdata.param_v        = &amgnsparam->param_v;
    precdata.param_p        = &amgnsparam->param_p;
    precdata.ITS_param_v    = &ITS_param_v;
    precdata.ITS_param_p    = &ITS_param_p;
    precdata.mgl_data_v     = mgl_v;
    precdata.mgl_data_p     = mgl_p;
    precdata.ILU_p          = &LU_p;
    
    precdata.max_levels     = mgl_v[0].num_levels;
    precdata.print_level    = amgnsparam->param_v.print_level;
    precdata.maxit          = amgnsparam->param_v.maxit;
    precdata.amg_tol        = amgnsparam->param_v.tol;
    precdata.cycle_type     = amgnsparam->param_v.cycle_type;
    precdata.smoother       = amgnsparam->param_v.smoother;
    precdata.presmooth_iter = amgnsparam->param_v.presmooth_iter;
    precdata.postsmooth_iter= amgnsparam->param_v.postsmooth_iter;
    precdata.relaxation     = amgnsparam->param_v.relaxation;
    precdata.coarse_scaling = amgnsparam->param_v.coarse_scaling;
    
    
    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;
    
    precdata.w = (REAL *)fasp_mem_calloc(precdata.col,sizeof(double));
    
    switch (precond_type) {
        case 1:
            prec.fct = fasp_precond_ns_bdiag;
            break;
        case 2:
            prec.fct = fasp_precond_ns_low_btri;
            break;
        case 3:
            prec.fct = fasp_precond_ns_up_btri;
            break;
        case 4:
            prec.fct = fasp_precond_ns_blu;
            break;
        default:
            printf("### ERROR: Unknown preconditioner type!\n");
            exit(ERROR_SOLVER_PRECTYPE);
    }
    
    setup_end = clock();
    
    if (PrtLvl>0) {
        setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
        printf("Setup costs %f.\n", setup_duration);
    }
    
    //------ solver phase ------//
    solver_start=clock();
    status=fasp_ns_solver_itsolver(Mat,b,x,&prec,itparam);
    solver_end=clock();
    
    if (PrtLvl>0) {
        solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
        printf(COLOR_RESET);
        printf("Solver costs %f seconds.\n", solver_duration);
        printf("Total costs %f seconds.\n", setup_duration + solver_duration);
    }
    
    //FINISHED:
    // clean up memory
    if (mgl_v) fasp_amg_data_free(mgl_v,&amgnsparam->param_v);
    if (itparam->precond_type_p == 1){fasp_dvec_free(&diag_S);}
    if (itparam->precond_type_p == 2) fasp_amg_data_free(mgl_p,&amgnsparam->param_p);
    
    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dcsr_free (&S);
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
/**
 * \fn static inline void get_schur_diagA (dCSRmat *B,dCSRmat *Bt,dCSRmat *A,dCSRmat *C,dCSRmat *S)
 * \brief Compute S = C+ B*diag(A)^{-1}*Bt
 *
 * \param *B:	       pointer to the Apu block
 * \param *Bt:	       pointer to the Aup block
 * \param *A:	       pointer to the Auu block
 * \param *C:          pointer to the App block
 * \param *S:          pointer to Schur complement
 *
 * \author Xiaozhe Hu and Lu Wang
 * \date 07/20/2014
 *
 * \note Assume B*diag(A)^{-1}*Bt has the right sign, i.e, it approximates -\Delta_p !!!
 */
static inline void get_schur_diagA (dCSRmat *B,
                                    dCSRmat *Bt,
                                    dCSRmat *A,
                                    dCSRmat *C,
                                    dCSRmat *S)
{
    INT colA = A->row;
    dvector diag_A;
    INT i;
    
    dCSRmat invA = fasp_dcsr_create(colA,colA,colA);
    
    fasp_dcsr_getdiag(A->row,A,&diag_A);
    
    // TODO: this part should be rewritten so that we do not form matrix invA explicitly!  -- Xiaozhe
    for (i=0;i<colA;i++)
    {
        invA.IA[i] = i;
        invA.JA[i] = i;
        if (diag_A.val[i] > SMALLREAL) invA.val[i]   = 1.0/diag_A.val[i];
        else invA.val[i] = 1.0;
    }
    invA.IA[colA] = colA;
    
    if (C) {
        dCSRmat tempA;
        fasp_blas_dcsr_rap(B, &invA, Bt, &tempA);
        fasp_blas_dcsr_add(C,1.0,&tempA,-1.0,S);
        fasp_dcsr_free(&tempA);
    }
    else {
        fasp_blas_dcsr_rap(B, &invA, Bt, S);
        fasp_blas_dcsr_axm(S, -1.0);
        
        /*
         // make it non-singular
         INT k,j,ibegin,iend;
         
         for (i=0;i<S->row;++i) {
         ibegin=S->IA[i]; iend=S->IA[i+1];
         for (k=ibegin;k<iend;++k) {
         j=S->JA[k];
         if ((j-i)==0) {
         S->val[k] = S->val[k] + 1e-8; break;
         } // end if
         } // end for k
         } // end for i
         */
    }
    
    fasp_dvec_free(&diag_A);
    fasp_dcsr_free(&invA);
}

/**
 * \fn static inline void get_schur_massP (dCSRmat *B,dCSRmat *Bt,dCSRmat *A, dCSRmat *Mp, REAL alpha, dCSRmat *S)
 * \brief Compute S = Mp + alpha * B*diag(A)^{-1}*Bt
 *
 * \param *B:	       pointer to the Apu block
 * \param *Bt:	       pointer to the Aup block
 * \param *A:	       pointer to the Auu block
 * \param *Mp:         pointer to the pressure mass
 * \param *alpha:      constant
 * \param *S:          pointer to Schur complement
 *
 * \author Xiaozhe Hu and Lu Wang
 * \date 07/20/2014
 *
 * \note Assume pressure block App=0 !!!
 * \note Assume B*diag(A)^{-1}*Bt has the right sign, i.e, it approximates -\Delta_p !!!
 * \note Still test the code !!!
 */
static inline void get_schur_massP (dCSRmat *B,
                                    dCSRmat *Bt,
                                    dCSRmat *A,
                                    dCSRmat *Mp,
                                    REAL     alpha,
                                    dCSRmat *S)
{
    INT colA = A->row;
    dvector diag_A;
    INT i;
    
    // get diagonal of A
    dCSRmat tempA;
    dCSRmat invA = fasp_dcsr_create(colA,colA,colA);
    fasp_dcsr_getdiag(A->row,A,&diag_A);
    
    // compute inverse of diagonal of A
    for (i=0;i<colA;i++)
    {
        invA.IA[i] = i;
        invA.JA[i] = i;
        if (diag_A.val[i] > SMALLREAL) invA.val[i]   = 1.0/diag_A.val[i];
        else invA.val[i] = 1.0;
    }
    invA.IA[colA] = colA;
    
    // compute B*diag(A)^{-1}*Bt
    fasp_blas_dcsr_rap(B, &invA, Bt, &tempA);
    
    // compute S = Mp + alpha * B*diag(A)^{-1}*Bt
    fasp_blas_dcsr_add(Mp, 1.0, &tempA, alpha, S);
    
    fasp_dvec_free(&diag_A);
    fasp_dcsr_free(&invA);
    fasp_dcsr_free(&tempA);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/