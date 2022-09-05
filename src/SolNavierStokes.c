/*! \file  SolNavierStokes.c
 *
 *  \brief Iterative solvers for Navier-Stokes matrices (main file)
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         PreNavierStokes.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2018 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  TODO: Clean up these functions! --Chensong
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

static inline void get_schur_diagA(dCSRmat *, dCSRmat *, dCSRmat *,
                                   dCSRmat *, dCSRmat *);
static inline void get_schur_pmass(dCSRmat *, dCSRmat *, dCSRmat *,
                                   dCSRmat *, REAL, dCSRmat *);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_ns_solver_itsolver (dCSRmat *A, dvector *b, dvector *x,
 *                                    precond *prec, itsolver_ns_param *itsparam)
 * \brief Solve Ax=b by standard Krylov methods
 *
 * \param A        pointer to the block dCSRmat matrix
 * \param b        pointer to the dvector of right hand side
 * \param x        pointer to the dvector of dofs
 * \param prec     pointer to the preconditioner data
 * \param itsparam pointer to parameters for iterative solvers
 *
 * \return         the number of iterations
 *
 * \author Lu Wang
 * \date   03/02/2012
 *
 * Modified by Xiaozhe Hu on 10/20/2013
 * Modified by Chensong Zhang on 03/18/2018
 */
SHORT fasp_ns_solver_itsolver(dBLCmat *A,
                              dvector *b,
                              dvector *x,
                              precond *prec,
                              itsolver_ns_param *itsparam)
{
    const SHORT PrtLvl = itsparam->print_level;
    const SHORT SolverType = itsparam->itsolver_type;
    const SHORT StopType = itsparam->stop_type;
    const SHORT MaxIt = itsparam->maxit;
    const SHORT restart = itsparam->restart;
    const REAL tol = itsparam->tol;
    const REAL abstol = itsparam->abstol;

    INT iter = ERROR_SOLVER_TYPE;
    REAL solver_start, solver_end, solver_duration;

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);
#endif

    fasp_gettime(&solver_start);

    switch (SolverType)
    {
    case SOLVER_BiCGstab:
        iter = fasp_solver_dblc_pbcgs(A, b, x, prec, tol,abstol, MaxIt, StopType, PrtLvl);
        break;
    case SOLVER_MinRes:
        iter = fasp_solver_dblc_pminres(A, b, x, prec, tol, abstol,MaxIt, StopType, PrtLvl);
        break;
    case SOLVER_VGMRES:
        iter = fasp_solver_dblc_pvgmres(A, b, x, prec, tol, abstol,MaxIt, restart, StopType, PrtLvl);
        break;
    case SOLVER_VFGMRES:
        iter = fasp_solver_dblc_pvfgmres(A, b, x, prec, tol,abstol, MaxIt, restart, StopType, PrtLvl);
        break;
    case SOLVER_GCR:
        iter = fasp_solver_dblc_pgcr(A, b, x, prec, tol, abstol,MaxIt, restart, StopType, PrtLvl);
        break;

    default:
        printf("### ERROR: Wrong itertive solver type %d!\n", SolverType);
        iter = ERROR_SOLVER_TYPE;
    }

    if ((PrtLvl > PRINT_MIN) && (iter >= 0))
    {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        fasp_cputime("Iterative method", solver_duration);
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return iter;
}

/**
 * \fn SHORT fasp_ns_solver_IR(dBLCmat *Mat, dvector *b, dvector *x,
 *                             precond *prec, itsolver_ns_param *itsparam)
 * \brief Solve Ax=b by IR methods for NS equations
 *
 * \param Mat       pointer to the dBLCmat matrix
 * \param b	        pointer to the dvector of right hand side
 * \param x	        pointer to the dvector of dofs
 * \param prec      pointer to parameters for  preconditioner
 * \param itsparam  pointer to parameters for iterative solvers
 *
 * \return          number of iterations
 *
 * \author  Ting Lai
 * \date    12/08/2022
 */
SHORT fasp_ns_solver_IR(dBLCmat *Mat,
                        dvector *b,
                        dvector *x,
                        precond *prec,
                        itsolver_ns_param *itsparam)
{
    INT i, j;
    INT col;
    SHORT iter = 0;
    SHORT num_IRiter = 0;
    SHORT status = FASP_SUCCESS;

    REAL rrn;
    REAL bnorm, rnorm;
    LONGREAL *ser_x;

    REAL *r_data;
    REAL *x_data;
    REAL *ser_b;

    const SHORT precond_type = itsparam->precond_type;
    const INT maxit_IR = itsparam->IRmaxit;
    const REAL IRtol = itsparam->tol;
    itsparam->tol = itsparam->IRtol;

    col = b->row;
    r_data = (REAL *)malloc(sizeof(REAL) * col);
    ser_x = (LONGREAL *)malloc(sizeof(LONGREAL) * col);
    
    // bnorm
    ser_b = b->val;
    bnorm = fasp_blas_darray_norm2(col, ser_b);
    rrn = 1.0;

    fasp_ldarray_set(col,ser_x,0.0); // TODO: x是输入的初值

    // IR Solve
    while (rrn > IRtol)
    {
        switch (precond_type)
        {
        case 0:
            status = fasp_ns_solver_itsolver(Mat, b, x, NULL, itsparam);
            break;
        default:
            status = fasp_ns_solver_itsolver(Mat, b, x, prec, itsparam);
        }

        // x=x+y;
        x_data = x->val;
        fasp_blas_ldarray_axpy(col,1,x_data,ser_x); // TODO: Add functions for LONGREAL

        // b-Ax long double
        fasp_darray_cp(col, ser_b, r_data);
        fasp_blas_ldblc_aAxpy(-1.0, Mat, ser_x, r_data);

        // rrn
        rnorm = fasp_blas_darray_norm2(col, r_data);
        rrn = rnorm / bnorm;
        
        
        if (status<0) {
            printf("\n### ERROR: IR Solver failed! Exit status = %d.\n\n", status);
            break;
        }else{
            iter += status; // TODO: What if status < 0 (error in fasp_ns_solver_itsolver)
            num_IRiter++;
        }    


        if (num_IRiter >= maxit_IR || rrn < IRtol)
        {
            printf("\nNumber of IR iterations = %d with IR tol relative residual %le . \n", iter, rrn);
            break;
        }
        // TODO: 算例901，IR没有作用

        fasp_darray_set(col,x_data,0.0);
        x->val = x_data;
        b->val = r_data;
    }

    for (i = 0; i < col; i++)
    {
        x_data[i] = ser_x[i];
    }
    x->val = x_data;
    return iter;
}

/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes (dBLCmat *Mat, dvector *b, dvector *x,
 *                                                  itsolver_ns_param *itsparam,
 *                                                  AMG_ns_param *amgparam,
 *                                                  ILU_param *iluparam,
 *                                                  SWZ_param *swzparam)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param Mat       pointer to the dBLCmat matrix
 * \param b         pointer to the dvector of right hand side
 * \param x         pointer to the dvector of dofs
 * \param itsparam   pointer to parameters for iterative solvers
 * \param amgparam  pionter to AMG parameters for N-S
 * \param iluparam  pionter to ILU parameters
 * \param swzparam  pionter to Schwarz parameters
 *
 * \return          number of iterations
 *
 * \author Lu Wang
 * \date   03/02/2012
 *
 * Modified by Xiaozhe Hu on 05/31/2016
 * Modified by Chensong Zhang on 03/15/2018
 */
SHORT fasp_solver_dblc_krylov_navier_stokes(dBLCmat *Mat,
                                            dvector *b,
                                            dvector *x,
                                            itsolver_ns_param *itsparam,
                                            AMG_ns_param *amgparam,
                                            ILU_param *iluparam,
                                            SWZ_param *swzparam)
{
    // parameters
    const SHORT PrtLvl = itsparam->print_level;
    const SHORT precond_type = itsparam->precond_type;
    const SHORT IR_type = itsparam->IR_type;
    const INT schwarz_mmsize = swzparam->SWZ_mmsize;
    const INT schwarz_maxlvl = swzparam->SWZ_maxlvl;
    const INT schwarz_type = swzparam->SWZ_type;

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    // Navier-Stokes 2 by 2 matrix
    dCSRmat *A = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B = Mat->blocks[2];
    dCSRmat *C = Mat->blocks[3];

    const INT n = A->row, m = B->row, nnzA = A->nnz;

    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S;
    SWZ_data schwarz_data;
    dvector diag_A;
    dvector diag_S;
    dCSRmat BABt;

    //------ setup phase ------//
    SHORT status = FASP_SUCCESS;
    REAL solver_start, solver_end, setup_start, setup_end;
    fasp_gettime(&setup_start);

#if DEBUG_MODE > 0
    printf("### DEBUG: n = %d, m = %d, nnz = %d\n", n, m, nnzA);
#endif

    //-------------------------//
    // setup AMG for velocity  //
    //-------------------------//
    AMG_data *mgl_v = fasp_amg_data_create(amgparam->param_v.max_levels);
    mgl_v[0].A = fasp_dcsr_create(n, n, nnzA);

    if (precond_type > 10)
    {

        dCSRmat BtB;
        fasp_blas_dcsr_mxm(Bt, B, &BtB);

        REAL gamma = 10;
        fasp_blas_dcsr_add(A, 1.0, &BtB, gamma, &mgl_v[0].A);

        fasp_dcsr_free(&BtB);
    }
    else
    {
        fasp_dcsr_cp(A, &mgl_v[0].A);
    }

    mgl_v[0].b = fasp_dvec_create(n);
    mgl_v[0].x = fasp_dvec_create(n);

    // setup AMG
    switch (amgparam->param_v.AMG_type)
    {
    case CLASSIC_AMG:
        fasp_amg_setup_rs(mgl_v, &amgparam->param_v);
        break;
    case SA_AMG:
        fasp_amg_setup_sa(mgl_v, &amgparam->param_v);
        break;
    case UA_AMG:
        fasp_amg_setup_ua(mgl_v, &amgparam->param_v);
        break;
    default:
        printf("### ERROR: Wrong AMG type %d!\n", amgparam->param_v.AMG_type);
        exit(ERROR_INPUT_PAR);
    }

    // get diagonal of A
    fasp_dcsr_getdiag(n, &mgl_v[0].A, &diag_A);

    //---------------------------//
    // setup Schur complement S  //
    //---------------------------//

    if (precond_type == 8 || precond_type == 9 ||
        precond_type == 18 || precond_type == 19)
    {

        fasp_blas_dcsr_mxm(B, Bt, &S);

        // change the sign of the BB^T
        fasp_blas_dcsr_axm(&S, -1.0);

        // make it non-singular
        INT i, k, j, ibegin, iend;

        for (i = 0; i < S.row; ++i)
        {
            ibegin = S.IA[i];
            iend = S.IA[i + 1];
            for (k = ibegin; k < iend; ++k)
            {
                j = S.JA[k];
                if (j == i)
                {
                    S.val[k] = S.val[k] + 1e-8;
                    break;
                } // end if
            }     // end for k
        }         // end for i
    }
    else if (precond_type == 10 || precond_type == 20)
    {

        fasp_blas_dcsr_mxm(B, Bt, &S);
        fasp_blas_dcsr_rap(B, A, Bt, &BABt);

        // change the sign of the BB^T
        fasp_blas_dcsr_axm(&S, -1.0);

        // make it non-singular
        INT i, k, j, ibegin, iend;

        for (i = 0; i < S.row; ++i)
        {
            ibegin = S.IA[i];
            iend = S.IA[i + 1];
            for (k = ibegin; k < iend; ++k)
            {
                j = S.JA[k];
                if (j == i)
                {
                    S.val[k] = S.val[k] + 1e-8;
                    break;
                } // end if
            }     // end for k
        }         // end for i
    }
    else
    {
        get_schur_diagA(B, Bt, A, C, &S);
    }

    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);

    AMG_data *mgl_p;
    ILU_data LU_p;

    if (itsparam->precond_type_p == 1)
    {
        fasp_dcsr_getdiag(0, &S, &diag_S);
    }

    else if (itsparam->precond_type_p == 2)
    {
        // Setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p = fasp_amg_data_create(amgparam->param_p.max_levels);
        mgl_p[0].A = fasp_dcsr_create(m, m, nnzS);
        fasp_dcsr_cp(As, &mgl_p[0].A);
        mgl_p[0].b = fasp_dvec_create(m);
        mgl_p[0].x = fasp_dvec_create(m);

        switch (amgparam->param_p.AMG_type)
        {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_p, &amgparam->param_p);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_p, &amgparam->param_p);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_p, &amgparam->param_p);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                   amgparam->param_p.AMG_type);
            exit(ERROR_INPUT_PAR);
        }
    }

    else if (itsparam->precond_type_p == 4)
    {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }

    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itsparam->print_level_v;
    ITS_param_v.itsolver_type = itsparam->itsolver_type_v;
    ITS_param_v.restart = itsparam->pre_restart_v;
    ITS_param_v.tol = itsparam->pre_tol_v;
    ITS_param_v.abstol = itsparam->pre_abstol_v;
    ITS_param_v.maxit = itsparam->pre_maxit_v;
    ITS_param_v.precond_type = itsparam->precond_type_v;

    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itsparam->print_level_p;
    ITS_param_p.itsolver_type = itsparam->itsolver_type_p;
    ITS_param_p.restart = itsparam->pre_restart_p;
    ITS_param_p.tol = itsparam->pre_tol_p;
    ITS_param_p.abstol = itsparam->pre_abstol_p;
    ITS_param_p.maxit = itsparam->pre_maxit_p;
    ITS_param_p.precond_type = itsparam->precond_type_p;

    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;

    precdata.colA = n;
    precdata.colB = m;
    precdata.col = n + m;
    precdata.M = M;
    precdata.B = B;
    precdata.Bt = Bt;
    precdata.C = C;
    precdata.BABt = &BABt;

    precdata.param_v = &amgparam->param_v;
    precdata.param_p = &amgparam->param_p;
    precdata.ITS_param_v = &ITS_param_v;
    precdata.ITS_param_p = &ITS_param_p;
    precdata.mgl_data_v = mgl_v;
    precdata.mgl_data_p = mgl_p;
    precdata.ILU_p = &LU_p;

    precdata.max_levels = mgl_v[0].num_levels;
    precdata.print_level = amgparam->param_v.print_level;
    precdata.maxit = amgparam->param_v.maxit;
    precdata.amg_tol = amgparam->param_v.tol;
    precdata.cycle_type = amgparam->param_v.cycle_type;
    precdata.smoother = amgparam->param_v.smoother;
    precdata.presmooth_iter = amgparam->param_v.presmooth_iter;
    precdata.postsmooth_iter = amgparam->param_v.postsmooth_iter;
    precdata.relaxation = amgparam->param_v.relaxation;
    precdata.coarse_scaling = amgparam->param_v.coarse_scaling;

    precdata.diag_A = &diag_A;
    precdata.diag_S = &diag_S;
    precdata.S = &S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;
    precdata.w = (REAL *)fasp_mem_calloc(precdata.col, sizeof(double));

    switch (precond_type)
    {
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
        printf("### WARNING: Unknown preconditioner type!\n");
    }

    setup_end = clock();

    if (PrtLvl > 0)
    {
        fasp_gettime(&setup_end);
        fasp_cputime("NS Setup", setup_end - setup_start);
    }

    //------ solve phase ------//
    fasp_gettime(&solver_start);

    if (IR_type > 0)
    {
        switch (precond_type)
        {
        case 0:
            status = fasp_ns_solver_IR(Mat, b, x, NULL, itsparam);
            break;
        default:
            status = fasp_ns_solver_IR(Mat, b, x, &prec, itsparam);
        }
    }
    else
    {
        switch (precond_type)
        {
        case 0:
            status = fasp_ns_solver_itsolver(Mat, b, x, NULL, itsparam);
            break;
        default:
            status = fasp_ns_solver_itsolver(Mat, b, x, &prec, itsparam);
        }
    }

    if (PrtLvl > 0)
    {
        fasp_gettime(&solver_end);
        printf(COLOR_RESET);
        fasp_cputime("NS Solve", solver_end - solver_start);
        fasp_cputime("NS Total", solver_end - setup_start);
    }

    // clean up memory
    if (mgl_v)
        fasp_amg_data_free(mgl_v, &amgparam->param_v);
    if (itsparam->precond_type_p == 1)
        fasp_dvec_free(&diag_S);
    if (itsparam->precond_type_p == 2)
        fasp_amg_data_free(mgl_p, &amgparam->param_p);

    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dvec_free(&diag_A);
    fasp_dcsr_free(&S);
    if (precond_type == 10 || precond_type == 20)
        fasp_dcsr_free(&BABt);

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes_pmass (dBLCmat *Mat, dvector *b,
 *                                                        dvector *x,
 *                                                        itsolver_ns_param *itsparam,
 *                                                        AMG_ns_param *amgparam,
 *                                                        ILU_param *iluparam,
 *                                                        SWZ_param *swzparam,
 *                                                        dCSRmat *Mp)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param Mat       pointer to the dBLCmat matrix
 * \param b	        pointer to the dvector of right hand side
 * \param x	        pointer to the dvector of dofs
 * \param itsparam  pointer to parameters for iterative solvers
 * \param amgparam  AMG parameters for NS
 * \param iluparam  ILU parameters
 * \param swzparam  Schwarz parameters
 * \param Mp        pointer to dCSRmat of the pressure mass matrix
 *
 * \return          number of iterations
 *
 * \author Xiaozhe Hu
 * \date   017/07/2014
 *
 * \note In general, this is for purely Stokes problem, NS problem with div-div
 *       stablization -- Xiaozhe
 */
SHORT fasp_solver_dblc_krylov_navier_stokes_pmass(dBLCmat *Mat,
                                                  dvector *b,
                                                  dvector *x,
                                                  itsolver_ns_param *itsparam,
                                                  AMG_ns_param *amgparam,
                                                  ILU_param *iluparam,
                                                  SWZ_param *swzparam,
                                                  dCSRmat *Mp)
{
    // parameters
    const SHORT PrtLvl = itsparam->print_level;
    const SHORT precond_type = itsparam->precond_type;
    const INT schwarz_mmsize = swzparam->SWZ_mmsize;
    const INT schwarz_maxlvl = swzparam->SWZ_maxlvl;
    const INT schwarz_type = swzparam->SWZ_type;

    // Navier-Stokes 4 by 4 matrix
    dCSRmat *A = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B = Mat->blocks[2];
    dCSRmat *C = Mat->blocks[3];
    const INT n = A->row, m = B->row, nnzA = A->nnz;

    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S, P;
    SWZ_data schwarz_data;
    dvector diag_S;

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    //------ setup phase ------//
    SHORT status = FASP_SUCCESS;
    REAL solver_start, solver_end, setup_start, setup_end;
    fasp_gettime(&setup_start);

    //-----------------------//
    // setup AMG for velocity
    //-----------------------//
    AMG_data *mgl_v = fasp_amg_data_create(amgparam->param_v.max_levels);
    mgl_v[0].A = fasp_dcsr_create(n, n, nnzA);
    fasp_dcsr_cp(A, &mgl_v[0].A);
    mgl_v[0].b = fasp_dvec_create(n);
    mgl_v[0].x = fasp_dvec_create(n);

    // setup AMG
    switch (amgparam->param_v.AMG_type)
    {
    case CLASSIC_AMG:
        fasp_amg_setup_rs(mgl_v, &amgparam->param_v);
        break;
    case SA_AMG:
        fasp_amg_setup_sa(mgl_v, &amgparam->param_v);
        break;
    case UA_AMG:
        fasp_amg_setup_ua(mgl_v, &amgparam->param_v);
        break;
    default:
        printf("### ERROR: Wrong AMG type %d!\n", amgparam->param_v.AMG_type);
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

    if (itsparam->precond_type_p == 1)
    {
        fasp_dcsr_getdiag(0, &S, &diag_S);
    }
    else if (itsparam->precond_type_p == 2)
    {
        // setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p = fasp_amg_data_create(amgparam->param_p.max_levels);
        mgl_p[0].A = fasp_dcsr_create(m, m, nnzS);
        fasp_dcsr_cp(As, &mgl_p[0].A);
        mgl_p[0].b = fasp_dvec_create(m);
        mgl_p[0].x = fasp_dvec_create(m);

        switch (amgparam->param_p.AMG_type)
        {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_p, &amgparam->param_p);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_p, &amgparam->param_p);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_p, &amgparam->param_p);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                   amgparam->param_p.AMG_type);
            exit(ERROR_INPUT_PAR);
        }
    }
    else if (itsparam->precond_type_p == 4)
    {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }

    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itsparam->print_level_v;
    ITS_param_v.itsolver_type = itsparam->itsolver_type_v;
    ITS_param_v.restart = itsparam->pre_restart_v;
    ITS_param_v.tol = itsparam->pre_tol_v;
    ITS_param_v.maxit = itsparam->pre_maxit_v;
    ITS_param_v.precond_type = itsparam->precond_type_v;

    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itsparam->print_level_p;
    ITS_param_p.itsolver_type = itsparam->itsolver_type_p;
    ITS_param_p.restart = itsparam->pre_restart_p;
    ITS_param_p.tol = itsparam->pre_tol_p;
    ITS_param_p.maxit = itsparam->pre_maxit_p;
    ITS_param_p.precond_type = itsparam->precond_type_p;

    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;

    precdata.colA = n;
    precdata.colB = m;
    precdata.col = n + m;
    precdata.M = M;
    precdata.B = B;
    precdata.Bt = Bt;
    precdata.C = C;

    precdata.param_v = &amgparam->param_v;
    precdata.param_p = &amgparam->param_p;
    precdata.ITS_param_v = &ITS_param_v;
    precdata.ITS_param_p = &ITS_param_p;
    precdata.mgl_data_v = mgl_v;
    precdata.mgl_data_p = mgl_p;
    precdata.ILU_p = &LU_p;

    precdata.max_levels = mgl_v[0].num_levels;
    precdata.print_level = amgparam->param_v.print_level;
    precdata.maxit = amgparam->param_v.maxit;
    precdata.amg_tol = amgparam->param_v.tol;
    precdata.cycle_type = amgparam->param_v.cycle_type;
    precdata.smoother = amgparam->param_v.smoother;
    precdata.presmooth_iter = amgparam->param_v.presmooth_iter;
    precdata.postsmooth_iter = amgparam->param_v.postsmooth_iter;
    precdata.relaxation = amgparam->param_v.relaxation;
    precdata.coarse_scaling = amgparam->param_v.coarse_scaling;

    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;

    precdata.w = (REAL *)fasp_mem_calloc(precdata.col, sizeof(double));

    switch (precond_type)
    {
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

    if (PrtLvl > 0)
    {
        fasp_gettime(&setup_end);
        fasp_cputime("NS Setup", setup_end - setup_start);
    }

    //------ solver phase ------//
    fasp_gettime(&solver_start);
    status = fasp_ns_solver_itsolver(Mat, b, x, &prec, itsparam);

    if (PrtLvl > 0)
    {
        fasp_gettime(&solver_end);
        printf(COLOR_RESET);
        fasp_cputime("NS Solve", solver_end - solver_start);
        fasp_cputime("NS Total", solver_end - setup_start);
    }

    // clean up memory
    if (mgl_v)
        fasp_amg_data_free(mgl_v, &amgparam->param_v);
    if (itsparam->precond_type_p == 1)
        fasp_dvec_free(&diag_S);
    if (itsparam->precond_type_p == 2)
        fasp_amg_data_free(mgl_p, &amgparam->param_p);

    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dcsr_free(&S);

#if DEBUG_MODE > 0
    printf("### DEBUG: [-End-] %s ...\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn SHORT fasp_solver_dblc_krylov_navier_stokes_schur_pmass (dBLCmat *Mat, dvector *b,
 *                                                              dvector *x,
 *                                                              itsolver_ns_param *itsparam,
 *                                                              AMG_ns_param *amgparam,
 *                                                              ILU_param *iluparam,
 *                                                              SWZ_param *swzparam,
 *                                                              dCSRmat *Mp)
 * \brief Solve Ax=b by standard Krylov methods for NS equations
 *
 * \param Mat       pointer to the dBLCmat matrix
 * \param b	        pointer to the dvector of right hand side
 * \param x	        pointer to the dvector of dofs
 * \param itsparam  pointer to parameters for iterative solvers
 * \param amgparam  AMG parameters for NS
 * \param iluparam  ILU parameters
 * \param swzparam  Schwarz parameters
 * \param Mp        pointer to dCSRmat of the pressure mass matrix
 *
 * \return          number of iterations
 *
 * \author Xiaozhe Hu
 * \date   017/07/2014
 *
 * \note In general, this is for NS problems without div-div stablization and
 *       pressure stablization (pressure block is zero), moreover, pressure mass
 *       matrix is provided.
 */
SHORT fasp_solver_dblc_krylov_navier_stokes_schur_pmass(dBLCmat *Mat,
                                                        dvector *b,
                                                        dvector *x,
                                                        itsolver_ns_param *itsparam,
                                                        AMG_ns_param *amgparam,
                                                        ILU_param *iluparam,
                                                        SWZ_param *swzparam,
                                                        dCSRmat *Mp)
{
    // parameters
    const SHORT PrtLvl = itsparam->print_level;
    const SHORT precond_type = itsparam->precond_type;
    const INT schwarz_mmsize = swzparam->SWZ_mmsize;
    const INT schwarz_maxlvl = swzparam->SWZ_maxlvl;
    const INT schwarz_type = swzparam->SWZ_type;

    // Navier-Stokes 4 by 4 matrix
    dCSRmat *A = Mat->blocks[0];
    dCSRmat *Bt = Mat->blocks[1];
    dCSRmat *B = Mat->blocks[2];
    dCSRmat *C = Mat->blocks[3];
    const INT n = A->row, m = B->row, nnzA = A->nnz;

    // preconditioner data
    dCSRmat *M = Mat->blocks[3];
    dCSRmat S, P;
    SWZ_data schwarz_data;
    dvector diag_S;

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    //------ setup phase ------//
    SHORT status = FASP_SUCCESS;
    REAL solver_start, solver_end, setup_start, setup_end;
    fasp_gettime(&setup_start);

    //-----------------------//
    // setup AMG for velocity
    //-----------------------//

    AMG_data *mgl_v = fasp_amg_data_create(amgparam->param_v.max_levels);
    mgl_v[0].A = fasp_dcsr_create(n, n, nnzA);
    fasp_dcsr_cp(A, &mgl_v[0].A);
    mgl_v[0].b = fasp_dvec_create(n);
    mgl_v[0].x = fasp_dvec_create(n);

    // setup AMG
    switch (amgparam->param_v.AMG_type)
    {
    case CLASSIC_AMG:
        fasp_amg_setup_rs(mgl_v, &amgparam->param_v);
        break;
    case SA_AMG:
        fasp_amg_setup_sa(mgl_v, &amgparam->param_v);
        break;
    case UA_AMG:
        fasp_amg_setup_ua(mgl_v, &amgparam->param_v);
        break;
    default:
        printf("### ERROR: Wrong AMG type %d!\n", amgparam->param_v.AMG_type);
        exit(ERROR_INPUT_PAR);
    }

    //-------------------------//
    // setup Schur complement S using pressure mass
    //-------------------------//

    get_schur_pmass(B, Bt, &mgl_v[0].A, Mp, 1e5, &S);
    // TODO: 1e5 is a parameter can be tuned, not sure how to tune now -- Xiaozhe

    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);

    AMG_data *mgl_p;
    ILU_data LU_p;

    if (itsparam->precond_type_p == 1)
    {
        fasp_dcsr_getdiag(0, &S, &diag_S);
    }

    else if (itsparam->precond_type_p == 2)
    {
        // Setup AMG for Schur Complement
        dCSRmat *As = &S;
        const INT nnzS = As->nnz;
        mgl_p = fasp_amg_data_create(amgparam->param_p.max_levels);
        mgl_p[0].A = fasp_dcsr_create(m, m, nnzS);
        fasp_dcsr_cp(As, &mgl_p[0].A);
        mgl_p[0].b = fasp_dvec_create(m);
        mgl_p[0].x = fasp_dvec_create(m);

        switch (amgparam->param_p.AMG_type)
        {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl_p, &amgparam->param_p);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl_p, &amgparam->param_p);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl_p, &amgparam->param_p);
            break;
        default:
            printf("### ERROR: Wrong AMG type %d for Schur Complement!\n",
                   amgparam->param_p.AMG_type);
            exit(ERROR_INPUT_PAR);
        }
    }

    else if (itsparam->precond_type_p == 4)
    {
        // setup ILU for Schur Complement
        fasp_ilu_dcsr_setup(&S, &LU_p, iluparam);
        fasp_mem_iludata_check(&LU_p);
    }

    //---------------------------------------//
    // Setup itsolver parameter for subblocks
    //---------------------------------------//
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itsparam->print_level_v;
    ITS_param_v.itsolver_type = itsparam->itsolver_type_v;
    ITS_param_v.restart = itsparam->pre_restart_v;
    ITS_param_v.tol = itsparam->pre_tol_v;
    ITS_param_v.maxit = itsparam->pre_maxit_v;
    ITS_param_v.precond_type = itsparam->precond_type_v;

    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itsparam->print_level_p;
    ITS_param_p.itsolver_type = itsparam->itsolver_type_p;
    ITS_param_p.restart = itsparam->pre_restart_p;
    ITS_param_p.tol = itsparam->pre_tol_p;
    ITS_param_p.maxit = itsparam->pre_maxit_p;
    ITS_param_p.precond_type = itsparam->precond_type_p;

    //-------------------------//
    // setup preconditioner
    //-------------------------//
    precond prec;
    precond_ns_data precdata;
    prec.data = &precdata;

    precdata.colA = n;
    precdata.colB = m;
    precdata.col = n + m;
    precdata.M = M;
    precdata.B = B;
    precdata.Bt = Bt;
    precdata.C = C;

    precdata.param_v = &amgparam->param_v;
    precdata.param_p = &amgparam->param_p;
    precdata.ITS_param_v = &ITS_param_v;
    precdata.ITS_param_p = &ITS_param_p;
    precdata.mgl_data_v = mgl_v;
    precdata.mgl_data_p = mgl_p;
    precdata.ILU_p = &LU_p;

    precdata.max_levels = mgl_v[0].num_levels;
    precdata.print_level = amgparam->param_v.print_level;
    precdata.maxit = amgparam->param_v.maxit;
    precdata.amg_tol = amgparam->param_v.tol;
    precdata.cycle_type = amgparam->param_v.cycle_type;
    precdata.smoother = amgparam->param_v.smoother;
    precdata.presmooth_iter = amgparam->param_v.presmooth_iter;
    precdata.postsmooth_iter = amgparam->param_v.postsmooth_iter;
    precdata.relaxation = amgparam->param_v.relaxation;
    precdata.coarse_scaling = amgparam->param_v.coarse_scaling;

    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;

    precdata.w = (REAL *)fasp_mem_calloc(precdata.col, sizeof(double));

    switch (precond_type)
    {
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

    if (PrtLvl > 0)
    {
        fasp_gettime(&setup_end);
        fasp_cputime("NS Setup", setup_end - setup_start);
    }

    //------ solver phase ------//
    fasp_gettime(&solver_start);

    status = fasp_ns_solver_itsolver(Mat, b, x, &prec, itsparam);

    if (PrtLvl > 0)
    {
        fasp_gettime(&solver_end);
        printf(COLOR_RESET);
        fasp_cputime("NS Solve", solver_end - solver_start);
        fasp_cputime("NS Total", solver_end - setup_start);
    }

    // clean up memory
    if (mgl_v)
        fasp_amg_data_free(mgl_v, &amgparam->param_v);
    if (itsparam->precond_type_p == 1)
        fasp_dvec_free(&diag_S);
    if (itsparam->precond_type_p == 2)
        fasp_amg_data_free(mgl_p, &amgparam->param_p);

    fasp_mem_free(precdata.w);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
    fasp_dcsr_free(&S);

#if DEBUG_MODE > 0
    printf("### DEBUG: [-End-] %s ...\n", __FUNCTION__);
#endif

    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
/**
 * \fn static inline void get_schur_diagA (dCSRmat *B,dCSRmat *Bt,dCSRmat *A,
 *                                         dCSRmat *C,dCSRmat *S)
 * \brief Compute S = C+ B*diag(A)^{-1}*Bt
 *
 * \param B   pointer to the Apu block
 * \param Bt  pointer to the Aup block
 * \param A   pointer to the Auu block
 * \param C   pointer to the App block
 * \param S   pointer to Schur complement
 *
 * \author Xiaozhe Hu and Lu Wang
 * \date 07/20/2014
 *
 * \note Assume B*diag(A)^{-1}*Bt has the right sign, i.e, it approximates -\Delta_p !!!
 */
static inline void get_schur_diagA(dCSRmat *B,
                                   dCSRmat *Bt,
                                   dCSRmat *A,
                                   dCSRmat *C,
                                   dCSRmat *S)
{
    INT colA = A->row;
    dvector diag_A;
    INT i;

    dCSRmat invA = fasp_dcsr_create(colA, colA, colA);

    fasp_dcsr_getdiag(A->row, A, &diag_A);

    // TODO: this part should be rewritten so that we do not form matrix invA explicitly!  -- Xiaozhe
    for (i = 0; i < colA; i++)
    {
        invA.IA[i] = i;
        invA.JA[i] = i;
        if (diag_A.val[i] > SMALLREAL)
            invA.val[i] = 1.0 / diag_A.val[i];
        else
            invA.val[i] = 1.0;
    }
    invA.IA[colA] = colA;

    if (C)
    {
        dCSRmat tempA;
        fasp_blas_dcsr_rap(B, &invA, Bt, &tempA);
        fasp_blas_dcsr_add(C, 1.0, &tempA, -1.0, S);
        fasp_dcsr_free(&tempA);
    }
    else
    {
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
 * \fn static inline void get_schur_pmass (dCSRmat *B,dCSRmat *Bt,dCSRmat *A,
 *                                         dCSRmat *Mp, REAL alpha, dCSRmat *S)
 * \brief Compute S = Mp + alpha * B*diag(A)^{-1}*Bt
 *
 * \param B      pointer to the Apu block
 * \param Bt     pointer to the Aup block
 * \param A      pointer to the Auu block
 * \param Mp     pointer to the pressure mass
 * \param alpha  constant
 * \param S      pointer to Schur complement
 *
 * \author Xiaozhe Hu and Lu Wang
 * \date   07/20/2014
 *
 * \note Assume pressure block App=0 !!!
 * \note Assume B*diag(A)^{-1}*Bt has the right sign \approx -\Delta_p !!!
 *
 * TODO: Still need to test the code !!! --Chensong
 */
static inline void get_schur_pmass(dCSRmat *B,
                                   dCSRmat *Bt,
                                   dCSRmat *A,
                                   dCSRmat *Mp,
                                   REAL alpha,
                                   dCSRmat *S)
{
    INT colA = A->row;
    dvector diag_A;
    INT i;

    // get diagonal of A
    dCSRmat tempA;
    dCSRmat invA = fasp_dcsr_create(colA, colA, colA);
    fasp_dcsr_getdiag(A->row, A, &diag_A);

    // compute inverse of diagonal of A
    for (i = 0; i < colA; i++)
    {
        invA.IA[i] = i;
        invA.JA[i] = i;
        if (diag_A.val[i] > SMALLREAL)
            invA.val[i] = 1.0 / diag_A.val[i];
        else
            invA.val[i] = 1.0;
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
