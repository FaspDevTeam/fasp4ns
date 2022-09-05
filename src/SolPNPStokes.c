/*! \file  SolPNPStokes.c
 *
 *  \brief Iterative solvers for PNP-Stokes system (main file)
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         PreNavierStokes.c and PrePNPStokes.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--Present by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dblc_krylov_pnp_stokes (dBLCmat *A, dvector *b, dvector *x,
 *                                             ITS_param *itparam, 
 *                                             ITS_param *itparam_pnp,
 *                                             AMG_param *amgparam_pnp,
 *                                             itsolver_ns_param *itparam_stokes,
 *                                             AMG_ns_param *amgparam_stokes,
 *                                             const int num_velocity,
 *                                             const int num_pressure)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A               Pointer to the coeff matrix in dBLCmat format
 * \param b               Pointer to the right hand side in dvector format
 * \param x               Pointer to the approx solution in dvector format
 * \param itparam         Pointer to parameters for iterative solvers
 * \param itparam_pnp     Pointer to parameters for iterative solvers for PNP
 * \param amgparam_pnp    Pointer to parameters for AMG solvers for PNP
 * \param itparam_stokes  Pointer to parameters for iterative solvers for Stokes
 * \param amgparam_stokes Pointer to parameters for AMG solvers for Stokes
 * \param num_velocity    Size of velocity vector
 * \param num_pressure    Size of pressure vector
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/12/2016
 */
INT fasp_solver_dblc_krylov_pnp_stokes (dBLCmat *A,
                                        dvector *b,
                                        dvector *x,
                                        ITS_param *itparam,
                                        ITS_param *itparam_pnp,
                                        AMG_param *amgparam_pnp,
                                        itsolver_ns_param *itparam_stokes,
                                        AMG_ns_param *amgparam_stokes,
                                        const int num_velocity,
                                        const int num_pressure)
{
    const SHORT prtlvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
    
    INT m, n, nnz, i, k;
    
    // local variables
    dCSRmat A_pnp_csr;
    dBSRmat A_pnp_bsr;
    
    dCSRmat A_stokes_csr;
    dBLCmat A_stokes_bcsr;
    dCSRmat S;
    dCSRmat BABt;
    
    // data for pnp
    AMG_data_bsr *mgl_pnp=fasp_amg_data_bsr_create(amgparam_pnp->max_levels);
    ILU_param iluparam_pnp;
    iluparam_pnp.print_level = amgparam_pnp->print_level;
    iluparam_pnp.ILU_lfil    = amgparam_pnp->ILU_lfil;
    iluparam_pnp.ILU_droptol = amgparam_pnp->ILU_droptol;
    iluparam_pnp.ILU_relax   = amgparam_pnp->ILU_relax;
    iluparam_pnp.ILU_type    = amgparam_pnp->ILU_type;
    ILU_data ILU_pnp;
    
    // data for stokes
    AMG_data *mgl_v = fasp_amg_data_create(amgparam_stokes->param_v.max_levels);
    AMG_data *mgl_p;
    dvector   res_p = fasp_dvec_create(num_pressure);
    dvector   sol_p = fasp_dvec_create(num_pressure);
    
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
        
        printf("Factorization for stokes diagonal block: \n");
        LU_diag[1] = fasp_umfpack_factorize(&A_stokes_csr, prtlvl);
#endif
    }
    
    /* diagonal blocks are solved inexactly */
    else if ( precond_type > 30 && precond_type < 40 ) {
        
        // pnp block
        {
            A_pnp_bsr = fasp_format_dcsr_dbsr(A->blocks[0], 3);
            
            // AMG for pnp
            /*
             // initialize A, b, x for mgl_pnp[0]
             mgl_pnp[0].A = fasp_dbsr_create(A_pnp_bsr.ROW, A_pnp_bsr.COL, A_pnp_bsr.NNZ, A_pnp_bsr.nb, A_pnp_bsr.storage_manner);
             mgl_pnp[0].b = fasp_dvec_create(mgl_pnp[0].A.ROW*mgl_pnp[0].A.nb);
             mgl_pnp[0].x = fasp_dvec_create(mgl_pnp[0].A.COL*mgl_pnp[0].A.nb);
             
             fasp_dbsr_cp(&A_pnp_bsr, &(mgl_pnp[0].A));
             
             switch (amgparam_pnp->AMG_type) {
             
             case SA_AMG: // Smoothed Aggregation AMG
             status = fasp_amg_setup_sa_bsr(mgl_pnp, amgparam_pnp); break;
             
             default:
             status = fasp_amg_setup_ua_bsr(mgl_pnp, amgparam_pnp); break;
             
             }
             
             if (status < 0) goto FINISHED;
             */
            
            // diagonal preconditioner for pnp
            /*
             // diag of the pnp matrix
             fasp_dvec_alloc(A_pnp_bsr.ROW*A_pnp_bsr.nb*A_pnp_bsr.nb, &diag_pnp);
             for (i = 0; i < A_pnp_bsr.ROW; ++i) {
             for (k = A_pnp_bsr.IA[i]; k < A_pnp_bsr.IA[i+1]; ++k) {
             if (A_pnp_bsr.JA[k] == i)
             memcpy(diag_pnp.val+i*A_pnp_bsr.nb*A_pnp_bsr.nb, A_pnp_bsr.val+k*A_pnp_bsr.nb*A_pnp_bsr.nb, A_pnp_bsr.nb*A_pnp_bsr.nb*sizeof(REAL));
             }
             }
             
             for (i=0; i<A_pnp_bsr.ROW; ++i){
             fasp_blas_smat_inv(&(diag_pnp.val[i*A_pnp_bsr.nb*A_pnp_bsr.nb]), A_pnp_bsr.nb);
             }
             */
            
            // BSR ILU for pnp
            /*
             // ILU setup
             if ( (status = fasp_ilu_dbsr_setup(&A_pnp_bsr, &ILU_pnp, &iluparam_pnp)) < 0 ) goto FINISHED;
             
             // check iludata
             if ( (status = fasp_mem_iludata_check(&ILU_pnp)) < 0 ) goto FINISHED;
             */
            
            // CSR ILU for pnp
            /*
             // ILU setup for whole matrix
             if ( (status = fasp_ilu_dcsr_setup(A->blocks[0],&ILU_pnp,&iluparam_pnp)) < 0 ) goto FINISHED;
             
             // check iludata
             if ( (status = fasp_mem_iludata_check(&ILU_pnp)) < 0 ) goto FINISHED;
             */
            
        }
        
        // stokes block
        {
            A_stokes_bcsr.brow = 2;
            A_stokes_bcsr.bcol = 2;
            A_stokes_bcsr.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
            for (i=0; i<4 ;i++) {
                A_stokes_bcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            }
            
            ivector velocity_idx;
            ivector pressure_idx;
            fasp_ivec_alloc(num_velocity, &velocity_idx);
            fasp_ivec_alloc(num_pressure, &pressure_idx);
            for (i=0; i<num_velocity; i++) velocity_idx.val[i] = i;
            for (i=0; i<num_pressure; i++) pressure_idx.val[i] = num_velocity + i;
            
            fasp_dcsr_getblk(A->blocks[3], velocity_idx.val, velocity_idx.val, velocity_idx.row, velocity_idx.row, A_stokes_bcsr.blocks[0]);
            fasp_dcsr_getblk(A->blocks[3], velocity_idx.val, pressure_idx.val, velocity_idx.row, pressure_idx.row, A_stokes_bcsr.blocks[1]);
            fasp_dcsr_getblk(A->blocks[3], pressure_idx.val, velocity_idx.val, pressure_idx.row, velocity_idx.row, A_stokes_bcsr.blocks[2]);
            fasp_dcsr_getblk(A->blocks[3], pressure_idx.val, pressure_idx.val, pressure_idx.row, pressure_idx.row, A_stokes_bcsr.blocks[3]);
            
            fasp_ivec_free(&velocity_idx);
            fasp_ivec_free(&pressure_idx);
            
            
            // AMG for velocity
            mgl_v[0].A=fasp_dcsr_create(A_stokes_bcsr.blocks[0]->row,A_stokes_bcsr.blocks[0]->col,A_stokes_bcsr.blocks[0]->nnz);
            fasp_dcsr_cp(A_stokes_bcsr.blocks[0], &mgl_v[0].A);
            mgl_v[0].b=fasp_dvec_create(A_stokes_bcsr.blocks[0]->row); mgl_v[0].x=fasp_dvec_create(A_stokes_bcsr.blocks[0]->col);
            
            switch (amgparam_stokes->param_v.AMG_type) {
                case CLASSIC_AMG:
                    fasp_amg_setup_rs(mgl_v, &amgparam_stokes->param_v);
                    break;
                case SA_AMG:
                    fasp_amg_setup_sa(mgl_v, &amgparam_stokes->param_v);
                    break;
                case UA_AMG:
                    fasp_amg_setup_ua(mgl_v, &amgparam_stokes->param_v);
                    break;
                default:
                    printf("Error: Wrong AMG type %d!\n",amgparam_stokes->param_v.AMG_type);
                    exit(ERROR_INPUT_PAR);
            }
            
            
            //-------------------------//
            // setup Schur complement S
            //-------------------------//
            fasp_blas_dcsr_mxm(A_stokes_bcsr.blocks[2], A_stokes_bcsr.blocks[1], &S);
            fasp_blas_dcsr_rap(A_stokes_bcsr.blocks[2], A_stokes_bcsr.blocks[0], A_stokes_bcsr.blocks[1], &BABt);
            
            // change the sign of the BB^T
            fasp_blas_dcsr_axm(&S, -1.0);
            
            // make it non-singular
            INT k,j,ibegin,iend;
            
            for (i=0;i<S.row;++i) {
                ibegin=S.IA[i]; iend=S.IA[i+1];
                for (k=ibegin;k<iend;++k) {
                    j=S.JA[k];
                    if ((j-i)==0) {
                        S.val[k] = S.val[k] + 1e-8; break;
                    } // end if
                } // end for k
            } // end for i
            
            dCSRmat *As = &S;
            const int nnzS = As->nnz;
            mgl_p=fasp_amg_data_create(amgparam_stokes->param_p.max_levels);
            mgl_p[0].A=fasp_dcsr_create(num_pressure,num_pressure,nnzS); fasp_dcsr_cp(As,&mgl_p[0].A);
            mgl_p[0].b=fasp_dvec_create(num_pressure); mgl_p[0].x=fasp_dvec_create(num_pressure);
            // setup AMG
            switch (amgparam_stokes->param_p.AMG_type) {
                case CLASSIC_AMG:
                    fasp_amg_setup_rs(mgl_p, &amgparam_stokes->param_p);
                    break;
                case SA_AMG:
                    fasp_amg_setup_sa(mgl_p, &amgparam_stokes->param_p);
                    break;
                case UA_AMG:
                    fasp_amg_setup_ua(mgl_p, &amgparam_stokes->param_p);
                    break;
                default:
                    printf("Error: Wrong AMG type %d for Schur Complement!\n",amgparam_stokes->param_p.AMG_type);
                    exit(ERROR_INPUT_PAR);
            }
            
        }
        
    }
    
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    // generate data for preconditioners
    // data for pnp
    precond_data_bsr precdata_pnp;
    precdata_pnp.print_level = amgparam_pnp->print_level;
    precdata_pnp.maxit = amgparam_pnp->maxit;
    precdata_pnp.tol = amgparam_pnp->tol;
    precdata_pnp.cycle_type = amgparam_pnp->cycle_type;
    precdata_pnp.smoother = amgparam_pnp->smoother;
    precdata_pnp.presmooth_iter = amgparam_pnp->presmooth_iter;
    precdata_pnp.postsmooth_iter = amgparam_pnp->postsmooth_iter;
    precdata_pnp.coarsening_type = amgparam_pnp->coarsening_type;
    precdata_pnp.relaxation = amgparam_pnp->relaxation;
    precdata_pnp.coarse_scaling = amgparam_pnp->coarse_scaling;
    precdata_pnp.amli_degree = amgparam_pnp->amli_degree;
    precdata_pnp.amli_coef = amgparam_pnp->amli_coef;
    precdata_pnp.tentative_smooth = amgparam_pnp->tentative_smooth;
    precdata_pnp.max_levels = mgl_pnp[0].num_levels;
    precdata_pnp.mgl_data = mgl_pnp;
    precdata_pnp.A = &A_pnp_bsr;
    
    // data for stokes
    // Setup itsolver parameters
    ITS_param ITS_param_v;
    fasp_param_solver_init(&ITS_param_v);
    ITS_param_v.print_level = itparam_stokes->print_level_v;
    ITS_param_v.itsolver_type = itparam_stokes->itsolver_type_v;
    ITS_param_v.restart = itparam_stokes->pre_restart_v;
    ITS_param_v.tol = itparam_stokes->pre_tol_v;
    ITS_param_v.maxit = itparam_stokes->pre_maxit_v;
    ITS_param_v.precond_type = itparam_stokes->precond_type_v;
    
    ITS_param ITS_param_p;
    fasp_param_solver_init(&ITS_param_p);
    ITS_param_p.print_level = itparam_stokes->print_level_p;
    ITS_param_p.itsolver_type = itparam_stokes->itsolver_type_p;
    ITS_param_p.restart = itparam_stokes->pre_restart_p;
    ITS_param_p.tol = itparam_stokes->pre_tol_p;
    ITS_param_p.maxit = itparam_stokes->pre_maxit_p;
    ITS_param_p.precond_type = itparam_stokes->precond_type_p;
    
    // data for stokes
    precond_ns_data  precdata_stokes;
    if ( precond_type > 30 && precond_type < 40 ) {
        precdata_stokes.colA = A_stokes_bcsr.blocks[0]->row;
        precdata_stokes.colB = A_stokes_bcsr.blocks[2]->row;
        precdata_stokes.col  = A_stokes_bcsr.blocks[0]->row + A_stokes_bcsr.blocks[2]->row;
        precdata_stokes.B  = A_stokes_bcsr.blocks[2];
        precdata_stokes.Bt = A_stokes_bcsr.blocks[1];
        precdata_stokes.C = A_stokes_bcsr.blocks[3];
        precdata_stokes.BABt  = &BABt;
    }
    
    precdata_stokes.param_v         = &amgparam_stokes->param_v;
    precdata_stokes.param_p         = &amgparam_stokes->param_p;
    precdata_stokes.ITS_param_v     = &ITS_param_v;
    precdata_stokes.ITS_param_p     = &ITS_param_p;
    precdata_stokes.mgl_data_v      = mgl_v;
    precdata_stokes.mgl_data_p      = mgl_p;
    
    precdata_stokes.max_levels      = mgl_v[0].num_levels;
    precdata_stokes.print_level     = amgparam_stokes->param_v.print_level;
    precdata_stokes.maxit           = amgparam_stokes->param_v.maxit;
    precdata_stokes.amg_tol         = amgparam_stokes->param_v.tol;
    precdata_stokes.cycle_type      = amgparam_stokes->param_v.cycle_type;
    precdata_stokes.smoother        = amgparam_stokes->param_v.smoother;
    precdata_stokes.presmooth_iter  = amgparam_stokes->param_v.presmooth_iter;
    precdata_stokes.postsmooth_iter = amgparam_stokes->param_v.postsmooth_iter;
    precdata_stokes.relaxation      = amgparam_stokes->param_v.relaxation;
    precdata_stokes.coarse_scaling  = amgparam_stokes->param_v.coarse_scaling;
    
    precdata_stokes.S  = &S;
    precdata_stokes.rp = &res_p;
    precdata_stokes.sp = &sol_p;
    
    precdata_stokes.w = (double *)fasp_mem_calloc(precdata_stokes.col,sizeof(double));
    
    // data for overall
    precond_pnp_stokes_data precdata;
    precdata.Abcsr = A;
    
#if WITH_UMFPACK
    // LU if exact solve
    precdata.LU_diag = LU_diag;
#endif

    // pnp part
    precdata.A_pnp_csr = &A_pnp_csr;
    precdata.A_pnp_bsr = &A_pnp_bsr;
    precdata.precdata_pnp = &precdata_pnp;
    precdata.pnp_fct = fasp_precond_dbsr_amg;
    precdata.ILU_pnp = &ILU_pnp;
    
    // stokes part
    precdata.A_stokes_csr = &A_stokes_csr;
    precdata.A_stokes_bcsr = &A_stokes_bcsr;
    precdata.precdata_stokes = &precdata_stokes;
    precdata.stokes_fct = fasp_precond_ns_LSCDGS;
    
    precdata.r = fasp_dvec_create(b->row);
    
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
            
        case 31:
            prec.fct = fasp_precond_pnp_stokes_diag_inexact;
            break;
            
        case 32:
            prec.fct = fasp_precond_pnp_stokes_lower_inexact;
            break;
            
        case 33:
            prec.fct = fasp_precond_pnp_stokes_upper_inexact;
            break;
            
        default:
            fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
            break;
    }
    
    if ( prtlvl >= PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        fasp_cputime("Setup totally", setup_duration);
    }
    
    
    // solver part
    fasp_gettime(&solver_start);
    
    status=fasp_solver_dblc_itsolver(A,b,x, &prec,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( prtlvl >= PRINT_MIN )
        fasp_cputime("Krylov method totally", solver_duration);
    
FINISHED:
    
    // clean
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
#if WITH_UMFPACK
        for (i=0; i<2; i++) fasp_umfpack_free_numeric(LU_diag[i]);
        
        fasp_dcsr_free(&A_pnp_csr);
        fasp_dcsr_free(&A_stokes_csr);
        
        fasp_dvec_free(&precdata.r);
#endif
    }
    /* diagonal blocks are solved by AMG */
    else if (precond_type > 30 && precond_type < 40) {
#if WITH_UMFPACK
        for (i=0; i<2; i++) fasp_umfpack_free_numeric(LU_diag[i]);
#endif
        fasp_dbsr_free(&A_pnp_bsr);
        fasp_amg_data_bsr_free(mgl_pnp);
        //if (&ILU_pnp) fasp_ilu_data_free(&ILU_pnp);
        
        fasp_dblc_free(&A_stokes_bcsr);
        fasp_dcsr_free(&S);
        fasp_dcsr_free(&BABt);
        fasp_amg_data_free(mgl_v, &amgparam_stokes->param_v);
        fasp_amg_data_free(mgl_p, &amgparam_stokes->param_p);
        fasp_dvec_free(&res_p);
        fasp_dvec_free(&sol_p);
        fasp_mem_free(precdata_stokes.w);
        
        fasp_dvec_free(&precdata.r);
        
    }
    
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
