/*! \file precond_ns.c
 *  \brief Preconditioners for (Navier-)Stokes problems
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"
/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_ns_bdiag (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning for ns equation
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 *
 */
void fasp_precond_ns_bdiag (REAL *r, 
                            REAL *z,
                            void *data)
{
	precond_ns_data *predata=(precond_ns_data *)data;
	
	const INT col = predata->col, colA = predata->colA, colB = predata->colB;
	
	//! prepare	AMG preconditioner 
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
	
	//! setup z;
	fasp_array_set(col, z, 0.0);
	
    //-------------------
	//! Solve velocity
    //-------------------
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,&amgparam);
	pcdata.max_levels = mgl[0].num_levels;
	pcdata.mgl_data = predata->mgl_data;
	precond pc; pc.data = &pcdata;
	pc.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);

    //-------------------------
	//! Solve Schur complement
    //-------------------------
    precond pc_s;
    pc_s.data = predata->diag_S;
    pc_s.fct  = fasp_precond_diag;

    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);

}

/**
 * \fn void fasp_precond_ns_low_btri (double *r, double *z, void *data)
 * \brief upper block diagonal preconditioning for ns equation
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 *
 */
void fasp_precond_ns_low_btri (REAL *r,
                                  REAL *z, 
                                  void *data)
{
	precond_ns_data *predata=(precond_ns_data *)data;
	const int col = predata->col, colA = predata->colA, colB = predata->colB;
	//const int maxit = predata->maxit;
	//double *diagptr=predata->diag_S->val;
    
	// local variables
	double	*tempr = predata->w;		
    
	//! prepare	AMG preconditioner
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
    
	//! back up r, setup z;
	fasp_array_cp(col, r, tempr);
	fasp_array_set(col, z, 0.0);
    
    //-------------------
	//! Solve velocity
    //-------------------
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,&amgparam);
	pcdata.max_levels = mgl[0].num_levels;
	pcdata.mgl_data = predata->mgl_data;
	precond pc; pc.data = &pcdata;
	pc.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    
    /*
    dCSRmat tmpA;
    dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(&mgl[0].A,ptrA);
    fasp_solver_umfpack(ptrA, &rv, &zv, 0);
    fasp_dcsr_free(ptrA);
     */

    //-------------------
    //! Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
	//! Solve Schur complement
    //-------------------------
    precond pc_s;
    pc_s.data = predata->diag_S;
    pc_s.fct  = fasp_precond_diag;
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);
    
    /*
    //dCSRmat tmpA;
    //dCSRmat *ptrA = &tmpA;
    fasp_dcsr_trans(predata->S,ptrA);
    fasp_solver_umfpack(ptrA, &rs, &zs, 0);
    fasp_dcsr_free(ptrA);
     */
	
    //! restore r
	fasp_array_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_up_btri (double *r, double *z, void *data)
 * \brief upper block diagonal preconditioning for ns equation
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 *
 */
void fasp_precond_ns_up_btri (REAL *r,
                              REAL *z,
                              void *data)
{
	precond_ns_data *predata=(precond_ns_data *)data;
	const int col = predata->col, colA = predata->colA, colB = predata->colB;
	//const int maxit = predata->maxit;
	//double *diagptr=predata->diag_S->val;
    
	// local variables
	double	*tempr = predata->w;
    
	//! prepare	AMG preconditioner
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
    
	//! back up r, setup z;
	fasp_array_cp(col, r, tempr);
	fasp_array_set(col, z, 0.0);
    
    //-------------------------
	//! Solve Schur complement
    //-------------------------
    precond pc_s;
    pc_s.data = predata->diag_S;
    pc_s.fct  = fasp_precond_diag;
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);
    
    //-------------------
    //! Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
	//! Solve velocity
    //-------------------
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,&amgparam);
	pcdata.max_levels = mgl[0].num_levels;
	pcdata.mgl_data = predata->mgl_data;
	precond pc; pc.data = &pcdata;
	pc.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    
	//! restore r
	fasp_array_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_sym_btri (double *r, double *z, void *data)
 * \brief upper block diagonal preconditioning for ns equation
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiozhe Hu, Lu Wang
 * \date 10/20/2013
 *
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
    
	//! prepare	AMG preconditioner
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
    
	//! back up r, setup z;
	fasp_array_cp(col, r, tempr);
	fasp_array_set(col, z, 0.0);
    //-------------------
	//! Solve velocity
    //-------------------
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,&amgparam);
	pcdata.max_levels = mgl[0].num_levels;
	pcdata.mgl_data = predata->mgl_data;
	precond pc; pc.data = &pcdata;
	pc.fct = fasp_precond_amg;
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    
    //-------------------
    //! Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->B, zv.val, rs.val);
    
    //-------------------------
	//! Solve Schur complement
    //-------------------------
    precond pc_s;
    pc_s.data = predata->diag_S;
    pc_s.fct  = fasp_precond_diag;
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);

    //-------------------
    //! Compute residule
    //-------------------
    //fasp_blas_dcsr_aAxpy(1.0, predata->C, zs.val, rs.val);
    
    //-------------------------
	//! Solve Schur complement
    //-------------------------
    
    fasp_solver_dcsr_pvfgmres(predata->S, &rs, &zs, &pc_s, 1.0e-2, 20, 20, 1, 0);
    
    //-------------------
    //! Compute residule
    //-------------------
    fasp_blas_dcsr_aAxpy(-1.0, predata->Bt, zs.val, rv.val);
    
    //-------------------
	//! Solve velocity
    //-------------------
    
    fasp_solver_dcsr_pvfgmres(&mgl[0].A, &rv, &zv, &pc, 1.0e-2, 20, 20, 1, 0);
    
    //fasp_blas_dvec_axpy(1.0,&z1v,&zv);
    //fasp_blas_dvec_axpy(1.0,&z1s,&zs);
    
	//! restore r
	fasp_array_cp(col, tempr, r);
}

/**
 * \fn void fasp_precond_ns_lsc (REAL *r, REAL *z, void *data)
 * \brief LSC preconditioning for ns equation
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
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
	//! prepare	 AMG preconditioner for A
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
    
    //! prepare iterative parameter
    itsolver_param  itparam;fasp_param_solver_init (&itparam);
    itparam.print_level = 0;
    itparam.itsolver_type = SOLVER_VFGMRES;
    itparam.precond_type   = 3;
    itparam.maxit          = 100;
    itparam.restart        = 100;
    itparam.tol            = predata->amg_tol;
    
	//! prepare	 AMG preconditioner for S
	//AMG_data *mgl_s = predata->mgl_data_p;
	AMG_param amgparam_s; fasp_param_amg_init(&amgparam_s);
	amgparam_s.cycle_type = predata->cycle_type;
	amgparam_s.smoother   = predata->smoother;
	amgparam_s.presmooth_iter  = predata->presmooth_iter;
	amgparam_s.postsmooth_iter = predata->postsmooth_iter;
	amgparam_s.relaxation      = predata->relaxation;
	amgparam_s.coarse_scaling  = predata->coarse_scaling;
	amgparam_s.ILU_levels      = predata->mgl_data_p->ILU_levels;
    
	//! back up r, setup z;
	fasp_array_cp(col, r, tempr);
	fasp_array_set(col, z, 0.0);
    
    //! Solve M
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
    
    //tmp  = fasp_blas_array_norm2(col,z);
    //printf("norm of z2=%e\n",tmp);
    
	//! Solve A by AMG
    mgl->b.row=colA; fasp_array_cp(colA,r,mgl->b.val); // residual is an input 
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
    
	//! restore r
	fasp_array_cp(col, tempr, r);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
