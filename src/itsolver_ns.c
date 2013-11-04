/*! \file itsolver_ns.c
 *  \brief Iterative solvers for (Navier-)ns-type matrices (main file)
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

void fasp_get_schur_somplement(dCSRmat *B,dCSRmat *Bt,dCSRmat *A,dCSRmat *C,dCSRmat *S,dCSRmat *P);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_ns_solver_itsolver(dCSRmat *A, dvector *b, dvector *x, 
 *                                   precond *prec, itsolver_param *itparam)
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
 * \note modified by Xiaozhe on 10/20/2013
 */
int fasp_ns_solver_itsolver(block_dCSRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *prec, 
                               itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	const int itsolver_type = itparam->itsolver_type;
	const int stop_type = itparam->stop_type;
	const double tol = itparam->tol; 
	const int MaxIt = itparam->maxit;
	const int restart = itparam->restart;

	clock_t solver_start=clock();
	int iter;
	switch (itsolver_type) {
        case SOLVER_BiCGstab:
            if (print_level>0) printf("BiCGstab method (Block CSR format) ...\n");
            iter=fasp_solver_bdcsr_pbcgs(A, b, x, prec, tol, MaxIt, stop_type, print_level);
            break;
            
		case SOLVER_MinRes:
			if (print_level>0) printf("Calling MinRes solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pminres(A, b, x, prec, tol, MaxIt,
                                           stop_type,print_level); break;
            
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRES solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pvgmres(A, b, x, prec, tol, MaxIt, restart, stop_type, print_level); break;
        case SOLVER_FGMRES:
			if (print_level>0) printf("Calling FGMRES solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pvfgmres(A, b, x, prec, tol, MaxIt, restart, stop_type, print_level); break;
        case SOLVER_GCR:
            if (print_level>0) printf("Calling GCR solver (Block CSR format) ...\n");
            iter=fasp_solver_bdcsr_pgcr(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart); break;
		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			iter = ERROR_SOLVER_TYPE;
            
	}
	
	if ((print_level>1) && (iter >= 0)) {
		clock_t solver_end=clock();	
		double solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
	
	return iter;
}	

/**
 * \fn int fasp_solver_bdcsr_krylov_navier_stokes (block_dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam, ns_param *precdata)
 * \brief Solve Ax=b by standard Krylov methods for NS equations 
 *
 * \param *A:	       pointer to the block_dCSRmat matrix
 * \param *b:	       pointer to the dvector of right hand side
 * \param *x:	       pointer to the dvector of dofs
 * \param *itparam:  pointer to parameters for iterative solvers
 * \param *precdata: pionter to preconditioner data for ns
 * \return           number of iterations
 *
 * \author Lu Wang
 * \date 03/02/2012
 *
 * \note modified by Xiaozhe on 10/20/2013
 */
int fasp_solver_bdcsr_krylov_navier_stokes (block_dCSRmat *Mat,
                                        dvector *b,
                                        dvector *x,
                                        itsolver_param *itparam,
                                        AMG_param *amgparam,
                                        ILU_param *iluparam,
                                        Schwarz_param *schparam)
{
	// parameters
	const int print_level = itparam->print_level;
	const int precond_type = itparam->precond_type;
    const INT schwarz_mmsize = schparam->schwarz_mmsize;
    const INT schwarz_maxlvl = schparam->schwarz_maxlvl;
    const INT schwarz_type   = schparam->schwarz_type;
	
	// Navier-Stokes 4 by 4 matrix
	dCSRmat *A  = Mat->blocks[0];
	dCSRmat *Bt = Mat->blocks[1];
	dCSRmat *B  = Mat->blocks[2];
	dCSRmat *C  = Mat->blocks[3];
	const int n = A->row, m = B->row, nnzA = A->nnz;
    
    // preconditioner data
	dCSRmat *M = Mat->blocks[3];
	dCSRmat S,P;
    Schwarz_data schwarz_data;
	dvector diag_S;
    
	// local variable
	clock_t solver_start, solver_end, setup_start, setup_end;
	double solver_duration, setup_duration;
	int status=SUCCESS;
        
    //------ setup phase ------//
    setup_start = clock();
    
    //-----------------------//
    // setup AMG for velocity
    //-----------------------//
    AMG_data *mgl=fasp_amg_data_create(amgparam->max_levels);
    mgl[0].A=fasp_dcsr_create(n,n,nnzA); fasp_dcsr_cp(A,&mgl[0].A);
    mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n);
    
    // setup AMG
    switch (amgparam->AMG_type) {
        case CLASSIC_AMG:
            fasp_amg_setup_rs(mgl, amgparam);
            break;
        case SA_AMG:
            fasp_amg_setup_sa(mgl, amgparam);
            break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl, amgparam);
            break;
        default:
            printf("Error: Wrong AMG type %d!\n",amgparam->AMG_type);
            exit(ERROR_INPUT_PAR);
    }
    
    //-------------------------//
    // setup Schur complement S
    //-------------------------//
    //printf("Create chur Complement S\n");
    fasp_get_schur_complement(B,Bt,A,C,&S,&P);
    //printf("Info of S:m = %d, n=%d nnz = %d\n",S.row,S.col,S.nnz);
    //dCSRmat *As = &S;
    dvector res_p = fasp_dvec_create(m);
    dvector sol_p = fasp_dvec_create(m);
    
    fasp_dcsr_getdiag(0,&S,&diag_S);
    
    
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
    
    precdata.max_levels     = mgl[0].num_levels;
	precdata.mgl_data       = mgl;
    precdata.print_level    = amgparam->print_level;
    precdata.maxit          = amgparam->maxit;
    precdata.amg_tol        = amgparam->tol;
    precdata.cycle_type     = amgparam->cycle_type;
    precdata.smoother       = amgparam->smoother;
    precdata.presmooth_iter = amgparam->presmooth_iter;
    precdata.postsmooth_iter= amgparam->postsmooth_iter;
    precdata.relaxation     = amgparam->relaxation;
    precdata.coarse_scaling = amgparam->coarse_scaling;
    
    precdata.S = &S;
    precdata.diag_S = &diag_S;
    precdata.rp = &res_p;
    precdata.sp = &sol_p;

    precdata.w = (double *)fasp_mem_calloc(precdata.col,sizeof(double));
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
            prec.fct = fasp_precond_ns_sym_btri;
            //precdata.w1 = (double *)fasp_mem_calloc(precdata.col,sizeof(double));
            break;
        default:
            printf("Error: Unknown preconditioner type!\n");
            exit(ERROR_SOLVER_PRECTYPE);
	}

	setup_end = clock();
	
	if (print_level>0) {
		setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
		printf("Setup costs %f.\n", setup_duration);
	}
    
	//------ solver phase ------// 
	solver_start=clock();
        status=fasp_ns_solver_itsolver(Mat,b,x,&prec,itparam);
        //status=fasp_ns_solver_itsolver(Mat,b,x,NULL,itparam);
	solver_end=clock();
	
	if (print_level>0) {
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Solver costs %f seconds.\n", solver_duration);	
		printf("Total costs %f seconds.\n", setup_duration + solver_duration);
	}
    
//FINISHED:
	// clean up memory
	if (mgl) fasp_amg_data_free(mgl,amgparam);
	fasp_mem_free(precdata.w);
	fasp_dvec_free(&diag_S);
    fasp_dvec_free(&res_p);
    fasp_dvec_free(&sol_p);
	fasp_dcsr_free (&S);
    //fasp_dcsr_free (&P);
	return status;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov_ns (block_dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam, ns_param *precdata)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param *A:	       pointer to the block_dCSRmat matrix
 * \param *b:	       pointer to the dvector of right hand side
 * \param *x:	       pointer to the dvector of dofs
 * \param *itparam:    pointer to parameters for iterative solvers
 * \param *precdata:   pionter to preconditioner data for ns
 *
 * \return             number of iterations
 *
 * \author Chensong Zhang
 * \date 11/25/2010
 */
INT fasp_solver_bdcsr_krylov_ns (block_dCSRmat *Mat, 
                                     dvector *b, 
                                     dvector *x, 
                                     itsolver_param *itparam,
                                     precond_ns_param *param, 
                                     precond_ns_data *precdata)
{
	// parameters
	const SHORT print_level  = itparam->print_level;
	const SHORT precond_type = itparam->precond_type;
	
	// ns matrix 
	dCSRmat *A = Mat->blocks[0];
	dCSRmat *B = Mat->blocks[1];
	const INT n = A->row, nnzA = A->nnz, m = B->row;	
	
	// preconditioner data
	dCSRmat *M = precdata->M;
	precond prec;
	AMG_param amgparam;
	dvector diag_M;	
	
	// local variable
	clock_t solver_start, solver_end, setup_start, setup_end;
	REAL solver_duration, setup_duration;
	INT status=SUCCESS;
	
	// initialize preconditioner 
	prec.data = &precdata; 
	switch (precond_type) {
		case 1:
			prec.fct = fasp_precond_ns_bdiag;
			break;
		default:
			printf("### ERROR: Unknown preconditioner type!\n");
			exit(ERROR_SOLVER_PRECTYPE);
	}
	
	// AMG parameters
	amgparam.print_level = param->print_level;
	amgparam.max_levels = param->max_levels;
	amgparam.AMG_type = param->AMG_type;
	
	//------ setup phase ------//
	setup_start = clock();
	
	precdata->colA = n;
	precdata->colB = m;
	precdata->col  = n+m;	
	
	// setup work array space
	precdata->w = (REAL *)fasp_mem_calloc(precdata->col,sizeof(REAL));
	// initialize AMG for A
	AMG_data *mgl=fasp_amg_data_create(amgparam.max_levels);
    mgl[0].A=fasp_dcsr_create(n,n,nnzA); fasp_dcsr_cp(A,&mgl[0].A);
	mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n);
	
	// setup AMG
	switch (amgparam.AMG_type) {
		case CLASSIC_AMG:
			fasp_amg_setup_rs(mgl, &amgparam);
			break;
		case SA_AMG:
			fasp_amg_setup_sa(mgl, &amgparam);
			break;
        case UA_AMG:
            fasp_amg_setup_ua(mgl, &amgparam);
		default:
			printf("### ERROR2: Wrong AMG type %d!\n",amgparam.AMG_type);
			exit(ERROR_INPUT_PAR);
	}	
	precdata->max_levels = mgl[0].num_levels;
	precdata->mgl_data = mgl;
	
	// setup diagonal for M
	fasp_dcsr_getdiag(0,M,&diag_M);	
	precdata->diag_M = &diag_M;
	
	setup_end = clock();
	
	if (print_level>PRINT_NONE) {
		setup_duration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Setup costs %f.\n", setup_duration);
	}
	
	//------ solver phase ------// 
	solver_start=clock();
	status=fasp_solver_bdcsr_itsolver(Mat,b,x,&prec,itparam);
	solver_end=clock();
	
	if (print_level>PRINT_NONE) {
		solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Solver costs %f seconds.\n", solver_duration);	
		printf("Total costs %f seconds.\n", setup_duration + solver_duration);
	}
	
	// clean up memory
	if (mgl) fasp_amg_data_free(mgl,&amgparam);
	fasp_mem_free(precdata->w);

	return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
void fasp_get_schur_complement(dCSRmat *B,dCSRmat *Bt,dCSRmat *A,dCSRmat *C,dCSRmat *S,dCSRmat *P)
{
    int colA = A->row;
    //int colB = B->row;
    dvector diag_A;
    int i;
    
    dCSRmat tempA;
    dCSRmat invA = fasp_dcsr_create(colA,colA,colA);
    fasp_dcsr_getdiag(A->row,A,&diag_A);
    for (i=0;i<colA;i++)
    {
        invA.IA[i] = i;
        invA.JA[i] = i;
        if (diag_A.val[i] > SMALLREAL) invA.val[i]   = 1.0/diag_A.val[i];
        else invA.val[i] = 1.0;
        //printf("A[%d] = %e, ",i,diag_A.val[i]);
        //invA.val[i] = 1.0;///diag_A.val[i];
    }
    invA.IA[colA] = colA;
    
    fasp_blas_dcsr_rap (B, &invA, Bt, &tempA);
    //fasp_blas_dcsr_mxm(B,Bt,&tempA);
    fasp_blas_dcsr_add(C,1.0,&tempA,-1.0,S);
    
    //
    //fasp_blas_dcsr_rap (B, &invA, S,  P);
    
    fasp_dvec_free(&diag_A);
    fasp_dcsr_free(&invA);
    fasp_dcsr_free(&tempA);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
