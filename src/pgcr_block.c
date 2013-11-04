/*! \file pgcr.c
 *  \brief Krylov subspace methods -- Preconditioned GCR.
 *
 *
 */  

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"
/**
 *	\fn int fasp_solver_bdcsr_pgcr (block_Reservoir *A, dvector *b, dvector *u, const int MaxIt, const double tol, \
 *        			 precond *pre, const int print_level, const int stop_type, const int restart)
 *	\brief A preconditioned GCR method for solving Au=b 
 *	\param *A	 pointer to the coefficient matrix
 *	\param *b	 pointer to the dvector of right hand side
 *	\param *u	 pointer to the dvector of dofs
 *	\param MaxIt integer, maximal number of iterations
 *	\param tol double float, the tolerance for stopage
 *	\param *pre pointer to the structure of precondition (precond) 
 *\param print_level how much information to print out
 *	\return the number of iterations
 *
 * \author Lu Wang
 * \date 03/02/2012
 *
 */
int fasp_solver_bdcsr_pgcr (block_dCSRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            const int MaxIt, 
                            const double tol,
                            precond *pre, 
                            const int print_level, 
                            const int stop_type, 
                            const int restart)
{
	INT i, j, j1, index;
    INT iter=0;
    INT m=restart;
	const INT nrow=b->row, nrow_1=nrow-1;
	REAL alpha, beta, gamma, tempr, tempe, tempb, tempu,temp2;
	REAL absres0=BIGREAL, absres, relres1, infnormu, factor;
	const REAL sol_inf_tol = 1e-16; // infinity norm tolrance
	
	if ((m<1)||(m>nrow)||(m>150)) m=10; // default restart level
    
	dvector *v=(dvector *)fasp_mem_calloc(m,sizeof(dvector));
	
	for (i=0;i<=m-1;++i) {
		v[i].row=nrow;
		v[i].val=(double*)fasp_mem_calloc(nrow,sizeof(double));
	}
    
    dvector *s=(dvector *)fasp_mem_calloc(m,sizeof(dvector));
	
	for (i=0;i<=m-1;++i) {
		s[i].row=nrow;
		s[i].val=(double*)fasp_mem_calloc(nrow,sizeof(double));
	}
	
    double *r=(double*)fasp_mem_calloc(nrow,sizeof(double));
    
	// compute norm for right hand side
	switch (stop_type) {
		case STOP_REL_RES:
			tempb=fasp_blas_array_norm2(nrow,b->val); // norm(b)
			break;
		case STOP_REL_PRECRES:
			if (pre != NULL) 
				pre->fct(b->val,s[0].val,pre->data); /* Preconditioning */
			else 
				fasp_array_cp(nrow,b->val,s[0].val); /* No preconditioner, B=I */
			tempb=sqrt(ABS(fasp_blas_array_dotprod(nrow,b->val,s[0].val)));
			break;
		case STOP_MOD_REL_RES:
			break;
		default:
			printf("Error: Unknown stopping criteria!\n");
			iter = ERROR_INPUT_PAR;
			goto FINISHED;
	}
	tempb=MAX(SMALLREAL,tempb);
	tempu=MAX(SMALLREAL,fasp_blas_array_norm2(nrow,u->val));
	
	// r = b-A*u
	fasp_array_cp(nrow,b->val,r); 
	fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
	tempe=fasp_blas_array_norm2(nrow,r);
	tempb=MAX(SMALLREAL,tempe);
	switch (stop_type) {
		case STOP_REL_RES:
			relres1=tempe/tempb; 
			break;
		case STOP_REL_PRECRES:
			if (pre == NULL)
				fasp_array_cp(nrow,r,s[0].val);
			else
				pre->fct(r,s[0].val,pre->data);
			temp2=sqrt(ABS(fasp_blas_array_dotprod(nrow,r,s[0].val)));
			relres1=temp2/tempb; 
			break;
		case STOP_MOD_REL_RES:
			relres1=tempe/tempu; 
			break;
	}
	
	if (relres1<tol) { fasp_mem_free(r); goto FINISHED; }
	
	if  (iter < 0) goto FINISHED;
    
	while (iter++<MaxIt)
	{      
		for (j=0;j<m;++j)
        {
            if (pre == NULL) {
                fasp_array_cp(nrow,r,s[j].val);  /* No preconditioner, B=I */
            }
            else {
                pre->fct(r,s[j].val,pre->data); /* Preconditioning */
            }
		
            fasp_blas_bdcsr_aAxpy(1.0,A,s[j].val,v[j].val);			
			for (i=0;i<j;++i)
			{
				alpha = fasp_blas_array_dotprod(nrow,v[j].val,v[i].val);
				fasp_blas_array_axpy(nrow, -alpha, v[i].val, v[j].val);
                fasp_blas_array_axpy(nrow, -alpha, s[i].val, s[j].val);
			}
			
			beta = fasp_blas_array_norm2(nrow,v[j].val);
            fasp_blas_array_ax(nrow,1.0/beta,v[j].val);
            fasp_blas_array_ax(nrow,1.0/beta,s[j].val);
            
            gamma=fasp_blas_array_dotprod(nrow,v[j].val,r);
            fasp_blas_array_axpy(nrow, gamma, s[j].val, u->val);
            //fasp_blas_array_axpy(nrow, -gamma, v[j].val, r);
            
            // r = b-A*u
            fasp_array_cp(nrow,b->val,r); 
            fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
            
            // absolute and relative residuals
            absres=sqrt(fasp_blas_array_dotprod(nrow,r,r));
            tempu=sqrt(fasp_blas_dvec_dotprod(u,u));		
            switch (stop_type) {
                case STOP_REL_RES:
                    relres1=absres/tempb; 
                    break;
                case STOP_REL_PRECRES:
                    if (pre == NULL)
                        fasp_array_cp(nrow,r,s[j].val);
                    else
                        pre->fct(r,s[j].val,pre->data);
                    temp2=sqrt(ABS(fasp_blas_array_dotprod(nrow,r,s[j].val)));
                    relres1=temp2/tempb; 
                    break;
                case STOP_MOD_REL_RES:
                    relres1=absres/tempu; 
                    break;
            }
            
            // contraction factor
            factor=absres/absres0;
            
            // output iteration information if needed	
            print_itinfo(print_level,stop_type,(iter-1)*m+j+1,relres1,absres,factor);
            
            absres0=absres;
            
            // solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
            infnormu = fasp_blas_array_norminf(nrow, u->val); 
            if (infnormu <= sol_inf_tol)
            {
                print_message(print_level, "GMRes stops: infinity norm of the solution is too small!\n");
                iter = ERROR_SOLVER_SOLSTAG;
                goto FINISHED;
            }
            
            if (relres1<tol) goto FINISHED;
		}		
	}
	
FINISHED:
	if (print_level>0) {
		if (iter>MaxIt){
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
			iter = ERROR_SOLVER_MAXIT;
		}
		else
			printf("Number of iterations = inner:%d(outer:%d) with relative residual %e.\n", j+1,iter, relres1);
	}
	
	fasp_mem_free(r);
	for (i=0;i<=m-1;++i) fasp_mem_free(v[i].val);
	fasp_mem_free(v);
    for (i=0;i<=m-1;++i) fasp_mem_free(s[i].val);
	fasp_mem_free(s);
	
	return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
