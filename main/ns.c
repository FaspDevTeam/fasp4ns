/*! \file ns.c
 *  \brief Test solvers for (Navier-)Stokes problems
 */

#include <time.h>
#include <math.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test purpose.
 * Lu Wang
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 */
int main (int argc, const char * argv[]) 
{
	block_dCSRmat A;
	dCSRmat RR,RW,WR,WW, Acsr, Mp;
	dvector b, uh, bcsr; //xapp;
    ivector u_idx, p_idx;
	int i,flag=0;
	
    /** initialize block_dCSRmat **/
    A.brow = 2;
    A.bcol = 2;
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    fasp_mem_check((void *)A.blocks, "block matrix:cannot allocate memory!\n", ERROR_ALLOC_MEM);
    A.blocks[0] = &RR;
    A.blocks[1] = &RW;
    A.blocks[2] = &WR;
    A.blocks[3] = &WW;
    
	/** Step 0. Read input parameters */
	char *inputfile = "ini/ns.dat";
	input_ns_param     inparam;  // parameters from input files
	itsolver_ns_param  itparam;  // parameters for itsolver
	AMG_ns_param      amgparam; // parameters for AMG
	ILU_param         iluparam; // parameters for ILU
    Schwarz_param     schparam; // parameters for Schwarz
    
    fasp_ns_param_input(inputfile,&inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &schparam);
	
	// Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int itsolver_type = inparam.solver_type;
	const int precond_type  = inparam.precond_type;
	//const int output_type   = inparam.output_type;
	
	printf("Test Problem %d\n", problem_num);
	
	/** Step 1. Assemble matrix and right-hand side */ 
	char filename1[512];// *datafile1;
	char filename2[512];// *datafile2;
	char filename3[512];// *datafile3;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
	strncpy(filename3,inparam.workdir,128);
	
	/* Assemble A and b. P1 FE discretization for Poisson. */
	if (problem_num == 1) {
        char *fileA = "data/test_1/Matrix_A";
        char *fileB = "data/test_1/Matrix_B";
        char *fileC = "data/test_1/Matrix_C";
        char *filerhs= "data/test_1/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    
	else {
		printf("Error: No such a problem number %d\n", problem_num);
		return -1;
	}
	
	/** Step 2. Check matrix properties */
	
	/** Step 3. Select a solver */
	
	/** Step 4. Solve the system */ 
    fasp_mem_usage();
    
    if (print_level>0) {
		printf("Max it num = %d\n", inparam.itsolver_maxit);
		printf("Tolerance  = %e\n", inparam.itsolver_tol);
	}
    
    // initial guess
	fasp_dvec_alloc(b.row, &uh); fasp_dvec_set(b.row,&uh,0.0);
    
    if (itsolver_type == 0) {
        dvector b_test, uh_test;
        int nrow = A.blocks[0]->row;
        fasp_dvec_alloc(nrow, &b_test); fasp_dvec_set(nrow,&b_test,0.0);
        for (i = 0;i < nrow;i++)
            b_test.val[i] = b.val[i];
        fasp_dvec_alloc(nrow, &uh_test); fasp_dvec_set(nrow,&uh_test,0.0);
		fasp_solver_amg(A.blocks[0], &b_test, &uh_test, &(amgparam.param_v));
	}
    else if (precond_type == PREC_DIAG){
        /*
        if (problem_num > 100 && problem_num < 100 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        /*
        else if ( problem_num > 200 && problem_num < 300 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        //else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        //}
    }
    else if ((precond_type == PREC_UP_TRI)|| (precond_type == PREC_UP_TRI)) {
        /*
        if (problem_num > 100 && problem_num < 100 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        /*
        else if ( problem_num > 200 && problem_num < 300 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        //else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        //}
    }
    else {
        /*
        if (problem_num > 100 && problem_num < 100 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        /*
        else if ( problem_num > 200 && problem_num < 300 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        //else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        //}
    }
	//for (i = 0;i < 4;i++)
    //    fasp_dcsr_free(&A.blocks[i]);
    fasp_dcsr_free(&RR);
    fasp_dcsr_free(&RW);
    fasp_dcsr_free(&WR);
    //fasp_dcsr_free(&WW);
	fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    free(A.blocks);
	printf("Finish\n");
	return flag;
}
