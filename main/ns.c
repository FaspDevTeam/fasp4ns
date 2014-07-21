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
	else if (problem_num == 2) {
        char *fileA = "data/test_2/Matrix_A";
        char *fileB = "data/test_2/Matrix_B";
        char *fileC = "data/test_2/Matrix_C";
        char *filerhs= "data/test_2/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    else if (problem_num == 3) {
        char *fileA = "data/test_3/Matrix_A";
        char *fileB = "data/test_3/Matrix_B";
        char *fileC = "data/test_3/Matrix_C";
        char *filerhs= "data/test_3/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    else if (problem_num == 4) {
        char *fileA = "data/test_4/Matrix_A";
        char *fileB = "data/test_4/Matrix_B";
        char *fileC = "data/test_4/Matrix_C";
        char *filerhs= "data/test_4/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    else if (problem_num == 5) {
        char *fileA = "data/test_5/Matrix_A";
        char *fileB = "data/test_5/Matrix_B";
        char *fileC = "data/test_5/Matrix_C";
        char *filerhs= "data/test_5/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    else if (problem_num == 6) {
        char *fileA = "data/test_6/Matrix_A";
        char *fileB = "data/test_6/Matrix_B";
        char *fileC = "data/test_6/Matrix_C";
        char *filerhs= "data/test_6/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
	
    else if (problem_num == 10) {
        
        char *fileA = "data/test_1/A.dat";
        char *fileB = "data/test_1/B.dat";
        char *fileBt = "data/test_1/Bt.dat";
        char *fileC = "data/test_1/C.dat";
        char *filerhs= "data/test_1/rhs.dat";
        
        fasp_dcoo_read(fileA,&RR);
        fasp_dcoo_read(fileBt,&RW);
        fasp_dcoo_read(fileB,&WR);
        fasp_dcoo_read(fileC,&WW);
        
        fasp_blas_dcsr_axm(&WW,1e+3);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    else if (problem_num == 11) {
        
        char *fileA = "data/data4Xiaozhe/A1.mat";
        char *fileB = "data/data4Xiaozhe/B1.mat";
        char *fileBt = "data/data4Xiaozhe/Bt1.mat";
        char *fileC = "data/data4Xiaozhe/C1.mat";
        char *filerhs= "data/data4Xiaozhe/rhs1.dat";
        
        fasp_dcoo_read(fileA,&RR);
        fasp_dcoo_read(fileBt,&RW);
        fasp_dcoo_read(fileB,&WR);
        fasp_dcoo_read(fileC,&WW);
        
        
        fasp_dvec_read(filerhs,&b);
        
    }
    else if (problem_num == 12) {
        
        char *fileA = "data/data4Xiaozhe/A2.mat";
        char *fileB = "data/data4Xiaozhe/B2.mat";
        char *fileBt = "data/data4Xiaozhe/Bt2.mat";
        char *fileC = "data/data4Xiaozhe/C2.mat";
        char *filerhs= "data/data4Xiaozhe/rhs2.dat";
        
        fasp_dcoo_read(fileA,&RR);
        fasp_dcoo_read(fileBt,&RW);
        fasp_dcoo_read(fileB,&WR);
        fasp_dcoo_read(fileC,&WW);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    else if (problem_num == 13) {
        
        char *fileA = "data/data4Xiaozhe/A3.mat";
        char *fileB = "data/data4Xiaozhe/B3.mat";
        char *fileBt = "data/data4Xiaozhe/Bt3.mat";
        char *fileC = "data/data4Xiaozhe/C3.mat";
        char *filerhs= "data/data4Xiaozhe/rhs3.dat";
        
        fasp_dcoo_read(fileA,&RR);
        fasp_dcoo_read(fileBt,&RW);
        fasp_dcoo_read(fileB,&WR);
        fasp_dcoo_read(fileC,&WW);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    else if (problem_num == 14) {
        
        char *fileA = "data/data4Xiaozhe/A4.mat";
        char *fileB = "data/data4Xiaozhe/B4.mat";
        char *fileBt = "data/data4Xiaozhe/Bt4.mat";
        char *fileC = "data/data4Xiaozhe/C4.mat";
        char *filerhs= "data/data4Xiaozhe/rhs4.dat";
        
        fasp_dcoo_read(fileA,&RR);
        fasp_dcoo_read(fileBt,&RW);
        fasp_dcoo_read(fileB,&WR);
        fasp_dcoo_read(fileC,&WW);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    
    else if (problem_num == 20) {
        char *fileA = "data/out/data_0/A11.dat";
        char *fileB = "data/out/data_0/A12.dat";
        char *fileC = "data/out/data_0/A21.dat";
        char *fileD = "data/out/data_0/A22.dat";
        char *filerhs= "data/out/data_0/b.dat";
        char *filex0= "data/out/data_0/x.dat";
        fasp_bdcsr_read_ruth(fileA,fileB,fileC,fileD,filerhs,filex0,&A,&b,&uh);
    }
    
    else if (problem_num == 47) {
        char *fileA = "data/3DUnstructuredFSI/Matrix_A";
        char *fileB = "data/3DUnstructuredFSI/Matrix_B";
        char *fileC = "data/3DUnstructuredFSI/Matrix_C";
        char *filerhs= "data/3DUnstructuredFSI/RHS";
        fasp_bdcsr_read(fileA,fileB,fileC,filerhs,&A,&b);
	}
    
    // --------------------------- //
    //  Examples wih pressure mass
    // --------------------------- //

    else if (problem_num == 101) {
        
        char *fileA   = "data/P2P0-with-divdiv/set-1/Auu.dat";
        char *fileB   = "data/P2P0-with-divdiv/set-1/Apu.dat";
        char *fileBt  = "data/P2P0-with-divdiv/set-1/Aup.dat";
        char *fileC   = "data/P2P0-with-divdiv/set-1/App.dat";
        char *filerhs = "data/P2P0-with-divdiv/set-1/b.dat";
        
        char *fileMp   = "data/P2P0-with-divdiv/set-1/Mp.dat";
        
        fasp_dcoo_read(fileA, A.blocks[0]);
        fasp_dcoo_read(fileBt,A.blocks[1]);
        fasp_dcoo_read(fileB, A.blocks[2]);
        fasp_dcoo_read(fileC, A.blocks[3]);
        
        fasp_dcoo_read(fileMp, &Mp);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    
    else if (problem_num == 102) {
        
        char *fileA   = "data/P2P0-with-divdiv/set-2/Auu.dat";
        char *fileB   = "data/P2P0-with-divdiv/set-2/Apu.dat";
        char *fileBt  = "data/P2P0-with-divdiv/set-2/Aup.dat";
        char *fileC   = "data/P2P0-with-divdiv/set-2/App.dat";
        char *filerhs = "data/P2P0-with-divdiv/set-2/b.dat";
        
        char *fileMp   = "data/P2P0-with-divdiv/set-2/Mp.dat";
        
        fasp_dcoo_read(fileA, A.blocks[0]);
        fasp_dcoo_read(fileBt,A.blocks[1]);
        fasp_dcoo_read(fileB, A.blocks[2]);
        fasp_dcoo_read(fileC, A.blocks[3]);
        
        fasp_dcoo_read(fileMp, &Mp);
        
        fasp_dvec_read(filerhs,&b);
        
    }
    
    else if (problem_num == 201) {
        
        // read in data
        char *fileA = "data/P2P0-no-divdiv/set-1/A.dat";
        char *filerhs = "data/P2P0-no-divdiv/set-1/b.dat";
        char *fileMp   = "data/P2P0-no-divdiv/set-1/Mp.dat";

        fasp_dcoo_shift_read(fileA, &Acsr);
        fasp_dvec_read(filerhs,&bcsr);
        fasp_dcoo_shift_read(fileMp, &Mp);

        // read in index for velocity and pressure
        char *file_u_idx = "data/P2P0-no-divdiv/set-1/ISu.dat";
        char *file_p_idx = "data/P2P0-no-divdiv/set-1/ISp.dat";

        fasp_ivecind_read(file_u_idx, &u_idx);
        fasp_ivecind_read(file_p_idx, &p_idx);
        
        // get Auu block
        fasp_dcsr_getblk(&Acsr, u_idx.val, u_idx.val, u_idx.row, u_idx.row, A.blocks[0]);
        // get Aup block
        fasp_dcsr_getblk(&Acsr, u_idx.val, p_idx.val, u_idx.row, p_idx.row, A.blocks[1]);
        // get Apu block
        fasp_dcsr_getblk(&Acsr, p_idx.val, u_idx.val, p_idx.row, u_idx.row, A.blocks[2]);
        // get App block
        fasp_dcsr_getblk(&Acsr, p_idx.val, p_idx.val, p_idx.row, p_idx.row, A.blocks[3]);

        // get right hand side
        fasp_dvec_alloc(bcsr.row, &b);
        
        for (i=0; i<u_idx.row; i++) b.val[i] = bcsr.val[u_idx.val[i]];
        for (i=u_idx.row; i<bcsr.row; i++) b.val[i] = bcsr.val[p_idx.val[i-u_idx.row]];

        
    }
    
    else if (problem_num == 202) {
        
        // read in data
        char *fileA = "data/P2P0-no-divdiv/set-2/A.dat";
        char *filerhs = "data/P2P0-no-divdiv/set-2/b.dat";
        char *fileMp   = "data/P2P0-no-divdiv/set-2/Mp.dat";
        
        fasp_dcoo_shift_read(fileA, &Acsr);
        fasp_dvec_read(filerhs,&bcsr);
        fasp_dcoo_shift_read(fileMp, &Mp);
        
        // read in index for velocity and pressure
        char *file_u_idx = "data/P2P0-no-divdiv/set-2/ISu.dat";
        char *file_p_idx = "data/P2P0-no-divdiv/set-2/ISp.dat";
        
        fasp_ivecind_read(file_u_idx, &u_idx);
        fasp_ivecind_read(file_p_idx, &p_idx);
        
        // get Auu block
        fasp_dcsr_getblk(&Acsr, u_idx.val, u_idx.val, u_idx.row, u_idx.row, A.blocks[0]);
        // get Aup block
        fasp_dcsr_getblk(&Acsr, u_idx.val, p_idx.val, u_idx.row, p_idx.row, A.blocks[1]);
        // get Apu block
        fasp_dcsr_getblk(&Acsr, p_idx.val, u_idx.val, p_idx.row, u_idx.row, A.blocks[2]);
        // get App block
        fasp_dcsr_getblk(&Acsr, p_idx.val, p_idx.val, p_idx.row, p_idx.row, A.blocks[3]);
        
        // get right hand side
        fasp_dvec_alloc(bcsr.row, &b);
        
        for (i=0; i<u_idx.row; i++) b.val[i] = bcsr.val[u_idx.val[i]];
        for (i=u_idx.row; i<bcsr.row; i++) b.val[i] = bcsr.val[p_idx.val[i-u_idx.row]];
        
        
    }
    
    else if (problem_num == 203) {
        
        // read in data
        char *fileA = "data/P2P0-no-divdiv/set-3/A.dat";
        char *filerhs = "data/P2P0-no-divdiv/set-3/b.dat";
        char *fileMp   = "data/P2P0-no-divdiv/set-3/Mp.dat";
        
        fasp_dcoo_shift_read(fileA, &Acsr);
        fasp_dvec_read(filerhs,&bcsr);
        fasp_dcoo_shift_read(fileMp, &Mp);
        
        // read in index for velocity and pressure
        char *file_u_idx = "data/P2P0-no-divdiv/set-3/ISu.dat";
        char *file_p_idx = "data/P2P0-no-divdiv/set-3/ISp.dat";
        
        fasp_ivecind_read(file_u_idx, &u_idx);
        fasp_ivecind_read(file_p_idx, &p_idx);
        
        // get Auu block
        fasp_dcsr_getblk(&Acsr, u_idx.val, u_idx.val, u_idx.row, u_idx.row, A.blocks[0]);
        // get Aup block
        fasp_dcsr_getblk(&Acsr, u_idx.val, p_idx.val, u_idx.row, p_idx.row, A.blocks[1]);
        // get Apu block
        fasp_dcsr_getblk(&Acsr, p_idx.val, u_idx.val, p_idx.row, u_idx.row, A.blocks[2]);
        // get App block
        fasp_dcsr_getblk(&Acsr, p_idx.val, p_idx.val, p_idx.row, p_idx.row, A.blocks[3]);
        
        // get right hand side
        fasp_dvec_alloc(bcsr.row, &b);
        
        for (i=0; i<u_idx.row; i++) b.val[i] = bcsr.val[u_idx.val[i]];
        for (i=u_idx.row; i<bcsr.row; i++) b.val[i] = bcsr.val[p_idx.val[i-u_idx.row]];
        
        
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
        if (problem_num > 100 && problem_num < 200 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else if ( problem_num > 200 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        }
    }
    else if ((precond_type == PREC_UP_TRI)|| (precond_type == PREC_UP_TRI)) {
        if (problem_num > 100 && problem_num < 200 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else if ( problem_num > 200 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        }
    }
    else {
        if (problem_num > 100 && problem_num < 200 ){
            flag = fasp_solver_bdcsr_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else if ( problem_num > 200 ) {
            flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
        else {
            flag = fasp_solver_bdcsr_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        }
    }
	//for (i = 0;i < 4;i++)
    //    fasp_dcsr_free(&A.blocks[i]);
    fasp_dcsr_free(&RR);
    fasp_dcsr_free(&RW);
    fasp_dcsr_free(&WR);
    fasp_dcsr_free(&WW);
	fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    free(A.blocks);
	printf("Finish\n");
	return flag;
}
