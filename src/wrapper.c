/*! \file wrapper.c
 *  \brief Wrappers for accessing functions by advanced users.
 *
 *  \note Input variables shall not need fasp.h at all!
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"
/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

#if 0

/**
 * \fn void fasp_krylov_ns_ (INT *nA, INT *nnzA, INT *ia, INT *ja, REAL *aval,
 *                               INT *nB, INT *nnzB, INT *ib, INT *jb, REAL *bval,
 *                               INT *nM, INT *nnzM, INT *im, INT *jm, REAL *mval,
 *                               REAL *b, REAL *u, REAL *beta,
 *                               REAL *tol, INT *maxit, INT *ptrlvl)
 *
 * \brief Solve [A B;B' O]x=b by Krylov method with block diagonal preconditioner
 *
 * \param nA       num of cols of A
 * \param nnzA     num of nonzeros of A
 * \param ia       IA of A in CSR format
 * \param ja       JA of A in CSR format
 * \param aval     VAL of A in CSR format
 * \param nB       num of cols of B
 * \param nnzB     num of nonzeros of B
 * \param ib       IA of B in CSR format
 * \param jb       JA of B in CSR format
 * \param bval     VAL of B in CSR format
 * \param nM       num of cols of M
 * \param nnzM     num of nonzeros of M
 * \param im       IA of M in CSR format
 * \param jm       JA of M in CSR format
 * \param mval     VAL of M in CSR format
 * \param nP       num of cols of P
 * \param nnzP     num of nonzeros of P
 * \param ip       IA of P in CSR format
 * \param jp       JA of P in CSR format
 * \param pval     VAL of P in CSR format
 * \param b        rhs vector
 * \param u        solution vector
 * \param tol      tolerance for iterative solvers
 * \param maxit    max num of iterations
 * \param ptrlvl   print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date 11/16/2010
 */
void fasp_fwrapper_krylov_ns_ (INT *nA,
                               INT *nnzA,
                               INT *ia,
                               INT *ja,
                               REAL *aval,
                               INT *nB,
                               INT *nnzB,
                               INT *ib,
                               INT *jb,
                               REAL *bval,
                               INT *nM,
                               INT *nnzM,
                               INT *im,
                               INT *jm,
                               REAL *mval,
                               INT *nP,
                               INT *nnzP,
                               INT *ip,
                               INT *jp,
                               REAL *pval,
                               REAL *b,
                               REAL *u,
                               REAL *beta,
                               REAL *tol,
                               INT *maxit,
                               INT *ptrlvl)
{
    dCSRmat matA, matB, matBt, matM, matP, zero;
    dBLCmat mat; // coefficient matrix
    dvector rhs, sol; // right-hand-side, solution
    itsolver_param  itparam;  // parameters for itsolver
    precond_ns_param psparam; // parameters for ns precond
    precond_ns_data  psdata; // data for ns precond
    
    // initialize itsolver parameters
    fasp_param_solver_init(&itparam);
    itparam.itsolver_type = SOLVER_MinRes;
    itparam.tol           = *tol;
    itparam.print_level   = *ptrlvl;
    itparam.maxit         = *maxit;
    
    // initialize precond parameters
    psparam.AMG_type = CLASSIC_AMG;
    psparam.print_level = *ptrlvl;
    psparam.max_levels = 15;
    
    // initialize matrix
    matA.row = *nA; matA.col = *nA; matA.nnz = *nnzA;
    matA.IA  = ia;  matA.JA  = ja; matA.val = aval;
    
    matB.row = *nB; matB.col = *nB; matB.nnz = *nnzB;
    matB.IA  = ib;  matB.JA  = jb; matB.val = bval;
    
    fasp_dcsr_trans(&matB, &matBt);
    mat.brow = mat.bcol = 2;
    mat.blocks[0] = &matA;
    mat.blocks[1] = &matB;
    mat.blocks[2] = &matBt;
    mat.blocks[3] = &zero;
    
    matM.row = *nB; matM.col = *nB; matM.nnz = *nnzB;
    matM.IA  = ib;  matM.JA  = jb;  matM.val = bval;
    
    matP.row = *nP; matP.col = *nP; matP.nnz = *nnzP;
    matP.IA  = ip;  matP.JA  = jp;  matP.val = pval;
    
    rhs.row = *nA+*nB; rhs.val = b;
    sol.row = *nA+*nB; sol.val = u;
    
    // initialize precond data
    psdata.colA = *nA;
    psdata.colB = *nB;
    psdata.col = *nA+*nB;
    psdata.beta = *beta;
    psdata.M = &matM;
    psdata.P = &matP;
    
    fasp_solver_dblc_krylov_ns(&mat, &rhs, &sol, &itparam, &psparam, &psdata);
}

#endif

/**
 * \fn void fasp_fwrapper_krylov_navier_stokes_ (INT *nA, INT *nnzA, INT *ia, INT *ja, REAL *aval,
 *                               INT *nB, INT *nnzB, INT *ib, INT *jb, REAL *bval,
 *                               INT *nM, INT *nnzM, INT *im, INT *jm, REAL *mval,
 *                               REAL *b, REAL *u, REAL *beta,
 *                               REAL *tol, INT *maxit, INT *ptrlvl)
 *
 * \brief Solve [A B;B' O]x=b by Krylov method with block diagonal preconditioner
 *
 * \param nA       num of cols of A
 * \param nnzA     num of nonzeros of A
 * \param ia       IA of A in CSR format
 * \param ja       JA of A in CSR format
 * \param aval     VAL of A in CSR format
 * \param nB       num of cols of B
 * \param nnzB     num of nonzeros of B
 * \param ib      IA of B in CSR format
 * \param jb      JA of B in CSR format
 * \param bval     VAL of B in CSR format
 * \param nC       num of cols of C
 * \param nnzC     num of nonzeros of C
 * \param ic       IA of C in CSR format
 * \param jc       JA of C in CSR format
 * \param cval     VAL of C in CSR format
 * \param b        rhs vector
 * \param u        solution vector
 *
 * \author Lu Wang
 * \date 03/14/2012
 */
void fasp_fwrapper_krylov_navier_stokes_ (INT *nA,
                                          INT *nnzA,
                                          INT *ia,
                                          INT *ja,
                                          REAL *aval,
                                          INT *nB,
                                          INT *nnzB,
                                          INT *ib,
                                          INT *jb,
                                          REAL *bval,
                                          INT *nC,
                                          INT *nnzC,
                                          INT *ic,
                                          INT *jc,
                                          REAL *cval,
                                          REAL *b,
                                          REAL *u)
{
    dBLCmat A; // coefficient matrix
    dCSRmat matA, matB, matBt, matC;
    dvector rhs, sol; // right-hand-side, solution
    precond_ns_param psparam; // parameters for ns precond
    precond_ns_data  psdata; // data for ns precond
    int i,flag;
    
    /** initialize dBLCmat **/
    A.brow = 2;
    A.bcol = 2;
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    fasp_mem_check((void *)A.blocks, "block matrix:cannot allocate memory!\n", ERROR_ALLOC_MEM);
    A.blocks[0] = &matA;
    A.blocks[1] = &matBt;
    A.blocks[2] = &matB;
    A.blocks[3] = &matC;
    
    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    //char *inputfile = "fasp4ns.dat";
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
    
    // initialize matrix
    matA.row = *nA; matA.col = *nA; matA.nnz = *nnzA;
    matA.IA  = ia;  matA.JA  = ja; matA.val = aval;
    matB.row = *nB; matB.col = *nA; matB.nnz = *nnzB;
    matB.IA  = ib;  matB.JA  = jb; matB.val = bval;
    matC.row = *nC; matC.col = *nC; matC.nnz = *nnzC;
    matC.IA  = ic;  matC.JA  = jc; matC.val = cval;
    
    //  Shift the index to start from 0 (for C routines)
    for (i=0;i<*nnzA;i++) matA.JA[i] --;
    for (i=0;i<*nnzB;i++) matB.JA[i] --;
    for (i=0;i<*nnzC;i++) matC.JA[i] --;
    
    fasp_dcsr_trans(&matB, &matBt);
    
    rhs.row = *nA+*nB; rhs.val = b;
    sol.row = *nA+*nB; sol.val = u;
    printf(" Finish Construct A\n");
    fasp_mem_usage();
    
    if (print_level>0) {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }
    
    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam, &amgparam, &iluparam, &schparam);
}


/**
 * \fn void fasp_fwrapper_krylov_navier_stokes_nsym_ (INT *nA, INT *nnzA, INT *ia, INT *ja, REAL *aval,
 *                               INT *nB, INT *nnzB, INT *ib, INT *jb, REAL *bval,
 *                               INT *nM, INT *nnzM, INT *im, INT *jm, REAL *mval,
 *                               REAL *b, REAL *u, REAL *beta,
 *                               REAL *tol, INT *maxit, INT *ptrlvl)
 *
 * \brief Solve [A B;B' O]x=b by Krylov method with block diagonal preconditioner
 *
 * \param nA       num of cols of A
 * \param nnzA     num of nonzeros of A
 * \param ia       IA of A in CSR format
 * \param ja       JA of A in CSR format
 * \param aval     VAL of A in CSR format
 * \param nB       num of cols of B
 * \param nnzB     num of nonzeros of B
 * \param ib      IA of B in CSR format
 * \param jb      JA of B in CSR format
 * \param bval     VAL of B in CSR format
 * \param nC       num of cols of C
 * \param nnzC     num of nonzeros of C
 * \param ic       IA of C in CSR format
 * \param jc       JA of C in CSR format
 * \param cval     VAL of C in CSR format
 * \param b        rhs vector
 * \param u        solution vector
 *
 * \author Lu Wang
 * \date 03/20/2014
 */
void fasp_fwrapper_krylov_navier_stokes_nsym_ (INT *nA,
                                               INT *nnzA,
                                               INT *ia,
                                               INT *ja,
                                               REAL *aval,
                                               INT *nB,
                                               INT *mB,
                                               INT *nnzB,
                                               INT *ib,
                                               INT *jb,
                                               REAL *bval,
                                               INT *nC,
                                               INT *mC,
                                               INT *nnzC,
                                               INT *ic,
                                               INT *jc,
                                               REAL *cval,
                                               REAL *b,
                                               REAL *u)
{
    dBLCmat A; // coefficient matrix
    dCSRmat matA, matB, matBt, matC;
    dvector rhs, sol; // right-hand-side, solution
    precond_ns_param psparam; // parameters for ns precond
    precond_ns_data  psdata; // data for ns precond
    int i,flag;
    
    /** initialize dBLCmat **/
    A.brow = 2;
    A.bcol = 2;
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    fasp_mem_check((void *)A.blocks, "block matrix:cannot allocate memory!\n", ERROR_ALLOC_MEM);
    A.blocks[0] = &matA;
    A.blocks[1] = &matBt;
    A.blocks[2] = &matB;
    A.blocks[3] = &matC;
    
    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    //char *inputfile = "fasp4ns.dat";
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
    
    // initialize matrix
    matA.row = *nA; matA.col = *nA; matA.nnz = *nnzA;
    matA.IA  = ia;  matA.JA  = ja; matA.val = aval;
    matBt.row = *nB; matBt.col = *mB; matBt.nnz = *nnzB;
    matBt.IA  = ib;  matBt.JA  = jb; matBt.val = bval;
    matB.row = *nC; matB.col = *mC; matB.nnz = *nnzC;
    matB.IA  = ic;  matB.JA  = jc; matB.val = cval;
    fasp_dcsr_alloc(*nC,*nC,1,&matC);
    
    //fasp_check_dCSRmat(A.blocks[0]);
    //fasp_check_dCSRmat (&matB);
    //fasp_check_dCSRmat (&matBt);
    //printf("A = %d by %d\n", matA.row,matA.col);
    //printf("B = %d by %d\n", matBt.row,matBt.col);
    //printf("C = %d by %d\n", matB.row,matB.col);
    //  Shift the index to start from 0 (for C routines)
    //for (i=0;i<*nnzA;i++) matA.JA[i] --;
    //for (i=0;i<*nnzB;i++) matB.JA[i] --;
    //for (i=0;i<*nnzC;i++) matC.JA[i] --;
    
    //fasp_dcsr_trans(&matB, &matBt);
    
    rhs.row = *nA+*nC; rhs.val = b;
    sol.row = *nA+*nC; sol.val = u;
    
    if (print_level>8) {
        fasp_dcoo_write("A.dat",  &matA);
        fasp_dcoo_write("Bt.dat", &matBt);
        fasp_dcoo_write("B.dat",  &matB);
        fasp_dcoo_write("C.dat",  &matC);
        fasp_dvec_write("rhs.dat", &rhs);
    }
    
    //printf(" Finish Construct A\n");
    //fasp_mem_usage();
    
    if (print_level>0) {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }
    
    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam, &amgparam, &iluparam, &schparam);
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
