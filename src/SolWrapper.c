/*! \file  SolWrapper.c
 *
 *  \brief Wrappers for accessing functions for advanced users
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         AuxInput.c, AuxParam.c, and SolNavierStokes.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  // TODO: Fix Doxygen. --Chensong
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_fwrapper_krylov_navier_stokes_sym_ (INT *nA, INT *nnzA, INT *ia,
 *                                                   INT *ja, REAL *aval, INT *nB,
 *                                                   INT *nnzB, INT *ib, INT *jb,
 *                                                   REAL *bval, INT *nC, INT *nnzC,
 *                                                   INT *ic, INT *jc, REAL *cval,
 *                                                   REAL *b, REAL *u)
 * \brief Solve [A B'; B C] u = b by Krylov method with block preconditioners
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
 * \param nC       num of cols of C
 * \param nnzC     num of nonzeros of C
 * \param ic       IA of C in CSR format
 * \param jc       JA of C in CSR format
 * \param cval     VAL of C in CSR format
 * \param b        rhs vector
 * \param u        solution vector
 *
 * \author Lu Wang
 * \date   03/14/2012
 *
 * Modified by Chensong Zhang on 03/13/2018
 */
void fasp_fwrapper_krylov_navier_stokes_sym_ (INT *nA,
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
    int i, flag;
    
    /** initialize dBLCmat **/
    A.brow = 2;
    A.bcol = 2;
    
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    if ( A.blocks == NULL ) {
        printf("### ERROR: Cannot allocate memory %s!\n", __FUNCTION__);
        exit(ERROR_ALLOC_MEM);
    }
    A.blocks[0] = &matA;
    A.blocks[1] = &matBt;
    A.blocks[2] = &matB;
    A.blocks[3] = &matC;
    
    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    input_ns_param     inparam;  // parameters from input files
    itsolver_ns_param  itparam;  // parameters for itsolver
    AMG_ns_param       amgparam; // parameters for AMG
    ILU_param          iluparam; // parameters for ILU
    SWZ_param          swzparam; // parameters for Schwarz
    
    fasp_ns_param_input(inputfile, &inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &swzparam);
    
    // set local parameters
    const int print_level   = inparam.print_level;
    const int problem_num   = inparam.problem_num;
    const int itsolver_type = inparam.solver_type;
    const int precond_type  = inparam.precond_type;
    
    // initialize matrix
    matA.row = *nA; matA.col = *nA; matA.nnz = *nnzA;
    matA.IA  = ia;  matA.JA  = ja;  matA.val = aval;
    matB.row = *nB; matB.col = *nA; matB.nnz = *nnzB;
    matB.IA  = ib;  matB.JA  = jb;  matB.val = bval;
    matC.row = *nC; matC.col = *nC; matC.nnz = *nnzC;
    matC.IA  = ic;  matC.JA  = jc;  matC.val = cval;
    
    // shift the index to start from 0 (for C routines)
    for (i=0;i<*nnzA;i++) matA.JA[i] --;
    for (i=0;i<*nnzB;i++) matB.JA[i] --;
    for (i=0;i<*nnzC;i++) matC.JA[i] --;

    // get transform of B
    fasp_dcsr_trans(&matB, &matBt);
    
    rhs.row = *nA + *nB; rhs.val = b;
    sol.row = *nA + *nB; sol.val = u;

    if (print_level>0) {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }
    
    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam,
                                                 &amgparam, &iluparam, &swzparam);
}

/**
 * \fn void fasp_fwrapper_krylov_navier_stokes_nsym_ (INT *nA, INT *nnzA, INT *ia,
 *                                                    INT *ja, REAL *aval, INT *nB,
 *                                                    INT *mB, INT *nnzB, INT *ib,
 *                                                    INT *jb, REAL *bval, INT *nC,
 *                                                    INT *mC, INT *nnzC, INT *ic,
 *                                                    INT *jc, REAL *cval, REAL *b,
 *                                                    REAL *u)
 * \brief Solve [A B; C O] u = b by Krylov method with block preconditioners
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
 * \date   03/20/2014
 *
 * Modified by Chensong Zhang on 03/13/2018
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
    if ( A.blocks == NULL ) {
        printf("### ERROR: Cannot allocate memory %s!\n", __FUNCTION__);
        exit(ERROR_ALLOC_MEM);
    }
    A.blocks[0] = &matA;
    A.blocks[1] = &matBt;
    A.blocks[2] = &matB;
    A.blocks[3] = &matC;
    
    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    input_ns_param     inparam;  // parameters from input files
    itsolver_ns_param  itparam;  // parameters for itsolver
    AMG_ns_param       amgparam; // parameters for AMG
    ILU_param          iluparam; // parameters for ILU
    SWZ_param          swzparam; // parameters for Schwarz
    
    fasp_ns_param_input(inputfile, &inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &swzparam);
    
    // Set local parameters
    const int print_level   = inparam.print_level;
    const int problem_num   = inparam.problem_num;
    const int itsolver_type = inparam.solver_type;
    const int precond_type  = inparam.precond_type;
    
    // initialize matrix
    matA.row  = *nA; matA.col  = *nA; matA.nnz  = *nnzA;
    matA.IA   = ia;  matA.JA   = ja;  matA.val  = aval;
    matBt.row = *nB; matBt.col = *mB; matBt.nnz = *nnzB;
    matBt.IA  = ib;  matBt.JA  = jb;  matBt.val = bval;
    matB.row  = *nC; matB.col  = *mC; matB.nnz  = *nnzC;
    matB.IA   = ic;  matB.JA   = jc;  matB.val  = cval;

    // generate an empty matrix
    fasp_dcsr_alloc(*nC,*nC,1,&matC);

    // shift the index to start from 0 (for C routines)
    for (i=0;i<*nnzA;i++) matA.JA[i] --;
    for (i=0;i<*nnzB;i++) matB.JA[i] --;
    for (i=0;i<*nnzC;i++) matC.JA[i] --;

    rhs.row = *nA+*nC; rhs.val = b;
    sol.row = *nA+*nC; sol.val = u;

    if (print_level>0) {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }
    
    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam,
                                                 &amgparam, &iluparam, &swzparam);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
