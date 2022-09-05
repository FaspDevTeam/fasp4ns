/*! \file  SolWrapper.c
 *
 *  \brief Wrappers for accessing functions for advanced users
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         AuxInput.c, AuxParam.c, and SolNavierStokes.c
 *
 *  \note  IMPORTANT: The wrappers DO NOT change the original matrix data. Users
 *         should shift the matrix indices in order to make the IA and JA to start
 *         from 0 instead of 1.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2018 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_fwrapper_dblc_krylov_nstokes_ (INT *nA, INT *nnzA, INT *ia,
 *                                              INT *ja, REAL *aval, INT *nB,
 *                                              INT *mB, INT *nnzB, INT *ib,
 *                                              INT *jb, REAL *bval, INT *nC,
 *                                              INT *mC, INT *nnzC, INT *ic,
 *                                              INT *jc, REAL *cval,
 *                                              REAL *b, REAL *u)
 *
 * \brief Solve [A B; C O] u = b by Krylov method with block preconditioners
 *
 * \param nA       num of rows/cols of A
 * \param nnzA     num of nonzeros of A
 * \param ia       IA of A in CSR format
 * \param ja       JA of A in CSR format
 * \param aval     VAL of A in CSR format
 * \param nB       num of rows of B
 * \param mB       num of cols of B
 * \param nnzB     num of nonzeros of B
 * \param ib       IA of B in CSR format
 * \param jb       JA of B in CSR format
 * \param bval     VAL of B in CSR format
 * \param nC       num of rows of C
 * \param mC       num of cols of C
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
 * Modified by Chensong Zhang on 03/16/2018: Fixed format problem
 * Modified by Chensong Zhang on 03/26/2018: Change file name
 */
void fasp_fwrapper_dblc_krylov_nstokes_(INT *nA,
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
    dCSRmat matA11, matA21, matA12, matA22;
    dvector rhs, sol;         // right-hand-side, solution
    precond_ns_param psparam; // parameters for ns precond
    precond_ns_data psdata;   // data for ns precond
    int i, flag;

    char *inputfile = "ini/ns.dat";
    input_ns_param inparam;    // parameters from input files
    itsolver_ns_param itparam; // parameters for itsolver
    AMG_ns_param amgparam;     // parameters for AMG
    ILU_param iluparam;        // parameters for ILU
    SWZ_param swzparam;        // parameters for Schwarz

    /** Step 0. Read input parameters */
    fasp_ns_param_input(inputfile, &inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &swzparam);

    // Set local parameters
    const int print_level = inparam.print_level;
    const int problem_num = inparam.problem_num;
    const int itsolver_type = inparam.solver_type;
    const int precond_type = inparam.precond_type;

#if DEBUG_MODE > 0
    printf("### DEBUG: nA = %d\n", *nA);
    printf("### DEBUG: nB = %d, mB = %d\n", *nB, *mB);
    printf("### DEBUG: nC = %d, mc = %d\n", *nC, *mC);
#endif

    // initialize dBLCmat pointer
    A.brow = 2;
    A.bcol = 2;
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    if (A.blocks == NULL)
    {
        printf("### ERROR: Cannot allocate memory %s!\n", __FUNCTION__);
        exit(ERROR_ALLOC_MEM);
    }
    A.blocks[0] = &matA11;
    A.blocks[1] = &matA12;
    A.blocks[2] = &matA21;
    A.blocks[3] = &matA22;

    // initialize matrix
    matA11.row = *nA;
    matA11.col = *nA;
    matA11.nnz = *nnzA;
    matA11.IA = ia;
    matA11.JA = ja;
    matA11.val = aval;

    matA12.row = *nB;
    matA12.col = *mB;
    matA12.nnz = *nnzB;
    matA12.IA = ib;
    matA12.JA = jb;
    matA12.val = bval;

    matA21.row = *nC;
    matA21.col = *mC;
    matA21.nnz = *nnzC;
    matA21.IA = ic;
    matA21.JA = jc;
    matA21.val = cval;

    if (print_level > 9)
    {
        fasp_dcsr_write_coo("A11.coo", &matA11);
        fasp_dcsr_write_coo("A12.coo", &matA12);
        fasp_dcsr_write_coo("A21.coo", &matA21);
        fasp_dvec_write("rhs.vec", &rhs);
    }

    // generate an empty matrix
    fasp_dcsr_alloc(*nC, *nC, 1, &matA22);

    // initialize rhs and sol vectors
    rhs.row = *nA + *nC;
    rhs.val = b;
    sol.row = *nA + *nC;
    sol.val = u;

    if (print_level > 0)
    {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }

    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam,
                                                 &amgparam, &iluparam, &swzparam);

    if (print_level > 9)
    {
        fasp_dvec_write("sol.vec", &sol);
        printf("Press ENTER to continue...");
        getchar();
    }
}

/**
 * \fn void fasp_fwrapper_dblc_krylov_sstokes_ (INT *nA, INT *nnzA, INT *ia,
 *                                              INT *ja, REAL *aval, INT *nB,
 *                                              INT *nnzB, INT *ib, INT *jb,
 *                                              REAL *bval, INT *nC, INT *nnzC,
 *                                              INT *ic, INT *jc, REAL *cval,
 *                                              REAL *b, REAL *u)
 *
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
 * Modified by Chensong Zhang on 03/16/2018: Fixed format problem
 * Modified by Chensong Zhang on 03/26/2018: Change file name
 */
void fasp_fwrapper_dblc_krylov_sstokes_(INT *nA,
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
    dCSRmat matA11, matA21, matA12, matA22;
    dvector rhs, sol;         // right-hand-side, solution
    precond_ns_param psparam; // parameters for ns precond
    precond_ns_data psdata;   // data for ns precond
    int i, flag;

    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    input_ns_param inparam;    // parameters from input files
    itsolver_ns_param itparam; // parameters for itsolver
    AMG_ns_param amgparam;     // parameters for AMG
    ILU_param iluparam;        // parameters for ILU
    SWZ_param swzparam;        // parameters for Schwarz

    fasp_ns_param_input(inputfile, &inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &swzparam);

    // set local parameters
    const int print_level = inparam.print_level;
    const int problem_num = inparam.problem_num;
    const int itsolver_type = inparam.solver_type;
    const int precond_type = inparam.precond_type;

#if DEBUG_MODE > 0
    printf("### DEBUG: nA = %d, nB = %d, nC = %d\n", *nA, *nB, *nC);
#endif

    // initialize dBLCmat pointer
    A.brow = 2;
    A.bcol = 2;
    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    if (A.blocks == NULL)
    {
        printf("### ERROR: Cannot allocate memory %s!\n", __FUNCTION__);
        exit(ERROR_ALLOC_MEM);
    }
    A.blocks[0] = &matA11;
    A.blocks[1] = &matA12;
    A.blocks[2] = &matA21;
    A.blocks[3] = &matA22;

    // initialize matrix
    matA11.row = *nA;
    matA11.col = *nA;
    matA11.nnz = *nnzA;
    matA11.IA = ia;
    matA11.JA = ja;
    matA11.val = aval;

    matA21.row = *nB;
    matA21.col = *nA;
    matA21.nnz = *nnzB;
    matA21.IA = ib;
    matA21.JA = jb;
    matA21.val = bval;

    matA22.row = *nC;
    matA22.col = *nC;
    matA22.nnz = *nnzC;
    matA22.IA = ic;
    matA22.JA = jc;
    matA22.val = cval;

    if (print_level > 9)
    {
        fasp_dcsr_write_coo("A11.coo", &matA11);
        fasp_dcsr_write_coo("A21.coo", &matA21);
        fasp_dcsr_write_coo("A22.coo", &matA22);
        fasp_dvec_write("rhs.vec", &rhs);
    }

    // get transform of B
    fasp_dcsr_trans(&matA21, &matA12);

    rhs.row = *nA + *nB;
    rhs.val = b;
    sol.row = *nA + *nB;
    sol.val = u;

    if (print_level > 0)
    {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }

    flag = fasp_solver_dblc_krylov_navier_stokes(&A, &rhs, &sol, &itparam,
                                                 &amgparam, &iluparam, &swzparam);

    if (print_level > 9)
    {
        fasp_dvec_write("sol.vec", &sol);
        printf("Press ENTER to continue...");
        getchar();
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
