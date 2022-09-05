/*! \file NavierStokes.c
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
 * \author Lu Wang
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 * Modified by Chensong Zhang on 03/13/2018
 */
int main(int argc, const char *argv[])
{
    dBLCmat A;
    dCSRmat RR, RW, WR, WW, Acsr, Mp;
    dvector b, uh, bcsr;
    ivector u_idx, p_idx;
    int i, flag = 0;

    /** initialize dBLCmat **/
    A.brow = 2;
    A.bcol = 2;

    A.blocks = (dCSRmat **)calloc(4, sizeof(dCSRmat *));
    if (A.blocks == NULL)
    {
        printf("### ERROR: Cannot allocate memory %s!\n", __FUNCTION__);
        return ERROR_ALLOC_MEM;
    }
    A.blocks[0] = &RR;
    A.blocks[1] = &RW;
    A.blocks[2] = &WR;
    A.blocks[3] = &WW;

    /** Step 0. Read input parameters */
    char *inputfile = "ini/ns.dat";
    input_ns_param inparam;    // parameters from input files
    itsolver_ns_param itparam; // parameters for itsolver
    AMG_ns_param amgparam;     // parameters for AMG
    ILU_param iluparam;        // parameters for ILU
    SWZ_param schparam;        // parameters for Schwarz

    fasp_ns_param_input(inputfile, &inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &schparam);

    // Set local parameters
    const int print_level = inparam.print_level;
    const int problem_num = inparam.problem_num;
    const int itsolver_type = inparam.solver_type;
    const int precond_type = inparam.precond_type;

    printf("Test Problem %d\n", problem_num);

    /** Step 1. Assemble matrix and right-hand side */
    char filename1[512]; // *datafile1;
    char filename2[512]; // *datafile2;
    char filename3[512]; // *datafile3;

    strncpy(filename1, inparam.workdir, 128);
    strncpy(filename2, inparam.workdir, 128);
    strncpy(filename3, inparam.workdir, 128);

    /* Assemble A and b. P1 FE discretization for Poisson. */
    if (problem_num == 1)
    {
        char *fileA11 = "data/test_1/Matrix_A";
        char *fileA21 = "data/test_1/Matrix_B";
        char *fileA22 = "data/test_1/Matrix_C";
        char *fileRhs = "data/test_1/RHS";
        fasp_dblc_read(fileA11, fileA21, fileA22, fileRhs, &A, &b);
    }

    // --------------------------- //
    //  Examples wih pressure mass
    // --------------------------- //

    else if (problem_num == 100)
    {
        char *fileA11 = "data/Stokes-P2P0/1/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/1/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/1/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/1/b.dat";
        char *fileMPr = "data/Stokes-P2P0/1/Mp.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);
        // fasp_dcoo_read(fileA22, A.blocks[3]);
        A.blocks[3] = NULL;

        fasp_dcoo_read(fileMPr, &Mp);

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 101)
    {

        char *fileA11 = "data/Stokes-P2P0/2/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/2/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/2/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/2/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);
        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 102)
    {
        char *fileA11 = "data/Stokes-P2P0/3/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/3/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/3/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/3/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);
        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 103)
    {

        char *fileA11 = "data/Stokes-P2P0/4/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/4/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/4/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/4/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 104)
    {
        char *fileA11 = "data/Stokes-P2P0/5/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/5/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/5/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/5/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 105)
    {
        char *fileA11 = "data/Stokes-P2P0/6/A11.dat";
        char *fileA21 = "data/Stokes-P2P0/6/A12.dat";
        char *fileA12 = "data/Stokes-P2P0/6/A21.dat";
        char *fileRhs = "data/Stokes-P2P0/6/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 201)
    {
        char *fileA11 = "data/Stokes-P2P1/1/A11.dat";
        char *fileA21 = "data/Stokes-P2P1/1/A12.dat";
        char *fileA12 = "data/Stokes-P2P1/1/A21.dat";
        char *fileRhs = "data/Stokes-P2P1/1/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 202)
    {
        char *fileA11 = "data/Stokes-P2P1/2/A11.dat";
        char *fileA21 = "data/Stokes-P2P1/2/A12.dat";
        char *fileA12 = "data/Stokes-P2P1/2/A21.dat";
        char *fileRhs = "data/Stokes-P2P1/2/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 203)
    {
        char *fileA11 = "data/Stokes-P2P1/3/A11.dat";
        char *fileA21 = "data/Stokes-P2P1/3/A12.dat";
        char *fileA12 = "data/Stokes-P2P1/3/A21.dat";
        char *fileRhs = "data/Stokes-P2P1/3/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 204)
    {
        char *fileA11 = "data/Stokes-P2P1/4/A11.dat";
        char *fileA21 = "data/Stokes-P2P1/4/A12.dat";
        char *fileA12 = "data/Stokes-P2P1/4/A21.dat";
        char *fileRhs = "data/Stokes-P2P1/4/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 205)
    {

        char *fileA11 = "data/Stokes-P2P1/5/A11.dat";
        char *fileA21 = "data/Stokes-P2P1/5/A12.dat";
        char *fileA12 = "data/Stokes-P2P1/5/A21.dat";
        char *fileRhs = "data/Stokes-P2P1/5/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 301)
    {
        char *fileA11 = "data/Stokes-RTP0/1/A11.dat";
        char *fileA21 = "data/Stokes-RTP0/1/A12.dat";
        char *fileA12 = "data/Stokes-RTP0/1/A21.dat";
        char *fileRhs = "data/Stokes-RTP0/1/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 302)
    {
        char *fileA11 = "data/Stokes-RTP0/2/A11.dat";
        char *fileA21 = "data/Stokes-RTP0/2/A12.dat";
        char *fileA12 = "data/Stokes-RTP0/2/A21.dat";
        char *fileRhs = "data/Stokes-RTP0/2/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 303)
    {
        char *fileA11 = "data/Stokes-RTP0/3/A11.dat";
        char *fileA21 = "data/Stokes-RTP0/3/A12.dat";
        char *fileA12 = "data/Stokes-RTP0/3/A21.dat";
        char *fileRhs = "data/Stokes-RTP0/3/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 304)
    {
        char *fileA11 = "data/Stokes-RTP0/4/A11.dat";
        char *fileA21 = "data/Stokes-RTP0/4/A12.dat";
        char *fileA12 = "data/Stokes-RTP0/4/A21.dat";
        char *fileRhs = "data/Stokes-RTP0/4/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 305)
    {
        char *fileA11 = "data/Stokes-RTP0/5/A11.dat";
        char *fileA21 = "data/Stokes-RTP0/5/A12.dat";
        char *fileA12 = "data/Stokes-RTP0/5/A21.dat";
        char *fileRhs = "data/Stokes-RTP0/5/b.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA21, A.blocks[1]);
        fasp_dcoo_read(fileA12, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 401)
    {
        // read in data
        char *fileA11 = "data/P2P0-no-divdiv/set-1/A.dat";
        char *fileRhs = "data/P2P0-no-divdiv/set-1/b.dat";
        char *fileMPr = "data/P2P0-no-divdiv/set-1/Mp.dat";

        fasp_dcoo_shift_read(fileA11, &Acsr);
        fasp_dvec_read(fileRhs, &bcsr);
        fasp_dcoo_shift_read(fileMPr, &Mp);

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

        for (i = 0; i < u_idx.row; i++)
            b.val[i] = bcsr.val[u_idx.val[i]];
        for (i = u_idx.row; i < bcsr.row; i++)
            b.val[i] = bcsr.val[p_idx.val[i - u_idx.row]];
    }

    else if (problem_num == 402)
    {
        // read in data
        char *fileA11 = "data/P2P0-no-divdiv/set-2/A.dat";
        char *fileRhs = "data/P2P0-no-divdiv/set-2/b.dat";
        char *fileMPr = "data/P2P0-no-divdiv/set-2/Mp.dat";

        fasp_dcoo_shift_read(fileA11, &Acsr);
        fasp_dvec_read(fileRhs, &bcsr);
        fasp_dcoo_shift_read(fileMPr, &Mp);

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

        for (i = 0; i < u_idx.row; i++)
            b.val[i] = bcsr.val[u_idx.val[i]];
        for (i = u_idx.row; i < bcsr.row; i++)
            b.val[i] = bcsr.val[p_idx.val[i - u_idx.row]];
    }

    else if (problem_num == 403)
    {
        // read in data
        char *fileA11 = "data/P2P0-no-divdiv/set-11/A.dat";
        char *fileRhs = "data/P2P0-no-divdiv/set-11/b.dat";
        char *fileMPr = "data/P2P0-no-divdiv/set-11/Mp.dat";

        fasp_dcoo_shift_read(fileA11, &Acsr);
        fasp_dvec_read(fileRhs, &bcsr);
        fasp_dcoo_shift_read(fileMPr, &Mp);

        // read in index for velocity and pressure
        char *file_u_idx = "data/P2P0-no-divdiv/set-11/ISu.dat";
        char *file_p_idx = "data/P2P0-no-divdiv/set-11/ISp.dat";

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

        for (i = 0; i < u_idx.row; i++)
            b.val[i] = bcsr.val[u_idx.val[i]];
        for (i = u_idx.row; i < bcsr.row; i++)
            b.val[i] = bcsr.val[p_idx.val[i - u_idx.row]];
    }

    else if (problem_num == 501)
    {
        char *fileA11 = "data/P2P0-with-divdiv/set-1/Auu.dat";
        char *fileA21 = "data/P2P0-with-divdiv/set-1/Apu.dat";
        char *fileA12 = "data/P2P0-with-divdiv/set-1/Aup.dat";
        char *fileA22 = "data/P2P0-with-divdiv/set-1/App.dat";
        char *fileRhs = "data/P2P0-with-divdiv/set-1/b.dat";
        char *fileMPr = "data/P2P0-with-divdiv/set-1/Mp.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);
        fasp_dcoo_read(fileA22, A.blocks[3]);

        fasp_dcoo_read(fileMPr, &Mp);
        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 502)
    {
        char *fileA11 = "data/P2P0-with-divdiv/set-2/Auu.dat";
        char *fileA21 = "data/P2P0-with-divdiv/set-2/Apu.dat";
        char *fileA12 = "data/P2P0-with-divdiv/set-2/Aup.dat";
        char *fileA22 = "data/P2P0-with-divdiv/set-2/App.dat";
        char *fileRhs = "data/P2P0-with-divdiv/set-2/b.dat";
        char *fileMPr = "data/P2P0-with-divdiv/set-2/Mp.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);
        fasp_dcoo_read(fileA22, A.blocks[3]);

        fasp_dcoo_read(fileMPr, &Mp);
        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 601)
    {
        char *fileA11 = "data/FSIMatrices/Matrix_A";
        char *fileA21 = "data/FSIMatrices/Matrix_B";
        char *fileA22 = "data/FSIMatrices/Matrix_C";
        char *fileRhs = "data/FSIMatrices/RHS";
        fasp_dblc_read(fileA11, fileA21, fileA22, fileRhs, &A, &b);
    }

    else if (problem_num == 701)
    {
        char *fileA11 = "data/saed/W10H025-NI50NJ50/A.dat";
        char *fileA21 = "data/saed/W10H025-NI50NJ50/B.dat";
        char *fileA12 = "data/saed/W10H025-NI50NJ50/Bt.dat";
        char *fileRhs = "data/saed/W10H025-NI50NJ50/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 702)
    {
        char *fileA11 = "data/saed/W10H025-NI50NJ50-C/A.dat";
        char *fileA21 = "data/saed/W10H025-NI50NJ50-C/B.dat";
        char *fileA12 = "data/saed/W10H025-NI50NJ50-C/Bt.dat";
        char *fileRhs = "data/saed/W10H025-NI50NJ50-C/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 703)
    {
        char *fileA11 = "data/saed/W100H025-NI50NJ50-C/A.dat";
        char *fileA21 = "data/saed/W100H025-NI50NJ50-C/B.dat";
        char *fileA12 = "data/saed/W100H025-NI50NJ50-C/Bt.dat";
        char *fileRhs = "data/saed/W100H025-NI50NJ50-C/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 704)
    {
        char *fileA11 = "data/saed/W100H025-NI100NJ100-C/A.dat";
        char *fileA21 = "data/saed/W100H025-NI100NJ100-C/B.dat";
        char *fileA12 = "data/saed/W100H025-NI100NJ100-C/Bt.dat";
        char *fileRhs = "data/saed/W100H025-NI100NJ100-C/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 705)
    {
        char *fileA11 = "data/saed/W200H025-NI100NJ100-C/A.dat";
        char *fileA21 = "data/saed/W200H025-NI100NJ100-C/B.dat";
        char *fileA12 = "data/saed/W200H025-NI100NJ100-C/Bt.dat";
        char *fileRhs = "data/saed/W200H025-NI100NJ100-C/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);

        A.blocks[3] = NULL;

        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 801)
    {
        char *fileA11 = "data/pnp-stokes/A.dat";
        char *fileA21 = "data/pnp-stokes/B.dat";
        char *fileA12 = "data/pnp-stokes/Bt.dat";
        char *fileA22 = "data/pnp-stokes/C.dat";
        char *fileRhs = "data/pnp-stokes/rhs.dat";

        fasp_dcoo_read(fileA11, A.blocks[0]);
        fasp_dcoo_read(fileA12, A.blocks[1]);
        fasp_dcoo_read(fileA21, A.blocks[2]);
        fasp_dcoo_read(fileA22, A.blocks[3]);
        fasp_dvec_read(fileRhs, &b);
    }

    else if (problem_num == 901)
    {
        char *fileA11 = "../data/DLD/incore/A11.coo";
        char *fileA21 = "../data/DLD/incore/A21.coo";
        char *fileA22 = "../data/DLD/incore/A22.coo";
        char *fileRhs = "../data/DLD/incore/rhs.vec";
        fasp_dblc_read3(fileA11, fileA21, fileA22, fileRhs, &A, &b);
    }

    else if (problem_num == 902)
    {
        char *fileA11 = "../data/DLD/outcore/A11.coo";
        char *fileA21 = "../data/DLD/outcore/A21.coo";
        char *fileA22 = "../data/DLD/outcore/A22.coo";
        char *fileRhs = "../data/DLD/outcore/rhs.vec";
        fasp_dblc_read3(fileA11, fileA21, fileA22, fileRhs, &A, &b);
    }

    else
    {
        printf("### ERROR: Unrecognised problem number %d\n", problem_num);
        return ERROR_INPUT_PAR;
    }

    /** Step 2. Check matrix properties */

    /** Step 3. Select a solver */

    /** Step 4. Solve the system */
    fasp_mem_usage();

    if (print_level > 0)
    {
        printf("Max it num = %d\n", inparam.itsolver_maxit);
        printf("Tolerance  = %e\n", inparam.itsolver_tol);
    }

    // initial guess
    fasp_dvec_alloc(b.row, &uh);
    fasp_dvec_set(b.row, &uh, 0.0);

    if (itsolver_type == 0)
    {
        dvector b_test, uh_test;
        int nrow = A.blocks[0]->row;
        fasp_dvec_alloc(nrow, &b_test);
        fasp_dvec_set(nrow, &b_test, 0.0);
        for (i = 0; i < nrow; i++)
            b_test.val[i] = b.val[i];
        fasp_dvec_alloc(nrow, &uh_test);
        fasp_dvec_set(nrow, &uh_test, 0.0);
        fasp_solver_amg(A.blocks[0], &b_test, &uh_test, &(amgparam.param_v));
    }
    else
    {
        /*
        if (problem_num > 100 && problem_num < 100 ){
            flag = fasp_solver_dblc_krylov_navier_stokes_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        /*
        else if ( problem_num > 200 && problem_num < 300 ) {
            flag = fasp_solver_dblc_krylov_navier_stokes_schur_complement_with_pressure_mass(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam, &Mp);
        }
         */
        // else {
        flag = fasp_solver_dblc_krylov_navier_stokes(&A, &b, &uh, &itparam, &amgparam, &iluparam, &schparam);
        //}
    }

    fasp_dvec_write("solu.vec", &uh);

    fasp_dcsr_free(&RR);
    fasp_dcsr_free(&RW);
    fasp_dcsr_free(&WR);
    fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    free(A.blocks);

    return flag;
}
