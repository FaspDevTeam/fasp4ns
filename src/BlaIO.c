/*! \file  BlaIO.c
 *
 *  \brief I/O functions for NS solvers
 *
 *  \note  This file contains Level-1 (Bla) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--Present by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  TODO: Need to simplify these functions!!! --Chensong
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_dblc_read (char *fileA11, char *fileA21, char *fileA22,
 *                          char *filerhs, dBLCmat *A, dvector *r)
 *
 * \brief Read the dBLCmat and right-hand side from files
 *
 * \param fileA11    file name of A in CSR format
 * \param fileA21    file name of B in CSR format
 * \param fileA22    file name of C in CSR format
 * \param filerhs    file name of right hand side
 * \param A          pointer to the dBLCmat matrix
 * \param r          pointer to the right-hand side vector
 *
 * \note
 *      E = (A B^T)
 *          (B C  )
 * File format:
 *   This routine reads three dCSRmat matrices from files!
 *
 * \author Lu WANG
 * \date   02/24/2012
 */
void fasp_dblc_read(char *fileA11,
                    char *fileA21,
                    char *fileA22,
                    char *filerhs,
                    dBLCmat *A,
                    dvector *r)
{
    int numA, nnz, numB, nnzb, nnzc;
    int i, ivalue, wall;
    double value;

    // read matrix A11
    FILE *fp = fopen(fileA11, "r");
    if (fp == NULL)
    {
        printf("### ERROR: Opening file %s failed!\n", fileA11);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s ...\n", __FUNCTION__, fileA11);

    wall = fscanf(fp, "%d %d", &numA, &nnz); // read dimension of the problem
    fasp_dcsr_alloc(numA, numA, nnz, A->blocks[0]);

    for (i = 0; i < numA + 1; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->IA[i] = ivalue;
    }
    for (i = 0; i < nnz; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->JA[i] = ivalue - 1;
    }
    for (i = 0; i < nnz; ++i)
    {
        wall = fscanf(fp, "%le", &value);
        A->blocks[0]->val[i] = value;
    }
    fclose(fp);

    // read matrix A21
    fp = fopen(fileA21, "r");
    if (fp == NULL)
    {
        printf("### ERROR: Opening file %s failed!\n", fileA21);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s ...\n", __FUNCTION__, fileA21);

    wall = fscanf(fp, "%d %d", &numB, &nnzb); // read dimension of the problem
    fasp_dcsr_alloc(numB, numA, nnzb, A->blocks[2]);

    for (i = 0; i < numB + 1; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->IA[i] = ivalue;
    }

    for (i = 0; i < nnzb; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->JA[i] = ivalue - 1;
    }

    for (i = 0; i < nnzb; ++i)
    {
        wall = fscanf(fp, "%le", &value);
        A->blocks[2]->val[i] = value;
    }
    fclose(fp);

    // get matrix A12
    fasp_dcsr_trans(A->blocks[2], A->blocks[1]);

    // read matrix A22
    fp = fopen(fileA22, "r");
    if (fp == NULL)
    {
        printf("### ERROR: Opening file %s failed!\n", fileA22);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s ...\n", __FUNCTION__, fileA22);

    wall = fscanf(fp, "%d %d", &numB, &nnzc); // read dimension of the problem
    fasp_dcsr_alloc(numB, numB, nnzc, A->blocks[3]);

    for (i = 0; i < numB + 1; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->IA[i] = ivalue;
    }

    for (i = 0; i < nnzc; ++i)
    {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->JA[i] = ivalue - 1;
    }

    for (i = 0; i < nnzc; ++i)
    {
        wall = fscanf(fp, "%le", &value);
        A->blocks[3]->val[i] = value;
    }
    fclose(fp);

    // read rhs
    fp = fopen(filerhs, "r");
    if (fp == NULL)
    {
        printf("### ERROR: Opening file %s failed!\n", filerhs);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s ...\n", __FUNCTION__, filerhs);

    fasp_dvec_alloc(numA + numB, r);
    for (i = 0; i < numA + numB; ++i)
    {
        wall = fscanf(fp, "%le", &value);
        r->val[i] = value;
    }
    if (wall < 0)
        printf("### ERROR: Something wrong when reading vaules!\n");
    fclose(fp);
}

/**
 * \fn void fasp_dblc_read3 (char *fileA11, char *fileA21, char *fileA22,
 *                           char *filerhs, dBLCmat *A, dvector *r)
 *
 * \brief Read A and rhs from files
 *
 * \param fileA11  file name of A in COO format
 * \param fileA21  file name of B in COO format
 * \param fileA22  file name of C in COO format
 * \param filerhs  file name of right-hand side
 * \param A        pointer to the dBLCmat
 * \param r        pointer to the right-hand side
 *
 * \note
 *      A = (A11 A12)
 *          (A21 A22)
 *
 * \author Chensong Zhang
 * \date   09/03/2022
 */
void fasp_dblc_read3(char *fileA11,
                     char *fileA21,
                     char *fileA22,
                     char *filerhs,
                     dBLCmat *A,
                     dvector *r)
{
    int num, i, wall;
    double value;

    // read matrices
    fasp_dcoo_read1(fileA11, A->blocks[0]);
    fasp_dcoo_read1(fileA21, A->blocks[2]);
    fasp_dcoo_read1(fileA22, A->blocks[3]);
    fasp_dcsr_trans(A->blocks[2], A->blocks[1]);

    num = A->blocks[0]->row + A->blocks[3]->row;

    // read rhs
    FILE *fp = fopen(filerhs, "r");
    if (fp == NULL)
    {
        printf("### ERROR: Opening file %s failed!\n", filerhs);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s ...\n", __FUNCTION__, filerhs);

    fasp_dvec_alloc(num, r);
    for (i = 0; i < num; ++i)
    {
        wall = fscanf(fp, "%le", &value);
        if (wall < 0)
        {
            printf("### ERROR: Something wrong when reading vaules!\n");
            exit(ERROR_READ_FILE);
        }
        r->val[i] = value;
    }

    fclose(fp);
}

/**
 * \fn void fasp_dblc_read_ruth (char *fileA11, char *fileA12, char *fileA21,
 *                               char *fileA22, char *filerhs, char *filex0,
 *                               dBLCmat *A, dvector *r, dvector *x0)
 *
 * \brief Read the dBLCmat matrix and rhs from files
 *
 * \param fileA11  file name of A
 * \param fileA12  file name of B
 * \param fileA21  file name of C
 * \param fileA22  file name of D
 * \param filerhs  file name of the right-hand side
 * \param filex0   file name of the initial guess
 * \param A        pointer to the dBLCmat matrix
 * \param r        pointer to the right-hand side vector
 * \param x0       pointer to the initial guess
 *
 * \note
 *      A = (A11 A12)
 *          (A21 A22)
 *
 * \author Lu WANG
 * \date   02/24/2012
 */
void fasp_dblc_read_ruth(char *fileA11,
                         char *fileA12,
                         char *fileA21,
                         char *fileA22,
                         char *filerhs,
                         char *filex0,
                         dBLCmat *A,
                         dvector *r,
                         dvector *x0)
{
    fasp_dcoo_read(fileA11, A->blocks[0]);
    fasp_dcoo_read(fileA12, A->blocks[1]);
    fasp_dcoo_read(fileA21, A->blocks[2]);
    fasp_dcoo_read(fileA22, A->blocks[3]);

    fasp_dvec_read(filerhs, r);
    fasp_dvec_read(filex0, x0);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
