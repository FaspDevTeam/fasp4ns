/*! \file  BlaIO.c
 *
 *  \brief I/O functions for NS solvers
 *
 *  \note  This file contains Level-1 (Bla) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  // TODO: Fix Doxygen. --Chensong
 *  // TODO: Remove unused functions. --Chensong
 */

#include "fasp.h"
#include "fasp_functs.h"
#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_dblc_read (char *fileA, char *fileB, char *fileC,
 *                          char *filerhs, dBLCmat *A, dvector *r)
 * \brief Read E and rhs from file in block_dSTRmat format
 *
 * \param fileA      file name of A
 * \param fileB      file name of B
 * \param fileC      file name of C
 * \param fileArhs   file name of right hand side
 * \param A          pointer to the dBLCmat
 *
 * \note
 * E = (A B^T)
 *     (B C)
 * File format:
 *   This routine reads a dCSRmat matrix from files in the following format:
 *
 * \author Lu WANG
 * \date   02/24/2012
 */
void fasp_dblc_read (char *fileA,
                     char *fileB,
                     char *fileC,
                     char *filerhs,
                     dBLCmat *A,
                     dvector *r)
{
    int numA,nnz,numB,nnzb,nnzc;
    int i, k, n;
    int ivalue,wall;
    double value;
    
    // read file A
    FILE *fp=fopen(fileA,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n", fileA);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileA);
    
    wall = fscanf(fp,"%d %d",&numA,&nnz); // read dimension of the problem
    fasp_dcsr_alloc (numA,numA,nnz,A->blocks[0]);
    // read matrix A
    for (i=0;i<numA+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->IA[i]=ivalue;
    }
    for (i=0;i<nnz;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->JA[i]=ivalue-1;
    }
    for (i=0;i<nnz;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[0]->val[i]=value;
    }
    fclose(fp);
    
    fp=fopen(fileB,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",fileB);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileB);
    
    wall = fscanf(fp,"%d %d",&numB,&nnzb); // read dimension of the problem
    fasp_dcsr_alloc (numB,numA,nnzb,A->blocks[2]);
    
    // read matrix B
    for (i=0;i<numB+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->IA[i]=ivalue;
    }
    
    for (i=0;i<nnzb;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->JA[i]=ivalue-1;
    }
    
    for (i=0;i<nnzb;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[2]->val[i]=value;
    }
    fclose(fp);
    fasp_dcsr_trans(A->blocks[2],A->blocks[1]);
    
    fp=fopen(fileC,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",fileC);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileC);
    
    wall = fscanf(fp,"%d %d",&numB,&nnzc); // read dimension of the problem
    fasp_dcsr_alloc (numB,numB,nnzc,A->blocks[3]);
    
    // read matrix B
    for (i=0;i<numB+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->IA[i]=ivalue;
    }
    
    for (i=0;i<nnzc;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->JA[i]=ivalue-1;
    }
    
    for (i=0;i<nnzc;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[3]->val[i]=value;
    }
    fclose(fp);
    
    fp=fopen(filerhs,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",filerhs);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, filerhs);
    
    fasp_dvec_alloc (numA+numB,r);
    for (i=0;i<numA+numB;++i) {
        wall = fscanf(fp, "%le", &value);
        r->val[i]=value;
    }
    fclose(fp);
}

/**
 * \fn void fasp_dblc_read1 (char *fileA, char *fileB, char *fileC, char *filerhs,
 *                           dBLCmat *A, dvector *r)
 * \brief Read E and rhs from file in block_dSTRmat format
 *
 * \param fileA      file name of A
 * \param fileB      file name of B
 * \param fileC      file name of C
 * \param filerhs    file name of right-hand side
 * \param A          pointer to the dBLCmat
 * \param r          pointer to the right-hand side
 * \note
 * E = (A B^T)
 *     (B C)
 * File format:
 *   This routine reads a dCSRmat matrix from files in the following format:
 *
 * \author Xiaozhe Hu
 * \date   11/01/2013
 */
void fasp_dblc_read_1 (char *fileA,
                       char *fileB,
                       char *fileC,
                       char *filerhs,
                       dBLCmat *A,
                       dvector *r)
{
    int numA,nnz,numB,nnzb,nnzc;
    int i, k, n;
    int ivalue,wall;
    double value;
    
    // read file A
    FILE *fp=fopen(fileA,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",fileA);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileA);
    
    wall = fscanf(fp,"%d %d",&numA,&nnz); // read dimension of the problem
    fasp_dcsr_alloc (numA,numA,nnz,A->blocks[0]);
    // read matrix A
    for (i=0;i<numA+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->IA[i]=ivalue;
    }
    for (i=0;i<nnz;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[0]->JA[i]=ivalue;
    }
    for (i=0;i<nnz;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[0]->val[i]=value;
    }
    fclose(fp);
    
    fp=fopen(fileB,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",fileB);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileB);
    
    wall = fscanf(fp,"%d %d",&numB,&nnzb); // read dimension of the problem
    fasp_dcsr_alloc (numB,numA,nnzb,A->blocks[2]);
    
    // read matrix B
    for (i=0;i<numB+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->IA[i]=ivalue;
    }
    
    for (i=0;i<nnzb;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[2]->JA[i]=ivalue;
    }
    
    for (i=0;i<nnzb;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[2]->val[i]=value;
    }
    fclose(fp);
    fasp_dcsr_trans(A->blocks[2],A->blocks[1]);
    
    fp=fopen(fileC,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",fileC);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, fileC);
    
    wall = fscanf(fp,"%d %d",&numB,&nnzc); // read dimension of the problem
    fasp_dcsr_alloc (numB,numB,nnzc,A->blocks[3]);
    
    // read matrix B
    for (i=0;i<numB+1;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->IA[i]=ivalue;
    }
    
    for (i=0;i<nnzc;++i) {
        wall = fscanf(fp, "%d", &ivalue);
        A->blocks[3]->JA[i]=ivalue;
    }
    
    for (i=0;i<nnzc;++i) {
        wall = fscanf(fp, "%le", &value);
        A->blocks[3]->val[i]=value;
    }
    fclose(fp);
    
    fp=fopen(filerhs,"r");
    if ( fp == NULL ) {
        printf("### ERROR: Opening file %s failed!\n",filerhs);
        exit(ERROR_OPEN_FILE);
    }
    printf("%s: reading file %s...\n", __FUNCTION__, filerhs);
    
    fasp_dvec_alloc (numA+numB,r);
    for (i=0;i<numA+numB;++i) {
        wall = fscanf(fp, "%le", &value);
        r->val[i]=value;
    }
    fclose(fp);
}

/**
 * \fn void fasp_dblc_read_ruth (char *fileA, char *fileB, char *fileC, 
 *                               char *fileD, char *filerhs, char *filex0, 
 *                               dBLCmat *A, dvector *r, dvector *x0)
 *
 * \brief Read E and rhs from file in block_dSTRmat format
 *
 * \param fileA      file name of A
 * \param fileB      file name of B
 * \param fileC      file name of C
 * \param fileArhs   file name of right hand side
 * \param A          pointer to the dBLCmat
 *
 * \note
 * E = (A B^T)
 *     (B C)
 * File format:
 *   This routine reads a dCSRmat matrix from files in the following format:
 *
 * \author Lu WANG
 * \date   02/24/2012
 */
void fasp_dblc_read_ruth (char *fileA,
                          char *fileB,
                          char *fileC,
                          char *fileD,
                          char *filerhs,
                          char *filex0,
                          dBLCmat *A,
                          dvector *r,
                          dvector *x0)
{
    fasp_dcoo_read (fileA,A->blocks[0]);
    fasp_dcoo_read (fileB,A->blocks[1]);
    fasp_dcoo_read (fileC,A->blocks[2]);
    fasp_dcoo_read (fileD,A->blocks[3]);
    fasp_dvec_read (filerhs,r);
    fasp_dvec_read (filex0,x0);
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
