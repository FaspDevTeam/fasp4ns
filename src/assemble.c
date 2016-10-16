/*! \file  assemble.c
 *  \brief Assembling for the N-S equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

//#include "basisP2.inl"
//#include "basisP0.inl"
//#include "assemble_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int assemble (block_dCSRmat *ptr_A, dvector *ptr_b, int gridnum)
 * \brief Assemble stiffness matrix *ptr_A and righ-hand side *ptr_b for uniform 
 *        grid (gridnum by gridnum).
 *
 * \param *ptr_A    pointer to stiffness matrix
 * \param *ptr_b    pointer to right hand side
 * \param gridnum   number of node per edges
 * \return          1 if succeed 0 if fail
 *
 * \author Lu Wang
 * \date 12/1/2011
 * 
 * TODO: Complete this code!!! --Chensong
 */
int assemble (block_dCSRmat *ptr_A,
              dvector *ptr_b,
              int gridnum)
{
	ddenmat nodes; // the first column stores the x coordinate of points
                   // the second column stores the y coordinate of points
	idenmat elem;  // triangulation: store 3 points and 3 edges corresponding
                   // to the element in each row
	idenmat bdnodes; // boundary node indeces
	idenmat bdedges; // boundary edge indeces
	int i,j,k,l;
	
	// get coarse grid from a disk file
	int IsExist = getGridInfo(&nodes, &elem, &bdnodes, &bdedges,gridnum);
	if (IsExist==0) {
		printf("### ERROR: Constructing grid fails!\n");
		return ERROR_UNKNOWN;
	}
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
