/*! \file init.c
 *  \brief init useful parameters
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_ns_param_init (precond_ns_param *sparam)
 * \brief Initialize itsolver_param
 * \param Input: pointer to itsolver_param
 * \date 02/29/2012 
 */
void fasp_precond_ns_param_init (precond_ns_param *sparam)
{
    sparam->AMG_type = CLASSIC_AMG;
    sparam->print_level = 0;
    sparam->max_levels = 20;
}

/**
 * \fn void fasp_param_solver_init (itsolver_param *pdata)
 * \brief Initialize itsolver_param
 * \param Input: pointer to itsolver_param
 * \date 2010/03/23 
 */
void fasp_precond_ns_data_init (precond_ns_data *sprecdata)
{
	sprecdata->colA = 0; 
    sprecdata->colB = 0;
    sprecdata->col  = 0;
	
	sprecdata->beta = 1.0;
	
	sprecdata->mgl_data    = NULL;
    sprecdata->mgl_data_p  = NULL;
	sprecdata->print_level = 0;
	sprecdata->max_levels  = 20;
	sprecdata->maxit    = 0;
	sprecdata->amg_tol     = 1e-8;
	sprecdata->cycle_type  = V_CYCLE;
	sprecdata->smoother    = SMOOTHER_GS;
	sprecdata->presmooth_iter  = 1;
	sprecdata->postsmooth_iter = 1;
	sprecdata->coarsening_type = 1;
	sprecdata->relaxation      = 1.0;
	sprecdata->coarse_scaling  = OFF;
    
	sprecdata->M      = NULL;
	sprecdata->diag_M = NULL;
	sprecdata->B      = NULL;
	sprecdata->Bt     = NULL;
    sprecdata->S      = NULL;
    sprecdata->diag_S = NULL;
    sprecdata->P      = NULL;
    sprecdata->w      = NULL;
}

/**
 * \fn void fasp_precond_ns_param_set (precond_ns_param *sparam)
 * \brief Initialize itsolver_param
 * \param Input: pointer to itsolver_param
 * \date 02/29/2012 
 */
void fasp_precond_ns_param_set(precond_ns_param *sparam,
                                   input_param *Input)
{
    sparam->AMG_type = Input->AMG_type;
    sparam->print_level = Input->print_level;
    sparam->max_levels = Input->AMG_levels;
}

/**
 * \fn void fasp_precond_ns_data_init (itsolver_param *pdata)
 * \brief Initialize itsolver_param
 * \param Input: pointer to itsolver_param
 * \author Lu Wang
 * \date 02/29/2012
 */
void fasp_precond_ns_data_set (precond_ns_data *sprecdata,
                                    input_param *Input)
{
    sprecdata->print_level = Input->AMG_type;
	sprecdata->max_levels  = Input->AMG_levels;
	sprecdata->maxit    = Input->AMG_maxit;
	sprecdata->amg_tol     = Input->AMG_tol;
	sprecdata->cycle_type  = Input->AMG_cycle_type;
	sprecdata->smoother    = Input->AMG_smoother;
	sprecdata->presmooth_iter  = Input->AMG_presmooth_iter;
	sprecdata->postsmooth_iter = Input->AMG_postsmooth_iter;
	sprecdata->coarsening_type = Input->AMG_coarsening_type;;
	sprecdata->relaxation      = Input->AMG_relaxation;
	sprecdata->coarse_scaling  = Input->AMG_coarse_scaling;
}

/**
 * \fn void fasp_ns_param_init(input_param *Input,precond_ns_param *sparam,precond_ns_data  *sprecdata)
 * \brief Initialize parameters from *inparam
 *
 * \param *inparam     input parameters
 * \param *itparam     iterative solver parameters
 * \param *amgparam    AMG parameters
 * \param *iluparam    ILU parameters
 *
 * \author Lu Wang
 * \date 02/29/2012 
 *
 */
void fasp_ns_param_init(input_param *Input,
                            precond_ns_param *sparam,
                            precond_ns_data  *sprecdata)
{
    fasp_precond_ns_param_init (sparam);
    fasp_precond_ns_data_init (sprecdata);
    
    fasp_precond_ns_param_set(sparam,Input);
    fasp_precond_ns_data_set (sprecdata,Input);
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
