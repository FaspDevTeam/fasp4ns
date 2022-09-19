/*! \file  AuxParam.c
 *
 *  \brief Initialize, set, or print input data and parameters
 *
 *  \note  This file contains Level-0 (Aux) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2018 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <stdio.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_ns_param_init (input_ns_param *inparam,
 *                              ITS_param *itsparam,
 *                              AMG_ns_param *amgparam,
 *                              ILU_param *iluparam,
 *                              SWZ_param *swzparam)
 *
 * \brief Initialize parameters, global variables, etc
 *
 * \param inparam       Input parameters
 * \param itsparam      Iterative solver parameters
 * \param amgparam      AMG parameters
 * \param iluparam      ILU parameters
 * \param swzparam      Schwarz parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * Modified by Xiaozhe Hu on 05/27/2014
 */
void fasp_ns_param_init (input_ns_param *inparam,
                         itsolver_ns_param *itsparam,
                         AMG_ns_param *amgparam,
                         ILU_param *iluparam,
                         SWZ_param *swzparam)
{
#if CHMEM_MODE
    total_alloc_mem   = 0; // initialize total memeory amount
    total_alloc_count = 0; // initialize alloc count
#endif
    
    if (itsparam) fasp_ns_param_solver_init(itsparam);
    
    if (amgparam) fasp_ns_param_amg_init(amgparam);
    
    if (iluparam) fasp_param_ilu_init(iluparam);
    
    if (swzparam) fasp_param_swz_init(swzparam);
    
    if (inparam) {
        if (itsparam) fasp_ns_param_solver_set(itsparam,inparam);
        if (amgparam) fasp_ns_param_amg_set(amgparam,inparam);
        if (iluparam) fasp_ns_param_ilu_set(iluparam,inparam);
        if (swzparam) fasp_ns_param_swz_set(swzparam,inparam);
    }
    else {
        printf("### WARNING: No input specified. Use default values instead!\n");
    }
}

/**
 * \fn void fasp_ns_param_input_init (input_ns_param *inparam)
 *
 * \brief Initialize input parameters
 *
 * \param inparam    Input parameters
 *
 * \author Lu Wang
 * \date   2014/01/01
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 */
void fasp_ns_param_input_init (input_ns_param *inparam)
{
    strcpy(inparam->workdir,"../data/");
    
    // Input/output
    inparam->print_level                = PRINT_MIN;
    inparam->output_type                = 0;
    
    // Problem information
    inparam->problem_num                = 10;
    inparam->stop_type                  = STOP_REL_RES;
    
    // Solver parameters
    inparam->solver_type                = SOLVER_CG;
    inparam->precond_type               = PREC_AMG;
    inparam->itsolver_tol               = 1e-6;
    inparam->itsolver_maxit             = 500;
    inparam->restart                    = 25;
    inparam->itsolver_abstol            = 1e-8;
    
    //IR Solver parameters
    inparam->IR_type                    = 0;      
    inparam->IRsolver_tol               = 1e-8; 
    inparam->IRsolver_maxit             = 100;    
    
    // Velocity block parameters
    inparam->itsolver_type_v            = SOLVER_CG;
    inparam->precond_type_v             = PREC_AMG;
    inparam->pre_tol_v                  = 1e-2;
    inparam->pre_maxit_v                = 20;
    inparam->pre_restart_v              = 20;
    
    // Pressure block parameters
    inparam->itsolver_type_p            = SOLVER_CG;
    inparam->precond_type_p             = PREC_AMG;
    inparam->pre_tol_p                  = 1e-2;
    inparam->pre_maxit_p                = 20;
    inparam->pre_restart_p              = 20;
    
    // ILU method parameters
    inparam->ILU_type                   = ILUk;
    inparam->ILU_lfil                   = 0;
    inparam->ILU_droptol                = 0.001;
    inparam->ILU_relax                  = 0;
    inparam->ILU_permtol                = 0.0;
    
    // Schwarz method parameters
    inparam->SWZ_mmsize                 = 200;
    inparam->SWZ_maxlvl                 = 2;
    inparam->SWZ_type                   = 1;
    
    // AMG method parameters
    inparam->AMG_type_v                 = UA_AMG;
    inparam->AMG_levels_v               = 20;
    inparam->AMG_cycle_type_v           = V_CYCLE;
    inparam->AMG_smoother_v             = SMOOTHER_GS;
    inparam->AMG_smooth_order_v         = CF_ORDER;
    inparam->AMG_presmooth_iter_v       = 2;
    inparam->AMG_postsmooth_iter_v      = 2;
    inparam->AMG_relaxation_v           = 1.0;
    inparam->AMG_coarse_dof_v           = 500;
    inparam->AMG_tol_v                  = 1e-4*inparam->itsolver_tol;
    inparam->AMG_maxit_v                = 1;
    inparam->AMG_ILU_levels_v           = 0;
    inparam->AMG_schwarz_levels_v       = 0;
    inparam->AMG_coarse_scaling_v       = OFF; // Require investigation --Chensong
    inparam->AMG_amli_degree_v          = 1;
    inparam->AMG_nl_amli_krylov_type_v  = 2;
    
    // Classical AMG specific
    inparam->AMG_coarsening_type_v      = 1;
    inparam->AMG_interpolation_type_v   = 1;
    inparam->AMG_max_row_sum_v          = 0.9;
    inparam->AMG_strong_threshold_v     = 0.5;
    inparam->AMG_truncation_threshold_v = 0.4;
    inparam->AMG_aggressive_level_v     = 0;
    inparam->AMG_aggressive_path_v      = 1;
    
    // Aggregation AMG specific
    inparam->AMG_strong_coupled_v       = 0.08;
    inparam->AMG_max_aggregation_v      = 9;
    inparam->AMG_tentative_smooth_v     = 0.67;
    inparam->AMG_smooth_filter_v        = ON;
    
    // AMG method parameters
    inparam->AMG_type_p                 = UA_AMG;
    inparam->AMG_levels_p               = 20;
    inparam->AMG_cycle_type_p           = V_CYCLE;
    inparam->AMG_smoother_p             = SMOOTHER_GS;
    inparam->AMG_smooth_order_p         = CF_ORDER;
    inparam->AMG_presmooth_iter_p       = 2;
    inparam->AMG_postsmooth_iter_p      = 2;
    inparam->AMG_relaxation_p           = 1.0;
    inparam->AMG_coarse_dof_p           = 500;
    inparam->AMG_tol_p                  = 1e-4*inparam->itsolver_tol;
    inparam->AMG_maxit_p                = 1;
    inparam->AMG_ILU_levels_p           = 0;
    inparam->AMG_schwarz_levels_p       = 0;
    inparam->AMG_coarse_scaling_p       = OFF; // Require investigation --Chensong
    inparam->AMG_amli_degree_p          = 1;
    inparam->AMG_nl_amli_krylov_type_p  = 2;
    
    // Classical AMG specific
    inparam->AMG_coarsening_type_p      = 1;
    inparam->AMG_interpolation_type_p   = 1;
    inparam->AMG_max_row_sum_p          = 0.9;
    inparam->AMG_strong_threshold_p     = 0.5;
    inparam->AMG_truncation_threshold_p = 0.4;
    inparam->AMG_aggressive_level_p     = 0;
    inparam->AMG_aggressive_path_p      = 1;
    
    // Aggregation AMG specific
    inparam->AMG_strong_coupled_p       = 0.08;
    inparam->AMG_max_aggregation_p      = 9;
    inparam->AMG_tentative_smooth_p     = 0.67;
    inparam->AMG_smooth_filter_p        = ON;
}

/**
 * \fn void fasp_ns_param_amg_init (AMG_param *amgparam)
 *
 * \brief Initialize AMG parameters
 *
 * \param amgparam    Parameters for AMG
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 */
void fasp_ns_param_amg_init (AMG_ns_param *amgparam)
{
    fasp_param_amg_init(&(amgparam->param_v));
    fasp_param_amg_init(&(amgparam->param_p));
}

/**
 * \fn void fasp_ns_param_amg_set (AMG_ns_param *param, input_ns_param *inparam)
 *
 * \brief Set AMG_param from INPUT
 *
 * \param param     Parameters for AMG
 * \param inparam   Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 */
void fasp_ns_param_amg_set (AMG_ns_param   *param,
                            input_ns_param *inparam)
{
    // iterative solver parameter for the velocity block
    param->param_v.AMG_type    = inparam->AMG_type_v;
    param->param_v.print_level = inparam->print_level;
    
    if (inparam->itsolver_type_v == SOLVER_AMG) {
        param->param_v.maxit = inparam->pre_maxit_v;
        param->param_v.tol   = inparam->pre_tol_v;
    }
    else if (inparam->itsolver_type_v == SOLVER_FMG) {
        param->param_v.maxit = inparam->pre_maxit_v;
        param->param_v.tol   = inparam->pre_tol_v;
    }
    else {
        param->param_v.maxit = inparam->AMG_maxit_v;
        param->param_v.tol   = inparam->AMG_tol_v;
    }
    
    param->param_v.max_levels           = inparam->AMG_levels_v;
    param->param_v.cycle_type           = inparam->AMG_cycle_type_v;
    param->param_v.smoother             = inparam->AMG_smoother_v;
    param->param_v.smooth_order         = inparam->AMG_smooth_order_v;
    param->param_v.relaxation           = inparam->AMG_relaxation_v;
    param->param_v.polynomial_degree    = inparam->AMG_polynomial_degree_v;
    param->param_v.presmooth_iter       = inparam->AMG_presmooth_iter_v;
    param->param_v.postsmooth_iter      = inparam->AMG_postsmooth_iter_v;
    param->param_v.coarse_dof           = inparam->AMG_coarse_dof_v;
    param->param_v.coarse_solver        = inparam->AMG_coarse_solver_v;
    param->param_v.coarse_scaling       = inparam->AMG_coarse_scaling_v;
    param->param_v.amli_degree          = inparam->AMG_amli_degree_v;
    param->param_v.amli_coef            = NULL;
    param->param_v.nl_amli_krylov_type  = inparam->AMG_nl_amli_krylov_type_v;
    
    param->param_v.coarsening_type      = inparam->AMG_coarsening_type_v;
    param->param_v.interpolation_type   = inparam->AMG_interpolation_type_v;
    param->param_v.strong_threshold     = inparam->AMG_strong_threshold_v;
    param->param_v.truncation_threshold = inparam->AMG_truncation_threshold_v;
    param->param_v.max_row_sum          = inparam->AMG_max_row_sum_v;
    param->param_v.aggressive_level     = inparam->AMG_aggressive_level_v;
    param->param_v.aggressive_path      = inparam->AMG_aggressive_path_v;
    
    param->param_v.aggregation_type     = inparam->AMG_aggregation_type_v;
    param->param_v.pair_number          = inparam->AMG_pair_number_v;
    param->param_v.quality_bound        = inparam->AMG_quality_bound_v;
    
    param->param_v.strong_coupled       = inparam->AMG_strong_coupled_v;
    param->param_v.max_aggregation      = inparam->AMG_max_aggregation_v;
    param->param_v.tentative_smooth     = inparam->AMG_tentative_smooth_v;
    param->param_v.smooth_filter        = inparam->AMG_smooth_filter_v;
    
    param->param_v.ILU_levels           = inparam->AMG_ILU_levels_v;
    param->param_v.ILU_type             = inparam->ILU_type;
    param->param_v.ILU_lfil             = inparam->ILU_lfil;
    param->param_v.ILU_droptol          = inparam->ILU_droptol;
    param->param_v.ILU_relax            = inparam->ILU_relax;
    param->param_v.ILU_permtol          = inparam->ILU_permtol;
    param->param_v.SWZ_levels           = inparam->AMG_schwarz_levels_v;
    param->param_v.SWZ_mmsize           = inparam->SWZ_mmsize;
    param->param_v.SWZ_maxlvl           = inparam->SWZ_maxlvl;
    param->param_v.SWZ_type             = inparam->SWZ_type;
    
    // iterative solver parameter for the pressure block
    param->param_p.AMG_type    = inparam->AMG_type_p;
    param->param_p.print_level = inparam->print_level;
    
    if (inparam->itsolver_type_p == SOLVER_AMG) {
        param->param_p.maxit   = inparam->pre_maxit_p;
        param->param_p.tol     = inparam->pre_tol_p;
    }
    else if (inparam->itsolver_type_p == SOLVER_FMG) {
        param->param_p.maxit   = inparam->pre_maxit_p;
        param->param_p.tol     = inparam->pre_tol_p;
    }
    else {
        param->param_p.maxit   = inparam->AMG_maxit_p;
        param->param_p.tol     = inparam->AMG_tol_p;
    }
    
    param->param_p.max_levels           = inparam->AMG_levels_p;
    param->param_p.cycle_type           = inparam->AMG_cycle_type_p;
    param->param_p.smoother             = inparam->AMG_smoother_p;
    param->param_p.smooth_order         = inparam->AMG_smooth_order_p;
    param->param_p.relaxation           = inparam->AMG_relaxation_p;
    param->param_p.polynomial_degree    = inparam->AMG_polynomial_degree_p;
    param->param_p.presmooth_iter       = inparam->AMG_presmooth_iter_p;
    param->param_p.postsmooth_iter      = inparam->AMG_postsmooth_iter_p;
    param->param_p.coarse_dof           = inparam->AMG_coarse_dof_p;
    param->param_p.coarse_solver        = inparam->AMG_coarse_solver_p;
    param->param_p.coarse_scaling       = inparam->AMG_coarse_scaling_p;
    param->param_p.amli_degree          = inparam->AMG_amli_degree_p;
    param->param_p.amli_coef            = NULL;
    param->param_p.nl_amli_krylov_type  = inparam->AMG_nl_amli_krylov_type_p;
    
    param->param_p.coarsening_type      = inparam->AMG_coarsening_type_p;
    param->param_p.interpolation_type   = inparam->AMG_interpolation_type_p;
    param->param_p.strong_threshold     = inparam->AMG_strong_threshold_p;
    param->param_p.truncation_threshold = inparam->AMG_truncation_threshold_p;
    param->param_p.max_row_sum          = inparam->AMG_max_row_sum_p;
    param->param_p.aggressive_level     = inparam->AMG_aggressive_level_p;
    param->param_p.aggressive_path      = inparam->AMG_aggressive_path_p;
    
    param->param_p.aggregation_type     = inparam->AMG_aggregation_type_p;
    param->param_p.pair_number          = inparam->AMG_pair_number_p;
    param->param_p.quality_bound        = inparam->AMG_quality_bound_p;
    
    param->param_p.strong_coupled       = inparam->AMG_strong_coupled_p;
    param->param_p.max_aggregation      = inparam->AMG_max_aggregation_p;
    param->param_p.tentative_smooth     = inparam->AMG_tentative_smooth_p;
    param->param_p.smooth_filter        = inparam->AMG_smooth_filter_p;
    
    param->param_p.ILU_levels           = inparam->AMG_ILU_levels_p;
    param->param_p.ILU_type             = inparam->ILU_type;
    param->param_p.ILU_lfil             = inparam->ILU_lfil;
    param->param_p.ILU_droptol          = inparam->ILU_droptol;
    param->param_p.ILU_relax            = inparam->ILU_relax;
    param->param_p.ILU_permtol          = inparam->ILU_permtol;
    param->param_p.SWZ_levels           = inparam->AMG_schwarz_levels_p;
    param->param_p.SWZ_mmsize           = inparam->SWZ_mmsize;
    param->param_p.SWZ_maxlvl           = inparam->SWZ_maxlvl;
    param->param_p.SWZ_type             = inparam->SWZ_type;
}

/**
 * \fn void fasp_ns_param_solver_init (itsolver_ns_param *itsparam)
 *
 * \brief Initialize parameters for iterative solvers
 *
 * \param itsparam    Parameters for AMG
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 */
void fasp_ns_param_solver_init (itsolver_ns_param *itsparam)
{
    itsparam->itsolver_type   = SOLVER_CG;
    itsparam->precond_type    = PREC_AMG;
    itsparam->stop_type       = STOP_REL_RES;
    itsparam->maxit           = 100;
    itsparam->tol             = 1e-8;
    itsparam->abstol          = 1e-8;
    itsparam->restart         = 20;
    itsparam->print_level     = 0;
    
    itsparam->IR_type         = 0;
    itsparam->IRtol           = 1e-8;
    itsparam->IRmaxit         = 100;
    // iterative solver parameter for the velocity block
    itsparam->itsolver_type_v = SOLVER_CG;
    itsparam->precond_type_v  = PREC_AMG;
    itsparam->pre_maxit_v     = 20;
    itsparam->pre_tol_v       = 1e-2;
    itsparam->pre_restart_v   = 20;
    itsparam->print_level_v   = 0;
    
    // iterative solver parameter for the pressure block
    itsparam->itsolver_type_p = SOLVER_CG;
    itsparam->precond_type_p  = PREC_AMG;
    itsparam->pre_maxit_p     = 20;
    itsparam->pre_tol_p       = 1e-2;
    itsparam->pre_restart_p   = 20;
    itsparam->print_level_p   = 0;
}

/**
 * \fn void fasp_ns_param_solver_set (itsolver_ns_param *itsparam,
 *                                    input_ns_param *inparam)
 *
 * \brief Set ITS_param with INPUT
 *
 * \param itsparam   Parameters for iterative solvers
 * \param inparam    Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * Modified by Xiaozhe Hu on 02/21/2014
 */
void fasp_ns_param_solver_set (itsolver_ns_param *itsparam,
                               input_ns_param *inparam)
{
    itsparam->print_level    = inparam->print_level;
    itsparam->itsolver_type  = inparam->solver_type;
    itsparam->precond_type   = inparam->precond_type;
    itsparam->stop_type      = inparam->stop_type;
    itsparam->restart        = inparam->restart;
    
    itsparam->tol   = inparam->itsolver_tol;
    itsparam->maxit = inparam->itsolver_maxit;  
    itsparam->abstol =inparam->itsolver_abstol;
    itsparam->IR_type   = inparam->IR_type;
    itsparam->IRtol   = inparam->IRsolver_tol;
    itsparam->IRmaxit = inparam->IRsolver_maxit;


    // iterative solver parameter for the velocity block
    itsparam->itsolver_type_v  = inparam->itsolver_type_v;
    itsparam->precond_type_v   = inparam->precond_type_v;
    itsparam->pre_restart_v    = inparam->pre_restart_v;
    
    if (itsparam->itsolver_type_v == SOLVER_AMG) {
        itsparam->pre_tol_v    = inparam->AMG_tol_v;
        itsparam->pre_abstol_v    = inparam->pre_abstol_v;
        itsparam->pre_maxit_v  = inparam->AMG_maxit_v;
    }
    else {
        itsparam->pre_tol_v    = inparam->pre_tol_v;
        itsparam->pre_abstol_v    = inparam->pre_abstol_v;
        itsparam->pre_maxit_v  = inparam->pre_maxit_v;
    }
    
    // iterative solver parameter for the pressure block
    itsparam->itsolver_type_p  = inparam->itsolver_type_p;
    itsparam->precond_type_p   = inparam->precond_type_p;
    itsparam->pre_restart_p    = inparam->pre_restart_p;
    
    if (itsparam->itsolver_type_p == SOLVER_AMG) {
        itsparam->pre_tol_p    = inparam->AMG_tol_p;
        itsparam->pre_abstol_p    = inparam->pre_abstol_p;
        itsparam->pre_maxit_p  = inparam->AMG_maxit_p;
    }
    else {
        itsparam->pre_tol_p    = inparam->pre_tol_p;
        itsparam->pre_abstol_p    = inparam->pre_abstol_p;
        itsparam->pre_maxit_p  = inparam->pre_maxit_p;
    }
    
    if (itsparam->print_level > 2) {
        itsparam->print_level_v = itsparam->print_level-3;
        itsparam->print_level_p = itsparam->print_level-3;
    }
    else {
        itsparam->print_level_v = 0;
        itsparam->print_level_p = 0;
    }
    
}

/**
 * \fn void fasp_ns_param_ilu_set (ILU_param *iluparam, input_ns_param *inparam)
 *
 * \brief Set ILU_param with INPUT
 *
 * \param iluparam    Parameters for ILU
 * \param inparam     Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 */
void fasp_ns_param_ilu_set (ILU_param *iluparam,
                            input_ns_param *inparam)
{
    iluparam->print_level = inparam->print_level;
    iluparam->ILU_type    = inparam->ILU_type;
    iluparam->ILU_lfil    = inparam->ILU_lfil;
    iluparam->ILU_droptol = inparam->ILU_droptol;
    iluparam->ILU_relax   = inparam->ILU_relax;
    iluparam->ILU_permtol = inparam->ILU_permtol;
}

/**
 * \fn void fasp_ns_param_swz_set (SWZ_param *swzparam, input_ns_param *inparam)
 *
 * \brief Set SWZ_param with INPUT
 *
 * \param swzparam    Parameters for Schwarz method
 * \param inparam     Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 */
void fasp_ns_param_swz_set (SWZ_param *swzparam,
                            input_ns_param *inparam)
{
    swzparam->print_level = inparam->print_level;
    swzparam->SWZ_type    = inparam->SWZ_type;
    swzparam->SWZ_maxlvl  = inparam->SWZ_maxlvl;
    swzparam->SWZ_mmsize  = inparam->SWZ_mmsize;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
