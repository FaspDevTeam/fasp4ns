/*! \file parameters.c
 *
 *  \brief Initialize, set, or print input data and parameters
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
 *                           itsolver_param *itsparam,
 *                           AMG_param *amgparam,
 *                           ILU_param *iluparam,
 *                           Schwarz_param *schparam)
 *
 * \brief Initialize parameters, global variables, etc
 *
 * \param inparam       Input parameters
 * \param itsparam      Iterative solver parameters
 * \param amgparam      AMG parameters
 * \param iluparam      ILU parameters
 * \param schparam      Schwarz parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 *
 */
void fasp_ns_param_init (input_ns_param *inparam,
                      itsolver_ns_param *itsparam,
                      AMG_ns_param *amgparam,
                      ILU_param *iluparam,
                      Schwarz_param *schparam)
{
#if CHMEM_MODE
    total_alloc_mem   = 0; // initialize total memeory amount
    total_alloc_count = 0; // initialize alloc count
#endif
    
    if (itsparam) fasp_ns_param_solver_init(itsparam);
    
    if (amgparam) fasp_ns_param_amg_init(amgparam);
    
    if (iluparam) fasp_param_ilu_init(iluparam);
    
    if (schparam) fasp_param_schwarz_init(schparam);
    
    if (inparam) {
        if (itsparam) fasp_ns_param_solver_set(itsparam,inparam);
        if (amgparam) fasp_ns_param_amg_set(amgparam,inparam);
        if (iluparam) fasp_ns_param_ilu_set(iluparam,inparam);
        if (schparam) fasp_ns_param_schwarz_set(schparam,inparam);
    }
    else {
        printf("### WARNING: No input specified. Use default values instead!\n");
    }
}

/**
 * \fn void fasp_param_input_init (input_ns_param *inparam)
 *
 * \brief Initialize input parameters
 *
 * \param inparam    Input parameters
 *
 * \author Lu Wang
 * \date   2014/01/1
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 *
 */
void fasp_ns_param_input_init (input_ns_param *inparam)
{
    strcpy(inparam->workdir,"../data/");
    
    // Input/output
    inparam->print_level              = PRINT_MIN;
    inparam->output_type              = 0;
    
    // Problem information
    inparam->problem_num              = 10;
    inparam->stop_type                = STOP_REL_RES;
    
    // Solver parameters
    inparam->solver_type              = SOLVER_CG;
    inparam->precond_type             = PREC_AMG;
    inparam->itsolver_tol             = 1e-6;
    inparam->itsolver_maxit           = 500;
    inparam->restart                  = 25;
    
    // Velocity block parameters
    inparam->itsolver_type_v             = SOLVER_CG;
    inparam->precond_type_v            = PREC_AMG;
    inparam->pre_tol_v                 = 1e-2;
    inparam->pre_maxit_v               = 20;
    inparam->pre_restart_v             = 20;
    
    // Pressure block parameters
    inparam->itsolver_type_p             = SOLVER_CG;
    inparam->precond_type_p            = PREC_AMG;
    inparam->pre_tol_p                 = 1e-2;
    inparam->pre_maxit_p               = 20;
    inparam->pre_restart_p             = 20;
    
    // ILU method parameters
    inparam->ILU_type                 = ILUk;
    inparam->ILU_lfil                 = 0;
    inparam->ILU_droptol              = 0.001;
    inparam->ILU_relax                = 0;
    inparam->ILU_permtol              = 0.0;
    
    // Schwarz method parameters
    inparam->Schwarz_mmsize           = 200;
	inparam->Schwarz_maxlvl           = 2;
	inparam->Schwarz_type             = 1;
    
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
    inparam->AMG_type_p                = UA_AMG;
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
 * \fn void fasp_param_amg_init (AMG_param *amgparam)
 *
 * \brief Initialize AMG parameters
 *
 * \param amgparam    Parameters for AMG
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 *
 */
void fasp_ns_param_amg_init (AMG_ns_param *amgparam)
{
    fasp_param_amg_init(&(amgparam->param_v));
    fasp_param_amg_init(&(amgparam->param_p));
    /*
    // AMG type 
    amgparam->AMG_type             = UA_AMG;
    amgparam->print_level          = PRINT_NONE;
    amgparam->maxit                = 1;
    amgparam->tol                  = 1e-8;
    
    // AMG method parameters
    amgparam->max_levels           = 15;
    amgparam->coarse_dof           = 500;
    amgparam->cycle_type           = V_CYCLE;
    amgparam->smoother             = SMOOTHER_GS;
    amgparam->smooth_order         = CF_ORDER;
    amgparam->presmooth_iter       = 2;
    amgparam->postsmooth_iter      = 2;
    amgparam->relaxation           = 1.0;
    amgparam->polynomial_degree    = 3;
    amgparam->coarse_scaling       = OFF; // Require investigation --Chensong
    amgparam->amli_degree          = 1;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = 2;
    
    // Classical AMG specific
    amgparam->coarsening_type      = 1;
    amgparam->interpolation_type   = 1;
    amgparam->max_row_sum          = 0.9;
    amgparam->strong_threshold     = 0.5;
    amgparam->truncation_threshold = 0.4;
    amgparam->aggressive_level     = 0;
    amgparam->aggressive_path      = 1;
    
    // Aggregation AMG specific
    amgparam->strong_coupled       = 0.08;
    amgparam->max_aggregation      = 9;
    amgparam->tentative_smooth     = 0.0;
    amgparam->smooth_filter        = OFF;
    
    // AMG type
    amgparam->AMG_p_type             = UA_AMG;
    amgparam->p_print_level          = PRINT_NONE;
    amgparam->p_maxit                = 1;
    amgparam->p_tol                  = 1e-2;
    
    // AMG method parameters
    amgparam->p_max_levels           = 15;
    amgparam->p_coarse_dof           = 500;
    amgparam->p_cycle_type           = V_CYCLE;
    amgparam->p_smoother             = SMOOTHER_GS;
    amgparam->p_smooth_order         = CF_ORDER;
    amgparam->p_presmooth_iter       = 2;
    amgparam->p_postsmooth_iter      = 2;
    amgparam->p_relaxation           = 1.0;
    amgparam->p_polynomial_degree    = 3;
    amgparam->p_coarse_scaling       = OFF; // Require investigation --Chensong
    amgparam->p_amli_degree          = 1;
    amgparam->p_amli_coef            = NULL;
    amgparam->p_nl_amli_krylov_type  = 2;
    
    // Classical AMG specific
    amgparam->p_coarsening_type      = 1;
    amgparam->p_interpolation_type   = 1;
    amgparam->p_max_row_sum          = 0.9;
    amgparam->p_strong_threshold     = 0.5;
    amgparam->p_truncation_threshold = 0.4;
    amgparam->p_aggressive_level     = 0;
    amgparam->p_aggressive_path      = 1;
    
    // Aggregation AMG specific
    amgparam->p_strong_coupled       = 0.08;
    amgparam->p_max_aggregation      = 9;
    amgparam->p_tentative_smooth     = 0.0;
    amgparam->p_smooth_filter        = OFF;
    
    // ILU smoother parameters
    amgparam->p_ILU_levels           = 0;
    amgparam->ILU_type             = ILUk;
    amgparam->ILU_levels           = 0;
    amgparam->ILU_lfil             = 0;
    amgparam->ILU_droptol          = 0.001;
    amgparam->ILU_relax            = 0;
    
    // Schwarz smoother parameters
    amgparam->p_schwarz_levels       = 0;
    amgparam->schwarz_levels       = 0;
    amgparam->schwarz_mmsize       = 200;
    amgparam->schwarz_maxlvl       = 2;
    amgparam->schwarz_type         = 1;
     */
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
 * \note Xiaozhe Hu modified on 02/21/2014
 *
 */
void fasp_ns_param_amg_set (AMG_ns_param *param,
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
    param->param_v.schwarz_levels       = inparam->AMG_schwarz_levels_v;
	param->param_v.schwarz_mmsize       = inparam->Schwarz_mmsize;
	param->param_v.schwarz_maxlvl       = inparam->Schwarz_maxlvl;
	param->param_v.schwarz_type         = inparam->Schwarz_type;

    
    // iterative solver parameter for the pressure block
    param->param_p.AMG_type    = inparam->AMG_type_p;
    param->param_p.print_level = inparam->print_level;
    
    if (inparam->itsolver_type_p == SOLVER_AMG) {
        param->param_p.maxit = inparam->pre_maxit_p;
        param->param_p.tol   = inparam->pre_tol_p;
    }
    else if (inparam->itsolver_type_p == SOLVER_FMG) {
        param->param_p.maxit = inparam->pre_maxit_p;
        param->param_p.tol   = inparam->pre_tol_p;
    }
    else {
        param->param_p.maxit = inparam->AMG_maxit_p;
        param->param_p.tol   = inparam->AMG_tol_p;
    }
    
    param->param_p.max_levels          = inparam->AMG_levels_p;
    param->param_p.cycle_type          = inparam->AMG_cycle_type_p;
    param->param_p.smoother            = inparam->AMG_smoother_p;
    param->param_p.smooth_order        = inparam->AMG_smooth_order_p;
    param->param_p.relaxation          = inparam->AMG_relaxation_p;
    param->param_p.polynomial_degree   = inparam->AMG_polynomial_degree_p;
    param->param_p.presmooth_iter      = inparam->AMG_presmooth_iter_p;
    param->param_p.postsmooth_iter     = inparam->AMG_postsmooth_iter_p;
    param->param_p.coarse_dof          = inparam->AMG_coarse_dof_p;
    param->param_p.coarse_scaling      = inparam->AMG_coarse_scaling_p;
    param->param_p.amli_degree         = inparam->AMG_amli_degree_p;
    param->param_p.amli_coef           = NULL;
    param->param_p.nl_amli_krylov_type = inparam->AMG_nl_amli_krylov_type_p;
    
    param->param_p.coarsening_type     = inparam->AMG_coarsening_type_p;
    param->param_p.interpolation_type  = inparam->AMG_interpolation_type_p;
    param->param_p.strong_threshold    = inparam->AMG_strong_threshold_p;
    param->param_p.truncation_threshold= inparam->AMG_truncation_threshold_p;
    param->param_p.max_row_sum         = inparam->AMG_max_row_sum_p;
    param->param_p.aggressive_level    = inparam->AMG_aggressive_level_p;
    param->param_p.aggressive_path     = inparam->AMG_aggressive_path_p;
    
    param->param_p.strong_coupled      = inparam->AMG_strong_coupled_p;
    param->param_p.max_aggregation     = inparam->AMG_max_aggregation_p;
    param->param_p.tentative_smooth    = inparam->AMG_tentative_smooth_p;
    param->param_p.smooth_filter       = inparam->AMG_smooth_filter_p;
    
    
    param->param_p.ILU_levels           = inparam->AMG_ILU_levels_p;
    param->param_p.ILU_type             = inparam->ILU_type;
    param->param_p.ILU_lfil             = inparam->ILU_lfil;
    param->param_p.ILU_droptol          = inparam->ILU_droptol;
    param->param_p.ILU_relax            = inparam->ILU_relax;
    param->param_p.ILU_permtol          = inparam->ILU_permtol;
    param->param_p.schwarz_levels       = inparam->AMG_schwarz_levels_p;
	param->param_p.schwarz_mmsize       = inparam->Schwarz_mmsize;
	param->param_p.schwarz_maxlvl       = inparam->Schwarz_maxlvl;
	param->param_p.schwarz_type         = inparam->Schwarz_type;
    /*
     param->AMG_type    = inparam->AMG_type;
     param->print_level = inparam->print_level;
     
     if (inparam->solver_v_type == SOLVER_AMG) {
     param->maxit = inparam->itsolver_maxit;
     param->tol   = inparam->itsolver_tol;
     }
     else if (inparam->solver_v_type == SOLVER_FMG) {
     param->maxit = inparam->itsolver_maxit;
     param->tol   = inparam->itsolver_tol;
     }
     else {
     param->maxit = inparam->AMG_maxit;
     param->tol   = inparam->AMG_tol;
     }
     
     param->max_levels           = inparam->AMG_levels;
     param->cycle_type           = inparam->AMG_cycle_type;
     param->smoother             = inparam->AMG_smoother;
     param->smooth_order         = inparam->AMG_smooth_order;
     param->relaxation           = inparam->AMG_relaxation;
     param->polynomial_degree    = inparam->AMG_polynomial_degree;
     param->presmooth_iter       = inparam->AMG_presmooth_iter;
     param->postsmooth_iter      = inparam->AMG_postsmooth_iter;
     param->coarse_dof           = inparam->AMG_coarse_dof;
     param->coarse_scaling       = inparam->AMG_coarse_scaling;
     param->amli_degree          = inparam->AMG_amli_degree;
     param->amli_coef            = NULL;
     param->nl_amli_krylov_type  = inparam->AMG_nl_amli_krylov_type;
     
     param->coarsening_type      = inparam->AMG_coarsening_type;
     param->interpolation_type   = inparam->AMG_interpolation_type;
     param->strong_threshold     = inparam->AMG_strong_threshold;
     param->truncation_threshold = inparam->AMG_truncation_threshold;
     param->max_row_sum          = inparam->AMG_max_row_sum;
     param->aggressive_level     = inparam->AMG_aggressive_level;
     param->aggressive_path      = inparam->AMG_aggressive_path;
     
     param->strong_coupled       = inparam->AMG_strong_coupled;
     param->max_aggregation      = inparam->AMG_max_aggregation;
     param->tentative_smooth     = inparam->AMG_tentative_smooth;
     param->smooth_filter        = inparam->AMG_smooth_filter;
     
     // pressure block AMG parameters
     param->AMG_p_type    = inparam->AMG_p_type;
     param->p_print_level = inparam->print_level;
     
     if (inparam->solver_p_type == SOLVER_AMG) {
     param->p_maxit = inparam->pre_p_maxit;
     param->p_tol   = inparam->pre_p_tol;
     }
     else if (inparam->solver_p_type == SOLVER_FMG) {
     param->maxit = inparam->pre_p_maxit;
     param->tol   = inparam->pre_p_tol;
     }
     else {
     param->maxit = inparam->AMG_p_maxit;
     param->tol   = inparam->AMG_p_tol;
     }
     
     param->p_max_levels          = inparam->AMG_p_levels;
     param->p_cycle_type          = inparam->AMG_p_cycle_type;
     param->p_smoother            = inparam->AMG_p_smoother;
     param->p_smooth_order        = inparam->AMG_p_smooth_order;
     param->p_relaxation          = inparam->AMG_p_relaxation;
     param->p_polynomial_degree   = inparam->AMG_p_polynomial_degree;
     param->p_presmooth_iter      = inparam->AMG_p_presmooth_iter;
     param->p_postsmooth_iter     = inparam->AMG_p_postsmooth_iter;
     param->p_coarse_dof          = inparam->AMG_p_coarse_dof;
     param->p_coarse_scaling      = inparam->AMG_p_coarse_scaling;
     param->p_amli_degree         = inparam->AMG_p_amli_degree;
     param->p_amli_coef           = NULL;
     param->p_nl_amli_krylov_type = inparam->AMG_p_nl_amli_krylov_type;
     
     param->p_coarsening_type     = inparam->AMG_p_coarsening_type;
     param->p_interpolation_type  = inparam->AMG_p_interpolation_type;
     param->p_strong_threshold    = inparam->AMG_p_strong_threshold;
     param->p_truncation_threshold= inparam->AMG_p_truncation_threshold;
     param->p_max_row_sum         = inparam->AMG_p_max_row_sum;
     param->p_aggressive_level    = inparam->AMG_p_aggressive_level;
     param->p_aggressive_path     = inparam->AMG_p_aggressive_path;
     
     param->p_strong_coupled      = inparam->AMG_p_strong_coupled;
     param->p_max_aggregation     = inparam->AMG_p_max_aggregation;
     param->p_tentative_smooth    = inparam->AMG_p_tentative_smooth;
     param->p_smooth_filter       = inparam->AMG_p_smooth_filter;
     
     param->ILU_levels           = inparam->AMG_ILU_levels;
     param->p_ILU_levels         = inparam->AMG_p_ILU_levels;
     param->ILU_type             = inparam->ILU_type;
     param->ILU_lfil             = inparam->ILU_lfil;
     param->ILU_droptol          = inparam->ILU_droptol;
     param->ILU_relax            = inparam->ILU_relax;
     param->ILU_permtol          = inparam->ILU_permtol;
     
     param->schwarz_levels       = inparam->AMG_schwarz_levels;
     param->p_schwarz_levels     = inparam->AMG_p_schwarz_levels;
     param->schwarz_mmsize       = inparam->Schwarz_mmsize;
     param->schwarz_maxlvl       = inparam->Schwarz_maxlvl;
     param->schwarz_type         = inparam->Schwarz_type;
     */
}

/**
 * \fn void fasp_ns_param_solver_init (itsolver_ns_param *itsparam)
 *
 * \brief Initialize AMG parameters
 *
 * \param amgparam    Parameters for AMG
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 */

void fasp_ns_param_solver_init(itsolver_ns_param *itsparam)
{
    itsparam->itsolver_type = SOLVER_CG;
	itsparam->precond_type  = PREC_AMG;
	itsparam->stop_type     = STOP_REL_RES;
	itsparam->maxit         = 100;
	itsparam->tol           = 1e-8;
	itsparam->restart       = 20;
	itsparam->print_level   = 0;
    
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
 * \fn void fasp_param_solver_set (itsolver_ns_param *itsparam, input_ns_param *inparam)
 *
 * \brief Set itsolver_param with INPUT
 *
 * \param itsparam   Parameters for iterative solvers
 * \param inparam    Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 *
 * \note Xiaozhe Hu modified on 02/21/2014
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
    
    // iterative solver parameter for the velocity block
    itsparam->itsolver_type_v  = inparam->itsolver_type_v;
    itsparam->precond_type_v   = inparam->precond_type_v;
    itsparam->pre_restart_v        = inparam->pre_restart_v;
    
    if (itsparam->itsolver_type_v == SOLVER_AMG) {
        itsparam->pre_tol_v   = inparam->AMG_tol_v;
        itsparam->pre_maxit_v = inparam->AMG_maxit_v;
    }
    else {
        itsparam->pre_tol_v   = inparam->pre_tol_v;
        itsparam->pre_maxit_v = inparam->pre_maxit_v;
    }
    
    // iterative solver parameter for the pressure block
    itsparam->itsolver_type_p  = inparam->itsolver_type_p;
    itsparam->precond_type_p   = inparam->precond_type_p;
    itsparam->pre_restart_p        = inparam->pre_restart_p;
    
    if (itsparam->itsolver_type_p == SOLVER_AMG) {
        itsparam->pre_tol_p   = inparam->AMG_tol_p;
        itsparam->pre_maxit_p = inparam->AMG_maxit_p;
    }
    else {
        itsparam->pre_tol_p   = inparam->pre_tol_p;
        itsparam->pre_maxit_p = inparam->pre_maxit_p;
    }
    if (itsparam->print_level > 2)
    { itsparam->print_level_v = itsparam->print_level-3;}
    else {itsparam->print_level_v = 0;}
    if (itsparam->print_level > 3)
    { itsparam->print_level_p = itsparam->print_level-4;}
    else {itsparam->print_level_p = 0;}
    
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
 * \fn void fasp_ns_param_schwarz_set (Schwarz_param *schparam, input_ns_param *inparam)
 *
 * \brief Set Schwarz_param with INPUT
 *
 * \param schparam    Parameters for Schwarz method
 * \param inparam     Input parameters
 *
 * \author Lu Wang
 * \date   2014/02/11
 */
void fasp_ns_param_schwarz_set (Schwarz_param *schparam,
                             input_ns_param *inparam)
{
    schparam->print_level    = inparam->print_level;
    schparam->schwarz_type   = inparam->Schwarz_type;
    schparam->schwarz_maxlvl = inparam->Schwarz_maxlvl;
    schparam->schwarz_mmsize = inparam->Schwarz_mmsize;
}


/**
 * \fn fasp_ns_amg_to_amg_param (AMG_param *amgparam, AMG_ns_param *inparam)
 *
 * \brief Set AMG_param with AMG_ns_param
 *
 * \param amgparam    Parameters for AMG
 * \param pcdata      Parameters for AMG for NS
 *
 * \author Lu Wang
 * \date   2014/02/11
 */
/*
void fasp_ns_amg_to_amg_param (AMG_param *param,
                             AMG_ns_param *inparam)
{
    param->AMG_type    = inparam->AMG_type;
    param->print_level = inparam->print_level;
    
    param->maxit = inparam->maxit;
    param->tol   = inparam->tol;
    
    param->max_levels           = inparam->max_levels;
    param->cycle_type           = inparam->cycle_type;
    param->smoother             = inparam->smoother;
    param->smooth_order         = inparam->smooth_order;
    param->relaxation           = inparam->relaxation;
    param->polynomial_degree    = inparam->polynomial_degree;
    param->presmooth_iter       = inparam->presmooth_iter;
    param->postsmooth_iter      = inparam->postsmooth_iter;
    param->coarse_dof           = inparam->coarse_dof;
    param->coarse_scaling       = inparam->coarse_scaling;
    param->amli_degree          = inparam->amli_degree;
    param->amli_coef            = NULL;
    param->nl_amli_krylov_type  = inparam->nl_amli_krylov_type;
    
    param->coarsening_type      = inparam->coarsening_type;
    param->interpolation_type   = inparam->interpolation_type;
    param->strong_threshold     = inparam->strong_threshold;
    param->truncation_threshold = inparam->truncation_threshold;
    param->max_row_sum          = inparam->max_row_sum;
    param->aggressive_level     = inparam->aggressive_level;
    param->aggressive_path      = inparam->aggressive_path;
    
    param->strong_coupled       = inparam->strong_coupled;
    param->max_aggregation      = inparam->max_aggregation;
    param->tentative_smooth     = inparam->tentative_smooth;
    param->smooth_filter        = inparam->smooth_filter;
    
}
 */
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
