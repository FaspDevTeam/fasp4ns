/*! \file input.c
 *  \brief Read input parameters
 */

#include "fasp.h"
#include "fasp_functs.h"

#include "fasp4ns.h"
#include "fasp4ns_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_param_check (input_ns_param *inparam)
 *
 * \brief Simple check on input parameters
 *
 * \param inparam     Input parameters
 *
 * \author Chensong Zhang
 * \date   09/29/2013
 *
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 * \note: modified by Xiaozhe Hu on May. 27, 2014
 *
 */
SHORT fasp_ns_param_check (input_ns_param *inparam)
{
    SHORT status = FASP_SUCCESS;
    
    if ( inparam->problem_num<0
        || inparam->solver_type<0
        || inparam->solver_type>50
        || inparam->precond_type<0
        || inparam->itsolver_tol<=0
        || inparam->itsolver_maxit<=0
        || inparam->stop_type<=0
        || inparam->stop_type>3
        || inparam->restart<0
        || inparam->ILU_type<=0
        || inparam->ILU_type>3
        || inparam->ILU_lfil<0
        || inparam->ILU_droptol<=0
        || inparam->ILU_relax<0
        || inparam->ILU_permtol<0
        || inparam->SWZ_mmsize<0
        || inparam->SWZ_maxlvl<0
        || inparam->SWZ_type<0
        || inparam->AMG_type_v<=0
        || inparam->AMG_type_v>3
        || inparam->AMG_cycle_type_v<=0
        || inparam->AMG_cycle_type_v>4
        || inparam->AMG_levels_v<0
        || inparam->AMG_ILU_levels_v<0
        || inparam->AMG_coarse_dof_v<=0
        || inparam->AMG_tol_v<0
        || inparam->AMG_maxit_v<0
        || inparam->AMG_coarsening_type_v<=0
        || inparam->AMG_coarsening_type_v>4
        || inparam->AMG_interpolation_type_v<0
        || inparam->AMG_interpolation_type_v>5
        || inparam->AMG_smoother_v<0
        || inparam->AMG_smoother_v>20
        || inparam->AMG_strong_threshold_v<0.0
        || inparam->AMG_strong_threshold_v>0.9999
        || inparam->AMG_truncation_threshold_v<0.0
        || inparam->AMG_truncation_threshold_v>0.9999
        || inparam->AMG_max_row_sum_v<0.0
        || inparam->AMG_presmooth_iter_v<0
        || inparam->AMG_postsmooth_iter_v<0
        || inparam->AMG_amli_degree_v<0
        || inparam->AMG_aggressive_level_v<0
        || inparam->AMG_aggressive_path_v<0
        || inparam->AMG_strong_coupled_v<0
        || inparam->AMG_max_aggregation_v<=0
        || inparam->AMG_tentative_smooth_v<0
        || inparam->AMG_smooth_filter_v<0
        || inparam->AMG_type_p<=0
        || inparam->AMG_type_p>3
        || inparam->AMG_cycle_type_p<=0
        || inparam->AMG_cycle_type_p>4
        || inparam->AMG_levels_p<0
        || inparam->AMG_ILU_levels_p<0
        || inparam->AMG_coarse_dof_p<=0
        || inparam->AMG_tol_p<0
        || inparam->AMG_maxit_p<0
        || inparam->AMG_coarsening_type_p<=0
        || inparam->AMG_coarsening_type_p>4
        || inparam->AMG_interpolation_type_p<0
        || inparam->AMG_interpolation_type_p>5
        || inparam->AMG_smoother_p<0
        || inparam->AMG_smoother_p>20
        || inparam->AMG_strong_threshold_p<0.0
        || inparam->AMG_strong_threshold_p>0.9999
        || inparam->AMG_truncation_threshold_p<0.0
        || inparam->AMG_truncation_threshold_p>0.9999
        || inparam->AMG_max_row_sum_p<0.0
        || inparam->AMG_presmooth_iter_p<0
        || inparam->AMG_postsmooth_iter_p<0
        || inparam->AMG_amli_degree_p<0
        || inparam->AMG_aggressive_level_p<0
        || inparam->AMG_aggressive_path_p<0
        || inparam->AMG_strong_coupled_p<0
        || inparam->AMG_max_aggregation_p<=0
        || inparam->AMG_tentative_smooth_p<0
        || inparam->AMG_smooth_filter_p<0
        ) status = ERROR_INPUT_PAR;
    
    return status;
}

/**
 * \fn void fasp_ns_param_input (char *filenm, input_ns_param *Input)
 *
 * \brief Read input parameters for NS problem from disk file
 *
 * \param filenm    File name for input file
 * \param Input     Input parameters
 *
 * \author Lu Wang
 * \date   02/15/2012
 *
 * \note Xiaozhe Hu modified on 02/21/2014
 * \note: modified by Xiaozhe Hu on May. 27, 2014
 *
 */
void fasp_ns_param_input (char *filenm,
                          input_ns_param *Input)
{
    char     buffer[500]; // Note: max number of char for each line!
    char   * wall;
    INT      val;
    SHORT    status = FASP_SUCCESS;
    
    // set default input parameters
    fasp_ns_param_input_init(Input);
    
    // if input file is not specified, use the default values
    if (filenm==NULL) return;
    
    FILE *fp = fopen(filenm,"r");
    if (fp==NULL) {
        printf("### ERROR: Could not open file %s...\n", filenm);
        exit(ERROR_OPEN_FILE);
    }
    
    while ( status == FASP_SUCCESS ) {
        INT   ibuff;
        REAL  dbuff;
        char  sbuff[500];
        
        val = fscanf(fp,"%s",buffer);
        if (val==EOF) break;
        if (val!=1){ status = ERROR_INPUT_PAR; break; }
        if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
            wall = fgets(buffer,500,fp); // skip rest of line
            continue;
        }
        
        // match keyword and scan for value
        if (strcmp(buffer,"workdir")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(Input->workdir,sbuff,128);
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"problem_num")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->problem_num=ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"print_level")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->print_level = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"output_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->output_type = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"solver_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->solver_type = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"precond_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->precond_type = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"stop_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->stop_type = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_tol")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_tol = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_maxit")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_maxit = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"solver_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_type_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"precond_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->precond_type_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_tol_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_tol_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_maxit_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_maxit_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_restart_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_restart_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"solver_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_type_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"precond_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->precond_type_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_tol_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_tol_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_maxit_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_maxit_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_restart_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->pre_restart_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"itsolver_restart")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->restart = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_ILU_levels_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_ILU_levels_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_schwarz_levels_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = FASP_SUCCESS; break; }
            Input->AMG_schwarz_levels_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        
        else if (strcmp(buffer,"AMG_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
                Input->AMG_type_v = CLASSIC_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                Input->AMG_type_v = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                Input->AMG_type_v = UA_AMG;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggregation_type_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggregation_type_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_pair_number_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_pair_number_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_quality_bound_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_quality_bound_v = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_strong_coupled_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_strong_coupled_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_max_aggregation_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_aggregation_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_tentative_smooth_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tentative_smooth_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smooth_filter_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                Input->AMG_smooth_filter_v = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_smooth_filter_v = OFF;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_scaling_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                Input->AMG_coarse_scaling_v = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_coarse_scaling_v = OFF;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_levels_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_levels_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_tol_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tol_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_maxit_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_maxit_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_dof_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarse_dof_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_solver_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarse_solver_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_cycle_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
                Input->AMG_cycle_type_v = V_CYCLE;
            else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
                Input->AMG_cycle_type_v = W_CYCLE;
            else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
                Input->AMG_cycle_type_v = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
                Input->AMG_cycle_type_v = NL_AMLI_CYCLE;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smoother_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
                Input->AMG_smoother_v = SMOOTHER_JACOBI;
            else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
                Input->AMG_smoother_v = SMOOTHER_GS;
            else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
                Input->AMG_smoother_v = SMOOTHER_SGS;
            else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
                Input->AMG_smoother_v = SMOOTHER_CG;
            else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
                Input->AMG_smoother_v = SMOOTHER_SOR;
            else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
                Input->AMG_smoother_v = SMOOTHER_SSOR;
            else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
                Input->AMG_smoother_v = SMOOTHER_GSOR;
            else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
                Input->AMG_smoother_v = SMOOTHER_SGSOR;
            else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
                Input->AMG_smoother_v = SMOOTHER_POLY;
            else if ((strcmp(buffer,"L1_DIAG")==0)||(strcmp(buffer,"l1_diag")==0))
                Input->AMG_smoother_v = SMOOTHER_L1DIAG;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smooth_order_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
                Input->AMG_smooth_order_v = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
                Input->AMG_smooth_order_v = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarsening_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarsening_type_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_interpolation_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_interpolation_type_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_level_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_level_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_path_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_path_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_presmooth_iter_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_presmooth_iter_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_postsmooth_iter_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_postsmooth_iter_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_relaxation_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_relaxation_v=dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_polynomial_degree_v")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_polynomial_degree_v = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_strong_threshold_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_strong_threshold_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_truncation_threshold_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_truncation_threshold_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_max_row_sum_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_row_sum_v = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_amli_degree_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_amli_degree_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_nl_amli_krylov_type_v")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_nl_amli_krylov_type_v = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_ILU_levels_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_ILU_levels_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_schwarz_levels_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = FASP_SUCCESS; break; }
            Input->AMG_schwarz_levels_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
                Input->AMG_type_p = CLASSIC_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                Input->AMG_type_p = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                Input->AMG_type_p = UA_AMG;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggregation_type_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggregation_type_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_pair_number_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_pair_number_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_quality_bound_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_quality_bound_p = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_strong_coupled_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_strong_coupled_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_max_aggregation_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_aggregation_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_tentative_smooth_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tentative_smooth_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smooth_filter_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                Input->AMG_smooth_filter_p = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_smooth_filter_p = OFF;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_scaling_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                Input->AMG_coarse_scaling_p = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_coarse_scaling_p = OFF;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_levels_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_levels_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_tol_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tol_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_maxit_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_maxit_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_dof_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarse_dof_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_solver_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarse_solver_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_cycle_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
                Input->AMG_cycle_type_p = V_CYCLE;
            else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
                Input->AMG_cycle_type_p = W_CYCLE;
            else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
                Input->AMG_cycle_type_p = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
                Input->AMG_cycle_type_p = NL_AMLI_CYCLE;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smoother_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
                Input->AMG_smoother_p = SMOOTHER_JACOBI;
            else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
                Input->AMG_smoother_p = SMOOTHER_GS;
            else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
                Input->AMG_smoother_p = SMOOTHER_SGS;
            else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
                Input->AMG_smoother_p = SMOOTHER_CG;
            else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
                Input->AMG_smoother_p = SMOOTHER_SOR;
            else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
                Input->AMG_smoother_p = SMOOTHER_SSOR;
            else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
                Input->AMG_smoother_p = SMOOTHER_GSOR;
            else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
                Input->AMG_smoother_p = SMOOTHER_SGSOR;
            else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
                Input->AMG_smoother_p = SMOOTHER_POLY;
            else if ((strcmp(buffer,"L1_DIAG")==0)||(strcmp(buffer,"l1_diag")==0))
                Input->AMG_smoother_p = SMOOTHER_L1DIAG;
            else
            { status = ERROR_INPUT_PAR; break; }
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smooth_order_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
                Input->AMG_smooth_order_p = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
                Input->AMG_smooth_order_p = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarsening_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarsening_type_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_interpolation_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_interpolation_type_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_level_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_level_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_path_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_path_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_presmooth_iter_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_presmooth_iter_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_postsmooth_iter_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_postsmooth_iter_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_relaxation_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_relaxation_p=dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_polynomial_degree_p")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_polynomial_degree_p = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_strong_threshold_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_strong_threshold_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_truncation_threshold_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_truncation_threshold_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_max_row_sum_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_row_sum_p = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_amli_degree_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_amli_degree_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_nl_amli_krylov_type_p")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_nl_amli_krylov_type_p = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ILU_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_type = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ILU_lfil")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_lfil = ibuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ILU_droptol")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_droptol = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ILU_relax")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_relax = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ILU_permtol")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_permtol = dbuff;
            wall = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"SWZ_mmsize")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->SWZ_mmsize = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"SWZ_maxlvl")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) {status = ERROR_INPUT_PAR; break; }
            Input->SWZ_maxlvl = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"SWZ_type")==0)
        {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->SWZ_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else
        {
            status = ERROR_INPUT_PAR;
            break;
        }		
    }
    
    fclose(fp);
    
    // sanity checks
    status = fasp_ns_param_check(Input);
    
#if DEBUG_MODE
    printf("### DEBUG: Reading input status = %d\n", status);
#endif
    
    // if meet unexpected input, stop the program
    fasp_chkerr(status,"fasp_param_input");
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
