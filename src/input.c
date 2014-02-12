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
 */
SHORT fasp_ns_param_check (input_ns_param *inparam)
{
    SHORT status = SUCCESS;
    
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
        || inparam->Schwarz_mmsize<0
        || inparam->Schwarz_maxlvl<0
        || inparam->Schwarz_type<0
        || inparam->AMG_type<=0
        || inparam->AMG_type>3
        || inparam->AMG_cycle_type<=0
        || inparam->AMG_cycle_type>4
        || inparam->AMG_levels<0
        || inparam->AMG_ILU_levels<0
        || inparam->AMG_coarse_dof<=0
        || inparam->AMG_tol<0
        || inparam->AMG_maxit<0
        || inparam->AMG_coarsening_type<=0
        || inparam->AMG_coarsening_type>4
        || inparam->AMG_interpolation_type<0
        || inparam->AMG_interpolation_type>5
        || inparam->AMG_smoother<0
        || inparam->AMG_smoother>20
        || inparam->AMG_strong_threshold<0.0
        || inparam->AMG_strong_threshold>0.9999
        || inparam->AMG_truncation_threshold<0.0
        || inparam->AMG_truncation_threshold>0.9999
        || inparam->AMG_max_row_sum<0.0
        || inparam->AMG_presmooth_iter<0
        || inparam->AMG_postsmooth_iter<0
        || inparam->AMG_amli_degree<0
        || inparam->AMG_aggressive_level<0
        || inparam->AMG_aggressive_path<0
        || inparam->AMG_strong_coupled<0
        || inparam->AMG_max_aggregation<=0
        || inparam->AMG_tentative_smooth<0
        || inparam->AMG_smooth_filter<0
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
 * \date   03/15/2012
 */
void fasp_ns_param_input (char *filenm, 
                          input_ns_param *Input)
{
	char     buffer[500]; // Note: max number of char for each line!
    char   * wall;
	INT      val; 
	SHORT    status = SUCCESS;
    
	// set default input parameters
	fasp_ns_param_input_init(Input);

	// if input file is not specified, use the default values
	if (filenm==NULL) return;
	
	FILE *fp = fopen(filenm,"r");
	if (fp==NULL) {
		printf("### ERROR: Could not open file %s...\n", filenm);
		exit(ERROR_OPEN_FILE);
	}
    
	while ( status == SUCCESS ) {
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
		
        else if (strcmp(buffer,"solver_v_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->solver_v_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"precond_v_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->precond_v_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"itsolver_v_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_v_tol = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"itsolver_v_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_v_maxit = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}

        else if (strcmp(buffer,"itsolver_v_restart")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_v_restart = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
        
        else if (strcmp(buffer,"solver_p_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->solver_p_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"precond_p_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->precond_p_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"itsolver_p_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_p_tol = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"itsolver_p_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_p_maxit = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
        
        else if (strcmp(buffer,"itsolver_p_restart")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->pre_p_restart = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
        
		else if (strcmp(buffer,"AMG_ILU_levels")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_ILU_levels = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		else if (strcmp(buffer,"AMG_schwarz_levels")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = SUCCESS; break; }
			Input->AMG_schwarz_levels = ibuff;
			fgets(buffer,500,fp); // skip rest of line
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
		
		else if (strcmp(buffer,"AMG_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
				Input->AMG_type = CLASSIC_AMG;
			else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
				Input->AMG_type = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
				Input->AMG_type = UA_AMG;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_strong_coupled")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_strong_coupled = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_max_aggregation")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_max_aggregation = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_tentative_smooth")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_tentative_smooth = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_smooth_filter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
				Input->AMG_smooth_filter = ON;
			else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
				Input->AMG_smooth_filter = OFF;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_coarse_scaling")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
				Input->AMG_coarse_scaling = ON;
			else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
				Input->AMG_coarse_scaling = OFF;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_levels")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_levels = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_tol = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_maxit = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_coarse_dof")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_coarse_dof = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_cycle_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
				Input->AMG_cycle_type = V_CYCLE;
			else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
				Input->AMG_cycle_type = W_CYCLE;
			else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
				Input->AMG_cycle_type = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
				Input->AMG_cycle_type = NL_AMLI_CYCLE; 
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_smoother")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
				Input->AMG_smoother = SMOOTHER_JACOBI;
			else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
				Input->AMG_smoother = SMOOTHER_GS;
			else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
				Input->AMG_smoother = SMOOTHER_SGS;
			else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
				Input->AMG_smoother = SMOOTHER_CG;				
			else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
				Input->AMG_smoother = SMOOTHER_SOR;
			else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
				Input->AMG_smoother = SMOOTHER_SSOR;
			else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
				Input->AMG_smoother = SMOOTHER_GSOR;
			else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
				Input->AMG_smoother = SMOOTHER_SGSOR;
			else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
				Input->AMG_smoother = SMOOTHER_POLY;
			else if ((strcmp(buffer,"L1_DIAG")==0)||(strcmp(buffer,"l1_diag")==0))
				Input->AMG_smoother = SMOOTHER_L1DIAG;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_smooth_order")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
                Input->AMG_smooth_order = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
                Input->AMG_smooth_order = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_coarsening_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_coarsening_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_interpolation_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_interpolation_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_aggressive_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_path")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_aggressive_path = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_presmooth_iter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_presmooth_iter = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_postsmooth_iter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_postsmooth_iter = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_relaxation")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_relaxation=dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_polynomial_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_polynomial_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_strong_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_strong_threshold = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_truncation_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_truncation_threshold = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_max_row_sum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_max_row_sum = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_amli_degree")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_amli_degree = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
        
        else if (strcmp(buffer,"AMG_nl_amli_krylov_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_nl_amli_krylov_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_p_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
            Input->AMG_p_type = CLASSIC_AMG;
			else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
            Input->AMG_p_type = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
            Input->AMG_p_type = UA_AMG;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_strong_coupled")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_strong_coupled = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_max_aggregation")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_max_aggregation = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_tentative_smooth")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_tentative_smooth = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_smooth_filter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
            Input->AMG_p_smooth_filter = ON;
			else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
            Input->AMG_p_smooth_filter = OFF;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_coarse_scaling")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
            Input->AMG_p_coarse_scaling = ON;
			else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
            Input->AMG_p_coarse_scaling = OFF;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_levels")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_levels = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_tol = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_maxit = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_coarse_dof")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_coarse_dof = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_cycle_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
            Input->AMG_p_cycle_type = V_CYCLE;
			else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
            Input->AMG_p_cycle_type = W_CYCLE;
			else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
            Input->AMG_p_cycle_type = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
            Input->AMG_p_cycle_type = NL_AMLI_CYCLE;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_smoother")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%s",buffer);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			
			if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
            Input->AMG_p_smoother = SMOOTHER_JACOBI;
			else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
            Input->AMG_p_smoother = SMOOTHER_GS;
			else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
            Input->AMG_p_smoother = SMOOTHER_SGS;
			else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
            Input->AMG_p_smoother = SMOOTHER_CG;
			else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
            Input->AMG_p_smoother = SMOOTHER_SOR;
			else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
            Input->AMG_p_smoother = SMOOTHER_SSOR;
			else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
            Input->AMG_p_smoother = SMOOTHER_GSOR;
			else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
            Input->AMG_p_smoother = SMOOTHER_SGSOR;
			else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
            Input->AMG_p_smoother = SMOOTHER_POLY;
			else if ((strcmp(buffer,"L1_DIAG")==0)||(strcmp(buffer,"l1_diag")==0))
            Input->AMG_p_smoother = SMOOTHER_L1DIAG;
			else
			{ status = ERROR_INPUT_PAR; break; }
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_p_smooth_order")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
            Input->AMG_p_smooth_order = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
            Input->AMG_p_smooth_order = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_p_coarsening_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_coarsening_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_interpolation_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_interpolation_type = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_p_aggressive_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_p_aggressive_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_p_aggressive_path")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_p_aggressive_path = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_p_presmooth_iter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_presmooth_iter = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_postsmooth_iter")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_postsmooth_iter = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_relaxation")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_relaxation=dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
        else if (strcmp(buffer,"AMG_p_polynomial_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_p_polynomial_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
		else if (strcmp(buffer,"AMG_p_strong_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_strong_threshold = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_truncation_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_truncation_threshold = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_max_row_sum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_max_row_sum = dbuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_p_amli_degree")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_amli_degree = ibuff;
			wall = fgets(buffer,500,fp); // skip rest of line
		}
        
        else if (strcmp(buffer,"AMG_p_nl_amli_krylov_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->AMG_p_nl_amli_krylov_type = ibuff;
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
		
        else if (strcmp(buffer,"Schwarz_mmsize")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->Schwarz_mmsize = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"Schwarz_maxlvl")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) {status = ERROR_INPUT_PAR; break; }
			Input->Schwarz_maxlvl = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"Schwarz_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			Input->Schwarz_type = ibuff;
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
