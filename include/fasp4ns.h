/*
 * fasp4ns.h
 *
 *
 */

/*! \file fasp4ns.h
 *  \brief Main header file for FASP4NS package
 */ 
#include "messages_ns.h"
#ifndef __FASP4NS_HEADER__/*-- allow multiple inclusions --*/
#define __FASP4NS_HEADER__

/*-----------------------*/
/*--- Data structures ---*/
/*-----------------------*/
/**
 * \brief Parameters passed to the preconditioner for generalized Navier-Stokes problems
 *
 */
typedef struct precond_ns_param {
    
	//! AMG type
	int AMG_type;
	//! print level in AMG preconditioner
	int print_level;
	//! max number of AMG levels
	int max_levels; 
	
} precond_ns_param;

/**
 * \brief Data passed to the preconditioner for generalized Navier-Stokes problems
 *
 */
typedef struct precond_ns_data {
    
	//! size of A, B, and whole matrix
	int colA, colB, col;
	
	double beta;
	
	AMG_data *mgl_data; /**< AMG data for velocity block */
    AMG_data *mgl_data_p; /**< AMG data for presure block */
    
	//! print level in AMG preconditioner
	int print_level;
	//! max number of AMG levels
	int max_levels;
	//! max number of iterations of AMG preconditioner
	int maxit;
	//! tolerance for AMG preconditioner
	double amg_tol;
	//! AMG cycle type
	int cycle_type;
	//! AMG smoother type
	int smoother;
	//! number of presmoothing
	int presmooth_iter;
	//! number of postsmoothing
	int postsmooth_iter;
	//! coarsening type
	int coarsening_type;
	//! relaxation parameter for SOR smoother
	double relaxation;
	//! switch of scaling of coarse grid correction
	int coarse_scaling;
    
    int AMG_tol;
    
	dCSRmat *M; /**< mass matrix for pressure*/
	dvector *diag_M; /**< diagonal of mass matrix M */ 
    dCSRmat *B; /**<  matrix B*/	
	dCSRmat *Bt; /**< matrix of transpose of B*/
	dCSRmat *C; /**<  matrix B*/
	dCSRmat *S;  /**< Schur Complement matrix*/
    dvector *diag_S; /**< diagonal of Schur Complement matrix S */
    dCSRmat *P;  /**< Poisson matrix of pressure*/
    dvector *rp; /**< residual for pressure */
    dvector *sp; /**< sol for pressure */
    ILU_data *LU_S; /**< LU date for schur */
    Schwarz_data *schwarz_S;
	//! temporary work space
	double *w; /**<  temporary work space for other usage */
	double *w1;
    
} precond_ns_data;

/**
 * \struct input_param
 * \brief Input parameters 
 *
 * Input parameters, reading from disk file
 */
typedef struct {
	
	// output flags
	SHORT print_level; /**< print level */
	SHORT output_type; /**< type of output stream */
	
	// problem parameters
	char workdir[256]; /**< working directory for data files */
	INT  problem_num; /**< problem number to solve */
	
	// parameters for iterative solvers
	SHORT itsolver_type; /**< type of iterative solvers */
	SHORT precond_type; /**< type of preconditioner for iterative solvers */
	SHORT stop_type; /**< type of stopping criteria for iterative solvers */
	REAL itsolver_tol; /**< tolerance for iterative linear solver */
	INT itsolver_maxit; /**< maximal number of iterations for iterative solvers */
	INT restart; /**< restart number used in GMRES */
	
	//pamameters for ILU
	SHORT ILU_type; /**< ILU type for decomposition*/
	INT ILU_lfil; /**< level of fill-in */
	REAL ILU_droptol; /**< drop tolerence */
	REAL ILU_relax; /**< add the sum of dropped elements to diagnal element in proportion relax */
	REAL ILU_permtol; /**< permutation tol */
	
	// parameters for AMG
	SHORT AMG_type; /**< Type of AMG */
	SHORT AMG_levels; /**< maximal number of levels */
	SHORT AMG_cycle_type; /**< type of cycle*/
	SHORT AMG_smoother; /**< type of smoother */
	REAL AMG_relaxation; /**< over-relaxation parameter for SOR */
	SHORT AMG_presmooth_iter; /**< number of presmoothing */
	SHORT AMG_postsmooth_iter; /**< number of postsmoothing */
	INT AMG_coarse_dof;	/**< minimal coarsest level dof */
	REAL AMG_tol; /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit; /**< max number of iterations for AMG if used as preconditioner */
	SHORT AMG_ILU_levels; /**< how many levels use ILU smoother */	
	SHORT AMG_coarse_scaling; /**< switch of scaling of the coarse grid correction */
	SHORT AMG_amli_degree; /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type; /**< type of krylov method used by nonlinear AMLI cycle */
	
	// parameters for classical AMG
	SHORT AMG_coarsening_type; /**< coarsening type */
	SHORT AMG_interpolation_type; /**< interpolation type */
	REAL AMG_strong_threshold; /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum; /**< maximal row sum */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled; /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation; /**< max size of each aggregate */
	REAL AMG_tentative_smooth; /**< relaxation parameter for smoothing the tentative prolongation */
	SHORT AMG_smooth_filter; /**< switch for filtered matrix for smoothing the tentative prolongation */
	
} input_ns_param; /**< Input parameters for NS problem */

#endif /* end if for __FASP4NS_HEADER__ */

/* Ene of fasp4ns.h */
