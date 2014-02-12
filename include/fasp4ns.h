/*
 * fasp4ns.h
 *
 *
 */

/*! \file fasp4ns.h
 *  \brief Main header file for FASP4NS package
 */ 
#include "messages_ns.h"
#include "fasp.h"

#ifndef __FASP4NS_HEADER__/*-- allow multiple inclusions --*/
#define __FASP4NS_HEADER__

/*-----------------------*/
/*--- Data structures ---*/
/*-----------------------*/
/**
 * \struct AMG_param
 * \brief Parameters for AMG solver
 *
 * \note This is needed for the AMG solver/preconditioner.
 */
typedef struct {
    AMG_param param_v;
    AMG_param param_p;
} AMG_ns_param;

/*
typedef struct {
	
	//! type of AMG method
	SHORT AMG_type;
	
    //! print level for AMG
	SHORT print_level;
	
    //! max number of iterations of AMG
	INT maxit;
	
    //! stopping tolerance for AMG solver
	REAL tol;
	
	//! max number of levels of AMG
	SHORT max_levels;
	
    //! max coarsest level dof
	INT coarse_dof;
	
    //! type of AMG cycle
	SHORT cycle_type;
	
    //! smoother type
	SHORT smoother;
	
    //! smoother order
	SHORT smooth_order;  // 1: nature order 2: C/F order (both are symmetric)
	
    //! number of presmoothers
	SHORT presmooth_iter;
	
    //! number of postsmoothers
	SHORT postsmooth_iter;
	
    //! relaxation parameter for SOR smoother
	REAL relaxation;
    
    //! degree of the polynomial smoother
    SHORT polynomial_degree;
	
    //! switch of scaling of the coarse grid correction
	SHORT coarse_scaling;
	
    //! degree of the polynomial used by AMLI cycle
	SHORT amli_degree;
	
    //! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
    
    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
	
	//! coarsening type
	SHORT coarsening_type;
	
    //! interpolation type
	SHORT interpolation_type;
	
	//! strong connection threshold for coarsening
	REAL strong_threshold;
	
    //! maximal row sum parameter
	REAL max_row_sum;
	
    //! truncation threshold
	REAL truncation_threshold;
    
    //! number of levels use aggressive coarsening
    INT aggressive_level;
    
    //! numebr of paths use to determin stongly coupled C points
    INT aggressive_path;
	
	//! strong coupled threshold for aggregate
	REAL strong_coupled;
	
    //! max size of each aggregate
	INT max_aggregation;
	
    //! relaxation parameter for smoothing the tentative prolongation
	REAL tentative_smooth;
	
    //! switch for filtered matrix used for smoothing the tentative prolongation
	SHORT smooth_filter;

	//! type of AMG method
	SHORT AMG_p_type;
	
    //! print level for AMG
	SHORT p_print_level;
	
    //! max number of iterations of AMG
	INT p_maxit;
	
    //! stopping tolerance for AMG solver
	REAL p_tol;
	
	//! max number of levels of AMG
	SHORT p_max_levels;
	
    //! max coarsest level dof
	INT p_coarse_dof;
	
    //! type of AMG cycle
	SHORT p_cycle_type;
	
    //! smoother type
	SHORT p_smoother;
	
    //! smoother order
	SHORT p_smooth_order;  // 1: nature order 2: C/F order (both are symmetric)
	
    //! number of presmoothers
	SHORT p_presmooth_iter;
	
    //! number of postsmoothers
	SHORT p_postsmooth_iter;
	
    //! relaxation parameter for SOR smoother
	REAL p_relaxation;
    
    //! degree of the polynomial smoother
    SHORT p_polynomial_degree;
	
    //! switch of scaling of the coarse grid correction
	SHORT p_coarse_scaling;
	
    //! degree of the polynomial used by AMLI cycle
	SHORT p_amli_degree;
	
    //! coefficients of the polynomial used by AMLI cycle
	REAL *p_amli_coef;
    
    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT p_nl_amli_krylov_type;
	
	//! coarsening type
	SHORT p_coarsening_type;
	
    //! interpolation type
	SHORT p_interpolation_type;
	
	//! strong connection threshold for coarsening
	REAL p_strong_threshold;
	
    //! maximal row sum parameter
	REAL p_max_row_sum;
	
    //! truncation threshold
	REAL p_truncation_threshold;
    
    //! number of levels use aggressive coarsening
    INT p_aggressive_level;
    
    //! numebr of paths use to determin stongly coupled C points
    INT p_aggressive_path;
	
	//! strong coupled threshold for aggregate
	REAL p_strong_coupled;
	
    //! max size of each aggregate
	INT p_max_aggregation;
	
    //! relaxation parameter for smoothing the tentative prolongation
	REAL p_tentative_smooth;
	
    //! switch for filtered matrix used for smoothing the tentative prolongation
	SHORT p_smooth_filter;

    //! number of levels use ILU smoother
	SHORT p_ILU_levels;
    //! number of levels use schwarz smoother
	INT p_schwarz_levels;
    
	//! number of levels use ILU smoother
	SHORT ILU_levels;
	
    //! ILU type for smoothing
	SHORT ILU_type;
	
    //! level of fill-in for ILUs and ILUk
	INT ILU_lfil;
	
    //! drop tolerence for ILUt
	REAL ILU_droptol;
	
    //! relaxiation for ILUs
	REAL ILU_relax;
	
    //! permuted if permtol*|a(i,j)| > |a(i,i)|
	REAL ILU_permtol;
    
	//! number of levels use schwarz smoother
	INT schwarz_levels;
	
    //! maximal block size
	INT schwarz_mmsize;
	
    //! maximal levels
	INT schwarz_maxlvl;
	
    //! type of schwarz method
	INT schwarz_type;
	
} AMG_ns_param; *//**< Parameters for AMG */

/**
 * \struct itsolver_param
 * \brief Parameters passed to iterative solvers
 */
typedef struct {
	
	SHORT itsolver_type; /**< solver type: see message.h */
	SHORT precond_type;  /**< preconditioner type: see message.h */
	SHORT stop_type;     /**< stopping criteria type */
	INT   maxit;         /**< max number of iterations */
	REAL  tol;           /**< convergence tolerance */
	INT   restart;       /**< number of steps for restarting: for GMRES etc */
	SHORT print_level;   /**< print level: 0--10 */
    
    SHORT solver_v_type; /**< solver type: see message.h */
	SHORT precond_v_type;  /**< preconditioner type: see message.h */
	INT   pre_v_maxit;         /**< max number of iterations */
	REAL  pre_v_tol;           /**< convergence tolerance */
	INT   pre_v_restart;       /**< number of steps for restarting: for GMRES etc */
	SHORT v_print_level;   /**< print level: 0--10 */
    
    SHORT solver_p_type; /**< solver type: see message.h */
	SHORT precond_p_type;  /**< preconditioner type: see message.h */
	INT   pre_p_maxit;         /**< max number of iterations */
	REAL  pre_p_tol;           /**< convergence tolerance */
	INT   pre_p_restart;       /**< number of steps for restarting: for GMRES etc */
	SHORT p_print_level;   /**< print level: 0--10 */
    
	
} itsolver_ns_param; /**< Parameters for iterative solvers */


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
    AMG_param *param_v;
    AMG_param *param_p;
    itsolver_param *itsolver_param_v;
    itsolver_param *itsolver_param_p;
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
	SHORT print_level;   /**< print level */
	SHORT output_type;   /**< type of output stream */
	
	// problem parameters
    char inifile[256];   /**< ini file name */
	char workdir[256];   /**< working directory for data files */
	INT  problem_num;    /**< problem number to solve */
	
	// parameters for iterative solvers
	SHORT solver_type;   /**< type of iterative solvers */
	SHORT precond_type;  /**< type of preconditioner for iterative solvers */
	SHORT stop_type;     /**< type of stopping criteria for iterative solvers */
	REAL itsolver_tol;   /**< tolerance for iterative linear solver */
	INT itsolver_maxit;  /**< maximal number of iterations for iterative solvers */
	INT restart;         /**< restart number used in GMRES */
	
    //parameters for the block preconditioner
    SHORT solver_v_type;   /**< type of velocity block */
    SHORT precond_v_type;  /**< type of preconditioner for velocity block */
    REAL pre_v_tol;   /**< tolerance for velocity block */
	INT  pre_v_maxit;  /**< maximal number for velocity block */
	INT pre_v_restart;         /**< restart number for velocity block */

    SHORT solver_p_type;   /**< type of pressure block */
    SHORT precond_p_type;  /**< type of preconditioner for pressure block */
    REAL pre_p_tol;   /**< tolerance for pressure block */
	INT  pre_p_maxit;  /**< maximal number for pressure block */
	INT pre_p_restart; /**< restart number for pressure block */
    
	//pamameters for ILU
	SHORT ILU_type;      /**< ILU type for decomposition*/
	INT ILU_lfil;        /**< level of fill-in */
	REAL ILU_droptol;    /**< drop tolerance */
	REAL ILU_relax;      /**< scaling factor: add the sum of dropped entries to diagnal */
	REAL ILU_permtol;    /**< permutation tolerance */
    
    // parameter for Schwarz
	INT Schwarz_mmsize;  /**< maximal block size */
	INT Schwarz_maxlvl;  /**< maximal levels */
	INT Schwarz_type;    /**< type of schwarz method */
	
	// parameters for AMG
	SHORT AMG_type;                /**< Type of AMG */
	SHORT AMG_levels;              /**< maximal number of levels */
	SHORT AMG_cycle_type;          /**< type of cycle */
	SHORT AMG_smoother;            /**< type of smoother */
	SHORT AMG_smooth_order;        /**< order for smoothers */
	REAL AMG_relaxation;           /**< over-relaxation parameter for SOR */
    SHORT AMG_polynomial_degree;   /**< degree of the polynomial smoother */
	SHORT AMG_presmooth_iter;      /**< number of presmoothing */
	SHORT AMG_postsmooth_iter;     /**< number of postsmoothing */
	INT AMG_coarse_dof;	           /**< minimal coarsest level dof */
	REAL AMG_tol;                  /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit;                 /**< number of iterations for AMG used as preconditioner */
	SHORT AMG_ILU_levels;          /**< how many levels use ILU smoother */
	SHORT AMG_coarse_scaling;      /**< switch of scaling of the coarse grid correction */
	SHORT AMG_amli_degree;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type; /**< type of krylov method used by nonlinear AMLI cycle */
    INT AMG_schwarz_levels;        /**< number of levels use schwarz smoother */
	
	// parameters for classical AMG
	SHORT AMG_coarsening_type;     /**< coarsening type */
	SHORT AMG_interpolation_type;  /**< interpolation type */
	REAL AMG_strong_threshold;     /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum;          /**< maximal row sum */
    INT AMG_aggressive_level;      /**< number of levels use aggressive coarsening */
    INT AMG_aggressive_path;       /**< number of paths used to determine strongly coupled C-set */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled;       /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation;       /**< max size of each aggregate */
	REAL AMG_tentative_smooth;     /**< relaxation factor for smoothing the tentative prolongation */
	SHORT AMG_smooth_filter;       /**< use filterfor smoothing the tentative prolongation or not */
    
    // parameters for AMG of pressure part
	SHORT AMG_p_type;                /**< Type of AMG */
	SHORT AMG_p_levels;              /**< maximal number of levels */
	SHORT AMG_p_cycle_type;          /**< type of cycle */
	SHORT AMG_p_smoother;            /**< type of smoother */
	SHORT AMG_p_smooth_order;        /**< order for smoothers */
	REAL AMG_p_relaxation;           /**< over-relaxation parameter for SOR */
    SHORT AMG_p_polynomial_degree;   /**< degree of the polynomial smoother */
	SHORT AMG_p_presmooth_iter;      /**< number of presmoothing */
	SHORT AMG_p_postsmooth_iter;     /**< number of postsmoothing */
	INT AMG_p_coarse_dof;	           /**< minimal coarsest level dof */
	REAL AMG_p_tol;                  /**< tolerance for AMG if used as preconditioner */
	INT AMG_p_maxit;                 /**< number of iterations for AMG used as preconditioner */
	SHORT AMG_p_ILU_levels;          /**< how many levels use ILU smoother */
	SHORT AMG_p_coarse_scaling;      /**< switch of scaling of the coarse grid correction */
	SHORT AMG_p_amli_degree;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_p_nl_amli_krylov_type; /**< type of krylov method used by nonlinear AMLI cycle */
    INT AMG_p_schwarz_levels;        /**< number of levels use schwarz smoother */
	
	// parameters for classical AMG
	SHORT AMG_p_coarsening_type;     /**< coarsening type */
	SHORT AMG_p_interpolation_type;  /**< interpolation type */
	REAL AMG_p_strong_threshold;     /**< strong threshold for coarsening */
	REAL AMG_p_truncation_threshold; /**< truncation factor for interpolation */
	REAL AMG_p_max_row_sum;          /**< maximal row sum */
    INT AMG_p_aggressive_level;      /**< number of levels use aggressive coarsening */
    INT AMG_p_aggressive_path;       /**< number of paths used to determine strongly coupled C-set */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_p_strong_coupled;       /**< strong coupled threshold for aggregate */
	INT AMG_p_max_aggregation;       /**< max size of each aggregate */
	REAL AMG_p_tentative_smooth;     /**< relaxation factor for smoothing the tentative prolongation */
	SHORT AMG_p_smooth_filter;       /**< use filterfor smoothing the tentative prolongation or not */
    
} input_ns_param; /**< Input parameters */

#endif /* end if for __FASP4NS_HEADER__ */

/* Ene of fasp4ns.h */
