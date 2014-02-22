/*
 * fasp4ns.h
 *
 *
 */

/*! \file fasp4ns.h
 *  \brief Main header file for FASP4NS package
 *
 * \note: modified by Xiaozhe Hu on Feb. 21, 2014
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
    /*---------------------------------------------------*/
    /* AMG parameters for the velocity block */
    /*---------------------------------------------------*/
    AMG_param param_v;
    
    /*---------------------------------------------------*/
    /* AMG parameters for the pressure block  */
    /*---------------------------------------------------*/
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
	
    /*---------------------------------------------------*/
    /* iteartive solver parameters for the whole system  */
    /*---------------------------------------------------*/
	SHORT itsolver_type; /**< solver type: see message.h */
	SHORT precond_type;  /**< preconditioner type: see message.h */
	SHORT stop_type;     /**< stopping criteria type */
	INT   maxit;         /**< max number of iterations */
	REAL  tol;           /**< convergence tolerance */
	INT   restart;       /**< number of steps for restarting: for GMRES etc */
	SHORT print_level;   /**< print level: 0--10 */
    
    /*------------------------------------------------------*/
    /* iteartive solver parameters for the velocity block  */
    /*------------------------------------------------------*/
    SHORT itsolver_type_v; /**< solver type: see message.h */
	SHORT precond_type_v;  /**< preconditioner type: see message.h */
	INT   pre_maxit_v;         /**< max number of iterations */
	REAL  pre_tol_v;           /**< convergence tolerance */
	INT   pre_restart_v;       /**< number of steps for restarting: for GMRES etc */
	SHORT print_level_v;   /**< print level: 0--10 */
    
    /*-----------------------------------------------------*/
    /* iteartive solver parameters for the pressure block  */
    /*-----------------------------------------------------*/
    SHORT itsolver_type_p; /**< solver type: see message.h */
	SHORT precond_type_p;  /**< preconditioner type: see message.h */
	INT   pre_maxit_p;         /**< max number of iterations */
	REAL  pre_tol_p;           /**< convergence tolerance */
	INT   pre_restart_p;       /**< number of steps for restarting: for GMRES etc */
	SHORT print_level_p;   /**< print level: 0--10 */
    
	
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
	
	double beta;  /** what is this? -- Xiaozhe */
	
    /*---------------------------------------------------*/
    /* AMG data for the velocity block  */
    /*---------------------------------------------------*/
	AMG_data *mgl_data_v; /**< AMG data for velocity block */
    
    /*---------------------------------------------------*/
    /* AMG data for the pressure block  */
    /*---------------------------------------------------*/
    AMG_data *mgl_data_p; /**< AMG data for presure block */
    
    /*---------------------------------------------------*/
    /* AMG parameters for the velocity block */
    /*---------------------------------------------------*/
    AMG_param *param_v;
    
    /*---------------------------------------------------*/
    /* AMG paramters for the pressure block  */
    /*---------------------------------------------------*/
    AMG_param *param_p;
    
    /*---------------------------------------------------*/
    /* Iterative solver paramters for the velocity block  */
    /*---------------------------------------------------*/
    itsolver_param *itsolver_param_v;
    
    /*---------------------------------------------------*/
    /* Iterative solver parameters for the pressure block */
    /*---------------------------------------------------*/
    itsolver_param *itsolver_param_p;
    
    /*---------------------------------------------------*/
    /* AMG param  (Do we need this anymore? --Xiaozhe) */
    /*---------------------------------------------------*/
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
    
    /*---------------------------------------------------*/
    /* Extra data */
    /*---------------------------------------------------*/
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
	
    /*---------------------------------------------------*/
	// parameters for iterative solvers
    /*---------------------------------------------------*/
	SHORT solver_type;   /**< type of iterative solvers */
	SHORT precond_type;  /**< type of preconditioner for iterative solvers */
	SHORT stop_type;     /**< type of stopping criteria for iterative solvers */
	REAL itsolver_tol;   /**< tolerance for iterative linear solver */
	INT itsolver_maxit;  /**< maximal number of iterations for iterative solvers */
	INT restart;         /**< restart number used in GMRES */
	
    /*---------------------------------------------------*/
    //parameters for the velocity block
    /*---------------------------------------------------*/
    SHORT itsolver_type_v;   /**< type of velocity block */
    SHORT precond_type_v;  /**< type of preconditioner for velocity block */
    REAL pre_tol_v;   /**< tolerance for velocity block */
	INT  pre_maxit_v;  /**< maximal number for velocity block */
	INT pre_restart_v;         /**< restart number for velocity block */

    /*---------------------------------------------------*/
	// parameters for AMG for the velocity block
    /*---------------------------------------------------*/
	SHORT AMG_type_v;                /**< Type of AMG */
	SHORT AMG_levels_v;              /**< maximal number of levels */
	SHORT AMG_cycle_type_v;          /**< type of cycle */
	SHORT AMG_smoother_v;            /**< type of smoother */
	SHORT AMG_smooth_order_v;        /**< order for smoothers */
	REAL AMG_relaxation_v;           /**< over-relaxation parameter for SOR */
    SHORT AMG_polynomial_degree_v;   /**< degree of the polynomial smoother */
	SHORT AMG_presmooth_iter_v;      /**< number of presmoothing */
	SHORT AMG_postsmooth_iter_v;     /**< number of postsmoothing */
	INT AMG_coarse_dof_v;	           /**< minimal coarsest level dof */
	REAL AMG_tol_v;                  /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit_v;                 /**< number of iterations for AMG used as preconditioner */
	SHORT AMG_ILU_levels_v;          /**< how many levels use ILU smoother */
	SHORT AMG_coarse_scaling_v;      /**< switch of scaling of the coarse grid correction */
	SHORT AMG_amli_degree_v;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type_v; /**< type of krylov method used by nonlinear AMLI cycle */
    INT AMG_schwarz_levels_v;        /**< number of levels use schwarz smoother */
	
	// parameters for classical AMG
	SHORT AMG_coarsening_type_v;     /**< coarsening type */
	SHORT AMG_interpolation_type_v;  /**< interpolation type */
	REAL AMG_strong_threshold_v;     /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold_v; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum_v;          /**< maximal row sum */
    INT AMG_aggressive_level_v;      /**< number of levels use aggressive coarsening */
    INT AMG_aggressive_path_v;       /**< number of paths used to determine strongly coupled C-set */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled_v;       /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation_v;       /**< max size of each aggregate */
	REAL AMG_tentative_smooth_v;     /**< relaxation factor for smoothing the tentative prolongation */
	SHORT AMG_smooth_filter_v;       /**< use filterfor smoothing the tentative prolongation or not */
    
    /*---------------------------------------------------*/
    //parameters for the pressure block
    /*---------------------------------------------------*/
    SHORT itsolver_type_p;   /**< type of pressure block */
    SHORT precond_type_p;  /**< type of preconditioner for pressure block */
    REAL pre_tol_p;   /**< tolerance for pressure block */
	INT  pre_maxit_p;  /**< maximal number for pressure block */
	INT pre_restart_p; /**< restart number for pressure block */
    
    /*---------------------------------------------------*/
	// parameters for AMG for the pressure block
    /*---------------------------------------------------*/
	SHORT AMG_type_p;                /**< Type of AMG */
	SHORT AMG_levels_p;              /**< maximal number of levels */
	SHORT AMG_cycle_type_p;          /**< type of cycle */
	SHORT AMG_smoother_p;            /**< type of smoother */
	SHORT AMG_smooth_order_p;        /**< order for smoothers */
	REAL AMG_relaxation_p;           /**< over-relaxation parameter for SOR */
    SHORT AMG_polynomial_degree_p;   /**< degree of the polynomial smoother */
	SHORT AMG_presmooth_iter_p;      /**< number of presmoothing */
	SHORT AMG_postsmooth_iter_p;     /**< number of postsmoothing */
	INT AMG_coarse_dof_p;	           /**< minimal coarsest level dof */
	REAL AMG_tol_p;                  /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit_p;                 /**< number of iterations for AMG used as preconditioner */
	SHORT AMG_ILU_levels_p;          /**< how many levels use ILU smoother */
	SHORT AMG_coarse_scaling_p;      /**< switch of scaling of the coarse grid correction */
	SHORT AMG_amli_degree_p;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type_p; /**< type of krylov method used by nonlinear AMLI cycle */
    INT AMG_schwarz_levels_p;        /**< number of levels use schwarz smoother */
	
	// parameters for classical AMG
	SHORT AMG_coarsening_type_p;     /**< coarsening type */
	SHORT AMG_interpolation_type_p;  /**< interpolation type */
	REAL AMG_strong_threshold_p;     /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold_p; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum_p;          /**< maximal row sum */
    INT AMG_aggressive_level_p;      /**< number of levels use aggressive coarsening */
    INT AMG_aggressive_path_p;       /**< number of paths used to determine strongly coupled C-set */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled_p;       /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation_p;       /**< max size of each aggregate */
	REAL AMG_tentative_smooth_p;     /**< relaxation factor for smoothing the tentative prolongation */
	SHORT AMG_smooth_filter_p;       /**< use filterfor smoothing the tentative prolongation or not */
    
    /*---------------------------------------------------*/
	//pamameters for ILU
    /*---------------------------------------------------*/
	SHORT ILU_type;      /**< ILU type for decomposition*/
	INT ILU_lfil;        /**< level of fill-in */
	REAL ILU_droptol;    /**< drop tolerance */
	REAL ILU_relax;      /**< scaling factor: add the sum of dropped entries to diagnal */
	REAL ILU_permtol;    /**< permutation tolerance */
    
    /*---------------------------------------------------*/
    // parameter for Schwarz
    /*---------------------------------------------------*/
	INT Schwarz_mmsize;  /**< maximal block size */
	INT Schwarz_maxlvl;  /**< maximal levels */
	INT Schwarz_type;    /**< type of schwarz method */
    
} input_ns_param; /**< Input parameters */

#endif /* end if for __FASP4NS_HEADER__ */

/* Ene of fasp4ns.h */
