/*
 * fasp4ns.h
 *
 *
 */

/*! \file fasp4ns.h
 *  \brief Main header file for FASP4NS package
 *
 * \note: modified by Xiaozhe Hu on Feb. 21, 2014
 * \note: modified by Xiaozhe Hu on May. 27, 2014
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
    /* AMG param for the whole system                        */
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
    
    /*---------------------------------------------------*/
    /* AMG parameters for the velocity block */
    /*---------------------------------------------------*/
    AMG_param param_v;
    
    /*---------------------------------------------------*/
    /* AMG parameters for the pressure block  */
    /*---------------------------------------------------*/
    AMG_param param_p;
    
} AMG_ns_param;

/**
 * \struct AMG_ns_data
 * \brief Data for AMG solvers for Navier-Stokes problems
 *
 * \note This is needed for the AMG solver/preconditioner for Navier-Stokes problems
 */
typedef struct {
    
    /* Level information */
    
    //! max number of levels
    SHORT max_levels;
    
    //! number of levels in use <= max_levels
    SHORT num_levels;
    
    /* Problem information */
    
    //! pointer to the matrix at level level_num
    block_dCSRmat A;
    
    //! restriction operator at level level_num
    block_dCSRmat R;
    
    //! prolongation operator at level level_num
    block_dCSRmat P;
    
    //! pointer to the right-hand side at level level_num
    dvector b;
    
    //! pointer to the iterative solution at level level_num
    dvector x;
    
    //! cycle type
    INT cycle_type;
    
    /* Coarsest grid information */
    
    //! pointer to the numerical factorization from UMFPACK
    void *Numeric;
    
    //! data for Intel MKL PARDISO
    Pardiso_data pdata;
    
    //! data for MUMPS
    Mumps_data mumps;
    
    /* Extra information */
    
    //! Temporary work space
    dvector w;
    
    //! SA information
    
    //! dimension of the near kernel for SAMG for velocity
    INT near_kernel_dim_v;
    
    //! dimension of the near kernel for SAMG for pressure
    INT near_kernel_dim_p;
    
    //! basis of near kernel space for SAMG for velocity
    REAL **near_kernel_basis_v;
    
    //! basis of near kernel space for SAMG for pressure
    REAL **near_kernel_basis_p;
    
    // Smoother order information
    
    //! pointer to the CF marker for velocity at level level_num
    ivector cfmark_v;
    
    //! pointer to the CF marker for pressure at level level_num
    ivector cfmark_p;
    
    // Advanced Smoother information
    
    //! number of levels use ILU smoother for velocity
    INT ILU_levels_v;
    
    //! number of levels use ILU smoother for pressure
    INT ILU_levels_p;
    
    //! ILU matrix for ILU smoother for velocity
    ILU_data LU_v;
    
    //! ILU matrix for ILU smoother for pressure
    ILU_data LU_p;
    
    //! number of levels use Schwarz smoother for velocity
    INT Schwarz_levels_v;
    
    //! number of levels use Schwarz smoother for pressure
    INT Schwarz_levels_p;
    
    //! data of Schwarz smoother for velocity
    Schwarz_data Schwarz_v;
    
    //! data of Schwarz smoother for pressure
    Schwarz_data Schwarz_p;
    
} AMG_ns_data; /**< Data for AMG */

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
	AMG_data *mgl_data_v;
    
    /*---------------------------------------------------*/
    /* AMG data for the pressure block  */
    /*---------------------------------------------------*/
    AMG_data *mgl_data_p;
    
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
    /* AMG data for the whole system  */
    /*---------------------------------------------------*/
    AMG_ns_data *mgl_ns_data; 
    
    /*---------------------------------------------------*/
    /* AMG param for the whole system                        */
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
    /* ILU data for the pressure block  */
    /*---------------------------------------------------*/
    ILU_data *ILU_p; /**< ILU data for presure block */
    
    /*---------------------------------------------------*/
    /* Extra data */
    /*---------------------------------------------------*/
	dCSRmat *M; /**< mass matrix for pressure*/
	dvector *diag_M; /**< diagonal of mass matrix M */
    dvector *diag_A; /**< diagonal of velocity block A */
    dCSRmat *B; /**<  matrix B*/	
	dCSRmat *Bt; /**< matrix of transpose of B*/
	dCSRmat *C; /**<  matrix C*/
	dCSRmat *S;  /**< Schur Complement matrix*/
    dCSRmat *BABt; /**<  matrix BABt*/
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
    INT AMG_coarse_solver_v;            /**< coarse grid solver */
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
    
    // parameters for aggregation based AMG
    SHORT AMG_aggregation_type_v;    /**< aggregation type */
    INT AMG_pair_number_v;           /**< number of pairs in matching algorithm */
    REAL AMG_quality_bound_v;        /**< threshold for pair wise aggregation */
    INT AMG_max_aggregation_v;       /**< max size of each aggregate */

	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled_v;       /**< strong coupled threshold for aggregate */
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
    INT AMG_coarse_solver_p;            /**< coarse grid solver */
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
    
    // parameters for aggregation based AMG
    SHORT AMG_aggregation_type_p;    /**< aggregation type */
    INT AMG_pair_number_p;           /**< number of pairs in matching algorithm */
    REAL AMG_quality_bound_p;        /**< threshold for pair wise aggregation */
    INT AMG_max_aggregation_p;       /**< max size of each aggregate */

	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled_p;       /**< strong coupled threshold for aggregate */
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

/**
 * \brief Data passed to the preconditioner for block preconditioning for block_dCSRmat format
 *
 * This is needed for the block preconditioner for pnp+stokes system.
 */
typedef struct {
    
    /*-------------------------------------*/
    /* Basic data for block preconditioner */
    /*-------------------------------------*/
    block_dCSRmat *Abcsr; /**< problem data, the blocks */
    
    dCSRmat *A_pnp_csr;      /**< data for pnp diagonal block in csr format*/
    dBSRmat *A_pnp_bsr;      /**< data for pnp diagonal block in bsr format*/
    
    dCSRmat *A_stokes_csr;      /**< data for pnp diagonal block in csr format*/
    block_dCSRmat *A_stokes_bcsr;      /**< data for pnp diagonal block in bsr format*/
    
    dvector r;            /**< temp work space */
    
    /*------------------------------*/
    /* Data for the diagonal blocks */
    /*------------------------------*/
    /*--- solve by direct solver ---*/
    void **LU_diag;       /**< LU decomposition for the diagonal blocks (for UMFpack) */
    
    /*---  solve by inexact solver ---*/
    ILU_data *ILU_pnp;
    dvector *diag_pnp;
    precond_data_bsr *precdata_pnp;          /**< data for pnp diagonal block */
    void (*pnp_fct)(REAL *, REAL *, void *);
    
    precond_ns_data  *precdata_stokes;      /**< data for stokes diagonal block */
    void (*stokes_fct)(REAL *, REAL *, void *);
    
} precond_pnp_stokes_data; /**< Precond data for block matrices */


#endif /* end if for __FASP4NS_HEADER__ */

/* Ene of fasp4ns.h */
