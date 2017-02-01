/*! \file  fasp4ns_functs.h
 *
 *  \brief Function decoration for the FASP package
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2008--2017 by the FASP team. All rights reserved.                
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  \warning DO NOT EDIT!!! This file is automatically generated!
 */ 

#include "fasp.h" 
#include "fasp_block.h" 

/*-------- In file: AuxInput.c --------*/

SHORT fasp_ns_param_check (const input_ns_param *inparam);

void fasp_ns_param_input (char *filenm,
                          input_ns_param *Input);


/*-------- In file: AuxParam.c --------*/

void fasp_ns_param_init (input_ns_param *inparam,
                         itsolver_ns_param *itsparam,
                         AMG_ns_param *amgparam,
                         ILU_param *iluparam,
                         SWZ_param *swzparam);

void fasp_ns_param_input_init (input_ns_param *inparam);

void fasp_ns_param_amg_init (AMG_ns_param *amgparam);

void fasp_ns_param_amg_set (AMG_ns_param *param,
                            input_ns_param *inparam);

void fasp_ns_param_solver_init(itsolver_ns_param *itsparam);

void fasp_ns_param_solver_set (itsolver_ns_param *itsparam,
                               input_ns_param *inparam);

void fasp_ns_param_ilu_set (ILU_param *iluparam,
                            input_ns_param *inparam);

void fasp_ns_param_swz_set (SWZ_param *swzparam,
                            input_ns_param *inparam);


/*-------- In file: BlaIO.c --------*/

void fasp_dblc_read (char *fileA,
                     char *fileB,
                     char *fileC,
                     char *filerhs,
                     dBLCmat *A,
                     dvector *r);

void fasp_dblc_read_1 (char *fileA,
                       char *fileB,
                       char *fileC,
                       char *filerhs,
                       dBLCmat *A,
                       dvector *r);

void fasp_dblc_read_ruth (char *fileA,
                          char *fileB,
                          char *fileC,
                          char *fileD,
                          char *filerhs,
                          char *filex0,
                          dBLCmat *A,
                          dvector *r,
                          dvector *x0);


/*-------- In file: PreNavierStokes.c --------*/

void fasp_precond_ns_bdiag (REAL *r,
                            REAL *z,
                            void *data);

void fasp_precond_ns_low_btri (REAL *r,
                               REAL *z,
                               void *data);

void fasp_precond_ns_up_btri (REAL *r,
                              REAL *z,
                              void *data);

void fasp_precond_ns_blu (REAL *r,
                          REAL *z,
                          void *data);

void fasp_precond_ns_simple (REAL *r,
                             REAL *z,
                             void *data);

void fasp_precond_ns_simpler (REAL *r,
                              REAL *z,
                              void *data);

void fasp_precond_ns_uzawa (REAL *r,
                            REAL *z,
                            void *data);

void fasp_precond_ns_projection (REAL *r,
                                 REAL *z,
                                 void *data);

void fasp_precond_ns_DGS (REAL *r,
                          REAL *z,
                          void *data);

void fasp_precond_ns_LSCDGS (REAL *r,
                             REAL *z,
                             void *data);

void fasp_precond_ns_sym_btri (REAL *r,
                               REAL *z,
                               void *data);

void fasp_precond_ns_lsc (REAL *r,
                          REAL *z,
                          void *data);


/*-------- In file: PrePNPStokes.c --------*/

void fasp_precond_pnp_stokes_diag (REAL *r,
                                   REAL *z,
                                   void *data);

void fasp_precond_pnp_stokes_lower (REAL *r,
                                    REAL *z,
                                    void *data);

void fasp_precond_pnp_stokes_upper (REAL *r,
                                    REAL *z,
                                    void *data);

void fasp_precond_pnp_stokes_diag_inexact (REAL *r,
                                           REAL *z,
                                           void *data);

void fasp_precond_pnp_stokes_lower_inexact (REAL *r,
                                            REAL *z,
                                            void *data);

void fasp_precond_pnp_stokes_upper_inexact (REAL *r,
                                            REAL *z,
                                            void *data);


/*-------- In file: SolNavierStokes.c --------*/

SHORT fasp_ns_solver_itsolver (dBLCmat *A,
                               dvector *b,
                               dvector *x,
                               precond *prec,
                               itsolver_ns_param *itparam);

SHORT fasp_solver_dblc_krylov_navier_stokes (dBLCmat *Mat,
                                             dvector *b,
                                             dvector *x,
                                             itsolver_ns_param *itparam,
                                             AMG_ns_param *amgnsparam,
                                             ILU_param *iluparam,
                                             SWZ_param *schparam);

SHORT fasp_solver_dblc_krylov_navier_stokes_with_pressure_mass (dBLCmat *Mat,
                                                                dvector *b,
                                                                dvector *x,
                                                                itsolver_ns_param *itparam,
                                                                AMG_ns_param *amgnsparam,
                                                                ILU_param *iluparam,
                                                                SWZ_param *schparam,
                                                                dCSRmat *Mp);

SHORT fasp_solver_dblc_krylov_navier_stokes_schur_complement_with_pressure_mass (dBLCmat *Mat,
                                                                                 dvector *b,
                                                                                 dvector *x,
                                                                                 itsolver_ns_param *itparam,
                                                                                 AMG_ns_param *amgnsparam,
                                                                                 ILU_param *iluparam,
                                                                                 SWZ_param *schparam,
                                                                                 dCSRmat *Mp);


/*-------- In file: SolPNPStokes.c --------*/

INT fasp_solver_dblc_krylov_pnp_stokes (dBLCmat *A,
                                        dvector *b,
                                        dvector *x,
                                        ITS_param *itparam,
                                        ITS_param *itparam_pnp,
                                        AMG_param *amgparam_pnp,
                                        itsolver_ns_param *itparam_stokes,
                                        AMG_ns_param *amgparam_stokes,
                                        const int num_velocity,
                                        const int num_pressure);


/*-------- In file: SolWrapper.c --------*/

void fasp_fwrapper_krylov_navier_stokes_ (INT *nA,
                                          INT *nnzA,
                                          INT *ia,
                                          INT *ja,
                                          REAL *aval,
                                          INT *nB,
                                          INT *nnzB,
                                          INT *ib,
                                          INT *jb,
                                          REAL *bval,
                                          INT *nC,
                                          INT *nnzC,
                                          INT *ic,
                                          INT *jc,
                                          REAL *cval,
                                          REAL *b,
                                          REAL *u);

void fasp_fwrapper_krylov_navier_stokes_nsym_ (INT *nA,
                                               INT *nnzA,
                                               INT *ia,
                                               INT *ja,
                                               REAL *aval,
                                               INT *nB,
                                               INT *mB,
                                               INT *nnzB,
                                               INT *ib,
                                               INT *jb,
                                               REAL *bval,
                                               INT *nC,
                                               INT *mC,
                                               INT *nnzC,
                                               INT *ic,
                                               INT *jc,
                                               REAL *cval,
                                               REAL *b,
                                               REAL *u);

/* Ene of fasp4ns_functs.h */
