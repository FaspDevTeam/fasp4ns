/*
 * message_ns.h
 *
 */

/*! \file message_ns.h
 *  \brief Main header file for FASP4NS package
 */ 
#ifndef __FASP_NS_MESSAGES__		/*-- allow multiple inclusions --*/
#define __FASP_NS_MESSAGES__

/** 
 * \brief Definition of iterative solver types
 */
#define SOLVER_MinRes       3  /**< Minimal Residual */
#define SOLVER_GMRES        4  /**< GMRES Method */
#define SOLVER_FGMRES       5  /**< FGMRES Method */
#define SOLVER_GCR          6  /**< GCR Method */

/** 
 * \brief Definition of preconditioner types
 */
#define PREC_NULL           0  /**< with no preconditioner */
#define PREC_DIAG           1  /**< with block diagonal preconditioner */
#define PREC_LOW_TRI        2  /**< with lower block triangle preconditioner */
#define PREC_UP_TRI         3  /**< with upper block triangle preconditioner */
#define PREC_SYM_TRI        4 /**< with symmetric block triangle preconditioner */


#endif
/* Ene of message_ns.h */
