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
//#define SOLVER_GCR          6  /**< GCR Method */

/** 
 * \brief Definition of preconditioner types
 */
#define PREC_NULL           0  /**< with no preconditioner */
#define PREC_DIAG           1  /**< with block diagonal preconditioner */
#define PREC_LOW_TRI        2  /**< with lower block triangle preconditioner */
#define PREC_UP_TRI         3  /**< with upper block triangle preconditioner */
#define PREC_SYM_TRI        4 /**< with symmetric block triangle preconditioner */

#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"
#define COLOR_RESET   "\x1b[0m"


#endif
/* Ene of message_ns.h */
