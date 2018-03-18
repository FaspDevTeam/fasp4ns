/*! \file  fasp4ns_const.h
 *
 *  \brief Constants used in the FASP4NS package
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2008--2018 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 */

#ifndef __FASP_NS_CONST__ /*-- allow multiple inclusions --*/
#define __FASP_NS_CONST__

/** 
 * \brief Definition of preconditioner types
 */
#define PREC_NULL           0  /**< with no preconditioner */
#define PREC_DIAG           1  /**< with block diagonal preconditioner */
#define PREC_LOW_TRI        2  /**< with lower block triangle preconditioner */
#define PREC_UP_TRI         3  /**< with upper block triangle preconditioner */
#define PREC_SYM_TRI        4 /**< with symmetric block triangle preconditioner */

/**
 * \brief Definition of colors
 */
#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"
#define COLOR_RESET   "\x1b[0m"

#endif

/* Ene of fasp4ns_const.h */
