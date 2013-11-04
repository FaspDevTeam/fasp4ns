/*! \file basis.c
 *  \brief Basis functions and problem information
 */
 
#include <stdio.h>
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn double f(double x, double y)
 * \brief load x direction of f of 2D Stokes equation.
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value
 *
 * \note Right hand side when the true solution.
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double f_u(double x, double y)
{
	return -2*cos(x)*sin(y);
}

/**
 * \fn double f(double x, double y)
 * \brief load y direction of f of 2D Stokes equation.
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value
 *
 * \note Right hand side when the true solution.
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double f_v(double x, double y)
{
	return -2*cos(x)*sin(y);
}

/**
 * \fn double f(double x, double y)
 * \brief load g of 2D Stokes equation.
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value
 *
 * \note Right hand side when the true solution.
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double g(double x, double y)
{
	return -2*cos(x)*sin(y);
}

/**
 * \fn double u(double x, double y)
 * \brief true solution u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value 
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double u(double x, double y)
{
	return -cos(x)*sin(y);
}

/**
 * \fn double u(double x, double y)
 * \brief true solution v
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value 
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double v(double x, double y)
{
	return -cos(x)*sin(y);
}

/**
 * \fn double u(double x, double y)
 * \brief true solution p
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value 
 *
 * \author Lu Wang
 * \date 11/30/2011
 */
double p(double x, double y)
{
	return 0;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
