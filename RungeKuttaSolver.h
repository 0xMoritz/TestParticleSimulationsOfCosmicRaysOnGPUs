/*!
 * @file RungeKuttaSolver.h
 * @brief ODE-solver using Runge Kutta method of second or fourth order
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "Engine.h"

/*! \brief Integration Engine working by the Runge Kutta method
 *
 * The Runge Kutta method evalutes the right hand side of the system of linear ODEs at multiple points to
 * increase the order of accuracy. Methods for different accuracies exist. This class implements a second
 * \f$\mathcal{O}(dt^2)\f$ and fourth order accurate \f$\mathcal{O}(dt^4)\f$ method with conventional coefficients.
 */
class RungeKuttaSolver : public Engine
{
private:
	int order;

	/*!
	 * Second order Runge-Kutta > Numerical Recipes Eq. 17.1.2
	 *
	 * \f$ k_1 = dt*f(t_n, y_n) \f$								(1)
	 * \f$ k_2 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1) \f$ 			(2)
	 * \f$ y_{n+1} = y_n + k_2 \f$ 								(3)
	 */
	void RK2(Vec& q, Vec& dqdt, Vec& q_out);

	/*!
	 * Fourth order Runge-Kutta > Numerical Recipes Eq. 17.1.3
	 *
	 * \f$ k_1 = dt*f(t_n, y_n)	\f$ 							(1)
	 * \f$ k_2 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1) \f$ 			(2)
	 * \f$ k_3 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2) \f$ 			(3)
	 * \f$ k_4 = dt*f(t_n + dt, y_n + k_3) \f$ 					(4)
	 * \f$ y_{n+1} = y_n + k_1/6 + k_2/3 + k_3/3 + k_4/6 \f$ 	(5)
	 */
	void RK4(Vec & q, Vec &dqdt, Vec &q_out);

public:
	RungeKuttaSolver(const int particleCount_, const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const FieldGenerator& field, const T omega_, const int order_);
	RungeKuttaSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, const FieldGenerator& field, const T omega_, const int order_);

	/*!
	 * One step of simulation using the Runge Kutta Method of the specified order
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
