/*!
 * @file BackgroundAnalyticSolution.h
 * @brief Solver like class to calculate analytical solution to background field
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner Geßner
 */

#pragma once

#include "Engine.h"


/*! \brief Engine class to calculate the analytical trajectory of a background field
 *
 * Engine-style pseudo-integrator that uses the analytical solution of a particle's trajectory
 * in a background homogenous magnetic field
 */
class BackgroundAnalyticSolution : public Engine
{
private:
	T v_perp; //!< velocity perpendicular to magnetic field
	T v_par;  //!< velocity parallel to magnetic field
	T OMEGA;  //!< Larmor frequency

public:
	/*!
	 * Assuming the initial y-velocity component is zero, and the initial position is 0. Additionally @field is
	 * assumed to only contain a background component
	 */
	BackgroundAnalyticSolution(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_);

	/*!
	 * One "Step" of the analytic solutions outputs the analytic solution at the time @t
	 * the @dqdt is redundant and no valid derivative will be returned.
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
