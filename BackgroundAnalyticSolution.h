/*
 * BackgroundAnalyticSolution.h
 *
 *  Created on: Jun 9, 2022
 *      Author: moritz
 */

#pragma once

#include "Engine.h"


/*!
 * Analytic solution of the coupled ODEs as provided by Kuhlen eq.185 ff. for comparison
 * v_perp and v_par are the perpendicular and parallel components of the velocity to the magnetic
 * field
 */
class BackgroundAnalyticSolution : public Engine
{
private:
	T v_perp;
	T v_par;
	T OMEGA;

public:
	/*!
	 * Assuming the initial y-velocity component is zero, and the initial position is 0
	 */
	BackgroundAnalyticSolution(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_);

	/*!
	 * One "Step" of the analytic solutions outputs the analytic solution at the time @t
	 * the @dqdt is redundant and no valid derivative will be returned.
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
