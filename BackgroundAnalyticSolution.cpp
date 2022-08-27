/*
 * BackgroundAnalyticSolution.cpp
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner Geßner
 */

#include "BackgroundAnalyticSolution.h"

using namespace std;


BackgroundAnalyticSolution::BackgroundAnalyticSolution(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, q0_, field, omega_)
{
	assert(q0[0] == 0 && q0[1] == 0 && q0[2] == 0 && q0[4] == 0);
	v_perp = q0[3];
	v_par = q0[5];
	T B = field.BField(q0)[2]; // z component of magnetic field, should be the same everywhere if the field is homogenous
	assert(field.BField(q0)[0] == 0 && field.BField(q0)[1] == 0);
	OMEGA = B/R;

	name = "Analytical calculation";
}

void BackgroundAnalyticSolution::Step(const T t, Vec& q, Vec& dqdt, Vec& q_out)
{
	q_out[0] = v_perp/OMEGA*sin(OMEGA*t);
	q_out[1] = v_perp/OMEGA*(cos(OMEGA*t) - 1);
	q_out[2] = v_par*t;
	// the speeds are not calculated because they are not needed for integration
}
