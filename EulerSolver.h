/*
 * EulerSolver.h
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "Engine.h"

/*!
 * Integration solver working by the Euler method, error in  First order of dt: \f$  \mathcal{O}(dt^1)
 */
class EulerSolver : public Engine
{
private:

public:
	EulerSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, T dt, T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_);

	/*!
	 * One Euler integration step of the coupled ODE of a particle moving in a background magnetic field
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
