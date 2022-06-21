/*
 * EulerSolver.cpp
 *
 *  Created on: Jun 9, 2022
 *      Author: moritz
 */

#include "EulerSolver.h"

using namespace std;


EulerSolver::EulerSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, T dt, T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, q0_, field, omega_)
{
	name = "Euler Simulation";
}

void EulerSolver::Step(const T t, Vec& q, Vec& dqdt, Vec& q_out)
{
	Scale(dqdt, dt); VecAdd(q, dqdt); // q_out = q + dt*dqdt
	Derivatives(q, dqdt, field.BField(q)); // Calculate the derivative for next iteration
	q_out = q;
}



