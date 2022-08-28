/*!
 * @file EulerSolver.h
 * @brief ODE-solver using the Euler method
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
	/*!
	 * Constructor for the ODE solver using the Euler method
	 * \param totalSimulationTime_ total time of the simulation
	 * \param outputPoints_ points written to the output file
	 * \param stepsPerOutput_ integration steps per written point
	 * \param dt time step
	 * \param R rigidity
	 * \param V_ initial absolute velocity
	 * \param q0_ initial state vector
	 * \param field reference to the magnetic field
	 * \param omega_ gyro-frequency
	 */
	EulerSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, T dt, T R, const T V_, const Vec q0_, FieldGenerator& field, const T omega_);

	/*!
	 * One Euler integration step of the coupled ODE of a particle moving in a background magnetic field
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
