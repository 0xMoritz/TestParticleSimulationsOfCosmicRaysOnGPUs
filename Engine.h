/*!
 * @file Engine.h
 * @brief Base class for ODE solvers
 *
 *  Created on: May 6, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "utility.h"
#include "Printer.h"
#include "FieldGenerator.h"

class Engine
{
protected:
	const T totalSimulationTime; /*! total simulation time */
	const int outputPoints;
	const int stepsPerOutput;
	const T dt; 				// [pc·c⁻¹]
	const T R_inverse; 		// [µG⁻¹·pc⁻¹] Charge divided by mass and Lorentz factor Q_MGamma = Q/m/γ
	const T V; 				// [c] abs velocity
	Vec q0 = {0,0,0,0,0,0};
	std::string name = ">>>MISSING<<<";
	const FieldGenerator& field;
	T firstWriteTime;
	T tWrite;
	int numWrites = 0;
    T timePerOutputIncrease;

	/*! @brief Calculates the derivatives at given state \f$q\f$.
	 *
	 * Calculates the Derivatives @dqdt dependent on the velocity for the differential equation
	 * the @q vector has Components \f$  q=(x, y, z, v_x, v_y, v_z)^T \f$ , @B is the magnetic field,
	 * The equations of motions are determined by the Lorentz-force
	 * \f$ \frac{d^2 r}{dt^2} = \frac{q}{m\gamma} v \times B = \frac{1}{R} v \times B \f$
	 * Where \f$\gamma\f$ is assumed to be constant (this is weakly violated by the Runge Kutta integration).
	 *
	 * \f$ dv_x/dt = 1/R(v2B3-v3B2) \f$ 			(1) \n
	 * \f$ dv_y/dt = 1/R(v3B1-v1B3) \f$ 			(2) \n
	 * \f$ dv_z/dt = 1/R(v1B2-v2B1) \f$ 			(3) \n
	 * \f$ dx/dt = v_x \f$ 							(4) \n
	 * \f$ dy/dt = v_y \f$ 							(5) \n
	 * \f$ dz/dt = v_z \f$ 							(6) \n
	 */
	void Derivatives(const Vec& q, Vec& dqdt, const Vec B);

public:
	/*!
	 * All properties of the particle are specified in this constructor
	 */
	Engine(const T totalSimulationTime_, int outputPoints_, int stepsPerOutput_, T dt_, T R, T V_, const FieldGenerator& field_, const T omega_);
	/// Brief description.
	/*!
	 * Detailed description.
	 */
	Engine(const T totalSimulationTime_, int outputPoints_, int stepsPerOutput_, T dt_, T R, T V_, Vec q0_, const FieldGenerator& field_, const T omega_);

	/*!
	 * calculation of one integration step. t is the absolute time, @q is the starting vector,
	 * @dqdt is the derivative of @q, that is expected as a valid evaluation at the point @q
	 * as an input, the evaluation will yield an output point @q_out. @dqdt will be updated and
	 * return the derivative at the point @q_out.
	 */
	virtual void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) = 0;

	/*!
	 * Will perform the simulation of the particle with the properties specified in the constructor.
	 * And return a \f$ q=(x,y,z, v_x, v_y, v_z) \f$  states vector @trajectory and time vector @time
	 * the return parameter is the computing time for the simulation in seconds
	 */
	virtual float Simulation(std::vector<Vec>& trajectory, Vec& time, int batchNo);
};
