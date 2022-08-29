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
	const T totalSimulationTime; 			//!< total time the particle should be simulated
	const int outputPoints;					//!< amount of positions recorded in the output file
	const int stepsPerOutput;				//!< integration steps per output point of the file
	const T dt; 							//!< \f$\left[\mathrm{pc}\cdot c^{-1}\right]\f$ time step
	const T R_inverse; 						//!< \f$\left[\mu\mathrm{G}^{-1}\cdot \mathrm{pc}^{-1}\right]\f$ Inverse rigidity \f$ \frac{1}{R}\f$
	const T V; 								//!< \f$\left[c\right]\f$ absolute velocity
	Vec q0 = {0,0,0,0,0,0};					//!< Initial state vector
	std::string name = ">>>MISSING<<<";		//!< Name of the Numerical solver used
	const FieldGenerator& field;			//!< Reference to the magnetic field
	// members for logarithic output times
	T firstWriteTime;						//!< time after which the first output point is written
	T tWrite;								//!< time after which the next output point can be written
	int numWrites = 0;						//!< recorded output points so far
	T timePerOutputIncrease;				//!< time that @p tWrite increases after a point has been recorded

	/*! @brief Calculates the derivatives at given state \f$q\f$.
	 *
	 * Calculates the Derivatives @dqdt dependent on the velocity for the differential equation
	 * the @q vector has Components \f$  q=(x, y, z, v_x, v_y, v_z)^T \f$ , @p B is the magnetic field,
	 * The equations of motions are determined by the Lorentz-force
	 * \f$ \frac{d^2 r}{dt^2} = \frac{q}{m\gamma} v \times B = \frac{1}{R} v \times B \f$
	 * Where \f$\gamma\f$ is assumed to be constant. The changing velocity due to integration errors will only be * considered classically.
	 *
	 * \f$ \frac{dv_x}{dt} = \frac{1}{R}(v_y B_z - v_z B_y) \f$ 			(1) \n
	 * \f$ \frac{dv_y}{dt} = \frac{1}{R}(v_z B_x - v_x B_z) \f$ 			(2) \n
	 * \f$ \frac{dv_z}{dt} = \frac{1}{R}(v_x B_y - v_y B_x) \f$ 			(3) \n
	 * \f$ \frac{dx}{dt} = v_x \f$ 							(4) \n
	 * \f$ \frac{dy}{dt} = v_y \f$ 							(5) \n
	 * \f$ \frac{dz}{dt} = v_z \f$ 							(6) \n
	 */
	void Derivatives(const Vec& q, Vec& dqdt, const Vec B);

public:
	/*!
	 * All properties of the particle are specified in this constructor
	 */
	Engine(const T totalSimulationTime_, int outputPoints_, int stepsPerOutput_, T dt_, T R, T V_, const FieldGenerator& field_, const T omega_);
	/*!
	 * Constructor with the additional specification of the starting state @p q0.
	 */
	Engine(const T totalSimulationTime_, int outputPoints_, int stepsPerOutput_, T dt_, T R, T V_, Vec q0_, const FieldGenerator& field_, const T omega_);

	/*!
	 * calculation of one integration step.
	 * \param t absolute time,
	 * \param q the starting vector,
	 * \param dqdt is the derivative of @p q, this is expected to be a valid evaluation at the point @p q
	 * as an input. @p dqdt will be updated to the derivative at the point @p q_out.
	 * \param the evaluation will yield an output point @p q_out.
	 */
	virtual void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) = 0;

	/*!
	 * Performs the simulation of the particle with the properties specified in the constructor.
	 * It will return a \f$ q=(x,y,z, v_x, v_y, v_z) \f$  states vector @p trajectory and time vector @p time.
	 * \return the computing time for the simulation in seconds
	 */
	virtual float Simulation(std::vector<Vec>& trajectory, Vec& time, int batchNo);
};
