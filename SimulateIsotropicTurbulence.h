/*
 * IsotropicTurbulenceSimulation.h
 *
 *  Created on: May 18, 2022
 *      Author: moritz
 */

#pragma once

#include "RungeKuttaSolver.h"
#include "BoostSolver.h"

/*!
 * Simulate the propagation of a single particle through a general magnetic @field
 */

/*!
 * Generates a isotropic turbulent field and simulates a number of particles with random
 * initial directions to diffuse through the field
 * The simulation properties are specified by the dimensionless parameters @R_norm, @eta and @modeCount:
 *	@R_norm is the normalized Rigidity, i.e. rigidity divided by coherence length
 * \f$ \tilde{R} = R_{norm} = \frac{R}{L_C \sqrt{B_0^2+\delta B^2}} \f$ ;
 * Rigidity: \f$ R = B ρ = \frac{p}{Q} \f$ , where ρ is the Larmor radius, p momentum, q charge
 * 	@eta is a dimensionless constant that determines the share of turbulence in the magnetic field:
 * \f$ \eta = \frac{ \langle \delta\vec{B}^2 \rangle }
 * { \langle \delta\vec{B}^2 \rangle + \vec{B}_0^2 } \f$
 * 	@modeCount is the number of modes, that are generated and summed over for calculating
 * the magnetic field. The modes will scale logarithmically from 0.04 pc⁻¹ to 40 pc⁻¹
 * 	@particleCount is the amount of particles simulated,
 * 	the trajectory for each particle will be stored in a separate file.
 * 	@SimulationTimeInGyrationCycles: time simulated in the simulation given in gyration periods
 * 	@MINSIMSTEPS: minimum number of integration steps per particle
 * 	@OUTPUTPOINTS: number of output lines generated
 */
void IsotropicTurbulenceSimulation(const T R_norm, const T eta, const T gamma, const T Lmin, const T Lmax, const int modeCount, const int fieldCount, const int particlesPerFieldCount, const T SimulationTimeInGyrationCycles, const int minSimSteps, const int outputPoints, const bool useBoost, const int seed);

