/*!
 * @file SimulateBackgroundField.h
 * @brief Prepares simulatinga particle moving through a background field using all solvers
 *
 *  Created on: May 18, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "BackgroundAnalyticSolution.h"
#include "EulerSolver.h"
#include "RungeKuttaSolver.h"
#include "BoostSolver.h"

/*!
 * Initializing the numerical simulation for a particle with rigidity R in a
 * homogeneous magnetic field parallel to the z-direction
 * 	@SimulationTimeInGyrationCycles: time simulated in the simulation given in gyration periods
 * 	@MINSIMSTEPS: minimum number of integration steps per particle
 * 	@OUTPUTPOINTS: number of output lines generated
 */
void SimulateBackgroundField(const T R_norm, const T SimulationTimeInGyrationCycles, const int minSimSteps, const int outputPoints);
