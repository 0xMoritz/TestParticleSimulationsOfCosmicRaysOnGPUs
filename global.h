/*!
 * @file global.h
 * @brief Contains Library includes and defines the StateType
 *
 *  Created on: May 4, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>		// for saving data
#include <cmath>      	// sin, cos, log
#include <limits>  		// for retrieving numeric limits
#include <ctime> 		// for clock and time
#include <cassert> 		// for assertions in debugging
#include <random>		// To generate the random field and random particle start positions
#include <vector>		// As a generic storage medium
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <boost/numeric/odeint.hpp>  // numerical integration
#include <boost/numeric/odeint/external/thrust/thrust.hpp>  // thrust support
#include <thrust/device_vector.h>    // library for interfacing the GPU
#include <thrust/host_vector.h>      // support of the CPU as well


// Value type for all calculations
using T = float;

// Container for particle state (pos and vel) as well as for magnetic field
using Vec = std::vector<T>;

// type for the numerical integration. device=GPU, host=CPU
using StateType = thrust::device_vector<T>;
using StateTypeVec = thrust::device_vector<StateType>;
//using StateType = thrust::host_vector<T>;
//using StateTypeVec = thrust::host_vector<StateType>;

using DeviceVector = thrust::device_vector<T>;
using HostVector = thrust::host_vector<T>;

/*! \mainpage Test particle simulation of cosmic rays on GPUs
 *
 * \section intro_sec Introduction
 *
 * This program's main purpose is the simulation of test particles in isotropically turbulent magnetic fields,
 * by numerically integrating the equations of motion using a Nvidia GPU. Moreover it can be used to perfom said
 * simulation using the CPU, run a test of different numerical integrators and generate isotropic turbulent
 * fields as well as evaluate those fields in a grid for further analysis.
 *
 * \section usage_sec Usage
 *
 * \subsection prerequisites Prerequisites
 * 		The integration of the particle's equations of motion are based on the boost Library ODEint,
 * 		which utilizes thrust. Both must be installed.
 * 		For compilation the NVCC cuda compiler is required. To yield an advantage with the execution on the GPU,
 * 		the amount of particles should be about \f$\gsim 10^3\f$.
 * 		Generating the documentation requires latex, ghostscript and Doxygen.
 * 		For running on the CPU, the comments in the file @p global.h have to be modified.
 *
 * \subsection running Running the program
 * 		The program has three modes of execution:
 * 		"i" for simulating particle trajectories in an isotropic
 * 		turbulent magnetic field. The amount of different field realizations and particles can be specified.
 * 		The particle trajectories are saved in individual files in either logarithmic or linear time steps.
 *
 * 		"B" simulates the motion of a single particle in a background magnetic field without turbulence. The
 *		numerical integration is performed with the Euler method, a Runge Kutta method 2. order, a Runge
 *		Kutta method 4. order using the conventional choice of coefficients and a Runge Kutta method as
 *		implemented by the boost library ODEint.
 *
 * 		"F" can be used to generate a turbulent magnetic field and evaluate it in a grid box. The data is
 * 		saved to a file and can be used for further analysis e.g. using a discrete Fourier transform... *
 *
 */
