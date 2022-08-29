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
