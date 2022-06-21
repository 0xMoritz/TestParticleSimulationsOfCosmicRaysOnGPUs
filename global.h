/*
 * global.h
 *
 *  Created on: May 4, 2022
 *      Author: moritz
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>		// for saving data
#include <cmath>
#include <limits>  		// for retrieving numeric limits
#include <ctime> 		// for clock and time
#include <cassert> 		// for assertions in debugging
#include <random>		// To generate the random field and random particle start positions
#include <vector>		// As a generic storage medium
#include <sys/stat.h>
#include <stdio.h>


#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

//#include <thrust/reduce.h> //?
//#include <thrust/functional.h> //?

// Value type for all calculations
using T = double;

// Container for particle state (pos and vel) as well as for magnetic field
using Vec = std::vector<T>;
using StateType = thrust::device_vector<T>;
using StateTypeVec = thrust::device_vector<StateType>;
//using StateType = thrust::host_vector<T>;
//using StateTypeVec = thrust::host_vector<StateType>;
using DeviceVector = thrust::device_vector<T>;
using HostVector = thrust::host_vector<T>;

//TODO: compare double, float, long double

