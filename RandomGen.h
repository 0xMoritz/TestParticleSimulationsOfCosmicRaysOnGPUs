/*!
 * @file RandomGen.h
 * @brief Auxillary class for random number generation
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "global.h"


class RandomGen
/*!
 * Singleton class for generating random numbers
 */
{
// https://gist.github.com/pazdera/1098119
private:
	// std::random_device rd;  // Could be used to obtain random number for seed
	// https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	static RandomGen* instance; 									//!< single instance of this class
	std::mt19937 RdmEngine; 											//!< Mersenne twister engine
	std::uniform_real_distribution<T> dist_0_2PI; //!< uniform distribution between 0 and 2π
	std::uniform_real_distribution<T> dist_m1_1;  //!< uniform distribution between -1 and +1

	/*
	 * Private Constructor to create @ref instance.
	 */
	RandomGen();

public:
	/*!
	 * An new instance will be created with the first call of this function.
	 * \return the single @ref instance of this class
	 */
	static RandomGen* GetInstance();

	/*!
	 * Change the seed of the random generator
	 * \param seed set new seed for the random number engine
	 */
	void ChangeSeed(const int seed);

	/*!
	 * Returns a random floating point value, uniformally distributed between \f$0\f$ and \f$2 \pi\f$
	 */
	T RandomFloat_0_2PI();

	/*!
	 * Returns a random floating point value, uniformally distributed between \f$ -1\f$ and \f$ +1\f$
	 */
	T RandomFloat_m1_1();
};
