/*
 * RandomGen.h
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
	// For Randomization:
	//std::random_device rd;  // Will be used to obtain a seed for the random number engine
	// https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	static RandomGen* instance;
	std::mt19937 RdmEngine;
	std::uniform_real_distribution<T> dist_0_2PI;
	std::uniform_real_distribution<T> dist_m1_1;

	/*
	 * Private Constructor to create @instance.
	 */
	RandomGen();

public:

	/*!
	 * An new instance will be created with the first call of this function.
	 */
	static RandomGen* GetInstance();

	/*!
	 * Change the seed of the random generator
	 */
	void ChangeSeed(const int seed);

	/*!
	 * Returns a random floating point value between 0 and 2 \f$ \pi\f$
	 */
	T RandomFloat_0_2PI();

	/*!
	 * Returns a random floating point value between -1 and 1
	 */
	T RandomFloat_m1_1();
};
