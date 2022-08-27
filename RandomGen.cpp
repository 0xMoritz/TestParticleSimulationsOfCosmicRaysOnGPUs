/*
 * RandomGen.cpp
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#include "RandomGen.h"

using namespace std;


// class RandomGen

/* Null, because instance will be initialized on demand. */
RandomGen* RandomGen::instance = 0;

RandomGen* RandomGen::GetInstance()
{
    if (instance == 0)
    {
        instance = new RandomGen();
    }

    return instance;
}
RandomGen::RandomGen() : RdmEngine(0), dist_0_2PI(0., 2.0*M_PI), dist_m1_1(-1., 1.)
{

}
void RandomGen::ChangeSeed(const int seed)
{
	RdmEngine.seed(seed);
}
T RandomGen::RandomFloat_0_2PI()
{
	return dist_0_2PI(RdmEngine);
};
T RandomGen::RandomFloat_m1_1()
{
	return dist_m1_1(RdmEngine);
};
