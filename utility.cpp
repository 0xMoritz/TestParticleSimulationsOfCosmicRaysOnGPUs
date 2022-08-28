/*
 * utility.cpp
 *
 *  Created on: May 4, 2022
 *      Author: Moritz Geßner
 */

#include "utility.h"

using namespace std;

// Vector Operations
T SqNorm(Vec v)
{
	T s = 0;
	for (size_t i=0;i<v.size();i++)
	{
		s += v[i]*v[i];
	}
	return s;
}


T ScalarProd(const Vec& v, const Vec& w)
{
	assert(v.size() <= w.size());
	T s = 0;
	for (size_t i=0;i<v.size();i++)
	{
		s += v[i]*w[i];
	}
	return s;
}


void Scale(Vec& v, T a)
{
	for (size_t i=0;i<v.size();i++)
	{
		v[i] = v[i]*a;
	}
}


Vec Sum(const Vec& v, const Vec& w)
{
	assert(v.size() == w.size());
	Vec u(v.size());
	for(size_t i=0;i<v.size();i++)
	{
		u[i] = v[i] + w[i];
	}
	return u;
}


void VecAdd(Vec& v, const Vec& w)
{
	assert(v.size() == w.size());
	Vec u(v.size());
	for(size_t i=0;i<v.size();i++)
	{
		v[i] = v[i] + w[i];
	}
}


void AppendVector(Vec& a, const Vec& b)
{
	// https://stackoverflow.com/questions/2551775/appending-a-vector-to-a-vector
	a.insert(a.end(), b.begin(), b.end());
}


T VSquared(const Vec& q)
{
	return (q[3]*q[3] + q[4]*q[4] + q[5]*q[5]); // classical kinetic energy
}


// Time
void PrintTime()
{
	time_t currentTime;
	struct tm *localTime;
	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	// https://stackoverflow.com/questions/530614/print-leading-zeros-with-c-output-operator
	// https://stackoverflow.com/questions/43493794/how-to-get-local-time-c
	cout << "[" << setw(2) << setfill('0') << localTime->tm_hour << ":"
				<< setw(2) << setfill('0') << localTime->tm_min << ":"
				<< setw(2) << setfill('0') << localTime->tm_sec << "] " << flush;
}


// Unitvec
Vec UnitVec(const T zeta, const T phi)
{
	Vec q(3);
	q[0] = sqrt(1-zeta*zeta) * cos(phi);
	q[1] = sqrt(1-zeta*zeta) * sin(phi);
	q[2] = zeta;
	return q;
}
