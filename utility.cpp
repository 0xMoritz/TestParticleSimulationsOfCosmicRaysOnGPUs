/*
 * utility.cpp
 *
 *  Created on: May 4, 2022
 *      Author: moritz
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
//template <typename Type>
//void AppendVector(vector<Type>& a, const vector<Type>& b)
void AppendVector(Vec& a, const Vec& b)
{
	a.insert(a.end(), b.begin(), b.end()); // https://stackoverflow.com/questions/2551775/appending-a-vector-to-a-vector
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
Vec UnitVec(const T eta, const T phi)
{
	Vec q(3);
	q[0] = sqrt(1-eta*eta) * cos(phi);
	q[1] = sqrt(1-eta*eta) * sin(phi);
	q[2] = eta;
	return q;
}


