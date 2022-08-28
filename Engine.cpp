/*
 * Engine.cpp
 *
 *  Created on: May 6, 2022
 *      Author: Moritz Geßner
 */

#include "Engine.h"

using namespace std;

/* Derivatives Function for a constant field in the z direction
void Engine::Derivatives(const Vec& q, Vec& dqdt, const T& OMEGA)
**
 * Calculates the Derivatives @dqdt dependent on the velocity
 * the @q vector has Components q=(x, y, z, v_x, v_y, v_z)^T
 * Here is the actual differential equation:
 *
 * dv_x/dt =  OMEGA*v_y			(1)
 * dv_y/dt = -OMEGA*v_x			(2)
 * dx/dt = v_x					(3)
 * dy/dt = v_y					(4)
 * dz/dt = v_z					(5)
 *
{
	dqdt[0] = q[3]; //		 	(3)
	dqdt[1] = q[4]; //		 	(4)
	dqdt[2] = q[5]; //		 	(5)
	dqdt[3] =  OMEGA*q[4]; //	(1)
	dqdt[4] = -OMEGA*q[3]; //	(2)
	dqdt[5] = 0;
}*/

void Engine::Derivatives(const Vec& q, Vec& dqdt, const Vec B)
{
	// (See declaration for differential equation)
	dqdt[0] = q[3]; //								(1)
	dqdt[1] = q[4]; //								(2)
	dqdt[2] = q[5]; //  							(3)
	dqdt[3] = R_inverse*(q[4]*B[2] - q[5]*B[1]); //	(4)
	dqdt[4] = R_inverse*(q[5]*B[0] - q[3]*B[2]); //	(5)
	dqdt[5] = R_inverse*(q[3]*B[1] - q[4]*B[0]); //	(6)
}

Engine::Engine(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt_, const T R, const T V_, const FieldGenerator& field_, const T omega_)
: totalSimulationTime(totalSimulationTime_), outputPoints(outputPoints_), stepsPerOutput(stepsPerOutput_), dt(dt_), V(V_), field(field_), R_inverse(1/R)
{
	firstWriteTime = 0.1*omega_;
	tWrite = firstWriteTime;
	timePerOutputIncrease = pow(totalSimulationTime/firstWriteTime, 1./outputPoints);	// [pc·c⁻¹]
}

Engine::Engine(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt_, const T R, const T V_, const Vec q0_, const FieldGenerator& field_, const T omega_)
: totalSimulationTime(totalSimulationTime_), outputPoints(outputPoints_), stepsPerOutput(stepsPerOutput_), dt(dt_), V(V_), q0(q0_), field(field_), R_inverse(1/R)
{
	firstWriteTime = 0.1*omega_;
	tWrite = firstWriteTime;
	timePerOutputIncrease = pow(totalSimulationTime/firstWriteTime, 1./outputPoints);	// [pc·c⁻¹]
}
float Engine::Simulation(vector<Vec>& trajectory, Vec& time, int batchNo)
{
	Vec q(q0); // running state variable TODO:use copy-constructor
	Vec dqdt(6);
	Derivatives(q, dqdt, field.BField(q)); // Stepper functions needs derivative to start
	Vec q_out(6);

	// Initiate Simulation
	PrintTime();
	const clock_t beginComputingTime = clock(); // https://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
	cout << name << " start" << flush;
	T t = 0; 	// [pc·c⁻¹]
	// Linear:
	trajectory.push_back(q); // Write zeroth entry
	time.push_back(t);
	for (int i = 0; i < outputPoints-1; i++)
	{
		if (i % (outputPoints/10) == 0)
			cout << "." << flush;
		for (int j = 0; j < stepsPerOutput; j++) // Do more simulation steps than are added to the output file
		{
			// Simulation step
			this->Step(t, q, dqdt, q_out);
			t += dt;
			q = q_out;
		}
		//Write point to storage
		trajectory.push_back(q);
		time.push_back(t);
	}
	// Logarithmic:
//	while (t <= totalSimulationTime)
//	{
//		// Simulation step
//		this->Step(t, q, dqdt, q_out);
//		t += dt;
//		q = q_out;
//		while (t > tWrite && t <= totalSimulationTime) // "While" in case one simulation step is smaller then the writing step
//		{
//			//cout << "tWrite:" << tWrite;
//			numWrites++;
//			if (t > totalSimulationTime*(numWrites/10))
//			{
//				cout << "." << flush;
//				numWrites++;
//			}
//			tWrite *= timePerOutputIncrease;
//			//Write point to storage
//			trajectory.push_back(q);
//			time.push_back(t);
//		}
//	}
	float timeElapsedInSeconds = float(clock() - beginComputingTime) /  CLOCKS_PER_SEC;
	cout << " finished. Time elapsed: " << timeElapsedInSeconds << " s" << endl;
	return timeElapsedInSeconds;
}
