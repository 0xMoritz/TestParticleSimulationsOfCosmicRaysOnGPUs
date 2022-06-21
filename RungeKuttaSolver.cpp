/*
 * RungeKuttaSolver.cpp
 *
 *  Created on: Jun 9, 2022
 *      Author: moritz
 */

#include "RungeKuttaSolver.h"

using namespace std;



RungeKuttaSolver::RungeKuttaSolver(const int particleCount_, const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const FieldGenerator& field, const T omega_, const int order_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, field, omega_)
{
	//TODO: take particleCount into account
	name = "Runge-Kutta Simulation";
	order = order_;
}
RungeKuttaSolver::RungeKuttaSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, const FieldGenerator& field, const T omega_, const int order_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, q0_, field, omega_)
{
	name = "Runge-Kutta Simulation";
	order = order_;
}

void RungeKuttaSolver::Step(const T t, Vec& q, Vec& dqdt, Vec& q_out)
{
	if (order == 2)
		RK2(q, dqdt, q_out);
	else if (order == 4)
		RK4(q, dqdt, q_out);
	//else if (order == 5)
		//RK5(q, dqdt, q_out);
}

void RungeKuttaSolver::RK2(Vec& q, Vec& dqdt, Vec& q_out)
/**
 * Second order Runge-Kutta > Numerical Recipes Eq. 17.1.2
 *
 * k_1 = dt*f(t_n, y_n)								(1)
 * k_2 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)			(2)
 * y_{n+1} = y_n + k_2								(3)
 */
{
	int n = q.size();
	Vec qt(n), dqt(n);
	for(int i=0;i<n;i++)
	{
		qt[i] = q[i] + 0.5*dt*dqdt[i]; // 			(1)
	}
	Derivatives(qt, dqt, field.BField(q));
	for(int i=0;i<n;i++)
	{
		q_out[i] = q[i] + dt*dqt[i]; // 			(2) & (3)
	}
	dqdt=dqt; // Derivative for next iteration
}

void RungeKuttaSolver::RK4(Vec& q, Vec& dqdt, Vec& q_out)
/**
 * Fourth order Runge-Kutta > Numerical Recipes Eq. 17.1.3
 *
 * k_1 = dt*f(t_n, y_n)								(1)
 * k_2 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)			(2)
 * k_3 = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)			(3)
 * k_4 = dt*f(t_n + dt, y_n + k_3)					(4)
 * y_{n+1} = y_n + k_1/6 + k_2/3 + k_3/3 + k_4/6	(5)
 */
{
	size_t n = q.size(); // 									(1)
	Vec dqm(n), dqt(n), qt(n); // temporary variables
	T dt_2 = 0.5*dt;
	T dt_6 = dt/6.0;
	for(size_t i=0;i<n;i++)
	{
		qt[i] = q[i] + dt_2*dqdt[i];
	}
	Derivatives(qt, dqt, field.BField(q)); // 					(2)
	for(size_t i=0;i<n;i++)
	{
		qt[i] = q[i] + dt_2*dqt[i];
	}
	Derivatives(qt, dqm, field.BField(q)); // 					(3)
	for(size_t i=0;i<n;i++)
	{
		qt[i]=q[i]+dt*dqm[i];
		dqm[i]+=dqt[i];
	}
	Derivatives(qt, dqt, field.BField(q)); //					(4)
	for(size_t i=0;i<n;i++)
	{
		q_out[i] = q[i] + dt_6*(dqdt[i]+dqt[i]+2.*dqm[i]); // 	(5)
	}
	dqdt=dqt; // Derivative for the next iteration
}


