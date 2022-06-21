/*
 * BoostintSolver.h
 *
 *  Created on: Jun 9, 2022
 *      Author: moritz
 */

#pragma once

#include "Engine.h"

struct DerivativeFunctor
{
public:
    const T R_inverse;
    const T B0;
    const int modeCount;
    const T* ptr_beta;
    const T* ptr_Ax;
    const T* ptr_Ay;
    const T* ptr_Az;
    const T* ptr_kx;
    const T* ptr_ky;
    const T* ptr_kz;

    DerivativeFunctor(T R_inverse_, T B0_, int modeCount_, const StateType& beta, const StateType& Ax, const StateType& Ay, const StateType& Az, const StateType& kx, const StateType& ky, const StateType& kz);

	template<class Tuple>
	__host__ __device__
	void operator()(Tuple tuple);  // this functor works on tuples of values
};

/*!
 *  Object to bring the Derivative function of Engine to the desired form for boost
 */
class LorentzForce
{
private:
	//FieldGenerator& field;
    const T R_inverse; // 1/R = Q/m/γ, [µG⁻¹·pc⁻¹]
    const int particleCount;
    const T B0;
    const int modeCount;
    StateType beta;
    StateType Ax;
    StateType Ay;
    StateType Az;
    StateType kx;
    StateType ky;
    StateType kz;

public:
    LorentzForce(const T R_inverse_, int particleCount_, T B0_, std::vector<Mode> modes_);
    void operator() (const StateType& q, StateType& dqdt, const T t);// const;
};

/*!
 * Observer class for writing points to a file (logarithmically with time) during
 * the integration process.
 */
struct Observer
{
public:
    std::vector<Vec>& states;
    Vec& times;
    Observer(std::vector<Vec>& states_, Vec& times_);

	/*!
	 * writes @x and @t to @m_states and @m_times after
	 * @timePerOutputPoint has passed.
	 */
    void operator()(const Vec& q, T t);
};

/*!
 * Engine to solve the differential equation by making use of the Boost ODEint library.
 * The differential equation an are given with the @DerivativeOperator class and
 * the data is recorded with the @Observer class. The Engine uses a constant
 * Runge-Kutta method of 4th order in @dt and runs on the CPU.
 */
class BoostSolver : public Engine
{
private:
	const int particleCount;
	void Randomize_q0();
	const T omega;
	const T Lc;

public:
	BoostSolver(const int particleCount_, const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V, const FieldGenerator& field, const T omega_, const T Lc_);
	BoostSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V, const Vec q0_, const FieldGenerator& field, const T omega_);

	/*!
	 * Override of the Simulation function, as the integration is triggered by a single
	 * call instead of looping the stepper manually.
	 */
	float Simulation(std::vector<Vec>& trajectory, Vec& time, int batchNo) override;

	/*!
	 * Dummy function
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
