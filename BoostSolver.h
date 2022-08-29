/*!
 * @file BoostSolver.h
 * @brief Solver using boost ODEint's Runge Kutta implementation
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "Engine.h"
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>

/*! \brief Functor for calculating the right hand side of the equations of motion on GPU
 *
 * Executing code on the GPU with thrust works via functors. These must contain a member function @p
 * operator() that is called in a @p thrust::for_each loop. Formatting the derivative in this way allows
 * thrust to perform the operation on the GPU, the order of execution is not defined. But since there is no
 * interaction between the particles this is not our concern. For calculating the magnetic field we can
 * unfortunately not use the @ref FieldGenerator class directly, properties of the Derivative Functor must be
 * stored as pointer, they are passed as @ref Statevector, i.e. as device_vector or host_vector depending on
 * the used processing unit.
 */
struct DerivativeFunctor
{
public:
    const T R_inverse;      //!< inverse rigidity
    const T B0;             //!< background component of the magnetic field
    const int modeCount;    //!< number of modes \f$n\f$
    const T* ptr_beta;      //!< pointer to the first element of an array containing the phases \f$\beta_i\f$
    const T* ptr_Ax;        //!< pointer to the first element of an array containing the Amplitude times the polarization in x direction
    const T* ptr_Ay;        //!< pointer to the first element of an array containing the Amplitude times the polarization in y direction
    const T* ptr_Az;        //!< pointer to the first element of an array containing the Amplitude times the polarization in z direction
    const T* ptr_kx;        //!< pointer to the first element of an array containing the x component of the wavevector
    const T* ptr_ky;        //!< pointer to the first element of an array containing the y component of the wavevector
    const T* ptr_kz;        //!< pointer to the first element of an array containing the z component of the wavevector

    /*! \brief constructor takes field properties and rigidity
     *
     * The constructor takes all field information (background field @p B0_) and modes as reference of device- or host-vector
     */
    DerivativeFunctor(T R_inverse_, T B0_, int modeCount_, const StateType& beta, const StateType& Ax, const StateType& Ay, const StateType& Az, const StateType& kx, const StateType& ky, const StateType& kz);

    /*! \brief calculates the derivative
     *
     * The operator satisfies the format required to be used in a @p thrust::for_each loop with
     * thrust::zip_iterator.
     * \param tuple Contains all particle states in a nested format
     */
	template<class Tuple>
	__host__ __device__
	void operator()(Tuple tuple);  // this functor works on tuples of values
};

/*!
 *  @brief integrator object for boost's integrate functions
 *
 * This object satisfies the requirements for the boost's integrator functions. It incorporates a thrust::tuple to integrate all particles in unison using the @ref DerivativeFunctor
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
    /*! \brief constructor supplies @ref DerivativeFunctor with @ref mode information
     *
     * This constructor initializes the @ref DerivativeFunctor and provides it with the pointers to mode properties
     */
    LorentzForce(const T R_inverse_, int particleCount_, T B0_, std::vector<Mode> modes_);
    /*! \brief operator for a single integration step of all particles
     *
     * Calling this function calculates the derivatives of all particles for propagation of the equations of
     * motion, in a boost integrator routine @q contains all particle states and @dqdt its derivatives in a
     * nested format
     */
    void operator() (const StateType& q, StateType& dqdt, const T t);// const;
};

/*! \brief Object to extract time and position of the integration process
 *
 * Observer class for writing points to a file (logarithmically with time) during
 * the integration process.
 */
struct Observer
{
public:
    std::vector<Vec>& states; //!< reference to a list of states to write the observed trajectory into
    Vec& times;				  //!< referente to a list of the time for each corresponding state
    /*! \brief Constructor that connects the internal @p states and @p times references to an external variable
     *
     * Given the vectors @p states_ and @p times_, @p Observer will create references to write the observed
     * trajectories into these containers
     */
    Observer(std::vector<Vec>& states_, Vec& times_);

	/*!
	 * writes @x and @t to @m_states and @m_times after
	 * @timePerOutputPoint has passed.
	 */
    void operator()(const Vec& q, T t);
};

/*! \brief  Engine to solve the differential equation by making use of the Boost ODEint library.
 *
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
	const bool logTime;

public:
    /*! @brief Constructor for simulation in isotropic turbulence
     *
     * Prepare the solver for the integration of test particles in a isotropically turbulent magnetic field
     * @param V the initial velocity of the particle
     * @param logTime_ determines, whether the times at which the trajectory is observed and written to files
     * is logarithmic or linear
     */
	BoostSolver(const int particleCount_, const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V, const FieldGenerator& field, const T omega_, const T Lc_, const bool logTime_);
    /*! @brief Constructor for test in background field
     *
     * Only one particle will be simulated in the background field.
     * @param field is expected to consist only of a background field
     * @param q0_ is the initial state of the particle
     */
	BoostSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V, const Vec q0_, const FieldGenerator& field, const T omega_);

	/*!
	 * Override of the Simulation function, as the integration is triggered by a single
	 * call instead of looping the stepper manually.
	 */
	float Simulation(std::vector<Vec>& trajectory, Vec& time, int batchNo) override;

	/*!
	 * Dummy function. Because @ref Simulation is overritten this will not be called for @p BoostSolver
	 */
	void Step(const T t, Vec& q, Vec& dqdt, Vec& q_out) override;
};
