/*
 * BoostintSolver.cu
 *
 *  Created on: Jun 16, 2022
 *      Author: moritz
 */

#include "BoostSolver.h"

using namespace std;

// https://www.boost.org/doc/libs/1_64_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/short_example.html

//stateType A = {}
//ptrA = raw_pointer_cast(&A[0])
//for (int i = 0; i < nModes; i++)
//{
//
//}

template<class Tuple>
__host__ __device__
void DerivativeFunctor::operator()(Tuple tuple) // this functor works on tuples of values
{
	// first, unpack the tuple into position values
//	thrust::tuple<T, T, T, T, T, T> q = thrust::get<0>(tuple);
	const T x   = thrust::get<0>(thrust::get<0>(tuple));
	const T y   = thrust::get<1>(thrust::get<0>(tuple));
	const T z   = thrust::get<2>(thrust::get<0>(tuple));
	const T v_x = thrust::get<3>(thrust::get<0>(tuple));
	const T v_y = thrust::get<4>(thrust::get<0>(tuple));
	const T v_z = thrust::get<5>(thrust::get<0>(tuple));
	// Evaluate magnetic field at the position
	T Bx = 0;
	T By = 0;
	T Bz = B0;
	for(int i=0; i<modeCount; i++)
	{
		T kr = ptr_kx[i] * x +  ptr_ky[i] * y + ptr_kz[i] * z; // scalar product of k and pos
		Bx += ptr_Ax[i] * cos(kr + ptr_beta[i]);
		By += ptr_Ay[i] * cos(kr + ptr_beta[i]);
		Bz += ptr_Az[i] * cos(kr + ptr_beta[i]);
	}
	// the differential equation
//	thrust::tuple<T&, T&, T&, T&, T&, T&> dqdt = thrust::get<1>(tuple);
	thrust::get<0>(thrust::get<1>(tuple)) = v_x;							// dx/dt
	thrust::get<1>(thrust::get<1>(tuple)) = v_y;							// dy/dt
	thrust::get<2>(thrust::get<1>(tuple)) = v_z;							// dz/dt
	thrust::get<3>(thrust::get<1>(tuple)) = R_inverse * (v_y*Bz - v_z*By);	// dv_x/dt
	thrust::get<4>(thrust::get<1>(tuple)) = R_inverse * (v_z*Bx - v_x*Bz);	// dv_y/dt
	thrust::get<5>(thrust::get<1>(tuple)) = R_inverse * (v_x*By - v_y*Bx);	// dv_z/dt
};

DerivativeFunctor::DerivativeFunctor(T R_inverse_, T B0_, int modeCount_, const StateType& beta, const StateType& Ax, const StateType& Ay, const StateType& Az, const StateType& kx, const StateType& ky, const StateType& kz)
: R_inverse(R_inverse_), B0(B0_),  modeCount(modeCount_), ptr_beta(thrust::raw_pointer_cast(&beta[0])),
  ptr_Ax(thrust::raw_pointer_cast(&Ax[0])), ptr_Ay(thrust::raw_pointer_cast(&Ay[0])), ptr_Az(thrust::raw_pointer_cast(&Az[0])),
  ptr_kx(thrust::raw_pointer_cast(&kx[0])), ptr_ky(thrust::raw_pointer_cast(&ky[0])), ptr_kz(thrust::raw_pointer_cast(&kz[0]))
{

}

// Derivative Operator //
LorentzForce::LorentzForce(const T R_inverse_, int particleCount_, T B0_, vector<Mode> modes)
: R_inverse(R_inverse_), particleCount(particleCount_), B0(B0_), modeCount(modes.size()), beta(modes.size()), Ax(modes.size()), Ay(modes.size()), Az(modes.size()), kx(modes.size()), ky(modes.size()), kz(modes.size())
{
	// wrap mode parameters to vectors
	Vec betaVec;
	Vec AxVec;
	Vec AyVec;
	Vec AzVec;
	Vec kxVec;
	Vec kyVec;
	Vec kzVec;
	for (int i = 0; i < modes.size(); i++)
	{
		betaVec.push_back(modes[i].beta);
		AxVec.push_back(modes[i].Axi[0]);
		AyVec.push_back(modes[i].Axi[1]);
		AzVec.push_back(modes[i].Axi[2]);
		kxVec.push_back(modes[i].k[0]);
		kyVec.push_back(modes[i].k[1]);
		kzVec.push_back(modes[i].k[2]);
//		cout << "A^2=" <<  to_string(modes[i].Axi[0]*modes[i].Axi[0] + modes[i].Axi[1]*modes[i].Axi[1] + modes[i].Axi[2]*modes[i].Axi[2]) << "k^2=" << to_string(modes[i].k[0]*modes[i].k[0] + modes[i].k[1]*modes[i].k[1] + modes[i].k[2]*modes[i].k[2]) << endl;
	}
	thrust::copy(betaVec.begin(), betaVec.end(), beta.begin());
	thrust::copy(AxVec.begin(), AxVec.end(), Ax.begin());
	thrust::copy(AyVec.begin(), AyVec.end(), Ay.begin());
	thrust::copy(AzVec.begin(), AzVec.end(), Az.begin());
	thrust::copy(kxVec.begin(), kxVec.end(), kx.begin());
	thrust::copy(kyVec.begin(), kyVec.end(), ky.begin());
	thrust::copy(kzVec.begin(), kzVec.end(), kz.begin());
}

//template <class State, class Derivative>
void LorentzForce::operator() (const StateType& q, StateType& dqdt, const T t)
{
	thrust::for_each(
		thrust::make_zip_iterator(thrust::make_tuple(
			// thrust::tuple is limited to a maximum of 10 components.
			// Hence a nested tuple is constructed
			thrust::make_zip_iterator(thrust::make_tuple(
				boost::begin(q) + 0*particleCount, 			// x
				boost::begin(q) + 1*particleCount, 			// y
				boost::begin(q) + 2*particleCount, 			// z
				boost::begin(q) + 3*particleCount, 			// v_x
				boost::begin(q) + 4*particleCount, 			// v_y
				boost::begin(q) + 5*particleCount)), 		// v_z
			thrust::make_zip_iterator(thrust::make_tuple(
				boost::begin(dqdt) + 0*particleCount,		// dx/dt
				boost::begin(dqdt) + 1*particleCount, 		// dy/dt
				boost::begin(dqdt) + 2*particleCount,		// dz/dt
				boost::begin(dqdt) + 3*particleCount,		// dv_x/dt
				boost::begin(dqdt) + 4*particleCount, 		// dv_y/dt
				boost::begin(dqdt) + 5*particleCount)))),	// dv_z/dt
		thrust::make_zip_iterator(thrust::make_tuple(
			thrust::make_zip_iterator(thrust::make_tuple(
				boost::begin(q) + 1*particleCount, 			// x
				boost::begin(q) + 2*particleCount, 			// y
				boost::begin(q) + 3*particleCount, 			// z
				boost::begin(q) + 4*particleCount, 			// v_x
				boost::begin(q) + 5*particleCount, 			// v_y
				boost::begin(q) + 6*particleCount)), 		// v_z
			thrust::make_zip_iterator(thrust::make_tuple(
				boost::begin(dqdt) + 1*particleCount,		// dx/dt
				boost::begin(dqdt) + 2*particleCount, 		// dy/dt
				boost::begin(dqdt) + 3*particleCount,		// dz/dt
				boost::begin(dqdt) + 4*particleCount,		// dv_x/dt
				boost::begin(dqdt) + 5*particleCount, 		// dv_y/dt
				boost::begin(dqdt) + 6*particleCount)))),	// dv_z/dt
		DerivativeFunctor(R_inverse, B0, modeCount, beta, Ax, Ay, Az, kx, ky, kz));
}

// Observer //
Observer::Observer(vector<Vec>& states_, vector<T>& times_)
: states(states_), times(times_)
{ }
void Observer::operator()(const Vec& q, T t)
{
	states.push_back(q);
	times.push_back(t);
}


// Boost Solver //
BoostSolver::BoostSolver(const int particleCount_, const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const FieldGenerator& field, const T omega_, const T Lc_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, field, omega_), particleCount(particleCount_), omega(omega_), Lc(Lc_)
{
	name = "Runge-Kutta (boost ODEint) Simulation";
}
BoostSolver::BoostSolver(const T totalSimulationTime_, const int outputPoints_, const int stepsPerOutput_, const T dt, const T R, const T V_, const Vec q0_, const FieldGenerator& field, const T omega_)
: Engine(totalSimulationTime_, outputPoints_, stepsPerOutput_, dt, R, V_, q0_, field, omega_), particleCount(1), omega(omega_), Lc(0)
{
	name = "Runge-Kutta (boost ODEint) Simulation";
	cout << endl << "Boost: q0z, R, dt, omega=" << q0_[5] << ", " << R << ", " << dt << ", " << omega_ << endl;
}
void BoostSolver::Randomize_q0()
{
	// Random Engine for initial particle velocity direction
	RandomGen* rg = RandomGen::GetInstance();
	// Randomize starting direction
	T phi = rg->RandomFloat_0_2PI();
	T eta = rg->RandomFloat_m1_1();
	// Generate random direction with v=c
	Vec randomVec(UnitVec(eta, phi));
	Scale(randomVec, V); // Apply speed (scaling of position is irrelevant as position is zero)
	q0 = {0, 0, 0}; // [pc],[pc·Ω] Initial position and velocity of the particle
	AppendVector(q0, randomVec); // Append velocity to position in initial state vector q0
	Scale(q0, V);
}
float BoostSolver::Simulation(vector<Vec>& trajectory_, Vec& time_, int batchNo)
{
	// Initialize simulation
	PrintTime();
	const clock_t beginComputingTime = clock(); // https://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
	cout << name << " start" << flush;

	// prepare q0 vector
	Vec q0Vec(6*particleCount);
	//cout << "maxVecSize=" << std::vector<T>::max_size() << endl;
	//cout << "maxdoubleVecSize=" << q0Vec.max_size() << endl;
	//cout << "len(q0Vec)=" << q0Vec.size() << endl;
	for (int n = 0; n < particleCount; n++)
	{
		if (particleCount > 1) // if only one particle is simulated q0 will be set from outside
			Randomize_q0();
		q0Vec[0*particleCount + n] = 0; // x
		q0Vec[1*particleCount + n] = 0; // y
		q0Vec[2*particleCount + n] = 0; // z
		q0Vec[3*particleCount + n] = q0[3]; // v_x
		q0Vec[4*particleCount + n] = q0[4]; // v_y
		q0Vec[5*particleCount + n] = q0[5]; // v_z
	}
	cout << "q=" << q0[0] << ", " << q0[1] << ", " << q0[2] << ", " << q0[3] << ", " << q0[4] << ", " << q0[5] << endl;
	StateType q(6*particleCount);
	//cout << "len(q)=" << q.size() << endl;
	thrust::copy(q0Vec.begin(), q0Vec.end(), q.begin());
	//cout << "len(q)=" << q.size() << endl;
	//cout << "maxStateTypeSize=" << q.max_size() << endl;

	// integrate
	vector<Vec> trajectory;
	Vec time;
	boost::numeric::odeint::runge_kutta4<StateType> stepper;
	LorentzForce lorentzForce(R_inverse, particleCount, field.GetB0(), field.GetModes());
	Observer observer(trajectory, time);

	for (int i = 0; i < outputPoints; i++)
	{
		Vec qVec(6*particleCount);
		//qVec.reserve(6*particleCount);
		thrust::copy(q.begin(), q.end(), qVec.begin());
		observer(qVec, i*dt*stepsPerOutput);
		T time = boost::numeric::odeint::integrate_n_steps(stepper, lorentzForce, q, (T)0., dt, (size_t)stepsPerOutput);//, observer);
		//size_t steps = boost::numeric::odeint::integrate_const(stepper, lorentzForce, q, i*stepsPerOutput*dt, (i+1)*stepsPerOutput*dt, dt);
		//cout << "steps, stepsPerOutput" << steps << ", " << stepsPerOutput << endl;
	}

	// Validate
	assert(trajectory.size() == time.size());
	float timeElapsedInSeconds = float(clock() - beginComputingTime) /  CLOCKS_PER_SEC;
	cout << " finished. Time elapsed: " << timeElapsedInSeconds << " s" << endl;

	// Print
	//vector<Printer*> printers;
	//printers.reserve(particleCount);
	cout << "Writing to files in '" << Printer::GetOutputPath() << "'..." << flush;
	// Instantiate Printers
	for (int particle = 0; particle < particleCount; particle++)
	{
		string filename = "batch" + to_string(batchNo) + "_particle" + to_string(particle) + ".csv";
		string header = "t/(OMEGA^-1); x/Lc; y/Lc; z/Lc; v_x/c, v_y/c, v_z/c";
		T normFac = 1/Lc;
		if (Lc == 0) // This case is for the homogeneous background field, where the correlation length would be infinite
		{
			header = "t/(OMEGA^-1); x/pc; y/pc; z/pc; v_x/c, v_y/c, v_z/c";
			normFac = 1;
		}
		Printer printer(filename, header);
		for (int i = 0; i < outputPoints; i++)
		{
			T t = time[i];
			Vec qVec(trajectory[i]); // Copy construct
			Vec q7Vec;
			q7Vec.push_back(t*omega);
			q7Vec.push_back(qVec[0*particleCount + particle] * normFac);
			q7Vec.push_back(qVec[1*particleCount + particle] * normFac);
			q7Vec.push_back(qVec[2*particleCount + particle] * normFac);
			q7Vec.push_back(qVec[3*particleCount + particle]); // [c]
			q7Vec.push_back(qVec[4*particleCount + particle]); // [c]
			q7Vec.push_back(qVec[5*particleCount + particle]); // [c]
			printer.Write(q7Vec);
		}
	}
	/*for (int particle = 0; particle < particleCount; particle++)
	{
		string filename = "batch" + to_string(batchNo) + "_particle" + to_string(particle) + ".csv";
		string header = "t/(OMEGA^-1); x/Lc; y/Lc; z/Lc; v_x/c, v_y/c, v_z/c";
		printers.push_back(new Printer(filename, header)); // TODO: how to do this without new?
	}
	cout << "len(printers)=" << printers.size() << endl;
	cout << "maxPrinterVecSize=" << printers.max_size() << endl;
	// Print data loop first over points and in the inner loop over particles (and printers) such that all printers Write "in parallel"
	for (int i = 0; i < outputPoints; i++)
	{
		T t = time[i];
		Vec qVec(trajectory[i]); // Copy construct
		for (int particle = 0; particle < particleCount; particle++)
		{
			Vec q7Vec;
			q7Vec.push_back(t*omega);
			q7Vec.push_back(qVec[0*particleCount + particle] / Lc);
			q7Vec.push_back(qVec[1*particleCount + particle] / Lc);
			q7Vec.push_back(qVec[2*particleCount + particle] / Lc);
			q7Vec.push_back(qVec[3*particleCount + particle]); // [c]
			q7Vec.push_back(qVec[4*particleCount + particle]); // [c]
			q7Vec.push_back(qVec[5*particleCount + particle]); // [c]
			printers[particle]->Write(q7Vec);
		}
	}
	// Delete Printers
	for (int particle = 0; particle < particleCount; particle++)
	{
		delete printers[particle];
	}*/
	cout << "finished writing." << endl;
	return timeElapsedInSeconds;
}
void BoostSolver::Step(const T t, Vec& q, Vec& dqdt, Vec& q_out)
{

}

