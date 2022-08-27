/*
 * SimulateBackgroundField.cpp
 *
 *  Created on: May 18, 2022
 *      Author: Moritz Geßner
 */

#include "SimulateBackgroundField.h"

using namespace std;


void SimulateBackgroundField(const T R_norm, const T simulationTimeInGyroperiods, const int minSimSteps, const int outputPoints)
{
	// Parameter list (units pc,µG,c):
	const T B0 = 4; 					// [µG] strength of background magnetic field
	const T V = 1;						// [c] absolute velocity of particle
	// Magnetic field
	const T B_mean = B0;				// [µG] mean B field
	// Particle related parameters
	const T R = R_norm*B_mean*1;		// [µG·pc] rigidity of particle ??
	//const T M_Q = R/V; 					// [µG·pc·c⁻¹] mass divided by charge
	const T Q_M = V/R;					// [c·pc⁻¹·µG⁻¹] charge divided by mass
	const T omega = Q_M*B_mean;			// [c·pc⁻¹] gyration frequency Ω for the mean Field
//	cout << "OMEGA=" << OMEGA << endl;
	// Simulation related parameters
	const T totalSimulationTime = simulationTimeInGyroperiods*2*M_PI/omega;	// [pc·c⁻¹] time simulated in the simulation
	assert(minSimSteps >= outputPoints);

	// Additional parameters for the background field
	const Vec q0= {0, 0, 0, V/sqrt(2.), 0, V/sqrt(2.)}; // Initial position and velocity of the particle

	// Composite parameters
	const int stepsPerOutput = minSimSteps % outputPoints==0 ? minSimSteps/outputPoints : minSimSteps/outputPoints+1;
	const int actualSimSteps = stepsPerOutput*outputPoints;
	const T dt = totalSimulationTime/actualSimSteps; // [pc·c⁻¹]
	const T VV0 = VSquared(q0); // [c²]

	FieldGenerator field(B0); //TODO: Simulationtime in Gyroperiods doesn't work

	// Print parameters
	cout << "Parameters:"
			<< "\n\tR_norm = " << R_norm
			<< "\n\tSimulationTimeInGyroperiods = " << simulationTimeInGyroperiods << " OMEGA^-1"
				<< " = " << totalSimulationTime*omega << "pc/c"
			<< "\n\tdt = " << dt*omega <<  " OMEGA^-1"
				<< " = " << dt << " pc/c"
			<< "\n\tminSimSteps = " << minSimSteps
			<< "\n\toutputPoints = " << outputPoints
			<< "\n\tactualSimSteps = " << actualSimSteps << endl;

	Engine* engine;


	// Doing the calculation with the analytical result
	BackgroundAnalyticSolution theo(totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, q0, field, omega);
	vector<Vec> theoTrajectory; // store the positions
	theoTrajectory.reserve(outputPoints);
	Vec theoTime; // store the time
	theoTime.reserve(outputPoints);
	engine = &theo;
	float theoCompTime = engine->Simulation(theoTrajectory, theoTime, 0);

	// Simulation with Euler method
	EulerSolver eulr(totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, q0, field, omega);
	vector<Vec> eulrTrajectory;
	eulrTrajectory.reserve(outputPoints);
	Vec eulrTime;
	eulrTime.reserve(outputPoints);
	engine = &eulr;
	float eulrCompTime = engine->Simulation(eulrTrajectory, eulrTime, 0);

	// Simulation with Runge Kutta method second order
	RungeKuttaSolver ruKu2(totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, q0, field, omega, 2);
	vector<Vec> ruKu2Trajectory;
	ruKu2Trajectory.reserve(outputPoints);
	Vec ruKu2Time;
	ruKu2Time.reserve(outputPoints);
	engine = &ruKu2;
	float ruKu2CompTime = engine->Simulation(ruKu2Trajectory, ruKu2Time, 0);

	// Simulation with Runge Kutta method fourth order
	RungeKuttaSolver ruKu4(totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, q0, field, omega, 4);
	vector<Vec> ruKu4Trajectory;
	ruKu4Trajectory.reserve(outputPoints);
	Vec ruKu4Time;
	ruKu4Time.reserve(outputPoints);
	engine = &ruKu4;
	float ruKu4CompTime = engine->Simulation(ruKu4Trajectory, ruKu4Time, 0);

	// Simulation with Boost ODEint Runge Kutta
	BoostSolver boost(totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, q0, field, omega);
	vector<Vec> boostTrajectory;
	boostTrajectory.reserve(outputPoints);
	Vec boostTime;
	boostTime.reserve(outputPoints);
	engine = &boost;
	float boostCompTime = engine->Simulation(boostTrajectory, boostTime, 0);
	Scale(boostTime, omega); // correct boost Time to be [boostTime] = [t*omega] = [1]

	// Writing to file
	PrintTime();
	string filename = "homogeneousBackground.csv";
	cout << "Printing to file '" << filename << "'..." << flush;
	Printer printer(filename, "t*OMEGA; theory: x/pc; y/pc; z/pc; vx/c; vy/c; vz/c; V^2;"
			"Euler: x; y; z; vx; vy; vz; V^2/c^2;"
			"Runge-Kutta2: x/pc; y/pc; z/pc; vx/c; vy/c; vz/c; V^2/c^2;"
			"Runge-Kutta4: x/pc; y/pc; z/pc; vx/c; vy/c; vz/c; V^2/c^2;");
			//"boost ODEint: t*Omega; x/pc; y/pc; z/pc; vx/c; vy/c; vz/c; V^2/c^2");
	// Write header
	printer.Write("totalSimulationTime = " + to_string(totalSimulationTime*omega)
			+ " omega^-1 = " + to_string(totalSimulationTime)
			+ " pc/c, minSimSteps = " + to_string(minSimSteps)
			+ ", R = " + to_string(R) + " uG pc"
			+ ", OMEGA = " + to_string(omega));

	//TODO: this is suboptimal:
//	int boostPoints = boostTrajectory.size();
//	int outputWrites = min(outputPoints, boostPoints);
	int outputWrites = outputPoints;

//	Vec a = theoTrajectory[theoTrajectory.size()-1];
//	Vec b = boostTrajectory[boostTrajectory.size()-1];
//	cout << "boostPoints=" << boostPoints << ", outputWrites=" << outputWrites << "theo[-1]=" << a[0] << ", " << a[1] << ", " << a[2] << " boost[-1]=" << b[0] << ", " << b[1] << ", " << b[2] << endl;

//	if (outputPoints != boostPoints)
//	{
//		cout << "outputPoints:" << outputPoints << " size of boost trajectory:" << boostPoints << endl;
//	}
//	assert(outputPoints >= boostPoints);
	for(int i=0;i<outputWrites;i++)
	{
		// Init and add Time
		Vec trajectoryList;
		trajectoryList.push_back(theoTime[i]*omega); // [t*omega]=1

		// Add analytical position and energy
		AppendVector(trajectoryList, theoTrajectory[i]);
		trajectoryList.push_back(VV0);

		// Add euler position and energy
		T eulrKinE = VSquared(eulrTrajectory[i]);
		AppendVector(trajectoryList, eulrTrajectory[i]);
		trajectoryList.push_back(eulrKinE);

		// Add Runge-Kutta second order position and energy
		T ruKu2KinE = VSquared(ruKu2Trajectory[i]);
		AppendVector(trajectoryList, ruKu2Trajectory[i]);
		trajectoryList.push_back(ruKu2KinE);

		// Add Runge-Kutta fourth order position and energy
		T ruKu4KinE = VSquared(ruKu4Trajectory[i]);
		AppendVector(trajectoryList, ruKu4Trajectory[i]);
		trajectoryList.push_back(ruKu4KinE);

		/*
		// Add Runge-Kutta boost ODEint position and energy
		if (i < boostPoints)
		{
			T boostKinE = VSquared(boostTrajectory[i]);
			Scale(boostTrajectory[i], 1./Lc); // [r]=[Lc]
			AppendVector(trajectoryList, boostTrajectory[i]);
			trajectoryList.push_back(boostKinE);
		}
		else
		{
			AppendVector(trajectoryList, Vec{0,0,0,0,0,0,0,0});
		}*/

		// Write to file and advance time
		printer.Write(trajectoryList);
	}
	printer.CloseFile();
	cout << "finished." << endl;
}
