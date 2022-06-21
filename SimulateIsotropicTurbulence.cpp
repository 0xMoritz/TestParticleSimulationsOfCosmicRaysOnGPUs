/*
 * SimulateIsotropicTurbulence.cpp
 *
 *  Created on: May 18, 2022
 *      Author: moritz
 */

#include "SimulateIsotropicTurbulence.h"

using namespace std;


void IsotropicTurbulenceSimulation(const T R_norm, const T eta, const T gamma, const T Lmin, const T Lmax, const int modeCount, const int fieldCount, const int particlesPerFieldCount, const T simulationTimeInGyroperiods, const int minSimSteps, const int outputPoints, const bool useBoost, const int seed)
{
	// Parameter list (units pc,µG,c):
	// Field related parameters // TODO: merge this parameter set with FieldGenerator.cpp: SimulateBackgroundField
	const T kmin = 2*M_PI/Lmax; 		// [pc⁻¹]
	const T kmax = 2*M_PI/Lmin; 		// [pc⁻¹]
//	cout << "L_max=" << Lmax << endl;	// should be ~ 150pc Kuhlen p.66
	// Chosen units:
	const T Lc = Lmax/5; 				// [pc] approximation for correlation length Kuhlen p. 66, Harari et al. eq.2.4ff
	const T B0 = 4; 					// [µG] strength of background magnetic field
	const T V = 1;						// [c] velocity
	// Magnetic field
	const T dBvar = eta/(1.-eta)*B0*B0;	// [(µG)²] variance of turbulent field
	//cout << "dBvar = " << dBvar << endl;
	const T B_mean = sqrt(B0*B0+dBvar);	// [µG] mean B field
	// Particle related parameters
	const T R = R_norm*B_mean*Lc;		// [µG·pc] rigidity of particle
	const T Q_M = V/R;					// [c·pc⁻¹·µG⁻¹] charge divided by mass
	const T omega = Q_M*B_mean;			// [c·pc⁻¹] gyration frequency Ω for the mean Field
//	cout << "OMEGA=" << OMEGA << endl;

	// Composite parameters
	const int stepsPerOutput = minSimSteps % outputPoints==0 ? minSimSteps/outputPoints : minSimSteps/outputPoints+1;
	const int actualSimSteps = stepsPerOutput*outputPoints;
	const T totalSimulationTime = simulationTimeInGyroperiods*2*M_PI/omega;	// [pc·c⁻¹] time simulated in the simulation
	const T dt = totalSimulationTime/actualSimSteps; // [pc·c⁻¹]
	assert(minSimSteps >= outputPoints);

	// Print parameters to console and file
	string parameterText = "Input Parameters:\n"
				"\tseed = " + to_string(seed) + "\n"
				+ "\tR_norm = " + to_string(R_norm) + "\n"
				+ "\teta = " + to_string(eta) + "\n"
				+ "\tgamma = " + to_string(gamma) + "\n"
				+ "\tmodeCount = " + to_string(modeCount) + "\n"
				+ "\tfieldCount = " + to_string(fieldCount) + "\n"
				+ "\tparticlePerFieldCount = " + to_string(particlesPerFieldCount) + "\n"
				+ "\ttotalSimulationTime = " + to_string(totalSimulationTime*omega) + " omega^-1 = " + to_string(totalSimulationTime) + " pc c^-1\n"
				+ "\tminSimSteps = " + to_string(minSimSteps) + "\n"
				+ "\toutputPoints = " + to_string(outputPoints) + "\n"
				+ "\tuseBoost = " + to_string(useBoost) + "\n"
				+ "\tLmin = " + to_string(Lmin) + " pc\n"
				+ "\tLmax = " + to_string(Lmax) + " pc\n"
					+ "\n"
				+ "\tR = " + to_string(R) + " uG pc\n"
				+ "\tOMEGA = " + to_string(omega) + " c pc\n"
				+ "\tLc = " + to_string(Lc) + " pc\n"
				+ "\tdt = " + to_string(dt*omega) + " OMEGA^-1 = " +to_string(dt) + " pc/c\n"
				+ "\tactualSimSteps = " + to_string(actualSimSteps) + "\n";
	cout << parameterText;
	Printer parameterFilePrinter("_info.txt");
	parameterFilePrinter.Write(parameterText); // Not closing here, last entry will be computation time

	// Generate the static, isotropic, turbulent field
	FieldGenerator field(modeCount, kmin, kmax, dBvar, B0, eta, gamma, Lc);


	// Measure total computation time:
	const clock_t beginComputingTime = clock();

	// Simulate a number of test particles with velocity v
	for (int batchNo=0; batchNo<fieldCount; batchNo++)
	{
		// Generate new Field
		field.Generate();

		// Instantiate Engine
		Engine* engine;

		if (useBoost)
		{
			// Simulation using Boost ODEint
			engine = new BoostSolver(particlesPerFieldCount, totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, field, omega, Lc);
		}
		else
		{
			// Simulation with Runge Kutta method fourth order
			engine = new RungeKuttaSolver(particlesPerFieldCount, totalSimulationTime, outputPoints, stepsPerOutput, dt, R, V, field, omega, 4);
		}
		vector<Vec> trajectory; trajectory.reserve(outputPoints);
		Vec time; time.reserve(outputPoints);
		// Simulation Call
		float computingTime = engine->Simulation(trajectory, time, batchNo);
	}

	// output total computation time and write to file
	float timeElapsedInSeconds = float(clock() - beginComputingTime) /  CLOCKS_PER_SEC;
	cout << endl << "\nTotal computation time: " << timeElapsedInSeconds << " s" << endl;
	parameterFilePrinter.Write("Total computation time = " + to_string(timeElapsedInSeconds) + " s");
	parameterFilePrinter.CloseFile();
}
