//============================================================================
// Name        : CosmicRayTransport.cpp
// Author      : Moritz Geßner (moritz.gessner@rwth-aachen.de)
// Copyright   : ©2022 Moritz Geßner All rights reserved
// Description : Main File for the simulation of cosmic rays propagating
//				 through homogeneous or turbulent magnetic fields as well as
//				 isotropic turbulent field construction for CPUs. The state of
//				 particles is stored as a Vec q=(x, y, z, v_x, v_y, v_z)^T.
//============================================================================

#include "CosmicRayTransport.h"
#include "SimulateBackgroundField.h"
#include "SimulateIsotropicTurbulence.h"
//#include "FieldGenerator.h"
//#include "Engines.h"

using namespace std; //TODO: gamma input with division

string Printer::outputPath{"output/"};


int ParseInput(int argc, char*argv[])
/**
 * @argc is the number of command line parameters including the command that executes the program
 * @argv[] is the list of parameters, @argv[0] will be the command to execute the program
 */
{
	// Standard parameters:
	T R_norm = 0.9;
	T eta = 0.1;
	T gamma = 5./3;
	int modeCount = 100;
	int fieldCount = 30; // Different fields to simulate particles on
	int particlesPerFieldCount = 1 << 10; // 2^10
	int fieldGenerationCount = 100;
	T simulationTimeInGyroperiods = 1e3;
	int minSimSteps = (int)2e5;
	int outputPoints = (int)1e4;
	int gridLength = 32;
	bool useBoost = 1;
	int seed = 0;
	T Lmin = 1;							// [pc] inner scale of turbulence
	T Lmax = 150;						// [pc] outer scale of turbulence
	bool logTime = 1;
	T evalBoxLenByLmax = 4;

	// Help messages
	string errorMsg = "[ERROR] invalid mode input. Type 'h' as argument for help.";
	string I_help = "Mode 'i' simulates a series of particles in an isotropic field\n"
			"\t./command i [outputPath] [seed] [R_norm] [eta] [gamma] [modeCount] [fieldCount] "
			"\n\t[particlePerFieldCount] [SimulationTime/gyroperiods] [minSimSteps] "
			"\n\t[outputPoints] [useBoost (1 or 0)] [Lmin/pc] [Lmax/pc] [logTime (1 or 0)]\n";
	string B_help = "Mode 'b' simulates one particle in a background field using various methods\n"
			"\t./command b [outputPath] [R_norm] [SimulationTime/gyroperiods] [minSimSteps] [outputPoints]\n";
	string F_help = "Mode 'f' generates a number of isotropic fields and evaluates these at grid points as output\n"
			"\t./command f [outputPath] [seed] [eta] [gamma] [modeCount] [GenerationCount]\n"
			"\t[gridLength] [Lmin/pc] [Lmax/pc] [evalBoxLen/Lmax]";
	string H_help = "Simulation of cosmic rays propagating through homogeneous or turbulent magnetic fields as well\n"
			"as isotropic turbulent field construction for CPUs.\n"
			"Use './command h' to bring up this help menu, "
			"or '.command h [mode]' to get information about operation mode [mode]\n";

	// Empty line for better visibility:
	cout << endl;
	// Read input
	if (argc < 2)
	{
		cout << errorMsg << endl << endl;
		return 1;
	}

	if (argv[1][0] == 'H' || argv[1][0] == 'h')
	{
		cout << "Using C++" << __cplusplus << endl;

		if (argc > 2)
		{
			if (argv[2][0] == 'I' || argv[2][0] == 'i')
			{
				cout << I_help << endl;
				return 0;
			}
			if (argv[2][0] == 'B' || argv[2][0] == 'b')
			{
				cout << B_help << endl;
				return 0;
			}
			if (argv[2][0] == 'F' || argv[2][0] == 'f')
			{
				cout << F_help << endl;
				return 0;
			}
			cout << "Error: Mode '" << argv[2][0] << " is not available." << endl;
			return 1;
		}
		else
		{
			cout << H_help << I_help << B_help << F_help << endl;
			return 0;
		}
	}
	// Determine filepath
	if (argc > 2)
	{
		string path = (string)argv[2];
		int dirCheck = Printer::ChangeOutputPath(path);
		if (dirCheck == 1)
		{
			cout << "Error: Directory '" << path << "' access denied." << endl;
			return 1;
		}
	}
	// Choose mode
	if (argv[1][0] == 'I' || argv[1][0] == 'i')
	{
		if (argc > 3)
			seed = atoi(argv[3]);
		if (argc > 4)
			R_norm = strtod(argv[4], nullptr);
		if (argc > 5)
			eta = strtod(argv[5], nullptr);
		if (argc > 6)
			gamma = strtod(argv[6], nullptr);
		if (argc > 7)
			modeCount = (int) strtod(argv[7], nullptr);
		if (argc > 8)
			fieldCount = (int) strtod(argv[8], nullptr);
		if (argc > 9)
			particlesPerFieldCount = (int) strtod(argv[9], nullptr);
		if (argc > 10)
			simulationTimeInGyroperiods = strtod(argv[10], nullptr);
		if (argc > 11)
			minSimSteps = (int) strtod(argv[11], nullptr);
		if (argc > 12)
			outputPoints = (int) strtod(argv[12], nullptr);
		if (argc > 13)
			useBoost = atoi(argv[13]);
		if (argc > 14)
			Lmin = strtod(argv[14], nullptr);
		if (argc > 15)
			Lmax = strtod(argv[15], nullptr);
		if (argc > 16)
			logTime = atoi(argv[16]);

		RandomGen* rd = RandomGen::GetInstance();
		rd->ChangeSeed(seed);
		IsotropicTurbulenceSimulation(R_norm, eta, gamma, Lmin, Lmax, modeCount, fieldCount, particlesPerFieldCount, simulationTimeInGyroperiods, minSimSteps, outputPoints, useBoost, seed, logTime);
	}
	else if (argv[1][0] == 'B' || argv[1][0] == 'b')
	{
		if (argc > 3)
			R_norm = strtod(argv[3], nullptr);
		if (argc > 4)
			simulationTimeInGyroperiods = strtod(argv[4], nullptr);
		if (argc > 5)
			minSimSteps = atoi(argv[5]);
		if (argc > 6)
			outputPoints = atoi(argv[6]);

		SimulateBackgroundField(R_norm, simulationTimeInGyroperiods, minSimSteps, outputPoints);
	}
	else if (argv[1][0] == 'F' || argv[1][0] == 'f')
	{
		if (argc > 3)
			seed = atoi(argv[3]);
		if (argc > 4)
			eta = strtod(argv[4], nullptr);
		if (argc > 5)
			gamma = strtod(argv[5], nullptr);
		if (argc > 6)
			modeCount = strtod(argv[6], nullptr);
		if (argc > 7)
			fieldGenerationCount = atoi(argv[7]);
		if (argc > 8)
			gridLength = atoi(argv[8]);
		if (argc > 9)
			Lmin = strtod(argv[9], nullptr);
		if (argc > 10)
			Lmax = strtod(argv[10], nullptr);
		if (argc > 11)
			evalBoxLenByLmax = strtod(argv[11], nullptr);

		RandomGen* rd = RandomGen::GetInstance();
		rd->ChangeSeed(seed);
		GenerateFields(eta, gamma, Lmin, Lmax, modeCount, fieldGenerationCount, gridLength, evalBoxLenByLmax, seed);
	}
	else
	{
		cout << errorMsg << endl << endl;
		return 1;
	}

	return 0;
}

/**if (argc == 2 && ((string)argv[1]).find("h") != string::npos) // If first argument contains 'h'

				<< " if not specified standard input is " << mode << " " << R << " " << eta << " " << modeCount << " " << amount
				<< " " << simulationTimeInGyrationCycles << " " << minSimSteps << " " << outputPoints << endl;*/

int main(int argc, char* argv[])
{
	// Initialize RandomGen Singleton
	RandomGen::GetInstance();
	return ParseInput(argc, argv);
}

//TODO: compute power spectrum, use logarithmically scaled concentric binsx2go
