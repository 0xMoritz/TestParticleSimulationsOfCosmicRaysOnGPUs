/*
 * FieldGenerator.cpp
 *
 *  Created on: May 5, 2022
 *      Author: moritz
 */

#include "FieldGenerator.h"

using namespace std;


void GenerateFields(const T eta, const T gamma, const T Lmin, const T Lmax, const T modeCount, const int fieldCount, const int gridLength)
{
	// process parameters
	const T kmin = 2*M_PI/Lmax; 			// [pc⁻¹]
	const T kmax = 2*M_PI/Lmin; 			// [pc⁻¹]
	const T Lc = Lmax/5; 					// [pc] approximation for correlation length Kuhlen p. 66, Harari et al. eq.2.4ff
	const T B0 = 4; 						// [µG] strength of background magnetic field
	const T dBvar = eta/(1.-eta)*B0*B0; 	// [µG²]
	FieldGenerator field(modeCount, kmin, kmax, dBvar, B0, eta, gamma, Lc);

	// Output parameters
	cout << "Parameters:"
			<< "\n\teta = " << eta
			<< "\n\tmodeCount = " << modeCount
			<< "\n\tfieldCount = " << fieldCount
			<< "\n\tgridLength = " << gridLength
			<< "\n\tLmin = " << Lmin
			<< "\n\tLmax = " << Lmax
			<< "\n\tkmin = " << kmin
			<< "\n\tkmax = " << kmax
			<< "\n\tLc = " << Lc
			<< "\n\tB0 = " << B0
			<< "\n\tdeltaBvar = " << dBvar << endl;

	// Generate many fields for statistics
	for (int i = 0; i < fieldCount;i++)
	{
		field.Generate();
		field.EvaluateField(gridLength, "field_" + to_string(i) + ".csv");
	}
}

Mode::Mode(Vec Axi_, Vec k_, T beta_)
{
	Axi = Axi_;
	k = k_;
	beta = beta_;
}
Vec Mode::ModeBfield(Vec x) const
{
	Vec res(Axi);
	Scale(res, std::cos(ScalarProd(k, x) + beta));
	return res;
}

const vector<Mode>& FieldGenerator::GetModes() const
{
	return modes;
}

Vec FieldGenerator::BField(const Vec& x) const
{
	Vec res = {0,0,0};
	for(int i=0;i<n;i++)
	{
		VecAdd(res, modes[i].ModeBfield(x));
	}
	res[2] += B0; // Add homogeneous background field
	return res;
}
T FieldGenerator::GetB0() const
{
	return B0;
}
void FieldGenerator::GeneratePowerspectrum()
{
	// Normalize spectrum:
	// g(k)= gFac k^{-5/3} // power law
	// var = delta B^2 = int_RR g(k) dk

	// Logarithmic scaling with N steps
	int N = (int)1e6;
	T logFac = pow(kmax/kmin, 1./(N-1)); // [1]
	//T logFac = pow(kmax/kmin, 1./N); // [1]
	T k = kmax;//kmin;
	T gSum = 0;
	for (int i=0;i<N;i++) // Numerical integral
	{
		T dk = (logFac-1)*k;
		gSum += pow(k, -gamma) * dk; // [pc^(γ-1)]
		k *= logFac;
	}
	T gFac = dBvar/gSum; // [µG²/pc^(γ-1)]

	/* Logarithmic scaling of k:
	 * k(i) = A*logFac^i
	 * kmin = k(0) = A  => A = kmin
	 * kmax = k(n-1) = A logFac^(n-1)  => logFac = (kmax/kmin)^(1/(n-1))
	 */
	//logFac = pow(kmax/kmin, 1./(n-1)); // Introduce logarithmic scaling again with number of modes n
	//logFac = pow(kmax/kmin, 1./n); // Introduce logarithmic scaling again with number of modes n
	k = (kmin + kmax)/2.;
	Vec veck(n); // [pc⁻¹]
	Vec vecA(n); // [µG]
	Vec vecg(n); // [µG²·pc]
	modes.clear();
	modes.reserve(n);
	for(int i=0;i<n;i++)
	{
		T dk = kmax-kmin;//(logFac-1)*k;
		cout << "dk=" << dk << endl;
		veck[i] = k;
		vecg[i] = gFac * pow(k, -gamma);
		vecA[i] = sqrt( 2.*pow(k, -gamma)*gFac*dk); // <cos²> = 1/2
		// Apply Amplitude to mode
		modes.push_back(GenerateMode(k, vecA[i]));
		//k *= logFac;
	}

//	// Write powerspectrum:
//	Printer printer("powerspectrum.csv", "k*pc; g(k)/microG^2*pc; A(k)/microG");
//	for(int i=0; i<n; i++)
//	{
//		Vec v = {veck[i], vecg[i], vecA[i]};
//		printer.Write(v);
//	}
//	printer.CloseFile();
}
void FieldGenerator::EvaluateGrid(vector<Vec>& fieldPoints, const int gridLength) const
{
	T max = 2*M_PI/kmin; // choose 2pi/kmin to incorporate each mode
	T dr = max/gridLength;
	fieldPoints.reserve((size_t)max/dr*max/dr*max/dr);
	for (T x=0; x<max; x+=dr)
	{
		cout << "." << flush;
		for (T y=0; y<max; y+=dr)
		{
			for (T z=0; z<max; z+=dr)
			{
				Vec pos = {x, y, z};
				Vec B = BField(pos);
				AppendVector(pos, B);
				fieldPoints.push_back(pos);
			}
		}
	}
}
//void FieldGenerator::EvaluateRandomPoints(const int N, vector<Vec>& fieldPoints)
//{
//	random_device rd;  // Will be used to obtain a seed for the random number engine
//	mt19937 RdmGen(rd());
//	T max = 2*M_PI/kmin; // choose 2pi/kmin to incorporate each mode
//	uniform_real_distribution<T> dist(0., max);
//
//	fieldPoints.reserve((size_t)N);
//
//	for (int i=0;i<N;i++)
//	{
//		T x = dist(RdmGen);
//		T y = dist(RdmGen);
//		T z = dist(RdmGen);
//
//		Vec pos = {x, y, z};
//		Vec B = BField(pos);
//		AppendVector(pos, B);
//		fieldPoints.push_back(pos);
//	}
//}
Mode FieldGenerator::GenerateMode(T k, T A)
{
	// Generate a random direction distributed isotropically (adopted from Kuhlen Eq.2.45)
	// Generate random unitvector and scale by k:
	Vec vec_k(3);
	T phi = rg->RandomFloat_0_2PI();
	T eta = rg->RandomFloat_m1_1();
	vec_k[0] = sqrt(1-eta*eta) * cos(phi);
	vec_k[1] = sqrt(1-eta*eta) * sin(phi);
	vec_k[2] = eta;
	Scale(vec_k, k);
	// Generate the orthogonal polarization vector from k (adopted from Kuhlen Eq.2.48) and Amplitude A
	// polarization has to be orthogonal for field to be divergence free
	Vec vec_Axi(3);
	T alpha = rg->RandomFloat_0_2PI();
	vec_Axi[0] = A*(-sin(alpha)*sin(phi) + cos(alpha)*cos(phi)*eta);
	vec_Axi[1] = A*( sin(alpha)*cos(phi) + cos(alpha)*sin(phi)*eta);
	vec_Axi[2] = A*(-sqrt(1-eta*eta)*cos(alpha));
	// Generate a random Phase distributed uniformally in [0,2 PI]
	T beta = rg->RandomFloat_0_2PI();
	// DebugInfo (Check whether vectors are normed correctly): cout << "k=" << k << " Norm |Axi|^2/A^2=" << SqNorm(vec_Axi)/A/A << " Norm |k|^2/k^2=" << SqNorm(vec_k)/k/k << " Scalarproduct vec_Axi·vec_k=" << ScalarProd(vec_Axi, vec_k) << endl;
	// concatenate properties in mode object
	Mode mode(vec_Axi, vec_k, beta);
	return mode;
}
FieldGenerator::FieldGenerator(const int n_, const T kmin_, const T kmax_, const T dBvar_, const T B0_, const T eta_, const T gamma_, const T Lc_)
{
	rg = RandomGen::GetInstance();
	n = n_;
	kmin = kmin_;
	kmax = kmax_;
	dBvar = dBvar_;
	B0 = B0_;
	eta = eta_;
	gamma = gamma_;
	Lc = Lc_;
}
void FieldGenerator::Generate()
{
	// Generate the Field:
	PrintTime();
	cout << "Generating field ..." << flush;
	clock_t beginComputingTime = clock(); // https://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
	// Generate the power spectrum
	GeneratePowerspectrum();
	cout << " finished. Time elapsed: " << float(clock() - beginComputingTime) /  CLOCKS_PER_SEC << endl;
}

void FieldGenerator::EvaluateField(const int gridLength, const string filename) const
{
	// Evaluate field:
	PrintTime();
	cout << "Evaluating field at gridpoints" << flush;
	clock_t beginComputingTime = clock();
	// Read at grid points and save to file for evaluation of the power spectrum in python
	vector<Vec> fieldPoints; // Format (x, y, z, B_x, B_y, B_z)^T
	EvaluateGrid(fieldPoints, gridLength);
	//EvaluateRandomPoints((int)1e5, fieldPoints);
	float timeElapsedInSeconds = float(clock() - beginComputingTime) /  CLOCKS_PER_SEC;
	cout << " finished. Time elapsed: " << timeElapsedInSeconds << endl;

	// Print evaulated grid:
	PrintTime();
	//Printer printer(filename, "x/Lc; y/Lc; z/Lc; B_x/µG; B_y/µG; B_z/µG");
	Printer printer(filename, "x/pc; y/pc; z/pc; B_x/µG; B_y/µG; B_z/µG");
	printer.Write("eta = " + to_string(eta) + ", modeCount = " + to_string(n) + ", kmin = " + to_string(kmin)
			+ ", kmax = " + to_string(kmax) + ", gridLength = " + to_string(gridLength)
			+ ", Lc = " + to_string(Lc) + ", timeElapsedInSeconds = " + to_string(timeElapsedInSeconds));
	for(size_t i=0; i<fieldPoints.size(); i++)
	{
		/*// Norm position to coherence length and write TODO: norm or not?
		fieldPoints[i][0] = fieldPoints[i][0]/Lc;
		fieldPoints[i][1] = fieldPoints[i][1]/Lc;
		fieldPoints[i][2] = fieldPoints[i][2]/Lc;*/
		printer.Write(fieldPoints[i]);
	}
	printer.CloseFile();
}

FieldGenerator::FieldGenerator(const T B0_) // Generate constant
{
	n = 0;
	// The following variables will not be used
	kmin = 0;
	kmax = 0;
	dBvar = 0;
	B0 = B0_;
	Lc = 0;
	eta = 0;
	// modes stays empty but will not be called because n=0
}

