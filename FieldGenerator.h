/*
 * FieldGenerator.h
 *
 *  Created on: May 5, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "utility.h"
#include "Printer.h"
#include "RandomGen.h"

/*!
 * Function for generating and Evaluating a number of fields for statistical analysis.
 * \f$ \eta = \frac{ \langle \delta\vec{B}^2 \rangle }
 * { \langle \delta\vec{B}^2 \rangle + \vec{B}_0^2 } \f$
 */
void GenerateFields(const T eta, const T gamma, const T Lmin, const T Lmax, const T modeCount, const int fieldCount, const int gridLength, const T evalBoxLenByLmax, const int seed);


class Mode
{
public:
	Vec k; 		// [pc⁻¹] direction and mode
	Vec Axi; 	// [µG] Amplitude and polarization
	T beta; 	// phase
	Mode(Vec Axi_, Vec k_, T beta_);
	Vec ModeBfield(Vec x) const;
};

/*!
 * Class that creates a isotropic magnetic field using the harmonic method
 */
class FieldGenerator
{
private:
	RandomGen* rg;
	// stored field:
	int n; 				// number of modes (background field excluded)
	T kmin;				// [pc⁻¹]
	T kmax;				// [pc⁻¹]
	T dBvar;			// [µG]
	T B0;				// [µG] strength of the background field parallel to the z axis
	T eta;				// [1]
	T gamma;			// [1]
	T Lc;				// [pc]
	T evalBoxLenByLmax; // [1]
	std::vector<Mode> modes; // stores all the modes

	/*!
	 * Generates a powerspectrum \f$k^{-\gamma}\f$ with the Energy normalized to δB².
	 */
	void GeneratePowerspectrum();

	/*!
	 * Calculates the magnetic field values by evaluating the sum of the modes at lattice points
	 */
	void EvaluateGrid(std::vector<Vec>& s, const int gridLength, const T evalBoxLenByLmax) const;

//	/*!
//	 * Calculates the magnetic field values by evaluating the sum of the modes at randomly chosen points
//	 */
//	void EvaluateRandomPoints(const int N, std::vector<Vec>& fieldPoints);

	/*!
	 * Generating a single random mode of the field
	 * \f$ \delta B = \sum_n \sqrt{2} A_n \hat{\xi}_n \cos(k_n \hat{k}_n * x + \beta_n) \f$
	 * (Kuhlen Eq.2.44)
	 */
	Mode GenerateMode(T k, T A);

public:
	const std::vector<Mode>& GetModes() const;

	/*!
	 * Define a homogeneous constant field pointing in the z direction. Superposed with
	 * a turbulent magnetic fields with @n_ modes spanning from @kmin to @kmax logarithmically
	 * and having a variance @dBvar_. @filename specifies a filename for points evaluated to check
	 * the powerspectrum.
	 */
	FieldGenerator(const int n_, const T kmin_, const T kmax_, const T dBvar_, const T B0_, const T eta_, const T gamma_, const T Lc);

	/*!
	 * Define a homogeneous constant field pointing in the z direction.
	 */
	FieldGenerator(const T strength);

	/*!
	 * Generates or Regenerates the field
	 */
	void Generate();

	/*!
	 * calculates the field as a sum of all the modes at a given point @x
	 */
	Vec BField(const Vec& x) const;

	T GetB0() const;

	/*!
	 * Evaluate a grid and store values in a file
	 */
	void EvaluateField(const int gridLength, const T evalBoxLenByLmax, const std::string filename) const;
};
