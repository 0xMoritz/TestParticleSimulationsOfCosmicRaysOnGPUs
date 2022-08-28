/*!
 * @file FieldGenerator.h
 * @brief Generates an isotropic turbulent magnetic
 *
 *  Created on: May 5, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "utility.h"
#include "Printer.h"
#include "RandomGen.h"

/*! @brief Function for generating and Evaluating a number of fields for statistical analysis
 *
 * Generates a number of magnetic fields with a background component parallel to \f$z\f$ and an isotropic
 * turbulent component constructed by the harmonic methods. Evaluates the field at grid points and writes the
 * values to files.
 * \param eta turbulence level defined as \f$ \eta = \frac{ \langle \delta\vec{B}^2 \rangle } { \langle \delta\vec{B}^2 \rangle + \vec{B}_0^2 } \f$
 * \param gamma spectral index of the power spectrum \f$ \delta \vec{B}^2(\vec{k})\propto k^{-\gamma} \f$ with which the isotropic field is constructed
 * \param Lmin corresponds to the upper bound of the power spectrums validity \f$k_{\max} = \frac{2\pi}{L_{\min}}\f$
 * \param Lmax corresponds to the lower bound and injection scale of the power spectrums validity \f$k_{\min} = \frac{2\pi}{L_{\max}}\f$
 * \param modeCount amount of modes \f$n\f$ in the discretized powerspectrum
 * \param fieldCount count of different field realizations
 * \param gridLength number of evaluations \f$N\f$ along one axis. Total amount of evaluated grid points will be \f$N^3\f$
 * \param evalBoxLenByLmax length of one side of the evaluation box given in multiples of \p Lmax.
 *  choose multiples of \p Lmax to increase the resolution on a fourier transform of the evaluated grid
 * \param seed seed for the random numbers generated for the construction of the modes
 */
void GenerateFields(const T eta, const T gamma, const T Lmin, const T Lmax, const T modeCount, const int fieldCount, const int gridLength, const T evalBoxLenByLmax, const int seed);

/*!
 * @brief Container to hold information about one magnetic mode
 *
 * In the harmonic method, the field is constructed from plane waves \f$ B_i(\vec{r}) = A_i(k_i) \hat{\xi}_i
 * \cos\left(k_i\hat{k}_i\cdot\vec{r}+\beta_i\right) \f$, one of these modes has a wave vector
 * \f$\vec{k}_i=k_i\hat{k}_i \f$ with a random direction \f$\hat{k}_i\f$, a corresponding amplitude \f$A_i\f$
 * with a random polarizatiion \f$ \hat{\xi}_i \f$ and phase \f$\beta_i\f$.
 */
class Mode
{
public:
	Vec k; 	 	//!< \f$ \left[\mathrm{pc}^{-1}\right] \f$ Wave vector of the mode
	Vec Axi; 	//!< \f$ \left[\mu\mathrm{G}\right] \f$ Amplitude and polarization
	T beta;		//!< Phase of the mode
	/*! @brief Constructs a mode object given all properties
	 */
	Mode(Vec Axi_, Vec k_, T beta_);
	/*! @brief Calculates the partial magnetic field of this wave
	 *
	 * Calculates the value of the magnetic field for this wave according to \f$ B_i(\vec{r}) = A_i(k_i) \hat{\xi}_i \cos\left(k_i\hat{k}_i\cdot\vec{r}+\beta_i\right) \f$.
	 */
	Vec ModeBfield(Vec x) const;
};

/*! @brief Class that creates a isotropic magnetic field using the harmonic method supoerposed with a background field
 *
 * Using the harmonic method, @ref FieldGenerator can generate a discretized model of an isotropic turbulent magnetic field. The field can be evaluated by the formula \f$\delta \vec{B}(\vec{r}) = \sum_{j=0}^{n-1} A_j \hat{\xi}_j \cos\left(k_j\hat{k}_j\cdot\vec{r}+\beta_j\right)\f$.
 * Information about the modes is stored in an internal vector of @ref Mode objects.
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
	T B0;					// [µG] strength of the background field parallel to the z axis
	T eta;				// [1]
	T gamma;			// [1]
	T Lc;					// [pc]
	T evalBoxLenByLmax; // [1]
	std::vector<Mode> modes; // stores all the modes

	/*!
	 * Generates a powerspectrum \f$k^{-\gamma}\f$ with the Energy normalized to δB² i.e. the variance of the turbulence.
	 * Normalize the powerspectrum and calculate the list of amplitudes corresponding to the given list of wavenumbers
	 */
	void GeneratePowerspectrum();

	/*!
	 * Calculates the magnetic field values by evaluating the sum of the modes at lattice points.
	 * Writes values to @p s, in the format vector<(x, y, z, B_x, B_y, B_z)>
	 */
	void EvaluateGrid(std::vector<Vec>& s, const int gridLength, const T evalBoxLenByLmax) const;

	/*!
	 * Generating a single random mode of the field
	 */
	Mode GenerateMode(T k, T A);

public:
	/*! @brief Retrieve the private vector of modes that constitute the field
	 *
	 * \return a const reference to the vector of @Mode objects that make up the field. To feed the field to the numerical integrator on the GPU, it needs to be accessible via pointers. It would be nicer to have an interface over @ref FieldGenerator but extracting the modes like this and converting them to pointers in @ref LorentzForce is a quick work-around that also does the trick
	 */
	const std::vector<Mode>& GetModes() const;

	/*! @brief Define properties of the homogeneous and turbulent part of the magnetic field
	 *
	 * Define a homogeneous constant field pointing in the z direction. Superposed with
	 * a turbulent magnetic fields with @n_ modes spanning from @kmin to @kmax logarithmically
	 * and having a variance @dBvar_. @filename specifies a filename for points evaluated to check
	 * the powerspectrum.
	 */
	FieldGenerator(const int n_, const T kmin_, const T kmax_, const T dBvar_, const T B0_, const T eta_, const T gamma_, const T Lc);

	/*! @brief Define a homogeneous constant field pointing in the z direction.
	 *
	 * This constructor exists for the sole purpose of creating a purely homogeneous field that has no isotropic component. This is used to test the different numerical solvers in a controlled environment, thus the only call of this function is in @ref SimulateBackgroundField.
	 */
	FieldGenerator(const T strength);

	/*! @brief Generates or regenerates the field
	 *
	 * Using the parameters given in the constructor @ref FieldGenerator, this member creates a new set of modes that satisfy the given constrains.
	 */
	void Generate();

	/*! @brief Calculates the field as a sum of all the modes at a given point @p x
	 *
	 * Calculates the value of the magnetic field at the point @p x using \f$\delta \vec{B}(\vec{r}) = \sum_{j=0}^{n-1} A_j \hat{\xi}_j \cos\left(k_j\hat{k}_j\cdot\vec{r}+\beta_j\right)\f$
	 */
	Vec BField(const Vec& x) const;

	/*! @brief retrieve the strength of constant component of the magnetic field
	 *
	 */
	T GetB0() const;

	/*! @brief Evaluate a grid and store values in a file
 	 *
	 * Evaluate the magnetic field on a grid of size @p evalBoxLenByLmax with \p gridLength points per axis and store them to @p filename
	 *
	 * \param gridLength number of evaluations \f$N\f$ along one axis. Total amount of evaluated grid points will be \f$N^3\f$
	 * \param evalBoxLenByLmax length of one side of the evaluation box given in multiples of \p Lmax.
	 *  choose multiples of \p Lmax to increase the resolution on a fourier transform of the evaluated grid
	 * \param filename name of the file to save the field values in
	 */
	void EvaluateField(const int gridLength, const T evalBoxLenByLmax, const std::string filename) const;
};
