/*!
 * @file CosmicRayTransport.h
 * @brief Main source file, parses input to execute a specific command
 *
 *  Created on: May 2, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "global.h"

/*!
 * Parse the command line arguments \p argv to determine execution mode
 * and parameters to execute corresponding simulation or function.
 * \param argc is the number of command line parameters including the command that executes the program
 * \param argv is the list of parameters, \p argv[0] will be the command to execute the program
 * \return 1 if reading errors occur, 0 for nominal execution
 */
int ParseInput(int argc, char*argv[]); // Forward declaration not necessary



/*! \mainpage Test particle simulation of cosmic rays on GPUs
 *
 * \section intro_sec Introduction
 *
 * This program's main purpose is the simulation of test particles in isotropically turbulent magnetic fields,
 * by numerically integrating the equations of motion using a Nvidia GPU. Moreover it can be used to perfom said
 * simulation using the CPU, run a test of different numerical integrators and generate isotropic turbulent
 * fields as well as evaluate those fields in a grid for further analysis.
 *
 * The program was developed during my bachelor's thesis "test particle simulations of cosmic rays on GPUs"
 * supervised by Prof. Philipp Mertsch and Dr. Vo Hong Minh Phan at RWTH University.
 *
 * Github page: https://github.com/0xMoritz/CosmicRaysOnGPUs/tree/main/CosmicRayTransport
 *
 * For questions, errors or bugs concerning program or thesis: moritz.gessner@rwth-aachen.de
 *
 * \section usage_sec Usage
 *
 * \subsection prerequisites Prerequisites
 * 		The integration of the particle's equations of motion are based on the boost Library ODEint,
 * 		which utilizes thrust. Both must be installed.
 * 		For compilation the NVCC cuda compiler is required. To yield an advantage with the execution on the GPU,
 * 		the amount of particles should be at least \f$\sim 10^3\f$.
 * 		Generating the documentation requires latex, ghostscript and Doxygen.
 * 		For running on the CPU, the comments in the file @p global.h have to be modified.
 *
 * \subsection compiling Compiling the program
 *      Calling @p make @p CosmicRayTransport.bin should build the program. To change between using GPU and CPU the
 *      comments in @p global.h must be changed manually. I lacked the time to introduce fancy compiler flags
 *      for this. But it should be easy to catch up on this and make it easier to switch.
 *
 *      With @p make @p man, the doxygen manual will be generated as a latex pdf file and html site. They will be
 *      found in @p Documentation/latex/refman.pdf and @p Documentation/html/index.html respectively.
 *
 * \subsection running Running the program
 * 		The program has three modes of execution:
 *
 * 		"i" for simulating particle trajectories in an isotropic
 * 		turbulent magnetic field. The amount of different field realizations and particles can be specified.
 * 		The particle trajectories are saved in individual files in either logarithmic or linear time steps.
 *
 * 		"b" simulates the motion of a single particle in a background magnetic field without turbulence. The
 *		numerical integration is performed with the Euler method, a Runge Kutta method 2. order, a Runge
 *		Kutta method 4. order using the conventional choice of coefficients and a Runge Kutta method as
 *		implemented by the boost library ODEint.
 *
 *      "f" can be used to generate a turbulent magnetic field and evaluate it in a grid box. The data is
 *      saved to a file and can be used for further analysis e.g. using a discrete Fourier transform...
 *
 * \section state State of the project
 *      The generation of fields has been thoroughly tested in the thesis. The core function of simulating test
 *      particles in isotropic turbulence is working properly.
 *
 *      However, some aspects can be improved:
 *          - The conventional implementation of
 *      Runge Kutta fourth order performs worse than the implementation by the library. This might be because
 *      of the choice of coordinates but there might be more to it.
 *          - Still from the process of developement, some passing of parameters is inefficient or unnecessary. One
 *      might frown upon the parameters needed for some functions and is encouraged to do some optimization him/
 *      herself.
 *          - Changing between GPU and CPU is so far only possible by changing comments in global.h, this could be
 *      much nicer with a compiler directive.
 *          - Parallelization so far only takes places by simulating many particles in parallel on the same magnetic
 *      field. Analysis however shows, that for particle numbers of \f$ \sim 10^3 \f$ the GPU's resources are
 *      not exhausted. The performance could be increased by also parallelizing the computation of multiple
 *      fields. Yet this will be cumbersome, as the states for different magnetic fields must be nested in
 *      thrust tuples.
 *          - Currently the trajecoty of every particle is saved as an individual file, this requires a lot of
 *      space. For the calculation of
 *      statistical properties such as the diffusion coefficients, it is more interesting to consider average
 *      over many particles and field realizations. The averaging process is currently done with python scripts
 *      but it should be more effective to automate part of the process in this program, e.g. save the mean
 *      quadratic displacement only.
 */
