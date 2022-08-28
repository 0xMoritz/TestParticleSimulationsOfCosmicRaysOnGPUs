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
