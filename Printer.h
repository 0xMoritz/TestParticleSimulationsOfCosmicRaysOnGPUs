/*
 * Printer.h
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "global.h"


/*!
 * Object for writing content to a file. A Printer Object is bound to one file,
 * that should be closed before destruction.
 */
class Printer
{
private:
	static std::string outputPath;
	std::ofstream file;

public:
	/*!
	 * Constructor for Printer Object, opens a file with @filename. @header will
	 * be written to the first line of the file.
	 */
	Printer(std::string filename, std::string header="");

//	Printer(const Printer& p);

	/*!
	 * Changes the global output Path, if the directory doesn't exist or can't be accessed, 1 will be returned.
	 * If the directory exist and is writeable, 0 will be returned
	 */
	static int ChangeOutputPath(std::string newPath);

	/*!
	 * Returns the global output Path
	 */
	static std::string GetOutputPath();

	/*!
	 * Writes a custom message to the file. Use for Header etc.
	 */
	void Write(std::string str);

	/*!
	 * Writes a @Vec vector to the file in csv format with ";" delimiter.
	 */
	void Write(const Vec& vec );

	/*!
	 * Closes the file.
	 */
	void CloseFile();

	/*!
	 * Does not call @CloseFile()
	 */
	~Printer();
};
