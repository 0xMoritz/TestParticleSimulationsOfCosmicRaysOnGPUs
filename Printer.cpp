/*
 * Printer.cpp
 *
 *  Created on: Jun 9, 2022
 *      Author: Moritz Geßner
 */

#include "Printer.h"

using namespace std;


// class Printer
Printer::Printer(string filename, string header)
: file(outputPath + filename)
{
	//boost::filesystem::create_directories(outputPath); // works only with c++17
	file << header << endl;
}

int Printer::ChangeOutputPath(string newPath)
{
	if (newPath.back()=='/')
		outputPath = newPath;
	else
		outputPath = newPath + "/";

	// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
	struct stat info;

	if( stat( outputPath.c_str(), &info ) != 0 )
	{
		// No writing permission
		return 1;
	}
	else if( info.st_mode & S_IFDIR )
	{
		return 0;
	}
	else
	{
		// Directory does not exist
		return 1;
	}
}
string Printer::GetOutputPath()
{
	return outputPath;
};
void Printer::Write(string str)
{
	file << str << endl;
}
void Printer::Write(const Vec& vec)
{
	for (size_t i=0;i<vec.size()-1;i++)
	{
		file << vec[i] << ";";
	}
	file << vec[vec.size()-1] << endl;
}
void Printer::CloseFile()
{
	file.close();
}
Printer::~Printer()
{

}
