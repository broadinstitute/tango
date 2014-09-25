//========================================================================
// Project     : PrepDB
// Name        : Parameter.h
// Author      : Xiao Yang
// Created on  : Oct 22, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2012)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or
//				 guaranteed support whatsoever. The Broad Institute cannot
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================


#ifndef PARAMETER_H_
#define PARAMETER_H_


#include "xutil.h"

/* input parameters */
class Parameter{

public:
	std::string idbdir; // input DB dir
	std::string irdir; // input read data dir
	int p;
	int k;
	double perc;    // downsampling percentage
 	bool silent;

	Parameter (int argc, char** argv): argnum(argc), arg(argv){
		init ();
		for (int i = 1; i < argnum; i += 2) {
			std::string option = argv[i];
			if (option.compare("-idbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				idbdir = argv[i + 1];
			} else if (option.compare("-irdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				irdir = argv[i + 1];
			} else if (option.compare("-p") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				p = atoi(argv[i + 1]);
			} else if (option.compare("-k") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				k = atoi(argv[i + 1]);
			} else if (option.compare("-perc") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				perc = atof(argv[i + 1]);
				if (perc <= 0) {
					std::cout << "\n[ERR] perc <= 0\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-silent") == 0) {
				silent = true;
				-- i;
			} else if (option.compare("-h") == 0) {
				printUsage(argv[0]);
			} else {
				std::cout << "unrecognized option " << option << "\n";
				printUsage(argv[0]);
			}
		}

		if (!silent) { printSpec(argv[0]); }
	}

private:
	int argnum;
	char** arg;

	void init () {
		silent = false;
		p = 12;
		perc = 1;
		k = 32;
	}

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Parameters\n";
			std::cout << "-p: default 12; number of cores to use\n";
			std::cout << "-k: default 32; kmer size (0, 32]\n";
			std::cout << "-perc: default 1; percentage of sampling\n";
			std::cout << "-idbdir: input folder containing database files in .fa format\n";
			std::cout << "-irdir: input folder containing read data files in .fq format\n";
			std::cout << "-silent: default false; no screen print-out\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void printSpec(char* exe) {
		if (!silent) {
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Running command: \n" << exe;
			std::cout << " -p " << p;
			std::cout << " -k " << k;
			std::cout << " -perc " << perc;
			if (!idbdir.empty()) std::cout << " -idbdir " << idbdir;
			if (!irdir.empty()) std::cout << " -irdbdir " << irdir;
			std::cout << "\n--------------------------------------------------------\n\n";
		}
	} // printSpec
};

#endif /* PARAMETER_H_ */
