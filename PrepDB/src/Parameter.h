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
	std::string iginoderaw; // raw gi->nodeID dmp file
	std::string iginode; // input gi->nodeID file
	std::string oginode; // output gi->nodeID file
	std::string itree;	// input file specifying tree structure
	//std::string inodename; // input file specifying node annotations
	std::string idbdir; // input DB dir
	std::string odbdir; // output DB dir
	std::string level;		// level taxonomic level
	std::string ogilevel;	// gi--> level map file
	int p;
	bool silent;
	Parameter (int argc, char** argv): argnum(argc), arg(argv){
		init ();
		for (int i = 1; i < argnum; i += 2) {
			std::string option = argv[i];
			if (option.compare("-iginode") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				iginode = argv[i + 1];
			} else if (option.compare("-iginoderaw") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				iginoderaw = argv[i + 1];
			} else if (option.compare("-oginode") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				oginode = argv[i + 1];
			} else if (option.compare("-idbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				idbdir = argv[i + 1];
			} else if (option.compare("-odbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				odbdir = argv[i + 1];
			} else if (option.compare("-itree") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				itree = argv[i + 1];
			} //else if (option.compare("-inodename") == 0) {
			//	if (argc < i + 2) printUsage (argv[0]);
			//	inodename = argv[i + 1];
			//}
			else if (option.compare("-level") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				level = argv[i + 1];
			} else if (option.compare("-ogilevel") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				ogilevel = argv[i + 1];
			} else if (option.compare("-p") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				p = atoi(argv[i + 1]);
			} else if (option.compare("-silent") == 0) {
				silent = true;
				-- i;
			} else if (option.compare("-h") == 0) {
				printUsage(argv[0]);
			} else {
				std::cout << "unrecognized option: " << option << "\n";
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
		level = "species";
		p = 16;
	}

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Parameters\n";
			std::cout << "-iginoderaw: (opt) raw gi_node dmp file from ncbi\n";
			std::cout << "-iginode: (opt) input gi_node file\n";
			std::cout << "-oginode: (opt) output giRange_node file\n";
			std::cout << "-itree: (opt) input file specifies tree structure\n";
			std::cout << "-idbdir: (opt) input folder contains databases; all files in fa format are considered\n";
			std::cout << "-odbdir: (opt) output folder to store output db files\n";
			std::cout << "-p: default 16; # procs\n";
			std::cout << "-level: default species; lowest tax level to be considered\n";
			std::cout << "-ogilevel: (opt) output file recording [gi][nodeID][level]\n";
			std::cout << "-silent: default false; no screen print-out\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void printSpec(char* exe) {
		if (!silent) {
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Running command: \n" << exe;
			if (!iginoderaw.empty()) std::cout << " -iginoderaw " << iginoderaw;
			if (!iginode.empty()) std::cout << " -iginode " << iginode;
			if (!oginode.empty()) std::cout << " -oginode " << oginode;
			if (!ogilevel.empty()) std::cout << " -ogilevel " << ogilevel;
			if (!itree.empty()) std::cout << " -itree " << itree;
			//if (!inodename.empty())	std::cout << " -inodename " << inodename;
			if (!idbdir.empty()) std::cout << " -idbdir " << idbdir;
			if (!odbdir.empty()) std::cout << " -odbdir " << odbdir;
			std::cout << " -p " << p;
			std::cout << " -level " << level << "\n";
			std::cout << "\n--------------------------------------------------------\n\n";
		}
	} // printSpec
};

#endif /* PARAMETER_H_ */
