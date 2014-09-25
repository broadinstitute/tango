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
#include "boost/filesystem.hpp"

/* input parameters */
class Parameter{

private:
	std::string idbdir_; // input DB dir
	std::string odir_;   // output dir storing skeleton files
	std::string osf_; // output file recording skeleton filelist
	strvec_t idbfl_;  // input DB filelist
	int p_;
	int k_;
	double perc_;   // downsampling percentage
 	bool silent_;

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Parameters\n";
			std::cout << "-p: default 12; number of cores to use\n";
			std::cout << "-k: default 32; kmer size (0, 32]\n";
			std::cout << "-perc: default 1; percentage of sampling\n";
			std::cout << "-idbdir: input folder containing database files in .fa format\n";
			std::cout << "-idbfl: comma separated database files in .fa format\n";
			std::cout << "-osf: file recording skeleton file list & # of sampled kmers\n";
			std::cout << "-odir: directory for output files\n";
			std::cout << "-silent: default false; no screen print-out\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void printSpec_(char* exe) {
		if (!silent_) {
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Running command: \n" << exe;
			std::cout << " -p " << p_;
			std::cout << " -k " << k_;
			std::cout << " -perc " << perc_;
			if (!idbdir_.empty()) std::cout << " -idbdir " << idbdir_;
			if (!idbfl_.empty()) {
				std::cout << " -idbfl ";
				int sz = idbfl_.size();
				for (int i = 0; i < sz; ++ i) {
					if (i != sz - 1) std::cout << idbfl_[i] << "," ;
					else std::cout << idbfl_[i];
				}
			}
			if (!osf_.empty()) std::cout << " -osf " << osf_;

			if (!odir_.empty()) std::cout << " -odir " << odir_;
			std::cout << "\n\n--------------------------------------------------------\n\n";
		}
	} // printSpec

public:

	Parameter (): silent_(false), p_(12), perc_(1), k_(32) {}

	void procinput (int argc, char** argv) {
		for (int i = 1; i < argc; i += 2) {
			std::string option = argv[i];
			if (option.compare("-idbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				idbdir_ = argv[i + 1];
			} else if (option.compare("-odir") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				odir_ = argv[i + 1];
			} else if (option.compare("-osf") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				osf_ = argv[i + 1];
			} else if (option.compare("-idbfl") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string entry = argv[i + 1];
				split (',', entry, std::back_inserter(idbfl_));
			} else if (option.compare("-p") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				p_ = atoi(argv[i + 1]);
			} else if (option.compare("-k") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				k_ = atoi(argv[i + 1]);
			} else if (option.compare("-perc") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				perc_ = atof(argv[i + 1]);
				if (perc_ <= 0) {
					std::cout << "\n[ERR] perc <= 0\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-silent") == 0) {
				silent_ = true;
				-- i;
			} else if (option.compare("-h") == 0) {
				printUsage(argv[0]);
			} else {
				std::cout << "unrecognized option " << option << "\n";
				printUsage(argv[0]);
			}
		}
		if (!silent_) { printSpec_(argv[0]); }

		if (idbdir_.empty() && idbfl_.empty()) {
			std::cout << "at least one of -idbdir and -idbfl_ should be specified\n";
			printUsage(argv[0]);
		}

		if (osf_.empty()) {
			std::cout << "-osf should be specified\n";
			printUsage(argv[0]);
		}

		if (odir_.empty()) {
			std::cout << "-odir should be specified\n";
			printUsage(argv[0]);
		}
	}

	std::string get_db_path() { return idbdir_; }
	std::string get_odir() { return odir_; }
	std::string get_osf () { return osf_; }
	strvec_t get_dbfl () { return idbfl_; }
	int get_num_threads () { return p_; }
	int get_k () { return k_; }
	double get_perc () { return perc_; }
	bool is_no_prt () { if (silent_) return true; else return false; }
	bool is_idbfl_empty () { if (idbfl_.empty()) return true; else return false; }
};

#endif /* PARAMETER_H_ */
