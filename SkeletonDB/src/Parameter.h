//========================================================================
// Project     : SkeletonDB
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
	std::string idbdir_;  // input DB dir
	strvec_t idbfl_;  // input DB filelist
	std::string odbdir_;	 // output reduced db dir
	std::string osfdir_; // output dir for skeleton files
	std::string osf_;	 // output file recording skeleton filelist
	int p_;
	int k_;
	double perc_;   // downsampling percentage
	int wl_; 		// window length
	int wd_; 		// window density
	int minl_;      // min contig length
 	bool q_;

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Parameters\n";
			std::cout << "-p: default 16; number of cores to use\n";
			std::cout << "-k: default 32; kmer size (0, 32]\n";
			std::cout << "-perc: default 1; percentage of sampling\n";
			std::cout << "-wl: default 10, window length\n";
			std::cout << "-wd: default 80 (/100), min window density\n";
			std::cout << "-minl: default 2000bp, min contig length\n";
			std::cout << "-idbdir: input folder containing database files in .fa format\n";
			std::cout << "-idbfl: comma separated database files in .fa format\n";
			std::cout << "-odbdir:\n-osfdir: output dir recording reduced db and skeleton files"
					"\n\t[WARNING]: these folders will be overwritten\n";
			std::cout << "-osf: output file recording skeleton filelist\n";
			std::cout << "-q: default false; quiet\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void printSpec_(char* exe) {
		if (!q_) {
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Running command: \n" << exe;
			std::cout << " -p " << p_;
			std::cout << " -k " << k_;
			std::cout << " -perc " << perc_;
			std::cout << " -wl " << wl_;
			std::cout << " -wd " << wd_;
			std::cout << " -minl " << minl_;
			if (!idbdir_.empty()) std::cout << " -idbdir " << idbdir_;
			if (!idbfl_.empty()) {
				std::cout << " -idbfl ";
				int sz = idbfl_.size();
				for (int i = 0; i < sz; ++ i) {
					if (i != sz - 1) std::cout << idbfl_[i] << "," ;
					else std::cout << idbfl_[i];
				}
			}
			if (!odbdir_.empty()) std::cout << " -odbdir " << odbdir_;
			if (!osfdir_.empty()) std::cout << " -osfdir " << osfdir_;

			if (!osf_.empty()) std::cout << " -osf " << osf_;
			std::cout << "\n\n--------------------------------------------------------\n\n";
		}
	} // printSpec

	void createDIR(const std::string& dir) {
		if (! boost::filesystem::exists(dir)) {
			if (! boost::filesystem::create_directory(dir)) {
				abording ("Could not create " + dir);
			}
		} else {
			if (!boost::filesystem::remove_all(dir)) {
				 abording ("Could not remove " + dir);
			}
			if (! boost::filesystem::create_directory(dir)) {
				abording ("Could not create " + dir);
			}
		}
	}

public:

	Parameter (): p_(16), k_(32), perc_(1), q_(false),
		wl_(10), wd_(80), minl_(2000) {}

	void procinput (int argc, char** argv) {
		for (int i = 1; i < argc; i += 2) {
			std::string option = argv[i];
			if (option.compare("-idbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				idbdir_ = argv[i + 1];
			} else if (option.compare("-idbfl") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string entry = argv[i + 1];
				split (',', entry, std::back_inserter(idbfl_));
			} else if (option.compare("-odbdir") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				odbdir_ = argv[i + 1];
			} else if (option.compare("-osfdir") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				osfdir_ = argv[i + 1];
			} else if (option.compare("-osf") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				osf_ = argv[i + 1];
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
			} else if (option.compare("-wl") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				wl_ = atof(argv[i + 1]);
				if (wl_ <= 0) {
					std::cout << "\n[ERR] wl <= 0\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-wd") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				wd_ = atof(argv[i + 1]);
				if (wd_ <= 50 || wd_ > 100) {
					std::cout << "\n[ERR] wd should be [50, 100]\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-minl") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				minl_ = atoi(argv[i + 1]);
			} else if (option.compare("-q") == 0) {
				q_ = true;
				-- i;
			} else if (option.compare("-h") == 0) {
				printUsage(argv[0]);
			} else {
				std::cout << "unrecognized option " << option << "\n";
				printUsage(argv[0]);
			}
		}
		if (!q_) { printSpec_(argv[0]); }

		if (idbdir_.empty() && idbfl_.empty()) {
			std::cout << "at least one of -idbdir and -idbfl should be specified\n";
			printUsage(argv[0]);
		}

		if (!odbdir_.empty()) createDIR(odbdir_);

		if (osfdir_.empty()) {
			std::cout << "-osfdir_ should be specified\n";
			printUsage(argv[0]);
		} else createDIR(osfdir_);

		if (osf_.empty()) {
			std::cout << "-osf should be specified\n";
			printUsage(argv[0]);
		}
	}

	std::string get_db_path() { return idbdir_; }
	//std::string get_odir() { return odir_; }
	std::string get_osf() { return osf_; }
	std::string get_odbdir() { return odbdir_; }
	std::string get_osfdir() { return osfdir_; }
	strvec_t get_dbfl () { return idbfl_; }
	int get_num_threads () { return p_; }
	int get_k() { return k_; }
	int get_wl() { return wl_; }
	int get_wd() { return wd_; }
	int get_minl() { return minl_; }
	double get_perc () { return perc_; }
	bool quiet () { if (q_) return true; else return false; }
	bool is_idbfl_empty () { if (idbfl_.empty()) return true; else return false; }
};

#endif /* PARAMETER_H_ */
