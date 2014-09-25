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
	//std::string idbdir_; // input DB dir
	std::string odir_;   // output dir
	std::string oprefix_; // output file prefix -- associated with time stamp
	std::string inodename_; // input file specifies name of node
	std::string itree_;	//input taxonomy tree
	std::string iskeleton_; // input file recording skeleton filelist
	//strvec_t idbfl_;  // input DB filelist
	strvec_t ipfqs_; // input paird read files
	strvec_t isfqs_; // input singleton read files
	int p_;
	int k_;
	//double perc_;   // downsampling percentage
	int w_;		    // window size when generating closely located kmers on ref seq
	int mcov_;
	// determine valid blocks for each skeleton
	int min_blksz_; 	// min valid block size
	int max_blkdist_; // max distance within kmers of a block
	// determine spurious taxa
	double min_perc_span_;
	int	max_perc_rep_;

 	bool silent_;

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Parameters\n";
			std::cout << "-p: default 16; number of cores to use\n";
			std::cout << "-k: default 32; kmer size (0, 32]\n";
			std::cout << "-w: default 2; window sz when generating local kmer pairs on ref seq\n";
			std::cout << "-mcov: default 3; min coverage when generating kmer pairs in read data\n";
			std::cout << "-min_blksz: default 2; minimum sz of a valid skeleton block\n";
			std::cout << "-max_blkdist: default 10; max distance between adjacent kmers of a valid block\n";
			std::cout << "-min_perc_span: default 5.0; min percentage of span of a valid taxa\n";
			std::cout << "-max_perc_rep: default 95; max percentage of repetitiveness of a valid taxa\n";
			std::cout << "-ipfqs: comma separated paired end FQ read files\n";
			std::cout << "-isfqs: comma separated singleton FQ read files\n";
			std::cout << "-iskeleton: file listing paths of skeleton files and # of sampling\n";
			std::cout << "-itree: nodes.dmp file from NCBI specifying taxonomy tree\n";
			std::cout << "-inodename: names.clean.dmp file from NCBI specifying names of taxa node\n";
			std::cout << "-odir: directory for output files\n";
			std::cout << "-oprefix: output file prefix\n";
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
			//std::cout << " -perc " << perc_;
			std::cout << " -w " << w_;
			std::cout << " -mcov " << mcov_;

			std::cout << " -min_blksz " << min_blksz_;
			std::cout << " -max_blkdist " << max_blkdist_;

			std::cout << " -min_perc_span " << min_perc_span_;
			std::cout << " -max_perc_rep " << max_perc_rep_;


			//if (!idbdir_.empty()) std::cout << " -idbdir " << idbdir_;
			/*if (!idbfl_.empty()) {
				std::cout << " -idbfl ";
				int sz = idbfl_.size();
				for (int i = 0; i < sz; ++ i) {
					if (i != sz - 1) std::cout << idbfl_[i] << "," ;
					else std::cout << idbfl_[i];
				}
			}*/
			if (!ipfqs_.empty()) {
				std::cout << " -ipfqs ";
				int sz = ipfqs_.size();
				for (int i = 0; i < sz; ++ i) {
					if (i != sz - 1) std::cout << ipfqs_[i] << "," ;
					else std::cout << ipfqs_[i];
				}
			}
			if (!isfqs_.empty()) {
				std::cout << " -isfqs ";
				int sz = isfqs_.size();
				for (int i = 0; i < sz; ++ i) {
					if (i != sz - 1) std::cout << isfqs_[i] << "," ;
					else std::cout << isfqs_[i];
				}
			}
			if (!itree_.empty()) std::cout << " -itree " << itree_;
			if (!inodename_.empty()) std::cout << " -inodename " << inodename_;
			if (!iskeleton_.empty()) std::cout << " -iskeleton " << iskeleton_;

			if (!odir_.empty()) std::cout << " -odir " << odir_;
			if (!oprefix_.empty()) std::cout << " -oprefix " << oprefix_;
			std::cout << "\n\n--------------------------------------------------------\n\n";
		}
	} // printSpec

	// set up output file prefix which is associated with time
	void set_oprefix_() {
		 mytime_t mytime = get_time ();
		 oprefix_ = std::to_string (mytime.day);
		 oprefix_ += "_" + std::to_string (mytime.hour);
		 oprefix_ += "_" + std::to_string (mytime.min);
		 oprefix_ += "_" + std::to_string (mytime.sec);
	}

public:

	Parameter (): silent_(false), p_(16),// perc_(1),
			k_(32), w_(2), mcov_(3), min_blksz_(2), max_blkdist_(10),
			min_perc_span_(5.0), max_perc_rep_(95){}

	void procinput (int argc, char** argv) {
		for (int i = 1; i < argc; i += 2) {
			std::string option = argv[i];
			/*if (option.compare("-idbdir") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				idbdir_ = argv[i + 1];
			} else */
			if (option.compare("-odir") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				odir_ = argv[i + 1];
			} else if (option.compare("-oprefix") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				oprefix_ = argv[i + 1];
			} else if (option.compare("-itree") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				itree_ = argv[i + 1];
			} else if (option.compare("-iskeleton") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				iskeleton_ = argv[i + 1];
			} else if (option.compare("-inodename") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				inodename_ = argv[i + 1];
			} /*else if (option.compare("-idbfl") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string entry = argv[i + 1];
				split (',', entry, std::back_inserter(idbfl_));
			}*/ else if (option.compare("-ipfqs") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string entry = argv[i + 1];
				split (',', entry, std::back_inserter(ipfqs_));
				if (ipfqs_.size() % 2 != 0) {
					std::cout << "\tspecify even number of files with -ipfqs\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-isfqs") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string entry = argv[i + 1];
				split (',', entry, std::back_inserter(isfqs_));
			} else if (option.compare("-p") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				p_ = atoi(argv[i + 1]);
			} else if (option.compare("-k") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				k_ = atoi(argv[i + 1]);
			} else if (option.compare("-w") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				w_ = atoi(argv[i + 1]);
			} else if (option.compare("-mcov") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				mcov_ = atoi(argv[i + 1]);
			} else if (option.compare("-min_blksz") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				min_blksz_ = atoi(argv[i + 1]);
			} else if (option.compare("-max_blkdist") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				max_blkdist_ = atoi(argv[i + 1]);
			} else if (option.compare("-min_perc_span") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				min_perc_span_ = atof(argv[i + 1]);
			} else if (option.compare("-max_perc_rep") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				max_perc_rep_ = atoi(argv[i + 1]);
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

		if (oprefix_.empty()) set_oprefix_();

		if (ipfqs_.empty() && isfqs_.empty()) {
			std::cout << "at least one of -ipfqs and -isfqs should be specified\n";
			printUsage(argv[0]);
		}

		if (iskeleton_.empty()) {
			std::cout << "-iskeleton should be specified\n";
			printUsage(argv[0]);
		} /*else if (!boost::filesystem::exists(iskeleton_)) {
			std::cout << "[ERROR]" << iskeleton_ << " doesn't exist!\n";
			exit(1);
		} */ else {
			if (inodename_.empty()) {
				std::cout << "-inodename should be specified\n";
				printUsage(argv[0]);
			}
			if (itree_.empty()) {
				std::cout << "-itree should be specified\n";
				printUsage(argv[0]);
			}
		}

		if (odir_.empty()) {
			std::cout << "-odir should be specified\n";
			printUsage(argv[0]);
		}
	}

	//std::string get_db_path() { return idbdir_; }
	std::string get_odir() { return odir_; }
	std::string get_oprefix() { return oprefix_; }
	std::string get_inodename() { return inodename_; }
	std::string get_itreefile() { return itree_; }
	std::string get_iskeletonfile () { return iskeleton_; }
	//strvec_t get_dbfl () { return idbfl_; }
	strvec_t get_ipfqs () { return ipfqs_; }
	strvec_t get_isfqs () { return isfqs_; }
	int get_num_threads () { return p_; }
	int get_k () { return k_; }
	//double get_perc () { return perc_; }
	int get_w() { return w_; }
	int get_mcov() { return mcov_; }
	int get_min_blksz() { return min_blksz_; }
	int get_max_blkdist() { return max_blkdist_; }
	double get_min_perc_span() { return min_perc_span_; }
	int get_max_perc_rep() { return max_perc_rep_; }
	bool quiet () { if (silent_) return true; else return false; }
	//bool is_idbfl_empty () { if (idbfl_.empty()) return true; else return false; }
	bool is_ipfqs_empty () { if (ipfqs_.empty()) return true; else return false; }
	bool is_isfqs_empty () { if (isfqs_.empty()) return true; else return false; }

};

#endif /* PARAMETER_H_ */
