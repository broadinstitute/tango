//========================================================================
// Project     : RankAbundance
// Name        : RankAbundance.cpp
// Author      : Xiao Yang
// Created on  : Nov 1, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2013)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or 
//				 guaranteed support whatsoever. The Broad Institute cannot 
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================


#include "Parameter.h"
#include "proc_db.h"
#include "proc_seq.h"

int main (int argc, char** argv){

	/*
	u64vec_t test (400000000);
	for (int i = 0; i < test.size(); ++ i) test[i] = i;
	for (int j = 0; j < 1; ++ j) {
		for (int i = 0; i < test.size(); ++ i) {
			if (i % 1000000 == 0) std::cout << i << "\n";
			std::binary_search(test.begin(), test.end(), i);
		}
	}
	std::cout << "done\n"; exit(1); */

	Parameter myPara (argc, argv);

 	omp_set_num_threads(myPara.p);

 	std::vector<u64ipair_t> record; // <kmer, fileID>
	//u64vec_t record; // kmers for all
//	std::vector<ipair_t> range; // <range of index wrt record, fileID>
	std::vector<strvec_t> input_files;
	strvec_t nodeNames;
	if (!myPara.idbdir.empty()) {
		if (!myPara.silent) std::cout << "Read " << myPara.idbdir << "\n";
	//	proc_db (record, range, input_files, myPara.perc, myPara.k,
	//				myPara.idbdir, myPara.silent);
		proc_db (record, input_files, nodeNames, myPara.perc, myPara.k,
					myPara.idbdir, myPara.silent);
	}

	std::vector<abund_t> abund (nodeNames.size()); // fileID [cnt]
	std::vector<ipair_t> abund_fID;

	if (!myPara.irdir.empty()) {
		int batch = 1000000;

		proc_seq (abund, record, myPara.irdir, myPara.k,
				batch, myPara.silent);
		abund_fID.resize(abund.size());
		for (int i = 0; i < abund.size(); ++ i) {
			abund_fID[i] = ipair_t (abund[i].total, i);
		}
		std::sort (abund_fID.begin(), abund_fID.end());
	}

	for (auto& x: abund_fID) {
		int fID = x.second;
		std::cout << std::setw(10) << nodeNames[fID]
		          << std::setw(10) << abund[fID].clade_k_num
		          << std::setw(10) << abund[fID].read_k_num
		          << std::setw(10) << abund[fID].total;
		if (abund[fID].read_k_num == 0) std::cout << std::setw(10) << "-" << "\n";
		else std::cout << std::setw(10) << abund[fID].total/abund[fID].read_k_num << "\n";
	}

	return (EXIT_SUCCESS);
}




