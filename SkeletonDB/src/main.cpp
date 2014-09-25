//========================================================================
// Project     : SkeletonDB
// Name        : main.cpp
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
// Prerequisite:
// 1) install boost system, filesystem library
// 2) c++0x (gcc4.7.2+)
//========================================================================


#include "Parameter.h"
#include "proc_db.h"
#include "boost/filesystem.hpp"

int main (int argc, char** argv) {

	// test get_bitkmer_by_alphabet
	/*u64vec_t klist; // all kmers from contig
	bdeque_t is_valid;
	std::string contig = "TCTGTTGGGAGACGGTGGTTTCCCGCGCATTACTACAAGCGTTAGCGTCGAGTGGTAAGAGCGTTGCAGGGTATAAACCTGTCGCTAAAGGGAGCAAAGAGACTGCAGAAGGGATGCGCAACAAAGATGCGCTGGTATTGCAAAGCGTCTCTTCGCTGGAGTTGCCGTATGAGGCGATTAATCCGATCGCGTTAAGCGAAGAAGAAAGCAGCGTCGCGCATAGCTGTCCTATAAATTACACCTTGCTGTCCAACGGTCTGGCGAGCCTGAGCGATAAAGTGGACCACGTGGTTGTGGAGGGTACTGGCGGCTGGCGTAGTTTGATGAATGATTTACGTCCATTGTCTGAATGGGTGGTACAGGAACAGTTACCGGTATTGATGGTGGTAGGGATTCAGGAGGGGTGCATTAATCATGCGCTACTGACTGCGCAAGCTGTCGCTAACGATGGACTGCCATTAATTGGCTGGGTGGCGAATCGTATCAATCCTGGCCTGGCGCATTATGCTGAAATTATTGATGTGCTTGGGAAAAAACTGCCAGCGCCGCTGATTGGCGAGTTACCGTACCTGCCGCGCGCCGAGCAGCGTGAGCTGGGGCAGTATATTCGTTTATCAATGCTCGGCAGCGTGCTGGCGGTAGATAGAATCATGGCGTAACGTTCGCGAGAGCACTGACGCGACTACACAGGCAATCAGCAGACCGGGGAGTAACTGATACTCGCCGGTCATTTCACAAATCATCAGGGTGGACATTATCGGCGCATGGGTTGTCGCCGCCAGTAGCGTCGCCATACCCGCCAGCCCCAATAAAATCGCTATTTCATCCGAGCCTGGCAGCCAGAATCCCCATATCCGACCCAGAAACATTCCAATGGATAATCCGACAAATAATGTTGGCGTAAAGACGCCGCCCGGCGCGCCAGACCCGCTGCTGGCGAGCACCGCCAGGATTTTGCACACAAAAWTCCGCCAATCAGCGAG";
	xny::get_bitkmer_by_alphabet<std::back_insert_iterator<u64vec_t>, uint64_t>
	(contig, std::back_inserter(klist), is_valid, 32);
	int index = 0;
	std::cout << contig << "\n";
	std::cout << "is valid size = " << is_valid.size() << "\n";
	for (auto& x: klist) {
		if (is_valid[index]) std::cout << index << ":" << xny::ID2Str(x, 32) << "\n";
		++index;
	}
	exit(1); */
	Parameter myPara;
	myPara.procinput(argc, argv);

 	omp_set_num_threads(myPara.get_num_threads());
 	omp_set_nested(1);

 	int batch = 1000000;
 	ProcDB myDB (myPara.get_perc(), myPara.get_k(), myPara.get_wd(),
 			myPara.get_wl(), myPara.get_minl(), batch, myPara.quiet());

 	myDB.set_db_filelist(myPara.get_db_path(), myPara.get_dbfl(),
 			myPara.get_odbdir());

 	for (auto& x: myPara.get_dbfl()) std::cout << x <<"\n";

 	ivec_t taxaKCnt;
 	myDB.set_db_skeleton_filelist(myPara.get_osfdir());

 	myDB.generate_skeleton_per_taxa(taxaKCnt);

 	myDB.output_db_skeleton_fl(myPara.get_osf(), taxaKCnt);

 	return EXIT_SUCCESS;
}



