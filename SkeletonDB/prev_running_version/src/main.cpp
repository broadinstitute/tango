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

	Parameter myPara;
	myPara.procinput(argc, argv);

 	omp_set_num_threads(myPara.get_num_threads());

 	int batch = 1000000;
 	ProcDB myDB (myPara.get_perc(), myPara.get_k(), batch, myPara.is_no_prt());

 	myDB.set_db_filelist(myPara.get_db_path(), myPara.get_dbfl());

 	ivec_t taxaKCnt;

 	for (auto& x: myPara.get_dbfl()) std::cout << x <<"\n";
 	myDB.set_db_skeleton_filelist(myPara.get_odir());
 	myDB.generate_skeleton_per_taxa(taxaKCnt);
 	myDB.output_db_skeleton_fl(myPara.get_osf(), taxaKCnt);
 	return EXIT_SUCCESS;

}



