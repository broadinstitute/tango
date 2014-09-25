//========================================================================
// Project     : PrepDB
// Name        : proc_db.h
// Author      : Xiao Yang
// Created on  : Oct 28, 2013
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


#ifndef PROC_DB_H_
#define PROC_DB_H_

#include "xutil.h"
#include "xny/seq_manip.hpp"
#include "xny/file_manip.hpp"
#include <sys/resource.h>



void proc_db (std::vector<u64ipair_t>& record, std::vector<strvec_t>& input_files,
		strvec_t& nodeNames, const double& perc, int k,
		const std::string& idir, bool silent) ;

void sample_kmers (u64vec_t& kmers, int k, const double& perc,
		const std::string& file);

void proc_contig (u64vec_t& kmers, int k, const double& perc,
		const std::string& contig);


#endif /* PROC_DB_H_ */
