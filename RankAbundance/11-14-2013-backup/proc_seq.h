//========================================================================
// Project     : RankAbundance
// Name        : proc_seq.h
// Author      : Xiao Yang
// Created on  : Nov 2, 2013
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


#ifndef PROC_SEQ_H_
#define PROC_SEQ_H_

#include "xutil.h"
#include "xny/file_manip.hpp"
#include "xny/seq_manip.hpp"
#include "jaz/fastx_iterator.hpp"
#include "ReadBioFile.h"


void proc_seq (std::vector<abund_t>& abund, const std::vector<u64ipair_t>& record,
	const std::string& idir,	int k, int batch, bool silent);

void batch_proc_seq (std::vector<abund_t>& abund, bvec_t& is_found, bvec_t& is_k_found,
	const std::vector<u64ipair_t>& record, const strvec_t& reads, int k) ;

/* sort <uint64_t, int> pair based on the first field */
struct cmp_u64ipair{
public:
	bool operator () (const u64ipair_t& lhs, const u64ipair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

#endif /* PROC_SEQ_H_ */
