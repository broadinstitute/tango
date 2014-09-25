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
#include "debugcode.h"

void read_context (std::vector<u64ipair_t>& record, bvec_t& is_record_m,
		const strvec_t& ipfqs, const strvec_t& isfqs, int k, int batch,
		bool silent);


void batch_read_context (bvec_t& is_record_found, bvec_t& is_record_m,
		std::vector<u64ipair_t>& record, const strvec_t& frags, int k);

void clean_records (std::vector<u64ipair_t>& record, bvec_t& is_record_m,
		const bvec_t& is_record_found);

void read_count (std::vector<clade_t>& cladeinfo,
		const std::vector<u64ipair_t>& record, const bvec_t& is_record_m,
		const strvec_t& ipfqs, const strvec_t& isfqs,
		int k, int batch, bool silent);

void batch_read_count (std::vector<clade_t>& cladeinfo, bvec_t& is_frag_found,
		const bvec_t& is_record_m, const std::vector<u64ipair_t>& record,
		const strvec_t& frags, int k);

void proc_seq (std::vector<abund_t>& abund, 	bvec_t& is_k_found,
		const std::vector<u64ipair_t>& record,
	const strvec_t& ipfqs, const strvec_t& isfqs, int k, int batch, bool silent);

void batch_proc_seq (std::vector<abund_t>& abund, bvec_t& is_r_found,
	std::vector<iset_t>& r_nodes, bvec_t& is_k_found,
	const std::vector<u64ipair_t>& record, const strvec_t& reads, int k);



#endif /* PROC_SEQ_H_ */
