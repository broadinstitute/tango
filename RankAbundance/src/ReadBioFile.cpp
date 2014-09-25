//========================================================================
// Project     : M-Vicuna
// Name        : ReadBioFile.cpp
// Author      : Xiao Yang
// Created on  : Jun 4, 2013
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

#include "ReadBioFile.h"


/**	Retrieve num records (header, read, qual) from fastq file pointed by
 * iter
 */
void add_fq_reads (std::vector<fqtuple_t>& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (*iter);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} // add_fq_reads

void add_fq_reads_only (strvec_t& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (std::get<1>(*iter));
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
}

/* concatenate paired-end reads to be a single fragment separated by an 'n'
 */
void create_frags (strvec_t& frags) {
	int num_frags = frags.size();
	if (num_frags % 2 != 0) 	abording ("proc_seq SC failed: "
			"number of paired reads not an even number");
	num_frags /= 2;
	for (int fragidx = 0; fragidx < num_frags; ++ fragidx) {
		frags[fragidx] += 'N';
		frags[fragidx] += frags[fragidx + num_frags];
		frags[fragidx + num_frags].clear();
	}
	frags.resize(num_frags);
} // create_fragments

/* get unique kmer set from input fragment
 */
void get_unique_frag_kmers (u64set_t& frag_kval_set, int k, const std::string& frag) {
	int k_pos_idx = 0;
	int last_k_pos = frag.length() - k;
	while (k_pos_idx <= last_k_pos) {
		//look for current kmer
		uint64_t val;
		std::string kmer = get_kstr(frag.substr(k_pos_idx, k));

		if (xny::str2ID<uint64_t>(val, kmer)) {
			frag_kval_set.insert(val);
			++ k_pos_idx;
		} else { // skip next found 'n/N'
			// find the last position of N
			int pos_n = kmer.find_last_of("nN");
			if (pos_n != std::string::npos) {
				k_pos_idx += pos_n + 1;
			} else k_pos_idx += k;
		}
	} // while (k_pos_idx <= last_k_pos)
}

void add_fa_reads (std::vector<strpair_t>& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (*iter);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} //add_fa_reads

void add_fa_reads_only (strvec_t& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (iter->second);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} //add_fa_reads

