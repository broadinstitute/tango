//========================================================================
// Project     : RankAbundance
// Name        : debugcode.h
// Author      : Xiao Yang
// Created on  : Dec 11, 2013
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


#ifndef DEBUGCODE_H_
#define DEBUGCODE_H_

#include "xutil.h"

void output_histogram(const ivec_t& cnts);

/* @brief	Given <kmer, cnt, flag> sorted by kmer, generate kmer cnt
 * histogram
 */
template<typename iterType>
void debug_gen_kmercnt_histogram (iterType begin, iterType end) {
	ivec_t cnts;
	for (iterType it = begin; it != end; ++ it) {
		cnts.push_back(std::get<1>(*it));
	}
	std::sort(cnts.begin(), cnts.end());
	std::cout << "\n[ --- DEBUG \n\n";
	std::cout << "Histogram (col1) freq (col2) this many kmers have\n";
	std::cout << "Freq\t\t#kmers\n-----------------------------------\n";
	output_histogram(cnts);
	std::cout << "----------------------------------------------------\n";
	std::cout << "\n DEBUG ---] \n\n";
}

/* @brief	Given record <kmer, taxaID> sorted by kmer, generate kmer
 * histogram in terms of number of taxa any kmer belongs to
 */
template<typename iterType>
void debug_gen_kmertaxa_histogram (iterType begin, iterType end) {
	ivec_t cnts;
	uint64_t kmer;
	int cnt = 1;
	if (end != begin) kmer = begin->first;
	for (iterType it = ++ begin; it != end; ++ it) {
		if (it->first == kmer) ++ cnt;
		else {
			cnts.push_back(cnt);
			cnt = 1;
			kmer = it->first;
		}
	}
	cnts.push_back(cnt);
	std::sort (cnts.begin(), cnts.end());

	std::cout << "\n[ --- DEBUG \n\n";
	std::cout << "Histogram (col1) num_taxa (col2) this many kmers belong to\n";
	std::cout << "#taxa\t\t#kmers\n-----------------------------------\n";
	output_histogram(cnts);
	std::cout << "----------------------------------------------------\n";
	std::cout << "\n DEBUG ---] \n\n";
} //


void debug_print_record_info (const bvec_t& is_record_found, const bvec_t& is_record_m,
		const std::vector<u64ipair_t>& record) ;

void debug_print_contig_skeleton (const bvec_t& is_record_m,
		const std::vector<u64ipair_t>& record, const u64vec_t& kmers) ;

void debug_print_node_kmership (const std::vector<u64ipair_t>& record,
		const bvec_t& is_k_found, const iset_t& debug_nodes, const ivec_t& nodes) ;
#endif /* DEBUGCODE_H_ */
