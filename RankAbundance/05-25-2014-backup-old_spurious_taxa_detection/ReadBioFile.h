//========================================================================
// Project     : M-Vicuna
// Name        : ReadBioFile.h
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


#ifndef READBIOFILE_H_
#define READBIOFILE_H_

#include "xutil.h"
#include "jaz/fastx_iterator.hpp"
#include "xny/file_manip.hpp"

void add_fq_reads (std::vector<fqtuple_t>& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end);

void add_fq_reads_only (strvec_t& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end);

void create_frags (strvec_t& frags);
void get_unique_frag_kmers (u64set_t& frag_kval_set, int k, const std::string& frag);

void add_fa_reads (std::vector<strpair_t>& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end);

void add_fa_reads_only (strvec_t& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end);

// F: functor   D: data to be processed
template <typename F, typename D>
void process_fq_fl(F func, D& data, int batch, const strvec_t& fl,
		bool is_pe = true, bool silent = false) {

	for (int fidx = 0; fidx < fl.size(); ++ fidx) {

		if (!silent) {
			if (is_pe) {
				++ fidx;
				std::cout << "\tproc paired-end files: \n\t\t"
					<< fl[fidx - 1] << "\n\t\t" << fl[fidx] << "\n";
			} else std::cout << "\tproc single-end file: " << fl[fidx] << "\n";
		}

		std::ifstream fh0, fh1;
		if (is_pe) {
			xny::openfile<std::ifstream>(fh0, fl[fidx-1]);
			xny::openfile<std::ifstream>(fh1, fl[fidx]);
		} else xny::openfile<std::ifstream>(fh0, fl[fidx]);

		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch/2, fq0, end);
	 		if (is_pe) {
	 			while (fq1 != end) {
	 				add_fq_reads_only (frags, batch/2, fq1, end);
	 		 		create_frags (frags);
	 			}
	 		}

	 		func (data, frags);

	 		total_frags += frags.size();

	 		frags.clear();

		} // while

		if (!silent) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";

		xny::closefile(fh0);
		if (is_pe) xny::closefile(fh1);
	} // for (int fidx = 0; fidx < fl.size(); ++ fidx) {
} // process_fq_fl

#endif /* READBIOFILE_H_ */
