//========================================================================
// Project     : RankAbundance
// Name        : proc_seq.cpp
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


#include "proc_seq.h"


void proc_seq (std::vector<abund_t>& abund, const std::vector<u64ipair_t>& record,
	const std::string& idir,	int k, int batch, bool silent){

	int num_record = record.size();
	strvec_t files;
	strvec_t fq_files;
	if (!xny::dir_list(idir, std::back_inserter(files))) {
		abording ("proc_seq dir_list failed");
	}
	for (auto& f: files) {
		std::string suf = xny::get_suffix (f, true, ".");
		if (suf.compare("fq") == 0 || suf.compare("fastq") == 0) {
			std::string path = idir + "/";
			fq_files.push_back(path + f);
		}
	}

	bvec_t is_k_found (num_record, false); //

	for (auto& f: fq_files) {
		if (!silent) std::cout << "\tproc: " << f << "\n";

		std::ifstream fh;
		xny::openfile<std::ifstream>(fh, f);
		bio::fastq_input_iterator<> fq(fh), end;

		int total_reads = 0;
		strvec_t reads;
		int unfound = 0;
		while (fq != end) {

	 		add_fq_reads_only (reads, batch, fq, end);

	 		//std::cout << batch << "\n";

	 		// keep track of whether any kmer in a read is found
	 		bvec_t is_r_found (reads.size(), false);
			batch_proc_seq (abund, is_r_found, is_k_found, record, reads, k);
			for (int i = 0; i < is_r_found.size(); ++ i) if (!is_r_found[i]) ++ unfound;

			total_reads += reads.size();

			reads.clear();

		} // while

		if (!silent) {
			std::cout << "\t\ttotal reads: " << total_reads <<  ", "
					<< unfound << " not found\n";
		}

		xny::closefile(fh);
	}

	// update abund
	//ivec_t baseline (abund.size(), 0);
	for (int i = 0; i < num_record; ++ i) {
		int fID = record[i].second;
		if (fID >= abund.size()) abording ("proc_seq SC failed");

		if (is_k_found[i]) ++ abund[fID].read_k_num;
		++ abund[fID].clade_k_num;
	}

}

void batch_proc_seq (std::vector<abund_t>& abund, bvec_t& is_found, bvec_t& is_k_found,
	const std::vector<u64ipair_t>& record, const strvec_t& reads, int k) {

	int num_r = reads.size();
	int num_f = abund.size();

	//std::cout << "num_r, num_f, num_record " << num_r << ", "
	//		<< num_f << ", " << record.size() << "\n";

	#pragma omp parallel
	{
		ivec_t private_abund (num_f, 0);
		#pragma omp for
		for (int r = 0; r < num_r; ++ r) {

			//if (r % 1000 == 0) std::cout << r << "\n";

			u64vec_t kmers;
			xny::get_bitkmer<std::back_insert_iterator<u64vec_t>, uint64_t>
				(reads[r], std::back_inserter(kmers), k, 0);
			u64set_t kmerset (kmers.begin(), kmers.end());
			for (auto& x: kmerset) {
				// search record
				//for (int i = 0; i < range.size(); ++ i) {
				//	int start = 0, end = range[i].first;
				//	if (i > 0) start = range[i - 1].first;
				//	if (std::binary_search(record.begin() + start,
				//			record.begin() + end, x)) {
				auto lb = std::lower_bound(record.begin(), record.end(),
						u64ipair_t (x, 0), cmp_u64ipair());
				int dist = std::distance(record.begin(), lb);
				for (; ; ++ lb, ++ dist) {
					if (lb->first == x) {
						is_k_found[dist] = true;
						int fileID = lb->second;
						if (fileID >= num_f) abording ("batch_proc_seq SC failed");
						++ private_abund[fileID];
						is_found[r] = true;
					} else break;
				}
			}
		}
		#pragma omp critical
		{
			for (int i = 0; i < num_f; ++ i) {
				abund[i].total += private_abund[i];
			}
		}
	}

}
