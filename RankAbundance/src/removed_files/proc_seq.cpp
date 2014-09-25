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

/* @brief	Scan read data and update [record] by removing ones not found;
 *			assign record to be 'm'ultiple if it occurs in multiple taxa
 */
void read_context (std::vector<u64ipair_t>& record, bvec_t& is_record_m,
		const strvec_t& ipfqs, const strvec_t& isfqs, int k, int batch,
		bool silent) {

	int num_record = record.size();
	bvec_t is_record_found (num_record, false);

	if (!silent) std::cout << "\tInitial record size: " << num_record << "\n";

	//-------------------------------------------------------------
	//----------------- process paired-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {
		if (!silent) std::cout << "\tproc paired-end files: \n\t\t"
				<< ipfqs[fidx-1] << "\n\t\t" << ipfqs[fidx] << "\n";

		std::ifstream fh0, fh1;
		xny::openfile<std::ifstream>(fh0, ipfqs[fidx-1]);
		xny::openfile<std::ifstream>(fh1, ipfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end && fq1 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch/2, fq0, end);
	 		add_fq_reads_only (frags, batch/2, fq1, end);
	 		create_frags (frags);

 			batch_read_context (is_record_found, is_record_m, record, frags, k);

			total_frags += frags.size();

			frags.clear();

		} // while

		if (!silent) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";


		xny::closefile(fh0);
		xny::closefile(fh1);
	} // for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {

	//-------------------------------------------------------------
	//----------------- process single-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 0; fidx < isfqs.size(); ++ fidx) {
		if (!silent) std::cout << "\tproc single-end file: "	<< isfqs[fidx] << "\n";

		std::ifstream fh0;
		xny::openfile<std::ifstream>(fh0, isfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {

	 		add_fq_reads_only (frags, batch, fq0, end);

	 		batch_read_context (is_record_found, is_record_m, record, frags, k);

			total_frags += frags.size();

 			frags.clear();
		} // while

		if (!silent) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";
		xny::closefile(fh0);

	}

	//debug_print_record_info (is_record_found, is_record_m, record);

	clean_records (record, is_record_m, is_record_found);

	if (!silent) std::cout << "\t# record found in read data: "
			<< record.size() << "\n";

	//debug_print_record_info (is_record_found, is_record_m, record);

} // read_context


void batch_read_context (bvec_t& is_record_found, bvec_t& is_record_m,
		std::vector<u64ipair_t>& record, const strvec_t& frags, int k){

	int num_frag = frags.size();

 	#pragma omp parallel for
	for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

		int k_pos_idx = 0;
		int last_k_pos = frags[fragidx].length() - k;
		iset_t fIDs; // record fIDs of kmers present in the fragment
		ivec_t record_indices; // indices of records present in the fragment

		// --- go through each kmers in the frag ---
		while (k_pos_idx <= last_k_pos) {
			//look for current kmer
			uint64_t val;
			std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k));

			if (xny::str2ID<uint64_t>(val, kmer)) {
				// look for val in the record
				auto lb = std::lower_bound(record.begin(), record.end(),
						u64ipair_t (val, 0), cmp_u64ipair());
				int dist = std::distance(record.begin(), lb);

				for (; ; ++ lb, ++ dist) {
					if (lb->first == val) {
						fIDs.insert(lb->second);
						record_indices.push_back(dist);
						is_record_found [dist] = true;
					} else break;
				}
				++ k_pos_idx;
			} else { // skip next found 'n/N'
				// find the last position of N
				int pos_n = kmer.find_last_of("nN");
				if (pos_n != std::string::npos) {
					k_pos_idx += pos_n + 1;
				} else k_pos_idx += k;
			}
		} // while (k_pos_idx <= last_k_pos)

		// ---
		if (fIDs.size() > 1) { // frag contains kmers from multi-taxa
			for (auto& x: record_indices) {
				is_record_m[x] = true;

				/*if (record[x].first == 7109313158666873182) {
					std::cout << frags[fragidx] << "\nfound!";
					exit(1);
				}*/
			}
		}

	} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

} // batch_read_context

/* @brief	Given is_record_found vector, clean out records that have not
 * been found in the data and meanwhile update is_record_m correspondingly
 */
void clean_records (std::vector<u64ipair_t>& record, bvec_t& is_record_m,
		const bvec_t& is_record_found) {
	int num = record.size();
	int p0 = 0, p1;
	while (p0 < num) {
		if (!is_record_found[p0]) {
			p1 = p0 + 1;
			break;
		} else ++ p0;
	}
	while (p1 < num) {
		if (is_record_found[p1]) {
			record[p0] = record[p1];
			is_record_m[p0] = is_record_m[p1];
			++ p0;
			++ p1;
		} else ++ p1;
	}
	record.resize(p0);
	is_record_m.resize(p0);
} // clean_records

/* @brief	Counting clades via read data
 */
void read_count (std::vector<clade_t>& cladeinfo,
		const std::vector<u64ipair_t>& record, const bvec_t& is_record_m,
		const strvec_t& ipfqs, const strvec_t& isfqs,
		int k, int batch, bool silent) {

	int num_record = record.size();
	//-------------------------------------------------------------
	//----------------- process paired-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {
		if (!silent) std::cout << "\tproc paired-end files: \n\t\t"
				<< ipfqs[fidx-1] << "\n\t\t" << ipfqs[fidx] << "\n";

		std::ifstream fh0, fh1;
		xny::openfile<std::ifstream>(fh0, ipfqs[fidx-1]);
		xny::openfile<std::ifstream>(fh1, ipfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end && fq1 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch/2, fq0, end);
	 		add_fq_reads_only (frags, batch/2, fq1, end);
	 		create_frags (frags);

	 		bvec_t is_frag_found (frags.size(), false);
 			batch_read_count (cladeinfo, is_frag_found,
 					is_record_m, record, frags, k);

			total_frags += frags.size();
			frags.clear();
		} // while

		if (!silent) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";

		xny::closefile(fh0);
		xny::closefile(fh1);
	} // for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {

	//-------------------------------------------------------------
	//----------------- process single-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 0; fidx < isfqs.size(); ++ fidx) {
		if (!silent) std::cout << "\tproc single-end file: "	<< isfqs[fidx] << "\n";

		std::ifstream fh0;
		xny::openfile<std::ifstream>(fh0, isfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {

	 		add_fq_reads_only (frags, batch, fq0, end);

	 		bvec_t is_frag_found (frags.size(), false);
	 		batch_read_count (cladeinfo, is_frag_found,
	 				is_record_m, record, frags, k);

			total_frags += frags.size();

 			frags.clear();
		} // while

		if (!silent) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";
		xny::closefile(fh0);
	} // for

} // read_count


void batch_read_count (std::vector<clade_t>& cladeinfo, bvec_t& is_frag_found,
		const bvec_t& is_record_m, const std::vector<u64ipair_t>& record,
		const strvec_t& frags, int k) {

	int num_frag = frags.size();
	int num_node = cladeinfo.size();
	#pragma omp parallel
	{
		std::vector<ipair_t> counts (num_node, ipair_t(0, 0)); // unique, multi
		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {
			int k_pos_idx = 0;
			int last_k_pos = frags[fragidx].length() - k;
			imap_t fID_cnt; // record fIDs of kmers present in the fragment
			int total_k_instance = 0;
			// --- go through each kmers in the frag ---
			bool is_m = false;
			while (k_pos_idx <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k));

				if (xny::str2ID<uint64_t>(val, kmer)) {
					// look for val in the record
					auto lb = std::lower_bound(record.begin(), record.end(),
							u64ipair_t (val, 0), cmp_u64ipair());
					int dist = std::distance(record.begin(), lb);

					for (; ; ++ lb, ++ dist) {
						if (lb->first == val) {
							auto it_fc = fID_cnt.find(lb->second);
							if (it_fc != fID_cnt.end()) { // found
								++ it_fc->second;
							} else fID_cnt[lb->second] = 1;
							if (is_record_m[dist]) is_m = true;
							++ total_k_instance;
						} else break;
					}
					++ k_pos_idx;
				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos_idx += pos_n + 1;
					} else k_pos_idx += k;
				}
			} // while (k_pos_idx <= last_k_pos)

			// --- accumulate counts ---
			for (auto& x: fID_cnt) x.second *= 1000 / total_k_instance;
			if (fID_cnt.size() > 1) is_m = true;
			for (auto& x: fID_cnt) {
				if (is_m) counts[x.first].second += x.second;
				else counts[x.first].first += x.second;
			}
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

		#pragma omp critical
		{
			for (int i = 0; i < num_node; ++ i) {
				cladeinfo[i].u_cnt += counts[i].first;
				cladeinfo[i].m_cnt += counts[i].second;
			}
		}
	} // # pragma omp parallel
} // batch_read_count

void proc_seq (std::vector<abund_t>& abund, 	bvec_t& is_k_found,
		const std::vector<u64ipair_t>& record, const strvec_t& ipfqs,
		const strvec_t& isfqs, int k, int batch, bool silent){

	int num_record = record.size();

//
	std::vector<iset_t> all_frag_nodes;

	//-------------------------------------------------------------
	//----------------- process paired-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {
		if (!silent) std::cout << "\tproc paired-end files: \n\t\t"
				<< ipfqs[fidx-1] << "\n\t\t" << ipfqs[fidx] << "\n";

		std::ifstream fh0, fh1;
		xny::openfile<std::ifstream>(fh0, ipfqs[fidx-1]);
		xny::openfile<std::ifstream>(fh1, ipfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		int unfound = 0;
		while (fq0 != end && fq1 != end) {
	 		// cat paired-end reads to be a single fragment separated by an 'n'
	 		add_fq_reads_only (frags, batch/2, fq0, end);
	 		add_fq_reads_only (frags, batch/2, fq1, end);
	 		create_frags (frags);
	 		int num_frags = frags.size();

	 		// keep track of whether any kmer in a read is found
	 		bvec_t is_frag_found (num_frags, false);
	 		std::vector<iset_t> frag_nodes (num_frags);
			batch_proc_seq (abund, is_frag_found, frag_nodes, is_k_found,
					record, frags, k);

			all_frag_nodes.insert (all_frag_nodes.end(),
					frag_nodes.begin(), frag_nodes.end());

			for (int i = 0; i < is_frag_found.size(); ++ i)
				if (!is_frag_found[i]) ++ unfound;

			total_frags += frags.size();

			frags.clear();

		} // while

		if (!silent) {
			std::cout << "\t\ttotal fragments: " << total_frags <<  ", "
					<< unfound << " not found\n";
		}

		xny::closefile(fh0);
		xny::closefile(fh1);
	} // for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {

	//-------------------------------------------------------------
	//----------------- process single-end files ------------------
	//-------------------------------------------------------------
	for (int fidx = 0; fidx < isfqs.size(); ++ fidx) {
		if (!silent) std::cout << "\tproc single-end file: "	<< isfqs[fidx] << "\n";

		std::ifstream fh0;
		xny::openfile<std::ifstream>(fh0, isfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), end;

		int total_frags = 0;
		strvec_t frags;
		int unfound = 0;

		while (fq0 != end) {

	 		add_fq_reads_only (frags, batch, fq0, end);
	 		int num_frags = frags.size();

	 		// keep track of whether any kmer in a read is found
	 		bvec_t is_frag_found (num_frags, false);
	 		std::vector<iset_t> frag_nodes (num_frags);
			batch_proc_seq (abund, is_frag_found, frag_nodes, is_k_found,
					record, frags, k);

			all_frag_nodes.insert (all_frag_nodes.end(),
					frag_nodes.begin(), frag_nodes.end());

			for (int i = 0; i < is_frag_found.size(); ++ i)
				if (!is_frag_found[i]) ++ unfound;

			total_frags += frags.size();

			frags.clear();

		} // while

		if (!silent) {
			std::cout << "\t\ttotal fragments: " << total_frags <<  ", "
					<< unfound << " not found\n";
		}

		xny::closefile(fh0);
	}

	// update abund
	for (int i = 0; i < num_record; ++ i) {
		int fID = record[i].second;
		if (fID >= abund.size()) abording ("proc_seq SC failed");

		if (is_k_found[i]) {
			++ abund[fID].read_k_num;
			// debug here:
			{
			//	std::cout << xny::ID2Str(record[i].first, k) << "\n";
			}
		}
		++ abund[fID].clade_k_num;
	}

	{ // generate histogram of #kmers in frag
		ivec_t cnts;
		for (auto& x: all_frag_nodes) cnts.push_back(x.size());
		std::sort (cnts.begin(), cnts.end());

		std::vector<ipair_t> cnt_freq;
		int freq = 1;
		int cnt;
		if (cnts.size()) cnt = cnts[0];
		for (int i = 1; i < cnts.size(); ++ i) {
			if (cnts[i] == cnt) ++ freq;
			else {
				cnt_freq.push_back(ipair_t (cnt, freq));
				freq = 1;
				cnt = cnts[i];
			}
		}
		cnt_freq.push_back(ipair_t (cnt, freq));
		std::sort (cnt_freq.begin(), cnt_freq.end());
		std::cout << "#kmers in frags histogram\n#kfreq\t#frags\n------------\n";
		for (auto& x: cnt_freq) {
			std::cout << x.first << "\t" << x.second << "\n";
		}
	}
}

void batch_proc_seq (std::vector<abund_t>& abund, bvec_t& is_frag_found,
	std::vector<iset_t>& frag_nodes, bvec_t& is_k_found,
	const std::vector<u64ipair_t>& record, const strvec_t& frags, int k) {

	int num_frag = frags.size();
	int num_node = abund.size();

	#pragma omp parallel
	{
		std::vector<ipair_t> private_abund (num_node, ipair_t(0, 0)); // unique, multi
		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

			imap_t fileID_cnt;
			int last_k_pos = frags[fragidx].length() - k;
			int num_k_sampled = 0; // for the current fragment
			int total_fIDs = 0;
			int k_pos = 0;

			// --- go through each non-overlapping kmers in the frag ---
			while (k_pos <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos, k));

				if (xny::str2ID<uint64_t>(val, kmer)) {
					// look for this val in the record
					auto lb = std::lower_bound(record.begin(), record.end(),
							u64ipair_t (val, 0), cmp_u64ipair());
					int dist = std::distance(record.begin(), lb);
					bool is_val_found = false;
					for (; ; ++ lb, ++ dist) {
						if (lb->first == val) {
							is_k_found[dist] = true;
							int fileID = lb->second;
							if (fileID >= num_node) abording ("batch_proc_seq SC failed");
							//++ private_abund[fileID];
							//fileIDs.insert(fileID);
							++ total_fIDs;
							is_val_found = true;
							auto it_fc = fileID_cnt.find(fileID);
							if (it_fc != fileID_cnt.end()) ++ it_fc->second;
							else fileID_cnt[fileID] = 1;
						} else break;
					}
					if (is_val_found) {
						++ num_k_sampled;
						is_frag_found[fragidx] = true;
						//k_pos += k;
						++ k_pos;
					} else ++ k_pos;

				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos += pos_n + 1;
					}
					else k_pos += k;
				}
			} // while (k_pos <= last_k_pos) {

			// --- assign fragment to node according to fileID_cnt and num_k_sampled ---
			iset_t fileIDs;
			ivec_t best_map_ids;
			for (auto& x: fileID_cnt) {
				if (x.second == num_k_sampled) best_map_ids.push_back(x.first);
				x.second *= (1000/total_fIDs);
				fileIDs.insert(x.first);
			}

			if (best_map_ids.size()) {
				fileIDs.clear();
				for (auto& x: best_map_ids) {
					fileIDs.insert(x);
					if (best_map_ids.size() > 1) { // multiple taxa
						private_abund[x].second += 1000/best_map_ids.size();
					} else { // single taxa (unique counts)
						private_abund[x].first += 1000/best_map_ids.size();
					}
				}
			} else {
				for (auto& x: fileID_cnt) {
					private_abund[x.first].second += x.second;
				}
			}


			/*// accumulate every kmer along every frag to hits
			u64vec_t kmers;
			xny::get_bitkmer<std::back_insert_iterator<u64vec_t>, uint64_t>
				(frags[fragidx], std::back_inserter(kmers), k, 0);
			u64set_t kmerset (kmers.begin(), kmers.end());

			iset_t fileIDs;

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
						if (fileID >= num_node) abording ("batch_proc_seq SC failed");
						++ private_abund[fileID];
						fileIDs.insert(fileID);
						is_frag_found[fragidx] = true;
					} else break;
				}
			} */

			frag_nodes[fragidx] = fileIDs;
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

		#pragma omp critical
		{
			for (int i = 0; i < num_node; ++ i) {
				abund[i].unique += private_abund[i].first;
				abund[i].multi += private_abund[i].second;
			}
		}
	}

} // batch_proc_seq
