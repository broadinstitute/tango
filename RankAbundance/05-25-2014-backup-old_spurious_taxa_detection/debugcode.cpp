//========================================================================
// Project     : RankAbundance
// Name        : debugcode.cpp
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


#include "debugcode.h"

/* @brief	Generate histogram for numbers in the cnts
 */
void output_histogram(const ivec_t& cnts) {
	std::vector<ipair_t> cnt_freq;
	int freq = 1;
	int cnt = -1;
	if (cnts.size()) cnt = cnts[0];
	else return;
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
	for (auto& x: cnt_freq) {
		std::cout << x.first << "\t" << x.second << "\n";
	}
}



void debug_print_record_info (const bvec_t& is_record_found,
		const bvec_t& is_record_m,
		const std::vector<u64ipair_t>& record) {
	int num = record.size();
	uint64_t prev_val;
	if (num) {
		prev_val = record.front().first;
		std::cout << prev_val;
		if (is_record_found.size() &&
				is_record_found.front()) std::cout << "\t+";
		else std::cout << "\t-";
		if (is_record_m.front()) std::cout << "\tm";
		std::cout << "\n";
	} else return;
	for (int i = 1; i < num; ++ i) {
		if (record[i].first != prev_val) {
			prev_val = record[i].first;
			std::cout << prev_val;
			if (is_record_found.size() &&
					is_record_found[i]) std::cout << "\t+";
			else std::cout << "\t-";
			if (is_record_m[i]) std::cout << "\tm";
			std::cout << "\n";
		}
	}
}

void debug_print_contig_skeleton (const bvec_t& is_record_m,
		const std::vector<u64ipair_t>& record, const u64vec_t& kmers) {
	int sz = kmers.size();
	for (int i = 0; i < sz; ++ i) {
		auto lb = std::lower_bound(record.begin(), record.end(),
				u64ipair_t (kmers[i], 0), cmp_u64ipair());
		int dist = std::distance(record.begin(), lb);
		if (lb->first == kmers[i]) { // found
			if (is_record_m[dist]) std::cout << "m";
			else std::cout << "u";
		} else std::cout << "-";
	}
}

void debug_print_node_kmership (const std::vector<u64ipair_t>& record,
		const bvec_t& is_k_found, const iset_t& debug_nodes, const ivec_t& nodes) {
	std::cout << "\nNodes to be search for: ";
	for (auto& x: debug_nodes) std::cout << "\t" << x;
	std::cout << "\n";

	std::map<int, iivec_t> nn_map;
	for (int i = 0; i < record.size(); ++ i) {
		if (is_k_found[i]) {
			int node = nodes[record[i].second];
			uint64_t kmer = record[i].first;
			auto it_dn = debug_nodes.find(node);
			if (it_dn != debug_nodes.end()) { // found target
				// go backwards and forwards

				//std::cout << "here\n";

				ivec_t shared_nodes;
				for (int j = i - 1; j >= 0; --j) {
					if (record[j].first == kmer) {
						shared_nodes.push_back(nodes[record[j].second]);
					} else break;
				}
				for (int j = i + 1; j < record.size(); ++ j) {
					if (record[j].first == kmer) {
						shared_nodes.push_back(nodes[record[j].second]);
					} else break;
				}
				auto it_nn = nn_map.find(node);
				if (it_nn != nn_map.end()) { //found
					it_nn->second.push_back(shared_nodes);
				} else nn_map[node] = iivec_t (1, shared_nodes);
			}
		}
	}

	// print out
	std::cout << "\n-----node_kmership----\n";
	for (auto& x: nn_map) {
		std::cout << "\nnode: " << x.first << "\n";
		for (auto& y: x.second) {
			if (!y.size()) std::cout << "-";
			else for (auto& z: y) std::cout << "\t" << z;
			std::cout << std::endl;

		}
	}
} // debug_print_node_kmership

