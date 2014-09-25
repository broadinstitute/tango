//========================================================================
// Project     : PrepDB
// Name        : proc_gi_node.cpp
// Author      : Xiao Yang
// Created on  : Oct 22, 2013
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

#include "proc_gi_node.h"

/* @brief	Given file in the format of [gi]\t[node], obtain all gis and
 * sort by value
 */
void get_gis (ivec_t& gis, const std::string& ifile, int batch, bool silent) {
	batch *= 10;
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);

	int cnt = 0;
	/*int gi, node;
	while (ifh >> gi >> node) {
		gis.push_back(gi);
		if (cnt % 2000000 == 0) std::cout << cnt << "\n";
		++ cnt;
	}*/

	while (ifh.good()) {

		strvec_t rows;
		read_file_chunk (rows, ifh, batch);
		int sz = rows.size();
		cnt += sz;
		if (!silent) std::cout << "\t " << cnt << " read\n";
		#pragma omp parallel
		{
			ivec_t private_gis;
			#pragma omp for
			for (int i = 0; i < sz; ++ i) {
				// split
				strvec_t entries;
				split('\t', rows[i], std::back_inserter(entries));
				if (entries.size() < 1) abording ("get_gis SC failed");
				private_gis.push_back(atoi(entries[0].c_str()));
			}
			#pragma omp critical
			{
				gis.insert(gis.end(), private_gis.begin(), private_gis.end());
			}
		}

	}
	ifh.close();
	if (!silent) std::cout << "\t" << gis.size() << " gi found\n";
	std::sort (gis.begin(), gis.end());
} // get_gis

/* @brief	Read input file in the format of [gi]\t[nodeID] (assuming
 * it is sorted with respect to [gi], then
 * 1) if ofile is specified, write to ofile in the same format as input
 * ([gi_0][nodeID], [gi_1][nodeID], ...) such that any gi number in between
 * [gi_j, gi_j+1) have the same [nodeID]
 * 2) if ofile is not specified, read input
 */
void get_giRange_node (std::vector<ipair_t>& vec_gi_node,
		const std::string& ifile, const std::string& ofile,
		int batch, bool silent) {

	bool is_to_read = true;

	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);
	std::ofstream ofh;
	if (!ofile.empty()) {
		ofh.open (ofile.c_str(), std::ofstream::out);
		if (!ofh.good()) abording ("Can't open " + ofile);
		is_to_read = false;
	}

	int num_node = 0, num_row = 0;
	while (ifh.good()) {
		std::vector<strvec_t> rows;
		read_in_tab_file (rows, ifh, batch);
		batch_get_giRange_node (vec_gi_node, rows, is_to_read);
		num_node += vec_gi_node.size();
		num_row += rows.size();
		if (!is_to_read) {
			for (auto& x : vec_gi_node) {
				ofh << std::get<0>(x) << "\t" << std::get<1>(x) << "\n";
			}
			vec_gi_node.clear();
		}
	}

	if (!silent) std::cout << "\tnums of rows, ranges: " << num_row << ", ";
	if (is_to_read) std::cout << vec_gi_node.size() << std::endl;
	else	  std::cout << num_node << std::endl;
	ifh.close();
	if (!ofile.empty()) ofh.close();

} //get_giRange_node

void batch_get_giRange_node (std::vector<ipair_t>& vec_gi_node,
		const std::vector<strvec_t>& rows, bool is_to_read) {

	if (rows.size() < 1 || rows[0].size() < 2) return;
	if (is_to_read) {
		//vec_gi_node.reserve(rows.size());
		for (auto& x: rows) {
			vec_gi_node.push_back(ipair_t(std::atoi(x.front().c_str()),
										  std::atoi(x.back().c_str())));
		}
	} else {
		int prev_start_gi = std::atoi(rows[0][0].c_str()),
			prev_node = std::atoi(rows[0][1].c_str());

		for (int i = 1; i < (int) rows.size(); ++ i) {
			if (rows[i].size() < 2) {
				std::cout << "\n[ERR] proc_gi_node_batch: entries in row " << i
						<< " < 2\n";
				exit(1);
			}
			int cur_node = std::atoi(rows[i][1].c_str());
			if (cur_node != prev_node) {
				vec_gi_node.push_back(
						ipair_t (prev_start_gi, prev_node));
				prev_start_gi = std::atoi(rows[i][0].c_str());
				prev_node = cur_node;
			}
		}
		// the last entry
		vec_gi_node.push_back(ipair_t (prev_start_gi, prev_node));
	}
} // batch_get_giRange_node
