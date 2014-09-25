//========================================================================
// Project     : PrepDB
// Name        : proc_tax_tree.h
// Author      : Xiao Yang
// Created on  : Oct 23, 2013
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


#ifndef PROC_TAX_TREE_H_
#define PROC_TAX_TREE_H_

#include "xutil.h"

void proc_tax_tree (std::vector<istrpair_t>& tree, const std::string& ifile,
		const std::string& level, int batch, bool silent);

void get_tax_types (int& maxNodeID, strset_t& allRanks,
		const std::string& ifile, int batch);

void batch_get_tax_types (int& maxNodeID, strset_t& allTaxRanks,
		const std::vector<strvec_t>& rows);

void gen_tree (std::vector<istrpair_t>& tree, const std::string& ifile,
		int batch);

void batch_gen_tree (std::vector<istrpair_t>& tree, const std::vector<strvec_t>& rows);

void rank_ordering (std::vector<strvec_t>& orderedRanks,
		const strvec_t& ranks, const std::vector<istrpair_t>& tree);

void generate_order (std::vector<strvec_t>& orderedRanks,
		const strvec_t& ranks, const bbvec_t& array);

#endif /* PROC_TAX_TREE_H_ */
