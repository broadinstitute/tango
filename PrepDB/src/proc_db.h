//========================================================================
// Project     : PrepDB
// Name        : proc_db.h
// Author      : Xiao Yang
// Created on  : Oct 28, 2013
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


#ifndef PROC_DB_H_
#define PROC_DB_H_

#include "xutil.h"
#include <sys/resource.h>

void proc_db (const std::string& odir, const std::string& ogilevel,
	const std::string& idir, const ivec_t& gis,
	const std::vector<ipair_t>& vec_gi_node,
	const std::vector<istrpair_t>& tree,	bool silent);

void get_gi_len_fa (std::vector<ipair_t>& vec_gi_len, const std::string& file);

void get_gi (int& gi, const std::string& header);

std::string gen_filename (const std::string& dir, int id, int subid);

void output_db (const strvec_t& output_files, std::vector<iset_t>& giset,
		const imap_t& gi_clsId, const strvec_t& input_files);

/* sort <int, int> pair based on the first field */
struct cmp_ipair{
public:
	bool operator () (const ipair_t& lhs, const ipair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

#endif /* PROC_DB_H_ */
