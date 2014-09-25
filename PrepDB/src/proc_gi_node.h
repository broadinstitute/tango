//========================================================================
// Project     : PrepDB
// Name        : proc_gi_node.h
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


#ifndef PROC_GI_NODE_H_
#define PROC_GI_NODE_H_

#include "xutil.h"

void get_gis (ivec_t& gis, const std::string& ifile, int batch, bool silent);

void get_giRange_node (std::vector<ipair_t>& vec_gi_node,
		const std::string& ifile, const std::string& ofile,
		int batch, bool silent);

void batch_get_giRange_node (std::vector<ipair_t>& vec_gi_node,
		const std::vector<strvec_t>& rows, bool is_to_read) ;

#endif /* PROC_GI_NODE_H_ */
