//========================================================================
// Project     : RankAbundance
// Name        : RankAbundance.cpp
// Author      : Xiao Yang
// Created on  : Nov 1, 2013
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
// Prerequisite:
// 1) install boost system, filesystem library
// 2) c++0x (gcc4.7.2+)
//========================================================================


#include "Parameter.h"
#include "proc_db.h"
#include "debugcode.h"
#include "boost/filesystem.hpp"

int main (int argc, char** argv) {
	/*test of get_neigbhorset_in_pairs
	std::vector<uu64pair_t> test {uu64pair_t (1, 3), uu64pair_t (2, 3), uu64pair_t (2, 3),
		uu64pair_t (2, 4), uu64pair_t (5, 2)};
	u64set_t nbs = get_neighborset_in_pairs(uint64_t(2), test, cmp_uu64pair());

	for (auto& x: nbs) std::cout << x << "\n";
	exit(1); */
	/* test of cmp_uu64ituple ()
	std::vector<uu64ituple_t> test { std::make_tuple (2, 3, 0), std::make_tuple (2, 3, 0),
		std::make_tuple (2,1,4), std::make_tuple (1, 1, 44)};

	std::sort(test.begin(), test.end(), cmp_uu64ituple());
	auto it = std::unique_copy (test.begin(), test.end(), test.begin());
	test.resize (std::distance (it, test.begin()));
	for (auto& x: test) {
		std::cout << std::get<0> (x) << ", " << std::get<1> (x) << "\n" ;
	}
	exit (1); */
	/*u64vec_t test {1133239868604512267, 2579872146756650429,
		2763909787844936593, 12221200084362843037};
	for (auto& x: test) std::cout << xny::ID2Str(x, 32) << "\n";
	exit (1); */

	Parameter myPara;

	myPara.procinput(argc, argv);

 	omp_set_num_threads(myPara.get_num_threads());

 	int batch = 1000000;
 	ProcDB myDB (myPara.get_k(), batch, myPara.quiet());


 	// test proc_contig_context
 	/*
 	std::vector<ktuple_t> rhs {std::make_tuple (1, 0, false), std::make_tuple (2, 0, true) };
 	u64vec_t skeleton {3, 3, 1, 1, 3, 3, 2 , 3, 3, 3, 1, 3, 1, 3, 3};
	myDB.set_kmerinfo_ (rhs) ;
	uu64vec_t blocks;
	int min_block_sz = 2, gap = 3;
	myDB.identify_valid_skeleton_blocks_ (blocks, skeleton, min_block_sz, gap);
	exit (1);
	*/
 	/*
 	u64vec_t kvals { 1, 3, 2, 1, 2, 5, 6, 6, 0, 9, 9, 9};
 	std::vector<ktuple_t> rhs;
 	rhs.push_back(std::make_tuple (1, 3, false));
 	rhs.push_back(std::make_tuple (5, 2, false));
 	rhs.push_back(std::make_tuple (7, 2, false));
 	rhs.push_back(std::make_tuple (8, 2, false));
	myDB.set_kmerinfo (rhs);
 	myDB.accrue_kmercnts(kvals);
 	myDB.debug_prnt_kmerinfo();
 	exit(1); */

 	//myDB.set_db_filelist(myPara.get_db_path(), myPara.get_dbfl());

 	// directly load data structure from file instead of processing db again
 	//myDB.debug_load_kmerTaxaID_from_file ("/seq/viral/analysis/xyang/FUO/scripts/RankAbundance/kmertaxaid.txt");
 	ivec_t taxaKCnt;

 	/*if (! boost::filesystem::exists(myPara.get_iskeletonfile())) {
 	 	myDB.set_db_skeleton_filelist(myPara.get_odir(), myPara.get_oprefix());
 		myDB.generate_skeleton_per_taxa(taxaKCnt);
 	 	myDB.output_db_skeleton_fl(myPara.get_iskeletonfile(), taxaKCnt);
 		return EXIT_SUCCESS;
 	} else*/

 	myDB.load_db_skeleton_fl(taxaKCnt, myPara.get_iskeletonfile());
 	myDB.import_taxa_names(myPara.get_inodename());
 	myDB.import_taxonomy_tree(myPara.get_itreefile());

 	// identify skeleton kmers in the read data
 	myDB.count_sampled_kmers_in_reads(taxaKCnt, myPara.get_ipfqs(),
 			myPara.get_isfqs());

 	// generate kmer pairs in taxa
 	std::vector<uu64pair_t> taxaKPairs; //
 	myDB.generate_locality_kpairs_per_taxa(taxaKPairs, myPara.get_w());

  	std::vector<uu64pair_t> rKPairs;
 	myDB.generate_locality_kpairs_in_reads(rKPairs, myPara.get_mcov(),
 			myPara.get_ipfqs(), myPara.get_isfqs());

 	myDB.set_kmer_multiplicity(taxaKPairs, rKPairs);

	int max_perc_flagged_kmers = 95;
 	myDB.detect_spurious_taxa(max_perc_flagged_kmers, myPara.get_odir(),
 			myPara.get_oprefix());

 	myDB.assign_frags(myPara.get_ipfqs(), myPara.get_isfqs(),
 			myPara.get_odir(), myPara.get_oprefix());

 	// 33917: Microbacterium barkeri, Micrococcus luteus: 1270
 	// 562: ecoli, 294: Pseudomonas fluorescens, 28901: Salmonella enterica
 	std::string opath = myPara.get_odir() + "/";
 	opath += myPara.get_oprefix() + "_skelton.txt";
 	myDB.debug_print_skeleton(opath, {28901, 1270, 33917});
 	//myDB.debug_print_skeleton({562, 621, 623, 33169, 33917});
 	//myDB.set_sampled_kmercnt_per_taxa (); // [taxaInfo_]->k_sampled

	return (EXIT_SUCCESS);


 	//myDB.link_sampled_kmers_via_reads (myPara.get_ipfqs(),  myPara.get_isfqs());
 	//myDB.rm_sampled_kmers_absent_in_reads ();

 	//myDB.debug_output_kmerTaxaID_to_file ("/seq/viral/analysis/xyang/FUO/scripts/RankAbundance/kmertaxaid.txt");
 	//exit (1);

 	//myDB.rm_spurious_support_per_taxa ();

 	//std::cout << "debug done\n";
 	//exit (1);

	/*
 	std::vector<u64ipair_t> record; // <kmer, fileID>
	std::vector<strvec_t> input_files;
	ivec_t nodes;
	std::vector<clade_t> cladeinfo;
	if (!myPara.get_db_path().empty() || !myPara.is_idbfl_empty()) {
		if (!myPara.is_no_prt()) {
			std::cout << "Read DB...\n";
			if (!myPara.get_db_path().empty()) {
				std::cout << "\tRead dir: " << myPara.get_db_path().empty() << "\n";
			} else std::cout << "\tRead files: ";
			for (auto& x: myPara.get_dbfl()) std::cout << "\t" << x;
			std::cout << std::endl;
		}

		get_input_db_files (input_files, nodes, myPara.get_db_path(), myPara.get_dbfl());


		if (!myPara.is_no_prt()) 	std::cout << "Step-wise sampling DB ...\n";
		db_sample (record, myPara.get_perc(), myPara.get_k(),
				input_files, myPara.is_no_prt());

		// update cladeinfo->k_sampled
		cladeinfo.resize(nodes.size());
		for (auto& x: record) {
			if (x.second < nodes.size()) ++ cladeinfo[x.second].k_sampled;
			else abording ("main SC0 failed");
		}
	}

	// [debug] analyze record and generate kmer freq histogram
	//debug_gen_kmer_histogram (record);

	if (!myPara.is_ipfqs_empty() || !myPara.is_isfqs_empty()) {
		int batch = 1000000;
		bvec_t is_record_m (record.size(), false);

		if (!myPara.is_no_prt()) 	std::cout << "Proc read support ...\n";
		read_context (record, is_record_m, myPara.get_ipfqs(), myPara.get_isfqs(),
				myPara.get_k(), batch, myPara.is_no_prt());

		//[may not be necessary]
		if (!myPara.is_no_prt()) 	std::cout << "Proc DB context ...\n";
		db_context (is_record_m, record, myPara.get_perc(), myPara.get_k(),
				input_files, myPara.is_no_prt());

		// update cladeinfo->k_in_reads
		for (auto& x: record) {
			if (x.second < nodes.size()) ++ cladeinfo[x.second].k_in_reads;
			else abording ("main SC1 failed");
		}

		if (!myPara.is_no_prt()) 	std::cout << "Measure abundance ...\n";
		read_count (cladeinfo, record, is_record_m, myPara.get_ipfqs(),
				myPara.get_isfqs(), myPara.get_k(), batch, myPara.is_no_prt());
	}

	exit(1);


	std::vector<abund_t> abund (nodes.size()); // fileID [cnt]
	std::vector<ipair_t> abund_fID;
	bvec_t is_k_found (record.size(), false);
	if (!myPara.is_ipfqs_empty() || !myPara.is_isfqs_empty()) {
		int batch = 1000000;
		proc_seq (abund, is_k_found, record, myPara.get_ipfqs(), myPara.get_isfqs(),
				myPara.get_k(), batch, myPara.is_no_prt());
		abund_fID.resize(abund.size());
		for (int i = 0; i < abund.size(); ++ i) {
			abund_fID[i] = ipair_t (abund[i].unique, i);
		}
		std::sort (abund_fID.begin(), abund_fID.end());
	}

	std::vector<istrvecpair_t> node_names;
	if (!myPara.get_inodename().empty()) {
		if (!myPara.is_no_prt()) std::cout << "Read " << myPara.get_inodename() << "\n";
		get_node_names (node_names, myPara.get_inodename());
		std::sort (node_names.begin(), node_names.end(), cmp_istrvecpair());
	}

	for (auto& x: abund_fID) {
		int fID = x.second;
		int node = nodes[fID];
		int dist;
		if (node_names.size()) {
			auto lb = std::lower_bound(node_names.begin(), node_names.end(),
					istrvecpair_t (node, strvec_t()), cmp_istrvecpair());
			dist = std::distance(node_names.begin(), lb);
			if (dist >= node_names.size()) {
				abording ("main SC0 failed");
			} else if (node_names[dist].first != node) {
				std::cout << node_names[dist].first << ", " << node << "\n";
				abording ("main SC1 failed");
			}
		}
		std::cout << std::setw(10) << nodes[fID]
		          << std::setw(10) << abund[fID].clade_k_num
		          << std::setw(10) << abund[fID].read_k_num
		          << std::setw(10) << abund[fID].unique/1000
				  << std::setw(10) << abund[fID].multi/1000;
		if (node_names.size()) {
			for(auto& s: node_names[dist].second) std::cout << " -- " << s;
		}
		std::cout << "\n\n";



		//if (abund[fID].read_k_num == 0) std::cout << std::setw(10) << "-" << "\n";
		//else std::cout << std::setw(10) << abund[fID].total/abund[fID].read_k_num << "\n";
	}

	iset_t debug_nodes {33917}; //, 1029823, 205844};
	debug_print_node_kmership (record, is_k_found, debug_nodes, nodes);
	return (EXIT_SUCCESS);
	*/
}



