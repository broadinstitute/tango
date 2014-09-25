//========================================================================
// Project     : SkeletonDB
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
#include "xny/file_manip.hpp"

class ProcDB {

public:

	ProcDB (const double& perc, int k, int wd, int wl, int minl,
			int batch, bool quiet): perc_ (perc), k_(k), wd_(wd),
			wl_(wl), minl_(minl), batch_(batch), silent_ (quiet) {};

	// identify database file list
	void set_db_filelist(const std::string& DBDir, const strvec_t& DBFileList,
			const std::string& odir);

	// files to record skeleton for each taxa
	void set_db_skeleton_filelist(const std::string& odir);

	// step-wise kmer sampling for each taxa
	void generate_skeleton_per_taxa(ivec_t& taxaKCnt);
	void output_db_skeleton_fl(const std::string& ofile, const ivec_t& taxaKCnt) {
		std::ofstream fh;
		xny::openfile(fh, ofile);
		int index = 0;
		for (auto& x: skeletonFileList_) {
			fh << x << "\t" << taxaKCnt[index] << "\n";
			++ index;
		}
		xny::closefile (fh);
	}
	void load_db_skeleton_fl(ivec_t& taxaKCnt, const std::string& ifile) {
		skeletonFileList_.clear();
		std::ifstream fh;
		xny::openfile(fh, ifile,  std::ios::in);
		std::string fname;
		int cnt;
		while (fh >> fname >> cnt) {
			skeletonFileList_.push_back(fname);
			taxaKCnt.push_back(cnt);
		}
		xny::closefile (fh);
	}


	// ----------- debug code ------------------
	void debug_print_skeleton(const ivec_t& gis);


private:
	double perc_;		// percentage of kmers to be sampled
	int k_;				// length of kmer
	int wd_;
	int wl_;
	int minl_;			//min len of contig to be kept
	int batch_;
	bool silent_;
	std::vector<istrvecpair_t> taxa_names_; // vector<taxaID, names>
	std::vector<iistrtuple_t> tree_;		// vecotr<tID, parent_tID, myTaxaRank:species,genus...>
	std::vector<strvec_t> fileList_; // database file list
	strvec_t ofileList_; // output database file list
	ivec_t taxaIDs_;					// IDs of taxa
	strvec_t skeletonFileList_; // skeleton file list


	//--------------------- handle [taxaIDs_] -------------------------
	int get_num_taxa_() { return taxaIDs_.size(); }
	void add_taxaID (int id) { taxaIDs_.push_back(id); }

	//----- functions related to generate_skeleton_per_taxa() -----
	void step_wise_sampling_kmer_from_file_(ivec_t& GIs,
			uu64vec_t& skeletons, u64set_t& sampled_kmers,
			std::ofstream& ofh, const std::string& file);
	void sample_kmer_from_contig_(uu64vec_t& low_density_k_windows,
			iivec_t& index_windows, const u64set_t& sampled_kmers,
			const std::string& contig);
	void identify_additional_indices_(ivec_t& idx_addition,
		const ivec_t& sampled_index, const u64vec_t& klist,
		bdeque_t& is_valid, int step, int contig_len, xny::low_complexity& lc);
	void identify_low_density_windows(iivec_t& windows,
			const ivec_t& coordinates);
	void output_contig(std::ofstream& ofh, int gi,
			const iivec_t& low_density_index_windows,
			const std::string& contig);


	void get_gi_from_header_(int& gi, const std::string& header) {
		strvec_t entries;
		split('|', header, std::back_inserter(entries));
		if (entries.size() < 2) abording ("get_gi_from_header_ SC failed");
		gi = std::atoi (entries[1].c_str());
	}
	void output_skeletons_(const ivec_t& GIs, const uu64vec_t& skeletons,
			const std::string& ofile);

};

// comparator for absolute value
struct cmp_abs{
public:
	bool operator () (const int& lhs, const int& rhs)
	const {	return std::abs(lhs) < abs(rhs); }
};

#endif /* PROC_DB_H_ */
