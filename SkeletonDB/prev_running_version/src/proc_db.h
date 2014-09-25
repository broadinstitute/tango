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
#include "xny/file_manip.hpp"

class ProcDB {

public:

	ProcDB (const double& perc, int k, int batch, bool silent):
		perc_ (perc), k_(k), batch_(batch), silent_ (silent) {};

	// identify database file list
	void set_db_filelist(const std::string& DBDir, const strvec_t& DBFileList);

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

	// import names of taxa, create [taxa_names_] (vector<taxaID, names>)
	void import_taxa_names(const std::string& ifile);

	// create [tree_] -- vector<tID, parent_tID, myTaxaRank:species,genus...>
	void import_taxonomy_tree(const std::string& ifile);

	// count sampled kmers that occur in reads
	void count_sampled_kmers_in_reads (const ivec_t& taxaKCnt,
			const strvec_t& ipfqs, const strvec_t& isfqs);

	// generate kmer pairs for each taxa
	void generate_locality_kpairs_per_taxa(
			std::vector<uu64pair_t>& taxaKPairs, int w);

	// generate kmer pairs for each fragment
	void generate_locality_kpairs_in_reads(
			std::vector<uu64pair_t>& rKPairs,
			int min_freq,
			const strvec_t& ipfqs,
			const strvec_t& isfqs);

	//Identify connectivity of kmer pairs in [taxaKPairs] relying
	// on [rKPairs] and set flag field in [kmerInfo_]
	void set_kmer_multiplicity(const std::vector<uu64pair_t>& taxaKPairs,
	 			const std::vector<uu64pair_t>& rKPairs);

	void detect_spurious_taxa(int max_perc, const std::string& odir,
			const std::string& oprefix);

	// assign scores to fragments
	void	 assign_frags(const strvec_t& ipfqs,	const strvec_t& isfqs,
			const std::string& odir, const std::string& oprefix);

	// ----------- debug code ------------------
	void debug_print_skeleton(const ivec_t& gis);


private:
	double perc_;		// percentage of kmers to be sampled
	int k_;				// length of kmer
	int batch_;
	bool silent_;
	std::vector<istrvecpair_t> taxa_names_; // vector<taxaID, names>
	std::vector<iistrtuple_t> tree_;		// vecotr<tID, parent_tID, myTaxaRank:species,genus...>
	std::vector<strvec_t> fileList_; // database file list
	ivec_t taxaIDs_;					// IDs of taxa
	strvec_t skeletonFileList_; // skeleton file list
	std::vector<u64ipair_t> kmerTaxaID_; // sorted list of <kmer, taxaID> -- can be a class itself
	std::vector<u64ibtuple_t> kmerInfo_; //vector <kmer, cnt_in_reads, is_m_flag>


	//--------------------- handle [taxaIDs_] -------------------------
	int get_num_taxa_() { return taxaIDs_.size(); }
	void add_taxaID (int id) { taxaIDs_.push_back(id); }

	//-------------------- handle [kmerInfo_] -------------------------
	void rm_lowfreq_kmers_(int freq);
	// set the flag field to be true
	bool set_kmerinfo_flag_(const uint64_t& kval) {
		auto it = std::lower_bound(kmerInfo_.begin(), kmerInfo_.end(),
				std::make_tuple(kval, 0, false), cmp_u64ibtuple());
		if (std::get<0> (*it) == kval) std::get<2> (*it) = true;
		else return false;
		return true;
	}
	bool get_kmerinfo_flag_(const uint64_t& kval) {
		auto it = std::lower_bound(kmerInfo_.begin(), kmerInfo_.end(),
				std::make_tuple(kval, 0, false), cmp_u64ibtuple());
		if (std::get<0> (*it) == kval) return std::get<2> (*it);
		else return false;
	}
	int get_flagged_kmer_count_() {
		int cnt = 0;
		for (auto& x: kmerInfo_) if (std::get<2>(x)) ++ cnt;
		return cnt;
	}
	int get_kmerinfo_freq_(const uint64_t& kmer) {
		auto it = std::lower_bound(kmerInfo_.begin(), kmerInfo_.end(),
				std::make_tuple(kmer, 0, false), cmp_u64ibtuple());
		if (std::get<0> (*it) == kmer) return (std::get<1> (*it));
		else return -1;
	}
	int get_kmerinfo_sz_() { return kmerInfo_.size(); }

	//-------------------- handle [kmerTaxaID_] -------------------------
	void add_to_kmertaxaid_(const std::vector<u64ipair_t>& rhs) {
		kmerTaxaID_.insert(kmerTaxaID_.end(), rhs.begin(), rhs.end());
	}
	void rm_dupl_kmertaxaid_() {	rm_duplicated_elements (kmerTaxaID_); }
	int get_kmertaxaid_sz_() { return kmerTaxaID_.size(); }
	/* Given a kmer, identify all taxaIDs it belongs to */
	iset_t get_kmertaxaid_tids_(const uint64_t& val) {
		iset_t ids;
		auto it = std::lower_bound(kmerTaxaID_.begin(),
			kmerTaxaID_.end(), u64ipair_t(val, 0), cmp_u64ipair());
		for (; it != kmerTaxaID_.end(); ++ it) {
			if (it->first == val) ids.insert (it->second);
			else break;
		}
		return ids;
	}
	bool is_kmer_in_taxa_(const uint64_t& val, int index_taxa) {
		auto it = std::lower_bound(kmerTaxaID_.begin(), kmerTaxaID_.end(),
				u64ipair_t (val, 0), cmp_u64ipair());
		for (; it != kmerTaxaID_.end(); ++ it) {
			if (it->first == val) {
				if (it->second == index_taxa) return true;
			} else return false;
		}
		return false;
	}

	//----- functions related to generate_skeleton_per_taxa() -----
	void step_wise_sampling_kmer_from_file_(ivec_t& GIs, uu64vec_t& skeletons,
			u64set_t& sampled_kmers, const std::string& file);
	void sample_kmer_from_contig_(u64vec_t& kmers,
			const u64set_t& sampled_kmers, const std::string& contig);
	void get_gi_from_header_(int& gi, const std::string& header) {
		strvec_t entries;
		split('|', header, std::back_inserter(entries));
		if (entries.size() < 2) abording ("get_gi_from_header_ SC failed");
		gi = std::atoi (entries[1].c_str());
	}
	void output_skeletons_(const ivec_t& GIs, const uu64vec_t& skeletons,
			const std::string& ofile);


};

#endif /* PROC_DB_H_ */
