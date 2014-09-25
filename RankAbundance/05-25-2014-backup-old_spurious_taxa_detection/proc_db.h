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

	ProcDB (int k, int batch, bool silent):
		k_(k), batch_(batch), silent_ (silent) {};

	// read input skeleton file list and kmers per file, identify taxaIDs
	void load_db_skeleton_fl(ivec_t& taxaKCnt, const std::string& ifile);

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
	void debug_print_skeleton(const std::string& opath, const ivec_t& gis);

	// link sampled kmers that occur in reads
	//void link_sampled_kmers_via_reads (const strvec_t& ipfqs, const strvec_t& isfqs);

	// remove sampled kmers that do not occur in reads
	//void rm_sampled_kmers_absent_in_reads ();

	// using context info to identify/rm contig regions that are unlikely
	// supported by read data even though kmers have been observed
	//void rm_spurious_support_per_taxa ();

	//void accrue_read_support_per_taxa ();

	// count initial # of sampled kmers per taxa
	/*
	void set_sampled_kmercnt_per_taxa () {
		taxaInfo_.resize(get_num_taxa_());
		for (auto& x: kmerTaxaID_) {
			if (x.second < get_num_taxa_()) ++ taxaInfo_[x.second].k_sampled;
			else abording ("ProcDB::get_taxa_sampled_kmercnt SC failed");
		}
	}*/

	/*
 	void debug_output_db_taxaKCnt (const ivec_t& taxaKCnt, const std::string& ofile){
 		std::ofstream fh;
		xny::openfile(fh, ofile);
		for (auto& x: taxaKCnt) fh << x << "\n";
 		xny::closefile (fh);
 	}
 	void debug_load_db_taxaKCnt (ivec_t& taxaKCnt, const std::string& file){

 		std::ifstream fh;
		xny::openfile (fh, file, std::ios::in);
		int cnt;
		while (fh >> cnt) taxaKCnt.push_back(cnt);
 		xny::closefile (fh);
 	}*/

private:
	//double perc_;		// percentage of kmers to be sampled
	int k_;				// length of kmer
	int batch_;
	bool silent_;
	std::vector<istrvecpair_t> taxa_names_; // vector<taxaID, names>
	std::vector<iistrtuple_t> tree_;		// vector<tID, parent_tID, myTaxaRank:species,genus...>
	//std::vector<strvec_t> fileList_; // database file list
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
	//void step_wise_sampling_kmer_from_file_(ivec_t& GIs, uu64vec_t& skeletons,
	//		u64set_t& sampled_kmers, const std::string& file);
	//void sample_kmer_from_contig_(u64vec_t& kmers,
	//		const u64set_t& sampled_kmers, const std::string& contig);
	/*void get_gi_from_header_(int& gi, const std::string& header) {
		strvec_t entries;
		split('|', header, std::back_inserter(entries));
		if (entries.size() < 2) abording ("get_gi_from_header_ SC failed");
		gi = std::atoi (entries[1].c_str());
	}*/
	//void output_skeletons_(const ivec_t& GIs, const uu64vec_t& skeletons,
	//		const std::string& ofile);

	//----- functions related to count_sampled_kmers_in_reads() -----
	void read_skeletons_(ivec_t& GIs, uu64vec_t& skeletons,
			const std::string& ifile);
	void read_skeletons_(u64vec_t& kmers, const std::string& ifile);
	void substract_counted_kmers_(u64vec_t& kmers);

	//----- functions related to generate_locality_kpairs_per_taxa() -----
	void identify_valid_skeleton_blocks_(uu64vec_t& blocks,
			const u64vec_t& skeleton, int min_block_sz, int gap);
	void proc_blocks_(std::vector<u64ipair_t>& ktaxa,
			std::vector<uu64ituple_t>& kpairtaxa,
			const uu64vec_t& blocks,
			int w,
			int taxaID);
	void proc_kpairtaxa_(std::vector<uu64pair_t>& taxakpairs,
			std::vector<uu64ituple_t>& kpairtaxa);

	//--------------- functions related to assign_frags() ---------------
	void output_taxa_list_(const std::vector<i4tuple_t>& taxarank,
			const std::string& odir, const std::string& oprefix);
	void build_tree_(imap_t& edges, const std::vector<i4tuple_t>& taxarank);
	void output_taxa_tree_(const imap_t& tree,
			const std::vector<i4tuple_t>& taxarank,
			const std::string& odir, const std::string& oprefix);
	void accumulate_scores_(imap_t& node_score, const imap_t& edges,
			const std::vector<i4tuple_t>& taxarank);
	std::string get_taxa_names_(int taxaID);
	std::string get_taxa_rank_(int taxaID);


	// ---------- debug code ------------------------
	void debug_prnt_kmerinfo () {
		for (auto& x: kmerInfo_) {
			std::cout << std::get<0> (x) << "\t" <<  std::get<1> (x) << "\t"
					<< std::get<2> (x) << "\n";
		}
	}

	/*
	void debug_output_kmerTaxaID_to_file (const std::string& ofile) {
		std::ofstream ofh (ofile.c_str());
		if (!ofh.good()) {
			std::cout << "cannot open " << ofile << "\n";
			exit (1);
		}
		int cnt = 0;
		for (auto& x: kmerTaxaID_) {
			ofh << x.first << "\t" << x.second << "\n";
			++ cnt;
		}
		ofh.close();
		std::cout <<  "\t" << cnt << " entries written to " << ofile << "\n";
	}
	void debug_load_kmerTaxaID_from_file (const std::string& ifile) {
		std::ifstream ifh (ifile.c_str());
		if (!ifh.good()) {
			std::cout << "cannot open " << ifile << "\n";
			exit (1);
		}
		kmerTaxaID_.clear();
		uint64_t kval;
		int taxaID;
		int cnt = 0;
		while (ifh >> kval >> taxaID) {
			kmerTaxaID_.push_back(u64ipair_t (kval, taxaID));
			++ cnt;
		}
		ifh.close();
		std::cout <<  "\t" << cnt << " entries read from " << ifile << "\n";
	}*/

	//void set_kmerinfo_ (const std::vector<u64ibtuple_t>& rhs) { kmerInfo_ = rhs; }
	//std::vector<u64bpair_t> kmerStat_; // vector <kmer, is_m_flag>
	//std::vector<taxa_t> taxaInfo_;		// taxa count information
	//void count_kmers_in_read_files_(const u64vec_t& kmers,	const strvec_t& fl, bool is_pe = true);
	//void count_kmers_in_frags (const u64vec_t& kmers, const strvec_t& frags);
	//void count_sampled_kmer_in_frags (const strvec_t& frags); // obsolete
	//void get_taxaIDs (iset_t& taxaIDs, const uint64_t& kval);
	//void link_sampled_kmer_in_frags (const strvec_t& frags);
	//void accrue_kcnts_(u64vec_t& kvals);
	//void accrue_kstat (std::vector<u64bpair_t>& kstat);
	//void set_kstat (const std::vector<u64bpair_t>& rhs) { kmerStat_ = rhs; }
	/*
	void print_kstat () {
		int num_m_stat = 0;
		for (auto& x: kmerStat_) if (x.second) ++ num_m_stat;
		std::cout << "\t# kmers, m-labled, u-labled:" << kmerStat_.size() << ", " <<
			num_m_stat << ", " << kmerStat_.size() - num_m_stat << "\n";
	}*/
	//void proc_contig_context (std::vector<u64bpair_t>& valid_kmers,
	//	u64set_t& kmer_to_del, const u64vec_t& skeleton);
	//void confirm_deletion (std::vector<u64bpair_t>& valid_kmers,
	//	u64set_t& kmer_to_del);

	//void count_kpair_in_read_files_(std::vector<uu64ituple_t>& rkpcnt,
	//		int min_freq, const strvec_t& fl, bool is_pe = true);
	//void count_kpair_in_frags_(std::vector<uu64pair_t>& rkp,
	//		int min_freq, const strvec_t& frags);
	//void accrue_kpcnts_(std::vector<uu64ituple_t>& rkpcnt,
	//		std::vector<uu64pair_t>& rkp);


	// obsolete functions
	//void add_kmer_from_contig (u64vec_t& kmers, const std::string& contig);
	//void sample_kmer_from_file_(ivec_t& GIs, uu64vec_t& skeletons,
	//		const std::string& file);
};


/*
struct taxa_t {
	int k_sampled;		//# kmers selected for the taxa
	int k_in_reads;		//# kmers in clade present in read data
	int u_cnt;			// # unique-aligned fragment assignment
	taxa_t () {
		k_sampled = 0;
		k_in_reads = 0;
		u_cnt = 0;
	}
};*/
/*
void get_node_names (std::vector<istrvecpair_t>& node_names, const std::string& ifile);

void db_sample (std::vector<u64ipair_t>& record, const double& perc, int k,
		const std::vector<strvec_t>& input_files, bool silent);


void db_context (bvec_t& is_record_m, std::vector<u64ipair_t>& record,
		const double& perc, int k, const std::vector<strvec_t>& input_files,
		bool silent);

void proc_contig_skeleton (bvec_t& is_record_m, bvec_t& is_to_del,
		std::vector<u64ipair_t>& record, const u64vec_t& kmers, int id);

void get_input_db_files (std::vector<strvec_t>& input_files,
		ivec_t& nodes, const std::string& idir, const strvec_t& ifl);

void sample_kmers (u64vec_t& kmers, int k, const double& perc,
		const std::string& file);

void sample_contig (u64vec_t& kmers, int k, const double& perc,
		const std::string& contig);
*/

#endif /* PROC_DB_H_ */
