//========================================================================
// Project     : PrepDB
// Name        : proc_db.cpp
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


#include <sys/resource.h>
#include <dirent.h>
#include "boost/filesystem.hpp"
#include "proc_db.h"
#include "xny/seq_manip.hpp"
#include "ReadBioFile.h"
#include "CountKmers.h"
#include "CountKP.h"
#include "AssignFrag.h"
#include "debugcode.h"

/* Parse skeleton file in the format of
 * [filename][\t][# sampled kmers], where the filename is the path
 * to the skeleton file for a specific taxaID and is named in the
 * format of path/[taxaID].txt
 *
 * Store file path to [skeletonFileList_] and taxaID to [taxaIDs]
 */
void ProcDB::load_db_skeleton_fl(ivec_t& taxaKCnt, const std::string& ifile) {
	skeletonFileList_.clear();
	std::ifstream fh;
	xny::openfile(fh, ifile,  std::ios::in);
	std::string fname;
	int cnt;
	while (fh >> fname >> cnt) {
		skeletonFileList_.push_back(fname);
		taxaKCnt.push_back(cnt);

		std::string suffix = xny::get_suffix (fname, true, "/");
		add_taxaID (std::atoi(xny::get_prefix (suffix, true, ".").c_str()));
	}
	xny::closefile (fh);
}

/* @brief	Import to [taxa_names_] from a tsv file in the format of
 * [taxaID][\t][scientific name]
 */
void ProcDB::import_taxa_names(const std::string& ifile) {
	if (ifile.empty()) return;
	if (!silent_) std::cout << "\nImport taxa names from " << ifile << "...\n";
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("ProcDB::import_taxa_names: Can't open " + ifile);
	int cnt = 0;

 	strvec_t rows;
	read_file_chunk (rows, ifh, INT_MAX);

	int prev_taxaID;
	strvec_t names;
	for (int i = 0; i < rows.size(); ++ i) {
		strvec_t entries;
		split('|', rows[i], std::back_inserter(entries));
		std::string name = entries[1];
		if (name.length() > 2) name = name.substr(1, name.length() - 2);

		if (i == 0) { // first entry
			prev_taxaID = std::atoi(entries[0].c_str());
			names.push_back(name);
		} else {
			int taxaID = std::atoi(entries[0].c_str());
			if (taxaID != prev_taxaID) {
				taxa_names_.push_back(istrvecpair_t (prev_taxaID, names));
				prev_taxaID = taxaID;
				names = strvec_t (1, name);
			} else names.push_back(name);
		}
	}
	if (names.size()) taxa_names_.push_back(istrvecpair_t (prev_taxaID, names));
	std::sort (taxa_names_.begin(), taxa_names_.end(), cmp_istrvecpair());

	ifh.close();

	if (!silent_) std::cout << "\t" << taxa_names_.size() << " taxa read\n";

	/*for (auto& x: taxa_names_) {
		std::cout << x.first << "\n";
	}*/
} // get_taxa_names

/* @brief	Import tree structure from ifile and store the result in
 * [tree_] that has the structure of [taxaID] [parent_taxaID] [taxonomy_type]
 * sorted wrt to [taxaID]
 */
void ProcDB::import_taxonomy_tree(const std::string& ifile) {
	if (ifile.empty()) return;
	if (!silent_) std::cout << "\nImport taxonomy tree from " << ifile << "...\n";
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);
	int cnt = 0;

 	strvec_t rows;
	read_file_chunk (rows, ifh, INT_MAX);

	for (int i = 0; i < rows.size(); ++ i) {
		strvec_t entries;
		split('|', rows[i], std::back_inserter(entries));
		if (entries.size() < 3) abording ("ProcDB::import_taxonomy_tree SC failed");
		tree_.push_back(iistrtuple_t (atoi(entries[0].c_str()),
				atoi(entries[1].c_str()), entries[2]));
	}
	ifh.close();

	std::sort (tree_.begin(), tree_.end(), cmp_iistrtuple ());
	if (!silent_) std::cout << "\t" << tree_.size() << " entries read\n";

} // ProcDB::import_taxonomy_tree(const std::string& ifile)

/* @brief	Read skeletons from files then counting kmers in reads and
 * generate <kmer, cnt> data structure where cnt > 1
 *
 * [taxaKCnt]: vector stores number of kmers for each taxon; this vector
 * will be used	to calculate per batch how many skeleton files to read in
 * for omp processing
 */
void ProcDB::count_sampled_kmers_in_reads(const ivec_t& taxaKCnt,
		const strvec_t& ipfqs, const strvec_t& isfqs) {

	if (!silent_) std::cout << "\nCount sampled kmers in reads ...\n";

	// partition files to be processed in parallel
	int sz_per_batch = 500000000; // 500 Mil entries of kmers ~ 4 Gbyte

	int num_taxa = get_num_taxa_();
	if (taxaKCnt.size() != num_taxa)
		abording ("ProcDB::count_sampled_kmers_in_reads SC failed");
	ivec_t partition_idx;
	int cnt = 0;// kmer counts for visited files
	for (int i = 0; i < num_taxa; ++ i) {
		cnt += taxaKCnt[i];
		if (cnt > sz_per_batch) {
			partition_idx.push_back (i);
			cnt = 0;
		}
	}
	partition_idx.push_back(num_taxa - 1);

	if (!silent_) std::cout << "\t" << partition_idx.size() << " iterations, "
			<< sz_per_batch/1000000 << " mil kmers per iter\n";

	int start = 0;
	for (int p = 0; p < partition_idx.size(); ++ p) {

		if (!silent_) std::cout << "\n\titeration: " << p + 1 << "\n";

		int end = partition_idx[p];

		u64vec_t kmers;

		// process current batch of files in parallel and store record [kmers]
		#pragma omp parallel
		{
			u64vec_t local;
			#pragma omp for schedule (dynamic)
			for (int i = start; i <= end; ++ i) { // process current batch of files
				read_skeletons_(local, skeletonFileList_[i]);
			}
			#pragma omp critical
			{
				kmers.insert(kmers.end(), local.begin(), local.end());
			}
		}

		//substract kmers that already been accounted for in [kmerInfo_]
		substract_counted_kmers_(kmers);

		CountKmers kmercounter (kmers.begin(), kmers.end(), k_);
		process_fq_fl(kmercounter, kmerInfo_, batch_, ipfqs);
		process_fq_fl(kmercounter, kmerInfo_, batch_, isfqs, false);

//		count_kmers_in_read_files_(kmers, ipfqs);
//		count_kmers_in_read_files_(kmers, isfqs, false);

		start = end + 1;
	}

	if (!silent_) std::cout << "\t# kmers present in reads: "
					<< get_kmerinfo_sz_ () << "\n";

	if (!silent_) std::cout << "\t# after remove low freq (= 1) kmers: ";
	rm_lowfreq_kmers_ (1);
	if (!silent_) std::cout << get_kmerinfo_sz_ () << "\n";


//	debug_gen_kmercnt_histogram(kmerInfo_.begin(), kmerInfo_.end());

} // ProcDB::count_sampled_kmer_in_reads

/* @brief Read skeletons into a 1-d kmer vector
 */
void ProcDB::read_skeletons_(u64vec_t& kmers, const std::string& ifile) {
	std::ifstream ifh;
	xny::openfile (ifh, ifile, std::ios::in);
	int num_rows;
	ifh >> num_rows;
	for (int i = 0; i < num_rows; ++ i) {
		int gi, num_elem;
		ifh >> gi >> num_elem;
		uint64_t elem;
		for (int j = 0; j < num_elem; ++ j) {
			ifh >> elem;
			kmers.push_back(elem);
		}
	}
	xny::closefile(ifh);
} // ProcDB::read_skeletons_(u64vec_t& kmers, const std::string& ifile)

/* @brief Read skeletons to be 2-d vectors
 */
void ProcDB::read_skeletons_(ivec_t& GIs, uu64vec_t& skeletons,
		const std::string& ifile) {
	std::ifstream ifh;
	xny::openfile (ifh, ifile, std::ios::in);
	int num_rows;
	ifh >> num_rows;
	skeletons.resize(num_rows);
	GIs.resize(num_rows);
	for (int i = 0; i < num_rows; ++ i) {
		int gi;
		int num_elem;
		ifh >> gi >> num_elem;
		uint64_t elem;
		u64vec_t skeleton (num_elem);
		for (int j = 0; j < num_elem; ++ j) {
			ifh >> elem;
			skeleton[j] = elem;
		}
		skeletons[i] = skeleton;
		GIs[i] = gi;
	}
	xny::closefile(ifh);
} // ProcDB::read_skeletons_(uu64vec_t& skeletons, const std::string& ifile)

/* @brief Remove from [kmers] elements that already exist in [kmerInfo_]
 */
void ProcDB::substract_counted_kmers_(u64vec_t& kmers) {

	rm_duplicated_elements (kmers);

	bvec_t is_found (kmers.size(), false);
	int idx = 0;
 	for (auto& x: kmers) {
		if (std::binary_search(kmerInfo_.begin(), kmerInfo_.end(),
				std::make_tuple(x, 0, false), cmp_u64ibtuple())) {
			is_found[idx] = true;
		}
		++ idx;
	}
	rm_marked_elements (kmers, is_found);
}//

void ProcDB::rm_lowfreq_kmers_(int freq) {
	bvec_t to_rm (kmerInfo_.size(), false);
	int idx = 0;
	for (auto& x: kmerInfo_) {
		if (std::get<1> (x) <= freq) to_rm[idx] = true;
		++ idx;
	}
	rm_marked_elements (kmerInfo_, to_rm);
}

/* @brief	Generate <k1, k2> pairs that are adjacent on some contigs
 * 			Set [kmerTaxaID_]
 */
void ProcDB::generate_locality_kpairs_per_taxa(
		std::vector<uu64pair_t>& taxaKPairs, int w) {

	if (!silent_) std::cout << "\nGenerate locality kmer pairs for taxa...\n";
	// process each skeleton files in parallel
	int num_taxa = get_num_taxa_();
	int min_block_sz = 2;
	int gap = 10;

	std::vector<uu64ituple_t> kpairtaxa;

	int debug_block_sz = 0;

	#pragma omp parallel
	{
		std::vector<uu64ituple_t> local_kpairtaxa;
		std::vector<u64ipair_t> local_ktaxa;

		int debug_local_block_sz = 0;

		#pragma omp for schedule (dynamic)
		for (int i = 0; i < num_taxa; ++ i) {
			uu64vec_t skeletons;
			ivec_t GIs;
			read_skeletons_(GIs, skeletons, skeletonFileList_[i]);

			for (auto& s: skeletons) {
				uu64vec_t blocks;
				identify_valid_skeleton_blocks_(blocks, s, min_block_sz, gap);
				{ // debug
					for (auto& x: blocks) {
						debug_local_block_sz += x.size();
					}
				}
				proc_blocks_(local_ktaxa, local_kpairtaxa, blocks, w, i);
			}
		}
		#pragma omp critical
		{
			debug_block_sz += debug_local_block_sz;

			kpairtaxa.insert(kpairtaxa.end(), local_kpairtaxa.begin(),
					local_kpairtaxa.end());

			add_to_kmertaxaid_(local_ktaxa);
		}
	}// 	#pragma omp parallel

	std::cout << "[debug] block sz: " << debug_block_sz << "\n";

	std::cout << "[debug] kpairtaxa sz: " << kpairtaxa.size() << "\n";

	proc_kpairtaxa_(taxaKPairs, kpairtaxa);
	kpairtaxa.clear();
	rm_dupl_kmertaxaid_();

	if (!silent_) {
		std::cout << "\t# kmer pairs: " << taxaKPairs.size() << "\n";
		std::cout << "\t[kmerTaxa_] size: " << get_kmertaxaid_sz_ () << "\n";
		std::cout << "\t" << get_flagged_kmer_count_ () << " kmers have been flagged\n";
	}

	//debug_gen_kmertaxa_histogram(kmerTaxaID_.begin(), kmerTaxaID_.end());
} // ProcDB::generate_locality_kmer_pairs_per_taxa

/* @brief 	Diving a skeleton into blocks, where each block
 *  			a) is separated by adjacent blocks by distance >= [gap]
 *  			b) the distance between adjacent elements in the block < [gap]
 *  			c) the num of elements in the block is >= min_block_sz
 *
 *			// process skeleton to remove --- (x -- x) ---- structure
 */
void ProcDB::identify_valid_skeleton_blocks_(uu64vec_t& blocks,
		const u64vec_t& skeleton, int min_block_sz, int gap) {

	//min_block_sz = 1; // debug value
	//gap = 2; // debug value

	// (-------)[X-X]{-------}
	// [] : block
	// (--) :prefix
	// {--} : suffix
	int cur_block_sz = 0;
	int prefix_gap_sz = 0;
	int suffix_gap_sz = 0;

	u64vec_t block;
	for (int i = 0; i < skeleton.size(); ++ i) {
		//kmer in skeleton found
		if (std::binary_search(kmerInfo_.begin(), kmerInfo_.end(),
			std::make_tuple(skeleton[i], 0, false), cmp_u64ibtuple())) {
			block.push_back(skeleton[i]);
			++ cur_block_sz;
			if (cur_block_sz > min_block_sz) prefix_gap_sz = 0;
			else suffix_gap_sz = 0;
		} else {
			if (prefix_gap_sz >= gap) {
				if (!cur_block_sz) ++ prefix_gap_sz;
				else {
					++ suffix_gap_sz;
					if (suffix_gap_sz >= gap) { // found block
						block.clear();
						prefix_gap_sz = suffix_gap_sz;
						suffix_gap_sz = 0;
						cur_block_sz = 0;
					}
				}
			} else {
				if (!cur_block_sz) {
					++ prefix_gap_sz;
					if (prefix_gap_sz >= gap) {
						if (block.size()) {
							blocks.push_back(block);
							block.clear();
							cur_block_sz = 0;
						}
					}
				}
				else {
					cur_block_sz = 0;
					prefix_gap_sz = 1;
				}
			}
		}
	}
	if (block.size()) blocks.push_back(block);

	/*{ // debug print
		std::cout << "[debug]\n";
		for (auto& x: blocks) {
			for (auto& y: x) std::cout << "\t" << y;
			std::cout << "\n";
		}
	}*/
} // ProcDB::identify_valid_skeleton_blocks_

/* @brief	Generate <k1, k2, taxaID> and <k, taxaID> from blocks
 */
void ProcDB::proc_blocks_(std::vector<u64ipair_t>& ktaxa,
		std::vector<uu64ituple_t>& kpairtaxa,
		const uu64vec_t& blocks,
		int w,
		int taxaID) {

	for (auto& block: blocks) {
		u64set_t elements;
		/*iset_t test (block.begin(), block.end());
		if (test.size() < 2) {
			std::cout << "found unique elem\n";
			exit(1);
		}*/
		int sz = block.size();
		for (int i = 0; i < sz; ++ i) {
			for (int j = 0; j < w; ++ j) {
				if (i + j + 1 < sz) {
					uint64_t smaller = block[i], larger = block[i + j + 1];
					if (smaller > larger) {
						smaller = larger;
						larger = block[i];
					}
					if (smaller != larger) {
						kpairtaxa.push_back(
								std::make_tuple(smaller, larger, taxaID));
						elements.insert(smaller);
						elements.insert(larger);
					}
				} else break;
			}
		}
		for (auto& x: elements) ktaxa.push_back(u64ipair_t (x, taxaID));
	}
} //ProcDB::proc_blocks_

/* @brief	For kmer pairs that belong to multiple taxa, set flag field
 *  in [kmerInfo_] to true. Convert kpairtaxa to taxakpairs structure by
 *  remove duplicated kpairs
 */
void ProcDB::proc_kpairtaxa_(std::vector<uu64pair_t>& taxakpairs,
		std::vector<uu64ituple_t>& kpairtaxa) {

	int sz = kpairtaxa.size();
	if (!sz) return;

	std::sort (kpairtaxa.begin(), kpairtaxa.end(), cmp_uu64ituple ());

	uu64pair_t elem (std::get<0>(kpairtaxa.front()), std::get<1>(kpairtaxa.front()));
	iset_t taxaIDs = { std::get<2>(kpairtaxa.front()) };
	taxakpairs.push_back(elem);
	for (int i = 1; i < sz; ++ i) {
		uu64pair_t cur_elem (std::get<0>(kpairtaxa[i]), std::get<1>(kpairtaxa[i]));
		if (cur_elem == elem) taxaIDs.insert(std::get<2> (kpairtaxa[i]));
		else {
			if (taxaIDs.size() > 1) { // set kmerInfo_
				if (!set_kmerinfo_flag_ (elem.first) ||
					!set_kmerinfo_flag_ (elem.second) ) {
					abording ("ProcDB::proc_taxakpairs_: SC failed");
				}
			}
			elem = cur_elem;
			taxaIDs = iset_t { std::get<2> (kpairtaxa[i]) };
			taxakpairs.push_back(elem);
		}
	}
	// process the last unprocessed element
	if (taxaIDs.size() > 1) { // set kmerInfo_
		if (!set_kmerinfo_flag_ (elem.first) ||
			!set_kmerinfo_flag_ (elem.second) ) {
			abording ("ProcDB::proc_taxakpairs_: SC failed");
		}
	}

} // ProcDB::proc_kpairtaxa_

/* @brief	Scan read data, count <k1, k2> pairs in any read/frag
 * if k1 and k2 belong to different taxa. Output in [rKPairs] kmer pairs
 * with cnt >= min_freq
 */
void ProcDB::generate_locality_kpairs_in_reads(
		std::vector<uu64pair_t>& rKPairs,
		int min_freq,
		const strvec_t& ipfqs,
		const strvec_t& isfqs) {

	if (!silent_) std::cout << "\nGenerate locality kmer pairs in reads\n";

	std::vector<uu64ituple_t> rkpcnt;

	ivec_t debug_kfragfreq;
	std::back_insert_iterator<ivec_t> it_debug_kfragfreq (debug_kfragfreq);
	CountKP kpcounter (kmerInfo_.begin(), kmerInfo_.end(), kmerTaxaID_.begin(),
			kmerTaxaID_.end(), it_debug_kfragfreq, k_, min_freq);
	process_fq_fl(kpcounter, rkpcnt, batch_, ipfqs);
	//process_fq_fl(kpcounter, rkpcnt, batch_, ipfqs); // treat as non-pairs
	process_fq_fl(kpcounter, rkpcnt, batch_, isfqs, false);

	//count_kpair_in_read_files_ (rkpcnt, min_freq, ipfqs);
	//count_kpair_in_read_files_ (rkpcnt, min_freq, isfqs, false);

	if (!silent_) std::cout << "\t# kmer pairs: " << rkpcnt.size() << "\n";
	for (auto& x: rkpcnt) {
		if (std::get<2> (x) >= min_freq)
			rKPairs.push_back(uu64pair_t (std::get<0> (x), std::get<1> (x)));
	}
	if (!silent_) {
		std::cout << "\t# after rm pairs with freq < " << min_freq
				<< ": " << rKPairs.size() << "\n";
	}

	/*std::cout << "\n[ --- DEBUG \n\n";
	std::cout << "Histogram (col1) #kmers (col2) this many fragments each contain\n";
	std::cout << "#kmers\t\t#fragments\n-----------------------------------\n";
	std::sort(debug_kfragfreq.begin(), debug_kfragfreq.end());
	output_histogram(debug_kfragfreq);
	std::cout << "----------------------------------------------------\n";
	std::cout << "\n DEBUG ---] \n\n";
	*/

} // ProcDB::generate_locality_kpairs_in_reads

/* @brief	Identify connectivity of kmer pairs in [taxaKPairs] relying
 * on [rKPairs] and set flag field in [kmerInfo_]
 */
void ProcDB::set_kmer_multiplicity(const std::vector<uu64pair_t>& taxaKPairs,
 			const std::vector<uu64pair_t>& rKPairs) {

	if (!silent_) std::cout << "\nSet kmer multiplicity...\n";
	int num_tkp = taxaKPairs.size();

	#pragma omp parallel for
	for (int i = 0; i < num_tkp; ++ i) {
		uint64_t kl = taxaKPairs[i].first, kr = taxaKPairs[i].second;
		u64set_t nl = get_neighborset_in_pairs(kl, rKPairs, cmp_uu64pair());
		u64set_t nr = get_neighborset_in_pairs(kr, rKPairs, cmp_uu64pair());
		nl.erase (kr);
		nr.erase (kl);
		for (auto& l: nl) {
			for (auto& r: nr) {
				uu64pair_t target (std::min(l, r), std::max(l, r));
				if (std::binary_search(taxaKPairs.begin(), taxaKPairs.end(),
						target)) {
					// set flag
					set_kmerinfo_flag_(l);
					set_kmerinfo_flag_(r);
					set_kmerinfo_flag_(kl);
					set_kmerinfo_flag_(kr);
				}
			}
		}
	}

	if (!silent_) std::cout << "\t" << get_flagged_kmer_count_()
			<< " kmers have been flagged\n";

} // ProcDB::set_kmer_multiplicity

/* @brief	Identify taxa with
 * 			flagged_kmer_cnt/total_sampled_kmer_cnt >= max_perc
 * 			and output spurious taxa to file
 */
void ProcDB::detect_spurious_taxa(int max_perc,
		const std::string& odir, const std::string& oprefix) {

	if (!silent_) std::cout << "\nDetect spurious taxa...\n";
	// prepare output file
	std::string opath = odir + "/spurious_taxa_";
	opath += oprefix + ".txt";
	std::ofstream ofh;
	xny::openfile(ofh, opath);

	int num_taxa = get_num_taxa_();

	// <num flagged kmer, num total kmer, taxa_idx>
	std::vector<iiituple_t> knum_taxa (num_taxa, iiituple_t {0, 0, -1});
	for (auto& x: kmerTaxaID_) {
		if (x.second >= num_taxa)
			abording ("ProcDB::detect_spurious_taxa SC failed");
		++ std::get<1> (knum_taxa[x.second]);
		std::get<2> (knum_taxa[x.second]) = x.second;
		if (get_kmerinfo_flag_(x.first)) {
			++ std::get<0> (knum_taxa[x.second]);
		}
	}
	std::sort (knum_taxa.begin(), knum_taxa.end(), cmp_iiituple());

	int idx = 0;
	int num_spurious = 0, num_sampled = 0;

	ofh << "Taxa_Index\tTaxaID\tNum_flagged_kmer\tNum_total_kmer\tNames\n";
	ofh << "-----------------------------------------------------------\n";

	iset_t spurious_taxa_index_set;
	for (int i = 0; i < knum_taxa.size(); ++ i) {
		int num_flagged_kmer = std::get<0> (knum_taxa[i]),
			num_total_kmer = std::get<1> (knum_taxa[i]);
		if (num_total_kmer) {
			++ num_sampled;
			int taxa_idx = std::get<2> (knum_taxa[i]);
			int taxaID = taxaIDs_[taxa_idx];

			if (num_flagged_kmer * 100 / num_total_kmer >= max_perc) {
				++ num_spurious;
				spurious_taxa_index_set.insert(taxa_idx);
				ofh << taxa_idx << ": " << taxaID << ": " << num_flagged_kmer
					<< ": " <<  num_total_kmer << " - "
					<< get_taxa_names_(taxaID) << "\n\n";

			}
		}
	} // for (int i = 0

	xny::closefile(ofh);

	if (!silent_) std::cout << "\tTaxa Info: total, num_sampled, num_spurious: "
							<< knum_taxa.size() << ", "
							<< num_sampled << ", " << num_spurious << "\n";

	/* remove kmers that are associated with spurious taxa from [kmerTaxaID_]*/
	bvec_t is_found (get_kmertaxaid_sz_(), false);
	idx = 0;
	for (auto& x: kmerTaxaID_) {
		if (spurious_taxa_index_set.count(x.second)) is_found[idx] = true;
		++ idx;
	}
	rm_marked_elements(kmerTaxaID_, is_found);
	if (!silent_) std::cout << "\t[kmerTaxa_] size: "
			<< get_kmertaxaid_sz_() << "\n";

} //ProcDB::detect_spurious_taxa

std::string ProcDB::get_taxa_names_(int taxaID) {
	std::string names;
	auto lb = std::lower_bound(taxa_names_.begin(), taxa_names_.end(),
			istrvecpair_t (taxaID, strvec_t()), cmp_istrvecpair());
	int dist = std::distance(taxa_names_.begin(), lb);
	if (dist >= taxa_names_.size()) {
		abording ("ProcDB::detect_spurious_taxa SC0 failed");
	} else if (taxa_names_[dist].first != taxaID) {
		std::cout << taxa_names_[dist].first << ", " << taxaID << "\n";
		abording ("ProcDB::detect_spurious_taxa SC1 failed");
	}
	int num_names = taxa_names_[dist].second.size();
	if (num_names) names += taxa_names_[dist].second.front();
	for (int i = 1; i < num_names; ++ i) {
		names += " | " + taxa_names_[dist].second[i];
	}
	//for (auto& x: taxa_names_[dist].second) names += x + " | ";
	return names;
} // ProcDB::get_taxa_names_

std::string ProcDB::get_taxa_rank_(int taxaID) {
	std::string rank;
	auto it = std::lower_bound(tree_.begin(), tree_.end(),
			iistrtuple_t(taxaID, 0, std::string()), cmp_iistrtuple());
	if (it == tree_.end() || std::get<0>(*it) != taxaID) {
		std::cout << "\tNode:" << taxaID << " couldn't be found\n";
		warning ("ProcDB::get_taxa_rank_ SC Failed");
	} else rank = std::get<2>(*it);
	return rank;
}// ProcDB::get_taxa_rank_


//
void ProcDB::assign_frags(const strvec_t& ipfqs,	const strvec_t& isfqs,
		const std::string& odir, const std::string& oprefix) {

	if (!silent_) std::cout << "\nAssign fragments to taxa\n";

	int num_taxa = get_num_taxa_();
	ivec_t taxacnts (num_taxa, 0);
	AssignFrag fassigner (kmerTaxaID_.begin(), kmerTaxaID_.end(), k_);
	process_fq_fl(fassigner, taxacnts, batch_, ipfqs);
	process_fq_fl(fassigner, taxacnts, batch_, isfqs, false);

	for (auto& x: taxacnts) x /= 1000; // adjust back

	// <num flagged kmer, num total kmer> per taxa
	std::vector<ipair_t> knum_taxa (num_taxa, ipair_t{0, 0});
	for (auto& x: kmerTaxaID_) {
		if (x.second >= num_taxa)
			abording ("ProcDB::detect_spurious_taxa SC failed");
		++ knum_taxa[x.second].second;
		if (get_kmerinfo_flag_(x.first)) {
			++ knum_taxa[x.second].first;
		}
	}

	// given x_0 x_1 ... x_n-1 as sorted taxacnts from high to low
	// identify the max I (max_I) such that
	// 1) (\sum_i=0^I x_i) / total_cnt >= 50%   and
	// 2) \sum_i=0^(I-1) x_i / (\sum_i=0^I x_i) <= 90% (if I > 0)
	// the intuition is to identify prevailing dominant taxa and remove
	// those from considering which taxa need to be output
	ivec_t sorted_taxacnts = taxacnts;
	std::sort(sorted_taxacnts.begin(), sorted_taxacnts.end());
	std::reverse(sorted_taxacnts.begin(), sorted_taxacnts.end());

	ivec_t prefix_sum (num_taxa, 0);
	for (int i = 0; i < num_taxa; ++ i) {
		if (i > 0) prefix_sum[i] += prefix_sum[i - 1] + sorted_taxacnts[i];
		else prefix_sum[i] = sorted_taxacnts[i];
	}
	int max_I = -1;
	for (int i = 0; i < num_taxa; ++ i) {
		 // remove 50% of content minimum
		if (100 * prefix_sum[i]/prefix_sum.back() >= 50) {
			if (i > 0) {
				// addition at least additional 10%
				if (100 * prefix_sum[i - 1] / prefix_sum[i] <= 90) {
					max_I = i;
				} else {
					if (max_I == -1) max_I = i;
					//else break;
				}
			} else max_I = i;
		}
	}

	int total_adjusted_cnt = prefix_sum.back();
	if (max_I != -1) total_adjusted_cnt -= prefix_sum[max_I];

	// <fragcnt, num flagged kmer, num total kmer, taxa_idx>
	std::vector<i4tuple_t> taxarank;
	std::vector<i4tuple_t> all_taxarank;
	for (int i = 0; i < num_taxa; ++ i) {
		// only consider taxa that account for >= 1/1000
		// of total excluding the most abundance
		if (taxacnts[i] * 1000 / total_adjusted_cnt > 0) {
			taxarank.push_back(i4tuple_t{taxacnts[i], knum_taxa[i].first,
				knum_taxa[i].second, i});
		}
		if (taxacnts[i] > 0) all_taxarank.push_back(i4tuple_t{taxacnts[i],
			knum_taxa[i].first, knum_taxa[i].second, i});
	}
	std::sort(taxarank.begin(), taxarank.end(), cmp_i4tuple());
	std::sort(all_taxarank.begin(), all_taxarank.end(), cmp_i4tuple());
	// output list
	output_taxa_list_(all_taxarank, odir, oprefix);

	// building tree structure of identified nodes
 	imap_t edges;
	build_tree_(edges, taxarank);
	output_taxa_tree_(edges, taxarank, odir, oprefix);

	// building full tree structure
	edges.clear();
	if (all_taxarank.size() > 100) {
		int num_to_erase = all_taxarank.size() - 100;
		all_taxarank.erase(all_taxarank.begin(),
				all_taxarank.begin() + num_to_erase);
	}
	build_tree_(edges, all_taxarank);
	std::string full_tree_prefix = oprefix + ".100";
	output_taxa_tree_(edges, all_taxarank, odir, full_tree_prefix);

} // ProcDB::generate_locality_kpairs_in_reads

/* @brief Output full list of identified taxa that has scores > 0
 */
void ProcDB::output_taxa_list_(const std::vector<i4tuple_t>& taxarank,
		const std::string& odir, const std::string& oprefix) {
	// prepare output file
	std::string opath = odir + "/taxa_rank_";
	opath += oprefix + ".txt";
	std::ofstream ofh;
	xny::openfile(ofh, opath.c_str());
	ofh << "Taxa_Index\tTaxaID\tScore\tNum_flagged_kmer\tNum_total_kmer\tNames\n";
	ofh << "-----------------------------------------------------------\n";

	for (int i = 0; i < taxarank.size(); ++ i) {
		int score = std::get<0> (taxarank[i]),
			num_flagged_kmer = std::get<1> (taxarank[i]),
			num_total_kmer = std::get<2> (taxarank[i]),
			taxa_idx = std::get<3> (taxarank[i]);

		int taxaID = taxaIDs_[taxa_idx];

		ofh << taxa_idx << ": " << taxaID << ": " << score << ": "
				<< num_flagged_kmer << ": " <<  num_total_kmer << " - "
				<< get_taxa_names_(taxaID) << "\n\n";
	} // for (int i = 0
	xny::closefile(ofh);
	if (!silent_) std::cout << "\toutput " << taxarank.size()
			<< " ranked list in " << opath << "\n";
} // output_taxa_list_

/* @brief Build taxonomic tree
 */
void ProcDB::build_tree_(imap_t& edges, const std::vector<i4tuple_t>& taxarank){
	if (!silent_) std::cout << "\tconstruct taxonomy tree...";
	// identify path to root for each taxID
	int sz = taxarank.size();
	iset_t taxaIDs, multi_taxaIDs; //non-leaf taxaIDs

	// build all edges using tree_
	for (int i = 0; i < sz; ++ i) {
		int taxa_idx = std::get<3> (taxarank[i]);
		int taxaID = taxaIDs_[taxa_idx];
		while (true) {
			// look for parent ID
			auto it = std::lower_bound(tree_.begin(), tree_.end(),
					iistrtuple_t(taxaID, 0, std::string()), cmp_iistrtuple());
			if (it == tree_.end() || std::get<0>(*it) != taxaID) { // parent ID not found
				std::cout << "\tNode:" << taxaID << " couldn't be found\n";
				warning ("ProcDB::build_tree_ SC Failed");
				break;
			} else { // parent found
				int parent_taxaID = std::get<1>(*it);
				if (taxaIDs.count(parent_taxaID))
					multi_taxaIDs.insert(parent_taxaID);
				else taxaIDs.insert(parent_taxaID);
				if (taxaID != parent_taxaID)	{
					edges[taxaID] = parent_taxaID;
					taxaID = parent_taxaID;
				} else break;
			}
		}
	}

	// skip intermediate nodes that contain only single child
	for (auto it = edges.begin(); it != edges.end(); ++ it) {
		int parent = it->second;
		while (! multi_taxaIDs.count(parent)) {
			auto it_e = edges.find(parent);
			if (it_e != edges.end()) parent = it_e->second;
			else break;
		}
		it->second = parent;
	}

	// remove obsolete edges
	for (auto it = taxaIDs.begin(); it != taxaIDs.end(); ++ it){
		if (!multi_taxaIDs.count(*it)) edges.erase(*it);
	}

	if (!silent_) std::cout << " complete\n";

} // build_tree_

/* @brief Output taxonomic tree in Graphviz (dot readable) format
 */
void ProcDB::output_taxa_tree_(const imap_t& edges,
		const std::vector<i4tuple_t>& taxarank,
		const std::string& odir, const std::string& oprefix) {

	// ---- calculating score for each node ---
	imap_t node_score;
	accumulate_scores_(node_score, edges, taxarank);

	// --- identify taxonomy rank and name for each node ---

	// prepare output file
	std::string opath = odir + "/taxa_tree_";
	opath += oprefix + ".txt";
	std::ofstream ofh;
	xny::openfile(ofh, opath.c_str());

	ofh << "digraph G {\n";
	for (auto it = edges.begin(); it != edges.end(); ++ it) {

		// --- identify and mark any edge that carry > 30% of the score
		//under the same node
		auto it_parent_score = node_score.find(it->second),
			 it_kid_score = node_score.find(it->first);
		std::string edge_property = "[style=dotted];";
		if (it_parent_score != node_score.end()
				&& it_kid_score != node_score.end()) {
			if (it_kid_score->second * 100 / it_parent_score->second > 30) {
				edge_property = "[color=red, style=bold];";
			}
		}

		ofh << "\t" << it->second << " -> " << it->first << edge_property << "\n";
	}
	for (auto it = node_score.begin(); it != node_score.end(); ++ it) {
		ofh << "\t" << it->first << "[label=\"" << it->second << "\\n"
			<< get_taxa_rank_(it->first) << "(" << it->first << ")\\n"
			<< get_taxa_names_(it->first)
			<< "\"];" << std::endl;
	}
	ofh << "}\n";
	xny::closefile(ofh);

	if (!silent_) std::cout << "\ttree has been written to " << opath << "\n";
} // ProcDB::output_taxa_tree_

/* @brief	Propagate scores through the tree
 */
void ProcDB::accumulate_scores_(imap_t& node_score, const imap_t& edges,
		const std::vector<i4tuple_t>& taxarank) {
	for (auto it = edges.begin(); it != edges.end(); ++ it) { // init
		node_score[it->first] = 0;
		node_score[it->second] = 0;
	}
	for (int i = 0; i < (int) taxarank.size(); ++ i) {
		int taxa_idx = std::get<3> (taxarank[i]),
			score = std::get<0> (taxarank[i]);
		int taxaID = taxaIDs_[taxa_idx];
		auto it = node_score.find(taxaID);
		if (it != node_score.end()) it->second += score;
		//propagate the score through the tree
		while (true) {
			auto it_e = edges.find(taxaID); //
			if (it_e != edges.end()) {
				int parentID = it_e->second;
				it = node_score.find(parentID);
				if (it != node_score.end()) it->second += score;
				taxaID = parentID;
			} else break;
		}
	}
} // ProcDB::accumulate_scores_


/* ---------------------------  DEBUG ---------------------------------
 *
 */
/* @brief	Given a vector of taxaIDs, print out skeleton for each of the
 * taxa
 */
void ProcDB::debug_print_skeleton(const std::string& opath,
		const ivec_t& taxaIDs) {
	std::ofstream ofh;
	xny::openfile(ofh, opath);
	for (auto& taxaID: taxaIDs) {
 		// identify index
		int index = 0;
		for (; index < taxaIDs_.size(); ++ index) {
			if (taxaIDs_[index] == taxaID) break;
		}
		if (index < taxaIDs_.size()) {
			std::string skeletonfile = skeletonFileList_[index];
			uu64vec_t skeletons;
			ivec_t GIs;
			read_skeletons_(GIs, skeletons, skeletonfile);

			// for each GI identify its order
			ivec_t gi_ranks;
			if (GIs.size()){
				int prev_gi = GIs.front();
				int rank = 0;
				gi_ranks.push_back(rank);
				for (int i = 1; i < GIs.size(); ++ i){
					if (GIs[i] == prev_gi){
						rank ++;
						gi_ranks.push_back(rank);
					} else {
						rank = 0;
						gi_ranks.push_back(rank);
						prev_gi = GIs[i];
					}
				}
			}

			int idx = 0;
			ofh << "[[[ " << taxaID << "]]]\n";

			strvec_t marked_skeletons;
			ivec_t marked_gis;
			ivec_t marked_ranks;
			for (auto& skeleton: skeletons) {
				std::string marked_skeleton;
				bool has_content = false;
				std::vector<iset_t> tids;
				for (auto& elem: skeleton) {
					if (is_kmer_in_taxa_(elem, index)) {
						if (get_kmerinfo_flag_(elem)) {
							tids.push_back(get_kmertaxaid_tids_(elem));
							marked_skeleton += "m";
							marked_skeleton += std::to_string (tids.back().size()) + ",";
							marked_skeleton += std::to_string (get_kmerinfo_freq_(elem));
							has_content = true;
						} else {
							tids.push_back(get_kmertaxaid_tids_(elem));
							marked_skeleton += "u";
							marked_skeleton += std::to_string (tids.back().size()) + ",";
							marked_skeleton += std::to_string (get_kmerinfo_freq_(elem));
							has_content = true;
						}
					} else marked_skeleton += "-"; // ofh << "-";
				}
				if (has_content) {
					marked_skeletons.push_back(marked_skeleton);
					marked_gis.push_back(GIs[idx]);
					marked_ranks.push_back(gi_ranks[idx]);
				}

				// print out tids
				/*for (auto& x: tids) {
					for (auto& y: x) ofh << y << "\t";
					ofh << "\n";
				}*/

				++ idx;
			}
			// print out to file

			for (int i = 0; i < marked_skeletons.size(); ++ i) {
				ofh << "GI = " << marked_gis[i] << "." << marked_ranks[i]
				    << "\n" << marked_skeletons[i] << "\n";
			}
		}
	}
} // ProcDB::debug_print_skeleton
/* @brief	Construct taxonomy tree
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// obsolete functions
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////


/* @brief	Set up file list to store skeleton for each taxa
 */
/* @brief	Analyze context information for contig per taxa and remove
 * kmers that are likely spuriously support for this taxa even though
 * they occur in reads
 *
 * Tasks: 1) For each taxa (walk through contig via info of [kmerStat_])
 * 			a) identify and clear -----Xn----- pattern
 * 			(n is the user specified para, min number of Xs to be observed)
 * 			b) if (Ms/Xs > 95%) mark taxa as spurious
 * 		    c) record kmer TaxaID pairs to be removed
 * 		  2) update kmerTaxaID_
 */
/*
void ProcDB::rm_spurious_support_per_taxa () {

	if (!silent_) std::cout << "Using context info to remove spurious taxa support ...\n";

	int num_taxa = get_num_taxa_();
	if (!num_taxa) return;

	if (!silent_) std::cout << "\t" << num_taxa << " taxa to process\n";

	std::vector<u64ipair_t> kt_to_del;

	#pragma omp parallel
	{
		std::vector<u64ipair_t> local_kt_to_del;

		#pragma omp for schedule (dynamic)
		for (int i = 0; i < num_taxa; ++ i) { // for each taxa

			std::vector<u64bpair_t> valid_kmers; // valid kmers, accrued for all contigs
			u64set_t kmer_to_del; 				// kmers to del, accured for all contigs

			int fnum = fileList_[i].size();
			for (int j = 0; j < fnum; ++ j) {
				std::ifstream fh (fileList_[i][j].c_str());
				if (!fh.good())	abording ("ProcDB::rm_spurious_support_per_taxa: "
							"can't open " + fileList_[i][j]);

				std::string contig;
				std::string header;
				std::string line;
				u64vec_t skeleton;

				while (std::getline (fh, line)) {
					if (line.empty()) continue;
					if (line.at(0) == '>') {
						skeleton.clear();
						add_kmer_from_contig (skeleton, contig);
						//if (i == 6761) {
						//	std::cout << header << "\n";
						//}
						header = line;
						proc_contig_context (valid_kmers, kmer_to_del, skeleton);
						contig.clear();
					} else contig += line;
				}
				skeleton.clear();
				add_kmer_from_contig (skeleton, contig);

				//if (i == 6761) {
				//	std::cout << header << "\n";
				//}
				proc_contig_context (valid_kmers, kmer_to_del, skeleton);

				fh.close();

			} // for (int j = 0

			// check for current taxa kmers to be removed from [kmerTaxaID_]
			if (valid_kmers.size()) {
				confirm_deletion (valid_kmers, kmer_to_del);
				if (!valid_kmers.size()) {
					std::cout << "Spurious taxa: " << i << ", " << taxaIDs_[i] << "\n";
				}
				for (auto& x: kmer_to_del) {
					local_kt_to_del.push_back(u64ipair_t (x, i));
				}
			}

		} // for (int i = 0

		#pragma omp critical
		{
			kt_to_del.insert(kt_to_del.end(), local_kt_to_del.begin(),
					local_kt_to_del.end());
		}

	} // 	#pragma omp parallel

	std::cout << "[debug]: total kt_to_del: " << kt_to_del.size() << "\n";
} //void ProcDB::rm_spurious_support_per_taxa ()
*/
/*
void ProcDB::add_kmer_from_contig (u64vec_t& kmers, const std::string& contig) {
	xny::low_complexity lc (30, 50, 80);
	int len = contig.length();
	int num = (len - k_ + 1) * perc_/100;
	if (num <= 0) return;
	int step = 100 / perc_; // equal to (len - k_ + 1) / num;
	for (int i = 0; i < num; ++ i) { //first char a or c, use as is, otherwise rvc
		uint64_t val;
		std::string kmer = get_kstr (contig.substr(i*step, k_));
		if (xny::str2ID<uint64_t>(val, kmer) && !lc (kmer)) kmers.push_back(val);
	}
} //ProcDB::add_kmer_from_contig
*/
/* @brief	Measure the percentage of Ms/Xs and kmers to be removed
 * for current taxa
 */
/*
void ProcDB::confirm_deletion (std::vector<u64bpair_t>& valid_kmers,
	u64set_t& kmer_to_del) {

	if (!valid_kmers.size()) return;
	int num_m = 0;
	for (auto& x: valid_kmers) if (x.second) ++ num_m;

	if (100 * num_m/valid_kmers.size() >= 95) { // the taxa is spurious
		for (auto& x: valid_kmers) kmer_to_del.insert(x.first);

		std::cout << num_m << ", " << valid_kmers.size() << "\n";

		valid_kmers.clear();
	} else {
		std::sort (valid_kmers.begin(), valid_kmers.end());
		u64set_t delset;
		for (auto& x: kmer_to_del) {
			if (!std::binary_search(valid_kmers.begin(), valid_kmers.end(),
				u64bpair_t (x, false), cmp_u64bpair())) delset.insert(x);
		}
		kmer_to_del = delset;
	}

} // ProcDB::confirm_deletion
*/
/*

void ProcDB::proc_contig_context (std::vector<u64bpair_t>& valid_kmers,
	u64set_t& kmer_to_del, const u64vec_t& skeleton) {
	int min_block_sz = 2;
	int step = 10;

	//min_block_sz = 2; // debug value
	//step = 3; // debug value

	// (-------)[X-X]{-------}
	// [] : block
	// (--) :prefix
	// {--} : suffix
	int cur_block_sz = 0;
	int prefix_gap_sz = 0;
	int suffix_gap_sz = 0;

	for (int i = 0; i < skeleton.size(); ++ i) {
		auto lb = std::lower_bound(kmerStat_.begin(), kmerStat_.end(),
				u64ipair_t (skeleton[i], 0), cmp_u64bpair());
		int dist = std::distance(kmerStat_.begin(), lb);
		if (lb->first == skeleton[i]) { // found
			valid_kmers.push_back(*lb);
			++ cur_block_sz;
			if (cur_block_sz > min_block_sz) prefix_gap_sz = 0;
			else suffix_gap_sz = 0;
		} else {
			if (prefix_gap_sz >= step) {
				if (!cur_block_sz) ++ prefix_gap_sz;
				else {
					++ suffix_gap_sz;
					if (suffix_gap_sz >= step) { // found block

						// remove store last cur_block_sz element
						for (int cnt = 0; cnt < cur_block_sz; ++ cnt) {
							kmer_to_del.insert(valid_kmers.back().first);
							valid_kmers.pop_back();
						}

						prefix_gap_sz = suffix_gap_sz;
						suffix_gap_sz = 0;
						cur_block_sz = 0;
					}
				}
			} else {
				if (!cur_block_sz) ++ prefix_gap_sz;
				else {
					cur_block_sz = 0;
					prefix_gap_sz = 1;
				}
			}
		}
	}

} // ProcDB::proc_contig_context
*/

 /* @brief	Get vec (nodeID, name lists) from input file */
/*
void get_node_names (std::vector<istrvecpair_t>& node_names, const std::string& ifile) {
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);
	int cnt = 0;

 	strvec_t rows;
	read_file_chunk (rows, ifh, INT_MAX);

	int prev_node;
	strvec_t names;
	for (int i = 0; i < rows.size(); ++ i) {
		strvec_t entries;
		split('|', rows[i], std::back_inserter(entries));
		std::string name = entries[1];
		if (name.length() > 2) name = name.substr(1, name.length() - 2);

		if (i == 0) { // first entry
			prev_node = std::atoi(entries[0].c_str());
			names.push_back(name);
		} else {
			int node = std::atoi(entries[0].c_str());
			if (node != prev_node) {
				node_names.push_back(istrvecpair_t (prev_node, names));
				prev_node = node;
				names = strvec_t (1, name);
			} else names.push_back(name);
		}
	}
	if (names.size()) node_names.push_back(istrvecpair_t (prev_node, names));

	std::cout << node_names.size() << " found\n";
	ifh.close();

} // get_node_names


// replaced by ProcDB::set_db_filelist
void get_input_db_files (std::vector<strvec_t>& input_files,
		ivec_t& nodes, const std::string& idir, const strvec_t& ifl) {

	if (!idir.empty()) {
		DIR* dir = opendir (idir.c_str());
		if (!dir) abording ("cannot open " + idir);
		struct dirent *dirp;

		std::vector<ipair_t> vec_gi_len;
		std::map<std::string, int> filePrefix_fileIdx;

		while ((dirp = readdir(dir))) {
			std::string filename = dirp->d_name;
			unsigned pos = filename.find_last_of('.');
			std::string suffix = filename.substr (pos + 1);
			std::string prefix = xny::get_prefix (filename, true, ".");
			std::string path = idir + "/" + filename;
			if (suffix.compare("fa") == 0 || suffix.compare("fasta") == 0) {
				auto it = filePrefix_fileIdx.find(prefix);
				if (it != filePrefix_fileIdx.end()) {
					input_files[it->second].push_back(path);
				} else {
					input_files.push_back(strvec_t(1, path));
					filePrefix_fileIdx[prefix] = input_files.size() - 1;
					nodes.push_back(std::atoi(prefix.c_str()));
				}
			}
		}
		closedir (dir);
	}

	for (auto& f: ifl) {
		std::string suffix = xny::get_suffix (f, true, "/");
		nodes.push_back(std::atoi(xny::get_prefix (suffix, true, ".").c_str()));
		input_files.push_back (strvec_t (1, f));
	}
} // get_input_db_files
*/

/* @brief Given [kmers] which are kmer sampled from a contig, identify
 * which kmers exists in the reads (i.e. [record]), update 'u'nique and
 * 'm'ultiple feature of the kmer as well as remove from record 'm'ulti-kmers
 * that are likely sampled via chance.
 */
/*
void proc_contig_skeleton (bvec_t& is_record_m, bvec_t& is_to_del,
		std::vector<u64ipair_t>& record, const u64vec_t& kmers, int id) {
	int sz = kmers.size();
	std::vector<char> states (sz);
	for (int i = 0; i < sz; ++ i) {
		auto lb = std::lower_bound(record.begin(), record.end(),
				u64ipair_t (kmers[i], 0), cmp_u64ipair());
		int dist = std::distance(record.begin(), lb);
		if (lb->first == kmers[i]) { // found
			if (is_record_m[dist]) states[i] = 'm';
			else states[i] = 'u';
		} else states[i] = '-';
	}

	// -- m/u --   >>> remove from record     [case a]
	// -u m --   >>>  m -> u					 [case b]
	// mm u mm   >>>  u -> m                  [case c]
	int step = 10;
	for (int i = step; i < sz - step; ++ i) {
		if (states[i] != '-') {
			bool is_case_a = true; // satisfy the checking criteria
			for (int j = 1; j <= step; ++ j) {
				if (states[i-j] != '-' || states[i+j] != '-') {
					is_case_a = false;
					break;
				}
			}
			if (is_case_a) {
 				// search in record
				auto lb = std::lower_bound(record.begin(), record.end(),
						u64ipair_t (kmers[i], 0), cmp_u64ipair());
				int dist = std::distance(record.begin(), lb);
				int target = -1;

				for (; dist < record.size(); ++ dist) {
					if (record[dist].first == kmers[i]) { // found the kmer
						if (record[dist].second == id) {
							target = dist;
							break;
						}
					} else break;
				}
				if (target != -1) {
					if (is_case_a) {
						is_to_del[target] = true;
						states[i] = '-';
					}
				}
			}
		}
	}
} // proc_contig_skeleton
*/

/* @brief	Scan each taxa and record <kmer, fileID> information sorted by
 * 			kmer
 *
 * replaced by ProcDB::sampling()
 */
/*
void db_sample (std::vector<u64ipair_t>& record, const double& perc, int k,
		const std::vector<strvec_t>& input_files, bool silent) {

	int num_node = input_files.size();
	if (!num_node) return;

	if (!silent) std::cout << "\t" << num_node << " taxa nodes to process\n";

	#pragma omp parallel shared (input_files)
	{
		std::vector<u64ipair_t> private_record;
		#pragma omp for schedule (dynamic)
		for (int i = 0; i < num_node; ++ i) {
			//std::cout << i << "\n";
			int fnum = input_files[i].size();
			for (int j = 0; j < fnum; ++ j) {
				u64vec_t kmers;
				sample_kmers (kmers, k, perc, input_files[i][j]);
				for (auto& x: kmers) private_record.push_back(u64ipair_t (x, i));
			}
			//private_record.insert(private_record.end(), kmers.begin(), kmers.end());
			//private_range.push_back(ipair_t(kmers.size(), i));
		}

		#pragma omp critical
		{
			record.insert(record.end(), private_record.begin(),	private_record.end());
			//range.insert(range.end(), private_range.begin(), private_range.end());
		}
	}

	//std::cout << "debug done\n";
	//exit(1);
	// sort
	std::sort (record.begin(), record.end());
	//for (auto& x: record) sum += x.size();
	if (!silent) std::cout << "\tinit size of record: " << record.size() << "\n";

	auto it = std::unique_copy (record.begin(), record.end(), record.begin());
	record.resize(std::distance(record.begin(), it));

	if (!silent) std::cout << "\tfinal size of record (rm dupl entries): "
			<< record.size() << "\n";

	//for (int i = 1; i < range.size(); ++ i) range[i].first += range[i - 1].first;

} // db_sample
*/

/* @brief	Remove entries of [kmerTaxaID_] that do not occur in reads
 */
/*
void ProcDB::rm_sampled_kmers_absent_in_reads () {

 	if (!silent_) std::cout << "Remove sampled kmers absent in reads ...\n";

	// linear scan [kmerTaxaID_] in-place replacing non-occurrence entries
	// with entries that occur in reads
	int sz = kmerTaxaID_.size();
	int p0 = 0, p1;

	// set p0 pointing to an entry that is not found in reads
	// and p1 pointing to an entry having kval diff from p0
	while (p0 < sz) {
		if (! std::binary_search(kmerStat_.begin(), kmerStat_.end(),
			u64bpair_t (kmerTaxaID_[p0].first, false), cmp_u64bpair())) {

			for (p1 = p0 + 1; p1 < sz; ++ p1) {
				if (kmerTaxaID_[p1].first != kmerTaxaID_[p0].first) break;
			}

			//p1 = p0 + 1;
			break;
		} else ++ p0;
	}

	if (p0 >= sz) {
		if (!silent_) std::cout << "\tno update\n";
		return;
	}

	while (p1 < sz) {
		uint64_t val = kmerTaxaID_[p1].first;
		if (std::binary_search(kmerStat_.begin(), kmerStat_.end(),
				u64bpair_t (val, false), cmp_u64bpair())) {

			kmerTaxaID_[p0] = kmerTaxaID_[p1];
			++ p0;
			++ p1;

			for (; p1 < sz; ++ p1, ++ p0) {
				if (kmerTaxaID_[p1].first == val) {
					kmerTaxaID_[p0] = kmerTaxaID_[p1];
				} else break;
			}

		} else {
			++ p1;

			for (; p1 < sz; ++ p1) {
				if (kmerTaxaID_[p1].first != val)  break;
			}
		}
	}
	kmerTaxaID_.resize(p0);

	if (!silent_) std::cout << "\tresulting # kmer taxaID pairs: " << p0 << "\n";

} // ProcDB::rm_sampled_kmers_absent_in_reads ()
*/
/*
void db_context (bvec_t& is_record_m, std::vector<u64ipair_t>& record,
		const double& perc, int k, const std::vector<strvec_t>& input_files,
		bool silent) {

	int num_node = input_files.size();
	if (!num_node) return;

	if (!silent) std::cout << "\t" << num_node << " taxa nodes to process\n";

	bvec_t is_to_del (record.size(), false);
	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < num_node; ++ i) {
		//std::cout << "i = " << i << "\n";

		int fnum = input_files[i].size();
		for (int j = 0; j < fnum; ++ j) {

			std::ifstream fh (input_files[i][j].c_str());
			if (!fh.good()) abording ("sample_kmers: can't open " +
					input_files[i][j]);

			//std::cout << input_files[i][j] << "\n";

			std::string contig;
			std::string header;
			std::string line;
			u64vec_t kmers;

			while (std::getline (fh, line)) {
				if (line.empty()) continue;
				if (line.at(0) == '>') {
					kmers.clear();
					sample_contig (kmers, k, perc, contig);
					//if (i == 6761) {
					//	std::cout << header << "\n";
					//}
					header = line;
					proc_contig_skeleton (is_record_m, is_to_del, record,
							kmers, i);
					contig.clear();
				} else contig += line;
			}
			kmers.clear();
			sample_contig (kmers, k, perc, contig);

			//if (i == 6761) {
			//	std::cout << header << "\n";
			//}
			proc_contig_skeleton (is_record_m, is_to_del, record, kmers, i);

			fh.close();

		} // for (int j = 0

	} // for (int i = 0

	for (int i = 0; i < is_to_del.size(); ++ i) {
		if (is_to_del[i]) is_to_del[i] = false;
		else is_to_del[i] = true;
	}
	clean_records (record, is_record_m, is_to_del);

	std::cout << "record size = " << record.size() << "\n";
} // db_context
*/
/* @brief Step-wise sampling kmers for each file
 */
/*
void sample_kmers (u64vec_t& kmers, int k, const double& perc, const std::string& file) {

	//std::cout << file << "\n";
	std::ifstream fh (file.c_str());
	if (!fh.good()) abording ("sample_kmers: can't open " + file);

	std::string contig;
	std::string line;

	while (std::getline (fh, line)) {
		if (line.empty()) continue;
		if (line.at(0) == '>') {
			sample_contig (kmers, k, perc, contig);
			contig.clear();
		} else contig += line;
	}
	sample_contig (kmers, k, perc, contig);
	fh.close();

} //sample_kmers
*/
/*
void sample_contig (u64vec_t& kmers, int k, const double& perc,
		const std::string& contig) {
	int len = contig.length();
	int num = (len - k + 1) * perc/100;
	if (num <= 0) return;
	int step = (len - k + 1) / num;
	for (int i = 0; i < num; ++ i) { //first char a or c, use as is, otherwise rvc
		uint64_t val;
		std::string kmer = get_kstr (contig.substr(i*step, k));
		if (xny::str2ID<uint64_t>(val, kmer)) kmers.push_back(val);
	}
}
*/
/* @brief	Get all taxaIDs that a given kmer present in
 */
/*
void ProcDB::get_taxaIDs (iset_t& taxaIDs, const uint64_t& kval) {
	auto lb = std::lower_bound(kmerTaxaID_.begin(), kmerTaxaID_.end(),
			u64ipair_t (kval, 0), cmp_u64ipair());
	int dist = std::distance(kmerTaxaID_.begin(), lb);
	for (; lb < kmerTaxaID_.end(); ++ lb, ++ dist) {
		if (lb->first == kval) {
			taxaIDs.insert(lb->second);
		} else break;
	}
} //
*/

/*
void ProcDB::link_sampled_kmer_in_frags (const strvec_t& frags) {

	int num_frag = frags.size();

	std::vector<u64bpair_t> kstat; //sampled kmers [kmerTaxaID_] that occur in [frags]
	#pragma omp parallel
	{
		std::vector<u64bpair_t> local_kstat; // local to each thread

		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {
			int k_pos_idx = 0;
			int last_k_pos = frags[fragidx].length() - k_;

			u64set_t frag_kval_set; // identified kmers in the frag
			iset_t taxaIDs; // taxaIDs associated with current fragment
			//bool is_multi = false; // multiple taxaIDs associated with frag
			while (k_pos_idx <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k_));

				if (xny::str2ID<uint64_t>(val, kmer)) {

						iset_t tmpIDs;
						get_taxaIDs (tmpIDs, val);
						if (tmpIDs.size() > 0) {
							taxaIDs.insert(tmpIDs.begin(), tmpIDs.end());
							frag_kval_set.insert(val);

						}

					++ k_pos_idx;
				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos_idx += pos_n + 1;
					} else k_pos_idx += k_;
				}
			} // while (k_pos_idx <= last_k_pos)

			bool is_multi = false;
			if (taxaIDs.size() <= 6 & taxaIDs.size() >= 2) is_multi = true;

			for (auto& val: frag_kval_set) {
				local_kstat.push_back (u64bpair_t (val, is_multi));
			}
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

		#pragma omp critical
		{
			kstat.insert(kstat.end(), local_kstat.begin(), local_kstat.end());
		}
	}

	accrue_kstat (kstat);

} // ProcDB::link_sampled_kmer_in_frags


void ProcDB::accrue_kstat (std::vector<u64bpair_t>& rhs) {

	if (!rhs.size()) return;
	std::sort (rhs.begin(), rhs.end());

	std::vector<u64bpair_t> merged_kstat;

	int idx_kstat = 0; // [kmerStat_]
	int idx_rhs = 1; // [kstat]

	// element in rhs
	u64bpair_t rhsElem = rhs.front();

	while (idx_rhs < rhs.size() && idx_kstat < kmerStat_.size()) {

		if (kmerStat_[idx_kstat].first >= rhsElem.first) {
			for (; idx_rhs < rhs.size(); ++ idx_rhs) {
				if (rhs[idx_rhs].first == rhsElem.first)
					rhsElem.second |=  rhs[idx_rhs].second;
				else break;
			}
			if (kmerStat_[idx_kstat].first == rhsElem.first) {
				rhsElem.second |= kmerStat_[idx_kstat].second;
				merged_kstat.push_back (rhsElem);
				++ idx_kstat;
			} else merged_kstat.push_back (rhsElem);

			rhsElem = rhs[idx_rhs];
			++ idx_rhs;
		} else {
			merged_kstat.push_back (kmerStat_[idx_kstat]);
			++ idx_kstat;
		}
	}

	// adding the remaining array
	if (idx_kstat < kmerInfo_.size()) {
		merged_kstat.insert(merged_kstat.end(), kmerStat_.begin() + idx_kstat,
				kmerStat_.end());
	} else {
		while (idx_rhs < rhs.size()) {
			if (rhs[idx_rhs].first == rhsElem.first)
				rhsElem.second |= rhs[idx_rhs].second;
			else {
				merged_kstat.push_back (rhsElem);
				rhsElem = rhs[idx_rhs];
			}
			++ idx_rhs;
		}
		merged_kstat.push_back (rhsElem);
	}

	set_kstat (merged_kstat);

} // void ProcDB::accrue_kstat
*/
/*
void ProcDB::count_sampled_kmer_in_frags (const strvec_t& frags) {

	int num_frag = frags.size();

	u64vec_t kval_vec; // identified kmers of [kmerTaxaID_] in [frags]
	#pragma omp parallel
	{
		u64vec_t local_kval_vec; // local to each thread
		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {
			int k_pos_idx = 0;
			int last_k_pos = frags[fragidx].length() - k_;

			// --- identify unique kmers in the frag ---
			u64set_t frag_kval_set;
			while (k_pos_idx <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k_));

				if (xny::str2ID<uint64_t>(val, kmer)) {
					frag_kval_set.insert(val);
					++ k_pos_idx;
				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos_idx += pos_n + 1;
					} else k_pos_idx += k_;
				}
			} // while (k_pos_idx <= last_k_pos)

			for (auto& val: frag_kval_set) {
				if (std::binary_search(kmerTaxaID_.begin(), kmerTaxaID_.end(),
						u64ipair_t(val, 0), cmp_u64ipair ())) {
					local_kval_vec.push_back(val);
				}
			}
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

		#pragma omp critical
		{
			kval_vec.insert(kval_vec.end(), local_kval_vec.begin(),
					local_kval_vec.end());
		}
	}

	accrue_kcnts_ (kval_vec);

} // ProcDB::verify_kmer_in_frags
*/




/*
void ProcDB::link_sampled_kmers_via_reads (const strvec_t& ipfqs, const strvec_t& isfqs) {

	if (!silent_) std::cout << "Link sampled kmers via reads ...\n";

	// process paired-end files
	for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {
		if (!silent_) std::cout << "\tproc paired-end files: \n\t\t"
				<< ipfqs[fidx-1] << "\n\t\t" << ipfqs[fidx] << "\n";

		std::ifstream fh0, fh1;
		xny::openfile<std::ifstream>(fh0, ipfqs[fidx-1]);
		xny::openfile<std::ifstream>(fh1, ipfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end && fq1 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch_/2, fq0, end);
	 		add_fq_reads_only (frags, batch_/2, fq1, end);
	 		//create_frags (frags);

	 		link_sampled_kmer_in_frags (frags);

			total_frags += frags.size();

			frags.clear();

		} // while

		if (!silent_) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";

		xny::closefile(fh0);
		xny::closefile(fh1);
	} // for (int fidx = 1; fidx < ipfqs.size(); fidx += 2) {

	// process single-end files
	for (int fidx = 0; fidx < isfqs.size(); ++ fidx) {
		if (!silent_) std::cout << "\tproc single-end file: " << isfqs[fidx] << "\n";

		std::ifstream fh0;
		xny::openfile<std::ifstream>(fh0, isfqs[fidx]);
		bio::fastq_input_iterator<> fq0 (fh0), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {

	 		add_fq_reads_only (frags, batch_, fq0, end);

	 		link_sampled_kmer_in_frags (frags);

			total_frags += frags.size();

 			frags.clear();
		} // while

		if (!silent_) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";
		xny::closefile(fh0);
	}

	if (!silent_) print_kstat ();
}
*/
/*
void ProcDB::count_kmers_in_read_files_(const u64vec_t& kmers,
		const strvec_t& fl, bool is_pe) {

	for (int fidx = 0; fidx < fl.size(); ++ fidx) {

		if (!silent_) {
			if (is_pe) {
				++ fidx;
				std::cout << "\tproc paired-end files: \n\t\t"
					<< fl[fidx - 1] << "\n\t\t" << fl[fidx] << "\n";
			} else std::cout << "\tproc single-end file: " << fl[fidx] << "\n";
		}

		std::ifstream fh0, fh1;
		if (is_pe) {
			xny::openfile<std::ifstream>(fh0, fl[fidx-1]);
			xny::openfile<std::ifstream>(fh1, fl[fidx]);
		} else xny::openfile<std::ifstream>(fh0, fl[fidx]);

		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch_/2, fq0, end);
	 		if (is_pe) {
	 			while (fq1 != end) {
	 				add_fq_reads_only (frags, batch_/2, fq1, end);
	 		 		create_frags (frags);
	 			}
	 		}

	 		count_kmers_in_frags (kmers, frags);

	 		total_frags += frags.size();
			frags.clear();

		} // while

		if (!silent_) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";

		xny::closefile(fh0);
		if (is_pe) xny::closefile(fh1);
	} // for (int fidx = 0; fidx < fl.size(); ++ fidx) {
} // ProcDB::count_kmers_in_read_files_
*/

/* @brief	Identify [kmers] that exist in frags
 */
/*
void ProcDB::count_kmers_in_frags (const u64vec_t& kmers, const strvec_t& frags) {

	int num_frag = frags.size();

	// whenever a traversed kmer x in [frags] exists in [kmers]
	// x will be added to kval_vec
	u64vec_t kval_vec;

	#pragma omp parallel
	{
		u64vec_t local_kval_vec; // local to each thread
		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {
			int k_pos_idx = 0;
			int last_k_pos = frags[fragidx].length() - k_;

			// --- identify unique kmers in the frag ---
			u64set_t frag_kval_set;
			while (k_pos_idx <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k_));

				if (xny::str2ID<uint64_t>(val, kmer)) {
					frag_kval_set.insert(val);
					++ k_pos_idx;
				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos_idx += pos_n + 1;
					} else k_pos_idx += k_;
				}
			} // while (k_pos_idx <= last_k_pos)

			for (auto& val: frag_kval_set) {
				if (std::binary_search(kmers.begin(), kmers.end(), val)){
					local_kval_vec.push_back(val);
				}
			}
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

		#pragma omp critical
		{
			kval_vec.insert(kval_vec.end(), local_kval_vec.begin(),
					local_kval_vec.end());
		}
	}

	accrue_kcnts_ (kval_vec);

} // ProcDB::count_sampled_kmer_in_frags (const u64vec_t& kmers, const strvec_t& frags)
*/

/* @brief Given a vector of kmer values, add these counts to [kmerInfo_]
 */
/*
void ProcDB::accrue_kcnts_ (u64vec_t& kvals) {
	if (!kvals.size()) return;
	std::sort (kvals.begin(), kvals.end());

	// convert kvals to <kmer, cnt> pairs
	std::vector<u64ipair_t> kcnt;
	kcnt.push_back(u64ipair_t (kvals.front (), 0)); // first element
	for (auto& x: kvals) { // ignore first element
		if (kcnt.back().first == x) ++ kcnt.rbegin()->second;
		else kcnt.push_back(u64ipair_t (x, 1));
	}
	kvals.clear();

	std::vector<u64ibtuple_t> merged_kmerinfo;

	int p0 = 0, p1 = 0;
	while (p0 < kmerInfo_.size() && p1 < kcnt.size()) {
		uint64_t val0 = std::get<0>(kmerInfo_[p0]), val1 = kcnt[p1].first;
		if (val0 == val1) {
			merged_kmerinfo.push_back (std::make_tuple (kcnt[p1].first,
					kcnt[p1].second + std::get<1> (kmerInfo_[p0]), false));
			++ p0;
			++ p1;
		} else if (val0 < val1) {
			merged_kmerinfo.push_back (kmerInfo_[p0]);
			++ p0;
		} else {
			merged_kmerinfo.push_back (std::make_tuple (kcnt[p1].first,
					kcnt[p1].second, false));
			++ p1;
		}
	}

	// adding the remaining array
	if (p0 < kmerInfo_.size()) {
		merged_kmerinfo.insert(merged_kmerinfo.end(), kmerInfo_.begin() + p0,
				kmerInfo_.end());
	} else {
		for (; p1 < kcnt.size(); ++ p1) {
			merged_kmerinfo.push_back (std::make_tuple (kcnt[p1].first,
					kcnt[p1].second, false));
		}
	}
	set_kmerinfo_ (merged_kmerinfo);

} // ProcDB::accrue_kcnts_
*/

/*
void ProcDB::count_kpair_in_read_files_(std::vector<uu64ituple_t>& rkpcnt,
		int min_freq, const strvec_t& fl, bool is_pe) {

	for (int fidx = 0; fidx < fl.size(); ++ fidx) {

		if (!silent_) {
			if (is_pe) {
				++ fidx;
				std::cout << "\tproc paired-end files: \n\t\t"
					<< fl[fidx - 1] << "\n\t\t" << fl[fidx] << "\n";
			} else std::cout << "\tproc single-end file: " << fl[fidx] << "\n";
		}

		std::ifstream fh0, fh1;
		if (is_pe) {
			xny::openfile<std::ifstream>(fh0, fl[fidx-1]);
			xny::openfile<std::ifstream>(fh1, fl[fidx]);
		} else xny::openfile<std::ifstream>(fh0, fl[fidx]);

		bio::fastq_input_iterator<> fq0 (fh0), fq1 (fh1), end;

		int total_frags = 0;
		strvec_t frags;

		while (fq0 != end) {
	 		// cat paired-end reads to be a single frag separated by an 'n'
	 		add_fq_reads_only (frags, batch_/2, fq0, end);
	 		if (is_pe) {
	 			while (fq1 != end) {
	 				add_fq_reads_only (frags, batch_/2, fq1, end);
	 		 		create_frags (frags);
	 			}
	 		}

	 		std::vector<uu64pair_t> rkp;
	 		count_kpair_in_frags_(rkp, min_freq, frags);
	 		accrue_kpcnts_(rkpcnt, rkp);

			total_frags += frags.size();
			frags.clear();

		} // while

		if (!silent_) std::cout << "\t\ttotal frags: " << total_frags <<  "\n";

		xny::closefile(fh0);
		if (is_pe) xny::closefile(fh1);
	} // for (int fidx = 0; fidx < fl.size(); ++ fidx) {
} // ProcDB::count_kpair_in_read_files_
*/

/*	@brief	Record <k1,k2,cnt> in [rkpcnt] if k1 and k2 belong to the
 * same read and from different taxa
 */
/*
void ProcDB::count_kpair_in_frags_(std::vector<uu64pair_t>& rkp,
		int min_freq, const strvec_t& frags) {

	int num_frag = frags.size();

	#pragma omp parallel
	{
		std::vector<uu64pair_t> local_rkp; // local to each thread
		#pragma omp for
		for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {
			int k_pos_idx = 0;
			int last_k_pos = frags[fragidx].length() - k_;

			// --- identify unique kmers in the frag ---
			u64vec_t frag_kval_vec;
			while (k_pos_idx <= last_k_pos) {
				//look for current kmer
				uint64_t val;
				std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k_));

				if (xny::str2ID<uint64_t>(val, kmer)) { // contain no "n"
					if (get_kmerinfo_freq_(val) > min_freq) { // cnt > min_freq
						frag_kval_vec.push_back(val);
					}
					++ k_pos_idx;
				} else { // skip next found 'n/N'
					// find the last position of N
					int pos_n = kmer.find_last_of("nN");
					if (pos_n != std::string::npos) {
						k_pos_idx += pos_n + 1;
					} else k_pos_idx += k_;
				}
			} // while (k_pos_idx <= last_k_pos)

			rm_duplicated_elements (frag_kval_vec);

			// for each kmer, obtain list of taxaIDs they belong to
			int sz = frag_kval_vec.size();
			std::vector<iset_t> taxaIDs (sz);
			for (int i = 0; i < sz; ++ i) {
				taxaIDs[i] = get_kmertaxaid_tids_(frag_kval_vec[i]);
			}

			// pairwise linking kmers
			for (int i = 0; i < sz - 1; ++ i) {
				for (int j = i + 1; j < sz; ++ j) {
					iset_t combined = taxaIDs[i];
					combined.insert (taxaIDs[j].begin(), taxaIDs[j].end());
					if (combined.size() > 1) {
						local_rkp.push_back(	uu64pair_t (
							std::min(frag_kval_vec[i], frag_kval_vec[j]),
							std::max(frag_kval_vec[i], frag_kval_vec[j])));
					}
				}
			}
		} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

		#pragma omp critical
		{
			rkp.insert(rkp.end(), local_rkp.begin(), local_rkp.end());
		}
	}
} //ProcDB::count_kpair_in_frags_
*/
/* @brief
 *
 */
/*
void ProcDB::accrue_kpcnts_(std::vector<uu64ituple_t>& rkpcnt,
		std::vector<uu64pair_t>& rkp) {
	if (!rkp.size()) return;
	std::sort (rkp.begin(), rkp.end());

	// convert rkp to <k1, k2, cnt> structure
	std::vector<uu64ituple_t> rhs_rkpcnt;

	rhs_rkpcnt.push_back(std::make_tuple (rkp[0].first, rkp[0].second, 0)); // first element
	for (auto& x: rkp) { // ignore first element
		if (uu64pair_t(std::get<0>(rhs_rkpcnt.back()),
			std::get<1>(rhs_rkpcnt.back())) == x) ++ std::get<2>(rhs_rkpcnt.back());
		else rhs_rkpcnt.push_back(std::make_tuple (x.first, x.second, 1));
	}
	rkp.clear();

	std::vector<uu64ituple_t> merged_rkpcnt;

	int p0 = 0, p1 = 0;
	while (p0 < rkpcnt.size() && p1 < rhs_rkpcnt.size()) {
		uu64pair_t kp0 {std::get<0>(rkpcnt[p0]), std::get<1>(rkpcnt[p0])},
				   kp1 {std::get<0>(rhs_rkpcnt[p1]), std::get<1>(rhs_rkpcnt[p1])};

		if (kp0 == kp1) {
			merged_rkpcnt.push_back (std::make_tuple (kp0.first, kp0.second,
				std::get<2>(rkpcnt[p0]) + std::get<2>(rhs_rkpcnt[p1])));
			++ p0;
			++ p1;
		} else if (kp0 < kp1) {
			merged_rkpcnt.push_back (rkpcnt[p0]);
			++ p0;
		} else {
			merged_rkpcnt.push_back (rhs_rkpcnt[p1]);
			++ p1;
		}
	}
	// adding the remaining array
	if (p0 < rkpcnt.size()) {
		merged_rkpcnt.insert(merged_rkpcnt.end(), rkpcnt.begin() + p0,
				rkpcnt.end());
	} else {
		merged_rkpcnt.insert(merged_rkpcnt.end(), rhs_rkpcnt.begin() + p1,
				rhs_rkpcnt.end());
	}

	rkpcnt = merged_rkpcnt;

} // ProcDB::accrue_kpcnts_
*/

/* @brief	Step-wise sampling each taxa and write skeletons to files
 */
/*
void ProcDB::create_skeletons_per_taxa (ivec_t& taxaKCnt) {

	if (!silent_) std::cout << "Create skeletons by step-wise sampling each taxa ...\n";

	int num_taxa = get_num_taxa_();
	taxaKCnt.resize(num_taxa);
	if (!num_taxa) return;

	if (!silent_) std::cout << "\t" << num_taxa << " taxa to process\n";

	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < num_taxa; ++ i) {
		int fnum = fileList_[i].size();
		uu64vec_t skeletons;
		ivec_t GIs;
		for (int j = 0; j < fnum; ++ j) {
			sample_kmer_from_file_(GIs, skeletons, fileList_[i][j]);
		}
		output_skeletons_(GIs, skeletons, skeletonFileList_[i]);
		int num_kmers = 0;
		for (auto& x: skeletons) taxaKCnt[i] += x.size();
	}

	if (!silent_) std::cout << "\tcomplete.\n";

} // ProcDB::create_taxa_skeletons
*/
/*
void ProcDB::sample_kmer_from_file_(ivec_t& GIs, uu64vec_t& skeletons,
		const std::string& file) {
	std::ifstream fh (file.c_str());
	if (!fh.good()) abording ("ProcDB::sample_kmer_from_file: can't open " + file);

	std::string contig;
	std::string line;
	int gi = -1;
	u64vec_t skeleton;

	while (std::getline (fh, line)) {
		if (line.empty()) continue;
		if (line.at(0) == '>') {
			skeleton.clear();
			if (!contig.empty()) {
				add_kmer_from_contig (skeleton, contig);
				skeletons.push_back(skeleton);
				GIs.push_back(gi);
				contig.clear();
			}
			get_gi_from_header_(gi, line);
		} else contig += line;
	}
	skeleton.clear();
	if (!contig.empty()) {
		add_kmer_from_contig(skeleton, contig);
		skeletons.push_back(skeleton);
		GIs.push_back(gi);
	}
	fh.close();
} // ProcDB::sample_kmer_from_file
*/
