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


/* @brief	Given taxa file list folder, set taxa [fileList_], which
 * is a 2-d string vector, the 1st-D specifies the taxon and the 2nd-D
 * specifies the list of files of this taxon
 */
void ProcDB::set_db_filelist(const std::string& DBDir, const strvec_t& DBFileList) {

	if (!silent_) std::cout << "Read in DB files...\n";

	if (!DBDir.empty()) {
		DIR* dir = opendir (DBDir.c_str());
		if (!dir) abording ("cannot open " + DBDir);
		struct dirent *dirp;

		std::vector<ipair_t> vec_gi_len;
		std::map<std::string, int> filePrefix_fileIdx;

		while ((dirp = readdir(dir))) {
			std::string filename = dirp->d_name;
			unsigned pos = filename.find_last_of('.');
			std::string suffix = filename.substr (pos + 1);
			std::string prefix = xny::get_prefix (filename, true, ".");
			std::string path = DBDir + "/" + filename;
			if (suffix.compare("fa") == 0 || suffix.compare("fasta") == 0) {
				auto it = filePrefix_fileIdx.find(prefix);
				if (it != filePrefix_fileIdx.end()) {
					fileList_[it->second].push_back(path);
				} else {
					fileList_.push_back(strvec_t(1, path));
					filePrefix_fileIdx[prefix] = fileList_.size() - 1;
					add_taxaID (std::atoi(prefix.c_str()));
				}
			}
		}
		closedir (dir);
	} // 	if (!DBDir.empty()) {

	for (auto& f: DBFileList) {
		std::string suffix = xny::get_suffix (f, true, "/");
		add_taxaID (std::atoi(xny::get_prefix (suffix, true, ".").c_str()));

		fileList_.push_back (strvec_t (1, f));
	}

 	if (!silent_) std::cout << "\tnum of taxa: " <<  get_num_taxa_() << "\n";

} //ProcDB::set_db_filelist

/* @brief	Set up file list to store skeleton for each taxa
 */
void ProcDB::set_db_skeleton_filelist (const std::string& odir) {
	 if (!taxaIDs_.size()) abording ("ProcDB::set_db_skeleton_filelist SC"
			 " taxaIDs_ empty: cannot set filelist");
	 if (! boost::filesystem::exists(odir)) {
		 if (! boost::filesystem::create_directory(odir)) {
		     abording ("ProcDB::set_db_skeleton_filelist: could not create "
		    		 + odir);
		 }
	 }

	 for (auto& x: taxaIDs_) {
		 std::string path = odir + "/";
		 path += std::to_string (x) + ".txt";
		 skeletonFileList_.push_back(path);
	 }
} // ProcDB::set_db_skeleton_filelist

/* @brief	Step-wise sampling each taxa to generate skeletons, and write
 * these to files. Taking care of over-sampling issue where many highly
 * similar sequences have been included as references for some taxa
 */
void ProcDB::generate_skeleton_per_taxa(ivec_t& taxaKCnt) {

	if (!silent_) std::cout << "Generate skeletons by sampling each taxa ...\n";

	int num_taxa = get_num_taxa_();
	taxaKCnt.resize(num_taxa);
	if (!num_taxa) return;

	if (!silent_) std::cout << "\t" << num_taxa << " taxa to process\n";

	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < num_taxa; ++ i) {
		int fnum = fileList_[i].size();
		uu64vec_t skeletons;
		ivec_t GIs;
		u64set_t sampled_kmers;
		for (int j = 0; j < fnum; ++ j) {
			step_wise_sampling_kmer_from_file_(GIs, skeletons,
					sampled_kmers, fileList_[i][j]);
		}
		output_skeletons_(GIs, skeletons, skeletonFileList_[i]);
		int num_kmers = 0;
		for (auto& x: skeletons) taxaKCnt[i] += x.size();

		/*{ // debug print
			std::cout << i << " : " << sampled_kmers.size() << "\n";
		}*/
	}

	if (!silent_) std::cout << "\tcomplete.\n";

} // generate_skeleton_per_taxa

/* @brief Step-wise sample each contig by first create a skeleton of this
 * using kmers that have been sampled previously, then create skeleton for
 * regions that have not been previously sampled and add the corresponding
 * kmers to the sampled kmer list to be used for remaining contigs.
 *
 * $param GIs: gi number corresponding to each contig
 * $param skeletons: skeleton for each contig (1-1 match with GIs)
 * $param sampled_kmers: sorted kmer list that have already been sampled
 * $param file: contig file
 */
void ProcDB::step_wise_sampling_kmer_from_file_(ivec_t& GIs,
	uu64vec_t& skeletons, u64set_t& sampled_kmers, const std::string& file){

	std::ifstream fh (file.c_str());
	if (!fh.good()) abording (
		"ProcDB::step_wise_sampling_kmer_from_file_: can't open " + file);

	//int num_contig = 0;
	std::string contig;
	std::string line;
	int gi = -1;
	u64vec_t skeleton;
	bool is_debug = false;
	while (std::getline (fh, line)) {
		if (line.empty()) continue;
		if (line.at(0) == '>') {
			//std::cout << "contig, size, sampled_kmers: " << num_contig << ", " << contig.size() << ", " << sampled_kmers.size() << "\n";
			//++ num_contig;

			skeleton.clear();
			if (!contig.empty()) {
				sample_kmer_from_contig_(skeleton, sampled_kmers, contig);
				if (is_debug) std::cout << skeleton.size() << "\n";
				// adding skeleton into sampled_kmers
				/*sampled_kmers.insert(sampled_kmers.end(),
						skeleton.begin(), skeleton.end());
				rm_duplicated_elements(sampled_kmers); */
				if (skeleton.size()) {
					sampled_kmers.insert(skeleton.begin(), skeleton.end());
					skeletons.push_back(skeleton);
					GIs.push_back(gi);
				}
				contig.clear();
			}
			get_gi_from_header_(gi, line);
			{
				if (gi == 445351323) {
					std::cout << gi << "\n";
					is_debug = true;
				}
			}

		} else contig += line;
	}
	skeleton.clear();
	if (!contig.empty()) {
		sample_kmer_from_contig_(skeleton, sampled_kmers, contig);

		// adding skeleton into sampled_kmers
		/*sampled_kmers.insert(sampled_kmers.end(),
				skeleton.begin(), skeleton.end());
		rm_duplicated_elements(sampled_kmers);*/
		if (skeleton.size()) {
			sampled_kmers.insert(skeleton.begin(), skeleton.end());
			skeletons.push_back(skeleton);
			GIs.push_back(gi);
		}
	}
	fh.close();
} // ProcDB::step_wise_sampling_kmer_from_file_

/* @brief	Sample kmers from contig by trying to accomodate kmers that
 * have been previously sampled, newly identified kmers are added to
 * sampled_kmers
 *
 * $kmers: newly sampled kmers
 * $sampled_kmers: kmers that have been included in the taxa skeletons
 * $contig: contig under consideration
 */
void ProcDB::sample_kmer_from_contig_(u64vec_t& kmers,
			const u64set_t& sampled_kmers, const std::string& contig){

	xny::low_complexity lc (30, 50, 80);
	int len = contig.length();
	int num = (len - k_ + 1) * perc_/100 + 1;
	if (num <= 0) return;
	int step = 100 / perc_; // equal to (len - k_ + 1) / num;
	if (!sampled_kmers.size()) {
		for (int i = 0; i < num; ++ i) { //first char a or c, use as is, otherwise rvc
			uint64_t val;
			std::string kmer = get_kstr (contig.substr(i*step, k_));
			if (!lc (kmer) && xny::str2ID<uint64_t>(val, kmer)) kmers.push_back(val);
		}
	} else { //
		u64vec_t klist; // all kmers from contig
		xny::get_bitkmer_by_alphabet<std::back_insert_iterator<u64vec_t>, uint64_t>
		(contig, std::back_inserter(klist), k_);

 		// identify skeleton of contig by checking which kmers have been sampled
		ivec_t sampled_index;
		for (int i = 0; i < klist.size(); ++ i) {
			/*if (std::binary_search(sampled_kmers.begin(),
					sampled_kmers.end(), klist[i])){ // found
				sampled_index.push_back(i);
			}*/
			if (sampled_kmers.count(klist[i])) sampled_index.push_back(i);
		}

		//std::cout << "index:\n";
		//for (auto& x:sampled_index) std::cout << x << "\t";
		//std::cout << "\n";

		// add to [sampled_index] positions that should have been sampled
		if (!sampled_index.size()) {
			for (int i = 0; i < num; ++ i) {
				int idx = i * step;
				if (idx < klist.size() && !lc (xny::ID2Str(klist[idx], k_)))
					sampled_index.push_back(idx);
			}
		} else {
			ivec_t idx_addition;
			// deal with the first index
			int num_to_sample = sampled_index.front() / step;
			for (int i = 0; i < num_to_sample; ++ i) {
				int idx = i * step;
				if (!lc (xny::ID2Str(klist[idx], k_)))
					idx_addition.push_back(idx);
			}
			// deal with intermediates
			for (int i = 1; i < sampled_index.size() - 1; ++ i) {
				num_to_sample = (sampled_index[i + 1] - sampled_index[i]) / step - 1;
				for (int j = 1; j <= num_to_sample; ++ j) {
					int idx = sampled_index[i] + j * step;
					if (idx < klist.size() && !lc (xny::ID2Str(klist[idx], k_)))
						idx_addition.push_back(idx);
				}
			}
			// deal with the last index
			num_to_sample = ((len - k_ + 1) - sampled_index.back()) / step;
			for (int i = 1; i <= num_to_sample; ++ i) {
				int idx = sampled_index.back() + i * step;
				if (idx < klist.size() && !lc (xny::ID2Str(klist[idx], k_)))
					idx_addition.push_back(idx);
			}

			//--- insert additional index and sort ---
			int sampled_sz = sampled_index.size(),
				total_sz = sampled_sz + idx_addition.size();

			// dont' add anything if 95% of genome has already been sampled
			if (100 * sampled_sz / total_sz > 95) {
				sampled_index.clear();
			} else {
				sampled_index.insert(sampled_index.end(), idx_addition.begin(),
						idx_addition.end());
				std::sort(sampled_index.begin(), sampled_index.end());
			}
		} // else

		// now get kmers via sampled_index
		for (auto& idx: sampled_index) kmers.push_back(klist[idx]);

	} // else

} // ProcDB::sample_kmer_from_contig_

void ProcDB::output_skeletons_(const ivec_t& GIs,
		const uu64vec_t& skeletons, const std::string& ofile) {
	if (GIs.size() != skeletons.size()) {
		abording ("ProcDB::output_skeletons_ failed");
	}
	std::ofstream ofh;
	xny::openfile(ofh, ofile);
	ofh << skeletons.size() << "\n";

	for (int i = 0; i < GIs.size(); ++ i) {
		ofh << GIs[i] << "\t" << skeletons[i].size(); // gi num_elem e1 e2 ...
		for (auto& elem: skeletons[i]) ofh << "\t" << elem;
		ofh << "\n";
	}
	xny::closefile(ofh);
} // ProcDB::output_skeletons_


