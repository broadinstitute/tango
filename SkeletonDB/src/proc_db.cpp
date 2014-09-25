//========================================================================
// Project     : SkeletonDB
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
 * specifies the list of files of this taxon.
 * Meanwhile set output taxa [ofileList_], where each entry specifies
 * the output file path of the taxa.
 */
void ProcDB::set_db_filelist(const std::string& DBDir,
		const strvec_t& DBFileList, const std::string& odir) {

	if (!silent_) std::cout << "Read in DB files...\n";

	std::map<std::string, int> filePrefix_fileIdx;

	if (!DBDir.empty()) {
		DIR* dir = opendir (DBDir.c_str());
		if (!dir) abording ("cannot open " + DBDir);
		struct dirent *dirp;

		while ((dirp = readdir(dir))) {
			std::string filename = dirp->d_name;
			unsigned pos = filename.find_last_of('.');
			std::string suffix = filename.substr (pos + 1);
			std::string prefix = xny::get_prefix (filename, true, ".");
			std::string path = DBDir;
			path += "/";
			path += filename;
			if (suffix.compare("fa") == 0 || suffix.compare("fasta") == 0) {
				auto it = filePrefix_fileIdx.find(prefix);
				if (it != filePrefix_fileIdx.end()) {
					fileList_[it->second].push_back(path);
				} else {
					fileList_.push_back(strvec_t(1, path));
					filePrefix_fileIdx[prefix] = fileList_.size() - 1;
					add_taxaID (std::atoi(prefix.c_str()));

					if (!odir.empty()) {
						std::string opath = odir;
						opath += "/";
						opath += prefix;
						opath += ".fa";
						ofileList_.push_back(opath);
					}
				}
			}
		}
		closedir (dir);
	} // 	if (!DBDir.empty()) {

	for (auto& f: DBFileList) {

		std::string suffix = xny::get_suffix (f, true, "/");
		std::string prefix = xny::get_prefix (suffix, true, ".");

		auto it = filePrefix_fileIdx.find(prefix);
		if (it != filePrefix_fileIdx.end()) {
			fileList_[it->second].push_back(f);
		} else {
			fileList_.push_back (strvec_t (1, f));
			filePrefix_fileIdx.insert(std::pair<std::string, int>(prefix,
					fileList_.size() - 1));
			add_taxaID (std::atoi(prefix.c_str()));
			if (!odir.empty()) {
				std::string opath = odir;
				opath += "/";
				opath += prefix;
				opath += ".fa";
				ofileList_.push_back(opath);
			}
		}
	}

 	if (!silent_) std::cout << "\tnum of taxa: " <<  get_num_taxa_() << "\n";

} //ProcDB::set_db_filelist

/* @brief	Set up file list to store skeleton for each taxa
 */
void ProcDB::set_db_skeleton_filelist (const std::string& odir) {
	 if (!taxaIDs_.size()) abording ("ProcDB::set_db_skeleton_filelist SC"
			 " taxaIDs_ empty: cannot set filelist");

	 for (auto& x: taxaIDs_) {
		 std::string path = odir;
		 path += "/";
 		 path += std::to_string (x);
 		 path += ".txt";
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
		std::ofstream ofh;
		xny::openfile(ofh, ofileList_[i]);

		uu64vec_t skeletons;
		ivec_t GIs;
		u64set_t sampled_kmers;
		for (int j = 0; j < fnum; ++ j) {
			step_wise_sampling_kmer_from_file_(GIs, skeletons,
					sampled_kmers, ofh, fileList_[i][j]);
		}
		output_skeletons_(GIs, skeletons, skeletonFileList_[i]);
		int num_kmers = 0;
		for (auto& x: skeletons) taxaKCnt[i] += x.size();
		xny::closefile(ofh);
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
	uu64vec_t& skeletons, u64set_t& sampled_kmers, std::ofstream& ofh,
	const std::string& file){

	std::ifstream fh (file.c_str());
	if (!fh.good()) abording (
		"ProcDB::step_wise_sampling_kmer_from_file_: can't open " + file);

	int num_contig = 0;
	std::string contig;
	std::string line;
	int gi = -1;
	uu64vec_t low_density_k_windows;
	iivec_t low_density_index_windows;
	bool is_debug = false;
	while (std::getline (fh, line)) {
		if (line.empty()) continue;
		if (line.at(0) == '>') {
			if (contig.length() >= minl_) {

				sample_kmer_from_contig_(low_density_k_windows,
						low_density_index_windows, sampled_kmers, contig);

				for(auto& w: low_density_k_windows) {
					sampled_kmers.insert(w.begin(), w.end());
					skeletons.push_back(w);
					GIs.push_back(gi);
				}

				output_contig(ofh, gi, low_density_index_windows, contig);

			}
			contig.clear();

			get_gi_from_header_(gi, line);

		} else contig += line;
	}

	// the last contig
 	if (contig.length() > minl_) {
		sample_kmer_from_contig_(low_density_k_windows,
				low_density_index_windows,  sampled_kmers, contig);

		// adding skeleton into sampled_kmers
		//sampled_kmers.insert(sampled_kmers.end(),
		//		skeleton.begin(), skeleton.end());
		//rm_duplicated_elements(sampled_kmers);
		for(auto& w: low_density_k_windows) {
			sampled_kmers.insert(w.begin(), w.end());
			skeletons.push_back(w);
			GIs.push_back(gi);
			//{ // debug print
			//	for (auto& x: w) std::cout << x << "\t";
			//	std::cout << "\n";
			//}
		}
		output_contig(ofh, gi, low_density_index_windows, contig);


		//if (skeleton.size()) {
		//	sampled_kmers.insert(skeleton.begin(), skeleton.end());
		//	skeletons.push_back(skeleton);
		//	GIs.push_back(gi);
		//}
	}

	fh.close();
} // ProcDB::step_wise_sampling_kmer_from_file_

/* @brief	Sample kmers from contig by trying to accommodate kmers that
 * have been previously sampled, newly identified kmers are added to
 * sampled_kmers, report only low density kmer_windows
 *
 * $low_density_k_windows: low density kmer windows
 * $sampled_kmers: kmers that have been included in the taxa skeletons
 * $contig: contig under consideration
 */
void ProcDB::sample_kmer_from_contig_(uu64vec_t& low_density_k_windows,
		iivec_t& low_density_index_windows, const u64set_t& sampled_kmers,
		const std::string& contig){
	low_density_k_windows.clear();
	low_density_index_windows.clear();

	xny::low_complexity lc (30, 50, 80);
	int len = contig.length();
	int num = (len - k_ + 1) * perc_/100 + 1;
	if (num <= 0) return;
	int step = 100 / perc_; // equal to (len - k_ + 1) / num;
	u64vec_t ld_k_window;
	ivec_t ld_index_window;
	if (!sampled_kmers.size()) { // the first contig
		for (int i = 0; i < num; ++ i) { //first char a or c, use as is, otherwise rvc
			uint64_t val;
			std::string kmer = get_kstr (contig.substr(i*step, k_));
			if (!lc (kmer) && xny::str2ID<uint64_t>(val, kmer)) {
				ld_k_window.push_back(val);
				ld_index_window.push_back(i*step);
			}
		}
		if (ld_k_window.size()) {
			low_density_k_windows.push_back(ld_k_window);
			low_density_index_windows.push_back(ld_index_window);
			ld_k_window.clear();
		}
	} else { //
		/*{
			std::cout << "\tget_bit_kmer\n";
		}*/
		u64vec_t klist; // all kmers from contig
		bdeque_t is_valid; // corresponding kmer is valid or not
		xny::get_bitkmer_by_alphabet<std::back_insert_iterator<u64vec_t>, uint64_t>
		(contig, std::back_inserter(klist), is_valid, k_);
		/*{
			std::cout << "\tidentify_skeleton_contig, sampled_k_size: "
					<< sampled_kmers.size() << "\n";
		}*/
 		// identify skeleton of contig by checking which kmers have been sampled
		ivec_t sampled_index;

		// replace this section with nested omp
		/*for (int i = 0; i < klist.size(); ++ i) {
			if (sampled_kmers.count(klist[i])) sampled_index.push_back(i);
		}*/


		#pragma omp parallel // nested parallel region
		{
			ivec_t local_sampled_index;
			#pragma omp for
			for (int i = 0; i < klist.size(); ++ i) {
				if (sampled_kmers.count(klist[i])) local_sampled_index.push_back(i);
			}
			#pragma omp critical
			{
				sampled_index.insert(sampled_index.end(), local_sampled_index.begin(),
						local_sampled_index.end());
			}
		}

		std::sort (sampled_index.begin(), sampled_index.end());

		/*{
			std::cout << "\tid_low_compl_w\n";
		}*/

		// add to [sampled_index] positions that should have been sampled
		if (!sampled_index.size()) {
			for (int i = 0; i < num; ++ i) {
				int idx = i * step;
				if (idx < klist.size() && is_valid[idx] &&
						!lc (xny::ID2Str(klist[idx], k_)))
					sampled_index.push_back(idx);
			}
			if (sampled_index.size()) {
				low_density_index_windows.push_back(sampled_index);
			}

		} else {

			// additional indices are stored as negative values
			ivec_t idx_addition;
			identify_additional_indices_(idx_addition, sampled_index,
					klist, is_valid, step, len, lc);
			// merge and sort all indices
			sampled_index.insert(sampled_index.end(), idx_addition.begin(),
					idx_addition.end());
			std::sort(sampled_index.begin(), sampled_index.end(), cmp_abs());
			idx_addition.clear();

			//--- now identify low density windows via sampled_index ---
			/*{ // debug print
				std::cout << "sampled index:\n";
				for (auto& x: sampled_index) std::cout << x << "\t";
				std::cout << "\n";
			}*/

			identify_low_density_windows(low_density_index_windows, sampled_index);

		} // else

		// identify all kmers specified by index_windows and generate
		// low_density_k_windows
		for (auto& w: low_density_index_windows) {
			for (auto& idx: w) {
				ld_k_window.push_back(klist[std::abs(idx)]);
			}
			low_density_k_windows.push_back(ld_k_window);
			ld_k_window.clear();
		}

		// now get kmers via sampled_index
		//for (auto& idx: sampled_index) kmers.push_back(klist[std::abs(idx)]);

	} // else

} // ProcDB::sample_kmer_from_contig_


void ProcDB::identify_additional_indices_(ivec_t& idx_addition,
	const ivec_t& sampled_index, const u64vec_t& klist, bdeque_t& is_valid,
	int step, int contig_len, xny::low_complexity& lc){

	// deal with the first index
	int num_to_sample = sampled_index.front() / step;
	for (int i = 0; i < num_to_sample; ++ i) {
		int idx = i * step;
		if (is_valid[idx] && !lc (xny::ID2Str(klist[idx], k_)))
			idx_addition.push_back(-1 * idx);
	}
	// deal with intermediates
	for (int i = 0; i < sampled_index.size() - 1; ++ i) {
		num_to_sample = (sampled_index[i + 1] - sampled_index[i]) / step - 1;
		for (int j = 1; j <= num_to_sample; ++ j) {
			int idx = sampled_index[i] + j * step;
			if (idx < klist.size() && is_valid[idx]
			     && !lc (xny::ID2Str(klist[idx], k_)))
				idx_addition.push_back(-1 * idx);
		}
	}
	// deal with the last index
	num_to_sample = ((contig_len - k_ + 1) - sampled_index.back()) / step;
	for (int i = 1; i <= num_to_sample; ++ i) {
		int idx = sampled_index.back() + i * step;
		if (idx < klist.size() && is_valid[idx]
		    && !lc (xny::ID2Str(klist[idx], k_)))
			idx_addition.push_back(-1 * idx);
	}
}

/* Segment [coordinates] into smaller windows where each has low density
 * 1.Linear scan index and mark all positions that a low density window
 * starts
 * 2. merge all low density windows
 * 3. Identify the first and the last unsampled index (negative value)
 * for each merged window and
 */
void ProcDB::identify_low_density_windows(iivec_t& windows,
		const ivec_t& coordinates){
	if (!coordinates.size()) return;
	int num_w = (int) coordinates.size() - wl_ + 1;
	int num_holes = 0;
	if (num_w <= 0) { // window is < target minimum length
		for (auto& c: coordinates) if (c < 0) ++ num_holes;
		if (100 - 100 * num_holes / wl_ < wd_) windows.push_back(coordinates);
		return;
	}

	ivec_t marked_indices;
	// the first window
	for (int i = 0; i < wl_; ++ i) if (coordinates[i] < 0) ++ num_holes;
	if (100 - 100 * num_holes / wl_ < wd_) marked_indices.push_back(0);

	for (int i = 1; i < num_w; ++ i) {
		if (coordinates[i + wl_ - 1] < 0) {
			if(coordinates[i - 1] >= 0) ++ num_holes;
		} else if (coordinates[i - 1] < 0) -- num_holes;
		//std::cout << i << ": " << num_holes << "\n";
		if (100 - 100 * num_holes / wl_ < wd_) marked_indices.push_back(i);
	}
	/*{ // debug print
		std::cout << "marked indices: ";
		for (auto& x: marked_indices) std::cout << x << "\t";
		std::cout << "\n";
	}*/
	// linear scan marked_indices and merge overlapping windows
	std::vector<ipair_t> range;
	int sz = marked_indices.size();
	if (sz) {
		int start = marked_indices.front(),
			end = start;
		for (int i = 1; i < sz; ++ i) {
			if (marked_indices[i] - end + 1 <= wl_) {
				end = marked_indices[i];
			} else {
				range.push_back(ipair_t(start, end + wl_ - 1));
				start = end = marked_indices[i];
			}
		}
		range.push_back(ipair_t(start, end + wl_ - 1));
	}
	/*{ // debug print
		std::cout << "range: \n";
		for (auto& x: range) std::cout << x.first << ", " << x.second << "\n";
	}*/

	// according to range, generate windows
	for (auto& x: range) {
		ivec_t window;
		window.assign(coordinates.begin() + x.first,
				coordinates.begin() + x.second + 1);
		windows.push_back(window);
	}
	/*{ // debug print
		std::cout << "window: \n";
		for (auto& x: windows){
			for(auto & y: x) std::cout << y << "\t";
			std::cout << "\n";
		}
	}*/
}// ProcDB::identify_low_density_windows

//
void ProcDB::output_contig(std::ofstream& ofh, int gi,
	const iivec_t& low_density_index_windows, const std::string& contig){

	strvec_t subseqs;
	for (auto& idx_window: low_density_index_windows){
		if (!idx_window.size()){
			std::cout << "ProcDB: SC failed!\n";
			exit(1);
		}
		int start = std::abs(idx_window.front()),
			len = std::abs(idx_window.back()) - start + k_;
		subseqs.push_back(contig.substr(start, len));
	}

	for (int i = 0; i < subseqs.size(); ++ i) {
		ofh << ">" << gi << "." << i << "\n" << subseqs[i] << "\n";
	}

}// ProcDB::output_contig

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


