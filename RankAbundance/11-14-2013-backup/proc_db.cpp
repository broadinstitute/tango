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


#include "proc_db.h"

void proc_db (std::vector<u64ipair_t>& record, std::vector<strvec_t>& input_files,
		strvec_t& nodeNames, const double& perc, int k,
		const std::string& idir, bool silent) {

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
				nodeNames.push_back(prefix);
			}
		}
	}
	int num_node = nodeNames.size();
	if (!silent) std::cout << "\t" << num_node
					<< " num of nodes to process\n";
	closedir (dir);

	if (!input_files.size()) return;

	//srand (time(NULL));

	// process each file
	/*std::vector<std::pair<uint64_t, int> > record;

	#pragma omp parallel
	{
		std::vector<std::pair<uint64_t, int> > private_record;
		#pragma omp for schedule (dynamic)
		for (int i = 0; i < num_files; ++ i) {
			//if (!silent) std::cout << i << "\n";
			u64vec_t kmers;
			sample_kmers (kmers, k, perc, input_files[i]);
			std::vector<std::pair<uint64_t, int> > local (kmers.size());
			for (int j = 0; j < kmers.size(); ++ j) {
				local[j] = std::pair<uint64_t, int> (kmers[j], i);
			}
			private_record.insert(private_record.end(), local.begin(), local.end());
		}
		#pragma omp critical
		{
			record.insert(record.end(), private_record.begin(),
					private_record.end());
			//std::cout << "\tsize of record: " << record.size() << "\n";
		}
	}
	if (!silent) std::cout << "\tsize of record: " << record.size() << "\n";
	*/

	#pragma omp parallel shared (input_files)
	{
		std::vector<u64ipair_t> private_record;
		//std::vector<ipair_t> private_range;
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
			record.insert(record.end(), private_record.begin(),
					private_record.end());
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

	if (!silent) std::cout << "\tfinal size of record: " << record.size() << "\n";

	//for (int i = 1; i < range.size(); ++ i) range[i].first += range[i - 1].first;

} // proc_db


void sample_kmers (u64vec_t& kmers, int k, const double& perc,
		const std::string& file) {

	//std::cout << file << "\n";
	std::ifstream fh (file.c_str());
	if (!fh.good()) abording ("sample_kmers: can't open " + file);

	std::string contig;
	std::string line;

	while (std::getline (fh, line)) {
		if (line.at(0) == '>') {
			proc_contig (kmers, k, perc, contig);
			contig.clear();
		} else contig += line;
	}
	proc_contig (kmers, k, perc, contig);
	fh.close();

	/*std::sort(kmers.begin(), kmers.end());
	auto it = std::unique_copy (kmers.begin(), kmers.end(), kmers.begin());
	kmers.resize(std::distance(kmers.begin(), it));
	*/
} //sample_kmers

void proc_contig (u64vec_t& kmers, int k, const double& perc,
		const std::string& contig) {
	int len = contig.length();
	int num = (len - k + 1) * perc/100;
	if (num <= 0) return;
	int step = (len - k + 1) / num;
	for (int i = 0; i < num; ++ i) {
		uint64_t ID_f, ID_r;
		std::string kmer = contig.substr(i*step, k);
		if (xny::str2ID(ID_f, kmer)) {
			std::string rvc = xny::get_rvc_str(kmer);
			xny::str2ID(ID_r, rvc);
			kmers.push_back(std::min(ID_f, ID_r));
		}
	}
}

/*
void proc_contig (u64vec_t& kmers, int k, const double& perc,
		const std::string& contig) {
	int len = contig.length();
	int base = 1000/perc;
	if (len - k + 1 >= 100/perc) {
		// get kmer
//		ivec_t pos (len - k + 1);
//		for (int i = 0; i < len - k + 1; ++ i) pos[i] = i;
//		std::random_shuffle(pos.begin(), pos.end());
//		int num_sample = (len - k + 1) * perc / 100;
//		for (int i = 0; i < num_sample; ++ i) {
			// convert kmer in contig[pos[i]]
//			uint64_t ID;
//			if (	str2ID (ID, contig.substr(pos[i], k))) kmers.push_back(ID);
//			if (str2ID (ID, contig.c_str() + pos[i], k)) kmers.push_back(ID);
//		}
		for (int i = 0; i < len - k + 1; ++ i) {
			if (rand () % base + 1 <= 10) {
				uint64_t ID;
				if (	str2ID (ID, contig.substr(i, k))) kmers.push_back(ID);
			}
		}
	}
}*/
