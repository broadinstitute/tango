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

void proc_db (const std::string& odir, const std::string& ogilevel,
	const std::string& idir, const ivec_t& gis,
	const std::vector<ipair_t>& vec_gi_node,
	const std::vector<istrpair_t>& tree,	bool silent) {

	// obtain input files
	DIR* dir = opendir (idir.c_str());
	if (!dir) abording ("cannot open " + idir);
	struct dirent *dirp;

	strvec_t input_files;
 	std::vector<ipair_t> vec_gi_len;
	while ((dirp = readdir(dir))) {
		std::string filename = dirp->d_name;
		unsigned pos = filename.find_last_of('.');
		std::string suffix = filename.substr (pos + 1);
		std::string path = idir + "/" + filename;
		if (suffix.compare("fa") == 0 || suffix.compare("fasta") == 0) {
			if (!silent) std::cout << "\tprocess " << path << "\n";
			get_gi_len_fa (vec_gi_len, path);
			input_files.push_back(path);
		}
	}

	if (!silent) std::cout << "\tnum of sequences: "
			<< vec_gi_len.size() << "\n";
	closedir (dir);

	if (vec_gi_node.empty()) abording ("vec_gi_node is empty");

	if (tree.empty()) abording ("tree is empty");

	// clustering the gi numbers
	int counter = 0;

	std::map<int, ivec_t> cluster; // tree node --> indices of vec_gi_len

	for (size_t i = 0; i < vec_gi_len.size(); ++ i) {

		int gi = vec_gi_len[i].first;

		// undocumented gi put to cluster[-1]
		if (!std::binary_search(gis.begin(), gis.end(), gi)) {
			auto it_cls = cluster.find(-1);
			if (it_cls != cluster.end()) it_cls->second.push_back(i);
			else cluster[-1] = ivec_t(1, i);
		} else { // documented gi
			auto it_gn =	std::upper_bound(vec_gi_node.begin(), vec_gi_node.end(),
					ipair_t(gi, -1), cmp_ipair());
			if (it_gn != vec_gi_node.end()) {
				int idx = it_gn - vec_gi_node.begin() - 1;
				// debug
				/*{
					std::cout << "gi = " << gi << "\n";
					std::cout << vec_gi_node[idx].first << "\t" <<
						vec_gi_node[idx].second << "\n";
				}*/
				int node = vec_gi_node[idx].second;
				if (node >= tree.size()) abording ("proc_db SC0 failed");
				int parentID = tree[node].first; // to be used as cluster ID
				auto it_cls = cluster.find(parentID);
				if (it_cls != cluster.end()) it_cls->second.push_back(i);
				else cluster[parentID] = ivec_t(1, i);
			} else {
				std::cout << "\nErr in proc_db: cannot find gi: " << gi << "\n";
				exit(1);
			}
		}
	}

	// print out clusters
	if (!silent) std::cout << "\tnum of clusters: " << cluster.size() - 1
			<< "\n";
	/*
	for (auto& cls: cluster) {
		std::cout << cls.second.size() << "\n";
	}
	std::cout << "\n"; */

	int max_file_sz = 500000000; // 500MB
	imap_t gi_fileID; // gi to fileID (starting from 0)
	strvec_t output_files;
	std::vector<iset_t> giset; // clsId --> {gi}s  used for closing file handle
	int fileID = 0;
	std::ofstream ofh_gilevel;
	if (!ogilevel.empty()) {
		ofh_gilevel.open (ogilevel.c_str(), std::ofstream::out); // version 2
	}

	for (auto& cls: cluster) {
		if (cls.first == -1) {
			if (!silent) std::cout << "\t" << cls.second.size()
					<< " (likely) outdated records [ignored]\n";
			continue;
		}
		int splits = 0;
		output_files.push_back(gen_filename (odir, cls.first, splits));

		int total_len = 0;
		for (auto& idx: cls.second) {
			if (idx >= vec_gi_len.size()) abording ("proc_db SC1 failed");
			int gi = vec_gi_len[idx].first;
			total_len += vec_gi_len[idx].second;
			gi_fileID[gi] = fileID;
			if (ofh_gilevel.is_open()) {
				ofh_gilevel << gi << "\t" << cls.first << "\t";
				// identify the taxaRank of the cls.first
				if (cls.first >= tree.size()) {
					abording ("proc_db SC2.failed");
				}
				ofh_gilevel << tree[cls.first].second << "\n";
			}

			if (giset.size() < fileID + 1) {
				iset_t tmp;
				tmp.insert (gi);
				giset.push_back(tmp);
			} else giset[fileID].insert(gi);
			if (total_len > max_file_sz) {
				++ fileID;
				++ splits;
				output_files.push_back(gen_filename (odir, cls.first, splits));
				total_len = 0;
			}
		}
		if (total_len != 0) ++ fileID;
	}

	if (!ogilevel.empty()) ofh_gilevel.close();

	if (!silent) std::cout << "\treorganize cluster by max size, num cls: "
			<< giset.size() << "\n";


	if (!silent) std::cout << "\tOutput\n";

	/*strvec_t output_files (giset.size());
	int fidx = 0;
	for (auto& f: output_files) {
		std::string fID = std::to_string(fidx);
		++ fidx;
		f = odir + "/";
		f += fID;
		f += ".fa";
	}*/

	/*
	for (auto& f: filenames) std::cout << f << ", ";
	std::cout << "\n"; */
	if (!odir.empty()) output_db (output_files, giset, gi_fileID, input_files);
} // proc_db

std::string gen_filename (const std::string& dir, int id, int subid) {
	std::string name = dir + "/";
	std::string fID = std::to_string(id);
	if (id == -1) fID = "u-1";
	name += fID;
	name += ".";
	name += std::to_string(subid);
	name += ".fa";
	return name;
}

void output_db (const strvec_t& output_files, std::vector<iset_t>& giset,
		const imap_t& gi_clsId, const strvec_t& input_files) {

	/* doesn't work on recent mac
	struct rlimit limit;
	 limit.rlim_cur = 65535;
	 limit.rlim_max = 65535;
	 if (setrlimit(RLIMIT_NOFILE, &limit) != 0) {
	    printf("setrlimit() failed with errno=%d\n", errno);
	    exit(1);
	 }*/

	int num_files = output_files.size();
//	std::vector<std::unique_ptr<std::ofstream> > FHs (num_files); // version 1
	std::vector<std::ofstream> FHs(num_files); //  version 2
	bvec_t is_open (num_files, false);
	for (size_t i = 0; i < input_files.size(); ++ i) {
		std::ifstream ifh (input_files[i].c_str());
		std::string line;
		int gi = -1;
		int fileID = -1;
		while (std::getline(ifh, line)) {
			if (line.at(0) == '>') {
				if (fileID != -1 && giset[fileID].empty()) FHs[fileID].close(); // close file handle
				get_gi (gi, line);
				if (gi != -1) {
					auto it = gi_clsId.find(gi);
					if (it == gi_clsId.end()) {
						//abording ("output_db SC0 fail");
						fileID = -1;  // guarantee the unrecorded gi are not to be written
					}
					else {
						fileID = it->second;
						if (fileID >= num_files) abording ("output_db SC1 fail");
						if (!is_open[fileID]) {
//							FHs[fileID] = std::unique_ptr<std::ofstream>
//								(new std::ofstream (output_files[fileID].c_str())); // version 1
							FHs[fileID].open (output_files[fileID].c_str(), std::ofstream::out); // version 2
							if (!FHs[fileID]) {
								abording ("Failed to open " + output_files[fileID]);
							}
							is_open[fileID] = true;
						}
//						*FHs[fileID] << line << "\n"; // version 1
						FHs[fileID] << line << "\n"; // version 2
						giset[fileID].erase(gi);
					}

				} else {
					warning ("gi not found in entry: " + line);
					fileID = -1;
				}

			} else {
//				if (fileID != -1) *FHs[fileID] << line << "\n"; // version 1
				if (fileID != -1) FHs[fileID] << line << "\n"; // version 2

			}
		}
		ifh.close();
	}
} //

/* @ Parse fasta file and get gi length pair for each entry
 */
void get_gi_len_fa (std::vector<ipair_t>& vec_gi_len, const std::string& file) {
	std::ifstream fh (file.c_str());
	if (!fh.good()) abording ("get_gi_len_fa: can't open " + file);

	std::string line;
	int counter = 0;
	int len = 0;
	int gi = -1;
	while (std::getline (fh, line)) {
		if (line.at(0) == '>') {
			if (len == 0) {
				get_gi (gi, line);
				continue;
			} else {
				vec_gi_len.push_back(ipair_t(gi, len));
				len = 0;
				get_gi (gi, line);
			}
			++ counter;
		} else len += line.length();
	}
	if (len > 0) {
		vec_gi_len.push_back(ipair_t(gi, len));
		++ counter;
	}

	// debug print
	/*
	std::cout << "num of entries " << counter << "\n";
	counter = 0;
	for (auto& x: vec_gi_len) {
		std::cout << x.first << ": " << x.second << "\n";
		if (++counter > 20) exit(1);
	}*/

	fh.close();
} // get_gi_len_fa

/* @ get gi number from header  in the format of ">gi|number|blabla"
 */
void get_gi (int& gi, const std::string& header) {
	size_t start = header.find_first_of ('|');
	size_t end = header.find_first_of ('|', start + 1);
	if (end > start + 1) gi = atoi(header.substr(start + 1, end - start - 1).c_str());
//	std::cout << header << "\n" << start << "\n" << end << "\n" << gi <<"\n"; exit (1);
} // get_gi

