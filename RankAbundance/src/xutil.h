//========================================================================
// Project     : PrepDB
// Name        : xutil.h
// Author      : Xiao Yang
// Created on  : Oct 22, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2012)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or 
//				 guaranteed support whatsoever. The Broad Institute cannot 
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================


#ifndef XUTIL_H_
#define XUTIL_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <climits>
#include <memory>
#include <omp.h>
#include <ctime>

#include "typedef.hpp"
#include "xny/seq_manip.hpp"



// defined for a given clade
struct abund_t {
	int clade_k_num;			//# kmers selected in the clade
	int read_k_num;			//# kmers in clade supported in frags
	int unique;				// # contributed by uniquely mapped frags
	int multi;  				//# contributed by multi-mapped frags
	abund_t () {
		clade_k_num = 0;
		read_k_num = 0;
		unique = 0;
		multi = 0;
	}
};

struct clade_t {
	int k_sampled;		//# kmers selected in the clade
	int k_in_reads;		//# kmers in clade present in read data
	int u_cnt;			// # unique-aligned fragment assignment
	int m_cnt;			// # multi-aligned fragment assignment
	clade_t () {
		k_sampled = 0;
		k_in_reads = 0;
		u_cnt = 0;
		m_cnt = 0;
	}
};

struct record_t {
	uint64_t kmer;
	int start_pos;
	int node;
	std::string refname;
	record_t (uint64_t& a, int b, int c, std::string d):
		kmer (a), start_pos (b), node (c), refname (d) {};
};

struct mytime_t {
	int day;
	int hour;
	int min;
	int sec;
	mytime_t (int d, int h, int m, int s): day (d), hour (h), min(m), sec (s) {};
	void print () {
		std::cout << "\ttime = " << day << ":"
				<< hour << ":" << min << ":" << sec << "\n";
	}
};

/** Copyed from Jaz lib
 * Split a string into a list of strings.
 *  @param pat is a string separator.
 *  @param s is a string to split.
 *  @param out is a model of OutputIterator
 *  where the output should be stored.
 */
template <typename charT, typename traits, typename Alloc, typename Iter>
void split(charT pat,
	     const std::basic_string<charT, traits, Alloc>& s, Iter out) {
    unsigned int pos = 0;

    for (unsigned int i = 0; i < s.size(); ++i) {
	  if (s[i] == pat) {
	      if (i - pos > 0) {
		  *(out++) =
		      std::basic_string<charT, traits, Alloc>(s, pos, i - pos);
	      }
	      pos = i + 1;
	  }
    }

    *(out++) =
	  std::basic_string<charT, traits, Alloc>(s, pos, s.size() - pos);
} // split

/*
template <typename T, typename Comparator = std::less<T> >
void rm_duplicated_elements (std::vector<T>& array, Comparator cmp = Comparator ()) {

	std::sort (array.begin(), array.end(), cmp);
	auto it = std::unique_copy(array.begin(), array.end(), array.begin());
	array.resize(std::distance (array.begin(), it));
} */

template <typename T>
void rm_duplicated_elements (std::vector<T>& array) {

	std::sort (array.begin(), array.end());
	auto it = std::unique_copy(array.begin(), array.end(), array.begin());
	array.resize(std::distance (array.begin(), it));
}

/* @brief Removed marked elements flag = true from array
 */
template<typename T>
void rm_marked_elements (std::vector<T>& array, const bvec_t& marks) {
	if (array.size() != marks.size()) {
		std::cout << "[ERR]: rm_marked_elements SC failed\n"; exit(1);
	}
	int p0 = 0, p1;
	int sz = array.size();
	while (p0 < sz) {
		if (marks[p0]) {
			p1 = p0 + 1;
			break;
		} else ++ p0;
	}
	if (p0 >= sz) return;
	while (p1 < sz) {
		if (marks[p1] == false) {
			array[p0] = array[p1];
			++ p0;
			++ p1;
		} else ++ p1;
	}
	array.resize(p0);
} //rm_marked_elements


inline mytime_t get_time () {
	time_t t = time (0);
	struct tm * now = localtime( & t );
	return mytime_t (now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
} // get_time

inline void abording (const std::string& msg) {
	std::cout << "[EXIT]: " << msg << "\n";
	exit(1);
}

inline void warning(const std::string& msg) {
	std::cout << "\t[WARNING]: " << msg << "\n\n";
}

inline void print_strvec (const std::string& delim,
		const std::vector<std::string>& list) {
	std::cout << list.size() << " entries.\n";
	for (int i = 0; i < (int) list.size(); ++ i) {
		std::cout << delim << i << ": " << list[i] << "\n";
	}
	std::cout << "\n";
} //print_strvec

inline void print_strset (const strset_t& data){
	strset_t::const_iterator it = data.begin();
	for (; it != data.end(); ++ it) {
		std::cout << "\t" << *it;
	}
	if (data.size()) std::cout << "\n";
}



/* read a chunk of tab separated file */
inline void read_file_chunk (strvec_t& rows,
		std::ifstream& fh, int batch) {

	rows.clear();
	std::string line;
	int counter = 0;
	while (std::getline (fh, line)) {
		if (line.at(0) == '#') continue;
		if (!line.empty()) { // skip blank lines
			rows.push_back(line);
		}
		++ counter;
		if (counter >= batch) return;
	}
}

/* to be universal, treat everything as a string first */
inline void read_in_tab_file (std::vector<strvec_t>& rows,
		std::ifstream& fh, int batch) {

	rows.clear();
	std::string line;
	std::istringstream iss;
	int counter = 0;
	while (std::getline (fh, line)) {
		iss.clear();
		iss.str(line);
		strvec_t elems;
		if (line.at(0) == '#') continue;
		if (!line.empty()) { // skip blank lines
		    std::copy(std::istream_iterator<std::string>(iss),
		             std::istream_iterator<std::string>(),
		             std::back_inserter(elems));
		    rows.push_back(elems);
		}
		++ counter;
		if (counter >= batch) return;
	}
}

inline std::string get_kstr (const std::string& in) {
	std::string kmer = in;
	char c = std::toupper(kmer.at(0));
	if (c == 'G' || c == 'T') xny::rvc_str(kmer);
	return kmer;
}

template<typename T, typename Comparator = std::less<T> >
std::set<T> get_neighborset_in_pairs(const T& elem,
		const std::vector<std::pair<T, T> >& pairs,
		Comparator cmp = Comparator()) {

	std::set<T> nbs;
	auto it = std::lower_bound (pairs.begin(), pairs.end(),
			std::pair<T, T> (elem, 0), cmp);
	for (; it != pairs.end(); ++ it) {
		if (it->first == elem) nbs.insert (it->second);
		else break;
	}
	return nbs;
} // get_neighborset_in_pairs

#endif /* XUTIL_H_ */
