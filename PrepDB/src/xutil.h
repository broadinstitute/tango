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
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <string>
#include <climits>
#include <dirent.h>
#include <memory>
#include <omp.h>
#include <stdlib.h>

typedef std::vector<std::string> strvec_t;
typedef std::vector<int> ivec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<bvec_t> bbvec_t;

typedef std::pair<int, int> ipair_t;
typedef std::pair<double, double> dpair_t;
typedef std::pair<std::string, std::string> strpair_t;
typedef std::pair<int, std::string> istrpair_t;

typedef std::set<int> iset_t;
typedef std::set<std::string> strset_t;

typedef std::map<std::string, int> strimap_t;
typedef std::map<int, int> imap_t;

typedef std::tuple<int, int, int> iiitpl_t;

template <typename T>
inline std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}
template <typename T>
inline T string_to (const std::string& st)
{
	T t;
	std::stringstream (st) >> t;
	return t;
}

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

#endif /* XUTIL_H_ */
