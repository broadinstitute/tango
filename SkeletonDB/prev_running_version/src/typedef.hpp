//========================================================================
// Project     : RankAbundance
// Name        : typedef.hpp
// Author      : Xiao Yang
// Created on  : Dec 27, 2013
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


#ifndef TYPEDEF_HPP_
#define TYPEDEF_HPP_

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <inttypes.h>
#include <tuple>

typedef std::vector<std::string> strvec_t;
typedef std::vector<int> ivec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint64_t> u64vec_t;
typedef std::vector<u64vec_t> uu64vec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<bvec_t> bbvec_t;

typedef std::pair<int, int> ipair_t;
typedef std::pair<double, double> dpair_t;
typedef std::pair<std::string, std::string> strpair_t;
typedef std::pair<int, std::string> istrpair_t;
typedef std::pair<uint64_t, int> u64ipair_t;
typedef std::pair<uint64_t, bool> u64bpair_t;
typedef std::pair<uint64_t, uint64_t> uu64pair_t;

typedef std::set<int> iset_t;
typedef std::set<std::string> strset_t;
typedef std::set<uint64_t> u64set_t;

typedef std::map<std::string, int> strimap_t;
typedef std::map<int, int> imap_t;

typedef std::pair<int, strvec_t> istrvecpair_t;

typedef std::tuple<uint64_t, int, bool> u64ibtuple_t;
typedef std::tuple<uint64_t, uint64_t, int> uu64ituple_t;
typedef std::tuple<int, int, int> iiituple_t;
typedef std::tuple<int, int, std::string> iistrtuple_t;
typedef std::tuple<int, int, int, int> i4tuple_t;
typedef std::tuple<std::string, std::string, std::string> fqtuple_t;

struct cmp_istrvecpair{
public:
	bool operator () (const  istrvecpair_t& lhs, const  istrvecpair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

struct cmp_u64ibtuple {
public:
	bool operator () (const u64ibtuple_t& lhs, const u64ibtuple_t& rhs)
	const {	return std::get<0> (lhs) < std::get<0> (rhs); }
};

struct cmp_iiituple {
public:
	bool operator () (const iiituple_t& lhs, const iiituple_t& rhs)
	const {	return std::get<0> (lhs) < std::get<0> (rhs); }
};

struct cmp_iistrtuple {
public:
	bool operator () (const iistrtuple_t& lhs, const iistrtuple_t& rhs)
	const {	return std::get<0> (lhs) < std::get<0> (rhs); }
};

struct cmp_i4tuple {
public:
	bool operator () (const i4tuple_t& lhs, const i4tuple_t& rhs)
	const {	return std::get<0> (lhs) < std::get<0> (rhs); }
};

struct cmp_u64bpair {
public:
	bool operator () (const u64bpair_t& lhs, const u64bpair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

/*
struct cmp_ipair {
public:
	bool operator () (const ipair_t& lhs, const ipair_t& rhs)
	const {	return lhs.first < rhs.first; }
};*/

struct cmp_u64ipair {
public:
	bool operator () (const u64ipair_t& lhs, const u64ipair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

struct cmp_uu64ituple {
public:
	bool operator () (const uu64ituple_t& lhs, const uu64ituple_t& rhs) const {
		if (std::get<0> (lhs) < std::get<0> (rhs)) return true;
		else if (std::get<0> (lhs) == std::get<0> (rhs)) {
			if (std::get<1> (lhs) < std::get<1> (rhs)) return true;
			else return false;
		} else return false;
	}
};

struct cmp_uu64pair {
public:
	bool operator () (const uu64pair_t& lhs, const uu64pair_t& rhs)
	const {	return lhs.first < rhs.first; }
};

#endif /* TYPEDEF_HPP_ */
