//========================================================================
// Project     : RankAbundance
// Name        : CountKmers.h
// Author      : Xiao Yang
// Created on  : Mar 5, 2014
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2014)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or 
//				 guaranteed support whatsoever. The Broad Institute cannot 
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================


#ifndef COUNTKMERS_H_
#define COUNTKMERS_H_

#include "typedef.hpp"
#include "ReadBioFile.h"

/**
 * Generate KmerInfo_ information from a set of fragments.
 */
class CountKmers{
public:
	CountKmers(u64vec_t::const_iterator first,
			u64vec_t::const_iterator last, int k):
				first_(first), last_(last), k_(k){}

	void operator()(std::vector<u64ibtuple_t>& kmerInfo,
			const strvec_t& frags) {

		int num_frag = frags.size();

		// whenever a traversed kmer x in [frags] exists in
		// [kmers] (first_, last_), x is added to kval_vec
		u64vec_t kvals;

		#pragma omp parallel
		{
			u64vec_t local_kvals; // local to each thread
			#pragma omp for
			for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

				u64set_t frag_kval_set;

				get_unique_frag_kmers (frag_kval_set, k_, frags[fragidx]);

				for (auto& val: frag_kval_set) {
					if (std::binary_search(first_, last_, val)){
						local_kvals.push_back(val);
					}
				}
			} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

			#pragma omp critical
			{
				kvals.insert(kvals.end(), local_kvals.begin(), local_kvals.end());
			}
		}

		accrue_kcnts_(kvals, kmerInfo);
	};

private:
	u64vec_t::const_iterator first_, last_;
	int k_;

	void accrue_kcnts_(u64vec_t& kvals, std::vector<u64ibtuple_t>& kmerInfo) {
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

		// merge <kmer, cnt> to [kmerInfo]
		std::vector<u64ibtuple_t> merged_ki;

		int p0 = 0, p1 = 0;
		while (p0 < kmerInfo.size() && p1 < kcnt.size()) {
			uint64_t val0 = std::get<0>(kmerInfo[p0]), val1 = kcnt[p1].first;
			if (val0 == val1) {
				merged_ki.push_back (std::make_tuple (kcnt[p1].first,
						kcnt[p1].second + std::get<1> (kmerInfo[p0]), false));
				++ p0;
				++ p1;
			} else if (val0 < val1) {
				merged_ki.push_back (kmerInfo[p0]);
				++ p0;
			} else {
				merged_ki.push_back (std::make_tuple (kcnt[p1].first,
						kcnt[p1].second, false));
				++ p1;
			}
		}

		// adding the left-overs from either [kmerInfo] or [kcnt]
		if (p0 < kmerInfo.size()) merged_ki.insert(merged_ki.end(),
				kmerInfo.begin() + p0, kmerInfo.end());
		else {
			for (; p1 < kcnt.size(); ++ p1) {
				merged_ki.push_back (std::make_tuple (kcnt[p1].first,
						kcnt[p1].second, false));
			}
		}

		kmerInfo = merged_ki;
	}
};

#endif /* COUNTKMERS_H_ */
