//========================================================================
// Project     : RankAbundance
// Name        : CountKP.h
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


#ifndef COUNTKP_H_
#define COUNTKP_H_

#include "xutil.h"
#include "ReadBioFile.h"

/*
 * Generate KmerPair information from a set of fragments
 */
class CountKP{
public:
	/* @param ki_first/last: iterators to [kmerInfo_]
	 * @param kt_first/last: iterators to [kmerTaxaID_]
	 * @param k_: length of kmer
	 * @param min_kpfreq: minimum frequency of kmer pairs
	 */
	CountKP(std::vector<u64ibtuple_t>::const_iterator ki_first,
			std::vector<u64ibtuple_t>::const_iterator ki_last,
			std::vector<u64ipair_t>::const_iterator kt_first,
			std::vector<u64ipair_t>::const_iterator kt_last,
			std::back_insert_iterator<ivec_t> it_kfragfreq, int k, int min_kpfreq):
				ki_first_(ki_first), ki_last_(ki_last),
				kt_first_(kt_first), kt_last_(kt_last),
				it_kfragfreq_(it_kfragfreq),
				k_(k), min_kpfreq_(min_kpfreq) {}

	/*	@brief	Operator that count kmer pairs in fragments
	 */
	void operator ()(std::vector<uu64ituple_t>& kpcnt, const strvec_t& frags){
		int num_frag = frags.size();

		std::vector<uu64pair_t> kp;

		#pragma omp parallel
		{
			std::vector<uu64pair_t> local_kp; // local to each thread
			ivec_t local_kfragfreq;
			#pragma omp for
			for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

				u64set_t frag_kval_set;
				get_unique_frag_kmers (frag_kval_set, k_, frags[fragidx]);

				u64vec_t frag_kval_vec;
				for (auto& x: frag_kval_set) {
					int freq = -1;
					auto it = std::lower_bound(ki_first_, ki_last_,
							std::make_tuple(x, 0, false), cmp_u64ibtuple());
					if (it != ki_last_ &&
						std::get<0> (*it) == x &&
						std::get<1> (*it) > min_kpfreq_) {
						frag_kval_vec.push_back(x);
					}
				}

				/*
				int k_pos_idx = 0;
				int last_k_pos = frags[fragidx].length() - k_;

				// --- identify unique kmers in the frag ---
				u64vec_t frag_kval_vec;
				while (k_pos_idx <= last_k_pos) {
					//look for current kmer
					uint64_t val;
					std::string kmer = get_kstr(frags[fragidx].substr(k_pos_idx, k_));

					if (xny::str2ID<uint64_t>(val, kmer)) { // contain no "n"

						int freq = -1;
						auto it = std::lower_bound(ki_first_, ki_last_,
								std::make_tuple(val, 0, false), cmp_u64ibtuple());
						if (it != ki_last_ &&
							std::get<0> (*it) == val &&
							std::get<1> (*it) > min_kpfreq_) {
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
				*/

				rm_duplicated_elements (frag_kval_vec);

				// for each kmer, obtain list of taxaIDs they belong to

				int sz = frag_kval_vec.size();
				std::vector<iset_t> taxaIDs (sz);
				for (int i = 0; i < sz; ++ i) {
					auto it = std::lower_bound(kt_first_, kt_last_,
						u64ipair_t(frag_kval_vec[i], 0), cmp_u64ipair());
					for (; it != kt_last_; ++ it) {
						if (it->first == frag_kval_vec[i]){
							taxaIDs[i].insert (it->second);
						} else break;
					}
				}

				//--- debug: capture # kmers each frag contains ----
				local_kfragfreq.push_back(sz);

	//			#pragma omp critical
				/*{ // debugging
					 if (sz == 119) {
						 std::cout << "frag = \n" << frags[fragidx] << "\n";
						 for (int i = 0; i < sz; ++ i) {
							for (auto it = taxaIDs[i].begin();
									it != taxaIDs[i].end(); ++ it) {
								std::cout << (*it) << "\t";
							}
							std::cout << "\n";
						 }
						 std::cout << "debug done\n";
						 exit(1);
					 }
				}*/
				// pairwise linking kmers
				for (int i = 0; i < sz - 1; ++ i) {
					for (int j = i + 1; j < sz; ++ j) {
						iset_t combined = taxaIDs[i];
						combined.insert (taxaIDs[j].begin(), taxaIDs[j].end());
						if (combined.size() > 1) {
							local_kp.push_back(	uu64pair_t (
								std::min(frag_kval_vec[i], frag_kval_vec[j]),
								std::max(frag_kval_vec[i], frag_kval_vec[j])));
						}
					}
				}
			} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx)

			#pragma omp critical
			{
				kp.insert(kp.end(), local_kp.begin(), local_kp.end());
				std::copy(local_kfragfreq.begin(),
						local_kfragfreq.end(), it_kfragfreq_);
			}
		}

		accrue_kpcnts_(kpcnt, kp);
	};
private:
	std::vector<u64ibtuple_t>::const_iterator ki_first_, ki_last_;
	std::vector<u64ipair_t>::const_iterator kt_first_, kt_last_;
	std::back_insert_iterator<ivec_t> it_kfragfreq_;
	int k_;
	int min_kpfreq_;

	/* @brief Accrue count of kmer pairs [kp] to [kpcnt]
	 */
	void accrue_kpcnts_(std::vector<uu64ituple_t>& kpcnt,
			std::vector<uu64pair_t>& kp){
		static int debug_cnt = 0;
		//std::cout << "[accrue_kpcnts_]: " << debug_cnt ++ << "\t";
		//std::cout << "kpcnt.size() = " << kpcnt.size() << "\n";

		if (!kp.size()) return;
		std::sort (kp.begin(), kp.end());

		// convert kp to <k1, k2, cnt> structure, initialize with first elem
		std::vector<uu64ituple_t> rhs_kpcnt {
			std::make_tuple (kp[0].first, kp[0].second, 0)};

		for (auto& x: kp) { // ignore first element
			if (uu64pair_t(std::get<0>(rhs_kpcnt.back()),
					std::get<1>(rhs_kpcnt.back())) == x)
				++ std::get<2>(rhs_kpcnt.back());
			else rhs_kpcnt.push_back(std::make_tuple (x.first, x.second, 1));
		}
		kp.clear();

		std::vector<uu64ituple_t> merged_kpcnt;

		int p0 = 0, p1 = 0;
		while (p0 < kpcnt.size() && p1 < rhs_kpcnt.size()) {
			uu64pair_t kp0 {std::get<0>(kpcnt[p0]), std::get<1>(kpcnt[p0])},
					   kp1 {std::get<0>(rhs_kpcnt[p1]), std::get<1>(rhs_kpcnt[p1])};

			if (kp0 == kp1) {
				merged_kpcnt.push_back (std::make_tuple (kp0.first, kp0.second,
					std::get<2>(kpcnt[p0]) + std::get<2>(rhs_kpcnt[p1])));
				++ p0;
				++ p1;
			} else if (kp0 < kp1) {
				merged_kpcnt.push_back (kpcnt[p0]);
				++ p0;
			} else {
				merged_kpcnt.push_back (rhs_kpcnt[p1]);
				++ p1;
			}
		}
		// adding the remaining array
		if (p0 < kpcnt.size()) merged_kpcnt.insert(merged_kpcnt.end(),
				kpcnt.begin() + p0, kpcnt.end());
		else merged_kpcnt.insert(merged_kpcnt.end(), rhs_kpcnt.begin() + p1,
				rhs_kpcnt.end());

		kpcnt = merged_kpcnt;
	};
};


#endif /* COUNTKP_H_ */
