//========================================================================
// Project     : RankAbundance
// Name        : AssignFrag.h
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


#ifndef ASSIGNFRAG_H_
#define ASSIGNFRAG_H_

#include "ReadBioFile.h"


/* Assign a set of fragments to taxa
 */
class AssignFrag{
public:
	/* @param kt_first/last: iterators to [kmerTaxaID_]
	 * @param k_: length of kmer
	 */
	AssignFrag(std::vector<u64ipair_t>::const_iterator kt_first,
			std::vector<u64ipair_t>::const_iterator kt_last, int k):
			kt_first_(kt_first), kt_last_(kt_last), k_(k) {}

	/* @brief	Operator that assigns fragments to taxa
	 */
	void operator ()(ivec_t& taxacnts, const strvec_t& frags) {

		int num_frag = frags.size();
		int num_taxa = taxacnts.size();

		#pragma omp parallel
		{
			ivec_t local_taxacnts (num_taxa, 0);
			#pragma omp for
			for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

				imap_t fileID_cnt;
	 			int total_fIDs = 0;

				u64set_t frag_kval_set;
				get_unique_frag_kmers (frag_kval_set, k_, frags[fragidx]);

				for (auto& x: frag_kval_set) {
					auto it = std::lower_bound(kt_first_, kt_last_,
							u64ipair_t(x, 0), cmp_u64ipair());
					for (; it != kt_last_; ++ it) {
						if (it->first == x){
							int fID = it->second;
							if (fID >= num_taxa)
								abording ("AssignFrag::operator() SC failed");
							auto it_fc = fileID_cnt.find(fID);
							if (it_fc != fileID_cnt.end()) ++ it_fc->second;
							else fileID_cnt[fID] = 1;
							++ total_fIDs;
						} else break;
					}
				}

				for (auto& x: fileID_cnt) {
					local_taxacnts[x.first] += (x.second * 1000) / total_fIDs;
				}
			} // for (int fragidx = 0; fragidx < num_frag; ++ fragidx) {

			#pragma omp critical
			{
				for (int i = 0; i < num_taxa; ++ i) {
					taxacnts[i] += local_taxacnts[i];
				}
			}
		}
	};

private:
	std::vector<u64ipair_t>::const_iterator kt_first_, kt_last_;
	int k_;
};



#endif /* ASSIGNFRAG_H_ */
