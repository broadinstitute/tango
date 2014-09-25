//========================================================================
// Project     : PrepDB
// Name        : proc_tax_tree.cpp
// Author      : Xiao Yang
// Created on  : Oct 23, 2013
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

#include "proc_tax_tree.h"

/* @brief	Process taxonomic tree information
 *
 */
void proc_tax_tree (std::vector<istrpair_t>& tree, const std::string& ifile,
		const std::string& level, int batch, bool silent) {

	if (!silent) std::cout << "\tGet taxonomy types\n";
	strset_t allRanks;
	int maxNodeID = -1;
	get_tax_types (maxNodeID, allRanks, ifile, batch);

	if (!silent) std::cout << "\tGen taxonomy tree\n";
	// vector of ([parentID]\t[*this taxonomy rank])
	tree = std::vector<istrpair_t> (maxNodeID + 1, istrpair_t(-1, std::string()));
	gen_tree (tree, ifile, batch);

	// debug print
	//std::cout << tree[367775].first << ", " << tree[367775].second << "\n";
	//std::cout << tree[552467].first << ", " << tree[552467].second << "\n";


	if (!silent) std::cout << "\tOrdering taxonomy types\n";
	std::vector<strvec_t> orderedRanks; // all linear paths
	strvec_t ranks (allRanks.begin(), allRanks.end());
	rank_ordering (orderedRanks, ranks, tree);

	// fill each entry of the 2x2 array to specify pairwise relationships
	int sz = ranks.size();
	bbvec_t array (sz, bvec_t(sz, false));

	strimap_t rank_arrayidx; // map rank_name to index in the array
	for (size_t i = 0; i < sz; ++ i) rank_arrayidx[ranks[i]] = i;

	for (auto& x: orderedRanks) {
		// collect all indices that < idx_i by checking  column
		ivec_t init_le_list;
		int idx_init = rank_arrayidx[x.front()];
		for (size_t c = 0; c < sz; ++ c) {
			if (array[c][idx_init] == true) init_le_list.push_back(c);
		}
		ivec_t init_ge_list;
		idx_init = rank_arrayidx[x.back()];
		for (size_t c = 0; c < sz; ++ c) {
			if (array[idx_init][c] == true) init_ge_list.push_back(c);
		}


		for (size_t i = 0; i < x.size(); ++ i) {
			int idx_i = rank_arrayidx[x[i]];

			for (size_t j = i + 1; j < x.size(); ++ j) {
				int idx_j = rank_arrayidx[x[j]];
				array[idx_i][idx_j] = true;
				for (auto& y: init_le_list) array[y][idx_j] = true;
			}
			for (auto& y: init_ge_list) array[idx_i][y] = true;
		}
	}

	// print out array
	/*
	int cnt = 0;
	for (auto& x: ranks) std::cout << cnt ++ << ": " << x << "\t";
	std::cout << "\n";
	cnt = 0;
	for (auto& x: array) {
		std::cout << cnt++ << " : ";
		for (const auto& y: x) std::cout << y << "  ";
		std::cout << std::endl;
	}*/
	///* print out paths
	if (!silent) {
		std::cout << "\n\tLists of ordered ranks:\n\t\t";
		int num_list = 0;
		for (auto& x: orderedRanks) {
			std::cout << "\n\t" << num_list ++  << "\n\t\t";
			int cnt = 0;
			for (auto& y: x) {
				++ cnt;
				std::cout << y << ", ";
				if (cnt % 10 == 0) std::cout << "\n\t\t";
			}
			std::cout << std::endl;
		}
	}
	if (!silent) std::cout << "\n\tGenerate cluster at " << level << " rank\n";
	int idx_level;
	strimap_t::iterator it_ra = rank_arrayidx.find(level);
	if (it_ra != rank_arrayidx.end()) idx_level = it_ra->second;
	else abording ("proc_tax_tree SC0: Cannot find " + level);

	// given level, update tree so that the parent for each node is >= level
	for (size_t i = 0; i < tree.size(); ++ i) {
		// debug
		/*if (i == 367775) {
			std::cout << "here\n";
		}*/

		int node = i;
		while (tree[node].first != -1 && tree[node].first != node) {
			int idx_rank;
			strimap_t::iterator it = rank_arrayidx.find (tree[node].second);
			if (it != rank_arrayidx.end()) idx_rank = it->second;
			else abording ("proc_tax_tree SC1: Cannot find " + tree[node].second);
			// compare rank with level
			if (array[idx_rank][idx_level] == true) {
				node = tree[node].first;
			}
			else break;
		}
		tree[i].first = node;
	}

	// clear out content to save memory
	//for (auto& x: tree) x.second.clear();

	//debug print out tree
	/*
	int cnt = 0;
	for (auto& x: tree) {
		if (x.first != -1) {
			std::cout << cnt << ": " << x.first << "  " << tree[x.first].second << "\n";
		}
		++ cnt;
	}*/
} // proc_tax_tree

/* @brief	Ordering taxomonic ranks according to tree
 */
void rank_ordering (std::vector<strvec_t>& orderedRanks,
		const strvec_t& ranks, const std::vector<istrpair_t>& tree) {

	int sz = ranks.size();

	//debug print out ranks
	/*
	for (size_t i = 0; i < ranks.size(); ++ i) std::cout << i
		<< ":" << ranks[i] << "\t";
	std::cout << "\n"; */

	// 2x2 array registering < information between ranks
	bbvec_t array (sz, bvec_t(sz, false));

	// map rank_name to index in the array
	strimap_t rank_arrayidx;
	for (size_t i = 0; i < sz; ++ i) rank_arrayidx[ranks[i]] = i;

	// fill in [array]: go through tree structure and generate < relationships
	// among different ranks
	// e.g. species < genus, genus < superfamily etc.

 	for (auto& x: tree) {
		int pID = x.first;
		if (pID == -1) continue;

		// search until find parent to be different than "no" (no rank)
		while (tree[pID].second.compare("no") == 0 && pID != tree[pID].first) {
			pID = tree[pID].first;
		}

		int idx0, idx1;

		if (x.second != tree[pID].second && tree[pID].second.compare("no") != 0) {
			strimap_t::iterator it = rank_arrayidx.find(x.second);
			if (it != rank_arrayidx.end()) idx0 = it->second;
			else abording ("rank_ordering: SC failed 0");
			it = rank_arrayidx.find(tree[pID].second);
			if (it != rank_arrayidx.end()) idx1 = it->second;
			else abording ("rank_ordering: SC failed 1");

			// only set "<" as "true"
			if (array[idx1][idx0] == true) {
				std::cout << "\n" << idx1 << ", " << idx0 << "\n";
				abording ("rank_ordering: contradicting events happening");
			}

			array[idx0][idx1] = true;
		}
 	} // for (auto& x: tree)

	/// debug print out the array
	/*for (auto& x: array) {
		for (const auto& y: x) std::cout << y << "  ";
		std::cout << std::endl;
	}*/

 	// given 2x2 array, order the ranks
 	generate_order (orderedRanks, ranks, array);


} // rank_ordering

void generate_order (std::vector<strvec_t>& orderedRanks,
		const strvec_t& ranks, const bbvec_t& array) {
	int sz = ranks.size();
	// start to order according to array
	struct vertex{
		std::string name;
		int idx;
		iset_t parents; // indices wrt ranks
		iset_t children; // indices wrt ranks
	};
	std::vector<vertex> g (sz), g_copy;
	for (size_t i = 0; i < sz; ++ i) {
		g[i].name = ranks[i];
		g[i].idx = i;
	}

	// generate graph from array by keeping as few edges as possible
	for (size_t i = 0; i < sz; ++ i) {
		for (size_t j = 0; j < sz; ++ j) {
			if (array[i][j] == true) { // i < j

				// to add a parent j make sure j is not the parent of
				// any of i's parent x
				bool is_found = false;
				for (auto& x: g[i].parents) {
					if (g[x].parents.count(j)) {
						is_found = true;
						break;
					}
				}
				if (!is_found) {
					g[i].parents.insert(j);
					g[j].children.insert(i);

					// if j share a child x with i, remove link b/t j and x
					for (auto& x: g[i].children) {
						if (g[x].parents.count(j)) {
							g[x].parents.erase(j);
							g[j].children.erase(x);
						}
					}
					// also if j share a parent x with i, remove link b/t i and x
					for (auto& x: g[j].parents) {
						if (g[i].parents.count(x)) {
							g[x].children.erase(i);
							g[i].parents.erase(x);
						}
					}
				}
			}
		}
	}

	// debug print the graph
	/*
	for (auto& x: g) {
		std::cout << x.idx << ": " << x.name << " (";
		for (auto& y: x.parents) std::cout << y << "\t";
		std::cout << ")\t (";
		for (auto& y: x.children) std::cout << y << "\t";
		std::cout << ")\n";
	}*/

	// remove bubbles
	// iterate: for each node without any parents add itself and its
	// parent list to all of its children (by accumulating the counts
	// of parents). Remove its relationships w/ all its children;
	// stop until no changes made.
	// Then for each node, remove any parent if it occurs more than once.
 	g_copy = g;
	iivec_t parents (sz, ivec_t (sz, 0));
	bool is_action_taken = true;
	while (is_action_taken) {
		is_action_taken = false;
		for (auto& x: g_copy) {
			if (x.parents.empty()) {
				for (auto& c: x.children) {
					++ parents[c][x.idx];
					for (size_t i = 0; i < sz; ++ i) {
						parents[c][i] += parents[x.idx][i];
					}
					is_action_taken = true;
					g_copy[c].parents.erase(x.idx);
				}
				x.children.clear();
			}
		}
	}

	// remove non-immediate parents from each node
	for (auto& x: g) {
		iset_t to_remove;
		for (auto& y: x.parents) {
			if (parents[x.idx][y] > 1) to_remove.insert(y);
		}
		for (auto& y: to_remove) {
			x.parents.erase(y);
			g[y].children.erase(x.idx);
		}
	}

	// debug print out simplified parent-children relationship
	/*for (auto& x: g) {
		std::cout << x.idx << ": " << x.name << " (";
		for (auto& y: x.parents) std::cout << y << "\t";
		std::cout << ")\t (";
		for (auto& y: x.children) std::cout << y << "\t";
		std::cout << ")\n";
	}*/

	// finalize ordering by traversing the graph
	int num_visited = 0;
	bvec_t visits (sz, false);
	while (num_visited < sz) {
		strvec_t path;
		for (auto& x: g) {
			// find node with empty child but not empty parent
			if (x.children.empty() && !x.parents.empty()) {
				strvec_t path;
				int node = x.idx;
				path.push_back(x.name);
				while (1) {
					if (visits[node] == false) {
						++ num_visited;
						visits[node] = true;
 					}
					// can find a not yet visited parent
					bool is_found = false;
					for (auto& y: g[node].parents) {
						if (visits[y] == false) {
							is_found = true;
							path.push_back(g[y].name);
							g[node].parents.erase(y);
							g[y].children.erase(node);
							node = y;
							break;
						}
					}
					if (!is_found) {
						// just to include the parent if any
						if (!g[node].parents.empty()) {
							int pID = *g[node].parents.begin();
							path.push_back(g[pID].name);
						}
						break;
					}
				}

				// store path
				if (path.size() > 1) orderedRanks.push_back(path);

				break;
			}
		}
	}

} //generate_order


void gen_tree (std::vector<istrpair_t>& tree, const std::string& ifile,
		int batch) {
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);
	int num_row = 0;
	while (ifh.good()) {
		std::vector<strvec_t> rows;
		read_in_tab_file (rows, ifh, batch);
		num_row += rows.size();
		batch_gen_tree (tree, rows);
	}
	ifh.close();
}

void batch_gen_tree (std::vector<istrpair_t>& tree, const std::vector<strvec_t>& rows) {
	int tree_sz = tree.size();
	if (rows.size() < 1) return;
	for (auto& x: rows) {
		if (x.size() < 5) {
			warning ("\tbatch_proc_tax_tree: entry number < 5\n");
			continue;
		}
		int nodeID = std::atoi(x[0].c_str());

		int parentID = std::atoi(x[2].c_str());
		std::string taxRank = x[4];
		if (nodeID >= tree_sz) abording ("batch_gen_tree: SC failure");
		tree[nodeID] = istrpair_t (parentID, taxRank);
	}
}

void get_tax_types (int& maxNodeID, strset_t& allRanks,
		const std::string& ifile, int batch) {
	std::ifstream ifh (ifile.c_str());
	if (!ifh.good()) abording ("Can't open " + ifile);
	int num_row = 0;
	while (ifh.good()) {
		std::vector<strvec_t> rows;
		read_in_tab_file (rows, ifh, batch);
		num_row += rows.size();
		batch_get_tax_types (maxNodeID, allRanks, rows);
	}
	ifh.close();
}

void batch_get_tax_types (int& maxNodeID, strset_t& allTaxRanks,
		const std::vector<strvec_t>& rows) {
	if (rows.size() < 1) return;
	for (auto& x: rows) {
		if (x.size() < 5) {
			warning ("\tbatch_proc_tax_tree: entry number < 5\n");
			continue;
		}
		maxNodeID = std::max(maxNodeID, std::atoi(x[0].c_str()));
		//int parentID = std::atoi(x[2].c_str());
		std::string taxRank = x[4];
		allTaxRanks.insert(taxRank);
	}
}
