//============================================================================
// Name        : PrepDB.cpp
// Author      : Xiao Yang
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

/*
 * Note: need to check the number of files allowed to be opened for
 * the system you are running
 * 1) MAC
 * 	a) to check: launchctl limit
 * 	b) to change: launchctl limit maxfiles 1000
 * 2) linux
 * 	a) to check: ulimit -a
 * 	b) to change: ulimit -n 1000
 */

#include "Parameter.h"
#include "proc_gi_node.h"
#include "proc_tax_tree.h"
#include "proc_db.h"

int main (int argc, char** argv){

	Parameter myPara (argc, argv);
	int batch = 1000000;

	omp_set_num_threads(myPara.p);

	/* parse gi to node information from raw dmp file to get all recorded gis */
	ivec_t gis;
	if (!myPara.iginoderaw.empty()) {
		if (!myPara.silent) {
			std::cout << "Obtain all gis from " << myPara.iginoderaw << "\n";
		}
		get_gis (gis, myPara.iginoderaw, batch, myPara.silent);
	}

	/* parse gi to node information */
	std::vector<ipair_t> vec_gi_node;
	if (!myPara.iginode.empty()) {
		if (!myPara.oginode.empty()) {
			if (!myPara.silent) {
				std::cout << "Parse " << myPara.iginode
				<< " to generate giRange->node to " << myPara.oginode << "\n";
			}
		} else if (!myPara.silent) {
			std::cout << "Read giRange->node from " << myPara.iginode << "\n";
		}

		get_giRange_node (vec_gi_node, myPara.iginode, myPara.oginode,
				batch, myPara.silent);

		if (!myPara.silent) std::cout << "\tdone\n";
	}


	/* parse tree node to parent information */
	std::vector<istrpair_t> tree;
	if (!myPara.itree.empty()) {
		if (!myPara.silent) std::cout << "Parse " << myPara.itree
									  << " to get tree structure\n";

		proc_tax_tree (tree, myPara.itree, myPara.level,
				batch, myPara.silent);
	}

	/* recreate database */
	if (!myPara.idbdir.empty() && (!myPara.odbdir.empty() || !myPara.ogilevel.empty())) {
		if (!myPara.silent) std::cout << "Analyze database " << myPara.idbdir << "\n";
		if (!myPara.odbdir.empty()) std::cout << "\t ouput as " << myPara.odbdir;
		if (!myPara.ogilevel.empty()) std::cout << "\t create gi->level: " << myPara.ogilevel;
		std::cout << "\n";

		proc_db (myPara.odbdir, myPara.ogilevel, myPara.idbdir,
				gis, vec_gi_node, tree, myPara.silent);

	}

	return (EXIT_SUCCESS);
}
