#ifndef _NEWICK_H_
#define _NEWICK_H_

#include <vector>
#include <string>

#include "phylo.h"

using namespace std;

/* Parse a single tree in Newick format */
PhyEdgeSet NewickParse(const string &str, map<string,unsigned int> &name_to_pos);

/* Use a single tree to assign each leaf label an integer ID to be used for all
 * future parsing */
map<string,unsigned int> AssignLeafLabels(const string &str); 

#endif /* _NEWICK_H_ */

