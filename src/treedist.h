/* 
 * Polynomial Time Distance Algorithm for Binary Rooted Phylogenetic Trees
 *
 * Note: Algorithm based on presentation by Scott Provan / Megan Owen,
 *       "Computing Geodesic Distances in Tree Space in Polynomial Time"
 * 
 * John Chakerian
 * chakj@stanford.edu
 *
 */

#ifndef __TREEDIST_H_
#define __TREEDIST_H_

#include <string>

#include "phylo.h"
#include "newick.h"

double ConeDistance(PhyEdgeSet &a, PhyEdgeSet &b);
double TreeDistance(PhyEdgeSet &a, PhyEdgeSet &b, stl_bool* incompatible);
double DisjointTreeDistance(PhyEdgeSet &a, PhyEdgeSet &b, stl_bool* incompatible, unsigned int max_id);
void ClampNegativeWeights(PhyEdgeSet &a);

enum maxflow_enum { LP, EK, EKF, KOL, PR };
extern maxflow_enum mf_method;
enum metric_enum { GEODESIC, CONEPATH };
extern metric_enum metric_method;

string tolower_str(string in);

#endif
