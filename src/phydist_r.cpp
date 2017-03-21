/* 
 * Polynomial Time Distance Algorithm for Binary Rooted Phylogenetic Trees
 *
 * Note: Algorithm based on presentation by Scott Provan / Megan Owen,
 *       "Computing Geodesic Distances in Tree Space in Polynomial Time"
 * 
 * John Chakerian
 * chakj@stanford.edu
 *
 * Last edited: Feb 28 2010
 *
 */

#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>

#include "treedist.h"

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> // for fmax2

using namespace std;

int compute_phylo_distance_matrix(vector<string> trees_in, bool verbose, double *m);
void build_tree_list(std::vector<string> &trees_in, std::vector<PhyEdgeSet> &trees, bool verbose);
double gromov_graycode(double *m, size_t n, double* deltas, int scale);

inline bool edge_less_than(const PhyEdge &a, const PhyEdge &b)
{
    for(unsigned int i = 0; i < a.split.size(); i++)
    {
        if(a.split[i] < b.split[i])
            return true;
        if(b.split[i] < a.split[i])
            return false;
    }
    // if we're here, they're equal, so a is NOT less than b.
    return false;
}

inline int edgeset_difference(const PhyEdgeSet &a, const PhyEdgeSet &b)
{
    int sim = 0; 
    for(unsigned int i = 0; i < a.size(); i++)
    {
        for(unsigned int j = 0; j < a.size(); j++)
        {
            if(a[i] == b[j])
            {
                sim++;
                break;
            }
        }
    }

    return a.size() - sim;
}

extern "C"
{
    // returns a distance matrix of doubles
    SEXP phycpp_compute_tree_distance_set(SEXP trees, SEXP verbose)
    {
        SEXP distmatrix;

        bool be_verbose = asLogical(verbose);

        // get number of trees
        int len = length(trees);

        vector<string> treevec(len);

        // convert 'trees' to a vector of tree strings
        for(int i = 0; i < len; i++)
            treevec[i] = CHAR(STRING_ELT(VECTOR_ELT(trees,i),0));

        PROTECT(distmatrix = allocMatrix(REALSXP,len, len));

        compute_phylo_distance_matrix(treevec, be_verbose, REAL(distmatrix));

        // replace -1 values with R N/A values
        for(int i = 0; i < len*len; i++)
            if(REAL(distmatrix)[i] == -1)
                REAL(distmatrix)[i] = R_NaReal;

        UNPROTECT(1);

        return distmatrix;
    }

    SEXP multiset_diff_integer(SEXP set1, SEXP set2)
    {
        unsigned len1 = length(set1);
        int *s1 = INTEGER(set1);
        unsigned len2 = length(set2);
        int *s2 = INTEGER(set2);

        SEXP newset;
        PROTECT(newset = allocVector(INTSXP, len1));
        int *a = INTEGER(newset);

        unsigned p=0;
        for(unsigned i = 0; i < len1; i++)
        {
            bool keep = true;
            for(unsigned j = 0; j < len2; j++)
            {
                if(s1[i] == s2[j])
                {
                    keep = false;
                    break;
                }
            }

            if(keep)
                a[p++] = s1[i];
        }

        for(unsigned i = p; i < len1; i++)
        {
            a[i] = NA_INTEGER;
        }

        
        UNPROTECT(1);
        return newset;
    }

    // returns edge differences
    SEXP phycpp_bin_trees(SEXP treelist)
    {
        // get number of trees
        int len = length(treelist);

        vector<string> treevec(len);

        // convert 'trees' to a vector of tree strings
        for(int i = 0; i < len; i++)
            treevec[i] = CHAR(STRING_ELT(VECTOR_ELT(treelist,i),0));

        std::vector<PhyEdgeSet> trees;
        build_tree_list(treevec, trees, false);

        SEXP distmatrix;
        PROTECT(distmatrix = allocMatrix(REALSXP,len,len));
        double *d = REAL(distmatrix);

        unsigned int sz = trees.size();

        // zero out diagonal
        for(unsigned int i = 0; i < trees.size(); i++)
            d[sz*i + i] = 0;

        for(unsigned int i = 0; i < trees.size(); i++)
        {
            for(unsigned int j = i; j < trees.size(); j++)
            {
                int sim = edgeset_difference(trees[i], trees[j]);
                d[sz*i + j] = sim;
                d[sz*j + i] = sim;
            }
        }

        UNPROTECT(1);

        return distmatrix;
    }

    SEXP gromov_distmatrix(SEXP distmatrix, SEXP bDeltas, SEXP scale_method)
    {
        bool list_deltas = asLogical(bDeltas);
        int scaleM = asInteger(scale_method);

        unsigned len = length(distmatrix);
        unsigned n = sqrt(static_cast<double>(len));
        double *d = REAL(distmatrix);

        SEXP g;

        if(list_deltas)
        {
            PROTECT(g = allocVector(REALSXP, (n*(n-1)*(n-2)*(n-3))/(4*3*2)));
            gromov_graycode(d, n, REAL(g), scaleM);
            UNPROTECT(1);
        }
        else
        {
            PROTECT(g = allocVector(REALSXP, 1));
            REAL(g)[0] = gromov_graycode(d, n, NULL, scaleM); 
            UNPROTECT(1);
        }

        return g;
    }
}

void build_tree_list(std::vector<string> &trees_in, std::vector<PhyEdgeSet> &trees, bool verbose)
{
    string t;

    // build the leaf label LUT 
    t = trees_in[0];
    map<string,unsigned int> strtbl = AssignLeafLabels(t);

    for(unsigned int treeno = 0; treeno < trees_in.size(); treeno++)
    {
        t = trees_in[treeno];

        if(verbose)
            Rprintf("Parsing tree %d\n", treeno);

        PhyEdgeSet tr = NewickParse(t, strtbl);

        ClampNegativeWeights(tr);
        trees.push_back(tr);
    }
}

int compute_phylo_distance_matrix(vector<string> trees_in, bool verbose, double *distance_matrix)
{
    std::vector<PhyEdgeSet> trees;
    build_tree_list(trees_in, trees, verbose);

    int ctr = 0;
    int tot = 0.5 * (trees.size() * (trees.size() - 1));

    // figure out how big the trees are
    int k = -1;
    while(trees[++k].size() == 0) ;

    stl_bool *incompatibility_buffer = new stl_bool[trees[k].size()*trees[k].size()];

    for(unsigned int j = 0; j < trees.size(); j++) // cols
    {
        for(unsigned int i = 0; i < j; i++) // rows
        {
            ctr++;
            if(verbose)
                Rprintf("%d/%d\t\t[%3.2f%%]\n", ctr,tot, (ctr/(double)tot)*100.0);

            if(trees[i].size() == 0 || trees[j].size() == 0) // mark invalid trees with -1 so we can mark them as N/A later
            {
                distance_matrix[i*trees.size()+j] = -1;
                distance_matrix[j*trees.size()+i] = -1;
            }
            else
            {
                double d = TreeDistance(trees[i],trees[j], incompatibility_buffer);
                distance_matrix[i*trees.size()+j] = d;
                distance_matrix[j*trees.size()+i] = d;
            }

        }
    }

    delete [] incompatibility_buffer;

    // fill in the diagonal
    for(unsigned int i = 0; i < trees.size(); i++)
    {
        distance_matrix[i*trees.size()+i] = 0;
    }

    return 0;
}

double gromov_graycode(double *m, size_t n, double *deltas, int scale) 
{
    /* implements Knuth 7.2.1.3 Alg R (Revolving Door) */

    unsigned c[] = {(unsigned)-1,0,1,2,3,(unsigned)n}; 
    unsigned t = 4; 
    unsigned j; 

    double s[3];
    double max = 0;

    double a,b;
    double raw_delta;

    unsigned i = 0;

R2:
    s[0] = m[c[1]*n + c[2]] + m[c[3]*n + c[4]];
    s[1] = m[c[1]*n + c[3]] + m[c[2]*n + c[4]];
    s[2] = m[c[1]*n + c[4]] + m[c[2]*n + c[3]];
    
    /* get the two largest */
    if(s[0] >= s[1])
    {
        a = s[0];
        if(s[1] >= s[2])
            b = s[1];
        else
            b = s[1];
    }
    else
    {
        a = s[1];
        if(s[2] >= s[0])
            b = s[2];
        else
            b = s[1];
    }

    raw_delta = fabs(a-b);
    switch(scale)
    {
        case 2:
            raw_delta /= fmax2(a,b);
            break;
        case 3:
            double dd[4];
            dd[0] = m[c[1]*n + c[4]] + m[c[1]*n + c[3]] + m[c[3]*n + c[4]];
            dd[1] = m[c[1]*n + c[4]] + m[c[1]*n + c[2]] + m[c[2]*n + c[4]];
            dd[2] = m[c[2]*n + c[3]] + m[c[3]*n + c[4]] + m[c[2]*n + c[4]];
            dd[3] = m[c[1]*n + c[2]] + m[c[1]*n + c[3]] + m[c[2]*n + c[3]];

            if(dd[0] >= dd[1] && dd[0] >= dd[2] && dd[0] >= dd[3])
                raw_delta /= dd[0];
            else if(dd[1] >= dd[0] && dd[1] >= dd[2] && dd[1] >= dd[3])
                raw_delta /= dd[1];
            else if(dd[2] >= dd[0] && dd[2] >= dd[1] && dd[2] >= dd[3])
                raw_delta /= dd[2];
            else
                raw_delta /= dd[3];
            break;
    }

    if(deltas)
    {
        deltas[i] = raw_delta;
    }

    if(raw_delta > max)
    {
        max = raw_delta;
    }

    i++;
//R3: /* N.B. assumes t is even */
    if(c[1] > 0) 
    {
        c[1]--;
        goto R2;
    }
    else
    {
        j = 2;
        goto R5;
    }
R4: 
    if(c[j] >= j)
    {
        c[j] = c[j-1];
        c[j-1] = j - 2;
        goto R2;
    }
    else
    {
        j++;
    }
R5: 
    if(c[j] + 1 < c[j+1])
    {
        c[j-1] = c[j];
        c[j]++;
        goto R2;
    }
    else
    {
        j++;
        if(j <= t)
        {
            goto R4; 
        }
    }

    return max/2.0;
}

