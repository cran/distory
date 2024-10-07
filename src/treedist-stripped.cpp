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
 * Note that this is a stripped version because it includes only a built-in
 * Edmonds-Karp network flow implementation; linear programming and Boost.Graph
 * implementations are not provided, since they would introduce library 
 * dependencies. 
 *
 */

#define R_NO_REMAP

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>

#include "treedist.h"

#include "phylo.h"
#include "newick.h"

using namespace std;

// Callback for clamping 0 or negative weights at a Very Small Value (that can
// still be squared)
void ClampNegativeWeights(PhyEdgeSet &a)
{
    for(unsigned int i = 0; i < a.size(); i++)
        if(a[i].weight < sqrt(DBL_MIN))
            a[i].weight = sqrt(DBL_MIN);
}

// Two edges are compatible if there is an intersection between the splits
// that is empty.
bool EdgesCompatible(PhyEdge &e1, PhyEdge &e2)
{
    // n.b. - left == false, right == true
    bool n1 = true;
    bool n2 = true;
    bool n3 = true;
    bool n4 = true;

    for(unsigned int i = 0; i < e1.split.size(); i++)
    {
        // left/left
        if(e1.split[i] == false && e2.split[i] == false)
            n1 = false;
        // right/right
        if(e1.split[i] == true && e2.split[i] == true)
            n2 = false;
        // right/left
        if(e1.split[i] == true && e2.split[i] == false)
            n3 = false;
        // left/right
        if(e1.split[i] == false && e2.split[i] == true)
            n4 = false;
    }

    return (n1 || n2 || n3 || n4);
}

double TreeDistance(PhyEdgeSet &a, PhyEdgeSet &b, stl_bool* incompatible)
{
    // 1. Separate out into shared edges & unique edges:
    //     a) foreach edge e in a:
    //          i) search b for any edge e' with a leaf set equivalent to e
    //          ii) if found, 'remove' e and e' and store as a pair in
    //              shared_edges
    //     b) all remaining edges are unique to each
    
    assert(a.size() == b.size());
 
    std::vector<pair<PhyEdge,PhyEdge> > shared_edges; 
    shared_edges.reserve(a.size());

    vector<int> a_shared_idxs;
    a_shared_idxs.reserve(a.size());

    vector<int> b_shared_idxs;
    b_shared_idxs.reserve(b.size());

    for(unsigned int i = 0; i < a.size(); i++)
    {
        for(unsigned int j = 0; j < b.size(); j++)
        {
            if(a[i] == b[j])
            {
                shared_edges.push_back( make_pair(a[i], b[j]) );
                a_shared_idxs.push_back(i);
                b_shared_idxs.push_back(j);
                break; // stop checking for this i
            }
        }
    }

    // The not-shared edges
    PhyEdgeSet p,q;

    for(unsigned int i = 0; i < a.size(); i++)
        if(find(a_shared_idxs.begin(), a_shared_idxs.end(), i) == a_shared_idxs.end())
            p.push_back(a[i]);

    for(unsigned int j = 0; j < b.size(); j++)
        if(find(b_shared_idxs.begin(), b_shared_idxs.end(), j) == b_shared_idxs.end())
            q.push_back(b[j]);

    std::vector<pair<PhyEdgeSet,PhyEdgeSet> > bins(shared_edges.size()+1); // last one is generic bin

    // 1.5. For every unique edge, classify under the tightest shared edge,

    // Classify left-tree edges
    for(unsigned int i = 0; i < p.size(); i++)
    {
        unsigned int smallest_subset = UINT_MAX;
        unsigned int edge_id = UINT_MAX;

        for(unsigned int j = 0; j < shared_edges.size(); j++)
        {
            unsigned int r = shared_edges[j].first.SubsetRemainder(p[i]);
            if(r < smallest_subset)
            {
                smallest_subset = r;
                edge_id = j;
                if(smallest_subset == 1) break;
            }
        }

        if(edge_id == UINT_MAX) // dump it in the generic bin
            bins[shared_edges.size()].first.push_back(p[i]);
        else
            bins[edge_id].first.push_back(p[i]);
    }

    // Classify right-tree edges
    for(unsigned int i = 0; i < q.size(); i++)
    {
        unsigned int smallest_subset = UINT_MAX;
        unsigned int edge_id = UINT_MAX;

        for(unsigned int j = 0; j < shared_edges.size(); j++)
        {
            unsigned int r = shared_edges[j].second.SubsetRemainder(q[i]);
            if(r < smallest_subset)
            {
                smallest_subset = r;
                edge_id = j;
                if(smallest_subset == 1) break;
            }
        }

        if(edge_id == UINT_MAX) // dump it in the generic bin
            bins[shared_edges.size()].second.push_back(q[i]);
        else
            bins[edge_id].second.push_back(q[i]);
    }

    // Precache all incompatibilities in a 'huge' logical array
    // (this could be made into a bitmatrix later, if space is in demand)
    // (this could also be upper-triangular, but again, space isn't a big deal
    // at the moment)
   
    // Keep in mind here that the IDs are _only_ valid within a certain tree.
    // Thus the ordering always must be (left-hand tree edge , right-hand tree
    // edge), and the matrix is NOT symmetric

    unsigned int max_id = a.size(); 

    for(unsigned int binid = 0; binid < bins.size(); binid++)
        for(unsigned int i = 0; i < bins[binid].first.size(); i++)
            for(unsigned int j = 0; j < bins[binid].second.size(); j++)
                incompatible[bins[binid].first[i].id*max_id + bins[binid].second[j].id] = EdgesCompatible(bins[binid].first[i], bins[binid].second[j]);

    // 2. For every shared edge bin, compute distance for all unique edges under
    // it.
    
    double unique_edge_factor = 0;
    for(unsigned int i = 0; i < bins.size(); i++)
        if(bins[i].first.size() > 0)
            unique_edge_factor += DisjointTreeDistance(bins[i].first, bins[i].second, incompatible, max_id);

    // 3. Calculate the final distance from both shared & unique edges
    // a) the shared edges have Euclidean distance:

    double shared_edge_factor = 0;
    for(unsigned int i = 0; i < shared_edges.size(); i++)
    {
        double t = shared_edges[i].first.weight - shared_edges[i].second.weight;
        shared_edge_factor += t*t;
    }

    // c) combine them together and return
    return sqrt(unique_edge_factor + shared_edge_factor);
}

struct NetworkFlowResult
{
    double flow;
    PhyEdgeSet A_i, B_i, A_i_prime, B_i_prime;
};

NetworkFlowResult EKNetworkFlow(PhyEdgeSet &A_i, PhyEdgeSet& B_i, stl_bool* incompatible, unsigned int max_id) 
{
    // make these convenient, since they won't be changing
    unsigned int ac = A_i.size(); 
    unsigned int bc = B_i.size(); 

    // compute square sums of A_i and B_i
    double sumsq_A_i = 0;
    for(unsigned int i = 0; i < ac; i++)
        sumsq_A_i += A_i[i].weight*A_i[i].weight;

    double sumsq_B_i = 0;
    for(unsigned int j = 0; j < bc; j++)
        sumsq_B_i += B_i[j].weight*B_i[j].weight;

    // A_i first, then B_i, then source, then sink
    unsigned int N = ac + bc + 2; // for offsets
    unsigned int s = ac+bc;     // source index
    unsigned int t = ac+bc+1;   // sink index

    // set up the variables for Edmonds-Karp

    // total flow
    double f = 0;

    double *F = new double[N*N]; // flow
    double *C = new double[N*N]; // capacity

    // set all the flow to zero
    fill(F,F + N*N,0);
    fill(C,C + N*N,0); // necessary?

    // add capacities from source to A_i
    for(unsigned int i = 0; i < A_i.size(); i++) // for (s,i)
    {
        C[N*s + i] = (A_i[i].weight*A_i[i].weight)/sumsq_A_i;
        C[N*i + s] = 0;
    }

    // add capacities from B_i to sink
    for(unsigned int j = 0; j < B_i.size(); j++) // for (j, t)
    {
        C[N*(ac+j) + t] = (B_i[j].weight*B_i[j].weight)/sumsq_B_i;
        C[N*t + (ac+j)] = 0;
    }

    // add in incompatible edges
    vector<pair<unsigned int,unsigned int> > incompatible_edge_pairs;
    for(unsigned int i = 0; i < A_i.size(); i++)
    {
        for(unsigned int j = 0; j < B_i.size(); j++)
        {
            if(incompatible[A_i[i].id*max_id + B_i[j].id ]) // add (i,j) and (j,i)
            {
                C[N*i+(ac+j)] = DBL_MAX;
                C[N*(ac+j)+i] = 0;
            }
        }
    }

    int *P = new int[N]; // path (reused)
    double *M = new double[N]; // capacity to node

    while(true)
    {
        // clean up old values, init the path
        fill(P, P+N, -1);
        fill(M, M+N, DBL_MAX);
        P[s] = -2;

        // do BFS to get a new path

        std::queue<unsigned int> tc; // to-check
        tc.push(s);

        bool res = false; // make sure we get a result
        while(!tc.empty())
        {
            unsigned int u = tc.front();
            tc.pop();

            // we need to get OUT-edges
            // decide if we're on a A_i, B_i, s, or t
            vector<unsigned int> E;
            if(u == s)
            {
                // add the nodes in A_i
                E.resize(ac,0);
                for(unsigned int i = 0; i < ac; i++)
                    E[i] = i;
            }
            else if(u == t)
            {
                // add the nodes in B_i
                E.resize(bc,0);
                for(unsigned int i = 0; i < bc; i++)
                    E[i] = ac + i;
            }
            else if(u < ac) // then it is in A_i
            {
                // add the nodes in A_i
                for(unsigned int j = 0; j < bc; j++)
                    if(incompatible[A_i[u].id*max_id + B_i[j].id ])
                        E.push_back(ac + j);
                // add the source - don't think we actually need to do this
                E.push_back(s);
            }
            else // then it is in B_i
            {
                // add the nodes in A_i
                for(unsigned int j = 0; j < ac; j++)
                    if(incompatible[A_i[j].id*max_id + B_i[u-ac].id ])
                        E.push_back(j);
                // add the sink
                E.push_back(t);
            }

            // now iterate through E 
            for(unsigned int i = 0; i < E.size(); i++)
            {
                if(C[ N*u+E[i] ] - F[ N*u+E[i] ] > 0 && P[E[i]] == -1)
                {
                    P[E[i]] = u;
                    M[E[i]] = min(M[u],C[ N*u+E[i] ] - F[ N*u+E[i] ]);
                    if(E[i] != t)
                        tc.push(E[i]);           
                    else
                    {
                        // we're done with the BFS - empty the queue out
                        while(!tc.empty())
                            tc.pop();
                        res = true;
                    }
                }
            }
        }
        
        double m = 0;
        if(res) m = M[t];

        // check capacity of the path, break if 0 (no new path found)
        if(m == 0)
            break; 

        // add capacity of the path to flow
        f += m;

        // travel from sink to source along the path, updating edges with +/- cap
        unsigned int v = t;
        while(v != s)
        {
            unsigned int u = P[v];
            F[N*u + v] += m;
            F[N*v + u] -= m;
            v = u;
        }
    }

    NetworkFlowResult r;
    r.flow = f;

    if(f < 1 - DBL_EPSILON*100)
    {
        // do a BFS to get splits

        // reset P
        fill(P, P+N, -1);

        std::queue<unsigned int> Q;
        Q.push(s);

        while(!Q.empty())
        {
            unsigned int u = Q.front();
            Q.pop();

            // we need to get OUT-edges
            // decide if we're on a A_i, B_i, s, or t
            vector<unsigned int> E;
            if(u == s)
            {
                // add the nodes in A_i
                E.resize(ac,0);
                for(unsigned int i = 0; i < ac; i++)
                    E[i] = i;
            }
            else if(u == t)
            {
                // add the nodes in B_i
                E.resize(bc,0);
                for(unsigned int i = 0; i < bc; i++)
                    E[i] = ac + i;
            }
            else if(u < ac) // then it is in A_i
            {
                // add the nodes in A_i
                for(unsigned int j = 0; j < bc; j++)
                    if(incompatible[ A_i[u].id*max_id + B_i[j].id ])
                        E.push_back(ac + j);
                // add the source - don't think we actually need to do this
                E.push_back(s);
            }
            else // then it is in B_i
            {
                // add the nodes in A_i
                for(unsigned int j = 0; j < ac; j++)
                    if(incompatible[A_i[j].id*max_id + B_i[u-ac].id ])
                        E.push_back(j);
                // add the sink
                E.push_back(t);
            }

            // now iterate through E 
            for(unsigned int i = 0; i < E.size(); i++)
            {
                if(C[ N*u+E[i] ] - F[ N*u+E[i] ] > 0 && P[E[i]] == -1)
                {
                    P[E[i]] = 1;
                    Q.push(E[i]); 
                }
            }
        }

        // form our split vectors
        for(unsigned int i = 0; i < N-2; i++) // ignore the last 2, since source/sink
        {
            if(i < ac) // in A_i
            {
                if(P[i] == 1)
                    r.A_i_prime.push_back(A_i[i]);
                else
                    r.A_i.push_back(A_i[i]);
            }
            else // in B_i
            {
                if(P[i] == 1)
                    r.B_i_prime.push_back(B_i[i-ac]);
                else
                    r.B_i.push_back(B_i[i-ac]);
            }
        }
    }

    delete [] C;
    delete [] F;
    delete [] P;
    delete [] M;

    return r;
}

double DisjointTreeDistance(PhyEdgeSet &a, PhyEdgeSet &b, stl_bool *incompatible, unsigned int max_id)
{
    // 2. For unique edges:
    //      a) initialize A_1 to all unique edges in a, B_1 to all unique edges in b
    //      b) for each edge e in A_i find all incompatible edges in b
    //      c) construct a bipartite graph based on these incompatible edges
    //      d) run max flow
    //      e) do BFS to find the min-weight vertex cover
    //      f) if min-weight vertex cover has weight < 1, split based on
    //          min-weight vertex cover and push blocks into the right places, else move to next A_i

    stack<pair<PhyEdgeSet, PhyEdgeSet> > unchecked_blocks;
    std::vector<pair<PhyEdgeSet, PhyEdgeSet> > finished_blocks; // NB blocks in reverse order

    pair<PhyEdgeSet, PhyEdgeSet> p;
    for(unsigned int i = 0; i < a.size(); i++)
        p.first.push_back(a[i]);

    for(unsigned int j = 0; j < b.size(); j++)
        p.second.push_back(b[j]);

    unchecked_blocks.push(p);

    while(!unchecked_blocks.empty())
    {
        PhyEdgeSet A_i = unchecked_blocks.top().first;
        PhyEdgeSet B_i = unchecked_blocks.top().second;
        unchecked_blocks.pop();

        // Some sanity checks
        if(A_i.size() == 0 || B_i.size() == 0)
        {
            finished_blocks.push_back( make_pair(A_i, B_i) ); 
            continue;
        }

        if(A_i.size() == 1 && B_i.size() == 1)
        {
            finished_blocks.push_back( make_pair(A_i, B_i) ); 
            continue;
        }

        NetworkFlowResult r = EKNetworkFlow(A_i, B_i, incompatible, max_id);

        if(r.flow < 1 - DBL_EPSILON*100) 
        {
            unchecked_blocks.push( make_pair(r.A_i_prime, r.B_i_prime) );
            unchecked_blocks.push( make_pair(r.A_i, r.B_i) );
        }
        else
        {
            finished_blocks.push_back( make_pair(A_i, B_i) );
        }
    }

    // the unique edges have distance sqrt( sum_i (||A_i|| + ||B_i||)^2 )
    double unique_edge_factor = 0;
    for(unsigned int i = 0; i < finished_blocks.size(); i++)
    {
        double accum = 0;

        PhyEdgeSet a = finished_blocks[i].first;
        PhyEdgeSet b = finished_blocks[i].second;

        // do the A_i block
        double t = 0;
        for(unsigned j = 0; j < a.size(); j++)
            t += a[j].weight*a[j].weight;
        t = sqrt(t);

        accum += t;

        // do the B_i block
        t = 0;
        for(unsigned j = 0; j < b.size(); j++)
            t += b[j].weight*b[j].weight;
        t = sqrt(t);

        accum += t;

        unique_edge_factor += accum*accum;
    }

    return unique_edge_factor; // n.b., not square rooted
}

