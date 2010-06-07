#ifndef _PHYLO_H_
#define _PHYLO_H_

#include <string>
#include <set>
#include <map>
#include <vector>
#include <climits>

using namespace std;

typedef unsigned char stl_bool; // to prevent bit packing, which slows access time down considerably due to non-addressability of bits

struct PhyEdge
{
    double weight;

    unsigned int id; // used so we can identify edges without doing the set comparisons

    vector<stl_bool> split; // downward == false, upward == true

    PhyEdge(unsigned short leafcount) 
    { 
        /* this is a hack to add in the root node implicitly. it is added as an
         * upward edge that is never assigned to. */
        split.resize(leafcount+1, true);
    }

    bool operator==(const PhyEdge &rhs) const  // for inter-tree comparisons
    {
        return (split == rhs.split); /* two edges are the same if their splits of leaves are the same */
    }

    bool operator<(const PhyEdge &rhs) const // FOR ORDERING WITHIN A SINGLE TREE ONLY
    {
        return id < rhs.id;
    }

    inline unsigned int SubsetRemainder(const PhyEdge &rhs) const
    {
        bool subset = true;

        // check left subset
        for(unsigned int i = 0; i < rhs.split.size(); i++)
        {
            if(rhs.split[i] == false && split[i] != false)
            {
                subset = false;
                break;
            }
        }

        if(subset)
        {
            // count # of lefts in this 
            int nl = 0;
            for(unsigned int i = 0; i < split.size(); i++)
                if(split[i] == false) nl++;
            // count # of lefts in RHS
            int nr = 0;
            for(unsigned int i = 0; i < rhs.split.size(); i++)
                if(rhs.split[i] == false) nr++;

            return nl - nr;
        }

        // check right subset

        subset = true; // reset flag
        for(unsigned int i = 0; i < rhs.split.size(); i++)
        {
            if(rhs.split[i] == true && split[i] != false)
            {
                subset = false;
                break;
            }
        }

        if(subset)
        {
            // count # of lefts in this 
            int nl = 0;
            for(unsigned int i = 0; i < split.size(); i++)
                if(split[i] == false) nl++;
            // count # of lefts in RHS
            int nr = 0;
            for(unsigned int i = 0; i < rhs.split.size(); i++)
                if(rhs.split[i] == true) nr++;

            return nl - nr;
        }

        return UINT_MAX;
    }
};

typedef vector<PhyEdge> PhyEdgeSet;

#endif /* _PHYLO_H_ */

