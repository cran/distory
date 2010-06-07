/* 
 * Polynomial Time Distance Algorithm for Binary Rooted Phylogenetic Trees
 *
 * Note: Algorithm based on presentation by Scott Provan / Megan Owen,
 *       "Computing Geodesic Distances in Tree Space in Polynomial Time"
 * 
 * John Chakerian
 * chakj@stanford.edu
 * 
 * Last edited: Apr 22 2010
 *
 */

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cctype>
#include <set>
#include <string>
#include <vector>
#include <stack>

#include <R.h>

#include "newick.h"

// Optionally parse a floating point value after a ':' to indicate edge weight.
// Sets the weight to 1 if no ':' is found.

double ParseWeight(const string &str, unsigned int pos, unsigned int *off)
{
    double w = 1;
    if(str[pos] == ':')
    {
        char* e;

        string s = str.substr(pos+1);
        const char* cs = s.c_str();

        w = strtod(cs, &e);
        if(cs == e)
        {
            w = 1;
        }
        else
        {
            pos += (e - cs)+1;
        }
    }

    if(off)
        *off = pos;

    return w;
}

map<string,unsigned int> AssignLeafLabels(const string &str)
{
    unsigned int ctr = 0;
    map<string, unsigned int> ret;

    bool rec = false;
    string t = "";
    for(unsigned int i = 0; i < str.size(); i++)
    {
        if(str[i] == ' ') continue;

        if(str[i] == '(' || str[i] == ',')
        {
            rec = true;
            continue;
        }

        if(rec)
        {
            if(isalpha(str[i]) || isdigit(str[i]) || str[i] == '_' || str[i] == '-')
                t += str[i];
            else
            {
                ret[t] = ctr;
                t = "";
                ctr++;
                rec = false;
            }
        }
    }
    
    return ret;
}

// We'll use a simple state machine to parse over the trees. Technically this
// isn't the best solution - a PDA would be more appropriate, however we're not
// interested in a full parse, only the partitions caused by each edge. Because
// of this, we can avoid recursion or any sort of tree data structure and just
// go right to the partitions. 

PhyEdgeSet NewickParse(const string &str, map<string,unsigned int> &name_to_pos)
{
    unsigned int pos = 0; // position in string
    unsigned int ctr = 0; // nesting level
    unsigned int id_ctr = 0; // for assigning each a unique id
    string errstr;

    vector< set<string> > in_leaves(30);
    stack<int> commas; // for ensuring trees are bifurcating

    vector<PhyEdge> edges;

    if(str.length() == 0)
    {
        return vector<PhyEdge>();
    }

    // move past leading whitespace
    while(str[pos] == ' ' || str[pos] == '\t')
        pos++;

    // we can bail out early here if we have an empty string or a ; character
    if(str[pos] == ';')
    {
        return vector<PhyEdge>();
    }

    // if the last non-whitespace character isn't a ;, bail out
    size_t tpos = str.length()-1;
    while(tpos > 0 && isspace(str[tpos])) 
        tpos--;

    if(tpos == 0) // empty string
    {
        return vector<PhyEdge>();
    }

    if(tpos > 0 && str[tpos] != ';')
    {
        errstr = "Tree not terminated with ';'";
        goto error;
    }

    goto start_state;

start_state:
    if(str[pos] == ';')
        goto finish;
    else if(str[pos] == '(')
        goto new_nesting;
    else
    {
        errstr = "Parse error at first character. A tree must have at least 2 taxa.";
        goto error;
    }

new_nesting:
    ctr++;
    in_leaves.push_back(set<string>());
    commas.push(0);

    pos++; // move to first character of new nesting

    // move past whitespace
    while(str[pos] == ' ' || str[pos] == '\t')
        pos++;

    if(str[pos] == '(')
        goto new_nesting; // start the process over again
    else if(isalpha(str[pos]) || isdigit(str[pos]) || str[pos] == '_' || str[pos] == '-')
        goto leaf_entry;
    else if(str[pos] == ':')
    {
        errstr = "Leaf nodes must have string labels.";
        goto error;
    }
    else if(str[pos] == ')')
    {
        errstr = "Empty blocks of form () are not allowed.";
        goto error;
    }
    else
    {
        errstr = "Parse error after '('";
        goto error;
    }

leaf_entry:
    // move past whitespace
    while(str[pos] == ' ' || str[pos] == '\t')
        pos++;

    // add the edge to the leaf
    {
        // get the leaf name
        unsigned int spos = pos;
        while(isalpha(str[pos]) || isdigit(str[pos]) || str[pos] == '_' || str[pos] == '-')
        {
            pos++;
        }

        string leaflabel = str.substr(spos,pos-spos);

        PhyEdge e(name_to_pos.size());
        e.id = id_ctr++;
        double w = ParseWeight(str,pos,&pos);
        e.weight = w;
        e.split[ name_to_pos[leaflabel] ] = false;

        edges.push_back(e);

        for(unsigned int i = 0; i < in_leaves.size(); i++)
        {
            in_leaves[i].insert(leaflabel);
        }
    }

    // move past extra whitespace
    while(str[pos] == ' ' || str[pos] == '\t')
        pos++;

    // decide what to do next
    if(str[pos] == ',')
    {
        commas.top()++;
        if(commas.top() > 1)
        {
            errstr = "More than one comma detected in a nesting.";
            goto error;
        }
        pos++;
        // move past whitespace
        while(str[pos] == ' ' || str[pos] == '\t')
            pos++;

        if(str[pos] == '(')
            goto new_nesting;
        else if(isdigit(str[pos]) || isalpha(str[pos]) || str[pos] == '_' || str[pos] == '-')
            goto leaf_entry;
        else
        {
            errstr = "Couldn't figure out what to do when moving past a ','.";
            goto error;
        }
    }
    else if(str[pos] == ')')
        goto close_nesting;
    else
    {
        errstr = "Couldn't figure out what to do after a first leaf entry ";
        goto error;
    }

close_nesting:
    if(str[pos] != ')')
    {
        errstr = "Close-nesting expected but no end-bracket found.";
        goto error;
    }
    pos++;
    if(str[pos] == ';')
        goto finish;

    // add the edge
    {
        PhyEdge e(name_to_pos.size());
        e.id = id_ctr++;
        e.weight = ParseWeight(str,pos,&pos);
        
        for(set<string>::iterator it = in_leaves[in_leaves.size()-1].begin(); it != in_leaves[in_leaves.size()-1].end(); ++it)
            e.split[ name_to_pos[*it] ] = false;

        edges.push_back(e);
    }
    in_leaves.pop_back();
    commas.pop();

    if(str[pos] == ';')
        goto finish;

    // decide what to do next
    if(str[pos] == ',')
    {
        commas.top()++;
        if(commas.top() > 1)
        {
            errstr = "More than one comma detected in a nesting. Trees MUST be binary.";
            goto error;
        }
        pos++;
        // move past whitespace
        while(str[pos] == ' ' || str[pos] == '\t')
            pos++;

        if(str[pos] == '(')
            goto new_nesting;
        else if(isdigit(str[pos]) || isalpha(str[pos]) || str[pos] == '_' || str[pos] == '-')
            goto leaf_entry;
        else
        {
            errstr = "Couldn't figure out what to do when moving past a ','.";
            goto error;
        }
    }
    else if(str[pos] == ')')
        goto close_nesting;
    else
    {
        errstr = "Couldn't figure out what to do after a second leaf entry";
        goto error;
    }

error:
    if(errstr == "") errstr = "Parser ran off the edge.";
    error("An error was encountered in parsing near position %d: %s\n", pos, errstr.c_str());

    return vector<PhyEdge>();
    
finish:
    return edges;
}

