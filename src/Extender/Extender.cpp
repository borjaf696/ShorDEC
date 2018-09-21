#include "Extender.h"

void Extension(Kmer kmer, DBG &dbg, vector<Kmer> & unitig, vector<Kmer> & out, vector<Kmer> & in, unordered_set<Kmer> added)
{
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    size_t in_ = dbg.in_degree(kmer);
    if (in_ == 1 && neighbors.size() == 1) {
        unitig.push_back(kmer);
        Extension(neighbors[0],dbg,unitig,out,in,added);
    }else if (added.find(kmer) == added.end())
    {
            added.emplace(kmer);
            if (neighbors.size() > 1) {
                out.push_back(kmer);
            }else
                in.push_back(kmer);
            unitig.push_back(kmer);
    }
}

vector<vector<Kmer>> UnitigExtender::Extend(Kmer kmer, DBG &dbg, vector<Kmer> & out, vector<Kmer> & in,unordered_set<Kmer> added)
{
    vector<vector<Kmer>> unitigs;
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    for (auto &k: neighbors)
    {
        vector<Kmer> unitig;
        unitig.push_back(kmer);
        Extension(k, dbg, unitig, out, in, added);
        unitigs.push_back(unitig);
    }
    return unitigs;
}