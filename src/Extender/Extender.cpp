#include "Extender.h"

void Extension(Kmer kmer, DBG &dbg, vector<Kmer> & unitig)
{
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    if (dbg.in_degree(kmer) == 1 && neighbors.size() == 1) {
        unitig.push_back(kmer);
        Extension(neighbors[0],dbg,unitig);
    }
}

vector<Kmer> UnitigExtender::Extend(Kmer kmer, DBG &dbg)
{
    vector<Kmer> unitig;
    unitig.push_back(kmer);
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    if (neighbors.size()==1)
        Extension(neighbors[0], dbg, unitig);
    return unitig;
}