#include "Extender.h"

void UnitigExtender::Extend(Kmer kmer, vector <Kmer> &unitig, bool first)
{
    vector<Kmer> neighbors = _dbg.getKmerNeighbors(kmer);
    if ((_dbg.in_degree(kmer) == 1 || first) && neighbors.size()==1) {
        unitig.push_back(kmer);
        Extend(neighbors[0],unitig, false);
    }else
        unitig.push_back(kmer);
}