#include "DBG.h"
bool NaiveDBG::is_solid(Kmer kmer) const
{
    std::unordered_set<Kmer>::const_iterator place
            = _dbg_naive.find(kmer);
    return (place != _dbg_naive.end());
}

vector<DnaSequence::NuclType> NaiveDBG::getNeighbors
        (const Kmer& kmer) const
{
    vector<DnaSequence::NuclType> nts;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        Kmer kmer_aux = kmer;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back( i);
    }
    return nts;
}

vector<Kmer> NaiveDBG::getKmerNeighbors
        (const Kmer & kmer) const
{
    vector<Kmer> nts;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        Kmer kmer_aux = kmer;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back( kmer_aux);
    }
    return nts;
}

size_t NaiveDBG::in_degree(Kmer k)
{
    size_t out = 0;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        Kmer kmer_aux = k;
        kmer_aux.appendLeft(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}

size_t NaiveDBG::out_degree(Kmer k)
{
    size_t out = 0;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i){
        kmer_aux = k;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
