#include "DBG.h"
bool NaiveDGB::is_solid(Kmer kmer) const {
    std::unordered_set<Kmer>::const_iterator place
            = _dbg_naive.find(kmer);
    return (place != _dbg_naive.end());
}

std::vector<DnaSequence::NuclType> NaiveDGB::getNeighbors
        (const Kmer& kmer) const{
    std::vector<DnaSequence::NuclType> nts;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        kmer_aux = kmer;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back( i);
    }
    return nts;
}

size_t NaiveDGB::in_degree(Kmer k){
    size_t out = 0;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = k;
        kmer_aux.appendLeft(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}

size_t NaiveDGB::out_degree(Kmer k){
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

void NaiveDGB::check_path(Kmer kmer, size_t &len) const {
    std::vector<DnaSequence::NuclType> neigh = getNeighbors(kmer);
    if (neigh.size() == 1) {
        len++;
        if (len < MIN_PATH_LEN){
            Kmer kmer_aux = kmer;
            kmer_aux.appendRight(neigh[0]);
            check_path(kmer_aux,len);
        }
    }if (neigh.size() > 1)
        len += MIN_PATH_LEN;
    else
        len += 0;
}