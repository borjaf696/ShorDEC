#include "DBG.h"
bool NaiveDGB::is_solid(Kmer kmer) const {
    /*std::unordered_set<std::string>::const_iterator place
            = _dbg_naive.find(kmer.str());
    return (place != _dbg_naive.end());*/
    std::unordered_set<Kmer>::const_iterator place
            = _dbg_naive.find(kmer);
    return (place != _dbg_naive.end());
}

std::vector<std::pair<Kmer,DnaSequence::NuclType>> NaiveDGB::getNeighbors
        (const Kmer& kmer) const{
    std::vector<std::pair<Kmer,DnaSequence::NuclType>> nts;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        Kmer kmer_aux = kmer;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(std::pair<Kmer, DnaSequence::NuclType>(kmer_aux, i));
    }
    return nts;
}

size_t NaiveDGB::in_degree(Kmer k){
    size_t out = 0;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        Kmer kmer_aux = k;
        kmer_aux.appendLeft(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}

size_t NaiveDGB::out_degree(Kmer k){
    size_t out = 0;
    for (DnaSequence::NuclType i = 0; i < 4; ++i){
        Kmer kmer_aux = k;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}

void NaiveDGB::check_path(Kmer kmer, size_t &len) const {
    std::vector<std::pair<Kmer,DnaSequence::NuclType>> neigh = getNeighbors(kmer);
    if (neigh.size() == 1) {
        len++;
        check_path(neigh[0].first, len);
    }if (neigh.size() > 1)
        len += MIN_PATH_LEN;
    else
        len += 0;
}