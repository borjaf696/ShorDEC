#include "DBG.h"
bool NaiveDBG::is_solid(Kmer kmer) const
{
    return (_dbg_naive.find(kmer) != _dbg_naive.end());
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
    //std::cout << "Original Kmer: "<<kmer.str()<<"\n";
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        Kmer kmer_aux = kmer;
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back( kmer_aux);
    }
    /*for (uint i = 0; i < nts.size(); ++i)
        std::cout <<"Neighbor: " <<nts[i].str() << "\n";*/
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

void NaiveDBG::show_info()
{
    for (auto p_k: _kmers_map)
    {
        std::cout << p_k.first.str() <<" -> "<< p_k.second.first<<", Correct? "
                  <<(_dbg_naive.find(p_k.first) == _dbg_naive.end())<<"\n";
    }
}
