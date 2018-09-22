#include "DBG.h"
/*
 * First false template -> single_end reads
 */
template<>
bool NaiveDBG<false>::is_solid(typename NodeType::DBGNode kmer) const
{
    return (_dbg_naive.find(kmer) != _dbg_naive.end());
}

template<>
vector<DnaSequence::NuclType> NaiveDBG<false>::getNeighbors
        (const typename NodeType::DBGNode& kmer) const
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

template<>
vector<typename NodeType<false>::DBGNode> NaiveDBG<false>::getKmerNeighbors
        (const typename NodeType<false>::DBGNode & kmer) const
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

template<>
size_t NaiveDBG<false>::in_degree(typename NodeType<false>::DBGNode k)
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

template<>
size_t NaiveDBG<false>::out_degree(typename NodeType<false>::DBGNode k)
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

template<>
void NaiveDBG<false>::show_info()
{
    for (auto p_k: _kmers_map)
    {
        std::cout << p_k.first.str() <<" -> "<< p_k.second.first<<", Correct? "
                  <<(_dbg_naive.find(p_k.first) == _dbg_naive.end())<<"\n";
    }
}

/*
 * Paired_end reads
 */

template<>
bool NaiveDBG<true>::is_solid(typename NodeType::DBGNode kmer) const
{
    return (_dbg_naive.find(kmer) != _dbg_naive.end());
}

template<>
vector<DnaSequence::NuclType> NaiveDBG<true>::getNeighbors
        (const typename NodeType::DBGNode& kmer) const
{
    vector<DnaSequence::NuclType> nts;
    return nts;
}

template<>
vector<typename NodeType<true>::DBGNode> NaiveDBG<true>::getKmerNeighbors
        (const typename NodeType::DBGNode & kmer) const
{
    vector<Kmer> nts;
    return nts;
}

template<>
size_t NaiveDBG<true>::in_degree(typename NodeType<false>::DBGNode k)
{
    size_t out = 0;
    return out;
}

template<>
size_t NaiveDBG<true>::out_degree(typename NodeType<false>::DBGNode k)
{
    size_t out = 0;
    return out;
}

template<>
void NaiveDBG<true>::show_info()
{

}
