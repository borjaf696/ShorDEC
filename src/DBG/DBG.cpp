#include "DBG.h"
/*
 * First false template -> single_end reads
 */
template<>
bool NaiveDBG<false>::is_solid(typename NodeType::DBGNode& kmer) const
{
    if (_is_standard) {
        Kmer kmer_aux = kmer;
        kmer_aux.standard();
        return (_dbg_naive.find(kmer_aux) != _dbg_naive.end());
    }else
        return (_dbg_naive.find(kmer) != _dbg_naive.end());
}

template<>
vector<DnaSequence::NuclType> NaiveDBG<false>::getNeighbors
        (typename NodeType::DBGNode kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<DnaSequence::NuclType> nts;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(i);
    }
    return nts;
}

template<>
vector<typename NodeType<false>::DBGNode> NaiveDBG<false>::getKmerNeighbors
        (typename NodeType<false>::DBGNode kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<Kmer> nts;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(kmer_aux.substr(1,Parameters::get().kmerSize));
    }
    return nts;
}

template<>
size_t NaiveDBG<false>::in_degree(typename NodeType<false>::DBGNode k)
{
    size_t out = 0;
    Kmer kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(k.str());
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
        kmer_aux = Kmer(k.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
template<>
void NaiveDBG<false>::_kmerCount() {
        Kmer kmer;
        KmerInfo<false> tail;
        bool first = false;
        for (auto &read:_sc.getIndex()){
            if (first)
                _tails.emplace(tail);
            Progress::update(read.first.getId());
            first = false;
            for (auto kmer_r: IterKmers<false>(read.second.sequence)) {
                /*
                 * TODO: Para hacer la parte de almancenar solo el kmer de menor tamanho chequear cada kmer con standar.
                 * recordar que falta testear el metodo y que por tanto antes de nada habra que asegurar esto.
                 */
                kmer = kmer_r.kmer;
                /*
                 * Lets change into standard form
                 */
                if (_is_standard)
                    kmer.standard();
                /*std::cout << "Kmer: "<<kmer.str()<<" Kmer(rc): " << kmer.rc().str()<<" ¿Es Menor? "<<(kmer < kmer.rc())<<"\n";
                if (kmer < kmer.rc())
                    exit(1);*/
                unordered_map<Kmer, pair<size_t,size_t>>::const_iterator place =
                                                                                 _kmers_map.find(kmer);
                if (place != _kmers_map.end()) {
                    _kmers_map[kmer].first++;
                    _kmers_map[kmer].second = min(_kmers_map[kmer].second,kmer_r.kmer_pos);
                    if (_kmers_map[kmer].first == Parameters::get().accumulative_h)
                    {
                        if (_kmers_map[kmer].second < Parameters::get().kmerSize / 2) {
                            first = true;
                            _heads.emplace(kmer_r);
                        }
                        tail = kmer_r;
                        _dbg_naive.emplace(kmer);
                        _dbg_nodes.emplace(kmer.substr(0,Parameters::get().kmerSize-1));
                        _dbg_nodes.emplace(kmer.substr(1,Parameters::get().kmerSize));
                    }
                } else
                    _kmers_map[kmer] = pair<size_t,size_t>(1,kmer_r.kmer_pos);
                if (_kmers_map[kmer].first == Parameters::get().accumulative_h) {
                    _dbg_naive.emplace(kmer);
                    _dbg_nodes.emplace(kmer.substr(0,Parameters::get().kmerSize-1));
                    _dbg_nodes.emplace(kmer.substr(1,Parameters::get().kmerSize));
                }
            }
        }
        std::cout << "Sizes: "<<_kmers_map.size()<<" "<<_dbg_naive.size()<<" "<<_dbg_nodes.size()<<"\n";
        _kmers_map.clear();
        /*for (auto &k: _dbg_naive)
            cout << "KmerNaive: "<<k.str()<<"\n";*/
        /*for (auto &k: _dbg_nodes) {
            cout << "KmerNodes: " << k.str() << "\n";
            vector<Kmer> neigh = getKmerNeighbors(k);
            for (auto n: neigh)
                cout << "Vecinos: "<<n.str() << "\n";
        }*/
        Progress::update(_sc.getIndex().size());
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
void NaiveDBG<true>::_kmerCount() {
}

template<>
bool NaiveDBG<true>::is_solid(typename NodeType::DBGNode& kmer) const
{
    return (_dbg_naive.find(kmer) != _dbg_naive.end());
}

template<>
vector<DnaSequence::NuclType> NaiveDBG<true>::getNeighbors
        (typename NodeType::DBGNode kmer) const
{
    vector<DnaSequence::NuclType> nts;
    return nts;
}

template<>
vector<typename NodeType<true>::DBGNode> NaiveDBG<true>::getKmerNeighbors
        (typename NodeType::DBGNode  kmer) const
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
