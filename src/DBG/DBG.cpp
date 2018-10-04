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

template <>
void NaiveDBG<false>::_insert(Kmer kmer, Kmer rc)
{
    _dbg_naive.emplace(kmer);
    _dbg_nodes.emplace(kmer.substr(0,Parameters::get().kmerSize-1));
    _dbg_nodes.emplace(kmer.substr(1,Parameters::get().kmerSize));
    _dbg_nodes.emplace(rc.substr(0,Parameters::get().kmerSize-1));
    _dbg_nodes.emplace(rc.substr(1,Parameters::get().kmerSize));
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
                kmer = kmer_r.kmer;
                /*
                 * Lets change into standard form
                 */
                if (_is_standard)
                    kmer.standard();
                unordered_map<Kmer, pair<size_t,size_t>>::const_iterator place =
                                                                                 _kmers_map.find(kmer);
                if (place != _kmers_map.end()) {
                    _kmers_map[kmer].first++;
                    _kmers_map[kmer].second = min(_kmers_map[kmer].second,kmer_r.kmer_pos);
                    if (_kmers_map[kmer].first == Parameters::get().accumulative_h)
                    {
                        /*
                         * First Version adding both forward and revComp
                         */
                        if (_kmers_map[kmer].second < Parameters::get().kmerSize / 2) {
                            first = true;
                            _heads.emplace(kmer_r);
                        }
                        tail = kmer_r;
                        _insert(kmer, kmer_r.kmer.rc());
                    }
                } else
                    _kmers_map[kmer] = pair<size_t,size_t>(1,kmer_r.kmer_pos);
                if (_kmers_map[kmer].first == Parameters::get().accumulative_h) {
                    /*
                     * First Version adding both forward and revComp
                     */
                    _insert(kmer, kmer_r.kmer.rc());
                }
            }
        }
        Progress::update(_sc.getIndex().size());
        std::cout << "Size Map: "<<_kmers_map.size()<<" Size Solid Kmers(as Edges): "<<_dbg_naive.size()
                  <<" Size Nodes Graph: "<<_dbg_nodes.size()<<"\n";
        _kmers_map.clear();
        /*for (auto &k: _dbg_naive)
            cout << "KmerNaive: "<<k.str()<<"\n";
        for (auto &k: _dbg_nodes) {
            cout << "KmerNodes: " << k.str() << "\n";
            vector<Kmer> neigh = getKmerNeighbors(k);
            for (auto n: neigh)
                cout << "Vecinos: "<<n.str() << "\n";
        }*/
    //sleep(10000);
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
template <>
void NaiveDBG<true>::_insert(Kmer kmer, Kmer rc)
{
    _dbg_naive.emplace(kmer);
    _dbg_nodes.emplace(kmer.substr(0,Parameters::get().kmerSize-1));
    _dbg_nodes.emplace(kmer.substr(1,Parameters::get().kmerSize));
    _dbg_nodes.emplace(rc.substr(0,Parameters::get().kmerSize-1));
    _dbg_nodes.emplace(rc.substr(1,Parameters::get().kmerSize));
}

template<>
void NaiveDBG<true>::_kmerCount() {
    std::cout << "A iterar\n";
    Pair_Kmer p_kmer;
    for (auto &read:_sc.getIndex())
    {
        Progress::update(read.first.getId());
        if ((read.first.getId() % 4) < 2)
        {
            for (auto k: IterKmers<true>(_sc.getSeq(read.second.getId()),_sc.getSeq(read.second.getPairId())))
            {
                /*
                 * Kmers pre_standar for dbg
                 */
                pair<Kmer,Kmer> nonstd_kmers = k.pair_kmer.getKmers();
                if (_is_standard)
                    k.pair_kmer.standard();
                std::cout << k.str() <<"\n";
                pair<Kmer,Kmer> kmers = k.pair_kmer.getKmers();
                std::cout << nonstd_kmers.first.str()<<"|"<<nonstd_kmers.second.str()<<"\n";
                if (_kmers_map.find(kmers.first) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.first];
                    /*
                     * Checking if we are above the threshold
                     */
                    if (++local_pair.first == Parameters::get().accumulative_h)
                        _insert(nonstd_kmers.first, nonstd_kmers.first.rc());
                    local_pair.second = min(local_pair.second,k.kmer_pos);
                    _kmers_map[kmers.first] = local_pair;
                }else {
                    _kmers_map[kmers.first].first = 0;
                    _kmers_map[kmers.first].second = k.kmer_pos;
                }
                if (_kmers_map.find(kmers.second) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.second];
                    /*
                     * Checking if we are above the threshold
                     */
                    if (++local_pair.first == Parameters::get().accumulative_h)
                        _insert(nonstd_kmers.second, nonstd_kmers.second.rc());
                    local_pair.second = min(local_pair.second, k.kmer_pos);
                }else {
                    _kmers_map[kmers.second].first = 0;
                    _kmers_map[kmers.second].second = k.kmer_pos;
                }
            }
        }
    }
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
