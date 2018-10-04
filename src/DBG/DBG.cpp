#include "DBG.h"
/*
 * Single_end reads
 */
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
                        _insert(kmer,kmer_r.kmer);
                    }
                } else
                    _kmers_map[kmer] = pair<size_t,size_t>(1,kmer_r.kmer_pos);
                if (_kmers_map[kmer].first == Parameters::get().accumulative_h) {
                    /*
                     * First Version adding both forward and revComp
                     */
                    _insert(kmer, kmer_r.kmer);
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
template<>
void NaiveDBG<true>::_kmerCount() {
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
                pair<Kmer,Kmer> kmers = k.pair_kmer.getKmers();
                if (_kmers_map.find(kmers.first) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.first];
                    /*
                     * Checking if we are above the threshold
                     */
                    if (++local_pair.first == Parameters::get().accumulative_h)
                        _insert(kmers.first, nonstd_kmers.first);
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
                        _insert(kmers.second, nonstd_kmers.second);
                    local_pair.second = min(local_pair.second, k.kmer_pos);
                }else {
                    _kmers_map[kmers.second].first = 0;
                    _kmers_map[kmers.second].second = k.kmer_pos;
                }
            }
        }
    }
    std::cout << "Size Map: "<<_kmers_map.size()<<" Size Solid Kmers(as Edges): "<<_dbg_naive.size()
              <<" Size Nodes Graph: "<<_dbg_nodes.size()<<"\n";
    /*for (auto &k: _dbg_naive)
        cout << "KmerNaive: "<<k.str()<<"\n";
    for (auto &k: _dbg_nodes) {
        cout << "KmerNodes: " << k.str() << "\n";
        vector <Kmer> neigh = getKmerNeighbors(k);
        for (auto n: neigh)
            cout << "Vecinos: " << n.str() << "\n";
    }*/
}

template<>
void NaiveDBG<true>::show_info()
{

}