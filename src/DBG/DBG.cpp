#include "DBG.h"
/*
 * Single_end reads
 */
/*
 * Number of in_edges
 */
template<>
size_t NaiveDBG<false>::in_degree(Node k)
{
    size_t out = 0;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(k.str());
        kmer_aux.appendLeft(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
/*
 * Number of out_neighbors
 */
template<>
size_t NaiveDBG<false>::out_degree(Node k)
{
    size_t out = 0;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i){
        kmer_aux = Kmer(k.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
/*
 * Get the Nt in the edges
 */
template<>
vector<DnaSequence::NuclType> NaiveDBG<false>::getNeighbors
        (Node kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<DnaSequence::NuclType> nts;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(i);
    }
    return nts;
}
/*
 * Get the neighbor k-mers
 */
template<>
vector<NaiveDBG<false>::Node> NaiveDBG<false>::getKmerNeighbors
        (Node kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<Kmer> nts;
    Node kmer_aux;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());

        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(kmer_aux.substr(1,Parameters::get().kmerSize));
    }
    return nts;
}

template<>
void NaiveDBG<false>::_kmerCount() {
        Node kmer;
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

//TODO: Check Standard
template<>
void NaiveDBG<false>::_remove_isolated_nodes()
{
    bool change = false, in_0_erase = true;
    vector<Node> erase;
    size_t cont_2 = 0;
    for (auto kmer:_dbg_nodes) {
        size_t cont = 0;
        size_t in_nodes_fw = in_degree(kmer),out_nodes_fw = out_degree(kmer);
        size_t in_nodes_total = in_nodes_fw;
        size_t len = 1;
        /*
         * InDegree = 0
         */
        if (!in_nodes_total) {
            std::cout << "Kmers con cero indegree: "<<kmer.str() << "\n";
            cont ++;
            cont_2++;
            vector<Node> aux = {kmer};
            _check_forward_path(len,aux);
            in_0_erase = _asses(erase,aux,len);
        }
        if (!in_0_erase)
            _in_0.push_back(kmer);
        /*
         * Unbalanced Nodes
         */
        if (!cont) {
            if (out_nodes_fw > in_nodes_fw)
                _in_0.push_back(kmer);
        }
        in_0_erase = true;
        /*
         * Check FWNeighbors and RCNeighbors (if proceeds)
         */
        vector<Node> neighbors = getKmerNeighbors(kmer);
        size_t num_neighbors = neighbors.size();
        size_t cont_fake_branches = 0;
        if ( num_neighbors > 1)
        {
            /*
             * Check branches
             */
            for (auto sibling:neighbors)
            {
                len = 1;
                vector<Node> aux = {sibling};
                _check_forward_path(len,aux);
                if (_asses(erase,aux,len)) {
                    cont_fake_branches++;
                }
            }
        }
    }
    if (erase.size() > 0) {
        change = true;
        _erase(erase);
    }
    /*
     * We have to iterate until convergence
     */
    if (change) {
        _in_0.clear();
        _remove_isolated_nodes();
    }else {
        cout << "KmerSolids: " << _dbg_nodes.size() << "; Suspicious Starts: " << _in_0.size() << "\n";
        cout << "Extra info:\n";
        _extra_info.show_info();
        for (auto k:_in_0)
            cout << "KmerSuspicious: " << k.str() << "\n";
    }
    /*for (auto k:_dbg_nodes)
        cout << "KmerNodes: "<<k.str()<<"\n";
    for (auto k:_dbg_naive)
        cout << "KmerSolidos: "<<k.str()<<"\n";*/
}

/*
 * Paired_end reads
 */
/*
 * Process ExtraInfo
 */
template<>
void NaiveDBG<true>::_insert_extra_info()
{
    for (auto &read:_sc.getIndex())
    {
        Progress::update(read.first.getId());
        if ((read.first.getId() % 4) < 2)
        {
            for (auto k: IterKmers<true>(_sc.getSeq(read.second.getId()), _sc.getSeq(read.second.getPairId())))
            {
                FuncNode nonstd_pk = k.pair_kmer;
                pair<Node, Node> sep_nodes = nonstd_pk.getKmers();
                if (is_solid(sep_nodes.first))
                {
                    Node k_left_rc = sep_nodes.first.rc(), k_right_rc = sep_nodes.second.rc();
                    Node origin = sep_nodes.first.substr(0, Parameters::get().kmerSize-1),
                            target = sep_nodes.first.substr(1,Parameters::get().kmerSize),
                            origin_right = sep_nodes.second.substr(0, Parameters::get().kmerSize-1),
                            target_right = sep_nodes.second.substr(1,Parameters::get().kmerSize);
                    Node origin_rc = k_left_rc.substr(0, Parameters::get().kmerSize-1),
                            target_rc = k_left_rc.substr(1, Parameters::get().kmerSize),
                            origin_right_rc = k_right_rc.substr(0, Parameters::get().kmerSize-1),
                            target_right_rc = k_right_rc.substr(1, Parameters::get().kmerSize);
                    /*
                     * TODO: Change to chunck insertion
                     */
                    /*std::cout << "Insertions: "<<"\n";
                    std::cout << "From: "<<origin.str()<<" To:"<<origin_right.str()<<"\n";
                    std::cout << "From: "<<target.str()<<" To:"<<target_right.str()<<"\n";
                    std::cout << "From: "<<target_right_rc.str()<<" To:"<<target_rc.str()<<"\n";
                    std::cout << "From: "<<origin_right_rc.str()<<" To:"<<origin_rc.str()<<"\n";*/
                    _extra_info.insert(origin, origin_right);
                    _extra_info.insert(target, target_right);
                    _extra_info.insert(target_right_rc, target_rc);
                    _extra_info.insert(origin_right_rc, origin_rc);
                }
            }
        }
    }
}
/*
  * Number of in_edges
  */
template<>
size_t NaiveDBG<true>::in_degree(Node k)
{
    size_t out = 0;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(k.str());
        kmer_aux.appendLeft(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
/*
 * Number of out_neighbors
 */
template<>
size_t NaiveDBG<true>::out_degree(Node k)
{
    size_t out = 0;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i){
        kmer_aux = Kmer(k.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            out++;
    }
    return out;
}
/*
 * Get the Nt in the edges: using pair_end constrains
 */
template<>
vector<DnaSequence::NuclType> NaiveDBG<true>::getNeighbors
        (Node kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<DnaSequence::NuclType> nts;
    Node kmer_aux;
    for (DnaSequence::NuclType i = 0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(i);

    }
    return nts;
}
/*
 * Get the neighbor k-mers: using pair_end constrains
 */
template<>
vector<NaiveDBG<true>::Node> NaiveDBG<true>::getKmerNeighbors
        (Node kmer) const
{
    if (kmer.length() == Parameters::get().kmerSize)
        kmer = kmer.substr(1, Parameters::get().kmerSize);
    vector<Node> nts;
    Node kmer_aux;
    for (DnaSequence::NuclType i=0; i < 4; ++i) {
        kmer_aux = Kmer(kmer.str());
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(kmer_aux.substr(1,Parameters::get().kmerSize));
    }
    return nts;
}

template<>
void NaiveDBG<true>::_kmerCount()
{
    Pair_Kmer p_kmer;
    for (auto &read:_sc.getIndex())
    {
        Progress::update(read.first.getId());
        if ((read.first.getId() % 4) < 2)
        {
            for (auto k: IterKmers<true>(_sc.getSeq(read.second.getId()),_sc.getSeq(read.second.getPairId())))
            {
                pair<Node,Node> nonstd_pk = k.pair_kmer.getKmers();
                /*
                 * Kmers pre_standar for dbg
                 */
                if (_is_standard)
                    k.pair_kmer.standard();
                pair<Node,Node> kmers = k.pair_kmer.getKmers();
                if (_kmers_map.find(kmers.first) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.first];
                    /*
                     * Checking if we are above the threshold
                     */
                    if (++local_pair.first == Parameters::get().accumulative_h)
                        _insert(kmers.first, nonstd_pk.first);
                    local_pair.second = min(local_pair.second,k.kmer_pos);
                    _kmers_map[kmers.first] = local_pair;
                }else
                    _kmers_map[kmers.first]= pair<size_t,size_t>(0,k.kmer_pos);
                if (_kmers_map.find(kmers.second) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.second];
                    /*
                     * Checking if we are above the threshold
                     */
                    if (++local_pair.first == Parameters::get().accumulative_h)
                        _insert(kmers.second, nonstd_pk.second);
                    local_pair.second = min(local_pair.second, k.kmer_pos);
                    _kmers_map[kmers.second] = local_pair;
                }else
                    _kmers_map[kmers.second] = pair<size_t,size_t>(0, k.kmer_pos);
            }
        }
    }
    std::cout << "Size Map: "<<_kmers_map.size()<<" Size Solid Kmers(as Edges): "<<_dbg_naive.size()
              <<" Size Nodes Graph: "<<_dbg_nodes.size()<<"\n";
    _kmers_map.clear();
    /*
     * Insert pair_end info from reads
     */
    _insert_extra_info();
    //_extra_info.show_info();
    /*for (auto &k: _dbg_naive)
        cout << "KmerNaive: "<<k.str()<<"\n";
    for (auto &k: _dbg_nodes) {
        cout << "KmerNodes: " << k.str() << "\n";
        vector <Kmer> neigh = getKmerNeighbors(k);
        for (auto n: neigh)
            cout << "Vecinos: " << n.str() << "\n";
    }*/
}

/*
 * Nodes updating:
 *      * Pruning the graph
 *      * Updating using pair_end information
 */
template<>
void NaiveDBG<true>::_to_pair_end()
{
    std::cout << "Not supported yet\n";
}
/*
 * Pruning methods
 */
template<>
void NaiveDBG<true>::_erase(vector<Node>& kmer_to_erase)
{
    for (auto kmer_erase:kmer_to_erase) {
        _dbg_nodes.erase(kmer_erase);
        for (uint i = 0; i < 8; i++) {
            Node new_kmer = kmer_erase;
            (i/4)?new_kmer.appendRight(i%4):new_kmer.appendLeft(i%4);
            if (_is_standard)
                new_kmer.standard();
            _dbg_naive.erase(new_kmer);
        }
        _extra_info.erase(kmer_erase);
    }
    /*cout << "NaiveSizePost: "<<_dbg_naive.size() << "\n";
    cout << "NodesSize: "<<_dbg_nodes.size()<<"\n";*/
    kmer_to_erase.clear();
}
template<>
void NaiveDBG<true>::_remove_isolated_nodes()
{
    bool change = false, in_0_erase = true;
    vector<Node> erase;
    size_t cont_2 = 0;
    for (auto kmer:_dbg_nodes) {
        size_t cont = 0;
        size_t in_nodes_fw = in_degree(kmer),out_nodes_fw = out_degree(kmer);
        size_t in_nodes_total = in_nodes_fw;
        size_t len = 1;
        /*
         * InDegree = 0
         */
        if (!in_nodes_total) {
            std::cout << "Kmers con cero indegree: "<<kmer.str() << "\n";
            cont ++;
            cont_2++;
            vector<Node> aux = {kmer};
            _check_forward_path(len,aux);
            in_0_erase = _asses(erase,aux,len);
        }
        if (!in_0_erase)
            _in_0.push_back(kmer);
        /*
         * Unbalanced Nodes
         */
        if (!cont) {
            if (out_nodes_fw > in_nodes_fw)
                _in_0.push_back(kmer);
        }
        in_0_erase = true;
        /*
         * Check FWNeighbors and RCNeighbors (if proceeds)
         */
        vector<Node> neighbors = getKmerNeighbors(kmer);
        size_t num_neighbors = neighbors.size();
        size_t cont_fake_branches = 0;
        if ( num_neighbors > 1)
        {
            /*
             * Check branches
             */
            for (auto sibling:neighbors)
            {
                len = 1;
                vector<Node> aux = {sibling};
                _check_forward_path(len,aux);
                if (_asses(erase,aux,len)) {
                    cont_fake_branches++;
                }
            }
        }
    }
    if (erase.size() > 0) {
        change = true;
        _erase(erase);
    }
    /*
     * We have to iterate until convergence
     */
    if (change) {
        _in_0.clear();
        _remove_isolated_nodes();
    }else {
        cout << "Size Solid kmers: "<<_dbg_naive.size()<< " Num nodes: "
             << _dbg_nodes.size() << " Suspicious Starts: " << _in_0.size() << "\n";
        cout << "Extra info:\n";
        _extra_info.show_info();
        for (auto k:_in_0)
            cout << "KmerSuspicious: " << k.str() << "\n";
        _to_pair_end();
    }
    /*for (auto k:_dbg_nodes)
        cout << "KmerNodes: "<<k.str()<<"\n";
    for (auto k:_dbg_naive)
        cout << "KmerSolidos: "<<k.str()<<"\n";*/
}

template<>
void NaiveDBG<true>::show_info()
{

}