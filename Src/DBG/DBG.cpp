#include "DBG.h"
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
        kmer_aux = Kmer(kmer);
        kmer_aux.appendRight(i);
        if (is_solid(kmer_aux))
            nts.push_back(i);
    }
    return nts;
}


template<>
void NaiveDBG<false>::show_info()
{
    for (auto p_k: _dbg_nodes)
    {
        vector<Node> neigh = getKmerNeighbors(p_k);
        cout << "Kmer: "<<p_k.str()<<"\n";
        for (auto n: neigh)
        {
            cout << "Neighbor: "<<n.str()<<"\n";
        }
    }
}

//Own->KmerCounting
template<>
void NaiveDBG<false>::_kmerCount() {
    Node kmer;
    KmerInfo<false> tail;
    bool first = false;
    size_t max_freq = 1;
    for (auto &read:_sc.getIndex()){
        if (read.first.getId() % 2)
            continue;
        if (first)
            _tails.emplace(tail);
        Progress::update(read.first.getId());
        first = false;
        for (auto kmer_r: IterKmers<false>(read.second.sequence)) {
            kmer = kmer_r.kmer;
            if (_is_standard)
                kmer.standard();
            if (_kmers_map.find(kmer) != _kmers_map.end()){
                max_freq = (++_kmers_map[kmer].first > max_freq)?_kmers_map[kmer].first:max_freq;
                _kmers_map[kmer].second = min(_kmers_map[kmer].second,kmer_r.kmer_pos);
            } else
                _kmers_map[kmer] = pair<size_t,size_t>(1,kmer_r.kmer_pos);
        }
    }
    Progress::update(_sc.getIndex().size());

    _buildGraphRepresentation(max_freq);
}
/*
 * Process ExtraInfo
 */
template<>
pair<bool, DBG<true>::Parent_Freq_Map > NaiveDBG<true>::getExtra(Node node)
{
    unordered_set<Node>::iterator it = _dbg_nodes.find(node);
    pair<bool, FreqMap> full_info = _extra_info.getInfoNode(it);
    return pair<bool, FreqMap>(full_info.first,full_info.second);
}

template<>
void NaiveDBG<true>::_insert_extra_info()
{
    Progress::get().size_total = _sc_paired.size();
    cout << "STAGE: Adding extra info\n";
    size_t adds = 0, num_read = 0;
    //Parameters::get().kmerSize -= 1;
    #pragma omp parallel
    {
        #pragma omp single
        for (auto &read:_sc_paired.getIndex()) {
            //Progress::update(read.first.getId());
            #pragma omp task shared(read)
            {
                Progress::update(num_read++);
                bool check = ((read.first.getId() % 4) == 0);
                if (check && ((read.first.getId() % 4) < 2)) {
                    unordered_set<Node>::iterator prev_node, prev_node_pair;
                    for (auto k: IterKmers<true>(read.second.sequence, _sc_paired.getSeq(read.second.getPairId()))) {
                        FuncNode nonstd_pk = k.pair_kmer;
                        pair <Node, Node> sep_nodes = nonstd_pk.getKmers();
                        if (!is_solid(sep_nodes.first) || !is_solid(sep_nodes.second))
                            continue;
                        pair <Kmer, Kmer> nodes_left = sep_nodes.first.preffixsuffix(),
                                nodes_right = sep_nodes.second.preffixsuffix();
                        unordered_set<Node>::iterator node_f_it = _dbg_nodes.find(nodes_left.first),
                                node_s_it = _dbg_nodes.find(nodes_right.first);
                        if (node_f_it == node_s_it)
                            continue;
                        if ((node_f_it != _dbg_nodes.end()) && (node_s_it != _dbg_nodes.end())) {
                            #pragma omp critical(update)
                            {
                                _extra_info.insert(node_f_it, node_s_it, adds);
                                _node_reads[(*node_s_it)].emplace(read.first.getId());
                            };
                        }
                    }
                }
            };
        }
    };
    //Parameters::get().kmerSize+=1;
    //_extra_info.show(10);
    Progress::update(_sc_paired.size());
    cout << "STAGE (inner): Plotting pairs distribution\n";
    _extra_info.show_distribution();
    cout << "Number of adds: "<<adds<<"\n";
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

template<>
void NaiveDBG<true>::_kmerCount()
{
    size_t max_freq = 1;
    Pair_Kmer p_kmer;
    for (auto &read:_sc.getIndex())
    {
        Progress::update(read.first.getId());
        if ((read.first.getId() % 4) < 2)
        {
            //Not complement info yet
            if (read.first.getId() % 4 == 1)
                continue;
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
                    local_pair.second = min(local_pair.second,k.kmer_pos);
                    max_freq = (local_pair.first++ > max_freq)?local_pair.first:max_freq;
                    _kmers_map[kmers.first] = local_pair;
                }else {
                    _kmers_map[kmers.first] = pair<size_t, size_t>(1, k.kmer_pos);
                }
                /*
                 * Second read -> kmers.first | ->kmers.second<-
                 */
                if (_kmers_map.find(kmers.second) != _kmers_map.end()) {
                    pair<size_t,size_t> local_pair = _kmers_map[kmers.second];
                    local_pair.second = min(local_pair.second, k.kmer_pos);
                    max_freq = (local_pair.first++ > max_freq)?local_pair.first:max_freq;
                    _kmers_map[kmers.second] = local_pair;
                }else {
                    _kmers_map[kmers.second] = pair<size_t, size_t>(1, k.kmer_pos);
                }
            }
        }
    }
    Progress::update(_sc.getIndex().size());
    std::cout << "Total Number of Bases: "<<_sc.getTotalBases()<<"\n";
    std::cout << "Average length read: "<<_sc.getAvLength()<<"\n";
    std::cout << "Total Kmers in all Reads: "<<_kmers_map.size()<<"\n";
    vector<size_t> histogram = getHistogram<Node,size_t>(_kmers_map, max_freq);
    Parameters::get().accumulative_h = Parameters::calculateAccumulativeParam(histogram, _sc.getTotalBases(),_sc.getAvLength());
    std::cout << "Threshold: "<<Parameters::get().accumulative_h<<"\n";
    for (auto kmer:_kmers_map)
    {
        if (kmer.second.first >= Parameters::get().accumulative_h)
        {
            _insert(kmer.first, kmer.first, kmer.second.first);
        }
    }
    std::cout<<"Total Solid K-mers(Graph Edges): "<<_dbg_naive.size()
             <<" Total Graph Nodes: "<<_dbg_nodes.size()<<"\n";
    _kmers_map.clear();
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
    /*for (auto node_kmer: _dbg_nodes)
    {
        vector<Node> neighbors = getKmerNeighbors(node_kmer);
        vector<unordered_set<Node>> neighbors_couples = getNeighborsCouples(node_kmer);
        if (_extra_info.find(node))
        {
            unordered_set<Node> couples = _extra_info[node];

        }
    }*/
}

template<>
void NaiveDBG<true>::show_info()
{
    for (auto p_k: _dbg_nodes)
    {
        vector<Node> neigh = getKmerNeighbors(p_k);
        cout << "Kmer: "<<p_k.str()<<"\n";
        for (auto n: neigh)
        {
            cout << "Neighbor: "<<n.str()<<"\n";
        }
    }
}

/*
 * Boost graphs -> Pair-end Reads
 */
template<>
void  boostDBG<true>::_insert_extra(SequenceContainer * sc, unordered_set<Node> * solids, unordered_set<Node> * nodes)
{
    Progress::get().size_total = sc->size();
    cout << "STAGE: Adding extra info (boost-time)\n";
    //Greedy solution :S
    Parameters::get().kmerSize-=1;
    size_t adds = 0, num_read = 0;
    unordered_set<uint32_t> ids;
    if (Parameters::get().remove_duplicates)
    {
        size_t numThreads = Parameters::get().numThreads;
        //Fichero1 -> left_reads.fa
        sc->write_left_chains();
        System::execute("bash -c \"sga index -t "+std::to_string(numThreads)+" left_reads.fa\"");
        System::execute("bash -c \"sga rmdup -t "+std::to_string(numThreads)+" left_reads.fa\"");
        ids = Parse::getIdsFromFiles("left_reads.rmdup.fa");
        //Fichero2 -> rigth_reads.fa
        sc->write_right_chains();
        System::execute("bash -c \"sga index -t "+std::to_string(numThreads)+" right_reads.fa\"");
        System::execute("bash -c \"sga rmdup -t "+std::to_string(numThreads)+" right_reads.fa\"");
        for (auto id: Parse::getIdsFromFiles("right_reads.rmdup.fa"))
        {
            ids.emplace(id);
        }
        cout << "number of reads: "<<sc->size()<<" NumberIds: "<<ids.size()<<endl;
    }
    #pragma omp parallel
    {
        #pragma omp single
        for (auto &read:sc->getIndex()) {
            //Progress::update(read.first.getId());
            #pragma omp task shared(read)
            {
                Progress::update(num_read++);
                bool check = ((read.first.getId() % 4) == 0);
                if (check && ((read.first.getId() % 4) < 2)) {
                    bool dup_check = (Parameters::get().remove_duplicates)?(ids.find(read.first.getId())!=ids.end()):true;
                    if (dup_check) {
                        unordered_set<Node>::iterator prev_node, prev_node_pair;
                        for (auto k: IterKmers<true>(read.second.sequence, sc->getSeq(read.second.getPairId()))) {
                            FuncNode nonstd_pk = k.pair_kmer;
                            pair <Node, Node> sep_nodes = nonstd_pk.getKmers();
                            unordered_set<Node>::iterator node_f_it = nodes->find(sep_nodes.first),
                                    node_s_it = nodes->find(sep_nodes.second);
                            if (node_f_it == node_s_it)
                                continue;
                            /*
                             * Forward
                             */
                            if ((node_f_it != nodes->end()) && (node_s_it != nodes->end())) {
                                graphBU node_bu_l = local_map[(*node_f_it)];
                                graphBU rep_node = _representants_map[local_map[(*node_s_it)]];
                                unordered_set<Node>::iterator representant = nodes->find(_g[rep_node].node);
                                #pragma omp critical(update)
                                {
                                    if (_representants_hits[rep_node].find(node_bu_l) ==
                                        _representants_hits[rep_node].end())
                                        _representants_hits[rep_node][node_bu_l] = 1;
                                    else
                                        _representants_hits[rep_node][node_bu_l]++;
                                    if (_representants_hits[rep_node][node_bu_l] == LIMIT_TO_REP) {
                                        _extra_info.insert(node_f_it, representant, adds);
                                        _node_reads[(*node_s_it)].emplace(read.first.getId());
                                    }
                                };
                            }

                            /*
                             * Reverse
                             */
                            if (Parameters::get().full_info) {
                                unordered_set<Node>::iterator node_f_it_rc = nodes->find(sep_nodes.first.rc()),
                                        node_s_it_rc = nodes->find(sep_nodes.second.rc());
                                if ((node_f_it_rc != nodes->end()) && (node_s_it_rc != nodes->end())) {
                                    graphBU node_bu_r = local_map[(*node_s_it_rc)];
                                    graphBU rep_node = _representants_map[local_map[(*node_f_it_rc)]];
                                    unordered_set<Node>::iterator representant = nodes->find(_g[rep_node].node);
                                    #pragma omp critical(update)
                                    {
                                        if (_representants_hits[rep_node].find(node_bu_r) ==
                                            _representants_hits[rep_node].end())
                                            _representants_hits[rep_node][node_bu_r] = 1;
                                        else
                                            _representants_hits[rep_node][node_bu_r]++;
                                        if (_representants_hits[rep_node][node_bu_r] == LIMIT_TO_REP) {
                                            _extra_info.insert(node_s_it_rc, representant, adds);
                                            _node_reads[(*node_s_it_rc)].emplace(read.first.rc().getId());
                                        }
                                    };
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Progress::update(sc->size());
    Parameters::get().kmerSize+=1;
}

template<>
int* boostDBG<true>::_floyds_warshall()
{
    std::cout << "STAGE: Building warshall matrix\n";
    size_t num_vertex = boost::num_vertices(_g);
    /*
     * Free!!!
     */
    int *dist = (int*) malloc(num_vertex*num_vertex*sizeof(int));
    /*
     * Initialize the dist_matrix
     */
    vertex_iterator v, vend;
    for (boost::tie(v,vend) = boost::vertices(_g); v != vend; ++v)
    {
        for (size_t i = 0; i < num_vertex; ++i)
        {
            if (_g[*v].id != (int)i)
                dist[_g[*v].id*num_vertex+(int)i] = INF;
            else
                dist[_g[*v].id*num_vertex+(int)i] = 1;
        }
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges((*v), _g);
        for(; neighbors.first != neighbors.second; ++neighbors.first)
        {
            dist[_g[*v].id*num_vertex+_g[boost::target(*neighbors.first,_g)].id] = 1;
        }
    }
    /*
     * Floyd Warshall
     */
    for (size_t k = 0; k < num_vertex; ++k) {
        for (size_t i = 0; i < num_vertex; ++i) {
            if (dist[i * num_vertex + k] > 2*DELTA_PATH_LEN)
                continue;
            for (size_t j = 0; j < num_vertex; ++j) {
                if ((dist[i * num_vertex + k] + dist[k * num_vertex + j]) < dist[i * num_vertex + j])
                    dist[i * num_vertex + j] = dist[i * num_vertex + k] + dist[k * num_vertex + j];
            }
        }
    }
    /*
     * Little print
     */
    for (size_t k = 0; k < num_vertex; ++k)
    {
        for (size_t i = 0; i < num_vertex; ++i)
            std::cout << dist[k * num_vertex + i] << " ";
        std::cout << "\n";
    }

    return dist;
}
template<>
void boostDBG<true>::_full_fil_matrix(size_t distance, size_t num_vertex, size_t executions,
                               graphBU source, graphBU cur_node, size_t mode, vector<bool> & reach_matrix,
                                      bool * token, unordered_set<graphBU> & check_store)
{
    /*
     * Recursion peak if we add all the inter-calls
     */
    size_t max_distance = (mode)?(mode==1)?(2*DELTA_PATH_LEN):(2*SHORT_LENGTH):(2*POLISH_PATH_LEN);
    if (distance >= max_distance){
        return;
    }
    distance++;
    pair<out_iterator, out_iterator> neighbors =
            boost::out_edges(cur_node, _g);
    for(; neighbors.first != neighbors.second; ++neighbors.first)
    {
        executions++;
        graphBU tmp = boost::target(*neighbors.first,_g);
        if (tmp == source || check_store.find(tmp) != check_store.end())
            continue;
        check_store.emplace(tmp);
        if (translator_vector[_g[tmp].id] < INF) {
            if ((translator_vector[_g[source].id] * num_vertex + translator_vector[_g[tmp].id]) > reach_matrix.size()) {
                cout << "Source: "<<translator_vector[_g[source].id]<<" "<<"Tmp: "<<translator_vector[_g[tmp].id]<<endl;
                cout << translator_vector[_g[source].id] * num_vertex + translator_vector[_g[tmp].id]<< "SizeMatriz: "<<reach_matrix.size()<<endl;
                cout << "Fallo del tamanho" << endl;
                exit(1);
            }
            reach_matrix[translator_vector[_g[source].id] * num_vertex + translator_vector[_g[tmp].id]] = true;
        }
        _full_fil_matrix(distance, num_vertex,executions ,source, tmp, mode, reach_matrix, token, check_store);
    }
}
template<>
bool boostDBG<true>::_reachable(int * dm, size_t row, size_t col)
{
    return (dm[row*boost::num_vertices(_g)+col] < (2*DELTA_PATH_LEN));
}
template<>
bool boostDBG<true>::_reachable(graphBU oriSource,
                                graphBU source,
                                graphBU target,
                                size_t distance,
                                size_t * branches,
                                vector<bool> & reach_vect,
                                vector<bool> & check_vect,
                                size_t num_vertex)
{
    if (check_vect[_g[target].id*num_vertex+_g[source].id] || check_vect[_g[source].id*num_vertex+_g[target].id] )
        return reach_vect[_g[source].id*num_vertex+_g[target].id];
    if (!(*branches))
        return false;
    if (distance >= (2*DELTA_PATH_LEN)) {
        (*branches)--;
        return false;
    }
    if (target == source) {
        (*branches)--;
        return true;
    }
    vector<graphBU> neigh = getKmerNeighbors(source);
    bool out = false;
    distance++;
    for (auto n: neigh)
    {
        reach_vect[_g[oriSource].id*num_vertex+_g[n].id] = true;
        check_vect[_g[oriSource].id*num_vertex+_g[n].id] = true;
        check_vect[_g[n].id*num_vertex+_g[oriSource].id] = true;
        reach_vect[_g[source].id*num_vertex+_g[target].id] = true;
        check_vect[_g[target].id*num_vertex+_g[source].id] = true;
        check_vect[_g[source].id*num_vertex+_g[target].id] = true;
        out |= _reachable(oriSource, n, target, distance, branches, reach_vect, check_vect, num_vertex);
        if (out)
            return out;
        (*branches)--;
    }
    return out;
}

template<>
void boostDBG<true>::_modify_info()
{
    /*
     * New graph definition:
     *      - BidirectionalS/UndirectedS
     */
    struct NodeInfoLocal {
        NodeInfoLocal():id(-1){}
        NodeInfoLocal(graphBU node, int32_t id):node(node), id(id){}
        NodeInfoLocal& operator=(const NodeInfoLocal& other)
        {
            this->node = other.node;
            this->id = other.id;
            return *this;
        }
        bool empty()
        {
            return (id == -1);
        }
        bool equal(graphBU node) const
        {
            return node == this->node;
        }
        bool operator ==(const NodeInfoLocal &other) const
        {
            return equal(other.node) && (this->id == other.id);
        }
        graphBU node;
        int32_t id;
    };
    typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, NodeInfo> Graph_l;
    typedef Graph_l::vertex_descriptor vertex_graph;
    typedef Graph_l::vertex_iterator vertex_it;
    /*
     * We are going to use pair_info (ExtraInfoNodes) to "modify" the graph including this new information
     */
    vertex_iterator v, vend;
    size_t num_vertex = getVertex(), num_representants = getNumRepresentants();
    //cout << "Num vertex: "<<num_vertex<<" Num representants: "<<getNumRepresentants()<<endl;
    int * distance_matrix;
    size_t option = 1;
    if (FLOYD)
        distance_matrix = _floyds_warshall();
    else
    {
        cout << "Inner Stage: Building distance matrix (representative to representative)"<<endl;
        reach = vector<bool>(num_representants*num_representants, false);
        //reach_by_read = vector<bool>(num_vertex * num_vertex, false);
        size_t distance, executions;
        unordered_set<graphBU> check_store;
        //Representative - Long path
        Progress::get().size_total = num_representants;
        size_t progress = 0;
        for (auto r: _first_last){
            distance = 0;
            executions = 0;
            bool token = false;
            check_store.emplace(r);
            _full_fil_matrix(distance, num_representants, executions, r, r, 1, reach, &token, check_store);
            check_store.clear();
            Progress::update(progress++);
        }
        Progress::update(num_representants);
    }
    //Progression
    cout << "STAGE: Processing cliques\n";
    Progress::get().size_total = num_vertex;
    //Incremental information:
    size_t cont = 0;
    //#pragma omp parallel for
    for (size_t v = 0; v < num_vertex; ++v){
            NodeInfo node_info = _g[v];
            Progress::update(cont++);
            //unordered_set<graphBU> node_extra = _map_extra_info[node_info.id];
            auto node_extra = getExtraInfoNode(v);
            /*cout << "Nodo (extra size): "<<node_extra.size()<<endl;
            for (auto i:node_extra)
                cout<<i<<", ";
            cout << endl;*/
            if (node_extra.empty())
                continue;
            auto endpoints = getKmerNeighbors(v);
            for (auto endpoint: endpoints)
            {
            /*pair <out_iterator, out_iterator> neighbors =
                        boost::out_edges((v), _g);
            for (; neighbors.first != neighbors.second; ++neighbors.first) {
                auto endpoint = boost::target(*neighbors.first, _g);*/
                NodeInfo neigh_info = _g[endpoint];
                size_t curr_node = 0;
                //unordered_set <graphBU> neigh_extra = _map_extra_info[neigh_info.id];
                auto neigh_extra = getExtraInfoNode(endpoint);
                /*for (auto i:neigh_extra)
                    cout<<i<<", ";
                cout << endl;*/
                if (neigh_extra.empty())
                    continue;
                unordered_set <graphBU> local_vect = getUnion(node_extra,neigh_extra),last_set, removed_representant;
                unordered_map<graphBU, vector<graphBU>> store_map;
                unordered_map<Node, graphBU> translate_map;
                unordered_map<graphBU, size_t> representant_hits;
                set<graphBU> local_set;
                size_t max_size = 0;
                if (option == 1)
                {
                    unordered_set<graphBU> to_rep;
                    unordered_set<graphBU> endpoints;
                    unordered_map<graphBU, size_t> representant_siege;
                    for (auto s:local_vect) {
                        //graphBU representative = _representants_map[s];
                        graphBU representative = getRepresentant(s);
                        //cout << "S: "<<s<<"Representant: "<<representative << endl;
                        if (store_map.find(representative) == store_map.end()) {
                            to_rep.emplace(representative);
                            (in(node_extra, s)?
                                    representant_siege[representative] = 0:representant_siege[representative] = 1);
                            removed_representant.emplace(representative);
                            store_map[representative] = vector<vertex_t>(1, s);
                            //representant_hits[representative] = _representants_hits[representative][v]+_representants_hits[representative][endpoint];
                            representant_hits[representative] = getRepresentantHitsNode(representative, v) + getRepresentantHitsNode(representative, endpoint);
                            //cout << "Hits: "<<representative<<" "<<representant_hits[representative]<<endl;
                            translate_map[_g[representative].node] = representative;
                            local_set.emplace(representative);
                        }else {
                            store_map[representative].push_back(s);
                            if (removed_representant.find(representative) != removed_representant.end())
                            {
                                if ((in(node_extra, s) && representant_siege[representative] == 1) ||
                                        (in(neigh_extra, s) && representant_siege[representative] == 0))
                                    removed_representant.erase(representative);
                            }
                        }
                    }
                    for (auto k:getKmerNeighbors(endpoint)) {
                        //for (auto s: _map_extra_info[_g[k].id]){
                        for (auto s: getExtraInfoNode(k)){
                            //graphBU representative = _representants_map[s];
                            graphBU representative = getRepresentant(s);
                            if (representant_hits.find(representative) != representant_hits.end())
                                //representant_hits[representative]+=_representants_hits[representative][k];
                                representant_hits[representative] += getRepresentantHitsNode(representative, k);
                            else{
                                to_rep.emplace(representative);
                                store_map[representative] = vector<graphBU>(1,s);
                                //representant_hits[representative] = _representants_hits[representative][k];
                                representant_hits[representative] = getRepresentantHitsNode(representative, k);
                                translate_map[_g[representative].node] = representative;
                                local_set.emplace(representative);
                            }
                        }
                    }
                    for (auto k:getInKmerNeighbors(v)) {
                        //for (auto s: _map_extra_info[_g[k].id]){
                        for (auto s: getExtraInfoNode(k)){
                            //graphBU representative = _representants_map[s];
                            graphBU representative = getRepresentant(s);
                            if (representant_hits.find(representative) != representant_hits.end())
                                //representant_hits[representative]+=_representants_hits[representative][k];
                                representant_hits[representative] += getRepresentantHitsNode(representative, k);
                            else{
                                store_map[representative] = vector<graphBU>(1,s);
                                //representant_hits[representative] = _representants_hits[representative][k];
                                translate_map[_g[representative].node] = representative;
                                representant_hits[representative] = getRepresentantHitsNode(representative, k);
                                local_set.emplace(representative);
                            }
                        }
                    }
                    //Remove representants with no more than MIN_REP element (pointless)
                    removed_representant.clear();
                    for (auto p:store_map)
                        (p.second.size() > max_size)?max_size=p.second.size():max_size=max_size;
                    for (auto p:store_map) {
                        if (p.second.size() < floor(max_size*MIN_SIZE_REP))
                            removed_representant.emplace(p.first);
                    }
                    local_vect = to_rep;
                }
                for (auto s: removed_representant) {
                    local_vect.erase(s);
                    local_set.erase(s);
                }
                //cout << "LocalVect: "<<local_vect.size()<<"\n";
                Graph_l local_graph;
                //unordered_map<Node, vertex_graph> local_node_map;
                unordered_map<size_t, vertex_graph> local_node_map;
                //cout << "NumIterations: "<<((size_t)last_set.size()*(last_set.size()-1)/2)<<endl;
                last_set = local_vect;
                unordered_map<vertex_graph, size_t> representative_size;
                for (auto s:local_set)
                {
                    NodeInfo s_info = _g[s];
                    //cout << "S: "<<s<<" "<<getRepresentantTranslation(s)<<endl;
                    vertex_graph source = boost::add_vertex(NodeInfo(s_info.node, curr_node++), local_graph);
                    representative_size[source] = representant_hits[s];
                    //local_node_map[s_info.node] = source;
                    local_node_map[getRepresentantTranslation(s)] = source;
                }
                for (auto s:last_set) {
                    NodeInfo s_info = _g[s];
                    //vertex_graph source = local_node_map[s_info.node];
                    vertex_graph source = local_node_map[getRepresentantTranslation(s)];
                    for (auto t:local_vect) {
                        if (s == t)
                            continue;
                        NodeInfo t_info = _g[t];
                        bool reached = false;
                        pair<graphBU, graphBU> s_pair, t_pair;
                        if (option == 1)
                        {
                            s_pair = _represent_first_last[s], t_pair = _represent_first_last[t];
                            /*reached = (reach[translator_vector[_g[s_pair.second].id]*num_representants+translator_vector[_g[t_pair.first].id]]
                                       || reach[translator_vector[_g[t_pair.second].id]*num_representants+translator_vector[_g[s_pair.first].id]]);*/
                            reached = (reach[getRepresentantTranslation(s_pair.second)*num_representants+getRepresentantTranslation(t_pair.first)]
                                       || reach[getRepresentantTranslation(t_pair.second)*num_representants+getRepresentantTranslation(s_pair.first)]);
                            /*if (node_info.id < 20)
                            {
                                cout << "S_pair (s): " << s_pair.first << " " << s_pair.second << " Translation: "
                                     << getRepresentantTranslation(s_pair.second) << " "
                                     << getRepresentantTranslation(s_pair.first)
                                     << " Tpair (t):" << t_pair.second << " " << t_pair.first << " Translation: "
                                     << getRepresentantTranslation(t_pair.first)
                                     << " " << getRepresentantTranslation(t_pair.second) << " " << reached << endl;
                            }*/
                        }
                        bool add = true;
                        if (reached) {
                            if (add) {
                                //vertex_graph target = local_node_map[t_info.node];
                                vertex_graph target = local_node_map[getRepresentantTranslation(t)];
                                boost::add_edge(source, target, local_graph);
                                //cout << "Adding vertex from: "<<source<<" to "<<target<<endl;
                            }
                        }
                    }
                    local_vect.erase(s);
                }
                /*if (node_info.id < 20) {
                    cout << "Node: "<<node_info.id<<endl;
                    show_graph<Graph_l, vertex_it>(local_graph);
                }*/

                /*
                 * Polish Graph
                 */
                map <vertex_graph, vector<Node>> clique_representation;
                unordered_set <vertex_graph> nodes_erase;
                /*
                 * Calculate all cliques in the local graph
                 */
                std::set<size_t> idCliques;
                priority_queue < pair < size_t, vector < vertex_graph >> > output;
                size_t num_vertex_local = boost::num_vertices(local_graph),
                        num_edges_local = boost::num_edges(local_graph);
                if (((num_vertex_local * (num_vertex_local - 1)) / 2) == (num_edges_local)) {
                    if ((node_extra.size()*PARENT_SON_LIMIT) > neigh_extra.size())
                        continue;
                    PairedInfoNode uniqueHaplotype;
                    for (auto s:last_set) {
                        if (option == 0)
                            uniqueHaplotype.emplace(_g[s].node);
                        if (option == 1) {
                            for (auto p:store_map[s])
                                uniqueHaplotype.emplace(_g[p].node);
                        }
                    }
                    #pragma omp critical(update)
                    {
                        //_g[endpoint].parent_cliques[node_info.node].push_back(uniqueHaplotype);
                        addHaplotype(endpoint, v, uniqueHaplotype);
                    }
                } else {
                    //cout << "V: "<<node_info.id<<endl;
                    output = findMaxClique < Graph_l, vertex_graph, vertex_it > (local_graph, idCliques, representative_size);
                    //cout << "Output: "<<output.size()<<endl;
                    unordered_set <vertex_graph> edge_transversed;
                    float total_size = 0.0, average_size = 0.0;
                    size_t num_clicks = 0;
                    if (output.empty()) {
                        vertex_it v2, v_end2;
                        for (tie(v2, v_end2) = boost::vertices(local_graph); v2 != v_end2; ++v2) {
                            PairedInfoNode local_haplotype;
                            if (option == 0) {
                                local_haplotype.emplace(local_graph[*v2].node);
                                auto p = clique_representation.find(*v2);
                                if (p != clique_representation.end()) {
                                    for (auto n: (*p).second)
                                        local_haplotype.emplace(n);
                                }
                            }
                            if (option == 1){
                                for (auto n: store_map[translate_map[local_graph[*v2].node]]) {
                                    local_haplotype.emplace(_g[n].node);
                                }
                            }
                            #pragma omp critical(update)
                            {
                                //_g[endpoint].parent_cliques[node_info.node].push_back(local_haplotype);
                                addHaplotype(endpoint, v, local_haplotype);
                            }
                            total_size += local_haplotype.size();
                            num_clicks++;
                        }
                    } else {
                        while (!output.empty() & (edge_transversed.size() != boost::num_vertices(local_graph))) {
                            pair <size_t, vector<vertex_graph>> top_click = output.top();
                            output.pop();
                            vector <vertex_graph> clique = top_click.second;
                            //cout << "Clique size: "<<clique.size()<<endl;
                            PairedInfoNode local_haplotype;
                            for (size_t i = 0; i < clique.size(); ++i) {
                                Node node_clique = local_graph[clique[i]].node;
                                if (option == 0)
                                {
                                    local_haplotype.emplace(node_clique);
                                    auto p = clique_representation.find(clique[i]);
                                    if (p != clique_representation.end())
                                        for (auto n: (*p).second) {
                                            local_haplotype.emplace(n);
                                        }
                                }
                                if (option == 1)
                                {
                                    for (auto n: store_map[translate_map[local_graph[clique[i]].node]]){
                                        local_haplotype.emplace(_g[n].node);
                                    }
                                }
                                edge_transversed.emplace(clique[i]);
                            }
                            #pragma omp critical(update)
                            {
                                //_g[endpoint].parent_cliques[node_info.node].push_back(local_haplotype);
                                addHaplotype(endpoint, v, local_haplotype);
                            }
                            total_size += local_haplotype.size();
                            num_clicks++;
                        }
                        if (output.empty() & (edge_transversed.size() != boost::num_vertices(local_graph))) {
                            vertex_it v2, v_end2;
                            for (tie(v2, v_end2) = boost::vertices(local_graph); v2 != v_end2; ++v2) {
                                if (edge_transversed.find(*v2) == edge_transversed.end()) {
                                    PairedInfoNode local_haplotype;
                                    if (option == 0)
                                    {
                                        auto p = clique_representation.find(*v2);
                                        local_haplotype.emplace(local_graph[*v2].node);
                                        if (p != clique_representation.end()) {
                                            for (auto n: (*p).second)
                                                local_haplotype.emplace(n);
                                        }
                                    }
                                    if (option == 1)
                                    {
                                        for (auto n: store_map[translate_map[local_graph[*v2].node]]) {
                                            local_haplotype.emplace(_g[n].node);
                                        }
                                    }
                                    #pragma omp critical(update)
                                    {
                                        //_g[endpoint].parent_cliques[node_info.node].push_back(local_haplotype);
                                        addHaplotype(endpoint, v, local_haplotype);
                                    }
                                    total_size += local_haplotype.size();
                                    num_clicks++;
                                }
                            }
                        }
                    }
                    average_size = total_size / num_clicks;
                    vector<PairedInfoNode> cliquesAvail = getCliquesWithParent(endpoint,v);
                    PairedInfoNode pairedInfoNeigh = getPairedInfoNode(endpoint),
                        pairedInfoNode = getPairedInfoNode(v);
                    size_t removes = 0, numberCliques = cliquesAvail.size();
                    vector<bool> fake_cliques(numberCliques/*_g[endpoint].parent_cliques[node_info.node].size()*/, false);
                    vector<size_t> caso(numberCliques/*_g[endpoint].parent_cliques[node_info.node].size()*/,6);
                    vector<size_t> join_cliques;
                    vector<unordered_set<Node>> intersections;
                    if (numberCliques != _g[endpoint].parent_cliques[node_info.node].size()){
                        cout << "ERROR! "<<numberCliques<<endl;
                        exit(1);
                    }
                    /*for (auto clique: _g[endpoint].parent_cliques[node_info.node])
                        intersections.push_back(getIntersection(clique, neigh_info.node_set));*/
                    for (auto clique: cliquesAvail) {
                        if (getIntersection(clique, pairedInfoNeigh) != getIntersection(clique,neigh_info.node_set))
                        {
                            cout << "ERROR! Intersecciones"<<endl;
                            exit(1);
                        }
                        intersections.push_back(getIntersection(clique, pairedInfoNeigh));
                    }
                    for (uint16_t i = 0; i < numberCliques/*_g[endpoint].parent_cliques[node_info.node].size()*/; ++i)
                    {
                        if (fake_cliques[i])
                            continue;
                        size_t num_neighs = getKmerNeighbors(v).size(),
                                clique_size = cliquesAvail[i].size()/*_g[endpoint].parent_cliques[node_info.node][i].size()*/;
                        /*unordered_set<Node> parent_none = getIntersection(_g[endpoint].parent_cliques[node_info.node][i], node_info.node_set),
                                neigh_none = intersections[i];*/
                        unordered_set<Node> parent_none = getIntersection(cliquesAvail[i],pairedInfoNode), neigh_none = intersections[i];
                        if (parent_none != getIntersection(_g[endpoint].parent_cliques[node_info.node][i], node_info.node_set)){
                            cout <<"ERROR! Parent None"<<endl;
                            exit(1);
                        }
                        if (neigh_none.empty() || parent_none.empty()) {
                            fake_cliques[i] = true;
                            if (parent_none.empty())
                                caso[i] = 7;
                            else {
                                caso[i] = 1;
                                total_size -= clique_size;
                            }
                            removes++;
                        }else {
                            if (clique_size < (average_size * CLICK_RATIO)) {
                                fake_cliques[i] = true;
                                caso[i] = 4;
                            } else {
                                    /*bool neigh_full = isSubset(_g[endpoint].parent_cliques[node_info.node][i],
                                                               neigh_info.node_set, SIMILARITY_RATIO)
                                    , parent_full = isSubset(_g[endpoint].parent_cliques[node_info.node][i],
                                                             node_info.node_set, SIMILARITY_RATIO),
                                            parent_same = isSame(parent_none, node_info.node_set, IDENTITY_RATIO);*/
                                    if (cliquesAvail[i] != _g[endpoint].parent_cliques[node_info.node][i]){
                                        cout << "ERROR! CliquesAvail"<<endl;
                                        exit(1);
                                    }
                                    bool neigh_full = isSubset(cliquesAvail[i],pairedInfoNeigh, SIMILARITY_RATIO),
                                        parent_full = isSubset(cliquesAvail[i],pairedInfoNode, SIMILARITY_RATIO),
                                        parent_same = isSame(parent_none,pairedInfoNode,IDENTITY_RATIO);
                                    if ((!neigh_full && parent_full) && !parent_same && (num_neighs > 1)) {
                                        fake_cliques[i] = true;
                                        caso[i] = 2;
                                        total_size -= clique_size;
                                        removes++;
                                    } else {
                                        for (uint16_t j = 0;
                                             j < numberCliques/*_g[endpoint].parent_cliques[node_info.node].size()*/; ++j) {
                                            if (i == j || fake_cliques[j])
                                                continue;
                                            if (isSubset(intersections[i], intersections[j]) &&
                                                (intersections[i].size() != intersections[j].size())) {
                                                fake_cliques[i] = true;
                                                caso[i] = 3;
                                                total_size -= clique_size;
                                                removes++;
                                                break;
                                            }
                                        }
                                    }
                            }
                        }
                    }
                    average_size = total_size / (num_clicks-removes);
                    for (size_t i = 0; i < fake_cliques.size(); i++)
                    {
                        if (fake_cliques[i] && caso[i] == 7){
                            if (in_degree(endpoint) == 1)
                                caso[i] = 4;
                        }
                        if (fake_cliques[i] && caso[i] == 4){
                            if ( cliquesAvail[i].size()/*_g[endpoint].parent_cliques[node_info.node][i].size()*/ > (average_size*CLICK_RATIO)) {
                                fake_cliques[i] = false;
                                caso[i] = 6;
                            }
                        }
                    }
                    size_t count = 0, count2 = 0;
                    for (auto i: fake_cliques) {
                        if (i) {
                            #pragma omp critical(update)
                            {
                                _g[endpoint].parent_cliques[node_info.node].
                                        erase(_g[endpoint].parent_cliques[node_info.node].begin() + count - count2++);
                            }
                        }
                        count++;
                    }
                    fake_cliques.clear();
                }
            }
        }
    //Free Floyd
    if (FLOYD)
        free(distance_matrix);
    else
        reach.clear();
    Progress::update(num_vertex);
    std::cout << "\n";
}
template<>
void boostDBG<true>::show_info(size_t max_it)
{
    vertex_iterator v, vend;
    size_t cont = 0;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        if (_g[*v].node_set.size() == 0)
            continue;
        if (cont++ == max_it)
            break;
        cout << " Kmer:"     << _g[*v].node.str()
                  << " id:"  << _g[*v].id
                  << " Puntero: " << (*v)
                  << " Coverage: "<<_g[*v].coverage
                  << endl;
        vector<Node> neighbors = getKmerNeighbors(_g[*v].node);
        cout << "Neighbors: ";
        for (auto n:neighbors)
            std::cout  << n.str()<<" ";
        cout << "\n Couples: ";
        for (auto n:_g[*v].node_set)
            std::cout << n.str()<<" ";
        cout << "NumCouples: "<<_g[*v].node_set.size()<<"\n";
        cout << "\n Haplotypes (with parent): ";
        for (auto n:_g[*v].parent_cliques){
            cout << "Parent: "<<n.first.str()<<"\n";
            for (auto n2:n.second) {
                for (auto k:n2)
                    cout << " " << k.str() << " ";
                cout << "\n";
            }
        }
        std::cout<<"\n";
    }
}

template<>
void boostDBG<true>::_transverse(const it_node &n,
                                 map<size_t, vector<size_t>> & graphUnitigs,
                                 map<size_t, graphBU> & map_seq_nodo_end,
                                 map<size_t, DnaSequence> & map_seqs,
                                 unordered_set<graphBU> & nodes_extended,
                                 DnaSequence & first_seq,
                                 size_t & seq_id)
{
    //Mucho ojo con este controlador
    unordered_set<it_node_set, hash_fn> rep_controller;
    stack<it_node> stack_node;
    stack<DnaSequence> stack_seq;
    stack<int> stack_prev_seqs;
    stack_node.emplace(n);
    stack_seq.emplace(first_seq);
    stack_prev_seqs.emplace(-1);
    size_t mas = 0, masSame = 0, igual = 0;
    bool show = false;
    while (!stack_node.empty()) {
        //cout << "Stack size: "<<stack_node.size()<<endl;
        it_node node_assay = stack_node.top();
        nodes_extended.emplace(node_assay.first);
        stack_node.pop();
        DnaSequence sequence = stack_seq.top();
        stack_seq.pop();
        NodeInfo nodeInfo = _g[node_assay.first];
        if (sequence.length() >= 100000)
            show = true;
        else
            show = false;
        size_t kmer_size = Parameters::get().kmerSize;
        pair <vector<it_node>, vector<it_node>> pair_neigh;
        if (!node_assay.second.empty()) {
            pair_neigh = getOutKmerNeighbors(node_assay);
        } else {
            vector <graphBU> neighs = getKmerNeighbors(node_assay.first);
            for (auto n_l:neighs)
                pair_neigh.first.push_back(it_node(n_l, PairedInfoNode()));
        }
        vector <it_node> neighbors = pair_neigh.first;
        if ((neighbors.size() == 1) /*(in_nodes == 1)*/) {
            if (rep_controller.find(it_node_set(node_assay)) == rep_controller.end()) {
                rep_controller.emplace(it_node_set(node_assay));
                sequence.append_nuc_right(nodeInfo.node.at(0));
                stack_seq.emplace(sequence);
                stack_node.emplace(neighbors[0]);
            } else {
                for (uint i = 0; i < kmer_size - 1; ++i)
                    sequence.append_nuc_right(nodeInfo.node.at(i));
                map_seq_nodo_end[seg] = node_assay.first;
                map_seqs[seg] = sequence;
                size_t seqTmp = stack_prev_seqs.top();stack_prev_seqs.pop();
                graphUnitigs[seqTmp].push_back(seg);
                graphUnitigs[seg] = vector<size_t>();
                seg++;
            }
        } else {
            if (neighbors.size() == 0) {
                igual++;
                if (show) {
                    cout << "NO NEIGHBORS\n";
                    cout << "Kmer: " << nodeInfo.node.str() << "\n";
                    cout << "Haplotype: " << endl;
                    show_set(node_assay.second);
                    cout << "Intersection: " << endl;
                    show_set(getIntersection(node_assay.second, nodeInfo.node_set));
                    cout << "Paired_info: " << endl;
                    show_set(nodeInfo.node_set);
                    cout << "Own Cliques:" << endl;
                    for (auto clique: nodeInfo.parent_cliques) {
                        cout << "Parent: " << clique.first.str() << endl;
                        for (auto c: clique.second) {
                            cout << "Clique: " << endl;
                            show_set(c);
                            cout << "Intersection: " << endl;
                            show_set(getIntersection(c, nodeInfo.node_set));
                        }
                    }
                    cout << "Neighbors Reales:" << endl;
                    for (auto n2: getKmerNeighbors(node_assay.first)) {
                        cout << "Neighbor: " << _g[n2].node.str() << endl;
                        cout << "Cliques with parent: " << endl;
                        for (auto cliques: _g[n2].parent_cliques) {
                            cout << "Parent: " << cliques.first.str() << endl;
                            if (cliques.first == nodeInfo.node.str()) {
                                for (auto click:cliques.second) {
                                    cout << "Clique: " << endl;
                                    show_set(click);
                                    cout << "Intersection: " << endl;
                                    show_set(getIntersection(click, nodeInfo.node_set));
                                }
                            }
                        }
                    }
                    cout << "Neighbors: " << endl;
                    for (auto n2: pair_neigh.second) {
                        cout << "Neighbor: " << _g[n2.first].node.str() << endl;
                        cout << "PairedInfo: " << endl;
                        show_set(_g[n2.first].node_set);
                        cout << "Cliques with parent: " << endl;
                        for (auto cliques: _g[n2.first].parent_cliques) {
                            cout << "Parent: " << cliques.first.str() << endl;
                            if (cliques.first == nodeInfo.node.str()) {
                                for (auto click:cliques.second) {
                                    cout << "Clique: " << endl;
                                    show_set(click);
                                    cout << "Intersection: " << endl;
                                    show_set(getIntersection(click, nodeInfo.node_set));
                                }
                            }
                        }
                        cout << "Intersection: " << endl;
                        show_set(getIntersection(n2.second, nodeInfo.node_set));
                    }
                    cin.get();
                }
                    for (auto n2: pair_neigh.second) {
                        if (node_launched.find(it_node_set(n2)) == node_launched.end()) {
                            node_launched.emplace(it_node_set(n2));
                            DnaSequence seq("");
                            seq.append_nuc_right(nodeInfo.node.at(0));
                            stack_prev_seqs.emplace(seg);
                            stack_seq.emplace(seq);
                            stack_node.emplace(n2);
                        }
                    }
            }
            if (neighbors.size() > 1) {
                if (show) {
                    cout << "NEIGHBORS HIGHER " << endl;
                    cout << "Kmer: " << nodeInfo.node.str() << "\n";
                    cout << "Haplotype: " << endl;
                    show_set(node_assay.second);
                    cout << "Intersection: " << endl;
                    show_set(getIntersection(node_assay.second, nodeInfo.node_set));
                    cout << "Paired_info: " << endl;
                    show_set(nodeInfo.node_set);
                    cout << "Neighbors Reales:" << endl;
                    for (auto n2: neighbors) {
                        cout << "Neighbor: " << _g[n2.first].node.str() << endl;
                        cout << "Cliques with parent: " << endl;
                        for (auto cliques: _g[n2.first].parent_cliques) {
                            cout << "Parent: " << cliques.first.str() << endl;
                            if (cliques.first == nodeInfo.node.str()) {
                                for (auto click:cliques.second) {
                                    cout << "Clique: " << endl;
                                    show_set(click);
                                    cout << "Intersection: " << endl;
                                    show_set(getIntersection(click, nodeInfo.node_set));
                                }
                            }
                        }
                    }
                    cout << "Neighbors: " << endl;
                    for (auto n2: pair_neigh.second) {
                        cout << "Neighbor: " << _g[n2.first].node.str() << endl;
                        cout << "Cliques with parent: " << endl;
                        for (auto cliques: _g[n2.first].parent_cliques) {
                            cout << "Parent: " << cliques.first.str() << endl;
                            if (cliques.first == nodeInfo.node.str()) {
                                for (auto click:cliques.second) {
                                    cout << "Clique: " << endl;
                                    show_set(click);
                                    cout << "Intersection: " << endl;
                                    show_set(getIntersection(click, nodeInfo.node_set));
                                }
                            }
                        }
                        cout << "Intersection: " << endl;
                        show_set(getIntersection(n2.second, nodeInfo.node_set));
                    }
                    cin.get();
                }
                //Checking
                bool same_neigh = true;
                Node pivot = _g[neighbors[0].first].node;
                for (auto it:neighbors) {
                    if (pivot != _g[it.first].node) {
                        same_neigh = false;
                        break;
                    }
                }
                if (same_neigh) {
                    masSame++;
                    //Launching new sequences:
                    if (node_launched.find(it_node_set(node_assay)) == node_launched.end())
                    {
                        node_launched.emplace(it_node_set(node_assay));
                        for (auto n2:neighbors)
                        {
                            DnaSequence seq(sequence);
                            seq.append_nuc_right(nodeInfo.node.at(0));
                            stack_seq.emplace(seq);
                            stack_node.emplace(n2);
                        }
                        for (auto n2: pair_neigh.second)
                        {
                            if (node_launched.find(it_node_set(n2)) == node_launched.end()) {
                                node_launched.emplace(it_node_set(n2));
                                DnaSequence seq("");
                                seq.append_nuc_right(nodeInfo.node.at(0));
                                stack_prev_seqs.emplace(seg);
                                stack_seq.emplace(seq);
                                stack_node.emplace(n2);
                            }
                        }
                    }
                    continue;
                }
            }
            if (neighbors.size() > 1) {
                mas++;
                if (node_launched.find(it_node_set(node_assay)) == node_launched.end()) {
                    node_launched.emplace(it_node_set(node_assay));
                    for (auto n2: neighbors) {
                        DnaSequence seq("");
                        seq.append_nuc_right(nodeInfo.node.at(0));
                        stack_prev_seqs.emplace(seg);
                        stack_seq.emplace(seq);
                        stack_node.emplace(n2);
                    }
                    //Launching new sequences:
                    for (auto n2: pair_neigh.second) {
                        if (node_launched.find(it_node_set(n2)) == node_launched.end()) {
                            node_launched.emplace(it_node_set(n2));
                            DnaSequence seq("");
                            seq.append_nuc_right(nodeInfo.node.at(0));
                            stack_prev_seqs.emplace(seg);
                            stack_seq.emplace(seq);
                            stack_node.emplace(n2);
                        }
                    }
                }
            }
            for (uint i = 0; i < kmer_size - 1; ++i)
                sequence.append_nuc_right(nodeInfo.node.at(i));
            map_seq_nodo_end[seg] = node_assay.first;
            map_seqs[seg] = sequence;
            size_t seqTmp = stack_prev_seqs.top();stack_prev_seqs.pop();
            graphUnitigs[seqTmp].push_back(seg);
            graphUnitigs[seg] = vector<size_t>();
            seg++;
        }
    }
}

template<>
void boostDBG<true>::extension(vector <Node> in_0, string path_to_write) {
    cout << "--------------------------------\nBoostDBG extension\n--------------------------------\n";
    map <size_t, vector<size_t>> map_nodo_seq_start;
    map <size_t, graphBU> map_seq_nodo_end;
    map <size_t, DnaSequence> map_seqs;
    unordered_set<graphBU> nodes_extended;
    size_t seqId = 0;
    for (auto n:_in_0_pairs)
    {
        if (node_launched.find(n) == node_launched.end()) {
            DnaSequence sequence("");
            sequence.append_nuc_right(_g[n.first].node.at(0));
            map_nodo_seq_start[n.first].push_back(seg);
            _transverse(n, map_nodo_seq_start, map_seq_nodo_end, map_seqs, nodes_extended, sequence, seqId);
        }
    }
    cout << "PercentageTransversed: "<<endl;
    cout << (float)(nodes_extended.size() /  boost::num_vertices(_g))<<"%"<< endl;
    cout << "Writing unitigs\n";
    _writeUnitigs(map_nodo_seq_start,map_seq_nodo_end, map_seqs, path_to_write);
}

template<>
void boostDBG<true>::_remove_outliers(Extra<true> * paired, DBG<true> * dbg, unordered_set<Node> * node_set)
{
    vertex_iterator v, vend;
    size_t num_vertex = boost::num_vertices(_g), num_representants = _first_last.size();
    reach = vector<bool>(num_representants*num_representants, false);
    size_t distance, executions;
    unordered_set<graphBU> check_store;
    //Polishing
    cout << "Inner Stage: Filling out reach matrix for polishing"<<endl;
    Progress::get().size_total = num_representants;
    bool option = false;
    size_t progress = 0;
    for (auto r:_first_last)
    {
        distance = 0;executions = 0;
        bool token = false;
        check_store.emplace(r);
        _full_fil_matrix(distance, num_representants, executions, r, r, 0, reach, &token, check_store);
        check_store.clear();
        Progress::update(progress++);
    }
    Progress::update(num_representants);
    //Checking outliers: not naive way
    Progress::get().size_total = num_vertex;
    size_t cont = 0;
    if (Parameters::get().polish) {
        cout << "Inner Stage: Polishing paired-end information" << endl;
        //#pragma omp parallel for
        for (size_t i = 0; i < num_vertex; ++i) {
            Progress::update(cont++);
            pair<bool, FreqMap> pair_end_info;
            if (option)
                pair_end_info = dbg->getExtra(_g[i].node);
            else
                _extra_info.getInfoNode(node_set->find(_g[i].node));
            if (pair_end_info.first) {
                FreqMap tmp;
                auto set_pair_end = pair_end_info.second;
                for (auto p: pair_end_info.second) {
                    if (p.second > 4 * POLISH_PATH_LEN) {
                        tmp[p.first] = 4 * POLISH_PATH_LEN;
                        set_pair_end.erase(p.first);
                    }
                }
                auto set_pair_end_2 = set_pair_end;
                for (auto p: set_pair_end) {
                    vertex_t node1 = local_map[(*p.first)];
                    //(tmp.find(p.first) == tmp.end())?tmp[p.first] = p.second:tmp[p.first]+=p.second;
                    if (tmp.find(p.first) == tmp.end())
                        tmp[p.first] = 0;
                    for (auto p2: set_pair_end_2) {
                        if (p.first == p2.first)
                            continue;
                        vertex_t node2 = local_map[(*p2.first)];
                        bool reachable = (reach[translator_vector[_g[_represent_first_last[node1].second].id] * num_representants
                                                + translator_vector[_g[_represent_first_last[node2].first].id]] ||
                                          reach[translator_vector[_g[_represent_first_last[node2].second].id] * num_representants
                                                + translator_vector[_g[_represent_first_last[node1].first].id]]);
                        if (translator_vector[_g[node1].id] == INF || translator_vector[_g[node2].id] == INF)
                            cout << "FAIL mayusculo"<<endl;
                        if (reachable) {
                            tmp[p.first]++;
                            (tmp.find(p2.first) == tmp.end()) ? tmp[p2.first] = 1 : tmp[p2.first]++;
                        }
                    }
                    set_pair_end_2.erase(p.first);
                }
                paired->mod_freqs(_g[i].node, tmp);
            }
        }
        Progress::update(num_vertex);
        cout << "Plotting Histogram\n";
        paired->show_distribution("_post_processing");
    }
    cout << "Removing outliers" << endl;
    paired->polish();
    cont = 0;
    Progress::get().size_total = num_vertex;
    cout << "Adding paired-end info"<<endl;
    //#pragma omp parallel for
    for (size_t i = 0; i < num_vertex; ++i)
    {
        Progress::update(cont++);
        pair<bool, FreqMap> pair_end_info;
        if (option)
            pair_end_info = dbg->getExtra(_g[i].node);
        else
            pair_end_info = _extra_info.getInfoNode(node_set->find(_g[i].node));
        if (pair_end_info.first)
        {
            for (auto p: paired->getPairs(_g[i].node)) {
                _g[i].node_set.emplace((*p));
            }
        }
    }
    Progress::update(num_vertex);
}

template<>
void boostDBG<true>::_get_representatives()
{
    cout << "Getting representatives\n";
    size_t num_vertex = boost::num_vertices(_g);
    vector<vertex_t> starts;
    vertex_iterator v, vend;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        if ((boost::in_degree(*v,_g) > 1) || (boost::out_degree(*v,_g) > 1) || (boost::in_degree(*v,_g) == 0))
            starts.push_back(*v);
    }
    cout << "Starts Size: "<<starts.size()<<endl;
    vector<bool> visited(num_vertex, false);
    for (auto p: starts){
        stack<vertex_t> kmer_check;
        vector<vertex_t> unitig_tmp;
        kmer_check.push(p);
        while(!kmer_check.empty())
        {
            vertex_t curr_kmer = kmer_check.top();
            size_t in_degree = boost::in_degree(curr_kmer,_g), out_degree=boost::out_degree(curr_kmer,_g);
            kmer_check.pop();
            bool prev_state = visited[_g[curr_kmer].id];
            if (!visited[_g[curr_kmer].id])
            {
                for (auto neigh:getKmerNeighbors(curr_kmer)) {
                    kmer_check.push(neigh);
                }
                visited[_g[curr_kmer].id] = true;
            }
            if ((in_degree > 1) || (out_degree > 1) || (out_degree == 0))
            {
                if (unitig_tmp.size())
                {
                    if (((out_degree > 1) || (out_degree == 0)) && (in_degree <= 1))
                        unitig_tmp.push_back(curr_kmer);
                    graphBU representative = unitig_tmp[floor((float) unitig_tmp.size() / 2)];
                    _representants.push_back(representative);
                    _represent_first_last[representative] = pair<graphBU, graphBU>(unitig_tmp[0],
                                                                                   unitig_tmp[unitig_tmp.size() - 1]);
                    _representants_hits[representative] = unordered_map<graphBU, size_t>();
                    for (auto k:unitig_tmp)
                        _representants_map[k] = representative;
                }
                unitig_tmp.clear();
                if ((in_degree > 1) && !prev_state && (out_degree <= 1))
                    unitig_tmp.push_back(curr_kmer);
                if ((in_degree > 1) && (out_degree > 1) && !prev_state){
                    _representants.push_back(curr_kmer);
                    _represent_first_last[curr_kmer] = pair<graphBU, graphBU>(curr_kmer, curr_kmer);
                    _representants_map[curr_kmer] = curr_kmer;
                }
            }else
                unitig_tmp.push_back(curr_kmer);
        }
    }
    /*
     * Fill vector
     */
    size_t num_representant = 0, num_representants_real = 0;
    for (auto r:_represent_first_last) {
        if (_first_last.find(r.second.first) == _first_last.end()) {
            translator_vector[_g[r.second.first].id] = num_representant++;
            _first_last.emplace(r.second.first);
        }
        if (_first_last.find(r.second.second) == _first_last.end()) {
            translator_vector[_g[r.second.second].id] = num_representant++;
            _first_last.emplace(r.second.second);
        }
    }
    num_representants_real = num_representant;
    num_representant = INF+1;
    for (auto r:_represent_first_last)
        if (_first_last.find(r.first) == _first_last.end())
            translator_vector[_g[r.first].id] = num_representant++;
    /*
     * Remove
     */
    size_t nodes = 0;
    for (auto b:visited) {
        if (!b) {
            cout << "Nodes still not visited: "<<_g[nodes].node.str() << endl;
            cout << "Indegree: "<<boost::in_degree(nodes, _g)<<" Outdegree: "<<boost::out_degree(nodes, _g)<<endl;
            for (auto k: getKmerNeighbors(nodes)){
                cout << "Neighbor: "<<_g[k].node.str()<<endl;
            }
            //exit(1);
        }
        nodes++;
    }
    cout << "End: "<<_representants.size()<<" representantes. Numero de representantes: "<<num_representants_real<<endl;
    //cin.get();
    //exit(1);
}

template<>
boostDBG<true>::boostDBG(string path_to_file, string dir_path, SequenceContainer * sc)
{
    cout << "STAGE: Building boost graph from file"<<endl;
    size_t retries = 0;
    _thirdPartyKmerCounting(path_to_file, dir_path, &retries);
    cout << "Initial information: "<<endl;
    _printInfo();
    //_show_internal_info();
    _polishing();
}
template<>
boostDBG<true>::boostDBG(DBG<true> * dbg)
{
    std::cout << "STAGE: Building boost graph (standard)"<<endl;
    pair<unordered_set<Node> * , unordered_set<Node>* > graph_struct = dbg->getNodes();
    for (auto k: (*graph_struct.second))
    {
        vector <Node> neigh = dbg->getKmerNeighbors(k);
        if (local_map.find(k) == local_map.end())
        {
            local_map[k] = boost::add_vertex(NodeInfo(k, _node_id++),_g);
        }
        for (auto k2: neigh) {
            if (graph_struct.second->find(k2)!=graph_struct.second->end())
            {
                if (local_map.find(k2) == local_map.end())
                {
                    local_map[k2] = boost::add_vertex(NodeInfo(k2, _node_id++),_g);
                    _g[local_map[k2]].parent_cliques[k] = vector<PairedInfoNode>();
                }else {
                    _g[local_map[k2]].parent_cliques[k] = vector<PairedInfoNode>();
                }
                boost::add_edge(local_map[k], local_map[k2], _g);
            }
        }
    }
    _node_reads = dbg->getNodeReads();
    //Fill_vector
    translator_vector = vector<size_t>(boost::num_vertices(_g), INF);
    cout << "STAGE: Getting representatives\n";
    _get_representatives();
    cout << "STAGE: Getting and polishing paired info\n";
    SequenceContainer * sc = dbg->getSequenceContainer();
    _insert_extra(sc, graph_struct.first, graph_struct.second);
    //Extra<true> * paired_info_pointer = dbg->getPairedInfo();
    Extra<true> * paired_info_pointer = &(_extra_info);
    _remove_outliers(paired_info_pointer, dbg,graph_struct.second);
    cout << "STAGE: Adding paired-end information\n";
    _insertExtraInfo(local_map, dbg);
    _modify_info();
    //show_info(20);
    vertex_iterator v, vend;
    cout << "STAGE: Getting starting points"<<endl;
    for (boost::tie(v,vend) = boost::vertices(_g);v!=vend;++v)
    {
        if (in_degree(*v) == 0) {
            //_in_0.push_back(*v);
            for (auto neigh: getKmerNeighbors(*v))
            {
                for (auto cliques: _g[neigh].parent_cliques)
                {
                    if (cliques.first == _g[*v].node) {
                        for (auto clique:cliques.second)
                            _in_0_pairs.push_back(pair<vertex_t, PairedInfoNode>(*v, clique));
                    }
                }
            }
        } else {
            for (auto parent: _g[*v].parent_cliques) {
                for (auto clique: parent.second) {
                    size_t in_nodes = in_degree(pair<vertex_t, PairedInfoNode>(*v, clique), false);
                    if (in_nodes == 0) {
                        //_in_0.push_back(*v);
                        _in_0_pairs.push_back(pair<vertex_t, PairedInfoNode>(*v, clique));
                        //break;
                    }
                }
            }
        }
    }
    cout << "Suspicious Starts: "<<_in_0.size()<<endl;
    cout << "Suspicious Real Starts: "<<_in_0_pairs.size()<<endl;
    //cin.get();
}

//ListDBG
template<>
void listDBG<false>::_transverse(Node n,
                                 map <Node, vector<size_t>> & map_nodo_seq_start,
                                 map <size_t, Node> & map_seq_nodo_end,
                                 map <size_t, DnaSequence> & map_seqs,
                                 DnaSequence & sequence)
{
    size_t kmer_size = Parameters::get().kmerSize;
    vector<Node> neighbors = getKmerNeighbors(n);
    if ((neighbors.size()==1) && (in_degree(n)==1))
    {
        sequence.append_nuc_right(n.at(0));
        _transverse(neighbors[0], map_nodo_seq_start, map_seq_nodo_end,map_seqs,sequence);
    }else
    {
        map_seq_nodo_end[seg] = n;
        if (neighbors.size()==0)
        {
            for (uint i = 0; i < kmer_size-1; ++i)
                sequence.append_nuc_right(n.at(i));
        }
        map_seqs[seg++] = sequence;
        if (neighbors.size() > 1)
        {
            if (map_nodo_seq_start.find(n) == map_nodo_seq_start.end())
            {
                map_nodo_seq_start[n] = vector<size_t>();
                for (auto n2: neighbors)
                {
                    DnaSequence seq("");
                    seq.append_nuc_right(n.at(0));
                    map_nodo_seq_start[n].push_back(seg);
                    _transverse(n2,map_nodo_seq_start, map_seq_nodo_end,map_seqs,seq);
                }
            }else{
                map_seq_nodo_end[seg-1] = n;
            }
        }
        if (in_degree(n) != 1)
        {
            if (map_nodo_seq_start.find(n) != map_nodo_seq_start.end()){
                map_seq_nodo_end[seg-1] = n;
            }
            else
            {
                for (auto n2: neighbors)
                {
                    DnaSequence seq("");
                    seq.append_nuc_right(n.at(0));
                    map_nodo_seq_start[n] = vector<size_t>(1,seg);
                    _transverse(n2, map_nodo_seq_start, map_seq_nodo_end, map_seqs,seq);
                }
            }
        }
    }
}

template<>
void listDBG<false>::_writeUnitigs(map <Node, vector<size_t>> map_nodo_seq_start,
                                   map <size_t, Node> map_seq_nodo_end,
                                   map <size_t, DnaSequence> map_seqs,
                                   string file_to_path)
{
    unordered_set<DnaSequence> seqs;
    vector<pair<size_t,size_t>> links;
    vector<bool> strand;
    for (auto s:map_seqs)
        seqs.emplace(s.second);
    for (auto s:map_nodo_seq_start)
        for (auto origin: s.second)
            for (auto end: map_nodo_seq_start[map_seq_nodo_end[origin]])
                links.push_back(pair<size_t,size_t>(end, origin));
    UnitigExtender<false>::_write_gfa(file_to_path, map_seqs,links);
}

template<>
void listDBG<false>::extension(vector <Node> in_0, string path_to_write)
{
    map<Node, vector<size_t>> map_nodo_seq_start;
    map<size_t, Node> map_seq_nodo_end;
    map<size_t, DnaSequence> map_seqs;
    for (auto n: _in_0)
    {
        map_nodo_seq_start[n] = vector<size_t>();
        for (auto n2: getKmerNeighbors(n))
        {
            DnaSequence sequence("");
            sequence.append_nuc_right(n.at(0));
            map_nodo_seq_start[n].push_back(seg);
            _transverse(n2, map_nodo_seq_start,map_seq_nodo_end,map_seqs,sequence);
        }
    }
    _writeUnitigs(map_nodo_seq_start,map_seq_nodo_end,map_seqs, path_to_write);

}
template<>
void listDBG<false>::_buildNewGraph(DBG<false> * dbg)
{
    cout << "STAGE: IntoAdjacentLists\n";
    /*
     * Check correction
     */
    pair<unordered_set<Node> *,unordered_set<Node> *> graph_shape = dbg->getNodes();
    _solid_kmers = (*graph_shape.first);
    for (auto n: (*graph_shape.second))
    {
        vector<Node> neighbors = dbg->getKmerNeighbors(n);
        _g[n].second = neighbors;
        for (auto n_2: neighbors)
        {
            if (_g.find(n_2)!=_g.end())
                _g[n_2].first.push_back(n);
            else
                _g[n_2].first = vector<Node>(1,n);
        }
    }
    //GetSuspicious Starts as unbalanced nodes
    for (auto n: _g)
    {
        if (n.second.first.size() < n.second.second.size())
        {
            _in_0.push_back(n.first);
        }
    }
    cout << "End Translation\n";
}

template<>
listDBG<false>::listDBG(DBG<false> * dbg)
{
    _buildNewGraph(dbg);
}
