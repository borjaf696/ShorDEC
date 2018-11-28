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
        kmer_aux = Kmer(kmer.str());
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
        /*
         * Particularidad mia de como he implementado el sequence container
         */
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
/*
 * Process ExtraInfo
 */
template<>
pair<bool, DBG<true>::Parent_Paired_Info> NaiveDBG<true>::getExtra(Node node)
{
    PairedInfoNode paired_info;
    unordered_set<Node>::iterator it = _dbg_nodes.find(node);
    cout << "Puntero: "<<&(it) << "; Nodo: "<<node.str()<<"\n";
    pair<bool, ExtraInfoNode> full_info = _extra_info.getInfoNode((it));
    for (auto e: full_info.second)
        paired_info.emplace(*e);
    return pair<bool, PairedInfoNode>(full_info.first,paired_info);
}

template<>
void NaiveDBG<true>::_insert_extra_info(unordered_map<Node, unordered_set<Node>> & paired_info)
{
    cout << "STAGE: Adding extra info\n";
    Parameters::get().kmerSize = Parameters::get().kmerSize - 1;
    for (auto &read:_sc.getIndex())
    {
        Progress::update(read.first.getId());
        bool check = (_is_standard)?_is_standard:(read.first.getId() % 2 == 0);
        if (check && ((read.first.getId() % 4) < 2))
        {
            for (auto k: IterKmers<true>(_sc.getSeq(read.second.getId()), _sc.getSeq(read.second.getPairId())))
            {
                FuncNode nonstd_pk = k.pair_kmer;
                pair<Node, Node> sep_nodes = nonstd_pk.getKmers();
                //auto start_build = std::chrono::high_resolution_clock::now();
                unordered_set<Node>::const_iterator node_f_it = _dbg_nodes.find(sep_nodes.first),
                        node_s_it = _dbg_nodes.find(sep_nodes.second);
                if (( node_f_it != _dbg_nodes.end()) &&
                        (node_s_it != _dbg_nodes.end()))
                {
                    _extra_info.insert(node_f_it, node_s_it);
                    cout << "Node from reads: "<<sep_nodes.first.str()<<"; Puntero: \n";
                    cout << "Node from iterator: "<<(*node_f_it).str()<<"\n";
                    cout << "Pair from reads: "<<sep_nodes.second.str()<<"\n";
                    cout << "Pair from iterator: "<<(*node_s_it).str()<<"\n";
                }
                /*if (is_solid(sep_nodes.first) && is_solid(sep_nodes.second))
                {
                    unordered_map<Node, unordered_set<Node>>::iterator key_it = paired_info.find(sep_nodes.first);
                    if ((key_it != paired_info.end()) && !(paired_info[sep_nodes.first].find(sep_nodes.second)!=paired_info[sep_nodes.first].end()))
                        paired_info[sep_nodes.first].emplace(sep_nodes.second);
                    else if (key_it == paired_info.end())
                    {
                        paired_info[sep_nodes.first] = unordered_set<Node>({sep_nodes.second});
                    }
                }*/
                /*auto finish_build = std::chrono::high_resolution_clock::now();
                cout << "Time out: "<<(finish_build - start_build).count()<<"\n";*/
            }
        }
    }
    Parameters::get().kmerSize = Parameters::get().kmerSize + 1;
    Progress::update(_sc.size());
}
template<>
void NaiveDBG<true>::_set_extra_info(unordered_map<Node, unordered_set<Node>> paired_info)
{
    cout << "STAGE: Setting extra info\n";
    for (auto duo: paired_info)
    {
        for (auto pair:duo.second)
        {
            /*vector<DnaSequence> all_seqs_left = duo.first.firstLastSubstr(Parameters::get().kmerSize-1, _is_standard),
                    all_seqs_right = pair.firstLastSubstr(Parameters::get().kmerSize-1, _is_standard);
            _extra_info.insert(all_seqs_left[0], all_seqs_right[0]);
            _extra_info.insert(all_seqs_left[1], all_seqs_right[1]);*/
            /*auto finish = std::chrono::high_resolution_clock::now();
            cout << "Time within: "<<(finish-start).count()<<"\n";*/
            /*if (_is_standard)
            {
                _extra_info.insert(all_seqs_right[2], all_seqs_left[2]);
                _extra_info.insert(all_seqs_right[3], all_seqs_left[3]);
            }*/
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
            _insert(kmer.first, kmer.first);
        }
    }
    std::cout<<"Total Solid K-mers(Graph Edges): "<<_dbg_naive.size()
             <<" Total Graph Nodes: "<<_dbg_nodes.size()<<"\n";
    _kmers_map.clear();
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

}

/*
 * Boost graphs -> Pair-end Reads
 */
template<>
int* boostDBG<true>::_floyds_warshall()
{
    std::cout << "Lets compute floyds warshall\n";
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
        vector<size_t> neigh;
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
                                vector<bool> & check_vect)
{
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
    size_t num_vertex = boost::num_vertices(_g);
    bool out = false;
    for (auto n: neigh)
    {
        reach_vect[_g[oriSource].id*num_vertex+_g[n].id] = true;
        check_vect[_g[oriSource].id*num_vertex+_g[n].id] = true;
        check_vect[_g[n].id*num_vertex+_g[oriSource].id] = true;
        out |= _reachable(oriSource, n, target, distance++, branches, reach_vect, check_vect);
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
    typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, NodeInfo> Graph_l;
    typedef Graph_l::vertex_descriptor vertex_graph;
    typedef Graph_l::out_edge_iterator out_it_local;
    typedef Graph_l::vertex_iterator vertex_it;
    /*
     * We are going to use pair_info (ExtraInfoNodes) to "modify" the graph including this new information
     */
    cout << "STAGE: Processing cliques\n";
    vertex_iterator v, vend;
    size_t num_vertex = boost::num_vertices(_g);
    int * distance_matrix;
    if (FLOYD)
        distance_matrix = _floyds_warshall();
    //TODO: Smaller
    vector<bool> reach_vect(num_vertex*num_vertex, false);
    vector<bool> check_vect(num_vertex*num_vertex, false);
    //New_graph container
    Graph tmp_graph;
    //Progression
    Progress::get().size_total = num_vertex;
    //Incremental information:
    unordered_set<graphBU> last_set = unordered_set<graphBU>();
    Graph_l local_graph;
    map<Node, vertex_graph> local_node_map;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges((*v), _g);
        NodeInfo node_info = _g[*v];
        Progress::update(node_info.id);
        if (_map_extra_info[node_info.id].empty())
            continue;
        for (; neighbors.first != neighbors.second; ++neighbors.first)
        {
            /*
             * Local Graph which joins evey single reachable node (from the others)
             * Reminder:
             *      - There are several nodes which are not reachable from me but can belong to my cliques
             *      - Ocurre cuando toda la informacion contenida en un click por parte del hijo ya ha sido introducida
             *      en otro click con anterioridad. En este caso probablemente estemos hablando de un puente. (punto
             *      de conexion entre multiples bifurcaciones)
             */
            auto endpoint = boost::target(*neighbors.first,_g);
            NodeInfo neigh_info = _g[endpoint];
            vector<bool> visited(num_vertex, false);
            size_t curr_node = 0;
            if (_map_extra_info[neigh_info.id].empty())
            {
                continue;
            }
            unordered_set<graphBU> extra_parent = _map_extra_info[node_info.id],
                    extra_neigh = _map_extra_info[neigh_info.id];
            unordered_set<vertex_t> local_vect = getUnion(extra_parent, extra_neigh);
            //Change last graph:
            unordered_set<graphBU> s_set = sustract(last_set,local_vect),
                    t_set = sustract(local_vect, last_set);
            for (auto s: s_set)
            {
                NodeInfo s_info = _g[s];
                if (local_node_map.find(s_info.node) != local_node_map.end())
                {
                    boost::clear_vertex(local_node_map[s_info.node], local_graph);
                    boost::remove_vertex(local_node_map[s_info.node], local_graph);
                    local_node_map.erase(s_info.node);
                }
            }
            last_set = local_vect;
            //cout << "LocalVect: "<<local_vect.size()<<"\n";
            auto start_build = std::chrono::high_resolution_clock::now();
            for (auto s:last_set)
            {
                NodeInfo s_info = _g[s];
                vertex_graph source;
                if (local_node_map.find(s_info.node) == local_node_map.end())
                {
                    source = boost::add_vertex(NodeInfo(s_info.node,curr_node++), local_graph);
                    local_node_map[s_info.node] = source;
                }else {
                    source = local_node_map[s_info.node];
                }
                for (auto t:local_vect)
                {
                    NodeInfo t_info = _g[t];
                    bool reached;
                    if (check_vect[s_info.id*num_vertex+t_info.id])
                        reached = reach_vect[s_info.id*num_vertex+t_info.id]||reach_vect[t_info.id*num_vertex+s_info.id];
                    else
                    {
                        size_t  branches = MAX_BRANCHES_CHECK;
                        reached = ((FLOYD)?
                                   _reachable(distance_matrix, s_info.id, t_info.id):
                                   _reachable(s,s,t,0, &branches, reach_vect, check_vect));
                        branches = MAX_BRANCHES_CHECK;
                        reached = ((reached)?
                                   reached:
                                   ((FLOYD)?_reachable(distance_matrix, t_info.id, s_info.id):
                                                      _reachable(t,t,s,0, &branches, reach_vect, check_vect)));
                    }
                    if (reached)
                    {
                        vertex_graph target;
                        if (local_node_map.find(t_info.node) == local_node_map.end())
                        {
                            target = boost::add_vertex(NodeInfo(t_info.node, curr_node++), local_graph);
                            local_node_map[t_info.node] = target;
                        }else
                            target = local_node_map[t_info.node];
                        boost::add_edge(source, target, local_graph);
                    }
                }
                local_vect.erase(s);
            }
            auto finish_build = std::chrono::high_resolution_clock::now();
            //cout << "TIME: "<<(finish_build-start_build).count()<<"\n";
            /*
             * Polish Graph
             */
            map<vertex_graph, vector<Node>> clique_representation;
            unordered_set<vertex_graph> nodes_erase;
            auto _mod_clique_graph = [&clique_representation, &local_graph, &nodes_erase] ()
            {
                vertex_it v, vend;
                for (boost::tie(v, vend) = boost::vertices(local_graph); v != vend; ++v)
                {
                    if (nodes_erase.find(*v) == nodes_erase.end())
                    {
                        clique_representation[(*v)] = vector<Node>();
                        pair<out_it_local, out_it_local> neigh = boost::out_edges((*v), local_graph);
                        for (; neigh.first != neigh.second; ++neigh.first)
                        {
                            auto endpoint = boost::target(*neigh.first,local_graph);
                            if (endpoint != (*v))
                            {
                                if (boost::degree(endpoint, local_graph) == boost::degree((*v), local_graph))
                                {
                                    nodes_erase.emplace(endpoint);
                                    clique_representation[(*v)].push_back(local_graph[endpoint].node);
                                }
                            }
                        }
                        if (clique_representation[(*v)].empty())
                            clique_representation.erase((*v));
                    }
                }
            };
            /*
             * Calculate all cliques in the local graph
             */
            auto start = std::chrono::high_resolution_clock::now();
            std::set<size_t> idCliques;
            priority_queue<pair<size_t,vector<vertex_graph>>> output;
            size_t num_vertex_local = boost::num_vertices(local_graph),
                    num_edges_local = boost::num_edges(local_graph);
            if (((num_vertex_local*(num_vertex_local-1))/2) == (num_edges_local-num_vertex_local))
            {
                PairedInfoNode uniqueHaplotype;
                for (auto s:last_set)
                    uniqueHaplotype.emplace(_g[s].node);
                _g[endpoint].parent_cliques[node_info.node].push_back(uniqueHaplotype);
            }else {
                if (num_vertex_local > 2)
                    _mod_clique_graph();
                for (auto n: nodes_erase) {
                    if (clique_representation.find(n) != clique_representation.end()) {
                        clique_representation[n].clear();
                        clique_representation.erase(n);
                    }
                    local_node_map.erase(local_graph[n].node);
                    boost::clear_vertex(n, local_graph);
                    boost::remove_vertex(n, local_graph);
                }
                output = findMaxClique < Graph_l, vertex_graph, vertex_it > (local_graph, idCliques);

                auto finish = std::chrono::high_resolution_clock::now();
                unordered_set <vertex_graph> edge_transversed;
                auto min = neigh_info.node_set;
                size_t num_clicks = 0;
                /*cout << "Local Graph edges: " << boost::num_edges(local_graph) << " Local Graph nodes: "
                     << boost::num_vertices(local_graph) << "\n";*/
                vertex_it v2, vend2;
                /*for (boost::tie(v2, vend2) = boost::vertices(local_graph); v2 != vend2; ++v2) {
                    cout << "Vertex: "<<local_graph[*v2].node.str()<<"\n";
                    pair <out_it_local, out_it_local> neigh2 = boost::out_edges((*v2), local_graph);
                    for (; neigh2.first != neigh2.second; ++neigh2.first) {
                        auto endpoint = boost::target(*neigh2.first,local_graph);
                        cout << "Neigh: "<<local_graph[endpoint].node.str()<<"\n";
                    }
                }*/
                if (boost::num_edges(local_graph) == 0)
                    _g[endpoint].parent_cliques[node_info.node].push_back(PairedInfoNode());
                while (!output.empty() & (edge_transversed.size() != boost::num_edges(local_graph))) {
                    pair <size_t, vector<vertex_graph>> top_click = output.top();
                    output.pop();
                    vector <vertex_graph> clique = top_click.second;
                    PairedInfoNode local_haplotype;
                    for (size_t i = 0; i < clique.size(); ++i) {
                        local_haplotype.emplace(local_graph[clique[i]].node);
                        edge_transversed.emplace(clique[i]);
                    }
                    /*
                     * Add cliques (presumed haplotypes which contains both nodes)
                     */
                    for (auto p:clique_representation) {
                        Node key = local_graph[p.first].node;
                        if (local_haplotype.find(key) != local_haplotype.end()) {
                            for (auto n:p.second)
                                local_haplotype.emplace(n);
                        }
                    }
                    if (!getIntersection(local_haplotype, node_info.node_set).empty() &
                        !getIntersection(local_haplotype, neigh_info.node_set).empty()) {
                        _g[endpoint].parent_cliques[node_info.node].push_back(local_haplotype);
                        min = getIntersection(min, local_haplotype);
                        num_clicks++;
                    }
                }
                /*
                 * Clean cliques
                 *      * min = max min set from neighbor -> but we need to have some neigh info.
                 */
                vector <vector<PairedInfoNode>::iterator> erase;
                for (auto i = _g[endpoint].parent_cliques[node_info.node].begin();
                     i != _g[endpoint].parent_cliques[node_info.node].end(); ++i) {
                    auto content = getIntersection(neigh_info.node_set, *i);
                    if (isSame(content, min)) {
                        if (((*i).size() != min.size()) & (num_clicks > 1))
                            erase.push_back(i);
                    }
                }
                /*
                 * Cleaning:
                 */
                for (auto s: erase)
                    _g[endpoint].parent_cliques[node_info.node].erase(s);
            }
        }
    }
    //Free Floyd
    if (FLOYD)
        free(distance_matrix);
    Progress::update(num_vertex);
    std::cout << "\n";
}
template<>
void boostDBG<true>::show_info()
{
    vertex_iterator v, vend;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        cout << " Kmer:"     << _g[*v].node.str()
                  << " id:"  << _g[*v].id
                  << " Puntero: " << (*v)
                  << "\n";
        vector<Node> neighbors = getKmerNeighbors(_g[*v].node);
        cout << "Neighbors: ";
        for (auto n:neighbors)
            std::cout  << n.str()<<" ";
        cout << "\n Couples: ";
        for (auto n:_g[*v].node_set)
            std::cout << n.str()<<" ";
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
void boostDBG<true>::_transverse(it_node n,
                                 map<graphBU, vector<size_t>> & map_nodo_seq_start,
                                 map<size_t, graphBU> & map_seq_nodo_end,
                                 map<size_t, DnaSequence> & map_seqs,
                                 DnaSequence & sequence)
{
    NodeInfo nodeInfo = _g[n.first];
    size_t kmer_size = Parameters::get().kmerSize, in_nodes = in_degree(n);
    cout << "Kmer: "<<nodeInfo.node.str()<<"\n HaplotipoACC: \n";
    for( auto t: n.second)
        cout << t.str()<<"\n";
    cout << "Kmer Extra Info: \n";
    for (auto t: _g[n.first].node_set)
        cout << t.str()<<"\n";
    cout << "Sequence: "<<sequence.str()<<"\n";
    vector<it_node> neighbors = getOutKmerNeighbors(n);
    cout << "NumNeighbors: "<<neighbors.size()<<":"<<in_nodes<<"\n";
    for (auto it: neighbors)
    {
        cout <<"VECINO: "<< _g[it.first].node.str() <<"\n";
        for (auto hap1: it.second)
            cout <<hap1.str()<<"\n";
        cout << "\n";
    }
    if ((neighbors.size()==1) && (in_nodes==1))
    {
        sequence.append_nuc_right(nodeInfo.node.at(0));
        _transverse(neighbors[0], map_nodo_seq_start, map_seq_nodo_end,map_seqs,sequence);
    }else
    {
        map_seq_nodo_end[seg] = n.first;
        if (!neighbors.size())
        {
            for (uint i = 0; i < kmer_size-1; ++i)
                sequence.append_nuc_right(nodeInfo.node.at(i));
        }
        map_seqs[seg++] = sequence;
        if (neighbors.size() > 1)
        {
            cout << "OutDegree: "<<_g[n.first].node.str()<<"\n";
            if (map_nodo_seq_start.find(n.first) == map_nodo_seq_start.end())
            {
                map_nodo_seq_start[n.first] = vector<size_t>();
                for (auto n2: neighbors)
                {
                    DnaSequence seq("");
                    seq.append_nuc_right(nodeInfo.node.at(0));
                    map_nodo_seq_start[n.first].push_back(seg);
                    _transverse(n2,map_nodo_seq_start, map_seq_nodo_end,map_seqs,seq);
                }
            }else
                map_seq_nodo_end[seg-1] = n.first;
        }
        if (in_nodes > 1)
        {
            cout << "Indegree: "<<_g[n.first].node.str()<<"\n";
            if (map_nodo_seq_start.find(n.first) != map_nodo_seq_start.end())
                map_seq_nodo_end[seg-1] = n.first;
            else
            {
                cout << "Launching: "<<neighbors.size()<<"\n";
                /*
                 * Cabe la posibilidad de 1 vecino con este haplotipo pero varios vecinos con todos los haplotipos posibles.
                 */
                if (neighbors.size())
                {
                    map_nodo_seq_start[n.first] = vector<size_t>();
                    for (auto n2:neighbors)
                    {
                        DnaSequence seq("");
                        seq.append_nuc_right(nodeInfo.node.at(0));
                        map_nodo_seq_start[n.first] = vector<size_t>(1,seg);
                        _transverse(n2, map_nodo_seq_start, map_seq_nodo_end, map_seqs,seq);
                    }
                }
            }
        }
    }
}

template<>
void boostDBG<true>::extension(vector <Node> in_0, string path_to_write) {
    cout << "--------------------------------\nBoostDBG extension\n--------------------------------\n";
    map <graphBU, vector<size_t>> map_nodo_seq_start;
    map <size_t, graphBU> map_seq_nodo_end;
    map <size_t, DnaSequence> map_seqs;
    for (auto n:_in_0) {
        NodeInfo node_info = _g[n];
        map_nodo_seq_start[n] = vector<size_t>();
        for (auto n2: getKmerNeighbors(n)) {
            NodeInfo neigh_info = _g[n2];
            //Cambiar neigh_info.node_set -> por haplotipo del padre
            for (auto hap: neigh_info.parent_cliques[node_info.node])
            {
                cout << "Hap 1\n";
                DnaSequence sequence("");
                cout << "New launching" << "\n";
                sequence.append_nuc_right(node_info.node.at(0));
                map_nodo_seq_start[n].push_back(seg);
                _transverse(pair<graphBU, PairedInfoNode>(n2, hap), map_nodo_seq_start, map_seq_nodo_end, map_seqs,
                            sequence);
            }
        }
    }
    /*
    cout << "Seq_start: \n";
    for (auto s:map_nodo_seq_start){
        cout << _g[s.first].node.str()<<": ";
        for (auto t:s.second) {
            cout << t << "\t";
        }
        cout << "\n";
    }
    cout << "Seq end: \n";
    for (auto s:map_seq_nodo_end){
        cout << s.first <<":"<<_g[s.second].node.str();
        cout << "\n";
    }*/
    cout << "Writing unitigs\n";
    _writeUnitigs(map_nodo_seq_start,map_seq_nodo_end, map_seqs, path_to_write);
}

template<>
boostDBG<true>::boostDBG(DBG<true> * dbg)
{
    std::cout << "STAGE: Filling boost graph\n";
    map<Node, vertex_t > local_map;
    pair<unordered_set<Node>, unordered_set<Node>> graph_struct = dbg->getNodes();
    for (auto k: graph_struct.second)
    {
        vector <Node> neigh = dbg->getKmerNeighbors(k);
        vertex_t origin, target;
        cout << "ExtraInfo: "<<k.str()<<"\n";
        if (local_map.find(k) == local_map.end())
        {
            pair<bool, PairedInfoNode> extraInfo = dbg->getExtra(k);
            cout << "ExtraInfo first "<<extraInfo.first<<"\n";
            for (auto t: extraInfo.second)
                cout << t.str()<<"\n";
            if (extraInfo.first)
                origin = boost::add_vertex(NodeInfo(k, _node_id++, extraInfo.second),_g);
            else
                origin = boost::add_vertex(NodeInfo(k, _node_id++),_g);
            local_map[k] = origin;
        }else
            origin = local_map[k];
        for (auto k2: neigh) {
            if (local_map.find(k2) == local_map.end())
            {
                pair<bool, PairedInfoNode> extraInfo = dbg->getExtra(k2);
                if (extraInfo.first)
                {
                    target = boost::add_vertex(NodeInfo(k2,_node_id++,extraInfo.second),_g);
                }
                else
                    target = boost::add_vertex(NodeInfo(k2, _node_id++),_g);
                _g[target].parent_cliques[k] = vector<PairedInfoNode>();
                local_map[k2] = target;
            }else {
                target = local_map[k2];
                _g[target].parent_cliques[k] = vector<PairedInfoNode>();
            }
            edge_t e = boost::add_edge(origin, target, _g).first;
            _g[e] = EdgeInfo(k2.at(Parameters::get().kmerSize-2), _edge_id++);
        }
    }
    _insertExtraInfo(_g, local_map);
    _modify_info();
    //GetInDegree 0
    //show_info();
    vertex_iterator v, vend;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        if (in_degree(*v) == 0)
        {
            _in_0.push_back(*v);
        }
    }
    cout << "Suspicious Starts: "<<_in_0.size()<<"\n";
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
    vector<DnaSequence> seqs;
    vector<pair<size_t,size_t>> links;
    vector<bool> strand;
    for (auto s:map_seqs)
        seqs.push_back(s.second);
    for (auto s:map_nodo_seq_start)
        for (auto origin: s.second)
            for (auto end: map_nodo_seq_start[map_seq_nodo_end[origin]])
                links.push_back(pair<size_t,size_t>(end, origin));
    UnitigExtender<false>::_write_gfa(file_to_path, seqs,links);
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
    pair<unordered_set<Node>,unordered_set<Node>> graph_shape = dbg->getNodes();
    _solid_kmers = graph_shape.first;
    for (auto n: graph_shape.second)
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
    /*for (auto n:_solid_kmers)
        cout << "Kmer: "<<n.str()<<":"<<n.hash()<<"\n";
    cout << "Start Kmers: "<<_in_0.size()<<"\n";*/
    //show_info();

    /*for (auto n: _in_0)
        cout << "StartPoint: "<<n.str()<<"\n";*/
}

template<>
listDBG<false>::listDBG(DBG<false> * dbg)
{
    _buildNewGraph(dbg);
}
