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
                if (is_solid(sep_nodes.first) && is_solid(sep_nodes.second))
                {
                    Node origin = sep_nodes.first.substr(0, Parameters::get().kmerSize-1),
                            target = sep_nodes.first.substr(1,Parameters::get().kmerSize),
                            origin_right = sep_nodes.second.substr(0, Parameters::get().kmerSize-1),
                            target_right = sep_nodes.second.substr(1,Parameters::get().kmerSize);
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
                    if (_is_standard)
                    {
                        Node k_left_rc = sep_nodes.first.rc(), k_right_rc = sep_nodes.second.rc();
                        Node origin_rc = k_left_rc.substr(0, Parameters::get().kmerSize-1),
                                target_rc = k_left_rc.substr(1, Parameters::get().kmerSize),
                                origin_right_rc = k_right_rc.substr(0, Parameters::get().kmerSize-1),
                                target_right_rc = k_right_rc.substr(1, Parameters::get().kmerSize);
                        _extra_info.insert(target_right_rc, target_rc);
                        _extra_info.insert(origin_right_rc, origin_rc);
                    }
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
            dist[_g[*v].id*num_vertex+i] = INF;
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
                if (dist[i * num_vertex + k] + dist[k * num_vertex + j] < dist[i * num_vertex + j])
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
void boostDBG<true>::_modify_info()
{
    /*
     * We are going to use pair_info (ExtraInfoNodes) to "modify" the graph including this new information
     */
    vertex_iterator v, vend;
    size_t num_vertex = boost::num_vertices(_g);
    int * distance_matrix = _floyds_warshall();
    //New_graph container
    Graph tmp_graph;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges((*v), _g);
        NodeInfo node_info = _g[*v];
        /*
         * Mirar como corregir
         */
        std::cout << "Kmer: "<<_g[*v].node.str()<<" ";
        for (; neighbors.first != neighbors.second; ++neighbors.first)
        {
            /*
             * New graph definition:
             *      - BidirectionalS/UndirectedS
             */
            typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, NodeInfo> Graph_l;
            typedef Graph_l::vertex_descriptor vertex_graph;
            typedef Graph_l::vertex_iterator vertex_it;
            /*
             * Local Graph which joins evey single reachable node (from the others)
             */
            auto endpoint = boost::target(*neighbors.first,_g);

            NodeInfo neigh_info = _g[endpoint];
            unordered_set<vertex_t> local_vect;
            vector<bool> visited(num_vertex, false);
            Graph_l local_graph;
            size_t curr_node = 0;
            map<Node, vertex_graph> local_node_map;
            std::cout << "Neighbor: "<<neigh_info.node.str()<<"\n";
            if (_map_extra_info[node_info.id].empty())
            {
                break;
            }
            if (_map_extra_info[neigh_info.id].empty())
            {
                continue;
            }
            local_vect = getUnion(_map_extra_info[node_info.id], _map_extra_info[neigh_info.id]);
            for (auto s:local_vect)
            {
                NodeInfo s_info = _g[s];
                vertex_graph source;
                if (local_node_map.find(s_info.node) == local_node_map.end())
                {
                    source = boost::add_vertex(NodeInfo(s_info.node,curr_node++), local_graph);
                    local_node_map[s_info.node] = source;
                }else
                    source = local_node_map[s_info.node];
                for (auto t:local_vect)
                {
                    NodeInfo t_info = _g[t];
                    if (visited[t_info.id] || t_info.id == s_info.id)
                        continue;
                    std::cout << " " << s_info.node.str() << " "<<t_info.node.str();
                    std::cout << " " << _reachable(distance_matrix, s_info.id, t_info.id)<<"\n";
                    if (_reachable(distance_matrix, s_info.id, t_info.id))
                    {
                        vertex_graph target;
                        if (local_node_map.find(t_info.node) == local_node_map.end())
                        {
                            target = boost::add_vertex(NodeInfo(t_info.node, curr_node++), local_graph);
                            local_node_map[t_info.node] = target;
                        }else
                            target = local_node_map[t_info.node];
                        boost::add_edge(source, target, local_graph);
                        visited[s_info.id] = 1;
                    }
                }
            }
            /*
             * Calculate all maximal cliques in the local graph
             */
            vector<vector<vertex_graph>> output = findMaxClique<Graph_l, vertex_graph,vertex_it>(local_graph);
            while (boost::num_edges(local_graph)) {
                for (auto clique: output)
                {
                    ExtraInfoNode local_haplotype;
                    for (size_t i = 0; i < clique.size(); ++i)
                    {
                        local_haplotype.emplace(local_graph[clique[i]].node);
                        for (size_t j = i + 1; j < clique.size(); ++j) {
                            local_haplotype.emplace(local_graph[clique[j]].node);
                            cout << local_graph[clique[i]].id << " " << local_graph[clique[j]].id << "\n";
                            boost::remove_edge(clique[i], clique[j], local_graph);
                        }
                    }
                    /*
                     * Add cliques (presumed haplotypes which contains both nodes)
                     */
                    _g[endpoint].parent_cliques[node_info.node].push_back(local_haplotype);
                }
                output = findMaxClique<Graph_l, vertex_graph, vertex_it>(local_graph);
            }
            /*
             * Get in_neighbors
             */
            pair<in_iterator,in_iterator> in_neighbors =boost::in_edges((*v), _g);
            cout << "InNeighbors: "<<_g[*v].node.str()<<"\n";
            for (; in_neighbors.first != in_neighbors.second; ++in_neighbors.first)
            {
                cout << _g[boost::source(*in_neighbors.first, _g)].node.str() << "\n";
            }
        }
    }
    //Free Floyd
    free(distance_matrix);
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
    size_t kmer_size = Parameters::get().kmerSize;
    cout << "Kmer: "<<nodeInfo.node.str()<<"\n";
    cout << "InDegree size: "<< in_degree(n)<< " Esperable 1\n";
    vector<it_node> neighbors = getOutKmerNeighbors(n);
    cout << "Neighbor size: "<< neighbors.size()<< " Esperable 1\n";
    for (auto t:n.second)
        cout << t.str()<<"\n";
    if ((neighbors.size()==1) && (in_degree(n)==0))
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
                map_nodo_seq_start[n.first].push_back(seg);
        }
        if (in_degree(n) > 0)
        {
            if (map_nodo_seq_start.find(n.first) != map_nodo_seq_start.end())
                map_nodo_seq_start[n.first].push_back(seg);
            else
            {
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
                        map_nodo_seq_start[n.first].push_back(seg);
                        _transverse(n2, map_nodo_seq_start, map_seq_nodo_end, map_seqs,seq);
                    }
                }
            }
        }
    }
}

template<>
void boostDBG<true>::extension(vector <Node> in_0, string path_to_write)
{
    cout << "--------------------------------\nBoostDBG extension\n--------------------------------\n";
    map<graphBU, vector<size_t>> map_nodo_seq_start;
    map<size_t, graphBU> map_seq_nodo_end;
    map<size_t, DnaSequence> map_seqs;
    for(auto n:_in_0)
    {
        NodeInfo node_info = _g[n];
        map_nodo_seq_start[n] = vector<size_t>();
        for (auto n2: getKmerNeighbors(n))
        {
            NodeInfo neigh_info = _g[n2];
            for (auto p: neigh_info.parent_cliques[node_info.node])
            {
                DnaSequence sequence("");
                sequence.append_nuc_right(node_info.node.at(0));
                map_nodo_seq_start[n].push_back(seg);
                _transverse(pair<graphBU, ExtraInfoNode>(n2,p)
                        , map_nodo_seq_start,map_seq_nodo_end,map_seqs,sequence);
            }
        }
    }
    for (auto s:map_seqs)
        cout << "Sequence: "<<s.second.str()<<"\n";
}

template<>
boostDBG<true>::boostDBG(DBG<true> * dbg)
{
    std::cout << "Trying to fill the graph: \n";
    map<Node, vertex_t > local_map;
    pair<unordered_set<Node>, unordered_set<Node>> graph_struct = dbg->getNodes();
    for (auto k: graph_struct.second)
    {
        vector <Node> neigh = dbg->getKmerNeighbors(k);
        vertex_t origin, target;
        if (local_map.find(k) == local_map.end())
        {
            pair<bool, ExtraInfoNode> extraInfo = dbg->getExtra(k);
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
                pair<bool, ExtraInfoNode> extraInfo = dbg->getExtra(k2);
                if (extraInfo.first)
                    target = boost::add_vertex(NodeInfo(k2,_node_id++,extraInfo.second),_g);
                else
                    target = boost::add_vertex(NodeInfo(k2, _node_id++),_g);
                _g[target].parent_cliques[k] = vector<ExtraInfoNode>();
                local_map[k2] = target;
            }else {
                target = local_map[k2];
                _g[target].parent_cliques[k] = vector<ExtraInfoNode>();
            }
            edge_t e = boost::add_edge(origin, target, _g).first;
            _g[e] = EdgeInfo(k2.at(Parameters::get().kmerSize-2), _edge_id++);
        }
    }
    _insertExtraInfo(_g);
    _modify_info();
    //GetInDegree 0
    show_info();
    vertex_iterator v, vend;
    for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
    {
        if (!in_degree(_g[*v].node))
        {
            cout << "Kmer: "<<_g[*v].node.str()<<"\n";
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
        cout << "Sequence: "<<sequence.str()<<"\n";
        _transverse(neighbors[0], map_nodo_seq_start, map_seq_nodo_end,map_seqs,sequence);
    }else
    {
        cout << "SEG: "<<seg<< " NODE: "<<n.str()<<"\n";
        map_seq_nodo_end[seg] = n;
        if (neighbors.size()==0)
        {
            for (uint i = 0; i < kmer_size; ++i)
                sequence.append_nuc_right(n.at(i));
            cout <<"Finalizo en :"<<n.str()<<"-----> Sequence: " <<sequence.str()<<"\n";
        }
        map_seqs[seg++] = sequence;
        if (neighbors.size() > 1)
        {
            cout << "Kmer with more than 1 outDegree: "<<n.str()<<"\n";
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
            cout << "Kmer with more than 1 inDegree: "<<n.str()<<"\n";
            if (map_nodo_seq_start.find(n) != map_nodo_seq_start.end()){
                map_seq_nodo_end[seg-1] = n;
            }
            else
            {
                cout << "Soy: "<<n.str() << "---> Extiendo\n";
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
        cout << "KMER START: "<<n.str()<<"\n";
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
    /*for (auto n:_solid_kmers)
        cout << "Kmer: "<<n.str()<<"\n";*/
    cout << "Start Kmers: "<<_in_0.size()<<"\n";
    //show_info();
    cout << "SOLID INFORMATION\n";

    /*for (auto n: _in_0)
        cout << "StartPoint: "<<n.str()<<"\n";*/
}

template<>
listDBG<false>::listDBG(DBG<false> * dbg)
{
    _buildNewGraph(dbg);
}
