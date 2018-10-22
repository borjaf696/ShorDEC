#include <unordered_map>
#include <map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <boost/version.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include "../Extender/Extender.h"

//Constants
#define MIN_PATH_LEN 1
#define DELTA_PATH_LEN 4

using namespace std;
template<typename T>
class Internal
{
    Internal(T n):_n(n){}
    T getNode(){return _n;}
private:
    T _n;
};

template<bool P>
class NaiveDBG: public DBG<P>
{
public:
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
    typedef typename DBG<P>::Parent_Extra ExtraInfoNode;
    typedef typename BUgraph<Node>::graphBU graphBU;
    typedef typename BUgraph<FuncNode>::graphBU graphBU_Func;

    size_t in_degree(Node);
    size_t out_degree(Node);
    vector<DnaSequence::NuclType> getNeighbors (Node) const;
    vector<Node> getKmerNeighbors
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

    /*
     * Define how to check if rc or just forward
     */
    NaiveDBG(SequenceContainer& sc):_sc(sc)
    {
        Progress::get().size_total = _sc.getIndex().size();
        Progress::get().show = true;
        _kmerCount();
        _cleaning();
    }
    /*
     * Check whether a kmer is solid or not
     */
    bool is_solid(Node& kmer) const
    {
        if (_is_standard) {
            Node kmer_aux = kmer;
            kmer_aux.standard();
            return (_dbg_naive.find(kmer_aux) != _dbg_naive.end());
        }else
            return (_dbg_naive.find(kmer) != _dbg_naive.end());
    }
    /*
     * "Length" of the DBG
     */
    size_t length() const
    {
        return _dbg_naive.size();
    }

    typename DBG<P>::Heads get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }

    void ProcessTigs(string path_to_write)
    {
        std::cout << "Lets start: Naive\n";
        extension(_in_0,path_to_write);
        std::cout << "End Unitigs: Naive\n";
    }

    void show_info();

    //Extension
    void extension(vector<Node> in_0, string path_to_write)
    {
        size_t _curr_segment = 0;
        unordered_map<Node, vector<size_t>> _fin_segs;
        vector<pair<size_t,size_t>> _links;
        //Sequences
        vector<DnaSequence> _seqs;
        std::cout << "Voy a extender como un pro!\n";
        /*
         * Set of already assesed heads
         */
        unordered_set<Node> added;
        /*
         * New heads with out and in > 1
         */
        stack<graphBU> in, out;
        /*
         * Unitigs: entran por orden de curr_segment, ¡Nos ahorramos indexar las secuencias!
         */
        vector<vector<Node>> unitigs;
        for (auto k: in_0)
        {
            for (auto &p: this->extend(k,out,in,added,_curr_segment,_fin_segs))
                unitigs.push_back(p);
        }
        //Remove the kmers from in_0
        in_0.clear();
        while (!in.empty() && !out.empty())
        {
            /*
             * Lets check kmers with out_degree > 1
             */
            while (!out.empty())
            {
                Node k = out.top();
                out.pop();
                Node n_k = getNode(k);
                if (_fin_segs.find(n_k) != _fin_segs.end())
                    for (uint i = 0; i < _fin_segs[n_k].size(); i++) {
                        _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[n_k][i]));
                    }
                for (auto &p: this->extend(k,out,in,added, _curr_segment, _fin_segs))
                    unitigs.push_back(p);
            }
            /*
             * Lets check kmers with in_degree > 1
             */
            while(!in.empty())
            {
                Node k = in.top();
                in.pop();
                Node n_k = getNode(k);
                if (_fin_segs.find(n_k) != _fin_segs.end())
                    for (uint i = 0; i < _fin_segs[n_k].size(); i++) {
                        _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[n_k][i]));
                    }
                for (auto &p: this->extend(k,out,in,added,_curr_segment,_fin_segs))
                    unitigs.push_back(p);
            }
        }
        /*
         * Construct the DnaSequences
         */
        UnitigExtender<P>::_construct_sequences(unitigs,_fin_segs,_seqs);
        /*
         * Lets write them
         */
        UnitigExtender<P>::_write_gfa(path_to_write,_seqs,_links);
    }

    //Get
    vector<Node> getEngagers()
    {
        return _in_0;
    }
    pair<unordered_set<Node>, unordered_set<Node>> getNodes()
    {
        return pair<unordered_set<Node>, unordered_set<Node>>(_dbg_naive, _dbg_nodes);
    };

    pair<bool,ExtraInfoNode> getExtra(Node node)
    {
        return _extra_info.getInfoNode(node);
    }

    Node getNode(Node t)
    {
        return t;
    }

    //Operators
    NaiveDBG& operator=(const NaiveDBG& other)
    {
        if (this->length()!= other.length())
        {
            this->_sc = other._sc;
            _dbg_naive = other._dbg_naive;
            _heads = other._heads;
            _tails = other._tails;
            _in_0 = other._in_0;
        }
        return *this;
    }

private:
    /*
     * Kmer Counting
     * Naive DBG construction + heads + tails
     */
    void _kmerCount();
    void _remove_isolated_nodes();
    /*
     * Only paired_version gives impl
     */
    void _insert_extra_info();
    void _to_pair_end();

    void _cleaning()
    {
        _remove_isolated_nodes();
    }
    /*
     * Insertion into the graph_nodes and graph_edges
     */
    void _insert(Node k, Node kmer)
    {
        Node origin = kmer.substr(0, Parameters::get().kmerSize-1),
                target = kmer.substr(1, Parameters::get().kmerSize);
        _dbg_naive.emplace(k);
        _dbg_nodes.emplace(origin);
        _dbg_nodes.emplace(target);
        if (_is_standard)
        {
            Node rc = kmer.rc();
            Node origin_rc = rc.substr(0, Parameters::get().kmerSize-1),
                    target_rc = rc.substr(1, Parameters::get().kmerSize);
            _dbg_nodes.emplace(origin_rc);
            _dbg_nodes.emplace(target_rc);
        }
    }
    /*
     * Check real neighbors:
     */
    void _build_pair_neighs(Node node, unordered_set<Node>& set_full) const
    {
        for (uint i = 0; i < 4; i++)
        {
            Node new_node = node;
            new_node.appendRightReplace(i);
            set_full.emplace(new_node);
        }
    }

    vector<Node> _check_valid(vector<Node> neighbors,Node node) const
    {
        vector<Node> real_neighbors;
        if (!_extra_info.find(node))
            return real_neighbors;
        unordered_set<Node> node_pairs = _extra_info[node], set_full;
        for (auto s:node_pairs)
            _build_pair_neighs(s, set_full);
        for (auto k: neighbors)
        {
            unordered_set<Node> result = getIntersection(set_full, _extra_info[k]);
            if (!result.empty())
                real_neighbors.push_back(k);
        }
        return real_neighbors;
    }

    vector<DnaSequence> _get_sequences(vector<vector<Node>> unitigs)
    {
        /*
         * From the vector of unitigs get all the unitigs
         */
        cout << "Unitigs\n";
        vector<DnaSequence> dna_vect;
        for (auto &vect: unitigs){
            cout << "New Unitig: \n";
            size_t cont = 0;
            DnaSequence seq_build(vect[0].str());
            for (auto &k: vect) {
                if (cont)
                    seq_build.append_nuc_right(k.at(Parameters::get().kmerSize-1));
                cont++;
            }
            dna_vect.push_back(seq_build);
            cout << "FinalSequence: "<<seq_build.str() <<"\n";
        }
        return dna_vect;
    }
    /*
     * k1->k2->k3 (To standard post append) -> Only length matters :P
     */
    void _check_forward_path(size_t& len_fw, vector<Node>& k_vec) const
    {
        /*
         * Check neighbors of the Kmer
         */
        Node aux = k_vec.back();
        std::vector<DnaSequence::NuclType> neigh_fw = getNeighbors(aux);
        if (neigh_fw.size() == 1) {
            len_fw++;
            if (len_fw < MIN_PATH_LEN) {
                Node node_aux = aux;
                node_aux.appendRightReplace(neigh_fw[0]);
                k_vec.push_back(node_aux);
                _check_forward_path(len_fw, k_vec);
            }
        }
        if (neigh_fw.size() > 1)
            len_fw += MIN_PATH_LEN;
    }

    bool _asses(vector<Node> &erase,vector<Node> aux, size_t len)
    {
        if (len < MIN_PATH_LEN)
            for (auto k:aux) {
                erase.push_back(k);
            }
        return (len < MIN_PATH_LEN);
    }

    void _erase(vector<Node>& kmer_to_erase)
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
        }
        /*cout << "NaiveSizePost: "<<_dbg_naive.size() << "\n";
        cout << "NodesSize: "<<_dbg_nodes.size()<<"\n";*/
        kmer_to_erase.clear();
    }

    /*
     * I/O
     */
    void _write_unitigs(vector<DnaSequence> dnaSequences, string filename)
    {
        FILE* fout = fopen(filename.c_str(), "w");
        if (!fout)
            throw std::runtime_error("Can't open " + filename);
        //size_t num_unitig = 0;
        for (auto& seq : dnaSequences)
        {
            std::string seq_ =seq.str()+"\n";
            std::string header = ">NumUnitig\n";
            fwrite(header.data(), sizeof(header.data()[0]),
                   header.size(), fout);
            fwrite(seq_.data(), sizeof(seq_.data()[0]),
                   seq_.size(), fout);
        }
    }
    /*
     * First Counter
     */
    unordered_map<Node, pair<size_t,size_t>> _kmers_map;
    /*
     * PairedInfo
     */
    Extra<P> _extra_info;
    /*
     * DBG_naive -> stores the set of solid Kmers
     * DBG_nodes -> stores the set of (K-1)mers
     * TODO: Pack (Extra, _dbg_naive, _dbg_nodes) under the same struct and define template
     */
    unordered_set<Node> _dbg_naive, _dbg_nodes;
    unordered_set<KmerInfo<P>> _heads,_tails;
    //Extension points
    vector<Node> _in_0;
    //Extend
    SequenceContainer& _sc;
    //Standard
    bool _is_standard = false;
};
/*
 * Boost implementation
 */
template <bool P> class boostDBG:public DBG<P>
{
public:
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
    typedef typename NodeType<P>::set_couples ExtraInfoNode;
    /*
     * Boost structure
     */
    /*
     * SingleEnd
     */
    //TODO: Check how to manage this thing
    /*
     * PairEnd
     */
    struct NodeInfo {
        NodeInfo():id(-1){}
        NodeInfo(Node node, int32_t id):node(node), id(id){}
        NodeInfo(Node node, int32_t id, ExtraInfoNode extra):node(node), id(id),node_set(extra){}
        NodeInfo(const NodeInfo & nodeInfo):node(nodeInfo.node),id(nodeInfo.id),node_set(nodeInfo.node_set){}
        NodeInfo& operator=(const NodeInfo& other)
        {
            node = other.node;
            return *this;
        }
        bool equal(Node node) const
        {
            return node == this->node;
        }
        bool operator ==(const NodeInfo &other) const
        {
            return equal(other.node) && (this->id == other.id);
        }
        Node node;
        int32_t id;
        ExtraInfoNode node_set;
        map<Node,vector<ExtraInfoNode>> parent_cliques;
    };

    struct EdgeInfo {
        EdgeInfo(){}
        EdgeInfo(int8_t nt_received, size_t id):nt(nt_received),id(id){}
        EdgeInfo& operator=(const EdgeInfo& other)
        {
            nt = other.nt;
            return *this;
        }
        int8_t nt;
        size_t id;
    };
    struct GraphInfo{};
    typedef typename boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, NodeInfo,
            EdgeInfo> Graph;
    /*
     * Bundles for properties
     */
    typedef typename Graph::vertex_descriptor vertex_t;
    typedef typename Graph::edge_descriptor edge_t;
    typedef typename Graph::adjacency_iterator adjacency_iterator;
    typedef typename Graph::out_edge_iterator out_iterator;
    typedef typename Graph::in_edge_iterator in_iterator;
    typedef typename Graph::vertex_iterator vertex_iterator;
    typedef typename BUgraph<vertex_t>::graphBU graphBU;
    /*
     * TODO: Needs implementation
     */
    boostDBG(DBG<P> *);
    bool is_solid(Node&) const
    {
        return true;
    }
    size_t length() const
    {
        return 0;
    }
    vector<typename DnaSequence::NuclType> getNeighbors
            (Node) const
    {
        vector<DnaSequence::NuclType> neigh;
        return vector<typename DnaSequence::NuclType>();
    }

    vector<Node> getKmerNeighbors(Node node) const
    {
        vector<Node> neigh;
        vertex_t vertex = _getNode(node);
        if (vertex == nullptr)
            return neigh;
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges(vertex, _g);
        for(; neighbors.first != neighbors.second; ++neighbors.first)
        {
            neigh.push_back(_g[boost::target(*neighbors.first,_g)].node);
        }
        return neigh;
    }

    vector<graphBU> getKmerNeighbors(graphBU node) const
    {
        vector<graphBU> neigh;
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges(node, _g);
        for (; neighbors.first != neighbors.second; ++neighbors.first)
        {
            neigh.push_back(boost::target(*neighbors.first,_g));
        }
        return neigh;
    }

    size_t in_degree(Node node)
    {
        graphBU v_node = this->_getNode(node);
        return boost::in_degree(v_node, _g);
    }

    size_t in_degree(graphBU node)
    {
        return boost::in_degree(node, _g);
    }

    size_t out_degree(Node node)
    {
        graphBU v_node = this->_getNode(node);
        return boost::degree(v_node, _g);
    }

    size_t out_degree(graphBU node)
    {
        return boost::degree(node, _g);
    }

    //Extesion ->TODO: move
    vector<vector<Node>> extend(graphBU kmer,
                                       stack<graphBU > &out,
                                       stack<graphBU > &in,
                                       unordered_set<graphBU > added,
                                       size_t &curr_segment,
                                       unordered_map<Node, vector<size_t>> &fin_segs)
    {
        vector<vector<Node>> unitigs;
        Node k_node = getNode(kmer);
        vector<graphBU> neighbors = getKmerNeighbors(kmer);
        for (auto &k: neighbors)
        {
            vector<Node> unitig;
            unitig.push_back(k_node);
            pair<size_t,Node> result =  _Extension(k,unitig,out, in, added);
            if (result.first == 1 || result.first == 2)
            {
                fin_segs[result.second].push_back(curr_segment);
            }
            curr_segment++;
            unitigs.push_back(unitig);
        }
        return unitigs;
    }

    void extension(vector<Node> in_0, string path_to_write)
    {
        //Extension info:
        size_t _curr_segment = 0;
        unordered_map<Node, vector<size_t>> _fin_segs;
        vector<pair<size_t,size_t>> _links;
        vector<DnaSequence> _seqs;

        std::cout << "Pair_End extension!\n";
        /*
         * Set of already assesed heads
         */
        unordered_set<graphBU> added;
        /*
         * New heads with out and in > 1
         */
        stack<graphBU> in, out;
        /*
         * Unitigs: entran por orden de curr_segment, ¡Nos ahorramos indexar las secuencias!
         */
        vector<vector<Node>> unitigs;
        for (auto k: in_0)
        {
            //Searching the node -> TODO: ¿Change?
            graphBU k_node = _getNode(k);
            for (auto &p: extend(k_node,out,in,added,_curr_segment,_fin_segs))
                unitigs.push_back(p);
        }
        //Remove the kmers from in_0
        in_0.clear();
        while (!in.empty() && !out.empty())
        {
            /*
             * Lets check kmers with out_degree > 1
             */
            while (!out.empty())
            {
                graphBU k = out.top();
                out.pop();
                Node n_k = getNode(k);
                if (_fin_segs.find(n_k) != _fin_segs.end())
                    for (uint i = 0; i < _fin_segs[n_k].size(); i++) {
                        _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[n_k][i]));
                    }
                for (auto &p: extend(k,out,in,added, _curr_segment, _fin_segs))
                    unitigs.push_back(p);
            }
            /*
             * Lets check kmers with in_degree > 1
             */
            while(!in.empty())
            {
                graphBU k = in.top();
                in.pop();
                Node n_k = getNode(k);
                if (_fin_segs.find(n_k) != _fin_segs.end())
                    for (uint i = 0; i < _fin_segs[n_k].size(); i++) {
                        _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[n_k][i]));
                    }
                for (auto &p: extend(k,out,in,added,_curr_segment,_fin_segs))
                    unitigs.push_back(p);
            }
        }
        /*
         * Construct the DnaSequences
         */
        UnitigExtender<P>::_construct_sequences(unitigs,_fin_segs,_seqs);
        /*
         * Lets write them
         */
        UnitigExtender<P>::_write_gfa(path_to_write,_seqs,_links);
    }

    //Getter
    vector<Node> getEngagers()
    {
        return _in_0;
    }

    typename DBG<P>::Heads get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }

    pair<unordered_set<Node>, unordered_set<Node>> getNodes()
    {
        return pair<unordered_set<Node>, unordered_set<Node>>(unordered_set<Node>(), unordered_set<Node>());
    }

    pair<bool,ExtraInfoNode> getExtra(Node node)
    {
        /*
         * Iterate over all nodes and return the paired info for the Node = node
         */
        return pair<bool,ExtraInfoNode>(false,ExtraInfoNode ());
    }

    Node getNode(Node t)
    {
        return t;
    }

    Node getNode(graphBU k)
    {
        return _g[k].node;
    }

    void ProcessTigs(string path_to_write)
    {
        cout << "Lets try Unitigs Pair_End graph\n";
        //UnitigExtender<P, graphBU>::full_extension(this, vector<graphBU>(), path_to_write);
        cout << "Unitigs donette\n";
    }
    //Show methods
    void show_info()
    {
        typename Graph::vertex_iterator v, vend;
        for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v) {
            std::cout << " Kmer:"     << _g[*v].node.str()
                      << " id:"  << _g[*v].id
                      << "\n";
            vector<Node> neighbors = getKmerNeighbors(_g[*v].node);
            std::cout << "Neighbors: ";
            for (auto n:neighbors)
                std::cout  << n.str()<<" ";
            std::cout << "\n";
        }
    }
private:
    //Extension:
    pair<size_t, Node> _Extension(graphBU kmer,vector<Node> & unitig,
                                         stack<graphBU> & out, stack<graphBU> & in,
                                         unordered_set<graphBU> added)
    {
        vector<graphBU> neighbors = getKmerNeighbors(kmer);
        size_t in_ = in_degree(kmer);
        Node kmer_node = getNode(kmer);
        if (in_ == 1 && neighbors.size() == 1) {
            unitig.push_back(kmer_node);
            return _Extension(neighbors[0],unitig,out,in,added);
        }else if (neighbors.size () == 0) {
            unitig.push_back(kmer_node);
            return pair<size_t, Node>(0, kmer_node);
        }else if (added.find(kmer) == added.end())
        {
            added.emplace(kmer);
            if (neighbors.size() > 1) {
                out.push(kmer);
                unitig.push_back(kmer_node);
                return pair<size_t,Node>(1,kmer_node);
            }
            in.push(kmer);
            unitig.push_back(kmer_node);
            return pair<size_t, Node>(2,kmer_node);
        }
        return pair<size_t, Node>(0,kmer_node);
    }
    //Boost Get Node:
    vector<graphBU> _getRealNeighbors(graphBU node){

        return vector<graphBU>();
    }
    graphBU _getNode(Node node) const
    {
        vertex_iterator v, vend;
        for (boost::tie(v, vend) = boost::vertices(_g); v != vend; ++v)
        {
            if (_g[*v].equal(node))
            {
                return (*v);
            }
        }
        return nullptr;
    }

    template<typename T>
    void _insertExtraInfo(const T &g)
    {
        std::cout << "Import extrainfo as vertex_t: "<<_node_id<<"\n";
        vertex_iterator v, vend;
        _map_extra_info.resize(_node_id);
        for (boost::tie(v,vend) = boost::vertices(g); v != vend; ++v)
        {
            std::cout << g[*v].node.str()<<"\n";
            NodeInfo node_info = g[*v];
            if (!node_info.node_set.empty())
            {
                unordered_set<vertex_t> local;
                for (auto s:g[*v].node_set)
                {
                    std::cout << s.str()<<"\n";
                    vertex_t pointerToExtraInfo = _getNode(s);
                    std::cout << g[pointerToExtraInfo].node.str()<<"\n";
                    local.emplace(pointerToExtraInfo);
                }
                std::cout <<"Cualquier cosa: "<<g[*v].id<<" InfoSize: "<<_map_extra_info.size()
                          <<"\n";
                _map_extra_info[g[*v].id] = local;
            }
        }
        std::cout << "ExtraInfo imported as vertex_t\n";
    }
    /*
     * Completing info
     */
    int* _floyds_warshall();
    bool _reachable(int*, size_t, size_t);
    void _modify_info();
    /*auto _add_vertex(int32_t id, Node node);
    auto _add_edge(int32_t, int8_t);*/
    void _kmerCount()
    {}
    void _cleaning()
    {}
    //Graph
    Graph _g;
    //Properties + nodes (this should be fixed with vecS)
    vector<unordered_set<vertex_t>> _map_extra_info;
    vector<Node> _in_0;
    //Node_id
    int32_t _node_id = 0, edge_id = 0;
    //Need fix
    unordered_set<KmerInfo<P>> _heads,_tails;
};