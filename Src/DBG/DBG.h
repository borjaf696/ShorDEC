#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <queue>
#include <chrono>
#include <boost/version.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include "../Extender/Extender.h"

//Constants
#define MAX_RETRIES 2
#define TIME_WAIT 2
#define MAX_SC_SIZE 2000000

#define MIN_PATH_LEN 100

#define POLISH_PATH_LEN 20
#define SHORT_LENGTH 65
#define DELTA_PATH_LEN 600
#define FLOYD 0
#define MAX_BRANCHES_CHECK 5

#define MAX_LOCAL_GRAPH 40
#define SAME_RATIO 0.2
#define DO_POLISH 0

#define MIN_SIZE_REP 0.0

#define PARENT_SON_LIMIT 0.1
#define LIMIT_TO_REP 5

#define CLICK_RATIO 0.1
#define SIMILARITY_RATIO 1.0
#define IDENTITY_RATIO 1.0
#define IDENTITY_RATIO_UNION 1.0

#define COV_RATIO 5


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
    typedef typename DBG<P>::Parent_Paired_Info PairedInfoNode;
    typedef typename DBG<P>::Parent_Freq_Map FreqMap;
    typedef typename BUgraph<Node>::graphBU graphBU;
    typedef typename BUgraph<FuncNode>::graphBU graphBU_Func;
    void clear()
    {
        _heads.clear();
        _tails.clear();
        _in_0.clear();
        _dbg_nodes.clear();
        _dbg_naive.clear();
        _extra_info.clear();
        _node_reads.clear();
    }

    size_t in_degree(Node);
    size_t out_degree(Node);
    vector<DnaSequence::NuclType> getNeighbors (Node) const;
    vector<Node> getKmerNeighbors_adhoc(Node kmer) const
    {
        vector<Kmer> nts;
        Node kmer_aux;
        for (DnaSequence::NuclType i=0; i < 4; ++i) {
            kmer_aux = Kmer(kmer);
            kmer_aux.appendRight(i);
            if (is_solid(kmer_aux)){
                nts.push_back(kmer_aux.substr(1,Parameters::get().kmerSize));
            }
        }
        return nts;
    }
    vector<Node> getKmerNeighbors
            (Node kmer) const
    {
        if (kmer.length() == Parameters::get().kmerSize)
        {
            kmer = kmer.substr(1, Parameters::get().kmerSize);
        }
        vector<Kmer> nts;
        Node kmer_aux;
        for (DnaSequence::NuclType i=0; i < 4; ++i) {
            kmer_aux = Kmer(kmer);
            kmer_aux.appendRight(i);
            if (is_solid(kmer_aux)){
                nts.push_back(kmer_aux.substr(1,Parameters::get().kmerSize));
            }
        }
        return nts;
    }

    vector<Node> getInKmerNeighbors(Node node) const
    {
        //TODO: Complete
        vector<Node> neigh;
        return neigh;
    }

    /*
     * Define how to check if rc or just forward
     */
    NaiveDBG(SequenceContainer& sc, SequenceContainer& sc_paired, bool thirdPartyCount, string path_to_file, string pair_dir, string program, bool last):_sc(sc), _sc_paired(sc_paired)
    {
        size_t trials = 0;
        Progress::get().size_total = _sc.getIndex().size();
        Progress::get().show = true;
        if (!thirdPartyCount)
            _kmerCount();
        else
            _thirdPartyKmerCounting(path_to_file, pair_dir, program, &trials, last);
        _cleaning();
        //Export solid information
        _write_solid();
        //show_info();
        /*if (P)
            _insert_extra_info();*/
    }
    /*
     * Write fasta solids
     */
    void _write_solid(){
        string path = "solids.fasta";
        ofstream solid_file;
        solid_file.open ("solids",std::fstream::in | std::fstream::out | std::fstream::trunc);
        for (auto k:_dbg_naive) {
            solid_file<<">ACCC"<<endl;
            solid_file << k.str() << endl;
        }
        solid_file.close();
        cout << "END FILE"<<endl;
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
        unordered_set<DnaSequence> _seqs;
        cout << "STAGE: Extension\n";
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
    unordered_set<Node> getSolidKmers()
    {
        return _dbg_naive;
    }

    vector<Node> getEngagers()
    {
        return _in_0;
    }

    pair<unordered_set<Node> * , unordered_set<Node> * > getNodes()
    {
        return pair<unordered_set<Node> * , unordered_set<Node> *>(&_dbg_naive, &_dbg_nodes);
    };

    unordered_map<Node, unordered_set<size_t>> getNodeReads()
    {
        return _node_reads;
    };

    pair<bool,FreqMap> getExtra(Node node)
    {
        return {false, FreqMap()};
    };

    Extra<P> * getPairedInfo(){
        return &_extra_info;
    }

    void insert(Node node, size_t freq)
    {
        if (_is_standard)
            node.standard();
        _insert(node, node, freq);
    }

    Node getNode(Node t)
    {
        return t;
    }

    SequenceContainer * getSequenceContainer()
    {
        return &(_sc_paired);
    }

    //Operators
    NaiveDBG& operator=(const NaiveDBG& other)
    {
        this->_sc = other._sc;
        _dbg_naive = other._dbg_naive;
        _dbg_nodes = other._dbg_nodes;
        _heads = other._heads;
        _tails = other._tails;
        _in_0 = other._in_0;
        return *this;
    }

private:
    /*
     * Kmer Counting
     * Naive DBG construction + heads + tails
     */
    void _kmerCount();
    //ThirdPartyCounting
    void _thirdPartyKmerCounting(string path_to_file_count, string dir_pairs, string third, size_t * retries, bool last)
    {
        vector<string> files;
        if (path_to_file_count != "")
            for (auto f: System::getAllFaFqFiles(path_to_file_count))
                files.push_back(f);
        if (dir_pairs != "")
            for (auto f: System::getAllFaFqFiles(dir_pairs))
                files.push_back(f);
        if (files.size() > 1)
        {
            cout << "Appending files:"<<std::endl;
            for (auto f:files)
                cout << "F: "<<f<<std::endl;
            path_to_file_count = System::appendFiles(files, "newFile_tmp.fasta");
        }
        else if (files.size() == 1)
            path_to_file_count = files[0];
        string instruction = "";
        if (third == "jelly")
            instruction += "bash -c \"./Utils/script/jellyfish_script ";
        else
            instruction += "bash -c \"./Utils/script/dsk_script ";
        instruction += path_to_file_count+" ";
        instruction+= to_string(Parameters::get().kmerSize)+" ";
        instruction+="output output.txt >/dev/null 2>&1\"";
        cout << "INSTRUCTION: "<<instruction<<"\n";
        while ((system(instruction.c_str())))
        {
            cout << "Problem executing: "<<instruction<<"\n";
            if ((*retries) > MAX_RETRIES)
                exit(1);
            cout << "Retrying with an alternative approach (first wait "<<TIME_WAIT<<"s)\n";
            sleep(TIME_WAIT);
            if (third == "jelly")
                third = "dsk";
            else
                third = "jelly";
            (*retries)++;
            _thirdPartyKmerCounting(path_to_file_count, dir_pairs, third, retries, last);
            return;
        }
        size_t max_freq = 0;
        if (third == "jelly")
            max_freq = createCountMapJelly<Node,size_t>(_kmers_map, "output.txt");
        else
        {
            vector<size_t> histogram = getHistogramFromFile<size_t>("histogram.txt");
            cout << "Histogram\n";
            createCountMapDSK<DBG<P>,Node>(this,"output.txt", histogram, _sc.getTotalBases()+_sc_paired.getTotalBases(),
                                           std::max(_sc.getAvLength(),_sc_paired.getAvLength()), last);
            _printInfo();
            return;
        }
        _buildGraphRepresentation(max_freq);
    }

    bool _asses(vector<Node> &erase, unordered_set<Node> & erase_nodes,
                vector<Node> aux, vector<Node> aux_nodes,size_t len)
    {
        if (len < MIN_PATH_LEN) {
            for (auto k:aux) {
                erase.push_back(k);
            }
            for (auto k:aux_nodes) {
                erase_nodes.emplace(k);
            }
        }
        return (len < MIN_PATH_LEN);
    }

    void _erase(vector<Node>& kmer_to_erase, unordered_set<Node> & nodes_erase)
    {
        /*
         * Ojo a la hebra no hay que borrar los dos llevar cuenta de que k-mer hay que borrar en cada momento... Esto puede ser peliagudo
         */
        for (auto kmer_erase:kmer_to_erase) {
            _dbg_naive.erase(kmer_erase);
        }
        for (auto node_erase:nodes_erase) {
            _dbg_nodes.erase(node_erase);
        }
        kmer_to_erase.clear();
        nodes_erase.clear();
    }

    /*
     * k1->k2->k3 (To standard post append) -> Only length matters :P
     */
    void _check_forward_path(size_t& len_fw, vector<Node>& k_vec, vector<Node> & k_node,Node node, bool show = false)
    {
        /*
         * Check neighbors of the Kmer
         */
        vector<DnaSequence::NuclType> neigh_fw = getNeighbors(node);
        if (neigh_fw.size() > 1)
            len_fw +=  MIN_PATH_LEN;
        if (neigh_fw.size() == 1) {
            len_fw++;
            if (len_fw <  MIN_PATH_LEN) {
                Kmer kmer_solido = node,
                        kmer_neigh = node;
                kmer_solido.appendRight(neigh_fw[0]);
                if (_is_standard)
                    kmer_solido.standard();
                kmer_neigh.appendRightReplace(neigh_fw[0]);
                k_vec.push_back(kmer_solido);
                k_node.push_back(kmer_neigh);
                _check_forward_path(len_fw, k_vec, k_node, kmer_neigh, show);
            }
        }
    }

    void _remove_isolated_nodes()
    {
        bool change = false;
        vector<Node> erase;
        unordered_set<Node> erase_nodes;
        for (auto kmer:_dbg_nodes)
        {
            size_t cont = 0;
            size_t in_nodes_fw = in_degree(kmer);
            size_t in_nodes_total = in_nodes_fw;
            size_t len = 1;
            /*
             * InDegree = 0
             */
            if (!in_nodes_total) {
                cont ++;
                if (!out_degree(kmer))
                {
                    erase_nodes.emplace(kmer);
                    continue;
                }
                vector<Node> aux = {};
                vector<Node> aux_nodes = {kmer};
                _check_forward_path(len,aux, aux_nodes, kmer);
                _asses(erase,erase_nodes,aux,aux_nodes,len);
            }
            /*
             * Check FWNeighbors and RCNeighbors (if proceeds)
             */
            vector<DnaSequence::NuclType> neighbors = getNeighbors(kmer);
            size_t num_neighbors = neighbors.size();
            if ( num_neighbors > 1)
            {
                /*
                 * Check branches
                 */

                for (auto nt:neighbors)
                {
                    Node k_aux = kmer, sibling = kmer;
                    k_aux.appendRight(nt);
                    if (_is_standard)
                        k_aux.standard();
                    sibling.appendRightReplace(nt);
                    len = 1;
                    vector<Node> aux = {k_aux};
                    vector<Node> aux_nodes = {sibling};
                    _check_forward_path(len,aux, aux_nodes, sibling);
                    _asses(erase,erase_nodes,aux,aux_nodes,len);
                }
            }
        }
        if (erase.size() > 0) {
            change = true;
            _erase(erase, erase_nodes);
        }
        /*
         * We have to iterate until convergence
         */
        if (change) {
            _remove_isolated_nodes();
        }else {
            cout << "Update DBG status:\n";
            _updateInfo();
            /*for (auto k:_in_0)
                cout << "KmerSuspicious: " << k.str() << "\n";*/
        }
        /*for (auto k:_dbg_nodes)
            cout << "KmerNodes: "<<k.str()<<"\n";
        for (auto k:_dbg_naive)
            cout << "KmerSolidos: "<<k.str()<<"\n";*/
    }

    void _updateInfo()
    {
        std::cout << "Total Solid k-Mers after cleaning: "<<_dbg_naive.size()
                                                          << " Total Graph Nodes: "<<_dbg_nodes.size()<<"\n";
    }

    void _printInfo()
    {
        std::cout << "Total Number of Bases: "<<_sc.getTotalBases()<<"\n";
        std::cout << "Average length read: "<<_sc.getAvLength()<<"\n";
        std::cout << "Total Kmers in all Reads: "<<_kmers_map.size()<<"\n";
        std::cout << "Threshold: "<<Parameters::get().accumulative_h<<"\n";
        std::cout<<"Total Solid K-mers(Graph Edges): "<<_dbg_naive.size()
                 <<" Total Graph Nodes: "<<_dbg_nodes.size()<<"\n";
    }

    void _buildGraphRepresentation(size_t max_freq)
    {
        vector<size_t> histogram = getHistogram<Node,size_t>(_kmers_map, max_freq);
        Parameters::get().accumulative_h = Parameters::calculateAccumulativeParam(histogram, _sc.getTotalBases(),_sc.getAvLength());
        cout << "Graph Representation\n";
        for (auto kmer:_kmers_map)
        {
            if (kmer.second.first >= Parameters::get().accumulative_h)
            {
                if (_is_standard)
                {
                    Node node = kmer.first;
                    node.standard();
                    _insert(node, node, kmer.second.first);
                }else{
                    Node node = kmer.first;
                    _insert(kmer.first, kmer.first, kmer.second.first);
                }
            }
        }
        cout << "END GRAPH\n";
        _printInfo();
        _kmers_map.clear();
    }
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
    void _insert(Node k, Node kmer, size_t freq)
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
     * AdHoc Methods
     */
    size_t _getScSize()
    {
        return (_is_standard)?_sc.size():_sc.size()/2;
    }
    /*
     * First Counter
     */
    unordered_map<Node, pair<size_t,size_t>> _kmers_map, _nodes_tmp_map;
    /*
     * PairedInfo
     */
    Extra<P> _extra_info;
    /*
     * DBG_naive -> stores the set of solid Kmers
     * DBG_nodes -> stores the set of (K-1)mers
     */
    unordered_set<Node> _dbg_naive, _dbg_nodes;
    unordered_map<Node, unordered_set<size_t>> _node_reads;
    unordered_set<KmerInfo<P>> _heads,_tails;
    //Extension points
    vector<Node> _in_0;
    //Extend
    SequenceContainer& _sc, _sc_paired;
    //Canonical Representation
    bool _is_standard = true;
};
/*
 * Boost implementation
 */
template <bool P> class boostDBG:public DBG<P>
{
public:
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
    typedef typename DBG<P>::Parent_Extra ExtraInfoNode;
    typedef typename DBG<P>::Parent_Paired_Info PairedInfoNode;
    typedef typename DBG<P>::Parent_Freq_Map FreqMap;
    /*
     * Boost structure
     */
    /*
     * PairEnd
     */
    struct NodeInfo {
        NodeInfo():id(-1){}
        NodeInfo(Node node, int32_t id):node(node), id(id)
        {
            coverage = 0;
            node_set = PairedInfoNode ();
            parent_cliques = map<Node, vector<PairedInfoNode>>();
        }
        NodeInfo(Node node, int32_t id, PairedInfoNode extra):node(node), id(id),node_set(extra){}
        NodeInfo(const NodeInfo & nodeInfo):node(nodeInfo.node),coverage(nodeInfo.coverage),id(nodeInfo.id),node_set(nodeInfo.node_set)
                ,parent_cliques(nodeInfo.parent_cliques){}
        NodeInfo(Node node, int32_t id, size_t ex_coverage):node(node),coverage(ex_coverage), id(id){
        }
        NodeInfo& operator=(const NodeInfo& other)
        {
            node = other.node;
            id = other.id;
            node_set = other.node_set;
            coverage = other.coverage;
            parent_cliques = other.parent_cliques;
            return *this;
        }
        bool empty()
        {
            return (id == -1);
        }
        bool equal(Node node) const
        {
            return node == this->node;
        }
        bool operator ==(const NodeInfo &other) const
        {
            return equal(other.node) && (this->id == other.id) && (node_set==other.node_set) && (coverage == other.coverage)
                   && (parent_cliques == other.parent_cliques);
        }
        Node node;
        size_t coverage;
        int32_t id;
        PairedInfoNode node_set;
        map<Node,vector<PairedInfoNode>> parent_cliques;
    };

    struct EdgeInfo {
        EdgeInfo():active(true){}
        EdgeInfo(size_t coverage): coverage(coverage),active(true){}
        EdgeInfo(int8_t nt_received, size_t id):nt(nt_received),id(id), active(true){}
        EdgeInfo& operator=(const EdgeInfo& other)
        {
            nt = other.nt;
            id = other.id;
            coverage = other.coverage;
            active = other.active;
            return *this;
        }
        int8_t nt;
        size_t id;
        size_t coverage;
        bool active;
    };
    struct GraphInfo{};
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, NodeInfo,
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
    typedef typename Graph::edge_iterator edge_iterator;
    typedef typename BUgraph<vertex_t>::graphBU graphBU;
    typedef unordered_set<graphBU> Extratmp;
    typedef pair<graphBU, PairedInfoNode> it_node;

    struct it_node_set{
        it_node_set(it_node node){
            _node.first = node.first;
            _node.second = node.second;
        }
        it_node_set(const it_node_set & other){
            _node.first = other._node.first;
            _node.second = other._node.second;
        }
        it_node_set(it_node_set && other){
            _node.first = std::move(other._node.first);
            _node.second = std::move(other._node.second);
        }
        it_node_set& operator=(const it_node_set & other)
        {
            _node.first = other._node.first;
            _node.second = other._node.second;
            return *this;
        }
        it_node_set& operator=(it_node_set && other)
        {
            _node = std::move(other._node);
            return *this;
        }
        bool operator==(const it_node_set &other) const
        {
            return (other._node.first == _node.first) && same(other._node.second, _node.second);
        }
        bool operator!=(const it_node_set &other) const
        {
            return (other._node.first != _node.first) || (!same(other._node.second, _node.second));
        }
        void show() const
        {
            cout << "Nodo: "<<_node.first<<" Size: "<<_node.second.size()<<endl;
            show_set(_node.second);
        }
        it_node _node;
    };
    struct hash_fn{
        size_t operator() (const it_node_set &other) const
        {
            size_t res = 0;
            boost::hash_combine(res, other._node.first);
            for (auto p: other._node.second)
                boost::hash_combine(res, p.hash());
            return res;
        }
    };
    /*
     * BoostDBG
     */
    boostDBG(DBG<P> *);
    boostDBG(string, string, SequenceContainer*);

    void clear()
    {
        _heads.clear();
        _tails.clear();
        _in_0.clear();
    }

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
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges(vertex, _g);
        for(; neighbors.first != neighbors.second; ++neighbors.first)
        {
            if (_g[*neighbors.first].active)
                neigh.push_back(_g[boost::target(*neighbors.first,_g)].node);
        }
        return neigh;
    }

    vector<graphBU> getKmerNeighbors(graphBU node) const
    {
        vector<graphBU> neigh = vector<graphBU>();
        pair<out_iterator, out_iterator> neighbors =
                boost::out_edges(node, _g);
        for (; neighbors.first != neighbors.second; ++neighbors.first)
        {
            if (_g[*neighbors.first].active)
                neigh.push_back(boost::target(*neighbors.first,_g));
        }
        return neigh;
    }

    vector<graphBU> getInKmerNeighbors(graphBU node) const
    {
        vector<graphBU> neigh = vector<graphBU>();
        pair<in_iterator,in_iterator> in_neighbors =
                boost::in_edges(node, _g);
        for (; in_neighbors.first != in_neighbors.second; ++in_neighbors.first) {
            if (_g[*in_neighbors.first].active)
                neigh.push_back(boost::source(*in_neighbors.first, _g));
        }
        return neigh;
    }

    pair<vector<it_node>, vector<it_node>> getOutKmerNeighbors(const it_node& node, bool full = true)
    {
        vector<it_node> neigh, neigh_alter;
        NodeInfo nodeInfo = _g[node.first], neighInfo;
        for (auto n:getKmerNeighbors(node.first))
        {
            neighInfo = _g[n];
            vector<ExtraInfoNode> rejected;
            vector<bool> append;
            if (neighInfo.parent_cliques[nodeInfo.node].empty() && node.second.empty())
            {
                neigh.push_back(pair<vertex_t, PairedInfoNode>(n, PairedInfoNode()));
                continue;
            }
            if (neighInfo.parent_cliques[nodeInfo.node].empty())
            {
                neigh_alter.push_back(pair<vertex_t, PairedInfoNode>(n, PairedInfoNode()));
                continue;
            }
            for (auto s:neighInfo.parent_cliques[nodeInfo.node])
            {
                if (node.second.empty()){
                    neigh.push_back(pair<vertex_t,PairedInfoNode>(n,s));
                    continue;
                }
                /*
                 * The parent node paired-end info has to be in every son haplotype -> Redundant
                 */
                if (isSame(getIntersection(node.second, nodeInfo.node_set),
                           getIntersection(s, nodeInfo.node_set), IDENTITY_RATIO)) {
                    neigh.push_back(pair<vertex_t, PairedInfoNode>(n, s));
                }
            }
        }
        if (neigh.empty() && full)
        {
            neigh =vector<it_node>(neigh_alter);
            neigh_alter.clear();
        }
        return pair<vector<it_node>,vector<it_node>>(neigh,neigh_alter);
    }

    vector<Node> getInKmerNeighbors(Node node) const
    {
        //TODO: Complete
        vector<Node> neigh;
        return neigh;
    }

    vector<it_node> getInNeighbors(const it_node& node, bool check_all)
    {
        vector<it_node> neigh, neigh_aux;
        unordered_map<Node,vertex_t> in_neighs;
        NodeInfo nodeInfo = _g[node.first];
        pair<in_iterator,in_iterator> in_neighbors =
                boost::in_edges(node.first, _g);
        for (; in_neighbors.first != in_neighbors.second; ++in_neighbors.first)
        {
            //in_neighs[_g[boost::source(*in_neighbors.first,_g)].node] = boost::source(*in_neighbors.first, _g);
            if (!_g[*in_neighbors.first].active)
                continue;
            graphBU endpoint_first = boost::source(*in_neighbors.first,_g);
            for (auto cliques:_g[endpoint_first].parent_cliques)
            {
                if (cliques.second.empty())
                {
                    neigh.push_back(pair<graphBU, PairedInfoNode>(endpoint_first, PairedInfoNode()));
                    if (!check_all)
                        return neigh;
                }
                for (auto clique:cliques.second){
                    if (isSame(getIntersection(clique,_g[endpoint_first].node_set),
                               getIntersection(node.second, _g[endpoint_first].node_set), IDENTITY_RATIO)) {
                        neigh.push_back(pair<graphBU, PairedInfoNode>(endpoint_first, clique));
                        if (!check_all)
                            return neigh;
                    }
                }
            }
        }
        if (check_all)
            return (neigh.size() > 0)?neigh:neigh_aux;
        else
            return neigh;
    }

    Extra<P> * getPairedInfo()
    {
        return nullptr;
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

    size_t in_degree(it_node node, bool check_all)
    {
        return getInNeighbors(node, check_all).size();
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

    size_t out_degree(it_node node)
    {
        return getKmerNeighbors(node).size();
    }

    void extension(vector<Node> in_0, string path_to_write);

    //Getter
    vector<Node> getEngagers()
    {
        vector<Node> in_0;
        for (auto n: _in_0)
            in_0.push_back(_g[n].node);
        return in_0;
    }

    typename DBG<P>::Heads get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }

    pair<unordered_set<Node> *, unordered_set<Node> *> getNodes()
    {
        return pair<unordered_set<Node> *, unordered_set<Node> * >(new unordered_set<Node>(), new unordered_set<Node>());
    }

    unordered_map<Node, unordered_set<size_t>> getNodeReads()
    {
        return _node_reads;
    };

    pair<bool, FreqMap> getExtra(Node node)
    {
        return {false, FreqMap()};
    };

    unordered_set<graphBU> getExtraInfoNode(graphBU node)
    {
        return _map_extra_info[_g[node].id];
    }

    PairedInfoNode getPairedInfoNode(graphBU node)
    {
        return _g[node].node_set;
    }

    graphBU getRepresentant(graphBU node)
    {
        return _representants_map[node];
    }

    size_t getRepresentantHitsNode(graphBU representant, graphBU node)
    {
        return _representants_hits[representant][node];
    }

    size_t getRepresentantTranslation(graphBU representant)
    {
        return translator_vector[_g[representant].id];
    }


    size_t getNumberCliques(graphBU endpoint,graphBU parent)
    {
        return _g[endpoint].parent_cliques[_g[parent].node].size();
    }

    vector<PairedInfoNode> getCliquesWithParent(graphBU endpoint,graphBU parent)
    {
        return _g[endpoint].parent_cliques[_g[parent].node];
    }

    void addHaplotype(graphBU endpoint, graphBU parent, PairedInfoNode haplotype)
    {
        _g[endpoint].parent_cliques[_g[parent].node].push_back(haplotype);
    }

    unordered_set<Node> getSolidKmers()
    {
        return unordered_set<Node>();
    }

    Node getNode(Node t)
    {
        return t;
    }

    Node getNode(graphBU k)
    {
        return _g[k].node;
    }

    SequenceContainer * getSequenceContainer()
    {
        return nullptr;
    }

    void ProcessTigs(string path_to_write)
    {
        cout << "Unitigs in APDB"<<endl;
        //UnitigExtender<P, graphBU>::full_extension(this, vector<graphBU>(), path_to_write);
        extension(vector<Node>(),path_to_write);
        cout << "Reported in : "<<path_to_write<<endl;
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

    void show_info(size_t);
private:
    void _show_internal_info()
    {
        cout << "Showing internal information: "<<endl;
        edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(_g); ei != ei_end; ++ei) {
            vertex_t  s = source(*ei, _g), t = target(*ei, _g);
            cout << "Source: "<<_g[s].node.str()<< " Coverage: "<<_g[s].coverage << endl;
            cout << "Target: "<<_g[t].node.str()<<" Coverage: "<<_g[t].coverage << endl;
            cout << "Edge Information (coverage): "<<_g[*ei].coverage <<endl;
        }
        std::cout << std::endl;
    }
    void _insert(Node kmer, size_t coverage)
    {
        graphBU origin_node, target_node;
        pair<Node,Node> suffix_preffix(kmer.substr(0, Parameters::get().kmerSize-1),kmer.substr(1, Parameters::get().kmerSize));
        pair<Node,Node> suffix_preffix_2 = kmer.preffixsuffix();
        typename unordered_map<Node, graphBU>::const_iterator origin = local_map.find(suffix_preffix.first),
                target = local_map.find(suffix_preffix.second);
        if (origin == local_map.end())
        {
            origin_node = boost::add_vertex(NodeInfo(suffix_preffix.first, _node_id++, coverage), _g);
            local_map[suffix_preffix.first] = origin_node;
        }else {
            origin_node = (*origin).second;
            _g[origin_node].coverage = std::max(_g[origin_node].coverage, coverage);
        }
        if (target == local_map.end())
        {
            target_node = boost::add_vertex(NodeInfo(suffix_preffix.second, _node_id++, coverage), _g);
            local_map[suffix_preffix.second] = target_node;
        }else {
            target_node = (*target).second;
            _g[target_node].coverage = std::max(_g[target_node].coverage, coverage);
        }
        boost::add_edge(origin_node, target_node, EdgeInfo(coverage), _g).first;
    }
    void _thirdPartyKmerCounting(string path_to_file_count, string dir_pairs, size_t * retries)
    {
        vector <string> files;
        if (path_to_file_count != "")
            for (auto f: System::getAllFaFqFiles(path_to_file_count))
                files.push_back(f);
        if (dir_pairs != "")
            for (auto f: System::getAllFaFqFiles(dir_pairs))
                files.push_back(f);
        if (files.size() > 1)
            path_to_file_count = System::appendFiles(files, "newFile_tmp.fasta");
        else if (files.size() == 1)
            path_to_file_count = files[0];
        string instruction = "";
        instruction += "bash -c \"./Utils/script/dsk_script ";
        instruction += path_to_file_count + " ";
        instruction += to_string(Parameters::get().kmerSize) + " ";
        instruction += "output output.txt >/dev/null 2>&1\"";
        cout << "INSTRUCTION: " << instruction << "\n";
        while ((system(instruction.c_str()))) {
            cout << "Problem executing: " << instruction << "\n";
            if ((*retries) > MAX_RETRIES)
                exit(1);
            cout << "Retrying with an alternative approach (first wait " << TIME_WAIT << "s)\n";
            sleep(TIME_WAIT);
            (*retries)++;
            _thirdPartyKmerCounting(path_to_file_count, dir_pairs, retries);
            return;
        }
        _directConstruction("output.txt");
    }
    void _checkStrain(graphBU node, size_t * edges_disconnected, bool forward = true, size_t coverage = 0)
    {
        vector<graphBU> neighbors = (forward)?getKmerNeighbors(node):getInKmerNeighbors(node);
        while ( neighbors.size() )
        {
            edges_disconnected++;
            for(auto n:(forward)?getKmerNeighbors(neighbors[0]):getInKmerNeighbors(neighbors[0]))
                if (_g[n].coverage < coverage*COV_RATIO)
                    neighbors.push_back(n);
            if (forward) {
                out_iterator ei, ei_end;
                for (tie(ei, ei_end) = out_edges(node, _g); ei != ei_end; ++ei)
                    _g[*ei].active = false;
            }else{
                in_iterator ei, ei_end;
                for (tie(ei, ei_end) = in_edges(node, _g); ei != ei_end; ++ei)
                    _g[*ei].active = false;
            }
        }
    }
    void _polishing()
    {
        size_t edges_disconnected = 0;
        cout << "Starting polish"<<endl;
        edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(_g); ei != ei_end; ++ei)
        {
            if (_g[*ei].active) {
                if (_g[source(*ei, _g)].coverage > (_g[*ei].coverage*COV_RATIO) ||
                    _g[target(*ei, _g)].coverage > (_g[*ei].coverage*COV_RATIO)) {
                    edges_disconnected++;
                    _g[*ei].active = false;
                    (_g[source(*ei, _g)].coverage > _g[*ei].coverage)?_checkStrain(target(*ei,_g),&edges_disconnected,true, _g[*ei].coverage)
                        :_checkStrain(source(*ei,_g),&edges_disconnected,false,_g[*ei].coverage);
                }
            }
        }
        cout << "Number of removals: "<<edges_disconnected<<endl;
    }
    void _directConstruction(string inFile)
    {
        int lineNo = 1, count = 0;
        Node y;
        std::ifstream infile(inFile);
        for( std::string line; getline( infile, line ); )
        {
            if (line.back() == '\n' or line.back() == '\r')
            {
                line.pop_back();
            }
            count = stoi(line.substr(Parameters::get().kmerSize+1,line.size()-(Parameters::get().kmerSize+1)));
            y = Node(line.substr(0, Parameters::get().kmerSize));
            _insert(y, count);
            //map[t] = pair<Y,Y>(count,count);
            lineNo++;
        }
        infile.close();
        cout << "End building"<<endl;
        cout << "Initial properties: "<<endl;
    }
    /*
     * Check Information
     */
    void _printInfo()
    {
        size_t num_vertex = boost::num_vertices(_g), num_edges = boost::num_edges(_g);
        cout << "Number of vertex: "<<num_vertex<<endl;
        cout << "Number of edges: "<<num_edges<<endl;
    }
    /*
     * Boost insert extra info
     */
    void _insert_extra(SequenceContainer*, unordered_set<Node>*, unordered_set<Node>*);
    /*
     * actives_haplotypes contains haplotypes in source actives (but only pairs of source)
     */
    vector<ExtraInfoNode> _getSharedHaplotypes(Node source, Node target, vector<ExtraInfoNode> active_haplotypes)
    {
        vector<ExtraInfoNode> newHaplotypes;
        for (auto k: active_haplotypes)
        {
            NodeInfo target_info = _g[*target];
            for (auto k2: target_info.parent_cliques[source])
            {
                if (isSubset(k,k2))
                {
                    newHaplotypes.push_back(getIntersection(k2,target_info.node_set));
                }
            }
        }
        return newHaplotypes;
    }
    vector<ExtraInfoNode> _getSharedHaplotypes(NodeInfo source, NodeInfo target, vector<ExtraInfoNode> active_haplotypes)
    {
        vector<ExtraInfoNode> newHaplotypes;
        for (auto k: active_haplotypes)
        {
            for (auto k2: target.parent_cliques[source.node])
            {
                if (isSubset(k,k2))
                {
                    newHaplotypes.push_back(getIntersection(k2,target.node_set));
                }
            }
        }
        return newHaplotypes;
    }
    void _transverse(const it_node &,
                     map<size_t, vector<size_t>> &,
                     map<size_t, graphBU> &,
                     map<size_t, DnaSequence> &,
                     unordered_set<graphBU> &,
                     DnaSequence &,
                     size_t &);
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
        return graphBU();
    }

    void _insertExtraInfo(unordered_map<Node, vertex_t> map, DBG<P> * dbg)
    {
        vertex_iterator v, vend;
        size_t num_vertex = boost::num_vertices(_g);
        _map_extra_info.resize(num_vertex);
        for(size_t i = 0; i < num_vertex; ++i)
        {
            for (auto s:_g[i].node_set){
                vertex_t pointerToExtraInfo = map[s];
                _map_extra_info[_g[i].id].emplace(pointerToExtraInfo);
            }
        }
    }

    void _writeUnitigs(map <size_t , vector<size_t>> graphUnitigs,
                                       map <size_t, graphBU > map_seq_nodo_end,
                                       map <size_t, DnaSequence> map_seqs,
                                       string file_to_path)
    {
        unordered_set<DnaSequence> seqs;
        vector<pair<size_t,size_t>> links;
        vector<bool> strand;
        for (auto s:map_seqs)
            seqs.emplace(s.second);
        for (auto p:graphUnitigs)
        {
            size_t from = p.first;
            for (auto to: p.second)
                links.push_back(pair<size_t, size_t>(to, from));
        }
        UnitigExtender<false>::_write_gfa(file_to_path, map_seqs,links);
        UnitigExtender<false>::_write_fasta(file_to_path+".fasta", seqs);
    }
    /*
     * Get Representants
     */
    void _get_representatives();
    size_t getVertex()
    {
        return boost::num_vertices(_g);
    }
    size_t getNumRepresentants()
    {
        return _first_last.size();
    }
    /*
     * Completing info
     */
    int* _floyds_warshall();
    void _full_fil_matrix(size_t, size_t, size_t, graphBU, graphBU, size_t, vector<bool> &, bool*, unordered_set<graphBU>&);
    bool _reachable(int*, size_t, size_t);
    bool _reachable(graphBU, graphBU , graphBU , size_t, size_t *, vector<bool> &, vector<bool> &, size_t);
    void _remove_outliers(Extra<P> *, DBG<P> *,unordered_set<Node> *);
    void _modify_info();
    //Graph
    Graph _g;
    vector<bool> reach;
    vector<size_t> translator_vector;
    unordered_map<Node, unordered_set<size_t>> _node_reads;
    //Representants
    vector<graphBU> _representants;
    unordered_set<graphBU> _first_last;
    unordered_map<graphBU, graphBU> _representants_map;
    unordered_map<graphBU, unordered_map<graphBU,size_t>> _representants_hits;
    unordered_map<graphBU, pair<graphBU, graphBU>> _represent_first_last;
    unordered_set<it_node_set, hash_fn> node_launched;
    unordered_set<graphBU> extended;

    //Properties + nodes (this should be fixed with vecS)
    unordered_map<Node, graphBU > local_map;
    Extra<P> _extra_info;
    vector<unordered_set<graphBU>> _map_extra_info;
    vector<map<graphBU, vector<Extratmp>>> _map_click_parents;
    vector<graphBU> _in_0;
    vector<it_node> _in_0_pairs;
    //Node_id
    int32_t _node_id = 0, _edge_id = 0;
    size_t seg = 0;
    //Need fix
    unordered_set<KmerInfo<P>> _heads,_tails;
};


//Adjacency list
template <bool P> class listDBG:public DBG<P>
{
public:
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
    typedef typename DBG<P>::Parent_Extra ExtraInfoNode;
    typedef typename DBG<P>::Parent_Paired_Info PairedInfoNode;
    typedef typename DBG<P>::Parent_Freq_Map FreqMap;
    typedef typename BUgraph<Node>::graphBU graphBU;
    typedef DnaSequence Unitig;
    typedef unordered_map<Node,pair<vector<Node>,vector<Node>>> Graph;
    listDBG(DBG<P> *);

    void clear()
    {
        _g.clear();
        _solid_kmers.clear();
        _heads.clear();
        _tails.clear();
        _in_0.clear();
    }

    size_t in_degree(Node node)
    {
        return _g[node].first.size();
    }

    size_t out_degree(Node node)
    {
        return _g[node].second.size();
    }
    vector<DnaSequence::NuclType> getNeighbors (Node node) const
    {
        return vector<DnaSequence::NuclType>();
    }

    vector<Node> getKmerNeighbors
            (Node node) const
    {
        if (node.length() == Parameters::get().kmerSize)
            node = node.substr(1, Parameters::get().kmerSize);
        return _g.at(node).second;
    }

    vector<Node> getInKmerNeighbors(Node node) const
    {
        if (node.length() == Parameters::get().kmerSize)
            node = node.substr(0,Parameters::get().kmerSize-1);
        return _g.at(node).first;
    }

    bool is_solid(Node & node) const
    {
        if (_is_standard) {
            Node node_aux = node;
            node_aux.standard();
            return (_solid_kmers.find(node_aux) != _solid_kmers.end());
        }else
            return (_solid_kmers.find(node) != _solid_kmers.end());
    }

    size_t length() const
    {
        return _g.size();
    }

    void extension(vector<Node> in_0, string path_to_write);

    pair<unordered_set<Node> * ,unordered_set<Node> *> getNodes()
    {
        return pair<unordered_set<Node> * , unordered_set<Node> * >(new unordered_set<Node>(), new unordered_set<Node>());
    }

    unordered_set<Node> getSolidKmers()
    {
        return _solid_kmers;
    }
    pair<bool,FreqMap> getExtra(Node node)
    {
        return pair<bool,FreqMap>(false,FreqMap());
    }
    unordered_map<Node, unordered_set<size_t>> getNodeReads()
    {
        return _node_reads;
    };
    vector<Node> getEngagers()
    {
        return _in_0;
    }

    Extra<P> * getPairedInfo()
    {
        return nullptr;
    }

    SequenceContainer * getSequenceContainer()
    {
        return nullptr;
    }

    void show_info()
    {
        for (auto n:_g)
        {
            cout << "Kmer: "<< n.first.str()<<"\n";
            for (auto in_:n.second.first)
                cout << "InNode: "<<in_.str()<<"\n";
            cout << "InDegree: "<<in_degree(n.first)<<"\n";
            for (auto out_:getKmerNeighbors(n.first))
                cout << "OutNode2: "<<out_.str()<<"\n";
            cout << "OutDegree: "<<out_degree(n.first)<<"\n";
        }
    }

    void ProcessTigs(string path_to_write)
    {
        cout << "Go for Unitigs\n";
        extension(_in_0, path_to_write);
        cout << "End Unitigs\n";
    }

    typename DBG<P>::Heads  get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }
private:
    void _transverse(Node,
                     map <Node, vector<size_t>> &,
                     map <size_t, Node> &,
                     map <size_t, DnaSequence> &,
                     DnaSequence &);
    void _writeUnitigs(map <Node, vector<size_t>> ,
                       map <size_t, Node> ,
                       map <size_t, DnaSequence>,
                       string);
    void _buildNewGraph(DBG<P> *);
    Graph _g;
    unordered_set<Node> _solid_kmers;
    vector<Node> _in_0;
    size_t seg = 0;
    //InfoHeads
    typename DBG<P>::Heads _heads,_tails;
    unordered_map<Node, unordered_set<size_t>> _node_reads;
    //Canonical Representation
    bool _is_standard = true;
};

template <bool P> class gatbDBG:public DBG<P>
{
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
    typedef typename DBG<P>::Parent_Extra ExtraInfoNode;
    typedef typename DBG<P>::Parent_Paired_Info PairedInfoNode;
    typedef typename DBG<P>::Parent_Freq_Map FreqMap;
    typedef typename BUgraph<Node>::graphBU graphBU;
    typedef DnaSequence Unitig;
    typedef unordered_map<Node,pair<vector<Node>,vector<Node>>> Graph;
    gatbDBG(string, string, SequenceContainer*);
    void clear()
    {
    }

    size_t in_degree(Node node)
    {
        return 0;
    }

    size_t out_degree(Node node)
    {
        return 0;
    }
    vector<DnaSequence::NuclType> getNeighbors (Node node) const
    {
        return vector<DnaSequence::NuclType>();
    }

    vector<Node> getKmerNeighbors
            (Node node) const
    {
        return vector<Node>();
    }

    vector<Node> getInKmerNeighbors(Node node) const
    {
        return vector<Node>();
    }

    bool is_solid(Node & node) const
    {
       return false;
    }

    size_t length() const
    {
        return 0;
    }

    void extension(vector<Node> in_0, string path_to_write);

    pair<unordered_set<Node> * ,unordered_set<Node> *> getNodes()
    {
        return pair<unordered_set<Node> * , unordered_set<Node> * >(new unordered_set<Node>(), new unordered_set<Node>());
    }

    unordered_set<Node> getSolidKmers()
    {
        return unordered_set<Node>();
    }
    pair<bool,FreqMap> getExtra(Node node)
    {
        return pair<bool,FreqMap>(false,FreqMap());
    }
    unordered_map<Node, unordered_set<size_t>> getNodeReads()
    {
        return _node_reads;
    };
    vector<Node> getEngagers()
    {
        return vector<Node>();
    }

    Extra<P> * getPairedInfo()
    {
        return nullptr;
    }

    SequenceContainer * getSequenceContainer()
    {
        return nullptr;
    }

    void show_info()
    {

    }

    void ProcessTigs(string path_to_write)
    {

    }

    typename DBG<P>::Heads  get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }
private:
    void _transverse(Node,
                     map <Node, vector<size_t>> &,
                     map <size_t, Node> &,
                     map <size_t, DnaSequence> &,
                     DnaSequence &);
    void _writeUnitigs(map <Node, vector<size_t>> ,
                       map <size_t, Node> ,
                       map <size_t, DnaSequence>,
                       string);
    void _buildNewGraph(DBG<P> *);
    Graph _g;
    unordered_set<Node> _solid_kmers;
    vector<Node> _in_0;
    size_t seg = 0;
    //InfoHeads
    typename DBG<P>::Heads _heads,_tails;
    unordered_map<Node, unordered_set<size_t>> _node_reads;
    //Canonical Representation
    bool _is_standard = true;
};