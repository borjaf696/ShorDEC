#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include "../Extender/Extender.h"


//Constants
#define MIN_PATH_LEN 10

template<bool P>
class NaiveDBG: public DBG<P>
{
public:
    typedef typename DBG<P>::Parent_Node Node;
    typedef typename DBG<P>::Parent_FuncNode FuncNode;
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
     * Get the Nt in the edges
     */
    vector<DnaSequence::NuclType> getNeighbors
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
     * "Length" of the DBG
     */
    size_t length() const
    {
        return _dbg_naive.size();
    }
    /*
     * Number of in_edges
     */
    size_t in_degree(Node k)
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
    size_t out_degree(Node k)
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

    typename DBG<P>::Heads get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }

    void ProcessTigs(string path_to_write)
    {
        std::cout << "Lets start\n";
        UnitigExtender<P>::full_extension(*this,_in_0,path_to_write);
        std::cout << "End Unitigs\n";
    }

    void show_info();

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

    void _cleaning()
    {
        _remove_isolated_nodes();
    }

    /*
     * Insertion into the graph_nodes and graph_edges
     */
    void _insert(Node, FuncNode, bool = false);

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

    //TODO: Revisar todo lo asociado con los Standard, pensar en hacer 2 instancias 1 para fw y otra para rc
    void _remove_isolated_nodes()
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
     */
    unordered_set<Node> _dbg_naive, _dbg_nodes;
    unordered_set<KmerInfo<P>> _heads,_tails;
    //Extension points
    vector<Node> _in_0;
    //Extend
    SequenceContainer& _sc;
    //Standard
    bool _is_standard = true;
};