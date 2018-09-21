#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include "../Extender/Extender.h"


//Constants
#define MIN_PATH_LEN 1

class NaiveDBG: public DBG
{
public:
    NaiveDBG(SequenceContainer& sc):_sc(sc)
    {
        Progress::get().size_total = _sc.getIndex().size();
        Progress::get().show = true;
        _kmerCount();
        _cleaning();
    }
    bool is_solid(Kmer) const;
    vector<DnaSequence::NuclType> getNeighbors
            (const Kmer &) const;
    vector<Kmer> getKmerNeighbors
            (const Kmer &) const;
    size_t length() const
    {
        return _dbg_naive.size();
    }
    size_t in_degree(Kmer);
    size_t out_degree(Kmer);

    Heads get(bool behaviour) const
    {
        return (behaviour)?_heads:_tails;
    }

    void ProcessTigs(string path_to_write)
    {
        std::cout << "Lets start\n";
        _getTigs(path_to_write);
        std::cout << "End Unitigs\n";
    }

    void show_info();

    //Operators
    /*
     * Check if this is working properly
     */
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
    void _kmerCount()
    {
        Kmer kmer;
        KmerInfo tail;
        bool first = false;
        for (auto &read:_sc.getIndex()){
            if (first)
                _tails.emplace(tail);
            Progress::update(read.first.getId());
            first = false;
            for (auto kmer_r: IterKmers(read.second.sequence)) {
                kmer = kmer_r.kmer;
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
                    }
                } else
                    _kmers_map[kmer] = pair<size_t,size_t>(1,kmer_r.kmer_pos);
                if (Parameters::get().accumulative_h == 1)
                    if (_kmers_map[kmer].first == Parameters::get().accumulative_h)
                        _dbg_naive.emplace(kmer);

            }
        }
        _kmers_map.clear();
        /*for (auto &k: _dbg_naive)
            cout << "Kmer: "<<k.str()<<"\n";*/
        Progress::update(_sc.getIndex().size());
    }

    void _cleaning()
    {
        _remove_isolated_nodes();
    }

    void _getTigs(string path_to_unitigs)
    {
        unordered_set<Kmer> kmers_added;
        /*
         * Lists of kmers which have more than one out degree and in degree
         */
        vector<Kmer> out, in;
        vector<vector<Kmer>> unitigs;
        /*
         * Paths started because 0 in_degree
         */
        for (auto &k: _in_0)
            kmers_added.emplace(k);
        for (auto &k:_in_0)
        {
            vector<vector<Kmer>> inter_unitigs = UnitigExtender::Extend(k,*this,out,in,kmers_added);
            for (auto &vect:inter_unitigs)
                unitigs.push_back(vect);
        }
        /*
         * Paths started because > 1 out_degree
         */
        for (auto &k:out)
        {
            vector<vector<Kmer>> inter_unitigs = UnitigExtender::Extend(k,*this,out,in,kmers_added);
            for (auto &vect:inter_unitigs)
                unitigs.push_back(vect);
        }
        /*
         * Paths started because > 1 in_degree
         */
        for (auto &k:in)
        {
            vector<vector<Kmer>> inter_unitigs = UnitigExtender::Extend(k,*this,out,in,kmers_added);
            for (auto &vect:inter_unitigs)
                unitigs.push_back(vect);
        }
        vector<DnaSequence> sequence_for_unitig = _get_sequences(unitigs);
        _write_unitigs(sequence_for_unitig, path_to_unitigs);

    }

    vector<DnaSequence> _get_sequences(vector<vector<Kmer>> unitigs)
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

    void _check_path(size_t& len, vector<Kmer>& k_vec) const
    {
        Kmer aux = k_vec.back();
        std::vector<DnaSequence::NuclType> neigh = getNeighbors(aux);
        if (neigh.size() == 1) {
            len++;
            if (len < MIN_PATH_LEN){
                Kmer kmer_aux = aux;
                kmer_aux.appendRight(neigh[0]);
                k_vec.push_back(kmer_aux);
                _check_path(len,k_vec);
            }
        }if (neigh.size() > 1)
            len += MIN_PATH_LEN;
        else
            len += 0;
    };

    bool _asses(vector<Kmer> &erase,vector<Kmer> aux, size_t len)
    {
        if (len < MIN_PATH_LEN)
            for (auto k:aux)
                erase.push_back(k);
        return len < MIN_PATH_LEN;
    }

    void _erase(vector<Kmer>& kmer_to_erase)
    {
        for (auto kmer_erase:kmer_to_erase)
            _dbg_naive.erase(kmer_erase);
        kmer_to_erase.clear();
    }

    void _remove_isolated_nodes()
    {
        bool change = false, in_0_erase = true;
        vector<Kmer> erase;
        for (auto kmer:_dbg_naive) {
            size_t cont = 0;
            size_t len = 1, in_nodes = in_degree(kmer);
            if (!in_nodes) {
                cont ++;
                /*
                 * Check isolated nodes
                 */
                vector<Kmer> aux;
                aux.push_back(kmer);
                _check_path(len,aux);
                in_0_erase = _asses(erase,aux,len);
            }
            if (!in_0_erase)
                _in_0.push_back(kmer);
            in_0_erase = true;
            //Implementar algo tipo cache para evitar pasar varias veces por un kmer
            vector<Kmer> neighbors = getKmerNeighbors(kmer);
            size_t cont_fake_branches = 0;
            if ( neighbors.size() > 1)
            {
                /*
                 * Check branches
                 */
                for (auto sibling:neighbors)
                {
                    len = 1;
                    vector<Kmer> aux;
                    aux.push_back(sibling);
                    _check_path(len,aux);
                    if (_asses(erase,aux,len))
                        cont_fake_branches++;
                }
            }
        }
        if (erase.size() > 0)
            change = true;
        _erase(erase);
        /*
         * We have to iterate until convergence
         */
        if (change)
            _remove_isolated_nodes();
        cout << "KmerSolids: "<<_dbg_naive.size() << "; Suspicious Starts: "<<_in_0.size()<< "\n";
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

    unordered_map<Kmer, pair<size_t,size_t>> _kmers_map;
    //_dbg_naive graph, set of first solid k-mers
    unordered_set<Kmer> _dbg_naive;
    unordered_set<KmerInfo> _heads,_tails;
    //Extension points
    vector<Kmer> _in_0, _out_0, _in_1;

    //Extend
    SequenceContainer& _sc;
};