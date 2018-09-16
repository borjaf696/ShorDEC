#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include "../Extender/Extender.h"


//Constants
#define MIN_PATH_LEN 5

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

    void ProcessTigs()
    {
        _getTigs();
    }

    void show_info();


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
            }
        }
        _kmers_map.clear();
        Progress::update(_sc.getIndex().size());
    }

    void _cleaning()
    {
        _remove_isolated_nodes();
        std::cout << "Sizes: "<< _in_0.size() << "-"<<_in_1.size() << "-"<<_out_0.size() << "\n";
    }

    void _getTigs()
    {}

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
        bool in_erase = false, already_in = false;
        vector<Kmer> erase;
        for (auto kmer:_dbg_naive) {
            size_t len = 0, in_nodes = in_degree(kmer);
            if (!in_nodes) {
                /*
                 * Check isolated nodes
                 */
                vector<Kmer> aux;
                aux.push_back(kmer);
                _check_path(len,aux);
                in_erase = _asses(erase,aux,len);
            }
            /*
             * Path origin: in_degree == 0 and path > min_path
             */
            if (!in_erase && !in_nodes) {
                already_in = true;
                _in_0.push_back(kmer);
            }
            if (in_nodes > 1){
                already_in = true;
                _in_1.push_back(kmer);
            }
            //Implementar algo tipo cache para evitar pasar varias veces por un kmer
            vector<Kmer> neighbors = getKmerNeighbors(kmer);
            size_t cont_fake_branches = 0;
            if ( neighbors.size() > 1)
            {
                std::cout << "Neighbors \n";
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
            if (cont_fake_branches == neighbors.size())
                erase.push_back(kmer);
            else if ((neighbors.size() - cont_fake_branches > 1) && !already_in)
                /*
                 * Path origin: out_degree > 0 and not inserted yet
                 */
                _out_0.push_back(kmer);
        }
        _erase(erase);
        cout << "KmerSolids: "<<_dbg_naive.size() << "\n";
    }

    unordered_map<Kmer, pair<size_t,size_t>> _kmers_map;

    //_dbg_naive graph, set of first solid k-mers
    unordered_set<Kmer> _dbg_naive;
    unordered_set<Kmer> _extenders;
    unordered_set<KmerInfo> _heads,_tails;

    //Extension points
    vector<Kmer> _in_0, _out_0, _in_1;

    //Extend
    const SequenceContainer& _sc;
};