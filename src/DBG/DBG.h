#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <unistd.h>
#include "../ReadData/kmer.h"
#include "../Utils/utils.h"

//Constants
#define MIN_PATH_LEN 5

using namespace std;
class DGB;

class DGB
{
public:
    DGB(){}
    virtual bool is_solid(Kmer kmer) const = 0;
    virtual size_t length() const = 0;
    virtual std::vector<DnaSequence::NuclType> getNeighbors
            (const Kmer &) const = 0;
    virtual size_t in_degree(Kmer) = 0;
    virtual size_t out_degree(Kmer) = 0;
private:
    virtual void _kmerCount() = 0;
};

class NaiveDGB: public DGB
{
public:
    NaiveDGB(SequenceContainer& sc):_sc(sc){
        Progress::get().size_total = _sc.getIndex().size();
        Progress::get().show = true;
        _kmerCount();
        //TODO: Better cleaning approach. Remove short branches (Leenas idea)
        _remove_isolated_nodes();
    }

    bool is_solid(Kmer) const;
    void check_path(Kmer,size_t&) const;

    std::vector<DnaSequence::NuclType> getNeighbors
            (const Kmer &) const;

    size_t length() const{
        return _dbg_naive.size();
    }

    size_t in_degree(Kmer);
    size_t out_degree(Kmer);

private:
    void _kmerCount()
    {
        Kmer kmer;
        for (auto &read:_sc.getIndex()){
            Progress::update(read.first.getId());
            for (auto kmer_r: IterKmers(read.second.sequence)) {
                kmer = kmer_r.kmer;//kmer_r.kmer;
                std::unordered_map<Kmer, size_t>::const_iterator place =
                        _kmers_map.find(kmer);
                if (place != _kmers_map.end()) {
                    _kmers_map[kmer]++;
                    if (_kmers_map[kmer] == Parameters::get().accumulative_h)
                        _dbg_naive.emplace(kmer);
                } else
                    _kmers_map[kmer] = 1;
            }
        }
        Progress::update(_sc.getIndex().size());
    }

    void _remove_isolated_nodes()
    {
        std::vector<Kmer> erase;
        for (auto kmer:_dbg_naive) {
            size_t in = in_degree(kmer);
            size_t len = 0;
            if (!in) {
                check_path(kmer,len);
                if ( len < MIN_PATH_LEN)
                    erase.push_back(kmer);
            }
        }
        for (auto kmer:erase)
            _dbg_naive.erase(kmer);
        std::cout << _dbg_naive.size() << "\n";
    }
    std::unordered_map<Kmer, size_t> _kmers_map;
    std::unordered_set<Kmer> _dbg_naive;
    const SequenceContainer& _sc;
};