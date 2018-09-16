#include <unistd.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include "../ReadData/kmer.h"
#include "../Utils/utils.h"

using namespace std;

typedef unordered_set<KmerInfo> Heads;
struct Node_ext
{
    Kmer kmer;
    size_t _in = 0, _out = 0;
};

class DBG
{
public:
    DBG(){}
    virtual bool is_solid(Kmer kmer) const = 0;
    virtual size_t length() const = 0;
    virtual vector<DnaSequence::NuclType> getNeighbors
            (const Kmer &) const = 0;
    virtual vector<Kmer> getKmerNeighbors
            (const Kmer &) const = 0;
    virtual size_t in_degree(Kmer) = 0;
    virtual size_t out_degree(Kmer) = 0;
    virtual Heads  get(bool) const = 0;
    virtual void ProcessTigs() = 0;
    //Show methods
    virtual void show_info() = 0;
private:
    virtual void _kmerCount() = 0;
    virtual void _cleaning() = 0;
    virtual void _getTigs() = 0;

    //All DBG can handle uni/omnitigs
    vector<DnaSequence> _tigs;
};
