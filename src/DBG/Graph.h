#include <unistd.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include "../ReadData/kmer.h"
#include "../Utils/utils.h"

using namespace std;
template <bool> struct NodeType;
template<> struct NodeType<false>{
    typedef Kmer DBGNode;
};
template<> struct NodeType<true>{
    typedef Kmer DBGNode;
};

struct Node_ext
{
    Kmer kmer;
    size_t _in = 0, _out = 0;
};
template<bool P>
class DBG:public NodeType<P>
{
public:
    typedef unordered_set<KmerInfo<P>> Heads;
    DBG(){}
    virtual bool is_solid(typename NodeType<P>::DBGNode&) const = 0;
    virtual size_t length() const = 0;
    virtual vector<typename DnaSequence::NuclType> getNeighbors
            (const typename NodeType<P>::DBGNode &) const = 0;
    virtual vector<typename NodeType<P>::DBGNode> getKmerNeighbors
            (const typename NodeType<P>::DBGNode &) const = 0;
    virtual size_t in_degree(typename NodeType<P>::DBGNode) = 0;
    virtual size_t out_degree(typename NodeType<P>::DBGNode) = 0;
    virtual Heads  get(bool) const = 0;
    virtual void ProcessTigs(string) = 0;
    //Show methods
    virtual void show_info() = 0;
private:
    virtual void _kmerCount() = 0;
    virtual void _cleaning() = 0;

    //All DBG can handle uni/omnitigs
    vector<DnaSequence> _tigs;
};
