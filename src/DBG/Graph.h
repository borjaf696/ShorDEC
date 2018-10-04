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
template<> struct NodeType<false>
{
    typedef Kmer DBGNode;
};
template<> struct NodeType<true>
{
    typedef Kmer DBGNode;
};

template <bool> struct ExtraType;
template<> struct ExtraType<false>{};
template<> struct ExtraType<true>
{
    Pair_Kmer pair;
};

namespace std
{
    template <>
    struct hash<ExtraType<true>>{
        size_t operator()(const ExtraType<true>& extraType) const {
            return extraType.pair.hash();
        }
    };

    template <>
    struct hash<ExtraType<false>>{
        size_t operator()(const ExtraType<false>& extraType) const {
            return 0;
        }
    };

    template <>
    struct hash<NodeType<true>>{
        size_t operator()(const NodeType<true>& nodeType) const {
            return 0;
        }
    };

    template <>
    struct hash<NodeType<false>>{
        size_t operator()(const NodeType<false>& nodeType) const {
            return 0;
        }
    };
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
            (typename NodeType<P>::DBGNode) const = 0;
    virtual vector<typename NodeType<P>::DBGNode> getKmerNeighbors
            (typename NodeType<P>::DBGNode) const = 0;
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
