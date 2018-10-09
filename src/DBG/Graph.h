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
    typedef Kmer DBGFuncInp;
};
template<> struct NodeType<true>
{
    typedef Kmer DBGNode;
    typedef Pair_Kmer DBGFuncInp;
};

template <bool P> struct Extra{
    virtual void insert(typename NodeType<P>::DBGNode, typename NodeType<P>::DBGNode);
    virtual void show_info();
    virtual void clear();
    virtual bool find(typename NodeType<true>::DBGNode);
    virtual unordered_set<typename NodeType<P>::DBGNode> operator[] (typename NodeType<P>::DBGNode) const;
};
template <> struct Extra<false>
{
    void show_info() {}
    void clear() {}
    bool find(typename NodeType<true>::DBGNode){
        return true;
    }
};
template <> struct Extra<true>
{
    unordered_map<typename NodeType<true>::DBGNode, unordered_set<typename NodeType<true>::DBGNode>> mapPair;
    virtual unordered_set<typename NodeType<true>::DBGNode> operator[] (typename NodeType<true>::DBGNode node) const
    {
        return mapPair.at(node);
    }
    void clear()
    {
        mapPair.clear();
    }
    bool find(typename NodeType<true>::DBGNode node) const
    {
        return mapPair.find(node)!=mapPair.end();
    }
    void insert(typename NodeType<true>::DBGNode key, typename NodeType<true>::DBGNode val)
    {
        mapPair[key].emplace(val);
    }
    void show_info()
    {
        for (auto k:mapPair)
        {
            cout << k.first.str()<<":";
            for (auto k_:k.second)
                cout << k_.str() << ", ";
            cout << "\n";
        }
    }
};

namespace std
{
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
    typedef typename NodeType<P>::DBGNode Parent_Node;
    typedef typename NodeType<P>::DBGFuncInp Parent_FuncNode;
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
