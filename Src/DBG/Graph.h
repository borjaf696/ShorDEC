#include <unistd.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <boost/graph/adjacency_list.hpp>
#include "../ReadData/kmer.h"
#include "../Utils/utils.h"

using namespace std;
template <bool> struct NodeType;
template<> struct NodeType<false>
{
    typedef Kmer DBGNode;
    typedef Kmer DBGFuncInp;
    //TODO: Corregir
    typedef int set_couples;
};
template<> struct NodeType<true>
{
    typedef Kmer DBGNode;
    typedef Pair_Kmer DBGFuncInp;
    typedef unordered_set<DBGNode> set_couples;
};

template<typename T>
struct BUgraph{
    typedef T graphBU;
};

template <bool P> struct Extra{
    virtual void insert(typename NodeType<P>::DBGNode, typename NodeType<P>::DBGNode);
    virtual void show_info();
    virtual void clear();
    virtual bool find(typename NodeType<true>::DBGNode);
    virtual void erase(typename NodeType<true>::DBGNode);
    virtual unordered_set<typename NodeType<P>::DBGNode> operator[] (typename NodeType<P>::DBGNode) const;
    virtual pair<bool,typename NodeType<false>::set_couples> getInfoNode(typename NodeType<P>::DBGNode);
};
template <> struct Extra<false>
{
    void show_info() {}
    void clear() {}
    bool find(typename NodeType<false>::DBGNode){
        return true;
    }
    pair<bool,typename NodeType<false>::set_couples> getInfoNode(typename NodeType<false>::DBGNode)
    {
        return pair<bool, typename NodeType<false>::set_couples >(false,0);
    }
};
/*
 * TODO: Optimize
 */
template <> struct Extra<true>
{
    unordered_map<typename NodeType<true>::DBGNode, NodeType<true>::set_couples> mapPair;
    unordered_set<typename NodeType<true>::DBGNode> key_and_value;
    pair<bool,typename NodeType<true>::set_couples> getInfoNode(typename NodeType<true>::DBGNode node)
    {
        if (find(node))
            return pair<bool,typename NodeType<true>::set_couples>(true,mapPair[node]);
        else
            return pair<bool,typename NodeType<true>::set_couples>(false,NodeType<true>::set_couples());
    }
    void erase(typename NodeType<true>::DBGNode node)
    {
        unordered_map<typename NodeType<true>::DBGNode, NodeType<true>::set_couples>::iterator map_iterator =
                mapPair.find(node);
        if (map_iterator != mapPair.end())
            mapPair.erase(node);
        /*
         * How to know places where a key is stored as value? -> This is painful
         */
        if (key_and_value.find(node) != key_and_value.end())
        {
            for (auto & pk: mapPair)
            {
                if (pk.second.find(node) != pk.second.end())
                    pk.second.erase(node);
            }
        }
    }
    unordered_set<typename NodeType<true>::DBGNode> operator[] (typename NodeType<true>::DBGNode node) const
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
        if (mapPair.find(val) != mapPair.end())
            key_and_value.emplace(val);
    }
    void show_info()
    {
        for (auto k:mapPair)
        {
            cout <<"\t" <<k.first.str()<<":";
            for (auto k_:k.second)
                cout << k_.str() << " ";
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
    typedef typename NodeType<P>::set_couples Parent_Extra;

    DBG(){}
    virtual void clear() = 0;
    virtual bool is_solid(Parent_Node&) const = 0;
    virtual size_t length() const = 0;
    virtual vector<typename DnaSequence::NuclType> getNeighbors
            (Parent_Node) const = 0;
    virtual size_t out_degree(Parent_Node) = 0;
    //Getter
    virtual Heads  get(bool) const = 0;
    virtual pair<unordered_set<Parent_Node>,
            unordered_set<Parent_Node>> getNodes() = 0;
    //Todo: think better (or not)
    virtual pair<bool,typename NodeType<P>::set_couples> getExtra(Parent_Node) = 0;
    virtual vector<Parent_Node> getEngagers() = 0;
    /*
     * This is like a hot fix but however->polymorphism :(
     */
    virtual unordered_set<Parent_Node> getSolidKmers() = 0;
    virtual Parent_Node getNode(Parent_Node k)
    {return k;}
    virtual void insert(Parent_Node node){};
    virtual size_t in_degree(Parent_Node) = 0;

    virtual vector<Parent_Node> getKmerNeighbors
            (Parent_Node) const = 0;
    //Show methods
    virtual void show_info() = 0;
    //Extension -> Todo: Move
    virtual void ProcessTigs(string) = 0;
    vector<vector<Parent_Node>> extend(Parent_Node kmer,
                                   stack<Parent_Node> &out,
                                   stack<Parent_Node> &in,
                                   unordered_set<Parent_Node> added,
                                   size_t &curr_segment,
                                   unordered_map<Parent_Node, vector<size_t>> &fin_segs)
    {
        vector<vector<Parent_Node>> unitigs;
        Parent_Node k_node = getNode(kmer);
        vector<Parent_Node> neighbors = getKmerNeighbors(kmer);
        for (auto &k: neighbors)
        {
            vector<Parent_Node> unitig;
            unitig.push_back(k_node);
            pair<size_t,Parent_Node> result =  _Extension(k,unitig,out, in, added);
            if (result.first == 1 || result.first == 2)
            {
                fin_segs[getNode(result.second)].push_back(curr_segment);
            }
            curr_segment++;
            unitigs.push_back(unitig);
        }
        return unitigs;
    }
    virtual void extension(vector<Parent_Node> in_0, string path_to_write) = 0;
private:

    pair<size_t, Parent_Node> _Extension(Parent_Node kmer,vector<Parent_Node> & unitig,
                                  stack<Parent_Node> & out, stack<Parent_Node> & in,
                                  unordered_set<Parent_Node> added)
    {
        vector<Parent_Node> neighbors = getKmerNeighbors(kmer);
        size_t in_ = in_degree(kmer);
        if (in_ == 1 && neighbors.size() == 1) {
            unitig.push_back(kmer);
            return _Extension(neighbors[0],unitig,out,in,added);
        }else if (neighbors.size () == 0) {
            unitig.push_back(kmer);
            return pair<size_t, Parent_Node>(0, kmer);
        }else if (added.find(kmer) == added.end())
        {
            added.emplace(kmer);
            if (neighbors.size() > 1) {
                out.push(kmer);
                unitig.push_back(kmer);
                return pair<size_t,Parent_Node>(1,kmer);
            }
            in.push(kmer);
            unitig.push_back(kmer);
            return pair<size_t, Parent_Node>(2,kmer);
        }
        return pair<size_t, Parent_Node>(0,kmer);
    }

    //All DBG can handle uni/omnitigs
    vector<DnaSequence> _tigs;
};
