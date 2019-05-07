#include <unistd.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include "../ReadData/kmer.h"

using namespace std;
template<typename T>
struct BUgraph{
    typedef T graphBU;
};
template <bool> struct NodeType;
template<> struct NodeType<false>
{
    typedef Kmer DBGNode;
    typedef Kmer DBGFuncInp;
    //TODO: Corregir
    typedef int * key;
    typedef unordered_set<key> value;
    typedef int static_paired_info;
    typedef int frequency_map;
};
template<> struct NodeType<true>
{
    typedef Kmer DBGNode;
    typedef Pair_Kmer DBGFuncInp;
    struct iterator_hash
    {
        size_t operator()(unordered_set<DBGNode>::iterator it) const
        {
            return (*it).hash();
        }
    };
    typedef unordered_set<DBGNode>::iterator key;
    typedef unordered_set<key, iterator_hash> value;
    typedef unordered_set<DBGNode> static_paired_info;
    typedef unordered_map<key,size_t,iterator_hash> frequency_map;
};

template <bool P> struct Extra{
    virtual void insert(typename NodeType<P>::key , typename NodeType<P>::key, size_t &);
    virtual void clear();
    virtual bool find(typename NodeType<true>::key );
    virtual void erase(typename NodeType<true>::key );
    virtual typename NodeType<P>::value operator[] (typename NodeType<P>::key ) const;
    virtual pair<bool,typename NodeType<P>::frequency_map> getInfoNode(typename NodeType<P>::key );
};
template <> struct Extra<false>
{
    void clear() {}
    bool find(typename NodeType<false>::key ){
        return true;
    }
    pair<bool,typename NodeType<false>::frequency_map> getInfoNode(typename NodeType<false>::key )
    {
        return pair<bool, typename NodeType<false>::frequency_map>(false,typename NodeType<false>::frequency_map());
    }
};
/*
 * TODO: Optimize
 */
template <> struct Extra<true>
{
    unordered_map<NodeType<true>::key, NodeType<true>::value, NodeType<true>::iterator_hash> mapPair;
    unordered_map<NodeType<true>::key, NodeType<true>::frequency_map, NodeType<true>::iterator_hash> freq_map;
    unordered_map<NodeType<true>::DBGNode, NodeType<true>::key> node_to_key;

    pair<bool,typename NodeType<true>::frequency_map> getInfoNode(typename NodeType<true>::key node)
    {
        if (find(node)) {
            return pair<bool, typename NodeType<true>::frequency_map>(true, freq_map[node]);
        }else
            return pair<bool, typename NodeType<true>::frequency_map>(false,NodeType<true>::frequency_map());
    }
    void show(size_t num_shows)
    {
        size_t cont = 0;
        for (auto p:mapPair){
            cout << "Kmer: "<<(*p.first).str()<<endl;
            for (auto s: p.second)
                cout << " "<<(*s).str();
            cout << endl;
            if (cont++ == num_shows)
                return;
        }
    }
    typename NodeType<true>::value getPairs(typename NodeType<true>::DBGNode node)
    {
        return mapPair[node_to_key[node]];
    }
    typename NodeType<true>::value operator[] (typename NodeType<true>::key node) const
    {
        return mapPair.at(node);
    }
    void clear()
    {
        mapPair.clear();
    }
    bool find(typename NodeType<true>::key node) const
    {
        return mapPair.find(node)!=mapPair.end();
    }
    bool in(typename NodeType<true>::key key, typename NodeType<true>::key val)
    {
        return mapPair[key].find(val) != mapPair[key].end();
    }
    void mod_freqs(typename NodeType<true>::DBGNode node, typename NodeType<true>::frequency_map freqs_change)
    {
        freq_map[node_to_key[node]] = freqs_change;
        /*for (auto val: freqs_change)
            (freq_map[node_to_key[node]])[val.first] = val.second;*/
    }
    void polish()
    {
        vector<size_t> thresholds;
        for (auto k: freq_map)
        {
            vector<size_t> histogram(1000,0);
            size_t max = 0;
            for (auto kk:k.second)
            {
                histogram[kk.second]++;
                max = (kk.second > max)?kk.second:max;
            }
            size_t threshold = floor(k.second.size()*1.0), sum = 0, i = 0;
            for (i = max; i > 0; --i)
            {
                sum += histogram[i];
                if (sum >= threshold)
                    break;
            }
            //thresholds.push_back(i);
            //Remove everything under threshold
            for (auto kk:k.second)
            {
                if (kk.second < i /*thresholds[it]*/)
                    (mapPair[k.first]).erase(kk.first);
            }
        }

        /*size_t it = 0;
        for (auto k: freq_map)
        {
            for (auto kk:k.second)
            {
                if (kk.second < thresholds[it])
                    (mapPair[k.first]).erase(kk.first);
            }
            it++;
        }*/
    };
    void show_distribution(string extra = "_pre_processing")
    {
        string histogram_file = "histogram_pairs"+extra;
        unordered_map<NodeType<true>::DBGNode, vector<size_t>> histogram_map;
        for (auto k: freq_map)
        {
            vector<size_t> histogram(1000,0);
            for (auto kk:k.second) {
                histogram[kk.second]++;
            }
            histogram_map[(*k.first)] = histogram;
        }
        System::write_histogram<NodeType<true>::DBGNode>(histogram_map, histogram_file);
    }
    void insert(typename NodeType<true>::key key, typename NodeType<true>::key val, size_t & adds)
    {
        if (mapPair.find(key)!=mapPair.end()) {
            pair<unordered_set<NodeType<true>::key, NodeType<true>::iterator_hash>::iterator, bool> insertion =  mapPair[key].emplace(val);
            if (insertion.second)
                adds++;
            if (freq_map[key].find(val) == freq_map[key].end())
                (freq_map[key])[val] = 1;
            else
                ((freq_map[key])[val])++;
        }else{
            mapPair[key] = {val};
            node_to_key[(*key)] = key;
            (freq_map[key])[val] = 1;
            adds++;
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
    typedef typename NodeType<P>::value Parent_Extra;
    typedef typename NodeType<P>::static_paired_info Parent_Paired_Info;
    typedef typename NodeType<P>::frequency_map Parent_Freq_Map;

    DBG(){}
    virtual void clear() = 0;
    virtual bool is_solid(Parent_Node&) const = 0;
    virtual size_t length() const = 0;
    virtual vector<typename DnaSequence::NuclType> getNeighbors
            (Parent_Node) const = 0;
    virtual size_t out_degree(Parent_Node) = 0;
    //Getter
    virtual Heads get(bool) const = 0;
    virtual pair<unordered_set<Parent_Node>*,
            unordered_set<Parent_Node>*> getNodes() = 0;
    //Todo: think better (or not)
    virtual pair<bool,Parent_Freq_Map> getExtra(Parent_Node) = 0;
    virtual SequenceContainer * getSequenceContainer() = 0;
    virtual Extra<P>* getPairedInfo() = 0;
    virtual vector<Parent_Node> getEngagers() = 0;
    virtual unordered_map<Parent_Node, unordered_set<size_t>> getNodeReads() = 0;
    /*
     * This is like a hot fix but however->polymorphism :(
     */
    virtual unordered_set<Parent_Node> getSolidKmers() = 0;
    virtual Parent_Node getNode(Parent_Node k)
    {return k;}
    virtual void insert(Parent_Node, size_t){};
    virtual size_t in_degree(Parent_Node) = 0;

    virtual vector<Parent_Node> getKmerNeighbors(Parent_Node) const = 0;
    virtual vector<Parent_Node> getInKmerNeighbors(Parent_Node) const = 0;
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
