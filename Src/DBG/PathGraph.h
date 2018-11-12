#include <unordered_map>
#include <queue>
#include <vector>
#include <stdlib.h>
#include <unistd.h>
#include "lemon/list_graph.h"
#include "../DBG/DBG.h"
/*
 * Path graph where:
 *  - Vertes are solid k-mers.
 *  - Edges are paths labelled with the edit distance and the path between two k-mer solids.
 * Implementation
 *  - PathGraphAdj: Adjacent lists
 *  - PathGraphLemon: Lemon Implementation
 * */

/*Algoritmo de Dijstra para buscar el camino optimo.*/
struct Edge{
    Edge(DnaSequence seq1, size_t ed1)
            :seq(seq1),ed(ed1),vis(false) {}
    DnaSequence seq;
    size_t ed;
    //For check transversed
    bool vis;
};
template<bool P>
class PathGrap{
public:
    typedef KmerInfo<P> Node;
    PathGrap(){};

    //Constructor methods
    virtual void add_edge(const Node&,const Node &,size_t,DnaSequence) = 0;
    //Transverse methods
    virtual DnaSequence shortest_path(const Node&, const Node&) = 0;
    virtual bool covered(const Node&) = 0;
    //Show Methods
    virtual void show() = 0;
    //Re-build method
    virtual DnaSequence build_optimal_read(std::vector<Node>) = 0;
    //Check methods
    virtual size_t num_vertex() = 0;
    virtual size_t num_edges() = 0;
};

template <bool P>
class PathGraphAdj: public PathGrap<P>
{
public:
    PathGraphAdj(){};
    ~PathGraphAdj()
    {
        _adj_list.clear();
    }

    void add_edge(const typename PathGrap<P>::Node&,const typename PathGrap<P>::Node&, size_t = 0, DnaSequence = DnaSequence()) ;
    DnaSequence shortest_path(const typename PathGrap<P>::Node&, const typename PathGrap<P>::Node&);
    bool covered(const typename PathGrap<P>::Node&);
    void show();

    //TODO:Change to pvt
    DnaSequence build_optimal_read(std::vector<typename PathGrap<P>::Node>);

    //Operations over the graph
    size_t num_vertex();
    size_t num_edges();
    bool check_isolated()
    {
        for (auto k:_adj_list) {
            if (is_isolated(k.first))
                return true;
        }
        return false;
    }
private:
    typedef std::pair<std::vector<typename PathGrap<P>::Node>,std::vector<Edge>> AdjList;
    class CompareDist
    {
    public:
        bool operator()(std::pair<int,std::vector<typename PathGrap<P>::Node>> p1,
                        std::pair<int,std::vector<typename PathGrap<P>::Node>> p2)
        {
            return p1.first > p2.first;
        }
    };
    size_t out_degree(const typename PathGrap<P>::Node & k)
    {
        return _adj_list[k].first.size();
    }
    size_t in_degree(const typename PathGrap<P>::Node &k)
    {
        size_t out = 0;
        for (auto al:_adj_list)
        {
            typename std::vector<typename PathGrap<P>::Node>::iterator node_it = std::find(al.second.first.begin()
                    ,al.second.first.end(), k);
            if (node_it != al.second.first.end())
                out++;
        }
        return out;
    }

    bool is_isolated(const typename PathGrap<P>::Node &k)
    {
        return (!in_degree(k))&&(!out_degree(k));
    }

    std::vector<std::pair<int,std::vector<typename PathGrap<P>::Node>>> push_in_queue
            (std::pair<int,std::vector<typename PathGrap<P>::Node>> q_el)
    {
        //List of all possible outcomes
        std::vector<std::pair<int,std::vector<typename PathGrap<P>::Node>>> output;

        //Evaluated Node
        typename PathGrap<P>::Node last_path_node = q_el.second.back();
        //Check all neighbors availables and add distances.
        int distance = q_el.first;
        //std::cout << last_path_kmer.str()<<"\n";
        for (uint i = 0; i < _adj_list[last_path_node].first.size(); ++i)
        {
            /*std::cout << "Kmer Vecino: "<<_adj_list[last_path_node].first[i].str()<<"\n";
            std::cout << "PathTo: "<<_adj_list[last_path_node].second[i].ed<<"\n";*/
            //Use the same instance I think is crazy
            std::vector<typename PathGrap<P>::Node> k_list = q_el.second;
            //If the edge has been already visited just skip it
            if (!_adj_list[last_path_node].second[i].vis) {
                Kmer kmer_aux = _adj_list[last_path_node].first[i].kmer;
                k_list.push_back(typename PathGrap<P>::Node(kmer_aux,_adj_list[last_path_node].first[i].kmer_pos));
                _adj_list[last_path_node].second[i].vis = true;
            }else
                continue;
            //Add the elements which are going to go to the prior_queue!!!
            output.push_back(std::pair<int,std::vector<typename PathGrap<P>::Node>>(distance+_adj_list[last_path_node].second[i].ed,
                                       k_list));
        }
        /*for (uint i = 0; i < output.size(); ++i)
            std::cout << output[i].first<<" "<<output[i].second.back().str()<<"\n";*/
        return output;
    }
    std::unordered_map<typename PathGrap<P>::Node,AdjList> _adj_list;
};

template <bool P>
class PathGraphBoost: public PathGrap<P>
{
    typedef KmerInfo<P> Node;
    struct Edge {
        Edge(float weight, DnaSequence seq) : _weight(weight),_seq(seq) {}
        float _weight;
        DnaSequence _seq;
    };
    struct NodeGraph {
        NodeGraph():_id(-1){}
        NodeGraph (size_t id, Node node):_id(id),_node(node){}
        size_t _id;
        Node _node;
    };
    typedef typename boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS, NodeGraph, Edge> Graph;
    typedef typename boost::graph_traits <Graph>::vertex_descriptor vertex_t;

    PathGraphBoost(){};
    //Constructor methods
    void add_edge(const Node& source,const Node & target,size_t ed,DnaSequence seq)
    {
        vertex_t s, t;
        if (_map_node_vertex.find(source) != _map_node_vertex.end()){
            s = _map_node_vertex[source];
        }else{
            s = boost::add_vertex(NodeGraph(_curr_id++, source),_path_graph);
            _map_node_vertex[source] = s;
        }
        if (_map_node_vertex.find(target) != _map_node_vertex.end())
            t = _map_node_vertex[target];
        else{
            t = boost::add_vertex(NodeGraph(_curr_id++, target),_path_graph);
            _map_node_vertex[target] = t;
        }
        boost::add_edge(s,t,Edge(ed, seq),_path_graph);
    }
    //Transverse methods -> Dijkstra algorithm wont work for boost -> not faster at least
    DnaSequence shortest_path(const Node&, const Node&)
    {
        auto weights = boost::get(&Edge::_weight, _path_graph);
        return DnaSequence();
    }
    bool covered(const Node&)
    {
        return false;
    }
    //Show Methods
    void show()
    {

    }
    //Re-build method
    DnaSequence build_optimal_read(std::vector<Node>)
    {
        return DnaSequence();
    }
    //Check methods
    size_t num_vertex()
    {
        return 0;
    }
    size_t num_edges()
    {
        return 0;
    }
    std::unordered_map<Node, vertex_t> _map_node_vertex;
    size_t _curr_id = 0;
    Graph _path_graph;
};

template <bool P>
class PathGraphLemon: public PathGrap<P>
{
public:
    void add_edge(const typename PathGrap<P>::Node& ,
                  const typename PathGrap<P>::Node&, size_t = 0,
                  DnaSequence = DnaSequence());
    size_t num_vertex();
    size_t num_edges();
private:
    lemon::ListDigraph g;
};