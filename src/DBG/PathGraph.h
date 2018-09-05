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
typedef Kmer Node;

class PathGrap{
public:
    PathGrap(){};

    //Constructor methods
    virtual void add_edge(const Node&,const Node &,size_t,DnaSequence) = 0;

    //Transverse methods
    virtual DnaSequence shortest_path(const Node&, const Node&) = 0;
    virtual bool covered(const Node&) = 0;

    //Re-build method
    virtual DnaSequence build_optimal_read(std::vector<Node>) = 0;

    //Check methods
    virtual size_t num_vertex() = 0;
    virtual size_t num_edges() = 0;
};

class PathGraphAdj: public PathGrap
{
public:
    PathGraphAdj(){};

    void add_edge(const Node&,const Node&, size_t = 0, DnaSequence = DnaSequence()) ;
    DnaSequence shortest_path(const Node&, const Node&);
    bool covered(const Node&);

    //TODO:Change to pvt
    DnaSequence build_optimal_read(std::vector<Node>);

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
    typedef std::pair<std::vector<Node>,std::vector<Edge>> AdjList;
    class CompareDist
    {
    public:
        bool operator()(std::pair<int,std::vector<Node>> p1, std::pair<int,std::vector<Node>> p2)
        {
            return p1.first > p2.first;
        }
    };
    size_t out_degree(const Node & k)
    {
        return _adj_list[k].first.size();
    }
    size_t in_degree(const Node &k)
    {
        size_t out = 0;
        for (auto al:_adj_list)
        {
            std::vector<Node>::iterator node_it = std::find(al.second.first.begin()
                    ,al.second.first.end(), k);
            if (node_it != al.second.first.end())
                out++;
        }
        return out;
    }

    bool is_isolated(const Node &k)
    {
        return (!in_degree(k))&&(!out_degree(k));
    }

    std::vector<std::pair<int,std::vector<Node>>> push_in_queue(std::pair<int,std::vector<Node>> q_el)
    {
        //List of all possible outcomes
        std::vector<std::pair<int,std::vector<Node>>> output;

        //Evaluated Node
        Node last_path_node = q_el.second.back();
        //Check all neighbors availables and add distances.
        int distance = q_el.first;
        //std::cout << last_path_kmer.str()<<"\n";
        for (uint i = 0; i < _adj_list[last_path_node].first.size(); ++i)
        {
            //Use the same instance I think is crazy
            std::vector<Node> k_list = q_el.second;
            //If the edge has been already visited just skip it
            if (!_adj_list[last_path_node].second[i].vis) {
                Node node_aux;
                node_aux = _adj_list[last_path_node].first[i];
                _adj_list[last_path_node].second[i].vis = true;
                k_list.push_back(node_aux);
            }else
                continue;
            //Add the elements which are going to go to the prior_queue!!!
            output.push_back(std::pair<int,std::vector<Node>>(distance+_adj_list[last_path_node].second[i].ed,
                                       k_list));
        }
        /*for (uint i = 0; i < output.size(); ++i)
            std::cout << output[i].first<<" "<<output[i].second.back().str()<<"\n";*/
        return output;
    }
    std::unordered_map<Node,AdjList> _adj_list;
};

class PathGraphLemon: public PathGrap
{
public:
    void add_edge(const Node& , const Node&, size_t = 0, DnaSequence = DnaSequence());
    size_t num_vertex();
    size_t num_edges();
private:
    lemon::ListDigraph g;
};