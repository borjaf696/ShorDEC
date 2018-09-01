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
    virtual void add_edge(const Kmer&,const Kmer &,size_t,DnaSequence) = 0;

    //Transverse methods
    virtual std::vector<Kmer> shortest_path(const Kmer&, const Kmer&) = 0;

    //Check methods
    virtual size_t num_vertex() = 0;
    virtual size_t num_edges() = 0;
};

class PathGraphAdj: public PathGrap
{
public:
    PathGraphAdj(){};

    void add_edge(const Kmer&,const Kmer&, size_t = 0, DnaSequence = DnaSequence());
    std::vector<Kmer> shortest_path(const Kmer&, const Kmer&);

    //Operations over the graph
    size_t num_vertex();
    size_t num_edges();
    bool check_isolated()
    {
        for (auto k:_adj_list) {
            std::cout << k.first.str()<<"1\n";
            if (is_isolated(k.first))
                return true;
            std::cout << k.first.str()<<"2\n";
        }
        return false;
    }
private:
    typedef std::pair<std::vector<Node>,std::vector<Edge>> AdjList;
    class CompareDist
    {
    public:
        bool operator()(std::pair<int,std::vector<Kmer>> p1, std::pair<int,std::vector<Kmer>> p2)
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

    std::vector<std::pair<int,std::vector<Kmer>>> push_in_queue(std::pair<int,std::vector<Kmer>> q_el)
    {
        //List of all possible outcomes
        std::vector<std::pair<int,std::vector<Kmer>>> output;

        //Evaluated Kmer
        Kmer last_path_kmer = q_el.second.back();
        //Check all neighbors availables and add distances.
        int distance = q_el.first;
        //std::cout << last_path_kmer.str()<<"\n";
        for (uint i = 0; i < _adj_list[last_path_kmer].first.size(); ++i)
        {
            //Use the same instance I think is crazy
            std::vector<Kmer> k_list = q_el.second;
            Kmer kmer_aux;
            //If the edge has been already visited just skip it
            if (!_adj_list[last_path_kmer].second[i].vis) {
                kmer_aux = _adj_list[last_path_kmer].first[i];
                _adj_list[last_path_kmer].second[i].vis = true;
            }else
                continue;
            k_list.push_back(kmer_aux);
            //Add the elements which are going to go to the prior_queue!!!
            output.push_back(std::pair<int,std::vector<Kmer>>(distance+_adj_list[last_path_kmer].second[i].ed,
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
    void add_edge(const Kmer& , const Kmer&, size_t = 0, DnaSequence = DnaSequence());
    size_t num_vertex();
    size_t num_edges();
private:
    lemon::ListDigraph g;
};