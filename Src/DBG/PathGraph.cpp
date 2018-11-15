#include "PathGraph.h"

using namespace std;
//PathGraphAdj -> Single_end(false)
template<>
void PathGraphAdj<false>::add_edge(const typename PathGrap<false>::Node &source,
                                  const typename PathGrap<false>::Node &target,
                                  size_t edit, DnaSequence path)
{
    //std::cout << "Edge->From: "<<source.str()<<" to "<<target.str()<<"\n";
    unordered_map<Node,AdjList>::const_iterator it = _adj_list.find(source);
    if (it == _adj_list.end()) {
        _adj_list[source] = AdjList({target}, {Edge(path,edit)});
    }else{
        _adj_list[source].first.push_back(target);
        _adj_list[source].second.push_back(Edge(path,edit));
    }
}

template<>
size_t PathGraphAdj<false>::num_vertex()
{
    return _adj_list.size();
}

template<>
size_t PathGraphAdj<false>::num_edges()
{
    size_t edges_count = 0;
    for (auto vertex:_adj_list)
        edges_count += vertex.second.first.size();
    return edges_count;
}

//TODO:Remove this method
template<>
DnaSequence PathGraphAdj<false>::build_optimal_read(vector<typename PathGrap<false>::Node> nodes)
{
    /*
     * From optimal Kmers re-build the optimal Read
     * */
    DnaSequence optimal_read;
    AdjList adj;
    uint i = 0, pos_int = 0;
    for (i = 0; i < nodes.size()-1; ++i)
    {
        adj = _adj_list[nodes[i]];
        //std::cout << nodes[i].str()<<"\n";
        vector<typename PathGrap<false>::Node>::iterator pos = find(adj.first.begin(),adj.first.end(),nodes[i+1]);
        pos_int = distance(adj.first.begin(),pos);
        for (uint j = 0; j < adj.second[pos_int].seq.length(); ++j) {
            optimal_read.append_nuc_right(adj.second[pos_int].seq.atRaw(j));
        }
    }
    uint j;
    /*
     * If ed=0 the path already written belong to the first kmer
     */
    for (j=0;j < Parameters::get().kmerSize;++j) {
        optimal_read.append_nuc_right(nodes[i].kmer.at(j));
    }

    return optimal_read;
}

/*Dijkstra approach with priority queue
 * TODO: Use a better approach than re-build the path
 * */
template<>
DnaSequence PathGraphAdj<false>::shortest_path(const typename PathGrap<false>::Node &source,
                                               const typename PathGrap<false>::Node &target){
    //First lets check the nodes degree
    if (source == target)
        return source.kmer.getSeq();
    vector<typename PathGrap<false>::Node> optimal_path;
    priority_queue<pair<int,vector<typename PathGrap<false>::Node>>
            , vector<pair<int,vector<typename PathGrap<false>::Node>>>, CompareDist> prior_q;
    if (check_isolated())
    {
        cout << ":(\n";
        exit(1);
    }
    prior_q.push(pair<int,vector<typename PathGrap<false>::Node>>(0,{source}));
    while (!prior_q.empty())
    {
        //Extract from the queue the top element
        pair<int,vector<typename PathGrap<false>::Node>> head = prior_q.top();
        prior_q.pop();

        /*cout << "Execute! "<<head.second.back().str()<<"\n";
        cout << "Size: "<<prior_q.size()<<" Empty: "<<prior_q.empty()<<"\n";*/


        vector<pair<int,vector<typename PathGrap<false>::Node>>> to_queue = push_in_queue(head);
        for (auto i:to_queue)
        {
            if (i.second.back() != target) {
                //cout <<"Kmer working over: " <<i.second.back().str() << "\n";
                prior_q.push(i);
            }else {
                optimal_path = i.second;
                return build_optimal_read(optimal_path);
            }
        }
    }
    std::cout << "Oups\n";
    return build_optimal_read(optimal_path);
}

/*
 * Check if one solid kmer is being covered or not. Remember every kmer solid should be in the path graph
 * */
template<>
bool PathGraphAdj<false>::covered(const Node &kmer)
{
    return !(_adj_list.end() == _adj_list.find(kmer));
}

/*
 * Show methods
 */
template<>
void PathGraphAdj<false>::show()
{
    for (auto k_list: _adj_list)
    {
        std::cout << " Kmer: "<<k_list.first.kmer.str() << " Pos: "<< k_list.first.kmer_pos << "->\n";
        for (uint i = 0; i < k_list.second.first.size(); ++i){
            std::cout << k_list.second.first[i].kmer.str() << " "<<k_list.second.second[i].seq.str()
                      <<":"<<k_list.second.second[i].ed<<"\t";
        }
        std::cout << "\n";
    }
}


/*
 * Pathgraph -> Pair_end (true)
 */
template<>
void PathGraphAdj<true>::add_edge(const typename PathGrap<true>::Node &source,
                                   const typename PathGrap<true>::Node &target,
                                   size_t edit, DnaSequence path)
{

}

template<>
size_t PathGraphAdj<true>::num_vertex()
{
    return _adj_list.size();
}

template<>
size_t PathGraphAdj<true>::num_edges()
{
    size_t edges_count = 0;
    return edges_count;
}

//TODO:Remove this method
template<>
DnaSequence PathGraphAdj<true>::build_optimal_read(vector<typename PathGrap<true>::Node> nodes)
{
    /*
     * From optimal Kmers re-build the optimal Read
     * */
    DnaSequence optimal_read;

    return optimal_read;
}

/*Dijkstra approach with priority queue
 * TODO: Use a better approach than re-build the path
 * */
template<>
DnaSequence PathGraphAdj<true>::shortest_path(const typename PathGrap<true>::Node &source,
                                               const typename PathGrap<true>::Node &target){
    vector<typename PathGrap<true>::Node> optimal_path;
    return build_optimal_read(optimal_path);
}

/*
 * Check if one solid kmer is being covered or not. Remember every kmer solid should be in the path graph
 * */
template<>
bool PathGraphAdj<true>::covered(const Node &kmer)
{
    return !(_adj_list.end() == _adj_list.find(kmer));
}

/*
 * Show methods
 */
template<>
void PathGraphAdj<true>::show()
{
}
