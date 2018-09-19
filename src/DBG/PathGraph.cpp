#include "PathGraph.h"

using namespace std;
//PathGraphAdj

void PathGraphAdj::add_edge(const Node &source,const Node &target, size_t edit, DnaSequence path)
{
    unordered_map<Node,AdjList>::const_iterator it = _adj_list.find(source);
    if (it == _adj_list.end()) {
        _adj_list[source] = AdjList({target}, {Edge(path,edit)});
    }else{
        _adj_list[source].first.push_back(target);
        _adj_list[source].second.push_back(Edge(path,edit));
    }
}

size_t PathGraphAdj::num_vertex()
{
    return _adj_list.size();
}

size_t PathGraphAdj::num_edges()
{
    size_t edges_count = 0;
    for (auto vertex:_adj_list)
        edges_count += vertex.second.first.size();
    return edges_count;
}

//TODO:Remove this method
DnaSequence PathGraphAdj::build_optimal_read(vector<Node> nodes)
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
        vector<Node>::iterator pos = find(adj.first.begin(),adj.first.end(),nodes[i+1]);
        pos_int = distance(adj.first.begin(),pos);
        if (adj.second[pos_int].ed == 0) {
            for (uint j = 0; j < adj.second[pos_int].seq.length(); ++j) {
                optimal_read.append_nuc_right(adj.second[pos_int].seq.atRaw(j));
            }
        }else
        {
            for (uint j = 0; j < Parameters::get().kmerSize; ++j) {
                optimal_read.append_nuc_right(nodes[i].kmer.at(j));
            }
            for (uint j = 0; j < adj.second[pos_int].seq.length(); ++j) {
                optimal_read.append_nuc_right(adj.second[pos_int].seq.atRaw(j));
            }
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
DnaSequence PathGraphAdj::shortest_path(const Node &source, const Node &target){
    //First lets check the nodes degree
    if (source == target)
        return source.kmer.getSeq();
    vector<Node> optimal_path;
    priority_queue<pair<int,vector<Node>>, vector<pair<int,vector<Node>>>, CompareDist> prior_q;
    if (check_isolated())
    {
        cout << ":(\n";
        exit(1);
    }
    prior_q.push(pair<int,vector<Node>>(0,{source}));
    while (!prior_q.empty())
    {
        //Extract from the queue the top element
        pair<int,vector<Node>> head = prior_q.top();
        prior_q.pop();

        /*cout << "Execute! "<<head.second.back().str()<<"\n";
        cout << "Size: "<<prior_q.size()<<" Empty: "<<prior_q.empty()<<"\n";*/

        vector<pair<int,vector<Node>>> to_queue = push_in_queue(head);
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
bool PathGraphAdj::covered(const Node &kmer)
{
    return !(_adj_list.end() == _adj_list.find(kmer));
}

/*
 * Show methods
 */
void PathGraphAdj::show()
{
    for (auto k_list: _adj_list)
    {
        std::cout << " Kmer: "<<k_list.first.kmer.str() << " Pos: "<< k_list.first.kmer_pos << "\n";
    }
}