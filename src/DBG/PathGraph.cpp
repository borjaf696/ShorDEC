#include "PathGraph.h"

using namespace std;
//PathGraphAdj

void PathGraphAdj::add_edge(const Kmer &source,const Kmer &target, size_t edit, DnaSequence path)
{
    unordered_map<Node,AdjList>::const_iterator it = _adj_list.find(source);
    if (it == _adj_list.end()) {
        _adj_list[source] = AdjList({target}, {Edge(path, edit)});
    }else{
        _adj_list[source].first.push_back(target);
        _adj_list[source].second.push_back(Edge((!path.length())?path
                                                :path.substr(Parameters::get().kmerSize,path.length()),edit));
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

/*Dijkstra approach with priority queue*/
vector<Kmer> PathGraphAdj::shortest_path(const Kmer &source, const Kmer &target){
    //First lets check the nodes degree
    vector<Kmer> optimal_path;
    priority_queue<pair<int,vector<Kmer>>, vector<pair<int,vector<Kmer>>>, CompareDist> prior_q;
    if (check_isolated()) {
        cout << ":(\n";
        exit(1);
    }
    prior_q.push(pair<int,vector<Kmer>>(0,{source}));
    while (!prior_q.empty())
    {
        //Extract from the queue the top element
        pair<int,vector<Kmer>> head = prior_q.top();
        prior_q.pop();

        /*cout << "Execute! "<<head.second.back().str()<<"\n";
        cout << "Size: "<<prior_q.size()<<" Empty: "<<prior_q.empty()<<"\n";*/

        vector<pair<int,vector<Kmer>>> to_queue = push_in_queue(head);
        for (auto i:to_queue)
        {
            if (i.second.back() != target) {
                //cout <<"Kmer working over: " <<i.second.back().str() << "\n";
                prior_q.push(i);
            }else {
                optimal_path = i.second;
                return optimal_path;
            }
        }
    }

    return optimal_path;
}