#include "PathGraph.h"

using namespace std;
//PathGraphAdj

void PathGraphAdj::add_edge(const Kmer &source,const Kmer &target, size_t edit, DnaSequence path)
{
    size_t kmer_size = Parameters::get().kmerSize;
    if (path.length() == kmer_size)
        std::cout << "Path: "<<path.str() << "\n";
    unordered_map<Node,AdjList>::const_iterator it = _adj_list.find(source);
    if (it == _adj_list.end()) {
        _adj_list[source] = AdjList({target}, {Edge((path.length() < kmer_size)?
                                                    path :path.substr(kmer_size,path.length()-kmer_size),edit)});
    }else{
        _adj_list[source].first.push_back(target);
        _adj_list[source].second.push_back(Edge((path.length() < kmer_size)?path
                                                :path.substr(kmer_size,path.length()-kmer_size),edit));
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
    for (uint i = 0; i < nodes.size()-2; ++i)
    {
        AdjList adj = _adj_list[nodes[i]];
        vector<Node>::iterator pos = find(adj.first.begin(),adj.first.end(),nodes[i+1]);
        uint pos_int = distance(adj.first.begin(),pos);
        std::cout << nodes[i].str()<<"\n";
        if (adj.second[pos_int].ed == 0) {
            for (uint j = 0; j < adj.second[pos_int].seq.length(); ++j) {
                optimal_read.append_nuc_right(adj.second[pos_int].seq.atRaw(j));
                std::cout <<"1" <<adj.second[pos_int].seq.at(j);
            }
            std::cout << "Arriba: "<<optimal_read.str() << "\n";
        }else
        {

            for (uint j = 0; j < Parameters::get().kmerSize; ++j) {
                optimal_read.append_nuc_right(nodes[i].at(j));
            }
            for (uint j = 0; j < adj.second[pos_int].seq.length(); ++j) {
                optimal_read.append_nuc_right(adj.second[pos_int].seq.atRaw(j));
                std::cout << adj.second[pos_int].seq.at(j);
            }
            std::cout << "\n";
            std::cout << optimal_read.str() << "\n";
        }
    }
    std::cout << "A dormir "<<optimal_read.str() <<"\n";
    sleep(10000);
    return optimal_read;
}

/*Dijkstra approach with priority queue
 * TODO: Use a better approach than re-build the path
 * */
DnaSequence PathGraphAdj::shortest_path(const Kmer &source, const Kmer &target){
    //First lets check the nodes degree
    vector<Kmer> optimal_path;
    priority_queue<pair<int,vector<Kmer>>, vector<pair<int,vector<Kmer>>>, CompareDist> prior_q;
    if (check_isolated())
    {
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
                return build_optimal_read(optimal_path);
            }
        }
    }

    return build_optimal_read(optimal_path);
}