/*
 * Path graph where:
 *  - Vertes are solid k-mers.
 *  - Edges are paths labelled with edit distance.
 * Implementation
 *  - Adjacent lists.
 * */
#include <unordered_map>

struct Edge{
    DnaSequence seq;
    size_t ed;
};

typedef Kmer Node;

class PathGrap{
public:
    PathGrap(){};
};

class PathGraphAdj: public PathGrap
{
public:
    PathGraphAdj(){};
private:
    std::unordered_map<Kmer,std::pair<std::vector<Node>,std::vector<Edge>>> _adj_list;
};