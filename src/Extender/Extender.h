#include "../DBG/Graph.h"

using namespace std;

class UnitigExtender
{
public:
    static vector<vector<Kmer>> Extend(Kmer,DBG&,vector<Kmer>&,vector<Kmer>&,unordered_set<Kmer>);
};