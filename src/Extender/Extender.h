#include "../DBG/Graph.h"

using namespace std;

class UnitigExtender
{
public:
    //TODO: Esto hay que mirarlo
    static vector<vector<Kmer>> Extend(Kmer,DBG&,vector<Kmer>&,vector<Kmer>&,unordered_set<Kmer>);
    static void full_extension(DBG&, vector<Kmer>, string);

private:
    /*
    * Info unitigs
    */
    static void _construct_sequences(vector<vector<Kmer>>);
    static void _write_gfa(string);

    static size_t _curr_segment;
    static unordered_map<Kmer, vector<size_t>> _fin_segs;
    static vector<pair<size_t,size_t>> _links;
    //Sequences
    static vector<DnaSequence> _seqs;
};