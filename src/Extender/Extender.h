#include "../DBG/Graph.h"
#include <stack>
/*
 * UnitigExtender depends on Graph type and NodeType.
 */
using namespace std;

template<bool P>
class UnitigExtender
{
public:
    //TODO: Esto hay que mirarlo
    static vector<vector<Kmer>> Extend(Kmer,DBG<P>&,stack<Kmer>&,stack<Kmer>&,unordered_set<Kmer>);
    static void full_extension(DBG<P>&, vector<Kmer>, string);

private:
    /*
    * Info unitigs
    */
    static void _construct_sequences(vector<vector<Kmer>> unitigs)
    {
        for (auto & vect:unitigs) {
            size_t cont = 0;
            DnaSequence seq_local;
            if (_fin_segs.find(vect[0]) == _fin_segs.end())
                seq_local = DnaSequence(vect[0].str());
            for (auto k: vect)
                if (cont++)
                    seq_local.append_nuc_right(k.at(Parameters::get().kmerSize-1));
            _seqs.push_back(seq_local);
        }
        /*
         * Lets check the results
         */
        /*for (auto seq : _seqs)
            cout << "Seqs: "<<seq.str() << "\n";*/
    }
    static void _write_gfa(string filename)
    {
        FILE* fout = fopen(filename.c_str(), "w");
        if (!fout)
            throw std::runtime_error("Can't open " + filename);
        //size_t num_unitig = 0;
        string header = "H\tVN:Z:1\n";
        fwrite(header.data(), sizeof(header.data()[0])
                ,header.size(), fout);
        int num_seq = 0;
        /*
         * Write seqs -> gfa
         */
        for (auto& seq : _seqs)
        {
            string s_line = "S\t"+to_string(num_seq++)+"\t";
            s_line+= seq.str()+"\n";
            fwrite(s_line.data(), sizeof(s_line.data()[0]),
                   s_line.size(), fout);
        }
        /*
         * Write Links
         */
        for (auto & link:_links)
        {
            string l_line = "L\t"+to_string(link.second)+"\t"+to_string(link.first)+"\t*\n";
            fwrite(l_line.data(), sizeof(l_line.data()[0])
                    ,l_line.size(), fout);
        }
    }

    static size_t _curr_segment;
    static unordered_map<Kmer, vector<size_t>> _fin_segs;
    static vector<pair<size_t,size_t>> _links;
    //Sequences
    static vector<DnaSequence> _seqs;
};