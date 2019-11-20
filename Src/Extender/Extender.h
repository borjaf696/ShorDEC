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
    /*
     * Typedefs
     */
    typedef typename NodeType<P>::DBGNode Node;
    /*
     * Info unitigs -> Build at the same time we transverse the graph
     */
    static void _construct_sequences(vector<vector<Node>> unitigs,unordered_map<Node, vector<size_t>> _fin_segs
            , unordered_set<DnaSequence> & _seqs)
    {
        for (auto & vect:unitigs) {
            size_t cont = 0;
            DnaSequence seq_local;
            if (_fin_segs.find(vect[0]) == _fin_segs.end())
                seq_local = DnaSequence(vect[0].str());
            for (auto k: vect)
                if (cont++)
                    seq_local.append_nuc_right(k.at(Parameters::get().kmerSize-2));
            _seqs.emplace(seq_local);
        }
    }

    static void _write_gfa(string filename, unordered_set<DnaSequence> _seqs, vector<pair<size_t,size_t>> _links)
    {
        FILE* fout = fopen(filename.c_str(), "w");
        std::cout << "FileName: "<<filename<<"\n";
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
        unordered_set<DnaSequence> removedSequences;
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
            string l_line = "L\t"+to_string(link.second)+"\t"+"+\t"+to_string(link.first)+"\t+\t*\n";
            fwrite(l_line.data(), sizeof(l_line.data()[0])
                    ,l_line.size(), fout);
        }
        fclose(fout);
    }

    static void _write_gfa(string filename, map<size_t,DnaSequence> _seqs, vector<pair<size_t,size_t>> _links)
    {
        cout << "Exporting GFA"<<endl;
        unordered_set<size_t> removed;
        FILE* fout = fopen(filename.c_str(), "w");
        std::cout << "FileName: "<<filename<<"\n";
        if (!fout)
            throw std::runtime_error("Can't open " + filename);
        //size_t num_unitig = 0;
        string header = "H\tVN:Z:1\n";
        fwrite(header.data(), sizeof(header.data()[0])
                ,header.size(), fout);
        /*
         * Write seqs -> gfa
         */
        for (auto& mseq : _seqs)
        {
            DnaSequence seq = mseq.second;
            if (seq.length() < 500)
            {
                removed.emplace(mseq.first);
                continue;
            }
            string s_line = "S\t"+to_string(mseq.first)+"\t";
            s_line+= seq.str()+"\n";
            fwrite(s_line.data(), sizeof(s_line.data()[0]),
                   s_line.size(), fout);
        }
        /*
         * Write Links
         */
        for (auto & link:_links)
        {
            if (removed.find(link.first) != removed.end() || removed.find(link.second) != removed.end())
                continue;
            string l_line = "L\t"+to_string(link.second)+"\t"+"+\t"+to_string(link.first)+"\t+\t*\n";
            fwrite(l_line.data(), sizeof(l_line.data()[0])
                    ,l_line.size(), fout);
        }
        fclose(fout);
    }

    static void _write_fasta(string filename, unordered_set<DnaSequence> _seqs, size_t minLength = 500, std::string execFile = "Utils/script/PostProcessing/removeDup.py")
    {
        std::string execFile_alt = "/Utils/script/PostProcessing/removeDup.py"
        FILE* fout = fopen(filename.c_str(), "w");
        std::cout << "FileName: "<<filename<<endl;
        if (!fout)
            throw std::runtime_error("Can't open "+filename);
        int num_seq = 0;
        for (auto &seq: _seqs)
        {
            if (seq.length() < 500)
                continue;
            string s_line = ">"+std::to_string(num_seq)+"\n";
            s_line += seq.str()+"\n";
            fwrite(s_line.data(), sizeof(s_line.data()[0]), s_line.size(), fout);
            num_seq++;
        }
        fclose(fout);
        if (Parameters::get().postProcess)
        {
            std::cout << "Removing redundant unitigs"<<std::endl;
            if (System::execute("bash -c \"python "+execFile+" "+filename+"\""))
            {
                System::execute("bash -c \"python "+execFile_alt+" "+filename+"\"")
            }
            std::cout << "Done!"<<std::endl;
        }
    }
};