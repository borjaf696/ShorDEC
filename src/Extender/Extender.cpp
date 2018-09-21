#include "Extender.h"


size_t UnitigExtender::_curr_segment = 0;
unordered_map<Kmer, vector<size_t>> UnitigExtender::_fin_segs;
vector<pair<size_t,size_t>> UnitigExtender::_links;
vector<DnaSequence> UnitigExtender::_seqs;

pair<size_t, Kmer> Extension(Kmer kmer, DBG &dbg, vector<Kmer> & unitig, vector<Kmer> & out, vector<Kmer> & in, unordered_set<Kmer> added)
{
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    size_t in_ = dbg.in_degree(kmer);
    if (in_ == 1 && neighbors.size() == 1) {
        unitig.push_back(kmer);
        Extension(neighbors[0],dbg,unitig,out,in,added);
    }else if (added.find(kmer) == added.end())
    {
        added.emplace(kmer);
        if (neighbors.size() > 1) {
            out.push_back(kmer);
            unitig.push_back(kmer);
            return {1,kmer};
        }
        in.push_back(kmer);
        unitig.push_back(kmer);
        return {2,kmer};
    }
    return {0,kmer};
}

vector<vector<Kmer>> UnitigExtender::Extend(Kmer kmer, DBG &dbg, vector<Kmer> & out, vector<Kmer> & in,
                            unordered_set<Kmer> added)
{
    vector<vector<Kmer>> unitigs;
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    for (auto &k: neighbors)
    {
        vector<Kmer> unitig;
        unitig.push_back(kmer);
        pair<size_t,Kmer> result = Extension(k, dbg, unitig, out, in, added);
        if (result.first == 1 || result.first == 2)
                _fin_segs[result.second].push_back(_curr_segment++);
        unitigs.push_back(unitig);
    }
    return unitigs;
}



void UnitigExtender::full_extension(DBG & dbg, vector <Kmer> in_0, string path_to_write)
{
    /*
     * Set of already assesed heads
     */
    unordered_set<Kmer> added;
    /*
     * New heads with out and in > 1
     */
    vector<Kmer> in, out;
    /*
     * Unitigs: entran por orden de curr_segment, ¡Nos ahorramos indexar las secuencias!
     */
    vector<vector<Kmer>> unitigs;
    for (auto k: in_0)
    {
        for (auto &p: Extend(k,dbg,out,in,added))
            unitigs.push_back(p);
    }
    //Remove the kmers from in_0
    in_0.clear();
    while (in.size() && out.size())
    {
        /*
         * Lets check kmers with out_degree > 1
         */
        for (auto k: out)
        {
            if (_fin_segs.find(k) != _fin_segs.end())
                for (uint i = 0; i < _fin_segs[k].size(); i++)
                    _links.push_back(pair<size_t,size_t>(_curr_segment,_fin_segs[k][i]));
            for (auto &p: Extend(k,dbg,out,in,added))
                unitigs.push_back(p);
        }
        /*
         * Lets check kmers with in_degree > 1
         */
        for (auto k: in)
        {
            if (_fin_segs.find(k) != _fin_segs.end())
                for (uint i = 0; i < _fin_segs[k].size(); i++)
                     _links.push_back(pair<size_t,size_t>(_curr_segment,_fin_segs[k][i]));
            for (auto &p: Extend(k,dbg,out,in,added))
                unitigs.push_back(p);
        }
    }
    /*
     * Construct the DnaSequences
     */
    _construct_sequences(unitigs);
    /*
     * Lets write them
     */
    _write_gfa(path_to_write);
}

void UnitigExtender::_construct_sequences(vector<vector<Kmer>> unitigs)
{
    for (auto & vect:unitigs) {
        size_t cont = 0;
        DnaSequence seq_local(vect[0].str());
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

void UnitigExtender::_write_gfa(string filename)
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
        string l_line = "L\t"+to_string(link.first)+"\t"+to_string(link.second)+"\t*\n";
        fwrite(l_line.data(), sizeof(l_line.data()[0])
                ,l_line.size(), fout);
    }
}