#include "Extender.h"

/*
 * SingleEnd -> false
 */
template<bool P>
size_t UnitigExtender<P>::_curr_segment = 0;
template<bool P>
unordered_map <Kmer, vector<size_t>> UnitigExtender<P>::_fin_segs;
template<bool P >
vector <pair<size_t, size_t>> UnitigExtender<P>::_links;
template<bool P>
vector <DnaSequence> UnitigExtender<P>::_seqs;

pair<size_t, Kmer> Extension(Kmer kmer, DBG<false> &dbg, vector<Kmer> & unitig, stack<Kmer> & out, stack<Kmer> & in,
                             unordered_set<Kmer> added)
{
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    size_t in_ = dbg.in_degree(kmer);
    if (in_ == 1 && neighbors.size() == 1) {
        unitig.push_back(kmer);
        return Extension(neighbors[0],dbg,unitig,out,in,added);
    }else if (neighbors.size () == 0) {
        unitig.push_back(kmer);
        return pair<size_t, Kmer>(0, kmer);
    }else if (added.find(kmer) == added.end())
    {
        added.emplace(kmer);
        if (neighbors.size() > 1) {
            out.push(kmer);
            unitig.push_back(kmer);
            return pair<size_t,Kmer>(1,kmer);
        }
        in.push(kmer);
        unitig.push_back(kmer);
        return pair<size_t, Kmer>(2,kmer);
    }
    return pair<size_t, Kmer>(0,kmer);
}
template<>
vector<vector<Kmer>> UnitigExtender<false>::Extend(Kmer kmer, DBG<false> &dbg, stack<Kmer> & out, stack<Kmer> & in,
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
        {
            _fin_segs[result.second].push_back(_curr_segment);
        }
        _curr_segment++;
        unitigs.push_back(unitig);
    }
    return unitigs;
}


template<>
void UnitigExtender<false>::full_extension(DBG<false> & dbg, vector <Kmer> in_0, string path_to_write)
{
    /*
     * Set of already assesed heads
     */
    unordered_set<Kmer> added;
    /*
     * New heads with out and in > 1
     */
    stack<Kmer> in, out;
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
    while (!in.empty() && !out.empty())
    {
        /*
         * Lets check kmers with out_degree > 1
         */
        while (!out.empty())
        {
            Kmer k = out.top();
            out.pop();
            if (_fin_segs.find(k) != _fin_segs.end())
                for (uint i = 0; i < _fin_segs[k].size(); i++) {
                    _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[k][i]));
                }
            for (auto &p: Extend(k,dbg,out,in,added))
                unitigs.push_back(p);
        }
        /*
         * Lets check kmers with in_degree > 1
         */
        while(!in.empty())
        {
            Kmer k = in.top();
            in.pop();
            if (_fin_segs.find(k) != _fin_segs.end())
                for (uint i = 0; i < _fin_segs[k].size(); i++) {
                    _links.push_back(pair<size_t, size_t>(_curr_segment, _fin_segs[k][i]));
                }
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

/*
 * Pair_end
 */
/*
 * SingleEnd -> false
 */

pair<size_t, Kmer> Extension(Kmer kmer, DBG<true> &dbg, vector<Kmer> & unitig, stack<Kmer> & out, stack<Kmer> & in,
                             unordered_set<Kmer> added)
{
    vector<Kmer> neighbors = dbg.getKmerNeighbors(kmer);
    return pair<size_t, Kmer>(0,Kmer());
}
template<>
vector<vector<Kmer>> UnitigExtender<true>::Extend(Kmer kmer, DBG<true> &dbg, stack<Kmer> & out, stack<Kmer> & in,
                                                   unordered_set<Kmer> added)
{
    vector<vector<Kmer>> unitigs;
    return unitigs;
}

template<>
void UnitigExtender<true>::full_extension(DBG<true> & dbg, vector <Kmer> in_0, string path_to_write)
{

}
