#include "path.h"

void insert_queue(std::queue<Kmer>& t, std::queue<size_t> &p_t, Kmer kmer, size_t pos)
{
    if (t.size() < NUM_PREV_POST_KMER) {
        t.push(kmer);
        p_t.push(pos);
    }else {
        t.pop();p_t.pop();
        t.push(kmer);p_t.push(pos);
    }
}

/*
 * Backtrack from the first solid k-mer until the beginning of the read.
 * TODO: Need test
 */
template<>
size_t Path<false>::extend_head(const DnaSequence &sub_sequence
        ,KmerInfo<false> t
        ,const DBG<false> &dbg
        ,size_t * score_ed
        ,size_t * max_pos
        ,char * expected_head
        ,size_t &branches
        ,bool behaviour
        ,KmerInfo<false> &which)
{
    size_t max_allowed_len = (t.kmer_pos*EXTRA_LEN);
    size_t best_len = MAX_PATH_LEN, count_heads = 0;
    char path[MAX_PATH_LEN+1];
    std::stack<stack_el*> neighbors, stack_swap;
    //For each possible head we are going to look for a path
    for (auto head_info: dbg.get(behaviour))
    {
        /*
         * If the kmers are the same -> skip
         * Topological info -> skip
         */
        if (head_info.kmer.str() == t.kmer.str())
            continue;
        if (head_info.kmer_pos > t.kmer_pos)
            continue;
        Kmer head = head_info.kmer;
        std::vector<DnaSequence::NuclType> nts
                = dbg.getNeighbors(head);
        for (uint i = 0; i < nts.size(); ++i)
        {
            Kmer kmer_aux = head;
            kmer_aux.appendRight(nts[i]);
            neighbors.push(new stack_el(kmer_aux,1,nts[i]));
        }

        size_t fail_len = t.kmer_pos;
        while (neighbors.size() > 0 && branches > 0)
        {
            stack_el * el_kmer = neighbors.top();
            neighbors.pop();
            Kmer cur_kmer = el_kmer->kmer;
            size_t pos = el_kmer->pos;
            DnaSequence::NuclType nuc = DnaSequence::nfi(el_kmer->nuc);

            //Deallocate the stack_el
            delete el_kmer;

            if (pos >= max_allowed_len)
                break;

            path[pos-1] = nuc;
            _DP[pos*(MAX_PATH_LEN+1)] = pos;
            size_t min = MAX_PATH_LEN;

            for (uint k = 1; k <= (fail_len); ++k)
            {
                size_t ev_pos = pos*(MAX_PATH_LEN+1)+k;
                (path[pos-1] == sub_sequence.at(k-1))?
                        _DP[ev_pos] = std::min(
                                _DP[ev_pos-(MAX_PATH_LEN+1)-1]
                                ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                                        ,_DP[ev_pos-1])+1):
                        _DP[ev_pos] = std::min(
                                _DP[ev_pos-(MAX_PATH_LEN+1)-1]+1
                                ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                                        ,_DP[ev_pos-1])+1);
                if (_DP[ev_pos] < min)
                    min = _DP[ev_pos];
            }
            /*
            * If edit_distance es min y además es la pos más alta aceptado
            */
            if (min < *(score_ed))
            {
                if (cur_kmer == t.kmer)
                {
                    branches--;
                    size_t edit_distance = _DP[pos * (MAX_PATH_LEN + 1) + fail_len];
                    if (edit_distance < *(score_ed) && pos > *(max_pos))
                    {
                        (*score_ed) = edit_distance;
                        (*max_pos) = pos;
                        path[pos] = '\0';
                        std::memcpy(expected_head, path, pos);
                        best_len = pos;
                        which = head_info;
                    }
                }else{
                    nts = dbg.getNeighbors(cur_kmer);
                    for (uint i = 0; i < nts.size(); i++){
                        Kmer kmer_aux = Kmer(cur_kmer.str());
                        kmer_aux.appendRight(nts[i]);
                        neighbors.push(new stack_el(Kmer(kmer_aux.str()),pos+1,nts[i]));
                    }
                    if (!nts.size())
                        branches--;
                }
            }else
                branches--;
        }
        for (uint i = 0; i > neighbors.size(); ++i)
        {
            stack_el * top_el = neighbors.top();
            neighbors.pop();
            delete top_el;
        }
        count_heads++;
        //std::cout << "NumOfNeighbors "<<neighbors.size()<< " "<<best_len<< "\n";
    }
    expected_head[best_len] = '\0';
    return best_len;
}

/*
 * Sub_sequence: part of the read to cover with the DBG
 * h: vector de posibles head-sources
 * t: vector de posibles tails
 * DBG: grafo de bruijn
 * score_ed: mejor score encontrado hasta el momento
 * length_cover: parte de la lectura cubierta
 * */
template<>
size_t Path<false>::extend(const DnaSequence &sub_sequence
        ,KmerInfo<false> h
        ,KmerInfo<false> t
        ,const DBG<false> &dbg
        ,size_t * score_ed
        ,char *expected_path
        ,size_t & branches)
{
    size_t best_len = MAX_PATH_LEN, kmerSize = Parameters::get().kmerSize;
    Kmer kmer_objective = Kmer(t.kmer.substr(0,kmerSize-1));
    unordered_set<Kmer> processed;
    char path[MAX_PATH_LEN+1];
    stack<stack_el*> neighbors;
    /*
     * Extract (k-1)mer from header to start the extension
     */
    vector<Kmer> nts
            = dbg.getKmerNeighbors(h.kmer);
    //Continuacion del kmer actual
    for (auto k: nts)
    {
        neighbors.push(new stack_el(k,1,k.at(kmerSize-2)));
        //Avoid repeat Kmers
        processed.emplace(k);
    }
    //Expected fail length (maximum fail length)
    size_t fail_len = t.kmer_pos - h.kmer_pos;
    while (neighbors.size() > 0 && branches > 0)
    {
        stack_el * el_kmer = neighbors.top();
        neighbors.pop();
        Kmer cur_kmer = el_kmer->kmer;
        size_t pos = el_kmer->pos;
        DnaSequence::NuclType nuc = DnaSequence::nfi(el_kmer->nuc);

        //Deallocate the stack_el
        delete el_kmer;

        //This can me improved->Why leave
        if (pos >= MAX_PATH_LEN)
            return MAX_PATH_LEN;

        path[pos-1] = nuc;
        _DP[pos*(MAX_PATH_LEN+1)] = pos;
        size_t min = MAX_PATH_LEN;
        /*
         * Different penalization over indels and subs
         * */
        for (uint k = 1; k <= (fail_len); ++k)
        {
            size_t ev_pos = pos*(MAX_PATH_LEN+1)+k;
            char nt = (sub_sequence.at(k-1)=='N')?path[pos-1]:sub_sequence.at(k-1);
            (path[pos-1] == nt)?
                    _DP[ev_pos] = std::min(
                            _DP[ev_pos-(MAX_PATH_LEN+1)-1]
                            ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                                    ,_DP[ev_pos-1])+1):
                    _DP[ev_pos] = std::min(
                    _DP[ev_pos-(MAX_PATH_LEN+1)-1]+1
                    ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                            ,_DP[ev_pos-1])+1);
            if (_DP[ev_pos] < min)
                min = _DP[ev_pos];
        }
        //Comparative between current path and alternative paths already visited
        if (min < (*score_ed)) {
            //Reach the expected node
            if (cur_kmer == kmer_objective) {
                branches--;
                size_t edit_distance = _DP[pos * (MAX_PATH_LEN + 1) + fail_len];
                if (edit_distance < (*score_ed)) {
                    (*score_ed) = edit_distance;
                    path[pos] = '\0';
                    memcpy(expected_path, path, pos+1);
                    best_len = pos;
                }
            } else {
                //We havent reached our objective we extend the path again
                nts = dbg.getKmerNeighbors(cur_kmer);
                for (auto k: nts)
                {
                    if (processed.find(k) == processed.end())
                    {
                        processed.emplace(k);
                        neighbors.push(new stack_el(k,pos+1,k.at(kmerSize-2)));
                    }
                }
                if (!nts.size())
                    branches--;
            }
        }else
            branches--;
    }
    while (neighbors.size())
    {
        stack_el * rem_el = neighbors.top();
        neighbors.pop();
        //Delete
        delete rem_el;
    }
    //std::cout << "NumOfNeighbors "<<neighbors.size()<< " "<<best_len<< "\n";
    return best_len;
}

//Container path
//TODO: Posiciones salvarlas y empezar la correccion
template<>
size_t PathContainer<false>::check_read()
{
    for (auto cur_kmer: IterKmers<false>(_seq)){
        //Optimizar esto
        /*Kmer kmer(cur_kmer.kmer.getSeq().substr(0
                ,cur_kmer.kmer.getSeq().length()));*/
        Kmer kmer = cur_kmer.kmer;
        if (_dbg.is_solid(kmer))
        {
            _solid.push_back(KmerInfo<false>(kmer, cur_kmer.kmer_pos));
        }
    }
    /*for (uint i = 0; i < _solid.size(); ++i)
        std::cout << _solid[i].kmer.str()<<" "<<_solid[i].kmer_pos<<"\n";*/
    return 1;
}

/*Mirar para cada Kmer solido todos los kmers siguientes a este susceptibles de trazar
 *un camino entre ellos:
 *  Si son adyacentes -> No se busca -> -1 trial
 *  Si estan solapados -> No se busca
 *  Si estan demasiado lejos no se busca
 * En cualquier otro caso se busca -> -1 trial
 * TODO: Check for more accuracy!
 * */
template<>
int PathContainer<false>::check_solids(size_t pos_i, size_t pos_j, size_t i,size_t j
        ,size_t &num_trials, Kmer& kmer)
{
    //Kmer are all adjacent.
    if ((pos_j-pos_i) == j-i) {
        num_trials++;
        return 3;
    }
    //Kmers far away from eachother
    if (((pos_j-pos_i)*(1.0+ERROR_RATE)) >= MAX_PATH_LEN) {
        num_trials = MAX_NUM_TRIALS;
        return 2;
    }
    //Both Kmers are overlapped.
    if ((pos_j - pos_i) < ((1-ERROR_RATE)*Parameters::get().kmerSize)) {
        return 1;
    }
    //They are neither too far nor too close -> True
    num_trials++;
    return 0;
}

//TODO: Correct tail->Optimize all process
template<>
DnaSequence PathContainer<false>::correct_read() {
    //For correction:
    PathGraphAdj<false> path_graph = PathGraphAdj<false>();
    Path<false> path;
    DnaSequence seq_head, seq_tail;
    KmerInfo<false> first_kmer, last_kmer;
    //Kmer_size
    size_t kmer_size = Parameters::get().kmerSize;
    //Empty solids
    bool no_solids = false;
    if (!_solid.size()){
        no_solids = true;
    }
    else {
        if (_solid.size() == _seq.length()-kmer_size+1)
        {
            return _seq;
        }
        //When first k-mer is not solid
        if (_solid[0].kmer_pos > 0)
        {
            size_t selected_pos = _solid[0].kmer_pos;
            seq_head = _seq.substr(0,selected_pos);
        }
        //When the last k-mer is not solid
        if (_solid[_solid.size()-1].kmer_pos < _seq.length())
        {
            seq_tail = _seq.substr(_solid[_solid.size()-1].kmer_pos+kmer_size,_seq.length());
        }
        for (uint i = 0; i < _solid.size() - 1; ++i) {
            /*Check the head/tail of the read
             * Planned approach: copy head/tail (Naive so far)
             * */
            size_t num_trials = 0;
            for (uint j = i + 1; j < _solid.size() && num_trials < MAX_NUM_TRIALS; ++j) {
                int outcome = check_solids(_solid[i].kmer_pos,
                                           _solid[j].kmer_pos, i, j, num_trials, (j == (i + 1) ? _solid[i].kmer
                                                                                             : _solid[j].kmer));
                if (outcome == 0) {
                    //Extension path
                    char way[MAX_PATH_LEN + 1];
                    size_t ed_score = MAX_PATH_LEN, max_branch = MAX_BRANCH;
                    /*
                     * The path starts at the end of the first k-mer or in the middle, but never at the beginning
                     */
                    size_t start_path = std::min(_solid[i].kmer_pos+kmer_size
                                                 ,_solid[j].kmer_pos);
                    size_t distance = _solid[j].kmer_pos-_solid[i].kmer_pos;
                    size_t len = path.extend(_seq.substr(start_path, distance), _solid[i],
                                             _solid[j], _dbg, &ed_score, way, max_branch);
                    //Add minimum edit path to optimal paths.
                    if (len < MAX_PATH_LEN) {
                        std::string way_string(way);
                        way_string = way_string.substr(0,len-kmer_size+1);
                        /*std::cout << "KmerBegin: "<<_solid[i].kmer.str()<<" "<<_solid[j].kmer.str()<<"\n";
                        std::cout << "Path buscado: "<<_seq.substr(start_path,distance).str()<<"\n";
                        std::cout << "Original seq: "<<_seq.str()<<" "<<start_path<<"\n";
                        std::cout << "Path Found: "<<way<<" "<<way_string<<"\n"<<"SCORE: "<<ed_score<<"\n";*/
                        path_graph.add_edge(_solid[i], _solid[j], ed_score
                                ,(!ed_score)?DnaSequence(_solid[i].kmer.str()+way_string)
                                            :DnaSequence(way_string));
                    } else {
                        //TODO: Crear logs y meter estas salidas ahi
                        /*if (max_branch == 0)
                            std::cout << "Maximum branches reached\n";
                        else
                            std::cout << "No path from source to target " << _solid[i].kmer_pos << "-"
                                      << _solid[j].kmer_pos<< "\n";*/
                        size_t start = _solid[i].kmer_pos+Parameters::get().kmerSize;
                        if (_solid[j].kmer_pos-start <= 0)
                            path_graph.add_edge(_solid[i],_solid[j]
                                    ,0,_seq.substr(_solid[i].kmer_pos,_solid[j].kmer_pos-_solid[i].kmer_pos));
                        else
                            path_graph.add_edge(_solid[i], _solid[j]
                                    ,_solid[j].kmer_pos - start
                                    ,_seq.substr(start,_solid[j].kmer_pos - start));
                    }
                } else {
                    if (outcome == 1) {
                        /*
                         * When there is a jump between 2 nodes at the end of the seq (we need to have at least one copy
                         * of each solid_kmer in the path_graph)
                         */
                        if (j == _solid.size()-1)
                        {
                            //Fix!!!
                            /*std::cout << _solid[i].kmer.str()<<" "<<_solid[i+1].kmer.str()<<" "<<_solid[i+1].kmer_pos << " " <<
                                      _solid[i].kmer.substr(0, _solid[i+1].kmer_pos - _solid[i].kmer_pos).str()<<"\n";*/
                            path_graph.add_edge(_solid[i],_solid[i+1],0
                                    ,_solid[i].kmer.substr(0, _solid[i+1].kmer_pos - _solid[i].kmer_pos));
                            /*
                             * If not it will connect with the last kmer
                             */
                            num_trials = MAX_NUM_TRIALS;
                        }
                        continue;
                    }
                    else if (outcome == 2)
                        break;
                    else {
                        /*std::cout << _solid[i].kmer.str()<<" "<<_solid[j].kmer.str()<<" "<<_solid[j].kmer_pos <<" "<<
                                  _solid[i].kmer.substr(0, _solid[j].kmer_pos - _solid[i].kmer_pos).str()<<"\n";*/
                        path_graph.add_edge(_solid[i], _solid[j], 0,
                                            _solid[i].kmer.substr(0, _solid[j].kmer_pos - _solid[i].kmer_pos));
                        num_trials = MAX_NUM_TRIALS;
                        break;
                    }
                }
            }
            if (num_trials < MAX_NUM_TRIALS)
            {
                size_t fail_len = _solid[_solid.size()-1].kmer_pos-_solid[i].kmer_pos;
                path_graph.add_edge(_solid[i],_solid[_solid.size()-1]
                        ,fail_len,_seq.substr(_solid[i].kmer_pos,_seq.length()));
            }
            if (!path_graph.covered(_solid[i]))
            {
                /*Agregas 5 aristas a los nodos adjacentes*/
                for (uint d = 0; d < (uint)std::min(MAX_NUM_TRIALS,(int)(_solid.size()-i)); ++d) {
                    size_t fail_len = _solid[i + d].kmer_pos - _solid[i].kmer_pos;
                    //TODO: Needs correction
                    path_graph.add_edge(_solid[i], _solid[i + d], 0, _seq.substr(_solid[i].kmer_pos, fail_len));
                }
            }
        }
    }
    /*std::cout << "Vertex of the path graph: "<<path_graph.num_vertex()<<" Edges in the path graph: "
                                                                     << path_graph.num_edges() << "\n";
    //path_graph.show();*/
    if (_solid.size() > 0)
    {
        if (_solid[0].kmer_pos == seq_head.length())
            first_kmer = _solid[0];
        last_kmer = _solid[_solid.size()-1];

    }
    if (!no_solids)
    {
        //std::cout << "Source: "<<first_kmer.kmer.str() << " Target: "<<last_kmer.kmer.str() << "\n";
        DnaSequence path_found = path_graph.shortest_path(first_kmer, last_kmer);
        /*if (_solid[0].kmer_pos != seq_head.length()) {
            std::cout << "Secuencia original: "<<_seq.str() <<"\n";
            std::cout << "Cabeza/Cola: " << seq_head.str() << "-" << seq_tail.str() << "\n";
            std::cout << "Optimal Read: " << path_found.str() << "\n";
        }*/
        DnaSequence full_path(seq_head.str() + path_found.str() + seq_tail.str());
        /*std::cout << "Head: "<<seq_head.str()<<"\nCola: "<<seq_tail.str()<<"\n";
        std::cout << "Secuencia Recuperada: "<<full_path.str()<<" "<<"\n";
        std::cout << "Secuencia Original: "<<_seq.str()<<"\n";*/
        /*
         * PathGraph "Visualization"
         */
        //path_graph.show();
        /*if (_seq.str() != full_path.str())
            std::cout <<"\n"<< _seq.str() << "\n"
                      << full_path.str() << "\n";*/
        return full_path;
    }
    return _seq;
}

//ReadsCorrector
template<>
void ReadCorrector<false>::correct_reads() {
    std::cout << "STAGE: Reads Correction\n";
    #pragma omp parallel
    {
        #pragma omp single
        for (auto &read: _sc.getIndex())
        {
            #pragma omp task shared(read)
            {
                //std::cout << read.first.getId() <<" " << read.second.sequence.length()<<"\n";
                /*if (read.first.getId() != 416)
                    continue;*/
                Progress::update(read.first.getId());
                PathContainer<false> pc(read.first,_dbg,read.second.sequence);
                DnaSequence seq = pc.correct_read();
                _sc.setRead(read.first.getId(),seq);
                //std::cout << "Chain: "<<_sc.getSeq(read.first.getId()).str()<<"\n";
                //std::cout << read.first.getId() << " "<<pc.getSolidLength()<<"\n"
            };
        }
    };
    Progress::update(_sc.getIndex().size());
}

/*
 * Pair_end
 */

template<>
size_t Path<true>::extend_head(const DnaSequence &sub_sequence
        ,KmerInfo<true> t
        ,const DBG<true> &dbg
        ,size_t * score_ed
        ,size_t * max_pos
        ,char * expected_head
        ,size_t &branches
        ,bool behaviour
        ,KmerInfo<true> &which)
{
    size_t best_len = MAX_PATH_LEN;
    return best_len;
}

/*
 * Sub_sequence: part of the read to cover with the DBG
 * h: vector de posibles head-sources
 * t: vector de posibles tails
 * DBG: grafo de bruijn
 * score_ed: mejor score encontrado hasta el momento
 * length_cover: parte de la lectura cubierta
 * */
template<>
size_t Path<true>::extend(const DnaSequence &sub_sequence
        ,KmerInfo<true> h
        ,KmerInfo<true> t
        ,const DBG<true> &dbg
        ,size_t * score_ed
        ,char *expected_path
        ,size_t & branches)
{
    size_t best_len = MAX_PATH_LEN;
    return best_len;
}

//Container path
//TODO: Posiciones salvarlas y empezar la correccion
template<>
size_t PathContainer<true>::check_read()
{
    return 1;
}

/*Mirar para cada Kmer solido todos los kmers siguientes a este susceptibles de trazar
 *un camino entre ellos:
 *  Si son adyacentes -> No se busca -> -1 trial
 *  Si estan solapados -> No se busca
 *  Si estan demasiado lejos no se busca
 * En cualquier otro caso se busca -> -1 trial
 *
 * */
template<>
int PathContainer<true>::check_solids(size_t pos_i, size_t pos_j, size_t i,size_t j
        ,size_t &num_trials, Kmer& kmer)
{
    return 0;
}

//TODO: Correct tail->Optimize all process
template<>
DnaSequence PathContainer<true>::correct_read() {
    return _seq;
}

//ReadsCorrector
template<>
void ReadCorrector<true>::correct_reads() {
    std::cout << "PAIR_END READS CORRECTION!\n";

}