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
 * Sub_sequence: part of the read to cover with the DBG
 * h: vector de posibles head-sources
 * t: vector de posibles tails
 * DBG: grafo de bruijn
 * score_ed: mejor score encontrado hasta el momento
 * length_cover: parte de la lectura cubierta
 * */
size_t Path::extend(const DnaSequence &sub_sequence
        ,std::pair<Kmer, size_t> h
        ,std::pair<Kmer, size_t> t
        ,const DGB &dbg
        ,size_t * score_ed
        ,char *expected_path
        ,size_t & branches)
{
    //Unicamente nos centraremos en el primer Head/Tail -> preguntar manhana
    size_t best_len = MAX_PATH_LEN;
    char path[MAX_PATH_LEN+1];
    std::stack<stack_el*> neighbors;

    std::vector<DnaSequence::NuclType> nts
            = dbg.getNeighbors(h.first);
    std::cout<<"NUM NEIGHBORS: "<<nts.size()<<"\n";
    //Continuacion del kmer actual
    for (uint i = 0; i < nts.size(); ++i)
    {
        Kmer kmer_aux = h.first;
        kmer_aux.appendRight(nts[i]);
        std::cout << "Vecinos "<<nts[i]<<"\n";
        neighbors.push(new stack_el(Kmer(kmer_aux.str()),1));
    }
    //Expected fail length
    size_t fail_len = t.second - h.second;
    while (neighbors.size() > 0 && branches > 0)
    {
        stack_el * el_kmer = neighbors.top();
        neighbors.pop();
        Kmer cur_kmer = el_kmer->kmer;
        size_t pos = el_kmer->pos;

        //Deallocate the stack_el
        //delete el_kmer;

        if (pos >= MAX_PATH_LEN)
            return MAX_PATH_LEN;

        path[pos-1] = cur_kmer.at(0);
        _DP[pos*(MAX_PATH_LEN+1)] = pos;
        size_t min = MAX_PATH_LEN;
        /*
         * Different penalization over indels and subs
         * */
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
        //Comparative between current path and alternative paths already visited
        if (min < (*score_ed)) {
            //Reach the expected node
            if (cur_kmer.str() == t.first.str()) {
                branches--;
                size_t edit_distance = _DP[pos * (MAX_PATH_LEN + 1) + fail_len];
                if (edit_distance < (*score_ed)) {
                    (*score_ed) = edit_distance;
                    std::memcpy(expected_path, path, pos);
                    best_len = pos;
                }
            } else {
                //We havent reached our objective we extend the path again
                nts = dbg.getNeighbors(cur_kmer);
                for (uint i = 0; i < nts.size(); i++){
                    Kmer kmer_aux = Kmer(cur_kmer.str());
                    kmer_aux.appendRight(nts[i]);
                    neighbors.push(new stack_el(Kmer(kmer_aux.str()),pos+1));
                }
                if (!nts.size())
                    branches--;
            }
        }else
            branches--;
    }
    std::cout << "NumOfNeighbors "<<neighbors.size()<< " "<<best_len<< "\n";
    return best_len;
}

//Container path
//TODO: Posiciones salvarlas y empezar la correccion
size_t PathContainer::check_read()
{
    for (auto cur_kmer: IterKmers(_seq)){
        //Optimizar esto
        /*Kmer kmer(cur_kmer.kmer.getSeq().substr(0
                ,cur_kmer.kmer.getSeq().length()));*/
        Kmer kmer = cur_kmer.kmer;
        size_t cur_pos = cur_kmer.kmer_pos;
        bool solid = _dbg.is_solid(kmer);
        if (solid)
            _solid.push_back(std::pair<Kmer, size_t>(kmer, cur_pos));
    }
    /*for (uint i = 0; i < _solid.size(); ++i)
        std::cout << _solid[i].first.str()<<"\n";*/
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
int PathContainer::check_solids(size_t pos_i, size_t pos_j, size_t i,size_t j
        ,size_t &num_trials, Kmer& kmer)
{
    //Kmer are all adjacent.
    if ((pos_j-pos_i) == j-i) {
        num_trials++;
        return 3;
    }
    //Kmers far away from eachother
    if (((pos_j-pos_i)*(1.0+ERROR_RATE)) >= MAX_PATH_LEN) {
        return 2;
    }
    //Both Kmers are overlapped: PREGUNTAR
    if ((pos_j - pos_i) < (1-ERROR_RATE)*Parameters::get().kmerSize) {
        return 1;
    }
    //They are neither too far nor too close -> True
    num_trials++;
    return 0;
}
//TODO: Glue Head/Tails
DnaSequence PathContainer::correct_read() {
    //For correction:
    PathGraphAdj path_graph = PathGraphAdj();
    Path path;
    DnaSequence seq_head, seq_tail;
    Kmer first_kmer, last_kmer;
    //Kmer_size
    size_t kmer_size = Parameters::get().kmerSize;
    std::cout << "Number of solid kmers: "<< _solid.size()<<"\n";
    if (_solid.size() == 0) {
        Kmer source, target;
        bool sb = false;
        for (auto kmer_info: IterKmers(_seq)) {
            if (!sb) {
                source = kmer_info.kmer;
                first_kmer = source;
                sb = true;
                continue;
            }
            target = kmer_info.kmer;
            path_graph.add_edge(source,target,0,source.substr(0,1));
            source = kmer_info.kmer;
        }
        last_kmer = source;
    } else {
        //When first k-mer is not solid
        if (_solid[0].second > 0)
        {
            seq_head = _seq.substr(0,_solid[0].second);
        }
        //When the last k-mer is not solid
        if (_solid[_solid.size()-1].second < _seq.length())
        {
            seq_tail = _seq.substr(_solid[_solid.size()-1].second+kmer_size,_seq.length());
        }
        for (uint i = 0; i < _solid.size() - 1; ++i) {
            /*Check the head/tail of the read
             * Planned approach: copy head/tail (Naive so far)
             * */
            size_t num_trials = 0;
            for (uint j = i + 1; j < _solid.size() && num_trials < MAX_NUM_TRIALS; ++j) {
                int outcome = check_solids(_solid[i].second,
                                           _solid[j].second, i, j, num_trials, (j == (i + 1) ? _solid[i].first
                                                                                             : _solid[j].first));
                if (outcome == 0) {
                    //Extension path
                    char way[MAX_PATH_LEN + 1];
                    size_t ed_score = MAX_PATH_LEN, max_branch = MAX_BRANCH;
                    size_t len = path.extend(_seq.substr(_solid[i].second, _solid[j].second + kmer_size), _solid[i],
                                             _solid[j], _dbg, &ed_score, way, max_branch);
                    //Add minimum edit path to optimal paths.
                    if (len < MAX_PATH_LEN) {
                        std::string way_string(way);
                        //std::cout << "Path Found: "<<way<<" "<<way_string<<"\n";
                        path_graph.add_edge(_solid[i].first, _solid[j].first, len - kmer_size, DnaSequence(way_string));
                    } else {
                        //TODO: Crear logs y meter estas salidas ahi
                        if (max_branch == 0)
                            std::cout << "Maximum branches reached\n";
                        else
                            std::cout << "No path from source to target " << _solid[i].second << "-" << _solid[j].second
                                      << "\n";
                        path_graph.add_edge(_solid[i].first, _solid[j].first,
                                            _solid[j].second - _solid[i].second - Parameters::get().kmerSize,
                                            _seq.substr(_solid[i].second, _solid[j].second - _solid[i].second));
                    }
                } else {
                    //Preguntar por los K-mers que estan muy cercanos
                    if (outcome == 1)
                        continue;
                    else if (outcome == 2)
                        break;
                    else {
                        path_graph.add_edge(_solid[i].first, _solid[j].first, 0,
                                            _solid[i].first.substr(0, _solid[j].second - _solid[i].second));
                    }
                }
            }
            if (num_trials < MAX_NUM_TRIALS)
            {
                size_t fail_len = _solid[_solid.size()-1].second-_solid[i].second;
                path_graph.add_edge(_solid[i].first,_solid[_solid.size()-1].first
                        ,fail_len,_seq.substr(_solid[i].second,_seq.length()));
            }
            if (!path_graph.covered(_solid[i].first))
            {
                /*Esto solo puede ocurrir cuando los kmers estan muy cerca asi que no hace falta if*/
                size_t fail_len = _solid[i+1].second-_solid[i].second;
                path_graph.add_edge(_solid[i].first, _solid[i + 1].first
                        ,0, _seq.substr(_solid[i].second,fail_len));
            }
        }
    }
    std::cout << "Vertex of the path graph: "<<path_graph.num_vertex()<<" Edges in the path graph: "
                                                                     << path_graph.num_edges() << "\n";
    std::cout << "Secuencia original: "<<_seq.str() <<"\n";
    if (_solid.size() > 0)
    {
        first_kmer = _solid[0].first;
        last_kmer = _solid[_solid.size()-1].first;
    }
    std::cout << "Source: "<<first_kmer.str() << " Target: "<<last_kmer.str() << "\n";
    DnaSequence path_found = path_graph.shortest_path(first_kmer,last_kmer);
    std::cout << "Cabeza/Cola: "<<seq_head.str()<<"-"<<seq_tail.str()<<"\n";
    std::cout << "Optimal Read: "<<path_found.str()<<"\n";
    return path_found;
}

//ReadsCorrector
void ReadCorrector::correct_reads() {
    std::cout << "Correcting reads\n";
    for (auto &read: _sc.getIndex())
    {
        if (read.first.getId() > 10)
            continue;
        std::cout << read.first.getId() <<" " << read.second.sequence.length()<<"\n";
        Progress::update(read.first.getId());
        PathContainer pc(read.first,_dbg,read.second.sequence);
        DnaSequence seq = pc.correct_read();
        _sc.setRead(read.first.getId(),seq);
        std::cout << "Chain: "<<_sc.getSeq(read.first.getId()).str()<<"\n";
        //std::cout << read.first.getId() << " "<<pc.getSolidLength()<<"\n";
    }
    Progress::update(_sc.getIndex().size());
}