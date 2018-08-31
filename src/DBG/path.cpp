#include "path.h"
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stack>

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
        ,size_t & score_ed
        ,size_t & length_cover
        ,char *expected_path
        ,size_t & branches)
{
    //Unicamente nos centraremos en el primer Head/Tail -> preguntar manhana
    size_t best_len = -1;
    char path[MAX_PATH_LEN+1];
    std::stack<stack_el*> neighbors;

    std::vector<std::pair<Kmer,DnaSequence::NuclType>> nts
            = dbg.getNeighbors(h.first);

    //Continuacion del kmer actual
    for (int i = 0; i < (int) nts.size(); ++i)
        neighbors.push(new stack_el(nts[i].first,1));

    std::cout << sub_sequence.str() << "\n";
    std::cout << t.first.str() << " " <<" "<<h.first.str() <<"\n";
    //Expected fail length
    size_t fail_len = t.second - h.second;
    while (neighbors.size() > 0 && branches > 0)
    {
        stack_el * el_kmer = neighbors.top();
        neighbors.pop();
        Kmer cur_kmer = el_kmer->kmer;
        size_t pos = el_kmer->pos;

        //Deallocate the stack_el
        delete el_kmer;

        if (pos >= MAX_PATH_LEN){
            std::cout << "Path extremely large\n";
            return 0;
        }

        std::cout << cur_kmer.str()<< " "<<pos << "\n";
        path[pos-1] = cur_kmer.at(0);
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
        //Comparative between current path and alternative paths already visited
        if (min < score_ed) {
            //Reach the expected node
            if (cur_kmer.str() == t.first.str()) {
                branches--;
                size_t edit_distance = _DP[pos * (MAX_PATH_LEN + 1) + fail_len];
                if (edit_distance < score_ed) {
                    score_ed = edit_distance;
                    std::cout << "TE ENCONTREEEE" << score_ed << "\n";
                    std::memcpy(expected_path, path, pos);
                    best_len = pos;
                }
            } else {
                //We havent reached our objective we extend the path again
                nts = dbg.getNeighbors(cur_kmer);
                for (uint i = 0; i < nts.size(); i++)
                    neighbors.push(new stack_el(nts[i].first, pos + 1));
                std::cout << "Num neighbors " << nts.size() << "\n";
                if (!nts.size())
                    branches--;
                std::cout << "Neighbors size: " << neighbors.size() << "\n";
            }
        }else
            branches--;
    }
    return best_len;
}

//Container path
//TODO: Posiciones salvarlas y empezar la correccion
size_t PathContainer::check_read(){
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
 * Mirar lo de como copiar Kmers eso tengo que resolverlo cuanto antes...
 * */
int PathContainer::check_solids(size_t pos_i, size_t pos_j, size_t i,size_t j
        ,size_t &num_trials, std::vector<Kmer> &checked_solid, Kmer& kmer)
{
    //Kmer are all adjacent.
    if ((pos_j-pos_i) == j-i) {
        checked_solid.push_back(kmer);
        num_trials++;
        return 1;
    }
    //Kmers far away from eachother
    if ((pos_j-pos_i)*(1.0+ERROR_RATE) >= MAX_DISTANCE_READ) {
        std::cout << pos_j<< " " << pos_i<<" Demasiado lejos\n";
        return 2;
    }
    //Both Kmers are overlapped: PREGUNTAR
    if ((pos_j - pos_i) < (1-ERROR_RATE)*Parameters::get().kmerSize) {
        std::cout << "Demasiado cerca\n";
        return 1;
    }
    //They are neither too far nor too close -> True
    num_trials++;
    return 0;
}

void PathContainer::correct_read() {
    PathGraphAdj pathGraphAdj();
    std::vector<Kmer> checked_solids;
    std::vector<DnaSequence> corrected_paths;
    Path path(_seq);
    if (_solid.size() == 0)
        for (auto kmer_info: IterKmers(_seq))
            checked_solids.push_back(kmer_info.kmer);
    else
        for (int i = 0; i < (int) _solid.size()-1; ++i) {
            //Check the head of the read
            size_t num_trials = 0;
            for (int j = i + 1; j < (int)_solid.size() && num_trials < MAX_NUM_TRIALS; ++j){
                int outcome =check_solids(_solid[i].second,
                                          _solid[j].second,(size_t)i,(size_t)j
                        ,num_trials,checked_solids, (j == (i+1)?_solid[i].first
                                                               :_solid[j].first));
                if (!outcome)
                {
                    //Extension path
                    char way[MAX_PATH_LEN+1] = "None\0";
                    size_t ed_score = MAX_PATH_LEN, length_cover, max_branch = MAX_BRANCH;
                    size_t kmer_size = Parameters::get().kmerSize;
                    size_t len = path.extend(_seq.substr(_solid[i].second
                            ,_solid[j].second+kmer_size)
                            ,_solid[i], _solid[j],_dbg
                            ,ed_score, length_cover
                            ,way, max_branch);
                    if (len != 0){
                        std::cout << "Path Found: "<<way<<"\n";
                        corrected_paths.push_back(DnaSequence(std::string(way)));
                    }
                    else {
                        if (max_branch == 0)
                            std::cout << "Maximum branches reached\n";
                        else
                            std::cout << "No path from source to target\n";
                    }
                }else{
                    if (outcome == 1)
                        continue;
                    else
                        break;
                }
            }
        }
    sleep(10000);
    std::cout << "Size of checked solids "<<checked_solids.size() << " "
              << " Size of kmer_solids " << _solid.size()<< "\n";
}

//ReadsCorrector
void ReadCorrector::correct_reads() {
    std::cout << "Correcting reads\n";
    for (auto &read: _sc.getIndex()){
        Progress::update(read.first.getId());
        PathContainer pc(read.first,_dbg,read.second.sequence);
        pc.correct_read();
        //std::cout << read.first.getId() << " "<<pc.getSolidLength()<<"\n";
    }
    Progress::update(_sc.getIndex().size());
}