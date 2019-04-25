#include "path.h"
#include <math.h>

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
 */
template<>
size_t Path<false>::extend_head(const DnaSequence &sub_sequence
        ,KmerInfo<false> t
        ,const DBG<false> * dbg
        ,size_t * score_ed
        ,char * expected_head
        ,size_t &branches
        ,size_t * distance
        ,bool behave)
{
    size_t best_len = MAX_PATH_LEN,
            kmerSize = Parameters::get().kmerSize,
            startPos = t.kmer_pos+kmerSize-2;
    char * path = (char*) malloc ( sizeof(char) * (MAX_PATH_LEN+1));
    stack<stack_el*> neighbors;
    vector<Kmer> nts = (behave)?dbg->getInKmerNeighbors(t.kmer)
                               :dbg->getKmerNeighbors(t.kmer);
    //Continuacion del kmer actual
    for (auto k: nts)
    {
        neighbors.push(new stack_el(k,1,k.at((behave)?0:kmerSize-2)));
    }
    //Expected fail length (maximum fail length)
    size_t fail_len = (behave)?t.kmer_pos:(sub_sequence.length()-startPos-2);
    size_t offset = ceil(DELTA_EXTENSION*fail_len);
    size_t maximal_allowed_extension = fail_len+offset;
    if (maximal_allowed_extension >= MAX_PATH_LEN)
    {
        free(path);
        return MAX_PATH_LEN;
    }
    //For each possible head we are going to look for a path
    while (neighbors.size() > 0 && branches > 0)
    {
        stack_el * el_kmer = neighbors.top();
        neighbors.pop();
        Kmer cur_kmer = el_kmer->kmer;
        size_t pos = el_kmer->pos;
        DnaSequence::NuclType nuc = DnaSequence::nfi(el_kmer->nuc);
        //Deallocate the stack_el
        delete el_kmer;
        if (pos >= maximal_allowed_extension)
        {
            branches--;
            continue;
        }

        size_t pointInPath = (behave)?(maximal_allowed_extension-pos):(pos-1), cover = 0;
        path[pointInPath] = nuc;
        _DP[pos*(MAX_PATH_LEN+1)] = pos;
        size_t min = MAX_PATH_LEN+1;
        int64_t pivote = fail_len;
        if (behave) {
            /*cout << "Curr: "<<cur_kmer.str()<<" Diagonal: "<<_DP[pos * (MAX_PATH_LEN + 1) + pivote - (MAX_PATH_LEN + 1) + 1]
                 <<" Right: "<<_DP[pos * (MAX_PATH_LEN + 1) + pivote + 1]<<" Izquierda: "<<_DP[pos * (MAX_PATH_LEN + 1) + pivote - (MAX_PATH_LEN + 1)]<<endl;*/
            size_t extra_i = 0, extra_d = 0, extra_r = 0;
            for (int k = pivote; k > 0; --k) {
                size_t ev_pos = pos * (MAX_PATH_LEN + 1) + k;
                if ((pos - 1) == 0) {
                    extra_i = (pivote+1) - _DP[ev_pos - (MAX_PATH_LEN + 1)];
                    extra_d = (pivote+1) - _DP[ev_pos - (MAX_PATH_LEN + 1) + 1];
                }else{
                    extra_i = _DP[ev_pos - (MAX_PATH_LEN + 1)];
                    extra_d = _DP[ev_pos - (MAX_PATH_LEN + 1) + 1];
                }
                if ((pivote == k)) {
                    extra_r = pos;
                    extra_d += (pos-1);
                }else
                    extra_r = _DP[ev_pos + 1];
                (path[pointInPath] == sub_sequence.at(k-1))?
                        _DP[ev_pos] = std::min(extra_d,std::min(extra_i + 1,extra_r + 1)):
                        _DP[ev_pos] = std::min(extra_d + 1,std::min(extra_i + 1,extra_r + 1));
                if (_DP[ev_pos] <= min) {
                    min = _DP[ev_pos];
                    cover = pivote - k + 1;
                }
            }
        }else {
            for (int64_t k = 1; k <= (int) pivote; k++) {
                size_t ev_pos = pos * (MAX_PATH_LEN + 1) + k;
                (path[pointInPath] == sub_sequence.at(startPos + k + 1)) ?
                        _DP[ev_pos] = std::min(
                                _DP[ev_pos - (MAX_PATH_LEN + 1) - 1],
                                std::min(_DP[ev_pos - (MAX_PATH_LEN + 1)], _DP[ev_pos - 1]) + 1) :
                        _DP[ev_pos] = std::min(
                                _DP[ev_pos - (MAX_PATH_LEN + 1) - 1] + 1,
                                std::min(_DP[ev_pos - (MAX_PATH_LEN + 1)], _DP[ev_pos - 1]) + 1);
                if (_DP[ev_pos] <= min)
                {
                    min = _DP[ev_pos];
                    cover = k;
                }
            }
        }
        /*
        * If edit_distance es min y además es la pos más alta aceptado
        */
        if (fail_len - cover + min <= (*score_ed))
        {
            (*score_ed) = fail_len - cover + min;
            path[(behave) ? maximal_allowed_extension : pos] = '\0';
            (behave) ? memcpy(expected_head, &path[maximal_allowed_extension - pos], pos + 1)
                     : memcpy(expected_head, path, pos + 1);
            best_len = 0;
            (*distance) = cover;
        }
        /*
         * Lets continue extension:
         */
        nts = (behave) ? dbg->getInKmerNeighbors(cur_kmer)
                       : dbg->getKmerNeighbors(cur_kmer);
        for (auto k:nts) {
            neighbors.push(new stack_el(k, pos + 1, k.at((behave) ? 0 : kmerSize - 2)));
        }
        if (!nts.size()) {
            branches--;
        }
            /*if (fail_len <= pos)
            {
                branches--;
                size_t edit_distance = (behave)?_DP[pos * (MAX_PATH_LEN + 1) + 1]
                                               :_DP[pos*(MAX_PATH_LEN+1)+pivote];
                if (edit_distance < *(score_ed))
                {
                    (*score_ed) = edit_distance;
                    path[(behave)?maximal_allowed_extension:pos] = '\0';
                    (behave)?memcpy(expected_head, &path[maximal_allowed_extension-pos], pos+1)
                            :memcpy(expected_head,path,pos+1);
                    best_len = 0;
                    (*distance) = pos;
                }
            }else{
                nts = (behave)?dbg->getInKmerNeighbors(cur_kmer)
                              :dbg->getKmerNeighbors(cur_kmer);
                for (auto k:nts){
                    if (processed.find(k) == processed.end())
                    {
                        processed.emplace(k);
                        neighbors.push(new stack_el(k,pos+1,k.at((behave)?0:kmerSize-2)));
                    }
                }
                if (!nts.size())
                {
                    branches--;
                    if (pos < fail_len)
                    {
                        size_t edit_distance = pivote - (behave)?(_DP[pos * (MAX_PATH_LEN + 1) + pivote - pos]-1)
                                                       :_DP[pos*(MAX_PATH_LEN + 1)+pos];
                        if (edit_distance < (*score_ed) && pos > (*distance))
                        {
                            (*score_ed) = edit_distance;
                            path[(behave)?maximal_allowed_extension:pos] = '\0';
                            (behave)?memcpy(expected_head, &path[maximal_allowed_extension - pos], pos +1)
                                    :memcpy(expected_head, path, pos+1);
                            best_len = 0;
                            (*distance) = pos;
                        }
                    }
                }
            }*/
    }
    for (uint i = 0; i > neighbors.size(); ++i)
    {
        stack_el * top_el = neighbors.top();
        neighbors.pop();
        delete top_el;
    }
    free(path);
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
        ,const DBG<false> * dbg
        ,size_t * score_ed
        ,char *expected_path
        ,size_t & branches
        ,char * alternative_path
        ,size_t * distance
        ,size_t * alternative_ed)
{
    size_t best_len = MAX_PATH_LEN, kmerSize = Parameters::get().kmerSize;
    size_t startPos = h.kmer_pos + kmerSize;
    Kmer kmer_objective = Kmer(t.kmer.substr(1,kmerSize-1));
    char * path = (char*) malloc ( sizeof(char) * (MAX_PATH_LEN+1));
    stack<stack_el*> neighbors;
    /*cout << "Kmer partida: "<<h.kmer.str()<<endl;
    cout << "Kmer final: "<<kmer_objective.str()<<endl;*/
    /*
     * Extract (k-1)mer from header to start the extension
     */
    vector<Kmer> nts;
    //Continuacion del kmer actual
    nts = dbg->getKmerNeighbors(h.kmer);
    for (auto k: nts) {
        neighbors.push(new stack_el(k, 1, k.at(kmerSize - 2)));
    }

    //Expected fail length (maximum fail length)
    size_t fail_len = t.kmer_pos+kmerSize - startPos;
    //size_t maximal_ed = ceil(fail_len*DELTA_EXTENSION)+1;
    /*cout << "SOURCE K-mer: "<<h.kmer.str()<<" TARGET: "<<kmer_objective.str()<<endl;
    cout << "Pos1: "<<h.kmer_pos<<" Pos2: "<<t.kmer_pos<<endl;
    cout << "OVERLAP: "<<overlaped<<endl;
    cout << "FAIL: "<<fail_len<< " Maximal_allowed: "<<maximal_ed<<endl;*/

    while ((neighbors.size() > 0) && (branches > 0))
    {
        stack_el * el_kmer = neighbors.top();
        neighbors.pop();
        Kmer cur_kmer = el_kmer->kmer;
        size_t pos = el_kmer->pos;
        DnaSequence::NuclType nuc = DnaSequence::nfi(el_kmer->nuc);
        //Deallocate the stack_el
        delete el_kmer;

        if (pos >= MAX_PATH_LEN)
        {
            best_len = MAX_PATH_LEN;
            branches--;
            continue;
        }

        size_t pointInPath = pos-1;
        path[pointInPath] = nuc;
        _DP[pos*(MAX_PATH_LEN+1)] = pos;
        size_t min = MAX_PATH_LEN, cover = 0;
        for (uint k = 1; k <= (fail_len); ++k)
        {
            size_t ev_pos = pos*(MAX_PATH_LEN+1)+k;
            (path[pos-1] == sub_sequence.at(startPos+k-1))?
                    _DP[ev_pos] = std::min(
                            _DP[ev_pos-(MAX_PATH_LEN+1)-1]
                            ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                                    ,_DP[ev_pos-1])+1):
                    _DP[ev_pos] = std::min(
                            _DP[ev_pos-(MAX_PATH_LEN+1)-1]+1
                            ,std::min(_DP[ev_pos-(MAX_PATH_LEN+1)]
                                    ,_DP[ev_pos-1])+1);
            if (_DP[ev_pos] <= min)
            {
                min = _DP[ev_pos];
                cover = k;
            }
            /*if (overlaped)
                cout << sub_sequence.at(startPos+k-1)<<" "<<path[pos-1]<<endl<<"TARGET: "<<kmer_objective.str()<< " CUR: "<<cur_kmer.str()<<endl;*/
        }
        if (min > (*score_ed) || min > (*alternative_ed)) {
            branches--;
            continue;
        }
        if (fail_len - cover + min <= (*alternative_ed))
        {
            (*alternative_ed) = fail_len - cover + min;
            path[pos] = '\0';
            /*cout << "Path: "<<path<<" ED: "<<(*alternative_ed)<<" Cover: "<<cover<<endl;
            cout << " "<<sub_sequence.at(startPos+cover-1)<<endl;*/
            memcpy(alternative_path, path, pos+1);
            (*distance) = cover;
        }
        //Comparative between current path and alternative paths already visited
        if (min <= (*score_ed)) {
            //Reach the expected node
            if (cur_kmer == kmer_objective) {
                size_t edit_distance = _DP[pos * (MAX_PATH_LEN + 1) + fail_len];
                if (edit_distance < (*score_ed)) {
                    (*score_ed) = edit_distance;
                    path[pos] = '\0';
                    memcpy(expected_path,path,pos+1);
                    best_len = pos;
                }
            } else {
                //We have not reached our objective we extend the path again
                nts = dbg->getKmerNeighbors(cur_kmer);
                for (auto k: nts)
                {
                    neighbors.push(new stack_el(k,pos+1,k.at(kmerSize-2)));
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
    free(path);
    //(*alternative_ed) = (*alternative_ed) - kmerSize + (fail_len-kmerSize-((*distance)-1));
    return best_len;
}

//Container path
template<>
size_t PathContainer<false>::check_read()
{
    bool occupied = false, origin = false;
    size_t remaining_solids = 0, seq_length = _seq.length(), followed_solid = 0;
    KmerInfo<false> lka, pks;
    for (auto cur_kmer: IterKmers<false>(_seq))
    {
        if (_dbg->is_solid(cur_kmer.kmer))
        {
            if (cur_kmer.kmer_pos == 0)
                _head = false;
            if (cur_kmer.kmer_pos == (seq_length-Parameters::get().kmerSize))
                _tail = false;
            if (remaining_solids)
            {
                lka = cur_kmer;
                --remaining_solids;
                _solid.push_back(cur_kmer);
            }else if (!origin) {
                origin = true;
            }
            pks = cur_kmer;
            occupied = true;
            followed_solid++;
        }else{
            if (occupied)
            {
                if (_solid.size() == 0)
                {
                    if (followed_solid > 1)
                        _solid.push_back(pks);
                }
                else if (lka != pks)
                {
                    _solid.push_back(pks);
                    if (origin){
                        origin = false;
                        _origin.emplace(pks);
                    }
                }else if (followed_solid == 1)
                {
                    _solid.pop_back();
                }
                occupied = false;
            }
            followed_solid = 0;
            remaining_solids = MAX_NUM_TRIALS;
        }
    }
    /*if (_readId == 27472){
    for (uint i = 0; i < _solid.size(); ++i)
        std::cout << _solid[i].kmer.str()<<" "<<_solid[i].kmer_pos<<"\n";
    for (auto k:_origin)
        cout << k.kmer.str()<<" "<<k.kmer_pos<<endl;
    std::cout << "NumSolids: "<<_solid.size()<<"\n";
    std::cout << "HEAD: "<< _head<<endl;
    std::cout << "TAIL: "<<_tail<<endl;}*/
    return 1;
}

/*Mirar para cada Kmer solido todos los kmers siguientes a este susceptibles de trazar
 *un camino entre ellos:
 *  Si son adyacentes -> No se busca -> -1 trial
 *  Si estan solapados -> No se busca
 *  Si estan demasiado lejos no se busca
 * En cualquier otro caso se busca -> -1 trial
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
    //They are neither too far nor too close -> True
    if ((pos_j-pos_i) < (Parameters::get().kmerSize*(1-MAX_ERROR_RATE))){
        return 1;
    }
    num_trials++;
    return 0;
}

template<>
DnaSequence PathContainer<false>::correct_read() {
    //For correction:
    PathGraphAdj<false> path_graph = PathGraphAdj<false>();
    Path<false> path;
    DnaSequence seq_head, seq_tail, full_path;
    KmerInfo<false> first_kmer, last_kmer;

    size_t kmer_size = Parameters::get().kmerSize;
    bool no_solids = false;
    if (!_solid.size()){
        no_solids = true;
    }
    else {
        if (_solid.size() == _seq.length()-kmer_size+1)
        {
            return _seq;
        }
        //First k-mer - not solid
        if (_head && kmer_size >= 15)
        {
            char * way = (char*) malloc (sizeof(char) * (MAX_PATH_LEN+1));
            size_t score_ed = MAX_PATH_LEN,selected = 0, distance = 0;
            for (size_t i = 0; i < (size_t)std::min(MAX_NUM_TRIALS, (int)_solid.size()); ++i)
            {
                size_t branches = MAX_BRANCH,
                        new_len = path.extend_head(_seq,_solid[i],_dbg,&score_ed,way,branches,&distance,true);
                if (new_len == MAX_PATH_LEN)
                    continue;
                if (new_len == 0)
                {
                    selected = i;
                }
            }
            string way_string(way);
            if (distance < _solid[selected].kmer_pos)
                seq_head = DnaSequence(_seq.substr(0, _solid[selected].kmer_pos - (distance)).str()+way_string);
            else
                seq_head = DnaSequence(way_string);
            /*size_t max_allowed_ed = ceil(_solid[selected].kmer_pos*ERROR_HEADTAIL_RATE)+1;
            string way_string(way);
            seq_head = ((score_ed <= max_allowed_ed))?
                       ((distance < _solid[selected].kmer_pos)?
                       DnaSequence(_seq.substr(0,_solid[selected].kmer_pos-distance).str()+way_string):DnaSequence(way_string))
                                                     :_seq.substr(0, _solid[0].kmer_pos);*/
            /*cout <<"BASIC HEAD: "<<_seq.substr(0,_solid[selected].kmer_pos).str()<< "\n HEAD: "
                 <<seq_head.str()<<" "<<_solid[selected].kmer_pos<<"\n"<<way<<" "<<score_ed<<endl
                 <<"Distance: "<<distance<<" Origin: "<<_solid[selected].kmer_pos<<endl;*/
            _solid.erase(_solid.begin(), _solid.begin()+selected);
            free(way);
        }else
        {
            if (_solid[0].kmer_pos > 0)
                seq_head = _seq.substr(0,_solid[0].kmer_pos);
        }
        //Last k-mer - not solid
        if (_tail && kmer_size >= 15)
        {
            char * way = (char*) malloc (sizeof(char) * (MAX_PATH_LEN+1));
            size_t score_ed = MAX_PATH_LEN,selected = 0, distance = 0;
            for (size_t i = 0; i < (size_t)std::min(MAX_NUM_TRIALS,(int)_solid.size()); ++i)
            {
                size_t branches = MAX_BRANCH,
                        new_len = path.extend_head(_seq,_solid[_solid.size()-1-i],_dbg,&score_ed,way,branches, &distance,false);
                if (new_len == MAX_PATH_LEN)
                    continue;
                if (new_len == 0)
                {
                    selected = i;
                }
            }
            string way_string(way);
            size_t index = _solid.size()-1-selected;
            if (distance < (_seq.length() - _solid[selected].kmer_pos))
                seq_tail = DnaSequence(way_string+_seq.substr(_solid[index].kmer_pos+kmer_size+distance, _seq.length()).str());
            else
                seq_tail = DnaSequence(way_string);
            /*size_t max_allowed_ed = ceil(_solid[index].kmer_pos*ERROR_HEADTAIL_RATE)+1;
            seq_tail = ((score_ed <= max_allowed_ed))?
                       ((distance < (_seq.length()-_solid[index].kmer_pos-kmer_size))?
                       DnaSequence(way_string+_seq.substr(_solid[index].kmer_pos+kmer_size+distance,_seq.length()).str()):DnaSequence(way_string))
                                                     :_seq.substr(_solid[_solid.size()-1].kmer_pos+kmer_size,_seq.length());*/
            /*cout <<"BASIC TAIL: "<< _seq.substr(_solid[index].kmer_pos+kmer_size,_seq.length()).str()
                 <<"\n TAIL: "<<seq_tail.str()<<" "<<_solid[index].kmer_pos<<"\nWAY: "<<way<<" SE:"<<score_ed<<endl;*/
            _solid.resize(_solid.size()-selected);
            free(way);
        }else
        {
            if (_solid[_solid.size()-1].kmer_pos+kmer_size < _seq.length())
                seq_tail = _seq.substr(_solid[_solid.size()-1].kmer_pos+kmer_size,_seq.length());
        }
        for (uint i = 0; i < _solid.size() - 1; ++i) {
            /*Check the head/tail of the read
             * Planned approach: copy head/tail (Naive so far)
             * */
            size_t num_trials = 0;
            for (uint j = i + 1; (j < _solid.size()) && (num_trials < MAX_NUM_TRIALS); ++j) {
                //cout << "Solid: "<<_solid[i].kmer.str()<<":"<<_solid[i].kmer_pos<<" "<<_solid[j].kmer.str()<<":"<<_solid[j].kmer_pos<<endl;
                int outcome = check_solids(_solid[i].kmer_pos,
                                           _solid[j].kmer_pos, i, j, num_trials, (j == (i + 1) ? _solid[i].kmer
                                                                                             : _solid[j].kmer));
                if (_origin.find(_solid[j]) != _origin.end()){
                    num_trials = MAX_NUM_TRIALS;
                    size_t separation = _solid[j].kmer_pos-_solid[i].kmer_pos;
                    path_graph.add_edge(_solid[i], _solid[j], separation,
                                        _seq.substr(_solid[i].kmer_pos+kmer_size, separation));
                    continue;
                }
                if (outcome == 0) {
                    //Extension path
                    char * way = (char*) malloc (sizeof(char) * (MAX_PATH_LEN+1)),
                            * alternative_way = (char*) malloc (sizeof(char) * (MAX_PATH_LEN+1));
                    size_t ed_score = MAX_PATH_LEN, max_branch = MAX_BRANCH, distance = 0, alternative_ed = MAX_PATH_LEN;
                    /*
                     * The path starts at the end of the first k-mer or in the middle, but never at the beginning
                     */
                    size_t len = MAX_PATH_LEN;
                    len = path.extend(_seq, _solid[i],
                                  _solid[j], _dbg, &ed_score, way,
                                  max_branch, alternative_way, &distance, &alternative_ed);
                    //Add minimum edit path to optimal paths.
                    if ((len < MAX_PATH_LEN)) {
                        string way_string(way);
                        if ((_solid[i].kmer_pos + kmer_size) <= _solid[j].kmer_pos)
                        {
                            path_graph.add_edge(_solid[i], _solid[j], ed_score
                                    ,DnaSequence(way_string));
                        }else{
                            path_graph.add_edge(_solid[i], _solid[j], ed_score
                                    ,DnaSequence(way_string));
                        }
                        /*size_t startPath =_solid[i].kmer_pos + kmer_size,
                                distance = _solid[j].kmer_pos - _solid[i].kmer_pos;
                        cout << "ReadId: "<<_readId<<endl;
                        std::cout << "Original seq: " << _seq.str() << " " << startPath << " " << " "
                                  << _solid[j].kmer_pos << " " << distance << "\n";
                        std::cout << "Path buscado: " << _seq.substr(startPath, distance).str() << "\n";
                        std::cout << "Original seq: " << _seq.str() << " " << startPath << "\n";
                        std::cout << "Path Found: " << way << " " << way_string << "\n" << "SCORE: " << ed_score<< "\n";*/
                    } else {
                        string way_string(alternative_way);
                        /*if (max_branch == 0)
                            std::cout << "Maximum branches reached"<<endl<<" Num_Read: "<<_readId
                                      <<endl<<"Path: "<<way_string<<endl<<"Distance: "<<distance<<endl<<"EditScore: "<<alternative_ed<<endl;
                        else if (distance != 0)
                            std::cout << "No path from source to target " << _solid[i].kmer_pos << "-"
                                      << _solid[j].kmer_pos<< "\n"<<_readId<<endl<<" Num_Read: "<<_readId
                                      <<endl<<"Path: "<<way_string<<endl<<"Distance: "<<distance<<endl<<"EditScore: "<<alternative_ed<<endl;*/
                        /*else
                            cout << "Case weird: "<<max_branch<<" Distance: "<<distance<<_readId
                                 <<" Places: "<<_solid[i].kmer_pos<<"-"<<_solid[j].kmer_pos<<endl<<"ORIGIN: "<<(_origin.find(_solid[j]) != _origin.end())<<endl;*/
                        size_t start = _solid[i].kmer_pos+kmer_size;
                        if (_solid[j].kmer_pos < start)
                        {
                            size_t separation = _solid[j].kmer_pos -_solid[i].kmer_pos;
                            if (distance != 0) {
                                string final_chain = way_string +
                                                     ((distance >= separation) ? "" :
                                                      _seq.substr(_solid[i].kmer_pos+kmer_size + distance,separation-distance).str());
                                path_graph.add_edge(_solid[i], _solid[j], alternative_ed, DnaSequence(final_chain));
                                /*cout << "Final(1): "<<final_chain<<endl<<alternative_ed<<endl;
                                cout << "Final(2): "<<_seq.substr(_solid[i].kmer_pos+kmer_size, _solid[j].kmer_pos-_solid[i].kmer_pos).str()<<endl;*/
                            } else {
                                path_graph.add_edge(_solid[i], _solid[j],0,
                                                    _seq.substr(_solid[i].kmer_pos+kmer_size,(_solid[j].kmer_pos - _solid[i].kmer_pos)));
                            }
                        }
                        else
                        {
                            if (distance != 0)
                            {
                                size_t cover = _solid[j].kmer_pos-_solid[i].kmer_pos;
                                string final_chain = way_string+((distance >= cover)?"":
                                                     _seq.substr(_solid[i].kmer_pos+kmer_size+distance,cover-distance).str());
                                /*string final_chain = _solid[i].kmer.str()+way_string.substr(0, (distance >= (_solid[j].kmer_pos - start))?_solid[j].kmer_pos - start:way_string.length())
                                                    + ((distance >= (_solid[j].kmer_pos-start))?"":_seq.substr(start+distance, _solid[j].kmer_pos - start - distance).str());*/
                                path_graph.add_edge(_solid[i], _solid[j], alternative_ed, DnaSequence(final_chain));
                                /*cout << "Final-(1): "<<DnaSequence(final_chain).str()<<endl<<alternative_ed<<endl<<"Way: "<<way_string<<endl;
                                cout << "Distance: "<<distance<<" Max: "<<(_solid[j].kmer_pos-_solid[i].kmer_pos)<<endl;
                                cout << "Final-(2): "<<_seq.substr(_solid[i].kmer_pos+kmer_size, _solid[j].kmer_pos-_solid[i].kmer_pos).str()<<endl;*/
                            }else{
                                path_graph.add_edge(_solid[i], _solid[j]
                                        ,0,_seq.substr(_solid[i].kmer_pos, _solid[j].kmer_pos-_solid[i].kmer_pos));
                            }
                        }
                    }
                    free(alternative_way);
                    free(way);
                } else {
                    if (outcome == 1) {
                        /*
                         * When there is a jump between 2 nodes at the end of the seq (we need to have at least one copy
                         * of each solid_kmer in the path_graph)
                         */
                        if (j == _solid.size() - 1) {
                            //Fix!!!
                            /*std::cout << _solid[i].kmer.str()<<" "<<_solid[i+1].kmer.str()<<" "<<_solid[i+1].kmer_pos << " " <<
                                      _solid[i].kmer.substr(0, _solid[i+1].kmer_pos - _solid[i].kmer_pos).str()<<"\n";*/
                            path_graph.add_edge(_solid[i], _solid[i + 1], _solid[i + 1].kmer_pos - _solid[i].kmer_pos,
                                                _solid[i].kmer.substr(0, _solid[i + 1].kmer_pos - _solid[i].kmer_pos));
                            num_trials = MAX_NUM_TRIALS;
                        }
                        continue;
                    } else if (outcome == 2){
                        path_graph.add_edge(_solid[i], _solid[j], _solid[j].kmer_pos-_solid[i].kmer_pos
                                ,_seq.substr(_solid[i].kmer_pos+kmer_size, (_solid[j].kmer_pos-_solid[i].kmer_pos)));
                        break;
                    }else {
                        /*std::cout << _solid[i].kmer.str()<<" "<<_solid[j].kmer.str()<<" "<<_solid[j].kmer_pos <<" "<<
                                  _solid[i].kmer.substr(0, _solid[j].kmer_pos - _solid[i].kmer_pos).str()<<"\n";*/
                        path_graph.add_edge(_solid[i], _solid[j], 0,
                                            _solid[i].kmer.substr(0, _solid[j].kmer_pos - _solid[i].kmer_pos));
                    }
                }
            }
            if (num_trials < MAX_NUM_TRIALS)
            {
                size_t fail_len = _solid[_solid.size()-1].kmer_pos-_solid[i].kmer_pos;
                path_graph.add_edge(_solid[i],_solid[_solid.size()-1]
                        ,fail_len,_seq.substr(_solid[i].kmer_pos+kmer_size,fail_len));
            }
            if (!path_graph.covered(_solid[i]))
            {
                size_t start = _solid[i].kmer_pos+Parameters::get().kmerSize;
                for (uint d = 1; d < (uint)std::min(MAX_NUM_TRIALS,(int)(_solid.size()-i)); ++d) {
                    if (_solid[i+d].kmer_pos-start <= 0)
                    {
                        path_graph.add_edge(_solid[i], _solid[i + d], _solid[i + d].kmer_pos - _solid[i].kmer_pos,
                                            _seq.substr(_solid[i].kmer_pos+kmer_size, (_solid[i+d].kmer_pos - _solid[i].kmer_pos)));
                    }else
                    {
                        path_graph.add_edge(_solid[i], _solid[i+d]
                                ,_solid[i+d].kmer_pos - start
                                ,_seq.substr(_solid[i].kmer_pos+kmer_size,(_solid[i+d].kmer_pos-_solid[i].kmer_pos)));
                    }
                }
            }
        }
    }
    /*std::cout << "Vertex of the path graph: "<<path_graph.num_vertex()<<" Edges in the path graph: "
                                                                     << path_graph.num_edges() << "\n";
    path_graph.show();*/
    if (_solid.size() > 0)
    {
        first_kmer = _solid[0];
        last_kmer = _solid[_solid.size()-1];
    }
    if (!no_solids)
    {
        //std::cout << "Source: "<<first_kmer.kmer.str() << " Target: "<<last_kmer.kmer.str() << "\n";
        //DnaSequence path_found = path_graph.shortest_path(first_kmer, last_kmer);
        /*if (_solid[0].kmer_pos != seq_head.length()) {
            std::cout << "Secuencia original: "<<_seq.str() <<"\n";
            std::cout << "Cabeza/Cola: " << seq_head.str() << "-" << seq_tail.str() << "\n";
            std::cout << "Optimal Read: " << path_found.str() << "\n";
        }*/
        full_path = DnaSequence(seq_head.str() + path_graph.shortest_path(first_kmer, last_kmer).str() + seq_tail.str());
        /*if (full_path==DnaSequence("AAAAGTGAATCAGAGTTAGTCAGTCAAATAATAAAGCAGTTAATAAAAAAGGAAAAGGTCTACCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGAAGAACATGAGAAATATCACACTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAG"))
        {
            cout << "Read: "<<_readId<<endl;
            exit(1);
        }*/
        /*std::cout << "Head: " << seq_head.str() << "\nCola: " << seq_tail.str() << "\n";
        std::cout << "Secuencia Recuperada: " << full_path.str() << " " << "\n";
        std::cout << "Secuencia Original: " << _seq.str() << "\n";
        cout << " Lectura: " << _readId << "\n";
        exit(1);*/
        /*
         * PathGraph "Visualization"
         */
        //path_graph.show();
        /*if (_seq.str() != full_path.str())
            std::cout <<"\n"<< _seq.str() << "\n"
                      << full_path.str() << "\n";*/
        return full_path;
    }
    return DnaSequence("");
}

//ReadsCorrector
template<>
void ReadCorrector<false>::correct_reads() {
    cout << "STAGE: Reads Correction: "<<_sc.size()<<"\n";
    Progress::get().size_total = _sc.size();
    size_t cont = 0;
    #pragma omp parallel
    {
        #pragma omp single
        for (auto &read: _sc.getIndex())
        {
            #pragma omp task shared(read)
            {
                /*if (read.first.getId() != 27472)
                {
                    continue;
                }
                cout << "Read: " << read.first.getId() << "\n";
                cout << "Sequence: " << read.second.sequence.str() << "\n";
                if (!(read.first.getId() % 10000))
                    cout << "Read: "<<read.first.getId()<<"\n";*/
                //cout << "Read: "<<read.first.getId()<<endl;
                //Progress::update(read.first.getId());
                Progress::update(cont++);
                if (read.first.getId() % 2 == 0)
                {
                    PathContainer<false> pc(read.first,(&_dbg),read.second.sequence);
                    DnaSequence full_seq = pc.correct_read();
                    if (full_seq.length() != 0) {
                       /* if (read.second.sequence.length() > (full_seq.length()+10)) {
                            cout << "ReadId: "<<read.first.getId()<<endl;
                            cout << "Sequence: " << read.second.sequence.str() << endl;
                            cout << "NewSequence: " << full_seq.str() << endl;
                            exit(1);
                        }*/
                        _sc.setRead(read.first.getId(), full_seq);
                        //_sc.setRead(read.first.rc().getId(), full_seq.complement());
                    }
                }
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
        ,const DBG<true> * dbg
        ,size_t * score_ed
        ,char * expected_head
        ,size_t &branches
        ,size_t * distance
        ,bool behave)
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
        ,const DBG<true> * dbg
        ,size_t * score_ed
        ,char *expected_path
        ,size_t & branches
        ,char * alternative_path
        ,size_t * distance
        ,size_t * alternative_ed)
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