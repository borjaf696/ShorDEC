#include <vector>
#include <queue>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stack>
#include "../DBG/PathGraph.h"

//Constants for the path extension
#define MAX_PATH_LEN 300
#define MAX_BRANCH 5
#define ERROR_RATE 0.01
#define NUM_PREV_POST_KMER 5
#define MAX_DISTANCE_READ 100
#define MAX_NUM_TRIALS 5

/*El camino se extiende desde un k-mer source (head) hasta un k-mer source (tail).
 * De todos los caminos extendidos se selecciona aquel con distancia de edicion minima con
 * la secuencia que se esta corrigiendo. Como algunos k-mers solidos pueden ser erroneos
 * se busca mas alla de dicho k-mer atendiendo a mas k-mer sources posibles.
 * */

struct stack_el{
    stack_el(Kmer kmer1, size_t pos1)
            :kmer(kmer1),pos(pos1)
            {}
    Kmer kmer;
    size_t pos;
};

class Path
{
public:
    //Constructor
    Path()
    {
        _DP = std::vector<size_t>((MAX_PATH_LEN+1)*(MAX_PATH_LEN+1),0);
        for (uint i = 0; i < (1+ERROR_RATE)*MAX_PATH_LEN+1; ++i)
            _DP[i] = i;
    }
    //Extender
    size_t extend(const DnaSequence&
            ,std::pair<Kmer,size_t>
            ,std::pair<Kmer,size_t>
            ,const DGB&
            ,size_t*
            ,char*
            ,size_t&);
private:
    std::vector<size_t> _DP;
};

class PathContainer
{
public:
    PathContainer(FastaRecord::Id id, const DGB& dbg, const DnaSequence& seq):
            _readId(id),_dbg(dbg),_seq(seq)
    {
        check_read();
    };

    void correct_read();
    size_t getSolidLength(){
        return _solid.size();
    }
    //Mover PVT
    size_t check_read();
    int check_solids(size_t,size_t,size_t,size_t
            ,size_t &,Kmer&);
private:
    std::vector<Path> _path_extended;
    std::vector<std::pair<Kmer,size_t>> _solid;
    FastaRecord::Id _readId;
    const DGB &_dbg;
    const DnaSequence &_seq;
};

class ReadCorrector
{
public:
    ReadCorrector(const SequenceContainer &sc, const DGB &dbg):
            _sc(sc),_dbg(dbg)
    {
        correct_reads();
    }

    void correct_reads();
private:
    const SequenceContainer &_sc;
    const DGB &_dbg;
};