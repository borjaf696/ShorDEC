#include <vector>
#include <queue>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stack>
#include "../DBG/PathGraph.h"
#include <omp.h>

//Constants for the path extension
#define MAX_PATH_LEN 300
#define EXTRA_LEN 1.15
#define MAX_BRANCH 5
#define ERROR_RATE 0.01
#define NUM_PREV_POST_KMER 5
#define MAX_DISTANCE_READ 100
#define MAX_NUM_TRIALS 3

/*El camino se extiende desde un k-mer source (head) hasta un k-mer source (tail).
 * De todos los caminos extendidos se selecciona aquel con distancia de edicion minima con
 * la secuencia que se esta corrigiendo. Como algunos k-mers solidos pueden ser erroneos
 * se busca mas alla de dicho k-mer atendiendo a mas k-mer sources posibles.
 * TODO: Change name because path is not exactly what is this.
 * */

using namespace std;

struct stack_el{
    stack_el(Kmer kmer1, size_t pos1, DnaSequence::NuclType nuc1)
            :kmer(kmer1),pos(pos1),nuc(nuc1)
            {}
    Kmer kmer;
    size_t pos;
    DnaSequence::NuclType nuc;
};
template<bool P>
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
    ~Path()
    {
        _DP.clear();
    }
    //Extender
    size_t extend(const DnaSequence&
            ,KmerInfo<P>
            ,KmerInfo<P>
            ,const DBG<P> *
            ,size_t*
            ,char*
            ,size_t&);
    size_t extend_head(const DnaSequence&
            ,KmerInfo<P>
            ,const DBG<P>&
            ,size_t*
            ,size_t*
            ,char*
            ,size_t&
            ,bool
            ,KmerInfo<P>&);
private:
    std::vector<size_t> _DP;
};

template<bool P>
class PathContainer
{
public:
    PathContainer(FastaRecord::Id id, const DBG<P> * dbg, DnaSequence seq):
            _readId(id),_dbg(dbg),_seq(seq)
    {
        check_read();
    };

    ~PathContainer()
    {
        _path_extended.clear();
        _solid.clear();
        _solid.shrink_to_fit();
    }

    DnaSequence correct_read();
    size_t getSolidLength(){
        return _solid.size();
    }
    //Mover PVT
    size_t check_read();
    int check_solids(size_t,size_t,size_t,size_t
            ,size_t &,Kmer&);
private:
    std::vector<Path<P>> _path_extended;
    std::vector<KmerInfo<P>> _solid;
    FastaRecord::Id _readId;
    const DBG<P> * _dbg;
    DnaSequence _seq;
};

template<bool P>
class ReadCorrector
{
public:
    ReadCorrector(SequenceContainer &sc, const DBG<P> &dbg):
            _sc(sc),_dbg(dbg)
    {
        correct_reads();
    }

    void correct_reads();
private:
    SequenceContainer &_sc;
    const DBG<P> &_dbg;
};