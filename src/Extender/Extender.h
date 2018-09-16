#include "../DBG/DBG.h"

using namespace std;

class Extender {
public:
    //Constructor estandar
    Extender(DBG& dbg):_dbg(dbg) {};
    //Extender
    virtual void Extend(Kmer, vector<Kmer>&, bool) = 0;
protected:
    DBG& _dbg;
};

class UnitigExtender: public Extender
{
public:
    void Extend(Kmer, vector<Kmer>&, bool);
};