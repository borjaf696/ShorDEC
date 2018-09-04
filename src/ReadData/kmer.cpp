#include "kmer.h"
#include "../Utils/utils.h"

//Kmer
Kmer::Kmer(const DnaSequence &ds, size_t start, size_t length):exist(true) {
    _seq = *(ds.substr(start,length));
}

void Kmer::appendRight(DnaSequence::NuclType symbol) {
    _seq.append_with_replace_right(symbol);
}

void Kmer::appendLeft(DnaSequence::NuclType symbol) {
    _seq.append_with_replace_left(symbol);
}

DnaSequence::NuclType Kmer::at(size_t index) const
{
    return _seq.atRaw(index);
}

Kmer& Kmer::operator=(const Kmer &other) {
    _seq = other._seq;
    return *this;
}

bool Kmer::operator==(const Kmer& other) const {
    return (_seq == other._seq);
}

bool Kmer::operator!=(const Kmer& other) const{
    return _seq != other._seq;
}

//Kmer iterator
KmerIt::KmerIt(const DnaSequence *seq, size_t pos)
        :_own_seq(seq),_pos(pos)
{
    if (pos != seq->length() - Parameters::get().kmerSize )
        _kmer = Kmer(*seq,0,Parameters::get().kmerSize);
}

bool KmerIt::operator==(const KmerIt &kmerIt) const {
    return _own_seq == kmerIt._own_seq && _pos == kmerIt._pos;
}

bool KmerIt::operator!=(const KmerIt &kmerIt) const {
    return !(*this == kmerIt);
}

KmerIt& KmerIt::operator++() {
    size_t newPos = _pos + Parameters::get().kmerSize;
    _kmer.appendRight(_own_seq->atRaw(newPos));
    ++_pos;
    return *this;
}

KmerInfo KmerIt::operator*() const {
    return KmerInfo(_kmer,_pos);
}

//IterKmers
KmerIt IterKmers::begin() {
    if (_seq.length() < Parameters::get().kmerSize)
        return this->end();
    return KmerIt(&_seq,0);
}

KmerIt IterKmers::end() {
    return KmerIt(&_seq, _seq.length()-Parameters::get().kmerSize);
}