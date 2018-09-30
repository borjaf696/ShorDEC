#include "kmer.h"
#include "../Utils/utils.h"

/*
 * Single_end kmer
 */
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

bool Kmer::operator<(const Kmer &other) const {
    return (_seq < other._seq);
}

bool Kmer::operator>(const Kmer &other) const {
    return (_seq > other._seq);
}

bool Kmer::operator==(const Kmer& other) const {
    return (_seq == other._seq);
}

bool Kmer::operator!=(const Kmer& other) const{
    return _seq != other._seq;
}

//Kmer iterator -> Single_end (false)
KmerIt<false>::KmerIt(const DnaSequence *seq, size_t pos)
        :_own_seq(seq),_pos(pos)
{
    if (pos != seq->length() - Parameters::get().kmerSize+1)
        _kmer = Kmer(*seq,0,Parameters::get().kmerSize);
}

bool KmerIt<false>::operator==(const KmerIt &kmerIt) const {
    return _own_seq == kmerIt._own_seq && _pos == kmerIt._pos;
}

bool KmerIt<false>::operator!=(const KmerIt &kmerIt) const {
    return !(*this == kmerIt);
}

KmerIt<false>& KmerIt<false>::operator++() {
    size_t newPos = _pos + Parameters::get().kmerSize;
    if (newPos < _own_seq->length())
        _kmer.appendRight(_own_seq->atRaw(newPos));
    ++_pos;
    return *this;
}

KmerInfo<false> KmerIt<false>::operator*() const {
    return KmerInfo<false>(_kmer,_pos);
}

//IterKmers
KmerIt<false> IterKmers<false>::begin() {
    if (_seq.length() < Parameters::get().kmerSize)
        return this->end();
    return KmerIt<false>(&_seq,0);
}

KmerIt<false> IterKmers<false>::end() {
    return KmerIt<false>(&_seq, _seq.length()-Parameters::get().kmerSize+1);
}

/*
 * Pair_end kmer
 */

Pair_Kmer::Pair_Kmer(const DnaSequence &ds_1, size_t s1, size_t l1,
        const DnaSequence &ds_2,size_t s2, size_t l2):exist(true) {
    _seq_left = *(ds_1.substr(s1,l1));
    _seq_right = *(ds_2.substr(s2,l2));
}

void Pair_Kmer::appendRight(DnaSequence::NuclType symbol_left, DnaSequence::NuclType symbol_right) {
    _seq_left.append_with_replace_right(symbol_left);
    _seq_right.append_with_replace_right(symbol_right);
}

void Pair_Kmer::appendLeft(DnaSequence::NuclType symbol_left, DnaSequence::NuclType symbol_right) {
    _seq_left.append_with_replace_left(symbol_left);
    _seq_right.append_with_replace_left(symbol_right);
}

pair<DnaSequence::NuclType,DnaSequence::NuclType> Pair_Kmer::at(size_t index) const
{
    return {_seq_left.atRaw(index),_seq_right.atRaw(index)};
}

Pair_Kmer& Pair_Kmer::operator=(const Pair_Kmer &other) {
    _seq_left = other._seq_left;
    _seq_right = other._seq_right;
    return *this;
}

bool Pair_Kmer::operator==(const Pair_Kmer& other) const {
    return (_seq_left == other._seq_left && _seq_right == other._seq_right);
}

bool Pair_Kmer::operator!=(const Pair_Kmer& other) const{
    return (_seq_left != other._seq_left || _seq_right != other._seq_right);
}
