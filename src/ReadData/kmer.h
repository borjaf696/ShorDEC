#include <unordered_map>
#include <memory>

#include "sequence_container.h"

class Kmer 
{
public:
	Kmer(){}
	~Kmer(){}
	Kmer(const DnaSequence&,size_t,size_t);
	Kmer(DnaSequence seq):_seq(seq){}
	Kmer(const Kmer & kmer){
		_seq = kmer.getSeq();
	}
	Kmer(const Kmer *kmer){
		_seq = kmer->getSeq();
	}

	Kmer(const std::string &string_own)
	{
		_seq = DnaSequence(string_own);
	}

	DnaSequence rc(){
		return _seq.complement();
	}

	void appendRight(DnaSequence::NuclType symbol);
	void appendLeft(DnaSequence::NuclType symbol);

	const DnaSequence getSeq() const {
		return _seq;
	}

	const DnaSequence* getSeq_ref() const{
		return &_seq;
	}

	DnaSequence::NuclType at(size_t index) const;

	//Chequear con algo mas de calma
	bool operator==(const Kmer&) const;
	bool operator!=(const Kmer&) const;
	Kmer& operator=(const Kmer&);

	std::size_t hash() const{
		return std::hash<std::string>()(_seq.str());
	}
	std::string str() const
	{
		return _seq.str();
	}
private:
	DnaSequence _seq;
};

namespace std{
	template <>
	struct hash<Kmer>{
		std::size_t operator()(const Kmer& kmer) const{
			return kmer.hash();
		}
	};
};
struct KmerInfo{
	KmerInfo(Kmer kmer1, size_t pos)
			:kmer(kmer1),kmer_pos(pos)
	{}
	Kmer kmer;
	size_t kmer_pos;
};

class KmerIt{
public:
	typedef std::forward_iterator_tag iterator_category;
	KmerIt(const DnaSequence*,size_t);

	bool operator==(const KmerIt&) const;
	bool operator!=(const KmerIt&) const;

	//Return the kmer which is being read right now
	KmerInfo operator*() const;
	KmerIt& operator++();
protected:
	const DnaSequence* _own_seq;
	size_t _pos;
	Kmer _kmer;
};

class IterKmers
{
public:
	IterKmers(const DnaSequence& seq):_seq(seq){}
	KmerIt begin();
	KmerIt end();

private:
	const DnaSequence& _seq;
};