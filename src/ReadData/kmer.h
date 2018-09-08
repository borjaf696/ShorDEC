#include <unordered_map>
#include <memory>

#include "sequence_container.h"

class Kmer 
{
public:
	bool exist = true;
	Kmer():exist(false){}
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

	DnaSequence substr(size_t start, size_t end) const
	{
		return _seq.substr(start,end);
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
    KmerInfo()
	{}
	KmerInfo(Kmer kmer1, size_t pos)
			:kmer(kmer1),kmer_pos(pos)
	{}
	size_t hash() const
	{
		return std::hash<std::string>()(kmer.str()+std::to_string(kmer_pos));
	}
	KmerInfo(const KmerInfo & k_info)
			:kmer(k_info.kmer),kmer_pos(k_info.kmer_pos)
	{
	}
    KmerInfo& operator=(const KmerInfo& k_info)
    {
        kmer = k_info.kmer;
        kmer_pos = k_info.kmer_pos;
        return *this;
    }
	bool operator==(const KmerInfo & k_info) const
	{
		return(kmer == k_info.kmer) && (kmer_pos == k_info.kmer_pos);
	}
	bool operator!=(const KmerInfo & k_info) const
	{
	    return !((kmer == k_info.kmer) && (kmer_pos == k_info.kmer_pos));
	}
	Kmer kmer;
	size_t kmer_pos;
};

namespace std
{
	template <>
	struct hash<KmerInfo>
	{
		size_t operator()(const KmerInfo & k_info) const
		{
			return k_info.hash();
		}
	};
}

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