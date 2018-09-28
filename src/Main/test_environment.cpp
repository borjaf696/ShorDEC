#include <getopt.h>
#include <iostream>
#include <unordered_set>
#include "../ReadData/sequence_container.h"
#include "../DBG/path.h"

auto print_use = [](char ** argc)
{printf("Usage: %s "
                "-f [path_to_file][path_to_dir] -k [kmer_size] -u [path_to_unitigs.gfa] -p [pair_end] -t [num_threads]"
        ,argc[0]);};

bool parse_args(int argv, char **argc, std::string& path_to_file, std::string& path_to_write
        ,std::string& path_unitigs, bool &pair_end)
{
    int opt = 0;
    char optString[] = "f:k:o:u:t:p";
    while ((opt = getopt(argv,argc,optString))!=-1){
        switch (opt)
        {
            case 'f':
                path_to_file = optarg;
                break;
            case 'k':
                Parameters::get().kmerSize = atoi(optarg);
                //Cambiar
                Parameters::get().accumulative_h = 4;
                break;
            case 'o':
                path_to_write = optarg;
                break;
            case 'u':
                path_unitigs = optarg;
                break;
            case 'p':
                pair_end = true;
                break;
            case 't':
                break;
        }
    }
    std::cout << "Number args: "<<argv<<"\n";
    if (argv < 8) {
        print_use(argc);
        return false;
    }
    return true;

}

int main(int argv, char ** argc){

	std::string path_to_file(argc[1]),path_to_write, path_unitigs;
    bool pair_end = false;
    if (!parse_args(argv,argc,path_to_file,path_to_write, path_unitigs, pair_end))
        exit(0);
	std::cout << path_to_file << "\n";
	SequenceContainer sc;
    sc.load(path_to_file, pair_end);


	//sc.loadFromFile(path_to_file);

    /*sc.getSeq(FastaRecord::Id(0)).set(3,11);
    std::cout << sc.getSeq(FastaRecord::Id(0).getId()).str() << "\n";
    std::cout << "Kmersize "<< Parameters::get().kmerSize << "\n";
    Kmer kmer(sc.getSeq(FastaRecord::Id(0).getId()),0,10);
    std::cout << "Kmer "<<kmer.str()<<"\n";
    std::cout << kmer.at(0)<<"\n";
    for (auto kmer:IterKmers(sc.getSeq(FastaRecord::Id(0))))
        std::cout << kmer.str() << "\n";
    Kmer kmer(sc.getSeq(FastaRecord::Id(0).getId()),0,15);
    DnaSequence s1 = sc.getSeq(FastaRecord::Id(0).getId());
    std::cout << s1.str() << " \n";
    DnaSequence s2 = s1;
    s2.append_nuc_right(0);
    std::cout << s2.str() << " \n";
    Kmer kmer2 = kmer;
    std::cout << "Kmer "<< kmer.str() << "\n";
    kmer2.appendRight(0);
    std::cout << "Kmer " << kmer2.str() << "\n";
    std::unordered_set<Kmer> s;
    s.insert(kmer);
    s.insert(kmer2);
    std::unordered_set<Kmer>::const_iterator it = s.find(kmer);
    if (it != s.end())
        std::cout<< "Correcto" << "\n";
    for (auto k: s)
        std::cout << k.str() << "\n";
    Kmer kmer(sc.getSeq(FastaRecord::Id(0).getId()),0,10);
    std::vector<Kmer> k_vect;
    for (uint i = 0; i < 4 ; ++i)
    {
        Kmer kmer_aux = kmer;
        kmer_aux.appendRight(i);
        k_vect.push_back(kmer_aux);
    }
    for (uint i = 0; i < 4; ++i)
        std::cout << k_vect[i].str() << "\n";*/
}