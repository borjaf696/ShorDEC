#include <getopt.h>
#include <iostream>
#include <unordered_set>
#include "sequence_container.h"
#include "../DBG/path.h"
auto print_use = [](char ** argc)
{printf("Usage: %s "
                "-f [path_to_file][path_to_dir] -k [kmer_size] -t [num_threads]\n"
        ,argc[0]);};

bool parse_args(int argv, char **argc, std::string& path_to_file, std::string& path_to_write)
{
    int opt = 0;
    char optString[] = "f:k:o:t";
    while ((opt = getopt(argv,argc,optString))!=-1){
        switch (opt)
        {
            case 'f':
                path_to_file = optarg;
                break;
            case 'k':
                Parameters::get().kmerSize = atoi(optarg);
                //Cambiar
                Parameters::get().accumulative_h = 100;
                break;
            case 'o':
                path_to_write = optarg;
                break;
            case 't':
                break;
        }
    }
    if (argv < 5) {
        print_use(argc);
        return false;
    }
    return true;

}

int main(int argv, char ** argc){

	std::string path_to_file(argc[1]),path_to_write;
    if (!parse_args(argv,argc,path_to_file,path_to_write))
        exit(0);
	std::cout << path_to_file << "\n";
	SequenceContainer sc;
	sc.loadFromFile(path_to_file);

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

    size_t kmer_sizes[] = {1,2,2};
    /*
     * Iteratively we are going to correct the reads
     */
    std::cout << "Size sequence container: "<<sc.getIndex().size() << "\n";
    for (auto i: kmer_sizes) {
        Parameters::get().kmerSize = Parameters::get().kmerSize * i;
        std::cout << "Building DBG, Kmer Size: "<<Parameters::get().kmerSize << "\n";
        NaiveDBG naiveDBG(sc);
        //naiveDBG.show_info();
        ReadCorrector read(sc, naiveDBG);
    }
    std::cout <<"Size sequence container: " <<sc.getIndex().size() << "\n";
    std::cout << "Writing new reads in: "<<path_to_write << "\n";
    sc.writeSequenceContainer(path_to_write);
    std::cout << "Creating unitigs\n";
}
