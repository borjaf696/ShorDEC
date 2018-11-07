#include <getopt.h>
#include <iostream>
#include <unordered_set>
#include <vector>
#include "../ReadData/sequence_container.h"
#include "../DBG/path.h"

auto print_blanks = [](size_t num_blanks)
{
    for (size_t i = 0; i < num_blanks; ++i)
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
};
auto print_use = [](char ** argc)
{printf("Usage: %s "
                "-f [path_to_file][path_to_dir] -k [kmer_size] -u [path_to_unitigs.gfa] -p [pair_end] -t [num_threads]"
        ,argc[0]);};

bool parse_args(int argv, char **argc, std::string& path_to_file, std::string& path_to_write
        ,std::string& path_unitigs, std::string& program, bool &pair_end, bool & thirdPartyCount)
{
    int opt = 0;
    char optString[] = "f:k:o:u:h:t:r:c:p";
    std::vector<bool> mandatory(3,false);
    while ((opt = getopt(argv,argc,optString))!=-1){
        switch (opt)
        {
            case 'f':
                path_to_file = optarg;
                mandatory[0] = true;
                break;
            case 'k':
                Parameters::get().kmerSize = atoi(optarg);
                mandatory[1] = true;
                break;
            case 'h':
                Parameters::get().accumulative_h = atoi(optarg);
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
            case 'c':
                thirdPartyCount = true;
                program = optarg;
                break;
            case 'r':
                Parameters::get().missmatches = atof(optarg);
                mandatory[2] = true;
                break;
            case 't':
                break;
        }
    }
    std::cout << "Number args: "<<argv<<"\n";
    for (auto val : mandatory)
        if (!val){
            print_use(argc);
            return false;
        }
    return true;

}

int main(int argv, char ** argc){

	std::string path_to_file(argc[1]),path_to_write, path_unitigs, program;
    bool pair_end = false, thirdPartyCount = false;
    if (!parse_args(argv,argc,path_to_file,path_to_write, path_unitigs, program, pair_end, thirdPartyCount))
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

    size_t kmer_sizes[] = {1,2,2};
    /*
     * Iteratively we are going to correct the reads
     */
    NaiveDBG<false> naiveDBG(sc, thirdPartyCount, path_to_file, program);
    std::cout << "Number of Reads: " << sc.getIndex().size() << "\n";
    for (auto i: kmer_sizes) {
        Parameters::get().kmerSize = Parameters::get().kmerSize * i;
        std::cout << "Building DBG\n Kmer Size: " << Parameters::get().kmerSize << "\n";
        if (i > 1)
            naiveDBG = NaiveDBG<false>(sc,thirdPartyCount, path_to_write, program);
        listDBG<false> listDBG(&naiveDBG);
        ReadCorrector<false> read(sc, listDBG);
        std::cout << "Writing new reads in: " << path_to_write << "\n";
        sc.writeSequenceContainer(path_to_write);
        print_blanks(5);
        //sc.ShowInfo();
    }
    /*
     * Depending on pair_end info one approach or another
     */
    std::cout << "Creating unitigs and writing unitigs: " << path_unitigs << "\n";
    if (pair_end)
    {
        DBG<true> * dbg = new NaiveDBG<true>(sc,thirdPartyCount, path_to_write, program);
        boostDBG<true> boostDBG1(dbg);
        boostDBG1.ProcessTigs(path_unitigs);
        exit(1);
    }
    NaiveDBG<false> dbg = NaiveDBG<false>(sc,thirdPartyCount, path_to_write, program);
    listDBG<false> listDBG(&dbg);
    listDBG.ProcessTigs(path_unitigs);
    //dbg.ProcessTigs(path_unitigs);
}
