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

void iteraKmers(const SequenceContainer & sc)
{
    std::vector<KmerInfo<false>> v_info;
    for (auto &read: sc.getIndex())
    {
        std::cout << "Read Number: "<<read.first.getId()<<"\n";
        for (auto k: IterKmers<false>(read.second.sequence))
        {
            v_info.push_back(k);
        }
    }
    std::cout << "Ahora que llame a los destructores de secuencias\n";
    v_info.clear();
}

int main(int argv, char ** argc){
    std::cout << "Test Environment\n";
    std::string path_to_file(argc[1]),path_to_write, path_unitigs, program;
    bool pair_end = false, thirdPartyCount = false;
    if (!parse_args(argv,argc,path_to_file,path_to_write, path_unitigs, program, pair_end, thirdPartyCount))
        exit(0);
	std::cout << path_to_file << " "<<path_unitigs<<"\n";
	SequenceContainer sc;
    sc.load(path_to_file, pair_end);
    //NaiveDBG<false> graph(sc);
    DBG<true> * graph = new NaiveDBG<true>(sc, thirdPartyCount,path_to_file, program);
    //Lets try boost
    boostDBG<true> boostDBG1(graph);
    boostDBG1.ProcessTigs(path_unitigs);
    std::cout << "END!\n";
}
