#include <getopt.h>
#include <iostream>
#include <unordered_set>
#include "../ReadData/sequence_container.h"
#include "../DBG/path.h"

auto print_use = [](char ** argc)
{printf("Usage: %s "
                "-f [path_to_file][path_to_dir] -k [kmer_size] -u [path_to_unitigs.gfa] -p [pair_end] -t [num_threads]\n"
        ,argc[0]);};

bool parse_args(int argv, char **argc, std::string& path_to_file, std::string& path_to_dir, std::string& path_to_write
        ,std::string& path_unitigs, std::string& program, bool &pair_end, bool & thirdPartyCount)
{
    int opt = 0;
    char optString[] = "s:k:o:u:h:t:r:c:g:p:";
    std::vector<bool> mandatory(3,false);
    while ((opt = getopt(argv,argc,optString))!=-1){
        switch (opt)
        {
            case 's':
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
                mandatory[0] = true;
                path_to_dir = optarg;
                break;
            case 'c':
                thirdPartyCount = true;
                program = optarg;
                break;
            case 'r':
                Parameters::get().missmatches = atof(optarg);
                mandatory[2] = true;
                break;
            case 'g':
                Parameters::get().genome_size = atoi(optarg);
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
    v_info.clear();
}

int main(int argv, char ** argc){
    std::cout << "Test Environment\n";
    std::string path_to_file(""),dir_pairs(""), path_to_write(""), path_unitigs(""), program("");
    bool pair_end = false, thirdPartyCount = false;
    if (!parse_args(argv,argc,path_to_file,dir_pairs,path_to_write, path_unitigs, program, pair_end, thirdPartyCount))
        exit(0);
    SequenceContainer sc_single, sc_paired;
    if (path_to_file != "")
    {
        std::cout <<"SingleEnd: "<<path_to_file<<"\n";
        sc_single.load(path_to_file, false);
    }
    if (dir_pairs != "")
    {
        std::cout <<"PairedDir: "<<dir_pairs<<"\n";
        sc_paired.load(dir_pairs, true);
    }
    //NaiveDBG<false> graph(sc);
    if (pair_end)
    {
        DBG<true> * graph = new NaiveDBG<true>(sc_single, sc_paired, thirdPartyCount,path_to_file, dir_pairs, program, false);
        //Lets try boost
        sc_single.clear();
        sc_paired.clear();
        boostDBG<true> boostDBG1(graph);
        boostDBG1.ProcessTigs(path_unitigs);
        exit(1);
    }
    NaiveDBG<false> * graph_se = new NaiveDBG<false>(sc_single, sc_paired,thirdPartyCount, path_to_write, dir_pairs, program, true);
    listDBG<false> listDBG(graph_se);
    //Clean naiveDBG
    listDBG.ProcessTigs(path_unitigs);
    //Clean listDBG
    listDBG.clear();
    std::cout << "END!\n";
}
