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

    NaiveDBG<false> graph(sc);
}
