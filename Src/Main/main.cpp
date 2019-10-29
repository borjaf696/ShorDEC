#include <getopt.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <boost/config.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
#include "../Utils/OptionPrinter.hpp"
#include "../ReadData/sequence_container.h"
#include "../DBG/path.h"

namespace po = boost::program_options;
namespace
{
    const size_t HELP = 2;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
}

auto print_blanks = [](size_t num_blanks)
{
    for (size_t i = 0; i < num_blanks; ++i)
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
};
auto print_use = [](char ** argc)
{printf("Usage: %s "
                "-s [single_end_file] -p [paired_end_file] -k [kmer_size] -u [path_to_unitigs.gfa] -o [output file]"
                " -c [dsk/jelly] -b [Do_Error_Correction (Default - FALSE)] -h [Do_polish (Default - FALSE)] -d [0 (ReverseComplement paired end/otherwise]"
                "-f [use bot forward and reverse hits (Default - False) -t [Number of threads (not fully parallelized)] -n [Remove duplicates (for speed)] \n"
        ,argc[0]);};

bool parse_args(int argv, char **argc, std::string& path_to_file, std::string& dir_pairs, std::string& path_to_write
        ,std::string& path_unitigs, std::string& program, bool &pair_end, bool & thirdPartyCount, bool & do_correction
        ,bool & do_polish, bool & meta)
{
    int opt = 0;
    char optString[] = "s:k:o:u:h:t:r:c:g:p:bhfmn";
    std::vector<bool> mandatory(3,false);
    while ((opt = getopt(argv,argc,optString))!=-1){
        switch (opt)
        {
            case 'm':
                meta = true;
                break;
            case 'b':
                do_correction = true;
                break;
            case 's':
                path_to_file = optarg;
                mandatory[0] = true;
                break;
            case 'n':
                Parameters::get().remove_duplicates = true;
                break;
            case 'f':
                Parameters::get().full_info = true;
                break;
            case 'k':
                Parameters::get().kmerSize = atoi(optarg);
                mandatory[1] = true;
                break;
            case 'h':
                do_polish = true;
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
                dir_pairs = optarg;
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
                Parameters::get().numThreads = atoi(optarg);
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

int parse_args_boost(int argc, char **argv, std::string * path_to_file, std::string * dir_pairs, std::string * path_to_write
        ,std::string * path_unitigs, std::string * program, std::string * reference, std::string * metaquastPath
        , bool * pair_end, bool * thirdPartyCount, bool * do_correction, bool * do_polish, bool * meta)
{
    po::options_description des("Options");
    des.add_options()
            ("help,h","Help message")
            ("meta,m",po::bool_switch(meta),"Metagenomic analysis")
            ("correct,b",po::bool_switch(do_correction),"Do correction step (high impact in time)")
            ("single,s",po::value<std::string>(path_to_file),"Path to single-end file")
            ("paired,p",po::value<std::string>(dir_pairs)->required(),"Path to directory with paired-end reads")
            ("duplicated,n",po::bool_switch(&(Parameters::get().remove_duplicates)),"Remove duplicated reads")
            ("fullInfo,f",po::bool_switch(&(Parameters::get().full_info)),"Use all information available")
            ("kmerSize,k",po::value<size_t>(&(Parameters::get().kmerSize))->required(),"Size of the k-mers")
            ("polish,h",po::bool_switch(&(Parameters::get().polish)),"Do polishing step (high impact in time)")
            ("output,o",po::value<std::string>(path_to_write)->required(),"Output sequence file")
            ("unitigs,u",po::value<std::string>(path_unitigs)->required(),"Output unitigs file")
            ("third,c",po::value<std::string>(program),"Third party to use on k-mer counting (dsk)")
            ("error,r",po::value<double>(&(Parameters::get().missmatches)),"Estimated error ratio")
            (",t",po::value<size_t>(&(Parameters::get().numThreads)),"Number of threads to use")
            ("postprocess,",po::bool_switch(&(Parameters::get().postProcess)),"Post process output removing full duplicates")
            (",g",po::value<size_t>(&(Parameters::get().genome_size))," Estimated number of nodes (only for testing purposes)")
            ("outputfmt,",po::bool_switch(&(Parameters::get().gfa)),"GFA format activate")
            ("reference,",po::value<std::string>(reference),"Using this reference for metaquast")
            ("metaquastpath,",po::value<std::string>(metaquastPath),"Metaquast path (default path ../quast/metaquast.py)");
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,argv).options(des).run(),vm);
        if ( vm.count("help")  )
        {
            std::cout << "viaDBG" << endl << endl;
            rad::OptionPrinter::printStandardAppDesc("viaDBG",
                                                     cout,
                                                     des);
            return HELP;
        }
        po::notify(vm);
    }catch(boost::program_options::required_option& e)
    {
        rad::OptionPrinter::formatRequiredOptionError(e);
        std::cerr << "ERROR: " << e.what() << endl << endl;
        rad::OptionPrinter::printStandardAppDesc("viaDBG",
                                                 cout,
                                                 des);
        return ERROR_IN_COMMAND_LINE;
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << "ERROR: " << e.what() << endl << endl;
        rad::OptionPrinter::printStandardAppDesc("viaDBG",
                                                 cout,
                                                 des);
        return ERROR_IN_COMMAND_LINE;
    }
    return SUCCESS;
}

void basicReport(std::string pathFile, std::string dirPairs, bool doCorrection)
{
    std::cout << "**************************************CORE PARAMETERS*******************************"<<std::endl;
    std::cout << "Single end file: "<<pathFile<<std::endl;
    std::cout << "Paired end dir: "<<dirPairs<<std::endl;
    std::cout << "Correction: "<<doCorrection<<std::endl;
    std::cout << "Polishing: "<<Parameters::get().polish<<std::endl;
    std::cout << "Remove duplicates: "<<Parameters::get().remove_duplicates<<std::endl;
    std::cout << "Use revcomp information: "<<Parameters::get().full_info<<std::endl;
    std::cout << "Missmatches: "<<Parameters::get().missmatches<<std::endl;
    std::cout << "PostProcessing: "<<Parameters::get().postProcess<<std::endl;
    std::cout << "************************************************************************************"<<std::endl;
}

int main(int argv, char ** argc){
	std::string path_to_file(""), dir_pairs(""), output_path, path_unitigs, program="dsk", reference = ""
	        , metaquastPath = '../quast/metaquast.py';
    bool pair_end = false, thirdPartyCount = true, do_correction = false, do_polish = false, meta = false;
    /*if (!parse_args(argv,argc,path_to_file,dir_pairs,
                    output_path, path_unitigs, program, pair_end, thirdPartyCount, do_correction, do_polish, meta))
        exit(0);*/
    if (parse_args_boost(argv,argc, &path_to_file, &dir_pairs, &output_path,&path_unitigs,&program, &reference,
            metaquastPath, &pair_end, &thirdPartyCount,&do_correction,&do_polish,&meta))
        exit(0);
    basicReport(path_to_file, dir_pairs, do_correction);
    //cout << "Correct: "<<do_correction<<endl<<" Polish: "<<do_polish<<endl<<" Meta: "<<meta<<endl<<" ThirdParty: "<<
    SequenceContainer sc_single, sc_paired;
    pair_end = (dir_pairs!="");
    if (path_to_file != "")
    {
        std::cout <<"SingleEnd: "<<path_to_file<<"\n";
        sc_single.load(path_to_file, false);
        std::cout << "Number of single-end reads: " << sc_single.getIndex().size() << "\n";
    }
    if (dir_pairs != "")
    {
        std::cout <<"PairedDir: "<<dir_pairs<<"\n";
        sc_paired.load(dir_pairs, true);
        std::cout << "Number of paired-end reads: "<< sc_paired.getIndex().size()<<"\n";
    }
    if (meta)
    {
        boostDBG<true> boostDBG(path_to_file, dir_pairs, &sc_paired);
        exit(1);
    }
    if (do_correction)
    {
        double kmer_sizes[] = {1, 2, 2};
        /*
         * Iteratively we are going to correct the reads
         */
        NaiveDBG<false> naiveDBG(sc_single, sc_paired, thirdPartyCount, path_to_file, dir_pairs, program, false);
        for (auto i: kmer_sizes) {
            Parameters::get().kmerSize = std::min(Parameters::get().kmerSize * i, (double) 120);
            std::cout << "Building DBG\n Kmer Size: " << Parameters::get().kmerSize << "\n";
            if (i > 1)
                naiveDBG = NaiveDBG<false>(sc_single, sc_paired, thirdPartyCount, output_path, "", program, false);
            listDBG<false> listDBG(&naiveDBG);
            //Clean naiveDBG:
            naiveDBG.clear();

            if (path_to_file != "")
                ReadCorrector<false> read(sc_single, listDBG);
            if (dir_pairs != "")
                ReadCorrector<false> read(sc_paired, listDBG);
            //Clean listDBG:
            listDBG.clear();

            std::cout << "Writing new reads in: " << output_path << "\n";
            if (path_to_file != "")
                sc_single.writeSequenceContainer(output_path + "SequenceContainer_single.fasta");
            if (dir_pairs != "")
                sc_paired.writeSequenceContainer(output_path + "SequenceContainer_paired.fasta");
            Parameters::reestimateParams(program, output_path);
            print_blanks(5);
            //sc.ShowInfo();
        }
    }
    /*
     * Depending on pair_end info one approach or another
     */
    //Parameters::get().kmerSize = 120;
    if (pair_end)
    {
        DBG<true> * graph = new NaiveDBG<true>(sc_single, sc_paired, thirdPartyCount,path_to_file, dir_pairs, program, false);
        if (path_to_file != "")
            sc_single.clear();
        if (dir_pairs != "")
            sc_paired.clear();
        boostDBG<true> boostDBG1(graph);
        std::cout << "Writing unitigs(core): " << path_unitigs << "\n";
        boostDBG1.ProcessTigs(path_unitigs);
        graph->clear();
    }else {
        NaiveDBG<false> naiveDBG = NaiveDBG<false>(sc_single, sc_paired, thirdPartyCount, output_path, "", program, true);
        listDBG<false> listDBG(&naiveDBG);
        //Clean naiveDBG
        naiveDBG.clear();

        std::cout << "Creating unitigs and writing unitigs: " << path_unitigs << "\n";
        listDBG.ProcessTigs(path_unitigs);
        //Clean listDBG
        listDBG.clear();
    }

    /*
     * Aligning reference
     */
    if (reference != "")
    {
        /*
         * Alinear!
         */
    }
    cout << "Thanks for trying our tool!"<<endl;

    //dbg.ProcessTigs(path_unitigs);
}
