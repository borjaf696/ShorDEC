#include <math.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <boost/version.hpp>
#include <boost/regex.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>

#define INF 99999

#define CLICK_LIMIT 2
//Smooth factor to take into accout legitimate k-mers at erroneous positions
#define SMOOTH_FACTOR 1.0

using namespace std;
/*
 * Map -> histogram<K,Y>, histogram<K,container>
 */
template<typename T, typename Y>
vector<size_t> getHistogram(const unordered_map<T,Y> map, size_t max)
{
    vector<size_t> histogram(max,0);
    for (auto key:map)
        histogram[key.second]++;
    return histogram;
};

template<typename T, typename Y>
vector<size_t> getHistogram(const unordered_map<T,pair<Y,Y>> map, size_t max)
{
    vector<size_t> histogram(max+1,0);
    for (auto key:map)
        histogram[key.second.first]++;
    return histogram;
};
template<typename T>
vector<T> getHistogramFromFile(std::string path_to_file)
{
    int lineNo = 0;
    std::vector<T> output = {0};
    std::ifstream infile(path_to_file);
    bool prev_num = false;
    for( std::string line; getline( infile, line ); )
    {
        lineNo++;
        if (lineNo < 9)
            continue;
        //Abundance
        std::string s;
        size_t count = 0;
        while (!prev_num)
        {
            getline(infile,line);
            s = boost::regex_replace(line,boost::regex("[^0-9]*([0-9]+).*"),string("\\1"));
            if (std::isdigit(s[0]))
                prev_num = true;
            count ++;
            if (count > 5)
                return output;
        }
        //Number which matters
        prev_num = false;
        getline(infile,line);
        s = boost::regex_replace(line,boost::regex("[^0-9]*([0-9]+).*"),string("\\1"));
        if (line.back() == '\n' or line.back() == '\r')
        {
            line.pop_back();
        }
        output.push_back(stoi(s));
    }
    infile.close();
    return output;
};
/*
 * Set operations
 */
template<typename T>
void show_set(const unordered_set<T> & set1)
{
    for (auto p: set1)
        cout << " "<<p.str();
    cout << endl;
}
template<typename T>
unordered_set<T> sustract(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    unordered_set<T> s_set;
    for (auto i = set1.begin(); i != set1.end(); ++i)
        if (set2.find(*i) == set2.end())
            s_set.emplace(*i);
    return s_set;
}
template<typename T>
unordered_set<T> getIntersection(const unordered_set<T>& set_full,const unordered_set<T>& nn_pair)
{
    unordered_set<T> intersection;
    for (auto i = nn_pair.begin(); i != nn_pair.end(); i++) {
        if (set_full.find(*i) != set_full.end())
            intersection.insert(*i);
    }
    return intersection;
}

template<typename T>
unordered_set<T> getUnion(const unordered_set<T>& set1, const unordered_set<T>& set2)
{
    unordered_set<T> result = set1;
    result.insert(set2.begin(),set2.end());
    return result;
}

template<typename T>
bool in(const unordered_set<T> & set1, T el)
{
    return (set1.find(el)!=set1.end());
}

template<typename T>
bool isSubset(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    for (auto el: set1)
        if (set2.find(el) == set2.end())
            return false;
    return true;
}

template <typename T>
bool isSubset(const unordered_set<T> & set1, const unordered_set<T> & set2, double percent)
{
    int32_t num_fails = ceil(set1.size()*(1-percent)), fails = 0;
    if (set1.size() > set2.size()+num_fails)
    {
        for (auto el:set2)
            if (set1.find(el) == set1.end())
                if (++fails >= num_fails) return false;
        return true;
    }
    for (auto el: set1)
        if (set2.find(el) == set2.end())
            if (++fails>=num_fails) return false;
    return true;
}

template<typename T>
unordered_set<T> isGetSubset(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    unordered_set<T> result;
    for (auto el: set1) {
        if (set2.find(el) == set2.end())
            return unordered_set<T>();
        result.emplace(el);
    }
    return result;
}

template<typename T>
bool isSame(const unordered_set<T> & set1, const unordered_set<T> & set2, double percent)
{
    size_t max = std::max(set1.size(), set2.size());
    int32_t num_fails = ceil(max*(1-percent)), fails = 0;
    if (set1.empty() && !set2.empty())
        return false;
    if (set2.empty() && !set1.empty())
        return false;
    for (auto el:(max ==set1.size())?set1:set2)
        if ((max==set1.size())?set2.find(el) == set2.end():set1.find(el)==set2.end())
            if (++fails >= num_fails) return false;
    return true;
}

template<typename T>
bool isSame(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    if (set1.empty() && set2.empty())
        return false;
    for (auto el: set1)
    {
        if (set2.find(el)==set2.end())
            return false;
    }
    for (auto el:set2)
    {
        if (set1.find(el) == set1.end())
            return false;
    }
    return true;
}

template<typename T>
bool same(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    for (auto el: set1)
    {
        if (set2.find(el) == set2.end())
            return false;
    }
    for (auto el: set2)
    {
        if (set2.find(el) == set2.end())
            return false;
    }
    return true;
}
/*
 * Parsing
 */
struct Parse{
    /*
     * Get fasta ids from sga files
     */
    static std::unordered_set<uint32_t> getIdsFromFiles(std::string path)
    {
        std::unordered_set<uint32_t> ids;
        std::ifstream infile(path);
        for( std::string line; getline( infile, line ); )
        {
            if (line[0] == '>')
            {
                std::stringstream ss(line);
                string token;
                getline(ss,token, ' ');
                token.erase(0,1);
                ids.emplace(stoi(token));
            }
        }
        return ids;
    }
};
/*
 * System operations
 */
struct System{
    /*
     * System utils
     */
    static std::vector<string> getAllFaFqFiles(std::string path)
    {
        namespace fs = boost::filesystem;
        std::vector<string> output;
        fs::path fs_path(path);
        if (fs::is_regular_file(fs_path))
            output.push_back(path);
        else
        {
            for (auto & p : fs::directory_iterator(path))
            {
                std::ostringstream oss;
                oss << p;
                std::string converted_path = oss.str().substr(1, oss.str().size() - 2);
                output.push_back(converted_path);
            }
        }
        return output;
    }

    static void execute(std::string instruction)
    {
        if (system(instruction.c_str()))
        {
            cout << "Fail on: "<<instruction<<"\n";
            exit(1);
        }
    }

    static std::string appendFiles(std::vector<string> files, std::string newFile)
    {
        std::string instruction = "cat ";
        for (auto s: files)
            instruction += s+" ";
        instruction += ">"+newFile;
        if (system(instruction.c_str())){
            cout << "Problem executing: "<<instruction<<"\n";
            cout << "Exiting\n";
            exit(1);
        }
        return newFile;
    }
    template<typename T>
    static void write_histogram(std::unordered_map<T,std::vector<size_t>> histogram, std::string file_name)
    {
        ofstream file;
        file.open(file_name);
        for (auto k: histogram)
        {
            file << (k.first.str()+"\n");
            for (auto t: k.second)
                file << (std::to_string(t)+" ");
            file << endl;
        }
        file.close();
    }
};

/*
 * Parameters
 */
struct Parameters
{
    static Parameters& get(){
        static Parameters param;
        return param;
    }
    void set(size_t param, std::string param_read)
    {
        if (param_read == "KmerSize")
            Parameters::get().kmerSize = param;
        if (param_read == "Accumulative")
            Parameters::get().accumulative_h = param;
        if (param_read == "missmatches")
            Parameters::get().missmatches = param;
        if (param_read == "numKmers")
            Parameters::get().num_unique_kmers = param;
        if (param_read == "genome_size")
            Parameters::get().genome_size = param;
        if (param_read == "full_info")
            Parameters::get().full_info = param;
        if (param_read == "metagenomic")
            Parameters::get().metagenomic = param;
        if (param_read == "remove_duplicates")
            Parameters::get().remove_duplicates = param;
        if (param_read == "numThreads")
            Parameters::get().numThreads = param;
    }
    /*
    * Calculate H.
    */
    static size_t calculateAccumulativeParam(vector<size_t> histogram, size_t totalBases, size_t avgLength)
    {
        size_t kmersAccumulated = 0, h = 0;
        float kmer_size = (float)Parameters::get().kmerSize-1, avgLeng = (float)avgLength+1, tBases = (float)totalBases;
        //float threshold = tBases*Parameters::get().missmatches*((avgLeng-2*kmer_size)*kmer_size/avgLeng+(kmer_size*(kmer_size-1))/avgLeng)+1;
        float threshold = (tBases*Parameters::get().missmatches*(kmer_size*(avgLeng-2*kmer_size)/avgLeng + (kmer_size-1)/2*(2*kmer_size)/avgLeng))*SMOOTH_FACTOR;
        std::cout << "Fail Expected Kmers: "<<threshold<<"\n";
        size_t total_kmers = 0;
        for (auto k:histogram)
            total_kmers+=k;
        cout << "Total kmers: "<<total_kmers<<" "<<Parameters::get().genome_size<<"\n";
        for (auto k:histogram)
        {
            total_kmers-=k;
            if (!h)
                kmersAccumulated+=k;
            else
                kmersAccumulated+=(k*h);
            if (((float)kmersAccumulated > threshold) || (total_kmers<Parameters::get().genome_size))
            {
                break;
            }
            h++;
        }
        cout << "KmersErased: "<<kmersAccumulated<<endl;
        cout << "Remaining k-mers: "<<total_kmers<<endl;
        return h;
    }

    static void reestimateParams(std::string third, std::string path_to_file)
    {
        vector<string> files;
        for (auto f: System::getAllFaFqFiles(path_to_file))
            files.push_back(f);
        if (files.size() > 1)
            path_to_file = System::appendFiles(files, "newFile_tmp.fasta");
        else if (files.size() == 1)
            path_to_file = files[0];
        string instruction = "";
        if (third != "dsk")
        {
            cout << "Option only available with \"DSK\"\n";
            exit(1);
        }
        else
            instruction += "bash -c \"./Utils/script/dsk_script ";
        instruction += path_to_file+" ";
        instruction+= to_string(Parameters::get().kmerSize)+" ";
        instruction+="output output.txt >/dev/null 2>&1\"";
        if (system(instruction.c_str()))
        {
            cout << "Fail on: "<<instruction<<"\n";
            exit(1);
        }
        std::vector<size_t> histogram = getHistogramFromFile<size_t>("histogram.txt");
        size_t total_kmer = 0;
        for (auto k:histogram)
            total_kmer+=k;
        cout << "First val: "<<Parameters::get().num_unique_kmers<<"; Total_kmers: "<<total_kmer<<"\n";
        Parameters::get().num_unique_kmers -= total_kmer;
        cout << "Fails - total: "<<Parameters::get().num_unique_kmers<<"\n";
        Parameters::get().num_unique_kmers /= (Parameters::get().kmerSize-1)/2;
        cout << "Fails corrected: "<<Parameters::get().num_unique_kmers<<"\n";
    }
    double missmatches;
    size_t genome_size = 0;
    double num_unique_kmers = 0;
    size_t accumulative_h = 0;
    bool full_info = false, metagenomic = false, remove_duplicates = false;
    size_t kmerSize;
    size_t numThreads;
    bool show = false;
};

struct Progress
{
    static void update(size_t num_actual){
        if (Progress::get().show) {
            int val=(int) (((float) num_actual) / ((float) Progress::get().size_total) * 100);
            if (val % 10 == 0)
                if (!Progress::get().f[val/10]) {
                    Progress::get().f[val/10] = true;
                    show_progress(num_actual);
                }
        }
    }

    static Progress& get(){
        static Progress progress;
        return progress;
    }
    std::vector<bool> f = std::vector<bool>(10,false);
    size_t size_total;
    bool show = false;
private:
    static void show_progress(size_t num_actual){
        size_t val = (size_t)((float)num_actual/
                              (float)Progress::get().size_total * 100);
        std::cout << val << ((val==100)?"%\n":"% ") << std::flush;
        if (val == 100)
            _prepare_next();
    }
    static void _prepare_next(){
        Progress::get().f = std::vector<bool>(10,false);
    }
};

/*
 * Show Graph
 */
template<typename G, typename Vi>
void show_graph(const G & graph)
{
    Vi vertex, v_end;
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex)
    {
        cout <<graph[*vertex].id <<"->";
        typename G::adjacency_iterator adjVertex, adjVertEnd;
        for (tie(adjVertex, adjVertEnd) = boost::adjacent_vertices(*vertex, graph); adjVertex != adjVertEnd; ++adjVertex)
            cout << graph[*adjVertex].id<<", ";
        cout << endl;
    }
    cout <<endl;
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex)
    {
        cout <<graph[*vertex].id <<"->"<<graph[*vertex].node.str();
        cout << endl;
    }

};
/*
 * Check -> maximal_clique (using heuristic approach)
 */
template<typename G, typename Vt, typename Vi>
std::priority_queue<pair<size_t,vector<Vt>>> findMaxClique(const G & graph, std::set<size_t> & idCliques,
                                                           unordered_map<Vt,size_t> & store_map ) {
    if (!boost::num_edges(graph))
        return std::priority_queue<pair<size_t,vector<Vt>>>();
    std::priority_queue<pair<size_t,vector<Vt>>> endCliques;
    std::vector<Vt> maxClique;
    std::vector<Vt> tmpClique;
    std::unordered_set<Vt> already_checked;
    /*
     * Sum -> shift id to the left. Example:
     *      - 0 1 4 -> 1 0 0 1 1 -> 19
     *      - 0 3 2 -> 1 1 0 0 -> 10
     * Both cases sum = 5 but id_different. Problem with cliques larger than 64 Nodes.
     */
    int64_t sum;
    size_t min_flow;

    auto findMaxCliqueWithVertex = [&sum, &already_checked, &min_flow](const Vt vertex, const int maxCliqueSize, const G &graph,
                                                            unordered_map<Vt,size_t> store_map, size_t * score)
    {
        std::vector<Vt> clique;
        clique.emplace_back(vertex);
        sum |= 1 << graph[vertex].id;

        std::set<Vt> candidateNeighbors;

        std::unordered_set<Vt> visited;
        visited.emplace(vertex);

        typename G::adjacency_iterator adjVertex, adjVertEnd;
        for (tie(adjVertex, adjVertEnd) = boost::adjacent_vertices(vertex, graph); adjVertex != adjVertEnd; ++adjVertex) {
            candidateNeighbors.emplace(*adjVertex);
        }

        std::set<Vt> tmp;

        while (!candidateNeighbors.empty()) {
            const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph, store_map](const Vt &lhs, const Vt &rhs) {
                if ((store_map.at(lhs)*boost::degree(lhs,graph)) == (store_map.at(rhs)*boost::degree(rhs,graph)))
                    return graph[lhs].id > graph[rhs].id;
                return (store_map.at(lhs)*boost::degree(lhs,graph)) < (store_map.at(rhs)*boost::degree(rhs,graph));
            });
            /*const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph, store_map](const Vt &lhs, const Vt &rhs) {
                if ((boost::degree(lhs,graph)) == (boost::degree(rhs,graph)))
                    return graph[lhs].id > graph[rhs].id;
                return (boost::degree(lhs,graph)) < (boost::degree(rhs,graph));
            });*/

            const auto highestDegVert = *highestDegNeighborIt;
            min_flow = (store_map.at(highestDegVert) < min_flow)?store_map.at(highestDegVert):min_flow;
            clique.emplace_back(highestDegVert);
            (*score) += store_map.at(highestDegVert);
            //already_checked.emplace(highestDegVert);
            /*
             * Questionable
             */
            sum |= 1 << graph[highestDegVert].id;
            visited.emplace(highestDegVert);

            for (tie(adjVertex, adjVertEnd) = boost::adjacent_vertices(highestDegVert, graph); adjVertex != adjVertEnd; ++adjVertex) {
                if (candidateNeighbors.find(*adjVertex)!=candidateNeighbors.end() && visited.find(*adjVertex) == visited.end()) {
                    tmp.insert(*adjVertex);
                }
            }
            candidateNeighbors = std::move(tmp);
        }
        return clique;
    };
    /*typedef std::function<bool(std::pair<Vt, size_t>, std::pair<Vt, size_t>)> Comparator;
    Comparator compFunctor =
            [](std::pair<Vt, size_t> elem1 ,std::pair<Vt, size_t> elem2)
            {
                return elem1.second > elem2.second;
            };

    // Declaring a set that will store the pairs using above comparision logic
    std::set<std::pair<Vt, size_t>, Comparator> sortedNodes(store_map.begin(), store_map.end(), compFunctor);
    Vt vertex;*/
    Vi vertex, v_end;
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex)
    //for (auto v: sortedNodes)
    {
        //vertex = v.first;
        sum = 0;
        min_flow = 0;
        size_t score = 0;
        /*
         * Coger cada nodo como nodo de arranque
         */
        if (already_checked.find(*vertex) == already_checked.end())
        {
            score += store_map.at(*vertex);
            tmpClique = findMaxCliqueWithVertex(*vertex, maxClique.size(), graph, store_map, &score);
            if ((tmpClique.size() >= CLICK_LIMIT) && (idCliques.find(sum) == idCliques.end())) {
                idCliques.emplace(sum);
                for (auto c: tmpClique)
                    store_map[c]-=min_flow;
                endCliques.push(pair < size_t, vector < Vt >> (score, tmpClique));
            }
        }
    }
    return endCliques;
}

/*
 * Process Counts:
 *      * First JellyFish output
 *      * Second DSK output
 */

template<typename T, typename Y>
size_t createCountMapJelly(unordered_map<T,pair<Y,Y>> & map, string path_to_file)
{
    int lineNo = 1;
    size_t max_freq = 0;
    Y count;
    std::string node;
    std::ifstream infile(path_to_file);
    for( std::string line; getline( infile, line ); )
    {
        if (line[0] == '>')
        {
            count = stoi(line.substr(1,line.size()-1));
            max_freq = (count > max_freq)?count:max_freq;
        }
        else
        {
            if (line.back() == '\n' or line.back() == '\r')
            {
                line.pop_back();
            }
            node = line;
            T t(node);
            map[t] = pair<Y,Y>(count,count);
        }
        lineNo++;
    }
    infile.close();
    //TODO: Change para que reciba el previous max_size y busque los menores directos de este
    return max_freq;
}

template<typename T, typename Y>
void createCountMapDSK(T * dbg, string path_to_file, vector<size_t> histogram, size_t tb, size_t avL, bool last)
{
    int lineNo = 1;
    size_t count;
    Y y;
    std::string node;
    std::ifstream infile(path_to_file);
    for( std::string line; getline( infile, line ); )
        lineNo++;
    infile.close();
    std::cout << "Previous Kmers: "<<Parameters::get().num_unique_kmers<<"\n";
    if (Parameters::get().num_unique_kmers == 0)
        Parameters::get().num_unique_kmers = lineNo-1;
    else
    {
        float totalFails = tb*Parameters::get().missmatches, basesCorrected = Parameters::get().num_unique_kmers;
        Parameters::get().missmatches = (totalFails-basesCorrected) / tb;
        /*Parameters::get().missmatches /= ((float)Parameters::get().num_unique_kmers /(float)(lineNo-1))*2;*/
        Parameters::get().num_unique_kmers = lineNo-1;
    }
    std::cout << "Number of unique Kmers: "<<Parameters::get().num_unique_kmers<<"\n";
    std::cout << "MissMatches: "<<Parameters::get().missmatches<<"\n";
    if (!last)
    {
        Parameters::get().accumulative_h = (Parameters::calculateAccumulativeParam(histogram, tb, avL)==1)?
                                           2:Parameters::calculateAccumulativeParam(histogram, tb, avL);
    }
    infile.close();
    infile.open(path_to_file);
    for( std::string line; getline( infile, line ); )
    {
        if (line.back() == '\n' or line.back() == '\r')
        {
            line.pop_back();
        }
        count = stoi(line.substr(Parameters::get().kmerSize+1,line.size()-(Parameters::get().kmerSize+1)));
        if (count >= Parameters::get().accumulative_h)
        {
            y = Y(line.substr(0,Parameters::get().kmerSize));
            dbg->insert(y, count);
        }
        //map[t] = pair<Y,Y>(count,count);
        lineNo++;
    }
    infile.close();
}
