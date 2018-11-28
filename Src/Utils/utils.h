#include <fstream>
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
        for (auto k:histogram)
        {
            if (!h)
                kmersAccumulated+=k;
            else
                kmersAccumulated+=(k*h);
            if ((float)kmersAccumulated >= threshold)
                break;
            h++;
        }
        return h;
    }
    double missmatches;
    size_t accumulative_h = 0;
    size_t kmerSize;
    size_t numThreads;
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
};

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
bool isSubset(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
    for (auto el: set1)
        if (set2.find(el) == set2.end())
            return false;
    return true;
}

template<typename T>
bool isSame(const unordered_set<T> & set1, const unordered_set<T> & set2)
{
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

/*
 * Translate
 */

/*
 * Common graph operations
 */
/*
 * Check -> maximal_clique (using heuristic approach)
 */
template<typename G, typename Vt, typename Vi>
std::priority_queue<pair<size_t,vector<Vt>>> findMaxClique(const G & graph, std::set<size_t> & idCliques) {
    if (!boost::num_edges(graph))
        return std::priority_queue<pair<size_t,vector<Vt>>>();
    std::priority_queue<pair<size_t,vector<Vt>>> endCliques;
    std::vector<Vt> maxClique;
    std::vector<Vt> tmpClique;
    /*
     * Sum -> shift id to the left. Example:
     *      - 0 1 4 -> 1 0 0 1 1 -> 19
     *      - 0 3 2 -> 1 1 0 0 -> 10
     * Both cases sum = 5 but id_different. Problem with cliques larger than 64 Nodes.
     */
    int64_t sum;

    auto findMaxCliqueWithVertex = [&sum](const Vt vertex, const int maxCliqueSize, const G &graph)
    {
        std::vector<Vt> clique;
        clique.emplace_back(vertex);
        sum |= 1 << graph[vertex].id;

        std::unordered_set<Vt> candidateNeighbors;

        std::unordered_set<Vt> visited;
        visited.emplace(vertex);

        typename G::adjacency_iterator adjVertex, adjVertEnd;
        for (tie(adjVertex, adjVertEnd) = boost::adjacent_vertices(vertex, graph); adjVertex != adjVertEnd; ++adjVertex) {
            candidateNeighbors.emplace(*adjVertex);
        }

        std::unordered_set<Vt> tmp;

        while (!candidateNeighbors.empty()) {
            const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph](const Vt &lhs, const Vt &rhs) {
                return boost::degree(lhs,graph) < boost::degree(rhs,graph);
            });

            const auto highestDegVert = *highestDegNeighborIt;
            clique.emplace_back(highestDegVert);
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
    Vi vertex, v_end;
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex)
    {
        sum = 0;
        tmpClique = findMaxCliqueWithVertex(*vertex, maxClique.size(), graph);
        if ((tmpClique.size() >= CLICK_LIMIT) && (idCliques.find(sum) == idCliques.end()))
        {
            idCliques.emplace(sum);
            endCliques.push(pair<size_t,vector<Vt>>(tmpClique.size(), tmpClique));
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
void createCountMapDSK(T * dbg, string path_to_file)
{
    int lineNo = 1;
    size_t count;
    Y y;
    std::string node;
    std::ifstream infile(path_to_file);
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
            dbg->insert(y);
        }
        //map[t] = pair<Y,Y>(count,count);
        lineNo++;
    }
    infile.close();
}