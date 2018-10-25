#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <boost/version.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>

#define INF 99999
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
    static size_t calculateAccumulativeParam(vector<size_t> histogram, size_t totalBases)
    {
        size_t kmersAccumulated = 0, h = 0;
        size_t threshold = totalBases*Parameters::get().missmatches*(Parameters::get().kmerSize-1)+1;
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

/*
 * Set operations
 */
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
std::vector<vector<Vt>> findMaxClique(const G & graph) {
    if (!boost::num_edges(graph))
        return std::vector<vector<Vt>>();
    std::vector<vector<Vt>> endCliques;
    std::vector<Vt> maxClique;
    std::vector<Vt> tmpClique;
    std::set<size_t> idCliques;
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
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex) {
        sum = 0;
        tmpClique = findMaxCliqueWithVertex(*vertex, maxClique.size(), graph);
        if (tmpClique.size() >= 2 && tmpClique.size() > maxClique.size()) {
            endCliques.clear();
            idCliques.clear();
            maxClique = std::move(tmpClique);
            endCliques.push_back(maxClique);
            idCliques.emplace(sum);
        }else if (tmpClique.size() == maxClique.size() && idCliques.find(sum) == idCliques.end())
        {
            idCliques.emplace(sum);
            endCliques.push_back(tmpClique);
        }
    }
    return endCliques;
}
