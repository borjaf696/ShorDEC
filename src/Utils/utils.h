#include <vector>
#include <unordered_set>
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
    }
    size_t accumulative_h;
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

/*
 * Common graph operations
 */
/*
 * Check -> maximal_clique (using heuristic approach)
 */
template<typename G, typename Vt, typename Vi>
std::vector<Vt> findMaxClique(const G & graph) {
    if (!boost::num_edges(graph))
        return std::vector<Vt>();
    std::vector<Vt> maxClique;
    std::vector<Vt> tmpClique;

    auto findMaxCliqueWithVertex = [](const Vt vertex, const int maxCliqueSize, const G &graph)
    {
        std::vector<Vt> clique;
        clique.emplace_back(vertex);

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
        if (boost::degree(*vertex, graph) >= maxClique.size()) {
            tmpClique = findMaxCliqueWithVertex(*vertex, maxClique.size(), graph);

            if (tmpClique.size() >= 2 && tmpClique.size() > maxClique.size()) {
                maxClique = std::move(tmpClique);
            }
        }
    }
    return maxClique;
}
