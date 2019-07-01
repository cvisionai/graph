#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_GREEDY_ADDITIVE_HXX
#define ANDRES_GRAPH_MULTICUT_GREEDY_ADDITIVE_HXX

#include <cstddef>
#include <iterator>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>

#include "andres/partition.hxx"

namespace andres {
namespace graph {
namespace multicut {

/// Greedy agglomerative decomposition of a graph.
///
template<typename GRAPH, typename EVA, typename ELA>
void greedyAdditiveEdgeContraction(
    const GRAPH& graph,
    EVA const& edge_values,
    ELA& edge_labels
)
{
    class DynamicGraph
    {
    public:
        DynamicGraph(size_t n) :
            vertices_(n)
        {}

        bool edgeExists(size_t a, size_t b) const
        {
            return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
        }

        std::map<size_t, typename EVA::value_type> const& getAdjacentVertices(size_t v) const
        {
            return vertices_[v];
        }

        typename EVA::value_type getEdgeWeight(size_t a, size_t b) const
        {
            return vertices_[a].at(b);
        }

        void removeVertex(size_t v)
        {
            for (auto& p : vertices_[v])
                vertices_[p.first].erase(v);

            vertices_[v].clear();
        }

        void updateEdgeWeight(size_t a, size_t b, typename EVA::value_type w)
        {
            vertices_[a][b] += w;
            vertices_[b][a] += w;
        }

    private:
        std::vector<std::map<size_t, typename EVA::value_type>> vertices_;
    };

    struct Edge
    {
        Edge(size_t _a, size_t _b, typename EVA::value_type _w)
        {
            if (_a > _b)
                std::swap(_a, _b);

            a = _a;
            b = _b;

            w = _w;
        }

        size_t a;
        size_t b;
        size_t edition;
        typename EVA::value_type w;

        bool operator <(Edge const& other) const
        {
            return w < other.w;
        }
    };

    std::vector<std::map<size_t, size_t>> edge_editions(graph.numberOfVertices());
    DynamicGraph original_graph_cp(graph.numberOfVertices());
    std::priority_queue<Edge> Q;

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
    {
        auto a = graph.vertexOfEdge(i, 0);
        auto b = graph.vertexOfEdge(i, 1);

        original_graph_cp.updateEdgeWeight(a, b, edge_values[i]);

        auto e = Edge(a, b, edge_values[i]);
        e.edition = ++edge_editions[e.a][e.b];

        Q.push(e);
    }

    andres::Partition<size_t> partition(graph.numberOfVertices());

    while (!Q.empty())
    {
        auto edge = Q.top();
        Q.pop();

        if (!original_graph_cp.edgeExists(edge.a, edge.b) || edge.edition < edge_editions[edge.a][edge.b])
            continue;

        if (edge.w < typename EVA::value_type())
            break;

        auto stable_vertex = edge.a;
        auto merge_vertex = edge.b;

        if (original_graph_cp.getAdjacentVertices(stable_vertex).size() < original_graph_cp.getAdjacentVertices(merge_vertex).size())
            std::swap(stable_vertex, merge_vertex);

        partition.merge(stable_vertex, merge_vertex);

        for (auto& p : original_graph_cp.getAdjacentVertices(merge_vertex))
        {
            if (p.first == stable_vertex)
                continue;

            original_graph_cp.updateEdgeWeight(stable_vertex, p.first, p.second);

            auto e = Edge(stable_vertex, p.first, original_graph_cp.getEdgeWeight(stable_vertex, p.first));
            e.edition = ++edge_editions[e.a][e.b];

            Q.push(e);
        }

        original_graph_cp.removeVertex(merge_vertex);
    }

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
        edge_labels[i] = partition.find(graph.vertexOfEdge(i, 0)) == partition.find(graph.vertexOfEdge(i, 1)) ? 0 : 1;
}

/// Greedy agglomerative decomposition of a graph with constraints.
///
template<typename GRAPH, typename EVA, typename ECA, typename ELA>
void constrainedGreedyAdditiveEdgeContraction(
    const GRAPH& graph,
    EVA const& edge_values,
    ECA const& constraints,
    ELA& vertex_labels
)
{
    class DynamicGraph
    {
    public:
        DynamicGraph(size_t n) :
            vertices_(n)
        {}

        bool edgeExists(size_t a, size_t b) const
        {
            return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
        }

        std::map<size_t, typename EVA::value_type> const& getAdjacentVertices(size_t v) const
        {
            return vertices_[v];
        }

        typename EVA::value_type getEdgeWeight(size_t a, size_t b) const
        {
            return vertices_[a].at(b);
        }

        void removeVertex(size_t v)
        {
            for (auto& p : vertices_[v])
                vertices_[p.first].erase(v);

            vertices_[v].clear();
        }

        void updateEdgeWeight(size_t a, size_t b, typename EVA::value_type w)
        {
#pragma omp atomic update
            vertices_[a][b] += w;
#pragma omp atomic update
            vertices_[b][a] += w;
        }

    private:
        std::vector<std::map<size_t, typename EVA::value_type>> vertices_;
    };

    struct Edge
    {
        Edge(size_t _a, size_t _b, typename EVA::value_type _w)
        {
            if (_a > _b)
                std::swap(_a, _b);

            a = _a;
            b = _b;

            w = _w;
        }

        size_t a;
        size_t b;
        size_t edition;
        typename EVA::value_type w;

        bool operator <(Edge const& other) const
        {
            return w < other.w;
        }
    };

    std::vector<std::map<size_t, size_t>> edge_editions(graph.numberOfVertices());
    DynamicGraph original_graph_cp(graph.numberOfVertices());
    std::priority_queue<Edge> Q;

    #pragma omp parallel for
    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
    {
        auto a = graph.vertexOfEdge(i, 0);
        auto b = graph.vertexOfEdge(i, 1);

        original_graph_cp.updateEdgeWeight(a, b, edge_values[i]);

        auto e = Edge(a, b, edge_values[i]);
        e.edition = ++edge_editions[e.a][e.b];

#pragma omp critical(UpdatePriorityQueue)
        Q.push(e);
    }

    andres::Partition<size_t> partition(graph.numberOfVertices());
    auto constraint_cp = constraints;

    while (!Q.empty())
    {
        auto edge = Q.top();
        Q.pop();

        if (!original_graph_cp.edgeExists(edge.a, edge.b) || edge.edition < edge_editions[edge.a][edge.b])
            continue;

        if (edge.w < typename EVA::value_type())
            break;

        auto stable_vertex = edge.a;
        auto merge_vertex = edge.b;

        if (original_graph_cp.getAdjacentVertices(stable_vertex).size() < original_graph_cp.getAdjacentVertices(merge_vertex).size())
            std::swap(stable_vertex, merge_vertex);

        auto it_merge = std::find_if(constraint_cp.begin(), constraint_cp.end(),
            [&merge_vertex](const std::pair<long, long>& element) {
                return element.first == merge_vertex || element.second == merge_vertex;
            }
        );
        auto it = it_merge;
        bool constrained = false;
        auto stable_set = partition.find(stable_vertex);
        while (it != constraint_cp.end())
        {
            if (it->first == merge_vertex) {
                if (partition.find(it->second) == stable_set) {
                    constrained = true;
                }
            }
            if (it->second == merge_vertex) {
                if (partition.find(it->first) == stable_set) {
                    constrained = true;
                }
            }
            if (constrained) {
                break;
            }
            ++it;
        }
        if (constrained) {
            continue;
        }

        partition.merge(stable_vertex, merge_vertex);

        while (it_merge != constraint_cp.end()) {
            if (it_merge->first == merge_vertex) {
                it_merge->first = stable_vertex;
            }
            if (it_merge->second == merge_vertex) {
                it_merge->second = stable_vertex;
            }
            ++it_merge;
        }

        using MapType=std::map<size_t, typename EVA::value_type>;
        MapType const & vertexMap = original_graph_cp.getAdjacentVertices(merge_vertex);
        std::vector<size_t> keys;
        keys.reserve(vertexMap.size());
        for (auto & vertex : vertexMap)
        {
          keys.push_back(vertex.first);
        }

#pragma omp parallel for
        for (size_t idx = 0; idx < keys.size(); idx++ )
        {
            size_t key = keys[idx];
            typename EVA::value_type value = vertexMap.at(key);
            if (key == stable_vertex)
                continue;

            original_graph_cp.updateEdgeWeight(stable_vertex, key, value);

            auto e = Edge(stable_vertex,key, original_graph_cp.getEdgeWeight(stable_vertex, key));
            e.edition = ++edge_editions[e.a][e.b];

#pragma omp critical(UpdatePriorityQueue)
            Q.push(e);
        }

        original_graph_cp.removeVertex(merge_vertex);
    }

    vertex_labels.resize(partition.numberOfElements());
    partition.elementLabeling(vertex_labels.begin());
}
} // namespace multicut
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_GREEDY_ADDITIVE_HXX
