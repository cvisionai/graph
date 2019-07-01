#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

void testMulticut() {
    // a simple weighted graph in which an optimal multicut is non-trivial

    andres::graph::Graph<> graph;
    graph.insertVertices(6);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 3); // 1
    graph.insertEdge(1, 2); // 2
    graph.insertEdge(1, 4); // 3
    graph.insertEdge(2, 5); // 4
    graph.insertEdge(3, 4); // 5
    graph.insertEdge(4, 5); // 6

    std::vector<double> weights(7);
    weights[0] = 5;
    weights[1] = -20;
    weights[2] = 5;
    weights[3] = 5;
    weights[4] = -20;
    weights[5] = 5;
    weights[6] = 5;

    std::vector<char> edge_labels(graph.numberOfEdges());
    andres::graph::multicut::greedyAdditiveEdgeContraction(graph, weights, edge_labels);

    test(edge_labels[0] == 0);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 0);
    test(edge_labels[3] == 1);
    test(edge_labels[4] == 1);
    test(edge_labels[5] == 0);
    test(edge_labels[6] == 0);
}

void testContraintedGreedyEdgeContraction()
{
  andres::graph::Graph<> graph;
    graph.insertVertices(8);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 2); // 1
    graph.insertEdge(1, 3); // 2
    graph.insertEdge(1, 4); // 3
    graph.insertEdge(1, 5); // 4
    graph.insertEdge(2, 3); // 5
    graph.insertEdge(2, 4); // 6
    graph.insertEdge(2, 5); // 7

    std::vector<double> weights(8);
    weights[0] = 0.9;
    weights[1] = 0.7;
    weights[2] = 2.5;
    weights[3] = 7.5;
    weights[4] = 0.8;
    weights[5] = 0.8;
    weights[6] = 2.5;
    weights[7] = 3.0;

    std::vector<std::pair<long,long>> constraints(4);
    constraints[0] = std::pair<long,long>(1,2);
    constraints[1] = std::pair<long,long>(3,4);
    constraints[2] = std::pair<long,long>(4,5);
    constraints[3] = std::pair<long,long>(3,5);

    std::vector<long> arg;
    andres::graph::multicut::constrainedGreedyAdditiveEdgeContraction(graph, weights, constraints, arg);


    test(arg[1] != arg[2]);
    test(arg[3] != arg[4]);
    test(arg[4] != arg[5]);
    test(arg[3] != arg[5]);
    test(arg[1] == arg[4]);
    test(arg[0] == arg[1]);
    test(arg[2] == arg[5]);

}

int main()
{
    testMulticut();
    testContraintedGreedyEdgeContraction();
    return 0;
}
