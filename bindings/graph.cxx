#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"

using namespace boost::python;
using namespace andres::graph;

BOOST_PYTHON_MODULE(graph) {
    class_<std::vector<double> >("DoubleVec")
        .def(vector_indexing_suite<std::vector<double> >())
    ;
    class_<std::vector<long> >("LongVec")
        .def(vector_indexing_suite<std::vector<long> >())
    ;
    class_<Graph<>>("Graph")
        .def("insertVertices", &Graph<>::insertVertices)
        .def("insertEdge", &Graph<>::insertEdge)
    ;
    def("greedyAdditiveEdgeContraction",
        &multicut::greedyAdditiveEdgeContraction<
            Graph<>,
            std::vector<double>,
            std::vector<long> 
        >
    );
}
