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
    class_<std::vector<float> >("FloatVec")
        .def(vector_indexing_suite<std::vector<float> >())
    ;
    class_<std::vector<long> >("LongVec")
        .def(vector_indexing_suite<std::vector<long> >())
    ;
    class_<std::pair<long, long> >("LongPair")
        .def(init<long, long>())
        .def_readwrite("first", &std::pair<long, long>::first)
        .def_readwrite("second", &std::pair<long, long>::second)
    ;
    class_<std::vector<std::pair<long, long> > >("PairVec")
        .def(vector_indexing_suite<std::vector<std::pair<long, long> > >())
    ;
    class_<Graph<>>("Graph")
        .def("insertVertices", &Graph<>::insertVertices)
        .def("insertEdge", &Graph<>::insertEdge)
    ;
    def("greedyAdditiveEdgeContraction",
        &multicut::greedyAdditiveEdgeContraction<
            Graph<>,
            std::vector<float>,
            std::vector<long> 
        >
    );
    def("constrainedGreedyAdditiveEdgeContraction",
        &multicut::constrainedGreedyAdditiveEdgeContraction<
            Graph<>,
            std::vector<float>,
            std::vector<std::pair<long, long> >,
            std::vector<long> 
        >
    );
}
