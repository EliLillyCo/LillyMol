// Python bindings for TSubstructure class.
#include <iostream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"  

#include "tsubstructure.h"

namespace py = pybind11;
using pybind_substructure::TSubstructure;

PYBIND11_MODULE(lillymol_tsubstructure, s) {
	py::class_<pybind_substructure::TSubstructure>(s, "TSubstructure")
	.def(py::init<>())
	.def("read_queries", &TSubstructure::ReadQueries)
	.def("add_query_from_smarts", &TSubstructure::AddQueryFromSmarts)
        .def("substructure_search",static_cast<bool (pybind_substructure::TSubstructure::*)(const std::string&)>(&TSubstructure::SubstructureSearch))
	.def("substructure_search",static_cast<bool (pybind_substructure::TSubstructure::*)(Molecule&)>(&TSubstructure::SubstructureSearch))
	.def("num_matches",
          [](TSubstructure& ts, const std::string& smi)->std::vector<int> {
            return ts.NumofMatches(smi);
          },
          "Number of times each query matches"
        )
	.def("num_matches",
          [](TSubstructure& ts, Molecule& m)->std::vector<int> {
            return ts.NumofMatches(m);
          },
          "Number of times each query matches"
        )
	.def("label_matched_atoms",
          [](TSubstructure& ts, const std::string& smi)->std::string {
            return ts.LabelMatchedAtoms(smi);
          },
          "Number of times each query matches"
        )
	.def("label_matched_atoms",
          [](TSubstructure& ts, Molecule& m)->std::string {
            return ts.LabelMatchedAtoms(m);
          },
          "Number of times each query matches"
        )
	.def_readwrite("reduce_to_largest_fragment", &TSubstructure::reduce_to_largest_fragment)
	.def_readwrite("make_implicit_hydrogens_explicit", &TSubstructure::make_implicit_hydrogens_explicit)
	.def_readwrite("unique_embeddings_only", &TSubstructure::unique_embeddings_only)
	.def_readwrite("find_one_embedding_per_root_atom", &TSubstructure::find_one_embedding_per_root_atom)
	.def_readwrite("perceive_symmetry_equivalent_matches", &TSubstructure::perceive_symmetry_equivalent_matches)
	.def_readwrite("max_matches_to_find", &TSubstructure::max_matches_to_find)
	.def_readwrite("query", &TSubstructure::query) 
	.def_readwrite("must_match_all_queries", &TSubstructure::must_match_all_queries)
	.def_readwrite("isotope", &TSubstructure::isotope)
	;
}
