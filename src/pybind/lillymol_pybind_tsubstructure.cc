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
        .def("add_queries_from_smarts", &TSubstructure::AddQueriesFromSmarts, "Add multiple queries from smarts")
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
          "Smiles with matched atoms assigned an isotope"
        )
	.def("label_matched_atoms",
          [](TSubstructure& ts, Molecule& m)->int {
            return ts.LabelMatchedAtoms(m);
          },
          "Matched atoms in `m` are labelled with isotope"
        )
#ifdef THIS_DOES_NOT_WORK
        // due to how things are passed through pybind
        .def("label_matched_atoms",
          [](TSubstructure& ts, std::vector<Molecule>& mols)->uint32_t {
            int rc = 0;
            for (Molecule& m : mols) {
              if (ts.LabelMatchedAtoms(m)) {
                std::cerr << "After label " << m.unique_smiles() << '\n';
                ++rc;
              }
            }
            return rc;
          },
          "Apply isotopic labels to matched atoms"
        )
#endif
        .def("substructure_search",
          [](TSubstructure& ts, std::vector<Molecule>& mols) ->std::vector<bool> {
            const uint32_t number_molecules = mols.size();
            std::vector<bool> results(number_molecules);
            for (uint32_t i = 0; i < number_molecules; ++i) {
              results[i] = ts.SubstructureSearch(mols[i]);
            }
            return results;
          },
          "For each molecule, the number of queries matching"
        )
        .def("substructure_search",
          [](TSubstructure& ts, const std::vector<std::string>& smiles)->std::vector<bool>{
            const uint32_t nmols = smiles.size();
            std::vector<bool> result(nmols);
            Molecule m;
            for (uint32_t i = 0; i < nmols; ++i) {
              if (! m.build_from_smiles(smiles[i])) {
                std::cerr << "Invalid smiles '" << smiles[i] << "' ignored\n";
                result[i] = false;
              } else if (ts.SubstructureSearch(m)) {
                result[i] = true;
              } else {
                result[i] = false;
              }
            }

            return result;
          },
          "Perform substructure search on list of smiles"
        )
        .def("num_matches",
          [](TSubstructure& ts, std::vector<Molecule>& mols)->std::vector<std::vector<uint32_t>> {
            return ts.NumberMatches(mols);
          },
          "For each molecule, return the per query number of matches"
        )
        .def("all_queries_match",
          [](TSubstructure& ts, Molecule& m)->bool {
            ts.must_match_all_queries = true;
            return ts.SubstructureSearch(m);
          },
          ""
        )


        .def("number_queries", &TSubstructure::number_queries, "Number of queries defined")
        // These behaviour modifiers should be set before any queries are read.
	.def("set_reduce_to_largest_fragment", &TSubstructure::set_reduce_to_largest_fragment)
	.def("set_make_implicit_hydrogens_explicit", &TSubstructure::set_make_implicit_hydrogens_explicit)
        .def("set_label_by_query_number", &TSubstructure::set_label_by_query_number, "Label matched atoms by query number")
	.def("set_unique_embeddings_only", &TSubstructure::set_unique_embeddings_only, "Only find unique embeddings")
	.def("set_find_one_embedding_per_root_atom", &TSubstructure::set_find_one_embedding_per_root_atom, "for each first atom matched, find one embedding")
	.def("set_perceive_symmetry_equivalent_matches", &TSubstructure::set_perceive_symmetry_equivalent_matches, "Find symmetry related matches")
	.def("set_max_matches_to_find", &TSubstructure::set_max_matches_to_find, "The max number of embeddings to find")
        .def("set_labeled_smiles_are_unique", &TSubstructure::set_labeled_smiles_are_unique, "When labelled smiles are returned should they be unique")
	.def_readwrite("query", &TSubstructure::query) 
	.def_readwrite("must_match_all_queries", &TSubstructure::must_match_all_queries)
	.def_readwrite("isotope", &TSubstructure::isotope)
        .def("query_matched", &TSubstructure::QueryMatched, "Vector of how many times each query matched")
        .def("query_names", &TSubstructure::query_names, "Names of each query")
	;
}
