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
        // These release the GIL for the duration of the search, so that a python
        // program calling them from several threads gets real concurrency. See
        // the note above the batch entry points for the safety contract: a
        // TSubstructure is thread compatible, not thread safe, and a search
        // writes to the molecule it searches.
        .def("substructure_search",
          [](TSubstructure& ts, const std::string& smi)->bool {
            py::gil_scoped_release release;
            return ts.SubstructureSearch(smi);
          },
          "Substructure search a smiles"
        )
	.def("substructure_search",
          [](TSubstructure& ts, Molecule& m)->bool {
            py::gil_scoped_release release;
            return ts.SubstructureSearch(m);
          },
          "Substructure search a molecule"
        )
	.def("num_matches",
          [](TSubstructure& ts, const std::string& smi)->std::vector<int> {
            py::gil_scoped_release release;
            return ts.NumofMatches(smi);
          },
          "Number of times each query matches"
        )
	.def("num_matches",
          [](TSubstructure& ts, Molecule& m)->std::vector<int> {
            py::gil_scoped_release release;
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
        // This could not work with std::vector<Molecule>, because pybind copies
        // the list at the boundary and the labels were applied to the copies and
        // then discarded. std::vector<Molecule*> labels the caller's molecules,
        // which is the whole point of the method.
        .def("label_matched_atoms",
          [](TSubstructure& ts, std::vector<Molecule*> mols)->uint32_t {
            for (uint32_t i = 0; i < mols.size(); ++i) {
              if (mols[i] == nullptr) {
                throw py::value_error("label_matched_atoms:None at index " + std::to_string(i));
              }
            }

            uint32_t rc = 0;
            py::gil_scoped_release release;
            for (Molecule* m : mols) {
              if (ts.LabelMatchedAtoms(*m)) {
                ++rc;
              }
            }
            return rc;
          },
          "Apply isotopic labels to matched atoms, in place. Returns how many molecules were labelled"
        )
        // The batch entry points below release the GIL for the whole loop, so
        // other python threads run while the search proceeds. This matters more
        // here than for a single search: holding the GIL across a batch of a few
        // thousand molecules stalls every other python thread for tens of
        // milliseconds. Nothing in the released region touches a python object.
        //
        // A TSubstructure is thread compatible, not thread safe - it holds the
        // queries, which accumulate match state - so give each thread its own.
        // A search also writes to the molecule it searches, forcing ring and
        // aromaticity perception, so molecules must not be shared between
        // threads either. For the Molecule overload that is automatic, since
        // pybind copies the list into a std::vector<Molecule> (the same reason
        // the label_matched_atoms batch overload above cannot work); for the
        // smiles overload the molecules are built inside C++.
        //
        // Note the Molecule* below. It is not decoration.
        //
        // With std::vector<Molecule> pybind11 has to convert the python list into
        // a vector of Molecule before the search starts, copy constructing every
        // one, about 5.4us each. list_caster reserves, so no reallocation is
        // involved - it is the element copies. Worse, that conversion happens
        // while the GIL is held, so it cannot overlap between threads, and it
        // capped the scaling at under 3x. Measured that way this batch call was
        // SLOWER than simply calling once per molecule, at every thread count,
        // which made the obvious looking optimisation a pessimisation.
        //
        // std::vector<Molecule*> takes pointers to the molecules the caller
        // already has, so nothing is copied. Same python signature - callers
        // still pass a list of Molecule. Measured on 32 cores, one query, 4000
        // drug sized molecules per thread:
        //
        //                                1 thread   4 threads     16 threads
        //   one call per molecule           54k      181k (3.4x)   199k (3.7x)
        //   batch, vector<Molecule>         36k       83k (2.3x)    99k (2.7x)
        //   batch, vector<Molecule*>        57k      218k (3.8x)   570k (10.0x)
        //   batch, vector<string>           45k      163k (3.7x)   517k (11.6x)
        //
        // So keep the star. If a batch entry point here is ever slower than the
        // loop it replaces, look at what pybind is copying at the boundary before
        // looking anywhere else.
        //
        // One consequence of not copying. Preprocess() may reduce a molecule to
        // its largest fragment or make its implicit Hydrogens explicit, if those
        // options are set, and those now happen to the caller's molecules rather
        // than to throwaway copies. That matches what the single molecule
        // overload has always done - the copy was the anomaly, not the mutation -
        // and with the default options Preprocess does nothing.
        .def("substructure_search",
          [](TSubstructure& ts, std::vector<Molecule*> mols) ->std::vector<bool> {
            const uint32_t number_molecules = mols.size();
            // Taking pointers means pybind converts a python None to nullptr,
            // where the by value form would have raised. Check before releasing
            // the GIL, since throwing needs it. Without this, a None in the list
            // dereferences nullptr and takes the interpreter down.
            for (uint32_t i = 0; i < number_molecules; ++i) {
              if (mols[i] == nullptr) {
                throw py::value_error("substructure_search:None at index " + std::to_string(i));
              }
            }

            std::vector<bool> results(number_molecules);
            {
              py::gil_scoped_release release;
              for (uint32_t i = 0; i < number_molecules; ++i) {
                results[i] = ts.SubstructureSearch(*mols[i]);
              }
            }
            return results;
          },
          "For each molecule, whether the queries matched"
        )
        .def("substructure_search",
          [](TSubstructure& ts, const std::vector<std::string>& smiles)->std::vector<bool>{
            const uint32_t nmols = smiles.size();
            std::vector<bool> result(nmols);
            {
              py::gil_scoped_release release;
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
            }

            return result;
          },
          "Perform substructure search on list of smiles"
        )
        // Molecule*, for the reason given above the batch substructure_search.
        .def("num_matches",
          [](TSubstructure& ts, std::vector<Molecule*> mols)->std::vector<std::vector<uint32_t>> {
            for (uint32_t i = 0; i < mols.size(); ++i) {
              if (mols[i] == nullptr) {
                throw py::value_error("num_matches:None at index " + std::to_string(i));
              }
            }
            py::gil_scoped_release release;
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
	.def("set_reduce_to_largest_fragment", &TSubstructure::set_reduce_to_largest_fragment, "Transform molecule before matching")
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
