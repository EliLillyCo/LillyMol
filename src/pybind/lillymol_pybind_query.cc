#include <iostream>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "Molecule_Lib/substructure.h"

#include "Molecule_Lib/substructure.pb.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_query, q)
{
  py::class_<Substructure_Results>(q, "SubstructureResults")
    .def(py::init<>())
    .def("number_embeddings", &Substructure_Results::number_embeddings, "Number embeddings")
    .def("each_embedding_set_vector",
      [](const Substructure_Results& sresults, int matoms, int value)->std::vector<int> {
        return sresults.EachEmbeddingSetVector(matoms, value);
      },
      "vector of all matched atoms"
    )
    .def("max_query_atoms_matched_in_search", &Substructure_Results::max_query_atoms_matched_in_search,
        "max atoms matched during last search"
    )
    .def("__iter__",
      [](const Substructure_Results& sresults) {
        return py::make_iterator(sresults.embeddings().begin(), sresults.embeddings().end());
      },
      py::keep_alive<0, 1>()
    )
    .def("__getitem__",
      [](const Substructure_Results& sresults, int ndx) {
        return sresults.embedding(ndx);
      },
      py::return_value_policy::reference
    )
  ;

  py::class_<Substructure_Query>(q, "SubstructureQuery")
    .def(py::init<>())
    .def("build_from_smarts", static_cast<int (Substructure_Query::*)(const std::string&)>(&Substructure_Query::CreateFromSmarts),
        "build from smarts")
    .def("set_only_keep_matches_in_largest_fragment", &Substructure_Query::set_only_keep_matches_in_largest_fragment, "set_only_keep_matches_in_largest_fragment")
    .def("set_embeddings_do_not_overlap", &Substructure_Query::set_embeddings_do_not_overlap, "set_embeddings_do_not_overlap")
    .def("set_find_one_embedding_per_atom", &Substructure_Query::set_find_one_embedding_per_atom, "set_find_one_embedding_per_atom")
    .def("set_find_unique_embeddings_only", &Substructure_Query::set_find_unique_embeddings_only, "set_find_unique_embeddings_only")
    .def("set_max_matches_to_find", &Substructure_Query::set_max_matches_to_find, "set_max_matches_to_find")
    .def("set_perceive_symmetry_equivalent_matches", &Substructure_Query::set_perceive_symmetry_equivalent_matches, "set_perceive_symmetry_equivalent_matches")
    .def("set_min_atoms_to_match", &Substructure_Query::set_min_atoms_to_match, "set_min_atoms_to_match")
    .def("set_max_atoms_to_match", &Substructure_Query::set_max_atoms_to_match, "set_max_atoms_to_match")
    .def("max_query_atoms_matched_in_search", &Substructure_Query::max_query_atoms_matched_in_search, "max_query_atoms_matched_in_search")
    .def("substructure_search", static_cast<int (Substructure_Query::*)(Molecule*)>(&Substructure_Query::substructure_search), "substructure_search")
    .def("substructure_search",
      [](Substructure_Query& qry, Molecule& m, Substructure_Results& sresults) {
        return qry.substructure_search(m, sresults);
      },
      "Substructure search"
    )
    .def("substructure_search_matches",
      [](Substructure_Query& qry, Molecule& m)->std::vector<Set_of_Atoms>{
        std::vector<Set_of_Atoms> results;
        Substructure_Results query_results;
        if (! qry.substructure_search(&m, query_results)) {
          return results;
        }

        results.reserve(query_results.number_embeddings());
        for (const Set_of_Atoms* s : query_results.embeddings()) {
          results.push_back(Set_of_Atoms(*s));
        }
        return results;
      },
      "Substructure search with matches"
    )
    .def("construct_from_proto", static_cast<int (Substructure_Query::*)(const SubstructureSearch::SubstructureQuery& proto)>(&Substructure_Query::ConstructFromProto), "Construct from proto")
    .def("__repr__",
      [](const Substructure_Query &q) {
        IWString rc;
        rc << "<SubstructureQuery " << q.comment() << '>';
        return std::string(rc.data(), rc.length());
      })
      .def("read_proto",
        [](Substructure_Query& qry, const std::string& fname)->bool{
          IWString tmp(fname.data(), fname.size());
          return qry.ReadProto(tmp);
        },
        "read from textproto file"
      )
      // Maybe one day when we get move constructures for a query
      /*.def("QueryFromSmarts",
        [](const std::string& smarts){
          return QueryFromSmarts(smarts);
        },
        "Query created from smarts - or None"
      )*/
    ;

}
