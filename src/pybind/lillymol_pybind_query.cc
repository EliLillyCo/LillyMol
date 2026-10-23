#include <iostream>
#include <memory>
#include <optional>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Lib/substructure.pb.h"

namespace py = pybind11;

namespace {

struct SmartsSearchOptions {
  std::optional<int> max_matches_to_find;
  std::optional<bool> unique_embeddings_only;
  std::optional<bool> one_embedding_per_start_atom;
  std::optional<bool> perceive_symmetry_equivalent_matches;
};

std::unique_ptr<Substructure_Query>
BuildQueryFromSmarts(const std::string& smarts,
                     const SmartsSearchOptions& options) {
  auto query = std::make_unique<Substructure_Query>();
  if (! query->CreateFromSmarts(smarts)) {
    throw py::value_error("Invalid smarts");
  }

  if (options.max_matches_to_find) {
    query->set_max_matches_to_find(*options.max_matches_to_find);
  }
  if (options.unique_embeddings_only) {
    query->set_find_unique_embeddings_only(*options.unique_embeddings_only);
  }
  if (options.one_embedding_per_start_atom) {
    query->set_find_one_embedding_per_atom(*options.one_embedding_per_start_atom);
  }
  if (options.perceive_symmetry_equivalent_matches) {
    query->set_perceive_symmetry_equivalent_matches(
        *options.perceive_symmetry_equivalent_matches);
  }

  return query;
}

std::vector<std::vector<int>>
SubstructureSearchMatches(Substructure_Query& query, Molecule& molecule) {
  Substructure_Results query_results;
  if (! query.substructure_search(&molecule, query_results)) {
    return std::vector<std::vector<int>>();
  }

  std::vector<std::vector<int>> results;
  results.reserve(query_results.number_embeddings());
  for (const Set_of_Atoms* embedding : query_results.embeddings()) {
    std::vector<int> atoms;
    atoms.reserve(embedding->number_elements());
    for (atom_number_t atom : *embedding) {
      atoms.push_back(atom);
    }
    results.push_back(std::move(atoms));
  }

  return results;
}

}  // namespace

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
    .def("build_from_molecule", [](Substructure_Query& q, Molecule& m)->bool {
        Molecule_to_Query_Specifications mqs;
        mqs.set_make_embedding(1);
        return q.create_from_molecule(m, mqs);
      },
      "Convert `m` to a query with default conditions"
    )
    .def("set_only_keep_matches_in_largest_fragment", &Substructure_Query::set_only_keep_matches_in_largest_fragment, "set_only_keep_matches_in_largest_fragment")
    .def("set_embeddings_do_not_overlap", &Substructure_Query::set_embeddings_do_not_overlap, "set_embeddings_do_not_overlap")
    .def("set_find_one_embedding_per_atom", &Substructure_Query::set_find_one_embedding_per_atom, "set_find_one_embedding_per_atom")
    .def("set_find_unique_embeddings_only", &Substructure_Query::set_find_unique_embeddings_only, "set_find_unique_embeddings_only")
    .def("set_max_matches_to_find", &Substructure_Query::set_max_matches_to_find, "set_max_matches_to_find")
    .def("set_perceive_symmetry_equivalent_matches", &Substructure_Query::set_perceive_symmetry_equivalent_matches, "set_perceive_symmetry_equivalent_matches")
    .def("set_min_atoms_to_match", &Substructure_Query::set_min_atoms_to_match, "set_min_atoms_to_match")
    .def("set_max_atoms_to_match", &Substructure_Query::set_max_atoms_to_match, "set_max_atoms_to_match")
    .def("max_query_atoms_matched_in_search", &Substructure_Query::max_query_atoms_matched_in_search, "max_query_atoms_matched_in_search")
    // The searches below release the GIL for the duration of the C++ work, so
    // that other python threads can run. A search is 6 to 18 microseconds, which
    // is long enough for that to be worth the ~1 microsecond a release and
    // reacquire pair costs. Nothing inside the released region touches a python
    // object; the arguments were converted before it and the return value is
    // converted after it.
    //
    // This does NOT make searching thread safe, and cannot. A Substructure_Query
    // is thread compatible, not thread safe - it accumulates match state such as
    // _max_query_atoms_matched_in_search - so give each thread its own query.
    // Building one from smarts costs about 7 microseconds, so that is cheap.
    // Less obviously, a search also WRITES to the molecule it searches, forcing
    // ring and aromaticity perception on it, so two threads must not search the
    // same Molecule either.
    .def("substructure_search",
      [](Substructure_Query& qry, Molecule& m)->uint32_t {
        py::gil_scoped_release release;
        return qry.substructure_search(&m);
      },
      "substructure_search"
    )
    .def("substructure_search",
      [](Substructure_Query& qry, Molecule& m, Substructure_Results& sresults)->uint32_t {
        py::gil_scoped_release release;
        return qry.substructure_search(m, sresults);
      },
      "Substructure search"
    )
    .def("substructure_search_matches",
      [](Substructure_Query& qry, Molecule& m)->std::optional<std::vector<Set_of_Atoms>>{
        std::vector<Set_of_Atoms> results;
        {
          py::gil_scoped_release release;
          Substructure_Results query_results;
          if (! qry.substructure_search(&m, query_results)) {
            return std::nullopt;
          }

          results.reserve(query_results.number_embeddings());
          for (const Set_of_Atoms* s : query_results.embeddings()) {
            results.push_back(Set_of_Atoms(*s));
          }
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
    .def("read_msi",
      [](Substructure_Query& qry, const std::string& fname)->bool {
        const IWString myfname(fname);
        return qry.read(myfname);
      },
      "read MSI style query from 'fname'"
    )
  ;

  q.def("QueryFromSmarts",
    [](const std::string& smarts)->std::unique_ptr<Substructure_Query> {
      auto result = std::make_unique<Substructure_Query>();
      if (! result->CreateFromSmarts(smarts)) {
        return nullptr;
      }
      return result;
    },
    py::arg("smarts"),
    "Return a SubstructureQuery built from `smarts`, or None if parsing fails."
  );

  q.def("HasSubstructMatch",
    [](Molecule& molecule, const std::string& smarts,
       std::optional<int> max_matches_to_find,
       std::optional<bool> unique_embeddings_only,
       std::optional<bool> one_embedding_per_start_atom,
       std::optional<bool> perceive_symmetry_equivalent_matches)->bool {
      SmartsSearchOptions options;
      options.max_matches_to_find = max_matches_to_find;
      options.unique_embeddings_only = unique_embeddings_only;
      options.one_embedding_per_start_atom = one_embedding_per_start_atom;
      options.perceive_symmetry_equivalent_matches = perceive_symmetry_equivalent_matches;

      std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(smarts, options);
      py::gil_scoped_release release;
      return query->substructure_search(&molecule);
    },
    py::arg("molecule"),
    py::arg("smarts"),
    py::kw_only(),
    py::arg("max_matches_to_find") = py::none(),
    py::arg("unique_embeddings_only") = py::none(),
    py::arg("one_embedding_per_start_atom") = py::none(),
    py::arg("perceive_symmetry_equivalent_matches") = py::none(),
    "Return True if `smarts` matches `molecule`."
  );

  q.def("CountSubstructMatches",
    [](Molecule& molecule, const std::string& smarts,
       std::optional<int> max_matches_to_find,
       std::optional<bool> unique_embeddings_only,
       std::optional<bool> one_embedding_per_start_atom,
       std::optional<bool> perceive_symmetry_equivalent_matches)->uint32_t {
      SmartsSearchOptions options;
      options.max_matches_to_find = max_matches_to_find;
      options.unique_embeddings_only = unique_embeddings_only;
      options.one_embedding_per_start_atom = one_embedding_per_start_atom;
      options.perceive_symmetry_equivalent_matches = perceive_symmetry_equivalent_matches;

      std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(smarts, options);
      py::gil_scoped_release release;
      return query->substructure_search(&molecule);
    },
    py::arg("molecule"),
    py::arg("smarts"),
    py::kw_only(),
    py::arg("max_matches_to_find") = py::none(),
    py::arg("unique_embeddings_only") = py::none(),
    py::arg("one_embedding_per_start_atom") = py::none(),
    py::arg("perceive_symmetry_equivalent_matches") = py::none(),
    "Return the number of embeddings for `smarts` in `molecule`."
  );

  q.def("GetSubstructMatches",
    [](Molecule& molecule, const std::string& smarts,
       std::optional<int> max_matches_to_find,
       std::optional<bool> unique_embeddings_only,
       std::optional<bool> one_embedding_per_start_atom,
       std::optional<bool> perceive_symmetry_equivalent_matches)->std::vector<std::vector<int>> {
      SmartsSearchOptions options;
      options.max_matches_to_find = max_matches_to_find;
      options.unique_embeddings_only = unique_embeddings_only;
      options.one_embedding_per_start_atom = one_embedding_per_start_atom;
      options.perceive_symmetry_equivalent_matches = perceive_symmetry_equivalent_matches;

      std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(smarts, options);
      py::gil_scoped_release release;
      return SubstructureSearchMatches(*query, molecule);
    },
    py::arg("molecule"),
    py::arg("smarts"),
    py::kw_only(),
    py::arg("max_matches_to_find") = py::none(),
    py::arg("unique_embeddings_only") = py::none(),
    py::arg("one_embedding_per_start_atom") = py::none(),
    py::arg("perceive_symmetry_equivalent_matches") = py::none(),
    "Return embeddings as lists of atom numbers, in query atom order."
  );

}
