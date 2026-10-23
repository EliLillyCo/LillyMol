#include "Molecule_Lib/target.h"
#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

using pybind_substructure::TSubstructure;

void
BindTSubstructure(nb::module_& m) {
  // How aromatic bonds in a molecule are matched by bonds in a query. This is
  // process wide state in the C++ library, not a property of a TSubstructure, and
  // it is the same setting as -M kekule, -M nokekule and -M fkekule in the
  // tsubstructure command line tool.
  //
  // The range is checked here because the C++ setter takes any int, and only
  // 0, 1 and 2 have defined meanings.
  m.def("set_aromatic_bonds_lose_kekule_identity",
        [](int mode) {
          if (mode < 0 || mode > 2) {
            throw nb::value_error(
                "set_aromatic_bonds_lose_kekule_identity:mode must be 0, 1 or 2");
          }
          set_aromatic_bonds_lose_kekule_identity(mode);
        },
        nb::arg("mode"),
        "Set how aromatic bonds match the bonds in substructure queries. "
        "0 (the default): an aromatic bond matches an aromatic query bond, and also "
        "its Kekule single or double form, so 'a-a' matches the bonds inside a "
        "benzene ring. "
        "1: an aromatic bond matches only aromatic query bonds, so 'a-a' means two "
        "aromatic atoms joined by a non aromatic single bond. "
        "2: as 1, but aromatic rings that have no alternative Kekule form keep their "
        "Kekule bond types. "
        "This is process wide, not per TSubstructure. It is applied when a search is "
        "done, and when queries are built from molecules, so set it before either. "
        "Creating an iwdescr descriptor calculator resets it to 0.");
  m.def("aromatic_bonds_lose_kekule_identity",
        []() { return aromatic_bonds_lose_kekule_identity(); },
        "The current value of the aromatic bond matching mode. See "
        "set_aromatic_bonds_lose_kekule_identity.");

  nb::class_<TSubstructure>(m, "TSubstructure")
      .def(nb::init<>())
      .def("read_queries", &TSubstructure::ReadQueries, nb::arg("directive"))
      .def("add_query_from_smarts", &TSubstructure::AddQueryFromSmarts,
           nb::arg("smarts"))
      .def("add_queries_from_smarts", &TSubstructure::AddQueriesFromSmarts,
           nb::arg("smarts"))
      .def("substructure_search",
           [](TSubstructure& ts, const std::string& smiles) {
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(smiles);
           },
           nb::arg("smiles"))
      .def("substructure_search",
           [](TSubstructure& ts, Molecule& mol) {
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(mol);
           },
           nb::arg("mol"))
      // Note the Molecule*. With std::vector<Molecule> nanobind copy constructs
      // every molecule at the boundary, about 5.4us each, and does it while
      // holding the GIL so it cannot overlap between threads. Measured on the
      // pybind equivalent, that made this batch call slower than a python loop
      // over the single molecule overload, and capped scaling at under 3x on 32
      // cores; with pointers it reached 10x. The python signature is the same.
      //
      // Taking pointers means None arrives as nullptr, so the list is checked
      // before the GIL is released - throwing needs the GIL, and dereferencing
      // nullptr takes the interpreter down.
      .def("substructure_search",
           [](TSubstructure& ts, std::vector<Molecule*> mols) {
             const uint32_t number_molecules = mols.size();
             for (uint32_t i = 0; i < number_molecules; ++i) {
               if (mols[i] == nullptr) {
                 throw nb::value_error("substructure_search:None in list of molecules");
               }
             }

             std::vector<bool> results(number_molecules);
             {
               nb::gil_scoped_release release;
               for (uint32_t i = 0; i < number_molecules; ++i) {
                 results[i] = ts.SubstructureSearch(*mols[i]);
               }
             }
             return results;
           },
           nb::arg("mols"))
      .def("substructure_search",
           [](TSubstructure& ts, const std::vector<std::string>& smiles) {
             const uint32_t number_molecules = smiles.size();
             std::vector<bool> results(number_molecules);
             {
               nb::gil_scoped_release release;
               Molecule mol;
               for (uint32_t i = 0; i < number_molecules; ++i) {
                 if (!mol.build_from_smiles(smiles[i])) {
                   std::cerr << "Invalid smiles '" << smiles[i] << "' ignored\n";
                   results[i] = false;
                 } else {
                   results[i] = ts.SubstructureSearch(mol);
                 }
               }
             }
             return results;
           },
           nb::arg("smiles"))
      .def("num_matches",
           [](TSubstructure& ts, const std::string& smiles) {
             nb::gil_scoped_release release;
             return ts.NumofMatches(smiles);
           },
           nb::arg("smiles"))
      .def("num_matches",
           [](TSubstructure& ts, Molecule& mol) {
             nb::gil_scoped_release release;
             return ts.NumofMatches(mol);
           },
           nb::arg("mol"))
      // Molecule*, for the reason given above the batch substructure_search.
      .def("num_matches",
           [](TSubstructure& ts, std::vector<Molecule*> mols) {
             for (uint32_t i = 0; i < mols.size(); ++i) {
               if (mols[i] == nullptr) {
                 throw nb::value_error("num_matches:None in list of molecules");
               }
             }
             nb::gil_scoped_release release;
             return ts.NumberMatches(mols);
           },
           nb::arg("mols"))
      .def("number_matches",
           [](TSubstructure& ts, const std::vector<std::string>& smiles) {
             nb::gil_scoped_release release;
             return ts.NumberMatches(smiles);
           },
           nb::arg("smiles"))
      .def("label_matched_atoms",
           [](TSubstructure& ts, const std::string& smiles) {
             return ts.LabelMatchedAtoms(smiles);
           },
           nb::arg("smiles"))
      .def("label_matched_atoms",
           [](TSubstructure& ts, Molecule& mol) {
             return ts.LabelMatchedAtoms(mol);
           },
           nb::arg("mol"))
      // Only possible taking pointers. With std::vector<Molecule> the labels are
      // applied to the copies nanobind made at the boundary and then discarded,
      // which is why the pybind bindings carried this under an ifdef for years.
      .def("label_matched_atoms",
           [](TSubstructure& ts, std::vector<Molecule*> mols) {
             for (uint32_t i = 0; i < mols.size(); ++i) {
               if (mols[i] == nullptr) {
                 throw nb::value_error("label_matched_atoms:None in list of molecules");
               }
             }

             uint32_t rc = 0;
             nb::gil_scoped_release release;
             for (Molecule* m : mols) {
               if (ts.LabelMatchedAtoms(*m)) {
                 ++rc;
               }
             }
             return rc;
           },
           nb::arg("mols"))
      .def("matched_atoms",
           [](TSubstructure& ts, const std::string& smiles) {
             return ts.MatchedAtoms(smiles);
           },
           nb::arg("smiles"))
      .def("matched_atoms",
           [](TSubstructure& ts, Molecule& mol) {
             return ts.MatchedAtoms(mol);
           },
           nb::arg("mol"))
      .def("all_queries_match",
           [](TSubstructure& ts, Molecule& mol) {
             ts.must_match_all_queries = true;
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(mol);
           },
           nb::arg("mol"))
      .def("number_queries", &TSubstructure::number_queries)
      .def("set_reduce_to_largest_fragment", &TSubstructure::set_reduce_to_largest_fragment,
           nb::arg("value"))
      .def("set_make_implicit_hydrogens_explicit",
           &TSubstructure::set_make_implicit_hydrogens_explicit, nb::arg("value"))
      .def("set_label_by_query_number", &TSubstructure::set_label_by_query_number,
           nb::arg("value"))
      .def("set_unique_embeddings_only", &TSubstructure::set_unique_embeddings_only,
           nb::arg("value"))
      .def("set_find_one_embedding_per_root_atom",
           &TSubstructure::set_find_one_embedding_per_root_atom, nb::arg("value"))
      .def("set_perceive_symmetry_equivalent_matches",
           &TSubstructure::set_perceive_symmetry_equivalent_matches, nb::arg("value"))
      .def("set_max_matches_to_find", &TSubstructure::set_max_matches_to_find,
           nb::arg("value"))
      .def("set_labeled_smiles_are_unique", &TSubstructure::set_labeled_smiles_are_unique,
           nb::arg("value"))
      .def_rw("query", &TSubstructure::query)
      .def_rw("must_match_all_queries", &TSubstructure::must_match_all_queries)
      .def_rw("isotope", &TSubstructure::isotope)
      .def("query_matched", &TSubstructure::QueryMatched)
      .def("query_names", &TSubstructure::query_names);
}

}  // namespace lillymol_nb
