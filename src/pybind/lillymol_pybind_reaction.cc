#include <iostream>
#include <memory>
#include <optional>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"
// #include "pybind11/operators.h"

#include "Foundational/iwmisc/proto_support.h"
#include "Molecule_Lib/iwreaction.h"

namespace py = pybind11;

// The reaction entry points below take std::vector<Molecule*> rather than
// std::vector<Molecule>. With the by value form pybind copy constructs every
// molecule in the list at the boundary, about 5.4us each, which is a third of the
// cost of a whole reaction - measured 10.50us passing one sidechain by reference
// against 13.92us passing it in a list. Over an enumeration driven by a million
// python calls that is seconds of pure copying. Taking pointers copies nothing and
// does not change the python signature; callers pass the same list.
//
// Nothing is exposed by not copying: add_sidechain_reagent builds a
// Molecule_and_Embedding from a copy of what it is given, so the caller's
// molecules are not touched. Verified - a reaction leaves its sidechain molecule
// byte identical.
//
// What pointers do expose is None, which pybind converts to nullptr, and every
// use below dereferences. One check, so the eight call sites cannot disagree.
static void
CheckNoNullMolecules(const std::vector<Molecule*>& mols, const char* caller) {
  for (uint32_t i = 0; i < mols.size(); ++i) {
    if (mols[i] == nullptr) {
      throw py::value_error(std::string(caller) + ":None at index " + std::to_string(i));
    }
  }
}

static bool
ReadReaction(const std::string& fname, IWReaction& rxn) {
  IWString tmp(fname.data(), fname.size());
  std::optional<ReactionProto::Reaction> proto =
                iwmisc::ReadTextProto<ReactionProto::Reaction>(tmp);
  if (! proto) {
    return false;
  }

  Sidechain_Match_Conditions smc;

  if (! rxn.ConstructFromProto(*proto, tmp, smc)) {
    return false;
  }

  return true;
}

PYBIND11_MODULE(lillymol_reaction, rxn) 
{
  py::class_<Sidechain_Match_Conditions>(rxn, "SidechainMatchConditions")
    .def(py::init<>())
    .def("set_make_new_reagent_for_each_hit", &Sidechain_Match_Conditions::set_make_new_reagent_for_each_hit, "Each match generates regioisomer")
    .def("set_max_matches_to_find", &Sidechain_Match_Conditions::set_max_matches_to_find, "max matches")
    .def("set_strip_reagents_to_largest_fragment", &Sidechain_Match_Conditions::set_strip_reagents_to_largest_fragment, "use largest fragment")
    .def("set_ignore_not_reacting", &Sidechain_Match_Conditions::set_ignore_not_reacting, "ignore non matches")
    .def("set_find_unique_embeddings_only", &Sidechain_Match_Conditions::set_find_unique_embeddings_only, "unique embeddings")
    .def("set_one_embedding_per_start_atom", &Sidechain_Match_Conditions::set_one_embedding_per_start_atom, "One embedding per start atom")
    .def("set_ignore_symmetry_related_matches", &Sidechain_Match_Conditions::set_ignore_symmetry_related_matches, "Ignore symmetry")
  ;

  py::class_<Reaction_Iterator>(rxn, "ReactionIterator")
    .def(py::init<>())
    .def(py::init<const IWReaction&>())
    .def("initialise", &Reaction_Iterator::initialise, "Initialise for reaction")
    .def("active",
      [](const Reaction_Iterator& rxnit)->bool{
        return rxnit.active();
      },
      "True if still active"
    )
    .def("increment",
      [](Reaction_Iterator& rxnit) {
        rxnit++;
      },
      "Move to next reagent"
    )
    .def("reagent", &Reaction_Iterator::reagent, "Get reagent for sidechain")
    .def("reset", &Reaction_Iterator::reset, "Reset the iterator")
    .def("debug_print", [](const Reaction_Iterator& rxnit) {
      rxnit.debug_print(std::cerr);
      },
      "debugging info"
    )
  ;

  py::class_<IWReaction, Substructure_Query>(rxn, "Reaction")
    .def(py::init<>())
    .def("name",
      [](const IWReaction& rxn)->std::string{
        return rxn.comment().AsString();
      },
      "Name"
    )
    .def("set_name",
      [](IWReaction& rxn, const std::string& name) {
        rxn.set_comment(name);
      },
      "Assign new name to reaction"
    )
    .def("read",
      [](IWReaction& rxn, const std::string& fname)->bool{
        return ReadReaction(fname, rxn);
      },
      "read textproto reaction"
    )
    // IWReaction does not have a move operator, maybe that is why this does not work.
    //.def("rxn_from_file",
    //  [](const std::string& fname)->std::optional<IWReaction>{
    //    IWReaction result;
    //    if (ReadReaction(fname, result)) {
    //      return result;
    //    }
    //    return std::nullopt;
    //  },
    //  "Read textproto reaction"
    //)

    .def("construct_from_smirks",
      [](IWReaction& rxn, const std::string& smirks)->bool{
        const const_IWSubstring tmp(smirks.data(), smirks.size());
        return rxn.construct_from_smirks(tmp);
      },
      "From smirks"
    )
    .def("construct_from_textproto",
      [](IWReaction& rxn, const std::string& textproto)->bool {
        IWString dirname(".");  // maybe make an argument
        Sidechain_Match_Conditions smc;
        if (! rxn.ConstructFromTextProto(textproto, dirname, smc))  {
          return false;
        }

        return true;
      },
      ""
    )
    .def("write_msi",
      [](IWReaction& rxn, const std::string& fname) {
        IWString tmp(fname);
        return rxn.write_msi(tmp);
      }
    )
    .def("sidechain_name",
      [](const IWReaction& rxn, int sidechain_number, int reagent_number)->std::optional<std::string> {
        const Sidechain_Reaction_Site * sidechain = rxn.sidechain(sidechain_number);
        if (sidechain == nullptr) {
          return std::nullopt;
        }
        const Molecule_and_Embedding* reagent = sidechain->reagent(reagent_number);
        if (reagent == nullptr) {
          return std::nullopt;
        }
        const IWString& s = reagent->name();

        std::string rc(s.data(), s.length());
        return rc;
      },
      "sidechain_name(sidechain_number, reagent_number)   The name of reagent `reagent_number` in sidechain `sidechain_number`"
    )
    .def("number_sidechains", &IWReaction::number_sidechains, "Number of sidechains")
    .def("number_sidechains_with_reagents", &IWReaction::number_sidechains_with_reagents, "number_sidechains_with_reagents")
    .def("set_one_embedding_per_start_atom", &IWReaction::set_one_embedding_per_start_atom, "one embedding per start atom")
    .def("add_sidechain_reagents",
      [](IWReaction& rxn, int sidechain, const char* fname, FileType file_type, Sidechain_Match_Conditions& smc)->bool{
        return rxn.add_sidechain_reagents(sidechain, fname, file_type, smc);
      },
      "Add reagents to a sidechain"
    )
    .def("add_sidechain_reagent",
      [](IWReaction& rxn, int sidechain, Molecule& m, const Sidechain_Match_Conditions& smc)->bool {
        return rxn.add_sidechain_reagent(sidechain, m, smc);
      },
      ""
    )
    // Deliberately not binding remove_no_delete_all_reagents. Reagents added from
    // python go through add_sidechain_reagent, which copies, so the reaction owns
    // them and not deleting them leaks one Molecule_and_Embedding per reagent.
    .def("remove_all_reagents", &IWReaction::remove_all_reagents,
         "Destroy all reagents in all sidechains. Returns the number of sidechains"
    )
    .def("reagent_names",
      [](const IWReaction& rxn, const Reaction_Iterator& iter)->std::vector<std::string> {
        const int n = iter.number_sidechains();
        std::vector<std::string> result(n);
        for (int i = 0; i < n; ++i) {
          int j = iter.reagent(i);
          const Sidechain_Reaction_Site* r = rxn.sidechain(i);
          const Molecule_and_Embedding* s = r->reagent(j);
          std::string tmp(s->name().data(), s->name().length());
          result.emplace_back(std::move(tmp));
        }
        return result;
      },
      "return the sidechains at the current position"
    )
    .def("substructure_search",
      [](IWReaction& rxn, Molecule& m, Substructure_Results& sresults) {
        return rxn.substructure_search(m, sresults);
      }
    )
    .def("substructure_search_matches",
      [](IWReaction& rxn, Molecule& m)->std::optional<std::vector<Set_of_Atoms>>{
        Substructure_Results sresults;
        if (! rxn.substructure_search(m, sresults)) {
          return std::nullopt;
        }

        std::vector<Set_of_Atoms> rc;
        rc.reserve(sresults.number_embeddings());

        for (const Set_of_Atoms* s : sresults.embeddings()) {
          rc.push_back(Set_of_Atoms(*s));
        }

        return rc;
      },
      "perform substructure search against scaffold queries"
    )
    .def("substructure_search",
      [](IWReaction& rxn, Molecule& m, Substructure_Results& sresults)->bool {
        if (! rxn.substructure_search(m, sresults)) {
          return false;
        }

        return true;
      },
      "perform substructure_search against scaffold queries"
    )
    .def("in_place_transformations",
      [](IWReaction& rxn, Molecule& m)->bool{
        return rxn.in_place_transformations(m);
      },
      "apply reaction to 'm'"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, const Set_of_Atoms* embedding, Molecule& product)->bool{
        return rxn.perform_reaction(&scaffold, embedding, product);
      },
      "perform reaction"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, const Set_of_Atoms* embedding)->std::optional<Molecule>{
        Molecule product;
        if (! rxn.perform_reaction(&scaffold, embedding, product)) {
          return std::nullopt;
        }
        return product;
      },
      "Perform reaction with particular set of matched atoms"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, const Molecule& scaffold, const Set_of_Atoms* embedding,
         const Reaction_Iterator& iter)->std::optional<Molecule>{
        Molecule product;
        if (! rxn.perform_reaction(&scaffold, embedding, iter, product)) {
          return std::nullopt;
        }
        return product;
      },
      "generate product based on embedding and iter"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, Molecule& sidechain)->std::optional<std::vector<Molecule>> {
        return rxn.perform_reaction(scaffold, sidechain);
      },
      ""
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, const Set_of_Atoms& scaffold_embedding,
         std::vector<Molecule*> sidechain)->std::optional<Molecule> {
        CheckNoNullMolecules(sidechain, "perform_reaction");
        // Just use default match conditions.
        Sidechain_Match_Conditions smc;
        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, *sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i]->name() << '\n';
            rxn.remove_all_reagents();
            return std::nullopt;
          }
        }
        Molecule result;
        int rc = rxn.perform_reaction(&scaffold, &scaffold_embedding, result);

        rxn.remove_all_reagents();
        if (rc) {
          return result;
        }
        std::cerr << "Cannot react " << scaffold.name() << '\n';
        return std::nullopt;
      },
      "React scaffold with the sidechains - assumes 1 query match per sidechain"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, std::vector<Molecule*> sidechain)->std::optional<Molecule> {
        CheckNoNullMolecules(sidechain, "perform_reaction");

        // Default conditions, multiple matches not allowed.
        Sidechain_Match_Conditions smc;

        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, *sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i]->name() << '\n';
            rxn.remove_all_reagents();
            return std::nullopt;
          }
        }

        Substructure_Results sresults;
        if (rxn.substructure_search(scaffold, sresults) != 1) {
          std::cerr << "perform_reaction:not 1 match to scaffold " << scaffold.name() << '\n';
          rxn.remove_all_reagents();
          return std::nullopt;
        }

        Molecule result;
        int rc = rxn.perform_reaction(&scaffold, sresults.embedding(0), result);

        rxn.remove_all_reagents();
        if (rc) {
          return result;
        }
        std::cerr << "Cannot react " << scaffold.name() << '\n';
        return std::nullopt;

      },
      "React scaffold with sidechains, assuming one substructure match all round"
    )
    .def("perform_reaction_to_list",
      [](IWReaction& rxn, Molecule& scaffold, std::vector<Molecule*> sidechain)->std::vector<Molecule> {
        CheckNoNullMolecules(sidechain, "perform_reaction_to_list");
        std::vector<Molecule> result;

        Sidechain_Match_Conditions smc;
        // Multiple sidechain matches enumerated.
        smc.set_make_new_reagent_for_each_hit(1);

        int number_reagents = 0;
        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, *sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i]->name() << '\n';
            rxn.remove_all_reagents();
            return result;
          }
          number_reagents += rxn.sidechain(0)->number_reagents();
        }

        // Make allowances for 2 scaffold matches. Resizing is expected to be expensive.
        result.reserve(2 * number_reagents);

        Substructure_Results sresults;
        if (rxn.substructure_search(scaffold, sresults) == 0) {
          std::cerr << "perform_reaction:no match to scaffold " << scaffold.name() << '\n';
          rxn.remove_all_reagents();
          return result;
        }

        Reaction_Iterator iter;
        for (iter.initialise(rxn); iter.active(); iter++) {
          Molecule product;
          if (! rxn.perform_reaction(&scaffold, sresults, iter, product)) {
            std::cerr << "Reaction involving " << scaffold.name() << " failed, returning partial result\n";
            rxn.remove_all_reagents();
            return result;
          }
          result.push_back(product);
        }

        rxn.remove_all_reagents();
        return result;
      },
      "For each scaffold embedding, generate list of products"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, std::vector<Molecule*> reagents)->std::optional<Molecule> {
        CheckNoNullMolecules(reagents, "perform_reaction");
        // reagents[0] is the scaffold, so an empty list would read off the end.
        if (reagents.empty()) {
          throw py::value_error("perform_reaction:empty list, reagents[0] must be the scaffold");
        }

        Substructure_Results scaffold_sresults;
        if (rxn.substructure_search(*reagents[0], scaffold_sresults) == 0) {
          std::cerr << "perform_reaction::no match to scaffold " << reagents[0]->name() << '\n';
          return std::nullopt;
        }

        rxn.remove_all_reagents();

        Sidechain_Match_Conditions smc;
        smc.set_ignore_multiple_substucture_matches(1);

        for (size_t i = 1; i < reagents.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i - 1, *reagents[i], smc)) {
            std::cerr << "perform_reaction:cannot add " << reagents[i]->name() << '\n';
            return std::nullopt;
          }
        }

        Molecule product;
        if (! rxn.perform_reaction(reagents[0], scaffold_sresults.embedding(0), product)) {
          std::cerr << "perform_reaction:reaction failed\n";
          return std::nullopt;
        }
          
        return product;
      },
      "Perform the reaction if there is one embedding in the scaffold and in each sidechain"
    )

  ;

  // These do not work very well.
  rxn.def("set_smirks_lost_atom_means_remove_frgment", &set_smirks_lost_atom_means_remove_frgment, "atoms lost in a smirks are removed");
  rxn.def("set_smirks_remove_elements_in_lhs_but_missing_in_rhs", &set_smirks_remove_elements_in_lhs_but_missing_in_rhs,
        "unmapped reagent atoms with elements not in RHS are removed");
}
