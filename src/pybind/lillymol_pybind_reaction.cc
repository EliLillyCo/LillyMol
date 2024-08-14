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

static bool
ReadReaction(const std::string& fname, IWReaction& rxn) {
  IWString tmp(fname.data(), fname.size());
  std::optional<ReactionProto::Reaction> proto =
                iwmisc::ReadTextProto<ReactionProto::Reaction>(tmp);
  if (! proto) {
    return false;
  }

  if (! rxn.ConstructFromProto(*proto, tmp)) {
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

  ;

  py::class_<IWReaction, Substructure_Query>(rxn, "Reaction")
    .def(py::init<>())
    .def("name",
      [](const IWReaction& rxn)->std::string{
        return rxn.comment().AsString();
      },
      "Name"
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
        if (! rxn.ConstructFromTextProto(textproto, dirname))  {
          return false;
        }

        return true;
      },
      ""
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
    .def("remove_no_delete_all_reagents", &IWReaction::remove_no_delete_all_reagents,
         "remove, without destroying, all sidechain reagents"
    )
    .def("substructure_search",
      [](IWReaction& rxn, Molecule& m, Substructure_Results& sresults) {
        return rxn.substructure_search(m, sresults);
      }
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
         std::vector<Molecule>& sidechain)->std::optional<Molecule> {
        // Just use default match conditions.
        Sidechain_Match_Conditions smc;
        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i].name() << '\n';
            rxn.remove_no_delete_all_reagents();
            return std::nullopt;
          }
        }
        Molecule result;
        int rc = rxn.perform_reaction(&scaffold, &scaffold_embedding, result);

        rxn.remove_no_delete_all_reagents();
        if (rc) {
          return result;
        }
        std::cerr << "Cannot react " << scaffold.name() << '\n';
        return std::nullopt;
      },
      "React scaffold with the sidechains - assumes 1 query match per sidechain"
    )
    .def("perform_reaction",
      [](IWReaction& rxn, Molecule& scaffold, std::vector<Molecule>& sidechain)->std::optional<Molecule> {

        // Default conditions, multiple matches not allowed.
        Sidechain_Match_Conditions smc;

        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i].name() << '\n';
            rxn.remove_no_delete_all_reagents();
            return std::nullopt;
          }
        }

        Substructure_Results sresults;
        if (rxn.substructure_search(scaffold, sresults) != 1) {
          std::cerr << "perform_reaction:not 1 match to scaffold " << scaffold.name() << '\n';
          rxn.remove_no_delete_all_reagents();
          return std::nullopt;
        }

        Molecule result;
        int rc = rxn.perform_reaction(&scaffold, sresults.embedding(0), result);

        rxn.remove_no_delete_all_reagents();
        if (rc) {
          return result;
        }
        std::cerr << "Cannot react " << scaffold.name() << '\n';
        return std::nullopt;

      },
      "React scaffold with sidechains, assuming one substructure match all round"
    )
    .def("perform_reaction_to_list",
      [](IWReaction& rxn, Molecule& scaffold, std::vector<Molecule>& sidechain)->std::vector<Molecule> {
        std::vector<Molecule> result;

        Sidechain_Match_Conditions smc;
        // Multiple sidechain matches enumerated.
        smc.set_make_new_reagent_for_each_hit(1);

        int number_reagents = 0;
        for (uint32_t i = 0; i < sidechain.size(); ++i) {
          if (! rxn.add_sidechain_reagent(i, sidechain[i], smc)) {
            std::cerr << "perform_reaction:cannot add sidechain reagent " << sidechain[i].name() << '\n';
            rxn.remove_no_delete_all_reagents();
            return result;
          }
          number_reagents += rxn.sidechain(0)->number_reagents();
        }

        // Make allowances for 2 scaffold matches. Resizing is expected to be expensive.
        result.reserve(2 * number_reagents);

        Substructure_Results sresults;
        if (rxn.substructure_search(scaffold, sresults) == 0) {
          std::cerr << "perform_reaction:no match to scaffold " << scaffold.name() << '\n';
          rxn.remove_no_delete_all_reagents();
          return result;
        }

        Reaction_Iterator iter;
        for (iter.initialise(rxn); iter.active(); iter++) {
          Molecule product;
          if (! rxn.perform_reaction(&scaffold, sresults, iter, product)) {
            std::cerr << "Reaction involving " << scaffold.name() << " failed, returning partial result\n";
            rxn.remove_no_delete_all_reagents();
            return result;
          }
          result.push_back(product);
        }

        rxn.remove_no_delete_all_reagents();
        return result;
      },
      "For each scaffold embedding, generate list of products"
    )
      
  ;

  rxn.def("set_smirks_lost_atom_means_remove_frgment", &set_smirks_lost_atom_means_remove_frgment, "atoms lost in a smirks are removed");
}
