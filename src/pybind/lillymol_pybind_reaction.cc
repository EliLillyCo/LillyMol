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
    .def("number_sidechains", &IWReaction::number_sidechains, "Number of sidechains")
    .def("number_sidechains_with_reagents", &IWReaction::number_sidechains_with_reagents, "number_sidechains_with_reagents")
    .def("set_one_embedding_per_start_atom", &IWReaction::set_one_embedding_per_start_atom, "one embedding per start atom")
    .def("add_sicechain_reagents",
      [](IWReaction& rxn, int sidechain, const char* fname, FileType file_type, Sidechain_Match_Conditions& smc)->bool{
        return rxn.add_sidechain_reagents(sidechain, fname, file_type, smc);
      },
      "Add reagents to a sidechain"
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
      
  ;
}
