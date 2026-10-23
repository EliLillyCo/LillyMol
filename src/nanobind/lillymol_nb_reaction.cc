// Reaction bindings for the nanobind LillyMol pilot.

#include "nanobind/lillymol_nb_internal.h"

#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Foundational/iwmisc/proto_support.h"
#include "Molecule_Lib/iwreaction.h"

namespace lillymol_nb {
namespace {

bool
ReadReaction(const std::string& fname, IWReaction& reaction) {
  IWString tmp(fname.data(), fname.size());
  std::optional<ReactionProto::Reaction> proto =
      iwmisc::ReadTextProto<ReactionProto::Reaction>(tmp);
  if (!proto) {
    return false;
  }

  Sidechain_Match_Conditions smc;
  return reaction.ConstructFromProto(*proto, tmp, smc);
}

std::string
ReactionName(const IWReaction& reaction) {
  return reaction.comment().AsString();
}

void
SetReactionName(IWReaction& reaction, const std::string& name) {
  reaction.set_comment(name);
}

bool
ConstructReactionFromSmirks(IWReaction& reaction, const std::string& smirks) {
  const_IWSubstring tmp(smirks.data(), smirks.size());
  return reaction.construct_from_smirks(tmp);
}

bool
ConstructReactionFromTextProto(IWReaction& reaction, const std::string& textproto) {
  IWString dirname(".");
  Sidechain_Match_Conditions smc;
  return reaction.ConstructFromTextProto(textproto, dirname, smc);
}

bool
WriteReactionMsi(IWReaction& reaction, const std::string& fname) {
  IWString tmp(fname);
  return reaction.write_msi(tmp.null_terminated_chars());
}

std::optional<std::string>
SidechainName(const IWReaction& reaction, int sidechain_number, int reagent_number) {
  if (sidechain_number < 0 || sidechain_number >= reaction.number_sidechains()) {
    return std::nullopt;
  }

  const Sidechain_Reaction_Site* sidechain = reaction.sidechain(sidechain_number);
  if (sidechain == nullptr) {
    return std::nullopt;
  }

  const Molecule_and_Embedding* reagent = sidechain->reagent(reagent_number);
  if (reagent == nullptr) {
    return std::nullopt;
  }

  return reagent->name().AsString();
}

bool
AddSidechainReagents(IWReaction& reaction, int sidechain, const char* fname,
                     FileType file_type, const Sidechain_Match_Conditions& smc) {
  return reaction.add_sidechain_reagents(sidechain, fname, file_type, smc);
}

bool
AddSidechainReagent(IWReaction& reaction, int sidechain, Molecule& mol,
                    const Sidechain_Match_Conditions& smc) {
  return reaction.add_sidechain_reagent(sidechain, mol, smc);
}

std::vector<std::string>
ReagentNames(const IWReaction& reaction, const Reaction_Iterator& iter) {
  const int n = iter.number_sidechains();
  std::vector<std::string> result;
  result.reserve(n);

  for (int i = 0; i < n; ++i) {
    const int reagent_number = iter.reagent(i);
    const Sidechain_Reaction_Site* sidechain = reaction.sidechain(i);
    if (sidechain == nullptr) {
      continue;
    }
    const Molecule_and_Embedding* reagent = sidechain->reagent(reagent_number);
    if (reagent == nullptr) {
      continue;
    }
    result.emplace_back(reagent->name().AsString());
  }

  return result;
}

std::optional<std::vector<std::vector<int>>>
ReactionSubstructureSearchMatches(IWReaction& reaction, Molecule& mol) {
  Substructure_Results sresults;
  if (!reaction.substructure_search(mol, sresults)) {
    return std::nullopt;
  }

  return SubstructureResultsAsVectors(sresults);
}

bool
ReactionSubstructureSearch(IWReaction& reaction, Molecule& mol,
                           Substructure_Results& sresults) {
  return reaction.substructure_search(mol, sresults);
}

std::optional<Molecule>
PerformReactionWithEmbedding(IWReaction& reaction, Molecule& scaffold,
                             const Set_of_Atoms& embedding) {
  Molecule product;
  if (!reaction.perform_reaction(&scaffold, &embedding, product)) {
    return std::nullopt;
  }
  return product;
}

bool
PerformReactionWithProduct(IWReaction& reaction, Molecule& scaffold,
                           const Set_of_Atoms& embedding, Molecule& product) {
  return reaction.perform_reaction(&scaffold, &embedding, product);
}

std::optional<Molecule>
PerformReactionWithIterator(IWReaction& reaction, const Molecule& scaffold,
                            const Set_of_Atoms& embedding,
                            const Reaction_Iterator& iter) {
  Molecule product;
  if (!reaction.perform_reaction(&scaffold, &embedding, iter, product)) {
    return std::nullopt;
  }
  return product;
}

std::optional<std::vector<Molecule>>
PerformReactionWithSidechain(IWReaction& reaction, Molecule& scaffold, Molecule& sidechain) {
  return reaction.perform_reaction(scaffold, sidechain);
}

// The reaction entry points take std::vector<Molecule*> rather than
// std::vector<Molecule>. With the by value form nanobind copy constructs every
// molecule in the list at the boundary, about 5.4us each, which is a third of the
// cost of a whole reaction - measured on the pybind equivalent, 10.50us passing a
// sidechain by reference against 13.92us passing it in a list. Taking pointers
// copies nothing and does not change the python signature.
//
// Nothing is exposed by not copying: add_sidechain_reagent builds a
// Molecule_and_Embedding from a copy of what it is given, so the caller's
// molecules are not modified.
//
// What pointers do expose is None, which nanobind converts to nullptr, and every
// use below dereferences. One check, so the call sites cannot disagree.
void
CheckNoNullMolecules(const std::vector<Molecule*>& mols, const char* caller) {
  for (uint32_t i = 0; i < mols.size(); ++i) {
    if (mols[i] == nullptr) {
      throw nb::value_error(
          (std::string(caller) + ":None at index " + std::to_string(i)).c_str());
    }
  }
}

std::optional<Molecule>
PerformReactionWithSidechainVectorAndEmbedding(IWReaction& reaction, Molecule& scaffold,
                                               const Set_of_Atoms& scaffold_embedding,
                                               std::vector<Molecule*> sidechains) {
  CheckNoNullMolecules(sidechains, "perform_reaction");
  Sidechain_Match_Conditions smc;
  for (uint32_t i = 0; i < sidechains.size(); ++i) {
    if (!reaction.add_sidechain_reagent(i, *sidechains[i], smc)) {
      std::cerr << "perform_reaction:cannot add sidechain reagent "
                << sidechains[i]->name() << '\n';
      reaction.remove_all_reagents();
      return std::nullopt;
    }
  }

  Molecule result;
  const int rc = reaction.perform_reaction(&scaffold, &scaffold_embedding, result);
  reaction.remove_all_reagents();
  if (rc) {
    return result;
  }

  std::cerr << "Cannot react " << scaffold.name() << '\n';
  return std::nullopt;
}

std::optional<Molecule>
PerformReactionWithSidechainVector(IWReaction& reaction, Molecule& scaffold,
                                   std::vector<Molecule*> sidechains) {
  CheckNoNullMolecules(sidechains, "perform_reaction");
  Sidechain_Match_Conditions smc;
  for (uint32_t i = 0; i < sidechains.size(); ++i) {
    if (!reaction.add_sidechain_reagent(i, *sidechains[i], smc)) {
      std::cerr << "perform_reaction:cannot add sidechain reagent "
                << sidechains[i]->name() << '\n';
      reaction.remove_all_reagents();
      return std::nullopt;
    }
  }

  Substructure_Results sresults;
  if (reaction.substructure_search(scaffold, sresults) != 1) {
    std::cerr << "perform_reaction:not 1 match to scaffold " << scaffold.name() << '\n';
    reaction.remove_all_reagents();
    return std::nullopt;
  }

  Molecule result;
  const int rc = reaction.perform_reaction(&scaffold, sresults.embedding(0), result);
  reaction.remove_all_reagents();
  if (rc) {
    return result;
  }

  std::cerr << "Cannot react " << scaffold.name() << '\n';
  return std::nullopt;
}

std::vector<Molecule>
PerformReactionToList(IWReaction& reaction, Molecule& scaffold,
                      std::vector<Molecule*> sidechains) {
  CheckNoNullMolecules(sidechains, "perform_reaction_to_list");
  std::vector<Molecule> result;

  Sidechain_Match_Conditions smc;
  smc.set_make_new_reagent_for_each_hit(1);

  int number_reagents = 0;
  for (uint32_t i = 0; i < sidechains.size(); ++i) {
    if (!reaction.add_sidechain_reagent(i, *sidechains[i], smc)) {
      std::cerr << "perform_reaction:cannot add sidechain reagent "
                << sidechains[i]->name() << '\n';
      reaction.remove_all_reagents();
      return result;
    }
    const Sidechain_Reaction_Site* sidechain = reaction.sidechain(i);
    if (sidechain != nullptr) {
      number_reagents += sidechain->number_reagents();
    }
  }

  result.reserve(2 * number_reagents);

  Substructure_Results sresults;
  if (reaction.substructure_search(scaffold, sresults) == 0) {
    std::cerr << "perform_reaction:no match to scaffold " << scaffold.name() << '\n';
    reaction.remove_all_reagents();
    return result;
  }

  Reaction_Iterator iter;
  for (iter.initialise(reaction); iter.active(); iter++) {
    Molecule product;
    if (!reaction.perform_reaction(&scaffold, sresults, iter, product)) {
      std::cerr << "Reaction involving " << scaffold.name()
                << " failed, returning partial result\n";
      reaction.remove_all_reagents();
      return result;
    }
    result.push_back(product);
  }

  reaction.remove_all_reagents();
  return result;
}

std::optional<Molecule>
PerformReactionWithReagents(IWReaction& reaction, std::vector<Molecule*> reagents) {
  CheckNoNullMolecules(reagents, "perform_reaction");
  if (reagents.empty()) {
    throw std::invalid_argument("perform_reaction requires at least one reagent");
  }

  Substructure_Results scaffold_sresults;
  if (reaction.substructure_search(*reagents[0], scaffold_sresults) == 0) {
    std::cerr << "perform_reaction::no match to scaffold " << reagents[0]->name() << '\n';
    return std::nullopt;
  }

  reaction.remove_all_reagents();

  Sidechain_Match_Conditions smc;
  smc.set_ignore_multiple_substucture_matches(1);

  for (size_t i = 1; i < reagents.size(); ++i) {
    if (!reaction.add_sidechain_reagent(i - 1, *reagents[i], smc)) {
      std::cerr << "perform_reaction:cannot add " << reagents[i]->name() << '\n';
      return std::nullopt;
    }
  }

  Molecule product;
  if (!reaction.perform_reaction(reagents[0], scaffold_sresults.embedding(0), product)) {
    std::cerr << "perform_reaction:reaction failed\n";
    reaction.remove_all_reagents();
    return std::nullopt;
  }

  reaction.remove_all_reagents();
  return product;
}

}  // namespace

void
BindReaction(nb::module_& m) {
  nb::class_<Sidechain_Match_Conditions>(m, "SidechainMatchConditions")
      .def(nb::init<>())
      .def("set_make_new_reagent_for_each_hit",
           &Sidechain_Match_Conditions::set_make_new_reagent_for_each_hit,
           nb::arg("value"), "Each match generates a regioisomer")
      .def("set_max_matches_to_find",
           &Sidechain_Match_Conditions::set_max_matches_to_find,
           nb::arg("value"), "Set maximum matches")
      .def("set_strip_reagents_to_largest_fragment",
           &Sidechain_Match_Conditions::set_strip_reagents_to_largest_fragment,
           nb::arg("value"), "Use largest fragment")
      .def("set_ignore_not_reacting",
           &Sidechain_Match_Conditions::set_ignore_not_reacting,
           nb::arg("value"), "Ignore non-matching reagents")
      .def("set_find_unique_embeddings_only",
           &Sidechain_Match_Conditions::set_find_unique_embeddings_only,
           nb::arg("value"), "Find unique embeddings")
      .def("set_one_embedding_per_start_atom",
           &Sidechain_Match_Conditions::set_one_embedding_per_start_atom,
           nb::arg("value"), "One embedding per start atom")
      .def("set_ignore_symmetry_related_matches",
           &Sidechain_Match_Conditions::set_ignore_symmetry_related_matches,
           nb::arg("value"), "Ignore symmetry-related matches");

  nb::class_<Reaction_Iterator>(m, "ReactionIterator")
      .def(nb::init<>())
      .def(nb::init<const IWReaction&>(), nb::arg("reaction"))
      .def("initialise", &Reaction_Iterator::initialise, nb::arg("reaction"),
           "Initialise for reaction")
      .def("active", [](const Reaction_Iterator& iter) -> bool { return iter.active(); },
           "True if still active")
      .def("increment", [](Reaction_Iterator& iter) { iter++; }, "Move to next reagent")
      .def("reagent", &Reaction_Iterator::reagent, nb::arg("sidechain"),
           "Get reagent for sidechain")
      .def("reset", &Reaction_Iterator::reset, "Reset the iterator")
      .def("debug_print", [](const Reaction_Iterator& iter) {
             iter.debug_print(std::cerr);
           }, "debugging info");

  nb::class_<IWReaction, Substructure_Query>(m, "Reaction")
      .def(nb::init<>())
      .def("name", &ReactionName, "Name")
      .def("set_name", &SetReactionName, nb::arg("name"), "Assign new reaction name")
      .def("read", [](IWReaction& reaction, const std::string& fname) -> bool {
             return ReadReaction(fname, reaction);
           }, nb::arg("fname"), "Read textproto reaction")
      .def("construct_from_smirks", &ConstructReactionFromSmirks, nb::arg("smirks"),
           "Build from SMIRKS")
      .def("construct_from_textproto", &ConstructReactionFromTextProto,
           nb::arg("textproto"), "Build from a reaction textproto string")
      .def("write_msi", &WriteReactionMsi, nb::arg("fname"))
      .def("sidechain_name", &SidechainName, nb::arg("sidechain_number"),
           nb::arg("reagent_number"),
           "Return the name of a reagent in a sidechain")
      .def("number_sidechains", &IWReaction::number_sidechains, "Number of sidechains")
      .def("number_sidechains_with_reagents",
           &IWReaction::number_sidechains_with_reagents,
           "Number of sidechains with reagents")
      .def("set_one_embedding_per_start_atom", &IWReaction::set_one_embedding_per_start_atom,
           nb::arg("value"), "One embedding per start atom")
      .def("add_sidechain_reagents", &AddSidechainReagents, nb::arg("sidechain"),
           nb::arg("fname"), nb::arg("file_type"), nb::arg("match_conditions"),
           "Add reagent molecules from a file to a sidechain")
      .def("add_sidechain_reagent", &AddSidechainReagent, nb::arg("sidechain"),
           nb::arg("mol"), nb::arg("match_conditions"),
           "Add a reagent molecule to a sidechain")
      // Deliberately not binding remove_no_delete_all_reagents. Reagents added
      // from python go through add_sidechain_reagent, which copies, so the
      // reaction owns them and not deleting them leaks one Molecule_and_Embedding
      // per reagent.
      .def("remove_all_reagents", &IWReaction::remove_all_reagents,
           "Destroy all reagents in all sidechains. Returns the number of sidechains")
      .def("reagent_names", &ReagentNames, nb::arg("iterator"),
           "Return sidechain reagent names at the iterator position")
      .def("substructure_search", &ReactionSubstructureSearch, nb::arg("mol"),
           nb::arg("results"))
      .def("substructure_search_matches", &ReactionSubstructureSearchMatches,
           nb::arg("mol"), "Search scaffold queries and return embeddings")
      .def("in_place_transformations",
           [](IWReaction& reaction, Molecule& mol) -> bool {
             return reaction.in_place_transformations(mol);
           },
           nb::arg("mol"), "Apply reaction to mol")
      .def("perform_reaction", &PerformReactionWithProduct, nb::arg("scaffold"),
           nb::arg("embedding"), nb::arg("product"), "Perform reaction into product")
      .def("perform_reaction", &PerformReactionWithEmbedding, nb::arg("scaffold"),
           nb::arg("embedding"), "Perform reaction with a scaffold embedding")
      .def("perform_reaction", &PerformReactionWithIterator, nb::arg("scaffold"),
           nb::arg("embedding"), nb::arg("iterator"),
           "Generate product based on embedding and iterator")
      .def("perform_reaction", &PerformReactionWithSidechain, nb::arg("scaffold"),
           nb::arg("sidechain"), "React scaffold with one sidechain reagent")
      .def("perform_reaction", &PerformReactionWithSidechainVectorAndEmbedding,
           nb::arg("scaffold"), nb::arg("scaffold_embedding"), nb::arg("sidechains"),
           "React scaffold with sidechains and a specified scaffold embedding")
      .def("perform_reaction", &PerformReactionWithSidechainVector, nb::arg("scaffold"),
           nb::arg("sidechains"),
           "React scaffold with sidechains, requiring one scaffold embedding")
      .def("perform_reaction_to_list", &PerformReactionToList, nb::arg("scaffold"),
           nb::arg("sidechains"), "Generate products for each scaffold embedding")
      .def("perform_reaction", &PerformReactionWithReagents, nb::arg("reagents"),
           "React a list whose first item is scaffold and remaining items are sidechains");

  m.def("set_smirks_lost_atom_means_remove_frgment",
        &set_smirks_lost_atom_means_remove_frgment,
        nb::arg("value"), "Atoms lost in a SMIRKS are removed");
  m.def("set_smirks_remove_elements_in_lhs_but_missing_in_rhs",
        &set_smirks_remove_elements_in_lhs_but_missing_in_rhs,
        nb::arg("value"),
        "Unmapped reagent atoms with elements not in RHS are removed");
}

}  // namespace lillymol_nb
