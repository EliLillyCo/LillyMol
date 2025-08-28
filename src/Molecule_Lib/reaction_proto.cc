#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_set>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "iwreaction.h"

using std::cerr;

namespace fs = std::filesystem;

// A set of atoms which are being removed during the reaction.
struct ChangingAtoms {
  resizable_array<uint32_t> atom_being_removed;
};

// Applies to ScaffoldReactionSite and SidechainReactionSite protos.
// Mostly looks for atoms that are associated with some directive that
// are in the set of atoms being deleted.
// Would be nice to have an estimate of the number of atoms in the
// smarts or query so we could check matched atom numbers...
template <typename P>
bool
InternallyConsistent(const P& proto, const ChangingAtoms& changing_atoms) {
  for (const auto& make_bond : proto.make_bond()) {
    if (changing_atoms.atom_being_removed.contains(make_bond.a1())) {
      cerr << "Matched atom " << make_bond.a1() << " in make bond being deleted\n";
      return false;
    }
    if (changing_atoms.atom_being_removed.contains(make_bond.a2())) {
      cerr << "Matched atom " << make_bond.a2() << " in make bond being deleted\n";
      return false;
    }
  }

  for (const auto& change_bond : proto.change_bond()) {
    if (changing_atoms.atom_being_removed.contains(change_bond.a1())) {
      cerr << "Matched atom " << change_bond.a1() << " in change_bond being deleted\n";
      return false;
    }
    if (changing_atoms.atom_being_removed.contains(change_bond.a2())) {
      cerr << "Matched atom " << change_bond.a2() << " in change_bond being deleted\n";
      return false;
    }
  }

  for (const auto& isotope : proto.isotope()) {
    for (uint32_t a : isotope.atom()) {
      if (changing_atoms.atom_being_removed.contains(a)) {
        cerr << "Matched atom " << a << " isotope being deleted\n";
        return false;
      }
    }
  }

  for (const auto& change_isotope : proto.change_isotope()) {
    for (uint32_t a : change_isotope.atom()) {
      if (changing_atoms.atom_being_removed.contains(a)) {
        cerr << "Change isotope atom " << a << " being deleted\n";
        return false;
      }
    }
  }

  for (const auto& change_isotope : proto.change_isotope()) {
    for (uint32_t a : change_isotope.atom()) {
      if (changing_atoms.atom_being_removed.contains(a)) {
        cerr << "Change isotope atom " << a << " being deleted\n";
        return false;
      }
    }
  }

  for (uint32_t a : proto.single_bond()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "Single bond atom " << a << " being deleted\n";
      return false;
    }
  }

  for (uint32_t a : proto.double_bond()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "double bond atom " << a << " being deleted\n";
      return false;
    }
  }

  for (uint32_t a : proto.triple_bond()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "triple bond atom " << a << " being deleted\n";
      return false;
    }
  }

  for (const auto& change_element : proto.change_element()) {
    if (changing_atoms.atom_being_removed.contains(change_element.atom())) {
      cerr << "Change element atom " << change_element.atom() << " being deleted\n";
      return false;
    }
  }

  for (const auto& formal_charge : proto.formal_charge()) {
    if (changing_atoms.atom_being_removed.contains(formal_charge.atom())) {
      cerr << "Formal charge atom " << formal_charge.atom() << " being deleted\n";
      return false;
    }
  }

  for (const auto& change_formal_charge : proto.change_formal_charge()) {
    if (changing_atoms.atom_being_removed.contains(change_formal_charge.atom())) {
      cerr << "Change formal charge atom " << change_formal_charge.atom() << " being deleted\n";
      return false;
    }
  }

  for (uint32_t a : proto.unfix_implicit_hydrogens()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "Unfix implicit hydrogens atom " << a << " being deleted\n";
      return false;
    }
  }

  for (uint32_t a : proto.invert_chirality()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "invert_chirality atom " << a << " being deleted\n";
      return false;
    }
  }

  for (uint32_t a : proto.remove_chirality()) {
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "remove_chirality atom " << a << " being deleted\n";
      return false;
    }
  }

  for (int i = 0; i < proto.atom_and_isotope_size(); i += 2) {
    uint32_t a = proto.atom_and_isotope(i);
    if (changing_atoms.atom_being_removed.contains(a)) {
      cerr << "atom_and_isotope atom " << a << " being deleted\n";
      return false;
    }
  }

  return true;
}

template <typename P>
bool
IdentifyAtomsBeingRemoved(const P& proto,
                          ChangingAtoms& changing_atoms) {
  for (uint32_t a : proto.remove_atom()) {
    if (! changing_atoms.atom_being_removed.add_if_not_already_present(a)) {
      cerr << "Duplicate atom being removed " << a << '\n';
      return false;
    }
  }

  for (uint32_t a : proto.remove_atoms()) {
    if (! changing_atoms.atom_being_removed.add_if_not_already_present(a)) {
      cerr << "Duplicate atom being removed " << a << '\n';
      return false;
    }
  }

  for (uint32_t a : proto.remove_fragment()) {
    if (! changing_atoms.atom_being_removed.add_if_not_already_present(a)) {
      cerr << "Duplicate remove fragment " << a << '\n';
      return false;
    }
  }

  return true;
}

static bool
InternallyConsistent(const ReactionProto::ScaffoldReactionSite& proto,
                ChangingAtoms& changing_atoms) {
  if (!IdentifyAtomsBeingRemoved(proto, changing_atoms)) {
    cerr << proto.ShortDebugString() << '\n';
    return false;
  }

  return InternallyConsistent<ReactionProto::ScaffoldReactionSite>(proto, changing_atoms);
}

static bool
SidechainInternallyConsistent(const ReactionProto::SidechainReactionSite& proto,
                const ChangingAtoms& changing_in_scaffold) {
  ChangingAtoms changing_in_sidechain;

  if (! IdentifyAtomsBeingRemoved(proto, changing_in_sidechain)) {
    cerr << proto.ShortDebugString() << '\n';
    return false;
  }

  if (! InternallyConsistent(proto, changing_in_sidechain)) {
    return false;
  }

  for (const auto& join : proto.join()) {
    if (join.has_a1()) {  // TOo hard to check component specifications.
      uint32_t a1 = join.a1();
      if (changing_in_scaffold.atom_being_removed.contains(a1)) {
        cerr << "Scaffold join atom " << a1 << " being removed\n";
        return false;
      }
    }

    if (join.has_a2()) {
      uint32_t a2 = join.a2();
      if (changing_in_sidechain.atom_being_removed.contains(a2)) {
        cerr << "Sidechain join atom " << a2 << " being removed\n";
        return false;
      }
    }
  }

  // Replace atoms are different. They MUST be removed in the scaffold.
  for (const auto& replace_atom : proto.replace_atom()) {
    if (!replace_atom.has_a1()) {
      uint32_t a1 = replace_atom.a1();
      if (changing_in_scaffold.atom_being_removed.contains(a1)) {
        cerr << "Scaffold replace atom " << a1 << " NOT being removed\n";
        return false;
      }
    }
    if (replace_atom.has_a2()) {
      uint32_t a2 = replace_atom.a2();
      if (changing_in_sidechain.atom_being_removed.contains(a2)) {
        cerr << "Sidechain replace atom " << a2 << " being removed\n";
        return false;
      }
    }
  }

  return true;
}

static bool
InternallyConsistent(const ReactionProto::Reaction& proto) {
  if (! proto.has_scaffold()) {
    cerr << "InternallyConsistent:no scaffold\n";
    cerr << proto.ShortDebugString() << '\n';
    return false;
  }

  ChangingAtoms changing_atoms;

  if (! InternallyConsistent(proto.scaffold(), changing_atoms)) {
    cerr << proto.scaffold().ShortDebugString() << '\n';
    return false;
  }

  for (const auto& sidechain : proto.sidechain()) {
    if (! SidechainInternallyConsistent(sidechain, changing_atoms)) {
      cerr << sidechain.ShortDebugString() << '\n';
      return false;
    }
  }

  return true;
}

template <typename P>
int
WriteError(const char * message, const P & proto)
{
  cerr << message << ' ' << proto.ShortDebugString() << "\n";

  return 0;
}

template <typename P>
int
Match_Conditions::ConstructFromProto(const P& proto) {
  if (proto.has_ignore_not_reacting())
    _ignore_not_reacting = proto.ignore_not_reacting();
  if (proto.has_find_unique_embeddings())
    _find_unique_embeddings = proto.find_unique_embeddings();
  if (proto.has_process_hit_number())
    _process_hit_number = proto.process_hit_number();
  if (proto.has_one_embedding_per_start_atom())
    _one_embedding_per_start_atom = proto.one_embedding_per_start_atom();
  if (proto.has_ignore_symmetry_related_matches())
    _ignore_symmetry_related_matches = proto.ignore_symmetry_related_matches();
  if (proto.has_multiple_match_string())
    _multiple_match_string = proto.multiple_match_string();
  if (proto.has_suppress_if_more_than_this_many_substructure_search_hits())
    _suppress_if_more_than_this_many_substructure_search_hits = proto.suppress_if_more_than_this_many_substructure_search_hits();
  if (proto.has_embeddings_can_overlap()) {
    _embeddings_can_overlap = proto.embeddings_can_overlap();
  }
  if (proto.has_max_matches_to_find()) {
    _max_matches_to_find = proto.max_matches_to_find();
  }

  return 1;
}

int
Scaffold_Match_Conditions::ConstructFromProto(const ReactionProto::ScaffoldMatchConditions& proto) {
  Match_Conditions::ConstructFromProto(proto);

  _enumerate_scaffold_hits_individually  = proto.enumerate_scaffold_hits_individually();
  _combinatorial_expansion_of_scaffold_hits = proto.combinatorial_expansion_of_scaffold_hits();
  _max_matches_to_find = proto.max_matches_to_find();

  return 1;
}

int
Sidechain_Match_Conditions::ConstructFromProto(const ReactionProto::SidechainMatchConditions& proto) {
  Match_Conditions::ConstructFromProto(proto);

  _make_new_reagent_for_each_hit = proto.make_new_reagent_for_each_hit();
  _max_matches_to_find = proto.max_matches_to_find();
  _strip_reagents_to_largest_fragment = proto.strip_reagents_to_largest_fragment();

  return 1;
}

template <typename P>
int
BondFromProto(Bond& b, const P& proto, const bool self_bonds_ok = false,
              const bool must_have_btype = true)
{
  if (! proto.has_a1())
    return WriteError("BondFromProto:no a1", proto);
  if (! proto.has_a2())
    return WriteError("BondFromProto:no a2", proto);
  if (must_have_btype && ! proto.has_btype())
    return WriteError("BondFromProto:no btype", proto);
  if (! self_bonds_ok && proto.a1() == proto.a2())
    return WriteError("BondFromProto:self bonds not allowed", proto);
  b.set_a1(proto.a1());
  b.set_a2(proto.a2());
  if (proto.has_btype())
  {
    switch (proto.btype())
    {
      case SubstructureSearch::SS_SINGLE_BOND:
        b.set_bond_type(SINGLE_BOND);
        break;
      case SubstructureSearch::SS_DOUBLE_BOND:
        b.set_bond_type(DOUBLE_BOND);
        break;
      case SubstructureSearch::SS_TRIPLE_BOND:
        b.set_bond_type(TRIPLE_BOND);
        break;
      default:
        return WriteError("BondFromProto:unrecognized bond type", proto);
    }
  }

//cerr << "BondFromProto "; b.debug_print(cerr);

  return 1;
}

int
BondFromPairOfAtoms(Bond& b, const ReactionProto::PairOfAtoms& proto, const bool self_bonds_ok = false)
{
  if (! proto.has_a1())
    return WriteError("BondFromPairOfAtoms:no a1", proto);
  if (! proto.has_a2())
    return WriteError("BondFromPairOfAtoms:no a2", proto);
  if (! self_bonds_ok && proto.a1() == proto.a2())
    return WriteError("BondFromPairOfAtoms:self bonds not allowed", proto);

  b.set_a1(proto.a1());
  b.set_a2(proto.a2());

  return 1;
}

int
Reaction_Stereo_Centre::ConstructFromProto(const ReactionProto::StereoCenter& proto)
{
  if (! proto.has_center())
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:no center", proto);
  if (! proto.has_top_front())
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:no top_front", proto);
  if (! proto.has_top_back())
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:no top_back", proto);
  if (! proto.has_left_down())
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:no left_down", proto);
  if (! proto.has_right_down())
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:no right_down", proto);

  if (! _ssc[0].ConstructFromProto(proto.center()) ||
      ! _ssc[1].ConstructFromProto(proto.top_front()) ||
      ! _ssc[2].ConstructFromProto(proto.top_back()) ||
      ! _ssc[3].ConstructFromProto(proto.left_down()) ||
      ! _ssc[4].ConstructFromProto(proto.right_down()))
    return WriteError("Reaction_Stereo_Centre::ConstructFromProto:invalid info", proto);

  return 1;
}

int
Stereo_Centre_Component::ConstructFromProto(const ReactionProto::StereoCenterComponent& proto)
{
  if (proto.has_atom()) {
    _implicit_hydrogen = false;
    if (! Matched_Atom_in_Component::ConstructFromProto(proto.atom()))
      return WriteError("Stereo_Centre_Component::ConstructFromProto:invalid atom", proto);
  } else {  // We do not guard against _implicit_hydrogen: false
    _implicit_hydrogen = true;
  }

  return 1;
}

int
Matched_Atom_in_Component::ConstructFromProto(const ReactionProto::MatchedAtomInComponent& proto)
{
  if (proto.has_component_and_atom())
  {
    const const_IWSubstring tmp(proto.component_and_atom().data(), proto.component_and_atom().size());
    return _construct_from_string(proto.component_and_atom());
  }

  if (!proto.has_atom())
    return WriteError("Reaction_Change_Element::ConstructFromProto:no atom", proto);
  if (!proto.has_component())
    return WriteError("Reaction_Change_Element::ConstructFromProto:no component", proto);
  if (proto.component() < 0)
    return WriteError("Reaction_Change_Element::ConstructFromProto:invalid component", proto);

  _matched_atom = proto.atom();
  _component = proto.component() - 1;  // Components are offset, so -1 is the scaffold.

  return 1;
}

// Parse the component.atom specification.
// Note that historical usage was with a '.' separator, but it should be ':'. Allow both.

int
Matched_Atom_in_Component::_construct_from_string(const const_IWSubstring& s)
{
  char sep = '.';
  if (s.contains(':')) {
    sep = ':';
  }

  int i = 0;
  const_IWSubstring c, a;
  if (! s.nextword(c, i, sep) || c.empty() ||
      ! s.nextword(a, i, sep) || a.empty())
  {
    cerr << "Matched_Atom_in_Component::_construct_from_string:invalid spec '" << s << "'\n";
    return 0;
  }

  if (! c.numeric_value(_component) || _component < 0 ||
      ! a.numeric_value(_matched_atom) || _matched_atom < 0)
  {
    cerr << "Matched_Atom_in_Component::_construct_from_string:invalid vaues '" << s << "'\n";
    return 0;
  }

  _component--;

  return 1;
}

int
Reaction_Change_Element::ConstructFromProto(const ReactionProto::ChangeElement& proto)
{
  if (!proto.has_atom())
    return WriteError("Reaction_Change_Element::ConstructFromProto:no atom", proto);
  if (!proto.has_element())
    return WriteError("Reaction_Change_Element::ConstructFromProto:no element", proto);
  const IWString mystring(proto.element());

  const Element * e = get_element_from_symbol_no_case_conversion(mystring);
  if (nullptr == e)
    return WriteError("Reaction_Change_Element::ConstructFromProto:invalid element", proto);

  _atom = proto.atom();
  _element = e;

  return 1;
}

int
Reaction_Formal_Charge::ConstructFromProto(const ReactionProto::FormalCharge& proto)
{
  if (! proto.has_atom())
    return WriteError("Reaction_Formal_Charge::ConstructFromProto:no atom", proto);
  if (! proto.has_formal_charge())
    return WriteError("Reaction_Formal_Charge::ConstructFromProto:no formal_charge", proto);
  _atom = proto.atom();
  _fc = proto.formal_charge();  // Should we check for reasonablness?

  return 1;
}

int
Reaction_Change_Formal_Charge::ConstructFromProto(const ReactionProto::ChangeFormalCharge& proto)
{
  if (! proto.has_atom())
    return WriteError("Reaction_Change_Formal_Charge::ConstructFromProto:no atom", proto);
  if (! proto.has_delta())
    return WriteError("Reaction_Change_Formal_Charge::ConstructFromProto:no delta", proto);
  _atom = proto.atom();
  _delta = proto.delta();  // Should we check for reasonablness?

  return 1;
}

int
Reaction_Place_Isotope::ConstructFromProto(const ReactionProto::PlaceIsotope& proto)
{
  if (proto.atom_size() == 0)
    return WriteError("Reaction_Place_Isotope::ConstructFromProto:no atom", proto);
  if (! proto.has_isotope())
    return WriteError("Reaction_Place_Isotope::ConstructFromProto:no isotope", proto);

  for (int a : proto.atom()) {
    _atom << a;
  }

  _isotope = proto.isotope();

  return 1;
}

int
Reaction_Increment_Isotope::ConstructFromProto(const ReactionProto::IncrementIsotope& proto)
{
  if (proto.atom().empty()) {
    return WriteError("Reaction_Increment_Isotope::ConstructFromProto:no atom", proto);
  }

  int properties_found = 0;

  for (int m : proto.atom()) {
    _atom << m;
  }

  if (proto.has_delta()) {
    _isotope = proto.delta();
    ++properties_found;
  }
  if (proto.has_multiply()) {
    _multiply = proto.multiply();
    ++properties_found;
  }
  if (proto.has_integer_divide()) {
    _integer_divide = proto.integer_divide();
    ++properties_found;
  }
  if (proto.has_modulo()) {
    _modulo = proto.modulo();
    ++properties_found;
  }

  if (properties_found == 0) {
    cerr << "Reaction_Increment_Isotope:ConstructFromProto:no increments\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_only_change_if()) {
    if (proto.only_change_if().has_values()) {
      if (! BuildOnlyChangeIf(proto.only_change_if().values())) {
        cerr << "Reaction_Increment_Isotope::ConstructFromProto:invalid only change if\n";
        cerr << proto.ShortDebugString() << '\n';
        return 0;
      }
    } else if (proto.only_change_if().has_range()) {
      if (! BuildOnlyChangeIf(proto.only_change_if().range())) {
        cerr << "Reaction_Increment_Isotope::ConstructFromProto:invalid only change if\n";
        cerr << proto.ShortDebugString() << '\n';
        return 0;
      }
    } else {
      return 0;
    }
  }

  return 1;
}

int
Reaction_Increment_Isotope::BuildOnlyChangeIf(const ReactionProto::Range& proto) {
  if (! proto.has_min() || ! proto.has_max() || proto.min() > proto.max()) {
    cerr << "ReactionProto::BuildOnlyChangeIf:invalid proto " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _only_change_if.set_min(proto.min());
  _only_change_if.set_max(proto.max());

  return 1;
}

int
Reaction_Increment_Isotope::BuildOnlyChangeIf(const ReactionProto::Values& proto) {
  if (proto.values().empty()) {
    cerr << "ReactionProto::BuildOnlyChangeIf:invalid proto " << proto.ShortDebugString() << '\n';
    return 0;
  }

  for (int a : proto.values()) {
    if (! _only_change_if.add_if_not_already_present(a)) {
      cerr << "ReactionProto::BuildOnlyChangeIf:duplicate value " << a << '\n';
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Reaction_Invert_Isotope::ConstructFromProto(const ReactionProto::PlaceIsotope& proto)
{
  if (proto.atom().empty())
    return WriteError("Reaction_Invert_Isotope::ConstructFromProto:no atom", proto);
  if (! proto.has_isotope())
    return WriteError("Reaction_Invert_Isotope::ConstructFromProto:no isotope", proto);

  for (int m : proto.atom()) {
    _atom << m;
  }

  _isotope = proto.isotope();

  return 1;
}

int
Reaction_Bond_Length::ConstructFromProto(const ReactionProto::BondLength& proto,
     const int component)
{
  if (! proto.has_a1() && !proto.has_c1())
    return WriteError("Reaction_Bond_Length::ConstructFromProto:no a1/c1", proto);
  if (! proto.has_a2() && !proto.has_c2())
    return WriteError("Reaction_Bond_Length::ConstructFromProto:no a2/c2", proto);
  if (! proto.has_distance())
    return WriteError("Reaction_Bond_Length::ConstructFromProto:no distance", proto);
  if (proto.distance() < 0.0f)
    return WriteError("Reaction_Bond_Length::ConstructFromProto:invalid distance", proto);

  if (proto.has_a1())
  {
    _atom[0].set_matched_atom(proto.a1());
    _atom[0].set_in_component(component - 1);
  }
  else
  {
    if (! _atom[0].ConstructFromProto(proto.c1())) {
      return WriteError("Reaction_Bond_Length::ConstructFromProto:invalid c1", proto);
    }
  }

  if (proto.has_a2()) {
    _atom[1].set_matched_atom(proto.a2());
    _atom[1].set_in_component(component - 1);
  } else if (proto.has_c2()) {
    if (! _atom[1].ConstructFromProto(proto.c2())) {
      return WriteError("Reaction_Bond_Length::ConstructFromProto:invalid c2", proto);
    }
  }

  _desired_length = proto.distance();

  return 1;
}

int
Reaction_Bond_Length::BuildProto(ReactionProto::BondLength& proto) const {
  if (_atom[0].in_scaffold()) {
    proto.set_a1(_atom[0].atom());
  } else {
    auto* c = proto.mutable_c1();
    c->set_component(_atom[0].in_component());
    c->set_atom(_atom[0].atom());
  }

  if (_atom[1].in_scaffold()) {
    proto.set_a2(_atom[1].atom());
  } else {
    auto* c = proto.mutable_c2();
    c->set_component(_atom[1].in_component());
    c->set_atom(_atom[1].atom());
  }

  proto.set_distance(_desired_length);

  return 1;
}

int
Reaction_Bond_Angle::ConstructFromProto(const ReactionProto::BondAngle& proto,
    const int component)
{
  if (! proto.has_a1() && ! proto.has_c1())
    return WriteError("Reaction_Bond_Angle::ConstructFromProto:no a1/c1", proto);
  if (! proto.has_a2() && ! proto.has_c2())
    return WriteError("Reaction_Bond_Angle::ConstructFromProto:no a2/c2", proto);
  if (! proto.has_a3() && ! proto.has_c3())
    return WriteError("Reaction_Bond_Angle::ConstructFromProto:no a3/c3", proto);
  if (! proto.has_angle())
    return WriteError("Reaction_Bond_Angle::ConstructFromProto:no angle", proto);

  if (proto.has_a1()) {
    _atom[0].set_matched_atom(proto.a1());
    _atom[0].set_in_component(component - 1);
  } else if (proto.has_c1()) {
    if (! _atom[0].ConstructFromProto(proto.c1()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c1", proto);
  }

  if (proto.has_a2()) {
    _atom[1].set_matched_atom(proto.a2());
    _atom[1].set_in_component(component - 1);
  } else if (proto.has_c2()) {
    if (! _atom[1].ConstructFromProto(proto.c2()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c2", proto);
  }

  if (proto.has_a3()) {
    _atom[2].set_matched_atom(proto.a3());
    _atom[2].set_in_component(component - 1);
  } else if (proto.has_c3()) {
    if (! _atom[2].ConstructFromProto(proto.c3()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c3", proto);
  }

  _desired_angle = proto.angle() * DEG2RAD;

  return 1;
}

int
Reaction_Bond_Angle::BuildProto(ReactionProto::BondAngle& proto) const {
  if (_atom[0].in_scaffold()) {
    proto.set_a1(_atom[0].atom());
  } else {
    auto* c = proto.mutable_c1();
    c->set_component(_atom[0].in_component());
    c->set_atom(_atom[0].atom());
  }

  if (_atom[1].in_scaffold()) {
    proto.set_a2(_atom[1].atom());
  } else {
    auto* c = proto.mutable_c2();
    c->set_component(_atom[1].in_component());
    c->set_atom(_atom[1].atom());
  }

  if (_atom[2].in_scaffold()) {
    proto.set_a3(_atom[2].atom());
  } else {
    auto* c = proto.mutable_c3();
    c->set_component(_atom[2].in_component());
    c->set_atom(_atom[2].atom());
  }

  proto.set_angle(_desired_angle);

  return 1;
}

int
Reaction_Dihedral_Angle::ConstructFromProto(const ReactionProto::DihedralAngle& proto,
    const int component)
{
  if (! proto.has_a1() && ! proto.has_c1())
    return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:no a1/c1", proto);
  if (! proto.has_a2() && ! proto.has_c2())
    return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:no a2/c2", proto);
  if (! proto.has_a3() && ! proto.has_c3())
    return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:no a3/c3", proto);
  if (! proto.has_a4() && ! proto.has_c4())
    return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:no a4/c4", proto);
  if (! proto.has_angle())
    return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:no angle", proto);

  if (proto.has_a1()) {
    _atom[0].set_matched_atom(proto.a1());
    _atom[0].set_in_component(component - 1);
  } else if (proto.has_c1()) {
    if (! _atom[0].ConstructFromProto(proto.c1()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c1", proto);
  }

  if (proto.has_a2()) {
    _atom[1].set_matched_atom(proto.a2());
    _atom[1].set_in_component(component - 1);
  } else if (proto.has_c2()) {
    if (! _atom[1].ConstructFromProto(proto.c2()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c2", proto);
  }

  if (proto.has_a3()) {
    _atom[2].set_matched_atom(proto.a3());
    _atom[2].set_in_component(component - 1);
  } else if (proto.has_c3()) {
    if (! _atom[2].ConstructFromProto(proto.c3()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c3", proto);
  }

  if (proto.has_a4()) {
    _atom[3].set_matched_atom(proto.a4());
    _atom[3].set_in_component(component - 1);
  } else if (proto.has_c4()) {
    if (! _atom[3].ConstructFromProto(proto.c4()))
      return WriteError("Reaction_Dihedral_Angle::ConstructFromProto:invalid c4", proto);
  }

  _desired_angle = proto.angle() * DEG2RAD;

  return 1;
}

int
Reaction_Dihedral_Angle::BuildProto(ReactionProto::DihedralAngle& proto) const {
  if (_atom[0].in_scaffold()) {
    proto.set_a1(_atom[0].atom());
  } else {
    auto* c = proto.mutable_c1();
    c->set_component(_atom[0].in_component());
    c->set_atom(_atom[0].atom());
  }

  if (_atom[1].in_scaffold()) {
    proto.set_a2(_atom[1].atom());
  } else {
    auto* c = proto.mutable_c2();
    c->set_component(_atom[1].in_component());
    c->set_atom(_atom[1].atom());
  }

  if (_atom[2].in_scaffold()) {
    proto.set_a3(_atom[2].atom());
  } else {
    auto* c = proto.mutable_c3();
    c->set_component(_atom[2].in_component());
    c->set_atom(_atom[2].atom());
  }

  if (_atom[3].in_scaffold()) {
    proto.set_a4(_atom[3].atom());
  } else {
    auto* c = proto.mutable_c4();
    c->set_component(_atom[3].in_component());
    c->set_atom(_atom[3].atom());
  }

  proto.set_angle(_desired_angle);

  return 1;
}

int
Reaction_3D_Replace::ConstructFromProto(const ReactionProto::CoordinateTransfer& proto,
                int component)
{
  if (proto.atoms().empty()) {
    return WriteError("Reaction_3D_Replace::ConstructFromProto:no atoms", proto);
  }

  _n = proto.atoms_size();
  _a1 = new Matched_Atom_in_Component[_n];
  _a2 = new Matched_Atom_in_Component[_n];

  for (int i = 0; i < _n; ++i)
  {
    const ReactionProto::ComponentAndAtom& atom = proto.atoms(i);
    if (! _a1[i].ConstructFromProto(atom.c1())) {
      return WriteError("Reaction_3D_Replace::ConstructFromProto:invalid a1", proto);
    }
    _a2[i].set_in_component(component - 1);
    _a2[i].set_matched_atom(atom.a2());
    cerr << "Reaction_3D_Replace:i " << i << " component " << component << " atom " << atom.a2() << '\n';
  }

  _weight = new double[_n];

  _weight[0] = 1.0;

  for (int i = 1; i < _n; i++)
  {
    _weight[i] = 0.1;
  }

  return 1;
}

int
Reaction_3D_Replace::BuildProto(ReactionProto::CoordinateTransfer& proto) const {
  cerr << "ReactionProto::BuildProto:not implemented\n";
  return 1;
}

int
Reaction_Wedge_Bond::ConstructFromProto(const ReactionProto::WedgeBond& proto)
{
  if (! proto.has_a1())
    return WriteError("Reaction_Wedge_Bond::ConstructFromProto:no a1", proto);
  if (! proto.has_a2())
    return WriteError("Reaction_Wedge_Bond::ConstructFromProto:no a2", proto);
  if (! proto.has_direction())
    return WriteError("Reaction_Wedge_Bond::ConstructFromProto:no direction", proto);

  _a1 = proto.a1();
  _a2 = proto.a2();
  _direction = proto.direction();

  return 1;
}

int
Replace_Atom::ConstructFromProto(const ReactionProto::ReplaceAtom& proto,
    const int component)
{
  if (! proto.has_a1() && ! proto.has_c1())
    return WriteError("Replace_Atom::ConstructFromProto:no a1/c1", proto);
  if (! proto.has_a2() && ! proto.has_c2())
    return WriteError("Replace_Atom::ConstructFromProto:no a2/c2", proto);

  if (proto.has_a1()) {
    _a1.set_matched_atom(proto.a1());
    _a1.set_in_component(-1);  // In scaffold.
  } else if (proto.has_c1()) {
    if (! _a1.ConstructFromProto(proto.c1()))
      return WriteError("Replace_Atom::ConstructFromProto:invalid c1", proto);
  }

  if (proto.has_a2()) {
    _a2.set_matched_atom(proto.a2());
    _a2.set_in_component(component - 1);
  } else if (proto.has_c2()) {
    if (! _a2.ConstructFromProto(proto.c2()))
      return WriteError("Replace_Atom::ConstructFromProto:invalid c2", proto);
  }

  if (proto.has_invert_chirality()) {
    _invert_chirality = proto.invert_chirality();
  }

  return 1;
}

int
Inter_Particle_Bond::ConstructFromProto(const ReactionProto::InterParticleBond& proto,
    const int component)
{
  if (! proto.has_a1() && ! proto.has_c1())
    return WriteError("Inter_Particle_Bond::ConstructFromProto:no a1/c1", proto);
  if (! proto.has_a2() && ! proto.has_c2())
    return WriteError("Inter_Particle_Bond::ConstructFromProto:no a2/c2", proto);

  if (proto.has_a1()) {
    _a1.set_matched_atom(proto.a1());
    _a1.set_in_component(-1);  // Scaffold
  } else if (proto.has_c1()) {
    if (! _a1.ConstructFromProto(proto.c1()))
      return WriteError("Inter_Particle_Bond::ConstructFromProto:invalid c1", proto);
  }

  if (proto.has_a2()) {
    _a2.set_matched_atom(proto.a2());
    _a2.set_in_component(component - 1);
  } else if (proto.has_c2()) {
    if (! _a2.ConstructFromProto(proto.c2()))
      return WriteError("Inter_Particle_Bond::ConstructFromProto:invalid c2", proto);
  }

  // The default bond type is a single bond.
  if (! proto.has_btype()) {
    _bt = SINGLE_BOND;
  } else {
    switch (proto.btype())
    {
      case SubstructureSearch::SS_SINGLE_BOND:
        _bt = SINGLE_BOND;
        break;
      case SubstructureSearch::SS_DOUBLE_BOND:
        _bt = DOUBLE_BOND;
        break;
      case SubstructureSearch::SS_TRIPLE_BOND:
        _bt = TRIPLE_BOND;
        break;
      default:
        return WriteError("Inter_Particle_Bond::ConstructFromProto:unrecognized bond type", proto);
    }
  }

  if (proto.has_max_distance()) {
    _distance = proto.max_distance();
  }

  if (proto.isotope_join_requirement() == ReactionProto::JoinRequirement::NONE) {
  } else if (proto.isotope_join_requirement() == ReactionProto::JoinRequirement::ATOMIC_NUMBER) {
    _sidechain_isotope_requirement = SidechainIsotopeRequirement::kAtomicNumber;
  } else if (proto.isotope_join_requirement() == ReactionProto::JoinRequirement::ISOTOPE) {
    _sidechain_isotope_requirement = SidechainIsotopeRequirement::kIsotope;
  } else if (proto.isotope_join_requirement() == ReactionProto::JoinRequirement::ATYPE) {
    _sidechain_isotope_requirement = SidechainIsotopeRequirement::kAtype;
  } else {
    cerr << "Inter_Particle_Bond::ConstructFromProto:unrecognised sidechain isotope requirement\n";
    return 0;
  }

  if (proto.has_align_3d()) {
    _align_3d = proto.align_3d();
  }

  return 1;
}

// If the proto has a comment field, copy it to `destination`.
template <typename P>
void
MaybeCopyComment(const P& proto,
                 IWString& destination) {
  if (proto.comment_size() == 0) {
    return;
  }

  // Should we insert a separator between entries?
  if (proto.comment_size()) {
    for (const std::string& s : proto.comment()) {
      destination << s;
    }
  }
}

int
No_Reaction::ConstructFromProto(const ReactionProto::NoReaction& proto)
{
  MaybeCopyComment(proto, _comment);

  if (proto.has_scaffold_no_reaction() &&
      ! _scaffold_no_reaction.ConstructFromProto(proto.scaffold_no_reaction()))
     return WriteError("No_Reaction::ConstructFromProto:invalid scaffold no reaction query", proto);
  if (proto.has_sidechain_no_reaction() &&
      ! _sidechain_no_reaction.ConstructFromProto(proto.sidechain_no_reaction()))
     return WriteError("No_Reaction::ConstructFromProto:invalid sidechain no reaction query", proto);

  return 1;
}

// Return a file name to use for the query_file directive.
// If a file with the name `query_file_name` exists already, return it.
// Otherwise, check for existence relative to the directory of `proto_file_name`.
std::optional<IWString>
PathNameOfQueryFile(const IWString& proto_file_name,
                   const IWString& query_file_name_iws) {
  const std::string query_file_name(query_file_name_iws.data(), query_file_name_iws.length());

  if (fs::exists(fs::path(query_file_name))) {
    // cerr << "File name " << query_file_name << " alread exists, good\n";
    return query_file_name_iws;
  }

  // Try query_file_name in same directory as proto_file_name.

  // cerr << "proto_file_name '" << proto_file_name << "'\n";
  fs::path new_name = fs::path(std::string(proto_file_name.data(), proto_file_name.length())).replace_filename(query_file_name);
  // cerr << "With directory '" << new_name.string() << "'\n";

  if (fs::exists(new_name)) {
    return new_name.string();
  }

  return std::nullopt;
}

// Return a proto specified by a query_file directive in `owning_proto_fname`.
// Shell variables are expanded, and PathNameOfQueryFile is called to possibly
// find the file.
std::optional<SubstructureSearch::SubstructureQuery>
ProtoFromQueryFile(const std::string& qfile,
                   const IWString& owning_proto_fname) {
  IWString query_file(qfile);
  std::optional<IWString> maybe_expanded = query_file.ExpandEnvironmentVariables();
  if (! maybe_expanded) {
    cerr << "ProtoFromQueryFile:shell variable expansion failed " << query_file << '\n';
    return std::nullopt;
  }

  std::optional<IWString> maybe_path = PathNameOfQueryFile(owning_proto_fname, *maybe_expanded);
  if (!maybe_path) {
    cerr << "ProtoFromQueryFile:cannot find query_file " << *maybe_expanded << '\n';
    return std::nullopt;
  }

  std::optional<SubstructureSearch::SubstructureQuery> maybe_query =
            iwmisc::ReadTextProto<SubstructureSearch::SubstructureQuery>(*maybe_path);
  if (! maybe_query) {
    cerr << "ProtoFromQueryFile:cannot parse query file\n";
    return std::nullopt;
  }

  return *maybe_query;
}

// Query constraints and queries have been read, and the _match_conditions
// variable has been filled. Some values from match_conditions need
// to be propagated to the Substructure_Query objects.
template <typename M>
int
Reaction_Site::InitialiseQueryConstraints(const M& match_conditions) {
  for (Substructure_Query* q : _queries) {
    if (match_conditions.find_unique_embeddings_only()) {
      q->set_find_unique_embeddings_only(1);
    }
    if (match_conditions.one_embedding_per_start_atom()) {
      q->set_find_one_embedding_per_atom(1);
    }
    if (match_conditions.ignore_symmetry_related_matches()) {
      q->set_do_not_perceive_symmetry_equivalent_matches(1);
    }
    if (! match_conditions.embeddings_can_overlap()) {
      q->set_embeddings_do_not_overlap(1);
    }
    if (auto s = match_conditions.max_matches_to_find(); s) {
      q->set_max_matches_to_find(s);
    }
  }

  return 1;
}

template <typename P>
int
Reaction_Site::ConstructFromProto(const P& proto, const IWString& fname)
{
  MaybeCopyComment(proto, _comment);

  if (proto.query_size() > 0) {
    for (const SubstructureSearch::SubstructureQuery& query : proto.query()) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->ConstructFromProto(query)) {
        return WriteError("Reaction_Site::ConstructFromProto:invalid query", proto);
      }
      _queries << q.release();
    }
  }

  if (proto.query_file_size() > 0) {
    for (const std::string& qfile : proto.query_file()) {
      std::optional<SubstructureSearch::SubstructureQuery> maybe_query = 
                ProtoFromQueryFile(qfile, fname);
      if (! maybe_query) {
        return WriteError("Reaction_Site::ConstructFromProto:cannot get query_file proto ", proto);

      }

      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (!q->ConstructFromProto(*maybe_query)) {
        return WriteError("Reaction_Site::ConstructFromProto:cannot parse query ", *maybe_query);
      }
      _queries.add(q.release());
    }
  }

  if (proto.smarts().size() > 0) {
    for (const std::string& smt : proto.smarts()) {
      IWString smarts(smt);
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        return WriteError("Scaffold_Reaction_Site::ConstructFromProto:invalid smarts", proto);
      }
      _queries << q.release();
    }
  }

  for (const auto& bond_to_be_made : proto.make_bond()) {
    std::unique_ptr<Bond> b(new Bond());
    if (! BondFromProto(*b, bond_to_be_made))
      return WriteError("Reaction_Site::ConstructFromProto:invalid bond to be made ", bond_to_be_made);
    _bonds_to_be_made.add(b.release());
  }

  for (const auto& bond_to_be_changed : proto.change_bond()) {
    std::unique_ptr<Bond> b(new Bond());
    if (! BondFromProto(*b, bond_to_be_changed))
      return WriteError("Reaction_Site::ConstructFromProto:invalid bond to be changed ", bond_to_be_changed);
    _bonds_to_be_changed.add(b.release());
  }

  for (const auto& bond_to_be_broken : proto.break_bond()) {
    std::unique_ptr<Bond> b(new Bond());
    b->set_bond_type(SINGLE_BOND);   // Seems invalid if not set.
    if (! BondFromPairOfAtoms(*b, bond_to_be_broken, /*no self bonds*/false))
      return WriteError("Reaction_Site::ConstructFromProto:invalid bond to be made ", bond_to_be_broken);
    _bonds_to_be_broken.add(b.release());
  }

  if (proto.break_bond_between_atoms_size() > 0) {
    if (proto.break_bond_between_atoms_size() % 2 != 0) {
      return WriteError("Reaction_Site::ConstructFromProto:must have even number of break_bond_between_atoms ",
                         proto);
    }

    for (int i = 0; i < proto.break_bond_between_atoms_size(); i += 2) {
      uint32_t a1 = proto.break_bond_between_atoms(i);
      uint32_t a2 = proto.break_bond_between_atoms(i + 1);
      // The bond type is not used.
      _bonds_to_be_broken << new Bond(a1, a2, SINGLE_BOND);
    }
  }

  for (const auto atom :proto.remove_atom()) {
    _atoms_to_be_removed.add_if_not_already_present(atom);
  }

  for (const auto atom: proto.remove_atoms()) {
    _atoms_to_be_removed.add_if_not_already_present(atom);
  }

  for (const auto atom :proto.remove_fragment()) {
    _fragments_to_be_removed.add_if_not_already_present(atom);
  }

  for (const auto atom :proto.keep_fragment()) {
    _fragments_to_be_kept.add_if_not_already_present(atom);
  }

  for (const auto& element_to_be_changed : proto.change_element()) {
    std::unique_ptr<Reaction_Change_Element> rce(new Reaction_Change_Element);
    if (! rce->ConstructFromProto(element_to_be_changed))
       return WriteError("Reaction_Site::ConstructFromProto:invalid element to change", proto);
    _elements_to_change.add(rce.release());
  }

  for (const auto& formal_charge_to_assign : proto.formal_charge()) {
    std::unique_ptr<Reaction_Formal_Charge> rfc(new Reaction_Formal_Charge);
    if (! rfc->ConstructFromProto(formal_charge_to_assign))
       return WriteError("Reaction_Site::ConstructFromProto:invalid formal charge", proto);
    _formal_charges_to_assign.add(rfc.release());
  }

  for (const auto& formal_charge_to_be_changed : proto.change_formal_charge()) {
    std::unique_ptr<Reaction_Change_Formal_Charge> rcfc(new Reaction_Change_Formal_Charge);
    if (! rcfc->ConstructFromProto(formal_charge_to_be_changed))
       return WriteError("Reaction_Site::ConstructFromProto:invalid formal charge to change", proto);
    _formal_charges_to_change.add(rcfc.release());
  }

  for (const auto& isotope_to_assign : proto.isotope()) {
    std::unique_ptr<Reaction_Place_Isotope> rpi(new Reaction_Place_Isotope);
    if (! rpi->ConstructFromProto(isotope_to_assign))
       return WriteError("Reaction_Site::ConstructFromProto:invalid isotope", proto);
    _isotopes_to_assign.add(rpi.release());
  }

  if (proto.atom_and_isotope_size() > 0) {
    if (proto.atom_and_isotope_size() % 2 != 0) {
      return WriteError("Reaction_Site::ConstructFromProto:atom_and_isotope must have even number of entries",
                proto);
    }
    for (int i = 0; i < proto.atom_and_isotope_size(); i += 2) {
      uint32_t atom = proto.atom_and_isotope(i);
      uint32_t iso = proto.atom_and_isotope(i + 1);
      _isotopes_to_assign << new Reaction_Place_Isotope(atom, iso);
    }
  }

  for (const auto& isotope_to_be_changed : proto.change_isotope()) {
    std::unique_ptr<Reaction_Increment_Isotope> rii(new Reaction_Increment_Isotope);
    if (! rii->ConstructFromProto(isotope_to_be_changed))
       return WriteError("Reaction_Site::ConstructFromProto:invalid isotope to be changed", proto);
    _isotopes_to_increment.add(rii.release());
  }

  for (const auto& isotope_to_be_inverted : proto.invert_isotope()) {
    std::unique_ptr<Reaction_Invert_Isotope> rii(new Reaction_Invert_Isotope);
    if (! rii->ConstructFromProto(isotope_to_be_inverted))
       return WriteError("Reaction_Site::ConstructFromProto:invalid isotope to be inverted", proto);
    _isotopes_to_invert.add(rii.release());
  }

  // 3D related

  for (const auto& bond_length : proto.bond_length()) {
    std::unique_ptr<Reaction_Bond_Length> rbl(new Reaction_Bond_Length);
    if (! rbl->ConstructFromProto(bond_length, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid bond length", proto);
    _reaction_bond_length.add(rbl.release());
  }

  for (const auto& bond_angle : proto.bond_angle()) {
    std::unique_ptr<Reaction_Bond_Angle> rba(new Reaction_Bond_Angle);
    if (! rba->ConstructFromProto(bond_angle, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid bond angle", proto);
    _reaction_bond_angle.add(rba.release());
  }

  for (const auto& dihedral_angle : proto.dihedral_angle()) {
    std::unique_ptr<Reaction_Dihedral_Angle> rda(new Reaction_Dihedral_Angle);
    if (! rda->ConstructFromProto(dihedral_angle, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid dihedral_angle", proto);
    _reaction_dihedral_angle.add(rda.release());
  }

  for (const auto& coordinate_transfer : proto.coordinate_transfer()) {
    std::unique_ptr<Reaction_3D_Replace> r3d(new Reaction_3D_Replace);
    if (! r3d->ConstructFromProto(coordinate_transfer, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid reaction 3d replace", proto);
    _reaction_3d_replace.add(r3d.release());
  }

  for (const auto& wedge_bond : proto.wedge_bonds()) {
    std::unique_ptr<Reaction_Wedge_Bond> rwb(new Reaction_Wedge_Bond);
    if (! rwb->ConstructFromProto(wedge_bond))
       return WriteError("Reaction_Site::ConstructFromProto:invalid reaction wedge bond", proto);
    _wedge_bonds_to_place.add(rwb.release());
  }

  for (const auto& replace_atom : proto.replace_atom()) {
    std::unique_ptr<Replace_Atom> rpla(new Replace_Atom);
    if (! rpla->ConstructFromProto(replace_atom, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid replace atom", proto);
    _replace_atom.add(rpla.release());
  }

  for (const auto& data : proto.copy_isotope()) {
    std::unique_ptr<CopyIsotope> ci = std::make_unique<CopyIsotope>();
    if (! ci->ConstructFromProto(data, _unique_id - 1)) {
      cerr << "SidechainReactionSite::ConstructFromProto:invalid copy isotope\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
    _copy_isotope << ci.release();
  }

  for (const auto x : proto.unfix_implicit_hydrogens()) {
    _unfix_implicit_hydrogens.add_if_not_already_present(x);
  }

  for (const auto& inactive : proto.inactive()) {
    std::unique_ptr<Substructure_Query> ss(new Substructure_Query);
    if (! ss->ConstructFromProto(inactive))
       return WriteError("Reaction_Site::ConstructFromProto:invalid inactive query", proto);
    _inactive.add(ss.release());
  }

  for (const auto x : proto.invert_chirality())
    _stereo_centres_to_invert.add_if_not_already_present(x);

  for (const auto x : proto.remove_chirality())
    _chiral_centres_to_remove.add_if_not_already_present(x);

  for (const auto& cip : proto.cip_stereo()) {
    std::unique_ptr<SiteCipStereo> tmp = std::make_unique<SiteCipStereo>();
    if (! tmp->ConstructFromProto(cip)) {
      return WriteError("Reaction_Site::ConstructFromProto:invalid CIP stereo ", cip);
    }
    _cip_stereo << tmp.release();
  }

  if (proto.has_toggle_kekule_form())
  {
    if (! _toggle_kekule_form.ConstructFromProto(proto.toggle_kekule_form()))
      return WriteError("Reaction_Site::ConstructFromProto:invalid toggle kekule form", proto);
  }

  _ignore_multiple_matches_involving_atoms_not_changing = proto.ignore_multiple_matches_involving_atoms_not_changing();
  _ignore_multiple_matches_involving_changing_atoms = proto.ignore_multiple_matches_involving_changing_atoms();

  _noop_reaction = proto.noop_reaction();

  return 1;
}

int
Scaffold_Reaction_Site::ConstructFromProto(const ReactionProto::ScaffoldReactionSite& proto,
                                           const IWString& fname)
{
  if (proto.has_id()) {
    _unique_id = proto.id();
  } else {
    _unique_id = 0;
  }

  // Jan 2024.
  // The initial implementation of match_conditions with textproto input did not
  // properly process the match_conditions, so the default was never seen. Set it.
  if (! proto.has_match_conditions()) {
    set_find_unique_embeddings_only(0);
    _match_conditions.set_find_unique_embeddings_only(0);
  } else if (! _match_conditions.ConstructFromProto(proto.match_conditions())) {
    return WriteError("Reaction_Site::ConstructFromProto:invalid match conditions", proto);
  }

  if (!Reaction_Site::ConstructFromProto(proto, fname)) {
    return WriteError("ScaffoldReactionSite::ConstructFromProto:invalid Reaction_Site", proto);
  }

  // Transfer settings from _match_conditions to queries.
  InitialiseQueryConstraints<Scaffold_Match_Conditions>(_match_conditions);

  if (proto.single_bond_size() > 0) {
    if (proto.single_bond_size() % 2 != 0) {
      return WriteError("Reaction_Site::ConstructFromProto:single_bond requires even number of atoms ",
                        proto);
    }
    for (int i = 0; i < proto.single_bond_size(); i += 2) {
      uint32_t a1 = proto.single_bond(i);
      uint32_t a2 = proto.single_bond(i + 1);
      _bonds_to_be_made << new Bond(a1, a2, SINGLE_BOND);
    }
  }

  if (proto.double_bond_size() > 0) {
    if (proto.double_bond_size() % 2 != 0) {
      return WriteError("Reaction_Site::ConstructFromProto:double_bond requires even number of atoms ",
                        proto);
    }
    for (int i = 0; i < proto.double_bond_size(); i += 2) {
      uint32_t a1 = proto.double_bond(i);
      uint32_t a2 = proto.double_bond(i + 1);
      _bonds_to_be_made << new Bond(a1, a2, DOUBLE_BOND);
    }
  }

  if (proto.triple_bond_size() > 0) {
    if (proto.triple_bond_size() % 2 != 0) {
      return WriteError("Reaction_Site::ConstructFromProto:triple_bond requires even number of atoms ",
                        proto);
    }
    for (int i = 0; i < proto.triple_bond_size(); i += 2) {
      uint32_t a1 = proto.triple_bond(i);
      uint32_t a2 = proto.triple_bond(i + 1);
      _bonds_to_be_made << new Bond(a1, a2, TRIPLE_BOND);
    }
  }

  return 1;
}

// This is complicated by the fact that there may be match conditions specified in the
// proto and there might also be conditions in `smc` - which has typically come from the
// -M option to trxn. 
int
Sidechain_Reaction_Site::ConstructFromProto(const ReactionProto::SidechainReactionSite& proto,
                        int uid,
                        const IWString& fname,
                        const Sidechain_Match_Conditions& smc) {
  if (!proto.has_id()) {
    cerr << "SidechainReactionSite::ConstructFromProto:no id, setting " << uid << '\n';
    cerr << proto.ShortDebugString() << '\n';
    _unique_id = uid;
  } else if (uid == proto.id()) {
    _unique_id = proto.id();
  } else {
    cerr << "SidechainReactionSite::ConstructFromProto:non sequential id got " <<
             proto.id() << " expected " << uid << " expect trouble\n";
    cerr << proto.ShortDebugString() << '\n';
  }

  // First priority is what might have come from the command line.
  if (smc.active()) {
    _match_conditions = smc;
  } else if (proto.has_match_conditions()) {
    if (! _match_conditions.ConstructFromProto(proto.match_conditions())) {
      return WriteError("Reaction_Site::ConstructFromProto:invalid match conditions", proto);
    }
  } else {
    // In sidechains it is OK to suppress symmetry related matches by default.
    // Monitor this, symmetry has been an ongoing problem.
    //set_do_not_perceive_symmetry_equivalent_matches(1);
    _match_conditions.set_find_unique_embeddings_only(1);
  }

  if (!Reaction_Site::ConstructFromProto(proto, fname)) {
    return WriteError("ScaffoldReactionSite::ConstructFromProto:invalid Reaction_Site", proto);
  }

  // Transfer settings from _match_conditions to queries.
  InitialiseQueryConstraints<Sidechain_Match_Conditions>(_match_conditions);

  for (const auto& inter_particle_bond : proto.join()) {
    std::unique_ptr<Inter_Particle_Bond> ipb(new Inter_Particle_Bond);
    if (! ipb->ConstructFromProto(inter_particle_bond, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid inter particle bond", proto);
    _inter_particle_bonds.add(ipb.release());
  }

  if (proto.single_bond_size() > 0) {
    if (proto.single_bond_size() % 2 != 0) {
      cerr << "SidechainReactionSite::ConstructFromProto:number single bonds must be even\n";
      return 0;
    }
    for (int i = 0; i < proto.single_bond_size(); i += 2) {
      std::unique_ptr<Inter_Particle_Bond> ipb = std::make_unique<Inter_Particle_Bond>();
      ipb->SetDefaultJoin(proto.single_bond(i), _unique_id - 1, proto.single_bond(i + 1), SINGLE_BOND);
      _inter_particle_bonds << ipb.release();
    }
  }

  if (proto.double_bond_size() > 0) {
    if (proto.double_bond_size() % 2 != 0) {
      cerr << "SidechainReactionSite::ConstructFromProto:number double bonds must be even\n";
      return 0;
    }
    for (int i = 0; i < proto.double_bond_size(); i += 2) {
      std::unique_ptr<Inter_Particle_Bond> ipb = std::make_unique<Inter_Particle_Bond>();
      ipb->SetDefaultJoin(proto.double_bond(i), _unique_id - 1, proto.double_bond(i + 1), DOUBLE_BOND);
      _inter_particle_bonds << ipb.release();
    }
  }

  if (proto.triple_bond_size() > 0) {
    if (proto.triple_bond_size() % 2 != 0) {
      cerr << "SidechainReactionSite::ConstructFromProto:number triple bonds must be even\n";
      return 0;
    }
    for (int i = 0; i < proto.triple_bond_size(); i += 2) {
      std::unique_ptr<Inter_Particle_Bond> ipb = std::make_unique<Inter_Particle_Bond>();
      ipb->SetDefaultJoin(proto.triple_bond(i), _unique_id - 1, proto.triple_bond(i + 1), TRIPLE_BOND);
      _inter_particle_bonds << ipb.release();
    }
  }

  // TODO:ianwatson
  // Need to check to make sure that any _inter_particle_bonds bonds do NOT
  // involve atoms that are being removed. That is a silent bug right now.

  for (const auto& no_reaction : proto.no_reaction()) {
    std::unique_ptr<No_Reaction> nrxn(new No_Reaction);
    if (! nrxn->ConstructFromProto(no_reaction))
       return WriteError("Reaction_Site::ConstructFromProto:invalid no reaction", proto);
    _no_reaction.add(nrxn.release());
  }

  if (proto.has_make_implicit_hydrogens_explicit()) {
    _make_implicit_hydrogens_explicit = proto.make_implicit_hydrogens_explicit();
  }

  // Do this last after all query modifiers have been applied.
  if (proto.reagent_size() > 0) {
    if (_queries.size()) {
      _copy_match_conditions_to_query();
    } else {
      Substructure_Query::create_from_smarts("[*]");
      Substructure_Query::set_max_matches_to_find(1);
    }

    for (const auto& reagent : proto.reagent()) {
      if (! _add_reagent(reagent)) {
        cerr << "Sidechain_Reaction_Site::ConstructFromProto:invalid reagent " << reagent << "\n";
        return WriteError("Sidechain_Reaction_Site::ConstructFromProto:invalid reagent", proto);
      }
    }
  }


  return 1;
}

// Add an inter particle join involving matched atom `a1` in the scaffold and matched
// atoms `a2` in our sidechain.
int
Sidechain_Reaction_Site::AddInterParticleBond(atom_number_t a1, atom_number_t a2,
                bond_type_t btype) {
  static constexpr int kScaffold = -1;

  _inter_particle_bonds << new Inter_Particle_Bond(kScaffold, a1, a2, btype);

  // A little dangerous, we should be using the index into the _sidechain array.
  _inter_particle_bonds.back()->is_part_of_component(_unique_id - 1);

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagent(const std::string& smiles)
{
  const_IWSubstring mysmiles(smiles.data(), smiles.length());
  cerr << "Addomg reagent " << mysmiles << '\n';

  Molecule_and_Embedding * mae = new Molecule_and_Embedding();
  if (! mae->build_from_smiles(mysmiles)) {
    cerr << "Sidechain_Reaction_Site::_add_reagent:invalid reagent smiles '" << smiles << "'\n";
    return 0;
  }

  if (_make_implicit_hydrogens_explicit) {
    Make_Implicit_Hydrogens_Explicit mihe;
    mihe.set_isotope(_make_implicit_hydrogens_explicit);
    mae->make_implicit_hydrogens_explicit(mihe);
  }

  if (! add_reagent(mae, _match_conditions)) {
    delete mae;
    cerr << "Sidechain_Reaction_Site::_add_reagent:cannot add reagent " << smiles << "\n";
    return 0;
  }

  return 1;
}

int
IWReaction::Read(IWString& fname) {
  std::optional<ReactionProto::Reaction> maybe_proto =
      iwmisc::ReadTextProtoCommentsOK<ReactionProto::Reaction>(fname);
  if (! maybe_proto) {
    cerr << "IWReaction::Read:cannot read reaction proto from '" << fname << "'\n";
    return 0;
  }

  Sidechain_Match_Conditions smc;

  return ConstructFromProto(*maybe_proto, fname, smc);
}

int
IWReaction::ConstructFromProto(const ReactionProto::Reaction& proto,
                               const IWString& file_name,
                               const Sidechain_Match_Conditions& smc) {
  if (! InternallyConsistent(proto)) {
    return 0;
  }

  if (!proto.has_scaffold()) {
    return WriteError("IWReaction::ConstructFromProto:no scaffold", proto);
  }

  if (! Scaffold_Reaction_Site::ConstructFromProto(proto.scaffold(), file_name)) {
    return WriteError("IWReaction::ConstructFromProto:invalid scaffold", proto);
  }

  if (proto.has_name()) {
    _comment = proto.name();
  } else  {
    MaybeCopyComment(proto, _comment);
  }

  if (proto.has_scaffold_match_conditions()) {
    cerr << "IWReaction::ConstructFromProto:deprecation error, reaction has scaffold_match_conditions\n";
    cerr << "Move to ScaffoldMatchConditions message in the scaffold\n";
    return 0;
  }
#ifdef OBSOLETE
  if (proto.has_scaffold_match_conditions() &&
      ! _match_conditions.ConstructFromProto(proto.scaffold_match_conditions()))
    return WriteError("IWReaction::ConstructFromProto:invalid scaffold match conditions", proto);

  // Copy match conditions to query.
  if (_match_conditions.find_unique_embeddings_only()) {
    set_find_unique_embeddings_only(1);
  }
  if (_match_conditions.one_embedding_per_start_atom()) {
    set_one_embedding_per_start_atom(1);
  }
  if (! _match_conditions.embeddings_can_overlap()) {
    set_embeddings_do_not_overlap(1);
  }
#endif

  if (proto.sidechain_size() > 0) {
    for (int i = 0; i < proto.sidechain_size(); ++i) {
      std::unique_ptr<Sidechain_Reaction_Site> sc(new Sidechain_Reaction_Site);
      if (! sc->ConstructFromProto(proto.sidechain(i), i+1, file_name, smc)) {
        return WriteError("IWReaction::ConstructFromProto:invalid sidechain", proto);
      }
      _sidechains.add(sc.release());
    }
  }

//cerr << "Read " << _sidechains.number_elements() << " sidechains\n";

  for (const auto& stereo_center : proto.reaction_stereo_center()) {
    std::unique_ptr<Reaction_Stereo_Centre> sc(new Reaction_Stereo_Centre);
    if (! sc->ConstructFromProto(stereo_center)) {
      return WriteError("IWReaction::ConstructFromProto:invalid reacton stereo center", proto);
    }
    _reaction_stereo_centre.add(sc.release());
  }

  _make_implicit_hydrogens_explicit = proto.make_implicit_hydrogens_explicit();
  if (_make_implicit_hydrogens_explicit) {
    for (Sidechain_Reaction_Site* r : _sidechains) {
      r->set_make_implicit_hydrogens_explicit(true);
    }
  }

  if (proto.has_append_reagent_name()) {
    _append_names = proto.append_reagent_name();
  }

  if (proto.has_append_to_name()) {
    _append_to_name = proto.append_to_name();
  }

  _query_files_in_current_directory = proto.query_files_in_current_directory();

  _find_kekule_forms_for_bad_valence = proto.find_kekule_forms_for_bad_valence();

  if (proto.has_noop_reaction()) {
    _noop_reaction = proto.noop_reaction();
  }

  for (const auto& cip:  proto.cip_stereo()) {
    std::unique_ptr<ReactionCipStereo> tmp = std::make_unique<ReactionCipStereo>();
    if (! tmp->ConstructFromProto(cip)) {
      return WriteError("IWReaction::ConstructFromProto:invalid cip stereo center", proto);
    }
    _cip_stereo << tmp.release();
  }

  _determine_has_sidechain_isotope_requirement();

  if (! check_internal_consistency()) {
    return WriteError("IWReaction::ConstructFromProto:not internally consistent", proto);
  }

  return 1;
}

int
IWReaction::_determine_has_sidechain_isotope_requirement() {
  _has_sidechain_isotope_requirement = 0;

  for (const Sidechain_Reaction_Site* s : _sidechains) {
    if (s->has_sidechain_isotope_requirement()) {
      ++_has_sidechain_isotope_requirement;
    }
  }

  return _has_sidechain_isotope_requirement;
}


int
Sidechain_Reaction_Site::has_sidechain_isotope_requirement() const {
  for (Inter_Particle_Bond* b : _inter_particle_bonds) {
    if (b->has_sidechain_isotope_requirement()) {
      return 1;
    }
  }

  return 0;
}

int
SiteCipStereo::ConstructFromProto(ReactionProto::CipStereoAtom const& proto) {
  if (! proto.has_atom()) {
    cerr << "SiteCipStereo::ConstructFromProto:no atom " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_rs()) {
    cerr << "SiteCipStereo::ConstructFromProto:no rs " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _atom = proto.atom();
  switch (proto.rs()) {
    case SubstructureSearch::CIP_R:
      _rs = CahnIngoldPrelog::R;
      break;
    case SubstructureSearch::CIP_S:
      _rs = CahnIngoldPrelog::S;
      break;
    default:
      break;
  }

  return 1;
}

int
SiteCipStereo::BuildProto(ReactionProto::CipStereoAtom& proto) const {
  proto.set_atom(_atom);
  switch (_rs) {
    case CahnIngoldPrelog::R:
      proto.set_rs(SubstructureSearch::CIP_R);
      break;
    case CahnIngoldPrelog::S:
      proto.set_rs(SubstructureSearch::CIP_S);
      break;
    default:
      cerr << "SiteCipStereo::BuildProto:what stereo " << _rs << '\n';
  }

  return 1;
}

int
ReactionCipStereo::ConstructFromProto(ReactionProto::CipStereoReaction const& proto) {
  if (! proto.has_atom()) {
    cerr << "ReactionCipStereo::ConstructFromProto:no atom " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_rs()) {
    cerr << "ReactionCipStereo::ConstructFromProto:no rs " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _atom.ConstructFromProto(proto.atom())) {
    cerr << "ReactionCipStereo::ConstructFromProto:cannot build mached atom " << proto.ShortDebugString() << '\n';
    return 0;
  }

  switch (proto.rs()) {
    case SubstructureSearch::CIP_R:
      _rs = CahnIngoldPrelog::R;
      break;
    case SubstructureSearch::CIP_S:
      _rs = CahnIngoldPrelog::S;
      break;
    default:
      break;
  }

  return 1;
}

int
IWReaction::ConstructFromTextProto(const std::string& textproto, const IWString& file_name,
                const Sidechain_Match_Conditions& smc) {
  ReactionProto::Reaction proto;
  if (!google::protobuf::TextFormat::ParseFromString(textproto, &proto)) {
    cerr << "IWReaction:ConstructFromTextProto:cannot parse text proto " << textproto << '\n';
    return 0;
  }

  return ConstructFromProto(proto, file_name, smc);
}

int
IWReaction::BuildProto(ReactionProto::Reaction& proto) const {
  if (! _comment.empty()) {
    proto.set_name(_comment.data(), _comment.size());
  }

  proto.mutable_scaffold()->set_id(0);

  // TODO:ianwatson implement these sometime
  for (const Reaction_Stereo_Centre* rsc : _reaction_stereo_centre) {
    (void) rsc;
  }
    
  for (const ReactionCipStereo* rcs : _cip_stereo) {
    (void) rcs;
  }

  if (_append_names) {
    proto.set_append_reagent_name(true);
  }
  if (! _append_to_name.empty()) {
    proto.set_append_to_name(_append_to_name.data(), _append_to_name.size());
  }
  if (_query_files_in_current_directory) {
    proto.set_query_files_in_current_directory(true);
  }
  if (! _reaction_directory.empty()) {
    proto.set_reaction_directory(_reaction_directory.data(), _reaction_directory.size());
  }
  if (_find_kekule_forms_for_bad_valence) {
    proto.set_find_kekule_forms_for_bad_valence(true);
  }
  if (_make_implicit_hydrogens_explicit) {
    proto.set_make_implicit_hydrogens_explicit(true);
  }

  if (_smarts.size()) {
    for (const IWString* s : _smarts) {
      std::string tmp(s->data(), s->size());
      proto.mutable_scaffold()->add_smarts(tmp);
    }
  }

  cerr << "Reaction has " << _sidechains.size() << " sidechains\n";

  int id = 1;
  for (const Sidechain_Reaction_Site* sidechain : _sidechains) {
    auto* s = proto.add_sidechain();
    s->set_id(id);
    ++id;
    sidechain->BuildProto(*s);
  }

  Reaction_Site::BuildProto(*proto.mutable_scaffold());

  // This should be deprecated, match conditions are no longer part of a reaction
  // but instead belong to scaffolds and sidechains.
  if (! _match_conditions.IsDefault()) {
    AddMatchConditions(*proto.mutable_scaffold());
  }

  return 1;
}

// BOth the Scaffold_Reaction_Site and Sidechain_Reaction_Site need to set
// match conditions kinds of proto messages. 
// The underlying Match_Conditions object can set the common values into
// whichever type P is.
template <typename P>
int
Match_Conditions::BuildProto(P& proto) const {

  if (_ignore_not_reacting) {
    proto.set_ignore_not_reacting(_ignore_not_reacting);
  }

  if (_find_unique_embeddings) {
    proto.set_find_unique_embeddings(_find_unique_embeddings);
  }

  if (_process_hit_number >= 0) {
    proto.set_process_hit_number(_process_hit_number);
  }

  if (_one_embedding_per_start_atom) {
    proto.set_one_embedding_per_start_atom(_one_embedding_per_start_atom);
  }

  if (_ignore_symmetry_related_matches) {
    proto.set_ignore_symmetry_related_matches(_ignore_symmetry_related_matches);
  }

  if (_embeddings_can_overlap) {
    proto.set_embeddings_can_overlap(true);
  } else {
    proto.set_embeddings_can_overlap(false);
  }

  if (_multiple_match_string.size() > 0) {
    const IWString& s = _multiple_match_string;
    proto.set_multiple_match_string(s.data(), s.length());
  }

  if (uint32_t s = _suppress_if_more_than_this_many_substructure_search_hits;
        s < std::numeric_limits<int32_t>::max()) {
    proto.set_suppress_if_more_than_this_many_substructure_search_hits(s);
  }

  if (_max_matches_to_find > 0) {
    proto.set_max_matches_to_find(_max_matches_to_find);
  }

  return 1;
}

// Call the underlying Reaction_Site method for all the inherited attributes.
int
Scaffold_Reaction_Site::AddMatchConditions(ReactionProto::ScaffoldReactionSite& proto) const {
  if (_match_conditions.IsDefault()) {
    return 1;
  }

  ReactionProto::ScaffoldMatchConditions* smc = proto.mutable_match_conditions();

  _match_conditions.BuildProto(*smc);

  return 1;
}

template <typename P>
int
Reaction_Site::BuildProto(P& proto) const {
  for (const int a : _atoms_to_be_removed) {
    proto.add_remove_atom(a);
  }

  for (const Bond* bond : _bonds_to_be_made) {
    auto* b = proto.add_make_bond();
    b->set_a1(bond->a1());
    b->set_a2(bond->a2());
    if (bond->is_single_bond()) {
      b->set_btype(SubstructureSearch::SS_SINGLE_BOND);
    } else if (bond->is_double_bond()) {
      b->set_btype(SubstructureSearch::SS_DOUBLE_BOND);
    } else if (bond->is_triple_bond()) {
      b->set_btype(SubstructureSearch::SS_TRIPLE_BOND);
    }
  }

  for (const Bond* bond : _bonds_to_be_changed) {
    auto* b = proto.add_change_bond();
    b->set_a1(bond->a1());
    b->set_a2(bond->a2());
    if (bond->is_single_bond()) {
      b->set_btype(SubstructureSearch::SS_SINGLE_BOND);
    } else if (bond->is_double_bond()) {
      b->set_btype(SubstructureSearch::SS_DOUBLE_BOND);
    } else if (bond->is_triple_bond()) {
      b->set_btype(SubstructureSearch::SS_TRIPLE_BOND);
    }
  }

  for (const Bond* bond : _bonds_to_be_broken) {
    auto* b = proto.add_break_bond();
    b->set_a1(bond->a1());
    b->set_a2(bond->a2());
  }

  for (int f : _fragments_to_be_removed) {
    proto.add_remove_fragment(f);
  }

  for (int f : _fragments_to_be_kept) {
    proto.add_keep_fragment(f);
  }

  for (const Reaction_Change_Element* ce :  _elements_to_change) {
    auto* q = proto.add_change_element();
    q->set_atom(ce->atom());
    const IWString& s = ce->element()->symbol();
    q->set_element(s.data(), s.size());
  }

  for (const Reaction_Formal_Charge* rfc:  _formal_charges_to_assign) {
    auto* a = proto.add_formal_charge();
    a->set_atom(rfc->atom());
    a->set_formal_charge(rfc->charge());
  }
    
  for (const Reaction_Change_Formal_Charge* rcfc : _formal_charges_to_change) {
    auto* q = proto.add_change_formal_charge();
    q->set_atom(rcfc->atom());
    q->set_delta(rcfc->delta());
  }

  for (const Reaction_Place_Isotope* rpi :  _isotopes_to_assign) {
    auto* q = proto.add_isotope();
    rpi->BuildProto(*q);
  }
    
  for (const Reaction_Increment_Isotope* rii : _isotopes_to_increment) {
    auto* q = proto.add_change_isotope();
    rii->BuildProto(*q);
  }
    
  for (const Reaction_Invert_Isotope* rii :  _isotopes_to_invert) {
    auto* q = proto.add_invert_isotope();
    rii->BuildProto(*q);
  }

  // TODO:ianwatson implement these....
  for (const Reaction_Dihedral_Angle* rda : _reaction_dihedral_angle) {
    auto * q = proto.add_dihedral_angle();
    rda->BuildProto(*q);
  }

  for (const Reaction_Bond_Length* rbl : _reaction_bond_length) {
    auto* q  = proto.add_bond_length();
    rbl->BuildProto(*q);
  }

  for (const Reaction_Bond_Angle* rba : _reaction_bond_angle) {
    auto* q = proto.add_bond_angle();
    rba->BuildProto(*q);
  }

  for (const Reaction_3D_Replace* r3dr : _reaction_3d_replace) {
    auto* q = proto.add_coordinate_transfer();
    r3dr->BuildProto(*q);
  }

  for (const Replace_Atom* ra : _replace_atom) {
    auto* q = proto.add_replace_atom();
    ra->BuildProto(*q);
  }
    
  if (_inactive.size() > 0) {
    cerr << "Reaction_Site::BuildProto:_inactive queries not handled\n";
  }
  for (const Substructure_Query* q : _inactive) {
    (void) q;
  }
    
  for (int a : _stereo_centres_to_invert) {
    proto.add_invert_chirality(a);
  }

  for (int a : _chiral_centres_to_remove) {
    proto.add_remove_chirality(a);
  }

  if (_toggle_kekule_form.active()) {
    _toggle_kekule_form.BuildProto(*proto.mutable_toggle_kekule_form());
  }

  if (_ignore_multiple_matches_involving_atoms_not_changing) {
    proto.set_ignore_multiple_matches_involving_atoms_not_changing(true);
  }

  if (_ignore_multiple_matches_involving_changing_atoms) {
    proto.set_ignore_multiple_matches_involving_changing_atoms(true);
  }

  if (_noop_reaction) {
    proto.set_noop_reaction(true);
  }

  for (const SiteCipStereo* scp : _cip_stereo) {
    auto* q = proto.add_cip_stereo();
    scp->BuildProto(*q);
  }

  return 1;
}

int
Sidechain_Reaction_Site::BuildProto(ReactionProto::SidechainReactionSite& proto) const {
  for (const IWString* smarts : _smarts) {
    const std::string tmp(smarts->data(), smarts->size());
    proto.add_smarts(tmp);
  }

  for (const Inter_Particle_Bond* bond : _inter_particle_bonds) {
    bond->BuildProto(*proto.add_join());
  }

  Reaction_Site::BuildProto(proto);

  for (Molecule_and_Embedding* r : _reagents) {
    const IWString& s = r->smiles();
    if (r->name().empty()) {
      proto.add_reagent(s.data(), s.size());
    } else {
      IWString tmp(s);
      tmp << ' ' << r->name();
      proto.add_reagent(tmp.data(), tmp.size());
    }
  }

  Reaction_Site::BuildProto(proto);

  if (! _match_conditions.IsDefault()) {
    _match_conditions.BuildProto(*proto.mutable_match_conditions());
  }

  return 1;
}

int
Sidechain_Match_Conditions::BuildProto(ReactionProto::SidechainMatchConditions& proto) const {
  Match_Conditions::BuildProto(proto);

  if (_make_new_reagent_for_each_hit) {
    proto.set_make_new_reagent_for_each_hit(_make_new_reagent_for_each_hit);
  }

  if (_strip_reagents_to_largest_fragment) {
    proto.set_strip_reagents_to_largest_fragment(true);
  }

  if (_make_new_reagent_for_each_hit) {
    proto.set_make_new_reagent_for_each_hit(true);
  }

  return 1;
}

int
Scaffold_Match_Conditions::BuildProto(ReactionProto::ScaffoldMatchConditions& proto) const {
  Match_Conditions::BuildProto(proto);

  if (_combinatorial_expansion_of_scaffold_hits) {
    proto.set_combinatorial_expansion_of_scaffold_hits(_combinatorial_expansion_of_scaffold_hits);
  }

  if (_enumerate_scaffold_hits_individually) {
    proto.set_enumerate_scaffold_hits_individually(true);
  }

  return 1;
}

int
Inter_Particle_Bond::BuildProto(ReactionProto::InterParticleBond& proto) const {
  if (_bt == SINGLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_SINGLE_BOND);
  } else if (_bt == DOUBLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_DOUBLE_BOND);
  } else if (_bt == TRIPLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_TRIPLE_BOND);
  } else {
    cerr << "InterParticleBond::BuildProto:unrecognised btype " << _bt << '\n';
    // should we return?
  }

  if (_a1.in_scaffold()) {
    proto.set_a1(_a1.atom());
  } else {
    auto* c = proto.mutable_c1();
    c->set_component(_a1.in_component() + 1);
    c->set_atom(_a1.atom());
  }

  proto.set_a2(_a2.atom());

  return 1;
}

int
Reaction_Place_Isotope::BuildProto(ReactionProto::PlaceIsotope& proto) const {
  for (int a : _atom) {
    proto.add_atom(a);
  }

  proto.set_isotope(_isotope);

  return 1;
}

int
Reaction_Increment_Isotope::BuildProto(ReactionProto::IncrementIsotope& proto) const {
  for (int a : _atom) {
    proto.add_atom(a);
  }

  if (_multiply > 0) {
    proto.set_multiply(_multiply);
  }
  if (_isotope != 0) {
    proto.set_delta(_isotope);
  }
  if (_integer_divide > 0) {
    proto.set_integer_divide(_integer_divide);
  }
  if (_modulo > 0) {
    proto.set_modulo(_modulo);
  }

  if (! _only_change_if.is_set())  {
  } else if (_only_change_if.empty()) {  // must be min and max
    isotope_t tmp = 0;
    _only_change_if.min(tmp);
    proto.mutable_only_change_if()->mutable_range()->set_min(tmp);
    _only_change_if.max(tmp);
    proto.mutable_only_change_if()->mutable_range()->set_max(tmp);
  } else {
    for (int a : _only_change_if.ValuesMatched()) {
      proto.mutable_only_change_if()->mutable_values()->mutable_values()->Add(a);
    }
  }

  return 1;
}

int
Reaction_Invert_Isotope::BuildProto(ReactionProto::PlaceIsotope& proto) const {
  for (int a : _atom) {
    proto.add_atom(a);
  }

  return 1;
}

int
Replace_Atom::BuildProto(ReactionProto::ReplaceAtom& proto) const {
  if (_a1.in_scaffold()) {
    proto.set_a1(_a1.atom());
  } else {
    auto* c = proto.mutable_c1();
    c->set_component(_a1.in_component());
    c->set_atom(_a1.atom());
  }

  proto.set_a2(_a2.atom());

  if (_invert_chirality) {
    proto.set_invert_chirality(true);
  }

  return 1;
}

int
CopyIsotope::ConstructFromProto(const ReactionProto::CopyIsotopeData& proto,
                                int component) {
  if (! proto.has_from() && ! proto.has_cfrom())
    return WriteError("CopyIsotope::ConstructFromProto:no from1/cfrom", proto);
  if (! proto.has_to() && ! proto.has_cto())
    return WriteError("CopyIsotope::ConstructFromProto:no from2/cto", proto);

  // cerr << "CopyIsotope::ConstructFromProto component " << component << '\n';
  if (proto.has_from()) {
    _from.set_matched_atom(proto.from());
    _from.set_in_component(component);
  } else if (proto.has_cfrom()) {
    if (! _from.ConstructFromProto(proto.cfrom()))
      return WriteError("Inter_Particle_Bond::ConstructFromProto:invalid cfrom", proto);
  }

  if (proto.has_to()) {
    _to.set_matched_atom(proto.to());
    _to.set_in_component(component);
  } else if (proto.has_cto()) {
    if (! _to.ConstructFromProto(proto.cto()))
      return WriteError("CopyIsotope::ConstructFromProto:invalid cto", proto);
  }

  if (_from.in_component() == _to.in_component() && _from.atom() == _to.atom()) {
    cerr << "CopyIsotope::ConstructFromProto:from and to identical\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
CopyIsotope::BuildProto(ReactionProto::CopyIsotopeData& proto) const {
  cerr << "CopyIsotopeData::BuildProto:not implemented, see Ian\n";
  return 1;
}
