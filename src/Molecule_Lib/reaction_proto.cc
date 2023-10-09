#include <filesystem>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "iwreaction.h"

using std::cerr;
using std::endl;

namespace fs = std::filesystem;

template <typename P>
int
WriteError(const char * message, const P & proto)
{
  cerr << message << ' ' << proto.ShortDebugString() << "\n";

  return 0;
}

int
Match_Conditions::ConstructFromProto(const ReactionProto::MatchConditions& proto)
{
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

  return 1;
}

int
Scaffold_Match_Conditions::ConstructFromProto(const ReactionProto::ScaffoldMatchConditions& proto)
{
  if (! proto.has_match_conditions())
    return WriteError("Scaffold_Match_Conditions::ConstructFromProto:no match conditions", proto);

  if (!Match_Conditions::ConstructFromProto(proto.match_conditions()))
    return WriteError("Scaffold_Match_Conditions::ConstructFromProto:invalid match conditions", proto);

  _enumerate_scaffold_hits_individually  = proto.enumerate_scaffold_hits_individually();
  _combinatorial_expansion_of_scaffold_hits = proto.combinatorial_expansion_of_scaffold_hits();

  return 1;
}

int
Sidechain_Match_Conditions::ConstructFromProto(const ReactionProto::SidechainMatchConditions& proto)
{
  if (! proto.has_match_conditions())
    return WriteError("Sidechain_Match_Conditions::ConstructFromProto:no match conditions", proto);

  if (!Match_Conditions::ConstructFromProto(proto.match_conditions()))
    return WriteError("Sidechain_Match_Conditions::ConstructFromProto:invalid match conditions", proto);

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
  if (! proto.has_atom())
    return WriteError("Reaction_Place_Isotope::ConstructFromProto:no atom", proto);
  if (! proto.has_isotope())
    return WriteError("Reaction_Place_Isotope::ConstructFromProto:no isotope", proto);

  _atom = proto.atom();
  _isotope = proto.isotope();

  return 1;
}

int
Reaction_Increment_Isotope::ConstructFromProto(const ReactionProto::IncrementIsotope& proto)
{
  if (! proto.has_atom())
    return WriteError("Reaction_Increment_Isotope::ConstructFromProto:no atom", proto);
  if (! proto.has_delta())
    return WriteError("Reaction_Increment_Isotope::ConstructFromProto:no delta", proto);

  _atom = proto.atom();
  _isotope = proto.delta();

  return 1;
}

int
Reaction_Invert_Isotope::ConstructFromProto(const ReactionProto::PlaceIsotope& proto)
{
  if (! proto.has_atom())
    return WriteError("Reaction_Invert_Isotope::ConstructFromProto:no atom", proto);
  if (! proto.has_isotope())
    return WriteError("Reaction_Invert_Isotope::ConstructFromProto:no isotope", proto);

  _atom = proto.atom();
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

#ifdef OBSOLETE_AASD
int
Reaction_3D_Replace::ConstructFromProto(const ReactionProto::ThreeDReplace& proto)
{
  if (proto.a1().empty())
    return WriteError("Reaction_3D_Replace::ConstructFromProto:no a1", proto);
  if (proto.a2().empty())
    return WriteError("Reaction_3D_Replace::ConstructFromProto:no a2", proto);

  if (proto.a1_size() != proto.a2_size())
    return WriteError("Reaction_3D_Replace::ConstructFromProto:inconsistent a1/a2 sizes", proto);

  _n = proto.a1_size();
  _a1 = new Matched_Atom_in_Component[_n];
  _a2 = new Matched_Atom_in_Component[_n];

  for (int i = 0; i < _n; ++i)
  {
    if (! _a1[i].ConstructFromProto(proto.a1(i)) ||
        ! _a2[i].ConstructFromProto(proto.a2(i)))
      return WriteError("Reaction_3D_Replace::ConstructFromProto:invalid a1/a2", proto);
  }

  _weight = new double[_n];

  _weight[0] = 1.0;

  for (int i = 1; i < _n; i++)
  {
    _weight[i] = 0.1;
  }

  return 1;
}
#endif

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
    _a1.set_in_component(component - 1);
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


  return 1;
}

int
No_Reaction::ConstructFromProto(const ReactionProto::NoReaction& proto)
{
  if (proto.has_comment())
    _comment = proto.comment();

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

template <typename P>
int
Reaction_Site::ConstructFromProto(const P& proto, const IWString& fname)
{
  if (proto.has_comment())
    _comment = proto.comment();

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

  for (const auto& bond_to_be_broken : proto.break_bond()) {
    std::unique_ptr<Bond> b(new Bond());
    b->set_bond_type(SINGLE_BOND);   // Seems invalid if not set.
    if (! BondFromPairOfAtoms(*b, bond_to_be_broken, /*no self bonds*/false))
      return WriteError("Reaction_Site::ConstructFromProto:invalid bond to be made ", bond_to_be_broken);
    _bonds_to_be_broken.add(b.release());
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

  for (const auto x : proto.unfix_implicit_hydrogens())
    _unfix_implicit_hydrogens.add_if_not_already_present(x);

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

  if (!Reaction_Site::ConstructFromProto(proto, fname)) {
    return WriteError("ScaffoldReactionSite::ConstructFromProto:invalid Reaction_Site", proto);
  }

  return 1;
}

int
Sidechain_Reaction_Site::ConstructFromProto(const ReactionProto::SidechainReactionSite& proto,
                        const IWString& fname)
{
  if (!proto.has_id())
    return WriteError("Sidechain_Reaction_Site::ConstructFromProto:no id", proto);

  _unique_id = proto.id();

  if (!Reaction_Site::ConstructFromProto(proto, fname))
      return WriteError("ScaffoldReactionSite::ConstructFromProto:invalid Reaction_Site", proto);

  for (const auto& inter_particle_bond : proto.join()) {
    std::unique_ptr<Inter_Particle_Bond> ipb(new Inter_Particle_Bond);
    if (! ipb->ConstructFromProto(inter_particle_bond, _unique_id))
       return WriteError("Reaction_Site::ConstructFromProto:invalid inter particle bond", proto);
    _inter_particle_bonds.add(ipb.release());
  }

  // TODO:ianwatson
  // Need to check to make sure that any _inter_particle_bonds bonds do NOT
  // involve atoms that are being removed. That is a silent bug right now.

  if (proto.has_match_conditions() &&
      ! _match_conditions.ConstructFromProto(proto.match_conditions()))
    return WriteError("Reaction_Site::ConstructFromProto:invalid match conditions", proto);

  for (const auto& no_reaction : proto.no_reaction()) {
    std::unique_ptr<No_Reaction> nrxn(new No_Reaction);
    if (! nrxn->ConstructFromProto(no_reaction))
       return WriteError("Reaction_Site::ConstructFromProto:invalid no reaction", proto);
    _no_reaction.add(nrxn.release());
  }

  if (proto.has_make_implicit_hydrogens_explicit())
    _make_implicit_hydrogens_explicit = proto.make_implicit_hydrogens_explicit();

  // Do this last after all query modifiers have been applied.
  if (proto.reagent_size() > 0)
  {
    set_do_not_perceive_symmetry_equivalent_matches(1);   // let's just make this the default

    _copy_match_conditions_to_query();

    for (const auto& reagent : proto.reagent()) {
      if (! _add_reagent(reagent))
      {
        cerr << "Sidechain_Reaction_Site::ConstructFromProto:invalid reagent " << reagent << "\n";
        return WriteError("Sidechain_Reaction_Site::ConstructFromProto:invalid reagent", proto);
      }
    }
  }

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagent(const std::string& smiles)
{
  const_IWSubstring mysmiles(smiles.data(), smiles.length());

  Molecule_and_Embedding * mae = new Molecule_and_Embedding();
  if (! mae->build_from_smiles(mysmiles))
  {
    cerr << "Sidechain_Reaction_Site::_add_reagent:invalid reagent smiles '" << smiles << "'\n";
    return 0;
  }

  if (_make_implicit_hydrogens_explicit)
  {
    Make_Implicit_Hydrogens_Explicit mihe;
    mihe.set_isotope(_make_implicit_hydrogens_explicit);
    mae->make_implicit_hydrogens_explicit(mihe);
  }

  if (! add_reagent(mae, _match_conditions))
  {
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

  return ConstructFromProto(*maybe_proto, fname);
}

int
IWReaction::ConstructFromProto(const ReactionProto::Reaction& proto,
                               const IWString& file_name)
{
  if (!proto.has_scaffold()) {
    return WriteError("IWReaction::ConstructFromProto:no scaffold", proto);
  }

  if (! Scaffold_Reaction_Site::ConstructFromProto(proto.scaffold(), file_name)) {
    return WriteError("IWReaction::ConstructFromProto:invalid scaffold", proto);
  }

  if (proto.has_comment()) {
    _comment = proto.comment();
  } else if (proto.has_name()) {
    _comment = proto.name();
  }

  if (proto.has_scaffold_match_conditions() &&
      ! _match_conditions.ConstructFromProto(proto.scaffold_match_conditions()))
    return WriteError("IWReaction::ConstructFromProto:invalid scaffold match conditions", proto);

  for (const auto& sidechain : proto.sidechain()) {
    std::unique_ptr<Sidechain_Reaction_Site> sc(new Sidechain_Reaction_Site);
    if (! sc->ConstructFromProto(sidechain, file_name)) {
      return WriteError("Reaction_Site::ConstructFromProto:invalid sidechain", proto);
    }
    _sidechains.add(sc.release());
  }

//cerr << "Read " << _sidechains.number_elements() << " sidechains\n";

  for (const auto& stereo_center : proto.reaction_stereo_center()) {
    std::unique_ptr<Reaction_Stereo_Centre> sc(new Reaction_Stereo_Centre);
    if (! sc->ConstructFromProto(stereo_center)) {
      return WriteError("Reaction_Site::ConstructFromProto:invalid reacton stereo center", proto);
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
      return WriteError("Reaction_Site::ConstructFromProto:invalid cip stereo center", proto);
    }
    _cip_stereo << tmp.release();
  }

  _determine_has_sidechain_isotope_requirement();

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
