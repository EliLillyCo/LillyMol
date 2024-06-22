#include <cstdint>

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>

/*
  Functions associated with reading and writing substructure query objects and
  such.
*/

#include <google/protobuf/message.h>
#include <google/protobuf/text_format.h>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/parse_smarts_tmp.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"
#include "Molecule_Lib/target.h"

using std::cerr;

using down_the_bond::DownTheBond;

constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

// Having unbalanced braces in the code messes up matching in the editor.
constexpr char kOpenBrace = '{';
constexpr char kCloseBrace = '}';

using google::protobuf::Descriptor;
using google::protobuf::FieldDescriptor;
using google::protobuf::Reflection;

template <typename P>
void
SetValue(uint32_t value,
         const char * min_max,
         const char * name_stem,
         const Descriptor* descriptor,
         const Reflection* reflection,
         P& proto) {
  IWString name;
  name << min_max << name_stem;
  const FieldDescriptor* fd = descriptor->FindFieldByName(name.null_terminated_chars());
  reflection->SetUInt32(&proto, fd, value);
}

template <typename P>
void
SetValueInt(int32_t value,
         const char * min_max,
         const char * name_stem,
         const Descriptor* descriptor,
         const Reflection* reflection,
         P& proto) {
  IWString name;
  name << min_max << name_stem;
  const FieldDescriptor* fd = descriptor->FindFieldByName(name.null_terminated_chars());
  reflection->SetInt32(&proto, fd, value);
}

// Convert `values` to `proto`.
// The name for the proto field is `name_stem`, with
// min_name_stem and max_name_stem possibly being populated.
template <typename TheMatcher, typename P>
int
SetProtoValues(const TheMatcher& values,
               const char * name_stem,
               P& proto) {
  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  if (values.size() > 0) {
    const FieldDescriptor * fd = descriptor->FindFieldByName(name_stem);
    for (int v : values) {
      reflection->AddUInt32(&proto, fd, static_cast<uint32_t>(v));
    }
  }

  int v;
  if (values.min(v)) {
    SetValue<P>(v, "min_", name_stem, descriptor, reflection, proto);
  }
  if (values.max(v)) {
    SetValue<P>(v, "max_", name_stem, descriptor, reflection, proto);
  }

  return 1;
}
template <typename TheMatcher, typename P>
int
SetProtoValuesInt(const TheMatcher & values,
               const char * name_stem,
               P& proto) {
  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  if (values.size() > 0) {
    const FieldDescriptor * fd = descriptor->FindFieldByName(name_stem);
    for (int v : values) {
      reflection->AddInt32(&proto, fd, static_cast<uint32_t>(v));
    }
  }

  int v;
  if (values.min(v)) {
    SetValueInt<P>(v, "min_", name_stem, descriptor, reflection, proto);
  }
  if (values.max(v)) {
    SetValueInt<P>(v, "max_", name_stem, descriptor, reflection, proto);
  }

  return 1;
}

template <typename P>
int
SetProtoIfSet(int value,
              int minval,
              const char * name,
              P& proto) {
  if (value < minval) {
    return 0;
  }

  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  const FieldDescriptor* fd = descriptor->FindFieldByName(name);
  reflection->SetUInt32(&proto, fd, value);

  return 1;
}

#define MATCHER_FROM_PROTO(proto, name, ztype, destination) { \
  if (proto.name ## _size() > 0) { \
    for (ztype v : proto.name()) { \
      destination.add(v); \
    } \
  } \
  if (proto.has_min_ ## name()) { \
    destination.set_min(proto.min_ ## name()); \
  } \
  if (proto.has_max_ ## name()) { \
    destination.set_max(proto.max_ ## name()); \
  } \
}

// Build a proto from a Matcher object.
#define PROTO_FROM_MATCHER(matcher, name, ztype, proto) { \
  for (const auto v : matcher.ValuesMatched()) { \
    proto.add_ ## name(v); \
  } \
  ztype s; \
  if (matcher.min(s)) { \
    proto.set_min_ ## name(s); \
  } \
  if (matcher.max(s)) { \
    proto.set_max_ ## name(s); \
  } \
}


// Convert between the enumeration in the proto and the operators
// needed by IW_Logical_Expression.
int
AddOperator(const SubstructureSearch::Operator op,
            IW_Logical_Expression& destination)
{
  switch(op)
  {
    case SubstructureSearch::SS_OR:
      destination.add_operator(IW_LOGEXP_OR);
      break;
    case SubstructureSearch::SS_AND:
      destination.add_operator(IW_LOGEXP_AND);
      break;
    case SubstructureSearch::SS_XOR:
      destination.add_operator(IW_LOGEXP_XOR);
      break;
    case SubstructureSearch::SS_LP_AND:
      destination.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
      break;
    default:
      cerr << "AddOperator:unrecognized operator " << op << '\n';
      return 0;
  }

  return 1;
}

//static int set_element_hits_needed_during_molecule_to_query = 1;

/*
  In several places we need to quickly ascertain whether or not an msi_object
  specifies a rejection or not
*/


// Elements can be specified as either numeric (atomic number) or string (atomic symbol)
// forms, or both.
static int
fetch_elements(const SubstructureSearch::SubstructureAtomSpecifier & proto,
               resizable_array<const Element*>& ele,
               resizable_array<int>& element_unique_id,
               extending_resizable_array<int>& element_uid)
{
  if (proto.atomic_symbol_size() > 0) {
    for (const auto& s : proto.atomic_symbol()) 
    {
      const_IWSubstring atomic_symbol(s);
      const Element * e = get_element_from_symbol_no_case_conversion(atomic_symbol);
      if (nullptr == e && auto_create_new_elements())
        e = create_element_with_symbol(atomic_symbol);

      if (nullptr == e) {
        cerr << "fetch_elements:no element for symbol '" << atomic_symbol << "'\n";
        return 0;
      }
      
      ele.add_if_not_already_present(e);
      element_unique_id.add_if_not_already_present(e->unique_id());
      element_uid[e->unique_id()] = 1;
    }
  }

  if (proto.atomic_number_size() > 0) {
    for (const auto z : proto.atomic_number()) {
      const Element * e = get_element_from_atomic_number(z);
      if (nullptr == e) {
        cerr << "fetch_elements:no element for atomic number " << z << '\n';
        return 0;
      }
      
      ele.add_if_not_already_present(e);
      element_unique_id.add_if_not_already_present(e->unique_id());
      element_uid[e->unique_id()] = 1;
    }
  }

  assert (ele.ok());

  return 1;
}


#define MIN_NOT_SPECIFIED 9393
#define MAX_NOT_SPECIFIED 14163

/*
  Possibly add `value` to `specifier` if it is within
  the range [min_val,max_val].
*/

template <typename T>
int
AppendIntValue(const T value,
               T min_val, T max_val,
               Min_Max_Specifier<T> & specifier)
{
  assert (specifier.ok());

  if (value < min_val || value > max_val)
  {
    cerr << "AppendIntValue:out of range " << value << " must bt btw " << min_val << " and " << max_val << '\n';
    return 0;
  }

  specifier.add(value);

  return 1;
}

// The objective is to fill `target` with the numbers in
// `values`, `minval` and `maxval`.
// Currently `attribute` is not used, nor are the max/min allowed values.
template <typename T>
bool
ReallyGruesome(const ::google::protobuf::RepeatedField<T>& values,
               const std::optional<T>& minval, const std::optional<T>& maxval,
               const T min_allowed, const T max_allowed,
               const char * attribute,
               Min_Max_Specifier<int>& target)
{
  if (minval.has_value())
    target.set_min(minval.value());

  if (maxval.has_value())
    target.set_max(maxval.value());

  for (const auto v : values) {
    target.add(v);
  }

  return true;
}

#define MAYBEMIN(p, attribute, T) (p.has_min_ ## attribute() ? std::optional<T>(p.min_ ## attribute()) : std::optional<T>())
#define MAYBEMAX(p, attribute, T) (p.has_max_ ## attribute() ? std::optional<T>(p.max_ ## attribute()) : std::optional<T>())

#define GETVALUES(p, attribute, lowest_allowed, max_allowed)    ReallyGruesome<uint32_t>(p.attribute(), MAYBEMIN(p, attribute, uint32_t), MAYBEMAX(p, attribute, uint32_t), lowest_allowed, max_allowed, #attribute, _ ## attribute)
#define GETVALUESINT(p, attribute, lowest_allowed, max_allowed) ReallyGruesome<int32_t>(p.attribute(),  MAYBEMIN(p, attribute, int32_t), MAYBEMAX(p, attribute, int32_t), lowest_allowed, max_allowed, #attribute, _ ## attribute)
#define GETFLOATVALUES(p, attribute, lowest_allowed, max_allowed) ReallyGruesome<float>(p.attribute(),  MAYBEMIN(p, attribute, float32), MAYBEMAX(p, attribute, float), lowest_allowed, max_allowed, #attribute, _ ## attribute)

using ::google::protobuf::RepeatedField;

template <typename T>
bool
really_gruesome(const RepeatedField<T> values, const T minval, const T maxval,
                const T min_value_allowed, const T max_value_allowed,
                Min_Max_Specifier<T>& destination, const char * name)
{
  if (values.size() > 0)
  {
    for (const T v : values)
    {
      if (v < min_value_allowed || v > max_value_allowed)
      {
        cerr << "really_gruesome::value of " << name << " " << v << " out of range, must be btw " << min_value_allowed << " and " << max_value_allowed << "\n";
        return false;
      }

      destination.add_if_not_already_present(v);
    }

    return true;
  }

  if (minval < min_value_allowed)
  {
    cerr << "really_gruesome::value of " << name << " " << minval << " out of range, must be less than " << min_value_allowed << "\n";
    return false;
  }

  if (minval > max_value_allowed)
  {
    cerr << "really_gruesome::value of " << name << " " << maxval << " out of range, must be more than " << max_value_allowed << "\n";
    return false;
  }

  destination.set_min(minval);
  destination.set_min(maxval);

  return true;
}

/*
  Process a MinMaxSpecifierUint or MinMaxSpecifierint proto (type P).
*/



int
Substructure_Atom_Specifier::construct_from_proto(const SubstructureSearch::SubstructureAtomSpecifier & proto)
{
  assert (ok());

//cerr << "SubstructureAtomSpecifier::construct_from_proto:from " << proto.ShortDebugString() << '\n';
  if (! fetch_elements(proto, _element, _element_unique_id, _element_uid)) {
    return 0;
  }

  if (proto.has_preference_value())
    _preference_value = proto.preference_value();

  // Should relate these to the min and max reasonable formal charge values.
  if (!GETVALUESINT(proto, formal_charge, -12, 12))
    return 0;

  MATCHER_FROM_PROTO(proto, ncon, int, _ncon);
//if (! GETVALUES(proto, ncon, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, nbonds, int, _nbonds);
//if (!GETVALUES(proto, nbonds, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, valence, int, _valence);
//if (!GETVALUES(proto, valence, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, nrings, int, _nrings);
//if (!GETVALUES(proto, nrings, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, ring_bond_count, int, _ring_bond_count);
//if (!GETVALUES(proto, ring_bond_count, 0, no_limit))
//  return 0;

  if (!GETVALUES(proto, ncon2, 0, no_limit))
    return 0;
  MATCHER_FROM_PROTO(proto, hcount, int, _hcount);
//if (!GETVALUES(proto, hcount, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, ring_size, int, _ring_size);
//if (!GETVALUES(proto, ring_size, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, attached_heteroatom_count, int, _attached_heteroatom_count);

  MATCHER_FROM_PROTO(proto, lone_pair_count, int, _lone_pair_count);
//if (!GETVALUES(proto, lone_pair_count, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, unsaturation, int, _unsaturation);
//if (!GETVALUES(proto, unsaturation, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, daylight_x, int, _daylight_x);
//if (!GETVALUES(proto, daylight_x, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, isotope, uint32_t, _isotope);
//if (!GETVALUES(proto, isotope, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, aryl, int, _aryl);

  MATCHER_FROM_PROTO(proto, vinyl, int, _vinyl);

  MATCHER_FROM_PROTO(proto, fused_system_size, int, _fused_system_size);
//if (!GETVALUES(proto, fused_system_size, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, symmetry_degree, int, _symmetry_degree);
//if (!GETVALUES(proto, symmetry_degree, 0, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, scaffold_bonds_attached_to_ring, int, _scaffold_bonds_attached_to_ring);
//if (!GETVALUES(proto, scaffold_bonds_attached_to_ring, 0, no_limit))
//  return 0;

  if (proto.has_symmetry_group())
  {
    _symmetry_group = proto.symmetry_group();
    if (_symmetry_group == 0) {
      cerr << "SubstructureAtomSpecifier::construct_from_proto:_symmetry_group must be > 0\n";
      return 0;
    }
  }

  if (proto.has_aromatic()) {
    if (proto.aromatic())
      _aromaticity = AROMATIC;
    else
      _aromaticity = NOT_AROMATIC;
  }

  if (proto.has_match_spinach_only())
    _match_spinach_only = proto.match_spinach_only();

  if (proto.has_chirality())
    _chirality = proto.chirality();

  MATCHER_FROM_PROTO(proto, heteroatoms_in_ring, int, _heteroatoms_in_ring);
//if (!GETVALUES(proto, heteroatoms_in_ring, 0, no_limit))
//  return 0;

  // Min should be related to min aromatic ring size.
  MATCHER_FROM_PROTO(proto, aromatic_ring_size, int, _aromatic_ring_size);
//if (!GETVALUES(proto, aromatic_ring_size, 4, no_limit))
//  return 0;

  MATCHER_FROM_PROTO(proto, aliphatic_ring_size, int, _aliphatic_ring_size);
//if (!GETVALUES(proto, aliphatic_ring_size, 3, no_limit))
//  return 0;

  if (proto.has_all_rings_kekule()) 
    _all_rings_kekule = proto.all_rings_kekule();

  if (proto.has_user_atom_type())
    _userAtomType = proto.user_atom_type();

  if (proto.has_atom_type())
    _atom_type = proto.atom_type();

  if (! proto.has_spiro()) {
    _spiro = -1;
  } else if (proto.spiro()) {
    _spiro = 1;
  } else {
    _spiro = 0;
  }

  if (! proto.has_cip()) {
  } else if (proto.cip() == SubstructureSearch::CahnIngoldPrelog::CIP_R) {
    _cip_chirality = CahnIngoldPrelog::R;
  } else if (proto.cip() == SubstructureSearch::CahnIngoldPrelog::CIP_S) {
    _cip_chirality = CahnIngoldPrelog::S;
  } else if (proto.cip() == SubstructureSearch::CahnIngoldPrelog::CIP_NEITHER) {
    _cip_chirality = CahnIngoldPrelog::kNeither;
  }
 
  (void) count_attributes_specified();

  assert (ok());

  return 1;
}

int
Substructure_Atom::_create_preference_from_proto(const SubstructureSearch::SubstructureAtomSpecifier& proto)
{
  if (! proto.has_preference_value()) {
    cerr << "Substructure_Atom_Specifier::_create_preference_from_proto:no preference value\n";
    return 0;
  }
//cerr << "_create_preference_from_proto " << proto.ShortDebugString() << '\n';

  std::unique_ptr<Substructure_Atom_Specifier> a = std::make_unique<Substructure_Atom_Specifier>();

  if (! a->construct_from_proto(proto)) {
    cerr << "Substructure_Atom_Specifier::_create_preference_from_proto:invalid preference specification\n";
    return 0;
  }

  _preferences.add(a.release());

  return 1;
}

/*
  There is a common function of examining an object and extracting its
  OPERATOR attribute and then adding the appropriate operator to a
  logical_expression
*/

static int
ExtractOperator(const google::protobuf::RepeatedField<int>& ops,
                const int ops_needed,
                IW_Logical_Expression& logexp,
                int default_operator,
                const char * caller)
{
  if (ops.empty()) {
    for (int i = 0; i < ops_needed; ++i) {
      logexp.add_operator(default_operator);
    }

    return 1;
  }

  if (ops.size() != ops_needed) {
    cerr << "ExtractOperator::operator count mismatch, needed " << ops_needed << " got " << ops.size() << "\n";
    return 0;
  }

  for (const auto op : ops)
  {
    switch (op) {
      case SubstructureSearch::Operator::SS_OR:
        logexp.add_operator(IW_LOGEXP_OR);
        break;
      case SubstructureSearch::Operator::SS_AND:
        logexp.add_operator(IW_LOGEXP_AND);
        break;
      case SubstructureSearch::Operator::SS_XOR:
        logexp.add_operator(IW_LOGEXP_XOR);
        break;
      case SubstructureSearch::Operator::SS_LP_AND:
        logexp.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
        break;
      default:
        cerr << "extract_operator:unrecognized operator " << op << "\n";
        return 0;
    };
  }

#ifdef DEBUG_EXTRACT_OPERATOR
  cerr << "extract_operator built\n";
  logexp.debug_print(cerr);
#endif

  return 1;
}

/*
  Parse input that looks like:

    - 1
    -!@ 1 2 3

    The bond, then a set of atoms.
*/

static int
parse_smarts_bond_attribute(const IWString& attribute,
                            const atom_number_t must_not_be,
                            Set_of_Atoms & other_end,
                            IWString & bond_smarts)
{
  const int nw = attribute.nwords();

  if (nw < 2)
  {
    cerr << "parse_smarts_bond_attribute: bond smarts must have at least two tokens '" << attribute << "'\n";
    return 0;
  }

  const int natoms = nw-1;  // first token is the bond

  other_end.resize_keep_storage(natoms);

  int i = 0;
  attribute.nextword(bond_smarts, i);

  const_IWSubstring s;
  while (attribute.nextword(s, i))
  {
    atom_number_t zatom;
    if (! s.numeric_value(zatom) || zatom < 0)
    {
      cerr << "parse_smarts_bond_attribute: invalid atom number '" << s << "'\n";
      return 0;
    }

    if (zatom == must_not_be) 
    {
      cerr << "parse_smarts_bond_attribute:self bonds not allowed\n";
      return 0;
    }

    other_end.add_if_not_already_present(zatom);
  }

  return bond_smarts.length();
}

bond_type_t
BondTypeFromSSBond(const google::protobuf::RepeatedField<int>& bond_type)
{
  bond_type_t btype = 0;

  for (const auto b : bond_type)
  {
    switch (b) {
      case SubstructureSearch::BondType::SS_SINGLE_BOND:
        btype |= SINGLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_DOUBLE_BOND:
        btype |= DOUBLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_TRIPLE_BOND:
        btype |= TRIPLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_AROMATIC_BOND:
        btype |= AROMATIC_BOND;
        break;
      case SubstructureSearch::BondType::SS_BOND:
        btype |= (SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND);
        break;
      default:
        cerr << "Substructure_Environment::BondTypeFromSSBond:unrecognised bond type " << b << "\n";
        return 0;
    };
  }

  return btype;
}

int
Substructure_Atom::_process_substructure_bond(const SubstructureSearch::SubstructureBond& bond,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 == bond.btype_size()) {
    cerr << "Substructure_Atom::_process_substructure_bond:no bonds\n";
    return 0;
  }

  if (! bond.has_other_end()) {
    cerr << "Substructure_Atom::_process_substructure_bond:no other end\n";
    return 0;
  }

  const atom_number_t other_end = bond.other_end();
  if (nullptr == completed[other_end]) {
    cerr << "Substructure_Atom::_process_substructure_bond:other end atom " << other_end << " not defined\n";
    return 0;
  }

  const bond_type_t btype = BondTypeFromSSBond(bond.btype());

  if (0 == btype) {
    cerr << "Substructure_Atom::_process_substructure_bond:invalid bond type specification\n";
    cerr << bond.ShortDebugString() << '\n';
    return 0;
  }

  // cerr << "Substructure_Atom::_process_substructure_bond:adding bond from " << _unique_id << " to " << other_end << " completed " << completed[other_end]->unique_id() << '\n';

  Substructure_Bond * b = new Substructure_Bond;

  b->set_bond_type(btype);
  b->set_atom(completed[other_end]);

  return _add_bond(b);
}

int
Substructure_Atom::_process_smarts_bond(const IWString& input,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
  Set_of_Atoms other_end;
  IWString bond_smarts;

  if (! parse_smarts_bond_attribute(input, _unique_id, other_end, bond_smarts))
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond:invalid bond smarts " << input << '\n';
    return 0;
  }

  for (const atom_number_t o : other_end)
  {
    if (nullptr == completed[o])
    {
      cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom " << o << " has not been defined\n";
      return 0;
    }

    if (this == completed[o])
    {
      cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom bonded to itself\n";
      return 0;
    }

    Substructure_Bond * b = new Substructure_Bond;

    int characters_processed;
    if (! b->construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
    {
      cerr << "parse_smarts_bond_attribute:cannot parse bond smarts " << input << "\n";
      delete b;
      return 0;
    }

    if (characters_processed != bond_smarts.length())
    {
      cerr << "parse_smarts_bond_attribute:extra junk at end of bond smarts " << input << "\n";
      cerr << characters_processed << " characters processed\n";
      delete b;
      return 0;
    }

    b->set_atom(completed[o]);

    _add_bond(b);
    if (_parent == nullptr) {
      _parent = completed[o];
    }
  }

  return 1;
}

/*
  An attribute bond will have the other end as the first int value
*/

// substructure_bond: "-!@ 1 2 3"

int
Substructure_Environment::_process_attachment_via_substructure_bond(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  Set_of_Atoms attachments;
  IWString bond_smarts;
  if (!parse_smarts_bond_attribute(proto.substructure_bond(), -1, attachments, bond_smarts)) {
    cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid substructure bond " << proto.substructure_bond() << "'\n";
    return 0;
  }

  for (const auto zatom : attachments) {
    if (nullptr == completed[zatom]) {
      cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid atom " << zatom <<"\n";
      return 0;
    }
    _possible_parents.add(completed[zatom]);
  }

  int characters_processed;
  if (! _bond.construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
  {
    cerr << "Substructure_Environment::_process_attribute_bond:invalid bond smarts " << bond_smarts << '\n';
    return 0;
  }

  if (characters_processed != bond_smarts.length()) {
    cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid bond smarts '" << bond_smarts << "'\n";
    return 0;
  }

  return 1;
}

/*
  The environment has just one bond. If you want the environment
  multiply attached, add a ring bond to the first atom in the env
*/

int
Substructure_Environment::_process_how_connected(const SubstructureSearch::SubstructureEnvironment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (!proto.has_attachment()) {
    cerr << "Substructure_Environment::_process_how_connected:no attachment\n";
    return 0;
  }

  const SubstructureSearch::EnvironmentAttachment& attachment = proto.attachment();

  if (attachment.attachment_point_size() > 0)
    return _process_attachment_bonds(attachment, completed);
  else if (attachment.has_substructure_bond())
    return _process_attachment_via_substructure_bond(attachment, completed);

  cerr << "Substructure_Atom::_process_how_connected:environment must have attachment\n";
  return 0;
}

int 
Substructure_Environment::_process_attachment_bonds(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 == proto.attachment_point_size() || 0 == proto.btype_size()) {
    cerr << "Substructure_Environment::_process_attachment_bonds:must have both attachment and bond attributes\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  Set_of_Atoms attachments;
  for (const auto a : proto.attachment_point())
  {
    if (attachments.add_if_not_already_present(a)) {
      if (nullptr == completed[a]) {
        cerr << "Substructure_Environment::_process_attachment_bonds:unrecognized atom number " << a << "\n";
        return 0;
      }
      _possible_parents.add(completed[a]);
    }
  }

  const bond_type_t btype = BondTypeFromSSBond(proto.btype());

  _bond.set_bond_type(btype);

  return 1;
}

/*
  When creating an environment, the environment components will
  register themselves as children. They are really environment components,
  so after they are done, we move any extra children to components.
*/

int
Substructure_Atom::_create_environment_from_proto(const SubstructureSearch::SubstructureAtomEnvironment& proto,
                                           extending_resizable_array<Substructure_Atom *> & completed)
{
  assert (_environment.empty());

  int initial_children = _children.number_elements();

  for (const auto& substructure_atom : proto.substructure_atom())
  {
    Substructure_Atom * a = new Substructure_Atom;
    if (! a->construct_from_proto(substructure_atom, completed))
    {
      cerr << "_create_environment_from_proto: cannot create environment component " << substructure_atom.ShortDebugString() << '\n';
      delete a;
      return 0;
    }
  }

  int final_children = _children.number_elements();
  if (initial_children == final_children)
  {
    cerr << "Substructure_Atom::create_env_from_proto: children not attached to root!\n";
    cerr << proto.ShortDebugString();
    return 0;
  }

  while (_children.number_elements() > initial_children)
    _environment.transfer_in(_children, initial_children);

  return _environment.number_elements();
}


// Make the assignment target=proposed if fn(proposed) is OK.
template <typename T, typename F>
bool
AssignValue(T proposed, T & target, F fn, const char * variable_name)
{
  if (! fn(proposed)) {
    cerr << "AssignValue:invalid " << variable_name  << ' ' << proposed << "\n";
    return false;
  }

  target = proposed;

  return true;
}
 
/*
  The sub-objects of a Substructure_Atom will be either bonds, or
  the Substructure_Atoms which are its components.
  Note that components are not allowed to have bonds.
*/

int
Substructure_Atom::construct_from_proto(const SubstructureSearch::SubstructureAtom& proto,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
//#define DEBUG_COMPLTED
#ifdef DEBUG_COMPLTED
  cerr << "Building Substructure_Atom " << proto.ShortDebugString() << '\n';
  for (int i = 0; i < 6; i++)
  {
    cerr << "Completed[" << i << "] = " << completed[i] << '\n';
  }
#endif

  if (!proto.has_id()) {
    cerr << "Substructure_Atom::_construct_from_proto:no id " << proto.ShortDebugString() << "\n";
    return 0;
  }

  _unique_id = proto.id();

  if (proto.has_match_as_match())
    _match_as_match_or_rejection = proto.match_as_match();

  if (nullptr != completed[_unique_id]) {
    cerr << "Substructure_Atom::construct_from_msi_object: Yipes, atom " << _unique_id << " already allocated\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_initial_atom_number())
    _initial_atom_number = proto.initial_atom_number();
  else
    _initial_atom_number = _unique_id;   // Not sure if this is the right thing to do or not.

  if (proto.has_ring_id() && ! AssignValue(proto.ring_id(), _ring_id, [](int r) { return r > 0;}, "ring_id"))
    return 0;

  if (proto.has_fragment_id() && ! AssignValue(proto.fragment_id(), _fragment_id, [](int f) { return f > 0;}, "fragment_id"))
    return 0;

  if (proto.has_global_match_id()) {
    _global_match_id = proto.global_match_id();
  }

  if (proto.has_fused_system_id())
    _fused_system_id = proto.fused_system_id();

  if (proto.has_text_identifier())
    _text_identifier = proto.text_identifier();

  if (proto.has_numeric_value())
    _numeric_value = proto.numeric_value();

  if (proto.has_include_in_embedding())
    _include_in_embedding = proto.include_in_embedding();

  if (proto.has_sum_all_preference_hits())
    _sum_all_preference_hits = proto.sum_all_preference_hits();

  if (proto.atom_properties_size() == 1) {
    const auto& spec = proto.atom_properties(0);
    Substructure_Atom_Specifier* tmp = this;
    if (! tmp->construct_from_proto(spec)) {
      cerr << "Substructure_Atom::_construct_from_proto:cannot parse Substructure_Atom_Specifier\n";
      cerr << spec.ShortDebugString() << "\n";
      return 0;
    }
  } else {
    for (int i = 0; i < proto.atom_properties_size(); ++i) {
      const auto& spec = proto.atom_properties(i);

      if (i > 0 && ! spec.has_logical_operator()) {
        cerr << "Substructure_Atom::_construct_from_proto:second Substructure_Atom_Specifier's must have operator\n";
        cerr << spec.ShortDebugString();
        return 0;
      }

      std::unique_ptr<Substructure_Atom_Specifier> tmp = std::make_unique<Substructure_Atom_Specifier>();
      if (!tmp->construct_from_proto(spec)) {
        cerr << "Substructure_Atom::_construct_from_proto:cannot parse Substructure_Atom_Specifier\n";
        cerr << spec.ShortDebugString() << "\n";
        return 0;
      }

      _components.add(tmp.release());
      if (i > 0)
        AddOperator(spec.logical_operator(), _operator);
    }
  }

// For simplicity, one can specify only one of atom smarts, smiles and smarts

  const int smiles_and_smarts = proto.has_atom_smarts() +
                                proto.has_smarts() +
                                proto.has_smiles();

  if (_components.number_elements() > 0 && smiles_and_smarts) {
    cerr << "Substructure_Atom::_construct_from_proto:cannot mix Substructure_Atom_Specifier with smiles/smarts specifications\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (smiles_and_smarts > 1) {
    cerr << "Substructure_Atom::_construct_from_proto:can have only one of 'atom_smarts', 'smarts' or 'smiles'\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_atom_smarts()) {
    IWString csmarts = proto.atom_smarts();

    if (! (csmarts.starts_with('[') && csmarts.ends_with(']'))) {
      cerr << "Substructure_Atom::_construct_from_proto: atomic smarts must be within [] '" << csmarts << "'\n";
      return 0;
    }

    if (! construct_from_smarts_token(csmarts)) {
      cerr << "Substructure_Atom::_construct_from_proto: cannot interpret smarts '" << csmarts << "'\n";
      return 0;
    }
  }

  if (proto.has_smiles()) {
    const IWString smiles = proto.smiles();
    if (! parse_smiles_specifier(smiles)) {
      cerr << "Substructure_Atom::_construct_from_proto:invalid smiles '" << proto.smiles() << "'\n";
      return 0;
    }
  }

  if (proto.has_smarts()) {
    const_IWSubstring smarts = proto.smarts();
    if (! parse_smarts_specifier(smarts)) {
      cerr << "Substructure_Atom::_construct_from_proto:invalid smarts '" << proto.smarts() << "'\n";
      return 0;
    }
  }

  completed[_unique_id] = this;

  for (const auto& env : proto.environment()) {
    if (! _create_environment_from_proto(env, completed)) {
      cerr << "Substructure_Atom::_construct_from_proto:invalid environment\n";
      cerr << env.ShortDebugString() << "\n";
      return 0;
    }
  }

  for (const auto & pref : proto.preference()) {
    if (! _create_preference_from_proto(pref)) {
      cerr << "Substructure_Atom::_construct_from_proto:invalid preference\n";
      cerr << pref.ShortDebugString() << "\n";
      return 0;
    }
  }

  if (proto.single_bond_size() > 0) {
    if (! _add_bonds(proto.single_bond(), SINGLE_BOND, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid single bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.double_bond_size() > 0) {
    if (! _add_bonds(proto.double_bond(), DOUBLE_BOND, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid double bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.triple_bond_size() > 0) {
    if (! _add_bonds(proto.triple_bond(), TRIPLE_BOND, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid triple bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.aromatic_bond_size() > 0) {
    if (! _add_bonds(proto.aromatic_bond(), AROMATIC_BOND, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid aromatic bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.bond_size() > 0) {
    if (! _add_bonds(proto.bond(), SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid any bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.bond_smarts_size() > 0) {
    for (const auto& smt : proto.bond_smarts()) {
      if (! _process_smarts_bond(smt, completed)) {
        cerr << "Substructure_Atom::construct_from_proto:invalid bond smarts bond\n";
        cerr << proto.ShortDebugString() << '\n';
        return 0;
      }
    }
  }

#ifdef QUERY_BONDS_PROCESSED_GREEDY
  // These are not being done after all atoms have been instantiated.
  if (proto.query_bond_size() > 0) {
    for (const auto & bond : proto.query_bond()) {
      if (! _process_substructure_bond(bond, completed)) {
        cerr << "Substructure_Atom::construct_from_proto:invalid substructure bond\n";
        cerr << proto.ShortDebugString() << '\n';
        return 0;
      }
    }
  }
#endif

  if (!GETVALUES(proto, unmatched_atoms_attached, 0, no_limit))
    return 0;

  for (const SubstructureSearch::SubstructureAtom& c : proto.children()) {
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();

    // Each env gets its own unique numbers and ids.
    extending_resizable_array<Substructure_Atom *> completed;
    if (! a->construct_from_proto(c, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid env " << c.ShortDebugString() << '\n';
      return 0;
    }
    _children << a.release();
  }

  if (proto.has_atom_type_group())
    _atom_type_group = proto.atom_type_group();

  assert (ok());
  return 1;
}

// Bonds are discerned after all the atoms have been constructed.
int
Substructure_Atom::FormBonds(const SubstructureSearch::SubstructureAtom& proto,
                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (proto.query_bond_size() == 0) {
    return 1;
  }

//#define DEBUG_FORM_BONDS
#ifdef DEBUG_FORM_BONDS
  cerr << "Substructure_Atom::FormBonds:adding " << proto.query_bond_size() << " query bonds\n";
#endif
  for (const auto & bond : proto.query_bond()) {
    if (! _process_substructure_bond(bond, completed)) {
      cerr << "Substructure_Atom::construct_from_proto:invalid substructure bond\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

// Add bond(s) of type `btype` to the query atoms in `atoms`.
int
Substructure_Atom::_add_bonds(const google::protobuf::RepeatedField<uint32_t>& atoms,
    const bond_type_t btype,
    extending_resizable_array<Substructure_Atom *> & completed)
{
  for (const auto a : atoms) {
    if (nullptr == completed[a]) {
      cerr << "Substructure_Atom::_add_bonds:non existent atom " << a << '\n';
      return 0;
    }

    if (static_cast<int>(a) == _unique_id) {
      cerr << "Substructure_Atom::_add_bonds:self bonds not allowed, atom " << a << '\n';
      return 0;
    }

    Substructure_Bond * b = new Substructure_Bond;
    b->set_atom(completed[a]);
    b->set_type(btype);
    _add_bond(b);
  }

  return 1;
}

int
Single_Substructure_Query::WriteProto(const char * fname)
{
  assert (nullptr != fname);

  std::ofstream os(fname, std::ios::out);
  if (! os.good()) {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteProto(os);
}

int
Single_Substructure_Query::WriteProto(std::ostream& output) 
{
  cerr << "Single_Substructure_Query::WriteProto:not implemented yet\n";
  return 0;
}

// Some fairly nasty preprocessor things. Later done via proto reflection.
// Not really clear what is the best approach.
// The preprocessor catches errors at compile time, whereas reflection errors
// are run-time.
#define ADDREPEATEDFIELD(p, field) { \
  for (auto value : _ ## field) { \
    proto.add_ ## field(value); \
  } \
}

// Populates a Min_Max_Specifier.
#define SETPROTOVALUES(p, attribute, type) { \
  ADDREPEATEDFIELD(p, attribute) \
  type tmp; \
  if (_ ## attribute.min(tmp)) \
    p.set_min_ ## attribute(tmp); \
  if (_ ## attribute.max(tmp)) \
    p.set_max_ ## attribute(tmp); \
}

// Convert from the operator types in IW_Logical_Expression to the
// SubstructureSearch::Operator
SubstructureSearch::Operator
OperatorTranslate(int op) {
  switch (op) {
    case IW_LOGEXP_AND:
      return SubstructureSearch::SS_AND;
    case IW_LOGEXP_OR:
      return SubstructureSearch::SS_OR;
    case IW_LOGEXP_XOR:
      return SubstructureSearch::SS_XOR;
    case IW_LOGEXP_LOW_PRIORITY_AND:
      return SubstructureSearch::SS_LP_AND;
    default:
      cerr << "Unrecognised operator " << op << '\n';
      return SubstructureSearch::SS_AND;
  }
}

int
Single_Substructure_Query::BuildProto(SubstructureSearch::SingleSubstructureQuery& proto) const
{
  if (_comment.length() > 0)
    proto.set_name(_comment.rawchars(), _comment.length());

  for (const auto* env : _environment) {
    env->BuildProto(*proto.add_environment());
  }

  if (! _originating_smarts.empty()) {
    proto.set_smarts(_originating_smarts.data(), _originating_smarts.size());
    return 1;
  }

  if (_find_one_embedding_per_start_atom) {
    proto.set_one_embedding_per_start_atom(_find_one_embedding_per_start_atom);
  }
  if (_normalise_rc_per_hits_needed > 0) {
    proto.set_normalise_rc_per_hits_needed(_normalise_rc_per_hits_needed);
  }
  if (_subtract_from_rc > 0) {
    proto.set_subtract_from_rc(_subtract_from_rc);
  }
  if (_max_matches_to_find > 0) {
    proto.set_max_matches_to_find(_max_matches_to_find);
  }
  if (!_save_matched_atoms) {
    proto.set_save_matched_atoms(_save_matched_atoms);
  }
  if (_ncon_ignore_singly_connected) {
    proto.set_ncon_ignore_singly_connected(_ncon_ignore_singly_connected);
  }
  if (_do_not_perceive_symmetry_equivalent_matches) {
    proto.set_perceive_symmetric_equivalents(false);
  } else {
    proto.set_perceive_symmetric_equivalents(true);
  }
  if (_implicit_ring_condition >= 0) {
    proto.set_implicit_ring_condition(_implicit_ring_condition);
  }
  if (_all_hits_in_same_fragment) {
    proto.set_all_hits_in_same_fragment(_all_hits_in_same_fragment);
  }
  if (_only_keep_matches_in_largest_fragment) {
    proto.set_only_match_largest_fragment(_only_keep_matches_in_largest_fragment);
  }
  if (_embeddings_do_not_overlap) {
    proto.set_embeddings_do_not_overlap(_embeddings_do_not_overlap);
  }
  if (_sort_by_preference_value) {
    proto.set_sort_by_preference_value(_sort_by_preference_value);
  }
  if (_fail_if_embeddings_too_close) {
    proto.set_fail_if_embeddings_too_close(_fail_if_embeddings_too_close);
  }

  ADDREPEATEDFIELD(proto, numeric_value);

  // TODO: implement this
#ifdef IMPLEMENT_NO_MATCHED_ATOMS_BETWEEN
  for (const auto * nmb : _no_matched_atoms_between) {
    SubstructureSearch::NoMatchedAtomsBetween * p = proto.add_no_matched_atoms_between();
    nmb->BuildProtoNoBond(*p);
  }
#endif
  if (_no_matched_atoms_between_exhaustive) {
    proto.set_no_matched_atoms_between_exhaustive(_no_matched_atoms_between_exhaustive);
  }
  if (! _environment_must_match_unmatched_atoms) {
    proto.set_environment_must_match_unmatched_atoms(_environment_must_match_unmatched_atoms);
  }
  if (_find_unique_embeddings_only) {
    proto.set_unique_embeddings_only(_find_unique_embeddings_only);
  }
  if (! _respect_initial_atom_numbering) {
    proto.set_respect_initial_atom_numbering(_respect_initial_atom_numbering);
  }
  if (_compress_embeddings) {
    proto.set_compress_embeddings(_compress_embeddings);
  }


  for (const auto * lnk : _link_atom) {
    SubstructureSearch::LinkAtoms* p = proto.add_link_atoms();
    lnk->BuildProto(*p);
  }

  if (_matched_atoms_to_check_for_hits_too_close > 0) {
    proto.set_distance_between_hits_ncheck(_matched_atoms_to_check_for_hits_too_close);
  }

  SetProtoValues(_attached_heteroatom_count, "attached_heteroatom_count", proto);
  SetProtoValues(_hits_needed, "hits_needed", proto);
  SetProtoValues(_ring_atoms_matched, "ring_atoms_matched", proto);
  SetProtoValues(_heteroatoms_matched, "heteroatoms_matched", proto);
  SetProtoValues(_heteroatoms_in_molecule, "heteroatoms_in_molecule", proto);
  SetProtoValues(_natoms, "natoms", proto);
  SetProtoValues(_nrings, "nrings", proto);
  SetProtoValues(_ncon, "ncon", proto);
  SetProtoValues(_fused_rings, "fused_rings", proto);
  SetProtoValues(_strongly_fused_rings, "strongly_fused_rings", proto);
  SetProtoValues(_isolated_rings, "isolated_rings", proto);
  SetProtoValues(_isolated_ring_objects, "isolated_ring_objects", proto);
  SetProtoValues(_aromatic_rings, "aromatic_rings", proto);
  SetProtoValues(_non_aromatic_rings, "non_aromatic_rings", proto);
  SetProtoValues(_distance_between_hits, "distance_between_hits", proto);
  SetProtoValues(_number_isotopic_atoms, "number_isotopic_atoms", proto);
  SetProtoValues(_number_fragments, "number_fragments", proto);
  SetProtoValues(_distance_between_root_atoms, "distance_between_root_atoms", proto);
  SetProtoValues(_atoms_in_spinach, "atoms_in_spinach", proto);
  SetProtoValues(_inter_ring_atoms, "inter_ring_atoms", proto);
  SetProtoValues(_unmatched_atoms, "unmatched_atoms", proto);
  SetProtoValues(_net_formal_charge, "net_formal_charge", proto);
  proto.set_environment_must_match_unmatched_atoms(_environment_must_match_unmatched_atoms);

  if (_min_fraction_atoms_matched > 0.0f) {
    proto.set_min_fraction_atoms_matched(_min_fraction_atoms_matched);
  }
  if (_max_fraction_atoms_matched > 0.0f) {
    proto.set_max_fraction_atoms_matched(_max_fraction_atoms_matched);
  }
  if (_min_all_matches_fraction_atoms_matched > 0.0f) {
    proto.set_min_all_matches_fraction_atoms_matched(_min_all_matches_fraction_atoms_matched);
  }
  if (_max_all_matches_fraction_atoms_matched > 0.0f) {
    proto.set_max_all_matches_fraction_atoms_matched(_max_all_matches_fraction_atoms_matched);

  }

  // TODO: Do something with _environments_can_share_attachment_points
  for (const auto* env : _environment) {
    env->BuildProto(*proto.add_environment());
  }
  for (const auto* env : _environment_rejections) {
    env->BuildProto(*proto.add_environment_no_match());
  }

  for (const auto * ring_sys : _ring_system_specification) {
    SubstructureSearch::SubstructureRingSystemSpecification * p = proto.add_ring_system_specifier();
    ring_sys->BuildProto(*p);
  }

  for (const auto op : _ring_specification_logexp.Operators()) {
    proto.add_ring_specification_logexp(OperatorTranslate(op));
  }

  for (const auto * ring : _ring_specification) {
    SubstructureSearch::SubstructureRingSpecification* p = proto.add_ring_specifier();
    ring->BuildProto(*p);
  }
  for (const auto op : _ring_system_specification_logexp.Operators()) {
    proto.add_ring_system_specifier_logexp(OperatorTranslate(op));
  }

  if (_find_unique_embeddings_only) {
    proto.set_unique_embeddings_only(_find_unique_embeddings_only);
  }

  _required_molecular_properties.BuildProto(*proto.mutable_required_molecular_properties());

  SetProtoValues(_aromatic_atoms, "aromatic_atoms", proto);
  ADDREPEATEDFIELD(proto, heteroatoms);
  if (_respect_initial_atom_numbering) {
    proto.set_respect_initial_atom_numbering(_respect_initial_atom_numbering);
  }
  if (_compress_embeddings) {
    proto.set_compress_embeddings(_compress_embeddings);
  }

  for (const auto * c: _chirality) {
    c->BuildProto(*proto.add_chiral_centre());
  }

  // atom type

  for (const auto* geom : _geometric_constraints) {
    geom->BuildProto(*proto.add_geometric_constraints());
  }

  for (const auto* sep : _separated_atoms) {
    sep->BuildProto(*proto.add_separated_atoms());
  }

  for (const InterRingAtoms* ira : _inter_ring_region) {
    ira->BuildProto(*proto.add_inter_ring_region());
  }

  for (const Substituent* substituent : _substituent) {
    substituent->BuildProto(*proto.add_substituent());
  }

  for (const DownTheBond* dtb : _down_the_bond) {
    dtb->BuildProto(*proto.add_down_the_bond());
  }

  for (const Region* region : _region) {
    SubstructureSearch::Region* r = proto.add_region();
    region->BuildProto(*r);
  }

  extending_resizable_array<Substructure_Atom*> atoms;
  for (Substructure_Atom * a : _root_atoms) {
    a->collect_all_atoms(atoms);
  }

  for (const Substructure_Atom * atom : atoms) {
    if (atom != nullptr) {
      atom->BuildProto(*proto.add_query_atom());
    }
  }

  return 1;
}

// Many operations involve checking if a value is > `minval` and if so
// then set the proto attribute.
#define SET_PROTO_IF_SET(p, attribute, minval) { \
  if (_ ## attribute > minval) { \
    p.set_ ## attribute(_ ## attribute); \
  } \
}

int
Substructure_Atom::BuildProto(SubstructureSearch::SubstructureAtom& proto) const {
  proto.set_id(_unique_id);  // 1

  if (! _match_as_match_or_rejection) {  // Default is true.
    proto.set_match_as_match(_match_as_match_or_rejection);  // 2
  }

  if (_text_identifier.length() > 0) {
    proto.set_text_identifier(_text_identifier.data(), _text_identifier.length());  // 3
  }

  SetProtoIfSet(_atom_map_number, 0, "atom_map_number", proto);  // 4
  SetProtoIfSet(_initial_atom_number, 0, "initial_atom_number", proto);  // 4
#ifdef OR_ID_NO_LONGER_PROCESSED
  SET_PROTO_IF_SET(proto, or_id, 0);  // 6
#endif
  SET_PROTO_IF_SET(proto, ring_id, 0);  // 9
  SET_PROTO_IF_SET(proto, fused_system_id, 0);  // 10
  SET_PROTO_IF_SET(proto, fragment_id, 0);  // 11
  SetProtoIfSet(_global_match_id, 0, "global_match_id", proto);  // 34

  double nv;
  if (_numeric_value.value(nv)) {  // 12
    proto.set_numeric_value(nv);
  }
  if (! _include_in_embedding) {
    proto.set_include_in_embedding(false);  // 13
  }

  //  smarts // 14
  //  atom_smarts // 15
  //  smiles // 16

  // environment looks hard.  // 17

  // The bonds are hard because sometimes they can be represented
  // as a SubstructureBond, and sometimes they need to be
  // represented as a bond smarts.
  for (const Substructure_Bond* b : _bonds) {  // 21
    b->BuildProto(proto);
  }

  // bond smarts 22

  // single, double triple aromatic, bond 25-29

  for (const Substructure_Atom_Specifier* a : _preferences) {  // 23
    a->BuildProto(*proto.add_preference());
  }

  if (_sum_all_preference_hits) {
    proto.set_sum_all_preference_hits(true);  // 24
  }

  if (_unmatched_atoms_attached.is_set()) {
    SETPROTOVALUES(proto, unmatched_atoms_attached, int);
  }

//SET_PROTO_IF_SET(proto, atom_type_group, 0);
  SetProtoIfSet(_atom_map_number, 0, "atom_type_group", proto);

  if (_components.number_elements() == 0) {
    const Substructure_Atom_Specifier* me = this;
    me->BuildProto(*proto.add_atom_properties());
  } else {
    for (const Substructure_Atom_Specifier* c : _components) {
      c->BuildProto(*proto.add_atom_properties());
    }
  }

  if (_components.number_elements() > 1) {
    for (int i = 0; i < _components.number_elements(); ++i) {
      SubstructureSearch::SubstructureAtomSpecifier* atom_prop = proto.add_atom_properties();
      _components[i]->BuildProto(*atom_prop);
      if (i == 0) {
        continue;
      }
      const int op = _operator.op(i - 1);
      // TODO(ianwatson) consolidate this logic into a separate function.
      if (op == IW_LOGEXP_AND) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_AND);
      } else if (op == IW_LOGEXP_OR) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_OR);
      } else if (op == IW_LOGEXP_XOR) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_XOR);
      } else if (op == IW_LOGEXP_LOW_PRIORITY_AND) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_LP_AND);
      } else {
        cerr << "Substructure_Atom::BuildProto:unrecognized operator type " << op << '\n';
        return 0;
      }
    }
  }

  _environment.BuildProto(proto);

  for (const Substructure_Atom* c : _children) {
    SubstructureSearch::SubstructureAtom* child = proto.add_children();
    c->BuildProto(*child);
  }

  return 1;
}

// If the bond looks simple, _b is null, then update the
// bonds in `proto`. Otherwise use a smarts string.
int
Substructure_Bond::BuildProto(SubstructureSearch::SubstructureAtom& proto) const {

  if (_b == nullptr) {
    SubstructureSearch::SubstructureBond* b = proto.add_query_bond();

    b->set_other_end(_a1->unique_id());

    if (_bond_types & SINGLE_BOND) {
      b->add_btype(SubstructureSearch::SS_SINGLE_BOND);
    }
    if (_bond_types & DOUBLE_BOND) {
      b->add_btype(SubstructureSearch::SS_DOUBLE_BOND);
    }
    if (_bond_types & TRIPLE_BOND) {
      b->add_btype(SubstructureSearch::SS_TRIPLE_BOND);
    }
    if (_bond_types & AROMATIC_BOND) {
      b->add_btype(SubstructureSearch::SS_AROMATIC_BOND);
    }

    return 1;
  }

  IWString smt = Smarts();
  smt << ' ' << _a1->unique_id();
  proto.add_bond_smarts(smt.data(), smt.size());

  return 1;
}

int
Substructure_Atom_Specifier::BuildProto(SubstructureSearch::SubstructureAtomSpecifier& proto) const {
  for (const Element* e : _element) {
    if (e->is_in_periodic_table()) {
      proto.add_atomic_number(e->atomic_number());
    } else {
      proto.add_atomic_symbol(e->symbol().data(), e->symbol().length());
    }
  }

  PROTO_FROM_MATCHER(_ncon, ncon, int, proto); // 3
  SetProtoValues(_ncon2, "ncon2", proto);  // 6
  PROTO_FROM_MATCHER(_nbonds, nbonds, int, proto); // 9
  PROTO_FROM_MATCHER(_valence, valence, int, proto); // 79
  SetProtoValuesInt(_formal_charge, "formal_charge", proto);  // 12
  PROTO_FROM_MATCHER(_nrings, nrings, int, proto); // 15
  PROTO_FROM_MATCHER(_ring_bond_count, ring_bond_count, int, proto);  // 15
  PROTO_FROM_MATCHER(_ring_size, ring_size, int, proto);  // 21
  PROTO_FROM_MATCHER(_hcount, hcount, int, proto); // 24
  if (_aromaticity == SUBSTRUCTURE_NOT_SPECIFIED) {
  }  else if (_aromaticity == AROMATIC) {  // 27
    proto.set_aromatic(true);
  }  else if (_aromaticity == NOT_AROMATIC) {
    proto.set_aromatic(false);
  }
  if (_chirality != SUBSTRUCTURE_NOT_SPECIFIED)
    proto.set_chirality(true);  // 28
  PROTO_FROM_MATCHER(_aromatic_ring_size, aromatic_ring_size, int, proto);  // 30
  PROTO_FROM_MATCHER(_aliphatic_ring_size, aliphatic_ring_size, int, proto);  // 33
  PROTO_FROM_MATCHER(_attached_heteroatom_count, attached_heteroatom_count, int, proto);  // 36
  PROTO_FROM_MATCHER(_lone_pair_count, lone_pair_count, int, proto);  // 39
  PROTO_FROM_MATCHER(_unsaturation, unsaturation, int, proto);  // 42
  PROTO_FROM_MATCHER(_daylight_x, daylight_x, int, proto);  // 45
  PROTO_FROM_MATCHER(_isotope, isotope, uint32_t, proto);  // 48
  PROTO_FROM_MATCHER(_aryl, aryl, int, proto);  // 51
  PROTO_FROM_MATCHER(_vinyl, vinyl, int, proto);  // 54
  PROTO_FROM_MATCHER(_fused_system_size, fused_system_size, int, proto);  // 58
  if (_all_rings_kekule != SUBSTRUCTURE_NOT_SPECIFIED) {
    proto.set_all_rings_kekule(true);  // 60
  }
  PROTO_FROM_MATCHER(_heteroatoms_in_ring, heteroatoms_in_ring, int, proto);  // 61
  SET_PROTO_IF_SET(proto, match_spinach_only, -1);  // 64
  PROTO_FROM_MATCHER(_scaffold_bonds_attached_to_ring, scaffold_bonds_attached_to_ring, int, proto);  // 65
  SET_PROTO_IF_SET(proto, preference_value, 0);  // 68
  PROTO_FROM_MATCHER(_symmetry_degree, symmetry_degree, int, proto);  // 69
  SET_PROTO_IF_SET(proto, symmetry_group, 0);  // 72
  if (_spiro < 0) {
  } else if (_spiro == 0) {
    proto.set_spiro(false);
  } else {
    proto.set_spiro(true);
  }

  switch (_cip_chirality) {
    case CahnIngoldPrelog::kUnspecified:
      break;
    case CahnIngoldPrelog::kNeither:
      proto.set_cip(SubstructureSearch::CIP_NEITHER);
      break;
    case CahnIngoldPrelog::R:
      proto.set_cip(SubstructureSearch::CIP_R);
      break;
    case CahnIngoldPrelog::S:
      proto.set_cip(SubstructureSearch::CIP_S);
      break;
    default:
      break;
  }

  // Not sure what to do with operator...
  // atom typing is not implemented, seems ambiguous.

  return 1;
}

int
Single_Substructure_Query::_parse_ring_specifier(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  for (const auto& ring_spec : proto.ring_specifier()) {
    std::unique_ptr<Substructure_Ring_Specification> r = std::make_unique<Substructure_Ring_Specification>();
    if (! r->ConstructFromProto(ring_spec)) {
      cerr << "SingleSubstructureQuery::_parse_ring_specifier:could not create ring specifier\n";
      cerr << ring_spec.ShortDebugString();
      return 0;
    }

    _ring_specification << r.release();
  }

  const int ring_specification_count = _ring_specification.number_elements();

  if (1 == ring_specification_count)    // no operators with first object
  {
    if (!proto.ring_specification_logexp().empty()) {
      cerr << "Single_Substructure_Query::_parse_ring_specifier_object: spurious operator\n";
      return 0;
    }
    return 1;
  }

  if (0 == proto.ring_specification_logexp_size())
  {
    for (int i = 0; i < ring_specification_count - 1; ++i)
    {
      _ring_specification_logexp.add_operator(IW_LOGEXP_AND);
    }

    return 1;
  }

  if (ring_specification_count != proto.ring_specification_logexp_size() + 1) {
    cerr << "Single_Substructure_Query::_parse_ring_specifier_object:inconsistent ring specifications " <<
             ring_specification_count << " and operators " << proto.ring_specification_logexp_size() << '\n';
    return 0;
  }

  if (! ExtractOperator(proto.ring_specification_logexp(), ring_specification_count - 1,
          _ring_specification_logexp, IW_LOGEXP_AND,
          "Single_Substructure_Query::_parse_ring_specifier_object"))
    return 0;

  return 1;
}

int
Single_Substructure_Query::_parse_ring_system_specifier(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  if (0 == proto.ring_system_specifier_size())
    return 1;

  for (const auto& spec : proto.ring_system_specifier())
  {
    Substructure_Ring_System_Specification * r = new Substructure_Ring_System_Specification();
    if (! r->ConstructFromProto(spec))
    {
      delete r;

      cerr << "SingleSubstructureQuery::_parse_ring_system_specifier:could not create ring system specifier object from proto\n";
      cerr << spec.ShortDebugString() << "\n";
      return 0;
    }

    _ring_system_specification.add(r);
  }

  const int ring_system_spec_count = _ring_system_specification.number_elements();

  if (1 == ring_system_spec_count)  // No operators to worry about.
  {
    if (! proto.ring_system_specifier_logexp().empty()) {
      cerr << "Single_Substructure_Query::_parse_ring_system_specifier:spurious ring_system_specifier_logexp\n";
      return 0;
    }
    return 1;
  }

  if (proto.ring_system_specifier_logexp().empty())
  {
    for (int i = 0; i < ring_system_spec_count - 1; ++i) {
      _ring_system_specification_logexp.add_operator(IW_LOGEXP_AND);
    }

    return 1;
  }

  if (ring_system_spec_count != proto.ring_system_specifier_logexp_size() + 1) {
    cerr << "Single_Substructure_Query::_parse_ring_system_specifier:inconsistent operator count " <<
       ring_system_spec_count << " vs " << proto.ring_system_specifier_logexp_size() << '\n';
    return 0;
  }

  if (! ExtractOperator(proto.ring_system_specifier_logexp(), ring_system_spec_count - 1,
        _ring_system_specification_logexp, IW_LOGEXP_AND, "Single_Substructure_Query::_parse_ring_system_specifier_object"))
    return 0;

  return 1;
}

/*
  Is an msi object a root atom or not.  Basically, if it has an
  attribute which ends in "bond" (so as to not also match "nbonds"),
  then it is not a root
*/

static bool
IsRootSubstructureAtom(const SubstructureSearch::SubstructureAtom& proto)
{
  if (proto.query_bond_size() > 0)
    return false;

  if (proto.bond_smarts_size() > 0)
    return false;

  if (proto.single_bond_size() > 0)
    return false;

  if (proto.double_bond_size() > 0)
    return false;

  if (proto.triple_bond_size() > 0)
    return false;

  if (proto.aromatic_bond_size() > 0)
    return false;

  if (proto.bond_size() > 0)
    return false;

  return true;
}

/*
  The environment object 
  It can specify query_atoms, smiles or smarts
*/

int
Single_Substructure_Query::_construct_environment_from_proto(
    const google::protobuf::RepeatedPtrField<SubstructureSearch::SubstructureEnvironment>& env,
    extending_resizable_array<Substructure_Atom *> & completed,
    resizable_array_p<Substructure_Environment> & destination)
{
  for (const auto& e : env)
  {
    std::unique_ptr<Substructure_Environment> a = std::make_unique<Substructure_Environment>();

    _collect_all_atoms(completed);

    int rc = a->construct_from_proto(e, completed);

    if (0 == rc) {
      cerr << "Single_Substructure_Query::_construct_environment_from_proto:invalid environment " << e.ShortDebugString() << "\n";
      return 0;
    }

    destination.add(a.release());
  }

  return 1;
}

int
Single_Substructure_Query::_construct_matched_atoms_match(
    const google::protobuf::RepeatedPtrField<SubstructureSearch::MatchedAtomMatch>& matched_atoms_match) {

  for (const auto & mam : matched_atoms_match) {
    std::unique_ptr<MatchedAtomMatch> m = std::make_unique<MatchedAtomMatch>();
    if (! m->ConstructFromProto(mam)) {
      cerr << "Single_Substructure_Query::_construct_matched_atoms_match:invalid " << mam.ShortDebugString() << '\n';
      return 0;
    }
    if (mam.has_logexp()) {
      AddOperator(mam.logexp(), _matched_atom_match_operator);
    } else if (_matched_atom_match.size() > 0) {
      AddOperator(SubstructureSearch::SS_AND, _matched_atom_match_operator);
    }
    _matched_atom_match << m.release();
  }

  return 1;
}

int
Substructure_Environment::BuildProto(SubstructureSearch::SubstructureEnvironment& proto) const
{
  proto.set_id(_unique_id);

  SET_PROTO_IF_SET(proto, or_id, 0);  // 8
  SET_PROTO_IF_SET(proto, and_id, 0);  // 9

  SETPROTOVALUES(proto, hits_needed, int);  // 10
  if (_no_other_substituents_allowed) {
    proto.set_no_other_substituents_allowed(true);  // 13
  }
  if (_environments_can_share_attachment_points) {
    proto.set_env_matches_can_share_attachment_points(true);  // 15
  }
  SET_PROTO_IF_SET(proto, max_matches_to_find, 0);  // 16
  if (_hydrogen_ok_as_environment_match) {
    proto.set_hydrogen_ok(true);  // 17
  }
  if (_max_environment_matches_per_attachment_point > 0) {
    proto.set_max_env_matches_per_anchor(_max_environment_matches_per_attachment_point);
  }

  // The rest of this is too hard.
  for (const Substructure_Atom* a : _possible_parents) {
    a->BuildProto(*proto.add_query_atom());
  }
  cerr << "Substructure_Environment::BuildProto:implement this sometime\n";

  return 1;
}

int
Substructure_Environment::construct_from_proto(const SubstructureSearch::SubstructureEnvironment& proto,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               atom_number_t possible_parent,
                               bond_type_t possible_parent_bond_type)
{
// Process any bonds which have been specified as attributes

  if (! _process_how_connected(proto, completed))
  {
    cerr << "Substructure_Environment::construct_from_proto:cannot connect " << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER != possible_parent)
    _add_possible_parent(possible_parent, possible_parent_bond_type, completed);

  if (_possible_parents.empty()) {
    cerr << "Substructure_Environment::construct_from_proto: environment not connected\n";
    return 0;
  }

  if (proto.hits_needed_size() > 0) {
    for (const auto x : proto.hits_needed()) {
      _hits_needed.add_if_not_already_present(x);
    }
  }
  if (proto.has_min_hits_needed()) {
    _hits_needed.set_min(proto.min_hits_needed());
  }
  if (proto.has_max_hits_needed()) {
    _hits_needed.set_max(proto.max_hits_needed());
  }

  if (!_hits_needed.ok()) {
    cerr << "Substructure_Environment::construct_from_proto:invalid hits needed\n";
    _hits_needed.debug_print(cerr);
  }

  if (proto.has_no_other_substituents_allowed())
    _no_other_substituents_allowed = proto.no_other_substituents_allowed();

  if (proto.has_hydrogen_ok())
    _hydrogen_ok_as_environment_match = proto.hydrogen_ok();

  if (proto.has_max_env_matches_per_anchor())
    _max_environment_matches_per_attachment_point = proto.max_env_matches_per_anchor();

  if (proto.has_env_matches_can_share_attachment_points())
    _environments_can_share_attachment_points = proto.env_matches_can_share_attachment_points();

  if (proto.has_max_matches_to_find())
    _max_matches_to_find = proto.max_matches_to_find();

  if (proto.has_or_id()) {
    _or_id = proto.or_id();

    if (_or_id <= 0) {
      cerr << "Substructure_Environment::construct_from_proto: or_id values must be whole +ve numbers " << _or_id << " invalid\n";
      return 0;
    }
  }

  if (proto.has_and_id()) {
    _and_id = proto.and_id();
    if (_and_id <= 0) {
      cerr << "Substructure_Environment::construct_from_proto: and_id values must be whole +ve numbers " << _and_id << " invalid\n";
      return 0;
    }
  }

  if (_and_id && _or_id) {
    cerr << "Substructure_Environment::construct_from_proto:cannot have both OR and AND specifications\n";
    return 0;
  }

// Now construct the structural specification. There can be any number

  for (const std::string& smarts: proto.smarts()) {
    if (smarts.empty()) {
      cerr << "Substructure_Environment::construct_from_proto:empty smarts\n";
      return 0;
    }
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();

    const_IWSubstring x = smarts;
    const char * s = smarts.data();

    if (x.length() > 0 && (isdigit(s[0]) || '>' == s[0] || '<' == s[0] || s[0] == kOpenBrace)) {
      int chars_consumed = substructure_spec::SmartsNumericQualifier(s, smarts.length(), _hits_needed);
      if (chars_consumed == 0) {
        cerr << "Substructure_Environment::construct_from_proto:invalid numeric qualifier '" << x << "'\n";
        return 0;
      }

      x += chars_consumed;
    }

    if (! a->parse_smarts_specifier(x)) {
      return 0;
    }

//  cerr << "Build query from smarts " << (*att) << '\n';
    add(a.release());
  }

  for (const auto& smiles : proto.smiles()) {
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();

    if (! a->parse_smiles_specifier(smiles)) {
      return 0;
    }

//  cerr << "Build query from smiles " << (*att) << '\n';
    add(a.release());
  }

  for (const auto& query_atom: proto.query_atom()) {
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();
    if (! a->construct_from_proto(query_atom, completed)) {
      return 0;
    }

    if (IsRootSubstructureAtom(query_atom)) {
      add(a.release());
    } else {  // Will be given to its parent below.
      a.release();
    }
  }

  // Now that all atoms are available in the completed array, form bonds.
  // Note that this only handles query_bond specifications, not single_bond, double_bond...
  for (const auto& atom : proto.query_atom()) {
    Substructure_Atom * a = completed[atom.id()];
    if (a == nullptr) {
      cerr << "NULL completed " << atom.id() << '\n';
      return 0;
    }
    if (! a->FormBonds(atom, completed)) {
      cerr << "Substructure_Environment::construct_from_proto:FormBonds failed\n";
      return 0;
    }
    if (!IsRootSubstructureAtom(atom) && a->nbonds() == 0) {
      cerr << "Non root Substructure_Atom not bonded, id " << a->unique_id() << "\n";
      cerr << atom.ShortDebugString() << "\n";
      return 0;
    }
  }

  return 1;
}


template <typename T>
void
transfer_to_our_array(resizable_array_p<T> & to,
                      const resizable_array<T *> & from)
{
  for (T* x : from)
  {
    to.add(x);
  }

  return;
}

template void transfer_to_our_array(resizable_array_p<Substructure_Atom> &, const resizable_array<Substructure_Atom *> &);
template void transfer_to_our_array(resizable_array_p<Bond> &, const resizable_array<Bond *> &);
template void transfer_to_our_array(resizable_array_p<Link_Atom> &, const resizable_array<Link_Atom *> &);

SeparatedAtoms::SeparatedAtoms() {
  _a1 = INVALID_ATOM_NUMBER;
  _a2 = INVALID_ATOM_NUMBER;
}

int
SeparatedAtoms::Build(const SubstructureSearch::SeparatedAtoms& proto) {
  _a1 = proto.a1();
  _a2 = proto.a2();
  if (_a1 == _a2) {
    cerr << "SeparatedAtoms::Build:separated atoms must be distinct\n";
    return 0;
  }
  if (proto.bonds_between().size() > 0) {
    for (auto s : proto.bonds_between()) {
      _separation.add_if_not_already_present(s);
    }
  }
  if (proto.min_bonds_between() > 0) {
    _separation.set_min(proto.min_bonds_between());
  }
  if (proto.max_bonds_between() > 0) {
    if (! _separation.set_max(proto.max_bonds_between())) {
      cerr << "SeparatedAtoms::Build:invalid max bonds between\n";
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER

/*
  Parsing a smarts is complicated by the possible presence of the '.' character,
  and the need to process the '...' construct (no_matched_atoms_between)

  Valid inputs are

  'C-O-C'
  'C-O-C.[NH2]'
  'Br...C(=O)Cl'

*/




template <typename FROM, typename TO>
int
FetchRepeatedField(const google::protobuf::RepeatedField<FROM> & values,
                   resizable_array<TO> & destination)
{
  for (const FROM z : values)
  {
    destination.add(static_cast<TO>(z));
  }

  return destination.number_elements();
}

int
Single_Substructure_Query::_construct_from_proto(const SubstructureSearch::SingleSubstructureQuery & proto)
{
  assert (ok());
  assert (_root_atoms.empty());

  _numeric_value.resize(0);

  int environments_can_share_attachment_points = -1;    // indicates not set

  if (proto.has_name())
    _comment = proto.name();
  else if (proto.has_label())
    _comment = proto.label();

  if (proto.has_subtract_from_rc())
    _subtract_from_rc = proto.subtract_from_rc();

  if (proto.has_respect_initial_atom_numbering())
    _respect_initial_atom_numbering = proto.respect_initial_atom_numbering();

  if (proto.has_max_matches_to_find())
    _max_matches_to_find = proto.max_matches_to_find();

  if (proto.has_save_matched_atoms())
    _save_matched_atoms = proto.save_matched_atoms();

  if (proto.has_ncon_ignore_singly_connected())
    _ncon_ignore_singly_connected = proto.ncon_ignore_singly_connected();

  if (proto.has_perceive_symmetric_equivalents())
    _do_not_perceive_symmetry_equivalent_matches = !proto.perceive_symmetric_equivalents();

  if (proto.has_embeddings_do_not_overlap())
    _embeddings_do_not_overlap = proto.embeddings_do_not_overlap();

  if (proto.has_implicit_ring_condition())
    _implicit_ring_condition = proto.implicit_ring_condition();

  if (proto.has_all_hits_in_same_fragment())
    _all_hits_in_same_fragment = proto.all_hits_in_same_fragment();

  if (proto.has_only_match_largest_fragment())
    _only_keep_matches_in_largest_fragment = proto.only_match_largest_fragment();

  if (proto.has_one_embedding_per_start_atom())
    _find_one_embedding_per_start_atom = proto.one_embedding_per_start_atom();

  if (proto.has_sort_by_preference_value())
    _sort_by_preference_value = proto.sort_by_preference_value();

  if (proto.numeric_value_size() > 0)
    FetchRepeatedField( proto.numeric_value(), _numeric_value);

  if (proto.has_fail_if_embeddings_too_close())
    _fail_if_embeddings_too_close = proto.fail_if_embeddings_too_close();

  if (proto.has_distance_between_hits_ncheck())
  {
    _matched_atoms_to_check_for_hits_too_close = proto.distance_between_hits_ncheck();
    if (_matched_atoms_to_check_for_hits_too_close == 0)
    {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid distance_between_hits_ncheck\n";
      return 0;
    }
  }

  if (proto.has_environment_must_match_unmatched_atoms())
    _environment_must_match_unmatched_atoms = proto.environment_must_match_unmatched_atoms();

  if (proto.has_environments_can_share_attachment_points())
    environments_can_share_attachment_points = proto.environments_can_share_attachment_points();

  if (proto.heteroatoms_size() > 0)
  {
    if (!FetchRepeatedField(proto.heteroatoms(), _heteroatoms))
    {
      cerr << "SingleSubstructureQuery::_construct_from_proto:invalid heteroatoms specification\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.has_required_molecular_properties()) {
    if (! _required_molecular_properties.ConstructFromProto(proto.required_molecular_properties())) {
      cerr << "SingleSubstructureQuery::_construct_from_proto:cannot build required molecular properties\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  // The smiles and smarts directives scramble the _unique_id.
  if (_respect_initial_atom_numbering &&
      (proto.has_smiles() || proto.has_smarts()))
  {
    cerr << "SingleSubstructureQuery::construct_from_proto:cannot use respect_initial_atom_numbering with smiles or smarts directives\n";
    return 0;
  }

  if (proto.has_smiles())
  {
    const IWString smiles = proto.smiles();
    Substructure_Atom * a = new Substructure_Atom;
    if (! a->parse_smiles_specifier(smiles))
    {
      delete a;
      cerr << "Single_Substructure_Query::_construct_from_proto: cannot parse '" << smiles << "'\n";
      return 0;
    }

    _root_atoms.add(a);
  }

  if (proto.has_smarts())
  {
    const_IWSubstring smarts = proto.smarts();

    if (! create_from_smarts(smarts))
    {
      cerr << "Single_Substructure_Query::_construct_from_proto: cannot parse '" << smarts << "'\n";
      return 0;
    }
  }

  if (proto.has_one_embedding_per_start_atom())
    set_find_one_embedding_per_atom(proto.one_embedding_per_start_atom());

  if (proto.has_unique_embeddings_only())
    set_find_unique_embeddings_only(proto.unique_embeddings_only());

  if (proto.has_normalise_rc_per_hits_needed())
    set_normalise_rc_per_hits_needed(proto.normalise_rc_per_hits_needed());

  if (proto.no_matched_atoms_between_size() > 0)
  {
    for (const auto& nma : proto.no_matched_atoms_between())
    {
      if (nma.a1() == nma.a2())
      {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid no matched atoms between " << nma.ShortDebugString() << "'\n";
        return 0;
      }

      Bond * b = new Bond(nma.a1(), nma.a2(), SINGLE_BOND);

      _no_matched_atoms_between.add(b);
    }
  }

  if (proto.link_atoms_size() > 0)
  {
    for (const auto& link : proto.link_atoms())
    {
      Link_Atom * l = new Link_Atom;

      if (! l->ConstructFromProto(link))
      {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid link atom specification\n";
        cerr << link.ShortDebugString() << '\n';
        delete l;
        return 0;
      }

      _link_atom.add(l);
    }
  }

  if (proto.down_the_bond_size() > 0) {
    for (const auto& dtb : proto.down_the_bond()) {
      std::unique_ptr<DownTheBond> down_the_bond(std::make_unique<DownTheBond>());
      if (! down_the_bond->ConstructFromProto(dtb)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid DownTheBond\n";
        cerr << dtb.ShortDebugString() << '\n';
        return 0;
      }
      _down_the_bond.add(down_the_bond.release());
    }
  }

  if (proto.substituent_size() > 0) {
    for (const auto& subst : proto.substituent()) {
      std::unique_ptr<Substituent> substituent(std::make_unique<Substituent>());
      if (! substituent->ConstructFromProto(subst)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid Substituent\n";
        cerr << subst.ShortDebugString() << '\n';
        return 0;
      }
      _substituent.add(substituent.release());
    }
  }

  if (proto.substituent_no_match_size() > 0) {
    for (const auto& subst : proto.substituent_no_match()) {
      std::unique_ptr<Substituent> substituent(std::make_unique<Substituent>());
      if (! substituent->ConstructFromProto(subst)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid Substituent\n";
        cerr << subst.ShortDebugString() << '\n';
        return 0;
      }
      _substituent_no_match.add(substituent.release());
    }
  }

  if (proto.region_size() > 0) {
    for (const auto& r : proto.region()) {
      std::unique_ptr<Region> region(std::make_unique<Region>());
      if (! region->ConstructFromProto(r)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid Region\n";
        cerr << r.ShortDebugString() << '\n';
        return 0;
      }
      _region.add(region.release());
    }
  }

  if (proto.nearby_atoms_size() > 0) {
    for (const auto& r : proto.nearby_atoms()) {
      std::unique_ptr<NearbyAtoms> nearby(std::make_unique<NearbyAtoms>());
      if (! nearby->ConstructFromProto(r)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid NearbyAtoms\n";
        cerr << r.ShortDebugString() << '\n';
        return 0;
      }
      _nearby_atoms.add(nearby.release());
    }
  }

  if (proto.has_sort_matches() > 0)
  {
    _sort_matches_by = 0;

    const IWString s(proto.sort_matches());

    int j = 0;
    const_IWSubstring token;
    while (s.nextword(token, j))
    {
      int direction = 1;
      if (token.ends_with('+'))
        token.chop();
      else if (token.ends_with('-'))
      {
        token.chop();
        direction = -1;
      }
      if ("ncon" == token)
      {
        if (1 == direction)
          _sort_matches_by |= SORT_MATCHES_BY_INCREASING_NCON;
        else
          _sort_matches_by |= SORT_MATCHES_BY_DECREASING_NCON;
      }
      else if ("maxd" == token)
      {
        if (1 == direction)
          _sort_matches_by |= SORT_MATCHES_BY_INCREASING_MAXD;
        else
          _sort_matches_by |= SORT_MATCHES_BY_DECREASING_MAXD;
      }
      else
      {
        cerr << "Single_Substructure_Query:construct_from_proto:unrecognised sort attribute '" << s << "'\n";
        return 0;
      }
    }
  }

  if (_only_keep_matches_in_largest_fragment && _all_hits_in_same_fragment)
  {
    cerr << "Single_Substructure_Query::_construct_from_proto:the _only_keep_matches_in_largest_fragment and _all_hits_in_same_fragment attributes are mutually inconsistent\n";
    return 0;
  }

  constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, attached_heteroatom_count, 0, no_limit))
    return 0;

  if (_hits_needed.is_set())   // from a numeric qualifier on a smarts
    ;
  else if (! GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, ring_atoms_matched, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, heteroatoms_matched, 0, no_limit))
    return 0;

  if (_heteroatoms.number_elements() && 
      (! _attached_heteroatom_count.is_set() && ! _heteroatoms_in_molecule.is_set()))
  {
    cerr << "Single_Substructure_Query::_construct_from_proto: " << _heteroatoms.number_elements() <<
            " heteroatoms defined, but no attached_heteroatom_count\n";
    return 0;
  }

  if (! GETVALUES(proto, distance_between_hits, 0, no_limit)) {
    return 0;
  }

  if (! GETVALUES(proto, ncon, 0, no_limit)) {
    return 0;
  }

#ifdef NOW_IN_REQUIRED_MOLECULAR_PROPERTIES
  if (! GETVALUES(proto, natoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, nrings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, heteroatoms_in_molecule, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, fused_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, strongly_fused_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, isolated_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, isolated_ring_objects, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, aromatic_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, aromatic_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, non_aromatic_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, number_isotopic_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, number_fragments, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, atoms_in_spinach, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, inter_ring_atoms, 0, no_limit))
    return 0;

  // Need to sync up with reasonable formal charge
  if (! GETVALUESINT(proto, net_formal_charge, -12, 12))
    return 0;
#endif

  if (! GETVALUES(proto, distance_between_root_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, unmatched_atoms, 0, no_limit))
    return 0;

  if (proto.has_min_fraction_atoms_matched())
    _min_fraction_atoms_matched = proto.min_fraction_atoms_matched();

  if (proto.has_max_fraction_atoms_matched())
    _max_fraction_atoms_matched = proto.max_fraction_atoms_matched();

  if (proto.has_min_fraction_atoms_matched() && !proto.has_max_fraction_atoms_matched())
    _max_fraction_atoms_matched = 1.0f;
  else if (_min_fraction_atoms_matched > _max_fraction_atoms_matched)
  {
    cerr << "Single_Substructure_Query::_construct_from_proto:invalid min/max fraction atoms matched\n";
    cerr << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (proto.has_min_all_matches_fraction_atoms_matched()) {
    _min_all_matches_fraction_atoms_matched = proto.min_all_matches_fraction_atoms_matched();
  }

  if (proto.has_max_all_matches_fraction_atoms_matched()) {
    _max_all_matches_fraction_atoms_matched = proto.max_all_matches_fraction_atoms_matched();
  }

  if (_min_all_matches_fraction_atoms_matched > 0.0f &&
      _max_all_matches_fraction_atoms_matched == 0.0f) {
    _max_all_matches_fraction_atoms_matched = 1.0f;
  }

  if (_min_all_matches_fraction_atoms_matched >
      _max_all_matches_fraction_atoms_matched) {
    cerr << "Single_Substructure_Query::_construct_from_proto:invalid min/max all atoms\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.ring_specifier_size() > 0)
  {
    if ( ! _parse_ring_specifier(proto))
      return 0;
  }

  if (proto.ring_system_specifier_size() > 0)
  {
    if (! _parse_ring_system_specifier(proto))
      return 0;
  }

  if (proto.element_hits_needed_size() > 0)
  {
    if (! _parse_element_hits_needed(proto))
      return 0;
  }

  if (proto.has_compress_embeddings())
    _compress_embeddings = proto.compress_embeddings();

  extending_resizable_array<Substructure_Atom *> completed;

  // Form root atoms first.
  for (const auto& atom : proto.query_atom())
  {
    if (IsRootSubstructureAtom(atom))
    {
      Substructure_Atom * r = new Substructure_Atom;

      if (! r->construct_from_proto(atom, completed))
      {
        delete r;
        return 0;
      }

      _root_atoms.add(r);
    }
  }

  // And now non-root atoms.
  for (const auto & atom : proto.query_atom()) {
    if (IsRootSubstructureAtom(atom))
      continue;

    Substructure_Atom* a = new Substructure_Atom();
    if (! a->construct_from_proto(atom, completed)) {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid query atom " << atom.ShortDebugString() << "\n";
      delete a;
      return 0;
    }
  }

  // Now that all atoms are available in the completed array, form bonds.
  // Note that this only handles query_bond specifications, not single_bond, double_bond...
  for (const auto& atom : proto.query_atom()) {
    Substructure_Atom * a = completed[atom.id()];
    if (! a->FormBonds(atom, completed)) {
      cerr << "Single_Substructure_Query::_construct_from_proto:FormBonds failed\n";
      return 0;
    }
    if (!IsRootSubstructureAtom(atom) && a->nbonds() == 0) {
      cerr << "Non root Substructure_Atom not bonded, id " << a->unique_id() << "\n";
      cerr << atom.ShortDebugString() << "\n";
      return 0;
    }
  }

// OK if the root is not built, but let's make sure there are some global
// attributes available in that case

  if (_root_atoms.empty()) {
//  cerr << "Single_Substructure_Query::_construct_from_proto: Root atom not found\n";

//  apr 99. There are some global conditions which don't make sense with no root atom

//  TODO add check for atom types here.
    if (_element_hits_needed.number_elements() || _distance_between_root_atoms.is_set() || 
        _no_matched_atoms_between.number_elements() || _attached_heteroatom_count.is_set() ||
        _ncon.is_set() || _heteroatoms_matched.is_set() || _ring_atoms_matched.is_set())
    {
      cerr << "Single_Substructure_Query::_construct_from_proto: global conditions incompatible with no root atom\n";
      return 0;
    }

//  Should list all the possible match conditions

    const int nat = _compute_attribute_counts();

    if (nat > 0) {
      // Great, got things to match.
    } else if (_natoms.is_set() || _nrings.is_set() || _aromatic_rings.is_set() ||
             _atoms_in_spinach.is_set() || _inter_ring_atoms.is_set() ||
             _net_formal_charge.is_set() || 
             _ring_specification.size() > 0 ||
             _ring_system_specification.size() > 0) {
    } else {
      cerr << "No root atoms and no global attributes specified, rejected\n";
      return 0;
    }
  }

  if (! _construct_environment_from_proto(proto.environment(), completed, _environment)) {
    cerr << "Single_Substructure_Query::_construct_from_proto: environment interpretation failed\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _construct_environment_from_proto(proto.environment_no_match(), completed, _environment_rejections)) {
    cerr << "Single_Substructure_Query::_construct_from_proto: environment interpretation failed\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (environments_can_share_attachment_points >= 0) {
    for (int i = 0; i < _environment.number_elements(); ++i) {
      _environment[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
    }
    for (int i = 0; i < _environment_rejections.number_elements(); ++i) {
      _environment_rejections[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
    }
  }

  if (proto.has_atom_type()) {
    const auto t1 = proto.atom_type();
    const_IWSubstring t2(t1.data(), t1.size());
    _atom_typing = new Atom_Typing_Specification();
    if (! _atom_typing->build(t2)) {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid atom type specification " << t2 << "'\n";
      return 0;
    }
  }

  if (proto.no_matched_atoms_between().size() > 0) {
    // TODO(ianwatson) implement
  }

  if (proto.region().size() > 0) {
    // TODO(ianwatson) impossible
  }

  if (proto.geometric_constraints().size() > 0) {
    for (const auto& constraint : proto.geometric_constraints()) {
      std::unique_ptr<geometric_constraints::SetOfGeometricConstraints> c =
              std::make_unique<geometric_constraints::SetOfGeometricConstraints>();
      if (! c->BuildFromProto(constraint)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid distance constraint " << constraint.ShortDebugString() << '\n';
        return 0;
      }
      if (! _atom_numbers_in_geometric_constraints_ok(c.get())) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid atom numbers in geometric constraint " << constraint.ShortDebugString() << '\n';
        return 0;
      }
      _geometric_constraints.add(c.release());
    }
  }

  if (proto.separated_atoms().size() > 0) {
    for (const auto& separated_atoms : proto.separated_atoms()) {
      std::unique_ptr<SeparatedAtoms> s(new SeparatedAtoms);
      if (! s->Build(separated_atoms)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid separated atoms proto " << separated_atoms.ShortDebugString() << '\n';
        return 0;
      }
      _separated_atoms.add(s.release());
    }
  }

  if (proto.matched_atom_must_be_size() > 0) {
    if (! _construct_matched_atoms_match(proto.matched_atom_must_be())) {
      cerr << "Single_Substructure_Query::_construct_matched_atoms_match:cannot parse matched atoms must be\n";
      return 0;
    }
  }

  for (const auto& inter_ring : proto.inter_ring_region()) {
    std::unique_ptr<InterRingAtoms> s = std::make_unique<InterRingAtoms>();
    if (! s->ConstructFromProto(inter_ring)) {
      cerr << "Single_Substructure_Query::ConstructFromProto:invalid inter ring " << inter_ring.ShortDebugString() << '\n';
      return 0;
    }
    _inter_ring_region.add(s.release());
  }

//if (_root_atoms.number_elements())
//  _adjust_for_internal_consistency();

  return 1;
}

int
Single_Substructure_Query::_atom_numbers_in_geometric_constraints_ok(const geometric_constraints::SetOfGeometricConstraints* constraints) const {
  const resizable_array<int> atom_numbers_present = constraints->AtomNumbersPresent();
  for (int a : atom_numbers_present) {
    const Substructure_Atom* q = query_atom_with_initial_atom_number(a);
    if (q == nullptr) {
      cerr << "Single_Substructure_Query::_atom_numbers_in_geometric_constraints_ok:no query atom " << a << '\n';
      return 0;
    }
  }

  return 1;
}

// Note that this silently does the wrong thing if called multiple times.
int
Single_Substructure_Query::ConstructFromProto(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  assert (ok());

  if (! _construct_from_proto(proto))
  {
    cerr << "SingleSubstructureQuery::ConstructFromProto:cannot build from proto\n";
    cerr << proto.ShortDebugString() << "\n";
    return 0;
  }

// Now that all the atoms have been defined, we can read in any chirality

  for (const auto& chiral : proto.chiral_centre())
  {
    if (! _build_chirality_specification_from_proto(chiral))
    {
      cerr << "Single_Substructure_Query::ConstructFromProto:invalid chirality '" << chiral.ShortDebugString() << "'\n";
      return 0;
    }
  }

  _preferences_present = 0;
  for (const auto* root_atom : _root_atoms) {
    if (root_atom->preferences_present()) {
      _preferences_present = 1;
      break;
    }
  }

  _fragment_ids_present = 0;
  for (const auto * root_atom : _root_atoms)
  {
    if (root_atom->fragment_ids_present())
    {
      _fragment_ids_present = 1;
      break;
    }
  }

  assert (ok());

  (void) min_atoms_in_query();

  for (int i = 0; i < _no_matched_atoms_between.number_elements(); i++)
  {
    const Bond * b = _no_matched_atoms_between[i];
    if (b->a1() >= _min_atoms_in_query || b->a2() >= _min_atoms_in_query)
    {
      cerr << "Single_Substructure_Query::ConstructFromProto: illegal no_matched_atoms_between specifier\n";
      cerr << "There are as few as " << _min_atoms_in_query << " query atoms\n";
      cerr << "No matched atoms between " << i << " specifies atoms " << b->a1() << " and " << b->a2() << '\n';
      cerr << "this is impossible\n";
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::ReadProto(iwstring_data_source & input)
{
  assert (input.ok());

  input.set_ignore_pattern("^#");
  input.set_skip_blank_lines(1);

  IWString contents;
  const_IWSubstring line;
  while (input.next_record(line))
  {
    if (contents.length() > 0)
      contents << "\n";
    contents << line;
    if (line.starts_with(kCloseBrace))
      break;
  }

  const std::string s(contents.rawchars(), contents.length());

  SubstructureSearch::SingleSubstructureQuery proto;

  if (!google::protobuf::TextFormat::ParseFromString(s, &proto))
  {
    cerr << "SingleSubstructureQuery::ReadProto:cannot parse proto data\n";
    cerr << s << '\n';
    return 0;
  }

  return ConstructFromProto(proto);
}

int
Single_Substructure_Query::ReadProto(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadProto(input);
}

int
Single_Substructure_Query::ReadProto(const IWString & fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadProto(input);
}

SubstructureSearch::SubstructureQuery
Substructure_Query::BuildProto() const
{
  SubstructureSearch::SubstructureQuery to_be_returned;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->assign_unique_numbers();
  }

  to_be_returned.set_name(_comment.rawchars(), _comment.length());

  for (const auto* query : *this)
  {
    query->BuildProto(*to_be_returned.add_query());
  }

  if (_each_component_search)
    to_be_returned.set_match_each_component(true);

  return to_be_returned;
}

int
Substructure_Query::WriteProto(IWString & fname) const
{
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars()))
  {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteProto(output);
}

int
Substructure_Query::WriteProto(IWString_and_File_Descriptor& output) const
{
  const SubstructureSearch::SubstructureQuery proto = BuildProto();

  std::string s;
  if (! google::protobuf::TextFormat::PrintToString(proto, &s))
  {
    cerr << "Substructure_Query::WriteProto:PrintToString failed\n";
    return 0;
  }

  return output.write(s.data(), s.length());
}

int
Substructure_Ring_Base::ConstructFromProto(const SubstructureSearch::SubstructureRingBase& proto)
{
  if (proto.has_match_as_match()) {
    _match_as_match_or_rejection = ! proto.match_as_match();
  }

  if (proto.has_all_hits_in_same_fragment()) {
    _all_hits_in_same_fragment = proto.all_hits_in_same_fragment();
  }
    
  if (proto.has_environment()) {
    if (_environment_atom.number_elements()) {
      cerr << "Substructure_Ring_Base::ConstructFromProto:environment already specified '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    const IWString env = proto.environment();
    
    if (! _construct_environment(env)) {
      cerr << "Substructure_Ring_Base::ConstructFromProto:invalid environment '" << env << "'\n";
      return 0;
    }
  }

  if (proto.has_environment_can_match_in_ring_atoms())
    _environment_can_match_in_ring_atoms = proto.environment_can_match_in_ring_atoms();

  if (proto.has_set_global_id()) {
    _set_global_id = proto.set_global_id();
  }

  if (proto.has_ring_extends_to_carbonyl()) {
    _ring_extends_to_carbonyl = proto.ring_extends_to_carbonyl();
  }

  if (!GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ncon, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, heteroatom_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, attached_heteroatom_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, within_ring_unsaturation, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_with_pi_electrons, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fully_saturated_atoms, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, largest_number_of_bonds_shared_with_another_ring, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, strongly_fused_ring_neighbours, 0, no_limit))
    return 0;

  for (const auto& substituent : proto.substituent()) {
    std::unique_ptr<Substituent> s = std::make_unique<Substituent>();
    if (! s->ConstructFromProto(substituent)) {
      cerr << "Substructure_Ring_Base::ConstructFromProto:invalid substituent " << substituent.ShortDebugString() << '\n';
      return 0;
    }
    _substituent.add(s.release());
  }

  return 1;
}


int
Substructure_Ring_Specification::ConstructFromProto(SubstructureSearch::SubstructureRingSpecification const& proto)
{
  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (! Substructure_Ring_Base::ConstructFromProto(proto.base()))
    return 0;

  if (!GETVALUES(proto, ring_size, 3, no_limit))
    return 0;

  if (!GETVALUES(proto, fused, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_aromatic_neighbours, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_non_aromatic_neighbours, 0, no_limit))
    return 0;

  if (proto.has_aromatic())
  {
    if (proto.aromatic())
      _aromatic = 1;
    else
      _aromatic = 0;
  }

  if (proto.has_spiro_fusion_count()) {
    _spiro_fusion_count = proto.spiro_fusion_count();
  }

  return 1;
}

// Ring_size is handled specially, because min_ring_size has a meaning of its own,
// and so cannot be handled by GETVALUES.
int
Substructure_Ring_System_Specification::GetRingSize(const SubstructureSearch::SubstructureRingSystemSpecification& proto) {
  for (int rsize : proto.ring_size()) {
    _ring_sizes.add(rsize);
  }

  if (proto.min_ring_size_count() > 0) {
    _ring_sizes.set_min(proto.min_ring_size_count());
  }
  if (proto.max_ring_size_count() > 0) {
    _ring_sizes.set_max(proto.max_ring_size_count());
  }

  return 1;
}

int
Substructure_Ring_System_Specification::ConstructFromProto(const SubstructureSearch::SubstructureRingSystemSpecification& proto)
{
  if (! Substructure_Ring_Base::ConstructFromProto(proto.base()))
    return 0;

  if (_atoms_with_pi_electrons.is_set() || _fully_saturated_atoms.is_set())
    _need_per_atom_array = 1;

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, rings_in_system, 0, no_limit))
    return 0;

  if (! GetRingSize(proto)) {
    return 0;
  }

  if (!GETVALUES(proto, aromatic_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, non_aromatic_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, degree_of_fusion, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_in_system, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, strongly_fused_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, number_spinach_groups, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, number_non_spinach_groups, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_in_spinach_group, 0, no_limit))
    return 0;

  if (proto.every_group_matches_atoms_in_spinach()) {
    _every_group_matches_atoms_in_spinach_group = true;
  }
  if (proto.every_group_matches_length_of_spinach()) {
    _every_group_matches_length_of_spinach_group = true;
  }

  if (!GETVALUES(proto, length_of_spinach_group, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, distance_to_another_ring, 0, no_limit))
    return 0;

  if (proto.has_ring_systems_extend_across_spiro_fusions()) {
    _ring_systems_extend_across_spiro_fusions = proto.ring_systems_extend_across_spiro_fusions();
  }

  extending_resizable_array<int> ring_sizes_encountered;

  for (const auto& ring_size_requirement : proto.ring_size_requirement()) {
    std::unique_ptr<RingSizeRequirement> r = std::make_unique<RingSizeRequirement>();
    if (! r->ConstructFromProto(ring_size_requirement)) {
      cerr << "Substructure_Ring_System_Specification::ConstructFromProto:invalid ring size\n";
      cerr << ring_size_requirement.ShortDebugString() << '\n';
      return 0;
    }

    // Can only check ring sizes across RingSizeRequirement's if they specify a single value
    if (ring_size_requirement.ring_size_size() != 1) {
      // Cannot check across values.
    } else if (ring_sizes_encountered[ring_size_requirement.ring_size(0)]) {
      cerr << "Substructure_Ring_System_Specification::ConstructFromProto:duplicate ring size specificiation in ring size requirement\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    } else {
      ring_sizes_encountered[ring_size_requirement.ring_size(0)] = 1;
    }

    _ring_size_requirement.add(r.release());
  }

  return 1;
}

RingSizeRequirement::RingSizeRequirement() {
}

int
RingSizeRequirement::ConstructFromProto(const SubstructureSearch::RingSizeRequirement& proto) 
{
  for (uint32_t rsize : proto.ring_size()) {
    _ring_size.add(rsize);
  }

  if (proto.has_min_ring_size()) {
    _ring_size.set_min(proto.min_ring_size());
  }
  if (proto.has_max_ring_size()) {
    _ring_size.set_max(proto.min_ring_size());
  }

  if (!_ring_size.is_set()) {
    cerr << "RingSizeRequirement::ConstructFromProto:no specifications " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (!GETVALUES(proto, count, 0, no_limit)) {
    return 0;
  }

  // A requirement with a ring size, but no count is interpreted as
  // a need for at least one instance of that ring size.
  if (! _count.is_set()) {
    _count.set_min(1);
  }

  return 1;
}

int
Single_Substructure_Query::_parse_element_hits_needed(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  for (const auto& spec : proto.element_hits_needed())
  {
    Elements_Needed * e = new Elements_Needed();
    if (! e->ConstructFromProto(spec))
    {
      cerr << "SingleSubstructureQuery::_parse_element_hits_needed:invalid input " << spec.ShortDebugString() << "\n";
      delete e;
      return 0;
    }

    _element_hits_needed.add(e);
  }

  return 1;
}

int
RequiredMolecularProperties::ParseElementsNeeded(const SubstructureSearch::RequiredMolecularProperties& proto)
{
  for (const auto & spec : proto.elements_needed())
  {
    std::unique_ptr<Elements_Needed> e = std::make_unique<Elements_Needed>();
    if (! e->ConstructFromProto(spec)) {
      cerr << "RequiredMolecularProperties::ParseElementsNeeded:invalid " << spec.ShortDebugString() << "\n";
      return 0;
    }

    _elements_needed << e.release();
  }

  return 1;
}

int
Elements_Needed::ConstructFromProto(const SubstructureSearch::ElementsNeeded& proto)
{
  for (int z : proto.atomic_number()) {
    const Element * e = get_element_from_atomic_number(z);
    if (e == nullptr) {
      cerr << "Elements_Needed::construct_from_proto:no element for " << z << '\n';
      return 0;
    }
    _elements.add(e);
  }
  for (std::string symbol : proto.atomic_symbol()) {
    const Element * e = get_element_from_symbol_no_case_conversion(symbol.data(), symbol.length());
    if (nullptr == e) {
      cerr << "Elements_Needed::construct_from_proto:invalid atomic symbol " << proto.ShortDebugString() << "\n";
      return 0;
    }
    _elements.add(e);
  }

  if (_elements.empty()) {
    cerr << "Elements_Needed::construct_from_proto:no elements\n";
    return 0;
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (! GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  if (proto.multiple_values_operator_size() > 0) {
    if (_elements.number_elements() == 1) {
      cerr << "Elements_Needed::construct_from_proto:operator with just one value " << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  // The simple case of just one element, and no operator.
  if (_elements.number_elements() == 1 && proto.multiple_values_operator_size() == 0) {
    return 1;
  }

  // Default operator when multiple elements present is AND.
  if (_elements.number_elements() > 1 && proto.multiple_values_operator_size() == 0) {
    _operator = Op::AND;
    return 1;
  }

  // Deliberate design decision to limit complexity.
  const auto op = proto.multiple_values_operator(0);
  switch (op) {
    case SubstructureSearch::SS_AND:
      _operator = Op::AND;
      break;
    case SubstructureSearch::SS_OR:
      _operator = Op::OR;
      break;
    case SubstructureSearch::SS_XOR:
      _operator = Op::XOR;
      break;
    default:
      cerr << "Elements_Needed::construct_from_proto:unrecognised operator " << proto.ShortDebugString() << '\n';
      return 0;
  }

  return 1;
}

int
Link_Atom::ConstructFromProto(const SubstructureSearch::LinkAtoms & proto)
{
  if (!proto.has_a1() || !proto.has_a2())
  {
    cerr << "LinkAtoms::ConstructFromProto:incomplete atom spec " << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (proto.distance_size() > 0)
    ;
  else if (proto.has_min_distance() || proto.has_max_distance())
    ;
  else
  {
    cerr << "LinkAtoms::ConstructFromProto:incomplete distance spec " << proto.ShortDebugString() << "\n";
    return 0;
  }

  _a1 = proto.a1();
  _a2 = proto.a2();
  if (_a1 == _a2)
  {
    cerr << "Link_Atom::ConstructFromProto:atoms must be distinct " << proto.ShortDebugString() << "\n";
    return 0;
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, distance, 0, no_limit))
    return 0;

  return 1;
}

int
DownTheBond::ConstructFromProto(const SubstructureSearch::DownTheBond& proto) {
  if (! proto.has_a1() || ! proto.has_a2()) {
    cerr << "DownTheBond::ConstructFromProto:a1 or a2 missing " << proto.ShortDebugString() << '\n';
    return 0;
  }
  _a1 = proto.a1();
  _a2 = proto.a2();

  if (_a1 == _a2) {
    cerr << "DownTheBond::ConstructFromProto:a1 a2 the same " << proto.ShortDebugString() << '\n';
    return 0;
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (! GETVALUES(proto, natoms, 0, no_limit)) {
    return 0;
  }

  return 1;
}

int
DownTheBond::BuildProto(SubstructureSearch::DownTheBond& proto) const {
  proto.set_a1(_a1);
  proto.set_a2(_a2);
  SETPROTOVALUES(proto, natoms, int);

  return 1;
}

#ifdef IMPLEMENT_THIS
TODO
int
Single_Substructure_Query::_add_component(const SubstructureSearch::SubstructureChiralCenter::AtomNumberHydrogenLonePair&& atom_or,
    void (Substructure_Chiral_Centre::*)(atom_number_t) setter,
    resizable_array<int>& numbers_encountered)
{
}
#endif

int
Single_Substructure_Query::_build_chirality_component(const SubstructureSearch::SubstructureChiralCenter& proto,
                                                      void (Substructure_Chiral_Centre::* ptrmfn)(const Substructure_Atom*),
                                                      const SubstructureSearch::AtomNumberHydrogenLonePair& atom_or,
                                                      Substructure_Chiral_Centre& c,
                                                      resizable_array<int>& numbers_encountered)
{
  static constexpr int hydrogen = -9;

  if (atom_or.has_atom_number())
  {
    const int uid = atom_or.atom_number();
    const Substructure_Atom * a = query_atom_with_initial_atom_number(uid);
    if (nullptr == a)
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:no top front atom '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    if (!numbers_encountered.add_if_not_already_present(uid))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:duplicate atom number " << proto.ShortDebugString() << '\n';
      return 0;
    }
    (c.*ptrmfn)(a);
  }
  else  // hydrogen or lone pair 
  { 
    const auto h_or_lp = atom_or.h_or_lp();
    if (h_or_lp == SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN)
      ;
    else if (h_or_lp == SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR) 
    {
      cerr <<"Single_Substructure_Query::_build_chirality_specification_from_proto:lp directive not supported yet " << proto.ShortDebugString() << '\n';
      return 0;
    }
    else
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:unrecognised chirality component " << proto.ShortDebugString() << '\n';
      return 0;
    }

    if (!numbers_encountered.add_if_not_already_present(hydrogen))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:multiple H " << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_build_chirality_specification_from_proto(const SubstructureSearch::SubstructureChiralCenter& proto)
{
  if (! proto.has_center())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:chirality with no center atom " << proto.ShortDebugString() << '\n';
    return 0;
  }
#ifdef NOW_CHECKED_BELOW
  if (! proto.has_top_front() || ! proto.has_top_back() ||
      ! proto.has_left_down() || ! proto.has_right_down())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:all components of a chiral center must be specified\n";
    cerr << proto.ShortDebugString() << '\n';
    return 9;
  }
#endif

  const int uid = proto.center();

  const Substructure_Atom * a = query_atom_with_initial_atom_number(uid);
  if (nullptr == a) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:no centre atom id "
         << uid << " in '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  std::unique_ptr<Substructure_Chiral_Centre> c = std::make_unique<Substructure_Chiral_Centre>();

  c->set_centre(a);

  resizable_array<int> numbers_encountered;   // make sure no duplicates

  numbers_encountered.add(uid);

  int unspecified = 0;

  // Top front.

  auto atom_or = proto.top_front();

  if (! proto.has_top_front()) {
    ++unspecified;
  } else if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_top_front,
                                          proto.top_front(), *c, numbers_encountered)) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid top front " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_top_back()) {
    ++unspecified;
  } else if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_top_back,
                                          proto.top_back(), *c, numbers_encountered)) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid top back " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_left_down()) {
    ++unspecified;
  } else if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_left_down,
                                          proto.left_down(), *c, numbers_encountered)) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid left down back " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_right_down()) {
    ++unspecified;
  } if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_right_down,
                                     proto.right_down(), *c, numbers_encountered)) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid right down back " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (unspecified > 1) {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto: " << unspecified << " attributes unspecified\n";
    return 0;
  }

  for (int i = 0; i < numbers_encountered.number_elements(); ++i) {
    for (int j = i + 1; j < numbers_encountered.number_elements(); ++j) {
      if (numbers_encountered[i] == numbers_encountered[j]) {
        cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:duplicate atom numbers '" << proto.ShortDebugString() << "'\n";
        return 0;
      }
    }
  }

  _chirality.add(c.release());

  return 1;
}

int
Substructure_Query::ConstructFromProto(const SubstructureSearch::SubstructureQuery& proto)
{
  if (proto.has_name()) {
    _comment = proto.name();
  } else if (proto.comment().size() > 0) {
    _comment = proto.comment(0);
  }

  if (proto.query().empty() && proto.query_file().empty()) {
    cerr << "Substructure_Query::ConstructFromProto:no components in query\n";
    return 0;
  }

  for (const auto& query : proto.query())
  {
    std::unique_ptr<Single_Substructure_Query> q = std::make_unique<Single_Substructure_Query>();

    if (! q->ConstructFromProto(query)) {
      cerr << "Substructure_Query::ConstructFromProto:cannot parse component\n";
      cerr << query.ShortDebugString() << '\n';
      return 0;
    }
    add(q.release());
  }

  for (const std::string& fname : proto.query_file()) {
    IWString tmp(fname);
    std::optional<SubstructureSearch::SubstructureQuery> maybe_qry =
       iwmisc::ReadTextProtoCommentsOK<SubstructureSearch::SubstructureQuery>(tmp);
    if (! maybe_qry) {
      cerr << "Substructure_Query::ConstructFromProto:cannot read query from '" << fname << "'\n";
      return 0;
    }
    std::unique_ptr<Single_Substructure_Query> qry = std::make_unique<Single_Substructure_Query>();
    if (maybe_qry->query().empty()) {
      cerr << "Substructure_Query::ConstructFromProto:no query in proto '" << fname << "'\n";
      return 0;
    }
    if (! qry->ConstructFromProto(maybe_qry->query(0))) {
      cerr << "Substructure_Query::ConstructFromProto:cannot parse " << maybe_qry->ShortDebugString() << '\n';
      return 0;
    }
    add(qry.release());
  }

  if (proto.has_match_each_component())
    _each_component_search = proto.match_each_component();

  if (1 == _number_elements) {
    if (0 == proto.logexp_size())  // Great.
      return 1;

    cerr << "Substructure_Query::ConstructFromProto:operator(s) with one component\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (0 == proto.logexp_size())
    return _number_elements;

  if (_each_component_search && proto.logexp_size() > 0) {
    cerr << "Substructure_Query::ConstructFromProto:cannot combine match_each_component with operators\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  _operator.RemoveAllOperators();

  if (! ExtractOperator(proto.logexp(), _number_elements - 1, _operator, IW_LOGEXP_OR, "Substructure_Query::ConstructFromProto"))
  {
    cerr << "Substructure_Query::ConstructFromProto:cannot process operators\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  return _number_elements;
}

int
Substructure_Query::ReadProto(const IWString& fname) {
  IWString tmp(fname);  // we need a non const version
  std::optional<SubstructureSearch::SubstructureQuery> maybe_proto =
    iwmisc::ReadTextProto<SubstructureSearch::SubstructureQuery>(tmp);
  if (! maybe_proto) {
    cerr << "Substructure_Query::ReadProto:cannot read proto from '" << fname << "'\n";
    return 0;
  }

  return ConstructFromProto(*maybe_proto);
}

int
Substructure_Ring_Base::BuildProto(SubstructureSearch::SubstructureRingBase& proto) const {
  if (! _match_as_match_or_rejection) {  // Default is true.
    proto.set_match_as_match(_match_as_match_or_rejection);  // 1
  }

  SETPROTOVALUES(proto, hits_needed, int);  // 2
  SETPROTOVALUES(proto, attached_heteroatom_count, int);  // 5
  SETPROTOVALUES(proto, heteroatom_count, int);  // 8
  SETPROTOVALUES(proto, ncon, int);  // 11
  if (_all_hits_in_same_fragment) {  // 14
    proto.set_all_hits_in_same_fragment(true);
  }
  SETPROTOVALUES(proto, within_ring_unsaturation, int);  // 16
  SETPROTOVALUES(proto, largest_number_of_bonds_shared_with_another_ring, int);  // 19
  SETPROTOVALUES(proto, atoms_with_pi_electrons, int);  // 26
  SETPROTOVALUES(proto, fully_saturated_atoms, int);  // 29
  SETPROTOVALUES(proto, strongly_fused_ring_neighbours, int);  // 32
  // environment as string....
  if (_environment_can_match_in_ring_atoms) {
    proto.set_environment_can_match_in_ring_atoms(true);  // 23
  }
  if (_set_global_id) {
    proto.set_set_global_id(_set_global_id);
  }
  if (_ring_extends_to_carbonyl) {
    proto.set_ring_extends_to_carbonyl(true);
  }

  for (const Substituent* substituent : _substituent) {
    substituent->BuildProto(*proto.mutable_substituent()->Add());
  }

  return 1;
}

int
Substructure_Ring_Specification::BuildProto(SubstructureSearch::SubstructureRingSpecification& proto) const {

  ::Substructure_Ring_Base::BuildProto(*proto.mutable_base());

  SETPROTOVALUES(proto, ring_size, int);  // 2
  if (_aromatic == AROMATIC) {
    proto.set_aromatic(true);   // 5
  } else if (_aromatic == NOT_AROMATIC) {
    proto.set_aromatic(false);   // 5
  }
  SETPROTOVALUES(proto, fused, int);   // 6
  SETPROTOVALUES(proto, fused_aromatic_neighbours, int);   // 9
  SETPROTOVALUES(proto, fused_non_aromatic_neighbours, int);   // 12

  return 1;
}

int
Substructure_Ring_System_Specification::BuildProto(SubstructureSearch::SubstructureRingSystemSpecification& proto) const {
  SETPROTOVALUES(proto, rings_in_system, int);  // 2

  if (_ring_sizes.is_set()) {
    for (int rsize : _ring_sizes) {
      proto.add_ring_size(rsize);
    }
    int m;
    if (_ring_sizes.min(m)) {
      proto.set_min_ring_size_count(m);
    }
    if (_ring_sizes.max(m)) {
      proto.set_max_ring_size_count(m);
    }
  }

  for (const RingSizeRequirement * rsc : _ring_size_requirement) {
    SubstructureSearch::RingSizeRequirement* s = proto.add_ring_size_requirement();
    rsc->BuildProto(*s);
  }
  SETPROTOVALUES(proto, aromatic_ring_count, int);  // 11
  SETPROTOVALUES(proto, non_aromatic_ring_count, int);  // 14
  SETPROTOVALUES(proto, degree_of_fusion, int);  // 17
  SETPROTOVALUES(proto, atoms_in_system, int);  // 20
  SETPROTOVALUES(proto, number_spinach_groups, int);  // 23
  SETPROTOVALUES(proto, number_non_spinach_groups, int);  // 26
  SETPROTOVALUES(proto, atoms_in_spinach_group, int);  // 29
  SETPROTOVALUES(proto, length_of_spinach_group, int);  // 32
  SETPROTOVALUES(proto, distance_to_another_ring, int);  // 35
  SETPROTOVALUES(proto, strongly_fused_ring_count, int);  // 38

  if (_ring_systems_extend_across_spiro_fusions) {
    proto.set_ring_systems_extend_across_spiro_fusions(true);
  }

  if (_every_group_matches_atoms_in_spinach_group) {
    proto.set_every_group_matches_atoms_in_spinach(true);
  }

  if (_every_group_matches_length_of_spinach_group) {
    proto.set_every_group_matches_length_of_spinach(true);
  }

  return 1;
}

int
RingSizeRequirement::BuildProto(SubstructureSearch::RingSizeRequirement & proto) const {
  SETPROTOVALUES(proto, ring_size, int);
  SETPROTOVALUES(proto, count, int);
  return 1;
}

int
Elements_Needed::BuildProto(SubstructureSearch::ElementsNeeded& proto) const {
  for (const Element * e : _elements) {
    const IWString & s = e->symbol();
    const std::string as_string(s.data(), s.length());
    proto.add_atomic_symbol(as_string);
  }
  SETPROTOVALUES(proto, hits_needed, int);
  switch (_operator) {
    case Op::NOT_SET:
      break;
    case Op::AND:
      proto.add_multiple_values_operator(SubstructureSearch::SS_AND);
      break;
    case Op::OR:
      proto.add_multiple_values_operator(SubstructureSearch::SS_OR);
      break;
    case Op::XOR:
      proto.add_multiple_values_operator(SubstructureSearch::SS_XOR);
      break;
    default:
      cerr << "Elements_Needed::BuildProto:unrecognised operator " << static_cast<int>(_operator) << '\n';
      return 0;
  }
  return 1;
}

int
Link_Atom::BuildProto(SubstructureSearch::LinkAtoms& proto) const {
  proto.set_a1(_a1);
  proto.set_a2(_a2);
  SETPROTOVALUES(proto, distance, int);
  SetProtoValues(_distance, "distance", proto);
  return 1;
}

int
Substructure_Chiral_Centre::BuildProto(SubstructureSearch::SubstructureChiralCenter& proto) const {
  if (_numeric == nullptr) {
    proto.set_center(_centre->initial_atom_number());
  } else {
    proto.set_center(_numeric->a());
  }

  int unspecified = 0;

  int tf;
  if (_numeric != nullptr) {
    tf = _numeric->top_front();
  } else if (_top_front != nullptr) {
    tf = _top_front->initial_atom_number();
  } else {
    tf = INVALID_ATOM_NUMBER;
    ++unspecified;
  }

  if (tf == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_top_front()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (tf == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_top_front()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_top_front()->set_atom_number(tf);
  }

  int tb;
  if (_numeric != nullptr) {
    tb = _numeric->top_back();
  } else if (_top_back != nullptr) {
    tb = _top_back->initial_atom_number();
  } else {
    tb = INVALID_ATOM_NUMBER;
    ++unspecified;
  }

  if (tb == INVALID_ATOM_NUMBER) {
  } else if (tb == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_top_back()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (tb == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_top_back()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_top_back()->set_atom_number(tb);
  }

  int ld;
  if (_numeric != nullptr) {
    ld = _numeric->left_down();
  } else if (_left_down != nullptr) {
    ld = _left_down->initial_atom_number();
  } else {
    ld = INVALID_ATOM_NUMBER;
    ++unspecified;
  }

  if (ld == INVALID_ATOM_NUMBER) {
  } else if (ld == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_left_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (ld == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_left_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_left_down()->set_atom_number(ld);
  }

  int rd;
  if (_numeric != nullptr) {
    rd = _numeric->right_down();
  } else if (_right_down != nullptr) {
    rd = _right_down->initial_atom_number();
  } else {
    rd = INVALID_ATOM_NUMBER;
    ++unspecified;
  }

  if (rd == INVALID_ATOM_NUMBER) {
  } else if (rd == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_right_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (rd == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_right_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_right_down()->set_atom_number(rd);
  }

  if (unspecified > 1) {
    cerr << "SubstructureChiralCenter::BuildProto: " << unspecified << " unspecified aspects\n";
  }

  return 1;
}

int
SeparatedAtoms::BuildProto(SubstructureSearch::SeparatedAtoms& proto) const {
  proto.set_a1(_a1);
  proto.set_a2(_a2);
  SetProtoValues(_separation, "bonds_between", proto);
  return 1;
}

int
Substituent::ConstructFromProto(const SubstructureSearch::Substituent& proto) {
  if (!GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;
  if (!GETVALUES(proto, natoms, 1, no_limit))
    return 0;
  if (!GETVALUES(proto, nrings, 0, no_limit))
    return 0;
  if (!GETVALUES(proto, length, 1, no_limit))
    return 0;

  if (proto.has_set_global_id()) {
    _set_global_id = proto.set_global_id();
  }

  // This code copied to InterRingAtoms.
  // Could have set up a template function, but not worth the bother.
  for (const std::string& smarts : proto.required_smarts()) {
    const IWString smt(smarts.data(), smarts.size());
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->create_from_smarts(smt)) {
      cerr << "Substituent::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      return 0;
    }
    _required.add(qry.release());
    _required_smarts.add(new IWString(smarts));
  }

  for (const std::string& smarts : proto.disqualifying_smarts()) {
    const IWString smt(smarts.data(), smarts.size());
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->create_from_smarts(smt)) {
      cerr << "Substituent::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      return 0;
    }
    _disqualifying.add(qry.release());
    _disqualifying_smarts.add(new IWString(smarts));
  }

  return 1;
}

int
Substituent::BuildProto(SubstructureSearch::Substituent& proto) const {
  SetProtoValues(_hits_needed, "hits_needed", proto);
  SetProtoValues(_natoms, "natoms", proto);
  SetProtoValues(_nrings, "nrings", proto);
  SetProtoValues(_length, "length", proto);

  if (_set_global_id > 0) {
    proto.set_set_global_id(_set_global_id);
  }

  for (const IWString* smt : _required_smarts) {
    proto.add_required_smarts(smt->AsString());
  }

  for (const IWString* smt : _disqualifying_smarts) {
    proto.add_disqualifying_smarts(smt->AsString());
  }

  cerr << "Substituent::BuildProto:implement this something\n";
  return 0;
}

int
InterRingAtoms::ConstructFromProto(const SubstructureSearch::InterRingAtoms& proto) {
  if (!GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;
  if (!GETVALUES(proto, natoms, 1, no_limit))
    return 0;
  if (!GETVALUES(proto, ring_connections, 2, no_limit))
    return 0;
  if (!GETVALUES(proto, length, 2, no_limit))
    return 0;

  if (proto.has_set_global_id()) {
    _set_global_id = proto.set_global_id();
  }

  for (const auto length : proto.required_length()) {
    _required_length.add(length);
  }

  if (_required_length.size() > 1) {
    Int_Comparator_Larger icl;
    _required_length.iwqsort(icl);
  }

  // Code copied from Substituent.
  for (const std::string& smarts : proto.required_smarts()) {
    const IWString smt(smarts.data(), smarts.size());
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->create_from_smarts(smt)) {
      cerr << "InterRingAtoms::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      return 0;
    }
    _required.add(qry.release());
    _required_smarts.add(new IWString(smarts));
  }

  for (const std::string& smarts : proto.disqualifying_smarts()) {
    const IWString smt(smarts.data(), smarts.size());
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->create_from_smarts(smt)) {
      cerr << "InterRingAtoms::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      return 0;
    }
    _disqualifying.add(qry.release());
    _disqualifying_smarts.add(new IWString(smarts));
  }

  if (!proto.defined_by_ring_global_id().empty()) {
    if (proto.defined_by_ring_global_id_size() < 2) {
      cerr << "InterRingAtoms::ConstructFromProto:must have at least 2 defined_by_ring_global_id " << proto.ShortDebugString() << '\n';
      return 0;
    }
    for (auto def : proto.defined_by_ring_global_id()) {
      if (! _defined_by_ring_global_id.add_if_not_already_present(def)) {
        cerr << "InterRingAtoms::ConstructFromProto:duplicate defined_by_ring_global_id " << proto.ShortDebugString() << '\n';
        return 0;
      }
    }
  }

  return 1;
}

int
InterRingAtoms::BuildProto(SubstructureSearch::InterRingAtoms& proto) const {
  SetProtoValues(_hits_needed, "hits_needed", proto);
  SetProtoValues(_natoms, "natoms", proto);
  SetProtoValues(_ring_connections, "ring_connections", proto);
  SetProtoValues(_length, "length", proto);

  if (_set_global_id > 0) {
    proto.set_set_global_id(_set_global_id);
  }

  for (int length : _required_length) {
    proto.add_required_length(length);
  }

  for (const IWString* smt : _required_smarts) {
    proto.add_required_smarts(smt->AsString());
  }

  for (const IWString* smt : _disqualifying_smarts) {
    proto.add_disqualifying_smarts(smt->AsString());
  }

  for (int& g : _defined_by_ring_global_id) {
    proto.add_defined_by_ring_global_id(g);
  }

  cerr << "InterRingAtoms::BuildProto:implement this something\n";
  return 0;
}

int
MatchedAtomMatch::ConstructFromProto(const SubstructureSearch::MatchedAtomMatch& proto) {
  if (proto.atom_size() == 0) {
    cerr << "MatchedAtomMatch::ConstructFromProto:no atom " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.smarts_size() == 0) {
    cerr << "MatchedAtomMatch::ConstructFromProto:no smarts " << proto.ShortDebugString() << '\n';
    return 0;
  }

  for (int atom : proto.atom()) {
    _atoms << atom;
  }

  for (const std::string& smarts : proto.smarts()) {
    if (smarts.empty()) {
      cerr << "MatchedAtomMatch::ConstructFromProto:skipping empty smarts\n";
      continue;
    }
    bool match_as_match;
    const_IWSubstring s;
    if (smarts[0] == '!') {
      match_as_match = false;
      s = smarts;
      s += 1;
    } else {
      match_as_match = true;
      s = smarts;
    }

    std::unique_ptr<Single_Substructure_Query> a = std::make_unique<Single_Substructure_Query>();
    if (! a->create_from_smarts(s)) {
      cerr << "MatchedAtomMatch::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }

    if (match_as_match) {
      _positive_matches << a.release();
    } else {
      _negative_matches << a.release();
    }
  }

  return 1;
}

int
MatchedAtomMatch::BuildProto(SubstructureSearch::MatchedAtomMatch& proto) const {
  for (int atom : _atoms) {
    proto.add_atom(atom);
  }

  cerr << "MatchedAtomMatch::BuildProto:cannot emit queries\n";

  return 1;
}

int
Substructure_Atom_Environment::BuildProto(SubstructureSearch::SubstructureAtom& proto) const {
  cerr << "SubstructureAtom::BuildProto:has " << _number_elements << " components\n";
  if (_number_elements == 0) {
    return 1;
  }

  SubstructureSearch::SubstructureAtomEnvironment* env = proto.add_environment();
  for (Substructure_Atom* a : *this) {
    SubstructureSearch::SubstructureAtom* in_proto = env->add_substructure_atom();
    a->BuildProto(*in_proto);
  }

  // If just one component, no operator.
  if (_number_elements == 1) {
    return 1;
  }

  const std::string op = _operator.ToString().AsString();
  cerr << "Encoding logexp as '" << op << "'\n";
  env->set_op(op);

  return 1;
}

int
Substructure_Atom_Environment::BuildFromProto(const SubstructureSearch::SubstructureAtomEnvironment& proto) {
  cerr << "SubstructureAtomEnvironment::BuildFromProto:implement this\n";
  for (const SubstructureSearch::SubstructureAtom& atom_proto : proto.substructure_atom()) {
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();
    extending_resizable_array<Substructure_Atom *> completed;
    if (! a->construct_from_proto(atom_proto, completed)) {
      cerr << "Substructure_Atom::BuildFromProto:cannot parse atom " << atom_proto.ShortDebugString() << '\n';
      return 0;
    }
    this->add(a.release());
  }

  // Cannot be empty.
  if (_number_elements == 0) {
    cerr << "Substructure_Atom_Environment::BuildFromProto:no components read " << proto.ShortDebugString() << '\n';
    return 0;
  }

  // If there is just one component, no operator needed.
  if (! proto.has_op()) {
    if (_number_elements == 1) {
      return 1;
    }
    cerr << "Substructure_Atom_Environment::BuildFromProto:one component, but operator " << proto.op() << "'\n";
    return 0;
  }

  if (! proto.has_op() || proto.op().empty()) {
    cerr << "Substructure_Atom_Environment::BuildFromProto:got " << _number_elements << " components, but no op " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _operator.BuildFromString(proto.op())) {
    cerr << "Substructure_Atom_Environment::BuildFromProto:cannot parse operator '" << proto.op() << "'\n";
    return 0;
  }

  return 1;
}

int
RequiredBond::ConstructFromProto(const SubstructureSearch::RequiredBond& proto) {
  if (proto.has_atomic_number_1()) {
    _atomic_number_1 = proto.atomic_number_1();
  } else if (proto.has_atomic_symbol_1()) {
    const IWString s = proto.atomic_symbol_1();
    isotope_t notused;
    const Element* e = get_element_from_symbol(s, notused);
    if (e == nullptr) {
      cerr << "RequiredBond::ConstructFromProto:unrecognised symbol '" << s << "'\n";
      return 0;
    }
    _atomic_number_1 = e->atomic_number();
  } else {
    cerr << "RequiredBond::ConstructFromProto:must specify first atom type\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_btype()) {
    if (proto.btype() == SubstructureSearch::SS_SINGLE_BOND) {
      _btype = SINGLE_BOND;
    } else if (proto.btype() == SubstructureSearch::SS_DOUBLE_BOND) {
      _btype = DOUBLE_BOND;
    } else if (proto.btype() == SubstructureSearch::SS_TRIPLE_BOND) {
      _btype = TRIPLE_BOND;
    } else {
      cerr << "RequiredBond::ConstructFromProto:Unrecognised bond type\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  } else {
    _btype = SINGLE_BOND;
  }

  if (proto.has_atomic_number_2()) {
    _atomic_number_2 = proto.atomic_number_2();
  } else if (proto.has_atomic_symbol_2()) {
    const IWString s = proto.atomic_symbol_2();
    isotope_t notused;
    const Element* e = get_element_from_symbol(s, notused);
    if (e == nullptr) {
      cerr << "RequiredBond::ConstructFromProto:unrecognised symbol '2 " << s << "'\n";
      return 0;
    }
    _atomic_number_2 = e->atomic_number();
  } else {
    cerr << "RequiredBond::ConstructFromProto:must specify first atom type\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.has_min_count()) {
    _min_count = proto.min_count();
    if (_min_count <= 0) {
      cerr << "RequiredBond::ConstructFromProto:invalid min count\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
RequiredBond::BuildProto(SubstructureSearch::RequiredBond& proto) const {
  proto.set_atomic_number_1(_atomic_number_1);
  if (_btype == SINGLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_SINGLE_BOND);
  } else if (_btype == DOUBLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_DOUBLE_BOND);
  } else if (_btype == TRIPLE_BOND) {
    proto.set_btype(SubstructureSearch::SS_TRIPLE_BOND);
  } else {
    cerr << "RequiredBond::BuildProto:what kind of bond do I have " << _btype << '\n';
    return 0;
  }
  proto.set_atomic_number_2(_atomic_number_2);

  return 1;
}


int
Region::ConstructFromProto(const SubstructureSearch::Region& proto) {
  for (auto a : proto.atom()) {
    _matched_atom.add_if_not_already_present(a);
  }

  if (_matched_atom.size() < 2) {
    cerr << "Region::ConstructFromProto:too few atoms specified\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (_matched_atom.size() > 2) {
    cerr << "Region::ConstructFromProto:regions other than 2 atoms not implemented, see Ian\n";
    return 0;
  }

  if (! GETVALUES(proto, natoms, 0, no_limit)) {
    return 0;
  }

  if (! GETVALUES(proto, nrings, 0, no_limit)) {
    return 0;
  }

  if (! GETVALUES(proto, atoms_not_on_shortest_path, 0, no_limit)) {
    return 0;
  }

  if (_natoms.is_set()) {
  } else if (_nrings.is_set()) {
  } else if (_atoms_not_on_shortest_path.is_set()) {
  } else {
    cerr << "Region::ConstructFromProto:nothing specified\n";
    return 0;
  }

  return 1;
}

int
Region::BuildProto(SubstructureSearch::Region& proto) const {
  for (int a : _matched_atom) {
    proto.add_atom(a);
  }

  if (_natoms.is_set()) {
    SetProtoValues(_natoms, "natoms", proto);
  }

  if (_nrings.is_set()) {
    SetProtoValues(_natoms, "nrings", proto);
  }

  if (_atoms_not_on_shortest_path.is_set()) {
    SetProtoValues(_atoms_not_on_shortest_path, "atoms_not_on_shortest_path", proto);
  }

  return 1;
}

int
RequiredMolecularProperties::ConstructFromProto(const SubstructureSearch::RequiredMolecularProperties& proto) {
  _attributes_specified = 0;

  MATCHER_FROM_PROTO(proto, natoms, uint32_t, _natoms);
  MATCHER_FROM_PROTO(proto, nrings, uint32_t, _nrings);
  MATCHER_FROM_PROTO(proto, heteroatoms_in_molecule, uint32_t, _heteroatoms_in_molecule);
  MATCHER_FROM_PROTO(proto, fused_rings, uint32_t, _fused_rings);
  MATCHER_FROM_PROTO(proto, strongly_fused_rings, uint32_t, _strongly_fused_rings);
  MATCHER_FROM_PROTO(proto, isolated_rings, uint32_t, _isolated_rings);
  MATCHER_FROM_PROTO(proto, ring_systems, uint32_t, _ring_systems);
  MATCHER_FROM_PROTO(proto, aromatic_rings, uint32_t, _aromatic_rings);
  MATCHER_FROM_PROTO(proto, aromatic_atoms, uint32_t, _aromatic_atoms);
  MATCHER_FROM_PROTO(proto, non_aromatic_rings, uint32_t, _non_aromatic_rings);
  MATCHER_FROM_PROTO(proto, number_isotopic_atoms, uint32_t, _number_isotopic_atoms);
  MATCHER_FROM_PROTO(proto, number_fragments, uint32_t, _number_fragments);
  MATCHER_FROM_PROTO(proto, atoms_in_spinach, uint32_t, _atoms_in_spinach);
  MATCHER_FROM_PROTO(proto, inter_ring_atoms, uint32_t, _inter_ring_atoms);
  MATCHER_FROM_PROTO(proto, net_formal_charge, uint32_t, _net_formal_charge);

  _any_net_formal_charge = proto.any_net_formal_charge();

  if (proto.elements_needed_size() > 0) {
    if (! ParseElementsNeeded(proto)) {
      return 0;
    }
  }

  DiscernAttributesSpecified();

  return 1;
}

int
RequiredMolecularProperties::BuildProto(SubstructureSearch::RequiredMolecularProperties& proto) const {
  PROTO_FROM_MATCHER(_natoms, natoms, int, proto);
  PROTO_FROM_MATCHER(_nrings, nrings, int, proto);
  PROTO_FROM_MATCHER(_heteroatoms_in_molecule, heteroatoms_in_molecule, int, proto);
  PROTO_FROM_MATCHER(_fused_rings, fused_rings, int, proto);
  PROTO_FROM_MATCHER(_strongly_fused_rings, strongly_fused_rings, int, proto);
  PROTO_FROM_MATCHER(_isolated_rings, isolated_rings, int, proto);
  PROTO_FROM_MATCHER(_ring_systems, ring_systems, int, proto);
  PROTO_FROM_MATCHER(_aromatic_rings, aromatic_rings, int, proto);
  PROTO_FROM_MATCHER(_aromatic_atoms, aromatic_atoms, int, proto);
  PROTO_FROM_MATCHER(_non_aromatic_rings, non_aromatic_rings, int, proto);
  PROTO_FROM_MATCHER(_number_isotopic_atoms, number_isotopic_atoms, int, proto);
  PROTO_FROM_MATCHER(_number_fragments, number_fragments, int, proto);
  PROTO_FROM_MATCHER(_atoms_in_spinach, atoms_in_spinach, int, proto);
  PROTO_FROM_MATCHER(_inter_ring_atoms, inter_ring_atoms, int, proto);
  PROTO_FROM_MATCHER(_net_formal_charge, net_formal_charge, int, proto);

  if (_any_net_formal_charge) {
    proto.set_any_net_formal_charge(true);
  }

  for (const auto * ele : _elements_needed) {
    ele->BuildProto(*proto.add_elements_needed());
  }

  return 1;
}

#define INCREMENT_IF_SET(matcher, target) { \
  if (matcher.is_set()) { \
    ++target; \
  }\
}

int
RequiredMolecularProperties::DiscernAttributesSpecified() {
  _attributes_specified = 0;
  INCREMENT_IF_SET(_natoms, _attributes_specified);
  INCREMENT_IF_SET(_nrings, _attributes_specified);
  INCREMENT_IF_SET(_heteroatoms_in_molecule, _attributes_specified);
  INCREMENT_IF_SET(_fused_rings, _attributes_specified);
  INCREMENT_IF_SET(_strongly_fused_rings, _attributes_specified);
  INCREMENT_IF_SET(_isolated_rings, _attributes_specified);
  INCREMENT_IF_SET(_ring_systems, _attributes_specified);
  INCREMENT_IF_SET(_aromatic_rings, _attributes_specified);
  INCREMENT_IF_SET(_aromatic_atoms, _attributes_specified);
  INCREMENT_IF_SET(_non_aromatic_rings, _attributes_specified);
  INCREMENT_IF_SET(_number_isotopic_atoms, _attributes_specified);
  INCREMENT_IF_SET(_number_fragments, _attributes_specified);
  INCREMENT_IF_SET(_atoms_in_spinach, _attributes_specified);
  INCREMENT_IF_SET(_inter_ring_atoms, _attributes_specified);
  INCREMENT_IF_SET(_net_formal_charge, _attributes_specified);
  if (_elements_needed.size() > 0) {
    ++_attributes_specified;
  }
  if (_any_net_formal_charge) {
    ++_attributes_specified;
  }

  return _attributes_specified;
}

int
NearbyAtoms::ConstructFromProto(const SubstructureSearch::NearbyAtoms& proto) {
#ifdef DEBUG_NEARBY_ATOMS_CONSTRUCT_FROM_PROTO
  cerr << "NearbyAtoms::ConstructFromProto building from\n";
  cerr << proto.ShortDebugString() << '\n';
#endif

  int nquery = proto.smarts_size() + proto.query_size() + proto.query_file_size();
  if (nquery != 1) {
    cerr << "NearbyAtoms::ConstructFromProto:must specify exactly 1 query source " << nquery << " invalid\n";
    return 0;
  }

  for (const std::string& smt : proto.smarts()) {
    std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
    if (! q->create_from_smarts(smt)) {
      cerr << "NearbyAtoms::ConstructFromProto:invalid smarts '" << smt << "'\n";
      return 0;
    }
    _query << q.release();
  }

  for (const std::string& query_string : proto.query()) {
    SubstructureSearch::SubstructureQuery query_proto;
    if (!google::protobuf::TextFormat::ParseFromString(query_string, &query_proto)) {
      cerr << "NearbyAtoms::ConstructFromProto:cannot parse proto data\n";
      cerr << query_string << '\n';
      return 0;
    }

    std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
    if (! q->ConstructFromProto(query_proto)) {
      cerr << "NearbyAtoms::ConstructFromProto:invalid query '" << query_string << "'\n";
      return 0;
    }
    _query << q.release();
  }

  // We do not make allowance for directory path here...
  for (const std::string& fname : proto.query_file()) {
    IWString tmp(fname);
    std::optional<SubstructureSearch::SubstructureQuery> maybe_proto = 
      iwmisc::ReadTextProto<SubstructureSearch::SubstructureQuery>(tmp);
    if (! maybe_proto) {
      cerr << "NearbyAtoms::ConstructFromProto:cannot read '" << fname << "'\n";
      return 0;
    }
    std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
    if (! q->ConstructFromProto(*maybe_proto)) {
      cerr << "NearbyAtoms::ConstructFromProto:Invalid query '" << fname << "'\n";
      return 0;
    }
    _query << q.release();
  }

  if (proto.has_can_overlap_matched_atoms()) {
    _can_overlap_matched_atoms = proto.can_overlap_matched_atoms();
  }

  if (!GETVALUES(proto, hits_needed, 0, no_limit)) {
    return 0;
  }

  if (! _hits_needed.is_set()) {
    _hits_needed.set_min(1);
  }

  for (uint32_t a : proto.matched_atom()) {
    _matched_atom << a;
  }

  if (! GETVALUES(proto, bonds_between, 0, no_limit)) {
    return 0;
  }

  if (proto.has_rejection()) {
    _rejection = proto.rejection();
  }

  return 1;
}

int
NearbyAtoms::BuildProto(SubstructureSearch::NearbyAtoms& proto) const {
  cerr << "NearbyAtoms::BuildProto:queries not implemented\n";

  SetProtoValues(_hits_needed, "hits_needed", proto);
  SetProtoValues(_hits_needed, "bonds_between", proto);

  for (const uint32_t a : _matched_atom) {
    proto.add_matched_atom(a);
  }

  if (_can_overlap_matched_atoms) {
    proto.set_can_overlap_matched_atoms(_can_overlap_matched_atoms);
  }

  if (_rejection) {
    proto.set_rejection(true);
  }

  return 1;
}

namespace iwsubstructure {
constexpr char kCloseBrace = '}';

std::optional<std::string>
GetNextQueryTextProto(iwstring_data_source& input) {
  std::string result;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (result.empty()) {
      if (! buffer.starts_with("query")) {
        cerr << "GetNextQueryTextProto:First line in textproto must be 'query' '" << buffer << "' invalid\n";
        return std::nullopt;
      }
    }

    result.append(buffer.data(), buffer.length());
    result.push_back('\n');
    if (buffer == kCloseBrace) {
      return result;
    }
  }

  if (result.empty()) {
    return std::nullopt;
  }

  return result;
}


}  // namespace iwsubstructure
