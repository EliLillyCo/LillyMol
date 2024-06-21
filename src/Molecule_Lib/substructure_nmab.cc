#include <algorithm>
#include <iostream>
#include <memory>

// Processing associated with no matched atoms between.

#include "Foundational/iwmisc/misc.h"

#include "mdl_molecule.h"
#include "misc2.h"
#include "parse_smarts_tmp.h"
#include "path.h"
#include "substructure.h"
#include "target.h"

using std::cerr;

static constexpr char kOpenBrace = '{';
static constexpr char kCloseBrace = '}';

namespace nmab
{
void
fill_matched_atoms_array(const Query_Atoms_Matched& matched_atoms, int* matched)
{
  for (const Substructure_Atom* a : matched_atoms) {
    if (!a->include_in_embedding()) {
      continue;
    }

    matched[a->atom_number_matched()] = 1;
  }

  return;
}

int
AtomsInRegion(const Molecule& m, atom_number_t zatom, atom_number_t prev, int* visited,
              atom_number_t destination)
{
  visited[zatom] = 1;
  int rc = 1;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    // Checking visited might also cover this case.
    if (j == prev) {
      continue;
    }

    if (j == destination) {
      ++visited[destination];
      continue;
    }

    if (visited[j]) {
      continue;
    }

    rc += AtomsInRegion(m, j, zatom, visited, destination);
  }

  return rc;
}

// Given two atoms, how many atoms are in between them?
// This will fail if `a1` and `a2` are in the same ring system.
// `visited` is working storage, and `region` will be populated with
// the matched atoms if we find a region.
int
AtomsInRegion(const Molecule& m, atom_number_t a1, atom_number_t a2,
              int* visited,
              int* region) {
  // Make sure a1 has the lowest number of connections.
  if (m.ncon(a1) > m.ncon(a2)) {
    std::swap(a1, a2);
  }

  const int matoms = m.natoms();

  const Atom& a = m.atom(a1);

  // For each atom attached to a1, try to find a2. There should be one non-zero result.
  int atoms_in_region = 0;

  for (const Bond* b : a) {
    atom_number_t j = b->other(a1);
    // If a1 and a2 are bonded, zero atoms in region.
    if (j == a2) {
      return 0;
    }

    std::fill_n(visited, matoms, 0);
    visited[a1] = 1;

    int tmp = AtomsInRegion(m, j, a1, visited, a2);
    if (! visited[a2]) {
      continue;
    }

    // There can only be one region. If we have found a prevous region,
    // that is a failure.
    if (atoms_in_region > 0) {
      return 0;
    }

    atoms_in_region = tmp;
    std::copy_n(visited, matoms, region);
  }

  return atoms_in_region;
}

}  // namespace nmab

// #define DEBUG_NO_MATCHED_ATOM_PATH

/*
  Can we find a path with no matched atoms between two atoms.
  Note that the path must be a shortest path between the two
  atoms. Note that such a path will not be unique.
*/

static int
no_matched_atom_path(const Atom** atoms, atom_number_t destination,
                     int* already_covered_this_path, const int* matched,
                     atom_number_t my_atom, const int current_path_length,
                     const int needed_path_length)
{
  // assert (0 == matched[my_atom]);   NOT TRUE, the first atom will be on the matched
  // atom list

  assert(0 == already_covered_this_path[my_atom]);

  already_covered_this_path[my_atom] = 1;  // put ourselves on the path

  const Atom* a = atoms[my_atom];

#ifdef DEBUG_NO_MATCHED_ATOM_PATH
  cerr << "no_matched_atom_path continues with atom " << my_atom << " type "
       << a->atomic_symbol() << " ncon " << a->ncon() << " destination " << destination
       << '\n';
#endif

  for (const Bond* b : *a) {
    const atom_number_t j = b->other(my_atom);

    if (destination == j &&
        current_path_length + 1 ==
            needed_path_length) {  // bingo. Do this test first as destination will be set
                                   // in matched[]
      return 1;
    }

    if (already_covered_this_path[j]) {
      continue;
    }

    if (matched[j]) {
      continue;
    }

    if (1 == atoms[j]->ncon()) {  // terminal atom, no sense following it
      continue;
    }

    if (no_matched_atom_path(atoms, destination, already_covered_this_path, matched, j,
                             current_path_length + 1, needed_path_length)) {
      return 1;
    }
  }

  already_covered_this_path[my_atom] = 0;  // take ourselves off the path

  return 0;
}

//#define DEBUG_NMAB_SATISFIED

int
Single_Substructure_Query::_nmab_satisfied(Molecule_to_Match& target,
                                           Query_Atoms_Matched& matched_atoms) const
{
#ifdef DEBUG_NMAB_SATISFIED
  cerr << "XJXJXJ checking " << _nmab.size() << " nmab\n";
#endif

  Molecule* m = matched_atoms[0]->current_hold_atom()->m();

  const int matoms = m->natoms();

  int* tmp = new_int(matoms + matoms);
  std::unique_ptr<int[]> free_tmp(tmp);

  nmab::fill_matched_atoms_array(matched_atoms, tmp);

#ifdef DEBUG_NMAB_SATISFIED
  cerr << " matched atoms";
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i]) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
  cerr << "Checking " << _nmab.size() << " nmab conditions\n";
#endif

  if (_no_matched_atoms_between_exhaustive) {
    for (No_Matched_Atoms_Between* nmab : _nmab) {
      if (!nmab->MatchesExhaustive(target, matched_atoms, tmp, tmp + matoms)) {
        return 0;
      }
    }
  } else {
    for (No_Matched_Atoms_Between* nmab : _nmab) {
      if (!nmab->Matches(target, matched_atoms, tmp, tmp + matoms)) {
        return 0;
      }
    }
  }

  return 1;
}

int
Single_Substructure_Query::_no_matched_atoms_between_satisfied(
    Query_Atoms_Matched& matched_atoms) const
{
#ifdef DEBUG_NO_MATCHED_ATOMS_BETWEEN_SATISFIED
  cerr << "Testing " << _nmab.size() << " matched atoms between constraints\n";
#endif

  Molecule* m = matched_atoms[0]->current_hold_atom()->m();

  const int matoms = m->natoms();

  int* matched = new_int(matoms);
  std::unique_ptr<int[]> free_matched(matched);

  nmab::fill_matched_atoms_array(matched_atoms, matched);

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);
  m->atoms(atoms);

  int* tmp = new_int(matoms);  // used in constructing the path
  std::unique_ptr<int[]> free_tmp(tmp);

  for (const Bond* b : _no_matched_atoms_between) {
    //  the numbers in the Bond object refer to query atoms. Need to convert these to atom
    //  numbers in the molecule

    //  const Substructure_Atom * qa1 = matched_atoms[b->a1()];
    //  const Substructure_Atom * qa2 = matched_atoms[b->a2()];

    assert(matched_atoms[b->a1()]->is_matched());
    assert(matched_atoms[b->a2()]->is_matched());

    const atom_number_t a1 = matched_atoms[b->a1()]->current_hold_atom()->atom_number();
    const atom_number_t a2 = matched_atoms[b->a2()]->current_hold_atom()->atom_number();

#ifdef DEBUG_NO_MATCHED_ATOMS_BETWEEN_SATISFIED
    cerr << "Testing constraint between matched atoms " << b->a1() << " and " << b->a2()
         << '\n';
    cerr << "Corresponding to atoms " << a1 << " and " << a2 << " in the molecule\n";
#endif

    if (m->fragment_membership(a1) != m->fragment_membership(a2)) {
      return 1;
    }

    if (!no_matched_atom_path(atoms, a1, tmp, matched, a2, 0, m->bonds_between(a1, a2))) {
      return 0;
    }
  }

  return 1;  // all constraints satisfied
}

int
Single_Substructure_Query::_distance_between_root_atoms_satisfied(
    Query_Atoms_Matched& matched_atoms) const
{
  assert(_distance_between_root_atoms.is_set());

  int nr = _root_atoms.number_elements();

  assert(nr > 1);  // doesn't make sense otherwise

  Molecule* m = matched_atoms[0]->current_hold_atom()->m();

  for (int i = 0; i < nr; i++) {
    const Substructure_Atom* ri = _root_atoms[i];
    assert(ri->is_matched());

    atom_number_t ai = ri->current_hold_atom()->atom_number();

    for (int j = i + 1; j < nr; j++) {
      const Substructure_Atom* rj = _root_atoms[j];
      assert(rj->is_matched());

      atom_number_t aj = rj->current_hold_atom()->atom_number();
      if (m->fragment_membership(ai) != m->fragment_membership(aj)) {
        continue;
      }

      int d = m->bonds_between(ai, aj);
      if (!_distance_between_root_atoms.matches(d)) {
        return 0;
      }
    }
  }

  return 1;
}

Link_Atom::Link_Atom()
{
  _a1 = INVALID_ATOM_NUMBER;
  _bt = INVALID_BOND_TYPE;
  _a2 = INVALID_ATOM_NUMBER;

  _bond_topology = -1;

  _mdl_atom_data = nullptr;

  return;
}

Link_Atom::~Link_Atom()
{
  if (_mdl_atom_data != nullptr) {
    delete _mdl_atom_data;
  }
}

Link_Atom::Link_Atom(const Link_Atom& rhs)
{
  _a1 = rhs._a1;
  _bt = rhs._bt;
  _a2 = rhs._a2;

  _bond_topology = rhs._bond_topology;

  _symbol = rhs._symbol;

  _e = nullptr;

  if (nullptr != _mdl_atom_data) {
    delete _mdl_atom_data;
  }

  _mdl_atom_data = rhs._mdl_atom_data;

  _distance = rhs._distance;

  assert(_distance.ok());

  return;
}

int
Link_Atom::debug_print(std::ostream& os) const
{
  if (_a1 < 0) {
    os << "Link_Atom::debug_print:not initialised\n";
    return os.good();
  }

  os << "Link_Atom::debug_print:between matched atoms " << _a1 << " and " << _a2;
  if (_bond_topology >= 0) {
    os << " bond topology " << _bond_topology;
  }

  if (_symbol.length()) {
    os << ", symbol '" << _symbol << "'";
  }

  os << '\n';

  return _distance.debug_print(os);
}

Link_Atom_Current_State::Link_Atom_Current_State()
{
  _lhs = INVALID_ATOM_NUMBER;
  _rhs = INVALID_ATOM_NUMBER;
  _bt = INVALID_BOND_TYPE;

  return;
}

int
Link_Atom_Current_State::initialise(Link_Atom const* l)
{
  _lhs = l->a1();
  _rhs = l->a2();

  const bond_type_t b = l->btype();

  if (IS_SINGLE_BOND(b)) {
    _bt = SINGLE_BOND;
  } else if (IS_DOUBLE_BOND(b)) {
    _bt = DOUBLE_BOND;
  } else if (IS_TRIPLE_BOND(b)) {
    _bt = TRIPLE_BOND;
  } else {
    _bt = SINGLE_BOND;
  }

  return 1;
}

// #define DEBUG_LINK_ATOMS_SATISFIED

int
Single_Substructure_Query::_link_atoms_satisfied(Query_Atoms_Matched& matched_atoms) const
{
#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Testing " << _link_atom.number_elements() << " link atoms\n";
#endif

  for (const Link_Atom* l : _link_atom) {
    if (!_link_atom_satisfied(*l, matched_atoms)) {
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_link_atom_satisfied(const Link_Atom& l,
                                                Query_Atoms_Matched& matched_atoms) const
{
#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Matched atoms";
  for (int i = 0; i < matched_atoms.number_elements(); i++) {
    const Substructure_Atom* a = matched_atoms[i];

    cerr << ' ' << a->current_hold_atom()->atom_number();
  }
  cerr << '\n';
#endif

  Molecule* m = matched_atoms[0]->current_hold_atom()->m();

  const atom_number_t a1 = matched_atoms[l.a1()]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[l.a2()]->current_hold_atom()->atom_number();

  if (m->fragment_membership(a1) != m->fragment_membership(a2)) {
    //  cerr << "Single_Substructure_Query::_link_atom_satisfied:atoms in different
    //  fragments\n";
    return 0;
  }

  int d = m->bonds_between(a1, a2) - 1;

#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Link atoms " << l.a1() << " and " << l.a2() << '\n';

  cerr << "matched atoms " << a1 << " and " << a2 << " are " << d
       << " bonds apart, satisfies separation " << l.satisfies_separation(d) << '\n';
#endif

  return l.satisfies_separation(d);
}

int
Link_Atom::swap_atoms(atom_number_t n1, atom_number_t n2)
{
  if (_a1 == n1) {
    _a1 = n2;
    return 1;
  }

  if (_a2 == n1) {
    _a2 = n2;
    return 1;
  }

  return 1;
}

int
Link_Atom::atom_is_being_removed(atom_number_t zatom)
{
  assert(zatom != _a1);
  assert(zatom != _a2);

  if (_a1 > zatom) {
    _a1--;
  }

  if (_a2 > zatom) {
    _a2--;
  }

  return 1;
}

int
Link_Atom::initialise_from_mdl_record(const const_IWSubstring& buffer, int matoms,
                                      atom_number_t& a)
{
  if (buffer.nwords() < 6) {
    cerr << "Link_Atom::initialise_from_mdl_record:must have at least 6 words\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);  // M
  buffer.nextword(token, i);  // LIN
  buffer.nextword(token, i);  // sqeuential number, not of interest

  buffer.nextword(token, i);  // atom number

  if (!token.numeric_value(a) || a < 1 || a > matoms) {
    cerr << "MDL_Molecule::_parse_link_record:invalid central atom\n";
    return 0;
  }

  a--;  // convert to our numbering

  buffer.nextword(token, i);  //   range

  int range;
  if (token.numeric_value(range) && range > 0) {
    _distance.set_max(range);
  } else {
    _distance.reset();

    if (_distance.initialise(token)) {
      // convert from atom count to bond count
      for (int j = 0; j < _distance.number_elements(); j++) {
        _distance[j] = _distance[j] + 1;
      }
    } else {
      cerr << "Link_Atom::initialise_from_mdl_record:invalid range '" << token << "'\n";
      return 0;
    }
  }

  assert(_distance.ok());

  buffer.nextword(token, i);
  if (!token.numeric_value(_a1) || _a1 < 1) {
    cerr << "Link_Atom::initialise_from_mdl_record:invalid a1 '" << token << "'\n";
    return 0;
  }

  _a1--;

  buffer.nextword(token, i);
  if (!token.numeric_value(_a2) || _a2 < 1) {
    cerr << "Link_Atom::initialise_from_mdl_record:invalid a2 '" << token << "'\n";
    return 0;
  }

  _a2--;

  // cerr << "Link_Atom:initialised, between " << _a1 << " and " << _a2 << '\n';

  return 1;
}

int
Link_Atom::set_symbol(const const_IWSubstring& s)
{
  _symbol = s;

  _e = get_element_from_symbol_no_case_conversion(s);  // may return NULL, A, Q, [CD4]

  return 1;
}

int
Single_Substructure_Query::add_no_matched_atoms_between_initial_atom_numbers(int n1,
                                                                             int n2)
{
  Substructure_Atom* a1 = query_atom_with_initial_atom_number(n1);
  Substructure_Atom* a2 = query_atom_with_initial_atom_number(n2);

  if (nullptr == a1 || nullptr == a2) {
    cerr << "Single_Substructure_Query::add_no_matched_atoms_between_initial_atom_"
            "numbers:no query atoms found\n";
    cerr << "a1 = " << a1 << " or a2 = " << a2 << '\n';
    return 0;
  }

  Bond* b = new Bond(a1->unique_id(), a2->unique_id(), SINGLE_BOND);

  _no_matched_atoms_between.add(b);

  return _no_matched_atoms_between.number_elements();
}

int
Link_Atom::create_next_variant(MDL_Molecule& m, Link_Atom_Current_State& lacs) const
{
  int maxdist = 0;  // keep the compiler quiet.
  (void)_distance.max(maxdist);

  if (lacs.number_placed() == maxdist) {
    return 0;
  }

  const int natoms = m.natoms();

  // If this is a clean element, just use it, otherwise transfer query specifications

  if (nullptr != _e) {  // atom was C, O, etc...
    m.add(_e);
  } else if (_symbol.length()) {  // atom was *, A, Q, [SD2], etc...
    m.add_atom_based_on_symbol(_symbol);
  } else {
    m.add_atom_based_on_symbol("*");  // how could this happen?
  }

#ifdef DEBUG_CREATE_NEXT_VARIANT
  cerr << "Start with '" << m.smiles() << "'\n";
  cerr << " adding bond " << lacs.lhs() << " and " << natoms << '\n';
#endif

  // When nothing has been placed, there will be a gap between lhs and rhs.
  // Otherwise we need to remove the existing bond to rhs and insert our new atom

  if (0 == lacs.number_placed()) {
    m.add_bond(lacs.lhs(), natoms, lacs.btype_for_molecule(), _bt, _bond_topology);
  } else {
    m.remove_bond_between_atoms(lacs.last_atom_placed(), lacs.rhs());
    m.add_bond(lacs.last_atom_placed(), natoms, lacs.btype_for_molecule(), _bt,
               _bond_topology);
  }

  m.add_bond(natoms, lacs.rhs(), lacs.btype_for_molecule(), _bt, _bond_topology);

#ifdef DEBUG_CREATE_NEXT_VARIANT
  cerr << "create_next_variant created '" << m.smiles() << '\n';
#endif

  lacs.add(natoms);

  if (!m.arrays_allocated()) {
    m.build(m);
  }

  if (nullptr != _mdl_atom_data) {
    m.set_mdl_atom_data(natoms, _mdl_atom_data);
  }

  return 1;
}

/*
  We are deleting the atom which was the link atom. We need to transfer
  any query specifications that may have been associated with it.
*/

int
Link_Atom::set_mdl_atom_data(const MDL_Atom_Data* a)
{
  if (nullptr == _mdl_atom_data) {
    _mdl_atom_data = new MDL_Atom_Data(*a);
  } else {
    (*_mdl_atom_data) = *a;
  }

  return 1;
}

No_Matched_Atoms_Between::No_Matched_Atoms_Between(const int a1, const int a2)
    : _a1(a1), _a2(a2)
{
}

// Starting at s[i], which must be 'open', return the
// index of the matching 'close' character.
// Returns a negative number if not found.
int
MatchingOpenCloseChar(const const_IWSubstring& s, int i, const char open,
                      const char close)
{
  assert(s[i] == open);

  int level = 0;
  for (; i < s.length(); ++i) {
    if (s[i] == close) {
      level--;
      if (level == 0) {
        return i;
      }
    } else if (s[i] == open) {
      level++;
    }
  }

  return -1;  // Did not find matching char.
}

NMAB_Token::NMAB_Token(int a1, int a2) : _a1(a1), _a2(a2)
{
  _op = IW_LOGEXP_UNDEFINED;
  _relational = 0;
}

// Consume characters from 's[i]' into a NMAB_Token.
// 'i' is incremented as characters are recognised.
//  optional leading &,^;
//  optional <>
//  optional digit
//  [smarts]

int
NMAB_Token::Parse(const const_IWSubstring& s, int& i)
{
  assert(s.length() > 0);

  if (s[i] == '&') {
    _op = IW_LOGEXP_AND;
    i++;
  } else if (s[i] == ',') {
    _op = IW_LOGEXP_OR;
    i++;
  } else if (s[i] == '^') {
    _op = IW_LOGEXP_XOR;
    i++;
  } else if (s[i] == ';') {
    _op = IW_LOGEXP_LOW_PRIORITY_AND;
    i++;
  }

  if (i == s.length()) {
    cerr << "NMAB_Token::Parse:incomplete specification '" << s << "'\n";
    return 0;
  }

  // After the leading operator, there may be a relational.

  if (s[i] == '<') {
    _relational = -1;
    i++;
  } else if (s[i] == '>') {
    _relational = 1;
    i++;
  }

  if (i == s.length()) {
    cerr << "NMAB_Token::Parse:incomplete specification '" << s << "'\n";
    return 0;
  }

  // Then a number or a smarts.
  if (s[i] == '[' || isdigit(s[i])) {
  } else {
    cerr << "NMAB_Token::Parse:Unrecognised termination '" << s << " '" << s[i] << "'\n";
    return 0;
  }

  if (isdigit(s[i])) {
    if (!GetNumberOrRange(s, i)) {
      return 0;
    }
  }

  if (_relational != 0 && _number.empty()) {
    cerr << "NMAB_Token::Parse:incomplete relational specification '" << s << "'\n";
    return 0;
  }

  if (i == s.length()) {
    return 1;
  }

  // If no smarts, we are done.
  if (s[i] != '[') {
    return IsComplete();
  }

  const int closing_sq_bracket = misc2::MatchingOpenCloseChar(s, i, '[', ']');
  if (closing_sq_bracket < 0) {
    cerr << "NMAB_Token::Parse:incomplete [smarts] '" << s << "'\n";
    return 0;
  }

  if (i + 1 == closing_sq_bracket) {
    cerr << "NMAB_Token::Parse:empty [smarts] '" << s << "'\n";
    return 0;
  }

  s.from_to(i + 1, closing_sq_bracket - 1, _smarts);

  i = closing_sq_bracket + 1;

  return 1;
}

int
NMAB_Token::GetNumberOrRange(const const_IWSubstring& s, int& i)
{
  assert(isdigit(s[i]));
  auto maybe_a_number = FetchNumeric(s, i);
  if (!maybe_a_number) {
    if (_relational != 0) {
      cerr << "NMAB_Token::Parse:number must follow relational '" << s << "'\n";
      return 0;
    }
  }

  _number.add(*maybe_a_number);
  if (i == s.length()) {
    return 1;
  }

  if (s[i] == '-') {  // A range specification.
    i += 1;
    if (i == s.length()) {
      cerr << "NMAB_Token::Parse:incomplete range specification '" << s << "'\n";
      return 0;
    }
    maybe_a_number = FetchNumeric(s, i);
    if (!maybe_a_number) {
      cerr << "NMAB_Token::Parse:invalid range specification '" << s << "'\n";
      return 0;
    }
    const int number2 = maybe_a_number.value();
    if (number2 < _number.last_item()) {
      cerr << "NMAB_Token::Parse:invalid range values '" << s << "'\n";
      return 0;
    }

    for (int i = _number.last_item() + 1; i <= number2; ++i) {
      _number.add_if_not_already_present(i);
    }
  }

  return 1;
}

bool
NMAB_Token::IsComplete() const
{
  if (_number.size() > 0) {
    return true;
  }
  if (_smarts.length() > 0) {
    return true;
  }

  return false;
}

// A testable function for parsing the {} directive.
// Note argument is a local copy.
int
TokeniseNMABSpecification(const_IWSubstring s, resizable_array_p<NMAB_Token>& tokens)
{
  if (!s.starts_with(kOpenBrace) || !s.ends_with(kCloseBrace)) {
    cerr << "No_Matched_Atoms_Between::Initialise:specification must be of the form "
            "{..}, '"
         << s << "' invalid\n";
    return 0;
  }

  s++;
  s.chop();

  // Or should we silently ignore this??
  if (s.empty()) {
    cerr << "No_Matched_Atoms_Between::Initialise:empty specification\n";
    return 0;
  }

  int i = 0;
  while (i < s.length()) {
    std::unique_ptr<NMAB_Token> token =
        std::make_unique<NMAB_Token>(0, 1);  // arbitrary query atoms numbers.
    if (!token->Parse(s, i)) {
      cerr << "No_Matched_Atoms_Between::Initialise:invalid spec '" << s << "'\n";
      return 0;
    }
    tokens.add(token.release());
  }

  if (tokens.empty()) {
    return 0;
  }

  if (tokens[0]->op() != 0) {
    cerr << "No_Matched_Atoms_Between::Initialise:cannot have leading operator '" << s
         << "'\n";
    return 0;
  }

  for (int i = 1; i < tokens.number_elements(); ++i) {
    if (tokens[i]->op() == 0) {
      cerr << "No_Matched_Atoms_Between::Initialise:non-first tokens must have operator '"
           << s << "'\n";
      return 0;
    }
  }

  return tokens.number_elements();
}

// Parse the {} specification after the three dots.
int
No_Matched_Atoms_Between::Initialise(const ThreeDots& three_dots)
{
  const_IWSubstring s = three_dots.qualifier();
  if (s.empty()) {
    return 1;
  }

  resizable_array_p<NMAB_Token> tokens;
  if (!TokeniseNMABSpecification(s, tokens)) {
    cerr << "No_Matched_Atoms_Between::Initialise:invalid {} directive '" << s << "'\n";
    return 0;
  }

  // If there are consecutive numbers, separated by OR operators, compress those.
  for (int i = 0; i < tokens.number_elements() - 1; ++i) {
    if (tokens[i]->CombineAsRange(*tokens[i + 1])) {
      tokens.remove_item(i + 1);
      i--;
    }
  }

  for (const NMAB_Token* token : tokens) {
    if (token->op()) {
      _logexp.add_operator(token->op());
    }

    std::unique_ptr<NMAB_Operator> op = std::make_unique<NMAB_Operator>();

    if (!op->Build(*token)) {
      cerr << "No_Matched_Atoms_Between::Initialise:invalid token\n";
      return 0;
    }

    _specs.add(op.release());
  }

  // cerr << "No_Matched_Atoms_Between::Initialise: got " << _specs.number_elements() << "items\n";

  return 1;
}

int
NMAB_Token::CombineAsRange(const NMAB_Token& rhs)
{
  if (rhs.op() != IW_LOGEXP_OR) {  // Must have an OR operator.
    return 0;
  }

  if (rhs.relational() != 0) {  // Cannot be combined.
    return 0;
  }

  if (!rhs._smarts.empty()) {  // rhs is not another bond separation.
    return 0;
  }

  if (rhs._number.empty()) {  // No numbers there, cannot combine.
    return 0;
  }

  for (auto i : rhs._number) {
    _number.add_if_not_already_present(i);
  }

  return 1;
}

NMAB_Operator::NMAB_Operator()
{
  _relational = 0;
}

// This is mostly straightforward, except that we recognise some
// non standard 'smarts'.
int
NMAB_Operator::Build(const NMAB_Token& token)
{
  _relational = token.relational();

  _number = token.numbers();

  if (token.smarts().empty()) {
    return 1;
  }

  _query = std::make_unique<Substructure_Atom>();

  // cerr << "Building query from " << token.smarts() << '\n';
  IWString tmp;
  tmp << '[' << token.smarts() << ']';
  if (!_query->construct_from_smarts_token(tmp)) {
    cerr << "NMAB_Token::Build:invalid smarts '" << tmp << "'\n";
    return 0;
  }

  _query->count_attributes_specified();  // Ensure calculated.

  return 1;
}

int
NMAB_Operator::_relational_matches(const int s) const
{
  if (_relational > 0) {
    return s > _number[0];
  }

  if (_relational < 0) {
    return s < _number[0];
  }

  return _number.contains(s);
}

int
NMAB_Operator::Matches(Molecule_to_Match& target, const Set_of_Atoms& unmatched_atoms)
{
#ifdef DEBUG_NMAB_SATISFIED
  cerr << "NMAB_Operator::matches: query? " << (_query != nullptr) << '\n';
  cerr << "Examing " << unmatched_atoms.size() << " unmatched atoms\n";
#endif

  // No substructurre present, numbers are atom separations.
  if (_query == nullptr) { 
    // cerr << "NO query _relational_matches " <<  _relational_matches(unmatched_atoms.number_elements()) << '\n';
    return _relational_matches(unmatched_atoms.number_elements());
  }

  // If no unmatched atoms in the reguion, then all substructure matches will return 0 matches.
  if (unmatched_atoms.empty()) {
    return _relational_matches(0);
  }

  int* already_matched = new_int(target.natoms(), 1);
  std::unique_ptr<int[]> free_already_matched(already_matched);
  unmatched_atoms.set_vector(already_matched, 0);

  // If no number constraint, then all atoms must match.
  const bool all_must_match = _number.empty();

  int matches = 0;
  for (const atom_number_t a : unmatched_atoms) {
    const int tmp = _query->matches(target[a], already_matched);
    //  cerr << " At matched atom " << a << ' ' <<
    //  target.molecule()->smarts_equivalent_for_atom(a) << " match is " << tmp << '\n';
    if (!tmp && all_must_match) {
      return 0;
    }
    matches += tmp;
  }

  // If all must match and we got here, then all did in fact match.
  if (all_must_match) {
    return 1;
  }

#ifdef DEBUG_NMAB_SATISFIED
  cerr << "NMAB_Operator::Matches:Number of matches " << matches << " will return " << _relational_matches(matches) << '\n';
#endif

  return _relational_matches(matches);
}

//#define DEBUG_NO_MATCHED_ATOMS_BETWEEN_MATCHES

int
No_Matched_Atoms_Between::Matches(Molecule_to_Match& target,
                                  Query_Atoms_Matched& matched_atoms, const int* matched,
                                  int* tmp)
{
  Molecule* m = target.molecule();

  const int matoms = m->natoms();
  const atom_number_t a1 = matched_atoms[_a1]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[_a2]->current_hold_atom()->atom_number();

#ifdef DEBUG_NMAB_SATISFIED
  cerr << " matched atoms " << a1 << " " << m->smarts_equivalent_for_atom(a1) << " and "
       << a2 << " " << m->smarts_equivalent_for_atom(a2) << '\n';
#endif

  // probably should check if a1 and a2 are ok atom numbers...

  if (m->fragment_membership(a1) != m->fragment_membership(a2)) {
    return 0;
  }

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);
  m->atoms(atoms);

  std::fill_n(tmp, matoms, 0);
  if (!no_matched_atom_path(atoms, a1, tmp, matched, a2, 0, m->bonds_between(a1, a2))) {
    return 0;
  }

#ifdef DEBUG_NMAB_SATISFIED
  cerr << "Initial nmatched atoms";
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i]) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
#endif

  tmp[a2] = 0;  // temporary kludge

  Set_of_Atoms unmatched_atoms;
  unmatched_atoms.resize(matoms / 2);  // Just a guess.

  for (int i = 0; i < matoms; ++i) {
    if (tmp[i]) {
      unmatched_atoms.add(i);
    }
  }

#ifdef DEBUG_NMAB_SATISFIED
  cerr << " unmatched_atoms " << unmatched_atoms << '\n';
  cerr << "Check " << _specs.number_elements() << " conditions\n";
#endif

  if (_specs.empty()) {  // We have found a path, no constraints on it.
    return 1;
  }

  _logexp.reset();
  for (int i = 0; i < _specs.number_elements(); ++i) {
    cerr << " test " << i << " needed " << _logexp.result_needed(i) << '\n';
    if (!_logexp.result_needed(i)) {
      continue;
    }

    const int m = _specs[i]->Matches(target, unmatched_atoms);
#ifdef DEBUG_NMAB_SATISFIED
    cerr << "    condition " << i << " match " << m << '\n';
#endif
    _logexp.set_result(i, m);

    int rc;
    if (_logexp.evaluate(rc)) {
      return rc;
    }
  }

  return 0;  // Should never come here.
}

int
No_Matched_Atoms_Between::MatchesExhaustive(Molecule_to_Match& target,
                                            Query_Atoms_Matched& matched_atoms,
                                            const int* matched, int* tmp)
{
  Molecule* m = target.molecule();

  const int matoms = m->natoms();
  // cerr << "No_Matched_Atoms_Between::MatchesExhaustive, matched atoms " << _a1 << " and " << _a2 << '\n';

  const atom_number_t a1 = matched_atoms[_a1]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[_a2]->current_hold_atom()->atom_number();

#ifdef DEBUG_NMAB_SATISFIED
  cerr << " matched atoms " << a1 << " " << m->smarts_equivalent_for_atom(a1) << " and "
       << a2 << " " << m->smarts_equivalent_for_atom(a2) << '\n';
#endif

  // probably should check if a1 and a2 are ok atom numbers...

  if (m->fragment_membership(a1) != m->fragment_membership(a2)) {
    return 0;
  }

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);
  m->atoms(atoms);

  std::fill_n(tmp, matoms, 0);
  return _matches_exhaustive(target, matched_atoms, matched, tmp, a1, a2,
                             m->bonds_between(a1, a2));
}

int
No_Matched_Atoms_Between::_matches_exhaustive(Molecule_to_Match& target,
                                              const Query_Atoms_Matched& matched_atoms,
                                              const int* matched, int* in_path,
                                              const atom_number_t my_atom,
                                              const atom_number_t destination,
                                              const int current_distance)
{
  Molecule* m = target.molecule();

#ifdef DEBUG_NMAB_SATISFIED
  cerr << "_matches_exhaustive atom " << my_atom << " heading to " << destination
       << " current_distance " << current_distance << '\n';
#endif

  for (const Bond* b : *m->atomi(my_atom)) {
    const atom_number_t j = b->other(my_atom);
    if (in_path[j]) {
      continue;
    }

#ifdef DEBUG_NMAB_SATISFIED
    cerr << "  to atom " << j << " dist " << m->bonds_between(j, destination)
         << " in path " << in_path[j] << '\n';
#endif

    if (j == destination) {
      if (_is_a_match(target, matched_atoms, in_path)) {
        return 1;
      }
      break;
    } else if (matched[j]) {
      ;
    } else if (m->bonds_between(j, destination) <
               current_distance)  // Must be moving towards destination.
    {
      in_path[j] = 1;
      if (_matches_exhaustive(target, matched_atoms, matched, in_path, j, destination,
                              current_distance - 1)) {
        return 1;
      }
      in_path[j] = 0;
    }
  }

  return 0;
}

int
No_Matched_Atoms_Between::_is_a_match(Molecule_to_Match& target,
                                      const Query_Atoms_Matched& matched_atoms,
                                      const int* in_path)
{
#ifdef DEBUG_NMAB_SATISFIED
  cerr << "_is_a_match\n";
#endif
  if (_specs.empty()) {  // We have found a path, no constraints on it.
    return 1;
  }

  Molecule* m = target.molecule();

  const int matoms = m->natoms();

  Set_of_Atoms unmatched_atoms;
  unmatched_atoms.reserve(matoms / 2);  // Just a guess.

  for (int i = 0; i < matoms; ++i) {
    if (in_path[i]) {
      unmatched_atoms.add(i);
    }
  }

#ifdef DEBUG_NMAB_SATISFIED
  cerr << " unmatched_atoms " << unmatched_atoms << '\n';
  cerr << "Check " << _specs.number_elements() << " conditions\n";
#endif

  _logexp.reset();
  for (int i = 0; i < _specs.number_elements(); ++i) {
#ifdef DEBUG_NMAB_SATISFIED
    cerr << " test " << i << " needed " << _logexp.result_needed(i) << '\n';
#endif
    if (!_logexp.result_needed(i)) {
      continue;
    }

    const int m = _specs[i]->Matches(target, unmatched_atoms);
#ifdef DEBUG_NMAB_SATISFIED
    cerr << "    condition " << i << " match " << m << '\n';
#endif
    _logexp.set_result(i, m);

    int rc;
    if (_logexp.evaluate(rc)) {
      return rc;
    }
  }

  return 0;  // Should never come here.
}

// `region` contain the atoms in the region between atoms `e1` and `e2`.
// Some of those atoms have their `region` value 2, which means on the
// shorest path. Scan those atoms, and add any singly connected, doubly
// bonded atoms to the shortest path set.
int
AddDoublyBonded(Molecule& m, 
                atom_number_t e1, atom_number_t e2,
                int* region) {
  static constexpr int kShortestPath = 2;

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (region[i] != kShortestPath) {
      continue;
    }
    if (i == e1 || i == e2) {
      continue;
    }
    const Atom& a = m[i];
    if (a.ncon() < 3) {
      continue;
    }
    for (const Bond* b : a) {
      if (! b->is_double_bond()) {
        continue;
      }
      const atom_number_t o = b->other(i);
      if (m.ncon(o) == 1) {
        region[o] = kShortestPath;
        ++rc;
      }
    }
  }

  return rc;
}

// We have identified a region between `e1` and `e2` and the atoms in the
// region are set to 1 in `region`.
// For those atoms that are on a shortest path between e1 and e2, set their
// values in `region` to 2 and return the total.
int
IdentifyShortestPath(Molecule& m,
                     atom_number_t e1, atom_number_t e2,
                     int* region) {
  const int d12 = m.bonds_between(e1, e2);
  const int matoms = m.natoms();

  static constexpr int kShortestPath = 2;

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (region[i] == 0) {
      continue;
    }
    if (i == e1 || i == e2) {
      continue;
    }

    if (m.bonds_between(e1, i) + m.bonds_between(i, e2) == d12) {
      region[i] = kShortestPath;
      ++rc;
    }
  }

  return rc;
}

// #define DEBUG_REGION_MATCHES

int
Region::Matches(Molecule_to_Match& target,
                Query_Atoms_Matched & matched_atoms,
                const int* already_matched) const {
  Molecule* m = target.molecule();
  const int matoms = m->natoms();

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(2 * matoms);

  Set_of_Atoms edges;
  for (int i : _matched_atom) {
    edges << matched_atoms[i]->current_hold_atom()->atom_number();
  }

#ifdef DEBUG_REGION_MATCHES
  cerr << "Edges are atoms " << edges << '\n';
#endif

  // Only 2 point supported currently
  if (edges.size() != 2) {
    cerr << "Region::Matches:cannot match " << edges.size() << " regions\n";
    return 0;
  }

  int atoms_in_region = nmab::AtomsInRegion(*m, edges[0], edges[1], tmp.get() + matoms, tmp.get());

#ifdef DEBUG_REGION_MATCHES
  cerr << "Compute atoms_in_region " << atoms_in_region << '\n';
#endif

  // If the region does not exist, we cannot match it.
  if (atoms_in_region == 0) {
    return 0;
  }

  if (_natoms.is_set()) {
    if (! _natoms.matches(atoms_in_region)) {
      return 0;
    }
  }

#ifdef DEBUG_REGION_MATCHES
  for (int i = 0; i < m->natoms(); ++i) {
    cerr << " atm " << i << " region " << tmp[i] << '\n';
  }
#endif

  if (_nrings.is_set()) {
    int matched = 0;
    for (const Ring* r : m->sssr_rings()) {
      // cerr << *r << " set " << r->all_members_non_zero_in_array(tmp.get()) << '\n';
      if (r->all_members_non_zero_in_array(tmp.get())) {
        ++matched;
      }
    }

    if (! _nrings.matches(matched)) {
      return 0;
    }
  }

  if (_atoms_not_on_shortest_path.is_set()) {
    int on_shortest_path = IdentifyShortestPath(*m, edges[0], edges[1], tmp.get());
    // Or should this be an optional behaviour???
    on_shortest_path += AddDoublyBonded(*m, edges[0], edges[1], tmp.get());

    int not_in_shortest_path = atoms_in_region - on_shortest_path;
    if (! _atoms_not_on_shortest_path.matches(not_in_shortest_path)) {
      return 0;
    }
  }

  return 1;
}

NearbyAtoms::NearbyAtoms() {
  _rejection = 0;
  _can_overlap_matched_atoms = false;
}

//#define DEBUG_NEARBY_ATOMS_MATCHES

int
NearbyAtoms::Matches(Molecule_to_Match& target, Query_Atoms_Matched& qam,
                const int* already_matched) {
#ifdef DEBUG_NEARBY_ATOMS_MATCHES
  const int matoms = target.natoms();
  Molecule& m = *target.molecule();
  cerr << "NearbyAtoms::Matches " << m.smiles() << '\n';
  for (int i = 0; i < matoms; ++i) {
    if (already_matched[i]) {
      cerr << " atom " << i << " " << m.smarts_equivalent_for_atom(i) << '\n';
    }
  }
#endif // DEBUG_NEARBY_ATOMS_MATCHES

  Substructure_Results sresults;
  int got_match = 0;

  for (Substructure_Query* q : _query) {
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
#ifdef DEBUG_NEARBY_ATOMS_MATCHES
      cerr << "_can_overlap_matched_atoms " << _can_overlap_matched_atoms << " overlap " << 
          e->any_members_set_in_array(already_matched) << '\n';
#endif
      if (! _can_overlap_matched_atoms && e->any_members_set_in_array(already_matched)) {
        continue;
      }

      if (_matched_atom.empty()) {
        if (! AnyMatchedAtomInRange(target, already_matched, *e)) {
          return 0;
        }
      } else if (! SpecifiedmatchedAtomInRange(target, qam, *e)) {
        return 0;
      }

      ++got_match;
    }
  }

#ifdef DEBUG_NEARBY_ATOMS_MATCHES
  cerr << "Got " << got_match << " mathes to querie, matches " << _hits_needed.matches(got_match) << " rej " << _rejection << '\n';
#endif
  if (_hits_needed.matches(got_match)) {
    return ! _rejection;
  }

  return _rejection;
}

int
NearbyAtoms::AnyMatchedAtomInRange(Molecule_to_Match& target,
                const int* already_matched, const Set_of_Atoms& embedding) const {
  Molecule& m = *target.molecule();
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (! already_matched[i]) {
      continue;
    }

    for (const atom_number_t& a : embedding) {
      const int d = m.bonds_between(i, a);
      if (_bonds_between.matches(d)) {
        return 1;
      }
    }
  }

  return 0;
}

int
NearbyAtoms::SpecifiedmatchedAtomInRange(Molecule_to_Match& target,
                const Query_Atoms_Matched& qam,
                const Set_of_Atoms& embedding) const {
  Molecule& m = *target.molecule();

  for (uint32_t a : _matched_atom) {
    if (! qam.ok_index(a)) {
      cerr << "NearbyAtoms::SpecifiedmatchedAtomInRange:matched atom " << a << " invalid\n";
      cerr << "Query matches only " << qam.size() << " atoms\n";
      return 0;
    }

    const atom_number_t matched = qam[a]->atom_number_matched();

    for (const atom_number_t e : embedding) {
      const int d = m.bonds_between(e, matched);
      if (_bonds_between.matches(d)) {
        return 1;
      }
    }
  }

  return 0;
}
