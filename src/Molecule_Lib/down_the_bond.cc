#include <algorithm>
#include <iostream>
#include <memory>


#define RESIZABLE_ARRAY_IMPLEMENTATION 1

#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "substructure.h"
#include "target.h"

namespace down_the_bond {

using std::cerr;

constexpr char kOpenBrace = '{';
constexpr char kCloseBrace = '}';
constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

void
DownTheBond::DefaultValues() {
  _match_as_match = 1;
  _no_other_substituents_allowed = false;
  _match_individual_substituent = false;
  _all_queries_require_zero_hits = 0;
}

DownTheBond::DownTheBond() {
  _a1 = -1;
  _a2 = -1;
  DefaultValues();
}

DownTheBond::DownTheBond(int a1) {
  _a1 = a1;
  _a2 = -1;
  DefaultValues();
}

// Hold something like a>5;
// This means that the first token holds any operator specification,
// so parsing a>3;h[0] gets parsed into two tokens, with the first one
// holding the low priority and operator ';'.
// Note that currently operators are not supported. There is an implicit
// and operation for all specifications.
class DTB_Token {
  public:
    enum class Token {
      kUndefined = 0,
      kNatoms = 1,
      kheteroatoms = 2,
      kRingAtoms = 3,
      kUnsaturation = 4,
      kAromatic = 5,
      kMaxDistance = 6,
      kSmarts = 7
    };
    enum class Operator {
      kUndefined = 0,
      kHighPriorityAnd = 1,
      kOr = 2,
      kXor = 3,
      kLowPriorityAnd = 4
    };

    IWString _smarts;

  private:
    // The compiler will probably pad this for alignment.
    const Token _type;
    iwmatcher::Matcher<uint32_t> _numeric;
    Operator _op;

  // Private functions
    int ParseRange(const const_IWSubstring& buffer, int& ndx);
    int ParseMinMax(const const_IWSubstring& buffer, int& ndx);

  public:
    DTB_Token(const Token token);

    int Build(const const_IWSubstring& buffer, int& ndx);

    Token type() const {
      return _type;
    }

    void set_smarts(const const_IWSubstring& buffer, int istart, int nchars) {
      _smarts.set(buffer.data() + istart, nchars);
    }
    const IWString& smarts() const {
      return _smarts;
    }

    iwmatcher::Matcher<uint32_t> numeric() const {
      return _numeric;
    }

    void set_operator(Operator s) {
      _op = s;
    }

    // Initialise a matcher based on the values held.
    template <typename T> void set_matcher(iwmatcher::Matcher<T>& matcher) const;
};

DTB_Token::DTB_Token(const Token token) : _type(token) {
  _op = Operator::kUndefined;
}

template <typename T>
void
DTB_Token::set_matcher(iwmatcher::Matcher<T>& matcher) const {
  matcher = _numeric;
}

int
GetNumber(const const_IWSubstring& buffer, int& ndx, uint32_t& value) {
  value = 0;
  int rc = 0;
  while (ndx < buffer.length() && isdigit(buffer[ndx])) {
    value = value * 10 + buffer[ndx] - '0';
    ++ndx;
    ++rc;
  }

  return rc;
}

// The letter part of the token has already been set, parse the numeric
// qualifier
//   <3
//   >4
//   6
//   {3-5}

int
DTB_Token::Build(const const_IWSubstring& buffer, int& ndx) {
  if (ndx >= buffer.length()) {
    return 0;
  }

  const char c = buffer[ndx];
  if (isdigit(c)) {
    uint32_t value = 0;
    if (! GetNumber(buffer, ndx, value)) {
      return 0;
    }
    _numeric.add(value);

    return 1;
  }

  //  {3-5}, <4 >3
  if (c == kOpenBrace) {
    return ParseRange(buffer, ndx);
  } else if (c == '<') {
    return ParseMinMax(buffer, ndx);
  } else if (c == '>') {
    return ParseMinMax(buffer, ndx);
  } else {
    cerr << "DTB_Token::Build:unrecognised numeric qualifier '" << c << "'\n";
    return 0;
  }
}

// ndx will be pointing at '>' or '<'
int
DTB_Token::ParseMinMax(const const_IWSubstring& buffer, int& ndx) {
  if (ndx == buffer.length() - 1) {
    cerr << "DTB_Token::ParseMinMax:no numeric qualifier\n";
    return 0;
  }

  const char minmax = buffer[ndx];
  ++ndx;
  uint32_t value;
  if (! GetNumber(buffer, ndx, value)) {
    cerr << "DTB_Token::ParseMinMax:invalid numeric\n";
    return 0;
  }

  if (minmax == '<') {
    _numeric.set_max(value - 1);
  } else  if (minmax == '>') {
    _numeric.set_min(value + 1);
  } else {
    cerr << "DTB_Token::ParseMinMax:invalid operator '" << minmax << "'\n";
    return 0;
  }

  return 1;
}

// THis is not a robust parser. TODO:ianwatson see if we can use the
// code implemented in the substructure search code.
int
DTB_Token::ParseRange(const const_IWSubstring& buffer, int& ndx) {
  assert(buffer[ndx] == kOpenBrace);
  ++ndx;
  if (ndx >= buffer.length() - 1) {
    cerr << "DTB_Token::ParseRange:invalid range\n";
    return 0;
  }

  uint32_t minval = 0;
  if (! GetNumber(buffer, ndx, minval)) {
    cerr << "DTB_Token::ParseRange:invalid minval\n";
    return 0;
  }

  if (ndx >= buffer.length()) {
    return 0;
  }
  
  // {n} is OK, interpreted as a single value.
  if (buffer[ndx] == kCloseBrace) {
    ++ndx;
    _numeric.add(minval);
    return 1;
  }

  if (buffer[ndx] != '-') {
    cerr << "DTB_Token::ParseRange:not - separating values, got '" << buffer[ndx] << "'\n";
    return 0;
  }
  ++ndx;
  if (ndx >= buffer.length()) {
    return 0;
  }
  if (buffer[ndx] == kCloseBrace) {
    _numeric.set_min(minval);
    ++ndx;
    return 1;
  }

  uint32_t maxval = 0;
  if (! GetNumber(buffer, ndx, maxval)) {
    cerr << "DTB_Token::ParseRange:invalid maxval\n";
    return 0;
  }

  if (minval > maxval) {
    cerr << "DTB_Token::ParseRange:invalid range " << minval << ' ' << maxval << '\n';
    return 0;
  }

  if (buffer[ndx] != kCloseBrace) {
    cerr << "DTB_Token::ParseRange:not closing brace\n";
    return 0;
  }

  ++ndx;

  _numeric.set_min(minval);
  _numeric.set_max(maxval);

  return 1;
}
int
QueryMatches::Build(const IWString& smarts, const iwmatcher::Matcher<uint32_t>& numeric) {
  _query = std::make_unique<Substructure_Atom>();
  if (! _query->construct_from_smarts_token(smarts.data(), smarts.size())) {
    cerr << "QueryMatches::Build:invalid smarts '" << smarts << "'\n";
    return 0;
  }
  _query->count_attributes_specified();

  _hits_needed = numeric;

  return 1;
}

int
QueryMatches::RequiresZeroHits() const {
  if (! _hits_needed.is_set()) {
    return 0;
  }
  if (_hits_needed.match_any()) {
    return 0;
  }

  if (! _hits_needed.matches(0)) {
    return 0;
  }

  if (_hits_needed.size() != 1) {
    return 0;
  }

  uint32_t tmp;
  if (_hits_needed.min(tmp) || _hits_needed.max(tmp)) {
    return 0;
  }

  return 1;
}

int
FindClosingBracket(const const_IWSubstring& buffer, int ndx) {
  assert(buffer[ndx] == kOpenSquareBracket);

  const int nchars = buffer.length();
  int bracket_level = 1;

  for (int i = ndx + 1; i < nchars; ++i) {
    if (buffer[i] == kOpenSquareBracket) {
      ++bracket_level;
    }
    if (buffer[i] == kCloseSquareBracket) {
      --bracket_level;
      if (bracket_level == 0) {
        return i - ndx + 1;
      }
    }
  }

  return -1;
}

int
ParseToDTB(const const_IWSubstring& buffer,
            resizable_array_p<DTB_Token>& tokens) {
  int nchars = buffer.length();
  for (int i = 0; i < nchars; ++i) {
    std::unique_ptr<DTB_Token> token;
    int chars_consumed = 1;
    switch (buffer[i]) {
      case 'a':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kNatoms);
        break;
      case 'h':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kheteroatoms);
        break;
      case 'r':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kRingAtoms);
        break;
      case 'u':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kUnsaturation);
        break;
      case 'm':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kAromatic);
        break;
      case 'd':
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kMaxDistance);
        break;
      case kOpenSquareBracket:
        token = std::make_unique<DTB_Token>(DTB_Token::Token::kSmarts);
        chars_consumed = FindClosingBracket(buffer, i);
        if (chars_consumed < 0) {
          cerr << "ParseToDTB:no closing brace '" << buffer << "'\n";
          return 0;
        }
        token->set_smarts(buffer, i, i + chars_consumed);
        break;
      default:
        cerr << "ParseDTB:unrecognised property '" << buffer[i] << "'\n";
        return 0;
    }
    i += chars_consumed;
    int isave = i;
    if (! token->Build(buffer, i)) {
      cerr << "ParseDTB:invalid specification '" << buffer << "'\n";
      cerr << "                                 ";
      for (int i = 0; i < isave; ++i) {
        cerr << ' ';
      }
      cerr << "^\n";
      return 0;
    }

    tokens << token.release();

    if (i == buffer.length() - 1) {
      break;
    }
    if (buffer[i] == ';') {
      tokens.back()->set_operator(DTB_Token::Operator::kLowPriorityAnd);
    } else if (buffer[i] == '&') {
      tokens.back()->set_operator(DTB_Token::Operator::kHighPriorityAnd);
    } else if (buffer[i] == '|') {
      tokens.back()->set_operator(DTB_Token::Operator::kOr);
    } else if (buffer[i] == '^') {
      tokens.back()->set_operator(DTB_Token::Operator::kXor);
    }
  }

  return tokens.size();
}

// #define DEBUG_DOWN_THE_BOND_BUILD

// A DownTheBond directive consists of a token followed
// by a numeric qualifier.
// a2, a>3, a<4, a[4-9]

// The following tokens are recognised.
//   a  total number of atoms seen down the bond.
//   h  total number of heteroatoms seen down the bond.
//   r  total number of ring atoms seen down the bond.
//   m  total number of aromatic atoms seen down the bond.
//   u  total number of usaturated atoms seen down the bond.
//   d  longest distance of any atom from the _a2.
//   [smarts] smarts for atoms seen down the bond.
// For now, logical operators are not supported. All tokens
// must be separated by a ';' operator.
// a[3-9];r0;[n]1;d<5
int
DownTheBond::Build(const const_IWSubstring& buffer) {
  const_IWSubstring mybuffer(buffer);
  if (mybuffer[0] == '!') {
    _match_as_match = 0;
  }

  resizable_array_p<DTB_Token> tokens;
  if (! ParseToDTB(mybuffer, tokens)) {
    cerr << "DownTheBond::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }
#ifdef DEBUG_DOWN_THE_BOND_BUILD
  cerr << "Got " << tokens.size() << " Down the bond tokens\n";
#endif

  // We cannot multiply specify a condition.
  resizable_array<int> seen;
  seen.reserve(tokens.size());

  for (const DTB_Token* dtb : tokens) {
    if (dtb->type() == DTB_Token::Token::kUndefined) {
      return 0;
    }
    // OK to have multiple smarts
    if (dtb->type() == DTB_Token::Token::kSmarts) {
      continue;
    }
    if (! seen.add_if_not_already_present(static_cast<int>(dtb->type()))) {
      cerr << "DownTheBond::Build:duplicate specification\n";
      return 0;
    }
  }

  for (const DTB_Token* dtb : tokens) {
    switch (dtb->type()) {
      case DTB_Token::Token::kUndefined:
        return 0;
      case DTB_Token::Token::kNatoms:
        dtb->set_matcher(_natoms);
        break;
      case DTB_Token::Token::kRingAtoms:
        dtb->set_matcher(_ring_atom_count);
        break;
      case DTB_Token::Token::kheteroatoms:
        dtb->set_matcher(_heteroatom_count);
        break;
      case DTB_Token::Token::kUnsaturation:
        dtb->set_matcher(_unsaturation_count);
        break;
      case DTB_Token::Token::kAromatic:
        dtb->set_matcher(_aromatic_count);
        break;
      case DTB_Token::Token::kMaxDistance:
        dtb->set_matcher(_max_distance);
        break;
      case DTB_Token::Token::kSmarts:
        std::unique_ptr<QueryMatches> qm = std::make_unique<QueryMatches>();
        if (! qm->Build(dtb->smarts(), dtb->numeric())) {
          cerr << "DownTheBond::Build:bad smarts '" << dtb->smarts() << "'\n";
          return 0;
        }
        _query << qm.release();
        break;
    }
  }

  _all_queries_require_zero_hits = AllQueriesRequireZeroHits();

#ifdef DEBUG_DOWN_THE_BOND_BUILD
  cerr << "DownTheBond::Build:added " << _query.size() << " query specifications, _all_queries_require_zero_hits " << _all_queries_require_zero_hits << '\n';
#endif

  return 1;
}

static std::optional<int>
IdentifyDTB(const Molecule& m,
            atom_number_t zatom,
            atom_number_t previous_atom,
            atom_number_t avoid,
            int* visited) {
  int rc = 1;
  visited[zatom] = 1;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    // cerr << "From " << zatom << " to " << o << " previous_atom " << previous_atom << " visited " << visited[o] << " avoid " << avoid << '\n';
    if (o == previous_atom) {
      continue;
    }
    if (visited[o]) {
      continue;
    }
    if (o == avoid) {
      return std::nullopt;
    }

    std::optional<int> tmp = IdentifyDTB(m, o, zatom, avoid, visited);
    if (! tmp) {
      return std::nullopt;
    }
    rc += *tmp;
  }

  return rc;
}

// #define DEBUG_MATCHES_INDIVIDUAL_SUBSTITUENT

// THis is somewhat inefficient.
// If _no_other_substituents_allowed is set, then once we find a non-match
// we should return 0. But instead we continue looking for matches.
// This is just about code complexity.
int
DownTheBond::MatchesIndividualSubstituent(Molecule& m,
                atom_number_t a1,
                atom_number_t a2,
                int* visited) {

#ifdef DEBUG_MATCHES_INDIVIDUAL_SUBSTITUENT
  cerr << "DownTheBond::MatchesIndividualSubstituent:from " << a1 << " to " << a2 << '\n';
#endif
  int attachments_found = 0;
  int attachments_matched = 0;
  for (const Bond* b : m[a2]) {
    const atom_number_t o = b->other(a2);
    if (o == a1) {
      continue;
    }

    std::fill_n(visited, m.natoms(), 0);
    visited[a1] = 1;

    std::optional<int> natoms = IdentifyDTB(m, o, a2, a1, visited);
    if (! natoms) {
      // cerr << "No fragment found to atom " << o << ' ' << m.smarts_equivalent_for_atom(o) << '\n';
      continue;
    }
#ifdef DEBUG_MATCHES_INDIVIDUAL_SUBSTITUENT
    for (int i = 0; i < m.natoms(); ++i) {
      cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " visited " << visited[i] << '\n';
    }
#endif
    visited[a1] = 0;
    visited[a2] = 1;

#ifdef DEBUG_MATCHES_INDIVIDUAL_SUBSTITUENT
    cerr << "Found " << *natoms << " atoms in attachment\n";
#endif

    ++attachments_found;
    int matched_here = 0;
    // Need to add 1 to natoms because a2 was not counted above.
    // All other properties are computed from the visited array.
    if (_natoms.is_set() && ! _natoms.matches(*natoms + 1)) {
      continue;
    } else {
      matched_here = 1;
    }
    if (_heteroatom_count.is_set() && ! OkHeteratomCount(m, visited)) {
      continue;
    } else {
      matched_here = 1;
    }
    if (_ring_atom_count.is_set() && ! OkRingAtomCount(m, visited)) {
      continue;
    } else {
      matched_here = 1;
    }
    if (_unsaturation_count.is_set() && ! OkUnsaturationCount(m, visited)) {
      continue;
    } else {
      matched_here = 1;
    }
    if (_aromatic_count.is_set() && ! OkAromaticCount(m, visited)) {
      continue;
    } else {
      matched_here = 1;
    }
    if (_max_distance.is_set() && ! OkMaxDistance(m, a2, visited)) {
      continue;
    } else {
      matched_here = 1;
    }

    if (matched_here) {
      ++attachments_matched;
    }
  }

#ifdef DEBUG_MATCHES_INDIVIDUAL_SUBSTITUENT
  cerr << "Found " << attachments_found << " attachments, " << attachments_matched << " matched\n";
#endif
  if (attachments_matched == 0) {
    return ! _match_as_match;
  }

  if (_no_other_substituents_allowed && attachments_matched < attachments_found) {
    return ! _match_as_match;
  }

  return 1;
}

int
DownTheBond::OkHeteratomCount(const Molecule& m, const int* visited) const {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (visited[i] == 0) {
      continue;
    }

    if (m.atomic_number(i) != 6) {
      ++rc;
    }
  }

  return _heteroatom_count.matches(rc);
}

int
DownTheBond::OkUnsaturationCount(Molecule& m, const int* visited) const {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (visited[i] == 0) {
      continue;
    }

    if (m.is_aromatic(i)) {
      continue;
    }

    if (m.unsaturation(i)) {
      ++rc;
    }
  }

  return _unsaturation_count.matches(rc);
}

int
DownTheBond::OkAromaticCount(Molecule& m, const int* visited) const {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (visited[i] == 0) {
      continue;
    }

    if (m.is_aromatic(i)) {
      ++rc;
    }
  }

  return _aromatic_count.matches(rc);
}

int
DownTheBond::OkMaxDistance(Molecule& m, atom_number_t a2, const int* visited) const {
  const int matoms = m.natoms();
  int maxdist = 0;
  for (int i = 0; i < matoms; ++i) {
    if (! visited[i]) {
      continue;
    }
    const int d = m.bonds_between(a2, i);
    if (d > maxdist) {
      maxdist = d;
    }
  }

  return _max_distance.matches(maxdist);
}

int
DownTheBond::OkRingAtomCount(Molecule& m, const int* visited) const {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (visited[i] == 0) {
      continue;
    }

    if (m.ring_bond_count(i) > 0) {
      ++rc;
    }
  }

  return _ring_atom_count.matches(rc);
}

// #define DEBUG_DOWN_THE_BOND_MATCHES

int
DownTheBond::Matches(Molecule_to_Match& target,
                      Query_Atoms_Matched & matched_atoms,
                      int * visited) {
  Molecule& m = *target.molecule();

  if (! matched_atoms.ok_index(_a1) || ! matched_atoms.ok_index(_a2)) {
    cerr << "DownTheBond:Matches:invalid matched query atom number " << _a1 << " or " << _a2 << '\n';
    cerr << "Have " << matched_atoms.size() << " matched query atoms\n";
    return 0;
  }

  const atom_number_t a1 = matched_atoms[_a1]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[_a2]->current_hold_atom()->atom_number();
#ifdef DEBUG_DOWN_THE_BOND_MATCHES
  cerr << "atoms " << a1 << " and " << a2 << '\n';
#endif

  if (! m.ok_atom_number(a1) || ! m.ok_atom_number(a2)) {
    cerr << "DownTheBond::Matches:invalid atom number " << a1 << " or " << a2 << '\n';
    return 0;
  }

  // cerr << "_match_individual_substituent " << _match_individual_substituent << '\n';
  if (_match_individual_substituent) {
    return MatchesIndividualSubstituent(m, a1, a2, visited);
  }

  const bool compute_natoms = _natoms.is_set();
  const bool compute_heteroatoms = _heteroatom_count.is_set();
  const bool compute_ring_atoms = _ring_atom_count.is_set();
  const bool compute_unsaturation = _unsaturation_count.is_set();
  const bool compute_aromatic = _aromatic_count.is_set();

  resizable_array<atom_number_t> atom_stack;
  int heteroatoms = m.atomic_number(a2) != 6;
  int ring_atoms = m.is_ring_atom(a2);
  int aromatic = 0;
  int unsaturation = 0;
  if (m.is_aromatic(a2)) {
    aromatic = 1;
  } else if (! m.saturated(a2)) {
    unsaturation = 1;
  }

  const Atom& atom2 = m.atom(a2);
  for (const Bond* b : atom2) {
    const atom_number_t j = b->other(a2);
    if (j == a1) {
      continue;
    }
    atom_stack << j;
  }

  if (atom_stack.empty()) {
    return NoAtomsDownTheBond(m, a1, a2);
  }

  std::fill_n(visited, m.natoms(), 0);

  visited[a1] = 1;
  visited[a2] = 1;
  int number_visited = 1;
  while (! atom_stack.empty()) {
    const atom_number_t i = atom_stack.pop();
    if (visited[i]) {
      continue;
    }
    visited[i] = 1;
    ++number_visited;
    const Atom& atom = m.atom(i);
    if (compute_heteroatoms && atom.atomic_number() != 6) {
      ++heteroatoms;
    }
    if (compute_ring_atoms && m.ring_bond_count(i)) {
      ++ring_atoms;
    }
    if (compute_unsaturation && ! m.is_aromatic(i) && atom.unsaturated()) {
      ++unsaturation;
    }
    if (compute_aromatic && m.is_aromatic(i)) {
      ++aromatic;
    }

    for (const Bond* bond : atom) {
      atom_number_t j = bond->other(i);
      if (j == a1) {  // Must be part of a loop, must fail.
        return ! _match_as_match;
      }
      if (visited[j]) {
        continue;
      }
      atom_stack << j;
    }
  }

#ifdef DEBUG_DOWN_THE_BOND_MATCHES
  cerr << "From atom " << a1 << " visited " << number_visited << " atoms\n";
  cerr << "heteroatoms " << compute_heteroatoms << " value " << heteroatoms << '\n';
  cerr << "_heteroatom_count.matches " << _heteroatom_count.matches(heteroatoms) << '\n';
#endif

  if (compute_heteroatoms && ! _heteroatom_count.matches(heteroatoms)) {
    return ! _match_as_match;
  }
  if (compute_ring_atoms && ! _ring_atom_count.matches(ring_atoms)) {
    return ! _match_as_match;
  }
  if (compute_unsaturation && ! _unsaturation_count.matches(unsaturation)) {
    return ! _match_as_match;
  }
  if (compute_aromatic && ! _aromatic_count.matches(aromatic)) {
    return ! _match_as_match;
  }
  if (compute_natoms && ! _natoms.matches(number_visited)) {
    return ! _match_as_match;
  }

#ifdef DEBUG_DOWN_THE_BOND_MATCHES
  cerr << "Have " << _query.size() << " query constraints\n";
#endif
  if (_query.empty()) {
    return _match_as_match;
  }

  const int matoms = m.natoms();
  // If any atoms have been visited, we run the queries.
  // If no atoms have been visited, and all queries are for 0 occurrences, that is ok.
  // If no atoms have been visited, and all queries require a positive match, that is a fail.

  // First case, some atoms visited.
  if (std::any_of(visited, visited + matoms, [](int v) {
      return v == 1;
    })) {
  } else if (_all_queries_require_zero_hits) {
    return _match_as_match;
  } else {
    return !_match_as_match;;
  }

  std::unique_ptr<int[]> already_matched = std::make_unique<int[]>(matoms);

  for (QueryMatches* qm : _query) {
    uint32_t nhits = 0;
    for (int i = 0; i < matoms; ++i) {
      // cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " visited " << visited[i] << '\n';
      if (!visited[i]) {
        continue;
      }
      std::fill_n(already_matched.get(), matoms, 0);
      if (qm->ss_atom().matches(target[i], already_matched.get())) {
        ++nhits;
      }
    }
#ifdef DEBUG_DOWN_THE_BOND_MATCHES
    cerr << "Got " << nhits << " hits\n";
#endif
    if (! qm->Matches(nhits)) {
      return ! _match_as_match;
    }
  }

  return _match_as_match;
}

// If _query is empty, it does not matter what this function does.
int
DownTheBond::AllQueriesRequireZeroHits() const {
  for (const QueryMatches* q : _query) {
    if (!q->RequiresZeroHits()) {
      return 0;
    }
  }

  return 1;
}

// Atom a2 is terminal. There are no further atoms. By definition, atom a2 is
// included in what is downt he bond, so we have 1 matched atom.
int
DownTheBond::NoAtomsDownTheBond(Molecule& m, atom_number_t a1, atom_number_t a2) {
  int positive_match_found = 0;

  if (! _natoms.is_set()) {
  } else if (_natoms.matches(1)) {
    ++positive_match_found;
  } else {
    return ! _match_as_match;
  }

  if (! _heteroatom_count.is_set()) {
  } else {
    int h = (m.atomic_number(a2) != 6);
    if (_heteroatom_count.matches(h)) {
      ++positive_match_found;
    } else {
      return ! _match_as_match;
    }
  }

  if (! _ring_atom_count.is_set()) {
  } else if (_ring_atom_count.matches(0)) {
    ++positive_match_found;
  } else {
    return ! _match_as_match;
  }

  if (! _unsaturation_count.is_set()) {
  } else {
    int u = 0;
    if (m.is_aromatic(a2)) {
      u = 0;
    } else if (m.saturated(a2)) {
      u = 0;
    } else {
      u = 1;
    }
    if (_unsaturation_count.matches(u)) {
      ++positive_match_found;
    } else {
      return ! _match_as_match;
    }
  }

  if (! _aromatic_count.is_set()) {
  } else if (_aromatic_count.matches(0)) {
    ++positive_match_found;
  } else {
    return ! _match_as_match;
  }

  if (! _max_distance.is_set()) {
  } else if (_max_distance.matches(0)) {
    ++positive_match_found;
  } else {
    return ! _match_as_match;
  }

  // Nothing has rejected the value. If there was nothing specified, that is a match.
  if (positive_match_found == 0) {
    return 1;
  }

  return _match_as_match;
}

}  // namespace down_the_bond

using down_the_bond::DownTheBond;

int
Single_Substructure_Query::_down_the_bond_satisfied(Molecule_to_Match& target,
                Query_Atoms_Matched& matched_atoms) const {
  std::unique_ptr<int[]> visited = std::make_unique<int[]>(target.natoms());

  for (DownTheBond * dtb : _down_the_bond) {
    if (! dtb->Matches(target, matched_atoms, visited.get())) {
      return 0;
    }
  }

  return 1;
}
