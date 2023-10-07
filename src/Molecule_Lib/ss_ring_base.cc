#include <ctype.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "path.h"
#include "substructure.h"
#include "target.h"

using std::cerr;

Substructure_Ring_Base::Substructure_Ring_Base()
{
  _match_as_match_or_rejection = 1;

  _all_hits_in_same_fragment = 0;

  _only_keep_matches_in_largest_fragment = 0;

  _is_heteroatom = new_int(HIGHEST_ATOMIC_NUMBER + 1, 1);
  _is_heteroatom[6] = 0;

//  Should we include Hydrogen?

  _environment_can_match_in_ring_atoms = 0;

  _set_global_id = 0;

  _ring_extends_to_carbonyl = false;

  return;
}

Substructure_Ring_Base::~Substructure_Ring_Base()
{
  if (nullptr != _is_heteroatom)
    delete [] _is_heteroatom;

  return;
}

//#define DEBUG_ENVIRONMENT_MATCHES

int
Substructure_Ring_Base::_environment_matches(Molecule_to_Match & target,
                                             Query_Atoms_Matched & matched_query_atoms,
                                             int * previously_matched_atoms)
{
  int rc = 0;

  int atom_to_process = 0;

#ifdef DEBUG_ENVIRONMENT_MATCHES
  cerr << "Processing " << matched_query_atoms.number_elements() << " environment children\n";
#endif

  while (atom_to_process >= 0)
  {
#ifdef DEBUG_ENVIRONMENT_MATCHES
    cerr << "Processing child " << atom_to_process << '\n';
#endif

    Substructure_Atom * a = const_cast<Substructure_Atom *>(matched_query_atoms[atom_to_process]);
#ifdef DEBUG_ENVIRONMENT_MATCHES
    cerr << "Matched " << a->is_matched() << " anchor " << a->anchor()->atom_number() << " Children= " << a->number_children() << '\n';
#endif

    if (! a->move_to_next_match_from_current_anchor(previously_matched_atoms, matched_query_atoms)) {
#ifdef DEBUG_ENVIRONMENT_MATCHES
      cerr << "Move to next failed for query environment atom " << a->unique_id() << '\n';
#endif

      a->remove_your_children(matched_query_atoms, previously_matched_atoms);

      atom_to_process--;
      if (atom_to_process < 0) {
        return rc;
      }
    }
    else
    {
#ifdef DEBUG_ENVIRONMENT_MATCHES
      cerr << "Move to next match succeeded " << a->unique_id() << "(" << a->atom_number_matched() <<
              "), or = " << a->or_id() << " atom to process = " << atom_to_process << " matched = " << matched_query_atoms.size() << '\n';
#endif

      a->add_your_children(matched_query_atoms);   // does nothing if already added

      atom_to_process++;

      if (atom_to_process >= matched_query_atoms.number_elements())      // the == condition is the only one which should ever happen
      {
        rc++;
        break;        // break from while (atom_to_process >= 0) loop, we are only interested in one embedding per start atom
      }
    }
  }

#ifdef DEBUG_ENVIRONMENT_MATCHES
  cerr << "Substructure_Ring_Base::_environment_matches: returning " << rc << '\n';
#endif

  return rc;
}

/*
  Dec 2008. Initially I used ring as a non-const array. Once I found a
  match at a given point in ring[] I set that entry to zero. But that has
  the undesirable effect of preventing disubstituted things from counting
  properly. Previously these recorded one spinach group, now they report
  two. We could theoretically change this to an optional behaviour
  OP1(=O)Oc2c(O1)cccc2 PBCHM78561
  OP1(=O)COc2c(OC1)cccc2 PBCHM746738
*/

int
Substructure_Ring_Base::_environment_matches(Molecule_to_Match & target,
                                             Substructure_Atom & root_atom,
                                             const int * ring,
                                             int * already_matched)
{
#ifdef DEBUG_ENVIRONMENT_MATCHES
  cerr << "Start environment match, root atom has " << root_atom.attributes_specified() << " attributes specified\n";
#endif

  Query_Atoms_Matched qam;    // scope here just for efficiency
  qam.resize(20);            // 20 seems pretty large

  const int matoms = target.natoms();

  int * copy_ring = new int[matoms]; std::unique_ptr<int[]> free_copy_ring(copy_ring);
  copy_vector(copy_ring, ring, matoms);

  int nhits = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == copy_ring[i]) {     // not in the ring or already matched
      continue;
    }

    Target_Atom & a = target[i];

    copy_vector(already_matched, copy_ring, matoms);

    already_matched[i] = 0;    // the matches() function needs it this way

    if (! root_atom.matches(a, already_matched)) {
      continue;
    }

    root_atom.set_hold(&a, already_matched);

    if (qam.number_elements()) {
      qam.resize_keep_storage(0);
    }

    if (0 == root_atom.add_your_children(qam))    // root atom only, no children
    {
#ifdef DEBUG_ENVIRONMENT_SEARCH
      cerr << "Root atom hit, nhits = " << (nhits + 1) << '\n';
#endif

      nhits++;
      root_atom.recursive_release_hold();
      continue;
    }

    if (_environment_can_match_in_ring_atoms)
    {
      set_vector(already_matched, matoms, 0);
      already_matched[i] = 1;
    }

    int tmp = _environment_matches(target, qam, already_matched);

#ifdef DEBUG_ENVIRONMENT_SEARCH
    cerr << "At atom " << i << " matches = " << tmp << '\n';
#endif

    root_atom.recursive_release_hold();

    already_matched[i] = 1;

    if (tmp)
      copy_ring[i] = 0;

    nhits += tmp;
  }

#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "Substructure_Ring_Base::_environment_matches: returning " << nhits << '\n';
#endif

  return nhits;
}

/*
  Note some important design issues here. We do a greedy match. The first environment
  component finds a match. Then the second environment component tries to find a
  different position to match. This could create problems in cases like

  *-N||*-NO2

  It would be random as to whether or not such a query would match nitro-aniline - would
  depend on which one was encountered first. Make sure you are specific when
  multiple environments are present
*/

//#define DEBUG_L1_ENVIRONMENT_MATCHES

int
Substructure_Ring_Base::_environment_matches(Molecule_to_Match & target,
                                             int * ring,
                                             int * already_matched)
{
  int ne = _environment_atom.number_elements();
#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
  cerr << "Substructure_Ring_Base::_environment_matches: have " << ne << " components to match\n";
#endif

  for (int i = 0; i < ne; i++) {
    // should check to see if this result is actually needed by _environment_logexp. risky, maybe TODO...
    Substructure_Atom * r = _environment_atom[i];

    int nhits = _environment_matches(target, *r, ring, already_matched);

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    cerr << "Result for component " << i << " is " << nhits << '\n';
#endif

    if (_environment_numerical_requirement[i]->is_set()) {
      nhits = _environment_numerical_requirement[i]->matches(nhits);
    }

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    if (_environment_numerical_requirement[i]->is_set())
      cerr << "After filtering by numerical requirement result is " << nhits << '\n';
#endif

    _environment_logexp.set_result(i, nhits);

    int zresult;

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    if (_environment_logexp.evaluate(zresult)) {
      cerr << "Result available, will return " << zresult << '\n';
    }
#endif

    if (_environment_logexp.evaluate(zresult)) {
      return zresult;
    }
  }

  return 0;
}

/*
  RING is an array over the molecule. Atoms which are in either the single ring,
  or the ring system, will be set.
*/

int
Substructure_Ring_Base::_environment_matches(Molecule_to_Match & target,
                                             int * ring)
{
#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
  cerr << "Starting _environment_matches, have " << _environment_atom.size() << " components\n";
#endif

  _environment_logexp.reset();

  if (_environment_atom.empty()) {
    return 1;
  }

  int * already_matched = new int[target.natoms()]; std::unique_ptr<int[]> free_already_matched(already_matched);

  return _environment_matches(target, ring, already_matched);
}

/*
  renv will look like
  1c-F&&<3C-[OH]...
*/

static int
is_operator_character(char c)
{
  if ('&' == c)
    return 1;

  if ('|' == c)
    return 1;

  if ('^' == c)
    return 1;

  if (';' == c)
    return 1;
   
  return 0;
}

static int
is_logexp_operator(const const_IWSubstring & s,
                   int offset,
                   int & zop)
{
  if (s.matches_at_position(offset, "&&"))
    zop = IW_LOGEXP_AND;
  else if (s.matches_at_position(offset, "||"))
    zop = IW_LOGEXP_OR;
  else if (s.matches_at_position(offset, "^^"))
    zop = IW_LOGEXP_XOR;
  else if (s.matches_at_position(offset, ";;"))
    zop = IW_LOGEXP_LOW_PRIORITY_AND;
  else
    return 0;

  return 1;
}

static int
count_components(const const_IWSubstring & renv)
{
  int nchars = renv.length();

  if (0 == nchars)
    return 0;

  nchars--;     // so we can always check the next character

  int rc = 1;

  for (int i = 0; i < nchars; i++)
  {
    if (renv[i + 1] != renv[i])
      continue;

    if (is_operator_character(renv[i]))
      rc++;
  }

  return rc;
}

/*
  examine RENV and pick out an optional operator at the front, and the smarts that follows
*/

static int
get_next_token(const const_IWSubstring & renv,
               int & renv_ndx,
               int & zop,
               IWString & token)
{
  zop = 0;
  token.resize_keep_storage(0);

  if (is_logexp_operator(renv, renv_ndx, zop))
    renv_ndx += 2;

  while (renv_ndx < renv.length())
  {
    int notused;

    if (is_logexp_operator(renv, renv_ndx, notused))
      return token.length();

    token.add(renv[renv_ndx]);
    renv_ndx++;
  }

  return token.length();
}

/*
  Strip leading digits from S and form the number in C
*/

static int
consume_digits(IWString & s,
               int & c)
{
  c = 0;

  while (s.length())
  {
    char c0 = s[0];
    
    int tmp = c0 - '0';
    if (tmp < 0 || tmp > 9)
      return s.length();

    c = c * 10 + tmp;
    s.remove_leading_chars(1);
  }

  return s.length();
}

static int
discern_leading_numerical_qualifier(Min_Max_Specifier<int> & m,
                                    IWString & token)
{
  char c0 = token[0];

  int relational = 0;

  if ('<' == c0)
  {
    relational = -1;
    token.remove_leading_chars(1);
  }
  else if ('>' == c0)
  {
    relational = 1;
    token.remove_leading_chars(1);
  }
  else if (isdigit(c0))
    ;
  else
    return 1;

  int c;
  if (! consume_digits(token, c))
  {
    cerr << "discern_leading_numerical_qualifier: invalid numeric specifier\n";
    return 0;
  }

  if (0 == relational)
    m.add(c);
  else if (relational < 0)
    m.set_max(c - 1);
  else if (relational > 0)
    m.set_min(c + 1);
    
  return 1;
}

// The ring environment is built from Substructure_Atoms, and is not a full Substructure_Query
// and so certain functionality is not present. For now, the only one being checked
// is the ... directive, there may be others.
int
Substructure_Ring_Base::_construct_environment(const const_IWSubstring & renv)
{
  if (renv.contains("...{")) {
    cerr << "Substructure_Ring_Base::_construct_environment:cannot contain no matched atoms directive '" << renv << "'\n";
    return 0;
  }

  int ne = count_components(renv);

  if (0 == ne) {
    cerr << "Substructure_Ring_Base::_construct_environment: no components in '" << renv << "'\n";
    return 0;
  }

  if (renv.ends_with('&') || renv.ends_with('|') || renv.ends_with('^') || renv.ends_with(';'))
  {
    cerr << "Substructure_Ring_Base::_construct_environment: invalid environment '" << renv << "'\n";
    return 0;
  }

  _environment_numerical_requirement.make_room_for_extra_items(ne);
  _environment_atom.make_room_for_extra_items(ne);

  int renv_ndx = 0;

  for (int i = 0; i < ne; i++)
  {
    int zop;
    IWString token;
    if (! get_next_token(renv, renv_ndx, zop, token)) {
      cerr << "Substructure_Atom::_construct_environment: invalid specification, renv_ndx = " << renv_ndx << '\n';
      return 0;
    }

//#define DEBUG_BUILD_ENV
#ifdef DEBUG_BUILD_ENV
    cerr << "i = " << i << " zop '" << zop << "' and smarts '" << token << "'\n";
#endif

    if (0 == i && 0 != zop) {
      cerr << "Substructure_Ring_Base::_construct_environment: first token cannot have an operator\n";
      return 0;
    }

    assert (0 == i || 0 != zop);

    Min_Max_Specifier<int> * mms = new Min_Max_Specifier<int>;

    if (! discern_leading_numerical_qualifier(*mms, token)) {
      cerr << "Cannot discern leading numerical qualifier\n";
      return 0;
    }

    const_IWSubstring smarts(token);

    Substructure_Atom * a = new Substructure_Atom;

    if (! a->parse_smarts_specifier(smarts)) {
      cerr << "Substructure_Ring_Base::_construct_environment: invalid smarts '" << token << "'\n";
      return 0;
    }

    a->count_attributes_specified();    // needed to initialise some things

    _environment_numerical_requirement.add(mms);
    _environment_atom.add(a);

    if (1 == _environment_atom.number_elements()) {  // no operator for the first one
      ;
    } else if (zop) {
      if (! _environment_logexp.add_operator(zop)) {
        cerr << "Substructure_Ring_Base::_construct_environment: huh, operator not recognised '" << zop << "'\n";
        return 0;
      }
    } else {
      _environment_logexp.add_operator('&');
    }
  }

#ifdef DEBUG_BUILD_ENV
  cerr << "Substructure_Ring_Base::_construct_environment: after building, have " << _environment_atom.size() << " components\n";
  _environment_logexp.debug_print(cerr);
#endif

  return 1;
}

// Add =C and =N atoms that are attached to the ring system defined by
// `in_system[_set_global_id]`.
int
Substructure_Ring_Base::ExtendToCarbonyl(const Molecule& m,
        int * in_system) const {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (in_system[i] != _set_global_id) {
      continue;
    }
    const Atom * a = m.atomi(i);
    if (a->ncon() == 2) {
      continue;
    }
    const Atom& atom = m.atom(i);
    for (const Bond * b : atom) {
      if (b->is_single_bond()) {
        continue;
      }

      const atom_number_t j = b->other(i);
      if (in_system[j]) {
        continue;
      }

      const Atom& aj = m.atom(j);
      if (aj.ncon() != 1) {
        continue;
      }

      const atomic_number_t zj = aj.atomic_number();
      if (zj == 6) {
        continue;
      } else if (zj == 8 || zj == 7 || zj == 16) {
      } else {
        continue;
      }

      in_system[j] = in_system[i];
      ++rc;
    }
  }

  return rc;
}

int
IdentifySubstituent(const Molecule& m,
                    atom_number_t zatom,
                    atom_number_t previous_atom,
                    int* in_substituent) {
  in_substituent[zatom] = 2;
  int rc = 1;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    if (j == previous_atom) {
      continue;
    }

    if (in_substituent[j] == 1 || in_substituent[j] == 3) {
      return -1;
    }
    if (in_substituent[j] == 2) {
      continue;
    }
    int tmp = IdentifySubstituent(m, j, zatom, in_substituent);
    if (tmp < 0) {
      return -1;
    }
    rc += tmp;
  }

  return rc;
}

int
IdentifySubstituent(const Molecule& m,
                    atom_number_t zatom,
                    int* in_substituent) {
  int rc = 0;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    if (in_substituent[j] == 1) {
      continue;
    }
    int tmp = IdentifySubstituent(m, j, zatom, in_substituent);
    if (tmp < 0) {
      return -1;
    }
    rc += tmp;
  }

  return rc;
}

void
TranslateNumbers(int * storage,
                 int n,
                 int from,
                 int to) {
  for (int i = 0; i < n; ++i) {
    if (storage[i] == from) {
      storage[i] = to;
    }
  }
}

Substituent::Substituent() {
  _set_global_id = -1;
}

// for all atoms for which storage[i] == flag, set the global_id
// attribute for the corresponding atom in `target`.
void
Substituent::SetTargetGlobalIds(int * storage,
                                int flag,
                                Molecule_to_Match& target) const {
  const int matoms = target.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (storage[i] != flag) {
      continue;
    }
    target[i].set_global_id(_set_global_id);
  }
}

void
Substituent::FillMatchedAtomsArray(std::unique_ptr<int[]>& matched_by_global_specs,
                        const int matoms,
                        const int * storage,
                        int flag) const {
  if (! matched_by_global_specs) {
    matched_by_global_specs.reset(new_int(matoms));
  }

  for (int i = 0; i < matoms; ++i) {
    if (storage[i] == flag) {
      matched_by_global_specs.get()[i] = _set_global_id;
    }
  }
}

// We will be doing a substructure search on `destination`, but want
// to restrict the search to only those atoms for which
// subset[i] == flag. For any atom NOT in that subset, invalidate
// the atoms in `destination` so they cannot match.
// Returns the number of atoms suppressed.
int
InvalidateOtherAtoms(const int* subset,
               int flag,
               Molecule_to_Match& destination) {
  const int matoms = destination.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (subset[i] != flag) {
      destination[i].invalidate();
      ++rc;
    }
  }

  return rc;
}

//#define DEBUG_SUBSTITUENT_MATCHES

// A wrapper for MatchesInner, that handles _set_global_id if needed.
// If queries are present, things are more complex. We only want to
// match a subset of the atoms in the molecule.
// to just a subset of the atoms in `target`. The `save_atomic_numbers`
// array will be set if queries are present, and we copy the atomic
// numbers from `target`. Then, before a search, `target` will be poisoned
// and atoms not to be searched will be set to an unreaslistic atomic number.
int
Substituent::Matches(Molecule_to_Match& target, const int * ring_atoms,
                     int * storage,
                     std::unique_ptr<int[]>& matched_by_global_specs) {
  const int rc = MatchesInner(target, ring_atoms, storage);

#ifdef DEBUG_SUBSTITUENT_MATCHES
  cerr << "Substituent::Matches:from Inner " << rc << '\n';
  for (int i = 0; i < target.natoms(); ++i) {
    cerr << " storage[" << i << "] " << storage[i] << '\n';
  }
#endif
  if (! rc) {
    return 0;
  }

#ifdef DEBUG_SUBSTITUENT_MATCHES
  cerr << "Substituent::Matches: global id " << _set_global_id << '\n';
#endif
  if (_set_global_id < 0) {
    return rc;
  }

  // Turn off the ring atoms, only the substituent(s) atoms.
  TranslateNumbers(storage, target.natoms(), 1, 0);

  // All the atoms that were marked as having a successful substituent match.
  SetTargetGlobalIds(storage, 3, target);
  FillMatchedAtomsArray(matched_by_global_specs, target.natoms(), storage, 3);

#ifdef DEBUG_SUBSTITUENT_MATCHES
  cerr << "Substituent::Matches:global ids assigned\n";
  for (int i = 0; i < target.natoms(); ++i) {
    cerr << " atom " << i << " gid " << target[i].global_id() << '\n';
  }
#endif

  return rc;
}

// Does all the work of matching.
// The storage array is initialized to 1 by copying from `ring_atoms`.
// A value of 2 is used when identifying and testing a substituent.
// Once a valid substituent has been identified, the value will be
// adjusted to 3.
int
Substituent::MatchesInner(Molecule_to_Match& target, const int * ring_atoms,
                          int * storage) {
  // A slightly risky cast, should be OK.
  Molecule& m = const_cast<Molecule&>(*target.molecule());
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_atoms[i]) {
      storage[i] = 1;
    } else {
      storage[i] = 0;
    }
  }

#ifdef DEBUG_SUBSTITUENT_MATCHES
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << " storage " << storage[i] << '\n';
  }
#endif

  int substituents_found = 0;

  // for each ring atom, look for a substituent.
  int matches_found = 0;
  for (int i = 0; i < matoms; ++i) {
    if (storage[i] != 1) {
      continue;
    }
    const Atom& a = m.atom(i);
    // ring atoms will not have substituents if they are 2 connected.
    // If we ever process things other than ring atoms, this will need to be changed.
    if (a.ncon() == 2) {
      continue;
    }

    // Turn off atoms left by a previously failed attempt.
    TranslateNumbers(storage, matoms, 2, 0);
    const int atoms_in_substituent = IdentifySubstituent(m, i, storage);
#ifdef DEBUG_SUBSTITUENT_MATCHES
    cerr << "from atom " << i << " " << atoms_in_substituent << " atoms_in_substituent\n";
#endif
    if (atoms_in_substituent <= 0) {
      continue;
    }

    ++substituents_found;

    if (! _natoms.is_set()) {
    } else if (_natoms.matches(atoms_in_substituent)) {
    } else {
      continue;
    }

    if (! _nrings.is_set()) {
    } else if (OkNrings(m, storage, 2)) {
    } else {
      continue;
    }

    if (! _length.is_set()) {
    } else if (OkLength(m, storage, i, 2)) {
    } else {
      continue;
    }

    if (_required.size() > 0 || _disqualifying.size() > 0) {
      int got_required_match = 0;
      int got_rejected_match = 0;
      RunQueries(target, storage, 2, got_required_match, got_rejected_match);
      if (got_rejected_match) {
        continue;
      }
      if (_required.size() > 0 && ! got_required_match) {
        continue;
      }
    }

    TranslateNumbers(storage, matoms, 2, 3);  // Mark as having a successful substituent match.
    ++matches_found;
  }

#ifdef DEBUG_SUBSTITUENT_MATCHES
  cerr << "Examined " << substituents_found << " substituent, found " << matches_found << " matches\n";
#endif
  if (substituents_found == 0) {
    return 0;
  }

  if (_hits_needed.is_set()) {
    return _hits_needed.matches(matches_found);
  }

  return matches_found;
}

//#define DEBUG_OKNRINGS
int
Substituent::OkNrings(Molecule& m, const int* storage, int flag) const {
#ifdef DEBUG_OKNRINGS
  cerr << "Substituent::OkNrings:checking " << m.nrings() << " rings, flag " << flag << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " storage " << storage[i] << '\n';
  }
#endif

  const int nrings = m.nrings();
  int rings_in_substituent = 0;

  for (int i = 0; i < nrings; ++i) {
    const Ring* r = m.ringi(i);
    if (r->any_members_set_in_array(storage, flag)) {
      ++rings_in_substituent;
    }
  }

#ifdef DEBUG_OKNRINGS
  cerr << "rings_in_substituent " << rings_in_substituent << " matches " << _nrings.matches(rings_in_substituent) << '\n';
#endif

  return _nrings.matches(rings_in_substituent);
}

int
Substituent::OkLength(Molecule& m, const int* storage, atom_number_t anchor, int flag) const {
  const int matoms = m.natoms();
  int longest_length = 0;
  for (int i = 0; i < matoms; ++i) {
    if (storage[i] != flag) {
      continue;
    }
    const int d = m.bonds_between(i, anchor);
#ifdef DEBUG_OKLENGTH
    cerr << d << " bonds between anchor " << anchor << " and " << i << '\n';
#endif
    if (d > longest_length) {
      longest_length = d;
    }
  }

#ifdef DEBUG_OKLENGTH
  cerr << longest_length << " longest_length, matches? " << _length.matches(longest_length) << '\n';
#endif
  return _length.matches(longest_length);
}

int
AllMatchedAtomsInStorage(const Substructure_Results& sresults,
                         const int* storage,
                         int flag) {
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    if (e->all_members_set_in_array(storage, flag)) {
      return 1;
    }
  }

  return 0;
}

// Return true if all atoms in `embedding` have their value
// in `storage` == `flag`.
// Could use std::find_if_not, but I find this clearer and 
// more compact.
bool
AllAtomsMatchFlag(const Set_of_Atoms& embedding,
                  const int* storage,
                  int flag) {
  for (atom_number_t a : embedding) {
    if (storage[a] != flag) {
      return false;
    }
  }

  // All atoms in `embeddings` have storage[i] == flag.
  return true;
}

// Among the embeddings in `sresults` is there one that hits just
// atoms for which storage[i] == flag
int
AnEmbeddingAllInRegion(const Substructure_Results& sresults,
                         const int* storage,
                         int flag) {
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    if (AllAtomsMatchFlag(*e, storage, flag)) {
      return 1;
    }
  }

  return 0;
}

//#define DEBUG_RUNQUERIES

// There is a lot to NOT like with how I am doing these subset queries.
// The most logical way of doing this would be to use the whole molecule
// `target` and then just look for matches that are all in the region of
// interest. This does not work because of smarts that have numeric
// qualifiers. 2[CD1] for example. That may not match the whole molecule,
// but it might match a substituent (or inter-ring region). So, we have
// this less than satisfactory subsetting.
int
Substituent::RunQueries(Molecule_to_Match& target,
                       const int * storage,
                       int flag,
                       int& got_required_match,
                       int& got_rejected_match) {
  if (_required.empty() && _disqualifying.empty()) {
    return 1;
  }

#ifdef DEBUG_RUNQUERIES
  cerr << "Substituent::RunQueries: checking " << _required.size() << " required and " << _disqualifying.size() << " disqualifying queries\n";
  for (int i = 0; i < target.natoms(); ++i) {
    if (flag == storage[i]) {
      cerr << i << ' ' << target.molecule()->smarts_equivalent_for_atom(i) << " being searched\n";
    }
  }
#endif
  Molecule mcopy(*target.molecule());
  Molecule_to_Match subset_target(&mcopy);
  InvalidateOtherAtoms(storage, flag, subset_target);

  if (!_required.empty()) {
    for (Substructure_Query* qry : _required) {
      Substructure_Results sresults;
      if (!qry->substructure_search(subset_target, sresults)) {
        continue;
      }
      if (! AllMatchedAtomsInStorage(sresults, storage, flag)) {
        continue;
      }
      ++got_required_match;
      break;
    }
  }

  for (Substructure_Query* qry : _disqualifying) {
    Substructure_Results sresults;
    if (!qry->substructure_search(subset_target, sresults)) {
      continue;
    }
    if (AllMatchedAtomsInStorage(sresults, storage, flag)) {
      return 0;
    }
  }

  return 1;
}

InterRingAtoms::InterRingAtoms() {
  _set_global_id = -1;
}

int
InterRingAtoms::Matches(Molecule_to_Match& target,
                        const int region_number,
                        InterRingRegionData& data,
                        std::unique_ptr<int[]>& matched_by_global_specs) {
#ifdef DEBUG_INTER_RING_ATOMS_MATCHES
  cerr << "InterRingAtoms::Matches check region number " << region_number << '\n';
  for (int i = 0; i < target.natoms(); ++i) {
    cerr << " atom " << i << " region " << data.region[i]  << '\n';
  }
  cerr << "atoms_in_region " << data.atoms_in_region << " ring_connections " << data.ring_connections << '\n';
#endif

  if (_natoms.is_set() && ! _natoms.matches(data.atoms_in_region)) {
    return 0;
  }

  if (_ring_connections.is_set() && 
      ! _ring_connections.matches(data.ring_connections.number_elements())) {
    return 0;
  }

  if (_required_length.size() > 0) {
    if (! OkRequiredLengths(*target.molecule(), data, region_number)) {
      return 0;
    }
  }

  if (_length.is_set() && ! OkLengths(*target.molecule(), data, region_number)) {
    return 0;
  }

  return RunQueries(target, data, region_number);
}

// Given a set of ring connections that define an inter-ring region
// place into `result` all the inter-atom distances between all pairs
// in `ring_connections`.
// Returns result.size();
int
GetInterRingLengths(Molecule& m,
                    const Set_of_Atoms& ring_connections,
                    std::vector<int>& result) {
  const int n = ring_connections.number_elements();
  if (n == 2) {
    result.push_back(m.bonds_between(ring_connections[0], ring_connections[1]));
    return 1;
  }

  for (int i = 0; i < n; ++i) {
    const atom_number_t a1 = ring_connections[i];
    for (int j = i + 1; j < n; ++j) {
      int d = m.bonds_between(a1, ring_connections[j]);
      result.push_back(d);
    }
  }

  return result.size();
}

// For the required distances to match, there must be an exact and
// full match between what is in _required_length and what is found.
int
InterRingAtoms::OkRequiredLengths(Molecule& m,
                          InterRingRegionData& data,
                          int region_number) {
  if (_required_length.empty()) {
    return 1;
  }

  std::vector<int> distances;
  const int number_distances = GetInterRingLengths(m, data.ring_connections, distances);
#ifdef DEBUG_INTER_RING_ATOMS_OK_REQUIRED_LENGTHS
  cerr << "From " << data.ring_connections << " get " << number_distances << " distances, cmp " << _required_length.size() << " dist";
  for (int d : distances) {
    cerr << ' ' << d;
  }
  cerr << '\n';
#endif
  if (number_distances != _required_length.number_elements()) {
    return 0;
  }

  if (number_distances > 1) {
    std::sort(distances.begin(), distances.end());
  }

  for (int i = 0; i < number_distances; ++i) {
    if (_required_length[i] != distances[i]) {
      return 0;
    }
  }

  return 1;
}

int
InterRingAtoms::OkLengths(Molecule& m,
                          InterRingRegionData& data,
                          int region_number) {
  const int n = data.ring_connections.number_elements();

  // Handle the common case of a linker.
  if (n == 2) {
    return _length.matches(m.bonds_between(data.ring_connections[0], data.ring_connections[1]));
  }

  std::vector<int> distances;
  GetInterRingLengths(m, data.ring_connections, distances);
  for (int d : distances) {
    if (_length.matches(d)) {
      return 1;
    }
  }

  return 0;
}

int
InterRingAtoms::MatchesNhits(int nhits) const {
  // Not set, means any number of occurrences is a match.
  if (! _hits_needed.is_set()) {
    return nhits > 0;
  }

  // Zero matches might indicate a positive match.
  return _hits_needed.matches(nhits);
}

// The inter-ring region query is complicated by the fact that we want
// the ring atoms to be included in what gets searched.
// This function sets up a Molecule_to_Match for the subset of atoms
// to be searched, which includes the ring atoms.
int
InterRingAtoms::RunQueries(Molecule_to_Match& target,
                          InterRingRegionData& data,
                          int region_number) {
  if (_required.empty() && _disqualifying.empty()) {
    return 1;
  }

  Molecule mcopy(*target.molecule());
  Molecule_to_Match subset_target(&mcopy);

  // Temporarily add the ring atoms to this subset.
  data.ring_connections.set_vector(data.region, region_number);

  InvalidateOtherAtoms(data.region, region_number, subset_target);

  const int rc = RunQueriesInner(subset_target, data, region_number);

  // Ring atoms back to where they belong.
  data.ring_connections.set_vector(data.region, 1);

  return rc;
}

// Perform the searches against `subset_target`.
// The check for AllMatchedAtomsInStorage is probably redundant since
// the only atoms that can match are those still active in `subset_target`.
// Remove if this becomes a performance issue.
int
InterRingAtoms::RunQueriesInner(Molecule_to_Match& subset_target,
                          InterRingRegionData& data,
                          int region_number) {
  if (!_required.empty()) {
    bool got_match = false;
    for (Substructure_Query* qry : _required) {
      IWString fname("/tmp/query");
      qry->write_msi(fname);
      Substructure_Results sresults;
      if (!qry->substructure_search(subset_target, sresults)) {
        continue;
      }
      if (! AllMatchedAtomsInStorage(sresults, data.region, region_number)) {
        continue;
      }
      got_match = true;
      break;
    }
    if (! got_match) {
      return 0;
    }
  }

  for (Substructure_Query* qry : _disqualifying) {
    Substructure_Results sresults;
    if (!qry->substructure_search(subset_target, sresults)) {
      continue;
    }
    if (AllMatchedAtomsInStorage(sresults, data.region, region_number)) {
      return 0;
    }
  }

  return 1;
}
