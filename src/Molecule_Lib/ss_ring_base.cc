#include <stdlib.h>
#include <ctype.h>
#include <memory>

#include "misc.h"

#include "substructure.h"
#include "target.h"

Substructure_Ring_Base::Substructure_Ring_Base ()
{
  _match_as_match_or_rejection = 1;

  _all_hits_in_same_fragment = 0;

  _only_keep_matches_in_largest_fragment = 0;

  _is_heteroatom = new_int(HIGHEST_ATOMIC_NUMBER + 1, 1);
  _is_heteroatom[6] = 0;

//  Should we include Hydrogen?

  _environment_can_match_in_ring_atoms = 0;

  return;
}

Substructure_Ring_Base::~Substructure_Ring_Base ()
{
  if (NULL != _is_heteroatom)
    delete [] _is_heteroatom;

  return;
}

//#define DEBUG_ENVIRONMENT_MATCHES

int
Substructure_Ring_Base::_environment_matches (Molecule_to_Match & target,
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
    cerr << "Processing child " << atom_to_process << endl;
#endif

    Substructure_Atom * a = const_cast<Substructure_Atom *>(matched_query_atoms[atom_to_process]);

    if (! a->move_to_next_match_from_current_anchor(previously_matched_atoms, matched_query_atoms))
    {
#ifdef DEBUG_ENVIRONMENT_MATCHES
      cerr << "Move to next failed for query environment atom " << a->unique_id() << endl;
#endif

      a->remove_your_children(matched_query_atoms, previously_matched_atoms);

      atom_to_process--;
      if (atom_to_process < 0)
        return rc;
    }
    else
    {
#ifdef DEBUG_ENVIRONMENT_SEARCH
      cerr << "Move to next match succeeded " << a->unique_id() << "(" << a->atom_number_matched() <<
              "), or = " << a->or_id() << " atom to process = " << atom_to_process << " matched = " << matched_query_atoms.number_elements() << endl;
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
Substructure_Ring_Base::_environment_matches (Molecule_to_Match & target,
                                              Substructure_Atom & root_atom,
                                              const int * ring,
                                              int * already_matched)
{
#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "Start environment match, root atom has " << root_atom.attributes_specified() << " attributes specified\n";
#endif

  Query_Atoms_Matched qam;    // scope here just for efficiency
  qam.resize(20);            // 20 seems pretty large

  int matoms = target.natoms();

  int * copy_ring = new int[matoms]; std::unique_ptr<int[]> free_copy_ring(copy_ring);
  copy_vector(copy_ring, ring, matoms);

  int nhits = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == copy_ring[i])     // not in the ring or already matched
      continue;

    Target_Atom & a = target[i];

    copy_vector(already_matched, copy_ring, matoms);

    already_matched[i] = 0;    // the matches() function needs it this way

    if (! root_atom.matches(a, already_matched))
      continue;

    root_atom.set_hold(&a, already_matched);

    if (qam.number_elements())
      qam.resize_keep_storage(0);

    if (0 == root_atom.add_your_children(qam))    // root atom only, no children
    {
#ifdef DEBUG_ENVIRONMENT_SEARCH
      cerr << "Root atom hit, nhits = " << (nhits + 1) << endl;
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
    cerr << "At atom " << i << " matches = " << tmp << endl;
#endif

    root_atom.recursive_release_hold();

    already_matched[i] = 1;

    if (tmp)
      copy_ring[i] = 0;

    nhits += tmp;
  }

#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "Substructure_Ring_Base::_environment_matches: returning " << nhits << endl;
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
Substructure_Ring_Base::_environment_matches (Molecule_to_Match & target,
                                              int * ring,
                                              int * already_matched)
{
  int ne = _environment_atom.number_elements();

  for (int i = 0; i < ne; i++)
  {
    Substructure_Atom * r = _environment_atom[i];

    int nhits = _environment_matches(target, *r, ring, already_matched);

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    cerr << "Result for atom " << i << " is " << nhits << endl;
#endif

    if (_environment_numerical_requirement[i]->is_set())
      nhits = _environment_numerical_requirement[i]->matches(nhits);

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    if (_environment_numerical_requirement[i]->is_set())
      cerr << "After filtering by numerical requirement result is " << nhits << endl;
#endif

    _environment_logexp.set_result(i, nhits);

    int zresult;

#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
    if (_environment_logexp.evaluate(zresult))
      cerr << "Result available, will return " << zresult << endl;
#endif

    if (_environment_logexp.evaluate(zresult))
      return zresult;
  }

  return 0;
}

/*
  RING is an array over the molecule. Atoms which are in either the single ring,
  or the ring system, will be set.
*/

int
Substructure_Ring_Base::_environment_matches (Molecule_to_Match & target,
                                              int * ring)
{
#ifdef DEBUG_L1_ENVIRONMENT_MATCHES
  cerr << "Starting _environment_matches\n";
#endif

  _environment_logexp.reset();

  if (0 == _environment_atom.number_elements())
    return 1;

  int matoms = target.natoms();

  int * already_matched = new int[matoms]; std::unique_ptr<int[]> free_already_matched(already_matched);

  return _environment_matches(target, ring, already_matched);
}

/*
  renv will look like
  1c-F&&<3C-[OH]...
*/

static int
is_operator_character (char c)
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
is_logexp_operator (const const_IWSubstring & s,
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
count_components (const const_IWSubstring & renv)
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
get_next_token (const const_IWSubstring & renv,
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
consume_digits (IWString & s,
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
discern_leading_numerical_qualifier (Min_Max_Specifier<int> & m,
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

int
Substructure_Ring_Base::_construct_environment (const const_IWSubstring & renv)
{
  int ne = count_components(renv);

  if (0 == ne)
  {
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
    if (! get_next_token(renv, renv_ndx, zop, token))
    {
      cerr << "Substructure_Atom::_construct_environment: invalid specification, renv_ndx = " << renv_ndx << endl;
      return 0;
    }

//#define DEBUG_BUILD_ENV
#ifdef DEBUG_BUILD_ENV
    cerr << "i = " << i << " zop '" << zop << "' and smarts '" << token << "'\n";
#endif

    if (0 == i && 0 != zop)
    {
      cerr << "Substructure_Ring_Base::_construct_environment: first token cannot have an operator\n";
      return 0;
    }

    assert (0 == i || 0 != zop);

    Min_Max_Specifier<int> * mms = new Min_Max_Specifier<int>;

    if (! discern_leading_numerical_qualifier(*mms, token))
    {
      cerr << "Cannot discern leading numerical qualifier\n";
      return 0;
    }

    const_IWSubstring smarts(token);

    Substructure_Atom * a = new Substructure_Atom;

    if (! a->parse_smarts_specifier(smarts))
    {
      cerr << "Substructure_Ring_Base::_construct_environment: invalid smarts '" << token << "'\n";
      return 0;
    }

    a->attributes_specified();    // needed to initialise some things

    _environment_numerical_requirement.add(mms);
    _environment_atom.add(a);

    if (1 == _environment_atom.number_elements())  // no operator for the first one
      ;
    else if (zop)
    {
      if (! _environment_logexp.add_operator(zop))
      {
        cerr << "Substructure_Ring_Base::_construct_environment: huh, operator not recognised '" << zop << "'\n";
        return 0;
      }
    }
    else
      _environment_logexp.add_operator('&');
  }

#ifdef DEBUG_BUILD_ENV
  cerr << "Substructure_Ring_Base::_construct_environment: after building\n";
  _environment_logexp.debug_print(cerr);
#endif

  return 1;
}
