#include <stdlib.h>
#include <memory>
using namespace std;

#include "cmdline.h"
#include "misc.h"

#include "molecule.h"
#include "qry_wstats.h"
#include "charge_assigner.h"
#include "target.h"

/*
  The simple flag is an attempt at efficiency.

  Basically it is an attempt to determine, one time, the matched atom
  within each query that has the numeric value.

  while this sounds good, it might break if any of the queries were OR
  queries, with different matched atoms
*/

Charge_Assigner::Charge_Assigner ()
{
  _verbose = 0;

  _molecules_examined = 0;
  _molecules_changed = 0;
  _negative_charges_assigned = 0;
  _positive_charges_assigned = 0;

  _min_distance_between_charges = 3;

  _positive_element = NULL;
  _positive_isotope = 0;
  _negative_element = NULL;
  _negative_isotope = 0;

  _remove_chiral_centres_from_changed_atoms = 0;

  _preserve_implicit_hydrogen_count = 0;

  _attach_explicit_hydrogen_to_positive_atom = 0;

  _apply_isotopic_labels = 0;
    
  _molecules_receiving_negative_charges = 0;
  _molecules_receiving_positive_charges = 0;

  _overwrite_existing_formal_charges = 0;

// For those queries with just one charged atom (probably all) these can speed things up

  _which_atom = NULL;
  _charge_to_assign = NULL;

  _apply_charges_to_molecule = 1;

  _assigned_charge_multiplier = 0;

  return;
}

Charge_Assigner::~Charge_Assigner ()
{
  DELETE_IF_NOT_NULL(_which_atom);
  DELETE_IF_NOT_NULL(_charge_to_assign);

  _molecules_examined = -5;

  return;
}

int
Charge_Assigner::report (ostream & os) const
{
  os << "Report on Charge_Assigner with " << _number_elements << " queries\n";
  os << "Changed " << _molecules_changed << " of " << _molecules_examined << " molecules examined\n";
  if (0 == _molecules_examined)
    return os.good();

  os << "Assigned " << _negative_charges_assigned << " negative charges to " << _molecules_receiving_negative_charges << " molecules\n";
  os << "Assigned " << _positive_charges_assigned << " positive charges to " << _molecules_receiving_positive_charges << " molecules\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const Substructure_Hit_Statistics * q = _things[i];
    os << "Charge Assigner Query " << i << ": ";
    q->report(os, _verbose);
  }

  return os.good();
}

void
display_standard_charge_assigner_options (ostream & os,
                                          char cflag)
{
  os << "  -" << cflag << " <...>       Charge assigner specification, enter \"-" << cflag << " help\" for usage\n";

  return;
}

void
display_all_charge_assigner_options (ostream & os,
                                          char cflag)
{
  os << "  -" << cflag << "             charge assigner flags\n";
  os << "  -" << cflag << " over        allow overwriting existing charges\n";
  os << "  -" << cflag << " <file>      read query from <file>\n";
  os << "  -" << cflag << " <F:file>    read queries from <file>\n";
  os << "  -" << cflag << " <P:symbol>  change positive atoms to element <symbol>\n";
  os << "  -" << cflag << " <N:symbol>  change negative atoms to element <symbol>\n";
  os << "  -" << cflag << " isotope     when changing atom types, isotopically label them\n";
  os << "  -" << cflag << " rmchiral    when changing atom types, remove chiral centres\n";
  os << "  -" << cflag << " hcount      when changing atom types, preserve implicit hydrogens\n";
  os << "  -" << cflag << " nodetach    don't temporarily remove explicit hydrogens\n";
  os << "  -" << cflag << " noremove    don't delete explicit Hydrogens no longer needed\n";
  os << "  -" << cflag << " athpos      attach an explicit Hydrogen to positive atoms\n";
  os << "  -" << cflag << " minsep=<n>  minimum separation between assigned charges\n";

  return;
}

int
Charge_Assigner::construct_from_command_line (Command_Line & cl,
                      int verbose,
                      char cflag)
{
  _verbose = verbose;

  resizable_array_p<Substructure_Hit_Statistics> & tmp = *this;

  int simple_flag = 0;

  const_IWSubstring opt;
  int i = 0;
  while (cl.value(cflag, opt, i++))
  {
    if (_temp_detach_hydrogens.recognised_directive(opt))
      continue;

    if (opt.starts_with("over"))
    {
      _overwrite_existing_formal_charges = 1;
      if (verbose)
        cerr << "Charge_Assigner::construct_from_command_line: existing formal charges can be overwritten\n";
    }
    else if (opt == "simple")
    {
      simple_flag = 1;
    }
    else if (opt.starts_with("P:"))
    {
      opt.remove_leading_chars(2);
      if (NULL == (_positive_element = get_element_from_symbol(opt, _positive_isotope)))
      {
        cerr << "Cannot discern positive element '" << opt << "'\n";
        return 0;
      }
      if (verbose)
        cerr << "Positively charged atoms will be replace by element '" << _positive_element->symbol() << "'\n";
    }
    else if (opt.starts_with("N:"))
    {
      opt.remove_leading_chars(2);
      if (NULL == (_negative_element = get_element_from_symbol(opt, _negative_isotope)))
      {
        cerr << "Cannot discern negative element '" << opt << "'\n";
        return 0;
      }
      if (verbose)
        cerr << "Negatively charged atoms will be replace by element '" << _negative_element->symbol() << "'\n";
    }
    else if ("isotope" == opt)
    {
      _apply_isotopic_labels = 1;
      if (verbose)
        cerr << "Changed atoms will be isotopically labeled\n";
    }
    else if (opt.starts_with("F:"))
    {
      opt.remove_leading_chars(2);
      if (! queries_from_file(opt, tmp, 1, verbose))   // queries always in same directory as control file
      {
        cerr << "Charge_Assigner: cannot read queries from file specifier 'F:" << opt << "'\n";
        return 0;
      }
    }
    else if ("rmchiral" == opt)
    {
      _remove_chiral_centres_from_changed_atoms = 1;
      if (verbose)
        cerr << "Chiral centres on transformed atoms removed\n";
    }
    else if ("hcount" == opt)
    {
      _preserve_implicit_hydrogen_count = 1;
      if (verbose)
        cerr << "Will preserve implicit hydrogens on transformed atoms\n";
    }
    else if ("athpos" == opt)
    {
      _attach_explicit_hydrogen_to_positive_atom = 1;

      if (verbose)
        cerr << "Will attach an explicit Hydrogen atom to atoms assigned positive charges\n";
    }
    else if (opt.starts_with("minsep="))
    {
      opt.remove_leading_chars(7);
      if (! opt.numeric_value(_min_distance_between_charges) || _min_distance_between_charges < 1)
      {
        cerr << "Charge_Assigner::invalid min distance '" << opt << "'\n";
        return 0;
      }
    }
    else if ("help" == opt)
    {
      display_all_charge_assigner_options(cerr, cflag);
      exit(7);    // note exit!
    }
    else 
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics(opt);
      if (! q->read(opt))
      {
        cerr << "Charge_Assigner: cannot read query from '" << opt << "'\n";
        delete q;
        return 0;
      }

      assert (q->ok());

      add(q);

      if (verbose)
        cerr << "Charge_Assigner::construct_from_command_line: made query from '" << opt << "'\n";
    }
  }

  if (0 == _number_elements)
  {
    cerr << "Charge_Assigner::construct_from_command_line: no queries specified\n";
    return 0;
  }

  if (NULL == _negative_element && NULL == _positive_element && _apply_isotopic_labels)
  {
    cerr << "You have asked to apply isotopic labels, but have not specified atoms to change\n";
    cerr << "Use 'P:<symbol>' and/or 'N:<symbol>' together with 'isotope'\n";
    return 0;
  }

  if (simple_flag)
  {
    _which_atom = new_int(_number_elements, -99);    // -99 means not initialised
    _charge_to_assign = new_int(_number_elements, -99);
  }

  return _all_queries_have_numeric_value();
}

void
Charge_Assigner::_increment_global_counters (const Set_of_Atoms & positive_charges_assigned,
                                             const Set_of_Atoms & negative_charges_assigned)
{
  if (positive_charges_assigned.number_elements())
  {
    _molecules_receiving_positive_charges++;
    _positive_charges_assigned += positive_charges_assigned.number_elements();
  }

  if (negative_charges_assigned.number_elements())
  {
    _molecules_receiving_negative_charges++;
    _negative_charges_assigned += negative_charges_assigned.number_elements();
  }

  return;
}

int
Charge_Assigner::_all_queries_have_numeric_value () const
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Substructure_Query * cq = _things[i];   // the (possibly) composite substructure query

    int nq = cq->number_elements();
    for (int j = 0; j < nq; j++)
    {
      const Single_Substructure_Query * q = cq->item(j);
      if (! _numeric_value_present(q))
      {
        cerr << "Charge_Assigner::_all_queries_have_numeric_value: no numeric value found, query " << i << ", component " << j << endl;
        cerr << cq->comment() << endl;
        return 0;
      }
    }
  }

  return _number_elements;
}

int
Charge_Assigner::_numeric_value_present (const Single_Substructure_Query * q) const
{
  int nr = q->root_atoms();
  for (int i = 0; i < nr; i++)
  {
    const Substructure_Atom * a = q->root_atom(i);

    if (_numeric_value_present(a))
      return 1;
  }

  return 0;
}

int
Charge_Assigner::_numeric_value_present (const Substructure_Atom * a) const
{
  double unused;
  if (a->numeric_value(unused))
    return 1;

  int nc = a->number_children();
  for (int i = 0; i < nc; i++)
  {
    Substructure_Atom * child = a->child(i);
    if (_numeric_value_present(child))
      return 1;
  }

  return 0;
}

/*
  common code for putting a specific charge on a given atom
*/

int
Charge_Assigner::_add_charge_to_atom (Molecule_to_Match & target,
                           Set_of_Atoms & positive_charges_assigned,
                           Set_of_Atoms & negative_charges_assigned,
                           formal_charge_t * charges_assigned,
                           atom_number_t zatom,
                           formal_charge_t fc) const
{
  if (charges_assigned[zatom])    // hmm, already assigned a charge to this one
    return 0;

  if (target[zatom].formal_charge() && 0 == _overwrite_existing_formal_charges)
    return 0;

  if (fc > 0)
    positive_charges_assigned.add(zatom);
  else
    negative_charges_assigned.add(zatom);

  charges_assigned[zatom] = fc;

  return 1;
}

static atom_number_t
atom_with_numeric_value (const Query_Atoms_Matched & qami,
                         const Set_of_Atoms * e)
{
  int n = qami.number_elements();

  double tmp;

  for (int i = 0; i < n; i++)
  {
    const Substructure_Atom * a = qami[i];

    if (a->numeric_value(tmp))
      return e->item(i);
  }

  cerr << "No matched atom has a numeric value\n";

  return INVALID_ATOM_NUMBER;
}

/*
  The matches may have preference values associated with them. If we
  have two matches that are close, but one is lower preference value
  than the other, we can remove the lower preference match
*/

int
Charge_Assigner::_remove_lower_preference_hits(Molecule & m,
                                Substructure_Results & sresults) const
{
  int n = sresults.number_embeddings();

  if (1 == n)
    return 1;

//#define DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
#ifdef DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
  cerr << "Removing low priority hits from " << n << " embeddings\n";
#endif

  for (int i = 0; i < n; i++)
  {
    const Query_Atoms_Matched * qami = sresults.query_atoms_matching(i);

    if (0 == qami->preference_value())
      continue;

    atom_number_t ai = atom_with_numeric_value(*qami, sresults.embedding(i));

#ifdef DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
    cerr << " i = " << i << " atom " << ai << endl;
#endif

    if (INVALID_ATOM_NUMBER == ai)    // huh
      continue;

    for (int j = n - 1; j > i; j--)
    {
      const Query_Atoms_Matched * qamj = sresults.query_atoms_matching(j);

      if (0 == qamj->preference_value())
        continue;

      if (qamj->preference_value() == qami->preference_value())
        continue;

#ifdef DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
      cerr << "Preference values " << qami->preference_value() << " and " << qamj->preference_value() << endl;
#endif

      atom_number_t aj = atom_with_numeric_value(*qamj, sresults.embedding(j));

      if (INVALID_ATOM_NUMBER == aj)
        continue;

      if (m.fragment_membership(ai) != m.fragment_membership(aj))
        continue;

      if (ai == aj)
        ;
      else if (m.bonds_between(ai, aj) > _min_distance_between_charges)
        continue;

#ifdef DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
      if (ai != aj)
        cerr << "Atoms " << ai << " and " << aj << " dist " << m.bonds_between(ai, aj) << endl;
      else
        cerr << "Two queries hit atom " << ai << endl;
#endif
      sresults.remove_embedding(j);
      n--;
    }
  }

#ifdef DEBUG_REMOVE_LOWER_PREFERENCE_MATCHES
  cerr << "After removal of lower preference values " << sresults.number_embeddings() << " embeddings\n";
#endif

  return 1;
}

//#define DEBUG_CHARGE_ASSIGNER_PROCESS

int
Charge_Assigner::_process (Molecule_to_Match & target,
                           Set_of_Atoms & positive_charges_assigned,
                           Set_of_Atoms & negative_charges_assigned,
                           formal_charge_t * charges_assigned)
{
  _molecules_examined++;

  int rc = 0;    // we return the number of charges assigned

  for (int i = 0; i < _number_elements; i++)
  {
    Substructure_Results sresults;

    int nhits = _things[i]->substructure_search(target, sresults);
    if (0 == nhits)
      continue;

#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
    cerr << "Charge assigner query " << i << " (" << _things[i]->comment() << ") hits " << nhits << " times\n";
#endif

    if (sresults.number_embeddings() > 1)
      _remove_lower_preference_hits(*(target.molecule()), sresults);

    nhits = sresults.number_embeddings();

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * s = sresults.embedding(j);
#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
      cerr << " hit " << j << " is " << (*s) << endl;
#endif

//    If we already know which member of the embedding carries the charge, just do it

      if (_which_atom && _which_atom[i] >= 0)
      {
        atom_number_t l = s->item(_which_atom[i]);
        if (_add_charge_to_atom(target, positive_charges_assigned, negative_charges_assigned, charges_assigned, l, _charge_to_assign[l]))
          rc++;
        continue;
      }

      const Query_Atoms_Matched * qam = sresults.query_atoms_matching(j);

      assert (s->number_elements() == qam->number_elements());

      int found_non_zero_charge = 0;   // make sure every query has a charge to assign

      int ns = s->number_elements();
      for (int k = 0; k < ns; k++)
      {
        const Substructure_Atom * a = qam->item(k);

        double tmp;
        if (! a->numeric_value(tmp))
          continue;

        formal_charge_t fc = static_cast<formal_charge_t>(tmp);   // double -> int conversion
        if (0 == fc)
          continue;

        if (fc > 0)
          fc = _assigned_charge_multiplier * i + fc;
        else
          fc = - (_assigned_charge_multiplier * i) + fc;

        if (_which_atom && _which_atom[i] < 0)   // not initialised
        {
          _which_atom[i] = k;
          _charge_to_assign[i] = fc;
        }

        found_non_zero_charge++;

        atom_number_t l = s->item(k);

#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
        cerr << "k = " << k << " formal charge of " << fc << " to be placed on atom " << l << ", type " << target[l].atomic_number() << endl;
#endif

        if (_add_charge_to_atom(target, positive_charges_assigned, negative_charges_assigned, charges_assigned, l, fc))
          rc++;
      }

      if (0 == found_non_zero_charge)
        cerr << "Hmmm, query '" << _things[i]->comment() << "' has no non-zero charges\n";
    }
  }

  _increment_global_counters(positive_charges_assigned, negative_charges_assigned);

#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
  cerr << "_process returning " << rc << endl;
#endif

  return rc;
}

void
Charge_Assigner::_remove_positive_charge_hits_on_chiral_atoms (Molecule & m,
                        Set_of_Atoms & positive_charges_assigned,
                        formal_charge_t * charges_assigned) const
{
  for (int i = 0; i < positive_charges_assigned.number_elements(); i++)
  {
    atom_number_t ai = positive_charges_assigned[i];

    if (NULL == m.chiral_centre_at_atom(ai))
      continue;

    if (0 == m.hcount(ai))
      continue;

    positive_charges_assigned.remove_item(i);
    charges_assigned[ai] = 0;
    i--;
  }

  return;
}

static int
advance_shell (const Molecule & m,
               int * tmp,
               int flag)
{
  int matoms = m.natoms();

  int frontier_flag = (flag << 1) | flag;    // two bits

  int expanded = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == (tmp[i] & flag))   // not part of this expansion
      continue;

    if (frontier_flag == (tmp[i] & frontier_flag))  // already part of frontier
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

//  cerr << "Growing frontier from atom " << i << endl;
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (flag & tmp[k])   // already part of this expansion
        continue;

      tmp[k] |= frontier_flag;
      expanded = 1;
//    cerr << "Atom " << k << " added to frontier\n";
    }
  }

  if (0 == expanded)
    return 0;

  frontier_flag = (flag << 1);   // we now need to turn off this bit

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (tmp[i] & frontier_flag)
    {
//    cerr << "Atom " << i << " was on frontier\n";
      tmp[i] ^= frontier_flag;
      rc++;
    }
  }

  return rc;
}

static atom_number_t
identify_survivor (Molecule & m,
                   atom_number_t a1,
                   atom_number_t a2,
                   int * tmp)
{
  set_vector(tmp, m.natoms(), 0);

  tmp[a1] = 1;   // some non overlapping bits
  tmp[a2] = 8;

//cerr << "Checking atoms " << a1 << " and " << a2 << endl;

  if (m.ncon(a1) < m.ncon(a2))
    return a1;
  else if (m.ncon(a1) > m.ncon(a2))
    return a2;

  while (1)
  {
//  cerr << "Advancing shell for atom " << a1 << endl;
    int c1 = advance_shell(m, tmp, 1);
//  cerr << "Advancing shell for atom " << a2 << endl;
    int c2 = advance_shell(m, tmp, 8);
//  cerr << "Shells expanded " << c1 << " and " << c2 << endl;

    if (c1 < c2)
      return a1;
    if (c1 > c2)
      return a2;

//  scores are identical

    if (0 == c1)        // exhausted, return a random choice
      return a1;
  }
}

//#define DEBUG_IDENTIFY_CHARGED_ATOMS_TOO_CLOSE

int
Charge_Assigner::_identify_charged_atoms_too_close (Molecule & m,
                                                    const Set_of_Atoms & s,
                                                    int * times_too_close) const
{
  Set_of_Atoms too_close;

#ifdef DEBUG_IDENTIFY_CHARGED_ATOMS_TOO_CLOSE
  cerr << "Looking for charges too close in " << s.number_elements() << " matches\n";
#endif

  for (int i = 0; i < s.number_elements(); i++)
  {
    atom_number_t ai = s[i];

    int fi = m.fragment_membership(ai);

    for (int j = i + 1; j < s.number_elements(); j++)
    {
      atom_number_t aj = s[j];

      if (m.fragment_membership(aj) != fi)
        continue;

#ifdef DEBUG_IDENTIFY_CHARGED_ATOMS_TOO_CLOSE
      cerr << m.bonds_between(ai, aj) << " bonds between " << ai << " and " << aj << endl;
#endif

      if (ai == aj)    // can this happen?
        ;
      else if (m.bonds_between(ai, aj) > _min_distance_between_charges)
        continue;

#ifdef DEBUG_IDENTIFY_CHARGED_ATOMS_TOO_CLOSE
      cerr << "Atoms " << ai << " '" << m.smarts_equivalent_for_atom(ai) << "' and " << aj << " '" << m.smarts_equivalent_for_atom(aj) << "' too close " << m.bonds_between(ai, aj) << endl;
#endif

      too_close.add(ai);
      too_close.add(aj);
    }
  }

  if (0 == too_close.number_elements())   // no atoms too close
    return 0;

  set_vector(times_too_close, m.natoms(), 0);

  too_close.increment_vector(times_too_close);

  return too_close.number_elements() / 2;
}

void
Charge_Assigner::_remove_hits_too_close_isolation_score (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const
{
  if (s.number_elements() < 2)    // no close matches to eliminate
    return;

  int matoms = m.natoms();

  int * times_too_close = new_int(matoms); std::unique_ptr<int[]> free_times_too_close(times_too_close);

  if (0 == _identify_charged_atoms_too_close(m, s, times_too_close))  // great, nothing too close
    return;

// We have atoms too close. If there is an atom that is in a too-close
// relationship multiple times, remove the charge from it

  int istart = iwmax_of_array(times_too_close, matoms);   // most likely 2 or maybe 3
//cerr << "istart " << istart << endl;

  int charges_removed = 0;

  for (int i = istart; i>= 2; i--)
  {
    int charges_removed_this_loop = 0;

    for (int j = s.number_elements() - 1; j >= 0; j--)
    {
      atom_number_t k = s[j];

      if (times_too_close[k] < i)
        continue;

      charges_assigned[k] = 0;
      s.remove_item(j);
      charges_removed_this_loop++;
    }

    if (0 == charges_removed_this_loop)
      continue;

    charges_removed = 1;
    if (0 == _identify_charged_atoms_too_close(m, s, times_too_close))
      return;
  }

// All we are left with now are pairs, no triples or higher groups

  Set_of_Atoms survivors;

  for (int i = 0; i < s.number_elements(); i++)
  {
    atom_number_t ai = s[i];

    for (int j = i + 1; j < s.number_elements(); j++)
    {
      atom_number_t aj = s[j];

      if (m.fragment_membership(ai) != m.fragment_membership(aj))
        continue;

      if (m.bonds_between(ai, aj) > _min_distance_between_charges)
        continue;

      atom_number_t survivor = identify_survivor(m, ai, aj, times_too_close);

      if (ai == survivor)
        charges_assigned[aj] = 0;
      else
        charges_assigned[ai] = 0;

      survivors.add(survivor);
    }
  }

  s = survivors;

  return;
}

void
Charge_Assigner::_remove_hits_too_close (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const
{
  if (s.number_elements() < 2)    // no close matches to eliminate
    return;

#ifdef DEBUG_REMOVE_HITS_TOO_CLOSE
  cerr << "Atoms are " << s << endl;
#endif

  for (int i = 0; i < s.number_elements(); i++)
  {
    atom_number_t ai = s[i];

    int fi = m.fragment_membership(ai);

    for (int j = i + 1; j < s.number_elements(); j++)
    {
      atom_number_t aj = s[j];

      if (m.fragment_membership(aj) != fi)
        continue;

#ifdef DEBUG_REMOVE_HITS_TOO_CLOSE
      cerr << m.bonds_between(ai, aj) << " bonds between " << ai << " and " << aj << endl;
#endif

      if (ai == aj)
        ;
      else if (m.bonds_between(ai, aj) > _min_distance_between_charges)
        continue;

#ifdef DEBUG_REMOVE_HITS_TOO_CLOSE
      cerr << "Removing atom " << aj << ", '" << m.smarts_equivalent_for_atom(aj) << "'\n";
#endif

      s.remove_item(j);
      charges_assigned[aj] = 0;
      j--;
    }
  }

  return;
}

int
Charge_Assigner::_enumerate_possibilities0 (Molecule & m,
                                          const Set_of_Atoms & s,
                                          resizable_array_p<Set_of_Atoms> & possibilities) const
{
  int n = s.number_elements();

  cerr << "Enumerating " << s << endl;

  if (n < 2)
  {
    Set_of_Atoms * t = new Set_of_Atoms(s);
    possibilities.add(t);

    return 1;
  }

  return _enumerate_possibilities1(m, s, 0, possibilities);
}

int
Charge_Assigner::_enumerate_possibilities1 (Molecule & m,
                                          const Set_of_Atoms & s,
                                          int istart,
                                          resizable_array_p<Set_of_Atoms> & possibilities) const
{
  Set_of_Atoms equivalent;

  int n = s.number_elements();

  for (int i = istart; i < n; i++)
  {
    atom_number_t ai = s[i];

    int fi = m.fragment_membership(ai);

    for (int j = i + 1; j < n; j++)
    {
      atom_number_t aj = s[j];

      if (m.fragment_membership(aj) != fi)
        continue;

      if (m.bonds_between(ai, aj) > _min_distance_between_charges)
        continue;

      equivalent.add_if_not_already_present(ai);
      equivalent.add_if_not_already_present(aj);
    }
  }

  if (0 == equivalent.number_elements())   // just one item
  {
    if (0 == possibilities.number_elements())
    {
      Set_of_Atoms * t = new Set_of_Atoms;
      t->add(s[istart]);
      possibilities.add(t);
    }
    else
    {
      int n = possibilities.number_elements();
      for (int i = 0; i < n; i++)
      {
        Set_of_Atoms * t = new Set_of_Atoms(*(possibilities[i]));
        t->add(s[istart]);
        possibilities.add(t);
      }
    }
  }
  else    // multiple possibilities
  {
    if (0 == possibilities.number_elements())
    {
      int n = equivalent.number_elements();
      for (int i = 0; i < n; i++)
      {
        Set_of_Atoms * t = new Set_of_Atoms;
        t->add(equivalent[i]);
        possibilities.add(t);
      }
    }
    else
    {
      int np = possibilities.number_elements();
      int ne = equivalent.number_elements();

      for (int i = 0; i < ne; i++)
      {
        atom_number_t ei = equivalent[i];
        for (int j = 0; j < np; j++)
        {
          Set_of_Atoms * t = new Set_of_Atoms(*(possibilities[j]));
          t->add(ei);
          possibilities.add(t);
        }
      }
    }
  }

  if (istart < s.number_elements() - 1)
    return _enumerate_possibilities1(m, s, istart + 1, possibilities);

  return possibilities.number_elements();
}

/*
  The HITS argument is optional. If present, it will record the number of
  hits for each query

  Make sure we don't put a +ve charge on any chiral nitrogen

  Jan 2005. Since the queries are prioritised, we need to record the atoms hit
  in the order they are identified
*/

int
Charge_Assigner::_process (Molecule & m,
                           formal_charge_t * charges_assigned)
{
  Set_of_Atoms positive_charges_assigned;
  Set_of_Atoms negative_charges_assigned;

  int matoms = m.natoms();

  Molecule_to_Match target(&m);

  int rc = _process(target, positive_charges_assigned, negative_charges_assigned, charges_assigned);

  if (m.chiral_centres())    // one or more chiral centres present
    _remove_positive_charge_hits_on_chiral_atoms(m, positive_charges_assigned, charges_assigned);

//cerr << "_min_distance_between_charges " << _min_distance_between_charges << ", npos " << positive_charges_assigned.number_elements() << endl;

  if (_min_distance_between_charges > 0)
  {
    _remove_hits_too_close_isolation_score(m, positive_charges_assigned, charges_assigned);
//  _remove_hits_too_close(m, negative_charges_assigned, charges_assigned);  this messes up putting 2 charges on a phosphoric acid
  }

// Should we do anything about positive and negative being too close to each other?

  for (int i = 0; i < matoms; i++)
  {
    if (0 == charges_assigned[i])
      continue;

    int ih;
    if (_preserve_implicit_hydrogen_count && (_positive_element || _negative_element))
      ih = m.implicit_hydrogens(i);
    else
      ih = -1;

    if (charges_assigned[i] > 0)
    {
      if (NULL == _positive_element)      // just apply a charge
      {
        if (_apply_charges_to_molecule)
          m.set_formal_charge(i, charges_assigned[i]);
      }
      else if (_apply_isotopic_labels)
      {
        m.set_element(i, _positive_element);
        m.set_isotope(i, m.atomic_number(i));
      }
      else
      {
        m.set_element(i, _positive_element);
        if (_positive_isotope)
          m.set_isotope(i, _positive_isotope);
      }
    }
    else if (charges_assigned[i] < 0)
    {
      if (NULL == _negative_element)      // just apply a charge
      {
        if (_apply_charges_to_molecule)
          m.set_formal_charge(i, charges_assigned[i]);
      }
      else if (_apply_isotopic_labels)
      {
        m.set_element(i, _negative_element);
        m.set_isotope(i, m.atomic_number(i));
      }
      else
      {
        m.set_element(i, _negative_element);
        if (_negative_isotope)
          m.set_isotope(i, _negative_isotope);
      }
    }

    if (ih >= 0)
      m.set_implicit_hydrogens(i, ih, 1);    // set the "sticky" bit

    if (_remove_chiral_centres_from_changed_atoms && (_positive_element || _negative_element) &&
        NULL != m.chiral_centre_at_atom(i))
      m.remove_chiral_centre_at_atom(i);

    if (_verbose > 1)
      cerr << "Charge_Assigner::_process: atom " << i << " (" << m.atomic_symbol(i) << ") assigned " << charges_assigned[i] << endl;
  }

  if (rc)
    _molecules_changed++;

  return rc;
}

/*
  October 2013. Nasty problems with a 3D molecule that had an explicit hydrogen early in the connection table.
  We use the temp detach hydrogens, process the molecule, which adds a hydrogen atom.
  But when temp detach does the re-attach, if notices that the atom is no longer needed and so drops it.
  That totally messes up the ordering in the charges_assigned array.
  To be safe, we always move explicit hydrogens to the end of the connection table
*/

int
Charge_Assigner::process (Molecule & m,
                          formal_charge_t * charges_assigned)
{
  int matoms = m.natoms();
  if (0 == matoms)
  {
    cerr << "Charge_Assigner::process: cannot process no-struct\n";
    return 0;
  }

  m.move_hydrogens_to_end_of_connection_table();   // temp detach may remove hydrogens that are found to be superfluous, and that will mess up our array

  const int hydrogens_present = _temp_detach_hydrogens.detach_atoms(m);

  int i_need_to_delete_charges_assigned = 0;

// Should we zero the array if it has been passed to us???

  if (NULL == charges_assigned)
  {
    charges_assigned = new_int(matoms);
    i_need_to_delete_charges_assigned = 1;
  }

  int rc = _process(m, charges_assigned);

  if (hydrogens_present)
    _temp_detach_hydrogens.reattach_atoms(m);

#ifdef DEBUG_CHARGE_ASSIGNER
  for (auto i = 0; i < matoms; ++i)
  {
    if (charges_assigned[i] != 0)
      cerr << "Atom " << i << " " << m.smarts_equivalent_for_atom(i) << " assigned " << charges_assigned[i] << endl;
  }
#endif

  if (_attach_explicit_hydrogen_to_positive_atom)
    _do_make_implicit_hydrogens_explicit(m, charges_assigned);

  if (i_need_to_delete_charges_assigned)
    delete [] charges_assigned;

  return rc;
}

int
Charge_Assigner::process (Molecule & m,
                          resizable_array_p<Molecule> & charged_forms)
{
  int matoms = m.natoms();
  if (0 == matoms)
  {
    cerr << "Charge_Assigner::process: cannot process no-struct\n";
    return 0;
  }

  _temp_detach_hydrogens.detach_atoms(m);

  int rc = _process(m, charged_forms);

  _temp_detach_hydrogens.reattach_atoms(m);

  if (charged_forms.number_elements() > 0)
    _molecules_changed++;

  return rc;
}

static void
transfer_charges (Molecule & m,
                  const formal_charge_t * charges_assigned)
{
  int matoms = m.natoms();

  for (int j = 0; j < matoms; j++)
  {
    if (0 != charges_assigned[j])
      m.set_formal_charge(j, charges_assigned[j]);
  }

  return;
}

int
Charge_Assigner::_process (Molecule & m,
                           resizable_array_p<Molecule> & charged_forms)
{
  int matoms = m.natoms();

  formal_charge_t * charges_assigned = new formal_charge_t[matoms]; std::unique_ptr<formal_charge_t[]> free_charges_assigned(charges_assigned);
  set_vector(charges_assigned, matoms, static_cast<formal_charge_t>(0));

  Set_of_Atoms negative_charges_assigned, positive_charges_assigned;

  Molecule_to_Match target(&m);

#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
  cerr << "Processing '" << m.name() << "'\n";
#endif

  if (0 == _process(target, positive_charges_assigned, negative_charges_assigned, charges_assigned))
    return 0;

#ifdef DEBUG_CHARGE_ASSIGNER_PROCESS
  cerr << "To '" << m.name() << " assign " << positive_charges_assigned.number_elements() << " +ve and " << negative_charges_assigned.number_elements() << " -ve charges\n";
#endif

  if (positive_charges_assigned.number_elements() < 2)    // no enumeration possible
  {
    Molecule * t = new Molecule(m);
    t->set_name(m.name());

    negative_charges_assigned.set_vector(charges_assigned, -1);
    transfer_charges(*t, charges_assigned);
    charged_forms.add(t);

    return 1;
  }

  resizable_array_p<Set_of_Atoms> set_of_positive_charges;

  _enumerate_possibilities0(m, positive_charges_assigned, set_of_positive_charges);

  int n = set_of_positive_charges.number_elements();

  for (int i = 0; i < n; i++)
  {
    set_vector(charges_assigned, matoms, static_cast<formal_charge_t>(0));

    const Set_of_Atoms * p = set_of_positive_charges[i];

    p->set_vector(charges_assigned, 1);

    negative_charges_assigned.set_vector(charges_assigned, -1);

    Molecule * t = new Molecule(m);
    t->set_name(m.name());

    transfer_charges(*t, charges_assigned);

    charged_forms.add(t);
  }

  if (positive_charges_assigned.number_elements())
  {
    _positive_charges_assigned += positive_charges_assigned.number_elements();
    _molecules_receiving_positive_charges++;
  }

  if (negative_charges_assigned.number_elements())
  {
    _negative_charges_assigned += negative_charges_assigned.number_elements();
    _molecules_receiving_negative_charges++;
  }

  return charged_forms.number_elements();
}

int
Charge_Assigner::_do_make_implicit_hydrogens_explicit (Molecule & m,
                                                const int * charges_assigned) const
{
  int rc = 0;

  Make_Implicit_Hydrogens_Explicit mihe;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
//  if (0 != charges_assigned[i])
//    cerr << "Atom " << m.smarts_equivalent_for_atom(i) << " has charge " << charges_assigned[i] << ", ih " << m.implicit_hydrogens(i) << endl;

    if (charges_assigned[i] <= 0)
      continue;

    if (1 != m.implicit_hydrogens(i))
      continue;

    mihe.set_atom(i);

    m.make_implicit_hydrogens_explicit(mihe);
    rc++;
  }

  return rc;
}

int
Charge_Assigner::build (const const_IWSubstring & s)
{
  resizable_array_p<Substructure_Hit_Statistics> & tmp = *this;

  const_IWSubstring token;

  int verbose = 0;

  for (auto i = 0; s.nextword(token, i); )
  {
    if (token.starts_with("F:"))
    {
      token.remove_leading_chars(2);
      if (! queries_from_file(token, tmp, 1, verbose))   // queries always in same directory as control file
      {
        cerr << "Charge_Assigner: cannot read queries from file specifier 'F:" << token << "'\n";
        return 0;
      }
    }
    else if ("verbose" == token)
      verbose = 1;
    else
    {
      cerr << "Charge_Assigner::build:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  if (0 == _number_elements)
  {
    cerr << "Charge_Assigner::build:no queries\n";
    return 0;
  }

  return 1;
}
