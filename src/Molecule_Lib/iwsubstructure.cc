#include <iostream>
#include <iomanip>
#include <memory>
using std::cerr;
using std::endl;

#include "misc.h"

#include "path.h"
#include "substructure.h"
#include "target.h"
#include "misc2.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#define SUBSTRUCTURE_MAGIC_NUMBER -1010

static int _use_fingerprints_for_screening_substructure_searches = 0;

int
use_fingerprints_for_screening_substructure_searches ()
{
  return _use_fingerprints_for_screening_substructure_searches;
}

void
set_use_fingerprints_for_screening_substructure_searches (int s)
{
  _use_fingerprints_for_screening_substructure_searches = s;
}

static unsigned int report_multiple_hits_threshold = 10000;

void 
set_report_multiple_hits_threshold (unsigned int s)
{
  report_multiple_hits_threshold = s;
}

static int file_scope_ignore_chirality_in_smarts_input = 0;

void
set_ignore_chirality_in_smarts_input (int s)
{
  file_scope_ignore_chirality_in_smarts_input = s;
}

int
ignore_chirality_in_smarts_input ()
{
  return file_scope_ignore_chirality_in_smarts_input;
}

static int file_scope_environment_must_match_unmatched_atoms = 1;

void
set_query_environment_must_match_unmatched_atoms(const int s)
{
  file_scope_environment_must_match_unmatched_atoms = s;
}

void
Single_Substructure_Query::_default_values ()
{
  _magic = SUBSTRUCTURE_MAGIC_NUMBER;

  _find_one_embedding_per_start_atom = 0;
  _rings_in_query = 0;
  _find_unique_embeddings_only = 0;
  _embeddings_do_not_overlap = 0;

  _normalise_rc_per_hits_needed = 0;
  _subtract_from_rc = 0;

  _max_matches_to_find = 0;

  _save_matched_atoms = 1;

  _environment_must_match_unmatched_atoms = file_scope_environment_must_match_unmatched_atoms;

  _ncon_ignore_singly_connected = 0;

  _all_hits_in_same_fragment = 0;

  _consistency_count = -1;

  _max_atoms_in_query = 0;
  _min_atoms_in_query = -1;    // must be negative

  _do_not_perceive_symmetry_equivalent_matches = 0;

  _preferences_present = 0;

  _fragment_ids_present = 0;

  _need_to_compute_ring_membership = -1;
  _need_to_compute_aromaticity = -1;

  _matches_before_checking_environment = 0;
  _no_match_to_environment = 0;
  _match_to_environemt_rejection = 0;

  _sort_by_preference_value = 0;

  _embeddings_violating_distance_constraints = 0;

  _rejection = 0;

  _respect_initial_atom_numbering = 0;
  _highest_initial_atom_number = -1;

  _implicit_ring_condition = -1;     // -1 means not specified

  _matched_atoms_to_check_for_hits_too_close = 0;

  _fail_if_embeddings_too_close = 0;

  _only_keep_matches_in_largest_fragment = 0;

  _min_fraction_atoms_matched = static_cast<float>(0.0);
  _max_fraction_atoms_matched = static_cast<float>(0.0);

  _fingerprint = NULL;

  _ring_ids_present = -1;    // -1 means not known
  _fused_system_ids_present = -1;

  _spinach_match_requirements_present = 0;

  _sort_matches_by = 0;

  _global_conditions_present = -1;    // special flag to indicate unknown. Resolved during first use

  _first_root_atom_with_symmetry_group = -1;   // none of our root atoms have symmetry grouping info

  return;
}

Single_Substructure_Query::Single_Substructure_Query () :
                            _max_query_atoms_matched(0)
{
  _default_values();

  return;
}

Single_Substructure_Query::Single_Substructure_Query (const char * comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::Single_Substructure_Query (const IWString & comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::Single_Substructure_Query (const const_IWSubstring & comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::~Single_Substructure_Query ()
{
  assert (ok());

  _magic = 0;

  if (-2 == _rings_in_query)
    cerr << "Deleting already deleted query\n";

  _rings_in_query = -2;    // keep consistent with check above.

  if (NULL != _fingerprint)
    delete _fingerprint;

  return;
}

//#define SHOW_SSQ_OK

int
Single_Substructure_Query::ok () const
{
#ifdef SHOW_SSQ_OK
  cerr << "Checking ok for query '" << _comment << "' & " << hex << this << dec << endl;
  if (SUBSTRUCTURE_MAGIC_NUMBER != _magic)
    iwabort();
#endif

  if (SUBSTRUCTURE_MAGIC_NUMBER != _magic)
    return 0;

  int nr = _root_atoms.number_elements();
  for (int i = 0; i < nr; i++)
  {
    const Substructure_Atom * a = _root_atoms[i];
    if (! a->ok_recursive())
      return 0;
  }

#ifdef SHOW_SSQ_OK
  cerr << "Atoms are OK, check embeddings " << _embeddings.number_elements() << " and " << _query_atoms_matched.number_elements() << endl;
#endif

  if (_max_matches_to_find > 0 && ! _hits_needed.matches(_max_matches_to_find))
  {
    cerr << "Inconsistency on max_matches_to_find " << _max_matches_to_find << endl;
    _hits_needed.debug_print(cerr);
    return 0;
  }

  if (0 == _save_matched_atoms && _find_unique_embeddings_only)
    return 0;

#ifdef SHOW_SSQ_OK
  cerr << "Substructure_Query is OK\n";
#endif

  return 1;
}

int
Single_Substructure_Query::debug_print (std::ostream & os, const IWString & indentation)
{
  assert (os.good());

  assign_unique_numbers();

  os << indentation << "Query";

  if (_comment.length())
    os << " '" << _comment << "'";

  os << " with " << _root_atoms.number_elements() << " root atoms\n";

  if (_environment.number_elements())
    os << indentation << "Query contains " << _environment.number_elements() << " environment members\n";

  if (_ncon.is_set())
    os << indentation << "_ncon " << _ncon << endl;
  if (_natoms.is_set())
    os << indentation << "_natoms " << _natoms << endl;
  if (_nrings.is_set())
    os << indentation << "_nrings " << _nrings << endl;
  if (_aromatic_rings.is_set())
    os << indentation << "_aromatic_rings " << _aromatic_rings << endl;
  if (_non_aromatic_rings.is_set())
    os << indentation << "_non_aromatic_rings " << _non_aromatic_rings << endl;
  if (_fused_rings.is_set())
    os << indentation << "_fused_rings " << _fused_rings << endl;
  if (_strongly_fused_rings.is_set())
    os << indentation << "_strongly_fused_rings " << _strongly_fused_rings << endl;
  if (_isolated_rings.is_set())
    os << indentation << "_isolated_rings " << _isolated_rings << endl;

  if (_max_matches_to_find)
    os << indentation << " max matches to find " << _max_matches_to_find << endl;
  if (_find_one_embedding_per_start_atom)
    os << indentation << " find one embedding per start atom\n";
  if (_embeddings_do_not_overlap)
    os << indentation << " embeddings do not overlap\n";
  if (_find_unique_embeddings_only)
    os << indentation << " find unique embeddings only\n";
  if (_hits_needed.is_set())
    os << indentation << " hits_needed " << _hits_needed << endl;
  if (_all_hits_in_same_fragment)
    os << indentation << " all hits in same fragment\n";
  if (_subtract_from_rc)
    os << indentation << " substract from hits " << _subtract_from_rc << endl;
  if (_preferences_present)
    os << indentation << " preferences present\n";
  if (_fragment_ids_present)
    os << indentation << " fragment ids present\n";
  if (0 == _save_matched_atoms)
    os << indentation << " matched atoms not saved\n";
  if (_only_keep_matches_in_largest_fragment)
    os << indentation << " only keep matches in largest fragment\n";
  if (_chirality.number_elements())
    os << indentation << _chirality.number_elements() << " chirality specifications\n";
  if (_first_root_atom_with_symmetry_group >= 0)
    os << _first_root_atom_with_symmetry_group << " root atom has symmetry group info\n";

  IWString ind = indentation;
  ind += "  ";

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    os << indentation << "Root atom " << i << endl;
    const Substructure_Atom * r = _root_atoms[i];
    r->recursive_debug_print(os, ind);
  }

  int ne = _environment.number_elements();
  if (ne)
    os << indentation << ne << " environment components\n";
  for (int i = 0; i < ne; i++)
  {
    os << indentation << "Environment  " << i << endl;

    const Substructure_Environment * e = _environment[i];

    e->debug_print(os, ind);
  }

  int nr = _environment_rejections.number_elements();
  if (nr)
    os << indentation << nr << " environment rejections\n";
  for (int i = 0; i < nr; i++)
  {
    os << indentation << "Environment Rejection  " << i << endl;

    const Substructure_Environment * e = _environment_rejections[i];

    e->debug_print(os, ind);
  }

  nr = _ring_specification.number_elements();
  if (nr)
  {
    os << ind << nr << " ring specifications\n";
    for (int i = 0; i < nr; i++)
    {
      const Substructure_Ring_Specification * r = _ring_specification[i];
      r->debug_print(os, ind);
    }
  }

  int nrs = _ring_system_specification.number_elements();
  if (nrs)
  {
    os << ind << nrs << " ring system specifications\n";
    for (int i = 0; i < nrs; i++)
    {
      const Substructure_Ring_System_Specification * r = _ring_system_specification[i];
      r->debug_print(os, ind);

    }
  }

  return 1;
}

int
Single_Substructure_Query::terse_details (std::ostream & os, const IWString & indentation) 
{
  assert (os.good());

  assign_unique_numbers();

  if (_comment.length())
    os << indentation << "Comment '" << _comment << "'\n";

  if (_environment.number_elements())
    os << indentation << "Query contains " << _environment.number_elements() << " environment members\n";

  IWString ind = indentation + "  ";

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    const Substructure_Atom * r = _root_atoms[i];

    os << indentation << "  Root atom " << i << ' ';
    r->terse_details(os, ind);
  }

  int ne = _environment.number_elements();
  if (ne)
    os << indentation << ne << " environment components\n";
  for (int i = 0; i < ne; i++)
  {
    os << indentation << "Environment  " << i << endl;

    const Substructure_Environment * e = _environment[i];

    e->terse_details(os, ind);
  }

  int nr = _environment_rejections.number_elements();
  if (nr)
    os << indentation << nr << " environment rejections\n";
  for (int i = 0; i < nr; i++)
  {
    os << indentation << "Environment Rejection " << i << endl;

    const Substructure_Environment * r = _environment_rejections[i];

    r->terse_details(os, ind);
  }

  return 1;
}

int
Single_Substructure_Query::numeric_value (double & d, int ndx) const
{
  if (! _numeric_value.ok_index(ndx))
    return 0;

  d = _numeric_value[ndx];

  return 1;
}

void
Single_Substructure_Query::_collect_all_atoms (extending_resizable_array<Substructure_Atom *> & completed) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    _root_atoms[i]->collect_all_atoms(completed);
  }

  return;
}

void
Single_Substructure_Query::_determine_if_ring_ids_are_present()
{
  _ring_ids_present = 0;

  int n = _root_atoms.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (_root_atoms[i]->ring_ids_present())
    {
      _ring_ids_present = 1;
      return;
    }
  }

  return;
}

void
Single_Substructure_Query::_determine_if_fused_system_ids_are_present()
{
  _fused_system_ids_present = 0;

  int n = _root_atoms.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (_root_atoms[i]->fused_system_ids_present())
    {
      _fused_system_ids_present = 1;
      return;
    }
  }

  return;
}

void
Single_Substructure_Query::_determine_if_spinach_specifications_are_present ()
{
  _spinach_match_requirements_present = 0;

  int n = _root_atoms.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (_root_atoms[i]->spinach_match_specified())
    {
      _spinach_match_requirements_present = 1;
      return;
    }
  }

  return;
}

void
Single_Substructure_Query::_determine_if_symmetry_groups_present()
{
  const int n = _root_atoms.number_elements();

  _first_root_atom_with_symmetry_group = -1;

//cerr << "Single_Substructure_Query::_determine_if_symmetry_groups_present:checking " << n << " root atoms\n";

  for (int i = 0; i < n; ++i)
  {
    if (_root_atoms[i]->symmetry_groups_present())
    {
      _first_root_atom_with_symmetry_group = i;
      break;
    }
  }

//cerr << "Single_Substructure_Query::_determine_if_symmetry_groups_present " << _first_root_atom_with_symmetry_group << endl;

  return;
}

int
Single_Substructure_Query::_compute_attribute_counts ()
{
  int rc = 0;

  int n = _root_atoms.number_elements();
  for (int i = 0; i < n; i++)
  {
    rc += _root_atoms[i]->attributes_specified();
  }

  n = _environment.number_elements();
  for (int i = 0; i < n; i++)
  {
    rc += _environment[i]->attributes_specified();
  }

  n = _environment_rejections.number_elements();
  for (int i = 0; i < n; i++)
  {
    rc += _environment_rejections[i]->attributes_specified();
  }

  return rc;
}

int
Single_Substructure_Query::set_find_one_embedding_per_atom (int i)
{
  _find_one_embedding_per_start_atom = i;

  return 1;
}

int
Single_Substructure_Query::set_find_unique_embeddings_only (int s)
{
  if (0 == _save_matched_atoms && s > 0)
  {
    cerr << "Single_Substructure_Query::set_find_unique_embeddings_only: must save embeddings\n";
    return 0;
  }

  _find_unique_embeddings_only = s;

  return 1;
}

int
Single_Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches (int s)
{
  if (0 == _save_matched_atoms && s)
  {
    cerr << "Single_Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches: must save matched atoms\n";
    return 0;
  }

  _do_not_perceive_symmetry_equivalent_matches = s;

  return 1;
}

int
Single_Substructure_Query::set_save_matched_atoms (int s)
{
  if (_do_not_perceive_symmetry_equivalent_matches && 0 == s)
  {
    cerr << "Single_Substructure_Query::set_save_matched_atoms: must save embeddings to determine symmetry equivalency\n";
    return 0;
  }

  if (_find_unique_embeddings_only && 0 == s)
  {
    cerr << "Single_Substructure_Query::set_save_matched_atoms: must save embeddings to determine uniqueness\n";
    return 0;
  }

  if ((_sort_matches_by || _sort_by_preference_value) && 0 == s)
  {
    cerr << "Single_Substructure_Query::set_save_matched_atoms: must have embeddings to order preferences\n";
    return 0;
  }

  _save_matched_atoms = s;

  return 1;
}

int
Single_Substructure_Query::set_normalise_rc_per_hits_needed (int i)
{
  _normalise_rc_per_hits_needed = i;

  return 1;
}

int
Single_Substructure_Query::set_max_matches_to_find (int m)
{
  assert (m > 0);

  if (_hits_needed.is_set())
    cerr << "Single_Substructure_Query::set_max_matches_to_find:not set due to _hits_needed\n";
  else
    _max_matches_to_find = m;

  return 1;
}

int
Single_Substructure_Query::set_min_matches_to_find (int m)
{
  assert (m > 0);

  _hits_needed.set_min(m);

  return 1;
}

int
Single_Substructure_Query::set_min_atoms_to_match (int s)
{
  _natoms.set_min(s);

  return _natoms.ok();
}

int
Single_Substructure_Query::set_max_atoms_to_match (int s)
{
  _natoms.set_max(s);

  return _natoms.ok();
}

/*
  This looks wrong right now, seems as if it should do a recursive search...
*/

int
Single_Substructure_Query::involves_rings () const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    const Substructure_Atom * r = _root_atoms[i];
    if (r->involves_rings())
      return 1;
  }

// Should put something in here if any of the global conditions specifies a ring

  return 0;
}

void
Single_Substructure_Query::set_ncon (int s)
{
  _ncon.reset();

  _ncon.add(s);

  return;
}

/*
  September 2000. Bob Coner wanted to be able to perceive implicit rings
  in a query match. Note that in general, this is very hard, but if we
  restrict ourselves to linear queries, it should work OK
*/

int
Single_Substructure_Query::_has_implicit_rings (Query_Atoms_Matched & matched_query_atoms,
                            const int * already_matched) const
{
  int nm = matched_query_atoms.number_elements();

  if (nm < 3)     // rings are at least 3
    return 0;

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * ai = matched_query_atoms[i];

    const Atom * mai = ai->current_hold_atom()->atom();

    for (int j = i + 1; j < nm; j++)
    {
      const Substructure_Atom * aj = matched_query_atoms[j];

      if (ai == aj->parent())    // for linear chains this won't happen, but leave the condition here in case we generalise this
        continue;

      atom_number_t mj = aj->current_hold_atom()->atom_number();

      if (! mai->is_bonded_to(mj))
        continue;

      if (! aj->has_ring_closure_bond_to(ai))    // atoms are bonded, but no ring closure bond. Therefore an implicit ring
        return 1;
    }
  }

  return 0;     // no implicit rings found
}

int
Single_Substructure_Query::_fused_system_id_conditions_satisfied (Query_Atoms_Matched & matched_query_atoms) const
{
  int nm = matched_query_atoms.number_elements();

  Molecule * m = matched_query_atoms[0]->current_hold_atom()->m();

//int matoms = m->natoms();

//#define DEBUG_CHECK_FUSED_SYSTEM_IDS
#ifdef DEBUG_CHECK_FUSED_SYSTEM_IDS
  cerr << "Checking " << nm << " matched atoms for fused system ids '" << m->name() << "', " << m->natoms() << " atoms\n";
  for (auto i = 0; i < m->nrings(); ++i)
  {
    const Ring * ri = m->ringi(i);
    cerr << " ring " << i << " " << (*ri) << endl;
  }
#endif

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * sai = matched_query_atoms[i];

    int fsidi = sai->fused_system_id();
    
    if (0 == fsidi)    // not set
      continue;

    atom_number_t ai = sai->current_hold_atom()->atom_number();

#ifdef DEBUG_CHECK_FUSED_SYSTEM_IDS
    cerr << "Atom " << ai << " " << m->atomic_symbol(ai) << " fsid " << fsidi << " in a ring? " << m->is_ring_atom(ai) << endl;
#endif

//  if (! m->is_ring_atom(ai))     feb 2016. Very important change. fsids are useful to specify an atom that is not in the same ring system as another
//    return 0;

    for (int j = i + 1; j < nm; j++)
    {
      const Substructure_Atom * saj = matched_query_atoms[j];

      int fsidj = saj->fused_system_id();

      if (0 == fsidj)
        continue;

      atom_number_t aj = saj->current_hold_atom()->atom_number();

//    if (! m->is_ring_atom(aj))
//      return 0;

#ifdef DEBUG_CHECK_FUSED_SYSTEM_IDS
//    cerr << "atoms " << ai << ' ' << m->atomic_symbol(ai) << " fsidi " << fsidi << " and " << aj << ' ' << m->atomic_symbol(aj) << " fsidj " << fsidj << " checking same ring system\n";
      cerr << "atoms " << ai << ' ' << m->atomic_symbol(ai) << " fsidi " << fsidi << " and " << aj << ' ' << m->atomic_symbol(aj) << " fsidj " << fsidj << " same ring system " << m->in_same_ring_system(ai, aj) << endl;
#endif

      if (fsidi == fsidj)   // must be in the same ring system
      {
        if (! m->in_same_ring_system(ai, aj))
          return 0;
      }
      else
      {
        if (m->in_same_ring_system(ai, aj))
          return 0;
      }
    }
  }

#ifdef DEBUG_CHECK_RING_IDS
  cerr << "Ring IDS match\n";
#endif

  return 1;
}

int
Single_Substructure_Query::_ring_id_conditions_satisfied (Query_Atoms_Matched & matched_query_atoms) const
{
  int nm = matched_query_atoms.number_elements();

  Molecule * m = matched_query_atoms[0]->current_hold_atom()->m();

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * sai = matched_query_atoms[i];

    int ridi = sai->ring_id();

    if (0 == ridi)
      continue;

    atom_number_t ai = sai->current_hold_atom()->atom_number();

    for (int j = i + 1; j < nm; j++)
    {
      const Substructure_Atom * saj = matched_query_atoms[j];

      int ridj = saj->ring_id();

      if (0 == ridj)
        continue;

      atom_number_t aj = saj->current_hold_atom()->atom_number();

#ifdef DEBUG_CHECK_RING_IDS
      cerr << "atoms " << ai << " ridi " << ridi << " and " << aj << " ridj " << ridj << " same ring " << m->in_same_ring(ai, aj) << endl;
#endif

      if (ridi == ridj)   // must be in the same ring
      {
        if (! m->in_same_ring(ai, aj))
          return 0;
      }
      else
      {
        if (m->in_same_ring(ai, aj))
          return 0;
      }
    }
  }

#ifdef DEBUG_CHECK_RING_IDS
  cerr << "Ring IDS match\n";
#endif

  return 1;
}

int
Single_Substructure_Query::_spinach_match_requirements_satisfied (Query_Atoms_Matched & matched_query_atoms,
                                        Molecule_to_Match & target_molecule) const
{
  int nm = matched_query_atoms.number_elements();

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * sai = matched_query_atoms[i];

    int s = sai->spinach_match_value();
    if (s < 0)    // not specified
      continue;

    Target_Atom * a = sai->current_hold_atom();

    if (NULL == a)
      continue;

    atom_number_t matched_to = a->atom_number();

//#define DEBUG_SPINACH_MATCH_REQUIREMENTS_SATISFIED
#ifdef DEBUG_SPINACH_MATCH_REQUIREMENTS_SATISFIED
    cerr << "Query atom " << sai->unique_id() << " matched with atom " << matched_to << " spinach? " << target_molecule.is_spinach(matched_to) << endl;
#endif

    if (target_molecule.is_spinach(matched_to))
    {
      if (0 == s)
        return 0;
    }
    else    // matched atom not in spinach
    {
      if (s > 0)
        return 0;
    }
  }

  return 1;
}

//#define DEBUG_FRAGMENT_ID_CONDITIONS_SATISFIED

int
Single_Substructure_Query::_fragment_id_conditions_satisfied (Query_Atoms_Matched & matched_query_atoms) const
{
  int nm = matched_query_atoms.number_elements();

  Molecule * m = matched_query_atoms[0]->current_hold_atom()->m();

#ifdef DEBUG_FRAGMENT_ID_CONDITIONS_SATISFIED
  cerr << "_fragment_id_conditions_satisfied:examining " << nm << " matched query atoms\n";
#endif

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * ai = matched_query_atoms[i];

    int fi;
    if (! ai->fragment_id(fi))
      continue;

    int mfi = m->fragment_membership(ai->current_hold_atom()->atom_number());

    for (int j = i + 1; j < nm; j++)
    {
      const Substructure_Atom * aj = matched_query_atoms[j];

      int fj;
      if (! aj->fragment_id(fj))
        continue;

      int mfj = m->fragment_membership(aj->current_hold_atom()->atom_number());

#ifdef DEBUG_FRAGMENT_ID_CONDITIONS_SATISFIED
      cerr << " fragment match conditions i = " << i << " atom " << ai->current_hold_atom()->atom_number() << " mfi " << mfi << " fi " << fi << " j = " << j << " atom " << aj->current_hold_atom()->atom_number() << " mfj " << mfj << " fj " << fj << endl;
#endif

      if (fi == fj && mfi == mfj)
        continue;

      if (fi != fj && mfi != mfj)
        continue;

      return 0;
    }
  }

  return 1;
}

/*
  We have a match - before the environment is checked.
  Check whether or not any specification for ncon met.
*/

int
Single_Substructure_Query::_embedding_ncon_satisfied (Query_Atoms_Matched & matched_query_atoms,
                            const int * already_matched) const
{
  if (! _ncon.is_set())
    return 1;

// Determine the number of attached atoms.

  int attached_atoms = 0;

  int nm = matched_query_atoms.number_elements();
  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * q = matched_query_atoms[i];
    Target_Atom * a = q->current_hold_atom();

    int acon = a->ncon();
    for (int j = 0; j < acon; j++)
    {
      const Target_Atom * aj = a->other(j).other();   // Target_Atom::other() returns a Bond_and_Target_Atom object
      if (already_matched[aj->atom_number()])
        continue;

      if (_ncon_ignore_singly_connected && 1 == aj->ncon())
        continue;

      attached_atoms++;
    }
  }

  return _ncon.matches(attached_atoms);
}

//#define DEBUG_ATTACHED_HETEROATOM_COUNT

/*
  We have a match - before the environment is checked.
  Check whether or not any specification for attached heteroatom count is met.
*/

int
Single_Substructure_Query::_attached_heteroatom_count_satisfied (Query_Atoms_Matched & matched_query_atoms,
                            const int * already_matched) const
{
  if (! _attached_heteroatom_count.is_set())
    return 1;

// Determine the number of attached heteroatoms.

  int ahc = 0;    // the count of attached heteroatoms

  int nm = matched_query_atoms.number_elements();

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
  cerr << "Checking " << nm << " matched query atoms for attached heteroatoms\n";
#endif

  for (int i = 0; i < nm; i++)
  {
    const Substructure_Atom * q = matched_query_atoms[i];
    Target_Atom * a = q->current_hold_atom();

    int acon = a->ncon();

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
    cerr << "Query atom " << q->unique_id() << " matched with atom " << a->atom_number() << " acon = " << acon << endl;
#endif

    for (int j = 0; j < acon; j++)
    {
      const Target_Atom * aj = a->other(j).other();   // Target_Atom::other() returns a Bond_and_Target_Atom object

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
      cerr << "  atom " << aj->atom_number() << " z = " << aj->atomic_number() << " is attached";
      if( already_matched[aj->atom_number()])
        cerr << " already matched";
      cerr << endl;
#endif

      if ( already_matched[aj->atom_number()])
        continue;

      atomic_number_t z = aj->atomic_number();

      if (0 == _heteroatoms.number_elements())
      {
        if (6 != z)
          ahc++;
      }
      else if (_heteroatoms.contains(z))
        ahc++;
    }
  }

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
  cerr << "Final attached heteroatom count is " << ahc << " match is " << _attached_heteroatom_count.matches(ahc) << endl;
#endif

  return _attached_heteroatom_count.matches(ahc);
}

int
Single_Substructure_Query:: _ring_atoms_matched_satisfied (Query_Atoms_Matched & matched_atoms) const
{
  int r = 0;    // the number of ring atoms in the embedding

  for (int i = 0; i < matched_atoms.number_elements(); i++)
  {
    if (matched_atoms[i]->current_hold_atom()->is_ring_atom())
      r++;
  }

  return _ring_atoms_matched.matches(r);
}

int
Single_Substructure_Query:: _heteroatoms_matched_satisfied (Query_Atoms_Matched & matched_atoms) const
{
  int h = 0;    // the number of heteroatoms in the embedding

  for (int i = 0; i < matched_atoms.number_elements(); i++)
  {
    if (6 != matched_atoms[i]->current_hold_atom()->atomic_number())
      h++;
  }

  return _heteroatoms_matched.matches(h);
}

//#define DEBUG_MATCHED_ATOM_POST_CHECK

int
Single_Substructure_Query::_global_query_conditions_also_matched(Query_Atoms_Matched & matched_query_atoms,
                            const int * already_matched,
                            Molecule_to_Match & target_molecule) const
{
#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Checking global conditions\n";
#endif

  if (_root_atoms.number_elements() > 1 && _distance_between_root_atoms.is_set())
  {
    if (! _distance_between_root_atoms_satisfied(matched_query_atoms))
      return 0;
  }

  if (_no_matched_atoms_between.number_elements())
  {
    if (! _no_matched_atoms_between_satisfied(matched_query_atoms))
      return 0;
  }

  if (_link_atom.number_elements())
  {
    if (! _link_atoms_satisfied(matched_query_atoms))
      return 0;
  }

  if (! _attached_heteroatom_count_satisfied(matched_query_atoms, already_matched))
    return 0;

  if (! _embedding_ncon_satisfied(matched_query_atoms, already_matched))
    return 0;

//  cerr << "_fragment_ids_present " << _fragment_ids_present << endl;
  if (_fragment_ids_present && ! _fragment_id_conditions_satisfied(matched_query_atoms))
    return 0;

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Check ring ids? " << _ring_ids_present << endl;
#endif

  if (_ring_ids_present && ! _ring_id_conditions_satisfied(matched_query_atoms))
    return 0;

  if (_fused_system_ids_present && ! _fused_system_id_conditions_satisfied(matched_query_atoms))
    return 0;

  if (_spinach_match_requirements_present && ! _spinach_match_requirements_satisfied(matched_query_atoms, target_molecule))
    return 0;

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Checking " << _element_hits_needed.number_elements() << " element hits needed\n";
#endif

  for (int i = 0; i < _element_hits_needed.number_elements(); i++)
  {
    if (! _element_hits_needed[i]->matches(matched_query_atoms))
      return 0;
  }

  if (_heteroatoms_matched.is_set() && ! _heteroatoms_matched_satisfied(matched_query_atoms))
    return 0;

  if (_ring_atoms_matched.is_set() && ! _ring_atoms_matched_satisfied(matched_query_atoms))
    return 0;

  if (0 == _implicit_ring_condition)
  {
    if (_has_implicit_rings(matched_query_atoms, already_matched))
      return 0;
  }
  else if (1 == _implicit_ring_condition)
  {
    if (! _has_implicit_rings(matched_query_atoms, already_matched))
      return 0;
  }

  if (0 == matched_query_atoms.number_elements())   // don't know how to handle the unmatched atoms stuff
    return 1;

  if (_unmatched_atoms.is_set())
  {
    const Molecule * m = matched_query_atoms[0]->current_hold_atom()->m();

    int unma = m->natoms() - matched_query_atoms.number_elements();

    if (! _unmatched_atoms.matches(unma))
      return 0;
  }

  if (static_cast<float>(0.0) != _min_fraction_atoms_matched ||
      static_cast<float>(0.0) != _max_fraction_atoms_matched)
  {
    const Molecule * m = matched_query_atoms[0]->current_hold_atom()->m();

    float f = static_cast<float>(matched_query_atoms.number_elements()) / static_cast<float>(m->natoms());

//  cerr << "Fraction atoms matched " << f << endl;

    if (static_cast<float>(0.0) == _min_fraction_atoms_matched)
      ;
    else if (f < _min_fraction_atoms_matched)
      return 0;
      
    if(static_cast<float>(0.0) == _max_fraction_atoms_matched)
      ;
    else if (f > _max_fraction_atoms_matched)
      return 0;
  }

  if (_first_root_atom_with_symmetry_group >= 0)
  {
    if (! _symmetry_group_specifications_matches(matched_query_atoms, target_molecule))
      return 0;
  }

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Global conditions OK\n";
#endif

  return 1;
}

/*
  This will be expensive. We have a list of matched atoms and their associated symmetry group specification.
  first task is to collect all that info from the atoms. Note that this could be done once and stored in the
  Single_Substructure_Query object. but since this is not used a lot, we do not worry about doing that. 
  fix if it ever becomes a problem...
*/
 
int
Single_Substructure_Query::_symmetry_group_specifications_matches(const Query_Atoms_Matched & matched_query_atoms,
                                        Molecule_to_Match & target) const
{
  Molecule & m = *(target.molecule());

  const int * msim = m.symmetry_classes();

  const int n = matched_query_atoms.number_elements();

#ifdef DEBUG_SYMMETRY_GROUP_SPECIFICATIONS_MATCHES
  cerr << "Single_Substructure_Query::_symmetry_group_specifications_matches: " << n << " matched atoms\n";
#endif

  for (int i = 0; i < n; ++i)
  {
    const Substructure_Atom * ai = matched_query_atoms[i];

    const auto si = ai->first_symmetry_group();

    if (si <= 0)
      continue;

     const atom_number_t mi = ai->atom_number_matched();

#ifdef DEBUG_SYMMETRY_GROUP_SPECIFICATIONS_MATCHES
    cerr << " i = " << i << " symmetry group " << si << ", atom " << mi << " symm " << msim[mi] << endl;
#endif

    for (int j = i + 1; j < n; ++j)
    {
      const Substructure_Atom * aj = matched_query_atoms[j];

      const auto sj = aj->first_symmetry_group();

      if (sj <= 0)
        continue;

      const atom_number_t mj = aj->atom_number_matched();

      if (si == sj)
      {
        if (msim[mi] != msim[mj])
          return 0;
      }
      else if (msim[mi] == msim[mj])    // si and sj are different
        return 0;
    }
  }

  return 1;    // no conflicts found
}

static int
number_valid_atoms_in_set (const Set_of_Atoms & s)
{
  int rc = 0;

  int n = s.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (s[i] >= 0)
      rc++;
  }

  return rc;
}

/*
  The embedding will be added to the list. Note that some Substructure_Atom
  objects may be excluded from the embedding, so there may be fewer atoms
  in NEW_EMBEDDING than there are in MATCHED_ATOMS.
*/

int
Single_Substructure_Query::_add_embedding (Query_Atoms_Matched & matched_atoms,
                                    Set_of_Atoms * new_embedding,
                                    Substructure_Results & results)
{
  assert (_root_atoms[0] == matched_atoms[0]);

  if (matched_atoms.number_elements() >= new_embedding->number_elements())
    ;
  else if (number_valid_atoms_in_set(*new_embedding) >= matched_atoms.number_elements())
    ;
  else   // too many atoms in the embedding for the number of query atoms matched!
  {
    cerr << "Adding embedding matched atoms contains " << matched_atoms.number_elements() << " embedding " << new_embedding->number_elements() << endl;
    cerr << (*new_embedding) << endl;
    return 0;
  }

  results.add_embedding(new_embedding, matched_atoms);

  if (_max_matches_to_find > 0 &&
      static_cast<int>(results.hits_found()) >= _max_matches_to_find)
    results.set_complete();

  return 1;
}

//#define DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING

/*
  We are dealing with a multi-root query, and a match has been found
  for the first root atom (and all the atoms attached to it).
*/

int
Single_Substructure_Query::_find_next_root_atom_embedding (Query_Atoms_Matched & matched_atoms,
                                    Molecule_to_Match & target_molecule,
                                    int * already_matched,
                                    Substructure_Results & results)
{
  _iroot++;
  assert (_root_atoms.ok_index(_iroot));

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
  cerr << "_find_next_root_atom_embedding: iroot = " << _iroot << endl;
#endif

  Substructure_Atom * r = _root_atoms[_iroot];

  assert (NULL == r->current_hold_atom());

  int istart, istop;
  if (! r->determine_start_stop(target_molecule, istart, istop))
  {
    _iroot--;
    return 0;
  }

  assert (istart >= 0 && istart <= istop && istop <= target_molecule.natoms());

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
  cerr << "Looking between atoms " << istart << " and " << istop << endl;
#endif

  int rc = 0;
  for (int i = istart; i < istop; i++)
  {
#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Testing atom " << i << " already matched = " << already_matched[i] << endl;
#endif

    if (already_matched[i])
      continue;

    Target_Atom & a = target_molecule[i];

    if (! r->matches(a, already_matched))
      continue;

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Looking for embedding, root atom " << _iroot << " matched atom " << i << endl;
#endif

    int tmp = _find_embedding(target_molecule, a, matched_atoms, already_matched, r, results);

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Found " << tmp << " matches at matched atom " << i << endl;
#endif

    if (0 == tmp)    // no matches found this atom
      continue;

    rc += tmp;

    if (_find_one_embedding_per_start_atom)
      break;

    if (_max_matches_to_find > 0 && rc >= _max_matches_to_find)
      break;
  }

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
  cerr << "_find_next_root_atom_embedding, root " << _iroot << " returning " << rc << endl;
#endif

  _iroot--;

  return rc;
}

//#define DEBUG_GOT_EMBEDDING

/*
  We have a match.
  First task is to check whether or not more root atoms need to be matched.

  We refer to atoms in the target to determine whether or not a new
  embedding is unique.
*/

int
Single_Substructure_Query::_got_embedding(Query_Atoms_Matched & matched_atoms,
                                    Molecule_to_Match & target_molecule,
                                    int * already_matched,
                                    Substructure_Results & results)
{
#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Got embedding, " << matched_atoms.number_elements () << " atoms,  iroot = " << _iroot << " nroot = " << _root_atoms.number_elements() << endl;
#endif

  if (_iroot < _root_atoms.number_elements() - 1)
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "iroot " << _iroot << ", " << matched_atoms.number_elements() << " atoms matched\n";
#endif

    return _find_next_root_atom_embedding(matched_atoms, target_molecule, already_matched, results);
  }

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Processing embedding with " << matched_atoms.number_elements() << " atoms matched: ";
  for (int i = 0; i < matched_atoms.number_elements(); i++)
  {
    cerr << ' ' << matched_atoms[i]->atom_number_matched();
  }
  cerr << endl;
#endif

  if (! _global_query_conditions_also_matched(matched_atoms, already_matched, target_molecule))
    return 0;

  if (_chirality.number_elements() && ! _chiral_atoms_matched(matched_atoms, target_molecule))
    return 0;

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "_got_embedding: match OK before environment\n";
#endif

  _matches_before_checking_environment++;

  if (! _query_environment_also_matched(matched_atoms, target_molecule.natoms()))
    return 0;

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "_got_embedding: environment also matched\n";
#endif

  if (! _save_matched_atoms)
  {
    results.got_embedding();

    if (0 == results.hits_found() % report_multiple_hits_threshold)
      cerr << "Query got " << results.hits_found() << " hits\n";

//  if (results.hits_found() > 40000000)
//    exit(0);

    if (_max_matches_to_find > 0 && results.number_embeddings() >= _max_matches_to_find)
      return _max_matches_to_find;

    return 1;
  }

  Set_of_Atoms * new_embedding = new Set_of_Atoms;
  int na = matched_atoms.number_elements();

// Note that if any of the matched atoms are to be excluded from the embedding, we may
// end up with INVALID_ATOM_NUMBER atoms in the final embedding.

  if (_respect_initial_atom_numbering)
  {
    assert (_highest_initial_atom_number >= 0);
    new_embedding->extend(_highest_initial_atom_number + 1, INVALID_ATOM_NUMBER);
  }
  else
    new_embedding->resize(na);

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Scanning " << na << " matched atoms, _respect_initial_atom_numbering " << _respect_initial_atom_numbering << endl;
#endif

  for (int i = 0; i < na; i++)
  {
    const Substructure_Atom * a = matched_atoms[i];

    if (! a->include_in_embedding())
      continue;

    atom_number_t ma = a->atom_number_matched();
//  assert (! new_embedding->contains(ma));

    if (_respect_initial_atom_numbering)
      new_embedding->seti(a->initial_atom_number(), ma);
    else
      new_embedding->add(ma);

#ifdef DEBUG_GOT_EMBEDDING
//  cerr << "Processed atom " << ma << " initial " << a->initial_atom_number() << endl;
#endif
  }

// Not sure what it would mean if all atoms had been excluded from the embedding
// so for now, let's prohibit that

  assert (new_embedding->number_elements() > 0);

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Got embedding " << *new_embedding << ", find unique? " << _find_unique_embeddings_only << endl;
#endif

  if (_find_unique_embeddings_only && ! results.embedding_is_unique(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding is not unique, deleting it\n";
#endif
    delete new_embedding;
    return 0;
  }

  if (_embeddings_do_not_overlap && results.embedding_overlaps_previous_embeddings(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding overlaps previous embeddings, deleting it\n";
#endif

    delete new_embedding;
    return 0;
  }

  if (_do_not_perceive_symmetry_equivalent_matches && 
       results.embedding_is_symmetry_related(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding symmetry related, discarded\n";
#endif
    delete new_embedding;
    return 0;
  }

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Finally, got an embedding\n";
#endif

  return _add_embedding(matched_atoms, new_embedding, results);
}

int
remove_atoms_with_same_or(Query_Atoms_Matched & matched_atoms,
                          const int istart,
                          const int orid)
{
  int rc = 0;
  int na = matched_atoms.number_elements();
  for (int i = istart; i < na; i++)
  {
    const Substructure_Atom * b = matched_atoms[i];
    if (b->or_id() == orid)
    {
      matched_atoms.remove_item(i);
      i--;
      na--;
      rc++;
    }
  }

  return rc;
}

/*
  ROOT_ATOM has been matched. Let's see if a complete embedding can
  be obtained starting with that initial match

  Note that there is inefficiency in the or removal step, because as
  an atom moves to multiple atoms, each time we try to remove all its
  OR'd atoms.
*/

//#define DEBUG_FIND_EMBEDDING

int
Single_Substructure_Query::_find_embedding(Molecule_to_Match & target_molecule,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results)
{
  _max_query_atoms_matched.extra(matched_atoms.number_elements());
  results.matched_this_many_query_atoms(matched_atoms.number_elements());

// In multi-root queries, there may already be atoms in MATCHED_ATOMS.

  int number_initially_matched = matched_atoms.number_elements();

#ifdef DEBUG_FIND_EMBEDDING
  cerr << "_find_embedding: " << number_initially_matched << " atoms initially matched, have " << results.number_embeddings() << " embeddings\n";
#endif

  matched_atoms.add(root_atom);

  int atom_to_process = number_initially_matched + 1;    // before we do the add - C++ arrays start with 0

  if (0 == root_atom->add_your_children(matched_atoms))    // single atom query
  {
    already_matched[root_atom->atom_number_matched()] = 1;
#ifdef DEBUG_FIND_EMBEDDING
    cerr << "_find_embedding:single atom query, matched, have " << results.number_embeddings() << " embeddings so far\n";
#endif
    return _got_embedding(matched_atoms, target_molecule, already_matched, results);
  }

  int rc = 0;

  while (atom_to_process > number_initially_matched)    // we ignore any which may have previously been matched
  {
    Substructure_Atom * a = (Substructure_Atom *) matched_atoms[atom_to_process];

#ifdef DEBUG_FIND_EMBEDDING
    cerr << "Atom to process = " << atom_to_process << " id " << a->unique_id() << 
            " array contains " << matched_atoms.number_elements() << endl;
    for (int i = 0; i < matched_atoms.number_elements(); i++)
    {
      const Substructure_Atom * a = matched_atoms[i];
      cerr << " " << a->unique_id();
      if (a->is_matched())
        cerr << "(" << a->atom_number_matched() << ")";
    }
    cerr << endl;
    if (NULL == a->parent())
      cerr << "Returning " << rc << endl;
#endif

    if (NULL == a->parent())    // must be a root atom, done
      return rc;

    if (! a->move_to_next_match_from_current_anchor(already_matched, matched_atoms))
    {
#ifdef DEBUG_FIND_EMBEDDING
      cerr << "Move to next failed for atom " << a->unique_id() << endl;
#endif

      a->remove_your_children(matched_atoms, already_matched);
      if (a->or_id() && atom_to_process < matched_atoms.number_elements() - 1 &&
          a->or_id() == matched_atoms[atom_to_process + 1]->or_id())
        matched_atoms.remove_item(atom_to_process);
      else
        atom_to_process--;
    }
    else
    {
      _max_query_atoms_matched.extra(atom_to_process + 1);
      results.matched_this_many_query_atoms(atom_to_process + 1);

#ifdef DEBUG_FIND_EMBEDDING
      cerr << "Move to next match succeeded " << 
              a->unique_id() << "(" << a->atom_number_matched() <<
              "), or = " << a->or_id() <<
              " atom to process = " << atom_to_process << " matched = " << matched_atoms.number_elements() << endl;
#endif
      int orid = a->or_id();
      if (orid)
        remove_atoms_with_same_or(matched_atoms, atom_to_process + 1, orid);

      a->add_your_children(matched_atoms);   // does nothing if already added

      if (atom_to_process == matched_atoms.number_elements() - 1)
      {
#ifdef DEBUG_FIND_EMBEDDING
        cerr << "AlL query atoms matched, calling _got_embedding, B4 have " << results.number_embeddings() << " rc " << rc << endl;
#endif
        if (_got_embedding(matched_atoms, target_molecule, already_matched, results))
        {
          rc++;
  
#ifdef DEBUG_FIND_EMBEDDING
          cerr << "Rc incremented to " << rc << " have " << results.number_embeddings() << " embeddings\n";
          if (results.matching_complete() || _find_one_embedding_per_start_atom)
            cerr << "Returning " << rc << endl;
#endif       
          if (results.matching_complete() || _find_one_embedding_per_start_atom)
            return rc;
        }
      }
      else
        atom_to_process++;
    }
  }


#ifdef DEBUG_FIND_EMBEDDING
  cerr << "Single_Substructure_Query::_find_embedding: returning rc = " << rc << endl;
#endif

  return rc;
}

int
Single_Substructure_Query::_find_embedding (Molecule_to_Match & target_molecule,
                                     Target_Atom & a,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results)
{
  int initial_size = matched_atoms.number_elements();

  root_atom->set_hold(&a, already_matched);

#ifdef DEBUG_FIND_EMBEDDING
  cerr << "Single_Substructure_Query::_find_embedding:atom " << a.atom_number() << " begin have " << results.number_embeddings() << " embdddings\n";
#endif

  int rc = _find_embedding(target_molecule, matched_atoms, already_matched, root_atom, results);

#ifdef DEBUG_FIND_EMBEDDING
  cerr << "Single_Substructure_Query::_find_embedding:atom " << a.atom_number() << " returning " << rc << " have " << results.number_embeddings() << " embdddings\n";
#endif

  root_atom->release_hold(already_matched);

  matched_atoms.resize_keep_storage(initial_size);

  return rc;
}

/*
  Right now, we just check to see if there are rings in the query itself.
  If there are none, we return. We should probably check all our components
  to make sure that none of them specify an nrings value of anything but 0
*/

int
Single_Substructure_Query::_adjust_for_internal_consistency()
{
  Molecule tmp ;
  if (! create_molecule(tmp))
  {
    cerr << "Single_Substructure_Query::_adjust_for_internal_consistency: create molecule failed\n";
    return 0;
  }

  if (0 == tmp.nrings())
    return 1;

  int rc = 1;
//int rc = _root.set_ring_membership(tmp);    IMPLEMENT THIS SOME TIME

  _consistency_count++;

  return rc;
}

//#define DEBUG_SUBSTRUCTURE_QUERY

int
Single_Substructure_Query::_substructure_search (Molecule_to_Match & target_molecule,
                                                 int * already_matched,
                                                 Substructure_Results & results)
{
//int matoms = target_molecule.natoms();

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "Beginning atom matching over " << target_molecule.natoms() << " atoms\n";
#endif

  _iroot = 0;

  Substructure_Atom * r = _root_atoms[_iroot];

  int jstart, jstop;
  if (INVALID_ATOM_NUMBER != target_molecule.start_matching_at())
  {
    jstart = target_molecule.start_matching_at();
    jstop = jstart + 1;
  }
  else if (! r->determine_start_stop(target_molecule, jstart, jstop))
    return 0;

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "Start atoms " << jstart << " and " << jstop << endl;
#endif

//assert (jstart >= 0 && jstop >= jstart && jstop <= matoms);

  Query_Atoms_Matched matched_atoms;

  for (int j = jstart; j < jstop; j++)   // loop over atoms in molecule
  {
    if (already_matched[j])   // will only be the case when embeddings are not allowed to overlap
      continue;

    Target_Atom & target_atom = target_molecule[j];

    if (! r->matches(target_atom, already_matched))
      continue;

    _max_query_atoms_matched.extra(1);
    results.matched_this_many_query_atoms(1);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
    cerr << "Atom " << j << " matches query atom 0, " << results.number_embeddings() << " embeddings so far, root " << _iroot << endl;
#endif

    const int tmp = _find_embedding(target_molecule, target_atom, matched_atoms, already_matched, r, results);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
    if (tmp)
    {
      cerr << "Found embedding starting with atom " << j << ", tmp = " << tmp;
      if (results.number_embeddings())
        cerr << " last embedding " << *results.embedding(results.number_embeddings() - 1);
      cerr << endl;
    }
    else
      cerr << "No embedding found, tmp = " << tmp << endl;
    for (int i = 0; i < target_molecule.natoms(); ++i)
    {
      cerr << " alread matched " << i << ' ' << already_matched[i] << target_molecule.molecule()->smarts_equivalent_for_atom(i) << endl;
    }
#endif

    if (tmp && _all_hits_in_same_fragment)
    {
      int k = target_atom.fragment_membership();
#ifdef DEBUG_SUBSTRUCTURE_QUERY
      cerr << "Hit with atom " << target_atom.atom_number() << " in frag " << k << endl;
#endif
      results.got_hit_in_fragment(k);
    }

    if (results.matching_complete())
      return results.return_code();

    if (0 == _embeddings_do_not_overlap)   // that is, embeddings can overlap..
      set_vector(already_matched, target_molecule.natoms(), 0);
  }

  assert (0 == _iroot);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "_substructure_search: returning " << results.number_embeddings() << endl;
#endif

  return results.return_code();
}

int
Single_Substructure_Query::_match_elements_needed (Molecule_to_Match & target_molecule) const
{
  int ne = _elements_needed.number_elements();

  for (int i = 0; i < ne; i++)
  {
    if (! _elements_needed[i]->matches(target_molecule))
      return 0;
  }

  return 1;
}

int
Single_Substructure_Query::_match_heteroatom_specifications (Molecule_to_Match & target_molecule)
{
  if (0 == _heteroatoms.number_elements())
    return _heteroatoms_in_molecule.matches(target_molecule.heteroatom_count());

// Now this gets ugly (and slow)

  int heteroatoms_in_target = 0;
  int natoms = target_molecule.natoms();

  for (int i = 0; i < natoms; i++)
  {
    atomic_number_t z = target_molecule[i].atomic_number();

    if (_heteroatoms.contains(z))
      heteroatoms_in_target++;
  }

  return  _heteroatoms_in_molecule.matches(heteroatoms_in_target);
}

/*
  Note that each of the ring_specifications could match the same ring.
  Would not be too hard to change if ever that becomes an issue
*/

int
Single_Substructure_Query::_match_ring_system_specifications (Molecule_to_Match & target_molecule)
{
#ifdef DEBUG_MATCH_RING_SYSTEM_SPECIFICATIONS
  cerr << "Single_Substructure_Query::_match_ring_system_specifications:checking " << _ring_system_specification.number_elements() << " ring system specifications\n";
#endif

  int nr = _ring_system_specification.number_elements();

  if (1 == nr)
    return _ring_system_specification[0]->matches(target_molecule);

  _ring_system_specification_logexp.reset();

  for (int i = 0; i < nr; i++)
  {
    if (! _ring_system_specification_logexp.result_needed(i))
      continue;

    int m = _ring_system_specification[i]->matches(target_molecule);

    _ring_system_specification_logexp.set_result(i, m);

    int rc;

    if (_ring_system_specification_logexp.evaluate(rc))
      return rc;
  }

  return 1;
}

/*
  Note that each of the ring_specifications could match the same ring.
  Would not be too hard to change if ever that becomes an issue
*/

int
Single_Substructure_Query::_match_ring_specifications (Molecule_to_Match & target_molecule)
{
  int nr = _ring_specification.number_elements();

  if (1 == nr)
    return _ring_specification[0]->matches(target_molecule);

  _ring_specification_logexp.reset();

  for (int i = 0; i < nr; i++)
  {
    if (! _ring_specification_logexp.result_needed(i))
      continue;

    Substructure_Ring_Specification * ri = _ring_specification[i];

    int m = ri->matches(target_molecule);

    _ring_specification_logexp.set_result(i, m > 0);

    int rc;

    if (_ring_specification_logexp.evaluate(rc))
      return rc;
  }

  return 1;
}

int
Single_Substructure_Query::_match_nrings_specifications (Molecule_to_Match & target_molecule)
{
  return _nrings.matches (target_molecule.nrings());
}

int
Single_Substructure_Query::_match_ring_type_specifications(Molecule_to_Match & target_molecule)
{
  if (_aromatic_rings.is_set() && 
      ! _aromatic_rings.matches(target_molecule.aromatic_ring_count()))
    return 0;

  if (_non_aromatic_rings.is_set() &&
      ! _non_aromatic_rings.matches(target_molecule.non_aromatic_ring_count()))
    return 0;

  if (_isolated_rings.is_set() &&
      ! _isolated_rings.matches(target_molecule.isolated_ring_count()))
    return 0;

  if (_isolated_ring_objects.is_set() &&
      ! _isolated_ring_objects.matches(target_molecule.ring_object_count()))
    return 0;

  if (_fused_rings.is_set() &&
      ! _fused_rings.matches(target_molecule.fused_ring_count()))
    return 0;

  if (_strongly_fused_rings.is_set() &&
      ! _strongly_fused_rings.matches(target_molecule.strongly_fused_ring_count()))
    return 0;

  return 1;
}

int
Single_Substructure_Query::_spinach_atoms_match (Molecule_to_Match & target) const
{
  Molecule & m = *(target.molecule());

  int matoms = m.natoms();

  int * spinach = new int[matoms]; std::unique_ptr<int[]> free_spinach(spinach);

  int ais = m.identify_spinach(spinach);

  if (! _atoms_in_spinach.matches(ais))
    return 0;

  if (! _inter_ring_atoms.is_set())
    return 1;

  int number_inter_ring_atoms = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (spinach[i])    // inter ring atoms are not in the spinach
      continue;

    if (target[i].is_non_ring_atom())
      number_inter_ring_atoms++;
  }

  return _inter_ring_atoms.matches(number_inter_ring_atoms);
}

/*
  If is very common for there to be no global conditions present, so we can often save time by detecting that once
  and then not checking on subsequent calls. Make sure you update this function as new global conditions are added!
*/

int
Single_Substructure_Query::_discern_global_conditions_present ()
{
  if (_natoms.is_set() || _heteroatoms_in_molecule.is_set() || _nrings.is_set() || _number_isotopic_atoms.is_set() ||
      _number_fragments.is_set() || _ring_specification.number_elements() || _ring_system_specification.number_elements() || 
      _elements_needed.number_elements() || _atoms_in_spinach.is_set() || _inter_ring_atoms.is_set() || 
      _aromatic_rings.number_elements() || _non_aromatic_rings.is_set() || _isolated_rings.is_set() || _isolated_ring_objects.is_set() || 
      _fused_rings.is_set() || _strongly_fused_rings.is_set() || _net_formal_charge.is_set() ||
      _first_root_atom_with_symmetry_group >= 0)
  {
    _global_conditions_present = 1;
    return 1;
  }

  return 0;
}

//#define DEBUG_CHECK_GLOBAL_CONDITIONS

/*
  At the beginning of a substructure search, check any global specifications
  about the query - fused system sizes, number of rings, etc...

  Make sure that _nrings is checked first!
*/

int
Single_Substructure_Query::_match_global_specifications (Molecule_to_Match & target_molecule)
{
#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Checking global specifications: _natoms " << _natoms.is_set() << ", matches? " << _natoms.matches(target_molecule.natoms()) << endl;
  cerr << "Target contains " << target_molecule.natoms() << " atoms\n";
#endif

  if (_global_conditions_present < 0)
    _discern_global_conditions_present();

  if (0 == _global_conditions_present)
    return 1;

  if (_natoms.is_set())
  {
    if (! _natoms.matches(target_molecule.natoms()))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check heteroatoms " << _heteroatoms_in_molecule.is_set() << endl;
#endif

  if (_heteroatoms_in_molecule.is_set())
  {
    if (! _match_heteroatom_specifications(target_molecule))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check nrings " << _nrings.is_set() << endl;
#endif

  if (_nrings.is_set())
  {
    if (! _nrings.matches(target_molecule.nrings()))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check isotopes " << _number_isotopic_atoms.is_set() << endl;
#endif

  if (_number_isotopic_atoms.is_set())
  {
    if (! _number_isotopic_atoms.matches(target_molecule.number_isotopic_atoms()))
      return 0;
  }

  if (_number_fragments.is_set())
  {
    if (! _number_fragments.matches(target_molecule.number_fragments()))
      return 0;
  }

//if (_nrings.is_set())
//{
//  if (! _match_nrings_specifications(target_molecule))
//    return 0;
//}

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check ring specification(s) " << _ring_specification.number_elements() << endl;
#endif

  if (_ring_specification.number_elements())
  {
    if (! _match_ring_specifications(target_molecule))
      return 0;
  }

  if (_ring_system_specification.number_elements())
  {
    if (! _match_ring_system_specifications(target_molecule))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check elements needed? " << _elements_needed.number_elements() << endl;
#endif

  if (_elements_needed.number_elements())
  {
    if (! _match_elements_needed(target_molecule))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check ring type specifications\n";
#endif

  if (! _match_ring_type_specifications(target_molecule))
    return 0;

  if (_atoms_in_spinach.is_set() || _inter_ring_atoms.is_set())
  {
    if (! _spinach_atoms_match(target_molecule))
      return 0;
  }

  if (_net_formal_charge.is_set() && ! _net_formal_charge.matches(target_molecule.net_formal_charge()))
    return 0;

  return 1;
}

//#define DEBUG_SUBSTRUCTURE_SEARCH

/*
  All queries come through this point.
  We release all holds on the way out. This is because if we don't,
  some of the atoms may still be pointing to target atoms which
  have been deleted.
*/

int
Single_Substructure_Query::_substructure_search (Molecule_to_Match & target_molecule,
                                                 Substructure_Results & results)
{
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Begin common _substructure_search code, " << _root_atoms.number_elements() << " root atoms\n";
#endif

  if (! _match_global_specifications(target_molecule))
  {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
    cerr << "Global specifications not matched\n";
#endif
    return 0;
  }

  if (0 == _root_atoms.number_elements())
    return 1;

  if (_do_not_perceive_symmetry_equivalent_matches)
    results.set_symmetry_class(target_molecule.molecule()->symmetry_classes());
  else
    results.set_symmetry_class(NULL);

  int * tmp = new_int(target_molecule.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  const int rc = _substructure_search(target_molecule, tmp, results);

  const int nr = _root_atoms.number_elements();
  for (int i = 0; i < nr; i++)
  {
    _root_atoms[i]->recursive_release_hold();
  }

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Return code from _substructure_search " << rc << " results contains " << results.number_embeddings() << " embeddings\n";
#endif

  if (_rejection)
    return ! rc;

  return rc;
}

/*
  The return code from substructure search is complicated by the _hits_needed object,
  and especially when there is just a min_hits_needed specified.
  Let's examine some cases.
    _min_hits_needed      hits           rc
       1                   8             8
       2                   8             4
       3                   8             2
*/

int
Single_Substructure_Query::substructure_search (Molecule_to_Match & target_molecule,
                                                Substructure_Results & results)
{
  assert (target_molecule.ok());

  int matoms = target_molecule.natoms();

  results.initialise(matoms);     // no returns before this.

  if (NULL == _fingerprint)
    ;
  else if (! target_molecule.is_superset(*_fingerprint))
    return 0;

  if (_need_to_compute_aromaticity < 0)    // first time this query has been invoked
  {
    _examine_bond_specifications();

    _compute_attribute_counts();
    min_atoms_in_query();
    if (_max_atoms_in_query <= 0)
      assign_unique_numbers();    // these are sequential, 0-

    _determine_if_ring_ids_are_present();

    _determine_if_spinach_specifications_are_present();

    _determine_if_symmetry_groups_present();

    _determine_if_fused_system_ids_are_present();

    if (_respect_initial_atom_numbering)
    {
      for (int i = 0; i < _root_atoms.number_elements(); i++)
      {
        int tmp = _root_atoms[i]->highest_initial_atom_number();
        if (tmp > _highest_initial_atom_number)
          _highest_initial_atom_number = tmp;
      }
    }
  }

  if (matoms < _min_atoms_in_query)
    return 0;     

// If 2 == aromatic_bonds_lose_kekule_identity(), we need a temporary array
// of bond types. Even if the query doesn't specify aromatic bonds, we
// need to convert aromatic rings so they match single bonds

  bond_type_t * save_bt = NULL;

  if (2 == aromatic_bonds_lose_kekule_identity())
  {
    Molecule * m = target_molecule.molecule();

    save_bt = new bond_type_t[m->nedges()];

    m->get_bond_types(save_bt);

    if (0 == m->set_bond_types_for_isis_aromaticity_matching())
    {
      delete [] save_bt;
      save_bt = NULL;
    }
  }
  else if (_need_to_compute_aromaticity)
  {
    target_molecule.establish_aromatic_bonds();
  }
  else if (_need_to_compute_ring_membership)
    (void) target_molecule.molecule()->ring_membership();

  results.set_save_matched_atoms(_save_matched_atoms);

  int nf = 1;
  if (_all_hits_in_same_fragment)
  {
    nf = target_molecule.molecule()->number_fragments();
    results.size_hits_per_fragment_array(nf);
  }
  else if (_only_keep_matches_in_largest_fragment)
    nf = target_molecule.molecule()->number_fragments();

  int rc = _substructure_search(target_molecule, results);

//target_molecule.debug_print(cerr);

  if (NULL != save_bt)
  {
    target_molecule.molecule()->set_bond_types_no_set_modified(save_bt);
    delete [] save_bt;
  }

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Return code from ss = " << rc;
  if (_hits_needed.is_set())
    cerr << " Hits needed is set. match " << _hits_needed.matches(rc) << endl;
  else
    cerr << " Hits needed not set, result " << rc << endl;
#endif

  if (0 == rc)
  {
    if (_hits_needed.is_set())
      return _hits_needed.matches(0);
    else
      return 0;
  }

//cerr << _preferences_present << " _preferences_present\n";
  if (1 == rc || 0 == _preferences_present)   // if only one embedding, or no preferences, do nothing
    ;
  else if (_sort_by_preference_value)
    (void) results.sort_by_preference_value();
  else
    rc = results.remove_low_preference_hits();

  if (rc > 1 && _sort_matches_by)
    results.sort_matches(target_molecule, _sort_matches_by);

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "rc = " << rc << " _distance_between_hits " << _distance_between_hits.is_set() << endl;
#endif

  if (rc > 1 && _distance_between_hits.is_set())
  {
    int tmp = results.remove_hits_violating_distance(target_molecule, _distance_between_hits, _matched_atoms_to_check_for_hits_too_close);   // the number of hits removed

    if (tmp && _fail_if_embeddings_too_close)
      return 0;

    rc = rc - tmp;
    assert (rc > 0);   // it cannot have removed all the hits
  }

// Do any adjustments for fragment stuff

  if (0 == rc)    // no embeddings to worry about
    ;
  else if (_all_hits_in_same_fragment)
  {
    if (nf > 1 && rc > 1)
    {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
      cerr << "nf = " << nf << " and all hits in same frag\n";
#endif

      const extending_resizable_array<int> & hits_per_fragment = results.hits_per_fragment();

      int found_match = 0;
      for (int i = 0; i < nf; i++)
      {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
        cerr << hits_per_fragment[i] << " hits in fragment " << i << " matches " << _hits_needed.matches(hits_per_fragment[i]) << endl;
#endif
        if (_hits_needed.matches(hits_per_fragment[i]))
        {
          found_match = 1;
          rc = hits_per_fragment[i];
          break;
        }
      }
  
      if (! found_match)
        return 0;
    }
  }
  else if (_only_keep_matches_in_largest_fragment)
  {
    if (rc > 0 && nf > 1)
    {
      results.remove_hits_not_in_largest_fragment(target_molecule);
      rc = results.number_embeddings();
    }
  }

  if (! _hits_needed.matches(rc))
    return 0;

// check on any modifications to the return code.

// This is somewhat arbitrary. By definition, we honour _subtract_from_rc if
// it is specified, and in that case, ignore all other modifications.

  if (_subtract_from_rc)
  {
    rc -= _subtract_from_rc;
    if (rc < 0)
      rc = 0;

    return rc;
  }

// If we are not to modify our return code, return here.

  if (0 == _normalise_rc_per_hits_needed)
    return rc;

// If it specifies more than one values for nhits, just return

  if (_hits_needed.number_elements() > 1)
    return rc;

// If we specified 3 hits and we got 3, then we just return 1

  if (1 == _hits_needed.number_elements())
    return 1;

// If a max value is specified, just return.

  int maxv;
  if (_hits_needed.max(maxv))
    return rc;

// Now comes the complex part of perhaps modifying rc

  int minv;
  if (! _hits_needed.min(minv))     // no minimum specified (Hmmm, due to the logic above, returning here should never happend)
    return rc;

  assert (minv > 0);

  rc /= minv;

  assert (rc > 0);

  return rc;
}

/*
  There is a problem with bonds with specifications for ring membership
  or aromaticity. Whereas if we ask an atom for its ring membership,
  it can compute it, a Bond object cannot do that. Therefore if we
  will be asking for ring bonds or aromatic bonds, we need to ask
  the molecule to do that in advance
*/

int
Single_Substructure_Query::_examine_bond_specifications ()
{
  int nr = _root_atoms.number_elements();

  _need_to_compute_ring_membership = 0;

  for (int i = 0; i < nr; i++)
  {
    Substructure_Atom * r = _root_atoms[i];
    if (r->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
    {
      _need_to_compute_aromaticity = 1;
      break;
    }
  }

  int ne = _environment.number_elements();
  for (int i = 0; i < ne; i++)
  {
    Substructure_Environment * e = _environment[i];
    if (e->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
    {
      _need_to_compute_aromaticity = 1;
      break;
    }
  }

  nr = _environment_rejections.number_elements();
  for (int i = 0; i < nr; i++)
  {
    Substructure_Environment * r = _environment_rejections[i];
    if (r->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
    {
      _need_to_compute_aromaticity = 1;
      break;
    }
  }

  if (_need_to_compute_aromaticity < 0)
    _need_to_compute_aromaticity = 0;

//cerr << "At end of _examine_bond_specifications, ring = " << _need_to_compute_ring_membership <<
//        " arom = " << _need_to_compute_aromaticity << endl;

  return 1;
}

int
Single_Substructure_Query::substructure_search(Molecule * m,
                                               Substructure_Results & results)
{
//assert (ok ());
  assert (OK_MOLECULE(m));

  int matoms = m->natoms();

// If the atom has fewer atoms than the "atoms" in the query, this cannot work

  if (matoms < _min_atoms_in_query)
    return 0;

// If the query contains ring(s), but the molecule contains fewer, there can be no match

  if (_rings_in_query > 0 && m->nrings() < _rings_in_query)
    return 0;

  Molecule_to_Match target(m);

  int rc = substructure_search(target, results);

  return rc;
}

/*
  Thought about doing this, but with the code as it is now, it is just too inefficient.
  In cycling through the atoms in a molecule, we would need to allocate all the
  temporary stuff for each atom. Well, I suppose we could keep track of which
  molecule we last searched, but for now, I'm taking the easy way out. See routine
  matched_atoms
*/

/*int
Single_Substructure_Query::match_atom (Molecule * m, atom_number_t a)
{
  assert (ok());
  assert (OK_ATOM_NUMBER(m, a));

  cerr << "This is not implemented\n";
  iwabort ();
  return 0;
}*/

/*int
Single_Substructure_Query::is_atom_matched (atom_number_t a) const
{
  for (int i = 0; i < _embeddings.number_elements(); i++)
  {
    const Set_of_Atoms * p = _embeddings[i];
    if (p->contains(a))
      return 1;
  }

  return 0;      // not found in any of the embeddings.
}*/

/*
  Does this query specify values for ring sizes anywhere?
*/

/*int
Single_Substructure_Query::_ring_sizes_specified (resizable_array<int> & ring_sizes) const
{
  return _root.ring_sizes_specified(ring_sizes);
}*/

int
Single_Substructure_Query::assign_unique_numbers()
{
  int n = 0;

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    _root_atoms[i]->assign_unique_atom_numbers(n);
  }

  int ne = _environment.number_elements();
  for (int i = 0; i < ne; i++)
  {
    _environment[i]->assign_unique_atom_numbers(n);
  }

  int nr = _environment_rejections.number_elements();
  for (int i = 0; i < nr; i++)
  {
    _environment_rejections[i]->assign_unique_atom_numbers(n);
  }

  _max_atoms_in_query = n;

  return _max_atoms_in_query;   // used for generating a new unused number
}

int
Single_Substructure_Query::unique_numbers_from_initial_atom_numbers ()
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    if (! _root_atoms[i]->unique_id_from_initial_atom_number())
    {
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::max_atoms_in_query ()
{
  assert (ok());

  if (_max_atoms_in_query <= 0)
    assign_unique_numbers();

  return _max_atoms_in_query;
}

int
Single_Substructure_Query::min_atoms_in_query ()
{
  if (_min_atoms_in_query > 0)
    return _min_atoms_in_query;

  _min_atoms_in_query = 0;
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    _min_atoms_in_query += _root_atoms[i]->min_atoms_for_match();
  }

  return _min_atoms_in_query;
}

/*
  For the I'th embedding, how many non-matched atoms are bonded to
  the embedding
*/

/*int
Single_Substructure_Query::connections_to_embedding (const Molecule & m,
                                              int embedding) const
{
  assert (ok ());
  assert (_embeddings.ok_index(embedding));

  const Set_of_Atoms * p = _embeddings[embedding];

  int rc = 0;
  int np = p->number_elements();
  for (int i = 0; i < np; i++)
  {
    atom_number_t j = p->item(i);
    int jcon = m.ncon(j);
    for (int k = 0; k < jcon; k++)
    {
      atom_number_t l = m.other(j, k);
      if (! p->contains(l))
        rc++;
    }
  }

  return rc;
}*/

Query_Atoms_Matched::Query_Atoms_Matched ()
{
  _preference_value = 0;

  return;
}

Query_Atoms_Matched::~Query_Atoms_Matched ()
{
  if (-12345 == _preference_value)
    cerr << "Deleting already deleted Query_Atoms_Matched object\n";

  _preference_value = -12345;
}

int
Single_Substructure_Query::print_environment_matches (std::ostream & os) const
{
  os << "  Environment matches for single query '" << _comment << "'\n";

  int ne = _environment.number_elements();
  int nr = _environment_rejections.number_elements();
  if (0 == ne && 0 == nr)
  {
    os << "  No environment(s)\n";
    return os.good();
  }
  
  os << "  " << _matches_before_checking_environment << " matches before checking environment(s)\n";
  if (ne)
    os << "  " << _no_match_to_environment << " embeddings failed to match the environment(s)\n";
  if (nr)
    os << "  " << _match_to_environemt_rejection << " embeddings matched an environment rejection\n";

  if (_embeddings_violating_distance_constraints)
    os << "  " << _embeddings_violating_distance_constraints << " embeddings violated a distance constraint\n";

  if (_matches_before_checking_environment)
  {
    for (int i = 0; i < ne; i++)
    {
      os << "  environment component " << i << endl;
      _environment[i]->print_environment_matches(os);
    }
  }

  if (_match_to_environemt_rejection)
  {
    for (int i = 0; i < nr; i++)
    {
      os << "  environment rejection " << i << endl;
      _environment_rejections[i]->print_environment_matches(os);
    }
  }

  return os.good();
}

Substructure_Atom *
Single_Substructure_Query::query_atom_with_initial_atom_number (atom_number_t a) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    Substructure_Atom * rc = _root_atoms[i]->query_atom_with_initial_atom_number(a);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

Substructure_Atom *
Single_Substructure_Query::query_atom_with_atom_map_number(atom_number_t a) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    Substructure_Atom * rc = _root_atoms[i]->query_atom_with_atom_map_number(a);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

int
Single_Substructure_Query::highest_initial_atom_number () const
{
  int rc = -1;

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    int h = _root_atoms[i]->highest_initial_atom_number();

    if (h > rc)
      rc = h;
  }

  return rc;
}

int
Single_Substructure_Query::highest_atom_map_number () const
{
  int rc = -1;

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    int h = _root_atoms[i]->highest_atom_map_number();

    if (h > rc)
      rc = h;
  }

  return rc;
}

void
Single_Substructure_Query::identify_atom_numbers (extending_resizable_array<int> & a) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    _root_atoms[i]->identify_atom_numbers(a);
  }

  return;
}

void
Single_Substructure_Query::identify_atom_map_numbers (extending_resizable_array<int> & a) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    _root_atoms[i]->identify_atom_map_numbers(a);
  }

  return;
}

const Substructure_Bond *
Single_Substructure_Query::bond_between_atoms (int a1, int a2) const
{
  const Substructure_Bond * rc;

  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    rc = _root_atoms[i]->bond_between_atoms(a1, a2);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

const Substructure_Bond *
Single_Substructure_Query::bond_between_atom_map_numbers(int a1, int a2) const
{
  const Substructure_Bond * rc = NULL;

  for (int i = 0; i < _root_atoms.number_elements(); ++i)
  {
    rc = _root_atoms[i]->bond_between_atom_map_numbers(a1, a2);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

int
Single_Substructure_Query::query_atom_with_isotope (int iso) const
{
  for (auto i = 0; i < _root_atoms.number_elements(); ++i)
  {
    auto rc = _root_atoms[i]->query_atom_with_isotope(iso);

    if (rc >= 0)
      return rc;
  }

  return -1;
}
void
Single_Substructure_Query::assign_unique_id_from_atom_number_if_set (extending_resizable_array<int> & numbers_in_use)
{
  for (unsigned int i = 0; i < _root_atoms.size(); ++i)
  {
    _root_atoms[i]->assign_unique_id_from_atom_number_if_set(numbers_in_use);
  }

  return;
}

int
Single_Substructure_Query::assign_atom_map_numbers(int & amap)
{
  int rc = 0;

  for (unsigned int i = 0; i < _root_atoms.size(); ++i)
  {
    _root_atoms[i]->assign_atom_map_numbers(amap);
  }

  return rc;
}

int
Single_Substructure_Query::print_connectivity_graph(std::ostream & output) const
{
  output << "Single_Substructure_Query::print_connectivity_graph:has " << _root_atoms.number_elements() << " root atoms\n";
  for (int i = 0; i < _root_atoms.number_elements(); ++i)
  {
    output << "ROOT " << i << endl;
    _root_atoms[i]->print_connectivity_graph(output);
  }

  return 1;
}
