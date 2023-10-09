#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;

#define SUBSTRUCTURE_MAGIC_NUMBER -1010

static int _use_fingerprints_for_screening_substructure_searches = 0;

int
use_fingerprints_for_screening_substructure_searches()
{
  return _use_fingerprints_for_screening_substructure_searches;
}

void
set_use_fingerprints_for_screening_substructure_searches(int s)
{
  _use_fingerprints_for_screening_substructure_searches = s;
}

static unsigned int report_multiple_hits_threshold = 10000;

void 
set_report_multiple_hits_threshold(unsigned int s)
{
  report_multiple_hits_threshold = s;
}

static int file_scope_ignore_chirality_in_smarts_input = 0;

void
set_ignore_chirality_in_smarts_input(int s)
{
  file_scope_ignore_chirality_in_smarts_input = s;
}

int
ignore_chirality_in_smarts_input()
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
Single_Substructure_Query::_default_values()
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
  _min_all_matches_fraction_atoms_matched = static_cast<float>(0.0);
  _max_all_matches_fraction_atoms_matched = static_cast<float>(0.0);

  _fingerprint = nullptr;

  _ring_ids_present = -1;    // -1 means not known
  _fused_system_ids_present = -1;

  _spinach_match_requirements_present = 0;

  _sort_matches_by = 0;

  _global_conditions_present = -1;    // special flag to indicate unknown. Resolved during first use

  _first_root_atom_with_symmetry_group = -1;   // none of our root atoms have symmetry grouping info

  _no_match_unless_nhits = 0;

  _unmatched_atoms_attached_specified = 0;

  _compress_embeddings = 0;

  _atom_typing = nullptr;

  _atom_type_groups_present = 0;
    
  _no_matched_atoms_between_exhaustive = 1;

  return;
}

Single_Substructure_Query::Single_Substructure_Query() :
                            _max_query_atoms_matched(0)
{
  _default_values();

  return;
}

Single_Substructure_Query::Single_Substructure_Query(const char * comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::Single_Substructure_Query(const IWString & comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::Single_Substructure_Query(const const_IWSubstring & comment) :
                            _max_query_atoms_matched(0)
{
  _default_values();

  _comment = comment;

  return;
}

Single_Substructure_Query::~Single_Substructure_Query()
{
  assert(ok());

  _magic = 0;

  if (-2 == _rings_in_query)
    cerr << "Deleting already deleted query\n";

  _rings_in_query = -2;    // keep consistent with check above.

  if (nullptr != _fingerprint)
    delete _fingerprint;

  if (_atom_typing != nullptr)
    delete _atom_typing;

  return;
}

//#define SHOW_SSQ_OK

int
Single_Substructure_Query::ok() const
{
#ifdef SHOW_SSQ_OK
  cerr << "Checking ok for query '" << _comment << "' & " << hex << this << dec << '\n';
  if (SUBSTRUCTURE_MAGIC_NUMBER != _magic)
    iwabort();
#endif

  if (SUBSTRUCTURE_MAGIC_NUMBER != _magic)
    return 0;

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (! r->ok_recursive())
      return 0;
  }

#ifdef SHOW_SSQ_OK
  cerr << "Atoms are OK, check embeddings " << _embeddings.number_elements() << " and " << _query_atoms_matched.number_elements() << '\n';
#endif

  if (_max_matches_to_find > 0 && ! _hits_needed.matches(_max_matches_to_find))
  {
    cerr << "Inconsistency on max_matches_to_find " << _max_matches_to_find << '\n';
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
Single_Substructure_Query::debug_print(std::ostream & os, const IWString & indentation)
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
    os << indentation << "_ncon " << _ncon << '\n';
  if (_natoms.is_set())
    os << indentation << "_natoms " << _natoms << '\n';
  if (_nrings.is_set())
    os << indentation << "_nrings " << _nrings << '\n';
  if (_aromatic_rings.is_set())
    os << indentation << "_aromatic_rings " << _aromatic_rings << '\n';
  if (_non_aromatic_rings.is_set())
    os << indentation << "_non_aromatic_rings " << _non_aromatic_rings << '\n';
  if (_fused_rings.is_set())
    os << indentation << "_fused_rings " << _fused_rings << '\n';
  if (_strongly_fused_rings.is_set())
    os << indentation << "_strongly_fused_rings " << _strongly_fused_rings << '\n';
  if (_isolated_rings.is_set())
    os << indentation << "_isolated_rings " << _isolated_rings << '\n';

  if (_max_matches_to_find)
    os << indentation << " max matches to find " << _max_matches_to_find << '\n';
  if (_find_one_embedding_per_start_atom)
    os << indentation << " find one embedding per start atom\n";
  if (_embeddings_do_not_overlap)
    os << indentation << " embeddings do not overlap\n";
  if (_find_unique_embeddings_only)
    os << indentation << " find unique embeddings only\n";
  if (_hits_needed.is_set())
    os << indentation << " hits_needed " << _hits_needed << '\n';
  if (_all_hits_in_same_fragment)
    os << indentation << " all hits in same fragment\n";
  if (_subtract_from_rc)
    os << indentation << " substract from hits " << _subtract_from_rc << '\n';
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
    os << indentation << "Root atom " << i << '\n';
    const Substructure_Atom * r = _root_atoms[i];
    r->recursive_debug_print(os, ind);
  }

  int ne = _environment.number_elements();
  if (ne)
    os << indentation << ne << " environment components\n";
  for (int i = 0; i < ne; i++)
  {
    os << indentation << "Environment  " << i << '\n';

    const Substructure_Environment * e = _environment[i];

    e->debug_print(os, ind);
  }

  int nr = _environment_rejections.number_elements();
  if (nr)
    os << indentation << nr << " environment rejections\n";
  for (int i = 0; i < nr; i++)
  {
    os << indentation << "Environment Rejection  " << i << '\n';

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
Single_Substructure_Query::terse_details(std::ostream & os, const IWString & indentation) 
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
    os << indentation << "Environment  " << i << '\n';

    const Substructure_Environment * e = _environment[i];

    e->terse_details(os, ind);
  }

  int nr = _environment_rejections.number_elements();
  if (nr)
    os << indentation << nr << " environment rejections\n";
  for (int i = 0; i < nr; i++)
  {
    os << indentation << "Environment Rejection " << i << '\n';

    const Substructure_Environment * r = _environment_rejections[i];

    r->terse_details(os, ind);
  }

  return 1;
}

int
Single_Substructure_Query::numeric_value(double & d, int ndx) const
{
  if (! _numeric_value.ok_index(ndx))
    return 0;

  d = _numeric_value[ndx];

  return 1;
}

void
Single_Substructure_Query::_collect_all_atoms(extending_resizable_array<Substructure_Atom *> & completed) const
{
  for (Substructure_Atom* r : _root_atoms)
  {
    r->collect_all_atoms(completed);
  }

  return;
}

void
Single_Substructure_Query::_determine_if_ring_ids_are_present()
{
  _ring_ids_present = 0;

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->ring_ids_present())
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

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->fused_system_ids_present())
    {
      _fused_system_ids_present = 1;
      return;
    }
  }

  return;
}

void 
Single_Substructure_Query::_determine_if_unmatched_atom_counts_are_present()
{
  _unmatched_atoms_attached_specified = 0;

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->unmatched_atoms_attached_specified()) {
      _unmatched_atoms_attached_specified++;
    }
  }

//cerr << "_determine_if_unmatched_atom_counts_are_present found " << _unmatched_atoms_attached_specified << '\n';

  return;
}

void
Single_Substructure_Query::_determine_if_atom_type_groups_are_present()
{
  _atom_type_groups_present = 0;

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->atom_type_groups_present())
      _atom_type_groups_present++;
  }

  return;
}

void
Single_Substructure_Query::_determine_if_spinach_specifications_are_present()
{
  _spinach_match_requirements_present = 0;

  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->spinach_match_specified())
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

//cerr << "Single_Substructure_Query::_determine_if_symmetry_groups_present " << _first_root_atom_with_symmetry_group << '\n';

  return;
}

int
Single_Substructure_Query::_compute_attribute_counts()
{
  int rc = 0;

  for (Substructure_Atom* r : _root_atoms)
  {
    rc += r->count_attributes_specified();
  }

  for (Substructure_Environment * e : _environment)
  {
    rc += e->count_attributes_specified();
  }

  for (Substructure_Environment* e : _environment_rejections)
  {
    rc += e->count_attributes_specified();
  }

  return rc;
}

int
Single_Substructure_Query::set_find_one_embedding_per_atom(int i)
{
  _find_one_embedding_per_start_atom = i;

  return 1;
}

int
Single_Substructure_Query::set_find_unique_embeddings_only(int s)
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
Single_Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches(int s)
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
Single_Substructure_Query::set_save_matched_atoms(int s)
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
Single_Substructure_Query::set_normalise_rc_per_hits_needed(int i)
{
  _normalise_rc_per_hits_needed = i;

  return 1;
}

int
Single_Substructure_Query::set_max_matches_to_find(int m)
{
  assert (m > 0);

  if (_hits_needed.is_set())
    cerr << "Single_Substructure_Query::set_max_matches_to_find:not set due to _hits_needed\n";
  else
    _max_matches_to_find = m;

  return 1;
}

int
Single_Substructure_Query::set_min_matches_to_find(int m)
{
  assert (m > 0);

  _hits_needed.set_min(m);

  return 1;
}

int
Single_Substructure_Query::set_min_atoms_to_match(int s)
{
  _natoms.set_min(s);

  return _natoms.ok();
}

int
Single_Substructure_Query::set_max_atoms_to_match(int s)
{
  _natoms.set_max(s);

  return _natoms.ok();
}

/*
  This looks wrong right now, seems as if it should do a recursive search...
*/

int
Single_Substructure_Query::involves_rings() const
{
  for (const Substructure_Atom* r : _root_atoms)
  {
    if (r->involves_rings())
      return 1;
  }

// Should put something in here if any of the global conditions specifies a ring

  return 0;
}

void
Single_Substructure_Query::set_ncon(int s)
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
Single_Substructure_Query::_has_implicit_rings(Query_Atoms_Matched & matched_query_atoms,
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
Single_Substructure_Query::_fused_system_id_conditions_satisfied(Query_Atoms_Matched & matched_query_atoms) const
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
    cerr << " ring " << i << " " << (*ri) << '\n';
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
    cerr << "Atom " << ai << " " << m->atomic_symbol(ai) << " fsid " << fsidi << " in a ring? " << m->is_ring_atom(ai) << '\n';
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
      cerr << "atoms " << ai << ' ' << m->atomic_symbol(ai) << " fsidi " << fsidi << " and " << aj << ' ' << m->atomic_symbol(aj) << " fsidj " << fsidj << " same ring system " << m->in_same_ring_system(ai, aj) << '\n';
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
Single_Substructure_Query::_ring_id_conditions_satisfied(Query_Atoms_Matched & matched_query_atoms) const
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
      cerr << "atoms " << ai << " ridi " << ridi << " and " << aj << " ridj " << ridj << " same ring " << m->in_same_ring(ai, aj) << '\n';
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
Single_Substructure_Query::_spinach_match_requirements_satisfied(Query_Atoms_Matched & matched_query_atoms,
                                        Molecule_to_Match & target_molecule) const
{
  for (const Substructure_Atom * sai : matched_query_atoms)
  {
    int s = sai->spinach_match_value();
    if (s < 0)    // not specified
      continue;

    Target_Atom * a = sai->current_hold_atom();

    if (nullptr == a)
      continue;

    atom_number_t matched_to = a->atom_number();

//#define DEBUG_SPINACH_MATCH_REQUIREMENTS_SATISFIED
#ifdef DEBUG_SPINACH_MATCH_REQUIREMENTS_SATISFIED
    cerr << "Query atom " << sai->unique_id() << " matched with atom " << matched_to << " spinach? " << target_molecule.is_spinach(matched_to) << '\n';
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
      cerr << " fragment match conditions i = " << i << " atom " << ai->current_hold_atom()->atom_number() << " mfi " << mfi << " fi " << fi << " j = " << j << " atom " << aj->current_hold_atom()->atom_number() << " mfj " << mfj << " fj " << fj << '\n';
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
    cerr << "Query atom " << q->unique_id() << " matched with atom " << a->atom_number() << " acon = " << acon << '\n';
#endif

    for (int j = 0; j < acon; j++)
    {
      const Target_Atom * aj = a->other(j).other();   // Target_Atom::other() returns a Bond_and_Target_Atom object

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
      cerr << "  atom " << aj->atom_number() << " z = " << aj->atomic_number() << " is attached";
      if( already_matched[aj->atom_number()])
        cerr << " already matched";
      cerr << '\n';
#endif

      if ( already_matched[aj->atom_number()])
        continue;

      atomic_number_t z = aj->atomic_number();

      if (_heteroatoms.empty())
      {
        if (6 != z)
          ahc++;
      }
      else if (_heteroatoms.contains(z))
        ahc++;
    }
  }

#ifdef DEBUG_ATTACHED_HETEROATOM_COUNT
  cerr << "Final attached heteroatom count is " << ahc << " match is " << _attached_heteroatom_count.matches(ahc) << '\n';
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

// #define DEBUG_MATCHED_ATOM_POST_CHECK

// An embedding has been found. Do the atoms in that embedding satisfy
// any conditions that may be present that impose conditions on an
// embedding.
int
Single_Substructure_Query::_global_query_conditions_also_matched(Query_Atoms_Matched & matched_query_atoms,
                            const int * already_matched,
                            Molecule_to_Match & target_molecule) const
{
#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Checking global conditions\n";
#endif

  if (_root_atoms.number_elements() > 1 && _distance_between_root_atoms.is_set()) {
    if (! _distance_between_root_atoms_satisfied(matched_query_atoms))
      return 0;
  }

  // cerr << "_no_matched_atoms_between " << _no_matched_atoms_between.size() << '\n';
  if (_no_matched_atoms_between.number_elements()) {
    if (! _no_matched_atoms_between_satisfied(matched_query_atoms))
      return 0;
  }

  // cerr << "_link_atom " << _link_atom.size() << '\n';
  if (_link_atom.number_elements()) {
    if (! _link_atoms_satisfied(matched_query_atoms))
      return 0;
  }
  
  if (_down_the_bond.number_elements()) {
    if (! _down_the_bond_satisfied(target_molecule, matched_query_atoms)) {
      return 0;
    }
  }

  if (_substituent.size()) {
    if (! _substituent_satisfied(target_molecule, matched_query_atoms)) {
      return 0;
    }
  }

  // cerr << "nmab " << _nmab.size() << " nmab items\n";
  if (_nmab.number_elements()) {
    if (! _nmab_satisfied(target_molecule, matched_query_atoms)) {
      return 0;
    }
  }

  if (! _attached_heteroatom_count_satisfied(matched_query_atoms, already_matched)) {
    return 0;
  }

  if (! _embedding_ncon_satisfied(matched_query_atoms, already_matched)) {
    return 0;
  }

//  cerr << "_fragment_ids_present " << _fragment_ids_present << '\n';
  if (_fragment_ids_present && ! _fragment_id_conditions_satisfied(matched_query_atoms)) {
    return 0;
  }

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Check ring ids? " << _ring_ids_present << '\n';
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

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "_unmatched_atoms_attached_specified value " << _unmatched_atoms_attached_specified << '\n';
#endif
  if (_unmatched_atoms_attached_specified)
  {
    if (! _unmatched_atoms_attached_matched(matched_query_atoms, already_matched, target_molecule))
      return 0;
  }

  if (matched_query_atoms.empty())   // don't know how to handle the unmatched atoms stuff
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

//  cerr << "Fraction atoms matched " << f << '\n';

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


  if (_geometric_constraints.number_elements() > 0) {
    Set_of_Atoms embedding;
    for (const auto * q : matched_query_atoms) {
      embedding << q->atom_number_matched();
    }
    for (const geometric_constraints::SetOfGeometricConstraints* constraint : _geometric_constraints) {
      if (!constraint->Matches(*target_molecule.molecule(), embedding)) {
        return 0;
      }
    }
  }

  if (_separated_atoms.number_elements() > 0) {
    std::unique_ptr<Set_of_Atoms> embedding = _make_new_embedding(matched_query_atoms);
    for (const SeparatedAtoms * separated_atoms : _separated_atoms) {
      if (! separated_atoms->Matches(*target_molecule.molecule(), *embedding)) {
        return 0;
      }
    }
  }

  if (_region.size() > 0) {
    for (const Region* r : _region) {
      if (! r->Matches(target_molecule, matched_query_atoms, already_matched)) {
        return 0;
      }
    }
  }

#ifdef DEBUG_MATCHED_ATOM_POST_CHECK
  cerr << "Global conditions OK\n";
#endif

  return 1;
}

/*
  This will be expensive.  We have a list of matched atoms and their associated
  symmetry group specification.  first task is to collect all that info from the
  atoms.  Note that this could be done once and stored in the
  Single_Substructure_Query object.  but since this is not used a lot, we do not
  worry about doing that.  fix if it ever becomes a problem...
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
    cerr << " i = " << i << " symmetry group " << si << ", atom " << mi << " symm " << msim[mi] << '\n';
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

// We have one or more embeddings in Query_Atoms_Matched.
int
Single_Substructure_Query::AllMatchesFractionAtomsMatched(const Substructure_Results & sresults,
                Molecule_to_Match& target_molecule) const {
  const int matoms = target_molecule.natoms();

  std::unique_ptr<int[]> hit_atom(new_int(matoms));
  constexpr int kOne = 1;
  sresults.each_embedding_set_vector(hit_atom.get(), kOne);
  int atoms_matched = std::count(hit_atom.get(), hit_atom.get() + matoms, kOne);

  const float f = iwmisc::Fraction<float>(atoms_matched, matoms);
  // cerr << "Matched " << atoms_matched << " of " << matoms << " f " << f << '\n';
  if (f < _min_all_matches_fraction_atoms_matched) {
    return 0;
  }
  if (f > _max_all_matches_fraction_atoms_matched) {
    return 0;
  }
  return 1;
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
    cerr << "Adding embedding matched atoms contains " << matched_atoms.number_elements() << " embedding " << new_embedding->number_elements() << '\n';
    cerr << (*new_embedding) << '\n';
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
Single_Substructure_Query::_find_next_root_atom_embedding(Query_Atoms_Matched & matched_atoms,
                                    Molecule_to_Match & target_molecule,
                                    int * already_matched,
                                    std::unique_ptr<int[]>& matched_by_global_specs,
                                    Substructure_Results & results)
{
  _iroot++;
  assert (_root_atoms.ok_index(_iroot));

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
  cerr << "_find_next_root_atom_embedding: iroot = " << _iroot << '\n';
#endif

  Substructure_Atom * root_atom = _root_atoms[_iroot];

  assert (nullptr == root_atom->current_hold_atom());

  int istart, istop;
  if (! root_atom->determine_start_stop(target_molecule, istart, istop))
  {
    _iroot--;
    return 0;
  }

  assert (istart >= 0 && istart <= istop && istop <= target_molecule.natoms());

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
  cerr << "Looking between atoms " << istart << " and " << istop << '\n';
#endif

  int rc = 0;
  for (int i = istart; i < istop; i++)
  {
#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Testing atom " << i << " already matched = " << already_matched[i] << '\n';
#endif

    if (already_matched[i])
      continue;

    Target_Atom & a = target_molecule[i];

    if (! root_atom->matches(a, already_matched))
      continue;

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Looking for embedding, root atom " << _iroot << " matched atom " << i << '\n';
#endif

    int tmp = _find_embedding(target_molecule, a, matched_atoms, already_matched, matched_by_global_specs, root_atom, results);

#ifdef DEBUG_FIND_NEXT_ROOT_ATOM_EMBEDDING
    cerr << "Found " << tmp << " matches at matched atom " << i << '\n';
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
  cerr << "_find_next_root_atom_embedding, root " << _iroot << " returning " << rc << '\n';
#endif

  _iroot--;

  return rc;
}


// Based on the matched atoms in `matched_atoms`, generate an embedding derived from
// those matched atoms.
// Must respect each atom's include_in_embedding attribute, as well as our own
// _respect_initial_atom_numbering.
// If an atom is not included, that position in `new_embedding` will be INVALID_ATOM_NUMBER.
std::unique_ptr<Set_of_Atoms>
Single_Substructure_Query::_make_new_embedding(const Query_Atoms_Matched& matched_atoms) const {
  std::unique_ptr<Set_of_Atoms> new_embedding(new Set_of_Atoms()); // Will be returned.
  int number_matched_atoms = matched_atoms.number_elements();

  if (_respect_initial_atom_numbering)
  {
    assert (_highest_initial_atom_number >= 0);
    new_embedding->extend(_highest_initial_atom_number + 1, INVALID_ATOM_NUMBER);
  }
  else
    new_embedding->resize(number_matched_atoms);

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Scanning " << number_matched_atoms << " matched atoms, _respect_initial_atom_numbering " << _respect_initial_atom_numbering << '\n';
#endif

  for (int i = 0; i < number_matched_atoms; i++)
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
//  cerr << "Processed atom " << ma << " initial " << a->initial_atom_number() << '\n';
#endif
  }

  return new_embedding;
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
                                    std::unique_ptr<int[]>& matched_by_global_specs,
                                    Substructure_Results & results)
{
#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Got embedding, " << matched_atoms.number_elements () << " atoms,  iroot = " << _iroot << " nroot = " << _root_atoms.number_elements() << '\n';
#endif

  if (_iroot < _root_atoms.number_elements() - 1)
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "iroot " << _iroot << ", " << matched_atoms.number_elements() << " atoms matched\n";
#endif

    return _find_next_root_atom_embedding(matched_atoms, target_molecule, already_matched, matched_by_global_specs, results);
  }

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Processing embedding with " << matched_atoms.number_elements() << " atoms matched: ";
  for (int i = 0; i < matched_atoms.number_elements(); i++)
  {
    cerr << ' ' << matched_atoms[i]->atom_number_matched();
  }
  cerr << '\n';
#endif

  if (! _global_query_conditions_also_matched(matched_atoms, already_matched, target_molecule))
    return 0;

  if (_chirality.number_elements() && ! _chiral_atoms_matched(matched_atoms, target_molecule))
    return 0;

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "_got_embedding: match OK before environment\n";
#endif

  if (! _matched_atoms_match_also_matched(matched_atoms, target_molecule)) {
    return 0;
  }

  _matches_before_checking_environment++;

  if (! _query_environment_also_matched(matched_atoms, target_molecule.natoms(), matched_by_global_specs)) {
    return 0;
  }

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "_got_embedding: environment also matched\n";
#endif

  if (! _save_matched_atoms) {
    results.got_embedding();

    if (0 == results.hits_found() % report_multiple_hits_threshold)
      cerr << "Query got " << results.hits_found() << " hits\n";

//  if (results.hits_found() > 40000000)
//    exit(0);

    if (_max_matches_to_find > 0 &&
        results.hits_found() >= static_cast<unsigned int>(_max_matches_to_find)) {
      results.set_complete();
      return _max_matches_to_find;
    }

    return 1;
  }

  std::unique_ptr<Set_of_Atoms> new_embedding = _make_new_embedding(matched_atoms);

// Not sure what it would mean if all atoms had been excluded from the embedding
// so for now, let's prohibit that

  assert (new_embedding->number_elements() > 0);

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Got embedding " << *new_embedding << ", find unique? " << _find_unique_embeddings_only << '\n';
#endif

  if (_find_unique_embeddings_only && ! results.embedding_is_unique(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding is not unique, deleting it\n";
#endif
    return 0;
  }

  if (_embeddings_do_not_overlap && results.embedding_overlaps_previous_embeddings(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding overlaps previous embeddings, deleting it\n";
#endif

    return 0;
  }

  if (_do_not_perceive_symmetry_equivalent_matches && 
       results.embedding_is_symmetry_related(*new_embedding))
  {
#ifdef DEBUG_GOT_EMBEDDING
    cerr << "Embedding symmetry related, discarded\n";
#endif
    return 0;
  }

  if (_atom_type_groups_present && ! _atom_type_groupings_matched(matched_atoms)) {
    return 0;
  }

#ifdef DEBUG_GOT_EMBEDDING
  cerr << "Finally, got an embedding\n";
#endif

  return _add_embedding(matched_atoms, new_embedding.release(), results);
}

int
remove_atoms_with_same_or(Query_Atoms_Matched & matched_atoms,
                          const int istart,
                          const int orid)
{
  int rc = 0;
  int number_matched_atoms = matched_atoms.number_elements();
  for (int i = istart; i < number_matched_atoms; i++)
  {
    const Substructure_Atom * b = matched_atoms[i];
    if (b->or_id() == orid)
    {
      matched_atoms.remove_item(i);
      i--;
      number_matched_atoms--;
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
                                     std::unique_ptr<int[]>& matched_by_global_specs,
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
    return _got_embedding(matched_atoms, target_molecule, already_matched, matched_by_global_specs, results);
  }

  int rc = 0;

  while (atom_to_process > number_initially_matched)    // we ignore any which may have previously been matched
  {
    Substructure_Atom * a = (Substructure_Atom *) matched_atoms[atom_to_process];

#ifdef DEBUG_FIND_EMBEDDING
    cerr << "Atom to process = " << atom_to_process << " id " << a->unique_id() << 
            " array contains " << matched_atoms.number_elements() << '\n';
    for (int i = 0; i < matched_atoms.number_elements(); i++)
    {
      const Substructure_Atom * a = matched_atoms[i];
      cerr << " " << a->unique_id();
      if (a->is_matched())
        cerr << "(" << a->atom_number_matched() << ")";
    }
    cerr << '\n';
    if (nullptr == a->parent())
      cerr << "Returning " << rc << '\n';
#endif

    if (nullptr == a->parent())    // must be a root atom, done
      return rc;

    if (! a->move_to_next_match_from_current_anchor(already_matched, matched_atoms))
    {
#ifdef DEBUG_FIND_EMBEDDING
      cerr << "Move to next failed for atom " << a->unique_id() << '\n';
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
              " atom to process = " << atom_to_process << " matched = " << matched_atoms.number_elements() << '\n';
#endif
      if (a->or_id())
        remove_atoms_with_same_or(matched_atoms, atom_to_process + 1, a->or_id());

      a->add_your_children(matched_atoms);   // does nothing if already added

      if (atom_to_process == matched_atoms.number_elements() - 1)
      {
#ifdef DEBUG_FIND_EMBEDDING
        cerr << "AlL query atoms matched, calling _got_embedding, B4 have " << results.number_embeddings() << " rc " << rc << '\n';
#endif
        if (_got_embedding(matched_atoms, target_molecule, already_matched, matched_by_global_specs, results)) {
          rc++;
  
#ifdef DEBUG_FIND_EMBEDDING
          cerr << "Rc incremented to " << rc << " have " << results.number_embeddings() << " embeddings\n";
          if (results.matching_complete() || _find_one_embedding_per_start_atom)
            cerr << "Returning " << rc << '\n';
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
  cerr << "Single_Substructure_Query::_find_embedding: returning rc = " << rc << '\n';
#endif

  return rc;
}

int
Single_Substructure_Query::find_embedding(Molecule_to_Match & target_molecule,
                                     Target_Atom & a,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results) {
  if (_need_to_compute_aromaticity < 0) {    // first time this query has been invoked
    _one_time_initialisations();
  }

  // Not set here, but needed to satisfy the args.
  std::unique_ptr<int[]> matched_by_global_specs;

  int save_iroot = _iroot;
  _iroot = 0;
  int rc = _find_embedding(target_molecule, a, matched_atoms, already_matched,
                        matched_by_global_specs, root_atom, results);
  _iroot = save_iroot;
  return rc;
}

int
Single_Substructure_Query::_find_embedding(Molecule_to_Match & target_molecule,
                                     Target_Atom & a,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     std::unique_ptr<int[]>& matched_by_global_specs,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results)
{
  int initial_size = matched_atoms.number_elements();

  root_atom->set_hold(&a, already_matched);

#ifdef DEBUG_FIND_EMBEDDING
  cerr << "Single_Substructure_Query::_find_embedding:atom " << a.atom_number() << " begin have " << results.number_embeddings() << " embdddings\n";
#endif

  int rc = _find_embedding(target_molecule, matched_atoms, already_matched, matched_by_global_specs, root_atom, results);

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
Single_Substructure_Query::_substructure_search(Molecule_to_Match & target_molecule,
                                                int * already_matched,
                                                std::unique_ptr<int[]>& matched_by_global_specs,
                                                Substructure_Results & results)
{
//int matoms = target_molecule.natoms();

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "Beginning atom matching over " << target_molecule.natoms() << " atoms\n";
  cerr << "_respect_initial_atom_numbering " << _respect_initial_atom_numbering << '\n';
#endif

  _iroot = 0;

  Substructure_Atom * root_atom = _root_atoms[_iroot];

  int jstart, jstop;
  if (INVALID_ATOM_NUMBER != target_molecule.start_matching_at())
  {
    jstart = target_molecule.start_matching_at();
    jstop = jstart + 1;
  }
  else if (! root_atom->determine_start_stop(target_molecule, jstart, jstop))
    return 0;

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "Start atoms " << jstart << " and " << jstop << '\n';
#endif

//assert (jstart >= 0 && jstop >= jstart && jstop <= matoms);

  Query_Atoms_Matched matched_atoms;

  for (int j = jstart; j < jstop; j++)   // loop over atoms in molecule
  {
    if (already_matched[j])   // will only be the case when embeddings are not allowed to overlap
      continue;

    Target_Atom & target_atom = target_molecule[j];

    if (! root_atom->matches(target_atom, already_matched)) {
      continue;
    }

    _max_query_atoms_matched.extra(1);
    results.matched_this_many_query_atoms(1);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
    cerr << "Atom " << j << " matches query atom 0, " << results.number_embeddings() << " embeddings so far, root " << _iroot << '\n';
#endif

    const int tmp = _find_embedding(target_molecule, target_atom, matched_atoms, already_matched, matched_by_global_specs, root_atom, results);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
    if (tmp)
    {
      cerr << "Found embedding starting with atom " << j << ", tmp = " << tmp;
      if (results.number_embeddings())
        cerr << " last embedding " << *results.embedding(results.number_embeddings() - 1);
      cerr << '\n';
    }
    else
      cerr << "No embedding found, tmp = " << tmp << '\n';
    for (int i = 0; i < target_molecule.natoms(); ++i)
    {
      cerr << " alread matched " << i << ' ' << already_matched[i] << target_molecule.molecule()->smarts_equivalent_for_atom(i) << '\n';
    }
#endif

    if (tmp && _all_hits_in_same_fragment)
    {
      int k = target_atom.fragment_membership();
#ifdef DEBUG_SUBSTRUCTURE_QUERY
      cerr << "Hit with atom " << target_atom.atom_number() << " in frag " << k << '\n';
#endif
      results.got_hit_in_fragment(k);
    }

    if (results.matching_complete())
      return results.return_code();

//  NOT correct
//  if (tmp && _find_one_embedding_per_start_atom)
//    return results.return_code();

    if (0 == _embeddings_do_not_overlap)   // that is, embeddings can overlap..
      set_vector(already_matched, target_molecule.natoms(), 0);
  }

  assert (0 == _iroot);

#ifdef DEBUG_SUBSTRUCTURE_QUERY
  cerr << "_substructure_search: returning " << results.number_embeddings() << '\n';
#endif

  return results.return_code();
}

int
Single_Substructure_Query::_match_elements_needed(Molecule_to_Match& target_molecule) const {
  for (const Elements_Needed* e : _elements_needed) {
    if (! e->matches(target_molecule)) {
      return 0;
    }
  }

  return 1;
}

int
RequiredMolecularProperties::MatchElementsNeeded(Molecule_to_Match & target_molecule) const {
  for (const Elements_Needed* e : _elements_needed) {
    if (! e->matches(target_molecule)) {
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_match_heteroatom_specifications(Molecule_to_Match & target_molecule)
{
  if (_heteroatoms.empty())
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

  return _heteroatoms_in_molecule.matches(heteroatoms_in_target);
}

/*
  Note that each of the ring_specifications could match the same ring.
  Would not be too hard to change if ever that becomes an issue
*/

int
Single_Substructure_Query::_match_ring_system_specifications(Molecule_to_Match & target_molecule,
                        std::unique_ptr<int[]>& matched_by_global_specs)
{
#ifdef DEBUG_MATCH_RING_SYSTEM_SPECIFICATIONS
  cerr << "Single_Substructure_Query::_match_ring_system_specifications:checking " << _ring_system_specification.number_elements() << " ring system specifications\n";
#endif

  int nr = _ring_system_specification.number_elements();

  if (1 == nr)
    return _ring_system_specification[0]->matches(target_molecule, matched_by_global_specs);

  _ring_system_specification_logexp.reset();

  for (int i = 0; i < nr; i++)
  {
    if (! _ring_system_specification_logexp.result_needed(i))
      continue;

    int m = _ring_system_specification[i]->matches(target_molecule, matched_by_global_specs);

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

//#define DEBUG_MATCH_RING_SPECIFICATIONS

int
Single_Substructure_Query::_match_ring_specifications(Molecule_to_Match & target_molecule,
                std::unique_ptr<int[]> & matched_by_global_specs)
{
  int nr = _ring_specification.number_elements();

#ifdef DEBUG_MATCH_RING_SPECIFICATIONS
  cerr << "Single_Substructure_Query::_match_ring_specifications:tsting " << nr << " ring specifications\n";
#endif

  if (1 == nr) {
    return _ring_specification[0]->matches(target_molecule, matched_by_global_specs);
  }

  _ring_specification_logexp.reset();

  for (int i = 0; i < nr; i++)
  {
    if (! _ring_specification_logexp.result_needed(i))
      continue;

    Substructure_Ring_Specification * ri = _ring_specification[i];

    int m = ri->matches(target_molecule, matched_by_global_specs);

    _ring_specification_logexp.set_result(i, m > 0);

    int rc;

    if (_ring_specification_logexp.evaluate(rc)) {
      return rc;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_match_nrings_specifications(Molecule_to_Match & target_molecule)
{
  return _nrings.matches(target_molecule.nrings());
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
Single_Substructure_Query::_spinach_atoms_match(Molecule_to_Match & target) const
{
  Molecule & m = *(target.molecule());

  int matoms = m.natoms();

  int * spinach = new int[matoms]; std::unique_ptr<int[]> free_spinach(spinach);

  int ais = m.identify_spinach(spinach);

//cerr << "Single_Substructure_Query::_spinach_atoms_match:ais " << ais << '\n';
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

//cerr << "Single_Substructure_Query::_spinach_atoms_match:number_inter_ring_atoms " << number_inter_ring_atoms << '\n';
  return _inter_ring_atoms.matches(number_inter_ring_atoms);
}

/*
  If is very common for there to be no global conditions present, so we can often save time by detecting that once
  and then not checking on subsequent calls. Make sure you update this function as new global conditions are added!
*/

int
Single_Substructure_Query::_discern_global_conditions_present ()
{
  if (_required_molecular_properties.attributes_specified()) {
    return 1;
  }

  if (_natoms.is_set() || _heteroatoms_in_molecule.is_set() || _nrings.is_set() || _number_isotopic_atoms.is_set() ||
      _number_fragments.is_set() || _ring_specification.number_elements() || _ring_system_specification.number_elements() || 
      _elements_needed.number_elements() || _atoms_in_spinach.is_set() || _inter_ring_atoms.is_set() || 
      _aromatic_rings.number_elements() || _non_aromatic_rings.is_set() || _isolated_rings.is_set() || _isolated_ring_objects.is_set() || 
      _fused_rings.is_set() || _strongly_fused_rings.is_set() || _net_formal_charge.is_set() ||
      _first_root_atom_with_symmetry_group >= 0 || _aromatic_atoms.is_set())
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

  `matched_by_global_specs` might be set by one of the functions called. Those
  entities might place a per-atom value in that array, which will later be
  used for checking Global ID values.

  Make sure that _nrings is checked first!
*/

int
Single_Substructure_Query::_match_global_specifications(Molecule_to_Match & target_molecule,
                std::unique_ptr<int[]>& matched_by_global_specs)
{
#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Checking global specifications: _natoms " << _natoms.is_set() << ", matches? " << _natoms.matches(target_molecule.natoms()) << '\n';
  cerr << "Target contains " << target_molecule.natoms() << " atoms\n";
#endif

  // This is a kludge, need to fix this to get clear separation and efficiency.
  // TODO:ianwatson clean this up
  if (_global_conditions_present < 0) {
    _discern_global_conditions_present();
    if (_global_conditions_present == 0 && _required_molecular_properties.attributes_specified()) {
      _global_conditions_present =1;
    }
  }

  if (0 == _global_conditions_present)
    return 1;

  // First check those properties that update global specs
  if (_ring_specification.number_elements()) {
    if (! _match_ring_specifications(target_molecule, matched_by_global_specs)) {
      return 0;
    }
  }

  if (_ring_system_specification.number_elements()) {
    if (! _match_ring_system_specifications(target_molecule, matched_by_global_specs)) {
      return 0;
    }
  }

  if (_inter_ring_region.size() > 0) {
    if (! InterRingRegionsMatch(target_molecule, matched_by_global_specs)) {
      return 0;
    }
  }

  if (_required_molecular_properties.attributes_specified()) {
    return _required_molecular_properties.Matches(target_molecule);
  }

  if (_natoms.is_set())
  {
    if (! _natoms.matches(target_molecule.natoms()))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check heteroatoms " << _heteroatoms_in_molecule.is_set() << '\n';
#endif

  if (_heteroatoms_in_molecule.is_set()) {
    if (! _match_heteroatom_specifications(target_molecule)) {
      return 0;
    }
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check nrings " << _nrings.is_set() << '\n';
#endif

  if (_nrings.is_set())
  {
    if (! _nrings.matches(target_molecule.nrings()))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check isotopes " << _number_isotopic_atoms.is_set() << '\n';
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
  cerr << "Do we need to check ring specification(s) " << _ring_specification.number_elements() << '\n';
#endif

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check elements needed? " << _elements_needed.number_elements() << '\n';
#endif

  if (_elements_needed.number_elements())
  {
    if (! _match_elements_needed(target_molecule))
      return 0;
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check ring type specifications\n";
#endif

  if (! _match_ring_type_specifications(target_molecule)) {
    return 0;
  }

  if (_atoms_in_spinach.is_set() || _inter_ring_atoms.is_set()) {
    if (! _spinach_atoms_match(target_molecule)) {
      return 0;
    }
  }

  if (_net_formal_charge.is_set() && ! _net_formal_charge.matches(target_molecule.net_formal_charge()))
    return 0;

  if (_aromatic_atoms.is_set() && ! _aromatic_atoms_matches(target_molecule)) {
    return 0;
  }

  return 1;
}

int
Single_Substructure_Query::_aromatic_atoms_matches(Molecule_to_Match& target_molecule) const
{
  const int matoms = target_molecule.natoms();
  int aromatic_atom_count = 0;
  for (int i = 0; i < matoms; ++i) {
    if (target_molecule[i].is_aromatic())
      aromatic_atom_count++;
  }

  return _aromatic_atoms.matches(aromatic_atom_count);
}


//#define DEBUG_SUBSTRUCTURE_SEARCH

/*
  All queries come through this point.
  We release all holds on the way out. This is because if we don't,
  some of the atoms may still be pointing to target atoms which
  have been deleted.
*/

int
Single_Substructure_Query::_substructure_search(Molecule_to_Match & target_molecule,
                                                Substructure_Results & results)
{
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Begin common _substructure_search code, " << _root_atoms.number_elements() << " root atoms\n";
  cerr << "_respect_initial_atom_numbering " << _respect_initial_atom_numbering << '\n';
#endif

  // Perhaps one of the global conditions will mark some atoms.
  std::unique_ptr<int[]> matched_by_global_specs;

  if (! _match_global_specifications(target_molecule, matched_by_global_specs))
  {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
    cerr << "Global specifications not matched\n";
#endif

    return 0;
  }
  if (_root_atoms.empty()) {
    return 1;
  }

  if (_do_not_perceive_symmetry_equivalent_matches)
    results.set_symmetry_class(target_molecule.molecule()->symmetry_classes());
  else
    results.set_symmetry_class(nullptr);

  if (_atom_typing != nullptr) {
    if (! target_molecule.AssignAtomTypes(*_atom_typing))
      return 0;
  }

  int * tmp = new_int(target_molecule.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  int rc = _substructure_search(target_molecule, tmp, matched_by_global_specs, results);

  for (auto* root_atom : _root_atoms) {
    root_atom->recursive_release_hold();
  }

  if (matched_by_global_specs && rc > 0) {
    rc = results.RemoveEmbeddingsNotSatisfyingGlobalId(matched_by_global_specs.get());
  }

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Return code from _substructure_search " << rc << " results contains " << results.number_embeddings() << " embeddings\n";
  if (results.number_embeddings())
    cerr << "     last embedding " << *results.embedding(results.number_embeddings() - 1) << '\n';
#endif

  if (_rejection) {
    return ! rc;
  }

  return rc;
}

// The first time a query is invoked, some initialisations are needed.

int
Single_Substructure_Query::_one_time_initialisations()
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

  _determine_if_unmatched_atom_counts_are_present();

  _determine_if_atom_type_groups_are_present();

  if (_respect_initial_atom_numbering)
  {
    for (int i = 0; i < _root_atoms.number_elements(); i++)
    {
      int tmp = _root_atoms[i]->highest_initial_atom_number();
      if (tmp > _highest_initial_atom_number)
        _highest_initial_atom_number = tmp;
    }
  }

  return 1;
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
Single_Substructure_Query::substructure_search(Molecule_to_Match & target_molecule,
                                               Substructure_Results & results)
{
  assert (target_molecule.ok());

  const int matoms = target_molecule.natoms();

  results.initialise(matoms);     // no returns before this.

  if (nullptr == _fingerprint)
    ;
  else if (! target_molecule.is_superset(*_fingerprint))
    return 0;

//cerr << " cmp " << matoms << " and _min_atoms_in_query " << _min_atoms_in_query << '\n';
     
  if (matoms < _min_atoms_in_query)
    return 0;     

  if (_need_to_compute_aromaticity < 0)    // first time this query has been invoked
  {
    _one_time_initialisations();
  }

// If 2 == aromatic_bonds_lose_kekule_identity(), we need a temporary array
// of bond types. Even if the query doesn't specify aromatic bonds, we
// need to convert aromatic rings so they match single bonds

  bond_type_t * save_bt = nullptr;

  if (2 == aromatic_bonds_lose_kekule_identity())
  {
    Molecule * m = target_molecule.molecule();

    save_bt = new bond_type_t[m->nedges()];

    m->get_bond_types(save_bt);

    if (0 == m->set_bond_types_for_isis_aromaticity_matching())
    {
      delete [] save_bt;
      save_bt = nullptr;
    }
  }
  else if (_need_to_compute_aromaticity)
    target_molecule.establish_aromatic_bonds();
  else if (_need_to_compute_ring_membership)
    (void) target_molecule.molecule()->ring_membership();

  results.set_save_matched_atoms(_save_matched_atoms);

  int nf = 1;
  if (_all_hits_in_same_fragment) {
    nf = target_molecule.molecule()->number_fragments();
    results.size_hits_per_fragment_array(nf);
  } else if (_only_keep_matches_in_largest_fragment) {
    nf = target_molecule.molecule()->number_fragments();
  }

  int rc = _substructure_search(target_molecule, results);

//target_molecule.debug_print(cerr);

  if (nullptr != save_bt) {
    target_molecule.molecule()->set_bond_types_no_set_modified(save_bt);
    delete [] save_bt;
  }

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Return code from ss = " << rc;
  if (_hits_needed.is_set())
    cerr << " Hits needed is set. match " << _hits_needed.matches(rc) << '\n';
  else
    cerr << " Hits needed not set, result " << rc << '\n';
#endif

  if (0 == rc) {
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
  cerr << "rc = " << rc << " _distance_between_hits " << _distance_between_hits.is_set() << '\n';
#endif

  if (rc > 1 && _distance_between_hits.is_set()) {
    int tmp = results.remove_hits_violating_distance(target_molecule, _distance_between_hits, _matched_atoms_to_check_for_hits_too_close);   // the number of hits removed

    if (tmp && _fail_if_embeddings_too_close)
      return 0;

    rc = rc - tmp;
    assert (rc > 0);   // it cannot have removed all the hits
  }

// Do any adjustments for fragment stuff

  if (0 == rc) {  // no embeddings to worry about
    ;
  } else if (_all_hits_in_same_fragment) {
    if (nf > 1 && rc > 1) {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
      cerr << "nf = " << nf << " and all hits in same frag\n";
#endif

      const extending_resizable_array<int> & hits_per_fragment = results.hits_per_fragment();

      int found_match = 0;
      for (int i = 0; i < nf; i++) {
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
        cerr << hits_per_fragment[i] << " hits in fragment " << i << " matches " << _hits_needed.matches(hits_per_fragment[i]) << '\n';
#endif
        if (_hits_needed.matches(hits_per_fragment[i])) {
          found_match = 1;
          rc = hits_per_fragment[i];
          break;
        }
      }
  
      if (! found_match)
        return 0;
    }
  } else if (_only_keep_matches_in_largest_fragment) {
    if (rc > 0 && nf > 1) {
      results.remove_hits_not_in_largest_fragment(target_molecule);
      rc = results.number_embeddings();
    }
  }

  if (! _hits_needed.matches(rc))
    return 0;

  if (_min_all_matches_fraction_atoms_matched > 0.0f ||
      _max_all_matches_fraction_atoms_matched > 0.0f) {
    if (!AllMatchesFractionAtomsMatched(results, target_molecule)) {
      return 0;
    }
  }

  if (_compress_embeddings && rc > 1)
    rc = results.CompressToSingleEmbedding();

// check on any modifications to the return code.

// This is somewhat arbitrary. By definition, we honour _subtract_from_rc if
// it is specified, and in that case, ignore all other modifications.

  if (_subtract_from_rc) {
    rc -= _subtract_from_rc;
    if (rc < 0)
      return 0;

    return rc;  // Arbitrary decision
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
Single_Substructure_Query::_examine_bond_specifications()
{
  int nr = _root_atoms.number_elements();

  _need_to_compute_ring_membership = 0;

  for (int i = 0; i < nr; i++)
  {
    Substructure_Atom * root_atom = _root_atoms[i];
    if (root_atom->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
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
    Substructure_Environment * rej = _environment_rejections[i];
    if (rej->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
    {
      _need_to_compute_aromaticity = 1;
      break;
    }
  }

  if (_need_to_compute_aromaticity < 0)
    _need_to_compute_aromaticity = 0;

//cerr << "At end of _examine_bond_specifications, ring = " << _need_to_compute_ring_membership <<
//        " arom = " << _need_to_compute_aromaticity << '\n';

  return 1;
}

int
Single_Substructure_Query::substructure_search(Molecule * m,
                                               Substructure_Results & results)
{
//assert (ok ());
  assert (OK_MOLECULE(m));

  const int matoms = m->natoms();

// If the atom has fewer atoms than the "atoms" in the query, this cannot work

  if (matoms < _min_atoms_in_query)
    return 0;

// If the query contains ring(s), but the molecule contains fewer, there can be no match

  if (_rings_in_query > 0 && m->nrings() < _rings_in_query)
    return 0;

  Molecule_to_Match target(m);

  return substructure_search(target, results);
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
      os << "  environment component " << i << '\n';
      _environment[i]->print_environment_matches(os);
    }
  }

  if (_match_to_environemt_rejection)
  {
    for (int i = 0; i < nr; i++)
    {
      os << "  environment rejection " << i << '\n';
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

    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

Substructure_Atom *
Single_Substructure_Query::query_atom_with_atom_map_number(atom_number_t a) const
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    Substructure_Atom * rc = _root_atoms[i]->query_atom_with_atom_map_number(a);

    if (nullptr != rc)
      return rc;
  }

  return nullptr;
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

    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

const Substructure_Bond *
Single_Substructure_Query::bond_between_atom_map_numbers(int a1, int a2) const
{
  const Substructure_Bond * rc = nullptr;

  for (int i = 0; i < _root_atoms.number_elements(); ++i)
  {
    rc = _root_atoms[i]->bond_between_atom_map_numbers(a1, a2);

    if (nullptr != rc)
      return rc;
  }

  return nullptr;
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
    output << "ROOT " << i << '\n';
    _root_atoms[i]->print_connectivity_graph(output);
  }

  return 1;
}

int
Single_Substructure_Query::_unmatched_atoms_attached_matched(const Query_Atoms_Matched& matched_query_atoms,
       const int * already_matched,
       Molecule_to_Match& target) const
{
  for (Substructure_Atom* a : matched_query_atoms)
  {
    if (!a->unmatched_atoms_attached_matches(matched_query_atoms, already_matched, target))
      return 0;
  }

  return 1;
}
int
Molecule_to_Match::AssignAtomTypes(Atom_Typing_Specification& atom_typing) 
{
  if (_m == nullptr)
  {
    cerr << "Molecule_to_Match::AssignAtomTypes:no molecule\n";
    return 0;
  }

  atom_type_t * atype = new atom_type_t[_natoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  if (! atom_typing.assign_atom_types(*_m, atype))
  {
    cerr << "Molecule_to_Match::AssignAtomTypes:cannot assign atom types\n";
    return 0;
  }

  for (int i = 0; i < _natoms; ++i)
  {
    _target_atom[i].set_atom_type(atype[i]);
  }

  return 1;
}

//#define DEBUG_ATOM_TYPE_GROUPINGS_MATCHED

int
Single_Substructure_Query::_atom_type_groupings_matched(const Query_Atoms_Matched& matched_query_atoms) const
{
  const int n = matched_query_atoms.number_elements();
#ifdef DEBUG_ATOM_TYPE_GROUPINGS_MATCHED
  cerr << "Single_Substructure_Query::_atom_type_groupings_matched:checking " << n << " matched atoms\n";
#endif
  for (int i = 0; i < n; ++i)
  {
    const Substructure_Atom* ai = matched_query_atoms[i];

    if (ai->atom_type_group() == 0)
      continue;

    const auto atype = ai->current_hold_atom()->atom_type();
#ifdef DEBUG_ATOM_TYPE_GROUPINGS_MATCHED
    cerr << " matched atom " << i << " has type " << atype << " group " << ai->atom_type_group() << '\n';
#endif

    for (int j = i + 1; j < n; ++j)
    {
      const Substructure_Atom* aj = matched_query_atoms[j];
      if (aj->atom_type_group() == 0)
        continue;

#ifdef DEBUG_ATOM_TYPE_GROUPINGS_MATCHED
      cerr << "   matched atom " << j << " has type " << aj->current_hold_atom()->atom_type() << " group " << aj->atom_type_group() << '\n';
#endif

      if (ai->atom_type_group() == aj->atom_type_group())  // Same group, types must match.
      {
        if (atype != aj->current_hold_atom()->atom_type())
          return 0;
      }
      else   // Different groups, types must not match.
      {
        if (atype == aj->current_hold_atom()->atom_type())
          return 0;
      }
    }
  }

  return 1;
}

int
SeparatedAtoms::Matches(Molecule& m,
        const Set_of_Atoms& embedding) const {
  atom_number_t a1 = embedding[_a1];
  atom_number_t a2 = embedding[_a2];
  // cerr << "Matched atoms are " << a1 << " and " << a2 << " betw " << m.bonds_between(a1, a2) << '\n';
  return _separation.matches(m.bonds_between(a1, a2));
}

int
Single_Substructure_Query::ForgetOriginatingSmarts() {
  if (_originating_smarts.empty()) {
    return 0;
  }
  _originating_smarts.resize(0);
  return 1;
}

void
InterRingRegionData::BeginNewRegion(int number) {
  region_number = number;
  atoms_in_region = 0;
  ring_connections.resize_keep_storage(0);
}

// Identify the atoms in an inter-ring region
void
IdentifyRegion(InterRingRegionData& irrd,
               atom_number_t zatom) {
  irrd.region[zatom] = irrd.region_number;
  const Atom& a = irrd.m->atom(zatom);
  ++irrd.atoms_in_region;
  for (const Bond * b : a) {
    atom_number_t j = b->other(zatom);
    if (irrd.region[j] == irrd.region_number) {
      continue;
    }
    if (irrd.m->ring_bond_count(j)) {
      irrd.ring_connections << j;
      continue;
    }
    IdentifyRegion(irrd, j);
  }
}

int
Single_Substructure_Query::InterRingRegionsMatch(Molecule_to_Match& target,
                std::unique_ptr<int[]>& matched_by_global_specs) {
  // Inter ring regions need a minimum of 2 rings.
  if (target.nrings() < 2) {
    return 0;
  }
  // Look for something like naphthalene, no inter-ring regions there.
  if (target.nrings() == 2 && target.molecule()->ringi(0)->is_fused()) {
    return 0;
  }

  const int matoms = target.natoms();
  std::unique_ptr<int[]> region(new int[matoms]);
  target.molecule()->identify_spinach(region.get());
  // Invert so we identify scaffold atoms. Scaffold==1, spinach==0
  for (int i = 0; i < matoms; ++i) {
    region[i] = ! region[i];
  }

  const int number_region_specifications = _inter_ring_region.number_elements();

  // for each _inter_ring_region we need to keep track of which regions
  // it matches.
  std::vector<std::vector<int>> region_matched(number_region_specifications);

  // Ring atoms are numbered 1, so start region numbering at 2, 
  int region_number = 2;
  InterRingRegionData irrd(target.molecule(), region.get());

  for (int i = 0; i < matoms; ++i) {
    if (region[i] == 0) {  // spinach atom, not interested.
      continue;
    }
    if (target[i].ring_bond_count() > 0) {  // scaffold atom in a ring.
      continue;
    }
    // We have an atom that is NOT spinach and is NOT in a ring.
    // Therefore part of an inter-ring region.

    irrd.BeginNewRegion(region_number);
    IdentifyRegion(irrd, i);
#ifdef DEBUG_INTER_RING_REGION_MATCH
    cerr << " begin inter ring region " << region_number << " with atom " << i << '\n';
#endif
    for (int j = 0; j < number_region_specifications; ++j) {
      if (! _inter_ring_region[j]->Matches(target, region_number, irrd, matched_by_global_specs)) {
#ifdef DEBUG_INTER_RING_REGION_MATCH
        cerr << "IRR " << j << " does not match region " << region_number << '\n';
#endif
        continue;
      }
      region_matched[j].push_back(region_number);
    }
    ++region_number;
  }

  // Now check to see that each inter region specification matched the required
  // number of times.
  // Also include the region numbers of matched regions.
  resizable_array<int> regions_matched;
  int found_set_global_id = 0;
  for (int i = 0; i < number_region_specifications; ++i) {
#ifdef DEBUG_INTER_RING_REGION_MATCH
    cerr << "IRR " << i << " matches " << region_matched[i].size() << " regions\n";
#endif
    if (! _inter_ring_region[i]->MatchesNhits(region_matched[i].size())) {
      return 0;
    }
    if (_inter_ring_region[i]->global_id() > 0) {
      found_set_global_id = 1;
    }
    // Update global list of regions matched.
    for (int region : region_matched[i]) {
      regions_matched.add_if_not_already_present(region);
    }
  }

#ifdef DEBUG_INTER_RING_REGION_MATCH
  cerr << "All inter ring regions match, found_set_global_id " << found_set_global_id << '\n';
#endif
  if (! found_set_global_id) {
    return 1;
  }

  // For each region matched, mark the atoms in matched_by_global_specs and
  // set the global id in each of the atoms in `target`.
  // Note that if the global id is already set, we do not update it, so only
  // the first inter-ring region gets to update the global id. Need to change
  // to a bit field.
  if (! matched_by_global_specs) {
    matched_by_global_specs.reset(new_int(matoms));
  }

  for (int i = 0; i < number_region_specifications; ++i) {
    if (region_matched[i].empty()) {
      continue;
    }
    const int gid = _inter_ring_region[i]->global_id();

    if (gid <= 0) {
      continue;
    }
    for (int r : region_matched[i]) {
      for (int j = 0; j < matoms; ++j) {
        if (region[j] != r) {
          continue;
        }
        if (matched_by_global_specs[j] > 0) {
          continue;
        }
        target[j].set_global_id(gid);
        matched_by_global_specs[j] = gid;
      }
    } 
  }
 
  return 1;
}

void
SetRingAtoms(Query_Atoms_Matched& matched_query_atoms,
             int * ring_atom) {
  for (const Substructure_Atom* a : matched_query_atoms) {
    atom_number_t j = a->current_hold_atom()->atom_number();
    ring_atom[j] = 1;
  }
}

int
Single_Substructure_Query::_substituent_satisfied(Molecule_to_Match& target_molecule,
                Query_Atoms_Matched& matched_query_atoms) const {
  const int matoms = target_molecule.natoms();

  std::unique_ptr<int[]> matched_by_global_specs(new_int(matoms));
  std::unique_ptr<int[]> ring_atoms(new_int(matoms));
  SetRingAtoms(matched_query_atoms, ring_atoms.get());
  std::unique_ptr<int[]> storage(new_int(matoms));

  int rc = 0;
  for (Substituent* subst : _substituent) {
    if (! subst->Matches(target_molecule, ring_atoms.get(), storage.get(),
                matched_by_global_specs)) {
      continue;
    }
    ++rc;
  }

  return rc;
}

MatchedAtomMatch::MatchedAtomMatch() {
  _unique_id = -1;
}

//#define DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED

int
Single_Substructure_Query::_matched_atoms_match_also_matched(Query_Atoms_Matched & matched_query_atoms, 
                Molecule_to_Match& target) {
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
  cerr << "Single_Substructure_Query::_matched_atoms_match_also_matched:checking " << _matched_atom_match.size() << " MatchedAtomMatch\n";
  if (_matched_atom_match.empty()) {
    return 1;
  }
#endif

  std::unique_ptr<int[]> already_matched(new_int(target.molecule()->natoms()));

  for (MatchedAtomMatch* mam : _matched_atom_match) {
    if (! mam->Matches(matched_query_atoms, target, already_matched.get())) {
      return 0;
    }
  }

  return 1;
}

int
MatchedAtomMatch::Matches(Query_Atoms_Matched& matched_query_atoms,
                        Molecule_to_Match& target,
                        int * already_matched) {
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
  cerr << "MatchedAtomMatch::Matches: " << _atoms.size() << " atoms " << _positive_matches.size() << " positive amd " << _negative_matches.size() << " negative matches\n";
#endif
  const int matoms = target.natoms();

  Query_Atoms_Matched matched_atoms;
  Substructure_Results sresults;

  for (int atom : _atoms) {
    if (! matched_query_atoms.ok_index(atom)) {
      cerr << "MatchedAtomMatch::Matches:invalid atom number " << atom << " embedding has " << matched_query_atoms.size() << " items\n";
      return 0;
    }
    const Substructure_Atom* matched_atom = matched_query_atoms[atom];
    Target_Atom& target_atom = target[matched_atom->atom_number_matched()];
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
    cerr << "Target atom tpe " << target_atom.atomic_number() << '\n';
#endif

    if (_positive_matches.size() > 0) {
      bool got_match = false;
      for (Single_Substructure_Query* q : _positive_matches) {
        std::fill_n(already_matched, matoms, 0);
        Substructure_Atom* root = const_cast<Substructure_Atom*>(q->root_atom(0));
        if (! root->matches(target_atom, already_matched)) {
          continue;
        }
        if (q->find_embedding(target, target_atom, matched_atoms,
                              already_matched, root, sresults)) {
          got_match = true;
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
          cerr << "MatchedAtomMatch::Matches:got positive match\n";
#endif
          break;
        }
      }
      if (! got_match) {
        return 0;
      }
    }

    for (Single_Substructure_Query* q : _negative_matches) {
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
      cerr << "Trying negative match\n";
#endif
      std::fill_n(already_matched, matoms, 0);
      Substructure_Atom* root = const_cast<Substructure_Atom*>(q->root_atom(0));
      if (! root->matches(target_atom, already_matched)) {
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
        cerr << "MatchedAtomMatch::Matches:No match to root atom\n";
#endif
        continue;
      }
      if (q->find_embedding(target, target_atom, matched_atoms,
                            already_matched, root, sresults)) {
#ifdef DEBUG_MATCHED_ATOMS_MATCH_ALSO_MATCHED
        cerr << "MatchedAtomMatch:Matches:got negative match\n";
#endif
        return 0;
      }
    }
  }

  return 1;
}

RequiredMolecularProperties::RequiredMolecularProperties() {
  _attributes_specified = 0;
  _any_net_formal_charge = 0;
}

#define INCREMENENT_RETURN_IF_DONE(attributes_checked, attributes_specified) { \
  ++attributes_checked; \
  if (attributes_checked == attributes_specified) { \
    return 1; \
  } \
}

// Code copied from Single_Substructure_Query.
int
RequiredMolecularProperties::Matches(Molecule_to_Match & target_molecule) {
  int attributes_checked = 0;

  if (_natoms.is_set()) {
    if (! _natoms.matches(target_molecule.natoms())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check heteroatoms " << _heteroatoms_in_molecule.is_set() << '\n';
#endif

  if (_heteroatoms_in_molecule.is_set()) {
    if (! _heteroatoms_in_molecule.matches(target_molecule.heteroatom_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Do we need to check nrings " << _nrings.is_set() << '\n';
#endif

  if (_nrings.is_set()) {
    if (! _nrings.matches(target_molecule.nrings())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check isotopes " << _number_isotopic_atoms.is_set() << '\n';
#endif

  if (_number_isotopic_atoms.is_set()) {
    if (! _number_isotopic_atoms.matches(target_molecule.number_isotopic_atoms())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_number_fragments.is_set()) {
    if (! _number_fragments.matches(target_molecule.number_fragments())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check elements needed? " << _elements_needed.number_elements() << '\n';
#endif

  if (_elements_needed.number_elements()) {
    if (! MatchElementsNeeded(target_molecule)) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

#ifdef DEBUG_CHECK_GLOBAL_CONDITIONS
  cerr << "Check ring type specifications\n";
#endif

  if (_atoms_in_spinach.is_set() || _inter_ring_atoms.is_set()) {
    if (! SpinachAtomsMatch(target_molecule)) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_net_formal_charge.is_set()) {
    if (! _net_formal_charge.matches(target_molecule.net_formal_charge())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_aromatic_atoms.is_set()) {
    if (! _aromatic_atoms.matches(target_molecule.molecule()->aromatic_atom_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_aromatic_rings.is_set()) {
    if (! _aromatic_rings.matches(target_molecule.aromatic_ring_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_non_aromatic_rings.is_set()) {
    if (! _non_aromatic_rings.matches(target_molecule.non_aromatic_ring_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_isolated_rings.is_set()) {
    if (! _isolated_rings.matches(target_molecule.isolated_ring_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_ring_systems.is_set()) {
    if (! _ring_systems.matches(target_molecule.ring_object_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_fused_rings.is_set()) {
    if (! _fused_rings.matches(target_molecule.fused_ring_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  if (_strongly_fused_rings.is_set()) {
    if (! _strongly_fused_rings.matches(target_molecule.strongly_fused_ring_count())) {
      return 0;
    }
    INCREMENENT_RETURN_IF_DONE(attributes_checked, _attributes_specified);
  }

  cerr << "RequiredMolecularProperties::Matches: default return - should not happen\n";
  return 1;
}

int
RequiredMolecularProperties::SpinachAtomsMatch(Molecule_to_Match & target) const
{
  Molecule & m = *(target.molecule());

  int matoms = m.natoms();

  int * spinach = new int[matoms]; std::unique_ptr<int[]> free_spinach(spinach);

  int ais = m.identify_spinach(spinach);

//cerr << "Single_Substructure_Query::_spinach_atoms_match:ais " << ais << '\n';
  if (! _atoms_in_spinach.matches(ais)) {
    return 0;
  }

  if (! _inter_ring_atoms.is_set()) {
    return 1;
  }

  int number_inter_ring_atoms = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (spinach[i])    // inter ring atoms are not in the spinach
      continue;

    if (target[i].is_non_ring_atom())
      number_inter_ring_atoms++;
  }

//cerr << "Single_Substructure_Query::_spinach_atoms_match:number_inter_ring_atoms " << number_inter_ring_atoms << '\n';
  return _inter_ring_atoms.matches(number_inter_ring_atoms);
}
