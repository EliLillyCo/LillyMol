#include <iostream>
#include <memory>
using std::cerr;
using std::endl;


#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "chiral_centre.h"
#include "iwreaction.h"
#include "mdl_molecule.h"

static void
add_changes(const resizable_array<int> & stuff,
            int * v)
{
  for (int i = 0; i < stuff.number_elements(); i++)
  {
    int j = stuff[i];

    v[j]++;
  }

  return;
}

int
Reaction_Site::determine_which_matched_atoms_are_changed ()
{
  assert (nullptr == _matched_atom_changed);

  int h = highest_initial_atom_number();

  if (h <= 0)
  {
    cerr << "Reaction_Site::determine_which_matched_atoms_are_changed:no query, assuming 100 atoms\n";
    h = 100;
  }

  _matched_atom_changed = new_int(h + 1);

  for (int i = 0; i < _bonds_to_be_made.number_elements(); i++)
  {
    const Bond * b = _bonds_to_be_made[i];

    _matched_atom_changed[b->a1()]++;
    _matched_atom_changed[b->a2()]++;
  }

  for (int i = 0; i < _bonds_to_be_broken.number_elements(); i++)
  {
    const Bond * b = _bonds_to_be_broken[i];

    _matched_atom_changed[b->a1()]++;
    _matched_atom_changed[b->a2()]++;
  }

  add_changes(_atoms_to_be_removed, _matched_atom_changed);

  add_changes(_fragments_to_be_removed, _matched_atom_changed);

  for (int i = 0; i < _elements_to_change.number_elements(); i++)
  {
    const Reaction_Change_Element * e = _elements_to_change[i];

    int j = e->atom();

    _matched_atom_changed[j]++;
  }

  for (int i = 0; i < _formal_charges_to_assign.number_elements(); i++)
  {
    const Reaction_Formal_Charge * e = _formal_charges_to_assign[i];

    int j = e->atom();

    _matched_atom_changed[j]++;
  }

  for (int i = 0; i < _formal_charges_to_change.number_elements(); i++)
  {
    const Reaction_Change_Formal_Charge * e = _formal_charges_to_change[i];

    int j = e->atom();

    _matched_atom_changed[j]++;
  }

  for (const Reaction_Place_Isotope* iso : _isotopes_to_assign) {
    for (int atom : iso->matched_atoms()) {
      _matched_atom_changed[atom]++;
    }
  }

  for (int i = 0; i < _replace_atom.number_elements(); i++)
  {
    const Replace_Atom * e = _replace_atom[i];

    const Matched_Atom_in_Component & ma1 = e->a1();

    int j = ma1.atom();

    _matched_atom_changed[j]++;

    const Matched_Atom_in_Component & ma2 = e->a2();

     j = ma2.atom();

    _matched_atom_changed[j]++;
  }

  for (int i = 0; i < _wedge_bonds_to_place.number_elements(); i++)
  {
    const Reaction_Wedge_Bond * e = _wedge_bonds_to_place[i];

    int j = e->a1();

    _matched_atom_changed[j]++;

    j = e->a2();

    _matched_atom_changed[j]++;
  }

  add_changes(_stereo_centres_to_invert, _matched_atom_changed);

  add_changes(_chiral_centres_to_remove, _matched_atom_changed);

  return 1;
}

/*
  Some attributes are not part of the inherited Reaction_Site class
*/

int
IWReaction::_determine_which_matched_atoms_are_changed ()
{
  for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++)
  {
    const Reaction_Stereo_Centre * e = _reaction_stereo_centre[i];

    int j = e->atom();

    _matched_atom_changed[j]++;
  }

  return 1;
}

int
Reaction_Site::another_reagent_changes_your_matched_atom (int a)
{
  assert (nullptr != _matched_atom_changed);

  assert (a >= 0);

  _matched_atom_changed[a]++;

  return 1;
}

int
Reaction_Site::_remove_multiple_hits_that_do_not_involve_changing_atoms (Molecule & m,
                                        Substructure_Results & sresults) const
{
  assert (nullptr != _matched_atom_changed);

//#define DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
  cerr << "Checking " << sresults.number_embeddings() << " embeddings for non changing atoms\n";
#endif

  const int * sym = m.symmetry_classes();

  for (int i = sresults.number_embeddings() - 1; i > 0; i--)
  {
    const Set_of_Atoms * ei = sresults.embedding(i);

    for (int j = 0; j < i; j++)
    {
      const Set_of_Atoms * ej = sresults.embedding(j);

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
      cerr << "What about " << (*ei) << " and " << (*ej) << endl;
#endif

      if (_common_matched_atoms_do_not_change(*ei, *ej, sym))
      {
        sresults.remove_embedding(i);
        break;
      }
    }
  }

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
  cerr << "Final " << sresults.number_embeddings() << endl;
#endif

  return 1;
}

int
Reaction_Site::_common_matched_atoms_do_not_change (const Set_of_Atoms & e1,
                                                    const Set_of_Atoms & e2,
                                                    const int * symm) const
{
#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
  cerr << "Examining embeddings " << e1 << " and " << e2 << endl;
#endif

  int n = e1.number_elements();

  if (n > e2.number_elements())    // different sized embeddings, hopefully unusual
    n = e2.number_elements();

  int same_atom_found = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j1 = e1[i];
    atom_number_t j2 = e2[i];

    if (j1 == j2)    // great, same atom in both embeddings
    {
      same_atom_found++;
      continue;
    }
      
//  Embeddings have a different atom at position I

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_DO_NOT_INVOLVE_CHANGING_ATOMS
    cerr << "Atoms " << j1 << " and " << j2 << " changed? " << _matched_atom_changed[i] << " symmetry " << symm[j1] << " and " << symm[j2] << endl;
#endif

    if (0 == _matched_atom_changed[i])    // atom not changed by the reaction, doesn't matter
      continue;

//  if (symm[j1] == symm[j2])    // case of any di-phenol (o, m or p), make sure we get just want one match
//  {
//    same_atom_found++;
//    continue;
//  }

    return 0;
  }

  if (0 == same_atom_found)     // embeddings don't overlap at all
    return 0;

  return 1;   // never found a different atom where things were changing
}

int
IWReaction::setup_to_skip_multiple_embeddings_involving_non_changing_atoms ()
{
  Reaction_Site::determine_which_matched_atoms_are_changed();

  _determine_which_matched_atoms_are_changed();

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    _sidechains[i]->determine_which_matched_atoms_are_changed();
  }

  for (int i = 0; i < ns; i++)
  {
    _sidechains[i]->notify_scaffold_of_atoms_involved_in_changes(*this);
  }

  for (int i = 0; i < ns; i++)
  {
    for (int j = 0; j < ns; j++)
    {
      if (i == j)
        continue;

      _sidechains[j]->notify_sidechain_of_atoms_involved_in_changes(_sidechains[i]);
    }
  }

//#define ECHO_CHANGED_ATOMS
#ifdef ECHO_CHANGED_ATOMS
  cerr << "IWReaction::setup_to_skip_multiple_embeddings_involving_non_changing_atoms\n";
  debug_print(cerr);
  for (int i = 0; i < ns; i++)
  {
    _sidechains[i]->debug_print(cerr, "  ");
  }
#endif

  return 1;
}

/*
  Doing Replace_Atom and Inter_Particle_Bond notifications look the same, so put
  them in a template
*/

template <typename T>
int
notify_scaffold_of_atom_changes (Reaction_Site & r,
                                 const resizable_array_p<T> & c)
{
#ifdef DEBUG_NOTIFY_SCAFFOLD_OF_ATOM_CHANGES
  cerr << "Notifying scaffold of " << c.number_elements() << " possible changes\n";
#endif

  for (int i = 0; i < c.number_elements(); i++)
  {
    const T * a = c[i];

    const Matched_Atom_in_Component & a1 = a->a1();

#ifdef DEBUG_NOTIFY_SCAFFOLD_OF_ATOM_CHANGES
    cerr << "A1 in scaffold " << a1.in_scaffold() << endl;
#endif

    if (a1.in_scaffold())
    {
      int j = a1.matched_atom();
      r.another_reagent_changes_your_matched_atom(j);
    }

    const Matched_Atom_in_Component & a2 = a->a2();

#ifdef DEBUG_NOTIFY_SCAFFOLD_OF_ATOM_CHANGES
    cerr << "A2 in scaffold " << a2.in_scaffold() << endl;
#endif

    if (a2.in_scaffold())
    {
      int j = a2.matched_atom();
      r.another_reagent_changes_your_matched_atom(j);
    }
  }

  return 1;
}
                        
template int notify_scaffold_of_atom_changes(Reaction_Site &, const resizable_array_p<Replace_Atom> &);
template int notify_scaffold_of_atom_changes(Reaction_Site &, const resizable_array_p<Inter_Particle_Bond> &);

int
Sidechain_Reaction_Site::notify_scaffold_of_atoms_involved_in_changes (Reaction_Site & r) const
{
  notify_scaffold_of_atom_changes(r, _replace_atom);
  notify_scaffold_of_atom_changes(r, _inter_particle_bonds);

  return 1;
}

template <typename T>
int
notify_sidechain_of_atom_changes (Sidechain_Reaction_Site * r,
                                  const resizable_array_p<T> & c)
{
  for (int i = 0; i < c.number_elements(); i++)
  {
    const T * a = c[i];

    const Matched_Atom_in_Component & a1 = a->a1();
    if (a1.in_component() == r->unique_id())
    {
      int j = a1.matched_atom();
      r->another_reagent_changes_your_matched_atom(j);
    }

    const Matched_Atom_in_Component & a2 = a->a2();
    if (a2.in_component() == r->unique_id())
    {
      int j = a2.matched_atom();
      r->another_reagent_changes_your_matched_atom(j);
    }
  }

  return 1;
}

template int notify_sidechain_of_atom_changes(Sidechain_Reaction_Site *, const resizable_array_p<Replace_Atom> &);
template int notify_sidechain_of_atom_changes(Sidechain_Reaction_Site *, const resizable_array_p<Inter_Particle_Bond> &);

int
Sidechain_Reaction_Site::notify_sidechain_of_atoms_involved_in_changes (Sidechain_Reaction_Site * r) const
{
  notify_sidechain_of_atom_changes(r, _replace_atom);
  notify_sidechain_of_atom_changes(r, _inter_particle_bonds);

  return 1;
}

int
Reaction_Site::_remove_multiple_hits_that_hit_atoms_being_changed (int matoms,
                                                Substructure_Results & sresults) const
{
  assert (nullptr != _matched_atom_changed);

  int * tmp = new int[matoms]; std::unique_ptr<int[]> freetmp(tmp);

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_HIT_ATOMS_BEING_CHANGED
  cerr << "Trimming " << sresults.number_embeddings() << " hits in molecule with " << matoms << " atoms\n";
#endif

  for (int i = sresults.number_embeddings() - 1; i > 0; i--)
  {
    const Set_of_Atoms * ei = sresults.embedding(i);

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_HIT_ATOMS_BEING_CHANGED
    cerr << "Checking embedding " << i << " " << (*ei) << endl;
#endif

    set_vector(tmp, matoms, 0);

    int n = ei->number_elements();

    for (int j = 0; j < n; j++)    // mark all the atoms that EI will change
    {
//    cerr << "Matched atom " << j << " changed? " << _matched_atom_changed[j] << endl;

      if (0 == _matched_atom_changed[j])     // the J'th member of the embedding is never changed
        continue;

      atom_number_t k = ei->item(j);

      tmp[k] = 1;

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_HIT_ATOMS_BEING_CHANGED
      cerr << "Atom " << k << " changes\n";
#endif
    }

    for (int j = 0; j < i; j++)
    {
      const Set_of_Atoms * ej = sresults.embedding(j);

#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_HIT_ATOMS_BEING_CHANGED
      cerr << "Compare " << (*ej) << endl;
#endif

      if (_common_matched_atoms_change(tmp, *ej))
      {
#ifdef DEBUG_REMOVE_MULTIPLE_HITS_THAT_HIT_ATOMS_BEING_CHANGED
        cerr << "Embedding " << i << " being removed\n";
#endif

        sresults.remove_embedding(i);
        break;
      }
    }
  }

  return 1;
}

/*
  Do these two embeddings hit any atoms that are changed
*/

int
Reaction_Site::_common_matched_atoms_change (const int * changed_by_other_embedding,
                                             const Set_of_Atoms & e) const
{
  int n = e.number_elements();

//cerr << "Embedding has " << n << " members\n";

  for (int i = 0; i < n; i++)
  {
//  cerr << " i = " << i << " changed? " << _matched_atom_changed[i] << endl;

    if (0 == _matched_atom_changed[i])    // we are looking for members of the embeddings that change
      continue;

    atom_number_t j = e[i];

    if (changed_by_other_embedding[j])
      return 1;
  }

  return 0;   // never found a case where the common atoms were being changed
}

/*
  We want to discard any matches in srtmp that are not in sresults
*/

static int
keep_only_hits_in_original_set (const Substructure_Results & sresults,
                                int matoms,
                                int * tmp,
                                Substructure_Results & srtmp)
{
  int ne1 = sresults.number_embeddings();
  int ne2 = srtmp.number_embeddings();

  for (int i = ne2 - 1; i >= 0; i--)
  {
    set_vector(tmp, matoms, 0);
    srtmp.embedding(i)->set_vector(tmp, 1);

    int found_match = 0;
    for (int j = 0; j < ne1; j++)
    {
      if (sresults.embedding(j)->all_members_set_in_array(tmp, 1))
      {
        found_match = 1;
        break;
      }
    }

    if (! found_match)
      srtmp.remove_embedding(i);
  }

  return srtmp.number_embeddings();
}

//#define DEBUG__PERFORM_REACTION_RECHECK_SCAFFOLD_REACTIVITY

/*
  This is very tricky because there may be atoms or fragments removed from the scaffold
*/

int
IWReaction::_perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold,
                                                           const Substructure_Results & sresults,
                                                           Molecule & result,
                                                           Enumeration_Temporaries & etmp)
{
  _add_molecule(result, *scaffold);

  if (! _perform_reaction(result, sresults.embedding(0), etmp))
    return 0;

  _break_bonds_to_atoms_and_fragments_to_be_removed(result, etmp);

  int initial_natoms = scaffold->natoms();

  int * tmp = new int[initial_natoms]; std::unique_ptr<int[]> free_tmp(tmp);

  int ne = sresults.number_embeddings();

#ifdef DEBUG__PERFORM_REACTION_RECHECK_SCAFFOLD_REACTIVITY
  cerr << "Checking " << ne << " embeddings\n";
#endif

  for (int i = 1; i < ne; i++)
  {
    Substructure_Results srtmp;

    int nhits = substructure_search(result, srtmp);

#ifdef DEBUG__PERFORM_REACTION_RECHECK_SCAFFOLD_REACTIVITY
    cerr << " iteration " << i << " found " << nhits << " hits in '" << result.smiles() << "'\n";
#endif

    if (0 == nhits)
      break;

    if (0 == keep_only_hits_in_original_set(sresults, initial_natoms, tmp, srtmp))
      continue;

    if (! _perform_reaction(result, srtmp.embedding(0), etmp))
      return 0;

    _break_bonds_to_atoms_and_fragments_to_be_removed(result, etmp);
  }

  return _do_atom_removals(result, etmp);
}

int
IWReaction::perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold,
                                                          const Substructure_Results & sresults,
                                                          const Reaction_Iterator & iter,
                                                          Molecule & result)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  _copy_iterator_data_to_etmp(iter, etmp);

  return _perform_reaction_recheck_scaffold_reactivity(scaffold, sresults, result, etmp);
}

int
IWReaction::perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold,
                                                          const Substructure_Results & sresults,
                                                          Molecule & result)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  _take_first_reagent_from_each_sidechain(etmp);

  return _perform_reaction_recheck_scaffold_reactivity(scaffold, sresults, result, etmp);
}

int
IWReaction::_break_bonds_to_atoms_and_fragments_to_be_removed (Molecule & result,
                                const Enumeration_Temporaries & etmp) const
{
  const Set_of_Atoms & e = etmp.atoms_to_be_removed();

  int n = e.number_elements();

  for (int i = 0; i < n; i++)
  {
    result.remove_bonds_to_atom(e[i], 1);    // 2nd arg means preserve chirality

//  cerr << "After removing bonds to '" << e[i] << "' result is " << result.smiles() << "'\n";
  }

  return 1;
  const resizable_array<const Atom *> & rmf = etmp.fragments_to_be_removed();

  n = rmf.number_elements();

  for (int i = 0; i < n; i++)
  {
    cerr << "Implement fragment breakage, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
  }

  return 1;
}

static int
read_record_from_file_of_reactions (const const_IWSubstring & fname, 
                        const const_IWSubstring & dirname,
                        resizable_array_p<IWReaction> & rxn,
                        Sidechain_Match_Conditions & smc)
{
  IWReaction * tmp = new IWReaction;

  int rc;

  if (fname.starts_with(IWDIRECTORY_SEPARATOR))
    rc = tmp->do_read(fname, smc);
  else
  {
    IWString full_path_name;
    full_path_name << dirname << IWDIRECTORY_SEPARATOR << fname;
    rc = tmp->do_read(full_path_name, smc);
  }

  if (0 == rc)
  {
    cerr << "Cannot read reaction '" << fname << "'\n";
    delete tmp;
    return 0;
  }

  rxn.add(tmp);

  return 1;
}

static int
read_file_of_reactions (iwstring_data_source & input,
                        const const_IWSubstring & dirname,
                        resizable_array_p<IWReaction> & rxn,
                        Sidechain_Match_Conditions & smc)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    buffer.strip_trailing_blanks();

    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    if (! read_record_from_file_of_reactions(buffer, dirname, rxn, smc))
    {
      cerr << "Invalid reaction file '" << buffer << "' line " << input.lines_read() << endl;
      return 0;
    }
  }

  return rxn.number_elements();
}

static int
read_file_of_reactions (const const_IWSubstring & fname, 
                        resizable_array_p<IWReaction> & rxn,
                        Sidechain_Match_Conditions & smc)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "read_file_of_reactions:cannot open '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring dirname(fname);
  int i = dirname.rindex(IWDIRECTORY_SEPARATOR);
  if (i > 0)
    dirname.iwtruncate(i);

  return read_file_of_reactions(input, dirname, rxn, smc);
}

int
read_reactions (const Command_Line & cl,
                resizable_array_p<IWReaction> & rxn,
                Sidechain_Match_Conditions & smc,
                char flag)
{
  int i = 0;
  const_IWSubstring r;

  while (cl.value(flag, r, i++))
  {
    if (r.starts_with("F:"))
    {
      r.remove_leading_chars(2);
      if (! read_file_of_reactions(r, rxn, smc))
      {
        cerr << "Cannot read file of reactons '" << r << "'\n";
        return 0;
      }
    }
    else
    {
      IWReaction * tmp = new IWReaction;
      if (! tmp->do_read(r, smc))
      {
        cerr << "Cannot read reaction '" << r << "'\n";
        return 0;
      }

      rxn.add(tmp);
    }
  }

  return rxn.number_elements();
}

static int
count_number_of_reactions (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return input.records_remaining();
}

int
count_number_of_reactions (const Command_Line & cl,
                           char flag)
{
  int rc = 0;

  int i = 0;
  const_IWSubstring r;

  while (cl.value(flag, r, i++))
  {
    if (r.starts_with("F:"))
    {
      r.remove_leading_chars(2);

      int tmp = count_number_of_reactions(r);
      if (0 == tmp)
      {
        cerr << "Missing or empty file '" << r << "'\n";
        return 0;
      }

      rc += tmp;
    }
    else if (r.starts_with("SMIRKS:"))
      rc++;
    else if (r.starts_with("K:"))
    {
      r.remove_leading_chars(2);


      int tmp = count_number_of_reactions(r);
      if (0 == tmp)
      {
        cerr << "Missing or empty smirks file '" << r << "'\n";
        return 0;
      }

      rc += tmp;
    }
    else
      rc++;
  }

  return rc;
}

int
write_isotopically_labelled_smiles(Molecule & m,
                                   std::ostream & os)
{
  int matoms = m.natoms();

  isotope_t * isosave = new isotope_t[matoms]; std::unique_ptr<isotope_t[]> free_isosave(isosave);

  m.get_isotopes(isosave);

  for (int i = 0; i < matoms; i++)
  {
    m.set_isotope(i, i);
  }

  os << m.smiles() << ' ' << m.name() << '\n';

  m.set_isotopes(isosave);

  return 1;
}

const Substructure_Atom *
IWReaction::query_atom_with_initial_atom_number(atom_number_t z) const 
{
  const Substructure_Query & q = *this;
  const Substructure_Atom * rc = q.query_atom_with_initial_atom_number(z);

  if (nullptr != rc)
    return rc;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    const Substructure_Atom * rc = _sidechains[i]->query_atom_with_initial_atom_number(z);
    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

Reaction_Site *
IWReaction::_reaction_site_with_initial_atom_number(atom_number_t z)
{
  const Substructure_Atom * a = this->query_atom_with_initial_atom_number(z);
  if (nullptr != a)
    return this;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    const Substructure_Atom * a = _sidechains[i]->query_atom_with_initial_atom_number(z);

    if (nullptr != a)
      return _sidechains[i];
  }

  return nullptr;
}

Reaction_Site *
IWReaction::_reaction_site_with_atom_map_number(const int amap)
{
  const Substructure_Atom * a = this->query_atom_with_atom_map_number(amap);
  if (nullptr != a)
    return this;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    const Substructure_Atom * a = _sidechains[i]->query_atom_with_atom_map_number(amap);

    if (nullptr != a)
      return _sidechains[i];
  }

  return nullptr;
}

const Substructure_Bond *
IWReaction::bond_between_atoms (atom_number_t a1, atom_number_t a2) const
{
//cerr << "IWReaction::bond_between_atoms:looking for bond between " << a1 << " and " << a2 << endl;
  const Substructure_Query & q = *this;

  const Substructure_Bond * rc = q.bond_between_atoms(a1, a2);

  if (nullptr != rc)
    return rc;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    rc = _sidechains[i]->bond_between_atoms(a1, a2);
    if (nullptr != rc)
      return rc;
  }

  return rc;
}

Reaction_Site *
IWReaction::_component_with_bond (const atom_number_t a1, const atom_number_t a2)
{
  const Substructure_Query & q = *this;

  if (nullptr != q.bond_between_atoms(a1, a2))
    return this;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    if (nullptr != _sidechains[i]->bond_between_atoms(a1, a2))
      return _sidechains[i];
  }

  return nullptr;
}

Reaction_Site *
IWReaction::_component_with_bond_between_mapped_atoms (const int a1, const int a2)
{
  const Substructure_Query & q = *this;

  if (nullptr != q.bond_between_atom_map_numbers(a1, a2))
    return this;

  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    if (nullptr != _sidechains[i]->bond_between_atom_map_numbers(a1, a2))
      return _sidechains[i];
  }

  return nullptr;
}

int
IWReaction::_sidechain_with_mapped_atom(const int x) const
{
  for (unsigned int i = 0; i < _sidechains.size(); ++i)
  {
    if (nullptr != _sidechains[i]->query_atom_with_atom_map_number(x))
      return i;
  }

  return -1;
}

int
IWReaction::add_sidechain_reagents(int sidechain, const char* fname, FileType input_type,
                         const Sidechain_Match_Conditions& smc) {
  return _sidechains[sidechain]->add_reagents(fname, input_type, smc); }
