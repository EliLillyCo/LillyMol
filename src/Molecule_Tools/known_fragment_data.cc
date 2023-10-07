#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/smiles.h"
#include "known_fragment_data.h"

//#define FD_MAX_NATOMS 100

Known_Fragment_Data::Known_Fragment_Data()
{
  _only_check_molecular_formula = 0;

  _remove_everything_if_all_fragments_match = 1;

  return;
}

Known_Fragment_Data::~Known_Fragment_Data()
{
  return;
}

int
Known_Fragment_Data::debug_print(std::ostream & os) const
{
  os << "Fragment data contains " << _known_salt_mf.size() << " known salts MF and " << _known_parent_mf.size() << " known parent MF fragments\n";

  os << "Known salts MF\n";
  _common_print_set(_known_salt_mf, os);
  os << "Known salts USMI\n";
  _common_print_set(_known_salt_usmi, os);

  if (_known_parent_mf.size() > 0)
  {
    os << "Known parents MF\n";
    _common_print_set(_known_salt_mf, os);
    os << "Known parents USMI\n";
    _common_print_set(_known_salt_usmi, os);
  }

  return os.good();
}

int
Known_Fragment_Data::_common_print_set(const _formula_usmi & h,
                                        std::ostream & os) const
{
  for (_formula_usmi::const_iterator i = h.begin(); i != h.end(); ++i)
  {
    os << " '" << (*i) << "'\n";
  }

  return os.good();
}

int
Known_Fragment_Data::read_known_salts(const const_IWSubstring & fname)
{
  return _common_read(fname, _known_salt_mf, _known_salt_usmi);
}

int
Known_Fragment_Data::read_known_parents(const const_IWSubstring & fname)
{
  return _common_read(fname, _known_parent_mf, _known_parent_usmi);
}

int
Known_Fragment_Data::_common_read(const const_IWSubstring & fname,
                             _formula_usmi & mf,
                             _formula_usmi & usmi)
{
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);

  if (! input.good())
  {
    cerr << "Known_Fragment_Data::read_known_salts:cannot open '" << fname << "'\n";
    return 0;
  }

  return _common_read(input, mf, usmi);
}

int
Known_Fragment_Data::_common_read(data_source_and_type<Molecule> & input,
                             _formula_usmi & mf,
                             _formula_usmi & usmi)
{
  Temporarily_Set_Include_Chiral_Info_in_Smiles tsics(0);

  Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    if (! m->valence_ok())
      cerr << "Warning, invalid valence '" << m->name() << "'\n";

    if (! _add_to_hash(*m, mf, usmi))
      return 0;
  }

  return static_cast<int>(mf.size());
}

int
Known_Fragment_Data::_add_to_hash(Molecule & m,
                             _formula_usmi & mf,
                             _formula_usmi & usmi)
{
  _natoms[m.natoms()]++;

  IWString tmp;
  m.isis_like_molecular_formula_dot_between_fragments(tmp);
  mf.insert(tmp);

  usmi.insert(m.unique_smiles());

  return 1;
}

static int
just_acid_and_carbon(Molecule & m,
                      const int * fragment_membership,
                      int f)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (f != fragment_membership[i])
      continue;

    if (m.is_ring_atom(i))
      return 0;

    atomic_number_t z = m.atomic_number(i);

    if (6 == z || 8 == z || 15 == z || 16 == z)
      ;
    else
      return 0;
  }

// First scan for an acid group. Ignore the carbons

  int * already_done = new_int(matoms);std::unique_ptr<int> free_already_done(already_done);

  int found_acid = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (f != fragment_membership[i])
      continue;

    if (already_done[i])
      continue;

    const Atom * a = m.atomi(i);

    if (1 != a->ncon())
      continue;

    if (2 != a->nbonds())
      continue;

    atomic_number_t z = a->atomic_number();

    if (8 == z || 16 == z)
      ;
    else
      continue;

    atom_number_t csp = a->other(i, 0);    // The Carbon, Sulphur or Phosphorus atom at the centre of the acid

    const Atom * acsp = m.atomi(csp);

    atomic_number_t zcsp = acsp->atomic_number();

    if (6 == zcsp || 15 == zcsp || 16 == zcsp)    // check probably reduncant
      ;
    else
      return 0;

    already_done[i] = 1;
    already_done[csp] = 1;

    for (int j = 0; j < acsp->ncon(); j++)
    {
      const Bond * b = acsp->item(j);

      atom_number_t k = b->other(csp);

      if (k == i)
        continue;

      const Atom * ak = m.atomi(k);

      if (6 == ak->atomic_number() && ak->ncon() == ak->nbonds() && ak->ncon() < 3)
        ;
      else if (8 == ak->atomic_number())
        ;
      else if (16 == ak->atomic_number())
        ;
      else
        return 0;

      already_done[k] = 1;
    }

    found_acid++;
  }

  if (0 == found_acid)
    return 0;

// The remaining atoms must be carbons

  for (int i = 0; i < matoms; i++)
  {
    if (f != fragment_membership[i])
      continue;

    if (already_done[i])
      continue;

    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      return 0;

    if (a->ncon() > 2)
      return 0;

    if (a->nbonds() > a->ncon())
      return 0;
  }

  return 1;
}

static int
just_nitrogen_and_carbon(Molecule & m,
                          const int * fragment_membership,
                          int f)
{
  int number_4_connected_nitrogen = 0;
  int number_3_connected_nitrogen = 0;
  int number_2_connected_nitrogen = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (f != fragment_membership[i])
      continue;

    const Atom * a = m.atomi(i);

    atomic_number_t z = a->atomic_number();

    if (a->nbonds() > a->ncon())
      return 0;

    if (6 == z)
    {
      if (a->ncon() > 2 || m.is_ring_atom(i))
        return 0;
    }
    else if (7 == z)
    {
      if (4 == a->ncon())
        number_4_connected_nitrogen++;
      else if (3 == a->ncon())
        number_3_connected_nitrogen++;
      else if (2 == a->ncon())
        number_2_connected_nitrogen++;
      else
        return 0;
    }
    else
      return 0;
  }

  if (number_4_connected_nitrogen > 0 && 0 == number_3_connected_nitrogen)
    return 1;

  if (0 == number_4_connected_nitrogen && 1 == number_3_connected_nitrogen)
    return 1;

  if (0 == number_4_connected_nitrogen && 0 == number_3_connected_nitrogen && 1 == number_2_connected_nitrogen)
    return 1;

  return 0;
}

int
Known_Fragment_Data::_remove_non_organic(Molecule & m) const
{
  resizable_array<int> fragments_to_be_removed;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Element * e = m.elementi(i);

    if (e->organic())
      continue;

    int f = m.fragment_membership(i);

    fragments_to_be_removed.add_if_not_already_present(f);
  }

  if (fragments_to_be_removed.empty())
    return 0;

  return _delete_set_of_fragments(m, fragments_to_be_removed);
}

int
Known_Fragment_Data::_remove_soap(Molecule & m) const
{
  int * fragment_membership = new int[m.natoms()]; std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  m.fragment_membership(fragment_membership);

  int nf = m.number_fragments();

  resizable_array<int> fragments_to_be_removed;

  for (int i = 0; i < nf; i++)
  {
    if (just_acid_and_carbon(m, fragment_membership, i))
      fragments_to_be_removed.add(i);
    else if (just_nitrogen_and_carbon(m, fragment_membership, i))
      fragments_to_be_removed.add(i);
  }

//if (fragments_to_be_removed.number_elements())
//  cerr << m.smiles() << ' ' << m.name() << " SOAP or QUAT\n";

  if (fragments_to_be_removed.empty())
    return 0;

  return _delete_set_of_fragments(m, fragments_to_be_removed);
}

int
Known_Fragment_Data::process(Molecule & m)
{
  _remove_soap(m);

  const int nf = m.number_fragments();

  if (nf <= 1)
    return 1;

  _remove_non_organic(m);

  if (m.number_fragments() <= 1)
    return 1;

  resizable_array_p<Molecule> fragments;

  IWString mf;
  m.isis_like_molecular_formula_dot_between_fragments(mf);

  resizable_array<int> fragments_to_be_removed;

  _scan_known_fragments(m, mf, fragments, _known_salt_mf, _known_salt_usmi, fragments_to_be_removed);

  if (fragments_to_be_removed.number_elements() + 1 == fragments.number_elements())   // only one fragment NOT being removed
    return m.delete_fragments(fragments_to_be_removed);

  if (fragments_to_be_removed.number_elements() == m.number_fragments())   // everything being removed
  {
    if (_remove_everything_if_all_fragments_match)
      return m.resize(0);

    return 0;
  }

  int atoms_being_removed = 0;
  for (int i = 0; i < nf; ++i)
  {
    if (fragments_to_be_removed.contains(i))
      atoms_being_removed += m.atoms_in_fragment(i);
  }

  if ((m.natoms() - atoms_being_removed) <= 2)    // arbitrary decision. If we will be left with two or fewer atoms, that cannot be right
  {
    return 0;
  }

  resizable_array<int> fragments_to_keep;

  _scan_known_fragments(m, mf, fragments, _known_parent_mf, _known_parent_usmi, fragments_to_keep);

  if (fragments_to_keep.empty())
  {
    if (fragments_to_be_removed.number_elements())
      return m.delete_fragments(fragments_to_be_removed);

    return 0;
  }

  return 1;
}

int
Known_Fragment_Data::_scan_known_fragments(Molecule & m,
                                            const IWString & molecular_formula,
                                            resizable_array_p<Molecule> & fragments,
                                            _formula_usmi & mf,
                                            _formula_usmi & usmi,
                                            resizable_array<int> & fragments_identified) const
{
  int i = 0;
  IWString token;
  int frag = 0;
  while (molecular_formula.nextword(token, i, '.'))    // scan tokens in the molecular formula
  {
//  cerr << "Scanning molecular formula token '" << token << "'\n";
    if (! mf.contains(token))
    {
      frag++;
      continue;
    }

    if (_only_check_molecular_formula)
    {
      fragments_identified.add(frag);
      frag++;
      continue;
    }

//  Need to check the unique smiles

    if (fragments.empty())   // something has the formula, need to do more...
      m.create_components(fragments);

    Molecule * fi = fragments[frag];

    if (1 == fi->natoms())    // just one atom, formula match OK
    {
      fragments_identified.add(frag);
      frag++;
      continue;
    }

    if (0 == _natoms[fi->natoms()])    // no molecules with this number of atoms in the set
      continue;

    Temporarily_Set_Include_Chiral_Info_in_Smiles tsics(0);

    const IWString & unique_smiles = fi->unique_smiles();

    if (usmi.contains(unique_smiles))
      fragments_identified.add(frag);
    
    frag++;
  }

  return fragments_identified.number_elements();
}

/*
  There is some ambiguity here if we have been instructed to NOT _remove_everything_if_all_fragments_match
  What should we do in that case? Right now, we leave the molecule alone
*/

int
Known_Fragment_Data::_delete_set_of_fragments(Molecule & m,
                                              const resizable_array<int> & fragments_to_be_removed) const
{
  if (fragments_to_be_removed.number_elements() == m.number_fragments())
  {
    if (_remove_everything_if_all_fragments_match)
      m.resize(0);
  }
  else
    m.delete_fragments(fragments_to_be_removed);

  return 1;
}

void
Known_Fragment_Data::deactivate()
{
  _only_check_molecular_formula = 0;

  _remove_everything_if_all_fragments_match = 1;

  _known_salt_mf.clear();
  _known_salt_usmi.clear();
  _known_parent_mf.clear();
  _known_parent_usmi.clear();

  return;
}
