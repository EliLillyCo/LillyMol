#include <stdlib.h>
#include <memory>

#define COMPILING_CAREFUL_FRAG

#include "Foundational/iwmisc/misc.h"

#include "molecule.h"

/*
  Fileconv has the ability to reduce to a fragment based on various criteria

  Need some concept of which elements are more desirable
*/

static int organic_desirability[HIGHEST_ATOMIC_NUMBER] = {
  0,
  0,    // 1 Hydrogen 
  0,
  0,
  0,
  0,
  3,    // 6 Carbon
  5,    // 7 Nitrogen
  4,    // 8 Oxygen
  3,    // 9 Fluorine
  0,
  0,
  0,
  0,
  0,
  1,    // 15 Phosphorus
  2,    // 16 Sulphur
  2,    // 17 Chlorine
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  2,    // 35 Bromine
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  2,    // 53 Iodine
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
};

int
Molecule::_identify_fragment_undesirable_groups(int * exclude) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (exclude[i])
      continue;

    const Atom * a = _things[i];

    if (7 == a->atomic_number())
    {
      _is_nitro(i, exclude);
    }
    else if (16 == a->atomic_number())
    {
      _is_sulphate_like(i, exclude);
    }
    else if (! a->element()->organic())
      exclude[i] = 2;
  }

  return 1;
}

int
Molecule::identify_largest_organic_fragment(Set_of_Atoms & atoms_to_be_removed,
                                            int & fragments_same_size_as_largest_organic)
{
  int nf = number_fragments();

  if (nf <= 1)
    return 1;

  int max_atoms_in_fragment = 0;
  int fragment_with_max_atoms = -1;
  int max_organic_atoms_in_fragment = 0;
  int fragment_with_most_organic_atoms = -1;

  fragments_same_size_as_largest_organic = 0;

  const int * fragment_membership = _fragment_information.fragment_membership();
  const resizable_array<int> & atoms_in_fragment = _fragment_information.atoms_in_fragment();

  int * exclude = new_int(_number_elements); std::unique_ptr<int[]> free_exclude(exclude);

  _identify_fragment_undesirable_groups(exclude);

  for (int i = 0; i < nf; i++)
  {
    if (atoms_in_fragment[i] > max_atoms_in_fragment)
    {
      max_atoms_in_fragment = atoms_in_fragment[i];
      fragment_with_max_atoms = i;
    }

    int organic_atoms_this_fragment = 0;
    int very_undesirable_this_fragment = 0;
    for (int j = 0; j < _number_elements; j++)
    {
//    cerr << "Atom " << j << " in fragment " << aj->fragment_membership() << endl;

      if (fragment_membership[j] != i)
        continue;

      if (exclude[j])
      {
        if (exclude[j] > 1)
          very_undesirable_this_fragment++;
        continue;
      }

      const Atom * aj = _things[j];

      if (! aj->element()->organic())
        continue;

      atomic_number_t z = aj->atomic_number();

      organic_atoms_this_fragment += organic_desirability[z];
    }

//  cerr << very_undesirable_this_fragment << " very_undesirable_this_fragment\n";

    if (very_undesirable_this_fragment)
      organic_atoms_this_fragment -= 10 * very_undesirable_this_fragment;

//  cerr << "organic_atoms_this_fragment " << organic_atoms_this_fragment << endl;

    if (organic_atoms_this_fragment > max_organic_atoms_in_fragment)
    {
      max_organic_atoms_in_fragment = organic_atoms_this_fragment;
      fragment_with_most_organic_atoms = i;
      fragments_same_size_as_largest_organic = 1;
    }
    else if (organic_atoms_this_fragment == max_organic_atoms_in_fragment)
      fragments_same_size_as_largest_organic++;
  }

//if (max_organic_atoms_in_fragment > 0 && fragments_same_size_as_largest_organic > 1)
//  cerr << fragments_same_size_as_largest_organic << " fragments size " << max_organic_atoms_in_fragment << ' ' << _molecule_name << endl;

  int fragment_to_keep;

  if (0 == max_organic_atoms_in_fragment)     // no fragments contain organic atoms
    fragment_to_keep = fragment_with_max_atoms;
  else
    fragment_to_keep = fragment_with_most_organic_atoms;

//cerr << "Keeping fragment " << fragment_to_keep << endl;
   
  atoms_to_be_removed.resize(max_atoms_in_fragment);

  for (int i = 0; i < _number_elements; i++)
  {
    if (fragment_membership[i] != fragment_to_keep)
      atoms_to_be_removed.add(i);
  }

//cerr << "Removing atoms " << atoms_to_be_removed << endl;

  return 1;
}

int
Molecule::reduce_to_largest_organic_fragment()
{
  Set_of_Atoms atoms_to_be_removed;

  int notused;

  (void) identify_largest_organic_fragment(atoms_to_be_removed, notused);

  if (atoms_to_be_removed.number_elements())
    remove_atoms(atoms_to_be_removed);

  return 1;
}

class Fragment_Data
{
  private:
    int _frag_id;
    int _organic_atoms;
    int _nitro;
    int _phosphorus;
    int _sulphur;
    int _sulphate;
    int _non_nitro_nitrogen;
    int _oxygen;
    int _hydrogens;
    int _non_organic_atoms;

  public:
    Fragment_Data();

    int debug_print(std::ostream &) const;

    int  frag_id() const { return _frag_id;}
    void set_frag_id(int f) { _frag_id = f;}

    void another_organic_atom() { _organic_atoms++;}
    int  organic_atoms() const { return _organic_atoms;}

    void another_nitro() { _nitro++;}
    void another_nitrogen() {_non_nitro_nitrogen++;}
    void another_sulphur() { _sulphur++;}
    void another_sulphate() { _sulphate++;}
    void another_oxygen() { _oxygen++;};
    void another_phosphorus() { _phosphorus++;};
    void extra_hydrogens(int h) { _hydrogens++;}
    void another_non_organic_atom() { _non_organic_atoms++;}

    int phosphorus_atoms() const { return _phosphorus;}
    int nitro_groups() const { return _nitro;}
    int nitrogen_atoms() const { return _non_nitro_nitrogen;}
    int oxygen_atoms() const { return _oxygen;}
    int sulphur_atoms() const { return _sulphur;}
    int hydrogen_atoms() const { return _hydrogens;}
    int non_organic_atoms() const { return _non_organic_atoms;}
};

Fragment_Data::Fragment_Data()
{
  _frag_id = -1;
  _organic_atoms = 0;
  _nitro = 0;
  _phosphorus = 0;
  _sulphur = 0;
  _sulphate = 0;
  _non_nitro_nitrogen = 0;
  _oxygen = 0;
  _hydrogens = 0;
  _non_organic_atoms = 0;

  return;
}

int
Fragment_Data::debug_print(std::ostream & os) const
{
  os << "Fragment data for fragment " << _frag_id << ", " << _organic_atoms << " organic atoms\n";
  os << "  " << _nitro << " nitro groups\n";
  os << "  " << _phosphorus << " phosphorus atoms\n";
  os << "  " << _sulphur << " sulphur atoms\n";
  os << "  " << _sulphate << " sulphate groups\n";
  os << "  " << _oxygen << " oxygen atoms\n";
  os << "  " << _non_nitro_nitrogen << " nitrogen atoms\n";
  os << "  " << _hydrogens << " hydrogen atoms\n";
  os << "  " << _non_organic_atoms << " non organic atoms\n";

  return os.good();
}

//#define DEBUG_CAREFUL_FRAG

int
Molecule::reduce_to_largest_fragment_carefully()
{
  int nf = number_fragments();
  
  if (1 == nf)
    return 1;

#ifdef DEBUG_CAREFUL_FRAG
  cerr << "Reducing '" << _molecule_name << "' from " << nf << " fragments\n";
#endif

  Set_of_Atoms atoms_to_be_removed;
  int number_instances_of_largest_fragment;

  identify_largest_organic_fragment(atoms_to_be_removed, number_instances_of_largest_fragment);

  if (atoms_to_be_removed.empty())
    return 1;

// If there are no organic fragments, NUMBER_INSTANCES_OF_LARGEST_FRAGMENT will be 0

  if (number_instances_of_largest_fragment <= 1)
  {
    remove_atoms(atoms_to_be_removed);

    return 1;
  }

// Now the more difficult case of multiple organic fragments with the
// same number of atoms.

  Fragment_Data * fd = new Fragment_Data[nf]; std::unique_ptr<Fragment_Data[]> free_fd(fd);

  int * already_counted = new_int(_number_elements); std::unique_ptr<int[]> free_already_counted(already_counted);

#ifdef DEBUG_CAREFUL_FRAG
  cerr << _molecule_name << " has " << number_instances_of_largest_fragment << " instances of largest fragment\n";
#endif

  return _reduce_to_largest_fragment_carefully(fd, already_counted);
}

static int
fragment_data_comparitor(const void * v1, const void * v2)
{
  const Fragment_Data * f1 = (const Fragment_Data *) v1;
  const Fragment_Data * f2 = (const Fragment_Data *) v2;

// We don't like non-organic atoms

  if (f1->non_organic_atoms() == f2->non_organic_atoms())   // most common case, both zero
    ;
  else if (f1->non_organic_atoms() > f2->non_organic_atoms())   
    return 1;
  else if (f1->non_organic_atoms() < f2->non_organic_atoms())   
    return -1;

  if (f1->organic_atoms() > f2->organic_atoms())
    return -1;

  if (f1->organic_atoms() < f2->organic_atoms())
    return 1;

// We don't like Phosphorus

  if (f1->phosphorus_atoms() == f2->phosphorus_atoms())   // most common case, both zero
    ;
  else if (f1->phosphorus_atoms() > f2->phosphorus_atoms())
    return 1;
  else if (f1->phosphorus_atoms() < f2->phosphorus_atoms())
    return -1;

// Definitely dis-favour nitro groups

  if (f1->nitro_groups() > f2->nitro_groups())
    return 1;

  if (f1->nitro_groups() < f2->nitro_groups())
    return -1;

// Sulphur isn't good either

  if (f1->sulphur_atoms() > f2->sulphur_atoms())
    return 1;

  if (f1->sulphur_atoms() < f2->sulphur_atoms())
    return -1;

// We like Nitrogen atoms

  if (f1->nitrogen_atoms() > f2->nitrogen_atoms())
    return -1;

  if (f1->nitrogen_atoms() < f2->nitrogen_atoms())
    return 1;

// We like oxygen atoms

  if (f1->oxygen_atoms() > f2->oxygen_atoms())
    return -1;

  if (f1->oxygen_atoms() < f2->oxygen_atoms())
    return 1;

// We like Hydrogens

  if (f1->hydrogen_atoms() > f2->hydrogen_atoms())
    return -1;

  if (f1->hydrogen_atoms() < f2->hydrogen_atoms())
    return 1;

// If still not resolved, resolve by fragment id

#ifdef DEBUG_CAREFUL_FRAG
  cerr << "Resolved only by fragment id\n";
#endif

  if (f1->frag_id() < f2->frag_id())
    return -1;

  return 1;
}

int
Molecule::_reduce_to_largest_fragment_carefully(Fragment_Data * fd,
                                                int * already_counted)
{
  int nf = number_fragments();

  for (int i = 0; i < nf; i++)
  {
    fd[i].set_frag_id(i);
  }

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    int f = fragment_membership[i];

    if (a->element()->organic())
      fd[f].another_organic_atom();
    else
      fd[f].another_non_organic_atom();
  }

  int highest_organic_atom_count = fd[0].organic_atoms();
  for (int i = 1; i < nf; i++)
  {
    const Fragment_Data & fdi = fd[i];

    if (fdi.organic_atoms() > highest_organic_atom_count)
      highest_organic_atom_count = fdi.organic_atoms();
  }

// Look for the features

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];     // not const because of implicit_hydrogens()

    int fid = fragment_membership[i];

    Fragment_Data & fdi = fd[fid];

    if (highest_organic_atom_count != fdi.organic_atoms())
      continue;

    if (7 == a->atomic_number())
    {
      if (_is_nitro(i, already_counted))
        fdi.another_nitro();
      else
      {
        fdi.another_nitrogen();
        already_counted[i] = 1;
        fdi.extra_hydrogens(a->implicit_hydrogens());
      }
    }
    else if (16 == a->atomic_number())
    {
      if (_is_sulphate_like(i, already_counted))
        fdi.another_sulphate();
      else
      {
        fdi.extra_hydrogens(a->implicit_hydrogens());
        fdi.another_sulphur();
      }
    }
    else if (15 == a->atomic_number())
      fdi.another_phosphorus();
    else if (1 == a->atomic_number())
      fdi.extra_hydrogens(1);
  }

// Now do oxygens, as some of them will have been consumed by Nitro and Sulphates

  for (int i = 0; i < _number_elements; i++)
  {
    if (already_counted[i])
      continue;

    const Atom * ai = _things[i];

    if (8 != ai->atomic_number())
      continue;

    int fid = fragment_membership[i];

    Fragment_Data & fdi = fd[fid];

    if (highest_organic_atom_count != fdi.organic_atoms())
      continue;

    fdi.another_oxygen();
  }

  qsort(fd, nf, sizeof(Fragment_Data), fragment_data_comparitor);

#ifdef DEBUG_CAREFUL_FRAG
  for (int i = 0; i < nf; i++)
  {
    fd[i].debug_print(cerr);
  }
#endif

  int fragment_to_keep = fd[0].frag_id();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    int f = fragment_membership[i];

    if (f != fragment_to_keep)
      atoms_to_be_removed.add(i);
  }

  return remove_atoms(atoms_to_be_removed);
}

/*
  Atom N is a nitrogen atom. Is it part of a Nitro group
*/

int
Molecule::_is_nitro(atom_number_t n, int * already_done) const
{
  const Atom * an = _things[n];

  assert (7 == an->atomic_number());

  if (3 != an->ncon())
    return 0;

  if (5 != an->nbonds())
    return 0;

  atom_number_t o1 = INVALID_ATOM_NUMBER;    // the first doubly bonded oxygen found

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = an->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(n);

    if (8 != _things[o]->atomic_number())
      continue;

    if (INVALID_ATOM_NUMBER == o1)
      o1 = o;
    else      // already got o1, so 'o' is the second doubly bonded oxygen
    {
      already_done[n] = 1;
      already_done[o1] = 1;
      already_done[o] = 1;

      return 1;
    }
  }

  return 0;
}

/*
  Atom S is a Sulphur. Is it part of one of the various SO combinations
*/

int
Molecule::_is_sulphate_like(atom_number_t s, int * already_done) const
{
  const Atom * as = _things[s];

  int scon = as->ncon();

  if (1 == scon)
    return 0;

  Set_of_Atoms oxygens;
  oxygens.resize(scon);

  int double_bonds_found = 0;

  for (int i = 0; i < scon; i++)
  {
    const Bond * b = as->item(i);

    if (b->is_double_bond())
      double_bonds_found++;
    else
      continue;

    atom_number_t o = b->other(s);

    if (8 == _things[o]->atomic_number())
      oxygens.add(o);
  }

  if (0 == double_bonds_found || oxygens.empty())
    return 0;

  already_done[s] = 1;

  oxygens.set_vector(already_done, 1);

  return 1;
}
