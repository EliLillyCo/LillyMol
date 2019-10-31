/*
  
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#define IW_MULTI_THREAD_OMP
#ifdef IW_MULTI_THREAD_OMP
#include <omp.h>
#endif

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "cmdline.h"
#include "report_progress.h"
#include "iw_stl_hash_map.h"
#include "iwbits.h"
#include "misc.h"

#define ISTREAM_AND_TYPE_IMPLEMENTATION
#include "istream_and_type.h"
#include "smiles.h"
#include "substructure.h"
#include "molecule_to_query.h"
#include "aromatic.h"
#include "target.h"
#include "path.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_joined = 0;

/*
  If we have a reaction with two join points, it might be forming an extra ring
*/

static int extra_rings_formed = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int remove_chirality = 1;

static int examine_ring_related_properties = 1;

static Report_Progress report_progress;

#define NPROPERTIES 22
#define MP_NATOMS 0
#define MP_CARBON 1
#define MP_NITROGEN 2
#define MP_OXYGEN 3
#define MP_FLUORINE 4
#define MP_SULPHUR 5
#define MP_CHLORINE 6
#define MP_BROMINE 7
#define MP_IODINE 8
#define MP_PHOSPHORUS 9
#define MP_HYDROGEN 10
#define MP_SINGLE_BONDS 11
#define MP_DOUBLE_BONDS 12
#define MP_TRIPLE_BONDS 13
#define MP_CSC 14
#define MP_CDC 15
#define MP_CSN 16
#define MP_CDN 17
#define MP_CSO 18
#define MP_CDO 19
#define MP_NRINGS 20
#define MP_OTHER_ELEMENT 21

static int atomic_number_to_index [] = {
  MP_OTHER_ELEMENT,    // 0
  MP_HYDROGEN,    // 1
  MP_OTHER_ELEMENT,    // 2
  MP_OTHER_ELEMENT,    // 3
  MP_OTHER_ELEMENT,    // 4
  MP_OTHER_ELEMENT,    // 5
  MP_CARBON,    // 6
  MP_NITROGEN,    // 7
  MP_OXYGEN,    // 8
  MP_FLUORINE,    // 9
  MP_OTHER_ELEMENT,    // 10
  MP_OTHER_ELEMENT,    // 11
  MP_OTHER_ELEMENT,    // 12
  MP_OTHER_ELEMENT,    // 13
  MP_OTHER_ELEMENT,    // 14
  MP_PHOSPHORUS,    // 15
  MP_SULPHUR,    // 16
  MP_CHLORINE,    // 17
  MP_OTHER_ELEMENT,    // 18
  MP_OTHER_ELEMENT,    // 19
  MP_OTHER_ELEMENT,    // 20
  MP_OTHER_ELEMENT,    // 21
  MP_OTHER_ELEMENT,    // 22
  MP_OTHER_ELEMENT,    // 23
  MP_OTHER_ELEMENT,    // 24
  MP_OTHER_ELEMENT,    // 25
  MP_OTHER_ELEMENT,    // 26
  MP_OTHER_ELEMENT,    // 27
  MP_OTHER_ELEMENT,    // 28
  MP_OTHER_ELEMENT,    // 29
  MP_OTHER_ELEMENT,    // 30
  MP_OTHER_ELEMENT,    // 31
  MP_OTHER_ELEMENT,    // 32
  MP_OTHER_ELEMENT,    // 33
  MP_OTHER_ELEMENT,    // 34
  MP_BROMINE,    // 35
  MP_OTHER_ELEMENT,    // 36
  MP_OTHER_ELEMENT,    // 37
  MP_OTHER_ELEMENT,    // 38
  MP_OTHER_ELEMENT,    // 39
  MP_OTHER_ELEMENT,    // 40
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 50
  MP_OTHER_ELEMENT,    // 51
  MP_OTHER_ELEMENT,    // 52
  MP_IODINE,    // 53
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT,    // 1
  MP_OTHER_ELEMENT     // 1
};

class Molecular_Properties
{
  public:
    int _mp[NPROPERTIES];

  private:
    void _increment_individual_bond (bond_type_t b, atomic_number_t z1, atomic_number_t z2);

  public:
    Molecular_Properties();

    template <typename T> void compute_molecular_properties (T & m);
    template <typename T> void increment_molecular_properties (T & m);
    template <typename T> void increment_bonds (const T & m1, const resizable_array<bond_type_t> & bonds, const T & m2) ;

    int  match (const Molecular_Properties & rhs) const;

    int natoms () const { return _mp[MP_NATOMS];}

    int is_subset (const Molecular_Properties & rhs) const;
};

Molecular_Properties::Molecular_Properties()
{
  std::fill_n(_mp, NPROPERTIES, 0);
}

void
Molecular_Properties::_increment_individual_bond (bond_type_t b,
                                                  atomic_number_t z1,
                                                  atomic_number_t z2)
{
  if (z1 > z2)
    std::swap(z1, z2);

  if (SINGLE_BOND == b)
  {
    _mp[MP_SINGLE_BONDS]++;
    if (6 == z1)
    {
      if (6 == z2)
        _mp[MP_CSC]++;
      else if (7 == z2)
        _mp[MP_CSN]++;
      else if (8 == z2)
        _mp[MP_CSO]++;
    }
  }
  else if (DOUBLE_BOND == b)
  {
    _mp[MP_DOUBLE_BONDS]++;
    if (6 == z1)
    {
      if (6 == z2)
        _mp[MP_CDC]++;
      else if (7 == z2)
        _mp[MP_CDN]++;
      else if (8 == z2)
        _mp[MP_CDO]++;
    }
  }
  else if (TRIPLE_BOND == b)
    _mp[MP_TRIPLE_BONDS]++;

  return;
}

template <typename T>
void
Molecular_Properties::compute_molecular_properties (T & m)
{
  std::fill_n(_mp, NPROPERTIES, 0);

  increment_molecular_properties(m);

  return;
}

template <typename T>
void
Molecular_Properties::increment_molecular_properties (T & m)
{
  const auto matoms = m.natoms();

  _mp[MP_NATOMS] += matoms;

  for (auto i = 0; i < matoms; ++i)
  {
    const auto z = m.atomic_number(i);
    _mp[atomic_number_to_index[z]]++;
  }

  const auto nb = m.nedges();

  for (auto i = 0; i < nb; ++i)
  {
    const auto b = m.bondi(i);

    const auto z1 = m.atomic_number(b->a1());
    const auto z2 = m.atomic_number(b->a2());

    if (b->is_single_bond())
      _increment_individual_bond(SINGLE_BOND, z1, z2);
    else if (b->is_double_bond())
      _increment_individual_bond(DOUBLE_BOND, z1, z2);
    else if (b->is_triple_bond())
      _increment_individual_bond(TRIPLE_BOND, z1, z2);
  }

  _mp[MP_NRINGS] += m.nrings() + extra_rings_formed;

  return;
}

template <typename T>
void
Molecular_Properties::increment_bonds (const T & m1,
                                       const resizable_array<bond_type_t> & bonds,
                                       const T & m2) 
{
  const auto nb = bonds.size();

  for (unsigned int i = 0; i < nb; ++i)
  {
    auto z1 = m1.atomic_number(i);
    auto z2 = m2.atomic_number(i);
//  cerr << "Bond type " << bonds[i] << " btw " << z1 << " and " << z2 << endl;

    _increment_individual_bond(bonds[i], z1, z2);
  }

  return;
}

int
Molecular_Properties::match (const Molecular_Properties & rhs) const
{
  for (auto i = 0; i < NPROPERTIES; ++i)
  {
//  cerr << "Molecular_Properties::match:i = " << i << " cmp " << _mp[i] << " and " << rhs._mp[i] << endl;

    if (_mp[i] != rhs._mp[i])
      return 0;
  }

  return 1;
}

int
Molecular_Properties::is_subset (const Molecular_Properties & rhs) const
{
  for (auto i = 0; i < NPROPERTIES; ++i)
  {
    if (_mp[i] > rhs._mp[i])
      return 0;
  }

  return 1;
}

#define NUMBER_AROMATIC_PROPERTIES 16
#define AMP_AROMATIC_ATOMS 0
#define AMP_AROMATIC_RINGS 1
#define AMP_ISOLATED_RINGS 2
#define AMP_IN_MORE_THAN_ONE_RING 3
#define AMP_AROMATIC_CARBON 4
#define AMP_AROMATIC_NITROGEN 5
#define AMP_AROMATIC_OXYGEN 6
#define AMP_RSIZE_3 7
#define AMP_RSIZE_4 8
#define AMP_RSIZE_5 9
#define AMP_RSIZE_6 10
#define AMP_RSIZE_7 11
#define AMP_RSIZE_X 12

class Aromatic_Properties
{
  private:
    int _mp[NUMBER_AROMATIC_PROPERTIES];

  public:

    template <typename T> void compute_ring_and_aromatic_properties (T &);
    template <typename T> void increment_ring_and_aromatic_properties (T & m);

    int matches (const Aromatic_Properties & rhs) const;

    int is_subset (const Aromatic_Properties & rhs) const;

    void debug_print (std::ostream &) const;
};

template <typename T>
void
Aromatic_Properties::compute_ring_and_aromatic_properties (T & m)
{
  std::fill_n(_mp, NUMBER_AROMATIC_PROPERTIES, 0);

  increment_ring_and_aromatic_properties(m);

  return;
}

template <typename T>
void
Aromatic_Properties::increment_ring_and_aromatic_properties (T & m)
{
  m.compute_aromaticity_if_needed();

  const auto nr = m.nrings();

  if (0 == nr)
    return;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (ri->is_aromatic())
      _mp[AMP_AROMATIC_RINGS]++;

    if (0 == ri->fused_ring_neighbours())
      _mp[AMP_ISOLATED_RINGS]++;

    const auto rsize = ri->size();

    if (rsize > 7)
      _mp[AMP_RSIZE_X]++;
    else
      _mp[AMP_RSIZE_3 + rsize - 3]++;
  }

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    const auto nri = m.nrings(i);

    if (0 == nri)
      continue;

    if (nri > 1)
      _mp[AMP_IN_MORE_THAN_ONE_RING]++;

    if (! m.is_aromatic(i))
      continue;

    const auto z = m.atomic_number(i);

    if (6 == z)
      _mp[AMP_AROMATIC_CARBON]++;
    else if (7 == z)
      _mp[AMP_AROMATIC_NITROGEN]++;
    else if (8 == z)
      _mp[AMP_AROMATIC_OXYGEN]++;
  }

  return;
}

int
Aromatic_Properties::matches (const Aromatic_Properties & rhs) const
{
  for (auto i = 0; i < NUMBER_AROMATIC_PROPERTIES; ++i)
  {
    if (_mp[i] != rhs._mp[i])
      return 0;
  }

  return 1;
}

int
Aromatic_Properties::is_subset (const Aromatic_Properties & rhs) const
{
  for (auto i = 0; i < NUMBER_AROMATIC_PROPERTIES; ++i)
  {
    if (_mp[i] > rhs._mp[i])
      return 0;
  }

  return 1;
}

void
Aromatic_Properties::debug_print (std::ostream & os) const
{
  os << "AP:";

  for (auto i = 0; i < NPROPERTIES; ++i)
  {
    os << ' ' << _mp[i];
  }
  os << '\n';

  return;
}

class Reagent : public Molecule
{
  private:
    Molecular_Properties _mp;
    Aromatic_Properties _amp;

//  We do a substructure search to figure out which target molecules we could possibly match

    IW_Bits_Base _can_match;
    int _nset;

  public:

    void compute_molecular_properties();

    int establish_possible_matches (Molecule_to_Match * target, int n);
    template <typename T> int establish_possible_matches_thread_safe (const resizable_array_p<T> & m,
                                                 Molecular_Properties * mpr, Aromatic_Properties * apr);

    int nset () const { return _nset;}
    int is_set (const int b) const { return _can_match.is_set(b);}

    int can_match (int s) const { return _can_match.is_set(s);}
    const IW_Bits_Base match_bits() const { return _can_match;}
};

void Reagent::compute_molecular_properties ()
{
  _mp.compute_molecular_properties(*this);

//cerr << Molecule::name() << " computing properties " << examine_ring_related_properties << endl;
  if (examine_ring_related_properties)
    _amp.compute_ring_and_aromatic_properties(*this);

  return;
}

int
Reagent::establish_possible_matches (Molecule_to_Match * target,
                                     int n)
{
  const auto matoms = Molecule::natoms();

  Substructure_Query q;

  Molecule_to_Query_Specifications mqs;

  q.create_from_molecule(*this, mqs);

  int nb = n / 32 * 32;
  if (nb != n)
    nb += 32;

  _can_match.allocate_space_for_bits(nb);

  int rc = 0;

  for (auto i = 0; i < n; ++i)
  {
    if (target[i].natoms() <= matoms)   // clearly this reagent cannot be used to create that target
      continue;

    if (q.substructure_search(target[i]))
    {
      _can_match.set(i, 1);
      rc++;
    }
  }

  _nset = _can_match.nset();

  return rc;
}

template <typename T>
int
Reagent::establish_possible_matches_thread_safe (const resizable_array_p<T> & m,
                                                 Molecular_Properties * mpr,
                                                 Aromatic_Properties * apr)
{
  const auto matoms = Molecule::natoms();

  Substructure_Query q;

  Molecule_to_Query_Specifications mqs;
//mqs.set_make_embedding(1);
  mqs.set_substituents_only_at_isotopic_atoms(1);

  q.create_from_molecule(*this, mqs);

  q.set_max_matches_to_find(1);
  q.set_save_matched_atoms(0);

//IWString foo("foo.qry");
//q.write_msi(foo);

  const int n = m.size();

  int nb = n / 32 * 32;
  if (nb != n)
    nb += 32;

  _can_match.allocate_space_for_bits(nb);

//int specialq = Molecule::name().contains("R(194)");
//if (specialq)
//  cerr << "Processing '" << Molecule::name() << "'\n";

  for (auto i = 0; i < n; ++i)
  {
    T & mi = *(m[i]);

    if (mi.natoms() <= matoms)   // clearly this reagent cannot be used to create that target
      continue;

    if (! _mp.is_subset(mpr[i]))
      continue;

    if (! _amp.is_subset(apr[i]))
      continue;

    Substructure_Results sresults;
    sresults.set_save_query_atoms_matched(0);

    if (q.substructure_search(mi, sresults))
      _can_match.set(i, 1);
//  else if (specialq)
//    cerr << "Zeo hits for '" << Molecule::name() << "'\n";
//  else
//    cerr << "Processing '" << mi.name() << " only matched " << sresults.max_query_atoms_matched_in_search() << endl;
  }

  _nset = _can_match.nset();

//if (specialq)
//  cerr << "Special molecule '" << Molecule::name() << " nset " << _nset << endl;

  return _nset;
}

class Molecule_to_Generate : public Molecule
{
  private:
    Molecular_Properties _mp;
    Aromatic_Properties _amp;
    IWString _unique_smiles;
    int _found;

  public:
    Molecule_to_Generate ();

    void compute_molecular_properties ();

    void set_found (int s) { _found = s;}
    int  found () const { return _found;}

    const Molecular_Properties & mpr () const { return _mp;}
    const Aromatic_Properties  & amp () const { return _amp;}

    const IWString & unique_smiles () const { return _unique_smiles;}

};

Molecule_to_Generate::Molecule_to_Generate ()
{
  _found = 0;

  return;
}

void
Molecule_to_Generate::compute_molecular_properties ()
{
  _unique_smiles = Molecule::unique_smiles();

  _mp.compute_molecular_properties(*this);

  if (examine_ring_related_properties)
    _amp.compute_ring_and_aromatic_properties(*this);

  return;
}

template <typename T>
class Set_of_Molecules_to_Generate
{
  private:
    resizable_array_p<T> _m;
    IW_STL_Hash_Map_int _unique_smiles;
    int * _atom_count_to_find;

  public:
    Set_of_Molecules_to_Generate();
    ~Set_of_Molecules_to_Generate();

    int read_targets (const char * fname, int input_type);
    int read_targets (data_source_and_type<T> & input);

    int initialise ();
};

template <typename T>
Set_of_Molecules_to_Generate<T>::Set_of_Molecules_to_Generate()
{
  _atom_count_to_find = nullptr;

  return;
}

template <typename T>
Set_of_Molecules_to_Generate<T>::~Set_of_Molecules_to_Generate ()
{
  if (nullptr != _atom_count_to_find)
    delete [] _atom_count_to_find;

  return;
}


template <typename T>
int
Set_of_Molecules_to_Generate<T>::read_targets (const char * fname,
                                        int input_type)
{
  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<T> input(input_type, fname);
  if (! input.good())
  {
    cerr << "read_reagents:cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return read_targets(input);
}

template <typename T>
int
Set_of_Molecules_to_Generate<T>::read_targets (data_source_and_type<T> & input)
{
  T * m;

  while (NULL != (m = input.next_molecule()))
  {
    preprocess(*m);

    _m.add(m);
  }

  return _m.size();
}

template <typename T>
int
Set_of_Molecules_to_Generate<T>::initialise ()
{
  const auto n = _m.size();

  for (auto i = 0; i < n; ++i)
  {
    const auto m = _m[i];

    m->compute_molecular_properties();
    if (examine_ring_related_properties)
      m->compute_ring_and_aromatic_properties();

    const IWString & s = m->unique_smiles();

    if (_unique_smiles.contains(s))
    {
      cerr << "Duplicate structure in target '" << m->name() << "'\n";
      return 0;
    }

    _unique_smiles[s] = i;
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Generates reaction products from labelled reagents\n";
  cerr << "  -R <fname>    files of isotopically labelled reagents (sorted by atom count)\n";
  cerr << "  -b <1,2,3>    bond types to use\n";
  cerr << "  -x <n>        the reaction forms <n> extra rings compared to reagents\n";
  cerr << "  -X <ele>      the reagents have extra <ele> atoms instead of isotopes\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  cerr << "  -r <n>        report progress every <n> products\n";
  cerr << "  -M <stem>     write reagents to files stem<filecount>\n";

  exit(rc);
}

/*
  We only process the first atom we find, don't bother checking for other atoms
*/

template <typename T>
int
replace_atoms_with_isotopes (T & m, 
                             const Element * e)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (e != m.elementi(i))
      continue;

    if (1 != m.ncon(i))
    {
      cerr << "Element " << e->symbol() << " not singly connected, cannot process\n";
      return 0;
    }

    auto j = m.other(i, 0);

    m.remove_atom(i);
    if (j > i)
      j--;

    if (j > 0)
      m.swap_atoms(j, 0);

    m.set_isotope(0, 1);

    return 1;
  }

  return 0;
}

template <typename T>
int
replace_atoms_with_isotopes(resizable_array_p<T> & m, 
                            const Element * e)
{
  const auto n = m.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    T & mi = *(m[i]);

    if (! replace_atoms_with_isotopes(*(m[i]), e))
    {
      cerr << "No atoms of type '" << e->symbol() << "' in '" << m[i]->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T>
void
remove_isotopes (resizable_array_p<T> & m, const int niso)
{
  const auto n = m.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    T & mi = *(m[i]);

    auto matoms = mi.natoms();
    if (matoms > niso)
      matoms = niso;

    for (auto j = 0; j < matoms; ++j)
    {
      mi.set_isotope(j, 0);
    }
  }

  return;
}

template <typename T>
int
write_set_of_molecules(const resizable_array_p<T> & m, 
                       IWString_and_File_Descriptor & output)
{
  const auto n = m.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    output << m[i]->smiles() << ' ' << m[i]->name() << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  output.flush();

  return output.good();
}

template <typename T>
int
write_set_of_molecules(const resizable_array_p<T> & m, IWString & fname)
{
  IWString_and_File_Descriptor output;

  if (! output.open(fname.null_terminated_chars()))
  {
    cerr << "write_set_of_molecules::cannot open '" << fname << "'\n";
    return 0;
  }

  return write_set_of_molecules(m, output);
}

void
assemble_non_zero_values (const extending_resizable_array<int> & v,
                          resizable_array<int> & nzero)
{
  const auto n = v.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    if (0 != v[i])
      nzero.add(i);
  }

  return;
}

int
are_the_same (const resizable_array<int> & a1,
              const resizable_array<int> & a2)
{
  return a1 == a2;   

#ifdef XUSE_OPERATOR
  const auto n = a1.size();

  if (n != a2.size())
    return 0; 

  for (unsigned int i = 0; i < n; ++i)
  {
    if (a1[i] != a2[i])
      return 0;
  }

  return 1;
#endif
}

template <typename T>
int
molecules_sorted_by_atom_count (const resizable_array_p<T> & mols)
{
  const auto n = mols.size();

  auto prev_natoms = mols[0]->natoms();

  for (unsigned int i = 1; i < n; ++i)
  {
    const auto a = mols[i]->natoms();

    if (a < prev_natoms)
    {
      cerr << "Atom count out of order, i = " << i << " molecule '" << mols[i]->name() << "'\n";
      return 0;
    }

    prev_natoms = a;
  }

  return 1;
}

template <typename T>
int
identify_isotopes_molecule (T & m,
                            extending_resizable_array<int> & isotopes)
{
  const auto matoms = m.natoms();

  Set_of_Atoms iso;

  for (auto i = 0; i < matoms; ++i)
  {
    const auto j = m.isotope(i);
    if (0 == j)
      continue;

    isotopes[j]++;
    iso.add(i);
  }

  const auto rc = iso.size();

  if (0 == rc)    // no isotopes
    return 0;

  if (1 == rc)
  {
    const auto j = iso[0];
    if (0 != j)
      m.swap_atoms(j, 0);
    m.unset_all_implicit_hydrogen_information(0);
    return 1;
  }

// We need to get the atoms ordered with the isotopic atoms first, and in order of increasing isotope
// First just move them to the front of the list

  for (unsigned int i = 0; i < rc; ++i)
  {
    const auto j = iso[i];
    if (i != j)
      m.swap_atoms(i, j);

    m.unset_all_implicit_hydrogen_information(i);
  }

  if (2 == rc)
  {
    if (m.isotope(0) > m.isotope(1))
      m.swap_atoms(0, 1);
    return 2;
  }

// We don't handle more than two atoms, but let's see if we are just lucky

  for (unsigned int i = 1; i < rc; ++i)
  {
    if (m.isotope(i-1) < m.isotope(i))
      continue;

    cerr << "More than two isotopes not implemented, contact Ian\n";
    return 0;
  }

  return rc;
}

template <typename T>
int
identify_isotopes (const resizable_array_p<T> & mols,
                   extending_resizable_array<int> & isotopes)
{
  const auto n = mols.size();

  int rc = 0;

  for (unsigned int i = 0; i < n; ++i)
  {
    const auto niso = identify_isotopes_molecule(*(mols[i]), isotopes);

    if (0 == niso)
      return 0;

    if (0 == rc)
      rc = niso;
    else if (niso != rc)
    {
      cerr << "Found " << niso << " isotopic atoms in '" << mols[i]->name() << " but found " << rc << " in other molecules in the file\n";
      return 0;
    }
  }

  return rc;
}

template <typename T>
void
preprocess (T & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (remove_chirality)
    m.remove_all_chiral_centres();

  return;
}


template <typename T>
int
read_reagents (data_source_and_type<T> & input,
               resizable_array_p<T> & mols)
{
  T * m;

  while (NULL != (m = input.next_molecule()))
  {
    preprocess(*m);

    mols.add(m);
  }

  return mols.size();
}

template <typename T>
int
read_reagents (const char * fname,
               int input_type,
               resizable_array_p<T> & mols)
{
  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<T> input(input_type, fname);
  if (! input.good())
  {
    cerr << "read_reagents:cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return read_reagents(input, mols);
}

static int
determine_highest_atom_count_to_find (const extending_resizable_array<int> & f)
{
  for (auto i = f.size() - 1; i >= 0; --i)
  {
    if (0 != f[i])
      return i;
  }

  cerr << "Huh, not looking for any particular atom counts!!!????\n";

  return 0;  
}

static int
compute_nfound (const resizable_array_p<Molecule_to_Generate> & m)
{
  const auto n = m.size();

  int rc = 0;

  for (unsigned int i = 0; i < n; ++i)
  {
    if (m[i]->found() > 0)
      rc++;
  }

  return rc;
}

static void
do_join (Molecule & m1,
         const Molecule & m2,
         const resizable_array<bond_type_t> & bonds)
{
  const auto initial_matoms = m1.natoms();

  m1.add_molecule(&m2);

  for (unsigned int i = 0; i < bonds.size(); ++i)
  {
    m1.add_bond(0, initial_matoms + i, bonds[i]);
  }

  return;
}

/*
  We return the index of the target molecule we find
*/

static int
find_a_matching_target_molecule (Molecule & m1,
                                 const IW_Bits_Base & m1_match_bits,
                                 Reagent & m2,
                                 const resizable_array<bond_type_t> & bonds,
                                 const resizable_array_p<Molecule_to_Generate> & to_generate,
                                 const IW_STL_Hash_Map_int & usmi_to_ndx)

{
  Molecular_Properties mpr;
  mpr.compute_molecular_properties(m1);
  mpr.increment_molecular_properties(m2);

  const Molecule & m2r = m2;
  mpr.increment_bonds(m1, bonds, m2r);

  Aromatic_Properties amp;
  if (examine_ring_related_properties)
  {
    amp.compute_ring_and_aromatic_properties(m1);
    amp.increment_ring_and_aromatic_properties(m2);
  }

//cerr << "From '" << m1.smiles() << " and " << m2.smiles() << " what do we get, nb " << bonds.size() << endl;

  const auto n = to_generate.size();

  int molecule_joined = 0;

  for (unsigned int i = 0; i < n; ++i)
  {
    if (! m1_match_bits.is_set(i) || ! m2.can_match(i))
      continue;

    const Molecule_to_Generate * t = to_generate[i];

    if (t->found())
      continue;

    const auto tmatoms = t->natoms();

    if (tmatoms < mpr.natoms())
      continue;

    if (tmatoms > mpr.natoms())    // done because to_generate is sorted
      return -1;

    if (! mpr.match(t->mpr()))
      continue;

    if (! molecule_joined)
    {
      do_join(m1, m2, bonds);
      molecule_joined = 1;
      molecules_joined++;
//    cerr << "Generated " << m1.smiles() << "\n";
    }

    if (examine_ring_related_properties && ! amp.matches(t->amp()))
      continue;

//  cerr << "Compare " << m1.unique_smiles() << " and " << t->unique_smiles() << endl;
//  cerr << "Found match for taret molecule " << i << endl;
    if (m1.unique_smiles() == t->unique_smiles())
    {
      to_generate[i]->set_found(1);

      return i;
    }

//  now that we have generated the unique smiles, check it against all

    const auto f = usmi_to_ndx.find(m1.unique_smiles());
    if (f == usmi_to_ndx.end())
      continue;

    const auto j = (*f).second;
    to_generate[j]->set_found(1);
    return j;
  }

  return -1;    // no molecule found
}

/*
*/


static int
molecules_from_reactions (resizable_array_p<Molecule_to_Generate> & to_generate,
                          const extending_resizable_array<int> & atom_count_to_find,
                          const int nr,
                          resizable_array_p<Reagent> * reagent,
                          const resizable_array<bond_type_t> & bonds,
                          IWString_and_File_Descriptor & output)
{
  assert (2 == nr);

  resizable_array_p<Reagent> & r1 = reagent[0];
  const auto nr1 = r1.size();
  resizable_array_p<Reagent> & r2 = reagent[1];
  const auto nr2 = r2.size();

  const auto highest_atom_count_to_find = determine_highest_atom_count_to_find(atom_count_to_find);

  if (verbose)
    cerr << "Highest atom count requested " << highest_atom_count_to_find << endl;

  IW_STL_Hash_Map_int usmi_to_ndx;

  for (unsigned int i = 0; i < to_generate.size(); ++i)
  {
    usmi_to_ndx[to_generate[i]->unique_smiles()] = i;
  }

  if (verbose)
    cerr << "Unique smiles hash filled\n";

  unsigned int nfound = 0;

  for (unsigned int i = 0; i < nr1; ++i)
  {
    Molecule m(*(r1[i]));

    const auto initial_matoms = m.natoms();

    if (initial_matoms >= highest_atom_count_to_find)
      break;

    for (unsigned int j = 0; j < nr2; ++j)
    {
      const Molecule * mj = r2[j];

//    cerr << "Product of " << m.name() << " and " << mj->name() << " will contain " << (initial_matoms + mj->natoms()) << " atoms\n";

      if (0 == atom_count_to_find[initial_matoms + mj->natoms()])
        continue;

      if (initial_matoms + mj->natoms() > highest_atom_count_to_find)
        break;

//    cerr << "Passed atom count test " << m.smiles() << " " << r2[j]->smiles() << endl;

      const auto f = find_a_matching_target_molecule(m, r1[i]->match_bits(), *(r2[j]), bonds, to_generate, usmi_to_ndx);
      if (f >= 0)
      {
        output << to_generate[f]->smiles() << ' ' << to_generate[f]->name();
        output << ' ' << r1[i]->smiles() << ' ' << r1[i]->name();
        output << ' ' << r2[j]->smiles() << ' ' << r2[j]->name();
        output << '\n';
        output.write_if_buffer_holds_more_than(4096);

        nfound++;

        if (nfound == to_generate.size())
          return nfound;
      }

      if (m.natoms() > initial_matoms)
        m.resize(initial_matoms);

      if (report_progress())
        cerr << "Examined " << report_progress.times_called() << " molecules, found " << compute_nfound(to_generate) << " of " << to_generate.size() << ", joined " << molecules_joined << " molecules, completed " << i << " (" << initial_matoms << " atoms) of " << nr1 << " reagents\n";
    }
  }

  return nfound;
}

static int
any_reagents_can_make (const resizable_array_p<Reagent> & reagents,
                       const int b)
{
  const auto n = reagents.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    if (reagents[i]->is_set(b))
      return 1;
  }

  return 0;
}

static int
molecules_from_reactions (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lR:b:r:x:M:X:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  set_include_isotopic_information_in_unique_smiles(0);

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress reporting mechanism (-r)\n";
      usage(2);
    }
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', extra_rings_formed) || extra_rings_formed < 0)
    {
      cerr << "The number of extra rings formed (-x) must be a non negative number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Reaction results in " << extra_rings_formed << " extra rings compared to reagents\n";
  }

  int nr =cl.option_count('R');

  if (0 == nr)
  {
    cerr << "Must specify at least two reagent files via the -R option\n";
    usage(1);
  }

  if (1 == nr)
  {
    cerr << "Must specify at least two reagent files via the -R option\n";
    usage(1);
  }

  resizable_array_p<Reagent> * reagents = new resizable_array_p<Reagent>[nr]; std::unique_ptr<resizable_array_p<Reagent>[]> free_reagents(reagents);

  unsigned long long combinations = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const char * r = cl.option_value('R', i);

    if (! read_reagents(r, input_type, reagents[i]))
    {
      cerr << "Cannot read reagents from '" << r << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << reagents[i].size() << " reagents from '" << r << "'\n";

    if (! molecules_sorted_by_atom_count(reagents[i]))
    {
      cerr << "Reagent sets must be sorted by atom count, use msort. '" << r << "'\n";
      return 2;
    }

    if (0 == combinations)
      combinations = reagents[i].size();
    else
      combinations *= reagents[i].size();
  }

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    const Element * e = get_element_from_symbol_no_case_conversion(x);
    if (NULL == e)
    {
      cerr << "NO element of type '" << x << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Will remove and label atoms of type '" << x << "'\n";

    for (auto i = 0; i < nr; ++i)
    {
      replace_atoms_with_isotopes(reagents[i], e);
    }
  }

  for (auto i = 0; i < nr; ++i)
  {
    const auto & ri = reagents[i];

    const auto n = ri.size();

    for (unsigned int j = 0; j < n; ++j)
    {
      ri[j]->compute_molecular_properties();
    }
  }

  extending_resizable_array<int> * isotopes_used = new extending_resizable_array<int>[nr]; std::unique_ptr<extending_resizable_array<int>[]> free_isotopes_used(isotopes_used);

  if (verbose)
    cerr << "IDentifying isotopes in " << nr << " sets of reagents\n";

  for (auto i = 0; i < nr; ++i)
  {
    if (! identify_isotopes(reagents[i], isotopes_used[i]))
    {
      cerr << "Reagent set " << i << " cannot identify isotopes used\n";
      return 1;
    }
  }

  int niso = 0;

  for (auto i = 0; i < nr; ++i)
  {
    resizable_array<int> i1;
    assemble_non_zero_values(isotopes_used[i], i1);

    for (auto j = i + 1; j < nr; ++j)
    {
      resizable_array<int> i2;
      assemble_non_zero_values(isotopes_used[j], i2);

      if (! are_the_same(i1, i2))
      {
        cerr << "Isotope count mismatch between reagent set " << i << " (" << i1.size() << ") and " << j << " (" << i2.size() << ")\n";
        return 2;
      }
    }

    niso = i1.size();
    if (verbose > 1)
    {
      cerr << "Found isotopes ";
      for (auto i = 0; i < niso; ++i)
      {  
        cerr << ' ' << i1[i];
      }  
      cerr << endl;
    }
  }

  if (verbose)
    cerr << nr << " reagent sets, using " << niso << " isotopes, can generate " << combinations << " product molecules\n";

  if (2 != nr)
  {
    cerr << "Sorry, don't know how to process " << nr << " reagent sets, see Ian\n";
    return 1;
  }

  resizable_array<bond_type_t> bonds_to_add;

  if (cl.option_present('b'))
  {
    if (niso != cl.option_count('b'))
    {
      cerr << "Reagents each have " << niso << " isotopes, but " << cl.option_count('b') << " bond making specifications, impossible\n";
      return 1;
    }

    int i = 0;
    int b = 0;
    while (cl.value('b', b, i++))
    {
      if (1 == b)
        bonds_to_add.add(SINGLE_BOND);
      else if (2 == b)
        bonds_to_add.add(DOUBLE_BOND);
      else if (3 == b)
        bonds_to_add.add(TRIPLE_BOND);
      else
      {
        cerr << "Unrecognised bond type " << b << "\n";
        usage(1);
      }
    }
  }
  else
  {
    for (auto i = 0; i < niso; ++i)
    {
      bonds_to_add.add(SINGLE_BOND);
    }
  }

  if (1 != niso)
  {
    examine_ring_related_properties = 0;

    cerr << "More than two attachment points, will not consider ring properties\n";
  }

//set_substituents_only_at_isotopic_atoms(1);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  resizable_array_p<Molecule_to_Generate> to_generate;

  for (auto i = 0; i < cl.number_elements(); ++i)
  {
    if (! read_reagents(cl[i], input_type, to_generate))
    {
      cerr << "Cannot read molecules from '" << cl[i] << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << to_generate.size() << " molecules to generate\n";
  }

  if (! molecules_sorted_by_atom_count(to_generate))
  {
    cerr << "The set of molecules to be generated must be sorted by atom count, use msort\n";
    return 2;
  }

  const auto ngenerate = to_generate.size();

  extending_resizable_array<int> atom_count_to_find;

  for (unsigned int i = 0; i < ngenerate; ++i)
  {
    to_generate[i]->compute_molecular_properties();

    atom_count_to_find[to_generate[i]->natoms()]++;
  }

  Molecular_Properties * mp = new Molecular_Properties[ngenerate]; std::unique_ptr<Molecular_Properties[]> free_mp(mp);
  Aromatic_Properties *  ap = new Aromatic_Properties[ngenerate];  std::unique_ptr<Aromatic_Properties[]>  free_ap(ap);
  for (unsigned int i = 0; i < ngenerate; ++i)
  {
    auto ti = to_generate[i];
    mp[i].compute_molecular_properties(*ti);
    ap[i].compute_ring_and_aromatic_properties(*ti);
  }

  if (verbose)
    cerr << "Computed molecular properties for " << ngenerate << " target molecules\n";

#if defined(IW_MULTI_THREAD_OMP)
#pragma omp parallel
  for (auto i = 0; i < nr; i++)
  {
    auto & ri = reagents[i];     // a set of reagents

    const auto n = ri.size();

#pragma omp for schedule(dynamic,512) nowait
    for (unsigned int j = 0; j < n; ++j)
    {
      ri[j]->establish_possible_matches_thread_safe(to_generate, mp, ap);
    }
  }
#else     // single threaded, can share Molecule_to_Match objects 

  if (1)
  {
    Molecule_to_Match * target = new Molecule_to_Match[ngenerate]; std::unique_ptr<Molecule_to_Match[]> free_target(target);
    for (unsigned int i = 0; i < ngenerate; ++i)
    {
      target[i].initialise_molecule(to_generate[i]);
    }
    for (auto i = 0; i < nr; ++i)
    {
      auto & ri = reagents[i];

      const auto n = ri.size();

      for (auto j = 0; j < n; ++j)
      {
        ri[j]->establish_possible_matches(target, static_cast<int>(ngenerate));
      }
    }
  }
#endif

  for (auto i = 0; i < nr; ++i)
  {
    auto & ri = reagents[i];

    const auto initial_size = ri.size();

    for (int j = initial_size - 1; j >= 0; j--)
    {
//    cerr << "Reagent set " << i << " reagent " << j << " '" << ri[j]->name() << "' nset " << ri[j]->nset() << endl;
      if (0 == ri[j]->nset())
        ri.remove_item(j);
    }

    const auto final_size = ri.size();

    if (verbose)
      cerr << "Reagent set " << i << " started with " << initial_size << " reagents, " << final_size << " could match a product\n";

    if (0 == final_size)
    {
      cerr << "Reagent set " << i << " does not match any of the proposed product molecules\n";
      return 2;
    }

    if (cl.option_present('M'))
    {
      IWString m = cl.string_value('M');
      m << i << ".smi";
      write_set_of_molecules(ri, m);
    }
  }

// At this stage, we no longer need the isotopic information

  for (auto i = 0; i < nr; ++i)
  {
    remove_isotopes(reagents[i], niso);
  }

// Any product molecules that cannot be made by anything should be dropped

  if (1)            // always
  {
    const auto initial_size = to_generate.size();

    int molecules_dropped = 0;
    for (int i = initial_size - 1; i >= 0; --i)
    {
      int can_be_made = 1;
      for (auto j = 0; j < nr; ++j)
      {
        if (0 == any_reagents_can_make(reagents[j], i))
        {
          can_be_made = 0;
          break;
        }
      }

      if (! can_be_made)
      {
        to_generate[i]->set_found(-1);    // special flag to indicate not attempted
        molecules_dropped++;
      }
    }

    if (verbose)
      cerr << "Dropped " << molecules_dropped << " of " << initial_size << " molecules due to lack of suitable reagents\n";
  }

  IWString_and_File_Descriptor output(1);

  molecules_from_reactions(to_generate, atom_count_to_find, nr, reagents, bonds_to_add, output);

  output.flush();

  if (verbose)
  {
    cerr << "Processed " << to_generate.size() << " molecules, found " << compute_nfound(to_generate) << endl;
    cerr << molecules_joined << " actual product molecules formed\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = molecules_from_reactions(argc, argv);

  return rc;
}
