#include <iostream>
#include <memory>
#include <limits>
using std::cerr;
using std::endl;


#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/md5.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/linear_path_fingerprint.h"

namespace LFP {

#ifdef COMPILE_WITH_WATCH_BITS

static std::unordered_set<unsigned int> watch_bit;
static const Molecule * current_molecule = nullptr;

//cerr << "formed " << b__.bit_number(q__) << ' ' << reason__ << endl;  // if we want to monitor bit formation

  #define BIT_FORMED(q__, b__, reason__) {\
      const auto f = watch_bit.find(b__.bit_number(q__)); \
      if (f != watch_bit.end()) \
        write_state(*current_molecule, b__.bit_number(q__), reason__, cerr);\
  }\

void
add_bit_to_watch(const unsigned int b)
{
  watch_bit.insert(b);
}

static int
read_bits_to_watch (iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    unsigned int b;
    if (! buffer.numeric_value(b))
    {
      cerr << "read_bits_to_watch:invalid bit " << buffer << "'\n";
      return 0;
    }

    watch_bit.insert(b);
  }

  return watch_bit.size();
}

int
read_bits_to_watch (const IWString & fname)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "LFP::read_bits_to_watch:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_bits_to_watch(input);
}

#else
//#define BIT_FORMED(q__, b__, r__) {cerr << "formed " << b__.bit_number(q__) << ' ' << r__ << endl;}
  #define BIT_FORMED(q__, b__, r__)
#endif

static int default_min_path_length = 0;

void
set_default_min_path_length (const int s)
{
  default_min_path_length = s;
}
static int default_max_path_length = 7;

void
set_default_max_path_length (const int s)
{
  default_max_path_length = s;
}

static int default_min_heteroatoms_at_path_ends = -1;
void
set_default_min_heteroatoms_at_path_ends (const int s)
{
  default_min_heteroatoms_at_path_ends = s;
}


/*
  There are some bits set that will not permit accurate perception of
  substructures. These can be enabled or disabled as needed.
*/

static int default_only_bits_preserving_substructure_perception = 0;

static int defalt_min_heteroatoms_at_path_ends = std::numeric_limits<int>::max();

void 
set_defalt_min_heteroatoms_at_path_ends (const int s)
{
  defalt_min_heteroatoms_at_path_ends = s;
}

static int default_omit_ch2 = 0;

void set_default_omit_ch2 (const int s)
{
  default_omit_ch2 = s;
}

static int default_set_formal_charge_bits = 0;

void
set_default_set_formal_charge_bits (const int s)
{
  default_set_formal_charge_bits = s;
}

static int default_formally_charged_atoms_lose_identity = 0;

void
set_default_formally_charged_atoms_lose_identity(const int s)
{
  default_formally_charged_atoms_lose_identity = s;
}

static int default_max_path_length_isotopic_bits = -1;

void
set_default_max_path_length_isotopic_bits (const int s)
{
  default_max_path_length_isotopic_bits = s;
}

static int default_place_bits_for_rings = 0;

void
set_default_place_bits_for_rings (const int s)
{
  default_place_bits_for_rings = s;
}

static int default_bits_for_atoms_in_multiple_rings = 0;

void
set_default_bits_for_atoms_in_multiple_rings (const int s)
{
  default_bits_for_atoms_in_multiple_rings = s;
}

static int default_terminal_path_bits = -1;

void
set_default_terminal_path_bits (const int s)
{
  default_terminal_path_bits = s;
}

static int default_maximal_daylight_compatibility = 0;

void
set_maximal_daylight_compatibility (const int s)
{
  default_maximal_daylight_compatibility = s;
}

static int default_bits_for_hydrogen_attachments = -1;

void
set_default_bits_for_hydrogen_attachments(int s)
{
  default_bits_for_hydrogen_attachments = s;
}

static int default_differentiate_ring_and_chain_bonds = 0;

void
set_default_differentiate_ring_and_chain_bonds(const int s)
{
  default_differentiate_ring_and_chain_bonds = s;
}

Linear_Fingerprint_Defaults::Linear_Fingerprint_Defaults()
{
  _min_path_length = default_min_path_length;
  _max_path_length = default_max_path_length;
  _min_heteroatoms_at_path_ends = -1;
  _only_bits_preserving_substructure_perception = default_only_bits_preserving_substructure_perception;
  _omit_ch2 = default_omit_ch2;
  _set_formal_charge_bits = default_set_formal_charge_bits;
  _formally_charged_atoms_lose_identity = default_formally_charged_atoms_lose_identity;
  _max_path_length_isotopic_bits = default_max_path_length_isotopic_bits;
  _isotopic_paths = 0;
  _place_bits_for_rings = default_place_bits_for_rings;
  _bits_for_atoms_in_multiple_rings = default_bits_for_atoms_in_multiple_rings;
  _terminal_path_bits = default_terminal_path_bits;
  _maximal_daylight_compatibility = default_maximal_daylight_compatibility;
  _formally_charged_atoms_lose_identity = default_formally_charged_atoms_lose_identity;
  _bits_for_hydrogen_attachments = default_bits_for_hydrogen_attachments;
  _differentiate_ring_and_chain_bonds = default_differentiate_ring_and_chain_bonds;

  return;
}

static Linear_Fingerprint_Defaults default_linear_fingerprint_defaults;

Linear_Fingerprint_Creator::Linear_Fingerprint_Creator()
{
  _lfpd = &default_linear_fingerprint_defaults;
}

Linear_Fingerprint_Creator::Linear_Fingerprint_Creator (const Linear_Fingerprint_Defaults * s)
{
  _lfpd = s;
}

class MFingerprint
{
  private:
    int _matoms;

    atom_number_t * _path;     // which atoms are in the path
    int * _bond;               // the bonds along the path
    int * _in_path;            // whether or not an atom is in the path

    Atom ** _atoms;            // the atoms of the molecule

    aromaticity_type_t * _arom;    // atom aromaticity

    atomic_number_t * _atomic_number;    // atomic number of each atom in the molecule

    int * _ncon;               // ncon() for each atom in the molecule

    int * _unsaturation;       // unsaturation value for each atom

//  resizable_array_p<List_of_Ring_Sizes> _lors;   // ring sizes for each atom

    int * _atom_hash_value;     // hash value for each atom

    int * _path_hash_value;     // atom hash values along the path

    int * _ch2;                 // those path members which are CH2's

    int * _nrings;

    int _molecule_contains_fused_rings;

    int * _hcount;

    isotope_t * _isotope;

    int _path_length;     // the length (in bonds) of the current path

//  We can specify a minimum number of heteroatoms at the ends of a path.
//  Typically this is used for hit followup work where we want to minimise
//  the influence of carbon atoms.
//  Valid values are 0 (no restrictions) 1 or 2

    int _min_heteroatoms_at_ends;

#ifdef COUNT_TIMES_ATOM_IN_PATH
    int * _times_in_path;
#endif

//  private functions

    int _compute_atom_hash_value (atom_number_t i) const;
    int _compute_atom_hash_value_non_periodic_table (atom_number_t, int) const;
    int _compute_atom_hash_value_include_unsaturation (atom_number_t) const;
    template <typename T> void _set_cluster_bit (atom_number_t centre,
                           int b1, atom_number_t a1, int b2, atom_number_t a2, int b3, atom_number_t a3, T & bits);
    template <typename T> void _build_branched_paths (const Linear_Fingerprint_Defaults *, const int * include_these_atoms, T & bits);
    template <typename T> int  _build (atom_number_t aprev, const int * include_these_atoms, const Linear_Fingerprint_Defaults * lfpd, T & bits);
    template <typename T> void _set_auxiliary_bits (unsigned int zbit1, unsigned int zbit2, int is_ring, int astart, int astop, const Linear_Fingerprint_Defaults *, T & bits) const;
    template <typename T> void _set_bit_forward (int, const Linear_Fingerprint_Defaults * lfpd, T & bits) const;
    template <typename T> void _set_bit_backward (int, const Linear_Fingerprint_Defaults * lfpd, T & bits) const;
    template <typename T> void _set_bit (int, const Linear_Fingerprint_Defaults * lfpd, T & bits) const;
    template <typename T> void _omit_ch2 (int, T & bits) const;
    template <typename T> void _omit_ch2 (int, int, T & bits) const;
    template <typename T> void _set_bit_omit (int, int, T & bits) const;
    template <typename T> void _do_formal_charge_bits (int, int, int, T & bits) const;
    template <typename T> void _do_formal_charge_bits (int, int, T & bits) const;
    template <typename T> void _do_isotope_bits (int istart, int istop, int istep, T & bits) const;
    template <typename T> void _do_isotope_bits (T & bits) const;
    template <typename T> void _do_ring_bits (const Linear_Fingerprint_Defaults * lfpd, T & bits) const;
    template <typename T> void _set_terminal_path_bits (const int direction, T & bits) const;
    template <typename T> void _set_terminal_path_bits (T & bits) const;
    template <typename T> void _set_bits_for_atoms_in_multiple_rings(int direction, const Linear_Fingerprint_Defaults * lfpd, T & bits) const;

    int write_state (const Molecule & m, const unsigned int b, const char * reason, std::ostream &) const;

  public:
    MFingerprint (Molecule &, const int * atype, const Linear_Fingerprint_Defaults *);
    ~MFingerprint();

    void set_min_heteroatoms_at_path_ends (int);

    template <typename T> int fingerprint_large_ring (const Ring &, const Linear_Fingerprint_Defaults * lfpd, T & bits);

    template <typename T> int build (atom_number_t, const int * inc, const Linear_Fingerprint_Defaults *, T & bits);

#ifdef COUNT_TIMES_ATOM_IN_PATH
    void place_labels_for_times_in_path(Molecule & m) const;
#endif
};

/*
*/

MFingerprint::MFingerprint (Molecule & m,
                            const int * atype,
                            const Linear_Fingerprint_Defaults * lfpd)
{
  const int max_path_length = lfpd->max_path_length();

  _matoms = m.natoms();

  _path = new atom_number_t[_matoms + 1];    // don't use max path length because we may be fingerprinting large rings
  _bond = new int[max_path_length];

  _in_path = new int[_matoms];

  _atom_hash_value = new int[_matoms];

  _path_hash_value = new int[max_path_length + 1];    // atom hash values along the path

  _ch2 = new int[max_path_length + 1];

  _ncon = new int[_matoms];

  m.ncon(_ncon);

  _atomic_number = new atomic_number_t[_matoms];

  m.atomic_numbers(_atomic_number);

  _arom = new aromaticity_type_t[_matoms];

  m.aromaticity(_arom);

  _unsaturation = new int[_matoms];

  _atoms = new Atom *[_matoms];
  m.atoms((const Atom **) _atoms);    // the a->nbonds() method is non const, that's why we need a cast

  _hcount = new int[_matoms];

  for (int i = 0; i < _matoms; ++i)
  {
    _hcount[i] = m.hcount(i);
  }

  if (lfpd->max_path_length_isotopic_bits() >= 0)
  {
    _isotope = new isotope_t[_matoms];
    m.get_isotopes(_isotope);
  }
  else
    _isotope = nullptr;

  _nrings = new int[_matoms];
  m.ring_membership(_nrings);

  copy_vector(_atom_hash_value, atype, _matoms);
  set_vector(_unsaturation, _matoms, 0);

  _molecule_contains_fused_rings = 0;

  if (lfpd->bits_for_atoms_in_multiple_rings())
  {
    for (int i = 0; i < _matoms; i++)
    {
      if (_nrings[i] > 1)
      {
        _molecule_contains_fused_rings = 1;
        break;
      }
    }
  }

  _path_length = -1;

#ifdef CHECK_COLLISIONS
  _first_path = new resizable_array<int> * [bits_per_iwmfingerprint];
  for (int i = 0; i < bits_per_iwmfingerprint; i++)
  {
    _first_path[i] = nullptr;
  }
#endif

  _min_heteroatoms_at_ends = 0;

#ifdef COUNT_TIMES_ATOM_IN_PATH
  _times_in_path = new_int(_matoms);
#endif

  return;
}

MFingerprint::~MFingerprint()
{
  if (-9 == _matoms)
    cerr << "Deleting already deleted MFingerprint\n";

  _matoms = -9;

  delete [] _atoms;
  delete [] _path;
  delete [] _bond;
  delete [] _in_path;
  if (nullptr != _arom)
    delete [] _arom;
  delete [] _ncon;
  if (nullptr != _atomic_number)
    delete [] _atomic_number;

  delete [] _unsaturation;
  delete [] _nrings;

  delete [] _atom_hash_value;
  delete [] _path_hash_value;
  delete [] _ch2;

  if (nullptr != _hcount)
    delete [] _hcount;

  if (nullptr != _isotope)
    delete [] _isotope;

  _path_length = -4;

#ifdef CHECK_COLLISIONS
  for (int i = 0; i < bits_per_iwmfingerprint; i++)
  {
    if (nullptr != _first_path[i])
      delete [] _first_path[i];
  }
  
  delete [] _first_path;
#endif

#ifdef COUNT_TIMES_ATOM_IN_PATH
  delete [] _times_in_path;
#endif

  return;
}

void
MFingerprint::set_min_heteroatoms_at_path_ends (int m)
{
  assert (m >= 0 && m <= 2);

  _min_heteroatoms_at_ends = m;

  return;
}

/*int
MFingerprint::_compute_topotorsion_atom_hash_value (atom_number_t zatom) const
{
  atomic_number_t z = _atomic_number[zatom];

  int ncon = _ncon[zatom];

  int pi;
  if (! _atoms[zatom]->pi_electrons(pi))
  {
    cerr << "MFingerprint::_compute_topotorsion_atom_hash_value: cannot compute pi electrons for atom " << zatom << ", type " << _atoms[zatom]->atomic_symbol() << ", ncon = " << _ncon[zatom] << endl;
    pi = 0;
  }

  return z * 104 + 7 * ncon + pi;
}*/

/*
  Compute the hash value for a single atom
  
*/

int
MFingerprint::_compute_atom_hash_value (atom_number_t zatom) const
{
//if (topotorsion_atom_types)
//  return _compute_topotorsion_atom_hash_value(zatom);

  atomic_number_t z = _atomic_number[zatom];

// An aliphatic carbon can have any kind of bond attached, whereas an
// aromatic carbon can have only arom, single or double

  if (6 == z)
  {
    if (_arom[zatom])
      return 0;
    else
      return 3;
  }

// Number the bond types 1 == aromatic, 0 == single, and we don't
// have to leave extra space for bonds to an unsaturated atom - all
// bonds to it are single bonds.

  if (7 == z)
  {
    if (_arom[zatom])
      return 7;
    else
      return 10;
  }

// only one bond type to unsaturated Nitrogen

  if (8 == z)
  {
    if (_arom[zatom])
      return 14;
    else
      return 17;
  }

// only one bond type to unsaturated Oxygen

  if (9 == z)
    return 20;

  if (17 == z)
    return 21;

  if (35 == z)
    return 22;

  if (53 == z)
    return 23;

  if (16 == z)      // sulphur
  {
    if (! _arom[zatom])
      return 24;
    else
      return 27;
  }

// leave room for 2 bond types to aliphatic sulphur

  if (15 == z)      // phosphorus
    return 28;

  if (5 == z)
    return 29;

  if (14 == z)
    return 30;

  return _compute_atom_hash_value_non_periodic_table(zatom, 31);
}

int
MFingerprint::_compute_atom_hash_value_include_unsaturation (atom_number_t zatom) const
{
//if (topotorsion_atom_types)
//  return _compute_topotorsion_atom_hash_value(zatom);

  atomic_number_t z = _atomic_number[zatom];

  const Atom * a = _atoms[zatom];

  int unsaturation;

  if (_arom[zatom])
    unsaturation = 0;
  else
    unsaturation = a->nbonds() - _ncon[zatom];

// aromatic carbon can have only arom, single or double (3)
// an unsaturated carbon can have single, double or triple (3)
// a saturated carbon can have only single bonds (1)

  if (6 == z)
  {
    if (_arom[zatom])
      return 0;
    else if (unsaturation)
      return 3;
    else
      return 7;
  }

// Aromatic N can have aromatic, single and double bonds (3)
// Unsaturated N can have single, double, triple (3)
// A saturated N can have only single bonds (1)

  if (7 == z)
  {
    if (_arom[zatom])
      return 8;
    else if (unsaturation)
      return 11;
    else
      return 15;
  }

// aromatic oxygen can have only aromatic bonds - remember, aromatic bonds are number 1
// unsaturated oxygen can have single or double (2)
// Saturated oxygen can have only single bonds (1)

  if (8 == z)
  {
    if (_arom[zatom])   // only aromatic bonds, but aromatic is 1
      return 16;
    else if (unsaturation)
      return 18;
    else
      return 20;
  }

  if (9 == z)
    return 21;

  if (17 == z)    // there is no unsaturated Cl in drugs
    return 22;

  if (35 == z)
    return 23;

  if (53 == z)
    return 24;

// Aromatic sulphur has only aromatic (numbered 1) (1)
// Unsaturated sulphur has single, double (2)
// Saturated sulphur has single bonds (1)

  if (16 == z)      // sulphur
  {
    if (_arom[zatom])
      return 25;
    else if (unsaturation)
      return 27;
    else
      return 29;
  }

  if (15 == z)      // phosphorus
    return 30;

  if (5 == z)
    return 31;

  if (14 == z)
    return 32;

  return _compute_atom_hash_value_non_periodic_table (zatom, 33);
}

int
MFingerprint::_compute_atom_hash_value_non_periodic_table (atom_number_t zatom,
                                int hstart) const
{
  atomic_number_t z = _atomic_number[zatom];

  if (kInvalidAtomicNumber == z)
  {
    int h = _atoms[zatom]->element()->atomic_symbol_hash_value();

    if (h < 0)    // not a one or two letter element!
      return hstart + HIGHEST_ATOMIC_NUMBER;

//  May 2004. Without the + 9, fingerprinting the 'a' atom didn't work properly - collided with something

    return hstart + 9 + HIGHEST_ATOMIC_NUMBER + 2 * h + _arom[zatom];
  }

  if (z <= 0)   // not sure what other cases there are
  {
    return hstart + HIGHEST_ATOMIC_NUMBER - z;
  }

  if (z <= 18)
    return hstart + 17 * z;

  return hstart + z + z;    // all others
}

/*
  Just like the regular path bit setting, except that when dealing with
  either of the end atoms, if they are charged, we key on their formal
  charge (if present);
*/

template <typename T>
void
MFingerprint::_do_formal_charge_bits (int istart, int istop, int istep,
                                      T & bits) const
{
  atom_number_t astart = _path[istart];
  atom_number_t astop  = _path[istop];

  unsigned int zbit = 931 + 7823 * _atoms[astart]->formal_charge() + 32771 * _atoms[astop]->formal_charge();

  istart += istep;    // we have already accounted for the first atom

  for (int i = istart; i != istop; i += istep)
  {
    zbit = 11517 * zbit + _path_hash_value[i] + _bond[i - 1];
  }

  bits[zbit]++;   BIT_FORMED(zbit, bits, "formal charge")

#ifdef DEBUG_DO_FORMAL_CHARGE_BITS
  cerr << "Formal charge bit " << (zbit % bits_per_iwmfingerprint) << endl;
#endif

  return;
}

/*
  We need to ensure that formal charge paths are treated in a canonical form
*/

template <typename T>
void
MFingerprint::_do_formal_charge_bits (atom_number_t astart,
                                      atom_number_t astop,
                                      T & bits) const
{
//cerr << " astart " << astart << " path[0] " << _path[0] << ", path length " << _path_length << endl;
//cerr << " astop " << astop << " path[_path_length] " << _path[_path_length] << endl;
  assert ((astart == _path[0] && astop == _path[_path_length]) || (astart == _path[_path_length] && astop == _path[0]));

  formal_charge_t fc0 = _atoms[astart]->formal_charge();

  if (0 == _path_length)
  {
    unsigned int b;

    if (fc0 >= 0)
      b = (767772 * fc0);
    else
      b = (-21118 * fc0);

#ifdef DEBUG_DO_FORMAL_CHARGE_BITS
    cerr << "Zero path length formal charge bit " << b << endl;
#endif

    bits[b]++;   BIT_FORMED(b, bits, "formal charge zero")

    return;
  }

  formal_charge_t fc1 = _atoms[astop]->formal_charge();

  if (fc0 > fc1)
    _do_formal_charge_bits(0, _path_length, 1, bits);
  else if (fc0 < fc1)
    _do_formal_charge_bits(_path_length, 0, -1, bits);
  else if (_path_hash_value[0] > _path_hash_value[_path_length])
    _do_formal_charge_bits(0, _path_length, 1, bits);
  else if (_path_hash_value[0] < _path_hash_value[_path_length])
    _do_formal_charge_bits(_path_length, 0, -1, bits);
  else
    _do_formal_charge_bits(0, _path_length, 1, bits);

  return;
}

//#define DEBUG_DO_ISOTOPE_BITS

template <typename T>
void 
MFingerprint::_do_isotope_bits (int istart,
                                int istop,
                                int istep,
                                T & bits) const
{
  atom_number_t astart = _path[istart];
  atom_number_t astop  = _path[istop];

  unsigned int zbit = 3313 * _isotope[astart] + 7 * _isotope[astop];

  istart += istep;    // we have already accounted for the first atom

  int d = 0;    // measure of distance into the vector, regardless of direction

  for (int i = istart; i != istop; i += istep)
  {
    zbit = 11517 * zbit + _path_hash_value[i] + _bond[i - 1];

    atom_number_t j = _path[i];

    d++;

    if (0 != _isotope[j])
      zbit += 17;
  }

#ifdef DEBUG_DO_ISOTOPE_BITS
  cerr << "MFingerprint::_do_isotope_bits:setting " << (zbit % bits_per_iwmfingerprint) << endl;
#endif

  bits[zbit]++;      BIT_FORMED(zbit, bits, "isotope")
}

template <typename T>
void
MFingerprint::_do_isotope_bits (T & bits) const
{
  assert (nullptr != _isotope);

  atom_number_t astart = _path[0];

  isotope_t i0 = _isotope[astart];
  if (0 == _path_length)
  {
    unsigned int b = 77265 + _path_hash_value[0] * 87 + i0 % 87;

    bits[b]++;      BIT_FORMED(b, bits, "isotope 0")

#ifdef DEBUG_DO_ISOTOPE_BITS
  cerr << "MFingerprint::_do_isotope_bits:setting " << (b % bits_per_iwmfingerprint) << endl;
#endif

    return;
  }

  atom_number_t astop = _path[_path_length];

  isotope_t i1 = _isotope[astop];

  if (i0 > i1)
    _do_isotope_bits(0, _path_length, 1, bits);
  else if (i0 < i1)
    _do_isotope_bits(_path_length, 0, -1, bits);
  else if (0 == i0 && 0 == i1)
    ;
  else if (_path_hash_value[0] < _path_hash_value[_path_length])
    _do_isotope_bits(0, _path_length, 1, bits);
  else
    _do_isotope_bits(_path_length, 0, -1, bits);

  return;
}

//#define DEBUG_SET_CLUSTER_BIT

template <typename T>
void
MFingerprint::_set_cluster_bit (atom_number_t centre,
                                int b1, atom_number_t a1,
                                int b2, atom_number_t a2,
                                int b3, atom_number_t a3,
                                T & bits)
{
  int h1 = _atom_hash_value[a1] + b1;
  int h2 = _atom_hash_value[a2] + b2;
  int h3 = _atom_hash_value[a3] + b3;

#ifdef DEBUG_SET_CLUSTER_BIT
  cerr << "Setting cluster bit " << centre << " a1 = " << a1 << " (" << h1 << ") a2 = " << a2 << " (" << h2 << ") a3 = " << a3 << " (" << h3 << ")\n";
#endif

// Enforce a canonical order on these.

  if (h1 > h3)
  {
    int tmp = h1;
    h1 = h3;
    h3 = tmp;
  }

  if (h1 > h2)
  {
    int tmp = h1;
    h1 = h2;
    h2 = tmp;
  }

  if (h2 > h3)
  {
    int tmp = h2;
    h2 = h3;
    h3 = tmp;
  }

  assert (h1 <= h2 && h2 <= h3);

  unsigned int zbit = _atom_hash_value[centre] + 15 * h1 + 911 * h2 + 3741 * h3;

  bits[zbit]++;      BIT_FORMED(zbit, bits, "cluster")

#ifdef DEBUG_SET_CLUSTER_BIT
  cerr << "Cluster bit " << (zbit % nbits) << endl;
#endif

  return;
}

int
Linear_Fingerprint_Defaults::numeric_bond_code(const Bond * b) const
{
  if (b->is_aromatic())
    return 1;

  int rc ;

  if (b->is_single_bond())
    rc = 0;
  else if (b->is_double_bond())
    rc = 2;
  else if (b->is_triple_bond())
    rc = 3;
  else if (IS_COORDINATION_BOND(b->btype()))
    rc = 4;
  else 
    rc = 5;

  if (_differentiate_ring_and_chain_bonds && b->nrings())
    rc += 18;

  return rc;
}

template <typename T>
void
MFingerprint::_build_branched_paths (const Linear_Fingerprint_Defaults * lfpd, 
                                     const int * include_these_atoms,
                                     T & bits)
{
  assert (2 == _path_length);

  const atom_number_t a1 = _path[1];
  if (2 == _ncon[a1])        // middle atom is not branched, cannot do anything
    return;

  const Atom * centre = _atoms[a1];

  const atom_number_t a0 = _path[0];
  const atom_number_t a2 = _path[2];

  for (int i = 0; i < _ncon[a1]; i++)
  {
    const Bond * b = centre->item(i);

    atom_number_t j = b->other(a1);
    if (a0 == j || a2 == j)     // we are looking for atoms outside the path
      continue;

    if (0 == include_these_atoms[j])
      continue;

    _set_cluster_bit(a1, _bond[0], a0, _bond[1], a2, lfpd->numeric_bond_code(b), j, bits);
  }

  return;
}

//#define DEBUG_SET_AUX_BITS

#ifdef DEBUG_SET_AUX_BITS
#define PRINT_BIT(b, reason) { cerr << "setting " << (reason) << " bit " << (b) << endl; }
#else
#define PRINT_BIT(b, reason)
#endif

#define INCREMENT_VECTOR(v, b, reason) \
  { \
    PRINT_BIT(b, reason) \
    (v)[(b)]++;\
  }

template <typename T>
void
MFingerprint::_set_auxiliary_bits (unsigned int zbit1, unsigned int zbit2,
                                   int is_ring,
                                   int astart, int astop,
                                   const Linear_Fingerprint_Defaults * lfpd,
                                   T & bits) const
{
#ifdef DEBUG_SET_AUX_BITS
  cerr << "Setting auxiliary bits\n";
#endif

  if (_path_length <= lfpd->bits_for_hydrogen_attachments())
  {
    if (_hcount[astart])
    {
      const unsigned int hbit = 6 + zbit2 + 3131;
      bits[hbit]++;       BIT_FORMED(hbit, bits, "hydrogen0")
    }
    if (_hcount[astop])
    {
      const unsigned int hbit = 6 + zbit1 + 8181;
      bits[hbit]++;     BIT_FORMED(hbit, bits, "hydrogen0")
    }
  }

  if (is_ring)
  {
    unsigned int rbit1 = (zbit1 * 93717 * _path_length);
    bits[rbit1]++;      BIT_FORMED(rbit1, bits, "ring path")
  }
//else if (_nrings[astart] || _nrings[astop])
  else if (1 == _path_length && (_nrings[astart] || _nrings[astop]))
  {
    unsigned int rbit = (zbit1 + zbit2 + 119);
    bits[rbit]++;       BIT_FORMED(rbit, bits, "ring path 1")

//  if (_nrings[astart] > 1 || _nrings[astop] > 1)    // could never get this to work
//  {
//    rbit = (rbit * 31) % nbits;
//    bits[rbit]++;
//  }
  }

// do something for branching. Make sure you keep the numbers the same in
// both places. Note that since these depend on ncon, they break the
// substructure relation

  if (lfpd->only_bits_preserving_substructure_perception())
    ;
  else if (_ncon[astart] > 2 && 6 != _atomic_number[astart])
  {
    unsigned int cbit = (zbit1 + 1024) >> 1;

    bits[cbit]++;     BIT_FORMED(cbit, bits, "ncon2")

    if (_ncon[astart] > 3)
    {
      const unsigned int cb = ((cbit + zbit2) & 8715);
      bits[cb]++;     BIT_FORMED(cb, bits, "ncon3")
    }

    if (0 == _arom[astart] && _unsaturation[astart])
    {
      const unsigned int bu = cbit + 31290;
      bits[bu]++;    BIT_FORMED(bu, bits, "unsat0")
    }
  }

  if (lfpd->only_bits_preserving_substructure_perception())
    ;
  else if (_path_length > 1 && _ncon[astop] > 2 && 6 != _atomic_number[astop])
  {
    unsigned int cbit = (zbit1 + 8723) >> 1;
    bits[cbit]++;  BIT_FORMED(cbit, bits, "ncon2")

    if (_ncon[astop] > 3)
    {
      const unsigned int x = (cbit + zbit1) & 8715;
      bits[x]++;         BIT_FORMED(x, bits, "ncon3")
    }

    if (0 == _arom[astop] && _unsaturation[astop])
    {
      const unsigned int x = cbit + 31290;
      bits[x]++;     BIT_FORMED(x, bits, "unsat")
    }
  }

// and unsaturation

  if (0 == _arom[astart] && _unsaturation[astart])
  {
    const unsigned int ubit = (zbit2 * 13);
    bits[ubit]++;    BIT_FORMED(ubit, bits, "unsat0")
  }

  if (_path_length > 0 && 0 == _arom[astop] && _unsaturation[astop])
  {
    unsigned int ubit = (zbit2 * 13);
    bits[ubit]++;    BIT_FORMED(ubit, bits, "unsat")
  }

  if (_path_length > 2 && (6 != _atomic_number[astart] || 6 != _atomic_number[astop]))
  {
    unsigned int unsat = 0;
    unsigned int branches = 0;

    for (int i = 1; i < _path_length; i++)
    {
      atom_number_t a = _path[i];
      branches += _ncon[a] - 2;
      if (0 == _arom[a])
        unsat += _unsaturation[a];
    }

    if (unsat)
    {
      const unsigned int x = unsat + 387712;
      bits[x]++;     BIT_FORMED(x, bits, "unsat")
    }

    if (lfpd->only_bits_preserving_substructure_perception())
      ;
    else if (branches)
    {
      const unsigned int x = static_cast<unsigned int>(_path_length) * branches + zbit2%6;
      bits[x]++;     BIT_FORMED(x, bits, "branches")
    }
  }

  if (_path_length > lfpd->set_formal_charge_bits())   // it gets initialised to -1
    ;
  else if (_atoms[astart]->formal_charge() || _atoms[astop]->formal_charge())
    _do_formal_charge_bits(astart, astop, bits);

  if (_path_length > lfpd->max_path_length_isotopic_bits())   // it gets initialised to -1
    ;
  else if (_isotope[astart] || _isotope[astop])
    _do_isotope_bits(bits);

  return;
}

//#define DEBUG_SET_BIT

/*
  Within both _set_bit_forward and _set_bit_backward, the coefficients
  used must be the same. Perhaps the functions could be combined...
*/

#define IWMFP_PBS1 31
#define IWMFP_PBS2 4703
#define IWMFP_XBS1 31
#define IWMFP_XBS2 4703

template <typename T>
void
MFingerprint::_set_bit_forward (int is_ring,
                                const Linear_Fingerprint_Defaults * lfpd,
                                T & bits) const
{
#ifdef DEBUG_SET_BIT
  cerr << "Set forward\n";
#endif

  unsigned int zbit1 = IWMFP_PBS1 * _path_length + _path_hash_value[0];
  unsigned int zbit2 = IWMFP_PBS2 * _path_length + _path_hash_value[0];

#ifdef COUNT_TIMES_ATOM_IN_PATH
//cerr << "Registering path of length " << _path_length << ", 0 " << _path[0] << endl;
  for (auto i = 0; i <= _path_length; ++i)
  {
//  cerr << " i = " << i << " item " << _path[i] << endl;
    _times_in_path[_path[i]]++;
  }
#endif

//cerr << "Initial " << zbit1 << " and " << zbit2 << endl;

  for (int i = 1; i <= _path_length; i++)
  {
    zbit1 = IWMFP_XBS1 * zbit1 + _path_hash_value[i] + _bond[i - 1];
    zbit2 = IWMFP_XBS2 * zbit2 + _path_hash_value[i] + _bond[i - 1];
//  cerr << "Now " << zbit1 << " and " << zbit2 << endl;
  }

  bits[zbit1]++;    BIT_FORMED(zbit1, bits, "forward")
  if (_path_length > 2 && zbit2 != zbit1)
  {
    bits[zbit2]++;   BIT_FORMED(zbit2, bits, "forward2")
  }

#ifdef DEBUG_SET_BIT
  cerr << "path of length " << _path_length << " sets bits " << zbit1 << " and " << zbit2 << endl;
#endif

  if (! lfpd->maximal_daylight_compatibility())
    _set_auxiliary_bits(zbit1, zbit2, is_ring, _path[0], _path[_path_length], lfpd, bits);

  return;
}

template <typename T>
void
MFingerprint::_set_bit_backward (int is_ring,
                                 const Linear_Fingerprint_Defaults * lfpd,
                                 T & bits) const
{
#ifdef DEBUG_SET_BUT
  cerr << "Set backward " << _path_hash_value[0] << " " << _path_hash_value[1] << endl;
#endif

#ifdef COUNT_TIMES_ATOM_IN_PATH
//cerr << "Registering path of length " << _path_length << ", 0 " << _path[0] << endl;
  for (auto i = 0; i <= _path_length; ++i)
  {
//  cerr << " i = " << i << " item " << _path[i] << endl;
    _times_in_path[_path[i]]++;
  }
#endif

  unsigned int zbit1 = IWMFP_PBS1 * _path_length + _path_hash_value[_path_length];
  unsigned int zbit2 = IWMFP_PBS2 * _path_length + _path_hash_value[_path_length];

  for (int i = _path_length - 1; i >= 0; i--)
  {
    zbit1 = IWMFP_XBS1 * zbit1 + _path_hash_value[i] + _bond[i];
    zbit2 = IWMFP_XBS2 * zbit2 + _path_hash_value[i] + _bond[i];
  }

  bits[zbit1]++;  BIT_FORMED(zbit1, bits, "backward")

  if (_path_length > 2 && zbit2 != zbit1)
  {
    bits[zbit2]++; BIT_FORMED(zbit2, bits, "backward2")
  }

#ifdef DEBUG_SET_BIT
  cerr << "path of length " << _path_length << " sets bits " << zbit1 << " and " << zbit2 << endl;
#endif

  if (! lfpd->maximal_daylight_compatibility())
    _set_auxiliary_bits(zbit1, zbit2, is_ring, _path[_path_length], _path[0], lfpd, bits);

  return;
}

/*
  A complete path is in the object. Set the appropriate bit
  The first thing to do is decide on the direction in order to
  get a canonical path
*/

template <typename T>
void
MFingerprint::_set_bit (int astart,
                        const Linear_Fingerprint_Defaults * lfpd,
                        T & bits) const
{
#ifdef DEBUG_SET_BIT
  cerr << "Setting bit for path length " << _path_length << endl;
  for (int i = 0; i < _path_length; i++)
  {
    cerr << "Atom " << _path[i] << " type " << _atoms[_path[i]]->atomic_symbol() << " hash " << _path_hash_value[i] << " bond " << _bond[i] << endl;
  }
  cerr << "Atom " << _path[_path_length] << " type " << _atoms[_path[_path_length]]->atomic_symbol() << " hash " << _path_hash_value[_path_length] << endl;
#endif

// If there is a specification on the minimum number of heteroatoms at ends
// of the path, check that now

  if (0 == _min_heteroatoms_at_ends)
    ;
  else if ((6 != _atomic_number[_path[0]]) + (6 != _atomic_number[_path[_path_length]]) < _min_heteroatoms_at_ends)
    return;

  if (lfpd->isotopic_paths() && ! _atoms[_path[_path_length]]->is_isotope())
    return;

  if (0 == _path_length)
  {
    _set_bit_forward(0, lfpd, bits);      // last arg 0 means not a ring path
    return;
  }

  int is_ring;
  if (_path_length < 2)
    is_ring = 0;
  else if (_in_path[astart])
    is_ring = 1;
  else
    is_ring = 0;

// Need to figure out whether this path gets recorded in the forward or
// backward direction.

  int i1 = 0;
  int i2 = _path_length;

  int direction = 1;

  while (i2 > i1)
  {
    if (_path_hash_value[i1] < _path_hash_value[i2])
    {
      direction = 1;
      break;
    }
    if (_path_hash_value[i1] > _path_hash_value[i2])
    {
      direction = -1;
      break;
    }

    if (_bond[i1] > _bond[i2 - 1])
    {
      direction = 1;
      break;
    }
    if (_bond[i1] < _bond[i2 - 1])
    {
      direction = -1;
      break;
    }

    i1++;
    i2--;
  }

  if (direction > 0)
  {
    _set_bit_forward(is_ring, lfpd, bits);
    if (lfpd->bits_for_atoms_in_multiple_rings() && _molecule_contains_fused_rings && _path_length > 2)
      _set_bits_for_atoms_in_multiple_rings(1, lfpd, bits);
  }
  else
  {
    _set_bit_backward(is_ring, lfpd, bits);
    if (lfpd->bits_for_atoms_in_multiple_rings() && _molecule_contains_fused_rings && _path_length > 2)
      _set_bits_for_atoms_in_multiple_rings(-1, lfpd, bits);
  }

// We only do the trick of leaving out CH2 groups when the path is length 3 or longer

  if (_path_length < 3 || 0 == lfpd->omit_ch2())
    return;

// If omit_ch2 is just 1, we only omit the CH2 groups if there are heteroatoms at the ends

  if (lfpd->omit_ch2() > 1)     // do it regardless
    ;
  else if (6 == _atomic_number[astart] || 6 == _atomic_number[_path[_path_length]])
    return;

// Look for first and last CH2 groups

  int first_ch2 = 0;
  int last_ch2;       // gcc complains about not initialised, that's OK.
  for (int i = 1; i < _path_length; i++)
  {
    if (_ch2[i])
    {
      if (0 == first_ch2)
        first_ch2 = i;

      last_ch2 = i;
    }
  }

#ifdef DEBUG_SET_BIT
  cerr << "First ch2 at " << first_ch2 << " last at " << last_ch2 << endl;
#endif

  if (first_ch2)
    _omit_ch2(first_ch2, last_ch2, bits);

  return;
}

/*
  Making the paths in which a CH2 group is omitted is very similar to making
  a regular path. 
  Our strategy is to scan inwards from the ends looking for the first CH2.
  Note that in order to preserve symmetry, we may need to generate two paths.
  For example

  X - CH2 - Y - Z - CH2 - W
  
  We need to generate X - Y - Z - CH2 -W and X - CH2 - Y - Z - W
*/

template <typename T>
void
MFingerprint::_omit_ch2 (int first_ch2, int last_ch2,
                         T & bits) const
{
  int d1 = first_ch2;
  int d2 = _path_length - last_ch2;

  if (d1 < d2)
    _omit_ch2(first_ch2, bits);
  else if (d1 > d2)
    _omit_ch2(last_ch2, bits);
  else
  {
    _omit_ch2(first_ch2, bits);
    _omit_ch2(last_ch2, bits);
  }

  return;
}

/*
  We are to set some bits based on a path which omits OMIT.
  First decide on a direction
*/

template <typename T>
void
MFingerprint::_omit_ch2 (int omit,
                         T & bits) const
{
  int direction = 1;

  int i1 = 1;
  int i2 = _path_length - 1;

  while (i2 > i1 && omit < 0)
  {
    if (omit == i1)
    {
      i1++;
      continue;
    }

    if (omit == i2)
    {
      i2--;
      continue;
    }

    if (_path_hash_value[i1] < _path_hash_value[i2])
    {
      break;
    }
    if (_path_hash_value[i1] > _path_hash_value[i2])
    {
      direction = -1;
      break;
    }

    if (_bond[i1] > _bond[i2 - 1])
    {
      break;
    }
    if (_bond[i1] < _bond[i2 - 1])
    {
      direction = -1;
      break;
    }

    i1++;
    i2--;
  }

  _set_bit_omit(direction, omit, bits);

  return;
}

/*
  Set a bit based on a path from which a member of the path is omitted.

  We need to keep this in-synch with _set_bit_forward and _set_bit_backward
  in order to set the same bit
*/

template <typename T>
void
MFingerprint::_set_bit_omit (int direction, int omit,
                             T & bits) const
{
#ifdef DEBUG_SET_BIT_OMIT
  cerr << "Omitting path member " << omit << endl;
#endif

  int istart, istop;

  if (direction > 0)
  {
    istart = 0;
    istop = _path_length + 1;
  }
  else
  {
    istart = _path_length;
    istop = -1;
  }

  unsigned int zbit1 = IWMFP_PBS1 * (_path_length - 1) + _path_hash_value[istart];
  unsigned int zbit2 = IWMFP_PBS2 * (_path_length - 1) + _path_hash_value[istart];

  for (int i = istart + direction; i != istop; i += direction)
  {
    if (omit == i)
      continue;

    zbit1 = IWMFP_XBS1 * zbit1 + _path_hash_value[i] + _bond[i - 1];
    zbit2 = IWMFP_XBS2 * zbit2 + _path_hash_value[i] + _bond[i - 1];
  }

  bits[zbit1]++;   BIT_FORMED(zbit1, bits, "omit")
  bits[zbit2]++;   BIT_FORMED(zbit2, bits, "omit")

#ifdef DEBUG_SET_BIT_OMIT
  cerr << "By omitting " << omit << " bits " << zbit1 << " and " << zbit2 << endl;
#endif

  return;
}

template <typename T>
void
MFingerprint::_set_bits_for_atoms_in_multiple_rings(int direction,
                                        const Linear_Fingerprint_Defaults * lfpd,
                                        T & bits) const
{
  int istart, istop, istep, bdelta;
  if (direction > 0)
  {
    istart = 0;
    istop = _path_length;
    istep = 1;
    bdelta = -1;
  }
  else
  {
    istart = _path_length - 1;
    istop = -1;
    istep = -1;
    bdelta = 0;
  }
  
  unsigned int atoms_in_multiple_rings = 0;

  unsigned int b = _path_hash_value[istart];

  for (int i = istart + istep; i != istop; i += istep)
  {
    atom_number_t j = _path[i];

    if (_nrings[j] > 1)
    {
      atoms_in_multiple_rings++;
      b += 152213 * _bond[i + bdelta] * (100 + _path_hash_value[i]);
    }
    else
      b += 28331 * _bond[i + bdelta] * (10 + _path_hash_value[i]);

//  cerr << " Atom " << j << " _nrings " << _nrings[j] << " bond " << _bond[i + bdelta] << " phv " << _path_hash_value[i] << "  b = " << b << endl;
  }

//if (atoms_in_multiple_rings)
//  cerr << "Path of length " << _path_length << " setting bit " << (b % nbits) << endl;
  if (atoms_in_multiple_rings)
  {
    const unsigned int x = b + 97 * atoms_in_multiple_rings * static_cast<unsigned int>(_path_length);
    bits[x]++;  BIT_FORMED(x, bits, "multiple rings")
  }

  return;
}

/*
  Build all fingerprints emanating from the current start atom
*/

template <typename T>
int
MFingerprint::_build (atom_number_t aprev, 
                      const int * include_these_atoms,
                      const Linear_Fingerprint_Defaults * lfpd,
                      T & bits)
{
  assert (_path_length >= 0);

  atom_number_t astart = _path[_path_length];

//cerr << "MFingerprint::_build:from " << aprev << " to " << astart << " _path_length " << _path_length << endl;

  assert(include_these_atoms[astart]);

  assert (astart >= 0 && astart < _matoms);     // must be a valid atom for our molecule

  if ((_path_length <= lfpd->set_formal_charge_bits() || lfpd->formally_charged_atoms_lose_identity()) &&
      _atoms[astart]->formal_charge())
  {
    _path_hash_value[_path_length] = 7900 + 977 * _atoms[astart]->formal_charge();
  }
  else
    _path_hash_value[_path_length] = _atom_hash_value[astart];

  const Atom * a = _atoms[astart];

  const int acon = a->ncon();

  _ch2[_path_length] = (6 == _atomic_number[astart] && 2 == acon && 2 == a->nbonds());

  const int min_path_length = lfpd->min_path_length();
  const int max_path_length = lfpd->max_path_length();

  if (_path_length >= min_path_length && _path_length <= max_path_length)   // in range, write the path
    _set_bit(astart, lfpd, bits);

  if (2 == _path_length && min_path_length <= 2)
    _build_branched_paths(lfpd, include_these_atoms, bits);

  if (_in_path[astart])     // we have formed a ring, don't continue from here
  {
    if (lfpd->place_bits_for_rings() && astart == _path[0])
      _do_ring_bits(lfpd, bits);
    return 1;
  }

  if (_path_length <= lfpd->terminal_path_bits())
    _set_terminal_path_bits(bits);

  if (_path_length == max_path_length)    // we are at the longest path allowed, cannot extend it
    return 1;

//cerr << "_build from " << aprev << " to " << astart << " nset " << bits.nset() << endl;

// Looks like this path is continuing...

  _in_path[astart] = 1;

  _path_length++;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(astart);

    if (aprev == j)      // cannot turn back on ourselves
      continue;

    if (0 == include_these_atoms[j])
      continue;

    _bond[_path_length - 1] = lfpd->numeric_bond_code(b);
    _path[_path_length] = j;

    _build(astart, include_these_atoms, lfpd, bits);
  }

  _path_length--;
  _in_path[astart] = 0;

  return 1;
}

template <typename T>
int
MFingerprint::build (atom_number_t astart, 
                     const int * include_these_atoms,
                     const Linear_Fingerprint_Defaults * lfpd, 
                     T & bits)
{
  if (2 == _min_heteroatoms_at_ends && 6 == _atomic_number[astart])
    return 1;

  if (lfpd->isotopic_paths() && ! _atoms[astart]->is_isotope())
    return 1;

  _path_length = 0;

  set_vector(_in_path, _matoms, 0);

  if (! include_these_atoms[astart])
  {
    cerr << "MFingerprint::build:start atom not included\n";
    abort();
  }

  _path[0] = astart;

  return _build(INVALID_ATOM_NUMBER, include_these_atoms, lfpd, bits);
}

template <typename T>
int
MFingerprint::fingerprint_large_ring (const Ring & r,
                                      const Linear_Fingerprint_Defaults * lfpd,
                                      T & bits)
{
  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    _path[i] = r[i];
  }

  _path_length = ring_size;

  _do_ring_bits(lfpd, bits);

  return 1;
}

/*
  We have formed a ring
*/

template <typename T>
void
MFingerprint::_do_ring_bits (const Linear_Fingerprint_Defaults * lfpd, 
                             T & bits) const
{
  int ring_size = _path_length;

//#define DEBUG_DO_RING_BITS
#ifdef DEBUG_DO_RING_BITS
  cerr << "MFingerprint::_do_ring_bits:found ring of size " << ring_size << endl;
  for (int i = 0; i < _path_length; i++)
  {
    cerr << "Atom " << _path[i] << endl;
  }
#endif

  const int n = lfpd->place_bits_for_rings();

  for (int i = 0; i < n; i++)
  {
    unsigned int b = ((i + 1) * (_path_length) * 327);

#ifdef DEBUG_DO_RING_BITS
    cerr << "ring_size " << ring_size << " setting bit " << b << " currently " << bits[b] << endl;
#endif

    bits[b]++; BIT_FORMED(b, bits, "ring bits")
  }

  int heteroatoms_in_ring = 0;

  for (int i = 0; i < _path_length; i++)
  {
    atom_number_t j = _path[i];

    if (6 != _atoms[j]->atomic_number())
      heteroatoms_in_ring++;
  }

  unsigned int b = (88 * heteroatoms_in_ring + ring_size) * 7691;

#ifdef DEBUG_DO_RING_BITS
  cerr << " b = " << b << " -> " << b << endl;
#endif

  bits[b]++;    BIT_FORMED(b, bits, "ring hetero")

  return;
}

/*
  We make a simplifying assumption by only fingerprinting terminal
  paths where one end is terminal. Saves an expensive canonicalisation
*/

template <typename T>
void
MFingerprint::_set_terminal_path_bits (T & bits)const
{
  atom_number_t a0 = _path[0];
  atom_number_t a1 = _path[_path_length];

  if (_ncon[a0] > 1 && _ncon[a1] > 1)    // neither end of the path is terminal
    return;

// At least one of the groups is terminal. Only do atoms which could be something other than terminal

  int t0 = 0;
  if (1 == _ncon[a0] && _atoms[a0]->element()->normal_valence() > 1)
    t0 = 1;

  int t1 = 0;
  if (1 == _ncon[a1] && _atoms[a1]->element()->normal_valence() > 1)
    t1 = 1;

// Only do it if just one end is terminal

  if (t0 && 0 == t1)
    _set_terminal_path_bits(1, bits);
  else if (0 == t0 && t1)
    _set_terminal_path_bits(-1, bits);

  return;
}

/*
  Convoluted logic to ensure that we step through the path
  in a consistent manner. To figure this out, just write out
  a couple of simple cases...

  Atom0 Bond0 Atom1 Bond1 Atom2 Bond2 Atom3

  Going left to right, we initialise with Atom0, and then
  take pairs (Bond0 Atom1) (Bond1 Atom2) (Bond2 Atom3)
  Going right to left, we initialise with Atom3 and then
  take pairs (Bond2 Atom2) (Bond1 Atom1) (Bond0 Atom0)
*/

template <typename T>
void
MFingerprint::_set_terminal_path_bits (const int direction,
                                       T & bits) const
{
  unsigned int b;
  int atom_offset;

  int istart, istop;
  if (direction > 0)
  {
    istart = 0;
    istop = _path_length;
    b = _path_hash_value[0];
    atom_offset = 1;
  }
  else
  {
    istart = _path_length - 1;
    istop = -1;
    b = _path_hash_value[_path_length];
    atom_offset = 0;
  }

  for (int i = istart; i != istop; i += direction)
  {
    b = b * (_bond[i] + 85) + _path_hash_value[i + atom_offset];
  }

  const unsigned int x = 3 * b + static_cast<unsigned int>(_path_length) + 1;
  bits[x]++;    BIT_FORMED(x, bits, "terminal")

  return;
}

#ifdef COMPILE_WITH_WATCH_BITS
int
MFingerprint::write_state (const Molecule & m,
                           const unsigned int b,
                           const char * reason,
                           std::ostream & output) const
{
  cerr << "MFingerprint::write_state:bit " << b << " path of length " << _path_length << " " << reason << ", in " << m.name() << endl;
  for (int i = 0; i < _path_length; ++i)
  {
    output << _path[i] << ' ' << m.smarts_equivalent_for_atom(_path[i]) << '\n';
  }

  return 1;
}
#endif

#ifdef COUNT_TIMES_ATOM_IN_PATH
void
MFingerprint::place_labels_for_times_in_path(Molecule & m) const
{
  for (auto i = 0; i < _matoms; ++i)
  {
    m.set_isotope(i, _times_in_path[i]);
  }

  return;
}
#endif

// Haven't implemented this yet because it isn't being used

#if 0
int
IWMFingerprint::construct_fingerprint (Molecule & m,
                                       Sparse_Fingerprint_Creator & sfp)
{
  unsigned int bsave = bits_per_iwmfingerprint;

  unsigned int mynbits = 1000007;    // just some large arbitrary number

  if (bits_per_iwmfingerprint <= 100000)
    set_iwmfingerprint_nbits(mynbits);

  int rc = construct_fingerprint(m);

  if (0 == rc)
    return 0;

  for (unsigned int i = 0; i < mynbits; i++)
  {
    sfp.hit_bit(i, _bvector[i]);
  }

  if (mynbits > bsave)
    bits_per_iwmfingerprint = bsave;

  return rc;
}
#endif

template <typename T, typename B>
int
Linear_Fingerprint_Creator::process (Molecule & m,
                                     const T * atype,
                                     const int * include_these_atoms,
                                     B & bits)
{
#ifdef COMPILE_WITH_WATCH_BITS
  current_molecule = &m;
#endif

//bits.clear();     is this a bug or a feature. Turned off for now because of rxn_fingerprint

  const int matoms = m.natoms();

  if (0 == matoms)
    return 1;

  MFingerprint mfp(m, reinterpret_cast<const int *>(atype), _lfpd);   // get the proper type matching sometime...

  if (_lfpd->min_heteroatoms_at_path_ends() >= 0)
    mfp.set_min_heteroatoms_at_path_ends(_lfpd->min_heteroatoms_at_path_ends());

// If we are fingerprinting rings, capture any that are longer than our max path length

  if (_lfpd->place_bits_for_rings())
  {
    const int nr = m.nrings();
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);

      if (ri->number_elements() >= _lfpd->max_path_length())
        mfp.fingerprint_large_ring(*ri, _lfpd, bits);
    }
  }

  for (int i = 0; i < matoms; i++)
  {
//  cerr << "atom " << i << " include_these_atoms " << include_these_atoms[i] << endl;
    if (include_these_atoms[i])
      (void) mfp.build(i, include_these_atoms, _lfpd, bits);
  }

#ifdef COUNT_TIMES_ATOM_IN_PATH
  mfp.place_labels_for_times_in_path(m);
#endif

  return 1;
}

template <typename T, typename B>
int
Linear_Fingerprint_Creator::process (Molecule & m,
                                     const T * atype,
                                     B & bits)
{
  const int matoms = m.natoms();

  if (0 == matoms)    
  {
    bits.clear();
    return 1;
  }

  int * include_these_atoms = new_int(matoms, 1); std::unique_ptr<int[]> free_include_these_atoms(include_these_atoms);

  return process(m, atype, include_these_atoms, bits);
}

int
display_misc_fingerprint_options (std::ostream & output,
                                  char flag)
{
  output << " -" << flag << " nb=<n>   number of bits in each fingerprint\n";
  output << " -" << flag << " unsinca  aromatic atoms are considered unsaturated\n";
  output << " -" << flag << " unsat    include unsaturation state in atom type\n";
//output << " -" << flag << " tt       use Topotorsion types for atom types\n";
  output << " -" << flag << " ap       set bits corresponding to atom pairs\n";
  output << " -" << flag << " iatfc    initial atom types forgotten for charged atoms\n";
  output << " -" << flag << " mdc      maximal Daylight compatibility\n";
  output << " -" << flag << " rb=<n>   extra bits set for rings of size <= <n>\n";
  output << " -" << flag << " hb=<l>   set bits for attached hydrogens through path length <l>\n";
  output << " -" << flag << " ip       only do paths anchored with an isotopic atom\n";
  output << " -" << flag << " ip=<l>   set extra bits for isotopic paths up to length <l>\n";
  output << " -" << flag << " tpb=<n>  set extra bits for paths involving terminal atoms\n";
  output << " -" << flag << " amr=<n>  set extra bits when there are atoms in multiple rings\n";
  output << " -" << flag << " ss       only set bits that fully preserve substructure perception\n";

  return output.good();
}

int
Linear_Fingerprint_Defaults::build (Command_Line & cl,
                                    char flag,
                                    const int verbose)
{
  if (! cl.option_present(flag))
    return 1;

  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y); ++i)
  {
    const_IWSubstring directive;    // several of our things look like XX=nn
    int dvalue;
    int dvalue_valid;

    if (y.contains('='))
      dvalue_valid = y.split_into_directive_and_value(directive, '=', dvalue);
    else
      dvalue_valid = 0;

    if ("och2" == directive)
    {
      _omit_ch2 = dvalue;
      if (verbose)
        cerr << "Omit ch2 set to " << dvalue << endl;
    }
    else if ("iatfc" == y)
    {
      _formally_charged_atoms_lose_identity = 1;
      if (verbose)
        cerr << "Formally charged atoms fingerprinted only as charged form\n";
    }
    else if ("mdc" == y)
    {
      _maximal_daylight_compatibility = 1;
      if (verbose)
        cerr << "Will generate fingerprints similar to Daylight\n";
    }
    else if (y.starts_with("rb="))
    {
      if (! dvalue_valid || dvalue < 1)
      {
        cerr << "The number of ring bits to set (rb=) must be a whole positive number\n";
        return 0;
      }

      _place_bits_for_rings = dvalue;
      if (verbose)
        cerr << "Will set " << dvalue << " bits for each ring\n";
    }
    else if (y.starts_with("hb="))
    {
      if (! dvalue_valid || dvalue < 1)
      {
        cerr << "The add Hydrogen bits value must be a whole +ve number\n";
        return 0;
      }

      _bits_for_hydrogen_attachments = dvalue;

      if (verbose)
        cerr << "Attached hydrogen bits set for paths to length " << dvalue << endl;
    }
    else if (y.starts_with("mpl="))
    {
      if (! dvalue_valid || dvalue < 2)
      {
        cerr << "The max path length directive must be followed by a whole non-negative number\n";
        return 0;
      }

      _max_path_length = dvalue;
      if (verbose)
        cerr << "Max path length " << dvalue << endl;
    }
    else if (y.starts_with("ip="))
    {
      if (! dvalue_valid || dvalue < 0)
      {
        cerr << "The isotopic bits directive must be followed by a whole non-negative number\n";
        return 0;
      }

      _max_path_length_isotopic_bits = dvalue;
      if (verbose)
        cerr << "Extra bits set for isotopic paths less than " << dvalue << " atoms\n";
    }
    else if (y.starts_with("tpb="))
    {
      if (! dvalue_valid || dvalue < 0)
      {
        cerr << "The terminal path bits directive must be followed by a whole non-negative number\n";
        return 0;
      }

      _terminal_path_bits = dvalue;
      if (verbose)
        cerr << "Extra bits for terminal paths up to length " << dvalue << " bonds\n";
    }
    else if (y.starts_with("amr="))
    {
      if (! dvalue_valid || dvalue < 0)
      {
        cerr << "The atoms in multiple rings directive must be followed by a whole non-negative number\n";
        return 0;
      }

      _bits_for_atoms_in_multiple_rings = dvalue;
      if (verbose)
        cerr << "Will set bits for multiple ring membership " << dvalue << " bonds\n";
    }
    else if ("ss" == y)
    {
      _only_bits_preserving_substructure_perception = 1;
      if (verbose)
        cerr << "Only bits preserving substructure relationships generated\n";
    }
    else if ("help" == y)
    {
      display_misc_fingerprint_options(cerr, 'Y');
      return 0;
    }
    else 
    {
      cerr << "Linear_Fingerprint_Defaults::parse_misc_fingerprint_options::unrecognised -" << flag << " qualifier '" << y << "'\n";
      return 0;
    }
  }

  return 1;
}

int
parse_misc_fingerprint_options (Command_Line & cl,
                                char flag,
                                int verbose)
{
  return default_linear_fingerprint_defaults.build(cl, flag, verbose);
}

static int default_linear_path_nbits = 2048;

void
set_default_linear_path_nbits (const int s)
{
  default_linear_path_nbits = s;
}

Fixed_Width_Fingerprint_Bits::Fixed_Width_Fingerprint_Bits () : _nbits(default_linear_path_nbits)
{
  _bits = new fwf_bits_t[_nbits];

  std::fill_n(_bits, _nbits, 0);

  return;
}

Fixed_Width_Fingerprint_Bits::Fixed_Width_Fingerprint_Bits(int s) : _nbits(s)
{
  _bits = new fwf_bits_t[_nbits];

  std::fill_n(_bits, _nbits, 0);

  return;
}

Fixed_Width_Fingerprint_Bits::~Fixed_Width_Fingerprint_Bits()
{
  delete [] _bits;

  return;
}

void
Fixed_Width_Fingerprint_Bits::clear ()
{
  std::fill_n(_bits, _nbits, 0);

  return;
}

int
Fixed_Width_Fingerprint_Bits::identify_differences (const Fixed_Width_Fingerprint_Bits & rhs,
                                                    std::ostream & output) const
{
  assert(_nbits == rhs._nbits);

  int rc = 0;    

  for (int i = 0; i < _nbits; ++i)
  {
    if (_bits[i] == rhs._bits[i])
      continue;

    output << "bit mismatch " << i << " values " << _bits[i] << " and " << rhs._bits[i] << endl;
    rc++;
  }

  return rc;
}

fwf_bits_t & 
Fixed_Width_Fingerprint_Bits::operator[](const unsigned int s)
{
//cerr << "Fixed_Width_Fingerprint_Bits::getting " << s << endl;
  return _bits[s % _nbits];
}

int
Sparse_Fingerprint_Bits::identify_differences (const Sparse_Fingerprint_Bits & rhs,
                                                  std::ostream & output) const
{
  int rc = 0;
  for (auto i : *this)
  {
    const auto f = rhs.find(i.first);
    if (f == rhs.end())
    {
      output << "Sparse_Fingerprint_Bits::identify_differences:bit " << i.first << " missing from rhs\n";
      rc++;
    }

    if (i.second == f->second)
      continue;

    output << "Sparse_Fingerprint_Bits::identify_differences:bit count mismatch. Bit " << i.first << " values " << i.second << ' ' << f->second << '\n';
    rc++;
  }

  if (size() == rhs.size())
    return rc;

  output << "Sparse_Fingerprint_Bits::identify_differences:warning different sizes " << size() << " vs " << rhs.size() << '\n';
  rc++;

  return rc;
}

int
Fixed_Width_Fingerprint_Bits::append_fingerprint (const IWString & tag, IWString & output) const
{
  output << tag;
  IW_Bits_Base b;
  b.construct_from_array_of_ints(_bits, _nbits);
  IWString tmp;
  b.daylight_ascii_representation_including_nset_info(tmp);
  output << tmp << ">\n";
  return 1;
}

int
Sparse_Fingerprint_Bits::append_fingerprint (const IWString & tag, IWString & output) const
{
  Sparse_Fingerprint_Creator sfp;

  for (auto i : *this)
  {
    sfp.hit_bit(i.first, i.second);
//  cerr << "bit " << i.first << ' ' << i.second << " times\n";
  }

  output << tag;

  IWString tmp;
  sfp.daylight_ascii_form_with_counts_encoded(tmp);
//cerr << "Encoding " << tmp.length() << " " << size() << " bits\n";

  output << tmp << ">\n";

  return 1;
}

template <typename T>
int
Fixed_Width_Fingerprint_Bits::write_count(const IWDigits& iwdigits, T &output) const
{
  int nset = 0;

  for (int i = 0; i < _nbits; ++i)
  {
    iwdigits.append_number(output, _bits[i]);
    if (_bits[i])
      nset++;
  }

  iwdigits.append_number(output, nset);

  return 1;
}

template <typename T>
int
Sparse_Fingerprint_Bits::write_count(const IWDigits& iwdigits, T &output) const
{
  iwdigits.append_number(output, nset());

  return 1;
}

template <typename T>
int
Fixed_Width_Fingerprint_Bits::write_as_zero_and_one (T & output) const
{
  int nset = 0;

  for (int i = 0; i < _nbits; ++i)
  {
    if (_bits[i])
    {
      output << " 1";
      nset++;
    }
    else
      output << " 0";
  }

  output << ' ' << nset << '\n';

  return 1;
}

template <typename T>
int
Sparse_Fingerprint_Bits::write_as_zero_and_one (T & output) const
{
  cerr << "Sparse_Fingerprint_Bits:write_as_zero_and_one:not implemented\n";
  return 0;
}

int
Fixed_Width_Fingerprint_Bits::nset() const
{
  int rc = 0;
  for (int i = 0; i < _nbits; ++i)
  {
    if (_bits[i] > 0)
      rc++;
  }

  return rc;
}

int
Fixed_Width_Fingerprint_Bits::truncate_to_max_hits(const int x)
{
  int rc = 0;

  for (int i = 0; i < _nbits; ++i)
  {
    if (_bits[i] <= static_cast<fwf_bits_t>(x))
      continue;

    _bits[i] = x;
    rc++;
  }

  return rc;
}

int 
Sparse_Fingerprint_Bits::truncate_to_max_hits (const int x)
{
  int rc = 0;

  for (auto b : *this)
  {
    if (b.second <= x)
      continue;

    b.second = x;
    rc++;
  }

  return rc;
}

int
Sparse_Fingerprint_Bits::write_as_md5_sum (IWString & output) const
{
  return unordered_map_to_md5(*this, output);
}

template <typename T>
int
Sparse_Fingerprint_Bits::write_as_feature_count(const IWString & sep, T & output) const
{
  const auto sz = size();

//cerr << "Sparse_Fingerprint_Bits::write_as_feature_count:have " << sz << " bits set\n";

  if (0 == sz)
    return 0;

  sfp_bits_type * u = new sfp_bits_type[sz]; std::unique_ptr<sfp_bits_type[]> free_u(u);

  unsigned int ndx = 0;

  for (auto f : *this)
  {
    u[ndx] = f.first;
    ndx++;
  }

  assert (ndx == sz);

  std::sort(u, u + ndx);

  for (unsigned int i = 0; i < ndx; ++i)
  {
    if (i > 0)
      output << sep;

    const auto f = find(u[i]);

    output << u[i] << ':' << f->second;
  }

  return 1;
}

int
Fixed_Width_Fingerprint_Bits::write_as_md5_sum(IWString & output) const
{
  fwf_bits_t * u = new fwf_bits_t[_nbits]; std::unique_ptr<fwf_bits_t[]> free_u(u);

  int ndx = 0;
  for (int i = 0; i < _nbits; ++i)
  {
    if (0 == _bits[i])
      continue;

    u[ndx] = _bits[i];
    ndx++;
  }

  if (0 == ndx)
    return 0;

  MD5_CTX ctx;
  MD5Init(&ctx);
  MD5Update(&ctx, reinterpret_cast<unsigned char *>(u), ndx * sizeof(fwf_bits_t));

  unsigned char digest[16];
  MD5Final(digest, &ctx);

  output.append_hex(digest, 16);

  return 1;
}

template <typename T>
int
Fixed_Width_Fingerprint_Bits::write_as_feature_count(const IWString & sep, T & output) const
{
  int first = 1;

  for (int i = 0; i < _nbits; ++i)
  {
    if (0 == _bits[i])
      continue;

    if (first)
      first = 0;
    else
      output << sep;

    output << i << ':' << _bits[i];
  }

  return 1;
}

//int
//Sparse_Fingerprint_Bits::write_as_feature_count(IWString const & sep, IWString_and_File_Descriptor& output) const
//{
//}

template int Fixed_Width_Fingerprint_Bits::write_as_feature_count(const IWString & sep, IWString_and_File_Descriptor & output) const;
template int Fixed_Width_Fingerprint_Bits::write_count<IWString>(const IWDigits&, IWString&) const;
template int Sparse_Fingerprint_Bits::write_as_zero_and_one(IWString &) const;
template int Fixed_Width_Fingerprint_Bits::write_as_zero_and_one<IWString>(IWString&) const;
//template int LFP::Linear_Fingerprint_Creator::process<LFP::Sparse_Fingerprint_Bits>(Molecule&, int const*, LFP::Sparse_Fingerprint_Bits&);
//template int LFP::Linear_Fingerprint_Creator::process<LFP::Sparse_Fingerprint_Bits>(Molecule&, int const*, int const *, LFP::Sparse_Fingerprint_Bits&);

template int LFP::Sparse_Fingerprint_Bits::write_count<IWString>(const IWDigits&, IWString&) const;
//template int LFP::Linear_Fingerprint_Creator::process(Molecule&, int const*, LFP::Fixed_Width_Fingerprint_Bits&);
//template int LFP::Linear_Fingerprint_Creator::process(Molecule&, int const*, int const *, LFP::Fixed_Width_Fingerprint_Bits&);
template int LFP::Sparse_Fingerprint_Bits::write_as_feature_count(const IWString & sep, IWString_and_File_Descriptor & output) const;

template int LFP::Linear_Fingerprint_Creator::process<unsigned int, LFP::Fixed_Width_Fingerprint_Bits>(Molecule&, unsigned int const*, int const*, LFP::Fixed_Width_Fingerprint_Bits&);
template int LFP::Linear_Fingerprint_Creator::process<unsigned int, LFP::Sparse_Fingerprint_Bits>(Molecule&, unsigned int const*, int const*, LFP::Sparse_Fingerprint_Bits&);
template int LFP::Linear_Fingerprint_Creator::process<int, LFP::Fixed_Width_Fingerprint_Bits>(Molecule&, int const*, LFP::Fixed_Width_Fingerprint_Bits&);

} // end namespace LFP
