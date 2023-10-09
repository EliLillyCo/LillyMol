#include <iostream>
#include <limits>
#include <ctype.h>

using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "molecule.h"
#include "iwmfingerprint.h"
#include "path.h"

static int bits_per_iwmfingerprint = 2048;

static IWDigits iwdigits(256);

/*
  Just need something with a constructor
*/

class NotUsed
{
  private:
  public:
    NotUsed();
};

NotUsed::NotUsed()
{
  iwdigits.set_include_leading_space(1);
}

static NotUsed notused;    // er, this really does not work, it is order dependent and should be fixed sometime

void
set_iwmfingerprint_nbits (int n)
{
  assert (n > 1);

  bits_per_iwmfingerprint = n;
  
  if (0 != n % 8)
    cerr << "\nWarning, you have requested " << n << " bits per fingerprint, but that is not divisible by 8. Things will fail!\n\n";

  return;
}

int
iwmfingerprint_nbits()
{
  return bits_per_iwmfingerprint;
}

static int form_bit_vector_during_fingerprint_creation = 1;

void
set_form_bit_vector_during_fingerprint_creation(int s)
{
  form_bit_vector_during_fingerprint_creation = s;
}

static int min_path_length = 0;
static int max_path_length = 7;

int
set_min_path_length (int l)
{
  assert (l >= 0 && l <= max_path_length);

  min_path_length = l;

  return 1;
}

int
set_max_path_length (int l)
{
  assert (l >= 0 && l >= min_path_length);

  max_path_length = l;

  return 1;
}

static int omit_ch2 = 0;

void 
set_omit_ch2 (int o)
{
  omit_ch2 = o;
}

/*
  There are some bits set that will not permit accurate perception of
  substructures. These can be enabled or disabled as needed.
*/

static int only_bits_preserving_substructure_perception = 0;

void 
set_only_bits_preserving_substructure_perception(int s)
{
  only_bits_preserving_substructure_perception = s;
}

/*
  set_formal_charge_bits is actually the longest path for which formal charge
  bits are applied
*/

static int set_formal_charge_bits = -1;

void
set_include_formal_charge_bits (int s)
{
  set_formal_charge_bits = s;
}

static int formally_charged_atoms_lose_identity = 0;

void
set_formally_charged_atoms_lose_identity (int s)
{
  formally_charged_atoms_lose_identity = s;
}

static int do_atom_pair_bits = 0;

void
set_do_atom_pair_bits (int s)
{
  do_atom_pair_bits = s;
}

static int isotopic_paths = 0;

void
set_isotopic_paths (int s)
{
  isotopic_paths = s;
}

static int unsaturated_includes_aromatic = 0;

void
set_unsaturated_includes_aromatic (int s)
{
  unsaturated_includes_aromatic = s;

  return;
}

static int maximal_daylight_compatibility = 0;

void
set_maximal_daylight_compatibility (int s)
{
  maximal_daylight_compatibility = s;
}

/*
  I saw cases where a hydroxy was found to be similar to an ether. That's
  great for substructure searching, but less so for similarity searching.
  Therefore, we can turn on extra bits if the ends of a path have hydrogen
  attachments
  All paths of length less than or equal to BITS_FOR_HYDROGEN_ATTACHMENTS 
  will have those bits turned on
*/

static int bits_for_hydrogen_attachments = -1;

void
set_bits_for_hydrogen_attachments (int s)
{
  bits_for_hydrogen_attachments = s;
}

/*
  Do we fingerprint an isotopic variant the same or differently

  We do extra isotopic bits for all paths of length <= include_isotopic_information,
  so
*/

static int max_path_length_isotopic_bits = -1;

void
set_max_path_length_isotopic_bits (int s)
{
  max_path_length_isotopic_bits = s;
}

static int place_bits_for_rings = 0;

void
set_include_bits_for_rings (int s)
{
  place_bits_for_rings = s;
}

static int bits_for_atoms_in_multiple_rings = 0;

void
set_bits_for_atoms_in_multiple_rings(int s)
{
  bits_for_atoms_in_multiple_rings = s;
}

static int include_unsaturation_in_atom_type = 0;

void
set_include_unsaturation_in_atom_type (int s)
{
  include_unsaturation_in_atom_type = s;
}

/*
  Saw cases where we have a functional group on the outside
  of a molecule but a neighbour gets matches with the group inside
  the molecule. Set extra bits if a path is terminal
*/

static int terminal_path_bits = -1;

void
set_terminal_path_bits (int s)
{
  terminal_path_bits = s;
}

static int differentiate_ring_and_chain_bonds = 0;

void
set_differentiate_ring_and_chain_bonds(const int s)
{
  differentiate_ring_and_chain_bonds = s;
}

void
IWMFingerprint::_default_values()
{
  _bvector = nullptr;
  _auxiliary_bvector = nullptr;

  _min_heteroatoms_at_path_ends = 0;

#ifdef COUNT_TIMES_ATOM_IN_PATH
  _cycles_adjust_times_used = 0;
#endif

  return;
}

IWMFingerprint::IWMFingerprint()
{
  _default_values();

  return;
}

IWMFingerprint::IWMFingerprint (int s) : IW_Bits_Base (s)
{
  _default_values();

  return;
}

IWMFingerprint::~IWMFingerprint()
{
  if (nullptr != _bvector)
    delete [] _bvector;
  if (nullptr != _auxiliary_bvector)
    delete [] _auxiliary_bvector;

  return;
}

void
IWMFingerprint::set_min_heteroatoms_at_path_ends (int m)
{
  assert (m >= 0 && m <= 2);

  _min_heteroatoms_at_path_ends = m;

  return;
}

int
IWMFingerprint::delete_vector_representation()
{
  if (nullptr != _bvector)
  {
    delete [] _bvector;
    _bvector = nullptr;
  }

  if (nullptr != _auxiliary_bvector)
  {
    delete [] _auxiliary_bvector;
    _auxiliary_bvector = nullptr;
  }
  
  return 1;
}

int
IWMFingerprint::write_as_zero_and_one (std::ostream & os) const
{
  IWString buffer;
  buffer.resize(bits_per_iwmfingerprint * 2);

  if (_bvector)
  {
    if (_bvector[0])
      buffer << '1';
    else
      buffer << '0';

    for (int i = 1; i < bits_per_iwmfingerprint; i++)
    {
      if (_bvector[i])
        buffer += " 1";
      else
        buffer += " 0";
    }
  }
  else
    (void) IW_Bits_Base::append_string_form(buffer);

  os << buffer;

  return os.good();
}

int
IWMFingerprint::write_as_zero_and_one (IWString & buffer) const
{
  if (_bvector)
  {
    if (_bvector[0])
      buffer << '1';
    else
      buffer << '0';

    for (int i = 1; i < bits_per_iwmfingerprint; i++)
    {
      if (_bvector[i])
        buffer += " 1";
      else
        buffer += " 0";
    }
  }
  else
    (void) IW_Bits_Base::append_string_form(buffer);

  return 1;
}

int
IWMFingerprint::write_count (std::ostream & os) const
{
  assert (nullptr != _bvector);

  IWString buffer;
  buffer.resize(bits_per_iwmfingerprint * 3);

  buffer += _bvector[0];

  for (int i = 1; i < bits_per_iwmfingerprint; i++)
  {
    iwdigits.append_number(buffer, _bvector[i]);
//  buffer << ' ' << _bvector[i];
  }

  os << buffer;

  return os.good();
}

int
IWMFingerprint::write_count (IWString & buffer) const
{
  assert (nullptr != _bvector);


  buffer += _bvector[0];

  for (int i = 1; i < bits_per_iwmfingerprint; i++)
  {
    iwdigits.append_number(buffer, _bvector[i]);
  }

  return 1;
}

int 
IWMFingerprint::truncate_to_max_hits (int m)
{
  int rc = 0;
  for (int i = 0; i < bits_per_iwmfingerprint; i++)
  {
    if (_bvector[i] > m)
    {
      _bvector[i] = m;
      rc++;
    }
  }

  return rc;
}

int
IWMFingerprint::nset() const
{
  assert (nullptr != _bvector);

  int rc = 0;

  for (int i = 0; i < bits_per_iwmfingerprint; i++)
  {
    if (0 != _bvector[i])
      rc++;
  }

  return rc;
}

/*
  This class is used as a helper
*/

//#define CHECK_COLLISIONS

#ifdef CHECK_COLLISIONS
static int collisions = 0;
#endif

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

#ifdef CHECK_COLLISIONS
    resizable_array<int> ** _first_path;
    void _report_collision (int zbit) const;
    void _check_collision_forward (int zbit, const int * bvector) const;
    void _check_collision_backward (int zbit, const int * bvector) const;
#endif

//  private functions

    int _compute_atom_hash_value (atom_number_t i) const;
    int _compute_atom_hash_value_non_periodic_table (atom_number_t, int) const;
    int _compute_atom_hash_value_include_unsaturation (atom_number_t) const;
    void _set_cluster_bit (int * bvector, atom_number_t centre,
                           int b1, atom_number_t a1, int b2, atom_number_t a2, int b3, atom_number_t a3);
    void _build_branched_paths (int * bvector, const int * include_these_atoms);
    int  _build (atom_number_t aprev, int * bvector, int * auxiliary_vector, const int * include_these_atoms);
    void _set_auxiliary_bits (int * bvector, unsigned int zbit1, unsigned int zbit2, int is_ring, int astart, int astop) const;
    void _set_bit_forward (int *, int *, int) const;
    void _set_bit_backward (int *, int *, int) const;
    void _set_bit (int *, int *, int) const;
    void _omit_ch2 (int *, int) const;
    void _omit_ch2 (int *, int, int) const;
    void _set_bit_omit (int *, int, int) const;
    void _do_formal_charge_bits (int, int, int, int *) const;
    void _do_formal_charge_bits (int *, int, int) const;
    void _do_isotope_bits (int istart, int istop, int istep, int * bvector) const;
    void _do_isotope_bits (int * bvector) const;
    void _do_ring_bits (int * bvector) const;
    void _set_terminal_path_bits (int * auxiliary_bvector) const;
    void _set_terminal_path_bits (int * auxiliary_bvector, int) const;
    void _set_bits_for_atoms_in_multiple_rings(int direction, int * bvector) const;

  public:
    MFingerprint (Molecule &, Atom_Typing_Specification &, const int * atype);
    ~MFingerprint();

    void set_min_heteroatoms_at_path_ends (int);

    int fingerprint_large_ring (const Ring &, int *);

    int build (atom_number_t, int *, int *, const int * inc);

#ifdef COUNT_TIMES_ATOM_IN_PATH
    void place_labels_for_times_in_path(Molecule & m) const;
    const int * times_in_path() const { return _times_in_path;}
    int generate_bits(const Set_of_Atoms & s, const resizable_array<int> & bonds, int * bvector, int * auxiliary_bvector);
#endif
};

/*
*/

MFingerprint::MFingerprint (Molecule & m,
                            Atom_Typing_Specification & ats,
                            const int * atype)
{
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

// Work hard to avoid aromaticity computation unless really needed

  _arom = new aromaticity_type_t[_matoms];

  if (ats.active() || nullptr != atype)
    set_vector(_arom, _matoms, 0);
  else
    m.aromaticity(_arom);

  _unsaturation = new int[_matoms];

  _atoms = new Atom *[_matoms];
  m.atoms((const Atom **) _atoms);    // the a->nbonds() method is non const, that's why we need a cast

  if (bits_for_hydrogen_attachments >= 0)
    _hcount = new int[_matoms];
  else
    _hcount = nullptr;

  if (max_path_length_isotopic_bits >= 0)
  {
    _isotope = new isotope_t[_matoms];
    m.get_isotopes(_isotope);
  }
  else
    _isotope = nullptr;

  _nrings = new int[_matoms];
  m.ring_membership(_nrings);

  if (nullptr != atype)
  {
    copy_vector(_atom_hash_value, atype, _matoms);
    set_vector(_unsaturation, _matoms, 0);
    if (nullptr != _hcount)
      set_vector(_hcount, _matoms, 0);
  }
  else if (ats.active())
  {
    ats.assign_atom_types(m, _atom_hash_value);
    set_vector(_unsaturation, _matoms, 0);
    if (nullptr != _hcount)
      set_vector(_hcount, _matoms, 0);
  }
  else
  {

// Note that (for efficiency) we don't consider aromatic atoms with a double bond outside the ring as being unsaturated

    for (int i = 0; i < _matoms; i++)
    {
      if (include_unsaturation_in_atom_type)
        _atom_hash_value[i] = _compute_atom_hash_value_include_unsaturation(i);
      else
        _atom_hash_value[i] = _compute_atom_hash_value(i);

      if (! unsaturated_includes_aromatic && _arom[i])
        _unsaturation[i] = 0;
      else
        _unsaturation[i] = _atoms[i]->nbonds() - _atoms[i]->ncon();

//   Don't set Hydrogen bits for aromatic Nitrogens. Too many problems with tautomers
//   C12=C(N=CN=C1O)NN=C2 
//   C12=CNN=C1N=CN=C2O 

      if (bits_for_hydrogen_attachments >= 0)
      {
        if (_arom[i] && 7 == _atomic_number[i])   // always 0
          _hcount[i] = 0;
        else
          _hcount[i] = m.hcount(i);
      }
    }
  }

  _molecule_contains_fused_rings = 0;

  if (bits_for_atoms_in_multiple_rings)
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

#ifdef CHECK_COLLISIONS

void
MFingerprint::_report_collision (int zbit) const
{
  collisions++;

  cerr << "Collision on bit " << zbit << endl;

  const resizable_array<int> * fp = _first_path[zbit];

  cerr << "First path (" << (fp->number_elements() / 2) << "):";
  for (int i = 0; i < fp->number_elements(); i++)
  {
    cerr << ' ' << fp->item(i);
  }

  cerr << endl;

  cerr << "New path (" << _path_length << "):";
  for (int i = 0; i <= _path_length; i++)
  {
    cerr << ' ' << _path_hash_value[i];
    if (i < _path_length)
      cerr << ' ' << _bond[i];
  }

  cerr << endl;

  return;
}

void
MFingerprint::_check_collision_forward (int zbit, const int * bvector) const
{
  if (1 == bvector[zbit])
  {
    assert (nullptr == _first_path[zbit]);
    _first_path[zbit] = new resizable_array<int>;
    _first_path[zbit]->resize(2 * _path_length);
    for (int i = 0; i <= _path_length; i++)
    {
      _first_path[zbit]->add(_path_hash_value[i]);
      if (i < _path_length)
        _first_path[zbit]->add(_bond[i]);
    }

    return;
  }

// We already have a path which hit this bit. Is this a duplicate?

  const resizable_array<int> * fp = _first_path[zbit];
  assert (nullptr != fp);

  int nfp = fp->number_elements();

  if (2 * _path_length + 1 != nfp)
  {
    _report_collision(zbit);
    return;
  }

// Lengths are the same, are they the same path?

  int j = 0;
  for (int i = 0; i <= _path_length; i++)
  {
    if (_path_hash_value[i] != fp->item(j))
    {
      _report_collision(zbit);
      return;
    }
    j++;
    if (i == _path_length)
      break;

    if (_bond[i] != fp->item(j))
    {
      _report_collision(zbit);
      return;
    }

    j++;
  }

  return;
}

void
MFingerprint::_check_collision_backward (int zbit, const int * bvector) const
{
  if (1 == bvector[zbit])
  {
    assert (nullptr == _first_path[zbit]);
    _first_path[zbit] = new resizable_array<int>;
    _first_path[zbit]->resize(2 * _path_length);
    _first_path[zbit]->add(_path_hash_value[_path_length]);
    for (int i = _path_length - 1; i >= 0; i--)
    {
      _first_path[zbit]->add(_bond[i]);
      _first_path[zbit]->add(_path_hash_value[i]);
    }

    return;
  }

// We already have a path which hit this bit. Is this a duplicate?

  const resizable_array<int> * fp = _first_path[zbit];
  assert (nullptr != fp);

  int nfp = fp->number_elements();

  if (2 * _path_length + 1 != nfp)
  {
    cerr << "Length collision\n";
    _report_collision(zbit);
    return;
  }

// Lengths are the same, are they the same path?

  if (_path_hash_value[_path_length] != fp->item(0))
  {
    _report_collision(zbit);
    return;
  }

  int j = 1;     // which item in FP to compare
  for (int i = _path_length - 1; i >= 0; i--)
  {
    if (_bond[i] != fp->item(j))
    {
      _report_collision(zbit);
      return;
    }

    j++;
    if (_path_hash_value[i] != fp->item(j))
    {
      cerr << "Backward collosion at i = " << i << " j = " << j << endl;
      _report_collision(zbit);
      return;
    }
    j++;
  }

  return;
}

#endif

//#define DEBUG_DO_FORMAL_CHARGE_BITS

/*
  Just like the regular path bit setting, except that when dealing with
  either of the end atoms, if they are charged, we key on their formal
  charge (if present);
*/

void
MFingerprint::_do_formal_charge_bits (int istart, int istop, int istep,
                                      int * bvector) const
{
  atom_number_t astart = _path[istart];
  atom_number_t astop  = _path[istop];

  unsigned int zbit = 931 + 7823 * _atoms[astart]->formal_charge() + 32771 * _atoms[astop]->formal_charge();

  istart += istep;    // we have already accounted for the first atom

  for (int i = istart; i != istop; i += istep)
  {
    zbit = 11517 * zbit + _path_hash_value[i] + _bond[i - 1];
  }

  bvector[zbit % bits_per_iwmfingerprint]++;

#ifdef DEBUG_DO_FORMAL_CHARGE_BITS
  cerr << "Formal charge bit " << (zbit % bits_per_iwmfingerprint) << endl;
#endif

  return;
}

/*
  We need to ensure that formal charge paths are treated in a canonical form
*/

void
MFingerprint::_do_formal_charge_bits (int * bvector,
                                      atom_number_t astart,
                                      atom_number_t astop) const
{
//cerr << " astart " << astart << " path[0] " << _path[0] << ", path length " << _path_length << endl;
//cerr << " astop " << astop << " path[_path_length] " << _path[_path_length] << endl;
  assert ((astart == _path[0] && astop == _path[_path_length]) || (astart == _path[_path_length] && astop == _path[0]));

  formal_charge_t fc0 = _atoms[astart]->formal_charge();

  if (0 == _path_length)
  {
    unsigned int b;

    if (fc0 >= 0)
      b = (767772 * fc0) % bits_per_iwmfingerprint;
    else
      b = (-21118 * fc0) % bits_per_iwmfingerprint;

#ifdef DEBUG_DO_FORMAL_CHARGE_BITS
    cerr << "Zero path length formal charge bit " << b << endl;
#endif

    bvector[b]++;

    return;
  }

  formal_charge_t fc1 = _atoms[astop]->formal_charge();

  if (fc0 > fc1)
    _do_formal_charge_bits(0, _path_length, 1, bvector);
  else if (fc0 < fc1)
    _do_formal_charge_bits(_path_length, 0, -1, bvector);
  else if (_path_hash_value[0] > _path_hash_value[_path_length])
    _do_formal_charge_bits(0, _path_length, 1, bvector);
  else if (_path_hash_value[0] < _path_hash_value[_path_length])
    _do_formal_charge_bits(_path_length, 0, -1, bvector);
  else
    _do_formal_charge_bits(0, _path_length, 1, bvector);

  return;
}

//#define DEBUG_DO_ISOTOPE_BITS

void 
MFingerprint::_do_isotope_bits (int istart,
                                int istop,
                                int istep,
                                int * bvector) const
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

  bvector[zbit % bits_per_iwmfingerprint]++;
}

void
MFingerprint::_do_isotope_bits (int * bvector) const
{
  assert (nullptr != _isotope);

  atom_number_t astart = _path[0];

  // Deliberately keep i0 and i1 as int's for backwards compatability.
  int i0 = _isotope[astart];
  if (0 == _path_length)
  {
    unsigned int b = 77265 + _path_hash_value[0] * 87 + i0 % 87;

    bvector[b % bits_per_iwmfingerprint]++;

#ifdef DEBUG_DO_ISOTOPE_BITS
  cerr << "MFingerprint::_do_isotope_bits:setting " << (b % bits_per_iwmfingerprint) << endl;
#endif

    return;
  }

  atom_number_t astop = _path[_path_length];

  int i1 = _isotope[astop];

  if (i0 > i1)
    _do_isotope_bits(0, _path_length, 1, bvector);
  else if (i0 < i1)
    _do_isotope_bits(_path_length, 0, -1, bvector);
  else if (0 == i0 && 0 == i1)
    ;
  else if (_path_hash_value[0] < _path_hash_value[_path_length])
    _do_isotope_bits(0, _path_length, 1, bvector);
  else
    _do_isotope_bits(_path_length, 0, -1, bvector);

  return;
}

//#define DEBUG_SET_CLUSTER_BIT

void
MFingerprint::_set_cluster_bit (int * bvector, atom_number_t centre,
                                int b1, atom_number_t a1,
                                int b2, atom_number_t a2,
                                int b3, atom_number_t a3)
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

  bvector[zbit % bits_per_iwmfingerprint]++;

#ifdef DEBUG_SET_CLUSTER_BIT
  cerr << "Cluster bit " << (zbit % bits_per_iwmfingerprint) << endl;
#endif

  return;
}

static inline int
numeric_bond_code (const Bond * b)
{
  if (b->is_aromatic())
    return 1;

  int rc;
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
  
  if (differentiate_ring_and_chain_bonds && b->nrings())
    rc += 18;

  return rc;
}

void
MFingerprint::_build_branched_paths (int * bvector,
                                     const int * include_these_atoms)
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

    if (nullptr != include_these_atoms && 0 == include_these_atoms[j])
      continue;

    _set_cluster_bit(bvector, a1, _bond[0], a0, _bond[1], a2, numeric_bond_code(b), j);
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

void
MFingerprint::_set_auxiliary_bits (int * bvector,
                                   unsigned int zbit1, unsigned int zbit2,
                                   int is_ring,
                                   int astart, int astop) const
{
#ifdef DEBUG_SET_AUX_BITS
  cerr << "Setting auxiliary bits\n";
#endif

  if (_path_length <= bits_for_hydrogen_attachments)
  {
    if (_hcount[astart])
    {
//    unsigned int hbit = ((zbit1 << 6) + zbit2 + 3131) % bits_per_iwmfingerprint;   what I intended
      unsigned int hbit = (zbit1 << (6 + zbit2 + 3131)) % bits_per_iwmfingerprint;   // what the language specifies
      INCREMENT_VECTOR(bvector, hbit, "hcount1");
    }
    if (_hcount[astop])
    {
//    unsigned int hbit = ((zbit2 << 6) + zbit1 + 8181) % bits_per_iwmfingerprint;   what I intended
      unsigned int hbit = (zbit2 << (6 + zbit1 + 8181)) % bits_per_iwmfingerprint;   // what the language specifies
      INCREMENT_VECTOR(bvector, hbit, "hcount2");
    }
  }

  if (is_ring)
  {
    unsigned int rbit1 = (zbit1 * 93717 * _path_length) % bits_per_iwmfingerprint;
    INCREMENT_VECTOR(bvector, rbit1, "ring");
  }
//else if (_nrings[astart] || _nrings[astop])
  else if (1 == _path_length && (_nrings[astart] || _nrings[astop]))
  {
    unsigned int rbit = (zbit1 + zbit2 + 119) % bits_per_iwmfingerprint;
    INCREMENT_VECTOR(bvector, rbit, "ring end");

//  if (_nrings[astart] > 1 || _nrings[astop] > 1)    // could never get this to work
//  {
//    rbit = (rbit * 31) % bits_per_iwmfingerprint;
//    bvector[rbit]++;
//  }
  }

// do something for branching. Make sure you keep the numbers the same in
// both places. Note that since these depend on ncon, they break the
// substructure relation

  if (only_bits_preserving_substructure_perception)
    ;
  else if (_ncon[astart] > 2 && 6 != _atomic_number[astart])
  {
    unsigned int cbit = (zbit1 + bits_per_iwmfingerprint) >> 1;

    INCREMENT_VECTOR(bvector, cbit, "branch0");

    if (_ncon[astart] > 3)
      INCREMENT_VECTOR(bvector, (((cbit + zbit2) & 8715) % bits_per_iwmfingerprint), "ncon3");    // parens added Jan 2017

    if (0 == _arom[astart] && _unsaturation[astart])
      bvector[(cbit + 31290) % bits_per_iwmfingerprint]++;
  }

  if (only_bits_preserving_substructure_perception)
    ;
  else if (_path_length > 1 && _ncon[astop] > 2 && 6 != _atomic_number[astop])
  {
    unsigned int cbit = (zbit1 + bits_per_iwmfingerprint) >> 1;

    INCREMENT_VECTOR(bvector, cbit, "branch0");

    if (_ncon[astop] > 3)
      INCREMENT_VECTOR(bvector, ((cbit + zbit1) & 8715) % bits_per_iwmfingerprint, "ncon3");    // parens added Jan 2017

    if (0 == _arom[astop] && _unsaturation[astop])
    INCREMENT_VECTOR(bvector, (cbit + 31290) % bits_per_iwmfingerprint, "arom");
  }

// and unsaturation

  if (0 == _arom[astart] && _unsaturation[astart])
  {
    unsigned int ubit = (zbit2 * 13) % bits_per_iwmfingerprint;
    INCREMENT_VECTOR(bvector, ubit, "unsaturation");
  }

  if (_path_length > 0 && 0 == _arom[astop] && _unsaturation[astop])
  {
    unsigned int ubit = (zbit2 * 13) % bits_per_iwmfingerprint;
    INCREMENT_VECTOR(bvector, ubit, "unsaturation");
  }

  if (_path_length > 2 && (6 != _atomic_number[astart] || 6 != _atomic_number[astop]))
  {
    int unsat = 0;
    int branches = 0;

    for (int i = 1; i < _path_length; i++)
    {
      atom_number_t a = _path[i];
      branches += _ncon[a] - 2;
      if (0 == _arom[a])
        unsat += _unsaturation[a];
    }

    if (unsat)
      INCREMENT_VECTOR(bvector, (unsat + 387712 + bits_per_iwmfingerprint) % bits_per_iwmfingerprint, "unsat");
//    bvector[(unsat + 387712 + bits_per_iwmfingerprint) % bits_per_iwmfingerprint]++;

    if (only_bits_preserving_substructure_perception)
      ;
    else if (branches)
      INCREMENT_VECTOR(bvector, (_path_length * branches + zbit2%6) % bits_per_iwmfingerprint, "branches");
  }

  if (_path_length > set_formal_charge_bits)   // it gets initialised to -1
    ;
  else if (_atoms[astart]->formal_charge() || _atoms[astop]->formal_charge())
    _do_formal_charge_bits(bvector, astart, astop);

  if (_path_length > max_path_length_isotopic_bits)   // it gets initialised to -1
    ;
  else if (_isotope[astart] || _isotope[astop])
    _do_isotope_bits(bvector);

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

void
MFingerprint::_set_bit_forward (int * bvector, int * auxiliary_bvector, int is_ring) const
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

  zbit1 = zbit1 % bits_per_iwmfingerprint;
  zbit2 = zbit2 % bits_per_iwmfingerprint;

  bvector[zbit1]++;
  if (_path_length > 2 && zbit2 != zbit1)
    bvector[zbit2]++;

#ifdef DEBUG_SET_BIT
  cerr << "path of length " << _path_length << " sets bits " << zbit1 << " and " << zbit2 << endl;
#endif

#ifdef CHECK_COLLISIONS
  _check_collision_forward(zbit1, bvector);
  if (zbit1 != zbit2)
    _check_collision_forward(zbit2, bvector);
#endif

  if (! maximal_daylight_compatibility)
    _set_auxiliary_bits(auxiliary_bvector, zbit1, zbit2, is_ring, _path[0], _path[_path_length]);

  return;
}

void
MFingerprint::_set_bit_backward (int * bvector, int * auxiliary_bvector, int is_ring) const
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

  zbit1 = zbit1 % bits_per_iwmfingerprint;
  zbit2 = zbit2 % bits_per_iwmfingerprint;

  bvector[zbit1]++;
  if (_path_length > 2 && zbit2 != zbit1)
    bvector[zbit2]++;

#ifdef DEBUG_SET_BIT
  cerr << "path of length " << _path_length << " sets bits " << zbit1 << " and " << zbit2 << endl;
#endif

#ifdef CHECK_COLLISIONS
  _check_collision_backward(zbit1, bvector);
  if (zbit1 != zbit2)
    _check_collision_backward(zbit2, bvector);
#endif

  if (! maximal_daylight_compatibility)
    _set_auxiliary_bits(auxiliary_bvector, zbit1, zbit2, is_ring, _path[_path_length], _path[0]);

  return;
}

/*
  A complete path is in the object. Set the appropriate bit
  The first thing to do is decide on the direction in order to
  get a canonical path
*/

void
MFingerprint::_set_bit(int * bvector, int * auxiliary_bvector, int astart) const
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

  if (isotopic_paths && ! _atoms[_path[_path_length]]->is_isotope())
    return;

  if (0 == _path_length)
  {
    _set_bit_forward(bvector, auxiliary_bvector, 0);      // last arg 0 means not a ring path
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
    _set_bit_forward(bvector, auxiliary_bvector, is_ring);
    if (bits_for_atoms_in_multiple_rings && _molecule_contains_fused_rings && _path_length > 2)
      _set_bits_for_atoms_in_multiple_rings(1, auxiliary_bvector);
  }
  else
  {
    _set_bit_backward(bvector, auxiliary_bvector, is_ring);
    if (bits_for_atoms_in_multiple_rings && _molecule_contains_fused_rings && _path_length > 2)
      _set_bits_for_atoms_in_multiple_rings(-1, auxiliary_bvector);
  }

// We only do the trick of leaving out CH2 groups when the path is length 3 or longer

  if (_path_length < 3 || 0 == omit_ch2)
    return;

// If omit_ch2 is just 1, we only omit the CH2 groups if there are heteroatoms at the ends

  if (omit_ch2 > 1)     // do it regardless
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
    _omit_ch2(bvector, first_ch2, last_ch2);

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

void
MFingerprint::_omit_ch2 (int * bvector, int first_ch2, int last_ch2) const
{
  int d1 = first_ch2;
  int d2 = _path_length - last_ch2;

  if (d1 < d2)
    _omit_ch2(bvector, first_ch2);
  else if (d1 > d2)
    _omit_ch2(bvector, last_ch2);
  else
  {
    _omit_ch2(bvector, first_ch2);
    _omit_ch2(bvector, last_ch2);
  }

  return;
}

/*
  We are to set some bits based on a path which omits OMIT.
  First decide on a direction
*/

void
MFingerprint::_omit_ch2 (int * bvector, int omit) const
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

  _set_bit_omit(bvector, direction, omit);

  return;
}

/*
  Set a bit based on a path from which a member of the path is omitted.

  We need to keep this in-synch with _set_bit_forward and _set_bit_backward
  in order to set the same bit
*/

void
MFingerprint::_set_bit_omit (int * bvector, int direction, int omit) const
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

  bvector[zbit1 % bits_per_iwmfingerprint]++;
  bvector[zbit2 % bits_per_iwmfingerprint]++;

#ifdef DEBUG_SET_BIT_OMIT
  cerr << "By omitting " << omit << " bits " << (zbit1 % bits_per_iwmfingerprint) << " and " << (zbit2 % bits_per_iwmfingerprint) << endl;
#endif

  return;
}

void
MFingerprint::_set_bits_for_atoms_in_multiple_rings(int direction,
                                        int * bvector) const
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

  int atoms_in_multiple_rings = 0;

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
//  cerr << "Path of length " << _path_length << " setting bit " << (b % bits_per_iwmfingerprint) << endl;
  if (atoms_in_multiple_rings)
    bvector[(b + 97 * atoms_in_multiple_rings * _path_length) % bits_per_iwmfingerprint]++;

  return;
}

/*
  Build all fingerprints emanating from the current start atom
*/

int
MFingerprint::_build(atom_number_t aprev, int * bvector, int * auxiliary_bvector,
                     const int * include_these_atoms)
{
  assert (_path_length >= 0 && _path_length <= max_path_length);

  atom_number_t astart = _path[_path_length];

  assert (astart >= 0 && astart < _matoms);     // must be a valid atom for our molecule

  if ((_path_length <= set_formal_charge_bits || formally_charged_atoms_lose_identity) &&
      _atoms[astart]->formal_charge())
  {
    _path_hash_value[_path_length] = 7900 + 977 * _atoms[astart]->formal_charge();
  }
  else
    _path_hash_value[_path_length] = _atom_hash_value[astart];

  Atom * a = _atoms[astart];

  int acon = a->ncon();

  _ch2[_path_length] = (6 == _atomic_number[astart] && 2 == acon && 2 == a->nbonds());

  if (_path_length >= min_path_length && _path_length <= max_path_length)   // in range, write the path
    _set_bit(bvector, auxiliary_bvector, astart);

  if (2 == _path_length && min_path_length <= 2)
    _build_branched_paths(bvector, include_these_atoms);

  if (_in_path[astart])     // we have formed a ring, don't continue from here
  {
    if (place_bits_for_rings && astart == _path[0])
      _do_ring_bits(auxiliary_bvector);
    return 1;
  }

  if (_path_length <= terminal_path_bits)
    _set_terminal_path_bits(auxiliary_bvector);

  if (_path_length == max_path_length)    // we are at the longest path allowed, cannot extend it
    return 1;

// Looks like this path is continuing...

  _in_path[astart] = 1;

  _path_length++;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(astart);

    if (aprev == j)      // cannot turn back on ourselves
      continue;

    if (nullptr != include_these_atoms && 0 == include_these_atoms[j])
      continue;

    _bond[_path_length - 1] = numeric_bond_code(b);
    _path[_path_length] = j;

    _build(astart, bvector, auxiliary_bvector, include_these_atoms);
  }

  _path_length--;
  _in_path[astart] = 0;

  return 1;
}

int
MFingerprint::build (atom_number_t astart, int * bvector, int * auxiliary_bvector,
                     const int * include_these_atoms)
{
  if (2 == _min_heteroatoms_at_ends && 6 == _atomic_number[astart])
    return 1;

  if (isotopic_paths && ! _atoms[astart]->is_isotope())
    return 1;

  _path_length = 0;

  set_vector(_in_path, _matoms, 0);

  if (nullptr != include_these_atoms)
  {
    if (! include_these_atoms[astart])
    {
      cerr << "MFingerprint::build:start atom not included\n";
      abort();
    }
  }

  _path[0] = astart;

  return _build(INVALID_ATOM_NUMBER, bvector, auxiliary_bvector, include_these_atoms);
}

int
MFingerprint::fingerprint_large_ring (const Ring & r,
                                      int * bvector)
{
  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    _path[i] = r[i];
  }

  _path_length = ring_size;

  _do_ring_bits(bvector);

  return 1;
}

/*
  We have formed a ring
*/

void
MFingerprint::_do_ring_bits (int * bvector) const
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

  for (int i = 0; i < place_bits_for_rings; i++)
  {
    unsigned int b = ((i + 1) * (_path_length) * 327) % bits_per_iwmfingerprint;

#ifdef DEBUG_DO_RING_BITS
    cerr << "ring_size " << ring_size << " setting bit " << b << " currently " << bvector[b] << endl;
#endif

    bvector[b]++;
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
  cerr << " b = " << b << " -> " << (b % bits_per_iwmfingerprint) << endl;
#endif

  bvector[b % bits_per_iwmfingerprint]++;

  return;
}

/*
  We make a simplifying assumption by only fingerprinting terminal
  paths where one end is terminal. Saves an expensive canonicalisation
*/

void
MFingerprint::_set_terminal_path_bits (int * bvector)const
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
    _set_terminal_path_bits(bvector, 1);
  else if (0 == t0 && t1)
    _set_terminal_path_bits(bvector, -1);

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

void
MFingerprint::_set_terminal_path_bits (int * bvector,
                                       int direction) const
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

  bvector[(3 * b + _path_length + 1) % bits_per_iwmfingerprint]++;

  return;
}

#ifdef COUNT_TIMES_ATOM_IN_PATH
void
MFingerprint::place_labels_for_times_in_path(Molecule & m) const
{
  for (auto i = 0; i < _matoms; ++i)
  {
//  m.set_isotope(i, _times_in_path[i]);
    m.set_atom_map_number(i, _times_in_path[i]);
  }

  return;
}

int
MFingerprint::generate_bits(const Set_of_Atoms & s,
                            const resizable_array<int> & bonds,
                            int * bvector,
                            int * auxiliary_bvector)
{
  const int n = s.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = s[i];

    _path[i] = j;
    _path_hash_value[i] = _atom_hash_value[j];
    _in_path[j] = 1;
  }

  assert (n == bonds.number_elements() + 1);

  _path_length = n - 1;

  for (int i = 0; i < _path_length; ++i)
  {
    _bond[i] = bonds[i];
  }

  _set_bit(bvector, auxiliary_bvector, INVALID_ATOM_NUMBER);   // this last arg could be problematic, watch...

  return 1;
}

#endif

//#define CHECK_AUXILIARY_COLLISIONS

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

/*
  Finally, the function which does all the work
*/

int
IWMFingerprint::construct_fingerprint (Molecule & m,
                                       Atom_Typing_Specification & ats,
                                       const int * inc)
{
  return _construct_fingerprint(m, ats, nullptr, inc);
}

int 
IWMFingerprint::construct_fingerprint (Molecule & m,
                                       const int * atype,
                                       const int * inc)
{
  return _construct_fingerprint(m, _atom_typing_specification, atype, inc);
}

int 
IWMFingerprint::construct_fingerprint (Molecule & m,
                                       const int * atype)
{
  return _construct_fingerprint(m, _atom_typing_specification, atype, nullptr);
}

int
IWMFingerprint::construct_fingerprint (Molecule & m)
{
  return _construct_fingerprint(m, _atom_typing_specification, nullptr, nullptr);
}

#ifdef COUNT_TIMES_ATOM_IN_PATH
int
IWMFingerprint::construct_fingerprint_label_by_times_in_path (Molecule & m)
{
  int rc = _construct_fingerprint(m, _atom_typing_specification, nullptr);


  return rc;
}
#endif

int
IWMFingerprint::_construct_fingerprint (Molecule & m,
                                        Atom_Typing_Specification & ats,
                                        const int * atype,
                                        const int * include_these_atoms)
{
  assert (nullptr == _bvector);
  assert (0 == _nbits);

  allocate_space_for_bits(bits_per_iwmfingerprint);

  _bvector = new_int(bits_per_iwmfingerprint);

  _auxiliary_bvector = new_int(bits_per_iwmfingerprint);

  MFingerprint mfp(m, ats, atype);

  mfp.set_min_heteroatoms_at_path_ends(_min_heteroatoms_at_path_ends);

  const int matoms = m.natoms();

  if (0 == matoms)
  {
    IW_Bits_Base::clear();
    return 1;
  }

// If we are fingerprinting rings, capture any that are longer than our max path length

  if (place_bits_for_rings)
  {
    const int nr = m.nrings();
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);

      if (ri->number_elements() >= max_path_length)
        mfp.fingerprint_large_ring(*ri, _bvector);
    }
  }

  if (nullptr == include_these_atoms)
  {
    for (int i = 0; i < matoms; i++)    // do all atoms
    {
      (void) mfp.build(i, _bvector, _auxiliary_bvector, include_these_atoms);
    }
  }
  else
  {
    for (int i = 0; i < matoms; i++)
    {
      if (include_these_atoms[i])
        (void) mfp.build(i, _bvector, _auxiliary_bvector, include_these_atoms);
    }
  }

#ifdef COUNT_TIMES_ATOM_IN_PATH
//cerr << _cycles_adjust_times_used << " _cycles_adjust_times_used\n";
  if (_cycles_adjust_times_used)
    _extra_bits_to_equalise_times_in_path_4(m, mfp, include_these_atoms, _cycles_adjust_times_used);
//  _extra_bits_to_equalise_times_in_path_3(m, mfp, include_these_atoms, _cycles_adjust_times_used);

  mfp.place_labels_for_times_in_path(m);
#endif

#ifdef CHECK_AUXILIARY_COLLISIONS
  int collisions = 0;
#endif

  for (int i = 0; i < bits_per_iwmfingerprint; i++)
  {
#ifdef CHECK_AUXILIARY_COLLISIONS
    if (_bvector[i] && _auxiliary_bvector[i])
    {
      cerr << "Primary/secondary collision on bit " << i << " values " << _bvector[i] << " and " << _auxiliary_bvector[i] << endl;
      collisions++;
    }
#endif

    _bvector[i] += _auxiliary_bvector[i];
  }

#ifdef CHECK_AUXILIARY_COLLISIONS
  cerr << collisions << " primary/secondary collisions\n";
#endif

  if (form_bit_vector_during_fingerprint_creation)
    return IW_Bits_Base::construct_from_array_of_ints(_bvector, bits_per_iwmfingerprint);

  return 1;
}

int 
operator == (const IWMFingerprint & lhs, const IWMFingerprint & rhs)
{
  const IW_Bits_Base & lhs_ref = lhs;
  const IW_Bits_Base & rhs_ref = rhs;

  return lhs_ref == rhs_ref;
}

int
parse_path_length_options (Command_Line & cl, char cflag, int verbose)
{
  const_IWSubstring s;
  cl.value(cflag, s);

  if (s.length() < 3)
  {
    cerr << "The -s option must be off the form 'min/max'\n";
    return 0;
  }

  int i = s.index('/');
  if (i <= 0 || i == s.length() - 1)
  {
    cerr << "The -s option must be off the form 'min/max'\n";
    return 0;
  }

  int min_path_length = 0;
  while (isdigit(s[0]))
  {
    min_path_length = 10 * min_path_length + s[0] - '0';
    s++;
  }

  if ('/' != s[0])
  {
    cerr << "The -s option must be off the form 'min/max'\n";
    return 0;
  }

  s++;    // skip over the /

  int max_path_length;
  if (! s.numeric_value(max_path_length) || max_path_length < min_path_length)
  {
    cerr << "The max path length must be longer than the min path length\n";
    return 0;
  }

  set_max_path_length(max_path_length);     // set max first, otherwise entering -s 8/8 would fail
  set_min_path_length(min_path_length);

  if (verbose)
    cerr << "Will generate paths between " << min_path_length << " and " << max_path_length << endl;

  return 1;
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
  output << " -" << flag << " xrcb     differentiate ring and chain bonds\n";

  return output.good();
}

int
parse_misc_fingerprint_options (Command_Line & cl,
                                char flag,
                                int verbose)
{
  if (cl.option_present(flag))
  {
    int i = 0;
    const_IWSubstring y;
    while (cl.value(flag, y, i++))
    {
      const_IWSubstring directive;    // several of our things look like XX=nn
      int dvalue;
      int dvalue_valid;

      if (y.contains('='))
        dvalue_valid = y.split_into_directive_and_value(directive, '=', dvalue);
      else
        dvalue_valid = 0;

      if ("unsat" == y)
      {
        set_include_unsaturation_in_atom_type(1);
        if (verbose)
          cerr << "Unsaturation included in atom type\n";
      }
      else if ("unsinca" == y)
      {
        set_unsaturated_includes_aromatic(1);
        if (verbose)
          cerr << "Aromatic atoms can be considered unsaturated\n";
      }
/*    else if ("tt" == y)
      {
        set_topotorsion_atom_types(1);
        if (verbose)
          cerr << "Will use Topological Torsion atom types\n";
      }*/
      else if ("ap" == y)
      {
        set_do_atom_pair_bits(1);
        if (verbose)
          cerr << "Will set bits based on atom pair information\n";
      }
      else if ("och2" == directive)
      {
        set_omit_ch2(dvalue);
        if (verbose)
          cerr << "Omit ch2 set to " << dvalue << endl;
      }
      else if ("iatfc" == y)
      {
        set_formally_charged_atoms_lose_identity(1);
        if (verbose)
          cerr << "Formally charged atoms fingerprinted only as charged form\n";
      }
      else if ("mdc" == y)
      {
        set_maximal_daylight_compatibility(1);
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

        set_include_bits_for_rings(dvalue);
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

        set_bits_for_hydrogen_attachments(dvalue);

        if (verbose)
          cerr << "Attached hydrogen bits set for paths to length " << dvalue << endl;
      }
      else if (y.starts_with("nb=") || y.starts_with("nbits="))
      {
        if (! dvalue_valid || dvalue < 8)
        {
          cerr << "The number of bits in each fingerprint must be a whole +ve number\n";
          return 0;
        }

        set_iwmfingerprint_nbits(dvalue);
        if (verbose)
          cerr << "Fingerprints created with " << dvalue << " bits\n";
      }
      else if ("ip" == y)
      {
        set_isotopic_paths(1);
        if (verbose)
          cerr << "Only isotopically anchored paths will be considered\n";
      }
      else if (y.starts_with("mpl="))
      {
        if (! dvalue_valid || dvalue < 2)
        {
          cerr << "The max path length directive must be followed by a whole non-negative number\n";
          return 0;
        }

        set_max_path_length(dvalue);
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

        set_max_path_length_isotopic_bits(dvalue);
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

        set_terminal_path_bits(dvalue);
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

        set_bits_for_atoms_in_multiple_rings(dvalue);
        if (verbose)
          cerr << "Will set bits for multiple ring membership " << dvalue << " bonds\n";
      }
      else if ("ss" == y)
      {
        only_bits_preserving_substructure_perception = 1;
        if (verbose)
          cerr << "Only bits preserving substructure relationships generated\n";
      }
      else if ("xrcb" == y)
      {
        differentiate_ring_and_chain_bonds = 1;
        if (verbose)
          cerr << "Will differentiate between ring and chain bonds\n";
      }
      else if ("help" == y)
      {
        display_misc_fingerprint_options(cerr, 'Y');
        return 0;
      }
      else 
      {
        cerr << "Unrecognised -" << flag << " qualifier '" << y << "'\n";
        return 0;
      }
    }
  }

  return 1;

}

int
IWMFingerprint::initialise_atom_typing_specification (const const_IWSubstring & p)
{
  return _atom_typing_specification.build(p);
}

#ifdef COUNT_TIMES_ATOM_IN_PATH

class Possible_Fingerprint_Addition
{
  private:
    Set_of_Atoms _atoms;
    resizable_array<int> _bonds;    // numeric bond code

    int _desirability;

  public:
    Possible_Fingerprint_Addition();

    void add_atom(const atom_number_t s) { _atoms.add(s);}
    void add_atom(const atom_number_t a, const int b) { _atoms.add(a);_bonds.add(b);}

    int desirability() const { return _desirability;}

    const Set_of_Atoms & atoms() const { return _atoms;}
    const resizable_array<int> & bonds() const { return _bonds;}

    int compute_desirability(const int * times_in_path);
};

Possible_Fingerprint_Addition::Possible_Fingerprint_Addition()
{
  _desirability = 0;
}

int
Possible_Fingerprint_Addition::compute_desirability(const int * times_in_path)
{
  _desirability = 0;

  const int n = _atoms.number_elements();

  int mintip = std::numeric_limits<int>::max();

  int sum = 0;

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = _atoms[i];

    const int x = times_in_path[j];

    sum += x;

    if (x < mintip)
      mintip = x;
  }

  _desirability = 100000 * mintip + sum;

  return _desirability;
}

class PFA_Comparator
{
  private:
  public:
    int operator() (const Possible_Fingerprint_Addition * p1, const Possible_Fingerprint_Addition * p2) const
                            {
                              if (p1->desirability() < p2->desirability())
                                return -1;
                              if (p1->desirability() > p2->desirability())
                                return 1;
                              return 0;
                            }
};


static void
determine_max_min_times_in_path(const int * times_in_path,
                                const int n,
                                int & mintip,
                                int & maxtip)
{
  mintip = times_in_path[0];
  maxtip = times_in_path[0];

  for (int i = 1; i < n; ++i)
  {
    if (times_in_path[i] < mintip)
      mintip = times_in_path[i];
    else if (times_in_path[i] > maxtip)
      maxtip = times_in_path[i];
  }

  return;
}

int
IWMFingerprint::_extra_bits_to_equalise_times_in_path_3(Molecule & m,
                                        MFingerprint & mfc,
                                        const int * include_these_atoms,
                                        const int niter)
{
  const int matoms = m.natoms();

  const int * times_in_path = mfc.times_in_path();

  resizable_array_p<Possible_Fingerprint_Addition> possible_additions;

  PFA_Comparator pfac;

  for (int i = 0; i < niter; ++i)
  {
    possible_additions.resize_keep_storage(0);

    int maxtip, mintip;
    determine_max_min_times_in_path(times_in_path, matoms, mintip, maxtip);
    int cutoff = (mintip + maxtip) / 2;

    for (int j = 0; j < matoms; ++j)
    {
      if (times_in_path[j] > cutoff)
        continue;

      if (nullptr != include_these_atoms && 0 == include_these_atoms[j])
        continue;
  
      _all_3_paths_with_centre(m, j, include_these_atoms, times_in_path, possible_additions);
    }

    const int n = possible_additions.number_elements();
//  cerr << "Iteration " << i << " found " << n << " possible additions\n";

    if (0 == n)   // seems unlikely
      return 0;

    possible_additions.iwqsort(pfac);

// We need to add new fingerprints carefully. All things that are equivalent must be done together.

    for (int j = 0; j < n; ++j)
    {
      if (j > 0 && pfac(possible_additions[0], possible_additions[j]))    // are different
        break;

      mfc.generate_bits(possible_additions[j]->atoms(), possible_additions[j]->bonds(), _bvector, _auxiliary_bvector);
    }
  }

  return 1;
}

int
IWMFingerprint::_all_3_paths_with_centre(const Molecule & m,
                                         const atom_number_t zatom,
                                         const int * include_these_atoms,
                                         const int * times_in_path,
                                         resizable_array_p<Possible_Fingerprint_Addition> & possible_additions) const
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b1 = a->item(i);

    const atom_number_t ci = b1->other(zatom);

    if (nullptr != include_these_atoms && 0 == include_these_atoms[ci])
      continue;

    for (int j = i + 1; j < acon; ++j)
    {
      const Bond * b2 = a->item(j);

      const atom_number_t cj = b2->other(zatom);

      if (nullptr != include_these_atoms && 0 == include_these_atoms[cj])
        continue;

      Possible_Fingerprint_Addition * p = new Possible_Fingerprint_Addition();

      p->add_atom(ci);
      p->add_atom(zatom, numeric_bond_code(b1));
      p->add_atom(cj, numeric_bond_code(b2));

      p->compute_desirability(times_in_path);

      possible_additions.add(p);
    }
  }

  return 1;
}

int
IWMFingerprint::_all_4_paths_from_bond(const Molecule & m,
                                       const Bond * b,
                                       const int * include_these_atoms,
                                       const int * times_in_path,
                                       resizable_array_p<Possible_Fingerprint_Addition> & possible_additions) const
{
  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  const int a1con = aa1->ncon();
  const int a2con = aa2->ncon();

  for (int i = 0; i < a1con; ++i)
  {
    const Bond * b1 = aa1->item(i);

    const atom_number_t x1 = b1->other(a1);

    if (a2 == x1)
      continue;

    if (nullptr != include_these_atoms && 0 == include_these_atoms[x1])
      continue;

    for (int j = 0; j < a2con; ++j)
    {
      const Bond * b2 = aa2->item(j);

      const atom_number_t x2 = b2->other(a2);

      if (a1 == x2)
        continue;

      if (nullptr != include_these_atoms && 0 == include_these_atoms[x2])
        continue;

      Possible_Fingerprint_Addition * p = new Possible_Fingerprint_Addition();
      p->add_atom(x1);
      p->add_atom(a1, numeric_bond_code(b1));
      p->add_atom(a2, numeric_bond_code(b));
      p->add_atom(x2, numeric_bond_code(b2));

      p->compute_desirability(times_in_path);

      possible_additions.add(p);
    }
  }

  return 1;
}

int
IWMFingerprint::_extra_bits_to_equalise_times_in_path_4(Molecule & m,
                                        MFingerprint & mfc,
                                        const int * include_these_atoms,
                                        const int niter)
{
  const int matoms = m.natoms();

  const int nedges = m.nedges();

  const int * times_in_path = mfc.times_in_path();

  resizable_array_p<Possible_Fingerprint_Addition> possible_additions;

  PFA_Comparator pfac;

  for (int i = 0; i < niter; ++i)
  {
    possible_additions.resize_keep_storage(0);

    int maxtip, mintip;
    determine_max_min_times_in_path(times_in_path, matoms, mintip, maxtip);
    int cutoff = (mintip + maxtip) / 2;

    for (int j = 0; j < nedges; ++j)
    {
      const Bond * b = m.bondi(j);

      const atom_number_t a1 = b->a1();
      const atom_number_t a2 = b->a2();

      if (times_in_path[a1] > cutoff && times_in_path[a2] > cutoff)
        continue;

      if (nullptr != include_these_atoms && (0 == include_these_atoms[a1] || 0 == include_these_atoms[a2]))
        continue;

      _all_4_paths_from_bond(m, b, include_these_atoms, times_in_path, possible_additions);
    }

    const int n = possible_additions.number_elements();
//  cerr << "Iteration " << i << " found " << n << " possible additions\n";

    if (0 == n)   // seems unlikely
      return 0;

    possible_additions.iwqsort(pfac);

// We need to add new fingerprints carefully. All things that are equivalent must be done together.

    for (int j = 0; j < n; ++j)
    {
      if (j > 0 && pfac(possible_additions[0], possible_additions[j]))    // are different
        break;

      mfc.generate_bits(possible_additions[j]->atoms(), possible_additions[j]->bonds(), _bvector, _auxiliary_bvector);
    }
  }

  return 1;
}

template resizable_array_p<Possible_Fingerprint_Addition>::resizable_array_p();
template resizable_array_p<Possible_Fingerprint_Addition>::~resizable_array_p();


#endif
