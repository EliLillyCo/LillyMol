#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"

#include "maccskeys_fn5.h"

using std::cerr;
using std::endl;

void
MACCSKeys::_default_values ()
{
  _aromatic_nitrogens_do_not_have_hydrogens = 0;

  _nbits = 192;

  return;
}

MACCSKeys::MACCSKeys ()
{
  _default_values();

  return;
}

int
MACCSKeys::set_nbits (int s)
{
  if (s <= 0)
  {
    cerr << "MACCSKeys::set_nbits:invalid value " << s << endl;
    return 0;
  }

  if (0 != s % 8)
  {
    cerr << "MACCSKeys::set_nbits:nbits must be a multiple of 8, " << s << " invalid\n";
    return 0;
  }

  if (s > 192)
  {
    cerr << "MACCSKeys::set_nbits:we have a max of 192 keys defined, " << s << " unsupported\n";
    return 0;
  }

  _nbits = s;

  return 1;
}

static int level2_threshold [] = {
  2,     // key 0   2.32637
  1,     // key 1   1.29712
  1,     // key 2   1.282324
  1,     // key 3   1.211719
  2,     // key 4   1.625305
  1,     // key 5   1.133073
  2,     // key 6   1.730573
  1,     // key 7   1.0334
  1,     // key 8   1.05614
  1,     // key 9   1.187549
  1,     // key 10   1.248802
  1,     // key 11   1.034698
  2,     // key 12   1.970081
  1,     // key 13   1.365192
  1,     // key 14   1.197764
  2,     // key 15   1.588901
  1,     // key 16   1.285887
  1,     // key 17   1.067911
  2,     // key 18   1.685646
  1,     // key 19   1.010039
  2,     // key 20   2.274425
  1,     // key 21   1.336127
  1,     // key 22   1.044099
  1,     // key 23   1.041522
  1,     // key 24   1.088208
  1,     // key 25   1.071293
  1,     // key 26   1.140006
  2,     // key 27   1.928823
  1,     // key 28   1.47359
  1,     // key 29   1.150772
  1,     // key 30   1.086388
  2,     // key 31   1.863367
  1,     // key 32   1.181033
  1,     // key 33   1.091568
  1,     // key 34   1.384407
  2,     // key 35   2.331103
  1,     // key 36   1.078796
  1,     // key 37   1.038228
  1,     // key 38   1.294023
  1,     // key 39   1.053052
  1,     // key 40   1.039636
  1,     // key 41   1.34122
  2,     // key 42   2.227006
  1,     // key 43   1.128487
  1,     // key 44   1.453928
  3,     // key 45   2.56371
  1,     // key 46   1.098094
  2,     // key 47   1.771517
  1,     // key 48   1.46953
  1,     // key 49   1.10763
  2,     // key 50   2.312624
  1,     // key 51   1.076154
  1,     // key 52   1.288975
  1,     // key 53   1.24555
  1,     // key 54   1.225561
  1,     // key 55   1.018173
  1,     // key 56   1.165394
  1,     // key 57   1.218056
  2,     // key 58   2.160829
  1,     // key 59   1.138902
  2,     // key 60   2.129522
  1,     // key 61   1.127266
  3,     // key 62   2.628532
  2,     // key 63   2.283072
  1,     // key 64   1.257803
  4,     // key 65   3.642178
  1,     // key 66   1.172121
  1,     // key 67   1.052929
  2,     // key 68   1.821157
  1,     // key 69   1.075863
  1,     // key 70   1.319902
  2,     // key 71   1.689758
  1,     // key 72   1.443498
  2,     // key 73   1.930519
  2,     // key 74   1.719574
  1,     // key 75   1.415125
  5,     // key 76   4.742134
  4,     // key 77   3.714824
  3,     // key 78   2.695983
  2,     // key 79   2.040346
  2,     // key 80   2.316072
  2,     // key 81   1.675143
  1,     // key 82   1.146073
  1,     // key 83   1.290493
  1,     // key 84   1.068638
  1,     // key 85   1.411451
  1,     // key 86   1.454847
  1,     // key 87   1.205534
  1,     // key 88   1.203188
  2,     // key 89   1.685763
  2,     // key 90   1.734786
  2,     // key 91   2.046442
  1,     // key 92   1.335854
  1,     // key 93   1.394207
  2,     // key 94   2.140344
  2,     // key 95   1.718768
  1,     // key 96   1.32123
  2,     // key 97   1.972328
  1,     // key 98   1.371072
  1,     // key 99   1.149938
  3,     // key 100   2.52238
  1,     // key 101   1.011411
  2,     // key 102   2.050475
  1,     // key 103   1.334123
  2,     // key 104   1.558624
  3,     // key 105   2.609813
  1,     // key 106   1.230884
  1,     // key 107   1.380271
  2,     // key 108   2.165668
  1,     // key 109   1.423979
  2,     // key 110   1.607156
  3,     // key 111   2.651621
  1,     // key 112   1.268281
  1,     // key 113   1.490815
  1,     // key 114   1.242277
  2,     // key 115   1.736712
  2,     // key 116   1.967522
  2,     // key 117   1.822076
  3,     // key 118   2.649653
  2,     // key 119   1.819278
  3,     // key 120   2.765679
  2,     // key 121   2.416838
  2,     // key 122   2.052335
  1,     // key 123   1.191346
  2,     // key 124   2.082804
  1,     // key 125   1.198947
  1,     // key 126   1.368367
  1,     // key 127   1.211217
  5,     // key 128   4.529893
  5,     // key 129   4.913997
  1,     // key 130   1.147156
  2,     // key 131   1.572001
  2,     // key 132   1.984907
  1,     // key 133   1.389492
  2,     // key 134   1.92625
  1,     // key 135   1.496362
  2,     // key 136   1.806896
  1,     // key 137   1
  3,     // key 138   2.714807
  1,     // key 139   1.206046
  3,     // key 140   2.6453
  2,     // key 141   2.104606
  3,     // key 142   2.964203
  2,     // key 143   1.576038
  2,     // key 144   1.78995
  2,     // key 145   2.392797
  1,     // key 146   1.084459
  1,     // key 147   1.113913
  2,     // key 148   1.589959
  1,     // key 149   1.197222
  2,     // key 150   1.840145
  1,     // key 151   1.439922
  2,     // key 152   1.688804
  3,     // key 153   3.228122
  1,     // key 154   1.471567
  2,     // key 155   2.197356
  4,     // key 156   3.936198
  3,     // key 157   2.81755
  3,     // key 158   3.427067
  1,     // key 159   1.102107
  2,     // key 160   1.676797
  1,     // key 161   1.110776
  13,     // key 162   13.48983
  1,     // key 163   1.044401
  1,     // key 164   1.365189
  3,     // key 165   3.193642
  1,     // key 166   1.405154
  14,     // key 167   14.44107
  1,     // key 168   1.308397
  2,     // key 169   2.26848
  1,     // key 170   1.076979
  1,     // key 171   1.059862
  1,     // key 172   1.244386
  1,     // key 173   1.0756
  2,     // key 174   2.304352
  2,     // key 175   1.769305
  1,     // key 176   1.152277
  2,     // key 177   1.698764
  1,     // key 178   1.380591
  1,     // key 179   1.04041
  1,     // key 180   1.401356
  1,     // key 181   1.147481
  1,     // key 182   1.014824
  1,     // key 183   1.040559
  1,     // key 184   1.222028
  1,     // key 185   1.0666
  4,     // key 186   3.757186
  1,     // key 187   1.068092
  1,     // key 188   1.455155
  1,     // key 189   1.010624
  1,     // key 190   1.020083
  2     // key 191   2.095454
};

class MK_Molecular_Properties
{
  private:
    int _natoms;
    int _nr;

    atomic_number_t * _z;
    int * _ncon;
    const Atom * * _atom;
    int * _ring_membership;
    int * _aromatic;
    int * _in_same_ring;
    int * _implicit_hydrogens;
    int * _attached_heteroatom_count;
    int * _spinach;
    int * _ring_bond_count;
    int * _ring_already_done;

    atom_number_t _first_carbon;
    atom_number_t _last_carbon;

    atom_number_t _first_heteroatom;
    atom_number_t _last_heteroatom;

    atom_number_t _first_ring_atom;
    atom_number_t _last_ring_atom;

    atom_number_t _first_non_ring_atom;
    atom_number_t _last_non_ring_atom;

    atom_number_t _first_oxygen;
    atom_number_t _last_oxygen;

    atom_number_t _first_nitrogen;
    atom_number_t _last_nitrogen;

//  private functions

    int _initialise_ring_arrays (Molecule & m);
    int _initialise_implicit_hydrogens (Molecule & m);

  public:
    MK_Molecular_Properties(Molecule &);
    ~MK_Molecular_Properties();

    inline int natoms () const { return _natoms;}
    inline int nrings () const { return _nr;}

    inline atom_number_t first_carbon() const { return _first_carbon;}
    inline atom_number_t last_carbon() const { return _last_carbon;}

    inline atom_number_t first_oxygen() const { return _first_oxygen;}
    inline atom_number_t last_oxygen() const { return _last_oxygen;}

    inline atom_number_t first_nitrogen() const { return _first_nitrogen;}
    inline atom_number_t last_nitrogen() const { return _last_nitrogen;}

    inline atom_number_t first_heteroatom() const { return _first_heteroatom;}
    inline atom_number_t last_heteroatom() const { return _last_heteroatom;}

    inline atom_number_t first_ring_atom() const { return _first_ring_atom;}
    inline atom_number_t last_ring_atom() const { return _last_ring_atom;}

    inline atom_number_t first_non_ring_atom() const { return _first_non_ring_atom;}
    inline atom_number_t last_non_ring_atom() const { return _last_non_ring_atom;}

    inline const atomic_number_t * z() const { return _z;}
    inline const atomic_number_t * ncon() const { return _ncon;}
    inline const Atom * const * atoms() const { return _atom;}
    inline const int * aromatic() const { return _aromatic;}
    inline const int * implicit_hydrogens() const { return _implicit_hydrogens;}
    inline const int * ring_membership() const { return _ring_membership;}
    inline const int * in_same_ring() const { return _in_same_ring;}
    inline const int * attached_heteroatom_count () const { return _attached_heteroatom_count;}
    inline const int * spinach () const { return _spinach;}
    inline const int * ring_bond_count () const { return _ring_bond_count;}
    inline int * ring_already_done () const { return _ring_already_done;}
};

static int
compute_ring_bond_count (const Atom & a, const int acon)
{
  int rc = 0;

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a[i];
    if (b->nrings() > 0)
      rc++;
  }

  return rc;
}

MK_Molecular_Properties::MK_Molecular_Properties(Molecule & m)
{
  _natoms = m.natoms();
  _nr = m.nrings();

  _z = new atomic_number_t[_natoms];
  m.atomic_numbers(_z);

  _ncon = new int[_natoms];
  m.ncon(_ncon);

  _atom = new const Atom *[_natoms];
  m.atoms(_atom);

  _implicit_hydrogens = new int[_natoms];

  _ring_membership = new int[_natoms];
  m.ring_membership(_ring_membership);

  _aromatic = new int[_natoms];

  if (_nr)
  {
    m.aromaticity(_aromatic);
    _ring_already_done = new int[_nr];
  }
  else
  {
    set_vector(_aromatic, _natoms, 0);
    _ring_already_done = nullptr;
  }

  if (_nr > 0)
    _initialise_ring_arrays(m);
  else
    _in_same_ring = nullptr;

  _attached_heteroatom_count = new int[_natoms];

  _spinach = new int[_natoms];
  m.identify_spinach(_spinach);
  m.compute_aromaticity_if_needed();

  _ring_bond_count = new int[_natoms];
  if (_nr)
  {
    for (int i = 0; i < _natoms; ++i)
    {
      _ring_bond_count[i] = compute_ring_bond_count(*_atom[i], _ncon[i]);
    }
  }
  else
    set_vector(_ring_bond_count, _natoms, 0);

  _initialise_implicit_hydrogens(m);

  _first_carbon = INVALID_ATOM_NUMBER;
  _first_heteroatom = INVALID_ATOM_NUMBER;
  _first_non_ring_atom = INVALID_ATOM_NUMBER;
  _first_ring_atom = INVALID_ATOM_NUMBER;
  _first_oxygen = INVALID_ATOM_NUMBER;
  _first_nitrogen = INVALID_ATOM_NUMBER;

  _last_carbon = INVALID_ATOM_NUMBER;
  _last_heteroatom = INVALID_ATOM_NUMBER;
  _last_ring_atom = INVALID_ATOM_NUMBER;
  _last_non_ring_atom = INVALID_ATOM_NUMBER;
  _last_oxygen = INVALID_ATOM_NUMBER;
  _last_nitrogen = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _natoms; i++)
  {
    _attached_heteroatom_count[i] = m.attached_heteroatom_count(i);

    if (6 == _z[i])
    {
      if (INVALID_ATOM_NUMBER == _first_carbon)
        _first_carbon = i;
      _last_carbon = i + 1;
    }
    else
    {
      if (INVALID_ATOM_NUMBER == _first_heteroatom)
        _first_heteroatom = i;

      _last_heteroatom = i + 1;

      if (7 == _z[i])
      {
        if (INVALID_ATOM_NUMBER == _first_nitrogen)
          _first_nitrogen = i;
        _last_nitrogen = i + 1;
      }
      else if (8 == _z[i])
      {
        if (INVALID_ATOM_NUMBER == _first_oxygen)
          _first_oxygen = i;
        _last_oxygen = i + 1;
      }
    }

    if (0 == _ring_membership[i])
    {
      if (INVALID_ATOM_NUMBER == _first_non_ring_atom)
        _first_non_ring_atom = i;
      _last_non_ring_atom = i + 1;
    }
    else
    {
      if (INVALID_ATOM_NUMBER == _first_ring_atom)
        _first_ring_atom = i;
      _last_ring_atom = i + 1;
    }
  }

  if (INVALID_ATOM_NUMBER == _first_ring_atom)
  {
    _first_ring_atom = 0;
    _last_ring_atom = 0;
  }

  if (INVALID_ATOM_NUMBER == _first_carbon)
  {
    _first_carbon = 0;
    _last_carbon = 0;
  }

  if (INVALID_ATOM_NUMBER == _first_oxygen)
  {
    _first_oxygen = 0;
    _last_oxygen = 0;
  }

  if (INVALID_ATOM_NUMBER == _first_nitrogen)
  {
    _first_nitrogen = 0;
    _last_nitrogen = 0;
  }

  if (INVALID_ATOM_NUMBER == _first_non_ring_atom)
  {
    _first_non_ring_atom = 0;
    _last_non_ring_atom = 0;
  }


  return;
}

MK_Molecular_Properties::~MK_Molecular_Properties()
{
  delete [] _z;
  delete [] _ncon;
  delete [] _atom;
  delete [] _ring_membership;
  delete [] _aromatic;
  if (nullptr != _in_same_ring)
    delete [] _in_same_ring;
  if (nullptr != _ring_already_done)
    delete [] _ring_already_done;
  delete [] _implicit_hydrogens;
  delete [] _attached_heteroatom_count;
  delete [] _spinach;
  delete [] _ring_bond_count;

  return;
}

int
MK_Molecular_Properties::_initialise_ring_arrays (Molecule & m)
{

  _in_same_ring = new_int(_natoms * _natoms);
  assert (nullptr != _in_same_ring);

  for (int i = 0; i < _nr; i++)
  {
    const Ring * r = m.ringi(i);
    int rs = r->number_elements();
    for (int j = 0; j < rs; j++)
    {
      atom_number_t aj = r->item(j);
      for (int k = j + 1; k < rs; k++)
      {
        atom_number_t ak = r->item(k);
        _in_same_ring[aj * _natoms + ak] = 1;
        _in_same_ring[ak * _natoms + aj] = 1;
      }
    }
  }

  return 1;
}

int
MK_Molecular_Properties::_initialise_implicit_hydrogens(Molecule & m)
{
  for (int i = 0; i < _natoms; i++)
  {
    _implicit_hydrogens[i] = const_cast<Atom *>(_atom[i])->implicit_hydrogens();
  }

  return 1;
}

/*
  '[CH3]A[CH3]'
*/

int
MACCSKeys::_key0 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  const auto istop = mpr.last_carbon();
  for (auto i = mpr.first_carbon(); i < istop; ++i)
  {
    if (6 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const auto b = atoms[i]->item(0);

    if (! b->is_single_bond())
      continue;

    const auto x = b->other(i);

    if (4 != ncon[x])
      continue;

//  cerr << "from " << i << " Found middle atom " << x << " ncon " << ncon[x] << endl;

    const auto ax = atoms[x];

    for (auto j = 0; j < ncon[x]; ++j)
    {
      const auto b = ax->item(j);

      if (! b->is_single_bond())
        continue;

      const auto k = b->other(x);

      if (k <= i)
        continue;

      if (6 != z[k])
        continue;

      if (1 != ncon[k])
        continue;

      rc++;
    }
  }

  return rc;
}

int
MACCSKeys::_key1 (Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  int rc = 0;

  const int * ring_membership = mpr.ring_membership();
  const auto matoms = mpr.natoms();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  for (int i = 0; i < matoms; i++)
  {
    if (1 != ncon[i])
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(i);

    if (ring_membership[j])
      rc++;
  }

  return rc;
}

/*
  amide
*/

int
MACCSKeys::_key2 (const Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  const auto matoms = mpr.natoms();

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != atoms[i]->ncon())
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(i);

    if (6 != z[j])
      continue;

    if (3 != ncon[j])
      continue;

    const Atom * c = atoms[j];

    int found_nitrogen = 0;
    int found_another_carbon = 0;

    for (int k = 0; k < 3; k++)
    {
      atom_number_t l = c->other(j, k);

      if (l == i)
        continue;

      if (6 == z[l])
        found_another_carbon++;
      else if (7 == z[l])
        found_nitrogen++;
      else
        break;
    }

    if (1 == found_another_carbon && 1 == found_nitrogen)
      rc++;
  }

  return rc;
}

/*
  Singly connected aromatic ring. Just one exocyclic bond
*/

int
MACCSKeys::_key3 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const int * ncon = mpr.ncon();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    int ring_size = ri->number_elements();

    int connected_outside_ring = 0;

    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      if (2 != ncon[k])
        connected_outside_ring++;
    }

    if (1 == connected_outside_ring)
      rc++;
  }

  return rc;
}

/*
  '[#7,#8]@[D2]@[D2]@[#7,#8]'
*/

int
MACCSKeys::_key4 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  const auto istop = mpr.last_heteroatom();

  for (int i = mpr.first_heteroatom(); i < istop; ++i)
  {
    if (6 == z[i])
      continue;

    if (0 == ring_membership[i])
      continue;

    const auto ai = atoms[i];

    for (auto i1 = 0; i1 < ncon[i]; ++i1)
    {
      const auto b1 = ai->item(i1);

      if (0 == b1->nrings())
        continue;

      const auto a2 = b1->other(i);

      if (2 != ncon[a2])
        continue;

      const auto aa2 = atoms[a2];

      for (auto i2 = 0; i2 < ncon[a2]; ++i2)
      {
        const Bond * b2 = aa2->item(i2);

        if (0 == b2->nrings())
          continue;

        const auto a3 = b2->other(a2);

        if (a3 == i)
          continue;

        if (2 != ncon[a3])
          continue;

        const auto aa3 = atoms[a3];

        for (auto i3 = 0; i3 < ncon[a3]; ++i3)
        {
          const Bond * b = aa3->item(i3);

          if (0 == b->nrings())
            continue;

          const auto a4 = b->other(a3);

          if (a4 == a2)
            continue;

          if (6 == z[a4])
            continue;

          if (a4 > i)
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  Para substituted aromatic rings
*/

static int
bonded_to_non_ring_atom (const Molecule & m,
                         const int * ring_membership,
                         atom_number_t zatom,
                         const Atom * a)
{
  int acon = a->ncon();   // should be 3

  for (int i = 0; i < acon; i++)
  {
    if (0 == ring_membership[a->other(zatom, i)])
      return 1;
  }

  return 0;
}

int
MACCSKeys::_key5 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const int * ring_membership = mpr.ring_membership();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (6 != ri->number_elements())
      continue;

    int first_3_substitued = -1;
    int second_3_substitued = -1;
    int third_three_substitued = -1;
    for (int j = 0; j < 6; j++)
    {
      atom_number_t k = ri->item(j);

      if (2 == ncon[k])
        continue;

      if (! bonded_to_non_ring_atom(m, ring_membership, k, atoms[k]))
        continue;

      if (-1 == first_3_substitued)
        first_3_substitued = j;
      else if (-1 == second_3_substitued)
        second_3_substitued = j;
      else
        third_three_substitued++;
    }

    if (-1 != third_three_substitued)
      continue;

    if (first_3_substitued + 3 == second_3_substitued)
      rc++;
  }

  return rc;
}

/*
  Key 6 is lanthanide
  Oct 2007, change to a-!@[!#6G0D>1T0R0]
*/

int
MACCSKeys::_key6 (const Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();
  const int * ncon = mpr.ncon();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * ai = atoms[i];
    for (int j = 0; j < 3; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atomic_number_t k = b->other(i);

      if (ring_membership[k])
        continue;

      if (1 == ncon[k])
        continue;

      const Atom * ak = atoms[k];

      if (ncon[k] < ak->nbonds())
        continue;

      if (attached_heteroatom_count[k] > 0)
        continue;

      rc++;
    }
  }

  return rc;
}

/*
  Key 7 is Group VB, VIB, VIIB   V, etc
  As I don't know what this is, just detect V
  Oct 2007. Change to [CD2H2R0]-[CD2H2R0]
  Change to count number of rings with multiple OH attachments
  Problems with SSSR stuff, just count number of OH's bonded to ring atoms
*/

#ifdef OLD_MK_NOT_USED_ANYMORE
int
MACCSKeys::_key7 (Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto * ri = m.ringi(i);

    auto jstop = ri->number_elements();

    if (jstop > 7)     // arbitrary decision
      continue;

    int oxygens = 0;
    for (auto j = 0; j < jstop; ++j)
    {
      atom_number_t k = ri->item(j);

      if (2 == ncon[k])
        continue;

      const Atom * ak = atoms[k];

      for (auto l = 0; l < ncon[k]; ++l)
      {
        const auto o = ak->other(k, l);

        if (8 != z[o])
          continue;

        if (1 == ncon[o])
          oxygens++;
      }
    }

    if (oxygens > 1)
      rc++;
  }

  return rc;
}
#endif

int
MACCSKeys::_key7 (Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_bond_count = mpr.ring_bond_count();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();
  const Atom * const * atoms = mpr.atoms();

  int istop = mpr.last_oxygen();

  int rc = 0;

  for (int i = mpr.first_oxygen(); i < istop; ++i)
  {
    if (8 != z[i])
      continue;
    if (1 != ncon[i])
      continue;
    if (1 != implicit_hydrogens[i])
      continue;

    const auto r = atoms[i]->other(i, 0);

    if (ring_bond_count[r])
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 8 is a four membered heterocycle

  too rare, discard
*/

int
MACCSKeys::_key8 (Molecule & m,
      const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);
    if (4 != ri->number_elements())
      continue;

    int heteroatoms_found = 0;
    for (int j = 0; j < 4; j++)
    {
      atom_number_t k = ri->item(j);
      if (6 != z[k])
        heteroatoms_found++;
    }

    if (1 == heteroatoms_found)
      rc++;
  }

  return rc;
}
#endif

/*
  Number of rings smaller than 5
*/

int
MACCSKeys::_key8 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (ri->number_elements() > 4)
      break;

    rc++;
  }

  return rc;
}

/*
  [!#6]([aR2])[aR2]
*/

int
MACCSKeys::_key9 (Molecule & m,
                  const MK_Molecular_Properties & mpr) const
{
  if (m.nrings() < 3)
    return 0;

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_bond_count = mpr.ring_bond_count();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  const auto istop = mpr.last_heteroatom();

  for (auto i = mpr.first_heteroatom(); i < istop; ++i)
  {
    if (6 == z[i])
      continue;

    if (ncon[i] < 2)
      continue;

    int aromatics = 0;

    const auto ai = atoms[i];

    for (auto j = 0; j < ncon[i]; ++j)
    {
      atom_number_t k = ai->other(i, j);

      if (! aromatic[k])
        continue;

      if (ring_bond_count[k] > 2)
        aromatics++;
    }

    if (aromatics > 1)
      rc++;
  }

  return rc;
}

/*
  Aromatic ring with > 1 heteroatom
*/

int
MACCSKeys::_key10 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const atomic_number_t * z = mpr.z();

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    int heteroatoms = 0;

    const auto jstop = ri->number_elements();

    for (auto j = 0; j < jstop; ++j)
    {
      auto k = ri->item(j);

      if (6 != z[k])
        heteroatoms++;
    }

    if (heteroatoms > 1)
      rc++;
  }

  return rc;
}

/*
  Key 11 is four membered ring
*/

int
MACCSKeys::_key11 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (4 == r->number_elements())
    {
      rc++;
    }
  }

  return rc;
}

/*
  CH3AACH2A
*/

int
MACCSKeys::_key12 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (1 != ncon[i])
      continue;

    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (a2 == i)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t ch2 = atom_a2->other(a2, l);

          if (ch2 == a1)
            continue;

          if (6 != z[ch2])
            continue;
  
          if (implicit_hydrogens[ch2] < 2)
            continue;
  
          if (ncon[ch2] > 1)
            rc++;
        }
      }
    }
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 13 is ON(C)C 
*/

int
MACCSKeys::_key13 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istop = mpr.last_nitrogen();
  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    int acon = a->ncon();

    if (acon < 3)
      continue;

    int carbon_attachments = 0;
    int oxygen_attachments = 0;
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (6 == z[k])
        carbon_attachments++;
      else if (8 == z[k])
        oxygen_attachments++;
    }

    if (1 == oxygen_attachments && 2 == carbon_attachments)
      rc++;
  }

  return rc;
}
#endif

/*
  key 13 is a-!@[G1R0]'
*/

int
MACCSKeys::_key13 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_bond_count = mpr.ring_bond_count();
  const Atom * const * atoms = mpr.atoms();

  const int matoms = mpr.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == aromatic[i])    // start with aromatic atoms
      continue;

    if (3 != ncon[i])
      continue;

    const auto ai = atoms[i];

    if (2 != ring_bond_count[i])
      continue;

    for (int j = 0; j < 3; ++j)
    {
      const Bond * b = ai->item(j);
      if (b->nrings())
        continue;

      const auto k = b->other(i);

      if (ring_bond_count[k])
        continue;

      if (ncon[k] < m.nbonds(k))
        rc++;
    }
  }

  return rc;
}

/*
  Key 14 is S-S, Disulphide

  There are too few of these so change
  Oct 2007 aromatic nitrogen
  Sep 2013 [CT3]
*/

int
MACCSKeys::_key14 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();

  int rc = 0;

  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (3 == attached_heteroatom_count[i])
      rc++;
  }

  return rc;
}

/*
  Key 15 is OC(O)O C attached to 3 O 

  Note an ambiguity, should we count C with 4 oxygens as a hit.  For
  now, we do not

  Too few molecules hit this. 
  Oct 2007 change to [aR2]:[aD2]:[aD3]-!@A
*/

int
MACCSKeys::_key15 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_bond_count = mpr.ring_bond_count();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (3 != ring_bond_count[i])
      continue;

    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    List_of_Ring_Sizes lors;
    m.ring_sizes_for_atom(i, lors);

    if (lors.last_item() > 7)
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * bj = ai->item(j);

      if (! bj->is_aromatic())
        continue;

      atom_number_t c1 = bj->other(i);

      if (2 != ncon[c1])
        continue;

      if (! aromatic[c1])
        continue;

      const Atom * ac = atoms[c1];
      for (int k = 0; k < 2; k++)
      {
        atom_number_t c2 = ac->other(c1, k);

        if (c2 == i)
          continue;

        if (! aromatic[c2])   // not sure if that could happen
          continue;

        if (2 != ring_bond_count[c2])
          continue;

        if (3 == ncon[c2])
          rc++;
      }
    }
  }

  return rc;
}

/*
  '[D1]-[D2]-a'
*/

int
MACCSKeys::_key16 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  if (0 == mpr.nrings())
    return 0;

  const auto matoms = mpr.natoms();
  const int * aromatic = mpr.aromatic();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (1 != ncon[i])
      continue;

    const auto a2 = atoms[i]->other(i, 0);

    if (2 != ncon[a2])
      continue;

    const auto aa1 = atoms[a2];

    for (auto j = 0; j < 2; ++j)
    {
      const auto a3 = aa1->other(a2, j);

      if (a3 == i)
        continue;

      if (aromatic[a3])
        rc++;
    }
  }

  return rc;
}

/*
  Key 17 is C#C, alkyne

  Does not hit enough molecules, change to any triple bond
*/

#ifdef OLD_MK_NOT_USED_ANYMORE

int
MACCSKeys::_key17 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] > 2 || 0 == ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (ai->nbonds() <= 2)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item(j);

      if (! b->is_triple_bond())
        continue;

      atom_number_t k = b->other(i);
      if (k > i && 6 == z[k])
        rc++;
    }
  }

  return rc;
}
#endif

int
MACCSKeys::_key17 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const int nb = m.nedges();

  int rc = 0;

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = m.bondi(i);

    if (b->is_triple_bond())
      rc++;
  }

  return rc;
}

/*
  Key 18 is Group IIIa, B, etc.
  Oct 2007, change to a-[D1]
*/

int
MACCSKeys::_key18 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      const Bond * bj = ai->item(j);
      if (bj->nrings())
        continue;

      atom_number_t k = bj->other(i);

      if (1 == ncon[k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 19 is seven membered ring
*/

int
MACCSKeys::_key19 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (7 == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  Key 20 is Si
  Oct 2007, change to [aD3R2](:a)(:a):a
*/

int
MACCSKeys::_key20 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (1 == ring_membership[i])
      continue;

    int aromatic_connections = 0;

    const Atom *ai = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      const Bond * bj = ai->item(j);

      if (bj->is_aromatic())
        aromatic_connections++;
    }

    if (3 == aromatic_connections)
      rc++;
  }

  return rc;
}

/*
  Key 21 is C=C(Q)Q, 1,1-dihetero vinyl, purine

  Dec 2004. Programme was sensitive to Kekule form. If I changed
  things to just C=C (aliphatic) then hardly anything hits. So,
  I allow double or aromatic bond.
  
  Special problems with:
  C12=C(NC(=O)NC3=CC=CC=C3)N=CC=C1C=CC=C2
*/

int
MACCSKeys::_key21 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;
    if (3 != ncon[i])
      continue;
     
    const Atom * a = atoms[i];

    if (4 != a->nbonds())
      continue;

    int found_vinyl = 0;
    int found_heteroatoms = 0;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = a->item(j);

      atom_number_t k = b->other(i);

      if (6 != z[k])
          found_heteroatoms++;
      else if (b->is_double_bond() || b->is_aromatic())
        found_vinyl++;
    }
    
    if (1 == found_vinyl && 2 == found_heteroatoms)
      rc++;
  }

  return rc;
}

/*
  Key 22 is a three membered ring
*/

int
MACCSKeys::_key22 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (3 == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  Key 23 is NC(O)O, C attached to 2 O, 1 N, carbamyl, BOC
*/

int
MACCSKeys::_key23 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])    // start search with N
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t c = ai->other(i, j);
      if (6 != z[c])
        continue;

      const Atom * ac = atoms[c];

      int oxygen_count = 0;
      for (int l = 0; l < ncon[c]; l++)
      {
        atom_number_t o = ac->other(c, l);
        if (8 == z[o])
          oxygen_count++;
      }

      if (2 == oxygen_count)
        rc++;
    }
  }

  return rc;
}

/*
  Key 24 is N-O, Hydroxylamine, oxime, isoxazole
*/

int
MACCSKeys::_key24 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;

  const int istop = mpr.last_nitrogen();
  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (8 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 25 is  NC(N)N
*/

int
MACCSKeys::_key25 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    int nitrogen_count = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (7 == z[k])
        nitrogen_count++;
    }

    if (3 == nitrogen_count)
    {
      rc++;
    }
  }

  return rc;
}

/*
  Key 26 is C$=C($A)($A)
  Dollar means that bond in a ring

  Note that I am somewhat nervous about this with aromatic systems.
*/

int
MACCSKeys::_key26 (Molecule & m,
        const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * in_same_ring = mpr.in_same_ring();

  if (nr < 2)        // system must have at least two rings
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  if (mpr.first_carbon() > istart)
    istart = mpr.first_carbon();

  int istop = mpr.last_ring_atom();
  if (mpr.last_carbon() < istop)
    istop = mpr.last_carbon();

  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (ring_membership[i] < 2)
      continue;

    const Atom * a = atoms[i];

    if (4 != a->nbonds())
      continue;

    int found_vinyl_carbon = 0;
    int found_other_ring_atoms = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      atom_number_t k = b->other(i);

      if (0 == ring_membership[k])
        continue;

      if (0 == in_same_ring[i * matoms + k])
        continue;

      if (6 == z[k] && (b->is_double_bond() || b->is_aromatic()))
        found_vinyl_carbon++;
      else
        found_other_ring_atoms++;
    }

    if (1 == found_vinyl_carbon && found_other_ring_atoms >= 2)
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 27 is Iodine
  Too rare
*/

int
MACCSKeys::_key27 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (53 == z[i])
    {
      rc++;
    }
  }

  return rc;
}
#endif

/*
  Key 27 is nitrogen atoms in a 2 aromatic ring fused system
*/

int
MACCSKeys::_key27 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * aromatic = mpr.aromatic();

  int rc = 0;

  const auto nr = m.nrings();

  const auto matoms = mpr.natoms();

  if (nr < 2)
    return 0;

  int * ring_already_done = mpr.ring_already_done();
  set_vector(ring_already_done, nr, 0);

  for (int i = 0; i < nr; ++i)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (1 != ri->fused_ring_neighbours())
      continue;

    for (int j = i + 1; j < nr; ++j)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);

      if (! rj->is_aromatic())
        continue;

      if (rj->fused_system_identifier() != ri->fused_system_identifier())
        continue;

      if (! ri->is_fused_to(rj))
        continue;

      if (1 != rj->fused_ring_neighbours())
        continue;

      ring_already_done[j] = 1;

      for (int k = 0; k < matoms; ++k)
      {
        if (7 != z[k])
          continue;

        if (! aromatic[k])
          continue;

        if (ri->contains(k) || rj->contains(k))
          rc++;
      }
    }
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 28 is QCH2Q
*/

int
MACCSKeys::_key28 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;
    if (2 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (implicit_hydrogens[i] < 2)
      continue;

    int heteroatom_count = 0;
    for (int j = 0; j < 2; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        heteroatom_count++;
    }

    if (2 == heteroatom_count)
      rc++;
  }

  return rc;
}
#endif

/*
  '[OD1,NH>0]~*~*~*~*~[OD1,NH>0]'

  with some ring bonds in there - that part not implemented
*/

//#define DEBUG_KEY28

static int
rotatable_bonds_between (Molecule & m,
                         const atom_number_t astart,
                         const atom_number_t zend,
                         const int d,
                         const MK_Molecular_Properties & mpr)
{
  const auto a = mpr.atoms()[astart];

  const int acon = mpr.ncon()[astart];

  int maxrotb = 0;

#ifdef DEBUG_KEY28
  cerr << " from atom " << astart << " distance " << d << endl;
#endif

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const auto j = b->other(astart);

#ifdef DEBUG_KEY28
    cerr << "   btw " << j << " and " << zend << " dist " << m.bonds_between(j, zend) << ", cmp " << d << endl;
#endif

    if (m.bonds_between(j, zend) != d)
      continue;

    int rotatable = b->is_single_bond() && 0 == b->nrings();

    if (j == zend)
      return rotatable;

    int tmp = rotatable + rotatable_bonds_between(m, j, zend, d-1, mpr);

#ifdef DEBUG_KEY28
    cerr << "  from " << astart << " to " << zend << " rbc " << tmp << endl;
#endif

    if (tmp > maxrotb)
      maxrotb = tmp;
  }

  return maxrotb;
}

int
MACCSKeys::_key28 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  const int matoms = mpr.natoms();

#ifdef DEBUG_KEY28
  if (1)
  {
    Molecule mcopy(m);
    for (int i = 0; i < matoms; ++i)
    {
      mcopy.set_isotope(i, i);
    }
    cerr << mcopy.smiles() << ' ' << m.name() << '\n';
  }
#endif
  
  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (7 == z[i] && implicit_hydrogens[i] > 0)
      ;
    else if (8 == z[i] && 1 == ncon[i])
      ;
    else
      continue;

    for (int j = i + 1; j < matoms; ++j)
    {
      if (6 == z[i])
        continue;

      if (15 != z[i] + z[j])   // nitrogen and oxygen only
        continue;

      if (7 == z[j] && implicit_hydrogens[j] > 0)
        ;
      else if (8 == z[j] && 1 == ncon[j])
        ;
      else
        continue;

      if (6 != m.bonds_between(i, j))
        continue;
      
      const auto rtb = rotatable_bonds_between(m, i, j, 5, mpr);
#ifdef DEBUG_KEY28
      cerr << " btw " << i << " and " << j << " rtb " << rtb << endl;
#endif
      if (rtb <= 3)
        rc++;
    }
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 29 is Phosphorus
*/

int
MACCSKeys::_key29 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (15 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

#endif

/*
  [CD4T0]
*/

int
MACCSKeys::_key29 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();
  const int * ncon = mpr.ncon();

  const int matoms = mpr.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i)
  {
    if (6 != z[i])
      continue;

    if (0 != attached_heteroatom_count[i])
      continue;

    if (4 != ncon[i])
      continue;

    rc++;
  }

  return rc;
}

/*
  Key 30 is CQ(C)(C)A
  Too few molecules hit this
  Oct 2007 [nH]
*/

int
MACCSKeys::_key30 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();
  const int * aromatic = mpr.aromatic();

  int rc = 0;

  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (! aromatic[i])
      continue;

   if (implicit_hydrogens[i])
     rc++;
  }

  return rc;
}

/*
  Key 31 is QX - halogen bonded to a heteroatom
  Never hits anything, change:
    Oct 2007 change to [CR0]=[CR0]
*/

int
MACCSKeys::_key31 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ring_membership[i])
      continue;

    int icon = ncon[i];
    if (1 == icon)
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < icon; j++)
    {
      const Bond * bj = ai->item(j);

      if (! bj->is_double_bond())
        continue;

      atom_number_t k = bj->other(i);

      if (6 != z[k])
        continue;

      if (ring_membership[k])
        continue;

      rc++;
    }
  }

  return rc;
}

/*
  '[D1]-[D3]-[D1]'
*/

int
MACCSKeys::_key32 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();

  const auto matoms = mpr.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (3 != ncon[i])
      continue;

    if (ring_membership[i])
      continue;

    const Atom * a = atoms[i];

    int d1 = 0;
    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = a->other(i, j);
      if (1 == ncon[k])
        d1++;
    }

    if (2 == d1)
      rc++;
  }

  return rc;
}

/*
  Key 33 is NS
*/

int
MACCSKeys::_key33 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (7 == z[k])
        rc++;
    }
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 34 is CH2=A
*/

int
MACCSKeys::_key34 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (implicit_hydrogens[i] < 2)
      continue;

    if (2 != a->nbonds())
      continue;

    rc++;
  }

  return rc;
}
#endif

/*
  key 34 is '[CD1,F]-[a]:[aD3R1]'
*/

int
MACCSKeys::_key34 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  const int matoms = mpr.natoms();
  if (matoms == 1) {
    return 0;
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (6 == z[i] && 1 == ncon[i])
      ;
    else if (9 == z[i])
      ;
    else
      continue;

    const auto j = atoms[i]->other(i, 0);

    if (! aromatic[j])
      continue;

    const auto aj = atoms[j];

    if (3 != ncon[j])    // very hard to imagine
      continue;

    int found_d3 = 0;
    for (int k = 0; k < 3; ++k)
    {
      const Bond * b = aj->item(k);

      if (! b->is_aromatic())
        continue;

      const auto l = b->other(j);

      if (aromatic[l] && 3 == ncon[l])
      {
        found_d3 = 1;
        break;
      }
    }

    if (found_d3)
      rc++;
  }

  return rc;
}

/*
  Key 35 is Group 1A (alkali metal)
  Oct 2007. change to [A]@-[aD3R2](:a):a
*/

int
MACCSKeys::_key35 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (1 == ring_membership[i])
      continue;

    const Atom * ai = atoms[i];

    int aromatic_attachments = 0;
    int aliphatic_attachments = 0;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = ai->item(j);

      if (0 == b->nrings())   // how could that happen?
        continue;

      if (b->is_aromatic())
        aromatic_attachments++;
      else
        aliphatic_attachments++;
    }

    if (2 == aromatic_attachments && 1 == aliphatic_attachments)
      rc++;
  }

  return rc;
}

/*
  Key 36 is S heterocycle
*/

int
MACCSKeys::_key36 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  if (mpr.first_heteroatom() > istart)
    istart = mpr.first_heteroatom();

  int istop = mpr.last_ring_atom();
  if (mpr.last_heteroatom() < istop)
    istop = mpr.last_heteroatom();

  for (int i = istart; i < istop; i++)
  {
    if (16 == z[i] && ring_membership[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  Key 37 is NC(O)N
*/

int
MACCSKeys::_key37 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * a = atoms[i];

    int nitrogen = 0;
    int oxygen = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (7 == z[k])
        nitrogen++;
      else if (8 == z[k])
        oxygen++;
    }

    if (nitrogen >= 2 && oxygen)
      rc++;
  }

  return rc;
}

/*
  Key 38 is NC(C)N
*/

int
MACCSKeys::_key38 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    int carbon = 0;
    int nitrogen = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        carbon++;
      else if (7 == z[k])
        nitrogen++;
    }

    if (nitrogen >= 2 && carbon)
      rc++;
  }


  return rc;
}

/*
  Key 39 is OS(O)O
  Too rare, get rid of the 3rd oxygen requirement
*/

int
MACCSKeys::_key39 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    int oxygen = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (8 == z[k])
        oxygen++;
    }

    if (oxygen >= 2)    // change to require two
      rc++;
  }

  return rc;
}

/*
  '[CD4]([D>1])([D>1])[D>1]'
*/

int
MACCSKeys::_key40 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (4 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    int dg1 = 0;

    for (int j = 0; j < 4; j++)
    {
      const Bond * b = a->item(j);

      if (! b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (ncon[k] > 1)
        dg1++;
    }

    if (4 == dg1)
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE

/*
  Key 41 is C#N
*/

int
MACCSKeys::_key41 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (3 != a->nbonds())
      continue;

    atomic_number_t c = a->other(i, 0);
    if (6 == z[c])
      rc++;
  }

  return rc;
}
#endif

/*
  Key41 is O ~ [G1] ~ N
*/

int
MACCSKeys::_key41 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;

  const int istop = mpr.last_oxygen();

  for (int i = mpr.first_oxygen(); i < istop; ++i)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Bond * b = atoms[i]->item(0);

    const auto c = b->other(i);

    if (6 != z[c])
      continue;

    if (3 != ncon[c])
      continue;

    const auto ac = atoms[c];

    if (4 != ac->nbonds())
      continue;

    int nitrogen_found = 0;

    for (int j = 0; j < 3; ++j)
    {
      const auto n = ac->other(c, j);
      if (7 == z[n])
      {
        nitrogen_found = 1;
        break;
      }
    }

    if (nitrogen_found)
      rc++;
  }

  return rc;
}
/*
  Key 42 is F
*/

int
MACCSKeys::_key42 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (9 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  Key 43 is QHAQH
*/

int
MACCSKeys::_key43 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      const Atom * ak = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        atom_number_t q2 = ak->other(k, l);
        if (i == q2)
          continue;

        if (6 != z[q2] && q2 > i && implicit_hydrogens[q2])
          rc++;
      }
    }
  }

  return rc;
}

/*
  Key 44 is OTHER
  Oct 2007 change to [CG1R0]=!@[G1R0]
*/

int
MACCSKeys::_key44 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ring_membership[i])
      continue;

    const Atom * ai = atoms[i];

    int icon = ncon[i];

    if (icon + 1 != ai->nbonds())
      continue;

    for (int j = 0; j < icon; j++)
    {
      const Bond * bj = ai->item(j);

      if (! bj->is_double_bond())
        continue;

      atom_number_t k = bj->other(i);

      if (ring_membership[k])
        continue;

      if (ncon[k] + 1 != atoms[k]->nbonds())
        continue;

      rc++;
      break;
    }
  }

  return rc;
}

/*
  Key 45 is C=CN
*/

int
MACCSKeys::_key45 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  const int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    if (a->nbonds() != a->ncon() + 1)
      continue;

    int vinyl = 0;
    int nitrogen = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      Bond * b = a->item(j);

      atom_number_t k = b->other(i);

      if (6 == z[k])
      {
        if (b->is_double_bond() || b->is_aromatic())
          vinyl++;
      }
      else if (7 == z[k])
        nitrogen++;
    }

    if (vinyl && nitrogen)
      rc++;
  }

  return rc;
}

/*
  Key 46 is Br
*/

int
MACCSKeys::_key46 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (35 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  Key 47 is SAN
*/

int
MACCSKeys::_key47 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = ai->other(i, j);

      const Atom * ak = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        atom_number_t n = ak->other(k, l);
        if (7 == z[n])
          rc++;
      }
    }
  }

  return rc;
}

/*
  Key 48 is OQ(O)O
  Too rare, 
  Oct 2007 change to [O,S]=[S,C;R]@N
*/

int
MACCSKeys::_key48 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (1 != ncon[i])
      continue;

    if (8 == z[i] || 16 == z[i])
      ;
    else
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t c = b->other(i);

    if (3 != ncon[c])
      continue;

    if (0 == ring_membership[c])
      continue;

    if (6 == z[c] || 16 == z[c])
      ;
    else
      continue;

    const Atom * ac = atoms[c];

    for (int j = 0; j < 3; j++)
    {
      const Bond * bj = ac->item(j);

      if (! bj->is_single_bond())
        continue;

      if (0 == bj->nrings())
        continue;

      atom_number_t n = bj->other(c);

      if (7 == z[n])
        rc++;
    }
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 49 is charge
*/

int
MACCSKeys::_key49 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  int rc = m.number_formally_charged_atoms();

  return rc;
}
#endif

/*
  Two fully saturated heteroatoms attached to the same ring
*/

static int
count_saturated_heteroatoms_attached (Molecule & m,
                                      const Ring & r,
                                      const MK_Molecular_Properties & mpr)
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; ++i)
  {
    const atom_number_t j = r[i];

    if (2 == ncon[j])    // no bonds outside the ring
      continue;

    const auto aj = atoms[j];

    for (int k = 0; k < ncon[j]; ++k)
    {
      const Bond * b = aj->item(k);
      if (b->nrings())
        continue;

      const auto l = b->other(j);

      if (6 == z[l])
        continue;

      if (ring_membership[l])
        continue;

      if (ncon[l] < m.nbonds(l))
        continue;

      rc++;
    }
  }

  return rc > 1;
}

int
MACCSKeys::_key49 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    rc += count_saturated_heteroatoms_attached(m, *ri, mpr);
  }

  return rc;
}

/*
  Key 50 is C=C(C)C
*/

int
MACCSKeys::_key50 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * a = atoms[i];

    int vinyl = 0;
    int carbon = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      atom_number_t k = b->other(i);

      if (6 == z[k] && b->is_double_bond())
        vinyl++;
      else if (6 == z[k])
        carbon++;
    }
    
    if (vinyl && carbon >= 2)
      rc++;
  }

  return rc;
}

/*
  'a[C,S](=[O,S])Na'
*/

int
MACCSKeys::_key51 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  if (mpr.nrings() < 2)
    return 0;

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (1 != ncon[i])
      continue;

    if (8 == z[i] || 16 == z[i])
      ;
    else
      continue;

    const auto c = atoms[i]->other(i, 0);

    if (3 != ncon[c])
      continue;

    if (6 == z[c] || 16 == z[c])
      ;
    else
      continue;

    const auto ac = atoms[c];

    int found_aromatic = 0;
    int found_nitrogen = 0;

    for (auto j = 0; j < ncon[c]; ++j)
    {
      const auto k = ac->other(c, j);

      if (aromatic[k])
        found_aromatic++;
      else if (7 == z[k])
      {
        const auto ak = atoms[k];
        for (auto l = 0; l < ncon[k]; ++l)
        {
          const auto aa = ak->other(k, l);

          if (aa == c)
            continue;

          if (aromatic[aa])
          {
            found_nitrogen = 1;
            break;
          }
        }
      }
    }

    if (found_aromatic && found_nitrogen)
      rc++;
  }

  return rc;
}

/*
  Key 52 is NN
*/

int
MACCSKeys::_key52 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (7 == z[k] && k > i)
        rc++;
    }
  }

  return rc;
}

/*
  Key 53 is QHAAAQH
*/

int
MACCSKeys::_key53 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * ai = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);
      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        const atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          const atom_number_t a3 = atom_a2->other(a2, l);
          if (a3 == a1)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int n = 0; n < ncon[a3]; n++)
          {
            const atom_number_t q2 = atom_a3->other(a3, n);
             
            if (a1 == q2 || a2 == q2)
              continue;

            if (q2 <= i)
              continue;

            if (6 == z[q2])
              continue;

            if (0 == implicit_hydrogens[q2])
              continue;

            if (_aromatic_nitrogens_do_not_have_hydrogens && 7 == z[q2] && 
                ring_membership[q2] && aromatic[q2])
              continue;

            rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  Key 54 is QHAAQH
*/

int
MACCSKeys::_key54 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * ai = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    if (_aromatic_nitrogens_do_not_have_hydrogens && 7 == z[i] &&
        aromatic[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t q2 = atom_a2->other(a2, l);
          if (q2 == a1)
            continue;

          if (q2 <= i)
            continue;

          if (6 == z[q2])
            continue;

          if (0 == implicit_hydrogens[q2])
            continue;

          rc++;
        }
      }
    }
  }

  return rc;
}

/*
  spiro
*/

static atom_number_t
identify_spiro_fusion_atom (const Ring & r1, const Ring & r2)
{
  const int r1size = r1.number_elements();

  atom_number_t rc = INVALID_ATOM_NUMBER;

  for (int i = 0; i < r1size; ++i)
  {
    int j = r2.index(r1[i]);
    if (j < 0)
      continue;

    if (INVALID_ATOM_NUMBER == rc)
      rc = r2[j];
    else                        // have found multiple atoms in common, not spiro
      return 0;
  }

  return rc;
}

int
MACCSKeys::_key55 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;

  Set_of_Atoms spiro_fusion_points;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    for (auto j = i + 1; j < nr; ++j)
    {
      const auto rj = m.ringi(j);

      if (ri->fused_system_identifier() == rj->fused_system_identifier())
        continue;

      const auto s = identify_spiro_fusion_atom(*ri, *rj);

      if (INVALID_ATOM_NUMBER == s)
        continue;

      if (spiro_fusion_points.contains(s))
        continue;

      spiro_fusion_points.add(s);
      rc++;
    }
  }

  return rc;
}

/*
  Key 56 is ON(O)C
*/

int
MACCSKeys::_key56 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    int oxygen = 0;
    int carbon = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        carbon++;
      else if (8 == z[k])
        oxygen++;
    }
    
    if (carbon && oxygen >= 2)
      rc++;
  }

  return rc;
}

/*
  Key 57 is O heterocycle
*/

int
MACCSKeys::_key57 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  if (mpr.first_heteroatom() > istart)
    istart = mpr.first_heteroatom();

  int istop = mpr.last_ring_atom();
  if (mpr.last_heteroatom() < istop)
    istop = mpr.last_heteroatom();

  for (int i = istart; i < istop; i++)
  {
    if (8 == z[i] && ring_membership[i])
      rc++;
  }

  return rc;
}

/*
  Key 58 is QSQ
  Too correlated with others
  Change to a-[CD>1]
*/

int
MACCSKeys::_key58 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (! aromatic[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      const auto b = a->item(j);

      if (b->nrings())
        continue;

      const auto k = b->other(i);
      if (6 != z[k])
        continue;

      if (ncon[k] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  Key 59 is Snot%A%A - we interpret this as S attached to an aromatic
*/

int
MACCSKeys::_key59 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      
      if (b->is_aromatic())
        continue;

      atom_number_t k = b->other(i);

      if (aromatic[k])
        rc++;
    }
  }

  return rc;
}

/*
  S~O
*/

int
MACCSKeys::_key60 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const auto k = a->other(i, j);

      if (8 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  '[D2R0H2](-a)-[R]'
*/

int
MACCSKeys::_key61 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  if (m.nrings() < 2)
    return 0;

  const auto matoms = mpr.natoms();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (2 != ncon[i])
      continue;

    if (ring_membership[i])
      continue;

    const auto a = atoms[i];

    int got_aromatic = 0;
    int got_ring = 0;

    for (auto j = 0; j < 2; ++j)
    {
      const auto k = a->other(i, j);

      if (0 == ring_membership[k])
        continue;

      if (aromatic[k])
      {
        if (got_aromatic)
          got_ring = 1;
        else
          got_aromatic = 1;
      }
      else
        got_ring++;
    }

    if (got_aromatic && got_ring)
      rc++;
  }

  return rc;
}

/*
  Key 62 is A$A!A$A
  Two ring atoms, not in the same ring
*/

int
MACCSKeys::_key62 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();
  const int matoms = mpr.natoms();

  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * in_same_ring = mpr.in_same_ring();

  if (nr < 2)
    return 0;

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (0 == ring_membership[i])
      continue;

    if (2 == ncon[i])    // cannot have any bonds outside the ring
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == ring_membership[k])
        continue;

      if (0 == in_same_ring[i * matoms + k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 63 is N=O
*/

int
MACCSKeys::_key63 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  const int istop = mpr.last_nitrogen();
  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(i);
      if (8 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 64 is A$A!S
  which is chain S attached to a ring
*/

int
MACCSKeys::_key64 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  if (mpr.first_heteroatom() > istart)
    istart = mpr.first_heteroatom();

  int istop = mpr.last_heteroatom();
  if (mpr.last_heteroatom() < istop)
    istop = mpr.last_heteroatom();

  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    if (ring_membership[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (ring_membership[k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 65 is C%N - aromatic bond
*/

int
MACCSKeys::_key65 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  if (mpr.first_ring_atom() > istart)
    istart = mpr.first_ring_atom();

  int istop = mpr.last_heteroatom();
  if (mpr.last_ring_atom() < istop)
    istop = mpr.last_ring_atom();

  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (0 == ring_membership[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (! b->is_aromatic())
        continue;

      atom_number_t k = b->other(i);

      assert (mpr.aromatic()[k]);

      if (6 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  Key 66 is CC(C)(C)A
  A carbon, with at least four connections, at least three of which are C
*/

int
MACCSKeys::_key66 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (4 != ncon[i])
      continue;

    const Atom * a = atoms[i];
    if (4 != a->nbonds())
      continue;

    int carbon = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        carbon++;
    }
    
    if (carbon >= 3)
      rc++;
  }

  return rc;
}

/*
  6 membered, para substituted aromatic, in the scaffold
*/

int
MACCSKeys::_is_para_substituted_within_scaffold (const Molecule & m,
                                             const Set_of_Atoms & r,
                                             const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * spinach = mpr.spinach();
  const Atom * const * atoms = mpr.atoms();

  resizable_array<int> scaffold_connection;

  for (auto i = 0; i < 6; ++i)
  {
    const auto j = r[i];

    if (2 == ncon[j])
      continue;

    const auto aj = atoms[j];

    for (auto k = 0; k < 3; ++k)
    {
      const auto b = aj->item(k);

      if (b->nrings())
        continue;

      const auto l = b->other(j);

      if (! spinach[l])
        scaffold_connection.add(i);
    }
  }

  const auto n = scaffold_connection.number_elements();

  if (2 != n)
    return 0;

  if (scaffold_connection[0] + 3 == scaffold_connection[1])   // hopefully the most common case
    return 1;
  else
    return 0;
}

int
MACCSKeys::_key67 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (6 != ri->number_elements())
      continue;

    if (_is_para_substituted_within_scaffold(m, *ri, mpr))
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Key 68 is QHQH (&...)
*/

int
MACCSKeys::_key68 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (k < i)
        continue;

      if (6 != z[k] && implicit_hydrogens[k])
        rc++;
    }
  }

  return rc;
}
#endif

/*
  Key 68 is 'a!@[CT0]!@[CT0]'
*/

int
MACCSKeys::_key68 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const auto matoms = mpr.natoms();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_bond_count = mpr.ring_bond_count();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();
  const int * aromatic = mpr.aromatic();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (2 != ring_bond_count[i])
      continue;

    const auto ai = atoms[i];

    atom_number_t c1 = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; ++j)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      const atom_number_t x = b->other(i);

      if (6 != z[x])
        continue;

      if (attached_heteroatom_count[x] > 0)
        continue;

      c1 = x;
      break;
    }

    if (INVALID_ATOM_NUMBER == c1)
      continue;

    const auto ac1 = atoms[c1];

    const auto nconc1 = ncon[c1];

    for (int j = 0; j < nconc1; ++j)
    {
      const Bond * b = ac1->item(j);

      if (b->nrings())
        continue;

      const atom_number_t c2 = b->other(c1);

      if (c2 == i)
        continue;

      if (6 != z[c2])
        continue;

      if (attached_heteroatom_count[c2] > 0)
        continue;

      rc++;     // deliberately no break so we can count them - possible aC(C)C
    }
  }

  return rc;
}

/*
  key69 is QQH

  this is one case where the definition is open to interpretation
*/

int
MACCSKeys::_key69 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        continue;

      if (0 == implicit_hydrogens[k])
        rc++;
      else if (k > i)
        rc++;
    }
  }

  return rc;
}

/*
  key70 is QNQ
*/

int
MACCSKeys::_key70 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    int heteroatoms = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        heteroatoms++;
    }
    
    if (heteroatoms >= 2)
      rc++;
  }

  return rc;
}

/*
  key71 is NO
*/

int
MACCSKeys::_key71 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (8 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key72 is OAAO
*/

int
MACCSKeys::_key72 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];
        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t o2 = atom_a2->other(a2, l);
          if (a1 == o2)
            continue;

          if (8 == z[o2] && o2 > i)
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key73 is S=
*/

int
MACCSKeys::_key73 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (b->is_double_bond())
        rc++;
    }
  }

  return rc;
}

/*
  key74 is CH3ACH3
*/

int
MACCSKeys::_key74 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (implicit_hydrogens[i] < 3)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t c2 = atom_a1->other(a1, k);
        if (c2 <= i)
          continue;

        if (6 != z[c2])
          continue;

        if (1 != ncon[i])
          continue;

        if (implicit_hydrogens[c2] >= 3)
          rc++;
      }
    }
  }

  return rc;
}

/*
  key75 is A!N$A which is Ring nitrogen with a branch outside the ring
*/

int
MACCSKeys::_key75 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * in_same_ring = mpr.in_same_ring();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  if (mpr.first_ring_atom() > istart)
    istart = mpr.first_ring_atom();

  int istop = mpr.last_heteroatom();
  if (mpr.last_ring_atom() < istop)
    istop = mpr.last_ring_atom();

  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (0 == ring_membership[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == ring_membership[k] || 0 == in_same_ring[i * matoms + k])
        rc++;
    }
  }

  return rc;
}

/*
  key76 is C=C(A)A
*/

int
MACCSKeys::_key76 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  const int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * a = atoms[i];
    if (4 != a->nbonds())
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (b->is_double_bond())
        ;
      else if (b->is_aromatic())
        ;
      else
        continue;

      atom_number_t k = b->other(i);

      if (6 == z[k])
      {
        rc++;
        break;
      }
    }
  }

  return rc;
}

/*
  key77 is NAN
*/

int
MACCSKeys::_key77 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t n = atom_a1->other(a1, k);
        if (n == i)
          continue;

        if (7 == z[n])
          rc++;
      }
    }
  }

  return rc;
}

/*
  key78 is C=N
*/

int
MACCSKeys::_key78 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    if (a->ncon() + 1 != a->nbonds())
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_double_bond())
        ;
      else if (b->is_aromatic())
        ;
      else
        continue;

      atom_number_t k = b->other(i);

      if (6 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key79 is NAAN
*/

int
MACCSKeys::_key79 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
     continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t n = atom_a2->other(a2, l);
          if (n == a1)
            continue;

          if (7 == z[n] && n > i)
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key80 is NAAAN
*/

int
MACCSKeys::_key80 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
     continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = ai->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a1 == a3)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int q = 0; q < ncon[a3]; q++)
          {
            atom_number_t n = atom_a3->other(a3, q);
            if (n == a2)
              continue;

            if (7 == z[n] && n > i)
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key81 is SA(A)A
  Which is a Sulphur bonded to an atom with at least three connections
*/

int
MACCSKeys::_key81 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (ncon[k] >= 3)
        rc++;
    }
  }

  return rc;
}

/*
  Ester 'O=[CT2D3]-[OD2T0]'
*/

int
MACCSKeys::_key82 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();

  int rc = 0;
  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (1 != ncon[i])
      continue;

    if (8 == z[i] || 16 == z[i])
      ;
    else
      continue;

    if (2 != atoms[i]->nbonds())
      continue;

    const auto c = atoms[i]->other(i, 0);

    if (6 != z[c])
      continue;

    if (3 != ncon[c])
      continue;

    const auto ac = atoms[c];

    int got_carbon = 0;
    int got_oxygen = 0;

    for (auto j = 0; j < 3; ++j)
    {
      const Bond * b = ac->item(j);

      if (! b->is_single_bond())
        continue;

      const auto k = b->other(c);

      if (6 == z[k])
        got_carbon = 1;
      else if ((8 == z[k] || 16 == z[k]) && 0 == attached_heteroatom_count[k])
        got_oxygen = 1;
    }

    if (got_oxygen && got_carbon)
      rc++;
  }

  return rc;
}

static int
is_heterocycle (const Ring * r,
                const atomic_number_t * z)
{
  int ring_size = r->number_elements();
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t a = r->item(i);
    if (6 != z[a])
      return 1;
  }

  return 0;
}

/*
  key83 is QAAAA@1 which is five membered heterocycle
*/

int
MACCSKeys::_key83 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (5 != r->number_elements())
      continue;

    if (is_heterocycle(r, z))
      rc++;
  }

  nr = m.non_sssr_rings();
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.non_sssr_ring(i);
    if (5 != r->number_elements())
      continue;

    if (is_heterocycle(r, z))
      rc++;
  }

  return rc;
}

/*
  key84 is NH2
*/

int
MACCSKeys::_key84 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 == z[i] && implicit_hydrogens[i] >= 2)
    {
      rc++;
    }
  }

  return rc;
}

/*
  key85 is CN(C)C
*/

int
MACCSKeys::_key85 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    int carbon = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        carbon++;
    }

    if (carbon >= 3)
      rc++;
  }

  return rc;
}

/*
  key86 is CH2QCH2
*/

int
MACCSKeys::_key86 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    int ch2 = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k] && implicit_hydrogens[k] >= 2)
        ch2++;
    }
    
    if (ch2 >= 2)
      rc++;
  }

  return rc;
}

/*
  '[N,O,S;D1]=[CRD3]-[ND2R]'
*/

int
MACCSKeys::_key87 (Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  if (0 == mpr.nrings())
    return 0;

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;

  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const auto b = atoms[i]->item(0);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(i);

    if (6 != z[j])
      continue;

    if (3 != ncon[j])
      continue;

    if (0 == ring_membership[j])
      continue;

    const auto aj = atoms[j];

    for (auto k = 0; k < 3; ++k)
    {
      const auto b = aj->item(k);

      if (! b->is_single_bond())
        continue;

      const auto n = b->other(j);

      if (7 != z[n])
        continue;

      if (2 != ncon[n])
        continue;

      if (0 == ring_membership[n])
        continue;

      rc++;
      break;
    }
  }

  return rc;
}

/*
  key88 is S
*/

int
MACCSKeys::_key88 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (16 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key89 is OAAAO
*/

int
MACCSKeys::_key89 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a1 == a3)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int q = 0; q < ncon[a3]; q++)
          {
            atom_number_t o = atom_a3->other(a3, q);
            if (o == a2)
              continue;

            if (8 == z[o] && o > i)
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key90 is QHAACH2A
*/

int
MACCSKeys::_key90 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (a2 == i)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t ch2 = atom_a2->other(a2, l);
          if (ch2 == a1)
            continue;

          if (6 == z[ch2] && implicit_hydrogens[ch2] >= 2 && ncon[ch2] > 1)
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key91 is QHAAACH2A
*/

int
MACCSKeys::_key91 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (a2 == i)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a3 == a1)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int q = 0; q < ncon[a3]; q++)
          {
            atom_number_t ch2 = atom_a3->other(a3, q);
            if (ch2 == a2)
              continue;

            if (6 == z[ch2] && implicit_hydrogens[ch2] >= 2 &&
                ncon[ch2] > 1)
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key92 is OC(N)C
*/

int
MACCSKeys::_key92 (const Molecule & m,
                   const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    int carbon = 0;
    int oxygen = 0;
    int nitrogen = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k])
        carbon++;
      else if (7 == z[k])
        nitrogen++;
      else if (8 == z[k])
        oxygen++;
    }
    
    if (carbon && nitrogen && oxygen)
      rc++;
  }

  return rc;
}

/*
  key93 is QCH3
*/

int
MACCSKeys::_key93 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k] && implicit_hydrogens[k] >= 3)
        rc++;
    }
  }

  return rc;
}

/*
  key94 is QN
*/

int
MACCSKeys::_key94 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (7 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key95 is NAAO
*/

int
MACCSKeys::_key95 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
     continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t o = atom_a2->other(a2, l);
          if (o == a2)
            continue;

          if (8 == z[o])
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key96 is 5 membered ring
*/

int
MACCSKeys::_key96 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (5 == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  key97 is NAAAO - we process as OAAAN
*/

int
MACCSKeys::_key97 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (i == a2)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a1 == a3)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int q = 0; q < ncon[a3]; q++)
          {
            atom_number_t n = atom_a3->other(a3, q);
            if (n == a2)
              continue;

            if (7 == z[n])
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key98 is QAAAAA@1, or 6 membered heterocycle
*/

int
MACCSKeys::_key98 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (6 == r->number_elements() && is_heterocycle(r, z))
      rc++;
  }

  nr = m.non_sssr_rings();
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.non_sssr_ring(i);
    if (6 == r->number_elements() && is_heterocycle(r, z))
      rc++;
  }

  return rc;
}

/*
  key99 is C=C
*/

int
MACCSKeys::_key99 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;

  const int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_aromatic())
        continue;

      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(i);

      if (k < i && 6 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key100 is ACH2N
*/

int
MACCSKeys::_key100 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 == z[k] && implicit_hydrogens[k] >= 2 && ncon[k] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  key101 is 8 membered ring or larger
  Too rare, change to 7 or larger
*/

int
MACCSKeys::_key101 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (r->number_elements() >= 7)    // changed
      rc++;
  }

  return rc;
}

/*
  key102 is QO
*/

int
MACCSKeys::_key102 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_oxygen();
  for (int i = mpr.first_oxygen(); i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key103 is CL
*/

int
MACCSKeys::_key103 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (17 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key104 is QHACH2A
*/

int
MACCSKeys::_key104 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

    const Atom * a = atoms[i];

    if (0 == implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t ch2 = atom_a1->other(a1, k);
        if (ch2 == i)
          continue;

        if (6 != z[ch2])
          continue;

        if (implicit_hydrogens[ch2] >= 2 && ncon[ch2] > 1)
          rc++;
      }
    }
  }

  return rc;
}

/*
  key105 is A$A($A)$A
    or an atom in 3 rings
*/

int
MACCSKeys::_key105 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();
  const int matoms = mpr.natoms();

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * in_same_ring = mpr.in_same_ring();

  if (nr < 2)
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  int istop  = mpr.last_ring_atom();
  for (int i = istart; i < istop; i++)
  {
    if (ring_membership[i] < 2)
      continue;

    const Atom * a = atoms[i];

    int ring_bonds = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == ring_membership[k])
        continue;

      if (in_same_ring[i * matoms + k])
        ring_bonds++;
    }

    if (3 == ring_bonds)
      rc++;
  }

  return rc;
}

/*
  key106 is QA(Q)Q
*/

int
MACCSKeys::_key106 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = atoms[i];

    int heteroatoms = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        heteroatoms++;
    }

    if (heteroatoms >= 3)
      rc++;
  }

  return rc;
}

/*
  aromatic attached halogen
*/

int
MACCSKeys::_key107 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (1 != ncon[i])
      continue;

    if (6 == z[i] || 7 == z[i])      // efficiency
      continue;

    if (9 == z[i] || 17 == z[i] || 35 == z[i] || 53 == z[i] || 85 == z[i])
    {
      const Atom * a = atoms[i];

      for (int j = 0; j < ncon[i]; j++)
      {
        atom_number_t k = a->other(i, j);
        if (aromatic[k])
          rc++;
      }
    }
  }

  return rc;
}

/*
  key108 is CH3AAACH2A
*/

int
MACCSKeys::_key108 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    const Atom * a = atoms[i];

    if (implicit_hydrogens[i] < 3)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (a2 == i)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a3 == a1)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int q = 0; q < ncon[a3]; q++)
          {
            atom_number_t ch2 = atom_a3->other(a3, q);
            if (ch2 == a2)
              continue;

            if (6 != z[ch2])
              continue;

            if (implicit_hydrogens[ch2] < 2)
              continue;

            if (ncon[ch2] > 1)
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key109 is ACH2O
*/

int
MACCSKeys::_key109 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        continue;

      if (implicit_hydrogens[k] < 2)
        continue;

      if (ncon[k] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  key110 is NCO
*/

int
MACCSKeys::_key110 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t c = a->other(i, j);
      if (6 != z[c])
        continue;

      const Atom * atom_c = atoms[c];

      for (int k = 0; k < ncon[c]; k++)
      {
        atom_number_t o = atom_c->other(c, k);
        if (8 == z[o])
          rc++;
      }
    }
  }

  return rc;
}

/*
  key111 is NACH2A
*/

int
MACCSKeys::_key111 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t ch2 = atom_a1->other(a1, k);
        if (6 != z[ch2])
          continue;

        if (implicit_hydrogens[ch2] < 2)
          continue;

        if (ncon[ch2] > 1)
          rc++;
      }
    }
  }

  return rc;
}

/*
  key112 is AA(A)(A)A
  which is a four connected atom
*/

int
MACCSKeys::_key112 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (4 == ncon[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key113 is Onot%A%A
  or an oxygen bonded via a non-aromatic bond to an aromatic ring
*/

int
MACCSKeys::_key113 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  if (0 == nr)
    return 0;

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;

  int istop = mpr.last_oxygen();
  for (int i = mpr.first_oxygen(); i < istop; i++)
  {
    if (8 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_aromatic())
        continue;

      atom_number_t k = b->other(i);

      if (aromatic[k])
        rc++;
    }
  }

  return rc;
}

/*
  key114 is CH3CH2A
*/

int
MACCSKeys::_key114 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

     const Atom *a = atoms[i];

    if (implicit_hydrogens[i] < 3)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        continue;

      if (implicit_hydrogens[k] < 2)
        continue;

      if (ncon[k] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  key115 is CH3ACH2A
*/

int
MACCSKeys::_key115 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (1 != ncon[i])     // this part had been omitted
      continue;

    if (1 != atoms[i]->nbonds())
      continue;

    atom_number_t A = atoms[i]->other(i, 0);

     const Atom *a = atoms[A];

    for (int j = 0; j < ncon[A]; j++)
    {
      atom_number_t ch2 = a->other(A, j);

      if (i == ch2)
        continue;

      if (6 != z[ch2])
        continue;

      if (2 != implicit_hydrogens[ch2])
        continue;

      if (ncon[ch2] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  key116 is CH3AACH2A
*/

int
MACCSKeys::_key116 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (6 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

     const Atom *a = atoms[i];

    if (3 != implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

       const Atom *atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);
        if (a2 == i)
          continue;

         const Atom *atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t ch2 = atom_a2->other(a2, l);

          if (ch2 == a1)
            continue;

          if (6 != z[ch2])
            continue;
  
          if (implicit_hydrogens[ch2] < 2)
            continue;
  
          if (ncon[ch2] > 1)
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key117 is NAO
*/

int
MACCSKeys::_key117 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

       const Atom *atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t o = atom_a1->other(a1, k);
        if (o == i)
          continue;

        if (8 == z[o])
          rc++;
      }
    }
  }

  return rc;
}

/*
  key118 is ACH2CH2A > 1
  We just produce the count
*/

int
MACCSKeys::_key118 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

     const Atom *a = atoms[i];

    if (implicit_hydrogens[i] < 2)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (k < i)
        continue;

      if (6 != z[k])
        continue;

      if (ncon[k] < 2)
        continue;

      if (implicit_hydrogens[k] >= 2)
        rc++;
    }
  }

  return rc;
}

/*
  key119 is N=A
  Note that we intentionally double count N=N
*/

int
MACCSKeys::_key119 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_double_bond())
        rc++;
    }
  }

  return rc;
}

/*
  key120 is Heterocyclic atom > 1
  We just count them
*/

int
MACCSKeys::_key120 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  if (mpr.first_ring_atom() > istart)
    istart = mpr.first_ring_atom();

  int istop = mpr.last_heteroatom();
  if (mpr.last_ring_atom() < istop)
    istop = mpr.last_ring_atom();

  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i] && ring_membership[i])
      rc++;
  }

  return rc;
}

/*
  key121 is N heterocycle
*/

int
MACCSKeys::_key121 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  if (mpr.first_ring_atom() > istart)
    istart = mpr.first_ring_atom();

  int istop = mpr.last_heteroatom();
  if (mpr.last_ring_atom() < istop)
    istop = mpr.last_ring_atom();

  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (ring_membership[i])
      rc++;
  }

  return rc;
}

/*
  'a!@[R]'
*/

int
MACCSKeys::_key122 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 2)
    return 0;

  const auto matoms = mpr.natoms();

  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (0 == aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const auto ai = atoms[i];

    for (auto j = 0; j < 3; ++j)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      const auto k = b->other(i);

      if (ring_membership[k])
        rc++;
    }
  }

  return rc;
}

/*
  key123 is OCO
*/

int
MACCSKeys::_key123 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t c = a->other(i, j);
      if (6 != z[c])
        continue;

       const Atom *atom_c = atoms[c];

      for (int k = 0; k < ncon[c]; k++)
      {
        atom_number_t o2 = atom_c->other(c, k);
        if (o2 == i)
          continue;

        if (o2 < i)
          continue;

        if (8 == z[o2])
          rc++;
      }
    }
  }

  return rc;
}

/*
  key124 is QQ
*/

int
MACCSKeys::_key124 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (k > i && 6 != z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key125 is Arom ring > 1
    We just count them

  Change to number of five membered aromatic rings
*/

int
MACCSKeys::_key125 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (! r->is_aromatic())
      continue;

    if (5 == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  key126 is A!O!A
*/

int
MACCSKeys::_key126 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();
  const int matoms = mpr.natoms();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * in_same_ring = mpr.in_same_ring();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

     const Atom *a = atoms[i];

    int chain_bonds = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == nr || 0 == in_same_ring[i * matoms + k])   // careful, in_same_ring may be null
        chain_bonds++;
    }

    if (chain_bonds >= 2)
      rc++;
  }

  return rc;
}

/*
  '[G>0R0]-a([aD2])[aD2]' 
*/

int
MACCSKeys::_key127 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  if (0 == mpr.nrings())
    return 0;

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();


  const auto matoms = mpr.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i])
      continue;

    const auto ai = atoms[i];

    if (ncon[i] == ai->nbonds())    // fully saturated
      continue;

    for (auto j = 0; j < ncon[i]; ++j)
    {
      const auto b = ai->item(j);

      if (! b->is_single_bond())
        continue;

      const auto k = b->other(i);

      if (! aromatic[k])
        continue;

      if (3 != ncon[k])
        continue;

      const auto ak = atoms[k];

      int aD2 = 0;
      for (auto l = 0; l < 3; ++l)
      {
        const auto n = ak->other(k, l);

        if (n == i)
          continue;

        if (2 == ncon[n])
          aD2++;
      }

      if (2 == aD2)
        rc++;
     }
  }

  return rc;
}

/*
  key128 is ACH2AAACH2A
*/

int
MACCSKeys::_key128 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

    if (2 != implicit_hydrogens[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

      const Atom * atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);

        if (a2 == i)
          continue;

        const Atom * atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t a3 = atom_a2->other(a2, l);
          if (a3 == a1 || a3 == i)
            continue;

          const Atom * atom_a3 = atoms[a3];

          for (int n = 0; n < ncon[a3]; n++)
          {
            atom_number_t ch2 = atom_a3->other(a3, n);
            if (ch2 == a2 || ch2 == a1 || ch2 == i)
              continue;

            if (ch2 < i)
              continue;

            if (6 != z[ch2])
              continue;

            if (ncon[ch2] < 2)
              continue;

            if (2 == implicit_hydrogens[ch2])
              rc++;
          }
        }
      }
    }
  }

  return rc;
}

/*
  key129 is ACH2AACH2A
*/

int
MACCSKeys::_key129 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

     const Atom *a = atoms[i];

    if (2 != implicit_hydrogens[i])
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

       const Atom *atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t a2 = atom_a1->other(a1, k);

        if (a2 == i)
          continue;

         const Atom *atom_a2 = atoms[a2];

        for (int l = 0; l < ncon[a2]; l++)
        {
          atom_number_t ch2 = atom_a2->other(a2, l);
          if (ch2 == a1 || ch2 == i)
            continue;

          if (ch2 < i)
            continue;

          if (6 != z[ch2])
            continue;

          if (ncon[ch2] < 2)
            continue;

          if (2 == implicit_hydrogens[ch2])
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  key130 is QQ > 1 (&..)
    this is the same as bit 124
  Change to c-!@c
*/

int
MACCSKeys::_key130 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (k > i)
        continue;

      if (6 != z[k])
        continue;

      if (! aromatic[k])
        continue;

      rc++;
      break;
    }
  }

  return rc;
}

/*
  key131 is QH > 1
    we just count them
*/

int
MACCSKeys::_key131 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i] && implicit_hydrogens[i])
      rc++;
  }

  return rc;
}

/*
  key132 is OACH2A
*/

int
MACCSKeys::_key132 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t a1 = a->other(i, j);

       const Atom *atom_a1 = atoms[a1];

      for (int k = 0; k < ncon[a1]; k++)
      {
        atom_number_t ch2 = atom_a1->other(a1, k);
        if (6 != z[ch2])
          continue;

        if (implicit_hydrogens[ch2] >= 2)
          rc++;
      }
    }
  }

  return rc;
}

/*
  key133 is A$A!N
    A nitrogen attached to an aromatic ring
*/

int
MACCSKeys::_key133 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      
      if (b->is_aromatic())
        continue;

      atom_number_t k = b->other(i);

      if (aromatic[k])
        rc++;
    }
  }

  return rc;
}

/*
  key134 is X (HALOGEN)
*/

int
MACCSKeys::_key134 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (z[i] < 9)
      continue;

    if (9 == z[i] || 17 == z[i] || 35 == z[i] || 53 == z[i] || 85 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  '[cD3](!@*)[cD3](!@*)[cD3]!@*'
    which is three adjacent aromatic atoms, all substituted outside the ring
*/

int
MACCSKeys::_key135 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  if (0 == nr)
    return 0;

  Set_of_Atoms cD3;    // first we identify the possible atoms

  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->nrings())
        continue;

      cD3.add(i);      // aromatic with exocyclic bond
      break;
    }
  }

// Now we need to see if there are any groups of three that are bonded

  const auto n = cD3.number_elements();

  if (n < 3)
    return 0;

  int rc = 0;

  for (auto i = 0; i < n; ++i)
  {
    const auto ai = cD3[i];
    for (auto j = 0; j < n; ++j)
    {
      if (j == i)
        continue;

      const auto aj = cD3[j];

      if (! m.are_bonded(ai, aj))
        continue;

      for (auto k = j + 1; k < n; ++k)
      {
        if (k == i)
          continue;

        const auto ak = cD3[k];

        if (m.are_bonded(ai, ak))
          rc++;
      }
    }
  }

  return rc;
}

/*
  key136 is O=A > 1
    we just count them
*/

int
MACCSKeys::_key136 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
 const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    if (1 != a->ncon())
      continue;

    const Bond * b = a->item(0);
    if (b->is_double_bond())
        rc++;
  }

  return rc;
}

/*
  key137 is Heterocycle
*/

int
MACCSKeys::_key137 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (is_heterocycle(r, z))
      rc++;
  }

  nr = m.non_sssr_rings();
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.non_sssr_ring(i);
    if (is_heterocycle(r, z))
      rc++;
  }

// Because of the problems counting rings (especially macrocycles) we only
// allow 0 and 1 values for this bit

  if (rc > 0)
    rc = 1;

  return rc;
}

/*
  key138 is QCH2A > 1
    we just count
*/

int
MACCSKeys::_key138 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (6 != z[k])
        continue;

      if (implicit_hydrogens[k] < 2)
        continue;

      if (ncon[k] > 1)
        rc++;
    }
  }

  return rc;
}

/*
  key139 is OH
*/

int
MACCSKeys::_key139 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 == z[i] && implicit_hydrogens[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key140 is O > 3 (&..)
*/

int
MACCSKeys::_key140 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key141 is CH3 > 2
    we just count them

  Since this is the same as key 149, we do something different for
  key 149 - see below
*/

int
MACCSKeys::_key141 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 == z[i] && implicit_hydrogens[i] >= 3)
    {
      rc++;
    }
  }

  return rc;
}

/*
  key142 is N > 1
    we just count them
*/

int
MACCSKeys::_key142 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 == z[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key143 is A$A!O
    which is an Oxygen connected to a ring via a non-ring bond
*/

int
MACCSKeys::_key143 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * in_same_ring = mpr.in_same_ring();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (0 == ring_membership[k])
        continue;

      if (0 == in_same_ring[i * matoms + k])
        rc++;
    }
  }

  return rc;
}

/*
  key144 is Anot%A%Anot%A
    which is an ortho substituted aromatic ring, with non aromatic bonds outside the ring
*/

int
MACCSKeys::_key144 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  int istop  = mpr.last_ring_atom();
  for (int i = istart; i < istop; i++)
  {
    if (0 == aromatic[i])
      continue;

    if (2 == ncon[i])    // must have at least 3 connections
      continue;

    const Atom *a = atoms[i];

    atom_number_t non_arom_outside_ring = INVALID_ATOM_NUMBER;
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (aromatic[k])
        continue;

      non_arom_outside_ring = k;
      break;
    }

    if (INVALID_ATOM_NUMBER == non_arom_outside_ring)
      continue;

//  Now identify the next atom in the aromatic ring

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (! b->is_aromatic())
        continue;

      atom_number_t k = b->other(i);

      assert (aromatic[k]);

      if (k < i)
        continue;

      if (ncon[k] <= 2)
        continue;

      const Atom *atom_k = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        const Bond * b = atom_k->item(l);

        if (b->nrings())
          continue;

        const atom_number_t x = b->other(k);

        if (aromatic[x])
          continue;

        rc++;
      }
    }
  }

  return rc;
}

/*
  key145 is 6M ring > 1
    we just count them
*/

int
MACCSKeys::_key145 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (6 == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  key146 is O > 2, this is the same as key 140
  change to be a-[R0]-a
*/

int
MACCSKeys::_key146 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const Atom * const * atoms = mpr.atoms();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;
  int istop = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (! aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    int found_match = 0;

    for (int j = 0; j < 3; j++)
    {
      atom_number_t r0 = ai->other(i, j);

      if (ring_membership[r0])
        continue;

      int r0con = ncon[r0];
      if (1 == r0con)
        continue;

      const Atom * ar0 = atoms[r0];

      for (int k = 0; k < r0con; k++)
      {
        atom_number_t c2 = ar0->other(r0, k);

        if (c2 >= i)   // avoid double counting and doubling back
          continue;

        if (aromatic[c2])
        {
          found_match = 1;
          rc++;
          break;
        }
      }

      if (found_match)
        break;
    }
  }

  return rc;
}

/*
  '[CR0]([R])[R]'
*/

int
MACCSKeys::_key147 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  if (m.nrings() < 2)
    return 0;

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (ring_membership[i])
      continue;

    if (ncon[i] < 2)
      continue;

    const Atom *a = atoms[i];

    int ring_atoms = 0;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (ring_membership[k])
        ring_atoms++;
    }

    if (ring_atoms > 1)
      rc++;
  }

  return rc;
}

/*
  key148 is AQ(A)A which is a 3 connected heteroatom
*/

int
MACCSKeys::_key148 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 != z[i] && ncon[i] >= 3)
    {
      rc++;
    }
  }

  return rc;
}

/*
  key149 is CH3 > 1
    which is the same as key 141

  Since it is the same as key 141, we define it to be fully
  saturated, 4 connected carbon atoms - except things like CF3
*/

int
MACCSKeys::_key149 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();
  const int * ncon = mpr.ncon();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();

  int rc = 0;
  int istop = mpr.last_carbon();
  for (int i = mpr.first_carbon(); i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (implicit_hydrogens[i] > 0)
      continue;

    if (4 != ncon[i])
      continue;

    if (attached_heteroatom_count[i] >= 3)
      continue;

    rc++;
  }

  return rc;
}

/*
  key150 is A!A$A!A
    which is an ortho substituted ring
*/

int
MACCSKeys::_key150 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const int nr = mpr.nrings();

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * in_same_ring = mpr.in_same_ring();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istart = mpr.first_ring_atom();
  int istop  = mpr.last_ring_atom();
  for (int i = istart; i < istop; i++)
  {
    if (0 == ring_membership[i])
      continue;

    if (ncon[i] < 3)
      continue;

     const Atom *a = atoms[i];

    atom_number_t outside_ring = INVALID_ATOM_NUMBER;
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (0 == ring_membership[k] || 0 == in_same_ring[i * matoms + k])
      {
        outside_ring = k;
        break;
      }
    }

    if (INVALID_ATOM_NUMBER == outside_ring)
      continue;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (k < i)
        continue;

      if (ncon[k] < 3)
        continue;

      if (0 == ring_membership[k] || 0 == in_same_ring[i * matoms + k])
        continue;

       const Atom *atom_k = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        atom_number_t n = atom_k->other(k, l);
        if (n == i)
          continue;

        if (0 == ring_membership[n] || 0 == in_same_ring[k * matoms + n])
        {
          rc++;
          break;
        }
      }
    }
  }

  return rc;
}

/*
  key151 is NH
*/

int
MACCSKeys::_key151 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 == z[i] && implicit_hydrogens[i])
    {
      rc++;
    }
  }

  return rc;
}

/*
  key152 is OC(C)C
*/

int
MACCSKeys::_key152 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t c = a->other(i, j);
      if (6 != z[c])
        continue;

       const Atom *atom_c = atoms[c];

      int carbons = 0;
      for (int k = 0; k < ncon[c]; k++)
      {
        atom_number_t l = atom_c->other(c, k);
        if (6 == z[l])
          carbons++;
      }

      if (carbons >= 2)
        rc++;
    }
  }

  return rc;
}

/*
  '[#7]!@*[#7]'
*/

int
MACCSKeys::_key153 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istop = mpr.last_nitrogen();
  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const auto * b = a->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (ncon[k] < 2)
        continue;

      const auto ak = atoms[k];

      for (auto l = 0; l < ncon[k]; ++l)
      {
        const auto n = ak->other(k, l);

        if (7 == z[n])    // will do some level of double counting
          rc++;
      }
    }
  }

  return rc;
}

/*
  key154 is C=O
*/

int
MACCSKeys::_key154 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    if (1 != a->ncon())
      continue;

    const Bond * b = a->item(0);
    if (! b->is_double_bond())
      continue;

    atom_number_t k = b->other(i);

    if (6 == z[k])
      rc++;
  }

  return rc;
}

/*
  key155 is A!CH2!A
    which is a CH2 with two chain bonds
*/

int
MACCSKeys::_key155 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const int nr = mpr.nrings();

  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * in_same_ring = mpr.in_same_ring();
  const int * implicit_hydrogens = mpr.implicit_hydrogens();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (implicit_hydrogens[i] < 2)
      continue;

    int chain_bonds = 0;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (0 == nr || 0 == in_same_ring[i * matoms + k])   // careful, in_same_ring may be null
        chain_bonds++;
    }
    
    if (chain_bonds >= 2)
      rc++;
  }

  return rc;
}

/*
  key156 is NA(A)A
    which is a Nitrogen attached to a 3 connected atom
*/

int
MACCSKeys::_key156 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

     const Atom *a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (ncon[k] >= 3)
        rc++;
    }
  }

  return rc;
}

/*
  key157 is C-O
*/

int
MACCSKeys::_key157 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (8 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);
      if (6 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key158 is C-N
*/

int
MACCSKeys::_key158 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  int istart = mpr.first_heteroatom();
  int istop = mpr.last_heteroatom();
  for (int i = istart; i < istop; i++)
  {
    if (7 != z[i])
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (nullptr == aromatic)
        ;
      else if (aromatic[k])
        continue;

      if (6 == z[k])
        rc++;
    }
  }

  return rc;
}

/*
  key159 is O > 1
    which is the same as key 140
  Change to acidic group
*/

int
MACCSKeys::_key159 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(i);

    int jcon = ncon[j];

    if (jcon < 3)
      continue;

    if (6 == z[j] || 16 == z[j] || 15 == z[j])
      ;
    else
      continue;

    const Atom * aj = atoms[j];

    if (jcon == aj->nbonds())
      continue;

    int found_doubly_bonded_oxygen = 0;

    for (int k = 0; k < jcon; k++)
    {
      const Bond * b = aj->item(k);

      if (! b->is_double_bond())
        continue;

      atom_number_t o = b->other(j);

      if (8 == z[o] || 16 == z[o])
      {
        found_doubly_bonded_oxygen = 1;
        break;
      }
    }

    if (found_doubly_bonded_oxygen)
      rc++;
  }

  return rc;
}

/*
  key160 is CH3
    which is the same as key 149. 
    Key 149 was redefined, so let's do somthing slightly different here.
    fully saturated Carbon, >2 connections, in a ring
*/

int
MACCSKeys::_key160 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  int istart = mpr.first_carbon();
  int istop = mpr.last_carbon();
  for (int i = istart; i < istop; i++)
  {
    if (6 != z[i])
      continue;

    if (0 == ring_membership[i])
      continue;

    if (2 == ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (a->nbonds() != ncon[i])
      continue;

    rc++;
  }

  return rc;
}

/*
  key161 is N
    which is the same as key 142

  We change the definition to a-[R0]-[R0]-a

  variables
    c1 r01 r02 c2
*/

int
MACCSKeys::_key161 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * aromatic = mpr.aromatic();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  int istop = mpr.last_ring_atom();
  for (int c1 = mpr.first_ring_atom(); c1 < istop; c1++)
  {
    if (! aromatic[c1])
      continue;

    if (3 != ncon[c1])
      continue;

    const Atom * ai = atoms[c1];

    for (int j = 0; j < 3; j++)
    {
      const Bond * bj = ai->item(j);

      if (bj->nrings())
        continue;

      atom_number_t r01 = bj->other(c1);

      if (ring_membership[r01])
        break;

      int r01con = ncon[r01];

      if (1 == r01con)
        break;

      const Atom * ar01 = atoms[r01];

      for (int k = 0; k < r01con; k++)
      {
        const Bond * bk = ar01->item(k);

        if (bk->nrings())
          continue;

        atom_number_t r02 = bk->other(r01);

        if (ring_membership[r02])   // this also covers r02 == c1 because c1 is in a ring
          continue;

        int r02con = ncon[r02];

        if (1 == r02con)
          continue;

        const Atom * ar02 = atoms[r02];

        for (int l = 0; l < r02con; l++)
        {
          atom_number_t c2 = ar02->other(r02, l);

          if (c2 > c1)    // don't double count
            continue;

          if (aromatic[c2])
            rc++;
        }
      }

      break;   // after we find 1 bond outside the ring we are done
    }
  }

  return rc;
}

/*
  key162 is AROM
    For a count, we count the number of aromatic atoms
*/

int
MACCSKeys::_key162 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  const int * aromatic = mpr.aromatic();

  if (0 == nr)
    return 0;

  int rc = 0;
  int istop  = mpr.last_ring_atom();
  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (aromatic[i])
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  key163 is 6M Ring
    which is the same as key 145
  Change to
    a-!@[R]@[R]-!@a

  Variables are
    a1-r11-r12-a2
*/

int
MACCSKeys::_key163 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  if (m.nrings() < 3)
    return 0;

  int rc = 0;
  int istop  = mpr.last_ring_atom();

  for (int a1 = mpr.first_ring_atom(); a1 < istop; a1++)
  {
    if (! aromatic[a1])
      continue;

    if (3 != ncon[a1])
      continue;

    const Atom * aa1 = atoms[a1];

    int found_match = 0;

    for (int i = 0; i < 3; i++)
    {
      const Bond * bi = aa1->item(i);

      if (bi->nrings())
        continue;

      atom_number_t r11 = bi->other(a1);

      if (0 == ring_membership[r11])
        continue;

      int r11con = ncon[r11];

      const Atom * ar11 = atoms[r11];

      for (int j = 0; j < r11con; j++)
      {
        const Bond * bj = ar11->item(j);

        if (0 == bj->nrings())
          continue;

        atom_number_t r12 = bj->other(r11);

        if (0 == ring_membership[r12])
          continue;

        int r12con = ncon[r12];
        if (2 == r12con)
          continue;

        const Atom * ar12 = atoms[r12];
        for (int k = 0; k < r12con; k++)
        {
          const Bond * bk = ar12->item(k);

          if (bk->nrings())
            continue;

          atom_number_t a2 = bk->other(r12);

          if (a2 > a1)    // don't double count
            continue;

          if (aromatic[a2])
          {
            found_match++;
            rc++;
          }
        }
      }

      if (found_match)
        break;
    }
  }

  return rc;
}
#endif

/*
  '[/IWrh1ND3R1]1[CD2][CD2][CD3][CD2][CD2]1'
*/

static int
is_piperidine (const Ring & r,
               const MK_Molecular_Properties & mpr)
{
  assert (6 == r.number_elements());

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();

  int n = -1;
  int cd3 = -1;

  for (int i = 0; i < 6; ++i)
  {
    const auto j = r[i];
    if (7 == z[j])
    {
      if (n < 0)
        n = i;
      else
        return 0;
    }
    else if (6 != z[j])
      return 0;
    else if (3 == ncon[j])
    {
      if (cd3 < 0)
        cd3 = i;
      else
        return 0;
    }
  }

  if (n < 0 || cd3 < 0)
    return 0;

  if (n + 3 == cd3)
    return 1;
  if (cd3 + 3 == n)
    return 1;

  return 0;
}

int
MACCSKeys::_key163 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int nr = mpr.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;
  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_aromatic())
      continue;

    if (6 != ri->number_elements())
      continue;

    rc += is_piperidine(*ri, mpr);
  }

  return rc;
}

/*
  key164 is O
    which is the same as key 140
  change to
  a-!@[R0G>0]
*/

int
MACCSKeys::_key164 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;

  int istop = mpr.last_ring_atom();

  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (! aromatic[i])
      continue;

    int icon = ncon[i];

    if (3 != icon)
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (ring_membership[k])
        continue;

      if (ncon[k] < atoms[k]->nbonds())
      {
        rc++;
        break;
      }
    }
  }

  return rc;
}

/*
  key165 is RING
*/

int
MACCSKeys::_key165 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{

  int rc = mpr.nrings();

  return rc;
}

/*
  key166 is Fragments
    Don't know what this is
*/

int
MACCSKeys::_key166 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;

  int istop = mpr.last_ring_atom();

  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (0 == ring_membership[i])
      continue;

    if (aromatic[i])
      continue;

    int icon = ncon[i];

    if (3 != icon)
      continue;

    const Atom * ai = atoms[i];

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (ring_membership[k])
        continue;

      if (ncon[k] < atoms[k]->nbonds())
      {
        rc++;
        break;
      }
    }
  }

  return rc;
}

/*
  Once upon a time I had key 119 incorrect. I've corrected it and this is the old
  version that had proven useful
*/

int
MACCSKeys::_key167 (const Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (7 == z[i])           // correct definition has '7 != z[i]'
      continue;

    const Atom * a = atoms[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (b->is_double_bond())
        rc++;
    }
  }

  return rc;
}

/*
  -s '[O,N,S;D1]=[C,S,N;D3;T2]-[O,S,N;D2]-C'
*/

int
MACCSKeys::_key168 (Molecule & m,
       const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();
  const int * aromatic = mpr.aromatic();

  int rc = 0;

  int istop = mpr.last_heteroatom();
  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(i);

    if (3 != ncon[j])
      continue;

    if (2 != attached_heteroatom_count[j])
      continue;

    if (aromatic[j])
      continue;

    const auto aj = atoms[j];

    for (auto k = 0; k < 3; ++k)
    {
      const Bond * b = aj->item(k);

      if (! b->is_single_bond())
        continue;

      const auto l = b->other(j);

      if (6 == z[l])    // the second heteroatom
        continue;

      if (1 == ncon[l])    // cannot be terminal
        continue;

      if (6 == z[j])
      {
        if (0 == attached_heteroatom_count[l])
          rc++;
      }
      else
      {
        if (1 == attached_heteroatom_count[l])
          rc++;
      }
    }
  }

  return rc;
}

/*
  '[r5a]!@[r6a]'
*/

int
MACCSKeys::_key169 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  if (m.nrings() < 2)
    return 0;

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_bond_count = mpr.ring_bond_count();
  const int * aromatic = mpr.aromatic();

  int rc = 0;
  const auto istop = mpr.last_ring_atom();

  for (int i = mpr.first_ring_atom(); i < istop; i++)
  {
    if (3 != ncon[i])
      continue;

    if (2 != ring_bond_count[i])
      continue;

    if (! aromatic[i])
      continue;

    const auto ai = atoms[i];

    for (auto j = 0; j < 3; ++j)
    {
      const auto b = ai->item(j);

      if (b->nrings())
        continue;

      const auto k = b->other(i);

      if (! aromatic[k])
        continue;

      const auto ri = m.ring_containing_atom(i);
      const auto rk = m.ring_containing_atom(k);

      if (30 == ri->number_elements() * rk->number_elements())    // want a 5 and a 6
        rc++;
    }
  }

  return rc;
}


/*
  '[AR]-[R0]-[R0]-a'
*/

int
MACCSKeys::_key170 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();

  int nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (0 == ring_membership[i])
      continue;

    if (3 != ncon[i])
      continue;

    atom_number_t chain1 = INVALID_ATOM_NUMBER;

    const Atom * ai = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (0 == ring_membership[k])
      {
        chain1 = k;
        break;
      }
    }

    if (INVALID_ATOM_NUMBER == chain1)
      continue;

    const Atom * ach1 = atoms[chain1];

    for (int j = 0; j < ach1->ncon(); j++)
    {
      atom_number_t k = ach1->other(chain1, j);

      if (k == i)
        continue;

      if (ring_membership[k])
        continue;

      const Atom * ak = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        atom_number_t n2 = ak->other(k, l);

        if (n2 == chain1)
          continue;

        if (n2 < i)
          continue;

        if (0 == ring_membership[n2])
          continue;

        if (3 != ncon[n2])
          continue;

        if (aromatic[i] && ! aromatic[n2])
          rc++;
        else if (! aromatic[i] && aromatic[n2])
          rc++;
      }
    }
  }

  return rc;
}

/*
  Key 171 is 'a-[R0]-[R0]-[R0]-a'
*/

int
MACCSKeys::_key171 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const int matoms = mpr.natoms();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();

  int nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (0 == ring_membership[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (! aromatic[i])
      continue;

    atom_number_t chain1 = INVALID_ATOM_NUMBER;

    const Atom * ai = atoms[i];

    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (0 == ring_membership[k])
      {
        chain1 = k;
        break;
      }
    }

    if (INVALID_ATOM_NUMBER == chain1)
      continue;

    const Atom * ach1 = atoms[chain1];

    for (int j = 0; j < ncon[chain1]; j++)
    {
      atom_number_t k = ach1->other(chain1, j);

      if (k == i)
        continue;

      if (ring_membership[k])
        continue;

      const Atom * ak = atoms[k];

      for (int l = 0; l < ncon[k]; l++)
      {
        atom_number_t x = ak->other(k, l);

        if (x == chain1)
          continue;

        if (ring_membership[x])
          continue;

        const Atom * ax = atoms[x];

        for (int y = 0; y < ncon[x]; y++)
        {
          atom_number_t w = ax->other(x, y);

          if (k == w)
            continue;

          if (w < i)
            continue;

          if (aromatic[w])
            rc++;
        }
      }
    }
  }

  return rc;
}

/*
  Key 172 is chain Nitrogen not in spinach
*/

int
MACCSKeys::_key172 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ring_membership = mpr.ring_membership();

  if (m.nrings() < 2)
    return 0;

  const int * spinach = mpr.spinach();

  int rc = 0;
  const auto istop = mpr.last_nitrogen();

  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (spinach[i])
      continue;

    if (0 == ring_membership[i])
      rc++;
  }

  return rc;
}

/*
  Key 173 is 'N(a)(C)C'
*/

int
MACCSKeys::_key173 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();

  int nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;

  const auto istop = mpr.last_nitrogen();
  for (int i = mpr.first_nitrogen(); i < istop; i++)
  {
    if (7 != z[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (aromatic[i])
      continue;

    int aromatic = 0;
    int carbon = 0;

    const Atom * ai = atoms[i];

    if (3 != ai->nbonds())
      continue;

    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (m.is_aromatic(k))
        aromatic++;
      else if (6 == z[k])
        carbon++;
    }

    if (2 == carbon && 1 == aromatic)
      rc++;
  }

  return rc;
}

int
shortest_distance(Molecule & m, 
                  const Ring & r1,
                  const Ring & r2)
{
  int ring1_size = r1.number_elements();
  int ring2_size = r2.number_elements();

  int rc = m.natoms();

  for (int i = 0; i < ring1_size; i++)
  {
    atom_number_t a1 = r1[i];

    for (int j = 0; j < ring2_size; j++)
    {
      atom_number_t a2 = r2[j];

      if (a1 == a2)   // don't want spiro fusions
        return 0;

      int b = m.bonds_between(a1, a2);

      if (b < rc)
        rc = b;
    }
  }

  return rc;
}

/*
  Key 174 is closeness of rings.
  Note that we set five bits

  Later note: the calling programme over-writes the last bit set, don't worry about it...
*/

void
MACCSKeys::_key174 (Molecule & m,
                    const MK_Molecular_Properties & mpr,
                    int * keys_i_set) const
{

  int nr = m.nrings();

  set_vector(keys_i_set, 5, 0);

  if (nr < 2)
    return;

// Simply do not process molecules with strongly fused rings

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);
    if (ri->strongly_fused_ring_neighbours())
      return;
  }

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->fused_ring_neighbours() > 1)
      continue;

    if (ri->strongly_fused_ring_neighbours())
      continue;

    int aromi = ri->is_aromatic();

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi(j);
      if (ri->fused_system_identifier() == rj->fused_system_identifier())
        continue;

      if (rj->strongly_fused_ring_neighbours())
        continue;

      if (rj->fused_ring_neighbours() > 1)
        continue;

      int aromj = rj->is_aromatic();

      int b = shortest_distance(m, *ri, *rj);

#ifdef DEBUG_KEY_174
      cerr << "Btw " << (*ri) << " and " << (*rj) << " dist " << b << endl;
#endif

      if (b <= 4)
      {
        keys_i_set[0]++;

        if (aromi && aromj)
          keys_i_set[1]++;
        else if (0 == aromi && 0 == aromj)
          keys_i_set[2]++;
        else
          keys_i_set[3]++;
      }
      else
        keys_i_set[4]++;
    }
  }

#ifdef DEBUG_KEY_174
  for (int i = 0; i < 5; ++i)
  {
    cerr << " i = " << i << " v " << keys_i_set[i] << endl;
  }
#endif

  return;
}


/*
  Key 178 is '[OD1]=A-a'
*/

int
MACCSKeys::_key178 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * ring_membership = mpr.ring_membership();
  const int * aromatic = mpr.aromatic();

  int nr = m.nrings();

  if (nr < 1)
    return 0;

  int rc = 0;

  const auto istop = mpr.last_oxygen();
  if (INVALID_ATOM_NUMBER == istop)
    return 0;

  for (int i = mpr.first_oxygen(); i < istop; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    atom_number_t x = atoms[i]->other(i, 0);

    if (ring_membership[x])
      continue;

    if (3 != ncon[x])
      continue;

    const Atom * ax = atoms[x];

    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = ax->other(x, j);

      if (0 == ring_membership[k])
        continue;

      if (aromatic[k])
      {
        rc++;
        break;
      }
    }
  }

  return rc;
}

/*
  Key 179 is para substituted aliphatic ring
*/

int
MACCSKeys::_key179 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();

  int nr = m.nrings();

  if (nr < 1)
    return 0;

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_aromatic())
      continue;

    if (6 != ri->number_elements())
      continue;

    int first_join = -1;
    int second_join = -1;

    for (int j = 0; j < 6; j++)
    {
      atom_number_t k = ri->item(j);

      if (2 == ncon[k])
      {
        if (6 != z[k])
        {
          second_join = -1;
          break;
        }
      }
      else if (3 != ncon[k])
      {
        second_join = -1;
        break;
      }
      else if (first_join < 0)
        first_join = j;
      else if (second_join < 0)
        second_join = j;
      else
      {
        second_join = -1;
        break;
      }
    }

    if (second_join < 0)
      continue;

    if (first_join + 3 == second_join)
      rc++;
  }

  return rc;
}

/*
  Key 180 is '[!#6R0D2]-a'
*/

int
MACCSKeys::_key180 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int nr = m.nrings();

  if (nr < 1)
    return 0;

  int rc = 0;

  const auto istop = mpr.last_heteroatom();

  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (ring_membership[i])
      continue;

    if (2 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (2 != ai->nbonds())
      continue;

    int aromatic_count = 0;

    for (int j = 0; j < 2; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (aromatic[k])
        aromatic_count++;
    }

    if (1 == aromatic_count)
      rc++;
  }

  return rc;
}


/*
  6 membmered aromatic fused to an aliphatic
*/

int
MACCSKeys::_key181 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;
  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

//  cerr << "Ring " << i << " is aromatic, frn " << ri->fused_ring_neighbours() << endl;

    if (1 != ri->fused_ring_neighbours())
      continue;

    const auto r0 = ri->fused_neighbour(0);

//  cerr << " second ring strongly_fused_ring_neighbours " << r0->strongly_fused_ring_neighbours() << " total nbrs " << r0->fused_ring_neighbours() << endl;

    if (r0->strongly_fused_ring_neighbours())
      continue;

    if (r0->is_aromatic())
      continue;

    if (r0->fused_ring_neighbours() > 2)
      continue;

    rc++;
  }

  return rc;
}

/*
  '[/IWfss1aD3]1[aD2][aD3][aD2][aD2]1'
  We do both ortho and para
*/

void
MACCSKeys::_key182 (Molecule & m,
                    const MK_Molecular_Properties & mpr,
                    int * keys_i_set) const
{
  const auto nr = m.nrings();

  set_vector(keys_i_set, 2, 0);

  if (0 == nr)
    return;

  resizable_array<int> substitutions;

  const int * ncon = mpr.ncon();

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (5 != ri->number_elements())
      continue;

    if (ri->is_fused())
      continue;

    substitutions.resize_keep_storage(0);

    for (auto j = 0; j < 5; ++j)
    {
      const auto k = ri->item(j);

      if (2 == ncon[k])
        continue;

      substitutions.add(j);
    }

    if (2 != substitutions.number_elements())
      continue;

//  cerr << substitutions[0] << " and " << substitutions[1] << endl;
    if ((substitutions[0] + 1 == substitutions[1]) || (substitutions[0] + 4 == substitutions[1]))
      keys_i_set[0]++;
    else
      keys_i_set[1]++;
  }

  return;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  Atom in 3 aromatic rings
*/

int
MACCSKeys::_key184 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 3)
    return 0;

  const auto matoms = mpr.natoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (3 != ncon[i])
      continue;

    if (! aromatic[i])
      continue;

    if (3 != ring_membership[i])
      continue;

    int aromatic_rings_containing = 0;

    for (auto j = 0; j < nr; ++j)
    {
      const auto rj = m.ringi(j);

      if (! rj->is_aromatic())
        continue;

      if (rj->fused_ring_neighbours() < 2)
        continue;

      if (! rj->contains(i))
        continue;

      if (rj->strongly_fused_ring_neighbours())
      {
        aromatic_rings_containing = 0;
        break;
      }

      aromatic_rings_containing++;
    }

    if (3 == aromatic_rings_containing)
      rc++;
  }

  return rc;
}
#endif

/*
  '[/IWVy0ND3R1T0]'
*/

int
MACCSKeys::_key184 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
//const int * ring_membership = mpr.ring_membership();
  const int * attached_heteroatom_count = mpr.attached_heteroatom_count();

  const int istop = mpr.last_nitrogen();

  int rc = 0;

  for (int i = mpr.first_nitrogen(); i < istop; ++i)
  {
    if (7 != z[i])
      continue;

    if (aromatic[i])
      continue;

    if (3 != ncon[i])
      continue;

    if (attached_heteroatom_count[i])
      continue;

    if (3 != m.nbonds(i))
      continue;

    int vinyl_found = 0;
    int ring_bonds_found = 0;

    const Atom * n = atoms[i];

    for (int j = 0; j < 3; ++j)
    {
      const Bond * b = n->item(j);
      if (b->nrings())
        ring_bonds_found++;

      const auto k = b->other(i);

      if (aromatic[k])
        continue;

      if (ncon[k] < m.nbonds(k))
        vinyl_found = 1;
    }

    if (2 == ring_bonds_found && 0 == vinyl_found)
      rc++;
  }

  return rc;
}

#ifdef OLD_MK_NOT_USED_ANYMORE
/*
  '[/IWfss1a]!@[R]@[R]-[/IWfss1a]'
*/

static int
identify_fusion_point (Molecule & m,
                       const Ring & ri,
                       atom_number_t & a1,
                       const MK_Molecular_Properties & mpr)
{
  a1 = INVALID_ATOM_NUMBER;

  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();

  const auto n = ri.number_elements();

  for (auto i = 0; i < n; ++i)
  {
    const auto j = ri[i];

    if (2 == ncon[j])
      continue;

    const auto aj = atoms[j];

    for (auto k = 0; k < 3; ++k)
    {
      const auto b = aj->item(k);

      if (b->nrings())
        continue;

      const auto l = b->other(j);

      if (1 == ncon[l])
        continue;

      if (INVALID_ATOM_NUMBER == a1)
        a1 = l;
      else
        return 0;
    }
  }

  if (INVALID_ATOM_NUMBER == a1)
    return 0;

  return 1;
}

/*
  Aromatic ring !@ ring - ring !@ aromatic ring

  Aromatic rings must be mostly terminal
*/

int
MACCSKeys::_key185 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 3)
    return 0;

  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (ri->fused_ring_neighbours())
      continue;

    atom_number_t a1;     // outside the aromatic ring
    if (! identify_fusion_point(m, *ri, a1, mpr))
      continue;

    if (0 == ring_membership[a1])
      continue;

    const auto aa1 = atoms[a1];

    for (auto j = 0; j < ncon[a1]; ++j)
    {
      const auto b = aa1->item(j);

      if (0 == b->nrings())
        continue;

      const auto a2 = b->other(a1);

      if (a2 < a1)
        continue;

      const auto aa2 = atoms[a2];

      for (auto k = 0; k < ncon[a2]; ++k)     // look for a pendant aromatic
      {
        const auto b = aa2->item(k);

        if (b->nrings())
          continue;

        const auto a3 = b->other(a2);

        if (! aromatic[a3])
          continue;

        const auto r3 = m.ring_containing_atom(a3);

        if (! r3->is_aromatic())   // should not happen
          continue;

        if (r3->is_fused())
          continue;

        atom_number_t notused;
        if (! identify_fusion_point(m, *r3, notused, mpr))
          continue;

        rc++;
      }
    }
  }

  return rc;
}

#endif

/*
  Aliphatic ring with two connections
*/


static int
count_branching (Molecule & m,
                 const Ring & r,
                 const MK_Molecular_Properties & mpr)
{
  const int * ncon = mpr.ncon();
  const Atom * const * atoms = mpr.atoms();
  const int * ring_membership = mpr.ring_membership();

  const int ring_size = r.number_elements();

  int rc = 0;

  for (int i = 0; i < ring_size; ++i)
  {
    const auto j = r[i];

    if (2 == ncon[j])
      continue;

    if (1 != ring_membership[j])
      return 0;

    const auto aj = atoms[j];

    for (int k = 0; k < ncon[j]; ++k)
    {
      const Bond * b = aj->item(k);

      if (b->nrings())
        continue;

      const auto l = b->other(j);

      if (1 == ncon[l])
        continue;

      rc++;
    }
  }

  return rc;
}

/*
  two branches off an aliphatic ring
*/

int
MACCSKeys::_key185 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_aromatic())
      continue;

    if (ri->strongly_fused_ring_neighbours())
      continue;

    if (2 == count_branching(m, *ri, mpr))
      rc++;
  }

  return rc;
}

/*
  '[!#6a]:[aD3]:[!#6a]'
*/

int
MACCSKeys::_key186 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const atomic_number_t * z = mpr.z();
  const Atom * const * atoms = mpr.atoms();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int nr = m.nrings();

  if (nr < 1)
    return 0;

  int rc = 0;

  const auto istop = mpr.last_heteroatom();

  for (int i = mpr.first_heteroatom(); i < istop; i++)
  {
    if (6 == z[i])
      continue;

    if (0 == aromatic[i])
      continue;

    const auto ai = atoms[i];

    for (auto j = 0; j < ncon[i]; ++j)
    {
      const auto b = ai->item(j);

      if (0 == b->nrings())
        continue;

      const auto k = b->other(i);

      if (3 != ncon[k])
        continue;

      if (1 != ring_membership[k])
        continue;

      const auto ak = atoms[k];

      for (auto l = 0; l < ncon[k]; ++l)
      {
        const auto b = ak->item(l);

        if (0 == b->nrings())
          continue;

        const auto n = b->other(k);

        if (6 == z[n])
          continue;

        rc++;
      }
    }
  }

  return rc;
}

/*
  Piperazine, piperadine - and morpholine and  similar
*/

int
MACCSKeys::_key187 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 0;

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (ri->is_aromatic())
      continue;

    if (6 != ri->number_elements())
      continue;

    resizable_array<int> h;

    for (auto j = 0; j < 6; ++j)
    {
      const auto k = ri->item(j);

      if (7 == z[k] || 8 == z[k])
        h.add(j);
      else if (6 == z[k] && ncon[k] > 2)
        h.add(j);
    }

    if (2 != h.number_elements())
      continue;

    if (h[0] + 3 == h[1])
      rc++;
  }

  return rc;
}

/*
  '[D1]a:a!@[D>1]'
*/

int
MACCSKeys::_key188 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto matoms = mpr.natoms();

  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;
  for (auto i = 0; i < matoms; ++i)
  {
    if (1 != ncon[i])
      continue;

    const auto j = atoms[i]->other(i, 0);

    if (! aromatic[j])
      continue;

    const auto aj = atoms[j];

    for (auto k = 0; k < 3; ++k)
    {
      const auto l = aj->other(j, k);

      if (i == l)
        continue;

      if (3 != ncon[l])
        continue;

      const auto al = atoms[l];

      for (auto k = 0; k < ncon[l]; ++k)
      {
        const auto b = al->item(k);

        if (b->nrings())
          continue;

        const auto n = b->other(l);

        if (ncon[n] > 1)
          rc++;
      }
    }
  }

  return rc;
}

/*
  Strongly fused ring systems
*/

int
MACCSKeys::_key189 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{
  const auto nr = m.nrings();

  if (nr < 2)
    return 0;

  int rc = 0;

  resizable_array<int> fsid_complete;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);
    if (! ri->is_fused())
      continue;

    const auto fsid = ri->fused_system_identifier();

    if (! fsid_complete.add_if_not_already_present(fsid))
      continue;

    if (ri->strongly_fused_ring_neighbours())
      rc++;
    else
    {
      int tmp = 0;
      for (auto j = i + 1; j < nr; ++j)
      {
        const auto rj = m.ringi(j);

        if (rj->fused_system_identifier() != fsid)
          continue;

        if (! rj->strongly_fused_ring_neighbours())
          continue;

        tmp = 1;
        break;
      }
      rc += tmp;
    }
  }

  return rc;
}

/*
  '[nD3R>1]'
*/

int
MACCSKeys::_key190 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{

  const int * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const int * ring_membership = mpr.ring_membership();

  int rc = 0;

  const auto istop = mpr.last_nitrogen();
  for (auto i = mpr.first_nitrogen(); i < istop; ++i)
  {
    if (7 != z[i])
      continue;

    if (! aromatic[i])
      continue;

    if (1 == ring_membership[i])
      continue;

    if (3 != ncon[i])
      continue;

    rc++;
  }

  return rc;
}

/*
  'n:[D3a](:a):[D3a]'
*/

int
MACCSKeys::_key191 (Molecule & m,
                    const MK_Molecular_Properties & mpr) const
{

  const atomic_number_t * z = mpr.z();
  const int * ncon = mpr.ncon();
  const int * aromatic = mpr.aromatic();
  const Atom * const * atoms = mpr.atoms();

  int rc = 0;

  const auto istop = mpr.last_nitrogen();
  for (auto i = mpr.first_nitrogen(); i < istop; ++i)
  {
    if (7 != z[i])
      continue;

    if (0 == aromatic[i])
      continue;

    const auto ai = atoms[i];

    for (auto j = 0; j < ncon[i]; ++j)
    {
      const auto b = ai->item(j);

      if (! b->is_aromatic())
        continue;

      const auto k = b->other(i);

      if (3 != ncon[k])
        continue;

      const auto ak = atoms[k];

      int ad3 = 0;
      for (auto l = 0; l < ncon[k]; ++l)
      {
        const auto b = ak->item(l);

        if (! b->is_aromatic())
          continue;

        const auto n = b->other(k);

        if (n == i)
          continue;

        if (3 == ncon[n])
          ad3++;
      }

      if (ad3)
        rc++;
    }
  }

  return rc;
}

int
MACCSKeys::operator() (Molecule & m,  int * keys) const
{
  MK_Molecular_Properties mpr(m);

#ifdef MK_TIME_TEST
  time_t tzero = time(NULL);
  int niter = 101411612;
  for (auto i = 0; i < niter; ++i)
  {
    (void) _key4(m, mpr);
  }
  time_t t = time(NULL);
  cerr << "Key 4 " << (t - tzero) << endl;
  tzero = t;
  for (auto i = 0; i < niter; ++i)
  {
    (void) _key18(m, mpr);
  }
  t = time(NULL);
  cerr << "Key 18 " << (t - tzero) << endl;
#endif

  if (m.number_fragments() > 1) {
    cerr << "MACCSKeys:Cannot process multi fragment molecules '" << m.name() << "'\n";
    return 0;
  }

  keys[0] = _key0(m, mpr);
  keys[1] = _key1(m, mpr);
  keys[2] = _key2(m, mpr);
  keys[3] = _key3(m, mpr);
  keys[4] = _key4(m, mpr);
  keys[5] = _key5(m, mpr);
  keys[6] = _key6(m, mpr);
  keys[7] = _key7(m, mpr);
  keys[8] = _key8(m, mpr);
  keys[9] = _key9(m, mpr);
  keys[10] = _key10(m, mpr);
  keys[11] = _key11(m, mpr);
  keys[12] = _key12(m, mpr);
  keys[13] = _key13(m, mpr);
  keys[14] = _key14(m, mpr);
  keys[15] = _key15(m, mpr);
  keys[16] = _key16(m, mpr);
  keys[17] = _key17(m, mpr);
  keys[18] = _key18(m, mpr);
  keys[19] = _key19(m, mpr);
  keys[20] = _key20(m, mpr);
  keys[21] = _key21(m, mpr);
  keys[22] = _key22(m, mpr);
  keys[23] = _key23(m, mpr);
  keys[24] = _key24(m, mpr);
  keys[25] = _key25(m, mpr);
  keys[26] = _key26(m, mpr);
  keys[27] = _key27(m, mpr);
  keys[28] = _key28(m, mpr);
  keys[29] = _key29(m, mpr);
  keys[30] = _key30(m, mpr);
  keys[31] = _key31(m, mpr);
  keys[32] = _key32(m, mpr);
  keys[33] = _key33(m, mpr);
  keys[34] = _key34(m, mpr);
  keys[35] = _key35(m, mpr);
  keys[36] = _key36(m, mpr);
  keys[37] = _key37(m, mpr);
  keys[38] = _key38(m, mpr);
  keys[39] = _key39(m, mpr);
  keys[40] = _key40(m, mpr);
  keys[41] = _key41(m, mpr);
  keys[42] = _key42(m, mpr);
  keys[43] = _key43(m, mpr);
  keys[44] = _key44(m, mpr);
  keys[45] = _key45(m, mpr);
  keys[46] = _key46(m, mpr);
  keys[47] = _key47(m, mpr);
  keys[48] = _key48(m, mpr);
  keys[49] = _key49(m, mpr);
  keys[50] = _key50(m, mpr);
  keys[51] = _key51(m, mpr);
  keys[52] = _key52(m, mpr);
  keys[53] = _key53(m, mpr);
  keys[54] = _key54(m, mpr);
  keys[55] = _key55(m, mpr);
  keys[56] = _key56(m, mpr);
  keys[57] = _key57(m, mpr);
  keys[58] = _key58(m, mpr);
  keys[59] = _key59(m, mpr);
  keys[60] = _key60(m, mpr);
  keys[61] = _key61(m, mpr);
  keys[62] = _key62(m, mpr);
  keys[63] = _key63(m, mpr);
  keys[64] = _key64(m, mpr);
  keys[65] = _key65(m, mpr);
  keys[66] = _key66(m, mpr);
  keys[67] = _key67(m, mpr);
  keys[68] = _key68(m, mpr);
  keys[69] = _key69(m, mpr);
  keys[70] = _key70(m, mpr);
  keys[71] = _key71(m, mpr);
  keys[72] = _key72(m, mpr);
  keys[73] = _key73(m, mpr);
  keys[74] = _key74(m, mpr);
  keys[75] = _key75(m, mpr);
  keys[76] = _key76(m, mpr);
  keys[77] = _key77(m, mpr);
  keys[78] = _key78(m, mpr);
  keys[79] = _key79(m, mpr);
  keys[80] = _key80(m, mpr);
  keys[81] = _key81(m, mpr);
  keys[82] = _key82(m, mpr);
  keys[83] = _key83(m, mpr);
  keys[84] = _key84(m, mpr);
  keys[85] = _key85(m, mpr);
  keys[86] = _key86(m, mpr);
  keys[87] = _key87(m, mpr);
  keys[88] = _key88(m, mpr);
  keys[89] = _key89(m, mpr);
  keys[90] = _key90(m, mpr);
  keys[91] = _key91(m, mpr);
  keys[92] = _key92(m, mpr);
  keys[93] = _key93(m, mpr);
  keys[94] = _key94(m, mpr);
  keys[95] = _key95(m, mpr);
  keys[96] = _key96(m, mpr);
  keys[97] = _key97(m, mpr);
  keys[98] = _key98(m, mpr);
  keys[99] = _key99(m, mpr);
  keys[100] = _key100(m, mpr);
  keys[101] = _key101(m, mpr);
  keys[102] = _key102(m, mpr);
  keys[103] = _key103(m, mpr);
  keys[104] = _key104(m, mpr);
  keys[105] = _key105(m, mpr);
  keys[106] = _key106(m, mpr);
  keys[107] = _key107(m, mpr);
  keys[108] = _key108(m, mpr);
  keys[109] = _key109(m, mpr);
  keys[110] = _key110(m, mpr);
  keys[111] = _key111(m, mpr);
  keys[112] = _key112(m, mpr);
  keys[113] = _key113(m, mpr);
  keys[114] = _key114(m, mpr);
  keys[115] = _key115(m, mpr);
  keys[116] = _key116(m, mpr);
  keys[117] = _key117(m, mpr);
  keys[118] = _key118(m, mpr);
  keys[119] = _key119(m, mpr);
  keys[120] = _key120(m, mpr);
  keys[121] = _key121(m, mpr);
  keys[122] = _key122(m, mpr);
  keys[123] = _key123(m, mpr);
  keys[124] = _key124(m, mpr);
  keys[125] = _key125(m, mpr);
  keys[126] = _key126(m, mpr);
  keys[127] = _key127(m, mpr);
  keys[128] = _key128(m, mpr);
  keys[129] = _key129(m, mpr);
  keys[130] = _key130(m, mpr);
  keys[131] = _key131(m, mpr);
  keys[132] = _key132(m, mpr);
  keys[133] = _key133(m, mpr);
  keys[134] = _key134(m, mpr);
  keys[135] = _key135(m, mpr);
  keys[136] = _key136(m, mpr);
  keys[137] = _key137(m, mpr);
  keys[138] = _key138(m, mpr);
  keys[139] = _key139(m, mpr);
  keys[140] = _key140(m, mpr);
  keys[141] = _key141(m, mpr);
  keys[142] = _key142(m, mpr);
  keys[143] = _key143(m, mpr);
  keys[144] = _key144(m, mpr);
  keys[145] = _key145(m, mpr);
  keys[146] = _key146(m, mpr);
  keys[147] = _key147(m, mpr);
  keys[148] = _key148(m, mpr);
  keys[149] = _key149(m, mpr);
  keys[150] = _key150(m, mpr);
  keys[151] = _key151(m, mpr);
  keys[152] = _key152(m, mpr);
  keys[153] = _key153(m, mpr);
  keys[154] = _key154(m, mpr);
  keys[155] = _key155(m, mpr);
  keys[156] = _key156(m, mpr);
  keys[157] = _key157(m, mpr);
  keys[158] = _key158(m, mpr);
  keys[159] = _key159(m, mpr);
  keys[160] = _key160(m, mpr);
  keys[161] = _key161(m, mpr);
  keys[162] = _key162(m, mpr);
  keys[163] = _key163(m, mpr);
  keys[164] = _key164(m, mpr);
  keys[165] = _key165(m, mpr);
  keys[166] = _key166(m, mpr);

  if (166 ==_nbits)
    return 1;

  keys[167] = _key167(m, mpr);
  keys[168] = _key168(m, mpr);
  keys[169] = _key169(m, mpr);
  keys[170] = _key170(m, mpr);
  keys[171] = _key171(m, mpr);
  keys[172] = _key172(m, mpr);
  keys[173] = _key173(m, mpr);

  _key174(m, mpr, keys + 174);

  keys[178] = _key178(m, mpr);
  keys[179] = _key179(m, mpr);
  keys[180] = _key180(m, mpr);
  keys[181] = _key181(m, mpr);

  _key182(m, mpr, keys + 182);

  keys[184] = _key184(m, mpr);
  keys[185] = _key185(m, mpr);
  keys[186] = _key186(m, mpr);
  keys[187] = _key187(m, mpr);
  keys[188] = _key188(m, mpr);
  keys[189] = _key189(m, mpr);
  keys[190] = _key190(m, mpr);
  keys[191] = _key191(m, mpr);
  
  return 1;
}

int
MACCSKeys::set_level_2_fingerprint (int * v) const
{
  int rc = 0;
  for (auto i = 0; i < _nbits; ++i)
  {
    if (v[i] >= level2_threshold[i])
    {
      v[i] = 1;
      rc++;
    }
    else
      v[i] = 0;
  }

  return rc;
}

