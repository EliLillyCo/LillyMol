#include <iostream>
#include <limits>

#include "Foundational/iwmisc/misc.h"

#include "atom_pair_fingerprint.h"

namespace atom_pair_fingerprint {

using std::cerr;
using std::endl;

using std::tie;

AtomPairFingerprint::AtomPairFingerprint() {
  _min_separation = 1;
  _max_separation = std::numeric_limits<int>::max();

  _fingerprint_bonded_atoms_with_btype = false;

  _check_for_collisions = false;

  _collisions_found = 0;

  _need_to_examine_bits_formed = false;

  _include_out_of_range_separations = false;

  _examine_coverage = false;

  _coverage = nullptr;

  // By setting the magic numbers to different values, different fingerprints can
  // be formed.

  _bond_magic1 = 52;
  _bond_magic2 = 4;
  _bond_magic3 = 202;
  _bond_magic4 = 312;

  _magic1 = 5092308;
  _magic2 = 389932;
  _magic3 = 755551;
  _magic4 = 613868;
}

AtomPairFingerprint::~AtomPairFingerprint()
{
  if (nullptr != _coverage)
    delete [] _coverage;
}

int
AtomPairFingerprint::OpenStreamForBitMeanings(const char * fname) {
  if (_stream_for_bit_meanings.is_open()) {
    cerr << "AtomPairFingerprint::OpenStreamForBitMeanings:already open\n";
    return 0;
  }

  if (! _stream_for_bit_meanings.open(fname)) {
    cerr << "AtomPairFingerprint::OpenStreamForBitMeanings:cannot open '" << fname << "'\n";
    return 0;
  }

  _need_to_examine_bits_formed = true;

  return 1;
}

int
AtomPairFingerprint::Fingerprint(Molecule& m,
        const int * include_atom,
        const atom_type_t* atom_type,
        Sparse_Fingerprint_Creator& sfc) {

  if (0 == m.natoms()) {
    return 0;
  }

  m.recompute_distance_matrix();

  if (nullptr == include_atom && 1 == m.number_fragments() && ! _examine_coverage) {
    return _FingerprintFast(m, atom_type, sfc);
  }

  // Construct a list of atoms to be processed

  const int matoms = m.natoms();

  atom_number_t * tmp = new atom_number_t[matoms]; std::unique_ptr<atom_number_t[]> free_tmp(tmp);

  int nprocess = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (nullptr == include_atom)
    {
      tmp[nprocess] = i;
      nprocess++;
    }
    else if (include_atom[i])
    {
      tmp[nprocess] = i;
      nprocess++;
    }
  }

  if (_examine_coverage)
    _coverage = new_int(matoms);
  else
    _coverage = nullptr;

  return _Fingerprint(m, tmp, nprocess, atom_type, sfc);
}

int
AtomPairFingerprint::_Fingerprint(Molecule& m,
        const atom_number_t * atoms_to_process,
        const int nprocess,
        const atom_type_t* atom_type,
        Sparse_Fingerprint_Creator& sfc) {

  if (0 == _min_separation) {
    for (int i = 0; i < nprocess; ++i) {
      sfc.hit_bit(atom_type[atoms_to_process[i]] + _magic1);
    }
  }

  APairOfAtoms a_pair(m, atom_type, sfc);

  const int nfrag = m.number_fragments();

  for (int i = 0; i < nprocess; ++i) {
    const atom_number_t a1 = atoms_to_process[i];
    a_pair.SetA1(a1);

    // Only set fragment membership if needed.
    const int frag1 = nfrag > 1 ? m.fragment_membership(a1) : 0;

    for (int j = i + 1; j < nprocess; ++j) {
      const atom_number_t a2 = atoms_to_process[j];

      if (1 == nfrag)  // no need to check
        ;
      else if (frag1 != m.fragment_membership(a2))
        continue;

      a_pair.SetA2(a2);
      _FormFingerprint(a_pair);
    }
  }

  if (_examine_coverage)
    _CoverageToMolecule(a_pair);

  return sfc.nbits();
}

int
AtomPairFingerprint::_FingerprintFast(Molecule& m,
        const atom_type_t* atom_type,
        Sparse_Fingerprint_Creator& sfc) {
  const int matoms = m.natoms();

  if (0 == _min_separation) {
    for (int i = 0; i < matoms; ++i) {
      sfc.hit_bit(atom_type[i] + _magic1);
    }
  }

  APairOfAtoms a_pair(m, atom_type, sfc);

  for (int i = 0; i < matoms; ++i) {
    a_pair.SetA1(i);
    for (int j = i + 1; j < matoms; ++j) {
      a_pair.SetA2(j);
      _FormFingerprint(a_pair);
    }
  }

  return sfc.nbits();
}

void
AtomPairFingerprint::_FormFingerprint(APairOfAtoms& a_pair) {

  int d = a_pair.d();
  if (0 == d)
    abort();

  if (d < _min_separation)
    return;

  if (d <= _max_separation)  // Great.
    ;
  else if (_include_out_of_range_separations)  // Adjust the distance.
  {
    d = _max_separation + 1;
    a_pair._d = d;
  }
  else  // Not processed.
    return;

  if (1 == d && _fingerprint_bonded_atoms_with_btype)
    _FingerprintBond(a_pair);
  else
    _FingerprintPair(a_pair);
}

uint64_t
AtomPairFingerprint::_BondConstant(const Bond & b) const
{
  if (b.is_aromatic())
    return _bond_magic4;
  
  if (b.is_single_bond())
    return _bond_magic1;
  if (b.is_double_bond())
    return _bond_magic2;
  if (b.is_triple_bond())
    return _bond_magic3;

  cerr << "AtomPairFingerprint::_BondConstant:what kind of bond ";
  b.debug_print(cerr);

  return 1;
}

void
AtomPairFingerprint::_FingerprintBond(APairOfAtoms& a_pair) {
  const Bond * b = a_pair.BondBetweenAtoms();
  assert(nullptr != b);

  const uint64_t btype = _BondConstant(*b);

  uint64_t t1, t2;
  tie (t1, t2) = a_pair.OrderedAtomTypes();

  const uint64_t bitnum = t1 * _magic1 + btype * (t2 + _magic2);
  a_pair.HitBit(bitnum);

  if (_need_to_examine_bits_formed)
    _ExamineBit(a_pair);
}


void
AtomPairFingerprint::_FingerprintPair(APairOfAtoms& a_pair) {
  uint64_t t1, t2;
  tie(t1, t2) = a_pair.OrderedAtomTypes();

  const uint64_t b = t1 * _magic3 + a_pair.d() * (_magic4 + t2);
  a_pair.HitBit(b);

  if (_need_to_examine_bits_formed)
    _ExamineBit(a_pair);

  if (_examine_coverage)
    _UpdateCoverage(a_pair);
}

void
AtomPairFingerprint::_ExamineBit(APairOfAtoms& a_pair) {
  if (_check_for_collisions)
    _DoCheckForCollisions(a_pair);

  if (_stream_for_bit_meanings.is_open())
    _DoWriteBitMeanings(a_pair);

  if (_examine_coverage)
    _UpdateCoverage(a_pair);

  return;
}

void
AtomPairFingerprint::_DoCheckForCollisions(const APairOfAtoms& a_pair) {

  uint64_t t1, t2;
  tie(t1, t2) = a_pair.OrderedAtomTypes();

  std::tuple<uint64_t, int, uint64_t> tmp{t1, a_pair.d(), t2};

  const uint64_t bitnum = a_pair.bitnum();

  const auto f = _bit_explanation.find(bitnum);
  if (f == _bit_explanation.end()) {  // First time.
    _bit_explanation.emplace(bitnum, std::move(tmp));
    return;
  }

  if (f->second == tmp)
    return;

//#define ECHO_BIT_COLLISIONS
#ifdef ECHO_BIT_COLLISIONS
  cerr << "Collision on bit " << bitnum << endl;
  cerr << std::get<0>(tmp) << ' ' << std::get<0>(f->second) << endl;
  cerr << std::get<1>(tmp) << ' ' << std::get<1>(f->second) << endl;
  cerr << std::get<2>(tmp) << ' ' << std::get<2>(f->second) << endl;
#endif

  _collisions_found++;

  return;
}

void
AtomPairFingerprint::_DoWriteBitMeanings(APairOfAtoms& a_pair) {
  uint64_t t1, t2;
  tie(t1, t2) = a_pair.OrderedAtomTypes();

  const atom_number_t a1 = a_pair.a1();
  const atom_number_t a2 = a_pair.a2();

  _stream_for_bit_meanings << a_pair.smiles() << ' ' << a_pair.name() << 
      ' ' << a1 << ' ' << a_pair.m().smarts_equivalent_for_atom(a1) << ' ' << a_pair.atype(a1) <<
      ' ' << a_pair.d() <<
      ' ' << a2 << ' ' << a_pair.m().smarts_equivalent_for_atom(a2) << ' ' << a_pair.atype(a2) <<
      ' ' << a_pair.bitnum() << '\n';

  _stream_for_bit_meanings.write_if_buffer_holds_more_than(8192);

  return;
}

void
AtomPairFingerprint::_UpdateCoverage(const APairOfAtoms& a_pair) {
  _coverage[a_pair.a1()]++;
  _coverage[a_pair.a2()]++;

  return;
}

void
AtomPairFingerprint::_CoverageToMolecule(APairOfAtoms& a_pair)
{
  assert (nullptr != _coverage);

  const Molecule& molref = a_pair.m();

  Molecule * m = const_cast<Molecule *>(&molref);

  const int matoms = m->natoms();
  for (int i = 0; i < matoms; ++i) {
    m->set_atom_map_number(i, _coverage[i]);
  }

  return;
}

int
AtomPairFingerprint::ReportCollisions(std::ostream& output) const {
  output << "Accumulated " << _bit_explanation.size() << " bits, found " << _collisions_found << " collisions\n";

  return output.good();
}

}  // namespace atom_pair_fingerprint
