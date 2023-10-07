#include <algorithm>
#include <iostream>

#include "linear_fingerprint.h"

#include "Foundational/iwmisc/misc.h"

namespace linear_fingerprint {

using std::cerr;
using std::endl;

namespace internal {

constexpr int exclude_atom = -1;

Options::Options () 
{
  _min_length = 0;
  _max_length = 7;

  _paths_can_cross = false;
  _fingerprint_ring_presence = false;
  _check_for_collisions = false;
  _check_coverage = false;

  return;
}

LinearFpStatus::LinearFpStatus(const Options& opt, const Molecule& m,
    const int * include_atom,
    const atom_type_t * atype,
    Sparse_Fingerprint_Creator& sfc) :
                              _options(opt),
                              _matoms(m.natoms()),
                              _nedges(m.nedges()),
                              _m(m),
                              _atype(atype),
                              _bond_list(m.bond_list()),
                              _sfc(sfc)
{
  // These should be const, but the constructor already has too many things being
  // initialised.
  _magic1 = 412333;
  _magic2 = 1929111;
  _magic3 = 698212;

  _bond_magic1 = 102;
  _bond_magic2 = 351;
  _bond_magic3 = 6001;
  _bond_magic4 = 52;

  if (_options._check_coverage)
    _coverage = new_int(_matoms);
  else
    _coverage = nullptr;

  _stream_for_bit_meanings = nullptr;

  _need_to_examine_bits_formed = false;

  _atom_in_path = new_int(_matoms);

  _path_length = 0;

  _path = new uint64_t[2 * _nedges + 1];
  _path_index = new int[2 * _nedges + 2];

  if (0 == _nedges)  // Single atom, no bond info.
  {
    _bond_constant = nullptr;
    _bond_in_path = nullptr;
    _atom = nullptr;

    return;
  }

  _atom = new const Atom*[_matoms];

  _m.atoms(_atom);

  _bond_constant = new uint64_t[_nedges];
  
  _bond_in_path = new_int(_nedges);

  if (nullptr != include_atom) {
    for (int i = 0; i < _matoms; ++i)
    {
      if (include_atom[i])
        _atom_in_path[i] = 0;
      else
        _atom_in_path[i] = exclude_atom;
    }

    for (int i = 0; i < _nedges; ++i) {
      const Bond * b = m.bondi(i);
      const int bond_number = b->bond_number();
      assert(bond_number < _nedges);

      if (include_atom[b->a1()] && include_atom[b->a2()]) {
        _bond_constant[bond_number] = _BondHash(*b);
      }
      else
        _bond_in_path[bond_number] = 1;
    }
  }
  else 
  {
    for (int i = 0; i < _nedges; ++i) {
      const Bond * b = _bond_list[i];
      _bond_constant[b->bond_number()] = _BondHash(*b);
    }
  }

  return;
}

LinearFpStatus::~LinearFpStatus()
{
  if (nullptr != _atom)
    delete [] _atom;
  if (nullptr != _bond_constant)
    delete [] _bond_constant;
  delete [] _atom_in_path;
  if (nullptr != _bond_in_path)
    delete [] _bond_in_path;
  delete [] _path;
  delete [] _path_index;

  if (nullptr != _coverage)
    delete [] _coverage;

  return;
}

uint64_t
LinearFpStatus::_BondHash(const Bond& b) const {
  if (b.is_aromatic())
    return _bond_magic4;

  if (b.is_single_bond())
    return _bond_magic1;

  if (b.is_double_bond())
    return _bond_magic2;

  if (b.is_triple_bond())
    return _bond_magic3;

  cerr << "LinearFpStatus:_BondHash:what kind of bond ";
  b.debug_print(cerr);

  return 1;
}

void
LinearFpStatus::_AddBondToPath(const Bond & b, const atom_number_t next_atom) {
  const int bond_number = b.bond_number();
  assert(0 == _bond_in_path[bond_number]);

#ifdef DEBUG_LINEAR_FP
  cerr << "At length " << _path_length << " adding bond number " << bond_number << " value " << _bond_constant[bond_number] << endl;
#endif

  _path[_path_length] = _bond_constant[bond_number];
  _path_index[_path_length] = bond_number;

  _bond_in_path[bond_number] = 1;

  _path_length++;

  _path[_path_length] = _atype[next_atom];
  _path_index[_path_length] = next_atom;

  _atom_in_path[next_atom]++;

  _path_length++;

}

void
LinearFpStatus::_PopPath() {

  _path_length--;

  const int atom_number = _path_index[_path_length];
  assert(_atom_in_path[atom_number]);
  _atom_in_path[atom_number]--;

  _path_length--;

  const int bond_number = _path_index[_path_length];
  assert(_bond_in_path[bond_number]);
  _bond_in_path[bond_number] = 0;

  return;
}

void
LinearFpStatus::_HandleMoleculeWithNoBonds()
{
  if (_options._min_length > 0)
    return;

  for (int i = 0; i < _matoms; ++i)
  {
    _StartPath(i);
    _MaybeFormBit();
  }

  return;
}

int
LinearFpStatus::Fingerprint() {
  if (0 == _matoms)
    return 0;

  _path_length = 1;  // This function places a single atom in _path

  if (0 == _nedges) {
    _HandleMoleculeWithNoBonds();
    return 1;
  }

  if (nullptr != _stream_for_bit_meanings)
    _WriteLabelledSmiles();

  for (int i = 0; i < _matoms; ++i)
  {
    if (exclude_atom == _atom_in_path[i])
      continue;

    _StartPath(i);
    _MaybeFormBit();
    if (_options._max_length > 0) {
      _Expand();
    }
    _path_length = 1;
    _atom_in_path[i] = 0;
  }

  if (_options._check_coverage)
    _AssignCoverageDataToMolecule();

  return _sfc.nbits();
}

void
LinearFpStatus::_Expand() 
{
  if (_path_length / 2 >= _options._max_length)
    return;

  const atom_number_t a1 = _path_index[_path_length - 1];

  const Atom * a = _atom[a1];

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i) {
    const Bond * b = a->item(i);
    if (_bond_in_path[b->bond_number()])
      continue;

    const atom_number_t a2 = b->other(a1);
    if (exclude_atom == _atom_in_path[a2])
      continue;

    bool a2_already_in_path;
    if (!_atom_in_path[a2])  // The easy case
      a2_already_in_path = false;
    else if (_options._fingerprint_ring_presence ||
             _options._paths_can_cross)
      a2_already_in_path = true;
    else  // Avoid placed atom.
      continue;

    _AddBondToPath(*b, a2);
    _MaybeFormBit();
    if (a2_already_in_path) 
    {
      if (_options._fingerprint_ring_presence)
        _FormRingBit();
      if (!_options._paths_can_cross)
        continue;
    }
    _Expand();
    _PopPath();
  }
}

void
LinearFpStatus::_StartPath(const int atom_number) 
{
  assert(1 == _path_length);  // Shortcut

  _path[0] = _atype[atom_number];
  _atom_in_path[atom_number] = 1;
  _path_index[0] = atom_number;

  return;
}

void
LinearFpStatus::_MaybeFormBit() {
  if (_path_length / 2 < _options._min_length)
    return;

  if (_path_length / 2 > _options._max_length)  // Should never happen.
    return;

#ifdef DEBUG_LINEAR_FP
  cerr << "LinearFpStatus::_MaybeFormBit\n";
  _PrintPath(cerr);
#endif

  if (1 == _path_length) {  // Nothing to canonicalise.
    _FormFingerprintForward();
    return;
  }

  if (_path_index[0] < _path_index[_path_length - 1])  // Do not double count paths
    return;

  // Need a canonical direction for the path

  const int half = _path_length / 2;

  for (int lhs = 0; lhs < half; ++lhs) {
    const int rhs = _path_length - lhs - 1;
    if (_path[lhs] == _path[rhs])  // Not resolved.
      continue;
    
    if (_path[lhs] < _path[rhs])
      _FormFingerprintForward();
    else
      _FormFingerprintBackward();

    return;
  }

  // Tie not broken, any direction will do.
  _FormFingerprintBackward();
}

// Atom at end of path occurs somewhere previously. Find it.
void
LinearFpStatus::_FormRingBit()
{
  const int target = _path_index[_path_length - 1];

  int first_index = -1;
  for (int i = 0; i < (_path_length - 1); i += 2)
  {
    if (_path_index[i] == target)
    {
      first_index = i;
      break;
    }
  }

  if (first_index < 0) 
  {
    cerr << "LinearFpStatus:_FormRingBit:first occurrence not found\n";
    _PrintPath(cerr);
    return;
  }

  uint64_t t1 = _path[first_index];
  uint64_t t2 = _path[_path_length - 1];

  if (t1 < t2)
    std::swap(t1, t2);

  _sfc.hit_bit(_magic1 * t1 + (_path_length - first_index) * (t2 + _magic2));
}

void
LinearFpStatus::_FormFingerprintForward() 
{
  uint64_t b = 0;
  for (int i = 0; i < _path_length; ++i)
  {
    b += (i + _magic3) * _path[i];
  }

  _sfc.hit_bit(b);

  _ExamineBit(b);

  return;
}

void
LinearFpStatus::_FormFingerprintBackward() 
{
  uint64_t b = 0;
  int factor = _magic3;
  for (int i = _path_length - 1; i >= 0; --i, ++factor)
  {
    b += factor * _path[i];
  }

  _sfc.hit_bit(b);

  _ExamineBit(b);

  return;
}

void
LinearFpStatus::_ExamineBit(const uint64_t b)
{
#ifdef DEBUG_LINEAR_FP
  cerr << "Formed bit " << b <<endl;
#endif

  if (! _need_to_examine_bits_formed)
    return;

  if (_options._check_coverage)
    _DoCheckCoverage();

  if (nullptr != _stream_for_bit_meanings)
    _WriteBit(b);

  return;
}

void
LinearFpStatus::_WriteLabelledSmiles() const {
  Molecule tmp(_m);
  for (int i = 0; i < _matoms; ++i) {
    tmp.set_atom_map_number(i, i);
  }

  *_stream_for_bit_meanings << tmp.smiles() << ' ' << _m.name() << "\n";
  _stream_for_bit_meanings->write_if_buffer_holds_more_than(8192);
}

void
LinearFpStatus::_AssignCoverageDataToMolecule() const {
  Molecule* tmp = const_cast<Molecule*>(&_m);   // loss of const OK.
  for (int i = 0; i < _matoms; ++i) {
    tmp->set_atom_map_number(i, _coverage[i]);
  }
}

void
LinearFpStatus::_DoCheckCoverage()
{
  for (int i = 0; i < _path_length; i += 2) {
    _coverage[_path_index[i]]++;
  }

  return;
}

void
LinearFpStatus::_PrintPath(std::ostream& output) const
{
  output << "Path with " << (_path_length/2) << " bonds\n";
  for (int i = 0; i < _path_length; ++i) 
  {
    if (0 == i % 2)   // an atom
    {
      output << " atom " << _path_index[i];
    }
    else  // a bond
    {
      output << " bond " << _path_index[i];
    }
  }
  output << "\n";

  for (int i = 0; i < _path_length; ++i) 
  {
    if (0 == i % 2)   // an atom
    {
      output << " atom " << _path[i];
    }
    else  // a bond
    {
      output << " bond " << _path[i];
    }
  }

  output << "\n";

  return;
}

void
LinearFpStatus::_WriteBit(const uint64_t b) {
  Molecule tmp(_m);  // Very inefficient.

  for (int i = 0; i < _path_length; i += 2) {
    tmp.set_atom_map_number(_path_index[i], i + 1);
  }

  *_stream_for_bit_meanings << tmp.smiles() << ' ' << _m.name() << ' ' << _path_length << ' ' << b << "\n";
  _stream_for_bit_meanings->write_if_buffer_holds_more_than(8192);
}

}  // namespace internal

LinearFingerprintGenerator::LinearFingerprintGenerator()
{
}

int
LinearFingerprintGenerator::Fingerprint(Molecule & m,
    const int * include_atom,
    const atom_type_t * atype,
    Sparse_Fingerprint_Creator& sfp)
{
  if (0 == m.natoms())
    return 0;

  m.assign_bond_numbers_to_bonds_if_needed();
  m.compute_aromaticity_if_needed();

  internal::LinearFpStatus lfp(_options, m, include_atom, atype, sfp);

  if (_stream_for_bit_meanings.is_open())
    lfp.SetStreamForBitMeanings(&_stream_for_bit_meanings);

  return lfp.Fingerprint();
}

int
LinearFingerprintGenerator::OpenStreamForBitMeanings(const char * fname) {
  if (_stream_for_bit_meanings.is_open()) {
    cerr << "LinearFingerprintGenerator::OpenStreamForBitMeanings:already open\n";
    return 0;
  }

  if (! _stream_for_bit_meanings.open(fname)) {
    cerr << "LinearFingerprintGenerator::OpenStreamForBitMeanings:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
LinearFingerprintGenerator::ReportCollisions(std::ostream& output) const 
{
  output << "LinearFingerprintGenerator::ReportStatus:not implemented\n";
  return output.good();
}

}  // namespace linear_fingerprint
