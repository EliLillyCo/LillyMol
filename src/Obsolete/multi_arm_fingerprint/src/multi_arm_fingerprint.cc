#include "multi_arm_fingerprint.h"

#include <algorithm>
#include <array>
#include <limits>
#include <vector>

namespace multi_arm_fingerprint {
namespace {

// MurmurHash3's 32-bit word mixing and avalanche, expressed on unsigned words.
// Fields are mixed separately, never via machine-dependent byte representations.
uint32_t
Mix(uint32_t h, uint32_t word) {
  word *= 0xcc9e2d51u;
  word = (word << 15) | (word >> 17);
  word *= 0x1b873593u;
  h ^= word;
  h = (h << 13) | (h >> 19);
  return h * 5u + 0xe6546b64u;
}

uint32_t
Finish(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6bu;
  h ^= h >> 13;
  h *= 0xc2b2ae35u;
  return h ^ (h >> 16);
}

uint32_t
BondCode(const Bond& bond) {
  if (bond.is_aromatic()) {
    return 1;
  }
  if (bond.is_single_bond()) {
    return 2;
  }
  if (bond.is_double_bond()) {
    return 3;
  }
  if (bond.is_triple_bond()) {
    return 4;
  }
  return 5;  // Other bond types deliberately share one category.
}

struct ArmPath {
  uint32_t _hash;
  bool _terminal;
};

using Arm = std::vector<std::vector<ArmPath>>;

class Enumeration {
 private:
  const Molecule& _m;
  const uint32_t* _types;
  const Options& _options;
  Sparse_Fingerprint_Creator& _result;
  Statistics& _statistics;
  // Path position, or -1 if absent. Shared scratch, restored after backtracking.
  std::vector<int> _position;
  std::vector<Arm> _arms;
  std::array<uint32_t, 4> _selected{};
  uint32_t _centre_type;
  int _limit;

  void
  Extend(int atom, int parent, int length, uint32_t hash, Arm& arm) {
    arm[length].push_back({Finish(Mix(hash, length)), false});
    ++_statistics._arm_paths;
    if (length == _limit) {
      return;
    }
    bool extended = false;
    for (const Bond* bond : _m[atom]) {
      const int next = bond->other(atom);
      if (next == parent) {
        continue;
      }
      uint32_t next_hash = Mix(hash, BondCode(*bond));
      if (_position[next] >= 0) {
        if (!_options._ring_closures) {
          continue;
        }
        next_hash = Mix(next_hash, 0x52494e47u);  // RING field marker.
        next_hash = Mix(next_hash, length - _position[next] + 1);
        arm[length + 1].push_back({Finish(Mix(next_hash, length + 1)), true});
        ++_statistics._arm_paths;
        extended = true;
        continue;
      }
      extended = true;
      next_hash = Mix(Mix(next_hash, 0x41544f4du), _types[next]);
      _position[next] = length + 1;
      Extend(next, atom, length + 1, next_hash, arm);
      _position[next] = -1;
    }
    // No permitted continuation: this path can persist at later radii.
    arm[length].back()._terminal = !extended;
  }

  // Canonicalise a complete selection and record its fingerprint occurrence.
  void
  RecordCombination() {
    const int narms = _arms.size();
    auto sorted = _selected;
    // Sort a fixed-size array. For three arms the sentinel sorts last and is
    // excluded from hashing; repeated arm hashes remain part of the feature.
    if (narms == 3) {
      sorted[3] = std::numeric_limits<uint32_t>::max();
    }
    std::sort(sorted.begin(), sorted.end());
    uint32_t hash = Mix(Mix(0x4d415231u, _centre_type), narms);  // MAR1 version seed.
    for (int i = 0; i < narms; ++i) {
      hash = Mix(hash, sorted[i]);
    }
    const uint32_t bit = Finish(hash);
    ++_statistics._combinations;
    if (_options._presence_only &&
        _result.bits_found().find(bit) != _result.bits_found().end()) {
      return;
    }
    // Feature counts cannot realistically overflow for the molecules processed;
    // omit an overflow check here.
    _result.hit_bit(bit);
  }

  // Phase two: select one previously computed path hash from every physical arm.
  // No molecular traversal occurs here. `arm` is the next arm to select, `total`
  // is the sum of lengths already selected, and _selected[0..arm-1] holds their
  // hashes. Trying each length bucket and each hash in it generates the Cartesian
  // product of the arm choices, with recursion depth only three or four.
  //
  // Reserve at least one bond for each remaining arm when calculating the length
  // available to this arm. For example, three arms with a total budget of five
  // allow the first arm up to length three, leaving one bond for each other arm.
  // Valid length tuples include (1,1,1), (1,2,2), and (3,1,1), where paths exist.
  //
  // Once all arms have been selected, reject totals below the minimum and pass
  // the complete selection to RecordCombination. That method sorts the arm hashes
  // to remove arm-order dependence, mixes them with the centre type and arm count,
  // and records the bit. Different selections can increment the same bit; in
  // presence-only mode only the first such occurrence is recorded.
  void
  Combine(int arm, int total) {
    if (arm == static_cast<int>(_arms.size())) {
      if (total >= _options._min_total_length) {
        RecordCombination();
      }
      return;
    }
    const int remaining = _arms.size() - arm - 1;
    const int available = _options._max_total_length - total - remaining;
    for (int length = 1; length <= std::min(_limit, available); ++length) {
      for (const ArmPath& path : _arms[arm][length]) {
        _selected[arm] = path._hash;
        Combine(arm + 1, total + length);
      }
    }
  }

  // At radius r use paths of length r or shorter paths that have terminated.
  // At least one arm must reach r: otherwise this exact selection was already
  // emitted at an earlier radius. Actual path lengths still obey the total budget.
  void
  CombineAtRadius(int arm, int total, int radius, bool reaches_radius) {
    if (arm == static_cast<int>(_arms.size())) {
      if (reaches_radius && total >= _options._min_total_length) {
        RecordCombination();
      }
      return;
    }
    const int remaining = _arms.size() - arm - 1;
    const int available = _options._max_total_length - total - remaining;
    for (int length = 1; length <= std::min(radius, available); ++length) {
      for (const ArmPath& path : _arms[arm][length]) {
        if (length < radius && !path._terminal) {
          continue;
        }
        _selected[arm] = path._hash;
        CombineAtRadius(arm + 1, total + length, radius,
                        reaches_radius || length == radius);
      }
    }
  }

 public:
  Enumeration(const Molecule& m, const uint32_t* types, const Options& options,
              Sparse_Fingerprint_Creator& result, Statistics& statistics)
      : _m(m),
        _types(types),
        _options(options),
        _result(result),
        _statistics(statistics),
        _position(m.natoms(), -1) {
  }

  void
  Process(atom_number_t root) {
    const Atom& atom = _m[root];
    const int degree = atom.ncon();
    if (degree != 3 && degree != 4) {
      return;
    }
    if (_options._centre_degree && degree != _options._centre_degree) {
      return;
    }
    if (degree > _options._max_total_length) {
      return;
    }
    ++_statistics._centres;
    _limit = _options._max_total_length - degree + 1;
    if (_options._max_arm_length) {
      _limit = std::min(_limit, _options._max_arm_length);
    }
    // A simple path can visit each heavy atom once and then close a ring.
    _limit = std::min(_limit, _m.natoms());
    _arms.assign(degree, Arm(_limit + 1));
    _centre_type = _types[root];
    _position[root] = 0;
    for (int i = 0; i < degree; ++i) {
      const Bond* bond = atom[i];
      const atom_number_t next = bond->other(root);
      uint32_t hash = Mix(0x41524d31u, BondCode(*bond));  // ARM1 version seed.
      hash = Mix(Mix(hash, 0x41544f4du), _types[next]);
      _position[next] = 1;
      Extend(next, root, 1, hash, _arms[i]);
      _position[next] = -1;
    }
    _position[root] = -1;
    if (_options._synchronised) {
      for (int radius = 1; radius <= _limit; ++radius) {
        CombineAtRadius(0, 0, radius, false);
      }
    } else {
      Combine(0, 0);
    }
  }
};
}  // namespace

bool
Options::Valid() const {
  if (_min_total_length < 3 || _max_total_length < _min_total_length ||
      _max_arm_length < 0 ||
      (_centre_degree != 0 && _centre_degree != 3 && _centre_degree != 4)) {
    return false;
  }
  for (int degree : {3, 4}) {
    if (_centre_degree && degree != _centre_degree) {
      continue;
    }
    if (degree > _max_total_length) {
      continue;
    }
    if (!_max_arm_length ||
        static_cast<int64_t>(degree) * _max_arm_length >= _min_total_length) {
      return true;
    }
  }
  return false;
}

Generator::Generator() {
  _atom_typing.build("UST:Y");
}

bool
Generator::SetOptions(const Options& options) {
  if (!options.Valid()) {
    return false;
  }
  _options = options;
  return true;
}

int
Generator::SetAtomTyping(const const_IWSubstring& specification) {
  return _atom_typing.build(specification);
}

int
Generator::Fingerprint(Molecule& m, Sparse_Fingerprint_Creator& result, const int* roots,
                       Statistics* statistics) {
  result.clear();
  if (statistics) {
    *statistics = Statistics{};
  }
  std::vector<uint32_t> types(m.natoms());
  if (m.natoms() && !_atom_typing.assign_atom_types(m, types.data())) {
    return 0;
  }
  return Fingerprint(m, types.data(), result, roots, statistics);
}

int
Generator::Fingerprint(Molecule& m, const uint32_t* types,
                       Sparse_Fingerprint_Creator& result, const int* roots,
                       Statistics* statistics) const {
  result.clear();
  Statistics local;
  Statistics& stats = statistics ? *statistics : local;
  stats = Statistics{};
  if (!_options.Valid() || (m.natoms() && !types) || m.natoms(1) != 0) {
    return 0;
  }
  m.compute_aromaticity_if_needed();
  Enumeration enumeration(m, types, _options, result, stats);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    if (roots && !roots[atom]) {
      continue;
    }
    enumeration.Process(atom);
  }
  stats._distinct_bits = result.nbits();
  return 1;
}
}  // namespace multi_arm_fingerprint
