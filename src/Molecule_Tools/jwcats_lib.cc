#include "Molecule_Tools/jwcats_lib.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Molecule_Lib/substructure.h"

namespace jwcats {

const char* const kPairName[kPairsPerDistance] = {"AA", "AD", "AP", "AN", "AH",
                                                  "DD", "DP", "DN", "DH", "PP",
                                                  "PN", "PH", "NN", "NH", "HH"};

int
PairNumber(int p1, int p2) {
  assert((p1 > -1) && (p1 < kNumberProperties) && (p2 > -1) && (p2 < kNumberProperties));

  int small_no = p1;
  int large_no = p2;

  if (small_no > large_no) {
    small_no = p2;
    large_no = p1;
  }

  return small_no * kNumberProperties - small_no * (small_no - 1) / 2 + large_no -
         small_no;
}

JWCats::JWCats() {
}

void
JWCats::BuildWriteArray() {
  _array_size = kPairsPerDistance * _max_bond_separation;
  _write_array_value.resize(_array_size);
  std::fill(_write_array_value.begin(), _write_array_value.end(), 1);

  if (_include_hydrophobic_pairs) {
    return;
  }

  for (int i = 0; i < _array_size; ++i) {
    if (i % kPairsPerDistance == kPairsPerDistance - 1) {
      _write_array_value[i] = 0;
    }
  }
}

void
JWCats::SetIncludeHydrophobicPairs(int s) {
  _include_hydrophobic_pairs = s;
  _initialised = 0;
}

void
JWCats::SetMaxBondSeparation(int s) {
  _max_bond_separation = s;
  _initialised = 0;
}

int
JWCats::Initialise() {
  if (_max_bond_separation < 1) {
    return 0;
  }
  if (_min_bond_separation < 0) {
    return 0;
  }
  if (_min_bond_separation > _max_bond_separation) {
    return 0;
  }
  if (_scaling_type < 0 || _scaling_type > 3) {
    return 0;
  }

  BuildWriteArray();
  _initialised = 1;

  return 1;
}

std::vector<std::string>
JWCats::FeatureNames() const {
  std::vector<std::string> result;
  if (!_initialised) {
    return result;
  }
  result.reserve(_array_size);

  for (int i = 0; i < _max_bond_separation; ++i) {
    for (int j = 0; j < kPairsPerDistance; ++j) {
      const int ndx = i * kPairsPerDistance + j;
      if (!_write_array_value[ndx]) {
        continue;
      }
      result.push_back("jwc_B" + std::to_string(i + 1) + "P" + kPairName[j]);
    }
  }

  return result;
}

ComputeStatus
JWCats::Compute(Molecule& m, Result& result) {
  if (!_initialised) {
    return ComputeStatus::kNotInitialised;
  }

  m.recompute_distance_matrix();

  if (_charge_assigner.active()) {
    _charge_assigner.process(m);
  }

  if (_make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  const int n_atoms = m.natoms();
  result.scaled_counts.resize(_array_size);
  std::fill(result.scaled_counts.begin(), result.scaled_counts.end(), 0.0);
  std::fill(result.property_count.begin(), result.property_count.end(), 0);
  result.number_heavy_atoms = 0;

  std::vector<int> donor_acceptor(n_atoms, 0);
  if (_donor_acceptor_assigner.active()) {
    _donor_acceptor_assigner.process(m, donor_acceptor.data());
  }

  std::array<std::vector<int>, kNumberProperties> properties;
  for (int i = 0; i < kNumberProperties; ++i) {
    properties[i].resize(n_atoms);
    std::fill(properties[i].begin(), properties[i].end(), 0);
  }

  for (int i = 0; i < n_atoms; ++i) {
    const Atom* atom_i = m.atomi(i);

    if (atom_i->atomic_number() == 1) {
      continue;
    }

    ++result.number_heavy_atoms;

    switch (donor_acceptor[i]) {
      case 0:
        break;
      case 1:
        properties[0][i] = 1;
        break;
      case 2:
        properties[0][i] = 1;
        properties[1][i] = 1;
        break;
      case 3:
        properties[1][i] = 1;
        break;
      default:
        break;
    }

    formal_charge_t fci = atom_i->formal_charge();
    if (fci == 1) {
      properties[2][i] = 1;
    } else if (fci == -1) {
      properties[3][i] = 1;
    }
  }

  if (_use_queries_to_determine_hydrophobicity) {
    const int nq = _queries.number_elements();
    for (int i = 0; i < nq; ++i) {
      Substructure_Results sresults;
      const int nhits = _queries[i]->substructure_search(m, sresults);
      if (nhits == 0) {
        continue;
      }

      for (int j = 0; j < nhits; ++j) {
        const Set_of_Atoms* embedding = sresults.embedding(j);
        embedding->set_vector(properties[4].data(), 1);
      }
    }
  } else {
    if (!m.compute_Gasteiger_partial_charges()) {
      return ComputeStatus::kMissingChargeData;
    }

    for (int i = 0; i < n_atoms; ++i) {
      if (std::fabs(m.charge_on_atom(i)) <= 0.2) {
        properties[4][i] = 1;
      }
    }
  }

  if (_scaling_type) {
    for (int i = 0; i < n_atoms; ++i) {
      for (int j = 0; j < kNumberProperties; ++j) {
        if (properties[j][i]) {
          ++result.property_count[j];
        }
      }
    }
  }

  for (int i = 0; i < n_atoms; ++i) {
    if (m.atomic_number(i) == 1) {
      continue;
    }

    for (int j = i + 1; j < n_atoms; ++j) {
      if (m.atomic_number(j) == 1) {
        continue;
      }

      const int bond_distance = m.bonds_between(i, j);

      if (bond_distance > _max_bond_separation) {
        continue;
      }
      if (bond_distance < _min_bond_separation) {
        continue;
      }

      for (int p1 = 0; p1 < kNumberProperties; ++p1) {
        if (properties[p1][i] == 0) {
          continue;
        }

        for (int p2 = 0; p2 < kNumberProperties; ++p2) {
          if (properties[p2][j] == 0) {
            continue;
          }

          const int property_pair_no = PairNumber(p1, p2);
          const int array_column =
              (bond_distance - 1) * kPairsPerDistance + property_pair_no;

          if (_scaling_type == 0) {
            result.scaled_counts[array_column] += 1.0;
          } else if (_scaling_type == 1) {
            result.scaled_counts[array_column] +=
                1.0 / static_cast<double>(result.number_heavy_atoms);
          } else if (_scaling_type == 2) {
            result.scaled_counts[array_column] +=
                1.0 / (result.property_count[p1] + result.property_count[p2]);
          } else if (_scaling_type == 3) {
            result.scaled_counts[array_column] +=
                static_cast<double>(result.number_heavy_atoms) /
                (result.property_count[p1] + result.property_count[p2]);
          } else {
            return ComputeStatus::kError;
          }
        }
      }
    }
  }

  return ComputeStatus::kOk;
}

}  // namespace jwcats
