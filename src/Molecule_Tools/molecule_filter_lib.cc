#include <iostream>
#include <optional>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/xlogp.h"
#include "molecule_filter_lib.h"

namespace molecule_filter_lib {

using std::cerr;

int
MoleculeFilter::Build(IWString& fname) {
  std::optional<MoleculeFilterData::Requirements> maybe_proto =
                iwmisc::ReadTextProto<MoleculeFilterData::Requirements>(fname);
  if (! maybe_proto) {
    cerr << "MoleculeFilter::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  _requirements = *maybe_proto;

  InitialiseOptionalFeatures();

  return 1;
}

int
MoleculeFilter::Build(const MoleculeFilterData::Requirements& proto) {
  _requirements = proto;

  InitialiseOptionalFeatures();

  return 1;
}

void
MoleculeFilter::InitialiseOptionalFeatures() {
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  set_display_psa_unclassified_atom_mesages(0);
  xlogp::SetIssueUnclassifiedAtomMessages(0);

  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_apply_zwitterion_correction(1);
}

int
CountHeteroatoms(const Molecule& m) {
  int rc = 0;
  m.each_atom_lambda([&rc](const Atom& a) {
    if (a.atomic_number() != 6) {
      ++rc;
    }
  });

  return rc;
}

int
AromaticRingCount(Molecule& m) {
  m.compute_aromaticity_if_needed();

  int rc = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      ++rc;
    }
  }

  return rc;
}

// Return true if we examine multiple fragments, which
// means that `largest_frag` will be different from `smiles`.
bool
LargestFragment(const const_IWSubstring& smiles,
                const_IWSubstring& largest_frag,
                int& natoms, int& nrings) {
  int max_atoms = 0;
  int i = 0;
  const_IWSubstring token;
  int fragments_examined = 0;
  while (smiles.nextword(token, i, '.')) {
    ++fragments_examined;
    int nri;
    int nat = count_atoms_in_smiles(token, nri);
    if (nat > max_atoms) {
      max_atoms = nat;
      nrings = nri;
      largest_frag = token;
    }
  }
  
  natoms = max_atoms;

  if (fragments_examined == 1) {
    return false;
  }
  return true;
}

std::tuple<int, int>
MaxRingSystemSize(Molecule& m, std::unique_ptr<int[]>& tmp) {
  const int matoms = m.natoms();

  m.compute_aromaticity_if_needed();

  if (! tmp) {
    tmp.reset(new int[matoms]);
  }
  std::fill_n(tmp.get(), matoms, 0);

  const int nrings = m.nrings();

  std::unique_ptr<int[]> ring_already_done = std::make_unique<int[]>(nrings);
  std::fill_n(ring_already_done.get(), nrings, 0);

  int max_system_size = 0;
  int max_aromatic_rings_in_system = 0;
  for (int i = 0; i < nrings; ++i) {
    if (ring_already_done[i]) {
      continue;
    }
    const Ring* ri = m.ringi(i);
    if (! ri->is_fused()) {
      continue;
    }

    int system_size = 1;
    int aromatic_rings_in_system;
    if (ri->is_aromatic()) {
      aromatic_rings_in_system = 1;
    } else {
      aromatic_rings_in_system = 0;
    }


    for (int j = i + 1; j < nrings; ++j) {
      if (ring_already_done[j]) {
        continue;
      }

      ring_already_done[j] = 1;
      const Ring* rj = m.ringi(j);
      if (ri->fused_system_identifier() == rj->fused_system_identifier()) {
        ++system_size;
        if (rj->is_aromatic()) {
          ++aromatic_rings_in_system;
        }
      }
    }
    if (system_size > max_system_size) {
      max_system_size = system_size;
    }
    if (aromatic_rings_in_system > max_aromatic_rings_in_system) {
      max_aromatic_rings_in_system = aromatic_rings_in_system;
    }
  }

  return std::make_tuple(max_system_size, max_aromatic_rings_in_system);
}

// Lifted from iwdescr.cc
void
RuleOfFive(Molecule & m, int& acceptor, int& donor) {
  acceptor = 0;
  donor = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    atomic_number_t z = m.atomic_number(i);
    // Intercept the most common case.
    if (z == 6) {
      continue;
    }

    if (z == 7 || z == 8) {
    } else {
      continue;
    }

    ++acceptor;

    const int h = m.hcount(i);

    // acceptor
    if (0 == h) {
      continue;
    }

    if (7 == z && h > 1) {
      donor += 2;
    } else {
      donor += 1;
    }
  }
}

int
HalogenCount(const Molecule& m) {
  static std::vector<int> halogen = {
    0,  // 0
    0,  // 1
    0,  // 2
    0,  // 3
    0,  // 4
    0,  // 5
    0,  // 6
    0,  // 7
    0,  // 8
    1,  // 9
    0,  // 10
    0,  // 11
    0,  // 12
    0,  // 13
    0,  // 14
    0,  // 15
    0,  // 16
    1,  // 17
    0,  // 18
    0,  // 19
    0,  // 20
    0,  // 21
    0,  // 22
    0,  // 23
    0,  // 24
    0,  // 25
    0,  // 26
    0,  // 27
    0,  // 28
    0,  // 29
    0,  // 30
    0,  // 31
    0,  // 32
    0,  // 33
    0,  // 34
    0,  // 35
    0,  // 36
    1,  // 37
    0,  // 38
    0,  // 39
    0,  // 40
    0,  // 41
    0,  // 42
    0,  // 43
    0,  // 44
    0,  // 45
    0,  // 46
    0,  // 47
    0,  // 48
    0,  // 49
    0,  // 50
    0,  // 51
    0,  // 52
    1   // 53
  };

  int rc = 0;

  for (const Atom* a : m) {
    const uint32_t z = a->atomic_number();
    if (z < halogen.size()) {
      rc += halogen[z];
    }
  }

  return rc;
}

int 
Sp3Carbon(Molecule & m) {
  int rc = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.saturated(i)) {
      ++rc;
    }
  }

  return rc;
}


int
MoleculeFilter::Ok(Molecule& m) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    return 0;
  }

  if (_requirements.has_min_natoms() && matoms < _requirements.min_natoms()) {
    return 0;
  }

  if (_requirements.has_max_natoms() && matoms > _requirements.max_natoms()) {
    return 0;
  }

  const int nrings = m.nrings();
  if (_requirements.has_min_nrings() && nrings < _requirements.min_nrings()) {
    return 0;
  }

  if (_requirements.has_max_nrings() && nrings > _requirements.max_nrings()) {
    return 0;
  }

  if (_requirements.has_min_heteroatom_count() ||
      _requirements.has_min_heteroatom_fraction() ||
      _requirements.has_max_heteroatom_fraction()) {
    const int hac = CountHeteroatoms(m);

    if (_requirements.has_min_heteroatom_count() && hac < _requirements.min_heteroatom_count()) {
      return 0;
    }

    if (_requirements.has_min_heteroatom_fraction() ||
        _requirements.has_max_heteroatom_fraction()) {
      const float haf = iwmisc::Fraction<float>(hac, matoms);
      if (_requirements.has_min_heteroatom_fraction() && haf < _requirements.min_heteroatom_fraction()) {
        return 0;
      }
      if (_requirements.has_max_heteroatom_fraction() && haf > _requirements.max_heteroatom_fraction()) {
        return 0;
      }
    }
  }

  if (_requirements.has_exclude_non_organic() && ! m.organic_only()) {
    return 0;
  }

  if (_requirements.has_exclude_isotopes() && m.number_isotopic_atoms() > 0) {
    return 0;
  }

  int arc = 0;
  int need_to_compute_aromatic_rings = 0;
  if (_requirements.has_min_aromatic_ring_count() ||
      _requirements.has_max_aromatic_ring_count() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    need_to_compute_aromatic_rings = 1;
  }

  if (need_to_compute_aromatic_rings) {
    // Check nrings first. If not enough rings if all were aromatic...
    if (_requirements.has_min_aromatic_ring_count() && nrings < _requirements.min_aromatic_ring_count()) {
      return 0;
    }

    arc = AromaticRingCount(m);
    if (_requirements.has_min_aromatic_ring_count() && arc < _requirements.min_aromatic_ring_count()) {
      return 0;
    }
    if (_requirements.has_max_aromatic_ring_count() && arc > _requirements.max_aromatic_ring_count()) {
      return 0;
    }
  }

  if (_requirements.has_min_aliphatic_ring_count() ||
      _requirements.has_max_aliphatic_ring_count()) {
    if (_requirements.has_min_aliphatic_ring_count() && nrings < _requirements.min_aliphatic_ring_count()) {
      return 0;
    }
    const int alring = nrings - AromaticRingCount(m);
    if (_requirements.has_min_aliphatic_ring_count() && alring < _requirements.min_aliphatic_ring_count()) {
      return 0;
    }
    if (_requirements.has_max_aliphatic_ring_count() && alring > _requirements.max_aliphatic_ring_count()) {
      return 0;
    }
  }

  if (_requirements.has_min_hba() || _requirements.has_max_hba() ||
      _requirements.has_min_hbd() || _requirements.has_max_hbd()) {
    int hba, hbd;
    RuleOfFive(m, hba, hbd);
    if (_requirements.has_min_hba() && hba < _requirements.min_hba()) {
      return 0;
    }
    if (_requirements.has_max_hba() && hba > _requirements.max_hba()) {
      return 0;
    }
    if (_requirements.has_min_hbd() && hbd < _requirements.min_hbd()) {
      return 0;
    }
    if (_requirements.has_max_hbd() && hbd > _requirements.max_hbd()) {
      return 0;
    }
  }

  if (_requirements.has_min_sp3_carbon()) {
    const int csp3 = Sp3Carbon(m);
    if (csp3 < _requirements.min_sp3_carbon()) {
      return 0;
    }
  }

  if (_requirements.has_max_halogen_count()) {
    const int h = HalogenCount(m);
    if (h > _requirements.max_halogen_count()) {
      return 0;
    }
  }

  if (_requirements.has_largest_ring_size()) {
    int rsze = m.ringi(nrings - 1)->number_elements();
    if (rsze > _requirements.largest_ring_size()) {
      return 0;
    }
  }

  if (_requirements.has_min_rotatable_bonds() || _requirements.has_max_rotatable_bonds()) {
    int rotb = _rotbond.Process(m);
    if (_requirements.has_min_rotatable_bonds() && rotb < _requirements.min_rotatable_bonds()) {
      return 0;
    }
    if (_requirements.has_max_rotatable_bonds() && rotb > _requirements.max_rotatable_bonds()) {
      return 0;
    }
  }

  if (_requirements.has_max_aromatic_density()) {
    const float aromdens = iwmisc::Fraction<float>(m.aromatic_atom_count(), matoms);
    if (aromdens > _requirements.max_aromatic_density()) {
      return 0;
    }
  }

  if (_requirements.has_max_distance() && matoms > _requirements.max_distance()) {
    const int d = m.longest_path();
    if (d > _requirements.max_distance()) {
      return 0;
    }
  }

  if (_requirements.has_min_tpsa() || _requirements.has_max_tpsa()) {
    float tpsa = novartis_polar_surface_area(m);
    if (_requirements.has_min_tpsa() && tpsa < _requirements.min_tpsa()) {
      return 0;
    }
    if (_requirements.has_max_tpsa() && tpsa > _requirements.max_tpsa()) {
      return 0;
    }
  }

  // A temporary array that some external functions might need.
  std::unique_ptr<int[]> tmp;

  if (_requirements.has_max_ring_system_size() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    if (nrings < _requirements.max_ring_system_size() &&
        nrings < _requirements.max_aromatic_rings_in_system()) {
      // no need to compute.
    }  else {
      const auto [max_ring_system_size, max_aromatic_rings_in_system] = MaxRingSystemSize(m, tmp);
      if (_requirements.has_max_ring_system_size() && max_ring_system_size > _requirements.max_ring_system_size()) {
        return 0;
      }
      if (_requirements.has_max_aromatic_rings_in_system() &&
          max_aromatic_rings_in_system > _requirements.max_aromatic_rings_in_system()) {
      }
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_alogp()) {
    std::optional<double> x = _alogp.LogP(m);
    if (! x) {
    } else if (_requirements.has_min_alogp() && *x < _requirements.min_alogp()) {
      return 0;
    }
    if (!x) {
    } else if (_requirements.has_max_alogp() && *x > _requirements.max_alogp()) {
      return 0;
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_xlogp()) {
    if (! tmp) {
      tmp.reset(new int[matoms]);
    }
    std::fill_n(tmp.get(), matoms, 0);
    std::optional<double> x = xlogp::XLogP(m, tmp.get());
    if (! x) {
    } else if (_requirements.has_min_xlogp() && *x < _requirements.min_xlogp()) {
      return 0;
    }
    if (!x) {
    } else if (_requirements.has_max_xlogp() && *x > _requirements.max_xlogp()) {
      return 0;
    }
  }

  return 1;
}


}  // namespace molecule_filter_lib
