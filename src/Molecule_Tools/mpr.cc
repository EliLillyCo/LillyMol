#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"

#include "mpr.h"

using std::cerr;

Molecular_Properties_Generator::Molecular_Properties_Generator ()
{
  _unsaturation_includes_aromatic = 0;

  return;
}

template <typename T>
int
Molecular_Properties_Generator::_molecular_properties_generation (Molecule & m,
                                 const int * ring_membership,
                                 T * properties) const
{
  int matoms = m.natoms();

  if (matoms > std::numeric_limits<T>::max())    // check for overflow. Note that if natoms overflows, so will many of the others. We don't check this!!! - too lazy
  {
    properties[0] = std::numeric_limits<T>::max();
  }
  else
    properties[0] = matoms;

  int nr = m.nrings();

  properties[1] = 0;
  properties[2] = nr;

  if (nr)
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * r = m.ringi(i);
      if (r->number_elements() > properties[1])
        properties[1] = r->number_elements();
    }
  }

  int ring_atoms = 0;
  int aromatic_atoms = 0;
  int fused_ring_atoms = 0;

  if (nr)    // the molecule has rings
  {
    assert (nullptr != ring_membership);

    for (int i = 0; i < matoms; i++)
    {
      if (0 == ring_membership[i])
        continue;

      ring_atoms++;
      if (m.is_aromatic(i))
        aromatic_atoms++;

      if (nr > 1 && ring_membership[i] > 1)
        fused_ring_atoms++;
    }
  }
  
  properties[3] = ring_atoms;
  properties[4] = aromatic_atoms;
  properties[5] = fused_ring_atoms;

  int heteroatom_count = 0;
  int unsaturation = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      heteroatom_count++;

    if (a->nbonds() > a->ncon())
    {
      if (_unsaturation_includes_aromatic)
        unsaturation++;
      else if (nr && ring_membership[i] && m.is_aromatic(i))
        ;
      else
        unsaturation++;
    }
  }

  properties[6] = heteroatom_count;
  properties[7] = unsaturation;

//#define DEBUG_PROPERTIES
#ifdef DEBUG_PROPERTIES
  cerr << "Molecule '" << m.name() << "' properties";
  for (int i = 0; i < NPROPERTIES; i++)
  {
    cerr << ' ' << int(properties[i]);
  }

  cerr << endl;
#endif

  return 1;
}


template <typename T>
int
Molecular_Properties_Generator::operator() (Molecule & m,
                                 T * properties) const
{
  int nr = m.nrings();

  if (0 == nr)
    return _molecular_properties_generation(m, nullptr, properties);

  int * ring_membership = new int[m.natoms()]; std::unique_ptr<int[]> free_ring_membership(ring_membership);

  m.ring_membership_including_non_sssr_rings(ring_membership);

#ifdef ECHO_RING_MEMBERSHIP
  for (int i = 0; i < m.natoms(); i++)
  {
    if (ring_membership[i])
      cerr << "Atom " << i << " in " << ring_membership[i] << " rings, '" << m.smarts_equivalent_for_atom(i) << "'\n";
  }
#endif

  return _molecular_properties_generation(m, ring_membership, properties);
}

template int Molecular_Properties_Generator::operator()(Molecule &, unsigned char *) const;
template int Molecular_Properties_Generator::operator()(Molecule &, int *) const;

namespace mpr {

int
LargestRingSize(Molecule& m) {
  if (m.nrings() == 0) {
    return 0;
  }

  int rc = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->number_elements() > rc) {
      rc = r->number_elements();
    }
  }

  return rc;
}

int
RingAtoms(Molecule& m) {
  if (m.nrings() == 0) {
    return 0;
  }

  int rc = 0;
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i)) {
      ++rc;
    }
  }
  
  return rc;
}

// We adopt a fairly loose definition here. Note that it
// includes spiro fused atoms.
int
FusedRingAtoms(Molecule& m) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) > 2) {
      ++rc;
    }
  }

  return rc;
}

int
HeteroatomCount(const Molecule& m) {
  int rc = 0;
  for (const Atom* a : m) {
    if (a->atomic_number() != 6) {
      ++rc;
    }
  }

  return rc;
}

int
Unsaturation(Molecule& m) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.is_aromatic(i)) {
      continue;
    }

    if (! m.saturated(i)) {
      ++rc;
    }
  }

  return rc;
}

void
AddDoublyBonded(Molecule& m,
                int* spinach) {
  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    if (! b->is_double_bond()) {
      continue;
    }
    if (b->is_aromatic()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (spinach[a1]) {
      if (! spinach[a2] && m.ncon(a1) == 1) {
        spinach[a1] = 0;
      }
    }
    if (spinach[a2]) {
      if (! spinach[a1] && m.ncon(a2) == 1) {
        spinach[a2] = 0;
      }
    }
  }
}

int
ScaffoldAtoms(Molecule& m,
              std::unique_ptr<int[]>& spinach) {
  const int matoms = m.natoms();

  if (! spinach) {
    spinach.reset(new int[matoms]);
    m.identify_spinach(spinach.get());
    AddDoublyBonded(m, spinach.get());
  }

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (! spinach[i]) {
      ++rc;
    }
  }

  return rc;
}

int
LinkerAtoms(Molecule& m,
            std::unique_ptr<int[]>& spinach) {
  const int matoms = m.natoms();

  if (! spinach) {
    spinach.reset(new int[matoms]);
    m.identify_spinach(spinach.get());
    AddDoublyBonded(m, spinach.get());
  }

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (spinach[i]) {
      continue;
    }

    if (m.ring_bond_count(i)) {
      continue;
    }

    ++rc;
  }

  return rc;
}

int
LargestRingSys(Molecule& m) {
  if (m.nrings() == 0) {
    return 0;
  }

  extending_resizable_array<int> count;

  for (const Ring* r : m.sssr_rings()) {
    const int fsid = r->fused_system_identifier();
    if (fsid >= 0) {
      ++count[r->fused_system_identifier()];
    }
  }

  return *std::max_element(count.begin(), count.end());
}

int
RotBonds(Molecule& m) {
  static quick_rotbond::QuickRotatableBonds simple_rotbond;
  simple_rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);

  return simple_rotbond.Process(m);
}

int
GenerateMolecularProperties(Molecule& m,
                const resizable_array<mpr::Feature>& features,
                float* result) {
  const int n = features.number_elements();

  const int matoms = m.natoms();
  if (matoms == 0) {
    std::fill_n(result, n, 0.0f);
    return 1;
  }

  // Multiple properties need the scaffold identified. The first
  // one will set this.
  std::unique_ptr<int[]> spinach;

  for (int i = 0; i < n; ++i) {
    switch (features[i]) {
      case mpr::kNatoms:
        result[i] = matoms;
        break;
      case mpr::kLargestRingSize:
        result[i] = LargestRingSize(m);
        break;
      case mpr::kNrings:
        result[i] = m.nrings();
        break;
      case mpr::kRingAtoms:
        result[i] = RingAtoms(m);
        break;
      case mpr::kAromaticAtoms:
        result[i] = m.aromatic_atom_count();
        break;
      case mpr::kFusedRingAtoms:
        result[i] = FusedRingAtoms(m);
        break;
      case mpr::kHeteroatomCount:
        result[i] = HeteroatomCount(m);
        break;
      case mpr::kUnsaturation:
        result[i] = Unsaturation(m);
        break;
      case mpr::kMaxDist:
        result[i] = m.longest_path();
        break;
      case mpr::kScaffoldAtoms:
        result[i] = ScaffoldAtoms(m, spinach);
        break;
      case mpr::kLinkerAtoms:
        result[i] = LinkerAtoms(m, spinach);
        break;
      case mpr::kLargestRingSys:
        result[i] = LargestRingSys(m);
        break;
      case mpr::kRotBonds:
        result[i] = RotBonds(m);
        break;
      default:
        cerr << "GenerateMolecularProperties:what feature " << features[i] << '\n';
        return 0;
    }
  }

  return 1;
}

void
DisplayOptions(char flag, std::ostream& output) {
  output << " -" << flag << " none      no optional properties\n";
  output << " -" << flag << " natoms    number of atoms\n";
  output << " -" << flag << " lrsz      largest ring size\n";
  output << " -" << flag << " nrings    number of rings\n";
  output << " -" << flag << " ringatom  number of ring atoms\n";
  output << " -" << flag << " aroma     number of aromatic atoms\n";
  output << " -" << flag << " frac      number of fused ring atoms\n";
  output << " -" << flag << " htroatom  number of heteroatoms\n";
  output << " -" << flag << " unsatura  number of unsaturated (non aromatic)\n";
  output << " -" << flag << " maxdist   longest path property (maximum distance)\n";
  output << " -" << flag << " scaffold  atoms in scaffold property\n";
  output << " -" << flag << " lrsysz    number of rings in the largest ring system\n";
  output << " -" << flag << " linkatom  number of linker atoms - chain atoms between rings\n";
  output << " -" << flag << " rotbond   number of rotatable bonds\n";

  ::exit(0);
}

int
MolecularPropertiesGenerator::Initialise(Command_Line& cl, char flag) {
  const_IWSubstring o;
  for (int i = 0; cl.value(flag, o, i); ++i) {
    int j = 0;
    const_IWSubstring token;
    while (o.nextword(token, j, ',')) {
      if (token == "none") {
        _features.resize_keep_storage(0);
      } else if (token == "natoms") {
        _features << mpr::kNatoms;
      } else if (token == "lrsz") {
        _features << mpr::kLargestRingSize;
      } else if (token == "nrings") {
        _features << mpr::kNrings;
      } else if (token == "ringatom") {
        _features << mpr::kRingAtoms;
      } else if (token == "aroma") {
        _features << mpr::kAromaticAtoms;
      } else if (token == "frac") {
        _features << mpr::kFusedRingAtoms;
      } else if (token == "htroatom") {
        _features << mpr::kHeteroatomCount;
      } else if (token == "unsatura") {
        _features << mpr::kUnsaturation;
      } else if (token == "mxdst") {
        _features << mpr::kMaxDist;
      } else if (token == "scaffold") {
        _features << mpr::kScaffoldAtoms;
      } else if (token == "lrsysz") {
        _features << mpr::kLargestRingSys;
      } else if (token == "linkatom") {
        _features << mpr::kLinkerAtoms;
      } else if (token == "rotbond") {
        _features << mpr::kRotBonds;
      } else if (token == "all") {
        _features << mpr::kNatoms;
        _features << mpr::kLargestRingSize;
        _features << mpr::kNrings;
        _features << mpr::kRingAtoms;
        _features << mpr::kAromaticAtoms;
        _features << mpr::kFusedRingAtoms;
        _features << mpr::kHeteroatomCount;
        _features << mpr::kUnsaturation;
        _features << mpr::kMaxDist;
        _features << mpr::kScaffoldAtoms;
        _features << mpr::kLargestRingSys;
        _features << mpr::kLinkerAtoms;
        _features << mpr::kRotBonds;
      } else {
        cerr << "Unrecognised -" << flag << " qualifier '" << token << "'\n";
        DisplayOptions(flag, cerr);
      }
    }
  }

  return 1;
}

void
MolecularPropertiesGenerator::AddDefaultFeatures() {
  _features << mpr::kNatoms;
  _features << mpr::kLargestRingSize;
  _features << mpr::kNrings;
  _features << mpr::kRingAtoms;
  _features << mpr::kAromaticAtoms;
  _features << mpr::kFusedRingAtoms;
  _features << mpr::kHeteroatomCount;
  _features << mpr::kUnsaturation;
}

template <typename T>
int
MolecularPropertiesGenerator::GenerateMolecularProperties(Molecule& m,
                T* result) {
  const int n = _features.number_elements();

  const int matoms = m.natoms();
  if (matoms == 0) {
    std::fill_n(result, n, 0.0f);
    return 1;
  }

  if (_features.empty()) {
    AddDefaultFeatures();
  }

  // Multiple properties need the scaffold identified. The first
  // one will set this.
  std::unique_ptr<int[]> spinach;

  for (int i = 0; i < n; ++i) {
    switch (_features[i]) {
      case mpr::kNatoms:
        result[i] = matoms;
        break;
      case mpr::kLargestRingSize:
        result[i] = LargestRingSize(m);
        break;
      case mpr::kNrings:
        result[i] = m.nrings();
        break;
      case mpr::kRingAtoms:
        result[i] = RingAtoms(m);
        break;
      case mpr::kAromaticAtoms:
        result[i] = m.aromatic_atom_count();
        break;
      case mpr::kFusedRingAtoms:
        result[i] = FusedRingAtoms(m);
        break;
      case mpr::kHeteroatomCount:
        result[i] = HeteroatomCount(m);
        break;
      case mpr::kUnsaturation:
        result[i] = Unsaturation(m);
        break;
      case mpr::kMaxDist:
        result[i] = m.longest_path();
        break;
      case mpr::kScaffoldAtoms:
        result[i] = ScaffoldAtoms(m, spinach);
        break;
      case mpr::kLinkerAtoms:
        result[i] = LinkerAtoms(m, spinach);
        break;
      case mpr::kLargestRingSys:
        result[i] = LargestRingSys(m);
        break;
      case mpr::kRotBonds:
        result[i] = RotBonds(m);
        break;
      default:
        cerr << "GenerateMolecularProperties:what feature " << _features[i] << '\n';
        return 0;
    }
  }

  return 1;
}

template
int MolecularPropertiesGenerator::GenerateMolecularProperties(Molecule& m, int*);
template
int MolecularPropertiesGenerator::GenerateMolecularProperties(Molecule& m, uint8_t*);

}  // namespace mpr
