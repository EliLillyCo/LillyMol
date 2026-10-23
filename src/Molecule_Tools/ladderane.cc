// Called from substructure_demerits, but also incoroporated into
// the Lilly Medchem Rules.
// Detects ladderane like structures as well as consecutive spiro
// fused 3 membered rings.

#include <algorithm>
#include <memory>
#include <vector>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

namespace ladderane {

// Return true if `r` is all carbon atoms.
bool
AllCarbon(const Molecule& m, const Ring& r) {
  const int ring_size = r.number_elements();
  for (int i = 0; i < ring_size; ++i) {
    atom_number_t a = r[i];
    if (m.atomic_number(a) != 6) {
      return false;
    }
  }

  return true;
}

// Count the number of ladderane fusions in `r4` where the
// fused system identifier of the ring is `fsid`
int
CountLadderane(Molecule& m,
               int fsid,
               const std::vector<const Ring*>& r4) {
  const int nr = r4.size();
  int rc = 0;
  for (int i = 0; i < nr; ++i) {
    const Ring* ri = r4[i];
    if (ri->fused_system_identifier() != fsid) {
      continue;
    }
    for (int j = i + 1; j < nr; ++j) {
      const Ring* rj = r4[j];
      if (rj->fused_system_identifier() != fsid) {
        continue;
      }
      if (ri->is_fused_to(rj)) {
        ++rc;
      }
    }
  }

  return rc;
}

// If a ladderane is not found, return 0;
// If a ladderane is found, return the number of 4 membered rings
// in the ladderane.
// `r4` is a list of four membered, fused, carbocycles.
int
CountLadderane(Molecule& m,
               const std::vector<const Ring*>& r4) {
  const int nr = r4.size();

  extending_resizable_array<int> fsid_done;

  int rc = 0;

  for (int i = 0; i < nr; ++i) {
    const Ring* r = r4[i];
    const int fsid = r->fused_system_identifier();

    if (fsid_done[fsid]) {
      continue;
    }

    int fusions = CountLadderane(m, fsid, r4);
    if (fusions > rc) {
      rc = fusions;
    }

    fsid_done[fsid] = 1;
  }

  return rc;
}

int
CountLadderane(Molecule& m) {
  const int nr = m.nrings();
  if (nr < 2) {
    return 0;
  }

  // Must be at least two four membered rings.

  std::vector<const Ring*> r4;

  for (int i = 0; i < nr; ++i) {
    const Ring* r = m.ringi(i);

    if (r->size() < 4) {
      continue;
    }

    if (r->size() > 4) {
      break;
    }

    if (! r->is_fused()) {
      continue;
    }

    // Eliminate likely cubanes.
    // Could eliminate other things too, but that's fine.
    if (r->fused_ring_neighbours() > 2) {
      continue;
    }

    if (r->largest_number_of_bonds_shared_with_another_ring() > 1) {
      continue;
    }

    if (! AllCarbon(m, *r)) {
      continue;
    }

    r4.push_back(r);
  }

  if (r4.size() < 2) {
    return 0;
  }

  return CountLadderane(m, r4);
}

// Recursively count the number if connected atoms that are
// set in `in_system`.
// For each atom found, set the corresponding value in `in_system` to 2.
int
SystemSize(const Molecule& m, atom_number_t zatom,
           int* in_system) {
  in_system[zatom] = 2;
  int rc = 1;

  const Atom* a = m.atomi(zatom);
  const int acon = a->ncon();
  for (int i = 0; i < acon; ++i) {
    atom_number_t o = a->other(zatom, i);
    if (in_system[o] != 1) {
      continue;
    }
    rc += SystemSize(m, o, in_system);
  }

  return rc;
}

int
CountPolySpiroCycloPropane(Molecule& m, const std::vector<const Ring*>& r3) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> in_system = std::make_unique<int[]>(matoms);

  const int nr = r3.size();
  for (int i = 0; i < nr; ++i) {
    r3[i]->set_vector(in_system.get(), 1);
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (in_system[i] != 1) {
      continue;
    }

    const int atoms_found = SystemSize(m, i, in_system.get());
    // Ignore two rings joined together. Or should that also be a demerit?
    if (atoms_found <= 5) {
      continue;
    }

    // 5 atoms-> 2 rings, 7 atoms-> 3 rings
    const int sys_size = atoms_found / 2;

    if (sys_size > rc) {
      rc = sys_size;
    }
  }

  return rc;
}

// REturn true if `r` contains a likely point of spiro fusion.
int
ContainsLikelySpiro(Molecule& m, const Ring& r) {
  for (int i = 0; i < r.number_elements(); ++i) {
    atom_number_t a = r[i];
    const int nr = m.nrings(a);
    // Not interested in more complex fusions.
    if (nr > 2) {
      return 0;
    }

    if (nr == 1) {
      continue;
    }

    // At this stage nr == 2
    if (m.ncon(a) == 4) {
      return 1;
    }
  }

  return 0;
}

int
CountPolySpiroCycloPropane(Molecule& m) {
  const int nrings = m.nrings();
  if (nrings < 2) {
    return 0;
  }

  std::vector<const Ring*> r3;

  for (int i = 0; i < nrings; ++i) {
    const Ring* r = m.ringi(i);
    if (r->size() > 3) {
      break;
    }

    // Spiro fusions are not detected as ring fusions.
    if (r->is_fused()) {
      continue;
    }

    if (! AllCarbon(m, *r)) {
      continue;
    }

    if (!ContainsLikelySpiro(m, *r)) {
      continue;
    }

    r3.push_back(r);
  }
  
  if (r3.size() < 2) {
    return 0;
  }

  return CountPolySpiroCycloPropane(m, r3);
}

}  // namespace ladderane
