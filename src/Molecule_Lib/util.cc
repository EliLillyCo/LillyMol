#include <iostream>
#include <queue>
#include <vector>

#include "util.h"

#include "atom.h"
#include "bond.h"
#include "istream_and_type.h"
#include "molecule.h"

namespace lillymol {

using std::cerr;

std::optional<uint64_t>
CountMoleculesInFile(const IWString& fname, FileType ftype) {
  if (ftype == FILE_TYPE_INVALID) {
    ftype = discern_file_type_from_name(fname);
    if (ftype == FILE_TYPE_INVALID) {
      cerr << "CountMoleculesInFile:cannot discern file type '" << fname << "'\n";
      return std::nullopt;
    }
  }

  data_source_and_type<Molecule> input(ftype, fname);
  if (! input.good()) {
    return std::nullopt;
  }

  return input.molecules_remaining();
}

int
SetAtomsWithinRadius(const Molecule& m, const Set_of_Atoms& seeds, int radius, int flag,
                     int* destination) {
  if (radius < 0) {
    return 0;
  }

  const int matoms = m.natoms();
  if (matoms == 0) {
    return 0;
  }

  std::vector<int> distance(matoms, -1);
  std::queue<atom_number_t> to_process;

  int changed = 0;
  for (atom_number_t atom : seeds) {
    if (atom < 0 || atom >= matoms) {
      continue;
    }

    if (distance[atom] < 0) {
      distance[atom] = 0;
      to_process.push(atom);
    }
    if (destination[atom] != flag) {
      destination[atom] = flag;
      ++changed;
    }
  }

  while (! to_process.empty()) {
    const atom_number_t atom = to_process.front();
    to_process.pop();

    if (distance[atom] >= radius) {
      continue;
    }

    const int next_distance = distance[atom] + 1;
    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (distance[other] >= 0) {
        continue;
      }

      distance[other] = next_distance;
      to_process.push(other);
      if (destination[other] != flag) {
        destination[other] = flag;
        ++changed;
      }
    }
  }

  return changed;
}

}  // namespace lillymol
