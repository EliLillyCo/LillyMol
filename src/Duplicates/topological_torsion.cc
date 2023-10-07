// Generate topological torsion fingerprints.

#include <algorithm>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "topological_torsion.h"

namespace topological_torsion {

constexpr int kMaxPathLength = 4;

// Return an int describing the number of bonds in `b`.
// Does not really matter what the numbers are.
int
NumberBonds(const Bond* b) {
  if (b->is_aromatic()) {
    return 4;
  }
  return b->number_of_bonds();
}

// Struct to avoid passing too many arguments around.
// Perhaps this should be doing more of the computations.
struct Args {
  // The current atom being added to the torsion.
  atom_number_t current_atom;

  // The atom in the path before `current_atom`.
  // Avoids backtracking.
  atom_number_t previous_atom;

  // First atom in the torsion - used to detect rings.
  atom_number_t first_atom;

  // Whether or not an atom can be included in a path.
  const int * include_atom;

  // The atom type for each atom
  const atom_type_t* atom_type;

  // Will be set to non zero for atoms already in the path
  int * visited;

  // The atom types of the atoms in the path.
  atom_type_t path[4];

  // The NumberBonds value for each Bond in the path.
  // Dimensioned 4 to handle 4 membered rings, normal
  // torsions only require 3 bonds.
  int bonds_in_path[4];

  // The current path length.
  int path_length = 0;

  // Private functions
  private:

  public:
    // Given 4 atoms and 3 bonds, return a canonical bit number for that.
    uint32_t CanonicalLinear() const;

    // The path length is 4, maybe it forms a 4 membered ring.
    // If atoms 0 and 3 are bonded, add that bond to bonds_in_path,
    // and return true.  Otherwise false.
    bool AddRingClosureBond(const Molecule& m);

    // Return canonical bit numbers for either a 3 membered ring or
    // a 4 membered ring.
    uint32_t CanonicalRingBit3() const;
    uint32_t CanonicalRingBit4() const;

  // For debugging
  void DebugPrint(std::ostream& output) const;
};

void
Args::DebugPrint(std::ostream& output) const {
  output << "Path len " << path_length;
  for (int i = 0; i < path_length; ++i) {
    if (i > 0) {
      output <<  " [" << bonds_in_path[i - 1] << "] ";
    }
    output << ' ' << path[i];
  }
  output << '\n';
  output << "current_atom " << current_atom << '\n';
}

// If there is a bond between the first and last atom in path,
// then add that bond number to bonds_in_path.
bool
Args::AddRingClosureBond(const Molecule& m) {
  const Bond * b = m.bond_between_atoms_if_present(first_atom, current_atom);
  if (b == nullptr) {
    return false;
  }

  bonds_in_path[3] = NumberBonds(b);
  return true;
}

// This should probably not exist, but should be inlined anyway.
void
HitBit(uint32_t bitnum, Sparse_Fingerprint_Creator& destination) {
  destination.hit_bit(bitnum);
}

//#define DEBUG_TOPOLOGICAL_TORSION

// Firm the bit in both directions, and just return the largest value.
uint32_t
Args::CanonicalLinear() const {
  static uint32_t magic[] {73, 46, 93};
  const uint32_t forward = path[0] + magic[0] * bonds_in_path[0] * path[1] +
                                     magic[1] * bonds_in_path[1] * path[2] +
                                     magic[2] * bonds_in_path[2] * path[3];
  const uint32_t reverse = path[3] + magic[0] * bonds_in_path[2] * path[2] +
                                     magic[1] * bonds_in_path[1] * path[1] +
                                     magic[2] * bonds_in_path[0] * path[0];
  return std::max(forward, reverse);
}

// Update `current_max` if `new_value` is larger.
void
MaybeUpdateHighest(uint32_t new_value, uint32_t& current_max) {
  if (new_value > current_max) {
    current_max = new_value;
  }
}

// A 3 membered ring is formed. There
// are 3 valid atoms and 3 valid bonds in `this`.
// There are 3 starting points, and two directions, so 6
// possible traversals of the ring. Return the max bit number.
uint32_t
Args::CanonicalRingBit3() const {
  static uint32_t magic[] {73, 46, 93};

  uint32_t highest = 0;
  // First try the clockwise directions
  uint32_t tmp;
  tmp = path[0] + magic[0] * bonds_in_path[0] * path[1] +
                  magic[1] * bonds_in_path[1] * path[2] +
                  magic[2] * bonds_in_path[2] * path[0];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[1] + magic[0] * bonds_in_path[1] * path[2] +
                  magic[1] * bonds_in_path[2] * path[0] +
                  magic[2] * bonds_in_path[0] * path[1];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[2] + magic[0] * bonds_in_path[2] * path[0] +
                  magic[1] * bonds_in_path[0] * path[1] +
                  magic[2] * bonds_in_path[1] * path[2];
  MaybeUpdateHighest(tmp, highest);
  // Now counterclockwise.
  tmp = path[0] + magic[0] * bonds_in_path[2] * path[2] +
                  magic[1] * bonds_in_path[1] * path[1] +
                  magic[2] * bonds_in_path[0] * path[0];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[1] + magic[0] * bonds_in_path[0] * path[0] +
                  magic[1] * bonds_in_path[2] * path[2] +
                  magic[2] * bonds_in_path[1] * path[1];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[2] + magic[0] * bonds_in_path[1] * path[1] +
                  magic[1] * bonds_in_path[0] * path[0] +
                  magic[2] * bonds_in_path[2] * path[2];
  MaybeUpdateHighest(tmp, highest);
  return highest;
}

// A 4 membered ring is formed. There are 4 atoms and 4 bonds.
// There are 4 starting atoms and 2 directions, so 8 possible
// paths. Return the highest bit number encountered.
uint32_t
Args::CanonicalRingBit4() const {
  static uint32_t magic[] {73, 46, 93, 27};

  uint32_t highest = 0;
  // First try the clockwise directions
  uint32_t tmp = path[0] + magic[0] * bonds_in_path[0] * path[1] +
                           magic[1] * bonds_in_path[1] * path[2] +
                           magic[2] * bonds_in_path[2] * path[3] +
                           magic[3] * bonds_in_path[3] * path[0];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[1] + magic[0] * bonds_in_path[1] * path[2] +
                  magic[1] * bonds_in_path[2] * path[3] +
                  magic[2] * bonds_in_path[3] * path[0] +
                  magic[3] * bonds_in_path[0] * path[1];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[2] + magic[0] * bonds_in_path[2] * path[3] +
                  magic[1] * bonds_in_path[3] * path[0] +
                  magic[2] * bonds_in_path[0] * path[1] +
                  magic[3] * bonds_in_path[1] * path[2];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[3] + magic[0] * bonds_in_path[3] * path[0] +
                  magic[1] * bonds_in_path[0] * path[1] +
                  magic[2] * bonds_in_path[1] * path[2] +
                  magic[3] * bonds_in_path[2] * path[3];
  MaybeUpdateHighest(tmp, highest);
  // Now counterclockwise.
  tmp = path[0] + magic[0] * bonds_in_path[3] * path[3] +
                  magic[1] * bonds_in_path[2] * path[2] +
                  magic[2] * bonds_in_path[1] * path[1] +
                  magic[3] * bonds_in_path[0] * path[0];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[1] + magic[0] * bonds_in_path[0] * path[0] +
                  magic[1] * bonds_in_path[3] * path[3] +
                  magic[2] * bonds_in_path[2] * path[2] +
                  magic[3] * bonds_in_path[1] * path[1];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[2] + magic[0] * bonds_in_path[1] * path[1] +
                  magic[1] * bonds_in_path[0] * path[0] +
                  magic[2] * bonds_in_path[3] * path[3] +
                  magic[3] * bonds_in_path[2] * path[2];
  MaybeUpdateHighest(tmp, highest);
  tmp = path[3] + magic[0] * bonds_in_path[2] * path[2] +
                  magic[1] * bonds_in_path[1] * path[1] +
                  magic[2] * bonds_in_path[0] * path[0] +
                  magic[3] * bonds_in_path[3] * path[3];
  MaybeUpdateHighest(tmp, highest);
  return highest;
}

// A full length linear path is loaded in `args`.
// Generate a bit if possible, and check for possible
// implicit 4 membered rings.
void
FormBit(Molecule& m,
        Args& args,
        const TorsionOptions& options,
        Sparse_Fingerprint_Creator& destination) {
  assert(args.path_length == kMaxPathLength);
#ifdef DEBUG_TOPOLOGICAL_TORSION
  cerr << "Forming bit, began " << args.first_atom << " terminated " << args.current_atom << '\n';
  args.DebugPrint(cerr);
#endif
  // Avoid forming all bits twice.
  if (args.path[0] <= args.path[3]) {
    HitBit(args.CanonicalLinear(), destination);
  }

  if (! options.fingerprint_4_membered_rings) {
    return;
  }

  if (! args.AddRingClosureBond(m)) {  // No 4 membered ring implied.
    return;
  }

  HitBit(args.CanonicalRingBit4(), destination);
}

// Forward declaration.
void
TopologicalTorsion(Molecule& m,
                   Args& args,
                   TorsionOptions& options,
                   Sparse_Fingerprint_Creator& destination);

// Torsion formation continues with atom `args.current_atom`.
void
TopologicalTorsionIdNext(Molecule& m,
                   Args& args,
                   TorsionOptions& options,
                   Sparse_Fingerprint_Creator& destination) {
  const atom_number_t current_atom = args.current_atom;
  const Atom& a = m.atom(current_atom); 
#ifdef DEBUG_TOPOLOGICAL_TORSION
  cerr << "TopologicalTorsion prev " << args.previous_atom << " current atom " << current_atom << " path_length " << args.path_length << '\n';
#endif
  for (const Bond * b : a) {
    const atom_number_t i = b->other(current_atom);
    if (! args.include_atom[i]) {
      continue;
    }
    if (! args.visited[i]) {
    } else if (i == args.previous_atom) {
      continue;
    } else if (options.fingerprint_3_membered_rings && args.path_length == 3) {
      args.bonds_in_path[2] = NumberBonds(b);
      HitBit(args.CanonicalRingBit3(), destination);
      continue;
    } else {
      continue;
    }
    const auto old_current_atom = args.current_atom;
    const auto old_previous_atom = args.previous_atom;
    args.current_atom = i;
    args.previous_atom = current_atom;
#ifdef DEBUG_TOPOLOGICAL_TORSION
    cerr << "Setting bond " << (args.path_length - 1) << " value " << NumberBonds(b) << " from " << m.atomic_symbol(current_atom) << " to " << m.atomic_symbol(i) << '\n';
#endif
    args.bonds_in_path[args.path_length - 1] = NumberBonds(b);
    TopologicalTorsion(m, args, options, destination);
    args.current_atom = old_current_atom;
    args.previous_atom = old_previous_atom;
  }
}

// This mostly manages the state of `args`, ensuring that
// state added to it is removed after the inner call to TopologicalTorsion.
// Place `current_atom` onto the path, then either form a bit, or continue
// forming the path.
void
TopologicalTorsion(Molecule& m,
                   Args& args,
                   TorsionOptions& options,
                   Sparse_Fingerprint_Creator& destination) {
  // Push current_atom onto the path.
  args.visited[args.current_atom] = 1;
  args.path[args.path_length] = args.atom_type[args.current_atom];
  args.path_length++;

  if (args.path_length == kMaxPathLength) {
    FormBit(m, args, options, destination);
  } else {
    TopologicalTorsionIdNext(m, args, options, destination);
  }

  // Reverse earlier changes.
  args.visited[args.current_atom] = 0;
  args.path_length--;
  args.path[args.path_length] = 0;  // Not really necessary.
}

// Does the actual fingerprinting. Not public, the public version
// makes sure `include_atom` is not null.
Sparse_Fingerprint_Creator
TopologicalTorsionInner(Molecule& m, 
                        const atom_type_t* atom_type,
                        const int * include_atom,
                        TorsionOptions& options) {
  const int natoms = m.natoms();

  Sparse_Fingerprint_Creator result;

  int * visited = new int[natoms]; std::unique_ptr<int[]> free_visited(visited);

  Args args;
  args.include_atom = include_atom;
  args.atom_type = atom_type;
  args.visited = visited;

#ifdef DEBUG_TOPOLOGICAL_TORSION
  for (int i = 0; i < natoms; ++i) {
    cerr << "atom " << i << " " << m.atomic_symbol(i) << ' ' << atom_type[i] << '\n';
  }
#endif

  for (int i = 0; i < natoms; ++i) {
    if (! include_atom[i]) {
      continue;
    }
    const Atom& a = m.atom(i);
    args.path[0] = atom_type[i];
    args.first_atom = i;
    for (const Bond * b : a) {
      const atom_number_t j = b->other(i);
      if (! include_atom[j]) {
        continue;
      }
      std::fill_n(visited, natoms, 0);
      visited[i] = 1;
      args.previous_atom = i;
      args.current_atom = j;
      args.path_length = 1;
      args.bonds_in_path[0] = NumberBonds(b);
#ifdef DEBUG_TOPOLOGICAL_TORSION
      cerr << "Begin torsion at atom " << i << " " << m.atomic_symbol(i) << " to " << j << " " << m.atomic_symbol(j) << " bond " << NumberBonds(b) << '\n';
#endif
      TopologicalTorsion(m, args, options, result);
    }
  }

  return result;
}

// The public function for topological torsion fingerprinting.
// Mostly this checks sizes and if OK, ensures invokes TopologicalTorsionInner
// with a non null include_atom.
Sparse_Fingerprint_Creator
TopologicalTorsion(Molecule& m, 
                   const atom_type_t* atom_type,
                   const int * include_atom,
                   TorsionOptions& options) {
  const int natoms = m.natoms();
  if (natoms < 3) {
    return Sparse_Fingerprint_Creator();
  }
  if (natoms == 3 && ! options.fingerprint_3_membered_rings) {
    return Sparse_Fingerprint_Creator();
  }

  m.compute_aromaticity_if_needed();

  if (include_atom != nullptr) {
    return TopologicalTorsionInner(m, atom_type, include_atom, options);
  }

  std::unique_ptr<int[]> my_include(new_int(natoms, 1));

  return TopologicalTorsionInner(m, atom_type, my_include.get(), options);
}

} // namespace topological_torsion
