#ifndef MOLECULE_TOOLS_SCAFFOLDS_H
#define MOLECULE_TOOLS_SCAFFOLDS_H

// Form all possible scaffold subsets for a given molecule.
// One of the criticisms of the Murcko scaffold is that the sequence
// methyl->ethyl->propyl->clclopropyl results in a different scaffold.
// This tool supports the idea of enumerating all possible scaffold
// subsets.

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds {

// Avoid passing around a bunch of separate arguments.
class PerMoleculeData {
  private:
    // Number of atoms in the molecule.
    int _matoms;

    // As returned by m.label_atoms_by_ring_system
    int* _ring_sys;

    // An indicator of whether or not an atom is in the current system
    // being assembled.
    int* _in_system;

    // The number of ring systems (returned by mlabel_atoms_by_ring_system).
    int _nsys;

    // For each non ring atom, a label which will be the same for all atoms
    // that are in the same inter-ring region.
    int* _region;

    // For linker regions, we keep track of the number of atoms in the region.
    int* _atoms_in_region;

    // When we assign ring_systems, we need to keep track of the region between
    // each pair of ring systems.
    int* _between;

    // When building the subsets, we need an array for that.
    int* _atoms_in_subset;

    // Avoid duplicates
    IW_STL_Hash_Set _seen;

    // The path tracing functions need a temporary array to keep track
    // of which atoms have been visited.
    int* _tmp;

  // private functions
    void SetBetween(int r1, int r2, int region_id);
    int ClassifyRegions(Molecule& m, const int* spinach, int classify_spinach);
    int PropagateRegion(Molecule& m, atom_number_t zatom, int flag);
    int ExtendRingSystem(const Molecule& m, atom_number_t zatom, int flag);
    int AssignLinkerRegion(const Molecule& m, atom_number_t previous,
                atom_number_t zatom, const int *spinach, int flag);
    int ApplySubstituentIsotope(Molecule& m, isotope_t substituent) const;
    int ApplyLinkerIsotope(Molecule& m, isotope_t substituent) const;

  public:
    // if `classify_spinach` is set, then we also identify how the spinach
    // regions are attached to the scaffold.
    PerMoleculeData(Molecule& m, int classify_spinach);
    ~PerMoleculeData();

    int nsys() const {
      return _nsys;
    }

    const int* ring_sys() const {
      return _ring_sys;
    }
    int* ring_sys() {
      return _ring_sys;
    }
    int ring_sys(int atom) const {
      return _ring_sys[atom];
    }
    const int* in_system() const {
      return _in_system;
    }

    int* in_system() {
      return _in_system;
    }
    int in_system(int atom) const {
      return _in_system[atom];
    }

    int* region() {
      return _region;
    }
    int region(int atom) const {
      return _region[atom];
    }

    int* atoms_in_subset() {
      return _atoms_in_subset;
    }

    // For ring atoms that are in the scaffold, and which have a non-ring
    // connection outside the ring, apply either `linker` or `substituent`
    // to the atom, depending on whether what is attached is an inter-ring
    // region or a ring substituent.
    int ApplyIsotopicLabels(Molecule& m, isotope_t linker,
                            isotope_t substituent) const;

    // Return true if the combination of ring systems in `state` is disconnected.
    int StateIsDisconnected(const std::vector<int>& state) const;

    // For every atom that is between ring systems `r1` and `r2`, set the
    // value in `atoms_in_subset` to `flag`.
    int AddInterRingAtoms(int r1, int r2, int* atoms_in_subset, int flag) const;

    // Return 1 if this molecule has been seen before.
    int Seen(Molecule& m);

    int* tmp() {
      return _tmp;
    }
};

class ScaffoldFinder {

  // All settable configuration options in the proto.
  scaffolds::ScaffoldsOptions _config;

  // Keep track of the number of ring systems in the incoming molecules.
  extending_resizable_array<int> _nsys;
  // The number of scaffolds generated.
  extending_resizable_array<int> _generated;

  // We can specify the maximum number of ring systems lost in any particular scaffold
  // combination
  int _max_systems_lost;
  // We can also specify a minumim number of ring systems in a combination.
  int _min_systems_in_subset;

   private:

  // private functions.

    int OkSubsetSize(uint32_t in_molecule, uint32_t in_subset) const;

    int Process(Molecule& m,
                 scaffolds::ScaffoldData& results);
    int Process(const Molecule& m,
                 const std::vector<int>& state,
                 PerMoleculeData& pmd,
                 scaffolds::ScaffoldData& results);

    int AddToResultsIfUnique(Molecule& candidate,
                 PerMoleculeData& pmd,
                 int nsys,
                 scaffolds::ScaffoldData& result);
    int OkLinkerSize(uint32_t atoms_in_linker) const;
    int MaybeApplyIsotopicLabels(Molecule& m, const int* spinach);
    int ApplyIsotopicLabels(Molecule& m, isotope_t linker,
                isotope_t substituent, const int* spinach) const;
    int ApplySubstituentIsotope(Molecule& m, isotope_t iso, const int* spinach) const;
    int ApplyLinkerIsotope(Molecule& m, isotope_t iso, const int* spinach) const;

  public:
    ScaffoldFinder();

    int Initialise(Command_Line& cl, char flag);
    int Initialise(const scaffolds::ScaffoldsOptions& proto);

    // Place in `results` all scaffold subsets of `m`.
    int MakeScaffolds(Molecule& m, scaffolds::ScaffoldData& results);

    int Report(std::ostream& output) const;
};

}  // namespace scaffolds

#endif // MOLECULE_TOOLS_SCAFFOLDS_H
