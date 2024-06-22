#ifndef MOLECULE_TOOLS_BDB_SMI2RINGS_BDB_LIB_H_
#define MOLECULE_TOOLS_BDB_SMI2RINGS_BDB_LIB_H_

#include <limits>
#include <memory>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"

namespace mol2rings {

//  For each ring system, a measure of uniqueness.

class Ring_System_Uniqueness {
 private:
  int _uid;
  int _uniqueness_measure;

  //  If we are looking up in a database with names, we can return
  //  the name of the molecule that first saw that ring in the database

  IWString _name_from_database;

  // If we are doing lookups across multiple databases, we keep track
  // of the number of databases in which a ring is found.
  int _ndb_found;

 public:
  Ring_System_Uniqueness();

  void set_uid(int s) {
    _uid = s;
  }

  int uid() const {
    return _uid;
  }

  void set_uniqueness_measure(int s) {
    _uniqueness_measure = s;
  }

  int uniqueness_measure() const {
    return _uniqueness_measure;
  }

  void set_name_from_database(const const_IWSubstring& s) {
    _name_from_database = s;
  }

  const IWString& name_from_database() const {
    return _name_from_database;
  }

  void set_ndb_found(int s) {
    _ndb_found = s;
  }

  int ndb_found() const {
    return _ndb_found;
  }
};


// A datastructure holding information about the molecule currently being
// processed. More about avoiding passing around long argument lists.
struct PerMoleculeData {
  Molecule& m;
  int* in_system;
  int* atype;

  // For each ring system discovered, we need to know the number of rings in
  // it. This array will only be partially filled if there are fused rings.
  int* ring_system_size;

  // For each ring system found, a measure of uniqueness.
  Ring_System_Uniqueness* rsu;

  // the unique smiles of the rings identified in this molecule.
  IW_STL_Hash_Set found_this_molecule;

  // If the molecule has isotopes, we capture them here.
  std::unique_ptr<isotope_t[]> starting_isotopes;

  // If we are adding functional groups.
  std::unique_ptr<isotope_t[]> functional_group;

  // We may need to adjust the smiles written for certain atom types.
  Set_of_Atoms positive_aromatic_nitrogen;
  // Three connected only;
  Set_of_Atoms aromatic_sulphur;

  // create_subset needs a cross reference array.
  std::unique_ptr<int[]> xref;

 public:
  PerMoleculeData(Molecule& m);
  ~PerMoleculeData();

  int natoms() const {
    return m.natoms();
  }

  int AssignAtomTypes(Atom_Typing_Specification& atom_typing);

  int IdentifyThreeConnectedAromaticSulphur(Molecule& m);
  int IdentifyPositiveAromaticNitrogen(Molecule& m);
  int IdentifyFunctionalGroups(Molecule& m);

  int MaybeAddOutsideRingAtoms(Molecule& m,
                         int ring_system_number);
  int MaybeConvertPositiveAromaticNitrogen(Molecule& m,
                int ring_system) const;

  int AddExtraRingSingleBondedNeighgours(Molecule& m, int uid);
  int AddConnectedAtomToPositiveAromaticNitrogen(Molecule& m, int uid) const;
  int ConvertPositiveAroamticNitrogenToNeutral(Molecule& m, int uid) const;
  int AddImplicitHydrogenToAromaticSulphur(Molecule& m, int ring_system) const;

  int CreateSubset(Molecule& m, int uid,
                IWString& unique_smiles) const;

  // Return true if found_this_molecule contains m.unique_smiles().
  bool AlreadySeen(const IWString& unique_smiles);

  // If starting_isotopes is set, set those isotopes in `m`.
  int MaybeRevertInitialIsotopes(Molecule& m);
};

// A class that can do lookups and stores with rings databases.
class Smi2Rings {
  private:
    int _verbose = 0;

    int label_attachment_points = 0;

    int add_double_bonded_neighbors_to_rings = 0;

    int add_single_bonded_neighbours_to_rings = 0;
    // This will add CH3, OH, NH2, etc groups
    int add_all_single_atom_ring_substituents = 0;

    int continue_chain_neighbours_through_multiple_bonds = 0;
    // If we are processing isotopically labelled attachment points, we can also
    // process the unsubstuted ring.
    int always_process_unsubstituted_ring = 0;

    int label_ring_atoms = 0;

    int label_attachment_points_with_environment = 0;

    int max_ring_size_to_consider = std::numeric_limits<int>::max();

    Atom_Typing_Specification _atom_typing;

  // private functions

  public:
    Smi2Rings();

    int Initialise(Command_Line& cl);
};

}  // namespace mol2rings

#endif // MOLECULE_TOOLS_BDB_SMI2RINGS_BDB_LIB_H_
