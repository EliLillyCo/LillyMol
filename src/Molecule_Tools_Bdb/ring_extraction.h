#ifndef MOLECULE_TOOLS_BDB_RING_EXTRACTION_H_
#define MOLECULE_TOOLS_BDB_RING_EXTRACTION_H_

#include <array>
#include <cstdint>
#include <limits>
#include <memory>

#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

namespace smi2rings_bdb {

// Key-generating settings shared by database building and lookup. Defaults
// preserve the CLI's existing behaviour; callers choose representations explicitly.
struct RingExtractionOptions {
  int _add_double_bonded_neighbors_to_rings = 0;
  int _add_single_bonded_neighbours_to_rings = 0;
  int _add_all_single_atom_ring_substituents = 0;
  int _continue_chain_neighbours_through_multiple_bonds = 0;
  int _label_attachment_points = 0;
  int _always_process_unsubstituted_ring = 0;
  int _label_ring_atoms = 0;
  int _label_attachment_points_with_environment = 0;
  int _label_attachment_points_with_existing_isotope = 0;
  int _max_ring_size_to_consider = std::numeric_limits<int>::max();
  int _discard_molecules_where_ring_is_whole_molecule = 0;
  int _max_ring_system_size_to_consider = std::numeric_limits<int>::max();
  int _convert_positive_aromatic_nitrogen_to_neutral = 0;
  int _add_connected_atom_to_positive_aromatic_nitrogen = 0;
  int _bridge_across_spiro_groups = 0;
  int _remove_invalid_chiral_centres = 0;
  int _write_3d_smiles_if_needed = 0;
  int _use_scaffolds = 0;

  void
  set_add_double_bonded_neighbors_to_rings(int value) {
    _add_double_bonded_neighbors_to_rings = value;
  }

  void
  set_add_single_bonded_neighbours_to_rings(int value) {
    _add_single_bonded_neighbours_to_rings = value;
  }

  void
  set_add_all_single_atom_ring_substituents(int value) {
    _add_all_single_atom_ring_substituents = value;
  }

  void
  set_continue_chain_neighbours_through_multiple_bonds(int value) {
    _continue_chain_neighbours_through_multiple_bonds = value;
  }

  void
  set_label_attachment_points(int value) {
    _label_attachment_points = value;
  }

  void
  set_always_process_unsubstituted_ring(int value) {
    _always_process_unsubstituted_ring = value;
  }

  void
  set_label_ring_atoms(int value) {
    _label_ring_atoms = value;
  }

  void
  set_label_attachment_points_with_environment(int value) {
    _label_attachment_points_with_environment = value;
  }

  void
  set_label_attachment_points_with_existing_isotope(int value) {
    _label_attachment_points_with_existing_isotope = value;
  }

  void
  set_max_ring_size_to_consider(int value) {
    _max_ring_size_to_consider = value;
  }

  void
  set_discard_molecules_where_ring_is_whole_molecule(int value) {
    _discard_molecules_where_ring_is_whole_molecule = value;
  }

  void
  set_max_ring_system_size_to_consider(int value) {
    _max_ring_system_size_to_consider = value;
  }

  void
  set_convert_positive_aromatic_nitrogen_to_neutral(int value) {
    _convert_positive_aromatic_nitrogen_to_neutral = value;
  }

  void
  set_add_connected_atom_to_positive_aromatic_nitrogen(int value) {
    _add_connected_atom_to_positive_aromatic_nitrogen = value;
  }

  void
  set_bridge_across_spiro_groups(int value) {
    _bridge_across_spiro_groups = value;
  }

  void
  set_remove_invalid_chiral_centres(int value) {
    _remove_invalid_chiral_centres = value;
  }

  void
  set_write_3d_smiles_if_needed(int value) {
    _write_3d_smiles_if_needed = value;
  }

  void
  set_use_scaffolds(int value) {
    _use_scaffolds = value;
  }
};

// Scratch data for one molecule. The molecule and options must outlive it.
// Not copyable: owns working arrays, reused across systems and representations.
struct RingExtractionData {
  Molecule& _m;
  const RingExtractionOptions& _options;
  int* _in_system;
  uint32_t* _atype;

  // For each ring system discovered, we need to know the number of rings in
  // it. This array will only be partially filled if there are fused rings.
  int* _ring_system_size;

  // the unique smiles of the rings identified in this molecule.
  IW_STL_Hash_Set _found_this_molecule;

  // If the molecule has isotopes, we capture them here.
  std::unique_ptr<isotope_t[]> _starting_isotopes;

  // If we are adding functional groups.
  std::unique_ptr<isotope_t[]> _functional_group;

  // We may need to adjust the smiles written for certain atom types.
  Set_of_Atoms _positive_aromatic_nitrogen;
  // Three connected only;
  Set_of_Atoms _aromatic_sulphur;

  // create_subset needs a cross reference array.
  std::unique_ptr<int[]> _xref;

 public:
  RingExtractionData(Molecule& m, const RingExtractionOptions& options);
  ~RingExtractionData();
  RingExtractionData(const RingExtractionData&) = delete;
  RingExtractionData& operator=(const RingExtractionData&) = delete;

  int
  natoms() const {
    return _m.natoms();
  }

  int AssignAtomTypes(Atom_Typing_Specification& atom_typing);

  int IdentifyThreeConnectedAromaticSulphur(Molecule& m);
  int IdentifyPositiveAromaticNitrogen(Molecule& m);
  int IdentifyFunctionalGroups(Molecule& m,
                               resizable_array_p<Substructure_Query>& queries);

  int MaybeAddOutsideRingAtoms(Molecule& m, int ring_system_number);
  int MaybeConvertPositiveAromaticNitrogen(Molecule& m, int ring_system) const;

  int AddExtraRingSingleBondedNeighgours(Molecule& m, int uid);
  int AddConnectedAtomToPositiveAromaticNitrogen(Molecule& m, int uid) const;
  int ConvertPositiveAroamticNitrogenToNeutral(Molecule& m, int uid) const;
  int AddImplicitHydrogenToAromaticSulphur(Molecule& m, int ring_system) const;

  int CreateSubset(Molecule& m, int uid, IWString& unique_smiles) const;

  // Return true if found_this_molecule contains m.unique_smiles().
  bool AlreadySeen(const IWString& unique_smiles);

  // If starting_isotopes is set, set those isotopes in `m`.
  int MaybeRevertInitialIsotopes(Molecule& m);
};

enum class RingKeyType { kContextual, kExistingIsotopes, kAttachment, kBare };

// A small iterable value, avoiding allocation when selecting key representations.
struct RingKeyTypes {
  std::array<RingKeyType, 4> _types{};
  int _size = 0;

  const RingKeyType*
  begin() const {
    return _types.data();
  }

  const RingKeyType*
  end() const {
    return _types.data() + _size;
  }
};

class RingExtractor {
 private:
  RingExtractionOptions _options;
  Atom_Typing_Specification _atom_typing;
  resizable_array_p<Substructure_Query> _queries;

 public:
  RingExtractionOptions&
  options() {
    return _options;
  }

  const RingExtractionOptions&
  options() const {
    return _options;
  }

  Atom_Typing_Specification&
  atom_typing() {
    return _atom_typing;
  }

  resizable_array_p<Substructure_Query>&
  functional_group_queries() {
    return _queries;
  }

  RingExtractor() = default;
  RingExtractor(const RingExtractor&) = delete;
  RingExtractor& operator=(const RingExtractor&) = delete;
  RingKeyTypes KeyTypes(bool building = false) const;
  int CreateKey(Molecule& m, RingExtractionData& data, int uid, RingKeyType type,
                IWString& key, bool building = false) const;

  // Operates on a caller-owned working molecule. Subsets are created as needed;
  // the full input is not copied per representation.
  int Prepare(Molecule& m, RingExtractionData& data);
  int IdentifySystems(Molecule& m, RingExtractionData& data) const;
  bool PrepareSystem(Molecule& m, RingExtractionData& data, int uid) const;
  void LabelEnvironment(Molecule& m, RingExtractionData& data, int uid) const;
  void LabelExistingIsotopes(Molecule& m, RingExtractionData& data, int uid) const;
  void LabelAttachments(Molecule& m, RingExtractionData& data, int uid) const;
};

}  // namespace smi2rings_bdb
#endif  // MOLECULE_TOOLS_BDB_RING_EXTRACTION_H_
