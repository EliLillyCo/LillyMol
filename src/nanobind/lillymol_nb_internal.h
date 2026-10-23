#ifndef LILLYMOL_NANOBIND_LILLYMOL_NB_INTERNAL_H_
#define LILLYMOL_NANOBIND_LILLYMOL_NB_INTERNAL_H_

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/hybridization.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/xlogp.h"
#include "nanobind/tsubstructure.h"

namespace nb = nanobind;

namespace lillymol_nb {

using molecule_processing::MoleculePreprocessing;

enum BondType {
  kUnknown = 0,
  kSingleBond = 1,
  kDoubleBond = 2,
  kTripleBond = 3,
  kAromaticBond = 4,
};

enum ChiralType {
  kChiUnspecified = 0,
  kChiTetrahedralCw = 1,
  kChiTetrahedralCcw = 2,
  kChiOther = 3,
};

struct SmartsSearchOptions {
  std::optional<int> max_matches_to_find;
  std::optional<bool> unique_embeddings_only;
  std::optional<bool> one_embedding_per_start_atom;
  std::optional<bool> perceive_symmetry_equivalent_matches;
};

struct RingInfo {
  Molecule* mol = nullptr;

  explicit RingInfo(Molecule* m) : mol(m) {}
};

struct ReaderOptions {
  bool largest_fragment = false;
  bool remove_chirality = false;
  bool remove_cis_trans_bonds = false;
  bool remove_isotopes = false;
  bool keep_sdf_tags = false;

  std::string sdf_identifier;
  bool sdf_tags_to_json = false;
  bool all_sdf_tags = false;
  bool first_sdf_tag = false;
  bool prepend_sdfid = true;

  bool mdlquiet = false;

  bool skip_invalid_valence = false;
  bool skip_non_organic = false;
  bool skip_non_periodic_table_elements = false;
};

struct ReaderContext {
  std::unique_ptr<data_source_and_type<Molecule>> reader;
  MoleculePreprocessing preprocessing;

  bool skip_invalid_valence = false;
  bool skip_non_organic = false;
  bool skip_non_periodic_table_elements = false;
  uint64_t skipped_invalid_valence = 0;
  uint64_t skipped_non_organic = 0;
  uint64_t skipped_non_periodic_table_elements = 0;

  bool element_creation_flags_saved = false;
  int saved_auto_create_new_elements = 0;
  int saved_atomic_symbols_can_have_arbitrary_length = 0;

  ReaderContext(const std::string& fname, FileType file_type);
  explicit ReaderContext(const std::string& fname);
  ~ReaderContext();

  bool ApplyOptions(const ReaderOptions& options);
  void SetPreprocessing(bool reduce_to_largest_fragment, bool remove_chirality,
                        bool remove_cis_trans_bonds, bool remove_isotopes);
  void SetSkipOptions(bool invalid_valence, bool non_organic,
                      bool non_periodic_table_elements);
  bool ShouldSkip(Molecule& mol);
  void EnableNonPeriodicTableElementReading();
  void RestoreElementCreationFlags();
  bool SetSdfOptions(const std::string& sdf_identifier, bool sdf_tags_to_json,
                     bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid,
                     bool mdlquiet);
  void ResetSdfOptions();
  std::optional<Molecule> Next();
};

struct ContextWriter {
  std::unique_ptr<Molecule_Output_Object> writer;

  explicit ContextWriter(const std::string& fname);
  ContextWriter(const std::string& fname, FileType file_type);

 private:
  void CommonFileOpen(const std::string& fname, FileType file_type);
};

BondType ToBondType(const Bond& bond);
bond_type_t ToLillyMolBondType(BondType bond_type);
void AppendChiralComponent(int atom, IWString& result);
std::string ChiralCentreRepr(const Chiral_Centre& centre);
std::optional<ChiralType> TetrahedralChirality(
    Molecule& mol, atom_number_t atom, bool check_is_chiral);
std::string AtomRepr(const Atom& atom);
std::string BondRepr(const Bond& bond);
exact_mass_t AtomExactMass(const Atom& atom);
int BondNrings(const Bond& bond);
std::vector<const Bond*> Bonds(const Molecule& mol);
std::string SetOfAtomsRepr(const Set_of_Atoms& atoms);
std::string SetOfAtomsStr(const Set_of_Atoms& atoms);
std::vector<int> GetRingMembership(Molecule& mol);
std::vector<int> RingBondCount(Molecule& mol);
std::vector<std::unique_ptr<Ring>> Rings(Molecule& mol);
std::vector<int> AtomicNumbers(const Molecule& mol);
std::vector<int> Isotopes(const Molecule& mol);
bool MoleculeContainsAtomicNumber(const Molecule& mol, atomic_number_t atomic_number);
bool MoleculeContainsElementSymbol(const Molecule& mol, const std::string& symbol);
bool MoleculeContainsSubstructureQuery(Molecule& mol, Substructure_Query& query);
std::optional<Molecule> ReaderNext(data_source_and_type<Molecule>& reader);
std::unique_ptr<Molecule> MolFromSmiles(const std::string& smiles);
std::vector<int> EmbeddingAsVector(const Set_of_Atoms& embedding);
std::vector<std::vector<int>> SubstructureResultsAsVectors(
    const Substructure_Results& sresults);
std::vector<Set_of_Atoms> SubstructureResultsAsSetOfAtoms(
    const Substructure_Results& sresults);
std::vector<std::vector<int>> SubstructureSearchMatches(
    Substructure_Query& query, Molecule& molecule);
std::unique_ptr<Substructure_Query> BuildQueryFromSmarts(
    const std::string& smarts, const SmartsSearchOptions& options);
SmartsSearchOptions MakeSmartsSearchOptions(
    std::optional<int> max_matches_to_find,
    std::optional<bool> unique_embeddings_only,
    std::optional<bool> one_embedding_per_start_atom,
    std::optional<bool> perceive_symmetry_equivalent_matches);
std::tuple<int, int> LipinskiHbaHbd(Molecule& mol);
std::optional<float> ALogP(Molecule& mol);
std::optional<float> XLogPValue(Molecule& mol);

void BindIo(nb::module_& m);
void BindAtomBond(nb::module_& m);
void BindChirality(nb::module_& m);
void BindSetOfAtomsAndRing(nb::module_& m);
void BindSubstructure(nb::module_& m);
void BindTSubstructure(nb::module_& m);
void BindMolecule(nb::module_& m);
void BindDescriptors(nb::module_& m);
void BindStandardise(nb::module_& m);
void BindFingerprint(nb::module_& m);
void BindReaction(nb::module_& m);
void BindTools(nb::module_& m);

}  // namespace lillymol_nb

#endif  // LILLYMOL_NANOBIND_LILLYMOL_NB_INTERNAL_H_
