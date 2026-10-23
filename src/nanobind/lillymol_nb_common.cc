#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {


ReaderContext::ReaderContext(const std::string& fname, FileType file_type) {
  IWString tmp(fname);
  reader = std::make_unique<data_source_and_type<Molecule>>(file_type, tmp);
}

ReaderContext::ReaderContext(const std::string& fname) {
  IWString tmp(fname);
  FileType file_type = discern_file_type_from_name(tmp);
  if (file_type == FILE_TYPE_INVALID) {
    std::cerr << "ReaderContext::ReaderContext:unrecognised type '" << fname << "'\n";
    return;
  }

  reader = std::make_unique<data_source_and_type<Molecule>>(file_type, tmp);
}

ReaderContext::~ReaderContext() {
  ResetSdfOptions();
  RestoreElementCreationFlags();
}

void
ReaderContext::SetPreprocessing(bool reduce_to_largest_fragment, bool remove_chirality,
                                bool remove_cis_trans_bonds, bool remove_isotopes) {
  preprocessing.set_reduce_to_largest_fragment(reduce_to_largest_fragment);
  preprocessing.set_remove_chirality(remove_chirality);
  preprocessing.set_remove_cis_trans_bonds(remove_cis_trans_bonds);
  preprocessing.set_remove_isotopes(remove_isotopes);
}

void
ReaderContext::EnableNonPeriodicTableElementReading() {
  if (!element_creation_flags_saved) {
    saved_auto_create_new_elements = auto_create_new_elements();
    saved_atomic_symbols_can_have_arbitrary_length =
        atomic_symbols_can_have_arbitrary_length();
    element_creation_flags_saved = true;
  }

  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);
}

void
ReaderContext::RestoreElementCreationFlags() {
  if (!element_creation_flags_saved) {
    return;
  }

  set_auto_create_new_elements(saved_auto_create_new_elements);
  set_atomic_symbols_can_have_arbitrary_length(
      saved_atomic_symbols_can_have_arbitrary_length);
  element_creation_flags_saved = false;
}

void
ReaderContext::SetSkipOptions(bool invalid_valence, bool non_organic,
                              bool non_periodic_table_elements) {
  skip_invalid_valence = invalid_valence;
  skip_non_organic = non_organic;
  skip_non_periodic_table_elements = non_periodic_table_elements;
  if (skip_non_periodic_table_elements) {
    EnableNonPeriodicTableElementReading();
  } else {
    RestoreElementCreationFlags();
  }
}

bool
ReaderContext::ShouldSkip(Molecule& mol) {
  if (skip_invalid_valence && !mol.valence_ok()) {
    ++skipped_invalid_valence;
    return true;
  }

  if (skip_non_organic && !mol.organic_only()) {
    ++skipped_non_organic;
    return true;
  }

  if (skip_non_periodic_table_elements) {
    const int matoms = mol.natoms();
    for (int i = 0; i < matoms; ++i) {
      if (!mol.elementi(i)->is_in_periodic_table()) {
        ++skipped_non_periodic_table_elements;
        return true;
      }
    }
  }

  return false;
}

bool
ReaderContext::SetSdfOptions(const std::string& sdf_identifier, bool sdf_tags_to_json,
                             bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid,
                             bool mdlquiet) {
  if (sdf_identifier.empty() && !sdf_tags_to_json && !all_sdf_tags && !first_sdf_tag &&
      prepend_sdfid && !mdlquiet) {
    return true;
  }

  if (!reader) {
    return false;
  }

  MDL_File_Supporting_Material& mdlfos =
      reader->mutable_molecule_read_options().mdl_file_supporting_material();
  const_IWSubstring tmp(sdf_identifier.data(), sdf_identifier.length());
  if (!mdlfos.set_sdf_identifier(tmp)) {
    return false;
  }

  mdlfos.set_sdf_tags_to_json(sdf_tags_to_json);
  mdlfos.set_fetch_all_sdf_identifiers(all_sdf_tags);
  mdlfos.set_take_first_tag_as_name(first_sdf_tag);
  mdlfos.set_prepend_sdfid(prepend_sdfid);
  if (mdlquiet) {
    mdlfos.set_report_unrecognised_records(false);
  }

  return true;
}

bool
ReaderContext::ApplyOptions(const ReaderOptions& options) {
  SetPreprocessing(options.largest_fragment, options.remove_chirality,
                   options.remove_cis_trans_bonds, options.remove_isotopes);
  SetSkipOptions(options.skip_invalid_valence, options.skip_non_organic,
                 options.skip_non_periodic_table_elements);
  if (reader) {
    reader->mutable_molecule_read_options()
        .mdl_file_supporting_material()
        .set_read_extra_text_info(options.keep_sdf_tags);
  }
  return SetSdfOptions(options.sdf_identifier, options.sdf_tags_to_json,
                       options.all_sdf_tags, options.first_sdf_tag,
                       options.prepend_sdfid, options.mdlquiet);
}

void
ReaderContext::ResetSdfOptions() {
  if (reader) {
    reader->reset_molecule_read_options();
  }
}

std::optional<Molecule>
ReaderContext::Next() {
  if (!reader) {
    return std::nullopt;
  }

  while (true) {
    std::unique_ptr<Molecule> mol(reader->next_molecule());
    if (!mol) {
      return std::nullopt;
    }

    if (preprocessing.active()) {
      preprocessing.Process(*mol);
    }

    if (ShouldSkip(*mol)) {
      continue;
    }

    return std::move(*mol);
  }
}

ContextWriter::ContextWriter(const std::string& fname) {
  CommonFileOpen(fname, FILE_TYPE_SMI);
}

ContextWriter::ContextWriter(const std::string& fname, FileType file_type) {
  CommonFileOpen(fname, file_type);
}

void
ContextWriter::CommonFileOpen(const std::string& fname, FileType file_type) {
  writer = std::make_unique<Molecule_Output_Object>();
  if (!writer->add_output_type(file_type)) {
    std::cerr << "MolWriterContext::CommonFileOpen:unrecognised type " << file_type << '\n';
    writer.reset(nullptr);
    return;
  }

  IWString tmp(fname);
  if (!writer->new_stem(tmp)) {
    std::cerr << "MolWriterContext::CommonFileOpen:cannot open '" << tmp << "'\n";
    writer.reset(nullptr);
    return;
  }
}

std::optional<Molecule>
ReaderNext(data_source_and_type<Molecule>& reader) {
  std::unique_ptr<Molecule> mol(reader.next_molecule());
  if (!mol) {
    return std::nullopt;
  }
  return std::move(*mol);
}

BondType ToBondType(const Bond& bond) {
  if (bond.is_aromatic()) {
    return BondType::kAromaticBond;
  }
  if (bond.is_single_bond()) {
    return BondType::kSingleBond;
  }
  if (bond.is_double_bond()) {
    return BondType::kDoubleBond;
  }
  if (bond.is_triple_bond()) {
    return BondType::kTripleBond;
  }
  throw std::invalid_argument("Unrecognised bond type");
}

bond_type_t ToLillyMolBondType(BondType bond_type) {
  switch (bond_type) {
    case BondType::kSingleBond:
      return SINGLE_BOND;
    case BondType::kDoubleBond:
      return DOUBLE_BOND;
    case BondType::kTripleBond:
      return TRIPLE_BOND;
    default:
      throw std::invalid_argument("Unrecognised bond type");
  }
}

void
AppendChiralComponent(int atom, IWString& result) {
  if (atom >= 0) {
    result << atom;
  } else if (atom == kChiralConnectionIsImplicitHydrogen) {
    result << 'H';
  } else if (atom == kChiralConnectionIsLonePair) {
    result << '^';
  } else {
    result << '?';
  }
}

std::string
ChiralCentreRepr(const Chiral_Centre& centre) {
  IWString result;
  result << "<Chiral_Centre atom " << centre.a();
  result << " tf ";
  AppendChiralComponent(centre.top_front(), result);
  result << " tb ";
  AppendChiralComponent(centre.top_back(), result);
  result << " ld ";
  AppendChiralComponent(centre.left_down(), result);
  result << " rd ";
  AppendChiralComponent(centre.right_down(), result);
  result << '>';
  return result.AsString();
}

std::optional<int>
IndexOf(const int* values, int needle) {
  for (int i = 0; i < 4; ++i) {
    if (values[i] == needle) {
      return i;
    }
  }

  return std::nullopt;
}

int
InversionCount(const std::vector<int>& values) {
  int result = 0;
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(values.size()); ++j) {
      if (values[i] > values[j]) {
        ++result;
      }
    }
  }

  return result;
}

int
ChiralConnectionSortKey(int member, int matoms) {
  if (member >= 0) {
    return member;
  }
  if (member == kChiralConnectionIsImplicitHydrogen) {
    return matoms;
  }
  if (member == kChiralConnectionIsLonePair) {
    return matoms + 1;
  }

  return matoms + 2;
}

std::optional<ChiralType>
TetrahedralChirality(Molecule& mol, atom_number_t atom, bool check_is_chiral) {
  if (atom < 0 || atom >= mol.natoms()) {
    throw std::invalid_argument("atom number outside molecule");
  }

  if (check_is_chiral && !::is_actually_chiral(mol, atom)) {
    return std::nullopt;
  }

  const Chiral_Centre* centre = mol.chiral_centre_at_atom(atom);
  if (centre == nullptr) {
    return std::nullopt;
  }

  if (!centre->chirality_known() || !centre->complete()) {
    return ChiralType::kChiUnspecified;
  }

  if (centre->lone_pair_count() > 0) {
    return ChiralType::kChiOther;
  }

  const int stored[4] = {
      centre->top_front(),
      centre->top_back(),
      centre->left_down(),
      centre->right_down()
  };

  std::vector<int> ordered_members(std::begin(stored), std::end(stored));
  const int matoms = mol.natoms();
  std::sort(ordered_members.begin(), ordered_members.end(),
            [matoms](int lhs, int rhs) {
              return ChiralConnectionSortKey(lhs, matoms) <
                     ChiralConnectionSortKey(rhs, matoms);
            });

  std::vector<int> permutation;
  permutation.reserve(4);
  for (int member : ordered_members) {
    std::optional<int> index = IndexOf(stored, member);
    if (!index) {
      return ChiralType::kChiUnspecified;
    }
    permutation.push_back(*index);
  }

  if (InversionCount(permutation) % 2 == 0) {
    return ChiralType::kChiTetrahedralCcw;
  }

  return ChiralType::kChiTetrahedralCw;
}


std::unique_ptr<Substructure_Query> BuildQueryFromSmarts(
    const std::string& smarts, const SmartsSearchOptions& options) {
  auto query = std::make_unique<Substructure_Query>();
  if (!query->CreateFromSmarts(smarts)) {
    throw std::invalid_argument("Invalid smarts");
  }

  if (options.max_matches_to_find) {
    query->set_max_matches_to_find(*options.max_matches_to_find);
  }
  if (options.unique_embeddings_only) {
    query->set_find_unique_embeddings_only(*options.unique_embeddings_only);
  }
  if (options.one_embedding_per_start_atom) {
    query->set_find_one_embedding_per_atom(*options.one_embedding_per_start_atom);
  }
  if (options.perceive_symmetry_equivalent_matches) {
    query->set_perceive_symmetry_equivalent_matches(
        *options.perceive_symmetry_equivalent_matches);
  }

  return query;
}

std::vector<int> EmbeddingAsVector(const Set_of_Atoms& embedding) {
  std::vector<int> result;
  result.reserve(embedding.number_elements());
  for (atom_number_t atom : embedding) {
    result.push_back(atom);
  }
  return result;
}

std::vector<std::vector<int>> SubstructureResultsAsVectors(
    const Substructure_Results& sresults) {
  std::vector<std::vector<int>> result;
  result.reserve(sresults.number_embeddings());
  for (const Set_of_Atoms* embedding : sresults.embeddings()) {
    result.push_back(EmbeddingAsVector(*embedding));
  }
  return result;
}

std::vector<Set_of_Atoms> SubstructureResultsAsSetOfAtoms(
    const Substructure_Results& sresults) {
  std::vector<Set_of_Atoms> result;
  result.reserve(sresults.number_embeddings());
  for (const Set_of_Atoms* embedding : sresults.embeddings()) {
    result.push_back(Set_of_Atoms(*embedding));
  }
  return result;
}

std::vector<std::vector<int>> SubstructureSearchMatches(
    Substructure_Query& query, Molecule& molecule) {
  Substructure_Results query_results;
  uint32_t matches;
  {
    nb::gil_scoped_release release;
    matches = query.substructure_search(&molecule, query_results);
  }
  if (matches == 0) {
    return std::vector<std::vector<int>>();
  }
  return SubstructureResultsAsVectors(query_results);
}

SmartsSearchOptions MakeSmartsSearchOptions(
    std::optional<int> max_matches_to_find,
    std::optional<bool> unique_embeddings_only,
    std::optional<bool> one_embedding_per_start_atom,
    std::optional<bool> perceive_symmetry_equivalent_matches) {
  SmartsSearchOptions options;
  options.max_matches_to_find = max_matches_to_find;
  options.unique_embeddings_only = unique_embeddings_only;
  options.one_embedding_per_start_atom = one_embedding_per_start_atom;
  options.perceive_symmetry_equivalent_matches = perceive_symmetry_equivalent_matches;
  return options;
}


bool MoleculeContainsAtomicNumber(const Molecule& mol, atomic_number_t atomic_number) {
  return mol.natoms(atomic_number) > 0;
}

bool MoleculeContainsElementSymbol(const Molecule& mol, const std::string& symbol) {
  const_IWSubstring tmp(symbol);
  const Element* element = get_element_from_symbol_no_case_conversion(tmp);
  if (element == nullptr) {
    throw std::invalid_argument("Unrecognised element type");
  }
  return mol.natoms(element) > 0;
}

bool MoleculeContainsSubstructureQuery(Molecule& mol, Substructure_Query& query) {
  nb::gil_scoped_release release;
  return query.substructure_search(&mol) > 0;
}

std::string AtomRepr(const Atom& atom) {
  IWString result;
  result << "<Atom " << atom.atomic_symbol() << " ncon " << atom.ncon();
  if (atom.isotope()) {
    result << " iso " << atom.isotope();
  }
  result << '>';
  return result.AsString();
}

std::string BondRepr(const Bond& bond) {
  IWString result;
  result << "<Bond " << bond.a1();
  if (bond.is_single_bond()) {
    result << '-';
  } else if (bond.is_double_bond()) {
    result << '=';
  } else if (bond.is_triple_bond()) {
    result << '#';
  }
  result << bond.a2() << '>';
  return result.AsString();
}

exact_mass_t AtomExactMass(const Atom& atom) {
  exact_mass_t result;
  atom.exact_mass(result);
  return result;
}

int BondNrings(const Bond& bond) {
  int nrings;
  if (bond.nrings(nrings)) {
    return nrings;
  }
  throw std::runtime_error("Bond.nrings:ring membership not computed");
}

std::vector<const Bond*> Bonds(const Molecule& mol) {
  std::vector<const Bond*> result;
  result.reserve(mol.nedges());
  for (const Bond* bond : mol.bond_list()) {
    result.push_back(bond);
  }
  return result;
}

std::string SetOfAtomsRepr(const Set_of_Atoms& atoms) {
  IWString result;
  result << "<Set_of_Atoms";
  for (atom_number_t atom : atoms) {
    result << ' ' << atom;
  }
  result << '>';
  return result.AsString();
}

std::string SetOfAtomsStr(const Set_of_Atoms& atoms) {
  IWString result;
  for (atom_number_t atom : atoms) {
    result << ' ' << atom;
  }
  return result.AsString();
}

std::vector<int> GetRingMembership(Molecule& mol) {
  std::vector<int> result(mol.natoms());
  mol.ring_membership(result.data());
  return result;
}

std::vector<int> RingBondCount(Molecule& mol) {
  const int* ring_bond_count = mol.ring_bond_count();
  return std::vector<int>(ring_bond_count, ring_bond_count + mol.natoms());
}

std::vector<std::unique_ptr<Ring>> Rings(Molecule& mol) {
  std::vector<std::unique_ptr<Ring>> result;
  for (const Ring* ring : mol.sssr_rings()) {
    result.push_back(std::make_unique<Ring>(*ring));
  }
  return result;
}

std::unique_ptr<Molecule> MolFromSmiles(const std::string& smiles) {
  auto result = std::make_unique<Molecule>();
  if (!result->build_from_smiles(smiles)) {
    return nullptr;
  }
  return result;
}

std::vector<int> AtomicNumbers(const Molecule& mol) {
  std::vector<int> result;
  const int matoms = mol.natoms();
  result.reserve(matoms);
  for (int i = 0; i < matoms; ++i) {
    result.push_back(mol.atomic_number(i));
  }
  return result;
}

std::vector<int> Isotopes(const Molecule& mol) {
  std::vector<int> result;
  const int matoms = mol.natoms();
  result.reserve(matoms);
  for (int i = 0; i < matoms; ++i) {
    result.push_back(mol.isotope(i));
  }
  return result;
}

std::tuple<int, int> LipinskiHbaHbd(Molecule& mol) {
  int hba = 0;
  int hbd = 0;
  mol.LipinskiHbaHbd(hba, hbd);
  return std::make_tuple(hba, hbd);
}

std::optional<float> ALogP(Molecule& mol) {
  alogp::ALogP calc;
  calc.set_use_alcohol_for_acid(1);
  calc.set_rdkit_phoshoric_acid_hydrogen(1);
  std::optional<double> result = calc.LogP(mol);
  if (!result) {
    return std::nullopt;
  }
  return static_cast<float>(*result);
}

std::optional<float> XLogPValue(Molecule& mol) {
  xlogp::XLogPCalc calc;
  std::optional<double> result = calc.LogP(mol);
  if (!result) {
    return std::nullopt;
  }
  return static_cast<float>(*result);
}


}  // namespace lillymol_nb
