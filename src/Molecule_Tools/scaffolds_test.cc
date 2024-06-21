// Tester for scaffold enumeration

#include "googletest/include/gtest/gtest.h"
#include "googlemock/include/gmock/gmock.h"
#include "google/protobuf/text_format.h"

#include "scaffolds.h"

namespace {

using testing::UnorderedElementsAreArray;

struct MoleculeResult {
  IWString smiles;
  resizable_array<IWString> results;
};

std::ostream&
operator<< (std::ostream& output, const MoleculeResult& mrt) {
  output << mrt.smiles << " ->";
  for (const IWString& s : mrt.results) {
    output << ' ' << s;
  }
  return output;
}

class TestScaffolds : public testing::TestWithParam<MoleculeResult> {
  protected:
    Molecule _m;

    scaffolds::ScaffoldFinder _scaffold_finder;

    scaffolds::ScaffoldData _scaffolds;

    // To do the final test, we convert the scaffold molecules to string form.
    resizable_array<IWString> _scaffold_smiles;
};

TEST_P(TestScaffolds, TestScaffolds) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Processing " << params.smiles << " results " << params.results.size() << " items\n";
  for (const IWString& s : params.results) {
    std::cerr << "Expecting '" << s << "'\n";
  }
#endif
  _scaffold_finder.MakeScaffolds(_m, _scaffolds);

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Got " << _scaffolds.size() << " values back\n";
  std::cerr << "params.results.size() " << params.results.size() << '\n';
#endif

  EXPECT_EQ(_scaffolds.subset_size(), params.results.number_elements());

  if (_scaffolds.subset().empty()) {
    return;
  }

  for (const auto& p : _scaffolds.subset()) {
    Molecule m;
    const IWString smiles = p.smi();
    m.build_from_smiles(smiles);
    std::cerr << "Generated " << m.smiles() << " from " << p.smi() << " usmi " << m.unique_smiles() << '\n';
    _scaffold_smiles << m.smiles();
  }

#ifdef DEBUG_TEST_ADASD
  for (const auto& s : params.results) {
    std::cerr << "Expecting " << s << '\n';
  }
#endif

  EXPECT_THAT(_scaffold_smiles, UnorderedElementsAreArray(params.results)) <<
        params;
}

INSTANTIATE_TEST_SUITE_P(TestScaffolds, TestScaffolds, testing::Values(
  MoleculeResult{"C", {}},
  MoleculeResult{"C1CC1", {{"C1CC1"}}},
  MoleculeResult{"C1CC1C", {{"C1CC1"}}},
  MoleculeResult{"C12CC1C2CC1CC1", {{"C1CC1", "C(C1CC1)C1C2CC12", "C1C2CC12"}}},
  MoleculeResult{"C1CC1CC1CC1", {{"C1CC1", "C(C1CC1)C1CC1"}}},
  MoleculeResult{"C1CC1CC1CC1CC1CC1", {{"C1CC1", "C(C1CC1)C1CC1", "C1C(CC2CC2)C1CC1CC1"}}},
  MoleculeResult{"C1CC1C(C1CC1)C1CC1", {{"C1CC1", "C(C1CC1)C1CC1", "C1C(C(C2CC2)C2CC2)C1"}}},
  MoleculeResult{"O=C1CC1CC(=O)CC1CC1", {{"C1CC1"}, {"O=C1CC1"}, {"O=C(CC1CC1=O)CC1CC1"}}},
  MoleculeResult{"O=C1CC1C(C)C(=O)C(C)(C)C1CC1", {{"C1CC1"}, {"O=C1CC1"}, {"O=C(CC1CC1=O)CC1CC1"}}},
  MoleculeResult{"O=C1CC1C(C)C(=O)C(C)(C)C1C(C)C1", {{"C1CC1"}, {"O=C1CC1"}, {"O=C(CC1CC1=O)CC1CC1"}}},
  MoleculeResult{"C1CC1CC(=O)CC1CC1", {{"C1CC1"}, {"O=C(CC1CC1)CC1CC1"}}}
));

struct ProtoSmilesResult {
  std::string textproto;
  IWString smiles;
  resizable_array<IWString> results;
};

std::ostream&
operator<< (std::ostream& output, const ProtoSmilesResult& mrt) {
  output << mrt.textproto << ' ' << mrt.smiles << " ->";
  for (const IWString& s : mrt.results) {
    output << ' ' << s;
  }
  return output;
}

class TestScaffoldsProto : public testing::TestWithParam<ProtoSmilesResult> {
  protected:
    scaffolds::ScaffoldsOptions _proto;

    Molecule _m;

    scaffolds::ScaffoldFinder _scaffold_finder;

    scaffolds::ScaffoldData _scaffolds;

    // To do the final test, we convert the scaffold molecules to string form.
    resizable_array<IWString> _scaffold_smiles;
};

TEST_P(TestScaffoldsProto, TestScaffolds) {
  const auto params = GetParam();

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.textproto, &_proto));
  ASSERT_TRUE(_scaffold_finder.Initialise(_proto));
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Processing " << params.smiles << " results " << params.results.size() << " items\n";
  for (const IWString& s : params.results) {
    std::cerr << "Expecting '" << s << "'\n";
  }
#endif
  _scaffold_finder.MakeScaffolds(_m, _scaffolds);

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Got " << _scaffolds.size() << " values back\n";
  std::cerr << "params.results.size() " << params.results.size() << '\n';
#endif

  EXPECT_EQ(_scaffolds.subset_size(), params.results.number_elements());

  if (_scaffolds.subset().empty()) {
    return;
  }

  for (const auto& p : _scaffolds.subset()) {
    Molecule m;
    const IWString smiles = p.smi();
    m.build_from_smiles(smiles);
    std::cerr << "Generated " << m.smiles() << " from " << p.smi() << " usmi " << m.unique_smiles() << '\n';
    _scaffold_smiles << m.smiles();
  }

#ifdef DEBUG_TEST_ADASD
  for (const auto& s : params.results) {
    std::cerr << "Expecting " << s << '\n';
  }
#endif

  EXPECT_THAT(_scaffold_smiles, UnorderedElementsAreArray(params.results)) <<
        params;
}
INSTANTIATE_TEST_SUITE_P(TestScaffoldsProto, TestScaffoldsProto, testing::Values(
  ProtoSmilesResult{R"pb()pb", "C", {}},
  ProtoSmilesResult{
    R"pb(
      remove_ring_based_non_scaffold_atoms: false
    )pb", "CC1CC1", {"CC1CC1"}},
  ProtoSmilesResult{
    R"pb(
      remove_ring_based_non_scaffold_atoms: false
    )pb", "CC1CC1C(C)(C)C1CC1", {"CC1CC1", "C1CC1", "CC1CC1CC1CC1"}},
  ProtoSmilesResult{
    R"pb(
      remove_linker_based_non_scaffold_atoms: false
    )pb", "CC1CC1", {"C1CC1"}},
  ProtoSmilesResult{
    R"pb(
      remove_linker_based_non_scaffold_atoms: false
    )pb", "C1CC1C(C)(C)C1CC1", {"C1CC1", "CC(C)(C1CC1)C1CC1"}},
  ProtoSmilesResult{
    R"pb(
      remove_linker_based_non_scaffold_atoms: false
    )pb", "C1CC1C(CC)C1CC1", {"C1CC1", "CCC(C1CC1)C1CC1"}},
  ProtoSmilesResult{
    R"pb(
      discard_cyclopropyl_ring: true
      remove_ring_based_non_scaffold_atoms: false
    )pb", "c1ccccc1CC1CC1", {"CC(Cc1ccccc1)C"}},
  ProtoSmilesResult{
    R"pb(
      linker_isotope: 1
      substituent_isotope: 2
    )pb", "c1c(C)cccc1CC(=O)NC(CCC)C1CC1", {"[1CH2]1CC1", "[2cH]1c[1cH]ccc1", "O=C(NC[1CH]1CC1)C[1c]1c[2cH]ccc1"}},
  ProtoSmilesResult{
    R"pb(
      max_systems_in_subset: 2
    )pb", "C1CC1C1CC1C1CC1", {"C1CC1", "C1CC1C1CC1"}},
  ProtoSmilesResult{
    R"pb(
      max_length_linker: 2
    )pb", "C1CC1CC1CC1CCCCC1CC1", {"C1CC1", "C(C1CC1)C1CC1"}}
));


}  // namespace
