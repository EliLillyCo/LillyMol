#include "Molecule_Tools_Bdb/substituent_identification_lib.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "gtest/gtest.h"


namespace {

std::string MakeTempDir() {
  std::string dirname = "/tmp/substituent_identification_lib_test.XXXXXX";
  char* tmpdir = mkdtemp(dirname.data());
  if (tmpdir == nullptr) {
    return "";
  }

  return dirname;
}

void WriteTextFile(const std::string& fname, const std::string& contents) {
  std::ofstream output(fname);
  output << contents;
}

}  // namespace

TEST(SubstituentIdentificationLookup, QuerySetupAndEmptyLookup) {
  SubstituentIdentificationLookup lookup;

  EXPECT_TRUE(lookup.AddQueryFromSmarts("[CD4]"));
  EXPECT_FALSE(lookup.AddQueryFromSmarts("["));

  lookup.set_min_shell_radius(0);
  lookup.set_break_molecule_at_first_two_matched_atoms(0);
  lookup.set_matched_atoms_to_process(1);
  lookup.set_min_substituent_size(1);
  lookup.set_max_substituent_size(4);
  lookup.set_max_atoms_in_product(20);
  lookup.set_min_examples_needed_for_addition(1);
  lookup.set_max_molecules_per_input_molecule(10);
  lookup.set_remove_isotopes_from_product(1);
  lookup.set_write_fragments_added(1);

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("C methane"));

  std::vector<SubstituentReplacement> replacements;
  EXPECT_TRUE(lookup.GenerateReplacements(mol, replacements));
  EXPECT_TRUE(replacements.empty());
  EXPECT_TRUE(lookup.GenerateReplacements(mol).empty());
}

TEST(SubstituentIdentificationLookup, DefaultStartingPointsWithoutDatabase) {
  SubstituentIdentificationLookup lookup;
  lookup.set_default_new_molecule_starting_points(1);

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CC ethane"));

  std::vector<SubstituentReplacement> replacements;
  EXPECT_TRUE(lookup.GenerateReplacements(mol, replacements));
  EXPECT_TRUE(replacements.empty());
}

TEST(SubstituentIdentificationLookup, ReadsCliBuiltDatabase) {
  const std::string tmpdir = MakeTempDir();
  ASSERT_FALSE(tmpdir.empty());

  const std::string build_fname = tmpdir + "/build.smi";
  const std::string dbname = tmpdir + "/subs.bdb";
  WriteTextFile(build_fname, "CC ethane\nCO methanol\nCN methylamine\n");

  std::vector<std::string> args = {
      "substituent_identification", "-d", dbname, "-B", "-R", "2", "-w", "1",
      "-M", "2", build_fname};
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (std::string& arg : args) {
    argv.push_back(arg.data());
  }

  ASSERT_EQ(SubstituentIdentificationMain(static_cast<int>(argv.size()), argv.data()), 0);

  SubstituentIdentificationLookup lookup;
  ASSERT_TRUE(lookup.AddDatabase(dbname));
  lookup.set_default_new_molecule_starting_points(1);
  lookup.set_max_substituent_size(2);

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("C methane"));

  std::vector<SubstituentReplacement> replacements = lookup.GenerateReplacements(mol);
  EXPECT_EQ(replacements.size(), 3u);

  std::vector<std::string> donors;
  for (const SubstituentReplacement& replacement : replacements) {
    donors.push_back(replacement.donor);
    EXPECT_EQ(replacement.name, "methane");
    EXPECT_EQ(replacement.radius, 1u);
    EXPECT_GE(replacement.examples, 1u);
    EXPECT_FALSE(replacement.smiles.empty());
    EXPECT_EQ(replacement.molecule.name(), "methane");
  }
  std::sort(donors.begin(), donors.end());
  EXPECT_EQ(donors, (std::vector<std::string>{"ethane", "methanol", "methylamine"}));
}
