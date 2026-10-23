#ifndef MOLECULE_TOOLS_BDB_SUBSTITUENT_IDENTIFICATION_LIB_H_
#define MOLECULE_TOOLS_BDB_SUBSTITUENT_IDENTIFICATION_LIB_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "Molecule_Lib/molecule.h"

struct SubstituentReplacement {
  Molecule molecule;
  std::string smiles;
  std::string name;
  std::string donor;
  uint32_t radius = 0;
  uint32_t examples = 0;
  std::string fragment_lost;
  std::string fragment_added;
};

class SubstituentIdentification;

class SubstituentIdentificationLookup {
 public:
  SubstituentIdentificationLookup();
  ~SubstituentIdentificationLookup();

  SubstituentIdentificationLookup(const SubstituentIdentificationLookup&) = delete;
  SubstituentIdentificationLookup& operator=(const SubstituentIdentificationLookup&) = delete;

  bool AddDatabase(const std::string& dbname);
  void close();
  bool AddQueryFromSmarts(const std::string& smarts);
  void set_default_new_molecule_starting_points(int s);
  void set_min_shell_radius(int s);
  void set_only_produce_molecules_at_biggest_radius(int s);
  void set_break_molecule_at_first_two_matched_atoms(int s);
  void set_matched_atoms_to_process(int s);
  void set_min_substituent_size(int s);
  void set_max_substituent_size(int s);
  void set_max_atoms_in_product(int s);
  void set_min_examples_needed_for_addition(uint32_t s);
  void set_max_molecules_per_input_molecule(uint32_t s);
  void set_remove_isotopes_from_product(int s);
  void set_write_fragments_added(int s);

  int GenerateReplacements(Molecule& m, std::vector<SubstituentReplacement>& results);
  std::vector<SubstituentReplacement> GenerateReplacements(Molecule& m);

 private:
  std::unique_ptr<SubstituentIdentification> _impl;
};

int SubstituentIdentificationMain(int argc, char** argv);

#endif  // MOLECULE_TOOLS_BDB_SUBSTITUENT_IDENTIFICATION_LIB_H_
