#include "unique_molecules_api.h"

namespace unique_molecules {
UniqueMolecules::UniqueMolecules() {
  _exclude_chiral_info = 0;
  _exclude_cis_trans_bonding_info = 1;
  _strip_to_largest_fragment = 0;
  _ignore_isotopes = 0;
  _constant_isotope = 0;

  _use_atom_hash = 1;

  _molecules_processed = 0;
  _duplicates_found = 0;
}

void
UniqueMolecules::AddToHashes(const Molecule& m) {
  InternalIsUnique(Molecule(m));
}

bool
UniqueMolecules::IsUnique(const Molecule& m) {
  ++_molecules_processed;

  // Create copy that may be altered.
  if (InternalIsUnique(Molecule(m))) {
    return true;
  }

  ++_duplicates_found;
  return false;
}

bool
UniqueMolecules::InternalIsUnique(Molecule&& m) {
  if (_strip_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
  if (_exclude_chiral_info) {
    m.remove_all_chiral_centres();
  }
  if (_exclude_cis_trans_bonding_info) {
    m.revert_all_directional_bonds_to_non_directional();
  }
  if (_ignore_isotopes) {
    m.unset_isotopes();
  }
  if (_etrans.active()) {
    _etrans.process(m);
  }
  if (_constant_isotope) {
  }

  return InternalIsUniqueInner(m);
}

bool
UniqueMolecules::InternalIsUniqueInner(Molecule& m) {

  IWString usmi;

  if (_mol2graph.active()) {
    IWString formula = m.molecular_formula();

    m.change_to_graph_form();

    usmi = m.unique_smiles();

    usmi << ':' << formula;
  } else {
    usmi = m.unique_smiles();
  }

  std::string usmi_string(usmi.data(), usmi.size());

  //cerr << "Testing unique smiles '" << usmi << "'\n";

  const int matoms = m.natoms();

  while (matoms >= _smiles_hash.number_elements())
  {
    _smiles_hash.add(new absl::flat_hash_set<std::string>);
    if (_use_atom_hash)
      _atom_hash.add(new std::unordered_set<uint64_t>());
  }

  //cerr << matoms << " smiles_hash " << _smiles_hash.number_elements() << " formula_hash " << formula_hash.number_elements() << endl;

  if (_use_atom_hash) {
    const uint64_t h = m.quick_bond_hash();

    std::unordered_set<uint64_t>* ha = _atom_hash[matoms];

    const auto f = ha->find(h);

    if (f == ha->end()) {  // new
      ha->insert(h);

      absl::flat_hash_set<std::string> * h = _smiles_hash[matoms];
      h->emplace(std::move(usmi_string));

      return true;
    }
  }

  absl::flat_hash_set<std::string>& h = *_smiles_hash[matoms];

  auto iter = h.find(usmi_string);

  if (iter == h.end()) {
    h.emplace(std::move(usmi_string));
    return true;
  }

  return false;
}

int
UniqueMolecules::Report(std::ostream& output) const {
  uint32_t tot = 0;
  for (const auto * s : _smiles_hash) {
    tot += s->size();
  }
  output << "UniqueMolecules:hash contains " << _smiles_hash.size() << " atom counts with " << tot << " unique smiles\n";
  output << _molecules_processed << " molecules " << _duplicates_found << " duplicates " << (_molecules_processed - _duplicates_found) << " unique structures\n";

  return 1;
}

}  // namespace unique_molecules
