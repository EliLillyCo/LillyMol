#ifndef PYBIND_MOLECULE_H_
#define PYBIND_MOLECULE_H_

#include <optional>
#include <string>

class Molecule;
class Substructure_Query;

extern std::optional<Molecule>
MolFromSmiles(const std::string& smi);

extern std::optional<Substructure_Query>
QueryFromSmarts(const std::string& smi);

#endif  // PYBIND_MOLECULE_H_
