// Molecule-preserving ring precedent queries for the BerkeleyDB module.
#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "Molecule_Tools_Bdb/ring_system_lookup.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

namespace {

class RingDatabase {
 private:
  std::unique_ptr<smi2rings_bdb::RingSystemLookup> _lookup;

 public:
  RingDatabase(nb::object databases, nb::object standardize, bool remove_isotopes,
               bool remove_chirality, bool largest_fragment, bool attachment_labels,
               bool bare_ring_fallback, int max_ring_size, int max_rings_per_system,
               const std::string& atom_typing)
      : _lookup(std::make_unique<smi2rings_bdb::RingSystemLookup>()) {
    if (max_ring_size < 3 || max_rings_per_system < 1) {
      throw nb::value_error("max_ring_size must be >= 3 and max_rings_per_system >= 1");
    }
    if (!attachment_labels && !bare_ring_fallback && atom_typing.empty()) {
      throw nb::value_error("At least one ring representation must be enabled");
    }
    auto& extraction = _lookup->extractor();
    auto& options = extraction.options();
    options.set_label_attachment_points(attachment_labels);
    options.set_always_process_unsubstituted_ring(bare_ring_fallback);
    options.set_max_ring_size_to_consider(max_ring_size);
    options.set_max_ring_system_size_to_consider(max_rings_per_system);
    if (!atom_typing.empty()) {
      if (!extraction.atom_typing().build(IWString(atom_typing))) {
        throw nb::value_error("Invalid atom_typing specification");
      }
      options.set_label_attachment_points_with_environment(1);
    }
    auto& preprocessing = _lookup->preprocessor();
    preprocessing.options().set_remove_isotopes(remove_isotopes);
    preprocessing.options().set_remove_chirality(remove_chirality);
    preprocessing.options().set_largest_fragment(largest_fragment);
    std::vector<std::string> directives;
    if (nb::isinstance<nb::bool_>(standardize)) {
      if (nb::cast<bool>(standardize)) {
        directives.emplace_back("all");
      }
    } else if (nb::isinstance<nb::str>(standardize)) {
      directives.push_back(nb::cast<std::string>(standardize));
    } else {
      directives = nb::cast<std::vector<std::string>>(standardize);
    }
    for (const auto& directive : directives) {
      if (!preprocessing.SetStandardisation(IWString(directive))) {
        throw nb::value_error("Invalid standardisation directive");
      }
    }
    std::vector<std::string> paths;
    if (nb::isinstance<nb::str>(databases)) {
      paths.push_back(nb::cast<std::string>(databases));
    } else {
      paths = nb::cast<std::vector<std::string>>(databases);
    }
    if (paths.empty()) {
      throw nb::value_error("At least one database path is required");
    }
    if (!_lookup->Open(paths)) {
      throw std::runtime_error("Cannot open ring databases for reading");
    }
  }

  std::optional<uint32_t>
  RarestRingCount(const Molecule& mol) {
    if (!_lookup) {
      throw std::runtime_error("RingDatabase is closed");
    }
    std::vector<smi2rings_bdb::RingSystemMatch> results;
    if (!_lookup->Lookup(mol, results)) {
      throw std::runtime_error("Ring database lookup failed");
    }
    if (results.empty()) {
      return std::nullopt;
    }
    uint32_t count = std::numeric_limits<uint32_t>::max();
    for (const auto& result : results) {
      count = std::min(count, result._found ? result._database._count : 0u);
    }
    return count;
  }

  void
  Close() {
    _lookup.reset();
  }
};
}  // namespace

void
RegisterRingDatabase(nb::module_& m) {
  nb::class_<RingDatabase>(m, "RingDatabase")
      .def(
          "__init__",
          [](RingDatabase* self, nb::object databases, nb::object standardize,
             bool remove_isotopes, bool remove_chirality, bool largest_fragment,
             bool attachment_labels, bool bare_ring_fallback, int max_ring_size,
             int max_rings_per_system, const std::string& atom_typing) {
            new (self)
                RingDatabase(databases, standardize, remove_isotopes, remove_chirality,
                             largest_fragment, attachment_labels, bare_ring_fallback,
                             max_ring_size, max_rings_per_system, atom_typing);
          },
          nb::arg("databases"), nb::kw_only(), nb::arg("standardize") = false,
          nb::arg("remove_isotopes") = false, nb::arg("remove_chirality") = false,
          nb::arg("largest_fragment") = false, nb::arg("attachment_labels") = true,
          nb::arg("bare_ring_fallback") = false,
          nb::arg("max_ring_size") = std::numeric_limits<int>::max(),
          nb::arg("max_rings_per_system") = std::numeric_limits<int>::max(),
          nb::arg("atom_typing") = "")
      .def("rarest_ring_count", &RingDatabase::RarestRingCount, nb::arg("mol"),
           "Minimum ring-system precedent count; zero for a missing ring, None "
           "when no eligible ring systems exist. Preserves the input molecule.")
      .def("close", &RingDatabase::Close, "Close databases. Idempotent.")
      .def("__enter__", [](nb::object self) { return self; })
      .def(
          "__exit__",
          [](RingDatabase& self, nb::object, nb::object, nb::object) {
            self.Close();
            return false;
          },
          nb::arg("exc_type").none(), nb::arg("exc_value").none(),
          nb::arg("traceback").none());
}
