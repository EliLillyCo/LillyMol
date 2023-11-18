// Fingerprint functions for LillyMol

#include <cstdint>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

namespace py = pybind11;

#ifdef LILLYMOL_VECTOR_OPAQUE
#endif

#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"

// Lifted from test_numpy_dtypes.cpp
template <typename T>
py::array mkarray_via_buffer(size_t n) {
    return py::array(py::buffer_info(
        nullptr, sizeof(T), py::format_descriptor<T>::format(), 1, {n}, {sizeof(T)}));
}

template <typename T>
void
LinearFingerprintToVector(IWMFingerprint& generator,
                          T* destination)  {
  const int nbits = generator.nbits();
  const int* bits = generator.vector();

  for (int i = 0; i < nbits; ++i) {
    destination[i] = bits[i];
  }
}

PYBIND11_MODULE(lillymol, m) 
{
  m.def("linear_fingerprint",
    [](Molecule& m)->py::array_t<uint8_t> {
      static constexpr int kNBits = 2048;

      py::array_t<uint8_t> result = mkarray_via_buffer<uint8_t>(kNBits);
      auto req = result.request();
      uint8_t* ptr = static_cast<uint8_t*>(req.ptr);

      IWMFingerprint generator(kNBits);
      generator.construct_fingerprint(m);

      LinearFingerprintToVector(generator, ptr);

      return result;
    },
    "Linear path based fingerprints with default atom type"
  );
  m.def("linear_fingerprint",
    [](Molecule& m, int nbits, const std::string& atype_specification)->std::optional<py::array_t<uint8_t>> {
      Atom_Typing_Specification atom_typing;
      if (atype_specification.empty()) {
        atom_typing.build("UST:AY");
      } else {
        const const_IWSubstring tmp(atype_specification);
        if (! atom_typing.build(tmp)) {
          std::cerr << "linear_fingerprint:invalid atom type '" << atype_specification << "'\n";
          return std::nullopt;
        }
      }

      const int matoms = m.natoms();
      std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);
      if (! atom_typing.assign_atom_types(m, atype.get())) {
        std::cerr << "linear_fingerprint::Cannot assign atom types\n";
        return std::nullopt;
      }
      
      if (nbits == 0) {
        nbits = 2048;
      }

      IWMFingerprint generator(nbits);

      generator.construct_fingerprint(m);

      py::array_t<uint8_t> result = mkarray_via_buffer<uint8_t>(nbits);
      auto req = result.request();
      uint8_t* ptr = static_cast<uint8_t*>(req.ptr);

      LinearFingerprintToVector(generator, ptr);

      return result;
    },
    py::arg("m"), py::kw_only(), py::arg("nbits"), py::arg("atype_specification"),
    "Linear fingerprint with atom type and number of bits"
  );

}
