// Fingerprint functions for LillyMol

#include <cstdint>

#include "Foundational/iwmisc/misc.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

namespace py = pybind11;

#ifdef LILLYMOL_VECTOR_OPAQUE
#endif

#include "Molecule_Lib/atom_pair_fingerprint.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/linear_fingerprint.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/ec_fingerprint.h"

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

class BaseFpGenerator {
  protected:
    Atom_Typing_Specification _atype;

    int _nbits;

  public:
    BaseFpGenerator(int nb);

    void set_nbits(int s) {
      _nbits = s;
    }

    bool SetAtomType(const std::string& atype);

    py::array_t<uint8_t> Fingerprint(Molecule& m) = delete;
};

BaseFpGenerator::BaseFpGenerator(int nb) {
  _nbits = nb;
  _atype.build("UST:ACHY");
}

bool
BaseFpGenerator::SetAtomType(const std::string& atype) {
  IWString tmp(atype);
  if (! _atype.build(tmp)) {
    std::cerr << "BaseFpGenerator::SetAtomType:invalid atom type '" << atype << "'\n";
    return false;
  }

  return true;
}

class LinearFingerprintByte : public BaseFpGenerator {
  private:
    linear_fingerprint::LinearFingerprintGenerator _fp;

  public:
    LinearFingerprintByte(int nb);

    py::array_t<uint8_t> Fingerprint(Molecule& m);
};

LinearFingerprintByte::LinearFingerprintByte(int nb) : BaseFpGenerator(nb) {
}

py::array_t<uint8_t> 
LinearFingerprintByte::Fingerprint(Molecule& m) {
  std::unique_ptr<linear_fingerprint::atom_type_t[]> atype = std::make_unique<linear_fingerprint::atom_type_t[]>(m.natoms());
  if (! _atype.assign_atom_types(m, atype.get())) {
    std::cerr << "ECFingerprintByte::Fingerprint:cannot assign atom types " << m.name() << '\n';
  }

  py::array_t<uint8_t> result = mkarray_via_buffer<uint8_t>(_nbits);
  auto req = result.request();
  uint8_t* ptr = static_cast<uint8_t*>(req.ptr);
  std::fill_n(ptr, _nbits, 0);

  Sparse_Fingerprint_Creator sfc;
  _fp.Fingerprint(m, nullptr, atype.get(), sfc);

  // No checking overflow
  for (const auto& [b, c] : sfc.bits_found()) {
    int bit = b % _nbits;
    ptr[bit] += c;
  }

  return result;
}

class ECFingerprintByte : public BaseFpGenerator {
  private:
    ec_fingerprint::ECFingerprint _fp;

  public:
    ECFingerprintByte(int nb);

    py::array_t<uint8_t> Fingerprint(Molecule& m);
};

ECFingerprintByte::ECFingerprintByte(int nb) : BaseFpGenerator(nb) {
}

py::array_t<uint8_t> 
ECFingerprintByte::Fingerprint(Molecule& m) {
  std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(m.natoms());
  if (! _atype.assign_atom_types(m, atype.get())) {
    std::cerr << "ECFingerprintByte::Fingerprint:cannot assign atom types " << m.name() << '\n';
  }

  py::array_t<uint8_t> result = mkarray_via_buffer<uint8_t>(_nbits);
  auto req = result.request();
  uint8_t* ptr = static_cast<uint8_t*>(req.ptr);
  std::fill_n(ptr, _nbits, 0);

  ec_fingerprint::ProduceFingerprint bits;
  _fp.Fingerprint(m, nullptr, atype.get(), bits);

  for (const auto& [bit, count] : bits.sfc().bits_found()) {
    const int b = bit % _nbits;

    if (count > std::numeric_limits<uint8_t>::max()) {
      ptr[b] = std::numeric_limits<uint8_t>::max();
    } else {
      ptr[b] += count;
    }
  }

  return result;
}

class AtomPairFingerprintByte : public BaseFpGenerator {
  private:
    atom_pair_fingerprint::AtomPairFingerprint _fp;

  public:
    AtomPairFingerprintByte(int nb);

    void set_min_separation(int s) {
      _fp.set_min_separation(s);
    }
    void set_max_separation(int s) {
      _fp.set_max_separation(s);
    }

    py::array_t<uint8_t> Fingerprint(Molecule& m);
};

AtomPairFingerprintByte::AtomPairFingerprintByte(int nb) : BaseFpGenerator(nb) {
}

py::array_t<uint8_t> 
AtomPairFingerprintByte::Fingerprint(Molecule& m) {
  std::unique_ptr<uint64_t[]> atype = std::make_unique<uint64_t[]>(m.natoms());
  if (! _atype.assign_atom_types(m, atype.get())) {
    std::cerr << "AtomPairFingerprintByte::Fingerprint:cannot assign atom types " << m.name() << '\n';
  }

  py::array_t<uint8_t> result = mkarray_via_buffer<uint8_t>(_nbits);
  auto req = result.request();
  uint8_t* ptr = static_cast<uint8_t*>(req.ptr);
  std::fill_n(ptr, _nbits, 0);

  Sparse_Fingerprint_Creator sfc;
  _fp.Fingerprint(m, nullptr, atype.get(), sfc);

  // No checking for overflow...
  for (const auto& [b, c]: sfc.bits_found()) {
    int bit = b % _nbits;
    ptr[bit] += c;
  }

  return result;
}

PYBIND11_MODULE(lillymol_fingerprint, m) 
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

  py::class_<ECFingerprintByte>(m, "ECFingerprintCreator")
    .def(py::init<int>())
    .def("set_nbits", &ECFingerprintByte::set_nbits, "Set the width of fingerprints generated")
    .def("set_atom_type",
      [](ECFingerprintByte& fp_creator, const std::string& s)->bool {
        return fp_creator.SetAtomType(s);
      },
      "Set atom type"
    )
    .def("fingerprint",
      [](ECFingerprintByte& fp_creator, Molecule& m) {
        return fp_creator.Fingerprint(m);
      },
      "Geneate a fixed width counted fingerprint for `m`"
    )
  ;

  py::class_<LinearFingerprintByte>(m, "LinearFingerprintCreator")
    .def(py::init<int>())
    .def("set_nbits", &LinearFingerprintByte::set_nbits, "Set the width of fingerprints generated")
    .def("set_atom_type",
      [](LinearFingerprintByte& fp_creator, const std::string& s)->bool {
        return fp_creator.SetAtomType(s);
      },
      "Set atom type"
    )
    .def("fingerprint",
      [](LinearFingerprintByte& fp_creator, Molecule& m) {
        return fp_creator.Fingerprint(m);
      },
      "Geneate a fixed width counted fingerprint for `m`"
    )
  ;

  py::class_<AtomPairFingerprintByte>(m, "AtomPairFingerprintCreator")
    .def(py::init<int>())
    .def("set_nbits", &AtomPairFingerprintByte::set_nbits, "Set the width of fingerprints generated")
    .def("set_atom_type",
      [](AtomPairFingerprintByte& fp_creator, const std::string& s)->bool {
        return fp_creator.SetAtomType(s);
      },
      "Set atom type"
    )
    .def("fingerprint",
      [](AtomPairFingerprintByte& fp_creator, Molecule& m) {
        return fp_creator.Fingerprint(m);
      },
      "Geneate a fixed width counted fingerprint for `m`"
    )
  ;

  m.def("tanimoto",
    [](const py::array_t<uint8_t>& fp1, const py::array_t<uint8_t>& fp2)->float {
      auto req1 = fp1.request();
      uint8_t* ptr1 = static_cast<uint8_t*>(req1.ptr);
      auto req2 = fp2.request();
      uint8_t* ptr2 = static_cast<uint8_t*>(req2.ptr);

      const int nbits = req1.size;
      int bits_in_common = 0;
      int all_bits = 0;
      for (int i = 0; i < nbits; ++i) {
        uint8_t c1 = ptr1[i];
        uint8_t c2 = ptr2[i];
        if (c1 <= c2) {
          bits_in_common += c1;
          all_bits += c2;
        } else {
          bits_in_common += c2;
          all_bits += c1;
        }
      }

      if (all_bits == 0) {
        return 0;
      }

      float result = iwmisc::Fraction<float>(bits_in_common, all_bits);

      return result;
    },
    "Compute Tanimito coefficient between two numpy vectors"
  );
}
