// Fingerprint bindings for the nanobind LillyMol pilot.

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

#include <nanobind/ndarray.h>

#include "Foundational/iwmisc/misc.h"
#include "Molecule_Lib/atom_pair_fingerprint.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/linear_fingerprint.h"
#include "Molecule_Tools/ec_fingerprint.h"
#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {
namespace {

using IntVector = std::vector<int>;
using UInt8NumpyArray = nb::ndarray<nb::numpy, uint8_t, nb::shape<-1>>;

IntVector
MakeIntVector(size_t n) {
  return IntVector(n, 0);
}

template <typename T>
void
LinearFingerprintToVector(IWMFingerprint& generator, T* destination) {
  const int nbits = generator.nbits();
  const int* bits = generator.vector();

  for (int i = 0; i < nbits; ++i) {
    destination[i] = bits[i];
  }
}

UInt8NumpyArray
MakeUInt8NumpyArray(size_t n) {
  uint8_t* data = new uint8_t[n];
  std::fill_n(data, n, uint8_t{0});
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<uint8_t*>(p); });
  return UInt8NumpyArray(data, {n}, owner);
}

void
IncrementUInt8(uint8_t& destination, uint32_t count) {
  const uint32_t available =
      static_cast<uint32_t>(std::numeric_limits<uint8_t>::max() - destination);
  if (count >= available) {
    destination = std::numeric_limits<uint8_t>::max();
  } else {
    destination += static_cast<uint8_t>(count);
  }
}

void
SparseFingerprintToUInt8Array(const Sparse_Fingerprint_Creator& sfc, UInt8NumpyArray& destination) {
  const int nbits = destination.shape(0);
  uint8_t* data = destination.data();
  for (const auto& [bit, count] : sfc.bits_found()) {
    IncrementUInt8(data[bit % nbits], count);
  }
}

void
LinearFingerprintToUInt8Array(IWMFingerprint& generator, UInt8NumpyArray& destination) {
  const int nbits = generator.nbits();
  const int* bits = generator.vector();
  uint8_t* data = destination.data();

  for (int i = 0; i < nbits; ++i) {
    data[i] = bits[i] > std::numeric_limits<uint8_t>::max()
                  ? std::numeric_limits<uint8_t>::max()
                  : static_cast<uint8_t>(bits[i]);
  }
}

class BaseFpGenerator {
 protected:
  Atom_Typing_Specification _atype;
  int _nbits;

 public:
  explicit BaseFpGenerator(int nbits) : _nbits(nbits) {
    _atype.build("UST:ACHY");
  }

  void set_nbits(int nbits) {
    _nbits = nbits;
  }

  bool SetAtomType(const std::string& atype) {
    IWString tmp(atype);
    if (!_atype.build(tmp)) {
      std::cerr << "BaseFpGenerator::SetAtomType:invalid atom type '" << atype << "'\n";
      return false;
    }

    return true;
  }
};

class LinearFingerprintByte : public BaseFpGenerator {
 private:
  linear_fingerprint::LinearFingerprintGenerator _fp;

 public:
  explicit LinearFingerprintByte(int nbits) : BaseFpGenerator(nbits) {}

  void set_max_length(uint32_t length) {
    _fp.set_max_length(length);
  }

  IntVector Fingerprint(Molecule& mol) {
    std::unique_ptr<linear_fingerprint::atom_type_t[]> atype =
        std::make_unique<linear_fingerprint::atom_type_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("LinearFingerprintCreator cannot assign atom types");
    }

    IntVector result = MakeIntVector(_nbits);
    int* ptr = result.data();

    Sparse_Fingerprint_Creator sfc;
    _fp.Fingerprint(mol, nullptr, atype.get(), sfc);

    for (const auto& [bit, count] : sfc.bits_found()) {
      ptr[bit % _nbits] += count;
    }

    return result;
  }

  UInt8NumpyArray FingerprintNumpy(Molecule& mol) {
    std::unique_ptr<linear_fingerprint::atom_type_t[]> atype =
        std::make_unique<linear_fingerprint::atom_type_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("LinearFingerprintCreator cannot assign atom types");
    }

    UInt8NumpyArray result = MakeUInt8NumpyArray(_nbits);
    Sparse_Fingerprint_Creator sfc;
    _fp.Fingerprint(mol, nullptr, atype.get(), sfc);
    SparseFingerprintToUInt8Array(sfc, result);

    return result;
  }
};

class ECFingerprintByte : public BaseFpGenerator {
 private:
  ec_fingerprint::ECFingerprint _fp;

 public:
  explicit ECFingerprintByte(int nbits) : BaseFpGenerator(nbits) {}

  void set_max_radius(uint32_t radius) {
    _fp.set_max_radius(radius);
  }

  IntVector Fingerprint(Molecule& mol) {
    std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("ECFingerprintCreator cannot assign atom types");
    }

    IntVector result = MakeIntVector(_nbits);
    int* ptr = result.data();

    ec_fingerprint::ProduceFingerprint bits;
    _fp.Fingerprint(mol, nullptr, atype.get(), bits);

    for (const auto& [bit, count] : bits.sfc().bits_found()) {
      const int b = bit % _nbits;
      if (count > std::numeric_limits<int>::max()) {
        ptr[b] = std::numeric_limits<int>::max();
      } else {
        ptr[b] += count;
      }
    }

    return result;
  }

  UInt8NumpyArray FingerprintNumpy(Molecule& mol) {
    std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("ECFingerprintCreator cannot assign atom types");
    }

    UInt8NumpyArray result = MakeUInt8NumpyArray(_nbits);
    ec_fingerprint::ProduceFingerprint bits;
    _fp.Fingerprint(mol, nullptr, atype.get(), bits);
    SparseFingerprintToUInt8Array(bits.sfc(), result);

    return result;
  }
};

class AtomPairFingerprintByte : public BaseFpGenerator {
 private:
  atom_pair_fingerprint::AtomPairFingerprint _fp;

 public:
  explicit AtomPairFingerprintByte(int nbits) : BaseFpGenerator(nbits) {}

  void set_min_separation(int separation) {
    _fp.set_min_separation(separation);
  }
  void set_max_separation(int separation) {
    _fp.set_max_separation(separation);
  }

  IntVector Fingerprint(Molecule& mol) {
    std::unique_ptr<uint64_t[]> atype = std::make_unique<uint64_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("AtomPairFingerprintCreator cannot assign atom types");
    }

    IntVector result = MakeIntVector(_nbits);
    int* ptr = result.data();

    Sparse_Fingerprint_Creator sfc;
    _fp.Fingerprint(mol, nullptr, atype.get(), sfc);

    for (const auto& [bit, count] : sfc.bits_found()) {
      ptr[bit % _nbits] += count;
    }

    return result;
  }

  UInt8NumpyArray FingerprintNumpy(Molecule& mol) {
    std::unique_ptr<uint64_t[]> atype = std::make_unique<uint64_t[]>(mol.natoms());
    if (!_atype.assign_atom_types(mol, atype.get())) {
      throw std::runtime_error("AtomPairFingerprintCreator cannot assign atom types");
    }

    UInt8NumpyArray result = MakeUInt8NumpyArray(_nbits);
    Sparse_Fingerprint_Creator sfc;
    _fp.Fingerprint(mol, nullptr, atype.get(), sfc);
    SparseFingerprintToUInt8Array(sfc, result);

    return result;
  }
};

IntVector
LinearFingerprintDefault(Molecule& mol) {
  static constexpr int kNBits = 2048;
  IntVector result = MakeIntVector(kNBits);
  int* ptr = result.data();

  IWMFingerprint generator(kNBits);
  generator.construct_fingerprint(mol);
  LinearFingerprintToVector(generator, ptr);

  return result;
}

UInt8NumpyArray
LinearFingerprintDefaultNumpy(Molecule& mol) {
  static constexpr int kNBits = 2048;
  UInt8NumpyArray result = MakeUInt8NumpyArray(kNBits);

  IWMFingerprint generator(kNBits);
  generator.construct_fingerprint(mol);
  LinearFingerprintToUInt8Array(generator, result);

  return result;
}

std::optional<IntVector>
LinearFingerprintWithOptions(Molecule& mol, int nbits, const std::string& atype_specification) {
  Atom_Typing_Specification atom_typing;
  if (atype_specification.empty()) {
    atom_typing.build("UST:AY");
  } else {
    const const_IWSubstring tmp(atype_specification);
    if (!atom_typing.build(tmp)) {
      std::cerr << "linear_fingerprint:invalid atom type '" << atype_specification << "'\n";
      return std::nullopt;
    }
  }

  const int matoms = mol.natoms();
  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);
  if (!atom_typing.assign_atom_types(mol, atype.get())) {
    std::cerr << "linear_fingerprint::Cannot assign atom types\n";
    return std::nullopt;
  }

  if (nbits == 0) {
    nbits = 2048;
  }

  IntVector result = MakeIntVector(nbits);
  int* ptr = result.data();

  IWMFingerprint generator(nbits);
  generator.construct_fingerprint(mol);
  LinearFingerprintToVector(generator, ptr);

  return result;
}


float
Tanimoto(const IntVector& fp1, const IntVector& fp2) {
  if (fp1.size() != fp2.size()) {
    throw std::invalid_argument("fingerprints must be arrays of equal length");
  }

  const int* ptr1 = fp1.data();
  const int* ptr2 = fp2.data();
  const size_t nbits = fp1.size();
  int bits_in_common = 0;
  int all_bits = 0;
  for (size_t i = 0; i < nbits; ++i) {
    const int c1 = ptr1[i];
    const int c2 = ptr2[i];
    if (c1 <= c2) {
      bits_in_common += c1;
      all_bits += c2;
    } else {
      bits_in_common += c2;
      all_bits += c1;
    }
  }

  if (all_bits == 0) {
    return 0.0f;
  }

  return iwmisc::Fraction<float>(bits_in_common, all_bits);
}

}  // namespace

void
BindFingerprint(nb::module_& m) {
  m.def("linear_fingerprint", &LinearFingerprintDefault,
        "Linear path based fingerprints with default atom type");
  m.def("linear_fingerprint_numpy", &LinearFingerprintDefaultNumpy,
        "Linear path based fingerprint as a NumPy uint8 array");
  m.def("linear_fingerprint", &LinearFingerprintWithOptions,
        nb::arg("m"), nb::kw_only(), nb::arg("nbits"), nb::arg("atype_specification"),
        "Linear fingerprint with atom type and number of bits");

  nb::class_<ECFingerprintByte>(m, "ECFingerprintCreator")
      .def(nb::init<int>())
      .def("set_nbits", &ECFingerprintByte::set_nbits,
           "Set the width of fingerprints generated")
      .def("set_atom_type", &ECFingerprintByte::SetAtomType,
           "Set atom type")
      .def("set_max_radius", &ECFingerprintByte::set_max_radius,
           "Max radius for fingerprints")
      .def("fingerprint", &ECFingerprintByte::Fingerprint,
           "Generate a fixed width counted fingerprint for a molecule")
      .def("fingerprint_numpy", &ECFingerprintByte::FingerprintNumpy,
           "Generate a fixed width counted fingerprint as a NumPy uint8 array");

  nb::class_<LinearFingerprintByte>(m, "LinearFingerprintCreator")
      .def(nb::init<int>())
      .def("set_nbits", &LinearFingerprintByte::set_nbits,
           "Set the width of fingerprints generated")
      .def("set_atom_type", &LinearFingerprintByte::SetAtomType,
           "Set atom type")
      .def("set_max_length", &LinearFingerprintByte::set_max_length,
           "Max path length")
      .def("fingerprint", &LinearFingerprintByte::Fingerprint,
           "Generate a fixed width counted fingerprint for a molecule")
      .def("fingerprint_numpy", &LinearFingerprintByte::FingerprintNumpy,
           "Generate a fixed width counted fingerprint as a NumPy uint8 array");

  nb::class_<AtomPairFingerprintByte>(m, "AtomPairFingerprintCreator")
      .def(nb::init<int>())
      .def("set_nbits", &AtomPairFingerprintByte::set_nbits,
           "Set the width of fingerprints generated")
      .def("set_atom_type", &AtomPairFingerprintByte::SetAtomType,
           "Set atom type")
      .def("set_min_separation", &AtomPairFingerprintByte::set_min_separation,
           "min separation between atoms")
      .def("set_max_separation", &AtomPairFingerprintByte::set_max_separation,
           "max separation between atoms")
      .def("fingerprint", &AtomPairFingerprintByte::Fingerprint,
           "Generate a fixed width counted fingerprint for a molecule")
      .def("fingerprint_numpy", &AtomPairFingerprintByte::FingerprintNumpy,
           "Generate a fixed width counted fingerprint as a NumPy uint8 array");

  m.def("tanimoto", &Tanimoto,
        "Compute Tanimoto coefficient between two counted integer sequences");
}

}  // namespace lillymol_nb
