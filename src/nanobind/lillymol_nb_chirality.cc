#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

namespace {

std::vector<int>
ChiralCentreAtoms(const Chiral_Centre& centre) {
  return {centre.top_front(), centre.top_back(), centre.left_down(), centre.right_down()};
}

int
ChiralCentreAtom(const Chiral_Centre& centre, int index) {
  switch (index) {
    case 0:
      return centre.top_front();
    case 1:
      return centre.top_back();
    case 2:
      return centre.left_down();
    case 3:
      return centre.right_down();
    default:
      throw std::out_of_range("Chiral_Centre atom index outside [0, 4)");
  }
}

}  // namespace

void
BindChirality(nb::module_& m) {
  m.attr("CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN") =
      nb::int_(kChiralConnectionIsImplicitHydrogen);
  m.attr("CHIRAL_CONNECTION_IS_LONE_PAIR") = nb::int_(kChiralConnectionIsLonePair);

  nb::enum_<ChiralType>(m, "ChiralType")
      .value("CHI_UNSPECIFIED", ChiralType::kChiUnspecified)
      .value("CHI_TETRAHEDRAL_CW", ChiralType::kChiTetrahedralCw)
      .value("CHI_TETRAHEDRAL_CCW", ChiralType::kChiTetrahedralCcw)
      .value("CHI_OTHER", ChiralType::kChiOther);

  nb::enum_<CahnIngoldPrelog>(m, "CIP")
      .value("R", CahnIngoldPrelog::R)
      .value("S", CahnIngoldPrelog::S)
      .value("Neither", CahnIngoldPrelog::kNeither)
      .value("Unspecified", CahnIngoldPrelog::kUnspecified);

  nb::class_<Chiral_Centre>(m, "Chiral_Centre")
      .def(nb::init<atom_number_t>(), nb::arg("atom"))
      .def("atom", &Chiral_Centre::a)
      .def("top_front", &Chiral_Centre::top_front)
      .def("top_back", &Chiral_Centre::top_back)
      .def("left_down", &Chiral_Centre::left_down)
      .def("right_down", &Chiral_Centre::right_down)
      .def("atoms", &ChiralCentreAtoms,
           "Return [top_front, top_back, left_down, right_down]")
      .def("__len__", [](const Chiral_Centre&) { return 4; })
      .def("__getitem__", &ChiralCentreAtom)
      .def("invert", &Chiral_Centre::invert)
      .def("involves",
           [](const Chiral_Centre& centre, atom_number_t atom) {
             return static_cast<bool>(centre.involves(atom));
           },
           nb::arg("atom"))
      .def("implicit_hydrogens", &Chiral_Centre::implicit_hydrogen_count)
      .def("lone_pairs", &Chiral_Centre::lone_pair_count)
      .def("__repr__", &ChiralCentreRepr);

  m.def("tetrahedral_chirality", &TetrahedralChirality,
        nb::arg("mol"), nb::arg("atom"), nb::arg("check_is_chiral") = false);
  m.def("is_actually_chiral",
        [](Molecule& mol, atom_number_t atom) {
          if (atom < 0 || atom >= mol.natoms()) {
            throw std::invalid_argument("atom number outside molecule");
          }
          return static_cast<bool>(::is_actually_chiral(mol, atom));
        },
        nb::arg("mol"), nb::arg("atom"));
  m.def("is_chiral_implicit_hydrogen",
        [](int connection) { return IsChiralImplicitHydrogen(connection); },
        nb::arg("connection"));
  m.def("is_chiral_lone_pair",
        [](int connection) { return IsChiralLonePair(connection); },
        nb::arg("connection"));
}

}  // namespace lillymol_nb
