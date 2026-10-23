#include <optional>
#include <tuple>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include "Depict/coords_2d.h"
#include "Depict/wedge_2d.h"
#include "Molecule_Lib/molecule.h"

namespace nb = nanobind;

namespace {

depict::Coords2DOptions
MakeOptions(float bond_length, depict::Precision precision,
            bool skip_minimization, bool even_angles, bool centre,
            bool honour_cis_trans) {
  depict::Coords2DOptions result;
  result.bond_length = bond_length;
  result.precision = precision;
  result.skip_minimization = skip_minimization;
  result.even_angles = even_angles;
  result.centre = centre;
  result.honour_cis_trans = honour_cis_trans;
  return result;
}

std::optional<depict::Coords2DResult>
PyGenerate2DCoordinates(Molecule& mol, const depict::Coords2DOptions& options) {
  depict::Coords2DResult result;
  if (!depict::Generate2DCoordinates(mol, options, result)) {
    return std::nullopt;
  }
  return result;
}

std::optional<std::tuple<Molecule, depict::Coords2DResult>>
PyGenerate2DCoordinatesCopy(const Molecule& mol,
                            const depict::Coords2DOptions& options) {
  Molecule copy(mol);
  std::optional<depict::Coords2DResult> result =
      PyGenerate2DCoordinates(copy, options);
  if (!result) {
    return std::nullopt;
  }
  return std::make_tuple(std::move(copy), *result);
}

depict::WedgeResult
AssignWedgeBonds(Molecule& mol) {
  depict::WedgeResult result;
  depict::AssignWedgeBonds(mol, result);
  return result;
}

}  // namespace

NB_MODULE(lillymol_depict, m) {
  // Molecule is registered by the main extension. Import it here so callers
  // can safely import lillymol_depict first.
  nb::module_::import_("lillymol");

  nb::enum_<depict::Precision>(m, "Precision")
      .value("QUICK", depict::Precision::kQuick)
      .value("STANDARD", depict::Precision::kStandard)
      .value("BEST", depict::Precision::kBest);

  nb::class_<depict::Coords2DOptions>(m, "Coords2DOptions")
      .def("__init__",
           [](depict::Coords2DOptions* options, float bond_length,
              depict::Precision precision, bool skip_minimization,
              bool even_angles, bool centre, bool honour_cis_trans) {
             new (options) depict::Coords2DOptions(
                 MakeOptions(bond_length, precision, skip_minimization,
                             even_angles, centre, honour_cis_trans));
           },
           nb::arg("bond_length") = 1.5f,
           nb::arg("precision") = depict::Precision::kStandard,
           nb::arg("skip_minimization") = false,
           nb::arg("even_angles") = false, nb::arg("centre") = true,
           nb::arg("honour_cis_trans") = true)
      .def_rw("bond_length", &depict::Coords2DOptions::bond_length)
      .def_rw("precision", &depict::Coords2DOptions::precision)
      .def_rw("skip_minimization",
              &depict::Coords2DOptions::skip_minimization)
      .def_rw("even_angles", &depict::Coords2DOptions::even_angles)
      .def_rw("centre", &depict::Coords2DOptions::centre)
      .def_rw("honour_cis_trans",
              &depict::Coords2DOptions::honour_cis_trans);

  nb::class_<depict::Coords2DResult>(m, "Coords2DResult")
      .def_ro("closest_nonbonded_approach",
              &depict::Coords2DResult::closest_nonbonded_approach)
      .def_ro("cis_trans_bonds", &depict::Coords2DResult::cis_trans_bonds)
      .def_ro("cis_trans_honoured",
              &depict::Coords2DResult::cis_trans_honoured);

  nb::class_<depict::WedgeResult>(m, "WedgeResult")
      .def_ro("chiral_centres", &depict::WedgeResult::chiral_centres)
      .def_ro("wedged", &depict::WedgeResult::wedged)
      .def_ro("not_stereogenic", &depict::WedgeResult::not_stereogenic)
      .def_ro("unresolved", &depict::WedgeResult::unresolved);

  m.def("generate_2d_coordinates", &PyGenerate2DCoordinates,
        nb::arg("mol"), nb::arg("options"),
        "Replace a molecule's coordinates with a 2D layout; return None on failure");
  m.def(
      "generate_2d_coordinates",
      [](Molecule& mol, float bond_length, depict::Precision precision,
         bool skip_minimization, bool even_angles, bool centre,
         bool honour_cis_trans) {
        return PyGenerate2DCoordinates(
            mol, MakeOptions(bond_length, precision, skip_minimization,
                             even_angles, centre, honour_cis_trans));
      },
      nb::arg("mol"), nb::arg("bond_length") = 1.5f,
      nb::arg("precision") = depict::Precision::kStandard,
      nb::arg("skip_minimization") = false,
      nb::arg("even_angles") = false, nb::arg("centre") = true,
      nb::arg("honour_cis_trans") = true,
      "Replace a molecule's coordinates with a configurable 2D layout");

  m.def("generate_2d_coordinates_copy", &PyGenerate2DCoordinatesCopy,
        nb::arg("mol"), nb::arg("options"),
        "Return (depicted molecule, result), leaving the input unchanged");
  m.def(
      "generate_2d_coordinates_copy",
      [](const Molecule& mol, float bond_length, depict::Precision precision,
         bool skip_minimization, bool even_angles, bool centre,
         bool honour_cis_trans) {
        return PyGenerate2DCoordinatesCopy(
            mol, MakeOptions(bond_length, precision, skip_minimization,
                             even_angles, centre, honour_cis_trans));
      },
      nb::arg("mol"), nb::arg("bond_length") = 1.5f,
      nb::arg("precision") = depict::Precision::kStandard,
      nb::arg("skip_minimization") = false,
      nb::arg("even_angles") = false, nb::arg("centre") = true,
      nb::arg("honour_cis_trans") = true,
      "Return a depicted copy and result without changing the input molecule");

  m.def("assign_wedge_bonds", &AssignWedgeBonds, nb::arg("mol"),
        "Assign verified wedge/hash bonds for existing 2D coordinates");
}
