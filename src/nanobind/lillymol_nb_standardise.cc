#include "nanobind/lillymol_nb_internal.h"

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"

namespace lillymol_nb {
namespace {

void
InitialiseChargeAssigner(Charge_Assigner& charge_assigner, const IWString& dirname) {
  if (!charge_assigner.BuildFromDir(dirname)) {
    throw std::runtime_error("ChargeAssigner:cannot initialise from '" + dirname.AsString() + "'");
  }
}

void
InitialiseChargeAssignerFromDefaultEnv(Charge_Assigner& charge_assigner) {
  if (!charge_assigner.BuildFromDefaultEnvs()) {
    throw std::runtime_error(
        "ChargeAssigner:cannot initialise from C3TK_DATA_PERSISTENT or LILLYMOL_HOME");
  }
}

int
AssignCharges(Charge_Assigner& charge_assigner, Molecule& mol) {
  const int matoms = mol.natoms();
  if (matoms == 0) {
    return 0;
  }

  std::unique_ptr<formal_charge_t[]> charges_assigned =
      std::make_unique<formal_charge_t[]>(matoms);
  if (!charge_assigner.process(mol, charges_assigned.get())) {
    return 0;
  }

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (charges_assigned[i] == 0) {
      continue;
    }
    mol.set_formal_charge(i, charges_assigned[i]);
    ++rc;
  }

  return rc;
}

void
InitialiseDonorAcceptor(Donor_Acceptor_Assigner& donor_acceptor, const IWString& dirname) {
  static constexpr int kVerbose = 0;
  if (!donor_acceptor.BuildFromDir(dirname, kVerbose)) {
    throw std::runtime_error("DonorAcceptor:cannot initialise from '" +
                             dirname.AsString() + "'");
  }
  donor_acceptor.set_apply_isotopic_labels(1);
}

void
InitialiseDonorAcceptorFromDefaultEnv(Donor_Acceptor_Assigner& donor_acceptor) {
  static constexpr int kVerbose = 0;
  if (!donor_acceptor.BuildFromDefaultEnv(kVerbose)) {
    throw std::runtime_error(
        "DonorAcceptor:cannot initialise from C3TK_DATA_PERSISTENT or LILLYMOL_HOME");
  }
  donor_acceptor.set_apply_isotopic_labels(1);
}

}  // namespace

void
BindStandardise(nb::module_& m) {
  nb::class_<Chemical_Standardisation>(m, "Standardise")
      .def(nb::init<>())
      .def("activate_all", &Chemical_Standardisation::activate_all,
           "Activate all transformations")
      .def("set_verbose", &Chemical_Standardisation::set_verbose,
           "Set verbosity")
      .def("process", &Chemical_Standardisation::process,
           "Apply active transformations to molecule");

  nb::class_<Element_Transformations>(m, "ElementTransformations")
      .def(nb::init<>())
      .def("active", [](const Element_Transformations& etrans) {
        return static_cast<bool>(etrans.active());
      }, "True if any transformations have been added")
      .def("add",
           [](Element_Transformations& etrans, const std::string& directive) {
             const IWString tmp(directive);
             return static_cast<bool>(etrans.Add(tmp));
           },
           nb::arg("directive"), "Add a transformation directive like 'I=Cl'")
      .def("process",
           [](Element_Transformations& etrans, Molecule& mol) { return etrans.process(mol); },
           nb::arg("mol"), "Apply transformations to molecule");

  nb::class_<Charge_Assigner>(m, "ChargeAssigner")
      .def("__init__",
           [](Charge_Assigner* charge_assigner) {
             new (charge_assigner) Charge_Assigner();
             InitialiseChargeAssignerFromDefaultEnv(*charge_assigner);
           },
           "Build charge queries from C3TK_DATA_PERSISTENT or LILLYMOL_HOME")
      .def("__init__",
           [](Charge_Assigner* charge_assigner, const std::string& query_dir) {
             new (charge_assigner) Charge_Assigner();
             InitialiseChargeAssigner(*charge_assigner, IWString(query_dir));
           },
           nb::arg("query_dir"),
           "Build charge queries from an explicit charges query directory")
      .def("active",
           [](const Charge_Assigner& charge_assigner) {
             return static_cast<bool>(charge_assigner.active());
           },
           "True if charge assignment queries are loaded")
      .def("set_min_distance_between_charges",
           &Charge_Assigner::set_min_distance_between_charges,
           nb::arg("distance"),
           "Specify minimum bond separation between formal charges assigned")
      .def("process", &AssignCharges, nb::arg("mol"),
           "Assign formal charges to a molecule and return the number of changed atoms");

  nb::class_<Donor_Acceptor_Assigner>(m, "DonorAcceptor")
      .def("__init__",
           [](Donor_Acceptor_Assigner* donor_acceptor) {
             new (donor_acceptor) Donor_Acceptor_Assigner();
             InitialiseDonorAcceptorFromDefaultEnv(*donor_acceptor);
           },
           "Build donor/acceptor queries from C3TK_DATA_PERSISTENT or LILLYMOL_HOME")
      .def("__init__",
           [](Donor_Acceptor_Assigner* donor_acceptor, const std::string& query_dir) {
             new (donor_acceptor) Donor_Acceptor_Assigner();
             InitialiseDonorAcceptor(*donor_acceptor, IWString(query_dir));
           },
           nb::arg("query_dir"),
           "Build donor/acceptor queries from an explicit hbonds query directory")
      .def("active",
           [](const Donor_Acceptor_Assigner& donor_acceptor) {
             return static_cast<bool>(donor_acceptor.active());
           },
           "True if donor or acceptor queries are loaded")
      .def("process",
           [](Donor_Acceptor_Assigner& donor_acceptor, Molecule& mol) {
             return donor_acceptor.process(mol);
           },
           nb::arg("mol"), "Assign donor/acceptor isotopic labels to a molecule");
}

}  // namespace lillymol_nb
