// Bindings for selected tools

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Tools/unique_molecules_api.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_tools, m) 
{
  using unique_molecules::UniqueMolecules;

  // This is a sub-optimal implementation. While functional, it is not efficient.
  // Much to my surprise I found that storing smiles in a python set() was much
  // faster than storing those same strings in the C++ map used in the current
  // implementation. TODO:ianwatson understand what is going on.
  // For now this is quite usable, just not efficient.
  py::class_<unique_molecules::UniqueMolecules>(m, "UniqueMolecules")
    .def(py::init<>())
    .def("set_exclude_chiral_info", &UniqueMolecules::set_exclude_chiral_info, "Ignore chirality")
    .def("set_exclude_cis_trans_bonding_info", &UniqueMolecules::set_exclude_cis_trans_bonding_info, "Ignore cis trans bonds")
    .def("set_strip_to_largest_fragment", &UniqueMolecules::set_strip_to_largest_fragment, "Strip to largest fragment")
    .def("set_ignore_isotopes", &UniqueMolecules::set_ignore_isotopes, "Ignore isotopes")
    .def("set_constant_isotope", &UniqueMolecules::set_constant_isotope, "Convert all isotopes to a constant")
    .def("element_transformations", 
      [](UniqueMolecules& u)->Element_Transformations&{
        return u.element_transformations();
      },
      py::return_value_policy::reference,
      "Element transformations"
    )
    // Element Transformations temporarily not working - need to figure out duplicate Element ptr issue
#ifdef ADD_ELEMENT_TRANSFORMATION_NOT_YET_WORKING
    .def("add_element_transformation",
      [](UniqueMolecules& u, const std::string& trans)->bool{
        const IWString s(trans);
        return u.element_transformations().Add(s);
      },
      "Add one element transformation"
    )
#endif
    .def("graph_specifications",
      [](UniqueMolecules& u)->Mol2Graph&{
        return u.graph_specifications();
      },
      py::return_value_policy::reference,
      "Tautomer matching specifications"
    )
    .def("add_to_hash", &UniqueMolecules::AddToHashes, "Add molecule to internal structures")
    .def("is_unique", &UniqueMolecules::IsUnique, "True if the molecule is unique")
    .def("report",
      [](const UniqueMolecules& u) {
        u.Report(std::cerr);
      },
      "Report"
    )
  ;
}
