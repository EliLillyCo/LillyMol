// Bindings for selected tools

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Tools/unique_molecules_api.h"
#include "Molecule_Tools_Bdb/iwecfp_database_lookup_lib.h"

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

//#ifdef THIS_WILL_REQUIRE_BERKELEY_DB_DYNAMIC_LIBRARIES
  py::class_<iwecfp_database_lookup::SP_Set_of_Databases>(m, "SyntheticPrecedentDatabases")
    .def(py::init<>())
    .def("add_database",
      [](iwecfp_database_lookup::SP_Set_of_Databases& databases, const std::string& dbname)->bool {
        return databases.AddDatabase(dbname);
      },
      "Add an existing synethetic precedent database"
    )
    .def("set_max_radius",
      [](iwecfp_database_lookup::SP_Set_of_Databases& databases, int max_radius)->bool {
        return databases.set_max_radius(max_radius);
      },
      "Set max radius for shells"
    )
    .def("per_shell_data",
      [](iwecfp_database_lookup::SP_Set_of_Databases& databases, Molecule& m)->std::vector<int> {
        std::vector<int> result(3);

        // Need to do something better here, but does it ever fail?
        if (! databases.PerShellData(m, result)) {
          std::cerr << "Failed\n";
        }
        return result;
      },
      "Report lowest bit prevalence at each radius"
    )
    .def("__repr__",
      [](const iwecfp_database_lookup::Set_of_Databases& databases){
        IWString result;
        result << "Synethetic Precedent database with " << databases.number_databases() << " databases";
        return result.AsString();
      }
    )
  ;
//#endif
}
