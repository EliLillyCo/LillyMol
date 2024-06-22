#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Tools_Bdb/iwecfp_database_lookup_lib.h"
#include "Molecule_Tools_Bdb/selimsteg.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_tools_bdb, m) 
{
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
    .def("slurp",
      [](iwecfp_database_lookup::SP_Set_of_Databases& databases, uint32_t min_examples)->bool {
        return databases.slurp(min_examples);
      },
      "Slurp database entries with min_examples or more examples to memory"
    )
    .def("__repr__",
      [](const iwecfp_database_lookup::Set_of_Databases& databases){
        IWString result;
        result << "Synethetic Precedent database with " << databases.number_databases() << " databases";
        return result.AsString();
      }
    )
  ;

  py::class_<selimsteg::Selimsteg>(m, "Selimsteg")
    .def(py::init<>())
    .def("open_database", &selimsteg::Selimsteg::OpenDatabase, "Open a BerkeleyDB datbase with Id->smiles mappints")
    .def("get_smiles", &selimsteg::Selimsteg::Lookup, "Fetch the smiles for an identifier")
    .def("get_molecule", &selimsteg::Selimsteg::GetMolecule, "Fetch a Molecule for an identifier")
    .def("get_molecules", &selimsteg::Selimsteg::GetMolecules, "Fetch a list of Molecules for list of identifiers")
  ;
}
