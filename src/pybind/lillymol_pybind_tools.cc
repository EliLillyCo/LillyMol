// Bindings for selected tools

#include <optional>
#include <string>
#include <unordered_map>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Tools/dicer_api.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/unique_molecules_api.h"
#include "Molecule_Tools_Bdb/iwecfp_database_lookup_lib.h"
#include "Molecule_Tools_Bdb/selimsteg.h"
#include "Molecule_Tools_Bdb/structure_database.h"

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
    .def("set_include_chiral_info", &UniqueMolecules::set_include_chiral_info, "Control whether chirality is considered")
    .def("set_include_cis_trans_bonding_info", &UniqueMolecules::set_include_cis_trans_bonding_info, "Control whether cis trans bonds are considered")
    .def("set_strip_to_largest_fragment", &UniqueMolecules::set_strip_to_largest_fragment, "Strip to largest fragment")
    .def("set_consider_isotopes", &UniqueMolecules::set_consider_isotopes, "Control whether isotopes are cosidered")
    .def("set_constant_isotope", &UniqueMolecules::set_constant_isotope, "Convert all isotopes to a constant")
    .def("set_standardize_molecules", &UniqueMolecules::set_standardize_molecules, "Control whether or not molecules are standardised")
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
//#endif

  py::class_<selimsteg::Selimsteg>(m, "Selimsteg")
    .def(py::init<>())
    .def("open_database", &selimsteg::Selimsteg::OpenDatabase, "Open a BerkeleyDB datbase with Id->smiles mappints")
    .def("get_smiles", &selimsteg::Selimsteg::Lookup, "Fetch the smiles for an identifier")
    .def("get_molecule", &selimsteg::Selimsteg::GetMolecule, "Fetch a Molecule for an identifier")
    .def("get_molecules", &selimsteg::Selimsteg::GetMolecules, "Fetch a list of Molecules for list of identifiers")
  ;

  py::enum_<structure_database::Lookup>(m, "LookupParams", py::arithmetic())
    .value("EXACT", structure_database::kExact)
    .value("STRIP", structure_database::kStrip)
    .value("NOCHIRAL", structure_database::kNoChiral)
    .value("GRAPH", structure_database::kGraph)
    .value("NOSTD", structure_database::kNoStandardise)
    .export_values();
  ;
  m.def("value", [](structure_database::Lookup e) { return static_cast<int>(e); });

  py::class_<structure_database::StructureDatabase>(m, "StructureDatabase")
    .def(py::init<>())
    .def("open_read",
      [](structure_database::StructureDatabase& db, const std::string& dbname)->bool {
        IWString tmp(dbname);
        return db.OpenForReading(tmp);
      },
      "Open a structure database for reading"
    )
    .def("lookup",
      [](structure_database::StructureDatabase& db, Molecule& m)->std::optional<std::string> {
        IWString tmp;
        const uint32_t mask = structure_database::Lookup::kExact;
        int rc = db.Lookup(m, mask, tmp);
        if (rc == 0) {
          return std::nullopt;
        }
        std::string result = tmp.AsString();
        return result;
      },
      ""
    )
    .def("lookup",
      [](structure_database::StructureDatabase& db, Molecule& m, const uint32_t params)-> std::optional<std::string> {
        IWString tmp;
        int rc = db.Lookup(m, params, tmp);
        if (rc == 0) {
          return std::nullopt;
        }
        std::string result = tmp.AsString();
        return result;
      },
      "Lookup molecule in database, returning id's of equivalent molecules"
    )
  ;
  py::class_<dicer_api::Dicer>(m, "Dicer")
    .def(py::init<>())
    .def("set_max_bonds_to_break", &dicer_api::Dicer::set_max_bonds_to_break, "set max bonds to break")
    .def("set_min_fragment_size", &dicer_api::Dicer::set_min_fragment_size, "set min fragment size")
    .def("set_max_fragment_size", &dicer_api::Dicer::set_max_fragment_size, "set max fragment size")
    .def("set_break_cc_bonds", &dicer_api::Dicer::set_break_cc_bonds, "set True of C-C bonds are to be broken")
    .def("set_label_join_points", &dicer_api::Dicer::set_label_join_points, "Isotope for join points")
    .def("set_accumulate_global_fragment_count", &dicer_api::Dicer::set_accumulate_global_fragment_count, "Set to True to accumulate fragments")
    .def("get_global_fragment_count", &dicer_api::Dicer::global_fragment_count, "Fetch the global fragment count")
    .def("set_perceive_symmetry_equivalent_matches", &dicer_api::Dicer::set_perceive_symmetry_equivalent_matches, "Set 0 to NOT perceive symmetry equivalent matches")
    .def("set_determine_fragment_counts", &dicer_api::Dicer::set_determine_fragment_counts, "If True determine the number of times a fragment appears in each starting molecule")
    .def("set_work_like_recap", &dicer_api::Dicer::set_work_like_recap, "Work like Recap - no recursion")
    .def("add_bond_break_smarts",
      [](dicer_api::Dicer& dicer, const std::string& smarts) {
        return dicer.AddBreakBondSmarts(smarts);
      },
      "Add bond breaking smarts"
    )
    .def("add_bond_break_query",
      [](dicer_api::Dicer& dicer, const std::string& fname) {
        return dicer.AddBreakBondQuery(fname);
      },
      "Add bond breaking query from textproto file"
    )
    .def("dice",
      [](dicer_api::Dicer& dicer, Molecule& m)->std::unordered_map<std::string, uint32_t> {
        std::unordered_map<std::string, uint32_t> result;
        dicer.Dice(m, result);
        return result;
      },
      "dice the molecule and return fragments and counts"
    )
  ;
}
