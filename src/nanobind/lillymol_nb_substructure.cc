#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindSubstructure(nb::module_& m) {
  nb::class_<Substructure_Results>(m, "SubstructureResults")
      .def(nb::init<>())
      .def("number_embeddings", &Substructure_Results::number_embeddings)
      .def("hits_found", &Substructure_Results::hits_found)
      .def("each_embedding_set_vector",
           [](const Substructure_Results& sresults, int matoms, int value) {
             return sresults.EachEmbeddingSetVector(matoms, value);
           },
           nb::arg("matoms"), nb::arg("value"))
      .def("max_query_atoms_matched_in_search",
           &Substructure_Results::max_query_atoms_matched_in_search)
      .def("reset", &Substructure_Results::reset)
      .def("embeddings", &SubstructureResultsAsVectors)
      .def("embeddings_as_set_of_atoms", &SubstructureResultsAsSetOfAtoms)
      .def("__len__", &Substructure_Results::number_embeddings)
      .def("__iter__",
           [](const Substructure_Results& sresults) {
             const SetOfEmbeddings& embeddings = sresults.embeddings();
             return nb::make_iterator<nb::rv_policy::reference_internal>(
                 nb::type<Substructure_Results>(), "EmbeddingIterator",
                 embeddings.begin(), embeddings.end());
           },
           nb::keep_alive<0, 1>())
      .def("__getitem__",
           [](const Substructure_Results& sresults, int index) {
             const Set_of_Atoms* embedding = sresults.embedding(index);
             if (embedding == nullptr) {
               throw std::out_of_range("SubstructureResults embedding index out of range");
             }
             return Set_of_Atoms(*embedding);
           });

  nb::class_<Substructure_Query>(m, "SubstructureQuery")
      .def(nb::init<>())
      .def("build_from_smarts", &Substructure_Query::CreateFromSmarts, nb::arg("smarts"))
      .def("build_from_molecule",
           [](Substructure_Query& query, Molecule& mol) {
             Molecule_to_Query_Specifications mqs;
             mqs.set_make_embedding(1);
             return static_cast<bool>(query.create_from_molecule(mol, mqs));
           },
           nb::arg("mol"))
      .def("number_elements", [](const Substructure_Query& query) { return query.number_elements(); })
      .def("active", &Substructure_Query::active)
      .def("set_only_keep_matches_in_largest_fragment",
           &Substructure_Query::set_only_keep_matches_in_largest_fragment)
      .def("set_embeddings_do_not_overlap", &Substructure_Query::set_embeddings_do_not_overlap)
      .def("set_find_one_embedding_per_atom", &Substructure_Query::set_find_one_embedding_per_atom)
      .def("set_find_unique_embeddings_only", &Substructure_Query::set_find_unique_embeddings_only)
      .def("set_max_matches_to_find", &Substructure_Query::set_max_matches_to_find)
      .def("set_perceive_symmetry_equivalent_matches",
           &Substructure_Query::set_perceive_symmetry_equivalent_matches)
      .def("set_min_atoms_to_match", &Substructure_Query::set_min_atoms_to_match)
      .def("set_max_atoms_to_match", &Substructure_Query::set_max_atoms_to_match)
      .def("max_query_atoms_matched_in_search",
           &Substructure_Query::max_query_atoms_matched_in_search)
      .def("substructure_search",
           [](Substructure_Query& query, Molecule& mol) {
             nb::gil_scoped_release release;
             return query.substructure_search(&mol);
           },
           nb::arg("mol"))
      .def("substructure_search",
           [](Substructure_Query& query, Molecule& mol, Substructure_Results& sresults) {
             nb::gil_scoped_release release;
             return query.substructure_search(mol, sresults);
           },
           nb::arg("mol"), nb::arg("sresults"))
      .def("substructure_search_matches",
           [](Substructure_Query& query, Molecule& mol) -> std::optional<std::vector<Set_of_Atoms>> {
             Substructure_Results sresults;
             uint32_t matches;
             {
               nb::gil_scoped_release release;
               matches = query.substructure_search(&mol, sresults);
             }
             if (matches == 0) {
               return std::nullopt;
             }
             return SubstructureResultsAsSetOfAtoms(sresults);
           },
           nb::arg("mol"))
      .def("substructure_search_match_lists",
           [](Substructure_Query& query, Molecule& mol) {
             return SubstructureSearchMatches(query, mol);
           },
           nb::arg("mol"))
      .def("read_proto",
           [](Substructure_Query& query, const std::string& fname) {
             return static_cast<bool>(query.ReadProto(IWString(fname)));
           },
           nb::arg("fname"))
      .def("construct_from_proto",
           [](Substructure_Query& query, nb::object proto) {
             if (!nb::hasattr(proto, "SerializeToString")) {
               throw nb::type_error("construct_from_proto requires a Python protobuf object");
             }
             nb::object serialized_obj = proto.attr("SerializeToString")();
             PyObject* bytes = serialized_obj.ptr();
             if (!PyBytes_Check(bytes)) {
               throw nb::type_error("SerializeToString() did not return bytes");
             }
             char* data = nullptr;
             Py_ssize_t size = 0;
             if (PyBytes_AsStringAndSize(bytes, &data, &size) != 0) {
               throw nb::python_error();
             }
             SubstructureSearch::SubstructureQuery cproto;
             if (!cproto.ParseFromArray(data, static_cast<int>(size))) {
               throw std::invalid_argument("construct_from_proto could not parse serialized SubstructureQuery proto");
             }
             return static_cast<bool>(query.ConstructFromProto(cproto));
           },
           nb::arg("proto"),
           "Construct from a Python SubstructureQuery protobuf object. The nanobind boundary serializes the proto before constructing the C++ query.")
      .def("read_msi",
           [](Substructure_Query& query, const std::string& fname) {
             return static_cast<bool>(query.read(IWString(fname)));
           },
           nb::arg("fname"))
      .def("__repr__",
           [](const Substructure_Query& query) {
             IWString result;
             result << "<SubstructureQuery " << query.comment() << '>';
             return result.AsString();
           });


}

}  // namespace lillymol_nb
