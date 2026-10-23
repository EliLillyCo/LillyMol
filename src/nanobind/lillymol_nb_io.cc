#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindIo(nb::module_& m) {
  nb::enum_<FileType>(m, "FileType")
      .value("SMI", FILE_TYPE_SMI)
      .value("USMI", FILE_TYPE_USMI)
      .value("SDF", FILE_TYPE_SDF);

  nb::enum_<BondType>(m, "BondType")
      .value("SINGLE_BOND", BondType::kSingleBond)
      .value("DOUBLE_BOND", BondType::kDoubleBond)
      .value("TRIPLE_BOND", BondType::kTripleBond)
      .value("AROMATIC_BOND", BondType::kAromaticBond);

  nb::class_<ReaderOptions>(m, "ReaderOptions")
      .def("__init__", [](ReaderOptions* options, bool largest_fragment,
                          bool remove_chirality, bool remove_cis_trans_bonds,
                          bool remove_isotopes, bool keep_sdf_tags,
                          const std::string& sdf_identifier, bool sdf_tags_to_json,
                          bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid,
                          bool skip_invalid_valence, bool skip_non_organic,
                          bool skip_non_periodic_table_elements) {
             new (options) ReaderOptions();
             options->largest_fragment = largest_fragment;
             options->remove_chirality = remove_chirality;
             options->remove_cis_trans_bonds = remove_cis_trans_bonds;
             options->remove_isotopes = remove_isotopes;
             options->keep_sdf_tags = keep_sdf_tags;
             options->sdf_identifier = sdf_identifier;
             options->sdf_tags_to_json = sdf_tags_to_json;
             options->all_sdf_tags = all_sdf_tags;
             options->first_sdf_tag = first_sdf_tag;
             options->prepend_sdfid = prepend_sdfid;
             options->skip_invalid_valence = skip_invalid_valence;
             options->skip_non_organic = skip_non_organic;
             options->skip_non_periodic_table_elements = skip_non_periodic_table_elements;
           },
           nb::kw_only(),
           nb::arg("largest_fragment") = false,
           nb::arg("remove_chirality") = false,
           nb::arg("remove_cis_trans_bonds") = false,
           nb::arg("remove_isotopes") = false,
           nb::arg("keep_sdf_tags") = false,
           nb::arg("sdf_identifier") = "",
           nb::arg("sdf_tags_to_json") = false,
           nb::arg("all_sdf_tags") = false,
           nb::arg("first_sdf_tag") = false,
           nb::arg("prepend_sdfid") = true,
           nb::arg("skip_invalid_valence") = false,
           nb::arg("skip_non_organic") = false,
           nb::arg("skip_non_periodic_table_elements") = false)
      .def_rw("largest_fragment", &ReaderOptions::largest_fragment)
      .def_rw("remove_chirality", &ReaderOptions::remove_chirality)
      .def_rw("remove_cis_trans_bonds", &ReaderOptions::remove_cis_trans_bonds)
      .def_rw("remove_isotopes", &ReaderOptions::remove_isotopes)
      .def_rw("keep_sdf_tags", &ReaderOptions::keep_sdf_tags)
      .def_rw("sdf_identifier", &ReaderOptions::sdf_identifier)
      .def_rw("sdf_tags_to_json", &ReaderOptions::sdf_tags_to_json)
      .def_rw("all_sdf_tags", &ReaderOptions::all_sdf_tags)
      .def_rw("first_sdf_tag", &ReaderOptions::first_sdf_tag)
      .def_rw("prepend_sdfid", &ReaderOptions::prepend_sdfid)
      .def_rw("mdlquiet", &ReaderOptions::mdlquiet)
      .def_rw("skip_invalid_valence", &ReaderOptions::skip_invalid_valence)
      .def_rw("skip_non_organic", &ReaderOptions::skip_non_organic)
      .def_rw("skip_non_periodic_table_elements",
              &ReaderOptions::skip_non_periodic_table_elements);

  nb::class_<MoleculePreprocessing>(m, "MoleculePreprocessing")
      .def("__init__", [](MoleculePreprocessing* preprocessing, bool largest_fragment,
                          bool remove_chirality, bool remove_cis_trans_bonds,
                          bool remove_isotopes) {
             new (preprocessing) MoleculePreprocessing();
             preprocessing->set_reduce_to_largest_fragment(largest_fragment);
             preprocessing->set_remove_chirality(remove_chirality);
             preprocessing->set_remove_cis_trans_bonds(remove_cis_trans_bonds);
             preprocessing->set_remove_isotopes(remove_isotopes);
           },
           nb::kw_only(),
           nb::arg("largest_fragment") = false,
           nb::arg("remove_chirality") = false,
           nb::arg("remove_cis_trans_bonds") = false,
           nb::arg("remove_isotopes") = false)
      .def("active", &MoleculePreprocessing::active)
      .def("set_largest_fragment", &MoleculePreprocessing::set_reduce_to_largest_fragment)
      .def("set_reduce_to_largest_fragment", &MoleculePreprocessing::set_reduce_to_largest_fragment)
      .def("set_remove_chirality", &MoleculePreprocessing::set_remove_chirality)
      .def("set_remove_cis_trans_bonds", &MoleculePreprocessing::set_remove_cis_trans_bonds)
      .def("set_remove_isotopes", &MoleculePreprocessing::set_remove_isotopes)
      .def("process", [](MoleculePreprocessing& preprocessing, Molecule& mol) {
        return preprocessing.Process(mol);
      }, nb::arg("mol"))
      .def("process_copy", [](MoleculePreprocessing& preprocessing, const Molecule& mol) {
        Molecule result(mol);
        preprocessing.Process(result);
        return result;
      }, nb::arg("mol"));

  nb::class_<data_source_and_type<Molecule>>(m, "Reader")
      .def(nb::init<>())
      .def("__init__", [](data_source_and_type<Molecule>* reader, FileType file_type,
                          const std::string& fname) {
        IWString tmp(fname);
        new (reader) data_source_and_type<Molecule>(file_type, tmp);
      }, nb::arg("file_type"), nb::arg("fname"))
      .def("set_verbose", &data_source_and_type<Molecule>::set_verbose, nb::arg("verbose"))
      .def("set_skip_first", &data_source_and_type<Molecule>::set_skip_first,
           nb::arg("skip_first"))
      .def("set_do_only", &data_source_and_type<Molecule>::set_do_only,
           nb::arg("do_only"))
      .def("set_ignore_connection_table_errors",
           &data_source_and_type<Molecule>::set_connection_table_errors_allowed,
           nb::arg("errors_allowed"))
      .def("connection_table_errors_encountered",
           &data_source_and_type<Molecule>::connection_table_errors_encountered)
      .def("open", [](data_source_and_type<Molecule>& reader, const std::string& fname,
                      FileType file_type) {
        return static_cast<bool>(reader.do_open(fname.c_str(), file_type));
      }, nb::arg("fname"), nb::arg("file_type"))
      .def("open", [](data_source_and_type<Molecule>& reader, const std::string& fname) {
        IWString tmp(fname);
        FileType file_type = discern_file_type_from_name(tmp);
        if (file_type == FILE_TYPE_INVALID) {
          return false;
        }
        return static_cast<bool>(reader.do_open(fname.c_str(), file_type));
      }, nb::arg("fname"))
      .def("molecules_remaining", &data_source_and_type<Molecule>::molecules_remaining)
      .def("next", &ReaderNext)
      .def("molecules_read", &data_source_and_type<Molecule>::molecules_read)
      .def("close", &data_source_and_type<Molecule>::do_close)
      .def("__iter__", [](data_source_and_type<Molecule>& reader) -> data_source_and_type<Molecule>& {
        return reader;
      }, nb::rv_policy::reference_internal)
      .def("__next__", [](data_source_and_type<Molecule>& reader) {
        std::optional<Molecule> mol = ReaderNext(reader);
        if (!mol) {
          throw nb::stop_iteration();
        }
        return std::move(*mol);
      })
      .def("__enter__", [](data_source_and_type<Molecule>& reader) -> data_source_and_type<Molecule>& {
        return reader;
      }, nb::rv_policy::reference_internal)
      .def("__exit__", [](data_source_and_type<Molecule>& reader, nb::args) {
        reader.do_close();
      })
      .def("__repr__", [](const data_source_and_type<Molecule>& reader) {
        IWString result;
        result << "<Reader read " << reader.molecules_read() << '>';
        return result.AsString();
      });

  nb::class_<Molecule_Output_Object>(m, "Writer")
      .def(nb::init<>())
      .def("set_verbose", &Molecule_Output_Object::set_verbose, nb::arg("verbose"))
      .def("add_output_type", &Molecule_Output_Object::add_output_type,
           nb::arg("file_type"))
      .def("new_stem", [](Molecule_Output_Object& writer, const std::string& stem) {
        const_IWSubstring tmp(stem.data(), stem.size());
        return static_cast<bool>(writer.new_stem(tmp));
      }, nb::arg("stem"))
      .def("set_molecules_per_file", &Molecule_Output_Object::set_molecules_per_file,
           nb::arg("molecules_per_file"))
      .def("molecules_written", &Molecule_Output_Object::molecules_written)
      .def("write", [](Molecule_Output_Object& writer, Molecule& mol) {
        return static_cast<bool>(writer.write(mol));
      }, nb::arg("mol"))
      .def("flush", &Molecule_Output_Object::do_flush)
      .def("close", &Molecule_Output_Object::do_close);

  nb::class_<ReaderContext>(m, "MolReaderContext")
      .def(nb::init<const std::string&, FileType>(), nb::arg("fname"), nb::arg("file_type"))
      .def(nb::init<const std::string&>(), nb::arg("fname"))
      .def("__init__", [](ReaderContext* context, const std::string& fname, FileType file_type,
                          const ReaderOptions& options) {
             new (context) ReaderContext(fname, file_type);
             if (!context->ApplyOptions(options)) {
               throw std::invalid_argument("Invalid sdf_identifier regular expression");
             }
           }, nb::arg("fname"), nb::arg("file_type"), nb::kw_only(), nb::arg("options"))
      .def("__init__", [](ReaderContext* context, const std::string& fname,
                          const ReaderOptions& options) {
             new (context) ReaderContext(fname);
             if (!context->ApplyOptions(options)) {
               throw std::invalid_argument("Invalid sdf_identifier regular expression");
             }
           }, nb::arg("fname"), nb::kw_only(), nb::arg("options"))
      .def("__init__", [](ReaderContext* context, const std::string& fname, FileType file_type,
                          bool largest_fragment, bool remove_chirality,
                          bool remove_cis_trans_bonds, bool remove_isotopes,
                          bool keep_sdf_tags, const std::string& sdf_identifier,
                          bool sdf_tags_to_json, bool all_sdf_tags, bool first_sdf_tag,
                          bool prepend_sdfid, bool mdlquiet,
                          bool skip_invalid_valence, bool skip_non_organic,
                          bool skip_non_periodic_table_elements) {
             new (context) ReaderContext(fname, file_type);
             context->SetPreprocessing(largest_fragment, remove_chirality,
                                       remove_cis_trans_bonds, remove_isotopes);
             context->SetSkipOptions(skip_invalid_valence, skip_non_organic,
                                     skip_non_periodic_table_elements);
             if (context->reader) {
               context->reader->mutable_molecule_read_options()
                   .mdl_file_supporting_material()
                   .set_read_extra_text_info(keep_sdf_tags);
             }
             if (!context->SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                                         first_sdf_tag, prepend_sdfid, mdlquiet)) {
               throw std::invalid_argument("Invalid sdf_identifier regular expression");
             }
           },
           nb::arg("fname"), nb::arg("file_type"), nb::kw_only(),
           nb::arg("largest_fragment") = false,
           nb::arg("remove_chirality") = false,
           nb::arg("remove_cis_trans_bonds") = false,
           nb::arg("remove_isotopes") = false,
           nb::arg("keep_sdf_tags") = false,
           nb::arg("sdf_identifier") = "",
           nb::arg("sdf_tags_to_json") = false,
           nb::arg("all_sdf_tags") = false,
           nb::arg("first_sdf_tag") = false,
           nb::arg("prepend_sdfid") = true,
           nb::arg("mdlquiet") = false,
           nb::arg("skip_invalid_valence") = false,
           nb::arg("skip_non_organic") = false,
           nb::arg("skip_non_periodic_table_elements") = false)
      .def("__init__", [](ReaderContext* context, const std::string& fname,
                          bool largest_fragment, bool remove_chirality,
                          bool remove_cis_trans_bonds, bool remove_isotopes,
                          bool keep_sdf_tags, const std::string& sdf_identifier,
                          bool sdf_tags_to_json, bool all_sdf_tags, bool first_sdf_tag,
                          bool prepend_sdfid, bool mdlquiet,
                          bool skip_invalid_valence, bool skip_non_organic,
                          bool skip_non_periodic_table_elements) {
             new (context) ReaderContext(fname);
             context->SetPreprocessing(largest_fragment, remove_chirality,
                                       remove_cis_trans_bonds, remove_isotopes);
             context->SetSkipOptions(skip_invalid_valence, skip_non_organic,
                                     skip_non_periodic_table_elements);
             if (context->reader) {
               context->reader->mutable_molecule_read_options()
                   .mdl_file_supporting_material()
                   .set_read_extra_text_info(keep_sdf_tags);
             }
             if (!context->SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                                         first_sdf_tag, prepend_sdfid, mdlquiet)) {
               throw std::invalid_argument("Invalid sdf_identifier regular expression");
             }
           },
           nb::arg("fname"), nb::kw_only(),
           nb::arg("largest_fragment") = false,
           nb::arg("remove_chirality") = false,
           nb::arg("remove_cis_trans_bonds") = false,
           nb::arg("remove_isotopes") = false,
           nb::arg("keep_sdf_tags") = false,
           nb::arg("sdf_identifier") = "",
           nb::arg("sdf_tags_to_json") = false,
           nb::arg("all_sdf_tags") = false,
           nb::arg("first_sdf_tag") = false,
           nb::arg("prepend_sdfid") = true,
           nb::arg("mdlquiet") = false,
           nb::arg("skip_invalid_valence") = false,
           nb::arg("skip_non_organic") = false,
           nb::arg("skip_non_periodic_table_elements") = false)
      .def("__enter__", [](ReaderContext& context) -> ReaderContext& {
        return context;
      }, nb::rv_policy::reference_internal)
      .def("__iter__", [](ReaderContext& context) -> ReaderContext& {
        return context;
      }, nb::rv_policy::reference_internal)
      .def("__next__", [](ReaderContext& context) {
        std::optional<Molecule> mol = context.Next();
        if (!mol) {
          throw nb::stop_iteration();
        }
        return std::move(*mol);
      })
      .def("next", &ReaderContext::Next)
      .def("apply_options", [](ReaderContext& context, const ReaderOptions& options) {
        if (!context.ApplyOptions(options)) {
          throw std::invalid_argument("Invalid sdf_identifier regular expression");
        }
      }, nb::arg("options"))
      .def("set_preprocessing", [](ReaderContext& context, bool largest_fragment,
                                   bool remove_chirality, bool remove_cis_trans_bonds,
                                   bool remove_isotopes) {
        context.SetPreprocessing(largest_fragment, remove_chirality,
                                 remove_cis_trans_bonds, remove_isotopes);
      }, nb::kw_only(),
         nb::arg("largest_fragment") = false,
         nb::arg("remove_chirality") = false,
         nb::arg("remove_cis_trans_bonds") = false,
         nb::arg("remove_isotopes") = false)
      .def("preprocessing_active", [](const ReaderContext& context) {
        return context.preprocessing.active();
      })
      .def("set_skip_options", [](ReaderContext& context, bool skip_invalid_valence,
                                   bool skip_non_organic,
                                   bool skip_non_periodic_table_elements) {
        context.SetSkipOptions(skip_invalid_valence, skip_non_organic,
                               skip_non_periodic_table_elements);
      }, nb::kw_only(),
         nb::arg("skip_invalid_valence") = false,
         nb::arg("skip_non_organic") = false,
         nb::arg("skip_non_periodic_table_elements") = false)
      .def("skipped_invalid_valence", [](const ReaderContext& context) {
        return context.skipped_invalid_valence;
      })
      .def("skipped_non_organic", [](const ReaderContext& context) {
        return context.skipped_non_organic;
      })
      .def("skipped_non_periodic_table_elements", [](const ReaderContext& context) {
        return context.skipped_non_periodic_table_elements;
      })
      .def("molecules_skipped", [](const ReaderContext& context) {
        return context.skipped_invalid_valence + context.skipped_non_organic +
               context.skipped_non_periodic_table_elements;
      })
      .def("set_sdf_options", [](ReaderContext& context, const std::string& sdf_identifier,
                                 bool sdf_tags_to_json, bool all_sdf_tags,
                                 bool first_sdf_tag, bool prepend_sdfid,
                                 bool mdlquiet) {
        if (!context.SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                                   first_sdf_tag, prepend_sdfid, mdlquiet)) {
          throw std::invalid_argument("Invalid sdf_identifier regular expression");
        }
      }, nb::kw_only(),
         nb::arg("sdf_identifier") = "",
         nb::arg("sdf_tags_to_json") = false,
         nb::arg("all_sdf_tags") = false,
         nb::arg("first_sdf_tag") = false,
         nb::arg("prepend_sdfid") = true,
         nb::arg("mdlquiet") = false)
      .def("reset_sdf_options", &ReaderContext::ResetSdfOptions)
      .def("set_ignore_connection_table_errors", [](ReaderContext& context, int value) {
        if (context.reader) {
          context.reader->set_connection_table_errors_allowed(value);
        }
      }, nb::arg("errors_allowed"))
      .def("connection_table_errors_encountered", [](const ReaderContext& context) {
        return context.reader ? context.reader->connection_table_errors_encountered() : 0;
      })
      .def("molecules_remaining", [](ReaderContext& context) -> uint64_t {
        return context.reader ? context.reader->records_remaining() : 0;
      })
      .def("molecules_read", [](const ReaderContext& context) -> uint64_t {
        return context.reader ? context.reader->molecules_read() : 0;
      })
      .def("close", [](ReaderContext& context) {
        if (context.reader) {
          context.reader->do_close();
        }
        context.ResetSdfOptions();
        context.RestoreElementCreationFlags();
      })
      .def("__exit__", [](ReaderContext& context, nb::args) {
        if (context.reader) {
          context.reader->do_close();
        }
        context.ResetSdfOptions();
        context.RestoreElementCreationFlags();
      });

  m.attr("ReaderContext") = m.attr("MolReaderContext");

  nb::class_<ContextWriter>(m, "MolWriterContext")
      .def(nb::init<const std::string&, FileType>(), nb::arg("fname"), nb::arg("file_type"))
      .def(nb::init<const std::string&>(), nb::arg("fname"))
      .def("__enter__", [](ContextWriter& context) -> ContextWriter& {
        return context;
      }, nb::rv_policy::reference_internal)
      .def("__exit__", [](ContextWriter& context, nb::args) {
        if (context.writer) {
          context.writer->do_close();
        }
      })
      .def("write", [](ContextWriter& context, Molecule& mol) {
        return context.writer && context.writer->write(mol);
      }, nb::arg("mol"))
      .def("set_molecules_per_file", [](ContextWriter& context, int value) {
        return context.writer ? context.writer->set_molecules_per_file(value) : 0;
      }, nb::arg("molecules_per_file"))
      .def("set_name_token_for_file_name", [](ContextWriter& context, int token) {
        if (context.writer) {
          context.writer->set_name_token_for_file_name(token);
        }
      }, nb::arg("token"));

  m.attr("ContextWriter") = m.attr("MolWriterContext");

  m.def("slurp",
        [](const std::string& fname) -> std::optional<std::vector<Molecule>> {
          IWString tmp(fname);
          FileType file_type = discern_file_type_from_name(tmp);
          if (file_type == FILE_TYPE_INVALID) {
            file_type = FILE_TYPE_SMI;
          }

          data_source_and_type<Molecule> input(file_type, tmp);
          if (!input.good()) {
            std::cerr << "slurp:cannot open '" << fname << "'\n";
            return std::nullopt;
          }

          const uint32_t number_molecules = input.molecules_remaining();
          std::vector<Molecule> result(number_molecules);
          for (uint32_t i = 0; i < number_molecules; ++i) {
            if (!input.next_molecule(result[i])) {
              return std::nullopt;
            }
          }

          return result;
        },
        nb::arg("fname"), "Read all molecules from fname");

  m.def("set_sdf_identifier",
        [](const std::string& sdf_identifier) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          const_IWSubstring tmp(sdf_identifier);
          return static_cast<bool>(mdlfos->set_sdf_identifier(tmp));
        },
        nb::arg("sdf_identifier"));
  m.def("set_ignore_bad_m",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_ignore_unrecognised_mdl_m_records(value);
        },
        nb::arg("value"));
  m.def("set_mdlquiet",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_report_unrecognised_records(!value);
        },
        nb::arg("value"));
  m.def("set_prepend_sdfid",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_prepend_sdfid(value);
        },
        nb::arg("value"));
  m.def("set_allsdfid",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_fetch_all_sdf_identifiers(value);
        },
        nb::arg("value"));
  m.def("set_sdf_tags_to_json",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_sdf_tags_to_json(value);
        },
        nb::arg("value"));
  m.def("set_firstsdftag",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_take_first_tag_as_name(value);
        },
        nb::arg("value"));
  m.def("set_allow_deuterium",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_allow_deuterium(value);
        },
        nb::arg("value"));
  m.def("set_allow_tritium",
        [](bool value) {
          MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
          mdlfos->set_allow_tritium(value);
        },
        nb::arg("value"));


}

}  // namespace lillymol_nb
