#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"
#include "Molecule_Lib/output.h"

namespace py = pybind11;

using molecule_processing::MoleculePreprocessing;

struct ReaderOptions {
  bool largest_fragment = false;
  bool remove_chirality = false;
  bool remove_cis_trans_bonds = false;
  bool remove_isotopes = false;
  bool keep_sdf_tags = false;

  std::string sdf_identifier;
  bool sdf_tags_to_json = false;
  bool all_sdf_tags = false;
  bool first_sdf_tag = false;
  bool prepend_sdfid = true;
};

std::unique_ptr<data_source_and_type<Molecule>>
MakeReader(const std::string& fname,
           FileType file_type) {
  std::unique_ptr<data_source_and_type<Molecule>> result = 
    std::make_unique<data_source_and_type<Molecule>>();

  IWString tmp(fname);
  if (! result->do_open(tmp, file_type)) {
    std::cerr << "Cannot open '" << fname << "'\n";
    return result;
  }

  return result;
}

struct ReaderContext {
  public:
    std::unique_ptr<data_source_and_type<Molecule>> _reader;
    MoleculePreprocessing _preprocessing;

  public:
    ReaderContext(const std::string& fname, FileType file_type);
    ReaderContext(const std::string& fname);
    ~ReaderContext();

    bool ApplyOptions(const ReaderOptions& options);
    void SetPreprocessing(bool reduce_to_largest_fragment, bool remove_chirality,
                          bool remove_cis_trans_bonds, bool remove_isotopes);
    bool SetSdfOptions(const std::string& sdf_identifier, bool sdf_tags_to_json,
                       bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid);
    void ResetSdfOptions();
    std::optional<Molecule> Next();
};

ReaderContext::ReaderContext(const std::string& fname, FileType file_type) {
  IWString tmp(fname);
  _reader = std::make_unique<data_source_and_type<Molecule>>(file_type, tmp);
}

ReaderContext::ReaderContext(const std::string& fname) {
  IWString tmp(fname);
  FileType file_type = discern_file_type_from_name(fname);
  if (file_type == FILE_TYPE_INVALID) {
    std::cerr << "ReaderContext::ReaderContext:unrecognised type '" << fname << "'\n";
    return;
  }

  _reader = std::make_unique<data_source_and_type<Molecule>>(file_type, tmp);
}

ReaderContext::~ReaderContext() {
  ResetSdfOptions();
}

bool
ReaderContext::ApplyOptions(const ReaderOptions& options) {
  SetPreprocessing(options.largest_fragment, options.remove_chirality,
                   options.remove_cis_trans_bonds, options.remove_isotopes);
  if (_reader) {
    _reader->mutable_molecule_read_options()
        .mdl_file_supporting_material()
        .set_read_extra_text_info(options.keep_sdf_tags);
  }
  return SetSdfOptions(options.sdf_identifier, options.sdf_tags_to_json,
                       options.all_sdf_tags, options.first_sdf_tag,
                       options.prepend_sdfid);
}

void
ReaderContext::SetPreprocessing(bool reduce_to_largest_fragment, bool remove_chirality,
                                bool remove_cis_trans_bonds, bool remove_isotopes) {
  _preprocessing.set_reduce_to_largest_fragment(reduce_to_largest_fragment);
  _preprocessing.set_remove_chirality(remove_chirality);
  _preprocessing.set_remove_cis_trans_bonds(remove_cis_trans_bonds);
  _preprocessing.set_remove_isotopes(remove_isotopes);
}

bool
ReaderContext::SetSdfOptions(const std::string& sdf_identifier, bool sdf_tags_to_json,
                             bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid) {
  if (sdf_identifier.empty() && !sdf_tags_to_json && !all_sdf_tags && !first_sdf_tag &&
      prepend_sdfid) {
    return true;
  }

  if (! _reader) {
    return false;
  }

  MDL_File_Supporting_Material& mdlfos =
      _reader->mutable_molecule_read_options().mdl_file_supporting_material();
  const_IWSubstring tmp(sdf_identifier.data(), sdf_identifier.length());
  if (! mdlfos.set_sdf_identifier(tmp)) {
    return false;
  }

  mdlfos.set_sdf_tags_to_json(sdf_tags_to_json);
  mdlfos.set_fetch_all_sdf_identifiers(all_sdf_tags);
  mdlfos.set_take_first_tag_as_name(first_sdf_tag);
  mdlfos.set_prepend_sdfid(prepend_sdfid);

  return true;
}

void
ReaderContext::ResetSdfOptions() {
  if (_reader) {
    _reader->reset_molecule_read_options();
  }
}

std::optional<Molecule>
ReaderContext::Next() {
  if (! _reader) {
    return std::nullopt;
  }

  std::unique_ptr<Molecule> m(_reader->next_molecule());
  if (! m) {
    return std::nullopt;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(*m);
  }

  return std::move(*m);
}

struct ContextWriter {
  public:
    std::unique_ptr<Molecule_Output_Object> _writer;

  private:
    void CommonFileOpen(const std::string& fname, FileType file_type);

  public:
    ContextWriter(const std::string& fname);
    ContextWriter(const std::string& fname, FileType file_type);
    // ContextWriter(const std::string& fname, const std::vector<FileType>& file_types);
};

ContextWriter::ContextWriter(const std::string& fname) {
  CommonFileOpen(fname, FILE_TYPE_SMI);
}

ContextWriter::ContextWriter(const std::string& fname, FileType file_type) {
  CommonFileOpen(fname, file_type);
}

void
ContextWriter::CommonFileOpen(const std::string& fname, FileType file_type) {
  _writer = std::make_unique<Molecule_Output_Object>();

#ifdef EERRQ
  switch (file_type) {
    case FileType::SMI:

  }
#endif

  if (! _writer->add_output_type(file_type)) {
    std::cerr << "WriterContext::CommonFileOpen:unrecognised type " << file_type << '\n';
    _writer.reset(nullptr);
    return;
  }

  IWString tmp(fname);
  if (! _writer->new_stem(tmp)) {
    std::cerr << "WriterContext::CommonFileOpen:cannot open '" << tmp << "'\n";
    _writer.reset(nullptr);
    return;
  }
}


PYBIND11_MODULE(lillymol_io, io)
{
  struct ReaderIterator {
    data_source_and_type<Molecule>& input;
    py::object ref;
    ReaderIterator(data_source_and_type<Molecule>& reader, py::object ref) :
                input(reader), ref(ref) {
    }
    Molecule* next() {
      Molecule* m = input.next_molecule();
      if (m == nullptr) {
        throw py::stop_iteration();
      }
      return m;
    }
  };

  py::class_<ReaderOptions>(io, "ReaderOptions")
    .def(py::init([](bool largest_fragment, bool remove_chirality,
                     bool remove_cis_trans_bonds, bool remove_isotopes,
                     bool keep_sdf_tags, const std::string& sdf_identifier,
                     bool sdf_tags_to_json, bool all_sdf_tags, bool first_sdf_tag,
                     bool prepend_sdfid) {
        ReaderOptions options;
        options.largest_fragment = largest_fragment;
        options.remove_chirality = remove_chirality;
        options.remove_cis_trans_bonds = remove_cis_trans_bonds;
        options.remove_isotopes = remove_isotopes;
        options.keep_sdf_tags = keep_sdf_tags;
        options.sdf_identifier = sdf_identifier;
        options.sdf_tags_to_json = sdf_tags_to_json;
        options.all_sdf_tags = all_sdf_tags;
        options.first_sdf_tag = first_sdf_tag;
        options.prepend_sdfid = prepend_sdfid;
        return options;
      }),
      py::kw_only(),
      py::arg("largest_fragment") = false,
      py::arg("remove_chirality") = false,
      py::arg("remove_cis_trans_bonds") = false,
      py::arg("remove_isotopes") = false,
      py::arg("keep_sdf_tags") = false,
      py::arg("sdf_identifier") = "",
      py::arg("sdf_tags_to_json") = false,
      py::arg("all_sdf_tags") = false,
      py::arg("first_sdf_tag") = false,
      py::arg("prepend_sdfid") = true)
    .def_readwrite("largest_fragment", &ReaderOptions::largest_fragment)
    .def_readwrite("remove_chirality", &ReaderOptions::remove_chirality)
    .def_readwrite("remove_cis_trans_bonds", &ReaderOptions::remove_cis_trans_bonds)
    .def_readwrite("remove_isotopes", &ReaderOptions::remove_isotopes)
    .def_readwrite("keep_sdf_tags", &ReaderOptions::keep_sdf_tags)
    .def_readwrite("sdf_identifier", &ReaderOptions::sdf_identifier)
    .def_readwrite("sdf_tags_to_json", &ReaderOptions::sdf_tags_to_json)
    .def_readwrite("all_sdf_tags", &ReaderOptions::all_sdf_tags)
    .def_readwrite("first_sdf_tag", &ReaderOptions::first_sdf_tag)
    .def_readwrite("prepend_sdfid", &ReaderOptions::prepend_sdfid)
  ;

  py::class_<MoleculePreprocessing>(io, "MoleculePreprocessing")
    .def(py::init([](bool largest_fragment, bool remove_chirality,
                     bool remove_cis_trans_bonds, bool remove_isotopes) {
        MoleculePreprocessing* result = new MoleculePreprocessing;
        result->set_reduce_to_largest_fragment(largest_fragment);
        result->set_remove_chirality(remove_chirality);
        result->set_remove_cis_trans_bonds(remove_cis_trans_bonds);
        result->set_remove_isotopes(remove_isotopes);
        return result;
      }),
      py::kw_only(),
      py::arg("largest_fragment") = false,
      py::arg("remove_chirality") = false,
      py::arg("remove_cis_trans_bonds") = false,
      py::arg("remove_isotopes") = false)
    .def("active", &MoleculePreprocessing::active)
    .def("set_largest_fragment", &MoleculePreprocessing::set_reduce_to_largest_fragment)
    .def("set_reduce_to_largest_fragment", &MoleculePreprocessing::set_reduce_to_largest_fragment)
    .def("set_remove_chirality", &MoleculePreprocessing::set_remove_chirality)
    .def("set_remove_cis_trans_bonds", &MoleculePreprocessing::set_remove_cis_trans_bonds)
    .def("set_remove_isotopes", &MoleculePreprocessing::set_remove_isotopes)
    .def("process", [](MoleculePreprocessing& preprocessing, Molecule& m)->int {
        return preprocessing.Process(m);
      })
    .def("process_copy", [](MoleculePreprocessing& preprocessing, const Molecule& m)->Molecule {
        Molecule result(m);
        preprocessing.Process(result);
        return result;
      })
  ;

  py::class_<ReaderIterator>(io, "Iterator")
    .def("__iter__", [](ReaderIterator &it) -> ReaderIterator& { return it; })
    .def("__next__", &ReaderIterator::next)
  ;

  py::class_<data_source_and_type<Molecule>>(io, "Reader")
    .def(py::init<>())
    .def(py::init<FileType, const std::string&>())
    .def("set_verbose",
      [](data_source_and_type<Molecule>& inp, int v) {
        return inp.set_verbose(v);
      },
      "verbosity"
    )
    .def("set_skip_first",
      [](data_source_and_type<Molecule>& inp, int s) {
        inp.set_skip_first(s);
      },
      "skip first n records"
    )
    .def("set_do_only",
      [](data_source_and_type<Molecule>& inp, int s) {
        inp.set_do_only(s);
      },
      "process only first n records"
    )
    .def("set_ignore_connection_table_errors",
      [](data_source_and_type<Molecule>& inp, int s) {
        inp.set_connection_table_errors_allowed(s);
      },
      "Number connection table errors allowed"
    )
    .def("connection_table_errors_encountered",
      [](data_source_and_type<Molecule>& me)->int {
        return me.connection_table_errors_encountered();
      }
    )
    .def("open",
      [](data_source_and_type<Molecule>& inp, const std::string& fname, FileType file_type)->bool{
        return inp.do_open(fname.c_str(), file_type);
      },
      "open file with file type"
    )
    .def("open",
      [](data_source_and_type<Molecule>& inp, const std::string& fname)->bool {
        IWString tmp(fname);
        FileType file_type = discern_file_type_from_name(tmp);
        if (file_type == FILE_TYPE_INVALID) {
          return false;
        }

        return inp.do_open(tmp, file_type);
      },
      "open file with inferred file type"
    )
    .def("molecules_remaining", &data_source_and_type<Molecule>::molecules_remaining,
      "Number of molecules yet to be read"
    )
    .def("next",
      [](data_source_and_type<Molecule>& inp)->std::optional<Molecule>{
        Molecule* m = inp.next_molecule();
        if (m == nullptr) {
          return std::nullopt;
        }
        return std::move(*m);
      },
      "Next molecule"
    )
    .def("molecules_read",
      [](const data_source_and_type<Molecule>&inp) {
        return inp.molecules_read();
      },
      "molecules read"
    )
    .def("__iter__", 
      [](py::object s) {
        return ReaderIterator(s.cast<data_source_and_type<Molecule> &>(), s);
      },
      py::keep_alive<0, 1>()
    )

    .def("__repr__",
      [](const data_source_and_type<Molecule>& inp) {
        IWString result;
        result << "<Reader from " << inp.fname() << " read " << inp.molecules_read() << '>';
        return std::string(result.data(), result.size());
      }
    )
    // Never did get these to work, maybe revisit sometime...
    .def("__enter__",
      [](data_source_and_type<Molecule>& inp) {
        std::cerr << "data_source_and_type __enter__\n";
        return;
      }
    )
    .def("__exit__",
      [](data_source_and_type<Molecule>& inp, 
        pybind11::object exc_type, pybind11::object exc_value, pybind11::object traceback) {
          inp.do_close();
      }
    )
  ;

  py::class_<Molecule_Output_Object>(io, "Writer")
    .def(py::init<>())
    .def("set_verbose",
      [](Molecule_Output_Object& writer, int v) {
        return writer.set_verbose(v);
      },
      "verbosity"
    )
    .def("add_output_type", &Molecule_Output_Object::add_output_type, "Output type")
    .def("new_stem",
      [](Molecule_Output_Object& writer, const std::string&stem)->bool{
        const const_IWSubstring s(stem.data(), stem.size());
        return writer.new_stem(s);
      },
      "Open new file"
    )
    .def("set_molecules_per_file", &Molecule_Output_Object::set_molecules_per_file, "Molecules per file")
    .def("molecules_written", &Molecule_Output_Object::molecules_written, "Molecules written")
    .def("write",
      [](Molecule_Output_Object& writer, Molecule& m)->bool{
        return writer.write(m);
      },
      "Write molecule"
    )
    .def("flush", &Molecule_Output_Object::do_flush, "flush")
    .def("close", &Molecule_Output_Object::do_close, "close")
  ;

  py::enum_<FileType>(io, "FileType")
    .value("SMI", FILE_TYPE_SMI)
    .value("USMI", FILE_TYPE_USMI)
    .value("SDF", FILE_TYPE_SDF)
    .export_values();
  ;

  py::class_<ReaderContext>(io, "MolReaderContext")
    .def(py::init([](const std::string& fname, FileType file_type,
                     const ReaderOptions& options) {
        ReaderContext* result = new ReaderContext(fname, file_type);
        if (! result->ApplyOptions(options)) {
          delete result;
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
        return result;
      }),
      py::arg("fname"), py::arg("file_type"), py::kw_only(),
      py::arg("options"))
    .def(py::init([](const std::string& fname, const ReaderOptions& options) {
        ReaderContext* result = new ReaderContext(fname);
        if (! result->ApplyOptions(options)) {
          delete result;
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
        return result;
      }),
      py::arg("fname"), py::kw_only(),
      py::arg("options"))
    .def(py::init([](const std::string& fname, FileType file_type,
                     bool largest_fragment, bool remove_chirality,
                     bool remove_cis_trans_bonds, bool remove_isotopes,
                     bool keep_sdf_tags, const std::string& sdf_identifier,
                     bool sdf_tags_to_json,
                     bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid) {
        ReaderContext* result = new ReaderContext(fname, file_type);
        result->SetPreprocessing(largest_fragment, remove_chirality,
                                 remove_cis_trans_bonds, remove_isotopes);
        result->_reader->mutable_molecule_read_options()
            .mdl_file_supporting_material()
            .set_read_extra_text_info(keep_sdf_tags);
        if (! result->SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                                    first_sdf_tag, prepend_sdfid)) {
          delete result;
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
        return result;
      }),
      py::arg("fname"), py::arg("file_type"), py::kw_only(),
      py::arg("largest_fragment") = false,
      py::arg("remove_chirality") = false,
      py::arg("remove_cis_trans_bonds") = false,
      py::arg("remove_isotopes") = false,
      py::arg("keep_sdf_tags") = false,
      py::arg("sdf_identifier") = "",
      py::arg("sdf_tags_to_json") = false,
      py::arg("all_sdf_tags") = false,
      py::arg("first_sdf_tag") = false,
      py::arg("prepend_sdfid") = true)
    .def(py::init([](const std::string& fname, bool largest_fragment,
                     bool remove_chirality, bool remove_cis_trans_bonds,
                     bool remove_isotopes, bool keep_sdf_tags,
                     const std::string& sdf_identifier, bool sdf_tags_to_json,
                     bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid) {
        ReaderContext* result = new ReaderContext(fname);
        result->SetPreprocessing(largest_fragment, remove_chirality,
                                 remove_cis_trans_bonds, remove_isotopes);
        result->_reader->mutable_molecule_read_options()
            .mdl_file_supporting_material()
            .set_read_extra_text_info(keep_sdf_tags);
        if (! result->SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                                    first_sdf_tag, prepend_sdfid)) {
          delete result;
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
        return result;
      }),
      py::arg("fname"), py::kw_only(),
      py::arg("largest_fragment") = false,
      py::arg("remove_chirality") = false,
      py::arg("remove_cis_trans_bonds") = false,
      py::arg("remove_isotopes") = false,
      py::arg("keep_sdf_tags") = false,
      py::arg("sdf_identifier") = "",
      py::arg("sdf_tags_to_json") = false,
      py::arg("all_sdf_tags") = false,
      py::arg("first_sdf_tag") = false,
      py::arg("prepend_sdfid") = true)
    .def("__enter__",
      [](ReaderContext& me)->ReaderContext& {
        return me;
      },
      py::return_value_policy::reference_internal
    )
    .def("__iter__", [](ReaderContext& me)->ReaderContext& { return me; },
      py::return_value_policy::reference_internal)
    .def("__next__",
      [](ReaderContext& me)->Molecule {
        std::optional<Molecule> m = me.Next();
        if (! m) {
          throw py::stop_iteration();
        }
        return std::move(*m);
      })
    .def("next", &ReaderContext::Next)
    .def("apply_options",
      [](ReaderContext& me, const ReaderOptions& options) {
        if (! me.ApplyOptions(options)) {
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
      })
    .def("set_preprocessing",
      [](ReaderContext& me, bool largest_fragment, bool remove_chirality,
         bool remove_cis_trans_bonds, bool remove_isotopes) {
        me.SetPreprocessing(largest_fragment, remove_chirality,
                            remove_cis_trans_bonds, remove_isotopes);
      },
      py::kw_only(),
      py::arg("largest_fragment") = false,
      py::arg("remove_chirality") = false,
      py::arg("remove_cis_trans_bonds") = false,
      py::arg("remove_isotopes") = false)
    .def("preprocessing_active",
      [](const ReaderContext& me)->bool {
        return me._preprocessing.active();
      })
    .def("set_sdf_options",
      [](ReaderContext& me, const std::string& sdf_identifier, bool sdf_tags_to_json,
         bool all_sdf_tags, bool first_sdf_tag, bool prepend_sdfid) {
        if (! me.SetSdfOptions(sdf_identifier, sdf_tags_to_json, all_sdf_tags,
                               first_sdf_tag, prepend_sdfid)) {
          throw py::value_error("Invalid sdf_identifier regular expression");
        }
      },
      py::kw_only(),
      py::arg("sdf_identifier") = "",
      py::arg("sdf_tags_to_json") = false,
      py::arg("all_sdf_tags") = false,
      py::arg("first_sdf_tag") = false,
      py::arg("prepend_sdfid") = true)
    .def("reset_sdf_options", &ReaderContext::ResetSdfOptions)
    .def("set_ignore_connection_table_errors",
      [](ReaderContext& me, int value) {
        me._reader->set_connection_table_errors_allowed(value);
      },
      "specify number of connection table errors to ignore"
    )
    .def("connection_table_errors_encountered",
      [](ReaderContext& me) {
        return me._reader->connection_table_errors_encountered();
      }
    )
    .def("molecules_remaining",
      [](ReaderContext& me)->uint64_t {
        if (! me._reader) {
          return 0;
        }
        return me._reader->records_remaining();
      },
      "number of molecules yet to be read - fails on pipes"
    )
    .def("molecules_read",
      [](const ReaderContext& me)->uint64_t {
        if (! me._reader) {
          return 0;
        }
        return me._reader->molecules_read();
      },
      "molecules read"
    )
    .def("__exit__",
      [](ReaderContext& me,
        pybind11::object exc_type, pybind11::object exc_value, pybind11::object traceback) {
          if (me._reader) {
            me._reader->do_close();
          }
          me.ResetSdfOptions();
      }
    )
  ;

  io.attr("ReaderContext") = io.attr("MolReaderContext");

  py::class_<ContextWriter>(io, "MolWriterContext")
    .def(py::init<const std::string&, FileType>())
    .def(py::init<const std::string&>())
    .def("__enter__",
      [](ContextWriter& me) {
        return me._writer.get();
      },
      py::return_value_policy::reference
    )
    .def("__exit__",
      [](ContextWriter& me, 
        pybind11::object exc_type, pybind11::object exc_value, pybind11::object traceback) {
          me._writer->do_close();
      }
    )
    .def("write",
      [](ContextWriter& me, Molecule& m)->bool {
        return me._writer->write(m);
      }
    )
    .def("set_molecules_per_file",
      [](ContextWriter& me, int value) {
        return me._writer->set_molecules_per_file(value);
      },
      "Constructs a sequence of files with N in each"
    )
    .def("set_name_token_for_file_name",
      [](ContextWriter& me, int token) {
        return me._writer->set_name_token_for_file_name(token);
      },
      "File names created from a token in the molecule name"
    )
  ;

  io.attr("ContextWriter") = io.attr("MolWriterContext");

  io.def("slurp",
    [](const std::string& fname)->std::optional<std::vector<Molecule>> {
      IWString tmp(fname);
      FileType itype = discern_file_type_from_name(tmp);
      if (itype == FILE_TYPE_INVALID) {
        itype = FILE_TYPE_SMI;  // give it a try.
      }
      data_source_and_type<Molecule> input(itype, tmp);
      if (! input.good()) {
        std::cerr << "slurp:cannot open '" << fname << "'\n";
        return std::nullopt;
      }

      uint32_t number_molecules = input.molecules_remaining();
      std::vector<Molecule> result(number_molecules);
      for (uint32_t i = 0; i < number_molecules; ++i) {
        if (! input.next_molecule(result[i])) {
          return std::nullopt;
        }
      }

      return result;
    },
    "Read all molecules from `fname`"
  );

  io.def("set_sdf_identifier",
    [](const std::string& sdfid)->bool {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      const_IWSubstring tmp(sdfid);
      return mdlfos->set_sdf_identifier(tmp);
    },
    "Set sdf tag containing identifier"
  );
  io.def("set_ignore_bad_m",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_ignore_unrecognised_mdl_m_records(v);
    },
    "ignore unrecognised M records while reading .sdf input"
  );
  io.def("set_mdlquiet",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_report_unrecognised_records(!v);
    },
    "quietly ignore unrecognised M records while reading .sdf input"
  );
  io.def("set_prepend_sdfid",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_prepend_sdfid(v);
    },
    "Prepend the sdf tag name to the identifier ->   tag:value"
  );
  io.def("set_allsdfid",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_fetch_all_sdf_identifiers(v);
    },
    "All sdf tags in the input concatenated to make the molecule name"
  );
  io.def("set_sdf_tags_to_json",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_sdf_tags_to_json(v);
    },
    "Molecule name is JSON representation of sdf tags"
  );
  io.def("set_firstsdftag",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_take_first_tag_as_name(v);
    },
    "The first sdf tag in the input is the molecule name"
  );
  io.def("set_allow_deuterium",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_allow_deuterium(v);
    },
    "D is interpreted as Dueterium [2H]"
  );
  io.def("set_allow_tritium",
    [](bool v) {
      MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
      mdlfos->set_allow_tritium(v);
    },
    "T is interpreted as Tritium [3H]"
  );


}
