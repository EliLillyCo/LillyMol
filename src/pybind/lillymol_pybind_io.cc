#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"

namespace py = pybind11;

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

  public:
    ReaderContext(const std::string& fname, FileType file_type);
    ReaderContext(const std::string& fname);
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

  py::class_<ReaderContext>(io, "ReaderContext")
    .def(py::init<const std::string&, FileType>())
    .def(py::init<const std::string&>())
    .def("__enter__",
      [](ReaderContext& me) {
        return me._reader.get();
      },
      py::return_value_policy::reference
    )
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
      [](ReaderContext& me) {
        if (! me._reader) {
          return 0;
        }
        return me._reader->records_remaining();
      },
      "number of molecules yet to be read - fails on pipes"
    )
    .def("__exit__",
      [](ReaderContext& me, 
        pybind11::object exc_type, pybind11::object exc_value, pybind11::object traceback) {
          me._reader->do_close();
      }
    )
  ;

  py::class_<ContextWriter>(io, "ContextWriter")
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

}
