// Generate a molecular formula fingerprint

#include <cstdint>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/mformula.h"

namespace formula_fingperprint {

using std::cerr;
using mformula::MFormula;

void
Usage() {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Fingerprint based on the molecular formula.
 -J <tag>       tag to use.
 -f             function as a TDT filter.
 -v             verbose output.
)";
  // clang-format off

  ::exit(0);
}

class Options {
  private:
    int _verbose;

    int _function_as_tdt_filter;

    IWString _smiles_tag;
    IWString _identifier_tag;
    IWString _tag;

  public:
    Options();

    int Initialise(Command_Line& cl);

    const IWString& smiles_tag() const {
      return _smiles_tag;
    }

    int WriteFingerprint(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _function_as_tdt_filter = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _tag = "NCFML<";
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (_verbose) {
      cerr << "Will function as a tdt filter\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    _tag.EnsureEndsWith('<');
  }

  return 1;
}

int
Options::WriteFingerprint(Molecule& m,
                          IWString_and_File_Descriptor& output) {
  
  MFormula mf;

  mf.Build(m);

  IWString fp;
  mf.ToSparseFingerprint(fp);

  if (_function_as_tdt_filter) {
    output << _tag << fp << ">\n";
    return 1;
  }

  output << _smiles_tag << m.smiles() << ">\n";
  output << _identifier_tag << m.name() << ">\n";
  output << _tag << fp << ">\n";
  output << "|\n";

  return 1;
}

int
Options::Report(std::ostream& output) const {
  return 1;
}

int
FormulaFingerprintFilterRecord(const_IWSubstring line,  // note local copy
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  line.remove_leading_chars(options.smiles_tag().length());
  line.chop();

  Molecule m;
  if (! m.build_from_smiles(line)) {
    cerr << "FormulaFingerprintFilterRecord:invalid smiles\n";
    return 0;
  }

  return options.WriteFingerprint(m, output);
}

int
FormulaFingerprintFilter(iwstring_data_source& input,
                         Options& options,
                         IWString_and_File_Descriptor& output) {

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
    if (! buffer.starts_with(options.smiles_tag())) {
      continue;
    }

    if (! FormulaFingerprintFilterRecord(buffer, options, output)) {
      cerr << "FormulaFingerprintFilter:error processing\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
FormulaFingerprintFilter(const char* fname,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "FormulaFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return FormulaFingerprintFilter(input, options, output);
}

int
FormulaFingerprint(Molecule& m,
                   Options& options,
                   IWString_and_File_Descriptor& output) {
  return options.WriteFingerprint(m, output);
}

int
FormulaFingerprint(data_source_and_type<Molecule>& input,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    FormulaFingerprint(*m, options, output);
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
FormulaFingerprint(const char* fname,
                   FileType input_type,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);
                          
  if (! input.good()) {
    cerr << "FormulaFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return FormulaFingerprint(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:EJ:fi:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage();
  }

  const int verbose = cl.option_present('v');

  Options options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage();
  }

  const int function_as_tdt_filter = cl.option_present('f');

  if (function_as_tdt_filter && cl.option_present('i')) {
    cerr << "Cannot combine TDT filter (-f) and input file (-i) options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage();
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficent arguments\n";
    Usage();
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    int rc;
    if (function_as_tdt_filter) {
      rc = FormulaFingerprintFilter(fname, options, output);
    } else {
      rc = FormulaFingerprint(fname, input_type, options, output);
    }

    if (rc == 0) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace formula_fingperprint

int 
main(int argc, char ** argv) {
  int rc = formula_fingperprint::Main(argc, argv);

  return rc;
}
