/*
  Computes ring substitution fingerprint
*/

#include <stdlib.h>

#include <cstdint>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/ring_substitution.h"

namespace ring_substitution_main {

using std::cerr;

using ring_substitution::RingSubstitutionGenerator;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(
  -J <tag>       tag for fingerprints
  -P <bonds>     maximum path length around edge of ring system
  -M ...         options for specifying just what kind of information to generate
  -l             reduce to largest fragment
  -g ...         chemical standardisation options
  -f             work as a TDT filter
  -i <type>      input specification
  -A ...         aromaticity specifications, enter '-A help' for info
  -E ...         element specifications, enter '-E help' for info
  -v             verbose output
)";

  exit(rc);
}

void
DisplayDashMOptions(std::ostream& os) {
  os << R"(
  -M posi    set bits to indicate positions on rings only
  -M simp    set bits to indicate positions and simple atom types
  -M full    set bits to indicate positions and complete atom types
  -M sfb     set bits to indicate presence or absence of single features
  -M 01      remove count information from fingerprints, just presence or absence
)";

  exit(2);
}

class Options {
 private:
  uint64_t _molecules_read = 0;
  uint64_t _molecules_with_no_rings = 0;

  IWString _fingerprint_tag = "NCRS<";

  Chemical_Standardisation _chemical_standardisation;

  bool _reduce_to_largest_fragment = false;
  int _verbose = 0;
  bool _flatten_fingerprints_to_01 = false;

  IWString _smiles_tag = "$SMI<";
  IWString _identifier_tag = "PCN<";

  bool _work_as_filter = false;

  Accumulator_Int<unsigned int> _acc_nset;

  RingSubstitutionGenerator _generator;

  void Preprocess(Molecule& m);
  bool Process(Molecule& m, IWString_and_File_Descriptor& output);
  bool Process(data_source_and_type<Molecule>& input,
               IWString_and_File_Descriptor& output);
  bool Process(const char* fname, FileType input_type,
               IWString_and_File_Descriptor& output);
  bool ProcessFilterRecord(const const_IWSubstring& buffer,
                           IWString_and_File_Descriptor& output);
  bool ProcessFilter(iwstring_data_source& input, IWString_and_File_Descriptor& output);
  bool ProcessFilter(const char* fname, IWString_and_File_Descriptor& output);
  void Report(std::ostream& output) const;

 public:
  int Main(int argc, char** argv);
};

void
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }
}

bool
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  Preprocess(m);

#ifdef DEBUG_RING_SUBSTITUTION
  for (int i = 0; i < m.natoms(); i++) {
    cerr << "Atom " << i << " assigned type " << atype[i] << '\n';
  }
#endif

  Sparse_Fingerprint_Creator sfpc;

  if (m.nrings() > 0) {
    _generator.Generate(m, sfpc);
  } else {
    ++_molecules_with_no_rings;
  }

  if (!_work_as_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  if (_flatten_fingerprints_to_01) {
    sfpc.flatten_to_01();
  }

  IWString tmp;
  sfpc.daylight_ascii_form_with_counts_encoded(_fingerprint_tag, tmp);

  output << tmp << '\n';

  if (_verbose > 2) {
    sfpc.debug_print(cerr);
  }

  if (!_work_as_filter) {
    output << "|\n";
  }

  if (_verbose) {
    _acc_nset.extra(sfpc.nset());
  }

  if (_verbose > 2) {
    cerr << _molecules_read << " '" << m.name() << "' hits " << sfpc.nbits() << " bits\n";
  }

  return true;
}

bool
Options::Process(data_source_and_type<Molecule>& input,
                 IWString_and_File_Descriptor& output) {
  Molecule* m;

  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!Process(*m, output)) {
      return false;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return true;
}

bool
Options::Process(const char* fname, FileType input_type,
                 IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(input_type != FILE_TYPE_INVALID);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return false;
  }

  if (_verbose > 2) {
    input.set_verbose(1);
  }

  return Process(input, output);
}

bool
Options::ProcessFilterRecord(const const_IWSubstring& buffer,
                             IWString_and_File_Descriptor& output) {
  const_IWSubstring tmp(buffer);
  tmp.remove_leading_chars(_smiles_tag.length());
  tmp.chop(1);

  Molecule m;
  if (!m.build_from_smiles(tmp)) {
    cerr << "Invalid smiles '" << tmp << "'\n";
    return false;
  }

  return Process(m, output);
}

bool
Options::ProcessFilter(iwstring_data_source& input,
                       IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(_smiles_tag)) {
      continue;
    }

    if (!ProcessFilterRecord(buffer, output)) {
      cerr << "Fatal error, line " << input.lines_read() << '\n';
      return false;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return true;
}

bool
Options::ProcessFilter(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return false;
  }

  return ProcessFilter(input, output);
}

void
Options::Report(std::ostream& output) const {
  output << "Processed " << _molecules_read << " molecules, " << _molecules_with_no_rings
         << " had no rings\n";
  if (_acc_nset.n() > 1) {
    output << "Fingerprints have between " << _acc_nset.minval() << " and "
           << _acc_nset.maxval() << " bits set, ave "
           << static_cast<float>(_acc_nset.average()) << '\n';
  }
}

int
Options::Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:J:lg:M:fP:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  _verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, _verbose > 1)) {
    cerr << "Cannot process -A option\n";
    Usage(11);
  }

  if (!process_elements(cl, _verbose > 1, 'E')) {
    cerr << "Cannot initialise elements\n";
    Usage(8);
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      Usage(32);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = true;

    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _fingerprint_tag);
    _fingerprint_tag.EnsureEndsWith('<');

    if (_verbose) {
      cerr << "Fingerprints written with tag '" << _fingerprint_tag << "'\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  } else if (cl.option_present('f')) {
    ;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('P')) {
    int max_path_length;
    if (!cl.value('P', max_path_length) || max_path_length < 2) {
      cerr << "The maximum path length (-P) must be a whole +ve number\n";
      Usage(7);
    }

    if (_verbose) {
      cerr << "Max path length " << max_path_length << '\n';
    }
    _generator.set_max_path_length(max_path_length);
  }

  if (cl.option_present('M')) {
    int nset = 0;

    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++)) {
      if (m.starts_with("posi")) {
        _generator.set_positional_information_only(true);
        nset++;
        if (_verbose) {
          cerr << "Only positional information will be recorded\n";
        }
      } else if (m.starts_with("simp")) {
        _generator.set_positional_information_only(false);
        _generator.set_simple_atom_types(true);
        nset++;
        if (_verbose) {
          cerr << "Simple atom types for substitutent atoms\n";
        }
      } else if ("full" == m) {
        _generator.set_positional_information_only(false);
        _generator.set_full_atom_types(true);
        nset++;
        if (_verbose) {
          cerr << "Full atom typing\n";
        }
      } else if ("sfb" == m) {
        _generator.set_place_single_feature_bits(true);
        if (_verbose) {
          cerr << "Will set bits according to presence or absence of features\n";
        }
      } else if ("01" == m) {
        _flatten_fingerprints_to_01 = true;
        if (_verbose) {
          cerr << "Fingerprints will be flattened to just 01 forms\n";
        }
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        DisplayDashMOptions(cerr);
        Usage(4);
      }
    }

    if (nset != 1) {
      cerr << "Sorry, must specify one, and exactly one of 'posi', 'simple' or 'full' "
              "for substituent atom types\n";
      Usage(6);
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f')) {
    _work_as_filter = true;
    if (cl.number_elements() != 1) {
      cerr << "Filter mode requires exactly one input file\n";
      return 3;
    }
    rc = ProcessFilter(cl[0], output) ? 0 : 1;
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!Process(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (_verbose) {
    Report(cerr);
  }

  return rc;
}

}  // namespace ring_substitution_main

int
main(int argc, char** argv) {
  ring_substitution_main::Options options;
  return options.Main(argc, argv);
}
