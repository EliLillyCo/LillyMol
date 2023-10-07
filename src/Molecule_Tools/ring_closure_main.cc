// Generate molecular variants created by adding ring closure bonds.
// Options are read from a RingClosure::ring_closure textproto.

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/ring_closure.h"

namespace ring_closure_main {

using std::cerr;

using ring_closure::RingClosureOptions;

void
Usage(int rc) {
  cerr << "Generates ring closure variants\n";
  cerr << " -C <proto>       RingClosure::ring_closure text proto with ring closure rules\n";
  cerr << " -p .             write the parent molecule\n";
  cerr << " -p <str>         write the parent molecule, appending <str> to name\n";
  cerr << " -u <str>         append <str> to the parent name in generated molecules\n";
  cerr << " -c               remove all chirality\n";
  cerr << " -l               reduce to largest fragment\n";
  cerr << " -S <stem>        file name stem for output\n";
  cerr << " -o <type>        type of output file\n";
  cerr << " -g ...           chemical standardisation\n";
  cerr << " -A ...           standard aromaticity options\n";
  cerr << " -E ...           standard element options\n";
  cerr << " -i <type>        input type\n";
  cerr << " -v               verbose output\n";
  ::exit(rc);
}

class LocalOptions {
  private:
    int _verbose;

    // If set, write the parent molecule to the output.
    int _write_parent;

    // If set, will append `_parent_suffix` to each parent written.
    IWString _parent_suffix;

    // If set, we will append _variant_suffix << ndx to each
    // newly generated molecule.
    IWString _variant_suffix;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions

    int MaybeWriteParent(Molecule& m, Molecule_Output_Object& output);
    int MaybeSetName(Molecule& m, int ndx) const;

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m,
                      RingClosureOptions& options,
                      Molecule_Output_Object& output);
};

LocalOptions::LocalOptions() {
  _verbose = 0;
  _write_parent = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
}

int
LocalOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      return 0;
    }
  }

  if (cl.option_present('p')) {
    IWString p = cl.string_value('p');
    _write_parent = 1;
    if (p != '.') {
      _parent_suffix << ' ' << p;
    }
  }

  if (cl.option_present('u')) {
    cl.value('u', _variant_suffix);
    if (_verbose) {
      cerr << "Wil append '" << _variant_suffix << "' to each variant\n";
    }
  }

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return 1;
}

int
LocalOptions::MaybeWriteParent(Molecule& m,
                Molecule_Output_Object& output) {
  if (! _write_parent) {
    return 0;
  }

  if (_parent_suffix.empty()) {
    return output.write(m);
  }

  const IWString save_name(m.name());

  IWString new_name(m.name());
  new_name << _parent_suffix;
  m.set_name(new_name);
  output.write(m);
  m.set_name(save_name);
  return 1;
}

// If this looks like a 3d molecule, setup 3d smiles.
void
Check3D(const Molecule& m) {
  if (m.highest_coordinate_dimensionality() == 3) {
    set_append_coordinates_after_each_atom(1);
  } else {
    set_append_coordinates_after_each_atom(0);
  }
}

int
LocalOptions::Process(Molecule& m,
                      RingClosureOptions& options,
                      Molecule_Output_Object& output) {
  Check3D(m);

  MaybeWriteParent(m, output);

  resizable_array_p<Molecule> variants = options.Process(m);

  for (int i = 0; i < variants.number_elements(); ++i) {
    Molecule* v = variants[i];
    MaybeSetName(*v, i);
    output.write(*v);
  }

  return 1;
}

// `m` is the `ndx` variant of a starting molecule.
// Set its name if needed.
int
LocalOptions::MaybeSetName(Molecule& m, int ndx) const  {
  if (_variant_suffix.empty()) {
    return 0;
  }

  static constexpr char kSep = '_';

  IWString append_to_name;
  append_to_name << kSep << _variant_suffix << kSep << ndx;
  m << append_to_name;

  return 1;
}

int
RingClosure(Molecule& m,
            RingClosureOptions& options,
            LocalOptions& local_options,
            Molecule_Output_Object& output) {
  return local_options.Process(m, options, output);
}
        
int
RingClosure(data_source_and_type<Molecule>& input,
            RingClosureOptions& options,
            LocalOptions& local_options,
            Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    local_options.Preprocess(*m);
    if (! RingClosure(*m, options, local_options, output)) {
      return 0;
    }
  }

  return 1;
}

int
RingClosure(const char* fname,
            FileType input_type,
            RingClosureOptions& options,
            LocalOptions& local_options,
            Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "Cannot discern input type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return RingClosure(input, options, local_options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:S:o:C:p:u:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(1);
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  set_copy_name_in_molecule_copy_constructor(1);
  set_display_no_kekule_form_message(0);
  set_allow_two_electron_systems_to_be_aromatic(1);
  set_warn_aromatic_chain_atoms(0);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  RingClosureOptions options;
  if (! options.Initialise(cl, 'C')) {
    cerr << "Cannot initialise ring closure options (-C)\n";
    return 1;
  }

  LocalOptions local_options;
  if (! local_options.Initialise(cl)) {
    cerr << "Cannot initialise job options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  Molecule_Output_Object output;
  if (! cl.option_present('S')) {
    cerr << "Must specify output file stem via the -S option\n";
    Usage(1);
  }

  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type (-o)\n";
    return 1;
  }

  if (cl.option_present('S')) {
    const IWString stem = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, stem)) {
      cerr << "Cannot overwrite input file(s) '" << stem << "'\n";
      return 1;
    }
    if (! output.new_stem(stem)) {
      cerr << "Cannot initialise output stem '" << stem << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written to '" << stem << "'\n";
    }
  }

  for (const char* fname : cl) {
    if (! RingClosure(fname, input_type, options, local_options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace ring_closure_main

int
main(int argc, char ** argv) {

  int rc = ring_closure_main::Main(argc, argv);

  return rc;
}
