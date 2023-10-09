// Computes the element density of a given atom type.

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace elemental_density {

using std::cerr;

void
Usage(int rc) {
  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    Substructure_Query _query;

    int _molecules_read = 0;

    extending_resizable_array<int> _acc_natoms; 
    Accumulator<double> _acc_density;
    extending_resizable_array<int> _acc_nhits; 
    extending_resizable_array<int> _acc_density_values;

    IWString_and_File_Descriptor _stream_for_individual_values;

    // If working as a filter, we can impose a range on what gets written.
    float _min_density = 0.0;
    float _max_density = 1.0;
    int _rejected_by_filters;

  // private functions.
    int WriteIndividualValue(Molecule& m, int nhits, float fraction);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }

    int Process(Molecule& m);
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _min_density = 0.0f;
  _max_density = 1.0f;
  _rejected_by_filters = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

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

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (! cl.option_present('s')) {
    cerr << "Must specify smarts\n";
    return 0;
  }

  if (cl.option_present('s')) {
    IWString smarts = cl.string_value('s');
    if (! _query.create_from_smarts(smarts)) {
      cerr << "Invalid smarts\n";
      return 0;
    }
  }

  if (cl.option_present('d')) {
    if (! cl.value('d', _min_density) || _min_density < 0.0f || _min_density > 1.0f) {
      cerr << "Invalid lower element density (-d)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will skip molecules with an elemental density < " << _min_density << '\n';
    }
  }

  if (cl.option_present('D')) {
    if (! cl.value('D', _max_density) || _max_density < _min_density || _max_density > 1.0f) {
      cerr << "Invalid upper element density (-D)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will skip molecules with an elemental density > " << _max_density << '\n';
    }
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }
    if (! _stream_for_individual_values.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for all values '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Individual values written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }
  for (int i = 0; i < _acc_natoms.number_elements(); ++i) {
    if (_acc_natoms[i] > 0) {
      output << _acc_natoms[i] << " molecules had " << i << " atoms\n";
    }
  }
  for (int i = 0; i  < _acc_nhits.number_elements(); ++i) {
    if (_acc_nhits[i]) {
      output << _acc_nhits[i] << " molecules had " << i << " hits\n";
    }
  }
  output << "Density btw " << _acc_density.minval() << " and " << _acc_density.maxval() << " ave " << _acc_density.average() << '\n';
  int n = 0;
  for (int i = 0; i < _acc_density_values.number_elements(); ++i) {
    n += _acc_density_values[i];
  }

  int sum = 0.0;
  for (int i = 0; i < _acc_density_values.number_elements(); ++i) {
    if (_acc_density_values[i] == 0) {
      continue;
    }

    const float f = static_cast<float>(i) / 100.0f;
    sum += _acc_density_values[i];
    output << i << " fraction " << f << " cumulative " << sum << " fraction " << iwmisc::Fraction<float>(sum, n) << '\n';
  }

  if (_min_density > 0.0f || _max_density < 1.0f) {
    output << _rejected_by_filters << " molecules rejected for out of range (" << _min_density << ',' << _max_density << ")\n";
  }

  return output.good();
}

int
Options::MaybeDiscernInputType(const char * fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }
  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m) {
  ++_molecules_read;
  const int matoms = m.natoms();
  ++_acc_natoms[matoms];
  Molecule_to_Match  target(&m);
  Substructure_Results sresults;
  const int nhits = _query.substructure_search(target, sresults);
  ++_acc_nhits[nhits];
  double fraction = iwmisc::Fraction<double>(nhits, matoms);
  _acc_density.extra(fraction);
  ++_acc_density_values[static_cast<int>(fraction * 100.0f)];

  if (_stream_for_individual_values.active()) {
    WriteIndividualValue(m, nhits, fraction);
  }
  return 1;
}

int
Options::WriteIndividualValue(Molecule& m,
                              int nhits,
                              float fraction) {
  if (fraction < _min_density || fraction > _max_density) {
    ++_rejected_by_filters;
    return 1;
  }
  _stream_for_individual_values << m.smiles();
  _stream_for_individual_values << ' ' << m.name();
  _stream_for_individual_values << ' ' << m.natoms();
  _stream_for_individual_values << ' ' << nhits;
  _stream_for_individual_values << ' ' << fraction;
  _stream_for_individual_values << '\n';
  _stream_for_individual_values.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
ElementalDensity(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! options.Process(*m)) {
      return 0;
    }
  }

  return 1;
}

int
ElementalDensity(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "ElementalDensity:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ElementalDensity(options, input, output);
}

int
ElementalDensity(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:Alcg:s:S:d:D:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! ElementalDensity(options, fname, output)) {
      cerr << "ElementalDensity::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.Report(std::cout);

  return 0;
}

}  // namespace elemental_density

int
main(int argc, char ** argv) {

  int rc = elemental_density::ElementalDensity(argc, argv);

  return rc;
}
