// Generate a fingerprint based on formal charges.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace charge_fingerprint {

using std::cerr;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(
Generates formal charge fingerprints.
Charges come either from existing formal charges in the molecule or from the charge assigner.
 -N ...                 charge assigner specifications
 -J <tag>               fingerprint tag
 -p <n>                 number of bit replicates
 -w                     also fingerprint the number of saturated O and N atoms that are neutral.
 -g ...                 chemical standardisation
 -v                     verbose output
)";
// clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    Charge_Assigner _charge_assigner;

    IWString _tag;

    int _bit_replicates;

    // We can optionally count the number of saturated O and N atoms
    // that are neutral.
    int _fingerprint_neutral_and_saturated;

    int _zwitterions;

    int _molecules_read = 0;

    // Statistics across all the molecules we process;
    int _neutral_molecules;
    extending_resizable_array<int> _molecules_with_negative_charges;
    extending_resizable_array<int> _molecules_with_positive_charges;
    extending_resizable_array<int> _molecules_with_neutral_nitrogen;
    extending_resizable_array<int> _molecules_with_neutral_oxygen;

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Fingerprint(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 1;
  _remove_chirality = 1;
  _bit_replicates = 1;
  _tag = "NCFC<";
  _fingerprint_neutral_and_saturated = 0;
  _molecules_read = 0;
  _zwitterions = 0;
  _neutral_molecules = 0;
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

  if (cl.option_present('N')) {
    if (! _charge_assigner.construct_from_command_line(cl, _verbose, 'N')) {
      cerr << "Options::Initialise:cannot initialise charge assigner (-N)\n";
      return 0;
    }
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _bit_replicates) || _bit_replicates < 1) {
      cerr << "The number of bit replicates (-p) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will generate " << _bit_replicates << " bit replicates\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    _tag.EnsureEndsWith('<');
    if (_verbose) {
      cerr << "Fingerprints generated with tag " << _tag << '\n';
    }
  }

  if (cl.option_present('w')) {
    _fingerprint_neutral_and_saturated = 1;
    if (_verbose) {
      cerr << "Will also fingerprint the number of saturated, neutral N and O atoms\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _neutral_molecules << " neutral molecules " <<
                iwmisc::Fraction<float>(_neutral_molecules, _molecules_read) << '\n';
  output << _zwitterions << " zwitterions\n";
  for (int i = 0; i < _molecules_with_positive_charges.number_elements(); ++i) {
    if (_molecules_with_positive_charges[i]) {
      output << _molecules_with_positive_charges[i] << " had " << i << " positive charges\n";
    }
  }
  for (int i = 0; i < _molecules_with_negative_charges.number_elements(); ++i) {
    if (_molecules_with_negative_charges[i]) {
      output << _molecules_with_negative_charges[i] << " had " << i << " negative charges\n";
    }
  }

  if (_fingerprint_neutral_and_saturated) {
    for (int i = 0; i < _molecules_with_neutral_nitrogen.number_elements(); ++i) {
      if (_molecules_with_neutral_nitrogen[i]) {
        output << _molecules_with_neutral_nitrogen[i] << " had " << i << " uncharged saturated nitrogen\n";
      }
    }
    for (int i = 0; i < _molecules_with_neutral_oxygen.number_elements(); ++i) {
      if (_molecules_with_neutral_oxygen[i]) {
        output << _molecules_with_neutral_oxygen[i] << " had " << i << " uncharged saturated oxygen\n";
      }
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
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

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Fingerprint(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  if (! Preprocess(m)) {
    return 0;
  }

  if (_charge_assigner.active()) {
    _charge_assigner.process(m);
  }

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  int npos = 0;
  int nneg = 0;
  // Count neutral, saturated nitrogen and oxygen atoms.
  int neutral_oxygen = 0;
  int neutral_nitrogen = 0;

  // Examine all the oxygen and nitrogen atoms. Ignore other charges.
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];

    atomic_number_t z = a.atomic_number();
    if (z == 7) {
    } else if (z == 8) {
    } else {
      continue;
    }

    const formal_charge_t q = a.formal_charge();
    if (q > 0) {
      ++npos;
    } else if (q < 0) {
      ++nneg;
    } else if (m.is_aromatic(i)) {
      continue;
    }  else if (a.unsaturation()) {
    } else if (7 == z) {
      ++neutral_nitrogen;
    } else if (8 == z) {
      ++neutral_oxygen;
    }
  }

  if (_verbose ) {
    if (nneg == 0 && npos == 0) {
      ++_neutral_molecules;
    }
    if (nneg > 0 && npos > 0) {
      ++_zwitterions;
    }

    ++_molecules_with_negative_charges[nneg];
    ++_molecules_with_positive_charges[npos];
    ++_molecules_with_neutral_nitrogen[neutral_nitrogen];
    ++_molecules_with_neutral_oxygen[neutral_oxygen];
  }

  // conditions for no output
  int empty = 0;
  if (_fingerprint_neutral_and_saturated) {
    if (npos == 0 && nneg == 0 &&
        neutral_nitrogen == 0 && neutral_oxygen == 0) {
      empty = 1;
    }
  } else {
    if (npos == 0 && nneg == 0) {
      empty = 1;
    }
  }

  if (empty) {
    output << _tag << ">\n";
    return 1;
  }

  Sparse_Fingerprint_Creator sfc;

  // If we are not fingerprinting the neutral, saturated heteroatoms,
  // we are done if there are no charges found.
  if (! _fingerprint_neutral_and_saturated) {
    for (int i = 0; i < _bit_replicates; ++i) {
      if (npos) {
        sfc.hit_bit(i * 4, npos);
      }
      if (nneg) {
        sfc.hit_bit(i * 4 + 1, npos);
      }
    }
  } else {
    for (int i = 0; i < _bit_replicates; ++i) {
      if (npos) {
        sfc.hit_bit(i * 4, npos);
      }
      if (nneg) {
        sfc.hit_bit(i * 4 + 1, npos);
      }
      if (neutral_oxygen) {
        sfc.hit_bit(i * 4 + 2, neutral_oxygen);
      }
      if (neutral_nitrogen) {
        sfc.hit_bit(i * 4 + 3, neutral_nitrogen);
      }
    }
  }

  IWString fp;
  sfc.daylight_ascii_form_with_counts_encoded(_tag, fp);
  output << fp << '\n';

  return 1;
}

int
ChargeFingerprint(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  output << smiles_tag << m.smiles() << ">\n";
  output << identifier_tag << m.name() << ">\n";
  options.Fingerprint(m, output);
  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
ChargeFingerprintFilterLine(Options& options,
                            const_IWSubstring& line,
                            IWString_and_File_Descriptor& output) {
  assert(line.starts_with(smiles_tag));
  assert(line.ends_with(">"));

  line.remove_leading_chars(smiles_tag.length());
  line.chop();

  Molecule m;
  if (! m.build_from_smiles(line)) {
    cerr << "ChargeFingerprintFilterLine:cannot parse '" << line << "'\n";
    return 0;
  }

  if (! options.Fingerprint(m, output)) {
    cerr << "ChargeFingerprintFilterLine:error fingerprinting " << line << '\n';
    return 0;
  }

  return 1;
}

int
ChargeFingerprintFilter(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    output << line << '\n';
    if (! line.starts_with(smiles_tag)) {
      continue;
    }

    if (! ChargeFingerprintFilterLine(options, line, output)) {
      cerr << "ChargeFingerprintFilter:error processing " << line << '\n';
      return 0;
    }
    output.write_if_buffer_holds_more_than(4092);
  }

  return 1;
}

int
ChargeFingerprintFilter(Options& options,
                const char* fname,
                IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ChargeFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return ChargeFingerprintFilter(options, input, output);
}

int
ChargeFingerprint(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);


    if (! ChargeFingerprint(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ChargeFingerprint(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ChargeFingerprint:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ChargeFingerprint(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:N:T:A:lcg:i:J:p:fw");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  int work_as_tdt_filter = 0;

  if (cl.option_present('f')) {
    work_as_tdt_filter = 1;
  } else if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (work_as_tdt_filter) {
    if (! ChargeFingerprintFilter(options, cl[0], output)) {
      cerr << "ChargeFingerprint:error processing '" << cl[0] << "'\n";
      return 1;
    }
  } else {
    for (const char * fname : cl) {
      if (! ChargeFingerprint(options, fname, input_type, output)) {
        cerr << "ChargeFingerprint::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace charge_fingerprint

int
main(int argc, char ** argv) {

  int rc = charge_fingerprint::Main(argc, argv);

  return rc;
}
