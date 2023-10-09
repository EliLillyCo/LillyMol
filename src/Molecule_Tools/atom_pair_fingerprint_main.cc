// Compute atom pair fingerprints.

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_pair_fingerprint.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

namespace atom_pair_fingerprint {

const char* prog_name = nullptr;

int verbose = 0;

int molecules_read = 0;

Chemical_Standardisation chemical_standardisation;

Atom_Typing_Specification atom_typing;

int reduce_to_largest_fragment = 0;

IWString tag;

bool function_as_tdt_filter = false;

bool check_for_collisions = false;

const IWString smiles_tag("$SMI<");
const IWString identifier_tag("PCN<");

// Potentially interesting statistics on molecules processed.

bool gather_molecule_statistics = false;
Accumulator_Int<uint> atom_count;
Accumulator_Int<uint> longest_path_acc;
extending_resizable_array<int> longest_path;

// If positive, the width of the fixed width fingerprint to produce.
int write_fixed_width_fingerprint = 0;

void
usage(int rc) {
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes atom pair fingerprints\n";
  cerr << "  -r <sep>      minimum bond separation (def 1)\n";
  cerr << "                set to 0 to get single atom type fingerprint included\n";
  cerr << "  -R <sep>      maximum bond separation (def none)\n";
  cerr << "  -P ...        atom type specification\n";
  cerr << "  -J <tag>      tag for fingerprints\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -X <fname>    look for bits in <fname> and provide explanations\n";
  cerr << "  -B <fname>    write all bits found to <fname>\n";
  cerr << "  -w <nbits>    generate fixed width binary fingerprints, `nbits` bits\n";
  cerr << "  -t            truncate all pairs beyond max_separation\n";
  cerr << "  -y            check for bit collisions\n";
  cerr << "  -c            produce labelled molecules with coverage\n";
  cerr << "  -s            gather statistics on molecules processed\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

using atom_pair_fingerprint::AtomPairFingerprint;

void
Preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

void
GatherStatistics(Molecule& m) {
  atom_count.extra(m.natoms());

  const int lp = m.longest_path();

  longest_path_acc.extra(lp);
  longest_path[lp]++;
}

int
DoAtomPairFingerprint(Molecule& m, const atom_type_t* atype,
                      AtomPairFingerprint& atom_pair_fp_gen,
                      IWString_and_File_Descriptor& output) {
  Sparse_Fingerprint_Creator sfc;

  // A single atom molecule would produce zero bits, or two atoms in different
  // fragments...
  if (!atom_pair_fp_gen.Fingerprint(m, nullptr, atype, sfc)) {
    cerr << "AtomPairFingerprint:fingerprinting failed " << m.smiles() << " ignored\n";
  }

  if (gather_molecule_statistics) {
    GatherStatistics(m);
  }

  if (!function_as_tdt_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  if (write_fixed_width_fingerprint > 0) {
    const IWString ascii = sfc.FixedWidthFingerprint(write_fixed_width_fingerprint);
    output << tag << ascii << ">\n";
  } else {
    IWString tmp;
    sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << "\n";
  }

  if (!function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

int
DoAtomPairFingerprint(Molecule& m, AtomPairFingerprint& atom_pair_fp_gen,
                      IWString_and_File_Descriptor& output) {
  atom_type_t* atype = new atom_type_t[m.natoms()];
  std::unique_ptr<atom_type_t[]> free_atype(atype);

  if (!atom_typing.assign_atom_types(m, atype)) {
    cerr << "AtomPairFingerprint::cannot assign atom types " << m.smiles() << ' '
         << m.name() << '\n';
    return 0;
  }

  return DoAtomPairFingerprint(m, atype, atom_pair_fp_gen, output);
}

int
DoAtomPairFingerprint(data_source_and_type<Molecule>& input,
                      AtomPairFingerprint& atom_pair_fp_gen,
                      IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    Preprocess(*m);

    if (!DoAtomPairFingerprint(*m, atom_pair_fp_gen, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
DoAtomPairFingerprintPipe(const_IWSubstring line,  // Note local copy
                          AtomPairFingerprint& atom_pair_fp_gen,
                          IWString_and_File_Descriptor& output) {
  assert(line.ends_with('>'));
  assert(line.starts_with(smiles_tag));

  line.remove_leading_chars(smiles_tag.length());
  line.chop();

  Molecule m;

  if (!m.build_from_smiles(line)) {
    cerr << "AtomPairFingerprintPipe:invalid smiles '" << line << "'\n";
    return 0;
  }

  Preprocess(m);

  return DoAtomPairFingerprint(m, atom_pair_fp_gen, output);
}

int
DoAtomPairFingerprintPipe(iwstring_data_source& input,
                          AtomPairFingerprint& atom_pair_fp_gen,
                          IWString_and_File_Descriptor& output) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    output << line << '\n';

    if (!line.starts_with(smiles_tag)) {
      continue;
    }

    if (!DoAtomPairFingerprintPipe(line, atom_pair_fp_gen, output)) {
      cerr << "AtomPairFingerprintPipe:invalid input " << line << "' line "
           << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
DoAtomPairFingerprint(const char* fname, FileType input_type,
                      AtomPairFingerprint& atom_pair_fp_gen,
                      IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (function_as_tdt_filter) {
    iwstring_data_source input(fname);
    if (!input.good()) {
      cerr << "AtomPairFingerprint::cannot open filter " << fname << endl;
      return 0;
    }

    return DoAtomPairFingerprintPipe(input, atom_pair_fp_gen, output);
  }

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return DoAtomPairFingerprint(input, atom_pair_fp_gen, output);
}

int
DoAtomPairFingerprint(int argc, char** argv) {
  Command_Line cl(argc, argv, "A:K:lg:i:J:P:bvftr:R:ysB:cw:");

  if (cl.unrecognised_options_encountered()) {
    usage(1);
  }

  verbose = cl.option_count('v');

  (void)process_elements(cl);

  if (!process_standard_smiles_options(cl, verbose)) {
    usage(7);
  }

  set_global_aromaticity_type(Daylight);

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(8);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  int min_separation = 0;
  if (cl.option_present('r')) {
    if (cl.value('r', min_separation) || min_separation < 1) {
      cerr << "The min atom pair separation (-r) must be > 0\n";
      return 1;
    }

    if (verbose) {
      cerr << "Will only fingerprint pairs at least " << min_separation
           << " bonds apart\n";
    }
  }

  int max_separation = std::numeric_limits<int>::max();
  if (cl.option_present('R')) {
    if (!cl.value('R', max_separation) || max_separation < min_separation) {
      cerr << "Max separation (-R) must be larger than min_separation " << min_separation
           << endl;
      return 1;
    }

    if (verbose) {
      cerr << "Will only fingerprint pairs separated by " << max_separation
           << " bonds or less\n";
    }
  }

  bool coverage = false;

  if (cl.option_present('c')) {
    coverage = true;
    if (verbose) {
      cerr << "Will generate coverage information\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (!atom_typing.build(p)) {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 1;
    }
  } else {
    atom_typing.build("UST:Y");
    if (verbose) {
      cerr << "Default compressed atomic number atom typing\n";
    }
  }

  if (cl.option_present('f')) {
    function_as_tdt_filter = 1;

    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }
  }

  if (cl.option_present('w')) {
    if (!cl.value('w', write_fixed_width_fingerprint) ||
        write_fixed_width_fingerprint < 8) {
      cerr << "The number of bits in a fixed fingerprint must be +ve\n";
      return 1;
    }

    if (verbose) {
      cerr << "Will generate fixed width fingerprints " << write_fixed_width_fingerprint
           << " bits\n";
    }
  }

  if (!cl.option_present('J')) {
    cerr << "Must specify tag via the -J option\n";
    usage(1);
  }

  cl.value('J', tag);
  if (!tag.ends_with('<')) {
    tag += '<';
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_tdt_filter) {
    if (cl.number_elements() > 1) {
      cerr << "When working as a filter, only one argument possible\n";
      usage(1);
    }
  } else if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot discern input type(s)\n";
      return 3;
    }
  } else if (!process_input_type(cl, input_type)) {
    cerr << prog_name << ": cannot discern input type\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    gather_molecule_statistics = true;
    if (verbose) {
      cerr << "Will gather statistics on molecules processed\n";
    }
  }

  AtomPairFingerprint atom_pair_fp_gen;
  atom_pair_fp_gen.set_min_separation(min_separation);
  atom_pair_fp_gen.set_max_separation(max_separation);
  atom_pair_fp_gen.set_examine_coverage(coverage);
  if (cl.option_present('y')) {
    atom_pair_fp_gen.set_check_for_collisions(true);
    check_for_collisions = true;
    if (verbose) {
      cerr << "WIll check for atom pair collisions\n";
    }
  }

  if (cl.option_present('t')) {
    atom_pair_fp_gen.set_include_out_of_range_separations(true);
    if (verbose) {
      cerr << "Will include out of range pairs in truncated form\n";
    }
  }

  if (cl.option_present('b')) {
    atom_pair_fp_gen.set_fingerprint_bonded_atoms_with_btype(true);
    if (verbose) {
      cerr << "Bonded atom pairs fingerprinted according to bond type\n";
    }
  }

  if (cl.option_present('B')) {
    const char* fname = cl.option_value('B');

    if (!atom_pair_fp_gen.OpenStreamForBitMeanings(fname)) {
      cerr << "Cannot open stream for bit meanings '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Bit meanings written to '" << fname << "'\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << prog_name << ": no inputs\n";
    usage(1);
  }

  IWString_and_File_Descriptor output(1);

  int rc;
  for (int i = 0; i < cl.number_elements(); ++i) {
    if (!DoAtomPairFingerprint(cl[i], input_type, atom_pair_fp_gen, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose || check_for_collisions) {
    cerr << molecules_read << " molecules read\n";
    if (check_for_collisions) {
      atom_pair_fp_gen.ReportCollisions(cerr);
    }
  }

  if (gather_molecule_statistics) {
    cerr << "Molecules had btw " << atom_count.minval() << " and " << atom_count.maxval()
         << " atoms, ave " << atom_count.average() << "\n";
    cerr << "Longest paths btw " << longest_path_acc.minval() << " and "
         << longest_path_acc.maxval() << " bonds, ave " << longest_path_acc.average()
         << "\n";
    for (int i = 0; i < longest_path.number_elements(); ++i) {
      if (longest_path[i]) {
        cerr << longest_path[i] << " molecules had a longest path of " << i << " bonds\n";
      }
    }
  }

  return rc;
}

}  // namespace atom_pair_fingerprint

int
main(int argc, char** argv) {
  atom_pair_fingerprint::prog_name = argv[0];

  int rc = atom_pair_fingerprint::DoAtomPairFingerprint(argc, argv);

  return rc;
}
