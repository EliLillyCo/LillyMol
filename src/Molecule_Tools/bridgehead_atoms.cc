// Places isotopes on bridghead atoms. Ultimately this might become
// part of substruture searches, but let's explore the idea here.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace bridgehead_atoms /* replace all occurrences */ {

using std::cerr;

using molecule_processing::MoleculePreprocessing;

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
  cerr << R"(Performs some task on a set of molecules.
 -a          what the -a option does
# commonly used LillyMol options.
 -I <isotope> isotope applied to matched atoms.
 -b           only put isotopic labels in bridghead atoms.
 -i <type>    input type, -i sdf -i ICTE -i mdlquiet -i SDFID:IDNUMBER
 -E ...       element options: -E autocreate -E anylength
 -A ...       aromaticity options: -A 2
 -g ...       chemical standardisation: -g all
 -T ...       element transformations: -T I=Cl -T Br=Cl
 -l           reduce to largest fragment
 -c           remove chirality
 -v           verbose output
)";
// clang-format on

  ::exit(rc);
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    MoleculePreprocessing _preprocessing;

    isotope_t _isotope = 0;

    resizable_array_p<Substructure_Query> _queries;

    // We can label all atoms involved in the strongly fused system or just what are likely to
    // be the bridgehead atoms.

    bool _label_bridgehead_only = true;

    extending_resizable_array<uint32_t> _matches;

    int _write_molecules_with_no_bridghead_atoms = 0;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint64_t _molecules_read = 0;

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

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _molecules_read = 0;
  _isotope = 0;
  _write_molecules_with_no_bridghead_atoms = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (! _preprocessing.Initialise(cl)) {
    cerr << "Options::Initialise:cannot initialise preprocessing\n";
    return 1;
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T')) {
        Usage(8);
    }
  }

  if (cl.option_present('q')) {
  }

  if (cl.option_present('I')) {
    if (! cl.value('I', _isotope)) {
      cerr << "Invalid isotope specification\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Bridgehead atoms labelled with " << _isotope << "\n";
    }
  }

  if (cl.option_present('b')) {
    _label_bridgehead_only = true;
    if (_verbose) {
      cerr << "Will only label bridghead atoms\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Processed " << _molecules_read << " molecules\n";
  for (int i = 0; i < _matches.number_elements(); ++i) {
    if (_matches[i] > 0) {
      cerr << _matches[i] << " molecules had " << i << " bridgehead atoms\n";
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  if (m.nrings() < 1) {
    if (_write_molecules_with_no_bridghead_atoms) {
      output << m.smiles() << ' ' << m.name() << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  std::unique_ptr<int[]> to_process(new_int(matoms));

  int got_matches = 0;

  for (int i = 0; i < m.nrings(); ++i) {
    const Ring* ri = m.ringi(i);
    if (ri->is_aromatic()) {
      continue;
    }

    if (ri->largest_number_of_bonds_shared_with_another_ring() < 2) {
      continue;
    }

    cerr << ri->ring_number() << " has " << ri->fused_ring_neighbours() << " fused neighbours\n";

    for (const Ring* rj : ri->fused_neighbours()) {
      if (ri->largest_number_of_bonds_shared_with_another_ring() < 2) {
        continue;
      }
      if (rj->is_aromatic()) {
        continue;
      }

      if (ri->compute_bonds_shared_with(*rj) == 1) {  // simple fused.
        continue;
      }

      cerr << "Strongly fused rings\n";
      ri->increment_vector(to_process.get(), 1);
      rj->increment_vector(to_process.get(), 1);
      ++got_matches;
    }
  }

  if (got_matches == 0) {
    if (_write_molecules_with_no_bridghead_atoms) {
      output << m.smiles() << ' ' << m.name() << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
    return 1;
  }

  int count = 0;
  for (int i = 0; i < matoms; ++i) {
    if (to_process[i] == 0) [[likely]] {
      continue;
    }

    if (_label_bridgehead_only) {
      if (m.ncon(i) > 2 && m.ring_bond_count(i) > 2)  {
        m.set_isotope(i, _isotope);
        ++count;
      }
    } else {
      m.set_isotope(i, _isotope);
      ++count;
    }
  }

  ++_matches[count];

  output << m.smiles() << ' ' << m.name() << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

// Replace all occurrences of BridgeHeadAtoms with a name appropriate
// for your application.
int
BridgeHeadAtoms(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
BridgeHeadAtoms(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! BridgeHeadAtoms(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
BridgeHeadAtoms(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "BridgeHeadAtoms:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return BridgeHeadAtoms(options, input, output);
}

int
BridgeHeadAtoms(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:s:q:I:b");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! BridgeHeadAtoms(options, fname, input_type, output)) {
      cerr << "BridgeHeadAtoms::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace bridgehead_atoms

int
main(int argc, char ** argv) {

  int rc = bridgehead_atoms::BridgeHeadAtoms(argc, argv);

  return rc;
}
