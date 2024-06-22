// Given a molecule, write different representations as tabular form.
// This is for dealing with the structure equivalence problem.
// The intent is tha the output from this program could be used to
// populate a table in a relational database, enabling joining on
// various structure abstractions.

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace representations {

using std::cerr;

void
Usage(int rc) {
  cerr << R"(Generates various unique smiles forms in tabular format
 -O ...          enable optional features, enter '-O help' for possibilities
 -g ...          apply chemical standardisation before any smiles formed
 -l              reduce to largest fragment before any smiles formed
 -c              remove chirality before any smiles formed
 -o <char>       output separator, space by default, try '-o tab' or '-o ,'
 -v              verbose output
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // These columns are optional.
    int _scaffold;
    int _scaffold_with_isotopes;

    IWString _output_separator;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions

    int AppendScaffold(Molecule& m, IWString_and_File_Descriptor& output);


  public:
    Options();

    int Initialise(Command_Line& cl);

    int WriteHeader(IWString_and_File_Descriptor& output) const;
    int Preprocess(Molecule& m);
    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

Options::Options() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _scaffold = 0;
  _scaffold_with_isotopes = 0;
  _output_separator = ' ';
  _input_type = FILE_TYPE_INVALID;
}

void
DisplayDashOOptions(std::ostream& output) {
  output << " -O scaffold       generate a column containing the molecular scaffold\n";
  output << " -O Iscaffold      generate a column containing the molecular scaffold with isotopes\n";

  exit(0);
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('o')) {
    cl.value('o', _output_separator);
    if (! char_name_to_char(_output_separator)) {
      cerr << "Invalid -o option, the following directives are understood\n";
      char_name_to_char_usage('o');
      return 1;
    }
  }

  if (cl.option_present('O')) {
    const_IWSubstring o;
    for (int i = 0; cl.value('O', o, i); ++i) {
      if (o == "scaffold") {
        _scaffold = 1;
        if (_verbose) {
          cerr << "Will generate a column for the molecular scaffold\n";
        }
      } else if (o == "Iscaffold") {
        _scaffold_with_isotopes = 1;
        if (_verbose) {
          cerr << "Will generate a column for the isotopicaly labelled scaffold\n";
        }
      } else if (o == "help") {
        DisplayDashOOptions(cerr);
      } else {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        DisplayDashOOptions(cerr);
      }
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  m.revert_all_directional_bonds_to_non_directional();

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "smiles";
  output << _output_separator << "id";
  output << _output_separator << "usmi";
  output << _output_separator << "mformula";

  if (_chemical_standardisation.active()) {
    output << _output_separator << "std";
  }

  output << _output_separator << "largest_frag";
  output << _output_separator << "largest_frag_mformula";

  output << _output_separator << "nochiral";
  output << _output_separator << "tautomer";

  if (_scaffold) {
    output << _output_separator << "scaffold";
  }
  if (_scaffold_with_isotopes) {
    output << _output_separator << "Iscaffold";
  }


  output << '\n';

  return 1;
}

// Place isotope 1 in `iso` for every atom that is in the scaffold
// and joins an atom in `spinach`.
void
StoreIsotopes(Molecule& m, const int *spinach, isotope_t* iso) {
  for (const Bond * b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (spinach[a1] == spinach[a2]) {
      continue;
    }
    if (spinach[a1]) {
      iso[a2] = 1;
    } else if (spinach[a2]) {
      iso[a1] = 1;
    }
  }
}

int
Options::Process(Molecule& m,
            IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  output << m.smiles() << _output_separator << m.name();
  output << _output_separator << m.unique_smiles();

  IWString mformula;
  m.formula_distinguishing_aromatic(mformula);
  output << _output_separator << mformula;

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
    output << _output_separator << m.unique_smiles();
  }

  Molecule m2(m);
  m2.reduce_to_largest_fragment_carefully();
  output << _output_separator << m2.unique_smiles();

  m2.formula_distinguishing_aromatic(mformula);
  output << _output_separator << mformula;

  m2.remove_all_chiral_centres();
  output << _output_separator << m2.unique_smiles();

  std::unique_ptr<Molecule> mcopy;
  if (_scaffold || _scaffold_with_isotopes) {
    mcopy.reset(new Molecule(m));
  }

  Mol2Graph mol2graph;
  mol2graph.TurnOnMostUsefulOptions();
  m2.change_to_graph_form(mol2graph);
  output << _output_separator << m2.unique_smiles() << ':' << mformula;

  if (_scaffold || _scaffold_with_isotopes) {
    AppendScaffold(*mcopy, output);
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::AppendScaffold(Molecule& m, IWString_and_File_Descriptor& output) {

  const int matoms = m.natoms();

  std::unique_ptr<int[]> spinach(new_int(matoms));
  m.identify_spinach(spinach.get());

  std::unique_ptr<isotope_t[]> isotope;
  if (_scaffold_with_isotopes) {
    isotope.reset(new isotope_t[matoms]);
    std::fill_n(isotope.get(), matoms, 0);
    StoreIsotopes(m, spinach.get(), isotope.get());
  }

  static constexpr char kNoScaffold = '.';

  m.remove_atoms(spinach.get(), 1);
  if (m.empty()) {
    if (_scaffold) {
      output << _output_separator << kNoScaffold;
    }
    if (_scaffold_with_isotopes) {
      output << _output_separator << kNoScaffold;
    }

    return 1;
  }

  if (_scaffold) {
    output << _output_separator << m.unique_smiles();
  }

  if (! _scaffold_with_isotopes) {
    return 1;
  }

  // Adjust the isotope array for the atoms that have been removed.
  int ndx = 0;
  for (int i = 0; i < matoms; ++i) {
    if (spinach[i]) {
      // Atom will be removed.
    } else {
      isotope[ndx] = isotope[i];
      ++ndx;
    }
  }

  // NOte that the isotope array is likely longer than `m` but that is OK.
  m.set_isotopes(isotope.get());
  output << _output_separator << m.unique_smiles();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Options:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  return 1;
}

int
MakeRepresentations(Options& options,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
MakeRepresentations(Options& options,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      return 0;
    }

    if (! MakeRepresentations(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
MakeRepresentations(Options& options,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = options.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return MakeRepresentations(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lco:O:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  options.WriteHeader(output);

  for (const char* fname : cl) {
    if (! MakeRepresentations(options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace representations

int
main(int argc, char** argv) {
  int rc = representations::Main(argc, argv);

  return rc;
}
