// Reads a file that contains previously sorted list of smiles.
// Typically this would be done by appending the molecular formula
// to the smiles, sorting that file, and then using this tool.

#include <iostream>

#include "absl/container/flat_hash_set.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace unique_molecules_sorted {

using std::cerr;

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
  cerr << R"(Detects duplicates in a file of sorted tokens. The assumption is that
if the tokens are different (for example molecular formula) then the molecules cannot
be the same. 
 -c <column>        column that contains the sort key.
 -v                 verbose output
)";
// clang-format on

  ::exit(rc);
}

// For each group of molecules that share a token, we need data to determine
// if a new molecule is unique or not.
class PreviousMolecule {
  private:
    // Unique smiles we have encountered so far.
    absl::flat_hash_set<IWString> _seen;

    int _remove_chirality;

    int _reduce_to_largest_fragment;

    // THe number of instances of the current molecule.
    uint32_t _number_instances;

    // THe token used for matching groups of molecules.
    IWString _token;

    // As a new group is started, we store the smiles of the first
    // molecule.
    // If we ever get a second molecule to compare, then we will need
    // to convert this smiles to a Molecule object and get the unique
    // smiles to check _seen.
    // If there is only one item in the group, we will never perform
    // smiles interpretation.
    IWString _first_molecule_smiles;

    IWString _id;  // is this needed?

    Molecule _m;

    // Private functions.

    bool BuildFromSmiles(const const_IWSubstring& smiles);

  public:
    PreviousMolecule();

    void set_remove_chirality(int s) {
      _remove_chirality = s;
    }
    void set_reduce_to_largest_fragment(int s) {
      _reduce_to_largest_fragment = s;
    }

    uint32_t number_instances() const {
      return _number_instances;
    }

    void FirstOfNewGroup(const const_IWSubstring& smiles,
                         const const_IWSubstring& token);

    int Reset(const const_IWSubstring& buffer,
              const const_IWSubstring& smiles,
              const const_IWSubstring& id,
              const const_IWSubstring& token);
    
    // Does a new token match what we have?
    bool SameToken(const const_IWSubstring& token) {
      return _token == token;
    }

    // Return true if the molecule represented by `smiles` has
    // not been encountered before.
    bool IsNew(const const_IWSubstring& smiles);
};

PreviousMolecule::PreviousMolecule() {
  _number_instances = 0;
  _remove_chirality = 0;
  _reduce_to_largest_fragment = 0;
}

// Do not do smiles interpretation at first, because there may be
// only one item in the group.
void
PreviousMolecule::FirstOfNewGroup(const const_IWSubstring& smiles,
                const const_IWSubstring& token) {
  _first_molecule_smiles = smiles;
  _token = token;
  _m.resize(0);
  _seen.clear();
  _number_instances = 1;
}

// build `_m` from `smiles`
bool
PreviousMolecule::BuildFromSmiles(const const_IWSubstring& smiles) {
  if (! _m.build_from_smiles(smiles)) {
    cerr << "PreviousMolecule::BuildFromSmiles:invalid smiles\n";
    cerr << smiles << '\n';
    return false;  // not sure what to return, maybe an optional<bool>..
  }

  if (_remove_chirality) {
    _m.remove_all_chiral_centres();
  }

  if (_reduce_to_largest_fragment) {
    _m.reduce_to_largest_fragment_carefully();
  }

  return true;
}

bool
PreviousMolecule::IsNew(const const_IWSubstring& smiles) {
  // If we have not done any smiles interepretation, we may need to now.
  if (_seen.empty()) {
    // If the smiles is the same as the previous smiles, this is a dup.
    // No need to do smiles interpretation yet.
    if (_first_molecule_smiles == smiles) {
      ++_number_instances;
      return 0;
    }

    // SMiles interpretation needed.
    if (! BuildFromSmiles(_first_molecule_smiles)) {
      return false;  // not sure what to return, maybe an optional<bool>..
    }

    // The unique smiles of the first molecule.
    _seen.insert(_m.unique_smiles());
    ++_number_instances;
  }

  // Now check this new molecule.

  if (! BuildFromSmiles(smiles)) {
    return false;  // not sure what to return, maybe an optional<bool>..
  }

  const IWString& usmi = _m.unique_smiles();

  if (auto iter = _seen.find(usmi); iter == _seen.end()) {
    _seen.insert(usmi);
    ++_number_instances;
    return true;
  }
  
  return false;
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    // The column that contains the sorted key.
    int _column = 0;

    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    uint64_t _molecules_read = 0;
    uint64_t _molecules_written = 0;

    PreviousMolecule _previous_molecule;

    extending_resizable_array<uint32_t> _duplicates;

    // Private functions.

    int GetTokenFromCol(const const_IWSubstring& buffer,
                const_IWSubstring& smiles,
                const_IWSubstring& id,
                const_IWSubstring& result);
    void DoWrite(const const_IWSubstring& buffer,
        IWString_and_File_Descriptor& output);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    // int Preprocess(Molecule& m);

    int UniqueMolecules(iwstring_data_source& input,
                IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _column = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

//if (cl.option_present('T')) {
//  if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
//    Usage(8);
//}

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
    _previous_molecule.set_reduce_to_largest_fragment(1);
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality\n";
    }
    _previous_molecule.set_remove_chirality(1);
  }

  if (cl.option_present('C')) {
    if (! cl.value('C', _column) || _column <= 0) {
      cerr << "Invalid column specification (-c)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Sorted keys in column " << _column << '\n';
    }
    --_column;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Wrote " << _molecules_written << " "
            << (_molecules_read - _molecules_written) << " duplicates\n";
  for (int i = 0; i < _duplicates.number_elements(); ++i) {
    if (_duplicates[i] == 0) {
      continue;
    }
    output << _duplicates[i] << " molecules had " << i << " variants\n";
  }

  return 1;
}

#ifdef NOT_USED_HERE
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

//if (_element_transformations.active()) {
//  _element_transformations.process(m);
//}

  return 1;
}
#endif

// `buffer` must contains smiles id and some number of other tokens.
// We place the first two tokens into `smiles` and `id`.
// Then look for the token in _column and place that in `result`.
int
Options::GetTokenFromCol(const const_IWSubstring& buffer,
                const_IWSubstring& smiles,
                const_IWSubstring& id,
                const_IWSubstring& result) {
  const_IWSubstring token;
  for (int col = 0, i = 0; buffer.nextword(token, i); ++col) {
    if (col == 0) {
      smiles = token;
    } else if (col == 1) {
      id = token;
    } else if (col == _column) {
      result = token;
      return 1;
    }
  }

  // Did not find _column - and maybe also not the others.
  cerr << "Cannot get column " << (_column + 1) << " from\n";
  cerr << buffer << '\n';
  return 0;
}

void
Options::DoWrite(const const_IWSubstring& buffer,
        IWString_and_File_Descriptor& output) {
  output << buffer << '\n';
  output.write_if_buffer_holds_more_than(4096);

  ++_molecules_written;
}

int
Options::UniqueMolecules(iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  // Within each group of identical tokens, all the unique
  // smiles that have been seen before. Will be reset every
  // time the token changes.

  // Scope here for efficiency, and to remove stuff from the loop.
  const_IWSubstring smiles; 
  const_IWSubstring id; 
  const_IWSubstring token;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! GetTokenFromCol(buffer, smiles, id, token)) {
      return 0;
    }

    ++_molecules_read;

    //cerr << smiles << ' ' << id << ' ' << token << " begin processing\n";

    // If the token is different from prev, we definitely have a new
    // molecule.
    if (! _previous_molecule.SameToken(token)) {
      static bool first_call = true;
      if (first_call) {
        first_call = false;
      } else {
        ++_duplicates[_previous_molecule.number_instances()];
      }

      DoWrite(buffer, output);

      _previous_molecule.FirstOfNewGroup(smiles, token);

      continue;
    }

    // Token is the same as the previous molecule we need to
    // generate a unique smiles to check.

    if (_previous_molecule.IsNew(smiles)) {
      DoWrite(buffer, output);
    }
  }

  ++_duplicates[_previous_molecule.number_instances()];

  return 1;
}

#ifdef NO_LONGER_USED
int
Options::CheckByUniqueSmiles(const const_IWSubstring& buffer,
                PreviousMolecule& previous_molecule,
                IWString_and_File_Descriptor& output) {
  Molecule newmol;

  const_IWSubstring smiles(buffer);
  smiles.truncate_at_first(' ');

  if (! newmol.build_from_smiles(smiles)) {
    cerr << "Options::CheckByUniqueSmiles:invalid smiles\n";
    cerr buffer << '\n';
    return 0;
  }

  if (! Preprocess(newmol)) {
      return 0;
    }
  }

  if (newmol.unique_smiles() == previous_molecule.unique_smiles) {
    ++previous_molecule.count;
    return 1;
  }

  return 0;
}
#endif

int
UniqueMolecules(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "UniqueMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  return options.UniqueMolecules(input, output);
}

int
UniqueMolecules(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE::A:lcC:");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! UniqueMolecules(options, fname, output)) {
      cerr << "UniqueMolecules::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace unique_molecules_sorted

int
main(int argc, char ** argv) {

  int rc = unique_molecules_sorted::UniqueMolecules(argc, argv);

  return rc;
}
