// A common task is to ask "is this/these molecule(s) in another set"

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecular_formula.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"


namespace grep_molecule {

using std::cerr;
using molecular_formula::MolecularFormula;

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
  cerr << R"(Look for one or molecules as exact matches in another set
Designed to work like grep, where the first argument is a patten and subsequent arguments
are files to be searched.
Multiple patterns (smiles) can be specified on the command line separated by commas
  'C methane,CC ethane,CCC propane'
will search for methane, ethane and propane in subsequent files. Quotes are essential.
 -f <fname>             like the -f option to grep, read the patterns (smiles) from <fname>
 -v                     verbose output
  )";

  ::exit(rc);
}

#ifdef NEED_TO_SEE_IF_THIS_WILL_HELP
  adding the -s option slowed things down, and this doubles down on that concept.
 -r <n>                 generate a hash containing <n> random smiles for each needle.
 -R <n>                 for each member of the haystack, check <n> random smiles
                        Suggest larger numbers for -r and smaller for -R
                        may work if the smiles are expected to be identical.
 -s                     store the starting smiles as text - do what real grep does.
                        Surprisingly, it usually slows things down.
#endif

struct MoleculeAndFormula {
  Molecule m;
  MolecularFormula<uint32_t> mf;
};

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    // We can also work exactly like grep.
    // As each needle is read, store the smiles and if there is a
    // match on a needle, we are done.
    int _store_starting_smiles = 0;
    IW_STL_Hash_Map_String _starting_smiles;

    // First check when a new molecule is read is whether or not any
    // of our needle molecules have a given atom count.
    extending_resizable_array<int> _atom_count;

    // the next step is a hash of molecular formula
    IW_STL_Hash_Map_uint _formula_to_find;

    // We can use the formula_distinguishing_aromatic() method to
    // generate a molecular formula that is more precise, but which
    // will still hopefully be cheaper than a unique smiles determination.
    IW_STL_Hash_Map_uint _aromatic_formula_to_find;

    // We load this hash with the unique smiles of the moleules to be found.
    // The value is the ID of the molecule needle with this unique smiles.
    // Note that we do not check for duplictes in the needles.
    IW_STL_Hash_Map_String _usmi_to_find;

    // Not clear if this is a good idea or not. If the smiles are written
    // by the same software, it might be just fine.

    // We can optionally read the needles and create some number of
    // random smiles for those molecules. Then we can do a string comparison
    // on the needle smiles and maybe we get lucky.
    uint32_t _number_random_needle_variants = 0;
    IW_STL_Hash_Map_String _random_needle_variants;
    // It might make sense to generate 100 random variants for each needle
    // molecule, but it will not make sense to do that for the haystack.
    // Likely generate a smaller number.
    uint32_t _number_random_haystack_variants = 0;

    int _molecules_read = 0;
    int _failed_formula_generation = 0;
    int _no_atom_count_match = 0;
    int _no_formula_match = 0;
    int _no_aromatic_formula_match = 0;
    int _unique_smiles_comparisons = 0;
    int _matched_as_smiles_string = 0;
    int _matches_found = 0;

  // Private functions

    int ReadNeedles(IWString& fname);
    int ReadNeedles(data_source_and_type<Molecule>& input);
    int ReadSmiles(const IWString& smiles);
    int AddToHashes(Molecule& m);
    int AddToHashes(const const_IWSubstring buffer, Molecule& m);
    int RandomVariantMatches(const const_IWSubstring& buffer, Molecule& m,
                IWString_and_File_Descriptor& output);
    int MatchesSmilesAsText(const_IWSubstring buffer, IWString_and_File_Descriptor& output);
    int GotMatch(const const_IWSubstring& buffer, 
                 const IWString& name_of_needle, IWString_and_File_Descriptor& output);

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

    // Helpful when the -i option is not given.
    int MaybeDiscernInputType(const char * fname);

    int Process(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
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

  if (cl.option_present('f')) {
    IWString fname;
    for (int i = 0; cl.value('f', fname, i); ++i) {
      if (! ReadNeedles(fname)) {
        cerr << "Cannot read '" << fname << "'\n";
        return 0;
      }
    }
  } else if (cl.empty()) {
    cerr << "Insufficient arguments, must specify pattern and targets\n";
    Usage(1);
  } else if (cl.size() == 1) {
    cerr << "Insufficient arguments, must specify both pattern and targets\n";
    Usage(1);
  } else {
    IWString smiles(cl[0]);
    if (! ReadSmiles(smiles)) {
      cerr << "Cannot interpret as smiles '" << smiles << "'\n";
      return 0;
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _number_random_needle_variants)) {
      cerr << "The number of random smiles needle variants (-r) must be a +ve whole number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "For each needle, will hash " << _number_random_needle_variants <<
          " random smiles\n";
    }
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _number_random_haystack_variants)) {
      cerr << "The number of random smiles haystack variants (-R) must be a +ve whole number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "For each member of the haystack, will generate " << _number_random_haystack_variants <<
          " random smiles\n";
    }
  }

  if (cl.option_present('s')) {
    _store_starting_smiles = 1;
    if (_verbose) {
      cerr << "Will store the starting smiles and compare smiles as text\n";
    }
  }

  return 1;
}

int
Options::ReadNeedles(IWString& fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "Options::ReadNeedles:cannot open '" << fname << '\n';
    return 0;
  }

  return ReadNeedles(input);
}

int
Options::ReadNeedles(data_source_and_type<Molecule>& input) {
  Molecule* m;
  for (int ndx = 0; (m = input.next_molecule()) != nullptr; ++ndx) {
    std::unique_ptr<Molecule> free_m(m);

    if (! Preprocess(*m)) {
      return 0;
    }

    if (! AddToHashes(*m)) {
      cerr << "Options::ReadNeedles:invalid needle " << m->smiles() <<
              ' ' << m->name() << '\n';
      return 0;
    }
  }

  return _usmi_to_find.size();
}

int
Options::ReadSmiles(const IWString& smiles) {
  IWString token;
  for (int i = 0; smiles.nextword(token, i, ','); ) {
    Molecule m;
    if (!m.build_from_smiles(token)) {
      cerr << "Options::ReadSmiles:invalid smiles '" << token << "'\n";
      return 0;
    }
    if (!Preprocess(m)) {
      return 0;
    }
    if (! AddToHashes(smiles, m)) {
      cerr << "Options::ReadSmiles:cannot add to hash '" << token << "'\n";
      return 0;
    }
  }

  return _usmi_to_find.size();
}

int
Options::AddToHashes(const_IWSubstring buffer, Molecule& m) {
  buffer.truncate_at_first(' ');

  _starting_smiles[buffer] = m.name();
  
  return AddToHashes(m);
}

int
Options::AddToHashes(Molecule& m) {
  static MolecularFormula<uint32_t> mfgen;
  static IWString mf;

  ++_atom_count[m.natoms()];

  if (! mfgen.Build(m)) {
    cerr << "Options::AddToHashes:cannot build formula\n";
    return 0;
  }

  mfgen.MakeFormula(mf);

  // We could be more careful about this and keep track of the number
  // of each formula requested, and then drop them when found...
  _formula_to_find[mf] = 1;

  // No checking for duplicates
  IWString name(m.name());

  IWString s;
  m.formula_distinguishing_aromatic(s);
  _aromatic_formula_to_find[s] = 1;

  _starting_smiles[m.smiles()] = name;

  _usmi_to_find[m.unique_smiles()] = name;

  if (_number_random_needle_variants) {
    _random_needle_variants[m.smiles()] = m.smiles();
    for (uint32_t i = 0; i < _number_random_needle_variants; ++i) {
      IWString s = m.random_smiles();
      _random_needle_variants[s] = m.name();
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _failed_formula_generation << " failed formula generation\n";
  output << _no_atom_count_match << " no atom count match\n";
  output << _no_formula_match << " no formula match\n";
  output << _no_aromatic_formula_match << " no aromatic formula match\n";
  output << _unique_smiles_comparisons << " unique smiles comparisons\n";
  if (_store_starting_smiles) {
    output <<  _matched_as_smiles_string << " matched as smiles string match\n";
  }
  output << _matches_found << " matches found\n";

  return 1;
}

// If the input type is known, return it.
// Otherwise examine the file name's suffix to 
// determine the type.
int
Options::MaybeDiscernInputType(const char * fname) {

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
Options::Process(const const_IWSubstring& buffer,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  static MolecularFormula<uint32_t> mfgen;

  const int matoms = mfgen.Build(buffer);
  if (matoms == 0) {
    ++_failed_formula_generation;

    cerr << "Options::Process:cannot process '" << buffer << "' ignored\n";
    return 0;
  }

  if (_atom_count[matoms] == 0) {
    ++_no_atom_count_match;
    return 0;
  }

  IWString mf;
  mfgen.MakeFormula(mf);

  auto iter_formula = _formula_to_find.find(mf);
  if (iter_formula == _formula_to_find.end()) {
    ++_no_formula_match;
    return 0;
  }

  if (_store_starting_smiles) {
    if (MatchesSmilesAsText(buffer, output)) {
      return 1;
    }
  }

  // Now we need to do things the expensive way.
  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "Options::Process:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  Preprocess(m);

  if (_number_random_needle_variants) {
    if (RandomVariantMatches(buffer, m, output)) {
      return 1;
    }
  }

  // Aromatic formula.
  if (1) {
    IWString s;
    m.formula_distinguishing_aromatic(s);
    auto iter = _aromatic_formula_to_find.find(s);
    if (iter == _aromatic_formula_to_find.end()) {
      ++_no_aromatic_formula_match;
      return 0;
    }
  }

  const IWString& usmi = m.unique_smiles();
  ++_unique_smiles_comparisons;

  auto iter_usmi = _usmi_to_find.find(usmi);
  if (iter_usmi == _usmi_to_find.end()) {
    return 0;
  }

  // cerr << "Going to GotMatch " << buffer << '\n';
  return GotMatch(buffer, iter_usmi->second, output);
}

int
Options::GotMatch(const const_IWSubstring& buffer,
                  const IWString& name_of_needle,
                  IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

  output << buffer << kSep << name_of_needle << '\n';

  ++_matches_found;

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
Options::RandomVariantMatches(const const_IWSubstring& buffer,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  for (uint32_t i = 0; i < _number_random_haystack_variants; ++i) {
    const IWString& s = m.random_smiles();
    auto iter = _random_needle_variants.find(s);
    if (iter == _random_needle_variants.end()) {
      continue;
    }
    // cerr << "Random variant GotMatch\n";
    return GotMatch(buffer, iter->second, output);
  }

  return 0;
}

int
Options::MatchesSmilesAsText(const_IWSubstring buffer, IWString_and_File_Descriptor& output) {
  IWString mybuffer(buffer);
  mybuffer.truncate_at_first(' ');
  auto iter = _starting_smiles.find(mybuffer);
  if (iter == _starting_smiles.end()) {
    return 0;
  }

  ++_matched_as_smiles_string;

  // cerr << "MatchesSmilesAsText GotMatch\n";
  return GotMatch(buffer, iter->second, output);
}

int
ApplicationName(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    options.Process(buffer, output);
  }


  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {

  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  return ApplicationName(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:f:clr:R:s");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  int istart;
  if (cl.option_present('f')) {
    istart = 0;
  } else {
    istart = 1;
  }
  for (int i = istart; i < cl.number_elements(); ++i) {
    const char* fname = cl[i];
    if (! ApplicationName(options, fname, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace grep_molecule

int
main(int argc, char ** argv) {

  int rc = grep_molecule::Main(argc, argv);

  return rc;
}

