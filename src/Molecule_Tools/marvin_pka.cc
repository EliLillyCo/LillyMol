// Take the output from Chemaxon's cxcalc pka and merge with
// a smiles. 
// Typical cxcalc output looks like
// id	apKa1	apKa2	bpKa1	bpKa2	atoms
// 1	9.74		10.78	-5.48	7,14,7
// 2	16.22		-1.65	-4.41	9,4,2
// 3			10.27		13
// 4	13.09	14.22	8.46	-3.03	30,28,26,6
// 5	19.59		-4.84	-7.28	17,5,18
// 6	3.85	15.11	-1.35		16,13,12
// 7	3.88	15.20	-1.36		17,13,12
// 8	7.93		0.10	-7.72	6,5,8
// 9					

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace marvin_pka {

using std::cerr;

void
Usage(int rc) {
  cerr << "Generates formally charged molecules based on a Chemaxon pKa calculation\n";
  cerr << "\n";
  cerr << " -M <fname>      output from 'cxcalc pka file.smi'\n";
  cerr << " -P <pH>         pH value to use\n";
  cerr << " -d <delta>      only assign charges that are > <delta> from pH\n";
  cerr << " -q <qry>        only write molecules that match <qry>\n";
  cerr << " -s <smarts>     only write molecules that match <smarts>\n";
  cerr << " -f q1,q2,...    only write molecules that have net formal charge values of q1 or q2 or...\n";
  cerr << " -S <fname>      output file name (default smiles to stdout)\n";
  cerr << " -o <type>       output type (default smiles)\n";
  cerr << " -l              reduce to largest fragment - messes with statistics generated\n";
  cerr << " -c              remove chirality\n";
  cerr << " -g <...>        chemical standardisation options - use with caution\n";
  cerr << " -v              verbose output\n";

  ::exit(rc);
}

constexpr float kMissing = std::numeric_limits<float>::max();

// Generated from parsing a cxcalc output record.
// There are varying numbers of tokens in an output from cxcalc.
struct CXResult {
  int id = -1;

  // Some of these may be missing.
  float calc[4];

  // For each item in `calc`, the corresponding atom number.
  int atoms[4];

  CXResult();

  int Build(const const_IWSubstring& buffer,
            const resizable_array<int>& columns);
  float apKa1() const {
    return calc[0];
  }
  float apKa2() const {
    return calc[1];
  }
  float bpKa1() const {
    return calc[2];
  }
  float bpKa2() const {
    return calc[3];
  }
};

CXResult::CXResult() {
  for (int i = 0; i < 4; ++i) {
    calc[i] = kMissing;
    atoms[i] = INVALID_ATOM_NUMBER;
  }
}

int
GetColumnStop(const const_IWSubstring& buffer,
              const resizable_array<int>& columns,
              int col) {
  if (col != columns.number_elements()) {
    return columns[col] + 1;
  }

  return buffer.nchars() - 1;
}

// this is complicated by the presence of missing values in
// the pK columns. Therefore we use the `atom` array to keep
// track of which columns have a numeric value, and use that
// to manage the atoms list.
// id	apKa1	apKa2	bpKa1	bpKa2	atoms
// 2	16.22		-1.65	-4.41	9,4,2
// 3			10.27		13
// 4	13.09	14.22	8.46	-3.03	30,28,26,6
int
CXResult::Build(const const_IWSubstring& buffer,
                 const resizable_array<int>& columns) {
  // Since the pKa values have gaps, we need to know which columns
  // the atom numbers refer to.
  resizable_array<int> atom;
  atom.resize(4);
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword_single_delimiter(token, i, '\t'); ++col) {
    if (col == 0) {
      if (! token.numeric_value(id) || id < 1) {
        cerr << "CXResult::Build:invalid identifier " << buffer << '\n';
        return 0;
      }
      continue;
    }

    if (token.empty()) {
      continue;
    }

    if (col >= 1 && col <= 4) {
      if (! token.numeric_value(calc[col - 1])) {
        cerr << "CXResult::Build:invalid numeric '" << token << "'\n";
        cerr << buffer << '\n';
        return 0;
      }
      atom << (col - 1);
      continue;
    }

    if (col == 5) {
      int i = 0;
      const_IWSubstring s;
      for (int ndx = 0; token.nextword(s, i, ','); ++ndx) {
        int atom_number;
        if (! s.numeric_value(atom_number)) {
          cerr << "CXResult::Build:invalid atom number " << token << "\n";
          cerr << buffer << '\n';
          return 0;
        }
        atoms[atom[ndx]] = (atom_number - 1);
        cerr << "ndx " << ndx << " is atom " << atoms[atom[ndx]] << '\n';
      }
    }
  }

  return 1;
}

class CXResults {
  private:
    // The starting columns for the data in the input file;

    resizable_array<int> _columns;

    // the sie of _results.
    int _nresults;

    // One for each line in the input file from cxcalc.
    CXResult * _results;

    // Need something to return if we are asked for a non
    // existent result. Should just abort instead...
    CXResult _empty_result;

  // private functions
    int Build(iwstring_data_source& input);

  public:
    CXResults();
    ~CXResults();

    int Build(const char * fname);

    int number_results() const {
      return _nresults;
    }

    const CXResult& Result(int i) const;
};

CXResults::CXResults() {
  _nresults = 0;
  _results = nullptr;
}

CXResults::~CXResults() {
  if (_results != nullptr) {
    delete [] _results;
  }
}

const CXResult&
CXResults::Result(int ndx) const {
  if (ndx < _nresults) {
    return _results[ndx];
  }
  cerr << "CXResult::Result:do not have " << ndx << " results\n";
  return _empty_result;
}

int
CXResults::Build(const char * fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "CXResults::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
CXResults::Build(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "CXResults::Build:cannot read header record\n";
    return 0;
  }

  _nresults = input.records_remaining();
  if (_nresults == 0) {
    cerr << "CXResults::Build:no data\n";
    return 0;
  }

  _results = new CXResult[_nresults];

  for (int i = 0; input.next_record(buffer), i < _nresults; ++i) {
    if (! _results[i].Build(buffer, _columns)) {
      cerr << "CXResults::Build:error processing line " << i << " '" << buffer << "'\n";
      return 0;
    }
  }

  return _nresults;
}


class Options {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    CXResults _marvin_results;

    float _pH;

    // Only assign if the computed pKa is > _delta from _pH
    float _delta;

    int _any_charge;
    int _zwitterion;
    int _first_positive_charge;
    int _second_positive_charge;
    int _first_negative_charge;
    int _second_negative_charge;

    // Only molecules matching this query are written.
    resizable_array_p<Substructure_Query> _only_write;

    resizable_array<int> _ok_net_formal_charge;

    int _rejected_for_net_formal_charge;

    int _molecules_written;

    // As a sanity check we can see what kinds of atoms get assigned + and - charges
    IW_STL_Hash_Map_int _negative_smarts;
    IW_STL_Hash_Map_int _positive_smarts;

  // private functions
    int Process(Molecule& m, const CXResult& result,
                Molecule_Output_Object& output);
    int MaybeSetFormalCharge(Molecule& m,
                     atom_number_t zatom,
                     int formal_charge);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);
    // Because reducing to the largest fragment perturbs atom numbers
    // that cannot be done in Preprocess.
    // Note that subsequently losing fragments meses up all the statistics
    // about what is assigned, because some of those may be in the fragments.
    // Too messy to worry about.
    int PostProcess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    void ReadAnotherMolecule() {
      ++_molecules_read;
    }

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }

    int Process(Molecule& m, Molecule_Output_Object& output);
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _pH = 7.0;
  _delta = 0.0;

  _rejected_for_net_formal_charge = 0;

  _any_charge = 0;
  _zwitterion = 0;
  _first_positive_charge = 0;
  _second_positive_charge = 0;
  _first_negative_charge = 0;
  _second_negative_charge = 0;

  _molecules_written = 0;
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

  if (! cl.option_present('M')) {
    cerr << "Options::Initialise:must specify file of cxcalc results via the -M option\n";
    Usage(1);
  } else {
    const char * fname = cl.option_value('M');
    if (! _marvin_results.Build(fname)) {
      cerr << "Options::Initialise:cannot read marvin results '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _marvin_results.number_results() << " cxcalc results from '" << fname << "'\n";
    }
  }

  if (cl.option_present('P')) {
    if (! cl.value('P', _pH) || _pH < 0.0 ||  _pH > 14.0) {
      cerr << "Options::Initialise:invalid pKa \n";
      return 0;
    }
    if (_verbose) {
      cerr << "Assignments assume pH " << _pH << '\n';
    }
  }

  if (cl.option_present('d')) {
    if (! cl.value('d', _delta) || _delta < 0.0)  {
      cerr << "Options::Initialise:invalid -d option\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only assign a formal charge if > " << _delta << " from pH\n";
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _only_write, _verbose, 'q')) {
      cerr << "Options::Initialise:cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _only_write << qry.release();
    }
  }

  if (cl.option_present('f')) {
    int f;
    for (int i = 0; cl.value('f', f, i); ++i) {
      _ok_net_formal_charge << f;
    }

    if (_verbose) {
      cerr << "Will only write molecules with net formal charge of";
      for (int f : _ok_net_formal_charge) {
        cerr << ' ' << f;
      }
      cerr << '\n';
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

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _any_charge << " molecules assigned any charge\n";
  output << _first_negative_charge << " assigned first negative charge\n";
  output << _second_negative_charge << " assigned second negative charge\n";
  output << _first_positive_charge << " assigned first positive charge\n";
  output << _second_positive_charge << " assigned second positive charge\n";
  output << _zwitterion << " assigned both +ve and -ve charges\n";
  output << _rejected_for_net_formal_charge << " rejected for invalid net formal charge\n";
  output << _molecules_written << " molecules written\n";

  output << "The following atom types were assigned positive charges\n";
  for (const auto& [smarts, count] : _positive_smarts) {
    output << " + " << smarts << ' ' << count << '\n';
  }
  output << "The following atom types were assigned negative charges\n";
  for (const auto& [smarts, count] : _negative_smarts) {
    output << " - " << smarts << ' ' << count << '\n';
  }

  return 1;
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
Options::PostProcess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 Molecule_Output_Object& output) {
  if (_molecules_read > _marvin_results.number_results()) {
    cerr << "Options::Process:not enough data " << _molecules_read << " molecules read "
         << _marvin_results.number_results() << " results\n";
    return 0;
  }

  const CXResult& result = _marvin_results.Result(_molecules_read - 1);

  return Process(m, result, output);
}

int
HitsAnyOf(Molecule& m,
          resizable_array_p<Substructure_Query>& queries) {
  Molecule_to_Match target(&m);
  for (Substructure_Query * q : queries) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }
  return 0;
}

int
Options::MaybeSetFormalCharge(Molecule& m,
                     atom_number_t zatom,
                     int formal_charge) {
  if (! m.ok_atom_number(zatom)) {
    cerr << "Invalid atom number " << zatom << " in " << m.smiles() << '\n';
    return 0;
  }

  if (formal_charge < 0) {
    _negative_smarts[m.smarts_equivalent_for_atom(zatom)]++;
  } else {
    _positive_smarts[m.smarts_equivalent_for_atom(zatom)]++;
  }

  m.set_formal_charge(zatom, formal_charge);

  return 1;
}

int
Options::Process(Molecule& m,
                 const CXResult& result,
                 Molecule_Output_Object& output) {
  float apKa1 = result.apKa1();
  float apKa2 = result.apKa2();
  float bpKa1 = result.bpKa1();
  float bpKa2 = result.bpKa2();

  int any_charge = 0;
  if (apKa1 != kMissing && apKa1 < (_pH - _delta)) {
    if (!MaybeSetFormalCharge(m, result.atoms[0], -1)) {
      return 0;
    }
    ++_first_negative_charge;
    ++any_charge;
  }
  if (apKa2 != kMissing && apKa2 < (_pH - _delta)) {
    if (!MaybeSetFormalCharge(m, result.atoms[1], -1)) {
      return 0;
    }
    ++_second_negative_charge;
  }
  if (bpKa1 != kMissing && bpKa1 > (_pH + _delta)) {
    if (!MaybeSetFormalCharge(m, result.atoms[2], 1)) {
      return 0;
    }
    ++_first_positive_charge;
    ++any_charge;
  }
  if (bpKa2 != kMissing && bpKa2 > (_pH + _delta)) {
    if (!MaybeSetFormalCharge(m, result.atoms[3], 1)) {
      return 0;
    }
    ++_second_positive_charge;
  }

  if (any_charge) {
    ++_any_charge;
    if (any_charge == 2) {
      ++_zwitterion;
    }
  }
  
  if (_only_write.empty()) {
  } else if (HitsAnyOf(m, _only_write)) {
  } else {
    return 1;
  }

  if (_ok_net_formal_charge.empty()) {
  } else if (! _ok_net_formal_charge.contains(m.net_formal_charge())) {
    ++_rejected_for_net_formal_charge;
    return 1;
  }

  ++_molecules_written;

  PostProcess(m);

  return output.write(m);
}

int
MarvinPka(Options& options,
                Molecule& m,
                Molecule_Output_Object& output) {
  return options.Process(m, output);
}

int
MarvinPka(Options& options,
                data_source_and_type<Molecule>& input,
                Molecule_Output_Object& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    options.ReadAnotherMolecule();

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! MarvinPka(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
MarvinPka(Options& options,
             const char * fname,
             Molecule_Output_Object& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "MarvinPka:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return MarvinPka(options, input, output);
}

int
MarvinPka(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:M:P:d:q:s:S:o:f:");

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

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (! output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output types\n";
      return 1;
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);
  }

  if (cl.option_present('S')) {
    const IWString stem = cl.string_value('S');
    if (! output.new_stem(stem)) {
      cerr << "Cannot set output file stem\n";
      return 1;
    }
  } else {
    output.new_stem("-");
  }

  for (const char * fname : cl) {
    if (! MarvinPka(options, fname, output)) {
      cerr << "MarvinPka::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.do_flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace marvin_pka

int
main(int argc, char ** argv) {

  int rc = marvin_pka::MarvinPka(argc, argv);

  return rc;
}
