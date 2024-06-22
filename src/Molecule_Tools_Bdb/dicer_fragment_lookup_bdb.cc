// Consumes the output from dicer and looks up the
// smiles in a Berkeley DB database containing
// Dicer.DicerFragment protos.

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

#include "google/protobuf/text_format.h"

#include "db_cxx.h"


#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_fragment_lookup_bdb {

using std::cerr;

struct SmilesNatoms {
  IWString smiles;
  int natoms;
  SmilesNatoms(const const_IWSubstring& s, int n) : smiles(s), natoms(n) {
  }
  SmilesNatoms(const std::string& s, int n) : smiles(s), natoms(n) {
  }
};

// As we look up fragments, we can do a substructure search to see how the
// fragment is embedded in the parent molecule. We set an isotope to the
// total number of times an atom is matched by a fragment.
class AtomCoverage {
  private:
    Molecule _m;

    // for each atom, the number of fragments matching the atom.
    int* _times_matched;

    // The number of atoms in the largest fragment that matches.
    int _largest_fragment_match;

    Molecule_to_Query_Specifications _mqs;

    std::unique_ptr<Molecule_to_Match> _target;

    Fraction_as_String _fraction;

  public:
    AtomCoverage(const IWString& smiles, const IWString& name);
    ~AtomCoverage();

    int UpdateCoverage(const IWString& smiles);
    int UpdateCoverage(Molecule& frag);

    int WriteCoverageSmiles(IWString_and_File_Descriptor& output);
};

AtomCoverage::AtomCoverage(const IWString& smiles, const IWString& name) {
  if (! _m.build_from_smiles(smiles)) {
    cerr << "AtomCoverage:invalid smiles '" << smiles << "'\n";
    return;
  }

  _m.set_name(name);

  _times_matched = new_int(_m.natoms());
  
  _mqs.set_substituents_only_at_isotopic_atoms(1);

  _target = std::make_unique<Molecule_to_Match>(&_m);

  _largest_fragment_match = 0;

  _fraction.initialise(0.0f, 1.0f, 3);
}

AtomCoverage::~AtomCoverage() {
  delete [] _times_matched;
}

int
AtomCoverage::UpdateCoverage(const IWString& smiles) {
  if (smiles.empty()) {
    cerr << "Ignoring empty smiles\n";
    return 0;
  }

  Molecule frag;
  if (! frag.build_from_smiles(smiles)) {
    cerr << "AtomCoverage::UpdateCoverage:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return UpdateCoverage(frag);
}

int
AtomCoverage::UpdateCoverage(Molecule& frag) {
  Substructure_Query qry;
  if (! qry.create_from_molecule(frag, _mqs)) {
    cerr << "AtomCoverage::UpdateCoverage:cannot convert " << frag.smiles() << " to query\n";
    return 0;
  }

  qry.set_find_unique_embeddings_only(1);

  Substructure_Results sresults;
  if (qry.substructure_search(*_target, sresults) == 0) {
    cerr << "AtomCoverage:UpdateCoverage:no query match for " << frag.smiles() <<
            " in " << _m.smiles() << '\n';
    return 0;
  }

  if (frag.natoms() > _largest_fragment_match) {
    _largest_fragment_match = frag.natoms();
  }

  for (const Set_of_Atoms* e : sresults.embeddings()) {
    e->increment_vector(_times_matched, 1);
  }

  return sresults.number_embeddings();
}

int
AtomCoverage::WriteCoverageSmiles(IWString_and_File_Descriptor& output) {
  _m.set_isotopes(_times_matched);

  const int niso = _m.number_isotopic_atoms();
  float f = iwmisc::Fraction<float>(niso, _m.natoms());

  static constexpr char kSep = ' ';

  output << _m.smiles() << kSep << _m.name() << kSep << niso << kSep;
  _fraction.append_number(output, f);
  output << kSep << _largest_fragment_match << kSep;

  f = iwmisc::Fraction<float>(_largest_fragment_match, _m.natoms());
  _fraction.append_number(output, f);


  return 1;
}

class DicerFragmentLookupImpl {
  private:
    int _verbose;

    resizable_array_p<Db> _database;

    // Molecules are considered OK if the count of the
    // rarest fragment is > _min_count_needed;
    // by default it is zero.
    int _min_count_needed;

    int _input_is_textproto;

    // the number of molecules passed to DoOutput.
    int _molecules_processed;

    // Some molecules will not have any fragments that pass
    // the rules for geneation.
    int _molecules_with_no_fragments;

    // The number of molecules above our threshold.
    int _above_min_count;

    // Optionally, molecules with fragments below the threshold
    // can be written.
    IWString_and_File_Descriptor _stream_for_below_threshold;

    // To make things easy to view with vf, we can write the rarest
    // fragment on a separate line.
    int _rarest_fragment_on_separate_line;

    // It can be interesting to record the atom count of
    // the rarest fragment.
    extending_resizable_array<int> _rarest_fragent_size;

    // If instructed, we can write per shell count data to cerr.
    int _write_per_natoms_data;

    // The initial implementation of this tool was looking
    // for the rarest fragment in the input. But quickly we
    // identified a need for looking for the most precedented
    // fragments.
    int _looking_for_highest_precedent;

    // We might have input that contains fragments smaller or larger
    // than what might be stored in the database, so have the ability
    // to filter fragments as they are read in.
    int _min_fragment_size;
    int _max_fragment_size;

    int _write_header_record;

    // The output separator in our output files.
    char _sep = ' ';

    // The values used to initialise per_natom_min.
    const int maxint = std::numeric_limits<int>::max();

    // As fragments are read, we can do a substructure search and increment the
    // number of times each atom is hit by one of the fragments fetched from the
    // database.
    int _compute_atom_coverage;

    IWString_and_File_Descriptor _stream_for_atom_coverage;

    // Private functions

    // Write result to `output`.
    int Write(const IWString& smiles,
                const IWString& name,
                const SmilesNatoms& rarest_frag,
                int count,
                IWString_and_File_Descriptor& output);
    int WritePerNatomsData(const IWString& smiles,
                const IWString&name, 
                const extending_resizable_array<int>& per_natom_min, 
                std::ostream& output) const;

    std::optional<dicer_data::DicerFragment> 
        FirstMatch(const const_IWSubstring& smiles);

    int HighestPrecedent(const IWString& smiles,
                                const IWString& name,
                                resizable_array_p<SmilesNatoms>& fragments,
                                IWString_and_File_Descriptor& output);

  public:
    DicerFragmentLookupImpl();

    int Initialise(Command_Line& cl);

    int input_is_textproto() const {
      return _input_is_textproto;
    }

    int OkAtomCount(int natoms) const;

    // Lookup `smiles` in each database and return the total
    // number of instances.
    int Lookup(const const_IWSubstring& smiles);

    int OkCount(int count) const {
      return count > _min_count_needed;
    }

    // A fragment has been looked up and a count determined.
    // Checks to see if this is to be written to `output` as
    // a passing molecule, or maybe to _stream_for_below_threshold
      // if below threshold.
    int DoOutput(const const_IWSubstring& smiles,
                 const const_IWSubstring& id,
                 int count,
                 const IWString& lowest_count_fragment,
                 IWString_and_File_Descriptor& output);

    int ProcessMolecule(const IWString& smiles,
                                const IWString& name,
                                resizable_array_p<SmilesNatoms>& fragments,
                                IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

DicerFragmentLookupImpl::DicerFragmentLookupImpl() {
  _verbose = 0;
  _input_is_textproto = 0;
  _write_header_record = 0;
  _molecules_processed = 0;
  _molecules_with_no_fragments = 0;
  _min_count_needed = 0;
  _above_min_count = 0;
  _min_fragment_size = 0;
  _max_fragment_size = 0;
  _rarest_fragment_on_separate_line = 0;
  _write_per_natoms_data = 0;
  _looking_for_highest_precedent = 0;
  _compute_atom_coverage = 0;
}

void
DisplayDashJOptions(std::ostream& output) {
  output << " -J def      default conditions\n";
  output << " -J          other qualifiers will be added\n";

  exit(0);
}

void
DisplayDashYOptions(std::ostream& output) {
  cerr << R"(
 -Y header              write a header record
 -Y vf2                 write the rarest fragment on a separate line, useful with vf
 -Y peratom             for each molecule, write number of failures as a function of size.
)";

  ::exit(0);
}

int
DicerFragmentLookupImpl::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');
  if (! cl.option_present('d')) {
    cerr << "Must specify one or more databases via the -d option\n";
    return 0;
  }

  int flags = 0;
  DBTYPE dbtype = DB_UNKNOWN;
  int mode = 0;
  DbEnv* env = NULL;

  IWString dbname;
  for (int i = 0; cl.value('d', dbname, i); ++i) {
    std::unique_ptr<Db> db = std::make_unique<Db>(env, DB_CXX_NO_EXCEPTIONS);
    if (int rc = db->open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);
        rc != 0) {
      cerr << "Cannot open '" << dbname << "' ";
      db->err(rc, "");
      return 0;
    }
    _database << db.release();
  }

  if (_verbose) {
    cerr << "Opened " << _database.size() << " dicer precedent databases\n";
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _min_count_needed)) {
      cerr << "DicerFragmentLookupImpl::Initialise:the -p option must be a whole number\n";
      return 0;
    }

    if (_min_count_needed < 0) {
      cerr << "Negative threshold, all molecules will pass\n";
    } else if (_verbose) {
      cerr << "Threshold for frequency set to " << _min_count_needed << '\n';
    }
  }

  if (cl.option_present('y')) {
    _rarest_fragment_on_separate_line = 1;
    if (_verbose) {
      cerr << "Rarest fragment written on separate line\n";
    }
  }

  if (cl.option_present('w')) {
    _write_per_natoms_data = 1;
    if (_verbose) {
      cerr << "Will write per natoms data for failed molecules\n";
    }
  }

  if (cl.option_present('J')) {
    IWString opt;
    for (int i = 0; cl.value('J', opt, i); ++i) {
      if (opt == "def") {
      } else if (opt.starts_with("")) {
      } else if (opt == "help") {
        DisplayDashJOptions(cerr);
      } else {
        cerr << "Unrecognised -J qualifier '" << opt << "'\n";
        DisplayDashJOptions(cerr);
      }
    }

    _looking_for_highest_precedent = 1;
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', _min_fragment_size) || _min_fragment_size < 1) {
      cerr << "The minimum fragment size (-c) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will discard fragments with fewer than " << _min_fragment_size << " atoms\n";
    }
  }

  if (cl.option_present('C')) {
    if (! cl.value('C', _max_fragment_size) || _max_fragment_size < 1) {
      cerr << "The maximum fragment size (-C) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will discard fragments with more than " << _max_fragment_size << " atoms\n";
    }
  }

  if (_min_fragment_size > _max_fragment_size) {
    cerr << "Inconsistent min " << _min_fragment_size << " and max " << 
            _max_fragment_size << " fragment size specifications\n";
    return 0;
  }

  if (cl.option_present('t')) {
    _input_is_textproto = 1;
    if (_verbose) {
      cerr << "Input is textproto\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "header") {
        _write_header_record = 1;
        if (_verbose) {
          cerr << "Will write a header record\n";
        }
      } else if (y == "vf2") {
        _rarest_fragment_on_separate_line = 1;
        if (_verbose) {
          cerr << "Rarest fragment written on separate line\n";
        }
      } else if (y == "peratom") {
        _write_per_natoms_data = 1;
        if (_verbose) {
          cerr << "Will write per natoms data for failed molecules\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    fname.EnsureEndsWith(".smi");
    if (!_stream_for_below_threshold.open(fname.null_terminated_chars())) {
      cerr << "DicerFragmentLookupImpl::Initialise:cannot open -X file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Molecules with fragments below " << _min_count_needed << " written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('F')) {
    _compute_atom_coverage = 1;

    IWString fname = cl.string_value('F');

    fname.EnsureEndsWith(".smi");

    if (!_stream_for_atom_coverage.open(fname.null_terminated_chars())) {
      cerr << "DicerFragmentLookupImpl::Initialise:cannot open -F file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Atom coverage data written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
DicerFragmentLookupImpl::Lookup(const const_IWSubstring& smiles) {
  Dbt dkey;
  dkey.set_data((void*)smiles.data());
  dkey.set_size(smiles.length());

  Dbt fromdb;
  int total_count = 0;
  for (Db* db : _database) {
    if (0 != db->get(NULL, &dkey, &fromdb, 0)) {
      // cerr << "'no match for '" << smiles << '\n';
      continue;
    }
    google::protobuf::io::ArrayInputStream input(fromdb.get_data(), fromdb.get_size());
    dicer_data::DicerFragment frag;
    if (! google::protobuf::TextFormat::Parse(&input, &frag)) {
      cerr << "DicerFragmentLookupImpl::Lookup:invalid db contents '";
      cerr.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
      cerr << "'\n";
      return 0;
    }
    // cerr << "smiles '" << smiles << " got " << frag.ShortDebugString() << '\n';
    total_count += frag.n();
  }

  return total_count;
}

// Return the proto of the first match to `smiles` from any
// of our databases.
std::optional<dicer_data::DicerFragment>
DicerFragmentLookupImpl::FirstMatch(const const_IWSubstring& smiles) {
  Dbt dkey;
  dkey.set_data((void*)smiles.data());
  dkey.set_size(smiles.length());

  Dbt fromdb;
  for (Db* db : _database) {
    if (0 != db->get(NULL, &dkey, &fromdb, 0)) {
      continue;
    }
    google::protobuf::io::ArrayInputStream input(fromdb.get_data(), fromdb.get_size());
    dicer_data::DicerFragment frag;
    if (! google::protobuf::TextFormat::Parse(&input, &frag)) {
      cerr << "DicerFragmentLookupImpl::Lookup:invalid db contents '";
      cerr.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
      cerr << "'\n";
      return std::nullopt;
    }

#ifdef DEBUG_FIRST_MATCH
    const_IWSubstring tmp((const char*)fromdb.get_data(), fromdb.get_size());
    cerr << "Key '" << smiles << "'\n";
    cerr << "From " << tmp << " build\n";
    cerr << frag.ShortDebugString() << '\n';
#endif
    return frag;
  }

  return std::nullopt;
}

int
DicerFragmentLookupImpl::DoOutput(const const_IWSubstring& smiles,
                                  const const_IWSubstring& id,
                                  int count,
                                  const IWString& lowest_count_fragment,
                                  IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  if (count > _min_count_needed) {
    output << smiles << _sep << id << _sep << count;
    if (_rarest_fragment_on_separate_line) {
      output << '\n';
    } else {
      output << _sep;
    }
    output << lowest_count_fragment << '\n';
    output.write_if_buffer_holds_more_than(4096);
    ++_above_min_count;
    return 1;
  }

  if (_stream_for_below_threshold.active()) {
    _stream_for_below_threshold << smiles << _sep <<id << _sep
         << count;
    if (_rarest_fragment_on_separate_line) {
      _stream_for_below_threshold << '\n';
    } else {
      _stream_for_below_threshold << _sep;
    }
    _stream_for_below_threshold << lowest_count_fragment << '\n';
    _stream_for_below_threshold.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
WriteAtomCoverage(AtomCoverage& atom_coverage, 
                  IWString_and_File_Descriptor& output) {
  atom_coverage.WriteCoverageSmiles(output);
  output << '\n';
    output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
DicerFragmentLookupImpl::OkAtomCount(int natoms) const {
  if (_max_fragment_size > 0 && natoms > _max_fragment_size) {
    return 0;
  }

  if (_min_fragment_size > 0 && natoms < _min_fragment_size) {
    return 0;
  }

  return 1;
}

int
DicerFragmentLookupImpl::ProcessMolecule(const IWString& smiles,
                                const IWString& name,
                                resizable_array_p<SmilesNatoms>& fragments,
                                IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  if (fragments.empty()) {
    ++_molecules_with_no_fragments;
    output << smiles << _sep << name << '\n';
    return 1;
  }

  if (_looking_for_highest_precedent) {
    return HighestPrecedent(smiles, name, fragments, output);
  }

  // Sort fragments so the smallest are first.
  fragments.iwqsort_lambda([] (const SmilesNatoms* f1, const SmilesNatoms* f2) {
    if (f1->natoms < f2->natoms) {
      return -1;
    }  else if (f1->natoms > f2->natoms) {
      return 1;
    } else {
      return 0;
    }
  });

  std::unique_ptr<AtomCoverage> atom_coverage;
  if (_compute_atom_coverage) {
    atom_coverage.reset(new AtomCoverage(smiles, name));
  }

  // For each fragment size, the lowest count found.
  extending_resizable_array<int> per_natom_min(maxint);

  // We could speed up this loop by breaking once we have a
  // fragment that falls below _min_count_needed.

  const int nfrag = fragments.number_elements();
  int lowest_count = maxint;
  int index_with_lowest_count = -1;

  for (int i = 0; i < nfrag; ++i) {

    int count = Lookup(fragments[i]->smiles);
    if (count == 0) {
      continue;
    }

    if (atom_coverage) {
      atom_coverage->UpdateCoverage(fragments[i]->smiles);
    }

    if (count < per_natom_min[fragments[i]->natoms]) {
      per_natom_min[fragments[i]->natoms] = count;
    }
    if (count < lowest_count) {
      lowest_count = count;
      index_with_lowest_count = i;
    }
  }

  if (atom_coverage) {
    WriteAtomCoverage(*atom_coverage, _stream_for_atom_coverage);
  }

  if (lowest_count > _min_count_needed) {
    ++_above_min_count;
    return Write(smiles, name, *fragments[index_with_lowest_count],
                 lowest_count, output);
  }

  if (_write_per_natoms_data) {
    WritePerNatomsData(smiles, name, per_natom_min, cerr);
  }

  if (! _stream_for_below_threshold.active()) {
    return 1;
  }

  return Write(smiles, name, *fragments[index_with_lowest_count],
               lowest_count, _stream_for_below_threshold);
}

int
DicerFragmentLookupImpl::HighestPrecedent(const IWString& smiles,
                                const IWString& name,
                                resizable_array_p<SmilesNatoms>& fragments,
                                IWString_and_File_Descriptor& output) {

  // Sort fragments so the largest are first.
  fragments.iwqsort_lambda([] (const SmilesNatoms* f1, const SmilesNatoms* f2) {
    if (f1->natoms < f2->natoms) {
      return 1;
    }  else if (f1->natoms > f2->natoms) {
      return -1;
    } else {
      return 0;
    }
  });

  std::unique_ptr<AtomCoverage> atom_coverage;
  if (_compute_atom_coverage) {
    atom_coverage.reset(new AtomCoverage(smiles, name));
  }

  for (const SmilesNatoms* iter : fragments) {
    std::optional<dicer_data::DicerFragment> maybe_frag = FirstMatch(iter->smiles);
    if (! maybe_frag) {
      continue;
    }

    if (atom_coverage) {
      atom_coverage->UpdateCoverage(maybe_frag->smi());
    }

    if (static_cast<int>(maybe_frag->n()) < _min_count_needed) {
      continue;
    }

    output << smiles << _sep << name << _sep;
    if (_rarest_fragment_on_separate_line) {
      output << '\n';
    } 

    output << iter->smiles << _sep << maybe_frag->n() << _sep
           << maybe_frag->par() << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  if (atom_coverage) {
    WriteAtomCoverage(*atom_coverage, _stream_for_atom_coverage);
  }

  return 1;
}

int
DicerFragmentLookupImpl::Write(const IWString& smiles,
                const IWString& name,
                const SmilesNatoms& rarest_frag,
                int count,
                IWString_and_File_Descriptor& output) {
  output << smiles << _sep << name << _sep << count;
  if (_rarest_fragment_on_separate_line) {
    output << '\n';
  } else {
    output << _sep;
  }
  output << rarest_frag.smiles << _sep << rarest_frag.natoms << '\n';
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
DicerFragmentLookupImpl::WritePerNatomsData(const IWString& smiles,
                const IWString& name,
                const extending_resizable_array<int>& per_natom_min,
                std::ostream& output) const {

  output << smiles << _sep << name << '\n';
  for (int i = 0;  i < per_natom_min.number_elements(); ++i) {
    if (per_natom_min[i] < maxint) {
      output << smiles << _sep << name << _sep << i << _sep << per_natom_min[i] << '\n';
    }
  }

  return 1;
}

int
DicerFragmentLookupImpl::Report(std::ostream& output) const {
  output << "Read " << _molecules_processed << " molecules\n";
  output << _molecules_with_no_fragments << " molecules generated no fragments\n";
  output << _above_min_count << " above " << _min_count_needed << '\n';
  return 1;
}

int
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Consumes the output from dicer and looks up fragments in a BerkeleyDB database.
If all fragments in the molecule are in the database, subject to the support requirement (-p)
then the molecule is considered to have 'passed'.
For example

 dicer -B proto ... -S example.textproto example.smi
 dicer_database_lookup_bdb -F atom_coverage -d collection.bdb -t example.textproto

 -d <dbname>            name of database(s) containing dicer fragments - built by dicer2bdb
 -t                     input is textproto
 -c <natoms>            minimum fragment size to process - discards these fragments upon input
 -C <natoms>            maximum fragment size to process - discards these fragments upon input
 -p <count>             minimum number of exemplars required for a molecule to pass.
                        enter a negative value and all molecules will pass
 -X <fname>             write failed molecules (one or more fragments below -p) to <fname>
 -J ...                 look up the most common fragments, enter '-J help' for info
 -Y ...                 various other options, enter '-Y help' for details.
 -F <fname>             write fragment coverage (labelled smiles) to <fname>
 -v                     verbose output
)";
  // clang-format on

  ::exit(rc);
}

int
SmilesAndId(const_IWSubstring& buffer,
            const_IWSubstring& smiles,
            const_IWSubstring& id) {
  int i = 0;
  if (! buffer.nextword(smiles, i) ||
      ! buffer.nextword(id, i)) {
    return 0;
  }
  if (buffer.empty() || smiles.empty()) {
    return 0;
  }

  return 1;
}

int
DicerFragmentLookup(iwstring_data_source& input,
                    DicerFragmentLookupImpl& dicer_fragment_lookup,
                    IWString_and_File_Descriptor& output) {
  std::unique_ptr<RE2> b_equals = std::make_unique<RE2>(" B=[0-9][0-9]*$");

  IWString current_smiles;
  IWString current_name;
  resizable_array_p<SmilesNatoms> fragments_this_molecule;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    const_IWSubstring smiles, id;
    if (! SmilesAndId(buffer, smiles, id)) {
      cerr << "Invalid input '" << buffer << "'\n";
      return 0;
    }

    // New molecule encountered, input ends with " B=\d+"
    if (iwre2::RE2PartialMatch(buffer, *b_equals)) {
      if (! current_smiles.empty()) {
        dicer_fragment_lookup.ProcessMolecule(current_smiles, current_name,
                  fragments_this_molecule, output);
      }
      current_smiles = smiles;
      current_name = id;
      fragments_this_molecule.resize_keep_storage(0);
      continue;
    }

    const int natoms = count_atoms_in_smiles(smiles);
    if (! dicer_fragment_lookup.OkAtomCount(natoms)) {
      continue;
    }
    fragments_this_molecule << new SmilesNatoms(smiles, natoms);
  }

  if (!fragments_this_molecule.empty()) {
    dicer_fragment_lookup.ProcessMolecule(current_smiles, current_name,
                fragments_this_molecule, output);
  }

  return 1;
}

int
DicerFragmentLookupTextproto(const dicer_data::DicedMolecule& proto,
                             DicerFragmentLookupImpl& dicer_fragment_lookup,
                             IWString_and_File_Descriptor& output) {
  if (proto.fragment_size() == 0) {
    return 1;

  }
  resizable_array_p<SmilesNatoms> fragments_this_molecule;
  fragments_this_molecule.reserve(proto.fragment_size());

  for (const dicer_data::DicerFragment& frag : proto.fragment()) {
    if (! dicer_fragment_lookup.OkAtomCount(frag.nat())) {
      continue;
    }

    fragments_this_molecule << new SmilesNatoms(frag.smi(), frag.nat());
  }

  const IWString current_smiles(proto.smiles());
  const IWString current_name(proto.name());

  return dicer_fragment_lookup.ProcessMolecule(current_smiles, current_name,
        fragments_this_molecule, output);
}

int
DicerFragmentLookupTextproto(const const_IWSubstring& buffer,
                             DicerFragmentLookupImpl& dicer_fragment_lookup,
                             IWString_and_File_Descriptor& output) {
  dicer_data::DicedMolecule mol;
  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());
  if (! google::protobuf::TextFormat::Parse(&input, &mol)) {
    cerr << "DicerFragmentLookupTextproto:cannot parse textproto\n";
    return 0;
  }

  return DicerFragmentLookupTextproto(mol, dicer_fragment_lookup, output);
}

int
DicerFragmentLookupTextproto(iwstring_data_source& input,
                             DicerFragmentLookupImpl& dicer_fragment_lookup,
                             IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! DicerFragmentLookupTextproto(buffer, dicer_fragment_lookup, output)) {
      cerr << "DicerFragmentLookupTextproto:invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
DicerFragmentLookup(const char * fname,
                    DicerFragmentLookupImpl& dicer_fragment_lookup,
                    IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DicerFragmentLookup:cannot open '" << fname << "'\n";
    return 0;
  }

  if (dicer_fragment_lookup.input_is_textproto()) {
    return DicerFragmentLookupTextproto(input, dicer_fragment_lookup, output);
  }

  return DicerFragmentLookup(input, dicer_fragment_lookup, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:c:C:X:yhwJ:tp:F:Y:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  DicerFragmentLookupImpl dicer_fragment_lookup;
  if (! dicer_fragment_lookup.Initialise(cl)) {
    cerr << "Cannot initialise dicer fragment lookup\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('h')) {
    if (cl.option_present('y')) {
      cerr << "Cannot combine -h and -y options\n";
    } else {
      output << "Smiles id count frag\n";
    }
  }

  for (const char * fname: cl) {
    if (! DicerFragmentLookup(fname, dicer_fragment_lookup, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    dicer_fragment_lookup.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragment_lookup_bdb

int
main(int argc, char ** argv) {

  int rc = dicer_fragment_lookup_bdb::Main(argc, argv);

  return rc;
}
