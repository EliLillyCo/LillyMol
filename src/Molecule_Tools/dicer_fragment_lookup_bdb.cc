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
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_fragment_lookup_bdb {

using std::cerr;

struct SmilesNatoms {
  IWString smiles;
  int natoms;
  SmilesNatoms(const const_IWSubstring& s, int n) : smiles(s), natoms(n) {
  }
};

class DicerFragmentLookupImpl {
  private:
    int _verbose;

    resizable_array_p<Db> _database;

    // Molecules are considered OK if the count of the
    // rarest fragment is > _min_count_needed;
    // by default it is zero.
    int _min_count_needed;

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

    // The output separator in our output files.
    char _sep = ' ';

    // The values used to initialise per_natom_min.
    int maxint = std::numeric_limits<int>::max();

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
  _molecules_processed = 0;
  _molecules_with_no_fragments = 0;
  _min_count_needed = 0;
  _above_min_count = 0;
  _rarest_fragment_on_separate_line = 0;
  _write_per_natoms_data = 0;
  _looking_for_highest_precedent = 0;
}

void
DisplayDashJOptions(std::ostream& output) {
  output << " -J def      default conditions\n";
  output << " -J          other qualifiers will be added\n";

  exit(0);
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
    }
    _database << db.release();
  }

  if (_verbose) {
    cerr << "Opened " << _database.size() << " dicer precedent databases\n";
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', _min_count_needed)) {
      cerr << "DicerFragmentLookupImpl::Initialise:the -c option must be a whole number\n";
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

  if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    if (!_stream_for_below_threshold.open(fname.null_terminated_chars())) {
      cerr << "DicerFragmentLookupImpl::Initialise:cannot open -X file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Molecules with fragments below " << _min_count_needed << " written to '" << fname << "'\n";
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


  // For each fragment size, the lowest count found.
  extending_resizable_array<int> per_natom_min(maxint);

  // We could speed up this loop by breaking once we have a
  // fragment that falls below _min_count_needed.

  const int nfrag = fragments.number_elements();
  int lowest_count = maxint;
  int index_with_lowest_count = -1;

  for (int i = 0; i < nfrag; ++i) {
    int count = Lookup(fragments[i]->smiles);

    if (count < per_natom_min[fragments[i]->natoms]) {
      per_natom_min[fragments[i]->natoms] = count;
    }
    if (count < lowest_count) {
      lowest_count = count;
      index_with_lowest_count = i;
    }
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

  for (const SmilesNatoms* iter : fragments) {
    std::optional<dicer_data::DicerFragment> maybe_frag = FirstMatch(iter->smiles);
    if (! maybe_frag) {
      continue;
    }

    if (static_cast<int>(maybe_frag->n()) < _min_count_needed) {
      continue;
    }

    output << smiles << _sep << name << _sep;
    if (_rarest_fragment_on_separate_line) {
      output << '\n';
    } else {
      output << _sep;
    }
    output << iter->smiles << _sep << maybe_frag->n() << _sep
           << maybe_frag->par() << '\n';
    output.write_if_buffer_holds_more_than(4096);
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
  cerr << "Consumes the output from dicer and looks up fragments in a database\n";
  cerr << " -d <dbname>     name of database(s) containing dicer fragments\n";
  cerr << " -c <count>      minimum number of exemplars required\n";
  cerr << "                 enter a negative value and all molecules will pass\n";
  cerr << " -X <fname>      write failed molecules (below -c) to <fname>\n";
  cerr << " -J ...          look up the most common fragments, enter '-r help' for info\n";
  cerr << " -y              write the rarest fragment on a separate line (easier with vf -c 2)\n";
  cerr << " -h              write a header record (cannot combine with -y)\n";
  cerr << " -w              write per natoms count for failed molecules\n";
  cerr << " -v              verbose output\n";
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
    fragments_this_molecule << new SmilesNatoms(smiles, natoms);
  }

  if (!fragments_this_molecule.empty()) {
    dicer_fragment_lookup.ProcessMolecule(current_smiles, current_name,
                fragments_this_molecule, output);
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

  return DicerFragmentLookup(input, dicer_fragment_lookup, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:c:X:yhwJ:");
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
