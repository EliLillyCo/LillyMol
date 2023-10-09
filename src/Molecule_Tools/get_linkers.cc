// Extract linker groups from a set of molecules.
// Break molecules at rings and collate the resulting fragments.
// This is conceptually somewhat like the inverse of ring_extraction.
// One flaw with this approach is that we only consider linkers
// that are chain only. More interesting would be to include
// rings. The scaffolds tool could perhaps do that.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace get_linkers {

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
  cerr << "Extracts inter-ring linkers to DicerFragment proto\n";
  cerr << " -m <natoms>   discard linkers with fewer than <natoms> atoms\n";
  cerr << " -M <natoms>   discard linkers with more  than <natoms> atoms\n";
  cerr << " -P ...        atom typing specification\n";
  cerr << " -c            discard chirality\n";
  cerr << " -l            reduce to largest fragment\n";
  cerr << " -v            verbose output\n";
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

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    Atom_Typing_Specification _atype;

    // We can impose atom count limits on the size of
    // linkers processed.
    int _min_natoms;
    int _max_natoms;

    // This will ultimately be written.
    IW_STL_Hash_Map<IWString, dicer_data::DicerFragment> _seen;

    int _molecules_read = 0;

  // private functions
    int OkAtomCount(const Molecule& m) const;
    int AddToHash(Molecule& m, const IWString& parent_name);

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

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After fragments have been accumulated
    int WriteFragments(IWString& fname);
    int WriteFragments(IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _min_natoms = 0;
  _max_natoms = std::numeric_limits<int>::max();

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

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_natoms) || _min_natoms < 1) {
      cerr << "The minimum number of atoms (-m) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only write fragments with at least " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_natoms) || _max_natoms < _min_natoms) {
      cerr << "The maximum number of atoms (-M) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only write fragments with fewer than " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('P')) {
    const IWString a = cl.string_value('P');
    if (! _atype.build(a)) {
      cerr << "Invalid atom typing specification '" << a << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::OkAtomCount(const Molecule& m) const {
  const int matoms = m.natoms();
  if (matoms < _min_natoms) {
    return 0;
  }
  if (matoms > _max_natoms) {
    return 0;
  }
  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Generated " << _seen.size() << " linkers\n";

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

  if (m.nrings() == 0) {
    return 0;
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

// Return true if `zatom` is doubly bonded to a ring atom.
bool
DoublyBondedToToRing(Molecule& m,
                     atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (m.ring_bond_count(j) > 0) {
      return true;
    }
  }

  return false;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  m.compute_aromaticity_if_needed();

  const int nedges = m.nedges();
  const int matoms = m.natoms();

  std::unique_ptr<int[]> break_bond(new_int(nedges));
  std::unique_ptr<uint32_t[]> atype;
  if (_atype.active()) {
    atype = std::make_unique<uint32_t[]>(matoms);
    std::fill_n(atype.get(), matoms, 0);
    _atype.assign_atom_types(m, atype.get());
  }

  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (b->nrings() > 0) {
      continue;
    }

    if (b->is_triple_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    int rbc1 = m.ring_bond_count(a1);
    int rbc2 = m.ring_bond_count(a2);

    // Skip biphenyls, no linker
    if (b->is_single_bond() && rbc1 > 0 && rbc2 > 0) {
      continue;
    }

    // If we are a single bond attached to a ring, definitely break this bond.
    // Note this test is more complex than it needs to be because we know
    // that this bond is not in a ring, so only one of rbc1, rbc2 can be nonzero.
    if (b->is_single_bond() && 
        ((rbc1 > 0 && rbc2 == 0) || (rbc1 == 0 && rbc2 > 0))) {
      break_bond[i] = 1;
      continue;
    }

    // A non ring bond and neither atom is in a ring, ignore.
    if (b->is_single_bond() && rbc1 == 0 && rbc2 == 0) {
      continue;
    }

    // If either a1 or a2 is doubly bonded to a ring, we break this bond.
    if (b->is_single_bond() && (DoublyBondedToToRing(m, a1) || DoublyBondedToToRing(m, a2))) {
      break_bond[i] = 1;
      continue;
    }
  }

  int rc = 0;
  for (int i = nedges - 1; i >= 0; --i) {
    if (! break_bond[i]) {
      continue;
    }

    const Bond* b = m.bondi(i);
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (atype) {
      m.set_isotope(a1, atype[a2]);
      m.set_isotope(a2, atype[a1]);
    } else {
      m.set_isotope(a1, m.atomic_number(a2));
      m.set_isotope(a2, m.atomic_number(a1));
    }
    m.remove_bond(i);
    ++rc;
  }

  if (rc == 0) {
    return 1;
  }

  resizable_array_p<Molecule> components;
  m.create_components(components);
  for (Molecule* c : components) {
    if (! OkAtomCount(*c)) {
      continue;
    }
    if (c->number_isotopic_atoms() < 2) {
      continue;
    }
    if (c->nrings()) {
      continue;
    }

    AddToHash(*c, m.name());
  }

  // The IWString_and_File_Descriptor object needs to be flushed.
  output.write_if_buffer_holds_more_than(4192);

  return 1;
}

int
Options::AddToHash(Molecule& m,
                   const IWString& parent_name) {
  const IWString& usmi = m.unique_smiles();

  const auto iter = _seen.find(usmi);
  if (iter != _seen.end()) {
    auto existing_count = iter->second.n();
    iter->second.set_n(existing_count + 1);
    return 1;
  }

  dicer_data::DicerFragment proto;
  proto.set_smi(usmi.AsString());
  proto.set_par(parent_name.AsString());

  proto.set_n(1);

  if (_atype.active()) {
    proto.set_iso(dicer_data::ATYPE);
  } else {
    proto.set_iso(dicer_data::Z);
  }
  _seen.emplace(std::make_pair(usmi, std::move(proto)));

  return 1;
}

int
Options::WriteFragments(IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "Options::WriteFragments:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteFragments(output);
}

int
Options::WriteFragments(IWString_and_File_Descriptor& output) {

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const auto& [_, proto] : _seen) {
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteFragments:cannot write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
GetLinkers(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
GetLinkers(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! GetLinkers(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetLinkers(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "GetLinkers:cannot determine input type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "GetLinkers:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return GetLinkers(options, input, output);
}

int
GetLinkers(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:P:m:M:");

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

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SDF;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GetLinkers(options, fname, input_type, output)) {
      cerr << "GetLinkers::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.WriteFragments(output);
  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace get_linkers

int
main(int argc, char ** argv) {

  int rc = get_linkers::GetLinkers(argc, argv);

  return rc;
}
