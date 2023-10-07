// Remove atoms from a molecule based on 3d coordinate matching.

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace excise_3d_overlap {

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
  cerr << "Performs some task on a set of molecules.\n";
  cerr << " -F <fname>   file containing the fragment\n";
  cerr << " -T <atol>    absolute tolerance for heavy atom matches\n";
  cerr << " -t <atol>    absolute tolerance for Hydrogen matches\n";
  cerr << " -X H         remove all Hydrogen atoms on input\n";
  cerr << " -S <fname>   write results to <fname>\n";
  cerr << " -v          verbose output\n";

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

    int _remove_hydrogens = 0;

    Chemical_Standardisation _chemical_standardisation;

    // One or more fragments read via the -F option.
    resizable_array_p<Molecule> _fragment;

    // The absolute difference in coordinates that is OK for a match.
    float _atol;

    // The absolute difference in coordinates when matching Hydrogen atoms.
    float _hydrogen_atol;

    // When we have a fragment that has > 1 connection, we enumerate all
    // possible labels. Query 0 will start placing isotopes at _delta, delta + 1...
    // and query 1 will start placing isotopes at (2*_delta), (2*_delta+1) ...
    int _delta;

    int _molecules_read = 0;

  // Private functions.
    int CloseEnough(const Atom& a1, const Atom& a2) const;
    int FragmentMatches(const Molecule& m,
                         const Molecule& frag,
                         int* matched,
                         int flag);
    int Process(Molecule& m,
                 const int * matched,
                 Molecule_Output_Object& output);
    int MakeVariants(const Molecule& m,
                      const int* matched,
                      int flag,
                      const Set_of_Atoms& adjacent,
                      Molecule_Output_Object& output);
    int PlaceOneIsotope(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent,
                      Molecule_Output_Object& output);
    int PlaceTwoIsotopes(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent1,
                      atom_number_t adjacent2,
                      Molecule_Output_Object& output);
    int PlaceThreeIsotopes(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent1,
                      atom_number_t adjacent2,
                      atom_number_t adjacent3,
                      Molecule_Output_Object& output);
    int MaybeProcessUnmatchedHydrogens(const Molecule& m,
                         const Molecule& frag,
                         int* matched,
                         int flag,
                         const int* xref,
                         const Set_of_Atoms& unmatched_hydrogens);

  public:
    Options();

    int verbose() const {
      return _verbose;
    }

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Process(Molecule& mol, Molecule_Output_Object& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _remove_hydrogens = 0;
  _molecules_read = 0;

  _atol = 1.0e-04;
  _hydrogen_atol = 0.1;

  _delta = 10;
}

std::unique_ptr<Molecule>
ReadFragment(IWString& fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_SDF, fname);
  if (! input.good()) {
    cerr << "ReadFragment:cannot open '" << fname << "'\n";
    return std::unique_ptr<Molecule>();
  }

  Molecule * m = input.next_molecule();
  if (m == nullptr) {
    cerr << "ReadFragment:no molecule\n";
    return std::unique_ptr<Molecule>();
  }

  return std::unique_ptr<Molecule>(m);
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  // We could make this more robust/general with an Element_Match object.
  if (cl.option_present('X')) {
    const_IWSubstring x = cl.option_value('X');
    if (x == 'H') {
      _remove_hydrogens = 1;
    } else {
      cerr << "Unrecognised -X qualifier '" << x << "'\n";
      return 0;
    }
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

  if (! cl.option_present('F')) {
    cerr << "Must specify one or more fragments via the -F option\n";
    Usage(1);
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      std::unique_ptr<Molecule> m = ReadFragment(fname);
      if (! m) {
        cerr << "Cannot read '" << fname << "'\n";
        return 0;
      }

      _fragment << m.release();
    }
  }

  if (cl.option_present('T')) {
    if (! cl.value('T', _atol) || _atol < 0.0) {
      cerr << "Invalid heavy atom absolute tolerance (-T)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Heavy atoms matced if within " << _atol << '\n';
    }
  }

  if (cl.option_present('t')) {
    if (! cl.value('T', _hydrogen_atol) || _hydrogen_atol < 0.0) {
      cerr << "Invalid hydrogen absolute tolerance (-t)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Hydrogens atoms matced if within " << _hydrogen_atol << '\n';
    }
  }

  if (_remove_hydrogens) {
    for (Molecule* m : _fragment) {
      m->remove_all(1);
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  // Other information about what has happened.

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

  if (_remove_hydrogens) {
    m.remove_all(1);
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }


  return 1;
}

int
AtomsCloseEnough(const Atom& a1,
            const Atom& a2,
            float atol) {
  return a1.distance(a2) <= atol;
}

int
Options::CloseEnough(const Atom& a1,
                     const Atom& a2) const {
  if (a1.atomic_number() == 1 && a2.atomic_number() == 1) {
    return AtomsCloseEnough(a1, a2, _hydrogen_atol);
  }

  return AtomsCloseEnough(a1, a2, _atol);
}

// For each atom in `frag`, identify the atom in `m` that is
// in the same spatial location.
int
Options::FragmentMatches(const Molecule& m,
                         const Molecule& frag,
                         int* matched,
                         int flag) {
  const int matoms = m.natoms();
  const int fatoms = frag.natoms();

  // Hydrogen atoms may not match, but they can perhaps be processed.
  Set_of_Atoms unmatched_hydrogens;

  // for each fragment atom, the atom number in `m` that matched.
  std::unique_ptr<int[]> matched_atom(new_int(fatoms, INVALID_ATOM_NUMBER));

  for (int i = 0; i < fatoms; ++i) {
    const Atom& f = frag.atom(i);

    int got_match = 0;
    for (int j = 0; j < matoms; ++j) {
      if (matched[j]) {
        continue;
      }

#ifdef DEBUG_EXCISE_3D_OVERLAP
      cerr << "btw " << i << ' ' <<  frag.smarts_equivalent_for_atom(i) << " and " << j <<
              m.smarts_equivalent_for_atom(j) << " dist " << f.distance(m.atom(j)) << '\n';
#endif
      if (! CloseEnough(f, m.atom(j))) {
        continue;
      }
      matched[j] = flag;
      got_match = 1;
      matched_atom[i] = j;
      if (_verbose > 1) {
        cerr << "match frag atom " << i << ' ' << frag.smarts_equivalent_for_atom(i) << " and " <<
             j << ' ' << m.smarts_equivalent_for_atom(j) << "\n";
      }
      break;
    }

    if (got_match) {  // great
      continue;
    }

    if (frag.atomic_number(i) == 1) {
      unmatched_hydrogens << i;
      continue;
    }

    if (! got_match) {
      cerr << "Options::FragmentMatches:spatial match for atom at " << f << " in " << m.name() << '\n';
      return 0;
    }
  }

  if (unmatched_hydrogens.empty()) {
    ;
  } else if (! MaybeProcessUnmatchedHydrogens(m, frag, matched, flag,
               matched_atom.get(), unmatched_hydrogens)) {
    cerr << "Cannot resolve unmatched Hydrogen atoms\n";
    return 0;
  }

  return 1;
}

// We need to sort a list of atoms and attached hydrogens.
struct AtomAndH {
  atom_number_t atom;
  atom_number_t hydrogen;
};

// We are attempting to match the atoms in `frag` with atoms in `m`.
// All the heavy atoms match, but some Hydrogen atoms do not.
// `xref` is a mapping from an atom number in `frag` to the atom that matched
// that atom in `m`.
int
Options::MaybeProcessUnmatchedHydrogens(const Molecule& m,
                         const Molecule& frag,
                         int* matched,
                         int flag,
                         const int* xref,
                         const Set_of_Atoms& unmatched_hydrogens) {
  // For each atom in unmatched_hydrogens, gather the heavy atom.
  Set_of_Atoms adjacent;
  for (atom_number_t h : unmatched_hydrogens) {
    // This should not happen.
    if (frag.ncon(h) == 0) {
      cerr << "Unattached atom in fragment '" << frag.name() << "'\n";
      return 0;
    }

    atom_number_t j = frag.other(h, 0);
    // If the attached atom was not matched to anything in `m`, we cannot proceed
    if (xref[j] == INVALID_ATOM_NUMBER) {
      return 0;
    }

    adjacent << j;
  }

  assert (adjacent.size() == unmatched_hydrogens.size());

  // We need to group things by the atom to which the hydrogen atoms are attached.
  std::unique_ptr<AtomAndH[]> atom_and_h = std::make_unique<AtomAndH[]>(adjacent.size());

  const int n = adjacent.size();
  for (int i = 0; i < n; ++i) {
    atom_and_h[i].atom = adjacent[i];
    atom_and_h[i].hydrogen = unmatched_hydrogens[i];
  }

  std::sort(atom_and_h.get(), atom_and_h.get() + n, [](const AtomAndH& ah1, const AtomAndH& ah2) {
    return ah1.atom < ah2.atom;
  });

#define DEBUG_UNMATCHED_HYDROGENS
#ifdef DEBUG_UNMATCHED_HYDROGENS
  for (int i = 0; i < n; ++i) {
    cerr << "Hydrogen " << atom_and_h[i].hydrogen << " attached to " << atom_and_h[i].atom << " matched to " << xref[atom_and_h[i].atom] << '\n';
  }
#endif

  // We now have a fragment atom and all its unmatched Hydrogens. We then need to look at `m`
  // and examine the matched atom, and see if there are matching H atoms.
  return 1;
}

int
Options::Process(Molecule& m,
                 Molecule_Output_Object& output) {
  ++_molecules_read;

  const int matoms = m.natoms();
  std::unique_ptr<int[]> matched(new_int(matoms));

  const int number_fragments = _fragment.number_elements();

  for (int i = 0; i < number_fragments; ++i) {
    if (! FragmentMatches(m, *_fragment[i], matched.get(), i + 1)) {
      cerr << "Options::Process:no geometry match to " << _fragment[i]->name() << " in " << m.name() << '\n';
      return 0;
    }
  }

  return Options::Process(m, matched.get(), output);
}

// For each atom in `m` for which matched[i] == flag, scan the
// adjacent atoms that have matched[j] == 0, and add them to 
// `destination`.
int
IdentifyUnmatchedJoining(const Molecule& m,
                         const int* matched,
                         int flag,
                         Set_of_Atoms& destination) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (matched[i] == 0) {
      continue;
    }
    if (matched[i] != flag) {
      continue;
    }

    for (const Bond* b : m.atom(i)) {
      atom_number_t j = b->other(i);
      if (matched[j] == 0) {
        destination << j;
      }
    }
  }

  return destination.number_elements();
}

int
Options::Process(Molecule& m,
                 const int * matched,
                 Molecule_Output_Object& output) {
  const int number_fragments = _fragment.number_elements();

  for (int i = 0; i < number_fragments; ++i) {
    Set_of_Atoms adjacent;
    IdentifyUnmatchedJoining(m, matched, i + 1, adjacent);
    if (adjacent.empty()) {
      cerr << "Options::Process:no unmatched adjacent atoms for fragment " << i << ' ' << _fragment[i]->name() << '\n';
    }
    MakeVariants(m, matched, i + 1, adjacent, output);
  }

  return 1;
}

int
Options::MakeVariants(const Molecule& m,
                      const int* matched,
                      int flag,
                      const Set_of_Atoms& adjacent,
                      Molecule_Output_Object& output) {
  if (adjacent.size() == 1) {
    return PlaceOneIsotope(m, matched, flag, adjacent[0], output);
  }

  if (adjacent.size() == 2) {
    PlaceTwoIsotopes(m, matched, flag, adjacent[0], adjacent[1], output);
    PlaceTwoIsotopes(m, matched, flag, adjacent[1], adjacent[0], output);
    return 1;
  }

  if (adjacent.size() == 3) {
    PlaceThreeIsotopes(m, matched, flag, adjacent[0], adjacent[1], adjacent[2], output);
    PlaceThreeIsotopes(m, matched, flag, adjacent[1], adjacent[0], adjacent[2], output);
    PlaceThreeIsotopes(m, matched, flag, adjacent[0], adjacent[2], adjacent[1], output);
    PlaceThreeIsotopes(m, matched, flag, adjacent[1], adjacent[2], adjacent[0], output);
    PlaceThreeIsotopes(m, matched, flag, adjacent[2], adjacent[0], adjacent[1], output);
    PlaceThreeIsotopes(m, matched, flag, adjacent[2], adjacent[1], adjacent[0], output);
    return 1;
  }

  cerr << "Options::MakeVariants:case of " << adjacent.size() << " attachments not covered, see Ian\n";
  return 0;
}


// Create a copy of `m`, place isotope `flag` on `adjacent`,
// remove all atoms for which matched[i]==flag and write.
int
Options::PlaceOneIsotope(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent,
                      Molecule_Output_Object& output) {
  Molecule mcopy(m);
  mcopy.set_isotope(adjacent, flag);
  mcopy.remove_atoms(matched, flag);
  mcopy.set_name(m.name());
  return output.write(mcopy);
}

// Create a copy of `m`, place isotope `flag` on `adjacent1` and
// `adjacent2.  Remove all atoms for which matched[i]==flag and write.
int
Options::PlaceTwoIsotopes(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent1,
                      atom_number_t adjacent2,
                      Molecule_Output_Object& output) {
  Molecule mcopy(m);
  mcopy.set_isotope(adjacent1, flag);
  mcopy.set_isotope(adjacent2, flag);
  mcopy.remove_atoms(matched, flag);
  mcopy.set_name(m.name());
  return output.write(mcopy);
}

// Create a copy of `m`, place isotope `flag` on `adjacent1`, `adjacent2,
// and `adjacent3.  Remove all atoms for which matched[i]==flag and write.
int
Options::PlaceThreeIsotopes(const Molecule& m,
                      const int* matched,
                      int flag,
                      atom_number_t adjacent1,
                      atom_number_t adjacent2,
                      atom_number_t adjacent3,
                      Molecule_Output_Object& output) {
  Molecule mcopy(m);
  mcopy.set_isotope(adjacent1, flag);
  mcopy.set_isotope(adjacent2, flag);
  mcopy.set_isotope(adjacent3, flag);
  mcopy.remove_atoms(matched, flag);
  mcopy.set_name(m.name());
  return output.write(mcopy);
}

int
Excise3DOverlap(Options& options,
                Molecule& m,
                Molecule_Output_Object& output) {
  return options.Process(m, output);
}

int
Excise3DOverlap(Options& options,
                data_source_and_type<Molecule>& input,
                Molecule_Output_Object& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! Excise3DOverlap(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Excise3DOverlap(Options& options,
             const char * fname,
             FileType input_type,
             Molecule_Output_Object& output) {

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "Excise3DOverlap:cannot determine input type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Excise3DOverlap:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return Excise3DOverlap(options, input, output);
}

int
Excise3DOverlap(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:o:S:F:t:T:X:");

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
    input_type = FILE_TYPE_SDF;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Molecule_Output_Object output;
  if (! cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    Usage(1);
  }

  if (cl.option_present('o')) {
    if (! output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 1;
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 1;
    }

    if (! output.new_stem(fname)) {
      cerr << "Cannot initialise output '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Results written to '" << fname << "'\n";
    }
  }

  for (const char * fname : cl) {
    if (! Excise3DOverlap(options, fname, input_type, output)) {
      cerr << "Excise3DOverlap::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace excise_3d_overlap

int
main(int argc, char ** argv) {

  int rc = excise_3d_overlap::Excise3DOverlap(argc, argv);

  return rc;
}
