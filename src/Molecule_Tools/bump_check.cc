// Bump checks on 3D molecules.

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace bump_check {

using std::cerr;

void
Usage(int rc) {
  cerr << "Performs 3D bump checks\n";
  cerr << " -t <distance>     the distance below which a bump check occurs\n";
  cerr << " -x                do not consider atom pairs bonded to a common atom\n";
  cerr << " -h                do NOT consider Hydrogen atoms\n";
  cerr << " -R                do NOT consider atoms in the same ring system\n";
  cerr << " -r                report all distances shorter than <distance>\n";
  cerr << " -s <smarts>       specify atom pairs to be checked (otherwise all non bonded atoms compared)\n";
  cerr << " -q <qury>         specify atom pairs as queries\n";
  cerr << " -S <fname>        write passing molecules to <fname>\n";
  cerr << " -U <fname>        write failing molecules to <fname>\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class CheckTooClose {
  private:
    int _verbose;

    int _molecules_read;

    int _ok_molecules;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // The distance for the bump check test.
    float _distance;

    int _skip_adjacent;

    int _skip_atoms_in_same_ring_system;

    int _exclude_hydrogen;

    resizable_array_p<Substructure_Query> _query;

    // Write molecules that do NOT violate bump checks.
    Molecule_Output_Object _output;

    // Write molecules that violate bump checks.
    Molecule_Output_Object _rejected;

    int _report_too_short;

    Accumulator<double> _acc_dist;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions.

    int IdentifyAtomsToConsider(Molecule& m, int* process_atom);

  public:
    CheckTooClose();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

CheckTooClose::CheckTooClose() {
  _verbose = 0;
  _molecules_read = 0;
  _ok_molecules = 0;
  _skip_adjacent = 0;
  _skip_atoms_in_same_ring_system = 0;
  _report_too_short = 0;
  _exclude_hydrogen = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
}

int
CheckTooClose::Initialise(Command_Line& cl) {
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

  if (! cl.option_present('t')) {
    cerr << "Must specify a distance threshold via the -t option\n";
    return 0;
  }

  if (cl.option_present('t')) {
    if (! cl.value('t', _distance) || _distance <= 0.0f) {
      cerr << "Invalid distance threshold (-t)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Distance threshold " << _distance << '\n';
    }
  }

  if (cl.option_present('x')) {
    _skip_adjacent = 1;
    if (_verbose) {
      cerr << "Will not consider adjacent atoms (bonded to common atom)\n";
    }
  }

  if (cl.option_present('r')) {
    _report_too_short = 1;
    if (_verbose) {
      cerr << "Will report all atom pairs closer than " << _distance << '\n';
    }
  }

  if (cl.option_present('R')) {
    _skip_atoms_in_same_ring_system = 1;
    if (_verbose) {
      cerr << "Will skip atoms in the same ring system\n";
    }
  }

  if (cl.option_present('h')) {
    _exclude_hydrogen = 1;
    if (_verbose) {
      cerr << "Will not consider Hydrogen atoms\n";
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

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }
      _query << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "Cannot process command line queries (-q)\n";
      return 0;
    }
  }

  for (Substructure_Query* q : _query) {
    q->set_find_unique_embeddings_only(1);
    q->set_perceive_symmetry_equivalent_matches(0);
  }

  set_append_coordinates_after_each_atom(1);

  if (cl.option_present('S')) {
    _output.add_output_type(FILE_TYPE_SMI);
    _output.add_output_type(FILE_TYPE_SDF);
    IWString fname = cl.string_value('S');
    if (_output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 0;
    }

    if (! _output.new_stem(fname)) {
      cerr << "Cannot open file for successful molecules (-S) '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "OK molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('U')) {
    _rejected.add_output_type(FILE_TYPE_SMI);
    _rejected.add_output_type(FILE_TYPE_SDF);
    IWString fname = cl.string_value('U');
    if (_rejected.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 0;
    }

    if (! _rejected.new_stem(fname)) {
      cerr << "Cannot open file for successful molecules (-U) '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Rejected molecules written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
CheckTooClose::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.natoms() == 0) {
    return 0;
  }

  return 1;
}

// Return true of `a1` and `a2` are both bonded to some atom.
int
JoinedToSameAtom(const Molecule& m,
                 atom_number_t a1,
                 atom_number_t a2) {
  const Set_of_Atoms c1 = m.connections(a1);

  for (const Bond* b : m.atom(a2)) {
    atom_number_t j = b->other(a2);
    if (c1.contains(j)) {
      return 1;
    }
  }
  return 0;
}

int
CheckTooClose::IdentifyAtomsToConsider(Molecule& m,
                        int* process_atom) {
  std::fill_n(process_atom, m.natoms(), 0);

  Molecule_to_Match target(&m);
  int got_match = 0;
  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }
    sresults.each_embedding_set_vector(process_atom, 1);
    ++got_match;
  }

  if (! got_match) {
    cerr << "CheckTooClose::IdentifyAtomsToConsider:no match to any of " << _query.size() <<
            " queries to " << m.name() << '\n';
    return 0;
  }

  return 1;
}

int
CheckTooClose::Process(Molecule& m) {
  ++_molecules_read;

  // if we are reporting distances, make the smarts as informative as possible.
  if (_report_too_short) {
    m.compute_aromaticity_if_needed();
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> process_atom(new_int(matoms, 1));
  if (! _query.empty()) {
    if (IdentifyAtomsToConsider(m, process_atom.get())) {
    } else {
      cerr << "Cannot identify atoms to process\n";
      return 0;
    }
  }

  int too_close = 0;

  for (int i = 0; i < matoms; ++i) {
    if (! process_atom[i]) {
      continue;
    }

    if (_exclude_hydrogen && m.atomic_number(i) == 1) {
      continue;
    }
    for (int j = i + 1; j < matoms; ++j) {
      if (! process_atom[j]) {
        continue;
      }

      if (_exclude_hydrogen && m.atomic_number(j) == 1) {
        continue;
      }

      if (m.are_bonded(i, j)) {
        continue;
      }

      if (_skip_atoms_in_same_ring_system && m.in_same_ring_system(i, j)) {
        continue;
      }

      if (_skip_adjacent && JoinedToSameAtom(m, i, j)) {
        continue;
      }

      const float d = m.distance_between_atoms(i, j);
      if (_verbose) {
        _acc_dist.extra(d);
      }

      if (d > _distance) {
        continue;
      }

      if (_report_too_short) {
        cerr << m.name() << " atoms " << i << ' ' << m.smarts_equivalent_for_atom(i) <<
                ' ' << j << ' ' << m.smarts_equivalent_for_atom(j) << " dist " << d << '\n';
      }
      ++too_close;
    }
  }

  if (too_close == 0) {
    ++_ok_molecules;

    if (_output.active()) {
      return _output.write(m);
    }

    return 1;
  }

  if (_rejected.active()) {
    IWString new_name(m.name());
    new_name << ' ' << too_close;
    m.set_name(new_name);
    return _rejected.write(m);
  }

  return 1;
}

int
CheckTooClose::Report(std::ostream& output) const {
  output << "CheckTooClose:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << _ok_molecules << " molecules did not violate distance " << _distance <<
            " fraction " << iwmisc::Fraction<float>(_ok_molecules, _molecules_read) << '\n';
  if (_acc_dist.n() > 0) {
    cerr << "Distances btw " << _acc_dist.minval() << " and " << _acc_dist.maxval() << " ave " << _acc_dist.average() << '\n';
  }

  return 1;
}

int
BumpCheck(CheckTooClose& check_too_close,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return check_too_close.Process(m);
}

int
BumpCheck(CheckTooClose& check_too_close,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! check_too_close.Preprocess(*m)) {
      return 0;
    }

    if (! BumpCheck(check_too_close, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
BumpCheck(CheckTooClose& check_too_close,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = check_too_close.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return BumpCheck(check_too_close, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcs:q:S:U:t:xrhR");
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

  CheckTooClose check_too_close;
  if (! check_too_close.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! BumpCheck(check_too_close, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    check_too_close.Report(cerr);
  }

  return 0;
}

}  // namespace bump_check

int
main(int argc, char** argv) {
  int rc = bump_check::Main(argc, argv);

  return rc;
}
