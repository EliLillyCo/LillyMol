// Molecule Filtering Tool.
// Designed as a first level filter for large collections.

#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/molecule_filter_lib.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/molecule_filter.pb.h"
#else
#include "molecule_filter.pb.h"
#endif

namespace molecule_filter {

using std::cerr;

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
  cerr << R"(Filters molecules based on easy to compute molecular properties.
Designed for rapidly triaging large collections of molecules in the most
computationally efficient way possible.
 -F <fname>     textproto describing constraints on filter - required.
 -c             remove chirality
 -l             reduce to largest fragment
 -B <fname>     write rejected molecules to <fname>
 -u             with utilities, write smiles id overall_utility instead of descriptor table
 -v             verbose output
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

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    int _write_smiles_with_utility = 0;

    Chemical_Standardisation _chemical_standardisation;

    uint64_t _molecules_read = 0;
    uint64_t _molecules_passed = 0;
    uint64_t _molecules_with_explicit_hydrogens = 0;
    uint64_t _molecules_with_isotopes_rehydrogenated = 0;

    molecule_filter_data::Requirements _requirements;

    molecule_filter_lib::MoleculeFilter _filter;

    IWString_and_File_Descriptor _reject_stream;

    // Values accumulated based on rejections
    uint64_t _too_few_atoms = 0;
    uint64_t _too_many_atoms = 0;
    uint64_t _too_few_rings = 0;
    uint64_t _too_many_rings = 0;
    uint64_t _too_few_heteroatoms = 0;
    uint64_t _too_many_heteroatoms = 0;
    uint64_t _min_heteroatom_fraction = 0;
    uint64_t _max_heteroatom_fraction = 0;
    uint64_t _too_few_aromatic_rings = 0;
    uint64_t _too_many_aromatic_rings = 0;
    uint64_t _too_few_aliphatic_rings = 0;
    uint64_t _too_many_aliphatic_rings = 0;
    uint64_t _ring_system_too_small = 0;
    uint64_t _ring_system_too_large = 0;
    uint64_t _too_few_aromatic_rings_in_system = 0;
    uint64_t _too_many_aromatic_rings_in_system = 0;
    uint64_t _ring_too_small = 0;
    uint64_t _ring_too_large = 0;
    uint64_t _non_organic = 0;
    uint64_t _isotope = 0;
    uint64_t _too_few_rotbond = 0;
    uint64_t _too_many_rotbond = 0;
    uint64_t _low_tpsa = 0;
    uint64_t _high_tpsa = 0;
    uint64_t _low_xlogp = 0;
    uint64_t _high_xlogp = 0;
    uint64_t _low_alogp = 0;
    uint64_t _high_alogp = 0;
    uint64_t _low_amw = 0;
    uint64_t _high_amw = 0;
    uint64_t _too_few_hba_rdkit = 0;
    uint64_t _too_many_hba_rdkit = 0;
    uint64_t _too_few_hbd_rdkit = 0;
    uint64_t _too_many_hbd_rdkit = 0;
    uint64_t _too_few_hba = 0;
    uint64_t _too_many_hba = 0;
    uint64_t _too_few_hbd = 0;
    uint64_t _too_many_hbd = 0;
    uint64_t _too_few_halogens = 0;
    uint64_t _too_many_halogens = 0;
    uint64_t _too_short = 0;
    uint64_t _too_long = 0;
    uint64_t _too_few_csp3 = 0;
    uint64_t _too_many_csp3 = 0;
    uint64_t _fcsp3_too_low = 0;
    uint64_t _fcsp3_too_high = 0;
    uint64_t _aromdens_too_low = 0;
    uint64_t _aromdens_too_high = 0;
    uint64_t _too_few_chiral = 0;
    uint64_t _too_many_chiral = 0;
    uint64_t _too_few_fragments = 0;
    uint64_t _too_many_fragments = 0;

    uint64_t _matches_exclusion_smarts = 0;
    uint64_t _no_match_required_smarts = 0;

    // for use with parallel processing.
    off_t _seek_to;
    off_t _stop_at;

  // Private functions
    int Process(Molecule& m, int matoms, int nrings,
                std::vector<double>& utility_values, double& overall_utility);
    int MaybeWriteToRejectStream(const const_IWSubstring& line);
    int WriteUtilityValues(const std::vector<double>& utility_values,
                           double overall_utility,
                           const const_IWSubstring& smiles,
                           const const_IWSubstring& id,
                           IWString_and_File_Descriptor& output);

    void NoteRejection(molecule_filter_lib::RejectionReason rejection_reason);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    int Process(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output);

    int has_utilities() const {
      return _filter.has_utilities();
    }

    int write_smiles_with_utility() const {
      return _write_smiles_with_utility;
    }

    int WriteUtilityHeader(IWString_and_File_Descriptor& output) const;

    // If we are seeking and stopping.
    // If _seek_to is set, seek to that offset in `input`.
    int SeekIfNeeded(iwstring_data_source& input) const;
    // If _stop_at is set, return OK of the current offset is less than _stop
    int OkContinue(iwstring_data_source& input) const;

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _write_smiles_with_utility = 0;
  _molecules_read = 0;

  _seek_to = 0;
  _stop_at = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
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

  if (cl.option_present('u')) {
    _write_smiles_with_utility = 1;
    if (_verbose) {
      cerr << "Will write smiles plus overall utility in utility mode\n";
    }
  }

  if (! cl.option_present('F')) {
    cerr << "Options::Initialise:must specify Requirements proto via the -F option\n";
    return 0;
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    std::optional<molecule_filter_data::Requirements> maybe_proto =
        iwmisc::ReadTextProtoCommentsOK<molecule_filter_data::Requirements>(fname);
    if (! maybe_proto) {
      cerr << "Options::Initialise:cannot read textproto '" << fname << "'\n";
      return 0;
    }

    _requirements = std::move(*maybe_proto);
  }

  if (_requirements.required_smarts_size() ||
      _requirements.must_not_have_smarts_size()) {
    cerr << "Options::Initialise:smarts not implemented\n";
    return 0;
  }

  if (_remove_chirality &&
      ((_requirements.has_min_chiral() && _requirements.min_chiral() > 0) ||
       _requirements.has_max_chiral())) {
    cerr << "Options::Initialise:removing chirality has been specified (-c)\n";
    cerr << "But the config file contains a chiral-centre filter. Impossible\n";
    return 0;
  }

  if ((_requirements.has_min_distance() || _requirements.has_max_distance()) && ! _reduce_to_largest_fragment) {
    cerr << "Options::Initialise:distance specified, but largest fragment not selected (-l)\n";
    cerr << "Automatically enabling largest fragment selection\n";
    _reduce_to_largest_fragment = 1;
  }

  if (! _filter.Build(_requirements)) {
    cerr << "Options::Initialise:cannot initialise MoleculeFilter\n";
    return 0;
  }

  if (cl.option_present('B')) {
    IWString fname = cl.string_value('B');
    fname.EnsureEndsWith(".smi");
    if (! _reject_stream.open(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot open rejection file '" << fname << "'\n";
      return 0;
    } 

    if (_verbose) {
      cerr << "Rejected molecules written to " << fname << "'\n";
    }
  }

  if (cl.option_present('i')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('i', s, i); ++i) {
      if (s.starts_with("seek=")) {
        s.remove_leading_chars(5);
        if (! s.numeric_value(_seek_to)) {
          cerr << "Invalid seek specification '" << s << "'\n";
          return 0;
        }
      } else if (s.starts_with("stop=")) {
        s.remove_leading_chars(5);
        if (! s.numeric_value(_stop_at)) {
          cerr << "Invalid stop specification '" << s << "'\n";
          return 0;
        }
      } else {
        cerr << "Unrecognised -i qualifier '" << s << "'\n";
        return 0;
      }
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules, passed " << _molecules_passed <<
                ' ' << iwmisc::Fraction<float>(_molecules_passed, _molecules_read) << '\n';
  if (_requirements.has_exclude_non_organic()) {
    output << _non_organic << " non organic " << '\n';
  }
  if (_requirements.has_exclude_isotopes()) {
    output << _isotope << " isotope " << '\n';
  }
  if (_requirements.has_min_natoms()) {
    output << _too_few_atoms << " too few atoms " << _requirements.min_natoms() << '\n';
  }
  if (_requirements.has_max_natoms()) {
    output << _too_many_atoms << " too many atoms " << _requirements.max_natoms() << '\n';
  }
  if (_requirements.has_min_nrings()) {
    output << _too_few_rings << " too few rings " << _requirements.min_nrings() << '\n';
  }
  if (_requirements.has_max_nrings()) {
    output << _too_many_rings << " too many rings " << _requirements.max_nrings() << '\n';
  }

  if (_requirements.has_min_heteroatom_count()) {
    output << _too_few_heteroatoms << " too few heteroatoms " << _requirements.min_heteroatom_count() << '\n';
  }
  if (_requirements.has_max_heteroatom_count()) {
    output << _too_many_heteroatoms << " too many heteroatoms " << _requirements.max_heteroatom_count() << '\n';
  }
  if (_requirements.has_min_heteroatom_fraction()) {
    output << _min_heteroatom_fraction << " min heteroatom fraction " << _requirements.min_heteroatom_fraction() << '\n';
  }
  if (_requirements.has_max_heteroatom_fraction()) {
    output << _max_heteroatom_fraction << " max heteroatom fraction " << _requirements.max_heteroatom_fraction() << '\n';
  }

  if (_requirements.has_min_aromatic_ring_count()) {
    output << _too_few_aromatic_rings << " too few aromatic rings " << _requirements.min_aromatic_ring_count() << '\n';
  }
  if (_requirements.has_max_aromatic_ring_count()) {
    output << _too_many_aromatic_rings << " too many aromatic rings " << _requirements.max_aromatic_ring_count() << '\n';
  }

  if (_requirements.has_min_aliphatic_ring_count()) {
    output << _too_few_aliphatic_rings << " too few aliphatic rings " << _requirements.min_aliphatic_ring_count() << '\n';
  }
  if (_requirements.has_max_aliphatic_ring_count()) {
    output << _too_many_aliphatic_rings << " too many aliphatic rings " << _requirements.max_aliphatic_ring_count() << '\n';
  }

  if (_requirements.has_min_ring_system_size()) {
    output << _ring_system_too_small << " ring systems too small " << _requirements.min_ring_system_size() << '\n';
  }
  if (_requirements.has_max_ring_system_size()) {
    output << _ring_system_too_large << " ring systems too large " << _requirements.max_ring_system_size() << '\n';
  }
  if (_requirements.has_min_aromatic_rings_in_system()) {
    output << _too_few_aromatic_rings_in_system << " too few aromatic rings in system " << _requirements.min_aromatic_rings_in_system() << '\n';
  }
  if (_requirements.has_max_aromatic_rings_in_system()) {
    output << _too_many_aromatic_rings_in_system << " too many aromatic rings in system " << _requirements.max_aromatic_rings_in_system() << '\n';
  }

  if (_requirements.has_min_largest_ring_size()) {
    output << _ring_too_small << " ring too small " << _requirements.min_largest_ring_size() << '\n';
  }
  if (_requirements.has_max_largest_ring_size()) {
    output << _ring_too_large << " ring too large " << _requirements.max_largest_ring_size() << '\n';
  }

  if (_requirements.has_min_rotatable_bonds()) {
    output << _too_few_rotbond << " too few rotatable bonds " << _requirements.min_rotatable_bonds() << '\n';
  }
  if (_requirements.has_max_rotatable_bonds()) {
    output << _too_many_rotbond << " too many rotatable bonds " << _requirements.max_rotatable_bonds() << '\n';
  }

  if (_requirements.has_min_tpsa()) {
    output << _low_tpsa << " low TPSA " << _requirements.min_tpsa() << '\n';
  }
  if (_requirements.has_max_tpsa()) {
    output << _high_tpsa << " high TPSA " << _requirements.max_tpsa() << '\n';
  }

  if (_molecules_with_explicit_hydrogens) {
    output << _molecules_with_explicit_hydrogens <<
              " molecules had explicit Hydrogens removed before computing properties\n";
  }

  if (_molecules_with_isotopes_rehydrogenated) {
    output << _molecules_with_isotopes_rehydrogenated <<
              " isotopic molecules had their Hydrogen count recomputed\n";
  }

  if (_requirements.has_min_hba_rdkit()) {
    output << _too_few_hba_rdkit << " too few HBA_RDKIT " << _requirements.min_hba_rdkit() << '\n';
  }
  if (_requirements.has_max_hba_rdkit()) {
    output << _too_many_hba_rdkit << " too many HBA_RDKIT " << _requirements.max_hba_rdkit() << '\n';
  }
  if (_requirements.has_min_hbd_rdkit()) {
    output << _too_few_hbd_rdkit << " too few HBD_RDKIT " << _requirements.min_hbd_rdkit() << '\n';
  }
  if (_requirements.has_max_hbd_rdkit()) {
    output << _too_many_hbd_rdkit << " too many HBD_RDKIT " << _requirements.max_hbd_rdkit() << '\n';
  }

  if (_requirements.has_min_amw()) {
    output << _low_amw << " low AMW " << _requirements.min_amw() << '\n';
  }
  if (_requirements.has_max_amw()) {
    output << _high_amw << " high AMW " << _requirements.max_amw() << '\n';
  }

  if (_requirements.has_min_alogp()) {
    output << _low_alogp << " low ALOGP " << _requirements.min_alogp() << '\n';
  }
  if (_requirements.has_max_alogp()) {
    output << _high_alogp << " high ALOGP " << _requirements.max_alogp() << '\n';
  }

  if (_requirements.has_min_xlogp()) {
    output << _low_xlogp << " low XLOGP " << _requirements.min_xlogp() << '\n';
  }
  if (_requirements.has_max_xlogp()) {
    output << _high_xlogp << " high XLOGP " << _requirements.max_xlogp() << '\n';
  }

  if (_requirements.has_min_hba()) {
    output << _too_few_hba << " too few HBA " << _requirements.min_hba() << '\n';
  }
  if (_requirements.has_max_hba()) {
    output << _too_many_hba << " too many HBA " << _requirements.max_hba() << '\n';
  }

  if (_requirements.has_min_hbd()) {
    output << _too_few_hbd << " too few HBD " << _requirements.min_hbd() << '\n';
  }
  if (_requirements.has_max_hbd()) {
    output << _too_many_hbd << " too many HBD " << _requirements.max_hbd() << '\n';
  }

  if (_requirements.has_min_halogen_count()) {
    output << _too_few_halogens << " too few halogens " << _requirements.min_halogen_count() << '\n';
  }
  if (_requirements.has_max_halogen_count()) {
    output << _too_many_halogens << " too many halogens " << _requirements.max_halogen_count() << '\n';
  }

  if (_requirements.has_min_distance()) {
    output << _too_short << " molecules too short " << _requirements.min_distance() << '\n';
  }
  if (_requirements.has_max_distance()) {
    output << _too_long << " molecules too long " << _requirements.max_distance() << '\n';
  }

  if (_requirements.has_min_sp3_carbon()) {
    output << _too_few_csp3 << " too few CSP3 " << _requirements.min_sp3_carbon() << '\n';
  }
  if (_requirements.has_max_sp3_carbon()) {
    output << _too_many_csp3 << " too many CSP3 " << _requirements.max_sp3_carbon() << '\n';
  }
  if (_requirements.has_min_sp3_carbon_fraction()) {
    output << _fcsp3_too_low << " fraction CSP3 too low " << _requirements.min_sp3_carbon_fraction() << '\n';
  }
  if (_requirements.has_max_sp3_carbon_fraction()) {
    output << _fcsp3_too_high << " fraction CSP3 too high " << _requirements.max_sp3_carbon_fraction() << '\n';
  }

  if (_requirements.has_min_aromatic_density()) {
    output << _aromdens_too_low << " aromatic density too low " << _requirements.min_aromatic_density() << '\n';
  }
  if (_requirements.has_max_aromatic_density()) {
    output << _aromdens_too_high << " aromatic density too high " << _requirements.max_aromatic_density() << '\n';
  }

  if (_requirements.has_min_chiral()) {
    output << _too_few_chiral << " too few chiral centres " << _requirements.min_chiral() << '\n';
  }
  if (_requirements.has_max_chiral()) {
    output << _too_many_chiral << " too many chiral centres " << _requirements.max_chiral() << '\n';
  }

  if (_requirements.has_min_number_fragments()) {
    output << _too_few_fragments << " too few fragments " << _requirements.min_number_fragments() << '\n';
  }
  if (_requirements.has_max_number_fragments()) {
    output << _too_many_fragments << " too many fragments " << _requirements.max_number_fragments() << '\n';
  }

  return 1;
}

int
Options::WriteUtilityHeader(IWString_and_File_Descriptor& output) const {
  if (! _filter.has_utilities() || _write_smiles_with_utility) {
    return 1;
  }

  output << "Name";
  for (int i = 0; i < _filter.number_utilities(); ++i) {
    output << ' ' << _filter.utility(i).name();
  }
  output << " overall_utility\n";

  return 1;
}

// LargestFragment used to be duplicated here, identical to the one in
// molecule_filter_lib. Use the library version.

// If chemical standardisation is in effect
int
Options::Process(const const_IWSubstring& line,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  const_IWSubstring smiles, id;
  line.split(smiles, ' ', id);

  bool smiles_changed = false;

  if (_requirements.has_min_number_fragments() || _requirements.has_max_number_fragments()) {
    const int nfrag = smiles.ccount('.') + 1;
    if (_requirements.has_min_number_fragments() && nfrag < _requirements.min_number_fragments()) {
      ++_too_few_fragments;
      MaybeWriteToRejectStream(line);
      return 0;
    }
    if (_requirements.has_max_number_fragments() && nfrag > _requirements.max_number_fragments()) {
      ++_too_many_fragments;
      MaybeWriteToRejectStream(line);
      return 0;
    }
  }

  // The atom counts below are HEAVY atom counts. The min_natoms/max_natoms
  // requirements are about heavy atoms, so an explicit Hydrogen must not
  // count towards them. Determining this costs nothing - the smiles scan that
  // counts atoms is happening anyway.
  const_IWSubstring largest_frag;
  int matoms = 0;
  int nrings = 0;
  int explicit_hydrogens = 0;
  if (_reduce_to_largest_fragment) {
    if (molecule_filter_lib::LargestFragment(smiles, largest_frag, matoms,
                                            nrings, explicit_hydrogens)) {
      smiles_changed = true;
    }
  } else {
    matoms = lillymol::count_atoms_in_smiles(smiles, nrings, explicit_hydrogens) -
             explicit_hydrogens;
    largest_frag = smiles;
  }

  if (matoms == 0) {
    return 0;
  }

  if (_requirements.has_min_natoms() && matoms < _requirements.min_natoms()) {
    ++_too_few_atoms;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_max_natoms() && matoms > _requirements.max_natoms()) {
    ++_too_many_atoms;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_min_nrings() && nrings < _requirements.min_nrings()) {
    ++_too_few_rings;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_max_nrings() && nrings > _requirements.max_nrings()) {
    ++_too_many_rings;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  Molecule m;
  if (! m.build_from_smiles(largest_frag)) {
    cerr << "MoleculeFilterLine:invalid smiles '" << line << "'\n";
    return 0;
  }

  if (m.empty()) {
    cerr << "MoleculeFilterLine:no atoms '" << line << "'\n";
    return 0;
  }

  // Every property here is intended to be about the heavy atom graph, and
  // several of the calculators behave differently, or quietly discard them,
  // when explicit Hydrogens are present. Remove them so the answer does not
  // depend on how the input happened to be written. The smiles scan above has
  // already told us whether there are any, so the common case costs nothing.
  // The output line is left as the caller supplied it - this is a filter, and
  // rewriting structures is reserved for the options that ask for it.
  if (explicit_hydrogens > 0) {
    m.remove_all(1);
    ++_molecules_with_explicit_hydrogens;
  }

  // In LillyMol an isotopic label is usually arbitrary - [1C] is a labelled
  // carbon, not a statement that the atom has no Hydrogens. A strict smiles
  // reading of a bracket atom says otherwise, so let the Hydrogen count be
  // recomputed, which is what transform_to_non_isotopic_form and set_isotopes
  // already do. Doing it here, once, means a property does not depend on
  // whether something else in the pipeline happened to trigger the recompute
  // first - notably ALogP::LogP, which round trips the isotopes and so used to
  // silently change the molecular weight of every isotopic molecule that
  // reached it.
  if (m.ContainsIsotopicAtoms()) {
    const int matoms_now = m.natoms();
    for (int i = 0; i < matoms_now; ++i) {
      if (m.isotope(i) > 0) {
        m.unset_all_implicit_hydrogen_information(i);
      }
    }
    ++_molecules_with_isotopes_rehydrogenated;
  }

  if (_requirements.has_min_chiral() || _requirements.has_max_chiral()) {
    const int chiral_centres = m.chiral_centres();
    if (_requirements.has_min_chiral() && chiral_centres < _requirements.min_chiral()) {
      ++_too_few_chiral;
      return 0;
    }
    if (_requirements.has_max_chiral() && chiral_centres > _requirements.max_chiral()) {
      ++_too_many_chiral;
      return 0;
    }
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    smiles_changed = true;
  }

  if (_chemical_standardisation.active()) {
    if (_chemical_standardisation.process(m)) {
      smiles_changed = true;
      matoms = m.natoms();
    }
  }

  std::vector<double> utility_values;
  double overall_utility = 0.0;

  if (Process(m, matoms, nrings, utility_values, overall_utility)) {
    ++_molecules_passed;

    if (_filter.has_utilities()) {
      const const_IWSubstring output_smiles = smiles_changed ? const_IWSubstring(m.smiles()) : smiles;
      if (! WriteUtilityValues(utility_values, overall_utility, output_smiles, id, output)) {
        MaybeWriteToRejectStream(line);
        return 0;
      }
    } else {
      // If the structure was altered, write the smiles, otherwise the original line
      if (smiles_changed) {
        output << m.smiles() << ' ' << id << '\n';
      } else {
        output << line << '\n';
      }
    }

    output.write_if_buffer_holds_more_than(4092);

    return 1;
  }

  MaybeWriteToRejectStream(line);

  return 0;
}

// The values are computed alongside the filtering, in Process, so that
// descriptors common to both are computed once per molecule.
int
Options::WriteUtilityValues(const std::vector<double>& utility_values,
                            const double overall_utility,
                            const const_IWSubstring& smiles,
                            const const_IWSubstring& id,
                            IWString_and_File_Descriptor& output) {
  if (_write_smiles_with_utility) {
    output << smiles << ' ' << id << ' ' << overall_utility << '\n';
    return 1;
  }

  output << id;
  for (double value : utility_values) {
    output << ' ' << value;
  }
  output << ' ' << overall_utility << '\n';

  return 1;
}

int
Options::MaybeWriteToRejectStream(const const_IWSubstring& line) {
  if (! _reject_stream.is_open()) {
    return 0;
  }

  _reject_stream << line << '\n';

  _reject_stream.write_if_buffer_holds_more_than(4096);

  return 1;
}

void
Options::NoteRejection(molecule_filter_lib::RejectionReason rejection_reason) {
  using molecule_filter_lib::RejectionReason;

  switch (rejection_reason) {
    case RejectionReason::kPass:
      return;
    case RejectionReason::kZeroAtoms:
      return;
    case RejectionReason::kTooFewAtoms:
      ++_too_few_atoms;
      return;
    case RejectionReason::kTooManyAtoms:
      ++_too_many_atoms;
      return;
    case RejectionReason::kTooFewRings:
      ++_too_few_rings;
      return;
    case RejectionReason::kTooManyRings:
      ++_too_many_rings;
      return;
    case RejectionReason::kTooFewHeteroatoms:
      ++_too_few_heteroatoms;
      return;
    case RejectionReason::kTooManyHeteroatoms:
      ++_too_many_heteroatoms;
      return;
    case RejectionReason::kMinHeteroatomFraction:
      ++_min_heteroatom_fraction;
      return;
    case RejectionReason::kMaxHeteroatomFraction:
      ++_max_heteroatom_fraction;
      return;
    case RejectionReason::kNonOrganic:
      ++_non_organic;
      return;
    case RejectionReason::kIsotope:
      ++_isotope;
      return;
    case RejectionReason::kTooManyChiral:
      ++_too_many_chiral;
      return;
    case RejectionReason::kTooFewAromaticRings:
      ++_too_few_aromatic_rings;
      return;
    case RejectionReason::kTooManyAromaticRings:
      ++_too_many_aromatic_rings;
      return;
    case RejectionReason::kTooFewAliphaticRings:
      ++_too_few_aliphatic_rings;
      return;
    case RejectionReason::kTooManyAliphaticRings:
      ++_too_many_aliphatic_rings;
      return;
    case RejectionReason::kTooFewHba:
      ++_too_few_hba;
      return;
    case RejectionReason::kTooManyHba:
      ++_too_many_hba;
      return;
    case RejectionReason::kTooFewHbd:
      ++_too_few_hbd;
      return;
    case RejectionReason::kTooManyHbd:
      ++_too_many_hbd;
      return;
    case RejectionReason::kTooFewSp3Carbon:
      ++_too_few_csp3;
      return;
    case RejectionReason::kTooManySp3Carbon:
      ++_too_many_csp3;
      return;
    case RejectionReason::kSp3CarbonFractionTooLow:
      ++_fcsp3_too_low;
      return;
    case RejectionReason::kSp3CarbonFractionTooHigh:
      ++_fcsp3_too_high;
      return;
    case RejectionReason::kTooFewHalogen:
      ++_too_few_halogens;
      return;
    case RejectionReason::kTooManyHalogen:
      ++_too_many_halogens;
      return;
    case RejectionReason::kRingTooSmall:
      ++_ring_too_small;
      return;
    case RejectionReason::kRingTooLarge:
      ++_ring_too_large;
      return;
    case RejectionReason::kTooFewRotatableBonds:
      ++_too_few_rotbond;
      return;
    case RejectionReason::kTooManyRotatableBonds:
      ++_too_many_rotbond;
      return;
    case RejectionReason::kAromaticDensityTooLow:
      ++_aromdens_too_low;
      return;
    case RejectionReason::kAromaticDensityTooHigh:
      ++_aromdens_too_high;
      return;
    case RejectionReason::kTooShort:
      ++_too_short;
      return;
    case RejectionReason::kTooLong:
      ++_too_long;
      return;
    case RejectionReason::kLowTpsa:
      ++_low_tpsa;
      return;
    case RejectionReason::kHighTpsa:
      ++_high_tpsa;
      return;
    case RejectionReason::kRingSystemTooSmall:
      ++_ring_system_too_small;
      return;
    case RejectionReason::kRingSystemTooLarge:
      ++_ring_system_too_large;
      return;
    case RejectionReason::kTooFewAromaticRingsInSystem:
      ++_too_few_aromatic_rings_in_system;
      return;
    case RejectionReason::kTooManyAromaticRingsInSystem:
      ++_too_many_aromatic_rings_in_system;
      return;
    case RejectionReason::kLowAlogp:
      ++_low_alogp;
      return;
    case RejectionReason::kHighAlogp:
      ++_high_alogp;
      return;
    case RejectionReason::kLowXlogp:
      ++_low_xlogp;
      return;
    case RejectionReason::kHighXlogp:
      ++_high_xlogp;
      return;
    case RejectionReason::kTooFewChiral:
      ++_too_few_chiral;
      return;
    case RejectionReason::kTooFewFragments:
      ++_too_few_fragments;
      return;
    case RejectionReason::kTooManyFragments:
      ++_too_many_fragments;
      return;
    case RejectionReason::kLowQed:
      return;
    case RejectionReason::kHighQed:
      return;
    case RejectionReason::kLowAmw:
      ++_low_amw;
      return;
    case RejectionReason::kHighAmw:
      ++_high_amw;
      return;
    case RejectionReason::kTooFewHbaRdkit:
      ++_too_few_hba_rdkit;
      return;
    case RejectionReason::kTooManyHbaRdkit:
      ++_too_many_hba_rdkit;
      return;
    case RejectionReason::kTooFewHbdRdkit:
      ++_too_few_hbd_rdkit;
      return;
    case RejectionReason::kTooManyHbdRdkit:
      ++_too_many_hbd_rdkit;
      return;
  }
}

// Returns 1 if `m` passes the filters, in which case `utility_values` and
// `overall_utility` have been filled in if any utilities are configured.
// Returns 0 if rejected by a filter, or if a utility could not be computed.
int
Options::Process(Molecule& m,
                 const int matoms,
                 const int nrings,
                 std::vector<double>& utility_values,
                 double& overall_utility) {
  if (m.natoms() != matoms) {
    cerr << "Atom count mismatch " << m.natoms() << " vs " << matoms << '\n';
  }

  molecule_filter_lib::RejectionReason rejection_reason;
  const int rc = _filter.OkAndUtilities(m, matoms, nrings, rejection_reason,
                                        utility_values, overall_utility);
  if (rc == 1) {
    return 1;
  }

  if (rc < 0) {
    cerr << "Options::Process:cannot evaluate utility values\n";
    return 0;
  }

  NoteRejection(rejection_reason);
  return 0;
}


int
Options::SeekIfNeeded(iwstring_data_source& input) const {
  if (_seek_to == 0) {
    return 1;
  }

  if (input.seekg(_seek_to)) {
    return 1;
  }

  cerr << "Options::SeekIfNeeded:cannot seek to " << _seek_to << '\n';
  return 0;
}

int
Options::OkContinue(iwstring_data_source& input) const {
  if (_stop_at == 0) {
    return 1;
  }

  off_t o = input.tellg();
  return o < _stop_at;
}

int
MoleculeFilterLine(Options& options,
                   const const_IWSubstring& line,
                   IWString_and_File_Descriptor& output) {

  options.Process(line, output);

  return 1;
}

int
MoleculeFilter(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  if (! options.SeekIfNeeded(input)) {
    return 0;
  }

  while (input.next_record(buffer)) {
    MoleculeFilterLine(options, buffer, output);

    if (! options.OkContinue(input)) {
      return 1;
    }
  }

  return 1;
}

int
MoleculeFilter(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "MoleculeFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return MoleculeFilter(options, input, output);
}

int
MoleculeFilter(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:F:B:i:u");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

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
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (options.has_utilities() && ! options.write_smiles_with_utility()) {
    options.WriteUtilityHeader(output);
  }

  for (const char * fname : cl) {
    if (! MoleculeFilter(options, fname, output)) {
      cerr << "MoleculeFilter::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace molecule_filter

int
main(int argc, char ** argv) {

  int rc = molecule_filter::MoleculeFilter(argc, argv);

  return rc;
}
