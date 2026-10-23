#include "nvrtspsa.h"

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

namespace nvrtspsa {

void
DisplayOptions(char flag) {
  cerr << " -" << flag << " rvstd        convert to charge separated forms - best compat with RDKit\n";
  cerr << " -" << flag << " zeroP        Phosphorus atoms make zero contribution\n";
  cerr << " -" << flag << " zeroS        Sulphur atoms make zero contribution\n";
  cerr << " -" << flag << " zeroU        all unclassified atoms make zero contribution\n";
  cerr << " -" << flag << " RDKit        as close as possible to the RDKit implementation\n";
}

}  // namespace nvrtspsa

static void
report_novartis_psa_unclassified_atom(const nvrtspsa::NovartisPolarSurfaceArea& options,
                                      Molecule &m, atom_number_t zatom, const Atom &a,
                                      int hcount, int is_aromatic) {
  if (!options.display_psa_unclassified_atom_messages()) {
    return;
  }

  cerr << m.name() << " unclassified '" << a.atomic_symbol();

  formal_charge_t fc = m.formal_charge(zatom);

  if (0 == fc) {
    ;
  } else if (fc > 0) {
    cerr << '+' << fc;
  } else {
    cerr << '-' << (-fc);
  }

  cerr << "' atom in Novartis Polar Surface area, atom " << zatom << '\n';;
  cerr << a.ncon() << " connections, ";
  if (a.formal_charge()) {
    cerr << "charge " << a.formal_charge() << ", ";
  }
  cerr << hcount << " hydrogens, aromatic " << is_aromatic << '\n';;

  return;
}

static double
novartis_polar_surface_area_nitrogen(const nvrtspsa::NovartisPolarSurfaceArea& options,
                                     Molecule &m, atom_number_t zatom, Atom &a,
                                     int is_aromatic) {
  // cerr << m.smiles() << " nitrogen atom " << m.smarts_equivalent_for_atom(zatom) << " atom " << zatom << '\n';

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (const Bond* b : a) {
    //b->debug_print(cerr);
    //cerr << "single " << b->is_single_bond() << " double " << b->is_double_bond() << " other " << m.smarts_equivalent_for_atom(b->other(zatom)) << '\n';

    if (b->is_aromatic()) {
      aromatic_bonds++;
    } else if (b->is_single_bond()) {
      single_bonds++;
    } else if (b->is_double_bond()) {
      double_bonds++;
    } else {
      triple_bonds++;
    }
  }

  const int hcount = a.implicit_hydrogens();
  const int ncon = a.ncon();
  formal_charge_t fc = a.formal_charge();
  // cerr << "Nitrogen atom " << m.smarts_equivalent_for_atom(zatom) << " single " << single_bonds << " double_bonds " << double_bonds << " fc " << fc << " arom " << is_aromatic << '\n';

  if (is_aromatic) {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds) {  // [n](:*):*
      return 12.89;
    }
    if (0 == fc && 0 == hcount && 3 == ncon && 3 == aromatic_bonds) {  // [n](:*)(:*):*
      return 4.41;
    }
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds &&
        1 == single_bonds) {  // [n](-*)(:*):*
      return 4.93;
    }
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds &&
        1 == double_bonds) {  // [n](=*)(:*):*
      return 8.39;
    }
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == aromatic_bonds) {  // [nH](:*):*
      return 15.79;
    }
    if (1 == fc && 0 == hcount && 3 == ncon && 3 == aromatic_bonds) {  // [n+](:*)(:*):*
      return 4.10;
    }
    if (1 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds &&
        1 == single_bonds) {  // [n+](-*):*):*
      return 3.88;
    }
    if (1 == fc && 1 == hcount && 2 == ncon && 2 == aromatic_bonds) {  // [nH+](:*):*
      return 14.14;
    }
  } else {
    if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds &&
        m.in_ring_of_given_size(zatom,
                                3)) {  // case of N1C=C1 maybe handled incorrectly //
                                       // [N]1(-*)-*-*-1     do ring queries first
      return 3.01;
    }
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == single_bonds &&
        m.in_ring_of_given_size(zatom, 3)) {  // [NH]1-*-*-1     do ring queries first
      return 21.94;
    }

    if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds) {  // [N](-*)-*
      return 3.24;
    }
    if (0 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds &&
        1 == double_bonds) {  // [N](-*)=*
      return 12.36;
    }
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == triple_bonds) {  // [N]#*
      return 23.79;
    }
    if (0 == fc && 0 == hcount && 3 == ncon && 1 == single_bonds &&
        2 == double_bonds) {  // [N](-*)(=*)(=*)
      return 11.68;
    }
    if (0 == fc && 0 == hcount && 2 == ncon && 1 == double_bonds &&
        1 == triple_bonds) {  // [N](=*)#*
      return 13.60;
    }
    if (0 == fc && 1 == hcount && 2 == ncon && 2 == single_bonds) {  // [NH](-*)-*
      return 12.03;
    }
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == double_bonds) {  // [NH]=*
      return 23.85;
    }
    if (0 == fc && 2 == hcount && 1 == ncon && 1 == single_bonds) {  // [NH2]-*
      return 26.02;
    }
    if (1 == fc && 0 == hcount && 4 == ncon && 4 == single_bonds) {  // [N+](-*)(-*)(-*)-*
      return 0.0;
    }
    if (1 == fc && 0 == hcount && 3 == ncon && 2 == single_bonds &&
        1 == double_bonds) {  // [N+](-*)(-*)=*
      return 3.01;
    }
    if (1 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds &&
        1 == triple_bonds) {  // [N+](-*)#*
      return 4.36;
    }
    if (1 == fc && 1 == hcount && 3 == ncon && 3 == single_bonds) {  // [NH+](-*)(-*)-*
      return 4.44;
    }
    if (1 == fc && 1 == hcount && 2 == ncon && 1 == single_bonds &&
        1 == double_bonds) {  // [NH+](-*)=*
      return 13.97;
    }
    if (1 == fc && 2 == hcount && 2 == ncon && 2 == single_bonds) {  // [NH2+](-*)-*
      return 16.61;
    }
    if (1 == fc && 2 == hcount && 1 == ncon && 1 == double_bonds) {  // [NH2+]=*
      return 25.59;
    }
    if (1 == fc && 3 == hcount && 1 == single_bonds) {  // [NH3+]-*
      return 27.64;
    }
  }

  // cerr << "unclassified nitrogen\n";
  report_novartis_psa_unclassified_atom(options, m, zatom, a, hcount, is_aromatic);

  if (options.return_zero_for_unclassified_atoms()) {
    return 0.0;
  } else {
    return 30.5 - ncon * 8.2 + hcount * 1.5;
  }
}

static double
novartis_polar_surface_area_oxygen(const nvrtspsa::NovartisPolarSurfaceArea& options,
                                   Molecule &m, atom_number_t zatom, Atom &a,
                                   int is_aromatic) {
  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++) {
    const Bond *b = a[i];

    if (b->is_aromatic()) {
      aromatic_bonds++;
    } else if (b->is_single_bond()) {
      single_bonds++;
    } else if (b->is_double_bond()) {
      double_bonds++;
    } else {
      triple_bonds++;
    }
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (is_aromatic) {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds) {  // [o](:*):*
      return 13.14;
    }
  } else {
    if (0 == fc && 0 == hcount && 2 == ncon &&
        m.in_ring_of_given_size(zatom, 3)) {  // [O]1-*-*-1     do ring queries first
      return 12.53;
    }
    if (0 == fc && 0 == hcount && 2 == ncon &&
        2 == single_bonds) {  // [O](-*)-*     do ring queries first
      return 9.23;
    }
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == double_bonds) {  // [O]=*
      return 17.07;
    }
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == single_bonds) {  // [OH]-*
      return 20.23;
    }
    if (-1 == fc && 0 == hcount && 1 == ncon && 1 == single_bonds) {  // [O-]-*
      return 23.06;
    }
  }

  report_novartis_psa_unclassified_atom(options, m, zatom, a, hcount, is_aromatic);

  if (options.return_zero_for_unclassified_atoms()) {
    return 0.0;
  } else {
    return 10.0;  // just a guess
  }
}

static double
novartis_polar_surface_area_sulphur(const nvrtspsa::NovartisPolarSurfaceArea& options,
                                    Molecule &m, atom_number_t zatom, Atom &a,
                                    int is_aromatic) {
  if (options.zero_for_all_sulphur_atoms()) {
    return 0.0;
  }

  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

//return 0.0;

  for (int i = 0; i < ncon; i++) {
    const Bond *b = a[i];

    if (b->is_aromatic()) {
      aromatic_bonds++;
    } else if (b->is_single_bond()) {
      single_bonds++;
    } else if (b->is_double_bond()) {
      double_bonds++;
    } else {
      triple_bonds++;
    }
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

// #define DEBUG_NOVARTIS_POLAR_SURFACE_AREA_SULPHUR
#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA_SULPHUR
  cerr << "atom " << zatom << " sulphur " << ncon << " connections, " << aromatic_bonds
       << " aromatic, " << single_bonds << " single bonds, " << double_bonds
       << " double bonds " << triple_bonds << " triple bonds\n";
#endif

  if (is_aromatic) {
    if (options.non_zero_contribution_for_SD2()) {
      if (0 == fc && 0 == hcount && 2 == ncon && 2 == aromatic_bonds) {  // [s](:*):*
        return 28.24;
      }
      if (0 == fc && 0 == hcount && 3 == ncon && 2 == aromatic_bonds &&
          1 == double_bonds) {  // [s](=*)(:*):*
        return 21.70;
      }
    }
    return 0.0;  // the paper doesn't seem to use its own data
  } else {
    if (0 == fc && 0 == hcount && 2 == ncon && 2 == single_bonds)  // [S](-*)-*
    {
      if (options.non_zero_contribution_for_SD2()) {
        return 25.30;
      } else {
        return 0.0;
      }
    }
    if (0 == fc && 0 == hcount && 1 == ncon && 1 == double_bonds) {  // [S]=*
      return 32.09;
    }
    if (0 == fc && 0 == hcount && 3 == ncon && 2 == single_bonds &&
        1 == double_bonds) {  // [S](-*)(-*)=*
      return 19.21;
    }
    if (0 == fc && 0 == hcount && 4 == ncon && 2 == single_bonds &&
        2 == double_bonds) {  // [S](-*)(-*)(=*)=*
      return 8.38;            // sometimes they don't use this, but I'll use it
    }
    if (0 == fc && 1 == hcount && 1 == ncon && 1 == single_bonds) {  // [SH]-*
      return 38.80;
    }
  }

  report_novartis_psa_unclassified_atom(options, m, zatom, a, hcount, is_aromatic);

  if (options.return_zero_for_unclassified_atoms()) {
    return 0.0;
  } else {
    return 20.0;  // just a guess
  }
}

static double
novartis_polar_surface_area_phosphorus(const nvrtspsa::NovartisPolarSurfaceArea& options,
                                       Molecule &m, atom_number_t zatom, Atom &a,
                                       int is_aromatic) {
  if (options.zero_for_all_phosphorus_atoms()) {
    return 0.0;
  }

  int ncon = a.ncon();

  int aromatic_bonds = 0;
  int single_bonds = 0;
  int double_bonds = 0;
  int triple_bonds = 0;

  for (int i = 0; i < ncon; i++) {
    const Bond *b = a[i];

    if (b->is_aromatic()) {
      aromatic_bonds++;
    } else if (b->is_single_bond()) {
      single_bonds++;
    } else if (b->is_double_bond()) {
      double_bonds++;
    } else {
      triple_bonds++;
    }
  }

  int hcount = a.implicit_hydrogens();

  formal_charge_t fc = a.formal_charge();

  if (0 == fc && 0 == hcount && 3 == ncon && 3 == single_bonds) {  // [P](-*)(-*)-*
    return 13.59;
  }
  if (0 == fc && 0 == hcount && 2 == ncon && 1 == single_bonds &&
      1 == double_bonds) {  // [P](-*)=*
    return 34.14;
  }
  if (0 == fc && 0 == hcount && 4 == ncon && 3 == single_bonds &&
      1 == double_bonds) {  // [P](-*)(-*)(-*)=*
    return 9.81;
  }
  if (0 == fc && 1 == hcount && 3 == ncon && 2 == single_bonds &&
      1 == double_bonds) {  // [PH](-*)(-*)=*
    return 23.47;
  }

  report_novartis_psa_unclassified_atom(options, m, zatom, a, hcount, is_aromatic);

  if (options.return_zero_for_unclassified_atoms()) {
    return 0.0;
  } else {
    return 10.0;  // just a guess
  }
}

// #define DEBUG_NOVARTIS_POLAR_SURFACE_AREA

/*
  Polar surface areas from
  J. Med Chem 2000, 43, 3714-3717
*/

double
nvrtspsa::NovartisPolarSurfaceArea::PolarSurfaceAreaInner(Molecule &m, const atomic_number_t *z,
                            const Atom **atom, const int *is_aromatic) const {
  double result = 0.0;

  const int matoms = m.natoms();

#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA
  cerr << "novartis_polar_surface_area:processing " << m.smiles() << ' ' << m.name()
       << " with " << matoms << " atoms\n";
#endif

  for (int i = 0; i < matoms; i++) {
    if (6 == z[i]) {
      continue;
    }

    Atom &ai = const_cast<Atom &>(*atom[i]);

    double delta;

    if (7 == z[i]) {
      delta = novartis_polar_surface_area_nitrogen(*this, m, i, ai, is_aromatic[i]);
    } else if (8 == z[i]) {
      delta = novartis_polar_surface_area_oxygen(*this, m, i, ai, is_aromatic[i]);
    } else if (16 == z[i]) {
      delta = novartis_polar_surface_area_sulphur(*this, m, i, ai, is_aromatic[i]);
    } else if (15 == z[i]) {
      delta = novartis_polar_surface_area_phosphorus(*this, m, i, ai, is_aromatic[i]);
    } else {
      delta = 0.0;
    }

    result += delta;

#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA
    if (0.0 != delta) {
      cerr << "atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " value " << delta
           << " result so far " << result << '\n';;
    }
#endif
  }

#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA
  cerr << "Total " << result << '\n';;
#endif

  return result;
}



/*
  The atom classification below assumes the Hydrogen suppressed graph. It reads
  implicit_hydrogens() and ncon(), so an explicit Hydrogen atom both hides a
  Hydrogen from the count and adds a connection. A primary amine written as
  N([H])[H] therefore matches the fully substituted pattern and scores 3.24
  instead of 26.02, an alcohol scores as an ether, and so on. The errors are not
  small; they are the values belonging to a different substitution pattern.

  So any explicit Hydrogens are removed first, on a copy, leaving the caller's
  molecule alone. RemoveAllHydrogenAtoms rather than remove_all(1), since for
  this calculation an isotopic or charged label says nothing about how many
  Hydrogens an atom has.

  Note also that this sets the global aromaticity type to Daylight for the
  duration, and restores it afterwards. A molecule that already had aromaticity
  computed under some other model keeps that perception, since the computation
  here is only "if needed". In practice LillyMol has used Daylight throughout
  for years, so this does not arise.
*/

std::optional<double>
nvrtspsa::NovartisPolarSurfaceArea::HydrogenSuppressedCopy(Molecule &m) const {
  Molecule mcopy(m);

  mcopy.RemoveAllHydrogenAtoms();

  return PolarSurfaceAreaHSuppressed(mcopy);
}

// `m` must have no explicit Hydrogens. Builds the arrays that the computation
// needs and evaluates.
std::optional<double>
nvrtspsa::NovartisPolarSurfaceArea::PolarSurfaceAreaHSuppressed(Molecule &m) const {
  const int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "Cannot compute polar surface area of molecule w/ no atoms\n";
    return std::nullopt;
  }

  int aromsave = global_aromaticity_type();

  set_global_aromaticity_type(Daylight);

  atomic_number_t *z = new atomic_number_t[matoms];
  std::unique_ptr<atomic_number_t[]> free_z(z);
  m.atomic_numbers(z);

  int *is_aromatic = new int[matoms];
  std::unique_ptr<int[]> free_is_aromatic(is_aromatic);

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < matoms; i++) {
    is_aromatic[i] = m.is_aromatic(i);
  }

  const Atom **atom = new const Atom *[matoms];
  std::unique_ptr<const Atom *[]> free_atom(atom);
  m.atoms(atom);

  std::optional<double> rc = PolarSurfaceAreaWithArrays(m, z, atom, is_aromatic);

  set_global_aromaticity_type(aromsave);

  return rc;
}

std::optional<double>
nvrtspsa::NovartisPolarSurfaceArea::PolarSurfaceAreaWithArrays(Molecule &m,
                            const atomic_number_t *z,
                            const Atom **atom, const int *is_aromatic) const {
  // This is a public API entry point, so make sure aromaticity is OK.
  // Because of what we are doing in iwdescr to reverse standardisations.
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();
  if (matoms == 0) {
    return std::nullopt;
  }

  if (! convert_to_charge_separated()) {
    return PolarSurfaceAreaInner(m, z, atom, is_aromatic);
  }

  Molecule mcopy(m);
  static std::unique_ptr<Chemical_Standardisation> chemical_standardisation;
  if (! chemical_standardisation) {
    chemical_standardisation = std::make_unique<Chemical_Standardisation>();
    chemical_standardisation->Activate(CS_REVERSE_NITRO, 0);
    chemical_standardisation->Activate(CS_REVERSE_NV5, 0);
  }

  chemical_standardisation->process(mcopy);

  mcopy.compute_aromaticity_if_needed();

  const Atom **my_atom = new const Atom *[matoms];
  std::unique_ptr<const Atom *[]> free_atom(my_atom);
  mcopy.atoms(my_atom);
  return PolarSurfaceAreaInner(mcopy, z, my_atom, is_aromatic);
}

std::optional<double>
nvrtspsa::NovartisPolarSurfaceArea::PolarSurfaceArea(Molecule &m, const atomic_number_t *z,
                            const Atom **atom, const int *is_aromatic) const {
  // The supplied arrays describe `m` as it stands. Removing Hydrogens renumbers
  // the atoms and invalidates them, so that case has to build its own.
  if (m.natoms(1) > 0) [[unlikely]] {
    return HydrogenSuppressedCopy(m);
  }

  return PolarSurfaceAreaWithArrays(m, z, atom, is_aromatic);
}

std::optional<double>
nvrtspsa::NovartisPolarSurfaceArea::PolarSurfaceArea(Molecule &m) const {
  if (m.natoms(1) > 0) [[unlikely]] {
    return HydrogenSuppressedCopy(m);
  }

  return PolarSurfaceAreaHSuppressed(m);
}

namespace nvrtspsa {

void
NovartisPolarSurfaceArea::SetRDKitCompatibility() {
  _convert_to_charge_separated = 1;
  _zero_for_all_phosphorus_atoms = 1;
  _zero_for_all_sulphur_atoms = 1;
}

int
NovartisPolarSurfaceArea::InitialiseOptions(const Command_Line& cl, char flag) {
  const int nw = cl.option_count(flag);
  for (int i = 0; i < nw; ++i) {
    const_IWSubstring w = cl.string_value(flag, i);
    if (w == "help") {
      DisplayOptions(flag);
      return 1;
    }

    if (w == "quiet") {
      set_display_psa_unclassified_atom_messages(0);
    } else if (w == "verbose") {
      set_display_psa_unclassified_atom_messages(1);
    } else if (w == "zeroU") {
      set_return_zero_for_unclassified_atoms(1);
    } else if (w == "zeroS") {
      set_zero_for_all_sulphur_atoms(1);
    } else if (w == "zeroP") {
      set_zero_for_all_phosphorus_atoms(1);
    } else if (w == "rvstd") {
      set_convert_to_charge_separated(1);
    } else if (w == "RDKit") {
      SetRDKitCompatibility();
    } else {
      cerr << "NovartisPolarSurfaceArea::InitialiseOptions:unrecognised option '" << w << "'\n";
      DisplayOptions(flag);
      return 0;
    }
  }

  return 1;
}

}  // namespace nvrtspsa
