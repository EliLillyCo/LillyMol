#include "nvrtspsa.h"

#include <iostream>
#include <memory>

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"

using std::cerr;
using std::endl;

static int display_psa_unclassified_atom_mesages = 1;

void
set_display_psa_unclassified_atom_mesages(int s) {
  display_psa_unclassified_atom_mesages = s;

  return;
}

static int return_zero_for_unclassified_atoms = 0;

void
set_return_zero_for_unclassified_atoms(int s) {
  return_zero_for_unclassified_atoms = s;

  return;
}

static int non_zero_constribution_for_SD2 = 1;

void
set_non_zero_constribution_for_SD2(int s) {
  non_zero_constribution_for_SD2 = s;
}

static void
report_novartis_psa_unclassified_atom(Molecule &m, atom_number_t zatom, const Atom &a,
                                      int hcount, int is_aromatic) {
  if (!display_psa_unclassified_atom_mesages) {
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

  cerr << "' atom in Novartis Polar Surface area, atom " << zatom << endl;
  cerr << a.ncon() << " connections, ";
  if (a.formal_charge()) {
    cerr << "charge " << a.formal_charge() << ", ";
  }
  cerr << hcount << " hydrogens, aromatic " << is_aromatic << endl;

  return;
}

static double
novartis_polar_surface_area_nitrogen(Molecule &m, atom_number_t zatom, Atom &a,
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

  report_novartis_psa_unclassified_atom(m, zatom, a, hcount, is_aromatic);

  if (return_zero_for_unclassified_atoms) {
    return 0.0;
  } else {
    return 5.0;  // just a guess
  }
}

static double
novartis_polar_surface_area_oxygen(Molecule &m, atom_number_t zatom, Atom &a,
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

  report_novartis_psa_unclassified_atom(m, zatom, a, hcount, is_aromatic);

  if (return_zero_for_unclassified_atoms) {
    return 0.0;
  } else {
    return 10.0;  // just a guess
  }
}

static double
novartis_polar_surface_area_sulphur(Molecule &m, atom_number_t zatom, Atom &a,
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

// #define DEBUG_NOVARTIS_POLAR_SURFACE_AREA_SULPHUR
#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA_SULPHUR
  cerr << "atom " << zatom << " sulphur " << ncon << " connections, " << aromatic_bonds
       << " aromatic, " << single_bonds << " single bonds, " << double_bonds
       << " double bonds " << triple_bonds << " triple bonds\n";
#endif

  if (is_aromatic) {
    if (non_zero_constribution_for_SD2) {
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
      if (non_zero_constribution_for_SD2) {
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

  report_novartis_psa_unclassified_atom(m, zatom, a, hcount, is_aromatic);

  if (return_zero_for_unclassified_atoms) {
    return 0.0;
  } else {
    return 20.0;  // just a guess
  }
}

static double
novartis_polar_surface_area_phosphorus(Molecule &m, atom_number_t zatom, Atom &a,
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

  report_novartis_psa_unclassified_atom(m, zatom, a, hcount, is_aromatic);

  if (return_zero_for_unclassified_atoms) {
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
novartis_polar_surface_area(Molecule &m, const atomic_number_t *z, const Atom **atom,
                            const int *is_aromatic) {
  const int matoms = m.natoms();

  double result = 0.0;

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
      delta = novartis_polar_surface_area_nitrogen(m, i, ai, is_aromatic[i]);
    } else if (8 == z[i]) {
      delta = novartis_polar_surface_area_oxygen(m, i, ai, is_aromatic[i]);
    } else if (16 == z[i]) {
      delta = novartis_polar_surface_area_sulphur(m, i, ai, is_aromatic[i]);
    } else if (15 == z[i]) {
      delta = novartis_polar_surface_area_phosphorus(m, i, ai, is_aromatic[i]);
    } else {
      delta = 0.0;
    }

    result += delta;

#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA
    if (0.0 != delta) {
      cerr << "atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " value " << delta
           << " result so far " << result << endl;
    }
#endif
  }

#ifdef DEBUG_NOVARTIS_POLAR_SURFACE_AREA
  cerr << "Total " << result << endl;
#endif

  return result;
}

double
novartis_polar_surface_area(Molecule &m) {
  int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "Cannot compute polar surface area of molecule w/ no atoms\n";
    return 0.0;
  }

  int aromsave = global_aromaticity_type();

  set_global_aromaticity_type(Daylight);

  atomic_number_t *z = new atomic_number_t[matoms];
  std::unique_ptr<atomic_number_t[]> free_z(z);
  m.atomic_numbers(z);

  const Atom **atom = new const Atom *[matoms];
  std::unique_ptr<const Atom *[]> free_atom(atom);
  m.atoms(atom);

  int *is_aromatic = new int[matoms];
  std::unique_ptr<int[]> free_is_aromatic(is_aromatic);

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < matoms; i++) {
    is_aromatic[i] = m.is_aromatic(i);
  }

  int rc = novartis_polar_surface_area(m, z, atom, is_aromatic);

  set_global_aromaticity_type(aromsave);

  return rc;
}

namespace nvrtspsa {

NovartisPolarSurfaceArea::NovartisPolarSurfaceArea() {
  _display_psa_unclassified_atom_mesages = 1;
  _return_zero_for_unclassified_atoms = 0;
  _non_zero_constribution_for_SD2 = 1;
}

std::optional<double>
NovartisPolarSurfaceArea::PolarSurfaceArea(Molecule& m) {
  if (m.empty()) {
    return std::nullopt;
  }

  return novartis_polar_surface_area(m);
}

}  // namespace nvrtspsa
