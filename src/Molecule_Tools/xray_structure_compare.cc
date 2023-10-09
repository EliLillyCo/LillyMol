/*
  Compare geometric features in 2 3D structures
*/

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/u3b.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int molecules_have_multiple_fragments = 0;

static int compute_initial_rms = 0;

static Molecule_Output_Object stream_for_superimposed_molecules;

static Molecule_Output_Object stream_for_isotopically_labelled_molecules;

static distance_t growth_range = 0.0;

static distance_t bonding_radius = 0.0;

static int produce_text_descriptions = 0;

static int relax_bonding_on_substructure_mismatch = 0;

static int remove_chirality_on_substructure_search_failure = 0;

/*
  We want to flag geometries that fall outside acceptable
  ranges.
  From Yong, we have data on sigma for bond angles and
  bond distances. We then multiply by a sigma factor
  and flag distances outside that tolerance.
*/

static double single_bond_distance_rms = 0.02;
static double non_single_bond_distance_rms = 0.015;

static double bond_angle_rms = 3.0;  // degrees

static double bond_distance_nsigma = 5.0;
static double bond_angle_nsigma = 6.0;

static int remove_explicit_hydrogens = 0;

Molecule_Output_Object directcolorfile_molecule;

static IWString_and_File_Descriptor directcolorfile_colour;

static resizable_array_p<IWString> nsigma_colour_array;

static resizable_array_p<Substructure_Query> delocalise_bonds;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Compares geometric freatures in 3D structures\n";
  cerr << "  -p <fname>    probe molecule\n";
  cerr << "  -f <dist>     atoms within <dist> will go in same fragment\n";
  cerr << "  -b <dist>     atoms within <dist> will be bonded\n";
  cerr << "  -S <fname>    write superimposed molecules to <fname>\n";
  cerr << "  -o <...>      output specification for -S output\n";
  cerr << "  -t <number>   write descriptive info on <n> worst fits\n";
  cerr << "  -x            relax bond type info on substructure search failure\n";
  cerr << "  -c            remove all chirality on substructure search failure\n";
  cerr << "  -r            compute initial RMS - before molecules moved\n";
  cerr << "  -N <fname>    provide atom names, one per line, | at end of each molecule\n";
  cerr << "  -I <fname>    write isotopically labelled molecules to <fname>\n";
  cerr << "  -G ...        specifications for tolerance reporting (rms, sigma), enter '-G help'\n";
  cerr << "  -D <stem>     name for creating csib directcolorfile files\n";
  cerr << "  -X <query>    query defining delocalised bonds (delocalised in input molecule)\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

class Atom_Names : public resizable_array_p<IWString>
{
 private:
  // private functions

  int
  _apply_isotopic_label(Molecule& m, atom_number_t zatom, const IWString& s) const;

 public:
  int
  build(iwstring_data_source&);

  int
  apply_isotopes(Molecule&) const;
};

int
Atom_Names::build(iwstring_data_source& input)
{
  resize(0);  // just in case

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    buffer.truncate_at_first(' ');

    if ('|' == buffer) {
      return _number_elements;
    }

    if (remove_explicit_hydrogens && buffer.starts_with('H')) {
      continue;
    }

    IWString* tmp = new IWString(buffer);
    add(tmp);
  }

  return _number_elements;  // end of file, no vertical bar separator
}

int
Atom_Names::apply_isotopes(Molecule& m) const
{
  if (m.natoms() != _number_elements) {
    cerr << "Atom_Names::apply_isotopes:have " << _number_elements
         << " atom names, but molecule has " << m.natoms() << " atoms, cannot process\n";
    return 0;
  }

  for (int i = 0; i < _number_elements; i++) {
    if (!_apply_isotopic_label(m, i, *(_things[i]))) {
      cerr << "Cannot apply isotope for atom " << i << endl;
      return 0;
    }
  }

  return 1;
}

int
Atom_Names::_apply_isotopic_label(Molecule& m, atom_number_t zatom,
                                  const IWString& s) const
{
  int nchars = s.length();

  if (0 == nchars) {
    return 1;
  }

  if (!isdigit(s.last_item())) {
    return 1;
  }

  int n = 0;
  int tenpower = 1;

  for (int ndx = s.length() - 1; ndx >= 0; ndx--) {
    int tmp = s[ndx] - '0';
    if (tmp < 0 || tmp > 9) {
      break;
    }

    n = tenpower * tmp + n;
    tenpower = tenpower * 10;
  }

  // cerr << "From '" << s << "' get isotope " << n << endl;

  if (0 != n) {
    m.set_isotope(zatom, n);
  }

  return 1;
}

class Geometric_Difference
{
 private:
  double _diff;

  IWString _description;

  IWString _append_text;

  atom_number_t _a1, _a2, _a3, _a4;

  //  private functions

  int
  _append_bond_symbol(const Bond* b);

  void
  _append_symbol_or_name(const Molecule& m, const Atom_Names& atom_names,
                         atom_number_t zatom);

 public:
  Geometric_Difference();

  double
  diff() const
  {
    return _diff;
  }

  void
  set_diff(double d)
  {
    _diff = d;
  }

  void
  set_append_text(const IWString& s)
  {
    _append_text = s;
  }

  void
  reset_append_text()
  {
    _append_text.resize(0);
  }

  int
  set_description(const Molecule&, const Atom_Names& atom_names, atom_number_t,
                  const Bond*, atom_number_t);
  int
  set_description(const Molecule& m, const Atom_Names& atom_names, atom_number_t a1,
                  const Bond* b1, atom_number_t a2, const Bond* b2, atom_number_t a3);
  int
  set_description(const Molecule& m, const Atom_Names& atom_names, atom_number_t a1,
                  const Bond* b1, atom_number_t a2, const Bond* b2, atom_number_t a3,
                  const Bond* b3, atom_number_t a4);

  int
  do_directfile_bond_distance(int nsigma, IWString_and_File_Descriptor& output) const;
  int
  do_directfile_bond_angle(int nsigma, IWString_and_File_Descriptor& output) const;

  int
  report(const char* s, IWString_and_File_Descriptor& output) const;
};

Geometric_Difference::Geometric_Difference()
{
  _diff = 0.0;

  _a1 = _a2 = _a3 = _a4 = INVALID_ATOM_NUMBER;

  return;
}

int
Geometric_Difference::_append_bond_symbol(const Bond* b)
{
  if (b->is_aromatic()) {
    _description << ':';
  } else if (b->is_single_bond()) {
    _description << '-';
  } else if (b->is_double_bond()) {
    _description << '=';
  } else if (b->is_triple_bond()) {
    _description << '#';
  } else {
    _description << '?';
  }

  return 1;
}

void
Geometric_Difference::_append_symbol_or_name(const Molecule& m,
                                             const Atom_Names& atom_names,
                                             atom_number_t zatom)
{
  if (0 == atom_names.number_elements()) {
    _description << m.atomic_symbol(zatom) << (zatom + 1);
  } else {
    _description << *(atom_names[zatom]);
  }

  return;
}

// clang-format off
int
Geometric_Difference::set_description (const Molecule & m,
                                       const Atom_Names & atom_names,
                                       atom_number_t a1,
                                       const Bond * b,
                                       atom_number_t a2)
// clang-format on
{
  _description.resize_keep_storage(0);

  _append_symbol_or_name(m, atom_names, a1);

  _append_bond_symbol(b);

  _append_symbol_or_name(m, atom_names, a2);

  _a1 = a1;
  _a2 = a2;

  return 1;
}

// clang-format off
int
Geometric_Difference::set_description(const Molecule & m,
                                      const Atom_Names & atom_names,
                                      atom_number_t a1,
                                      const Bond * b1,
                                      atom_number_t a2,
                                      const Bond * b2,
                                      atom_number_t a3)
// clang-format on
{
  set_description(m, atom_names, a1, b1, a2);

  _append_bond_symbol(b2);

  _append_symbol_or_name(m, atom_names, a3);

  _a3 = a3;

  return 1;
}

// clang-format off
int
Geometric_Difference::set_description(const Molecule & m,
                                      const Atom_Names & atom_names,
                                      atom_number_t a1,
                                      const Bond * b1,
                                      atom_number_t a2,
                                      const Bond * b2,
                                      atom_number_t a3,
                                      const Bond * b3,
                                      atom_number_t a4)
// clang-format on
{
  set_description(m, atom_names, a1, b1, a2, b2, a3);

  _append_bond_symbol(b3);

  _append_symbol_or_name(m, atom_names, a4);

  _a4 = a4;

  return 1;
}

int
Geometric_Difference::report(const char* s, IWString_and_File_Descriptor& output) const
{
  output << s << _description << ' ' << static_cast<float>(_diff);

  if (_append_text.length()) {
    output << _append_text;
  }

  output << '\n';

  return 1;
}

class Geometric_Difference_Comparitor
{
 private:
 public:
  int
  operator()(const Geometric_Difference&, const Geometric_Difference&) const;
};

int
Geometric_Difference_Comparitor::operator()(const Geometric_Difference& gd1,
                                            const Geometric_Difference& gd2) const
{
  if (gd1.diff() < gd2.diff()) {
    return 1;
  }

  if (gd1.diff() > gd2.diff()) {
    return -1;
  }

  return 0;
}

static int
superimpose_molecules(const Molecule& m1, const int* xref, Molecule& probe)
{
  int matoms = m1.natoms();

  double* c1 = new double[matoms * 3];
  std::unique_ptr<double[]> free_c1(c1);
  double* c2 = new double[matoms * 3];
  std::unique_ptr<double[]> free_c2(c2);

  double* weight = new double[matoms];
  std::unique_ptr<double[]> free_weight(weight);

  for (int i = 0; i < matoms; i++) {
    weight[i] = 1.0;

    int ix = xref[i];

    const Atom* a = m1.atomi(ix);

    c1[3 * i] = a->x();
    c1[3 * i + 1] = a->y();
    c1[3 * i + 2] = a->z();

    a = probe.atomi(i);

    c2[3 * i] = a->x();
    c2[3 * i + 1] = a->y();
    c2[3 * i + 2] = a->z();
  }

  int mode = 1;
  double u[9];
  double rms;
  double t[3];
  int ier = 0;

  u3b_(weight, c1, c2, &matoms, &mode, &rms, u, t, &ier);

  if (verbose > 1) {
    cerr << "Between " << m1.name() << " and " << probe.name() << " rms " << rms << '\n';
  }

  double rotmat11 = u[0];
  double rotmat12 = u[1];
  double rotmat13 = u[2];
  double rotmat21 = u[3];
  double rotmat22 = u[4];
  double rotmat23 = u[5];
  double rotmat31 = u[6];
  double rotmat32 = u[7];
  double rotmat33 = u[8];

  for (int i = 0; i < matoms; i++) {
#ifdef DEBUG_PROCESS_3D_REPLACE
    cerr << "Atom " << i << " '" << m.smarts_equivalent_for_atom(i) << "' is moving\n";
#endif

    const Atom* a = probe.atomi(i);

    double x0 = a->x() - t[0];
    double y0 = a->y() - t[1];
    double z0 = a->z() - t[2];

    double xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    double yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    double zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    probe.setxyz(i, xx, yy, zz);
  }

  return 1;
}

static void
write_both_identifiers(const IWString& m1name, const IWString& m2name,
                       IWString_and_File_Descriptor& output)
{
  output << m1name << ' ' << m2name << ' ';

  return;
}

/*
  Taking the difference between angles is difficult because things might wrap
  near 180 degrees
*/

static angle_t
angle_difference(angle_t a1, angle_t a2)
{
#ifdef DEBUG_ANGLE_DIFFERENCE
  if (fabs(a1 - a2) < 180.0) {
    cerr << "Angles " << a1 << " and " << a2 << " yield " << fabs(a1 - a2) << endl;
  }
#endif

  if (fabs(a1 - a2) < 180.0) {
    return fabs(a1 - a2);
  }

  // What happens if one is -178 and the other is +178. Result should be 4

#ifdef DEBUG_ANGLE_DIFFERENCE
  cerr << "*Angles " << a1 << " and " << a2 << " yield " << (360 - fabs(a1 - a2)) << endl;
#endif

  return 360 - fabs(a1 - a2);
}

static int
do_torsions(const Molecule& m1, const Atom_Names& atom_names, const Molecule& m2,
            const int* xref, Geometric_Difference* gd,
            IWString_and_File_Descriptor& output)
{
  int matoms = m1.natoms();

  Accumulator<angle_t> acc;

  int ndx = 0;

  for (int i = 0; i < matoms; i++) {
    int ix = xref[i];

    const Atom* ai = m1.atomi(i);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++) {
      const Bond* b1 = ai->item(j);

      atom_number_t k = b1->other(i);

      int kx = xref[k];

      const Atom* ak = m1.atomi(k);

      int kcon = ak->ncon();

      for (int l = 0; l < kcon; l++) {
        const Bond* b2 = ak->item(l);

        atom_number_t m = b2->other(k);

        if (m == i) {
          continue;
        }

        int mx = xref[m];

        const Atom* am = m1.atomi(m);

        int mcon = am->ncon();

        for (int n = 0; n < mcon; n++) {
          const Bond* b3 = am->item(n);

          atom_number_t o = b3->other(m);

          if (o == k || o <= i) {  // guard against 3 membered ring
            continue;
          }

          int ox = xref[o];

          //        cerr << "atoms " << i << ' ' << k << ' ' << m << ' ' << o << endl;

          angle_t a1 = m1.dihedral_angle(i, k, m, o) * RAD2DEG;
          angle_t a2 = m2.dihedral_angle(ix, kx, mx, ox) * RAD2DEG;

          double diff = angle_difference(a1, a2);

          acc.extra(diff);

          if (produce_text_descriptions) {
            gd[ndx].set_diff(diff);
            gd[ndx].set_description(m1, atom_names, i, b1, k, b2, m, b3, o);
            gd[ndx].reset_append_text();
            ndx++;
          }
        }
      }
    }
  }

  if (0 == acc.n()) {
    cerr << "No matching dihedrals found\n";
    return 0;
  }

  write_both_identifiers(m1.name(), m2.name(), output);

  output << "Dihedral " << acc.n() << " between " << acc.minval() << " and "
         << acc.maxval() << " ave " << static_cast<float>(acc.average()) << '\n';

  if (produce_text_descriptions) {
    Geometric_Difference_Comparitor gdc;
    iwqsort(gd, ndx, gdc);

    int istop = produce_text_descriptions;
    if (istop > ndx) {
      istop = ndx;
    }

    for (int i = 0; i < istop; i++) {
      gd[i].report(" DIHEDRAL ", output);
    }
  }

  return 1;
}

static int
append_bond_angle_out_of_bounds_message(Geometric_Difference& gd, const angle_t ideal,
                                        const angle_t actual, const double diff)
{
  // cerr << "BARMS " << bond_angle_rms << " DIFF " << diff << " bans " <<
  // bond_angle_nsigma << endl;

  if (bond_angle_rms < 0.0) {
    return 0;
  }

  if (diff < (bond_angle_rms * bond_angle_nsigma)) {
    return 0;
  }

  int nsigma = static_cast<int>(diff / bond_angle_rms);

  IWString tmp;

  tmp << " !!! " << nsigma << " sigma ideal " << ideal << " actual " << actual;

  gd.set_append_text(tmp);

  if (directcolorfile_colour.is_open()) {
    gd.do_directfile_bond_angle(nsigma, directcolorfile_colour);
  }

  // cerr << "Nsigma " << nsigma << endl;

  return 1;
}

static int
do_bond_angles(const Molecule& m1, const Atom_Names& atom_names,
               const Molecule& m2,  // probe
               const int* xref, Geometric_Difference* gd,
               IWString_and_File_Descriptor& output)
{
  int matoms = m1.natoms();

  Accumulator<angle_t> acc;

  int ndx = 0;

  for (int i = 0; i < matoms; i++) {
    int ix = xref[i];

    const Atom* ai = m1.atomi(i);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++) {
      const Bond* b1 = ai->item(j);

      atom_number_t k = b1->other(i);

      int kx = xref[k];

      for (int l = j + 1; l < acon; l++) {
        const Bond* b2 = ai->item(l);

        atom_number_t m = b2->other(i);

        int mx = xref[m];

        angle_t a1 = m1.bond_angle(k, i, m) * RAD2DEG;
        angle_t a2 = m2.bond_angle(kx, ix, mx) * RAD2DEG;

        double diff = angle_difference(a1, a2);
        acc.extra(diff);

        if (produce_text_descriptions) {
          gd[ndx].set_diff(diff);
          gd[ndx].set_description(m1, atom_names, k, b1, i, b2, m);
          gd[ndx].reset_append_text();
          if (bond_angle_nsigma > 0.0) {
            append_bond_angle_out_of_bounds_message(gd[ndx], a2, a1, diff);
          }
          ndx++;
        }
      }
    }
  }

  if (0 == acc.n()) {
    cerr << "No matching bond angles found\n";
    return 0;
  }

  write_both_identifiers(m1.name(), m2.name(), output);

  output << "Bond Angles " << acc.n() << " between " << acc.minval() << " and "
         << acc.maxval() << " ave " << static_cast<float>(acc.average()) << '\n';

  if (produce_text_descriptions) {
    Geometric_Difference_Comparitor gdc;
    iwqsort(gd, ndx, gdc);

    int istop = produce_text_descriptions;
    if (istop > ndx) {
      istop = ndx;
    }

    for (int i = 0; i < istop; i++) {
      gd[i].report(" ANGLE ", output);
    }
  }

  return 1;
}

int
Geometric_Difference::do_directfile_bond_distance(
    int nsigma, IWString_and_File_Descriptor& output) const
{
  directcolorfile_colour << _a1 << ' ' << _a2 << ' ';

  if (nsigma > nsigma_colour_array.number_elements()) {
    directcolorfile_colour << *(nsigma_colour_array.last_item());
  } else {
    directcolorfile_colour << *(nsigma_colour_array[nsigma]);
  }

  directcolorfile_colour << '\n';

  return 1;
}

int
Geometric_Difference::do_directfile_bond_angle(int nsigma,
                                               IWString_and_File_Descriptor& output) const
{
  const IWString* c;
  if (nsigma > nsigma_colour_array.number_elements()) {
    c = nsigma_colour_array.last_item();
  } else {
    c = nsigma_colour_array[nsigma];
  }

  directcolorfile_colour << _a1 << ' ' << _a2 << ' ' << (*c) << '\n';
  directcolorfile_colour << _a2 << ' ' << _a3 << ' ' << (*c) << '\n';

  return 1;
}

static int
append_bond_distance_out_of_bounds_message(Geometric_Difference& gd, double ideal,
                                           double actual, double diff, const Bond* b)
{
  // cerr << "Bond distanance discrepancy " << diff << " arom ? " << b->is_aromatic() << "
  // single " << b->is_single_bond() << endl;

  int nsigma = 0;

  if (b->is_aromatic()) {
    if (non_single_bond_distance_rms < 0.0) {
      return 0;
    }

    if (diff < (bond_distance_nsigma * non_single_bond_distance_rms)) {
      return 0;
    }

    nsigma = static_cast<int>(diff / non_single_bond_distance_rms);
  } else if (b->is_single_bond()) {
    if (single_bond_distance_rms < 0.0) {
      return 0;
    }

    if (diff < (bond_distance_nsigma * single_bond_distance_rms)) {
      return 0;
    }

    nsigma = static_cast<int>(diff / single_bond_distance_rms);
  } else {
    if (non_single_bond_distance_rms < 0.0) {
      return 0;
    }

    if (diff < (bond_distance_nsigma * non_single_bond_distance_rms)) {
      return 0;
    }

    nsigma = static_cast<int>(diff / non_single_bond_distance_rms);
  }

  // cerr << "Bond distance out of bopunds, nsigma " << nsigma << endl;

  IWString tmp;

  tmp << " !!! " << nsigma << " sigma ideal " << ideal << " actual " << actual;

  gd.set_append_text(tmp);

  if (directcolorfile_colour.is_open()) {
    gd.do_directfile_bond_distance(nsigma, directcolorfile_colour);
  }

  return 1;
}

static int
do_bond_distances(const Molecule& m1, const Atom_Names& atom_names,
                  const Molecule& m2,  // probe
                  const int* xref, Geometric_Difference* gd,
                  IWString_and_File_Descriptor& output)
{
  int matoms = m1.natoms();

  Accumulator<distance_t> acc;

  int ndx = 0;

  for (int i = 0; i < matoms; i++) {
    int ix = xref[i];

    const Atom* ai = m1.atomi(i);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++) {
      const Bond* b = ai->item(j);

      atom_number_t k = b->other(i);

      int kx = xref[k];

      if (kx < ix) {  // also covers the case of kx negative
        continue;
      }

      distance_t d1 = m1.distance_between_atoms(i, k);
      distance_t d2 = m2.distance_between_atoms(ix, kx);

      double d = fabs(d1 - d2);
      acc.extra(d);

      if (!produce_text_descriptions) {
        continue;
      }

      gd[ndx].set_diff(d);
      gd[ndx].set_description(m1, atom_names, i, b, k);
      gd[ndx].reset_append_text();

      if (bond_distance_nsigma > 0.0) {
        append_bond_distance_out_of_bounds_message(gd[ndx], d2, d1, d, b);
      }

      ndx++;
    }
  }

  if (0 == acc.n()) {
    cerr << "No bond distances encountered\n";
    return 0;
  }

  write_both_identifiers(m1.name(), m2.name(), output);

  output << "Bond Distances " << acc.n() << " between " << acc.minval() << " to "
         << acc.maxval() << " ave " << static_cast<float>(acc.average()) << '\n';

  if (produce_text_descriptions && ndx > 1) {
    Geometric_Difference_Comparitor gdc;
    iwqsort(gd, ndx, gdc);

    int istop = produce_text_descriptions;
    if (istop > ndx) {
      istop = ndx;
    }

    for (int i = 0; i < istop; i++) {
      gd[i].report(" BOND ", output);
    }
  }

  return 1;
}

struct Atom_and_Distance {
  atom_number_t _atom_number;
  distance_t _distance;
};

class Atom_and_Distance_Comparitor
{
 private:
 public:
  int
  operator()(const Atom_and_Distance&, const Atom_and_Distance&) const;
};

int
Atom_and_Distance_Comparitor::operator()(const Atom_and_Distance& ad1,
                                         const Atom_and_Distance& ad2) const
{
  if (ad1._distance < ad2._distance) {
    return -1;
  } else if (ad1._distance > ad2._distance) {
    return 1;
  }

  return 0;
}

template void
iwqsort(Atom_and_Distance*, int, Atom_and_Distance_Comparitor&);
template void
iwqsort<Atom_and_Distance, Atom_and_Distance_Comparitor>(Atom_and_Distance*, int,
                                                         Atom_and_Distance_Comparitor&,
                                                         void*);
template void
swap_elements<Atom_and_Distance>(Atom_and_Distance&, Atom_and_Distance&, void*);
template void
move_in_from_left<Atom_and_Distance, Atom_and_Distance_Comparitor>(
    Atom_and_Distance*, int&, int&, int, Atom_and_Distance_Comparitor&, void*);
template void
move_in_from_right<Atom_and_Distance, Atom_and_Distance_Comparitor>(
    Atom_and_Distance*, int&, int&, Atom_and_Distance_Comparitor&);
template void
compare_two_items<Atom_and_Distance, Atom_and_Distance_Comparitor>(
    Atom_and_Distance*, Atom_and_Distance_Comparitor&, void*);

template void
iwqsort(Geometric_Difference*, int, Geometric_Difference_Comparitor&);
template void
iwqsort<Geometric_Difference, Geometric_Difference_Comparitor>(
    Geometric_Difference*, int, Geometric_Difference_Comparitor&, void*);
template void
swap_elements<Geometric_Difference>(Geometric_Difference&, Geometric_Difference&, void*);
template void
move_in_from_left<Geometric_Difference, Geometric_Difference_Comparitor>(
    Geometric_Difference*, int&, int&, int, Geometric_Difference_Comparitor&, void*);
template void
move_in_from_right<Geometric_Difference, Geometric_Difference_Comparitor>(
    Geometric_Difference*, int&, int&, Geometric_Difference_Comparitor&);
template void
compare_two_items<Geometric_Difference, Geometric_Difference_Comparitor>(
    Geometric_Difference*, Geometric_Difference_Comparitor&, void*);

/*
  The maximum connectivity associated with each atomic number
*/

static extending_resizable_array<int> max_con(1);

static void
initialise_max_con()
{
  max_con.extend(HIGHEST_ATOMIC_NUMBER + 1);

  max_con[6] = 4;
  max_con[7] = 4;
  max_con[8] = 2;
  max_con[15] = 4;
  max_con[16] = 4;
  max_con[17] = 4;

  return;
}

static int
add_bonds(Molecule& m)
{
  int matoms = m.natoms();

  distance_t* d = new distance_t[matoms * matoms];
  std::unique_ptr<distance_t[]> free_d(d);

  Atom_and_Distance* ad = new Atom_and_Distance[matoms];
  std::unique_ptr<Atom_and_Distance[]> free_ad(ad);

  for (int i = 0; i < matoms; i++) {
    for (int j = i + 1; j < matoms; j++) {
      distance_t dij = m.distance_between_atoms(i, j);

      d[i * matoms + j] = d[j * matoms + i] = dij;
    }
  }

  Atom_and_Distance_Comparitor adc;

  for (int i = 0; i < matoms; i++) {
    int ndx = 0;
    for (int j = 0; j < matoms; j++) {
      if (i == j) {
        continue;
      }

      distance_t dij = d[i * matoms + j];

      if (dij > bonding_radius) {
        continue;
      }

      ad[ndx]._atom_number = j;
      ad[ndx]._distance = dij;
      ndx++;
    }

    if (ndx > 1) {
      iwqsort(ad, ndx, adc);
    }

    atomic_number_t zi = m.atomic_number(i);

    int ncon = max_con[zi];

    if (ncon > ndx) {
      ncon = ndx;
    }

    for (int j = 0; j < ncon; j++) {
      const Atom_and_Distance& adi = ad[j];

      if (i < adi._atom_number) {
        m.add_bond(i, adi._atom_number, SINGLE_BOND);
      }
    }
  }

  // Get rid of any bad valences by deleting the longest bond

  for (int i = 0; i < matoms; i++) {
    const Atom* ai = m.atomi(i);

    atomic_number_t z = ai->atomic_number();

    int acon = ai->ncon();

    if (acon <= max_con[z]) {
      continue;
    }

    for (int j = 0; j < acon; j++) {
      atom_number_t k = ai->other(i, j);

      ad[j]._atom_number = k;
      ad[j]._distance = d[i * matoms + k];
    }

    iwqsort(ad, acon, adc);

    while (ai->ncon() > max_con[z]) {
      atom_number_t j = ad[ai->ncon() - 1]._atom_number;

      m.remove_bond_between_atoms(i, j);
    }
  }

  if (verbose > 1) {
    cerr << "After filling in bonding '" << m.smiles() << "' " << m.name() << "\n";
  }

  return 1;
}

static void
fill_cross_reference_array(const Set_of_Atoms& e, const Query_Atoms_Matched& qam,
                           int* xref)
{
  int n = qam.number_elements();
  assert(n == e.number_elements());

  for (int i = 0; i < n; i++) {
    const Substructure_Atom* a = qam[i];

    atom_number_t j = a->initial_atom_number();

    atom_number_t k = e[i];

    xref[j] = k;

    //  cerr << "Atom " << k << " in molecule matched with probe atom " << j << "\n";
  }

  return;
}

static int
do_compute_initial_rms(const Molecule& m1, const Molecule& m2, const int* xref,
                       IWString_and_File_Descriptor& output)
{
  int matoms = m1.natoms();

  Accumulator<double> d2acc;

  for (int i = 0; i < matoms; i++) {
    const Atom* a1 = m1.atomi(i);

    const Atom* a2 = m2.atomi(xref[i]);

    coord_t d = a1->distance(*a2);

    d2acc.extra(d * d);
  }

  double rms = sqrt(d2acc.average());

  output << m1.name() << ' ' << m2.name() << " RMS " << static_cast<float>(rms) << '\n';

  return 1;
}

static int
xray_structure_compare(const Molecule& m1, const Atom_Names& atom_names,
                       const Molecule& m2,  // probe
                       int* xref, Geometric_Difference* gd,
                       IWString_and_File_Descriptor& output)
{
  if (0 == m1.nedges() || 0 == m2.nedges()) {
    cerr << "One or both of the molecules have zero bonds, cannot continue\n";
    return 0;
  }

#ifdef ECHO_CROSS_REFERENCE_ARRAY
  for (int i = 0; i < m1.natoms(); i++) {
    cerr << "xref " << i << " value " << xref[i] << endl;
  }
#endif

  if (m1.natoms() != m2.natoms()) {
    cerr << "Warning, molecule contains " << m1.natoms() << " atoms, ideal contains "
         << m2.natoms() << "\n";
  }

  if (compute_initial_rms) {
    do_compute_initial_rms(m1, m2, xref, output);
  }

  do_bond_distances(m1, atom_names, m2, xref, gd, output);
  do_bond_angles(m1, atom_names, m2, xref, gd, output);
  do_torsions(m1, atom_names, m2, xref, gd, output);

  return 1;
}

static void
fill_inter_atom_distance_array(const Molecule& m, distance_t* d)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    for (int j = i + 1; j < matoms; j++) {
      d[i * matoms + j] = d[j * matoms + i] = m.distance_between_atoms(i, j);
    }
  }

  return;
}

static double
compute_geometric_match(const Molecule& m1,
                        const Molecule& m2,  // probe
                        const Set_of_Atoms& e, const Query_Atoms_Matched& qam, int* xref,
                        distance_t* d)
{
  const int matoms = m1.natoms();

  fill_inter_atom_distance_array(m1, d);

  set_vector(xref, matoms, -1);

  fill_cross_reference_array(e, qam, xref);

  double rc = 0.0;

  for (int i = 0; i < matoms; i++) {
    int ix = xref[i];

    for (int j = i + 1; j < matoms; j++) {
      int jx = xref[j];

      distance_t dij = m2.distance_between_atoms(ix, jx);

// #define DEBUG_COMPUTE_GEOMETRIC_MATCH
#ifdef DEBUG_COMPUTE_GEOMETRIC_MATCH
      cerr << " i = " << i << " ix = " << ix << " j = " << j << " jx = " << jx << " dist "
           << dij << endl;
      cerr << "  other dist " << d[i * matoms + j] << " fabs "
           << fabs(dij - d[i * matoms + j]) << endl;
#endif

      rc += fabs(dij - d[i * matoms + j]);
    }
  }

#ifdef DEBUG_COMPUTE_GEOMETRIC_MATCH
  cerr << "Returning " << rc << endl;
#endif

  return rc;
}

static void
transfer_bonding_information(Molecule& m, const int* xref, const Molecule& probe)
{
  int nedges = m.nedges();

  for (int i = 0; i < nedges; i++) {
    const Bond* b1 = m.bondi(i);

    atom_number_t a1 = b1->a1();
    atom_number_t a2 = b1->a2();

    int x1 = xref[a1];
    int x2 = xref[a2];

    if (!probe.are_bonded(x1, x2)) {
      cerr << "Yipes, atoms " << a1 << " xref " << x1 << " and " << a2 << " xref " << x2
           << " not bonded in probe, see Ian\n";
      continue;
    }

    const Bond* b2 = probe.bond_between_atoms(x1, x2);

    if (b1->is_single_bond()) {
      if (!b2->is_single_bond()) {
        m.set_bond_type_between_atoms(a1, a2, b2->btype());
      }
    } else if (b1->is_double_bond()) {
      if (!b2->is_double_bond()) {
        m.set_bond_type_between_atoms(a1, a2, b2->btype());
      }
    } else if (b1->is_triple_bond()) {
      if (!b2->is_triple_bond()) {
        m.set_bond_type_between_atoms(a1, a2, b2->btype());
      }
    } else {
    }
  }

  return;
}

static int
do_delocalise_bonds(Molecule& m, Substructure_Query& q, const atom_number_t a1,
                    const atom_number_t a2)
{
  Substructure_Bond* bond_between_atoms =
      const_cast<Substructure_Bond*>(q.bond_between_atoms(a1, a2));
  if (nullptr == bond_between_atoms) {  // should not happen
    return 0;
  }

  if (verbose > 1) {
    cerr << "Setting match any bond " << m.smarts_equivalent_for_atom(a1) << " and "
         << m.smarts_equivalent_for_atom(a2) << endl;
  }

  bond_between_atoms->set_match_any();

  return 1;
}

static int
do_delocalise_bonds(Molecule& m, Substructure_Query& q, const Set_of_Atoms& e)
{
  const int n = e.number_elements();

  int rc = 0;
  for (int i = 0; i < n; ++i) {
    const atom_number_t a1 = e[i];
    for (int j = i + 1; j < n; ++j) {
      const atom_number_t a2 = e[j];
      if (!m.are_bonded(a1, a2)) {
        continue;
      }

      rc += do_delocalise_bonds(m, q, a1, a2);
    }
  }

  return 1;
}

static int
do_delocalise_bonds(Molecule& m, Substructure_Query& q,
                    Substructure_Query& delocalise_bonds)
{
  Substructure_Results sresults;
  const int nhits = delocalise_bonds.substructure_search(&m, sresults);
  if (0 == nhits) {
    return 0;
  }

  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);

    do_delocalise_bonds(m, q, *e);
  }

  return nhits;
}

static int
do_delocalise_bonds(Molecule& m, Substructure_Query& q,
                    resizable_array_p<Substructure_Query>& delocalise_bonds)
{
  int rc = 0;

  for (int i = 0; i < delocalise_bonds.number_elements(); ++i) {
    if (do_delocalise_bonds(m, q, *delocalise_bonds[i])) {
      rc++;
    }
  }

  IWString foo("foo.msi");
  q.write_msi(foo);

  return rc;
}

static int
do_substructure_search(Molecule& probe, Molecule& m, Substructure_Results& sresults)
{
  const int matoms = m.natoms();

  if (verbose > 1) {
    cerr << "matoms " << matoms << " probe " << probe.natoms() << endl;
  }

  if (matoms > probe.natoms()) {
    cerr << "Molecule has more atoms " << matoms << " than probe " << probe.natoms()
         << " cannot continue\n";
    return 0;
  }

  Molecule_to_Match target(&probe);

  Molecule_to_Query_Specifications mqs;

  mqs.set_make_embedding(1);
  mqs.set_all_ring_bonds_become_undefined(1);
  if (delocalise_bonds.number_elements()) {
    mqs.set_ignore_molecular_hydrogen_information(1);
  }

  Substructure_Query q;
  if (!q.create_from_molecule(m, mqs)) {
    cerr << "Cannot create query object from '" << m.name() << "'\n";
    return 0;
  }

  q.set_respect_initial_atom_numbering(0);

  if (delocalise_bonds.number_elements()) {
    do_delocalise_bonds(m, q, delocalise_bonds);
  }

  int nhits = q.substructure_search(target, sresults);

  if (nhits > 0) {
    return nhits;
  }

  if (remove_chirality_on_substructure_search_failure &&
      (m.chiral_centres() || probe.chiral_centres())) {
    cerr << "Substructure search failure, removing chirality\n";
    m.remove_all_chiral_centres();
    probe.remove_all_chiral_centres();
    return do_substructure_search(probe, m, sresults);
  }

  if (!relax_bonding_on_substructure_mismatch) {
    cerr << "OOps, no hits looking for '" << m.name() << "' in '" << probe.name()
         << "'\n";
    cerr << "Only matched " << sresults.max_query_atoms_matched_in_search() << " of "
         << matoms << " query atoms\n";
    return 0;
  }

  cerr << "No substructure match with '" << m.name() << "' "
       << sresults.max_query_atoms_matched_in_search()
       << " atoms matched, relaxing bonding information\n";

  mqs.set_just_atomic_number_and_connectivity(1);

  Substructure_Query q2;
  if (!q2.create_from_molecule(m, mqs)) {  // how could this happen?
    return 0;
  }

  q2.set_respect_initial_atom_numbering(0);

  if (q2.substructure_search(target, sresults)) {
    return 1;
  }

  // Start stripping off terminal groups...
  // No, cannot do this, the atom numbering needs to change, too hard for now,
  // come back later (if ever)

  return 0;

  for (int i = 0; i < matoms; i++) {
    if (1 != m.ncon(i)) {
      continue;
    }

    Molecule mcopy(m);
    mcopy.remove_atom(i);
    if (do_substructure_search(probe, mcopy, sresults)) {
      return 1;
    }
  }

  return 0;
}

static int
xray_structure_compare2(Molecule& m, const Atom_Names& atom_names,
                        const resizable_array_p<Molecule>& probe,
                        IWString_and_File_Descriptor& output)
{
  const int matoms = m.natoms();

  if (0 == atom_names.number_elements()) {
    ;
  } else if (atom_names.number_elements() != matoms) {
    cerr << "Size mismatch on atom names, molecule has " << matoms << " atoms, but "
         << atom_names.number_elements() << " atom names available\n";
    return 0;
  }

  if (3 != m.highest_coordinate_dimensionality()) {
    cerr << "Sorry, only works with 3D structures '" << m.name() << "'\n";
    return 0;
  }

  int atoms_in_probe = probe[0]->natoms();

  for (int i = 1; i < probe.number_elements(); i++) {
    if (probe[i]->natoms() > atoms_in_probe) {
      atoms_in_probe = probe[i]->natoms();
    }
  }

  int* xref = new int[matoms];
  std::unique_ptr<int[]> free_xref(xref);
  distance_t* d = new distance_t[matoms * matoms];
  std::unique_ptr<distance_t[]> free_d(d);

  double zbest = std::numeric_limits<double>::max();
  int best_embedding = -1;

  int np = probe.number_elements();

  for (int i = 0; i < np; i++) {
    Molecule* pi = probe[i];

    Substructure_Results sresults;

    int nhits = do_substructure_search(*pi, m, sresults);

    if (0 == nhits) {
      cerr << "Yipes, substructure search failure, cannot continue\n";
      return 0;
    }

    if (verbose > 1) {
      cerr << "Got " << nhits << " probe substructure matches\n";
    }

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);
      const Query_Atoms_Matched* qam = sresults.query_atoms_matching(j);

#ifdef DEBUG_COMPUTE_GEOMETRIC_MATCH
      cerr << "Testing " << (*e) << endl;
#endif

      double fit = compute_geometric_match(m, *pi, *e, *qam, xref, d);

#ifdef DEBUG_COMPUTE_GEOMETRIC_MATCH
      cerr << "Embedding " << i << " fit " << fit << endl;
#endif

      if (fit < zbest || 0 == j) {
        zbest = fit;
        best_embedding = j;
      }
    }

    if (verbose > 1) {
      cerr << m.name() << " best match number " << best_embedding << " diff " << zbest
           << endl;
    }

    const Set_of_Atoms* e = sresults.embedding(best_embedding);
    const Query_Atoms_Matched* qam = sresults.query_atoms_matching(best_embedding);

    set_vector(xref, matoms, -1);

    fill_cross_reference_array(*e, *qam, xref);

// #define ECHO_CROSS_REFERENCE_ARRAY
#ifdef ECHO_CROSS_REFERENCE_ARRAY
    cerr << "Cross reference array filled\n";
    for (int i = 0; i < matoms; i++) {
      cerr << " i = " << i << " xref[i] " << xref[i] << endl;
    }
#endif

    transfer_bonding_information(m, xref, *pi);

    Geometric_Difference* gd;
    if (produce_text_descriptions) {
      gd = new Geometric_Difference[matoms * 4];
    } else {
      gd = nullptr;
    }

    pi->compute_aromaticity_if_needed();

    if (!xray_structure_compare(m, atom_names, *pi, xref, gd,
                                output)) {  // hard to imagine what could be wrong
      return 0;
    }

    if (stream_for_superimposed_molecules.active()) {
      superimpose_molecules(m, xref, *pi);
      m.add_molecule(pi);
      stream_for_superimposed_molecules.write(m);
    }

    if (nullptr != gd) {
      delete[] gd;
    }
  }

  if (stream_for_isotopically_labelled_molecules.active()) {
    atom_names.apply_isotopes(m);
    stream_for_isotopically_labelled_molecules.write(m);
  }

  return 1;
}

static void
preprocess(Molecule& m)
{
  m.remove_all(1);  // this is Xray, no hydrogens!

  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
break_into_spatially_separated_molecules(Molecule& parent,
                                         resizable_array_p<Molecule>& components)
{
  int matoms = parent.natoms();

  distance_t* d = new distance_t[matoms * matoms];
  std::unique_ptr<distance_t[]> free_d(d);

  for (int i = 0; i < matoms; i++) {
    for (int j = i + 1; j < matoms; j++) {
      distance_t dij = parent.distance_between_atoms(i, j);

      d[i * matoms + j] = d[j * matoms + i] = dij;
    }
  }

  int* fragment_membership = new_int(matoms);
  std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  int atoms_classified = 0;

  int next_frag_id = 1;

  while (atoms_classified < matoms) {
    for (int i = 0; i < matoms; i++) {
      if (fragment_membership[i] > 0) {
        continue;
      }

      fragment_membership[i] = -next_frag_id;

      int atoms_added = 0;

      while (1) {
        for (int j = 0; j < matoms; j++) {
          if (0 != fragment_membership[j]) {
            continue;
          }

          distance_t dij = d[i * matoms + j];

          if (dij > growth_range) {
            continue;
          }

          fragment_membership[j] = -next_frag_id;
          atoms_added++;
        }
      }

      fragment_membership[i] = next_frag_id;

      atoms_classified += atoms_added;

      if (0 == atoms_added) {
        break;
      }

      next_frag_id++;
    }
  }

  for (int i = 0; i < matoms; i++) {
    if (fragment_membership[i] < 0) {
      fragment_membership[i] = -fragment_membership[i];
    }
  }

  if (verbose > 1) {
    cerr << parent.name() << " split into " << next_frag_id << " fragments at dist "
         << growth_range << endl;
  }

  return 1;
}

static int
process_components(const IWString& mname, const resizable_array_p<Molecule>& components,
                   const Atom_Names& atom_names, const resizable_array_p<Molecule>& probe,
                   IWString_and_File_Descriptor& output)
{
  if (verbose) {
    cerr << "Input molecule split into " << components.number_elements()
         << " components\n";
  }

  for (int i = 0; i < components.number_elements(); i++) {
    Molecule* ci = components[i];
    ci->set_name(mname);

    xray_structure_compare2(*ci, atom_names, probe, output);
  }

  return 1;
}

static int
convert_all_bonds_to_single_bonds(Molecule& m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* ai = m.atomi(i);

    int icon = ai->ncon();

    for (int j = 0; j < icon; j++) {
      const Bond* b = ai->item(j);

      if (b->is_single_bond()) {
        continue;
      }

      atom_number_t k = b->other(i);

      m.set_bond_type_between_atoms(i, k, SINGLE_BOND);
    }
  }

  return 1;
}

/*
 */

static int
xray_structure_compare(Molecule& m, iwstring_data_source& atom_name_stream,
                       resizable_array_p<Molecule>& probe,
                       IWString_and_File_Descriptor& output)
{
  if (bonding_radius > 0.0) {
    add_bonds(m);
  }

  Atom_Names atom_names;

  if (!atom_name_stream.is_open()) {
    ;
  } else if (!atom_names.build(atom_name_stream)) {
    cerr << "Cannot read atom names\n";
    return 0;
  }

  if (directcolorfile_molecule.active()) {
    directcolorfile_molecule.write(probe[0]);
  }

  if (molecules_have_multiple_fragments) {
    ;
  } else if (m.number_fragments() > 1) {
    resizable_array_p<Molecule> components;
    m.create_components(components);

    return process_components(m.name(), components, atom_names, probe, output);
  }

  if (0.0 == growth_range) {
    return xray_structure_compare2(m, atom_names, probe, output);
  }

  resizable_array_p<Molecule> components;
  break_into_spatially_separated_molecules(m, components);

  return process_components(m.name(), components, atom_names, probe, output);
}

static int
xray_structure_compare(data_source_and_type<Molecule>& input,
                       iwstring_data_source& atom_name_stream,
                       resizable_array_p<Molecule>& probe,
                       IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    //  cerr << "After preprocessing '" << m->smiles() << "' '" << m->unique_smiles() <<
    //  endl;

    if (!xray_structure_compare(*m, atom_name_stream, probe, output)) {
      return 0;
    }
  }

  return output.good();
}

static int
xray_structure_compare(const char* fname, FileType input_type,
                       iwstring_data_source& atom_name_stream,
                       resizable_array_p<Molecule>& probe,
                       IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return xray_structure_compare(input, atom_name_stream, probe, output);
}

static int
read_probe_molecule(data_source_and_type<Molecule>& input,
                    resizable_array_p<Molecule>& probe)
{
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    m->remove_all(1);

    //  cerr << "After removing hydrogen '" << m->smiles() << "'\n";

    probe.add(m);
  }

  return probe.number_elements();
}

static int
read_probe_molecule(const const_IWSubstring& fname, resizable_array_p<Molecule>& probe,
                    FileType input_type)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "Cannot discern input type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_probe_molecule(input, probe);
}

static void
display_dash_G_options(std::ostream& os)
{
  os << " -G sbrms=<f>  single   bond distance rms, default " << single_bond_distance_rms
     << endl;
  os << " -G mbrms=<f>  multiple bond distance rms, default "
     << non_single_bond_distance_rms << endl;
  os << " -G bnsg=<i>   bond distance nsigma, default " << bond_distance_nsigma << endl;
  os << " -G arms=<f>   bond angle rms, default " << bond_angle_rms << endl;
  os << " -G ansg=<i>   bond angle nsigma, default " << bond_angle_nsigma << endl;
  os << " -G NONE       turn off all aspects of this feature\n";

  exit(0);
}

static int
xray_structure_compare(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lS:o:xp:f:b:t:rN:I:G:jD:cHX:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  } else if (cl.option_present('j')) {
    molecules_have_multiple_fragments = 1;

    if (verbose) {
      cerr << "Molecules can have multiple fragments - like a PDB with no bonds\n";
    }
  }

  if (cl.option_present('f')) {
    if (!cl.value('f', growth_range) || growth_range <= 0.0) {
      cerr << "The within fragment growth radius (-f) option must be a positive number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Atoms within " << growth_range
           << " of other fragment atoms groupted together\n";
    }
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', bonding_radius) || bonding_radius <= 0.0) {
      cerr << "The bonding radius (-b) option must be a positive number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Atoms within " << bonding_radius << " of each other will be bonded\n";
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', produce_text_descriptions) || produce_text_descriptions < 1) {
      cerr << "The produce text descriptions option (-t) must have a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will produce text descriptions of the " << produce_text_descriptions
           << " worst matches\n";
    }
  }

  if (cl.option_present('r')) {
    compute_initial_rms = 1;

    if (verbose) {
      cerr << "Will compute initial RMS\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (!process_cmdline_token('X', x, delocalise_bonds, verbose)) {
        cerr << "Cannot interpret delocalised atom query '" << x << "'\n";
        return 1;
      }
    }

    if (verbose) {
      cerr << "Defined " << delocalise_bonds.number_elements()
           << " queries for delocalised atoms\n";
    }
  }

  if (cl.option_present('G')) {
    int i = 0;
    const_IWSubstring g;
    while (cl.value('G', g, i++)) {
      if (g.starts_with("sbrms=")) {
        g.remove_leading_chars(6);
        if (!g.numeric_value(single_bond_distance_rms) ||
            single_bond_distance_rms < 0.0) {
          cerr << "The single bond distance rms value must be a valid float, '" << g
               << "' invalid\n";
          return 3;
        }

        if (verbose) {
          cerr << "Single bond distance rms " << single_bond_distance_rms << endl;
        }
      } else if (g.starts_with("mbrms=")) {
        g.remove_leading_chars(6);
        if (!g.numeric_value(non_single_bond_distance_rms) ||
            non_single_bond_distance_rms < 0.0) {
          cerr << "The non-single bond distance rms value must be a valid float, '" << g
               << "' invalid\n";
          return 3;
        }

        if (verbose) {
          cerr << "Non single bond distance rms " << non_single_bond_distance_rms << endl;
        }
      } else if (g.starts_with("bnsg=")) {
        g.remove_leading_chars(5);
        if (!g.numeric_value(bond_distance_nsigma) || bond_distance_nsigma <= 0.0) {
          cerr << "The bond distance nsigma value must be a valid float, '" << g
               << "' invalid\n";
          return 3;
        }

        if (verbose) {
          cerr << "Bond distance nsigma " << bond_distance_nsigma << endl;
        }
      } else if (g.starts_with("arms=")) {
        g.remove_leading_chars(5);
        if (!g.numeric_value(bond_angle_rms) || bond_angle_rms < 0.0) {
          cerr << "The bond angle rms value must be a valid float, '" << g
               << "' is invalid\n";
          return 2;
        }

        if (verbose) {
          cerr << "Bond angle rms set to " << bond_angle_rms << endl;
        }
      } else if (g.starts_with("ansg=")) {
        g.remove_leading_chars(5);
        if (!g.numeric_value(bond_angle_nsigma) || bond_angle_nsigma <= 0.0) {
          cerr << "The bond angle nsigma value must be a valid float, '" << g
               << "' invalid\n";
          return 3;
        }

        if (verbose) {
          cerr << "Bond angle nsigma " << bond_angle_nsigma << endl;
        }
      } else if ("NONE" == g) {
        single_bond_distance_rms = -1.0;
        non_single_bond_distance_rms = -1.0;
        bond_distance_nsigma = 0.0;
        bond_angle_rms = -1.0;
        bond_angle_nsigma = 0.0;

        if (verbose) {
          cerr << "Comparison with sigma variations turned off\n";
        }
      } else if ("help" == g) {
        display_dash_G_options(cerr);
      } else {
        cerr << "Unrecognised -G directive '" << g << "'\n";
        display_dash_G_options(cerr);
      }
    }
  }

  if (cl.option_present('x')) {
    relax_bonding_on_substructure_mismatch = 1;

    if (verbose) {
      cerr << "Will relax bonding information on substructure mismatch\n";
    }
  }

  if (cl.option_present('c')) {
    remove_chirality_on_substructure_search_failure = 1;

    if (verbose) {
      cerr << "Will remove chirality on substructure mismatch\n";
    }
  }

  if (cl.option_present('H')) {
    remove_explicit_hydrogens = 1;

    if (verbose) {
      cerr << "Will remove explicit Hydrogen atoms\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (!cl.option_present('p')) {
    cerr << "Must specify the probe molecule via the -p option\n";
    usage(2);
  }

  resizable_array_p<Molecule> probe;

  if (cl.option_present('p')) {
    const_IWSubstring p = cl.string_value('p');

    if (0 == read_probe_molecule(p, probe, input_type)) {
      cerr << "Cannot read probe molecule(s) from '" << p << "'\n";
      return 3;
    }

    for (int i = 0; i < probe.number_elements(); i++) {
      Molecule* pi = probe[i];

      preprocess(*pi);

      if (3 != pi->highest_coordinate_dimensionality()) {
        cerr << "Sorry, only works with 3D structures '" << pi->name() << "'\n";
        return 0;
      }
      if (bonding_radius > 0.0) {
        convert_all_bonds_to_single_bonds(*pi);
        if (verbose > 2) {
          cerr << "Probe is '" << pi->smiles() << "'\n";
        }
      }
    }

    if (verbose) {
      cerr << "Read " << probe.number_elements() << " probe molecules\n";
    }
  }
  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (!cl.option_present('o')) {
      stream_for_superimposed_molecules.add_output_type(FILE_TYPE_SDF);
    } else if (!stream_for_superimposed_molecules.determine_output_types(cl)) {
      cerr << "Cannot determine output specifications for -S file '" << s << "'\n";
      return 4;
    }

    if (stream_for_superimposed_molecules.would_overwrite_input_files(cl, s)) {
      cerr << "The -S file cannot overwrite any input(s) '" << s << "'\n";
      return 4;
    }

    if (!stream_for_superimposed_molecules.new_stem(s)) {
      cerr << "Cannot initialise stream for superimposed molecules '" << s << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Superimposed molecules written to '" << s << "'\n";
    }
  }

  iwstring_data_source atom_name_stream;

  if (cl.option_present('N')) {
    const char* n = cl.option_value('N');

    if (!atom_name_stream.open(n)) {
      cerr << "Cannot open atom name file '" << n << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Atom names read from '" << n << "'\n";
    }
  }

  if (cl.option_present('I')) {
    if (!cl.option_present('N')) {
      cerr << "The stream for isotopically labeled molecule (-I) only makes sense with "
              "the -N option\n";
      usage(2);
    }

    const_IWSubstring s = cl.string_value('I');

    stream_for_isotopically_labelled_molecules.add_output_type(FILE_TYPE_SMI);

    if (stream_for_isotopically_labelled_molecules.would_overwrite_input_files(cl, s)) {
      cerr << "The -I file cannot overwrite any input(s) '" << s << "'\n";
      return 4;
    }

    if (!stream_for_isotopically_labelled_molecules.new_stem(s)) {
      cerr << "Cannot initialise stream for isotopically labelled molecules '" << s
           << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Isotopically labelled molecules written to '" << s << "'\n";
    }
  }

  initialise_max_con();

  set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present('D')) {
    const char* d = cl.option_value('D');

    directcolorfile_molecule.add_output_type(FILE_TYPE_SDF);
    if (directcolorfile_molecule.would_overwrite_input_files(cl, d)) {
      cerr << "Directcolorfile '" << d << "' cannot overwrite input file(s)\n";
      exit(2);
    }

    if (!directcolorfile_molecule.new_stem(d)) {
      cerr << "Cannot initialise directcolorfile_molecule '" << d << "'\n";
      return 2;
    }

    IWString dtmp(d);
    if (!dtmp.ends_with(".dcf")) {
      dtmp << ".dcf";
    }

    if (!directcolorfile_colour.open(dtmp.null_terminated_chars())) {
      cerr << "Cannot open directcolorfile_colour '" << d << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Initialised directcolorfiles '" << d << "'\n";
    }

    nsigma_colour_array.add(new IWString("blue"));
    nsigma_colour_array.add(new IWString("purple"));
    nsigma_colour_array.add(new IWString("orange"));
    nsigma_colour_array.add(new IWString("red"));
  }

  IWString_and_File_Descriptor output(1);

  int rc = xray_structure_compare(cl[0], input_type, atom_name_stream, probe, output);

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return 0 != rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = xray_structure_compare(argc, argv);

  return rc;
}
