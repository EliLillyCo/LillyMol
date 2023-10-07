/*
  Tester for the shadow code
*/
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/jama/jama_qr.h"
#include "Foundational/jama/jama_svd.h"
#include "Foundational/tnt/tnt_array2d_utils.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iw_vdw.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "shadow.h"

#include "Eigen/Core"
#include "Eigen/Dense"

using std::cerr;
using std::endl;
// #define DEBUG_RADIUS_OF_GYRATION_STUFF

#ifdef DEBUG_RADIUS_OF_GYRATION_STUFF
#include "output.h"
#endif

static int verbose = 0;

static int molecules_read = 0;

static float grid_resolution = 0.0;

static int compute_radius_of_gyration = 0;

static int write_ordered_areas = 0;

static int vdw_type = IW_VDW_SAVOL;

static double multiply_vdw_radii = 1.0;

static Accumulator<area_t> acc;

/*
  May 2003. Maybe we only want to consider some of the atoms when doing the alignment
*/

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*
  Or maybe we define two queries which define the two ends. First atom matches
*/

static resizable_array_p<Substructure_Hit_Statistics> end1queries;
static resizable_array_p<Substructure_Hit_Statistics> end2queries;

/*
  Normally we die if we cannot assign vdw radii. Optionally, we can just
  ignore the molecule
*/

static int ok_if_no_vdw_radii = 0;

/*
  We have one stream for each type of orientation requested
*/

static Molecule_Output_Object* stream_for_oriented_molecules = nullptr;

static int write_oriented_molecule[] = {0, 0, 0, 0};

/*
  print_grid will be 2 if we also want to see the perpendicular grid
*/

static int print_grid = 0;

static char on = '@';
static char off = '-';

static int rotate_to_longest_distance = 0;

#define OUTPUT_PRECISION 4

static IWString tag;

static IWString smiles_tag("$SMI<");

static int max_bit_count = 10;

typedef double tnt_float_type;

static int align_by_principal_components = 0;

static int pbf_only = 0;  // temporary capability, remove eventually

static int rdkit_only = 0;

static int check_rdkit = 0;

static int remove_hydrogens = 1;

static Accumulator<double> rdkit_diff;

static int position_by_extremeties = 1;

/*
  Run the computation across the Lilly collection and get extremeties
*/

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on

  // clang-format off
  cerr << "  -P <ext>       do alignment by principal components, extremeties equi-distant from origin\n";
  cerr << "  -P <ave>       do alignment by principal components, average spatial offsets equal around origin\n";
  cerr << "  -r <dist>      set the grid resolution to <dist> Angstroms\n";
  cerr << "  -p on=X        print the grid and use X as the on  character\n";
  cerr << "  -p off=x       print the grid and use x as the off character\n";
  cerr << "  -p def         print the grid and use the defaults: on='" << on << "', off='" << off << "'\n";
  cerr << "  -p Y           print the grid of the Y axis rotated orientation\n";
  cerr << "  -L             rotate to longest distance between atoms\n";
  cerr << "  -G             compute radius of gyration\n";
  cerr << "  -b             only compute PBF descriptor\n";
  cerr << "  -t             test mode, check internal computation against RDKIT implementation\n";
  cerr << "  -k             use only the RDKIT implementation\n";
  cerr << "  -h             do NOT remove hydrogen atoms\n";
  (void) display_standard_vdw_radius_types (cerr, 'V');
  cerr << "  -m <mult>      multiply vdw radii by factor\n";
  cerr << "  -u             just skip molecules for which VDW radius assignment fails\n";
  cerr << "  -d             also write shadow values as ordered set\n";
  cerr << "  -q <query>     substructure query. Only matched atoms considered when doing alignment\n";
  cerr << "  -s <smarts>    same as -q option\n";
  cerr << "  -e <query>     Define extremeties. Must be exactly two -e options\n";
  cerr << "  -J <tag>       work as a TDT filter\n";
  cerr << "  -y <count>     when producing fingerprints, max count\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -O <fname>     file name for oriented molecules - sdf output only\n";
  cerr << "  -E ...         standard element options, enter '-E help' for details\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
write_header(IWString_and_File_Descriptor& output) {
  output.resize(64000);

  if (pbf_only) {
    output << "Name PBF";
    if (check_rdkit) {
      output << " RDKIT Diff";
    }
    output << '\n';
    return 1;
  }

  output << "Name sh_maxdist sh_avedist sh_maxhdist sh_avehdist sh_shadow1";
  output << " sh_shadow2 sh_shadow3";

  if (write_ordered_areas) {
    output << " sh_shbig sh_shmid sh_shmin";
  }

  output << " sh_cigar1 sh_cigar13 sh_cigar23 sh_boxy sh_boxz sh_boxvol sh_sphere "
            "sh_nsphere sh_msphere sh_flat sh_nflat sh_mflat sh_tube sh_ntube sh_mtube";
  output << " sh_bow";
  output << " sh_zanisotropy sh_abszanisotropy";

  if (compute_radius_of_gyration) {
    output << " sh_rgspatl sh_rgmasswt";
  }

  output << '\n';

  return output.size();
}

static void
do_position_by_extremeties(Molecule& m, TNT::Array2D<tnt_float_type>& geometry) {
  coord_t xmin, xmax, ymin, ymax, zmin, zmax;

  m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);

  coord_t xmid = (xmax + xmin) * 0.5;
  coord_t ymid = (ymax + ymin) * 0.5;
  coord_t zmid = (zmax + zmin) * 0.5;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    Coordinates c;
    m.get_coords(i, c);

    coord_t newx = c.x() - xmid;
    coord_t newy = c.y() - ymid;
    coord_t newz = c.z() - zmid;

    m.setxyz(i, newx, newy, newz);

    geometry[i][0] = newx;
    geometry[i][1] = newy;
    geometry[i][2] = newz;

    //  cerr << newx << ' ' << newy << ' ' << newz << endl;
  }

  return;
}

static void
do_position_by_averages(Molecule& m, TNT::Array2D<tnt_float_type>& geometry) {
  Accumulator<double> x, y, z;

  const auto matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    Coordinates c;
    m.get_coords(i, c);

    x.extra(c.x());
    y.extra(c.y());
    z.extra(c.z());
  }

  double xave = x.average();
  double yave = y.average();
  double zave = z.average();

  for (int i = 0; i < matoms; i++) {
    Coordinates c;
    m.get_coords(i, c);

    coord_t newx = c.x() - xave;
    coord_t newy = c.y() - yave;
    coord_t newz = c.z() - zave;

    m.setxyz(i, newx, newy, newz);

    geometry[i][0] = newx;
    geometry[i][1] = newy;
    geometry[i][2] = newz;
  }

  return;
}

#ifdef DEBUG_RADIUS_OF_GYRATION_STUFF

static Molecule_Output_Object stream_for_centres;

#endif

static int
do_compute_radius_of_gyration(Molecule& m, const Atom** atom, coord_t& rg1,
                              coord_t& rg2) {
#ifdef DEBUG_RADIUS_OF_GYRATION_STUFF
  if (!stream_for_centres.active()) {
    stream_for_centres.add_output_type(SDF);

    if (!stream_for_centres.new_stem("radius_of_gyration")) {
      cerr << "Gack, cannot open radius of gyration file\n";
      return 0;
    }
  }
#endif

  int matoms = m.natoms();

  Coordinates spatial_centre;
  Coordinates mass_weighted_centre;
  atomic_mass_t total_mass = 0.0;

  for (int i = 0; i < matoms; i++) {
    const Atom& ai = *(atom[i]);

    spatial_centre += ai;

    //  mass weighted is more difficult

    atomic_mass_t a = ai.element()->atomic_mass();  // no isotopes considered

    Coordinates ci = ai;
    ci *= a;

    mass_weighted_centre += ci;
    total_mass += a;
  }

  spatial_centre /= static_cast<coord_t>(matoms);
  mass_weighted_centre /= static_cast<coord_t>(total_mass);

  double r1 = 0.0;
  double r2 = 0.0;

  for (int i = 0; i < matoms; i++) {
    const Coordinates& ai = *(atom[i]);

    double d = spatial_centre.distance_squared(ai);

    r1 += d;

    d = mass_weighted_centre.distance_squared(ai);

    r2 += d;
  }

  rg1 = sqrt(r1 / static_cast<double>(matoms));
  rg2 = sqrt(r2 / static_cast<double>(matoms));

#ifdef DEBUG_RADIUS_OF_GYRATION_STUFF
  Molecule mcopy = m;

  Atom* a1 = new Atom("Se");  // an element starting with S
  a1->setxyz(spatial_centre);
  mcopy.add(a1);

  Atom* a2 = new Atom("Mn");  // an element starting with M
  a2->setxyz(mass_weighted_centre);
  mcopy.add(a2);

  stream_for_centres.write(mcopy);
#endif

  return 1;
}

#ifdef USEFUL_DIAGNOSTIC
static void
report_average_position(const Molecule& m, const char* msg, std::ostream& output) {
  Accumulator<double> xx, yy;
  for (auto i = 0; i < m.natoms(); ++i) {
    xx.extra(m.x(i));
    yy.extra(m.y(i));
  }
  output << msg << ", x = " << xx.average() << " y " << yy.average() << endl;

  return;
}
#endif

static int
do_write_to_stream_for_oriented_molecules(Molecule& m, int sequence) {
  if (0 == write_oriented_molecule[sequence]) {
    return 1;
  }

  if (!stream_for_oriented_molecules[sequence].write(m)) {
    return 0;
  }

  return 1;
}

/*
  The molecule is aligned along the X axis. We need to make sure it is
  somewhat flat within the XY plane
*/

static void
do_rotate_into_xy_plane(Molecule& m, const Atom* const* atom,
                        const int* process_these_atoms) {
  int matoms = m.natoms();

  coord_t maxd = static_cast<coord_t>(0.0);
  atom_number_t furthest_atom = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; i++) {
    if (0 == process_these_atoms[i]) {
      continue;
    }

    const Atom* a = atom[i];

    coord_t d = a->y() * a->y() + a->z() * a->z();

    if (d > maxd) {
      maxd = d;
      furthest_atom = i;
    }
  }

  if (INVALID_ATOM_NUMBER == furthest_atom) {  // how could that happen?
    return;
  }

  // Now rotate so FURTHEST_ATOM is in the X/Y plane

  const Atom* a = atom[furthest_atom];

  Coordinates c(static_cast<coord_t>(0.0), a->y(), a->z());
  c.normalise();

  Coordinates yaxis(static_cast<coord_t>(0.0), static_cast<coord_t>(1.0),
                    static_cast<coord_t>(0.0));

  angle_t theta = yaxis.angle_between_unit_vectors(c);

#ifdef DEBUG_ROTATE_ABOUT_Y_AZIS
  cerr << "Furthest atom starts " << a->x() << " " << a->y() << " " << a->z() << endl;
  cerr << "Angle to Y axis " << (theta * RAD2DEG) << endl;
#endif

  Coordinates xaxis(static_cast<coord_t>(1.0), static_cast<coord_t>(0.0),
                    static_cast<coord_t>(0.0));

  if (a->z() < static_cast<coord_t>(0.0)) {
    m.rotate_atoms(xaxis, theta);
  } else {
    m.rotate_atoms(xaxis, -theta);
  }

#ifdef DEBUG_ROTATE_ABOUT_Y_AZIS
  cerr << "Furthest atom now " << a->x() << " " << a->y() << " " << a->z() << endl;
#endif

  // report_average_position(m, "After Y axis Rotation", cerr);

#ifdef DEBUG_ROTATE_ABOUT_Y_AZIS
  cerr << "Furthest atom now " << a->x() << " " << a->y() << " " << a->z() << endl;
#endif

  // assert (fabs(a->x()) < 0.01);

  return;
}

#ifdef NOT_USED
static void
do_rotate_about_z_axis(Molecule& m) {
  Coordinates zaxis(static_cast<coord_t>(0.0), static_cast<coord_t>(0.0),
                    static_cast<coord_t>(1.0));

  m.rotate_atoms(zaxis, M_PI * 0.5);

  return;
}
#endif

/*
  Are two atoms a pair of heteroatoms?
*/

static int
is_heteroatom_pair(const Atom* const* atom, int i1, int i2) {
  if (6 == atom[i1]->atomic_number() || 1 == atom[i1]->atomic_number()) {
    return 0;
  }

  if (6 == atom[i2]->atomic_number() || 1 == atom[i2]->atomic_number()) {
    return 0;
  }

  return 1;
}

/*
  To avoid passing around lots of arguments, we hold a bunch of information in an object
*/

class Shadow_Tmp {
 private:
  int* _process_these_atoms;

  vdw_radius_t* _vdw;

  Accumulator<coord_t> _inter_atomic_distance;

  Accumulator<coord_t> _inter_heteroatom_distance;

  //  Accumulator<float> _distance_to_bond_count_ratio;

 public:
  Shadow_Tmp(int);
  ~Shadow_Tmp();

  int*
  process_these_atoms() {
    return _process_these_atoms;
  }

  const int*
  process_these_atoms() const {
    return _process_these_atoms;
  }

  const Accumulator<coord_t>&
  inter_atomic_distance() const {
    return _inter_atomic_distance;
  }

  const Accumulator<coord_t>&
  inter_heteroatom_distance() const {
    return _inter_heteroatom_distance;
  }

  //  Accumulator<coord_t> & distance_to_bond_count_ratio() { return
  //  _distance_to_bond_count_ratio;} const Accumulator<coord_t> &
  //  distance_to_bond_count_ratio() const { return _distance_to_bond_count_ratio;}

  const vdw_radius_t*
  vdw() const {
    return _vdw;
  }

  int
  assign_vdw_radii(Molecule&, int);

  void
  determine_extremeties(const Molecule&, const Atom** const);
  void
  do_rotate_to_longest_distance(Molecule& m, const Atom** atom);
};

Shadow_Tmp::Shadow_Tmp(int matoms) {
  _process_these_atoms = new_int(matoms, 1);

  _vdw = new vdw_radius_t[matoms];

  return;
}

Shadow_Tmp::~Shadow_Tmp() {
  delete[] _process_these_atoms;

  delete[] _vdw;

  return;
}

int
Shadow_Tmp::assign_vdw_radii(Molecule& m, int vdw_type) {
  if (!::assign_vdw_radii(m, vdw_type, _vdw)) {
    return 0;
  }

  if (1.0 != multiply_vdw_radii) {
    int matoms = m.natoms();

    for (int i = 0; i < matoms; i++) {
      _vdw[i] *= multiply_vdw_radii;
    }
  }

  return 1;
}

#ifdef NOT_USED
static double
compute_pbf(const Molecule& m, const double eqa, const double eqb, const double eqc,
            const double eqd) {
  const double norm = 1.0 / sqrt(eqa * eqa + eqb * eqb + eqc * eqc);

  const auto matoms = m.natoms();

  double pbf = 0.0;

  for (auto i = 0; i < matoms; ++i) {
    Coordinates c;
    m.get_coords(i, c);

    double d = (eqa * c.x() + eqb * c.y() + eqc * c.z() + eqd) * norm;

    pbf += fabs(d);

#ifdef DEBUG_PBF
//  cerr << "Atom " << i << " at (" << c.x() << ',' << c.y() << ',' << c.z() << ") dist "
//  << d << endl;
#endif
  }

  return pbf / static_cast<double>(matoms);
}
#endif

#ifdef NOT_USED
static void
assemble_A_matrix_and_B_vector(const Molecule& m, TNT::Array2D<tnt_float_type>& a,
                               TNT::Array1D<tnt_float_type>& b) {
  const auto matoms = m.natoms();

  double x2 = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double y2 = 0.0;
  double yz = 0.0;
  double z2 = 0.0;

  Coordinates c;

  for (auto i = 0; i < matoms; ++i) {
    m.get_coords(i, c);
    double cx = c.x();
    double cy = c.y();
    double cz = c.z();
    x2 += cx * cx;
    xy += cx * cy;
    xz += cx * cz;
    y2 += cy * cy;
    yz += cy * cz;
    z2 += cz * cz;
  }

  a[0][0] = x2;
  a[0][1] = xy;
  a[0][2] = xz;
  a[1][0] = xy;
  a[1][1] = y2;
  a[1][2] = yz;
  a[2][0] = xz;
  a[2][1] = yz;
  a[2][2] = z2;

  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;

  return;
}

static void
assemble_A_matrix_and_B_vector(const Molecule& m, Eigen::Matrix3d& a,
                               Eigen::Vector3d& b) {
  const auto matoms = m.natoms();

  double x2 = 0.0;
  double xy = 0.0;
  double x = 0.0;
  double y2 = 0.0;
  double y = 0.0;
  double xz = 0.0;
  double yz = 0.0;
  double z = 0.0;

  Coordinates c;

  for (auto i = 0; i < matoms; ++i) {
    m.get_coords(i, c);
    double cx = c.x();
    double cy = c.y();
    double cz = c.z();
    x2 += cx * cx;
    xy += cx * cy;
    x += cx;
    y2 += cy * cy;
    y += cy;
    xz += cx * cz;
    yz += cy * cz;
    z += cz;
  }

  a(0, 0) = x2;
  a(0, 1) = xy;
  a(0, 2) = x;
  a(1, 0) = xy;
  a(1, 1) = y2;
  a(1, 2) = y;
  a(2, 0) = x;
  a(2, 1) = y;
  a(2, 2) = static_cast<double>(matoms);

  b(0) = xz;
  b(1) = yz;
  b(2) = z;

  return;
}
#endif

static bool
getBestFitPlane(const std::vector<Coordinates>& points, std::vector<double>& plane,
                const std::vector<double>* weights) {
  Coordinates origin(0, 0, 0);
  double wSum = 0.0;

  for (unsigned int i = 0; i < points.size(); ++i) {
    if (weights) {
      double w = (*weights)[i];
      wSum += w;
      origin += points[i] * w;
    } else {
      wSum += 1;
      origin += points[i];
    }
  }
  origin /= wSum;

  double sumXX = 0, sumXY = 0, sumXZ = 0, sumYY = 0, sumYZ = 0, sumZZ = 0;
  for (unsigned int i = 0; i < points.size(); ++i) {
    Coordinates delta = points[i] - origin;
    if (weights) {
      double w = (*weights)[i];
      delta *= w;
    }
    sumXX += delta.x() * delta.x();
    sumXY += delta.x() * delta.y();
    sumXZ += delta.x() * delta.z();
    sumYY += delta.y() * delta.y();
    sumYZ += delta.y() * delta.z();
    sumZZ += delta.z() * delta.z();
  }
  sumXX /= wSum;
  sumXY /= wSum;
  sumXZ /= wSum;
  sumYY /= wSum;
  sumYZ /= wSum;
  sumZZ /= wSum;

  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ, sumXY, sumYY, sumYZ, sumXZ, sumYZ, sumZZ;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
  if (eigensolver.info() != Eigen::Success) {
    cerr << "eigenvalue calculation did not converge" << std::endl;
    return 0.0;
  }
  Coordinates normal;
  normal.x() = eigensolver.eigenvectors()(0, 0);
  normal.y() = eigensolver.eigenvectors()(1, 0);
  normal.z() = eigensolver.eigenvectors()(2, 0);

  plane[0] = normal.x();
  plane[1] = normal.y();
  plane[2] = normal.z();
  plane[3] = -1 * normal.dot_product(origin);

  return true;
}

static double
distanceFromAPlane(const Coordinates& pt, const std::vector<double>& plane,
                   double denom) {
  double numer = 0.0;
  numer = std::abs(pt.x() * plane[0] + pt.y() * plane[1] + pt.z() * plane[2] + plane[3]);

  return numer / denom;
}

static double
PBFRD(const Molecule& mol) {
  const auto numAtoms = mol.natoms();
  if (numAtoms < 4) {
    return 0;
  }

  std::vector<Coordinates> points;
  points.reserve(numAtoms);
  for (int i = 0; i < numAtoms; ++i) {
    Coordinates c;
    mol.get_coords(i, c);
    points.push_back(c);
  }

  std::vector<double> plane(4);
  getBestFitPlane(points, plane, 0);

  double denom = 0.0;
  for (unsigned int i = 0; i < 3; ++i) {
    denom += plane[i] * plane[i];
  }
  denom = pow(denom, 0.5);

  double res = 0.0;
  for (int i = 0; i < numAtoms; ++i) {
    res += distanceFromAPlane(points[i], plane, denom);
  }
  res /= static_cast<double>(numAtoms);

  return res;
}

static int
do_align_by_principal_components(Molecule& m) {
  const auto matoms = m.natoms();

  // First centre the molecule

  TNT::Array2D<tnt_float_type> geometry(matoms, 3);

  if (position_by_extremeties) {
    do_position_by_extremeties(m, geometry);
  } else {
    do_position_by_averages(m, geometry);
  }

#ifdef ECHO_MATRICES
  cerr << "Starting matrix\n";
  cerr << geometry;
#endif

  JAMA::SVD<tnt_float_type> svd(geometry);

  // TNT::Array1D<tnt_float_type> singular_values;

  // svd.getSingularValues(singular_values);

#ifdef ECHO_MATRICES
  cerr << "TNT PCA singular values " << singular_values[0] << " " << singular_values[1]
       << ' ' << singular_values[2] << endl;
#endif

  TNT::Array2D<tnt_float_type> v;
  svd.getV(v);

  TNT::Array2D<tnt_float_type> t2 = matmult(geometry, v);

  for (int i = 0; i < matoms; i++) {
    m.setxyz(i, t2[i][0], t2[i][1], t2[i][2]);
    //  cerr << " atom " << i << " x " << t2[i][0] << " y " << t2[i][1] << " z " <<
    //  t2[i][2] << " sum " << pbf << endl;
  }

  return 1;
}

#ifdef NOT_USED
static void
move_to_origin(Molecule& m) {
  coord_t xmin, xmax, ymin, ymax, zmin, zmax;

  m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);

  m.translate_atoms(-(xmin + xmax) * 0.5, -(ymin + ymax) * 0.5, -(zmax - zmin) * 0.5);

  return;
}
#endif

#ifdef ANOTHER_WAY_OF_ALIGNING
static int
do_align_by_plane_of_best_fit_qr(Molecule& m) {
  move_to_origin(m);

  TNT::Array2D<tnt_float_type> a(3, 3);
  TNT::Array1D<tnt_float_type> b(3);

  assemble_A_matrix_and_B_vector(m, a, b);

#ifdef DEBUG_PBF
  cerr << a;
#endif

  JAMA::QR<tnt_float_type> qr(a);

  if (!qr.isFullRank()) {  // already aligned
    return 1;
  }

  TNT::Array1D<tnt_float_type> s = qr.solve(b);

  double eqa = s[0];
  double eqb = s[1];
  double eqc = -1.0;
  double eqd = s[2];

#ifdef DEBUG_PBF
  cerr << "Equation is A = " << eqa << " B = " << eqb << " C = " << eqc << " D = " << eqd
       << endl;
#endif

  const auto matoms = m.natoms();

  double pbf = compute_pbf(m, eqa, eqb, eqc, eqd);

  cerr << m.name() << ' ' << pbf / static_cast<double>(matoms) << endl;

  auto norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);

  // cerr << "Valyes are " << s[0] << ' ' << s[1] << ' ' << s[2] << ", norm " << norm <<
  // ", dim " << s.dim() << endl;

  if (norm < 1.0e-08) {
    cerr << "Best fit solution is zero\n";
    return 0;
  }

  double tmp = 0.0;
  for (auto i = 0; i < 3; ++i) {
    s[i] = s[i] / norm;
    tmp += s[i] * s[i];
  }

  // cerr << "Valyes are " << s[0] << ' ' << s[1] << ' ' << s[2] << " norm " << tmp <<
  // endl;

  return 1;
}
#endif

void
Shadow_Tmp::determine_extremeties(const Molecule& m, const Atom** atom) {
  const auto matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (0 == _process_these_atoms[i]) {
      continue;
    }

    const Atom* ai = atom[i];

    for (int j = i + 1; j < matoms; j++) {
      if (0 == _process_these_atoms[j]) {
        continue;
      }

      if (ai->is_bonded_to(j)) {
        continue;
      }

      const Atom* aj = atom[j];

      coord_t d = ai->distance(*aj) + _vdw[i] + _vdw[j];

      _inter_atomic_distance.extra(d);

      if (is_heteroatom_pair(atom, i, j)) {
        _inter_heteroatom_distance.extra(d);
      }
    }
  }

  return;
}

void
Shadow_Tmp::do_rotate_to_longest_distance(Molecule& m, const Atom** atom) {
  atom_number_t left = INVALID_ATOM_NUMBER;
  atom_number_t right = INVALID_ATOM_NUMBER;

  coord_t maxd = static_cast<coord_t>(0.0);

  const auto matoms = m.natoms();

  if (2 == matoms) {
    left = 0;
    right = 1;

    maxd = atom[0]->distance(*(atom[1]));

    _inter_atomic_distance.extra(maxd);

    if (is_heteroatom_pair(atom, 0, 1)) {
      _inter_heteroatom_distance.extra(maxd);
    }
  } else {
    for (auto i = 0; i < matoms; i++) {
      if (0 == _process_these_atoms[i]) {
        continue;
      }

      const Atom* ai = atom[i];

      for (int j = i + 1; j < matoms; j++) {
        if (0 == _process_these_atoms[j]) {
          continue;
        }

        //      if (ai->is_bonded_to(j))    // omit this test, Feb 2004. Why was it there?
        //        continue;

        const Atom* aj = atom[j];

        coord_t d = ai->distance(*aj) + _vdw[i] + _vdw[j];

        _inter_atomic_distance.extra(d);

        if (is_heteroatom_pair(atom, i, j)) {
          _inter_heteroatom_distance.extra(d);
        }

        if (d > maxd) {
          maxd = d;
          left = i;
          right = j;
        }
      }
    }
  }

  // In order to enforce more uniform behaviour, canonicalise left and right.

  if (m.x(left) > m.x(right)) {
    std::swap(left, right);
  }

  assert(INVALID_ATOM_NUMBER != left);
  assert(INVALID_ATOM_NUMBER != right);

  if (verbose > 2) {
    const Atom* a = m.atomi(left);

    cerr << m.name() << " extremity atoms " << left << " "
         << m.smarts_equivalent_for_atom(left) << ' ' << a->x() << ' ' << a->y() << ' '
         << a->z() << endl;

    a = m.atomi(right);
    cerr << " and " << right << ' ' << m.smarts_equivalent_for_atom(right) << ' '
         << a->x() << ' ' << a->y() << ' ' << a->z() << " dist " << maxd << endl;
  }

  // Translate the molecule so that atom LEFT is at the origin

  m.translate_atoms(-*(atom[left]));

  assert(static_cast<coord_t>(0.0) == atom[left]->x() &&
         static_cast<coord_t>(0.0) == atom[left]->y() &&
         static_cast<coord_t>(0.0) == atom[left]->z());

  // We now want to rotate the molecule so that RIGHT is along the X axis

  Coordinates x(static_cast<coord_t>(1.0), static_cast<coord_t>(0.0),
                static_cast<coord_t>(0.0));
  Coordinates r(atom[right]->x(), atom[right]->y(), atom[right]->z());
  r.normalise();
  angle_t theta = x.angle_between_unit_vectors(r);

  r.cross_product(x);
  r.normalise();

// #define DEBUG_ROTATE_TO_LOGNEST_DISTANCE
#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  cerr << "Right starts " << atom[right]->x() << ' ' << atom[right]->y() << ' '
       << atom[right]->z() << " angle " << (theta * RAD2DEG) << endl;
#endif

  m.rotate_atoms(r, theta);

#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  cerr << "Right now    " << atom[right]->x() << ' ' << atom[right]->y() << ' '
       << atom[right]->z() << " angle " << (theta * RAD2DEG) << endl;
#endif

  // At this stage, the molecule should have LEFT at the origin and RIGHT somwhere along
  // the X axis

  assert(static_cast<coord_t>(0.0) == atom[left]->x() &&
         static_cast<coord_t>(0.0) == atom[left]->y() &&
         static_cast<coord_t>(0.0) == atom[left]->z());
  assert(fabs(atom[right]->y()) < static_cast<coord_t>(0.02));
  assert(fabs(atom[right]->z()) < static_cast<coord_t>(0.02));

#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  report_average_position(m, "At end of rotate_2", cerr);
  cerr << "left is " << left << " (" << atom[left]->x() << ',' << atom[left]->y()
       << "), right " << right << " (" << atom[right]->x() << ',' << atom[right]->y()
       << ")\n";
#endif

  return;
}

/*
  Set the member of PROCESS_THESE_ATOMS that is the first matched atom of the queries
*/

static int
find_first_match(Molecule_to_Match& target,
                 resizable_array_p<Substructure_Hit_Statistics>& queries,
                 int* process_these_atoms) {
  for (int i = 0; i < queries.number_elements(); i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    const Set_of_Atoms* e = sresults.embedding(0);

    process_these_atoms[e->item(0)] = 1;
    return 1;
  }

  return 0;
}

static void
do_rotate_to_longest_distance_end_queries(Molecule& m, const Atom** atom,
                                          Shadow_Tmp& stmp) {
  int* process_these_atoms = stmp.process_these_atoms();
  set_vector(process_these_atoms, m.natoms(), 0);

  Molecule_to_Match target(&m);

  if (!find_first_match(target, end1queries, process_these_atoms)) {
    cerr << m.name() << " no match to end1 queries, default processing\n";
    set_vector(process_these_atoms, m.natoms(), 1);
    stmp.do_rotate_to_longest_distance(m, atom);
    return;
  }

  if (!find_first_match(target, end2queries, process_these_atoms)) {
    cerr << m.name() << " no match to end2 queries, default processing\n";
    set_vector(process_these_atoms, m.natoms(), 1);
  }

  stmp.do_rotate_to_longest_distance(m, atom);

  return;
}

static void
do_rotate_to_longest_distance(Molecule& m, const Atom** atom, Shadow_Tmp& stmp) {
  if (end1queries.number_elements() > 0) {
    do_rotate_to_longest_distance_end_queries(m, atom, stmp);
    return;
  }

  if (queries.empty()) {
    stmp.do_rotate_to_longest_distance(m, atom);
    return;
  }

  int matoms = m.natoms();

  int* process_these_atoms = stmp.process_these_atoms();
  set_vector(process_these_atoms, matoms, 0);

  Molecule_to_Match target(&m);

  int queries_matching = 0;

  for (int i = 0; i < queries.number_elements(); i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    queries_matching++;

    sresults.each_embedding_set_vector(process_these_atoms, 1);
  }

  if (0 == queries_matching) {
    cerr << "Warning, none of " << queries.number_elements()
         << " queries, match, default processing\n";
    set_vector(process_these_atoms, matoms, 1);
  }

  stmp.do_rotate_to_longest_distance(m, atom);

  return;
}

/*
  First rotate the molecule to lie in the XY plane.
  Remember that the shadow class does projections onto the YZ plane, so the order in which
  the values are computed might not necessarily correspond with the orientations here
*/

static int
tshadow(Molecule& m, const Atom** atom, Shadow_Tmp& stmp,
        resizable_array<float>& output) {
  if (!stmp.assign_vdw_radii(m, vdw_type)) {
    cerr << "Cannot assign VDW radii to " << molecules_read << " '" << m.name() << "'\n";
    return ok_if_no_vdw_radii;
  }

  if (align_by_principal_components) {
    do_align_by_principal_components(m);
    stmp.determine_extremeties(m, atom);
    do_write_to_stream_for_oriented_molecules(m, 1);
  } else if (!rotate_to_longest_distance) {
    stmp.determine_extremeties(m, atom);
  } else {
    do_rotate_to_longest_distance(m, atom, stmp);

    //  report_average_position(m, "After longest distance rotation", cerr);

    do_rotate_into_xy_plane(m, atom, stmp.process_these_atoms());

    do_write_to_stream_for_oriented_molecules(m, 1);
  }

  const vdw_radius_t* vdw = stmp.vdw();

  shadow::Shadow_Area sa;

  if (grid_resolution > 0.0) {
    sa.set_resolution(grid_resolution);
  }

  auto a1 = sa.shadow_area(m, vdw);

  if (verbose) {
    acc.extra(a1);
  }

  if (print_grid) {
    cerr << "Molecule has " << m.natoms() << " atoms, a1 " << a1 << endl;
    sa.print_grid(cerr, on, off);
  }

  const auto matoms = m.natoms();

  Accumulator<distance_t> out_of_plane, abs_out_of_plane;

  for (auto i = 0; i < matoms; ++i) {
    coord_t z = m.z(i);
    out_of_plane.extra(z);
    abs_out_of_plane.extra(fabs(z));
  }

  coord_t xmin, xmax, ymin, ymax, zmin, zmax;

  m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);

  float box_x = static_cast<float>(xmax - xmin);
  float box_y = static_cast<float>(ymax - ymin);
  float box_z = static_cast<float>(zmax - zmin);

  Accumulator<double> zleft, zright, abs_zleft, abs_zright, zz;

  coord_t xmid = (xmax + xmin) * static_cast<coord_t>(0.50);

  for (auto i = 0; i < matoms; ++i) {
    const Atom* ai = atom[i];

    //  cerr << "Examining atom " << i << " (" << ai->x() << ',' << ai->y() << ',' <<
    //  ai->z() << ")\n";

    zz.extra(fabs(ai->z()));

    if (ai->x() < xmid) {
      zleft.extra(ai->z());
      abs_zleft.extra(fabs(ai->z()));
    } else {
      zright.extra(ai->z());
      abs_zright.extra(fabs(ai->z()));
    }
  }

  Coordinates yaxis(static_cast<coord_t>(0.0), static_cast<coord_t>(1.0),
                    static_cast<coord_t>(0.0));

  m.rotate_atoms(yaxis, M_PI * 0.5);

  auto a2 = sa.shadow_area(m, vdw);

  do_write_to_stream_for_oriented_molecules(m, 2);

  if (print_grid > 1) {
    cerr << "After rotation about Y, area " << a2 << endl;
    sa.print_grid(cerr, on, off);
  }

  Coordinates zaxis(static_cast<coord_t>(0.0), static_cast<coord_t>(0.0),
                    static_cast<coord_t>(1.0));

  m.rotate_atoms(zaxis, M_PI * 0.5);

  do_write_to_stream_for_oriented_molecules(m, 3);

  auto a3 = sa.shadow_area(m, vdw);

  if (print_grid > 1) {
    cerr << "After rotation about Z, area " << a3 << endl;
    sa.print_grid(cerr, on, off);
  }

  coord_t rg1, rg2;

  if (compute_radius_of_gyration) {
    do_compute_radius_of_gyration(m, atom, rg1, rg2);
  } else {  // keep the compiler quiet.
    rg1 = 0;
    rg2 = 0;
  }

  // Write the descriptors

  const Accumulator<coord_t>& inter_atomic_distance = stmp.inter_atomic_distance();
  const Accumulator<coord_t>& inter_heteroatom_distance =
      stmp.inter_heteroatom_distance();

  output.add(inter_atomic_distance.maxval());
  output.add(inter_atomic_distance.average());
  output.add(inter_heteroatom_distance.maxval());

  if (0 == inter_heteroatom_distance.n()) {
    output.add(0.0f);
  } else {
    output.add(inter_heteroatom_distance.average_if_available_minval_if_not());
  }

  output.add(a1);  // shadow1
  output.add(a2);  // shadow2
  output.add(a3);  // shadow3

  float abig, amid, amin;
  if (a1 >= a2 && a2 >= a3) {
    abig = a1;
    amid = a2;
    amin = a3;
  } else if (a1 >= a3 && a3 >= a2) {
    abig = a1;
    amid = a3;
    amin = a2;
  } else if (a2 >= a1 && a1 >= a3) {
    abig = a2;
    amid = a1;
    amin = a3;
  } else if (a2 >= a3 && a3 >= a1) {
    abig = a2;
    amid = a3;
    amin = a1;
  } else if (a3 >= a2 && a2 >= a1) {
    abig = a3;
    amid = a2;
    amin = a1;
  } else if (a3 >= a1 && a1 >= a2) {
    abig = a3;
    amid = a1;
    amin = a2;
  } else {  // Will not happen, added to keep the compiler quiet.
    abig = 0;
    amid = 0;
    amin = 0;
  }

  if (write_ordered_areas) {
    output.add(abig);
    output.add(amid);
    output.add(amin);
  }

  assert(abig >= amid && amid >= amin);

  output.add(a1 / a2);                // cigar1
  output.add(a1 / a3);                // cigar13
  output.add(a2 / a3);                // cigar23
  output.add(box_y);                  // boxy
  output.add(box_z);                  // boxz
  output.add(box_x * box_y * box_z);  // boxvolx

  // To get the M* numbers, think of a unit sphere with X being the axis
  // of abig, y==amid and z==amin The perfect sphere is at (1,1,1), the
  // perfect flat shape at (1,0,0) and the perfect tube at (1,1,0) The
  // M* numbers are mostly the distances from these ideal locations

  output.add(abig - amin);           // sphere
  output.add((abig - amin) / abig);  // nsphere
  output.add(1.73205080756 - sqrt((1.0 - amid / abig) * (1.0 - amid / abig) +
                                  (1.0 - amin / abig) * (1.0 - amin / abig)));

  output.add(abig - amid);           // flat
  output.add((abig - amid) / abig);  // nflat
  output.add(1.73205080756 - sqrt((amid * amid) / (abig * abig) +
                                  (amin * amin) / (abig * abig)));  // mflat

  output.add(amid - amin);           // tube
  output.add((amid - amin) / abig);  // ntube
  output.add(1.73205080756 - sqrt((1.0 - amid / abig) * (1.0 - amid / abig) +
                                  (amin * amin) / (abig * abig)));  // mtube

  output.add(static_cast<float>(zz.average()));  // sh_bow, aka pbf

  output.add(static_cast<float>(
      fabs(zleft.average() -
           zright.average())));  // just abs diff between average Z values, sh_zanisotropy

  float zl = abs_zleft.average();
  float zr = abs_zright.average();

  if (fabs(zl) < 0.1 && fabs(zr) < 0.1) {  // sh_abszanisotropy"
    output.add(1.0f);
  } else if (zl < zr) {
    output.add(zl / zr);
  } else {
    output.add(zr / zl);
  }

  if (compute_radius_of_gyration) {
    output.add(rg1);
    output.add(rg2);
  }

  return 1;
}

static float shd_minval[] = {
    8.795,    6.187,       0,           0,          20.59,    41.03,  31.67,
    0.2741,   0.3433,      0.8302,      2.49,       0.006141, 0.1922, 9.83,
    0.145467, 0.930676,    0.280006,    0.0035504,  0.450689, 1.41,   0.0170617,
    0.873726, 0.000817891, 0.000267484, 0.00477638, 2.09,     2.114};

static float shd_maxval[] = {
    26.85,   12.42, 23.5,     12.97,   97.65,  157.4,    136,    0.8886, 0.9871,
    1.83,    13.32, 8.985,    1703,    91.15,  0.726361, 1.5719, 52.67,  0.453261,
    1.05822, 68.77, 0.593547, 1.38787, 1.9922, 2.22546,  1,      7.385,  7.393};

static int
do_fingerprint_output(const resizable_array<float>& shd,
                      IWString_and_File_Descriptor& output) {
  Sparse_Fingerprint_Creator sfc;

  const int n = shd.number_elements();

  for (auto i = 0; i < n; ++i) {
    float v = shd[i];
    if (v < shd_minval[i]) {
      sfc.hit_bit(i, 1);
    } else if (v > shd_maxval[i]) {
      sfc.hit_bit(i, max_bit_count);
    } else {
      const int c = static_cast<int>((v - shd_minval[i]) /
                                         (shd_maxval[i] - shd_minval[i]) * max_bit_count +
                                     0.4999) +
                    1;
      sfc.hit_bit(i, c);
    }
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, output);

  output << tmp << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
do_output(const resizable_array<float>& shd, IWString_and_File_Descriptor& output) {
  if (pbf_only) {
    if (rdkit_only) {
      return 1;
    }

    output.append_number(shd[1], OUTPUT_PRECISION);

    if (check_rdkit) {
      output << ' ';
      output.append_number(shd[2], OUTPUT_PRECISION);
      output << ' ';
      output.append_number(shd[3], OUTPUT_PRECISION);
    }

    return 1;
  }

  // normal output

  const auto n = shd.number_elements();

  for (auto i = 1; i < n; ++i) {
    output << ' ';
    output.append_number(shd[i], OUTPUT_PRECISION);
  }

  return 1;
}

static int
do_output(const IWString& mname, const resizable_array<float>& shd,
          IWString_and_File_Descriptor& output) {
  append_first_token_of_name(mname, output);
  output << ' ';

  output << shd[0];  // always gets written

  do_output(shd, output);

  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;

  if (pbf_only) {
    output << '\n';
  }

  return 1;
}

static int
tshadow(Molecule& m, resizable_array<float>& output) {
  const auto matoms = m.natoms();

  if (matoms < 2) {  // too hard otherwise
    return 1;
  }

  if (pbf_only) {
    if (rdkit_only) {
      const double pbf = PBFRD(m);
      output.add(pbf);
      return 1;
    }

    do_align_by_principal_components(m);

    double pbf = 0.0;

    for (auto i = 0; i < matoms; ++i) {
      pbf += fabs(m.z(i));
    }

    pbf /= static_cast<double>(matoms);

    output.add(pbf);

    if (check_rdkit) {
      const double rdkit_pbf = PBFRD(m);
      rdkit_diff.extra(abs(pbf - rdkit_pbf));
      output.add(static_cast<float>(rdkit_pbf));
      output.add(static_cast<float>(pbf - rdkit_pbf));
    }

    return 1;
  }

  Shadow_Tmp stmp(matoms);

  const Atom** atom = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atom(atom);
  m.atoms(atom);

  return tshadow(m, atom, stmp, output);
}

static void
preprocess(Molecule& m) {
  if (remove_hydrogens) {
    m.remove_all(1);
  }

  return;
}

static int
tshadow(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  int three_dimensional_molecule_encountered = 0;

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (three_dimensional_molecule_encountered) {  // assume everything OK
      ;
    } else if (3 == m->highest_coordinate_dimensionality()) {  // great
      three_dimensional_molecule_encountered = 1;
    } else if (1 == m->highest_coordinate_dimensionality())  // horrible
    {
      cerr << "Warning, no geometry '" << m->name() << "'\n";
      return 0;
    } else {  // benzene might be the first molecule in the dataset...
      cerr << "Warning, flat molecule '" << m->name() << "'\n";
    }

    preprocess(*m);

    resizable_array<float> shd(30);

    if (!tshadow(*m, shd)) {
      return 0;
    }

    do_output(m->name(), shd, output);
  }

  return 1;
}

static int
tshadow(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(input_type != FILE_TYPE_INVALID);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  return tshadow(input, output);
}

static int
tshadow_filter_record(const const_IWSubstring& buffer,
                      IWString_and_File_Descriptor& output) {
  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "Cannot interpret smiles '" << buffer << "'\n";
    return 0;
  }

  preprocess(m);  // we don't check for gettign a 3d structure

  resizable_array<float> shd(30);

  if (!tshadow(m, shd)) {
    return 0;
  }

  return do_fingerprint_output(shd, output);
}

static int
tshadow_filter(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    buffer.chop();
    buffer.remove_leading_chars(smiles_tag.length());

    if (!tshadow_filter_record(buffer, output)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
tshadow_filter(const char* fname,  // most likely '-'
               IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "tshadow_filter:cannot open '" << fname << "'\n";
    return 0;
  }

  return tshadow_filter(input, output);
}

static int
tshadow(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:E:r:V:p:LYGuq:s:e:O:m:dJ:P:btky:h");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot parse -A option\n";
    usage(5);
  }

  if (!process_elements(cl)) {
    cerr << "Cannot process element specifications (-E)\n";
    usage(2);
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', grid_resolution) || grid_resolution <= 0.0) {
      cerr << "The grid resolution (-r) option must be positive\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Grid resolution set to " << grid_resolution << endl;
    }
  }

  if (cl.option_present('V')) {
    if (!set_default_van_der_waals_radius_type(cl, 'V', vdw_type, verbose)) {
      cerr << "Cannot set default VDW radius type (-V option)\n";
      return 7;
    }
  }

  if (cl.option_present('b')) {
    pbf_only = 1;
    position_by_extremeties = 0;

    if (verbose) {
      cerr << "Will only produce PBF descriptor, position by average turned on\n";
    }

    if (cl.option_present('k')) {
      rdkit_only = 1;

      if (verbose) {
        cerr << "Will use the RDKIT algorithm\n";
      }
    }
  }

  if (cl.option_present('t')) {
    check_rdkit = 1;

    if (verbose) {
      cerr << "Will check results vs RDKIT implementation\n";
    }
  }

  if (cl.option_present('h')) {
    remove_hydrogens = 0;

    if (verbose) {
      cerr << "Will NOT remove Hydrogen atoms\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if ("ext" == p) {
      position_by_extremeties = 1;
      if (verbose) {
        cerr << "Alignment done by principal components, initial position by "
                "extremeties\n";
      }
    } else if ("ave" == p) {
      position_by_extremeties = 0;
      if (verbose) {
        cerr << "Alignment done by principal components, initial position averaged "
                "across origin\n";
      }
    } else {
      cerr << "Unrecognised -P qualifier '" << p << "'\n";
      usage(1);
    }

    align_by_principal_components = 1;
  }

  if (cl.option_present('p')) {
    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if (p.starts_with("def")) {
        on = '*';
        off = '.';
      } else if (p.starts_with("on=") && 4 == p.length()) {
        on = p[3];
      } else if (p.starts_with("off=") && 5 == p.length()) {
        off = p[4];
      } else if ('Y' == p) {
        print_grid = 2;
      } else {
        cerr << "Grid printing directives must be 'def', 'on=x' or 'off=y' only\n";
        cerr << "Directive '" << p << "' is invalid\n";
        usage(31);
      }
    }

    if (0 == print_grid) {
      print_grid = 1;
    }

    if (verbose) {
      cerr << "Will print the grid for each molecule, on = '" << on << "', off = '" << off
           << "'\n";
    }
  }

  if (cl.option_present('L')) {
    rotate_to_longest_distance = 1;

    if (verbose) {
      cerr << "Will rotate to longest inter-atom distance\n";
    }
  }

  if (cl.option_present('Y'))  // now the default
  {
    //  perpendicular_too = 1;

    if (verbose) {
      cerr << "Will rotate about Y axis for second computation\n";
    }
  }

  if (cl.option_present('G')) {
    compute_radius_of_gyration = 1;

    if (verbose) {
      cerr << "Spatial and mass weighted radii of gyration also computed\n";
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose > 1, 'q')) {
      cerr << "Cannot read queries\n";
      return 4;
    }

    if (verbose) {
      cerr << "Read " << queries.number_elements() << " queries\n";
    }
  }

  if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring s;
    while (cl.value('s', s, i++)) {
      Substructure_Hit_Statistics* tmp = new Substructure_Hit_Statistics;

      if (!tmp->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 6;
      }

      queries.add(tmp);
    }

    rotate_to_longest_distance = 1;
  }

  if (queries.number_elements() && cl.option_present('e')) {
    cerr << "Sorry, the -q/-s and -e options are mutually exclusive\n";
    usage(8);
  }

  if (!cl.option_present('e')) {
    ;
  } else if (2 != cl.option_count('e')) {
    cerr << "There must be exactly two -e options\n";
    usage(6);
  } else {
    const_IWSubstring e;
    cl.value('e', e, 0);

    if (!process_cmdline_token('e', e, end1queries, verbose > 1)) {
      cerr << "Cannot create end1 query '" << e << "'\n";
      return 6;
    }

    cl.value('e', e, 1);

    if (!process_cmdline_token('e', e, end2queries, verbose > 1)) {
      cerr << "Cannot create end2 query '" << e << "'\n";
      return 6;
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', tag);

    if (!tag.ends_with('<')) {
      tag << '<';
    }

    if (verbose) {
      cerr << "Will produce fingerprints with tag '" << tag << "'\n";
    }

    if (cl.option_present('y')) {
      if (!cl.value('y', max_bit_count) || max_bit_count < 1) {
        cerr << "The max bit count option (-y) must be a whole +ve number\n";
        usage(2);
      }

      if (verbose) {
        cerr << "Bits formed with a max count of " << max_bit_count << endl;
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (tag.length()) {
    ;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('O')) {
    int i = 0;
    IWString fname;
    const_IWSubstring o;

    while (cl.value('O', o, i++)) {
      if ('1' == o) {
        write_oriented_molecule[1] = 1;
      } else if ('2' == o) {
        write_oriented_molecule[2] = 1;
      } else if ('3' == o) {
        write_oriented_molecule[3] = 1;
      } else if (0 == fname.length()) {
        fname = o;
      } else {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        return 4;
      }
    }

    if (0 == fname.length()) {
      cerr << "Must specify file name in the -O options\n";
      usage(3);
    }

    //  If nothing specified for what to write, write everything

    if (0 == count_non_zero_occurrences_in_array(write_oriented_molecule, 4)) {
      set_vector(write_oriented_molecule, 4, 1);
    }

    set_write_isis_standard(1);
    set_write_mdl_charges_as_m_chg(1);

    stream_for_oriented_molecules = new Molecule_Output_Object[4];

    for (int i = 1; i <= 3; i++) {
      if (!write_oriented_molecule[i]) {
        continue;
      }

      stream_for_oriented_molecules[i].add_output_type(FILE_TYPE_SDF);

      IWString tmp;
      tmp << fname << i;

      if (!stream_for_oriented_molecules[i].new_stem(tmp)) {
        cerr << "Cannot initialise stem for oriented molecules '" << tmp << "'\n";
        return 3;
      }
    }

    if (verbose) {
      cerr << "Oriented molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('u')) {
    ok_if_no_vdw_radii = 1;

    if (verbose) {
      cerr << "Will skip molecules for which vdw radii assignment fails\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', multiply_vdw_radii) || multiply_vdw_radii <= 0.0) {
      cerr << "The multiply vdw radii factor option (-m) must be a non-negative value\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will multiply all vdw radii by " << multiply_vdw_radii << endl;
    }
  }

  if (cl.option_present('d')) {
    write_ordered_areas = 1;

    if (verbose) {
      cerr << "Will also write ordered shadow areas\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (tag.length() > 0) {
    if (!tshadow_filter(cl[0], output)) {
      rc = 1;
    }
  } else {
    if (!write_header(output)) {
      cerr << "Cannot write header\n";
      return 8;
    }

    for (int i = 0; i < cl.number_elements(); i++) {
      if (!tshadow(cl[i], input_type, output)) {
        cerr << "Error processing '" << cl[i] << "'\n";
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    if (acc.n() > 0) {
      cerr << "Shadow areas between " << acc.minval() << " and " << acc.maxval()
           << " average " << acc.average() << endl;
    }
  }

  if (check_rdkit) {
    cerr << "Performed " << rdkit_diff.n() << " comparisons with RDKIT implementation\n";
    if (rdkit_diff.n() > 0) {
      cerr << " differences between " << static_cast<float>(rdkit_diff.minval())
           << " and " << static_cast<float>(rdkit_diff.maxval());
      if (rdkit_diff.n() > 1) {
        cerr << " ave " << static_cast<float>(rdkit_diff.average());
      }
      cerr << endl;
    }
  }

  if (nullptr != stream_for_oriented_molecules) {
    delete[] stream_for_oriented_molecules;
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = tshadow(argc, argv);

  return rc;
}
