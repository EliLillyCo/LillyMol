/*
  Aligns molecules according to a substructure query
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"

#include "Molecule_Tools/spatially_common_matched_atoms.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static IWString stem_for_output;

static int max_matches_to_process = 0;

static int ignore_queries_not_hitting = 0;

static int write_molecules_not_matching_queries = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int embedding_to_process = -1;

static int matched_atoms_in_quadrant = 0;

/*
  When doing quadrant determinations, we can center the molecule based on
  either the extremeties of the coordinates, or the average coordinate
*/

static int use_spatial_extremeties_for_quadrant_determination = 1;

/*
  Matched atoms ORIGIN_ATOM, X_ATOM and Y_ATOM can be used to orient
  the remaining molecules
*/

static resizable_array<int> origin_atom;

static resizable_array<int> x_atom;

static resizable_array<int> y_atom;

/*
  Dec 99. Enhance to allow the coordinates of the "origin" and "axes" to be
  input on the command line
*/

static Space_Vector<coord_t> origin(0.0, 0.0, 0.0);

static Space_Vector<coord_t> xaxis(1.0, 0.0, 0.0);

static Space_Vector<coord_t> yaxis(0.0, 1.0, 0.0);

static extending_resizable_array<int> hits_found;

static Molecule_Output_Object stream_for_non_matches;

static int molecules_read = 0;
static int molecules_written = 0;

static int apply_isotopic_labels = 0;

static int align_so_first_matched_atom_is = 0;

static angle_t final_x_rotation = static_cast<angle_t>(0.0);
static angle_t final_y_rotation = static_cast<angle_t>(0.0);
static angle_t final_z_rotation = static_cast<angle_t>(0.0);

static int discern_embedding_via_spatial_location = 0;

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
  cerr << "Aligns molecules according to a substructure query\n";
  cerr << "Usage: " << prog_name << " <options> <input_file>...\n";
  cerr << "  -q <query>     specify substructure query\n";
  cerr << "  -s <smarts>    specify smarts for search\n";
  cerr << "  -O <number>    matched atom <number> is the origin\n";
  cerr << "  -X <number>    rotate so matched atom <number> is on the X axis\n";
  cerr << "  -Y <number>    rotate so matched atom <number> is on the Y axis\n";
//cerr << "  -r <coord>     coordinates for the oRigin (xxx.xxx,yyy.yyy,zzz.zzz)\n";
//cerr << "  -x <coord>     coordinates for the X axis (xxx.xxx,yyy.yyy,zzz.zzz)\n";
//cerr << "  -y <coord>     coordinates for the Y axis (xxx.xxx,yyy.yyy,zzz.zzz)\n";
//cerr << "  -c             align to centroid of matched atoms\n";
  cerr << "  -e <number>    process embedding number <number> (starts with 0)\n";
  cerr << "  -z i           ignore molecules not matching any query\n";
  cerr << "  -z w           write molecules not matching any query\n";
  cerr << "  -z f           take the first of multiple query matches\n";
  cerr << "  -w .           (3D) when multiple embeddings, use the one closest to most common location\n";
  cerr << "  -n <file>      write non-matching molecules to <file>\n";
  cerr << "  -Q <quad>      ensure the matched atoms are in a given quadrant\n";
  cerr << "  -Q geom        centre molecule at geometric centre\n";
  cerr << "  -m             remove atoms not matched by the query\n";
  cerr << "  -u             perceive unique embeddings only\n";
  cerr << "  -k             don't perceive symmetrically equivalent matches\n";
  cerr << "  -f <r,l>       align so first matched atom is rightmost\n";
  cerr << "  -F x=,y=,z=    apply one or more final rotations about specified axes\n";
  cerr << "  -p <iso>       apply isotopic labels to matched atoms\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -E autocreate  automatically create new elements when encountered\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
shift_to_average_coordinate(Molecule& m)
{
  const Atom* a = m.atomi(0);

  coord_t xsum = a->x();
  coord_t ysum = a->y();

  int matoms = m.natoms();

  for (int i = 1; i < matoms; i++) {
    a = m.atomi(i);

    xsum += a->x();
    ysum += a->y();
  }

  xsum = xsum / static_cast<coord_t>(matoms);
  ysum = ysum / static_cast<coord_t>(matoms);

  m.translate_atoms(-xsum, -ysum, 0.0);

  return 1;
}

static int
do_matched_atoms_in_quadrant(Molecule& m, const Set_of_Atoms& e, int dimensionality,
                             Molecule_Output_Object& output)
{
  if (3 == dimensionality) {
    cerr << "Sorry, quadrant processing set up for 2D only\n";
    return 0;
  }

  coord_t xmin, xmax, ymin, ymax;

  if (use_spatial_extremeties_for_quadrant_determination) {
    m.spatial_extremeties(xmin, xmax, ymin, ymax);
    coord_t dx = (xmin + xmax) * 0.5;
    coord_t dy = (ymin + ymax) * 0.5;

    m.translate_atoms(-dx, -dy, 0.0);
  } else {
    shift_to_average_coordinate(m);
  }

  // Now work out the average of the coordinates of the matched atoms
  // use xmin and ymin as accumulators

  xmin = static_cast<coord_t>(0.0);
  ymin = static_cast<coord_t>(0.0);

  int n = e.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = e[i];
    if (j < 0) {
      continue;
    }

    const Atom* a = m.atomi(j);

    xmin += a->x();
    ymin += a->y();
  }

  xmin = xmin / static_cast<coord_t>(n);
  ymin = ymin / static_cast<coord_t>(n);

  if (xmin < 0.0 && (1 == matched_atoms_in_quadrant || 4 == matched_atoms_in_quadrant)) {
    m.rotate_atoms(yaxis, M_PI);
  }

  if (ymin < 0.0 && (1 == matched_atoms_in_quadrant || 2 == matched_atoms_in_quadrant)) {
    m.rotate_atoms(xaxis, M_PI);
  }

  return output.write(m);
}

static int
align_molecule_single_atom(Molecule& m, atom_number_t zatom,
                           Molecule_Output_Object& output)
{
  Space_Vector<coord_t> o = *(m.atomi(zatom));

  m.translate_atoms(-o.x(), -o.y(), -o.z());

  int matoms = m.natoms();

  if (1 == matoms) {
    return 1;
  }

  float max_dist_from_origin = 0.0;
  atom_number_t furthest_atom = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; i++) {
    if (i == zatom) {
      continue;
    }

    Space_Vector<coord_t> ci(*(m.atomi(i)));

    float d = ci.norm();  // no need to take square root

    //  cerr << "Atom " << i << " at distance " << d << endl;
    if (d > max_dist_from_origin) {
      max_dist_from_origin = d;
      furthest_atom = i;
    }
  }

  if (INVALID_ATOM_NUMBER == furthest_atom) {
    cerr << "Huh, no furthest atom in '" << m.name() << "'\n";
    return 0;
  }

  const Atom* f = m.atomi(furthest_atom);

  Coordinates fx(f->x(), f->y(), f->z());

  fx.normalise();

  Coordinates negative_x(-1.0, 0.0, 0.0);

  angle_t theta = negative_x.angle_between_unit_vectors(fx);

  fx.cross_product(negative_x);

  fx.normalise();

  m.rotate_atoms(fx, theta);

  return output.write(m);
}

static int
handle_non_matching_molecules(Molecule& m, const int max_query_atoms_matched,
                              Molecule_Output_Object& output)
{
  if (verbose > 1) {
    cerr << "Zero hits to query, " << max_query_atoms_matched << " query atoms matched\n";
  }

  if (stream_for_non_matches.active()) {
    stream_for_non_matches.write(m);
  }

  if (!ignore_queries_not_hitting) {
    if (0 == verbose) {
      cerr << "No hits to query, molecule '" << m.name() << "', only matched "
           << max_query_atoms_matched << " query atoms\n";
    }
    return 0;
  }

  if (write_molecules_not_matching_queries) {
    molecules_written++;
    return output.write(m);
  }

  return 1;
}

static int
fetch_coordinates(const Molecule& m, const Set_of_Atoms& e,
                  const resizable_array<int>& ndx, Space_Vector<coord_t>& xyz)
{
  for (int i = 0; i < ndx.number_elements(); i++) {
    int j = ndx[i];

    if (!e.ok_index(j)) {
      cerr << "invalid matched atom number " << j << " query produced "
           << e.number_elements() << " items\n";
      return 0;
    }

    const Atom* a = m.atomi(e[j]);

    xyz += *a;
  }

  if (ndx.number_elements() > 1) {
    xyz /= static_cast<coord_t>(ndx.number_elements());
  }

  return 1;
}

// #define DEBUG_ALIGN_MOLECULE

static int
align_molecule(Molecule& m, const Set_of_Atoms& e)
{
  if (0 == x_atom.size()) {
    Space_Vector<coord_t> y;

    if (!fetch_coordinates(m, e, y_atom, y)) {
      return 0;
    }

    y.normalise();  // to a unit vector

    Coordinates yy = yaxis - origin;
    yy.normalise();

    coord_t theta = yy.angle_between_unit_vectors(y);

#ifdef DEBUG_ALIGN_MOLECULE
    cerr << "First rotation " << (theta * RAD2DEG) << endl;
#endif

    yy.cross_product(y);
    yy.normalise();  // probably not necessary

    m.rotate_atoms(yy, -theta);

    return 1;
  }

  if (x_atom.size() > 0) {
    Space_Vector<coord_t> x;

    if (!fetch_coordinates(m, e, x_atom, x)) {
      return 0;
    }

    x.normalise();  // to a unit vector

    Coordinates xx = xaxis - origin;
    xx.normalise();

    coord_t theta = xx.angle_between_unit_vectors(x);

#ifdef DEBUG_ALIGN_MOLECULE
    cerr << "First rotation " << (theta * RAD2DEG) << endl;
#endif

    xx.cross_product(x);
    xx.normalise();  // probably not necessary

    m.rotate_atoms(xx, -theta);
  }

  if (0 == y_atom.size()) {
    return 1;
  }

  Space_Vector<coord_t> y;

  if (!fetch_coordinates(m, e, y_atom, y)) {
    return 0;
  }

  // Note that we almost certainly won't be able to get the Y atom on the Y axis,
  // so the rotation will be around the X axis to get it in the X/Y plane

#ifdef DEBUG_ALIGN_MOLECULE
  cerr << "Y axis atom at " << y << endl;
#endif

  y -= origin;

  // In order to get the atom into the X,Y plane, we need to set to X coordinate to 0.0

  y.x() = 0.0;

  // what if the atom is positioned along the X axis

  if (y.norm() < static_cast<coord_t>(1.0e-04)) {
    return 1;
  }

  y.normalise();  // to a unit vector

  Coordinates yy = yaxis - origin;
  yy.normalise();

  angle_t theta = yy.angle_between_unit_vectors(y);

#ifdef DEBUG_ALIGN_MOLECULE
  cerr << "Between " << y << " and Y axis is " << theta << endl;
#endif

  Coordinates rotation_axis = xaxis - origin;  // we will rotate about here
  rotation_axis.normalise();

  if (y.z() > 0) {
    theta = -theta;
  }

  m.rotate_atoms(rotation_axis, theta + M_PI);

  // By convention, y_atom should have a positive Y coordinate

  atom_number_t first_y_atom = e[y_atom[0]];

  y = (*(m.atomi(first_y_atom))) - origin;

  if (y.y() < 0.0) {
    int matoms = m.natoms();
    for (int i = 0; i < matoms; i++) {
      coord_t y = m.y(i);
      m.sety(i, -y);
    }
  }

  return 1;
}

static void
set_all_z_coordinates_to_zero(Molecule& m)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    m.setz(i, static_cast<coord_t>(0.0));
  }

  return;
}

/*
  Translate the origin atom to the real (0.0, 0.0, 0.0) origin. Call another
  function to do the rotations, then move the molecule back to the actual origin
*/

static int
align_molecule_(Molecule& m, const Set_of_Atoms& e, int dimensionality,
                Molecule_Output_Object& output)
{
  assert(origin_atom.size() > 0);

  Space_Vector<coord_t> o;
  if (!fetch_coordinates(m, e, origin_atom, o)) {
    return 0;
  }

  m.translate_atoms(-(o.x()), -(o.y()), -(o.z()));

  int rc = align_molecule(m, e);
  if (0 == rc) {
    return 0;
  }

  if (0.0 != origin.x() || 0.0 != origin.y() || 0.0 != origin.z()) {
    m.translate_atoms(origin.x(), origin.y(), origin.z());
  }

  molecules_written++;

  if (2 == dimensionality) {
    set_all_z_coordinates_to_zero(m);
  }

  if (static_cast<angle_t>(0.0) != final_x_rotation) {
    Space_Vector<coord_t> x(1.0, 0.0, 0.0);
    m.rotate_atoms(x, final_x_rotation);
  }

  if (static_cast<angle_t>(0.0) != final_y_rotation) {
    Space_Vector<coord_t> y(0.0, 1.0, 0.0);
    m.rotate_atoms(y, final_y_rotation);
  }

  if (static_cast<angle_t>(0.0) != final_z_rotation) {
    Space_Vector<coord_t> z(0.0, 0.0, 1.0);
    m.rotate_atoms(z, final_z_rotation);
  }

  return output.write(m);
}

static int
do_apply_isotopic_labels(Molecule& m, const Set_of_Atoms& e, int iso)
{
  int n = e.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = e[i];

    m.set_isotope(j, iso);
  }

  return 1;
}

static int
align_molecule(Molecule& m, const Set_of_Atoms& e, int dimensionality,
               Molecule_Output_Object& output)
{
  if (apply_isotopic_labels) {
    do_apply_isotopic_labels(m, e, apply_isotopic_labels);
  }

  if (align_so_first_matched_atom_is) {
    return align_molecule_single_atom(m, e[0], output);
  }

  if (matched_atoms_in_quadrant > 0) {
    return do_matched_atoms_in_quadrant(m, e, dimensionality, output);
  }

  int rc = align_molecule_(m, e, dimensionality, output);

  if (apply_isotopic_labels) {
    do_apply_isotopic_labels(m, e, 0);
  }

  return rc;
}

static int
align_molecule(Molecule& m, Molecule_Output_Object& output)
{
  const int dimensionality = m.highest_coordinate_dimensionality();

  const int nq = queries.number_elements();

  int matches_found = 0;
  int max_query_atoms_matched = 0;

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(m, sresults);

    hits_found[nhits]++;

    if (verbose > 1) {
      cerr << "nhits = " << nhits << " '" << m.name() << "'\n";
    }

    if (0 == nhits) {
      if (max_query_atoms_matched < sresults.max_query_atoms_matched_in_search()) {
        max_query_atoms_matched = sresults.max_query_atoms_matched_in_search();
      }
      continue;
    }

    if (max_matches_to_process > 0 && nhits > max_matches_to_process) {
      nhits = max_matches_to_process;
    }

    for (int i = 0; i < nhits; i++) {
      if (embedding_to_process >= 0 && i != embedding_to_process) {
        continue;
      }

      const Set_of_Atoms* e = sresults.embedding(i);

      if (!align_molecule(m, *e, dimensionality, output)) {
        return 0;
      }
    }

    matches_found = 1;
    break;
  }

  if (0 == matches_found) {
    return handle_non_matching_molecules(m, max_matches_to_process, output);
  }

  return 1;
}

static int
align_molecule(resizable_array_p<Molecule>& molecules,
               const Substructure_Results* sresults, Molecule_Output_Object& output)
{
  const int n = molecules.number_elements();

  for (int i = 0; i < n; ++i) {
    Molecule& mi = *molecules[i];

    const int dimensionality = mi.highest_coordinate_dimensionality();

    if (!align_molecule(mi, *sresults[i].embedding(0), dimensionality, output)) {
      return 0;
    }
  }

  return 1;
}

static int
do_discern_embedding_via_spatial_location(data_source_and_type<Molecule>& input,
                                          Molecule_Output_Object& output)
{
  const int n = input.molecules_remaining();

  resizable_array_p<Molecule> molecules;
  molecules.resize(n);

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules.add(m);
  }

  assert(n == molecules.number_elements());

  molecules_read = n;

  Substructure_Results* sresults = new Substructure_Results[n];
  std::unique_ptr<Substructure_Results[]> free_sresults(sresults);

  for (int i = 0; i < n; ++i) {
    Molecule& mi = *molecules[i];
    const int nhits = queries[0]->substructure_search(mi, sresults[i]);
    if (0 == nhits) {
      if (ignore_queries_not_hitting) {
        continue;
      }

      if (write_molecules_not_matching_queries) {
        output.write(mi);
      }

      return 0;
    }
  }

  spatially_common_matched_atoms(molecules, sresults);

  return align_molecule(molecules, sresults, output);
}

static int
align_molecule(data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (!align_molecule(*m, output)) {
      return 0;
    }
  }

  return 1;
}

static int
align_molecule(const char* fname, FileType input_type, Molecule_Output_Object& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  static int first_call = 1;

  if (stem_for_output.length() && first_call)  // need to establish the output stream
  {
    if (!output.new_stem(stem_for_output, 1)) {
      cerr << "Cannot open output file(s) with stem '" << stem_for_output << "'\n";
      return 0;
    }

    first_call = 0;
  } else if (stem_for_output.length()) {  // should already be set
    ;
  } else  // need to open output file for this input
  {
    if (!output.new_stem(fname)) {
      cerr << "Cannot open output file(s) for input '" << fname << "'\n";
      return 0;
    }
  }

  if (discern_embedding_via_spatial_location) {
    return do_discern_embedding_via_spatial_location(input, output);
  }

  return align_molecule(input, output);
}

/*
  This can be many forms.

  A simple number

    3

  A set of numbers

    0,4,2

  One number with geometry
    '5=0.2122,3.114'

*/

static int
get_atom_and_coordinates(const const_IWSubstring& token,
                         resizable_array<int>& matched_atoms, Space_Vector<coord_t>& v)
{
  int a;

  if (token.numeric_value(a))  // just a number
  {
    matched_atoms.add(a);
    return 1;
  }

  if (token.contains('=')) {
    const_IWSubstring zatom = token.before('=');
    const_IWSubstring coords = token.after('=');

    if (0 == zatom.length() || 0 == coords.length()) {
      cerr << "Invalid atom and coordinate specification '" << token << "'\n";
      return 0;
    }

    if (!zatom.numeric_value(a) || a < 0) {
      cerr << "Invalid atom specifier '" << a << "' of '" << token << "'\n";
      return 0;
    }

    coord_t x, y, z;

    int token_number = 0;  // which token
    int j = 0;             // character position
    const_IWSubstring c;
    while (coords.nextword(c, j, ',')) {
      coord_t xyz;
      if (!c.numeric_value(xyz)) {
        cerr << "Invalid " << token_number << " coordinate specification '" << c
             << "' from '" << token << "'\n";
        return 0;
      }

      if (0 == token_number) {
        x = xyz;
      } else if (1 == token_number) {
        y = xyz;
      } else {
        z = xyz;
      }

      token_number++;
    }

    if (3 != token_number) {
      cerr << "Must specify 3 coordinates '" << token << "'\n";
      return 0;
    }

    v.setxyz(x, y, z);
    matched_atoms.add(a);

    return 1;
  }

  //  List of atoms

  int j = 0;
  const_IWSubstring anum;
  while (token.nextword(anum, j, ',')) {
    int k;
    if (!anum.numeric_value(k) || k < 0) {
      cerr << "Invalid matched atom number '" << anum << "'\n";
      return 0;
    }

    matched_atoms.add(k);
  }

  return matched_atoms.size();
}

static int
read_query_from_molecule(Substructure_Hit_Statistics& query, const_IWSubstring& fname)
{
  FileType input_type = discern_file_type_from_name(fname);
  if (input_type == FILE_TYPE_INVALID) {
    cerr << "Cannot discern file type '" << fname << "'\n";
    return 0;
  }

  data_source_and_type<MDL_Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  MDL_Molecule* m = input.next_molecule();
  if (nullptr == m) {
    cerr << "Cannot read structure from '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<MDL_Molecule> free_m(m);

  Molecule_to_Query_Specifications mqs;

  if (!query.create_from_molecule(*m, mqs)) {
    cerr << "Cannot create query from molecule in '" << fname << "'\n";
    return 0;
  }

  return 1;
}

static int
report_which_matched_atoms(std::ostream& os, const resizable_array<int>& matched_atoms,
                           const char* s)
{
  os << "Matched atom";
  if (matched_atoms.size() > 1) {
    os << 's';
  }

  for (int i = 0; i < matched_atoms.number_elements(); i++) {
    os << ' ' << matched_atoms[i];
  }

  os << " will define " << s << '\n';

  return 1;
}

static int
align_molecule(int argc, char** argv)
{
  Command_Line cl(argc, argv, "S:i:o:q:e:A:E:vO:X:Y:kus:z:n:r:x:y:p:f:2F:Q:w:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify one or more queries via the -q option\n";
    usage(8);
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose > 1, 'q')) {
      cerr << "Cannot process query specification(s) -q\n";
      return 3;
    }
  } else if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring s;
    while (cl.value('s', s, i++)) {
      Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;

      if (!q->create_from_smarts(s)) {
        cerr << "INvalid smarts '" << s << "'\n";
        return 12;
      }
      queries.add(q);
    }
  }

  int nq = queries.number_elements();

  if (0 == nq) {
    cerr << "No queries specified, cannot continue\n";
    usage(2);
  }

  if (cl.option_present('u')) {
    for (int i = 0; i < nq; i++) {
      queries[i]->set_find_unique_embeddings_only(1);
    }

    if (verbose) {
      cerr << "Only unique embeddings will be found\n";
    }
  }

  if (cl.option_present('k')) {
    for (int i = 0; i < nq; i++) {
      queries[i]->set_do_not_perceive_symmetry_equivalent_matches(1);
    }

    if (verbose) {
      cerr << "Symmetrically equivalent matches will not be found\n";
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', embedding_to_process) || embedding_to_process < 0) {
      cerr << "The embedding number (-e option) must be non-negative\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Will process embedding " << embedding_to_process << endl;
    }
  }

  if (cl.option_present('f')) {
    const_IWSubstring f = cl.string_value('f');

    if ('l' == f) {
      align_so_first_matched_atom_is = -1;
      if (verbose) {
        cerr << "Will align so first matched atom is leftmost\n";
      }
    } else if ('r' == f) {
      align_so_first_matched_atom_is = +1;
      if (verbose) {
        cerr << "Will align so first matched atom is rightmost\n";
      }
    } else {
      cerr << "Unrecognised -f qualifier '" << f << "'\n";
      usage(3);
    }
  }

  if (cl.option_present('F')) {
    const_IWSubstring f;
    int i = 0;
    while (cl.value('F', f, i++)) {
      const_IWSubstring token;
      int j = 0;
      while (f.nextword(token, j, ',')) {
        const_IWSubstring axis;
        angle_t angle;
        if (!token.split_into_directive_and_value(axis, '=', angle)) {
          cerr << "Invalid final rotation specification (-F) '" << f << "'\n";
          usage(4);
        }

        if ('x' == axis || 'X' == axis) {
          final_x_rotation = angle;
          if (verbose) {
            cerr << "Final X rotation of " << final_x_rotation << endl;
          }
          final_x_rotation = final_x_rotation * DEG2RAD;
        } else if ('y' == axis || 'Y' == axis) {
          final_y_rotation = angle;
          if (verbose) {
            cerr << "Final Y rotation of " << final_y_rotation << endl;
          }
          final_y_rotation = final_y_rotation * DEG2RAD;
        } else if ('z' == axis || 'Z' == axis) {
          final_z_rotation = angle;
          if (verbose) {
            cerr << "Final Z rotation of " << final_z_rotation << endl;
          }
          final_z_rotation = final_z_rotation * DEG2RAD;
        } else {
          cerr << "Unrecognised axis specifiecation '" << axis << "'\n";
          usage(9);
        }
      }
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', stem_for_output);
    if (verbose) {
      cerr << "Will use '" << stem_for_output << "' for output file(s)\n";
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    int i = 0;
    while (cl.value('z', z, i++)) {
      if ('i' == z) {
        ignore_queries_not_hitting = 1;
        if (verbose) {
          cerr << "Molecules not matching the query will be ignored\n";
        }
      } else if ('w' == z) {
        write_molecules_not_matching_queries = 1;
        if (verbose) {
          cerr << "MOlecules not matching the query will be written\n";
        }
      } else if ('f' == z) {
        max_matches_to_process = 1;
        if (verbose) {
          cerr << "Will take the first of multiple hits\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(15);
      }
    }
  }

  if (cl.option_present('Q')) {
    const_IWSubstring q;
    for (int i = 0; cl.value('Q', q, i); i++) {
      if ("geom" == q) {
        use_spatial_extremeties_for_quadrant_determination = 0;
        if (verbose) {
          cerr << "Will centre molecules to geometric centre for quadrant "
                  "determinations\n";
        }
      } else if (q.numeric_value(matched_atoms_in_quadrant) &&
                 matched_atoms_in_quadrant >= 1 && matched_atoms_in_quadrant <= 4) {
        if (verbose) {
          cerr << "Will try to put matched atoms in the " << matched_atoms_in_quadrant
               << " quadrant\n";
        }
      } else {
        cerr << "The matched_atoms_in_quadrant option (-Q) must be a valid quadrant\n";
        usage(3);
      }
    }

    if (0 == matched_atoms_in_quadrant) {
      cerr << "Must specify into which quadrant to place matched atoms (-Q)\n";
      usage(7);
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

  if (!cl.option_present('o')) {
    cerr << "Must specify output format by -o option\n";
    usage(18);
  }

  Molecule_Output_Object output;
  if (!output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    usage(28);
  }

  if (stem_for_output.length()) {  // don't check
    ;
  } else if (output.would_overwrite_input_files(cl, cl[0])) {
    cerr << "Cannot overwrite input\n";
    return 8;
  }

  if (stem_for_output.length()) {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (output.would_use_name(stem_for_output.null_terminated_chars(), cl[i])) {
        cerr << "Sorry, cannot overwrite my input '" << cl[i] << "'\n";
        return 71;
      }
    }
  }

  if (cl.option_present('n')) {
    if (!stream_for_non_matches.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s) for non-match file\n";
      usage(28);
    }

    IWString n = cl.string_value('n');

    if (!stream_for_non_matches.new_stem(n)) {
      cerr << "Cannot open non-match file '" << n << "'\n";
      return 18;
    }

    if (verbose) {
      cerr << "Molecules not matching the query written to '" << n << "'\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', apply_isotopic_labels) || apply_isotopic_labels < 1) {
      cerr << "The isotope to apply must be a non-negative number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will label matched atoms with isotope " << apply_isotopic_labels << endl;
    }
  }

  if (cl.option_present('f')) {
    ;
  } else if (cl.option_present('Q')) {
    ;
  } else if (!cl.option_present('O')) {
    cerr << "MUst specify origin atom via -O option\n";
    usage(27);
  }

  if (cl.option_present('O')) {
    if (!get_atom_and_coordinates(cl.string_value('O'), origin_atom, origin)) {
      cerr << "Invalid origin atom specification (-O)\n";
      usage(91);
    }

    if (verbose) {
      report_which_matched_atoms(cerr, origin_atom, "origin");
    }
  }

  if (cl.option_present('X')) {
    if (!get_atom_and_coordinates(cl.string_value('X'), x_atom, xaxis)) {
      cerr << "Invalid X axis specification (-X)\n";
      usage(91);
    }

    if (verbose) {
      report_which_matched_atoms(cerr, x_atom, "X axis");
    }
  }

  if (cl.option_present('Y')) {
    //  if (0 == x_atom.size())
    //  {
    //    cerr << "Cannot define the Y axis unless the X axis is also defined\n";
    //    usage(38);
    //  }

    if (!get_atom_and_coordinates(cl.string_value('Y'), y_atom, yaxis)) {
      cerr << "Invalid Y axis specification (-Y)\n";
      usage(91);
    }

    if (verbose) {
      report_which_matched_atoms(cerr, y_atom, "Y axis");
    }
  }

  // If we have single atoms, make sure they are all distinct

  if (1 == x_atom.number_elements() && 1 == origin_atom.number_elements() &&
      x_atom[0] == origin_atom[0]) {
    cerr << "Both the origin and X axis atoms are the same. Impossible\n";
    usage(27);
  }

  if (1 == y_atom.number_elements() && 1 == origin_atom.number_elements() &&
      y_atom[0] == origin_atom[0]) {
    cerr << "Both the origin and Y axis atoms are the same. Impossible\n";
    usage(27);
  }

  if (1 == y_atom.number_elements() && 1 == x_atom.number_elements() &&
      y_atom[0] == x_atom[0]) {
    cerr << "Both the Y and Y axis atoms are the same. Impossible\n";
    usage(27);
  }

  if (cl.option_present('r')) {
    const_IWSubstring r = cl.string_value('r');
    if (!origin.read(r, ',')) {
      cerr << "Cannot parse origin atom specification '" << r << "'\n";
      usage(21);
    }

    if (verbose) {
      cerr << "Origin atom at " << origin << endl;
    }
  }

  if (cl.option_present('x')) {
    const_IWSubstring x = cl.string_value('x');
    if (!xaxis.read(x, ',')) {
      cerr << "Cannot parse X axis atom specification '" << x << "'\n";
      usage(21);
    }

    if (verbose) {
      cerr << "X axis atom at " << xaxis << endl;
    }
  }

  if (cl.option_present('y')) {
    const_IWSubstring y = cl.string_value('y');
    if (!yaxis.read(y, ',')) {
      cerr << "Cannot parse Y axis atom specification '" << y << "'\n";
      usage(21);
    }

    if (verbose) {
      cerr << "Y axis atom at " << yaxis << endl;
    }
  }

  if (cl.option_present('w')) {
    if (queries.number_elements() > 1) {
      cerr << "Sorry, cannot use the -y option with multiple queries, (see Ian)\n";
      return 1;
    }

    discern_embedding_via_spatial_location = 1;

    if (verbose) {
      cerr << "When multiple matches, will remove matched atoms closest to most common "
              "location\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!align_molecule(cl[i], input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    for (int i = 0; i < hits_found.number_elements(); i++) {
      if (hits_found[i]) {
        cerr << hits_found[i] << " molecules had " << i << " hits to the query\n";
      }
    }

    cerr << molecules_read << " molecules read, " << molecules_written
         << " molecules written\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = align_molecule(argc, argv);

  return rc;
}
