/*
  Computes geometric principal components of a molecule
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/jama/jama_svd.h"
#include "Foundational/tnt/tnt_array2d_utils.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static Molecule_Output_Object stream_for_oriented_molecules;

static int position_by_extremeties = 0;

static const Element * axis_element = NULL;

typedef double tnt_float_type;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Computes geometric principal components\n";
  cerr << "  -e            do initial positioning by extremeties rather than averages\n";
  cerr << "  -S <fname>    stream for oriented molecules\n";
  cerr << "  -o <type>     output type for -S file\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static void
do_position_by_extremeties (Molecule & m,
                            TNT::Array2D<tnt_float_type> & geometry)
{
  coord_t xmin, xmax, ymin, ymax, zmin, zmax;

  m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);

  coord_t xmid = (xmax + xmin) * 0.5;
  coord_t ymid = (ymax + ymin) * 0.5;
  coord_t zmid = (zmax + zmin) * 0.5;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    Coordinates c;
    m.get_coords(i, c);

    coord_t newx = c.x() - xmid;
    coord_t newy = c.y() - ymid;
    coord_t newz = c.z() - zmid;

    m.setxyz(i, newx, newy, newz);

    geometry[i][0] = newx;
    geometry[i][1] = newy;
    geometry[i][2] = newz;

//  cerr << newx << ' ' << newy << ' ' << newz << '\n';
  }

  return;
}

static void
do_position_by_averages (Molecule & m,
                         TNT::Array2D<tnt_float_type> & geometry)
{
  Accumulator<float> x, y, z;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    Coordinates c;
    m.get_coords(i, c);

    x.extra(c.x());
    y.extra(c.y());
    z.extra(c.z());
  }

  double xave = x.average();
  double yave = y.average();
  double zave = z.average();

  for (int i = 0; i < matoms; i++)
  {
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

/*
  Make sure that the decomposition has worked

*/

#ifdef CHECK_COMPUTATION
static int
check_matrix_product (const Molecule & m,
                      const float * q,
                      const float * d,
                      const float * e,
                      const float * p)
{
  int matoms = m.natoms();

  cerr << matoms << " atoms being checked\n";

  float * tmp1 = new float[3 * 3]; std::unique_ptr<float[]>free_tmp1(tmp1);
  float * tmp2 = new float[3 * 3]; std::unique_ptr<float[]>free_tmp2(tmp2);
  float * tmp3 = new float[matoms * 3]; std::unique_ptr<float[]>free_tmp3(tmp3);

  set_vector(tmp1, 9, static_cast<float>(0.0));
  set_vector(tmp2, 9, static_cast<float>(0.0));
  set_vector(tmp3, matoms * 3, static_cast<float>(0.0));

  tmp1[0] = d[0];
  tmp1[3] = e[0];
  tmp1[4] = d[1];
  tmp1[7] = e[1];
  tmp1[8] = d[2];

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      float t = 0.0;

      for (int k = 0; k < 3; k++)
      {
        t += tmp1[f77(3, i, k)] * p[f77(3, k, j)];
      }

      tmp2[f77(3, i, j)] = t;
    }
  }

  cerr << matoms << " atoms\n";

  for (int i = 0; i < matoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      float t = 0.0;

      for (int k = 0; k < 3; k++)
      {
        t += q[f77(matoms, i, k)] * tmp2[f77(3, k, j)];
      }

      tmp3[f77(matoms, i, j)] = t;
    }
  }

  cerr << "Check transformation " << matoms << " atoms\n";
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    cerr << i << ' ' << a->x() << ' ' << tmp3[i] << ' ' << a->y() << ' ' << tmp3[matoms + i] << ' ' << a->z() << ' ' << tmp3[matoms + matoms + i] << '\n';
  }

  return 1;
}
#endif

/*
  Worked example from 
  http://www.puffinwarellc.com/index.php/news-and-articles/articles/30-singular-value-decomposition-tutorial.html?start=2
*/

#ifdef GOLF_EXAMPLE
static int
do_golf_example ()
{
  TNT::Array2D<tnt_float_type> geometry (9, 3);

  geometry[0][0] = 4.0;
  geometry[0][1] = 4.0;
  geometry[0][2] = 5.0;

  geometry[1][0] = 4.0;
  geometry[1][1] = 5.0;
  geometry[1][2] = 5.0;

  geometry[2][0] = 3.0;
  geometry[2][1] = 3.0;
  geometry[2][2] = 2.0;

  geometry[3][0] = 4.0;
  geometry[3][1] = 5.0;
  geometry[3][2] = 4.0;

  geometry[4][0] = 4.0;
  geometry[4][1] = 4.0;
  geometry[4][2] = 4.0;

  geometry[5][0] = 3.0;
  geometry[5][1] = 5.0;
  geometry[5][2] = 4.0;

  geometry[6][0] = 4.0;
  geometry[6][1] = 4.0;
  geometry[6][2] = 3.0;

  geometry[7][0] = 2.0;
  geometry[7][1] = 4.0;
  geometry[7][2] = 4.0;

  geometry[8][0] = 5.0;
  geometry[8][1] = 5.0;
  geometry[8][2] = 5.0;

  float xmid = 0.0;
  float ymid = 0.0;
  float zmid = 0.0;
  for (int i = 0; i < 9; i++)
  {
    xmid += geometry[i][0];
    ymid += geometry[i][1];
    zmid += geometry[i][2];
  }

  float xmin = xmid / 9.0;
  float ymin = ymid / 9.0;
  float zmin = zmid / 9.0;

  JAMA::SVD<tnt_float_type> svd(geometry);

  TNT::Array1D<tnt_float_type> singular_values;

  svd.getSingularValues(singular_values);

  TNT::Array2D<tnt_float_type> u;
  svd.getU(u);

  cerr << "U matrix\n";
  cerr << u;

  TNT::Array2D<tnt_float_type> s;
  svd.getS(s);

  cerr << "S matrix\n";
  cerr << s;

  TNT::Array2D<tnt_float_type> v;
  svd.getV(v);

// Need to transpose it to get vt

  for (int i = 0; i < 2; i++)
  {
    for (int j = i + 1; j < 3; j++)
    {
      tnt_float_type t = v[i][j];
      v[i][j] = v[j][i];
      v[j][i] = t;
    }
  }

  cerr << "VT matrix\n";
  cerr << v;

  TNT::Array2D<tnt_float_type> t1 = matmult(s, v);
  cerr << "T1 matrix\n";
  cerr << t1;

// Reproduce initial matrix

  TNT::Array2D<tnt_float_type> t2 = matmult(u, t1);
  cerr << "T2 matrix\n";
  cerr << t2;

  return 0;   // gets returned to shell level
}
#endif

static void
add_axes (Molecule & m,
          const Element * axis_element)
{
  assert (NULL != axis_element);

  coord_t xmin, xmax, ymin, ymax, zmin, zmax;

  m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);

  int initial_matoms = m.natoms();

  m.add(axis_element);
  m.setxyz(initial_matoms, 0.0, 0.0, 0.0);

  m.add(axis_element);
  m.setxyz(initial_matoms + 1, xmax * 1.10, 0.0, 0.0);
  m.add_bond(initial_matoms, initial_matoms + 1, SINGLE_BOND);

  m.add(axis_element);
  m.setxyz(initial_matoms + 2, 0.0, ymax * 1.10, 0.0);
  m.add_bond(initial_matoms, initial_matoms + 2, SINGLE_BOND);

  m.add(axis_element);
  m.setxyz(initial_matoms + 3, 0.0, 0.0, zmax * 1.10);
  m.add_bond(initial_matoms, initial_matoms + 3, SINGLE_BOND);

  return;
}

/*
  Compute principal components, optionally align
*/

static int
molecule_pca (Molecule & m,
              IWString_and_File_Descriptor & output)
{
  int matoms = m.natoms();

// First centre the molecule

  TNT::Array2D<tnt_float_type> geometry (matoms, 3);

  if (position_by_extremeties)
    do_position_by_extremeties(m, geometry);
  else
    do_position_by_averages(m, geometry);

#ifdef ECHO_MATRICES
  cerr << "Starting matrix\n";
  cerr << geometry;
#endif

  JAMA::SVD<tnt_float_type> svd(geometry);

  TNT::Array1D<tnt_float_type> singular_values;

  svd.getSingularValues(singular_values);

#ifdef ECHO_MATRICES
  for (int i = 0; i < singular_values.dim(); i++)
  {
    cerr << " sv " << singular_values[i] << '\n';
  }

  TNT::Array2D<tnt_float_type> u;
  svd.getU(u);

  cerr << "U matrix\n";
  cerr << u;
#endif

  TNT::Array2D<tnt_float_type> v;
  svd.getV(v);

#ifdef ECHO_MATRICES
  cerr << "V matrix\n";
  cerr << v;

  TNT::Array2D<tnt_float_type> s;
  svd.getS(s);

  cerr << "S matrix\n";
  cerr << s;
#endif

  if (stream_for_oriented_molecules.active())
  {
    TNT::Array2D<tnt_float_type> t2 = matmult(geometry, v);

    for (int i = 0; i < matoms; i++)
    {
#ifdef ECHO_MATRICES
      cerr << i << " move " << m.x(i) << ' ' << m.y(i) << ' ' << m.z(i) << " to " << t2[i][0] << ' ' << t2[i][1] << ' ' << t2[i][2] << '\n';
#endif
      m.setxyz(i, t2[i][0], t2[i][1], t2[i][2]);
    }

    if (NULL != axis_element)
      add_axes(m, axis_element);

    stream_for_oriented_molecules.write(m);
  }

  append_first_token_of_name(m.name(), output);
  for (int i = 0; i < 3; i++)
  {
    output << ' ' << singular_values[i];
  }

// Note that tubular molecules will have one large and two small
// singular values - different from how shadow areas are done

  float abig = singular_values[0];
  float amid = singular_values[1];
  float amin = singular_values[2];
  double sqrt3 = sqrt(3.0);

  output << ' ' << amin / static_cast<float>(matoms);
  output << ' ';
  output << (sqrt3 - sqrt((1.0-amid/abig)*(1.0-amid/abig) + (1.0-amin/abig)*(1.0-amin/abig)))/sqrt3; // msphere
  output << ' ';
  output << 1.0 - amin / abig;    // anywhere away from the xy plane is flat
  output << ' ';
  output << (sqrt3 - sqrt((amid*amid)/(abig*abig) + (amin*amin)/(abig*abig)))/sqrt3;  // mtube

  output << "\n";

  return output.good ();
}

static int
molecule_pca (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! molecule_pca(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
molecule_pca (const char * fname, FileType input_type, 
                IWString_and_File_Descriptor & output)
{
  assert (NULL != fname);

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return molecule_pca(input, output);
}

static int
molecule_pca(int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lS:o:eX:");

//#define GOLF_EXAMPLE
#ifdef GOLF_EXAMPLE
  return do_golf_example();
#endif

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('e'))
  {
    position_by_extremeties = 1;

    if (verbose)
      cerr << "Initial positions by extremeties\n";
  }

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    axis_element = get_element_from_symbol_no_case_conversion(x);

    if (NULL == axis_element)
    {
      create_element_with_symbol(x);
    }

    if (verbose)
      cerr << "Will add principal components axes, element '" << x << "'\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('S'))
  {
    if (! cl.option_present('o'))
      stream_for_oriented_molecules.add_output_type(FILE_TYPE_SDF);
    else if (! stream_for_oriented_molecules.determine_output_types(cl))
    {
      cerr << "Cannot determine output type(s)\n";
      return 2;
    }

    const_IWSubstring s = cl.string_value('S');

    if (stream_for_oriented_molecules.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 3;
    }

    if (! stream_for_oriented_molecules.new_stem(s))
    {
      cerr << "Cannot initialise stream for oriented molecules '" << s << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Oriented molecules written to '" << s << "'\n";
  }

  set_default_iwstring_float_concatenation_precision(4);
  set_default_iwstring_double_concatenation_precision(4);

#ifdef USE_INTEL_MKL
  mkl_set_num_threads ( 1 );
#endif

  IWString_and_File_Descriptor output(1);

  output << "ID PC1 PC2 PC3 NPC3 psphere pflat ptube\n";

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! molecule_pca(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = molecule_pca (argc, argv);

  return rc;
}
