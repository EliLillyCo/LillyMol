/*
 * This program provide the computation of 3D-Morse descriptors
 * Reference:
 * JCICS, 1996, 36, 334-344 by Shuur, Selzer and Gasteiger
 * Anal. Chem., 1997, 69, 2398 - 2406 by Schuur and Gasteiger
 */

#include <assert.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

#ifndef MOL2
#define MOL2  15
#endif
#ifndef IW_VDW_MOLVOL
#define IW_VDW_MOLVOL 3;
#endif

//DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

#define NATMTYPE 10 // type of atom
#define MAX_MORSE_DESCRIPTORS 32 // number of morse descirptors for each properties
#define NUMBER_OF_PROPERTY 6 // number of properties computed

//  When writing to a database we need chemical standardisation

static Chemical_Standardisation chemical_standardisation;

int verbose=0;

static int molecules_read = 0;

// keep track of number of warning and number of error

static int number_of_error = 0;

static int first_open = 1;

static int output_precision = 4;

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
  cerr << "  -p <n>         output precision (default " << output_precision << ")\n";
  cerr << "  -i <type>      input type, enter '-i help' for more info\n";
  (void) display_standard_aromaticity_options (cerr);
  (void) display_standard_chemical_standardisation_options (cerr, 'g');

#ifdef USE_IWMALLOC
  display_standard_debug_options (cerr, 'd');
#endif
  cerr << "  -v             verbose output\n";
  
  exit (rc);
}

static int
preprocess_molecule (Molecule & m)
{
  m.reduce_to_largest_fragment();  // always reduce to largest fragment
  
  return 1;
}

static int
output_result_header (IWString_and_File_Descriptor & output)
{
  output << "Name";

  for (int i=0; i< MAX_MORSE_DESCRIPTORS; i++)
  {
  /* commented out at this time because of the size for the name of descriptors
    buffer <<"JWMORS"<<i<<"ATNO ";
    buffer <<"JWMORS"<<i<<"ATMASS ";
    buffer <<"JWMORS"<<i<<"VDWV ";
    buffer <<"JWMORS"<<i<<"ELN ";
    buffer <<"JWMORS"<<i<<"POL ";
    buffer <<"JWMORS"<<i<<"CRG ";
    end of comment*/

  // at this time, use the short version
    output <<" morse_JWMO"<<i<<"AN";
    output <<" morse_JWMO"<<i<<"AM";
    output <<" morse_JWMO"<<i<<"VV";
    output <<" morse_JWMO"<<i<<"EN";
    output <<" morse_JWMO"<<i<<"PL";
    output <<" morse_JWMO"<<i<<"CG";
  }

  output << '\n';

  return 1;
}

/*
 * This procedure fills the bcut_matrix for the later computation
 * need to be changed to allow different ways of assigning the value
 */

#ifdef JIBOS_VERSION
static void
get_property_for_each_atom (Molecule &m, double ** property_array) 
{
  // this struct entry is used in assiging the weight for at_whim
  struct entry 
  {
    char    TYPE[6];
    double   ATNO;
    double   ATMASS;
    double   VDWV;
    double   ELN;
    double   POL;
  };

  // The value is adjusted so that the maximum is 1  
  struct entry AtmTypeList[NATMTYPE]=
  {
    // original -- before adjustment
    //    { "H", 1,  0.084, 0.299,  0.944, 0.379 },
    //    { "C", 6,  1.000, 1.000,  1.000, 1.000 },
    //    { "N", 7,  1.166, 0.695,  1.163, 0.625 },
    //    { "O", 8,  1.332, 0.512,  1.331, 0.456 },
    //    { "F", 9,  1.582, 0.410,  1.457, 0.316 },
    //    { "P", 15, 2.579, 1.181,  0.916, 2.063 },
    //    { "S", 16, 2.670, 1.088,  1.077, 1.648 },
    //    { "Cl",17, 2.952, 1.035,  1.265, 1.239 },
    //    { "Br",35, 6.653, 1.384,  1.172, 1.733 },
    //    { "I", 53, 10.556, 1.728,  1.012, 3.040 }
    //After adjustment
    { "H", 0.01887, 0.00796, 0.17303, 0.64791, 0.12467},
    { "C", 0.11321, 0.09437,0.57870, 0.68634, 0.32895},
    { "N", 0.13208, 0.11046, 0.40220, 0.79822, 0.20559},
    { "O", 0.15094, 0.12618, 0.29630, 0.91352, 0.15000},
    { "F", 0.16981, 0.14987, 0.23727, 1.00000, 0.10395},
    { "P", 0.28302, 0.24432, 0.68345, 0.62869, 0.67862},
    { "S", 0.30189, 0.25294, 0.62963, 0.73919, 0.54211},
    { "Cl",0.32076, 0.27965, 0.59896, 0.86822, 0.40757},
    { "Br",0.66038, 0.63026, 0.80093, 0.80439, 0.57007},
    { "I", 1.00000, 1.00000, 1.00000, 0.69458, 1.00000}
  };
  
  int n_atoms = m.natoms();

  // set up array to put the ghose crippen index value and hydrogen bond donor acceptor array

  // fill in the diagonal elements
  for (int i=0; i<n_atoms; i++)
    {
      int j = 0;
      while (j < NATMTYPE)
	{
	  if (AtmTypeList[j].TYPE == m.atomi(i)->atomic_symbol())
	    {
	      property_array[i][0] = AtmTypeList[j].ATNO;
	      property_array[i][1] = AtmTypeList[j].ATMASS;
	      property_array[i][2] = AtmTypeList[j].VDWV;
	      property_array[i][3] = AtmTypeList[j].ELN;
	      property_array[i][4] = AtmTypeList[j].POL;
	      break;
	    }
	  else j++;
	}
    }
}
#endif

/*
  Version from Ian that uses a switch on the atomic number. Runs a little
  faster by avoiding string comparisons
*/

static int
get_property_for_each_atom (const Molecule & m,
                            double ** property_array)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);
    switch (z)
    {
      case 1:
        property_array[i][0] = 0.01887;
        property_array[i][1] = 0.00796;
        property_array[i][2] = 0.17303;
        property_array[i][3] = 0.64791;
        property_array[i][4] = 0.12467;
        break;
      case 6:
        property_array[i][0] = 0.11321;
        property_array[i][1] = 0.09437;
        property_array[i][2] = 0.57870;
        property_array[i][3] = 0.68634;
        property_array[i][4] = 0.32895;
        break;
      case 7:
        property_array[i][0] = 0.13208;
        property_array[i][1] = 0.11046;
        property_array[i][2] = 0.40220;
        property_array[i][3] = 0.79822;
        property_array[i][4] = 0.20559;
        break;
      case 8:
        property_array[i][0] = 0.15094;
        property_array[i][1] = 0.12618;
        property_array[i][2] = 0.29630;
        property_array[i][3] = 0.91352;
        property_array[i][4] = 0.15000;
        break;
      case 9:
        property_array[i][0] = 0.16981;
        property_array[i][1] = 0.14987;
        property_array[i][2] = 0.23727;
        property_array[i][3] = 1.00000;
        property_array[i][4] = 0.10395;
        break;
      case 15:
        property_array[i][0] = 0.28302;
        property_array[i][1] = 0.24432;
        property_array[i][2] = 0.68345;
        property_array[i][3] = 0.62869;
        property_array[i][4] = 0.67862;
        break;
      case 16:
        property_array[i][0] = 0.30189;
        property_array[i][1] = 0.25294;
        property_array[i][2] = 0.62963;
        property_array[i][3] = 0.73919;
        property_array[i][4] = 0.54211;
        break;
      case 17:
        property_array[i][0] = 0.32076;
        property_array[i][1] = 0.27965;
        property_array[i][2] = 0.59896;
        property_array[i][3] = 0.86822;
        property_array[i][4] = 0.40757;
        break;
      case 35:
        property_array[i][0] = 0.66038;
        property_array[i][1] = 0.63026;
        property_array[i][2] = 0.80093;
        property_array[i][3] = 0.80439;
        property_array[i][4] = 0.57007;
        break;
      case 53:
        property_array[i][0] = 1.00000;
        property_array[i][1] = 1.00000;
        property_array[i][2] = 1.00000;
        property_array[i][3] = 0.69458;
        property_array[i][4] = 1.00000;
        break;
      default:
        cerr << "Unrecognised atom type '" << m.atomic_symbol(i) << "'\n";
        return 0;
    }
  }

  return 1;
}

static int
get_atomic_properties (Molecule & m,
                       double * p)
{
  const int matoms = m.natoms();

  if (! m.compute_Gasteiger_partial_charges())
  {
    cerr << "Cannot compute partial charges for '" << m.name() << "'\n";
    return 0;
  }

  for (int i = 0; i < matoms; ++i)
  {
    p[i * NUMBER_OF_PROPERTY + 5] = m.charge_on_atom(i);

    atomic_number_t z = m.atomic_number(i);
    switch (z)
    {
      case 1:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.01887;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.00796;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.17303;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.64791;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.12467;
        break;
      case 6:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.11321;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.09437;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.57870;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.68634;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.32895;
        break;
      case 7:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.13208;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.11046;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.40220;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.79822;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.20559;
        break;
      case 8:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.15094;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.12618;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.29630;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.91352;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.15000;
        break;
      case 9:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.16981;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.14987;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.23727;
        p[i * NUMBER_OF_PROPERTY + 3] = 1.00000;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.10395;
        break;
      case 15:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.28302;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.24432;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.68345;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.62869;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.67862;
        break;
      case 16:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.30189;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.25294;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.62963;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.73919;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.54211;
        break;
      case 17:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.32076;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.27965;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.59896;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.86822;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.40757;
        break;
      case 35:
        p[i * NUMBER_OF_PROPERTY + 0] = 0.66038;
        p[i * NUMBER_OF_PROPERTY + 1] = 0.63026;
        p[i * NUMBER_OF_PROPERTY + 2] = 0.80093;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.80439;
        p[i * NUMBER_OF_PROPERTY + 4] = 0.57007;
        break;
      case 53:
        p[i * NUMBER_OF_PROPERTY + 0] = 1.00000;
        p[i * NUMBER_OF_PROPERTY + 1] = 1.00000;
        p[i * NUMBER_OF_PROPERTY + 2] = 1.00000;
        p[i * NUMBER_OF_PROPERTY + 3] = 0.69458;
        p[i * NUMBER_OF_PROPERTY + 4] = 1.00000;
        break;
      default:
        cerr << "Unrecognised atom type '" << m.atomic_symbol(i) << "'\n";
        return 0;
    }
  }

  return 1;
}

static int
jwmorse (Molecule & m,
         const double * p,
         IWString_and_File_Descriptor & output)
{
  const int matoms = m.natoms();

  double * morse = new double[NUMBER_OF_PROPERTY * MAX_MORSE_DESCRIPTORS]; std::unique_ptr<double[]> free_morse(morse);
  std::fill_n(morse, NUMBER_OF_PROPERTY * MAX_MORSE_DESCRIPTORS, 0.0);

  double * sd = new double[MAX_MORSE_DESCRIPTORS]; std::unique_ptr<double[]> free_sd(sd);

  for (int i = 0; i < matoms; ++i)
  {
    const double * pi = p + (i * NUMBER_OF_PROPERTY);

    for (int j = i + 1; j < matoms; ++j)
    {
      const double * pj = p + (j * NUMBER_OF_PROPERTY);

      const double d = m.distance_between_atoms(i, j);
      sd[0] = 1.0;
      for (int s = 1; s < MAX_MORSE_DESCRIPTORS; ++s)            // very important to get a vectorised sin
      {
        sd[s] = sin(s * d) / (s * d);
      }

      for (int k = 0; k < NUMBER_OF_PROPERTY; ++k)
      {
        morse[k] += pi[k] * pj[k];
      }

      for (int s = 1; s < MAX_MORSE_DESCRIPTORS; ++s)
      {
        for (int k = 0; k < NUMBER_OF_PROPERTY; ++k)
        {
          morse[s * NUMBER_OF_PROPERTY + k] += sd[s] * pi[k] * pj[k];
        }
      }
    }
  }

  append_first_token_of_name(m.name(), output);

  for (int i = 0; i < (MAX_MORSE_DESCRIPTORS*NUMBER_OF_PROPERTY); ++i)
  {
    output << ' ' << static_cast<float>(morse[i]);
  }

  output << '\n';

  return 1;
}

#ifdef WORKING_VERSIONQWREQW
static int
jwmorse (Molecule & m,
         const double * p,
         IWString_and_File_Descriptor & output)
{
  const int matoms = m.natoms();

  double * morse = new double[NUMBER_OF_PROPERTY * MAX_MORSE_DESCRIPTORS]; unique_ptr<double[]> free_morse(morse);
  std::fill_n(morse, NUMBER_OF_PROPERTY * MAX_MORSE_DESCRIPTORS, 0.0);

  for (int i = 0; i < matoms; ++i)
  {
    const double * pi = p + (i * NUMBER_OF_PROPERTY);

    for (int j = i + 1; j < matoms; ++j)
    {
      const double d = m.distance_between_atoms(i, j);
    
      const double * pj = p + (j * NUMBER_OF_PROPERTY);

      for (int k = 0; k < NUMBER_OF_PROPERTY; ++k)
      {
        morse[k] += pi[k] * pj[k];
      }

      for (int s = 1; s < MAX_MORSE_DESCRIPTORS; ++s)
      {
        const double tmp = sin(s*d)/(s*d);

        for (int k = 0; k < NUMBER_OF_PROPERTY; ++k)
        {
          morse[s * NUMBER_OF_PROPERTY + k] += tmp * pi[k] * pj[k];
        }
      }
    }
  }

  append_first_token_of_name(m.name(), output);

  for (int i = 0; i < (MAX_MORSE_DESCRIPTORS*NUMBER_OF_PROPERTY); ++i)
  {
    output << ' ' << static_cast<float>(morse[i]);
  }

  output << '\n';

  return 1;
}
#endif

static int
jwmorse (Molecule & m,
         IWString_and_File_Descriptor & output)
{
  const int matoms = m.natoms();

  double * p = new double[matoms * NUMBER_OF_PROPERTY]; std::unique_ptr<double[]> free_p(p);

  if (! get_atomic_properties(m, p))
  {
    cerr << "jwmorse:cannot compute atomic properties for '" << m.name() << "'\n";
    return 0;
  }

  return jwmorse(m, p, output);
}

/* 
 */

static int
jwmorse (Molecule & m,
         IWString_and_File_Descriptor & output,
         double ** distance_matrix,
         double ** intermediate_matrix,
	 double ** property_array,
         double atomic_charge[],
         double morse_descriptors[])
{
  if (m.name().contains(' '))
  {
    const_IWSubstring tmp(m.name());
    tmp.truncate_at_first(' ');
    output << tmp;
  }
  else
    output << m.name();

  /* HERE is the main computation function */  

  int n_atoms = m.natoms();

  // fill the distance matrix
  for (int i=0; i<n_atoms; i++)
  {
    for (int j=i+1; j<n_atoms; j++)
      {
	distance_matrix[i][j] = m.distance_between_atoms (i, j);
      }
  }

  set_default_iwstring_float_concatenation_precision(output_precision);

  get_property_for_each_atom (m, property_array);

  for (int i=0; i<n_atoms; i++)
  {
    property_array[i][5] = atomic_charge[i];
  }

  for (int s=0; s<MAX_MORSE_DESCRIPTORS; s++)
  {
    if (0 == s) 
    {
      for (int i=0; i<n_atoms; i++)
      {
        for (int j=i+1; j<n_atoms; j++)
        {
          intermediate_matrix[i][j] = 1.0;
        }
      }
    }
    else
    {
      for (int i=0; i<n_atoms; i++)
      {
        for (int j=i+1; j<n_atoms; j++)
        {
          double temp = distance_matrix[i][j] * s;
          intermediate_matrix[i][j] = sin (temp)/temp;
        }
      }
    }
  
    for (int k=0; k<NUMBER_OF_PROPERTY; k++)
    {
      morse_descriptors[k] = 0.0;
      for (int i=0; i<n_atoms; i++)
      {
        for (int j=i+1; j<n_atoms; j++)
        {
          morse_descriptors[k] += intermediate_matrix[i][j] * property_array[i][k] * property_array[j][k];
        }
      }
    }
   
    for (int i=0; i<NUMBER_OF_PROPERTY; i++)
    {
      output << ' ' << static_cast<float>(morse_descriptors[i]);
    }
  }

  output << '\n';

  return 1;
}


#ifdef NOT_BEING_USED
static int
jwmorse (Molecule & m,
         IWString_and_File_Descriptor & output,
         std::ofstream &logfile )
{
  int rc=0;
  // calculate_partial_charge need to have all hydrogen added
  //m.make_implicit_hydrogens_explicit();

  int n_atoms = m.natoms();

  double * morse_descriptors = new double [NUMBER_OF_PROPERTY]; std::unique_ptr<double[]> free_morse_descriptors(morse_descriptors);

  double * atomic_charge = new double [n_atoms]; std::unique_ptr<double[]> free_atomic_charge(atomic_charge);

  double ** property_array = new double * [n_atoms];
  double ** distance_matrix = new double * [n_atoms];
  double ** intermediate_matrix = new double * [n_atoms];

  for (int i=0; i< n_atoms; i++)
    {
      distance_matrix[i] = new double [n_atoms];
      intermediate_matrix[i] = new double [n_atoms];
      property_array[i] = new double [NUMBER_OF_PROPERTY];
    }

  //  will switch to Gasteiger partial charge in the next release
  //  rc = m.compute_Gasteiger_partial_charges();
  //  if (0==rc) 
  //    {
  //      logfile<<m.name()<<"\tERROR in calculation of Gasteiger partial charge\n";
  //      if (verbose) cerr<<m.name()<<"\tERROR in calculation of Gasteiger partial charge\n";
  //    }
  
  //  rc = m.compute_Gasteiger_Huckel_partial_charges();

  rc = m.compute_Gasteiger_partial_charges();
  
  if (0==rc) 
  {
    if (verbose)
      { 
      cerr<<m.name()<<"\tERROR in calculation of Gasteiger partial charge\n";
      logfile<<m.name()<<"\tERROR in calculation of Gasteiger partial charge\n";
    }
  }
  else
  {
    for (int i=0; i<n_atoms; i++)
    {
      atomic_charge[i] = m.charge_on_atom (i);
    }

    //      for (int i=0; i<n_atoms; i++) cerr<<"atom i="<<i<<"\tcharge="<<atomic_charge[i]<<endl;

    rc=jwmorse (m, output, distance_matrix, intermediate_matrix, property_array, atomic_charge, morse_descriptors);
      
    if (0==rc) 
      number_of_error++;
  }

  for (int i=0; i<n_atoms; i++)
  {
    delete [] distance_matrix[i];
    delete [] intermediate_matrix[i];
    delete [] property_array[i];
  }
  delete [] intermediate_matrix;
  delete [] distance_matrix;
  delete [] property_array;

  return rc;
}
#endif

static int
jwmorse (data_source_and_type<Molecule> & input,
         IWString_and_File_Descriptor & output,
         std::ofstream & logfile)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose)
      cerr << molecules_read << " processing '" << m->name() << "'\n";

    (void) preprocess_molecule (*m);

    if (! jwmorse (*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
jwmorse (const char * fname, 
         FileType input_type,
         IWString_and_File_Descriptor & output,
         std::ofstream &logfile)
{
  data_source_and_type<Molecule> input (input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    if (verbose)
      logfile<<"\nCannot open "<<fname<<endl<<endl;
    return 0;
  }

  // print out the log file
  if (verbose)
    logfile<<"\nprocessing "<<fname<<endl<<endl;
  
  // print out the headerfile

  if (first_open)
  {
    if (! output_result_header (output))
      return 0;

    first_open =0;
  }

  return jwmorse (input, output, logfile);
}

static int
jwmorse (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:i:p:ng:");
  
  if (cl.unrecognised_options_encountered())
    {
      cerr << "Unrecognised options encountered\n";
      usage (1);
    }
  
  verbose = cl.option_count ('v');
  
  if (! process_elements (cl))
    {
      usage (2);
    }
  
  if (! process_standard_aromaticity_options (cl, verbose))
    {
      cerr << "Cannot process aromaticity options (-A)\n";
      usage (5);
    }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', output_precision) || output_precision < 1)
    {
      cerr << "The output precision must be a whole, +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Results written with " << output_precision << " precision\n";
  }
  
  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present ('i'))
    {
      if (! process_input_type (cl, input_type))
        {
          cerr << "Cannot determine input type\n";
          usage (6);
        }
    }
  else if (! all_files_recognised_by_suffix (cl))
    return 4;
  
  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }
  
  set_default_iwstring_float_concatenation_precision(output_precision);

  int rc = 0;
  
  std::ofstream logfile;
  if (verbose)
    logfile.open("jwmorse_descriptor.log", std::ios::out);

  if (!logfile) 
    {
      cerr << "jwmorse_descriptor.log file cannot be opened\n";
      return 0;
    }

  time_t current_time;
  if (verbose)
    {
      logfile<<"This file collect all the info about error during the calculation\n";
#ifdef GIT_HASH
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
      time(&current_time);
      logfile<<"calculation started at "<<ctime(&current_time);
      logfile<<"The complete Command was:"<<endl;
      for (int kk=0;kk<argc;kk++)
        {
          logfile<<argv[kk]<<" ";
        }
      logfile<<endl;
    }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! jwmorse (cl[i], input_type, output, logfile))
        {
          rc = i + 1;
          //      break;
        }
    }

    output.flush();

  // output info
  if (verbose)
    {
      time(&current_time);
      logfile<<"\n\ncalculation ended at "<<ctime(&current_time)<<endl;
      logfile<<"Total Molecules read: "<<molecules_read<<"\tTotal Error: "<<number_of_error<<endl;
      logfile<<"Error Rate: "<<(double) number_of_error/ (double) molecules_read<<endl;
      if (verbose)
        {
          cerr << "Read " << molecules_read << " molecules\n";
        }
    }
  
  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = jwmorse (argc, argv);
  return rc;
}
