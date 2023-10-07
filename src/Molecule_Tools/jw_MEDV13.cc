/*
 * This program compute the MEDV-13 descriptors
 * Reference:
 * 1. JCICS, 2001, 41, 321-329, Liu et al
 * 2. JCICS, 1999, 39, 951-957, Liu et al 
 */

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "e_state_computation_procedure_ori.h"

using std::cerr;
using std::endl;

#ifndef MOL2
#define MOL2  15
#endif
#ifndef IW_VDW_MOLVOL
#define IW_VDW_MOLVOL 3;
#endif

static int output_precision = 3;

//DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

// number of atom class
// 0 to 12 is reported in the paper while 13 refer to -I(=)(=)
#define NUMBER_OF_ATOM_CLASS 14
#define NUMBER_OF_ATOM_INTRINSIC_STATE_INDEX 44 // number of atom class
#define UNUSED_ATOM_CLASS -1 // class is from 0 to 12, -1 if not defined
// class is from 0 to 43, -1 is not defined
// 0 to 41 correspond to 1 to 42 reported in the paper 1, 42 refer to -I=(=), 43 refer to =S=
#define UNDEFINED_ATOM_ATTRIBUTE_INDEX -1

static const double intrinsic_state_from_intrinsic_state[NUMBER_OF_ATOM_INTRINSIC_STATE_INDEX] =
{2.0, 1.5, 1.3333, 1.25, 3.0, 2.0, 1.6667, 2.5, 4.0, 2.5, 2.5, 1.75, 1.5, 2.0, 1.6667, 1.8333, 2.4495, 1.8371, 3.6742, 2.2361, 1.6771, 1.4907, 3.3541, 2.2361, 4.4721, 2.7951, 1.9566, 2.2361, 2.2361, 1.7691, 1.1567, 2.3134, 1.134, 1.1227, 2.6458, 1.9108, 1.6536, 1.5345, 1.6149, 1.0559, 0.8696, 0.7764, 0.7937, 1.7010};

static Chemical_Standardisation chemical_standardisation;

int verbose=0;

static int molecules_read = 0;

// keep track of number of warning and number of error
static int number_of_error = 0;

static int ignore_unrecognized_atom_class = 0;

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

  cerr << "  -i <type>      input type\n";

  cerr << "  -f             ignore unrecognized atom class\n";

  (void) display_standard_aromaticity_options (cerr);

  (void) display_standard_chemical_standardisation_options (cerr, 'g');

  cerr << "  -v             verbose output\n";
  
  exit (rc);
}

static int
preprocess_molecule (Molecule & m)
{
  m.reduce_to_largest_fragment ();  // always reduce to largest fragment

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);
  
  return 1;
}

static int
output_result_header (IWString_and_File_Descriptor & output)
{
  output << "Name";

  // print out header for medv descriptor
  int sequence_no = 0;
  for (int i=0; i<NUMBER_OF_ATOM_CLASS; i++)
    for (int j=i; j<NUMBER_OF_ATOM_CLASS; j++)
      {
	if (sequence_no < 100)
	  output<<" medv_JWMEDV"<<sequence_no;
	else
	  output<<" medv_JWMEDVa"<<sequence_no % 10;
	sequence_no ++;
      }

  output << '\n';

  return 1;
}

/*
 * determine the atom class (atom type in the paper) needed
 * atom class from 0 to 12 total of 13, corresponds to 1 to 13 used in the paper 1
 * class 13 refer to -I(=)(=)
 */
static int
determine_atom_class (Molecule &m, int atom_class[],
                      const atomic_number_t * z)
{
  int n_atoms = m.natoms();

  set_vector (atom_class, n_atoms, UNUSED_ATOM_CLASS);

  for (int i=0; i<n_atoms; i++)
    {
      int connections_to_non_hydrogen_atom = m.ncon(i)-m.hcount(i);
      switch (z[i]) {
      case 6:
	switch (connections_to_non_hydrogen_atom) {
	case 1:
	  atom_class[i] = 0;
	  break;
	case 2:
	  atom_class[i] = 1;
	  break;
	case 3:
	  atom_class[i] = 2;
	  break;
	case 4:
	  atom_class[i] = 3;
	  break;
	}
	break;
	
      case 7:
	switch (connections_to_non_hydrogen_atom) {
	case 1:
	  atom_class[i] = 4;
	  break;
	case 2:
	  atom_class[i] = 5;
	  break;
	case 3:
	  atom_class[i] = 6;
	  break;
	case 4:
	  atom_class[i] = 7;
	  break;
	}
	break;
	
      case 8:
	switch (connections_to_non_hydrogen_atom) {
	case 1:
	  atom_class[i] = 8;
	  break;
	case 2:
	  atom_class[i] = 9;
	  break;
	}
	break;
      case 9:
	if (1== (connections_to_non_hydrogen_atom))
	  atom_class[i] = 12;
	break;

      case 15:
	switch (connections_to_non_hydrogen_atom) {
	case 1:
	  atom_class[i] = 4;
	  break;
	case 2:
	  atom_class[i] = 5;
	  break;
	case 3:
	  atom_class[i] = 6;
	  break;
	case 4:
	  atom_class[i] = 7;
	  break;
	}
	break;
	
      case 16:
	switch (connections_to_non_hydrogen_atom) {
	case 1:
	  atom_class[i] = 8;
	  break;
	case 2:
	  atom_class[i] = 9;
	  break;
	case 3:
	  atom_class[i] = 10;
	  break;
	case 4:
	  atom_class[i] = 11;
	  break;
	}
	break;

      case 17:
	if (1== (connections_to_non_hydrogen_atom))
	  atom_class[i] = 12;
	break;

      case 35:
	if (1== (connections_to_non_hydrogen_atom))
	  atom_class[i] = 12;
	break;

      case 53:
	if (1== (connections_to_non_hydrogen_atom))
	  atom_class[i] = 12;
	else if ((3== (connections_to_non_hydrogen_atom)) && (5==m.nbonds(i)))
	  atom_class[i] = 13;
	break;
      }
      if ((1 != m.atomic_number(i)) && (UNUSED_ATOM_CLASS == atom_class[i]))
	return 0;
    }
  return 1;
}


static void
determine_atom_intrinsic_state_index (Molecule &m, int a_index[],
                                      const atomic_number_t * z)
{
  int n_atoms = m.natoms();

  set_vector (a_index, n_atoms, UNDEFINED_ATOM_ATTRIBUTE_INDEX);

  for (int i=0; i<n_atoms; i++)
  {
    if ((m.formal_charge (i) > 0) && (m.hcount(i) < m.formal_charge(i)))
      ;
    else
    {
      int h_count = m.hcount(i) - m.formal_charge(i);
      int nconi = m.ncon(i) - m.formal_charge(i);
      int nbondsi = m.nbonds(i) - m.formal_charge(i);

      switch (z[i]) {
        case 6:
          if (nconi==nbondsi)
            switch (h_count) {
            case 3:
              a_index[i] = 0;
              break;
            case 2:
              a_index[i] = 1;
              break;
            case 1:
              a_index[i] = 2;
              break;
            case 0:
              a_index[i] = 3;
              break;
            }
          else if (nconi == nbondsi-1)
            switch (h_count) {
            case 2:
              if (number_of_bond_in_conjugated_system(m, i))
                a_index[i] = 10;
              else
                a_index[i] = 4;
              break;
            case 1:
              switch (number_of_bond_in_conjugated_system(m, i)) {
              case 2:
                a_index[i] = 13;
                break;
              case 1:
                a_index[i] = 11;
                break;
              case 0:
                a_index[i] = 5;
                break;
              }
              break;
            case 0:
              switch (number_of_bond_in_conjugated_system(m, i)) {
              case 3:
                a_index[i] = 15;
                break;
              case 2:
                a_index[i] = 14;
                break;
              case 1:
                a_index[i] = 12;
                break;
              case 0:
                a_index[i] = 6;
                break;
              }
              break;
            }
          else if (nconi == nbondsi-2)
            {
              if (is_triple_bonded (m, i))
                {
                  if (1==h_count)
                    a_index[i] = 8;
                  else if (0==h_count)
                    a_index[i] = 9;
                }
              else a_index[i] = 7;
            }
          
          break;
          
        case 7:
          if (nconi == nbondsi)
            switch (h_count) {
            case 2:
              a_index[i] = 19;
              break;
            case 1:
              a_index[i] = 20;
              break;
            case 0:
              a_index[i] = 21;
              break;
            }
          else if (nconi == nbondsi - 1)
            switch (h_count) {
            case 1:
              if (number_of_bond_in_conjugated_system(m, i))
                a_index[i] = 22;
              else
                a_index[i] = 25;
              break;
            case 0:
              switch (number_of_bond_in_conjugated_system(m, i)) {
              case 2:
                a_index[i] = 27;
                break;
              case 1:
                a_index[i] = 26;
                break;
              case 0:
                a_index[i] = 23;
                break;
              }
              break;
            }
          else if ((1 == nconi) && (3 == nbondsi))
            a_index[i] = 24;
          else if ((3 == nconi) && (5 == nbondsi))
            a_index[i] = 28;
          break;
          
        case 8:
          if (nconi == nbondsi)
            switch (h_count) {
            case 1:
              a_index[i] = 16;
              break;
            case 0:
              a_index[i] = 17;
              break;
            }
          else if (nconi == nbondsi-1)
            a_index[i] = 18;
          break;
          
        case 9:
          if (1== nconi)
            a_index[i] = 34;
          break;
          
        case 15:
          if (nconi == nbondsi)
            switch (nconi - h_count) {
            case 1:
              a_index[i] = 38;
              break;
            case 2:
              a_index[i] = 39;
              break;
            case 3:
              a_index[i] = 40;
              break;
            }
          else if (4 == nconi)
            a_index[i] = 41;
          
          break;
          
        case 16:
          if (nconi == nbondsi)
            switch (h_count) {
            case 1:
              a_index[i] = 29;
          break;
            case 0:
              a_index[i] = 30;
              break;
            }
          else if ((1== nconi) && (2== nbondsi))
            a_index[i] = 31;
          else if ((3== nconi) && (4== nbondsi))
            a_index[i] = 32;
          else if ((4== nconi) && (6== nbondsi))
            a_index[i] = 33;
          else if ((2== nconi) && (4== nbondsi))
            a_index[i] = 43;
          break;
          
        case 17:
          if (1== nconi)
            a_index[i] = 35;
          break;
          
        case 35:
          if (1== nconi)
            a_index[i] = 36;
          break;
          
        case 53:
          if (1== nconi)
            a_index[i] = 37;
          else if ((3==nconi) && (5==nbondsi))
            a_index[i] = 42;
          break;
        }
      }
  }
}

static int
value_of_atom_intrinsic_state (Molecule &m, const int i, double & value,
                               const atomic_number_t * z)
{
  double n =0.0;
  double v =0.0;
  switch (z[i]) {
  case 1:
    n=1.0;
    v=1.0;
    break;
  case 6:
    n=2.0;
    v=4.0;
    break;
  case 7:
    n=2.0;
    v=5.0;
    break;
  case 8:
    n=2.0;
    v=6.0;
    break;
  case 9:
    n=2.0;
    v=7.0;
    break;
  case 15:
    n=3.0;
    v=5.0;
    break;
  case 16:
    n=3.0;
    v=6.0;
    break;
  case 17:
    n=3.0;
    v=7.0;
    break;
  case 35:
    n=4.0;
    v=7.0;
    break;
  case 53:
    n=5.0;
    v=7.0;
    break;
  }

  double delta = m.ncon(i)-m.hcount(i);
  double delta_v = m.nbonds(i)-m.hcount(i);

  if (0.0==delta) return 0;
  else 
    value = sqrt(v/4.0) * ((2.0/n)*(2.0/n)*delta_v +1.0)/delta;
  return 1;
}
 
static int
determine_atom_intrinsic_state (Molecule &m, const int i_state_index[], double i_state[],
                                const atomic_number_t * z)
{
  int n_atoms = m.natoms();

  for (int i=0; i<n_atoms; i++)
    {
      if ((1!=z[i]) && (-1 == i_state_index[i]))
        {
          double value;
          if (!value_of_atom_intrinsic_state (m, i, value, z)) return 0;
          i_state[i] = value;
        }
      else
        i_state[i] = intrinsic_state_from_intrinsic_state[i_state_index[i]];
    }
  return 1;
}

/* 
 * compute the descriptors
 */
static int 
compute_medv_descriptors (Molecule &m, 
                          IWString_and_File_Descriptor & output,
                          const int atom_class[],
                          const double e_state_index[],
                          double * descriptors)
{
  std::fill_n(descriptors, NUMBER_OF_ATOM_CLASS * NUMBER_OF_ATOM_CLASS, 0.0);

  const int n_atoms = m.natoms();

#ifdef DEBUG_JWMEDV
  for (auto i = 0; i < n_atoms; ++i)
  {
    cerr << " atom " << i << " e_state_index " << e_state_index[i] << " class " << atom_class[i] << endl;
  }
#endif

  for (int i=0; i<n_atoms; i++)
  {
    const int aci = atom_class[i];

    if (UNUSED_ATOM_CLASS != aci)
    {
      for (int j=i+1; j<n_atoms; j++)
      {
        if (UNUSED_ATOM_CLASS != atom_class[j])
        {
          const int graph_distance = m.bonds_between(i, j);
          descriptors[aci * NUMBER_OF_ATOM_CLASS + atom_class[j]] += (e_state_index[i] * e_state_index[j]) / static_cast<double>(graph_distance * graph_distance); 
        }
      }
    }
  }
  
  append_first_token_of_name(m.name(), output);

  for (int i=0; i<NUMBER_OF_ATOM_CLASS; i++)
  {
    output << ' ' << static_cast<float>(descriptors[i * NUMBER_OF_ATOM_CLASS + i]);

    for (int j=i+1; j<NUMBER_OF_ATOM_CLASS; j++)
      output << ' ' << static_cast<float>(descriptors[i * NUMBER_OF_ATOM_CLASS + j] + descriptors[j * NUMBER_OF_ATOM_CLASS + i]); 
  }

  output << '\n';

  return 1;
}

/* 
 * This function is to calculate the MEDV descriptors for molecule
 */

static int
jwmedv (Molecule &m, IWString_and_File_Descriptor & output, 
        int atom_class[], int atom_intrinsic_state_index[], 
        double intrinsic_state[], double e_state_index[], 
        double * descriptors)
{
  const int matoms = m.natoms();

  atomic_number_t * z = new atomic_number_t[matoms]; std::unique_ptr<int[]> free_z(z);
  m.atomic_numbers(z);

  if (! determine_atom_class (m, atom_class, z))
  {
    if (ignore_unrecognized_atom_class) 
    {
      if (verbose)
        cerr<<"Warning"<<"\tMolecule "<<m.name()<<"\tUnrecognizable molecule structure, ignored"<<endl;
    }
    else return 0;
  }

  determine_atom_intrinsic_state_index (m, atom_intrinsic_state_index, z);

  if (! determine_atom_intrinsic_state (m, atom_intrinsic_state_index, intrinsic_state, z))
    return 0;

#ifdef DEBUG_JWMEDV
  for (auto i = 0; i < matoms; ++i)
  {
    cerr << "atom_intrinsic_state_index " << i << " value " << atom_intrinsic_state_index[i] << endl;
  }
#endif
  
  if (! determine_atom_e_state_index (m, e_state_index, intrinsic_state, z))
    return 0;

#ifdef DEBUG_JWMEDV
  for (auto i = 0; i < matoms; ++i)
  {
    cerr << " e_state_index " << i << " value " << e_state_index[i] << endl;
  }
#endif
  
  if (!compute_medv_descriptors (m, output, atom_class, e_state_index, descriptors)) 
    return 0;

  return 1;
}

/* 
 * Introduce and delete arrays to separate them from the main program
 */

static int
jwmedv (Molecule &m, IWString_and_File_Descriptor & output)
{
  // all hydrogen is needed for the computation of connections and bonds
  m.make_implicit_hydrogens_explicit();
  
  int n_atoms = m.natoms();
  int * atom_class = new int[n_atoms]; std::unique_ptr<int[]> free_atom_class(atom_class);

  int * atom_intrinsic_state_index = new int[n_atoms]; std::unique_ptr<int[]> free_atom_intrinsic_state_index(atom_intrinsic_state_index);

  double * atom_intrinsic_state = new double[n_atoms]; std::unique_ptr<double[]> free_atom_intrinsic_state(atom_intrinsic_state);
  double * electrotopological_state_index = new double[n_atoms]; std::unique_ptr<double[]> free_electrotopological_state_index(electrotopological_state_index);

  double medv_descriptors[NUMBER_OF_ATOM_CLASS * NUMBER_OF_ATOM_CLASS];

  auto rc=jwmedv (m, output, atom_class, atom_intrinsic_state_index, atom_intrinsic_state, electrotopological_state_index, medv_descriptors);

  if (0==rc) 
  {
    number_of_error++;
    
    if (verbose)
    {
      cerr<<"ERROR"<<"\tMolecule "<<m.name()<<"\tUnrecognizable molecule structure, notify Jibo"<<endl;
    }
  }

  return rc;
}

/* 
 * deal with each molecule
 */
static int
jwmedv (data_source_and_type<Molecule> & input, 
        IWString_and_File_Descriptor & output)
{
  set_default_iwstring_float_concatenation_precision(output_precision);

  Molecule * m;
  while (nullptr != (m = input.next_molecule ()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose)
      cerr << molecules_read << " processing '" << m->name () << "'\n";

    (void) preprocess_molecule (*m);

    if (!  jwmedv (*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

/*
 * deal with each files
 */

static int
jwmedv (const char * fname, FileType input_type, 
        IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input (input_type, fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jwmedv (input, output);
}

/*
 * deal with all the commandline options
 */

static int
jwmedv (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:fg:p:");
  
  if (cl.unrecognised_options_encountered ())
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
  
  set_global_aromaticity_type (Pearlman);
  
  ignore_unrecognized_atom_class = cl.option_count ('f');

  if (cl.option_present ('g'))   // only recognised with the -d option
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation (-g option)\n";
      usage (8);
    }
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

  if (0 == cl.number_elements ())
    {
      cerr << "Insufficient arguments\n";
      usage (2);
    }
  
  IWString_and_File_Descriptor output(1);

  if (! output_result_header (output))
    return 0;

  int rc = 0;
  
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! jwmedv (cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  // output info

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }
  
  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = jwmedv (argc, argv);
  return rc;
}
