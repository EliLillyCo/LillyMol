/*
 * This program is to compute the DiP (Distance Profiles) descriptors
 * Baumann, 507-519 QSAR 21 (2002)
 */

#include <stdlib.h>
#include <assert.h>

#include <iostream>
//#include <sstream>
#include <iomanip>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iw_vdw.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/surface_area_molvol.h"

using std::cerr;
using std::endl;

//DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

#define DIP_NATMTYPE 11 // type of DiP atom type

static int total_number_of_descriptors = 0;

static int number_of_regions = 10;
static double region_width = 1.0;

//  When populating a descriptor database, we may not want the usual
//  descriptor output (turn off with the -n option)

static int write_results_to_output_file = 1;

static IWDigits iwdigits;

int verbose=0;

static int molecules_read = 0;

// keep track of number of warning and number of error

static int number_of_error = 0;

static int compress_heavy_halogens = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Compute the JWDIP descriptors for molecules\n";
  cerr << "  -r <int>       number of region to do the computation (default 10) 1-20\n";
  cerr << "  -w <number>    region width (default 1.0) 0.2-2.0\n";
  cerr << "  -y             compress heavy halogens (I=Br=Cl)\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -A ...         standards aromaticity definitions\n";
  cerr << "  -E ...         standards element definitions\n";
  cerr << "  -v             verbose output\n";
  
  exit(rc);
}

static int
preprocess_molecule (Molecule & m)
{

  m.reduce_to_largest_fragment();  // always reduce to largest fragment

  m.remove_all(1);     // we do not do anything with Hydrogens
  
  return 1;
}

static char
character_type (int i)
{
  char c;
  if (0 == i) c = 'Y';
  else if (1==i) c = 'N';
  else if (2==i) c = 'O';
  else if (3==i) c = 'F';
  else if (4==i) c = 'P';
  else if (5==i) c = 'S';
  else if (6==i) c = 'L';
  else if (7==i) c = 'B';
  else if (8==i) c = 'I';
  else if (9==i) c = 'D';
  else if (10==i) c = 'T';
  else c= '0';

  return c;
}

static int
output_result_header (IWString_and_File_Descriptor & output)
{
  output << "Name";
  // print out header for bcut descriptor
  for (int i=0; i<DIP_NATMTYPE; i++)
    for (int j=i; j<DIP_NATMTYPE; j++)
      for (int k=0; k<number_of_regions +1; k++)
        output << " jwdp_JWDP"<<character_type(i)<<character_type(j)<<k;

  output << '\n';

  return 1;
}

static int find_DiP_atom_type_for_atoms (Molecule &m, int atom_type [])
{
  int n_atoms = m.natoms();
  for (int i=0; i<n_atoms; i++)
    {
      const Atom * atomi = m.atomi(i);

      atomic_number_t z = atomi->atomic_number();

      if (1 == z)
      {
        atom_type[i] = -1;
        continue;
      }

      int acon = atomi->ncon();
      int nbonds = atomi->nbonds();

      if (6 == z)
        {
          if (nbonds == acon+1)
            atom_type[i] = 9;
          else if (nbonds == acon +2)
            atom_type[i] = 18;
          else
            atom_type[i] = 0;
        }
      else if (7 == z)
        {
          if (nbonds == acon+1)
            atom_type[i] = 10;
          else if (nbonds == acon +2)
            atom_type[i] = 19;
          else
            atom_type[i] = 1;
        }
      else if (8 == z)
        {
          if (nbonds == acon+1)
            atom_type[i] = 11;
          else if (nbonds == acon +2)
            atom_type[i] = 20;
          else
            atom_type[i] = 2;
        }

      else if (9 == z)
        atom_type[i] = 3;
      else if (15 == z)
        atom_type[i] = 4;

      else if (16 == z)
        atom_type[i] = 5;

      else if (17 == z)
        atom_type[i] = 6;
      else if (35 == z)
      {
        if (compress_heavy_halogens)
          atom_type[i] = 6;
        else
          atom_type[i] = 7;
      }
      else if (53 == z)
      {
        if (compress_heavy_halogens)
          atom_type[i] = 6;
        else
          atom_type[i] = 8;
      }
      else
        atom_type[i] = -1;
    }
  return 1;
}

static int bit_number (int type_i, int type_j, int bin)
{
  int min = type_i;
  int max = type_j;
  if (min > max)
    {
      min = type_j;
      max = type_i;
    }

  int pair_number = min * DIP_NATMTYPE - min * (min-1) / 2 + max - min;

  return (pair_number * (number_of_regions + 1) + bin);
}

static int
DiP_procedure (Molecule &m, const int atom_type [], int descriptors [])
{
  int n_atoms = m.natoms();
  for (int i=0; i<n_atoms; i++)
    {
      if (-1 == atom_type[i])
        continue;
      for (int j=i+1; j<n_atoms; j++)
        {
          if (-1 == atom_type[j])
            continue;

          distance_t d = m.distance_between_atoms(i, j);
          int bin = (int) (d / region_width);
          if (bin > number_of_regions)
            bin = number_of_regions;

          int type_i = atom_type[i] % 9;
          int type_j = atom_type[j] % 9;

          // first, set the descriptor for type_i type_j
          descriptors[bit_number(type_i, type_j, bin)] ++;

          // second, if type_i != 0, set the descriptor for 0, type_j
          if (0 != type_i)
            descriptors[bit_number(0, type_j, bin)] ++;

          // third, if type_j != 0, set the descriptor for type_i, 0
          if (0 != type_j)
            descriptors[bit_number(type_i, 0, bin)] ++;

          // fourth, if type_i !=0 and type_j !=0, set the descriptor for 0, 0
          if ((0 != type_j) && (0 != type_i))
            descriptors[bit_number(0, 0, bin)] ++;

          if (atom_type[i] > 8)
            {
              int type_ii = 8 + atom_type[i] / 9;
              // first, set the descriptor for type_i, type_j
              descriptors[bit_number(type_ii, type_j, bin)] ++;
              // second, if type_j !=0, set the descriptor for type_i, 0
              if (0 != type_j)
                descriptors[bit_number(type_ii, 0, bin)] ++;
            }
          if (atom_type[j] > 8)
            {
              int type_jj = 8 + atom_type[j] / 9;
              // set the descriptor for type_i, type_j
              descriptors[bit_number(type_i, type_jj, bin)] ++;
              // second, if type_i !=0, set the descriptor for 0, type_j
              if (0 != type_i)
                descriptors[bit_number(0, type_jj, bin)] ++;
            }
          if ((atom_type[i] > 8) && (atom_type[j] > 8))
            {
              // set the descriptor for type_i, type_j
              int type_ii = 8 + atom_type[i] / 9;
              int type_jj = 8 + atom_type[j] / 9;
              descriptors[bit_number(type_ii, type_jj, bin)] ++;
            }
        }
    }
  return 1;
}

static int
compute_DiP_descriptor_for_molecule (Molecule &m, int atom_type[], int descriptors [])
{
  (void) find_DiP_atom_type_for_atoms(m, atom_type);

  return DiP_procedure(m, atom_type, descriptors);
}

/* 
 * This function is to calculate the DiP descriptors for molecule
 */

static int
compute_DiP_descriptor_for_molecule (Molecule &m, int * descriptors)
{
  int n_atoms = m.natoms();

  int * DiP_atom_type = new int [n_atoms]; std::unique_ptr<int[]> free_DiP_atom_type(DiP_atom_type);

  return compute_DiP_descriptor_for_molecule(m, DiP_atom_type, descriptors);
}


static int
jwDiP (Molecule &m, IWString_and_File_Descriptor & output)
{
// if we will need the unique smiles to store in the database, generate it now before the charge assigner and such.
// Since we don't know what state the input molecule is in, make a copy and apply chemical standardisation

  if (1==m.natoms())
  {
    cerr<<"Only one atom for this compounds, computation skipped"<<endl;
    return 1;
  }

  (void) preprocess_molecule(m);

  int * descriptors = new_int(total_number_of_descriptors); std::unique_ptr<int[]> free_descriptors(descriptors);
  
  compute_DiP_descriptor_for_molecule(m, descriptors);

  if (m.name().contains(' '))
  {
    const_IWSubstring tmp(m.name());
    tmp.truncate_at_first(' ');
    output << tmp;
  }
  else
    output << m.name();

  for (int i = 0; i < total_number_of_descriptors; i++)
  {
    iwdigits.append_number(output, descriptors[i]);
  }

  output << '\n';

  return 1;
}

static int
jwDiP (data_source_and_type<Molecule> & input, 
       IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1)
      cerr << molecules_read << " processing '" << m->name() << "'\n";

    if (! jwDiP(*m, output))
    {
      number_of_error++;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
jwDiP (const char * fname, FileType input_type,
       IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jwDiP(input, output);
}

static int
jwDiP (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:i:r:w:y");
  
  if (cl.unrecognised_options_encountered())
    {
      cerr << "Unrecognised options encountered\n";
      usage(1);
    }

  verbose = cl.option_count('v');
  
  if (! process_elements(cl))
    {
      usage(2);
    }
  
  if (! process_standard_aromaticity_options(cl, verbose))
    {
      cerr << "Cannot process aromaticity options (-A)\n";
      usage(5);
    }
  
  if (cl.option_present('r')) 
    {
      if (! cl.value('r', number_of_regions))
        {
          cerr << "Invalid value for number_of_regions (-r).  Has to be a integer between 1 and 20\n";
          usage(2);
        }
      
      if ((number_of_regions<1) || (number_of_regions>20))
        {
          cerr << "The value for number of regions is not valid, should be between 1 and 20 \n";
          cerr << "10 region is used instead\n";
          number_of_regions = 10;
        }
    }

  total_number_of_descriptors = DIP_NATMTYPE * (DIP_NATMTYPE + 1) / 2 * (number_of_regions + 1);

  if (write_results_to_output_file)
    iwdigits.initialise(total_number_of_descriptors + 1);

  if (cl.option_present('w')) 
  {
    if (! cl.value('w', region_width))
    {
      cerr << "Invalid value for region width (-w).  Has to be a value between 0.2 and 2.0\n";
      usage(2);
    }
      
    if ((region_width<0.2) || (region_width>2.0))
    {
      cerr << "The value for region width is not valid, should be between 0.2 and 2.0 \n";
      cerr << "1.0 region is used instead\n";
      region_width = 1.0;
    }
  }

  if (cl.option_present('y'))
  {
    compress_heavy_halogens = 1;

    if (verbose)
      cerr << "Will compress heavy halogens\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
    {
      if (! process_input_type(cl, input_type))
        {
          cerr << "Cannot determine input type\n";
          usage(6);
        }
    }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;
  
  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }
  
  iwdigits.set_include_leading_space(1);

  IWString_and_File_Descriptor output(1);
  output.resize(36000);

  if (! output_result_header(output))
    return 1;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! jwDiP(cl[i], input_type, output))
        {
          rc = i + 1;
          //      break;
        }
    }

    output.flush();

  // output info
  if (verbose)
  {
    if (verbose)
      cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = jwDiP(argc, argv);
  return rc;
}
