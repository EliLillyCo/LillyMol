/*
 * Calculate charge & surface area related descriptors
 * Part I   -- Jurs descriptors (total of 30)  -- Stanton & Jurs
 * Analytical Chemistry (1990) Vol 62, 2323-2329
 *
 * Part II  -- Electronic descriptors (total of 19) 
 * Chemical Review (1996) Vol 96, 1027-1043
 *
 * Part III -- Extended Jurs descriptors (total of 11)  -- Stanton et al.
 * JCICS (1992) Vol 32, 306-316
 *
 * Part IV -- Savol like descriptors (total of 6)
 * surface area is added here 
 *
 * Part V -- Aromaticity related descriptors (total of 2) 
 *
 * Partial charge calculation is based on Abraham's method
 * The last paper for partial charge is Journal of Computational Chemistry (1992) Vol 13 492-504
 *
 * Part VI -- COMMA descriptors (total of 24)  -- Silverman and Platt (exptended from original 13)
 *         -- COMMA2 descriptors (total of 28) -- Silverman four set of 7 descriptor (remove 1 --logP)
 * J. Med. Chem. (1996) 39, 2129 -2140 (COMMA)
 * JCICS (2000) 40, 1470-1476 Article (COMMA2)
 * J. Computational Chemistry 17 (3) 358-366 (1996)
 * In this case partial charge computation is based Gasteiger method
 */

#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "comma_related_procedures.h"
#include "jw_vector_array_operation.h"
#include "read_control_file.h"
#include "surface_area_molvol.h"

using std::cerr;
using std::endl;

#define NUMBER_OF_COMMA_DESCRIPTOR 53
#define NUMBER_OF_HYDROPHOBIC_JURS 28

//DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

#define NATMTYPE 10 // type of atom

/*
 * By default, the isotopic labels are just the line number in the file.
 * If the control file contains 'ISO=nnn' we set that for the isotopic label
 *
 * Each atom type can have a confidence level.
 * These range from 0 (just a guess) to 100 (unambiguous)
 *
 * We also offer the option of NOT marking the atoms as processed. this
 * allows the possibility of multiple queries hitting the same atoms
 *
 * We allow a variable number of atoms to match per query. 
 *
 * By default, we are an atom additivity method, and so just process
 * the first atom in any embedding
 */

int max_atoms_per_embedding_to_process = 1;

resizable_array<int> atoms_to_match;
resizable_array<int> mark_atoms;
resizable_array<int> label;
resizable_array<int> confidence;
resizable_array_p<Substructure_Hit_Statistics> queries;

static Donor_Acceptor_Assigner donor_acceptor_assigner;
static Charge_Assigner charge_assigner;

//  When writing to a database we need chemical standardisation

static Chemical_Standardisation chemical_standardisation;

int verbose=0;

static int molecules_read = 0;

static int use_gasteiger_partial_charge = 0;

static int compute_jurs_descriptor = 0;
static int compute_hydrophobic_jurs_descriptor = 0;
static int compute_electronic_descriptor = 0;
static int compute_savol_like_descriptor = 0;
static int compute_aromaticity_related_descriptor = 0;
static int compute_comma_related_descriptor = 0;

// allow different probe to be used here
static double probe_radius = 1.4;
static int vdw_radius_type = IW_VDW_MOLVOL;

// keep track of number of warning and number of error
static int number_of_error = 0;

static double min_value_comma_descriptor = -2000.0;
static double max_value_comma_descriptor =  2000.0;

static IWString tag;
static IWString smiles_tag("$SMI<");

static double bit_count_scaling_factor = 10.0;

/*
  Run this against the Lilly database and get the 1 and 99'th percentiles
*/

static double min_jurs [] = {
317.9,   // jurs_JWJSASA
112.1,   // jurs_JWJPPSA1
23.55,   // jurs_JWJPNSA1
-224.9,   // jurs_JWJDPSA1
53.26,   // jurs_JWJPPSA2
-939.4,   // jurs_JWJPNSA2
92.97,   // jurs_JWJDPSA2
3.48,   // jurs_JWJPPSA3
-64.83,   // jurs_JWJPNSA3
10.71,   // jurs_JWJDPSA3
0.266,   // jurs_JWJFPSA1
0.04376,   // jurs_JWJFNSA1
0.1457,   // jurs_JWJFPSA2
-1.23,   // jurs_JWJFNSA2
0.007342,   // jurs_JWJFPSA3
-0.1038,   // jurs_JWJFNSA3
41.5,   // jurs_JWJWPSA1
11.26,   // jurs_JWJWNSA1
18.4,   // jurs_JWJWPSA2
-789.5,   // jurs_JWJWNSA2
1.435,   // jurs_JWJWPSA3
-53.63,   // jurs_JWJWNSA3
0.0886,   // jurs_JWJRPCG
0.1065,   // jurs_JWJRNCG
0,   // jurs_JWJRPCS
0,   // jurs_JWJRNCS
0,   // jurs_JWJTPSA
283.1,   // jurs_JWJTASA
0,   // jurs_JWJRPSA
0.728,   // jurs_JWJRASA
0,   // jurs_JWJ2SSAH
0,   // jurs_JWJ2CHGD
0,   // jurs_JWJ2ACGD
0,   // jurs_JWJ2SSAA
0,   // jurs_JWJ2CNTH
0,   // jurs_JWJ2CNTA
1,   // jurs_JWJ2RHTA
0,   // jurs_JWJ2RSAH
0,   // jurs_JWJ2RSAA
0,   // jurs_JWJ2RSHM
0,   // jurs_JWJ2RSAM
0.4831,   // jurs_JWELED1
0.02885,   // jurs_JWELED2
0.2415,   // jurs_JWELED3
0.0192,   // jurs_JWELED4
-0.2477,   // jurs_JWELED5
0.07609,   // jurs_JWELED6
0.004256,   // jurs_JWELED7
-0.3744,   // jurs_JWELED8
-0.1015,   // jurs_JWELED9
0.1266,   // jurs_JWELED10
0.005383,   // jurs_JWELED11
0.01367,   // jurs_JWELED12
0.001033,   // jurs_JWELED13
0.02898,   // jurs_JWELED14
0.005551,   // jurs_JWELED15
0.002641,   // jurs_JWELED16
0.1944,   // jurs_JWELEDP1
0.0378,   // jurs_JWELEDP2
0.7071,   // jurs_JWELEDP3
-34.82,   // jurs_JWSSovEn
475.5,   // jurs_JWSVol
11.65,   // jurs_JWSPSA12
4.213,   // jurs_JWSPSA22
2.106,   // jurs_JWSPSA32
138.5,   // jurs_JWSHpbSA
0,   // jurs_JWAromSA
0,   // jurs_JWAromRa
116,   // jurs_JWPPHS1
12.37,   // jurs_JWPNHS1
87.81,   // jurs_JWPPHS2
-3106,   // jurs_JWPNHS2
17.77,   // jurs_JWPPHS3
-195.4,   // jurs_JWPNHS3
-254.3,   // jurs_JWDHS1
365.5,   // jurs_JWDHS2
55.89,   // jurs_JWDHS3
0.2637,   // jurs_JWFPHS1
0.02361,   // jurs_JWFNHS1
0.2337,   // jurs_JWFPHS2
-3.998,   // jurs_JWFNHS2
0.04222,   // jurs_JWFPHS3
-0.3376,   // jurs_JWFNHS3
42.97,   // jurs_JWWPHS1
5.973,   // jurs_JWWNHS1
31.19,   // jurs_JWWPHS2
-2709,   // jurs_JWWNHS2
6.582,   // jurs_JWWPHS3
-146.5,   // jurs_JWWNHS3
0.05284,   // jurs_JWRPH
0.09814,   // jurs_JWRNH
0.05038,   // jurs_JWRPHS
0,   // jurs_JWRNHS
0.02279,   // jurs_JWSURR1
-9.193,   // jurs_JWSURR2
-5.979   // jurs_JWSURR3
};

static double max_jurs [] = {
1000,   // jurs_JWJSASA
806.1,   // jurs_JWJPPSA1
435.7,   // jurs_JWJPNSA1
679,   // jurs_JWJDPSA1
1647,   // jurs_JWJPPSA2
-8.772,   // jurs_JWJPNSA2
2303,   // jurs_JWJDPSA2
53.29,   // jurs_JWJPPSA3
-2.581,   // jurs_JWJPNSA3
97.81,   // jurs_JWJDPSA3
0.9562,   // jurs_JWJFPSA1
0.734,   // jurs_JWJFNSA1
1.764,   // jurs_JWJFPSA2
-0.01894,   // jurs_JWJFNSA2
0.07893,   // jurs_JWJFPSA3
-0.00486,   // jurs_JWJFNSA3
785,   // jurs_JWJWPSA1
356.2,   // jurs_JWJWNSA1
1595,   // jurs_JWJWPSA2
-3.795,   // jurs_JWJWNSA2
45.78,   // jurs_JWJWPSA3
-1.263,   // jurs_JWJWNSA3
0.539,   // jurs_JWJRPCG
0.7521,   // jurs_JWJRNCG
17.26,   // jurs_JWJRPCS
26.88,   // jurs_JWJRNCS
161.9,   // jurs_JWJTPSA
945.8,   // jurs_JWJTASA
0.2719,   // jurs_JWJRPSA
1,   // jurs_JWJRASA
0,   // jurs_JWJ2SSAH
0,   // jurs_JWJ2CHGD
0,   // jurs_JWJ2ACGD
0,   // jurs_JWJ2SSAA
0,   // jurs_JWJ2CNTH
0,   // jurs_JWJ2CNTA
1,   // jurs_JWJ2RHTA
0,   // jurs_JWJ2RSAH
0,   // jurs_JWJ2RSAA
0,   // jurs_JWJ2RSHM
0,   // jurs_JWJ2RSAM
5.027,   // jurs_JWELED1
0.1618,   // jurs_JWELED2
2.514,   // jurs_JWELED3
0.1506,   // jurs_JWELED4
-0.05017,   // jurs_JWELED5
0.5086,   // jurs_JWELED6
0.04822,   // jurs_JWELED7
-0.1252,   // jurs_JWELED8
-0.01447,   // jurs_JWELED9
0.5086,   // jurs_JWELED10
0.02957,   // jurs_JWELED11
0.5296,   // jurs_JWELED12
0.03675,   // jurs_JWELED13
0.5553,   // jurs_JWELED14
0.06524,   // jurs_JWELED15
0.03599,   // jurs_JWELED16
4.956,   // jurs_JWELEDP1
24.56,   // jurs_JWELEDP2
10.86,   // jurs_JWELEDP3
6.387,   // jurs_JWSSovEn
1766,   // jurs_JWSVol
445,   // jurs_JWSPSA12
343.6,   // jurs_JWSPSA22
323.4,   // jurs_JWSPSA32
849.2,   // jurs_JWSHpbSA
555.6,   // jurs_JWAromSA
0.8011,   // jurs_JWAromRa
761,   // jurs_JWPPHS1
494.5,   // jurs_JWPNHS1
5021,   // jurs_JWPPHS2
-6.404,   // jurs_JWPNHS2
226.9,   // jurs_JWPPHS3
-2.798,   // jurs_JWPNHS3
648.8,   // jurs_JWDHS1
6157,   // jurs_JWDHS2
300.7,   // jurs_JWDHS3
0.9746,   // jurs_JWFPHS1
0.7361,   // jurs_JWFNHS1
5.976,   // jurs_JWFPHS2
-0.01383,   // jurs_JWFNHS2
0.3815,   // jurs_JWFPHS3
-0.0056,   // jurs_JWFNHS3
717,   // jurs_JWWPHS1
434.6,   // jurs_JWWNHS1
4599,   // jurs_JWWPHS2
-2.891,   // jurs_JWWNHS2
176.1,   // jurs_JWWPHS3
-1.374,   // jurs_JWWNHS3
0.5,   // jurs_JWRPH
1,   // jurs_JWRNH
33.63,   // jurs_JWRPHS
34.59,   // jurs_JWRNHS
2.723,   // jurs_JWSURR1
-0.003545,   // jurs_JWSURR2
-0.02291    // jurs_JWSURR3
};

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Calculate charge surface area related descriptors for molecules\n";

  cerr << "  -j             compute jurs and extended jurs descriptors\n";
  cerr << "  -h             compute hydrophobic jurs descriptors\n";
  cerr << "  -e             compute electronic property related descriptors\n";
  cerr << "  -s             compute savol like descriptors\n";
  cerr << "  -a             compute aromaticity related descriptors\n";
  cerr << "  -b             compute CoMMA related descriptors -- -F option has to be turned on\n";

  cerr << "                 by default, if none of -j -e -s -a -b is specified, all flags are turned on\n\n";
  cerr << "  -r <radius>    set solvent probe radius (default = 1.40 for water)\n";
  cerr << "  -J <tag>       tag for fingerprint output\n";
  cerr << "  -y <scale>     bit count scaling factor (default 10.0)\n";

  cerr << "  -i <type>      input type\n";

  cerr << "  -H <..>        donor acceptor assignment, enter '-H help' for info\n";

  cerr << "  -F file        file in which ghose-crippen parameters are saved\n";

  cerr << "  -G             use gasteiger partial charge in computation\n";
  cerr << "  -M min=<m>     specify minimum value for comma descriptors\n";
  cerr << "  -M max=<m>     specify maximum value for comma descriptors\n";

  (void) display_standard_vdw_radius_types (cerr, 'V');

  (void) display_standard_charge_assigner_options (cerr, 'N');

  (void) display_standard_aromaticity_options (cerr);
  
  (void) display_standard_chemical_standardisation_options (cerr, 'g');

#ifdef USE_IWMALLOC
  display_standard_debug_options (cerr, 'd');
#endif
  cerr << "  -v             verbose output\n";
  
  exit (rc);
}

static void
hit_bit (Sparse_Fingerprint_Creator & sfc,
         int b,
         double v,
         double minval, double maxval,
         double scale)
{
#ifdef CHECK_OUT_OF_RANGE
  if (v < minval || v > maxval)
    cerr << "Value out of range " << v << " min " << minval << " max " << maxval << " feature " << b << endl;
#endif

  if (v <= minval)
  {
    sfc.hit_bit(b, 1);
    return;
  }

  if (v >= maxval)
    v = maxval;

  int c = static_cast<int> ((v - minval) / (maxval - minval) * scale) + 1;

  sfc.hit_bit(b, c);

  return;
}

static int
ghose_crippin_procedure (Molecule & m, int * already_hit, 
			 double * hydrophobicity_index)
{
  Molecule_to_Match target (&m);
  
  int nq = queries.number_elements();
  
  for (int i = 0; i < nq; i++)
    {
      Substructure_Results sresults;
      
      int nhits = queries[i]->substructure_search (target, sresults);
      
      if (0 == nhits)
	continue;
      
      double nvi;
      (void) queries[i]->numeric_value (nvi);
      
      for (int j = 0; j < nhits; j++)
	{
	  Set_of_Atoms *e =  const_cast<Set_of_Atoms *> (sresults.embedding (j));   // loss of const OK
	  
	  int hit_atom = e->item(0);
	  if (already_hit[hit_atom] == 0)
	    {
	      already_hit[hit_atom] = 1;
	      hydrophobicity_index[hit_atom]= nvi;
	    }
	}
    }
  
  int atoms_not_hit = 0;
  int matoms = m.natoms();
  for (int i=0; i<matoms; i++)
  {
    if (already_hit[i] == 0)
      atoms_not_hit ++;
  }
  
  if (atoms_not_hit>1) 
    return 0;
  else return 1;
}

static int
ghose_crippen_procedure (Molecule &m, double hydrophobicity_index[])
{
  int n_atoms = m.natoms();

  int * already_hit = new_int(n_atoms); std::unique_ptr<int[]> free_already_hit(already_hit);

  set_vector(hydrophobicity_index, n_atoms, 0.0);

  int rc = ghose_crippin_procedure (m, already_hit, hydrophobicity_index);

  return rc;

}

static int
preprocess_molecule (Molecule & m)
{
  m.reduce_to_largest_fragment();  // always reduce to largest fragment
  
  return 1;
}

// compute three different definition of polar surface area
// method 1:  anything but carbon is considered polar
static int
compute_polar_surface_area_1 (Molecule &m, area_t atom_area[], double & polar_surface_area)
{
  int n_atom = m.natoms();
  polar_surface_area = 0.0;
  for (int i= 0; i<n_atom;i++)
    {
      if (m.atomi(i)->atomic_number() != 6)
	polar_surface_area +=atom_area[i];
    }
  return 1;
}

// method 2:  O & N (execpt in nitro -N(=O)(=O)) S (except in -S(=O)(=O)-) 
static int
compute_polar_surface_area_2 (Molecule &m, area_t atom_area[], double &polar_surface_area)
{
  int n_atom = m.natoms();
  polar_surface_area = 0.0;
  for (int i= 0; i<n_atom;i++)
    {
      const Atom * atom_i=m.atomi(i);
      int icon = atom_i ->ncon();
      if (icon>=4)
	continue; // anything has connection large than 4 -- cannot be seen
      else if ((7==atom_i->atomic_number()) && (3==icon) && (5 == atom_i->nbonds()))
	continue; // nitro group is taken care of here
      else if ((8 == atom_i->atomic_number()) || (7 == atom_i->atomic_number()) || (16 == atom_i->atomic_number())) 
	polar_surface_area +=atom_area[i];
      else continue;
    }
  return 1;
}

// method 3:  O, N and H linked to O N.
static int
compute_polar_surface_area_3 (Molecule &m, area_t atom_area[], double &polar_surface_area)
{
  int n_atom = m.natoms();
  polar_surface_area = 0;
  for (int i= 0; i<n_atom;i++)
    {
      const Atom * atom_i=m.atomi(i);
      if ((atom_i->atomic_number()==7) || (atom_i->atomic_number()==8) || ((atom_i->atomic_number()==1)&& ((m.atomi(m.other(i,0))->atomic_number()==7)||(m.atomi(m.other(i,0))->atomic_number()==8))))
	{
	  polar_surface_area +=atom_area[i];
	}
    }

  return 1;
}

// compute hydrophobic surface area
static int
compute_hydrophobic_surface_area (Molecule &m,
                                  area_t atom_area[],
                                  double &hydrophobic_surface_area)
{
  int n_atom = m.natoms();
  hydrophobic_surface_area = 0.0;
  for (int i = 0; i < n_atom; i++)
  {
    const Atom * a = m.atomi (i);

    atomic_number_t z = a->atomic_number();

    if (7 == z)
      continue;
    if (8 == z)
      continue;
    if (16 == z)
      continue;

    //  Note the strangeness of a pyridine here. Only one of the neighbouring carbons has a
    //  multiple bond to a heteroatom. Let's just skip over that complication.

    if (6 == z && m.multiple_bond_to_heteroatom(i))
      continue;
    
    int icon = a->ncon();
    
    if (icon >= 4)     // will be invisible - the C in CF3 for example
      continue;
    
    if (icon > 1 && (17 == z || 35 == z || 53 == z))    // only singly connected halogens are hydrophillic
      continue;
    
    hydrophobic_surface_area += atom_area[i];
  }
  return 1;
}

// solvation parameter
static double loncharich_carbon = 12.3 * 1.0e-03;
static double loncharich_nitrogen = -115.7 * 1.0e-03;
static double loncharich_oxygen   = -115.7 * 1.0e-03;
static double loncharich_sulphur  = -18.3 * 1.0e-03;
static double loncharich_ominus   = -175.4 * 1.0e-03;
static double loncharich_nplus    = -185.5 * 1.0e-03;
static double loncharich_fluorine = 21.0 * 1.0e-03;
static double loncharich_chlorine = -1.5 * 1.0e-03;
static double loncharich_bromine  = -5.2 * 1.0e-03;
static double loncharich_iodine   = -9.3 * 1.0e-03;

// compute savol like solvation energy
static int
compute_solvation_energy (Molecule &m, area_t atom_area[], double &energy)
{
  // calculate the solvation energy, result is saved in energy
  int n_atom = m.natoms();
  energy = 0.0;
  for (int i =0; i<n_atom; i++)
    {
      atomic_number_t z = m.atomic_number (i);

      // for the hydrogens, give its surface area to the bonded atom;
      if (1==z)
	z = m.atomic_number (m.other (i, 0));

      if (1==z)
	;
      if (53 == z)
	energy += atom_area[i] * loncharich_iodine;
      else if (35 == z)
	energy += atom_area[i] * loncharich_bromine;
      else if (17 == z)
	energy += atom_area[i] * loncharich_chlorine;
      else if (9 == z)
	energy += atom_area[i] * loncharich_fluorine;
      else if (7 == z && 1 == m.formal_charge (i))
	energy += atom_area[i] * loncharich_nplus;
      else if (16 == z)
	energy += atom_area[i] * loncharich_sulphur;
      else if (8 == z && -1 == m.formal_charge (i))
	energy += atom_area[i] * loncharich_ominus;
      else if (8 == z)
	energy += atom_area[i] * loncharich_oxygen;
      else if (7 == z)
	energy += atom_area[i] * loncharich_nitrogen;
      else if (6 == z)
	energy += atom_area[i] * loncharich_carbon;
      else if (15 == z) // might want to do something else here
	{;}
      else return 0;
    }
  return 1;
}

static int
molecule_surface_area (Molecule &m,
                       area_t atom_area[],
                       volume_t &total_molecule_volume)
{
  int rc;  // return status
  
  area_t total_molecule_area=0.0;
  
  Surface_Area_Molvol molecule_surface; 
  
  molecule_surface.set_vdw_type (vdw_radius_type);
  molecule_surface.set_probe_radius (probe_radius); 

  //  cout <<m.name()<<" ";
  Molecule mcopy = m;
  rc=molecule_surface.surface_area(mcopy, atom_area, total_molecule_area, total_molecule_volume);
  
  //  cout<<endl;
  if (0==rc) return 0;
  // checking to see if calculated atom surface area is good
#if (DEBUG_SWITCH) 
  cerr<<"Inside function molecule_surface_area"<<endl;
  cerr<<"Atom coordinates"<<endl<<endl;
  
  int n_atom=m.natoms();
  for (int i=0;i<n_atom;i++)
    {
      cerr<<"Atom number "<<i<<"   ";
      cerr<<m.x(i)<<"    "<<m.y(i)<<"    "<<m.z(i)<<endl;
    }
  
  cerr<<"Inside calculate_partial_charge, just exited from surface_area() of surface_area_molvol.cc"<<endl;
  for (int i=0;i<n_atom;i++)
    {
      // printout the surface_area for each atom of the molecule
      cerr<<atom_area[i]<< "      ";
    }
  cerr<<endl<<endl;
#endif
  
  return 1;
}

static int
output_result_header (IWString_and_File_Descriptor & output)
{
  output << "Name";

  if (compute_jurs_descriptor)
    { 
      // print out header for jurs descriptor
      output<<" jurs_JWJSASA jurs_JWJPPSA1 jurs_JWJPNSA1 jurs_JWJDPSA1 jurs_JWJPPSA2 jurs_JWJPNSA2 jurs_JWJDPSA2 jurs_JWJPPSA3 jurs_JWJPNSA3 jurs_JWJDPSA3";
      output<<" jurs_JWJFPSA1 jurs_JWJFNSA1 jurs_JWJFPSA2 jurs_JWJFNSA2 jurs_JWJFPSA3 jurs_JWJFNSA3";
      output<<" jurs_JWJWPSA1 jurs_JWJWNSA1 jurs_JWJWPSA2 jurs_JWJWNSA2 jurs_JWJWPSA3 jurs_JWJWNSA3";
      output<<" jurs_JWJRPCG jurs_JWJRNCG jurs_JWJRPCS jurs_JWJRNCS jurs_JWJTPSA jurs_JWJTASA jurs_JWJRPSA jurs_JWJRASA";


      // print out header for hydrogen bond extension of jurs descriptor
      output <<" jurs_JWJ2SSAH jurs_JWJ2CHGD jurs_JWJ2ACGD jurs_JWJ2SSAA jurs_JWJ2CNTH jurs_JWJ2CNTA";
      output <<" jurs_JWJ2RHTA jurs_JWJ2RSAH jurs_JWJ2RSAA jurs_JWJ2RSHM jurs_JWJ2RSAM";
    }

  if (compute_electronic_descriptor)
    {
      // print out header for electronic descriptor
      output<<" jurs_JWELED1 jurs_JWELED2 jurs_JWELED3 jurs_JWELED4 jurs_JWELED5 jurs_JWELED6 jurs_JWELED7 jurs_JWELED8 jurs_JWELED9";
      output<<" jurs_JWELED10 jurs_JWELED11 jurs_JWELED12 jurs_JWELED13 jurs_JWELED14 jurs_JWELED15 jurs_JWELED16";
      output<<" jurs_JWELEDP1 jurs_JWELEDP2 jurs_JWELEDP3";
    }

  if (compute_savol_like_descriptor)
    {
      // commented out long version of the descriptor name
      output<<" jurs_JWSSovEn jurs_JWSVol jurs_JWSPSA12 jurs_JWSPSA22 jurs_JWSPSA32 jurs_JWSHpbSA"; 
    }

  if (compute_aromaticity_related_descriptor)
    {
      output<<" jurs_JWAromSA jurs_JWAromRa";
    }

  if (compute_hydrophobic_jurs_descriptor)
    {
      output <<" jurs_JWPPHS1 jurs_JWPNHS1 jurs_JWPPHS2 jurs_JWPNHS2 jurs_JWPPHS3 jurs_JWPNHS3 jurs_JWDHS1 jurs_JWDHS2 jurs_JWDHS3";
      output <<" jurs_JWFPHS1 jurs_JWFNHS1 jurs_JWFPHS2 jurs_JWFNHS2 jurs_JWFPHS3 jurs_JWFNHS3";
      output <<" jurs_JWWPHS1 jurs_JWWNHS1 jurs_JWWPHS2 jurs_JWWNHS2 jurs_JWWPHS3 jurs_JWWNHS3";
      output <<" jurs_JWRPH jurs_JWRNH jurs_JWRPHS jurs_JWRNHS jurs_JWSURR1 jurs_JWSURR2 jurs_JWSURR3";
    }

  if (compute_comma_related_descriptor)
    {
      for (int i =0; i<NUMBER_OF_COMMA_DESCRIPTOR; i++)
	output<<" comma_JWCoMM"<<i;
    }

  output << '\n';

  return 1;
}

/* 
 * compute the relationship between the center of dipole and center of mass
 */
static void 
compute_electronic_property_comma_descriptor (double dipole_center[], double dipole_moment[], 
					      double quadrupole_moment[], double center[],
					      double ** reference_axis, double descriptor[])
{
  double axis[3];
  double displacement[3];

  for (int i=0; i<3; i++)
    {
      displacement[i] = dipole_center[i] - center[i];
      //      cerr<<"displacement["<<i<<"]="<<displacement[i]<<endl;
    }

  for (int i=0; i<3; i++)
    {
      //      cerr<<"axis i="<<i<<" ";

      for (int j=0; j<3; j++)
	{
	  axis[j] = reference_axis[j][i];
	  //	  cerr<<axis[j]<<" ";
	}

      //      cerr<<endl;

      descriptor[i] = fabs (dot (displacement, axis));
      //      cerr<<"descriptor="<<descriptor[i]<<endl;
      descriptor[3+i] = fabs (dot (dipole_moment, axis));
    }

  for (int i=0; i<2; i++)
    descriptor[6+i] = quadrupole_moment[i];

}

/* 
 * compute the quadrupole moment (principal quadrupole moment and Qxx, Qyy)
 */

static double 
compute_quadrupole_moment (Molecule &m,
                           double ** quadrupole_moment_tensor,
                           double quadrupole_moment[])
{
  double ** quadrupole_moment_axis = new double *[3];
  for (int i=0; i<3; i++)
    quadrupole_moment_axis[i] = new double[3];

//double quadrupole_moment_axis[9];

  (void) m.find_eigen_value_for_matrix (3, quadrupole_moment_tensor, quadrupole_moment, quadrupole_moment_axis);

  for (int i=0; i<3; i++)
    delete[] quadrupole_moment_axis[i];

  delete[] quadrupole_moment_axis;

  return quadrupole_moment[0];
}

/*
 * compute the center of the dipole and the corresponding quadrupole moment tensor
 */

static int 
compute_center_of_dipole (double dipole_moment[], double ** quadrupole_moment_tensor,
			  double dipole_center[])
{
  double tmp_Qp[3];
  
  // initialize the Q dot p vector
  for (int i=0; i<3; i++)
    tmp_Qp[i] = 0;

  for (int i=0; i<3; i++)
    dipole_center[i] = 0;

  // compute Q dot p vector
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      tmp_Qp[i] += quadrupole_moment_tensor[i][j] * dipole_moment[j];
 
  double tmp =0;
  double length_square =0;

  for (int i=0; i<3; i++)
    {
      tmp += dipole_moment[i] * tmp_Qp[i];
      length_square += dipole_moment[i] * dipole_moment[i];
    }

  if (length_square >1e-3)
    {
      tmp /= 4 * length_square;

      for (int i=0; i<3; i++)
	dipole_center[i] = (tmp_Qp[i] - tmp * dipole_moment[i]) * 2 / 3 / length_square;
      return 1;
    }
  return 0;
}

static void
compute_Qxx_Qyy_at_dipole_center_using_inertial_axis (int n_atoms, double atomic_partial_charge[], double **coordinate, 
						      double origin[], double ** axis, 
						      double quadrupole_moment[])
{
  for (int i=0; i<3; i++)
    quadrupole_moment[i] = 0;
  
  for (int i=0; i<n_atoms; i++)
    {
      double tmp_after_translation[3];
      double tmp_after_rotation[3];
      
      for (int j=0; j<3; j++)
	tmp_after_translation[j] = (coordinate[i][j] - origin[j]) ;

      for (int j=0; j<3; j++)
	{
	  tmp_after_rotation[j] = 0;
	  for (int k=0; k<3; k++) 
	    tmp_after_rotation[j] += tmp_after_translation[k] * axis[k][j];
	}

      double xx = tmp_after_rotation[0] * tmp_after_rotation[0];
      double yy = tmp_after_rotation[1] * tmp_after_rotation[1];
      double zz = tmp_after_rotation[2] * tmp_after_rotation[2];

      quadrupole_moment[0] +=  (2 * xx - yy - zz) * atomic_partial_charge[i];
      quadrupole_moment[1] += (-xx + 2 * yy - zz) * atomic_partial_charge[i];
    }
}
/*
 * Compute dipole moment and quadrupole moment tensor from the origin
 */

static void
compute_dipole_and_quadrupole_from_origin (int n_atoms, double atomic_partial_charge[], double dipole_moment[], 
					   double ** quadrupole_moment_tensor, double ** coordinate, double origin[])
{
  // initialize the dipole moment
  for (int i=0; i<3; i++)
    dipole_moment[i] =0;

  /*  // compute the dipole moment vector
  for (int i=0; i<n_atoms; i++)
    for (int j=0; j<3; j++)
      dipole_moment[j] += atomic_partial_charge[i] * coordinate[i][j];
  */

  // initialize the quadrupole_moment_tensor
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      quadrupole_moment_tensor[i][j] = 0;

  for (int i=0; i<n_atoms; i++)
    {
      // compute the origin-adjusted cooridinate
      double tmp_x = coordinate[i][0] - origin[0];
      double tmp_y = coordinate[i][1] - origin[1];
      double tmp_z = coordinate[i][2] - origin[2];

      // compute the dipole moment vector
      dipole_moment[0] += atomic_partial_charge[i] * tmp_x;
      dipole_moment[1] += atomic_partial_charge[i] * tmp_y;
      dipole_moment[2] += atomic_partial_charge[i] * tmp_z;

      // compute the quadrupole_moment_tensor
      double xx = tmp_x * tmp_x;
      double yy = tmp_y * tmp_y;
      double zz = tmp_z * tmp_z;
      double xy = tmp_x * tmp_y;
      double xz = tmp_x * tmp_z;
      double yz = tmp_y * tmp_z;
      quadrupole_moment_tensor[0][0] +=  (2 * xx - yy - zz) * atomic_partial_charge[i];
      quadrupole_moment_tensor[1][1] += (-xx + 2 * yy - zz) * atomic_partial_charge[i];
      quadrupole_moment_tensor[2][2] += (-xx - yy + 2 * zz) * atomic_partial_charge[i];
      quadrupole_moment_tensor[0][1] += 3 * xy * atomic_partial_charge[i];
      quadrupole_moment_tensor[0][2] += 3 * xz * atomic_partial_charge[i];
      quadrupole_moment_tensor[1][2] += 3 * yz * atomic_partial_charge[i];
    }
  quadrupole_moment_tensor[1][0] = quadrupole_moment_tensor[0][1];
  quadrupole_moment_tensor[2][0] = quadrupole_moment_tensor[0][2];
  quadrupole_moment_tensor[2][1] = quadrupole_moment_tensor[1][2];

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      quadrupole_moment_tensor[i][j] /=2;
}


/*
 * compute the geometry and gravity center
 */

static void
compute_geometry_and_gravity_center (Molecule &m, double ** coordinate, double geometry_center[], double gravity_center[])
{
  int n_atoms = m.natoms();
  double * geometry_field = new double[n_atoms]; std::unique_ptr<double[]> free_geometry_field (geometry_field);
  double * gravity_field = new double[n_atoms];  std::unique_ptr<double[]> free_gravity_field (gravity_field);

  for (int i=0; i<n_atoms; i++)
    {
      geometry_field[i] = 1;
      gravity_field[i] = m.atomic_mass (i);
    }

  compute_field_center (n_atoms, coordinate, geometry_center, geometry_field);
  compute_field_center (n_atoms, coordinate, gravity_center, gravity_field);

}
 
/*
 * This is the actual procedure for the computation of CoMMA descriptor
*/
static void
comma_related_descriptor_computation_procedure (Molecule &m,
                        double atom_partial_charge[],
                        double comma_descriptor[],
                        double ** coordinate, double & ratio)
{
  int n_atoms = m.natoms();

  double origin[3];
  for (int i=0; i<3; i++)
    origin[i] = 0.0;

  double geometry_center[3];
  double gravity_center[3];
  double dipole_center[3];

  double geometry_electronic_related[8];
  double gravity_electronic_related[8];

  // store the quadrupole_moment Qxx, Qyy, Qzz
  double quadrupole_moment[3];
  double principal_quadrupole_moment;

  double geometry_moment_of_inertia[3];
  double gravity_moment_of_inertia[3];

  double dipole_moment[3];

  double ** geometry_moment_of_inertia_axis = new double *[3];
  double ** gravity_moment_of_inertia_axis = new double *[3];
  double ** quadrupole_moment_tensor = new double *[3];

  for (int i=0; i<3; i++)
    {
      geometry_moment_of_inertia_axis[i] = new double[3];
      gravity_moment_of_inertia_axis[i] = new double[3];
      quadrupole_moment_tensor[i] = new double[3];
    }

  double * gravity_field = new double[n_atoms]; std::unique_ptr<double[]> free_gravity_field (gravity_field);
  double * geometry_field = new double[n_atoms]; std::unique_ptr<double[]> free_geometry_field (geometry_field);

  for (int i=0; i<n_atoms; i++)
    {
      gravity_field[i] = m.atomic_mass (i);
      geometry_field[i] = 1.0;
    }
  // compute the geometry center and the gravity center
  compute_geometry_and_gravity_center (m, coordinate, geometry_center, gravity_center);

  // compute moment of inertia for geometry and gravity around their corresponding center Ix, Iy, Iz for geometry and gravity
  compute_moment_of_inertia (m, n_atoms, coordinate, geometry_center, geometry_field, geometry_moment_of_inertia_axis, geometry_moment_of_inertia);

  compute_moment_of_inertia (m, n_atoms, coordinate, gravity_center, gravity_field, gravity_moment_of_inertia_axis, gravity_moment_of_inertia);

  for (int i=0; i<3; i++)
    comma_descriptor[2+i] = geometry_moment_of_inertia[i];

  for (int i=0; i<3; i++)
    comma_descriptor[13+i] = gravity_moment_of_inertia[i];
    
  if (geometry_moment_of_inertia[2] != 0.0)
    ratio = geometry_moment_of_inertia[0] / geometry_moment_of_inertia[2];

  // compute the dipole and quadrupole using (0, 0, 0) as the origin

  compute_dipole_and_quadrupole_from_origin (n_atoms, atom_partial_charge, dipole_moment, quadrupole_moment_tensor, coordinate, origin);

  // compute the center of dipole, 
  if (! compute_center_of_dipole (dipole_moment, quadrupole_moment_tensor, dipole_center))
    {
      if (verbose)
        cerr<<"molecule '"<<m.name()<<"' has very small dipole_moment in computation, dipole center procedure set center to (0,0,0)\n";
    }
  else
    {
      // compute the dipole and quadrupole using the center of dipole as the origin
      compute_dipole_and_quadrupole_from_origin (n_atoms, atom_partial_charge, dipole_moment, quadrupole_moment_tensor, coordinate, dipole_center);
      
      // at dipole center, the eigen value for the quadrupole_moment_tensor are p, -p, 0. p is the principal quadrupole moment
      principal_quadrupole_moment = compute_quadrupole_moment (m, quadrupole_moment_tensor, quadrupole_moment);
      
      // at the dipole center, using geometry inertia axis as the reference axis, eigen value Qxx, Qyy, Qzz Qzz = -(Qxx + Qyy)  
      compute_Qxx_Qyy_at_dipole_center_using_inertial_axis (n_atoms, atom_partial_charge, coordinate, dipole_center, geometry_moment_of_inertia_axis, quadrupole_moment);
      
      // compute the value along geometry moment of inertia axis, of dx, dy, dz (displacement), px, py, pz (dipole), and Qxx, Qyy
      compute_electronic_property_comma_descriptor (dipole_center, dipole_moment, quadrupole_moment, geometry_center, geometry_moment_of_inertia_axis, geometry_electronic_related);
      
      // at the dipole center, using gravity inertia axis as the reference axis, eigen value Qxx, Qyy, Qzz Qzz = -(Qxx + Qyy)  
      compute_Qxx_Qyy_at_dipole_center_using_inertial_axis (n_atoms, atom_partial_charge, coordinate, dipole_center, gravity_moment_of_inertia_axis, quadrupole_moment);
      
      // compute the value along gravity moment of inertia axis, of dx, dy, dz (displacement), px, py, pz (dipole), and Qxx, Qyy
      compute_electronic_property_comma_descriptor (dipole_center, dipole_moment, quadrupole_moment, gravity_center, gravity_moment_of_inertia_axis, gravity_electronic_related);
      
      // length of dipole moment
      comma_descriptor[0] = vector_length (dipole_moment);
      
      // length of principal quadrupole moment
      comma_descriptor[1] = fabs (principal_quadrupole_moment);
      
      for (int i=0; i<8; i++)
        comma_descriptor[5+i] = geometry_electronic_related[i];

      for (int i=0; i<8; i++)
        comma_descriptor[16+i] = gravity_electronic_related[i];
    }

  for (int i=0; i<3; i++)
    {
      delete[] geometry_moment_of_inertia_axis[i];
      delete[] gravity_moment_of_inertia_axis[i];
      delete[] quadrupole_moment_tensor[i];
    }
  
  delete[] geometry_moment_of_inertia_axis;
  delete[] gravity_moment_of_inertia_axis;
  delete[] quadrupole_moment_tensor;
}

/*
 * Assign field for COMMA computation
 */
struct entry 
{
  IWString TYPE;
  double field_parameter[3]; 
  //0 -- vdw volume; 
  //1 -- electronegativity;
  //2 -- polarizability;
};
  
static struct entry AtmTypeList[NATMTYPE]=
{
  { "H", {0.299,  0.944, 0.379} },
  { "C", {1.000,  1.000, 1.000} },
  { "N", {0.695,  1.163, 0.625} },
  { "O", {0.512,  1.331, 0.456} },
  { "F", {0.410,  1.457, 0.316} },
  { "P", {1.181,  0.916, 2.063} },
  { "S", {1.088,  1.077, 1.648} },
  { "Cl", {1.035,  1.265, 1.239} },
  { "Br", {1.384,  1.172, 1.733} },
  { "I", {1.728,  1.012, 3.040} }
};

static int
assign_field_for_at_comma_computation (Molecule &m,
                                double field[],
                                int field_index)
{
  
  int n_atoms = m.natoms();

  for (int i = 0; i < n_atoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    switch(z)
    {
      case 1:
        field[i] = AtmTypeList[0].field_parameter[field_index];
        break;
      case 6:
        field[i] = AtmTypeList[1].field_parameter[field_index];
        break;
      case 7:
        field[i] = AtmTypeList[2].field_parameter[field_index];
        break;
      case 8:
        field[i] = AtmTypeList[3].field_parameter[field_index];
        break;
      case 9:
        field[i] = AtmTypeList[4].field_parameter[field_index];
        break;
      case 15:
        field[i] = AtmTypeList[5].field_parameter[field_index];
        break;
      case 16:
        field[i] = AtmTypeList[6].field_parameter[field_index];
        break;
      case 17:
        field[i] = AtmTypeList[7].field_parameter[field_index];
        break;
      case 35:
        field[i] = AtmTypeList[8].field_parameter[field_index];
        break;
      case 53:
        field[i] = AtmTypeList[9].field_parameter[field_index];
        break;
      default:
        return 0;
    }
  }

  return 1;
}

/* 
 * The main procedure for the computation of comma2 descriptor
 */
static void
atomic_comma2_related_descriptor_computation_procedure (Molecule &m, double comma2_descriptor[], double ** coordinate)
{
  int n_atoms = m.natoms();

  double geometry_center[3];
  double * geometry_field = new double[n_atoms]; std::unique_ptr<double[]> free_geometry_field (geometry_field);
  double geometry_moment_of_inertia[3];

  double * hydrophobicity_index = new double[n_atoms]; std::unique_ptr<double[]> free_hydrophobicity_index (hydrophobicity_index);

  double * field = new double[n_atoms]; std::unique_ptr<double[]> free_field (field);
  double field_comma2_descriptor[7];

  double ** geometry_moment_of_inertia_axis = new double *[3];
  
  for (int i=0; i<3; i++)
    geometry_moment_of_inertia_axis[i] = new double[3];
  
  set_vector (geometry_field, n_atoms, 1.0);

  compute_field_center (n_atoms, coordinate, geometry_center, geometry_field);
  compute_moment_of_inertia (m, n_atoms, coordinate, geometry_center, geometry_field, geometry_moment_of_inertia_axis, geometry_moment_of_inertia);
  
  // compute the comma descriptors (COMMA2 using ghose and crippen hydrophobicity index as field
  ghose_crippen_procedure (m, hydrophobicity_index);
  
  comma2_related_descriptor_computation_procedure (m, n_atoms, hydrophobicity_index, field_comma2_descriptor, coordinate, geometry_center, geometry_moment_of_inertia_axis);

  for (int i=0; i<7; i++)
    comma2_descriptor[i] = field_comma2_descriptor[i];

  // for three other field (vdw volume, atomic electronegativity and atomic polarizability) compute COMMA2 descriptors
  for (int i=0; i<3; i++)
    {
      assign_field_for_at_comma_computation (m, field, i);
      comma2_related_descriptor_computation_procedure (m, n_atoms, field, field_comma2_descriptor, coordinate, geometry_center, geometry_moment_of_inertia_axis);

      for (int j=0; j<7; j++)
        comma2_descriptor[(i+1)*7 +j] = field_comma2_descriptor[j];
    }

  for (int i=0; i<3; i++)
    delete[] geometry_moment_of_inertia_axis[i];
  
  delete[] geometry_moment_of_inertia_axis;
}

/* 
 * This function compute the CoMMA related descriptors
 */

static void
comma_related_descriptor_computation_procedure (Molecule &m,
                                double atom_partial_charge[],
                                double comma_descriptor[])
{
  for (int i=0; i<NUMBER_OF_COMMA_DESCRIPTOR; i++)
    comma_descriptor[i] = 0.0;

  int n_atoms = m.natoms();
  double ** coordinate_array = new double *[n_atoms];
  double comma2_descriptor[28];

  for (int i=0; i<n_atoms; i++)
    {
      coordinate_array[i] = new double[3];
      coordinate_array[i][0] = m.x(i);
      coordinate_array[i][1] = m.y(i);
      coordinate_array[i][2] = m.z(i);
    }

  double ratio = 1;
  // compute the comma descriptors (COMMA using geometry center and gravity center)
  comma_related_descriptor_computation_procedure (m, atom_partial_charge, comma_descriptor, coordinate_array, ratio);

  // compute 28 atomic comma2 descriptors (geometry center and geometry_meoment_of_inertia_axis 
  // as the reference point and axis frame, respectively)
  atomic_comma2_related_descriptor_computation_procedure (m, comma2_descriptor, coordinate_array);

  for (int i=0; i<28; i++)
    comma_descriptor[24+i] = comma2_descriptor[i];

  comma_descriptor[52] = ratio;

  for (int i=0; i<n_atoms; i++)
    delete[] coordinate_array[i];

  delete[] coordinate_array;
}

static int
compute_hydrophobic_jurs_procedure (Molecule &m, area_t atom_area[], double hydrophobic_jurs[])
{
  int natoms = m.natoms();
  double * hydrophobicity_index = new double[natoms]; std::unique_ptr<double[]> free_hydrophobicity_index (hydrophobicity_index);
  ghose_crippen_procedure (m, hydrophobicity_index);

  //  cerr<<m.smiles()<<endl;
  //  for (int i=0; i<natoms; i++)
  // cerr<<"Atom "<<i<<"   atomic_symbol "<<m.atomi(i)->atomic_symbol()<<"  Surface area="<<atom_area[i]<<"  hydro="<<hydrophobicity_index[i]<<endl;

  set_vector (hydrophobic_jurs, NUMBER_OF_HYDROPHOBIC_JURS, 0.0);

  double total_hydrophobic_constant = 0.0;
  double total_hydrophilic_constant = 0.0;
  double total_surface_area = 0.0;
  double most_hydrophobic_constant = 0.0;
  double most_hydrophilic_constant = 0.0;

  double surface_area_most_hydrophobic = 0.0;
  double surface_area_most_hydrophilic = 0.0;

  for (int i=0; i<natoms; i++)
    {
      total_surface_area += atom_area[i];

      if (hydrophobicity_index[i] > 0.0)
        {
          // total of hydrophobic surface area  PPHS-1
          hydrophobic_jurs[0] += atom_area[i];
          // total hydrophobic constant
          total_hydrophobic_constant += hydrophobicity_index[i];
          // atomic constant weighted PPHS-3;
          hydrophobic_jurs[4] += atom_area[i] * hydrophobicity_index[i];

          if (hydrophobicity_index[i] > most_hydrophobic_constant)
            {
              most_hydrophobic_constant = hydrophobicity_index[i];
              surface_area_most_hydrophobic = atom_area[i];
            }
          else if (hydrophobicity_index[i] == most_hydrophobic_constant) 
            {
              if (atom_area[i] > surface_area_most_hydrophobic) 
                surface_area_most_hydrophobic = atom_area[i];
            }
        }
      else if (hydrophobicity_index[i] < 0.0)
        {
          // total of hydrophilic surface area PNHS-1
          hydrophobic_jurs[1] += atom_area[i];
          // total hydrophilic constant
          total_hydrophilic_constant += hydrophobicity_index[i];
          // atomic constant weighted PNHS-3;
          hydrophobic_jurs[5] += atom_area[i] * hydrophobicity_index[i];

          if (hydrophobicity_index[i] < most_hydrophilic_constant)
            {
              most_hydrophilic_constant = hydrophobicity_index[i];
              surface_area_most_hydrophilic = atom_area[i];
            }
          else if (hydrophobicity_index[i] == most_hydrophilic_constant) 
            {
              if (atom_area[i] > surface_area_most_hydrophilic) 
                surface_area_most_hydrophilic = atom_area[i];
            }
        }
    }
  // total hydrophobic weighted PPHS  PPHS-2
  hydrophobic_jurs[2] = total_hydrophobic_constant * hydrophobic_jurs[0];
  // total hydrophilic weighted PNHS  PNHS-2
  hydrophobic_jurs[3] = total_hydrophilic_constant * hydrophobic_jurs[1];
  // difference between PPHS-1 and PNHS-1
  hydrophobic_jurs[6] = hydrophobic_jurs[0] - hydrophobic_jurs[1];
  // difference between PPHS-2 and PNHS-2
  hydrophobic_jurs[7] = hydrophobic_jurs[2] - hydrophobic_jurs[3];
  // difference between PPHS-3 and PNHS-3
  hydrophobic_jurs[8] = hydrophobic_jurs[4] - hydrophobic_jurs[5];
  // fractional hydrophobic, hydrophilic surface area;
  for (int i=0; i<6; i++)
    hydrophobic_jurs[9+i] = hydrophobic_jurs[i] / total_surface_area;
  // surace weighted hydrophobic/hydrophilic surface area
  for (int i=0; i<6; i++)
    hydrophobic_jurs[15+i] = hydrophobic_jurs[i] * total_surface_area / 1000.0;

  // relative hydrophobicity
  if (0.0 == total_hydrophobic_constant)
    hydrophobic_jurs[21] = 1.0;
  else
    hydrophobic_jurs[21] = most_hydrophobic_constant / total_hydrophobic_constant;

  // relative hydrophilicity
  if (0 == total_hydrophilic_constant)
    hydrophobic_jurs[22] = 1.0;
  else
    hydrophobic_jurs[22] = most_hydrophilic_constant / total_hydrophilic_constant;

  // relative hydrophobic surface area
  hydrophobic_jurs[23] = surface_area_most_hydrophobic * hydrophobic_jurs[21];
  // relative hydrophilic surface area
  hydrophobic_jurs[24] = surface_area_most_hydrophilic * hydrophobic_jurs[22];

  // SURR1 SURR2 and SURR3
  if (0.0 != hydrophobic_jurs[0])
    hydrophobic_jurs[25] = hydrophobic_jurs[1] / hydrophobic_jurs[0];
  if (0.0 != hydrophobic_jurs[2])
    hydrophobic_jurs[26] = hydrophobic_jurs[3] / hydrophobic_jurs[2];
  if (0.0 != hydrophobic_jurs[4])
    hydrophobic_jurs[27] = hydrophobic_jurs[5] / hydrophobic_jurs[4];

  return 1;
}


/* 
 * This function is to calculate the charge surface area related descriptors for molecule
 */

static int

jwsa (Molecule &m, 
      IWString_and_File_Descriptor & output, 
      double atom_partial_charge[],
      double comma_partial_charge[],
      area_t atom_area[],
      int da_results[],
      double comma_descriptor[],
      double hydrophobic_jurs[])
{

  int rc;  // return status

  int n_atoms = m.natoms();
  
  // compute partial atomic charge -- only needed when jurs descriptor or electronic_descriptor are needed
  if ((compute_jurs_descriptor) || (compute_electronic_descriptor))
  {
    if (! use_gasteiger_partial_charge)
      rc = m.compute_Abraham_partial_charges();
    else
      rc = m.compute_Gasteiger_partial_charges();

    if (0==rc) 
    {
      if (verbose)
        cerr<<m.name()<<"\tERROR in calculation of partial charge\n";
      return rc;
    }

    for (int i=0; i<n_atoms; i++) {
      atom_partial_charge[i]= m.charge_on_atom (i);
    }
  }

  if (compute_comma_related_descriptor)
  {
    if (! use_gasteiger_partial_charge)
    {
      rc = m.compute_Gasteiger_partial_charges();
        
      if (0==rc) 
      {
        if (verbose)
          cerr<<m.name()<<"\tERROR in calculation of partial charge\n";
        return rc;
      }
    }
      
    for (int i=0; i<n_atoms; i++)
      comma_partial_charge[i]= m.charge_on_atom (i);

    comma_related_descriptor_computation_procedure (m, comma_partial_charge, comma_descriptor);
  }

  int name_written = 0;

  Sparse_Fingerprint_Creator sfc;

  // calculate the surface area, result stored in atom_area

  if (compute_jurs_descriptor || compute_savol_like_descriptor || compute_aromaticity_related_descriptor || compute_electronic_descriptor || compute_hydrophobic_jurs_descriptor)
  {
    volume_t molecule_volume = 0;
    rc=molecule_surface_area (m, atom_area, molecule_volume);
    
    if (0==rc) 
    {
      cerr<<m.name()<<"\tERROR in calculation of surface area\n";
      return rc;
    }

    (void) compute_hydrophobic_jurs_procedure (m, atom_area, hydrophobic_jurs);
    
    // used for jurs descriptor calculation & electronic descriptor calculation
    double total_surface_area = 0.0;
    double total_positive_charge = 0.0;  
    double most_positive_charge = 0.0;
    double total_negative_charge = 0.0;
    double most_negative_charge = 0.0;
    double surface_area_of_most_positive_atom = 0.0;
    double surface_area_of_most_negative_atom = 0.0;
    
    // used for dipolar moment calculation -- electronic descriptors
    double positive_weight_x=0;
    double positive_weight_y=0;
    double positive_weight_z=0;
    double negative_weight_x=0;
    double negative_weight_y=0;
    double negative_weight_z=0;
    
    // used for electronic descriptor calculation only
    int positive_atom_number = 0;
    int negative_atom_number = 0;
    double sum_square_positive_charge = 0.0;
    double sum_square_negative_charge = 0.0;
    
    // jurs descriptor, total of 30
    // six initial one (3 pairs) partial positive (negative) surface area
    // total charge weighted positive (negative) surface area
    // atomic charge weighted positive (negative) surface area -- 6 -- total 6
    double PPSA_1 = 0.0, PPSA_2 = 0.0, PPSA_3 =0.0;
    double PNSA_1 = 0.0, PNSA_2 = 0.0, PNSA_3 =0.0; 
    
    // total hydrophobic surface area (absolute partial charge<0.2) -- 1 -- total 26 
    // total polar surface area (absolute partial charge>=0.2)  -- 1 -- total 27
    double TASA = 0.0;
    double TPSA = 0.0;
    
    // do the calculation (mainly, summation)
    for (int i=0; i<n_atoms; i++)
    {
      total_surface_area += atom_area[i];
      if (atom_partial_charge[i] >0.0) 
      {
        total_positive_charge += atom_partial_charge[i];
        PPSA_1 += atom_area[i];
        PPSA_3 += atom_area[i] * atom_partial_charge[i];
        if (atom_partial_charge[i]>most_positive_charge) 
        {
          most_positive_charge = atom_partial_charge[i];
          surface_area_of_most_positive_atom = atom_area[i];
        }
        if (atom_partial_charge[i] <0.2) TASA += atom_area[i];
        else TPSA += atom_area[i];
        
        // for electronic descriptors
        positive_atom_number++;
        sum_square_positive_charge += atom_partial_charge[i]*atom_partial_charge[i];
        
        // for dipolar moment (related)
        positive_weight_x += atom_partial_charge[i]*m.x(i);
        positive_weight_y += atom_partial_charge[i]*m.y(i);
        positive_weight_z += atom_partial_charge[i]*m.z(i);
      }
      if (atom_partial_charge[i] <0.0) 
      {
        total_negative_charge += atom_partial_charge[i];
        PNSA_1 += atom_area[i];
        PNSA_3 += atom_area[i] * atom_partial_charge[i];
        if (atom_partial_charge[i]<most_negative_charge) 
        {
          most_negative_charge = atom_partial_charge[i];
          surface_area_of_most_negative_atom = atom_area[i];
        }
        if (atom_partial_charge[i] >-0.2) TASA += atom_area[i];
        else TPSA += atom_area[i];
        
        // for electronic descriptors
        negative_atom_number++;
        sum_square_negative_charge += atom_partial_charge[i]*atom_partial_charge[i];
        
        // for dipolar moment (related)
        negative_weight_x += atom_partial_charge[i]*m.x(i);
        negative_weight_y += atom_partial_charge[i]*m.y(i);
        negative_weight_z += atom_partial_charge[i]*m.z(i);
      }
    }
      
    // center of positive charge and negative charge
    double JW_DP_1 = 0.0;
    double JW_DP_2 = 0.0;
    
    if ((total_positive_charge != 0) && (total_negative_charge !=0))
    {
      double cp_x = positive_weight_x / total_positive_charge;
      double cp_y = positive_weight_y / total_positive_charge;
      double cp_z = positive_weight_z / total_positive_charge;
      
      double cn_x = negative_weight_x / total_negative_charge;
      double cn_y = negative_weight_y / total_negative_charge;
      double cn_z = negative_weight_z / total_negative_charge;
      // dipolar moment for the molecule & square of dipolar moment  -- electronic descriptor --3 -- total 19
      double pn_distance = sqrt((cp_x-cn_x)*(cp_x-cn_x)+(cp_y-cn_y)*(cp_y-cn_y)+(cp_z-cn_z)*(cp_z-cn_z));
      JW_DP_1 = total_positive_charge * pn_distance *4.8;
      JW_DP_2 = JW_DP_1*JW_DP_1;
    }
      
      // topological electronic index
    double JW_DP_3 = 0.0;
    for (int i=0; i<n_atoms; i++) 
      for (int j=i+1; j<n_atoms;j++)
      {
        if (m.distance_between_atoms(i,j) !=0.0)
        {
          if (atom_partial_charge[i]>atom_partial_charge[j])
            JW_DP_3 += (atom_partial_charge[i]-atom_partial_charge[j])/(m.distance_between_atoms(i,j)*m.distance_between_atoms(i,j));
          else
            JW_DP_3 += (atom_partial_charge[j]-atom_partial_charge[i])/(m.distance_between_atoms(i,j)*m.distance_between_atoms(i,j));
        }
      }
    
    PPSA_2 = PPSA_1 * total_positive_charge;
    PNSA_2 = PNSA_1 * total_negative_charge;
    // differential of the initial 3 pairs  -- 3 -- total 9
    double DPSA_1 = PPSA_1 - PNSA_1;
    double DPSA_2 = PPSA_2 - PNSA_2;
    double DPSA_3 = PPSA_3 - PNSA_3;
    
    // fractional of the initial 6 (3 pairs)  -- 6 -- total 15
    double FPSA_1 = 0.0;
    double FNSA_1 = 0.0;
    double FPSA_2 = 0.0;
    double FNSA_2 = 0.0;
    double FPSA_3 = 0.0;
    double FNSA_3 = 0.0;
    
    if (0 !=total_surface_area)
    {
      FPSA_1 = PPSA_1/total_surface_area;
      FNSA_1 = PNSA_1/total_surface_area;
      FPSA_2 = PPSA_2/total_surface_area;
      FNSA_2 = PNSA_2/total_surface_area;
      FPSA_3 = PPSA_3/total_surface_area;
      FNSA_3 = PNSA_3/total_surface_area;
    }
    
    // surface weighted charged partial surface area  -- 6 -- total 21
    double WPSA_1 = PPSA_1*total_surface_area/1000.0;
    double WNSA_1 = PNSA_1*total_surface_area/1000.0;
    double WPSA_2 = PPSA_2*total_surface_area/1000.0;
    double WNSA_2 = PNSA_2*total_surface_area/1000.0;
    double WPSA_3 = PPSA_3*total_surface_area/1000.0;
    double WNSA_3 = PNSA_3*total_surface_area/1000.0;
    
    // relative postive (negative) charge -- 2 -- total 23
    double RPCG = 0.0;
    double RNCG = 0.0;
    if (0.0 != total_positive_charge)
    {
      RPCG = most_positive_charge/total_positive_charge;
      RNCG = most_negative_charge/total_negative_charge;
    }
    // relative postive (negative) charge surface area -- 2 -- total 25
    double RPCS = surface_area_of_most_positive_atom*RPCG;
    double RNCS = surface_area_of_most_negative_atom*RNCG;
    
    // relative hydrophobic surface area (absolute partial charge<0.2) -- 1 -- total 28 
    // relative polar surface area (absolute partial charge>=0.2)  -- 1 -- total 29
    
    double RASA = 0.0;
    double RPSA = 0.0;
    if (0.0 != total_surface_area)
    {
      RASA = TASA / total_surface_area;
      RPSA = TPSA / total_surface_area;
    }
    // total molecular solvent-accessible surface area -- 1 -- total 30
    double SASA = total_surface_area;
    
    // Electronic descriptors -- Total of 19
    double JWED_1 = 0.0;
    double JWED_2 = 0.0;
    double JWED_3 = 0.0;
    double JWED_4 = 0.0;
    double JWED_5 = 0.0;
      
    // Most absolue, positive, negative charge & their corresponding average  -- 6 -- total 11
    double JWED_6 = 0.0;
    double JWED_7 = 0.0;
    double JWED_8 = 0.0;
    double JWED_9 = 0.0;
    double JWED_10 = 0.0;
    double JWED_11 = 0.0;
    
    // squares of absolute, positive and negative charge and their average -- 5 -- total 16
    double JWED_12 = 0.0;
    double JWED_13 = 0.0;
    double JWED_14 = 0.0;
    double JWED_15 = 0.0;
    double JWED_16 = 0.0;
    
    // Total absolute , positive, negative charge & their corresponding average  -- 5 -- total 5
    JWED_1 = total_positive_charge - total_negative_charge;
    if (0 != n_atoms)
      JWED_2 = JWED_1 / n_atoms;
    JWED_3 = total_positive_charge; // total_negative_charge is just the opposite
    if ( 0!= positive_atom_number)
      JWED_4 = JWED_3 / positive_atom_number ;
    if (0!= negative_atom_number)
      JWED_5 = total_negative_charge / negative_atom_number;
    
    // Most absolue, positive, negative charge & their corresponding average  -- 6 -- total 11
    JWED_6 = most_positive_charge;
    if (0 != positive_atom_number)
      JWED_7 = JWED_6 / positive_atom_number;
    JWED_8 = most_negative_charge;
    if (0!= negative_atom_number)
      JWED_9 = JWED_8 / negative_atom_number;
    if (most_positive_charge>-most_negative_charge) JWED_10 = most_positive_charge;
    else JWED_10 = -most_negative_charge;
    if (0 != n_atoms)
      JWED_11 = JWED_10 / n_atoms;
    
    // squares of absolute, positive and negative charge and their average -- 5 -- total 16
    JWED_12 = sum_square_positive_charge;
    if (0 != positive_atom_number)
      JWED_13 = JWED_12 / positive_atom_number;
    JWED_14 = sum_square_negative_charge;
    if (0 != negative_atom_number)
      JWED_15 = JWED_14 / negative_atom_number;
    if (0 != n_atoms)
      JWED_16 = (JWED_12 + JWED_14) / n_atoms;
    
    // compute the savol like descriptors
    double JW_SovEn = 0.0;
    double JW_Vol = molecule_volume;
    double JW_PSA1 = 0.0;
    double JW_PSA2 = 0.0;
    double JW_PSA3 = 0.0;
    double JW_HpbSA = 0.0;
    
    if (compute_savol_like_descriptor)
    {
      rc = compute_solvation_energy (m, atom_area, JW_SovEn);
      if (0==rc) { if (verbose) cerr<<m.name()<<"\tERROR in calculation of solvation energy\n"; return rc;}
      rc = compute_polar_surface_area_1 (m, atom_area, JW_PSA1);
      if (0==rc) { if (verbose) cerr<<m.name()<<"\tERROR in calculation of polar surface area 1\n";return rc;}
      rc = compute_polar_surface_area_2 (m, atom_area, JW_PSA2);
      if (0==rc) { if (verbose) cerr<<m.name()<<"\tERROR in calculation of polar surface area 2\n";return rc;}
      rc = compute_polar_surface_area_3 (m, atom_area, JW_PSA3);
      if (0==rc) { if (verbose) cerr<<m.name()<<"\tERROR in calculation of polar surface area 3\n";return rc;}
      
      rc = compute_hydrophobic_surface_area (m, atom_area, JW_HpbSA);
      if (0==rc) {if (verbose) cerr<<m.name()<<"\tERROR in calculation of hydrophobic surface area\n";return rc;}
    }
    
    // start computing the aromaticity related descriptors
    double JW_AromSA=0.0;
    double JW_AromRa=0.0;
    if (compute_aromaticity_related_descriptor)
    {
      m.compute_aromaticity();
      for (int i=0;i<n_atoms; i++)
      {
        if (1==m.atomi(i)->atomic_number())
        { 
          if (m.is_aromatic(m.other(i,0))) JW_AromSA += atom_area[i];
        }
        else
          if (m.is_aromatic(i)) JW_AromSA += atom_area[i];
      }
      if (0.0!= total_surface_area)
        JW_AromRa = JW_AromSA / total_surface_area;
    }
    
    // start computing the H-bond donor acceptor -- jurs extention -- descriptors  Total of 11
    int number_of_hydrogen_donated = 0;
    double total_difference_in_charge_for_h_donor_pair = 0;
    
    double SSAH=0.0; // sum of the surface areas of hydrogens which can be donated
    double CHGD=0.0; // max charge difference between hydrogen and its corresponding heteroatom
    double ACGD = 0.0; // average difference in charge between all pairs of H-bond donors
    double SSAA = 0.0; // sum of the surface areas of all H-bond acceptor groups
    int CNTH = 0; // simple count of all H-bond donor groups
    int CNTA = 0; // simple count of all H-bond acceptor groups
    double RHTA = 0.0; // ratio of number of donor groups to number of acceptor groups
    double RSAH = 0.0; // average surface area of hydrogens which can be donated
    double RSAA = 0.0; // average surface area of H-bond acceptor groups
    double RSHM = 0.0; // fraction of the total molecular surface area associated with hydrogens which can be donated
    double RSAM = 0.0; // fraction of the total molecular surface area associated with H-bond acceptor groups
    
    //    cerr << "Before charge assigner, Number of Atom"<<n_atoms<<endl;
    if (charge_assigner.active())
      (void) charge_assigner.process (m);
    //    cerr << "After charge assigner,  Number of Atom"<<m.natoms()<<endl;
    if (donor_acceptor_assigner.active())
      donor_acceptor_assigner.process (m, da_results);
    //    cerr << "After donor_acceptor_assigner,  Number of Atom"<<m.natoms()<<endl;

    if (compute_jurs_descriptor)  // if jurs descriptor needs to be calculated, assign_donor_acceptor
    {
      for (int i=0; i<n_atoms; i++)
      {
        int neighbors = m.atomi(i)->ncon();
        switch (da_results[i])
        {
          case 1:
            // H-bond acceptors
            CNTA++;
            SSAA += atom_area[i];
            break;
            
          case 2:
            // Both acceptors & donors
            // as acceptor
            CNTA++;
            SSAA += atom_area[i];
            
            // as donors
            CNTH++;
            for (int j=0; j<neighbors; j++)
            {
              // if hydrogen, then doing addition
              int jth_atom = m.other(i,j);
              const Atom * atomj = m.atomi(jth_atom);
              if (atomj->atomic_number()==1)
              {
                number_of_hydrogen_donated ++;
                SSAH += atom_area[jth_atom];
                double difference_in_charge_for_h_donor_pair = atom_partial_charge[jth_atom]-atom_partial_charge[i];
                total_difference_in_charge_for_h_donor_pair += difference_in_charge_for_h_donor_pair;
                if (difference_in_charge_for_h_donor_pair > CHGD)
                  CHGD = difference_in_charge_for_h_donor_pair;
              }
            }
            break;
            
          case 3:
            // H-bond donor
            CNTH++;
            for (int j=0; j<neighbors; j++)
            {
              // if hydrogen, then doing addition
              int jth_atom = m.other(i,j);
              const Atom * atomj = m.atomi(jth_atom);
              if (atomj->atomic_number()==1)
              {
                number_of_hydrogen_donated ++;
                SSAH += atom_area[jth_atom];
                double difference_in_charge_for_h_donor_pair = atom_partial_charge[jth_atom]-atom_partial_charge[i];
                total_difference_in_charge_for_h_donor_pair += difference_in_charge_for_h_donor_pair;
                if (difference_in_charge_for_h_donor_pair > CHGD)
                  CHGD = difference_in_charge_for_h_donor_pair;
              }
            }
            break;
          default:
            break;
        }
      }
      
      // calculation of H-bond extension of jurs descriptor
      if (number_of_hydrogen_donated == 0)
      { 
        ACGD = 0.0;
        RSAH = 0.0;
      }  
      else
      { 
        ACGD = total_difference_in_charge_for_h_donor_pair / number_of_hydrogen_donated;
        RSAH = SSAH / number_of_hydrogen_donated;
      }
      
      if (CNTA!=0)
      { 
        RHTA = (double) CNTH / (double) CNTA;
        RSAA = SSAA / CNTA;
      }
      else 
      {
        if (0==CNTH) RHTA = 1;
        else RHTA = 2 * CNTA;  // WHEN CNTA=0, CNTH!=0, RHTA is set as 2* CNTA -- do not want to set as infinite
        RSAA = 0;
      }

      if (0 != total_surface_area)
      {
        RSHM = SSAH / total_surface_area;
        RSAM = SSAA / total_surface_area;
      }
    }

//  computations done, perform output

    if (compute_jurs_descriptor)
    {
      if (tag.length())
      {
        hit_bit(sfc, 0, SASA, min_jurs[0], max_jurs[0], bit_count_scaling_factor);
        hit_bit(sfc, 1, PPSA_1, min_jurs[1], max_jurs[1], bit_count_scaling_factor);
        hit_bit(sfc, 2, PNSA_1, min_jurs[2], max_jurs[2], bit_count_scaling_factor);
        hit_bit(sfc, 3, DPSA_1, min_jurs[3], max_jurs[3], bit_count_scaling_factor);
        hit_bit(sfc, 4, PPSA_2, min_jurs[4], max_jurs[4], bit_count_scaling_factor);
        hit_bit(sfc, 5, PNSA_2, min_jurs[5], max_jurs[5], bit_count_scaling_factor);
        hit_bit(sfc, 6, DPSA_2, min_jurs[6], max_jurs[6], bit_count_scaling_factor);
        hit_bit(sfc, 7, PPSA_3, min_jurs[7], max_jurs[7], bit_count_scaling_factor);
        hit_bit(sfc, 8, PNSA_3, min_jurs[8], max_jurs[8], bit_count_scaling_factor);
        hit_bit(sfc, 9, DPSA_3, min_jurs[9], max_jurs[9], bit_count_scaling_factor);
        hit_bit(sfc, 10, FPSA_1, min_jurs[10], max_jurs[10], bit_count_scaling_factor);
        hit_bit(sfc, 11, FNSA_1, min_jurs[11], max_jurs[11], bit_count_scaling_factor);
        hit_bit(sfc, 12, FPSA_2, min_jurs[12], max_jurs[12], bit_count_scaling_factor);
        hit_bit(sfc, 13, FNSA_2, min_jurs[13], max_jurs[13], bit_count_scaling_factor);
        hit_bit(sfc, 14, FPSA_3, min_jurs[14], max_jurs[14], bit_count_scaling_factor);
        hit_bit(sfc, 15, FNSA_3, min_jurs[15], max_jurs[15], bit_count_scaling_factor);
        hit_bit(sfc, 16, WPSA_1, min_jurs[16], max_jurs[16], bit_count_scaling_factor);
        hit_bit(sfc, 17, WNSA_1, min_jurs[17], max_jurs[17], bit_count_scaling_factor);
        hit_bit(sfc, 18, WPSA_2, min_jurs[18], max_jurs[18], bit_count_scaling_factor);
        hit_bit(sfc, 19, WNSA_2, min_jurs[19], max_jurs[19], bit_count_scaling_factor);
        hit_bit(sfc, 20, WPSA_3, min_jurs[20], max_jurs[20], bit_count_scaling_factor);
        hit_bit(sfc, 21, WNSA_3, min_jurs[21], max_jurs[21], bit_count_scaling_factor);

        hit_bit(sfc, 22, RPCG, min_jurs[22], max_jurs[22], bit_count_scaling_factor);
        hit_bit(sfc, 23, RNCG, min_jurs[23], max_jurs[23], bit_count_scaling_factor);
        hit_bit(sfc, 24, RPCS, min_jurs[24], max_jurs[24], bit_count_scaling_factor);
        hit_bit(sfc, 25, RNCS, min_jurs[25], max_jurs[25], bit_count_scaling_factor);
  
        hit_bit(sfc, 26, TPSA, min_jurs[26], max_jurs[26], bit_count_scaling_factor);
        hit_bit(sfc, 27, TASA, min_jurs[27], max_jurs[27], bit_count_scaling_factor);
        hit_bit(sfc, 28, RPSA, min_jurs[28], max_jurs[28], bit_count_scaling_factor);
        hit_bit(sfc, 29, RASA, min_jurs[29], max_jurs[29], bit_count_scaling_factor);
        hit_bit(sfc, 30, SSAH, min_jurs[30], max_jurs[30], bit_count_scaling_factor);
        hit_bit(sfc, 31, CHGD, min_jurs[31], max_jurs[31], bit_count_scaling_factor);
        hit_bit(sfc, 32, ACGD, min_jurs[32], max_jurs[32], bit_count_scaling_factor);
        hit_bit(sfc, 33, SSAA, min_jurs[33], max_jurs[33], bit_count_scaling_factor);
        hit_bit(sfc, 34, CNTH, min_jurs[34], max_jurs[34], bit_count_scaling_factor);
        hit_bit(sfc, 35, CNTA, min_jurs[35], max_jurs[35], bit_count_scaling_factor);
  
        hit_bit(sfc, 36, RHTA, min_jurs[36], max_jurs[36], bit_count_scaling_factor);
        hit_bit(sfc, 37, RSAH, min_jurs[37], max_jurs[37], bit_count_scaling_factor);
        hit_bit(sfc, 38, RSAA, min_jurs[38], max_jurs[38], bit_count_scaling_factor);
        hit_bit(sfc, 39, RSHM, min_jurs[39], max_jurs[39], bit_count_scaling_factor);
        hit_bit(sfc, 40, RSAM, min_jurs[40], max_jurs[40], bit_count_scaling_factor);
      }
      else
      {
        // output original jurs descritor
        if (! name_written)
        {
          append_first_token_of_name(m.name(), output);
          name_written = 1;
        }
        output << ' ' << SASA;
        output << ' ' << PPSA_1 <<' '<< PNSA_1<<' '<<DPSA_1;
        output << ' ' << PPSA_2 <<' '<< PNSA_2<<' '<<DPSA_2;
        output << ' ' << PPSA_3 <<' '<< PNSA_3<<' '<<DPSA_3;
        
        output << ' ' << FPSA_1 <<' '<< FNSA_1;
        output << ' ' << FPSA_2 <<' '<< FNSA_2;
        output << ' ' << FPSA_3 <<' '<< FNSA_3;
        
        output << ' ' << WPSA_1 <<' '<< WNSA_1;
        output << ' ' << WPSA_2 <<' '<< WNSA_2;
        output << ' ' << WPSA_3 <<' '<< WNSA_3;
        
        output << ' ' << RPCG <<' '<< RNCG;
        
        output << ' ' << RPCS <<' '<< RNCS;
        
        output << ' ' << TPSA <<' '<< TASA;
        output << ' ' << RPSA <<' '<< RASA;
        
        // output hydrogen_bond extention of jurs descriptors
        output <<' ' << SSAH<<' '<<CHGD<<' '<<ACGD<<' '<<SSAA<<' '<<CNTH<<' '<<CNTA;
        output <<' ' << RHTA<<' '<<RSAH<<' '<<RSAA<<' '<<RSHM<<' '<<RSAM;
      }
    }
    
    if (compute_electronic_descriptor)
    {
      if (tag.length())
      {
        hit_bit(sfc, 41, JWED_1, min_jurs[41], max_jurs[41], bit_count_scaling_factor);
        hit_bit(sfc, 42, JWED_2, min_jurs[42], max_jurs[42], bit_count_scaling_factor);
        hit_bit(sfc, 43, JWED_3, min_jurs[43], max_jurs[43], bit_count_scaling_factor);
        hit_bit(sfc, 44, JWED_4, min_jurs[44], max_jurs[44], bit_count_scaling_factor);
        hit_bit(sfc, 45, JWED_5, min_jurs[45], max_jurs[45], bit_count_scaling_factor);
        hit_bit(sfc, 46, JWED_6, min_jurs[46], max_jurs[46], bit_count_scaling_factor);
        hit_bit(sfc, 47, JWED_7, min_jurs[47], max_jurs[47], bit_count_scaling_factor);
        hit_bit(sfc, 48, JWED_8, min_jurs[48], max_jurs[48], bit_count_scaling_factor);
        hit_bit(sfc, 49, JWED_9, min_jurs[49], max_jurs[49], bit_count_scaling_factor);
        hit_bit(sfc, 50, JWED_10, min_jurs[50], max_jurs[50], bit_count_scaling_factor);
        hit_bit(sfc, 51, JWED_11, min_jurs[51], max_jurs[51], bit_count_scaling_factor);
        hit_bit(sfc, 52, JWED_12, min_jurs[52], max_jurs[52], bit_count_scaling_factor);
        hit_bit(sfc, 53, JWED_13, min_jurs[53], max_jurs[53], bit_count_scaling_factor);
        hit_bit(sfc, 54, JWED_14, min_jurs[54], max_jurs[54], bit_count_scaling_factor);
        hit_bit(sfc, 55, JWED_15, min_jurs[55], max_jurs[55], bit_count_scaling_factor);
        hit_bit(sfc, 56, JWED_16, min_jurs[56], max_jurs[56], bit_count_scaling_factor);
        hit_bit(sfc, 57, JW_DP_1, min_jurs[57], max_jurs[57], bit_count_scaling_factor);
        hit_bit(sfc, 58, JW_DP_2, min_jurs[58], max_jurs[58], bit_count_scaling_factor);
        hit_bit(sfc, 59, JW_DP_3, min_jurs[59], max_jurs[59], bit_count_scaling_factor);
      }
      else
      {
        if (! name_written)
        {
          append_first_token_of_name(m.name(), output);
          name_written = 1;
        }
        output << ' ' << JWED_1 <<' '<<JWED_2 <<' '<<JWED_3 <<' '<<JWED_4 <<' '<<JWED_5 ;
        output << ' ' << JWED_6 <<' '<<JWED_7 <<' '<<JWED_8 <<' '<<JWED_9 <<' '<<JWED_10<<' '<<JWED_11;
        output << ' ' << JWED_12<<' '<<JWED_13<<' '<<JWED_14<<' '<<JWED_15<<' '<<JWED_16;
        output << ' ' << JW_DP_1<<' '<<JW_DP_2<<' '<<JW_DP_3;
      }
    }
    
    if (compute_savol_like_descriptor)
    {
      if (tag.length())
      {
        hit_bit(sfc, 60, JW_SovEn, min_jurs[60], max_jurs[60], bit_count_scaling_factor);
        hit_bit(sfc, 61, JW_Vol, min_jurs[61], max_jurs[61], bit_count_scaling_factor);
        hit_bit(sfc, 62, JW_PSA1, min_jurs[62], max_jurs[62], bit_count_scaling_factor);
        hit_bit(sfc, 63, JW_PSA2, min_jurs[63], max_jurs[63], bit_count_scaling_factor);
        hit_bit(sfc, 64, JW_PSA3, min_jurs[64], max_jurs[64], bit_count_scaling_factor);
        hit_bit(sfc, 65, JW_HpbSA, min_jurs[65], max_jurs[65], bit_count_scaling_factor);
      }
      else
      {
        if (! name_written)
        {
          append_first_token_of_name(m.name(), output);
          name_written = 1;
        }
        output<<' ' << JW_SovEn<<' '<<  JW_Vol<<' '<< JW_PSA1<<' '<<JW_PSA2 ;
        output<<' ' << JW_PSA3 <<' '<<JW_HpbSA;
      }
    }
    
    if (compute_aromaticity_related_descriptor)
    {
      if (tag.length())
      {
        hit_bit(sfc, 66, JW_AromSA, min_jurs[66], max_jurs[66], bit_count_scaling_factor);
        hit_bit(sfc, 67, JW_AromRa, min_jurs[67], max_jurs[67], bit_count_scaling_factor);
      }
      else
      {
        if (! name_written)
        {
          append_first_token_of_name(m.name(), output);
          name_written = 1;
        }
        output<<' ' << JW_AromSA<<' '<<JW_AromRa;
      }

      if (compute_hydrophobic_jurs_descriptor)
      {
        if (tag.length())
        {
          for (int i = 0; i < NUMBER_OF_HYDROPHOBIC_JURS; ++i)
          {
            hit_bit(sfc, 68 + i, hydrophobic_jurs[i], min_jurs[68 + i], max_jurs[68 + i], bit_count_scaling_factor + i);
          }
        }
        else
        {
          if (! name_written)
          {
            append_first_token_of_name(m.name(), output);
            name_written = 1;
          }
          for (int i=0; i<NUMBER_OF_HYDROPHOBIC_JURS; i++) {
            output<<' ' << hydrophobic_jurs[i];
          }
        }
      }
    }
  }

  if (compute_comma_related_descriptor)
  {
    if (tag.length())
    {
      for (int i = 0; i < NUMBER_OF_COMMA_DESCRIPTOR; ++i)
      {
        double c = comma_descriptor[i];
        if (c > max_value_comma_descriptor)
          c = max_value_comma_descriptor;
        else if (c < min_value_comma_descriptor)
          c = min_value_comma_descriptor;

        hit_bit(sfc, 96 + i, c, min_value_comma_descriptor, max_value_comma_descriptor, bit_count_scaling_factor + i);
      }
    }
    else
    {
      if (! name_written)
        append_first_token_of_name(m.name(), output);

      for (int i=0; i<NUMBER_OF_COMMA_DESCRIPTOR; i++)
      {
        output << ' ';

        if (comma_descriptor[i] > max_value_comma_descriptor)
          output << max_value_comma_descriptor;
        else if (comma_descriptor[i] < min_value_comma_descriptor)
          output << min_value_comma_descriptor;
        else
          output << static_cast<float>(comma_descriptor[i]);
      }
    }
  }

  if (tag.length())
  {
    IWString tmp;
    sfc.daylight_ascii_form_with_counts_encoded(tag, output);

    output << tmp << '\n';
  }
  else
    output << '\n';

  return 1;
}

static int
jwsa (Molecule &m, 
      IWString_and_File_Descriptor & output)
{
  // calculate_partial_charge need to have all hydrogen added
  // but this has to be controlled by the input to the corina.  Other wise, molecule might have hydrogen but no coordinate
  //  m.make_implicit_hydrogens_explicit();
  
  int n_atoms = m.natoms();

  double * partial_charge = new double[n_atoms]; std::unique_ptr<double[]> free_partial_charge (partial_charge); 
  double * comma_partial_charge = new double[n_atoms]; std::unique_ptr<double[]> free_comma_partial_charge (comma_partial_charge);
  double comma_descriptor[NUMBER_OF_COMMA_DESCRIPTOR];
  double hydrophobic_jurs[NUMBER_OF_HYDROPHOBIC_JURS];

  area_t * area = new area_t[n_atoms]; std::unique_ptr<double[]> free_area (area);
  int * results = new_int(n_atoms); std::unique_ptr<int[]> free_results (results);

  int rc=jwsa (m, output, partial_charge, comma_partial_charge, area, results, comma_descriptor, hydrophobic_jurs);

  if (0==rc) 
    number_of_error++;

  return rc;
}


static int
jwsa (data_source_and_type<Molecule> & input, 
      IWString_and_File_Descriptor & output)
{
  Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose)
      cerr << molecules_read << " processing '" << m->name() << "'\n";

    (void) preprocess_molecule (*m);

    if (! jwsa (*m, output))
      continue;
    
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
jwsa (const char * fname, FileType input_type, 
      IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input (input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jwsa (input, output);
}

static int
jwsa_filter_record (const const_IWSubstring & smiles,
                    IWString_and_File_Descriptor & output)
{
  Molecule m;

  if (! m.build_from_smiles(smiles))
  {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  return jwsa (m, output);
}

static int
jwsa_filter (iwstring_data_source & input,
             IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(smiles_tag))
      continue;

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    if (! jwsa_filter_record (buffer, output))
      return 0;
  }

  output.flush();

  return 1;
}


static int
jwsa_filter (const char * fname,
             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jwsa_filter (input, output);
}

static void
process_dash_F_option (Command_Line &cl)
{
  const_IWSubstring f;
  int i = 0;
  while (cl.value ('F', f, i++))
  {
    if (! read_control_file (f))
    {
      cerr << "Cannot read control file '" << f << "'\n";
      usage (19);
    }
  }
  
  if (verbose)
      cerr << "Read " << queries.number_elements() << " queries\n";
  
  int nq = queries.number_elements();
  assert (nq > 0);
  
  Accumulator<double> values;
  
  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_one_embedding_per_atom (1);
    double v;
    if (! queries[i]->numeric_value (v))
      {
        cerr << "Yipes, no numeric value for query " << i << " '" << queries[i]->comment() << "'\n";
        usage (19);
      }
    
    values.extra (v);
  }
}

static int
jwsa (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:ei:sH:ar:jV:bGd:hN:F:ng:M:J:y:");
  
  if (cl.unrecognised_options_encountered())
    {
      cerr << "Unrecognised options encountered\n";
      usage (1);
    }
  
  verbose = cl.option_count ('v');

  set_default_iwstring_float_concatenation_precision(3);
  set_default_iwstring_double_concatenation_precision(4);

  use_gasteiger_partial_charge = cl.option_count ('G');

  if (! process_elements (cl))
  {
    usage (2);
  }
  
  if (! process_standard_aromaticity_options (cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage (5);
  }
  
  if (cl.option_present ('r')) 
  {
    if (! cl.value('r', probe_radius))
    {
      cerr << "Invalid value for probe radius (-r) has to be a number\n";
      usage (2);
    }
    
    if ((probe_radius<0) || (probe_radius>10))
    {
      cerr << "The value for probe radius is not valid, should be between 0 and 10 \n";
      cerr << "1.4 is used instead (water as a probe)\n";
      probe_radius = 1.4;
    }
  }
  
  if(cl.option_present ('V'))
  {
    if (! set_default_van_der_waals_radius_type (cl, 'V', vdw_radius_type, verbose))
    {
      cerr << "Cannot determine Van der Waals specification (-V option)\n";
      usage (11);
    }
  }
  
  if (cl.option_present ('j')) compute_jurs_descriptor = 1;
  if (cl.option_present ('e')) compute_electronic_descriptor = 1;
  if (cl.option_present ('s')) compute_savol_like_descriptor = 1;
  if (cl.option_present ('a')) compute_aromaticity_related_descriptor = 1;
  if ((cl.option_present ('b')) || (cl.option_present ('h')))
  {  
    if (! cl.option_present ('F'))
    {
      cerr << "Must specify one or more control files via the -F option when comma or hydrophobic jurs descriptor need to be computed\n";
      usage (19);
    }
      
    process_dash_F_option (cl);

    if (cl.option_present ('b'))
      compute_comma_related_descriptor = 1;
    if (cl.option_present ('h'))
      compute_hydrophobic_jurs_descriptor = 1;
  }

  if (!(compute_jurs_descriptor || compute_electronic_descriptor || compute_savol_like_descriptor || compute_aromaticity_related_descriptor || compute_comma_related_descriptor || compute_hydrophobic_jurs_descriptor))
  {
    cerr << "Default all of -j, -e, -s and -a will be turned on\n";
    compute_jurs_descriptor = 1;
    compute_electronic_descriptor = 1;
    compute_savol_like_descriptor = 1;
    compute_aromaticity_related_descriptor = 1;

    cerr << "Since -F option is on, -b and -h will be turned on too\n";

    if (cl.option_present ('F'))
    {
      process_dash_F_option (cl);

      compute_comma_related_descriptor = 1;
      compute_hydrophobic_jurs_descriptor = 1;
    }
  }

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (m.starts_with("min="))
      {
        m.remove_leading_chars(4);
        if (! m.numeric_value(min_value_comma_descriptor))
        {
          cerr << "Invalid minimum comma descriptor value '" << m << "'\n";
          usage(4);
        }
      }
      else if (m.starts_with("max="))
      {
        m.remove_leading_chars(4);
        if (! m.numeric_value(max_value_comma_descriptor))
        {
          cerr << "Invalid maximum comma descriptor value '" << m << "'\n";
          usage(4);
        }
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        usage(3);
      }
    }

    if (verbose)
      cerr << "Comma descriptors restricted to the range " << min_value_comma_descriptor << " to " << max_value_comma_descriptor << endl;

    if (min_value_comma_descriptor > max_value_comma_descriptor)
    {
      cerr << "Invalid range for comma descriptors " << min_value_comma_descriptor << " to " << max_value_comma_descriptor << endl;
      usage(4);
    }
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
  else if (cl.option_present('J'))
  {
    cl.value('J', tag);

    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "Will produce a fingerprint with tag '" << tag << "'\n";

    if (cl.option_present('y'))
    {
      if (! cl.value('y', bit_count_scaling_factor) || bit_count_scaling_factor <= 0.0)
      {
        cerr << "The bit count scaling factor option (-y) must be a positive value\n";
        usage(2);
      }

      if (verbose)
        cerr << "Bit count scaling factor set to " << bit_count_scaling_factor << endl;
    }
  }
  else if (! all_files_recognised_by_suffix (cl))
    return 4;

  if (cl.option_present ('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line (cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor/acceptor assignment object\n";
      usage (4);
    }
  }
  
  if (cl.option_present ('N'))
  {
    if (! charge_assigner.construct_from_command_line (cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage (33);
    }
  }
  
  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);
  output.resize(36000);

  if (0 == tag.length())
    output_result_header (output);
  
  int rc = 0;
  
  if (tag.length())
  {
    if (! jwsa_filter (cl[0], output))
      rc = 1;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! jwsa (cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

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
  int rc = jwsa (argc, argv);
  return rc;
}
