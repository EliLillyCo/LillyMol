/*
  We often need to assing Van-der-Waals radii
*/

#include <stdlib.h>

#include "cmdline.h"

#include "molecule.h"
#include "iw_vdw.h"

/*
  VDW parameters from Shrake and Rupley, J Mol Bio, 79, 351-371 (1973)
*/

static vdw_radius_t sr_vdw_nitrogen = 1.50;
static vdw_radius_t sr_vdw_oxygen   = 1.40;
static vdw_radius_t sr_vdw_sulphur  = 1.85;
static vdw_radius_t sr_vdw_ch       = 2.00;
static vdw_radius_t sr_vdw_aromc    = 1.85;
static vdw_radius_t sr_vdw_c        = 1.50;

// Parameters from savol3

static vdw_radius_t savol_vdw_nitrogen = 1.55;
static vdw_radius_t savol_vdw_oxygen   = 1.52;
static vdw_radius_t savol_vdw_sulphur  = 1.80;
static vdw_radius_t savol_vdw_c        = 1.70;

// These next parameters are common and come from savol3

static vdw_radius_t vdw_p        = 1.80;
static vdw_radius_t vdw_f        = 1.40;
static vdw_radius_t vdw_cl       = 1.75;
static vdw_radius_t vdw_br       = 1.95;
static vdw_radius_t vdw_i        = 2.10;
static vdw_radius_t vdw_h        = 1.20;
static vdw_radius_t vdw_si       = 2.10;
static vdw_radius_t vdw_ca       = 0.99;

// Parameters from Wikipedia, Dec 2005, Jeff Sutherland

static vdw_radius_t wiki_hydrogen   = 1.20;
static vdw_radius_t wiki_carbon     = 1.70;
static vdw_radius_t wiki_nitrogen   = 1.55;
static vdw_radius_t wiki_oxygen     = 1.40;
static vdw_radius_t wiki_fluorine   = 1.35;
static vdw_radius_t wiki_phosphorus = 1.90;
static vdw_radius_t wiki_sulphur    = 1.85;
static vdw_radius_t wiki_chlorine   = 1.80;

// J Phys Chem A. 2009 May 14; 113(19): 5806-5812. 

static vdw_radius_t mantina_boron = 1.92;
static vdw_radius_t mantina_selenium = 1.95;    // kind of interpolate their table 7

static int
assign_wiki_vdw_radii(Molecule & m,
                      vdw_radius_t * vdw)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    if (1 == z)
      vdw[i] = wiki_hydrogen;
    else if (6 == z)
      vdw[i] = wiki_carbon;
    else if (7 == z)
      vdw[i] = wiki_nitrogen;
    else if (8 == z)
      vdw[i] = wiki_oxygen;
    else if (16 == z)
      vdw[i] = wiki_sulphur;
    else if (9 == z)
      vdw[i] = wiki_fluorine;
    else if (17 == z)
      vdw[i] = wiki_chlorine;
    else if (35 == z)
      vdw[i] = vdw_br;
    else if (53 == z)
      vdw[i] = vdw_i;
    else if (15 == z)
      vdw[i] = wiki_phosphorus;
    else if (20 == z)
      vdw[i] = vdw_ca;
    else if (14 == z)
      vdw[i] = vdw_si;
    else if (5 == z)
      vdw[i] = mantina_boron;
    else if (34 == z)
      vdw[i] = mantina_selenium;
    else
    {
      cerr << "assign_wiki_vdw_radii:: what kind of atom is this? i = " << i <<
              " z = " << z << " connections = " << m.ncon(i) << endl;
      return 0;
    }
  }

  return 1;
}

static int
assign_savol_vdw_radii (Molecule & m,
                        vdw_radius_t * vdw)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    if (1 == z)
      vdw[i] = vdw_h;
    else if (6 == z)
      vdw[i] = savol_vdw_c;
    else if (7 == z)
      vdw[i] = savol_vdw_nitrogen;
    else if (8 == z)
      vdw[i] = savol_vdw_oxygen;
    else if (16 == z)
      vdw[i] = savol_vdw_sulphur;
    else if (9 == z)
      vdw[i] = vdw_f;
    else if (17 == z)
      vdw[i] = vdw_cl;
    else if (35 == z)
      vdw[i] = vdw_br;
    else if (53 == z)
      vdw[i] = vdw_i;
    else if (15 == z)
      vdw[i] = vdw_p;
    else if (20 == z)
      vdw[i] = vdw_ca;
    else if (14 == z)
      vdw[i] = vdw_si;
    else if (5 == z)
      vdw[i] = mantina_boron;
    else if (34 == z)
      vdw[i] = mantina_selenium;
    else
    {
      cerr << "assign_savol_vdw_radii:: what kind of atom is this? i = " << i <<
              " z = " << z << " connections = " << m.ncon(i) << endl;
      return 0;
    }
  }

  return 1;
}


static int
assign_shrake_and_rupley_vdw_radii (Molecule & m,
                                    vdw_radius_t * vdw)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);
    if (1 == z)
      vdw[i] = vdw_h;
    else if (6 == z)
    {
      aromaticity_type_t arom;
      if (m.aromaticity(i, arom) && IS_AROMATIC_ATOM(arom))
        vdw[i] = sr_vdw_aromc;
      else if (m.hcount(i))
        vdw[i] = sr_vdw_ch;
      else
        vdw[i] = sr_vdw_c;
    }
    else if (7 == z)
      vdw[i] = sr_vdw_nitrogen;
    else if (8 == z)
      vdw[i] = sr_vdw_oxygen;
    else if (16 == z)
      vdw[i] = sr_vdw_sulphur;
    else if (9 == z)
      vdw[i] = vdw_f;
    else if (17 == z)
      vdw[i] = vdw_cl;
    else if (35 == z)
      vdw[i] = vdw_br;
    else if (53 == z)
      vdw[i] = vdw_i;
    else if (15 == z)
      vdw[i] = vdw_p;
    else if (5 == z)
      vdw[i] = mantina_boron;
    else if (34 == z)
      vdw[i] = mantina_selenium;
    else
    {
      cerr << "assign_vdw_radii:: what kind of atom is this? i = " << i <<
              " z = " << z << " connections = " << m.ncon(i) << endl;
      return 0;
    }
  }

  return 1;
}

static int
assign_molvol_vdw_radii (Molecule & m,
                         vdw_radius_t * vdw)
{
  int matoms = m.natoms();

  int rc = 1;      // will be set to 0 if anything is inclassified

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    if (1 == z)
    {
      atom_number_t on = m.other(i, 0);
      atomic_number_t z = m.atomic_number(on);

      if (7 == z)
        vdw[i] = 1.125;
      else if (8 == z)
        vdw[i] = 1.10;
      else
        vdw[i] = 1.5;
    }
    else if (6 == z)
      vdw[i] = 1.9;
    else if (7 == z)
      vdw[i] = 1.82;
    else if (8 == z)
      vdw[i] = 1.74;
    else if (9 == z)
      vdw[i] = 1.65;
    else if (16 == z)
      vdw[i] = 2.11;
    else if (17 == z)
      vdw[i] = 2.03;
    else if (35 == z)
      vdw[i] = 2.18;
    else if (53 == z)
      vdw[i] = 2.32;
    else if (15 == z)
      vdw[i] = 2.05;
    else if (5 == z)
      vdw[i] = 1.98;
    else if (34 == z)
      vdw[i] = mantina_selenium;
    else
    {
      cerr << "assign_molvol_vdw_radii::unknown atom type, atom " << i << " atomic number " << z << endl;
      vdw[i] = 0.0;
      rc = 0;
    }
  }

  return rc;
}

static int
assign_sybyl63_vdw_radii (Molecule & m, 
                          vdw_radius_t * vdw)
{
  int matoms = m.natoms ();

  int rc = 1;      // will be set to 0 if anything is inclassified

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number (i);

    if (1 == z)
      vdw[i] = 1.08;
    else if (6 == z)
      vdw[i] = 1.53;
    else if (7 == z)
      vdw[i] = 1.45;
    else if (8 == z)
      vdw[i] = 1.36;
    else if (9 == z)
      vdw[i] = 1.30;
    else if (15 == z)
      vdw[i] = 1.75;
    else if (16 == z)
      vdw[i] = 1.70;
    else if (17 == z)
      vdw[i] = 1.65;
    else if (35 == z)
      vdw[i] = 1.80;
    else if (53 == z)
      vdw[i] = 2.05;
    else if (5 == z)
      vdw[i] = mantina_boron;
    else if (34 == z)
      vdw[i] = mantina_selenium;
    else
    {
      cerr << "assign_sybyl63_vdw_radii:unknown atom type, atom " << i << " atomic number " << z << endl;
      vdw[i] = 0.0;
      rc = 0;
    }
  }

  return rc;
}

int
assign_vdw_radii (Molecule & m,
                  int vdw_type,
                  vdw_radius_t * vdw)
{
  if (IW_VDW_SHRAKE_AND_RUPLEY == vdw_type)
    return assign_shrake_and_rupley_vdw_radii(m, vdw);
  else if (IW_VDW_SAVOL == vdw_type)
    return assign_savol_vdw_radii(m, vdw);
  else if (IW_VDW_MOLVOL == vdw_type)
    return assign_molvol_vdw_radii(m, vdw);
  else if (IW_VDW_SYBYL63 == vdw_type)
    return assign_sybyl63_vdw_radii(m, vdw);
  else if (IW_VDW_WIKI == vdw_type)
    return assign_wiki_vdw_radii(m, vdw);

  cerr << "What kind of vdw type is this " << vdw_type << endl;

  return 0;
}

int
display_standard_vdw_radius_types (std::ostream & os, char flag, int full_details)
{
  os << "  -" << flag << " <type>      specify Van der Waals radius type\n";

  if (full_details)
  {
    os << "  -" << flag << " savol   Savol radii\n";
    os << "  -" << flag << " shrake  Shrake and Rupley\n";
    os << "  -" << flag << " molvol  MOLVOL radii\n";
    os << "  -" << flag << " sybyl63 Sybyl  radii\n";
    os << "  -" << flag << " wiki    wiki   radii (Dec 2005)\n";
  }

  return os.good();
}

int
set_default_van_der_waals_radius_type (Command_Line & cl,
                                       char flag,
                                       int & vdw_type,
                                       int verbose)
{
  IWString v = cl.string_value(flag);

  v.to_lowercase();

  if ("help" == v)
  {
    display_standard_vdw_radius_types(cerr, flag, 1);
    exit(verbose);
  }

  if ("savol" == v)
  {
    vdw_type = IW_VDW_SAVOL;
    if (verbose)
      cerr << "Will use Savol Van der Waals radii\n";
  }
  else if ("shrake" == v)
  {
    vdw_type = IW_VDW_SHRAKE_AND_RUPLEY;
    if (verbose)
      cerr << "Will use Shrake and Rupley Van der Waals radii\n";
  }
  else if ("molvol" == v)
  {
    vdw_type = IW_VDW_MOLVOL;
    if (verbose)
      cerr << "Will use MOLVOL Van der Waals radii\n";
  }
  else if ("sybyl63" == v)
  {
    vdw_type = IW_VDW_SYBYL63;

    if (verbose)
      cerr << "Will use Sybyl-6.3 Van der Waals radii\n";
  }
  else if ("wiki" == v)
  {
    vdw_type = IW_VDW_WIKI;

    if (verbose)
      cerr << "Will use Wiki (Dec 2005) Van der Waals radii\n";
  }
  else
  {
    cerr << "Unrecognised VDW radius type specifier '" << v << "'\n";
    cerr << "Choose one of these....\n";
    display_standard_vdw_radius_types(cerr, flag, 1);
    return 0;
  }

  return 1;
}
