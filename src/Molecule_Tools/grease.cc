/*
  Computes the grease demerit
*/

#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/smiles.h"

using std::cerr;
using std::cout;
using std::endl;

static Accumulator<float> hratio;
static Accumulator_Int<int>   largest_hydrophobic_size;

static Molecule_Output_Object molecules_rejected;
static Molecule_Output_Object molecules_not_rejected;

static int reduce_to_largest_fregmant = 0;

/*
  We can make decisions on presence or absence of formal
  charges.
*/

static Charge_Assigner charge_assigner;

/*
  We can run in a simple mode in which all we check is presence of
  only a single hydrophillic group.

  If the molecule has only a single group of hydrophillic atoms,
  and there are single_hydrophillic_max_size or fewer atoms in
  that group, and if the fraction of hydrophillic atoms in the
  molecule is less than single_hydrophillic_ratio, then reject
  the molecule.

  A typical usage might be to have single_hydrophillic_max_size == 2
  and single_hydrophillic_ratio == 0.2
  which would reject any molecule with 10 atoms and a single
  hydrophillic group of size two atoms.
*/

static int single_hydrophillic_max_size = 0;
static float single_hydrophillic_ratio = 0.0;

/*
  Jan 98, Steve wanted to see print-outs of the ratio1
*/

static Molecule_Output_Object shr_test_output;

static int molecules_with_formal_charges = 0;

/*
  For 3rd party processing we want to be able to write the inigial
  smiles to cout
*/

static int write_uncharged_smiles = 0;

/*
  May 2000. Richard Lewis wanted to just have the rejection annotated
  onto the output stream
*/

static int write_all_molecules_uncharged = 0;

const char *prog_name = nullptr;

void
usage (int rc = 1)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Programme for making the \"grease\" demerit\n";
  cerr << "usage: " << prog_name << " options file1 file2 file3 ....\n";
  cerr << "  -q <file>      specify one or more substructure queries\n";
  cerr << "  -m <number>    specify lower bound for size of grease group\n";
  cerr << "  -f <ratio>     apply simple test only (single hydrophillic section)\n";
  cerr << "  -z <min size>  specify max size for single hydrophillic section test\n";
  cerr << "  -s             collect statistics about grease qualities\n";
  cerr << "  -B <name>      specify stem for output of rejected molecules\n";
  cerr << "  -G <name>      specify stem for output of non rejected (good) molecules\n";
  cerr << "  -c             any molecule with a formal charge is NOT hydrophobic\n";
  cerr << "  -d <number>    add <number> to the hydrophillic side if formal charge present\n";
  cerr << "  -u             write uncharged smiles to stdout (for 3rd party processing)\n";
  cerr << "  -l             strip to largest fragment\n";
  display_standard_charge_assigner_options (cerr, 'N');
  cerr << "  -E <symbol>    standard element options, enter '-E help' for info\n";
  display_standard_aromaticity_options (cerr);
  cerr << "  -i ...         specify input type\n";
  cerr << "  -o  ...        specify output type\n";
  cerr << "  -v             verbose output\n";
// clang-format on

  exit (rc);
}

static int verbose = 0;

/*
  The smallest group of C atoms which is considered grease
*/

static int min_grease_size = 2;

/*
  We demerit based on the ratio of hydrophillic to phydrophobic atoms.
*/

static float hydrophillic_hydrophobic_ratio = 0.33;

static IWString stem_for_rejected;
static IWString stem_for_not_rejected;

static int collect_statistics = 0;

static int molecules_with_no_hydrophillic_atoms = 0;

/*
  One rule is that any molecule with a formal charge (we have a charge_assigner)
  is not considered hydrophobic
*/

static int formally_charged_molecules_are_hydrophillic = 0;

/*
  Alternatively, we can enlarge the size of any hydrophillic section if it
  contains a formally charged atom
*/

static int enlarge_formally_charged_hydrophillic_section = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static IWString demerit_string = "G";

static int
append_demerit_value (Molecule * m,
                      int demerit)
{
  if (0 == demerit)
    return 1;

  IWString cd;
  cd << m->name () << ' ' << demerit_string << '(' << demerit << ')';

  IWString new_name;

  new_name = cd;

  m->set_name (new_name);

  return 1;
}

/*
*/

static int
update_grease_id (int nmatches, 
                  const Substructure_Results & sresults,
                  int * grease_id)
{
  for (int i = 0; i < nmatches; i++)
  {
    const Set_of_Atoms * p = sresults.embedding (i);

    int np = p->number_elements ();
    for (int j = 0; j < np; j++)
    {
      atom_number_t a = p->item (j);
      grease_id[a]++;
    }
  }

  return 1;
}

class Grease : public Set_of_Atoms
{
  friend
    std::ostream & operator << (std::ostream &, const Grease &);

  private:
    int _non_ring_atoms;
    int _ring_atoms;
    int _hydrogens;

// To ensure we don't doubly count connections, keep track of all such atoms

  resizable_array<atom_number_t> _connections;

  public:
    Grease ();

    int ok () const;

    int identify_grease (Molecule & m, const atomic_number_t * z,
                         atom_number_t a, int id, int * grease_id);

    int ring_atoms () const { return _ring_atoms;}
    int non_ring_atoms () const { return _non_ring_atoms;}
    int ncon () const { return _connections.number_elements ();}
    int hydrogens () const { return _hydrogens;}
};

template class resizable_array_p<Grease>;
template class resizable_array_base<Grease *>;

Grease::Grease ()
{
  resize (8);    // 8 is some arbitrary number

  _connections.resize (8);
  _non_ring_atoms = 0;
  _ring_atoms = 0;
  _hydrogens = 0;

  return;
}

int
Grease::ok () const
{
  return 1;
}

#define IS_GREASY_ATOMIC_NUMBER(z) (6 == (z) || 9 == (z) || 17 == (z) || 35 == (z))

static int
is_greasy_atom(Molecule & m,
               const atomic_number_t * z,
               atom_number_t a)
{
  if (IS_GREASY_ATOMIC_NUMBER(z[a]))
    return 1;

// Jan 98, Steve Kaldor wants two connected S to be lipophillic

  if (16 == z[a] && 2 == m.ncon (a))
  {
    if (m.is_aromatic(a)) {
      return 0;
    }
    return 1;
  }


  return 0;
}

//#define DEBUG_IDENTIFY_GREASE

int
Grease::identify_grease (Molecule & m,
        const atomic_number_t * z,
        atom_number_t zatom,
        int id,
        int * grease_id)
{
  assert (0 == grease_id[zatom]);

#ifdef DEBUG_IDENTIFY_GREASE
  cerr << "Continuing grease expansion with atom " << zatom << endl;
#endif

  grease_id[zatom] = id;
  add (zatom);

  if (m.is_ring_atom (zatom))
    _ring_atoms++;
  else
    _non_ring_atoms++;

  Atom * a = const_cast<Atom *> (m.atomi (zatom));

  _hydrogens += a->implicit_hydrogens ();

  int rc = 1;

  int acon = m.ncon (zatom);
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type (zatom, i, j, bt);

    if (! is_greasy_atom (m, z, j))
    {
      _connections.add_if_not_already_present (j);
      continue;
    }

    if (grease_id[j])   // either already part of this group, or another group
    {
      if (id != grease_id[j])
        _connections.add_if_not_already_present (j);

      continue;
    }

    if (m.multiple_bond_to_heteroatom (j))
    {
      _connections.add_if_not_already_present (j);
      continue;
    }

    rc++;

    rc += identify_grease (m, z, j, id, grease_id);
  }

  return rc;
}

/*
  Counts of molecules rejected for various reasons
*/

static int ratio_reject = 0;

static int number_molecules_rejected = 0;

//#define DEBUG_HYDROPHILLIC_SECTION

static int
identify_hydrophillic_section (Molecule & m,
                               const atomic_number_t * z,
                               atom_number_t zatom,
                               int id,
                               int * grease_id,
                               Set_of_Atoms * hydrophillic)
{
#ifdef DEBUG_HYDROPHILLIC_SECTION
  cerr << "Continue hydrophillic section with atom " << a << endl;
#endif

  assert (grease_id[zatom] == 0);
  grease_id[zatom] = id;

  hydrophillic->add (zatom);

  const Atom * a = m.atomi (zatom);

  int rc = 1;

  if (enlarge_formally_charged_hydrophillic_section && a->formal_charge ())
    rc += enlarge_formally_charged_hydrophillic_section;

  int acon = a->ncon ();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t bt;
    a->other_and_type (zatom, i, j, bt);

    if (grease_id[j])
    {
      continue;
    }

    if (! is_greasy_atom (m, z, j) || m.multiple_bond_to_heteroatom (j))
    {
      rc += identify_hydrophillic_section (m, z, j, id, grease_id, hydrophillic);
    }
  }

  return rc;
}

std::ostream &
operator << (std::ostream & os, const Grease & g)
{
  os << "Grease with " << g._number_elements << " atoms, " <<
        g.ncon () << " connections. Atoms";
  for (int i = 0; i < g._number_elements; i++)
    os << ' ' << g._things[i];

  return os;
}

//#define DEBUG_GREASE

static int
grease (Molecule & m,
        int & demerit)
{
  demerit = 0;

  int matoms = m.natoms ();

  if (charge_assigner.active ())
    charge_assigner.process (m);

  if (formally_charged_molecules_are_hydrophillic && m.has_formal_charges ())
  {
    molecules_with_formal_charges++;
    if (verbose > 1)
      cerr << m.name () << " contains formal charges\n";
    return 1;
  }

  int * grease_id = new_int (matoms); std::unique_ptr<int[]> free_grease_id (grease_id);

  int nq = queries.number_elements ();
  for (int i = 0; i < nq; i++)
  {
    Substructure_Hit_Statistics * q = queries[i];

    Substructure_Results sresults;

    int nhits = q->substructure_search (m, sresults);
    if (nhits)
      update_grease_id (nhits, sresults, grease_id);
  }

  atomic_number_t * z = new atomic_number_t[matoms]; std::unique_ptr<atomic_number_t[]> free_z (z);

  m.atomic_numbers (z);

  resizable_array_p<Grease> grease;
  resizable_array_p<Set_of_Atoms> hydrophillic;

// At this stage, we should update the various counts for atoms identified by the search(es)
// Consider substructure stuff not really implemented yet....

  int rc = 0;

  int small_hydrophobic_sections = 0;  // those smaller than min size
  int hydrophillic_atoms = 0;
  int hydrophobic_atoms = 0;
  int atoms_in_hydrophobic_groups = 0;

  int id = 1;
  for (int i = 0; i < matoms; i++)
  {
    if (grease_id[i])
      continue;

    if (! is_greasy_atom (m, z, i) || m.multiple_bond_to_heteroatom (i))
    {
      Set_of_Atoms * tmp = new Set_of_Atoms;
      hydrophillic_atoms += identify_hydrophillic_section (m, z, i, id, grease_id, tmp);
      hydrophillic.add (tmp);
    }
    else
    {
      Grease * tmp = new Grease;

#ifdef DEBUG_IDENTIFY_GREASE
      cerr << "Begin grease identification with atom " << i << endl;
#endif

      rc += tmp->identify_grease (m, z, i, id, grease_id);
      hydrophobic_atoms += tmp->number_elements ();

      if (tmp->number_elements () < min_grease_size)
      {
        small_hydrophobic_sections++;
        delete tmp;
      }
      else
      {
        grease.add (tmp);
        atoms_in_hydrophobic_groups += tmp->number_elements ();
      }
    }

    id++;
  }

  int ng = grease.number_elements ();
  int nh = hydrophillic.number_elements ();

  if (verbose)
  {
    cerr << "Identified " << ng << " large grease groups,";
    if (small_hydrophobic_sections)
      cerr << " and " << small_hydrophobic_sections << " small regions,";
    cerr << ' ' << hydrophobic_atoms << " hydrophobic atoms\n";
    cerr << "And " << nh << " hydrophillic_sections containing " <<
            hydrophillic_atoms << " atoms\n";
  }

  if (verbose > 1)
  {
    for (int i = 0; i < ng; i++)
    {
      const Grease * g = grease[i];
      cerr << "Group " << i << ' ' << (*g) << endl;
    }
  }

  if (0 == ng)    // no grease in this molecule!
    return 1;

// Branch here if we are applying only the most loose test

  if (single_hydrophillic_ratio > 0.0)
  {
    if (0 == hydrophillic.number_elements ())
    {
      if (verbose)
        cerr << "No hydrophillic atoms\n";

      return 0;
    }

    if (hydrophillic.number_elements () > 1)    // more than one hydrophillic section
    {
      if (verbose > 1)
        cerr << "Not hydrophobic, " << hydrophillic.number_elements () << " hydrophillic groups\n";
      return 1;
    }

    const Set_of_Atoms * g = hydrophillic[0];
    if (single_hydrophillic_max_size > 0 && g->number_elements () > single_hydrophillic_max_size)
      return 1;

    if (shr_test_output.number_elements ())
    {
      IWString tmp = m.name ();
      tmp << " Ratio " << (float (g->number_elements ()) / float (matoms));
      m.set_name (tmp);

      shr_test_output.write (m);
    }

//  USE HYDROPHILLIC_ATOMS here because it includes the influence of formal
//  charges. Note that we are in the case of just one hydrophillic section,
//  so HYDROPHILLIC_ATOMS refers to the (possibly adjusted) score for that
//  single section

    if ( (float (hydrophillic_atoms) / float (matoms)) < single_hydrophillic_ratio)
    {
      if (verbose)
        cerr << "Fails single hydrophillic ratio test\n";

      return 0;
    }

    return 1;
  }

  int mkeep = 1;

  int lhs = 0;    // largest hydrophobic size
  for (int i = 0; i < ng; i++)
  {
    const Grease * g = grease[i];
    assert (g->ok ());

    int gne = g->number_elements ();

    if (gne > lhs)
      lhs = gne;

    if (gne > 6 && 1 == g->ncon ())
    {
      if (g->ring_atoms () <= 1)
        demerit += 100;
      else
        demerit += 7 * (gne - g->ring_atoms () / 2);

      if (verbose > 1)
        cerr << "Large (" << gne << ") singly connected grease, demerit = " << demerit << endl;

      continue;
    }

    if (gne > 6 && g->ncon () <= 2)
    {
      int tmp = gne - 2 * g->ncon ();
      if (tmp <= 0)
        demerit += 10;
      else
        demerit += 5 * tmp;

      if (verbose > 1)
        cerr << "More than 6 greasy atoms (" << gne << "), too few connections (" << g->ncon () << ") demerit = " << demerit << endl;

      continue;
    }

    if (gne - 2 * g->ncon () - g->non_ring_atoms () > 6)
    {
      int tmp = gne - 2 * g->ncon () - g->non_ring_atoms ();
      if (tmp <= 0)
        demerit += 10;
      else
        demerit += 5 * tmp;

      if (verbose > 1)
        cerr << "Too many greasy atoms: gne = " << gne << ", ncon = " << g->ncon () <<
                " non ring = " << g->non_ring_atoms () << " demerit = " << demerit << endl;
      continue;
    }

    if (6 == gne && 1 == g->ncon ())
    {
      demerit += 30;

      if (verbose)
        cerr << "Demerited for singly connected grease 6, demerit = " << demerit << endl;

      continue;
    }

    if (5 == gne && 1 == g->ncon ())
    {
      demerit += 20;

      if (verbose)
        cerr << "Demerited for singly connected grease 5, demerit = " << demerit << endl;

      continue;
    }

    if (4 == gne && 1 == g->ncon ())
    {
      demerit += 15;

      if (verbose)
        cerr << "Demerited for singly connected grease 4, demerit = " << demerit << endl;

      continue;
    }

    if (g->non_ring_atoms () >= 4)
    {
      demerit += 5 * (g->non_ring_atoms () - 4);

      if (verbose > 1)
        cerr << "Demerited for too many non ring atoms " << g->non_ring_atoms () << " demerit = " << demerit << endl;
    }

    if (gne >= 10)
    {
      int tmp = gne - (2 * g->ncon () + 1);
      if (tmp < 0)
        demerit += 5  * (gne - 10);
      else
        demerit += 8 * tmp;

      if (verbose > 1)
        cerr << "Demerited for very large grease " << gne << " ncon = " << g->ncon () << " demerit = " << demerit << endl;
    }
    else if (gne > 6 && g->ncon () < 2)
    {

      demerit += 8 * (gne - g->ncon () + 1);

      if (verbose > 1)
        cerr << "Demerited for grease too large " << gne << ", " <<
                g->ncon () << " connections, demerit = " << demerit << endl;
    }
    else if (gne >= 5)
    {
      int tmp = gne - (2 * g->ncon () + 1);
      if (g->ring_atoms ())
        tmp--;

      if (tmp <= 0)
        demerit += 5;
      else
        demerit += 5 * tmp;

      if (verbose > 1)
        cerr << "Demerited for grease = " << gne << " ncon = " << g->ncon () << ", demerit = " << demerit << endl;
    }
  }

  if (verbose > 1)
    cerr << "After examining individual grease groups, demerit = " << demerit << 
            " mkeep = " << mkeep << endl;

  if (collect_statistics)
    largest_hydrophobic_size.extra (lhs);

  if (0 == hydrophillic_atoms)
  {
    molecules_with_no_hydrophillic_atoms++;
    if (verbose > 1)
      cerr << "Rejected for no hydrophillic atoms\n";

    demerit += 100;

    return 0;
  }

// try to deal with large molecules with lots of small bits of grease

  if (ng > 4)
  {
    float tmp = 1.0 - float (ng - 4) / 10.0;
    demerit = int (demerit * tmp);
    if (verbose > 1)
      cerr << "Demerit adjusted for many groups " << ng << " demerit = " << demerit << endl;
  }

  float tmp = float (hydrophobic_atoms) / float (hydrophillic_atoms);
  if (collect_statistics)
    hratio.extra (tmp);

  if (nh > 2 * ng)
  {
    if (demerit >= 100)
      mkeep = 0;
    return mkeep;
  }

  if (tmp > 2.99 && 1 == nh)
  {
    ratio_reject++;
    if (verbose > 1)
      cerr << "Rejected for bad ratio " << tmp << endl;
    mkeep = 0;
  }
  else if (nh > ng + 3)
    ;
  else if ( (tmp = float (hydrophillic_atoms) / float (atoms_in_hydrophobic_groups)) < hydrophillic_hydrophobic_ratio)
  {
    if (verbose > 1)
      cerr << "Demerited for too few hydrophillic atoms " << tmp << endl;
    demerit += int (300.0 * (hydrophillic_hydrophobic_ratio - tmp));
  }

#ifdef DEBUG_GREASE
  cerr << "At end, demerit = " << demerit << " mkeep = " << mkeep << endl;
#endif

  if (demerit >= 100)
    mkeep = 0;

  return mkeep;
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fregmant)
    m.reduce_to_largest_fragment ();
  else if (1 != m.number_fragments ())
    cerr << "Warning, multi fragment molecule\n";

  m.compute_aromaticity ();
}

static int
grease (data_source_and_type<Molecule> & input)
{
  assert (input.good ());

  Molecule * m;
  while (nullptr != (m = input.next_molecule ()))
  {
    std::unique_ptr<Molecule> free_m (m);

    preprocess (*m);

    IWString smiles;
    if (write_uncharged_smiles || write_all_molecules_uncharged)
      smiles = m->smiles ();

    int demerit;
    int mkeep = grease (*m, demerit);

    if (write_all_molecules_uncharged)
    {
      cout << smiles << ' ' << m->name ();
      if (0 == mkeep)
        cout << " GREASY";
      cout << endl;
    }
    else if (0 == mkeep)
    {
      if (verbose)
        cerr << "Rejected\n";
      number_molecules_rejected++;

      if (molecules_rejected.active ())
        molecules_rejected.write (m);
    }
    else if (write_uncharged_smiles)
    {
      cout << smiles << ' ' << m->name () << endl;
    }
    else if (molecules_not_rejected.active ())
    {
      append_demerit_value (m, demerit);
      molecules_not_rejected.write (m);
    }

#ifdef USE_IWMALLOC
    check_all_malloced (stderr);
#endif
  }

  return 0;
}

static int
grease (const char * fname,
             FileType input_type)
{
  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name (fname);
    if (input_type == FILE_TYPE_INVALID)
    {
      cerr << "Cannot discern input file type from '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_verbose (verbose);

  return grease (input);
}

#include "Foundational/cmdline/cmdline.h"

int
grease (int argc, char **argv)
{
  Command_Line cl (argc, argv, "cN:f:z:sm:B:G:q:vE:A:i:o:lT:uU");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised options encountered\n";
    usage (4);
  }

  verbose = cl.option_count ('v');

  FileType input_type = FILE_TYPE_INVALID;
  if (! cl.option_present ('i'))
    ;
  else if (! process_input_type (cl, input_type))
    usage (1);

  if (! process_elements (cl))
    usage (3);

  if (! process_standard_aromaticity_options (cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage (8);
  }

  if (cl.option_present ('l'))
  {
    reduce_to_largest_fregmant = 1;
    if (verbose)
      cerr << "Multi fragment molecules will be trimmed to their largest form\n";
  }

  if (cl.option_present ('B'))
  {
    if (! cl.option_present ('o'))
    {
      cerr << "The -B option requires the -o option\n";
      usage (3);
    }

    cl.value ('B', stem_for_rejected);
    if (verbose)
      cerr << "Rejected molecules will be written to file stem '" << stem_for_rejected << "'\n";
  }

  if (cl.option_present ('G'))
  {
    if (! cl.option_present ('o'))
    {
      cerr << "The -G option requires the -o option\n";
      usage (4);
    }

    cl.value ('G', stem_for_not_rejected);
    if (verbose)
      cerr << "Non rejected molecules will be written to file stem '" << stem_for_not_rejected << "'\n";
  }

  if (cl.option_present ('u'))
  {
    if (cl.option_present ('G'))
    {
      cerr << "The -u and -G options are mutually exclusive\n";
      usage (7);
    }

    write_uncharged_smiles = 1;
    if (verbose)
      cerr << "Uncharged smiles written to cout\n";
  }

  if (cl.option_present ('U'))
  {
    if (cl.option_present ('u') || cl.option_present ('B') || cl.option_present ('G'))
    {
      cerr << "The -U option cannot be used with any of the -u, -B or -G options\n";
      usage (18);
    }

    write_all_molecules_uncharged = 1;
    if (verbose)
      cerr << "All molecules written to stdout\n";
  }

  if (cl.option_present ('o'))
  {
    if (0 == stem_for_rejected.length () && 0 == stem_for_not_rejected.length ())
    {
      cerr << "If you use the -o option, you must specify one of -G or -B\n";
      usage (8);
    }

    if (stem_for_rejected.length ())
    {
      if (! molecules_rejected.determine_output_types (cl))
      {
        cerr << "Could not determine output type(s), -o option\n";
        usage (4);
      }

      for (int i = 0; i < cl.number_elements (); i++)
      {
        if (molecules_rejected.would_use_name (stem_for_rejected, cl[i]))
        {
          cerr << "Cannot overwrite input file '" << cl[i] << "'\n";
          return 1;
        }
      }

      molecules_rejected.new_stem (stem_for_rejected);
    }

    if (stem_for_not_rejected.length ())
    {
      if (! molecules_not_rejected.determine_output_types (cl))
      {
        cerr << "Could not determine output type(s), -o option\n";
        usage (4);
      }

      for (int i = 0; i < cl.number_elements (); i++)
      {
        if (molecules_not_rejected.would_use_name (stem_for_rejected, cl[i]))
        {
          cerr << "Cannot overwrite input file '" << cl[i] << "'\n";
          return 1;
        }
      }
      molecules_not_rejected.new_stem (stem_for_not_rejected);
    }
  }

  if (cl.option_present ('m'))
  {
    if (! cl.value ('m', min_grease_size) || min_grease_size < 1)
    {
      cerr << "The -m switch requires a positive integer\n";
      usage (3);
    }

    if (verbose)
      cerr << "The smallest grease group is size " << min_grease_size << endl;
  }

  if (cl.option_present ('f'))
  {
    if (! cl.value ('f', single_hydrophillic_ratio) || single_hydrophillic_ratio <= 0.0 || single_hydrophillic_ratio >= 1.0)
    {
      cerr << "The -f option must be followed by a fraction between 0.0 and 1.0\n";
      usage (16);
    }

    if (cl.option_present ('z'))
    {
      if (! cl.value ('z', single_hydrophillic_max_size) || single_hydrophillic_max_size < 1)
      {
        cerr << "The -z option requires a whole positive number\n";
        usage (19);
      }
    }

    if (verbose)
    {
      cerr << "Will apply the single hydrophillic section test, ratio = " << single_hydrophillic_ratio;
      if (single_hydrophillic_max_size)
        cerr << " min size " << single_hydrophillic_max_size;
      cerr << endl;
    }
  }
  else if (cl.option_present ('z'))
  {
    cerr << "You can only specify a min size for the single hydrophillic section\n";
    cerr << "test if you also specify a ratio via the -f option\n";
    usage (17);
  }

  if (cl.option_present ('s'))
  {
    collect_statistics = 1;
    if (verbose)
      cerr << "Statistics about grease will be collected\n";
  }

  if (cl.option_present ('q'))
  {
    if (! process_queries (cl, queries, verbose, 'q'))
    {
      cerr << "Cannot process queries\n";
      return 8;
    }
  }

  if (cl.option_present ('N'))
  {
    if (! charge_assigner.construct_from_command_line (cl, verbose, 'N'))
    {
      cerr << "Cannot determine charge assigner from command line\n";
      usage (77);
    }
  }

  if (cl.option_present ('c'))
  {
    formally_charged_molecules_are_hydrophillic = 1;
    if (verbose)
      cerr << "Any molecule with a formal charge is NOT hydrophobic\n";

    if (! charge_assigner.active ())
      cerr << "Warning, -c option present, but charge assigner not active\n";
  }

  if (cl.option_present ('d'))
  {
    if (formally_charged_molecules_are_hydrophillic)
    {
      cerr << "The -c and -d options are mutually exclusive\n";
      usage (15);
    }

    if (! cl.value ('d', enlarge_formally_charged_hydrophillic_section))
    {
      cerr << "The -d option must be followed by a positive whole number\n";
      usage (19);
    }

    if (verbose)
      cerr << "Hydrophillic sections with formally charged atoms increased by " << enlarge_formally_charged_hydrophillic_section << " atoms\n";
  }

  if (cl.option_present ('T'))
  {
    IWString fname;
    cl.value ('T', fname);

    shr_test_output.add_output_type (FILE_TYPE_SMI);
    if (! shr_test_output.new_stem (fname))
    {
      cerr << "cannot open test output file '" << fname << "'\n";
      return 87;
    }

    if (verbose)
      cerr << "Will write single hydrophillic ratio values to '" << fname << "'\n";
  }

// By default we echo any chiral info present in the input

// There must be at least one token remaining on the command line

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (4);
  }

  for (int i = 0; i < cl.number_elements (); i++)
  {
    (void) grease (cl[i], input_type);
  }

  if (verbose)
  {
    int nq = queries.number_elements ();
    for (int i = 0; i < nq; i++)
    {
      const Substructure_Hit_Statistics * q = queries[i];
      q->report (cerr, verbose);
    }

    if (molecules_with_formal_charges)
      cerr << molecules_with_formal_charges << " molecules with formal charges not ignored\n";
  }

  if (collect_statistics)
  {
    cerr << "Hratio statistics for " << hratio.n () << " molecules\n";
    cerr << "Values between " << hratio.minval () << " and " << hratio.maxval ();
    if (hratio.n () > 1)
      cerr << " ave = " << hratio.average () << " var = " << hratio.variance ();
    cerr << endl;
    if (molecules_with_no_hydrophillic_atoms)
      cerr << molecules_with_no_hydrophillic_atoms << " molecules had zero hydrophillic atoms\n";

    if (number_molecules_rejected)
      cerr << number_molecules_rejected << " molecules rejected\n";

    cerr << "Largest hydrophobic sizes between " << largest_hydrophobic_size.minval () <<
            " and " << largest_hydrophobic_size.maxval ();
    if (largest_hydrophobic_size.n () > 1)
      cerr << " ave " << largest_hydrophobic_size.average () << " var = " << largest_hydrophobic_size.variance ();
    cerr << endl;

    if (charge_assigner.active ())
      charge_assigner.report (cerr);
  }

  return 0;
}

int 
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = grease (argc, argv);

  return rc;
}
