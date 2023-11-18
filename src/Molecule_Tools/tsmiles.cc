#include <iostream>
#include <memory>
#include <random>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int abort_on_error = 0;
static int stop_processing_molecule_on_first_error = 0;

static int test_subset_unique_smiles = 0;

static int molecules_with_errors = 0;
static int total_error_count = 0;

static int smiles_interpretation_errors = 0;
static int aromatic_smiles_interpretation_errors = 0;
static int amw_errors = 0;
static int unique_smiles_errors = 0;

static int test_swap_atoms = 0;
static int atom_swap_errors = 0;

static int permutations_per_molecule = 10;
static int permutations_performed = 0;

static int remove_non_chiral_chiral_centres = 0;

static int unset_unnecessary_implicit_hydrogens_known_values = 0;

static int write_sssr_upon_failure = 0;

static int test_build_from_aromatic_form = 0;

static int test_include_chiral_centres_in_smiles = 0;

static int test_same_ring_stuff = 0;

static int chiral_smiles_failures = 0;

static Report_Progress report_progress;

static int molecules_with_too_many_rings = 0;

static IWString_and_File_Descriptor stream_for_failed_molecules;

static int max_rings = 99;

static int destroy_aromaticity = 0;

static int molecules_ok_without_aromaticity = 0;

static int remove_all_chiral_centres = 0;

static std::random_device rd;

void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "usage: " << prog_name << " options file1 file2 file3....\n";
  cerr << "  -p <number>    do <number> permutations per molecule\n";
  cerr << "  -a             test the assignment operator\n";
  cerr << "  -b             test subset unique smiles\n";
  cerr << "  -c <0,1>       include chiral info in smiles (default 1)\n";
  cerr << "  -t <0,1>       include cis trans info in smiles (default 1)\n";
  cerr << "  -h             test set_include_chiral_info_in_smiles\n";
  cerr << "  -w <number>    do tests involving swapping atoms\n";
  cerr << "  -m <number>    test building from aromatic form\n";
  cerr << "  -x             terminate all processing on any error\n";
  cerr << "  -j             stop processing any molecule once it has an error\n";
  cerr << "  -u             remove marked chiral centres that aren't chiral\n";
  cerr << "  -y             remove any un-necessary square brackets (will get radicals wrong)\n";
  cerr << "  -q             test in_same_ring and in_same_ring_system functions\n";
  cerr << "  -R             write SSSR upon failure\n";
  cerr << "  -f             upon failure, destroy aromaticity and see if a non aromatic form would work\n";
  cerr << "  -r <number>    report progress every <n> molecules\n";
  cerr << "  -W             remove all chirality as the molecule is read\n";
  cerr << "  -s <number>    random number seed\n";
  cerr << "  -F <fname>     write smiles for failing molecules to <fname>\n";
  cerr << "  -i <type>      specify input type\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -k <nrings>    discard molecules with <nrings> or more (def 99)\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_smiles_options(cerr);
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int molecules_read = 0;

static int verbose = 0;

static int test_assignment = 0;

static int
do_destroy_aromaticity (Molecule & m)
{
  Set_of_Atoms double_bonds_to_be_removed;

  m.compute_aromaticity_if_needed();

  const int nb = m.nedges();

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = m.bondi(i);

    if (! b->is_aromatic())
      continue;

    double_bonds_to_be_removed.add(b->a1());
    double_bonds_to_be_removed.add(b->a2());
  }

  const int n = double_bonds_to_be_removed.size();

  if (0 == n)
    return 0;

  for (auto i = 0; i < n; i += 2)
  {
    m.set_bond_type_between_atoms(double_bonds_to_be_removed[i], double_bonds_to_be_removed[i+1], SINGLE_BOND);
  }

  m.append_to_name(" XAROM");

  return 1;
}

static void
do_remove_non_chiral_chiral_centres (Molecule & m)
{
  int nchiral = m.chiral_centres();
  if (0 == nchiral)
    return;

  int matoms = m.natoms();

//cerr << "Has " << nchiral << " chiral centres and " << matoms << " atoms\n";
  int centres_processed = 0;

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (nullptr == m.chiral_centre_at_atom(i))
      continue;

    if (! is_actually_chiral(m, i))
    {
      m.remove_chiral_centre_at_atom(i);
      rc++;
    }

    centres_processed++;

    if (centres_processed == nchiral)
      break;
  }

  if (verbose > 1 && rc)
    cerr << "removed " << rc << " false chiral centres\n";

  return;
}

/*
  Identify all connected atoms.
*/

static int
either_side (const Atom ** atom,
             atom_number_t a1,
             int n1,
             int * wside)
{
  const Atom * a = atom[a1];

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a1, i);
    if (n1 == wside[j])     // closing a ring or turning back on ourselves within our subset
      continue;

    wside[j] = n1;
    if (! either_side(atom, j, n1, wside))
      return 0;
  }

  return 1;
}

/*
  Classify a set of atoms according to which side of a given bond
  Those connected to A1 will be assigned N1, those connected via A2 are
  assigned N2
*/

static int
either_side (const Atom ** atom,
             atom_number_t a1,
             int n1,
             atom_number_t a2,
             int n2,
             int * wside)
{
  wside[a1] = n1;
  wside[a2] = n1;     // temporarily

  if (! either_side(atom, a1, n1, wside))
    return 0;

  wside[a1] = n2;     // temporarily
  wside[a2] = n2;

  if (! either_side(atom, a2, n2, wside))
    return 0;

  wside[a1] = n1;

  return 1;
}

/*
  Break the molecule at chain single bonds and form components
*/

static int
do_test_subset_unique_smiles (Molecule & m,
                              const Atom ** atom,
                              int * subset)
{
  (void) m.ring_membership();

  int matoms = m.natoms();

  int nb = m.nedges();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (b->nrings())
      continue;

    if (! b->is_single_bond())
      continue;

//  Divide the molecule across the bond

    (void) either_side(atom, b->a1(), 0, b->a2(), 1, subset);

    Smiles_Information smi_info;
    IWString usmi11 = m.unique_smiles(smi_info, subset);
    for (int i = 0; i < matoms; i++)
    {
      if (0 == subset[i])
        subset[i] = 1;
      else
        subset[i] = 0;

//    cerr << "Atom " << i << " subset " << subset[i] << endl;
    }

    smi_info.invalidate();

    IWString usmi12 = m.unique_smiles(smi_info, subset);

    Molecule tmp = m;

    tmp.remove_bond_between_atoms(b->a1(), b->a2());

    assert (tmp.number_fragments() > 1);

//  cerr << "Creating components\n";

    resizable_array_p<Molecule> components;
    (void) tmp.create_components(components);

//  cerr << "Created " << components.number_elements() << " components with " << components[0]->natoms() << " and " << components[1]->natoms() << " atoms\n";

    IWString usmi21 = components[0]->unique_smiles();
//  cerr << "Computing unique smiles for component with " << components[1]->natoms() << " atoms\n";
    IWString usmi22 = components[1]->unique_smiles();
//  cerr << "Unique smiles of components complete\n";
  }

  return 1;
}

static int
do_test_subset_unique_smiles (Molecule & m)
{
  const Atom ** atoms = new const Atom *[m.natoms()]; std::unique_ptr<const Atom *[]> free_atoms(atoms);
  
  m.atoms(atoms);

  int * tmp = new int[m.natoms()]; std::unique_ptr<int[]> free_tmp(tmp);

  int rc = do_test_subset_unique_smiles(m, atoms, tmp);

//delete atoms;

  return rc;
}

static int
do_test_assignment (Molecule & m, const IWString & usmi)
{
  molecular_weight_t mamw;
  if (! m.contains_non_periodic_table_elements())
    mamw = m.molecular_weight_ignore_isotopes();
  else
    mamw = 0.0;

  Molecule tmp(m);

  const IWString & tmpusmi = tmp.unique_smiles();

  if (tmpusmi != usmi)
  {
    cerr << "Yipes, unique smiles mismatch on constructor assignement\n";
    cerr << usmi << ' ' << m.name() << '\n';
    cerr << tmpusmi << ' ' << m.name() << '\n';
  }

  if (mamw > 0.0)
  {
    molecular_weight_t tmpamw = tmp.molecular_weight_ignore_isotopes();

    if (fabs(tmpamw - mamw) > 0.1)
    {
      cerr << "Molecular weight mismatch on constructor assignment, " << mamw << " vs " << tmpamw << '\n';
      return 0;
    }
  }

  tmp = m;

  if (tmp.unique_smiles() != usmi)
  {
    cerr << "Yipes, unique smiles mismatch on assignment operator\n";
    cerr << usmi << ' ' << m.name() << '\n';
    cerr << tmp.unique_smiles() << ' ' << m.name() << '\n';
  }

  return 1;
}

static int
report_canonical_ranking (std::ostream & os,
                          Molecule & m)
{
  os << "Canonical ordering:\n";
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    os << "Atom " << i << " (" << m.atomic_symbol(i) << ") canonical number " << 
          m.canonical_rank(i) << " symmetry class " << m.symmetry_class(i) << " bonds";
    int acon = m.ncon(i);
    for (int j = 0; j < acon; j++)
    {
      os << ' ' << m.other(i, j);
    }

    os << '\n';
  }

  if (write_sssr_upon_failure)
  {
    int nr = m.nrings();

    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);
      os << " ring " << i << ':';
      for (int j = 0; j < ri->number_elements(); j++)
      {
        atom_number_t k = ri->item(j);
        os << ' ' << k << " (" << m.atomic_symbol(k) << ')';
      }
      os << '\n';
    }
  }

  return os.good();
}

static int
report_aromatic_atom_count_if_different (Molecule & m1,
                                         Molecule & m2,
                                         std::ostream & os)
{
  int matoms = m1.natoms();

  if (matoms != m2.natoms())
  {
    cerr << "report_aromatic_atom_count_if_different:atom count mismatch " << matoms << " vs " << m2.natoms() << '\n';
    return 0;
  }

  if (m1.aromatic_atom_count() == m2.aromatic_atom_count())
    return os.good();

  os << "Aromatic atom count mismatch " << m1.aromatic_atom_count() << " vs " << m2.aromatic_atom_count() << '\n';

  return os.good();
}

static int
test_same_ring_system (Molecule & m,
                       const Ring & r1,
                       const Ring & r2)
{
  assert (r1.fused_system_identifier() == r2.fused_system_identifier());

  int r1size = r1.number_elements();
  int r2size = r2.number_elements();

  for (int i = 0; i < r1size; i++)
  {
    atom_number_t ai = r1[i];

    for (int j = 0; j < r2size; j++)
    {
      atom_number_t aj = r2[j];

      if (ai == aj)
        continue;

      if (! m.in_same_ring_system(ai, aj))
      {
        cerr << "In '" << m.name() << "' not in same ring system, atoms " << ai << " and " << r2[j] << '\n';
        cerr << r1 << '\n';
        cerr << r2 << '\n';
        return 0;
      }
    }
  }

  return 1;
}

static int
test_same_ring(Molecule & m,
               const Ring & r)
{
  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t ai = r[i];

    for (int j = i + 1; j < ring_size; j++)
    {
      if (! m.in_same_ring(ai, r[j]))
      {
        cerr << "In same ring failure '" << m.name() << "', atoms " << ai << " and " << r[j] << '\n';
        return 0;
      }
    }
  }

  return 1;
}

static int
test_same_ring_and_same_ring_system(Molecule & m)
{
  const int nr = m.nrings();

  if (0 == nr)
    return 1;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! test_same_ring(m, *ri))
      return 0;

    if (! ri->is_fused())
      continue;

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi(j);

      if (rj->fused_system_identifier() != ri->fused_system_identifier())
        continue;

      if (! test_same_ring_system(m, *ri, *rj))
        return 0;
    }
  }

  return 1;
}


static int
test_include_chirality_in_smiles(Molecule & m)
{
  if (0 == m.chiral_centres())
    return 1;

  set_include_chiral_info_in_smiles(0);

  m.invalidate_canonical_ordering_information();

  IWString usmi = m.unique_smiles();

  m.remove_all_chiral_centres();

  m.invalidate_smiles();

  IWString s2 = m.unique_smiles();

//cerr << "Compare '" << usmi << "'\n";
//cerr << "vs      '" << s2 << "'\n";

  if (usmi == s2)
    return 1;

  cerr << "set_include_chirality_in_smiles failure, " << m.name() << "\n";
  cerr << "smiles w/o chirality " << usmi << "\n";
  cerr << "smiles from achiral  " << s2 << "\n";

  chiral_smiles_failures++;

  return 0;
}
static int
choose_two_atom_numbers(int matoms,
                        atom_number_t & a1,
                        atom_number_t & a2)
{
  assert(matoms > 1);

  std::uniform_int_distribution<int> um(0, matoms - 1);

  a1 = um(rd);

  do
  {
    a2 = um(rd);
  }
  while(a2 == a1);

  return 1;   // will never come here
}


#ifdef NOT_NEEDED_NOW
static void
report_pi_electrons(Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    int pi;
    m.pi_electrons(i, pi);
    cerr << "Atom " << i << " " << m.smarts_equivalent_for_atom(i) << " has " << pi << " pi electrons\n";
  }
}
#endif

#ifdef MAYBE_IMPLEMENT_THIS_SOMETIME
class Counters_This_Molecle
{
  public:
    const IWString & _usmi;
    const molecular_weight_t _amw;

    int _smiles_interpretation_errors;
    int _aromatic_smiles_interpretation_errors;

  public:
    Counters_This_Molecle(const IWString & s, const molecular_weight_t & mw);

    void smiles_interpretation_error() { _smiles_interpretation_errors++;}
    void aromatic_smiles_interpretation_error() { _aromatic_smiles_interpretation_errors++;}
};

Counters_This_Molecle::Counters_This_Molecle(const IWString & s,
                                const molecular_weight_t & mw) : _usmi(s), _amw(mw)
{
  _smiles_interpretation_errors = 0;
  _aromatic_smiles_interpretation_errors = 0;
}

static int
tsmiles(Molecule & m,
        Counters_This_Molecle & c,
        const IWString & rsmi)
{
  Molecule x;
  if (! x.build_from_smiles(rsmi))
  {
    c.smiles_interpretation_error();
  }
}

static int
tsmiles(Molecule & m,
        const molecular_weight_t amw)
{
  const IWString & usmi = m.unique_smiles();

  Counters_This_Molecle c(m.unique_smiles(), amw);

  for (int i = 0; i < permutations_per_molecule; ++i)
  {
    permutations_performed++;

    IWString rsmi = m.random_smiles();
 
    if (verbose > 2)
      cerr << "Permutation " << i << " " << rsmi << "\n";

    if (tsmiles(m, c, rsmi))
      continue;
  }
}
#endif

/*
  We return the number of errors we encounter
*/

static int
tsmiles(Molecule & m,
        const IWString & usmiles,
        const molecular_weight_t amw)
{
  int errors_this_molecule = 0;

  IWString reason_for_failure;

//report_pi_electrons(m);

  for (int i = 0; i < permutations_per_molecule; i++)
  {
    if (errors_this_molecule && (abort_on_error || stop_processing_molecule_on_first_error))
      break;

    permutations_performed++;

    IWString rsmi = m.random_smiles();
 
    if (verbose > 2)
      cerr << "Permutation " << i << " " << rsmi << "\n";

    Molecule x;
    if (! x.build_from_smiles(rsmi.chars()))      // this should never happen
    {
      cerr << "Error: molecule " << molecules_read << " '" << m.name() << "'\n";
      cerr << "Permutation " << i << " rsmiles " << rsmi << "\n";
      cerr << "Cannot interpret smiles\n";
      smiles_interpretation_errors++;
      errors_this_molecule++;
      reason_for_failure << " Smiles Interpretation";
      continue;
    }

//  report_pi_electrons(x);

    IWString xusmi = x.unique_smiles();
    if (usmiles != xusmi)
    {
      cerr << "Error: molecule " << molecules_read << " '" << m.name() << "'\n";
      cerr << "Permutation " << i << " rsmiles " << rsmi << "\n";
      cerr << "Unique smiles mis-match\n";
      cerr << "Original " << usmiles << "\n";
      cerr << "Now      " << xusmi << "\n";
      report_aromatic_atom_count_if_different(m, x, cerr);
      if (verbose <= 1)
      {
        cerr << "Original molecule\n";
        report_canonical_ranking(cerr, m);
      }
      cerr << "Random variant\n";
      report_canonical_ranking(cerr, x);

      unique_smiles_errors++;
      errors_this_molecule++;
      reason_for_failure << " unique smiles mismatch";

      continue;
    }

    if (amw > 0.0)
    {
      molecular_weight_t xamw = x.molecular_weight_ignore_isotopes();
      if (fabs(amw - xamw) > 0.1)
      {
        cerr << "Error: molecule " << molecules_read << " '" << m.name() << "'\n";
        cerr << "Permutation " << i << " rsmiles " << rsmi << "\n";
        cerr << "AMW mismatch, originally " << amw << " now " << xamw << '\n';
        amw_errors++;
        errors_this_molecule++;
        reason_for_failure << " AMW mismatch";
        continue;
      }
    }
  }

  if (test_build_from_aromatic_form && m.contains_aromatic_atoms())
  {
    for (int i = 0; i < test_build_from_aromatic_form; i++)
    {
      m.compute_aromaticity_if_needed();
      set_include_aromaticity_in_smiles(1);
      IWString random_smiles = m.random_smiles();
      set_include_aromaticity_in_smiles(0);

//    cerr << "Random aromatic form " << i << ' ' << random_smiles << '\n';

      Molecule mtmp;
      if (! mtmp.build_from_smiles(random_smiles))
      {
        cerr << "Cannot parse aromatic smiles " << random_smiles << " '" << m.name() << "'\n";
        errors_this_molecule++;
        aromatic_smiles_interpretation_errors++;
        reason_for_failure << " aromatic smiles";
        if (stop_processing_molecule_on_first_error)
          break;
        continue;
      }

      IWString xusmi = mtmp.unique_smiles();
      if (usmiles != xusmi)
      {
        cerr << "Error: molecule " << molecules_read << " '" << m.name() << "'\n";
        cerr << "Permutation " << i << " rsmiles " << random_smiles << "\n";
        cerr << "Unique smiles mis-match (from aromatic)\n";
        cerr << "Original " << usmiles << "\n";
        cerr << "Random   " << random_smiles << "\n";
        cerr << "Now      " << xusmi << "\n";
        report_aromatic_atom_count_if_different(m, mtmp, cerr);

        unique_smiles_errors++;
        errors_this_molecule++;
        reason_for_failure << " aromatic smiles mismatch";
        if (stop_processing_molecule_on_first_error)
          break;
      }
    }
  }

  const int matoms = m.natoms();

  if (test_swap_atoms && matoms > 2)
  {
    Molecule x(m);
    atom_number_t a1, a2;   // scope here for efficiency

    for (int i = 0; i < test_swap_atoms; i++)
    {
//    cerr << "Testing atom swap " << i << '\n';

      choose_two_atom_numbers(matoms, a1, a2);
      x.swap_atoms(a1, a2);

      const IWString & xusmi = x.unique_smiles();

      if (xusmi != usmiles)
      {
        cerr << "Atom swap failure '" << m.name() << "'\n";
        cerr << "Parent unique smiles " << usmiles << "\n";
        cerr << "Atom swapped variant " << xusmi << "\n";
        atom_swap_errors++;
        errors_this_molecule++;
        reason_for_failure << " swap error";
        if (stop_processing_molecule_on_first_error)
          break;
      }
    }
  }

  if (test_include_chiral_centres_in_smiles)
  {
    if (! test_include_chirality_in_smiles(m))
    {
      errors_this_molecule++;
      reason_for_failure << " chirality";
    }
  }

  if (test_same_ring_stuff)
  {
    if (! test_same_ring_and_same_ring_system(m))
    {
      errors_this_molecule++;
      reason_for_failure << " same ring";
    }
  }

  if (errors_this_molecule)
  {
    total_error_count += errors_this_molecule;
    molecules_with_errors++;

    if (stream_for_failed_molecules.is_open())
    {
      m.invalidate_smiles();
      stream_for_failed_molecules << m.smiles() << ' ' << m.name() << reason_for_failure << '\n';
      stream_for_failed_molecules.write_if_buffer_holds_more_than(32768);
    }
  }
    
  return errors_this_molecule;
}

static int
tsmiles(Molecule & m)
{
  assert (m.ok());

  int matoms = m.natoms();

  if (verbose > 1)
    cerr << "Testing molecule " << molecules_read << " '" << m.name() << "' with " <<
            matoms << " atoms, " << m.nrings() << " rings\n";

  if (m.nrings() > max_rings)
  {
    cerr << "Too many rings " << m.nrings() << ", many potential problems, skipping '" << m.name() << "'\n";
    molecules_with_too_many_rings++;
    return 0;
  }

  if (remove_non_chiral_chiral_centres)
    do_remove_non_chiral_chiral_centres(m);

  if (unset_unnecessary_implicit_hydrogens_known_values)
    m.unset_unnecessary_implicit_hydrogens_known_values();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  const IWString usmiles = m.unique_smiles();

  molecular_weight_t amw;
  if (! m.contains_non_periodic_table_elements())
    amw = m.molecular_weight_ignore_isotopes();
  else
    amw = 0.0;

  if (verbose > 1)
    cerr << "Unique smiles " << usmiles << " amw = " << amw << '\n';

  if (verbose > 2)
    report_canonical_ranking(cerr, m);

  const resizable_array<atom_number_t> & atom_order_in_smiles = m.atom_order_in_smiles();
  if (matoms != atom_order_in_smiles.number_elements()) {
    cerr << "Size mismatch in atom order\n";
  }

  if (test_assignment)
    do_test_assignment(m, usmiles);

  if (test_subset_unique_smiles)
    do_test_subset_unique_smiles(m);

  int rc  = tsmiles(m, usmiles, amw);

  if (0 == rc)   // no failures
    return 0;

  if (! destroy_aromaticity)
    return rc;

  if (! m.contains_aromatic_atoms())
    return rc;

  if (! do_destroy_aromaticity(m))
    return rc;

  const IWString xusmi = m.unique_smiles();

  amw = m.molecular_weight_ignore_isotopes();

  rc = tsmiles(m, xusmi, amw);

  if (0 == rc)
  {
    molecules_ok_without_aromaticity++;
    if (verbose > 1)
      cerr << m.name() << " ok after aromaticity removal\n";
  }

  return rc;
}

/*
*/

static int
tsmiles(data_source_and_type<Molecule> & input)
{
  assert (input.good());

  int failures = 0;
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    if (verbose > 1)
      cerr << "Read " << input.molecules_read() << ", '" << m->name() << "'\n";

    molecules_read++;

    if (0 == m->natoms())
    {
      cerr << "Ignoring molecule with no atoms\n";
      continue;
    }

    failures += tsmiles(*m);

    if (failures && abort_on_error)
      break;

    if (report_progress())
      cerr << "Read " << molecules_read << " molecules, " << failures << " failures\n";
  }

  return failures;
}

static int
tsmiles(const char * fname, FileType input_type)
{
  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Processing '" << fname << "'\n";

  int failures = tsmiles(input);

  if (verbose)
    cerr << "Read " << input.molecules_read() << " molecules\n";

  return failures;
}

static int
tsmiles(int argc, char **argv)
{
  Command_Line cl(argc, argv, "jK:xp:s:A:vi:E:ac:t:bur:Rm:w:hqF:k:fWy");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(17);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl))
    usage(1);

  if (! process_standard_smiles_options(cl, verbose, 'K'))
    usage(2);

  FileType input_type = FILE_TYPE_INVALID;
  if (! cl.option_present('i'))
  {
    if (! all_files_recognised_by_suffix(cl))
    {
      cerr << "Cannot recognise file types by suffix(es)\n";
      return 3;
    }
  }
  else if (! process_input_type(cl, input_type))
    usage(2);

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose))
    {
      cerr << "cannot process aromaticit options\n";
      usage(3);
    }
  }
  else
  {
    cerr << "Warning, does not work well without Daylight aromaticity\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', permutations_per_molecule) || permutations_per_molecule <= 0)
    {
      cerr << "The -p option must be followed by a whole positive number\n";
      usage(4);
    }
  }
  else if (verbose)
    cerr << "Performing " << permutations_per_molecule << " by default\n";

  if (cl.option_present('x'))
  {
    abort_on_error = 1;
    if (verbose)
      cerr << "Programme will stop on encountering an error\n";
  }

  if (cl.option_present('y'))
  {
    unset_unnecessary_implicit_hydrogens_known_values = 1;

    if (verbose)
      cerr << "Will remove unnecessary square brackets\n";
  }

  if (cl.option_present('b'))
  {
    test_subset_unique_smiles = 1;
    if (verbose)
      cerr << "Will test unique smiles on subsets\n";
  }

  if (cl.option_present('j'))
  {
    stop_processing_molecule_on_first_error = 1;
    if (verbose)
      cerr << "No further permutations done a molecule for which an error has occurred\n";
  }

  if (cl.option_present('a'))
  {
    test_assignment = 1;
    if (verbose)
      cerr << "will test the assignment operator\n";
  }

  if (cl.option_present('f'))
  {
    destroy_aromaticity = 1;

    if (verbose)
      cerr << "Upon failure, will destroy aromaticity and see if molecule works without aromatic bonds\n";
  }

  if (cl.option_present('c'))
  {
    int include_chiral_info = 1;

    if (! cl.value('c', include_chiral_info))
    {
      cerr << "the -c option must be followed by a whole number\n";
      usage (7);
    }
    if (verbose)
      cerr << "Include chiral info set to " << include_chiral_info << '\n';

    set_include_chiral_info_in_smiles(include_chiral_info);
  }

  if (cl.option_present('W'))
  {
    remove_all_chiral_centres = 1;

    if (verbose)
      cerr << "Will remove all chiral centres before testing\n";
  }

  if (cl.option_present('h'))
  {
    test_include_chiral_centres_in_smiles = 1;
    if (verbose)
      cerr << "Will test set_include_chiral_info_in_smiles\n";
  }

  if (cl.option_present('q'))
  {
    test_same_ring_stuff = 1;

    if (verbose)
     cerr << "Will test in_same_ring and in_same_ring_system member functions\n";
  }

  if (cl.option_present('k'))
  {
    if (! cl.value('k', max_rings) || max_rings < 0)
    {
      cerr << "The maximum number of rings to process option (-k) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will not test molecules having more than " << max_rings << " rings\n";
  }

  if (cl.option_present('t'))
  {
    int t;
    if (! cl.value('t', t) || t < 0 || t > 1)
    {
      cerr << "The -t option must be either 0 or 1\n";
      usage(4);
    }

    set_include_cis_trans_in_smiles(t);
    set_include_directional_bonding_information_in_unique_smiles(t);

    if (verbose)
      cerr << "Include cis trans bonding flag set to " << t << '\n';

  }

  if (cl.option_present('u'))
  {
    remove_non_chiral_chiral_centres = 1;

    if (verbose)
      cerr << "False chiral centres will be removed\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress reporting\n";
      return 1;
    }
  }

  if (cl.option_present('R'))
  {
    write_sssr_upon_failure = 1;

    if (verbose)
      cerr << "Will write the SSSR upon failure\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', test_build_from_aromatic_form) || test_build_from_aromatic_form < 1)
    {
      cerr << "The number of steps for testing build from aromatic form (-m) must be +ve\n";
      usage(11);
    }

    set_input_aromatic_structures(1);

    if (verbose)
      cerr << "Will do " << test_build_from_aromatic_form << " tests of building from aromatic form\n";
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', test_swap_atoms) || test_swap_atoms < 1)
    {
      cerr << "The test swap atoms option (-w) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will swap " << test_swap_atoms << " atoms for each input molecule\n";
  }

  if (cl.option_present('s'))
  {
    random_number_seed_t s;
    if (! cl.value('s', s))
    {
      cerr << "Random number seeds (-s option) must be whole numbers\n";
      return 6;
    }

    set_smiles_random_number_seed(s);

    if (verbose)
      cerr << "Using random number seed " << s << '\n';
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments " << argc << "\n";
    usage(1);
  }

  if (cl.option_present('F'))
  {
    IWString f = cl.string_value('F');

    if (! f.ends_with(".smi"))
      f << ".smi";

    if (! stream_for_failed_molecules.open(f.null_terminated_chars()))
    {
      cerr << "Cannot open stream for failed molecules '" << f << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Failed molecules written to '" << f << "'\n";
  }

  int failures = 0;

  for (int i = 0; i < cl.number_elements(); i++)
    failures += tsmiles(cl[i], input_type);

  cerr << "After processing, " << failures << " failures\n";
  if (verbose || failures)
  {
    cerr << molecules_read << " molecules read " << permutations_performed <<
            " permutations performed\n";
    if (molecules_with_too_many_rings)
      cerr << "Skipped " << molecules_with_too_many_rings << " molecules_with_too_many_rings\n";
    if (molecules_with_errors)
    {
      cerr << total_error_count << " errors in " << molecules_with_errors << " molecules\n";
      cerr << smiles_interpretation_errors << " smiles interpretation errors\n";
      cerr << aromatic_smiles_interpretation_errors << " aromatic smiles interpretation errors\n";
      cerr << amw_errors << " molecular weight errors\n";
      cerr << unique_smiles_errors << " unique smiles errors\n";
      cerr << atom_swap_errors << " atom swap errors\n";

      if (destroy_aromaticity)
        cerr << molecules_ok_without_aromaticity << " molecules that failed were OK once aromaticity was removed\n";
    }
    else
      cerr << "All tests successful\n";
  }


  return (failures != 0);
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tsmiles(argc, argv);

  return rc;
}
