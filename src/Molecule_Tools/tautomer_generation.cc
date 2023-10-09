/*
  Enumerate tautomers.

  Each transformation gets a set of molecules to process.
  The transformation examines the first molecule in the
  set for a skeleton match. If that works, then all molecules
  in the set are examined more closely.

  The reason for that is that the first molecule in the set
  might not be in the right form for this transformation,
  but something that has been changed by some previous
  transformation might be OK
*/

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <limits>

#define RESIZABLE_ARRAY_IMPLEMENTATION 1

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#define ISTREAM_AND_TYPE_IMPLEMENTATION 1

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/toggle_kekule_form.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static extending_resizable_array<int> tautomer_stats;

static int molecules_written = 0;

static int tautomers_produced = 0;

static int process_acid = 0;
static int acids_flipped = 0;

static int process_keto_enol = 0;
static int keto_enol_flipped = 0;

static int allow_chain_amide_groups_to_flip = 0;

static int process_pyrrole = 0;
static int pyrrole_flipped = 0;

static int process_imidazole = 0;
static int imidazoles_flipped = 0;

static int process_tetrazole = 0;
static int tetrazoles_flipped = 0;

static int process_triazole_123 = 0;
static int triazole_123_flipped = 0;

static int process_triazole_124 = 0;
static int triazole_124_flipped = 0;

static int process_pyrazole = 0;
static int pyrazoles_flipped = 0;

static int five_membered_ring_kekule_toggled = 0;

static resizable_array_p<Substructure_Query> queries;
typedef resizable_array_p<Substructure_Query> Set_of_Queries;

static int reject_invalid_valence = 1;

static int invalid_valence_rejected = 0;

static int reject_for_aromaticity_loss = 0;

static int rejected_for_aromaticity_loss = 0;

static int must_preserve_hcount = 0;

static int rejected_for_hydrogen_change = 0;

static Molecule_Output_Object stream_for_invalid_valence;

static int max_tautomers_per_molecule = std::numeric_limits<int>::max();

static int terminated_by_max_tautomers = 0;

static int filter_to_unique_structures = 0;

static IWString prepend_before_sequence('.');

/*
  Useful for debugging, apply isotopic labels to all
  molecules as they are read
*/

static int apply_isotopic_labels = 0;

/*
  We need to keep track of molecules and the subset of atoms that can be
  processed
*/

class Molecule_Proc : public Molecule
{
  private:
    int * _process_these_atoms;

    int * _hcount;

  public:
    Molecule_Proc ();
    Molecule_Proc (const Molecule_Proc &);
    ~Molecule_Proc();

    int initialise (Set_of_Queries &);

    int process_atom(atom_number_t s) const { return _process_these_atoms[s];}

    int all_atoms_to_be_processed (const Set_of_Atoms & s) const;

    int determine_hcount ();

    int operator == (Molecule_Proc &);
};

Molecule_Proc::Molecule_Proc()
{
  _process_these_atoms = nullptr;
  _hcount = nullptr;

  return;
}

Molecule_Proc::Molecule_Proc (const Molecule_Proc & rhs) : Molecule (rhs)
{
  int matoms = rhs.natoms();

  assert (matoms > 0);

  _process_these_atoms = new int[matoms];

  copy_vector(_process_these_atoms, rhs._process_these_atoms, matoms);

  _hcount = nullptr;

  return;
}

Molecule_Proc::~Molecule_Proc()
{
  if (nullptr != _process_these_atoms)
    delete [] _process_these_atoms;

  if (nullptr != _hcount)
    delete [] _hcount;

  return;
}

int
Molecule_Proc::initialise (Set_of_Queries & queries)
{
  assert (nullptr == _process_these_atoms);

  int matoms = Molecule::natoms();

  _process_these_atoms = new int[matoms];

  int nq = queries.number_elements();

  if (0 == nq)
  {
    set_vector (_process_these_atoms, matoms, 1);
    return matoms;
  }

  set_vector(_process_these_atoms, matoms, 0);

  Molecule_to_Match target(this);

  for (int i = 0; i < queries.number_elements(); i++)
  {
    Substructure_Results sresults;

    if ( queries[i]->substructure_search(target, sresults))
      sresults.each_embedding_set_vector(_process_these_atoms, 1);
  }

  return count_non_zero_occurrences_in_array(_process_these_atoms, matoms);
}

int
Molecule_Proc::determine_hcount ()
{
  assert (nullptr == _hcount);

  int matoms = natoms();

  if (0 == matoms)
    return 0;

  _hcount = new int[matoms];

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    _hcount[i] = Molecule::implicit_hydrogens(i);

    rc += _hcount[i];
  }

  return rc;
}

int
Molecule_Proc::all_atoms_to_be_processed (const Set_of_Atoms & s) const
{
  return s.all_members_set_in_array(_process_these_atoms, 1);
}

int
Molecule_Proc::operator== (Molecule_Proc & rhs) 
{
  if (nullptr == rhs._hcount)
    rhs.determine_hcount();

  if (nullptr == _hcount)
    determine_hcount();

  int matoms = natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (_hcount[i] != rhs._hcount[i])
      return 0;
  }

  return 1;
}

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
  cerr << "  -T KE         enumerate Keto-Enol forms\n";
  cerr << "  -T ACID       enumerate acid tautomer forms\n";
  cerr << "  -T PY         enumerate pyrrole N tautomer forms\n";
  cerr << "  -s <smarts>   only process groups that match <smarts>\n";
  cerr << "  -q <query>    only process groups that match query <smarts>\n";
  cerr << "  -V            discard any tautomer generated with invalid valences\n";
  cerr << "  -a            discard any tautomer generated which loses aromaticity relative to parent\n";
  cerr << "  -h .          discard any tautomer generated that changes hcount\n";
  cerr << "  -X <n>        stop processing any molecules once <n> tautomers generated\n";
  cerr << "  -u 2          filter to unique structures - 2D consideration only\n";
  cerr << "  -u 3          filter to unique structures - 3D consideration\n";
  cerr << "  -Y <string>   prepend <string> before tautomer number\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
do_apply_isotopic_labels (Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    m.set_isotope(i, i);
  }

  return;
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (apply_isotopic_labels)
    do_apply_isotopic_labels(m);

  return;
}

/*
  Is there another aromatic ring that is fused to R.
  To make things easy, we only consider singly fused rings
  We return the ring number
*/

static const Ring *
identify_fused_aromatic_neighbour(Molecule_Proc & m,
                                  const Ring & r)
{
  m.compute_aromaticity_if_needed();

  for (int i = 0; i < r.fused_ring_neighbours(); i++)
  {
    const Ring * n = r.fused_neighbour(i);

    if (! n->is_aromatic())
      continue;

//  hard to know whether to include this or not
//  if (! m0.all_atoms_to_be_processed(*n))
//    continue;

    if (1 == n->fused_ring_neighbours())   // not fused to anything else, good
      return n;
  }

  return nullptr;
}

#ifdef NOT_USED_HERE
static int
contains_pyridine_and_pyrrole_nitrogens (Molecule_Proc & m,
                                         const Ring & r)
{
  int rs = r.number_elements();

  int pyrrole_forms = 0;
  int pyridine_forms = 0;

  for (int i = 0; i < rs; i++)
  {
    atom_number_t j = r[i];

    if (! m.process_atom(j))
      continue;

    const Atom * aj = m.atomi(j);

    if (7 != aj->atomic_number())
      continue;

    if (2 != aj->ncon())
      continue;

    if (0 == m.hcount(j))
      pyridine_forms++;
    else
      pyrrole_forms++;
  }

  if (0 == pyridine_forms)
    return 0;

  if (0 == pyrrole_forms)
    return 0;

  if (pyrrole_forms > 1)   // too difficult
    return 0;

  return 1;
}
#endif

/*
  The Ring class is handy for keeping track of the atoms around the
  ring system, so it will actually comprise two rings
*/

class Fused_Ring
{
  private:
    Ring _around_edge;

    resizable_array<int> _nitrogen_ndx;

//  private functions

    int _identify_atoms_around_edge (Molecule_Proc & m,
                                         atom_number_t zatom,
                                         int * atom_in_system);
  public:
    Fused_Ring();

    int is_fused_ring_with_nitrogens(Molecule_Proc & m,
                                     int * atom_in_system);
};

int
Fused_Ring::_identify_atoms_around_edge (Molecule_Proc & m,
                                         atom_number_t zatom,
                                         int * atom_in_system)
{
  const Atom * a = m.atomi(zatom);

  if (7 == a->atomic_number())
    _nitrogen_ndx.add(_around_edge.number_elements());  // before we add it

  _around_edge.add(zatom);

  atom_in_system[zatom] = 0;

  int acon = a->ncon();

  for (int j = 0; j < acon; j++)
  {
    const Bond * b = a->item(j);

    if (1 != b->nrings())
      continue;

    atom_number_t k = b->other(zatom);
    
    if (atom_in_system[k])
      _identify_atoms_around_edge(m, k, atom_in_system);
  }

  return 1;
}

int
Fused_Ring::is_fused_ring_with_nitrogens (Molecule_Proc & m,
                                          int * atom_in_system)
{
  int matoms = m.natoms();

  m.ring_membership();

  atom_number_t f = locate_item_in_array(1, matoms, atom_in_system);
  assert (f >= 0);

  if (! _identify_atoms_around_edge (m, f, atom_in_system))
    return 0;

  if (_nitrogen_ndx.number_elements() < 2)
    return 0;

// Now we need to figure out how far apart are the nitrogens

  return 1;
}

/*
  Two ring systems only.
*/

#ifdef NOT_USED_HERE
static int
walk_pyrrole_around_ring_system (resizable_array_p<Molecule_Proc> & m,
                     Molecule_Output_Object & output)
{
  Molecule_Proc * m0 = m[0];

  int nr = m0->nrings();

  if (0 == nr)
    return 1;

  m0->compute_aromaticity_if_needed();

  int * atom_in_system = new int[m0->natoms()]; std::unique_ptr<int[]> free_atom_in_system(atom_in_system);

  int * ring_in_system = new_int(nr); std::unique_ptr<int[]> free_ring_in_system(ring_in_system);

  resizable_array_p<Ring> sites_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (ri->fused_ring_neighbours() > 1)   // too hard
      continue;

    if (5 != ri->number_elements())
      continue;

    if (! contains_pyridine_and_pyrrole_nitrogens(*m0, *ri))
      continue;

    if (m0->fused_system_size(ri->item(0)) > 2)   // too hard, leave for later...
      continue;

    const Ring * fsdnbr = identify_fused_aromatic_neighbour(*m0, *ri);

    ri->set_vector(atom_in_system, 1);
    if (nullptr != fsdnbr)
      fsdnbr->set_vector(atom_in_system, 1);
  }

  return 1;
}
#endif

class Pyrazole
{
  private:

    const atom_number_t _n1;
    const atom_number_t _n2;
    atom_number_t _c3;
    atom_number_t _c4;
    atom_number_t _c5;

  public:
    Pyrazole (atom_number_t, atom_number_t);

//  After we have the nitrogens, we need to determine the carbons

    int determine_carbon_atoms (const Ring &, int, int);

    int process (Molecule_Proc &, resizable_array_p<Molecule_Proc> &) const;
};

Pyrazole::Pyrazole(atom_number_t a1, atom_number_t a2) : _n1(a1), _n2(a2)
{
  _c3 = INVALID_ATOM_NUMBER;
  _c4 = INVALID_ATOM_NUMBER;
  _c5 = INVALID_ATOM_NUMBER;

  return;
}

int
Pyrazole::determine_carbon_atoms (const Ring & r,
                                    int n1_ndx,
                                    int n2_ndx)
{
//cerr << "Processing ring " << r << endl;

  assert (_n1 == r[n1_ndx]);
  assert (_n2 == r[n2_ndx]);

  int direction = 1;
  if (n1_ndx + 1 == n2_ndx)   // great
    ;
  else if (n2_ndx == r.next_index_after_wrap(n1_ndx, -1))
    direction = -1;
  else     // should not happen
    return 0;

// note that n2_ndx is changed in each of these calls, but
// it is no longer needed after determining _c3

  _c3 = r.next_after_wrap(n2_ndx, direction);
  _c4 = r.next_after_wrap(n2_ndx, direction);
  _c5 = r.next_after_wrap(n2_ndx, direction);

//cerr << "Final " << _n1 << ", " << _n2 << ", " << _c3 << ", " << _c4 << ", " << _c5 << endl;

  return 1;
}

int
Pyrazole::process (Molecule_Proc & m,
             resizable_array_p<Molecule_Proc>& molecules) const
{
  const Bond * b23 = m.bond_between_atoms(_n2, _c3);

  Molecule_Proc * tmp = new Molecule_Proc(m);

  if (b23->is_single_bond())
  {
    tmp->set_bond_type_between_atoms(_n2, _c3, DOUBLE_BOND);
    tmp->set_bond_type_between_atoms(_c3, _c4, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_c4, _c5, DOUBLE_BOND);
    tmp->set_bond_type_between_atoms(_c5, _n1, SINGLE_BOND);
  }
  else
  {
    tmp->set_bond_type_between_atoms(_n2, _c3, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_c3, _c4, DOUBLE_BOND);
    tmp->set_bond_type_between_atoms(_c4, _c5, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_c5, _n1, DOUBLE_BOND);
  }

  molecules.add(tmp);

  pyrazoles_flipped++;

  return 1;
}

static int
do_pyrazoles (const Pyrazole & p,
              resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * m = molecules[i];

    p.process(*m, molecules);

    if (molecules.number_elements() >= max_tautomers_per_molecule)
    {
      terminated_by_max_tautomers++;
      return 1;
    }
  }

  return 1;
}

static int
looks_like_pyrazole (Molecule_Proc & m,
                     const Ring & r,
                     int & n1_ndx,
                     int & n2_ndx)
{
  n1_ndx = -1;
  n2_ndx = -1;

  int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    atomic_number_t z = m.atomic_number(j);

    if (6 == z)
      continue;

    if (7 != z)
      return 0;

    if (2 != m.ncon(j))   // impossible
      return 0;

    if (n1_ndx < 0)
      n1_ndx = i;
    else if (n2_ndx < 0)
      n2_ndx = i;
    else
      return 0;
  }

  if (n2_ndx < 0)
    return 0;

  if (n1_ndx + 1 == n2_ndx)
    ;
  else if (n2_ndx == r.next_index_after_wrap(n1_ndx, -1))
    ;
  else
    return 0;

  return 1;
}

static int
do_pyrazoles (resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  if (0 == n)
    return 0;

  Molecule_Proc * m0 = molecules[0];

  m0->compute_aromaticity_if_needed();

  int nr = m0->nrings();

  if (0 == nr)
    return 0;

  resizable_array_p<Pyrazole> rings_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (! ri->is_aromatic())
      continue;

    if (! m0->all_atoms_to_be_processed(*ri))
      continue;

    if (nullptr != identify_fused_aromatic_neighbour(*m0, *ri))
      continue;

    int n1_ndx, n2_ndx;
    if (! looks_like_pyrazole(*m0, *ri, n1_ndx, n2_ndx))
      continue;

    Pyrazole * tmp = new Pyrazole(ri->item(n1_ndx), ri->item(n2_ndx));

    tmp->determine_carbon_atoms(*ri, n1_ndx, n2_ndx);

    rings_to_process.add(tmp);
  }

  if (rings_to_process.empty())
    return 0;

  if (verbose > 1)
    cerr << "In '" << m0->name() << "' found " << rings_to_process.number_elements() << " possible pyrazoles\n";

  for (int i = 0; i < rings_to_process.number_elements(); i++)
  {
    const Pyrazole * p = rings_to_process[i];

    do_pyrazoles(*p, molecules);

    if (molecules.number_elements() >= max_tautomers_per_molecule)
      break;
  }

  return 1;
}

/*
  Place this kind of bonding pattern. The top atom is the pyrrole

     N               P
   /   \           /   \
  N     N         5     2
  ||    ||        |     |
  ||    ||        |     |
  N --- N         4 --- 3

*/

static int
add_tetrazole_pyrrole_on (Molecule_Proc & m,
                          atom_number_t pyrrole,
                          atom_number_t a2,
                          atom_number_t a3,
                          atom_number_t a4,
                          atom_number_t a5,
                          resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * tmp = new Molecule_Proc(m);

  tmp->set_bond_type_between_atoms(pyrrole, a2, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(a2, a3, DOUBLE_BOND);
  tmp->set_bond_type_between_atoms(a3, a4, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(a4, a5, DOUBLE_BOND);
  tmp->set_bond_type_between_atoms(a5, pyrrole, SINGLE_BOND);

  molecules.add(tmp);

  return 1;
}


/*
  There are four possible tautomers of a tetrazole

     C           C          C          C
   //  \       /  \\      //  \      /   \\
  N     N     N     N    N     N    N     N
  |     ||    |     |    |     |   ||     |
  |     ||    |     |    |     |   ||     |
  N --- N     N === N    N === N    N --- N

  case1        case2      case3      case4
*/

static int
do_tetrazole (Molecule_Proc & m,
              const Ring & r,
              resizable_array_p<Molecule_Proc> & molecules)
{
  int c_ndx = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (6 != m.atomic_number(j))
      continue;

    c_ndx = i;
    break;
  }

  if (c_ndx < 0)    // should not happen
    return 0;

  int n2_ndx = r.next_index_after_wrap(c_ndx);
  int n3_ndx = r.next_index_after_wrap(n2_ndx);
  int n4_ndx = r.next_index_after_wrap(n3_ndx);
  int n5_ndx = r.next_index_after_wrap(n4_ndx);

  atom_number_t c = r[c_ndx];
  atom_number_t n2 = r[n2_ndx];
  atom_number_t n3 = r[n3_ndx];
  atom_number_t n4 = r[n4_ndx];
  atom_number_t n5 = r[n5_ndx];

// we only need two bonds to get all cases

  const Bond * c1n2 = m.bond_between_atoms(c, n2);
  const Bond * n3n4 = m.bond_between_atoms(n3, n4);

/*
     C           C          C          C           1
   //  \       /  \\      //  \      /   \\      /   \
  N     N     N     N    N     N    N     N     5     2
  |     ||    |     |    |     |   ||     |     |     |
  |     ||    |     |    |     |   ||     |     |     |
  N --- N     N === N    N === N    N --- N     4 --- 3

  case1        case2      case3      case4
*/

  if (c1n2->is_single_bond() && n3n4->is_single_bond())   // case1
  {
    add_tetrazole_pyrrole_on(m, n2, n3, n4, n5, c, molecules);
    add_tetrazole_pyrrole_on(m, n3, n4, n5, c, n2, molecules);
    add_tetrazole_pyrrole_on(m, n5, c, n2, n3, n4, molecules);
  }
  else if (c1n2->is_double_bond() && n3n4->is_double_bond())  // case2
  {
    add_tetrazole_pyrrole_on(m, n2, n3, n4, n5, c, molecules);
    add_tetrazole_pyrrole_on(m, n3, n4, n5, c, n2, molecules);
    add_tetrazole_pyrrole_on(m, n4, n5, c, n2, n3, molecules);
  }
  else if (c1n2->is_single_bond() && n3n4->is_double_bond())   // case 3
  {
    add_tetrazole_pyrrole_on(m, n3, n4, n5, c, n2, molecules);
    add_tetrazole_pyrrole_on(m, n4, n5, c, n2, n3, molecules);
    add_tetrazole_pyrrole_on(m, n5, c, n2, n3, n4, molecules);
  }
  else if (c1n2->is_double_bond() && n3n4->is_single_bond())   // case4
  {
    add_tetrazole_pyrrole_on(m, n2, n3, n4, n5, c, molecules);
    add_tetrazole_pyrrole_on(m, n4, n5, c, n2, n3, molecules);
    add_tetrazole_pyrrole_on(m, n5, c, n2, n3, n4, molecules);
  }
  else             // should not happen
    return 0;

  tetrazoles_flipped += 3;

  return 1;
}

static int
do_tetrazoles(const Ring & r,
              resizable_array_p<Molecule_Proc> &molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * m = molecules[i];

    do_tetrazole (*m, r, molecules);
  }

  return 1;
}

static int
is_tetrazole (const Molecule & m,
              const Ring & r)
{
  int ncount = 0;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    atomic_number_t z = m.atomic_number(j);

    if (7 == z)
      ncount++;
    else if (6 == z)
      ;
    else
      return 0;
  }

  return 4 == ncount;
}

static int
do_tetrazoles (resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * m0 = molecules[0];

  int nr = m0->nrings();

  if (0 == nr)
    return 1;

  resizable_array<int> rings_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (! ri->is_aromatic())
      continue;

    if (ri->is_fused())
      continue;

    if (! is_tetrazole(*m0, *ri))
      continue;

    rings_to_process.add(i);
  }

  if (rings_to_process.empty())
    return 1;

  if (verbose > 1)
    cerr << "Found " << rings_to_process.number_elements() << " potential 124 triazoles in '" << m0->name() << "'\n";

  for (int i = 0; i < rings_to_process.number_elements(); i++)
  {
    int j = rings_to_process[i];

    const Ring * r = m0->ringi(j);

    Ring rcopy(*r);   // need a stable copy because m0 might change

    do_tetrazoles(rcopy, molecules);
  }

  return 1;
}

#ifdef NOT_USED_HERE
static int
is_imidazole (Molecule_Proc & m,
              const Ring & r,
              atom_number_t pyrrole,
              int center_index,
              atom_number_t & centre,
              int pyridine_index,
              atom_number_t & pyridine)
{
//#define DEBUG_IMIDAZOLE
#ifdef DEBUG_IMIDAZOLE
  cerr << "Looking for imidazole in " << r << endl;
  cerr << "Pyrrole " << pyrrole << " '" << m.smarts_equivalent_for_atom(pyrrole) << "'\n";
  cerr << "centre " << r[center_index] << " '" << m.smarts_equivalent_for_atom(r[center_index]) << "'\n";
  cerr << "Pyridine " << r[pyridine_index] << " '" << m.smarts_equivalent_for_atom(r[pyridine_index]) << "'\n";
#endif

  m.ring_membership();

  centre = r[center_index];

  if (! m.process_atom(centre))
    return 0;

  const Atom * c = m.atomi(centre);

  int ccon = c->ncon();

// Look for a double bond to the pyridine nitrogen

  for (int i = 0; i < ccon; i++)
  {
    const Bond * b = c->item(i);

    atom_number_t j = b->other(centre);

    if (! m.process_atom(j))
      continue;

    if (j == pyrrole)
      continue;

    const Atom * aj = m.atomi(j);

    if (8 == aj->atomic_number())    // will get picked up as keto-enol
      return 0;

    if (0 == b->nrings())    // we are now looking at the other bond in the ring
      continue;

    if (7 != m.atomic_number(j))
      return 0;

    if (! b->is_double_bond())
      return 0;

    if (0 != aj->formal_charge())
      return 0;

    if (3 != aj->nbonds())
      return 0;

    pyridine = j;
    return j == r[pyridine_index];
  }

  return 0;
}

static int
is_imidazole_form (Molecule_Proc & m,
                   const Ring & r,
                   atom_number_t & pyridine,
                   atom_number_t & something,
                   int pyrrole_index)
{
  int ring_size = r.number_elements();

  if (is_imidazole (m, r, r[pyrrole_index], (pyrrole_index + 1) % ring_size, something, (pyrrole_index + 2) % ring_size, pyridine))
    return 1;

  if (is_imidazole (m, r, r[pyrrole_index], (pyrrole_index - 1 + ring_size) % ring_size, something, (pyrrole_index - 2 + ring_size) % ring_size, pyridine))
    return 1;

  return 0;
}
#endif

class Triazole_124
{
  private:
    atom_number_t _n1;
    atom_number_t _n2;
    atom_number_t _c3;
    atom_number_t _n4;
    atom_number_t _c5;

  public:
    Triazole_124 ();

    int is_124_triazole (Molecule_Proc & m, const Ring &);

    int process (Molecule_Proc & m, resizable_array_p<Molecule_Proc> &) const;
};

Triazole_124::Triazole_124 ()
{
  _n1 = INVALID_ATOM_NUMBER;
  _n2 = INVALID_ATOM_NUMBER;
  _c3 = INVALID_ATOM_NUMBER;
  _n4 = INVALID_ATOM_NUMBER;
  _c5 = INVALID_ATOM_NUMBER;
}

int
Triazole_124::is_124_triazole (Molecule_Proc & m,
                               const Ring & r)
{
  int n = r.number_elements();

  int n1_ndx = -1;
  int n2_ndx = -1;
  int n4_ndx = -1;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    atomic_number_t z = m.atomic_number(j);

    if (6 == z)
      continue;

    if (7 != z)
      continue;

    if (2 != m.ncon(j))
      return 0;

    if (n1_ndx < 0)
      n1_ndx = i;
    else if (n2_ndx < 0)
      n2_ndx = i;
    else if (n4_ndx < 0)
      n4_ndx = i;
    else
      return 0;
  }

  if (n4_ndx < 0)
    return 0;

// Is this a 124 or a 123 triazole

//#define DEBUG_124_TRIAZOLE
#ifdef DEBUG_124_TRIAZOLE
  cerr << "Maybe 124 triazole indices: n1 " << _n1_ndx << " n2 " << _n2_ndx << " c3 " << _c3_ndx << " n4 " << _n4_ndx << " c5 " << _c5_ndx << endl;
#endif

  if (n1_ndx + 1 == n2_ndx && n4_ndx == r.next_index_after_wrap(n2_ndx))
    return 0;

// We have a 124 triazole. Now make sure th indices are correct
//    N1 N2 c  N4 c               no change needed
//    c  N1 N2 c  N4              no change needed
//    N1 c  N2 N4 c
//    c  N1 c  N2 N4
//  In all cases make sure N1-N2 are bonded and preserve directionality

  if (n1_ndx + 1 == n2_ndx && n2_ndx + 2 == n4_ndx)   // cases 1 and 2
    ;
  else if (n1_ndx + 2 == n2_ndx && n2_ndx + 1 == n4_ndx)  // cases 3 and 4
  {
    int tmp = n4_ndx;
    n4_ndx = n1_ndx;
    n1_ndx = n2_ndx;
    n2_ndx = tmp;
  }
  else
    return 0;

  _n1 = r[n1_ndx];
  _n2 = r[n2_ndx];
  _n4 = r[n4_ndx];

// n2_ndx and n4_ndx are no longer needed, so OK to destroy them here

  _c3 = r.next_after_wrap(n2_ndx);
  _c5 = r.next_after_wrap(n4_ndx);

  return 1;
}

/*
  We need to form three tautomers


      N              N             N
   //   \          /   \\        /   \
  C      C        C      C     C      C
  |      ||       ||     |     ||     ||
  N ---- N        N ---- N     N ---- N

    case1           case2        case3

  first task is to identify the NN bond
*/

int
Triazole_124::process (Molecule_Proc & m,
                       resizable_array_p<Molecule_Proc> & molecules) const
{
//cerr << "Examining 124 triazole " << r << endl;

// The N1-N2 bond must be a single bond

  const Bond * b = m.bond_between_atoms(_n1, _n2);

  if (! b->is_single_bond())
    return 0;

#ifdef DEBUG_124_TRIAZOLE
  cerr << "124 triazole indices: n1 " << _n1_ndx << " n2 " << _n2_ndx << " c3 " << _c3_ndx << " n4 " << _n4_ndx << " c5 " << _c5_ndx << endl;
#endif

  const Bond * n2c3 = m.bond_between_atoms(_n2, _c3);
  const Bond * c5n1 = m.bond_between_atoms(_c5, _n1);

/*
      N              N             N             4
   //   \          /   \\        /   \         /    \
  C      C        C      C     C      C       5      3
  |      ||       ||     |     ||     ||      |      |
  N ---- N        N ---- N     N ---- N       1 ---- 2

    case1           case2        case3
*/

  if (c5n1->is_single_bond() && n2c3->is_double_bond())   // case 1 above
  {
    Molecule_Proc * t2 = new Molecule_Proc(m);       // generate case 2
    t2->set_bond_type_between_atoms(_n2, _c3, SINGLE_BOND);
    t2->set_bond_type_between_atoms(_c3, _n4, DOUBLE_BOND);
    t2->set_bond_type_between_atoms(_n4, _c5, SINGLE_BOND);
    t2->set_bond_type_between_atoms(_c5, _n1, DOUBLE_BOND);
    molecules.add(t2);
    Molecule_Proc * t3 = new Molecule_Proc(m);       // generate case 3
    t3->set_bond_type_between_atoms(_n4, _c5, SINGLE_BOND);
    t3->set_bond_type_between_atoms(_c5, _n1, DOUBLE_BOND);
    molecules.add(t3);
  }
  else if (c5n1->is_double_bond() && n2c3->is_single_bond())   // case 2 above
  {
    Molecule_Proc * t1 = new Molecule_Proc(m);       // generate case 1
    t1->set_bond_type_between_atoms(_n2, _c3, DOUBLE_BOND);
    t1->set_bond_type_between_atoms(_c3, _n4, SINGLE_BOND);
    t1->set_bond_type_between_atoms(_n4, _c5, DOUBLE_BOND);
    t1->set_bond_type_between_atoms(_c5, _n1, SINGLE_BOND);
    molecules.add(t1);
    Molecule_Proc * t3 = new Molecule_Proc(m);       // generate case 3
    t3->set_bond_type_between_atoms(_n2, _c3, DOUBLE_BOND);
    t3->set_bond_type_between_atoms(_c3, _n4, SINGLE_BOND);
    molecules.add(t3);
  }
  else if (c5n1->is_double_bond() && n2c3->is_double_bond())   // case 3 above
  {
    Molecule_Proc * t1 = new Molecule_Proc(m);       // generate case 1
    t1->set_bond_type_between_atoms(_n4, _c5, DOUBLE_BOND);
    t1->set_bond_type_between_atoms(_c5, _n1, SINGLE_BOND);
    molecules.add(t1);
    Molecule_Proc * t2 = new Molecule_Proc(m);       // generate case 2
    t2->set_bond_type_between_atoms(_n2, _c3, SINGLE_BOND);
    t2->set_bond_type_between_atoms(_c3, _n4, DOUBLE_BOND);
    molecules.add(t2);
  }

  triazole_124_flipped++;

  return 1;
}

static int
do_isolated_124_triazoles(const Triazole_124 & t,
                          resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * mi = molecules[i];

    t.process(*mi, molecules);
  }

  return 1;
}

static int
do_isolated_124_triazoles (resizable_array_p<Molecule_Proc> & molecules)
{

  Molecule_Proc * m0 = molecules[0];

  int nr = m0->nrings();

  if (0 == nr)
    return 1;

  m0->compute_aromaticity_if_needed();

  resizable_array_p<Triazole_124> rings_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (! ri->is_aromatic())
      continue;

    if (ri->is_fused())
      continue;

    Triazole_124 * t = new Triazole_124();

    if (! t->is_124_triazole(*m0, *ri))
      delete t;
    else
      rings_to_process.add(t);
  }

  if (rings_to_process.empty())
    return 1;

  if (verbose > 1)
    cerr << "Found " << rings_to_process.number_elements() << " potential 124 triazoles in '" << m0->name() << "'\n";

  for (int i = 0; i < rings_to_process.number_elements(); i++)
  {
    const Triazole_124 * t = rings_to_process[i];

    do_isolated_124_triazoles (*t, molecules);
  }

  return 1;
}

/*
  
*/

#ifdef NOT_USED_HERE
static int
is_aromatic_fused (const Molecule_Proc & m,
                   const Ring & r)
{
  if (0 == r.fused_ring_neighbours())
    return 0;

  if (r.fused_ring_neighbours() > 1)
    return 1;

  const Ring * nbr = r.fused_neighbour(0);

  if (nbr->is_aromatic())
    return 0;

  return 1;
}
#endif

class Triazole_123
{
  private:

    atom_number_t _n1;
    atom_number_t _n2;
    atom_number_t _n3;
    atom_number_t _c4;
    atom_number_t _c5;

//  private functions

    int _do_123_triazole_from_case3 (Molecule_Proc & m,
                            resizable_array_p<Molecule_Proc> & molecules) const;

    int _do_123_triazole_to_case3 (Molecule_Proc & m,
                          resizable_array_p<Molecule_Proc> & molecules) const;

    int _do_123_triazole_pyrrole_end (const Molecule_Proc & m,
                             atom_number_t n1,
                             atom_number_t n2,
                             atom_number_t n3,
                             resizable_array_p<Molecule_Proc> & molecules) const;

  public:
    Triazole_123 ();

    int is_123_triazole (Molecule_Proc &, const Ring &);

    int process (Molecule_Proc &, resizable_array_p<Molecule_Proc> &) const;
};

Triazole_123::Triazole_123 ()
{
  _n1 = INVALID_ATOM_NUMBER;
  _n2 = INVALID_ATOM_NUMBER;
  _n3 = INVALID_ATOM_NUMBER;
  _c4 = INVALID_ATOM_NUMBER;
  _c5 = INVALID_ATOM_NUMBER;

  return;
}

/*
  If any of the N atoms in a 123 triazole are substituted, it cannot
  change
*/

int
Triazole_123::is_123_triazole (Molecule_Proc & m,
                               const Ring & r)
{
  int something_with_hcount = 0;

  int n = r.number_elements();

  int n1_ndx = -1;
  int n2_ndx = -1;
  int n3_ndx = -1;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    const Atom * aj = m.atomi(j);

    atomic_number_t z = aj->atomic_number();

    if (6 == z)
      continue;

    if (7 != z)
      return 0;

    if (n1_ndx < 0)
      n1_ndx = i;
    else if (n2_ndx < 0)
      n2_ndx = i;
    else if (n3_ndx < 0)
      n3_ndx = i;
    else            // must be a tetrazole
      return 0;

    if (2 != aj->ncon())
      return 0;

    if (something_with_hcount)
      ;
    else if (m.implicit_hydrogens(j))
      something_with_hcount++;
  }

//#define DEBUG_123_TRIAZOLE
#ifdef DEBUG_123_TRIAZOLE
  cerr << "123 triazole indices " << n1_ndx << ", " << n2_ndx << ", " << n3_ndx << endl;
#endif

  if (n3_ndx < 0)
    return 0;

  if (1 != something_with_hcount)
    return 0;

// the three nitrogen atoms must be consecutive.
// NNNCC   1
// CNNNC   2
// CCNNN   3
// NCCNN   4
// NNCCN   5

  if (n2_ndx != n1_ndx + 1)
    return 0;

  if (n1_ndx + 1 == n2_ndx && n2_ndx + 1 == n3_ndx)   // cases 1-3
    ;
  else if (0 == n1_ndx && 3 == n2_ndx && 4 == n3_ndx)  // case 4
  {
    n1_ndx = 3;
    n2_ndx = 4;
    n3_ndx = 0;
  }
  else if (0 == n1_ndx && 1 == n2_ndx && 4 == n3_ndx)  // case 5
  {
    n1_ndx = 4;
    n2_ndx = 0;
    n3_ndx = 1;
  }
  else
    return 0;

   _n1 = r[n1_ndx];
   _n2 = r[n2_ndx];
   _n3 = r[n3_ndx];

// n3_ndx finished with, so we can destroy it here

  _c4 = r.next_after_wrap(n3_ndx);
  _c5 = r.next_after_wrap(n3_ndx);

#ifdef DEBUG_123_TRIAZOLE
  cerr << "123 triazole returning true\n";
#endif

  return 1;
}

static int
in_more_than_one_aromatic_ring (Molecule_Proc & m,
                                atom_number_t zatom)
{
  int zatom_nrings = m.nrings(zatom);

  if (zatom_nrings < 2)
    return 0;

  int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (ri->contains(zatom))
      rc++;
  }

  if (rc > 1)
    return 1;

  return 0;
}

/*
  We are being called upon to convert case 3 into case1 and case2

     N          N          N             2
   //  \      /   \\     /   \         /   \
  N     N    N     N    N     N       1     3
  |     |    |     |    ||    ||      |     |
  C === C    C === C    C --- C       5 --- 4
  
   case1      case2      case3
*/

int
Triazole_123::_do_123_triazole_from_case3 (Molecule_Proc & m,
                            resizable_array_p<Molecule_Proc> & molecules) const
{
  if (in_more_than_one_aromatic_ring(m, _c4))
    return 0;

  Molecule_Proc * tmp1 = new Molecule_Proc(m);
  Molecule_Proc * tmp2 = new Molecule_Proc(m);

  tmp1->set_bond_type_between_atoms(_c4, _c5, DOUBLE_BOND);
  tmp2->set_bond_type_between_atoms(_c4, _c5, DOUBLE_BOND);

  tmp1->set_bond_type_between_atoms(_n3, _c4, SINGLE_BOND);
  tmp2->set_bond_type_between_atoms(_n3, _c4, SINGLE_BOND);

  tmp1->set_bond_type_between_atoms(_c5, _n1, SINGLE_BOND);
  tmp2->set_bond_type_between_atoms(_c5, _n1, SINGLE_BOND);

  tmp1->set_bond_type_between_atoms(_n1, _n2, DOUBLE_BOND);   // case1
  tmp2->set_bond_type_between_atoms(_n1, _n2, SINGLE_BOND);   // case2

  tmp1->set_bond_type_between_atoms(_n2, _n3, SINGLE_BOND);   // case1
  tmp2->set_bond_type_between_atoms(_n2, _n3, DOUBLE_BOND);   // case2

  molecules.add(tmp1);
  molecules.add(tmp2);

  return 1;
}


/*
  We are being called upon to convert case1 or case2 into case3

     N          N          N             2
   //  \      /   \\     /   \         /   \
  N     N    N     N    N     N       1     3
  |     |    |     |    ||    ||      |     |
  C === C    C === C    C --- C       5 --- 4
  
   case1      case2      case3
*/

int
Triazole_123::_do_123_triazole_to_case3 (Molecule_Proc & m,
                          resizable_array_p<Molecule_Proc> & molecules) const
{
  if (in_more_than_one_aromatic_ring(m, _c4))
    return 0;

  Molecule_Proc * tmp = new Molecule_Proc(m);

  tmp->set_bond_type_between_atoms(_n1, _n2, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(_n2, _n3, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(_n3, _c4, DOUBLE_BOND);
  tmp->set_bond_type_between_atoms(_c4, _c5, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(_c5, _n1, DOUBLE_BOND);

  molecules.add(tmp);

  return 1;
}

/*
  Just flip the bonds between n1 n2 and n3.
  N1 is the atom that is now the pyrrole
*/

int
Triazole_123::_do_123_triazole_pyrrole_end (const Molecule_Proc & m,
                             atom_number_t n1,
                             atom_number_t n2,
                             atom_number_t n3,
                             resizable_array_p<Molecule_Proc> & molecules) const
{
  Molecule_Proc * tmp = new Molecule_Proc(m);

  tmp->set_bond_type_between_atoms(n1, n2, SINGLE_BOND);
  tmp->set_bond_type_between_atoms(n2, n3, DOUBLE_BOND);

  molecules.add(tmp);

  return 1;
}

/*
  There are three tautomers for 123 triazoles

     N          N          N             1
   //  \      /   \\     /   \         /   \
  N     N    N     N    N     N       5     2
  |     |    |     |    ||    ||      |     |
  C === C    C === C    C --- C       4 --- 3
  
   case1      case2      case3
*/

int
Triazole_123::process (Molecule_Proc & m,
                 resizable_array_p<Molecule_Proc> & molecules) const
{
  const Bond * n1n2 = m.bond_between_atoms(_n1, _n2);
  const Bond * n2n3 = m.bond_between_atoms(_n2, _n3);

/*
     N          N          N             2
   //  \      /   \\     /   \         /   \
  N     N    N     N    N     N       1     3
  |     |    |     |    ||    ||      |     |
  C === C    C === C    C --- C       4 --- 3
  
   case1      case2      case3
*/

#ifdef DEBUG_123_TRIAZOLE
  cerr << "Processing 123 triazole indices: n1 " << n1_ndx << " n2 " << n2_ndx << " n3 " << n3_ndx << endl;
  cerr << "Processing 123 triazole atoms:   n1 " << n1 << " n2 " << n2 << " n3 " << n3 << endl;
#endif

  if (n1n2->is_double_bond() && n2n3->is_single_bond() && m.implicit_hydrogens(_n3))     // case1
  {
    _do_123_triazole_pyrrole_end(m, _n1, _n2, _n3, molecules);
    _do_123_triazole_to_case3 (m, molecules);
  }
  else if (n1n2->is_single_bond() && n2n3->is_double_bond() && m.implicit_hydrogens(_n1))    // case2
  {
    _do_123_triazole_pyrrole_end(m, _n3, _n2, _n1, molecules);
    _do_123_triazole_to_case3 (m, molecules);
  }
  else if (n1n2->is_single_bond() && n2n3->is_single_bond() && m.implicit_hydrogens(_n2))    // case3
    _do_123_triazole_from_case3(m, molecules);
  else           // should not happen
  {
    cerr << "Default case not covered\n";
    return 0;
  }

  triazole_123_flipped++;

  return 1;
}

static int
do_123_triazoles (const Triazole_123 & t,
                  resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * mi = molecules[i];

    t.process(*mi, molecules);
  }

  return 1;
}

static int
do_123_triazoles (resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * m0 = molecules[0];

  int nr = m0->nrings();

  if (0 == nr)
    return 1;

  m0->compute_aromaticity_if_needed();

  resizable_array_p<Triazole_123> rings_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (! ri->is_aromatic())
      continue;

    Triazole_123 * tmp = new Triazole_123();
    if (! tmp->is_123_triazole(*m0, *ri))
      delete tmp;
    else
      rings_to_process.add(tmp);
  }

  if (rings_to_process.empty())
    return 0;

  if (verbose > 1)
    cerr << "Found " << rings_to_process.number_elements() << " potential 123 triazoles in '" << m0->name() << "'\n";

  for (int i = 0; i < rings_to_process.number_elements(); i++)
  {
    const Triazole_123 * t = rings_to_process[i];

    do_123_triazoles(*t, molecules);
  }

  return 1;
}

/*
  We capture imidazoles only, just two nitrogens
*/

class Imidazole
{
  private:
    const atom_number_t _n1;
    const atom_number_t _c;
    const atom_number_t _n2;

  public:
    Imidazole (atom_number_t, atom_number_t, atom_number_t);

    int is_possible_imidazole (Molecule_Proc &);

    int process (Molecule_Proc &, resizable_array_p<Molecule_Proc> &) const;
};


Imidazole::Imidazole (atom_number_t a1, atom_number_t a2, atom_number_t a3) : _n1(a1), _c(a2), _n2(a3)
{
  return;
}

int
Imidazole::process (Molecule_Proc & m,
                    resizable_array_p<Molecule_Proc> & molecules) const
{
//const Atom * an1 = m.atomi(_n1);
//const Atom * ac  = m.atomi(_c);
//const Atom * an2 = m.atomi(_n2);

  int ih1 = m.implicit_hydrogens(_n1);
  int ih2 = m.implicit_hydrogens(_n2);

  if (0 == ih1 && 0 == ih2)
    return 0;

  if (ih1 && ih2)
    return 0;

  const Bond * b1 = m.bond_between_atoms(_c, _n1);
  const Bond * b2 = m.bond_between_atoms(_c, _n2);

  if (b1->is_single_bond() && b2->is_double_bond() && ih1)
  {
    Molecule_Proc * tmp = new Molecule_Proc(m);
    tmp->set_bond_type_between_atoms(_c, _n1, DOUBLE_BOND);
    tmp->set_bond_type_between_atoms(_c, _n2, SINGLE_BOND);
    molecules.add(tmp);
  }
  else if (b1->is_double_bond() && b2->is_single_bond() && ih2)
  {
    Molecule_Proc * tmp = new Molecule_Proc(m);
    tmp->set_bond_type_between_atoms(_c, _n1, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_c, _n2, DOUBLE_BOND);
    molecules.add(tmp);
  }
  else
    return 0;

  imidazoles_flipped++;

  return 1;
}

static int
do_imidazoles(const Imidazole & im,
              resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * m = molecules[i];

    im.process (*m, molecules);

    if (molecules.number_elements() >= max_tautomers_per_molecule)
    {
      terminated_by_max_tautomers++;
      return 1;
    }
  }

  return 1;
}

/*
  Any C-N-C-N-C site is a potential imidazole
*/

static int
identify_imidazoles_in_ring (Molecule_Proc & m,
                             const Ring & r,
                             resizable_array_p<Imidazole> & sites_to_process)
{
  int ring_size = r.number_elements();

  resizable_array<atomic_number_t> z;
  resizable_array<atom_number_t> anum;

  int free_nitrogens = 0;

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];
    anum.add(j);

    const Atom * aj = m.atomi(j);

    atomic_number_t zj = aj->atomic_number();

    if (! m.process_atom(j))
    {
      z.add(-1);
      continue;
    }

    if (6 == zj)
      z.add(6);
    else if (7 == zj)
    {
      if (2 == aj->ncon())
      {
        free_nitrogens++;
        z.add(7);
      }
      else
        z.add(0);
    }
    else
      z.add(0);
  }

#ifdef DEBUG_IMIDAZOLE
  cerr << "Atomic numbers";
  for (int i = 0; i < z.number_elements(); i++)
  {
    cerr << ' ' << z[i];
  }
  cerr << endl;
  cerr << free_nitrogens << " free nitrogens\n";
#endif

  if (2 != free_nitrogens)
    return 0;

// Now look for CNCNC sets

  int rc = 0;

  for (int i = 0; i < ring_size; i++)
  {
    atomic_number_t zi = z[i];

    if (7 == zi)    // we start on a non-nitrogen atom
      continue;

    int ndx = i;
    if (7 == z.next_after_wrap(ndx) && 6 == z.next_after_wrap(ndx) && 7 == z.next_after_wrap(ndx) && 7 != z.next_after_wrap(ndx))
    {
      int n1ndx = z.next_index_after_wrap(i);
      int cndx  = z.next_index_after_wrap(n1ndx);
      int n2ndx = z.next_index_after_wrap(cndx);

      Imidazole * id = new Imidazole(anum[n1ndx], anum[cndx], anum[n2ndx]);

      sites_to_process.add(id);
      rc++;
    }
  }

  return rc;
}

static int
do_imidazoles_nv (resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * m0 = molecules[0];

  int nr = m0->nrings();

  if (0 == nr)
    return 1;

  m0->compute_aromaticity_if_needed();

  resizable_array_p<Imidazole> sites_to_process;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m0->ringi(i);

//  if (5 != ri->number_elements())
//    continue;

    if (! ri->is_aromatic())
      continue;

    identify_imidazoles_in_ring (*m0, *ri, sites_to_process);
  }

  if (sites_to_process.empty())
    return 0;

  if (verbose > 1)
    cerr << "Found " << sites_to_process.number_elements() << " potential imidazoles in '" << m0->name() << "'\n";

  for (int i = 0; i < sites_to_process.number_elements(); i++)
  {
    const Imidazole * im = sites_to_process[i];

    do_imidazoles(*im, molecules);
  }

  return 1;
}

#ifdef NOT_USED_HERE
static int
connected_to_singly_attached_oxygen (const Molecule & m,
                                     atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t o = a->other(zatom, i);

    const Atom * oa = m.atomi(o);

    if (8 != oa->atomic_number())
      continue;

    if (1 == oa->ncon())
      return 1;
  }

  return 0;
}
#endif

/*
  anything that contains 2 or 3 nitrogens, with at least one pyrrole type
*/

#ifdef NOT_USED_HERE
static int
contains_pyrole_and_pyridine (Molecule_Proc & m,
                              const Ring & r)
{
  assert (5 == r.number_elements());

  int pyrrole_count = 0;
  int pyridine_count = 0;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    const Atom * aj = m.atomi(j);

    if (7 != aj->atomic_number())
      continue;

    if (0 != aj->formal_charge())
      continue;

    if (2 != aj->ncon())
      continue;

    if (0 == m.hcount(j))
    {
      pyridine_count++;
      continue;
    }

//  cerr << "Pyrrole is i = " << i << " atom " << j << " '" << m.smarts_equivalent_for_atom(j) << "'\n";

    atom_number_t notused1, notused2;
    if (is_imidazole_form (m, r, notused1, notused2, i))
    {
//    cerr << "Imidazole, skipping\n";
      return 0;
    }

    pyrrole_count++;
  }

  if (0 == pyrrole_count)
    return 0;

  if (0 == pyridine_count)
    return 0;

  return 1;
}

static int
do_5_membered_aromatics (resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * m0 = molecules[0];

  int nrings = m0->nrings();

  if (0 == nrings)
    return 0;

  m0->compute_aromaticity_if_needed();

  for (int i = 0; i < nrings; i++)
  {
    const Ring * ri = m0->ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (! ri->is_aromatic())
      continue;

    if (ri->is_fused())
      continue;

    if (! contains_pyrole_and_pyridine(*m0, *ri))
      continue;

    five_membered_ring_kekule_toggled++;

    Toggle_Kekule_Form tkf;
    const Bond * b = m0->bond_between_atoms (ri->item(0), ri->item(1));
    if (b->is_single_bond())
      tkf.add_bond(0, 1, DOUBLE_BOND);
    else
      tkf.add_bond(0, 1, SINGLE_BOND);

    int nmolecules = molecules.number_elements();

    for (int j = 0; j < nmolecules; j++)
    {
      Molecule_Proc * tmp = new Molecule_Proc(*(molecules[j]));

      int changed;
      tkf.process(*tmp, *ri, changed);
      if (! changed)
        delete tmp;
      else
      {
        molecules.add(tmp);
        if (molecules.number_elements() >= max_tautomers_per_molecule)
        {
          terminated_by_max_tautomers++;
          return 1;
        }
      }
    }
  }

  return molecules.number_elements();
}
#endif

/*
  We specifically make allowance for N-C(=O)-N
*/

class Keto_Enol_Form
{
  private:
    const atom_number_t _oxygen;
    const atom_number_t _carbon;
    atom_number_t _n1;
    atom_number_t _n2;

//  private functions

    int _do_flip (const Molecule_Proc & m, atom_number_t, resizable_array_p<Molecule_Proc> & molecules) const;

  public:
    Keto_Enol_Form(atom_number_t, atom_number_t);

    int is_possible_keto_enol_form (Molecule_Proc & m);

    int can_flip (const Molecule_Proc &) const;

    int do_flips (const Molecule_Proc & m, resizable_array_p<Molecule_Proc> & molecules) const;
};

Keto_Enol_Form::Keto_Enol_Form(atom_number_t o, atom_number_t c) :_oxygen(o), _carbon(c)
{
  _n1     = INVALID_ATOM_NUMBER;
  _n2     = INVALID_ATOM_NUMBER;

  return;
}

int
Keto_Enol_Form::can_flip (const Molecule_Proc & m) const
{
  if (4 != m.nbonds(_carbon))  // not sure how that could happen
    return 0;

  return 1;
}

int
Keto_Enol_Form::do_flips (const Molecule_Proc & m,
                          resizable_array_p<Molecule_Proc> & molecules) const
{
  assert (INVALID_ATOM_NUMBER != _n1);

  int rc = _do_flip(m, _n1, molecules);

  if (INVALID_ATOM_NUMBER != _n2)
    rc += _do_flip(m, _n2, molecules);

  return rc;
}

int
Keto_Enol_Form::_do_flip (const Molecule_Proc & m,
                          atom_number_t n,
                          resizable_array_p<Molecule_Proc> & molecules) const
{
  if (molecules.number_elements() >= max_tautomers_per_molecule)
  {
    terminated_by_max_tautomers++;
    return 1;
  }

  const Bond * bnc = m.bond_between_atoms(n, _carbon);
  const Bond * bco = m.bond_between_atoms(_carbon, _oxygen);

  if (bco->is_single_bond() && bnc->is_double_bond())
  {
    Molecule_Proc * tmp = new Molecule_Proc(m);
    tmp->set_bond_type_between_atoms(_carbon, n,       SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_carbon, _oxygen, DOUBLE_BOND);
    molecules.add(tmp);

    keto_enol_flipped++;

    return 1;
  }

  if (bco->is_double_bond() && bnc->is_single_bond() && 2 == m.nbonds(n))
  {
    Molecule_Proc * tmp = new Molecule_Proc(m);
    tmp->set_bond_type_between_atoms(_carbon, _oxygen, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(_carbon, n,       DOUBLE_BOND);
    molecules.add(tmp);

    keto_enol_flipped++;

    return 1;
  }

  return 0;
}

int
Keto_Enol_Form::is_possible_keto_enol_form (Molecule_Proc & m)
{
  const Atom * ac = m.atomi(_carbon);

  if (3 != ac->ncon())
    return 0;

  int rc = 0;

#ifdef DEBUG_KETO_ENOL_FORM
  cerr << "Looking for keto/enol from " << m.smarts_equivalent_for_atom(c) << endl;
#endif

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = ac->item(i);

    atom_number_t j = b->other(_carbon);

    if (_oxygen == j)
      continue;

    if (! m.process_atom(j))
      continue;

#ifdef DEBUG_KETO_ENOL_FORM
    cerr << "Examine atom " << j << " type " << m.smarts_equivalent_for_atom(j) << ", process " << m.process_atom(j) << endl;
#endif

    const Atom * aj = m.atomi(j);

    if (6 == aj->atomic_number())
      continue;

    if (7 != aj->atomic_number())     // some other heteroatom, definitely not
      return 0;

    if (2 != aj->ncon())     // no terminal or other types
      return 0;

    if (INVALID_ATOM_NUMBER == _n1)
      _n1 = j;
    else
      _n2 = j;

    rc += 1;
  }

#ifdef DEBUG_KETO_ENOL_FORM
  cerr << "looks_like_keto_enol_form returning " << rc << endl;
#endif

  return rc;
}

static int
do_keto_enol_form_nv (Molecule_Proc & m,
                      const Keto_Enol_Form & k,
                      resizable_array_p<Molecule_Proc> & molecules)
{
  return k.do_flips (m, molecules);
}

static int
do_keto_enol_forms_nv (const Keto_Enol_Form & k, 
                       resizable_array_p<Molecule_Proc> & molecules)
{
  int n = molecules.number_elements();

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * m = molecules[i];

    do_keto_enol_form_nv (*m, k, molecules);
  }

  return 1;
}

/*static int
looks_like_keto_enol_form (const Molecule_Proc & m,
                           const atom_number_t o,
                           const atom_number_t c,
                           atom_number_t & n1,    // by reference
                           atom_number_t & n2,    // by reference
                           bool & keto_to_enol)
{
  n1 = INVALID_ATOM_NUMBER;
  n2 = INVALID_ATOM_NUMBER;

  const Atom * ac = m.atomi(c);

  assert (3 == ac->ncon());

  bond_type_t bond_to_oxygen = SINGLE_BOND;   // not really using this

  int rc = 0;

#ifdef DEBUG_KETO_ENOL_FORM
  cerr << "Looking for keto/enol from " << m.smarts_equivalent_for_atom(c) << endl;
#endif

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = ac->item(i);

    atom_number_t j = b->other(c);

#ifdef DEBUG_KETO_ENOL_FORM
    cerr << "Examine atom " << j << " type " << m.smarts_equivalent_for_atom(j) << ", process " << m.process_atom(j) << endl;
#endif

    if (j == o)
    {
      if (b->is_single_bond())
        bond_to_oxygen = SINGLE_BOND;
      else
        bond_to_oxygen = DOUBLE_BOND;
      continue;
    }

    const Atom * aj = m.atomi(j);

    if (7 == aj->atomic_number())
      ;
    else if (6 == aj->atomic_number())
      continue;
    else               // some other heteroatom, definitely not
      return 0;

    if (! m.process_atom(j))
      continue;

    if (1 == aj->ncon())   // not sure, but this doesn't sound right
      return 0;

    if (2 != aj->ncon())
      return 0;

    if (b->is_double_bond())
      ;
    else if (2 != aj->nbonds())
      continue;

    if (b->is_single_bond())
      keto_to_enol = true;
    else
      keto_to_enol = false;

    if (INVALID_ATOM_NUMBER == n1)
      n1 = j;
    else
      n2 = j;

    rc += 1;
  }

#ifdef DEBUG_KETO_ENOL_FORM
  cerr << "looks_like_keto_enol_form returning " << rc << endl;
#endif

  return rc;
}*/

static int
do_keto_enol_forms_nv (resizable_array_p<Molecule_Proc> & molecules)
{
  resizable_array_p<Keto_Enol_Form> sites_to_process;

  Molecule_Proc * m0 = molecules[0];

  m0->compute_aromaticity_if_needed();

  int matoms = m0->natoms();

  for (int o = 0; o < matoms; o++)
  {
    if (! m0->process_atom(o))
      continue;

    const Atom * a = m0->atomi(o);

    if (8 != a->atomic_number())
      continue;

    if (1 != a->ncon())
      continue;

    if (0 != a->formal_charge())
      continue;

    const Bond * b1 = a->item(0);

    atom_number_t c = b1->other(o);

    if (! m0->process_atom(c))
      continue;

    if (allow_chain_amide_groups_to_flip)
      ;
    else if (m0->is_non_ring_atom(c))
      continue;

    const Atom * ac = m0->atomi(c);

    if (6 != ac->atomic_number())
      continue;

    if (3 != ac->ncon())
      continue;

    if (4 != ac->nbonds())
      continue;

    Keto_Enol_Form * tmp = new Keto_Enol_Form(o, c);
    if (! tmp->is_possible_keto_enol_form (*m0))
      delete tmp;
    else
      sites_to_process.add(tmp);
  }

  if (sites_to_process.empty())
    return 0;

  if (verbose > 1)
    cerr << "Found " << sites_to_process.number_elements() << " possible keto-enol forms\n";

  for (int i = 0; i < sites_to_process.number_elements(); i++)
  {
    const Keto_Enol_Form * k = sites_to_process[i];

    do_keto_enol_forms_nv(*k, molecules);

    if (molecules.number_elements() >= max_tautomers_per_molecule)
      break;
  }

  return 1;
}

/*static int
identify_keto_to_enol (Molecule_Proc & m,
                       atom_number_t & astart,
                       atom_number_t & o,
                       atom_number_t & c,
                       atom_number_t & n1,   // might be urea-like
                       atom_number_t & n2,
                       bool & keto_to_enol)
{
  for (int matoms = m.natoms(); astart < matoms; astart++)
  {
    if (! m.process_atom(astart))
      continue;

    const Atom * a = m.atomi(astart);

    if (8 != a->atomic_number())
      continue;

    if (1 != a->ncon())
      continue;

    if (0 != a->formal_charge())
      continue;

    const Bond * b1 = a->item(0);

    c = b1->other(astart);

    if (! m.process_atom(c))
      continue;

    if (allow_chain_amide_groups_to_flip)
      ;
    else if (m.is_non_ring_atom(c))
      continue;

    const Atom * ac = m.atomi(c);

    if (6 != ac->atomic_number())
      continue;

    if (3 != ac->ncon())
      continue;

    if (4 != ac->nbonds())
      continue;

    if (! looks_like_keto_enol_form(m, astart, c, n1, n2, keto_to_enol))
      continue;

    o = astart;
    astart++;
    return 1;
  }

  return 0;
}*/

/*static int
create_keto_enol_variant (const Molecule_Proc * m,
                          atom_number_t o,
                          atom_number_t c,
                          atom_number_t n,
                          int keto_to_enol,
                          resizable_array_p<Molecule_Proc> & molecules)
{
  Molecule_Proc * tmp = new Molecule_Proc(*m);

  if (keto_to_enol)
  {
    if (2 != m->nbonds(n))
      return 0;

    tmp->set_bond_type_between_atoms(o, c, SINGLE_BOND);
    tmp->set_bond_type_between_atoms(c, n, DOUBLE_BOND);
  }
  else
  {
    tmp->set_bond_type_between_atoms(o, c, DOUBLE_BOND);
    tmp->set_bond_type_between_atoms(c, n, SINGLE_BOND);
  }

  molecules.add(tmp);

  return 1;
}*/

/*static int
do_keto_enol_forms (resizable_array_p<Molecule_Proc> & molecules)
{
  int nmolecules = molecules.number_elements();
  assert (nmolecules > 0);

  Molecule_Proc * m = molecules[0];

  atom_number_t o, c, n1, n2;
  bool keto_to_enol;

  atom_number_t astart = 0;

  while (identify_keto_to_enol(*m, astart, o, c, n1, n2, keto_to_enol))
  {
    keto_enol_flipped++;

    for (int i = 0; i < nmolecules; i++)
    {
      create_keto_enol_variant (molecules[i], o, c, n1, keto_to_enol, molecules);

      if (INVALID_ATOM_NUMBER != n2)
        create_keto_enol_variant (molecules[i], o, c, n2, keto_to_enol, molecules);
    }

    nmolecules = molecules.number_elements();
  }

  return molecules.number_elements();
}*/

static int
identify_doubly_bonded_oxygen (const Molecule_Proc & m,
                               atom_number_t c,
                               atom_number_t & doubly_bonded_oxygen)
{
  doubly_bonded_oxygen = INVALID_ATOM_NUMBER;

  const Atom * a = m.atomi(c);

  if (3 != a->ncon())
    return 0;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(c);

    if (! m.process_atom(o))
      continue;

    if (8 != m.atomic_number(o))
      continue;

    doubly_bonded_oxygen = o;
    return 1;
  }

  return 0;
}

static int
identify_acid (const Molecule_Proc & m,
               atom_number_t & astart,
               atom_number_t & o1,
               atom_number_t & c,
               atom_number_t & o2)
{
  for (int matoms = m.natoms(); astart < matoms; astart++)
  {
    if (! m.process_atom(astart))
      continue;

    const Atom * a = m.atomi(astart);

    if (8 != a->atomic_number())
      continue;

    if (1 != a->ncon())
      continue;

    if (1 != a->nbonds())
      continue;

    c = a->other(astart, 0);

    if (! m.process_atom(c))
      continue;

    if (! identify_doubly_bonded_oxygen(m, c, o2))
      continue;

    o1 = astart;
    astart++;
    return 1;
  }

  return 0;
}

static int
do_acid_forms (resizable_array_p<Molecule_Proc> & molecules)
{
  int nmolecules = molecules.number_elements();
  assert (nmolecules > 0);

  const Molecule_Proc * m = molecules[0];

  atom_number_t o1, c, o2;

  atom_number_t astart = 0;

  while (identify_acid(*m, astart, o1, c, o2))
  {
    acids_flipped++;

    for (int i = 0; i < nmolecules; i++)
    {
      Molecule_Proc * tmp = new Molecule_Proc(*(molecules[i]));

      tmp->set_bond_type_between_atoms(o2, c, SINGLE_BOND);
      tmp->set_bond_type_between_atoms(o1, c, DOUBLE_BOND);
      molecules.add(tmp);
      
      if (molecules.number_elements() > max_tautomers_per_molecule)
      {
        terminated_by_max_tautomers++;
        break;
      }
    }

    nmolecules = molecules.number_elements();

    if (nmolecules >= max_tautomers_per_molecule)
      break;
  }

  return molecules.number_elements();
}

static void
append_sequence_number_to_name (Molecule & m,
                                const IWString & mname,
                                int seq)
{
  if (0 == seq)
    return;

  IWString tmp(mname);
  tmp << prepend_before_sequence << seq;

  m.set_name(tmp);

  return;
}

static int
handle_aromaticty_loss (Molecule & m,
                        int seq)
{
  rejected_for_aromaticity_loss++;
  if (stream_for_invalid_valence.active())
  {
    IWString tmp = m.name();
    tmp << " AROM";
    m.set_name(tmp);
    stream_for_invalid_valence.write(m);
  }

  return 1;
}

static int
handle_changed_hcount (Molecule & m,
                        int seq)
{
  rejected_for_hydrogen_change++;

  return 1;
}

static int
handle_bad_valence (Molecule & m,
                    int seq)
{
  invalid_valence_rejected++;

  if (stream_for_invalid_valence.active())
  {
    append_sequence_number_to_name(m, m.name(), seq);

    stream_for_invalid_valence.write(m);
  }

  return 1;
}

static int 
is_unique_2d (IW_STL_Hash_Set & produced,
              const IWString & usmi)
{
  if (produced.contains(usmi))
    return 0;

  produced.insert(usmi);

  return 1;
}

static int
is_unique_3d (const resizable_array_p<Molecule_Proc> & molecules,
              int n,
              Molecule_Proc & m)
{
  for (int i = 0; i < n; i++)
  {
    if (m == *(molecules[i]))
      cerr << "Discarding duplicate molecule '" << m.name() << "'\n";
    if (m == *(molecules[i]))
      return 0;
  }

  return 1;
}

static int
tautomer_generation (resizable_array_p<Molecule_Proc> & molecules,
                     Molecule_Output_Object & output)
{
  assert (1 == molecules.number_elements());

//if (keto_enol_form)
//  do_keto_enol_forms(molecules);

  if (process_keto_enol)
    do_keto_enol_forms_nv(molecules);

  if (! process_acid)
    ;
  else if (2 == filter_to_unique_structures)   // no sense doing acids in 2D
    ;
  else
    do_acid_forms(molecules);

  if (process_imidazole)
    do_imidazoles_nv (molecules);
  if (process_triazole_124)
    do_isolated_124_triazoles (molecules);
  if (process_triazole_123)
    do_123_triazoles (molecules);
  if (process_tetrazole)
    do_tetrazoles (molecules);
  if (process_pyrazole)
    do_pyrazoles (molecules);

// do_5_membered_aromatics (molecules);

  int n = molecules.number_elements();

  tautomers_produced += n;

  tautomer_stats[n]++;

  const IWString mname = molecules[0]->name();

  int arom;
  if (reject_for_aromaticity_loss)
    arom = molecules[0]->aromatic_atom_count();
  else
    arom = 0;

  int hcount;

  if (must_preserve_hcount)
    hcount = molecules[0]->implicit_hydrogens();
  else
    hcount = -1;

  IW_STL_Hash_Set produced;

  int seq = 0;

  for (int i = 0; i < n; i++)
  {
    Molecule_Proc * mi = molecules[i];

    if (0 == arom)   // either not checking or none in starting molecule
      ;
    else if (mi->aromatic_atom_count() < arom)
    {
      handle_aromaticty_loss(*mi, i);
      continue;
    }

    if (hcount < 0)
      ;
    else if (mi->implicit_hydrogens() != hcount)
    {
      handle_changed_hcount(*mi, i);
      continue;
    }

    if (reject_invalid_valence && ! mi->valence_ok())
    {
      handle_bad_valence(*mi, i);
      continue;
    }

    if (2 == filter_to_unique_structures)
    {
      if (! is_unique_2d(produced, mi->unique_smiles()))
        continue;

      mi->invalidate_smiles();
    }
    else if (3 == filter_to_unique_structures)
    {
      if (! is_unique_3d(molecules, i, *mi))
        continue;
    }

    append_sequence_number_to_name(*mi, mname, seq++);

    if (! output.write(mi))
      return 0;

    molecules_written++;
  }

  return output.good ();
}

static int
tautomer_generation_2 (resizable_array_p<Molecule_Proc> & molecules,
                       Molecule_Output_Object & output)
{
  Molecule_Proc * m = molecules[0];

  int matoms = m->natoms();

  if (0 == matoms)
  {
    cerr << "Skipping empty molecule '" << m->name() << "'\n";
    return 1;
  }

  if (0 == filter_to_unique_structures)
  {
    int d = m->highest_coordinate_dimensionality();
    if (3 == d)
      filter_to_unique_structures = 3;
    else
      filter_to_unique_structures = 2;
  }

  if (reject_invalid_valence && ! m->valence_ok())
    cerr << "Warning, rejecting invalid valences but initial molecule bad '" << m->name() << "'\n";

  if (0 == m->initialise (queries))
  {
    if (verbose)
      cerr << "No atoms to be processed in '" << m->name() << "'\n";

    return output.write(*m);
  }

  return tautomer_generation (molecules, output);
}

static int
tautomer_generation (data_source_and_type<Molecule_Proc> & input,
                     Molecule_Output_Object & output)
{
  Molecule_Proc * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    preprocess(*m);

    resizable_array_p<Molecule_Proc> molecules;
    molecules.add(m);

    if (! tautomer_generation_2 (molecules, output))
      return 0;
  }

  return 1;
}

static int
tautomer_generation (const char * fname, FileType input_type, 
                     Molecule_Output_Object & output)
{
  assert (nullptr != fname);

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule_Proc> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return tautomer_generation(input, output);
}

static int
tautomer_generation (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lT:s:q:o:V:au:X:mh:IS:Y:");

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
  else
    set_global_aromaticity_type(Daylight);

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

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring s;

    while (cl.value('s', s, i++))
    {
      Substructure_Query * q = new Substructure_Query;
      if (! q->create_from_smarts(s))
      {
        cerr << "Invalid smarts '" << s << "'\n";
        return 4;
      }

      queries.add(q);
    }
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose, 'q'))
    {
      cerr << "Cannot process queries (-q)\n";
      return 4;
    }
  }

  if (cl.option_present('h'))
  {
    must_preserve_hcount = 1;

    if (verbose)
      cerr << "Tautomers must preserve hydrogen count\n";
  }

  if (cl.option_present('T'))
  {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('T', t, i++))
    {
      if (t.starts_with("KE"))
      {
        process_keto_enol = 1;
      }
      else if (t.starts_with("AC"))
      {
        process_acid = 1;
      }
      else if ("nH" == t)
      {
        process_imidazole = 1;
        process_pyrazole = 1;
        process_triazole_123 = 1;
        process_triazole_124 = 1;
        process_tetrazole = 1;
      }
      else if (t.starts_with("IM"))
      {
        process_imidazole = 1;
      }
      else if (t.starts_with("PY"))
      {
        process_pyrazole = 1;
      }
      else if (t.starts_with("T123"))
      {
        process_triazole_123 = 1;
      }
      else if (t.starts_with("T124"))
      {
        process_triazole_124 = 1;
      }
      else
      {
        cerr << "Unrecognised -T qualifier '" << t << "'\n";
        usage(3);
      }
    }
  }
  else
  {
    process_acid = 1;
    process_keto_enol = 1;
    process_imidazole = 1;
    process_pyrazole = 1;
    process_triazole_123 = 1;
    process_triazole_124 = 1;
    process_tetrazole = 1;
  }

  if (cl.option_present('m'))
  {
    allow_chain_amide_groups_to_flip = 1;
    if (verbose)
      cerr << "Chain amide groups can flip\n";
  }

  if (cl.option_present('a'))
  {
    reject_for_aromaticity_loss = 1;

    if (verbose)
      cerr << "Will reject tautomers that lose aromatic form\n";
  }

  if (cl.option_present('I'))
  {
    apply_isotopic_labels = 1;

    if (verbose)
      cerr << "Will apply isotopic labels to molecules\n";
  }

  if (cl.option_present('V'))
  {
    const_IWSubstring v = cl.string_value('V');

    if ('.' == v)
      ;
    else
    {
      stream_for_invalid_valence.add_output_type(FILE_TYPE_SMI);
      if (stream_for_invalid_valence.would_overwrite_input_files(cl, v))
      {
        cerr << "Invalid valence rejection stream '" << v << "' cannot overwrite input file(s)\n";
        return 4;
      }

      if (! stream_for_invalid_valence.new_stem(v))
      {
        cerr << "Cannot initialise invalid valence stream to '" << v << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Invalid valence molecules written to '" << v << "'\n";
    }

    reject_invalid_valence = 1;

    if (verbose)
      cerr << "Will discard molecules with invalid valences\n";
  }

  if (cl.option_present('X'))
  {
    if (! cl.value('X', max_tautomers_per_molecule) || max_tautomers_per_molecule < 1)
    {
      cerr << "The maximum number of tautomers per molecule (-X) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will generate a max of " << max_tautomers_per_molecule << " tautomers per molecule\n";
  }

  if (cl.option_present('u'))
  { 
    const_IWSubstring u = cl.string_value('u');

    if ("3d" == u || '3' == u)
      filter_to_unique_structures = 3;
    else if ("2d" == u || '2' == u)
      filter_to_unique_structures = 2;
    else if ("none" == u)
      filter_to_unique_structures = -1;
    else
    {
      cerr << "The filter to unique structures option (-u) must be either 2 or 3\n";
      usage(3);
    }

    if (0 == verbose)
      ;
    else if (2 == filter_to_unique_structures)
      cerr << "Will filter to unique structures based on unique smiles\n";
    else if (3 == filter_to_unique_structures)
      cerr << "Will filter to unique structures based on atom number\n";
    else
      cerr << "No filtering to unique structures\n";
  }

  if (cl.option_present('Y'))
  {
    cl.value('Y', prepend_before_sequence);

    if (verbose)
      cerr << "Will insert '" << prepend_before_sequence << "' before sequence number\n";
  }

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SMI);
  else if (! output.determine_output_types(cl))
  {
    cerr << "Cannot determine output type(s), (-o)\n";
    usage(4);
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 4;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot initialise output streams to '" << s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "output written to '" << s << "'\n";
  }
  else if (! output.new_stem("-"))
  {
    cerr << "Gack, cannot redirect to stdout!!!\n";
    return 4;
  }

  set_copy_name_in_molecule_copy_constructor(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! tautomer_generation(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.do_flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, " << (tautomers_produced - molecules_read) << " tautomers produced\n";
    for (int i = 0; i < tautomer_stats.number_elements(); i++)
    {
      if (tautomer_stats[i])
        cerr << tautomer_stats[i] << " molecules had " << i << " tautomers\n";
    }
    cerr << "Wrote " << molecules_written << " molecules\n";

    cerr << keto_enol_flipped << " keto enol forms flipped\n";
    cerr << acids_flipped << " acids flipped\n";
    cerr << pyrazoles_flipped << " pyrazoles flipped\n";
    cerr << imidazoles_flipped << " imidazole type groups flipped\n";
    cerr << triazole_123_flipped << " triazoles 123 flipped\n";
    cerr << triazole_124_flipped << " triazoles 124 flipped\n";
    cerr << tetrazoles_flipped << " tetrazoles 124 flipped\n";
    cerr << five_membered_ring_kekule_toggled << " five membered aromatic rings kekule toggled\n";

    if (reject_invalid_valence)
      cerr << invalid_valence_rejected << " tautomers with invalid valences rejected\n";
    if (reject_for_aromaticity_loss)
      cerr << rejected_for_aromaticity_loss << " tautomers with aromaticity loss rejected\n";

    if (terminated_by_max_tautomers)
      cerr << terminated_by_max_tautomers << " terminations for having " << max_tautomers_per_molecule << " or more tautomers\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tautomer_generation (argc, argv);

  return rc;
}
