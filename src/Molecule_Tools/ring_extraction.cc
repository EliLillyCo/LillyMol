/*
  Extracts rings from molecules

  We produce the replacement smiles and the corresponding smarts. The atom ordering
  in these must be identical

  Sept 2013.

  There is a fatal flaw in here.
  We extract rings and ring systems. We also include doubly bonded exocyclic atoms.
  Because we need to keep the smiles and the smarts in sync, we have to put those atoms
  in he smarts too, but that means those rings will only match original rings
  that have that same arrangement.

  Bad, but too hard to fix
*/

#include <stdlib.h>
#include <memory>
using namespace std;

#include "cmdline.h"
#include "misc.h"

#define CHECK_RESDULT
#ifdef CHECK_RESDULT
#include "substructure.h"
#endif

#include "istream_and_type.h"
#include "path.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "reaction_duplicate.h"

#include "ring_ext_rep.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static extending_resizable_array<int> rings_per_molecule;

static int molecules_with_rings = 0;

static Reaction_Duplicate reaction_duplicate;

static int unique_rings_only = 0;

class Ring_Extraction_Replacement_Conditions rerc;

static int must_match_whole_ring = 1;

static int min_heteroatoms_needed_in_ring = 0;

static int extract_all_isolated_rings = 0;
static int extract_all_fused_systems = 0;

static int min_ring_size = 4;
static int max_ring_size = 7;

/*
  Must keep the numbers here consistent with what is in the header file
*/

static std::ofstream output_stream[RING_ARRAY_SIZE];
static IWString suffix[RING_ARRAY_SIZE];

static int check_smiles_formed = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Extracts rings from molecules\n";
  display_standard_ring_ext_rep_options(cerr);
  cerr << "  -h <number>   rings extracted must contain at least <number> heteroatoms\n";
  cerr << "  -N <stem>     get all rings, write to files starting with <stem>\n";
  cerr << "  -u            unique rings only\n";
  cerr << "  -V            check the validity of smiles formed\n";
  cerr << "  -R <size>     maximum ring size to process (default " << max_ring_size << ")\n";
  cerr << "  -4            allow four membered aromatics\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
append_chain_smiles (Molecule & m,
                     const int * include_atom,
                     atom_number_t zatom, 
                     const atom_number_t previous_atom,
                     IWString & s)
{
  s << m.atomic_symbol(zatom);

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int branches_so_far = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (! include_atom[j])
      continue;

    if (previous_atom == j)
      continue;

    if (branches_so_far)
      s << '(';

    if (b->is_double_bond())
      s << '=';
    else if (b->is_triple_bond())
      s << '#';

    append_chain_smiles(m, include_atom, j, zatom, s);

    if (branches_so_far)
      s << ')';

    branches_so_far++;
  }

  return 1;
}

static void
append_smiles_symbol (Molecule & m,
                      atom_number_t current_atom,
                      IWString & s)
{
  formal_charge_t fc = m.formal_charge(current_atom);

  if (m.isotope(current_atom))
  {
    s << '[' << m.isotope(current_atom) << m.atomic_symbol(current_atom);
    if (0 == fc)
      ;
    else if (fc > 0)
      s << '+';
    else
      s << '-';

    s << ']';
  }
  else if (0 == fc)
    s << m.atomic_symbol(current_atom);
  else if (fc < 0)
    s << '[' << m.atomic_symbol(current_atom) << "-]";
  else
    s << '[' << m.atomic_symbol(current_atom) << "+]";

  return;
}

static int
build_smiles (Molecule & m,
              atom_number_t first_ring_atom,
              const int * include_atom,
              int ring_is_aromatic,
              const Ring_Extraction_Replacement_Conditions & rerc,
              int build_smarts,
              IWString & s)
{
  atom_number_t current_atom = first_ring_atom;
  atom_number_t previous_atom = INVALID_ATOM_NUMBER;

  while (1)
  {
    if (build_smarts)
      rerc.append_connectivity_smarts(m, current_atom, include_atom, ring_is_aromatic, s);
    else
      append_smiles_symbol(m, current_atom, s);

    if (INVALID_ATOM_NUMBER != previous_atom)    // already written the first atom
      ;
    else if (! must_match_whole_ring)    // don't put in ring symbols
      ;
    else if (build_smarts)           // use a very unusual ring symbol
      s << "%99";
    else
      s << '1';

    const Atom * a = m.atomi(current_atom);

    int acon = a->ncon();

//  cerr << "Atom " << current_atom << " has " << acon << " neighbours, smiles so var '" << s << "'\n";

    Set_of_Atoms neighbours;    // make sure we take any double bonds first
    Set_of_Atoms doubly_bonded_outside_ring;

    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(i);

      atom_number_t j = b->other(current_atom);

      if (j == previous_atom)
        continue;

      if (! include_atom[j])
        continue;

      if (0 == b->nrings() && b->is_double_bond())
        doubly_bonded_outside_ring.add(i);
      else if (b->is_aromatic())
        neighbours.add(i);
      else if (0 == b->nrings())
        neighbours.add(i);
      else if (b->is_double_bond())    // ring double bond
        neighbours.insert_at_beginning(i);
      else
       neighbours.add(i);
    }

    for (int j = 0; j < doubly_bonded_outside_ring.number_elements(); j++)
    {
      neighbours.insert_at_beginning(doubly_bonded_outside_ring[j]);
    }

    atom_number_t next_atom = INVALID_ATOM_NUMBER;

    acon = neighbours.number_elements();

//  cerr << "Acon now " << acon << " atoms " << neighbours << endl;

    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(neighbours[i]);

      atom_number_t j = b->other(current_atom);

      if (b->nrings())
      {
        if (INVALID_ATOM_NUMBER != next_atom)
          continue;

        if (build_smarts)
        {
          if (b->is_aromatic())
            s << ':';
          else if (b->is_double_bond())
            s << '=';
        }
        else if (b->is_double_bond())
          s << '=';

        if (j == first_ring_atom)   // finished the ring, put any necessary ring closure info
        {
          if (! must_match_whole_ring)
            ;
          else if (build_smarts)
            s << "%99";
          else
            s << '1';
          continue;
        }

        next_atom = j;
      }
      else if (b->is_double_bond())
        s << "(=" << m.atomic_symbol(j) << ')';
      else
      {
        s << '(';
        if (b->is_double_bond())
          s << '=';

        append_chain_smiles(m, include_atom, j, current_atom, s);

        s << ')';
      }
    }

    if (INVALID_ATOM_NUMBER == next_atom)
      break;

    previous_atom = current_atom;
    current_atom = next_atom;
  }
  
  return 1;
}

static int
ring_is_duplicate (Molecule & m,
                   const int * include_atom)
{
  Molecule tmp;
  m.create_subset(tmp, include_atom);

  return reaction_duplicate.is_duplicate(tmp);
}

static int
write_ring (Molecule & m,
            atom_number_t first_ring_atom,
            const Ring_Extraction_Replacement_Conditions & rerc,
            const int * include_atom,
            const int ring_is_aromatic,
            std::ostream & output)
{
  if (unique_rings_only)
  {
    if (ring_is_duplicate(m, include_atom))
      return 1;
  }

//cerr << "write_ring, aromatic? " << ring_is_aromatic << endl;

  m.compute_aromaticity_if_needed();

  IWString smiles;
  build_smiles(m, first_ring_atom, include_atom, 0, rerc, 0, smiles);    // no aromaticity in the smiles

  if (check_smiles_formed)
  {
    Molecule notused;
    if (! notused.build_from_smiles(smiles))
    {
      cerr << "Yipes, '" << m.name() << " produced invalid smiles\n";
      cerr << m.smiles() << ' ' << m.name() << endl;
      cerr << smiles << endl;
      return 0;
    }
  }

  IWString smarts;
  build_smiles(m, first_ring_atom, include_atom, ring_is_aromatic, rerc, 1, smarts);

  output << smiles << ' ' << smarts << ' ' << m.name() << "\n";

  return output.good();
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (rerc.remove_chirality())
    m.remove_all_chiral_centres();

  return;
}

static int
contains_enough_heteroatoms (const Molecule & m,
                             const Set_of_Atoms & r,
                             int min_heteroatoms_needed_in_ring)
{
  int n = r.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    if (6 == m.atomic_number(r[i]))
      continue;

    rc++;
    if (rc >= min_heteroatoms_needed_in_ring)
      return 1;
  }

  return 0;
}

static void
append_smiles_bond_symbol (const Bond * b,
                           IWString & smiles)
{
  if (b->is_double_bond())
    smiles << '=';
  else if (b->is_single_bond())
    ;
  else if (b->is_triple_bond())    // hard to imagine
    smiles << '#';

  return;
}

static void
append_smarts_bond_symbol (const Bond * b,
                           IWString & smarts)
{
  if (b->is_aromatic())
    smarts << ':';
  else
    append_smiles_bond_symbol(b, smarts);
}


static int
write_ring_system (Molecule & m,
                   const int * aromatic,
                   Ring_Extraction_Replacement_Conditions & rerc,
                   int * completed,
                   const int stream_index)
{
  m.compute_aromaticity_if_needed();

  int istart = -1;

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (1 == m.nrings(i))
    {
      istart = i;
      break;
    }
  }

  if (istart < 0)
  {
    cerr << "write_ring_system:huh, did not find any atom in two rings!\n";
    return 0;
  }

  set_vector(completed, matoms, 0);

  IWString smiles, smarts;

  rerc.append_connectivity_smarts(m, istart, aromatic[istart], smarts);
  smiles << m.elementi(istart)->symbol();

  smarts << "%99";
  smiles << "2";

  atom_number_t prev = INVALID_ATOM_NUMBER;

  completed[istart] = 1;

//#define DEBUG_RING_SYSTEM_SMARTS
#ifdef DEBUG_RING_SYSTEM_SMARTS
  cerr << "Starting point is " << istart << ", starting smarts '" << smarts << "'\n";
#endif

  while (INVALID_ATOM_NUMBER != istart)
  {
    const Atom * a = m.atomi(istart);

    const auto acon = a->ncon();

#ifdef DEBUG_RING_SYSTEM_SMARTS
    cerr << "Current atom is " << istart << ", ncon " << acon << " smiles '" << smiles << "'\n";
#endif

//  first append any exocyclic or cross ring bonds

    if (acon > 2)
    {
      const Bond * across_ring = NULL;

      for (auto j = 0; j < acon; ++j)
      {
        const Bond * b = a->item(j);
  
        if (2 == b->nrings())
        {
          across_ring = b;
          continue;
        }

        if (b->nrings())
          continue;

        const auto k = b->other(istart);

        smarts << '(';
        smiles << '(';
        append_smiles_bond_symbol(b, smiles);
        append_smarts_bond_symbol(b, smarts);
        append_smiles_symbol(m, k, smiles);
        rerc.append_connectivity_smarts(m, k, 0, smarts);
        smarts << ')';
        smiles << ')';
      }

      if (across_ring)
      {
        append_smiles_bond_symbol(across_ring, smiles);
        append_smarts_bond_symbol(across_ring, smarts);
        smarts << "%98";
        smiles << '1';
      }
    }

    atom_number_t next_atom = INVALID_ATOM_NUMBER;

#ifdef DEBUG_RING_SYSTEM_SMARTS
    cerr << "Looking for next atom from " << istart << " smarts so far " << smarts << endl;
#endif

    for (auto j = 0; j < acon; ++j)    // continue around the outside of the ring
    {
      const Bond * b = a->item(j);

#ifdef DEBUG_RING_SYSTEM_SMARTS
      cerr << "Checking " << b->other(istart) << " rings " << b->nrings() << " complete " << completed[b->other(istart)] << " prev? " << (b->other(istart) == prev) << endl;
#endif

      if (1 != b->nrings())
        continue;

      const auto k = b->other(istart);

      if (k == prev)
        continue;


      if (completed[k])     // must be back to start atom
      {
        append_smiles_bond_symbol(b, smiles);
        append_smarts_bond_symbol(b, smarts);
        smiles << '2';
        smarts << "%99";
      }
      else if (INVALID_ATOM_NUMBER == next_atom)
      {
        append_smiles_bond_symbol(b, smiles);
        append_smarts_bond_symbol(b, smarts);
        rerc.append_connectivity_smarts(m, k, aromatic[k], smarts);
        append_smiles_symbol(m, k, smiles);

        next_atom = k;
      }
    }

#ifdef DEBUG_RING_SYSTEM_SMARTS
    cerr << "from " << istart << " next atom will be " << next_atom << ", smarts '" << smarts << "'\n";
#endif
    prev = istart;
    istart = next_atom;
  }

  std::ostream & output = output_stream[stream_index];

  output  << m.unique_smiles() << ' ' << smiles << ' ' << smarts << ' ' << m.name() << suffix[stream_index] << "\n";

  if (check_smiles_formed)
  {
    Molecule notused;

    if (! notused.build_from_smiles(smiles))
    {
      cerr << "Invalid smiles created '" << m.smiles() << ' ' << m.name() << endl;
      cerr << smiles << endl;
      return 0;
    }

    Substructure_Query q;
    if (! q.create_from_smarts(smarts))
    {
      cerr << "INvalid smarts created '" << smarts << "'\n";
      if (0 == q.substructure_search(&m))
      {
        cerr << "No substructure match '" << m.name() << "'\n";
        cerr << smarts << endl;
        return 0;
      }
    }
  }

  return 1;
}

static int
ok_fused_system_ring_size (const Ring & r,
                           extending_resizable_array<int> & arom,
                           extending_resizable_array<int> & aliph)
{
  const auto rsize = r.number_elements();

  if (6 == rsize)
    ;
  else if (5 == rsize)
    ;
  else if (4 == rsize)
    ;
  else if (7 == rsize)
    ;
  else
    return 0;

  if (r.is_aromatic())
    arom[rsize]++;
  else
    aliph[rsize]++;

  return 1;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
extract_fused_system (Molecule & m,
                      int * ring_already_done,
                      Ring_Extraction_Replacement_Conditions & rerc,
                      const int rstart,
                      int * include_atom,
                      std::ostream & output)
{
  extending_resizable_array<int> arom, aliph;

  set_vector(include_atom, m.natoms(), 0);

  int fsid = -1;
  int fss = 0;
  bool oksystem = true;

#ifdef DEBUG_RING_SYSTEM_SMARTS
  cerr << "Ring system starts with ring " << rstart << endl;
  iw_write_array(ring_already_done, m.nrings(), " ring already done", cerr);
#endif

  const auto nr = m.nrings();
  for (auto i = rstart; i < nr; ++i)
  {
    if (ring_already_done[i])
      continue;

    const auto ri = m.ringi(i);

    if (! ri->is_fused())
      continue;

#ifdef DEBUG_RING_SYSTEM_SMARTS
    cerr << " i = " << i << " ring " << (*ri) << endl;
#endif

    if (fsid < 0)
      fsid = ri->fused_system_identifier();
    else if (ri->fused_system_identifier() != fsid)
      continue;

    ring_already_done[i] = 1;
    fss++;

    if (ri->strongly_fused_ring_neighbours())
    {
      oksystem = false;
      continue;
    }

    if (! ok_fused_system_ring_size(*ri, arom, aliph))
    {
      oksystem = false;
      continue;
    }

    ri->set_vector(include_atom, 1);
  }

  if (2 != fss || ! oksystem)    // presumably the system contained rings of the wrong size
    return 0;

  int ndx = -1;
  if (2 == arom[6])
    ndx = FUSED_6a6a;
  else if (1 == arom[5] && 1 == arom[6])
    ndx = FUSED_5a6a;
  else if (2 == arom[5])
    ndx = FUSED_5a5a;
  else if (2 == aliph[5])
    ndx = FUSED_5A5A;
  else if (2 == aliph[6])
    ndx = FUSED_6A6A;
  else if (1 == arom[6]  && 1 == aliph[6])
    ndx = FUSED_6a6A;
  else if (1 == arom[5]  && 1 == aliph[6])
    ndx = FUSED_5a6A;
  else if (1 == arom[5] && 1 == aliph[5])
    ndx = FUSED_5a5A;
  else if (1 == aliph[5] && 1 == aliph[6])
    ndx = FUSED_5A6A;
  else if (1 == aliph[5] && 1 == arom[6])
    ndx = FUSED_5A6a;
  else if (1 == aliph[4] && 1 == aliph[6])
    ndx = FUSED_4A6A;
  else if (1 == aliph[4] && 1 == arom[6])
    ndx = FUSED_4A6a;
  else if (1 == aliph[4] && 1 == aliph[5])
    ndx = FUSED_4A5A;
  else if (1 == arom[6]  && 1 == aliph[7])
    ndx = FUSED_6a7A;
  else if (1 == aliph[5] && 1 == aliph[7])
    ndx = FUSED_5A7A;
  else if (1 == arom[6]  && 1 == arom[7])
    ndx = FUSED_6a7a;
  else if (1 == aliph[6] && 1 == aliph[7])
    ndx = FUSED_6A7A;
  else if (1 == arom[5]  && 1 == arom[7])
    ndx = FUSED_5a7a;
  else if (1 == arom[5]  && 1 == aliph[7])
    ndx = FUSED_5a7A;
  else if (1 == aliph[5]  && 1 == arom[7])
    ndx = FUSED_5A7a;
  else if (1 == arom[5] && 1 == aliph[8])
    ndx = FUSED_5a8A;
  else if (1 == arom[6] && 1 == aliph[8])
    ndx = FUSED_6a8A;
  else if (1 == aliph[4] && 1 == aliph[8])
    ndx = FUSED_4A8A;
  else if (1 == aliph[5] && 1 == aliph[8])
    ndx = FUSED_5A8A;
  else if (1 == aliph[6] && 1 == aliph[8])
    ndx = FUSED_6A8A;
  else if (1 == aliph[7] && 1 == aliph[8])
    ndx = FUSED_7A8A;
  else if (1 == aliph[4] && 1 == aliph[7])
    ndx = FUSED_4A7A;
  else if (1 == aliph[8] && 1 == aliph[8])
    ndx = FUSED_8A8A;
  else if (1 == aliph[6] && 1 == arom[7])
    ndx = FUSED_6A7a;
  else if (1 == arom[6] && 1 == aliph[7])
    ndx = FUSED_6a7A;
  else if (1 == arom[4] && 1 == arom[6])
    ndx = FUSED_4a6a;
  else if (2 == aliph[4])
    ndx = FUSED_4A4A;
  else if (1 == arom[7] && 1 == aliph[8])
    ndx = FUSED_7a8A;
  else if (1 == aliph[4] && 1 == arom[5])
    ndx = FUSED_4A5a;
  else if (1 == aliph[4] && 1 == arom[7])
    ndx = FUSED_4A7a;
  else if (1 == aliph[7] && 1 == arom[7])
    ndx = FUSED_7a7A;
  else if (2 == aliph[7])
    ndx = FUSED_7A7A;
  else if (2 == arom[7])
    ndx = FUSED_7a7a;
  else if (2 == arom[4])
    ndx = FUSED_4a4a;
  else if (1 == arom[4] && 1 == aliph[6])
    ndx = FUSED_4a6A;
  else if (1 == arom[4] && 1 == arom[7])
    ndx = FUSED_4a7a;
  else if (1 == arom[4] && 1 == arom[5])
    ndx = FUSED_4a5a;
  else if (1 == arom[4] && 1 == aliph[4])
    ndx = FUSED_4a4A;
  else if (1 == arom[4] && 1 == aliph[7])
    ndx = FUSED_4a7A;
  else
  {
    cerr << "Unknown ring size combination in " << m.smiles() << ' ' << m.name() << "\n";
    for (auto i = 0; i < arom.number_elements(); ++i)
    {
      if (arom[i])
        cerr << " arom[" << i << "] = " << arom[i] << endl;
    }
    for (auto i = 0; i < aliph.number_elements(); ++i)
    {
      if (aliph[i])
        cerr << " aliph[" << i << "] = " << aliph[i] << endl;
    }
    return 0;
  }

//cerr << "stream index is " << ndx << endl;

  if (! rerc.identify_atoms_associated_with_ring_system(m, include_atom))
    return 0;

  const auto matoms = m.natoms();

  int * aromatic = new int[matoms]; std::unique_ptr<int[]> free_aromatic(aromatic);
  for (auto i = 0; i < matoms; ++i)
  {
    if (m.is_aromatic(i))
      aromatic[i] = 1;
    else
      aromatic[i] = 0;
  }

  Molecule mcopy;
  m.create_subset(mcopy, include_atom);
  mcopy.set_name(m.name());

  int aromcount = 0;
  for (auto i = 0; i < matoms; ++i)
  {
    if (include_atom[i])
    {
      aromatic[aromcount] = aromatic[i];
      aromcount++;
    }
  }

  return write_ring_system(mcopy, aromatic, rerc, include_atom, ndx);   // the include atom array now used for different purposes
}

static int
extract_fused_systems (Molecule & m,
                       Ring_Extraction_Replacement_Conditions & rerc,
                       int * include_atom,
                       std::ostream & output)
{
  const auto nr = m.nrings();

  if (nr < 2)
    return 1;

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  for (auto i = 0; i < nr; ++i)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

//  cerr << "What about ring " << *ri << ", fused " << ri->is_fused() << endl;

    if (! ri->is_fused())
      continue;

    if (! extract_fused_system(m, ring_already_done, rerc, i, include_atom, output))
      return 1;    // we always return 1
  }

  return 1;
}

static int
write_ring (Molecule & m,
            const IWString & mname,
            const Ring & r,
            int * completed)
{
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  atom_number_t istart = INVALID_ATOM_NUMBER;

  for (auto i = 0; i < matoms; ++i)
  {
    if (2 == m.ncon(i))
    {
      istart = i;
      break;
    }
  }

  if (INVALID_ATOM_NUMBER == istart)    // too strange to worry about
    return 0;

//cerr << "Processing ring " << r << endl;

  const auto aromatic = r.is_aromatic();

  set_vector(completed, m.natoms(), 0);

  IWString smiles, smarts;

  append_smiles_symbol(m, istart, smiles);
  rerc.append_connectivity_smarts(m, istart, aromatic, smarts);

  smiles << '1';
  smarts << "%99";

  atom_number_t prev = INVALID_ATOM_NUMBER;

  while (INVALID_ATOM_NUMBER != istart)
  {
    completed[istart] = 1;

//  cerr << "Procesing atom " << istart << endl;

    const Atom * a = m.atomi(istart);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j)    // first any exocyclic bonds
    {
      const Bond * b = a->item(j);

      if (b->nrings())
        continue;

      const auto k = b->other(istart);

      smarts << '(';
      smiles << '(';
      append_smiles_bond_symbol(b, smiles);
      append_smarts_bond_symbol(b, smarts);
      rerc.append_connectivity_smarts(m, k, aromatic, smarts);
      append_smiles_symbol(m, k, smiles);
      smarts << ')';
      smiles << ')';
    }

    for (auto j = 0; j < acon; ++j)     // now continue around the ring
    {
      const Bond * b = a->item(j);

      if (0 == b->nrings())
        continue;

      const auto k = b->other(istart);

      if (k == prev)
        continue;

      if (completed[k])
      {
        append_smiles_bond_symbol(b, smiles);
        append_smarts_bond_symbol(b, smarts);
        istart = INVALID_ATOM_NUMBER;
        break;
      }

      append_smiles_bond_symbol(b, smiles);
      append_smarts_bond_symbol(b, smarts);
      append_smiles_symbol(m, k, smiles);
      rerc.append_connectivity_smarts(m, k, aromatic, smarts);

//    cerr << " from atom " << prev << " to " << istart << " smarts is " << smarts << endl;

      prev = istart;
      istart = k;
      break;
    }
  }

  smiles << '1';
  smarts << "%99";

  const auto rsize = r.number_elements();

  int ndx = 0;
  if (6 == rsize && aromatic)
    ndx = RING_6a;
  else if (6 == rsize && ! aromatic)
    ndx = RING_6A;
  else if (5 == rsize && aromatic)
    ndx = RING_5a;
  else if (5 == rsize && ! aromatic)
    ndx = RING_5A;
  else if (7 == rsize && aromatic)
    ndx = RING_7a;
  else if (7 == rsize && ! aromatic)
    ndx = RING_7A;
  else if (7 == rsize && aromatic)
    ndx = RING_7a;
  else if (4 == rsize && ! aromatic)
    ndx = RING_4A;
  else if (8 == rsize && ! aromatic)
    ndx = RING_8A;
  else
  {
    cerr << "Unprocessed ring form " << rsize << " aromatic? " << aromatic << ' ' << m.smiles() << ' ' << m.name() << endl;
    return 0;
  }
  
  std::ostream & output = output_stream[ndx];

  output << m.unique_smiles() << ' ' << smiles << ' ' << smarts << ' ' << mname << suffix[ndx] << '\n';

  return 1;
}

static int
do_extract_all_isolated_rings (Molecule & m)
{
  const auto nr = m.nrings();

  const auto matoms = m.natoms();

  int * include_these_atoms = new int[matoms]; std::unique_ptr<int[]> free_include_these_atoms(include_these_atoms);

  for (auto i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_fused())
      continue;

    if (ri->number_elements() < min_ring_size)
      continue;
    else if (ri->number_elements() > max_ring_size)
      continue;

//  cerr << "Processing ring " << (*ri) << endl;

    set_vector(include_these_atoms, matoms, 0);
    ri->set_vector(include_these_atoms, 1);
    if (! rerc.identify_atoms_associated_with_ring_system(m, include_these_atoms))
      continue;

    Molecule mcopy;
    m.create_subset(mcopy, include_these_atoms);
//  cerr << "Subset is " << mcopy.smiles() << endl;
    if (mcopy.natoms() < matoms)
      write_ring(mcopy, m.name(), *ri, include_these_atoms);
  }

  return 1;
}

static int
extract_rings (Molecule & m,
              std::ostream & output)
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 1;

  m.compute_aromaticity_if_needed();

  if (extract_all_isolated_rings)
    do_extract_all_isolated_rings(m);

  const auto matoms = m.natoms();

  int * in_same_ring = new_int(matoms * matoms); std::unique_ptr<int[]> free_in_same_ring(in_same_ring);

  initialise_in_same_ring_array(m, in_same_ring);

  int * include_atom = new int[m.natoms()]; std::unique_ptr<int[]> free_include_atom(include_atom);

  if (extract_all_fused_systems)
    return extract_fused_systems(m, rerc, include_atom, output);

  int rings_this_molecule = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! rerc.can_be_processed(m, *ri))
      continue;

    if (0 == min_heteroatoms_needed_in_ring)
      ;
    else if (! contains_enough_heteroatoms(m, *ri, min_heteroatoms_needed_in_ring))
      continue;

    Molecule mcopy(m);
    mcopy.set_name(m.name());
    set_vector(include_atom, matoms, 0);
    rerc.identify_atoms_associated_with_ring(mcopy, *ri, in_same_ring, include_atom);
    write_ring(mcopy, ri->item(0), rerc, include_atom, ri->is_aromatic(), output);
    rings_this_molecule++;
  }

  rings_per_molecule[rings_this_molecule]++;
  if (rings_this_molecule)
    molecules_with_rings++;


  return output.good();
}

static int
extract_rings (data_source_and_type<Molecule> & input,
               std::ostream & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! extract_rings(*m, output))
      return 0;
  }

  return output.good();
}

static int
extract_rings (const char * fname, int input_type, std::ostream & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return extract_rings(input, output);
}

static int
open_output_stream (const IWString & prefix,
                    const char * name,
                    ofstream & output,
                    IWString & suffix)
{
  IWString tmp(prefix);
  tmp << '_' << name << ".smi";
  output.open(tmp.null_terminated_chars(), std::ios::out);
  if (! output.good())
  {
    cerr << "open_output_stream:cannot open '" << tmp << "'\n";
    return 0;
  }

  suffix = '.';
  suffix << name;

  return 1;
}

static int
extract_rings (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lr:a:f:F:cub:wh:S:N:VR:4");

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

  if (cl.option_present('V'))
  {
    check_smiles_formed = 1;

    if (verbose)
      cerr << "Will check the smiles formed for validity\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_ring_size) || max_ring_size < min_ring_size || max_ring_size > 8)
    {
      cerr << "The maximum ring size option (-R) must be a valid ring size - internal limitation to size 8, see Ian\n";
      usage(1);
    }

    if (verbose)
      cerr << "Max ring size " << max_ring_size << endl;
  }

  if (! rerc.initialise(cl, verbose))
  {
    cerr << "Cannot initialise remove and extract option(s)\n";
    return 4;
  }

  if (cl.option_present('N'))
  {
    IWString n = cl.string_value('N');

    if (! open_output_stream(n, "4A", output_stream[RING_4A], suffix[RING_4A]) ||
        ! open_output_stream(n, "5A", output_stream[RING_5A], suffix[RING_5A]) ||
        ! open_output_stream(n, "5a", output_stream[RING_5a], suffix[RING_5a]) ||
        ! open_output_stream(n, "6a", output_stream[RING_6a], suffix[RING_6a]) ||
        ! open_output_stream(n, "6A", output_stream[RING_6A], suffix[RING_6A]) ||
        ! open_output_stream(n, "4A4A", output_stream[FUSED_4A4A], suffix[FUSED_4A4A]) ||
        ! open_output_stream(n, "4A5A", output_stream[FUSED_4A5A], suffix[FUSED_4A5A]) ||
        ! open_output_stream(n, "4A6a", output_stream[FUSED_4A6a], suffix[FUSED_4A6a]) ||
        ! open_output_stream(n, "4A6A", output_stream[FUSED_4A6A], suffix[FUSED_4A6A]) ||
        ! open_output_stream(n, "4A5a", output_stream[FUSED_4A5a], suffix[FUSED_4A5a]) ||
        ! open_output_stream(n, "5a5a", output_stream[FUSED_5a5a], suffix[FUSED_5a5a]) ||
        ! open_output_stream(n, "5a5A", output_stream[FUSED_5a5A], suffix[FUSED_5a5A]) ||
        ! open_output_stream(n, "5A5A", output_stream[FUSED_5A5A], suffix[FUSED_5A5A]) ||

        ! open_output_stream(n, "5a6a", output_stream[FUSED_5a6a], suffix[FUSED_5a6a]) ||
        ! open_output_stream(n, "5a6A", output_stream[FUSED_5a6A], suffix[FUSED_5a6A]) ||
        ! open_output_stream(n, "5A6a", output_stream[FUSED_5A6a], suffix[FUSED_5A6a]) ||
        ! open_output_stream(n, "5A6A", output_stream[FUSED_5A6A], suffix[FUSED_5A6A]) ||

        ! open_output_stream(n, "6A6A", output_stream[FUSED_6A6A], suffix[FUSED_6A6A]) ||
        ! open_output_stream(n, "6a6a", output_stream[FUSED_6a6a], suffix[FUSED_6a6a]) ||
        ! open_output_stream(n, "6a6A", output_stream[FUSED_6a6A], suffix[FUSED_6a6A]))
    {
      cerr << "Cannot open one or more output files with stem '" << n << "'\n";
      return 1;
    }

    if (max_ring_size > 6)
    {
      if (! open_output_stream(n,   "7A", output_stream[RING_7A],    suffix[RING_7A]) ||
          ! open_output_stream(n,   "7a", output_stream[RING_7a],    suffix[RING_7a]) ||
          ! open_output_stream(n, "4A7A", output_stream[FUSED_4A7A], suffix[FUSED_4A7A]) ||
          ! open_output_stream(n, "4A7a", output_stream[FUSED_4A7a], suffix[FUSED_4A7a]) ||
          ! open_output_stream(n, "5A7A", output_stream[FUSED_5A7A], suffix[FUSED_5A7A]) ||
          ! open_output_stream(n, "5a7A", output_stream[FUSED_5a7A], suffix[FUSED_5a7A]) ||
          ! open_output_stream(n, "5A7a", output_stream[FUSED_5A7a], suffix[FUSED_5A7a]) ||
          ! open_output_stream(n, "5a7a", output_stream[FUSED_5a7a], suffix[FUSED_5a7a]) ||
          ! open_output_stream(n, "6a7a", output_stream[FUSED_6a7a], suffix[FUSED_6a7a]) ||
          ! open_output_stream(n, "6a7A", output_stream[FUSED_6a7A], suffix[FUSED_6a7A]) ||
          ! open_output_stream(n, "6A7A", output_stream[FUSED_6A7A], suffix[FUSED_6A7A]) ||
          ! open_output_stream(n, "7A7A", output_stream[FUSED_7A7A], suffix[FUSED_7A7A]) ||
          ! open_output_stream(n, "7a7a", output_stream[FUSED_7a7a], suffix[FUSED_7a7a]) ||
          ! open_output_stream(n, "7a7A", output_stream[FUSED_7a7A], suffix[FUSED_7a7A]) ||
          ! open_output_stream(n, "6A7a", output_stream[FUSED_6A7a], suffix[FUSED_6A7a]))
      {
        cerr << "Cannot open one or more output files stem '" << n << "' ring size 7\n";
        return 1;
      }

      if (max_ring_size > 7)
      {
        if (! open_output_stream(n,   "8A", output_stream[RING_8A],    suffix[RING_8A]) ||
            ! open_output_stream(n, "5A8A", output_stream[FUSED_5A8A], suffix[FUSED_5A8A]) ||
            ! open_output_stream(n, "5a8A", output_stream[FUSED_5a8A], suffix[FUSED_5a8A]) ||
            ! open_output_stream(n, "6A8A", output_stream[FUSED_6A8A], suffix[FUSED_6A8A]) ||
            ! open_output_stream(n, "6a8A", output_stream[FUSED_6a8A], suffix[FUSED_6a8A]) ||
            ! open_output_stream(n, "7A8A", output_stream[FUSED_7A8A], suffix[FUSED_7A8A]) ||
            ! open_output_stream(n, "7a8A", output_stream[FUSED_7a8A], suffix[FUSED_7a8A]) ||
            ! open_output_stream(n, "8A8A", output_stream[FUSED_8A8A], suffix[FUSED_8A8A]) ||
            ! open_output_stream(n, "4A8A", output_stream[FUSED_4A8A], suffix[FUSED_4A8A]))
        {
          cerr << "Cannot open one or more output files stem '" << n << "' ring size 8\n";
          return 1;
        }
      }
    }

    if (cl.option_present('4'))
    {
      if (! open_output_stream(n, "4a4a", output_stream[FUSED_4a4a], suffix[FUSED_4a4a]) ||
          ! open_output_stream(n, "4a6A", output_stream[FUSED_4a6A], suffix[FUSED_4a6A]) ||
          ! open_output_stream(n, "4a7a", output_stream[FUSED_4a7a], suffix[FUSED_4a7a]) ||
          ! open_output_stream(n, "4a4A", output_stream[FUSED_4a4A], suffix[FUSED_4a4A]) ||
          ! open_output_stream(n, "4a7A", output_stream[FUSED_4a7A], suffix[FUSED_4a7A]) ||
          ! open_output_stream(n, "4a5a", output_stream[FUSED_4a5a], suffix[FUSED_4a5a]))
      {
        cerr << "Cannot open one or more 4a output files stem '" << n << "'\n";
        return 1;
      }
    }

    extract_all_isolated_rings = 1;
    extract_all_fused_systems = 1;
  }

  if (cl.option_present('u'))
  {
    unique_rings_only = 1;
    if (verbose)
      cerr << "Will only write unique rings\n";
  }

  if (cl.option_present('w'))
  {
    must_match_whole_ring = 0;

    if (verbose)
      cerr << "No requirement for matching whole rings\n";
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', min_heteroatoms_needed_in_ring) || min_heteroatoms_needed_in_ring < 0)
    {
      cerr << "The min heteroatoms in ring (-h) option must have a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Only rings containing at least " << min_heteroatoms_needed_in_ring << " heteroatoms will be processed\n";
  }

  int input_type = 0;

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

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! extract_rings(cl[i], input_type, std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (extract_all_fused_systems)
  {
    for (auto i = 0; i < 19; ++i)
    {
      output_stream[i].close();
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    if (unique_rings_only)
      reaction_duplicate.report(cerr);
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = extract_rings(argc, argv);

  return rc;
}
