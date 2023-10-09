/*
  For the Sample ID project
  Try to fix 3 connected nitrogens with 4 bonds
*/

#include <stdlib.h>

#include "cmdline.h"

#include "molecule.h"
#include "path.h"
#include "output.h"
#include "aromatic.h"
#include "istream_and_type.h"

static int verbose = 0;

static int molecules_read = 0;

static int molecules_changed = 0;

static extending_resizable_array<int> number_instances;

static int halfnitro_count = 0;

static int beta_hydrogen_count = 0;

static int gamma_hydrogen_count = 0;

static int gamma_nh2ohsh_count = 0;

static int n_oxide_count = 0;

static int ohsh_ortho_count = 0;

static int amine_ortho_count = 0;

static int carbocycle_count = 0;

static Molecule_Output_Object stream_for_not_changed;

static int append_change_to_name = 0;

static int write_original_molecule = 0;

static IW_Regular_Expression rx;

static int molecules_not_processed_because_of_rx = 0;

static int write_unprocessed_molecules = 0;

static int formal_charge_placed_on_halogen = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -U <file>      stream for molecules not changed\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -E autocreate  automatically create new elements when encountered\n";
  cerr << "  -w             write the unchanged molecule to the output\n";
  cerr << "  -a             append change reason to name\n";
  cerr << "  -R <rx>        only process molecules whose name matches <rx>\n";
  cerr << "  -v             verbose output\n";

  exit (rc);
}

static int
put_positive_charge_on_conterion (Molecule & m)
{
  int nf = m.number_fragments ();

  if (1 == nf)
    return 0;

  int matoms = m.natoms ();
  for (int i = 0; i < matoms; i++)
  {
    if (0 != m.ncon (i))
      continue;

    if (0 != m.formal_charge (i))
      continue;

    if (! m.elementi (i)->is_halogen ())
      continue;

    formal_charge_placed_on_halogen++;
    cerr << molecules_read << " formal charge on halogen " << m.name () << endl;
    m.set_formal_charge (i, -1);
    return 1;
  }

  return 0;
}

static int
add_xminus (Molecule & m,
            int n)
{
  int matoms = m.natoms ();

  const Atom * a = m.atomi (0);

  coord_t xmin = a->x ();
  coord_t xmax = a->x ();
  coord_t ymax = a->y ();
  coord_t ymin = a->y ();

  for (int i = 1; i < matoms; i++)
  {
    const Atom * ai = m.atomi (i);
    if (ai->x () < xmin)
      xmin = ai->x ();
    else if (ai->x () > xmax)
      xmax = ai->x ();

    if (ai->y () > ymax)
      ymax = ai->y ();
    else if (ai->y () < ymin)
      ymin = ai-> y ();
  }

  Atom * xminus = new Atom ("X");
  xminus->set_formal_charge (-1);

  coord_t min_separation = 2.0;

  coord_t ymid = (ymax + ymin) * 0.5;
  coord_t x = xmin - (xmax - xmin) * 0.10;
  if (xmin - x < min_separation)
    x = xmin - min_separation;

  xminus->setxyz (x, ymid, static_cast<coord_t> (0.0));

  m.add (xminus);

  if (1 == n)
    return 1;

  Atom * x2 = new Atom ("X");
  x2->set_formal_charge (-1);

  x = xmax + (xmax - xmin) * 0.10;
  if (x - xmax < min_separation)
    x = xmax + min_separation;

  x2->setxyz (x, ymid, static_cast<coord_t> (0.0));

  m.add (x2);

  assert (2 == n);

  return 1;
}

static void
append_to_name (Molecule & m,
                const char * to_append)
{
  IWString mname = m.name ();

  mname.append_with_spacer (to_append);

  m.set_name (mname);

  return;
}

/*
  Change an NH2, OH or SH outside a ring to the double bond form and shift the
  intervening bonds.
*/

static int
nh2ohsh_four_bonds_away (Molecule & m,
                        atom_number_t n)
{
  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t alpha = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    alpha = b->other (n);
    break;
  }

  if (INVALID_ATOM_NUMBER == alpha)
    return 0;

  const Atom * aa = m.atomi (alpha);

  atom_number_t beta = INVALID_ATOM_NUMBER;

  ncon = aa->ncon ();

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = aa->item (i);

    if (! b->is_single_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (alpha);

    if (7 == m.atomic_number (x) && 2 == m.ncon (x))
      ;
    else if (6 == m.atomic_number (x) && 3 == m.ncon (x))
      ;
    else
      continue;

    beta = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == beta)
    return 0;

  const Atom * ab = m.atomi (beta);

  ncon = ab->ncon ();

  atom_number_t gamma = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = ab->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (beta);

    if (3 != m.ncon (x))
      continue;

    if (6 != m.atomic_number (x))
      continue;

    gamma = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == gamma)
    return 0;

  const Atom * ag = m.atomi (gamma);

  ncon = ag->ncon ();

  atom_number_t delta = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = ag->item (i);

    if (! b->is_single_bond ())
      continue;

    if (b->nrings ())
      continue;

    atom_number_t x = b->other (gamma);

    if (0 == m.hcount (x))
      continue;

    if (1 != m.ncon (x))
      continue;

    if (m.nbonds (x) > m.ncon (x))
      continue;

    if (7 == m.atomic_number (x))
      ;
    else if (8 == m.atomic_number (x))
      ;
    else if (16 == m.atomic_number (x))
      ;
    else
      continue;

    delta = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == delta)
    return 0;

  m.set_bond_type_between_atoms (n, alpha, SINGLE_BOND);
  m.set_bond_type_between_atoms (alpha, beta, DOUBLE_BOND);
  m.set_bond_type_between_atoms (beta, gamma, SINGLE_BOND);
  m.set_bond_type_between_atoms (gamma, delta, DOUBLE_BOND);

  m.compute_aromaticity ();

  if (m.is_aromatic (n) && m.is_aromatic (alpha) && m.is_aromatic (beta) && m.is_aromatic (gamma))
    ;
  else
  {
    cerr << "Warning '" << m.name () << " gamma NH2OHSH aromaticity not achieved\n";
  }

  gamma_nh2ohsh_count++;

  if (append_change_to_name)
    append_to_name (m, "Gamma NH2OHSH");

  return 1;
}

static int
is_carbocycle (const Molecule & m,
               const Ring & r,
               atom_number_t n)
{
  int ring_size = r.number_elements ();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];
    
    if (j == n)
      continue;

    if (6 != m.atomic_number (j))
      return 0;
  }

  return 1;
}

/*
  One of these Nitrogens in a carbocycle is charged
*/

static int
is_carbocycle (Molecule & m,
               atom_number_t n)
{
  if (! m.is_ring_atom (n))
    return 0;

  if (0 != m.attached_heteroatom_count (n))
    return 0;

  int nr = m.nrings ();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi (i);

    if (! ri->contains (n))
      continue;

    if (! is_carbocycle (m, *ri, n))
      continue;

    m.set_formal_charge (n, 1);

    carbocycle_count++;

    if (append_change_to_name)
      append_to_name (m, "Carbocycle");

    return 1;
  }

  return 0;
}

/*
  Change [NH2]-C=N to [NH]=C-N
*/

static int
is_nh2_ortho_nitrogen (Molecule & m,
                       atom_number_t n)
{
  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t alpha = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (n);

    if (6 != m.atomic_number (x))
      continue;

    if (3 != m.ncon (x))
      continue;

    alpha = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == alpha)
    return 0;

  const Atom * aa = m.atomi (alpha);

  ncon = aa->ncon ();

  atom_number_t amine = INVALID_ATOM_NUMBER;
   
  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = aa->item (i);

    if (! b->is_single_bond ())
      continue;

    if (b->nrings ())
      continue;

    atom_number_t x = b->other (alpha);

    if (7 != m.atomic_number (x))
      continue;

    if (1 != m.ncon (x))
      continue;

    amine = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == amine)
    return 0;

  m.set_bond_type_between_atoms (n, alpha, SINGLE_BOND);
  m.set_bond_type_between_atoms (alpha, amine, DOUBLE_BOND);

  m.compute_aromaticity ();

  if (m.is_aromatic (n) && m.is_aromatic (alpha))
    ;
  else
  {
    cerr << "Warning, Amine ortho aromaticity not achieved '" << m.name () << "'\n";
  }

  if (append_change_to_name)
    append_to_name (m, "Amine Ortho");

  amine_ortho_count++;

  return 1;
}

/*
  Change [OH,SH]-C=N to [O,S]=C-N
*/

static int
is_ohsh_ortho_nitrogen (Molecule & m,
                        atom_number_t n)
{
  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t alpha = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (n);

    if (6 != m.atomic_number (x))
      continue;

    if (3 != m.ncon (x))
      continue;

    alpha = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == alpha)
    return 0;

  const Atom * aa = m.atomi (alpha);

  ncon = aa->ncon ();

  atom_number_t oxygen = INVALID_ATOM_NUMBER;
   
  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = aa->item (i);

    if (! b->is_single_bond ())
      continue;

    if (b->nrings ())
      continue;

    atom_number_t x = b->other (alpha);

    if (8 == m.atomic_number (x))
      ;
    else if (16 == m.atomic_number (x))
      ;
    else
      continue;

    if (1 != m.ncon (x))
      continue;

    oxygen = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == oxygen)
    return 0;

  m.set_bond_type_between_atoms (n, alpha, SINGLE_BOND);
  m.set_bond_type_between_atoms (alpha, oxygen, DOUBLE_BOND);

  m.compute_aromaticity ();

  if (m.is_aromatic (n) && m.is_aromatic (alpha))
    ;
  else
  {
    cerr << "Warning, OH,SH ortho aromaticity not achieved '" << m.name () << "'\n";
  }

  if (append_change_to_name)
    append_to_name (m, "OH,SH Ortho");

  ohsh_ortho_count++;

  return 1;
}

static int
is_n_oxide (Molecule & m,
            atom_number_t n)
{
  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_single_bond ())
      continue;

    if (b->nrings ())
      continue;

    atom_number_t x = b->other (n);

    if (8 != m.atomic_number (x))
      continue;

    if (1 != m.ncon (x))
      continue;

    singly_bonded_oxygen = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == singly_bonded_oxygen)
    return 0;

  m.set_formal_charge (n, 1);
  m.set_formal_charge (singly_bonded_oxygen, -1);

  m.compute_aromaticity ();

  if (m.is_aromatic (n))
    ;
  else
  {
    cerr << "Warning '" << m.name () << "' aromaticity not achieved with N-Oxide\n";
  }

  if (append_change_to_name)
    append_to_name (m, "N-Oxide");

  n_oxide_count++;

  return 1;
}

static int
hydrogen_four_bonds_away (Molecule & m,
                        atom_number_t n)
{
  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t alpha = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    alpha = b->other (n);
    break;
  }

  if (INVALID_ATOM_NUMBER == alpha)
    return 0;

  const Atom * aa = m.atomi (alpha);

  atom_number_t beta = INVALID_ATOM_NUMBER;

  ncon = aa->ncon ();

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = aa->item (i);

    if (! b->is_single_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (alpha);

    if (! m.in_same_ring (n, x))
      continue;

    beta = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == beta)
    return 0;

  const Atom * ab = m.atomi (beta);

  ncon = ab->ncon ();

  atom_number_t gamma = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = ab->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (beta);

    if (! m.in_same_ring (n, x))
      continue;

    gamma = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == gamma)
    return 0;

  const Atom * ag = m.atomi (gamma);

  ncon = ag->ncon ();

  atom_number_t delta = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = ag->item (i);

    if (! b->is_single_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (gamma);

    if (! m.in_same_ring (n, x))
      continue;

    if (1 != m.hcount (x))
      continue;

    if (m.nbonds (x) > m.ncon (x))
      continue;

    delta = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == delta)
    return 0;

  m.set_bond_type_between_atoms (n, alpha, SINGLE_BOND);
  m.set_bond_type_between_atoms (alpha, beta, DOUBLE_BOND);
  m.set_bond_type_between_atoms (beta, gamma, SINGLE_BOND);
  m.set_bond_type_between_atoms (gamma, delta, DOUBLE_BOND);

  m.compute_aromaticity ();

  if (m.is_aromatic (n) && m.is_aromatic (alpha) && m.is_aromatic (beta) && m.is_aromatic (gamma) && m.is_aromatic (delta))
    ;
  else
  {
    cerr << "Warning '" << m.name () << " Hydrogen 4 bonds, aromaticity not achieved\n";
  }

  gamma_hydrogen_count++;

  if (append_change_to_name)
    append_to_name (m, "Gamma Hydrogen");

  return 1;
}

/*
  See if we can move the double bond
*/

static int
hydrogen_one_bond_away (Molecule & m,
                        atom_number_t n)
{
#ifdef DEBUG_HYDROGEN_ONE_BOND_AWAY
  cerr << "Testing " << n << " for hydrogen_one_bond_away\n";
#endif

  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  atom_number_t doubly_bonded_connection = INVALID_ATOM_NUMBER;

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    if (! b->is_double_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (n);

    if (! m.is_ring_atom (x))    // redundant
      continue;

//  If the other atom is another ND3v4, don't process it - 480354

    if (7 == m.atomic_number (x) && 4 == m.nbonds (x) && 0 == m.formal_charge (x) && 3 == m.ncon (x))
      continue;

    doubly_bonded_connection = x;
    break;
  }

#ifdef DEBUG_HYDROGEN_ONE_BOND_AWAY
  cerr << "Doubly bonded connection is " << doubly_bonded_connection << endl;
#endif

  if (INVALID_ATOM_NUMBER == doubly_bonded_connection)
    return 0;

// Now look for an atom singly bonded

  const Atom * alpha = m.atomi (doubly_bonded_connection);

  atom_number_t beta = INVALID_ATOM_NUMBER;

  ncon = alpha->ncon ();

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = alpha->item (i);

    if (! b->is_single_bond ())
      continue;

    if (0 == b->nrings ())
      continue;

    atom_number_t x = b->other (doubly_bonded_connection);

    if (1 != m.hcount (x))
      continue;

    if (m.ncon (x) > 3)
      continue;

    if (m.nbonds (x) > m.ncon (x))
      continue;

    if (! m.in_same_ring (n, x))
      continue;

    beta = x;
    break;
  }

  if (INVALID_ATOM_NUMBER == beta)
    return 0;

  m.set_bond_type_between_atoms (n, doubly_bonded_connection, SINGLE_BOND);
  m.set_bond_type_between_atoms (doubly_bonded_connection, beta, DOUBLE_BOND);

  m.compute_aromaticity ();

  if (m.is_aromatic (n) && m.is_aromatic (doubly_bonded_connection) && m.is_aromatic (beta))
    ;
  else
  {
    cerr << "Warning, '" << m.name () << "' aromaticity not obtained\n";
//  m.set_bond_type_between_atoms (n, doubly_bonded_connection, DOUBLE_BOND);
//  m.set_bond_type_between_atoms (doubly_bonded_connection, beta, SINGLE_BOND);

//  m.compute_aromaticity ();
  }

  beta_hydrogen_count++;

  if (append_change_to_name)
    append_to_name (m, "Beta Hydrogen");

  return 1;
}

static int
is_half_nitro (Molecule & m,
               atom_number_t n)
{
  if (m.is_ring_atom (n))
    return 0;

  atom_number_t doubly_bonded_oxygen = INVALID_ATOM_NUMBER;
  atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;

  const Atom * an = m.atomi (n);

  int ncon = an->ncon ();

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = an->item (i);

    atom_number_t o = b->other (n);

    if (8 != m.atomic_number (o))
      continue;

    if (b->is_double_bond () && INVALID_ATOM_NUMBER == doubly_bonded_oxygen)
      doubly_bonded_oxygen = o;
    else if (b->is_single_bond () && INVALID_ATOM_NUMBER == singly_bonded_oxygen)
      singly_bonded_oxygen = o;
  }

  if (INVALID_ATOM_NUMBER == singly_bonded_oxygen || INVALID_ATOM_NUMBER == doubly_bonded_oxygen)
    return 0;

  m.set_formal_charge (n, 1);
  m.set_formal_charge (singly_bonded_oxygen, -1);

  halfnitro_count++;

  if (append_change_to_name)
    append_to_name (m, "Half Nitro");

  return 1;
}

static int
fix_nd3v4 (Molecule & m,
           atom_number_t n,
           int & changed,
           int & need_to_add_xminus)
{
  changed = 0;

  if (is_half_nitro (m, n))
  {
    changed = 1;
    return 1;
  }

  if (is_n_oxide (m, n))
  {
    changed = 1;
    return 1;
  }

  if (hydrogen_one_bond_away (m, n))
  {
    changed = 1;
    return 1;
  }

  if (hydrogen_four_bonds_away (m, n))
  {
    changed = 1;
    return 1;
  }

  if (nh2ohsh_four_bonds_away (m, n))
  {
    changed = 1;
    return 1;
  }

  if (is_ohsh_ortho_nitrogen (m, n))
  {
    changed = 1;
    return 1;
  }

  if (is_nh2_ortho_nitrogen (m, n))
  {
    changed = 1;
    return 1;
  }

  if (is_carbocycle (m, n))
  {
    if (put_positive_charge_on_conterion (m))
      ;
    else
      need_to_add_xminus++;
    changed = 1;
    return 1;
  }

  changed = 0;

  if (stream_for_not_changed.active ())
    stream_for_not_changed.write (m);

  return 1;
}

static int
fix_nd3v4 (Molecule & m,
           Molecule_Output_Object & output)
{
  IWString initial_molecular_formula;
  m.molecular_formula (initial_molecular_formula);

  int initial_aromatic_atom_count = m.aromatic_atom_count ();

  int matoms = m.natoms ();

  int cases_this_molecule = 0;
  int changes_this_molecule = 0;

  Molecule mcopy;
  if (write_original_molecule)
    mcopy.add_molecule (&m);

  int need_to_add_xminus = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi (i);

    if (7 != ai->atomic_number ())
      continue;

    if (3 != ai->ncon ())
      continue;

    if (4 != ai->nbonds ())
      continue;

    if (ai->formal_charge ())
      continue;

    cases_this_molecule++;

    int changed;
    if (fix_nd3v4 (m, i, changed, need_to_add_xminus) && changed)
    {
      changes_this_molecule++;
    }
  }

  if (changes_this_molecule)
  {
    if (write_original_molecule)
    {
      mcopy.set_name (m.name ());
      output.write (mcopy);
    }

    if (need_to_add_xminus)
      add_xminus (m, need_to_add_xminus);

    output.write (m);
    molecules_changed++;
  }
  else if (write_unprocessed_molecules)
  {
    output.write (m);
  }

  number_instances[cases_this_molecule]++;

  return 1;
}

static void
preprocess (Molecule & m)
{
  return;
}

static int
fix_nd3v4 (data_source_and_type<Molecule> & input,
                               Molecule_Output_Object & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule ()))
  {
    molecules_read++;

    if (rx.active () && ! rx.matches (m->name ()))
    {
      if (write_unprocessed_molecules)
        output.write (m);

      molecules_not_processed_because_of_rx++;

      delete m;
      continue;
    }

    if (verbose > 1)
      cerr << "Processing '" << m->molecule_name () << "'\n";

    preprocess (*m);

    if (! fix_nd3v4 (*m, output))
    {
      delete m;
      return 0;
    }

    delete m;
  }

  return 1;
}

static int
fix_nd3v4 (const char * fname, 
                               int input_type,
                               Molecule_Output_Object & output)
{
  if (0 == input_type)
  {
    input_type = discern_file_type_from_name (fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return fix_nd3v4 (input, output);
}

static int
fix_nd3v4 (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:U:i:o:S:awR:z");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! process_elements (cl, verbose > 1, 'E'))
  {
    cerr << "Cannot process -E option\n";
    usage (7);
  }

  set_auto_create_new_elements (1);

  if (! process_standard_aromaticity_options (cl, verbose > 1))
  {
    cerr << "Cannot process -A option\n";
    usage (11);
  }

  if (cl.option_present ('a'))
  {
    append_change_to_name = 1;

    if (verbose)
      cerr << "The reason for the change will be appended to the name\n";
  }

  if (cl.option_present ('w'))
  {
    write_original_molecule = 1;

    if (verbose)
      cerr << "Will write the unchanged molecule\n";
  }

  if (cl.option_present ('R'))
  {
    const IWString & r = cl.string_value ('R');

    if (! rx.set_pattern (r))
    {
      cerr << "Invalid regular expression specifier '" << r << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Will only process molecules whose name matches '" << rx.source () << "'\n";
  }

  if (cl.option_present ('z'))
  {
    if (! rx.active ())
    {
      cerr << "The write unprocessed molecules (-z) option only makes sense with a regular expression\n";
      usage (16);
    }

    write_unprocessed_molecules = 1;

    if (verbose)
      cerr << "Molecules that don't hit the regular expression will be written\n";
  }

  int input_type = 0;
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

  if (cl.option_present ('U'))
  {
    if (! cl.option_present ('o'))
      stream_for_not_changed.add_output_type (SMI);
    else if (! stream_for_not_changed.determine_output_types (cl))
    {
      cerr << "Cannot determine output type(s)\n";
      usage (18);
    }

    const_IWSubstring u = cl.string_value ('U');

    if (stream_for_not_changed.would_overwrite_input_files (cl, u))
    {
      cerr << "The -U option cannot overwrite its input file(s)\n";
      return 5;
    }

    if (! stream_for_not_changed.new_stem (u))
    {
      cerr << "Cannot initialise stem for unchanged molecules '" << u << "'\n";
      return 17;
    }
  }

  Molecule_Output_Object output;

  if (! cl.option_present ('o'))
    output.add_output_type (SMI);
  else if (! output.determine_output_types (cl))
  {
    cerr << "Cannot determine output type(s)\n";
    usage (18);
  }

  if (cl.option_present ('S'))
  {
    const_IWSubstring s = cl.string_value ('S');

    if (! output.new_stem (s))
    {
      cerr << "Cannot open output stream '" << s << "'\n";
      return 3;
    }
  }
  else if (! output.new_stem ("-"))
  {
    cerr << "Huh, cannot initialise stdout\n";
    return 2;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! fix_nd3v4 (cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Processed " << molecules_read << " molecules, changed " << molecules_changed << endl;
    if (rx.active ())
      cerr << molecules_not_processed_because_of_rx << " molecules not processed because of regular expression\n";

    for (int i = 0; i < number_instances.number_elements (); i++)
    {
      if (number_instances[i])
        cerr << number_instances[i] << " molecules had " << i << " instances of the group\n";
    }

    cerr << halfnitro_count << " half nitro groups\n";
    cerr << beta_hydrogen_count << " beta Hydrogen changes\n";
    cerr << gamma_hydrogen_count << " gamma Hydrogen changes\n";
    cerr << gamma_nh2ohsh_count << " gamma NH2OHSH changes\n";
    cerr << n_oxide_count << " N-Oxide changes\n";
    cerr << ohsh_ortho_count << " OH,SH Ortho changes\n";
    cerr << amine_ortho_count << " Amine Ortho changes\n";
    cerr << carbocycle_count << " carbocycle changes\n";
    cerr << formal_charge_placed_on_halogen << " formal charges placed on halogen atoms\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = fix_nd3v4 (argc, argv);

  return rc;
}
