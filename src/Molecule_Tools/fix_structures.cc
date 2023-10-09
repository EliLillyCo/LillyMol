#include <array>
#include <iostream>
#include <memory>

/*
  There are various structure errors that we can fix.
  These aren't really chemical standardisations, but are changes
*/

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"
#include "fix_structures.h"
#include "Molecule_Lib/substructure.h"

using std::cerr;
using std::endl;

Structure_Fixing::Structure_Fixing()
{
  _active = 0;

  _directive_fix_almost_nitro = 0;
  _directive_fix_substituted_pyridine = 0;
  _directive_fix_four_valent_nitrogens_charged = 0;
  _directive_remove_hydrogens_known_flag_to_fix_valence_errors = 0;
  _directive_try_removing_charges_on_atoms_with_bad_valences = 0;
  _directive_add_positive_charge_to_three_valent_sulphur = 0;
  _directive_add_positive_charge_to_two_valent_halogen = 0;
  _directive_fix_obviously_wrong_implicit_hydrogens = 0;

  _molecules_examined = 0;
  _molecules_changed = 0;

  _partial_nitros_changed = 0;
  _substituted_pyridines_changed = 0;
  _four_valent_nitrogens_charged = 0;
  _three_valent_sulphurs_changed = 0;
  _two_valent_halogens_changed = 0;
  _charges_removed = 0;
  _charged_carbons_neutralised = 0;
  _charged_halogens_neutralised = 0;

  _directive_remove_charges_from_carbon_atoms = 0;
  _directive_remove_charges_from_halogen_atoms = 0;

  return;
}

int 
Structure_Fixing::report (std::ostream & os) const
{
  os << "Structure changing examined " << _molecules_examined << " molecules, changed " << _molecules_changed << endl;

  if (_partial_nitros_changed)
    os << " fixed " << _partial_nitros_changed << " partial nitro groups\n";

  if (_substituted_pyridines_changed)
    os << " fixed " << _substituted_pyridines_changed << " substituted pyridines\n";

  if (_four_valent_nitrogens_charged)
    os << " put charges on " << _four_valent_nitrogens_charged << " 4 valent neutral nitrogens\n";

  if (_charges_removed)
    os << " removed " << _charges_removed << " formal charges on atoms with bad valences\n";

  if (_three_valent_sulphurs_changed)
    os << " applied " << _three_valent_sulphurs_changed << " positive charges to 3 valence Sulphur atoms\n";

  if (_two_valent_halogens_changed)
    os << " applied " << _two_valent_halogens_changed << " positive charges to 2 valent Halogen atoms\n";

  if (_charged_carbons_neutralised)
    os << " removed " << _charged_carbons_neutralised << " charges from carbon atoms\n";

  if (_charged_halogens_neutralised)
    os << " removed " << _charged_halogens_neutralised << " charges from halogen atoms\n";

  return os.good();
}

int
Structure_Fixing::usage (char c, std::ostream & os) const
{
  os << "  -" << c << " nitro     changes -N(=O)-O to a Nitro\n";
  os << "  -" << c << " pyridine  puts a + charge on substited Pyridines\n";
  os << "  -" << c << " pquatn    puts a + charge neutral Quaternary Nitrogens\n";
  os << "  -" << c << " SD3+      puts a + charge on 3 connected, saturated Sulphur\n";
  os << "  -" << c << " XD2+      puts a + charge on 2 connected, saturated Halogens (Cl, Br, I)\n";
  os << "  -" << c << " app=xxx   append xxx to each changed structure\n";
  os << "  -" << c << " rmhsqb    remove implicit hydrogens known flag from otherwise OK atoms\n";
  os << "  -" << c << " rmbcharge try removing charges on atoms with bad valences\n";
  os << "  -" << c << " Xnc       remove formal charges from halogen atoms\n";
  os << "  -" << c << " Cnc       remove formal charges from carbon atoms\n";
  os << "  -" << c << " all       turn on all known changes\n";

  return os.good();
}

int
Structure_Fixing::initialise (Command_Line & cl, char c, int verbose)
{
  int i = 0;
  const_IWSubstring f;
  while (cl.value(c, f, i++))
  {
    if ("help" == f)
    {
      usage(c, cerr);
      exit(1);
    }
    else if (f.starts_with("app="))
    {
      f.remove_leading_chars(4);
      _append_to_changed_structures = f;

      if (verbose)
        cerr << "Will append '" << _append_to_changed_structures << "' to each molecule changed\n";
    }
    else if ("all" == f)
    {
      _directive_fix_almost_nitro = 1;
      _directive_fix_substituted_pyridine = 1;
      _directive_fix_four_valent_nitrogens_charged = 1;
      _directive_remove_hydrogens_known_flag_to_fix_valence_errors = 1;
      _directive_try_removing_charges_on_atoms_with_bad_valences = 1;
      _directive_add_positive_charge_to_three_valent_sulphur = 1;
      _directive_add_positive_charge_to_two_valent_halogen = 1;
      _directive_remove_charges_from_carbon_atoms = 1;
      _directive_remove_charges_from_halogen_atoms = 1;
      _directive_fix_obviously_wrong_implicit_hydrogens = 1;
      _active = 1;

      if (verbose)
        cerr << "All structure fixes turned on\n";
    }
    else if ("nitro" == f)
    {
      _directive_fix_almost_nitro = 1;

      if (verbose)
        cerr << "Will fix partial Nitro groups\n";

      _active = 1;
    }
    else if ("pyridine" == f)
    {
      _directive_fix_substituted_pyridine = 1;

      if (verbose)
        cerr << "Will fix substituted Pyridine groups\n";

      _active = 1;
    }
    else if ("pquatn" == f)
    {
      _directive_fix_four_valent_nitrogens_charged = 1;

      if (verbose)
        cerr << "Will apply a positive charge to neutral quaternary Nitrogen atoms\n";

      _active = 1;
    }
    else if ("rmhsqb" == f)
    {
      _directive_remove_hydrogens_known_flag_to_fix_valence_errors = 1;

      if (verbose)
        cerr << "Will allow Hydrogens to fill in bad valences if possible\n";

      _active = 1;
    }
    else if ("rmbcharge" == f)
    {
      _directive_try_removing_charges_on_atoms_with_bad_valences = 1;

      if (verbose)
        cerr << "Will try removing formal charges on atoms with bad valences\n";

      _active = 1;
    }
    else if (f == "fcih") {
      _directive_fix_obviously_wrong_implicit_hydrogens = 1;
      if (verbose)
        cerr << "Will fix obviously incorrect implicit hydrogen specifications\n";

      _active = 1;
    }
    else if ("SD3+" == f)
    {
      _directive_add_positive_charge_to_three_valent_sulphur = 1;

      if (verbose)
        cerr << "Will add positive charges to fully saturated, 3 connected Sulphur atoms\n";

      _active = 1;
    }
    else if ("XD2+" == f)
    {
      _directive_add_positive_charge_to_two_valent_halogen = 1;

      if (verbose)
        cerr << "Will add a positive charge to fully saturated, 2 connected halogen atoms\n";
    }
    else if ("Xnc" == f)
    {
      _directive_remove_charges_from_halogen_atoms = 1;
      
      if (verbose)
        cerr << "Will remove formal charges from halogen atoms\n";
    }
    else if ("Cnc" == f)
    {
      _directive_remove_charges_from_carbon_atoms = 1;
      
      if (verbose)
        cerr << "Will remove formal charges from carbon atoms\n";
    }
    else
    {
      cerr << "Unrecognised structure fixing qualifier '" << f << "'\n";
      usage(c, cerr);
      return 0;
    }
  }

  if (_directive_fix_substituted_pyridine)
  {
    if (! _pyridine_query.create_from_smarts("[ND3+0r6]1=C-C=C-C=C1"))
    {
      cerr << "Gack, could not initialise Pyridine smarts\n";
      return 0;
    }
  }

  return 1;
}

int
Structure_Fixing::_do_try_removing_charges_on_atoms_with_bad_valences (Molecule & m,
                     const Atom ** atom)
{
  const int matoms = m.natoms();

  int rc = 0;

  // Try placing these formal charges in order
  static std::array<formal_charge_t, 3> to_try{0, -1, 1};

  for (int i = 0; i < matoms; i++)
  {
    if (m.valence_ok(i))
      continue;

    if (0 == atom[i]->formal_charge())
      continue;

    const formal_charge_t initial_formal_charge = atom[i]->formal_charge();

    bool fixed = false;
    for (formal_charge_t charge : to_try) {
      if (charge == initial_formal_charge)
        continue;

      m.set_formal_charge(i, charge);

      if (m.valence_ok(i)) {
        fixed = true;
        break;
      }
    }

    if (fixed) {
      _charges_removed++;
      rc++;
    }
    else {
      m.set_formal_charge(i, initial_formal_charge);
    }
  }

  return rc;
}

int
Structure_Fixing::_do_fix_four_valent_nitrogens_charged (Molecule & m,
                     const Atom ** atom)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = atom[i];

    if (7 != ai->atomic_number())
      continue;

    if (0 != ai->formal_charge())
      continue;

    if (4 != ai->ncon())
      continue;

    if (4 != ai->nbonds())
      continue;

    m.set_formal_charge(i, 1);

    rc++;

    _four_valent_nitrogens_charged++;
  }

  return rc;
}

int
Structure_Fixing::_do_fix_almost_nitro (Molecule & m,
                     const Atom ** atom,
                     int & has_3_connected_4_bonds_uncharged_nitrogen)
{
  int rc = 0;

  has_3_connected_4_bonds_uncharged_nitrogen = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = atom[i];

    if (7 != ai->atomic_number())
      continue;

    if (ai->formal_charge())
      continue;

    if (3 != ai->ncon())
      continue;

    if (4 != ai->nbonds())
      continue;

    has_3_connected_4_bonds_uncharged_nitrogen++;

    if (m.is_ring_atom(i))    // Nitro's are never in a ring
      continue;

    atom_number_t doubly_bonded_oxygen = INVALID_ATOM_NUMBER;
    atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; j++)
    {
      atom_number_t k = ai->other(i, j);

      const Atom * ak = atom[k];

      if (8 != ak->atomic_number())
        continue;

      const Bond * b = ai->item(j);

      if (b->is_double_bond())     // there can only be one double bond
        doubly_bonded_oxygen = k;
      else if (INVALID_ATOM_NUMBER != singly_bonded_oxygen)    // gack, got two singly bonded Oxygens, what is this?
      {
        singly_bonded_oxygen = INVALID_ATOM_NUMBER;
        break;
      }
      else
        singly_bonded_oxygen = k;
    }

    if (INVALID_ATOM_NUMBER == singly_bonded_oxygen || INVALID_ATOM_NUMBER == doubly_bonded_oxygen)
      continue;

    m.set_bond_type_between_atoms(i, singly_bonded_oxygen, DOUBLE_BOND);

    rc++;

    _partial_nitros_changed++;
  }

  return rc;
}

int
Structure_Fixing::_do_fix_substituted_pyridine (Molecule & m,
                             const Atom ** atom)
{
  int rc = 0;

  Substructure_Results sresults;

  int nhits = _pyridine_query.substructure_search(m, sresults);

//cerr << "Checking pyridine, nhits = " << nhits << endl;

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    assert (6 == e->number_elements());

    atom_number_t n = e->item(0);

    const Atom * an = atom[n];

    if (4 != an->nbonds())
      continue;

    assert (7 == an->atomic_number() && 3 == an->ncon());
    assert (0 == an->formal_charge());

    m.set_formal_charge(n, 1);
    rc++;

    _substituted_pyridines_changed++;
  }

  return rc;
}

int
Structure_Fixing::_do_add_positive_charge_to_three_valent_sulphur (Molecule & m,
                             const Atom ** atom)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (3 != atom[i]->ncon())
      continue;

    if (0 != atom[i]->formal_charge())
      continue;

    if (16 != atom[i]->atomic_number())
      continue;

    if (atom[i]->nbonds() > atom[i]->ncon())
      continue;

    m.set_formal_charge(i, 1);
    rc++;

    _three_valent_sulphurs_changed++;

  }

  return rc;
}

int
Structure_Fixing::_do_add_positive_charge_to_two_valent_halogen (Molecule & m,
                             const Atom ** atom)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (2 != atom[i]->ncon())
      continue;

    if (0 != atom[i]->formal_charge())
      continue;

    if (atom[i]->nbonds() > atom[i]->ncon())
      continue;

    atomic_number_t z = atom[i]->atomic_number();
    if (17 == z || 35 == z || 53 == z)
    {
      m.set_formal_charge(i, 1);
      rc++;

      _two_valent_halogens_changed++;
    }
  }

  return rc;
}

int
Structure_Fixing::_do_remove_charges_from_halogen_atoms(Molecule & m, const Atom ** atom)
{
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == atom[i]->formal_charge())
      continue;

    if (! atom[i]->is_halogen())
      continue;

    if (1 == atom[i]->formal_charge() && 2 == atom[i]->ncon())
      continue;

    m.set_formal_charge(i, 0);
    rc++;

    _charged_halogens_neutralised++;
  }

  return rc;
}

int
Structure_Fixing::_do_remove_charges_from_carbon_atoms(Molecule & m, const Atom ** atom)
{
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == atom[i]->formal_charge())
      continue;

    if (6 != atom[i]->atomic_number())
      continue;

    m.set_formal_charge(i, 0);
    rc++;

    _charged_carbons_neutralised++;
  }

  return rc;
}

/*
  We return the number of fixes we make
*/

int
Structure_Fixing::_process(Molecule & m, const Atom ** atom)
{
  int rc = 0;

  int has_3_connected_4_bonds_uncharged_nitrogen = -1;

  if (_directive_fix_almost_nitro)
    rc += _do_fix_almost_nitro(m, atom, has_3_connected_4_bonds_uncharged_nitrogen);

  if (_directive_fix_substituted_pyridine)
  {
    if (0 == has_3_connected_4_bonds_uncharged_nitrogen)   // no need to look
      ;
    else
      rc += _do_fix_substituted_pyridine(m, atom);
  }

  if (_directive_fix_four_valent_nitrogens_charged)
    rc += _do_fix_four_valent_nitrogens_charged(m, atom);

  if (_directive_try_removing_charges_on_atoms_with_bad_valences && ! m.valence_ok())
    _do_try_removing_charges_on_atoms_with_bad_valences(m, atom);

  if (_directive_remove_hydrogens_known_flag_to_fix_valence_errors)
    rc += m.remove_hydrogens_known_flag_to_fix_valence_errors();

  if (_directive_add_positive_charge_to_three_valent_sulphur)
    rc += _do_add_positive_charge_to_three_valent_sulphur(m, atom);

  if (_directive_add_positive_charge_to_two_valent_halogen)
    rc += _do_add_positive_charge_to_two_valent_halogen(m, atom);

  if (_directive_remove_charges_from_halogen_atoms)
    rc += _do_remove_charges_from_halogen_atoms(m, atom);

  if (_directive_remove_charges_from_carbon_atoms)
    rc += _do_remove_charges_from_carbon_atoms(m, atom);

  if (rc)
  {
    _molecules_changed++;
    if (_append_to_changed_structures.length())
    {
      IWString tmp(m.name());
      tmp.append_with_spacer(_append_to_changed_structures);
      m.set_name(tmp);
    }
  }

  return rc;
}

int
Structure_Fixing::process(Molecule & m)
{
//cerr << "Structure_Fixing::process " << _directive_fix_obviously_wrong_implicit_hydrogens << '\n';
  _molecules_examined++;

  const int matoms = m.natoms();

  if (matoms < 4)
    return 0;

  const Atom ** atoms = new const Atom *[matoms]; std::unique_ptr<const Atom *[]> free_atoms(atoms);

  m.atoms(atoms);

  return _process(m, atoms);
}
