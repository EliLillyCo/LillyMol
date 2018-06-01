#include <stdlib.h>

/*
  Define symbol so we get private molecule functions
*/

#define COMPILING_SMILES_CC
#define COMPILING_CTB

#include "cmdline.h"
#include "misc.h"

#include "molecule.h"
#include "misc2.h"
#include "smiles.h"
#include "chiral_centre.h"
#include "pearlman.h"

/*
  Some smiles parsers are picky about reusing numbers for ring closures.
  We can do it either way....
*/

static int reuse_ring_closure_numbers = 1;

void
set_smiles_reuse_ring_closure_numbers (int i)
{
  assert (i >= 0);

  reuse_ring_closure_numbers = i;

  return;
}

int
smiles_reuse_ring_closure_numbers()
{
  return reuse_ring_closure_numbers;
}


static int include_aromaticity_in_smiles = 0;

void
set_include_aromaticity_in_smiles (int i)
{
  if (i == include_aromaticity_in_smiles)
    return;

  include_aromaticity_in_smiles = i;

  return;
}

int
get_include_aromaticity_in_smiles()
{
  return include_aromaticity_in_smiles;
}

static int _include_cis_trans_in_smiles = 1;

void
set_include_cis_trans_in_smiles (int i)
{
  _include_cis_trans_in_smiles = i;

  return;
}

int
include_cis_trans_in_smiles()
{
  return _include_cis_trans_in_smiles;
}

static int _include_chiral_info_in_smiles = 1;

void
set_include_chiral_info_in_smiles (int i)
{
  _include_chiral_info_in_smiles = i;
}

int
include_chiral_info_in_smiles()
{
  return _include_chiral_info_in_smiles;
}

static int _ignore_chiral_info_on_input = 0;

void
set_ignore_chiral_info_on_input (int i)
{
  _ignore_chiral_info_on_input = i;

  return;
}

int
ignore_chiral_info_on_input()
{
  return _ignore_chiral_info_on_input;
}

static int append_coordinates_after_each_atom = 0;

void
set_append_coordinates_after_each_atom (int s)
{
  append_coordinates_after_each_atom = s;
}

static int _write_smiles_with_smarts_atoms = 0;

void
set_write_smiles_with_smarts_atoms (int s)
{
  _write_smiles_with_smarts_atoms = s;

  return;
}

int
write_smiles_with_smarts_atoms()
{
  return _write_smiles_with_smarts_atoms;
}

static int include_implicit_hydrogens_on_aromatic_n_and_p = 1;

void
set_include_implicit_hydrogens_on_aromatic_n_and_p (int s)
{
  include_implicit_hydrogens_on_aromatic_n_and_p = s;
}

static int _add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens = 0;

void
set_add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens (int s)
{
  _add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens = s;
}

int
add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens()
{
  return _add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens;
}

static int _write_single_bonds_in_smiles = 0;

int
write_single_bonds_in_smiles()
{
  return _write_single_bonds_in_smiles;
}

void
set_write_single_bonds_in_smiles (int s)
{
  _write_single_bonds_in_smiles = s;
}

static int include_hcount_in_smiles = 1;

void
set_include_hcount_in_smiles(int s)
{
  include_hcount_in_smiles = s;
}

static int file_scope_display_unusual_hcount_warning_messages = 1;

void
set_display_unusual_hcount_warning_messages(int s)
{
  file_scope_display_unusual_hcount_warning_messages = s;
}

int
display_unusual_hcount_warning_messages()
{
  return file_scope_display_unusual_hcount_warning_messages;
}

static int _include_atom_map_with_smiles = 1;

void
set_include_atom_map_with_smiles(const int s)
{
  _include_atom_map_with_smiles = s;
}

int
include_atom_map_with_smiles()
{
  return _include_atom_map_with_smiles;
}

/*
  Sept 2003. Ran into ring openings and closings from ChemDraw that
  went over 100!
  Jul 2014. These molecules generate all kinds of problems.
  How should %101 be parsed? Is it ring 10 followed by ring 1, or
  is it ring 101?

  Aug 2014. Implement the idea of %(nnn)
*/

/*
  Options for smiles are specified with the -SFLAG switch.
  We set it as a variable to make it easy to dynamically change it if
  ever that will be needed ???
*/

static char sflag = 'K';

int
display_standard_smiles_options (std::ostream & os)
{
  os << "  -" << sflag << " <...>       Enter \"-" << sflag << " help\" for SMILES options\n";

  return os.good();
}

int
display_all_smiles_options (char flag, std::ostream & os)
{
  assert (os.good());

  os << "  -" << flag << " nru         do not reuse ring closure numbers when outputting smiles\n";
  os << "  -" << flag << " ctb         include cis-trans bonds in the smiles\n";
  os << "  -" << flag << " <number>    use <number> as a seed for random smiles\n";
  os << "  -" << flag << " nochiral    exclude chiral info from all smiles produced\n";
  os << "  -" << flag << " coords      include coordinates in the smiles - '{{xx,yy,zz}}'\n";
  os << "  -" << flag << " Hiso        provide implicit Hydrogens to isotopic atoms needing H atoms\n";
  os << "  -" << flag << " nhis        exclude implicit Hydrogen information from smiles\n";
  os << "  -" << flag << " iusmi       include isotopic information when determining canonical order\n";
  os << "  -" << flag << " niusmi      exclude isotopic information when determining canonical order\n";
  os << "  -" << flag << " dbusmi      include directional bonding information when determining canonical order\n";
  os << "  -" << flag << " ndbusmi     exclude directional bonding information when determining canonical order\n";
  os << "  -" << flag << " rnoff=nn    add NN to all ring numbers\n";
  os << "  -" << flag << " nonH        do NOT write implicit Hydrogen info on pyrrole-like N atoms\n";
  os << "  -" << flag << " nihc        do NOT consider implicit Hydrogens during canonicalisation\n";
  os << "  -" << flag << " esssr       perceive the extended SSSR ring set\n";
  os << "  -" << flag << " noD         do NOT include the D operator in smarts generated\n";
  os << "  -" << flag << " iso01       when computing unique smiles, isotopes are zero or non zero\n";
  os << "  -" << flag << " rcsbd       include directionality in ring closure single bonds\n";
  os << "  -" << flag << " usv2        use faster, but incompatible unique smiles determination\n";
  os << "  -" << flag << " nd4h        four connected neutral Nitrogen atoms have a Hydrogen\n";

  return 1;
}

int
process_standard_smiles_options (Command_Line & cl, int verbose,
                                 const char sflag)
{
  int i = 0;
  IWString tmp;
  while (cl.value(sflag, tmp, i))
  {
    long seed;

    if ("nru" == tmp)
    {
      set_smiles_reuse_ring_closure_numbers(0);
      if (verbose)
        cerr << "Smiles will not re-use ring closure numbers\n";
    }
    else if ("ctb" == tmp)
    {
      set_include_cis_trans_in_smiles(1);
      if (verbose)
        cerr << "Cis trans bonds will be included in the smiles\n";
    }
    else if ("xtb" == tmp)
    {
      set_include_cis_trans_in_smiles(0);
      if (verbose)
        cerr << "Cis trans bonds will be excluded from smiles\n";
    }
    else if ("nochiral" == tmp)
    {
      set_include_chiral_info_in_smiles(0);
      if (verbose)
        cerr << "Will exclude chiral info from any smiles produced\n";
    }
    else if ("coords" == tmp)
    {
      set_append_coordinates_after_each_atom(1);
      if (verbose)
        cerr << "Will add coordinates to each atom in a smiles\n";
    }
    else if ("random" == tmp)
    {
      random_number_seed_t tmp = set_smiles_random_number_seed_random();
      if (verbose)
        cerr << "Smiles generated with a random seed " << tmp << endl;
    }
    else if ("iusmi" == tmp)
    {
       set_include_isotopic_information_in_unique_smiles(1);
    }
    else if ("niusmi" == tmp)
    {
      set_include_isotopic_information_in_unique_smiles(0);
    }
    else if ("dbusmi" == tmp)
    {
      set_include_directional_bonding_information_in_unique_smiles(1);
    }
    else if ("ndbusmi" == tmp)
    {
      set_include_directional_bonding_information_in_unique_smiles(0);
    }
    else if (tmp.starts_with("rnoff="))
    {
      tmp.remove_leading_chars(6);
      int o;
      if (! tmp.numeric_value(o) || o < 0)
      {
        cerr << "Invalid ring number offset '" << tmp << "'\n";
        return 0;
      }

      set_smiles_ring_number_offset(o);
    }
    else if ("nonH" == tmp)
    {
      include_implicit_hydrogens_on_aromatic_n_and_p = 0;
    }
    else if ("nihc" == tmp)
    {
      set_consider_implicit_hydrogens_in_unique_smiles(0);
    }
    else if ("esssr" == tmp)
    {
      set_perceive_sssr_rings(0);
    }
    else if ("noD" == tmp)
    {
      set_include_D_in_smarts(0);
    }
    else if ("Hiso" == tmp)
    {
      set_add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens(1);
    }
    else if ("nhis" == tmp)
    {
      set_include_hcount_in_smiles(0);
    }
    else if ("iso01" == tmp)
    {
      set_consider_isotopes_as_zero_and_non_zero(1);
    }
    else if ("rcsbd" == tmp)
    {
      set_include_directionality_in_ring_closure_bonds(1);
    }
    else if ("usv2" == tmp)
    {
      set_unique_determination_version(2);
    }
    else if ("nd4h" == tmp)
    {
      set_four_connected_neutral_nitrogen_has_h(1);
    }
    else if ("help" == tmp)
    {
      display_all_smiles_options(sflag, cerr);
      exit(5);    // note very different behaviour for help!
    }
    else if (cl.value(sflag, seed, i))
    {
      set_smiles_random_number_seed(random_number_seed_t(seed));
      if (verbose)
        cerr << "Random smiles generated with seed " << seed << endl;
    }
    else
    {
      cerr << "Unrecognised smiles specifier '" << tmp << "'\n";
      return 0;
    }

    i++;
  }

  return 1;
}

/*int
Molecule::_smiles_write_directional_bond (atom_number_t a,
                     atom_number_t anchor,
                     IWString & smiles)
{

  const Cis_Trans_Bond * ctb = part_of_cis_trans_bond (a, anchor);
  if (NULL == ctb)
  {
    cerr << "Molecule::process_atom_for_smiles: directional bond, but no CTB\n";
    cerr << "Atoms " << a << " and " << anchor << endl;
    debug_print (cerr);
    assert (NULL == "This should not happen");
  }

  if (anchor == ctb->left_up() && a == ctb->left_root())
    smiles += '\\';
  else if (anchor == ctb->left_down() && a == ctb->left_root())
    smiles += '/';
  else if (anchor == ctb->right_up() && a == ctb->right_root())
    smiles += '\\';
  else if (anchor == ctb->right_down() && a == ctb->right_root())
    smiles += '/';
  else
  {
    cerr << "Molecule::_smiles_write_directional_bond: anchor = " << anchor << " a = " << a << endl;
    ctb->debug_print(cerr);
    assert (NULL == "What's going on here");
  }

  return 1;
}*/

/*
  Append the appropriate bond type '/' or '\' to the smiles.
  We are passing through B on our way to NEXT_ATOM
*/

int
Molecule::_process_directional_bond_for_smiles (IWString & smiles,
                                                const Bond * b,
                                                atom_number_t next_atom)
{
#ifdef DEBUG_PROCESS_DIRECTIONAL_BOND_FOR_SMILES
  cerr << "Directional bond from " << b->a1() << " to " << b->a2() << endl;
  cerr << " up " << b->is_directional_up() << " down " << b->is_directional_down() << endl;
#endif

  if (b->is_directional_up() && next_atom == b->a2())
    smiles += '/';
  else if (b->is_directional_up() && next_atom == b->a1())
    smiles += '\\';
  else if (b->is_directional_down() && next_atom == b->a2())
    smiles += '\\';
  else if (b->is_directional_down() && next_atom == b->a1())
    smiles += '/';
  else
  {
    cerr << "Molecule::_process_directional_bond_for_smiles: bad directional bond to " << next_atom << endl;
    cerr << (*b) << endl;
    return 0;
  }

#ifdef DEBUG_PROCESS_DIRECTIONAL_BOND_FOR_SMILES
  cerr << "Smiles now '" << smiles << "'\n";
#endif

  return 1;
}

static void
do_append_coordinates (IWString & smiles, const Atom * a)
{
  smiles += "{{";
  smiles.append_number(a->x(), 5);
  smiles += ',';
  smiles.append_number(a->y(), 5);
  smiles += ',';
  smiles.append_number(a->z(), 5);
  smiles += "}}";

  return;
}

static void
finish_smiles_atom (IWString & smiles,
                    atomic_number_t z,
                    int hcount,
                    formal_charge_t fc)
{
//cerr << "finish_smiles_atom atomic number " << z << " with hcount " << hcount << endl;

  if (hcount > 0 && 1 != z)
  {
     smiles += 'H';
     if (hcount > 1)
       append_digit(smiles, hcount);
  }
  if (fc > 0)
  {
    for (int i = 0; i < fc; i++)
    {
      smiles += '+';
    }
  }
  else if (fc < 0)
  {
    for (int i = 0; i < -fc; i++)
    {
      smiles += '-';
    }
  }

  return;
}

/*
  this is unfortunately complex.  Ideally the information about
  embedding or not would be carried along in the arguments, but that
  would make the inferface excessively complex. So, we persevere
  with this rather awful manipulation of the global state
*/

int
Molecule::_append_smarts_equivalent (Smiles_Formation_Info & sfi,
                                     IWString & s)
{
  atom_number_t zatom = sfi.zatom();

  const IWString * const user_specified_atomic_smarts = sfi.user_specified_atomic_smarts();

  if (NULL == user_specified_atomic_smarts)   // ignore
    ;
  else if (0 != user_specified_atomic_smarts[zatom].length())
  {
    s << user_specified_atomic_smarts[zatom];
    return 1;
  }

  int restore_make_smarts_embedding = make_smarts_embedding();   // make sure we restore the default value when done

  set_make_smarts_embedding(sfi.make_smarts_embedding(zatom));

//append_smarts_equivalent_for_atom(zatom, s, sfi);

  if (NULL != sfi.include_atom())
    append_smarts_equivalent_for_atom(zatom, s, sfi.include_atom());
  else
    append_smarts_equivalent_for_atom(zatom, s);

  set_make_smarts_embedding(restore_make_smarts_embedding);

  return 1;
}

static int
append_permanent_aromatic (IWString & smiles,
                           const Atom * a)
{
  const IWString & s = a->element()->aromatic_symbol();

  int needs_square_brackets;
  if (a->isotope() > 0 || a->element()->needs_square_brackets())
  {
    needs_square_brackets = 1;
    smiles << '[';
  }
  else
   needs_square_brackets = 0;

  if (a->isotope() > 0)
    smiles << a->isotope();

  if (s.length())
    smiles << s;
  else
    smiles << a->element()->symbol();

  if (needs_square_brackets)
    smiles << ']';

  if (append_coordinates_after_each_atom)
    do_append_coordinates(smiles, a);

  return 1;
}

//#define DEBUG_PROCESS_ATOM_FOR_SMILES

/*
  This is really a const function, but the signature of in_same_aromatic_ring
  is non-const.

  When dealing with chiral atoms which have attached H atoms, we must
  make sure that the hcount be added to the atom primitive

  If a chiral atom has ring openings, then we need to know which
  atoms will be ring openings. This is array RING_OPENINGS. Similarly
  ring closures
*/

int
Molecule::_process_atom_for_smiles (Smiles_Formation_Info & sfi,
                     const int * zorder,
                     const resizable_array<const Bond *> & ring_opening_bonds,
                     const resizable_array<atom_number_t> & ring_closures,
                     const Chiral_Centre * c,
                     IWString & smiles)
{
  if (! sfi.write_smiles())
    return _append_smarts_equivalent(sfi, smiles);

  atom_number_t zatom = sfi.zatom();

  const IWString * const user_specified_atomic_smarts = sfi.user_specified_atomic_smarts();

  if (NULL == user_specified_atomic_smarts)   // ignore it
    ;
  else if (user_specified_atomic_smarts[zatom].length() > 0)
  {
    smiles << user_specified_atomic_smarts[zatom];
    return 1;
  }

  Atom * a = _things[zatom];

  const Element * e = a->element();

  if (a->permanent_aromatic())
    return append_permanent_aromatic(smiles, a);

  atom_number_t anchor = sfi.previous_atom();

#ifdef DEBUG_PROCESS_ATOM_FOR_SMILES
  cerr << "Processing atom " << zatom << "(" << a->atomic_symbol() << ") ncon " << a->ncon() << " nbonds " << nbonds(zatom);
  if (a->formal_charge())
    cerr << " charge " << a->formal_charge();
  cerr << " IH " << a->implicit_hydrogens();
  if (a->implicit_hydrogens_known())
    cerr << " ihknown";
  if (c)
    cerr << " chiral";
  cerr << " anchor is " << anchor << endl;
  a->debug_print(cerr);
#endif

//assert (ok_atom_number (zatom));

  int hcount = -1;     // initialised to a negative number

// Look for unusual H counts. Remember, we must specify the H count for any charged atoms

  const formal_charge_t fc = a->formal_charge();

  const int stored_hcount = include_hcount_in_smiles ? a->implicit_hydrogens() : 0;

#ifdef DEBUG_PROCESS_ATOM_FOR_SMILES
  cerr << "Atom type " << a->atomic_symbol() << " stored_hcount " << stored_hcount << " include_hcount_in_smiles " << include_hcount_in_smiles << endl;
  cerr << "has " << a->implicit_hydrogens() << " IH, known? " << a->implicit_hydrogens_known() << endl;
#endif

  if (! include_hcount_in_smiles)
    ;
  else if (0 == fc)
  {
    if (stored_hcount == _compute_implicit_hydrogens(zatom))    // the most common case - normal
      ;
    else if (a->implicit_hydrogens_known())     // probably a radical
      ;
    else
    {
      if (e->organic() && file_scope_display_unusual_hcount_warning_messages)
        cerr << "Unusual hcount, atom " << zatom << "(" << a->atomic_symbol() << ") ncon = " << a->ncon() << " nbonds " << a->nbonds() << " stored = " << stored_hcount << " implicit is " << _compute_implicit_hydrogens(zatom) << endl;
      hcount = stored_hcount;
    }
  }
  else if (stored_hcount)
    hcount = stored_hcount;

// Daylight needs to include non zero hcount with aromatic nitrogens and phosphorus.
// Actually, Daylight seems never to aromatise phosphorus...

//cerr << __LINE__ << " hcount " << hcount << endl;

  if (hcount >= 0)            // already set
    ;
  else if (! include_aromaticity_in_smiles)     // nothing to worry about
    ;
  else if (! include_implicit_hydrogens_on_aromatic_n_and_p)   // don't worry about anything
    ;
  else
  {
    const atomic_number_t z = e->atomic_number();

    if (7 == z || 15 == z)
    {
      aromaticity_type_t arom;
      if (aromaticity(zatom, arom) && IS_AROMATIC_ATOM(arom))
        hcount = stored_hcount;

#ifdef DEBUG_PROCESS_ATOM_FOR_SMILES
      cerr << "Aromatic N or P has " << hcount << " implicit hydrogens\n";
#endif
    }
  }

// Real problems if we are excluding chirality from the smiles. Argument C will be NULL if we are
// excluding chirality from smiles, so we need to fetch the value again

  int ihknown = a->implicit_hydrogens_known();
  if (! _include_chiral_info_in_smiles && ihknown && (NULL != chiral_centre_at_atom(zatom)) && hcount < 0 && 0 == a->formal_charge() && 0 == a->isotope() && e->organic())
    ihknown = 0;

  int need_to_close_square_bracket = 0;

// Whenever a square bracket atom is encountered, the H count must always be explicit. 

#ifdef DEBUG_PROCESS_ATOM_FOR_SMILES
  cerr << "Atom " << zatom << " z = " << a->atomic_number() << " ncon " << a->ncon() << " c is " << c << " hcount " << hcount << ", ihknown " << ihknown << " stored_hcount " << stored_hcount << endl;
  cerr << "Line " << __LINE__ << ", length " << smiles.length() << endl;
#endif

  if (c || 
      hcount >= 0 ||
      0 != fc ||
      e->needs_square_brackets() ||
      a->isotope() ||
      ihknown ||
      (a->atom_map() > 0 && _include_atom_map_with_smiles))       // implicit_hydrogens_known() - possibly adjusted above
  {
    smiles += '[';
//  cerr << "Appended square bracket '" << smiles << "' hcount " << hcount << " ihknown " << ihknown << endl;
    need_to_close_square_bracket = 1;
    if (hcount < 0)
      hcount = stored_hcount;
  }

  aromaticity_type_t arom = NOT_AROMATIC;
  if (include_aromaticity_in_smiles)
    (void) aromaticity(zatom, arom);

  int rc = e->append_smiles_symbol(smiles, arom, a->isotope());

  if (c)
    c->append_smiles_chirality_symbol(smiles, zorder, anchor, 
                                       ring_opening_bonds, ring_closures);

  finish_smiles_atom(smiles, a->atomic_number(), hcount, fc);

  if (_include_atom_map_with_smiles && a->atom_map() > 0)
    smiles << ':' << a->atom_map();

//cerr << "Square bracket needed " << need_to_close_square_bracket << endl;
  if (need_to_close_square_bracket)
    smiles += ']';

  if (append_coordinates_after_each_atom)
    do_append_coordinates(smiles, a);

//cerr << "After addition, length " << smiles.length() << endl;

  return rc;
}

/*
  When we have a zero or singly connected atom, building a smiles is
  easier
*/

int
Molecule::_process_atom_for_smiles (Smiles_Formation_Info & sfi,
                                    IWString & smiles)
{
  if (! sfi.write_smiles())
    return _append_smarts_equivalent(sfi, smiles);

  atom_number_t zatom = sfi.zatom();

  const IWString * user_specified_atomic_smarts = sfi.user_specified_atomic_smarts();

  if (NULL == user_specified_atomic_smarts)
    ;
  else if (user_specified_atomic_smarts[zatom].length() > 0)
  {
    smiles << user_specified_atomic_smarts[zatom];
    return 1;
  }

  Atom * a = _things[zatom];

  const Element * e = a->element();

  if (a->permanent_aromatic())
    return append_permanent_aromatic(smiles, a);

// Look for unusual H counts. Remember, we must specify the H count for any
// charged atoms

  int hcount = -1;     // initialised to a negative number

  const int stored_hcount = include_hcount_in_smiles ? a->implicit_hydrogens() : 0;

  const formal_charge_t fc = a->formal_charge();

  if (! include_hcount_in_smiles)
    ;
  else if (0 == fc)
  {
    if (stored_hcount == _compute_implicit_hydrogens(zatom))     // the most common and normal case
      ;
    else if (a->implicit_hydrogens_known())      // probably a radical or something
      ;
    else
    {
      if (e->organic() && a->ncon() && file_scope_display_unusual_hcount_warning_messages)
        cerr << "Unusual hcount, atom " << zatom << "(" << a->atomic_symbol() << ") ncon = " << a->ncon() << " stored = " << stored_hcount << " implicit is " << _compute_implicit_hydrogens(zatom) << endl;
      hcount = stored_hcount;
    }
  }
  else if (stored_hcount)
    hcount = stored_hcount;

  int need_to_close_square_bracket = 0;

#ifdef DEBUG_PROCESS_ATOM_FOR_SMILES
  cerr << "Smiles so far '" << smiles << "'\n";
#endif

// Whenever a square bracket atom is encountered, the H count must always be explicit. 

  if (hcount >= 0 || 0 != fc || e->needs_square_brackets() || a->isotope() || a->implicit_hydrogens_known() || (a->atom_map() && _include_atom_map_with_smiles))
  {
    smiles += '[';
    need_to_close_square_bracket = 1;
    if (hcount < 0)
      hcount = stored_hcount;
  }

  int rc = e->append_smiles_symbol(smiles, NOT_AROMATIC, a->isotope());

  finish_smiles_atom(smiles, a->atomic_number(), hcount, fc);

  if (a->atom_map() && _include_atom_map_with_smiles)
    smiles << ':' << a->atom_map();

  if (need_to_close_square_bracket)
    smiles += ']';

  if (append_coordinates_after_each_atom)
    do_append_coordinates(smiles, a);

  return rc;
}

/*int
include_atom_in_smiles (const Molecule * m, atom_number_t a)
{
  assert (OK_ATOM_NUMBER (m, a));

  const Element * e = m->elementi (a);

  if (1 != e->atomic_number())
    return 1;
  if (0 != m->formal_charge (a))
    return 1;
  if (m->isotope (a))
    return 1;

  return 0;
}*/

int
Molecule::_unset_implicit_hydrogens_known_if_computed_matches()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (! a->implicit_hydrogens_known())
      continue;

    int h;
    if (! a->compute_implicit_hydrogens(h))
      continue;

    if (a->implicit_hydrogens() != h)    // computed and actual differ, must remain a special case
      continue;

    a->set_implicit_hydrogens_known(0);
    rc++;
  }

  return rc;
}

int
Molecule::_unset_all_implicit_hydrogens_known_attributes()
{
  int rc = 0;

//cerr << "Searching " << _number_elements << " atoms\n";
  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (! a->implicit_hydrogens_known())
      continue;

    a->set_implicit_hydrogens_known(0);
    int notused;
    (void) a->recompute_implicit_hydrogens(notused);
    rc++;
  }

//cerr << "Processed " << rc << " atoms\n";

  return rc;
}

Smiles_Formation_Info::Smiles_Formation_Info (int na, int nr) : _rnm (nr)
{
  _natoms = na;

  _already_done = NULL;

  _previous_atom = INVALID_ATOM_NUMBER;

  _zatom = INVALID_ATOM_NUMBER;

  _include_atom = NULL;

  _write_smiles = 1;

  _make_smarts_embedding = NULL;

  _user_specified_atomic_smarts = NULL;

  return;
}

Smiles_Formation_Info::~Smiles_Formation_Info()
{
  return;
}

int
Smiles_Formation_Info::ok() const
{
  return _rnm.ok();
}

int
Smiles_Formation_Info::make_smarts_embedding (atom_number_t zatom) const
{
  if (NULL == _make_smarts_embedding)
    return 0;

  return _make_smarts_embedding[zatom];
}

/*void
Smiles_Formation_Info::set_create_embedding_smarts (atom_number_t zatom,
                                                    int s)
{
  if (NULL == _create_embedding)
    _create_embedding = new_int(_natoms);

  _create_embedding[zatom] = s;

  return;
*/

Temporarily_Set_Include_Chiral_Info_in_Smiles::Temporarily_Set_Include_Chiral_Info_in_Smiles (int c)
{
  _initial_value = include_chiral_info_in_smiles();

  set_include_chiral_info_in_smiles(c);

  return;
}

Temporarily_Set_Include_Chiral_Info_in_Smiles::~Temporarily_Set_Include_Chiral_Info_in_Smiles()
{
  set_include_chiral_info_in_smiles(_initial_value);

  return;
}

void
reset_smiles_support_file_scope_variables ()
{
  reuse_ring_closure_numbers = 1;
  include_aromaticity_in_smiles = 0;
  _include_cis_trans_in_smiles = 1;
  _include_chiral_info_in_smiles = 1;
  _ignore_chiral_info_on_input = 0;
  append_coordinates_after_each_atom = 0;
  _write_smiles_with_smarts_atoms = 0;
  include_implicit_hydrogens_on_aromatic_n_and_p = 1;
  _add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens = 0;
  _write_single_bonds_in_smiles = 0;
  include_hcount_in_smiles = 1;
  file_scope_display_unusual_hcount_warning_messages = 1;
  sflag = 'K';

  return;
}
