#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
using std::cerr;
using std::endl;

// Be sure to define this symbol so all the private functions get defined

#define COMPILING_MDL_CC
#define COMPILING_CTB

#include "misc.h"
#include "iwstring_data_source.h"
#include "string_data_source.h"
#include "iwcrex.h"

#include "readmdl.h"
#include "mdl.h"
#include "atom_alias.h"
#include "molecule.h"
#include "misc2.h"
#include "chiral_centre.h"
#include "rwmolecule.h"
#include "aromatic.h"
#include "mdl_atom_record.h"

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

void
MDL_File_Supporting_Material::_default_values()
{
  _isis_standard_records = 0;
  _ignore_unrecognised_m_records = 0;
  _report_unrecognised_records = 1;
  _die_on_erroneous_m_input = 1;
  _write_mdl_dollar_signs = 1;
  _write_mdl_chiral_flags = 1;
  _write_v30_mdl_files = 0;
  _ignore_self_bonds = 0;
  _write_mdl_chiral_flags = 1;
  _include_chiral_info_in_mdl_outputs = 1;
  _mdl_read_h_correct_chiral_centres = 1;
  _mdl_write_h_correct_chiral_centres = 1;
  _insert_between_sdf_name_tokens = ' ';
  _extract_isis_extregno = 0;
  _fetch_all_sdf_identifiers = 0;
  _take_first_tag_as_name = 0;
  _prepend_sdfid = 1;
  _discard_sdf_molecule_name = 0;
  _multi_record_tag_data_present = 1;
  _mdl_write_aromatic_bonds = 0;
  _mdl_write_aromatic_atoms = 0;
  _display_non_organic_chirality_messages = 1;
  _mdl_display_invalid_chiral_connectivity = 1;
  _truncate_long_symbols = 0;
  _change_long_symbols_to.resize(0);
  _discern_chirality_from_wedge_bonds = 0;
  _write_isotopes_as_numbers_rather_than_differences_from_normal = 0;
  _write_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;
  _read_isotopes_as_numbers_rather_than_differences_from_normal = 0;
  _read_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;
  _replace_first_sdf_tag.resize(0);
  _allow_deuterium = 0;
  _allow_tritium = 0;
  if (NULL != _input_bond_type_translation_table)
  {
    delete [] _input_bond_type_translation_table;
    _input_bond_type_translation_table = NULL;
  }
  _mdl_g_records_hold_atom_symbols = 0;
  _write_mdl_charges_as_m_chg = 1;
  _set_elements_based_on_atom_aliases = 0;
  _write_Rn_groups_as_element = 0;
  _convert_single_atom_sgroup_to_element = 0;
  _gsub_mdl_file_data = ' ';
}

MDL_File_Supporting_Material::MDL_File_Supporting_Material ()
{
  _input_bond_type_translation_table = NULL;

  _digits2 = NULL;
  _digits3 = NULL;

  _default_values();
}

MDL_File_Supporting_Material::~MDL_File_Supporting_Material ()
{
  if (NULL != _input_bond_type_translation_table)
    delete [] _input_bond_type_translation_table;

  if (NULL != _digits2)
  {
    delete [] _digits2;
    delete [] _digits3;
  }

  return;
}

int
MDL_File_Supporting_Material::set_mdl_input_bond_type_translation (int zfrom, int zto)
{
  if (NULL == _input_bond_type_translation_table)
  {
    _input_bond_type_translation_table = new int[20];
    for (int i = 0; i < 20; i++)
    {
      _input_bond_type_translation_table[i] = i;
    }
  }

  assert (zfrom >= 0 && zfrom < 20);

  _input_bond_type_translation_table[zfrom] = zto;

  return 1;
}

static MDL_File_Supporting_Material default_mdl_file_optional_settings;

MDL_File_Supporting_Material *
global_default_MDL_File_Supporting_Material ()
{
  return &default_mdl_file_optional_settings;
}

int
ignore_self_bonds ()
{
  return default_mdl_file_optional_settings.ignore_self_bonds();
}


int
interpret_d_as_deuterium ()
{
  return default_mdl_file_optional_settings.allow_deuterium();
}

int
interpret_t_as_tritium ()
{
  return default_mdl_file_optional_settings.allow_tritium();
}

int
set_display_non_organic_chirality_messages ()
{
  return default_mdl_file_optional_settings.display_non_organic_chirality_messages();
}

void
set_write_isis_standard (int s)
{
  default_mdl_file_optional_settings.set_write_isis_standard(s);
}

void
set_write_mdl_charges_as_m_chg (int s)
{
  default_mdl_file_optional_settings.set_write_mdl_charges_as_m_chg(s);
}

/*
  This function converts from MDL charge codes, to charge
  codes used by this programme.

  Note that we have a special number for radicals
*/

int
convert_from_mdl_charge (int chg)
{
  if (0 == chg)     // or perhaps it really means -4!!! curse mdl
    return 0;

  if (4 == chg)
    return MDL_RADICAL;

  return 4 - chg;
}

int
convert_to_mdl_charge (int chg)
{
  switch(chg)
  {
    case 0:
      return 0;
    case 1:
      return 3;
    case 2:
      return 2;
    case 3:
      return 1;
    case 4:
      return 0;
    case -1: 
      return 5;
    case -2:
      return 6;
    case -3:
      return 7;
  }

  cerr << "convert to mdl_charge: bad charge " << chg << "\n";

  return 0;
}

/*
  A function which fills the array based on a buffer
  A typical buffer might look like

  A, A3, I4, I4, I4,

  M  CHG  4   7  -1   8  -1  11   1  12   1
  12345678901234567890123456789012345678901234567890123456789012345678901234567890
           1         2         3         4         5         6         7         8

  BEWARE, THEY LIE!!

  I've seen several cases where the count is incorrect (more tokens than
  are present). Guard against this!

  May 2002. the desktop version produces a different number of spaces than the host
  version. Change to just extracting tokens - this will presumably break if atom
  numbers go above 999. 
*/

int
fill_atom_property_array (const IWString & buffer,
                          int & npairs,
                          Aprop * atom_properties)
{
  assert (buffer.starts_with('M'));

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  int nw = buffer.nwords();

  npairs = nw - 3;      // discount M  XXX c tokens

  if (npairs != npairs / 2 * 2)    // must be an even number of tokens
  {
    cerr << "mdl:must have even pairs '" << buffer << "'\n";

    if (mdlfos->die_on_erroneous_m_input())
      return 0;

    npairs--;    // convert to even number

    if (npairs <= 0)    // ignoring now empty record
      return 1;
  }

  npairs = npairs / 2;

  if (0 == npairs)
  {
    cerr << "mdl:very strange, M record with no info '" << buffer << "'\n";
    return 0;
  }

  if (npairs > MAX_PAIRS)
  {
    cerr << "MDL::fill_atom_property_array:too many pairs on record, cannot process. See Ian " << npairs << endl;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);    // M
  buffer.nextword(token, i);    // XXX
  buffer.nextword(token, i);    // count

  for (int j = 0; j < npairs; j++)
  {
    buffer.nextword(token, i);
    
    if (! token.numeric_value(atom_properties[j]._atom_number))
    {
      cerr << "mdl:invalid atom number '" << buffer << "'\n";
      return 0;
    }

    buffer.nextword(token, i);

    if (! token.numeric_value(atom_properties[j]._property))
    {
      cerr << "mdl: invalid numeric value '" << buffer << "'\n";
      return 0;
    }
  }
  
  return npairs;
}

/*
  Remember that the numbers in M  ISO records are differences from the
  normal mass
*/

int
Molecule::mdl_add_m_isotope (int ntokens,
                             const Aprop * atom_properties)
{
  assert (ntokens > 0);

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    int k = atom_properties[i]._property;

//  cerr << "Isotope for atom " << j << " value " << k << "\n";
    if (0 == j && 0 == k)     // have seen this a couple of times. Hard to know what they intended
      continue;

    if (0 == k)     // have seen this one too!
      continue;

    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_isotope: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }

    Atom * a = _things[j];             // the j'th atom of this molecule

    if (mdlfos->read_M_isotopes_as_numbers_rather_than_differences_from_normal())
      a->set_isotope(k);
    else
      a->set_isotope(a->element()->normal_isotope() + k);
  }
  
  return 1;
}

/*
  Warning, this implementation is not complete. The definition from
  MDL says that if the M  CHG record is present, 0 charge is then
  forced on all atoms not mentioned. I'm not doing that. Change if
  this becomes a problem.
*/

int
Molecule::mdl_add_m_formal_charge (int ntokens,
                                   const Aprop * atom_properties)
{
  assert (ntokens > 0);

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_charge: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }
    
    if (! reasonable_formal_charge_value(atom_properties[i]._property))
    {
      cerr << "mdl add formal charge: unreasonable charge value " << atom_properties[i]._property << endl;
      return 0;
    }

    _things[j]->set_formal_charge(atom_properties[i]._property);
  }
  
  return 1;
}

int
Molecule::mdl_add_m_radical (int ntokens,
                                             const Aprop * atom_properties)
{
  assert (ntokens > 0);

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_radical: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }

    _things[j]->set_implicit_hydrogens(0, 1);
    _things[j]->set_implicit_hydrogens_known(1);
  }
  
  return 1;
}



/*
  We have an atom with > 2 connections. We need to know whether there are two
  neighbouring atoms that could be delocalised.
  The two neighbours must be set in AROMATIC_ATOMS, and must share the
  same atomic number and connectivity

  I needed this because we may have the case of a delocalised Carboxyllic
  acid adjacent to an aromatic ring. The entry in AROMATIC_ATOMS for the
  carbon in the ring may be set.

  c1ccccc1C(=O)O

  Mar 2001. Reading charged sulphonic acids from mol2 files, you can have
  three neighbours

  Oct 2001. Encountered O--S(--O)(--O)--O

  March 2006. Allow PO2 and PO3
*/

int
Molecule::_has_delocalised_neighbours (atom_number_t zatom,
                                       const int * aromatic_atoms,
                                       const int * aromatic_bonds,
                                       Set_of_Atoms & s,
                                       const int * ab) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();
  s.resize(acon);

  int invariant = 0;    //  a quick and dirty computation to see if two atoms are the same

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (aromatic_atoms[j])    // good, will process
      ;
    else if (ab[zatom * _number_elements + j])   // aromatic bond btw zatom and j, process
      ;
    else
      continue;

    const Atom * aj = _things[j];

    int my_invariant;

    if (8 == aj->atomic_number())
      my_invariant = 10 * aj->atomic_number() + aj->ncon();   // breaks if we get more than 10 connections
    else if (7 == aj->atomic_number() && aj->ncon() < 3)   // guanidines
      my_invariant = 10 * aj->atomic_number();
    else
      continue;

//  cerr << "Atom " << j << " is an aromatic neighbour, invariant " << my_invariant << endl;

    if (0 == s.number_elements())    // first one
    {
      invariant = my_invariant;
      s.add(j);
      continue;
    }
    else if (my_invariant != invariant)    // different atom type,
      continue;                            // or should this be a failure

    if (0 != aj->formal_charge())    // make sure anything with an existing formal charge up front
      s.insert_at_beginning(j);
    else
      s.add(j);
  }

//cerr << "Found " << s.number_elements() << " aromatic neighbours\n";

  if (0 == s.number_elements())
    return 0;

  if (NULL == aromatic_bonds)
    return(s.number_elements() > 1);

//  Implement this sometime, maybe not needed

#ifdef IS_THIS_NEEDED
  int nb = nedges();

  for (int i = 0; i < nb; i++)
  {
    if (0 == aromatic_bonds[i])
      continue;

    const Bond * b = _bond_list[i];

    if (! b->involves(zatom))
      continue;
  }
#endif

  return 1;
}

int
Molecule::_unset_aromatic_bond_entry(atom_number_t a1, 
                                     atom_number_t a2,
                                     int * aromatic_bonds) const
{
  if (NULL == aromatic_bonds)
    return 1;

  int b = which_bond(a1, a2);

  assert (b >= 0);

  aromatic_bonds[b] = 0;

//cerr << " Bond between " << a1 << " and " << a2 << " is bond " << b << endl;

  return 1;
}

int
Molecule::_unset_aromatic_bond_entry(atom_number_t a1, 
                                     atom_number_t a2,
                                     int * aromatic_bonds,
                                     const int * ab) const
{
  const int bnumber = ab[a1 * _number_elements + a2] - 1;

  if (bnumber < 0)
    return 1;

  aromatic_bonds[bnumber] = 0;

//cerr << " Bond between " << a1 << " and " << a2 << " is bond " << b << endl;

  return 1;
}

/*
  For each bond that is aromatic, mark the corresponding atoms
  as having aromatic connections
*/

int
Molecule::_identify_atoms_at_ends_of_aromatic_bonds (const int * aromatic_bond,
                                                     int * aromatic_connections,
                                                     int * ab) const
{
  set_vector(aromatic_connections, _number_elements, 0);

  const int nb = nedges();

  for (int i = 0; i < nb; ++i)
  {
    if (! aromatic_bond[i])
      continue;

    const Bond * b = _bond_list[i];

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    aromatic_connections[a1]++;
    aromatic_connections[a2]++;

    ab[a1 * _number_elements + a2] = ab[a2 * _number_elements + a1] = i+1;    // we create an index into the bond list
  }

  return 1;
}

int
Molecule::_count_elements (const atomic_number_t z,
                           const Set_of_Atoms & e) const
{
  int rc = 0;

  for (int i = 0; i < e.number_elements(); ++i)
  {
    if (z == _things[e[i]]->atomic_number())
      rc++;
  }

  return rc;
}

//#define DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS

/*
  Sept 97, generalise this to be any two equivalent looking 'aromatic'
  bonds bonded to the same atom.

  Mar 2001, update to include sulph* acids with more than 2 delocalised
  oxygens

  Feb 2009. ALso do Nitro groups

  May 2015. Guanidine groups
*/

int
Molecule::process_delocalised_carbonyl_bonds (int * aromatic_atoms,
                                              int * aromatic_bonds)
{
  if (NULL == aromatic_bonds)    // cannot do anything
    return 1;

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
  if (NULL != aromatic_bonds)
  {
    cerr << "On entry to _process_delocalised_carbonyl_bonds\n";
    iw_write_array(aromatic_bonds, nedges(), "aromatic bonds", cerr);
  }
#endif

  int * rm = new int[_number_elements + _number_elements];std::unique_ptr<int[]> free_rm(rm);

  ring_membership(rm);

  int * number_aromatic_connections = rm + _number_elements;

// Transfer information from the aromatic bonds array to an array that allows us to quickly find if two atoms have an aromatic bond between them

  int * ab = new int[_number_elements * _number_elements]; std::unique_ptr<int[]> free_ab(ab);

  set_vector(ab, _number_elements * _number_elements, 0);

  _identify_atoms_at_ends_of_aromatic_bonds(aromatic_bonds, number_aromatic_connections, ab);

  for (int i = 0; i < _number_elements; i++)
  {
    if (rm[i])    // only process non ring atoms here
      continue;

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
    cerr << " atom " << i << " aromatic? " << aromatic_atoms[i] << " aromatic connections " << number_aromatic_connections[i] << endl;
#endif

    if (aromatic_atoms[i])
      ;
    else if (number_aromatic_connections[i] > 1)
      ;
    else
      continue;

    Atom * a = _things[i];

    int icon = a->ncon();

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
    cerr << "Molecule::_process_delocalised_carbonyl_bonds: atom " << i << " '" << a->atomic_symbol() << "' " << icon << " connections is non ring\n";
#endif

    if (icon < 2)
      continue;

    Set_of_Atoms s;

    if (! _has_delocalised_neighbours(i, aromatic_atoms, aromatic_bonds, s, ab))
      continue;

    if (1 == s.number_elements())   // deal with these later
      continue;

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
    cerr << "atom " << i << " atoms " << s << " delocalised\n";
#endif

    aromatic_atoms[i] = 0;
    s.set_vector(aromatic_atoms, 0);

    atom_number_t a0 = s[0];

#ifdef NEED_TO_PROCESS_AROMATIC_PO_FROM_TRIPOS
// Jan 2005.  Sybyl produces P:O bonds.  No, just one bond in a PO3
// grouping is aromatic, so this isn't needed

    if (1 == s.number_elements())
    {
      if (a->implicit_hydrogens() && _things[a0]->implicit_hydrogens())
      {
        set_bond_type_between_atoms(i, a0, DOUBLE_BOND);
        a->set_modified();
        _things[a0]->set_modified();
        if (NULL != aromatic_bonds)
        {
          int b = which_bond(i, a0);
          aromatic_bonds[b] = 0;
        }
        continue;
      }

      cerr << "Molecule::_process_delocalised_carbonyl_bonds:no available bonds\n";
      continue;
    }
#endif

    _unset_aromatic_bond_entry(i, a0, aromatic_bonds, ab);

//  Special case for the Nitro group

    if (7 == a->atomic_number() && 2 == s.number_elements() && 8 == atomic_number(a0) && 8 == atomic_number(s[1]))
    {
      set_bond_type_between_atoms(i, a0, DOUBLE_BOND);
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      _unset_aromatic_bond_entry(i, s[1], aromatic_bonds, ab);
      continue;
    }

    set_bond_type_between_atoms(i, a0, SINGLE_BOND);

    if (6 == a->atomic_number() && 3 == s.number_elements() && 3 == _count_elements(7, s))   // guanidine
    {
      _unset_aromatic_bond_entry(i, s[1], aromatic_bonds, ab);
      _unset_aromatic_bond_entry(i, s[2], aromatic_bonds, ab);
      if (1 == _things[s[1]]->ncon())
      {
        set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
        set_bond_type_between_atoms(i, s[2], SINGLE_BOND);
      }
      else if (1 == _things[s[2]]->ncon())
      {
        set_bond_type_between_atoms(i, s[1], SINGLE_BOND);
        set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      }
      else if (_things[s[1]]->nbonds() < _things[s[2]]->nbonds())
      {
        set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
        set_bond_type_between_atoms(i, s[2], SINGLE_BOND);
      }
      else     // could be dangerous, perhaps one of these atoms already has a double bond. ignore for now...
      {
        set_bond_type_between_atoms(i, s[1], SINGLE_BOND);
        set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      }
      continue;
    }

    if (-1 != _things[a0]->formal_charge())
      _things[a0]->set_formal_charge(-1);

    if (s.number_elements() > 1)
      _unset_aromatic_bond_entry(i, s[1], aromatic_bonds, ab);

    if (2 == s.number_elements())    // the most common case, carbonyl for example
    {
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      continue;
    }

    _unset_aromatic_bond_entry(i, s[2], aromatic_bonds, ab);

//  SO3 gets two double bonds, PO3 gets just one

    if (3 == s.number_elements())
    {
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      if (16 == a->atomic_number())
        set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      else
        set_bond_type_between_atoms(i, s[2], SINGLE_BOND);
      continue;
    }

//  SO3 gets two double bonds, PO3 gets just one

    if (4 == s.number_elements())
    {
      set_bond_type_between_atoms(i, s[1], SINGLE_BOND);
      _things[s[1]]->set_formal_charge(-1);
      set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      set_bond_type_between_atoms(i, s[3], DOUBLE_BOND);
      _unset_aromatic_bond_entry(i, s[3], aromatic_bonds, ab);
    }
    else
      cerr << "Molecule::_process_delocalised_carbonyl_bonds: This many neighbours " << s << endl;
  }

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
   cerr << "_process_delocalised_carbonyl_bonds returning\n";
#endif
  return 1;
}

//#define DEBUG_FINAL_PROCESSING_AROMATIC_MDL

int
Molecule::_final_processing_of_aromatic_mdl_input(int * aromatic_atoms,
                                int * aromatic_bonds)
{
#ifdef DEBUG_FINAL_PROCESSING_AROMATIC_MDL
  cerr << "Molecule::_final_processing_of_aromatic_mdl_input:\n";
#endif

  if (allow_delocalised_carbonyl_bonds())
  {
    if (! process_delocalised_carbonyl_bonds(aromatic_atoms, aromatic_bonds))
    {
      cerr << "Cannot process delocalised carbonyl bonds\n";
      return 0;
    }
  }

  if (locate_item_in_array(1, _number_elements, aromatic_atoms) < 0)    // there were no aromatic atoms in the input
    return 1;

  int rc = 1;

  if (! find_kekule_form(aromatic_atoms))
  {
    cerr << "read_molecule_mdl_ds: find kekule form failed";

    if (allow_input_without_valid_kekule_form())
    {
      cerr << ", using single bonds";
      _molecule_name += " (invalid KEKULE form)";
    }
    else
      rc = 0;

    cerr << ' ' << _molecule_name << "'\n";
  }

  return rc;
}

/*
  A frequent operation is to decide whether a problem with bad chirality
  should be fatal or not
*/

int
return_code_depending_on_ignore_incorrect_chiral_input()
{
  if (ignore_incorrect_chiral_input())
  {
    cerr << "Ignored\n";
    return 1;
  }

  return 0;
}

/*
  Shared between the V2 and V3 programmes
*/

Atom *
MDL_File_Supporting_Material::create_mdl_atom (const const_IWSubstring & ss,
                 int msdif,
                 int chg,
                 int is_radical) const
{
  const_IWSubstring zsymbol = ss;     // our own local copy

  if (ss.length() > 2)    // a residue or something
  {
    cerr << "Molecule::create_mdl_atom:warning: element '" << ss << "' encountered\n";

    if (_change_long_symbols_to.length())
      zsymbol = _change_long_symbols_to;
    else if (_truncate_long_symbols)
      zsymbol.iwtruncate(2);
    else if (atomic_symbols_can_have_arbitrary_length())
      ;
    else               // no way of processing it
      return NULL;
  }

  int iso = 0;     // never used with mdl files

  const Element * e = get_element_from_symbol(zsymbol, iso);
  if (NULL == e)
  {
    if (_allow_deuterium && 'D' == zsymbol)
    {
      e = get_element_from_symbol_no_case_conversion("H");
      iso = 2;
    }
    else if (_allow_tritium && 'T' == zsymbol)
    {
      e = get_element_from_symbol_no_case_conversion("H");
      iso = 3;
    }
    else if (! auto_create_new_elements())
    {
      cerr << "create_mdl_atom: unrecognised element '" << zsymbol << "'\n";
      return NULL;
    }
    else
      e = create_element_with_symbol(zsymbol);

    if (NULL == e)
    {
      cerr << "MDL_File_Optional_Settings::create_mdl_atom:cannot create element from '" << zsymbol << "'\n";
      return NULL;
    }
  }

  Atom * rc = new Atom(e);
  assert (NULL != rc);

  if (0 == msdif)
    ;
  else if (_read_isotopes_as_numbers_rather_than_differences_from_normal)
    iso = msdif;
  else
    iso = e->normal_isotope() + msdif;

  if (iso)
    rc->set_isotope(iso);

  if (chg)
  {
    if (! reasonable_formal_charge_value(chg))
    {
      cerr << "create mdl atom: bad charge value " << chg << endl;
      return NULL;
    }
    rc->set_formal_charge(chg);
  }

  if (is_radical)
  {
    rc->set_radical(1);
    rc->set_implicit_hydrogens(0, 1);
    rc->set_implicit_hydrogens_known(1);
  }

  return rc;
}

/*static int
parse_xyz_symbol_etc (const const_IWSubstring & buffer,
                      coord_t & x,
                      coord_t & y,
                      coord_t & z,
                      const_IWSubstring & atomic_symbol,
                      int & msdiff,
                      int & chg,
                      int & astere)
{
  int nw = buffer.nwords();

  if (nw < 4)
  {
    cerr << "parse_xyz_symbol_etc:too few tokens\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(x))
  {
    cerr << "Bad X value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(y))
  {
    cerr << "Bad Y value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(z))
  {
    cerr << "Bad Z value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (atomic_symbol, i);

  if (! buffer.nextword (token, i))
    return 1;

  if (! token.numeric_value(msdiff))
  {
    cerr << "Bad msdiff value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword (token, i))
    return 0;

  if (! token.numeric_value(chg))
  {
    cerr << "Bad chg value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword (token, i))
    return 1;

  if (! token.numeric_value(astere))
  {
    cerr << "Bad astere value '" << token << "'\n";
    return 0;
  }

  return 1;
}*/



/*
  When reading an atom record we have encountered a chirality flag
*/

int
Molecule::_mdl_atom_is_chiral_centre (atom_number_t zatom, int cfg,
                                      MDL_File_Supporting_Material & mdlfos)
{
  const Atom * a = _things[zatom];

  if (! a->element()->organic())
  {
    if (3 == cfg)     // unspecified, who cares
      return 1;

    if (mdlfos.display_non_organic_chirality_messages())
      cerr << "Ignoring non-organic chirality, atom " << zatom << " type '" << a->atomic_symbol() << "'\n";
    return 1;
  }

  if (1 == cfg || 2 == cfg)
  {
    Chiral_Centre * c = create_chiral_centre(zatom, 1);
    if (NULL == c)
    {
      cerr << "Cannot create chiral center at atom " << zatom << endl;
      return 0;
    }
    c->set_chirality_known(cfg);
  }
  else if (3 == cfg)    // maybe ignore these
    mdlfos.unspecified_chiral_atom_if_interested(zatom);
  else
  {
    cerr << "Molecule::_mdl_atom_is_chiral_centre: unrecognised chirality " << cfg << endl;
    return 0;
  }

  return 1;
}

/*
  We have read a bond record and DIRECTIONALITY is the 4'th number from the
  MDL bond record. Set the appropriate directionality flags in the appropriate
  bond
*/

int
Molecule::_mdl_set_bond_directionality (atom_number_t a1, atom_number_t a2,
                                        int directionality)
{
  assert (0 != directionality);

  for (int i = _bond_list.number_elements() - 1; i >= 0; i--)    // chances are this is the last bond in the list
  {
    Bond * b = _bond_list[i];
    if (! b->involves(a1, a2))
      continue;

    if (b->is_wedge_any())
    {
      cerr << "Molecule::_mdl_set_bond_directionality: bond between atoms " << a1 << " and " << a2 << " already directional\n";
      continue;
    }

    if (1 == directionality)
      b->set_wedge_up();
    else if (6 == directionality)
      b->set_wedge_down();
    else if (4 == directionality)
      b->set_wedge_either();
    else if (3 == directionality)
      b->set_cis_trans_either_double_bond();
    else
    {
      cerr << "Molecule::_mdl_set_bond_directionality: unknown directionality " << directionality << endl;
      return 0;
    }

    return 1;
  }

  cerr << "Molecule::_mdl_set_bond_directionality: HUH, no bond between atoms " << a1 << " and " << a2 << endl;
  return 0;
}

/*
  User visible numbering for bonds is
    single = 1
    double = 2
    triple = 3
    arom   = 4

  This converts those numbers into bond_type_t forms.
*/

int
convert_from_mdl_number_to_bond_type (int int_rep, bond_type_t & bt)
{
  switch(int_rep)
  {
    case 1:
      bt = SINGLE_BOND;
      return 1;

    case 2:
      bt = DOUBLE_BOND;
      return 1;

    case 3:
      bt = TRIPLE_BOND;
      return 1;

    case 4:
      bt = AROMATIC_BOND;
      return 1;

    default:
      cerr << "convert_to_bond_type: unrecognised type " << int_rep << endl;
      return 0;
  }

  assert (NULL == "Should not come here");
  return 0;
}

int
MDL_File_Supporting_Material::parse_bond_record (const_IWSubstring & buffer,
                   int na,
                   atom_number_t & a1, atom_number_t & a2,
                   int & bond_type_read_in,
                   int & directionality)
{
  if (buffer.length() < 9)
  {
    cerr << "Invalid bond record '" << buffer << "'\n";
    return 0;
  }

  directionality = 0;

  int ntokens = int3d(buffer, a1, a2, &bond_type_read_in);

  a1--;       // our atom numbers start at 0
  a2--;

  if ( (3 != ntokens) || a1 < 0 || a1 >= na || a2 < 0 || a2 >= na ||
       ((a1 == a2) && ! _ignore_self_bonds))
  {
    cerr << "parse_bond_record: error on bond record\n";

    if (3 != ntokens)
      cerr << "Bad token count on bond record " << ntokens << endl;
    else if (a1 < 0 || a1 > na)
      cerr << "Bad a1 value " << a1 << ", natoms = " << na << endl;
    else if (a2 < 0 || a2 > na)
      cerr << "Bad a2 value " << a2 << ", natoms = " << na << endl;
    else
      cerr << "a1 = " << a1 << " a2 = " << a2 << " natoms = " << na << endl;

    return 0;
  }

  if (buffer.length() < 12)       // no directionality present
    return 1;

  char direction = buffer[11];     // the 12'th column - probably should check columns 10 and 11....

  if ('0' == direction)            // no directionality
    return 1;

  directionality = direction - '0';    // the numeric form we return in our argument list

  if (0 == _accumulate_mdl_chirality_features)
    return 1;

  if ('1' == direction)
    _up_bonds_last_molecule.add(a1);
  else if ('6' == direction)
    _down_bonds_last_molecule.add(a1);
  else if ('4' == direction)
    _squiggle_bonds_last_molecule.add(a1);
  else if ('3' == direction)
  {
    assert (2 == bond_type_read_in);
    _unspecified_double_bond_atoms.add(a1);
    _unspecified_double_bond_atoms.add(a2);
  }
  else
  {
    cerr << "parse_bond_record: unrecognised bond directionality '" << direction << "'\n";
    return 0;
  }

  return 1;
}

/*
  We always attempt to read NB records.
*/

int
parse_m_sty_record (const const_IWSubstring & buffer,
                    int & sgroups_present)
{
  if (buffer.nwords() < 3)
  {
    cerr << "Molecule::_read_molecule_mdl_trailing_records:invalid record '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring n;
  buffer.word(2, n);

  int tmp;

  if (! n.numeric_value(tmp) || tmp < 1)
  {
    cerr << "Molecule::parse_m_sty_record:invalid count '" << buffer << "'\n";
    return 0;
  }

  sgroups_present += tmp;

  return 1;
}

int
extract_sdf_identifier(const IWString & buffer, IWString & id)
{
  int istart = buffer.index('<');
  int istop = buffer.rindex('>');

  if (istart > istop || istart < 0 || istop < 0)
  {
    cerr << "extract_sdf_identifier: HUH: '" << buffer << "'\n";
    return 0;
  }

  buffer.from_to(istart + 1, istop - 1, id);

  if (id.contains(' '))
    id.gsub(' ', '_');

  return 1;
}

/*
  Nov 99. GER needed to flag superatoms. These are things like

OMe
  -ISIS-  10279910532D

  4  3  0  0  0  0  0  0  0  0999 V2000
    7.5837   -2.1040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6320   -3.1626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9044   -2.4178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9458   -3.4454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
A    4
MeO
M  END
>  <REGNO>(2)
2

$$$$

  which is an example from Jeff Hanson, 27 Oct 99

  Note that we also flag records starting with a G
*/

static int a_records_found = 0;
static int g_records_found = 0;

static resizable_array_p<Atom_Alias> aliases;

int
a_records_found_most_recent_mdl_molecule()
{
  return a_records_found;
}

int
g_records_found_most_recent_mdl_molecule()
{
  return g_records_found;
}

const resizable_array_p<Atom_Alias> &
atom_aliases_most_recent_mdl_molecule()
{
  return aliases;
}

void
add_to_text_info (resizable_array_p<IWString> & text_info, const const_IWSubstring & zextra)
{
  IWString * s = new IWString(zextra);

  text_info.add(s);

  return;
}

/*
  Sept 2010. Turns out mdl is now putting digits in these records, so we can have. Hmmm,
  this showed up earlier, not sure why we did not find problems earlier.

> 2  <MODEL.EXTREG>

*/

#ifdef LOOKS_LIKE_SDF_TAG_OV
static int
looks_like_sdf_tag (const const_IWSubstring & buffer)
{
//cerr << "Checking for sdf tag '" << buffer << "'\n";

  if (! buffer.starts_with('>'))
    return 0;

  int n = buffer.length();

  if (n < 3)
    return 0;

  for (int i = 1; i < n; i++)
  {
    char c = buffer[i];

    if (' ' == c)
      continue;

    if ('<' != c)
      return 0;

    do
    {
      if ('>' == buffer[i])
        return 1;

      i++;
    } while (i < n);

    return 0;
  }

  return 0;
}
#endif

int
looks_like_sdf_tag (const const_IWSubstring & buffer)
{
//cerr << "Checking for sdf tag '" << buffer << "'\n";

  if (! buffer.starts_with('>'))   // cheapest test first
    return 0;

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  if (! buffer.nextword(token, i))   // strange, no second token on line!
    return 0;

  if (token.starts_with('<'))   // the most common case
    ;
  else         // maybe a digit in there
  {
    int notused;
    if (! token.numeric_value(notused) || notused < 0)
      return 0;

    if (! buffer.nextword(token, i))   // if we have a digit, there must be a following token
      return 0;

    if (! token.starts_with('<'))   // following the digit, there must be the 
      return 0;
  }

  if (token.ends_with('>'))
    return 1;

  while (buffer.nextword(token, i))
  {
    if (token.ends_with('>'))
      return 1;
  }

  return 0;
}


/*
  the MDL_Molecule class inherits from a Molecule. It needs to
  examine all the M records looking for substructural info.
  This function will process any records associated with the
  Molecule class
*/

int
Molecule::_common_parse_M_record (const const_IWSubstring & buffer,
                                  int & fatal)
{
  Aprop atom_properties[MAX_PAIRS];

#ifdef DEBUG_COMMON_PARSE_M_RECORD
  cerr << "_common_parse_M_record examining '" << buffer << "'\n";
#endif

  if (buffer.starts_with("M  CHG"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_formal_charge(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  ISO"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_isotope(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  RAD"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_radical(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  RGP"))
  {
    if (! _parse_M_RGP_record(buffer))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

#ifdef DEBUG_COMMON_PARSE_M_RECORD
  cerr << "_common_parse_M_record:did not recognise '" << buffer << "'\n";
#endif

// Not recognised here

  fatal = 0;
  return 0;
}


/*
  When reading RDFILES, the user can specify which fields to use for the name
  These are processed in order
*/

int
MDL_File_Supporting_Material::add_rdfile_identifier (const IWString & new_identifier)
{
  assert(new_identifier.length() > 0);

  return _rdfile_identifiers.add(new IWString(new_identifier));
}

int
MDL_File_Supporting_Material::set_rdfile_start_of_record (const const_IWSubstring & s)
{
  _rdfile_start_of_record = s;

//cerr << "rdfile_start_of_record set to '" << s << "'\n";

  return 1;
}

static int
skip_to_rdfile_start_of_record (iwstring_data_source & input,
                                const IWString & rdfile_start_of_record)
{
  const_IWSubstring buffer;
  
  int records_read_here = 0;

  while (input.next_record(buffer))
  {
    records_read_here++;

    if (buffer.starts_with(rdfile_start_of_record))
    {
      input.push_record();
      return 1;;
    }
  }

  if (0 == records_read_here)
  {
    cerr << "read mol rdf eof\n";
    return 0;
  }

  cerr << "EOF reading RDFILE, cannot find start record '" << rdfile_start_of_record << "', tried " << records_read_here << "\n";
  return 0;
}

//#define DEBUG_READ_MOLECULE_RDF_DS

/*
  Lots of ways an RDFILE identifier can match
*/

static int
rdfile_record_matches (const IWString & buffer,
                       const IWString & rdfile_identifier)
{
#ifdef DEBUG_READ_MOLECULE_RDF_DSQ
  cerr << "Checking rdfile match '" << buffer << "' vs '" << rdfile_identifier << "'\n";
#endif

  if (buffer.starts_with(rdfile_identifier))
    return 1;

  if (buffer.starts_with("$DTYPE MOL:") && buffer.matches_at_position(11, rdfile_identifier))
    return 1;

  if (buffer.starts_with("$DTYPE ") && buffer.matches_at_position(7, rdfile_identifier))
      return 1;

  return 0;
}

/*
  An RDFILE is somewhat different, in that it has a bunch of records starting with $

  The connection table starts with '$DATUM $MFMT'
  Also have seen files where it starts with '$MFMT $MIREG' sometimes.
  We must tell read_molecule_mdl_ds to return on encountering 'M  END'
*/

int
Molecule::read_molecule_rdf_ds (iwstring_data_source & input)
{
  assert(ok());
  assert(input.good());

  if (input.eof())
    return 0;

#ifdef DEBUG_READ_MOLECULE_RDF_DSQ
  cerr << "Molecule::read_molecule_rdf_ds:input.record_buffered? " << input.record_buffered() << endl;
#endif

  MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  const IWString & rdfile_start_of_record = mdlfos->rdfile_start_of_record();

  if (rdfile_start_of_record.length())   // scan through file till we get to header
  {
    if (! skip_to_rdfile_start_of_record(input, rdfile_start_of_record))
      return 0;
  }

  IWString possible_name;

  if (! _read_molecule_rdf_ds(input, possible_name, *mdlfos))
    return 0;

  if (0 == _molecule_name.length() && possible_name.length())
    set_name(possible_name);

  return 1;
}

int
Molecule::_read_molecule_rdf_ds (iwstring_data_source & input,
                                 IWString & possible_name,
                                 MDL_File_Supporting_Material & mdlfos)
{
  const resizable_array_p<IWString> & rdfile_identifiers = mdlfos.rdfile_identifiers();

  int nrdfid = rdfile_identifiers.number_elements();

  const IWString rdfile_start_of_record = mdlfos.rdfile_start_of_record();

// I've seen many kinds of RDF files. Some have the name in the name
// field, others have it as a data item of type LILLYNUMBER. Look for
// LILLYNUMBER things and save the next record.

  int records_read_here = 0;

  IWString buffer;    // must be an IWString, NOT a const_IWSubstring

  int have_read_structure = 0;

  while (input.next_record(buffer))
  {
    records_read_here++;

#ifdef DEBUG_READ_MOLECULE_RDF_DS
    cerr << "Just read '" << buffer << "'\n";
#endif

    if (buffer == "$DATUM $MFMT" || buffer.starts_with("$MFMT $MIREG"))
    {
      if (have_read_structure)
      {
        cerr << "Molecule::read_molecule_rdf_ds:multiple structures!!\n";
        return 0;
      }

#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Start reading mdl connection table at " << input.lines_read() << endl;
#endif

      if (! read_molecule_mdl_ds(input, 1))
        return 0;

      have_read_structure = 1;

#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Successfully read structure\n";
#endif

      if (0 == rdfile_start_of_record.length())
        return 1;
    }
    else if (nrdfid)
    {
      for (int i = 0; i < nrdfid; i++)
      {
        const IWString & rdfile_identifier = *(rdfile_identifiers[i]);

        if (rdfile_record_matches(buffer, rdfile_identifier))
        {
          const_IWSubstring tmp;
          if (! input.next_record(tmp))
          {
            cerr << "Molecule::read_molecule_rdf_ds: premature eof\n";
            return 0;
          }

          if (tmp.starts_with("$DATUM"))
            tmp.remove_leading_words(1);

          possible_name.append_with_spacer(rdfile_identifier);
          possible_name << ':' << tmp;
        }
      }
    }
    else if (buffer.contains("$DTYPE ") && buffer.contains(":LILLYNUMBER"))
    {
     (void) (input.next_record(possible_name));
      if (! possible_name.contains("$DATUM "))    // strip this some time.
      {
        cerr << "Molecule::read_molecule_rdf_ds: Yipes, no $datum '" << possible_name << "'\n";
        return 0;
      }
    }

#ifdef DEBUG_READ_MOLECULE_RDF_DS
    cerr << "Checking start '" << rdfile_start_of_record << "' and '" << buffer << "'\n";
#endif
    if (rdfile_start_of_record.length() && buffer.starts_with(rdfile_start_of_record) && records_read_here > 1)
    {
#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Found rdf start of record, structure? " << have_read_structure << endl;
#endif

      if (! have_read_structure)
        cerr << "Molecule::read_molecule_mdl_ds:no structure!\n";

      input.push_record();
      return 1;
    }
  }

  if (0 == records_read_here || have_read_structure)
  {
    cerr << "read mol rdf eof\n";
    return 0;
  }

  if (have_read_structure)
    return 1;

  if (possible_name.length())
    return 1;

  cerr << "EOF reading RDFILE\n";
  return 0;
}

int
Molecule::write_molecule_mdl (const char * fname, const char * comments) const
{
  assert(ok());
  assert (NULL != fname);

  std::ofstream output(fname);
  if (! output.good())
  {
    cerr << "Molecule::write_molecule_mdl: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_molecule_mdl(output, comments);
}

/*
  Function to write an MDL molfile.
*/

int
Molecule::write_extra_text_info (IWString & buffer) const
{
  int ne = _text_info.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const IWString * info = _text_info[i];

    buffer <<(*info) << newline_string();
  }

  return buffer.length();
}


/*
  Return the digit - 1,2,3 that is the stereo flag for a given atom
*/

int
Molecule::_mdl_atom_stereo_value (atom_number_t zatom,
                                  const MDL_File_Supporting_Material & mdlfos) const
{
  Chiral_Centre * c = chiral_centre_at_atom(zatom);

  if (NULL == c)
    return 0;

  if (! c->chirality_known())
    return 3;

// If there is an explicit Hydrogen, that atom must be the highest numbered atom
// Who knows what you are supposed to do if there is an implicit and an explicit Hydrogen

  if (c->implicit_hydrogen_count())        // we assume no explicit hydrogen also
    return c->mdl_stereo_centre_value();

  if (c->lone_pair_count())                // who knows how we are supposed to handle a lone pair and an explicit hydrogen
    return c->mdl_stereo_centre_value();

  if (! mdlfos.mdl_write_h_correct_chiral_centres())    // optional behaviour
    return c->mdl_stereo_centre_value();

// Look for an explicit Hydrogen atom. Remember, we pass the atoms SW N SE

  if (1 == _things[c->top_front()]->atomic_number())
    return c->mdl_stereo_centre_value(c->right_down(), c->top_back(), c->left_down());
  if (1 == _things[c->top_back()]->atomic_number())
    return c->mdl_stereo_centre_value(c->left_down(), c->top_front(), c->right_down());
  if (1 == _things[c->left_down()]->atomic_number())
    return c->mdl_stereo_centre_value(c->right_down(), c->top_front(), c->top_back());
  if (1 == _things[c->right_down()]->atomic_number())
    return c->mdl_stereo_centre_value(c->top_back(), c->top_front(), c->left_down());

  return c->mdl_stereo_centre_value();
}

const IWString &
MDL_File_Supporting_Material::digits2 (int n)
{
  if (NULL == _digits2)
    _initialise_digits();

  return _digits2[n];
}

const IWString &
MDL_File_Supporting_Material::digits3 (int n)
{
  if (NULL == _digits2)
    _initialise_digits();

  return _digits3[n];
}

void
MDL_File_Supporting_Material::_initialise_digits ()
{
  _digits2 = new IWString[100];
  _digits3 = new IWString[1000];

  assert (NULL != _digits3);

  for (int i = 0; i <= 999; i++)
  {
    IWString & d = _digits3[i];

    if (i >= 100)
      d << i;
    else if (i >= 10)
      d << ' ' << i;
    else
      d << "  " << i;
  }

  for (int i = 0; i <= 99; i++)
  {
    IWString & d = _digits2[i];
    if (i >= 10)
      d << i;
    else
      d << ' ' << i;
  }

  return;
}

/*
  One column of the SD file is the isotopic specification. The standard says
  only write as a number if it is in the range -3 +4. We don't really do it
  that way
*/

void
MDL_File_Supporting_Material::append_isotope_information (IWString & output_buffer,
                            int iso,
                            int normal_isotope) const
{
  int to_be_written;

  if (_write_isotopes_as_numbers_rather_than_differences_from_normal)
    to_be_written = iso;
  else
    to_be_written = iso - normal_isotope;

  if (to_be_written > 99)
    output_buffer += _digits2[0];
  else if (to_be_written >= 0)
    output_buffer += _digits2[to_be_written];
  else if (to_be_written > -10)
    output_buffer << to_be_written;      // just write it
  else       // will be done in the M ISO records
    output_buffer += _digits2[0];

  return;
}

/*
  Once reactions needed to be able to write part of an mdl file record,
  this becomes a stand-alone function
*/

int
Molecule::_write_mdl_atom_record_element_charge_and_chirality (atom_number_t i,
                                                   IWString & output_buffer,
                                                   MDL_File_Supporting_Material & mdlfos) const
{
  const Atom * a = _things[i];

  const Element * e = a->element();

  if (! e->is_in_periodic_table() && 'R' == e->symbol()[0] && 2 == e->symbol().length() && isdigit(e->symbol()[1]))
  {
    if (mdlfos.write_Rn_groups_as_element())
    {
      output_buffer = e->symbol();
      output_buffer << ' ';
    }
    else
      output_buffer = "R# ";
  }
  else if (mdlfos.mdl_write_aromatic_atoms() && NULL != _aromaticity && IS_AROMATIC_ATOM(_aromaticity[i]))
  {
    if (0 == e->aromatic_symbol().length())
      output_buffer = e->symbol();
    else
      output_buffer = e->aromatic_symbol();
    (void) output_buffer.extend(3, ' ');
  }
  else
  {
    output_buffer = e->symbol();
    (void) output_buffer.extend(3, ' ');
  }

  int iso = a->isotope();
  if (0 == iso)
    output_buffer += mdlfos.digits2(0);
  else
    mdlfos.append_isotope_information(output_buffer, iso, e->normal_isotope());
      
  output_buffer += mdlfos.digits3(convert_to_mdl_charge(a->formal_charge()));

  int s;
  if (mdlfos.include_chiral_info_in_mdl_outputs() && _chiral_centres.number_elements())
    s = _mdl_atom_stereo_value(i, mdlfos);
  else
    s = 0;

  if (s)
    output_buffer += mdlfos.digits3(s);
  else if (mdlfos.isis_standard_records())
    output_buffer += mdlfos.digits3(0);

  return 1;
}

/*
  This function has only been used for mdl files.
  Discern cis-trans centres from examining the depictions.

  The bond is

     a1          a5
       \        /
        a3 -- a4
       /       \
     a2         a6

  In order to get consistent placement, we need to "grow" cis-trans perceptions once we get one
*/

int
Molecule::discern_cis_trans_bonds_from_depiction()
{
  if (_number_elements < 4)
    return 1;

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

 (void) compute_aromaticity_if_needed();   // path_scoring needs aromaticity

  int nb = _bond_list.number_elements();

  int rc = 0;
  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->is_double_bond())
      continue;

    atom_number_t a3 = b->a1();
    atom_number_t a4 = b->a2();

    int acon = _things[a3]->ncon();
    if (1 == acon || acon > 3)
      continue;

    acon = _things[a4]->ncon();
    if (1 == acon || acon > 3)
      continue;

    if (6 != _things[a3]->atomic_number())    // for compatability with Daylight, we only do Carbon atoms
      continue;

    if (6 != _things[a4]->atomic_number())
      continue;

    if (mdlfos->unspecified_double_bond_atoms().contains(a3))
      continue;

    if (is_ring_atom(a3) && is_ring_atom(a4))   // if a3 and a4 are both ring atoms, ignore
      continue;

    if (_discern_cis_trans_bond_from_depiction(b))
      rc++;
  }

  return rc;
}

/*
  We need special treatment for the case where there are multiple
  wedge bonds to the same atom.
*/

int
Molecule::_multiple_wedge_bonds_to_atom (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_wedge_definitive())
      rc++;
  }

  return rc > 1;
}

//#define DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS

int
Molecule::discern_chirality_from_wedge_bonds()
{
#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Discerning chirality from wedge bonds if present\n";
#endif

  _chiral_centres.resize_keep_storage(0);

  int * wedge_bonds_at_atom = new_int(_number_elements); std::unique_ptr<int[]> free_wedge_bonds_at_atom(wedge_bonds_at_atom);

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_wedge_definitive())
      continue;

    wedge_bonds_at_atom[b->a1()]++;
    wedge_bonds_at_atom[b->a2()]++;
  }

  int rc = 1;     // success until found otherwise

  for (int i = 0; i < _number_elements; ++i)
  {
    if (0 == wedge_bonds_at_atom[i])
      continue;

    if (_things[i]->ncon() < 2)    // singly connected atoms cannot be chiral centres
    {
      wedge_bonds_at_atom[i] = 0;
      continue;
    }

    if (wedge_bonds_at_atom[i] < 2)
      continue;

    if (_discern_chirality_from_multiple_wedge_bonds(i))
      continue;

    const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

    if (mdlfos->mdl_display_invalid_chiral_connectivity())
      cerr << "Molecule::_discern_chirality_from_wedge_bonds: cannot determine chirality multiple wedge bonds at atom " << i << "\n";

    rc = 0;
  }

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_wedge_definitive())
      continue;

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    if (wedge_bonds_at_atom[a1] > 1 || wedge_bonds_at_atom[a2] > 1)    // processed above
      continue;

    int tmprc;

    if (b->is_wedge_up())
      tmprc = _discern_chirality_from_wedge_bond(a1, a2, 1);
    else if (b->is_wedge_down())
      tmprc = _discern_chirality_from_wedge_bond(a1, a2, -1);
    else if (b->is_wedge_any())
      tmprc = _create_unspecified_chirality_object(a1); 
    else
      continue;

    if (0 == tmprc)
    {
      const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

      if (mdlfos->mdl_display_invalid_chiral_connectivity())
      {
        cerr << "Molecule::_discern_chirality_from_wedge_bonds: cannot determine chirality\n";
        cerr << "Atoms " << b->a1() << " and " << b->a2() << ", molecule '" << _molecule_name << "'\n";
      }
      rc = 0;
    }
  }

  return rc;
}

/*
  If there are 4 atoms connected, we need to determine the directionality from
  the 3 other atoms - those not at the edge of the wedge bond.

  If there are 3 atoms connected, the other connection must be an implicit
  Hydrogen, and the directionality is determined by looking at the 3 connected
  atoms, including the one with the wedge bond
*/

int
Molecule::_discern_chirality_from_wedge_bond (atom_number_t a1, 
                          atom_number_t a2,
                          int direction)
{
  const Atom * aa1 = _things[a1];

#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Bond from atom " << a1 << ' ' << aa1->atomic_symbol() << " to " << a2 << " " << _things[a2]->atomic_symbol() << " direction " << direction << endl;
  cerr << "ncon = " << aa1->ncon() << endl;
#endif

  if (4 == aa1->ncon())
    return _discern_chirality_from_wedge_bond_4(a1, a2, direction);

// We assume that the Hydrogen or lone lair is on the opposite side of the
// page from the atom at the end of the wedge bond, so we reverse the direction

  if (3 == aa1->ncon())
    return _discern_chirality_from_wedge_bond_4(a1, INVALID_ATOM_NUMBER, -direction);

// Connectivity is low, maybe we are just at the wrong end of the bond

  const Atom * aa2 = _things[a2];

  if (4 == aa2->ncon())
    return _discern_chirality_from_wedge_bond_4(a2, a1, -direction);

  if (3 == aa2->ncon())
    return _discern_chirality_from_wedge_bond_4(a2, INVALID_ATOM_NUMBER, direction);

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  if (mdlfos->mdl_display_invalid_chiral_connectivity())
    cerr << "Molecule::_discern_chirality_from_wedge_bond: base atom " << aa1->atomic_symbol() << " has " << aa1->ncon() << " connections!\n";

  return return_code_depending_on_ignore_incorrect_chiral_input();
}

/*
  Atom ZATOM is 4 connected and has a wedge bond to atom A2
*/

int
Molecule::_discern_chirality_from_wedge_bond_4 (atom_number_t zatom,
                          atom_number_t a2,
                          int direction)
{
  const Atom * centre = _things[zatom];

// The other 3 atoms connected

  Coordinates c[3];
  atom_number_t a[3];

  int nfound = 0;
  for (int i = 0; i < centre->ncon(); i++)
  {
    atom_number_t j = centre->other(zatom, i);

    if (a2 == j)
      continue;

    a[nfound] = j;

    c[nfound] = *(_things[j]) - *centre;
    c[nfound].normalise();

    nfound++;
  }

  assert (3 == nfound);

// I observed that a positive dot product means anti-clockwise rotation from c[0]

  angle_t theta1 = c[0].angle_between_unit_vectors(c[1]);
  angle_t theta2 = c[0].angle_between_unit_vectors(c[2]);

  if (0.0 == theta1 && 0.0 == theta2)
  {
    cerr << "Molecule::_discern_chirality_from_wedge_bond_4: two zero angles encountered '" << _molecule_name << "'\n";
    return 1;     // ignore the problem
  }

  Coordinates x01(c[0]);
  x01.cross_product(c[1]);
  Coordinates x12(c[1]);
  x12.cross_product(c[2]);
  Coordinates x20(c[2]);
  x20.cross_product(c[0]);

  coord_t z01 = x01.z();
  coord_t z12 = x12.z();
  coord_t z20 = x20.z();

/*
  OK, there are bugs in this, but I don't have the time to chase them down. 
  If there are multiple wedge bonds to an atom, we will get multiple chiral
  centres. This example shows a case where the two chiral centre objects are
  incompatible - file was called t9b.mol


  -ISIS-  10300013162D

  5  4  0  0  0  0  0  0  0  0999 V2000
    1.4417   -3.2708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8292   -7.6583    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    5.8000  -13.8833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.1583   -1.6583    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   11.2000   -4.5250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  1  0  0  0
  1  2  1  0  0  0  0
  2  4  1  6  0  0  0
  2  5  1  0  0  0  0
M  END

  resolve this sometime if it ever matters.
*/

  int rotation;     // clockwise is positive - going 0 -> 1 -> 2

  if (z01 >= 0.0 && z12 >= 0.0)
  {
    rotation = 1;
  }
  else if (z01 < 0.0 && z12 < 0.0)
  {
    rotation = -1;
  }
  else if (z12 >= 0.0 && z20 >= 0.0)
  {
    rotation = 1;
  }
  else if (z12 < 0.0 && z20 < 0.0)
  {
    rotation = -1;
  }
  else if (z20 >= 0.0 && z01 >= 0.0)
  {
    rotation = 1;
  }
  else if (z20 < 0.0 && z01 < 0.0)
  {
    rotation = -1;
  }
  else
  {
    cerr << "Molecule::_discern_chirality_from_wedge_bond_4: unusual geometry z01 = " << z01 << " z12 = " << z12 << " z20 = " << z20 << endl;
    assert (NULL == "this should not happen");
    return 0;
  }

#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  for (int i = 0; i < 3; i++)
  {
    cerr << i << " is atom " << a[i] << ' ' << _things[a[i]]->atomic_symbol() << endl;
  }
  cerr << "z01 = " << z01 << " z12 = " << z12 << " z20 = " << z20 << " rotation " << rotation << endl;
#endif

  if (rotation < 0)
    return _create_chiral_centre(zatom, a[0], a[1], a[2], a2, direction);

  return _create_chiral_centre(zatom, a[0], a[2], a[1], a2, direction);
}

int
Molecule::_create_chiral_centre (atom_number_t zatom,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 atom_number_t a3,
                                 atom_number_t a4,
                                 int direction)
{
#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Molecule::_create_chiral_centre: " << zatom << " a1 = " << a1 << " a2 = " << a2 << " a3 = " << a3 << " a3 = " << a3 << " a4 = " << a4 << " direction " << direction << endl;
#endif

  Chiral_Centre * c = new Chiral_Centre(zatom);
  c->set_top_front(a1);
  c->set_chirality_known(1);

  if (INVALID_ATOM_NUMBER == a4)
  {
    int lp;

    if (1 == _things[zatom]->implicit_hydrogens())
      c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
    else if (_things[zatom]->lone_pair_count(lp) && 1 == lp)
      c->set_top_back(CHIRAL_CONNECTION_IS_LONE_PAIR);
    else
    {
      delete c;

      const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

      if (mdlfos->mdl_display_invalid_chiral_connectivity())
        cerr << "Molecule::_create_chiral_centre: atom " << zatom << " strange chemistry, '" << _molecule_name << "'\n";
      return return_code_depending_on_ignore_incorrect_chiral_input();
    }
  }
  else
    c->set_top_back(a4);

  if (direction > 0)
  {
    c->set_left_down(a2);
    c->set_right_down(a3);
  }
  else
  {
    c->set_left_down(a3);
    c->set_right_down(a2);
  }

  _chiral_centres.add(c);

  return 1;
}

int
Molecule::_create_unspecified_chirality_object (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  Chiral_Centre * c = new Chiral_Centre(zatom);

  if (3 == a->ncon())
    c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);

  for (int i = 0; i < a->ncon(); i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (0 == i)
      c->set_top_front(j);
    else if (1 == i)
      c->set_left_down(j);
    else if (2 == i)
      c->set_right_down(j);
    else 
      c->set_top_back(j);
  }

  c->set_chirality_known(0);

  _chiral_centres.add(c);

  return 1;
}

/*
  We have the case where there are multiple wedge bonds going to a given atom.
*/

int
Molecule::_discern_chirality_from_multiple_wedge_bonds (atom_number_t zatom)
{
//cerr << "Molecule::_discern_chirality_from_multiple_wedge_bonds:processing atom " << zatom << endl;

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  const Bond * wb1 = NULL;
  const Bond * wb2 = NULL;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_wedge_definitive())
      continue;

    if (NULL == wb1)
      wb1 = b;
    else if (NULL == wb2)
      wb2 = b;
    else
    {
      const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

      if (mdlfos->mdl_display_invalid_chiral_connectivity())
        cerr << "Molecule::_discern_chirality_from_multiple_wedge_bonds:warning, too many wedge bonds to atom " << zatom << ", ignored\n";
    }
  }

  assert (NULL != wb2);

// Once upon a time I had a check here trying to figure out if the directionalities are
// consistent, but I'm not sure there are any inconsistencies possible with just two
// bonds being processed - see above

  if (wb1->is_wedge_up())
    return _discern_chirality_from_wedge_bond(wb1->a1(), wb1->a2(), 1);
  else
    return _discern_chirality_from_wedge_bond(wb1->a1(), wb1->a2(), -1);
}

int
Molecule::_assign_strange_atomic_symbol_to_atom (atom_number_t zatom,
                                                 const_IWSubstring s)
{
  const Element * e = get_element_from_symbol_no_case_conversion(s);

  if (NULL != e)
    ;
  else if (! auto_create_new_elements())
  {
    cerr << "Molecule::_assign_strange_atomic_symbol_to_atom:autocreate elements needed for R# and G group\n";
    return 0;
  }
  else
    e = create_element_with_symbol(s);

  if (NULL == e)
  {
    cerr << "Molecule::_assign_strange_atomic_symbol_to_atom:cannot create element '" << s << "'\n";
    return 0;
  }

// cerr << "Created element '" << e->symbol() << "'\n";

  _things[zatom]->set_element(e);

  return 1;
}

int
Molecule::_parse_M_RGP_record (const const_IWSubstring & buffer)
{
  Aprop atom_properties[MAX_PAIRS];

  int npairs;
  if (! fill_atom_property_array(buffer, npairs, atom_properties))
    return 0;

  for (int i = 0; i < npairs; i++)
  {
    const Aprop & api = atom_properties[i];

    atom_number_t a = api._atom_number - 1;
    int rnumber = api._property;

    IWString tmp;
    tmp << "R" << rnumber;

    if (! _assign_strange_atomic_symbol_to_atom(a, tmp))
      return 0;
  }

  _set_modified();

  return 1;
}

int
Molecule::_write_M_RGP_records (const MDL_File_Supporting_Material & mdlfos, std::ostream & os) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Element * e = _things[i]->element();

    if (e->is_in_periodic_table())
      continue;

    const IWString & s = e->symbol();

    if (2 == s.length() && 'R' == s[0] && isdigit(s[1]))
    {
      if (mdlfos.write_Rn_groups_as_element())
        continue;
    }
    else
      continue;

    os << "M  RGP  1 " << std::setw(3) << (i + 1) << "   " << s[1] << newline_string();
  }

  return os.good();
}

int
Molecule::_set_elements_based_on_atom_aliases (const resizable_array_p<Atom_Alias> & a)
{
  set_atomic_symbols_can_have_arbitrary_length(1);
  set_auto_create_new_elements(1);

  for (int i = 0; i < a.number_elements(); i++)
  {
    const Atom_Alias * ai = a[i];

    atom_number_t j = ai->atom_number();

    assert (OK_ATOM_NUMBER (this, j));

    const IWString alias = ai->alias();

    const Element * e = get_element_from_symbol_no_case_conversion(alias);

    if (NULL == e)
      e = create_element_with_symbol(alias);

    assert (NULL != e);

    _things[j]->set_element(e);
  }

  return 1;
}

int
Molecule::_process_mdl_g_record (const IWString & g,
                                 const const_IWSubstring & buffer)
{
  if (0 == buffer.length())
  {
    cerr << "Molecule::_process_mdl_g_record:empty G record\n";
    return 0;
  }

  assert (g.starts_with("G "));

  if (3 != g.nwords())
  {
    cerr << "Molecule::_process_mdl_g_record:Invalid G record\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;
  g.nextword(token, i);

  g.nextword(token, i);

  atom_number_t a;
  if (! token.numeric_value(a) || a < 1 || a > _number_elements)
  {
    cerr << "Molecule::_process_mdl_g_record:invalid atom number in G record\n";
    return 0;
  }

  return _assign_strange_atomic_symbol_to_atom(a - 1, buffer);
}

int
MDL_File_Supporting_Material::set_sdf_identifier (const const_IWSubstring & sdfid)
{
  IWString mysdfid = sdfid;

  if (mysdfid.starts_with('^'))
    mysdfid.remove_item(0);

  IWString tmp = "^>.*<";
  tmp += mysdfid;

  mysdfid = tmp;

  if (mysdfid.ends_with('$'))
    mysdfid.chop();

  mysdfid += ">";

  if (! _sdf_identifier.set_pattern(mysdfid))
  {
    cerr << "Cannot set sdfid pattern to '" << mysdfid << "'\n";
    return 0;
  }

//cerr << "Pattern set to '" << mysdfid << "'\n";

  return 1;
}

void
MDL_File_Supporting_Material::set_mdl_name_in_m_tag(const const_IWSubstring & s)
{
  if (_name_in_m_tag.length() > 0)
    _name_in_m_tag.resize_keep_storage(0);

  _name_in_m_tag << ' ' << s << ' ';

  return;
}

#ifdef FILE_SCOPE_ACCUMULATION
static int accumulate_mdl_chirality_features = 0;

static Set_of_Atoms unspecified_chiral_atoms_last_molecule;

static Set_of_Atoms up_bonds_last_molecule, down_bonds_last_molecule, squiggle_bonds_last_molecule;

static Set_of_Atoms unspecified_double_bond_atoms;

void
set_mdl_accumulate_mdl_chirality_features (int s)
{
  accumulate_mdl_chirality_features = s;
}

const Set_of_Atoms &
mdl_unspecified_chiral_atoms()
{
  return unspecified_chiral_atoms_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_up_bonds()
{
  return up_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_squiggle_bonds()
{
  return squiggle_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_down_bonds()
{
  return down_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_unspecified_double_bond_atoms()
{
  return unspecified_double_bond_atoms;
}
#endif

void
MDL_File_Supporting_Material::set_mdl_change_long_symbols_to (const const_IWSubstring & s)
{
  assert(s.length() > 0 && s.length() <= 2);

  _change_long_symbols_to = s;

  const Element * e = get_element_from_symbol_no_case_conversion(s);

  if (NULL == e)
    e = create_element_with_symbol(s);

  return;
}

static MDL_File_Supporting_Material default_mdl_file_supporting_material;

/*
  This parses a string into int's, with each field being three digits long

  This would be trivial with FORTRAN, what am I missing?
*/

int
int3d (const const_IWSubstring & buffer, int & i1, int & i2, int * i3)
{
  if (buffer.length() < 6)
    return 0;

  i1 = 0;
  int tmp = -1;     // this is outside the loop to make sure that we detect blank fields.
  for (int i = 0; i < 3; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 0;

    i1 = i1 * 10 + tmp;
  }

  if (-1 == tmp)
    return 0;

  tmp = -1;

  i2 = 0;
  for (int i = 3; i < 6; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 1;

    i2 = i2 * 10 + tmp;
  }

  if (-1 == tmp)
    return 1;

  if (NULL == i3)
    return 2;

  tmp = -1;

  *i3 = 0;
  for (int i = 6; i < 9; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 1;

    *i3 = *i3 * 10 + tmp;
  }

  if (-1 == tmp)
    return 2;

  return 3;
}

void
MDL_File_Supporting_Material::reset_for_next_molecule ()
{
  _unspecified_chiral_atoms_last_molecule.resize_keep_storage(0);
  _down_bonds_last_molecule.resize_keep_storage(0);
  _up_bonds_last_molecule.resize_keep_storage(0);
  _squiggle_bonds_last_molecule.resize_keep_storage(0);
  _unspecified_double_bond_atoms.resize_keep_storage(0);

  _a_records_found = 0;
  _g_records_found = 0;
  _aliases.resize(0);

  if (NULL == _digits2)
    _initialise_digits();

  return;
}

int
MDL_File_Supporting_Material::process_mdl_bond_translation (const const_IWSubstring & mdlbt)
{
  const_IWSubstring f, t;
  if (! mdlbt.split(f, '=', t))
  {
    cerr << "An mdl bond translation must be of the form 'number=number'\n";
    return 0;
  }

  int zfrom;
  if (! f.numeric_value(zfrom) || zfrom < 0)
  {
    cerr << "process_mdl_bond_translation:invalid from bond type specification '" << f << "'\n";
    return 0;
  }

  int zto;
  if (! t.numeric_value(zto) || zto < 0)
  {
    cerr << "process_mdl_bond_translation:invalid to bond type specification '" << t << "'\n";
    return 0;
  }

  return set_mdl_input_bond_type_translation(zfrom, zto);
}

template <typename T>
int
MDL_File_Supporting_Material::write_atoms_and_bonds (int na, int nb, T & os)
{
  if (NULL == _digits2)
    _initialise_digits();

  os << _digits3[na];
  os << _digits3[nb];

  return 1;
}

int
MDL_File_Supporting_Material::sdf_identifier_matches (const IWString & buffer)
{
  if (! _sdf_identifier.active())
    return 0;

  return _sdf_identifier.matches(buffer);
}

int
MDL_File_Supporting_Material::unspecified_chiral_atom_if_interested (int a)
{
  if (! _accumulate_mdl_chirality_features)
    return 0;

  _unspecified_chiral_atoms_last_molecule.add(a);

  return 1;
}

int
MDL_File_Supporting_Material::translate_input_bond_type (int bt) const
{
  if (NULL == _input_bond_type_translation_table)
    return bt;
  else
    return _input_bond_type_translation_table[bt];
}

int
Molecule::_contains_isotope_above (int s) const
{
  for (auto i = 0; i < _number_elements; ++i)
  {
    if (_things[i]->isotope() > s)
      return 1;
  }

  return 0;
}

template int Molecule::write_molecule_mdl<std::ostream>(std::ostream&, IWString const&) const;
//template int Molecule::write_molecule_mdl<IWString_and_File_Descriptor>(IWString_and_File_Descriptor &, IWString const&) const;

template int Molecule::read_molecule_mdl_ds<String_Data_Source>(String_Data_Source&, int);
template int Molecule::read_molecule_mdl_ds<iwstring_data_source>(iwstring_data_source&, int);
template int Molecule::write_extra_text_info<IWString_and_File_Descriptor>(IWString_and_File_Descriptor&) const;
template int Molecule::_write_m_iso_records<std::ostream>(std::ostream&, int) const;
template int Molecule::_write_m_chg_records<std::ostream>(std::ostream&, int) const;
template int Molecule::_mdl_write_atoms_and_bonds_record<std::ostream>(std::ostream&, int, int, MDL_File_Supporting_Material&) const;

