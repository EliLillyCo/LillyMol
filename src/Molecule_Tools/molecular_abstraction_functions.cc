#include <stdlib.h>

#include <limits>

#include "misc.h"

//#include "iwmalloc.h"

#include "molecule.h"
#include "rwsubstructure.h"
#include "target.h"
#include "path.h"
#include "iwmfingerprint.h"
#include "molecular_abstraction_functions.h"
#include "is_actually_chiral.h"

static int default_append_count_to_tag = 0;

void
set_append_count_to_tag (int s)
{
  default_append_count_to_tag = s;
}

static int default_write_only_if_changes = 0;

void
set_write_only_if_changes (int s)
{
  default_write_only_if_changes = s;
}

static int remove_invalid_chiral_centres_before_writing = 0;

void
set_remove_invalid_chiral_centres_before_writing (int s)
{
  remove_invalid_chiral_centres_before_writing = s;
}

static Molecule_With_Info_About_Parent empty_molecule;

static int write_empty_molecule_on_no_match = 0;

void
set_write_empty_molecule_on_no_match(const int s)
{
  write_empty_molecule_on_no_match = s;
}

Molecular_Abstraction_Base_Class::Molecular_Abstraction_Base_Class()
{
  _molecules_processed = 0;
  _molecules_changed = 0;

  _nbits = 2048;

  _isotope = 0;

  _append_count_to_tag = default_append_count_to_tag;

  _write_only_if_changes = default_write_only_if_changes;

   _min_atoms_needed_for_write = 0;
   _max_atoms_allowed_for_write = numeric_limits<int>::max();

   _min_atom_ratio_needed_for_write = 0.0f;
   _max_atom_ratio_allowed_for_write = 0.0f;

  return;
}

Molecular_Abstraction_Base_Class::~Molecular_Abstraction_Base_Class()
{
  return;
}

int
Molecular_Abstraction_Base_Class::ok() const
{
  if (_write_tag.length() && _fingerprint_tag.length())
  {
    cerr << "Molecular_Abstraction_Base_Class::ok:cannot have both smiles and FP output\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Base_Class::report(ostream & os) const
{
  os << "Abstraction processed " << _molecules_processed << " changed " << _molecules_changed << '\n';

  return 1;
}

template <typename T>
int
numeric_value_after_equals (const_IWSubstring token,   // note local copy
                            T & v,
                            T minval,
                            T maxval)
{
  token.remove_up_to_first('=');

  if (! token.numeric_value(v) || v < minval || v > maxval)
  {
    cerr << "Molecular_Abstraction_Base_Class::_process:invalid specifier'" << token << "'\n";
    return 0;
  }

  return 1;
}

#ifdef __GNUG__
template int numeric_value_after_equals(const_IWSubstring, int &, int, int);
template int numeric_value_after_equals(const_IWSubstring, float &, float, float);
#endif

int
Molecular_Abstraction_Base_Class::_process(const const_IWSubstring & token,
                                           const char * caller,
                                           int & fatal)
{
  if ("WRITEC" == token)
  {
    _write_tag = caller;
    _append_count_to_tag = 1;
  }
  else if ("WRITEIF" == token)
  {
    _write_tag = caller;
    _write_only_if_changes = 1;
  }
  else if (token.starts_with("WRITE_MIN_ATOMS="))
  {
    if (! numeric_value_after_equals(token, _min_atoms_needed_for_write, 1, numeric_limits<int>::max()))
    {
      cerr << "Molecular_Abstraction_Base_Class::_process:invalid WRITE_MIN_ATOMS directive '" << token << "'\n";
      return 0;
    }
    _write_tag = caller;
  }
  else if (token.starts_with("WRITE_MAX_ATOMS="))
  {
    if (! numeric_value_after_equals(token, _max_atoms_allowed_for_write, 1, numeric_limits<int>::max()))
    {
      cerr << "Molecular_Abstraction_Base_Class::_process:invalid WRITE_MAX_ATOMS directive '" << token << "'\n";
      return 0;
    }
    _write_tag = caller;
  }
  else if (token.starts_with("WRITE_MIN_PARENT_ATOM_RATIO="))
  {
    if (! numeric_value_after_equals(token, _min_atom_ratio_needed_for_write, 0.0f, 1.0f))
    {
      cerr << "Molecular_Abstraction_Base_Class::_process:invalid WRITE_MIN_PARENT_ATOM_RATIO specification '" << token << "'\n";
      return 0;
    }
    _write_tag = caller;
  }
  else if (token.starts_with("WRITE_MAX_PARENT_ATOM_RATIO="))
  {
    if (! numeric_value_after_equals(token, _max_atom_ratio_allowed_for_write, 0.0f, 1.0f))
    {
      cerr << "Molecular_Abstraction_Base_Class::_process:invalid WRITE_MAX_PARENT_ATOM_RATIO specification '" << token << "'\n";
      return 0;
    }
    _write_tag = caller;
  }
  else if (token.starts_with("WRITE"))
  {
    if (! _parse_write_directive(token, caller))
    {
      fatal = 1;
      return 0;
    }
  }
  else if (token.starts_with("FP"))
  {
    if (! _parse_fp_directive(token, caller))
    {
      fatal = 1;
      return 0;
    }
  }
  else if (token.starts_with("AT="))
  {
    if (! _parse_at_directive(token, caller))
    {
      fatal = 1;
      return 0;
    }
  }
  else if (token.starts_with("ISO"))
  {
    if (! _parse_isotope_directive(token, caller))
    {
      fatal = 1;
      return 0;
    }
  }
  else
  {
    fatal = 0;
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Base_Class::_parse_write_directive(const const_IWSubstring & token,
                                        const char * caller)
{
  if ("WRITE" == token)
  {
    _write_tag = caller;
    return 1;
  }

  if (token.starts_with("WRITE=") && token.length() > 6)
  {
    token.from_to(6, token.length() - 1, _write_tag);
    return 1;
  }

  cerr << "Molecular_Abstraction_Base_Class::_parse_write_directive:cannot process '" << token << "'\n";
  return 0;
}

int
Molecular_Abstraction_Base_Class::_parse_fp_directive(const const_IWSubstring & token,
                                                      const char * caller)
{
  assert(0 == _fingerprint_tag.length());

  if ("FP" == token)
  {
    _smiles_tag << caller << '<';
    _fingerprint_tag << "FP" << caller << '<';
    return 1;
  }

// Jul 2015. Want to be able to specify width for each fp. This is not
// compatible with the FP= directive below, but that's ok, don't think it was being used

  if (token.starts_with("FP.") && token.length() > 3)
  {
    const_IWSubstring tmp(token);
    tmp.remove_leading_chars(3);

    if (! tmp.numeric_value(_nbits) || _nbits < 8)
    {
      cerr << "Molecular_Abstraction_Base_Class::_parse_fp_directive:invalid bits directive '" << token << "'\n";
      return 0;
    }
    _smiles_tag << caller << '<';
    _fingerprint_tag << "FP" << caller << '<';

    return 1;
  }

  if (token.starts_with("FP=") && token.length() > 3)
  {
    const_IWSubstring tmp(token);
    tmp.remove_leading_chars(3);

    _smiles_tag << tmp;
    if (! _smiles_tag.ends_with('<'))
      _smiles_tag << '<';

    if (! tmp.starts_with("FP"))
      _fingerprint_tag << "FP";

    _fingerprint_tag << tmp;

    if (! _fingerprint_tag.ends_with('<'))
      _fingerprint_tag << '<';

    return 1;
  }

  cerr << "Molecular_Abstraction_Base_Class::_parse_fp_directive:cannot process '" << token << "'\n";
  return 0;
}

int
Molecular_Abstraction_Base_Class::_parse_at_directive (const const_IWSubstring & token,
                                                       const char * caller)
{
  if (token.starts_with("AT=") && token.length() > 3)
  {
    const_IWSubstring tmp(token);
    tmp.remove_leading_chars(3);

    if (! _atom_typing_specification.build(tmp))
    {
      cerr << "Molecular_Abstraction_Base_Class::_parse_at_directive:invalid atom typing specification '" << token << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "Molecular_Abstraction_Base_Class::_parse_at_directive:unrecognised directory '" << token << "'\n";
  return 0;
}

int
Molecular_Abstraction_Base_Class::_parse_isotope_directive(const const_IWSubstring & token,
                                                const char * caller)
{
  if (0 != _isotope)
    cerr << "Molecular_Abstraction_Base_Class::_parse_isotope_directive:duplicate isotope directives\n";

  if ("ISO" == token || "isotope" == token)
  {
    _isotope = 1;
    return 1;
  }

  if (token.starts_with("ISO=") && token.length() > 4)
  {
    const_IWSubstring tmp(token);
    tmp.remove_leading_chars(4);

    if (! tmp.numeric_value(_isotope) || _isotope < 0)
    {
      cerr << "Molecular_Abstraction_Base_Class::_parse_isotope_directive:invalid isotope '" << token << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "Molecular_Abstraction_Base_Class::_parse_isotope_directive:cannot process '" << token << "'\n";
  return 0;
}

/*
  We have begun a request, but it cannot process the molecle. What should we do?
*/

int
Molecular_Abstraction_Base_Class::_handle_no_match_to_query(Molecule_With_Info_About_Parent & m,
                                                            IWString_and_File_Descriptor & output)
{
  if (write_empty_molecule_on_no_match)
    return _do_any_writing_needed(empty_molecule, 0, output);

  return _do_any_writing_needed(m, 0, output);
}

int
Molecular_Abstraction_Base_Class::_do_any_writing_needed (Molecule_With_Info_About_Parent & m,
                                                int rc,
                                                IWString_and_File_Descriptor & output)
{
  if (0 == rc && _write_only_if_changes)
    return 1;

  const int matoms = m.natoms();

  if (matoms < _min_atoms_needed_for_write)
    return 1;

  if (matoms > _max_atoms_allowed_for_write)
    return 1;

  if (_min_atom_ratio_needed_for_write > 0.0f || _max_atom_ratio_allowed_for_write > 0.0f)
  {
    float r = static_cast<float>(matoms) / static_cast<float>(m.parent_natoms());

    if (r < _min_atom_ratio_needed_for_write)
      return 1;

    if (_max_atom_ratio_allowed_for_write > 0.0f && r > _max_atom_ratio_allowed_for_write)
      return 1;
  }

  if (remove_invalid_chiral_centres_before_writing)
    do_remove_invalid_chiral_centres(m);

//cerr << "_do_any_writing_needed, _fingerprint_tag " << _fingerprint_tag << " _atom_typing_specification " << _atom_typing_specification.active() << endl;

  if (_write_tag.length())
  {
    output << m.unique_smiles() << ' ' << m.name() << ' ' << _write_tag;
    if (_append_count_to_tag)
      output << ':' << rc;
    output << '\n';
  }
  else if (_fingerprint_tag.length())
  {
    output << _smiles_tag << m.unique_smiles() << ">\n";
    IWMFingerprint fp;

    if (_atom_typing_specification.active())
    {
      int * atype = new int[matoms]; unique_ptr<int[]> free_atype(atype);
      _atom_typing_specification.assign_atom_types(m, atype);
      fp.construct_fingerprint(m, atype, NULL);
    }
    else
      fp.construct_fingerprint(m);

    output << _fingerprint_tag;

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);

    output << tmp << ">\n";
  }

  return 1;
}

Molecular_Abstraction_Transform::Molecular_Abstraction_Transform()
{
  return;
}

Molecular_Abstraction_Change_Bond_Type::Molecular_Abstraction_Change_Bond_Type()
{
  _bt = INVALID_BOND_TYPE;

  return;
}

int
Molecular_Abstraction_Transform::debug_print(ostream & os) const
{
  os << "Molecular_Abstraction_Transform::with " << _eto.number_elements() << " components\n";

  if (_eto.number_elements() != _smarts.number_elements())
  {
    cerr << "Molecular_Abstraction_Transform::debug_print:invalid counts\n";
    return 0;
  }

  for (int i = 0; i < _eto.number_elements(); i++)
  {
    os << " i = " << i << " transform to " << _eto[i]->symbol() << '\n';
  }

  return 1;
}

/*
  Will be of the form smarts->ele or smarts->ele
*/

static int
divide_into_smarts_and_element(const const_IWSubstring & token,
                               const_IWSubstring & smarts,
                               const_IWSubstring & ele)
{
  int ndx = token.index("->");
  if (ndx >= 0)
  {
    if (0 == ndx)
      return 0;

    if (ndx == token.length() - 2)
      return 0;

    token.from_to(0, ndx - 1, smarts);
    token.from_to(ndx + 2, token.length() - 1, ele);

    return 1;
  }

  ndx = token.rindex('=');

  if (ndx >= 0)
  {
    if (0 == ndx || ndx == token.length() - 1)
      return 0;

    token.from_to(0, ndx - 1, smarts);
    token.from_to(ndx + 1, token.length() - 1, ele);

    return 1;
  }

  return 0;
}

int
Molecular_Abstraction_Transform::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "TRANSFORM", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Base_Class::build:cannot process '" << token << "'\n";
      return 0;
    }

//  must be of the form smarts=ele or smarts->ele

    const_IWSubstring smarts, ele;
    if (! divide_into_smarts_and_element(token, smarts, ele))
    {
      cerr << "Molecular_Abstraction_Base_Class::build:must be of form 'smarts->ele', '" << token << "' not allowed\n";
      return 0;
    }

    Single_Substructure_Query * q = new Single_Substructure_Query;

    if (! q->create_from_smarts(smarts))
    {
      cerr << "Molecular_Abstraction_Transform::build:invalid smarts '" << smarts << "'\n";
      delete q;
      return 0;
    }

    q->set_find_unique_embeddings_only(1);

    _smarts.add(q);

    const Element * e = get_element_from_symbol_no_case_conversion(ele);
    if (NULL == e)
      e = create_element_with_symbol(ele);

    _eto.add(e);
  }

  if (0 == _smarts.number_elements())
  {
    cerr << "Molecular_Abstraction_Transform::build:no substructure\n";
    return 0;
  }

  if (_smarts.number_elements() != _eto.number_elements())
  {
    cerr << "Molecular_Abstraction_Transform::build:mismatch\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Transform::process (Molecule_With_Info_About_Parent & m,
                                          IWString_and_File_Descriptor & output)
{
  assert (_eto.number_elements() > 0);
  assert (_eto.number_elements() == _smarts.number_elements());

  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  int n = _eto.number_elements();

  Molecule_to_Match target(&m);

  const Element ** new_element = new const Element *[matoms]; std::unique_ptr<const Element *[]> free_new_element(new_element);
  for (int i = 0; i < matoms; i++)
  {
    new_element[i] = NULL;
  }

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    Substructure_Results sresults;

    int nhits = _smarts[i]->substructure_search(&m, sresults);

    if (0 == nhits)
      continue;

    rc++;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * s = sresults.embedding(j);

      atom_number_t k = s->item(0);

      new_element[k] = _eto[i];
    }
  }

  if (rc)
  {
    _molecules_changed++;
    for (int i = 0; i < matoms; i++)
    {
      if (NULL != new_element[i])
        m.set_element(i, new_element[i]);
    }
  }

  return _do_any_writing_needed(m, rc, output);
}

Molecular_Abstraction_All_Transform::Molecular_Abstraction_All_Transform()
{
  _eto = NULL;

  return;
}

int
Molecular_Abstraction_All_Transform::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "ALLTRANS", fatal))
      continue;
    else if (fatal)
      return 0;

//  Must be an element

    _eto = get_element_from_symbol_no_case_conversion(token);
    if (NULL == _eto)
    {
      _eto = create_element_with_symbol(token);
      if (NULL == _eto)
      {
        cerr << "Molecular_Abstraction_Base_Class::build:invalid element '" << token << "'\n";
        return 0;
      }
    }

  }

  if (NULL == _eto)
  {
    cerr << "Molecular_Abstraction_Base_Class::build:no element, assuming Carbon\n";
    _eto = get_element_from_atomic_number(6);
  }

  return 1;
}

int
Molecular_Abstraction_All_Transform::process(Molecule_With_Info_About_Parent & m,
                                      IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _handle_no_match_to_query(m, output);

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (_eto != m.elementi(i))
    {
      m.set_element(i, _eto);
      rc++;
    }
  }

  if (rc)
    _molecules_changed++;

  return _do_any_writing_needed(m, rc, output);
}

Molecular_Abstraction_Remove_Atom::Molecular_Abstraction_Remove_Atom()
{
  _rejoin_all = 0;

  return;
}

int
Molecular_Abstraction_Remove_Atom::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RMAT", fatal))
      continue;
    else if (fatal)
      return 0;

    if ("rejoin" == token)
    {
      _rejoin_all = 1;
      continue;
    }

    if (token.starts_with("RJ="))
    {
      token.remove_leading_chars(3);
      int j;
      if (! token.numeric_value(j) || j < 1)
      {
        cerr << "Molecular_Abstraction_Remove_Atom:build:invalid rejoin '" << token << "'\n";
        return 0;
      }

      _rejoin[j] = 1;
      continue;
    }

//  token must be smarts

//  cerr << "Building smarts from '" <<token << "'\n";
    if (! _smarts.create_from_smarts(token))
    {
      cerr << "Molecular_Abstraction_Remove_Atom::build:invalid smarts '" << token << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);
  }

  if (_smarts.highest_initial_atom_number() < 0)
  {
    cerr << "Molecular_Abstraction_Remove_Atom:build:no smarts\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Remove_Atom::process(Molecule_With_Info_About_Parent & m,
                                      IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;

  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  int * to_remove = new_int(matoms); std::unique_ptr<int[]> free_to_remove(to_remove);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    atom_number_t j = e->item(0);

    to_remove[j] = 1;
  }

  if (0 != _isotope)
  {
    for (int i = 0; i < matoms; i++)
    {
      if (! to_remove[i])
        continue;

      const Atom * ai = m.atomi(i);

      int icon = ai->ncon();
      for (int j = 0; j < icon; j++)
      {
        atom_number_t k = ai->other(j, i);
        if (! to_remove[k])
          m.set_isotope(k, _isotope);
      }
    }
  }

  if (_rejoin_all || _rejoin.number_elements())
    _do_removals_with_rejoins(m, to_remove);
  else
  {
    for (int i = matoms - 1; i >= 0; i--)
    {
      if (to_remove[i])
        m.remove_atom(i);
    }
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, nhits, output);
}

int
Molecular_Abstraction_Remove_Atom::_do_removals_with_rejoins(Molecule_With_Info_About_Parent & m,
                                        int * to_remove) const
{
// First check to see if any of the atoms to which we'll make connections
// will be removed

  int matoms = m.natoms();

#ifdef DEBUG_REMOVE_ATOM_PROCESS
  cerr << "Initial molecule has " << matoms << " atoms\n";
#endif

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (! to_remove[i])
      continue;

#ifdef DEBUG_REMOVE_ATOM_PROCESS
    cerr << "Will remove atom " << i << " type '" << m.smarts_equivalent_for_atom(i) << "'\n";
#endif

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    if (1 == acon)
      continue;

#ifdef DEBUG_REMOVE_ATOM_PROCESS
    if (_rejoin_all)
      cerr << "Will rejoin all\n";
    else
      cerr << "Rejoin value for " << acon << " is " << _rejoin[acon] << endl;
#endif

    if (_rejoin_all)
      ;
    else if (0 == _rejoin[acon])   // not reconnecting atoms with this number of connections
    {
      to_remove[i] = 2;     // just remove, no rejoin
      continue;
    }

    Set_of_Atoms connected;
    a->connections(i, connected);

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = connected[j];
      if (to_remove[k])
        to_remove[i] = 2;     // just remove, no rejoin
    }

#ifdef DEBUG_REMOVE_ATOM_PROCESS
    cerr << "After examining to_remove " << to_remove[i] << endl;
#endif
  }

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (0 == to_remove[i])
      continue;

    if (2 == to_remove[i])
    {
      m.remove_atom(i);
      continue;
    }

    Set_of_Atoms connected;
    int n = m.connections(i, connected);

    for (int j = 0; j < n; j++)
    {
      atom_number_t k = connected[j];
      if (k > i)
        connected[j] = k - 1;
    }

    m.remove_atom(i);
    for (int j = 0; j < n; j++)
    {
      atom_number_t a1 = connected[j];
      for (int k = j + 1; k < n; k++)
      {
        atom_number_t a2 = connected[k];
        if (! m.are_bonded(a1, a2))
          m.add_bond(a1, a2, SINGLE_BOND);
      }
    }
  }

  return 1;
}

Molecular_Abstraction_Delete_Atoms::Molecular_Abstraction_Delete_Atoms()
{
  return;
}

int
Molecular_Abstraction_Delete_Atoms::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RMATS", fatal))
      continue;
    else if (fatal)
      return 0;

//  token must be smarts

//  cerr << "Building smarts from '" <<token << "'\n";
    if (! _smarts.create_from_smarts(token))
    {
      cerr << "Molecular_Abstraction_Delete_Atoms::build:invalid smarts '" << token << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);
  }

  if (_smarts.highest_initial_atom_number() < 0)
  {
    cerr << "Molecular_Abstraction_Delete_Atoms:build:no smarts\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Delete_Atoms::process(Molecule_With_Info_About_Parent & m,
                                      IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  if (0 == m.natoms())
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;

  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  int matoms = m.natoms();

  int * to_remove = new_int(matoms); std::unique_ptr<int[]> free_to_remove(to_remove);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    e->set_vector(to_remove, 1);
  }

  if (_isotope)
  {
    for (int i = 0; i < matoms; i++)
    {
      if (! to_remove[i])
        continue;

      const Atom * ai = m.atomi(i);

      int icon = ai->ncon();
      for (int j = 0; j < icon; j++)
      {
        atom_number_t k = ai->other(i, j);
        if (! to_remove[k])
          m.set_isotope(k, _isotope);
      }
    }
  }

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (to_remove[i])
      m.remove_atom(i);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, nhits, output);
}

Molecular_Abstraction_Scaffold::Molecular_Abstraction_Scaffold()
{
  _keep_first_ring_attachment = 0;

  return;
}


int
Molecular_Abstraction_Scaffold::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "SCAFFOLD", fatal))
      continue;
    else if (fatal)
    {
      return 0;
    }
    else if ("keepfirst" == token)
    {
      _keep_first_ring_attachment = 1;
    }
  }

  return 1;
}

static int
apply_isotopic_labels(Molecule_With_Info_About_Parent & m,
                      const int * in_subset,
                      int isotope)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (! in_subset[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    if (2 == acon)
      continue;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (! in_subset[k])
      {
        m.set_isotope(i, isotope);
        rc++;
        break;
      }
    }
  }

  return rc;
}

int
Molecular_Abstraction_Scaffold::process(Molecule_With_Info_About_Parent & m,
                                IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  if (0 == m.nrings())
    return _handle_no_match_to_query(m, output);

  int * in_scaffold = new_int(matoms); std::unique_ptr<int[]> free_in_scaffold(in_scaffold);

  Set_of_Atoms add_at_end;

  m.ring_membership(in_scaffold);

  for (int i = 0; i < matoms; i++)
  {
    if (in_scaffold[i] <= 0)    // start with ring atoms
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    if (acon < 3)       // no branches to spinach here
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (in_scaffold[k])   // already known
        continue;

      if (1 == m.ncon(k))
      {
        if (b->is_double_bond())    // always add doubly bonded to ring
          in_scaffold[k] = 1;
        else if (_keep_first_ring_attachment)
          add_at_end.add(k);

        continue;
      }

//    at this stage, the bond from I to K may be the start of a group of spinach

      if (! _is_spinach(m, in_scaffold, i, k))
        ;
      else if (_keep_first_ring_attachment)
        add_at_end.add(k);
    }
  }

  if (add_at_end.number_elements())
    add_at_end.set_vector(in_scaffold, 1);

  if (_isotope)
    apply_isotopic_labels(m, in_scaffold, _isotope);

  int rc = 0;
  for (int i = matoms - 1; i >= 0; i--)
  {
    if (! in_scaffold[i])
    {
      m.remove_atom(i);
      rc++;
    }
  }

  if (rc)
    _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}

/*
  We are really asking do we encounter terminal atoms down this path
*/

int
Molecular_Abstraction_Base_Class::_is_spinach(Molecule_With_Info_About_Parent & m,
                                            int * in_scaffold,
                                            atom_number_t aprev,
                                            atom_number_t zatom) const
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

//Set_of_Atoms doubly_bonded_single_connection;
  Set_of_Atoms doubly_bonded;
  Set_of_Atoms to_be_checked;
  int found_singly_connected = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (j == aprev)
      continue;

    if (in_scaffold[j])   // we have reached a ring atom
    {
      in_scaffold[zatom] = 1;
      continue;
    }

    if (b->is_double_bond())
      doubly_bonded.add(j);

    int jcon = m.ncon(j);

    if (1 == jcon)
    {
      if (b->is_double_bond())
        ;
      else
        found_singly_connected++;

      continue;
    }
    else
      to_be_checked.add(j);
  }

// First down any other bonds to be checked. If any of them get to a ring,
// we are part of the scaffold

  int ntbc = to_be_checked.number_elements();
  if (ntbc)
  {
    int ring_encountered = 0;

    for (int i = 0; i < ntbc; i++)
    {
      atom_number_t j = to_be_checked[i];
  
      if (! _is_spinach(m, in_scaffold, zatom, j))
        ring_encountered = 1;
    }

    if (ring_encountered)
      in_scaffold[zatom] = 1;
  }

  if (! in_scaffold[zatom])
    return 1;    // yes, we are spinach

// If we are part of the scaffold, then anything doubly bonded must be included

  for (int i = 0; i < doubly_bonded.number_elements(); i++)
  {
    atom_number_t d = doubly_bonded[i];
    in_scaffold[d] = 1;
  }

  return 0;    // no, we are not spinach
}

int
Molecular_Abstraction_Base_Class::_identify_scaffold (Molecule_With_Info_About_Parent & m, 
                                                      int * in_scaffold,
                                                      int keep_first_ring_attachment) const
{
  int matoms = m.natoms();

  Set_of_Atoms add_at_end;

  m.ring_membership(in_scaffold);

  for (int i = 0; i < matoms; i++)
  {
    if (in_scaffold[i] <= 0)    // start with ring atoms
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    if (acon < 3)       // no branches to spinach here
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (in_scaffold[k])   // already known
        continue;

      if (1 == m.ncon(k))
      {
        if (b->is_double_bond())    // always add doubly bonded to ring
          in_scaffold[k] = 1;
        else if (keep_first_ring_attachment)
          add_at_end.add(k);

        continue;
      }

//    at this stage, the bond from I to K may be the start of a group of spinach

      if (! _is_spinach(m, in_scaffold, i, k))
        ;
      else if (keep_first_ring_attachment)
        add_at_end.add(k);
    }
  }

  if (add_at_end.number_elements())
    add_at_end.set_vector(in_scaffold, 1);

  return 1;
}

int
Molecular_Abstraction_Change_Bond_Type::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "CBT", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Change_Bond_Type::build:cannot process '" << token << "'\n";
      return 0;
    }

//  must be of the form smarts=btype or smarts->btype

    const_IWSubstring smarts, bt;
    if (! divide_into_smarts_and_element(token, smarts, bt))
    {
      cerr << "Molecular_Abstraction_Change_Bond_Type::build:must be of form 'smarts->bt', '" << token << "' not allowed\n";
      return 0;
    }

    if (! _smarts.create_from_smarts(smarts))
    {
      cerr << "Molecular_Abstraction_Change_Bond_Type::build:invalid smarts '" << smarts << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);

    if ('-' == bt || "SINGLE" == bt || '1' == bt)
      _bt = SINGLE_BOND;
    else if ('=' == bt || "DOUBLE" == bt || '2' == bt)
      _bt = DOUBLE_BOND;
    else if ('#' == bt || "TRIPLE" == bt || '3' == bt)
      _bt = TRIPLE_BOND;
    else if ('.' == bt || "NONE" == bt)
      _bt = NOT_A_BOND;
    else
    {
      cerr << "Molecular_Abstraction_Change_Bond_Type::build:unrecognised bond type '" << bt << "'\n";
      return 0;
    }
  }

  if (_smarts.highest_initial_atom_number() < 0)
  {
    cerr << "Molecular_Abstraction_Change_Bond_Type::build:no substructure\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Change_Bond_Type::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  if (0 == m.natoms())
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;

  int nhits = _smarts.substructure_search(&m, sresults);
//cerr << nhits << " hits to change bond query\n";

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    atom_number_t a1 = e->item(0);
    atom_number_t a2 = e->item(1);

    if (NOT_A_BOND == _bt)
      m.remove_bond_between_atoms(a1, a2);
    else
      m.set_bond_type_between_atoms(a1, a2, _bt);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, nhits, output);
}

Molecular_Abstraction_Change_All_Bonds::Molecular_Abstraction_Change_All_Bonds()
{
  _bt = SINGLE_BOND;

  return;
}

int
Molecular_Abstraction_Change_All_Bonds::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "ALLBONDS", fatal))  // great
      continue;
    else if (fatal)
      return 0;

    if ('-' == token || "SINGLE" == token || "single" == token)
      _bt = SINGLE_BOND;
    else if ('=' == token || "DOUBLE" == token || "double" == token)
      _bt = DOUBLE_BOND;
    else if ('#' == token || "TRIPLE" == token || "triple" == token)
      _bt = TRIPLE_BOND;
    else
    {
      cerr << "Molecular_Abstraction_Change_All_Bonds::build:unrecognised bond type '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Change_All_Bonds::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  int nb = m.nedges();

  int rc = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (0 == (b->btype() & _bt))
    {
      m.set_bond_type_between_atoms(b->a1(), b->a2(), _bt);
      rc++;
    }
  }

  if (rc)
    _molecules_changed++;

  return _do_any_writing_needed(m, rc, output);
}

Molecular_Abstraction_Rings::Molecular_Abstraction_Rings()
{
  return;
}

int
Molecular_Abstraction_Rings::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RINGS", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Rings::build:cannot process '" << token << "'\n";
      return 0;
    }
    else
    {
      cerr << "Molecular_Abstraction_Rings::build:directive not recognised '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
add_doubly_bonded_singly_connected(Molecule_With_Info_About_Parent & m,
                                   int * in_system)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (in_system[i])
      continue;

    const Atom * ai = m.atomi(i);

    if (1 != ai->ncon())
      continue;

    const Bond * b = ai->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(i);

    if (! in_system[j])
      continue;

    if (! m.is_aromatic(j))
      continue;

//  in_system[i] = 2;
    in_system[i] = in_system[j];
    rc++;
  }

  return rc;
}


int
Molecular_Abstraction_Rings::process(Molecule_With_Info_About_Parent & m,
                                     IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  _molecules_changed++;

  if (0 == m.nrings())
  {
    m.resize(0);
    return _do_any_writing_needed(m, 0, output);
  }

  int * ring_membership = new_int(matoms); std::unique_ptr<int[]> free_ring_membership(ring_membership);

  m.ring_membership(ring_membership);

// include any =O and such groups

  add_doubly_bonded_singly_connected(m, ring_membership);

  if (_isotope)
    apply_isotopic_labels(m, ring_membership, _isotope);

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (0 == ring_membership[i])
      m.remove_atom(i);
  }

  return _do_any_writing_needed(m, 1, output);
}

Set_of_Molecular_Abstractions::Set_of_Molecular_Abstractions()
{
  _n = 0;
  _a = NULL;

  return;
}

Molecular_Abstraction_Largest_Ring_System::Molecular_Abstraction_Largest_Ring_System()
{
  _spiro = 0;

  return;
}

int
Molecular_Abstraction_Largest_Ring_System::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "BIGRING", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Largest_Ring_System::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if ("spiro" == token)
      _spiro = 1;
    else
    {
      cerr << "Molecular_Abstraction_Largest_Ring_System::build:directive not recognised '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Largest_Ring_System::process(Molecule_With_Info_About_Parent & m,
                                                 IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  const int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  int nr = m.nrings();

  if (0 == nr)
    return _handle_no_match_to_query(m, output);

  int * ring_system_membership = new_int(matoms); std::unique_ptr<int[]> free_ring_system_membership(ring_system_membership);

  if (_spiro)
    m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership);
  else
    m.label_atoms_by_ring_system(ring_system_membership);

  int * atoms_in_system = new_int(nr + 1); std::unique_ptr<int[]> free_atoms_in_system(atoms_in_system);

  for (int i = 0; i < matoms; i++)
  {
    if (ring_system_membership[i])
      atoms_in_system[ring_system_membership[i]]++;
  }

  int atoms_in_largest_ring_system = 0;
  int largest_ring_system = -1;

  for (int i = 1; i <= nr; i++)
  {
    if (atoms_in_system[i] > atoms_in_largest_ring_system)
    {
      atoms_in_largest_ring_system = atoms_in_system[i];
      largest_ring_system = i;
    }
  }

  for (int i = 0; i < matoms; i++)
  {
    if (ring_system_membership[i] != largest_ring_system)
      ring_system_membership[i] = 0;
  }

  add_doubly_bonded_singly_connected(m, ring_system_membership);

  if (_isotope)
    apply_isotopic_labels(m, ring_system_membership, _isotope);

  int rc = 0;

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (0 == ring_system_membership[i])
    {
      m.remove_atom(i);
      rc++;
    }
  }

  if (rc)
    _molecules_changed++;

  return _do_any_writing_needed(m, rc, output);
}

Molecular_Abstraction_Abstract_Ring_Form::Molecular_Abstraction_Abstract_Ring_Form()
{
  _arom_ele = NULL;
  _aliph_ele = NULL;
  _label_by_ring_size = 0;

  return;
}

int
Molecular_Abstraction_Abstract_Ring_Form::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "ARF", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Abstract_Ring_Form::build:cannot process '" << token << "'\n";
      return 0;
    }

    if ("lrs" == token)
    {
      _label_by_ring_size = 1;
      continue;
    }

    if ("lhc" == token)
    {
      _label_by_heteroatom_count = 1;
      continue;
    }

    const_IWSubstring t1, t2;
    if (! token.split(t1, '=', t2) || 0 == t1.length() || 0 == t2.length())
    {
      cerr << "Molecular_Abstraction_Abstract_Ring_Form:build:invalid element specification '" << token << "'\n";
      return 0;
    }

    const Element * e = get_element_from_symbol_no_case_conversion(t2);
    if (NULL == e)
      e = create_element_with_symbol(t2);

    if ("ELE" == t1)
    {
      _arom_ele = e;
      _aliph_ele = e;
    }
    else if ("AROM" == t1 || "AR" == t1)
      _arom_ele = e;
    else if ("ALIPH" == t1 || "AL" == t1)
      _aliph_ele = e;
    else
    {
      cerr << "Molecular_Abstraction_Abstract_Ring_Form::build:directive not recognised '" << token << "'\n";
      return 0;
    }
  }

  if (NULL == _arom_ele && NULL == _aliph_ele)
  {
//  cerr << "Molecular_Abstraction_Abstract_Ring_Form::build:no element\n";
    _arom_ele = get_element_from_symbol_no_case_conversion("Ar");
    _aliph_ele = get_element_from_symbol_no_case_conversion("Al");
  }
  else if (NULL != _arom_ele && NULL != _aliph_ele)
    ;
  else
  {
    cerr << "Molecular_Abstraction_Abstract_Ring_Form::buid:must specify both aliphatic and aromatic elements\n";
    return 0;
  }

  return 1;
}

static int
looks_like_biphenyl (const Molecule_With_Info_About_Parent & m,
                     const Ring & r1,
                     const Ring & r2)
{
  int n1 = r1.number_elements();

  for (int i = 0; i < n1; i++)
  {
    atom_number_t j = r1[i];

    const Atom * aj = m.atomi(j);

    int jcon = aj->ncon();

   if (2 == jcon)
      continue;

    for (int k = 0; k < jcon; k++)
    {
      const Bond * b = aj->item (k);

     if (b->nrings())
        continue;

      atom_number_t l = b->other(j);

     if (r2.contains (l))
        return 1;
    }
  }

  return 0;
}

/*
  when placing abstract rings, each ring will have a set of atoms
  to which it will be joined, as well as a set of other rings
  to which it will be joined
*/

class Joins_to_Abstract_Ring
{
  private:
    int _aromatic;
    int _ring_size;
    int _heteroatom_count;

    resizable_array_p<Connection> _to_atom;
    resizable_array_p<Connection> _to_ring;

  public:
    Joins_to_Abstract_Ring();

    void set_aromatic(int s) {_aromatic = s;}
    int is_aromatic() const { return _aromatic;}

    void set_ring_size(int s) { _ring_size = s;}
    int ring_size() const { return _ring_size;}

    void bonded_to_atom(atom_number_t, bond_type_t);
    void bonded_to_ring(int, bond_type_t);

    const resizable_array_p<Connection> & bonds_to_atoms() const { return _to_atom;}
    resizable_array_p<Connection> & bonds_to_atoms() { return _to_atom;}

    const resizable_array_p<Connection> & bonds_to_rings() const { return _to_ring;}

    int count_heteroatoms (const Molecule & , const Ring &);
    int heteroatoms () const { return _heteroatom_count;}
};

Joins_to_Abstract_Ring::Joins_to_Abstract_Ring()
{
  _aromatic = 0;
  _ring_size = 0;
  _heteroatom_count = 0;

  return;
}

void
Joins_to_Abstract_Ring::bonded_to_atom(atom_number_t a, bond_type_t bt)
{
  Connection * c = new Connection(a, bt);

  _to_atom.add(c);

  return;
}

void
Joins_to_Abstract_Ring::bonded_to_ring(int a, bond_type_t bt)
{
  Connection * c = new Connection(a, bt);

  _to_ring.add(c);

  return;
}

int
Joins_to_Abstract_Ring::count_heteroatoms (const Molecule & m,
                                           const Ring & r)
{
  int ring_size = r.number_elements();

  _heteroatom_count = 0;

  for (int i = 0; i < ring_size; i++)
  {
    if (6 != m.atomic_number(r[i]))
      _heteroatom_count++;
  }

  return _heteroatom_count;
}

static int
eliminate_duplicates(const resizable_array_p<Connection> & btai,
                     resizable_array_p<Connection> & btaj)
{
  int rc = 0;

  int ni = btai.number_elements();

  for (int i = 0; i < ni; i++)
  {
    const Connection * ci = btai[i];

    atom_number_t a1 = ci->a2();

    for (int j = 0; j < btaj.number_elements(); j++)
    {
      const Connection * cj = btaj[j];

      if (cj->a2() == a1)
      {
        btaj.remove_item(j);
        rc++;
        break;
      }
    }
  }

  return rc;
}

static int
remove_duplicate_attachments_to_fused_rings(Molecule_With_Info_About_Parent & m,
                                            int n,
                                            Joins_to_Abstract_Ring * jar)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const Joins_to_Abstract_Ring & jari = jar[i];

    const Ring * ri = m.ringi(i);

    const resizable_array_p<Connection> & btai = jari.bonds_to_atoms();

    for (int j = i + 1; j < n; j++)
    {
      Joins_to_Abstract_Ring & jarj = jar[j];

      const Ring * rj = m.ringi(j);

      if (! ri->is_fused_to(rj))
        continue;

      resizable_array_p<Connection> & btaj = jarj.bonds_to_atoms();

      rc += eliminate_duplicates(btai, btaj);
    }
  }

  return rc;
}

void
Molecular_Abstraction_Abstract_Ring_Form::_do_apply_isotopic_label (Molecule_With_Info_About_Parent & m,
                                        atom_number_t zatom,
                                        const Joins_to_Abstract_Ring & jar) const
{
  int iso;

  if (_label_by_ring_size)
    iso = jar.ring_size();
  else
    iso = 0;

  if (_label_by_heteroatom_count)
    iso = 10 * iso + jar.heteroatoms();

  m.set_isotope(zatom, iso);

  return;
}

int
Molecular_Abstraction_Abstract_Ring_Form::_place_abstract_rings (Molecule_With_Info_About_Parent & m,
                                            const int * ring_membership,
                                            Joins_to_Abstract_Ring * jar)
{
  const int matoms = m.natoms();

  const int nr = m.nrings();

  if (_label_by_heteroatom_count)
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);

      int rn = ri->ring_number();

      jar[rn].count_heteroatoms(m, *ri);
    }
  }

// Identify the atoms bonded to each ring

#ifdef DEBUG_PLACE_ABSTRACT_RINGS
  cerr << "NR = " << nr << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    assert (i == ri->ring_number());

    Joins_to_Abstract_Ring & jari = jar[i];

    int n = ri->number_elements();

    jari.set_ring_size(n);

    if (ri->is_aromatic())
      jari.set_aromatic(1);

#ifdef DEBUG_PLACE_ABSTRACT_RINGS
    cerr << "Ring " << i << " has " << n << " atoms\n";
#endif

    for (int j = 0; j < n; j++)
    {
      atom_number_t k = ri->item(j);

      const Atom * ak = m.atomi(k);

      int kcon = ak->ncon();

#ifdef DEBUG_PLACE_ABSTRACT_RINGS
      cerr << "Atom " << k << " with " << kcon << " connections\n";
#endif

      for (int l = 0; l < kcon; l++)
      {
        const Bond * b = ak->item(l);

        atom_number_t n = b->other(k);

        if (ring_membership[n])
          continue;

#ifdef DEBUG_PLACE_ABSTRACT_RINGS
        cerr << "Ring " << i << " bonded to atom " << n << endl;
#endif
        jari.bonded_to_atom(n, b->btype());
      }
    }
  }

#ifdef DEBUG_PLACE_ABSTRACT_RINGS
  cerr << "Finished identifying atoms attached to rings\n";
#endif

// If there are cases where the same atom joins to two rings,
// and the rings are fused, remove one of those linkages

  if (nr > 1)
    remove_duplicate_attachments_to_fused_rings(m, nr, jar);

  int * tmp = new int[matoms]; std::unique_ptr<int[]> free_tmp (tmp);

// Now identify the bonds between the rings

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    set_vector(tmp, matoms, 0);

    ri->set_vector (tmp, 1);

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi (j);

      int c = rj->count_members_set_in_array (tmp, 1);

      bond_type_t bt = NOT_A_BOND;

     if (1 == c)    // spiro fused, probably very complex
        bt = TRIPLE_BOND;
      else if (2 == c)    // single bond shared
        bt = DOUBLE_BOND;
      else if (c > 2)
        bt = TRIPLE_BOND;
      else if (looks_like_biphenyl (m, *ri, *rj))
        bt = SINGLE_BOND;
      else
        continue;

      int rni = ri->ring_number();
      int rnj = rj->ring_number();

      jar[rni].bonded_to_ring(rnj, bt);
    }
  }

  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i])
      m.remove_bonds_to_atom(i);
  }

  int initial_matoms = matoms;

// Add bonds between newly created rings and atoms

  for (int i = 0; i < nr; i++)
  {
//  cerr << "Start ring " << i << endl;
    int matoms = m.natoms();

    if (jar[i].is_aromatic())
      m.add(_arom_ele);
    else
      m.add(_aliph_ele);

    if (_label_by_ring_size || _label_by_heteroatom_count)
      _do_apply_isotopic_label(m, matoms, jar[i]);

    const resizable_array_p<Connection> & s = jar[i].bonds_to_atoms();

    int n = s.number_elements();
    for (int j = 0; j < n; j++)
    {
      const Connection * c = s[j];

      atom_number_t k = c->a2();

//    cerr << "Adding bond between " << k << " and " << matoms << endl;
      m.add_bond(k, matoms, c->btype());
    }
  }

// Add bonds between rings

  for (int i = 0; i < nr; i++)
  {
    const Joins_to_Abstract_Ring & jari = jar[i];

    const resizable_array_p<Connection> & rr = jari.bonds_to_rings();

    int n = rr.number_elements();

    for (int j = 0; j < n; j++)
    {
      const Connection * c = rr[j];

      atom_number_t k = c->a2();

      m.add_bond(initial_matoms + i, initial_matoms + k, c->btype());
    }
  }

  return 1;
}

static int
eliminate_duplicate_fusions(Set_of_Atoms & spiro_joined,
                            Set_of_Atoms & conn1,
                            Set_of_Atoms & conn2)
{
  int n = spiro_joined.number_elements();
  assert (n == conn1.number_elements());
  assert (n == conn2.number_elements());

  for (int i = 0; i < n; i++)
  {
    const atom_number_t si = spiro_joined[i];
    const atom_number_t c1 = conn1[i];
    const atom_number_t c2 = conn2[i];

    for (int j = i + 1; j < n; j++)
    {
      if (si != spiro_joined[j])
        continue;

      atom_number_t c1j = conn1[j];
      atom_number_t c2j = conn2[j];

      if (c1 == c1j && c2 == c2j)
        ;
      else if (c1 == c2j && c2 == c1j)
        ;
      else                   // not sure how this could happen
        continue;

//    cerr << "Removing item " << j << endl;
      spiro_joined.remove_item(j);
      conn1.remove_item(j);
      conn2.remove_item(j);
      n--;
      j--;
//    cerr << "Decremented n to " << n << endl;
    }
  }

  assert (n == spiro_joined.number_elements());

  return n;
}

static atom_number_t
identify_spiro_atom_and_attachments(Molecule_With_Info_About_Parent & m,
                                    const Ring & r,
                                    const int * v,
                                    Set_of_Atoms & spiro_fused,
                                    Set_of_Atoms & a1,
                                    Set_of_Atoms & a2)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    if (0 == v[j])
      continue;

    if (m.nrings(j) > 2)  
      continue;

    spiro_fused.add(j);

    if (i == n - 1)
      a1.add(r[0]);
    else
      a1.add(r[i + 1]);

    if (0 == i)
      a2.add(r.last_item());
    else
      a2.add(r[i - 1]);

    return 1;
  }

  return 0;
}

//#define DEBUG_SPREAD_APART_ANY_SPIRO_RINGS

static const Element * spiro_element = NULL;

int
Molecular_Abstraction_Abstract_Ring_Form::_spread_apart_any_spiro_rings(Molecule_With_Info_About_Parent & m,
                                        int * tmp) const
{
  int nr = m.nrings();

  Set_of_Atoms spiro_joined;
  Set_of_Atoms conn1, conn2;

  int matoms = m.natoms();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    set_vector(tmp, matoms, 0);
    ri->set_vector(tmp, 1);

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi(j);

      if (1 != rj->count_members_set_in_array(tmp, 1))
        continue;

//    cerr << "Spiro fusion between rings " << i << " and " << j << endl;

      identify_spiro_atom_and_attachments(m, *rj, tmp, spiro_joined, conn1, conn2);
    }
  }

  int n = spiro_joined.number_elements();

  if (0 == n)
    return 0;

// Ran into problems with 
// where two rings share the same spiro fusion

  if (n > 1)
    n = eliminate_duplicate_fusions(spiro_joined, conn1, conn2);

  if (NULL == spiro_element)
    spiro_element = get_element_from_symbol_no_case_conversion("Sg");

#ifdef DEBUG_SPREAD_APART_ANY_SPIRO_RINGS
  cerr << "Found " << n << " spiro fusions\n";
#endif

  for (int i = 0; i < n; i++)
  {
    atom_number_t s = spiro_joined[i];
    atom_number_t a1 = conn1[i];
    atom_number_t a2 = conn2[i];

#ifdef DEBUG_SPREAD_APART_ANY_SPIRO_RINGS
    cerr << "Spiro atom " << s << " joined to " << a1 << " and " << a2 << endl;
#endif

    const Element * e = m.elementi(s);

    int matoms = m.natoms();

//  Break the ring

    m.remove_bond_between_atoms(s, a1);
    m.remove_bond_between_atoms(s, a2);

//  Add another atom and join to the open ring

    m.add(e);    // duplicate the spiro fused atom
    m.add_bond(matoms, a1, SINGLE_BOND);
    m.add_bond(matoms, a2, SINGLE_BOND);

//  Insert the spiro element between the two rings

    m.add(spiro_element);
    m.add_bond(matoms + 1, s, SINGLE_BOND);
    m.add_bond(matoms + 1, matoms, SINGLE_BOND);

//  If we have adjacent spiro fusions, we need to update the atom
//  number with the newly added ring atom

    for (int j = i + 1; j < n; j++)
    {
      if (conn1[j] == s)
        conn1[j] = m.natoms() - 2;
      else if (conn2[j] == s)
        conn2[j] = m.natoms() - 2;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Abstract_Ring_Form::process(Molecule_With_Info_About_Parent & m,
                                                  IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  int nr = m.nrings();

  if (0 == nr)
    return _handle_no_match_to_query(m, output);

  if (nr > 1)
  {
    int * tmp = new_int(matoms); std::unique_ptr<int[]> free_tmp(tmp);

    _spread_apart_any_spiro_rings(m, tmp);

    matoms = m.natoms();
  }

  m.compute_aromaticity_if_needed();

// First task is to figure out which rings are joined to each other and how

  int * ring_membership = new_int(matoms); std::unique_ptr<int[]> free_ring_membership(ring_membership);

  m.ring_membership(ring_membership);

  Joins_to_Abstract_Ring * jar = new Joins_to_Abstract_Ring[nr]; std::unique_ptr<Joins_to_Abstract_Ring[]> free_jar(jar);

  int rc = _place_abstract_rings (m, ring_membership, jar);

  if (rc)
  {
    for (int i = matoms - 1; i >= 0; i--)   // matoms is as in the original molecule
    {
      if (ring_membership[i])
        m.remove_atom(i);
    }

    _molecules_changed++;
  }

  return _do_any_writing_needed(m, rc, output);;
}

Molecular_Abstraction_Replace_Linker::Molecular_Abstraction_Replace_Linker()
{
  _linker_atom = NULL;

  return;
}

int
Molecular_Abstraction_Replace_Linker::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RPLINK", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Replace_Linker::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("ELE="))
    {
      token.remove_leading_chars(4);
      _linker_atom = get_element_from_symbol_no_case_conversion(token);
      if (NULL == _linker_atom)
        _linker_atom = create_element_with_symbol(token);

      assert (NULL != _linker_atom);
    }
    else
    {
      _linker_atom = get_element_from_symbol_no_case_conversion(token);
      if (NULL == _linker_atom)
        _linker_atom = create_element_with_symbol(token);

      assert (NULL != _linker_atom);
    }
//  else
//  {
//    cerr << "Molecular_Abstraction_Replace_Linker::build:directive not recognised '" << token << "'\n";
//    return 0;
//  }
  }

  if (NULL == _linker_atom)
    _linker_atom = get_element_from_atomic_number(3);  //  Lithium

  return 1;
}

/*
  Identify singly connected atoms that are doubly bonded to the chain
*/

static int
doubly_bonded_to_chain_atom(Molecule_With_Info_About_Parent & m,
                            atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  if (1 != a->ncon())
    return 0;

  if (2 != a->nbonds())
    return 0;

  atom_number_t j = a->other(zatom, 0);

  if (m.is_ring_atom(j))
    return 0;

  return 1;
}

static int
identify_linker_section(Molecule_With_Info_About_Parent & m,
                        int * ring_or_already_done,
                        atom_number_t zatom, 
                        Set_of_Atoms & ring_attachments)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (2 == ring_or_already_done[j])   // already done
      continue;

    if (1 == ring_or_already_done[j])
    {
      if (! doubly_bonded_to_chain_atom(m, j))
        ring_attachments.add(j);
      continue;
    }
      
    ring_or_already_done[j] = 2;
    identify_linker_section(m, ring_or_already_done, j, ring_attachments);
  }

  return ring_attachments.number_elements();
}

int
Molecular_Abstraction_Replace_Linker::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed ++;

  int nr = m.nrings();

  if (nr < 2)
    return _handle_no_match_to_query(m, output);

  int matoms = m.natoms();

// Identify those atoms that will not be part of linker groups
// Spinach and ring atoms are NOT linkers.Everything else is.

  int * non_linker = new int[matoms]; std::unique_ptr<int[]> free_non_linker(non_linker);

  int number_ring_systems = m.label_atoms_by_ring_system_including_spiro_fused(non_linker);

  if (1 == number_ring_systems)
    return _handle_no_match_to_query(m, output);

  for (int i = 0; i < matoms; i++)
  {
    if (non_linker[i])
      non_linker[i] = 4;         // any number but 1
  }

  m.identify_spinach_preset(non_linker);   // spinach is NOT a linker

  for (int i = 0; i < matoms; i++)
  {
    if (non_linker[i])
      non_linker[i] = 1;
  }

  resizable_array_p<Set_of_Atoms> ring_atoms_at_edges_of_linker_section;

  for (int i = 0; i < matoms; i++)
  {
    if (non_linker[i])   // ring atom, spinach or already processed
      continue;

    Set_of_Atoms * ring_attachments = new Set_of_Atoms;

    non_linker[i] = 2;   // indicates has been processed

    identify_linker_section(m, non_linker, i, *ring_attachments);

    ring_atoms_at_edges_of_linker_section.add(ring_attachments);
  }

  int nae = ring_atoms_at_edges_of_linker_section.number_elements();

  if (0 == nae)
    return _handle_no_match_to_query(m, output);

  int initial_matoms = m.natoms();

  for (int i = 0; i < ring_atoms_at_edges_of_linker_section.number_elements(); i++)
  {
    const Set_of_Atoms & s = *(ring_atoms_at_edges_of_linker_section[i]);

    int ns = s.number_elements();

    if (ns < 2)
    {
      cerr << "HUh, ns too small " << ns << " '" << m.smiles() << "'\n";
      continue;
    }

    assert (ns > 1);

    int atom_number_of_newly_added_linker = m.natoms();

    m.add(_linker_atom);

    for (int j = 0; j < ns; j++)
    {
      atom_number_t k = s[j];

      m.add_bond(atom_number_of_newly_added_linker, k, SINGLE_BOND);
    }
  }

  for (int i = initial_matoms - 1; i >= 0; i--)
  {
    if (2 == non_linker[i] || 1 == m.ncon(i))
      m.remove_atom(i);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}

Molecular_Abstraction_Fragment::Molecular_Abstraction_Fragment()
{
  return;
}

int
Molecular_Abstraction_Fragment::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "FRAGMENT", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Fragment::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("SMARTS="))
    {
      token.remove_leading_chars(7);
      if (! _fragment_to_keep.create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Fragment::build:invalid smarts '" << token << "'\n";
        return 0;
      }

      _fragment_to_keep.set_find_unique_embeddings_only(1);
    }
    else if (token.starts_with("KEEP="))
    {
      if (_fragment_to_keep.highest_initial_atom_number() >= 0)
      {
        cerr << "Molecular_Abstraction_Fragment::build:keep smarts alread specified\n";
        return 0;
      }
      token.remove_leading_chars(5);
      if (! _fragment_to_keep.create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Fragment::build:invalid keep smarts '" << token << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("REMOVE="))
    {
      if (_fragment_to_remove.highest_initial_atom_number() >= 0)
      {
        cerr << "Molecular_Abstraction_Fragment::build:remove smarts alread specified\n";
        return 0;
      }
      token.remove_leading_chars(7);
      if (! _fragment_to_remove.create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Fragment::build:invalid remove smarts '" << token << "'\n";
        return 0;
      }
    }
    else if (_fragment_to_keep.highest_initial_atom_number() >= 0)
    {
      cerr << "Molecular_Abstraction_Fragment::build:smarts already specified\n";
      return 0;
    }
    else
    {
      cerr << "Molecular_Abstraction_Fragment::build:directive not recognised '" << token << "'\n";
      return 0;
    }
  }

  if (_fragment_to_keep.highest_initial_atom_number() >= 0 &&
      _fragment_to_remove.highest_initial_atom_number() >= 0)
  {
    cerr << "Molecular_Abstraction_Fragment::build:cannot have both KEEP= and REMOVE= directives\n";
    return 0;
  }

  return 1;
}


int
Molecular_Abstraction_Fragment::_identify_fragments_hit (Molecule_With_Info_About_Parent & m,
                                        Single_Substructure_Query & q,
                                        int * hits_in_fragment)
{
  Substructure_Results sresults;

  int nhits = q.substructure_search(&m, sresults);

  if (0 == nhits)
    return 0;

  int rc = 0;
  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    int n = e->number_elements();

    for (int j = 0; j < n; j++)
    {
      atom_number_t k = e->item(j);

      int f = m.fragment_membership(k);

      if (hits_in_fragment[f] > 0)
        continue;

      hits_in_fragment[f] = 1;
      rc++;
    }
  }

  return rc;
}

int
Molecular_Abstraction_Fragment::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed ++;

  if (0 == m.natoms())
    return _do_any_writing_needed(m, 0, output);

  const int nf = m.number_fragments();

  if (nf <= 1)
    return _handle_no_match_to_query(m, output);

  if (0 == _fragment_to_keep.root_atoms() && 0 == _fragment_to_remove.root_atoms())
  {
    m.reduce_to_largest_fragment();
    _molecules_processed++;
    return _do_any_writing_needed(m, 1, output);
  }

  int * fragments_hit = new_int(nf); std::unique_ptr<int[]> free_fragments_hit(fragments_hit);

  int nhits;
  if (_fragment_to_keep.highest_initial_atom_number() >= 0)
    nhits = _identify_fragments_hit(m, _fragment_to_keep, fragments_hit);
  else if (_fragment_to_remove.highest_initial_atom_number() >= 0)
    nhits = _identify_fragments_hit(m, _fragment_to_remove, fragments_hit);
  else                  // should not happen
    return 0;

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  resizable_array<int> fragments_to_delete;
  if (_fragment_to_keep.highest_initial_atom_number())
  {
    for (int i = 0; i < nf; i++)
    {
      if (0 == fragments_hit[i])
        fragments_to_delete.add(i);
    }
  }
  else   // removing fragments
  {
    for (int i = 0; i < nf; i++)
    {
      if (fragments_hit[i])
        fragments_to_delete.add(i);
    }
  }

  if (fragments_to_delete.number_elements() == m.number_fragments())
    return _do_any_writing_needed(empty_molecule, 0, output);

  m.delete_fragments(fragments_to_delete);

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}

Molecular_Abstraction_Place_Isotope::Molecular_Abstraction_Place_Isotope()
{
  return;
}

void
Molecular_Abstraction_Place_Isotope::_initialise_isotope_on_first_matched_atom()
{
  _matched_atom.add(0);
  if (Molecular_Abstraction_Base_Class::_isotope > 0)
    _isotope.add(Molecular_Abstraction_Base_Class::_isotope);
  else
    _isotope.add(1);

  return;
}

void
Molecular_Abstraction_Place_Charge::_initialise_charge_on_first_matched_atom()
{
  _matched_atom.add(0);
  _charge.add(0);

  return;
}

static IW_Regular_Expression number_equals_number("^[0-9]+=[-]*[0-9]+$");

int
Molecular_Abstraction_Place_Isotope::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

//cerr << "Place isotope parsing '" << s << "'\n";

  if (0 == s.nwords())
  {
    _smarts.create_from_smarts("*");
    _matched_atom.add(0);
    _isotope.add(0);
    
    return 1;
  }

  if (1 == s.nwords())    // if just one token, put isotope 1 on first matched atom
  {
    if (! _smarts.create_from_smarts(s) || _smarts.highest_initial_atom_number() < 0)
    {
      cerr << "Molecular_Abstraction_Place_Isotope::build:invalid smarts '" << s << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);

    _initialise_isotope_on_first_matched_atom();
    
    return 1;
  }

// Severe confusion between the base class _isotope variable, and our array of isotopes.

  int iso_specified_here = 0;

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    if (token.starts_with("ISO="))
      iso_specified_here = 1;

    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "ISOTOPE", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Place_Isotope::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("SMARTS="))
    {
      token.remove_leading_chars(7);
      if (! _smarts.create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Place_Isotope::build:invalid smarts '" << token << "'\n";
        return 0;
      }

      _smarts.set_find_unique_embeddings_only(1);
      continue;
    }

    if (number_equals_number.matches(token))
    {
      const_IWSubstring sndx, siso;
      token.split(sndx, '=', siso);

      int ndx, iso;
      if (! sndx.numeric_value(ndx) || ndx < 0)
      {
        cerr << "Molecular_Abstraction_Place_Isotope::buid:invalid index '" << token << "'\n";
        return 0;
      }

      if (! siso.numeric_value(iso) || iso < 0)
      {
        cerr << "Molecular_Abstraction_Place_Isotope::buid:invalid isotope '" << token << "'\n";
        return 0;
      }

      _matched_atom.add(ndx);
      _isotope.add(iso);
      continue;
    }
        
    cerr << "Molecular_Abstraction_Place_Isotope::build:unrecognised directive '" << token << "'\n";
    return 0;
  }

// They must have entered just an isotope.

  if (_smarts.highest_initial_atom_number() < 0)
  {
    _smarts.create_from_smarts("*");
    _matched_atom.add(0);

    return 1;
  }

  assert (_matched_atom.number_elements() == _isotope.number_elements());

  if (_isotope.number_elements() > 0)
    ;
  else if (iso_specified_here && 0 == _matched_atom.number_elements())
  {
    _isotope.add(Molecular_Abstraction_Base_Class::_isotope);
    _matched_atom.add(0);
  }
  else
  {
    _initialise_isotope_on_first_matched_atom();
    return 1;
  }

  for (int i = 0; i < _matched_atom.number_elements(); i++)
  {
    atom_number_t j = _matched_atom[i];
    if (j > _smarts.highest_initial_atom_number())
    {
      cerr << "Molecular_Abstraction_Place_Isotope::build:matched atom " << j << " is invalid " << _smarts.highest_initial_atom_number() << endl;
      return 0;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Place_Isotope::process (Molecule_With_Info_About_Parent & m,
                                              IWString_and_File_Descriptor & output)
{
  _molecules_processed ++;

  if (0 == m.natoms())
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;
  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    for (int j = 0; j < _matched_atom.number_elements(); j++)
    {
      int k = _matched_atom[j];

      atom_number_t l = e->item(k);

//    cerr << "Setting isotope, atom " << l << " iso " << _isotope[j] << endl;
      m.set_isotope(l, _isotope[j]);
    }
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}

Molecular_Abstraction_Place_Charge::Molecular_Abstraction_Place_Charge()
{
  return;
}

int
Molecular_Abstraction_Place_Charge::build(const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  if (1 == s.nwords())    // if just one token, put isotope 1 on first matched atom
  {
    if (! _smarts.create_from_smarts(s) || _smarts.highest_initial_atom_number() < 0)
    {
      cerr << "Molecular_Abstraction_Place_Charge::build:invalid smarts '" << s << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);

    _initialise_charge_on_first_matched_atom();
    
    return 1;
  }

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "CHARGE", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Place_Charge::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("SMARTS="))
    {
      token.remove_leading_chars(7);
      if (! _smarts.create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Place_Charge::build:invalid smarts '" << token << "'\n";
        return 0;
      }

      _smarts.set_find_unique_embeddings_only(1);
      continue;
    }
    else if ('*' == token)
    {
      _smarts.create_from_smarts("*");
      _smarts.set_find_unique_embeddings_only(1);
      continue;
    }

    if (number_equals_number.matches(token))
    {
      const_IWSubstring sndx, siso;
      token.split(sndx, '=', siso);

      int ndx, charge;
      if (! sndx.numeric_value(ndx) ||
          ndx < 0)
      {
        cerr << "Molecular_Abstraction_Place_Charge::buid:invalid index '" << token << "'\n";
        cerr << "Must be of the form 'matched_atom=charge_value'\n";
        return 0;
      }

      if (! siso.numeric_value(charge))
      {
        cerr << "Molecular_Abstraction_Place_Charge::buid:invalid charge '" << token << "'\n";
        return 0;
      }

      _matched_atom.add(ndx);
      _charge.add(charge);
      continue;
    }
        
    cerr << "Molecular_Abstraction_Place_Charge::build:unrecognised directive '" << token << "'\n";
    return 0;
  }

  if (_smarts.highest_initial_atom_number() < 0)
    _smarts.create_from_smarts("*");

  assert (_matched_atom.number_elements() == _charge.number_elements());

  if (0 == _charge.number_elements())
  {
    _initialise_charge_on_first_matched_atom();
    return 1;
  }

  for (int i = 0; i < _matched_atom.number_elements(); i++)
  {
    atom_number_t j = _matched_atom[i];
    if (j > _smarts.highest_initial_atom_number())
    {
      cerr << "Molecular_Abstraction_Place_Charge::build:matched atom " << j << " is invalid " << _smarts.highest_initial_atom_number() << endl;
      return 0;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Place_Charge::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed ++;

  if (0 == m.natoms())
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;
  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  int this_molecule_changed = 0;

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    for (int j = 0; j < _matched_atom.number_elements(); j++)
    {
      int k = _matched_atom[j];

      atom_number_t l = e->item(k);

      if (m.set_formal_charge_if_different(l, _charge[j]))
        this_molecule_changed++;
    }
  }

  if (this_molecule_changed)
    _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}

Molecular_Abstraction_Compress_Consecutive::Molecular_Abstraction_Compress_Consecutive()
{
  return;
}

int
Molecular_Abstraction_Compress_Consecutive::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "COMPRCONSEC", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Compress_Consecutive::build:cannot process '" << token << "'\n";
      return 0;
    }

    if (_smarts.highest_initial_atom_number() >= 0)
    {
      cerr << "Molecular_Abstraction_Compress_Consecutive::build:only one smarts allowed\n";
      return 0;
    }

    if (! _smarts.create_from_smarts(token))
    {
      cerr << "Molecular_Abstraction_Compress_Consecutive::build:invalid smarts '" << token << "'\n";
      return 0;
    }

    _smarts.set_find_unique_embeddings_only(1);
  }

  if(_smarts.highest_initial_atom_number() < 0)
    _smarts.create_from_smarts("[CD2R0H2]");

  return 1;
}

static int
array_gsub(int * v, 
           int n, 
           int vfrom,
           int vto)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    if (vfrom != v[i])
      continue;

    v[i] = vto;
    rc++;
  }

  return rc;
}

/*
  We want to compress the group, leaving just one

  A-X-X-X-X-X-A -> A-X-A
*/

//#define DEBUG_COMPRESS_GROUP

int
Molecular_Abstraction_Compress_Consecutive::_compress_group(Molecule_With_Info_About_Parent & m,
               int * to_remove,
               int flag, 
               resizable_array<const Bond *> & join_points) const
{
  assert (join_points.number_elements() > 1);

  const Bond * b1 = join_points[0];
  const Bond * b2 = join_points[1];

#ifdef DEBUG_COMPRESS_GROUP
  cerr << "B1 is " << b1->a1() << " to " << b1->a2() << endl;
  cerr << "B2 is " << b2->a1() << " to " << b2->a2() << endl;
#endif

  atom_number_t rhs;  // we will make a bond between RETAINED and RHS

  if (flag == to_remove[b1->a1()])
    rhs = b1->a2();
  else
    rhs = b1->a1();

// At the other end, identify the atom that will be retained

  atom_number_t retained;
  atom_number_t lhs;
  if (flag == to_remove[b2->a1()])
  {
    lhs = b2->a2();
    retained = b2->a1();
  }
  else
  {
    lhs = b2->a1();
    retained = b2->a2();
  }

// Deal with cases when compressing consecutive ring CH2 groups

  if (lhs == rhs)
  {
    lhs = b1->other(lhs);   // get other end of bond

    if (m.are_bonded(lhs, retained))    // starting with 3 membered ring
      return 0;

    m.add_bond(lhs, retained, SINGLE_BOND);    // form 3 membered ring
    to_remove[lhs] = 0;
    to_remove[retained] = 0;

    return 1;
  }

#ifdef DEBUG_COMPRESS_GROUP
  cerr << "Retained atom " << retained << ", making bond to " << rhs << endl;
#endif

  if (rhs == retained || m.are_bonded(rhs, retained))   // happens when doing ring CH2's
  {
    array_gsub(to_remove, m.natoms(), flag, -1);
    return 0;
  }

  m.remove_bond_between_atoms(b1->a1(), b1->a2());   // break one end

  to_remove[retained] = 0;

  if (_isotope)
  {
    int c = count_occurrences_of_item_in_array(flag, m.natoms(), to_remove);
//  cerr << "Setting isotope " << c << endl;
    m.set_isotope(retained, c);
  }

  m.add_bond(retained, rhs, SINGLE_BOND);

  return 1;
}

static int
identify_group(Molecule_With_Info_About_Parent & m,
               atom_number_t zatom, 
               int * to_remove,
               int flag,
               resizable_array<const Bond *> & join_points)
{
  to_remove[zatom] = flag;

  int rc = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (flag == to_remove[j])   // already part of group
      continue;

    if (1 == to_remove[j])
      rc += identify_group(m, j, to_remove, flag, join_points);
    else if (to_remove[j] < 0)    // part of failed group - should not happen
      ;
    else
      join_points.add(b);
  }

  return rc;
}

int
Molecular_Abstraction_Compress_Consecutive::process(Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  _molecules_processed ++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  Substructure_Results sresults;

  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

// First identify all the candidate atoms

  int * to_remove = new_int(matoms); std::unique_ptr<int[]> free_to_remove(to_remove);

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    atom_number_t j = e->item(0);

    if (2 != m.ncon(j))
      continue;

    to_remove[j] = 1;
  }

  int next_number_to_assign = 2;

  int initial_matoms = m.natoms();

  for (int i = 0; i < initial_matoms; i++)
  {
    if (1 != to_remove[i])   // not to process, or already processed
      continue;

    resizable_array<const Bond *> join_points;
    if (identify_group(m, i, to_remove, next_number_to_assign, join_points) < 2)
      array_gsub(to_remove, matoms, next_number_to_assign, -1);
    else if (2 == join_points.number_elements())
      _compress_group(m, to_remove, next_number_to_assign, join_points);

    next_number_to_assign++;
  }

  if (count_non_zero_occurrences_in_array(to_remove, matoms) == matoms)
  {
    cerr << "Molecular_Abstraction_Compress_Consecutive::process:cannot remove all atoms '" << m.name() << "', not changing molecule\n";
    return _do_any_writing_needed(empty_molecule, 0, output);
  }

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (to_remove[i] > 1)
      m.remove_atom(i);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}


Set_of_Molecular_Abstractions::~Set_of_Molecular_Abstractions()
{
  if (NULL != _a)
  {
    for (int i = 0; i < _n; i++)
    {
      delete _a[i];
    }

    delete [] _a;
  }

  return;
}


Molecular_Abstraction_Ring_Systems::Molecular_Abstraction_Ring_Systems()
{
  _ele = get_element_from_atomic_number(86);    // Rn

  return;
}

int
Molecular_Abstraction_Ring_Systems::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RINGSYS", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Ring_Systems::build:cannot process '" << token << "'\n";
      return 0;
    }
    else
    {
      cerr << "Molecular_Abstraction_Ring_Systems::buld:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Molecular_Abstraction_Ring_Systems::_replace_ring_system (Molecule_With_Info_About_Parent & m, 
                     int matoms,
                     int flag,
                     const int * ring_systems) const
{
  m.add(_ele);

  int rc = 0;

  matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (flag != ring_systems[i])
      continue;

    rc++;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

//    if (k >= matoms)   // ring atom already added
//      continue;

      if (flag == ring_systems[k])
        continue;

      if (m.are_bonded(k, m.natoms() - 1))
        continue;

      m.add_bond(k, m.natoms() - 1, SINGLE_BOND);
    }
  }

  if (_isotope)
    m.set_isotope(m.natoms() - 1, rc);

  return rc;
}

int
Molecular_Abstraction_Ring_Systems::process(Molecule_With_Info_About_Parent & m,
                                            IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int nr = m.nrings();

  if (0 == nr)
    return _handle_no_match_to_query(m, output);

  int matoms = m.natoms();

  int * ring_systems = new int[matoms]; std::unique_ptr<int[]> free_ring_systems(ring_systems);

  int number_ring_systems = m.label_atoms_by_ring_system_including_spiro_fused(ring_systems);

  for (int i = 1; i <= number_ring_systems; i++)
  {
    _replace_ring_system(m, matoms, i, ring_systems);
  }

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (ring_systems[i])
      m.remove_atom(i);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}
int
Molecular_Abstraction_Remove_Bond::process(Molecule_With_Info_About_Parent & m,
                                           IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  Substructure_Results sresults;

  int nhits = _smarts.substructure_search(&m, sresults);

  if (0 == nhits)
    return _handle_no_match_to_query(m, output);

  int initial_number_fragments = m.number_fragments();

// We need to check the ring membership of the bonds being broken

  Set_of_Atoms b1, b2;

  Set_of_Atoms atoms_to_be_removed;

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    if (e->number_elements() < 2)
    {
      cerr << "Molecular_Abstraction_Remove_Bond::process:not enough matched atoms " << (*e) << endl;
      return 0;
    }

    atom_number_t a1 = e->item(0);
    atom_number_t a2 = e->item(1);

    if (! m.are_bonded(a1, a2))
    {
      cerr << "Molecular_Abstraction_Remove_Bond::process:in " << m.name() << "' atoms " << a1 << " and " << a2 << " not bnded\n";
      continue;
    }

    b1.add(a1);
    b2.add(a2);

    if (_remove_fragment < 0)
      continue;

    if (m.in_same_ring(a1, a2))   // does not make sense to remove fragments
      continue;

    if (0 == _remove_fragment)
      atoms_to_be_removed.add_if_not_already_present(a1);
    else if (1 == _remove_fragment)
      atoms_to_be_removed.add_if_not_already_present(a2);
  }

// Now remove the bonds

  nhits = b1.number_elements();

  for (int i = 0; i < nhits; i++)
  {
    atom_number_t a1 = b1[i];
    atom_number_t a2 = b2[i];

    if (! m.are_bonded(a1, a2))
      continue;

//  cerr << "Removing bond between " << a1 << " and " << a2 << endl;

    m.remove_bond_between_atoms(a1, a2);

    if (_isotope > 0)
    {
      m.set_isotope(a1, _isotope);
      m.set_isotope(a2, _isotope);
    }
  }

//cerr << "After bond removals '" << m.smiles() << "'\n";

  if (initial_number_fragments == m.number_fragments())   // cannot remove fragments, must have broken ring bonds
    ;
  else if (atoms_to_be_removed.number_elements() > 0)
  {
    resizable_array<int> fragments_to_be_removed;

    for (int i = 0; i < atoms_to_be_removed.number_elements(); i++)
    {
      int f = m.fragment_membership(atoms_to_be_removed[i]);

      fragments_to_be_removed.add_if_not_already_present(f);
    }

    m.delete_fragments(fragments_to_be_removed);
  }

  _molecules_changed++;

  return _do_any_writing_needed(m, 1, output);
}


int
Set_of_Molecular_Abstractions::build(const Molecular_Abstraction_Directives_Node & madn)
{
  _n = madn.number_abstractions();

  if (0 == _n)
  {
    cerr << "Set_of_Molecular_Abstractions::build:no abstractions\n";
    return 0;
  }

  _a = new Molecular_Abstraction_Base_Class *[_n];

  const Molecular_Abstraction_Directives_Node * m = &madn;

  for (int i = 0; i < _n; i++)
  {
    assert (NULL != m);

    int t = m->ztype();
//  cerr << "Type t is " << t <<endl;

    if (MAD_TYPE_SCAFFOLD == t)
      _a[i] = new Molecular_Abstraction_Scaffold();
    else if (MAD_TYPE_TRANSLATE == t)
      _a[i] = new Molecular_Abstraction_Transform();
    else if (MAD_TYPE_REMOVE_ATOM == t)
      _a[i] = new Molecular_Abstraction_Remove_Atom();
    else if (MAD_TYPE_RINGS == t)
      _a[i] = new Molecular_Abstraction_Rings();
    else if (MAD_TYPE_BIGRING == t)
      _a[i] = new Molecular_Abstraction_Largest_Ring_System();
    else if (MAD_TYPE_CHANGE_BOND_TYPE == t)
      _a[i] = new Molecular_Abstraction_Change_Bond_Type();
    else if (MAD_TYPE_REPLACE_LINKER == t)
      _a[i] = new Molecular_Abstraction_Replace_Linker();
    else if (MAD_TYPE_ABSTRACT_RING_FORM == t)
      _a[i] = new Molecular_Abstraction_Abstract_Ring_Form();
    else if (MAD_TYPE_REPLACE_LINKER == t)
      _a[i] = new Molecular_Abstraction_Replace_Linker();
    else if (MAD_TYPE_FRAGMENT == t)
      _a[i] = new Molecular_Abstraction_Fragment();
    else if (MAD_TYPE_REMOVE_ATOMS == t)
      _a[i] = new Molecular_Abstraction_Delete_Atoms();
    else if (MAD_TYPE_PLACE_ISOTOPE == t)
      _a[i] = new Molecular_Abstraction_Place_Isotope();
    else if (MAD_TYPE_PLACE_CHARGE == t)
      _a[i] = new Molecular_Abstraction_Place_Charge();
    else if (MAD_TYPE_ALL_ATOMS_TRANSFORM == t)
      _a[i] = new Molecular_Abstraction_All_Transform();
    else if (MAD_TYPE_ALL_BONDS_TRANSFORM == t)
      _a[i] = new Molecular_Abstraction_Change_All_Bonds();
    else if (MAD_TYPE_COMPRESS_CONSECUTIVE == t)
      _a[i] = new Molecular_Abstraction_Compress_Consecutive();
    else if (MAD_TYPE_RINGSYS == t)
      _a[i] = new Molecular_Abstraction_Ring_Systems();
    else if (MAD_TYPE_RMBOND == t)
      _a[i] = new Molecular_Abstraction_Remove_Bond();
    else if (MAD_TYPE_RMRD2 == t)
      _a[i] = new Molecular_Abstraction_Remove_Ring_CH2();
    else if (MAD_TYPE_INVSCAF == t)
      _a[i] = new Molecular_Abstraction_Inverse_Scaffold();
    else if (MAD_TYPE_SSS == t)
      _a[i] = new Molecular_Abstraction_Substructure_Search();
    else if (MAD_TYPE_SPINACH == t)
      _a[i] = new Molecular_Abstraction_Spinach();
    else if (MAD_TYPE_RMSCAFFOLD == t)
      _a[i] = new Molecular_Abstraction_Inverse_Scaffold();
    else
    {
      cerr << "Set_of_Molecular_Abstractions::build:unrecognised form " << t << endl;
      return 0;
    }

    if (! _a[i]->build(*m))
    {
      cerr << "Set_of_Molecular_Abstractions:build:cannot process\n";
      return 0;
    }

    m = m->next();
  }

  return 1;
}

int
Set_of_Molecular_Abstractions::process(Molecule_With_Info_About_Parent & m,
                                       IWString_and_File_Descriptor & output)
{
  if (0 == m.parent_natoms())
    m.compute_parent_natoms();

  for (int i = 0; i < _n; i++)
  {
    _a[i]->process(m, output);

    if (0 == m.natoms())
      break;
  }

  return 1;
}

int
Set_of_Molecular_Abstractions::what_is_being_written(int & writing_smiles,
                                                     int & writing_fingerprints) const
{
  for (int i = 0; i < _n; i++)
  {
    const Molecular_Abstraction_Base_Class * mabi = _a[i];

    if  (! mabi->ok())
      return 0;

    if (mabi->fingerprint_tag().length())
      writing_fingerprints++;

    if (mabi->write_tag().length())
      writing_smiles++;
  }

  if (writing_smiles && writing_fingerprints)
  {
    cerr << "Set_of_Molecular_Abstractions::what_is_being_written:cannot write both smiles and fingerprints\n";
    return 0;
  }

  if (! writing_smiles && ! writing_fingerprints)
  {
    cerr << "Set_of_Molecular_Abstractions::what_is_being_written:no output across " << _n << " transforms\n";
    return 0;
  }

  return 1;
}

int
Set_of_Molecular_Abstractions::report(ostream & os) const
{
  os << "Report on set of " << _n << " abstractions\n";

  for (int i = 0; i < _n; i++)
  {
    _a[i]->report(os);
  }

  return 1;
}

Molecular_Abstraction_Remove_Bond::Molecular_Abstraction_Remove_Bond()
{
  _remove_fragment = -1;

  return;
}

int
Molecular_Abstraction_Remove_Bond::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  const_IWSubstring smarts;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RMBOND", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Remove_Bond::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("REMOVE="))
    {
      token.remove_leading_chars(7);
      if (! token.numeric_value(_remove_fragment) || _remove_fragment < 0 || _remove_fragment > 1)
      {
        cerr << "Molecular_Abstraction_Remove_Bond::build:invalid fragment removal 'rmfrag=" << token << "'\n";
        return 0;
      }
    }
    else if (0 == smarts.length())
      smarts = token;
    else
    {
      cerr << "Molecular_Abstraction_Remove_Bond::buld:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  if (0 == smarts.length())
  {
    cerr << "Molecular_Abstraction_Remove_Bond::no smarts specified\n";
    return 0;
  }
  
  if (! _smarts.create_from_smarts(smarts))
  {
    cerr << "Molecular_Abstraction_Remove_Bond::invalid smarts '" << smarts << "'\n";
    return 0;
  }

  _smarts.set_find_unique_embeddings_only(1);

  return 1;
}

Molecular_Abstraction_Remove_Ring_CH2::Molecular_Abstraction_Remove_Ring_CH2()
{
  _element = get_element_from_atomic_number(6);

  return;
}

Molecular_Abstraction_Inverse_Scaffold::Molecular_Abstraction_Inverse_Scaffold ()
{
  _scaffold_chain_element = NULL;

  return;
}

Molecular_Abstraction_Substructure_Search::Molecular_Abstraction_Substructure_Search()
{
  _must_hit_at_least_one_query = 1;

  return;
}

Molecular_Abstraction_Spinach::Molecular_Abstraction_Spinach ()
{
  _remove_doubly_bonded_atoms_in_spinach = 0;

  _aromatic_element_replacement = NULL;
  _aliphatic_element_replacement = NULL;
  _chain_element_replacement = NULL;

  return;
}

int
Molecular_Abstraction_Remove_Ring_CH2::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  const_IWSubstring zsymbol;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "RMRD2", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Remove_Ring_CH2::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (0 == zsymbol.length())
      zsymbol = token;
    else
    {
      cerr << "Molecular_Abstraction_Remove_Ring_CH2::buld:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  if (0 == zsymbol.length())
    ;
  else if (NULL != (_element = get_element_from_symbol_no_case_conversion(zsymbol)))
    ;
  else
  {
    cerr << "Molecular_Abstraction_Remove_Ring_CH2::unrecognised element '" << zsymbol << "'\n";
    return 0;
  }
  
  return 1;
}

int
Molecular_Abstraction_Inverse_Scaffold::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  const_IWSubstring zsymbol;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "INVSCAF", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Inverse_Scaffold::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if (token.starts_with("SCH="))
    {
      token.remove_leading_chars(4);
      zsymbol = token;
    }
    else
    {
      cerr << "Molecular_Abstraction_Inverse_Scaffold::buld:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  if (0 == zsymbol.length())
    ;
  else if (NULL != (_scaffold_chain_element = get_element_from_symbol_no_case_conversion(zsymbol)))
    ;
  else
  {
    cerr << "Molecular_Abstraction_Inverse_Scaffold::unrecognised element '" << zsymbol << "'\n";
    return 0;
  }

  return 1;
}

int
Molecular_Abstraction_Substructure_Search::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

  const_IWSubstring zsymbol;

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "SSS", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Substructure_Search::build:cannot process '" << token << "'\n";
      return 0;
    }
    else if ("NONM" == token)
    {
      _must_hit_at_least_one_query = 0;
    }
    else if (token.starts_with("SMARTS=") || token.starts_with("SMARTS:"))
    {
      token.remove_leading_chars(7);
      Substructure_Query * q = new Substructure_Query;
      if (! q->create_from_smarts(token))
      {
        cerr << "Molecular_Abstraction_Substructure_Search::build:invalid smarts '" << token << "'\n";
        delete q;
        return 0;
      }

      _queries.add(q);
    }
    else if (process_cmdline_token('q', token, _queries, 0))
      ;
    else
    {
      cerr << "Molecular_Abstraction_Substructure_Search::build:invalid substructure specification '" << token << "'\n";
      return 0;
    }
  }

  if (0 == _queries.number_elements())
  {
    cerr << "Molecular_Abstraction_Substructure_Search::build:no queries\n";
    return 0;
  }

  return _queries.number_elements();
}

/*
  There are many options for processing to spinach form
*/

int
Molecular_Abstraction_Spinach::build (const Molecular_Abstraction_Directives_Node & madn)
{
  const IWString & s = madn.args();

  int i = 0;
  const_IWSubstring token;

//cerr << "Building spinach from '" << s << "'\n";

  while (s.nextword(token, i))
  {
    int fatal;
    if (Molecular_Abstraction_Base_Class::_process(token, "SPINACH", fatal))  // great
      continue;
    else if (fatal)
    {
      cerr << "Molecular_Abstraction_Substructure_Search::build:cannot process '" << token << "'\n";
      return 0;
    }

    if ("RMDBSC" == token)
      _remove_doubly_bonded_atoms_in_spinach = 1;
    else if (token.starts_with("AROM="))
    {
      token.remove_leading_chars(5);
      _aromatic_element_replacement = get_element_from_symbol_no_case_conversion(token);
    }
    else if (token.starts_with("ALIPH="))
    {
      token.remove_leading_chars(6);
      _aliphatic_element_replacement = get_element_from_symbol_no_case_conversion(token);
    }
    else if (token.starts_with("CHAIN="))
    {
      token.remove_leading_chars(6);
      _chain_element_replacement = get_element_from_symbol_no_case_conversion(token);
    }
    else
    {
      cerr << "Molecular_Abstraction_Substructure_Search::build:unrecognised directive '" << token << "' from '" << s << "'\n";
      return 0;
    }
  }

  return 1;
}
 
/*
  We need to make sure the stopping atom is the last item on the list
*/

int
Molecular_Abstraction_Remove_Ring_CH2::_identify_sequence (Molecule_With_Info_About_Parent & m,
                                                atom_number_t avoid,
                                                atom_number_t zatom,
                                                Set_of_Atoms & s) const
{
  if (s.contains(zatom))    // isolated ring
    return 0;

  s.add(zatom);

  if (1 != m.nrings(zatom))
    return 0;

  const Atom * a = m.atomi(zatom);

  if (_element != a->element())
    return 1;

  if (2 != a->ncon())
    return 1;

  if (avoid == a->other(zatom, 0))
    return 1 + _identify_sequence(m, zatom, a->other(zatom, 1), s);
  else
    return 1 + _identify_sequence(m, zatom, a->other(zatom, 0), s);
}

int
Molecular_Abstraction_Remove_Ring_CH2::_identify_ch2_to_remove (Molecule_With_Info_About_Parent & m,
                        atom_number_t zatom,
                        Set_of_Atoms & lhs,
                        Set_of_Atoms & rhs) const
{
  const Atom * a = m.atomi(zatom);

  int rc = _identify_sequence(m, zatom, a->other(zatom, 0), lhs);
  rc += _identify_sequence(m, zatom, a->other(zatom, 1), rhs);

  return rc;
}

int
Molecular_Abstraction_Remove_Ring_CH2::process (Molecule_With_Info_About_Parent & m,
                                IWString_and_File_Descriptor & output)
{
  if (0 == m.nrings())
    return _handle_no_match_to_query(m, output);

  int rc = 0;

  int matoms = m.natoms();

  int * to_remove = new int[matoms]; std::unique_ptr<int[]> free_to_remove(to_remove);

//#define DEBUG_RMRD2
#ifdef DEBUG_RMRD2
  cerr << "Molecule has " << matoms << " atoms, " << m.smiles() << endl;
#endif

  for (int i = 0; i < matoms; i++)
  {
    if (1 != m.nrings(i))
      continue;

    const Atom * a = m.atomi(i);

    if (_element != a->element())
      continue;

    if (2 != a->ncon())
      continue;

#ifdef DEBUG_RMRD2
    cerr << "Looking for removals either size of " << i << " " << m.smarts_equivalent_for_atom(i) << endl;
#endif

    Set_of_Atoms lhs, rhs;

    _identify_ch2_to_remove(m, i, lhs, rhs);

#ifdef DEBUG_RMRD2
    cerr << "Molecule with " << m.natoms() << " from " << i << " removing " << lhs << " and " << rhs << endl;
#endif

    set_vector(to_remove, matoms, 0);

//  If we have an isolated ring, we may need to shrink it
  
    if (rhs.number_elements() == lhs.number_elements() && i == lhs.last_item() && i == rhs.last_item())    // isolated ring, done
    {
      if (3 == rhs.number_elements())
        continue;

#ifdef DEBUG_RMRD2
      cerr << "Special processing for isolated ring " << lhs << ' ' << rhs << endl;
#endif

      lhs.pop();
      rhs.pop();
      atom_number_t a1 = lhs.pop();
      atom_number_t a2 = rhs.pop();
      lhs.remove_item(0);
      lhs.set_vector(to_remove, 1);
      m.add_bond(a1, a2, SINGLE_BOND);
    }
    else
    {
      atom_number_t a1 = lhs.pop();
      atom_number_t a2 = rhs.pop();
      if (a1 != a2)
        ;
      else if (lhs.number_elements() > rhs.number_elements())
        a2 = lhs.pop();
      else if (lhs.number_elements() < rhs.number_elements())
        a2 = rhs.pop();
      else if (0 == lhs.number_elements() && m.are_bonded(a1, a2))   // already 3 membered ring
        continue;
      else
        continue;

      lhs.set_vector(to_remove, 1);
      rhs.set_vector(to_remove, 1);
      if (m.are_bonded(a1, a2))
      {
        if (! m.are_bonded(i, a1))
          m.add_bond(i, a1, SINGLE_BOND);
        if (! m.are_bonded(i, a2))
          m.add_bond(i, a2, SINGLE_BOND);
      }
      else
      {
        m.add_bond(a1, a2, SINGLE_BOND);
        to_remove[i] = 1;
      }

#ifdef DEBUG_RMRD2
      cerr << "Determined " << a1 << " and " << a2 << endl;
#endif
    }

#ifdef DEBUG_RMRD2
    for (int j = 0; j < matoms; j++)
    {
      if (to_remove[j])
        cerr << "removing " << j << " " << m.smarts_equivalent_for_atom(j) << endl;
    }
#endif

    m.remove_atoms(to_remove);

//  if we removed atoms less than I, we need to adjust it down

    i = i - count_non_zero_occurrences_in_array(to_remove, i + 1);

#ifdef DEBUG_RMRD2
    cerr << "Molecule now '" << m.smiles() << "'\n";
#endif

    matoms = m.natoms();

    rc++;
  }

  if (0 == rc)
    return _handle_no_match_to_query(m, output);

  _molecules_changed++;

  return _do_any_writing_needed(m, rc, output);
}

static void
all_bonds_single (Molecule_With_Info_About_Parent & m,
                  atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    m.set_bond_type_between_atoms(zatom, j, SINGLE_BOND);
  }

  return;
}

#define SCAFFOLD_CHAIN 1
#define SCAFFOLD_AROMATIC 2
#define SCAFFOLD_ALIPHATIC 3

int
Molecular_Abstraction_Inverse_Scaffold::process (Molecule_With_Info_About_Parent & m,
                                IWString_and_File_Descriptor & output)
{
  _molecules_processed++;

  int matoms = m.natoms();

  if (0 == matoms)
    return _do_any_writing_needed(m, 0, output);

  int nr = m.nrings();

  if (0 == nr)
    return _handle_no_match_to_query(empty_molecule, output);

  _molecules_changed++;

  m.remove_all_chiral_centres();

  int * scaffold = new_int(matoms); std::unique_ptr<int[]> free_scaffold(scaffold);

  _identify_scaffold (m, scaffold, 0);    // 0 means exclude attachment point

  for (int i = 0; i < matoms; i++)
  {
    if (0 == scaffold[i])
      continue;

    if (0 == m.nrings(i))
      scaffold[i] = SCAFFOLD_CHAIN;
    else if (m.is_aromatic(i))
      scaffold[i] = SCAFFOLD_AROMATIC;
    else
      scaffold[i] = SCAFFOLD_ALIPHATIC;
  }

  for (int i = matoms - 1; i >= 0; i--)
  {
    if (0 == scaffold[i])          // no change to these
      continue;

    if (SCAFFOLD_AROMATIC == scaffold[i])
    {
      m.set_atomic_number(i, 18);      // Ar
      all_bonds_single(m, i);
    }
    else if (SCAFFOLD_ALIPHATIC == scaffold[i])
    {
      m.set_atomic_number(i, 13);     // Al
      all_bonds_single(m, i);
    }
    else if (SCAFFOLD_CHAIN == scaffold[i])
    {
      if (NULL != _scaffold_chain_element)
        m.set_element(i, _scaffold_chain_element);
    }
    else
      m.remove_atom(i);
  }

  if (NULL != _scaffold_chain_element)    // remove doubly bonded connections to scaffold
  {
    for (int i = m.natoms() - 1; i >= 0; i--)
    {
      if (0 == scaffold[i])    // we want to delete some scaffold atoms
        continue;

      const Atom * a = m.atomi(i);

      if (1 == a->ncon())
        m.remove_atom(i);
    }
  }

  return _do_any_writing_needed(m, 1, output);
}

int
Molecular_Abstraction_Substructure_Search::process (Molecule_With_Info_About_Parent & m,
                                IWString_and_File_Descriptor & output)
{
  Molecule_to_Match target(&m);

  for (int i = 0; i < _queries.number_elements(); i++)
  {
    if (_queries[i]->substructure_search(target))
    {
      if (_must_hit_at_least_one_query)
        return _do_any_writing_needed(m, 1, output);

//    matches are rejections, zorch the molecule

      m.resize(0);

      return 0;
    }
  }

// Did not hit any of the queries, zorch the molecule

  if (_must_hit_at_least_one_query)
  {
    m.resize(0);
    return 0;
  }

// did not hit any queries and they were rejections

  return _do_any_writing_needed(m, 1, output);
}

int
Molecular_Abstraction_Spinach::process (Molecule_With_Info_About_Parent & m,
                                        IWString_and_File_Descriptor & output)
{
  if (0 == m.nrings())
    return _handle_no_match_to_query(empty_molecule, output);

  const auto matoms = m.natoms();

  int * spch = new int[matoms + matoms + matoms]; std::unique_ptr<int[]> free_spch(spch);
  int * arom = spch + matoms;
  int * nrings = spch + matoms + matoms;

  m.ring_membership(nrings);
  m.aromaticity(arom);

  m.identify_spinach(spch);

  const auto ne = m.nedges();

  for (auto i = 0; i < ne; ++i)
  {
    const Bond * b = m.bondi(i);

    if (b->is_single_bond())
      continue;

    const auto a1 = b->a1();

    if (spch[a1])
      continue;

    const auto a2 = b->a2();

    if (spch[a2])
      continue;

    m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
  }

  if (NULL != _aromatic_element_replacement || NULL != _aliphatic_element_replacement || NULL != _chain_element_replacement || _remove_doubly_bonded_atoms_in_spinach)
  {
    for (auto i = matoms - 1; i >= 0; --i)
    {
      if (spch[i])    // part of spinach, we are processing scaffold atoms here
        continue;

      const auto a = m.atomi(i);

      if (_remove_doubly_bonded_atoms_in_spinach && 1 == a->ncon())
        m.remove_atom(i);
      else if (NULL != _aromatic_element_replacement && arom[i])
      {
        m.set_element(i, _aromatic_element_replacement);
        m.set_formal_charge_if_different(i, 0);
      }
      else if (NULL != _chain_element_replacement && 0 == nrings[i])
      {
        m.set_element(i, _chain_element_replacement);
        m.set_formal_charge_if_different(i, 0);
      }
      else if (NULL != _aliphatic_element_replacement && nrings[i] && 0 == arom[i])
      {
        m.set_element(i, _aliphatic_element_replacement);
        m.set_formal_charge_if_different(i, 0);
      }
    }
  }

  return _do_any_writing_needed(m, 1, output);
}
