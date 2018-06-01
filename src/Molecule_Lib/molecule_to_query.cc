#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "misc.h"

#include "molecule_to_query.h"

static int File_Scope_substituents_only_at_isotopic_atoms = 0;

void
set_substituents_only_at_isotopic_atoms (int s)
{
  File_Scope_substituents_only_at_isotopic_atoms = s;
}

int
substituents_only_at_isotopic_atoms()
{
  return File_Scope_substituents_only_at_isotopic_atoms;
}

static int File_Scope_must_have_substituent_at_every_isotopic_atom = 1;

void
set_must_have_substituent_at_every_isotopic_atom (int s)
{
  File_Scope_must_have_substituent_at_every_isotopic_atom = s;
}

int 
must_have_substituent_at_every_isotopic_atom()
{
  return File_Scope_must_have_substituent_at_every_isotopic_atom;
}

static int File_Scope_substitutions_only_at_non_isotopic_atoms = 0;

void
set_substitutions_only_at_non_isotopic_atoms (int s)
{
  File_Scope_substitutions_only_at_non_isotopic_atoms = s;
}

int
substitutions_only_at_non_isotopic_atoms ()
{
  return File_Scope_substitutions_only_at_non_isotopic_atoms;
}

/*
  2005. We are getting lots of scaffold offerings that have isotopic
  labels. We get an ortho substituted aniline. But we also get things
  that is hard to know what this means. Does it mean that there must
  be two connections at the Nitrogen, or does it mean at least
  one connection?
*/

static int File_Scope_isotope_count_means_extra_connections = 0;

void
set_isotope_count_means_extra_connections (int s)
{
  File_Scope_isotope_count_means_extra_connections = s;
}

int
isotope_count_means_extra_connections()
{
  return File_Scope_isotope_count_means_extra_connections;
}

static int respect_ring_membership = 0;

void
set_respect_ring_membership (int s)
{
  respect_ring_membership = s;
}

static int molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms = 0;

void
set_molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms (int s)
{
  molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms = s;
}

static int File_Scope_only_include_isotopically_labeled_atoms = 0;

int
only_include_isotopically_labeled_atoms()
{
  return File_Scope_only_include_isotopically_labeled_atoms;
}

void
set_only_include_isotopically_labeled_atoms(int s)
{
  File_Scope_only_include_isotopically_labeled_atoms = s;
}

static int File_Scope_only_aromatic_atoms_match_aromatic_atoms = 0;

void
set_only_aromatic_atoms_match_aromatic_atoms(int s)
{
  File_Scope_only_aromatic_atoms_match_aromatic_atoms = s;
}


Molecule_to_Query_Specifications::Molecule_to_Query_Specifications ()
{
  _make_embedding = 1;
  _built_from_isis_reaction_file = 0;
  _all_ring_bonds_become_undefined = 0;
  _non_ring_atoms_become_nrings_0 = respect_ring_membership;
  _atoms_conserve_ring_membership = 0;
  _ring_atoms_conserve_ring_membership = 0;
  _copy_bond_attributes = 0;
  _only_aromatic_atoms_match_aromatic_atoms = File_Scope_only_aromatic_atoms_match_aromatic_atoms;

  _atoms_in_molecule = 0;

  _ncon = -1;
  _min_ncon = -1;
  _max_ncon = -1;

  _condense_explicit_hydrogens_to_anchor_atoms = molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms;

  _just_atomic_number_and_connectivity = 0;

  _discern_hcount = 1;

  _nrings_is_ring_bond_count = 0;

  _min_extra_atoms_in_target = -1;
  _max_extra_atoms_in_target = -1;

  _min_fraction_atoms_matched = 0.0f;
  _max_fraction_atoms_matched = 1.0f;

  _use_preference_values_to_distinguish_symmetry = 0;

  _convert_explicit_hydrogens_to_match_any_atom = 0;
    
  _convert_all_aromatic_atoms_to_generic_aromatic = 0;

  _query_must_match_both_explicit_and_implicit_hydrogens = 0;

  _preserve_saturation = 0;

  _ignore_molecular_hydrogen_information = 0;

  _interpret_atom_alias_as_smarts = 1;

  _convert_explicit_hydrogens_to_match_any_atom_including_hydrogen = 0;

  _substituents_only_at_isotopic_atoms = 0;

  _isotopic_label_means = 0;

  _set_element_hits_needed_during_molecule_to_query = 1;

  _aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic = 0;

  _preserve_smallest_ring_size = 0;

  _bonds_preserve_ring_membership = 0;

  return;
}

Molecule_to_Query_Specifications::~Molecule_to_Query_Specifications ()
{
  assert (-5 != _atoms_in_molecule);

  _atoms_in_molecule = -5;

  return;
}


int
Molecule_to_Query_Specifications::smarts_for_element (const Element * e,
                                                      IWString & s) const
{
  for (int i = 0; i < _element_to_smarts.number_elements (); i++)
  {
    const Element_to_Smarts * ei = _element_to_smarts[i];

    if (e == ei->element ())
    {
      s = ei->smarts ();
      return 1;
    }
  }

  return 0;
}


int
Molecule_to_Query_Specifications::parse_directives (const const_IWSubstring & directives)
{
  int i = 0;
  const_IWSubstring token;

  while (directives.nextword (token, i, DIRECTIVE_SEPARATOR_TOKEN))
  {
    if (! _parse_directive (token))
    {
      cerr << "Molecule_to_Query_Specifications::parse_directives:invalid directive '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Molecule_to_Query_Specifications::_parse_directive (const_IWSubstring dir)    // passed by value, our own copy
{
  if (dir.starts_with ("onlysub="))
  {
    dir.remove_up_to_first ('=');
    if (! _parse_onlysub_directive (dir))
    {
      cerr << "Molecule_to_Query_Specifications::_parse_onlysub_directive:invalid onlysub directive\n";
      return 0;
    }
  }
  else if ("onlysubiso" == dir)
  {
    _substituents_only_at_isotopic_atoms = 1;

//  File_Scope_substituents_only_at_isotopic_atoms = 1;   // need to get this into the object

    _substitutions_only_at.create_from_smarts ("[!0*]");    // this may not be necessary now
  }
  else if ("oama" == dir)
  {
    _only_aromatic_atoms_match_aromatic_atoms = 1;
  }
  else
  {
    cerr << "Molecule_to_Query_Specifications::_parse_directive:unrecognised directive '" << dir << "'\n";
    return 0;
  }

  return 1;
}

int
Molecule_to_Query_Specifications::_parse_onlysub_directive (const const_IWSubstring smarts)
{
  if (_substitutions_only_at.active ())
  {
    cerr << "Molecule_to_Query_Specifications::_parse_onlysub_directive:query already active\n";
    return 0;
  }

  if (! _substitutions_only_at.create_from_smarts (smarts))
  {
    cerr << "Molecule_to_Query_Specifications::_parse_onlysub_directive:invalid smarts '" << smarts << "'\n";
    return 0;
  }
  
  return 1;
}






Element_to_Smarts::Element_to_Smarts ()
{
  _e = NULL;

  return;
}

int
Element_to_Smarts::build (const const_IWSubstring & buffer)
{
  if (buffer.contains (' '))
  {
    cerr << "Element_to_Smarts::build:whitespace not allowed\n";
    return 0;
  }

  const_IWSubstring e, s;

  if (! buffer.split (e, '=', s) || 0 == e.length () || 0 == s.length ())
  {
    cerr << "Element_to_Smarts::build:invalid specification '" << buffer << "'\n";
    return 0;
  }

  _e = get_element_from_symbol_no_case_conversion(e);
  if (NULL == _e)
  {
    set_auto_create_new_elements (1);
    set_atomic_symbols_can_have_arbitrary_length (1);

    _e = create_element_with_symbol (e);
  }

  _smarts = s;

  return 1;
}

template class resizable_array_p<Element_to_Smarts>;
template class resizable_array_base<Element_to_Smarts *>;

int
Molecule_to_Query_Specifications::set_smarts_for_atom (const const_IWSubstring & f)
{
  Element_to_Smarts * e = new Element_to_Smarts;

  if (! e->build (f))
  {
    cerr << "Molecule_to_Query_Specifications::set_smarts_for_atom:invalid element to smarts '" << f << "'\n";
    delete e;
    return 0;
  }

  _element_to_smarts.add (e);

  return 1;
}
