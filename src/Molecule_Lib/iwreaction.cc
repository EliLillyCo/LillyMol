#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>

using std::cerr;

#include "Foundational/iwmisc/misc.h"

#define ISTREAM_AND_TYPE_IMPLEMENTATION

#include "chiral_centre.h"
#include "istream_and_type.h"
#include "iwreaction.h"
#include "misc2.h"
#include "molecule.h"
#include "molecule_to_query.h"
#include "output.h"
#include "path.h"
#include "target.h"

template class data_source_and_type<Molecule_and_Embedding>;

static int strip_reagents_to_largest_fragment = 0;

static int transfer_text_info = 0;

void
set_reaction_transfer_text_info_to_products (int s)
{
  transfer_text_info = s;
}

static IWString component_separator = " + ";

void
set_component_separator(const const_IWSubstring & s)
{
  if ("NONE" == s)
    component_separator.resize(0);
  else
    component_separator = s;

  return;
}

void
set_strip_reagents_to_largest_fragment (int s)
{
  strip_reagents_to_largest_fragment = s;
}

static Molecule_Output_Object stream_for_sidechains_not_matching_query;

int
set_stream_for_sidechains_not_matching_query (const const_IWSubstring & fname)
{
  stream_for_sidechains_not_matching_query.add_output_type(FILE_TYPE_SMI);

  if (! stream_for_sidechains_not_matching_query.new_stem(fname))
  {
    cerr << "Cannot open stream for sidechains not reacting '" << fname << "'\n";
    return 0;
  }

  return 1;
}

static int iwreaction_display_no_atoms_in_query_message = 1;

void set_iwreaction_display_no_atoms_in_query_message(int s)
{
  iwreaction_display_no_atoms_in_query_message = s;
}

static int iwreaction_display_take_first_reagent_from_each_sidechain_warning_message = 1;

//#ifdef NOT_SURE_THIS_IS_A_GOOD_IDEA_DISALLOW_FOR_NOW  -
// in molecular_transformations, which uses this, we are now allowing the -Z option to ignore reactions that have no reactants,
// so we need to allow NOT showing this message  - Tad Hurst 20180209

void
set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(const int s)
{
  iwreaction_display_take_first_reagent_from_each_sidechain_warning_message = s;
}
//#endif

static int warn_duplicate_atom_deletion = 1;
void
set_warn_duplicate_atom_deletion(int s) {
  warn_duplicate_atom_deletion = s;
}

static int warn_non_bonded_breaks = 1;
void
set_warn_non_bonded_breaks(int s) {
  warn_non_bonded_breaks = s;
}

#define WRITE_NUMBERED_SMILES_ON_ERROR
#ifdef WRITE_NUMBERED_SMILES_ON_ERROR

static int
write_numbered_smiles (Molecule & m, 
                       std::ostream & os)
{
  int matoms = m.natoms();
  for (int i = 1; i < matoms; i++)
  {
    m.set_isotope(i, i);
  }

  os << m.smiles() << ' ' << m.name() << '\n';

  m.transform_to_non_isotopic_form();

  return os.good();
}

#endif

Enumeration_Temporaries::Enumeration_Temporaries (int s)
{
  assert (s >= 0);

  if (0 == s)
  {
    _atoms_in_growing_molecule = nullptr;
    _reagent = nullptr;

    return;
  }

  _atoms_in_growing_molecule = new int[s];

  _reagent = new const Molecule_and_Embedding *[s];

  for (int i = 0; i < s; i++)    // remove this loop once things are debugged
  {
    _reagent[i] = nullptr;
  }

  return;
}

Enumeration_Temporaries::~Enumeration_Temporaries()
{
  if (nullptr != _atoms_in_growing_molecule)
  {
    delete [] _atoms_in_growing_molecule;
    delete [] _reagent;
  }

  return;
}

IWReaction::IWReaction()
{
  _append_names = 0;

  _query_files_in_current_directory = 1;

  _find_kekule_forms_for_bad_valence = 0;

  _make_implicit_hydrogens_explicit = 0;

  _user_specified_void_ptr = nullptr;

  _match_via_atom_map = 0;

  _has_sidechain_isotope_requirement = 0;

  return;
}

IWReaction::~IWReaction()
{
}

int
IWReaction::number_reagents() const
{
  int rc = 0;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (s->single_reagent_only())    // HUH, why is this here?
      continue;

    rc += s->number_reagents();
  }

  return rc;
}

int
IWReaction::number_products_per_scaffold_embedding() const
{
  int rc = 1;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (0 == s->number_reagents())
      cerr << "IWReaction::number_products_per_scaffold_embedding: sidechain " << i << " with no reagents!\n";
    else
      rc = rc * s->number_reagents();
  }

  return rc;
}

/*
  Jun 2017. Not sure if this is a good idea or not. We may need to retain
  the ability to set the scaffold separately...
*/

void
IWReaction::set_one_embedding_per_start_atom(const int s)
{
//cerr << "IWReaction::set_one_embedding_per_start_atom:setting " << s << '\n';

  _match_conditions.set_one_embedding_per_start_atom(s);
  this->set_find_one_embedding_per_atom(s);

  for (int i = 0; i < _sidechains.number_elements(); ++i)
  {
    _sidechains[i]->set_find_one_embedding_per_atom(s);
  }

  return;
}

Molecule_and_Embedding *
IWReaction::reagent (const Reaction_Iterator & iterator) const
{
  assert (iterator.active());

  int ns = _sidechains.number_elements();
  if (0 == ns)
    return nullptr;

  if (ns > 1)      // ambiguous as to what is the "current" sidechain. Should be more careful and see if there is just one sidechain with varying reagents...
    return nullptr;

  int s = 0;

  const Sidechain_Reaction_Site * sc = _sidechains[s];

  int r = iterator.reagent(s);

  return sc->reagent(r);
}

/*
  Sometimes a reaction object already has all the reagents needed.
*/

int
IWReaction::all_sidechains_have_reagents() const
{
  assert (ok());

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    if (! _sidechains[i]->single_reagent_only())
      return 0;
  }

  return 1;    // yep, all the sidechains have their own reagent
}

int
IWReaction::number_sidechains_with_reagents() const
{
  int rc = 0;

  for (int i = 0; i < _sidechains.number_elements(); i++)
  {
    if (_sidechains[i]->number_reagents() > 0)
      rc++;
  }

  return rc;
}

Reaction_Change_Element::Reaction_Change_Element()
{
  _atom = -1;
  _element = nullptr;
}

/*
  An element transformation will look like

    (A C change_element "3 S")

  which means change atom 3 to a sulphur
*/

int
Reaction_Change_Element::construct_from_msi_attribute(const msi_attribute * att)
{
  const IWString & rce = att->stringval();

  if (2 != rce.nwords())
  {
    cerr << "Reaction_Change_Element::construct_from_msi_attribute: must have two string attributes\n";
    cerr << (*att) << '\n';
    cerr << rce << '\n';
    return 0;
  }

  const_IWSubstring token;
  rce.word(0, token);

  if (! token.numeric_value(_atom) || _atom < 0)
  {
    cerr << "Reaction_Change_Element::construct_from_msi_attribute: atom specifiers must be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  (void) rce.word(1, token);    // get the new element type

  _element = get_element_from_symbol_no_case_conversion(token);
  if (nullptr == _element)
  {
    cerr << "Reaction_Change_Element::construct_from_msi_attribute: unrecognised element '" << token << "'\n";
    cerr << (*att) << '\n';
    return 0;
  }

  return 1;
}

int
Reaction_Change_Element::write_msi(std::ostream & os, const const_IWSubstring & ind,
                                    const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"" << _atom << ' ' << _element->symbol() << "\")\n";

  return os.good();
}

int
Reaction_Change_Element::process (Molecule & result,
                                  const Set_of_Atoms & embedding,
                                  int offset)
{
  assert (embedding.ok_index(_atom));

  atom_number_t a = embedding[_atom];

  return result.set_element(a + offset, _element);   // will crash if (a + offset) illegal
}

Reaction_Formal_Charge::Reaction_Formal_Charge()
{
  _atom = -1;
  _fc = 0;
}

Reaction_Change_Formal_Charge::Reaction_Change_Formal_Charge()
{
  _atom = -1;
  _delta = 0;
}

int
Reaction_Formal_Charge::construct_from_msi_attribute (const msi_attribute * att)
{
  if (2 != att->number_int_values())
  {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: must have two int attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom = att->int_multi_value(0);

  if (_atom < 0)
  {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: atom specifiers must be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _fc = att->int_multi_value(1);

  return 1;
}

int
Reaction_Change_Formal_Charge::construct_from_msi_attribute (const msi_attribute * att)
{
  if (2 != att->number_int_values())
  {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: must have two int attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom = att->int_multi_value(0);

  if (_atom < 0)
  {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: atom specifiers must be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _delta = att->int_multi_value(1);

  return 1;
}

int
Reaction_Formal_Charge::write_msi(std::ostream & os, const const_IWSubstring & ind,
                           const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A I " << attribute_name << " \"" << _atom << ' ' << _fc << "\")\n";

  return os.good();
}

int
Reaction_Change_Formal_Charge::write_msi(std::ostream & os, const const_IWSubstring & ind,
                               const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A I " << attribute_name << " \"" << _atom << ' ' << _delta << "\")\n";

  return os.good();
}

int
Reaction_Formal_Charge::process (Molecule & result,
                                 const Set_of_Atoms & embedding,
                                 int offset)
{
  assert (embedding.ok_index(_atom));

  atom_number_t a = embedding[_atom];

  result.set_formal_charge(a + offset, _fc);   // will crash if (a + offset) illegal

  return 1;
}

int
Reaction_Change_Formal_Charge::process (Molecule & result,
                                        const Set_of_Atoms & embedding,
                                        int offset)
{
  assert (embedding.ok_index(_atom));

  atom_number_t a = embedding[_atom];

  formal_charge_t fc = result.formal_charge(a + offset) + _delta;

//cerr << "Initially formal charge is " << result.formal_charge (a + offset) << " will change to " << fc << '\n';

  result.set_formal_charge(a + offset, fc);   // will crash if (a + offset) illegal

  return 1;
}


Reaction_Place_Isotope::Reaction_Place_Isotope()
{
  _isotope = 0;
}

int
Reaction_Place_Isotope::construct_from_msi_attribute (const msi_attribute * att)
{
  if (2 != att->number_int_values())
  {
    cerr << "Reaction_Place_Isotope::construct_from_msi_attribute: must have two int attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom << att->int_multi_value(0);

  if (_atom[0] < 0)
  {
    cerr << "Reaction_Place_Isotope::construct_from_msi_attribute: atom specifiers must be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _isotope = att->int_multi_value(1);

  return 1;
}

int
Reaction_Place_Isotope::write_msi(std::ostream & os, const const_IWSubstring & ind,
                           const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A I " << attribute_name << " (" << _atom[0] << ' ' << _isotope << "))\n";

  return os.good();
}

int
Reaction_Place_Isotope::process (Molecule & result,
                                 const Set_of_Atoms & embedding,
                                 int offset)
{
  for (int matched_atom : _atom) {
    const atom_number_t a = embedding[matched_atom];
    assert (result.ok_atom_number(a));

    result.set_isotope(a + offset, _isotope);
  }

  return 1;
}

int
Reaction_Increment_Isotope::process (Molecule & result,
                                 const Set_of_Atoms & embedding,
                                 int offset)
{
  for (const int matched_atom : _atom) {
    const atom_number_t a = embedding[matched_atom];

    // cerr << "matched atom " << matched_atom << " current isotope is " << result.isotope(a+offset) << " delta " << _isotope << '\n';

    const isotope_t iso = result.isotope(a + offset) + _isotope;

    result.set_isotope(a + offset, iso);
  }

  return 1;
}

int
Reaction_Invert_Isotope::process(Molecule & result,
                                 const Set_of_Atoms & embedding,
                                 int offset)
{

  for (const int matched_atom : _atom) {
    assert (embedding.ok_index(matched_atom));

    const atom_number_t a = embedding[matched_atom];

    const isotope_t iso = result.isotope(a + offset);    // current value

  // cerr << "Inverting isotope to " << _isotope << '\n';

    if (iso == 0) {
      result.set_isotope(a + offset, _isotope);   // will crash if (a + offset) illegal
    } else {
      result.set_isotope(a + offset, 0);
    }
  }

  return 1;
}

/*
  We frequently need to see if an atom is mentioned in a set of bonds
*/

static int
atom_in_list_of_bonds (const resizable_array_p<Bond> & bonds,
                       int zatom)
{
  int nb = bonds.number_elements();
  for (int i = 0; i < nb; i++)
  {
    if (bonds[i]->involves(zatom))
      return 1;
  }

  return 0;
}

static int
check_atoms (const resizable_array<int> & atoms,
             int atoms_in_scaffold_query)
{
  int n = atoms.number_elements();
  for (int i = 0; i < n; i++)
  {
    if (atoms[i] >= atoms_in_scaffold_query)
    {
      cerr << "IWReaction::check_internal_consistency: invalid atom specifier " << atoms[i] << '\n';
      cerr << "Query has only " << atoms_in_scaffold_query << " atoms\n";
      return 0;
    }
  }

  return 1;
}

int
IWReaction::check_internal_consistency()
{
  if (! Reaction_Site::check_internal_consistency())
    return 0;

  int atoms_in_scaffold_query = max_atoms_in_query();

  int ns = _sidechains.number_elements();
  for (int i = 0; i < ns; i++)
  {
    if (! _sidechains[i]->check_internal_consistency(ns, atoms_in_scaffold_query, _atoms_to_be_removed, _fragments_to_be_removed))
    {
      cerr << "IWReaction::check_internal_consistency: sidechain " << i << " is invalid\n";
      return 0;
    }
  }

  return 1;
}

static int
check_bonds (const resizable_array_p<Bond> & bond,
             int atoms_in_query)
{
  int nb = bond.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = bond[i];

    if (b->a1() >= atoms_in_query || b->a2() >= atoms_in_query)
    {
      cerr << "check_bonds:invalid bond specifier\n";
      cerr << "query has " << atoms_in_query << " query atoms max\n";
      cerr << "Request for " << b->a1() << " or " << b->a2() << " is invalid\n";
      return 0;
    }
  }

  return 1;
}

void
Reaction_Site::set_find_unique_embeddings_only(int s) {
  if (_queries.number_elements() > 0) {
    for (Substructure_Query * q : _queries) {
      q->set_find_unique_embeddings_only(s);
    }
  } else {
    Substructure_Query::set_find_unique_embeddings_only(s);
  }
}

void
Reaction_Site::set_find_one_embedding_per_atom(int s) {
  if (_queries.number_elements() > 0) {
    for (Substructure_Query * q : _queries) {
      q->set_find_one_embedding_per_atom(s);
    }
  } else {
    Substructure_Query::set_find_one_embedding_per_atom(s);
  }
}

void
Reaction_Site::set_do_not_perceive_symmetry_equivalent_matches(int s) {
  if (_queries.number_elements() > 0) {
    for (Substructure_Query * q : _queries) {
      q->set_do_not_perceive_symmetry_equivalent_matches(s);
    }
  } else {
    Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches(s);
  }
}

void
Reaction_Site::set_max_matches_to_find(int s) {
  if (_queries.number_elements() > 0) {
    for (Substructure_Query * q : _queries) {
      q->set_max_matches_to_find(s);
    }
  } else {
    Substructure_Query::set_max_matches_to_find(s);
  }
}

void
Reaction_Site::set_embeddings_do_not_overlap(int s) {
  if (_queries.number_elements() > 0) {
    for (Substructure_Query * q : _queries) {
      q->set_embeddings_do_not_overlap(s);
    }
  } else {
    Substructure_Query::set_embeddings_do_not_overlap(s);
  }
}

int
Reaction_Site::check_internal_consistency()
{
  int atoms_in_query = 0;
  if (_queries.empty()) {
    atoms_in_query = Substructure_Query::max_atoms_in_query();
  } else {
    atoms_in_query = _queries[0]->max_atoms_in_query();
  }

  if (atoms_in_query > 0)
    ;
  else if (_noop_reaction)
    return 1;
  else
  {
    if (iwreaction_display_no_atoms_in_query_message)
      cerr << "Reaction_Site::check_internal_consistency:no atoms in query. Probably fixed reagent in ISIS file\n";
    return 1;
  }

//cerr << "Query has " << atoms_in_query << " atoms in the query\n";

  if (! check_bonds(_bonds_to_be_broken, atoms_in_query))
    return 0;

  if (! check_bonds(_bonds_to_be_made, atoms_in_query))
    return 0;

  if (! check_atoms(_atoms_to_be_removed, atoms_in_query))
    return 0;

// Are any of the atoms or fragments being removed involved in bonds?

  int n = _atoms_to_be_removed.number_elements();
  for (int i = 0; i < n; i++)
  {
    if (atom_in_list_of_bonds(_bonds_to_be_made, _atoms_to_be_removed[i]))
    {
      cerr << "Reaction_Site::check_internal_consistency: atom to be removed " << _atoms_to_be_removed[i] <<
              " is involved in a bond to be made\n";
      return 0;
    }
  }
  
  n = _fragments_to_be_removed.number_elements();

  if (n && _fragments_to_be_kept.number_elements())
  {
    cerr << "Reaction_Site::check_internal_consistency:cannot specify both fragments to remove and fragments to keep\n";
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    if (atom_in_list_of_bonds(_bonds_to_be_made, _fragments_to_be_removed[i]))
    {
      cerr << "Reaction_Site::check_internal_consistency: fragment to be removed " << _fragments_to_be_removed[i] <<
              " is involved in a bond to be made\n";
      return 0;
    }
  }

  return 1;
}

int
Sidechain_Reaction_Site::_check_inter_particle_atoms(const int atoms_in_scaffold_query,
                            const int scaffold_atom,
                            const resizable_array<int> & scaffold_atoms_to_be_removed,
                            const resizable_array<int> & scaffold_fragments_to_be_removed,
                            const int sidechain_atom) const
{
  if (scaffold_atom >= atoms_in_scaffold_query)
  {
    cerr << "Sidechain_Reaction_Site::check_internal_consistency: invalid inter particle bond\n";
    cerr << "Scaffold query contains " << atoms_in_scaffold_query << " atoms, so " << scaffold_atom << " is invalid\n";
    return 0;
  }

// Does the scaffold delete its end of an inter-particle bond

  if (scaffold_atoms_to_be_removed.contains(scaffold_atom) || scaffold_fragments_to_be_removed.contains(scaffold_atom))
  {
    cerr << "Sidechain_Reaction_Site::check_internal_consistency: invalid inter particle bond\n";
    cerr << "Atom " << scaffold_atom << " is slated for deletion in the scaffold\n";
    for (int i = 0; i < scaffold_atoms_to_be_removed.number_elements(); i++)
    {
      cerr << "Atom " << scaffold_atoms_to_be_removed[i] << " in the scaffold will be removed\n";
    }
    return 0;
  }

// Do we delete our end of the bond

  if (_atoms_to_be_removed.contains(sidechain_atom) || _fragments_to_be_removed.contains(sidechain_atom))
  {
    cerr << "Sidechain_Reaction_Site::check_internal_consistency: invalid inter particle bond\n";
    cerr << "Sidechain atom " << sidechain_atom << " is slated for deletion\n";
    return 0;
  }

  return 1;
}

int
Sidechain_Reaction_Site::check_internal_consistency(const int number_sidechains,
                                const int atoms_in_scaffold_query,
                                const resizable_array<int> & scaffold_atoms_to_be_removed,
                                const resizable_array<int> & scaffold_fragments_to_be_removed)
{
  if (! Reaction_Site::check_internal_consistency())
    return 0;

  for (int i = 0; i < _inter_particle_bonds.number_elements(); i++)
  {
    const Inter_Particle_Bond * b = _inter_particle_bonds[i];

    const Matched_Atom_in_Component & a1 = b->a1();
    if (! a1.in_scaffold())
      continue;

    const Matched_Atom_in_Component & a2 = b->a2();

    if (! _check_inter_particle_atoms(atoms_in_scaffold_query, a1.matched_atom(), scaffold_atoms_to_be_removed, scaffold_fragments_to_be_removed, a2.matched_atom()))
    {
      cerr << "Sidechain_Reaction_Site:check_internal_consistency:invalid inter particle bond\n";
      return 0;
    }
  }

  for (int i = 0; i < _replace_atom.number_elements(); i++)
  {
    const Replace_Atom * b = _replace_atom[i];

    const Matched_Atom_in_Component & a1 = b->a1();
    const Matched_Atom_in_Component & a2 = b->a2();

    if (a1.in_scaffold() && a1.matched_atom() >= atoms_in_scaffold_query)
    {
      cerr << "Sidechain_Reaction_Site::check_internal_consistency: scaffold atom " << a1.matched_atom() << " is invalid - max " << atoms_in_scaffold_query << '\n';
      return 0;
    }

    if (a2.in_scaffold() && a2.matched_atom() >= atoms_in_scaffold_query)
    {
      cerr << "Sidechain_Reaction_Site::check_internal_consistency: scaffold atom " << a2.matched_atom() << " is invalid - max " << atoms_in_scaffold_query << '\n';
      return 0;
    }

//  should do the inter particle checks here... Lots of possibilities however 

  }

  return 1;
}

Reaction_Site::Reaction_Site()
{
  _3d = 0;

  _ignore_multiple_matches_involving_atoms_not_changing = 0;

  _ignore_multiple_matches_involving_changing_atoms = 0;

  _matched_atom_changed = nullptr;

  _noop_reaction = 0;

  _match_via_atom_map = 0;

  return;
}

Reaction_Site::~Reaction_Site()
{
  if (_matched_atom_changed != nullptr) {
    delete [] _matched_atom_changed;
  }

  return;
}

Scaffold_Reaction_Site::Scaffold_Reaction_Site()
{
  return;
}

Scaffold_Reaction_Site::~Scaffold_Reaction_Site()
{
  return;
}

int
Reaction_Site::ok() const
{
  if (! Substructure_Query::ok())
  {
    cerr << "Reaction site substructure query is bad\n";
    return 0;
  }

  if (_bonds_to_be_made.number_elements())
  {
    int nb = _bonds_to_be_made.number_elements();

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_made[i];
      if (! b->ok())
      {
        cerr << "Reaction_Site::ok: bad bond to be made\n";
        return 0;
      }
    }
  }

  if (_bonds_to_be_broken.number_elements())
  {
    int nb = _bonds_to_be_broken.number_elements();

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_broken[i];
      if (! b->ok())
      {
        cerr << "Reaction_Site::ok: bad bond to be broken\n";
        return 0;
      }
    }
  }

  return 1;

// All of the atoms mentioned in atom removals and fragment removals must have also
// been mentioned in a bond breaking specification

  if (_atoms_to_be_removed.number_elements() || _fragments_to_be_removed.number_elements())
  {
    if (_fragments_to_be_removed.number_elements() && _bonds_to_be_broken.empty())
    {
      cerr << "Reaction_Site::ok: fragments to be removed, but no bonds broken\n";
      return 0;
    }

    if (_atoms_to_be_removed.number_elements())
    {
      int na = _atoms_to_be_removed.number_elements();
      for (int i = 0; i < na; i++)
      {
        if (! atom_in_list_of_bonds(_bonds_to_be_broken, _atoms_to_be_removed[i]))
        {
          cerr << "Reaction_Site::ok: atom to be removed " << _atoms_to_be_removed[i] << " not in bond breaking list\n";
          return 0;
        }
      }
    }

    if (_fragments_to_be_removed.number_elements())
    {
      int na = _fragments_to_be_removed.number_elements();
      for (int i = 0; i < na; i++)
      {
        if (! atom_in_list_of_bonds(_bonds_to_be_broken, _fragments_to_be_removed[i]))
        {
          cerr << "Reaction_Site::ok: fragment to be removed " << _fragments_to_be_removed[i] << " not in bond breaking list\n";
          return 0;
        }
      }
    }
  }

  return 1;
}

static int
common_bond_creation_thingy (int a1,
                             int a2,
                             bond_type_t bt,
                             resizable_array_p<Bond> & put_here)
{
  assert (a1 >= 0 && a2 >= 0 && a1 != a2);

  Bond * b = new Bond(a1, a2, bt);

  put_here.add(b);

  return 1;
}

int
Reaction_Site::add_bond_to_be_broken(const int a1, const int a2)
{
  return common_bond_creation_thingy(a1, a2, SINGLE_BOND, _bonds_to_be_broken);
}

int
Reaction_Site::add_bond_to_be_made(const int a1, const int a2, const bond_type_t bt)
{
  return common_bond_creation_thingy(a1, a2, bt, _bonds_to_be_made);
}

int
Reaction_Site::add_toggle_kekule_form(int a1,
                                      int a2, 
                                      bond_type_t bt)
{
  return _toggle_kekule_form.add_bond(a1, a2, bt);
}

int
Reaction_Site::add_fragment_to_be_removed (int f)
{
  _fragments_to_be_removed.add(f);

  return 1;
}

int
Sidechain_Reaction_Site::add_inter_particle_bond(const int frag1, const int a1,
                                                 const int a2, const bond_type_t bt)
{
  assert (a1 >= 0 && a2 >= 0);

  Inter_Particle_Bond * b = new Inter_Particle_Bond(frag1, a1, a2, bt);

  b->is_part_of_component(_sidechain_number);

  _inter_particle_bonds.add(b);

  return 1;
}

int
Reaction_Site::add_formal_charge_to_assign(const int zatom, const formal_charge_t fc)
{
  Reaction_Formal_Charge * tmp = new Reaction_Formal_Charge;
  tmp->set_atom(zatom);
  tmp->set_formal_charge(fc);

  _formal_charges_to_assign.add(tmp);

  return 1;
}

int
Reaction_Site::add_isotope_to_be_placed(const int a, const int iso)
{
  Reaction_Place_Isotope * tmp = new Reaction_Place_Isotope(a, iso);

  _isotopes_to_assign.add(tmp);

  return 1;
}

int
Reaction_Site::add_isotope_to_be_incremented (int a, int iso)
{
  Reaction_Increment_Isotope * tmp = new Reaction_Increment_Isotope(a, iso);

  _isotopes_to_increment.add(tmp);

  return 1;
}

int
Reaction_Site::_atom_slated_for_removal (int a) const
{
  if (_atoms_to_be_removed.contains(a))
  {
    cerr << "Sidechain_Reaction_Site::_atom_slated_for_removal: matched atom " << a << " slated for atom removal\n";
    return 1;
  }

  if (_fragments_to_be_removed.contains(a))
  {
    cerr << "Sidechain_Reaction_Site::_atom_slated_for_removal: matched atom " << a << " slated for fragment removal\n";
    return 1;
  }

  return 0;
}

int
Reaction_Site::_atom_slated_for_removal (const Matched_Atom_in_Component & a) const
{
  if (a.in_scaffold())
    return 0;

  int ma = a.matched_atom();

  return _atom_slated_for_removal(ma);
}

Sidechain_Reaction_Site::Sidechain_Reaction_Site()
{
  _sidechain_number = -1;

  _make_implicit_hydrogens_explicit = 0;

  return;
}

Sidechain_Reaction_Site::~Sidechain_Reaction_Site()
{
  _sidechain_number = -3;

  return;
}

int
Sidechain_Reaction_Site::ok() const
{
  if (! Reaction_Site::ok())
    return 0;

  if (_sidechain_number < 0)
    return 0;

  int nb = _inter_particle_bonds.number_elements();
  int nr = _replace_atom.number_elements();

// It is OK to have no bonds, maybe some other sidechain is using this one

//if (0 == nb && 0 == nr)
//{
//  cerr << "Reaction_Site::ok: no inter particle bonds\n";
//  return 0;
//}

// make sure that none of the atoms in the inter particle bonds are slated for
// deletion

  for (int i = 0; i < nb; i++)
  {
    const Inter_Particle_Bond * b = _inter_particle_bonds[i];

    if (_atom_slated_for_removal(b->a2()))
    {
      cerr << "Sidechain_Reaction_Site::ok: matched atom " << b->a2() << " slated for removal\n";
      return 0;
    }
  }

  for (int i = 0; i < nr; i++)
  {
    const Replace_Atom * b = _replace_atom[i];

    const Matched_Atom_in_Component & a1 = b->a1();

    if (_atom_slated_for_removal(a1))
      return 0;
  }

  return 1;
}

int
Sidechain_Reaction_Site::debug_print (std::ostream & os, const const_IWSubstring & ind) const
{
  int nj = _inter_particle_bonds.number_elements();
  int nr = _replace_atom.number_elements();

  os << ind << "Details on Sidechain reaction site ";
  if (nj)
    os << nj << " joins ";
  if (nr)
    os << nr << " replacements ";
  else
    os << _reagents.number_elements() << " reagents";
  os << '\n';

  if (_inactive.number_elements())
    os << _inactive.number_elements() << " inactive queries\n";

  Reaction_Site::debug_print(os, ind);

  os << ind << "  will make " << nj << " inter particle bonds\n";
  for (int i = 0; i < nj; i++)
  {
    const Inter_Particle_Bond * b = _inter_particle_bonds[i];
    os << ind << "    join " << b->a1() << " to " << b->a2() << " type " << b->btype() << '\n';
  }

  os << ind << "  will make " << nr << " substitutions\n";
  for (int i = 0; i < nr; i++)
  {
    const Replace_Atom * b = _replace_atom[i];

    os << ind << "    replace " << b->a1() << " with " << b->a2() << '\n';
  }

  return os.good();
}

int
Reaction_Site::debug_print(std::ostream & os,
                           const const_IWSubstring & ind) const
{
  os << ind << "Reaction Site '" << _comment << "'\n";
  if (Substructure_Query::comment().nchars())
    os << ind << "  Query '" << Substructure_Query::comment() << "'\n";
  for (const Substructure_Query * q : _queries) {
    if (q->comment().length() > 0) {
      os << q->comment() << '\n';
    }
  }

  if (_bonds_to_be_made.number_elements())
  {
    int nb = _bonds_to_be_made.number_elements();

    os << ind << "  will make " << nb << " bonds\n";

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_made[i];
      if (b->is_single_bond())
        os << ind << "    single";
      else if (b->is_double_bond())
        os << ind << "    double";
      else if (b->is_triple_bond())
        os << ind << "    triple";

      os << " bond between " << b->a1() << " and " << b->a2() << '\n';
    }
  }

  if (_bonds_to_be_broken.number_elements())
  {
    int nb = _bonds_to_be_broken.number_elements();

    os << ind << "  will remove " << nb << " bonds\n";

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_broken[i];

      os << ind << "    break bond between " << b->a1() << " and " << b->a2() << '\n';
    }
  }

  if (_atoms_to_be_removed.number_elements())
  {
    os << ind << " will remove atoms";

    int na = _atoms_to_be_removed.number_elements();
    for (int i = 0; i < na; i++)
      os << ' ' << _atoms_to_be_removed[i];

    os << '\n';
  }

  if (_fragments_to_be_removed.number_elements())
  {
    os << ind << "  will remove fragments starting with atoms";

    int na = _fragments_to_be_removed.number_elements();
    for (int i = 0; i < na; i++)
      os << ' ' << _fragments_to_be_removed[i];

    os << '\n';
  }
  else
    os << ind << "  no fragments removed\n";

  if (_stereo_centres_to_invert.number_elements())
  {
    os << ind << " will invert stereo centres";
    for (int i = 0; i < _stereo_centres_to_invert.number_elements(); i++)
    {
      os << ' ' << _stereo_centres_to_invert[i];
    }
    os << '\n';
  }

  if (nullptr != _matched_atom_changed)
  {
    cerr << "Atom changed values\n";
    int h = highest_initial_atom_number();
    for (int i = 0; i <= h; i++)
    {
      cerr << "Atom " << i << " changed " << _matched_atom_changed[i] << '\n';
    }
  }

  return os.good();
}

static Bond *
parse_bond_attribute (const msi_attribute * att,
                      bond_type_t bt = SINGLE_BOND)
{
  if (2 != att->number_int_values())
  {
    cerr << "A bond specification must have two int values\n";
    return 0;
  }

  int a1 = att->int_multi_value(0);
  int a2 = att->int_multi_value(1);

  if (a1 == a2)
  {
    cerr << "Identical atoms found at ends of bond " << a1 << " and " << a2 << '\n';
    cerr << (*att);
    return 0;
  }

  Bond * rc = new Bond(a1, a2, bt);

  assert (rc->ok());

  return rc;
}

static int
parse_atom_values (const msi_attribute * att,
                   resizable_array<int> & zatoms)
{
  if (1 == att->number_int_values())
  {
    int tmp;
    if (! att->value(tmp) || tmp < 0)
    {
      cerr << "Atom values must be whole non-negative numbers\n";
      return 0;
    }

    zatoms.add(tmp);

    return 1;
  }

  int nv = att->number_int_values();
  for (int i = 0; i < nv; i++)
  {
    int j = att->int_multi_value(i);

    if (zatoms.contains(j))
    {
      cerr << "Duplicate specification " << j << "\n";
      return 0;
    }

    zatoms.add(j);
  }

  return 1;
}

template <typename T>
int
process_attribute (const msi_attribute * att,
                   resizable_array_p<T> & things)
{
  T * t = new T;
  if (! t->construct_from_msi_attribute(att))
  {
    cerr << "Scaffold_Reaction_Site::construct_from_msi_object: cannot parse " << (*att);
    delete t;
    return 0;
  }

  things.add(t);

  return 1;
}

template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Change_Element> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Formal_Charge> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Change_Formal_Charge> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Place_Isotope> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Increment_Isotope> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Invert_Isotope> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Stereo_Centre> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Dihedral_Angle> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Bond_Length> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_Bond_Angle> &);
template int process_attribute (const msi_attribute *, resizable_array_p<Reaction_3D_Replace> &);

int
Reaction_Site::_parse_wedge_bonds (const msi_attribute * att)
{
  Reaction_Wedge_Bond * w = new Reaction_Wedge_Bond;

  if (! w->construct_from_msi_attribute(att))
  {
    cerr << "Reaction_Site::_parse_wedge_bonds: invalid attribute " << (*att) << '\n';
    delete w;
    return 0;
  }

  _wedge_bonds_to_place.add(w);

  return 1;
}

int
Reaction_Site::_parse_recompute_implicit_hydrogens (const msi_attribute * att)
{
  int n = att->number_int_values();
  if (0 == n)
  {
    cerr << "Reaction_Site::_parse_recompute_implicit_hydrogens: no values\n";
    return 0;
  }

  if (_unfix_implicit_hydrogens.number_elements() + n < _unfix_implicit_hydrogens.elements_allocated())
    _unfix_implicit_hydrogens.resize(_unfix_implicit_hydrogens.elements_allocated() + n);

  for (int i = 0; i < n; i++)
  {
    int j = att->int_multi_value(i);
    _unfix_implicit_hydrogens.add(j);
  }

  return 1;
}

int
Reaction_Site::_read_inactive_query (const msi_attribute * att)
{
  const IWString & smarts = att->stringval();

  return _read_inactive_query(smarts);
}

int
Reaction_Site::_read_inactive_query (const IWString & smarts)
{
  if (smarts.starts_with("Q:"))     // single query file
    return _read_inactive_from_query_file(smarts);

  const_IWSubstring tmp(smarts);

  if (smarts.starts_with("F:") || smarts.starts_with("Q:"))
    return process_cmdline_token(' ', tmp, _inactive, 0);

  Substructure_Query * q = new Substructure_Query;

  if (! q->create_from_smarts(smarts))
  {
    cerr << "Reaction_Site::_read_inactive_query: invalid smarts '" << smarts << "'\n";
    delete q;
    return 0;
  }

  _inactive.add(q);

  return 1;
}

int
Reaction_Site::_read_inactive_from_query_file (const IWString & fname)
{
  IWString tmp(fname);

  assert (tmp.starts_with("Q:"));

  tmp.remove_leading_chars(2);

  Substructure_Query * q = new Substructure_Query;

  if (! q->read(tmp))
  {
    cerr << "Reaction_Site::_read_inactive_from_query_file: cannot read query file '" << fname << "'\n";
    delete q;

    return 0;
  }

  _inactive.add(q);

  return 1;
}

int
Reaction_Site::_read_query_file_from_molecule(IWString & fname)
{
  FileType input_type = discern_file_type_from_name(fname);

  if (FILE_TYPE_INVALID == input_type)
    input_type = FILE_TYPE_SMI;

  data_source_and_type<Molecule> input(input_type, fname.null_terminated_chars());

  if (! input.good())
  {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot open '" << fname << "'\n";
    return 0;
  }

  Molecule * m = input.next_molecule();

  if (nullptr == m)
  {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot read molecule from '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<Molecule> free_m(m);

  Molecule_to_Query_Specifications mqs;

  // Duplicate this until we eliminate the inheritance.
  if (! Substructure_Query::create_from_molecule(*m, mqs))
  {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot build query from " << m->smiles() << '\n';
    return 0;
  }

  Substructure_Query * q = new Substructure_Query();
  q->create_from_molecule(*m, mqs);
  _queries.add(q);

  return 1;
}

int
Reaction_Site::_read_query_file(const IWString & fname,
                                const int query_files_in_current_directory,
                                const IWString & reaction_directory)
{
  if (fname.starts_with("M:"))
  {
    IWString tmp(fname);
    tmp.remove_leading_chars(2);
    return _read_query_file_from_molecule(tmp);
  }

  IWString query_fname;

   // cerr << "Parsing query file attribute, same " << query_files_in_current_directory << " rxn dir " << reaction_directory << " name '" << fname << "'\n";
  if (query_files_in_current_directory)
    query_fname = fname;
  else if (fname.starts_with('/'))
    query_fname = fname;
  else if (fname.contains("${")) {
    // shell variables will be expanded
    query_fname = fname;
  } else if (fname.starts_with("$")) {
    // Old style shell variable expansion, starts with $NAME
    IWString temp(fname);
    expand_environment_variables(temp.c_str(),query_fname);
  }
  else if (0 == reaction_directory.length())
    query_fname = fname;
  else
  {
    query_fname = reaction_directory;
    query_fname << '/' << fname;
  }

  if (! Substructure_Query::read(query_fname))
  {
    cerr << "Reaction_Site::_read_query_file_from_molecule: cannot process query '" << query_fname << "'\n";
    cerr << "In directory " << query_files_in_current_directory << " dir '" << reaction_directory << "'\n";
    return 0;
  }

  return 1;
}

/*
  Note that we'd really like things to be independent of the order of items 
  in the input msi_object, but we really haven't achieved that. For example
  the NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE attribute may fail if
  it appears AFTER the smiles. Fix someday if it becomes a problem
*/

int
Reaction_Site::construct_from_msi_object(const msi_object & msi,
                int query_files_in_current_directory,
                const IWString & reaction_directory)
{
  int nmsi = msi.number_elements();

  int found_query = 0;    // must be specified only once

  _unique_id = msi.object_id();

  std::unordered_map<int, int> atom_map_to_unique_id;

  for (int i = 0; i < nmsi; i++)
  {
    const msi_object * obj = msi[i];

    if ("Query" == obj->name() || "query" == obj->name())
    {
      if (! Substructure_Query::construct_from_msi_object(*obj))
      {
        cerr << "Reaction_Site::construct_from_msi_object: query interpretation failed\n";
        cerr << (*obj);
        return 0;
      }

      found_query = 1;
    }
    else if (NAME_OF_NO_REACTION_OBJECT == obj->name())    // will be parsed by the owning sidechain object
      ;
    else
    {
      cerr << "Reaction_Site::construct_from_msi_object: WHAT KIND OF OBJECT IS THIS '" << obj->name() << "'\n";
      return 0;
    }
  }

  if (_match_via_atom_map)
  {
    const int h = highest_initial_atom_number();

    for (int i = 0; i <= h; ++i)
    {
      const Substructure_Atom * a = query_atom_with_initial_atom_number(i);
      if (nullptr == a)
      {
        cerr << "Reaction_Site::construct_from_msi_object:no query atom with initial atom number " << i << '\n';
        return 0;
      }

      const int amap = a->atom_map_number();
      if (amap >= 0)
        atom_map_to_unique_id[amap] = i;

    }
  }

  int one_embedding = 0;    // NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE present or not
  int unique_embeddings = 0;
  int ignore_symmetry = 0;
  int max_matches_allowed = 0;
  int embeddings_do_not_overlap = 0;

  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++)
  {
    const msi_attribute * att = msi.attribute(i);

    if (NAME_OF_QUERY_FILE_ATTRIBUTE == att->name())
    {
      if (found_query)
      {
        cerr << "Reaction_Site::construct_from_msi_object: already found query\n";
        abort();
        return 0;
      }

      if (! _read_query_file(att->stringval(), query_files_in_current_directory, reaction_directory))
      {
        cerr << "Reaction_Site::construct_from_msi_object:cannot process query file '" << att->stringval() << "'\n";
        return 0;
      }

      found_query = 1;
    }
    else if (NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE == att->name())
    {
      att->value(one_embedding);
    }
    else if (NAME_OF_UNIQUE_EMBEDDINGS_ONLY_ATTRIBUTE == att->name())
    {
      att->value(unique_embeddings);
    }
    else if (NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE == att->name())
    {
      att->value(embeddings_do_not_overlap);
    }
    else if (NAME_OF_IGNORE_SYMMETRY_RELATED_MATCHES_ATTRIBUTE == att->name())
    {
      att->value(ignore_symmetry);
    }
    else if (NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE == att->name())
    {
      if (! att->value(max_matches_allowed) || max_matches_allowed < 1)
      {
        cerr << "IWReaction::_construct_from_msi_object: the '" << NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE << "' must be a whole positive number\n";
        return 0;
      }
    }
    else if (NAME_OF_QUERY_SMARTS_ATTRIBUTE == att->name())
    {
      if (found_query)
      {
        cerr << "Reaction_Site::construct_from_msi_object: already found query\n";
        abort();
        return 0;
      }

      const_IWSubstring s;
      att->value(s);

      if (! Substructure_Query::create_from_smarts(s))
      {
        cerr << "Reaction_Site::construct_from_msi_object: cannot parse smarts '" << s << "'\n";
        return 0;
      }

      _smarts << new IWString(s);

      found_query = 1;
    }
    else if (NAME_OF_INACTIVE_ATTRIBUTE == att->name())
    {
      if (! _read_inactive_query(att))
      {
        cerr << "Reaction_Site::construct_from_msi_object: invalid " << NAME_OF_INACTIVE_ATTRIBUTE << " attribite " << (*att) << '\n';
        return 0;
      }
    }
    else if (NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE == att->name())
    {
      int tmp;
      if (! att->value(tmp))
      {
        cerr << "Reaction_Site::construct_from_msi_object: invalid " << NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE << " attribute " << (*att) << '\n';
        return 0;
      }

      Substructure_Query::set_each_component_search(tmp);
    }
    else if (NAME_OF_BOND_TO_BREAK_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att);
      if (nullptr == b)
        return 0;

      _bonds_to_be_broken.add(b);
    }
    else if (NAME_OF_SINGLE_BOND_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att);
      if (nullptr == b)
        return 0;

      _bonds_to_be_made.add(b);
    }
    else if (NAME_OF_AROMATIC_BOND_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att, AROMATIC_BOND);
      if (nullptr == b)
        return 0;

      _bonds_to_be_made.add(b);
    }
    else if (NAME_OF_DOUBLE_BOND_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att, DOUBLE_BOND);
      if (nullptr == b)
        return 0;

      _bonds_to_be_made.add(b);
    }
    else if (NAME_OF_TRIPLE_BOND_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att, TRIPLE_BOND);
      if (nullptr == b)
        return 0;

      _bonds_to_be_made.add(b);
    }
    else if (NAME_OF_BOND_TO_BREAK_ATTRIBUTE == att->name())
    {
      Bond * b = parse_bond_attribute(att);
      if (nullptr == b)
        return 0;

      _bonds_to_be_broken.add(b);
    }
    else if (NAME_OF_REMOVE_ATOM_ATTRIBUTE == att->name())
    {
      if (! parse_atom_values(att, _atoms_to_be_removed))
      {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_REMOVE_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! parse_atom_values(att, _fragments_to_be_removed))
      {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_KEEP_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! parse_atom_values(att, _fragments_to_be_kept))
      {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_RECOMPUTE_IMPLICIT_HYDROGENS_ATTRIBUTE == att->name())
    {
      if (! _parse_recompute_implicit_hydrogens(att))
      {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_WEDGE_BOND_ATTRIBUTE == att->name())
    {
      if (! _parse_wedge_bonds(att))
      {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_NOOP_REACTION == att->name())
    {
      _noop_reaction = 1;
    }
    else if (NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE == att->name())   // ignore these, handled by derived object only
      ;
    else if (NAME_OF_INTER_PARTICLE_SUBSTITUTE_ATOM_ATTRIBUTE == att->name())   // ignore these, handled by derived object only
      ;
    else if (NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE == att->name())   // ignore these, handled by derived object only
      ;
    else if (NAME_OF_SIDECHAIN_SMILES == att->name())
      ;
    else if (NAME_OF_SIDECHAIN_MOLECULE == att->name())
      ;
    else if (NAME_OF_SIDECHAIN_FILE == att->name())
      ;
    else if (NAME_OF_CHANGE_ELEMENT_ATRRIBUTE == att->name())
    {
      if (! process_attribute(att, _elements_to_change))
        return 0;
    }
    else if (NAME_OF_FORMAL_CHARGE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _formal_charges_to_assign))
        return 0;
    }
    else if (NAME_OF_CHANGE_FORMAL_CHARGE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _formal_charges_to_change))
        return 0;
    }
    else if (NAME_OF_ISOTOPE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _isotopes_to_assign))
        return 0;
    }
    else if (NAME_OF_INCREMENT_ISOTOPE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _isotopes_to_increment))
        return 0;
    }
    else if (NAME_OF_INVERT_ISOTOPE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _isotopes_to_invert))
        return 0;
    }
    else if (NAME_OF_DIHEDRAL_ANGLE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_dihedral_angle))
        return 0;

      _3d = 1;
    }
    else if (NAME_OF_BOND_LENGTH_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_bond_length))
        return 0;

      _3d = 1;
    }
    else if (NAME_OF_BOND_ANGLE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_bond_angle))
        return 0;

      _3d = 1;
    }
    else if (NAME_OF_3D_REPLACE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_3d_replace))
        return 0;

      _3d = 1;
    }
    else if (NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE == att->name())
    {
      int tmp;
      if (! att->value(tmp) || tmp < 0)
      {
        cerr << "Reaction_Site::construct_from_msi_object: invalid " << NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE << " value\n";
        return 0;
      }

      _stereo_centres_to_invert.add(tmp);
    }
    else if (NAME_OF_REACTION_COMMENT == att->name())
    {
      att->value(_comment);
    }
    else if (NAME_OF_REMOVE_CHIRAL_CENTRE_ATTRIBUTE == att->name())
    {
      for (int i = 0; i < att->number_int_values(); i++)
      {
        _chiral_centres_to_remove.add(att->int_multi_value(i));
      }
    }
    else if (NAME_OF_TOGGLE_KEKULE_FORM_ATTRIBUTE == att->name())
    {
      if (! _toggle_kekule_form.add_bond_from_msi_attribute(*att))
      {
        cerr << "Reaction_Site::construct_from_msi_object:invalid toggke kekule form specification '" << (*att) << "'\n";
        return 0;
      }
    }
    else if (NAME_OF_UNLINK_UNMATCHED_ATOMS == att->name())
    {
      if (! parse_atom_values(att, _unlink_unmatched_atoms))
      {
        cerr << "Reaction_Site::construct_from_msi_object:invalid attribute " << (*att) << '\n';
        return 0;
      }
    }
    else
    {
      cerr << "Reaction_Site::construct_from_msi_object: unrecognised attribute '" << att->name() << "'\n";
      return 0;
    }
  }

/*
  We don't want things to be dependent on the order in the msi file. Otherwise
  these directives would only work if they preceeded the smarts...
*/

  if (found_query)
  {
    if (unique_embeddings)
      set_find_unique_embeddings_only(1);
    if (ignore_symmetry)
      set_do_not_perceive_symmetry_equivalent_matches(1);
    if (one_embedding)
      set_find_one_embedding_per_atom(1);
    if (max_matches_allowed)
      set_max_matches_to_find(max_matches_allowed);
    if (embeddings_do_not_overlap)
      set_embeddings_do_not_overlap(embeddings_do_not_overlap);
      
    return 1;
  }

  if (_noop_reaction)
    return 1;

  cerr << "Reaction_Site::construct_from_msi_object: no query found in object\n";
  cerr << "Using default '[*]', which will match the first atom regardless of type\n";

  Substructure_Query::create_from_smarts("[*]");
  Substructure_Query::set_max_matches_to_find(1);

  return 1;
}

int
Scaffold_Reaction_Site::construct_from_msi_object (const msi_object & msi,
                int query_files_in_current_directory,
                const IWString & reaction_directory)
{
  if (! Reaction_Site::construct_from_msi_object(msi, query_files_in_current_directory, reaction_directory))
    return 0;

//for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++)
//{
//  _reaction_stereo_centre[i]->all_atoms_in_scaffold();
//}

  for (int i = 0; i < _reaction_dihedral_angle.number_elements(); i++)
  {
    _reaction_dihedral_angle[i]->all_atoms_in_scaffold();
  }

  for (int i = 0; i < _reaction_bond_length.number_elements(); i++)
  {
    _reaction_bond_length[i]->all_atoms_in_scaffold();
  }

  for (int i = 0; i < _reaction_bond_angle.number_elements(); i++)
  {
    _reaction_bond_angle[i]->all_atoms_in_scaffold();
  }

  return 1;
}

void
Inter_Particle_Bond::_default_values() {
  _bt = SINGLE_BOND;

  _distance = 0.0f;

  _sidechain_isotope_requirement = SidechainIsotopeRequirement::kUndefined;

  _align_3d = 0.0f;
}

Inter_Particle_Bond::Inter_Particle_Bond()
{
  _default_values();

  _bt = SINGLE_BOND;

  _distance = 0.0f;

  _sidechain_isotope_requirement = SidechainIsotopeRequirement::kUndefined;

  _align_3d = 0.0f;

  return;
}

Inter_Particle_Bond::Inter_Particle_Bond(int frag1, atom_number_t z1,
                          atom_number_t z2,
                          bond_type_t b)
{
  _default_values();

  _a1.set_in_component(frag1);
  _a1.set_matched_atom(z1);

//_a2.set_in_component(qqqq)     // may need to set this sometime
  _a2.set_matched_atom(z2);

  _bt = b;

  _distance = 0.0f;

  return;
}

// If 3D constraints are in effect, return true if the distance
// between `a1` and `a2` is less than our distance.
int
Inter_Particle_Bond::OkDistance(const Molecule& m, atom_number_t a1,
                   atom_number_t a2) const {
  if (_distance <= 0.0f) {
    return 1;
  }

  float d = m.distance_between_atoms(a1, a2);

  return d <= _distance;
}

// a1 is in the scaffold, a2 in the sidechain.
int
Inter_Particle_Bond::OkSidechainIsotopeConstraint(Molecule& m,
                atom_number_t a1,
                atom_number_t a2) const {
  if (_sidechain_isotope_requirement == SidechainIsotopeRequirement::kUndefined) {
    return 1;
  }

  const isotope_t iso2 = m.isotope(a2);

  // If no isotope on a2, nothing to check.
  if (iso2 == 0) {
    return 1;
  }

  if (_sidechain_isotope_requirement == SidechainIsotopeRequirement::kAtomicNumber) {
    return iso2 == static_cast<isotope_t>(m.atomic_number(a1));
  }

  if (_sidechain_isotope_requirement == SidechainIsotopeRequirement::kIsotope) {
    return iso2 == m.isotope(a1);
  }

  return 1;
}

int
Inter_Particle_Bond::build(const IWString & s)
{
  const_IWSubstring token;
  int i = 0;

  if (! s.nextword(token, i))
  {
    cerr << "Inter_Particle_Bond::build:cannot extract first token '" << s << "'\n";
    return 0;
  }

  if (! _a1.construct(token))
  {
    cerr << "Inter_Particle_Bond::build:invalid A1 specification '" << s << "'\n";
    return 0;
  }

  if (! s.nextword(token, i))
  {
    cerr << "Inter_Particle_Bond::build:cannot extract second token '" << s << "'\n";
    return 0;
  }

  if (! _a2.construct(token))
  {
    cerr << "Inter_Particle_Bond::build:invalid A1 specification '" << s << "'\n";
    return 0;
  }

  int b;

  _bt = SINGLE_BOND;   // the default value

  if (! s.nextword(token, i))   // nothing specified, the default is good
    ;
  else if (! token.numeric_value(b) || b <= 0 || b > 3)
  {
    cerr << "Inter_Particle_Bond::build: invalid bond type '" << s << "'\n";
    return 0;
  }
  else if (1 == b)   // same as default
    ;
  else if (2 == b)
    _bt = DOUBLE_BOND;
  else if (3 == b)
    _bt = TRIPLE_BOND;
  else
  {
    cerr << "Inter_Particle_Bond::build: what kind of bond is this '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
Inter_Particle_Bond::debug_print(std::ostream & os) const
{
  os << "Inter particle bond between " << _a1 << " and " << _a2 << '\n';
  os << "Type " << _bt << '\n';

  return os.good();
}

/*
  Inter particle bonds will have two or three components.

    (A I join (0 0))          means make a single bond between atoms 0 and 0
    (A I join (0 0 1))        same
    (A I join (1 2 2))        make a double bond between atoms 1 and 2
    (A C join "S1.1 2")       make single bond between atom 1 in sidechain 1 and matched atom 2
    (A C join "1 S.0 2")
*/

static Inter_Particle_Bond *
parse_inter_particle_bond_attribute(const IWString & s)
{
  Inter_Particle_Bond * rc = new Inter_Particle_Bond;
  if (! rc->build(s))
  {
    delete rc;
    return nullptr;
  }

  return rc;
}

static Replace_Atom *
parse_replace_atom_attribute(const msi_attribute * att,
                             int default_component)
{
  Replace_Atom * r = new Replace_Atom;
  if (! r->construct_from_msi_object(att, default_component))
  {
    delete r;
    return nullptr;
  }

  return r;
}

int
Replace_Atom::write_msi(std::ostream & os, const const_IWSubstring & ind,
                         const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"" << _a1 << ' ' << _a2 << '"' << ")\n";

  return os.good();
}

std::ostream &
operator << (std::ostream & os, const Replace_Atom & r)
{
  os << "Replace atom " << r.a1() << " and " << r.a2();

  return os;
}

/*
  A replace_atom directive must specify two matched atoms
  We get the default component to use from the object that owns us

    (A C replace_atom "S.0 2")      replace scaffold atom 0 with sidechain atom 2
    (A C replace_atom "0 4")        replace sidechain atom 0 with sidechain atom 4
*/

int
Replace_Atom::construct_from_msi_object(const msi_attribute * att,
                                        int default_component)
{
  const IWString & s = att->stringval();

  if (2 != s.nwords())
  {
    cerr << "Replace_Atom::construct_from_msi_object: must have two tokens\n";
    return 0;
  }

  const_IWSubstring s1, s2;
  (void) s.split(s1, ' ', s2);

  if (! _a1.construct(s1, default_component))
  {
    cerr << "Replace_Atom::construct_from_msi_object: invalid a1 specification '" << s1 << "'\n";
    return 0;
  }

  if (! _a2.construct(s2, default_component))
  {
    cerr << "Replace_Atom::construct_from_msi_object: invalid a2 specification '" << s2 << "'\n";
    return 0;
  }

  return 1;
}

int
Replace_Atom::adjust_matched_atoms_in_component(const extending_resizable_array<int> & xref)
{
  _a1.adjust_matched_atoms_in_component(xref);
  _a2.adjust_matched_atoms_in_component(xref);

  return 1;
}

/*int
Replace_Atom::change_so_highest_component_number_is_a2()
{
  if (_a1.in_component() < _a2.in_component())
    return 0;                    // no need to change

  Matched_Atom_in_Component tmp (_a1);
  _a1 = _a2;
  _a2 = tmp;

  return 1;
}*/

// A substructure search of `target` has been done, and the results are in
// `sresults`. If any of those embeddings touch atoms matched by any of
// the queries in `inactive_qry`, remove those embeddings from `sresults`.
static int
remove_hits_that_touch_inactive_matches (Molecule_to_Match & target,
                                     const resizable_array_p<Substructure_Query> & inactive_qry,
                                     Substructure_Results & sresults)
{
//int * inactive = new_int(target.natoms()); std::unique_ptr<int[]> free_inactive(inactive);

  std::unique_ptr<int[]> inactive_matches;  // will be allocated if needed.

//cerr << "Checking " << n << " inactive queries\n";

  for (Substructure_Query* q : inactive_qry) {
    Substructure_Results inactive_result;

    const int nhits = q->substructure_search(target, inactive_result);

    if (0 == nhits)
      continue;

    if (! inactive_matches) {
      inactive_matches.reset(new_int(target.natoms()));
    }

    inactive_result.each_embedding_set_vector(inactive_matches.get(), 1);
  }

  int nhits = sresults.number_embeddings();

  if (! inactive_matches) {
    return nhits;
  }

  for (int i = nhits - 1; i >= 0; i--)
  {
    const Set_of_Atoms * e = sresults.embedding(i);
    if (e->any_members_set_in_array(inactive_matches.get()))
      sresults.remove_embedding(i);
  }

  return sresults.number_embeddings();
}

// Variations on the substructure search interface.

int
Reaction_Site::substructure_search(Molecule& m) {
  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  return substructure_search(target, sresults);
}

int
Reaction_Site::substructure_search(Molecule& m,
                                   Substructure_Results& sresults) {
  Molecule_to_Match target(&m);
  return _perform_substructure_search(target, sresults);
}

int
Reaction_Site::substructure_search(Molecule* m,
                                   Substructure_Results& sresults) {
  Molecule_to_Match target(m);
  return _perform_substructure_search(target, sresults);
}

int
Reaction_Site::substructure_search(Molecule_to_Match& target,
                                   Substructure_Results& sresults) {
  return _perform_substructure_search(target, sresults);
}

// Search through `_queries` until a match is found to `target`.
// Embeddings are stored in `sresults` and the number of matches
//is returned.
int
Reaction_Site::_perform_substructure_search(Molecule_to_Match& target,
                        Substructure_Results& sresults) {
  if (_queries.empty()) {
    return Substructure_Query::substructure_search(target, sresults);
  }

  for (Substructure_Query * q : _queries) {
    int nhits = q->substructure_search(target, sresults);
    if (nhits) {
      return nhits;
    }
  }

  return 0;
}

int
Reaction_Site::_determine_matched_atoms_checking_inactives(Molecule & m,
                                            Substructure_Results & sresults)
{
  Molecule_to_Match target(&m);

  int nhits = _perform_substructure_search(target, sresults);

//cerr << "Nhits " << nhits << " embeddings " << sresults.number_embeddings() << '\n';

  if (0 == nhits)
    return 0;

  if (_inactive.number_elements()) {
    remove_hits_that_touch_inactive_matches(target, _inactive, sresults);
  }

  int rc = sresults.number_embeddings();

//cerr << "From substructure_search, rc = " << rc << '\n';

  if (0 == rc)
    return 0;

  if (1 == rc)    // only one hit, nothing to worry about
    ;
  else if (nullptr == _matched_atom_changed)    // not doing any filtering
    ;
  else
  {
    if (_ignore_multiple_matches_involving_changing_atoms)
      _remove_multiple_hits_that_hit_atoms_being_changed(m.natoms(), sresults);

    if (_ignore_multiple_matches_involving_atoms_not_changing)
      _remove_multiple_hits_that_do_not_involve_changing_atoms(m, sresults);
  }

  return sresults.number_embeddings();
}

int
Sidechain_Reaction_Site::remove_first_reagent_no_delete()
{
  if (_reagents.empty()) {
    cerr << "Sidechain_Reaction_Site:remove_first_reagent_no_delete:no reagents\n";
    return 0;
  }

  _reagents.remove_no_delete(0);

  return 1;
}

int
Sidechain_Reaction_Site::remove_no_delete_all_reagents()
{
  if (_reagents.empty()) {
    cerr << "Sidechain_Reaction_Site:remove_first_reagent_no_delete:no reagents\n";
    return 0;
  }

  _reagents.resize_no_delete(0);

  return 1;
}

int
Sidechain_Reaction_Site::construct_from_msi_object(const msi_object & msi,
                int query_files_in_current_directory,
                const IWString & reaction_directory,
                const Sidechain_Match_Conditions & smc)
{
  _match_conditions = smc;

  if (! Reaction_Site::construct_from_msi_object(msi, query_files_in_current_directory, reaction_directory))
    return 0;

  _copy_match_conditions_to_query();

  const msi_attribute * att;
  int i = 0;
  while (nullptr != (att = msi.attribute(NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE, i++)))
  {
    Inter_Particle_Bond * b = parse_inter_particle_bond_attribute(att->stringval());
    if (nullptr == b)
      return 0;

    b->is_part_of_component(_unique_id);

    _inter_particle_bonds.add(b);
  }


  i = 0;
  while (nullptr != (att = msi.attribute(NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE, i++)))
  {
    Replace_Atom * b = parse_replace_atom_attribute(att, _sidechain_number);
    if (nullptr == b)
      return 0;

    _replace_atom.add(b);
  }

  if (_inter_particle_bonds.empty())
  {
    if (! _noop_reaction)
      cerr << "Sidechain_Reaction_Site::construct_from_msi_object: no joins specified\n";
  }

  const msi_object * nr;
  i = 0;
  while (nullptr != (nr = msi.component(NAME_OF_NO_REACTION_OBJECT, i++)))
  {
    No_Reaction * tmp = new No_Reaction;
    if (! tmp->construct_from_msi_object(*nr))
    {
      delete tmp;
      return 0;
    }

    _no_reaction.add(tmp);
  }

  att = msi.attribute(NAME_OF_MAKE_IMPLICIT_HYDROGEN_EXPLICIT_ATTRIBUTE);
  if (nullptr != att)
    att->value (_make_implicit_hydrogens_explicit);

// Is the fragment to add present?

  att = msi.attribute(NAME_OF_SIDECHAIN_SMILES);
  if (nullptr != att)
    return _add_reagents_from_smiles(*att);

  att = msi.attribute(NAME_OF_SIDECHAIN_FILE);
  if (nullptr != att)
    return _add_reagents_from_file(*att, smc);

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagents_from_smiles(const msi_attribute & att)
{
  const_IWSubstring smiles;
  att.value(smiles);

  Molecule_and_Embedding * reagent = new Molecule_and_Embedding;

  if (! reagent->build_from_smiles(smiles))
  {
    cerr << "Sidechain_Reaction_Site::construct_from_msi_object: cannot parse '" << att << '\n';
    return 0;
  }

// Now see if the query matches

  Substructure_Results tmp;

  const int nhits = _determine_matched_atoms_checking_inactives(*reagent, tmp);
//const Substructure_Query & qq = *this;
//cerr << qq.find_one_embedding_per_start_atom() << " nhits " << nhits << '\n';


  if (1 == nhits)
    ;
  else if (0 == nhits)
  {
    if (_noop_reaction)
      ;
    else
    {
      cerr << "Sidechain_Reaction_Site::construct_from_msi_object: yipes, no match to query\n";
      cerr << att;
      return 0;
    }
  }
  else if (nhits > 1)
  {
    cerr << "Sidechain_Reaction_Site::construct_from_msi_object: yipes, " << nhits << " matches to query\n";
    return 0;
  }

  if (! _noop_reaction)
  {
    reagent->collect_matched_atoms(tmp);

    if (_toggle_kekule_form.active())
    {
      if (! reagent->do_toggle_kekule_form(_toggle_kekule_form))
        return 0;
    }
  }

  _reagents.add(reagent);

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagents_from_file(const msi_attribute & att,
                                                 const Sidechain_Match_Conditions & smc)
{
  IWString fname;
  att.value(fname);

  data_source_and_type<Molecule_and_Embedding> input(FILE_TYPE_SMI, fname);

  if (! input.good())
  {
    cerr << "Sidechain_Reaction_Site::_add_reagents_from_file:cannot open '" << fname << "'\n";
    return 0;
  }

  return add_reagents(fname.null_terminated_chars(), FILE_TYPE_SMI, smc);
}

int
Sidechain_Reaction_Site::set_single_reagent(Molecule & m)
{
  if (_reagents.number_elements())
  {
    cerr << "Sidechain_Reaction_Site::set_single_reagent:already got reagents, impossible\n";
    return 0;
  }

  if (! Substructure_Query::active() && _queries.empty())    // no queries to perform
  {
    Set_of_Atoms s;

    const int matoms = m.natoms();
    s.extend(matoms);

    for (int i = 0; i < matoms; i++)
    {
      s[i] = i;
    }

    Molecule_and_Embedding * tmp = new Molecule_and_Embedding();
    tmp->add_molecule(&m);
    tmp->set_name(m.name());
    tmp->set_embedding(s);

//  cerr << "Sidechain_Reaction_Site::set_single_reagent:added " << tmp->smiles() << '\n';

    return _reagents.add(tmp);
  }

  Molecule_to_Match target(&m);
  Substructure_Results sresults;

  const int nhits = _perform_substructure_search(target, sresults);

  if (0 == nhits) {
    cerr << "Sidechain_Reaction_Site::set_single_reagent:no hits to query in " << m.smiles() << '\n';
    return 0;
  }

  if (nhits > 1) {
    cerr << "Sidechain_Reaction_Site::set_single_reagent: " << nhits << " hits to query in " << m.smiles() << '\n';
    return 0;
  }

  Molecule_and_Embedding * tmp = new Molecule_and_Embedding;
  tmp->add_molecule(&m);
  tmp->set_name(m.name());

  tmp->set_embedding(*sresults.embedding(0));

  return _reagents.add(tmp);
}

/*
  A newly created inter particle bond must be told its component
*/

int
Inter_Particle_Bond::is_part_of_component(int c)
{
  assert (c >= 0);

  _a2.set_in_component(c);

  return 1;
}

int
Inter_Particle_Bond::adjust_matched_atoms_in_component(const extending_resizable_array<int> & xref)
{
//cerr << "Adjusting inter particle bond\n";
  _a1.adjust_matched_atoms_in_component(xref);
  _a2.adjust_matched_atoms_in_component(xref);

  return 1;
}

Molecule_and_Embedding::Molecule_and_Embedding()
{
  _owning_sidechain = nullptr;

  return;
}

Molecule_and_Embedding::Molecule_and_Embedding(const Molecule& rhs) : Molecule(rhs) {
  _owning_sidechain = nullptr;
}

int
Molecule_and_Embedding::collect_matched_atoms(const Substructure_Results & sresults,
                                              int embedding)
{
  //assert (_embedding.empty());

  const Set_of_Atoms * s = sresults.embedding(embedding);
  assert (nullptr != s);

  _embedding += *s;

  return 1;
}

int
Molecule_and_Embedding::set_embedding (const Set_of_Atoms & s)
{
  _embedding = s;

  return (_embedding.number_elements() == s.number_elements());
}

int
Molecule_and_Embedding::do_toggle_kekule_form(Toggle_Kekule_Form & tkf) {
  int changed;

  return tkf.process(*this, _embedding, changed);
}

static int
write_array (std::ostream & os,
             const const_IWSubstring & ind,
             const resizable_array<int> & things,
             const const_IWSubstring & attribute_name)
{
  int nt = things.number_elements();

  os << ind << "  (A I " << attribute_name << ' ';

  if (nt > 1)
    os << '(';

  for (int i = 0; i < nt; i++)
  {
    if (i > 0)
      os << ' ';
    os << things[i];
  }

  if (nt > 1)
    os << ')';

  os << ")\n";

  return os.good();
}

template <typename T>
int
write_set (std::ostream & os,
           const const_IWSubstring & ind,
           const resizable_array_p<T> & things,
           const const_IWSubstring & attribute_name)
{
  for (int i = 0; i < things.number_elements(); i++)
  {
    (void) things[i]->write_msi(os, ind, attribute_name);
  }

  return os.good();
}

template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Change_Element> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Formal_Charge> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Change_Formal_Charge> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Place_Isotope> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Increment_Isotope> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Invert_Isotope> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Stereo_Centre> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Wedge_Bond> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Dihedral_Angle> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Bond_Length> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_Bond_Angle> &, const const_IWSubstring &);
template int write_set (std::ostream &, const const_IWSubstring &, const resizable_array_p<Reaction_3D_Replace> &, const const_IWSubstring &);

static int
write_bond ( std::ostream & os,
            const IWString & ind,
            const char * attribute_name,
            const Bond * b)
{
  os << ind << "  (A I " << attribute_name << " (" << b->a1() << ' ' << b->a2() << "))\n";

  return os.good();
}

static int
write_atoms (std::ostream & os,
             const IWString & ind,
             const char * attribute_name,
             const resizable_array<int> & atoms)
{
  os << ind << "  (A I " << attribute_name << ' ';

  int na = atoms.number_elements();
  if (1 == na)
  {
    os << atoms[0] << ")\n";
    return os.good();
  }

  os << '(';
  for (int i = 0; i < na; i++)
  {
    if (i > 0)
      os << ' ';
    os << atoms[i];
  }

  os << "))\n";

  return os.good();
}

static char
character_representation_of_bond (bond_type_t bt)
{
  if (IS_SINGLE_BOND(bt))
    return '1';
  if (IS_DOUBLE_BOND(bt))
    return '2';
  if (IS_TRIPLE_BOND(bt))
    return '3';

  cerr << "character_representation_of_bond:what kind of a bond is " << bt << '\n';
  return '1';
}

int
Inter_Particle_Bond::write_msi(std::ostream & os,
                               const IWString & ind) const
{
  os << ind << "  (A C " << NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE << " \"" << _a1 << ' ' << _a2 << ' ' << character_representation_of_bond(_bt) << "\")\n";

  return os.good();
}

std::ostream &
operator << (std::ostream & os, const Inter_Particle_Bond & b)
{
  os << "Inter particle bond between " << b.a1() << " and " << b.a2() << " type " << character_representation_of_bond(b.btype());

  return os;
}

/*
  writing a reaction_side object is complicated by the fact that the
  underlying object will never be written as a complete object
*/

int
Reaction_Site::_write_msi(std::ostream & os,
                           int & object_id,
                           const IWString & ind)
{
  if (_comment.length())
  {
    os << ind << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";
  }

  if (! Substructure_Query::write_msi(os, object_id, ind.nchars() + 2))
    return 0;

  if (_bonds_to_be_made.number_elements())
  {
    int nb = _bonds_to_be_made.number_elements();
    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_made[i];
      if (b->is_single_bond())
        write_bond(os, ind, NAME_OF_SINGLE_BOND_ATTRIBUTE, b);
      else if (b->is_double_bond())
        write_bond(os, ind, NAME_OF_DOUBLE_BOND_ATTRIBUTE, b);
      else if (b->is_triple_bond())
        write_bond(os, ind, NAME_OF_TRIPLE_BOND_ATTRIBUTE, b);
      else
      {
        cerr << "Reaction_Site::write_msi: What kind of bond is this " << (*b) << '\n';
        return 0;
      }
    }
  }

  if (_bonds_to_be_broken.number_elements())
  {
    int nb = _bonds_to_be_broken.number_elements();
    for (int i = 0; i < nb; i++)
    {
      const Bond * b = _bonds_to_be_broken[i];
      write_bond(os, ind, NAME_OF_BOND_TO_BREAK_ATTRIBUTE, b);
    }
  }

  if (_atoms_to_be_removed.number_elements())
    write_atoms(os, ind, NAME_OF_REMOVE_ATOM_ATTRIBUTE, _atoms_to_be_removed);

  if (_fragments_to_be_removed.number_elements())
    write_atoms(os, ind, NAME_OF_REMOVE_FRAGMENT_ATTRIBUTE, _fragments_to_be_removed);

  if (_fragments_to_be_kept.number_elements())
    write_atoms(os, ind, NAME_OF_KEEP_FRAGMENT_ATTRIBUTE, _fragments_to_be_kept);

  if (_elements_to_change.number_elements())
    write_set(os, ind, _elements_to_change, NAME_OF_CHANGE_ELEMENT_ATRRIBUTE);

  if (_formal_charges_to_assign.number_elements())
    write_set(os, ind, _formal_charges_to_assign, NAME_OF_FORMAL_CHARGE_ATTRIBUTE);
                         
  if (_formal_charges_to_change.number_elements())
    write_set(os, ind, _formal_charges_to_change, NAME_OF_CHANGE_FORMAL_CHARGE_ATTRIBUTE);
                         
  if (_isotopes_to_assign.number_elements())
    write_set(os, ind, _isotopes_to_assign, NAME_OF_ISOTOPE_ATTRIBUTE);

  if (_isotopes_to_increment.number_elements())
    write_set(os, ind, _isotopes_to_increment, NAME_OF_INCREMENT_ISOTOPE_ATTRIBUTE);

  if (_isotopes_to_invert.number_elements())
    write_set(os, ind, _isotopes_to_invert, NAME_OF_INVERT_ISOTOPE_ATTRIBUTE);

  if (_unfix_implicit_hydrogens.number_elements())
    write_array(os, ind, _unfix_implicit_hydrogens, NAME_OF_RECOMPUTE_IMPLICIT_HYDROGENS_ATTRIBUTE);

  if (_wedge_bonds_to_place.number_elements())
    write_set(os, ind, _wedge_bonds_to_place, NAME_OF_WEDGE_BOND_ATTRIBUTE);

  if (_reaction_dihedral_angle.number_elements())
    write_set(os, ind, _reaction_dihedral_angle, NAME_OF_DIHEDRAL_ANGLE_ATTRIBUTE);
                         
  if (_reaction_bond_length.number_elements())
    write_set(os, ind, _reaction_bond_length, NAME_OF_BOND_LENGTH_ATTRIBUTE);

  if (_reaction_bond_angle.number_elements())
    write_set(os, ind, _reaction_bond_angle, NAME_OF_BOND_ANGLE_ATTRIBUTE);

  if (_reaction_3d_replace.number_elements())
    write_set(os, ind, _reaction_3d_replace, NAME_OF_3D_REPLACE_ATTRIBUTE);

  if (_chiral_centres_to_remove.number_elements())
    write_array(os, ind, _chiral_centres_to_remove, NAME_OF_REMOVE_CHIRAL_CENTRE_ATTRIBUTE);

  if (_unlink_unmatched_atoms.number_elements())
    write_array(os, ind, _unlink_unmatched_atoms, NAME_OF_UNLINK_UNMATCHED_ATOMS);

  for (int i = 0; i < _stereo_centres_to_invert.number_elements(); i++)
  {
    os << ind << "  (A I " << NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE << ' ' << _stereo_centres_to_invert[i] << ")\n";
  }

  for (int i = 0; i < _replace_atom.number_elements(); i++)
  {
    _replace_atom[i]->write_msi(os, ind, NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE);
  }

  return os.good();
}

int
Sidechain_Reaction_Site::write_msi(std::ostream & os,
                                    int & object_id,
                                    int indentation)
{
  assert (ok());

  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << "(" << object_id ++ << ' ' << NAME_OF_SIDECHAIN_OBJECT << '\n';

  Reaction_Site::_write_msi(os, object_id, ind);

  if (_make_implicit_hydrogens_explicit)
    os << ind << "  (" << NAME_OF_MAKE_IMPLICIT_HYDROGEN_EXPLICIT_ATTRIBUTE << ' ' << _make_implicit_hydrogens_explicit << ")\n";

  for (int i = 0; i < _inter_particle_bonds.number_elements(); i++)
  {
    _inter_particle_bonds[i]->write_msi(os, ind);
  }

  for (int i = 0; i < _no_reaction.number_elements(); i++)
  {
    _no_reaction[i]->write_msi(os, object_id, indentation + 2);
  }

  os << ind << ")\n";

  return os.good();
}

int
Scaffold_Reaction_Site::write_msi(std::ostream & os,
                                   int & object_id,
                                   int indentation)
{
  assert (ok());

  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << "(" << object_id ++ << ' ' << NAME_OF_SCAFFOLD_OBJECT << '\n';

  Reaction_Site::_write_msi(os, object_id, ind);

  os << ind << ")\n";

  return os.good();
}

int
IWReaction::ok() const
{
  if (! Scaffold_Reaction_Site::ok())
  {
    cerr << "IWReaction::ok: scaffold is bad\n";
    return 0;
  }

  int ns = _sidechains.number_elements();
  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (! s->ok())
    {
      cerr << "IWReaction::ok: sidechain " << i << " is bad\n";
      return 0;
    }
  }

  return 1;
}

int
IWReaction::debug_print(std::ostream & os) const
{
  os << "Reaction";
  if (_comment.length())
    os << " '" << _comment << "'";
  os << '\n';

  os << "Scaffold\n";
  Reaction_Site::debug_print(os, "  ");

  int ns = _sidechains.number_elements();
  os << ns << " sidechains\n";
  for (int i = 0; i < ns; i++)
  {
    os << " Sidechain " << i << '\n';
    _sidechains[i]->debug_print(os, "  ");
  }

  return os.good();
}

int
IWReaction::set_number_sidechains(int s)
{
  assert (s >= 0);

  if (_sidechains.number_elements())
    _sidechains.resize_keep_storage(0);

  _sidechains.resize(s);

  for (int i = 0; i < s; i++)
  {
    Sidechain_Reaction_Site * r = new Sidechain_Reaction_Site();
    r->set_sidechain_number(i);

    if (_make_implicit_hydrogens_explicit)
      r->set_make_implicit_hydrogens_explicit(_make_implicit_hydrogens_explicit);

    _sidechains.add(r);
  }

  return 1;
}

int
IWReaction::_construct_from_msi_object(const msi_object & msi,
                                       const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  assert (_sidechains.empty());

  int nmsi = msi.number_elements();
//if (nmsi < 2)
//{
//  cerr << "IWReaction::construct_from_msi_object: a reaction object must have at least two components\n";
//  return 0;
//}

  const msi_attribute * att = msi.attribute(NAME_OF_MATCH_VIA_ATOM_MAP);
  if (nullptr != att)               // should probably check the value
    _match_via_atom_map = 1; 

  extending_resizable_array<int> unique_id_encountered;    // must have unique object ID's

  for (int i = 0; i < nmsi; i++)
  {
    const msi_object * obj = msi[i];

    int oid = obj->object_id();
    if (0 == unique_id_encountered[oid])
      unique_id_encountered[oid] = 1;
    else
    {
      cerr << "IWReaction::_construct_from_msi_object:duplicate object id " << oid << '\n';
      return 0;
    }

    if (NAME_OF_SCAFFOLD_OBJECT == obj->name())
    {
      Reaction_Site::set_match_via_atom_map(_match_via_atom_map);

      if (! Scaffold_Reaction_Site::construct_from_msi_object(*obj, _query_files_in_current_directory, _reaction_directory))
        return 0;
    }
    else if (NAME_OF_SIDECHAIN_OBJECT == obj->name() || OLD_NAME_OF_SIDECHAIN_OBJECT == obj->name())
    {
      Sidechain_Reaction_Site * s = new Sidechain_Reaction_Site;

      s->set_match_via_atom_map(_match_via_atom_map);

      s->set_sidechain_number(_sidechains.number_elements());

      if (_make_implicit_hydrogens_explicit)
        s->set_make_implicit_hydrogens_explicit(_make_implicit_hydrogens_explicit);

      if (! s->construct_from_msi_object(*obj, _query_files_in_current_directory, _reaction_directory, sidechain_match_conditions))
        return 0;

      _sidechains.add(s);
    }
    else
    {
      cerr << "IWReaction::_construct_from_msi_object: what kind of object is this?\n";
      cerr << (*obj);
      return 0;
    }
  }

  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++)
  {
    const msi_attribute * att = msi.attribute(i);

    if (NAME_OF_REACTION_COMMENT == att->name())
      att->value(_comment);
    else if (NAME_OF_APPEND_TO_NAME == att->name())
      att->value(_append_to_name);
    else if (NAME_OF_FIND_KEKULE_FORM_FOR_INVALID_VALENCES_ATTRIBUTE == att->name())
    {
      att->value(_find_kekule_forms_for_bad_valence);
    }
    else if (NAME_OF_STEREO_CENTRE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_stereo_centre))
        return 0;
    }
    else if (NAME_OF_OPTIONAL_STEREO_CENTRE_ATTRIBUTE == att->name())
    {
      if (! process_attribute(att, _reaction_stereo_centre))
        return 0;

      _reaction_stereo_centre.last_item()->set_optional(1);
    }
    else if (NAME_OF_NOOP_REACTION == att->name())
    {
      _noop_reaction = 1;
    }
    else
    {
      cerr << "IWReaction::_construct_from_msi_object: unrecognised attribute '" << (*att) << "'\n";
      return 0;
    }
  }

  return 1;
}

int
IWReaction::construct_from_msi_object (const msi_object & msi,
                                       const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  int rc = _construct_from_msi_object(msi, sidechain_match_conditions);

  _convert_from_msi_object_numbers_to_internal_indices();

  if (rc)
    rc = check_internal_consistency();

  if (rc)
    return rc;

  cerr << msi;

  return 0;
}

/*
  We have read the msi file. Many of the matched atom directives will be in terms of
  msi object numbers. We must convert those to our internal numbering
*/

int
IWReaction::_convert_from_msi_object_numbers_to_internal_indices()
{
  extending_resizable_array<int> xref(-9);

  int u = Reaction_Site::unique_id();

  xref[u] = -1;   // special index for scaffold

  for (int i = 0; i < _sidechains.number_elements(); i++)
  {
    u = _sidechains[i]->unique_id();

    xref[u] = i;
  }

//#define DEBUG_CONVERT_FROM_MSI_OBJECT_NUMBERS_TO_INTERNAL_INDICES
#ifdef DEBUG_CONVERT_FROM_MSI_OBJECT_NUMBERS_TO_INTERNAL_INDICES
  for (int i = 0; i < xref.number_elements(); i++)
  {
    if (-9 != xref[i])
      cerr << "Object ID " << i << " becomes " << xref[i] << '\n';
  }
#endif

  Reaction_Site::adjust_matched_atoms_in_component(xref);

  Reaction_Site::set_unique_id(0);

  for (int i = 0; i < _sidechains.number_elements(); i++)
  {
    _sidechains[i]->adjust_matched_atoms_in_component(xref);

    _sidechains[i]->set_unique_id(i + 1);
  }

  return 1;
}

/*
  We have read in some specifications involving matched atoms in component data. Right now
  the component numbers refer to object numbers. Convert those to our internal numbering
*/

template <typename T>
void
do_adjust_matched_atoms_in_component(resizable_array_p<T> & items,
                                     const extending_resizable_array<int> & xref)
{
  for (int i = 0; i < items.number_elements(); i++)
  {
    items[i]->adjust_matched_atoms_in_component(xref);
  }

  return;
}

template void do_adjust_matched_atoms_in_component<Reaction_3D_Replace>(resizable_array_p<Reaction_3D_Replace>&, extending_resizable_array<int> const&);
template void do_adjust_matched_atoms_in_component<Inter_Particle_Bond>(resizable_array_p<Inter_Particle_Bond>&, extending_resizable_array<int> const&);
template void do_adjust_matched_atoms_in_component<Reaction_Dihedral_Angle>(resizable_array_p<Reaction_Dihedral_Angle>&, extending_resizable_array<int> const&);
template void do_adjust_matched_atoms_in_component<Reaction_Bond_Length>(resizable_array_p<Reaction_Bond_Length>&, extending_resizable_array<int> const&);
template void do_adjust_matched_atoms_in_component<Reaction_Bond_Angle>(resizable_array_p<Reaction_Bond_Angle>&, extending_resizable_array<int> const&);
template void do_adjust_matched_atoms_in_component<Replace_Atom>(resizable_array_p<Replace_Atom>&, extending_resizable_array<int> const&);

int
Reaction_Site::adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref)
{
  do_adjust_matched_atoms_in_component(_reaction_dihedral_angle, xref);
  do_adjust_matched_atoms_in_component(_reaction_bond_length, xref);
  do_adjust_matched_atoms_in_component(_reaction_bond_angle, xref);
  do_adjust_matched_atoms_in_component(_reaction_3d_replace, xref);
  do_adjust_matched_atoms_in_component(_replace_atom, xref);

  return 1;
}

int
Sidechain_Reaction_Site::adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref)
{
  Reaction_Site::adjust_matched_atoms_in_component(xref);

  do_adjust_matched_atoms_in_component(_inter_particle_bonds, xref);

  return 1;
}

int
IWReaction::write_msi (const char * fname)
{
  std::ofstream output(fname, std::ios::out);
  if (! output.good())
  {
    cerr << "IWReaction::write_msi:cannot open '" << fname << "'\n";
    return 0;
  } 

  return write_msi(output);
}

int
IWReaction::write_msi(std::ostream & os)
{
  int object_id = 0;

  os << "(" << object_id++ << ' ' << NAME_OF_REACTION_OBJECT << '\n';

  if (_comment.length())
    os << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";

  if (_append_to_name.length())
    os << "  (A C " << NAME_OF_APPEND_TO_NAME << " \"" << _append_to_name << "\")\n";

  if (_find_kekule_forms_for_bad_valence)
    os << "  (A I " << NAME_OF_FIND_KEKULE_FORM_FOR_INVALID_VALENCES_ATTRIBUTE << ' ' << _find_kekule_forms_for_bad_valence << ")\n";

  if (_toggle_kekule_form.active())
    _toggle_kekule_form.write_msi(os, "  ", NAME_OF_TOGGLE_KEKULE_FORM_ATTRIBUTE);

// Because the reaction stereo centres may have different attribute names, we do them explicitly
// Hmmm, they should know their attribute name...

  IWString ind("  ");

  for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++)
  {
    Reaction_Stereo_Centre * r = _reaction_stereo_centre[i];
    if (r->optional())
      r->write_msi(os, ind, NAME_OF_OPTIONAL_STEREO_CENTRE_ATTRIBUTE);
    else
      r->write_msi(os, ind, NAME_OF_STEREO_CENTRE_ATTRIBUTE);
  }

  (void) Scaffold_Reaction_Site::write_msi(os, object_id, 2);

  int nsd = _sidechains.number_elements();

  for (int i = 0; i < nsd; i++)
  {
    _sidechains[i]->write_msi(os, object_id, 2);
  }

  os << ")\n";

  return os.good();
}

int
IWReaction::do_read (iwstring_data_source & input,
                     const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  msi_object msi;

  input.set_ignore_pattern("^#");    // these are comments

  if (! msi.read(input))
  {
    cerr << "IWReaction::read: cannot read msi object\n";
    return 0;
  }

  return construct_from_msi_object(msi, sidechain_match_conditions);
}

int
IWReaction::do_read (const const_IWSubstring & fname,
                     const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "IWReaction::read: cannot open '" << fname << "'\n";
    return 0;
  }

  _reaction_directory = fname;

// Now strip down to the directory part

  assert (! _reaction_directory.ends_with('/'));    // this wouldn't make sense

  int i = _reaction_directory.rindex('/');
  if (i < 0)
    _reaction_directory = "./";
  else
    _reaction_directory.iwtruncate(i + 1);

  int rc = do_read(input, sidechain_match_conditions);

  if (0 == rc)
    cerr << "IWReaction::read: cannot read reaction from '" << fname << "'\n";

  return rc;
}

int
IWReaction::do_read (const char * fname,
                  const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  const_IWSubstring tmp = fname;

  return do_read(tmp, sidechain_match_conditions);
}

int
IWReaction::do_read (const IWString & fname,
                     const Sidechain_Match_Conditions & sidechain_match_conditions)
{
  const_IWSubstring tmp = fname;

  return do_read(tmp, sidechain_match_conditions);
}

int
Reaction_Site::do_unfix_implicit_hydrogens(Molecule & result, 
                        const Set_of_Atoms & embedding,
                        int atom_offset) const
{
  const int n = _unfix_implicit_hydrogens.number_elements();

//#define DEBUG_DO_UNFIX_IMPLICIT_HYDROGENS
#ifdef DEBUG_DO_UNFIX_IMPLICIT_HYDROGENS
  cerr << "Reaction_Site::do_unfix_implicit_hydrogens:unfixing " << n << " implicit hydrogens in " << result.smiles() << '\n';
#endif

  if (0 == n)
    return 1;

  for (int i = 0; i < n; i++)
  {
    int j = _unfix_implicit_hydrogens[i];

    if (! embedding.ok_index(j))
    {
      cerr << "Reaction_Site::do_unfix_implicit_hydrogens:invalid index " << j << " contain " << embedding.number_elements() << '\n';
      abort();
    }

    atom_number_t k = atom_offset + embedding[j];

#ifdef DEBUG_DO_UNFIX_IMPLICIT_HYDROGENS
    cerr << "Unfixing atom " << k << ' ' << result.smarts_equivalent_for_atom(k) << '\n';
#endif

    if (! result.set_implicit_hydrogens_known(k, 0))
      return 0;

    result.recompute_implicit_hydrogens(k);
  }

#ifdef DEBUG_DO_UNFIX_IMPLICIT_HYDROGENS
  cerr << "After do_unfix_implicit_hydrogens " << result.smiles() << '\n';
  result.debug_print(cerr);
#endif

  return 1;
}

int
IWReaction::_add_molecule(Molecule & result, const Molecule & added)
{
  result.add_molecule(&added);

//cerr << "IWReaction::_add_molecule:_append_names " << _append_names << " name of added " << added.name() << '\n';

  if (_append_names)
  {
    IWString tmp = result.name();
    if (0 == tmp.length())
      ;
    else if (0 == component_separator.length())
      ;
    else
      tmp += component_separator;

    tmp << added.name();
    result.set_name(tmp);
  }

  if (transfer_text_info)
    added.copy_extra_text_info_to(result);

  return 1;
}

/*
  Transfer the info from the interator into an array of Molecule_and_Embedding
*/

/*int
IWReaction::_assemble_reagents (const Reaction_Iterator & riter,
                                Molecule_and_Embedding ** reagent) const
{
  int ns = _sidechains.number_elements();

  assert (ns > 0);

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];
    
    int j;     // which reagent for this sidechain

    if (s->single_reagent_only())
      j = 0;
    else
      j = riter.reagent(i);

    reagent[i] = s->reagent(j);
  }

  return 1;
}*/

void
IWReaction::_copy_iterator_data_to_etmp(const Reaction_Iterator & iter,
                                        Enumeration_Temporaries & etmp) const
{
  int ns = _sidechains.number_elements();

//cerr << "Copying " << ns << " sidechains to etmp\n";

  const Molecule_and_Embedding ** reagent = etmp.reagent();

  for (int i = 0; i < ns; i++)
  {
    int j = iter.reagent(i);

    assert (j >= 0 && j < _sidechains[i]->number_reagents());

    reagent[i] = _sidechains[i]->reagent(j);
  }

  return;
}

/*
  We are doing a reaction where we want the first reagent from each sidechain
*/

int
IWReaction::_take_first_reagent_from_each_sidechain(Enumeration_Temporaries & etmp) const
{
  const Molecule_and_Embedding ** reagent = etmp.reagent();

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    if (0 == _sidechains[i]->number_reagents())
    {
      if (iwreaction_display_take_first_reagent_from_each_sidechain_warning_message)
        cerr << "IWReaction::_take_first_reagent_from_each_sidechain:no reagents with sidechain " << i << " in " << _comment << '\n';
      return 0;
    }

//  cerr << " i = " << i << '\n';
//  cerr << _sidechains[i] << " has " << _sidechains[i]->number_reagents() << " reagents\n";
//  cerr << _sidechains[i]->reagent(0) << '\n';

    reagent[i] = _sidechains[i]->reagent(0);
  }

  return 1;
}

int
IWReaction::remove_no_delete_all_reagents() {
  for (Sidechain_Reaction_Site* s : _sidechains) {
    s->remove_no_delete_all_reagents();
  }

  return _sidechains.number_elements();
}

/*
  Process all matches in scaffold with a given reagent
*/

int
IWReaction::perform_reaction(const Molecule * scaffold,
                             const Substructure_Results & sresults,
                             const Reaction_Iterator & iter,
                             Molecule & result)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  _copy_iterator_data_to_etmp(iter, etmp);

  _add_molecule(result, *scaffold);

  const int ne = sresults.number_embeddings();
  for (int i = 0; i < ne; i++)
  {
    if (! _perform_reaction(result, sresults.embedding(i), etmp))
      return 0;
  }

  return _do_atom_removals(result, etmp);
}

int
IWReaction::perform_reaction(const Molecule * scaffold, 
                             const Set_of_Atoms * e, 
                             const Reaction_Iterator & iter,
                             Molecule & result)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  _copy_iterator_data_to_etmp(iter, etmp);

  _add_molecule(result, *scaffold);

  if (! _perform_reaction(result, e, etmp))
    return 0;

  return _do_atom_removals(result, etmp);
}

/*
  All sidechains must own their own fragments.
  Process all query matches in SCAFFOLD with all sidechain specifications
*/

int
IWReaction::perform_reaction(const Molecule * scaffold,
                             const Substructure_Results & sresults,
                             Molecule & result)
{
  assert (0 == result.natoms());

  _add_molecule(result, * scaffold);

  int ne = sresults.number_embeddings();
  if (0 == ne)
  {
    cerr << "IWReaction::perform_reaction: no embeddings\n";
    return 0;
  }

  Enumeration_Temporaries etmp(_sidechains.number_elements());

  if (! _take_first_reagent_from_each_sidechain(etmp))
    return 0;

  for (int i = 0; i < ne; i++)     // for each scaffold embedding
  {
    if (! _perform_reaction(result, sresults.embedding(i), etmp))
      return 0;
  }

  return _do_atom_removals(result, etmp);
}

int
IWReaction::perform_reaction(const Molecule * scaffold,
                             const Set_of_Atoms * e,
                             Molecule & result)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  if (! _take_first_reagent_from_each_sidechain(etmp))
    return 0;
//  cerr << "Reaction dump:\n";
//	write_msi(cerr);
//	cerr << "********************\n";
  
  _add_molecule(result, *scaffold);

  if (! _perform_reaction(result, e, etmp))
    return 0;

  return _do_atom_removals(result, etmp);
}

template <typename T>
int
do_transformations (Molecule & result,
                    const Set_of_Atoms & embedding,
                    int offset,
                    const resizable_array_p<T> & transformations)
{
  for (int i = 0; i < transformations.number_elements(); i++)
  {
    if (! transformations[i]->process(result, embedding, offset))
    {
      cerr << "Yipes, cannot process transformation " << i << '\n';
      return 0;
    }
  }

  return 1;
}

template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Change_Element> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Formal_Charge> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Change_Formal_Charge> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Place_Isotope> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Increment_Isotope> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Invert_Isotope> &);
template int do_transformations (Molecule &, const Set_of_Atoms &, int, const resizable_array_p<Reaction_Wedge_Bond> &);

//#define DEBUG_DO_ATOM_REMOVALS

int
IWReaction::_do_atom_removals(Molecule & result,
                              Enumeration_Temporaries & etmp) const
{
  Set_of_Atoms & atoms_to_be_removed = etmp.atoms_to_be_removed();

#ifdef DEBUG_DO_ATOM_REMOVALS
  cerr << "Removing atoms " << atoms_to_be_removed << '\n';
#endif

  const int nrm = atoms_to_be_removed.number_elements();

  if (nrm)
    result.remove_atoms(atoms_to_be_removed);

  result.recompute_implicit_hydrogens();    // do this even if no atoms being removed. We may have added an atom to something whose implicit H count will need to change

  const resizable_array<const Atom *> & fragments_to_be_removed = etmp.fragments_to_be_removed();

  const int nf = fragments_to_be_removed.number_elements();

#ifdef THIS_TEST_IS_NOT_NEEDED
  if (nf > result.number_fragments())
  {
    cerr << "IWReaction::_do_atom_removals: processing " << nf << " fragment removals\n";
    cerr << "BUT, molecule contains only " << result.number_fragments() << " fragments\n";
    cerr << "This is impossible\n";
    cerr << result.smiles() << '\n';
    return 0;
  }
#endif

// Since removing one fragment is so simple, do it by itself

  if (1 == nf)
  {
    atom_number_t j = result.which_atom(fragments_to_be_removed[0]);

    if (! result.ok_atom_number(j))
    {
      cerr << "IWReaction::_do_atom_removals: very bad news, removing fragment " << 0 << " of " << nf << " fragments to be removed\n";
      cerr << "Result molecule cannot find atom in " << result.smiles() << '\n';
      iwabort();
    }

    if (! result.remove_fragment_containing_atom(j))
    {
      cerr << "IWReaction::_do_atom_removals: bad news, could not remove fragment with atom " << j << '\n';
      return 0;
    }

    return 1;
  }

// Lots of potential problems with multiple fragments. Keep track of which fragments are
// to be removed

  resizable_array<int> fragments_removed;
  fragments_removed.resize(nf);

  atoms_to_be_removed.resize_keep_storage(0);

  for (int i = 0; i < nf; i++)
  {
    atom_number_t j = result.which_atom(fragments_to_be_removed[i]);

    if (! result.ok_atom_number(j))
    {
      cerr << "IWReaction::_do_atom_removals: very bad news, removing fragment " << i << " of " << nf << " fragments to be removed\n";
      cerr << "Result molecule cannot find atom\n";
      assert (nullptr == "This is very bad");
    }

    int jfm = result.fragment_membership(j);

#ifdef DEBUG_DO_ATOM_REMOVALS
    cerr << "Fragment " << i << " is atom " << j << " in fragment " << jfm << '\n';
#endif

    if (0 == i) {
      fragments_removed.add(jfm);
    } else if (fragments_removed.contains(jfm)) {
      cerr << "IWReaction::_do_atom_removals: fragment " << jfm << " being deleted multiple times, atom " << j << ", OK\n";
      cerr << _comment << '\n';
      continue;
    } else {
      fragments_removed.add(jfm);
    }

    result.add_atoms_in_fragment(atoms_to_be_removed, jfm);
  }

  result.remove_atoms(atoms_to_be_removed);

  return 1;
}

/*
  Have only implemented a limited subset of capabilities for in_place_transformation.
  Expand as needed
*/

int
IWReaction::_in_place_transformation(Molecule & m,
                                     const Set_of_Atoms * scaffold_embedding,
                                     Enumeration_Temporaries & etmp)
{
  int * atoms_in_growing_molecule = etmp.atoms_in_growing_molecule();

  const int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    Sidechain_Reaction_Site * s = _sidechains[i];
    //cerr << "Number_reagents= " << s->number_reagents() << '\n';
      
    if (1 != s->number_reagents())   /// cannot do in place transformation
    {
      cerr << "IWReaction::_in_place_transformation:sidechain has multiple reagents, cannot do in place\n";
      return 0;
    }

    atoms_in_growing_molecule[i] = m.natoms();

    _add_molecule(m, *( s->reagent(0)));
  }

  if (! _make_inter_particle_bonds(m, scaffold_embedding, etmp))
    return 0;

  if (!DoReplacements(m, scaffold_embedding, etmp))
    return 0;

  for (int i = 0; i < ns; ++i) {
    Sidechain_Reaction_Site * s = _sidechains[i];

    const int aap = atoms_in_growing_molecule[i];

    const Molecule_and_Embedding * x = s->reagent(0);

    if (! s->do_makes_breaks(m, x->embedding(), aap, etmp))
      return 0;

    if (! s->do_unfix_implicit_hydrogens(m, *(x->embedding()), aap))
      return 0;

    if (! s->DoReplacements(m, scaffold_embedding, etmp))
      return 0;
  }

  if (! _perform_scaffold_transformations(m, scaffold_embedding, etmp)) {
    return 0;
  }

  if (! _make_cip_stereo_centres(m, *scaffold_embedding, etmp)) {
    return 0;
  }

  return 1;
}

int
IWReaction::in_place_transformation(Molecule & m,
                                    const Set_of_Atoms * scaffold_embedding)
{
  assert (nullptr != scaffold_embedding);

  if (scaffold_embedding->max_val() >= m.natoms()) {
    cerr << "IWReaction::in_place_transformation:invalid embedding " << (*scaffold_embedding) << " for reagent with " << m.natoms() << " atoms\n";
    return 0;
  }

  Enumeration_Temporaries etmp(_sidechains.number_elements());
  
  if (! _take_first_reagent_from_each_sidechain(etmp)) {
    return 0;
  }

  if (! Reaction_Site::do_makes_breaks(m, *scaffold_embedding, 0, etmp))
    return 0;
  if (! _in_place_transformation(m, scaffold_embedding, etmp))
    return 0;

  if (_append_to_name.length())
    m.append_to_name(_append_to_name);

  return _do_atom_removals(m, etmp);
}

int
IWReaction::in_place_transformation(Molecule & m,
                                    const Set_of_Atoms * scaffold_embedding,
                                    const Reaction_Iterator & iter)
{
  assert (nullptr != scaffold_embedding);

  if (scaffold_embedding->max_val() >= m.natoms())
  {
    cerr << "IWReaction::in_place_transformation:invalid embedding " << (*scaffold_embedding) << " for reagent with " << m.natoms() << " atoms\n";
    return 0;
  }

  Enumeration_Temporaries etmp(_sidechains.number_elements());

  _copy_iterator_data_to_etmp(iter, etmp);

  int * atoms_in_growing_molecule = etmp.atoms_in_growing_molecule();

  for (int i = 0; i < _sidechains.number_elements(); ++i)
  {
    int ndx = iter.reagent(i);

    const Molecule_and_Embedding * r = _sidechains[i]->reagent(ndx);

    atoms_in_growing_molecule[i] = m.natoms();         /// before we add the sidechain

    m.add_molecule(r);
  }

  cerr << "Resnot now " << m.smiles() << " append " << _append_to_name << '\n';

  if (! Reaction_Site::do_makes_breaks(m, *scaffold_embedding, 0, etmp))
    return 0;

  if (! _make_inter_particle_bonds(m, scaffold_embedding, etmp))
    return 0;

  if (! _perform_scaffold_transformations(m, scaffold_embedding, etmp))
    return 0;

  if (_append_to_name.length())
    m.append_to_name(_append_to_name);

  return _do_atom_removals(m, etmp);
}

int
IWReaction::in_place_transformations(Molecule & m,
                                     const resizable_array<const Set_of_Atoms*> & embeddings)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  if (! _take_first_reagent_from_each_sidechain(etmp))
    return 0;

  int rc = 0;

  for (const Set_of_Atoms* e : embeddings) {
    if (_toggle_kekule_form.active()) {
      int changed;
      if (! _toggle_kekule_form.process(m, *e, changed))
        return 0;
    }

    if (! Reaction_Site::do_makes_breaks(m, *e, 0, etmp)) {
      return 0;
    }

    if (! _in_place_transformation(m, e, etmp)) {
      return 0;
    }

    rc++;
  }

  if (! _do_atom_removals(m, etmp)) {
    return 0;
  }

  if (_append_to_name.length() && rc) {
    m.append_to_name(_append_to_name);
  }

  return rc;
}

int
IWReaction::in_place_transformations(Molecule & m,
                                     const Substructure_Results & sresults)
{
  Enumeration_Temporaries etmp(_sidechains.number_elements());

  if (! _take_first_reagent_from_each_sidechain(etmp))
    return 0;

  int rc = 0;

  for (uint32_t i = 0; i < sresults.number_embeddings(); i++)
  {
    const Set_of_Atoms & e = *(sresults.embedding(i));

    if (_toggle_kekule_form.active())
    {
      int changed;
      if (! _toggle_kekule_form.process(m, e, changed))
        return 0;
    }

    if (! Reaction_Site::do_makes_breaks(m, e, 0, etmp))
      return 0;

    if (! _in_place_transformation(m, sresults.embedding(i), etmp))
      return 0;

    rc++;
  }

  if (! _do_atom_removals(m, etmp))
    return 0;

  if (_append_to_name.length() && rc)
    m.append_to_name(_append_to_name);

  return rc;
}

int
IWReaction::in_place_transformations(Molecule & m)
{
  Substructure_Results sresults;

  if (0 == _determine_matched_atoms_checking_inactives(m, sresults)) 
    return 0;

//cerr << "Query got " << sresults.number_embeddings() << " embeddings\n";

  return in_place_transformations(m, sresults);
}

/*
  The last things done during a reaction
*/

int
IWReaction::_perform_scaffold_transformations(Molecule & result,
                                              const Set_of_Atoms * scaffold_embedding,
                                              const Enumeration_Temporaries & etmp) const
{

//cerr << "IWReaction::_perform_scaffold_transformations begin with " << result.smiles() << '\n';

  do_transformations(result, *scaffold_embedding, 0, _elements_to_change);
//do_transformations(result, *scaffold_embedding, 0, _formal_charges_to_change);    // done in do_makes_breaks
  do_transformations(result, *scaffold_embedding, 0, _isotopes_to_assign);
  do_transformations(result, *scaffold_embedding, 0, _isotopes_to_increment);
  do_transformations(result, *scaffold_embedding, 0, _isotopes_to_invert);

  if (_unfix_implicit_hydrogens.number_elements())
    (void) do_unfix_implicit_hydrogens(result, *scaffold_embedding, 0);

  Reaction_Site::remove_and_invert_stereo_centres(result, *scaffold_embedding, 0);

  if (! _do_3d_replacements(result, scaffold_embedding, etmp))
    return 0;

  if (! _set_bond_lengths(result, scaffold_embedding, etmp))
    return 0;

  if (! _set_bond_angles(result, scaffold_embedding, etmp))
    return 0;

  if (! _set_dihedral_angles(result, scaffold_embedding, etmp))
    return 0;

  if (! _do_replacements(result, scaffold_embedding, etmp))
    return 0;

  if (_make_implicit_hydrogens_explicit)
    result.remove_all_atoms_with_isotope(_make_implicit_hydrogens_explicit);

  if (_reaction_stereo_centre.empty())
    ;
  else if (_make_stereo_centres(result, scaffold_embedding, etmp))
    ;
  else
  {
    cerr << "IWReaction::__perform_scaffold_transformations:reaction stereo centre problem\n";
    return 0;
  }

  return 1;
}

/*
  All calls to perform_reaction ultimately come down to this.
  In order to be able to process multiple hits in the same molecule,
  we do everything but add RESULT to SCAFFOLD

  It does no atom removals, so it can be called multiple times.
  The atoms and fragments to be removed are accumulated for
  subsequent processing
*/

// #define DEBUG_PERFORM_REACTION

int
IWReaction::_perform_reaction(Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              Enumeration_Temporaries & etmp)
{
  assert (nullptr != scaffold_embedding);

  Temporarily_Disable_Messages_About_Unable_to_Compute_Implicit_Hydrogens tdmaucih;

  if (_toggle_kekule_form.active()) {
    int changed;
    if (! _toggle_kekule_form.process(result, *scaffold_embedding, changed)) {
      cerr << "IWReaction::_perform_reaction:cannot toggle Kekule form " << result.smiles() << '\n';
      cerr << "LINE " << __LINE__ << ' ' << _comment << '\n';
      return 0;
    }
  }

  int chiral_centres_present = result.chiral_centres();

  int ns = _sidechains.number_elements();
//#define DEBUG_PERFORM_REACTION 1
#ifdef DEBUG_PERFORM_REACTION
  cerr << "Before _do_intra_particle_replacements, smiles is " << result.smiles() << '\n';
  cerr << "Scaffold embedding " << (*scaffold_embedding) << '\n';
#endif

  _do_intra_particle_replacements(result, scaffold_embedding);

#ifdef DEBUG_PERFORM_REACTION
  cerr << "Before adding " << ns << " sidechains, smiles is " << result.smiles() << '\n';
  cerr << "Scaffold embedding " << (*scaffold_embedding) << '\n';
#endif

// We must do any stereo preserving substitutions before any makes or breaks.
// Make sure we do these only for single reagent sidechains (always applied) and
// whichever of the other sidechains is being used - if any

  int * atoms_in_growing_molecule = etmp.atoms_in_growing_molecule();
  const Molecule_and_Embedding ** reagent = etmp.reagent();

  for (int i = 0; i < ns; i++) {
    atoms_in_growing_molecule[i] = result.natoms();

    const Molecule_and_Embedding * r = reagent[i];

    assert (nullptr != r);

    _add_molecule(result, *r);

#ifdef DEBUG_PERFORM_REACTION
    cerr << "After adding sidechain, smiles is '" << result.smiles() << '\n';
#endif
  }

//if (! _do_replacements(result, scaffold_embedding, etmp))
//  return 0;

// Now that all stereo preserving substitutions are done, start pulling things apart

  if (! Reaction_Site::do_makes_breaks(result, *scaffold_embedding, 0, etmp)) {
    return 0;
  }

  if (_append_to_name.length())
    result.append_to_name(_append_to_name);

#ifdef DEBUG_PERFORM_REACTION
  cerr << "After makes and breaks " << result.smiles() << "\n";
#endif

  if (! _make_inter_particle_bonds(result, scaffold_embedding, etmp)) {
    return 0;
  }

  for (int i = 0; i < ns; i++)
  {
    Sidechain_Reaction_Site * s = _sidechains[i];

    int aap = atoms_in_growing_molecule[i];

    if (! s->do_makes_breaks(result, reagent[i]->embedding(), aap , etmp))
      return 0;

    if (s->has_saved_chiral_centres())
      chiral_centres_present++;

    if (! s->do_unfix_implicit_hydrogens(result, *(reagent[i]->embedding()), aap))
      return 0;
  }

#ifdef DEBUG_PERFORM_REACTION
  cerr << "Are there any chiral centres? " << chiral_centres_present << " saved? " << _saved_chiral_centre.number_elements() << '\n';
#endif

  if (0 == chiral_centres_present)
    ;
  else if (_do_restore_saved_chiral_centres(result, *scaffold_embedding, etmp))
    ;
  else
  {
    cerr << "IWReaction::perform_reaction:problems with chiral centres\n";
    return 0;
  }

#ifdef DEBUG_PERFORM_REACTION
  cerr << "Just before doing transformations " << result.smiles() << "\n";
  cerr << _reaction_stereo_centre.number_elements() << " reaction stereo centres\n";
#endif

  if (_reaction_stereo_centre.empty())
    ;
  else if (_make_stereo_centres(result, scaffold_embedding, etmp))
    ;
  else
  {
    cerr << "IWReaction::perform_reaction:reaction stereo centre problem\n";
    return 0;
  }

  if (_cip_stereo.empty()) {
  } else if (_make_cip_stereo_centres(result, *scaffold_embedding, etmp)) {
  } else {
    cerr << "IWReaction::perform_reaction:CIP reaction stereo centre problem\n";
    return 0;
  }

// Now that stereo centres have been restored, make changes to them

//Reaction_Site::remove_and_invert_stereo_centres(result, *scaffold_embedding, 0);  done in _perform_scaffold_transformations
  for (int i = 0; i < ns; i++)
  {
    int aap = atoms_in_growing_molecule[i];

    _sidechains[i]->remove_and_invert_stereo_centres(result, *(reagent[i]->embedding()), aap);
  }

  if (! _perform_scaffold_transformations(result, scaffold_embedding, etmp)) {
    return 0;
  }

#ifdef DEBUG_PERFORM_REACTION
  cerr << "After performing scaffold transformations " << result.smiles() << '\n';
#endif

  if (_find_kekule_forms_for_bad_valence && ! result.valence_ok())
    _do_find_kekule_forms_for_bad_valence(result);

  return 1;
}

/*
  Only the sidechains with a single reagent will be applied
*/

/*int
IWReaction::_perform_reaction (Molecule & result,
                               const Set_of_Atoms * scaffold_embedding,
                               const Molecule * sidechain,
                               const Set_of_Atoms * sidechain_embedding)
{
  int atoms_already_present = result.natoms();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize (atoms_already_present);

  resizable_array<const Atom *> fragments_to_be_removed;
  fragments_to_be_removed.resize (atoms_already_present);

  if (! _perform_reaction (result, nullptr, scaffold_embedding, sidechain, sidechain_embedding,
                           atoms_to_be_removed, fragments_to_be_removed))
    return 0;

// All the transformations except atom removals have been done

  return _do_atom_removals(result, etmp);
}*/

/*
  We frequently need to know whether or not a matched atom number is consistent with an embedding
*/

/*static int
valid_member_of_embedding (const Set_of_Atoms & embedding,
                           int i,
                           const char * caller)
{
  if (embedding.ok_index(i))
    return 1;

  cerr << caller << " invalid index. Embedding contains " << embedding.number_elements() << " items, so " << i << " is invalid\n";
  return 0;
}*/

static int
valid_member_of_embedding(const Set_of_Atoms & embedding,
                          int i,
                          const char * caller,
                          atom_number_t & zatom)
{
  if (embedding.ok_index(i))
  {
    zatom = embedding[i];
    return 1;
  }

  if (i < 0)
  {
    int j = embedding.number_elements() + i;
    if (embedding.ok_index(j))
    {
//    cerr << "Matched atom " << i << " becomes atom " << j << " N = " << embedding.number_elements() << " atom " << embedding[j] << '\n';
      zatom = embedding[j];
      return 1;
    }
  }

  cerr << caller << " invalid index. Embedding contains " << embedding.number_elements() << " items, so " << i << " is invalid\n";
  return 0;
}

//#define DEBUG_DO_MAKES_AND_BREAKS

/*
  Well, this is quite complicated unfortunately.
  We have a function for deleting a series of atoms from a molecule,
  so that part is easy. The fragment deletion is more difficult,
  because we don't know how many, or which atoms will be removed
  during any fragment removal. For that reason, we store Atom
  pointers in an array, and use that to keep track of which
  atoms define the fragments to be removed. Note that this is
  only really necessary if there is more than one fragment to
  be deleted, but let's make this general

  EMBEDDING is the embedding in the sidechain which has been added.
  If we are a single reagent then our reagent should have been
  the one added and on input, EMBEDDING should be nullptr;
*/

int
Reaction_Site::do_makes_breaks(Molecule & result,
                               const Set_of_Atoms & embedding,
                               const int offset,
                               Enumeration_Temporaries & etmp)
{
#ifdef DEBUG_DO_MAKES_AND_BREAKS
  cerr << "do_makes_breaks:begin " << result.smiles() << ' ' << result.name() << '\n';
#endif

  const int chiral_centres_present = result.chiral_centres();

  if (chiral_centres_present)
  {
    _saved_chiral_centre.resize_keep_storage(0);
    _saved_chiral_centre.resize(chiral_centres_present);
  }

  int nb = _bonds_to_be_broken.number_elements();

#ifdef DEBUG_DO_MAKES_AND_BREAKS
  cerr << "Reaction_Site::do_makes_breaks: doing " << nb << " bond breaks, embedding " << embedding << '\n';
  cerr << result.smiles() << ' ' << result.natoms() << " atoms, offset " << offset << '\n';
#endif

// Make sure we never remove the same chiral centre twice!

  int * chiral_centre_available = new_int(result.natoms()); std::unique_ptr<int[]> free_chiral_centre_available(chiral_centre_available);

  if (chiral_centres_present)
  {
    result.at_centre_of_chiral_centre(chiral_centre_available, 1);

    if (_atoms_to_be_removed.number_elements() > 0)
      _discern_chiral_centres_to_be_saved_around_removed_atoms(result, embedding, offset, chiral_centre_available);
  }

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bonds_to_be_broken[i];

    assert (embedding.ok_index(b->a1()));
    assert (embedding.ok_index(b->a2()));

#ifdef DEBUG_DO_MAKES_AND_BREAKS
//  cerr << "getting items " << b->a1() << " and " << b->a2() << " from " << embedding << '\n';
#endif

    const atom_number_t a1 = embedding[b->a1()] + offset;
    const atom_number_t a2 = embedding[b->a2()] + offset;

//  cerr << a1 << " and " << a2 << '\n';

//  Save any chiral centres.

    if (chiral_centres_present)
    {
      if (1 == chiral_centre_available[a1])
      {
        chiral_centre_available[a1] = 0;
        _saved_chiral_centre.add(result.remove_no_delete_chiral_centre_at_atom(a1));
#ifdef DEBUG_DO_MAKES_AND_BREAKS
        cerr << "Saving chiral centre\n";
        _saved_chiral_centre.last_item()->debug_print(cerr);
#endif
      }

      if (1 == chiral_centre_available[a2])
      {
        chiral_centre_available[a2] = 0;
        _saved_chiral_centre.add(result.remove_no_delete_chiral_centre_at_atom(a2));
#ifdef DEBUG_DO_MAKES_AND_BREAKS
        cerr << "Saving chiral centre\n";
        _saved_chiral_centre.last_item()->debug_print(cerr);
#endif
      }
    }

#ifdef DEBUG_DO_MAKES_AND_BREAKS
    cerr << "Atoms " << a1 << " and " << a2 << '\n';
    cerr << result.smiles() << '\n';
#endif

//  Oct 99. Ran into a case where the query matched multiple times across the
//  same bond - and that bond was removed. The second time we tried to remove
//  the bond we got a fatal error. Therefore only try the removal if the atoms are bonded

    if (result.are_bonded(a1, a2)) {
      result.remove_bond_between_atoms(a1, a2);
    } else if (warn_non_bonded_breaks) {
      cerr << "Reaction_Site::do_makes_breaks: atoms " << a1 << " and " << a2 << " are not bonded. IGNORED. " << _comment << '\n';
      cerr << result.smiles() << '\n';
      cerr << result.smarts_equivalent_for_atom(a1) << ' ' << result.smarts_equivalent_for_atom(a2) << '\n';
      // abort();
    }
  }

#ifdef PRINT_SAVED_CHIRAL_CENTRES
  cerr << "Found " << _saved_chiral_centre.number_elements() << " chiral centres to be saved\n";
  if (_saved_chiral_centre.number_elements())
    _saved_chiral_centre[0]->debug_print(cerr);
#endif

  nb = _bonds_to_be_made.number_elements();

#ifdef DEBUG_DO_MAKES_AND_BREAKS
  cerr << "do_makes_breaks:making " << nb << " bonds to be made\n";
#endif

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bonds_to_be_made[i];

    atom_number_t a1;
    if (! valid_member_of_embedding(embedding, b->a1(), "Reaction_Site::do_makes_breaks: a1", a1))
      return 0;

    atom_number_t a2;
    if (! valid_member_of_embedding(embedding, b->a2(), "Reaction_Site::do_makes_breaks: a2", a2))
      return 0;

    a1 += offset;
    a2 += offset;

//  Ran into problems while making a double bond to an atom which had a chiral
//  centre. set_bond_type_between_atoms removes the chiral centre, but we must
//  now tell that atom that it's implicit hydrogens are no longer known

    result.set_implicit_hydrogens_known(a1, 0);
    result.set_implicit_hydrogens_known(a2, 0);

    if (! result.are_bonded(a1, a2))    // new bond - could enable the permanent aromatic thing if I wanted to
      result.add_bond(a1, a2, b->btype());
    else               // change bond
    {
      if (AROMATIC_BOND == b->btype())
      {
        Bond * b = const_cast<Bond *>(result.bond_between_atoms(a1, a2));
        b->set_permanent_aromatic(1);
      }
      else
        result.set_bond_type_between_atoms(a1, a2, b->btype());
    }

    result.recompute_implicit_hydrogens(a1);   // very important to recompute as soon as formed
    result.recompute_implicit_hydrogens(a2);
  }

#ifdef DEBUG_DO_MAKES_AND_BREAKS
  if (nb)
    cerr << "After making " << nb << " bonds '" << result.smiles() << "'\n";
#endif

// We don't actually remove any atoms (atom numbering may yet change), we just
// add things to the arrays of atoms/fragments awaiting deletion

  Set_of_Atoms & atoms_to_be_removed = etmp.atoms_to_be_removed();

  int na = _atoms_to_be_removed.number_elements();

  for (int i = 0; i < na; i++)
  {
//  cerr << "Embedding contains " << embedding.number_elements() << " items, remove atom " << _atoms_to_be_removed[i] << '\n';
    assert (embedding.ok_index(_atoms_to_be_removed[i]));

    const atom_number_t a = embedding[_atoms_to_be_removed[i]] + offset;

//  atoms_to_be_removed.add(a);
    if (atoms_to_be_removed.add_if_not_already_present(a))
      result.remove_bonds_to_atom(a);   // added Aug 2003
    else if (warn_duplicate_atom_deletion) {
      cerr << "Reaction_Site::do_makes_breaks:ignoring duplicate deletion, atom " << a << " in '" << result.name() << "'\n";
    }
  }

// The order in which we do the fragment related stuff is important

  int nf = result.number_fragments();

  if (nf > 1)
  {
    resizable_array<const Atom *> & fragments_to_be_removed = etmp.fragments_to_be_removed();

    na = _fragments_to_be_removed.number_elements();
    for (int i = 0; i < na; i++)
    {
      atom_number_t j = embedding[_fragments_to_be_removed[i]] + offset;

//    cerr << "Fragment being removed anchored at atom " << j << ' ' << result.smarts_equivalent_for_atom(j) << '\n';

      fragments_to_be_removed.add(result.atomi(j));
    }
  }

  na = _unlink_unmatched_atoms.number_elements();
  for (int i = 0; i < na; i++)
  {
    int j = _unlink_unmatched_atoms[i];
    atom_number_t a = embedding[j] + offset;

    int acon = result.ncon(a);
    Set_of_Atoms remove_bonds_to;

    for (int k = 0; k < acon; k++)
    {
      atom_number_t l = result.other(a, k);
      if (embedding.contains(l - offset))
        continue;

      remove_bonds_to.add(l);
    }

    for (int k = 0; k < remove_bonds_to.number_elements(); k++)
    {
      result.remove_bond_between_atoms(a, remove_bonds_to[k]);
    }

    nf = result.number_fragments();
  }

  na = _fragments_to_be_kept.number_elements();

  if (na > 0 && nf > 1)
  {
    resizable_array<const Atom *> & fragments_to_be_removed = etmp.fragments_to_be_removed();

    int * keep_frag = new_int(nf); std::unique_ptr<int[]> free_keep_frag(keep_frag);

    for (int i = 0; i < na; i++)
    {
      atom_number_t j = embedding[_fragments_to_be_kept[i]] + offset;

      int f = result.fragment_membership(j);

      keep_frag[f] = 1;
    }

    for (int i = 0; i < nf; i++)
    {
      if (keep_frag[i])
        continue;

      atom_number_t j = result.first_atom_in_fragment(i);

      fragments_to_be_removed.add_if_not_already_present(result.atomi(j));
    }
  }

  do_transformations(result, embedding, offset, _elements_to_change);
  do_transformations(result, embedding, offset, _formal_charges_to_assign);
  do_transformations(result, embedding, offset, _formal_charges_to_change);
  do_transformations(result, embedding, offset, _isotopes_to_assign);
//do_transformations(result, embedding, offset, _isotopes_to_increment);
  do_transformations(result, embedding, offset, _wedge_bonds_to_place);

  do_cahn_ingold_prelog_stereo(result, embedding, offset);

  return 1;
}

int
Reaction_Site::do_cahn_ingold_prelog_stereo(Molecule& result,
        const Set_of_Atoms& embedding,
        const int offset) const {
  int rc = 1;
  for (const SiteCipStereo* cips : _cip_stereo) {
    if (! do_cahn_ingold_prelog_stereo(result, embedding, offset, *cips)) {
      cerr << "Reaction_Site::do_cahn_ingold_prelog_stereo failed\n";
      rc = 0;
    }
  }

  return rc;
}

// Create a Chiral_Center at atom `zatom` and ensure that its CIP stereo
// value matches `rs`
static int
CreateChiralCentreMatchingCip(Molecule& m,
                atom_number_t zatom,
                CahnIngoldPrelog rs) {
  Chiral_Centre* c = m.create_chiral_centre(zatom);
  if (c == nullptr) {
    cerr << "CreateChiralCentreMatchingCip:chiral centre not created atom " << zatom << '\n';
    return 0;
  }

  // Note that we cannot just call CahnIngoldPrelog(c) because it is
  // changed by the invert_chirality_on_atom call - the CIP value
  // is in fact unchanged according to that measure. That would be
  // more efficient, this is clearer. Change if ever a problem.
  std::optional<CahnIngoldPrelog> cip = m.CahnIngoldPrelogValue(zatom);
  if (! cip) {
    return 0;
  }
  if (cip == CahnIngoldPrelog::kNeither) {
    return 0;
  }
  if (*cip == rs) {
    return 1;
  }
  m.invert_chirality_on_atom(zatom);
  cip = m.CahnIngoldPrelogValue(zatom);
  if (*cip == rs) {
    return 1;
  }

  cerr << "CreateChiralCentreMatchingCip:no match after switching chirality " << rs << '\n';
  cerr << "Current " << m.CahnIngoldPrelogValue(c) << '\n';
  return 0;
}

int
Reaction_Site::do_cahn_ingold_prelog_stereo(Molecule& result,
        const Set_of_Atoms& embedding,
        const int offset,
        const SiteCipStereo& cips) const {
  if (! embedding.ok_index(cips.atom())) {
    cerr << "Reaction_Site::do_cahn_ingold_prelog_stereo:invalid embedding member " <<
        cips.atom() << " in " << embedding << '\n';
    return 0;
  }

  atom_number_t a = embedding[cips.atom()];
  if (! result.ok_atom_number(a)) {
    cerr << "Reaction_Site::do_cahn_ingold_prelog_stereo:invalid atom number in result " <<
        cips.atom() << " in " << embedding << " result has " << result.natoms() << " atoms\n";
    return 0;
  }

  Set_of_Atoms connections;
  result.connections(a, connections);
  if (connections.size() < 3) {
    cerr << "Reaction_Site::do_cahn_ingold_prelog_stereo:insufficient connections " <<
        result.smarts_equivalent_for_atom(a) << '\n';
    return 0;
  }
  if (connections.size() > 4) {
    cerr << "Reaction_Site::do_cahn_ingold_prelog_stereo:excessive connections " <<
        result.smarts_equivalent_for_atom(a) << '\n';
    return 0;
  }

  return CreateChiralCentreMatchingCip(result, a, cips.cip());
}

int
Reaction_Site::remove_and_invert_stereo_centres(Molecule & result,
                                const Set_of_Atoms & embedding,
                                const int offset) const
{
  if (_stereo_centres_to_invert.number_elements())
    _do_invert_stereo_centres(result, embedding, offset);

  if (_chiral_centres_to_remove.number_elements())
    _do_remove_stereo_centres(result, embedding, offset);

  return 1;
}

/*
  Some atoms are going to be removed. Do we need to save any chiral centres that may be
  bonded
*/

int
Reaction_Site::_discern_chiral_centres_to_be_saved_around_removed_atoms(Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset,
                                          int * chiral_centre_available)
{
  const int n = _atoms_to_be_removed.number_elements();

//cerr << "Reaction_Site::_discern_chiral_centres_to_be_saved_around_removed_atoms:checking " << n << " atoms to be remove, offset " << offset << '\n';

  for (int i = 0; i < n; i++)
  {
    int j = _atoms_to_be_removed[i];

    atom_number_t k = embedding[j] + offset;

    assert (k >= 0 && k < result.natoms());

    const Atom * ak = result.atomi(k);

    int kcon = ak->ncon();

//  cerr << "Embedding member " << j << " becomes atom " << k << " ncon = " << kcon << '\n';

    for (int l = 0; l < kcon; l++)
    {
      atom_number_t j2 = ak->other(k, l);

      if (0 == chiral_centre_available[j2])
        continue;

      _saved_chiral_centre.add(result.remove_no_delete_chiral_centre_at_atom(j2));

#ifdef PRINT_SAVED_CHIRAL_CENTRES
      cerr << "Identified chiral centre to be removed\n";
      _saved_chiral_centre.last_item()->debug_print(cerr);
#endif

      chiral_centre_available[j2] = 0;
    }
  }

  return 1;
}

int
Reaction_Site::_do_invert_stereo_centres (Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset) const
{
//cerr << "Reaction_Site::_do_invert_stereo_centres: inverting " << _stereo_centres_to_invert.number_elements() << " atoms\n";

  if (0 == result.chiral_centres())    // inversions are always optional
    cerr << "No chiral centres to invert\n";
  if (0 == result.chiral_centres())    // inversions are always optional
    return 1;

  for (int i = 0; i < _stereo_centres_to_invert.number_elements(); i++)
  {
    int m = _stereo_centres_to_invert[i];

    atom_number_t a = embedding[m] + offset;

//  cerr << "Embedding member " << m << " is atom " << a << " present " << result.chiral_centre_at_atom(a) << " offset = " << offset << '\n';

    if (nullptr != result.chiral_centre_at_atom(a))
      result.invert_chirality_on_atom(a);
  }

  return 1;
}

int
Reaction_Site::_do_remove_stereo_centres (Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset) const
{
  if (0 == result.chiral_centres())    // none to remove, we are dong
    return 1;

  for (int i = 0; i < _chiral_centres_to_remove.number_elements(); i++)
  {
    int m = _chiral_centres_to_remove[i];

    atom_number_t a = embedding[m] + offset;

    if (nullptr != result.chiral_centre_at_atom(a))
      result.remove_chiral_centre_at_atom(a);
  }

  return 1;
}

int
Sidechain_Reaction_Site::add_reagents(const char * fname,
                                      FileType input_type,
                                      const Sidechain_Match_Conditions & smc)
{
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    if (0 == input_type)
    {
      cerr << "Sidechain_Reaction_Site::add_reagents: cannot discern input type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule_and_Embedding> input (input_type, fname);
  if (! input.ok())
  {
    cerr << "Sidechain_Reaction_Site::add_reagents: cannot open '" << fname << "'\n";
    return 0;
  }

  return add_reagents(input, smc);
}

int
Sidechain_Reaction_Site::_copy_match_conditions_to_query()
{
  if (_match_conditions.one_embedding_per_start_atom())
    Reaction_Site::set_find_one_embedding_per_atom(1);

  if (_match_conditions.ignore_symmetry_related_matches())
    Reaction_Site::set_do_not_perceive_symmetry_equivalent_matches(1);

  if (_match_conditions.find_unique_embeddings_only())
    Reaction_Site::set_find_unique_embeddings_only(1);

  int p = _match_conditions.process_hit_number();
  if (p >= 0)
    Reaction_Site::set_max_matches_to_find(p + 1);

  return 1;
}


int
Sidechain_Reaction_Site::add_reagents(data_source_and_type<Molecule_and_Embedding> & input,
                                      const Sidechain_Match_Conditions & smc)
{
  _match_conditions = smc;

  int nm = input.molecules_remaining();
  if (0 == nm)
  {
    cerr << "Sidechain_Reaction_Site::add_reagents: no molecules in file\n";
    return 0;
  }

  _reagents.resize(nm);

  Reaction_Site::set_do_not_perceive_symmetry_equivalent_matches(1);   // let's just make this the default

  _copy_match_conditions_to_query();

  Make_Implicit_Hydrogens_Explicit mihe;
  if (_make_implicit_hydrogens_explicit)
    mihe.set_isotope(_make_implicit_hydrogens_explicit);

  Molecule_and_Embedding * m;
  while (nullptr != (m = input.next_molecule()))
  {
    if (_make_implicit_hydrogens_explicit)
    {
      mihe.reset();
      m->make_implicit_hydrogens_explicit(mihe);
    }

//  cerr << "After making hydrogens explicit " << m->smiles() << "'\n";

    if (! add_reagent(m, smc))
    {
      cerr << "Sidechain_Reaction_Site::add_reagents: fatal error processing '" << m->name() << "'\n";
      delete m;
      return 0;
    }
  }

  return _reagents.number_elements();
}

int
Sidechain_Reaction_Site::add_reagent(Molecule_and_Embedding * m,
                                     const Sidechain_Match_Conditions & smc)
{
  if (smc.strip_reagents_to_largest_fragment()) {
    m->reduce_to_largest_fragment();
  }

  Substructure_Results sresults;

  int nhits = substructure_search(*m, sresults);  
  int reagent_is_ok = 1;

  if (0 == nhits) {     // cannot be used by this sidechain
    MaybeIssueNoHitsWarning(*m, sresults, smc);
    reagent_is_ok = 0;
  } else if (nhits > smc.suppress_if_more_than_this_many_substructure_search_hits()) {
    cerr << "Sidechain_Reaction_Site::add_reagent::too many hits " << nhits << " (max " << smc.suppress_if_more_than_this_many_substructure_search_hits() << ") for reagent " << m->name() << '\n';
    reagent_is_ok = 0;
  }

  if (! reagent_is_ok) {    
    if (stream_for_sidechains_not_matching_query.active())
      stream_for_sidechains_not_matching_query.write(*m);

    if (smc.ignore_not_reacting()) {
      delete m;
      return 1;
    }

    return 0;
  }

  //cerr << "Sidechain_Reaction_Site::add_reagent:_matched_atom_changed? " << _matched_atom_changed << '\n';

  if (nullptr != _matched_atom_changed)
  {
    _remove_multiple_hits_that_do_not_involve_changing_atoms(*m, sresults);
    nhits = sresults.number_embeddings();
  }

  if (nhits > 1 && smc.verbose() )
  {
    if (_comment.length())
      cerr << _comment;
    else
    {
      cerr << "Sidechain_Reaction_Site::add_reagent";
      cerr << ", " << nhits << " hits to reagent '" << m->name() << "'\n";
    }
  }

// cerr << "What about Kekule shifts " << _toggle_kekule_form.active() << '\n';;

  if (1 == nhits)    // hopefully the most common case
  {  
    m->collect_matched_atoms(sresults);

    if (_toggle_kekule_form.active() && ! m->do_toggle_kekule_form(_toggle_kekule_form)) {
        return 0;
    }

    _reagents.add(m);

    return 1;     // great, found a sidechain which can use this reagent
  }

  if (smc.process_hit_number() >= 0)
  {    
    if (nhits < smc.process_hit_number() + 1)     // cannot be used by this sidechain
    {
        cerr << "Sidechain_Reaction_Site::add_reagent: ";
        if (_comment.length())
          cerr << "'" << _comment << "' ";
        cerr << "hit number to process (" << smc.process_hit_number() << ") is out of range - #hits= " << nhits << '\n';

        return 0;
    }    
    m->collect_matched_atoms(sresults, smc.process_hit_number());

    if (_toggle_kekule_form.active() && ! m->do_toggle_kekule_form(_toggle_kekule_form))
        return 0;

    if (smc.multiple_match_string().length()) {
      IWString tmp = m->name();
      tmp.append_with_spacer(smc.multiple_match_string());

      m->set_name(tmp);
    }

    _reagents.add(m);

    return 1;
  }

// If we are making regio-isomers, make NHITS - 1 copies of the reagent

//cerr << "make_new_reagent_for_each_hit? " <<  (smc.make_new_reagent_for_each_hit()) << '\n';

  if (smc.make_new_reagent_for_each_hit())    
  {
    m->collect_matched_atoms(sresults);

    if (_toggle_kekule_form.active() && !m->do_toggle_kekule_form(_toggle_kekule_form))
      return 0;

    _reagents.add(m);

    for (int i = 1; i < nhits; i++)
    {
      Molecule_and_Embedding * mcopy = new Molecule_and_Embedding;
      mcopy->add_molecule(m);
      mcopy->set_name(m->name());
      mcopy->collect_matched_atoms(sresults, i);
      if (_toggle_kekule_form.active() && ! mcopy->do_toggle_kekule_form(_toggle_kekule_form))
        return 0;

      _reagents.add(mcopy);
    }

    return 1;
  }

// Looks like multiple hits, but we have no instructions on how to deal with that

  cerr << "Sidechain_Reaction_Site::add_reagent: " << nhits << " hits to query, fatal\n";
  cerr << "Molecule " << m->name() << '\n';

#ifdef ECHO_MULTIPLE_MATCHES
  for (int i = 0; i < sresults.number_embeddings(); i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);
    cerr << (*e) << '\n';
  }
#endif

  return 0;
}

void
Sidechain_Reaction_Site::MaybeIssueNoHitsWarning(const Molecule_and_Embedding& m,
                const Substructure_Results& sresults,
                const Sidechain_Match_Conditions& smc) const {
  if (! smc.issue_sidechain_no_match_warnings()) {
    return;
  }

  cerr << "Sidechain_Reaction_Site::add_reagent: ";
  if (_comment.length()) {
    cerr << "'" << _comment << "' ";
  }
  cerr << "no hits for reagent '" << m.name() << "', only matched " <<
              sresults.max_query_atoms_matched_in_search() << " query atoms\n";
}

int
Sidechain_Reaction_Site::add_reagent(const Molecule& m,
                                     const Sidechain_Match_Conditions & smc) {
  std::unique_ptr<Molecule_and_Embedding> tmp = std::make_unique<Molecule_and_Embedding>(m);
  if (add_reagent(tmp.get(), smc)) {
    tmp.release();
    return 1;
  }

  return 0;
}

int
Sidechain_Reaction_Site::add_reagent_embedding_identified (Molecule_and_Embedding * m,
                                const Sidechain_Match_Conditions & smc)
{
  if (_toggle_kekule_form.active() && ! m->do_toggle_kekule_form(_toggle_kekule_form))
      return 0;

  _reagents.add(m);

  return 1;
}

int
Sidechain_Reaction_Site::empty_reagents_array()
{
  _reagents.resize_keep_storage(0);

  return 1;
}

int
Sidechain_Reaction_Site::remove_last_reagent() {
  if (_reagents.empty()) {
    return 0;
  }

  _reagents.chop();

  return 1;
}

int
Sidechain_Reaction_Site::do_makes_breaks (Molecule & result,
                                const Set_of_Atoms * embedding,
                                int offset,
                                Enumeration_Temporaries & etmp)
{
  return Reaction_Site::do_makes_breaks(result, *embedding, offset, etmp);
}

// #define DEBUG_MAKE_INTER_PARTICLE_BONDS

int
determine_atom_number(const Set_of_Atoms & scaffold_embedding,
                      const Matched_Atom_in_Component & ma,
                      const Enumeration_Temporaries & etmp,
                      const char * caller,
                      atom_number_t & zresult)
{
  if (ma.in_scaffold())
    return valid_member_of_embedding(scaffold_embedding, ma.matched_atom(), caller, zresult);

  int c = ma.in_component();

  assert (c >= 0);

//const Molecule_and_Embedding * sc = etmp.reagent()[c];

#ifdef DEBUG_MAKE_INTER_PARTICLE_BONDS
  cerr << "The atom is in component " << c << '\n';
  cerr << "Reagent contains " << etmp.reagent()[c]->natoms() << " atoms\n";
#endif

  assert (etmp.reagent()[c]->ok());

  const Set_of_Atoms * e = etmp.reagent()[c]->embedding();

#ifdef DEBUG_MAKE_INTER_PARTICLE_BONDS
  cerr << "Embedding contains " << e->number_elements() << " atoms\n";
#endif

  if (! valid_member_of_embedding(*e, ma.matched_atom(), "Sidechain_Reaction_Site::make_inter_particle_bonds: a1", zresult))
    return 0;

  assert (zresult >= 0);

#ifdef DEBUG_DETERMINE_ATOM_NUMBER
  cerr << "Growing molecule " << c << " had " << aigm[c] << " atoms\n";
  assert (etmp.atoms_in_growing_molecule()[c] >= 0);
#endif

  zresult += etmp.atoms_in_growing_molecule()[c];

  return 1;
}

int
IWReaction::_make_inter_particle_bond(Molecule & result,
                                      const Set_of_Atoms * scaffold_embedding,
                                      const Inter_Particle_Bond & b,
                                      const Enumeration_Temporaries & etmp) const
{
#ifdef DEBUG_MAKE_INTER_PARTICLE_BONDS
  cerr << "IWReaction::_make_inter_particle_bond: bond " << b.a1() << " to " << b.a2() << " embedding " << *scaffold_embedding << '\n';
#endif

  atom_number_t a1;
  if (! determine_atom_number(*scaffold_embedding, b.a1(), etmp, "IWReaction::_make_inter_particle_bond: a1", a1))
    return 0;

  atom_number_t a2 = INVALID_ATOM_NUMBER;
  if (! determine_atom_number(*scaffold_embedding, b.a2(), etmp, "IWReaction::_make_inter_particle_bond: a1", a2))
    return 0;

#ifdef DEBUG_MAKE_INTER_PARTICLE_BONDS
  cerr << "Making " << b << " atoms " << a1 << " and " << a2 << '\n';
#endif

  if (! b.OkSidechainIsotopeConstraint(result, a1, a2)) {
    return 0;
  }

  // Align before bonds are made.
  if (b.align_3d() > 0.0f) {
    if (! lillymol::Position3D(result, a1, b.align_3d(), a2)) {
      cerr << "IWReaction::_make_inter_particle_bond:cannot 3d align atoms " <<
              a1 << result.smarts_equivalent_for_atom(a1) << " and " <<
              a2 << result.smarts_equivalent_for_atom(a2) << '\n';
    }
  } else if (! b.OkDistance(result, a1, a2)) {
    return 0;
  }
  if (! result.add_bond(a1, a2, b.btype())) {
    cerr << "IWReaction::_make_inter_particle_bond:cannot make bond between " << a1 << " and " << a2 << '\n';
    return 0;
  }
  result.set_implicit_hydrogens_known(a1, 0);
  result.set_implicit_hydrogens_known(a2, 0);

  result.recompute_implicit_hydrogens(a1);
  result.recompute_implicit_hydrogens(a2);

  return 1;
}

int
IWReaction::_make_inter_particle_bonds(Molecule & result,
                                       const Set_of_Atoms * scaffold_embedding,
                                       const Enumeration_Temporaries & etmp) const
{
//cerr << "IWReaction::_make_inter_particle_bonds:checking " << ns << " sidechains\n";

  for (const Sidechain_Reaction_Site* s : _sidechains) {
    int nb = s->number_inter_particle_bonds();

    // cerr << "SC has " << nb << " inter particle bonds\n";

    for (int j = 0; j < nb; j++) {
      const Inter_Particle_Bond * b = s->inter_particle_bond(j);

      if (! _make_inter_particle_bond(result, scaffold_embedding, *b, etmp)) {
        return 0;
      }
    }
  }

  return 1;
}

/*
  Note that we just do stereo centres associated with sidechains
*/

/*int
IWReaction::_make_stereo_centres (Molecule & result,
                                  const Set_of_Atoms * scaffold_embedding,
                                  const Enumeration_Temporaries & etmp) const
{
  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    int ns = s->number_stereo_centres();

    for (int j = 0; j < ns; j++)
    {
      const Reaction_Stereo_Centre * c = s->stereo_centre(j);

      if (! c->process (result, scaffold_embedding, etmp))
        return 0;
    }
  }

  return 1;
}*/

//#define DEBUG_SET_BOND_LENGTH

int
IWReaction::_set_bond_length(Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Reaction_Bond_Length & b,
                             const Enumeration_Temporaries & etmp) const
{
#ifdef DEBUG_SET_BOND_LENGTH
  cerr << "set bond length atoms are " << b.a1() << " and " << b.a2() << '\n';
  b.write_msi(cerr, "  ", "bond_process");
#endif

  atom_number_t a1;
  if (! determine_atom_number(*scaffold_embedding, b.a1(), etmp, "IWReaction::_set_bond_length: a1", a1))
    return 0;

  atom_number_t a2;
  if (! determine_atom_number(*scaffold_embedding, b.a2(), etmp, "IWReaction::_set_bond_length: a1", a2))
    return 0;

#ifdef DEBUG_SET_BOND_LENGTH
  cerr << "IWReaction::_set_bond_length:setting bond between " << a1 << " and " << a2 << " to length " << b.desired_length() << '\n';
  cerr << result.smarts_equivalent_for_atom(a1) << " and " << result.smarts_equivalent_for_atom(a2) << '\n';
  cerr << "Before setting length " << result.distance_between_atoms(a1, a2) << '\n';
#endif

  result.set_bond_length(a1, a2, b.desired_length());

#ifdef DEBUG_SET_BOND_LENGTH
  cerr << "After setting length " << result.distance_between_atoms(a1, a2) << '\n';
#endif

  return 1;
}

int
IWReaction::_set_bond_lengths(Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp,
                              const Reaction_Site & r) const
{
  int nb = r.number_bond_length_specifications();

//cerr << "IWReaction::_set_bond_lengths:setting " << nb << " reaction bond length specifications\n";

  for (int j = 0; j < nb; j++)
  {
    const Reaction_Bond_Length * b = r.bond_length_specification(j);

    if (! _set_bond_length(result, scaffold_embedding, *b, etmp))
      return 0;
  }

  return 1;
}

int
IWReaction::_set_bond_lengths(Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp) const
{
//cerr << "Setting scaffold bond lengths\n";
  if (! _set_bond_lengths(result, scaffold_embedding, etmp, *this))
    return 0;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (! _set_bond_lengths(result, scaffold_embedding, etmp, *s))
      return 0;
  }

  return 1;
}

int
IWReaction::_set_dihedral_angles(Molecule & result,
                                 const Set_of_Atoms * scaffold_embedding,
                                 const Enumeration_Temporaries & etmp,
                                 const Reaction_Site & r) const
{
  int nd = r.number_dihedral_angle_specifications();

  for (int j = 0; j < nd; j++)
  {
    const Reaction_Dihedral_Angle * b = r.dihedral_angle_specification(j);

    if (! b->process(result, scaffold_embedding, etmp))
      return 0;
  }

  return 1;
}

int
IWReaction::_set_dihedral_angles(Molecule & result,
                                 const Set_of_Atoms * scaffold_embedding,
                                 const Enumeration_Temporaries & etmp) const
{
  if (! _set_dihedral_angles(result, scaffold_embedding, etmp, *this))
    return 0;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (! _set_dihedral_angles(result, scaffold_embedding, etmp, *s))
      return 0;
  }

  return 1;
}

int
IWReaction::_set_bond_angles(Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Enumeration_Temporaries & etmp,
                             const Reaction_Site & r) const
{
  int nd = r.number_bond_angle_specifications();

  for (int j = 0; j < nd; j++)
  {
    const Reaction_Bond_Angle * b = r.bond_angle_specification(j);

    if (! b->process(result, scaffold_embedding, etmp))
      return 0;
  }

  return 1;
}

int
IWReaction::_set_bond_angles(Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Enumeration_Temporaries & etmp) const
{
  if (! _set_bond_angles(result, scaffold_embedding, etmp, *this))
    return 0;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (! _set_bond_angles(result, scaffold_embedding, etmp, *s))
      return 0;
  }

  return 1;
}


int
IWReaction::_do_3d_replacements(Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp,
                              const Reaction_Site & r) const
{
  int nr = r.number_3d_replace_specifications();

  for (int j = 0; j < nr; j++)
  {
    const Reaction_3D_Replace * b = r.r3d_replace_specification(j);

    if (! b->process(result, scaffold_embedding, etmp))
      return 0;
  }

  return 1;
}

int
IWReaction::_do_3d_replacements(Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp) const
{
  if (! _do_3d_replacements(result, scaffold_embedding, etmp, *this))
    return 0;

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Sidechain_Reaction_Site * s = _sidechains[i];

    if (! _do_3d_replacements(result, scaffold_embedding, etmp, *s))
      return 0;
  }

  return 1;
}

/*
  Some _replace_atom items may be within the current fragment
*/

int
Scaffold_Reaction_Site::_do_intra_particle_replacements(Molecule & result,
                                          const Set_of_Atoms * scaffold_embedding) const
{
//cerr << "Scaffold_Reaction_Site:have " << _replace_atom.number_elements() << " replacements\n";

  for (const Replace_Atom* r : _replace_atom)
  {
    const Matched_Atom_in_Component & a1 = r->a1();
    const Matched_Atom_in_Component & a2 = r->a2();

//  cerr << "in scaffold? " << a1.in_scaffold() << ' ' << a2.in_scaffold() << '\n';
    if (a1.in_scaffold() || a2.in_scaffold())    // HUH, appears exactly wrong!
      continue;

    int m1 = a1.matched_atom();
    int m2 = a2.matched_atom();

    atom_number_t ma1;
    if (! valid_member_of_embedding(*scaffold_embedding, m1, "Scaffold_Reaction_Site::_do_intra_particle_replacements:", ma1))
      return 0;

    atom_number_t ma2;
    if (! valid_member_of_embedding(*scaffold_embedding, m2, "Scaffold_Reaction_Site::_do_intra_particle_replacements:", ma2))
      return 0;

    if (! result.stereo_preserving_substitute(ma1, ma2))
      return 0;
  }

  return 1;
}

int
IWReaction::_do_replacement(Molecule & result,
                            const Set_of_Atoms * scaffold_embedding,
                            const Replace_Atom & r,
                            const Enumeration_Temporaries & etmp) const
{
  if (!r.DoReplacement(result, scaffold_embedding, etmp)) {
    cerr << "IWReaction::_do_replacement:replacement failed " << r << '\n';
    return 0;
  }

  if (_3d)
  {
    assert (nullptr == "implement this some time");
  }

  return 1;
}

int
Replace_Atom::DoReplacement(Molecule& result,
                            const Set_of_Atoms* scaffold_embedding,
                            const Enumeration_Temporaries& etmp) const
{
  atom_number_t a1;
  if (! determine_atom_number(*scaffold_embedding, _a1, etmp, "Replace_Atom::DoReplacement:a1", a1))
    return 0;

  atom_number_t a2;
  if (! determine_atom_number(*scaffold_embedding, _a2, etmp, "Replace_Atom::DoReplacement:a2", a2))
    return 0;

//cerr << "Atoms are " << a1 << " and " << a2 << '\n';
  if (! result.stereo_preserving_substitute(a1, a2))
  {
    cerr << "Replace_Atom::DoReplacement:failed, scaffold embedding " << (*scaffold_embedding) << '\n';
    cerr << *this << '\n';
    return 0;
  }

  // no 3d.

  return 1;
}

int
Reaction_Site::DoReplacements(Molecule& result,
                             const Set_of_Atoms * scaffold_embedding,
                             Enumeration_Temporaries & etmp) const
{
  for (const auto* r : _replace_atom)
  {
//  cerr << *r << " doing atom replacements\n";
    if (! r->DoReplacement(result, scaffold_embedding, etmp))
      return 0;
  }

  return 1;
}

int
IWReaction::_do_replacements(Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Enumeration_Temporaries & etmp) const
{
  for (const Sidechain_Reaction_Site* s : _sidechains)
  {
    const int nr = s->number_atom_replacements();

    for (int j = 0; j < nr; j++)
    {
      const Replace_Atom * r = s->replace_atom(j);

      if (! _do_replacement(result, scaffold_embedding, *r, etmp))
        return 0;
    }
  }

  return 1;
}

/*
  A No_Reaction object looks like

    (0 No_Reaction
      (A C smarts "xxx");
      (A C smarts "yyy")
    )

  where the first query specifies the scaffold and the 2nd specifies the sidechain
*/

int
No_Reaction::construct_from_msi_object(const msi_object & msi)
{
  if (0 != msi.number_elements() && 2 != msi.number_elements())
  {
    cerr << "No_Reaction::construct_from_msi_object: must specify 0 or 2 query objects\n";
    return 0;
  }

  if (2 == msi.number_elements())      // queries written in expanded form in the object
  {
    if (! _scaffold_no_reaction.construct_from_msi_object(*(msi[0])))
    {
      cerr << "No_Reaction::construct_from_msi_object: cannot read first query\n";
      return 0;
    }
    if (! _sidechain_no_reaction.construct_from_msi_object(*(msi[1])))
    {
      cerr << "No_Reaction::construct_from_msi_object: cannot read 2nd query\n";
      return 0;
    }
  }

// Note that we don't do any checking against the case of _scaffold_no_reaction
// being specified twice, once as an expanded object and once as a smarts

  int na = msi.number_attributes();

  for (int i = 0; i < na; i++)
  {
    const msi_attribute * att = msi.attribute(i);

    if (NAME_OF_REACTION_COMMENT == att->name())
    {
      _comment = att->stringval();
      continue;
    }

    int j = 0;       // processing scaffold or sidechain?
    if (_scaffold_no_reaction.empty())
      j = 0;
    else if (_sidechain_no_reaction.empty())
      j = 1;

    int success = 1;

    if (NAME_OF_QUERY_SMARTS_ATTRIBUTE == att->name())
    {
      if (0 == j)
        success = _scaffold_no_reaction.create_from_smarts(att->stringval());
      else if (1 == j)
        success = _sidechain_no_reaction.create_from_smarts(att->stringval());
    }
    else if (NAME_OF_QUERY_FILE_ATTRIBUTE == att->name())
    {
      if (0 == j)
        success = _scaffold_no_reaction.read(att->stringval());
      else if (1 == j)
        success = _sidechain_no_reaction.read(att->stringval());
    }
    else
    {
      cerr << "No_Reaction::construct_from_msi_object: unrecognised attribute '" << att->name() << "'\n";
      return 0;
    }

    if (0 == success)
    {
      cerr << "No_Reaction::construct_from_msi_object: cannot parse '" << (*att) << "'\n";
      return 0;
    }
  }

  if (_scaffold_no_reaction.empty() || _sidechain_no_reaction.empty())
  {
    cerr << "No_Reaction::construct_from_msi_object: one or both queries not specified\n";
    return 0;
  }

  return 1;
}

int
No_Reaction::write_msi(std::ostream & os,
                        int & object_id,
                        int indentation)
{
  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << '(' << object_id++ << ' ' << NAME_OF_NO_REACTION_OBJECT << '\n';

  if (_comment.length())
    os << ind << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";

  _scaffold_no_reaction.write_msi(os, object_id, indentation + 2);
  _sidechain_no_reaction.write_msi(os, object_id, indentation + 2);

  os << ind << ")\n";

  return os.good();
}

Reaction_Stereo_Centre::Reaction_Stereo_Centre()
{
  _optional = 0;

  return;
}

int
Reaction_Stereo_Centre::debug_print (std::ostream & os) const
{
  os << "Reaction Stereo Centre\n";
  os <<  " centre ";
  _ssc[0].debug_print(os);
  os <<  " top front ";
  _ssc[1].debug_print(os);
  os <<  " top back ";
  _ssc[2].debug_print(os);
  os <<  " left down ";
  _ssc[3].debug_print(os);
  os <<  " right down ";
  _ssc[4].debug_print(os);

  return os.good();
}

void
Reaction_Stereo_Centre::all_atoms_in_scaffold()
{
  for (int i = 0; i < 5; i++)
  {
    _ssc[i].set_in_component(0);
  }
}

int
Reaction_Stereo_Centre::construct_from_msi_attribute (const msi_attribute * msi)
{
  const_IWSubstring m;
  msi->value(m);

  int nw = m.nwords();

  if (nw < 4 || nw > 5)
  {
    cerr << "Reaction_Stereo_Centre::construct_from_msi_attribute: attribute must have 4 or 5 tokens\n";
    cerr << m << '\n';
    return 0;
  }

  for (int i = 0; i < nw; i++)
  {
    const_IWSubstring token;
    (void) m.word(i, token);

    if (! _ssc[i].construct(token))
    {
      cerr << "Reaction_Stereo_Centre::construct_from_msi_attribute: cannot parse attribute " << i << " '" << token << "'\n";
      cerr << m << '\n';
      return 0;
    }

//  cerr << "Built stereo component '" << _ssc[i] << "' from '" << token << "'\n";
  }

  if (4 == nw)
    _ssc[4].set_whatever();

  return 1;
}

int
Reaction_Stereo_Centre::write_msi(std::ostream & os, const const_IWSubstring & ind,
                           const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A C " << attribute_name << " \"";
  for (int i = 0; i < 5; i++)
  {
    if (i > 0)
      os << ' ';
    os << _ssc[i];
  }

  os << "\")\n";

  return os.good();
}

Stereo_Centre_Component::Stereo_Centre_Component()
{
  _implicit_hydrogen = -1;

  return;
}

int
Stereo_Centre_Component::active() const
{
  if (_matched_atom < 0 && _implicit_hydrogen < 0)
    return 0;

  return 1;
}

int
Stereo_Centre_Component::ok() const
{
  if (_matched_atom < 0 && _implicit_hydrogen < 0 && _component < 0)
    return 1;

  if (_matched_atom >= 0 && _implicit_hydrogen > 0)
    return 0;

  if (_implicit_hydrogen > 0 && _component > 0)
    return 0;

  return 1;
}

int
Stereo_Centre_Component::debug_print (std::ostream & os) const
{
  if (_implicit_hydrogen)
    os << "Implicit Hydrogen";
  else
    os << "Atom " << _matched_atom << " in component " << _component;

  os << '\n';

  return os.good();
}

void
Stereo_Centre_Component::set_matched_atom (int m)
{
  assert (m >= 0);

  _matched_atom = m;
  _implicit_hydrogen = 0;

  return;
}

void
Stereo_Centre_Component::set_implicit_hydrogen (int h)
{
  if (h)
  {
    _implicit_hydrogen = 1;
    _matched_atom = -1;
    _component = -1;
  }
  else
  {
    _implicit_hydrogen = 0;
  }

  return;
}

int
Stereo_Centre_Component::construct (const const_IWSubstring & token)
{
  if ('H' == token)
  {
    _implicit_hydrogen = 1;
    _matched_atom = -1;
    _component = -1;
    return 1;
  }

  _implicit_hydrogen = 0;

  return Matched_Atom_in_Component::construct(token);
}

std::ostream &
operator << (std::ostream & os, const Stereo_Centre_Component & scc)
{
  if (! scc.active())
  {
    os << "Inactive";
    return os;
  }

  assert (scc.ok());

  if (scc.in_scaffold())
    os << "S.";
  else
    os << "R" << scc.in_component() << '.';

  if (scc.implicit_hydrogen() > 0)
    os << "H";

  if (scc.matched_atom() >= 0)
    os << scc.matched_atom();

  return os;
}
   
int
IWReaction::_make_stereo_centres (Molecule & result,
                                  const Set_of_Atoms * scaffold_embedding,
                                  const Enumeration_Temporaries & etmp) const
{
//cerr << "Processing " << _reaction_stereo_centre.size() << " stereo centres\n";

  int rc = 1;
  for (Reaction_Stereo_Centre* r : _reaction_stereo_centre) {
    if (! r->process(result, scaffold_embedding, etmp)) {
      rc = 0;
    }
  }

  return rc;
}
 
int
IWReaction::_make_cip_stereo_centres(Molecule & result,
                                  const Set_of_Atoms&  scaffold_embedding,
                                  const Enumeration_Temporaries & etmp) const {
//  cerr << "Processing " << _cip_stereo.size() << " CIP stereo centres, scaffold_embedding " << scaffold_embedding << '\n';

  int rc = 1;
  for (const ReactionCipStereo* cip : _cip_stereo) {
    if (! _make_cip_stereo_centre(result, scaffold_embedding, etmp, *cip)) {
      rc = 0;
    }
  }

  return rc;
}

int
IWReaction::_make_cip_stereo_centre(Molecule & result,
                                  const Set_of_Atoms&  scaffold_embedding,
                                  const Enumeration_Temporaries & etmp,
                                  const ReactionCipStereo& cip) const {
  atom_number_t a;
  if (! determine_atom_number(scaffold_embedding, cip.matched_atom(), etmp,
      "_make_cip_stereo_centre", a)) {
  }

  if (! result.ok_atom_number(a)) {
    cerr << "IWReaction::_make_cip_stereo_centre:invalid atom " << a << '\n';
    return 0;
  }

  const Chiral_Centre* existing = result.chiral_centre_at_atom(a);
  if (existing != nullptr) {
    if (result.CahnIngoldPrelogValue(existing) == cip.cip()) {
      return 1;
    }
    result.invert_chirality_on_atom(a);
    if (result.CahnIngoldPrelogValue(existing) == cip.cip()) {
      return 1;
    }
    cerr << "IWReaction::_make_cip_stereo_centre:not matched " << a << '\n';
    return 0;
  }

  return CreateChiralCentreMatchingCip(result, a, cip.cip());
}

/*
  If we detect an invalid valence involving atoms we have added, try to find a different Kekule form.
  To make things easy, we only check for 5 bonded carbons in 6 membered rings
*/

/*int
Reaction_Site::_do_find_kekule_forms_for_bad_valence (Molecule & result,
                                                      const Set_of_Atoms & embedding,
                                                      int offset,
                                                      int * initially_aromatic_atoms) const
{
  cerr << "Trying different Kekule forms for fixing bad valences '" << result.smiles() << "'\n";

  set_display_abnormal_valence_messages(0);

  int matoms = result.natoms();

  int * aromatic_atoms = new int[matoms];

  std::unique_ptr<int[]> freeme (aromatic_atoms);

  (void) result.ring_membership();   // force SSSR and bonds to know ring membership

  for (int i = 0; i < matoms; i++)
  {
//  if (! initially_aromatic_atoms[i])
//    continue;

    const Atom * ai = result.atomi (i);

    if (6 != ai->atomic_number())     // we key on 5 valent carbons
      continue;

    if (ai->nbonds() <= 4)    // looking good
      continue;

    if (! result.in_ring_of_given_size (i, 6))     // only 6 membered aromatic rings have Kekule forms
      continue;

    set_vector (aromatic_atoms, matoms, 0);
    identify_atoms_in_ring_system (result, i, aromatic_atoms, initially_aromatic_atoms);
    for (int j = 0; j < matoms; j++)
    {
      cerr << "Atom " << j << " in ring system " << aromatic_atoms[j] << '\n';
    }
    int tmp = result.find_kekule_form (aromatic_atoms);
    cerr << "Found kekule form ? " << tmp << '\n';
  }

  set_display_abnormal_valence_messages (1);

  return 1;
}*/

/*
  We want to figure out whether a fused ring should be added to a possibly aromatic
  system. This is a bit of a kludge. Our test is that the new ring must bring
  at least one new double bond into the system. This is probably OK, because
  alternating Kekule systems usually don't cross such rings
*/

static int
could_extend_aromatic_system (Molecule & m,
                              const Ring & r,
                              const int * process_these)
{
  int new_double_bonds = 0;    // to possibly be aromatic, must have at least 1 double bond not already fully counted

  for (Ring_Bond_Iterator i(r); i != r.zend(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    if (process_these[a1] && process_these[a2])
      continue;

    const Bond * b = m.bond_between_atoms (a1, a2);
    if (b->is_double_bond())
      new_double_bonds++;

    if (m.hcount (a1) > 1)    // could imagine cases where this is wrong, but hopefully too rare to worry about
      return 0;
  }

  return new_double_bonds;
}

static int
identify_fused_system (Molecule & m,
                       const Ring & r,
                       int * process_these)
{
  r.set_vector(process_these, 1);

  int rc = 1;

  int nfused = r.fused_ring_neighbours();
  for (int i = 0; i < nfused; i++)
  {
    const Ring * ri = r.fused_neighbour(i);

    if (ri->all_members_set_in_array(process_these, 1))
      continue;

    if (! could_extend_aromatic_system(m, *ri, process_these))
      continue;

    rc += identify_fused_system(m, *ri, process_these);
  }

  return rc;
}

static int
find_kekule_forms_for_bad_valence (Molecule & m,
                                   atom_number_t c,
                                   int * process_these)
{
  const Ring * r = m.ring_containing_atom(c);
  assert (nullptr != r);

  if (r->is_fused())
    identify_fused_system(m, *r, process_these);

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == process_these[i])
      continue;

    int icon = m.ncon(i);

    if (icon > 3)     // a 4 connected atom is never aromatic
      return 0;

    const Atom * a = m.atomi(i);

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = a->other(i, j);

      if (process_these[k])
        m.set_bond_type_between_atoms(i, k, SINGLE_BOND);
    }
  }

#ifdef ECHO_ATOMS_IN_RING_SYSTEM
  for (int i = 0; i < matoms; i++)
  {
    if (process_these[i])
      cerr << " atom " << i << " in ring system\n";
  }
#endif

  return m.find_kekule_form(process_these);
}

static int
find_kekule_forms_for_bad_valence (Molecule & m,
                                   atom_number_t c)
{
  int * tmp = new_int(m.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  return find_kekule_forms_for_bad_valence(m, c, tmp);
}

/*
  We have a molecule with an invalid valence. We only process the case of 5 valent carbon
  in a 6 membered ring
*/

int
IWReaction::_do_find_kekule_forms_for_bad_valence (Molecule & m) const
{
//cerr << "Bad valence in '" << m.smiles() << ' ' << m.name() << "'\n";

  int matoms = m.natoms();

  int rc = 1;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      continue;

    if (a->nbonds() <= 4)    // looks like a happy Carbon
      continue;

    if (! m.in_ring_of_given_size(i, 6))    // we only process this case
      continue;

    if (! find_kekule_forms_for_bad_valence(m, i))
      rc = 0;
  }

  return rc;
}

static int
convert_implicit_hydrogen_to_explicit_connection (Molecule & m,
                                                  Chiral_Centre & c)
{
  const Atom * a = m.atomi(c.a());

  assert (4 == a->ncon());

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = a->other(c.a(), i);

    if (c.involves(j))
      continue;

    c.implicit_hydrogen_is_now_atom_number(j);
    return 1;
  }
  
  cerr << "convert_implicit_hydrogen_to_explicit_connection:huh, didn't find a match!\n";
  return 0;
}

/*
  We are building a chiral centre.  The last connection was left
  unspecified.  Find a missing connection.  We set c->right_down
*/

int
Reaction_Stereo_Centre::_find_missing_connection (Molecule & m,
                                                  atom_number_t zatom,
                                                  Chiral_Centre * c) const
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();
//cerr << "Chiral centre has " << acon << " connections\n";

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);
//  cerr << "Atom " << zatom << " bonded to atom " << j << " involves? " << c->involves(j) << '\n';
    if (c->involves(j))
      continue;

    c->set_right_down(j);
    return 1;
  }

// All attached atoms are already in the chiral centre object

  if (3 == acon)    // let's assume it is an implicit hydrogen (could be lone pair?)
  {
    c->set_right_down(kChiralConnectionIsImplicitHydrogen);
    return 1;
  }

  cerr << "Reaction_Stereo_Centre::_find_missing_connection: No missing connection!\n";
  c->debug_print(cerr);
  return 0;
}

static void
report_chiral_centre_status (const Molecule & m,
                             atom_number_t a1,
                             atom_number_t a2,
                             const char * origin)
{
  if (a2 < 0)
    return;

  if (m.are_bonded(a1, a2))
    return;

  cerr << "Chiral centre component atoms " << a1 << " and " << a2 << " not bonded " << origin << '\n';

  return;
}

static void
print_diagnostic_info_for_invalid_chiral_centre(const Molecule & m,
                                const Chiral_Centre * c)
{
  c->debug_print(cerr);

  atom_number_t a = c->a();

  report_chiral_centre_status(m, a, c->top_front(), "top front");
  report_chiral_centre_status(m, a, c->top_back(), "top back");
  report_chiral_centre_status(m, a, c->left_down(), "left down");
  report_chiral_centre_status(m, a, c->right_down(), "right down");
}

static int
chiral_component_bonded (const Atom * a,
                         atom_number_t k)
{
  if (k < 0)
    return 1;

  return a->is_bonded_to(k);
}

static int
all_parts_of_chiral_centre_bonded_to_centre (const Molecule & m, 
                                const Chiral_Centre * c)
{
  const Atom * a = m.atomi(c->a());

  if (! chiral_component_bonded(a, c->top_front()))
    return 0;
  if (! chiral_component_bonded(a, c->top_back()))
    return 0;
  if (! chiral_component_bonded(a, c->left_down()))
    return 0;
  if (! chiral_component_bonded(a, c->right_down()))
    return 0;

  return 1;
}

//#define DEBUG_REACTION_STEREO_CENTRE_PROCESS

int
Reaction_Stereo_Centre::process (Molecule & result,
                                 const Set_of_Atoms * scaffold_embedding,
                                 const Enumeration_Temporaries & etmp) const
{
  atom_number_t a;
  if (! determine_atom_number(*scaffold_embedding, _ssc[0], etmp, "Reaction_Stereo_Centre::process: a", a))
    return 0;

#ifdef DEBUG_REACTION_STEREO_CENTRE_PROCESS
  debug_print(cerr);
  cerr << "Reaction_Stereo_Centre::process: embedding " << (*scaffold_embedding) << '\n';
  cerr << _ssc[0] << " becomes atom " << a << '\n';
  cerr << result.smiles() << '\n';
#endif

  if (result.ncon(a) < 3)
  {
    cerr << "Reaction_Stereo_Centre::make_stereo_centre: atom " << a << " has " << result.ncon(a) << " connections. Impossible\n";

#ifdef WRITE_NUMBERED_SMILES_ON_ERROR
    write_numbered_smiles(result, cerr);
#endif

    if (_optional)
    {
      cerr << "Optional centre, ignored\n";
      return 1;
    }

    return 0;
  }

  Chiral_Centre * c = new Chiral_Centre(a);

  int acon = result.ncon(c->a());

  if (_ssc[1].implicit_hydrogen())
    c->set_top_front(kChiralConnectionIsImplicitHydrogen);
  else if (! determine_atom_number(*scaffold_embedding, _ssc[1], etmp, "Reaction_Stereo_Centre::process: top front", a))
    return 0;
  else
    c->set_top_front(a);

  if (_ssc[2].implicit_hydrogen())
    c->set_top_back(kChiralConnectionIsImplicitHydrogen);
  else if (! determine_atom_number(*scaffold_embedding, _ssc[2], etmp, "Reaction_Stereo_Centre::process: top back", a))
    return 0;
  else
    c->set_top_back(a);

  if (_ssc[3].implicit_hydrogen())
    c->set_left_down(kChiralConnectionIsImplicitHydrogen);
  else if (! determine_atom_number(*scaffold_embedding, _ssc[3], etmp, "Reaction_Stereo_Centre::process: left down", a))
    return 0;
  else
    c->set_left_down(a);


  if (! all_parts_of_chiral_centre_bonded_to_centre(result, c))
  {
    cerr << "Reaction_Stereo_Centre::process:invalid chiral centre formed\n";
    print_diagnostic_info_for_invalid_chiral_centre(result, c);
    return 0;
  }

// For the last connection, we may not know whether there is an atom or an implicit hydrogen

  if (_ssc[4].is_whatever())
  {
    if (! _find_missing_connection(result, c->a(), c))
    {
      cerr << "Reaction_Stereo_Centre::process:no unspecified connection\n";
      return 0;
    }
  }
  else if (_ssc[4].implicit_hydrogen())
    c->set_right_down(kChiralConnectionIsImplicitHydrogen);
  else if (! determine_atom_number(*scaffold_embedding, _ssc[4], etmp, "Reaction_Stereo_Centre::process: right down", a))
    return 0;
  else
    c->set_right_down(a);

// Make sure we can add the chiral centre

  int ih = c->implicit_hydrogen_count();

  if (0 == ih && 4 == acon)    // good
    ;
  else if (1 == ih && 3 == acon)    // good
    ;
  else if (1 == ih && 4 == acon)   // convert the implicit Hydrogen to the connection
    convert_implicit_hydrogen_to_explicit_connection(result, *c);
  else
  {
    cerr << "Reaction_Stereo_Centre::process:atom " << c->a() << " connectivity mismatch acon " << acon << " ih " << ih << "\n";
    c->debug_print(cerr);
    cerr << result.smiles() << '\n';
    delete c;
    return 1;     // make this a harmless error
  }

  if (! result.valid_chiral_centre(c))
  {
    cerr << "Reaction_Site::make_stereo_centre:: constructed chiral centre invalid\n";

    c->debug_print(cerr);
    atom_number_t centre = c->a();

    delete c;

    cerr << "In the molecule, atom " << centre << " bonded to atoms ";

    for (int i = 0; i < result.ncon(centre); i++)
    {
      atom_number_t j = result.other(centre, i);
      cerr << ' ' << j;
    }
    cerr << '\n';

#ifdef WRITE_NUMBERED_SMILES_ON_ERROR
    result.remove_all_chiral_centres();
    write_numbered_smiles(result, cerr);
#endif

    return 0;
  }

  c->set_chirality_known(1);

// Remove any existing chiral centre

  if (nullptr != result.chiral_centre_at_atom(c->a()))
    result.remove_chiral_centre_at_atom(c->a());

  if (! result.add_chiral_centre(c))
  {
    cerr << "Reaction_Site::make_stereo_centre:: cannot add chiral centre to result molecule\n";
    c->debug_print(cerr);

    return 0;
  }

  return 1;
}

Reaction_Wedge_Bond::Reaction_Wedge_Bond() : Pair_of_Atoms(INVALID_ATOM_NUMBER, INVALID_ATOM_NUMBER)
{
  _direction = 0;

  return;
}

/*
  A wedge bond must be of the form

    (A I wedge_bond a1 a2 dir)

  where dir is +1 or -1
*/

int
Reaction_Wedge_Bond::construct_from_msi_attribute(const msi_attribute * att)
{
  if (3 != att->number_int_values())
  {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: attribute must have 3 integer values\n";
    return 0;
  }

  int a = att->int_multi_value(0);
  Pair_of_Atoms::set_a1(a);

  a = att->int_multi_value(1);
  Pair_of_Atoms::set_a2(a);

  if (Pair_of_Atoms::a1() == Pair_of_Atoms::a2())
  {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: atoms the same - impossible\n";
    return 0;
  }

  int d = att->int_multi_value(2);

  if (1 == d)
    _direction = 1;
  else if (-1 == d)
    _direction = -1;
  else if (0 == d)
    _direction = 0;
  else
  {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: invalid direction " << d << '\n';
    return 0;
  }

  return 1;
}

int
Reaction_Wedge_Bond::write_msi(std::ostream & os, const IWString & ind, const const_IWSubstring & attribute_name) const
{
  os << ind << "  (A I " << attribute_name << " (" << _a1 << ' ' << _a2 << ' ' << _direction << "))\n";

  return os.good();
}

int
Reaction_Wedge_Bond::process (Molecule & result,
                              const Set_of_Atoms & embedding,
                              int offset) const
{
  atom_number_t a1 = embedding[_a1] + offset;
  atom_number_t a2 = embedding[_a2] + offset;

  if (! result.are_bonded(a1, a2))
  {
    cerr << "Reaction_Wedge_Bond::do_place_wedge_bonds: atoms " << a1 << " and " << a2 << " are not bonded\n";
    return 0;
  }

//cerr << "Placing " << _direction << " wedge bond between atoms " << a1 << " and " << a2 << '\n';

  if (! result.set_wedge_bond_between_atoms(a1, a2, _direction))
  {
    cerr << "Reaction_Wedge_Bond::do_place_wedge_bonds: cannot place wedge bond between atoms " << a1 << " and " << a2 << '\n';
    return 0;
  }

  return 1;
}

int
IWReaction::determine_matched_atoms(Molecule & m,
                                    Substructure_Results & sresults)
{
  return _determine_matched_atoms_checking_inactives(m, sresults);
}

int
IWReaction::add_sidechain_reagent(int sidechain, Molecule& reagent,
                        const Sidechain_Match_Conditions& smc) {
  if (sidechain >= _sidechains.number_elements()) {
    cerr << "IWReaction::add_sidechain_reagment:invalid sidechain " << sidechain << '\n';
    return 0;
  }

  return _sidechains[sidechain]->add_reagent(reagent, smc);
}

std::optional<std::vector<Molecule>>
IWReaction::perform_reaction(Molecule& scaffold, Molecule& sidechain) {
  std::vector<Molecule> products;

  if (_sidechains.size() != 1) {
    return std::nullopt;
  }

  Substructure_Results scaffold_results;
  if (this->substructure_search(scaffold, scaffold_results) == 0) {
    return std::nullopt;
  }

  int n = scaffold_results.number_embeddings();
  products.resize(n);

  // Add the reagent to the sidechain, will be removed below.
  Sidechain_Match_Conditions smc;
  if (! _sidechains[0]->add_reagent(sidechain, smc)) {
    return std::nullopt;
  }

  int ndx = 0;

  for (const Set_of_Atoms* scaffold_embedding : scaffold_results.embeddings()) {
    if (! perform_reaction(&scaffold, scaffold_embedding, products[ndx])) {
      return std::nullopt;
    }
    ++ndx;
  }

  _sidechains[0]->remove_last_reagent();

  return products;
}

int
IWReaction::perform_reaction(Molecule & scaffold, resizable_array_p<Molecule> & products) {
  Substructure_Results sresults;
  if (! this->substructure_search(scaffold, sresults)) {
    return 0;
  }

  for (const Set_of_Atoms* embedding : sresults.embeddings()) {
    std::unique_ptr<Molecule> product = std::make_unique<Molecule>();
    if (! perform_reaction(&scaffold, embedding, *product)) {
      return products.size();
    }
    products << product.release();
  }

  return products.size();
}

/*int
Reaction_Site::_do_toggle_kekule_form (Molecule & m,
                                       const Substructure_Results & sresults)
{
  int nhits = sresults.number_embeddings();

//#define DEBUG_DO_TOGGLE_KEKULE_FORM
#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
  cerr << "Doing toggle Kekule form on " << nhits << " hits\n";
#endif

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
    cerr << "Embedding is " << (*e) << '\n';
    _toggle_kekule_form.debug_print (cerr);
    IWString initial_smiles = m.smiles();
#endif

    int changed;
    if (! _toggle_kekule_form.process (m, *e, changed))
    {
      cerr << "Reaction_Site::_do_toggle_kekule_form:toggle kekule form failed for '" << m.name() << "'\n";
    }

#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
    if (changed)
      cerr << "Changed from '" << initial_smiles << "' to '" << m.smiles() << "'\n";
#endif
  }

  return 1;
}*/

/*
  Before we broke bonds, one of the atoms from which a bond was being
  removed had a chiral centre. We saved that chiral centre in _saved_chiral_centres.

  We don't worry about lone pairs...
*/

int 
IWReaction::_do_restore_saved_chiral_centres(Molecule & result,
                                             const Set_of_Atoms & scaffold_embedding,
                                             Enumeration_Temporaries & etmp)
{
  Reaction_Site::do_restore_saved_chiral_centres(result, scaffold_embedding, etmp);

  const Molecule_and_Embedding ** reagent = etmp.reagent();

  int ns = _sidechains.number_elements();

  for (int i = 0; i < ns; i++)
  {
    const Set_of_Atoms * e = reagent[i]->embedding();

    _sidechains[i]->do_restore_saved_chiral_centres(result, *e, etmp);
  }

  return 1;
}

int
Reaction_Site::do_restore_saved_chiral_centres (Molecule & result,
                                                const Set_of_Atoms & embedding,
                                                Enumeration_Temporaries & etmp)
{
  int n = _saved_chiral_centre.number_elements();

//#define DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
  cerr << "Reaction_Site::do_restore_saved_chiral_centres:restoring " << n << " saved chiral centres\n";
#endif

  for (int i = n - 1; i >= 0; i--)
  {
    Chiral_Centre * ci = _saved_chiral_centre[i];

    if (_do_restore_saved_chiral_centre(result, embedding, etmp, ci))
    {
      result.add_chiral_centre(ci);
      _saved_chiral_centre.remove_no_delete(i);
#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
      cerr << "Centre " << i << " restored\n";
#endif
    }
  }

  return 1;
}

int
Reaction_Site::_do_restore_saved_chiral_centre (Molecule & result,
                                                const Set_of_Atoms & scaffold_embedding,
                                                Enumeration_Temporaries & etmp,
                                                Chiral_Centre * c)
{
  atom_number_t zatom = c->a();

#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
  cerr << "Trying to restore chiral centre ";
  c->debug_print(cerr);
#endif

/* If the reaction is placing an explicit chiral centre on ZATOM, don't do anything

  for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++)
  {
    const Reaction_Stereo_Centre * r = _reaction_stereo_centre[i];

    const Stereo_Centre_Component & scci = r->centre();

    atom_number_t a;

    if (! determine_atom_number (scaffold_embedding, scci, etmp, "IWReaction::_do_restore_saved_chiral_centre", a))
      continue;

    if (a == zatom)
      return 0;
  }

// No explicit chiral centre on the atom, proceed*/

  const Atom * a = result.atomi(zatom);

  if (4 == a->ncon())    // could still be chiral
    ;
  else if (3 == a->ncon() && 3 == a->nbonds() && 1 == result.implicit_hydrogens(zatom) && 1 == result.hcount(zatom))    // could still be chiral
    ;
  else               // no longer a chiral centre
    return 0;

// Identify those atoms that were detached during the reaction

  Set_of_Atoms atoms_no_longer_attached;

  int atoms_in_chiral_centre = 0;

  int i = 0;
  atom_number_t j;
  while (INVALID_ATOM_NUMBER != (j = c->next_atom(i)))
  {
    atoms_in_chiral_centre++;

    if (! a->is_bonded_to(j))     // ZATOM still bonded to atom J in RESULT, no change
      atoms_no_longer_attached.add(j);    // may need to change to a new atom
  }

// The same atoms still attached, and there were 4 of them - should not happen - we know a bond was broken

  if (atoms_no_longer_attached.empty() && 4 == atoms_in_chiral_centre)    // should not happen, all atoms the same, but we know a bond was broken? Maybe they broke and re-formed something???
  {
    cerr << "IWReaction::_possibly_add_back_chiral_centre:strange, nothing has changed at atom " << zatom << '\n';
    return 1;
  }

  Set_of_Atoms atoms_newly_appeared;

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
    cerr << "Atom " << j << " is attached\n";
#endif

    if (! c->involves(j))
      atoms_newly_appeared.add(j);
  }

#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
  cerr << "After reaction " << atoms_no_longer_attached.number_elements() << " atoms no longer attached, " << atoms_newly_appeared.number_elements() << " atoms newly attached\n";
#endif

// The case where an atom is replacing an implicit Hydrogen. All atoms the same, only 3 explicit atoms in the chiral centre, but now 4 connections

  if (atoms_no_longer_attached.empty() && 3 == atoms_in_chiral_centre && 4 == a->ncon())
  {
    atom_number_t t = atoms_newly_appeared[0];

    c->lone_pair_is_now_atom_number(t);

    return 1;
  }

// One atom has obviously replaced another

  if (1 == atoms_newly_appeared.number_elements() && 1 == atoms_no_longer_attached.number_elements())
  {
    atom_number_t afrom = atoms_no_longer_attached[0];
    atom_number_t ato = atoms_newly_appeared[0];

#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
    cerr << "Changing atom " << afrom << " to " << ato << '\n';
#endif

    c->change_atom_number(afrom, ato);

    return 1;
  }

// If one atom has disappeared, we can turn that position into an implicit hydrogen

  if (atoms_newly_appeared.empty() && 1 == atoms_no_longer_attached.number_elements() &&
      0 == c->implicit_hydrogen_count())
  {
    atom_number_t afrom = atoms_no_longer_attached[0];

    c->atom_is_now_implicit_hydrogen(afrom);

    return 1;
  }

// Unless we have any replacement specifications, we can't do anything

  if (_replace_atom.empty())
    return 0;

// If there is an implicit Hydrogen on the chiral centre, we can have one atom that isn't
// specified as a replacement

//atom_number_t possible_implicit_hydrogen = INVALID_ATOM_NUMBER;

// Identify the replacements

  int replacements_identified = 0;

  for (int i = 0; i < _replace_atom.number_elements(); i++) {
    const Replace_Atom * r = _replace_atom[i];

    atom_number_t a1;
    if (! determine_atom_number(scaffold_embedding, r->a1(), etmp, "IWReaction::_do_restore_saved_chiral_centres: a1", a1))
      continue;

    if (! atoms_no_longer_attached.contains(a1))   // strange
      continue;

    atom_number_t a2;
    if (! determine_atom_number(scaffold_embedding, r->a2(), etmp, "IWReaction::_do_restore_saved_chiral_centres: a2", a2))
      continue;

    if (! atoms_newly_appeared.contains(a1))   // should not happen
      continue;

    c->change_atom_number(a1, a2);

    replacements_identified++;

    atoms_no_longer_attached.remove_first(a1);
    atoms_newly_appeared.remove_first(a2);
  }

// If we found a one-to-one correspondence, we are done

  if (atoms_no_longer_attached.empty() && atoms_newly_appeared.empty()) {
    return 1;
  }

// If there is one unmatched newly appeared atom, and the chiral centre had a lone pair, make that match

  if (atoms_no_longer_attached.empty() && 1 == atoms_newly_appeared.number_elements()
      && 1 == c->implicit_hydrogen_count()) {
    atom_number_t a = atoms_newly_appeared[0];

    c->implicit_hydrogen_is_now_atom_number(a);

    return 1;
  }

// Give up

#ifdef DEBUG_DO_RESTORE_SAVED_CHIRAL_CENTRES
  cerr << "Cannot restore chiral centre at atom " << a << '\n';
#endif

  return 0;
}

SiteCipStereo::SiteCipStereo() {
  _atom = -1;
  _rs = CahnIngoldPrelog::kUnspecified;
}

ReactionCipStereo::ReactionCipStereo() {
  _rs = CahnIngoldPrelog::kUnspecified;
}
