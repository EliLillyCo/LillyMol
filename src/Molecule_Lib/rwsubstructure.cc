#include <iostream>
#include <fstream>
#include <memory>
#include <limits>
using std::cerr;
using std::endl;

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

/*
  Functions associated with reading and writing substructure query objects and
  such.
*/

#include "msi_object.h"
#include "misc.h"

#include "molecule_to_query.h"
#include "istream_and_type.h"
#include "substructure.h"
#include "target.h"
#include "misc2.h"
#include "parse_smarts_tmp.h"
#include "mdl_molecule.h"
#include "path.h"
#include "smiles.h"

#ifdef FINGERPRINT_SUBSTRUCTURE_SEARCHES
#include "iwmfingerprint.h"
#endif

#define NAME_OF_RING_SPECIFIER_OBJECT "ring_specifier"
#define NAME_OF_RING_SYSTEM_SPECIFIER_OBJECT "ring_system_specifier"
#define NAME_OF_QUERY_BOND_OBJECT "query_bond"
#define NAME_OF_QUERY_ATOM_SPECIFIER_OBJECT "query_atom_specifier"
#define NAME_OF_ENVIRONMENT_OBJECT "environment"
#define NAME_OF_QUERY_ATOM_PREFERENCE_OBJECT "query_atom_preference"
#define NAME_OF_ENVIROMENT_REJECTION_OBJECT "environment_no_match"
#define NAME_OF_QUERY_ATOM_OBJECT "query_atom"
#define NAME_OF_ELEMENT_HITS_NEEDED_OBJECT "element_hits_needed"
#define NAME_OF_ELEMENTS_NEEDED_OBJECT "elements_needed"

#define NAME_OF_AND_OPERATOR "and"
#define NAME_OF_OR_OPERATOR "or"
#define NAME_OF_XOR_OPERATOR "xor"
#define NAME_OF_LOW_PRIORITY_AND_OPERATOR "low_priority_and"

#define NAME_OF_BOND_ATTRIBUTE "bond"
#define NAME_OF_SINGLE_BOND_ATTRIBUTE "single_bond"
#define NAME_OF_DOUBLE_BOND_ATTRIBUTE "double_bond"
#define NAME_OF_TRIPLE_BOND_ATTRIBUTE "triple_bond"
#define NAME_OF_AROMATIC_BOND_ATTRIBUTE "aromatic_bond"
#define NAME_OF_BOND_SMARTS_ATTRIBUTE "bond_smarts"

#define NAME_OF_AROMATICITY_ATTRIBUTE "aromatic"
#define NAME_OF_OR_AROMATICITY_SPECIFIER "or_aromatic"
#define NAME_OF_AND_AROMATICITY_SPECIFIER "and_aromatic"
#define NAME_OF_REJECTION_ATTRIBUTE "rejection"
#define NAME_OF_ATOMIC_NUMBER_ATTRIBUTE "atomic_number"
#define NAME_OF_ATOMIC_SYMBOL_ATTRIBUTE "atomic_symbol"
#define NAME_OF_NCON_ATTRIBUTE "ncon"
#define NAME_OF_FORMAL_CHARGE_ATTRIBUTE "formal_charge"
#define NAME_OF_CHARGED_ATTRIBUTE "charged"
#define NAME_OF_HCOUNT_ATTRIBUTE "hcount"
#define NAME_OF_NBONDS_ATTRIBUTE "nbonds"
#define NAME_OF_NRINGS_ATTRIBUTE "nrings"
#define NAME_OF_RING_BOND_COUNT_ATTRIBUTE "ring_bond_count"
#define NAME_OF_CHIRALITY_SPECIFIER "chirality"
#define NAME_OF_UNSATURATION_ATTRIBUTE "unsaturation"
#define NAME_OF_DAYLIGHT_X_ATTRIBUTE "daylight_x"
#define NAME_OF_LONE_PAIR_SPECIFIER "lone_pair"
#define NAME_OF_ISOTOPE_ATTRIBUTE "isotope"
#define NAME_OF_USER_ATOM_TYPE_ATTRIBUTE "user_atom_type"
#define NAME_OF_NUMBER_ISOTOPIC_ATOMS_ATTRIBUTE "isotopic_atoms"
#define NAME_OF_CHIRAL_CENTRE_ATTRIBUTE "chiral_center"
#define NAME_OF_SORT_MATCHES_BY "sort_matches"

// The isolated_ring attribute is deprecated, use fused system size

#define NAME_OF_ISOLATED_RING_ATTRIBUTE "isolated_ring"

#define NAME_OF_INITIAL_ATOM_NUMBER_ATTRIBUTE "initial_atom_number"
#define NAME_OF_ATOM_MAP_NUMBER_ATTRIBUTE "atom_map_number"
#define NAME_OF_FUSED_SYSTEM_SIZE_ATTRIBUTE "fused_system_size"
#define NAME_OF_RING_ID_SPECIFIER "ring_id"
#define NAME_OF_FRAGMENT_ID_SPECIFIER "fragment_id"
#define NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER "attached_heteroatom_count"
#define NAME_OF_COMMENT_ATTRIBUTE "comment"
#define NAME_OF_TEXT_IDENTIFIER_ATTRIBUTE "text_identifier"
#define NAME_OF_VERSION_ATTRIBUTE "version"
#define NAME_OF_VINYL_ATTRIBUTE "vinyl"
#define NAME_OF_ARYL_ATTRIBUTE  "aryl"
#define NAME_OF_DEFINE_HETEROATOMS_ATTRIBUTE "define_heteroatoms"
#define NAME_OF_HETEROATOMS_ATTRIBUTE "heteroatom_count"
#define NAME_OF_ENVIRONMENT_REJECTION "environment_rejection"
#define NAME_OF_ENVIRONMENT_HITS_NEEDED_ATTRIBUTE "environment_hits_needed"
#define NAME_OF_ENVIRONMENT_OPTIONAL_REJECTION "optional_rejection"
#define NAME_OF_MAX_ENVIRONMENT_MATCHES_PER_ATTACHMENT_POINT "max_env_matches_per_anchor"
#define NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS "env_matches_share_attachment_points"
#define NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE "hits_needed"
#define NAME_OF_FUSED_RING_SYSTEM_SIZE "fused_ring_system_size"
#define NAME_OF_ONE_EMBEDDING_PER_START_ATOM "one_embedding_per_start_atom"
#define NAME_OF_UNIQUE_EMBEDDINGS_ONLY "unique_embeddings_only"
#define NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE "embeddings_do_not_overlap"
#define NAME_OF_NORMALISE_RC_PER_HITS_NEEDED "normalise_rc_per_hits_needed"
#define NAME_OF_SUBTRACT_FROM_RC "subtract_from_rc"
#define NAME_OF_RESPECT_INITIAL_NUMBERING_ATTRIBUTE "respect_initial_numbering"
#define NAME_OF_MAX_MATCHES_TO_FIND "max_matches_to_find"
#define NAME_OF_NCON_IGNORE_SINGLY_CONNECTED "ncon_ignore_singly_connected"
#define NAME_OF_NUMERIC_VALUE_ATTRIBUTE "numeric_value"
#define NAME_OF_INCLUDE_IN_EMBEDDING_ATTRIBUTE "include_in_embedding"
#define NAME_OF_DO_NOT_PERCEIVE_SYMMETRIC_EQUIVALENTS "do_not_perceive_symmetric_equivalents"
#define NAME_OF_ATOM_SMARTS_ATTRIBUTE "atom_smarts"
#define NAME_OF_SMARTS_ATTRIBUTE "smarts"
#define NAME_OF_SMILES_ATTRIBUTE "smiles"
#define NAME_OF_RING_SIZE_ATTRIBUTE "ring_size"
#define NAME_OF_RINGS_THAT_MUST_MATCH_RING_SIZE_ATTRIBUTE "rings_that_must_match_ring_size"
#define NAME_OF_AROMATIC_RING_SIZE_ATTRIBUTE "aromatic_ring_size"
#define NAME_OF_ALIPHATIC_RING_SIZE_ATTRIBUTE "aliphatic_ring_size"
#define NAME_OF_RING_FUSED_ATTRIBUTE "fused"
#define NAME_OF_NATOMS_ATTRIBUTE "natoms"
#define NAME_OF_PREFERENCE_VALUE_ATTRIBUTE "preference_value"
#define NAME_OF_SUM_ALL_PREFERENCE_HITS_ATTRIBUTE "sum_all_preference_hits"
#define NAME_OF_DISTANCE_BETWEEN_HITS_ATTRIBUTE "distance_between_hits"
#define NAME_OF_DISTANCE_BETWEEN_HITS_NCHECK_ATTRIBUTE "distance_between_hits_ncheck"
#define NAME_OF_FAIL_IF_EMBEDDINGS_TOO_CLOSE_ATTRIBUTE "fail_if_embeddings_too_close"
#define NAME_OF_SORT_BY_PREFERENCE_VALUE_ATTRIBUTE "sort_by_preference_value"
#define NAME_OF_HETEROATOMS_MATCHED_ATTRIBUTE "heteroatoms_matched"
#define NAME_OF_RING_ATOMS_MATCHED_ATTRIBUTE "ring_atoms_matched"
#define NAME_OF_HETEROATOMS_IN_RING_ATTRIBUTE "heteroatoms_in_ring"
#define NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE "all_hits_in_same_fragment"
#define NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE "only_match_largest_fragment"
#define NAME_OF_NUMBER_FRAGMENTS_ATTRIBUTE "number_fragments"
#define NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE "match_each_component"
#define NAME_OF_NUMBER_SPINACH_GROUPS "number_spinach_groups"
#define NAME_OF_NUMBER_NON_SPINACH_GROUPS "number_non_spinach_groups"
#define NAME_OF_ATOMS_IN_SPINACH_ATTRIBUTE "atoms_in_spinach_group"
#define NAME_OF_LENGTH_OF_SPINACH_ATTRIBUTE "length_of_spinach_group"
#define NAME_OF_DISTANCE_TO_ANOTHER_RING_ATTRIBUTE "distance_to_another_ring"
#define NAME_OF_NET_FORMAL_CHARGE_ATTRIBUTE "net_formal_charge"

#define NAME_OF_ALL_RINGS_KEKULE "all_rings_kekule"

#define NAME_OF_SYMMETRY_DEGREE_ATTRIBUTE "symmetry_degree"
#define NAME_OF_SYMMETRY_GROUP_ATTRIBUTE "symmetry_group"

#define NAME_OF_ENVIRONMENT_MUST_MATCH_UNMATCHED_ATOMS_ATTRIBUTE "environment_must_match_unmatched_atoms"
#define NAME_OF_ENVIRONMENT_CAN_MATCH_IN_RING_ATOMS "environment_can_match_in_ring_atoms"

#define NAME_OF_SPINACH_ATOMS_ATTRIBUTE "spinach_atoms"
#define NAME_OF_INTER_RING_ATOMS_ATTRIBUTE "inter_ring_atoms"

#define NAME_OF_NO_OTHER_SUBSTITUENTS_ALLOWED_ATTRIBUTE "no_other_substituents_allowed"

#define NAME_OF_FUSED_RING_COUNT_SPECIFIER "fused_ring_count"
#define NAME_OF_ISOLATED_RING_COUNT_SPECIFIER "isolated_ring_count"
#define NAME_OF_AROMATIC_RING_COUNT_SPECIFIER "aromatic_ring_count"
#define NAME_OF_NON_AROMATIC_RING_COUNT_SPECIFIER "non_aromatic_ring_count"
#define NAME_OF_STRONGLY_FUSED_RING_COUNT_SPECIFIER "strongly_fused_ring_count"

#define NAME_OF_FUSED_AROMATIC_NEIGHBOURS_SPECIFIER "fused_aromatic_neighbours"
#define NAME_OF_FUSED_NON_AROMATIC_NEIGHBOURS_SPECIFIER "fused_non_aromatic_neighbours"
#define NAME_OF_LARGEST_NUMBER_SHARED_BONDS_SPECIFIER "shared_bonds"
#define NAME_OF_STRONGLY_FUSED_RING_NEIGHBOUR_SPECIFIER "strongly_fused_neighbours"
#define NAME_OF_DISTANCE_BETWEEN_ROOT_ATOMS_SPECIFIER "dist_between_root_atoms"
#define NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER "no_matched_atoms_between"
#define NAME_OF_LINK_ATOMS_SPECIFIER "link_atoms"

#define NAME_OF_ISOLATED_RING_OBJECTS "isolated_ring_objects"

#define NAME_OF_SAVE_HITS_ATTRIBUTE "save_hits"

#define NAME_OF_IMPLICIT_RING_CONDITION "implicit_ring_condition"

#define NAME_OF_ATOMS_WITH_PI_ELECTRONS "atoms_with_pi_electrons"

#define NAME_OF_UNMATCHED_ATOMS "unmatched_atoms"
#define NAME_OF_MIN_FRACTION_ATOMS_MATCHED "min_fraction_atoms_matched"
#define NAME_OF_MAX_FRACTION_ATOMS_MATCHED "max_fraction_atoms_matched"

#define NAME_OF_HYDROGEN_OK_AS_ENVIRONMENT_MATCH "hydrogen_ok"

static int set_element_hits_needed_during_molecule_to_query = 1;

/*
  In several places we need to quickly ascertain whether or not an msi_object
  specifies a rejection or not
*/

static int
is_rejection (const msi_object & msi)
{
  const msi_attribute * rj = msi.attribute(NAME_OF_REJECTION_ATTRIBUTE);

  if (NULL == rj)    // no rejection attribute present
   return 0;

  int rc;
  if (! rj->value(rc))
  {
    cerr << "The " << NAME_OF_REJECTION_ATTRIBUTE << " attribute must have an int value\n";
    cerr << msi << endl;
    iwabort();
    return 0;
  }

  return rc;
}

static int
fetch_aromaticity (const msi_object & msi, aromaticity_type_t & aromaticity)
{
  int rc = 0;

  const msi_attribute * arom;
  int i = 0;
  while (NULL != (arom = msi.attribute(NAME_OF_AROMATICITY_ATTRIBUTE, i++)))
  {
    int tmp;
    if (! arom->value(tmp))
    {
      cerr << "fetch_aromaticity: bad value '" << arom->stringval() << "'\n";
      assert (NULL == "This is not good");
    }

    if (0 == tmp)
    {
      if (SUBSTRUCTURE_NOT_SPECIFIED == aromaticity)
        aromaticity = NOT_AROMATIC;
      else
        SET_ALIPHATIC_ATOM(aromaticity);

      rc++;
    }
    else if (tmp > 0)
    {
      if (SUBSTRUCTURE_NOT_SPECIFIED == aromaticity)
        aromaticity = AROMATIC;
      else
        SET_AROMATIC_ATOM(aromaticity);

      rc++;
    }
  }

  return rc;
}

/*
  Atomic numbers are handled specially - they are not a min_max_specifier
  Note that regardless of how many attributes we process, we just increment
  ATTRIBUTES_SPECIFIED by 1
*/

static int
fetch_elements (const msi_object & msi,
                resizable_array<const Element *> & ele,
                resizable_array<int> & element_unique_id,
                int & attributes_specified)
{
  int ii = 0;

  const msi_attribute * attribute;
  while (NULL != (attribute = msi.attribute(NAME_OF_ATOMIC_NUMBER_ATTRIBUTE, ii++)))
  {
    int n = attribute->number_int_values();
    if (0 == n)
    {
      cerr << "Attribute is not of type int '" << (*attribute) << "'\n";
      return 0;
    }

    for (int i = 0; i < n; i++)
    {
      atomic_number_t z = attribute->int_multi_value(i);
      if (! REASONABLE_ATOMIC_NUMBER(z))
      {
        cerr << "fetch_elements::invalid atomic number " << z << endl;
        return 0;
      }

      const Element * e = get_element_from_atomic_number(z);
      if (NULL == e)
      {
        cerr << "fetch_elements:no element for atomic number " << z << endl;
        return 0;
      }
      
      ele.add(e);
      element_unique_id.add(e->unique_id());
    }
  }

  ii = 0;
  while (NULL != (attribute = msi.attribute(NAME_OF_ATOMIC_SYMBOL_ATTRIBUTE, ii++)))
  {
    int n = attribute->number_string_values();
    if (0 == n)
    {
      cerr << "Attribute has no string values '" << (*attribute) << "'\n";
      return 0;
    }

    for (int i = 0; i < n; i++)
    {
      const IWString * s = attribute->string_multi_value(i);
//    cerr << "ATOMIC SYMBOL " << *s << endl;

      const Element * e = get_element_from_symbol_no_case_conversion(*s);
      if (NULL == e && auto_create_new_elements())
        e = create_element_with_symbol(*s);

      if (NULL == e)
      {
        cerr << "fetch atomic numbers, no element for '" << (*attribute) << "'\n";
        return 0;
      }
      
      ele.add(e);
      element_unique_id.add(e->unique_id());
//    cerr << "Hash value " << e->atomic_symbol_hash_value() << endl;
    }
  }

  assert (ele.ok());

  if (ele.number_elements())
    attributes_specified++;

  return 1;
}

static int
get_float_attribute (const msi_object & msi,
                     const char * att_name,
                     float & v,
                     float minval, float maxval)
{
  const msi_attribute * att = msi.attribute(att_name);

  if (NULL == att)
    return 1;

  if (! att->value(v) || v < minval || v > maxval)
  {
    cerr << "invalid " << att_name << " specification '" << att->stringval() << "'\n";
    return 0;
  }


  return 1;
}

#define MIN_NOT_SPECIFIED 9393
#define MAX_NOT_SPECIFIED 14163

/*
  Process all the ATTRIBUTE_NAME attributes in an msi_object.
  Append resulting int values to SPECIFIER.
*/

static int
append_int_values (const msi_object & msi,
                   const IWString & attribute_name,
                   int min_val, int max_val,
                   Min_Max_Specifier<int> & specifier,
                   int & attributes_specified)
{
  assert (specifier.ok());
  assert (attributes_specified >=0);

  int ii = 0;

  const msi_attribute * attribute;
  while (NULL != (attribute = msi.attribute(attribute_name, ii++)))
  {
    if (0 == attribute->number_int_values())
    {
      cerr << "Attribute is not of type int '" << (*attribute) << "'\n";
      return 0;
    }

    int n = attribute->number_int_values();
    for (int i = 0; i < n; i++)
    {
      int j = attribute->int_multi_value(i);
      if (MIN_NOT_SPECIFIED != min_val && j < min_val)
      {
        cerr << "Value out of range, min is " << min_val << " '" << (*attribute) << "'\n";
        return 0;
      }
      
      if (MIN_NOT_SPECIFIED != max_val && j > max_val)
      {
        cerr << "Value out of range, min is " << min_val << " '" << (*attribute) << "'\n";
        return 0;
      }

      specifier.add(j);
      attributes_specified++;
    }
  }

  assert (specifier.ok());

  return 1;
}

static int
set_int_value (const msi_attribute & msi, int & value, int min_val, int max_val)
{
  int tmp;
  if (! msi.value(tmp))
  {
    cerr << "msi attribute value not interpretable as int\n";
    cerr << msi << endl;
    return 0;
  }

  if (MIN_NOT_SPECIFIED != min_val && tmp < min_val)
  {
    cerr << "msi attribute value out of range, got " << tmp << " min is " << min_val << endl;
    cerr << msi << endl;
    return 0;
  }

  if (MAX_NOT_SPECIFIED != max_val && tmp > max_val)
  {
    cerr << "msi attribute value out of range, got " << tmp << " max is " << max_val << endl;
    cerr << msi << endl;
    return 0;
  }

  value = tmp;
  return 1;
}

/*
  Unlike all the others, atomic numbers are an iwarchive, rather than an
  iwminmaxspc
*/

/*static int
really_gruesome (resizable_array<int> & specifier, const msi_object * msi,
                 const IWString & attribute_name, int & attributes_specified,
                 int min_value_allowed, int max_value_allowed)
{
  int ii = 0;

  const msi_attribute * attribute;
  while (NULL != (attribute = msi->attribute (attribute_name, ii++)))
  {
    int tmp;

    if (! set_int_value(*attribute, tmp, min_value_allowed, max_value_allowed))
      return 0;

    specifier.add(tmp);
    attributes_specified++;
  }

  return 1;
}*/

/*
  This extracts specifications for an attribute from an msi object.
  Note that even if we process multiple attributes, we only increment
  ATTRIBUTES_SPECIFIED by one, since it is really Substructure_Atom_Specifier::_attributes_specified
*/

template <typename T>
int
_really_gruesome (Min_Max_Specifier<T> & specifier, const msi_object & msi,
                  const IWString & attribute_name, int & attributes_specified,
                  T min_value_allowed, T max_value_allowed)
{
  int attributes_specified_here = 0;

  if (! append_int_values(msi, attribute_name, min_value_allowed, max_value_allowed,
                           specifier, attributes_specified_here))
    return 0;

// If we have some values specified, no need to go looking for specifications
// of min or max, that would be illegal (so will never happen right?)

  if (attributes_specified_here)
  {
    attributes_specified++;
    return 1;
  }

  IWString name = "min_" + attribute_name;
  const msi_attribute * attribute;
  if (NULL != (attribute = msi.attribute(name)))
  {
    int tmp;
    if (! set_int_value(*attribute, tmp, min_value_allowed, max_value_allowed))
      return 0;
    specifier.set_min(tmp);
    attributes_specified_here++;
  }

  name = "max_" + attribute_name;
  if (NULL != (attribute = msi.attribute(name)))
  {
    int tmp;
    if (! set_int_value(*attribute, tmp, min_value_allowed, max_value_allowed))
      return 0;
    specifier.set_max(tmp);
    attributes_specified_here++;
  }

  if (attributes_specified_here)
    attributes_specified++;
  else
    specifier.set_match_any(1);

  return 1;
}

template int _really_gruesome (Min_Max_Specifier<int> & specifier, const msi_object & msi,
                  const IWString & attribute_name, int & attributes_specified,
                  int min_value_allowed, int max_value_allowed);

template <typename T>
int
really_gruesome (Min_Max_Specifier<T> & specifier,
                 const msi_object & msi,
                 const IWString & attribute_name,
                 int & attributes_specified,
                 T min_value_allowed, T max_value_allowed)
{
  assert (specifier.ok());

  int rc = _really_gruesome(specifier, msi, attribute_name, attributes_specified,
                             min_value_allowed, max_value_allowed);

  if (rc)
    return rc;

  cerr << "Cannot parse '" << attribute_name << "' in msi object\n";

  cerr << msi << endl;

  return 0;
}

template int really_gruesome (Min_Max_Specifier<int> & specifier,
                 const msi_object & msi,
                 const IWString & attribute_name,
                 int & attributes_specified,
                 int min_value_allowed, int max_value_allowed);

/*static int
fetch_attribute (const msi_object & msi,
                 const char * attribute_name,
                 int & destination,
                 int minval,
                 int maxval,
                 int & attributes_specified,
                 int ok_not_present = 1)
{
  const msi_attribute * att = msi.attribute(attribute_name);
  if (NULL == att)
    return ok_not_present;       // attribute not present

  int tmp;
  if (! att->value(tmp))
  {
    cerr << "The " << attribute_name << " attribute must be a whole number\n";
    cerr << msi;
    return 0;
  }

  if (tmp < minval || tmp > maxval)
  {
    cerr << "The " << attribute_name << " attribute must be between " << minval << " and " << maxval << endl;
    cerr << msi;
    return 0;
  }

  destination = tmp;

  attributes_specified++;

  return 1;
}*/

int
Substructure_Atom_Specifier::construct_from_msi_object(const msi_object & msi)
{
  assert (ok());
  _attributes_specified = 0;

  const msi_attribute * att = msi.attribute(NAME_OF_PREFERENCE_VALUE_ATTRIBUTE);
  if (att)
  {
    if (! att->value(_preference_value))
    {
      cerr << "Bad preference value\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_CHARGED_ATTRIBUTE);
  if (att)
  {
    if(NULL != msi.attribute(NAME_OF_FORMAL_CHARGE_ATTRIBUTE))
      cerr << "Substructure_Atom_Specifier::construct_from_msi_object: warning, the " << NAME_OF_CHARGED_ATTRIBUTE << " and " << NAME_OF_FORMAL_CHARGE_ATTRIBUTE << " attributes may clash\n";

    _formal_charge.resize(8);
    for (int i= 1; i < 4; i++)
    {
      _formal_charge.add(i);
      _formal_charge.add(-i);
    }
  }

  if (! fetch_elements(msi, _element, _element_unique_id, _attributes_specified))
    return 0;

  if (! really_gruesome(_ncon, msi, NAME_OF_NCON_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_nbonds, msi, NAME_OF_NBONDS_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_nrings, msi, NAME_OF_NRINGS_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_ring_bond_count, msi, NAME_OF_RING_BOND_COUNT_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_ncon2, msi, "ncon2", _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_formal_charge, msi, NAME_OF_FORMAL_CHARGE_ATTRIBUTE, _attributes_specified, -4, 4))
    return 0;

  if (! really_gruesome(_hcount, msi, NAME_OF_HCOUNT_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_ring_size, msi, NAME_OF_RING_SIZE_ATTRIBUTE, _attributes_specified, 3, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_attached_heteroatom_count, msi, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, 
                              _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_lone_pair_count, msi, NAME_OF_LONE_PAIR_SPECIFIER, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_unsaturation, msi, NAME_OF_UNSATURATION_ATTRIBUTE, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_daylight_x, msi, NAME_OF_DAYLIGHT_X_ATTRIBUTE, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_isotope, msi, NAME_OF_ISOTOPE_ATTRIBUTE, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_aryl, msi, NAME_OF_ARYL_ATTRIBUTE, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_vinyl, msi, NAME_OF_VINYL_ATTRIBUTE, 
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_fused_system_size, msi, NAME_OF_FUSED_SYSTEM_SIZE_ATTRIBUTE,
                         _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_symmetry_degree, msi, NAME_OF_SYMMETRY_DEGREE_ATTRIBUTE, 
                        _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! _fetch_symmetry_group(msi))
    return 0;

  (void) fetch_aromaticity(msi, _aromaticity);

  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    _attributes_specified++;

  att = msi.attribute(NAME_OF_CHIRALITY_SPECIFIER);
  if (att)
  {
    if (! att->value(_chirality))
    {
      cerr << "Substructure_Atom_Specifier::construct_from_msi_object: bad chirality specifier\n";
      cerr << (*att) << endl;
      return 0;
    }

    _attributes_specified++;

//  Do something to audit the value
  }

  att = msi.attribute(NAME_OF_USER_ATOM_TYPE_ATTRIBUTE);
  if (att)
  {
    if (! att->value(_userAtomType))
    {
      cerr << "Substructure_Atom_Specifier::construct_from_msi_object: bad User Atom type\n";
      cerr << (*att) << endl;
      return 0;
    }

    _attributes_specified++;

//  Do something to audit the value
  }
// for backward compatibility, parse the deprecated "isolated_ring" attribute

  att = msi.attribute(NAME_OF_ISOLATED_RING_ATTRIBUTE);
  if (att)
  {
    int iso;
    if (! att->value(iso) || iso < 0)
    {
      cerr << "Substructure_Atom_Specifier::construct_from_msi_object: bad isolated ring specifier\n";
      cerr << (*att) << endl;
      return 0;
    } 

    if (0 == iso)
      _fused_system_size.set_min(2);
    else
      _fused_system_size.add(1);


    _attributes_specified++;
  }

  if (! really_gruesome(_heteroatoms_in_ring, msi, NAME_OF_HETEROATOMS_IN_RING_ATTRIBUTE, _attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_aromatic_ring_sizes, msi, NAME_OF_AROMATIC_RING_SIZE_ATTRIBUTE, _attributes_specified, 3, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_aliphatic_ring_sizes, msi, NAME_OF_ALIPHATIC_RING_SIZE_ATTRIBUTE, _attributes_specified, 3, MAX_NOT_SPECIFIED))
    return 0;

  att = msi.attribute(NAME_OF_ALL_RINGS_KEKULE);
  if (NULL != att)
  {
    att->value(_all_rings_kekule);    // should check the value...
    _attributes_specified++;
  }

//if (0 == _attributes_specified)
//  cerr << "Warning, no attributes specified\n";

  assert (ok());

  return 1;
}

int
Substructure_Atom_Specifier::_fetch_symmetry_group(const msi_object & msi)
{
  const msi_attribute * att = msi.attribute(NAME_OF_SYMMETRY_GROUP_ATTRIBUTE);
  if (NULL == att)
    return 1;

  if (1 != att->number_int_values())
  {
    cerr << "Substructure_Atom_Specifier::_fetch_symmetry_group:multiple values not supported, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 0;
  }

  if (! att->value(_symmetry_group) || _symmetry_group <= 0)
  {
    cerr << "Substructure_Atom_Specifier::_fetch_symmetry_group:invalid symmetry group specification\n";
    return 0;
  }

//cerr << "Substructure_Atom_Specifier::_fetch_symmetry_group:set " << _symmetry_group << endl;

  return 1;
}

/*
  Create a query bond.

  ATOMS must be an array of pointers to all query atoms found so far.
*/

/*int
Substructure_Bond::construct_from_msi_object (const msi_object * msi)
{
  assert (ok());

  int attributes_specified = 0;

  if (0 == msi->attribute_count ("type"))
  {
    attributes_specified = 1;
  }
  else if (! iwarchive<int>::add_values_from_msi_object (msi, "type"))
  {
    cerr << "Substructure_Bond::construct_from_msi_object: archive add failed\n";
    return 0;
  }
  else if (! _check_bond_types())
    return 0;
  else
    attributes_specified += _number_elements;

// Fetch the values for _nrings.

  if (! really_gruesome (_nrings, msi, NAME_OF_NRINGS_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  const msi_attribute * att = msi->attribute (NAME_OF_AROMATICITY_ATTRIBUTE);
  if (att)
  {
    int ar;
    if (! att->value (ar) || ar < 0)
    {
      cerr << "The " << NAME_OF_AROMATICITY_ATTRIBUTE << " specifier must be a non negative whole number\n";
      return 0;
    }

    set_and_aromatic (ar);
    attributes_specified++;
  }

  att = msi->attribute (NAME_OF_OR_AROMATICITY_SPECIFIER);
  if (att)
  {
    int ar;
    if (! att->value (ar) || ar < 0)
    {
      cerr << "The " << NAME_OF_OR_AROMATICITY_SPECIFIER << " specifier must be a non negative whole number\n";
      return 0;
    }

    set_or_aromatic (ar);
    attributes_specified++;
  }

  if (0 == attributes_specified)
  {
    cerr << "Substructure_Bond::construct_from_msi_object: no attributes specified\n";
    cerr << (*msi);
    cerr << "Warning only\n";
  }

  return 1;
}*/

/*int
Substructure_Bond::construct_from_msi_object (const msi_object * msi,
                                              extending_resizable_array<Substructure_Atom *> & completed)
{
  assert (completed);

  if (! Substructure_Bond_Specifier::construct_from_msi_object (msi))
    return 0;

  const msi_attribute * att = msi->attribute ("other");
  if (NULL == att)
  {
    cerr << "msi object for bond has no \"other\" attribute\n";
    return 0;
  }

  int i;
  if (! att->value (i) || i < 0 || NULL == completed[i] ||
      ! completed[i]->ok())
  {
    cerr << "Atom specifier at other end invalid '" << (*att) << "'\n";
    return 0;
  }

  _a1 = completed[i];

  return 1;
}*/

/*
  The multi-bond-type attribute looks like

    (A I bond (a t1 t2 t3...))

  where A is the atom at the other end of the bond and T1, T2, are the possible types
*/

int
Substructure_Bond::write_as_msi_attribute (std::ostream & os,
                                           int indentation) const
{
  assert (ok());

  IWString ind;
  ind.extend(indentation, ' ');

  os << ind;

  if (NULL != _b)
    return _write_as_smarts(os);

  os << "(A I ";

  int types_set = 0;

  if (_bond_types & SINGLE_BOND)
    types_set++;
  if (_bond_types & DOUBLE_BOND)
    types_set++;
  if (_bond_types & TRIPLE_BOND)
    types_set++;
  if (_bond_types & AROMATIC_BOND)
    types_set++;

  if (1 == types_set)
  {
    if (_bond_types & SINGLE_BOND)
      os << NAME_OF_SINGLE_BOND_ATTRIBUTE;
    else if (_bond_types & DOUBLE_BOND)
      os << NAME_OF_DOUBLE_BOND_ATTRIBUTE;
    else if (_bond_types & TRIPLE_BOND)
      os << NAME_OF_TRIPLE_BOND_ATTRIBUTE;
    else if (_bond_types & AROMATIC_BOND)
      os << NAME_OF_AROMATIC_BOND_ATTRIBUTE;

    os << ' ' << _a1->unique_id() << ")\n";

    return os.good();
  }

// The more complex case of multiple types

  os << NAME_OF_BOND_ATTRIBUTE << " (" << _a1->unique_id();

  if (SINGLE_BOND & _bond_types)
    os << ' ' << 1;
  if (DOUBLE_BOND & _bond_types)
    os << ' ' << 2;
  if (TRIPLE_BOND & _bond_types)
    os << ' ' << 3;
  if (AROMATIC_BOND & _bond_types)
    os << ' ' << 4;

  os << "))\n";

  return os.good();
}

int
Substructure_Bond::_write_as_smarts (std::ostream & os) const
{
  os << "(A C " << NAME_OF_BOND_SMARTS_ATTRIBUTE << " \"" << _a1->unique_id() << ' ';

  (void) __write_as_smarts(os);

  os << "\")\n";

  return os.good();
}

/*
  We have a bond type which may contain OR'd bond types
*/

static int
write_ord_bond_types (bond_type_t bt,
                      std::ostream & os)
{
  int written = 0;

  if (IS_SINGLE_BOND(bt))
  {
    os << '-';
    written++;
  }

  if (IS_DOUBLE_BOND(bt))
  {
    if (written)
      os << ',';
    os << '=';
    written++;
  }

  if (IS_TRIPLE_BOND(bt))
  {
    if (written)
      os << ',';
    os << '#';
    written++;
  }

  if (IS_AROMATIC_BOND(bt))
  {
    if (written)
      os << ',';
    os << ':';
    written++;
  }

  return os.good();
}

int
Substructure_Bond::__write_as_smarts (std::ostream & os) const
{
  if (0 != _bond_types)
    write_ord_bond_types(_bond_types, os);

  if (NULL == _b)
    return os.good();

  if (0 != _bond_types)
    os << ';';

  _b->smarts(os);

//cerr << "Substructure_Bond::__write_as_smarts: '" << _b->smarts (cerr) << "' operators = " << _logexp.number_operators() << endl;

  const Substructure_Bond_Specifier_Base * b = _b;
  for (int i = 0; i < _logexp.number_operators(); i++)
  {
    const_IWSubstring op;
    _logexp.op(i, op);

    if ('|' == op)     // we don't use | for OR here
      os << ',';
    else
      os << op;

    b = b->next();
    b->smarts(os);
  }

  return os.good();
}

/*
  The case when the Substructure_Bond object cannot be written as an attribute
*/

/*int
Substructure_Bond::_write_msi (std::ostream & os,
                               const IWString & ind) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    os << ind << "  (A I type  " << _things[i] << ")\n";
  }

  int j;
  if (_nrings.min (j))
    os << ind << "  (A I min_nrings " << j << ")\n";
  if (_nrings.max (j))
    os << ind << "  (A I max_nrings " << j << ")\n";

  _nrings.write_msi (os, NAME_OF_NRINGS_ATTRIBUTE, ind.nchars() + 2);

  if (SUBSTRUCTURE_NOT_SPECIFIED != _and_aromatic)
    os << ind << "  (A I " << NAME_OF_AND_AROMATICITY_SPECIFIER << ' ' << _and_aromatic << ")\n";

  if (_or_aromatic)
    os << ind << "  (A I " << NAME_OF_OR_AROMATICITY_SPECIFIER << ' ' << _or_aromatic << ")\n";

  return os.good();
}*/

/*int
Substructure_Bond::write_msi (std::ostream & os,
                              int & object_id, 
                              int indentation)
{
  IWString ind;
  if (indentation)
    ind.extend (indentation, ' ');

  os << ind << '(' << object_id++ << ' ' << NAME_OF_QUERY_BOND_OBJECT << endl;
  (void) _write_msi (os, ind);
  os << ind << ")\n";

  return os.good();
}*/

/*
*/

/*int
Substructure_Bond::write_as_msi_object (std::ostream & os,
                                        int & object_id, 
                                        int indentation) const
{
  IWString ind;
  if (indentation)
    ind.extend (indentation, ' ');

  os << ind << "(" << object_id++ << " Query_Bond\n";

  int other = _a1->unique_id();
  os << ind << "  (A I other " << other << ")\n";

  (void) _write_msi (os, ind);

  os << ind << ")\n";

  return os.good();
}*/

/*int
Substructure_Bond::write_as_msi_object_or_attribute (std::ostream & os,
                                                     int & object_id,
                                                     int indentation) const
{
  if (NULL == b)
    return write_as_msi_attribute (os, indentation);
  else
    return write_as_msi_object (os, object_id, indentation);
}*/

static int
all_atomic_numbers_positive (const resizable_array<const Element *> & element)
{
  for (int i = 0; i < element.number_elements(); ++i)
  {
    if (element[i]->atomic_number() <= 0)
      return 0;
  }

  return 1;
}

int
Substructure_Atom_Specifier::write_msi(std::ostream & os, int indentation,
                                        int unary_operator, int op)
{
  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  if (0 == unary_operator)
    os << ind << "(A I " << NAME_OF_REJECTION_ATTRIBUTE << " 1)\n";

  if (IW_LOGEXP_UNDEFINED == op)     // don't write the operator
    ;
  else if (IW_LOGEXP_AND == op)
    os << ind << "(A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"" << NAME_OF_AND_OPERATOR << "\")\n";
  else if (IW_LOGEXP_OR == op)
    os << ind << "(A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"" << NAME_OF_OR_OPERATOR << "\")\n";
  else if (IW_LOGEXP_XOR == op)
    os << ind << "(A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"" << NAME_OF_XOR_OPERATOR << "\")\n";
  else if (IW_LOGEXP_LOW_PRIORITY_AND == op)
    os << ind << "(A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"" << NAME_OF_LOW_PRIORITY_AND_OPERATOR << "\")\n";
  else
  {
    abort();
  }

  if (_preference_value)
    os << ind << "(A I " << NAME_OF_PREFERENCE_VALUE_ATTRIBUTE << ' ' << _preference_value << ")\n";

//for (int i = 0; i < _atomic_number.number_elements(); i++)
//  os << ind << "(A I atomic_number " << _atomic_number[i] << ")\n";

  int na = _element.number_elements();
  if (na > 0)
  {
    if (all_atomic_numbers_positive(_element))
    {
      os << ind << "(A I " << NAME_OF_ATOMIC_NUMBER_ATTRIBUTE << ' ';
      if (1 == na)
        os << _element[0]->atomic_number();
      else
      {
        os << '(';
        for (int i = 0; i < na; i++)
        {
          if (i)
            os << ' ';
          os << _element[i]->atomic_number();
        }
        os << ')';
      }
    }
    else
    {
      os << ind << "(A C " << NAME_OF_ATOMIC_SYMBOL_ATTRIBUTE << ' ';
      if (1 == na)
        os << '"' << _element[0]->symbol() << '"';
      else
      {
        os << '"';
        for (int i = 0; i < na; i++)
        {
          if (i)
            os << ' ';
          os << _element[i]->symbol();
        }
        os << '"';
      }
    }

    os << ")\n";
  }

  _ncon.write_msi(os, NAME_OF_NCON_ATTRIBUTE, indentation);
  _nbonds.write_msi(os, NAME_OF_NBONDS_ATTRIBUTE, indentation);
  _formal_charge.write_msi(os, NAME_OF_FORMAL_CHARGE_ATTRIBUTE, indentation);
  _hcount.write_msi(os, NAME_OF_HCOUNT_ATTRIBUTE, indentation);
  _nrings.write_msi(os, NAME_OF_NRINGS_ATTRIBUTE, indentation);
  _ring_bond_count.write_msi(os, NAME_OF_RING_BOND_COUNT_ATTRIBUTE, indentation);
  _ncon2.write_msi(os, "ncon2", indentation);
  _ring_size.write_msi(os, NAME_OF_RING_SIZE_ATTRIBUTE, indentation);
  _aromatic_ring_sizes.write_msi(os, NAME_OF_AROMATIC_RING_SIZE_ATTRIBUTE, indentation);
  _aliphatic_ring_sizes.write_msi(os, NAME_OF_ALIPHATIC_RING_SIZE_ATTRIBUTE, indentation);
  _attached_heteroatom_count.write_msi(os, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, indentation);
  _lone_pair_count.write_msi(os, NAME_OF_LONE_PAIR_SPECIFIER, indentation);
  _unsaturation.write_msi(os, NAME_OF_UNSATURATION_ATTRIBUTE, indentation);
  _daylight_x.write_msi(os, NAME_OF_DAYLIGHT_X_ATTRIBUTE, indentation);
  _isotope.write_msi(os, NAME_OF_ISOTOPE_ATTRIBUTE, indentation);
  if (_userAtomType != 0)
   	os << ind << "(A I " << NAME_OF_USER_ATOM_TYPE_ATTRIBUTE << " " << _userAtomType << ")\n";
  _aryl.write_msi(os, NAME_OF_ARYL_ATTRIBUTE, indentation);
  _vinyl.write_msi(os, NAME_OF_VINYL_ATTRIBUTE, indentation);
  _fused_system_size.write_msi(os, NAME_OF_FUSED_SYSTEM_SIZE_ATTRIBUTE, indentation);
  _heteroatoms_in_ring.write_msi(os, NAME_OF_HETEROATOMS_IN_RING_ATTRIBUTE, indentation);
  _symmetry_degree.write_msi(os, NAME_OF_SYMMETRY_DEGREE_ATTRIBUTE, indentation);
  
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
  {
    if (IS_AROMATIC_ATOM(_aromaticity))
      os << ind << "(A I " << NAME_OF_AROMATICITY_ATTRIBUTE << " 1)\n";
    if (IS_ALIPHATIC_ATOM(_aromaticity))
      os << ind << "(A I " << NAME_OF_AROMATICITY_ATTRIBUTE << " 0)\n";
  }

  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
    os << ind << "(A I " << NAME_OF_CHIRALITY_SPECIFIER << " " << _chirality << ")\n";

  if (_symmetry_group > 0)
    os << ind << "(A I " << NAME_OF_SYMMETRY_GROUP_ATTRIBUTE << ' ' << _symmetry_group << ")\n";

  if (SUBSTRUCTURE_NOT_SPECIFIED != _all_rings_kekule)
   os << ind << "(A I " << NAME_OF_ALL_RINGS_KEKULE << ' ' << _all_rings_kekule << ")\n";

  return os.good();
}

/*
  After a substructure_atom is written, it must write its children. 
*/

int
Substructure_Atom::_write_children_msi (std::ostream & os, 
                               int & object_id,
                               int indentation)
{
  int nc = _children.number_elements();
  for (int i = 0; i < nc; i++)
  {
    _children[i]->write_msi(os, object_id, NAME_OF_QUERY_ATOM_OBJECT, indentation);
  }

  return os.good();
}

/*
  The Substructure_Atom's get their _unique_id as their object id. We therefore
  need another series of numbers for the Substructure_Atom_Specifier's and 
  Substructure_Bond
*/

int
Substructure_Atom::write_msi (std::ostream & os, int & object_id,
                              const const_IWSubstring & stringid, int indentation)
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << _unique_id << ' ' << stringid << endl;

  if (_initial_atom_number >= 0)
    os << ind << "  (A I " << NAME_OF_INITIAL_ATOM_NUMBER_ATTRIBUTE << " " << _initial_atom_number << ")\n";

  if (_atom_map_number >= 0)
    os << ind << "  (A I " << NAME_OF_ATOM_MAP_NUMBER_ATTRIBUTE << " " << _atom_map_number << ")\n";

  if (_or_id)
    os << ind << "  (A I or " << _or_id << ")\n";

  if (0 == _match_as_match_or_rejection)
    os << ind << "  (A I " << NAME_OF_REJECTION_ATTRIBUTE << " 1)\n";

  if (! Substructure_Atom_Specifier::write_msi(os, indentation + 2))
    return 0;

  if (_ring_id)
    os << ind << "  (A I " << NAME_OF_RING_ID_SPECIFIER << ' ' << _ring_id << ")\n";

  if (_fragment_id)
    os << ind << "  (A I " << NAME_OF_FRAGMENT_ID_SPECIFIER << ' ' << _fragment_id << ")\n";

  if (_text_identifier.length())
    os << ind << "  (A C " << NAME_OF_TEXT_IDENTIFIER_ATTRIBUTE << " '" << _text_identifier << "')\n";

  if (0 == _include_in_embedding)
    os << ind << "  (A I " << NAME_OF_INCLUDE_IN_EMBEDDING_ATTRIBUTE << " 0)\n";

  if (_numeric_value.is_set())
  {
    double tmp;
    _numeric_value.value(tmp);
    os << ind << "  (A D " << NAME_OF_NUMERIC_VALUE_ATTRIBUTE << ' ' << tmp << ")\n";
  }

//cerr << "Writing " << _components.number_elements() << " components, operator ";
//_operator.debug_print (cerr);

  int nco = _components.number_elements();
  for (int i = 0; i < nco; i++)
  {
    os << ind << "  (" << object_id++ << " Query_Atom_Specifier\n";
    Substructure_Atom_Specifier * a = _components[i];
    if (i > 0)   // must include the operator
      a->write_msi(os, indentation + 4, _operator.unary_operator(i), _operator.op(i - 1));
    else                  // the last one does not have an operator
      a->write_msi(os, indentation + 4, _operator.unary_operator(i));

    os << ind << "  )\n";
  }

  if (_environment.active())
    _environment.write_msi(os, object_id, indentation + 2);

  for (int i = 0; i < _bonds.number_elements(); i++)
  {
    _bonds[i]->write_as_msi_attribute(os, indentation + 2);
  }

  if (_sum_all_preference_hits)
    os << ind << "  (A I " << NAME_OF_SUM_ALL_PREFERENCE_HITS_ATTRIBUTE << " 1)\n";

  int np = _preferences.number_elements();
  for (int i = 0; i < np; i++)
  {
    os << ind << "  (" << (object_id++) << ' ' << NAME_OF_QUERY_ATOM_PREFERENCE_OBJECT << endl;
    _preferences[i]->write_msi(os, indentation + 4);
    os << ind << "  )\n";
  }

  os << ind << ")\n";

  if (! _write_children_msi(os, object_id, indentation))
    return 0;

  return os.good();
}

int
Substructure_Environment::write_msi (std::ostream & os,
                                     int & object_id,
                                     int indentation)
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << ' ';
  if (_query_environment_match_as_rejection)
    os << NAME_OF_ENVIROMENT_REJECTION_OBJECT << endl;
  else
    os << NAME_OF_ENVIRONMENT_OBJECT << endl;

  if (_or_id)
    os << ind << "  (A I or " << _or_id << ")\n";

  if (_and_id)
    os << ind << "  (A I and " << _and_id << ")\n";

  if (_hits_needed.is_set())
    _hits_needed.write_msi(os, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, indentation + 2);

  if (_no_other_substituents_allowed)
    os << ind << "  (A I " << NAME_OF_NO_OTHER_SUBSTITUENTS_ALLOWED_ATTRIBUTE << ' ' << _no_other_substituents_allowed << ")\n";

  if (_max_environment_matches_per_attachment_point < std::numeric_limits<int>::max())
    os << ind << "  (A I " << NAME_OF_MAX_ENVIRONMENT_MATCHES_PER_ATTACHMENT_POINT << ' ' << _max_environment_matches_per_attachment_point << ")\n";

  if (std::numeric_limits<int>::max() != _max_matches_to_find)
    os << ind << "  (A I " << NAME_OF_MAX_MATCHES_TO_FIND << ' ' << _max_matches_to_find << ")\n";

  os << ind << "  (A I " << NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS << ' ' << _environments_can_share_attachment_points << ")\n";

  os << ind << "  (A I " << NAME_OF_HYDROGEN_OK_AS_ENVIRONMENT_MATCH << ' ' << _hydrogen_ok_as_environment_match << ")\n";

  IWString tmp;
  if (! _bond.bond_type_as_string(tmp))
  {
    cerr << "Substructure_Environment::write_msi: what kind of bond?\n";
    _bond.debug_print(cerr, ind);
  }

  os << ind << "  (A I " << tmp << "_bond";
  int np = _possible_parents.number_elements();
  assert (np);     // must have at least one parent
  if (1 == np)
    os << ' ' << _possible_parents[0]->unique_id();
  else
  {
    os << " (";
    for (int i = 0; i < np; i++)
    {
      Substructure_Atom * a = _possible_parents[i];
      if (i > 0)
        os << ' ';
      os << a->unique_id();
    }
    os << ')';
  }
  os << ')' << endl;

  for (int i = 0; i < _number_elements; i++)
  {
    Substructure_Atom * a = _things[i];
    if (! a->write_msi(os, object_id, NAME_OF_QUERY_ATOM_OBJECT, indentation + 2))
      break;
  }

  os << ind << ")\n";

  return os.good();
}

/*
  In order to deal with the NOT operator, we write each component as a separate
  environment object
*/

int
Substructure_Atom_Environment::write_msi (std::ostream & os, int & object_id, int indentation) const
{
  IWString ind;
  ind.extend(indentation, ' ');

  for (int i = 0; i < _number_elements; i++)
  {
    os << ind << '(' << object_id++ << ' ' << NAME_OF_ENVIRONMENT_OBJECT << endl;

    if (0 == _operator.unary_operator(i))
      os << ind << "  (A I " << NAME_OF_REJECTION_ATTRIBUTE << " 1)\n";
    if (_operator.number_operators() && i < _number_elements - 1)
    {
      int op = _operator.op (i);

      os << ind << "  (A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"";
      if (IW_LOGEXP_AND == op)
        os << NAME_OF_AND_OPERATOR;
      else if (IW_LOGEXP_OR == op)
        os << NAME_OF_OR_OPERATOR;
      else if (IW_LOGEXP_XOR == op)
        os << NAME_OF_XOR_OPERATOR;
      else if (IW_LOGEXP_LOW_PRIORITY_AND == op)
        os << NAME_OF_LOW_PRIORITY_AND_OPERATOR;

      os << "\")\n";
    }

    _things[i]->write_msi(os, object_id, NAME_OF_QUERY_ATOM_OBJECT, indentation + 2);

    os << ind << ")\n";
  }


  return os.good();
}

int
Substructure_Atom::_create_preference_from_msi_object (const msi_object & msi)
{
  Substructure_Atom_Specifier * a = new Substructure_Atom_Specifier;

  if (! a->construct_from_msi_object(msi))
  {
    delete a;
    return 0;
  }

  _preferences.add(a);

  return 1;
}

/*
  There is a common function of examining an object and extracting its
  OPERATOR attribute and then adding the appropriate operator to a
  logical_expression
*/

static int
extract_operator(const msi_object & msi,
                 IW_Logical_Expression & logexp,
                 int default_operator,
                 const char * caller)
{
  IWString op;
  if (! msi.string_value_for_attribute(NAME_OF_OPERATOR_ATTRIBUTE, op))
    logexp.add_operator(default_operator);
  else if (NAME_OF_AND_OPERATOR == op)
    logexp.add_operator(IW_LOGEXP_AND);
  else if (NAME_OF_OR_OPERATOR == op)
    logexp.add_operator(IW_LOGEXP_OR);
  else if (NAME_OF_XOR_OPERATOR == op)
    logexp.add_operator(IW_LOGEXP_XOR);
  else if (NAME_OF_LOW_PRIORITY_AND_OPERATOR == op)
    logexp.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
  else
  {
    cerr << caller << ": unrecognised operator '" << op << "'\n";
    return 0;
  }

  cerr << "Parsed '" << op << "' now ";
  logexp.debug_print(cerr);

  if (is_rejection(msi))
    logexp.set_unary_operator(logexp.number_operators(), 0);

  return 1;
}

int
Substructure_Atom::_add_component (const msi_object & msi)
{
  assert (ok());
  if (msi.number_elements())
  {
    cerr << "Substructure_Atom::_add_component: components for specifiers cannot have sub-objects\n";
    cerr << msi;
    return 0;
  }

  Substructure_Atom_Specifier * tmp = new Substructure_Atom_Specifier;
  if (! tmp->construct_from_msi_object(msi))
  {
    cerr << "Substructure_Atom::_add_component: component atom not formed\n";
    cerr << msi;
    delete tmp;
    return 0;
  }

  if (0 == _components.number_elements())    // no operator with first object
    ;
  else if (! extract_operator(msi, _operator, IW_LOGEXP_OR, "Substructure_Atom::_add_component"))
    return 0;

  if (is_rejection(msi))
    _operator.set_unary_operator(_operator.number_results() - 1, 0);

  _components.add(tmp);

  assert (ok());
  return 1;
}

/*
  Add a bond to a query atom
  When adding the first bond, that will be the bond to the parent.
  In that case, we must notify the parent that an extra child hs
  been created.
*/

int
Substructure_Atom::_add_bond (Substructure_Bond * b)
{
  _bonds.add(b);

  if (1 == _bonds.number_elements())     // first bond, must notify parent of extra child
  {
    _bond_to_parent = b;
    _parent = b->a();
    _parent->notify_extra_child(this);
  }
  else
  {
    if (b->a()->unique_id() > _parent->unique_id())
    {
      cerr << "Substructure_Atom::_add_bond: Warning, possible bond mis-ordering\n";
      cerr << "Atom " << _unique_id << ": parent is " << _parent->unique_id() <<
              " but ring closure bond to " << b->a()->unique_id() << endl;
      cerr << "Be prepared for possible trouble\n";
    }
  }

  return 1;
}

/*int
Substructure_Atom::_add_bond (const msi_object & msi,
                              extending_resizable_array<Substructure_Atom *> & completed)
{
  Substructure_Bond * b = new Substructure_Bond;
  if (! b->construct_from_msi_object (msi, completed))
  {
    delete b;
    return 0;
  }

  if (this == b->a())
  {
    cerr << "Substructure_Atom::_add_bond: atom bonded to itself\n";
    cerr << (*msi);
    return 0;
  }

  return _add_bond (b);
}*/

/*
  Look like

    1 -
    1 2 3 -!@

    A set of atoms and then the bond
*/

static int
parse_smarts_bond_attribute (const msi_attribute * att,
                             Set_of_Atoms & other_end,
                             IWString & bond_smarts)
{
  const IWString & Satt = att->stringval();    // 

  const int nw = Satt.nwords();

  if (nw < 2)
  {
    cerr << "parse_smarts_bond_attribute: bond smarts must have at least two tokens '" << Satt << "'\n";
    return 0;
  }

  const int natoms = nw-1;

  other_end.resize_keep_storage(natoms);

  int i = 0;
  const_IWSubstring Sanum;
  for (int j = 0; j < natoms; ++j)
  {
    Satt.nextword(Sanum, i);

    atom_number_t zatom;
    if (! Sanum.numeric_value(zatom) || zatom < 0)
    {
      cerr << "parse_smarts_bond_attribute: invalid atom number '" << Sanum << "'\n";
      return 0;
    }

    other_end.add_if_not_already_present(zatom);
  }

  Satt.nextword(bond_smarts, i);

  return bond_smarts.length();
}

int
Substructure_Atom::_process_attribute_smarts_bond (const msi_attribute * att,
                                                   extending_resizable_array<Substructure_Atom *> & completed)
{
  Set_of_Atoms other_end;
  IWString bond_smarts;

  if (! parse_smarts_bond_attribute(att, other_end, bond_smarts))
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond:invalid bond smarts " << *att << endl;
    return 0;
  }

  if (other_end.number_elements() > 1)
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond:cannot handle multi-valued attachment points\n";
    return 0;
  }

  const atom_number_t o = other_end[0];

  if (NULL == completed[o])
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom " << o << " has not been defined\n";
    return 0;
  }

  if (this == completed[o])
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom bonded to itself\n";
    return 0;
  }

  Substructure_Bond * b = new Substructure_Bond;

  int characters_processed;
  if (! b->construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
  {
    cerr << "parse_smarts_bond_attribute:cannot parse bond smarts " << *att << "\n";
    delete b;
    return 0;
  }

  if (characters_processed != bond_smarts.length())
  {
    cerr << "parse_smarts_bond_attribute:extra junk at end of bond smarts " << *att << "\n";
    cerr << characters_processed << " characters processed\n";
    delete b;
    return 0;
  }

  b->set_atom(completed[o]);

  return _add_bond(b);
}

/*
  An attribute bond will have the other end as the first int value
*/

int
Substructure_Atom::_process_attribute_bond (const msi_attribute * att,
                                            bond_type_t bond_type,
                                            extending_resizable_array<Substructure_Atom *> & completed)
{
//cerr << "Parsing bond from " << (*att) << " bt " << bond_type << endl;

  int other_end;
  if (1 == att->number_int_values())
    (void) att->value(other_end);
  else
    other_end = att->int_multi_value(0);

  if (other_end < 0)
  {
    cerr << "Invalid other end " << other_end << endl;
    return 0;
  }

  if (NULL == completed[other_end])
  {
    cerr << "Substructure_Atom::_process_attribute_bond: atom " << other_end << " has not been defined\n";
    return 0;
  }

  if (this == completed[other_end])
  {
    cerr << "Substructure_Atom::_add_attribute_bond: atom bonded to itself\n";
    return 0;
  }

  Substructure_Bond * b = new Substructure_Bond;
  b->set_atom(completed[other_end]);

  if (att->number_int_values() > 1)
  {
    bond_type = 0;

    for (int i = 1; i < att->number_int_values(); i++)
    {
      int j = att->int_multi_value(i);

      if (1 == j)
        bond_type = bond_type | SINGLE_BOND;
      else if (2 == j)
        bond_type = bond_type | DOUBLE_BOND;
      else if (3 == j)
        bond_type = bond_type | TRIPLE_BOND;
      else if (4 == j)
        bond_type = bond_type | AROMATIC_BOND;
      else
      {
        cerr << "Substructure_Atom::_add_attribute_bond: what kind of bond is " << j << endl;
        return 0;
      }
    }
  }

  b->set_type(bond_type);

  return _add_bond(b);
}

/*
  Attribute bonds for a Substructure_Environment object are
  different because they MAY have multiple points of attachment.
*/

int
Substructure_Environment::_process_attribute_bond (const msi_attribute * att,
                                            bond_type_t bond_type,
                                            extending_resizable_array<Substructure_Atom *> & completed)
{
//cerr << "Parsing '" << *att << "' btype " << bond_type << " as bond\n";

  if (0 != bond_type)     // only add if bond_type has defined a particular type
    _bond.set_type(bond_type);
  else
    _bond.set_match_any();

  int nb = att->number_int_values();
  if (1 == nb)         // NOT multi-valued
  {
    int a;
    if (! att->value(a) || a < 0)
    {
      cerr << "Substructure_Environment::_process_attribute_bond: invalid atom '" << (*att) << endl;
      return 0;
    }

    if (NULL == completed[a])
    {
      cerr << "Substructure_Environment::_process_attribute_bond: atom " << a << " not defined\n";
      return 0;
    }

    _possible_parents.add(completed[a]);

    return 1;
  }

// Multiply valued

  for (int i = 0; i < nb; i++)
  {
    int a = att->int_multi_value(i);

    if (a < 0)
    {
      cerr << "Invalid other end " << a << endl;
      return 0;
    }

    if (NULL == completed[a])
    {
      cerr << "Substructure_Environment::_process_attribute_bond: atom " << a << " has not been defined\n";
      return 0;
    }

    _possible_parents.add(completed[a]);
  }

  return 1;
}

/*
  Note that bonds as read in can be 0 1 2 3 or 4 only
  where 0 means unknown bond type
  and   4 means aromatic bond type

  Also allow things like

    (A I single|double|triple|aromatic_bond 0)
*/

static int
attribute_is_bond (const msi_attribute * att,
                   bond_type_t & result)
{
  if (NAME_OF_BOND_ATTRIBUTE == att->name())
  {
    result = SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND;
    return 1;
  }

  if (NAME_OF_SINGLE_BOND_ATTRIBUTE == att->name())
  {
    result = SINGLE_BOND;
    return 1;
  }

  if (NAME_OF_DOUBLE_BOND_ATTRIBUTE == att->name())
  {
    result = DOUBLE_BOND;
    return 1;
  }

  if (NAME_OF_TRIPLE_BOND_ATTRIBUTE == att->name())
  {
    result = TRIPLE_BOND;
    return 1;
  }

  if (NAME_OF_AROMATIC_BOND_ATTRIBUTE == att->name())
  {
    result = AROMATIC_BOND;
    return 1;
  }

  IWString a = att->name();
//cerr << "What about '" << a << "'\n";

  if (! a.ends_with("_bond"))
    return 0;

  a.chop(5);

  const_IWSubstring token;
  int i = 0;
  result = 0;
  while (a.nextword(token, i, '|'))
  {
//  cerr << "Check '" << token << "'\n";
    if ("single" == token)
      result |= SINGLE_BOND;
    else if ("double" == token)
      result |= DOUBLE_BOND;
    else if ("triple" == token)
      result |= TRIPLE_BOND;
    else if ("aromatic" == token)
      result |= AROMATIC_BOND;
    else
    {
      return 0;
    }
  }

//cerr << "Final value " << result << endl;
  return 1;
}

int
Substructure_Atom::_process_attribute_bonds (const msi_object & msi, 
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  const msi_attribute * att;
  int i = 0;
  while (NULL != (att = msi.attribute(i++)))
  {
    if (NAME_OF_BOND_SMARTS_ATTRIBUTE == att->name())
    {
      if (! _process_attribute_smarts_bond(att, completed))
        return 0;
    }

    bond_type_t bond_type;

    if (! attribute_is_bond(att, bond_type))
      continue;

//  cerr << "Bond '" << *att;

    if (! _process_attribute_bond(att, bond_type, completed))
      return 0;
  }

  return 1;
}

/*
  The environment has just one bond. If you want the environment
  multiply attached, add a ring bond to the first atom in the env
*/

int
Substructure_Environment::_process_attribute_bonds (const msi_object & msi,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  int rc = 0;

  const msi_attribute * att;
  int i = 0;
  while (NULL != (att = msi.attribute(i++)))
  {
//  cerr << "Substructure_Environment::_process_attribute_bonds: examining " << *att << endl;

    bond_type_t bond_type;

    if (! attribute_is_bond(att, bond_type))
      continue;

    if (rc)
    {
      cerr << "Substructure_Environment::_process_attribute_bonds: multiple bonds\n";
      cerr << "The environment can have just one bond attachment\n";
      return 0;
    }

//  cerr << "Processing '" << (*att) << "' as attribute bond\n";
    if (! _process_attribute_bond(att, bond_type, completed))
      return 0;

    rc++;
  }

  if (rc)
    return rc;

// If no bond directives, maybe a bond smarts directive

  i = 0;
  while (NULL != (att = msi.attribute(i++)))
  {
    if (NAME_OF_BOND_SMARTS_ATTRIBUTE != att->name())
      continue;

//  cerr << "Examining " << *att << endl;
    Set_of_Atoms other_end;
    IWString bond_smarts;

    if (! parse_smarts_bond_attribute(att, other_end, bond_smarts))
    {
      cerr << "Substructure_Environment::_process_attribute_bond:invalid bond smarts " << *att << endl;
      return 0;
    }

    int characters_processed;
    if (! _bond.construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
    {
      cerr << "Substructure_Environment::_process_attribute_bond:invalid bond smarts " << bond_smarts << endl;
      return 0;
    }

    for (int j = 0; j < other_end.number_elements(); ++j)
    {
      const atom_number_t a = other_end[j];

      if (NULL == completed[a])
      {
        cerr << "Substructure_Environment::_process_attribute_bonds: atom " << a << " has not been defined\n";
        return 0;
      }

      _possible_parents.add(completed[a]);
    }
  }

  return 1;
}

/*
  When creating an environment, the environment components will
  register themselves as children. They are really environment components,
  so after they are done, we move any extra children to components.
*/

int
Substructure_Atom::_create_environment_from_msi_object (const msi_object & msi,
                                           extending_resizable_array<Substructure_Atom *> & completed)
{
  assert (0 == _environment.number_elements());

  int initial_children = _children.number_elements();

  int nmsi = msi.number_elements();

  for (int i = 0; i < nmsi; i++)
  {
    const msi_object * m = msi[i];
    if (NAME_OF_QUERY_ATOM_OBJECT != m->name())
    {
      cerr << "_create_environment_from_msi_object: unknown type " << (*m);
      return 0;
    }

    Substructure_Atom * a = new Substructure_Atom;
    if (! a->construct_from_msi_object(*m, completed))
    {
      cerr << "_create_environment_from_msi_object: cannot create environment component " << i << endl;
      delete a;
      return 0;
    }
  }

  int final_children = _children.number_elements();
  if (initial_children == final_children)
  {
    cerr << "create_env_from_msi: children not attached to root!\n";
    cerr << msi;
    return 0;
  }

  while (_children.number_elements() > initial_children)
    _environment.transfer_in(_children, initial_children);

  return _environment.number_elements();
}

/*
  This does not do anything?
*/

int
Substructure_Atom_Environment::create_from_msi_object (msi_object & msi)
{
  IWString smarts;
  if (msi.string_value_for_attribute(NAME_OF_SMARTS_ATTRIBUTE, smarts))
  {
    if (msi.number_elements())
    {
      cerr << "Substructure_Atom_Environment::create_from_msi_object: smarts and multiple objects specified, impossible\n"; return 0;
    }
  }

  int n = msi.number_elements();
  for (int i = 0; i < n; i++)
  {
    const msi_object * oi = msi[i];
    if (NAME_OF_QUERY_ATOM_OBJECT != oi->name())
    {
      cerr << "Substructure_Atom_Environment::create_from_msi_object: unrecognised object\n";
      cerr << (*oi);
      return 0;
    }
  }

  return 1;
}
 
/*
  The sub-objects of a Substructure_Atom will be either bonds, or
  the Substructure_Atoms which are its components.
  Note that components are not allowed to have bonds.
*/

int
Substructure_Atom::_construct_from_msi_object (const msi_object & msi,
                                               extending_resizable_array<Substructure_Atom *> & completed)
{
  if (NAME_OF_QUERY_ATOM_OBJECT != msi.name())
  {
    cerr << "Substructure_Atom::_construct_from_msi_object: name must be '" << NAME_OF_QUERY_ATOM_OBJECT << "', '" << msi.name() << "' is invalid\n";
    return 0;
  }

#ifdef DEBUG_COMPLTED
  cerr << "Building atom\n";msi->print(cerr);
  for (int i = 0; i < 6; i++)
  {
    cerr << "Completed[" << i << "] = " << completed[i] << endl;
  }

  cerr << "Building atom from\n" << (*msi);
#endif

  if (! Substructure_Atom_Specifier::construct_from_msi_object(msi))
  {
    cerr << "Substructure_Atom::construct_from_msi_object: cannot get attributes\n";
    return 0;
  }

  _unique_id = msi.object_id();

  if (NULL != completed[_unique_id])
  {
    cerr << "Substructure_Atom::_construct_from_msi_object:Yipes, atom " << _unique_id << " already allocated\n";
    cerr << msi;
    return 0;
    iwabort();
  }

  const msi_attribute * att = msi.attribute(NAME_OF_INITIAL_ATOM_NUMBER_ATTRIBUTE);
  if (NULL != att)
  {
    if (! att->value(_initial_atom_number) || _initial_atom_number < 0)
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: initial atom number must be a non negative number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_OR_OPERATOR);
  if (att)
  {
    if (! att->value(_or_id) || 0 == _or_id)
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: or id must be a positive number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_RING_ID_SPECIFIER);
  if (att)
  {
    if (! att->value(_ring_id) || 0 == _ring_id)
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: ring id must be a positive number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_FRAGMENT_ID_SPECIFIER);
  if (att)
  {
    if (! att->value(_fragment_id) || 0 == _fragment_id)
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: fragment id must be a positive number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_TEXT_IDENTIFIER_ATTRIBUTE);
  if (att)
    att->value(_text_identifier);

  att = msi.attribute(NAME_OF_NUMERIC_VALUE_ATTRIBUTE);
  if (att)
  {
    double tmp;
    if (! att->value(tmp))
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: the " << NAME_OF_NUMERIC_VALUE_ATTRIBUTE << " attribute must be numeric\n";
      return 0;
    }

    _numeric_value.set(tmp);
  }

  att = msi.attribute(NAME_OF_INCLUDE_IN_EMBEDDING_ATTRIBUTE);
  if (att)
  {
    int tmp;
    if (! att->value(tmp))
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: the " <<
             NAME_OF_INCLUDE_IN_EMBEDDING_ATTRIBUTE << " attribute requires a whole number\n";
      return 0;
    }

    _include_in_embedding = tmp;
  }

// For simplicity, one can specify only one of atom smarts, smiles and smarts

  if ((NULL !=  msi.attribute(NAME_OF_ATOM_SMARTS_ATTRIBUTE)) +
      (NULL !=  msi.attribute(NAME_OF_SMARTS_ATTRIBUTE)) +
      (NULL !=  msi.attribute(NAME_OF_SMILES_ATTRIBUTE)) > 1)
  {
      cerr << "The '" << NAME_OF_ATOM_SMARTS_ATTRIBUTE << "', '" <<
              NAME_OF_SMILES_ATTRIBUTE << "' and '" <<
              NAME_OF_SMILES_ATTRIBUTE << "' attributes are mutually exclusive\n";
      return 0;
  }

  att = msi.attribute(NAME_OF_ATOM_SMARTS_ATTRIBUTE);
  if (att)
  {
    IWString csmarts;
    att->value(csmarts);

    if (! (csmarts.starts_with('[') && csmarts.ends_with(']')))
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: atomic smarts must be within [] '" << csmarts << "'\n";
      return 0;
    }

    if (! construct_from_smarts_token(csmarts))
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: cannot interpret smarts '" << csmarts << "'\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_SUM_ALL_PREFERENCE_HITS_ATTRIBUTE);
  if (att)
  {
    if (! att->value(_sum_all_preference_hits))
    {
      cerr << "Substructure_Atom::_construct_from_msi_object: the " << NAME_OF_SUM_ALL_PREFERENCE_HITS_ATTRIBUTE << " attribute must be a whole number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_SMILES_ATTRIBUTE);
  if (att)
  {
    if (! parse_smiles_specifier(att))
      return 0;
  }

  att = msi.attribute(NAME_OF_SMARTS_ATTRIBUTE);
  if (att)
  {
    if (! parse_smarts_specifier(att))
      return 0;
  }

  completed[_unique_id] = this;

  int nmsi = msi.number_elements();
  for (int i = 0; i < nmsi; i++)
  {
    msi_object & m = *(msi[i]);

    if (NAME_OF_QUERY_ATOM_OBJECT == m.name() || 
             NAME_OF_QUERY_ATOM_SPECIFIER_OBJECT == m.name())
    {
      if (! _add_component(m))
        return 0;
    }
    else if (NAME_OF_ENVIRONMENT_OBJECT == m.name())
    {
      if (! _create_environment_from_msi_object(m, completed))
        return 0;
    }
    else if (NAME_OF_QUERY_ATOM_PREFERENCE_OBJECT == m.name())
    {
      if (! _create_preference_from_msi_object(m))
        return 0;
    }
    else if (NAME_OF_QUERY_BOND_OBJECT == m.name())   // will be processed by calling function
      ;
    else
    {
      cerr << "Substructure_Atom::construct_from_msi_object: unrecognised type '" << 
             m.name() << "'\n";
      cerr << (NAME_OF_QUERY_BOND_OBJECT == m.name()) << endl;
      return 0;
    }
  }
      
  assert (ok());
  return 1;
}

/*
  This is mostly a wrapper to _construct_from_msi_object.
  If it fails, we print out the offending msi_object
  Because Substructure_Environment objects have different
  bonding types, the attribute bonds are processed here, leaving
  _construct_from_msi_object common to both classes.

  When processing query atoms, whenever a bond is found, the
  atom at the other end of the bond must already have been created
  and in the COMPLETED array. We notify the parent of the existence
  of another child for each first bond.
*/

int
Substructure_Atom::construct_from_msi_object(const msi_object & msi,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (is_rejection(msi))
    _match_as_match_or_rejection = 0;

  if (! _construct_from_msi_object(msi, completed))
  {
    cerr << msi;
    return 0;
  }

  if (_initial_atom_number < 0)
    _initial_atom_number = msi.object_id();

// Process any bonds which have been specified as attributes

  if (! _process_attribute_bonds(msi, completed))
  {
    cerr << msi;
    return 0;
  }

  return 1;
}

/*
  Does the ring consist of anything other than alternating single and double bonds
*/

int
is_fixed_kekule_form (const Molecule & m,
                      const Ring & r)
{

  bond_type_t previous_bond_was = 0;
  bond_type_t first_bond_encountered = 0;

//cerr << "Checking ring " << r << endl;
  for (Ring_Bond_Iterator j (r); j != r.end(); j++)
  {
    atom_number_t a1 = j.a1();
    atom_number_t a2 = j.a2();

    const Bond * b = m.bond_between_atoms(a1, a2);

    bond_type_t bt;
    if (b->is_single_bond())
      bt = SINGLE_BOND;
    else
      bt = DOUBLE_BOND;

//  cerr << "From " << a1 << " to " << a2 << " bt " << bt << " prev " << previous_bond_was << endl;
    if (0 == previous_bond_was)
      first_bond_encountered = bt;
    else if (previous_bond_was == bt)
      return 1;

    previous_bond_was = bt;
  }

// If the last bond we encountered is the same as the first we encountered,
// this must be a fixed system. Note we could work this out from the size
// of the ring, but this is easier

  return first_bond_encountered == previous_bond_was;   // if same, then fixed form
}

/*
  Is the bond B in a non-Kekule aromatic ring
*/

static int
determine_part_of_non_kekule_aromatic (Molecule & m,
                                       const Bond & b)
{
  m.compute_aromaticity_if_needed();

  int nr = m.nrings();

//cerr << "Checking '" << m.name() << "' for fixed kekule aromatic, nr = " << nr << endl;
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

//  cerr << "Ring " << i << " aromatic? " << ri->is_aromatic() << endl;
    if (! ri->is_aromatic())
      continue;

    int ring_size = ri->number_elements();

    if (ring_size > 6)
      return 0;

    if (! ri->contains_bond(b.a1(), b.a2()))
      continue;

    if (5 == ring_size)   // always
      return 1; 

    if (is_fixed_kekule_form(m, *ri))
    {
      cerr << "In '" << m.name() << " bond btw " << b.a1() << " and " << b.a2() << " in fixed kekule ring of size " << ring_size << endl;
      return 1;
    }
  }

  return 0;
}

static int
initialise_bond_attributes (Substructure_Bond & sb,
                            const MDL_Bond_Data * mdlbd,
                            const Bond * b,
                            atom_number_t j,
                            extending_resizable_array<Substructure_Atom *> & completed,
                            const Molecule_to_Query_Specifications & mqs,
                            int part_of_non_kekule_aromatic)
{
//#define DEBUG_INITIALISE_BOND_ATTRIBUTES
#ifdef DEBUG_INITIALISE_BOND_ATTRIBUTES
  cerr << "atoms " << b->a1() << " to " << b->a2() << " mdlbd->btype() " << mdlbd->btype() << " aromatic? " << b->is_aromatic() << endl;
#endif

  if (mqs.just_atomic_number_and_connectivity())
    sb.set_match_any();
  else if (0 == mdlbd->btype())
    sb.copy(b, mqs.copy_bond_attributes());
  else if (NULL == b)
    sb.set_match_any();
  else if (b->is_aromatic())
  {
    int ablki = aromatic_bonds_lose_kekule_identity();
    if (0 == ablki)
    {
      bond_type_t bt = b->btype();
      sb.set_type(bt & BOND_TYPE_ONLY_MASK);
    }
    else if (1 == ablki)
      sb.set_type(AROMATIC_BOND);
    else if (part_of_non_kekule_aromatic)
    {
      bond_type_t bt = b->btype();
      sb.set_type(bt & BOND_TYPE_ONLY_MASK);
    }
    else
      sb.set_type(AROMATIC_BOND);
  }
  else
    sb.set_type(mdlbd->btype());

#ifdef DEBUG_INITIALISE_BOND_ATTRIBUTES
  cerr << "Set bond type " << sb.types_matched() << endl;
#endif

  sb.set_atom(completed[j]);

  if (0 == mdlbd->bond_topology())   // nothing known
    ;
  else if (1 == mdlbd->bond_topology())
    sb.set_must_be_in_a_ring(1);
  else if (2 == mdlbd->bond_topology())
  {
    sb.set_must_be_in_a_ring(0);
//  cerr << "Bond topology 2 between atoms " << b->a1() << " and " << b->a2() << endl;
  }
  else
    cerr << "Unrecognised bond topology for atom " << b->other(j) << " " << mdlbd->bond_topology() << endl;

  if (mqs.bonds_preserve_ring_membership())   // we do not test for conflicting directives
  {
    if (b->nrings())
      sb.set_must_be_in_a_ring(1);
    else
      sb.set_must_be_in_a_ring(0);
  }

  return 1;
}

int
Substructure_Atom::_build_atomic_number_matches_from_atom_list(const ISIS_Atom_List & ali)
{
  if (ali.normal_list())
  {
    for (int j = 0; j < ali.number_elements(); j++)
    {
      const Element * ej = ali.elementi(j);

      set_element(ej);
    }

    return 1;
  }

  Substructure_Atom_Specifier * s = new Substructure_Atom_Specifier;

  for (int i = 0; i < ali.number_elements(); i++)
  {
    s->set_element(ali.elementi(i));
  }

  _components.add(s);

  int n = _components.number_elements();
  _operator.set_unary_operator(n - 1, 0);

  if (n > 1)
    _operator.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);

  return 1;
}

static int
count_heteroatoms_in_aromatic_ring (Molecule & m,
                                    atom_number_t zatom)
{
  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (! ri->contains(zatom))
      continue;

    int tmp = m.count_heteroatoms(*ri);

    if (tmp > rc)
      rc = tmp;
  }

  return rc;
}

/*
  We have an aromatic atom and the molecule to query specifier says that we
  should match any aromatic atom.
  The convention is that if convert_all_aromatic_atoms_to_generic_aromatic is
    1: any aromatic
    2: aromatic, heteroatom or carbon demepding on my_atom_number
    3. aromatic, same number of heteroatoms in ring
*/

int
Substructure_Atom::_set_match_any_aromatic (int convert_all_aromatic_atoms_to_generic_aromatic, 
                                            Molecule & m,
                                            atom_number_t my_atom_number)
{
  _element.resize(0);
  _aromaticity = AROMATIC;
  _hcount.reset();

  if (1 == convert_all_aromatic_atoms_to_generic_aromatic)
    return 1;

  if (2 == convert_all_aromatic_atoms_to_generic_aromatic)
  {
    if (6 == m.atomic_number(my_atom_number))
      _element.add(get_element_from_atomic_number(6));
    else
    {
      _element.add(get_element_from_atomic_number(7));
      _element.add(get_element_from_atomic_number(8));
      _element.add(get_element_from_atomic_number(16));
    }

    return 1;
  }

  else if (3 == convert_all_aromatic_atoms_to_generic_aromatic)
  {
    if (1 == m.nrings(my_atom_number))
      _heteroatoms_in_ring.add(count_heteroatoms_in_aromatic_ring(m, my_atom_number));
    return 1;
  }

  return 1;
}

//#define DEBUG_CREATE_FROM_MOLECULE

/*
  This is complicated by the need to set either ncon or min_ncon depending
  on whether there are implicit hydrogens in the molecule
*/

int
Substructure_Atom::create_from_molecule (Molecule & m,
                                         const MDL_File_Data & mdlfd,
                                         atom_number_t my_atom_number,
                                         Substructure_Atom * my_parent,
                                         const Molecule_to_Query_Specifications & mqs,
                                         extending_resizable_array<Substructure_Atom *> & completed,
                                         const int * include_these_atoms)
{
//cerr << "Substructure_Atom::creating from molecule with " << m.natoms() << " atoms, atom " << my_atom_number << endl;

  m.compute_aromaticity_if_needed();

//cerr << "At atom " << my_parent << endl;
//m.debug_print(cerr);

  assert (my_atom_number < m.natoms());

  _parent = my_parent;

  atom_number_t parent_atom_number;
  if (_parent)
    parent_atom_number = _parent->unique_id();
  else
    parent_atom_number = INVALID_ATOM_NUMBER;

  _unique_id = my_atom_number;
  _initial_atom_number = my_atom_number;

  assert (NULL == completed[my_atom_number]);

  completed[my_atom_number] = this;

#ifdef DEBUG_CREATE_FROM_MOLECULE
  cerr << "create from molecule atom " << my_atom_number << " unique id " << _unique_id;
  if (INVALID_ATOM_NUMBER != parent_atom_number)
    cerr << ". Parent " << parent_atom_number;
  cerr << endl;
#endif

  const Atom * a = m.atomi(my_atom_number);

  const Element * e = a->element();

  const MDL_Atom_Data * mdlad = mdlfd.mdl_atom_data(my_atom_number);

  const IWString initial_atomic_symbol = mdlad->initial_atomic_symbol();

  if (initial_atomic_symbol.length())
    e = get_element_from_symbol_no_case_conversion(initial_atomic_symbol);

  IWString smarts_specified;

//cerr << "Atom alias for " << my_atom_number << " " << mdlad->alias() << endl;
//cerr << "Examining atom " << mdlad->initial_atomic_symbol() << endl;
// Jan 2009. If the alias is a number, don't interpret it as a smarts

  if (mqs.interpret_atom_alias_as_smarts() && mdlad->alias().length())
  {
    const_IWSubstring s = mdlad->alias();
    int notused;
    if (s.numeric_value(notused) && notused >= 0)
      ;
    else
    {
      smarts_specified = mdlad->alias();

      if (! smarts_specified.starts_with('[') && smarts_specified.ends_with(']'))
      {
        cerr << "Substructure_Atom::create_from_molecule:invalid atomic smarts as alias '" << smarts_specified << "'\n";
        return 0;
      }
    }
  }
  else if (! mqs.interpret_atom_alias_as_smarts() && mdlad->alias().length())   // probably from Marvin, AH, AQ...
    _build_atomic_number_matches_from_atom_list(mdlad->atom_list());
  else if (mqs.smarts_for_element(e, smarts_specified))
  {
    if (! smarts_specified.starts_with('[') && smarts_specified.ends_with(']'))
    {
      cerr << "Substructure_Atom::create_from_molecule:invalid atomic smarts as alias '" << smarts_specified << "'\n";
      return 0;
    }
  }
  else if (e->is_in_periodic_table())
  {
    if (1 == e->atomic_number() && (mqs.convert_explicit_hydrogens_to_match_any_atom() || mqs.convert_explicit_hydrogens_to_match_any_atom_including_hydrogen()))
      ;
    else
      set_element(e);
  }
  else if ("A" == initial_atomic_symbol)    // match any atom
    ;
  else if (mdlad->atom_list().active())
    _build_atomic_number_matches_from_atom_list(mdlad->atom_list());
  else if ("L" == initial_atomic_symbol)   // looks like an ISIS atom list - should be handled previously
    ;
  else if ("Q" ==  initial_atomic_symbol)    // heteroatom
  {
    Substructure_Atom_Specifier * q = new Substructure_Atom_Specifier;
    q->set_element(get_element_from_atomic_number(6));
    _components.add(q);
    if (_components.number_elements() > 1)   // will never happen here
      _operator.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
    _operator.set_unary_operator(0, 0);
  }

  int acon = a->ncon();

#ifdef DEBUG_CREATE_FROM_MOLECULE
  cerr << "Atom " << my_atom_number << " type " << m.atomic_number(my_atom_number) << " has " << m.hcount(my_atom_number) << " hydrogens\n";
#endif

  if (smarts_specified.length() > 0)
    ;
  else if (1 == mdlad->exact_change())
  {
    _ncon.add(acon);
    _nbonds.add(a->nbonds());
  }
  else if (mqs.query_must_match_both_explicit_and_implicit_hydrogens())
    cerr << "query_must_match_both_explicit_and_implicit_hydrogens\n";
  else if ((MQS_ISOTOPE_MEANS_NCON & mqs.isotopic_label_means()) && a->isotope() > 0)
    _ncon.add(m.ncon(my_atom_number));
  else if (0 != mdlad->substitution())
  {
//  cerr << "Atom " << my_atom_number << " substitution " << mdlad->substitution() << endl;
    if (-2 == mdlad->substitution())
      _ncon.add(acon);
    else if (-1 == mdlad->substitution())
      _ncon.add(0);
    else if (mdlad->substitution() > 0)
      _ncon.add(mdlad->substitution());
    else if (-3 == mdlad->substitution())    // IAW extentions that means do nothing
      ;
    else
      set_min_ncon(acon);
  }
  else if (0 != mdlad->min_ncon())
    set_min_ncon(mdlad->min_ncon());
  else if (0 != mdlad->max_ncon())
    set_max_ncon(mdlad->max_ncon());
  else if (! mqs.make_embedding())
  {
    _ncon.add(acon);
    _nbonds.add(a->nbonds());
  }
  else
  {
    if (1 == acon && NULL != my_parent)    // doesn't make sense to specify anything
      ;
    else
      set_min_ncon(acon);
  }

//#define DEBUG_HCOUNT_STUFF
#ifdef DEBUG_HCOUNT_STUFF
  cerr << "Processing atom " << my_atom_number << endl;
  cerr << " hcount " << mdlad->hcount();
  cerr << " h0designator " << mdlad->h0designator();
  cerr << " min_hcount " << mdlad->min_hcount();
  cerr << endl;
#endif

  if (smarts_specified.length())
    ;
  else if (! mqs.discern_hcount())
    ;
  else if (1 == mdlad->exact_change())
    _hcount.add(const_cast<Atom *>(a)->implicit_hydrogens());
  else if (1 == mdlad->h0designator())
    _hcount.add(0);
  else if (mdlad->hcount() > 0)   // hcount 1 in ISIS means NO implicit Hydrogens
    _hcount.add(mdlad->hcount() - 1);
  else if (mdlad->min_hcount() > 0)
    _hcount.set_min(mdlad->min_hcount());
  else if (mdlad->explicit_hydrogen_atoms_removed())
    _hcount.set_min(mdlad->explicit_hydrogen_atoms_removed());
  else if (! a->element()->organic())    // cannot say anything about it
    ;
  else if (! e->is_in_periodic_table())    // can't say anything about hcount
    ;
  else if (1 == e->atomic_number())   // don't say anything about hcount
    ;
  else if (mqs.just_atomic_number_and_connectivity())
    ;
  else if (mqs.ignore_molecular_hydrogen_information())
    ;
  else if (_ncon.is_set() || -3 == mdlad->substitution())
    ;
  else
  {
    int exh = m.explicit_hydrogens(my_atom_number);
    int ih  = m.implicit_hydrogens(my_atom_number);

    if (0 == exh && 0 == ih)
      _hcount.add(0);
    else if (exh + ih < 3)       // specifying a max of 3 Hydrogens is useless
      _hcount.set_max(exh + ih);
  }

#ifdef DEBUG_HCOUNT_STUFF
  cerr << "atom " << my_atom_number << " hcount set to ";
  _hcount.debug_print(cerr);
#endif

//cerr << "preserve_saturation " << mqs.preserve_saturation() << ", uns " << mdlad->unsaturated() << endl;

  if (smarts_specified.length())
    ;
  else if (mqs.preserve_saturation())
  {
    const auto s = a->nbonds() - a->ncon();
    _unsaturation.add(s);
  }
  else if (0 == mdlad->unsaturated())
    ;
  else if (mqs.just_atomic_number_and_connectivity())
    ;
  else if (mdlad->unsaturated() > 0)
  {
    if (8 == e->atomic_number())   // Oxygen can only have 1 unsaturation
      _unsaturation.add(1);
    else
      _unsaturation.set_min(1);
  }
  else if (SPECIAL_MEANING_FULLY_SATURATED == mdlad->unsaturated())
    _unsaturation.add(0);

  if (smarts_specified.length())
    ;
  else if (mqs.just_atomic_number_and_connectivity())
    ;
  else if (a->formal_charge())
    set_formal_charge(a->formal_charge());

  if (smarts_specified.length())
    ;
  else if (only_include_isotopically_labeled_atoms())  // ignore isotopic labels
    ;
  else if (mqs.substituents_only_at_isotopic_atoms())   // ignore isotopes
    ;
  else if (0 != mqs.isotopic_label_means())
    ;
  else if (a->isotope())
    _isotope.add(a->isotope());
    
  if (a->userAtomType() > 0)      
    _userAtomType  = a->userAtomType();

  int nr = m.nrings(my_atom_number);

  int rbc = m.ring_bond_count(my_atom_number);

//cerr << "Atom " << my_atom_number << " nr = " << nr << " rbc " << rbc << endl;

  if (smarts_specified.length())
    ;
  else if (0 != mdlad->ring_bond())
  {
    int r = mdlad->ring_bond();
//  cerr << "Ring stuff for " << _initial_atom_number << " r = " << r << ", nr = " << nr << ", rbc " << rbc << endl;
    if (-1 == r)
      _ring_bond_count.add(0);
    else if (-2 == r)  // as drawn
      _ring_bond_count.add(rbc);
    else if (-3 == r)    // ignore
      ;
    else if (r > 1)
      _ring_bond_count.add(r);
  }
  else if (mqs.atoms_conserve_ring_membership() && nr > 0)
  {
    _ring_bond_count.add(rbc);
//  int ns = m.nrings_including_non_sssr_rings(my_atom_number);
//  cerr << "atoms_conserve_ring_membership, nr = " << nr << " ns = " << ns << endl;

//  for (int i = nr; i <= ns; i++)
//  {
//    _nrings.add(i);
//  }
  }
  else if (mqs.ring_atoms_conserve_ring_membership() && nr)
  {
    _nrings.set_min(nr);
  }
  else if (0 == nr && mqs.non_ring_atoms_become_nrings_0())
    _nrings.add(0);
  else if ((MQS_ISOTOPE_MEANS_RING_BOND_COUNT & mqs.isotopic_label_means()))
    _ring_bond_count.add(rbc);
  else if (rbc > 0)
    _ring_bond_count.set_min(rbc);
  else if (nr > 0)
    set_min_nrings(nr);

//cerr << "For atom " << my_atom_number << " nr = " << nr << " arom " << m.is_aromatic(my_atom_number) << endl;

  if (smarts_specified.length())
    ;
  else if (mqs.just_atomic_number_and_connectivity())
    ;
  else if (mqs.convert_all_aromatic_atoms_to_generic_aromatic() && 
           (m.is_aromatic(my_atom_number) || mdlad->aromatic()))
  {
    _set_match_any_aromatic(mqs.convert_all_aromatic_atoms_to_generic_aromatic(), m, my_atom_number);
  }
  else if (mqs.aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic())
  {
    if (m.is_aromatic(my_atom_number))
      update_aromaticity(AROMATIC);
    else
      update_aromaticity(NOT_AROMATIC);
  }
  else if (nr)
  {
    int known_aromaticity = mdlad->aromatic();

    if (known_aromaticity < 0)    // need to compute it
    {
      aromaticity_type_t arom;
      if (m.aromaticity(my_atom_number, arom))
      {
        if (IS_AROMATIC_ATOM(arom))
          update_aromaticity(arom);
        else if (mqs.only_aromatic_atoms_match_aromatic_atoms())
          update_aromaticity(NOT_AROMATIC);
//      if (IS_ALIPHATIC_ATOM(arom))     // only make positive decisions about aromaticity
//        update_aromaticity(arom);
      }
    }
    else if (AROMATIC == known_aromaticity)
      update_aromaticity(AROMATIC);
    else if (NOT_AROMATIC == known_aromaticity)
      update_aromaticity(NOT_AROMATIC);
  }
  else if (m.is_aromatic(my_atom_number))   // chain aromatic atom
    update_aromaticity(AROMATIC);
  else if (mqs.only_aromatic_atoms_match_aromatic_atoms())
  {
    if (m.is_aromatic(my_atom_number))
      update_aromaticity(AROMATIC);
    else
      update_aromaticity(NOT_AROMATIC);
  }

// We take the easy way out here and just use our component. Probably should do all of them if present

//cerr << "Building query " << mqs.use_preference_values_to_distinguish_symmetry() << endl;

  if (! mqs.use_preference_values_to_distinguish_symmetry())
    ;
  else if (1 == e->atomic_number())
    ;
  else
  {
    Set_of_Atoms s;
    m.symmetry_equivalents(my_atom_number, s);
    if (s.number_elements())
      _add_ncon_query_atom_preferences(m, my_atom_number);
  }

  if (smarts_specified.length())
    ;
  else if (mdlad->nbonds() >= 0)
    _nbonds.add(mdlad->nbonds());

// Because of cyclic structures, we need to be very careful how we build the query.
// Add the bond to the parent first

  if (INVALID_ATOM_NUMBER != parent_atom_number)
  {
    Substructure_Bond * sb = new Substructure_Bond;

    const Bond * b = m.bond_between_atoms(my_atom_number, parent_atom_number);
    assert (NULL != b);

    int bnumber = m.which_bond(my_atom_number, parent_atom_number);

//  When dealing with Hydrogens that can be anything, we need to allow the bond type to be anything. Bit of a kludge here I'm afraid

    if (1 == e->atomic_number() && (mqs.convert_explicit_hydrogens_to_match_any_atom() || mqs.convert_explicit_hydrogens_to_match_any_atom_including_hydrogen()))
    {
      initialise_bond_attributes(*sb, mdlfd.mdl_bond_data(bnumber), NULL, parent_atom_number, completed, mqs, 0);
    }
    else if (bnumber >= mdlfd.number_bonds())   // underlying mdl_molecule may have been expanded
      ;
    else if (aromatic_bonds_lose_kekule_identity())
      initialise_bond_attributes(*sb, mdlfd.mdl_bond_data(bnumber), b, parent_atom_number, completed, mqs, determine_part_of_non_kekule_aromatic(m, *b));
    else
      initialise_bond_attributes(*sb, mdlfd.mdl_bond_data(bnumber), b, parent_atom_number, completed, mqs, 0);

#ifdef DEBUG_CREATE_FROM_MOLECULE
    sb->debug_print(cerr, "From parent ");
#endif

    _bonds.add(sb);

    assert (completed[parent_atom_number] == _parent);

    _bond_to_parent = sb;
    _parent = completed[parent_atom_number];

    if (smarts_specified.length())
    {
//  cerr << "Building from smarts '" << smarts_specified << "' (parent)\n";
      if (! construct_from_smarts_token(smarts_specified))
      {
        cerr << "Substructure_Atom::create_from_molecule:invalid smarts in alias '" << smarts_specified << "'\n";
        return 0;
      }
      smarts_specified.resize(0);
    }
  }

  if (smarts_specified.length())
  {
//  cerr << "Building from smarts '" << smarts_specified << "'\n";
    if (! construct_from_smarts_token(smarts_specified))
    {
      cerr << "Substructure_Atom::create_from_molecule:invalid smarts in alias '" << smarts_specified << "'\n";
      return 0;
    }
  }

// Now that the parent is done, do any ring closures

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(my_atom_number);

    if (j == parent_atom_number)
      continue;

    if (NULL == completed[j])     // not yet done, further down the tree
      continue;

    if (only_include_isotopically_labeled_atoms() && 0 == m.isotope(j))
      continue;

#ifdef DEBUG_CREATE_FROM_MOLECULE
    cerr << "From atom " << my_atom_number << " go to atom " << j << endl;
#endif

    int bnumber = m.which_bond(b->a1(), b->a2());

    Substructure_Bond * sb = new Substructure_Bond();

    if (aromatic_bonds_lose_kekule_identity())
      initialise_bond_attributes(*sb, mdlfd.mdl_bond_data(bnumber), b, j, completed, mqs, determine_part_of_non_kekule_aromatic(m, *b));
    else
      initialise_bond_attributes(*sb, mdlfd.mdl_bond_data(bnumber), b, j, completed, mqs, 0);

    _bonds.add(sb);

#ifdef DEBUG_CREATE_FROM_MOLECULE
    cerr << "Ring closure bond back to atom " << j << endl;
#endif
    assert (NULL != completed[j]);
  }

// Do the children now

  for (int i = 0; i < acon; i++)
  {
    const atom_number_t j = a->other(my_atom_number, i);

#ifdef DEBUG_CREATE_FROM_MOLECULE
    cerr << "From atom " << my_atom_number << " child " << j << " inc " << (NULL == include_these_atoms ? '1' : include_these_atoms[j]) << " complete " << completed[j] << ' ' << m.smarts_equivalent_for_atom(j) << endl;
#endif

    if (NULL != include_these_atoms && ! include_these_atoms[j])
      continue;

#ifdef DEBUG_CREATE_FROM_MOLECULE
    if (NULL != include_these_atoms)
      cerr << "include_these_atoms[" << j << "] " << include_these_atoms[j] << endl;
    else
      cerr << "include_these_atoms not set\n";
#endif

    if (completed[j])    // must be a cyclic structure
      continue;

    if (only_include_isotopically_labeled_atoms() && 0 == m.isotope(j))
      continue;

    Substructure_Atom * child = new Substructure_Atom;

    child->create_from_molecule(m, mdlfd, j, this, mqs, completed, include_these_atoms);

    _children.add(child);
  }

// Jun 2002.  Wow, Beckman rearrangement with MFCD00046298 breaks. 
// The reaction drawing has a 6 membered ring, but this reagent
// doesn't have a 6 membered ring - at least according to SSSR stuff. 
// Therefore we need to leave out ring sizes...

#ifdef USE_RING_SIZES
  List_of_Ring_Sizes lors;
  m.ring_sizes_for_atom(my_atom_number, lors);
  for (int i = 0; i < lors.number_elements(); i++)
  {
    _ring_size.add(lors[i]);
  }
#endif

  return 1;
}

int
Substructure_Atom::_add_ncon_query_atom_preferences(const Molecule & m,
                                                    atom_number_t my_atom_number)
{
  for (int i = m.ncon(my_atom_number) + 1; i <= 4; i++)
  {
    add_ncon_preference_object(i, 5 * my_atom_number * m.natoms() + i);
  }

  return 1;
}

int
Substructure_Atom::_parse_smiles_specifier (Molecule & m,
                                            const MDL_File_Data & mdlfd,
                                            extending_resizable_array<Substructure_Atom *> & completed)
{
  Molecule_to_Query_Specifications mqs;

  if (! create_from_molecule(m, mdlfd, 0, _parent, mqs, completed))
    return 0;

// We need to make room to allow ourselves to be attached to the rest
// of the query. Because we don't know how we will be bonded (if at
// all), we need to be conservative

  if (_nbonds.number_elements())
  {
    int nb = _nbonds[0];
    _nbonds.resize(0);
    _nbonds.set_min(nb + 1);
  }

  if (_ncon.number_elements())
  {
    int nc = _ncon[0];
    _ncon.resize(0);
    _ncon.set_min(nc + 1);
  }

  if (_hcount.number_elements())
  {
    int hc = _hcount[0];
    if (0 == hc)
    {
      cerr << "Substructure_Atom::_parse_smiles_specifier: warning, no open valence in '" << m.smiles() << "'\n";
      _hcount.resize(0);
    }
    else 
    {
      _hcount.resize(0);
      if (1 == hc)
        _hcount.add(0);
      else
        _hcount.set_max(hc - 1);
    }
  }

  return 1;
}

/*
  We have found a smiles attribute. This begins a complete molecule specification
*/

int
Substructure_Atom::parse_smiles_specifier (const IWString & smiles)
{
  assert (smiles.nchars());

  Molecule m;
  if (! m.build_from_smiles(smiles))
  {
    cerr << "Substructure_Atom::_parse_smiles_specifier: bad smiles\n";
    return 0;
  }

  extending_resizable_array<Substructure_Atom *> completed(NULL);
  completed.resize(m.natoms());

// Unique numbers get scrambled in the process, so save it

  int uid = _unique_id;

  MDL_File_Data mdlfd;
  mdlfd.build(m);

  int rc = _parse_smiles_specifier(m, mdlfd, completed);

  assign_unique_atom_numbers(uid);

  return rc;
}

int
Substructure_Atom::parse_smiles_specifier (const msi_attribute * msi)
{
  IWString smiles;
  msi->value(smiles);

  return parse_smiles_specifier(smiles);
}

int
Substructure_Atom::parse_smarts_specifier (const msi_attribute * msi)
{
  const_IWSubstring smarts;
  msi->value(smarts);

  if (smarts.contains('.'))
  {
    cerr << "Substructure_Atom::parse_smarts_specifier: (msi) cannot contain '.'\n";
    return 0;
  }

// Unique numbers get scrambled in the process, so save it

  int uid = _unique_id;

//cerr << "Calling parse smarts specifier on '" << smarts << "'\n";
  int rc = parse_smarts_specifier(smarts);

  assign_unique_atom_numbers(uid);

  return rc;
}

std::ostream &
operator << (std::ostream & os, const Substructure_Atom & a)
{
  os << "Atom " << a.unique_id();
  if (NULL == a.current_hold_atom())
    os << ", not matched\n";
  else
    os << ", matched to " << a.current_hold_atom()->atom_number();

  return os;
}

std::ostream &
operator << (std::ostream & os, const Link_Atom & l)
{
  os << "Link atom between " << l.a1() << " and " << l.a2();
  return os;
}

/*
  Do the best we can to create a molecule from a query.
  We use the atomic number of the first Substructure_Atom_Specifier for each atom.
*/

int
Single_Substructure_Query::create_molecule (Molecule & m, int fill_min_ncon, 
                                     int set_implicit_hydrogens)
{
  assert (0 == m.natoms());
  assert (_root_atoms.number_elements());

  assign_unique_numbers();

  if (_min_atoms_in_query <= 0)
    (void) min_atoms_in_query();

//write ("b4createmol");
//cerr << "Seems the query contains a min of " << _min_atoms_in_query << " atoms\n";

  if (_comment.length())
    m.set_name(_comment);

  m.resize(_min_atoms_in_query);

  return _root_atoms[0]->create_molecule(m, fill_min_ncon, set_implicit_hydrogens);
}

/*int
Single_Substructure_Query::_initialise_sub_array (Molecule & m,
                                                  Substructure_Query & substitutions_only_at,
                                                  Molecule_to_Query_Specifications & mqs)
{
  Substructure_Results sresults;
  int nhits = substitutions_only_at.substructure_search(m, sresults);

  if (0 == nhits)
  {
    cerr << "Single_Substructure_Query::_initialise_sub_array:no hits to substitutions_only_at query\n";
    return 0;
  }
  
  int matoms = m.natoms();

  int * sub = new int[matoms]; std::unique_ptr<int[]> free_sub(sub);
  assert (NULL != sub);

  m.ncon(sub);     // sub[i] = m.ncon(i)

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);
    
    e->set_vector(sub, 0);
  }

  cerr << "Warning, removing isotopes\n";

  m.transform_to_non_isotopic_form();

  mqs.delete_substitution();

  return mqs.set_substitution(sub, matoms);
}*/

int
Single_Substructure_Query::create_from_molecule (MDL_Molecule & m, 
                                                 Molecule_to_Query_Specifications & mqs,
                                                 const int * include_these_atoms)
{
  assert (m.ok());
  assert (0 == _root_atoms.number_elements());

  if (! m.arrays_allocated())
    m.build(m);

#ifdef DEBUG_CREATE_FROM_MOLECULE
  cerr << "Single_Substructure_Query::create_from_molecule line " << __LINE__ << " " << m.name() << endl;
  for (int i = 0; i < m.natoms(); ++i)
  {
    if (include_these_atoms[i])
      cerr << "Will include atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << endl;
  }
#endif

//cerr << "Single_Substructure_Query::create_from_molecule:BEGIN   " << m.smiles() << endl;

  if (! _create_from_molecule(m, mqs, include_these_atoms))
    return 0;

  assign_unique_numbers();

#ifdef FINGERPRINT_SUBSTRUCTURE_SEARCHES
  if (use_fingerprints_for_screening_substructure_searches())
  {
    assert (NULL == _fingerprint);
    IWMFingerprint iwmfingerprint;
    iwmfingerprint.construct_fingerprint(m);

    _fingerprint = new IW_Bits_Base;
    (*_fingerprint) = iwmfingerprint;
  }
#endif

  if (mqs.set_element_hits_needed_during_molecule_to_query())
    _build_element_hits_needed(m, mqs, include_these_atoms);

  if (mqs.ncon() > 0)
    _ncon.add(mqs.ncon());
  if (mqs.min_ncon() > 0)
    _ncon.set_min(mqs.min_ncon());
  if (mqs.max_ncon() > 0)
    _ncon.set_max(mqs.max_ncon());

  return 1;
}

/*
  This is complicated by the fact that a query with A atoms in it
  may have Carbons put in there for aromaticity determinations and
  such...
*/

int
Single_Substructure_Query::_build_element_hits_needed (const MDL_Molecule & m,
                                const Molecule_to_Query_Specifications & mqs,
                                const int * include_these_atoms)
{
  int * count = new_int(HIGHEST_ATOMIC_NUMBER + 1); std::unique_ptr<int[]> free_count(count);

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (NULL != include_these_atoms && 0 == include_these_atoms[i])
      continue;

    if (only_include_isotopically_labeled_atoms() && 0 == m.isotope(i))
      continue;

    const Element * e = m.elementi(i);

    if (! e->is_in_periodic_table())
      continue;

    if (1 == e->atomic_number() && (mqs.convert_explicit_hydrogens_to_match_any_atom() || mqs.convert_explicit_hydrogens_to_match_any_atom_including_hydrogen()))
      continue;

    const MDL_Atom_Data * mad = m.mdl_atom(i);

    if ("A" == mad->initial_atomic_symbol())
      continue;
    
    if ("Q" == mad->initial_atomic_symbol())
      continue;

    if (mad->alias().length())
      continue;

    const ISIS_Atom_List & alist = mad->atom_list();

    if (alist.active())
      continue;

    atomic_number_t z = m.atomic_number(i);

    count[z]++;
  }

  for (int i = 0; i <= HIGHEST_ATOMIC_NUMBER; i++)
  {
    if (0 == count[i])
      continue;

    Elements_Needed * tmp = new Elements_Needed(i);
    tmp->set_min(count[i]);

    _elements_needed.add(tmp);
  }

  return 1;
}

static atom_number_t
first_isotopically_labelled_atom_in_fragment (Molecule & m,
                                              int f)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (f != m.fragment_membership(i))
      continue;

    if (0 != m.isotope(i))
      return i;
  }

  return INVALID_ATOM_NUMBER;
}

static atom_number_t
first_allowed_atom_in_fragment (Molecule & m,
                                const int f,
                                const extending_resizable_array<Substructure_Atom *> & atoms_already_done,
                                const int * include_these_atoms)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (NULL == include_these_atoms)
      ;
    else if (0 == include_these_atoms[i])
      continue;

    if (i >= atoms_already_done.number_elements())
      ;
    else if (NULL != atoms_already_done[i])
      continue;

    if (m.fragment_membership(i) != f)
      continue;

    return i;
  }

  return INVALID_ATOM_NUMBER;
}

int
Single_Substructure_Query::_create_from_molecule (MDL_Molecule & m,
                                                  Molecule_to_Query_Specifications & mqs,
                                                  const int * include_these_atoms)
{
  int matoms = m.natoms();
 
  if (0 == matoms)
    return 0;

  _respect_initial_atom_numbering = 1;

  set_comment(m.molecule_name());

  if (mqs.condense_explicit_hydrogens_to_anchor_atoms())
  {
    if (NULL != include_these_atoms)
    {
      cerr << "Single_Substructure_Query::_create_from_molecule:include_these_atoms not handled, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
      exit(1);
    }
    m.remove_explicit_hydrogens();
    m.compute_aromaticity_handle_atom_lists();
    matoms = m.natoms();
  }

  m.compute_aromaticity_handle_atom_lists();
  
#ifdef DEBUG_INITIALISE_SUB_ARRAY_STUFF
  cerr << "_create_from_molecule:before initialise_sub_array_from_only_sub_query\n";
  cerr << "hcount " << mqs.hcount() << endl;
  cerr << "h0designator " << mqs.h0designator() << endl;
  cerr << "min_hcount " << mqs.min_hcount() << endl;
#endif

  if (mqs.built_from_isis_reaction_file())
    ;
  else if (substituents_only_at_isotopic_atoms())
    m.only_allow_substitutions_at_isotopic_atoms(mqs);
  else if (mqs.substitutions_only_at().active())
    m.determine_attachment_points_by_query(mqs);
  else if (substitutions_only_at_non_isotopic_atoms())
    m.only_allow_substitutions_at_non_isotopic_atoms();

#ifdef DEBUG_INITIALISE_SUB_ARRAY_STUFF
  cerr << "_create_from_molecule:after initialise_sub_array_from_only_sub_query\n";
  cerr << "hcount " << mqs.hcount() << endl;
  cerr << "h0designator " << mqs.h0designator() << endl;
  cerr << "min_hcount " << mqs.min_hcount() << endl;
#endif

  if (mqs.min_extra_atoms_in_target() >= 0)
    _natoms.set_min(matoms + mqs.min_extra_atoms_in_target());
  else
    _natoms.set_min(matoms);

  if (mqs.max_extra_atoms_in_target() >= 0)
    _natoms.set_max(matoms + mqs.max_extra_atoms_in_target());

  if (mqs.min_fraction_atoms_matched() > 0.0F)
    _min_fraction_atoms_matched = mqs.min_fraction_atoms_matched();

  if (mqs.max_fraction_atoms_matched() > 0.0F)
    _max_fraction_atoms_matched = mqs.max_fraction_atoms_matched();

  if (_min_fraction_atoms_matched > _max_fraction_atoms_matched)
  {
    cerr << "Single_Substructure_Query::_create_from_molecule:inconsistent min " << _min_fraction_atoms_matched << " and max " << _max_fraction_atoms_matched << " fraction matched\n";
    return 0;
  }

  if (mqs.use_preference_values_to_distinguish_symmetry())
    _sort_by_preference_value = 1;

  extending_resizable_array<Substructure_Atom *> tmp;

  int nf = m.number_fragments();

  if (only_include_isotopically_labeled_atoms())
  {
    for (int i = 0; i < nf; i++)
    {
      atom_number_t astart = first_isotopically_labelled_atom_in_fragment(m, i);

      if (INVALID_ATOM_NUMBER == astart)
        continue;

      Substructure_Atom * r = new Substructure_Atom;

      if (! r->create_from_molecule(m, m, astart, NULL, mqs, tmp))
      {
        cerr << "Single_Substructure_Query::_create_from_molecule:failure for fragment " << i << endl;
        return 0;
      }
    
      _root_atoms.add(r);
    }
  }
  else if (NULL != include_these_atoms)
  {
//#define ECHO_SUBSET
#ifdef ECHO_SUBSET
    cerr << "_create_from_molecule:have subset of atoms\n";
    for (int i = 0; i < matoms; ++i)
    {
      cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " inc " << include_these_atoms[i] << endl;
    }
#endif
    for (int i = 0; i < nf; ++i)
    {
      atom_number_t astart;
      while (INVALID_ATOM_NUMBER != (astart = first_allowed_atom_in_fragment(m, i, tmp, include_these_atoms)))
      {
        Substructure_Atom * r = new Substructure_Atom;

        if (! r->create_from_molecule(m, m, astart, NULL, mqs, tmp, include_these_atoms))
        {
          cerr << "Single_Substructure_Query::_create_from_molecule:failure to build from atom " << astart << " in fragment " << i << endl;
          return 0;
        }

        _root_atoms.add(r);
      }
    }
  }
  else
  {
    for (int i = 0; i < nf; i++)
    {
      Substructure_Atom * r = new Substructure_Atom;

      atom_number_t astart = m.first_atom_in_fragment(i);
      
      if (! r->create_from_molecule(m, m, astart, NULL, mqs, tmp))
      {
        cerr << "Single_Substructure_Query::_create_from_molecule:failure for fragment " << i << endl;
        return 0;
      }
    
      _root_atoms.add(r);
    }
  }

  if (0 == _root_atoms.number_elements())
  {
    cerr << "Single_Substructure_Query::create_from_molecule:no root atoms found\n";
    return 0;
  }

//cerr << "Created query with " << _root_atoms.number_elements() << " root atoms, atoms in query " << _max_atoms_in_query << endl;

  int nr = m.number_sssr_rings();

  if (only_include_isotopically_labeled_atoms())
  {
    // some time put in a computation of the number of rings with isotopes - is complicated by fused rings...
  }
  else if (NULL != include_these_atoms)   // again, kind of hard
    ;
  else if (nr)
    _nrings.set_min(nr);

  for (int i = 0; i < m.number_link_atoms(); i++)
  {
    const Link_Atom * l = m.link_atom(i);;

    if (! add_link_atom(*l))
    {
      cerr << "ISIS Link Atom " << (*l) << endl;
      return 0;
    }
  }

  for (int i = 0; i < m.chiral_centres(); i++)
  {
    const Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);
    _add_chiral_centre(m, *c);
  }

  if (mqs.environment_near_substitution_points_specified() || mqs.environment_no_match_near_substitution_points_specified())
  {
    if (! _add_environment_according_to_matched_atoms(mqs))
    {
      cerr << "Single_Substructure_Query::_create_from_molecule:cannot build environment\n";
      return 0;
    }
  }

  if (NULL != include_these_atoms)
  {
    int * xref = new_int(matoms); std::unique_ptr<int[]> free_xref(xref);
    int ndx = 0;
    for (int i = 0; i < matoms; ++i)
    {
      if (include_these_atoms[i])
      {
        xref[i] = ndx;
        ndx++;
      }
    }

    _natoms.set_min(ndx);

    for (int i = 0; i < _root_atoms.number_elements(); ++i)
    {
      _root_atoms[i]->adjust_initial_atom_numbers(xref);
    }
  }

  return 1;
}

int
Substructure_Query::_create_query_and_add (MDL_Molecule & m,
                                           Molecule_to_Query_Specifications & mqs,
                                           const int * include_these_atoms)
{
//cerr << "Substructure_Query::_create_query_and_add:building from '" << m.smiles() << "'\n";

  Single_Substructure_Query * q = new Single_Substructure_Query;

  if (! q->create_from_molecule(m, mqs, include_these_atoms))
  {
    delete q;

    return 0;
  }

  add(q);

  _operator.add_operator('|');

  return 1;
}

/*
  The default behaviour for an ISIS atom list is to enumerate the possibilities
  as separate query components
*/

int
Substructure_Query::create_from_molecule (MDL_Molecule & m, 
                                          Molecule_to_Query_Specifications & mqs,
                                          const int * include_these_atoms)
{
  _comment = m.name();

#ifdef DEBUG_CREATE_FROM_MOLECULE
  cerr << "Substructure_Query::create_from_molecule:start with '" << m.smiles() << "'\n";
#endif

  const int nlink = m.number_link_atoms();

  if (0 == nlink)
    return _create_query_and_add(m, mqs, include_these_atoms);

  Link_Atom_Current_State * lacc = new Link_Atom_Current_State[nlink];
  std::unique_ptr<Link_Atom_Current_State[]> free_lacc(lacc);

  for (int i = 0; i < nlink; i++)
  {
    const Link_Atom * l = m.link_atom(i);

    lacc[i].initialise(l);
  }

// Transfer the link atoms out of the molecule and into this array

  resizable_array_p<Link_Atom> link_atoms;
  m.transfer_link_atoms(link_atoms);

  assert (nlink == link_atoms.number_elements());

  int rc = _enumerate_list_atom_possibilities(m, mqs, link_atoms, lacc, nlink, 0, include_these_atoms);

#ifdef DEBUG_CREATE_FROM_MOLECULE
  cerr << "Substructure_Query::create_from_molecule:created " << _number_elements << " molecules\n";
  cerr << "Final molecule is '" << m.smiles() << "'\n";
#endif

  return rc;
}

int
Substructure_Query::_enumerate_list_atom_possibilities (const MDL_Molecule & m,
                                                        Molecule_to_Query_Specifications & mqs,
                                                        const resizable_array_p<Link_Atom> & link_atoms,
                                                        Link_Atom_Current_State * lacc,
                                                        int nlink,
                                                        int ndx,
                                                        const int * include_these_atoms)
{
  lacc[ndx].reset();

  MDL_Molecule tmp(m);

  while (link_atoms[ndx]->create_next_variant(tmp, lacc[ndx]))
  {
    if (ndx == nlink - 1)
      _create_query_and_add(tmp, mqs, include_these_atoms);
    else
      _enumerate_list_atom_possibilities(tmp, mqs, link_atoms, lacc, nlink, ndx + 1, include_these_atoms);
  }

  return 1;
}

int
Substructure_Query::create_from_molecule (const Molecule & m,
                                          Molecule_to_Query_Specifications & mqs,
                                          const int * include_these_atoms)
{
  MDL_Molecule tmp(m);

  tmp.set_name(m.name());

  return create_from_molecule(tmp, mqs, include_these_atoms);

}

static int
add_environment (const msi_object & msi,
                 atom_number_t s,
                 extending_resizable_array<Substructure_Atom *> & completed,
                 resizable_array_p<Substructure_Environment> & destination)
{
  Substructure_Environment * e = new Substructure_Environment;
  if (! e->construct_from_msi_object(msi, completed, s, SINGLE_BOND))
  {
    cerr << "Single_Substructure_Query::_add_environment_according_to_matched_atoms:cannot parse environment specification\n";
    return 0;
  }

  destination.add(e);

  return 1;
}

int
Single_Substructure_Query::_add_environment_according_to_matched_atoms (Molecule_to_Query_Specifications & mqs)
{
  const msi_object & env = mqs.environment_near_substitution_points();
  const msi_object & env_no_match = mqs.environment_no_match_near_substitution_points();
  const Set_of_Atoms & substitution_points = mqs.substitution_points();

  int n = substitution_points.number_elements();

  if (0 == n)
  {
    cerr << "Single_Substructure_Query::_add_environment_according_to_matched_atoms:no substitution points\n";
    return 0;
  }

  extending_resizable_array<Substructure_Atom *> completed(NULL);

  _collect_all_atoms(completed);

#ifdef CHECK_ATOMS_DEFINED
  for (int i = 0; i < completed.number_elements(); i++)
  {
    if (NULL == completed[i])
      cerr << "No atom " << i << endl;
    else
      cerr << "Atom " << i << " is defined\n";
  }
#endif

  for (int i = 0; i < n; i++)
  {
    atom_number_t s = substitution_points[i];

    if (NULL == completed[s])
    {
      cerr << "Single_Substructure_Query::_add_environment_according_to_matched_atoms:no atom " << s << endl;
      return 0;
    }

    if (env.active())
    {
      if (! add_environment(env, s, completed, _environment))
        return 0;
    }

    if (env_no_match.active())
    {
      if (! add_environment(env_no_match, s, completed, _environment_rejections))
        return 0;
    }
  }

  return 1;
}

/*
  Each object in an msi object file needs an id
  The Substructure_Atoms get their own unique_id value. so we need some numbers for
  everyone else. This is variable OBJECT_ID
*/

int
Single_Substructure_Query::write_msi(std::ostream & os, int & object_id,
                                     const const_IWSubstring & logical_operator,
                                     int indentation)
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << " Query\n";

  os << ind << "  (A I " << NAME_OF_VERSION_ATTRIBUTE << " 2)\n";

  if (_comment.length())
    os << ind << "  (A C " << NAME_OF_COMMENT_ATTRIBUTE << " \"" << _comment << "\")\n";

  if (logical_operator.length())
    os << ind << "  (A C " << NAME_OF_OPERATOR_ATTRIBUTE << " \"" << logical_operator << "\")\n";

  for (int i = 0; i < _numeric_value.number_elements(); i++)
  {
    os << ind << "  (A D " NAME_OF_NUMERIC_VALUE_ATTRIBUTE << ' ' << _numeric_value[i] << ")\n";
  }

  if (_rejection)
    os << ind << "  (A I " << NAME_OF_REJECTION_ATTRIBUTE << " 1)\n";

  _attached_heteroatom_count.write_msi(os, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, indentation + 2);
  if (_heteroatoms.number_elements())
  {
    os << ind << "  (A I " << NAME_OF_DEFINE_HETEROATOMS_ATTRIBUTE << ' ';
    if (1 == _heteroatoms.number_elements())
      os << _heteroatoms[0];
    else
    {
      os << '(';
      for (int i = 0; i < _heteroatoms.number_elements(); i++)
      {
        if (i)
          os << ' ';
        os << _heteroatoms[i];
      }
      os << ')';
    }
    os << ")\n";
  }

  _heteroatoms_in_molecule.write_msi(os, NAME_OF_HETEROATOMS_ATTRIBUTE, indentation + 2);

  _heteroatoms_matched.write_msi(os, NAME_OF_HETEROATOMS_MATCHED_ATTRIBUTE, indentation + 2);

  if (_find_one_embedding_per_start_atom)
    os << ind << "  (A I " << NAME_OF_ONE_EMBEDDING_PER_START_ATOM << " 1)\n";
  if (_find_unique_embeddings_only)
    os << ind << "  (A I " << NAME_OF_UNIQUE_EMBEDDINGS_ONLY << " 1)\n";
  if (_embeddings_do_not_overlap)
    os << ind << "  (A I " << NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE << " 1)\n";
  if (_sort_by_preference_value)
    os << ind << "  (A I " << NAME_OF_SORT_BY_PREFERENCE_VALUE_ATTRIBUTE << " 1)\n";

  _hits_needed.write_msi(os, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, indentation + 2);

  if (_normalise_rc_per_hits_needed)
    os << ind << "  (A I " << NAME_OF_NORMALISE_RC_PER_HITS_NEEDED << " 1)\n";

  if (_subtract_from_rc)
    os << ind << "  (A I " << NAME_OF_SUBTRACT_FROM_RC << ' ' << _subtract_from_rc << ")\n";
    
  if (_all_hits_in_same_fragment)
    os << ind << "  (A I " << NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE << ' ' << _all_hits_in_same_fragment << ")\n";

  if (_only_keep_matches_in_largest_fragment)
    os << ind << "  (A I " << NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE << ' ' << _only_keep_matches_in_largest_fragment << ")\n";

  if (_respect_initial_atom_numbering)
    os << ind << "  (A I " << NAME_OF_RESPECT_INITIAL_NUMBERING_ATTRIBUTE << " 1)\n";

  if (_max_matches_to_find)
    os << ind << "  (A I " << NAME_OF_MAX_MATCHES_TO_FIND << ' ' << _max_matches_to_find << ")\n";

  if (0 == _save_matched_atoms)     // the default is 1
    os << ind << "  (A I " << NAME_OF_SAVE_HITS_ATTRIBUTE << " 0)\n";

// since this flag is rare, only write it if set to the non standard value

  if (0 == _environment_must_match_unmatched_atoms)
    os << ind << "  (A I " << NAME_OF_ENVIRONMENT_MUST_MATCH_UNMATCHED_ATOMS_ATTRIBUTE << " 0)\n";

  _natoms.write_msi(os, NAME_OF_NATOMS_ATTRIBUTE, indentation + 2);
  _nrings.write_msi(os, NAME_OF_NRINGS_ATTRIBUTE, indentation + 2);

  _isolated_rings.write_msi(os, NAME_OF_ISOLATED_RING_COUNT_SPECIFIER, indentation + 2);
  _fused_rings.write_msi(os, NAME_OF_FUSED_RING_COUNT_SPECIFIER, indentation + 2);
  _strongly_fused_rings.write_msi(os, NAME_OF_STRONGLY_FUSED_RING_COUNT_SPECIFIER, indentation + 2);

  _isolated_ring_objects.write_msi(os, NAME_OF_ISOLATED_RING_OBJECTS, indentation + 2);

  _aromatic_rings.write_msi(os, NAME_OF_AROMATIC_RING_COUNT_SPECIFIER, indentation + 2);
  _non_aromatic_rings.write_msi(os, NAME_OF_NON_AROMATIC_RING_COUNT_SPECIFIER, indentation + 2);

  _ncon.write_msi(os, NAME_OF_NCON_ATTRIBUTE, indentation + 2);

  _distance_between_hits.write_msi(os, NAME_OF_DISTANCE_BETWEEN_HITS_ATTRIBUTE, indentation + 2);

  if (_matched_atoms_to_check_for_hits_too_close > 0)
     os << ind << "  (A I " << NAME_OF_DISTANCE_BETWEEN_HITS_NCHECK_ATTRIBUTE << ' ' << _matched_atoms_to_check_for_hits_too_close << ")\n";

   if (_fail_if_embeddings_too_close)
     os << ind << "  (A I " << NAME_OF_FAIL_IF_EMBEDDINGS_TOO_CLOSE_ATTRIBUTE << ' ' << _fail_if_embeddings_too_close << ")\n";

  if (_sort_matches_by)
  {
    IWString tmp;
    if (SORT_MATCHES_BY_INCREASING_NCON & _sort_matches_by)
      tmp << "ncon+";
    else if (SORT_MATCHES_BY_DECREASING_NCON & _sort_matches_by)
      tmp << "ncon-";

    if (SORT_MATCHES_BY_INCREASING_MAXD & _sort_matches_by)
      tmp.append_with_spacer("maxd+", ',');
    else
      tmp.append_with_spacer("maxd-", ',');

    os << ind << "  (A C " << NAME_OF_SORT_MATCHES_BY << " \"" << tmp << "\")\n";
  }

  _net_formal_charge.write_msi(os, NAME_OF_NET_FORMAL_CHARGE_ATTRIBUTE, indentation + 2);

  _atoms_in_spinach.write_msi(os, NAME_OF_SPINACH_ATOMS_ATTRIBUTE, indentation + 2);
  _inter_ring_atoms.write_msi(os, NAME_OF_INTER_RING_ATOMS_ATTRIBUTE, indentation + 2);

  _unmatched_atoms.write_msi(os, NAME_OF_UNMATCHED_ATOMS, indentation + 2);

  if (static_cast<float>(0.0) != _min_fraction_atoms_matched)
    os << ind << "  (A D " << NAME_OF_MIN_FRACTION_ATOMS_MATCHED << ' ' << _min_fraction_atoms_matched << ")\n";

  if (static_cast<float>(0.0) != _max_fraction_atoms_matched)
    os << ind << "  (A D " << NAME_OF_MAX_FRACTION_ATOMS_MATCHED << ' ' << _max_fraction_atoms_matched << ")\n";


  _ring_atoms_matched.write_msi(os, NAME_OF_RING_ATOMS_MATCHED_ATTRIBUTE, indentation + 2);
  _number_isotopic_atoms.write_msi(os, NAME_OF_NUMBER_ISOTOPIC_ATOMS_ATTRIBUTE, indentation + 2);
  _number_fragments.write_msi(os, NAME_OF_NUMBER_FRAGMENTS_ATTRIBUTE, indentation + 2);

  if (_ncon_ignore_singly_connected)
    os << ind << "  (A I " << NAME_OF_NCON_IGNORE_SINGLY_CONNECTED << ' ' <<
                              _ncon_ignore_singly_connected << ")\n";

  if (_do_not_perceive_symmetry_equivalent_matches)
    os << ind << "  (A I " << NAME_OF_DO_NOT_PERCEIVE_SYMMETRIC_EQUIVALENTS << " 1)\n";

  if (0 == _implicit_ring_condition)
    os << ind << "  (A I " << NAME_OF_IMPLICIT_RING_CONDITION << " 0)\n";
  else if (1 == _implicit_ring_condition)
    os << ind << "  (A I " << NAME_OF_IMPLICIT_RING_CONDITION << " 1)\n";

  _distance_between_root_atoms.write_msi(os, NAME_OF_DISTANCE_BETWEEN_ROOT_ATOMS_SPECIFIER, indentation + 2);

  for (int i = 0; i < _elements_needed.number_elements(); i++)
  {
    _elements_needed[i]->write_msi(os, NAME_OF_ELEMENTS_NEEDED_OBJECT, object_id, indentation + 2);
  }

  for (int i = 0; i < _element_hits_needed.number_elements(); i++)
  {
    _element_hits_needed[i]->write_msi(os, NAME_OF_ELEMENT_HITS_NEEDED_OBJECT, object_id, indentation + 2);
  }

  int nr = _ring_specification.number_elements();
  for (int i = 0; i < nr; i++)
  {
    const Substructure_Ring_Specification * r = _ring_specification[i];
    r->write_msi(os, object_id, indentation + 2);
  }

  int nrs = _ring_system_specification.number_elements();
  for (int i = 0; i < nrs; i++)
  {
    const Substructure_Ring_System_Specification * r = _ring_system_specification[i];
    r->write_msi(os, object_id, indentation + 2);
  }

  if (_no_matched_atoms_between.number_elements())
  {
    for (int i = 0; i < _no_matched_atoms_between.number_elements(); i++)
    {
      const Bond * b = _no_matched_atoms_between[i];
      os << ind << "  (A I " << NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER << " (" << b->a1() << ' ' << b->a2() << "))\n";
    }
  }

  for (int i = 0; i < _link_atom.number_elements(); i++)
  {
    _link_atom[i]->write_msi(os, ind, NAME_OF_LINK_ATOMS_SPECIFIER);
  }

  for (int i = 0; i < _chirality.number_elements(); i++)
  {
    _chirality[i]->write_msi(os, ind, NAME_OF_CHIRALITY_SPECIFIER);
  }

  int rc = 1;
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    if (! _root_atoms[i]->write_msi(os, object_id, NAME_OF_QUERY_ATOM_OBJECT, indentation + 2))
      return 0;
  }

  int ne = _environment.number_elements();
  for (int i = 0; i < ne; i++)
  {
    Substructure_Environment * e = _environment[i];
    e->write_msi(os, object_id, indentation + 2);
  }

  nr = _environment_rejections.number_elements();
  for (int i = 0; i < nr; i++)
  {
    Substructure_Environment * r = _environment_rejections[i];
    r->write_msi(os, object_id, indentation + 2);
  }

  if (os.good())
    os << ind << ")\n";

  return rc;
}

int
Single_Substructure_Query::write_msi (std::ostream & os)
{
  if (0 == _max_atoms_in_query)
    (void) max_atoms_in_query();

  int i = _max_atoms_in_query + 1;
  return write_msi(os, i, "", 0);
}

int
Single_Substructure_Query::write (const char * fname)
{
  assert (NULL != fname);

  std::ofstream os(fname, std::ios::out);
  if (! os.good())
  {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_msi(os);
}

/*
  In constructing from an msi object, we have found a Substructure_Ring_Specification
  object. Parse it.
*/

int
Single_Substructure_Query::_parse_ring_specifier_object (const msi_object & msi)
{
  assert (NAME_OF_RING_SPECIFIER_OBJECT == msi.name());

  Substructure_Ring_Specification * r = new Substructure_Ring_Specification();
  if (! r->construct_from_msi_object(msi))
  {
    delete r;

    cerr << "Could not create ring specifier object from msi object\n";
    cerr << msi;
    return 0;
  }

  if (0 == _ring_specification.number_elements())    // no operators with first object
    ;
  else if (! extract_operator(msi, _ring_specification_logexp, IW_LOGEXP_AND, "Single_Substructure_Query::_parse_ring_specifier_object"))
    return 0;

  _ring_specification.add(r);

  return 1;
}

/*
  In constructing from an msi object, we have found a Substructure_Ring_System_Specification
  object. Parse it.
*/

int
Single_Substructure_Query::_parse_ring_system_specifier_object(const msi_object & msi)
{
  assert (NAME_OF_RING_SYSTEM_SPECIFIER_OBJECT == msi.name());

  Substructure_Ring_System_Specification * r = new Substructure_Ring_System_Specification();
  if (! r->construct_from_msi_object(msi))
  {
    delete r;

    cerr << "Could not create ring system specifier object from msi object\n";
    cerr << msi;
    return 0;
  }

  if (0 == _ring_system_specification.number_elements())
    ;
  else if (! extract_operator(msi, _ring_system_specification_logexp, IW_LOGEXP_AND, "Single_Substructure_Query::_parse_ring_system_specifier_object"))
  {
    delete r;
    return 0;
  }

  _ring_system_specification.add(r);

  return 1;
}

/*
  Is an msi object a root atom or not.  Basically, if it has an
  attribute which ends in "bond" (so as to not also match "nbonds"),
  then it is not a root
*/

static int
is_root_substructure_atom (const msi_object & msi)
{
  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++)
  {
    const msi_attribute * att = msi.attribute(i);
    if (att->name().ends_with("bond"))    // got a bond, not a root
      return 0;
    if (NAME_OF_BOND_SMARTS_ATTRIBUTE == att->name())   // a smarts bond, not a root
      return 0;
  }

  int nmsi = msi.number_elements();
  for (int i = 0; i < nmsi; i++)
  {
    const msi_object * m = msi.item(i);
    if (NAME_OF_QUERY_BOND_OBJECT == m->name())
      return 0;
  }

  return 1;
}

/*
  The environment object 
  It can specify query_atoms, smiles or smarts
*/

int
Single_Substructure_Query::_construct_environment_from_msi_object (const msi_object & msi,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               resizable_array_p<Substructure_Environment> & env)
{
//cerr << "Constructing env from\n" << (*msi);

  Substructure_Environment * a = new Substructure_Environment();

  _collect_all_atoms(completed);

  int rc = a->construct_from_msi_object(msi, completed);

  if (0 == rc)
  {
    cerr << msi;
    delete a;
    return 0;
  }

  env.add(a);

  return 1;
}

int
Substructure_Environment::_add_possible_parent (atom_number_t possible_parent,
                                                bond_type_t possible_parent_bond_type,
                                                extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 != possible_parent_bond_type)     // only add if bond_type has defined a particular type
    _bond.set_type(possible_parent_bond_type);
  else
    _bond.set_match_any();

  _possible_parents.add(completed[possible_parent]);

  return 1;
}

int
Substructure_Environment::construct_from_msi_object (const msi_object & msi,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               atom_number_t possible_parent,
                               bond_type_t possible_parent_bond_type)
{
  int nmsi = msi.number_elements();

  if (NAME_OF_ENVIROMENT_REJECTION_OBJECT == msi.name())
    _query_environment_match_as_rejection = 1;
  else if (NAME_OF_ENVIRONMENT_OBJECT == msi.name())
  {
    _query_environment_match_as_rejection = 0;

    const msi_attribute * att = msi.attribute(NAME_OF_ENVIRONMENT_REJECTION);
    if (att)
    {
      int i;
      if (! att->value(i))
      {
        cerr << "The rejection attribute must be numeric\n";
        return 0;
      }
      _query_environment_match_as_rejection = i;
    }
  }

// Process any bonds which have been specified as attributes

  if (! _process_attribute_bonds(msi, completed))
  {
    msi.print(cerr);
    return 0;
  }

  if (INVALID_ATOM_NUMBER != possible_parent)
    _add_possible_parent(possible_parent, possible_parent_bond_type, completed);

  int np = _possible_parents.number_elements();

  if (0 == np)
  {
    cerr << "Substructure_Environment::construct_from_msi_object: environment not connected\n";
    return 0;
  }

  int notused = 0;
  if (! really_gruesome(_hits_needed, msi, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, notused, 0, MAX_NOT_SPECIFIED))
    return 0;

  const msi_attribute * att = msi.attribute(NAME_OF_NO_OTHER_SUBSTITUENTS_ALLOWED_ATTRIBUTE);
  if (att)
  {
//  cerr << "Got " << NAME_OF_NO_OTHER_SUBSTITUENTS_ALLOWED_ATTRIBUTE << endl;
    if (! att->value(_no_other_substituents_allowed))
    {
      cerr << "Substructure_Environment::construct_from_msi_object: the " << NAME_OF_NO_OTHER_SUBSTITUENTS_ALLOWED_ATTRIBUTE <<
              " must be a whole number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_HYDROGEN_OK_AS_ENVIRONMENT_MATCH);
  if (att)
  {
    if (! att->value(_hydrogen_ok_as_environment_match) || _hydrogen_ok_as_environment_match < 0)
    {
      cerr << "Substructure_Environment::construct_from_msi_object: the " << NAME_OF_HYDROGEN_OK_AS_ENVIRONMENT_MATCH << " must be a whole number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_MAX_ENVIRONMENT_MATCHES_PER_ATTACHMENT_POINT);
  if (att)
  {
    if (! att->value(_max_environment_matches_per_attachment_point) || _max_environment_matches_per_attachment_point <= 0)
    {
      cerr << "Substructure_Environment::construct_from_msi_object: then " << NAME_OF_MAX_ENVIRONMENT_MATCHES_PER_ATTACHMENT_POINT << " must be a +ve whole number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS);
  if (att)
  {
    if (! att->value(_environments_can_share_attachment_points) || _environments_can_share_attachment_points < 0)
    {
      cerr << "Substructure_Environment::construct_from_msi_object: the " << NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS << " attribute must be a whole +ve number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_MAX_MATCHES_TO_FIND);
  if (att)
  {
    if (! att->value(_max_matches_to_find) || _max_matches_to_find < 1)
    {
      cerr << "Substructure_Environment::construct_from_msi_object:the " << NAME_OF_MAX_MATCHES_TO_FIND << " attribute must be a whole +ve number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_OR_OPERATOR);
  if (att)
  {
    if (! att->value(_or_id) || _or_id <= 0)
    {
      cerr << "Substructure_Environment::construct_from_msi_object: the " << NAME_OF_OR_OPERATOR <<
              " must be a positive whole number\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_AND_OPERATOR);
  if (att)
  {
    if (! att->value(_and_id) || _and_id <= 0)
    {
      cerr << "Substructure_Environment::construct_from_msi_object: the " << NAME_OF_AND_OPERATOR  <<
              " must be a positive whole number\n";
      return 0;
    }
  }

  if (_and_id && _or_id)
  {
    cerr << "Substructure_Environment::construct_from_msi_object:cannot have both OR and AND specifications\n";
    return 0;
  }

// Now construct the structural specification. There can be any number

  int i = 0;
  while (NULL != (att = msi.attribute(NAME_OF_SMARTS_ATTRIBUTE, i++)))
  {
    Substructure_Atom * a = new Substructure_Atom;
    a->set_unique_id(msi.object_id());

    const_IWSubstring x;
    att->value(x);
    const char * s = x.data();

    if (x.length() && (isdigit(s[0]) || '>' == s[0] || '<' == s[0]))
    {
      int value, qualifier;
      int chars_consumed = smarts_fetch_numeric(s, value, qualifier);
      if (0 == chars_consumed)
      {
        cerr << "Substructure_Environment::construct_from_msi_object:invalid numeric qualifier '" << x << "'\n";
        delete a;
        return 0;
      }
      if (qualifier > 0)
        _hits_needed.set_min(value+1);
      else if (qualifier < 0)
        _hits_needed.set_max(value-1);
      else
        _hits_needed.add(value);

      x += chars_consumed;
      if (! a->parse_smarts_specifier(x))
      {
        cerr << "Substructure_Environment::construct_from_msi_object:invalid smarts '" << x << "'\n";
        return 0;
      }
    }

    else if (! a->parse_smarts_specifier(att))
    {
      delete a;
      return 0;
    }

//  cerr << "Build query from smarts " << (*att) << endl;
    add(a);
  }

  i = 0;
  while (NULL != (att = msi.attribute(NAME_OF_SMILES_ATTRIBUTE, i++)))
  {
    Substructure_Atom * a = new Substructure_Atom;
    a->set_unique_id(msi.object_id());

    if (! a->parse_smiles_specifier(att))
    {
      delete a;
      return 0;
    }

//  cerr << "Build query from smiles " << (*att) << endl;
    add(a);
  }

  for (i = 0; i < nmsi; i++)
  {
    msi_object & m = *(msi[i]);

    if (NAME_OF_QUERY_BOND_OBJECT == m.name())
    {
    }
    else if (NAME_OF_QUERY_ATOM_OBJECT == m.name())
    {
      Substructure_Atom * a = new Substructure_Atom;
      if (! a->construct_from_msi_object(m, completed))
      {
        delete a;
        return 0;
      }

      if (is_root_substructure_atom(m))
        add(a);
    }
    else
    {
      cerr << "Unknown msi object type in query_environment\n" << m;
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_NEXT_CONNECTED

/*
  Parse disconnected smarts.
  Valid inputs include
  'CC'
  'C.O'
*/

static int
next_disconnected (const const_IWSubstring & smarts,
                   int & istart,
                   const_IWSubstring & disconnected)
{
  if (istart >= smarts.length())     // we finished on our previous call
    return 0;

#ifdef DEBUG_NEXT_CONNECTED
  cerr << "On entry, smarts is '" << smarts << "' istart = " << istart << endl;
#endif

  disconnected = smarts;

  if (istart)
    disconnected += istart;

  if (! disconnected.starts_with('.'))
    ;
  else if (disconnected.starts_with(".."))
  {
    cerr << "Invalid smarts no matched atoms specifier .... '" << smarts << "'\n";
    abort();
  }
  else               // starts with a period
  {
    disconnected++;
    istart++;
  }

  int nchars = disconnected.length();

  for (int i = 0; i < nchars; i++)
  {
    if ('.' != disconnected[i])
      continue;

    if (i == nchars - 1)    // gack, smarts ends in ., will be found illegal later
    {
      istart = smarts.length() + 9;
      return 1;
    }

    if ('.' != disconnected[i + 1])
    {
      disconnected.iwtruncate(i);
      istart += i;
      return 1;
    }

//  We have a period. The only valid thing is 3 consecutive periods

    i += 2;
  }

// If we come to here, there are no more disconnected pieces, we are done

  istart = smarts.nchars() + 9;     // we are done
  return 1;
}

template <typename T>
void
transfer_to_our_array (resizable_array_p<T> & to,
                       const resizable_array<T *> & from)
{
  for (int i = 0; i < from.number_elements(); i++)
  {
    to.add(from[i]);
  }

  return;
}

template void transfer_to_our_array (resizable_array_p<Substructure_Atom> &, const resizable_array<Substructure_Atom *> &);
template void transfer_to_our_array (resizable_array_p<Bond> &, const resizable_array<Bond *> &);
template void transfer_to_our_array (resizable_array_p<Link_Atom> &, const resizable_array<Link_Atom *> &);

//#define DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER

/*
  Parsing a smarts is complicated by the possible presence of the '.' character,
  and the need to process the '...' construct (no_matched_atoms_between)

  Valid inputs are

  'C-O-C'
  'C-O-C.[NH2]'
  'Br...C(=O)Cl'

*/

int
Single_Substructure_Query::_parse_smarts_specifier (const const_IWSubstring & smarts)
{
#ifdef DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER
  cerr << "Parsing smarts '" << smarts << "'\n";
#endif

  if (smarts.starts_with('.'))
  {
    cerr << "Single_Substructure_Query::_parse_smarts_specifier: smarts cannot start with '.'\n";
    cerr << "Bad smarts '" << smarts << "'\n";
    return 0;
  }

  const_IWSubstring disconnected;
  int i = 0;

// QUERY_ATOMS_CREATED helps keep track of atom numbers across disconnected smarts

  int query_atoms_created = 0;

  int component_grouping = 0;
  int component_id = 0;

  Parse_Smarts_Tmp pst;

  pst.set_natoms(smarts.length());

  while (next_disconnected(smarts, i, disconnected))
  {
#ifdef DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER
    cerr << "Smarts disconnected token '" << disconnected << "' being processed\n";
#endif

    int paren_balance = disconnected.balance('(', ')');

    if (disconnected.starts_with('('))
    {
      if (component_grouping)
      {
        cerr << "Single_Substructure_Query::_parse_smarts_specifier: component grouping cannot be nested\n";
        return 0;
      }

      disconnected++;
      component_id++;
      component_grouping = 1;

      paren_balance--;
      _fragment_ids_present = 1;
    }

    Substructure_Atom * a = new Substructure_Atom;

    if (component_grouping)
      a->set_fragment_id(component_id);

    if (-1 == paren_balance && disconnected.ends_with(')'))       // closing a component grouping
    {
      disconnected.chop();
      component_grouping = 0;
    }

    pst.add_root_atom(a);

    if (! a->parse_smarts_specifier(disconnected, pst, query_atoms_created))
    {
      cerr << "Single_Substructure_Query::_parse_smarts_specifier: cannot parse '" << disconnected << "'\n";
      delete a;
      return 0;
    }

    query_atoms_created = pst.last_query_atom_created() + 1;

//  cerr << "Query atoms created incremented to " << query_atoms_created << endl;
  }

  transfer_to_our_array(_root_atoms, pst.root_atoms());

  transfer_to_our_array(_no_matched_atoms_between, pst.no_matched_atoms_between());

  transfer_to_our_array(_link_atom, pst.link_atoms());

  if (! ignore_chirality_in_smarts_input())
    _build_chirality_specifications_from_atom_chiral_info();

#ifdef DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER
  cerr << "After parsing '" << smarts << "' we have " << _root_atoms.number_elements() << " root atoms\n";
#endif

  return 1;
}

static int
fetch_heteroatom_definitions(const msi_attribute * att,
                             resizable_array<atomic_number_t> & heteroatoms)
{
  int ndx = 0;
  if (att->valid_as_int())
  {
    int z;
    while (att->next_value(z, ndx))
    {
      if (z < 0 || z > HIGHEST_ATOMIC_NUMBER)
      {
        cerr << "Invalid atomic number in heteroatom spec '" << z << "'\n";
        return 0;
      }

      heteroatoms.add(z);
    }
  }
  else
  {
    IWString s = att->stringval();

    int i = 0;
    const_IWSubstring token;
    while (s.nextword(token, i))
    {
      const Element * e = get_element_from_symbol_no_case_conversion(token);
      if (NULL == e)
      {
        cerr << "fetch_heteroatom_definitions:unrecognised element specification '" << s << "'\n";
        return 0;
      }

      if (e->atomic_number () < 0)
      {
        cerr << "fetch_heteroatom_definitions::element " << s << " has no atomic number\n";
        return 0;
      }

      heteroatoms.add(e->atomic_number());
    }
  }

  return heteroatoms.number_elements();
}

static int
fetch_heteroatom_definitions(const msi_attribute & att,
                             int * heteroatm)
{
  set_vector(heteroatm, HIGHEST_ATOMIC_NUMBER + 1, 0);

  int z;
  int ndx = 0;

  while (att.next_value(z, ndx))
  {
    if (z < 0 || z > HIGHEST_ATOMIC_NUMBER)
    {
      cerr << "Invalid atomic number '" << z << "'\n";
      return 0;
    }

    heteroatm[z] = 1;
  }

  return 1;
}

int
Single_Substructure_Query::_construct_from_msi_object (const msi_object & msi,
                                                int * attribute_recognised)
{
  assert (ok());
  assert (0 == _root_atoms.number_elements());

  int version_attribute_found = 0;

  _numeric_value.resize(0);

  int environments_can_share_attachment_points = -1;    // indicates not set

  int natt = msi.attribute_count();
  for (int i = 0; i < natt; i++)
  {
    const msi_attribute * att = msi.attribute(i);
    if (NAME_OF_COMMENT_ATTRIBUTE == att->name())
    {
      att->value(_comment);
      attribute_recognised[i] = 1;
    }
    else if ("Label" == att->name())
    {
      att->value(_comment);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_VERSION_ATTRIBUTE == att->name())
    {
      int v;
      if (! att->value(v) || 2 != v)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: can only read version 2 queries\n";
        return 0;
      }
      version_attribute_found = 1;
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_DEFINE_HETEROATOMS_ATTRIBUTE == att->name())
    {
      if (! fetch_heteroatom_definitions(att, _heteroatoms))
        return 0;
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_ONE_EMBEDDING_PER_START_ATOM == att->name())
    {
      int oe;
      if (! att->value(oe))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_ONE_EMBEDDING_PER_START_ATOM <<
                "' attribute requires a whole number\n";
        return 0;
      }

      set_find_one_embedding_per_atom(oe);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_UNIQUE_EMBEDDINGS_ONLY == att->name())
    {
      int ue;
      if (! att->value(ue))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_UNIQUE_EMBEDDINGS_ONLY <<
                "' attribute requires a whole number\n";
        return 0;
      }

      set_find_unique_embeddings_only(ue);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_NORMALISE_RC_PER_HITS_NEEDED == att->name())
    {
      int nrc;
      if (! att->value(nrc) || nrc < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_NORMALISE_RC_PER_HITS_NEEDED <<
                "' attribute requires a non negative whole number\n";
        return 0;
      }

      set_normalise_rc_per_hits_needed(nrc);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_SUBTRACT_FROM_RC == att->name())
    {
      int nsr;
      if (! att->value(nsr))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_SUBTRACT_FROM_RC <<
                "' attribute requires a whole number\n";
        return 0;
      }

      _subtract_from_rc = nsr;
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_RESPECT_INITIAL_NUMBERING_ATTRIBUTE == att->name())
    {
      if (! att->value(_respect_initial_atom_numbering))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_RESPECT_INITIAL_NUMBERING_ATTRIBUTE <<
                "' attribute requires a whole number\n";
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_MAX_MATCHES_TO_FIND == att->name())
    {
      int mmf;
      if (! att->value(mmf) || mmf <= 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_MAX_MATCHES_TO_FIND <<
                "' attribute requires a whole positive number\n";
        return 0;
      }

      _max_matches_to_find = mmf;
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_SAVE_HITS_ATTRIBUTE == att->name())
    {
      int sha;
      if (! att->value(sha) || sha < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_SAVE_HITS_ATTRIBUTE <<
                "' attribute requires a whole non-negative number\n";
        return 0;
      }

      _save_matched_atoms = sha;
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_NCON_IGNORE_SINGLY_CONNECTED == att->name())
    {
      int nisc;
      if (! att->value(nisc) || nisc < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_NCON_IGNORE_SINGLY_CONNECTED <<
                "' attribute requires a whole number\n";
        return 0;
      }
      _ncon_ignore_singly_connected = nisc;

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_DO_NOT_PERCEIVE_SYMMETRIC_EQUIVALENTS == att->name())
    {
      if (! att->value(_do_not_perceive_symmetry_equivalent_matches))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_DO_NOT_PERCEIVE_SYMMETRIC_EQUIVALENTS <<
                "' attribute requires a whole number\n";
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_IMPLICIT_RING_CONDITION == att->name())
    {
      if (! att->value(_implicit_ring_condition))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_IMPLICIT_RING_CONDITION <<
                "' attribute requires a whole number\n";
        return 0;
      }

      if (_implicit_ring_condition > 1)
        _implicit_ring_condition = 1;

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! att->value(_all_hits_in_same_fragment))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE <<
                "' attribute requires a whole number\n";
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! att->value(_only_keep_matches_in_largest_fragment))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE <<
                "' attribute requires a whole number\n";
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE == att->name())
    {
      if (! att->value(_embeddings_do_not_overlap))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE <<
                "' attribute requires a whole number\n";
        return 0;
      }

      _find_one_embedding_per_start_atom = 1;    // by implication

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_SORT_BY_PREFERENCE_VALUE_ATTRIBUTE == att->name())
    {
      if (! att->value(_sort_by_preference_value))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_SORT_BY_PREFERENCE_VALUE_ATTRIBUTE <<
                "' attribute requires a whole number\n";
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_SMILES_ATTRIBUTE == att->name())
    {
      Substructure_Atom * a = new Substructure_Atom;
      if (! a->parse_smiles_specifier(att))
      {
        delete a;
        cerr << "Single_Substructure_Query::_construct_from_msi_object: cannot parse '" << NAME_OF_SMILES_ATTRIBUTE << "'\n";
        return 0;
      }

      _root_atoms.add (a);
    }
    else if (NAME_OF_SMARTS_ATTRIBUTE == att->name())
    {
      const_IWSubstring smarts;
      att->value(smarts);

      if (! create_from_smarts(smarts))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: cannot parse '" << NAME_OF_SMARTS_ATTRIBUTE << "'\n";
        return 0;
      }
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_NUMERIC_VALUE_ATTRIBUTE == att->name())
    {
      double d;
      if (! att->value(d))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the " << NAME_OF_NUMERIC_VALUE_ATTRIBUTE <<
                " must be numeric\n";
        cerr << (*att) << endl;
        return 0;
      }

      _numeric_value.add(d);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER == att->name())
    {
      if (2 != att->number_int_values())
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the " << NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER << " attribute must have two int values\n";
        cerr << (*att) << endl;
        return 0;
      }

      int a1 = att->int_multi_value(0);
      int a2 = att->int_multi_value(1);
      if (a1 < 0 || a2 < 0 || a1 == a2)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the " << NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER << " must have two valid numbers\n";
        cerr << (*att) << endl;
        return 0;
      }

      Bond * b = new Bond(a1, a2, SINGLE_BOND);

      _no_matched_atoms_between.add(b);

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_LINK_ATOMS_SPECIFIER == att->name())
    {
      Link_Atom * l = new Link_Atom;

      if (! l->construct_from_msi_attribute(att))
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object:invalid link atom specification\n";
        cerr << (*att) << endl;
        delete l;
        return 0;
      }

      _link_atom.add(l);
      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_DISTANCE_BETWEEN_HITS_NCHECK_ATTRIBUTE == att->name())
    {
      if (! att->value(_matched_atoms_to_check_for_hits_too_close) || _matched_atoms_to_check_for_hits_too_close < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object:the '" << NAME_OF_DISTANCE_BETWEEN_HITS_NCHECK_ATTRIBUTE << " must be a whole positive number\n";
        cerr << (*att) << endl;
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_FAIL_IF_EMBEDDINGS_TOO_CLOSE_ATTRIBUTE == att->name())
    {
      if (! att->value(_fail_if_embeddings_too_close) || _fail_if_embeddings_too_close < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object:the '" << NAME_OF_FAIL_IF_EMBEDDINGS_TOO_CLOSE_ATTRIBUTE << " must be a whole positive number\n";
        cerr << (*att) << endl;
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_ENVIRONMENT_MUST_MATCH_UNMATCHED_ATOMS_ATTRIBUTE == att->name())
    {
      if (! att->value(_environment_must_match_unmatched_atoms) || _environment_must_match_unmatched_atoms < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object:the '" << NAME_OF_ENVIRONMENT_MUST_MATCH_UNMATCHED_ATOMS_ATTRIBUTE << " must be non-negative\n";
        cerr << (*att) << endl;
        return 0;
      }

      attribute_recognised[i] = 1;
    }
    else if (NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS == att->name())
    {
      if (! att->value(environments_can_share_attachment_points) || environments_can_share_attachment_points < 0)
      {
        cerr << "Single_Substructure_Query::_construct_from_msi_object: the '" << NAME_OF_ENVIRONMENTS_CAN_SHARE_ATTACHMENT_POINTS << " attribute must be a whole non negative number\n";
        return 0;
      }
    }
    else if (NAME_OF_SORT_MATCHES_BY == att->name())
    {
      _sort_matches_by = 0;

      const_IWSubstring s;
      att->value(s);

      int j = 0;
      const_IWSubstring token;
      while (s.nextword(token, j))
      {
        int direction = 1;
        if (token.ends_with('+'))
          token.chop();
        else if (token.ends_with('-'))
        {
          token.chop();
          direction = -1;
        }
        if ("ncon" == token)
        {
          if (1 == direction)
            _sort_matches_by |= SORT_MATCHES_BY_INCREASING_NCON;
          else
            _sort_matches_by |= SORT_MATCHES_BY_DECREASING_NCON;
        }
        else if ("maxd" == token)
        {
          if (1 == direction)
            _sort_matches_by |= SORT_MATCHES_BY_INCREASING_MAXD;
          else
            _sort_matches_by |= SORT_MATCHES_BY_DECREASING_MAXD;
        }
        else
        {
          cerr << "Single_Substructure_Query:construct_from_msi_attribute:unrecognised sort attribute '" << s << "'\n";
          return 0;
        }
      }
    }
  }

  if (! _hits_needed.ok())
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object:inconsistent hits needed values\n";
    return 0;
  }

  if (_only_keep_matches_in_largest_fragment && _all_hits_in_same_fragment)
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object:the " << NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE << " and " << NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE << " attributes are mutually inconsistent\n";
    return 0;
  }

// For now, this is just a warning. Perhaps it should be fatal

  if (0 == version_attribute_found)
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object: no version attribute\n";
  }

  int nat = 0;
  if (! really_gruesome(_attached_heteroatom_count, msi, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, 
                          nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (_hits_needed.is_set())   // from a numeric qualifier on a smarts
    ;
  else if (! really_gruesome(_hits_needed, msi, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE,
                         nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_ring_atoms_matched, msi, NAME_OF_RING_ATOMS_MATCHED_ATTRIBUTE,
                         nat, 1, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_heteroatoms_matched, msi, NAME_OF_HETEROATOMS_MATCHED_ATTRIBUTE,
                         nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_heteroatoms_in_molecule, msi, NAME_OF_HETEROATOMS_ATTRIBUTE,
                         nat, 0, MAX_NOT_SPECIFIED))
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object: Cannot parse '" << NAME_OF_HETEROATOMS_MATCHED_ATTRIBUTE << "' attribute(s)\n";
    return 0;
  }

  if (_heteroatoms.number_elements() && 
      (! _attached_heteroatom_count.is_set() && ! _heteroatoms_in_molecule.is_set()))
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object: " << _heteroatoms.number_elements() <<
            " heteroatoms defined, but no attached_heteroatom_count\n";
    return 0;
  }

  if (! really_gruesome(_natoms, msi, NAME_OF_NATOMS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_nrings, msi, NAME_OF_NRINGS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_ncon, msi, NAME_OF_NCON_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_fused_rings, msi, NAME_OF_FUSED_RING_COUNT_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_strongly_fused_rings, msi, NAME_OF_STRONGLY_FUSED_RING_COUNT_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_isolated_rings, msi, NAME_OF_ISOLATED_RING_COUNT_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_isolated_ring_objects, msi, NAME_OF_ISOLATED_RING_OBJECTS, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_aromatic_rings, msi, NAME_OF_AROMATIC_RING_COUNT_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_non_aromatic_rings, msi, NAME_OF_NON_AROMATIC_RING_COUNT_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_distance_between_hits, msi, NAME_OF_DISTANCE_BETWEEN_HITS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_number_isotopic_atoms, msi, NAME_OF_NUMBER_ISOTOPIC_ATOMS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_number_fragments, msi, NAME_OF_NUMBER_FRAGMENTS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_distance_between_root_atoms, msi, NAME_OF_DISTANCE_BETWEEN_ROOT_ATOMS_SPECIFIER, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_atoms_in_spinach, msi, NAME_OF_SPINACH_ATOMS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_inter_ring_atoms, msi, NAME_OF_INTER_RING_ATOMS_ATTRIBUTE, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_unmatched_atoms, msi, NAME_OF_UNMATCHED_ATOMS, nat, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_net_formal_charge, msi, NAME_OF_NET_FORMAL_CHARGE_ATTRIBUTE, nat, MIN_NOT_SPECIFIED,  MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! get_float_attribute(msi, NAME_OF_MIN_FRACTION_ATOMS_MATCHED, _min_fraction_atoms_matched, 0.0, 1.0))
    return 0;

  if (! get_float_attribute(msi, NAME_OF_MAX_FRACTION_ATOMS_MATCHED, _max_fraction_atoms_matched, 0.0, 1.0))
    return 0;

  extending_resizable_array<Substructure_Atom *> completed;

// Since the environment may refer to query atoms, we cannot build it
// until the rest of the query is built

  int query_environment_present = 0;

  int nmsi = msi.number_elements();
  for (int i = 0; i < nmsi; i++)
  {
    msi_object & mi = *(msi[i]);

    if (NAME_OF_QUERY_BOND_OBJECT == mi.name())
    {
      cerr << "In version 2.0, bonds cannot be part of the query\n";
      return 0;
    }

    if (NAME_OF_ENVIRONMENT_OBJECT == mi.name())
    {
      query_environment_present = 1;
      continue;
    }

    if (NAME_OF_ENVIROMENT_REJECTION_OBJECT == mi.name())
    {
      query_environment_present = 1;
      continue;
    }

    if (NAME_OF_RING_SPECIFIER_OBJECT == mi.name())
    {
      if ( ! _parse_ring_specifier_object(mi))
        return 0;

      nat++;
      continue;
    }

    if (NAME_OF_RING_SYSTEM_SPECIFIER_OBJECT == mi.name())
    {
      if ( ! _parse_ring_system_specifier_object(mi))
        return 0;

      nat++;
      continue;
    }

    if (NAME_OF_ELEMENT_HITS_NEEDED_OBJECT == mi.name())
    {
      if (! _parse_element_hits_needed_object(mi))
        return 0;

      nat++;
      continue;
    }

    if (NAME_OF_ELEMENTS_NEEDED_OBJECT == mi.name())
    {
      if (! _parse_elements_needed_object(mi))
        return 0;

      nat++;
      continue;
    }

    if (NAME_OF_QUERY_ATOM_OBJECT != mi.name())
    {
      cerr << "Unrecognised type for msi object, ignored\n";
      cerr << mi;
      continue;
    }

    if (is_root_substructure_atom(mi))
    {
      Substructure_Atom * r = new Substructure_Atom;

      if (! r->construct_from_msi_object(mi, completed))
      {
        delete r;
        return 0;
      }

      nat++;
      _root_atoms.add(r);
    }
    else if (0 == _root_atoms.number_elements())
    {
      cerr << "Single_Substructure_Query::_construct_from_msi_object: no root atoms defined\n";
      return 0;
    }
    else
    {
      Substructure_Atom * a = new Substructure_Atom;
      if (! a->construct_from_msi_object(mi, completed))
      {
        delete a;
        return 0;
      }

      if (0 == a->nbonds())     // no bonds to anything already defined
      {
        cerr << "Non root Substructure_Atom not bonded\n";
        delete a;
        cerr << mi;
        return 0;
      }
    }
  }

// OK if the root is not built, but let's make sure there are some global
// attributes available in that case

  if (0 == _root_atoms.number_elements())
  {
    cerr << "Single_Substructure_Query::_construct_from_msi_object: Root atom not found\n";

//  apr 99. There are some global conditions which don't make sense with no root atom

    if (_element_hits_needed.number_elements() || _distance_between_root_atoms.is_set() || 
        _no_matched_atoms_between.number_elements() || _attached_heteroatom_count.is_set() ||
        _ncon.is_set() || _heteroatoms_matched.is_set() || _ring_atoms_matched.is_set())
    {
      cerr << "Single_Substructure_Query::_construct_from_msi_object: global conditions incompatible with no root atom\n";
      return 0;
    }

//  Should list all the possible match conditions

    if (nat > 0)
      ;
    else if (_natoms.is_set() || _nrings.is_set() || _aromatic_rings.is_set() || _atoms_in_spinach.is_set() || _inter_ring_atoms.is_set() || _net_formal_charge.is_set())
      ;
    else
    {
      cerr << "And no global attributes specified, rejected\n";
      return 0;
    }
  }

  if (query_environment_present)
  {
    for (int i = 0; i < nmsi; i++)
    {
      msi_object & mi = *(msi[i]);

      if (NAME_OF_ENVIRONMENT_OBJECT == mi.name())
      {
        if (! _construct_environment_from_msi_object(mi, completed, _environment))
        {
          cerr << "Single_Substructure_Query::_construct_from_msi_object: environment interpretation failed\n";
          cerr << mi << endl;
          return 0;
        }
      }

      else if (NAME_OF_ENVIROMENT_REJECTION_OBJECT == mi.name())
      {
        if (! _construct_environment_from_msi_object(mi, completed, _environment_rejections))
        {
          cerr << "Single_Substructure_Query::_construct_from_msi_object: environment interpretation failed\n";
          cerr << mi << endl;
          return 0;
        }
      }
    }

    if (environments_can_share_attachment_points >= 0)
    {
      for (int i = 0; i < _environment.number_elements(); ++i)
      {
        _environment[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
      }
      for (int i = 0; i < _environment_rejections.number_elements(); ++i)
      {
        _environment_rejections[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
      }
    }
  }

//if (_root_atoms.number_elements())
//  _adjust_for_internal_consistency();

  return 1;
}

int
Single_Substructure_Query::construct_from_msi_object(const msi_object & msi)
{
  assert (ok());

  if (NAME_OF_QUERY_OBJECT != msi.name())
  {
    cerr << "Single_Substructure_Query::construct_from_msi_object: msi object of incorrect type\n";
    return 0;
  }

// We need an array to mark each attribute as recognised or not.

  int * attribute_recognised = NULL;
  int nat = msi.attribute_count();
  if (nat)
    attribute_recognised = new_int(nat);

  int rc = _construct_from_msi_object(msi, attribute_recognised);

// Now that all the atoms have been defined, we can read in any chirality

  int n = msi.attribute_count(NAME_OF_CHIRALITY_SPECIFIER);
  if (n)
  {
    for (int i = 0; i < msi.number_attributes(); i++)
    {
      const msi_attribute * att = msi.attribute(i);
      if (NAME_OF_CHIRALITY_SPECIFIER != att->name())
        continue;

      attribute_recognised[i] = 1;
      if (! _build_chirality_specification_from_msi_attribute(att->stringval()))
      {
        cerr << "Single_Substructure_Query::construct_from_msi_object:invalid chirality '" << (*att) << "'\n";
        return 0;
      }
    }
  }

  if (rc && attribute_recognised)
  {
    for (int i = 0; i < nat; i++)
    {
      if (0 == attribute_recognised[i])
      {
//      implement this some time. The problem is things like 'min_....' are not recognised yet
//      const msi_attribute * att = msi.msi_attribute (i);
//      cerr << "Unrecognised attribute ";
//      rc = 0;
      }
    }
  }

  if (NULL != attribute_recognised)
    delete [] attribute_recognised;

  _preferences_present = 0;

  int nr = _root_atoms.number_elements();
  for (int i = 0; i < nr; i++)
  {
    if (_root_atoms[i]->preferences_present())
    {
      _preferences_present = 1;
      break;
    }
  }

  _fragment_ids_present = 0;
  for (int i = 0; i < nr; i++)
  {
    if (_root_atoms[i]->fragment_ids_present())
    {
      _fragment_ids_present = 1;
      break;
    }
  }

  assert (ok());

  (void) min_atoms_in_query();

  for (int i = 0; i < _no_matched_atoms_between.number_elements(); i++)
  {
    const Bond * b = _no_matched_atoms_between[i];
    if (b->a1() >= _min_atoms_in_query || b->a2() >= _min_atoms_in_query)
    {
      cerr << "Single_Substructure_Query::construct_from_msi_object: illegal " << NAME_OF_NO_MATCHED_ATOMS_BETWEEN_SPECIFIER << " specifier\n";
      cerr << "There are as few as " << _min_atoms_in_query << " query atoms\n";
      cerr << "No matched atoms between " << i << " specifies atoms " << b->a1() << " and " << b->a2() << endl;
      cerr << "this is impossible\n";
      return 0;
    }
  }

  return rc;
}

int
Single_Substructure_Query::read (iwstring_data_source & input)
{
  assert (input.ok());

  input.set_ignore_pattern("^#");
  input.set_skip_blank_lines(1);

  msi_object msi;

  if (! msi.read(input))
  {
    return 0;
  }

  input.set_ignore_pattern("");

  return construct_from_msi_object(msi);
}

int
Single_Substructure_Query::read (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return read(input);
}

int
Single_Substructure_Query::read (const IWString & fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return read(input);
}

/*
  Components are separated by the not-a-bond character, '.'
*/

int
Single_Substructure_Query::create_from_smiles (const IWString & smiles)
{
  if (! smiles.contains('.'))
  {
    Substructure_Atom * a = new Substructure_Atom;

    if (! a->parse_smiles_specifier(smiles))
    {
      delete a;
      return 0;
    }

    _root_atoms.add(a);

    return 1;
  }

// The more complex case of multiple components

  int i = 0;
  const_IWSubstring token;
  while (smiles.nextword(token, i))
  {
    Substructure_Atom * a = new Substructure_Atom;

    if (! a->parse_smiles_specifier(smiles))
    {
      delete a;
      return 0;
    }

    _root_atoms.add(a);
  }

  return _root_atoms.number_elements();
}

int
Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers (const_IWSubstring & smarts)
{
  int m;      // the numeric precondition
  int nchars;

  if ('>' == smarts[0])
  {
    smarts++;

    if ((0 == (nchars = fetch_numeric(smarts, m)) )|| (m < 0))
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers: the '>' directive must be followed by a whole number\n";
      return 0;
    }

    if (_hits_needed.is_set())
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers:_hits_needed already set\n";
      return 0;
    }

    if (! _hits_needed.set_min(m + 1))
      return 0;

    smarts += nchars;
    nchars++;     // to account for the > at the beginning
  }
  else if ('<' == smarts[0])
  {
    smarts++;

    if ((0 == (nchars = fetch_numeric(smarts, m))) || (m < 0))
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers: the '<' directive must be followed by a whole number\n";
      return 0;
    }

    if (_hits_needed.is_set())
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers:_hits_needed already set\n";
      return 0;
    }

    if (! _hits_needed.set_max(m - 1))
      return 0;

    smarts += nchars;
    nchars++;     // to account for the < at the beginning
  }
  else if ((nchars = fetch_numeric(smarts, m)))
  {
    if (m < 0)
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers: the count directive must be a whole positive number\n";
      return 0;
    }

    if (_hits_needed.is_set())
    {
      cerr << "Single_Substructure_Query::_parse_and_consume_optional_leading_numeric_qualifiers:_hits_needed already set\n";
      return 0;
    }

    if (! _hits_needed.add(m))
      return 0;

    smarts += nchars;
  }

  return 1;
}

//#define DEBUG_SSQ_CREATE_FROM_SMARTS

/*
  If the smarts starts with >nn or <nn these are interpreted as
  specifying the minimum or maximum number of hits needed.
  If it starts with a digit, that is a specific value for nhits.

  Valid inputs might be

  'CC'
  '0CC'
  '<3N-C.Cl'
  'Br...C(=O)Cl'

  Basically, this just strips off the initial numeric qualifiers and
  passes the rest to _parse_smarts_specifier
*/

int
Single_Substructure_Query::create_from_smarts (const IWString & smarts)
{
  if (smarts.nwords() > 1)
  {
    _comment = smarts;
    _comment.remove_leading_words(1);
  }

#ifdef DEBUG_SSQ_CREATE_FROM_SMARTS
  cerr << "Single_Substructure_Query::create_from_smarts: smarts is '" << smarts << "'\n";
#endif

  const_IWSubstring mysmarts = smarts;    // we may want to change it

  if (! _parse_and_consume_optional_leading_numeric_qualifiers(mysmarts))
  {
    cerr << "Single_Substructure_Query::create_from_smarts:invalid numeric specification '" << smarts << "'\n";
    return 0;
  }

#ifdef DEBUG_SSQ_CREATE_FROM_SMARTS
  cerr << "After parsing >< tokens, smarts is '" << mysmarts << "'\n";
#endif

  if (0 == mysmarts.nchars())
  {
    cerr << "Single_Substructure_Query::create_from_smarts: no characters to process '" << smarts << "'\n";
    return 0;
  }

  if ('!' == smarts[0])    // hmmm, why not use a leading 0. Is this correct?
  {
    _rejection = 1;
    mysmarts++;
  }

  if (0 == mysmarts.nchars())
  {
    cerr << "Single_Substructure_Query::create_from_smarts: no characters to process '" << smarts << "'\n";
    return 0;
  }

// Keeping the smarts order makes many things easier later on

  _respect_initial_atom_numbering = 1;

  return _parse_smarts_specifier(mysmarts);
}

int
Substructure_Query::write_msi (std::ostream & os, int & object_id, int indentation)
{
  assert (ok());
  assert (os.good());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->assign_unique_numbers();
  }

  if (1 == _number_elements)
    return _things[0]->write_msi(os, object_id, "", indentation);

// We have a composite query;

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << ' ' << NAME_OF_COMPOSITE_QUERY_OBJECT << endl;

  os << ind << "  (A I " << NAME_OF_VERSION_ATTRIBUTE << " 2)\n";

  if (_comment.length())
    os << ind << "  (A C " << NAME_OF_COMMENT_ATTRIBUTE << " \"" << _comment << "\")\n";

  if (_each_component_search)
    os << ind << "  (A I " << NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE << ' ' << _each_component_search << ")\n";
    
  if (_operator.all_operators_are('|'))
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (! _things[i]->write_msi(os, object_id, "", indentation + 2))
        return 0;
    }
  }
  else       // need to write operators
  {
    for (int i = 0; i < _number_elements; i++)
    {
      const_IWSubstring op;
      if (i < _number_elements - 1)
      {
        _operator.op(i, op);
        if ('&' == op)
          op = NAME_OF_AND_OPERATOR;
        else if ('|' == op)
          op = NAME_OF_OR_OPERATOR;
        else if ('^' == op)
          op = NAME_OF_XOR_OPERATOR;
        else if (';' == op)
          op = NAME_OF_LOW_PRIORITY_AND_OPERATOR;
        else
        {
          cerr << "Huh, what operator is this '" << op << "'\n";
          abort();
        }
      }

      if (! _things[i]->write_msi(os, object_id, op, indentation + 2))
        return 0;
    }
  }

  os << ind << ")\n";

  return os.good();
}

int
Substructure_Query::write_msi (std::ostream & os)
{
  assert (ok());
  assert (os.good());

  int i = 0;

  return write_msi(os, i, 0);
}

int
Substructure_Query::write_msi (IWString & fname)
{
  std::ofstream os(fname.null_terminated_chars(), std::ios::out);
  if (! os.good())
  {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_msi(os);
}

int
Substructure_Query::read (iwstring_data_source & input)
{
  assert (input.ok());

  input.set_ignore_pattern("^#");

  msi_object msi;

  if (! msi.read(input))
  {
    return 0;
  }

  input.set_ignore_pattern("");

  return construct_from_msi_object(msi);
}

static MDL_Molecule * 
read_first_molecule_from_file (const const_IWSubstring & fname)
{
  data_source_and_type<MDL_Molecule> input(MDL, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return input.next_molecule();
}

/*
  The argument to the ::read functions takes many roles.
  By default, it is the file name from which to read the query.
  If it starts with SMARTS: however, it is interpreted as a smarts

  To enforce this, we have a common entry point
*/

int
Substructure_Query::read (const const_IWSubstring & zname)
{
  if (zname.starts_with("SMARTS:"))
  {
    const_IWSubstring smarts = zname;    // our arg is const
    smarts.remove_leading_chars(7);

    if (! create_from_smarts(smarts))
    {
      cerr << "Substructure_Query::read: cannot parse 'SMARTS:" << smarts << "'\n";
      return 0;
    }

    return 1;
  }
  else if (zname.starts_with("SMILES:"))
  {
    const_IWSubstring smiles = zname;    // our arg is const
    smiles.remove_leading_chars(7);

    Single_Substructure_Query * q = new Single_Substructure_Query;
    if (! q->create_from_smiles(smiles))
    {
      cerr << "Substructure_Query::read: cannot parse '" << zname << "'\n";
      delete q;
      return 0;
    }

    add(q);

    if (_number_elements > 1)
      _operator.add_operator('|');

    return 1;
  }
  else if (zname.starts_with("I:"))
  {
    const_IWSubstring fname(zname);
    fname.remove_leading_chars(2);
    MDL_Molecule * m = read_first_molecule_from_file(fname);
    if (NULL == m)
    {
      cerr << "Cannot read query molecule from '" << fname << "'\n";
      return 0;
    }
    Molecule_to_Query_Specifications mqs;
    Single_Substructure_Query * q = new Single_Substructure_Query;
    if (! q->create_from_molecule(*m, mqs))
    {
      cerr << "Substructure_Query::read:cannot create query from molecule in '" << fname << "'\n";
      delete q;
      return 0;
    }

    add(q);
    if (_number_elements > 1)
      _operator.add_operator('|');

    return 1;
  }

  iwstring_data_source input(zname);
  if (! input.ok())
  {
    cerr << "Substructure_Query::read: cannot open '" << zname << "'\n";
    return 0;
  }

  return read(input);
}

int
Substructure_Query::read (const IWString & fname)
{
  const_IWSubstring tmp(fname);

  return read(tmp);
}

int
Substructure_Query::read (const char * fname)
{
  const_IWSubstring tmp(fname);

  return read(tmp);
}

int
Substructure_Query::create_from_smiles (const IWString & smiles)
{
  Single_Substructure_Query * q = new Single_Substructure_Query;

  if (! q->create_from_smiles(smiles))
  {
    delete q;
    return 0;
  }

  add(q);

  if (_number_elements > 1)
    _operator.add_operator('|');

  return 1;
}

//#define DEBUG_CREATE_FROM_SMARTS

/*
  Called from Substructure_Query::create_from_smarts to process the
  components of a smarts
*/

int
Substructure_Query::_parse_smarts_components (const IWString & smarts,
                                              char separator)
{
  int rc = 0;

  int i = 0;
  const_IWSubstring token;
  while (smarts.nextword(token, i, separator))
  {
    Single_Substructure_Query * q = new Single_Substructure_Query;

    if (! q->create_from_smarts(token))
    {
      delete q;
      return 0;
    }

    add(q);

    rc++;
  }

  return rc;
}

static int
fetch_next_smarts_component(const const_IWSubstring & smarts,
                            int smartslen,
                            int & ndx,
                            const_IWSubstring & s)
{
#ifdef DEBUG_CREATE_FROM_SMARTS
  cerr << "On entry " << smartslen << " chars to process, ndx = " << ndx << "\n";
#endif

  if (ndx >= smartslen)
    return 0;

  int istop = smartslen - 2;

  int istart = ndx;

  while (ndx <= istop)
  {
    char c1 = smarts[ndx];

    ndx++;

    if ('&' == c1 || '|' == c1 || '^' == c1 || ';' == c1 || ',' == c1)
      ;
    else 
      continue;

    char c2 = smarts[ndx];
    if (c1 != c2)
      continue;

    if (1 == ndx)
    {
      cerr << "smarts cannot start with operator\n";
      return 0;
    }

    if (c1 == c2)
    {
      ndx--;
      s = smarts.from_to(istart, ndx - 1);

#ifdef DEBUG_CREATE_FROM_SMARTS
      cerr << "Component is '" << s << "', ndx = " << ndx << "\n";
#endif
      return 1;
    }
  }

// If we come to here, then the whole rest of the string is the smarts token

  s = smarts.from_to(istart, smartslen - 1);

#ifdef DEBUG_CREATE_FROM_SMARTS
  cerr << "Whole lot '" << s << "'\n";
#endif

  return s.length();
}

static int
fetch_next_smarts_operator(const const_IWSubstring & smarts,
                           int & ndx,
                           int & op)
{
  if (ndx + 2 > smarts.length())   // all operators length 2
    return 0;

  char c1 = smarts[ndx];
  char c2 = smarts[ndx + 1];

#ifdef DEBUG_CREATE_FROM_SMARTS
  cerr << " Operator? '" << c1 << "' and '" << c2 << "'\n";
#endif

  if (c1 != c2)
    return 0;

  if ('&' == c1)
    op = IW_LOGEXP_AND;
  else if ('|' == c1 || ',' == c1)
    op = IW_LOGEXP_OR;
  else if ('^' == c1)
    op = IW_LOGEXP_XOR;
  else if (';' == c1)
    op = IW_LOGEXP_LOW_PRIORITY_AND;
  else
    return 0;

  ndx += 2;

  return 1;
}

int
Substructure_Query::_add_component_from_smarts (const const_IWSubstring & smarts,
                                                int op)
{
  Single_Substructure_Query * q = new Single_Substructure_Query;

  if (! q->create_from_smarts(smarts))
  {
    delete q;
    return 0;
  }

  resizable_array_p<Single_Substructure_Query>::add(q);

  if (IW_LOGEXP_UNDEFINED != op)
    _operator.add_operator(op);

#ifdef DEBUG_CREATE_FROM_SMARTS
    _operator.debug_print(cerr);
#endif

  return 1;
}

/*
  As an extension, a smarts can have multiple components separated
  by '&&' (and operator), '||' (or operator), '^^' (xor operator) or
  ';;' (low priority and operator).

  Note also that we don't check whether or not we have previously
  been initialised. Not sure if this is a bug or feature!
*/

int
Substructure_Query::create_from_smarts (const IWString & smarts)
{
  int smartslen;
  if (smarts.nwords() > 1)
  {
    _comment = smarts;
    _comment.remove_leading_words(1);
    smartslen = smarts.index(' ');
  }
  else
    smartslen = smarts.length();

#ifdef DEBUG_CREATE_FROM_SMARTS
  cerr << "Substructure_Query::create_from_smarts: smarts is '" << smarts << "'\n";
#endif

  int ndx = 0;
  const_IWSubstring s;
  if (! fetch_next_smarts_component(smarts, smartslen, ndx, s))
  {
    cerr << "Substructure_Query::create_from_smarts:no smarts\n";
    return 0;
  }

  if (! _add_component_from_smarts(s, IW_LOGEXP_UNDEFINED))
    return 0;

  int op = 0;    // initialise to keep the compiler quiet
  while (fetch_next_smarts_operator(smarts, ndx, op))
  {
    if (! fetch_next_smarts_component(smarts, smartslen, ndx, s))
    {
      cerr << "Substructure_Query::create_from_smarts:no smarts after operator\n";
      return 0;
    }

#ifdef DEBUG_CREATE_FROM_SMARTS
    cerr << "Smarts component is '" << s << "', operator " << op << "\n";
#endif

    if (! _add_component_from_smarts(s, op))
      return 0;
  }

#ifdef DEBUG_CREATE_FROM_SMARTS
  cerr << "After parsing " << number_elements() << " components, operator is\n";
  _operator.debug_print(cerr);
#endif

  return 1;
}

int
Substructure_Query::_single_query_construct_from_msi_object (const msi_object & msi)
{
  Single_Substructure_Query * tmp = new Single_Substructure_Query;
  if (! tmp->construct_from_msi_object(msi))
  {
    delete tmp;
    return 0;
  }

  add(tmp);

// Does this query have an operator associated with it

  IWString op;
  if (! msi.string_value_for_attribute(NAME_OF_OPERATOR_ATTRIBUTE, op))
  {
    if (_number_elements > 1)
      _operator.add_operator('|');     // the OR operator is the default

    return 1;
  }

  if (NAME_OF_AND_OPERATOR == op)
    _operator.add_operator('&');
  else if (NAME_OF_OR_OPERATOR == op)
    _operator.add_operator('|');
  else if (NAME_OF_XOR_OPERATOR == op)
    _operator.add_operator('^');
  else if (NAME_OF_LOW_PRIORITY_AND_OPERATOR == op)
    _operator.add_operator(';');
  else
  {
    cerr << "Substructure_Query::_single_query_construct_from_msi_object: which operator is this '" << op << "'\n";
    return 0;
  }

  return 1;
}

/*
  We are not entirely rigorous about parsing the msi object.
  For example, we would not detect the error in
    (A C smarts "O|N")
    (A C smarts "P")
    (A I operator and)
*/

int
Substructure_Query::_composite_query_construct_from_msi_object (const msi_object & msi)
{
  assert (0 == _operator.number_operators());

  int nmsi = msi.number_elements();
  const msi_attribute * att = msi.attribute(NAME_OF_SMARTS_ATTRIBUTE);

  if (0 == nmsi && NULL == att)
  {
    cerr << "Substructure_Query::_composite_query_construct_from_msi_object: no objects or smarts\n";
    return 0;
  }

// Looks like we should do something about operators here

  int aptr = 0;
  while (NULL != (att = msi.attribute(NAME_OF_SMARTS_ATTRIBUTE, aptr++)))
  {
    IWString smarts;
    att->value(smarts);

    if (! create_from_smarts(smarts))
    {
      cerr << "Substructure_Query::_composite_query_construct_from_msi_object: cannot parse smarts\n";
      return 0;
    }
  }

// Each query may have an operator associated with it. When we read in a
// query we need to add the operator that was written with the previous query

  int pending_operator = IW_LOGEXP_UNDEFINED;

#ifdef DEBUG_READ_COMPOSITE
  cerr << "Processing " << nmsi << " component query\n";
#endif

  for (int i = 0; i < nmsi; i++)
  {
    Single_Substructure_Query * tmp = new Single_Substructure_Query;
    if (! tmp->construct_from_msi_object(*(msi[i])))
    {
      delete tmp;
      return 0;
    }

//  The add() member function takes care of operators

#ifdef DEBUG_READ_COMPOSITE
    cerr << "Components " << _number_elements << " operators " << _operator.number_operators() << endl;
#endif

    if (IW_LOGEXP_UNDEFINED != pending_operator)
      add(tmp, pending_operator);
    else
      add(tmp);     // the OR operator is the default

#ifdef DEBUG_READ_COMPOSITE
    cerr << "Components " << _number_elements << " operators " << _operator.number_operators() << endl;
#endif

    pending_operator = IW_LOGEXP_UNDEFINED;

//  The last component doesn't include an operator.

    if (i == nmsi - 1)
      break;

    IWString op;
    if (! msi[i]->string_value_for_attribute(NAME_OF_OPERATOR_ATTRIBUTE, op))
      continue;        // pending_operator already set to undefined above

    if (NAME_OF_AND_OPERATOR == op)
      pending_operator = IW_LOGEXP_AND;
    else if (NAME_OF_OR_OPERATOR == op)
      pending_operator = IW_LOGEXP_OR;
    else if (NAME_OF_XOR_OPERATOR == op)
      pending_operator = IW_LOGEXP_XOR;
    else if (NAME_OF_LOW_PRIORITY_AND_OPERATOR == op)
      pending_operator = IW_LOGEXP_LOW_PRIORITY_AND;
    else
    {
      cerr << "Substructure_Query::_composite_query_construct_from_msi_object: what operator is this '" << op << "'\n";
      abort();
    }
  }

  assert (_number_elements == _operator.number_operators() + 1);

  att = msi.attribute(NAME_OF_COMMENT_ATTRIBUTE);
  if (NULL != att)
    att->value(_comment);

  att = msi.attribute(NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE);
  if (NULL != att)
  {
    if (! att->value(_each_component_search))
    {
      cerr << "Substructure_Query::_composite_query_construct_from_msi_object:invalid value for " << NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE << " attribute\n";
      return 0;
    }
  }

  att = msi.attribute(NAME_OF_OPERATOR_ATTRIBUTE);
  if (NULL == att)
    return 1;

// Any operator specified in the composite (silently) overrides any operators
// which may have come from the components

  IWString op;
  att->value(op);

  if (NAME_OF_AND_OPERATOR == op)
    _operator.set_all_operators('&');
  else if (NAME_OF_OR_OPERATOR == op)
    _operator.set_all_operators('|');
  else if (NAME_OF_XOR_OPERATOR == op)
    _operator.set_all_operators('^');
  else if (NAME_OF_LOW_PRIORITY_AND_OPERATOR == op)
    _operator.set_all_operators(';');
  else
  {
    cerr << "Substructure_Query::_composite_query_construct_from_msi_object: what operator is this '" << op << "'\n";
    return 0;
  }

  return 1;
}

/*
  We must be able to read either a set of queries, or just one.
*/

int
Substructure_Query::construct_from_msi_object (const msi_object & msi)
{
  int rc = 0;

  const msi_attribute * att;

  if (NAME_OF_QUERY_OBJECT == msi.name())
  {
    rc = _single_query_construct_from_msi_object(msi);
    if (_number_elements)
      _comment = _things[0]->comment();
  }
  else if (NAME_OF_COMPOSITE_QUERY_OBJECT == msi.name())
  {
    rc = _composite_query_construct_from_msi_object(msi);
  }
  else if (NULL != (att = msi.attribute(NAME_OF_SMARTS_ATTRIBUTE)))
  {
    IWString smarts;
    att->value(smarts);

    if (! create_from_smarts(smarts))
    {
      cerr << "Substructure_Query::construct_from_msi_object: cannot parse smarts\n";
      return 0;
    }
  }
  else
  {
    cerr << "Substructure_Query::construct_from_msi_object: unrecognised msi object type '" << msi.name() << "'\n";
    cerr << msi;
    return 0;
  }

  if (0 == rc)
    resize(0);

  return rc;
}

int
Substructure_Ring_Specification::construct_from_msi_object (const msi_object & msi)
{
  int attributes_specified = 0;

  if (! Substructure_Ring_Base::construct_from_msi_object(msi, attributes_specified))
  {
    cerr << "Substructure_Ring_Specification::construct_from_msi_object: base failed\n";
    return 0;
  }

  if (! really_gruesome(_ring_size, msi, NAME_OF_RING_SIZE_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_fused, msi, NAME_OF_RING_FUSED_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_fused_aromatic_neighbours, msi, NAME_OF_FUSED_AROMATIC_NEIGHBOURS_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_fused_non_aromatic_neighbours, msi, NAME_OF_FUSED_NON_AROMATIC_NEIGHBOURS_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  const msi_attribute * att = msi.attribute(NAME_OF_AROMATICITY_ATTRIBUTE);

  if (att)
  {
    if (! att->value(_aromatic) || _aromatic< 0)
    {
      cerr << "Invalid aromaticity value\n";
      return 0;
    }
  }

  return 1;
}

int
Substructure_Ring_Specification::write_msi (std::ostream & os, int & object_id, int indentation) const
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << '(' << object_id++ << ' ' << NAME_OF_RING_SPECIFIER_OBJECT << endl;

  (void) Substructure_Ring_Base::write_msi_attributes(os, object_id, ind);

  _ring_size.write_msi(os, NAME_OF_RING_SIZE_ATTRIBUTE, indentation + 2);

  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromatic)
    os << ind << "  (A I " << NAME_OF_AROMATICITY_ATTRIBUTE << ' ' << _aromatic << ")\n";

  _fused.write_msi(os, NAME_OF_RING_FUSED_ATTRIBUTE, indentation + 2);
  _fused_aromatic_neighbours.write_msi(os, NAME_OF_FUSED_AROMATIC_NEIGHBOURS_SPECIFIER, indentation + 2);
  _fused_non_aromatic_neighbours.write_msi(os, NAME_OF_FUSED_NON_AROMATIC_NEIGHBOURS_SPECIFIER, indentation + 2);
  _largest_number_of_bonds_shared_with_another_ring.write_msi(os, NAME_OF_LARGEST_NUMBER_SHARED_BONDS_SPECIFIER, indentation + 2);
  _strongly_fused_ring_neighbours.write_msi(os, NAME_OF_STRONGLY_FUSED_RING_NEIGHBOUR_SPECIFIER, indentation + 2);

  os << ind << ")\n";

  return os.good();
}

int
Substructure_Ring_Base::construct_from_msi_object(const msi_object & msi,
                                                  int & attributes_specified)
{
  if (msi.number_elements())
  {
    cerr << "Substructure_Ring_Base::construct_from_msi_object:no objects allowed in specification\n";
    return 0;
  }

  int nat = msi.attribute_count();

  for (int i = 0; i < nat; i++)
  {
    const msi_attribute * att = msi.attribute(i);

    if (NAME_OF_REJECTION_ATTRIBUTE == att->name())
    {
      int tmp;
      if (! att->value(tmp))
      {
        cerr << "Substructure_Ring_Base::construct_from_msi_object: the " << NAME_OF_REJECTION_ATTRIBUTE <<
                " attribute must be a whole number\n";
        return 0;
      }

      _match_as_match_or_rejection = ! tmp;
    }
    else if (NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! att->value(_all_hits_in_same_fragment))
      {
        cerr << "Substructure_Ring_Specification::construct_from_msi_object: the " << NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE <<
                " attribute must be a whole number\n";
        return 0;
      }
      attributes_specified++;
    }
    else if (NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE == att->name())
    {
      if (! att->value(_only_keep_matches_in_largest_fragment))
      {
        cerr << "Substructure_Ring_Specification::construct_from_msi_object: the " << NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE << " attribute must be a whole number\n";
        return 0;
      }
      attributes_specified++;
    }
    else if (NAME_OF_COMMENT_ATTRIBUTE == att->name())
    {
      (void) att->value(_comment);
    }
    else if (NAME_OF_ENVIRONMENT_OBJECT == att->name())
    {
      if (_environment_atom.number_elements())
      {
        cerr << "Substructure_Ring_Base::construct_from_msi_object:environment already specified '" << (*att) << "'\n";
        return 0;
      }

      const_IWSubstring env;
      att->value(env);

      if (! _construct_environment(env))
      {
        cerr << "Substructure_Ring_Base::construct_from_msi_object:invalid environment '" << (*att) << "'\n";
        return 0;
      }

      attributes_specified++;
    }
    else if (NAME_OF_ENVIRONMENT_CAN_MATCH_IN_RING_ATOMS == att->name())
    {
      if (! att->value(_environment_can_match_in_ring_atoms))
      {
        cerr << "Substructure_Ring_Base::construct_from_msi_object:invalid " << NAME_OF_ENVIRONMENT_CAN_MATCH_IN_RING_ATOMS << " specification\n";
        return 0;
      }

      attributes_specified++;
    }
    else if (NAME_OF_DEFINE_HETEROATOMS_ATTRIBUTE == att->name())
    {
      if (! fetch_heteroatom_definitions(*att, _is_heteroatom))
      {
        cerr << "Substructure_Ring_Base::construct_from_msi_object:invalid heteroatom definitions '" << (*att) << "'\n";
        return 0;
      }
    }
  }

  if (! really_gruesome(_hits_needed, msi, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_ncon, msi, NAME_OF_NCON_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_heteroatom_count, msi, NAME_OF_HETEROATOMS_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_attached_heteroatom_count, msi, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_within_ring_unsaturation, msi, NAME_OF_UNSATURATION_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_atoms_with_pi_electrons, msi, NAME_OF_ATOMS_WITH_PI_ELECTRONS, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_largest_number_of_bonds_shared_with_another_ring, msi, NAME_OF_LARGEST_NUMBER_SHARED_BONDS_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  if (! really_gruesome(_strongly_fused_ring_neighbours, msi, NAME_OF_STRONGLY_FUSED_RING_NEIGHBOUR_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
    return 0;

  return 1;
}

int
Substructure_Ring_System_Specification::construct_from_msi_object(const msi_object & msi)
{
  int attributes_specified = 0;

  if (! Substructure_Ring_Base::construct_from_msi_object(msi, attributes_specified))
    return 0;

  if (_atoms_with_pi_electrons.is_set())
    _need_per_atom_array = 1;

  if (! really_gruesome(_rings_in_system, msi, NAME_OF_NRINGS_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_ring_sizes, msi, NAME_OF_RING_SIZE_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_aromatic_ring_count, msi, NAME_OF_AROMATIC_RING_COUNT_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_non_aromatic_ring_count, msi, NAME_OF_NON_AROMATIC_RING_COUNT_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_degree_of_fusion, msi, NAME_OF_RING_FUSED_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_atoms_in_system, msi, NAME_OF_NATOMS_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_strongly_fused_ring_count, msi, NAME_OF_STRONGLY_FUSED_RING_COUNT_SPECIFIER, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_number_spinach_groups, msi, NAME_OF_NUMBER_SPINACH_GROUPS, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_number_non_spinach_groups, msi, NAME_OF_NUMBER_NON_SPINACH_GROUPS, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_atoms_in_spinach_group, msi, NAME_OF_ATOMS_IN_SPINACH_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_length_of_spinach_group, msi, NAME_OF_LENGTH_OF_SPINACH_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_distance_to_another_ring, msi, NAME_OF_DISTANCE_TO_ANOTHER_RING_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  if (! really_gruesome(_rings_that_must_match_ring_sizes, msi, NAME_OF_RINGS_THAT_MUST_MATCH_RING_SIZE_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    return 0;
  }

  return attributes_specified;
}

int
Substructure_Ring_Base::write_msi_attributes (std::ostream & os, int & object_id, const const_IWSubstring & ind) const
{
  if (_comment.length())
    os << ind << "  (A C " << NAME_OF_COMMENT_ATTRIBUTE << " \"" << _comment << "\")\n";

  if (0 == _match_as_match_or_rejection)
    os << ind << "  (A I " << NAME_OF_REJECTION_ATTRIBUTE << " 1)\n";

  int indentation = ind.length();

  _hits_needed.write_msi(os, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, indentation + 2);
  _attached_heteroatom_count.write_msi(os, NAME_OF_ATTACHED_HETEROATOM_COUNT_SPECIFIER, indentation + 2);
  _heteroatom_count.write_msi(os, NAME_OF_HETEROATOMS_ATTRIBUTE, indentation + 2);
  _ncon.write_msi(os, NAME_OF_NCON_ATTRIBUTE, indentation + 2);
  _within_ring_unsaturation.write_msi(os, NAME_OF_UNSATURATION_ATTRIBUTE, indentation + 2);
  _atoms_with_pi_electrons.write_msi(os, NAME_OF_ATOMS_WITH_PI_ELECTRONS, indentation + 2);

  if (_all_hits_in_same_fragment)
    os << ind << "  (A I " << NAME_OF_ALL_HITS_IN_SAME_FRAGMENT_ATTRIBUTE << " 1)\n";
  if (_only_keep_matches_in_largest_fragment)
    os << ind << "  (A I " << NAME_OF_ONLY_MATCH_LARGEST_FRAGMENT_ATTRIBUTE << " 1)\n";

  if (_environment_atom.number_elements())
  {
    os << ind << "  (A C " << NAME_OF_ENVIRONMENT_OBJECT << " ****)\n";
    cerr << "Substructure_Ring_Base::write_msi_attributes:cannot write ring environment\n";
  }

  return os.good();
}

int
Substructure_Ring_System_Specification::write_msi(std::ostream & os, int & object_id, int indentation) const
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << '(' << object_id++ << ' ' << NAME_OF_RING_SYSTEM_SPECIFIER_OBJECT << endl;

  (void) Substructure_Ring_Base::write_msi_attributes(os, object_id, ind);

  _rings_in_system.write_msi(os, NAME_OF_NRINGS_ATTRIBUTE, indentation + 2);
  _ring_sizes.write_msi(os, NAME_OF_RING_SIZE_ATTRIBUTE, indentation + 2);
  _rings_that_must_match_ring_sizes.write_msi(os, NAME_OF_RINGS_THAT_MUST_MATCH_RING_SIZE_ATTRIBUTE, indentation + 2);
  _aromatic_ring_count.write_msi(os, NAME_OF_AROMATIC_RING_COUNT_SPECIFIER, indentation + 2);
  _non_aromatic_ring_count.write_msi(os, NAME_OF_NON_AROMATIC_RING_COUNT_SPECIFIER, indentation + 2);
  _degree_of_fusion.write_msi(os, NAME_OF_RING_FUSED_ATTRIBUTE, indentation + 2);
  _atoms_in_system.write_msi(os, NAME_OF_NATOMS_ATTRIBUTE, indentation + 2);
  _strongly_fused_ring_count.write_msi(os, NAME_OF_STRONGLY_FUSED_RING_COUNT_SPECIFIER, indentation + 2);
  _number_spinach_groups.write_msi(os, NAME_OF_NUMBER_SPINACH_GROUPS, indentation + 2);
  _number_non_spinach_groups.write_msi(os, NAME_OF_NUMBER_NON_SPINACH_GROUPS, indentation + 2);
  _atoms_in_spinach_group.write_msi(os, NAME_OF_ATOMS_IN_SPINACH_ATTRIBUTE, indentation + 2);
  _length_of_spinach_group.write_msi(os, NAME_OF_LENGTH_OF_SPINACH_ATTRIBUTE, indentation + 2);
  _distance_to_another_ring.write_msi(os, NAME_OF_DISTANCE_TO_ANOTHER_RING_ATTRIBUTE, indentation + 2);

  os << ind << ")\n";

  return os.good();
}

int
Elements_Needed::write_msi (std::ostream & os, 
                            const char * object_name,
                            int & object_id,
                            int indentation) const
{
  if (_z == INVALID_ATOMIC_NUMBER)
  {
    cerr << "Elements_Needed::write_msi: cannot write invalid atomic number\n";
    return 0;
  }

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << ' ' << object_name << endl;

  os << ind << "  (A I " << NAME_OF_ATOMIC_NUMBER_ATTRIBUTE << ' ' << _z << ")\n";
  Min_Max_Specifier<int>::write_msi(os, NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, indentation + 2);

  os << ind << ")\n";

  return os.good();
}

int
Single_Substructure_Query::_parse_element_hits_needed_object (const msi_object & msi)
{

  Elements_Needed * e = new Elements_Needed();
  if (! e->construct_from_msi_object(msi))
  {
    delete e;
    return 0;
  }

  _element_hits_needed.add(e);

  return 1;
}

int
Single_Substructure_Query::_parse_elements_needed_object (const msi_object & msi)
{

  Elements_Needed * e = new Elements_Needed();
  if (! e->construct_from_msi_object(msi))
  {
    delete e;
    return 0;
  }

  _elements_needed.add(e);

  return 1;
}

int
Elements_Needed::construct_from_msi_object (const msi_object & msi)
{
  if (msi.number_elements())
  {
    cerr << "Elements_Needed::construct_from_msi_object: object cannot be complex\n";
    return 0;
  }

  const msi_attribute * att = msi.attribute(NAME_OF_ATOMIC_NUMBER_ATTRIBUTE);
  if (NULL == att)
  {
    cerr << "Elements_Needed::construct_from_msi_object: must give '" << NAME_OF_ATOMIC_NUMBER_ATTRIBUTE << "' attribute\n";
    return 0;
  }

  if (! att->value(_z) || ! REASONABLE_ATOMIC_NUMBER(_z))
  {
    cerr << "Elements_Needed::construct_from_msi_object: bad " << NAME_OF_ATOMIC_NUMBER_ATTRIBUTE << " attribute\n";
    return 0;
  }

  int attributes_specified = 0;
  if (! really_gruesome((*this), msi,
                 NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE, attributes_specified, 0, MAX_NOT_SPECIFIED))
  {
    cerr << "Elements_Needed::construct_from_msi_object: cannot parse " << NAME_OF_QUERY_HITS_NEEDED_ATTRIBUTE << " attribute\n";
    return 0;
  }

  return 1;
}

int
Link_Atom::construct_from_msi_attribute (const msi_attribute * msi)
{
  const IWString & s = msi->stringval();

  if (s.nwords() < 3)
  {
    cerr << "Link_Atom::construct_from_msi_object:invalid specification, must have at least 3 words '" << s << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  s.nextword(token, i);

  if (! token.numeric_value(_a1) || _a1 < 0)
  {
    cerr << "Link_Atom::construct_from_msi_object:invalid a1 specification, must be a positive number'" << s << "'\n";
    return 0;
  }

  s.nextword(token, i);

  if (! token.numeric_value(_a2) || _a2 < 0 || _a2 == _a1)
  {
    cerr << "Link_Atom::construct_from_msi_object:invalid a2 specification, must be a positive number'" << s << "'\n";
    return 0;
  }

  s.nextword(token, i);

  if (! _d.initialise(token))
  {
    cerr << "Link_Atom::construct_from_msi_object:invalid distance specification '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
Link_Atom::initialise_from_smarts (const_IWSubstring const & s)
{
  if (s.starts_with('{') && s.ends_with('}'))
    ;
  else 
  {
    cerr << "Link_Atom::initialise_from_smarts:specifier must be enclosed with {}\n";
    return 0;
  }

  const_IWSubstring tmp(s);
  tmp.remove_leading_chars(1);
  tmp.chop(1);

  if (0 == tmp.length())
  {
    cerr << "Link_Atom::initialise_from_smarts:empty specification, impossible\n";
    return 0;
  }

  if (! _d.initialise(tmp))
  {
    cerr << "Link_Atom::initialise_from_smarts:invalid range specification '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

int
Link_Atom::write_msi (std::ostream & os, 
                      const IWString & ind,
                      const char * attname) const
{
  os << ind << "  (A C " << attname << " \"" << _a1 << ' ' << _a2 << ' ';
  _d.write_compact_representation(os);

  os << "\")\n";

  return os.good();
}

/*
  This is kind of broken because we remove link atoms the moment they
  are read in - that is, they are no longer atoms, so zatom will be
  bogus
*/

int
Link_Atom::write_M_LIN(atom_number_t zatom, std::ostream & output) const
{
  output << "M  LIN " << zatom << ' ' << _a1 << ' ' << _a2 << " implement this";
  output << '\n';

  return output.good();
}

int
Single_Substructure_Query::add_link_atom (const Link_Atom & ila)
{
  Link_Atom * l = new Link_Atom(ila);

  const Substructure_Atom * a = query_atom_with_initial_atom_number(ila.a1());

  cerr << "Initial atom number " << ila.a1() << " becomes query atom " <<  a->unique_id() << endl;

  l->set_a1(a->unique_id());

  a = query_atom_with_initial_atom_number(ila.a2());

  cerr << "Initial atom number " << ila.a2() << " becomes query atom " <<  a->unique_id() << endl;

  l->set_a2(a->unique_id());

  _link_atom.add(l);

  return 1;
}
