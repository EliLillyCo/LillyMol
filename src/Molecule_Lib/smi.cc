#include <iostream>
#include <memory>
#include <tuple>
#include <ctype.h>

#define COMPILING_SMILES_CC
#define COMPILING_CTB

#include <iostream>

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "Foundational/iwmisc/misc.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "aromatic.h"
#include "chiral_centre.h"
#include "coordinate_box.h"
#include "element.h"
#include "misc2.h"
#include "molecule.h"
#include "moleculeio.h"
#include "parse_smarts_tmp.h"
#include "rwmolecule.h"
#include "smiles.h"
#include "substructure.h"

using std::cerr;
using std::endl;

using down_the_bond::DownTheBond;
using moleculeio::newline_string;

constexpr char kOparen = '(';
constexpr char kCparen = ')';
constexpr char kOpenBrace = '{';
constexpr char kCloseBrace = '}';

/*
  Aug 2002. 
   Aesthetic issue. If an atom comes in with [] around it, then it is assumed
   to have an explicitly known H count - even is the H count is quite normal.
*/

static int unset_implicit_hydrogens_known_if_possible = 0;

void
set_unset_implicit_hydrogens_known_if_possible (int s)
{
  unset_implicit_hydrogens_known_if_possible = s;
}

/*
  Of the case where the Hcount is just wrong

  C[C]C fix to CCC
*/

static int unset_all_implicit_hydrogens_known_attributes = 0;

void
set_unset_all_implicit_hydrogens_known_attributes (int s)
{
  unset_all_implicit_hydrogens_known_attributes = s;
}

static IWString identifier_dataitem ("PCN<");

void
set_tdt_identifier_dataitem (const const_IWSubstring & d)
{
  identifier_dataitem = d;

  if (! identifier_dataitem.ends_with('<'))
    identifier_dataitem += '<';

  return;
}

/*
  We have the ability to append one or more dataitems to the name field
*/

static resizable_array_p<IWString> dataitems_to_append;

void
set_tdt_append_dataitem (const const_IWSubstring & a)
{
  IWString * tmp = new IWString(a);

  if (! tmp->ends_with('<'))
    *tmp += '<';

  dataitems_to_append.add(tmp);

  return;
}

/*
  When we are concatenating things onto the name, we can add either the
  contents of the tdt record, or the whole record itself
*/

static int append_dataitem_content = 1;

void
set_tdt_append_dataitem_content (int s)
{
  append_dataitem_content = s;
}

static int file_scope_display_smiles_interpretation_error_messages = 1;

void
set_display_smiles_interpretation_error_messages(int s)
{
  file_scope_display_smiles_interpretation_error_messages = s;
}

int
display_smiles_interpretation_error_messages()
{
  return file_scope_display_smiles_interpretation_error_messages;
}

namespace smiles {
// The separator between smiles and id.
static char output_separator = ' ';

void
set_smiles_output_separator(char c) {
  output_separator = c;
}

}  // namespace smiles

static int
check_for_append (const const_IWSubstring & buffer,    // as read in
                  IWString & add_to_name)              // what will finally be appended
{
  int na = dataitems_to_append.number_elements();

  for (int i = 0; i < na; i++)
  {
    if (buffer.starts_with(*(dataitems_to_append[i])))
    {
      if (add_to_name.nchars())
        add_to_name += ' ';

      if (append_dataitem_content)    // need to remove tag and <>
      {
        add_to_name += buffer.substr(dataitems_to_append[i]->length());
        add_to_name.chop();    // remove trailing >
      }
      else
        add_to_name += buffer;
    }
  }

  return 1;
}

int
Molecule::write_molecule_usmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << unique_smiles();

  if (comment.length())
    os << smiles::output_separator << comment;
  
  os << newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

int
Molecule::write_molecule_nausmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << non_aromatic_unique_smiles();

  if (comment.length())
    os << smiles::output_separator << comment;
  
  os << newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

int
Molecule::write_molecule_smi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << smiles();

  if (comment.length())
    os << smiles::output_separator << comment;
  
  os << newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

int
Molecule::write_molecule_rsmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << random_smiles();

  if (comment.length())
    os << smiles::output_separator << comment;
  
  os << newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

// As part of smiles parsing a string is available as the molecule
// name. If parsing of Chemaxon extensions is enabled, examine the
// string to see if the first token can be consumed as Chemaxon
// directives.
// Arg `processing_quoted_smiles` is passed to MaybeParseAsChemaxonExtension.
int
Molecule::SmilesSetName(const char * s, int nchars,
                        int processing_quoted_smiles) {
  if (nchars == 0) {
    return 1;
  }

  const_IWSubstring name(s, nchars);
  // cerr << "Processing name '" << name << "'\n";
  name.strip_leading_blanks();
  name.strip_trailing_blanks();
  if (name.empty()) {
    _molecule_name.resize(0);
    return 1;
  }

  if (! smiles::DiscernChemaxonSmilesExtensions() ||
      ! name.starts_with('|')) {
    set_name(name);
    return 1;
  }

  // First token of name might be a Chemaxon extension
  return MaybeParseAsChemaxonExtension(name, processing_quoted_smiles);
}

// We are looking to see if `name` matches |...|.
// Return the index of the closing |.
int
ClosingVBar(const const_IWSubstring& name) {
  for (int i = 1; i < name.length(); ++i) {
    if (isspace(name[i])) {  // Chemaxon patterns do not include spaces it seems...
      return -1;
    }
    if (name[i] == '|') {
      return i;
    }
  }

  return -1;
}

// `name` starts with | and perhaps is a Chemaxon smiles extension
// of the form |chemaxon| actual_name
// If we can extract a valid Chemaxon directive, that is removed
// from `name`. The molecule name is set.
// If `processing_quoted_smiles` is set, that means the input would
// have been
// "smiles |chemaxon|" name
// so we need to look for the closing quote.

int
Molecule::MaybeParseAsChemaxonExtension(const_IWSubstring& name,
                        int processing_quoted_smiles) {
  assert(name.starts_with('|'));

  int found_closing_vbar = ClosingVBar(name);

  // No closing vertical bar is fine, name might be '|foo bar'. `name` is unchanged.
  // Or should this be an error??
  if (found_closing_vbar < 0) {
    set_name(name);
    return 1;
  }

  // But if we have a closing vertical bar, then what is in there must be a valid
  // Chemaxon directive. Disallow empty ||
  if (found_closing_vbar == 1) {
    cerr << "Molecule::MaybeParseAsChemaxonExtension:empty Chemaxon directive '" << name << "'\n";
    return 0;
  }

  const_IWSubstring chemaxon(name.rawchars() + 1, found_closing_vbar - 1);

  name.remove_leading_chars(found_closing_vbar + 1);
  if (processing_quoted_smiles) {
    if (name[0] != '"') {
      cerr << "Molecule::MaybeParseAsChemaxonExtension:no closing quote '" << name << "'\n";
      return 0;
    }
    name += 1;
  }

  if (! ParseChemaxonExtension(chemaxon)) {
    cerr << "Molecule::MaybeParseAsChemaxonExtension:invalid Chemaxon directive '" << name << "'\n";
    return 0;
  }

  name.strip_leading_blanks();
  set_name(name);
  return 1;
}

int
Molecule::ParseChemaxonExtension(const const_IWSubstring& chemaxon) {
  // Implement sometime...
  cerr << "Molecule::ParseChemaxonExtension:ignoring " << chemaxon << '\n';

  std::unique_ptr<int[]> claimed(new_int(chemaxon.length()));
  if (! ParseCoords(chemaxon, claimed.get())) {
    cerr << "Molecule::ParseChemaxonExtension:cannot process coordinates '" << chemaxon << "'\n";
    return 0;
  }

  if (! ParseSpecialAtoms(chemaxon, claimed.get())) {
    cerr << "Molecule::ParseChemaxonExtension:cannot process special atoms '" << chemaxon << "'\n";
    return 0;
  }
  return 1;
}

int
FirstUnclaimed(const const_IWSubstring& chemaxon,
               int istart,
               int * claimed,
               char to_find) {
  const int nchars = chemaxon.length();
  for (int i = istart; i < nchars; ++i) {
    if (claimed[i]) {
      continue;
    }
    if (chemaxon[i] == to_find) {
      return i;
    }
  }

  return -1;
}

// Examine the unclaimed characters in `chemaxon` and look for
// the pattern <open_group>..unclaimed...<close_group>
// This is used for finding things like
//  (...)
//  $...$
// In the Chemaxon smiles extensions.
// If no such pattern is found, return -1, -1, 1
// If the open group but no close group is found, return -1, -1, 0
// If both ends are found, return (open, close, 1).
std::tuple<int, int, int>
OpenAndCloseGroup(const const_IWSubstring& chemaxon,
                  char open_group,
                  char close_group,
                  int * claimed) {
  const int open_paren = FirstUnclaimed(chemaxon, 0, claimed, open_group);
  if (open_paren < 0) {
    return std::make_tuple(-1, -1, 1);  // No result, but OK.
  }

  const int close_paren = FirstUnclaimed(chemaxon, open_paren + 1, claimed, close_group);
  if (close_paren < 0) {
    cerr << "OpenAndCloseGroup:no close paren '" << chemaxon << "'\n";
    return std::make_tuple(-1, -1, 0);   // No result, fatal.
  }

  for (int i = open_paren; i <= close_paren; ++i) {
    claimed[i] = 1;
  }

  return std::make_tuple(open_paren, close_paren, 1);  // All good.
}

int
Molecule::ParseCoords(const const_IWSubstring& chemaxon, int * claimed) {
  const auto [open_paren, close_paren, rc] = OpenAndCloseGroup(chemaxon, kOparen, kCparen, claimed);
  if (open_paren < 0 || close_paren < 0) {
    return rc;
  }

  // const_IWSubstring coords(chemaxon.data() + open_paren, close_paren - open_paren + 1);  with open and close parens
  const_IWSubstring coords(chemaxon.data() + open_paren + 1, close_paren - open_paren - 1);
  const_IWSubstring atom_token;
  const_IWSubstring coord_token;
  resizable_array<coord_t> values;
  for (int i = 0, atom_number = 0; coords.nextword(atom_token, i, ';'); ++atom_number) {
    if (atom_number >= _number_elements) {
      cerr << "Molecule::ParseCoords:Too many atoms\n";
      return 0;
    }
    values.resize_keep_storage(0);
    for (int j = 0, xyz = 0; atom_token.nextword(coord_token, j, ','); ++xyz) {
      coord_t value;
      if (! coord_token.numeric_value(value)) {
        cerr << "Molecule::ParseCoords:invalid numeric '" << coord_token << "'\n";
        return 0;
      }
      values << value;
    }
    _things[atom_number]->Setxyz(values);
  }


  return 1;
}

// 
const Element *
GetElement(const char * symbol) {
  const Element * e = get_element_from_symbol_no_case_conversion(symbol);
  if (e != nullptr) {
    return e;
  }
  // Need to create it.
//const auto esave = auto_create_new_elements();
//set_auto_create_new_elements(1);
  e = create_element_with_symbol(symbol);
//set_auto_create_new_elements(esave);

  return e;
}

int
Molecule::ParseSpecialAtoms(const const_IWSubstring& chemaxon, int * claimed) {
  constexpr char dollar = '$';
  const auto [first_dollar, second_dollar, rc] = OpenAndCloseGroup(chemaxon, dollar, dollar, claimed);
  if (first_dollar < 0 || second_dollar < 0) {
    return rc;
  }

  const_IWSubstring special_atoms(chemaxon.data() + first_dollar + 1, second_dollar - first_dollar - 1);

  IWString token;
  for (int i = 0, atom_number = 0; special_atoms.nextword_single_delimiter(token, i, ';'); ++atom_number) {
    if (token.empty()) {
      if (_things[atom_number]->atomic_symbol() == "*") {
        _things[atom_number]->set_element(GetElement("A"));
      }
    } else if (token == "Q_e") {
      _things[atom_number]->set_element(GetElement("Q"));
    } else if (token == "star_e") {
      _things[atom_number]->set_element(GetElement("*"));
    } else if (token.ends_with("_p")) {
      token.chop(2);
      _things[atom_number]->set_element(GetElement(token.null_terminated_chars()));
    }
  }

  return 1;
}

/*
  This is used to keep track of ring openings and closings during 
  reading a smiles
*/

class Smiles_Ring_Status: public resizable_array<atom_number_t>
{
  private:
    resizable_array<bond_type_t> _bt;

    int _rings_encountered;

  public:
    Smiles_Ring_Status();

    int complete() const;
    int report_hanging_ring_closures(std::ostream &) const;

    int encounter(int, atom_number_t, atom_number_t &, bond_type_t &);

    int rings_encountered() const { return _rings_encountered;}
};

Smiles_Ring_Status::Smiles_Ring_Status()
{
  _rings_encountered = 0;
}

/*
  We will be complete when we are not storing any valid atom numbers
*/

int
Smiles_Ring_Status::complete() const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (INVALID_ATOM_NUMBER != _things[i])
      return 0;
  }

  return 1;
}

int
Smiles_Ring_Status::report_hanging_ring_closures(std::ostream & os) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (INVALID_ATOM_NUMBER != _things[i])
      os << "Hanging ring closure " << i << endl;
  }

  return 1;
}

/*
  For efficiency, we temporarily store any directionality of a
  bond in the bond type itself - this avoids having a 2nd stack
  of bond directionalities.

  Make sure any bits used don't collide with any already used by btype
  Both UP and DOWN use bit 1024 so we can check for non-directionality
  by making just one comparison
*/

#define SMI_DIRECTIONAL_UP 2048
#define SMI_DIRECTIONAL_DOWN 4096
#define SMI_DIRECTIONAL_EITHER (SMI_DIRECTIONAL_UP | SMI_DIRECTIONAL_DOWN)

/*
  At atom A we have encountered a ring opening/closing number.
  Analyse that ring number.
  If ring_number corresponds to closing a ring, we set OTHER_END and
  return 1.
  For ring openings we return 0;

  Note that we aren't using the bt variable. For this reason, we do not
  correctly process the smiles C=1CC1.
  Change this sometime if it ever becomes a problem
*/

int
Smiles_Ring_Status::encounter(int ring_number, atom_number_t a,
                              atom_number_t & other_end,
                              bond_type_t & bt)
{
  assert(ring_number >= 0);

  if (ring_number >= _number_elements)     // definitely a new ring.
  {
    extend(2 * ring_number + 1, INVALID_ATOM_NUMBER);
    _bt.extend(2 * ring_number + 1, SINGLE_BOND);
  }

// If ring_number is in range, and _things[ring_number] is an OK atom
//  number, then we have a ring closure

  else if (INVALID_ATOM_NUMBER != _things[ring_number])
  {
//  cerr << "Encounter ring " << ring_number << " making bond between atom " << _things[ring_number] << " type " << bt << " stored type " << _bt[ring_number] << endl;
    other_end = _things[ring_number];
    _things[ring_number] = INVALID_ATOM_NUMBER;

    if (NOT_A_BOND == _bt[ring_number])    // was not set during ring opening, use what just came in
    {
      if (NOT_A_BOND == bt)
        _bt[ring_number] = SINGLE_BOND | AROMATIC_BOND;
      else 
        _bt[ring_number] = bt;
    }
    else if (NOT_A_BOND == bt)      // was set during opening, and closing does not specify
      bt = _bt[ring_number];
    else if (bt == _bt[ring_number])    // great, same as on ring opening
      ;
    else if (SMI_DIRECTIONAL_UP & bt)
      _bt[ring_number] = bt ^ SMI_DIRECTIONAL_UP;
    else if (SMI_DIRECTIONAL_DOWN & bt)
      _bt[ring_number] = bt ^ SMI_DIRECTIONAL_DOWN;
    else
    {
      cerr << "Smiles_Ring_Status::encounter:inconsistent bonds at opening " << _bt[ring_number] << " and closing " << _bt[ring_number] << " of ring " << ring_number << " using last\n";
      _bt[ring_number] = bt;
    }

    bt = _bt[ring_number];
//  cerr << "Returned bond type " << bt << endl;
    return 1;    // ring closure
  }

// If we come to here, then this is a ring opening. Note that if we
// have a directional bond for a ring opening, we need to swap the 
// directionality because the closing will look down from the opposite
// direction

  if (bt & SMI_DIRECTIONAL_UP)
    bt = (bt ^ SMI_DIRECTIONAL_UP) | SMI_DIRECTIONAL_DOWN;
  else if (bt & SMI_DIRECTIONAL_DOWN)
    bt = (bt ^ SMI_DIRECTIONAL_DOWN) | SMI_DIRECTIONAL_UP;

  _things[ring_number] = a;
  _bt[ring_number] = bt;

//cerr << "_BT " << ring_number << ' ' << bt << endl;

  _rings_encountered++;

  return 0;      // ring opening
}

/*
  A hcount consists of an H followed by an optional digit. 
*/

static int
process_hcount(const char * smiles, int & hcount, int & nchars)
{
  hcount = 1;
  nchars = 1;
  smiles++;
  if (isdigit(*smiles))
  {
    hcount = *smiles - '0';
    nchars++;
  }

  return 1;
}

/*
  A charge specifier consists of a sequence of '-' or '+', OR
  a '-' or '+' followed by a digit.
*/

static int
process_charge_specifier(const char * smiles, int sign, formal_charge_t & fc, 
                         int & nchars)
{
  nchars = 1;
  smiles++;

  char c;
  if (1 == sign)
  {
    c = '+';
    fc = 1;
  }
  else if (-1 == sign)
  {
    c = '-';
    fc = -1;
  }
  else
  {
    cerr << "What is this!!!! " << sign << endl;
    return 0;
  }

  int count = 1;      // a count of the number of '-' or '+'
  while (c == *smiles)
  {
    count++;
    fc += sign;
    smiles++;
    nchars++;
  }

  if (count > 1)
  {
    if (reasonable_formal_charge_value(fc))
      ;
    else if (file_scope_display_smiles_interpretation_error_messages)
      cerr << "process_charge_specifier: possibly unreasonable formal charge " << fc << endl;

    return 1;
  }

// There was only one '-' or '+'. Do we now have a digit?

  if (isdigit(*smiles))
  {
    int tmp = *smiles - '0';
    fc *= tmp;
    nchars++;

    if (reasonable_formal_charge_value(fc))
      ;
    else if (file_scope_display_smiles_interpretation_error_messages)
      cerr << "process_charge_specifier: possibly unreasonable formal charge " << fc << endl;
  }

  return 1;
}

static int
process_chirality_specifier(const char * smiles, int & chiral_count, int & nchars)
{
  nchars = 1;
  chiral_count = 1;

  smiles++;
  if ('@' == *smiles)
  {
    nchars++;
    chiral_count++;
  }

  return 1;
}

static const Element *
fetch_or_create_R_element(const IWString & r)
{
  const Element * e = get_element_from_symbol_no_case_conversion(r);
  if (nullptr != e)
    return e;

  return create_element_with_symbol(r);
}

/*
  For the Organic Subset, we need rapid access to certain elements, so
  we invent this dummy class to enable us to initialise the element pointers
*/

static const Element * smi_element_star = nullptr;   // initialised for check in parse_smiles_token
static const Element * smi_element_b;
static const Element * smi_element_c;
static const Element * smi_element_n;
static const Element * smi_element_o;
static const Element * smi_element_f;
static const Element * smi_element_p;
static const Element * smi_element_s;
static const Element * smi_element_cl;
static const Element * smi_element_br;
static const Element * smi_element_i;
static const Element * smi_element_hydrogen;
static const Element * smi_element_a = nullptr;

static void
initialise_organic_subset()
{
  smi_element_star     = get_element_from_atomic_number(0);
  smi_element_hydrogen = get_element_from_atomic_number(1);
  smi_element_b    = get_element_from_atomic_number(5);
  smi_element_c    = get_element_from_atomic_number(6);
  smi_element_n    = get_element_from_atomic_number(7);
  smi_element_o    = get_element_from_atomic_number(8);
  smi_element_f    = get_element_from_atomic_number(9);
  smi_element_p    = get_element_from_atomic_number(15);
  smi_element_s    = get_element_from_atomic_number(16);
  smi_element_cl   = get_element_from_atomic_number(17);
  smi_element_br   = get_element_from_atomic_number(35);
  smi_element_i    = get_element_from_atomic_number(53);

  return;
}

/*
  Fundamental routine for parsing tokens in smiles.
  We return the number of characters we fully parse, all other
  info passed by reference.

  The caller must ensure that the token does NOT start with a '['
*/

int
parse_smiles_token(const char * smiles,
                   int characters_to_process,
                   const Element * &    e,
                   int & aromatic)
{
  int c = *smiles;

  aromatic = AROMATICITY_NOT_DETERMINED;

  if ('C' == *smiles)
  {
    if (characters_to_process > 1 && 'l' == smiles[1])
    {
      e = smi_element_cl;
      return 2;
    }

    e = smi_element_c;
    return 1;
  }

  if ('N' == c)
  {
    e = smi_element_n;
    return 1;
  }

  if ('O' == c)
  {
    e = smi_element_o;
    return 1;
  }

  if ('S' == c)
  {
    e = smi_element_s;
    return 1;
  }

  if ('F' == c)     // Aromatic F not recognised
  {
    e = smi_element_f;
    return 1;
  }

  if (characters_to_process > 1 && 'B' == c && 'r' == smiles[1])
  {
    e = smi_element_br;
    return 2;
  }

  if ('I' == c)    // Aromatic I not recognised
  {
    e = smi_element_i;
    return 1;
  }

  if ('P' == c)
  {
    e = smi_element_p;
    return 1;
  }

  if ('B' == c)
  {
    e = smi_element_b;
    return 1;
  }

  if ('c' == c)
  {
    e = smi_element_c;
    aromatic = AROMATIC;
    return 1;
  }

  if ('n' == c)
  {
    e = smi_element_n;
    aromatic = AROMATIC;
    return 1;
  }

  if ('o' == c)
  {
    e = smi_element_o;
    aromatic = AROMATIC;
    return 1;
  }

  if ('s' == c)
  {
    e = smi_element_s;
    aromatic = AROMATIC;
    return 1;
  }

  if ('p' == c)
  {
    e = smi_element_p;
    aromatic = AROMATIC;
    return 1;
  }

  if ('b' == c)     // aromatic Boron???
  {
    e = smi_element_b;
    aromatic = AROMATIC;
    return 1;
  }

  if ('H' == c)
  {
    e = smi_element_hydrogen;
    return 1;
  }

  if ('*' == c)
  {
    e = smi_element_star;
    return 1;
  }

  if ('a' == c)
  {
    if (nullptr == smi_element_a)
    {
      smi_element_a = get_element_from_symbol_no_case_conversion("a");
//    cerr << "smi_element_a now " << smi_element_a << endl;

      if (nullptr != smi_element_a)
        ;
      else if (! auto_create_new_elements())
      {
        cerr << "Cannot create 'a' element, use '-E autocreate'\n";
        return 0;
      }
      else
        smi_element_a = create_element_with_symbol('a');
    }

    e = smi_element_a;

//  cerr << "'a' atom found, arom " << e->permanent_aromatic() << endl;

    if (e->permanent_aromatic())
      aromatic = AROMATIC;

    return 1;
  }

  if (islower(c))    // some kind of aromatic other than anything recognised
  {
    e = get_element_from_symbol_no_case_conversion(smiles, 1);

    if (nullptr == e && auto_create_new_elements())
      e = create_element_with_symbol(smiles[0]);
      
    if (nullptr == e)
    {
      if (file_scope_display_smiles_interpretation_error_messages)
        cerr << "parse_smiles_token:invalid element specification '" << static_cast<char>(c) << "'\n";
      return 0;
    }

    if (e->permanent_aromatic())
      aromatic = AROMATIC;

    return 1;
  }

  if (file_scope_display_smiles_interpretation_error_messages)
  {
    cerr << "parse_smiles_token:cannot interpret smiles symbol '";
    cerr.write(smiles, characters_to_process);
    cerr << "'\n";
  }

  return 0;
}

static int
maybe_deuterium_or_tritium(const char * smiles,
                           const Element * & e,
                           int & atomic_mass)
{
  if ('D' == *smiles && element::interpret_d_as_deuterium())
    atomic_mass = 2;
  else if ('T' == *smiles && element::interpret_t_as_tritium())
    atomic_mass = 3;
  else
    return 0;

  e = get_element_from_atomic_number(1);

  return 1;
}

/*
  Parsing an atom enclosed in a square bracket
*/

//#define DEBUG_PARSE_SMILES_TOKEN

int
parse_smiles_token (const char * smiles,
                    int characters_to_process,
                    const Element * &    e,
                    aromaticity_type_t & aromatic,
                    formal_charge_t &    fc,
                    int             &    hcount,
                    int             &    chiral_count,
                    int             &    atomic_mass,
                    int             &    atom_map)
{
#ifdef DEBUG_PARSE_SMILES_TOKEN
  cerr << "parse_smiles_token ";

  for (int i = 0; i < characters_to_process; ++i)
  {
    cerr << smiles[i];
  }
  cerr << endl;
#endif

  aromatic = AROMATICITY_NOT_DETERMINED;
  fc = 0;
  hcount = 0;
  chiral_count = 0;
  atomic_mass = 0;
  atom_map = 0;

  int rc = 1;      //  1 char is the open square bracket.
  const char * original_smiles_ptr = smiles;
  smiles++;

// Is there an isotope specification

  while (isdigit(*smiles))
  {
    atomic_mass = 10 * atomic_mass + (*smiles - '0');
    rc++;
    smiles++;
  }

// Then we should get the element

  int tmp;    // number characters consumed by the element

  // cerr << "Begin check on smiles at " << *smiles << '\n';
  if (isalpha(*smiles))
  {
    if (atomic_symbols_can_have_arbitrary_length())     // note that if we have D and/or T in this case, we will not pick up the mass diff.... Fix if ever anyone cares...
      tmp = element_from_long_smiles_string(smiles, characters_to_process, e);
    else  
      tmp = element_from_smiles_string(smiles, characters_to_process, e);

    // cerr << "atomic_symbols_can_have_arbitrary_length " << atomic_symbols_can_have_arbitrary_length << " tmp " << tmp << '\n';

    if (0 == tmp && characters_to_process > 1 && ']' == smiles[1])     
      tmp = maybe_deuterium_or_tritium(smiles, e, atomic_mass);

    if (0 == tmp)
    {
      if (file_scope_display_smiles_interpretation_error_messages)
        cerr << "Cannot determine element '" << original_smiles_ptr << "'\n";
      return 0;
    }

    if (isupper(*smiles))   // definitely not aromatic
      ;
    else if (e->organic())
      aromatic = AROMATIC;
    else if (e->permanent_aromatic())
      aromatic = AROMATIC;
  }
  else if ('*' == *smiles)
  {
    e = smi_element_star;
    tmp = 1;
  }
  else if ('#' == *smiles && characters_to_process > 1 && isdigit(smiles[1]))
  {
    int z;
    tmp = fetch_numeric(smiles + 1, z);
//  cerr << "Got numeric atomic number, tmp = " << tmp << " z = " << z << endl;
    e = get_element_from_atomic_number(z);
    if (nullptr == e)
    {
      if (file_scope_display_smiles_interpretation_error_messages)
        cerr << "Unrecognised atomic number " << z << endl;
      return 0;
    }
    tmp++;     // don't forget the #
  }
  else
  {
    if (file_scope_display_smiles_interpretation_error_messages)
      cerr << "parse_smiles_token:must be alphabetic '" << *smiles << "'\n";
    return 0;
  }

  rc += tmp;
  smiles += tmp;

// There are a number of things which can follow the atomic symbol.
// H count, charges, chirality

// We have to make sure that we don't get some '+' and some '-', or more that
// two chirality things, or....

  int plus_encountered = 0;
  int minus_encountered = 0;
  int hcount_encountered = 0;
  int chiral_encountered = 0;

  while (']' != *smiles && rc < characters_to_process)
  {
    // cerr << "processing smiles char " << *smiles << " in " << smiles << '\n';
    if ('H' == *smiles)
    {
      if (hcount_encountered)
      {
        if (file_scope_display_smiles_interpretation_error_messages)
          cerr << "Duplicate H count specification '" << original_smiles_ptr << "'\n";
        return 0;
      }
      int tmp = 0;
      if (! process_hcount(smiles, hcount, tmp))
      {
        cerr << "Cannot process h count specification '" << original_smiles_ptr << "'\n";
        return 0;
      }
      rc += tmp;
      smiles += tmp;

      hcount_encountered++;
    }
    else if ('+' == *smiles)
    {
      if (plus_encountered)
      {
        cerr << "Duplicate '+' specification '" << original_smiles_ptr << "'\n";
        return 0;
      }
      int tmp = 0;
      if (! process_charge_specifier(smiles, 1, fc, tmp))
      {
        cerr << "Cannot process '+' specifier '" << original_smiles_ptr << "'\n";
        return 0;
      }

      rc += tmp;
      smiles += tmp;

      plus_encountered++;
    }
    else if ('-' == *smiles)
    {
      if (minus_encountered)
      {
        cerr << "Duplicate '+' specification '" << original_smiles_ptr << "'\n";
        return 0;
      }
      int tmp = 0;
      if (! process_charge_specifier(smiles, -1, fc, tmp))
      {
        cerr << "Cannot process '-' specifier '" << original_smiles_ptr << "'\n";
        return 0;
      }

      rc += tmp;
      smiles += tmp;

      minus_encountered++;
    }
    else if ('@' == *smiles)
    {
      if (chiral_encountered)
      {
        cerr << "Duplicate chirality specifier '" << original_smiles_ptr << "'\n";
        return 0;
      }

      int tmp = 0;
      if (! process_chirality_specifier(smiles, chiral_count, tmp))
      {
        cerr << "Cannot process chirality specifier '" << original_smiles_ptr << "'\n";
        return 0;
      }

      rc += tmp;
      smiles += tmp;

      chiral_encountered++;
    }
    else if (isdigit(*smiles) && 2 == rc && "R" == e->symbol())
    {
      int r = 0;
      while (isdigit(*smiles) && rc < characters_to_process)
      {
        r = 10 * r + (*smiles) - '0';
        rc++;
        smiles++;
      }

      IWString tmp;
      tmp << "R" << r;
      e = fetch_or_create_R_element(tmp);
      if (nullptr == e)
      {
        cerr << "Cannot build R# from smiles '" << tmp << "'\n";
        return 0;
      }
    }
    else if ('#' == *smiles && characters_to_process > 1 && isdigit(smiles[1]))
    {
      int z;
      int tmp = fetch_numeric(smiles + 1, z);
      e = get_element_from_atomic_number(z);
      if (nullptr == e)
      {
        if (file_scope_display_smiles_interpretation_error_messages)
          cerr << "Unrecognised atomic number " << z << endl;
        return 0;
      }
      smiles += (tmp + 1);
      rc += (tmp + 1);
    }
    else if (':' == *smiles && characters_to_process > 1 && isdigit(smiles[1]))    // an atom number specification
    {
      int tmp = fetch_numeric_char(smiles+1, atom_map, characters_to_process);
      smiles += (tmp + 1);
      rc += (tmp + 1);
    }
    else
    {
      if (file_scope_display_smiles_interpretation_error_messages)
        cerr << "Unrecognised smiles symbol '" << original_smiles_ptr << "', '" << *smiles << "'\n";
      return 0;
    }
  }

  if (chiral_encountered && moleculeio::ignore_all_chiral_information_on_input())
    chiral_encountered = 0;

  rc++;     // we got our closing bracket.
  return rc;
}

/*
*/

int
Molecule::read_molecule_smi_ds (iwstring_data_source & input)
{
  assert(input.good());

  if (input.eof())
    return 0;

  IWString buffer;

  EXTRA_STRING_RECORD(input, buffer, "read mol smi");

  if (buffer.length() < 1)
    return 0;

  if (! build_from_smiles(buffer.rawchars(), buffer.length()))
  {
    cerr << "Molecule::read_molecule_smi_ds: Cannot interpret smiles\n";
    cerr << buffer << endl;
    return 0;
  }

  return 1;
}

static int ignore_tdts_with_no_smiles = 0;

void
set_ignore_tdts_with_no_smiles (int s)
{
  ignore_tdts_with_no_smiles = s;
}

static IWString smiles_tag = "$SMI<";

void
set_smiles_tag (const const_IWSubstring & tag)
{
  smiles_tag = tag;

  if (! smiles_tag.ends_with('<'))
    smiles_tag += '<';

  return;
}

int
Molecule::read_molecule_tdt_ds (iwstring_data_source & input)
{
  assert(input.good());

  if (input.eof())
    return 0;

  input.set_skip_blank_lines(1);

  const_IWSubstring buffer;

  int got_structure = 0;
  int got_identifier = 0;

  IWString add_to_name;

  while (1)
  {
    EXTRA_STRING_RECORD(input, buffer, "read mol tdt");

    if (buffer.starts_with(smiles_tag) && 0 == got_structure)    // we only parse the first smiles in a TDT
    {
      buffer += smiles_tag.length();     // skip over '$SMI<'

      int close_angle_bracket = buffer.index('>');
      assert(close_angle_bracket >= 1);

      buffer.chop();
      if (! build_from_smiles(buffer.rawchars(), close_angle_bracket))
      {
        cerr << "Molecule::read_molecule_tdt_ds: Cannot interpret smiles\n";
        return 0;
      }

      got_structure = 1;

//    If this is a tdt in dump form, look for a PCN<> attribute.

      int i = buffer.find(">PCN<");
      if (i < 0)
        continue;

      const_IWSubstring zrest = buffer.substr(i + 5);

//    cerr << "The rest is '" << zrest << "'\n";

      if (zrest.nchars())
      {
        int j = zrest.find('>');
        zrest.iwtruncate(j);
        set_name(zrest);
      }

      return 1;     // molecule and name is all we can do in this case
    }

    if (buffer.starts_with(identifier_dataitem) && 0 == got_identifier) // PCN follows $SMI
    {
      IWString tmp = substr(buffer, identifier_dataitem.nchars());
      tmp.chop();
      set_name(tmp);

      got_identifier = 1;

      continue;
    }

    if ('|' == buffer && got_structure)
    {
      if (add_to_name.length())
        _molecule_name << ' ' << add_to_name;

      return 1;
    }

    if ('|' == buffer)
    {
      cerr << "Molecule::read_molecule_tdt_ds: no structure in TDT\n";
      if (0 == ignore_tdts_with_no_smiles)
        return 0;
      if (got_structure || got_identifier)
        return 0;

      continue;       // let's skip this TDT and look at the next one
    }

    if (moleculeio::read_extra_text_info())
    {
      IWString * tmp = new IWString(buffer);
      _text_info.add(tmp);
    }

    if (got_structure)
      check_for_append(buffer, add_to_name);
  }
}

int
ParseBoxedCoordinates(const_IWSubstring c,  // Local copy.
                      Coordinates& coords) {
  assert(c.starts_with('B'));
  c += 1;  // SKip over 'B'.
  coordinate_box::LayerPosition layer_position;
  const int consumed = coordinate_box::FromString(c, layer_position);
  if (consumed != c.length()) {
    return 0;
  }

  coordinate_box::ConcentricBox box;

  coords = box.CellToCoordinates<float>(layer_position);

  return c.length() + 1;
}

static int
parse_coordinates(const char * smiles,
                  int characters_to_process,
                  int characters_processed,
                  Coordinates & coords)
{
  smiles += characters_processed;

  assert(kOpenBrace == smiles[0] && kOpenBrace == smiles[1]);

  smiles += 2;

  characters_to_process -= 2;

// locate the closing braces

  int close_brace_pos = -1;
  for (int i = 0; i < characters_to_process - 1; i++)
  {
    if (kCloseBrace == smiles[i] && kCloseBrace == smiles[i + 1])
    {
      close_brace_pos = i;
      break;
    }
  }

  if (close_brace_pos <= 0)
  {
    cerr << "parse_coordinates: no closing braces found, or empty braces\n";
    return 0;
  }

  const_IWSubstring c(smiles, close_brace_pos);

  if (c[0] == 'B') {
    if (! ParseBoxedCoordinates(c, coords)) {
      cerr << "parse_coordinates::Cannot parse boxed coordinates '" << c << "'\n";
      return 0;
    }
  } else if (! coords.read(c, ',')) {
    cerr << "parse_coordinates: cannot parse coordinates '" << c << "'\n";
    return 0;
  }

  return 2 + close_brace_pos + 2;
}

/*
  We have read a smiles and are about to start looking for a Kekule form.
  Make sure we can deal with things like

  n1cc[20c]n[20c]1

  which is invalid because there are no H atoms on the two 20c atoms
*/

int
Molecule::_do_unfix_implicit_hydrogens_on_aromatic_isotopic_atoms(const int * aromatic_atoms)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == aromatic_atoms[i])
      continue;

    Atom * ai = _things[i];

    if (0 == ai->isotope())
      continue;

    if (! ai->implicit_hydrogens_known())
      continue;

    ai->set_implicit_hydrogens_known(0);
    rc++;
  }

  return rc;
}

int
Molecule::_smiles_add_bond (atom_number_t previous_atom,
                            atom_number_t current_atom,
                            bond_type_t bt)
{
  bond_type_t type_to_add;
  if (PERMANENT_AROMATIC_BOND == bt)
    type_to_add = bt;
  else
    type_to_add = BOND_TYPE_ONLY(bt);

  if (! add_bond(previous_atom, current_atom, type_to_add, 1))    // last arg, 1, means partial molecule
    return 0;

  if (0 == (bt & SMI_DIRECTIONAL_EITHER))       // non directional
    return 1;

//if (! discern_cis_trans_bonds())    Nov 03. Always discern cis-trans stuff in smiles
//  return 1;

// Fetch the last bond added to previous_atom and set its directionality

  Bond * b = _things[previous_atom]->last_item();

  assert(previous_atom == b->a1());

  if (bt & SMI_DIRECTIONAL_UP)
    b->set_directional_up();
  else 
    b->set_directional_down();

//cerr << "Bond " << (*b) << " added, " << hex << bt << dec << endl;

  return 1;
}

/*static int
fetch_ring_number (const char * s,
                   int & ring_number,
                   int & characters_processed)
{
  if (isdigit (*s))
  {
    ring_number = *s - '0';
    characters_processed = 1;
    return 1;
  }

  if ('%' != *s)
    return 0;

  s++;    // skip over the %

  int tmp = fetch_numeric (s, ring_number, 2);
  if (2 != tmp)
    return 0;
  
  characters_processed = 3;

  return 1;
}*/

static int
fetch_ring_number (const char * s,
                   int nchars,
                   int & ring_number,
                   int & characters_processed)
{
  if (isdigit(*s))
  {
    ring_number = *s - '0';
    characters_processed = 1;
    return 1;
  }

#ifdef DEBUG_FETCH_RING_NUMBER
  cerr << "fetch_ring_number:looking at '" << *s << "', nchars " << nchars << endl;
  for (int i = 0; i < nchars; ++i)
  {
    cerr << s[i];
  }
  cerr << endl;
#endif

  if ('%' != *s)
    return 0;

  if (nchars < 3)     // shortest would be %nn
    return 0;

  s++;    // skip over the %

  if (kOparen != *s)    // must be a simple %nn
  {
    if (2 != fetch_numeric_char(s, ring_number, 2))
      return 0;

    characters_processed = 3;

    return 1;
  }

// The more complex case of %(nnn)

  nchars -= 2;      // we have consumed %(
  s++;

  ring_number = 0;

  for (int i = 0; i < nchars; i++)
  {
    int j = *s - '0';

    if (j >= 0 && j <= 9)
    {
      ring_number = 10 * ring_number + j;
      s++;
    }
    else if (kCparen == *s)
    {
      characters_processed = i + 3;
      return 1;
    }
    else
      return 0;
  }

  return 0;
}

/*
  These are set up so that OR will work;
*/

#define PREVIOUS_TOKEN_NOT_SPECIFIED 0
#define PREVIOUS_TOKEN_RING_SPECIFIER 1
#define PREVIOUS_TOKEN_BOND_SPECIFIER 2
#define PREVIOUS_TOKEN_ATOM_SPECIFIER 4
#define PREVIOUS_TOKEN_OPEN_PAREN 8
#define PREVIOUS_TOKEN_CLOSE_PAREN 16
#define PREVIOUS_TOKEN_DOT 32

/*
  New way - used just in smarts for now
*/

#define PREVIOUS_TOKEN_WAS_NOT_SPECIFIED 0
#define PREVIOUS_TOKEN_WAS_RING 1
#define PREVIOUS_TOKEN_WAS_BOND 2
#define PREVIOUS_TOKEN_WAS_ATOM 3
#define PREVIOUS_TOKEN_WAS_OPEN_PAREN 4
#define PREVIOUS_TOKEN_WAS_CLOSE_PAREN 5
#define PREVIOUS_TOKEN_WAS_THREE_DOTS 6
#define PREVIOUS_TOKEN_WAS_DOWN_THE_BOND 7
// Don't forget to update CTDIM

#define CTDIM 8

static int * compatability_table = nullptr;     // never freed

static int
initialise_compatability_table()
{
  compatability_table = new_int(CTDIM * CTDIM);

// What can begin a smarts

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_NOT_SPECIFIED + PREVIOUS_TOKEN_WAS_ATOM] = 1;

// What can follow a ring

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow a bond

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_DOWN_THE_BOND] = 1;

// What can follow an atom

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow an open paren

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow a close paren

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow three dots

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_THREE_DOTS + PREVIOUS_TOKEN_WAS_ATOM] = 1;

// what can follow down the bond

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_DOWN_THE_BOND + PREVIOUS_TOKEN_WAS_ATOM] = 1;

  return 1;
}

static int
check_compatiability_table(int & previous_token_was, int nt)
{
  if (0 == compatability_table[CTDIM * previous_token_was + nt])
    return 0;

  previous_token_was = nt;

  return 1;
}

static int
previous_token_should_be(int previous_token_was, int possible_values)
{
  if (previous_token_was & possible_values)
    return 1;

  return 0;
}

static int
valid_end_of_smiles_character(int last_token)
{
  if (PREVIOUS_TOKEN_WAS_BOND == last_token)
    return 0;

  return 1;
}

// #define DEBUG_BUILD_FROM_SMILES

int
Molecule::_build_from_smiles(const char * smiles,
                             int characters_to_process,
                             Smiles_Ring_Status & ring_status,
                             int * aromatic_atoms)
{
  // initialise elements first time through
  if (nullptr == smi_element_star) {
    initialise_organic_subset();
  }

  int characters_processed = 0;

// We need to be somewhat careful about resizing, as we may be called recursively

  if (_elements_allocated - _number_elements < characters_to_process) {
    resize(_elements_allocated + characters_to_process);
    if (_bond_list.elements_allocated() - _bond_list.number_elements() < characters_to_process)
      _bond_list.resize(_bond_list.elements_allocated() + characters_to_process);
  }

  // If we are processing a quoted string, remove the leading quote and note the condtion.
  int processing_quoted_smiles = 0;
  if (*smiles == '"' && smiles::ProcessQuotedSmiles())  {
    ++smiles;
    --characters_to_process;
    processing_quoted_smiles = 1;
  }


// Various stack to keep track of branches. Probably should be combined into
// a single kind with a new object.

  resizable_array<atom_number_t>  atom_stack;
  resizable_array<int>            chirality_stack;

  int previous_token_was = PREVIOUS_TOKEN_NOT_SPECIFIED;

  bond_type_t previous_bond_type = INVALID_BOND_TYPE;
  atom_number_t previous_atom    = INVALID_ATOM_NUMBER;
  int previous_atom_chiral_count = 0;

  int paren_level = 0;

  int fragments_found = 0;

  while (characters_processed < characters_to_process)
  {
    const char * s = smiles + characters_processed;

    if (characters_processed && (isdigit(*s) || '%' == *s)) {
      int ring_number;
      int nchars;
      if (! fetch_ring_number(s, characters_to_process - characters_processed, ring_number, nchars))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Invalid ring specification");
        return 0;
      }
//    cerr << "NCHARS " << nchars << " for ring " << ring_number << endl;

      bond_type_t bt;
      if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
        bt = previous_bond_type;
      else if (PREVIOUS_TOKEN_DOT == previous_token_was)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Ring opening cannot begin smiles");
        return 0;
      }
      else
        bt = NOT_A_BOND;

      atom_number_t other_end;    // will be set when a ring closure
      if (ring_status.encounter(ring_number, previous_atom, other_end, bt))   // ring closing
      {
#ifdef DEBUG_BUILD_FROM_SMILES
        cerr << "Ring closing for ring " << ring_number << ", other end is " << other_end << ", bond " << bt << endl;
#endif

        if (previous_atom == other_end)
        {
          smiles_error_message(smiles, characters_to_process, characters_processed, "consecutive ring opening and closing");
          return 0;
        }

        if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was &&    // deal with things like CCC=1CCCC-1N
            NOT_A_BOND != bt &&
            bt != previous_bond_type)
        {
          if ((SMI_DIRECTIONAL_DOWN    ^ bt) == previous_bond_type)   // OK to have directional bond one end only, never really tested directional bonds both ends...
            bt = previous_bond_type;
          else if ((SMI_DIRECTIONAL_UP ^ bt) == previous_bond_type)
            bt = previous_bond_type;
          else
          {
            if (file_scope_display_smiles_interpretation_error_messages)
              cerr << "Inconsistent ring opening (" << bt << ") and closing (" << previous_bond_type << ") bond types\n";
            return 0;
          }
        }

        if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
          bt = previous_bond_type;
        else if (NOT_A_BOND == bt)
          bt = SINGLE_BOND;

//      If we are closing a ring with a directional bond, we are looking
//      down the bond in the opposite direction from when it was initially placed

        if (! _smiles_add_bond(previous_atom, other_end, bt))
        {
          smiles_error_message(smiles, characters_to_process, characters_processed, "invalid bonding arrangement");
          return 0;
        }

        Chiral_Centre * c = chiral_centre_at_atom(other_end);
        if (c)
          c->got_ring_closure_bond(ring_number, _number_elements - 1);

        if (previous_atom_chiral_count)    // the previous atom was chiral
        {
          if (! _smi_atom_bonded_to_chiral_centre(previous_atom, previous_atom_chiral_count, other_end))
          {
            smiles_error_message(smiles, characters_to_process, characters_processed,
                          "Error processing atom attached to chiral atom");
            return 0;
          }
        }
      }
      else if (previous_atom_chiral_count)
      {
#ifdef DEBUG_BUILD_FROM_SMILES
        cerr << "Ring opening for ring " << ring_number << " previous atom chiral " << previous_atom << endl;
#endif
        Chiral_Centre * c = chiral_centre_at_atom(previous_atom);
        assert(c);
        c->got_ring_opening_bond(ring_number, previous_atom_chiral_count);
      }
#ifdef DEBUG_BUILD_FROM_SMILES
      else
      {
        cerr << "Ring opening, " << ring_number << " previous atom not chiral, bt " << bt << endl;
      }
#endif

      characters_processed += nchars;
      previous_token_was = PREVIOUS_TOKEN_RING_SPECIFIER;
    }
    else if (characters_processed && SINGLE_BOND_SYMBOL == *s)
    {
      if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Consecutive bond symbols");
        return 0;
      }
      if (PREVIOUS_TOKEN_DOT == previous_token_was || PREVIOUS_TOKEN_NOT_SPECIFIED)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Fragment cannot start with bond");
        return 0;
      }
      characters_processed++;
      previous_bond_type = SINGLE_BOND;
      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
    }
    else if (characters_processed && DOUBLE_BOND_SYMBOL == *s)
    {
      if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Consecutive bond symbols");
        return 0;
      }
      if (PREVIOUS_TOKEN_DOT == previous_token_was || PREVIOUS_TOKEN_NOT_SPECIFIED)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Fragment cannot start with bond");
        return 0;
      }
      characters_processed++;
      previous_bond_type = DOUBLE_BOND;
      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
    }
    else if (characters_processed && TRIPLE_BOND_SYMBOL == *s)
    {
      if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Consecutive bond symbols");
        return 0;
      }
      if (PREVIOUS_TOKEN_DOT == previous_token_was || PREVIOUS_TOKEN_NOT_SPECIFIED)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Fragment cannot start with bond");
        return 0;
      }
      characters_processed++;
      previous_bond_type = TRIPLE_BOND;
      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
    }
    else if (characters_processed && AROMATIC_BOND_SYMBOL == *s)
    {
      if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Consecutive bond symbols");
        return 0;
      }
      if (PREVIOUS_TOKEN_DOT == previous_token_was || PREVIOUS_TOKEN_NOT_SPECIFIED)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Fragment cannot start with bond");
        return 0;
      }
      characters_processed++;
//    previous_bond_type = AROMATIC_BOND;
      previous_bond_type = PERMANENT_AROMATIC_BOND;
      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
    }
    else if (' ' == *s || '\t' == *s || ',' == *s)
    {
      if (0 != paren_level) {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Un-closed parenthesis");
        return 0;
      }
      if ( (PREVIOUS_TOKEN_ATOM_SPECIFIER != previous_token_was) &&
           (PREVIOUS_TOKEN_RING_SPECIFIER != previous_token_was)) {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Improperly terminated smiles");
        return 0;
      }

      if (characters_to_process - characters_processed > 1) {
        if (! SmilesSetName(s + 1, characters_to_process - characters_processed - 1, processing_quoted_smiles)) {
          return 0;
        }

        return 1;
      }
    }
    else if (kOparen == *s)
    {
      if (! previous_token_should_be(previous_token_was, 
                     PREVIOUS_TOKEN_ATOM_SPECIFIER | PREVIOUS_TOKEN_RING_SPECIFIER))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, 
                              "Branch specifier can only follow an atom or ring specifier");
        return 0;
      }

//    Push the various stacks.

      paren_level++;

      atom_stack.add(previous_atom);
      chirality_stack.add(previous_atom_chiral_count);

      characters_processed++;
      previous_token_was = PREVIOUS_TOKEN_OPEN_PAREN;
    }
    else if (characters_processed && kCparen == *s)
    {
      if (! previous_token_should_be(previous_token_was,
                PREVIOUS_TOKEN_ATOM_SPECIFIER | PREVIOUS_TOKEN_RING_SPECIFIER))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Close paren can only follow an atom");
        return 0;
      }
      if (atom_stack.empty())
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Parenthesis mismatch");
        return 0;
      }

      if (0 == paren_level)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Too many close paren's");
        return 0;
      }

      paren_level--;

      previous_atom = atom_stack.pop();
      previous_token_was = PREVIOUS_TOKEN_ATOM_SPECIFIER;
      previous_atom_chiral_count = chirality_stack.pop();
      characters_processed++;
#ifdef DEBUG_BUILD_FROM_SMILES
      cerr << "Closing paren found, previous atom now " << previous_atom << endl;
#endif
    }
    else if (characters_processed && '.' == *s)
    {
      if (0 != paren_level)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed,
                              "Mismatched parentheses");
        return 0;
      }

      if (! valid_end_of_smiles_character(previous_token_was))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed,
                              "Invalid end of smiles character");
        return 0;
      }

      if (PREVIOUS_TOKEN_DOT == previous_token_was)
        smiles_error_message(smiles, characters_to_process, characters_processed, "Consecutive '.' ignored");

      previous_atom = INVALID_ATOM_NUMBER;
      previous_token_was = PREVIOUS_TOKEN_DOT;
      characters_processed++;
      fragments_found++;
    }
    else if (characters_processed && ('/' == *s))
    {
      if (! previous_token_should_be(previous_token_was, 
                PREVIOUS_TOKEN_ATOM_SPECIFIER | PREVIOUS_TOKEN_OPEN_PAREN |
                PREVIOUS_TOKEN_RING_SPECIFIER))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Incorrectly placed / bond");
        return 0;
      }
       
      characters_processed++;

      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
      previous_bond_type = SINGLE_BOND;
      previous_bond_type |= SMI_DIRECTIONAL_UP;
    }
    else if (characters_processed && ('\\' == *s))
    {
      if (! previous_token_should_be(previous_token_was, 
                PREVIOUS_TOKEN_ATOM_SPECIFIER | PREVIOUS_TOKEN_OPEN_PAREN |
                PREVIOUS_TOKEN_RING_SPECIFIER))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Incorrectly placed \\ bond");
        return 0;
      }
       
      characters_processed++;

      previous_token_was = PREVIOUS_TOKEN_BOND_SPECIFIER;
      previous_bond_type = SINGLE_BOND;
      previous_bond_type |= SMI_DIRECTIONAL_DOWN;
    }
    else if (kOpenBrace == *s && characters_processed && (characters_to_process - characters_processed >=8) && kOpenBrace == s[1])
    {
      if (! previous_token_should_be(previous_token_was, PREVIOUS_TOKEN_ATOM_SPECIFIER))
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Coordinates must follow an atom");
        return 0;
      }

      Coordinates coords;

      int nchars = parse_coordinates(smiles, characters_to_process, characters_processed, coords);
      if (0 == nchars)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "cannot parse coordinates");
        return 0;
      }

      characters_processed += nchars;

      _things[_number_elements - 1]->setxyz(coords);
    }
    else
    {
      const Element * e;
      aromaticity_type_t aromatic;
      int fc;
      int hcount;
      int chiral_count;
      int atomic_mass = 0;
      int amap = 0;

      int tmp;
      if ('[' == *s)
        tmp = parse_smiles_token(s, characters_to_process - characters_processed,
                                  e, aromatic, fc, hcount, chiral_count, atomic_mass, amap);

      else
      {
        tmp = parse_smiles_token(s, characters_to_process - characters_processed, e, aromatic);
//      cerr << "After parsing '" << s[0] << "' arom " << aromatic << endl;

        hcount = 0;
        fc = 0;
        chiral_count = 0;
      }

#ifdef DEBUG_BUILD_FROM_SMILES
      cerr << "Parse smiles token '" << s << "'\n";
      cerr << "Element " << e->atomic_number() << " arom = " << aromatic <<
              " fc = " << fc << " hcount = " << hcount;
      if (chiral_count)
        cerr << " chiral = " << chiral_count;

      if (amap)
        cerr << " atom map " << amap;
      cerr << endl;
#endif

      if (0 == tmp)
      {
        smiles_error_message(smiles, characters_to_process, characters_processed, "Cannot parse");
        return 0;
      }

      Atom * a = new Atom(e);
      a->resize(3);

      if (atomic_mass)
        a->set_isotope(atomic_mass);
      if (fc)
        a->set_formal_charge(fc);

//    Oct 2000. Try the underlying ::add() call. Should be OK... 

      resizable_array_p<Atom>::add(a);

      if (AROMATICITY_NOT_DETERMINED != aromatic)
      {
        if (! input_aromatic_structures())
        {
          smiles_error_message(smiles, characters_to_process, characters_processed, "aromatic smiles input not enabled");
          return 0;
        }
        aromatic_atoms[_number_elements - 1] = 1;
      }

//    If any specific hcount value came in make sure to set it as a known value.

      if (hcount)
        a->set_implicit_hydrogens(hcount);

      if (amap > 0)
        a->set_atom_map(amap);

//    cerr << "Atom map " << a->atom_map() << endl;

      if ('[' == *s)
        a->set_implicit_hydrogens_known(1);

      if (INVALID_ATOM_NUMBER != previous_atom)
      {
        bond_type_t bt = SINGLE_BOND;
        if (PREVIOUS_TOKEN_BOND_SPECIFIER == previous_token_was)
          bt = previous_bond_type;

#ifdef DEBUG_BUILD_FROM_SMILES
        cerr << "Making bond between atom " << previous_atom << " and " << _number_elements - 1 << " type " << bt << endl;
        if (previous_atom_chiral_count)    // the previous atom was chiral
          cerr << "Previous atom was chiral\n";
        if (chiral_count)
          cerr << "New atom is chiral\n";
#endif

        _smiles_add_bond(previous_atom, _number_elements - 1, bt);

        if (previous_atom_chiral_count)    // the previous atom was chiral
        {
          if (! _smi_atom_bonded_to_chiral_centre(previous_atom, previous_atom_chiral_count, _number_elements - 1))
          {
            smiles_error_message(smiles, characters_to_process, characters_processed,
                          "Error processing atom attached to chiral atom");
            return 0;
          }
        }
      }

//    Check for chirality of the present atom after making any bonds
//    to the previous atom

      if (chiral_count)    // the present atom is chiral
      {
        if (hcount > 1)
        {
          smiles_error_message(smiles, characters_to_process, characters_processed,
                        "Too many Hydrogens on chiral atom");
          return 0;
        }

        Chiral_Centre * tmp = new Chiral_Centre(_number_elements - 1);
        if (! _smi_process_new_chiral_centre(tmp, hcount))
        {
          delete tmp;
          smiles_error_message(smiles, characters_to_process, characters_processed,
                        "Error processing chiral specifier");
          return 0;
        }
        _chiral_centres.add(tmp);
      }

      previous_atom = _number_elements - 1;
      previous_atom_chiral_count = chiral_count;
      characters_processed += tmp;
      previous_token_was = PREVIOUS_TOKEN_ATOM_SPECIFIER;
    }
//  cerr << "Character '" << *s << " interpreted as " << previous_token_was << " nchars = " << characters_processed << endl;
  }

  if (0 != paren_level)
  {
    smiles_error_message(smiles, characters_to_process, characters_processed, "mismatched parentheses");
    return 0;
  }

  if (! valid_end_of_smiles_character(previous_token_was))
  {
    smiles_error_message(smiles, characters_to_process, characters_processed, "invalid smiles end");
    return 0;
  }

//cerr << "LINE " << __LINE__ << " smiles\n";
//debug_print(cerr);

  return 1;
}

int
Molecule::_build_from_smiles (const char * smiles, int nchars,
                              int * aromatic_atoms)
{
  Smiles_Ring_Status ring_status;

  int rc = _build_from_smiles(smiles, nchars, ring_status, aromatic_atoms);

  if (rc && ! ring_status.complete())
  {
    if (file_scope_display_smiles_interpretation_error_messages)
    {
      cerr << "Hanging ring closures\n";
      ring_status.report_hanging_ring_closures(cerr);
    }
    return 0;
  }

  if (0 == rc)    // failed already
    return 0;

//(void) _assign_directional_bonds();

  if (moleculeio::ignore_all_chiral_information_on_input())
    _chiral_centres.resize(0);
  else if (! _check_for_incomplete_chiral_specifications())
    return 0;

// If we read any aromatic atoms, we must find a kekule form

  if (aromatic_atoms && locate_item_in_array(1, _number_elements, aromatic_atoms) >= 0)
  {
    if (add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens())
      _do_unfix_implicit_hydrogens_on_aromatic_isotopic_atoms(aromatic_atoms);
    else if (unset_all_implicit_hydrogens_known_attributes)
      _unset_all_implicit_hydrogens_known_attributes();

    if (! find_kekule_form(aromatic_atoms))
    {
      if (display_no_kekule_form_message())
        cerr << "Molecule::_build_from_smiles:find kekule form failed '" << _molecule_name << "'";

      if (allow_input_without_valid_kekule_form())
      {
        if (display_no_kekule_form_message())
          cerr << ", using single bonds";
        _molecule_name += " (invalid KEKULE form)";
      }
      else
        rc = 0;

      if (display_no_kekule_form_message())
        cerr << endl;
    }
  }

//if (! discern_cis_trans_bonds())    // no need to check anything
//  ;

  if (_finished_reading_smiles_assign_and_check_directional_bonds())    // did it successfully
    ;
  else if (moleculeio::ignore_bad_cis_trans_input())    // trouble, but we are ignoring errors
  {
    revert_all_directional_bonds_to_non_directional();
    _append_bad_cis_trans_input_text_to_name();
  }
  else            // death
    rc = 0;

  if (unset_implicit_hydrogens_known_if_possible)
    _unset_implicit_hydrogens_known_if_computed_matches();
  else if (unset_all_implicit_hydrogens_known_attributes)
    _unset_all_implicit_hydrogens_known_attributes();

  return rc;
}

int
Molecule::_build_from_smiles (const char * smiles, int nchars)
{
  assert(nchars >= 0);

  resize (0);

  if (nchars == 0) {
    _molecule_name.resize(0);
    return 1;
  }


  if ('.' != smiles[0])
    ;
  else if (1 == nchars)
    return 1;
  else if (' ' == smiles[1])
  {
    _molecule_name.set(smiles + 2, nchars - 2);
    return 1;
  }

  int * tmp;
  if (input_aromatic_structures())
    tmp = new_int(nchars);
  else 
    tmp = nullptr;

  int rc = _build_from_smiles(smiles, nchars, tmp);

  if (tmp)
    delete [] tmp;

  return rc;
}

/*
  Fundamental smiles parsing routine. Really just a front end for 
  _build_from_smiles
*/

int
Molecule::build_from_smiles (const char * smiles)
{
  return build_from_smiles(smiles, static_cast<int>(strlen(smiles)));
}

int
Molecule::build_from_smiles (const char * smiles,
                             int nchars)
{
  return _build_from_smiles(smiles, nchars);
}

int
Molecule::build_from_smiles (const const_IWSubstring & smiles)
{
  return _build_from_smiles(smiles.rawchars(), smiles.nchars());
}

int
Molecule::build_from_smiles (const IWString & smiles)
{
  return _build_from_smiles(smiles.rawchars(), smiles.nchars());
}

int
Molecule::build_from_smiles(const std::string& smiles) {
  return _build_from_smiles(smiles.data(), smiles.size());
}

/*
  We need to be able to specify the datatype name used when writing
  structures.
  This was first implemented to allow creation of the '$GRF' datatype
*/

static IWString datatype_name_for_structure_in_tdt_files = "$SMI<";

int 
set_datatype_name_for_structure_in_tdt_files (const char * new_dtname)
{
  assert(strlen(new_dtname));

  datatype_name_for_structure_in_tdt_files = new_dtname;

  if (! datatype_name_for_structure_in_tdt_files.ends_with('<'))
    datatype_name_for_structure_in_tdt_files += '<';

  return 1;
}

/*
  Writes all the data items except $SMI for a tdt.
*/

int
Molecule::_write_molecule_tdt_pcn (std::ostream & os,
                                    const IWString & comment) const
{
  if (_molecule_name.length())
  {
    if (_molecule_name.contains('>') || _molecule_name.contains('|'))
      os << "PCN<" << '"' << _molecule_name << "\">" << newline_string();
    else
      os << "PCN<" << _molecule_name << '>' << newline_string();
  }

  if (comment.length() && comment != _molecule_name)
    os << "REM<" << comment << '>' << newline_string();

  int nt = _text_info.number_elements();
  for (int i = 0; i < nt; i++)
  {
    os << (*_text_info[i]) << newline_string();
  }

  os << '|' << newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

static int
is_three_dots(const const_IWSubstring & smarts,
              int characters_processed,
              const_IWSubstring & three_dots_qualifier)
{
//cerr << "Checking for three dots, in '" << smarts << "', characters_processed = " << characters_processed << " examine '" << smarts[characters_processed] << "'\n";
  if ('.' != smarts[characters_processed])
    return 0;

  int dots = 1;

  int smarts_length = smarts.length();

  characters_processed++;

  while (characters_processed < smarts_length)
  {
    if ('.' == smarts[characters_processed])
    {
      dots++;
      characters_processed++;
    }
    else
      break;
  }

  if (3 != dots)
    return 0;

  three_dots_qualifier.make_empty();

  if (characters_processed == smarts_length)   // wrong, smarts cannot end with ...
    return 1;

  if (kOpenBrace != smarts[characters_processed])    // no qualifier
    return 1;

  for (int i = characters_processed + 1; i < smarts_length; i++)
  {
    if (kCloseBrace != smarts[i])
      continue;

    smarts.from_to(characters_processed, i, three_dots_qualifier);
    return 1;
  }

// Never found a closing brace, bad news

  smarts.from_to(characters_processed + 3, smarts_length - 1, three_dots_qualifier);

  return 1;
}

static int
is_down_the_bond(const const_IWSubstring & smarts,
              int characters_processed,
              const_IWSubstring & down_the_bond_qualifier) {
  if (smarts[characters_processed] != kOpenBrace) {
    return 0;
  }

  int close_brace = misc2::MatchingOpenCloseChar(smarts, characters_processed, kOpenBrace, kCloseBrace);
  if (close_brace < 0) {
    cerr << "is_down_the_bond:empty or unclosed {} directive\n";
    for (int i = 0; i < 10; ++i) {
      if (characters_processed + i == smarts.nchars()) {
        break;
      }
      cerr << smarts[characters_processed + i];
    }
    cerr << '\n';
    return 0;
  }

  smarts.from_to(characters_processed + 1, close_brace - 1, down_the_bond_qualifier);

  return close_brace - characters_processed;
}

int
Molecule::write_molecule_tdt(std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << smiles() << '>' << newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}

int
Molecule::write_molecule_tdt_unique (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << unique_smiles() << '>' << newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}

int
Molecule::write_molecule_tdt_nausmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << non_aromatic_unique_smiles() << '>' << newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}

/*
  In order to parse the '...' directive, we need to tell our caller what was
  the identity of the last atom created
*/

//#define DEBUG_BUILD_FROM_SMARTS

#define LOOKS_LIKE_SMARTS_BOND(b) ((SINGLE_BOND_SYMBOL == (b)) || (DOUBLE_BOND_SYMBOL == (b)) || (TRIPLE_BOND_SYMBOL == (b)) ||\
                                   (AROMATIC_BOND_SYMBOL == (b)) || ('~' == (b)) || ('@' == (b)) ||\
                                   ('!' == (b)) || ('&' == (b)) || (',' == (b)) ||\
                                   (';' == (b)) || ('^' == (b)) )

/*
  LAST_ATOM_CREATED keeps track of the atom numbers in the current fragment.
  ATOMS_IN_PREVIOUS_DISCONNECTED_SECTIONS is needed in order to be able to assign value
  initial atom numbers

*/

int
Substructure_Atom::_parse_smarts_specifier(const const_IWSubstring & qsmarts,
                                Parse_Smarts_Tmp & pst,
                                int atoms_in_previous_disconnected_sections,
                                Smiles_Ring_Status & ring_status)
{
#ifdef DEBUG_BUILD_FROM_SMARTS
  cerr << "Substructure_Atom::_parse_smarts_specifier: '" << qsmarts << "', " << atoms_in_previous_disconnected_sections << " atoms in previous fragments\n";
#endif

  extending_resizable_array<Substructure_Atom *> & completed = pst.completed();

  int characters_to_process = qsmarts.nchars();
  const char * smarts = qsmarts.rawchars();

  if (nullptr == smi_element_star)     // initialise elements first time through
    initialise_organic_subset();

  if (nullptr == compatability_table)
    initialise_compatability_table();

  int characters_processed = 0;

// We need to be somewhat careful about resizing, as we may be called recursively

// Various stacks to keep track of branches. Probably should be combined into
// a single kind with a new object.

  resizable_array<Substructure_Atom *>  atom_stack;
  resizable_array<int>                  chirality_stack;

  int previous_token_was = PREVIOUS_TOKEN_WAS_NOT_SPECIFIED;

  std::unique_ptr<Substructure_Bond> previous_bond;
//Substructure_Bond * previous_bond = nullptr;
  Substructure_Atom * previous_atom = nullptr;
  int previous_atom_chiral_count = 0;

  int atoms_this_fragment = 0;

  int paren_level = 0;

  const_IWSubstring three_dots_qualifier;    // any qualifier after a ... directive
  const_IWSubstring down_the_bond_qualifier;  // qualifier arter a -{} directive.

  while (characters_processed < characters_to_process)
  {
    const char * s = smarts + characters_processed;

#ifdef DEBUG_BUILD_FROM_SMARTS
    cerr << "Examining smarts token '" << *s << "', previous " << previous_token_was << endl;
#endif

    if (characters_processed && (isdigit(*s) || '%' == *s))    // we do not properly handle chirality here, fix sometime...
    {
      int ring_number;
      int nchars;
      if (! fetch_ring_number(s, characters_to_process - characters_processed, ring_number, nchars))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Invalid ring specification");
        return 0;
      }

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Is ring number " << ring_number << endl;
#endif

      std::unique_ptr<Substructure_Bond> b;
      bond_type_t bt = SINGLE_BOND; 
      if (PREVIOUS_TOKEN_WAS_BOND == previous_token_was)
      {
        b.reset(previous_bond.release());
        bt = b->types_matched();
      }
      else
      {
        b.reset(new Substructure_Bond);
        b->make_single_or_aromatic();
        bt = NOT_A_BOND;
      }

      atom_number_t other_end;    // will be set when a ring closure
      if (ring_status.encounter(ring_number, previous_atom->unique_id(), other_end, bt))   // ring closing.
      {
        assert(nullptr != completed[other_end]);

#ifdef DEBUG_BUILD_FROM_SMARTS
        cerr << "Ring closure, atom at other end is " << other_end << " prev " << previous_atom->unique_id() << endl;
#endif

        if (previous_atom->is_bonded_to(other_end))
        {
          cerr << "Substructure_Atom::_parse_smarts_specifier:two membered ring encountered\n";
          return 0;
        }

        b->set_atom(completed[other_end]);
        b->set_type(bt);
        previous_atom->_add_bond(b.release());
      }

      characters_processed += nchars;
      previous_token_was = PREVIOUS_TOKEN_WAS_RING;
    }
    else if (LOOKS_LIKE_SMARTS_BOND(*s) && check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
    {
      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Consecutive bonds");
        return 0;
      }

      previous_bond.reset(new Substructure_Bond);
      int ncb;
      if (! previous_bond->construct_from_smarts(s, characters_to_process - characters_processed, ncb))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot parse bond specifier");
        return 0;
      }

      characters_processed += ncb;
#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Recognised as bond\n";
#endif
    }
    else if (is_three_dots(qsmarts, characters_processed, three_dots_qualifier) && check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_THREE_DOTS))
    {
      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot have ... after bond");
        return 0;
      }

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Got three dots, qualifier '" << three_dots_qualifier << "'\n";
#endif

      characters_processed += 3 + three_dots_qualifier.length();

      ThreeDots * three_dots = new ThreeDots(previous_atom->unique_id(),
                                             pst.last_query_atom_created() + 1);

      if (three_dots_qualifier.length())
        three_dots->set_qualifier(three_dots_qualifier);

      pst.add_three_dots(three_dots);

      Substructure_Atom * a = new Substructure_Atom;
      pst.add_root_atom(a);

      const_IWSubstring newsmarts(qsmarts);
      newsmarts.remove_leading_chars(characters_processed);

      if (0 == newsmarts.length())
      {
        cerr << "Substructure_Atom::parse_smiles_token:smarts cannot end in ... operator\n";
        delete a;
        return 0;
      }

      return a->_parse_smarts_specifier(newsmarts, pst, atoms_in_previous_disconnected_sections + atoms_this_fragment, ring_status);
    }
    else if (is_down_the_bond(qsmarts, characters_processed, down_the_bond_qualifier) && check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_DOWN_THE_BOND))
    {
      std::unique_ptr<DownTheBond> dtb = std::make_unique<DownTheBond>(previous_atom->unique_id());
      if (! dtb->Build(down_the_bond_qualifier)) {
        cerr << "Substructure_Atom::parse_smiles_token:invalid down the bond {" << down_the_bond_qualifier << "}\n";
        return 0;
      }
      pst.add_down_the_bond(dtb.release());
      characters_processed += 1 + down_the_bond_qualifier.length() + 1;
      previous_token_was = PREVIOUS_TOKEN_WAS_DOWN_THE_BOND;
    }
    else if (' ' == *s)
    {
#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Found space, ending parsing\n";
#endif
      if (0 != paren_level)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Un-closed parenthesis");
        return 0;
      }

      return 1;
    }
    else if (kOparen == *s)
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_OPEN_PAREN))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, 
                              "Branch specifier can only follow an atom or ring specifier");
        return 0;
      }

//    assert(!previous_bond);

//    Push the various stacks.

      paren_level++;

      atom_stack.add(previous_atom);
      chirality_stack.add(previous_atom_chiral_count);

      characters_processed++;
    }
    else if (characters_processed && kCparen == *s)
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_CLOSE_PAREN))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Close paren can only follow an atom");
        return 0;
      }
      if (atom_stack.empty())
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Parenthesis mismatch");
        return 0;
      }

      if (0 == paren_level)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Too many close paren's");
        return 0;
      }

//    Jun 98. Handle something like '[CD3](:0)'   note it is a digit '0' rather than 'O' (oxygen)

      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Closing paren after bond??");
        return 0;
      }

      paren_level--;

      previous_atom = atom_stack.pop();
      previous_token_was = PREVIOUS_TOKEN_WAS_ATOM;
      previous_atom_chiral_count = chirality_stack.pop();
      characters_processed++;
    }
    else if ('.' == *s)
    {
      smiles_error_message(smarts, characters_to_process, characters_processed, 
                            "disconnect '.' not allowed in smarts");
      return 0;
    }
    else if (characters_processed && ('/' == *s))
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Incorrectly placed / bond");
        return 0;
      }
       
      characters_processed++;

      previous_bond.reset(new Substructure_Bond());
      previous_bond->make_single_or_aromatic();
    }
    else if (characters_processed && ('\\' == *s))
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Incorrectly placed \\ bond");
        return 0;
      }
       
      characters_processed++;

      previous_bond.reset(new Substructure_Bond());
      previous_bond->make_single_or_aromatic();
    }
    else
    {
      int save_previous_token_was = previous_token_was;
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_ATOM))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "incorrectly placed atom");
        return 0;
      }

      Substructure_Atom * a;
      if (0 == atoms_this_fragment)
        a = this;
      else
        a = new Substructure_Atom;

      int tmp;
      if ('[' == *s)
        tmp = a->construct_from_smarts_token(s, characters_to_process - characters_processed);
      else if ('*' == *s)
        tmp = 1;
      else
        tmp = a->construct_from_smiles_token(s, characters_to_process - characters_processed);

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Built atom processing '" << s << "' tmp = " << tmp << " previous_token_was " << previous_token_was << '\n';
#endif

      if (tmp == 0) {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot parse smarts");
        if (atoms_this_fragment > 0)      // only delete A if we created it
          delete a;
        return 0;
      }

      const int my_atom_number = atoms_in_previous_disconnected_sections + atoms_this_fragment;
      completed[my_atom_number] = a;
      pst.set_last_query_atom_created(my_atom_number);
      if (save_previous_token_was == PREVIOUS_TOKEN_WAS_DOWN_THE_BOND) {
        pst.down_the_bond().last_item()->set_a2(my_atom_number);
      }

      if (a->initial_atom_number() < 0)    // only change it if it is unset
        a->set_initial_atom_number(my_atom_number);

      a->set_unique_id(my_atom_number);

      atoms_this_fragment++;

      if (nullptr != previous_atom)
      {
        Substructure_Bond * b;

        if (previous_bond)
          b = previous_bond.release();
        else
        {
          b = new Substructure_Bond();
          b->make_single_or_aromatic();
        }

#ifdef DEBUG_BUILD_FROM_SMARTS
        cerr << "Making bond between atom " << previous_atom->unique_id() << " and " << atoms_this_fragment - 1 << endl;
#endif

        b->set_atom(previous_atom);
        a->_add_bond(b);
      }

      previous_atom = a;

//    No chirality in smarts

      characters_processed += tmp;
      previous_token_was = PREVIOUS_TOKEN_WAS_ATOM;
    }
#ifdef DEBUG_BUILD_FROM_SMARTS
    cerr << "Character '" << *s << " interpreted as " << previous_token_was << " nchars = " << characters_processed << endl;
#endif
  }

  if (0 != paren_level)
  {
    smiles_error_message(smarts, characters_to_process, characters_processed, "mismatched parentheses");
    return 0;
  }

  if (PREVIOUS_TOKEN_WAS_BOND == previous_token_was)
  {
    smiles_error_message(smarts, characters_to_process, characters_processed, "a smarts cannot end with a bond");
    return 0;
  }

  if (! pst.DownTheBondSpecificationsComplete()) {
    cerr << "Substructure_Atom::_parse_smarts_specifier:unclosed down the bond\n";
    return 0;
  }
#ifdef DEBUG_BUILD_FROM_SMARTS
  cerr << "Substructure_Atom::_parse_smarts_specifier: returning 1\n";
#endif

  return 1;
}

/*
  A single Substructure_Atom object can parse a complete smarts, because
  connected atoms become children
*/

int
Substructure_Atom::parse_smarts_specifier(const const_IWSubstring & smarts,
                                          Parse_Smarts_Tmp & pst,
                                          int atoms_in_previous_disconnected_sections)
{
  const_IWSubstring mysmarts(smarts);
  mysmarts.strip_leading_blanks();

  Smiles_Ring_Status ring_status;

  int rc = _parse_smarts_specifier(mysmarts, pst, atoms_in_previous_disconnected_sections, ring_status);

  Substructure_Atom::count_attributes_specified();

  ok_recursive();

  if (0 == rc)
    return 0;

  if (! ring_status.complete())
  {
    ring_status.report_hanging_ring_closures(cerr);
    return 0;
  }

  return rc;
}

int
Substructure_Atom::parse_smarts_specifier(const const_IWSubstring & smarts)
{
  Parse_Smarts_Tmp pst;

  pst.set_natoms(smarts.length());

  int notused = 0;

  return parse_smarts_specifier(smarts, pst, notused);
}

void
reset_smi_file_scope_variables()
{
  unset_implicit_hydrogens_known_if_possible = 0;
  unset_all_implicit_hydrogens_known_attributes = 0;
  identifier_dataitem  = "PCN<";
  dataitems_to_append.resize(0);
  append_dataitem_content = 1;
  file_scope_display_smiles_interpretation_error_messages = 1;
  smi_element_star = nullptr;   // initialised for check in parse_smiles_token
  ignore_tdts_with_no_smiles = 0;
  smiles_tag = "$SMI<";
  if (nullptr != compatability_table)
  {
    delete [] compatability_table;
    compatability_table = nullptr;
  }
  datatype_name_for_structure_in_tdt_files = "$SMI<";

  return;
}
