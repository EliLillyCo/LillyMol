/*
  This function initialises an array of elements.
*/

#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <limits>
#include <assert.h>

using std::cerr;
using std::endl;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "iwstring_data_source.h"
#include "element.h"
#include "misc2.h"

static resizable_array_p<Element> elements;

#define PLAUSIBLE_ATOMIC_NUMBER(q) ((q) >= 0 && (q) <= HIGHEST_ATOMIC_NUMBER)

#define OUTER_SHELL_ELECTRONS_NOT_KNOWN -18

static int include_isotopes_in_smiles = 1;

void
set_include_isotopes_in_smiles (int s)
{
  include_isotopes_in_smiles = s;
}

static int explicit_hydrogens_need_square_brackets_in_smiles = 1;

void
set_explicit_hydrogens_need_square_brackets_in_smiles (int s)
{
  explicit_hydrogens_need_square_brackets_in_smiles = s;

  if (elements.number_elements())    // will pretty much always be the case
  {
    Element * h = const_cast<Element *> (get_element_from_atomic_number (1));
    assert (NULL != h);

    h->set_needs_square_brackets (0);
  }

  return;
}

/*
  When atomic symbols have arbitrary length, we don't do any hashing on them
*/

static int _atomic_symbols_can_have_arbitrary_length = 0;

void
set_atomic_symbols_can_have_arbitrary_length (int s)
{
  _atomic_symbols_can_have_arbitrary_length = s;

  return;
}

int
atomic_symbols_can_have_arbitrary_length()
{
  return _atomic_symbols_can_have_arbitrary_length;
}

static int display_strange_chemistry_messages = 1;

void
set_display_strange_chemistry_messages(int s)
{
  display_strange_chemistry_messages = s;
}

/*
  Search the non periodic table elements. 
  If this ever becomes slow, use a hash
*/

const Element *
get_element_from_long_symbols (const char * asymbol,
                               int nchars)
{
//#define DEBUG_GET_ELEMENT_FROM_LONG_SYMBOLS
#ifdef DEBUG_GET_ELEMENT_FROM_LONG_SYMBOLS
  cerr << "get_element_from_long_symbols: getting element for '";
  cerr.write(asymbol, nchars);
  cerr << "' nchars " << nchars << "\n";
#endif

  const_IWSubstring tmp(asymbol, nchars);

  for (int i = HIGHEST_ATOMIC_NUMBER + 1; i < elements.number_elements(); i++)
  {
    if (tmp == elements[i]->symbol())
      return elements[i];
  }

  return NULL;
}


/*
  Make sure the hash table is correct.
  We need 26 for the single letter elements.
  For each two letter elements, we need 26 * 26 entries
  We also hash '*', and don't want it to collide with anything

  Nov 02. Allow elements like 'R1', so change all the 26 to 36
  May 03. Allow R# as an element
  Jan 04. Allow lowercase letters
*/

#define SIZE_OF_ELEMENT_HASH_TABLE (26 + 36 * 26 + 3 + 26)

static const Element * ehash[SIZE_OF_ELEMENT_HASH_TABLE];

static int
element_symbol_hash_function(const char * s,
                             const int nchars)
{
  if (nchars > 2)
  {
    if (_atomic_symbols_can_have_arbitrary_length)
      return -1;

    cerr << "element_symbol_hash_function:invalid atomic symbol '";
    cerr.write(s, nchars);
    cerr << "'\n";
    return -1;
  }

  int rc = *s - 'A';

  if (rc >= 0 && rc < 26)    // is an uppercase letter
  {
    if (1 == nchars)
      return rc;

    int tmp = s[1] - 'a';
    if (tmp >= 0 && tmp < 26)    // 2nd char is a lowercase letter
      return 26 + 36 * rc + tmp;

    tmp = s[1] - '0';
    if (tmp >= 0 && tmp <= 9)     // 2nd char is a digit
      return 26 + 36 * rc + 26 + tmp;

    tmp = s[1] - 'A';    // 2nd char is an uppercase letter, treat as lowercase (hmmm, is this such a good idea?)
    if (tmp >= 0 && tmp < 26)
      return 26 + 36 * rc + tmp;

    if ('R' == s[0] && '#' == s[1])
      return 26 + 36 * 26 + 2;
  }
  else if ('*' == *s && 1 == nchars)
    return 26 + 36 * 26;
  else if (2 == nchars && '?' == s[0] && '?' == s[1])
    return 26 + 36 * 26 + 1;

  rc = *s - 'a';

  if (rc >= 0 && rc < 26 && 1 == nchars)   // single lowercase letter
    return 26 + 36 * 26 + 3 + rc;          // at the end of the list

  cerr << "Illegal atomic symbol (len " << nchars << ") '";
  cerr.write(s, nchars);
  cerr << "' hash = " << rc << endl;
  return 26+36*26;
}

int
print_element_hash_table(std::ostream & os)
{
  for (int i = 0; i < SIZE_OF_ELEMENT_HASH_TABLE; i++)
  {
    const Element * e = ehash[i];

    if (NULL == e)
      continue;

    os << "Element hash " << i << " element '" << e->symbol() << "' hash value " << e->atomic_symbol_hash_value() << endl;
  }

  return os.good();
}

int
symbol_for_atomic_symbol_hash_value (int h,
                                     IWString & s)
{
  if (h < 0 || h >= SIZE_OF_ELEMENT_HASH_TABLE)
  {
    cerr << "symbol_for_atomic_symbol_hash_value:invalid hash value " << h << endl;
    return 0;
  }

//cerr << "What is the element for hash " << h << endl;

  if (h < 26)
  {
    s = 'A' + h;
    return 1;
  }

  if (26 + 36 * 36 == h)
  {
    s = '*';
    return 1;
  }

  if (26 + 36 * 36 + 1 == h)
  {
    s = "??";
    return 1;
  }

  if (h >= 26 + 36 * 26 + 3)
  {
    s.add('a' + h - (26 + 36 * 26 + 3));
    return 1;
  }

  h = h - 26;

  int c2 = h % 36;
  h = h / 36;

  assert (h >= 0 && h < 26);

  s.add('A' + h);

  if (c2 < 26)
    s.add('a' + c2);
  else
    s.add('0' + c2 - 26);

  return 1;
}

/*
  We are relying on the idea that most symbols will be short
*/

static int
element_symbol_hash_function_arbitrary_length (const char * s,
                                               const int nchars)
{
  uint64_t rc = SIZE_OF_ELEMENT_HASH_TABLE;

  for (int i = 0; i < nchars; ++i)
  {
    rc = 256 * rc + static_cast<unsigned int>(s[i]);
//  cerr << " element_symbol_hash_function_arbitrary_length i " << i << " rc " << rc << endl;
  }

  if (rc < SIZE_OF_ELEMENT_HASH_TABLE)    // make sure it is out of range
    return std::numeric_limits<int>::max() - rc;

  uint64_t m = static_cast<int>(std::numeric_limits<int>::max());

  while (rc > m)    // seems unlikely
  {
    rc = rc / 2;
  }

  return static_cast<int>(rc);
}

static void
init_elements (void)
{
  assert (0 == elements.number_elements());

  for (int i = 0; i < SIZE_OF_ELEMENT_HASH_TABLE; i++)
  {
    ehash[i] = NULL;
  }

  elements.resize(HIGHEST_ATOMIC_NUMBER + 20);

  for (int i = 0; i <= HIGHEST_ATOMIC_NUMBER; i++)
  {
    elements.add(new Element(i));
    elements.last_item()->set_unique_id(elements.number_elements());
  }

  return;
}

/*
  Common operation is to see if a symbol is a valid hashed element
*/

static int
check_valid_element_via_hash (const char * s, int nchars,
                              const Element * & result)
{
  const int hash = element_symbol_hash_function(s, nchars);

#ifdef DEBUG_CHECK_VALID_ELEMENT_VIA_HASH
  cerr << "check_valid_element_via_hash '";
  cerr.write(s, nchars);
  cerr << "' hash " << hash;
  if (hash < 0 || hash >= SIZE_OF_ELEMENT_HASH_TABLE)
    cerr << " out of range\n";
  else if (NULL == ehash[hash])
    cerr << " no element\n";
  else
    cerr << " element " << ehash[hash]->symbol() << endl;
#endif

  if (hash < 0)
    return 0;

  if (hash >= SIZE_OF_ELEMENT_HASH_TABLE)
    return 0;

  if (NULL == ehash[hash])
    return 0;

  result = ehash[hash];

  return nchars;
}

/*
  In order to get the elements created automatically, let's make a class
  which will get instantiated. 
*/

class element_creator
{
  private:
  public:
    element_creator();
};

element_creator::element_creator()
{
  init_elements();
}

static element_creator foo;      // the constructor will cause the elements to be created.

/*
  The variable _outer_shell_electrons. is initialised for only a few elements.
  It is really a count of the outer shell electrons (s and p)

  I haven't filled in all the metal flags.
*/

/*
  AMW data from the Handbook of Chemistry and Physics, on-line May 2001
*/

void
Element::_default_values (atomic_number_t i)
{
  assert (PLAUSIBLE_ATOMIC_NUMBER(i) || NOT_AN_ELEMENT == i); 

  _atomic_number  = i;
  _atomic_mass    = static_cast<atomic_mass_t>(0.0);
  _normal_isotope = 0;
  _organic        = 0;
  _metal          = 0;
  _normal_valence = VALENCE_NOT_DEFINED;
  _outer_shell_electrons = OUTER_SHELL_ELECTRONS_NOT_KNOWN;
  _exact_mass     = 0.0;

  _needs_square_brackets = 1;

  _atomic_symbol_hash_value = -1;

  _permanent_aromatic = 0;

  _natural_peptide = 0;

  switch (i)
  {
    case 0:
      _symbol = "*";
      _needs_square_brackets = 0;
      break;

    case 1:
      _symbol = "H";
      _normal_isotope = 1;
      _atomic_mass    = static_cast<atomic_mass_t>(1.00794);
      _needs_square_brackets = explicit_hydrogens_need_square_brackets_in_smiles;
      _normal_valence = 1;
      _organic = 1;
      _exact_mass = 1.007825035;
      break;

    case 2:
      _symbol = "He";
      _normal_isotope = 4;
      _atomic_mass    = static_cast<atomic_mass_t>(4.002602);
      _normal_valence = 0;
      _exact_mass     = 4.00260324;
      break;

    case 3:
      _symbol = "Li";
      _normal_isotope = 7;
      _atomic_mass    = static_cast<atomic_mass_t>(6.941);
      _normal_valence = 1;
      _exact_mass     = 7.0160030;
      _outer_shell_electrons = 1;
      _metal          = 1;
      break;

    case 4:
      _symbol = "Be";
      _normal_isotope = 9;
      _atomic_mass    = static_cast<atomic_mass_t>(9.012182);
      _exact_mass     = 9.0121822;
      break;

    case 5:
      _symbol = "B";
      _normal_isotope = 11;
      _atomic_mass    = static_cast<atomic_mass_t>(10.811);
      _needs_square_brackets = 0;
      _normal_valence = 3;
      _exact_mass     = 11.0093054;
      _outer_shell_electrons = 3;
      break;

    case 6:
      _symbol = "C";
      _aromatic_symbol = "c";
      _normal_isotope = 12;
      _atomic_mass    = static_cast<atomic_mass_t>(12.0107);
      _needs_square_brackets = 0;
      _normal_valence = 4;
      _organic = 1;
      _outer_shell_electrons = 4;
      _exact_mass     = 12.0000;
      break;

    case 7:
      _symbol = "N";
      _aromatic_symbol = "n";
      _normal_isotope = 14;
      _atomic_mass    = static_cast<atomic_mass_t>(14.00674);
      _needs_square_brackets = 0;
      _normal_valence = 3;
      add_alternate_valence(5);
      _organic = 1;
      _outer_shell_electrons = 5;
      _exact_mass            = 14.003074002;
      break;

    case 8:
      _symbol = "O";
      _aromatic_symbol = "o";
      _normal_isotope = 16;
      _atomic_mass    = static_cast<atomic_mass_t>(15.9994);
      _needs_square_brackets = 0;
      _normal_valence = 2;
      _organic = 1;
      _outer_shell_electrons = 6;
      _exact_mass     = 15.99491463;
      break;

    case 9:
      _symbol = "F";
      _normal_isotope = 19;
      _atomic_mass    = static_cast<atomic_mass_t>(18.9984032);
      _needs_square_brackets = 0;
      _normal_valence = 1;
      _organic = 1;
      _outer_shell_electrons = 7;
      _exact_mass     = 18.99840322;
      break;

    case 10:
      _symbol = "Ne";
      _normal_isotope = 20;
      _atomic_mass = static_cast<atomic_mass_t>(20.1797);
      _normal_valence = 0;
      _exact_mass     = 19.9924356;
      break;

    case 11:
      _symbol = "Na";
      _normal_isotope = 23;
      _atomic_mass = static_cast<atomic_mass_t>(22.989769);
      _normal_valence = 1;
      _exact_mass     = 22.9897677;
      _outer_shell_electrons = 1;
      _metal          = 1;
      break;

    case 12:
      _symbol = "Mg";
      _normal_isotope = 24;
      _atomic_mass = static_cast<atomic_mass_t>(24.3050);
      _normal_valence = 2;
      _exact_mass     = 23.9850423;
      _outer_shell_electrons = 2;
      _metal          = 1;
      break;

    case 13:
      _symbol = "Al";
      _normal_isotope = 27;
      _atomic_mass = static_cast<atomic_mass_t>(26.981539);
      _exact_mass     = 26.9815386;
      _metal          = 1;
      break;

    case 14:
      _symbol = "Si";
      _normal_isotope = 28;
      _atomic_mass = static_cast<atomic_mass_t>(28.0855);
      _normal_valence = 4;
      _outer_shell_electrons = 4;
      _exact_mass     = 27.9769271;
      break;

    case 15:
      _symbol = "P";
      _aromatic_symbol = "p";
      _normal_isotope = 31;
      _atomic_mass    = static_cast<atomic_mass_t>(30.973762);
      _needs_square_brackets = 0;
      _normal_valence = 3;
      add_alternate_valence(5);
      add_alternate_valence(6);   // added Jul 2000 for Natural Products
      _organic = 1;
      _outer_shell_electrons = 5;
      _exact_mass     = 30.9737620;
      break;


    case 16:
      _symbol = "S";
      _aromatic_symbol = "s";
      _normal_isotope = 32;
      _atomic_mass    = static_cast<atomic_mass_t>(32.065);
      _needs_square_brackets = 0;
      _normal_valence = 2;
      add_alternate_valence(4);   // added 11 Jun 96
      add_alternate_valence(6);
//    add_alternate_valence(8);   // omit, Feb 2004
      _outer_shell_electrons = 6;
      _organic = 1;
      _exact_mass     = 31.97207070;
      break;

    case 17:
      _symbol = "Cl";
      _normal_isotope = 35;
      _atomic_mass    = static_cast<atomic_mass_t>(35.4527);
      _normal_valence = 1;
      _needs_square_brackets = 0;
      add_alternate_valence(7);
      _outer_shell_electrons = 7;
      _organic = 1;
      _exact_mass     = 34.968852721;
      break;

    case 18:
      _symbol = "Ar";
      _normal_isotope = 40;
      _atomic_mass = static_cast<atomic_mass_t>(39.948);
      _normal_valence = 0;
      _exact_mass     = 39.9623837;
      break;

    case 19:
      _symbol = "K";
      _normal_isotope = 39;
      _atomic_mass    = static_cast<atomic_mass_t>(39.0983);
      _normal_valence = 1;
      _exact_mass     = 38.9637074;
      _outer_shell_electrons = 1;
      _metal          = 1;
      break;

    case 20:
      _symbol = "Ca";
      _normal_isotope = 40;
      _atomic_mass = static_cast<atomic_mass_t>(40.078);
      _normal_valence = 2;
      _exact_mass     = 39.9625906;
      _outer_shell_electrons = 2;
      _metal          = 1;
      break;

    case 21:
      _symbol = "Sc";
      _normal_isotope = 45;
      _atomic_mass = static_cast<atomic_mass_t>(44.955912);
      _exact_mass     = 44.95591000;
      break;

    case 22:
      _symbol = "Ti";
      _normal_isotope = 48;
      _atomic_mass = static_cast<atomic_mass_t>(47.867);
      _exact_mass     = 47.9479473;
      break;

    case 23:
      _symbol = "V";
      _normal_isotope = 51;
      _atomic_mass = static_cast<atomic_mass_t>(50.9415);
      _exact_mass     = 50.9439617;
      break;

    case 24:
      _symbol = "Cr";
      _normal_isotope = 52;
      _atomic_mass = static_cast<atomic_mass_t>(51.9961);
      _exact_mass     = 51.9405098;
      break;

    case 25:
      _symbol = "Mn";
      _normal_isotope = 55;
      _atomic_mass = static_cast<atomic_mass_t>(54.938045);
      _exact_mass     = 54.9380471;
      break;

    case 26:
      _symbol = "Fe";
      _normal_isotope = 56;
      _atomic_mass = static_cast<atomic_mass_t>(55.845);
      _exact_mass     = 55.9349393;
      break;

    case 27:
      _symbol = "Co";
      _normal_isotope = 59;
      _atomic_mass = static_cast<atomic_mass_t>(58.933195);
      _exact_mass     = 58.9331976;
      break;

    case 28:
      _symbol = "Ni";
      _normal_isotope = 59;
      _atomic_mass = static_cast<atomic_mass_t>(58.6934);
      _exact_mass     = 57.9353462;
      break;

    case 29:
      _symbol = "Cu";
      _normal_isotope = 64;
      _atomic_mass = static_cast<atomic_mass_t>(63.546);
      _exact_mass     = 62.9295989;
      break;

    case 30:
      _symbol = "Zn";
      _normal_isotope = 65;
      _atomic_mass = static_cast<atomic_mass_t>(65.38);
      _exact_mass     = 63.9291448;
      break;

    case 31:
      _symbol = "Ga";
      _normal_isotope = 70;
      _atomic_mass = static_cast<atomic_mass_t>(69.723);
      _exact_mass     = 68.925580;
      break;

    case 32:
      _symbol = "Ge";
      _normal_isotope = 73;
      _atomic_mass = static_cast<atomic_mass_t>(72.64);
      _exact_mass     = 73.9211774;
      break;

    case 33:
      _symbol = "As";
      _normal_isotope = 75;
      _atomic_mass = static_cast<atomic_mass_t>(74.92160);
      _exact_mass     = 74.9215942;
      break;

    case 34:
      _symbol = "Se";
      _aromatic_symbol = "se";
      _normal_isotope = 79;
      _atomic_mass = static_cast<atomic_mass_t>(78.96);
      _normal_valence = 2;
      _outer_shell_electrons = 6;
      _exact_mass     = 79.9165196;
      break;

    case 35:
      _symbol = "Br";
      _normal_isotope = 80;
      _atomic_mass    = static_cast<atomic_mass_t>(79.904);
      _normal_valence = 1;
      add_alternate_valence(5);   // 4 and 7 exist sometimes too
      add_alternate_valence(7);   // 4 and 7 exist sometimes too
      _needs_square_brackets = 0;
      _outer_shell_electrons = 7;
      _organic = 1;
      _exact_mass     = 78.9183361;
      break;
   
    case 36:
      _symbol = "Kr";
      _normal_isotope = 84;
      _atomic_mass = static_cast<atomic_mass_t>(83.798);
      _exact_mass     = 83.911507;
      break;
   
    case 37:
      _symbol = "Rb";
      _normal_isotope = 85;
      _atomic_mass = static_cast<atomic_mass_t>(85.4678);
      _exact_mass     = 84.911794;
      break;
   
    case 38:
      _symbol = "Sr";
      _normal_isotope = 88;
      _atomic_mass = static_cast<atomic_mass_t>(87.62);
      _exact_mass     = 87.9056188;
      break;
   
    case 39:
      _symbol = "Y";
      _normal_isotope = 89;
      _atomic_mass = static_cast<atomic_mass_t>(88.90585);
      _exact_mass     = 88.905849;    // possible problem at Trace Sciences
      break;
   
    case 40:
      _symbol = "Zr";
      _normal_isotope = 91;
      _atomic_mass = static_cast<atomic_mass_t>(91.224);
      _exact_mass     = 89.9047026;
      break;
   
    case 41:
      _symbol = "Nb";
      _normal_isotope = 93;
      _atomic_mass = static_cast<atomic_mass_t>(92.90638);
      _exact_mass     = 92.9063772;
      break;
   
    case 42:
      _symbol = "Mo";
      _normal_isotope = 96;
      _atomic_mass = static_cast<atomic_mass_t>(95.96);
      _exact_mass     = 97.9054073;
      break;
   
    case 43:
      _symbol = "Tc";
      _normal_isotope = 99;
      _atomic_mass = static_cast<atomic_mass_t>(98);     // not sure why these are different?
      break;
   
    case 44:
      _symbol = "Ru";
      _normal_isotope = 101;
      _atomic_mass = static_cast<atomic_mass_t>(101.07);
      _exact_mass     = 101.9043485;
      break;
   
    case 45:
      _symbol = "Rh";
      _normal_isotope = 103;
      _atomic_mass = static_cast<atomic_mass_t>(102.90550);
      _exact_mass     = 102.9055;
      break;
   
    case 46:
      _symbol = "Pd";
      _normal_isotope = 106;
      _atomic_mass = static_cast<atomic_mass_t>(106.42);
      _exact_mass     = 105.903478;
      break;
   
    case 47:
      _symbol = "Ag";
      _normal_isotope = 108;
      _atomic_mass = static_cast<atomic_mass_t>(107.8682);
      _exact_mass     = 106.905092;
      break;

    case 48:
      _symbol = "Cd";
      _normal_isotope = 112;
      _atomic_mass = static_cast<atomic_mass_t>(112.411);
      _exact_mass     = 113.903357;
      break;

    case 49:
      _symbol = "In";
      _normal_isotope = 115;
      _atomic_mass = static_cast<atomic_mass_t>(114.818);
      _exact_mass     = 114.903882;
      break;

    case 50:
      _symbol = "Sn";
      _normal_isotope = 119;
      _atomic_mass = static_cast<atomic_mass_t>(118.710);
      _exact_mass     = 119.9021991;
      break;
   
    case 51:
      _symbol = "Sb";
      _normal_isotope = 122;
      _atomic_mass = static_cast<atomic_mass_t>(121.760);
      _exact_mass     = 120.9038212;
      break;
   
    case 52:
      _symbol = "Te";
      _aromatic_symbol = "te";
      _normal_isotope = 128;
      _atomic_mass = static_cast<atomic_mass_t>(127.60);
      _normal_valence = 2;
      _outer_shell_electrons = 6;
      _exact_mass     = 129.906229;
      break;
   
    case 53:
      _symbol = "I";
      _normal_isotope = 127;
      _atomic_mass    = static_cast<atomic_mass_t>(126.90447);
      _needs_square_brackets = 0;
      _normal_valence = 1;
      add_alternate_valence(3);
      add_alternate_valence(5);
      add_alternate_valence(7);
      _outer_shell_electrons = 7;
      _organic = 1;
      _exact_mass     = 126.904473;
      break;

    case 54:
      _symbol = "Xe";
      _normal_isotope = 131;
      _atomic_mass = static_cast<atomic_mass_t>(131.293);
      _exact_mass     = 131.904144;
      break;

    case 55:
      _symbol = "Cs";
      _normal_isotope = 133;
      _atomic_mass    = static_cast<atomic_mass_t>(132.90545);
      _exact_mass     = 132.905429;    // possible problem at Trace Sciences
      break;

    case 56:
      _symbol = "Ba";
      _normal_isotope = 137;
      _atomic_mass    = static_cast<atomic_mass_t>(137.327);
      _exact_mass     = 137.905232;
      _metal          = 1;
      break;

    case 57:
      _symbol = "La";
      _normal_isotope = 139; 
      _atomic_mass = static_cast<atomic_mass_t>(138.9055);
      _exact_mass     = 138.906347;
      break;

    case 58:
      _symbol = "Ce";
      _normal_isotope = 140;
      _atomic_mass = static_cast<atomic_mass_t>(140.116);
      _exact_mass     = 139.905433;
      break;

    case 59:
      _symbol = "Pr";
      _normal_isotope = 147;
      _atomic_mass = static_cast<atomic_mass_t>(140.90765);
      _exact_mass     = 140.907647;
      break;

    case 60:
      _symbol = "Nd";
      _normal_isotope = 144;
      _atomic_mass = static_cast<atomic_mass_t>(144.242);
      _exact_mass     = 141.907719;   // Trace Sciences says this is the most abundant
      break;

    case 61:
      _symbol = "Pm";
      _normal_isotope = 147;
      _atomic_mass = static_cast<atomic_mass_t>(145.0);
      break;

    case 62:
      _symbol = "Sm";
      _normal_isotope = 150;
      _atomic_mass = static_cast<atomic_mass_t>(150.36);
      _exact_mass     = 151.919728;
      break;

    case 63:
      _symbol = "Eu";
      _normal_isotope = 152;
      _atomic_mass = static_cast<atomic_mass_t>(151.964);
      _exact_mass     = 152.921225;
      break;

    case 64:
      _symbol = "Gd";
      _normal_isotope = 157;
      _atomic_mass = static_cast<atomic_mass_t>(157.25);
      _exact_mass     = 157.924019;
      break;

    case 65:
      _symbol = "Tb";
      _normal_isotope = 159;
      _atomic_mass = static_cast<atomic_mass_t>(158.92534);
      _exact_mass     = 158.925342;
      break;

    case 66:
      _symbol = "Dy";
      _normal_isotope = 163;
      _atomic_mass = static_cast<atomic_mass_t>(162.50);
      _exact_mass     = 163.929171;
      break;

    case 67:
      _symbol = "Ho";
      _normal_isotope = 165;
      _atomic_mass = static_cast<atomic_mass_t>(164.93032);
      _exact_mass     = 164.930319;
      break;

    case 68:
      _symbol = "Er";
      _normal_isotope = 167;
      _atomic_mass = static_cast<atomic_mass_t>(167.259);
      _exact_mass     = 165.930290;
      break;

    case 69:
      _symbol = "Tm";
      _normal_isotope = 169;
      _atomic_mass = static_cast<atomic_mass_t>(168.93421);
      _exact_mass     = 168.934212;
      break;

    case 70:
      _symbol = "Yb";
      _normal_isotope = 173;
      _atomic_mass = static_cast<atomic_mass_t>(173.054);
      _exact_mass     = 171.936378;
      break;

    case 71:
      _symbol = "Lu";
      _normal_isotope = 175;
      _atomic_mass = static_cast<atomic_mass_t>(174.9668);
      _exact_mass     = 174.940770;
      break;

    case 72:
      _symbol = "Hf";
      _normal_isotope = 178;
      _atomic_mass = static_cast<atomic_mass_t>(178.49);
      _exact_mass     = 179.9465457;
      break;

    case 73:
      _symbol = "Ta";
      _normal_isotope = 181;
      _atomic_mass = static_cast<atomic_mass_t>(180.9479);
      _exact_mass     = 180.947992;
      break;

    case 74:
      _symbol = "W";
      _normal_isotope = 184;
      _atomic_mass = static_cast<atomic_mass_t>(183.84);
      _exact_mass     = 183.950928;
      break;

    case 75:
      _symbol = "Re";
      _normal_isotope = 186;
      _atomic_mass = static_cast<atomic_mass_t>(186.207);
      _exact_mass     = 186.955744;
      break;

    case 76:
      _symbol = "Os";
      _normal_isotope = 190;
      _atomic_mass = static_cast<atomic_mass_t>(190.23);
      _exact_mass     = 191.961467;
      break;

    case 77:
      _symbol = "Ir";
      _normal_isotope = 192;
      _atomic_mass = static_cast<atomic_mass_t>(192.217);
      _exact_mass     = 192.962917;
      break;

    case 78:
      _symbol = "Pt";
      _normal_isotope = 195;
      _atomic_mass = static_cast<atomic_mass_t>(195.084);
      _exact_mass     = 194.964766;
      break;

    case 79:
      _symbol = "Au";
      _normal_isotope = 197;
      _atomic_mass = static_cast<atomic_mass_t>(196.966569);
      _exact_mass     = 196.966543;
      break;

    case 80:
      _symbol = "Hg";
      _normal_isotope = 201;
      _atomic_mass = static_cast<atomic_mass_t>(200.59);
      _exact_mass     = 201.970617;
      break;

    case 81:
      _symbol = "Tl";
      _normal_isotope = 204;
      _atomic_mass = static_cast<atomic_mass_t>(204.3833);
      _exact_mass     = 204.974401;
      break;

    case 82:
      _symbol = "Pb";
      _normal_isotope = 207;
      _atomic_mass = static_cast<atomic_mass_t>(207.2);
      _exact_mass     = 207.976627;
      break;

    case 83:
      _symbol = "Bi";
      _normal_isotope = 207;
      _atomic_mass = static_cast<atomic_mass_t>(208.98040);
      _exact_mass     = 208.980374;
      break;

    case 84:
      _symbol = "Po";
      _normal_isotope = 210;
      _atomic_mass = static_cast<atomic_mass_t>(209);
      break;

    case 85:
      _symbol = "At";
      _normal_isotope = 210;
      _atomic_mass = static_cast<atomic_mass_t>(209.9871);
      break;

    case 86:
      _symbol = "Rn";
      _normal_isotope = 222;
      _atomic_mass = static_cast<atomic_mass_t>(222.0176);
      break;

    case 87:
      _symbol = "Fr";
      _normal_isotope = 223;
      _atomic_mass = static_cast<atomic_mass_t>(223.0197);
      break;

    case 88:
      _symbol = "Ra";
      _normal_isotope = 226;
      _atomic_mass = static_cast<atomic_mass_t>(226.0254);
      break;

    case 89:
      _symbol = "Ac";
      _normal_isotope = 227;
      _atomic_mass = static_cast<atomic_mass_t>(227.0277);
      break;

    case 90:
      _symbol = "Th";
      _normal_isotope = 232;
      _atomic_mass = static_cast<atomic_mass_t>(232.0381);
      _exact_mass     = 232.0381;
      break;

    case 91:
      _symbol = "Pa";
      _normal_isotope = 231;
      _atomic_mass = static_cast<atomic_mass_t>(231.03588);
      break;

    case 92:
      _symbol = "U";
      _normal_isotope = 238;
      _atomic_mass = static_cast<atomic_mass_t>(238.0289);
      _exact_mass     = 238.0507847;
      break;

    case 93:
      _symbol = "Np";
      _normal_isotope = 237;
      _atomic_mass = static_cast<atomic_mass_t>(237.0);
      break;

    case 94:
      _symbol = "Pu";
      _normal_isotope = 244;
      _atomic_mass = static_cast<atomic_mass_t>(244.0);
      break;

    case 95:
      _symbol = "Am";
      _normal_isotope = 243;
      _atomic_mass = static_cast<atomic_mass_t>(243.0614);
      break;

    case 96:
      _symbol = "Cm";
      _normal_isotope = 247;
      _atomic_mass = static_cast<atomic_mass_t>(247.0704);
      break;

    case 97:
      _symbol = "Bk";
      _normal_isotope = 247;
      _atomic_mass = static_cast<atomic_mass_t>(247.0703);
      break;

    case 98:
      _symbol = "Cf";
      _normal_isotope = 251;
      _atomic_mass = static_cast<atomic_mass_t>(251.0796);
      break;

    case 99:
      _symbol = "Es";
      _normal_isotope = 252;
      _atomic_mass = static_cast<atomic_mass_t>(252);
      break;

    case 100:
      _symbol = "Fm";
      _normal_isotope = 257;
      _atomic_mass = static_cast<atomic_mass_t>(257);
      break;

    case 101:
      _symbol = "Md";
      _normal_isotope = 258;
      _atomic_mass = static_cast<atomic_mass_t>(258);
      break;

    case 102:
      _symbol = "No";
      _normal_isotope = 259;
      _atomic_mass = static_cast<atomic_mass_t>(259);
      break;

    case 103:
      _symbol = "Lr";
      _normal_isotope = 262;
      _atomic_mass = static_cast<atomic_mass_t>(262);
      break;

    case 104:
      _symbol = "Rf";     // Rutherfordium
      _normal_isotope = 267;
      _atomic_mass = static_cast<atomic_mass_t>(267);
      break;

    case 105:
      _symbol = "Db";      // Dubnium
      _normal_isotope = 268;
      _atomic_mass = static_cast<atomic_mass_t>(262.1141);
      break;

    case 106:
      _symbol = "Sg";      //Seaborgium
      _normal_isotope = 271;
      _atomic_mass = static_cast<atomic_mass_t>(271);
      break;

    case 107:
      _symbol = "Bh";     // Borhium
      _normal_isotope = 272;
      _atomic_mass = static_cast<atomic_mass_t>(272);
      break;

    case 108:
      _symbol = "Hs";     // Hassium
      _normal_isotope = 270;
      _atomic_mass = static_cast<atomic_mass_t>(270);
      break;

    case 109:
      _symbol = "Mt";     // Meitnerium
      _normal_isotope = 276;
      _atomic_mass = static_cast<atomic_mass_t>(276);
      break;

    case 110:
      _symbol = "Ds";     // Darmstadtium
      _normal_isotope = 281;
      _atomic_mass = static_cast<atomic_mass_t>(281);
      break;

    case 111:
      _symbol = "Rg";     // Roentgenium
      _normal_isotope = 280;
      _atomic_mass = static_cast<atomic_mass_t>(280);
      break;

    case 112:
      _symbol = "Cn";     // Copernicium
      _normal_isotope = 285;
      _atomic_mass = static_cast<atomic_mass_t>(285);
      break;

    case 113:
      _symbol = "Nh";
      _normal_isotope = 284;
      _atomic_mass = static_cast<atomic_mass_t>(284);
      break;

    case 114:
      _symbol = "Fl";
      _normal_isotope = 289;
      _atomic_mass = static_cast<atomic_mass_t>(289);
      break;

    case 115:
      _symbol = "Mc";
      _normal_isotope = 288;
      _atomic_mass = static_cast<atomic_mass_t>(288);
      break;

    case 116:
      _symbol = "Lv";
      _normal_isotope = 293;
      _atomic_mass = static_cast<atomic_mass_t>(293);
      break;

    case 117:
      _symbol = "Ts";
      _normal_isotope = 294;
      _atomic_mass = static_cast<atomic_mass_t>(294);
      break;

    case 118:
      _symbol = "Og";
      _normal_isotope = 294;
      _atomic_mass = static_cast<atomic_mass_t>(294);
      break;

    case NOT_AN_ELEMENT:
      _atomic_number = NOT_AN_ELEMENT;
      break;

    default:
      _atomic_number = INVALID_ATOMIC_NUMBER;
      cerr << "new_element: element " << i << " cannot be initialised\n";
      iwabort();
  }

  if (0 == _symbol.length())    // happens sometimes when creating non-periodic-table elements
    return;

  (void) _symbol.null_terminated_chars();    // ensure null terminated

  if (0 == _aromatic_symbol.length())
    _aromatic_symbol = _symbol;

  if (_symbol.length() > 2)    // cannot hash these
    return;

  _atomic_symbol_hash_value = element_symbol_hash_function(_symbol.rawchars(), _symbol.length());

//cerr << "Hash value for '" << _symbol << "' is " << _atomic_symbol_hash_value << endl;

// Only write the hash once

  if (NULL == ehash[_atomic_symbol_hash_value])
    ehash[_atomic_symbol_hash_value] = this;

  return;
}

void
Element::_set_symbol (const char * symbol)
{
  if (NULL == symbol)
    _set_symbol(NULL, 0);
  else
    _set_symbol(symbol, static_cast<int>(strlen(symbol)));

  return;
}

void
Element::_set_symbol (const char * symbol, int nchars)
{
  if (NULL == symbol || 0 == nchars)
  {
    _symbol.resize(0);
    _atomic_symbol_hash_value = -1;

    return;
  }

  if (nchars <= 2)   // the most common case
  {
    _symbol.set(symbol, nchars);

    _atomic_symbol_hash_value = element_symbol_hash_function(_symbol.rawchars(), _symbol.length());

    if (_atomic_symbol_hash_value < 0)
      return;

//  cerr << "In set_symbol for '" << _symbol << "' myhashvalue " << _atomic_symbol_hash_value " existing value " << ehash[h] << endl;

    if (NULL == ehash[_atomic_symbol_hash_value])
      ehash[_atomic_symbol_hash_value] = this;

    return;
  }

  if (_atomic_symbols_can_have_arbitrary_length)
  {
    _symbol.set(symbol, nchars);
    _atomic_symbol_hash_value = element_symbol_hash_function_arbitrary_length(symbol, nchars);
    return;
  }

  cerr << "Element::_set_symbol:invalid symbol '";
  cerr.write(symbol, nchars);
  cerr << "'\n";

  _atomic_symbol_hash_value = -1;

  return;
}

void
Element::_set_symbol (const const_IWSubstring & symbol )
{
  _set_symbol(symbol.rawchars(), symbol.nchars());

  return;
}

Element::Element (int i)
{
  _default_values(i);

  return;
}

/*
  Common code for creating elements from different string specifications
*/

void
Element::_non_periodic_table_element_constructor (const char * s,
                               int nchars)
{
  assert (NULL != s);
  assert (nchars > 0);

  _default_values(NOT_AN_ELEMENT);
  _set_symbol(s, nchars);

  if (1 == nchars && islower(*s))
    _needs_square_brackets = 0;
  else
    _needs_square_brackets = 1;

#ifdef DEBUG_NPTEC
  cerr << "Element::_non_periodic_table_element_constructor: created ";
  cerr.write(s, nchars);
  cerr << " hash " << _atomic_symbol_hash_value << endl;
#endif

  return;
}

Element::Element (const char * symbol)
{
  _non_periodic_table_element_constructor(symbol, static_cast<int>(::strlen(symbol)));

  return;
}

Element::Element (const char * symbol, int nchars)
{
  _non_periodic_table_element_constructor(symbol, nchars);

  return;
}

Element::Element (const const_IWSubstring & symbol)
{
  _non_periodic_table_element_constructor(symbol.rawchars(), symbol.length());

  return;
}

/*
  As an element has only simple objects, its destructor is simple.
  Do I even need a destructor here?
*/

Element::~Element()
{
//assert(ok());    very important to leave this out, since some elements are not initialised - gaps in the sequence

  return;
}

int
Element::debug_print (std::ostream & os) const
{
  assert(os.good());

  os << "Element info (" << this << "), atomic number " << _atomic_number << " symbol '" << _symbol << "', unique id " << _unique_id << endl;

  if (! ok())
    os << "Element NOT ok\n";

  return 1;
}

void
debug_print_all_elements (std::ostream & os)
{
  assert (os.good());

  if (0 == elements.number_elements())
  {
    os << "Element array not initialised\n";
    return;
  }

  int nelements = elements.number_elements();

  for (int i = 0; i < nelements; i++)
  {
    os << "Element number " << i << endl;
    elements[i]->debug_print(os);
  }

  return;
}

/*
  Matches name with one of the atomic symbols. We need to do a
  case insensitive match, so we use tmp[4] to hold a case converted
  temporary copy of name. Some systems have stricmp and strcasecmp!!

  Feb 96:
    Modified to accept isotopes as leading digits
*/

const Element *
get_element_from_symbol(const char *name, int nchars, int & isotope)
{
  assert (NULL != name);
  assert (nchars > 0);

// Get any isotopic specification

  isotope = 0;
  while (isdigit(*name) && nchars > 0)
  {
    isotope = 10 * isotope + (*name - '0');
    name++;
    nchars--;
    assert(isotope <= 999);
  }

  if (0 == nchars)
    return NULL;

#ifdef DEBUG_GET_ELEMENT_FROM_SYMBOL
  cerr << "Creating element from '";
  cerr.write(name, nchars);
  cerr << "'\n";
#endif

  char tmp[4];

  if (1 == nchars || '@' == name[1] || 'H' == name[1] || '+' == name[1] || '-' == name[1])    // the most common cases
  {
    tmp[0] = toupper(*name);
    const int hash = element_symbol_hash_function(tmp, 1);
    if (hash < 0 || hash >= SIZE_OF_ELEMENT_HASH_TABLE)
      return NULL;

    return ehash[hash];
  }

// Two or more characters

  if (2 == nchars)
  {
    tmp[0] = toupper(*name);
    tmp[1] = name[1];
    if (isupper(tmp[1]))
      tmp[1] = tolower(tmp[1]);             // change if we ever want to have elements like qQ or QQ

    const int hash = element_symbol_hash_function(tmp, 2);
    if (hash < 0 || hash >= SIZE_OF_ELEMENT_HASH_TABLE)
      return 0;

    return ehash[hash];
  }

  if (_atomic_symbols_can_have_arbitrary_length)
    return get_element_from_long_symbols(name, nchars);

  return NULL;
}

/*
  Jan 2004. Ran into problems with the lowercase elements using get_element_from_symbol, 
  because it converts the first letter to uppercase
*/

//#define DEBUG_GET_ELEMENT_FROM_SYMBOL

const Element *
get_element_from_symbol_no_case_conversion(const char * s,
                                           const int nchars)
{
  if (0 == nchars)
  {
    cerr << "get_element_from_symbol_no_case_conversion:no symbol\n";
    return NULL;
  }

  if (nchars <= 2)
  {
    const int hash = element_symbol_hash_function(s, nchars);

#ifdef DEBUG_GET_ELEMENT_FROM_SYMBOL
    cerr << "Hash value " << hash;
    if (NULL == ehash[hash])
      cerr << ", no element defined\n";
    else
      cerr << " element '" << ehash[hash]->symbol() << endl;
#endif

    if (hash >= 0 && hash < SIZE_OF_ELEMENT_HASH_TABLE)
      return ehash[hash];
  }

  return get_element_from_long_symbols(s, nchars);
}

const Element *
get_element_from_symbol_no_case_conversion(const char * s)
{
  return get_element_from_symbol_no_case_conversion(s, static_cast<int>(::strlen(s)));
}

const Element *
get_element_from_symbol_no_case_conversion(const const_IWSubstring & s)
{
  return get_element_from_symbol_no_case_conversion(s.rawchars(), s.length());
}

const Element *
get_element_from_symbol_no_case_conversion(const IWString & s)
{
  return get_element_from_symbol_no_case_conversion(s.rawchars(), s.length());
}

const Element *
get_element_from_symbol (const char *name, int & isotope)
{
  return get_element_from_symbol(name, static_cast<int>(::strlen(name)), isotope);
}

const Element *
get_element_from_symbol (const const_IWSubstring & name, int & isotope)
{
  return get_element_from_symbol(name.rawchars(), name.nchars(), isotope);
}

const Element *
get_element_from_symbol (const IWString & name, int & isotope)
{
  return get_element_from_symbol(name.rawchars(), name.length(), isotope);
}

/*
  This is very inefficient - we need to malloc space for one character.
  But this function is only called from the smarts parsing routine, so
  who cares?
*/

const Element *
get_element_from_symbol (char name, int & isotope)
{
  IWString tmp = name;

  return get_element_from_symbol(tmp, isotope);
}

int
element_from_long_smiles_string (const char * asymbol,
                                 int nchars,
                                 const Element * & result)
{
//#define DEBUG_ELEMENT_FROM_LONG_SMILES_STRING
#ifdef DEBUG_ELEMENT_FROM_LONG_SMILES_STRING
  cerr << "Getting element '";
  cerr.write(asymbol, nchars);
  cerr << "'\n";
#endif

  int close_square_bracket = -1;

  for (int i = 0; i < nchars; i++)
  {
    if (']' == asymbol[i])
    {
      close_square_bracket = i;
      break;
    }
  }

  if (close_square_bracket < 0)
  {
    cerr << "get_element_from_symbol:no closing square bracket '";
    cerr.write(asymbol, nchars);
    cerr << "'\n";
    return 0;
  }

  if (0 == close_square_bracket)
  {
    cerr << "element_from_long_smiles_string:empty atom specification\n";
    return 0;
  }

// It could be a regular element. This is messy and not robust.

  if (1 == close_square_bracket && isalpha(asymbol[0]))
  {
    char tmp[4];
    tmp[0] = toupper(asymbol[0]);
    if (check_valid_element_via_hash(tmp, 1, result))
      return 1;
  }
  else if (close_square_bracket > 1 && isalpha(asymbol[0]) && ('@' == asymbol[1] || 'H' == asymbol[1] || isdigit(asymbol[1]) || '+' == asymbol[1] || '-' == asymbol[1] || ']' == asymbol[1]))
  {
    char tmp[4];
    tmp[0] = toupper(asymbol[0]);
    if (check_valid_element_via_hash(tmp, 1, result))
      return 1;
  }
  else if (2 == close_square_bracket && isupper(asymbol[0]) && islower(asymbol[1]))
  {
    if (check_valid_element_via_hash(asymbol, 2, result))
      return 2;
  }

#ifdef DEBUG_ELEMENT_FROM_LONG_SMILES_STRING
  cerr << "LINE " << __LINE__ << endl;
#endif

// Not recognised as a regular element

  if ('D' == asymbol[0] && close_square_bracket > 1 && ! isalpha(asymbol[1]) && interpret_d_as_deuterium())
  {
    cerr << "element_from_long_smiles_string:Deuterium detected, but isotope info not returned\n";
    result = get_element_from_atomic_number(1);
    return 1;
  }
  else if ('T' == asymbol[0] && close_square_bracket > 1 && ! isalpha(asymbol[1]) && interpret_t_as_tritium())
  {
    cerr << "element_from_long_smiles_string:Tritium detected, but isotope info not returned\n";
    result = get_element_from_atomic_number(1);
    return 1;
  }

//cerr << "LINE " << __LINE__ << " get_element_from_long_symbols\n";
  result = get_element_from_long_symbols(asymbol, close_square_bracket);
  if (NULL != result)
    return close_square_bracket;

  Element * e = new Element(asymbol, close_square_bracket);

  elements.add(e);
  elements.last_item()->set_unique_id(elements.number_elements());

  result = e;

  return close_square_bracket;
}

//#define DEBUG_ELEMENT_FROM_SMILES_STRING

/*
  Called to parse atomic symbols within square brackets of SMILES (not
  smarts). Because we know we are in square brackets, we know that
  there cannot be things like FP, but there could be nH
*/

int
element_from_smiles_string (const char * smiles,
                            const int nchars,
                            const Element * & result)
{
  assert (NULL != smiles);

  if (_atomic_symbols_can_have_arbitrary_length)
    return element_from_long_smiles_string(smiles, nchars, result);

#ifdef DEBUG_ELEMENT_FROM_SMILES_STRING
  cerr << "Fetching element '" << smiles[0];
  if (nchars > 1)
    cerr << smiles[1];
  cerr << "', nchars = " << nchars << endl;
#endif

  if (! isalpha(smiles[0]))
  {
    cerr << "Invalid element in smiles/smarts '" << smiles[0] << "'\n";
    return 0;
  }

  char ele[3];   // make our own copy so we can capitalise things

  ele[0] = toupper(smiles[0]);

  int nc;
  if (nchars > 1 && islower(smiles[1]))       // Cu for example
  {
    ele[1] = smiles[1];
    ele[2] = '\0';
    nc = 2;
  }
  else                 // U, nH N+ O2-
  {
    ele[1] = '\0';
    nc = 1;
  }

  if (check_valid_element_via_hash(ele, nc, result))
    return nc;

// Feb 2014. What if elements like R0 and A0 have been autocreated

  if (nchars > 1 && isdigit(smiles[1]))
  {
    ele[1] = smiles[1];
    ele[2] = '\0';
    int h = element_symbol_hash_function(ele, 2);
    if (ehash[h])
    {
      result = ehash[h];
      return 2;
    }

    ele[1] = '\0';
    nc = 1;
  }

//cerr << "Still no match, auto create = " << auto_create_new_elements() << " element is '" << ele << "', nc " << nc << endl;

  if (! auto_create_new_elements())
    return 0;    // no match found

// If we recognise D and T, return 0. The caller will then check to see if this
// is a hydrogen isotope. Not very efficient, but too much code to change
// otherwise

  if (1 == nc)
  {
    if ('D' == smiles[0] && interpret_d_as_deuterium())
      return 0;
    if ('T' == smiles[0] && interpret_t_as_tritium())
      return 0;
  }

// We can create any previously unknown elements
    
  result = create_element_with_symbol(ele);

//cerr << "Created element from '" << ele << "'\n";

  return nc;
}

/*
  There is a version for smarts because of things like 'CD' can be valid
  within a smarts.

  NOTE THAT WE NEVER AUTOCREATE ELEMENTS IN A SMARTS - just makes things easier

  Note that we don't need to know the number of characters to process
  because we know the atom is bounded by []
*/

int
element_from_smarts_string (const char * smiles, int characters_to_process,
                            const Element * & result)
{

// Try two character match on natural elements. 2nd char must be lowercase
// Note that we include element '*' as a "natural" element!!

  char ele[3];
  int nchars;

  ele[0] = toupper(smiles[0]);

  if (characters_to_process > 1 && islower(smiles[1]))
  {
    ele[1] = smiles[1];
    ele[2] = '\0';

    nchars = 2;
  }
  else          // single character matches
  {
    if ('D' == ele[0] || 'T' == ele[0])    // nov 97, these are invalid in a smarts
      return 0;

    ele[1] = '\0';
    nchars = 1;
  }

//cerr << "From smarts '" << smiles << "' (" << characters_to_process << " chars), ele is '" << ele << "'\n";

  int hash = element_symbol_hash_function(ele, nchars);

//cerr << "Hash value is " << hash << " ehash " << ehash[hash] << endl;

  if (NULL == ehash[hash])      // no element of that kind yet
    return 0;

  result = ehash[hash];
  return nchars;
}

void
check_elements_magic (void)
{
  int nelements = elements.number_elements();

  for (int i = 0; i < nelements; i++)
  {
    const Element *e = elements[i];
    e->ok();
  }

  return;
}

atomic_number_t
Element::atomic_number() const
{
  return _atomic_number;
}

int
Element::ok() const
{
  if (NOT_AN_ELEMENT == _atomic_number)
    ;
  else if (_atomic_number < 0)
    return 0;
  else if (_atomic_number > HIGHEST_ATOMIC_NUMBER)
    return 0;
  
  if (0 == _symbol.length())
    return 0;

  return 1;
}

const Element *
get_element_from_atomic_number (atomic_number_t z)
{
  if (z >= 0 && z <= HIGHEST_ATOMIC_NUMBER)
    return elements[z];

// Highly unlikely that this next section will work, as generally elements
// above HIGHEST_ATOMIC_NUMBER to not have atomic number values.

  for (int i = HIGHEST_ATOMIC_NUMBER + 1; i < elements.number_elements(); i++)
  {
    const Element * e = elements[i];

    if (z == e->atomic_number())
      return elements[i];
  }

  return NULL;
}

int
Element::append_smiles_symbol (IWString & smiles,
                               aromaticity_type_t arom, 
                               int isotope) const
{
  if (isotope && include_isotopes_in_smiles)
  {
//  cerr << "Adding isotope " << isotope << " to '" << smiles << "'\n";
    smiles.append_number(isotope);
  }

#ifdef DEBUG_APPEND_SMILES_SYMBOL
  if (IS_AROMATIC_ATOM(arom))
    cerr << "Adding '" << _aromatic_symbol << "' to smiles\n";
  else
    cerr << "Adding '" << _symbol << "'\n";
#endif

  if (IS_AROMATIC_ATOM(arom))
    smiles += _aromatic_symbol;
  else
    smiles += _symbol;

#ifdef DEBUG_APPEND_SMILES_SYMBOL
  cerr << "Resulting smiles '" << smiles << "'\n";
#endif

  return 1;
}

atomic_mass_t
Element::atomic_mass() const
{
  assert (ok());

  return _atomic_mass;
}

/*
  Even though this says create an element, we first check to see
  if the element already exists.
*/

const Element *
create_element_with_symbol (const char * symbol, int nchars)
{
#ifdef DEBUG_CREATE_ELEMENT_WITH_SYMBOL
  cerr << "create_element_with_symbol processing '";
  cerr.write(symbol,nchars);
  cerr << endl;
#endif

  if (nchars <= 0)
  {
    cerr << "create_element_with_symbol:no symbol\n";
    return NULL;
  }

  if (nchars <= 2)
    ;
  else if (! auto_create_new_elements())
    return NULL;
  else
  {
    const Element * e = get_element_from_long_symbols(symbol, nchars);
    if (NULL != e)
      return e;

    Element * rc = new Element(symbol, nchars);
    elements.add(rc);
    elements.last_item()->set_unique_id(elements.number_elements());

#ifdef DEBUG_CREATE_ELEMENT_WITH_SYMBOL
    cerr << "create_element_with_symbol: created from ";
    cerr.write(symbol, nchars);
    cerr << " uid " << elements.last_item()->unique_id() << endl;
#endif

    return rc;
  }

  const Element * e = get_element_from_symbol_no_case_conversion(symbol, nchars);
  if (NULL != e)
  {
    cerr << "create_element_with_symbol: cannot create new element with symbol '";
    cerr.write(symbol, nchars) << "', already present\n";

    e->debug_print(cerr);

    return NULL;
  }

  Element * rc = NULL;

  if (nchars > 1 && isupper(symbol[1]))    // symbols cannot have uppercase 2nd characters
  {
    IWString tmp(symbol, nchars);
    tmp[1] = tolower(tmp[1]);

    rc = new Element(tmp.rawchars(), nchars);
  }
  else
    rc = new Element(symbol, nchars);

  elements.add(rc);
  elements.last_item()->set_unique_id(elements.number_elements());

  return rc;
}

const Element *
create_element_with_symbol(const char * symbol)
{
  return create_element_with_symbol(symbol, static_cast<int>(::strlen(symbol)));
}

const Element *
create_element_with_symbol (const const_IWSubstring & symbol)
{
  return create_element_with_symbol(symbol.rawchars(), symbol.nchars());
}

const Element *
create_element_with_symbol (const IWString & symbol)
{
  return create_element_with_symbol(symbol.rawchars(), symbol.nchars());
}

int
Element::alternate_valence (int i) const
{
  assert (ok());

  assert (_alternate_valence.ok_index(i));

  return _alternate_valence[i];
}

int
Element::is_valid_valence (int v) const
{
  assert (ok());

  if (_normal_valence == v)
    return 1;

  return _alternate_valence.contains(v);
}

int
Element::is_halogen() const
{
  if (_atomic_number < 9)
    return 0;

  if (9 == _atomic_number)
    return 1;

  if (17 == _atomic_number)
    return 1;

  if (35 == _atomic_number)
    return 1;

  if (53 == _atomic_number)
    return 1;

  return 0;
}

int
Element::pi_electrons (int ncon, formal_charge_t fc, int & result) const
{
  if (OUTER_SHELL_ELECTRONS_NOT_KNOWN == _outer_shell_electrons)
  {
    result = 0;
    return 0;
  }

  int pe = _outer_shell_electrons - fc;    // a formal charge of -1 means one extra electron

  result = pe - ncon;             // each connection takes one electron.
  if (result > 2)
    result = result - 2;      // lone pairs
  else if (result < 0)
  {
    if (display_strange_chemistry_messages)
    {
      cerr << "Element::pi_electrons: strange chemistry, '" << _symbol << "' " <<
              ncon << " connections";
      if (fc)
        cerr << ',' << fc << " formal charge";
      cerr << ", " << result << " pi electrons (set to 0)\n";
    }

    result = 0;
  }

//cerr << "Atom type " << _atomic_number << " ncon = " << ncon << " fc = " << fc << " pi e's = " << result << endl;

  return 1;
}

int
Element::lone_pairs (int nbonds, formal_charge_t fc, int & result) const
{
  if (OUTER_SHELL_ELECTRONS_NOT_KNOWN == _outer_shell_electrons)
    return 0;

  result = _outer_shell_electrons - fc - nbonds;    // a formal charge of -1 means one extra electron

  if (result < 0)
  {
    if (display_strange_chemistry_messages)
    {
      cerr << "Element::lone_pairs: strange chemistry, '" << _symbol << "' " <<
              nbonds << " bonds";
      if (fc)
        cerr << ',' << fc << " formal charge";
      cerr << ", " << result << " pi electrons (set to 0)\n";
    }

    result = 0;
  }

  result = result / 2;

  return 1;
}

int
Element::outer_shell_electrons (int & result) const
{
  if (OUTER_SHELL_ELECTRONS_NOT_KNOWN == _outer_shell_electrons)
    return 0;

  result = _outer_shell_electrons;

  return 1;
}

/*
  If needed, we can automatically create new element types when encountered
*/

static int automatically_create_new_elements = 0;

int
set_auto_create_new_elements (int i)
{
  return automatically_create_new_elements = i;
}

int
auto_create_new_elements()
{
  return automatically_create_new_elements;
}


static int
read_ptable (const const_IWSubstring & buffer)
{
  int nw = buffer.nwords();
  if (0 == nw)     // ignore blank lines
    return 1;

  const_IWSubstring sym;
  buffer.word(0, sym);

  if (sym.length() > 2)
  {
    cerr << "read_ptable: sorry, can't create elements with length > 2 chars '" << sym << "', ignored\n";
    return 1;
  }

  Element * e = const_cast<Element *>(get_element_from_symbol_no_case_conversion(sym));
  if (NULL == e)
    e  = const_cast<Element *>(create_element_with_symbol(sym));

  if (NULL == e)
  {
    cerr << "read_ptable:sorry, cannot fetch or create element '" << sym << "'\n";
    return 0;
  }

  if (1 == nw)    // just a definition of an element
    return 1;

  return e->read_ptable_record(buffer);
}

static int
read_ptable (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! read_ptable(buffer))
    {
      cerr << "Fatal error processing PTABLE record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
read_ptable_file (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open ptable file '" << fname << "'\n";
    return 0;
  }

  return read_ptable(input);
}
      
int
display_standard_element_options (std::ostream & os)
{
  os << "  -E autocreate  automatically create new elements when encountered\n";
  os << "  -E PTABLE=file use 'file' for element data\n";
  os << "  -E O:Xx        element Xx is considered Organic\n";
#ifdef STUFF_FOR_SAMPLE_ID
  os << "  -E O:SAMID     SampleID makes some metals have valences\n";
#endif
  os << "  -E anylength   elements can be of arbitrary length\n";
  os << "  -E nsqb=xx     element <xx> needs square brackets\n";

  return os.good();
}

/*
  Process all -E options, creating a new element for each.
  If one of the elements is '*', then we set the auto create switch
*/

#include "cmdline.h"

int
process_elements(const Command_Line & cl,
                 int verbose,
                 char eflag)
{
  int j = 0;
  const_IWSubstring c;
  while (cl.value(eflag, c, j++))
  {
    if ('*' == c || "autocreate" == c)
    {
      set_auto_create_new_elements(1);
      continue;
    }

    if (c.starts_with("PTABLE="))
    {
      c.remove_leading_chars(7);
      if (! read_ptable_file(c))
      {
        cerr << "Cannot process ptable directive 'PTABLE=" << c << "'\n";
        return 0;
      }

      continue;
    }

    if (c.starts_with("O:"))
    {
      c.remove_leading_chars(2);

#ifdef STUFF_FOR_SAMPLE_ID
      if ("SAMID" == c)
      {
        make_organic(3);      // lithium
        make_organic(11);     // sodium
        make_organic(12);     // magnesium
        make_organic(14);     // silicon
        make_organic(19);     // potassium
        make_organic(20);     // calcium
        make_organic(34);     // selenium
        continue;
      }
#endif

      Element * e = const_cast<Element *>(get_element_from_symbol_no_case_conversion(c));
      if (NULL == e)
        e = const_cast<Element *>(create_element_with_symbol(c));

      e->set_organic(1);

      continue;
    }

    else if ("anylength" == c)
    {
      _atomic_symbols_can_have_arbitrary_length = 1;
      set_auto_create_new_elements(1);
      continue;
    }

    else if (c.starts_with("nsqb="))
    {
      c.remove_leading_chars(5);

      Element * e = const_cast<Element *>(get_element_from_symbol_no_case_conversion(c));
      if (NULL == e)
        e = const_cast<Element *>(create_element_with_symbol(c));

      e->set_needs_square_brackets(1);
      continue;
    }

    if ("help" == c)
    {
      display_standard_element_options(cerr);
      exit(1);
    }

    const Element * e = create_element_with_symbol(c);

    if (NULL == e)
    {
      cerr << "process_elements: yipes, could not create element '" << c << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Created element with symbol '" << e->symbol() << "'\n";

//#define DEBUG_EXTRA_ELEMENT_CREATION
#ifdef DEBUG_EXTRA_ELEMENT_CREATION
    IWString tmp;
    e->append_smiles_symbol(tmp, NOT_AROMATIC);
    cerr << "Smiles symbol for new element is '" << tmp << "', atomic number " <<
            e->atomic_number() << endl;
#endif
  }

  return 1;
}

int
Element::is_in_periodic_table() const
{
  if (_atomic_number >= 1 && _atomic_number <= HIGHEST_ATOMIC_NUMBER)
    return 1;

  return 0;
}

int
Element::read_ptable_record (const const_IWSubstring & buffer)
{
  int nw = buffer.nwords();

  assert (nw > 1);

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  int nchars = token.length();

  if (nchars > 2)
  {
    cerr << "Element::read_ptable_record: symbol too long '" << token << "'\n";
    return 0;
  }

  _symbol = token;

  if (1 == nw)
    return 1;

// Second token is a comma separated list of valences, with the first one being the normal valence - implement this some time

  (void) buffer.nextword(token, i);

  if (2 == nw)
    return 1;

// Next token is the atomic weight

  (void) buffer.nextword(token, i);

  if (! token.numeric_value(_atomic_mass) || _atomic_mass < 0.0)
  {
    cerr << "Element::read_ptable_record: invalid atomic mass '" << token << "'\n";
    return 0;
  }

  return 1;
}

/*
  Copy some - but not all element data. Remember, this is mostly intended for 
  letting the permanent aromatic elements know that they are actually carbon,
  nitrogen and oxygen atoms

  So, we don't bother copying all the different valences.
  We don't copy _permanent_aromatic
  We definitely don't change the atomic symbol hash value or unique_id
*/

void
Element::copy_element_data (const Element * rhs)
{
  _atomic_number  = rhs->_atomic_number;
  _atomic_mass    = rhs->_atomic_mass;
  _normal_isotope = rhs->_normal_isotope;
  _organic        = rhs->_organic;
  _metal          = rhs->_metal;
  _normal_valence = rhs->_normal_valence;
  _outer_shell_electrons = rhs->_outer_shell_electrons;
  _exact_mass     = rhs->_exact_mass;

  return;
}

void
de_allocate_periodic_table()
{
  elements.resize(0);

  return;
}

/*
  
*/

void
reset_element_file_scope_variables ()
{
  include_isotopes_in_smiles = 1;
  explicit_hydrogens_need_square_brackets_in_smiles = 1;
  _atomic_symbols_can_have_arbitrary_length = 0;
  display_strange_chemistry_messages = 1;
  automatically_create_new_elements = 0;

  for (int i = HIGHEST_ATOMIC_NUMBER+1; i < elements.number_elements(); i++)
  {
    int h = elements[i]->atomic_symbol_hash_value();
    if (h >= 0 && NULL != ehash[h])
    {
//    cerr << "Deleting hash for '" << elements[i]->symbol() << "'\n";
      ehash[h] = NULL;
    }
  }

  elements.resize(HIGHEST_ATOMIC_NUMBER+1);

  return;
}
