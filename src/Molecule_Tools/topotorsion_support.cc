#include <stdlib.h>
#include <iostream>

#include <limits>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "topotorsion_support.h"

#include "Foundational/iwmisc/iwdigits.h"

using std::cerr;
using std::endl;

static IWDigits iwdigits(10);

TopoTorsion::TopoTorsion ()
{
  IWString::resize (12);

  _count = 1;

  return;
}

template class resizable_array_p<TopoTorsion>;
template class resizable_array_base<TopoTorsion *>;

TopoTorsion::TopoTorsion (const char * s) : IWString (s)
{
  _count = 1;

  return;
}

static int
topological_torsion_comparitor (IWString * const * pps1, IWString * const * pps2)
{
  const IWString * s1 = *pps1;
  const IWString * s2 = *pps2;

  if (s1->length () == s2->length ())
    return ::strncmp (s1->rawchars (), s2->rawchars (), s1->length ());

  if (s1->length () < s2->length ())
    return -1;

  return 1;
}

#define TTQSSORT ( int (*) (TopoTorsion * const *, TopoTorsion * const *))

int
sort_and_count_tts (resizable_array_p<TopoTorsion> & tt)
{
  tt.sort (TTQSSORT topological_torsion_comparitor);

  int nt = tt.number_elements ();

  TopoTorsion * prev = tt[0];
  for (int i = 1; i < nt; i++)
  {
    TopoTorsion * t = tt[i];
    if (*t != *prev)
      prev = t;
    else
    {
      tt.remove_item (i);
      i--;
      nt--;
      prev->increment_count ();
    }
  }

  return 1;
}

/*
  Any single character atomic symbol is used in eight character torsions.
  Some two character elements are transformed.
*/

int
initialise_single_character_atom_symbols (const Molecule & m,
                                          IWString & asymbol)
{
  int matoms = m.natoms ();

  asymbol.extend (matoms, '_');

  for (int i = 0; i < matoms; i++)
  {
    const_IWSubstring as = m.atomic_symbol (i);

    if (1 == as.length ())
      asymbol[i] = as[0];
    else if ("Cl" == as)
      asymbol[i] = 'X';
    else if ("Br" == as)
      asymbol[i] = 'D';
    else if ("Se" == as)
      asymbol[i] = 'E';
    else
      cerr << "Unusual element '" << as << "' using wildcard\n";
  }

  return 1;
}


static void
invert (atom_number_t & a0,
        atom_number_t & a1,
        atom_number_t & a2,
        atom_number_t & a3)
{
  atomic_number_t tmp = a0;
  a0 = a3;
  a3 = tmp;

  tmp = a1;
  a1 = a2;
  a2 = tmp;

  return;
}

/*
  When forming the topological torsion, we must canonicalise the order
*/

void
form_canonical_order (Molecule & m,
                      const atomic_number_t * z,
                      const int * ncon,
                      atom_number_t & a0,
                      atom_number_t & a1,
                      atom_number_t & a2,
                      atom_number_t & a3)
{
  if (z[a1] < z[a2])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (z[a1] > z[a2])
    return;

  if (ncon[a1] < ncon[a2])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (ncon[a1] > ncon[a2])
    return;

  int pe1 = -1;
  (void) m.pi_electrons (a1, pe1);

  int pe2 = -1;
  (void) m.pi_electrons (a2, pe2);

  if (pe1 < pe2)
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (pe1 > pe2)
    return;

// Not resolved with a1 and a2, try a0 and a3

  if (z[a0] < z[a3])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (z[a0] > z[a3])
    return;

  if (ncon[a0] < ncon[a3])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (ncon[a0] > ncon[a3])
    return;

  int pe0 = -1;
  (void) m.pi_electrons (a0, pe0);

  int pe3 = -1;
  (void) m.pi_electrons (a3, pe3);

  if (pe0 < pe3)
  {
    invert (a0, a1, a2, a3);
    return;
  }

  return;
}

void
form_canonical_order (const int * invariant,
                      atom_number_t & a0,
                      atom_number_t & a1,
                      atom_number_t & a2,
                      atom_number_t & a3)
{
  if (invariant[a1] < invariant[a2])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (invariant[a1] > invariant[a2])
    return;

// Not resolved with a1 and a2, try a0 and a3

  if (invariant[a0] < invariant[a3])
  {
    invert (a0, a1, a2, a3);
    return;
  }

  if (invariant[a0] > invariant[a3])
    return;

  return;
}

/*
  We need to hash the ncon and the pi value to a unique number. 
*/

#define MAXNCP 62

static char eightcharcode[MAXNCP] = 
{
  '0',
  '1',
  '2',
  '3',
  '4',
  '5',
  '6',
  '7',
  '8',
  '9',
  'a',
  'b',
  'c',
  'd',
  'e',
  'f',
  'g',
  'h',
  'i',
  'j',
  'k',
  'l',
  'm',
  'n',
  'o',
  'p',
  'q',
  'r',
  's',
  't',
  'u',
  'v',
  'w',
  'x',
  'y',
  'z',
  'A',
  'B',
  'C',
  'D',
  'E',
  'F',
  'G',
  'H',
  'I',
  'J',
  'K',
  'L',
  'M',
  'N',
  'O',
  'P',
  'Q',
  'R',
  'S',
  'T',
  'U',
  'V',
  'W',
  'X',
  'Y',
  'Z'
};

void
append_numbers_eight_chars (IWString & tt,
             const Molecule & m,
             const IWString & asymbol,
             atom_number_t a,
             int pi,
             int ncon)
{
  tt += asymbol[a];

  ncon--;
  assert (ncon >= 0);

  assert (pi >= 0);

  int p = ncon * 8 + pi;

  if (p > MAXNCP)
  {
    cerr << "Yipes, cannot encode ncon " << ncon << " and pi " << pi << " value is " << p << 
            " max is " << MAXNCP << endl;
    tt << '*';
  }
  else
    tt << eightcharcode[p];

  return;
}


void
append_numbers_twelve_chars (IWString & tt,
             const Molecule & m,
             atom_number_t a,
             int pi,
             int ncon)
{
  const IWString & s = m.atomic_symbol (a);

  if (1 == s.length ())
    tt.add (s[0]);
  else if ("Cl" == s)
    tt.add ('X');
  else if ("Br" == s)
    tt.add ('Y');
  else
    tt.add ('Z');

  iwdigits.append_number (tt, ncon);
  iwdigits.append_number (tt, pi);

  return;
}

/*
  Make sure we use letters first as SAS requires variable names to start
  with a letter
*/

static char thirtysixchars [] = {
  'A',
  'B',
  'C',
  'D',
  'E',
  'F',
  'G',
  'H',
  'I',
  'J',
  'K',
  'L',
  'M',
  'N',
  'O',
  'P',
  'Q',
  'R',
  'S',
  'T',
  'U',
  'V',
  'W',
  'X',
  'Y',
  'Z',
  '0',
  '1',
  '2',
  '3',
  '4',
  '5',
  '6',
  '7',
  '8',
  '9'
};

/*
  We need a means of packing atomic numbers, pi and ncon values into a number
  less than 99. Each atomic number gets a value
*/

static int atomic_number_code [HIGHEST_ATOMIC_NUMBER + 1] =
{
  9,  //       atomic number 0
  9,  //       H
  9,  //       He
  9,  //  
  9,  //  
  9,  //  
  0,  //      Carbon
  1,  //      Nitrogen
  2,  //      Oxygen
  3,  //      Fluorine
  9,  //       Ne               10
  9,  //       Na
  9,  //       Mg
  9,  //       Al
  9,  //       Si
  4,  //      P
  5,  //      S
  6,  //      Cl
  9,  //       Ar
  9,  //       K
  9,  //                        20
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                        30
  9,  //      Ga
  9,  //      Ge
  9,  //      As
  9,  //      Se
  7,  //     Br
  9,  //       Kr
  9,  //       Rb
  9,  //      Sr
  9,  //      Y
  9,  //      Zr                40
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                        50
  9,  //  
  9,  //  
  8,  //     I
  9,  //      Xe
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                        60
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                         70
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                         80
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                         90
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //                         100
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9,  //  
  9
};

/*
  We set this up so that the first character in the torsion is never numeric.
*/

void
append_numeric_torsion_names (IWString & tt,
             const Molecule & m,
             atom_number_t a,
             int pi,
             int ncon)
{
  atomic_number_t z = m.atomic_number (a);

  int v;
  if (NOT_AN_ELEMENT == z)
    v = 9;
  else
    v = atomic_number_code[z];

  ncon--;

  int rc = 28 * v + 9 * ncon + pi;

  int d1 = rc / 36;
  int d2 = rc % 36;

  tt += thirtysixchars[d1];
  tt += thirtysixchars[d2];

  return;
}

unsigned char
TopoTorsion::count_as_unsigned_char () const
{
  if (_count <= std::numeric_limits<unsigned char>::max())
    return static_cast<unsigned char> (_count);

  return std::numeric_limits<unsigned char>::max();
}
