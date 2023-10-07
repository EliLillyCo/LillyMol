#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "Foundational/data_source/iwstring_data_source.h"
#include "molecule.h"
#include "rwmolecule.h"

/*
  Here is the file atom.typ, which was used to determine the atom types.

#################################################################
# *** $RCSfile: mmod.cc,v $
# *** $Revision: 1.1 $
# *** $Date: 2007/07/17 12:32:00 $
#################################################################

# atom.typ: information on individual atom types.

# This is a free-format file; fields may be separated by an arbitrary
#  number of blanks, but TABS MAY NOT BE USED. Blank lines are
#  gracefully ignored, as are leading blanks in a line.

# An unused, blank or unspecified field should be signified by an exclamation
#  point (!), which acts as a place-holder.

# Comments start with the character '#' and end with the end of the line.
#  A comment may appear after all valid entries on a line.

# Each line describes an atom type.  The entries on each line are as
#  follows.

# COLUMN 1:  Atom-type number.  This can be any number from 1 up to
#  a compiled-in maximum value, currently 300.  We would like to reserve
#  values up to 200 for expansion of MacroModel "standard" types;  therefore,
#  users wishing to add new types should restrict themselves to the range
#  201 through 300, inclusive.
# COLUMN 2: Atomic number. For a united atom, this is the atomic number of the 
#  underlying heavy atom.  Lone pairs get atomic number -1, dummy atoms -2.
#  The values -3 and 0 are reserved for use by MacroModel; user attempts
#  to specify these may or may not lead to desired effects.  I.e., please
#  avoid them.
# COLUMN 3: Atomic mass, in a.m.u.;  include mass of H's for united types.
# COLUMN 4: A unique two-character name for the type.  Care should be taken
#  not to conflict with either names of other types, or with names used for
#  aliases of atom equivalence classes in the force fields.
# COLUMN 5: The default color for MacroModel display.
# COLUMN 6: Van-der-Waals radius (used only in certain fall-back situations;
#  energetic calculations generally use values from .fld and .slv files).
# N.B.: Columns 7, 8, 9 and 10 are used in formal-charge delocalization.
# COLUMN 7: "T" if formal-charge delocalization is allowed over substituents, 
#  as in the central C in an allyl carbocation; otherwise, "F" or "!".
# COLUMN 8: "T" if, as a substituent atom, this atom type may participate
#  in delocalization of positive formal charges;  "F" or "!" otherwise.
#  Essentially, this column should get a "T" if, when initially single-bonded,
#   this atom type can become positivelyh charged in a resonance structure;  
#   for example, Br in:
#	  H-N+(-H)=C(-H)-Br  <---> H-N(-H)-C(-H)=Br+
# COLUMN 9: Two-letter name of atom type to which this type should be
#  considered equivalent if they are in 1-3 bonded disposition.  The only
#  obvious case is =O ('O2')and O- ('OM'), which are equivalent if they are in, 
#  say, a phosphate or a carboxylate.
# COLUMN 10: Pauling electronegativity of this atom type.
# COLUMN 11: "T" if this atom type can serve as a wild-card type in the 
#  force-field for other types of the same atomic number;  examples: 
#  C0, N0, etc.  Usually the "simplest" of a group of types sharing the 
#  same atomic number.
# COLUMN 12: For any united atom, the two-letter name of the underlying non-
#  united type with the same atomic number and hybridization; for example, all
#  sp3 united-atom carbons should have "C3" in this column.
# COLUMN 13: For any united atom, the number of H's that distinguishes it
#  from the root type given in COLUMN 12.  For example, CA is an
#  sp3 united CH;  it should have "C3" in column 12 and "1" in column 13.
# COLUMN 14: PDB atom name for any atom type that is associated with a unique
#  pdb atom name.  Since a PDB atom name is a 4-char field, we need a way of
#  encoding blanks within the name.  A dot ('.') in the name means a blank.
#  Thus, calcium might be encoded as 'CA..', which adheres to the PDB
#  convention, and distinguishes it from an alpha-carbon, which is '.CA.'.
#  If there are no dots in the name, it will be considered to represent a
#  left-aligned field; that is, 'CA' is equivalent to 'CA..'.
# For "standard" atom types, this field should be unspecified ('!'), since
#  for standard types there is not a one-to-one correspondence between PDB
#  atom names and MMOD atom names.  

#header lines follow:
#     1     2          3     4     5        6      7      8      9     10    11      12    13        14
#at_typ at_no      at_wt  name color  vdw_rad  deloc catdel 1-3eqv el_neg  wild UA_root  no_H  pdb_name

# Carbon, C:
1         6     12.01115    C1     2     1.78      T      F      !   2.50     F      C1     0         !
2         6     12.01115    C2     2     1.72      T      F      !   2.50     F      C2     0         !
3         6     12.01115    C3     2     1.70      F      F      !   2.50     F      C3     0         !
4         6     13.01912    CA     2     1.80      F      F      !   2.50     F      C3     1         !
5         6     14.02709    CB     2     1.90      F      F      !   2.50     F      C3     2         !
6         6     15.03506    CC     2     2.00      F      F      !   2.50     F      C3     3         !
7         6     13.01912    CD     2     1.80      T      F      !   2.50     F      C2     1         !
8         6     14.02709    CE     2     1.80      T      F      !   2.50     F      C2     2         !
9         6     13.01912    CF     2     1.70      T      F      !   2.50     F      C1     1         !
10        6     12.01115    CM     2     1.72      F      F      !   2.50     F       !     !         !
11        6     12.01115    CP     2     1.72      F      T      !   2.50     F       !     !         !
12        6     12.01115    CR     2     1.72      F      F      !   2.50     F       !     !         !
13        6      0.00000    !      2     1.72      F      F      !   2.50     F       !     !         !
14        6     12.01115    C0     2     1.70      F      F      !   2.50     T       !     !         !

# Oxygen, O:
15        8     15.99940    O2    16     1.50      T      T     OM   3.50     F       !     !         !
16        8     15.99940    O3    16     1.52      F      F      !   3.50     F      O3     0         !
17        8     17.00737    OA    16     1.60      F      F      !   3.50     F      O3     1         !
18        8     15.99940    OM    16     1.70      F      F     O2   3.50     F       !     !         !
19        8     18.01534    OW    16     1.80      F      F      !   3.50     F       !     !         !
20        8     15.99940    OP    16     1.52      F      F      !   3.50     F       !     !         !
21        8     15.99940    OQ    16     1.52      F      F      !   3.50     F       !     !         !
22        8      0.00000    !     16     1.52      F      F      !   3.50     F       !     !         !
23        8     15.99940    O0    16     1.60      F      F      !   3.50     T       !     !         !

# Nitrogen, N:
24        7     14.00670    N1     4     1.55      T      T      !   3.00     F      N1     0         !
25        7     14.00670    N2     4     1.55      T      T      !   3.00     F      N2     0         !
26        7     14.00670    N3     4     1.60      F      T      !   3.00     F      N3     0         !
27        7     15.01467    NA     4     1.70      F      T      !   3.00     F      N3     1         !
28        7     16.02264    NB     4     1.75      F      T      !   3.00     F      N3     2         !
29        7     15.01467    NC     4     1.70      T      T      !   3.00     F      N2     1         !
30        7     16.02264    ND     4     1.75      T      T      !   3.00     F      N2     2         !
31        7     14.00670    N4     4     1.55      T      F      !   3.00     F      N4     0         !
32        7     14.00670    N5     4     1.55      F      F      !   3.00     F      N5     0         !
33        7     15.01467    NE     4     1.60      F      F      !   3.00     F      N5     1         !
34        7     16.02264    NF     4     1.70      F      F      !   3.00     F      N5     2         !
35        7     17.03061    NG     4     1.80      F      F      !   3.00     F      N5     3         !
36        7     15.01467    NH     4     1.60      T      F      !   3.00     F      N4     1         !
37        7     16.02264    NI     4     1.70      T      F      !   3.00     F      N4     2         !
38        7            !    !      4        !      F      F      !   3.00     F       !     !         !
39        7            !    !      4        !      F      F      !   3.00     F       !     !         !
40        7     14.00670    N0     4     1.55      F      F      !   3.00     T       !     !         !

# Hydrogen, H:
41        1      1.00797    H1    21     1.20      F      F      !   2.00     F       !     !         !
42        1      1.00797    H2    17     1.00      F      F      !   2.00     F       !     !         !
43        1      1.00797    H3     5     1.10      F      F      !   2.00     F       !     !         !
44        1      1.00797    H4     5     1.20      F      F      !   2.00     F       !     !         !
45        1      1.00797    H5     5     1.20      F      F      !   2.00     F       !     !         !
46        1      1.00797    !      5        !      F      F      !   2.00     F       !     !         !
47        1      1.00797    !      5        !      F      F      !   2.00     F       !     !         !
48        1      1.00797    H0    21     1.20      F      F      !   2.00     T       !     !         !

# Sulfur, S:
49       16     32.06400    S1    13     1.80      T      T      !   2.50     F      S1     0         !
50       16     32.06400    SA    13     2.00      T      T      !   2.50     F      S1     1         !
51       16     32.06400    SM    13     1.80      F      F      !   2.50     F       !     !         !
52       16     32.06400    S0    13     1.80      T      F      !   2.50     T       !     !         !

# Phosphorus, P:
53       15     30.97380    P0     8     1.80      T      F      !   2.10     F       !     !         !

# Boron: B
54        5     11.00090    B2    10     1.75      F      F      !   0.00     F       !     !         !
55        5     11.00090    B3    10     1.75      F      F      !   0.00     F       !     !         !

# Halogens:
56        9     18.99840    F0    16     1.47      T      T      !   4.00     F       !     !         !
57       17     35.45300    Cl    13     1.75      T      T      !   3.00     F       !     !         !
58       35     79.90900    Br     8     1.85      T      T      !   2.80     F       !     !         !
59       53    126.90440    I0    19     1.98      T      T      !   2.50     F       !     !         !

# Silicon: Si:
60       14     28.08600    Si     4     2.10      T      F      !   1.80     F       !     !         !

# Dummy atom (for FEP): Du;  these are given atomic number -2:
61       -2     12.00000    Du    10     1.00      F      F      !   0.00     F       !     !         !

# Atom to be defined: Z0; perhaps no longer needed with increased atm. types:
62        !    100.00000    Z0     4     1.75      T      F      !   1.50     F       !     !         !

# Lone pair: Lp;  these get atomic number -1:
63       -1     12.00000    Lp    16     0.50      F      F      !   3.50     F       !     !         !

# Any atom: 00;  used as a wildcard in the force-field:
64        !     12.00000    00     4     0.00      T      F      !   0.00     F       !     !         !

# Metal ions;  vdw_rad are radii from .slv files; el_neg are from Levine, P.Chem., p655, but
#  are wrong, since the values refer to the elements, not the ions.  But they'll probably
#  never be used in the program.  Cs and Ba el_neg's are guesses.:
# Alkali:
65        3      6.94100    Li     4     1.43      F      F      !   1.00     F       !     !      LI..  # charge +1
66       11     22.99000    Na     4     1.92      F      F      !   0.90     F       !     !      NA..  # charge +1
67       19     39.09800    K0     4     2.12      F      F      !   0.80     F       !     !      K...  # charge +1
68       37     85.47000    Rb     4     2.26      F      F      !   0.80     F       !     !      RB..  # charge +1
69       55    132.91000    Cs     4     2.51      F      F      !   0.70     F       !     !      CS..  # charge +1

# Alkaline earth:
70       20     40.80000    Ca     4     1.73      F      F      !   1.00     F       !     !      CA..  # charge +2
71       56    137.34000    Ba     4     2.08      F      F      !   0.80     F       !     !      BA..  # charge +2

*/
static atomic_number_t
convert_from_mmod_atom_type(int mmod_type)
{
  if (mmod_type <= 0)
    return kInvalidAtomicNumber;

  if (mmod_type <= 14)
    return 6;

  if (mmod_type <= 23)
    return 8;

  if (mmod_type <= 40)
    return 7;

  if (mmod_type <= 48)
    return 1;

  if (mmod_type <= 52)
    return 16;

  if (mmod_type <= 53)
    return 15;

  if (mmod_type <= 55)
    return 5;

  if (mmod_type <= 56)
    return 9;

  if (mmod_type <= 57)
    return 17;

  if (mmod_type <= 58)
    return 35;

  if (mmod_type <= 59)
    return 53;

  if (mmod_type <= 60)
    return 14;

  if (65 == mmod_type)
    return 3;

  if (66 == mmod_type)
    return 11;

  if (67 == mmod_type)
    return 19;

  if (68 == mmod_type)
    return 37;

  if (69 == mmod_type)
    return 55;

  return kInvalidAtomicNumber;
}

int
Molecule::read_molecule_mmod_ds(iwstring_data_source& input)
{
  assert(ok());
  assert(input.good());

  IWString buffer;
  if (! input.next_record(buffer))
    return 0;

  set_name(buffer);

  int atoms_read = 0;
  while (input.next_record(buffer))
  {
    int atom_type;
    atom_number_t c1, c2, c3, c4, c5, c6;
    bond_type_t b1, b2, b3, b4, b5, b6;
    coord_t x, y, z;
    int tokens =
        IW_SSCANF(buffer.chars(), "%d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f", &atom_type,
                  &c1, &b1, &c2, &b2, &c3, &b3, &c4, &b4, &c5, &b5, &c6, &b6, &x, &y, &z);
    if (16 != tokens)
      return rwmolecule_error("read molecule mmod ds: incorrect token count", input);

    atomic_number_t zz = convert_from_mmod_atom_type(atom_type);
    if (kInvalidAtomicNumber == zz)
      return rwmolecule_error("read molecule mmod ds: Unknown atom type", input);

    Atom* a = new Atom(zz);
    assert(nullptr != a);
    a->setxyz(x, y, z);
    add(a);
    //  cerr << "Doing bonds " << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << " " << c6 << "\n";
    //  cerr << "atoms read = " << atoms_read << "\n";

    if (c1 > 0 && c1 > atoms_read)
      add_bond(atoms_read, c1 - 1, b1, 1);
    if (c2 > 0 && c2 > atoms_read)
      add_bond(atoms_read, c2 - 1, b2, 1);
    if (c3 > 0 && c3 > atoms_read)
      add_bond(atoms_read, c3 - 1, b3, 1);
    if (c4 > 0 && c4 > atoms_read)
      add_bond(atoms_read, c4 - 1, b4, 1);
    if (c5 > 0 && c5 > atoms_read)
      add_bond(atoms_read, c5 - 1, b5, 1);
    if (c6 > 0 && c6 > atoms_read)
      add_bond(atoms_read, c6 - 1, b6, 1);
    atoms_read++;
  }

  finished_bond_addition();
  check_bonding();

  if (natoms())
    return 1;
  else
    return 0;
}

// Not actually implemented.
#ifdef IMPLEMENT_THIS_SOMETIME
static int
macromodel_atom_type(const Molecule* m, atom_number_t i)
{
  assert(OK_ATOM_NUMBER(m, i));

  //atomic_number_t zi = m->atomic_number (i);

  return 0;
}
#endif

int
Molecule::write_molecule_mmod(std::ostream& output) const
{
#ifdef IMPLEMENT_THIS_SOMETIME
  assert(ok());
  if (! output.good())
    return 0;

  output << name() << "\n";

  int matoms = natoms();
  for (int i = 0; i < matoms && output.good(); i++)
  {
    int j = macromodel_atom_type (this, i);
  }

  return output.good() ? 1 : 0;
#endif
  return 1;
}
