#ifndef _IW_BASIC_MOLECULE_TYPES
#define _IW_BASIC_MOLECULE_TYPES

/*
  Various typedef's and such used for molecules and others.
*/

#ifdef UNIX
typedef int magic_number_t;
#else
typedef long magic_number_t;
#endif

//typedef long long          int64_t;
//typedef long unsigned long uint64_t;
typedef unsigned long long iw_uint64_t;

//typedef unsigned int uint32_t;

typedef float coord_t;

typedef float charge_t;
typedef int formal_charge_t;

typedef float energy_t;

typedef double area_t;
typedef double volume_t;

//#define REASONABLE_FORMAL_CHARGE(q) ( (q) >= -7 && (q) <= 7 )

typedef int atomic_number_t;     // proton count

typedef float atomic_mass_t;     // float to include isotopic weighting
typedef double exact_mass_t;     // these have lots of digits
typedef float molecular_weight_t;

#define INVALID_ATOMIC_NUMBER -1

#define INVALID_ATOM_TYPE -8

#define SINGLE_BOND_SYMBOL '-'
#define DOUBLE_BOND_SYMBOL '='
#define TRIPLE_BOND_SYMBOL '#'
#define AROMATIC_BOND_SYMBOL ':'

typedef float angle_t;

typedef float distance_t;

typedef float bond_length_t;

typedef int atom_number_t;     // number of each atom within a molecule, starts at 0 

// Value for atom numbers which are invalid.

#define INVALID_ATOM_NUMBER -1

/*
  Bond types are stored as by setting bits in a word.
  The reason for this is for aromaticity, so we can store not only whether
  a bond is single or double, but also whether it is aromatic.

  Bits are allocated according to the numbers below.

  A bond can be single, or single and aromatic.
  A bond can be double, or double and aromatic.

  Consult substructure.h before changing any of these!
*/

typedef unsigned int bond_type_t;

#define UNKNOWN_BOND_TYPE 0
#define SINGLE_BOND 1
#define DOUBLE_BOND 2
#define TRIPLE_BOND 4
#define AROMATIC_BOND 8
#define NOT_AROMATIC_BOND 0xfff7

#define ATOMS_NOT_BONDED -1

/*
  Jan 2004. The pre-existing aromatic concept is really a computed or perceived
  aromaticity. Now, we have the concept of an atom or bond that is purely
  aromatic, never anything else. We may or may not also perceive Kekule
  forms for such a bond
*/

#define PERMANENT_AROMATIC_BOND 16

#define NOT_A_BOND 32
#define INVALID_BOND_TYPE 64

#define COORDINATION_BOND 128

#define OK_BOND_TYPE(b) ((b) <= DOUBLE_BOND || \
                         (b) == TRIPLE_BOND || \
                         (IS_DOUBLE_BOND (b) && IS_AROMATIC_BOND (b)) || \
                         (IS_SINGLE_BOND (b) && IS_AROMATIC_BOND (b)) || \
                         (b) == NOT_A_BOND || IS_PERMANENT_AROMATIC_BOND (b) ||\
                         IS_AROMATIC_BOND(b) || \
                         COORDINATION_BOND==(b))

#define IS_SINGLE_BOND(b) ((b) & SINGLE_BOND)
#define IS_DOUBLE_BOND(b) ((b) & DOUBLE_BOND)
#define IS_TRIPLE_BOND(b) ((b) & TRIPLE_BOND)
#define IS_COMPUTED_AROMATIC_BOND(b) ((b) & AROMATIC_BOND)
#define IS_PERMANENT_AROMATIC_BOND(b) ((b) & PERMANENT_AROMATIC_BOND)
#define IS_AROMATIC_BOND(b) (((b) & AROMATIC_BOND) || ((b) & PERMANENT_AROMATIC_BOND))
#define IS_MULTIPLE_BOND(b) ((b) & DOUBLE_BOND || (b) & TRIPLE_BOND)
#define IS_NOT_A_BOND(b) (NOT_A_BOND == (b))
#define IS_INVALID_BOND_TYPE(b) (INVALID_BOND_TYPE == (b))
#define IS_COORDINATION_BOND(b) (COORDINATION_BOND & (b))

#define SET_AROMATIC_BOND(b) ((b) |= AROMATIC_BOND)                            
#define SET_NON_AROMATIC_BOND(b) ((b) &= NOT_AROMATIC_BOND)

/*
  Sometimes we need just type bond's type.
*/

#define BOND_TYPE_ONLY(b) (((b) & SINGLE_BOND) | ((b) & DOUBLE_BOND) | ((b) & TRIPLE_BOND))

#define BOND_TYPE_ONLY_MASK (SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND)

/*
  We can also indicate whether or not this is a ring bond
  Really, the only reason for putting this in the _btype is because
  of the way substructure searching works. Currently it works on the
  bond_type only. If we ever change that, we can remove this from
  _btype
*/

//#define RING_BOND 512
//#define IS_RING_BOND(qbb) (0 != ((qbb) & RING_BOND))
//#define IS_NON_RING_BOND(qbb) (0 == ((qbb) & RING_BOND))
//#define SET_RING_BOND(qbb) ((qbb) = (qbb) | RING_BOND)
//#define SET_NON_RING_BOND(qbb) ((qbb) = (qbb) ^ RING_BOND)


/*
  In ok_atom, we do impose a maximum connectivity. This is only a
  chemical reasonabless thing, and not a limitation of the programme.
  The programme imposes no maximum connectivity.
*/

#define MAX_CONNECTIVITY 8

#define DELETE_IF_NOT_NULL(p) if (NULL != (p)) { delete (p); (p) = NULL; }
#define DELETE_IF_NOT_NULL_ARRAY(p) if (NULL != (p)) { delete [] (p); (p) = NULL; }

typedef int boolean;

typedef int seed;      /* for random number generator random */

#ifdef UNIX
typedef int int32;
#define IW_BYTES_PER_INT 4
#define IW_BITS_PER_BYTE 8
#define IW_BITS_PER_INT (IW_BYTES_PER_INT * IW_BITS_PER_BYTE)
#else
typedef long int32;
#endif

#include <math.h>

#if ! defined (IWCYGWIN)
//#include <values.h>
#endif

#ifndef M_PI
#define M_PI 3.141592654
#endif

#define HALF_PI (M_PI * 0.50)

#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

#define NRINGS_UNKNOWN -91
#define IS_NON_RING_ATOM -74
#define IS_RING_ATOM  -73

/*
  After running into problems with fused ring systems (specifically
  where an aromatic ring is fused to a non aromatic system), I decided
  to add more complexity to aromaticity.

  An atom can be both aromatic and non-aromatic.
  
  All this is done by setting bits...
*/

// Not sure if this is a good idea or a bad idea, but we maintain consistency with
// the bit allocation in bond_type_t

#define _AROM_BIT 8
#define _ALIPH_BIT 128

typedef int aromaticity_type_t;

#define AROMATICITY_NOT_DETERMINED 0
#define NOT_AROMATIC _ALIPH_BIT
#define AROMATIC     _AROM_BIT

#define IS_AROMATIC_ATOM(c) (_AROM_BIT & (c))
#define IS_ALIPHATIC_ATOM(c) (_ALIPH_BIT & (c))
#define IS_AROMATIC_ONLY_ATOM(c) (_AROM_BIT == (_AROM_BIT & (c)))
#define IS_ALIPHATIC_ONLY_ATOM(c) (_ALIPH_BIT == (_ALIPH_BIT & (c)))
#define IS_AROMATIC_AND_ALIPHATIC_ATOM(c) ((_AROM_BIT & (c)) && (_ALIPH_BIT & (c)))

#define SET_AROMATIC_ATOM(c) ((c) = ((c) | _AROM_BIT))
#define SET_ALIPHATIC_ATOM(c) ((c) = ((c) |_ALIPH_BIT))

#define OK_ATOM_AROMATICITY(c) (((_AROM_BIT | _ALIPH_BIT) == (_AROM_BIT | _ALIPH_BIT | (c))) || \
                           (AROMATICITY_NOT_DETERMINED == (c)))

typedef float similarity_type_t;

/*
  We can run into problems in single precision if a value for arcos
  gets outside the correct range. This macro can fix that, and
  put the value back into [-1.0,1.0]
*/

#define OK_ACOS(x) ((x) > 1.0 ? 1.0 : (x) < -1.0 ? -1.0 : (x))

typedef unsigned int atom_invariant_type_t;

/*
  I often needed to identity the various file types.
*/

#define MDL    1
#define PDB    2
#define MMOD   4
#define SMI    5
#define USMI   6
#define MSI    7
#define TDT    8
#define RDF    9
#define QRY    10
/*#define BFILE  11*/
#define RSMI   12
#define UTDT   13
#define SDF    14
#define MOL2   15
#define IWMTYPE_PSF 16
#define IWMTYPE_CRD 17
#define IWMTYPE_CHM 18
#define IWMTYPE_MRK 19
#define IWMTYPE_WCHM 21

/*
  Karypis Graph Transaction file
*/

//#define IWMTYPE_GTF 22

// Unique smiles order, but without aromaticity

#define IWMTYPE_NAUSMI 20

// smarts are only written

#define IWMTYPE_SMT 23

#define IWMTYPE_MRV 24

#define IWMTYPE_INCHI 25

#define IWMTYPE_TDT_NAUSMI 26

#define IWMTYPE_CIF 27

/*enum 
{
  FILE_TYPE_INVALID,
  FILE_TYPE_MDL,
  FILE_TYPE_PDB,
  FILE_TYPE_MMOD,
  FILE_TYPE_SMI,
  FILE_TYPE_USMI,
  FILE_TYPE_MSI,
  FILE_TYPE_TDT,
  FILE_TYPE_RDF,
  FILE_TYPE_QRY,
  FILE_TYPE_BFILE,
} file_type_t;*/

typedef int chiral_type_t;

#define CHIRALITY_NOT_DETERMINED -99
#define NON_CHIRAL 0
#define CHIRAL_BUT_UNSPECIFIED 1
#define CHIRAL_R   2
#define CHIRAL_L   3

#define ATOM_TYPE_SYBYL "SYBYL"

#include "iwconfig.h"

#define NOT_COMPUTED -117

#endif
