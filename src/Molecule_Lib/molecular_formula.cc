
#include <iostream>

#include "Molecule_Lib/molecular_formula.h"

namespace molecular_formula {

constexpr char kOpenParen = '(';
constexpr char kCloseParen = '(';
constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';
constexpr char kSingleBond = '-';


// For each character that is ignored, set an entry in this array.
static int ignore[256] = {
  1,   // 0  
  1,   // 1  
  1,   // 2  
  1,   // 3  
  1,   // 4  
  1,   // 5  
  1,   // 6  
  1,   // 7  
  1,   // 8  
  1,   // 9  
  1,   // 10  
  1,   // 11  
  1,   // 12  
  1,   // 13  
  1,   // 14  
  1,   // 15  
  1,   // 16  
  1,   // 17  
  1,   // 18  
  1,   // 19  
  1,   // 20  
  1,   // 21  
  1,   // 22  
  1,   // 23  
  1,   // 24  
  1,   // 25  
  1,   // 26  
  1,   // 27  
  1,   // 28  
  1,   // 29  
  1,   // 30  
  1,   // 31  
  1,   // 32  
  1,   // 33 !
  1,   // 34 "
  1,   // 35 #
  1,   // 36 $
  1,   // 37 %
  1,   // 38 &
  1,   // 39 '
  1,   // 40 (
  1,   // 41 )
  1,   // 42 *
  1,   // 43 +
  1,   // 44 ,
  1,   // 45 -
  1,   // 46 .
  1,   // 47 /
  1,   // 48 0
  1,   // 49 1
  1,   // 50 2
  1,   // 51 3
  1,   // 52 4
  1,   // 53 5
  1,   // 54 6
  1,   // 55 7
  1,   // 56 8
  1,   // 57 9
  1,   // 58 :
  1,   // 59 ;
  1,   // 60 <
  1,   // 61 =
  1,   // 62 >
  1,   // 63 ?
  1,   // 64 @
  1,   // 65 A
  0,   // 66 B
  0,   // 67 C
  1,   // 68 D
  1,   // 69 E
  0,   // 70 F
  1,   // 71 G
  1,   // 72 H
  0,   // 73 I
  1,   // 74 J
  1,   // 75 K
  1,   // 76 L
  1,   // 77 M
  0,   // 78 N
  0,   // 79 O
  0,   // 80 P
  1,   // 81 Q
  1,   // 82 R
  0,   // 83 S
  1,   // 84 T
  1,   // 85 U
  1,   // 86 V
  1,   // 87 W
  1,   // 88 X
  1,   // 89 Y
  1,   // 90 Z
  0,   // 91 [
  1,   // 92 \
  0,   // 93 ]
  1,   // 94 ^
  1,   // 95 _
  1,   // 96 `
  1,   // 97 a
  1,   // 98 b
#ifdef NOT_SURE_WHAT_IS_GOING_ON
  I build this with ruby
256.times do |i|
  q = if i.chr.to_str =~ /[[:print:]]/
        i.chr.to_str
      else
        ' '
      end

  $stdout << "  0,   // #{i} #{q}\n"
end
But for some reason, the lowercase letters are offset. No idea
what went wrong, but ignoring it for now.
#endif
  1,   // 99 c
  0,   // 100 d
  1,   // 101 e
  1,   // 102 f
  1,   // 103 g
  1,   // 104 h
  1,   // 105 i
  1,   // 106 j
  1,   // 107 k
  1,   // 108 l
  1,   // 109 m
  1,   // 110 n
  0,   // 111 o
  0,   // 112 p
  0,   // 113 q
  1,   // 114 r
  1,   // 115 s
  0,   // 116 t
  1,   // 117 u
  1,   // 118 v
  1,   // 119 w
  1,   // 120 x
  1,   // 121 y
  1,   // 122 z
  1,   // 123 {
  1,   // 124 |
  1,   // 125 }
  1,   // 126 ~
  1,   // 127  
  1,   // 128  
  1,   // 129  
  1,   // 130  
  1,   // 131  
  1,   // 132  
  1,   // 133  
  1,   // 134  
  1,   // 135  
  1,   // 136  
  1,   // 137  
  1,   // 138  
  1,   // 139  
  1,   // 140  
  1,   // 141  
  1,   // 142  
  1,   // 143  
  1,   // 144  
  1,   // 145  
  1,   // 146  
  1,   // 147  
  1,   // 148  
  1,   // 149  
  1,   // 150  
  1,   // 151  
  1,   // 152  
  1,   // 153  
  1,   // 154  
  1,   // 155  
  1,   // 156  
  1,   // 157  
  1,   // 158  
  1,   // 159  
  1,   // 160  
  1,   // 161  
  1,   // 162  
  1,   // 163  
  1,   // 164  
  1,   // 165  
  1,   // 166  
  1,   // 167  
  1,   // 168  
  1,   // 169  
  1,   // 170  
  1,   // 171  
  1,   // 172  
  1,   // 173  
  1,   // 174  
  1,   // 175  
  1,   // 176  
  1,   // 177  
  1,   // 178  
  1,   // 179  
  1,   // 180  
  1,   // 181  
  1,   // 182  
  1,   // 183  
  1,   // 184  
  1,   // 185  
  1,   // 186  
  1,   // 187  
  1,   // 188  
  1,   // 189  
  1,   // 190  
  1,   // 191  
  1,   // 192  
  1,   // 193  
  1,   // 194  
  1,   // 195  
  1,   // 196  
  1,   // 197  
  1,   // 198  
  1,   // 199  
  1,   // 200  
  1,   // 201  
  1,   // 202  
  1,   // 203  
  1,   // 204  
  1,   // 205  
  1,   // 206  
  1,   // 207  
  1,   // 208  
  1,   // 209  
  1,   // 210  
  1,   // 211  
  1,   // 212  
  1,   // 213  
  1,   // 214  
  1,   // 215  
  1,   // 216  
  1,   // 217  
  1,   // 218  
  1,   // 219  
  1,   // 220  
  1,   // 221  
  1,   // 222  
  1,   // 223  
  1,   // 224  
  1,   // 225  
  1,   // 226  
  1,   // 227  
  1,   // 228  
  1,   // 229  
  1,   // 230  
  1,   // 231  
  1,   // 232  
  1,   // 233  
  1,   // 234  
  1,   // 235  
  1,   // 236  
  1,   // 237  
  1,   // 238  
  1,   // 239  
  1,   // 240  
  1,   // 241  
  1,   // 242  
  1,   // 243  
  1,   // 244  
  1,   // 245  
  1,   // 246  
  1,   // 247  
  1,   // 248  
  1,   // 249  
  1,   // 250  
  1,   // 251  
  1,   // 252  
  1,   // 253  
  1,   // 254  
  1    // 255  
};

// True if a character is a single letter element.
// C N O F S I
// We have a mapping from the letter to the kCarbon... symbols
static int is_single_letter_element[256] = {
  -1,   // 0  
  -1,   // 1  
  -1,   // 2  
  -1,   // 3  
  -1,   // 4  
  -1,   // 5  
  -1,   // 6  
  -1,   // 7  
  -1,   // 8  
  -1,   // 9  
  -1,   // 10  
  -1,   // 11  
  -1,   // 12  
  -1,   // 13  
  -1,   // 14  
  -1,   // 15  
  -1,   // 16  
  -1,   // 17  
  -1,   // 18  
  -1,   // 19  
  -1,   // 20  
  -1,   // 21  
  -1,   // 22  
  -1,   // 23  
  -1,   // 24  
  -1,   // 25  
  -1,   // 26  
  -1,   // 27  
  -1,   // 28  
  -1,   // 29  
  -1,   // 30  
  -1,   // 31  
  -1,   // 32  
  -1,   // 33 !
  -1,   // 34 "
  -1,   // 35 #
  -1,   // 36 $
  -1,   // 37 %
  -1,   // 38 &
  -1,   // 39 '
  -1,   // 40 (
  -1,   // 41 )
   9,   // 42 *  recognised as unspecified atom
  -1,   // 43 +
  -1,   // 44 ,
  -1,   // 45 -
  -1,   // 46 .
  -1,   // 47 /
  -1,   // 48 0
  -1,   // 49 1
  -1,   // 50 2
  -1,   // 51 3
  -1,   // 52 4
  -1,   // 53 5
  -1,   // 54 6
  -1,   // 55 7
  -1,   // 56 8
  -1,   // 57 9
  -1,   // 58 :
  -1,   // 59 ;
  -1,   // 60 <
  -1,   // 61 =
  -1,   // 62 >
  -1,   // 63 ?
  -1,   // 64 @
  -1,   // 65 A
   9,   // 66 B
   0,   // 67 C
  -1,   // 68 D
  -1,   // 69 E
   3,   // 70 F
  -1,   // 71 G
  -1,   // 72 H
   8,   // 73 I
  -1,   // 74 J
  -1,   // 75 K
  -1,   // 76 L
  -1,   // 77 M
   1,   // 78 N
   2,   // 79 O
   4,   // 80 P
  -1,   // 81 Q
  -1,   // 82 R
   5,   // 83 S
  -1,   // 84 T
  -1,   // 85 U
  -1,   // 86 V
  -1,   // 87 W
  -1,   // 88 X
  -1,   // 89 Y
  -1,   // 90 Z
  -1,   // 91 [
  -1,   // 92 \
  -1,   // 93 ]
  -1,   // 94 ^
  -1,   // 95 _
  -1,   // 96 `
  -1,   // 97 a
  -1,   // 98 b
#ifdef NOT_SURE_WHAT_IS_GOING_ON
  same thing here
#endif
  -1,   // 99 c
   0,   // 100 d
  -1,   // 101 e
  -1,   // 102 f
  -1,   // 103 g
  -1,   // 104 h
  -1,   // 105 i
  -1,   // 106 j
  -1,   // 107 k
  -1,   // 108 l
  -1,   // 109 m
  -1,   // 110 n
   1,   // 111 o
   2,   // 112 p
   4,   // 113 q
  -1,   // 114 r
  -1,   // 115 s
   5,   // 116 t
  -1,   // 117 u
  -1,   // 118 v
  -1,   // 119 w
  -1,   // 120 x
  -1,   // 121 y
  -1,   // 122 z
  -1,   // 123 {
  -1,   // 124 |
  -1,   // 125 }
  -1,   // 126 ~
  -1,   // 127  
  -1,   // 128  
  -1,   // 129  
  -1,   // 130  
  -1,   // 131  
  -1,   // 132  
  -1,   // 133  
  -1,   // 134  
  -1,   // 135  
  -1,   // 136  
  -1,   // 137  
  -1,   // 138  
  -1,   // 139  
  -1,   // 140  
  -1,   // 141  
  -1,   // 142  
  -1,   // 143  
  -1,   // 144  
  -1,   // 145  
  -1,   // 146  
  -1,   // 147  
  -1,   // 148  
  -1,   // 149  
  -1,   // 150  
  -1,   // 151  
  -1,   // 152  
  -1,   // 153  
  -1,   // 154  
  -1,   // 155  
  -1,   // 156  
  -1,   // 157  
  -1,   // 158  
  -1,   // 159  
  -1,   // 160  
  -1,   // 161  
  -1,   // 162  
  -1,   // 163  
  -1,   // 164  
  -1,   // 165  
  -1,   // 166  
  -1,   // 167  
  -1,   // 168  
  -1,   // 169  
  -1,   // 170  
  -1,   // 171  
  -1,   // 172  
  -1,   // 173  
  -1,   // 174  
  -1,   // 175  
  -1,   // 176  
  -1,   // 177  
  -1,   // 178  
  -1,   // 179  
  -1,   // 180  
  -1,   // 181  
  -1,   // 182  
  -1,   // 183  
  -1,   // 184  
  -1,   // 185  
  -1,   // 186  
  -1,   // 187  
  -1,   // 188  
  -1,   // 189  
  -1,   // 190  
  -1,   // 191  
  -1,   // 192  
  -1,   // 193  
  -1,   // 194  
  -1,   // 195  
  -1,   // 196  
  -1,   // 197  
  -1,   // 198  
  -1,   // 199  
  -1,   // 200  
  -1,   // 201  
  -1,   // 202  
  -1,   // 203  
  -1,   // 204  
  -1,   // 205  
  -1,   // 206  
  -1,   // 207  
  -1,   // 208  
  -1,   // 209  
  -1,   // 210  
  -1,   // 211  
  -1,   // 212  
  -1,   // 213  
  -1,   // 214  
  -1,   // 215  
  -1,   // 216  
  -1,   // 217  
  -1,   // 218  
  -1,   // 219  
  -1,   // 220  
  -1,   // 221  
  -1,   // 222  
  -1,   // 223  
  -1,   // 224  
  -1,   // 225  
  -1,   // 226  
  -1,   // 227  
  -1,   // 228  
  -1,   // 229  
  -1,   // 230  
  -1,   // 231  
  -1,   // 232  
  -1,   // 233  
  -1,   // 234  
  -1,   // 235  
  -1,   // 236  
  -1,   // 237  
  -1,   // 238  
  -1,   // 239  
  -1,   // 240  
  -1,   // 241  
  -1,   // 242  
  -1,   // 243  
  -1,   // 244  
  -1,   // 245  
  -1,   // 246  
  -1,   // 247  
  -1,   // 248  
  -1,   // 249  
  -1,   // 250  
  -1,   // 251  
  -1,   // 252  
  -1,   // 253  
  -1,   // 254  
  -1   // 255
};

template <typename T>
MolecularFormula<T>::MolecularFormula() :_digits(100) {
  _hash_valid = 0;
  _enable_hash_computation = 1;
}

template <typename T>
int
MolecularFormula<T>::Build(const const_IWSubstring& smiles) {
  const int nchars = smiles.length();

  std::fill_n(_count, kOther + 1, 0);

#ifdef NOT_SURE_WHAT_IS_GOING_ON
  // used to debug this problem.
  if (smiles == "Cr") {
    for (int i = 0; i < 256; ++i) {
      char c = static_cast<char>(i);
      std::cerr << i << " '" << c << "' ignore " << ignore[i] << '\n';
    }
  }
#endif

  int rc = 0;
  for (int i = 0; i < nchars; ++i) {
    unsigned char c = smiles[i];
    if (c == ' ' || c == '\t') {
      return rc;
    }

    // std::cerr << i << " '" << c << " ignore " << ignore[c] << '\n';
    if (ignore[c]) {
      continue;
    }

    // std::cerr << "MolecularFormula::Build:i " << i << " '" << c << "' " << is_single_letter_element[c] << '\n';
    char next_char = '\0';
    if (i != nchars - 1) {
      next_char = smiles[i + 1];
    }

    if (c == 'C' && next_char == 'l') {
      ++_count[kChlorine];
      ++i;
    } else if (c == 'B') {
      if (next_char == 'r') {
        ++_count[kBromine];
        ++i;
      } else {
        ++_count[kBoron];
      }
    } else if (is_single_letter_element[c] >= 0) {
     ++_count[is_single_letter_element[c]];
    } else if (c == kOpenSquareBracket) {
      if (! IdentifyElementSQB(smiles, i)) {
        return 0;
      }
    } else {  // No other elements are recognised outside []
      return 0;
    }

    ++rc;
  }

  return rc;
}

// It is a lot harder to identify the atom inside square brackets.
// [Cr] must not be interpreted as carbon.
// [Tc] must not be interpreted as aromatic carbom
// Don't forget Hydrogen [123CrH2@], [234Tc@@H3]
// [cH], [CH2], [C@H2], [12N+H3]
// Strategy is to skip ove any leading digits, and interpret
// the first letter(s) we see. then skip over everything else
// to the closing square braket.
template <typename T>
int
MolecularFormula<T>::IdentifyElementSQB(const const_IWSubstring& smiles,
                int& i) {
  assert(smiles[i] == kOpenSquareBracket);

  ++i;
  const int nchars = smiles.length();

  // skip over leading digits (isotopes)
  for ( ; i < nchars; ++i) {
    unsigned char c = smiles[i];
    if (c > '9' || c < '0') {
      break;
    }
  }

  unsigned char c = smiles[i];

  if (c == ' ' || c == '\t') {
    std::cerr << "MolecularFormula::IdentifyElementSQB no closing square bracket\n";
    return 0;
  }

  if (c == kCloseSquareBracket) {
    std::cerr << "MolecularFormula::IdentifyElementSQB:no element\n";
    return 0;
  }

  if (i == nchars - 1) {
    std::cerr << "MolecularFormula::IdentifyElementSQB:no closing square bracket\n";
    return 0;
  }

  char next_char = smiles[i + 1];

  if (c >= static_cast<unsigned char>('A') &&
      c <= static_cast<unsigned char>('Z'))  {
    if (next_char >= static_cast<unsigned char>('a') &&
        next_char <= static_cast<unsigned char>('z')) {
      if (c == static_cast<unsigned char>('C') &&
          next_char == static_cast<unsigned char>('l')) {
        ++_count[kChlorine];
        i += 2;
      } else if (c == static_cast<unsigned char>('B')) {
        if (next_char == static_cast<unsigned char>('r')) {
          ++_count[kBromine];
        }
        i += 2;
      } else {
        ++_count[kOther];
        ++i;
      }
    } else if (c == 'B') {
      ++_count[kBoron];
      ++i;
    } else if (is_single_letter_element[c] >= 0) {
      ++_count[is_single_letter_element[c]];
      ++i;
    } else {
      ++_count[kOther];
      ++i;
    }
  } else if (c >= static_cast<unsigned char>('a') &&
             c <= static_cast<unsigned char>('z')) {
    if (is_single_letter_element[c] >= 0) {
      ++_count[is_single_letter_element[c]];
      ++i;
    } else {
      return 0;   // There are no other aromatic elements.
    }
  } else {
    std::cerr << "MolecularFormula::IdentifyElementSQB:not an atomic symbol '" << c << '\n';
    return 0;
  }

  if (i == nchars) {
    return 0;
  }

  for ( ; i < nchars; ++i) {
    unsigned char c = smiles[i];
    if (c == kCloseSquareBracket) {
      return 1;
    }

    if (c == ' ' || c == '\t') {
      std::cerr << "MolecularFormula::IdentifyElementSQB no closing square bracket\n";
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
MolecularFormula<T>::MakeFormula(IWString& s) const {
  s.resize_keep_storage(0);

  return AppendFormula(s);
}

template <typename T>
int
MolecularFormula<T>::AppendFormula(IWString& s) const {
  static IWString asym[] = {
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "Cl",
    "Br",
    "I",
    "B",
    "*"
  };

  int rc = 0;
  for (int i = 0; i < kNTypes; ++i) {
    if (_count[i]) {
      s << asym[i] << _count[i];
      rc += _count[i];
    }
  }

  return rc;
}

template <typename T>
int
MolecularFormula<T>::Build(const Molecule& m) {
  std::fill_n(_count, kOther + 1, 0);

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);
    switch (z) {
      case 6:
        ++_count[kCarbon];
        break;
      case 7:
        ++_count[kNitrogen];
        break;
      case 8:
        ++_count[kOxygen];
        break;
      case 9:
        ++_count[kFluorine];
        break;
      case 15:
        ++_count[kPhosphorus];
        break;
      case 16:
        ++_count[kSulphur];
        break;
      case 17:
        ++_count[kChlorine];
        break;
      case 35:
        ++_count[kBromine];
        break;
      case 53:
        ++_count[kIodine];
        break;
      case 5:
        ++_count[kBoron];
        break;
      default:
        ++_count[kOther];
    }
  }

  return 1;
}


template class MolecularFormula<uint32_t>;

}  // namespace molecular_formula
