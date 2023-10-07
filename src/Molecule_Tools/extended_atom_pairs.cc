/*
  Implementation of Merck atom pairs
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

const char* prog_name = nullptr;

namespace extended_atom_pairs {

using std::cerr;

typedef unsigned int atype_t;
typedef unsigned int aptype_t;

// If writing fingerprints.
IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");

// Tag stem for fingerprints.
IWString tag("NCMAP");

// If writing fixed width fingerprints.
int fixed_width_size = 0;

// Write molecules with pairs labelled.
Molecule_Output_Object stream_for_labelled_molecules;

int verbose = 0;

// Read from a tdt file.
int function_as_filter = 0;

int molecules_read = 0;

// The min and max atom separation at which fingerprints are generated.
int min_distance = 0;
int max_distance = std::numeric_limits<int>::max();

// We can optionally fingerprint all distances that happen to be beyond
// max_distance, but assigning them all an arbitrary value.
int include_distances_beyond_max = 0;

// Sometimes it can be useful to just record the presence of pairs
// within the radius.
int all_distances_identical = 0;

/*
  With the -O option, we can find a subset of pairs. We assign a sequential number to
  each of the pairs.
  The value in the unordered_map will be the column number
*/

std::unordered_map<atype_t, unsigned int> pairs_to_find;

/*
  If we are creating a descriptor file, we need a cross reference between pair number and
  column number
*/

int write_descriptors = 0;

IWString descriptor_prefix("MAP");

/*
  Jul 2005. Need a file with
  id
  type1-dist-type2
  type1-dist-type2
  |
*/

IWString_and_File_Descriptor stream_with_all_pairs;

// Collect statistics on the number of bits set in each molecule.
Accumulator_Int<int> nbits;

// Used for numeric output.
IWDigits iwdigits;

/*
  We set this to 1 because the bond distance between atoms in different
  fragments is undefined.
  Feb 2010, relax this...
*/

int reduce_to_largest_fragment = 0;

int multi_fragment_molecules_encountered = 0;

Element_Transformations etrans;

Chemical_Standardisation chemical_standardisation;

Atom_Typing_Specification atom_typing_specification;

/*
  by default, we set one bit count per feature found. But when dealing
  with fuzziness, we will need to up the default
*/

int default_number_bits_to_set = 1;

/*
  Sometimes people want to figure out what certain bits are
*/

int explain_bits = 0;

int report_all_explain_bits = 1;

extending_resizable_array<int> separation_of_bits_explained;

std::unordered_map<int, IWString> explanation;

/*
  for each bit for which we are looking for an explanation, how many
  instances are found?
*/

std::unordered_map<unsigned int, int> bit_explained;

/*
  We can get some more specificity into the fingerprint by using the
  bonding information for connected atoms

  We can use any combination of
    1. Bond type
    2. Ring Membership
*/

#define AP_BOND_PROPERTY_RING_MEMBERSHIP 1
#define AP_BOND_PROPERTY_BOND_TYPE 2

int bonded_atom_type = 0;

resizable_array_p<Substructure_Query> only_process_queries;

int extra_bit_for_atoms_matched_by_queries = 0;

#define EXTRA_BIT_FOR_ATOMS_MATCHED_BY_QUERIES 411

int molecules_not_hit_by_any_queries = 0;

int molecules_with_only_one_matched_atom = 0;

int included_atoms_needed = 2;

/*
  We want to make fuzzy fingerprints.
  First we have a minimum separation for adding fuzziness
*/

int min_separation_for_fuzziness = std::numeric_limits<int>::max();

/*
  Sep 2016. What if the atom type was just "same" or "different"
*/

int same_or_different = 0;

int apply_isotopic_labels = 1;

static int flush_after_each_molecule = 0;

/*
  Jan 2017. Play with the idea of generating extra pairs to equalise times an atom is used
*/

/*
  A set of fuzziness is described by a set of integers, the largest
  gets applied to the actual distance.  So, 1,8,2 would mean the
  actual separation gets a value of 8, one less a value of 1, and one
  bond longer, a value of 2.
*/

class Fuzziness_Profile {
 private:
  int _minsep;
  int _central_separation;
  int _maxsep;

  int* _value;

 public:
  Fuzziness_Profile();
  ~Fuzziness_Profile();

  int debug_print(std::ostream& os) const;

  int build(const const_IWSubstring&);

  int
  minsep() const {
    return _minsep;
  }

  int
  maxsep() const {
    return _maxsep;
  }

  int
  central_separation() const {
    return _central_separation;
  }

  int
  value_at_separation(int s) const {
    return _value[s];
  }

  int set_bits(const int t1, const int t2, Sparse_Fingerprint_Creator& sfp) const;
};

Fuzziness_Profile** fuzziness_profile = nullptr;

Fuzziness_Profile::Fuzziness_Profile() {
  _minsep = 0;
  _central_separation = 0;
  _maxsep = 0;

  _value = nullptr;

  return;
}

Fuzziness_Profile::~Fuzziness_Profile() {
  if (nullptr != _value) {
    delete[] _value;
  }

  return;
}

/*
  Parse a text description of fuzziness.
  Must look something like

  o:a,b,c

  where `o` is the bond separation, and a,b,c correspond to the
  values associated with o-1,o,o+1 bond distances.
  TODO: swtich to RE2.
*/

int
Fuzziness_Profile::build(const const_IWSubstring& s) {
  const_IWSubstring bsep, fvalues;

  if (!s.split(bsep, ':', fvalues) || 0 == bsep.length() || bsep.length() > 2) {
    cerr << "Fuzziness_Profile::build:invalid separation specification '" << s << "'\n";
    return 0;
  }

  if (!bsep.numeric_value(_central_separation) || _central_separation < 1) {
    cerr << "Fuzziness_Profile::build:invalid separation '" << s << "'\n";
    return 0;
  }

  int n = fvalues.nwords(',');

  if (n < 2) {
    cerr << "Fuzziness_Profile::build:not enough values '" << fvalues << "'\n";
    return 0;
  }

  resizable_array<int> values;

  int i = 0;
  const_IWSubstring token;
  while (fvalues.nextword(token, i, ',')) {
    int v;
    if (!token.numeric_value(v) || v < 0) {
      cerr << "Fuzziness_Profile::build:invalid value '" << fvalues << "'\n";
      return 0;
    }
    values.add(v);
  }

  // Work out the highest number in there

  iwmaxid<int, int> maxval(values[0], 0);

  for (i = 1; i < n; i++) {
    maxval.extra(values[i]);
  }

  int maxwhere = maxval.which_is_max();

  _minsep = _central_separation - maxwhere - 1;
  _maxsep = _central_separation + values.number_elements() - maxwhere - 1;

  int number_to_allocate =
      _central_separation + 1 + values.number_elements();  // maybe waste a little space

  _value = new_int(number_to_allocate);

  for (int i = 0; i < n; i++) {
    int ndx = _central_separation - maxwhere + i - 1;
    _value[ndx] = values[i];
  }

  return 1;
}

int
Fuzziness_Profile::debug_print(std::ostream& os) const {
  os << "Fuzziness_Profile:between " << _minsep << " and " << (_maxsep - 1) << ", centre "
     << _central_separation << '\n';
  for (int i = _minsep; i < _maxsep; i++) {
    os << " at " << i << " bonds, increment " << _value[i] << '\n';
  }

  return 1;
}

/*
  At each separation, we have a fuzziness profile
*/

void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -c <number>    shortest bond distance to process\n";
  cerr << "  -C <number>    longest  bond distance to process\n";
  cerr << "  -D             use a fixed distance for all pairs beyond max distance\n";
  cerr << "  -P <type>      atom type specification, enter '-P help' for info\n";
  cerr << "  -f             work as a filter\n";
  cerr << "  -J <tag>       tag for fingerprints\n";
  cerr << "  -m <nbits>     if producing fixed width fingerprints (-J starts FP), nbits\n";
  cerr << "  -O <number>    only produce pair <number>\n";
  cerr << "  -O F:<file>    only produce pairs in <file>, '$ATYPE <type>' sets type\n";
  cerr << "  -d             produce a descriptor file (requires the -O option)\n";
//cerr << "  -x             write explanations for the bits in the -O file\n";
  cerr << "  -X all         write all explanations for the bits in the -O file\n";
  cerr << "  -X amap        use atom map labels instead of isotopes in the explanation output\n";
  cerr << "  -B <fname>     write raw pair values: type1-dist-type2\n";
  cerr << "  -F s:n,n,n     fuzziness specification, s = separation\n";
  cerr << "  -s <smarts>    only process atoms matching <smarts>\n";
  cerr << "  -q <query>     only process atoms matching <query>\n";
  cerr << "  -k             pairs can have just one matched atom (default is both ends must be matched)\n";
  cerr << "  -y             when using queries, set an extra bit for the number of atoms matched\n";
  cerr << "  -2             use version 2 augmented atom algorithm\n";
  cerr << "  -u             atom types are just 'same' or 'different'\n";
  cerr << "  -e ...         for bonded atoms, use 'type' and/or 'ring' for separation\n";
  cerr << "  -e ...         special processing for bonded atoms, enter '-e help' for infp\n";
  cerr << "  -G ...         obscure options, enter '-G help' for info\n";
  cerr << "  -b <n>         number of bit replicates to produce\n";
  cerr << "  -L <fname>     write atom type labelled atoms to <fname>\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -g ...         chemical standardisation options\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

// Apply either isotopic label `s`, or atom map number
// to atoms `a1` and `a2` in `m`.
void
label_two_atoms(Molecule& m, const atom_number_t a1, const atom_number_t a2,
                const int s) {
  if (apply_isotopic_labels) {
    m.set_isotope(a1, s);
    m.set_isotope(a2, s);
  } else {
    m.set_atom_map_number(a1, s);
    m.set_atom_map_number(a2, s);
  }
}

// If we are looking for bits of type `map`, write that information
// to `output`.
// Args:
//   m: Molecule
//   map: the value for this pair.
//   d: distance associated with this map.
//   a1,a2: the atom numbers.
//   t1,t2: the atom types of a1 and a2
//   ouptut: destination
int
do_explain_bits(const std::unordered_map<atype_t, unsigned int>& pairs_to_find,
                Molecule& m, const atype_t map, const int d, const atom_number_t a1,
                const atom_number_t a2, const atype_t t1, const atype_t t2,
                IWString_and_File_Descriptor& output) {
  if (pairs_to_find.find(map) == pairs_to_find.end()) {
    return 1;
  }

  constexpr char kSep = ' ';

  IWString tmp;
  tmp << map << kSep << d << kSep;
  if (t1 < t2) {
    tmp << t1 << kSep << t2;
  } else {
    tmp << t2 << kSep << t1;
  }

  const auto f = explanation.find(map);

  if (f == explanation.end()) {
    label_two_atoms(m, a1, a2, 1);
    output << m.smiles() << kSep << m.name() << kSep << tmp << kSep
           << m.smarts_equivalent_for_atom(a1) << kSep << m.smarts_equivalent_for_atom(a2)
           << '\n';
    label_two_atoms(m, a1, a2, 0);

    if (!report_all_explain_bits) {
      explanation[map] = tmp;
    }
    //  explanation.insert(std::make_pair<int,IWString>(map, tmp));

    separation_of_bits_explained[d]++;
    bit_explained[map]++;

    return 1;
  }

  if (tmp == (*f).second) {  // same as before
    return 1;
  }

  cerr << "Ambiguous bit, previous '" << (*f).second << ", now '" << tmp << "'\n";

  output << m.smiles() << kSep << m.name() << kSep << map << kSep << d << " *\n";

  return 1;
}

// Write an isotopically labelled version of `m` to `output`.
template <typename T>
int
write_labelled_molecule(Molecule& m, const T* atype, Molecule_Output_Object& output) {
  // If `m` has existing isotopic labels, need to save and restore them.
  std::unique_ptr<isotope_t[]> isosave;

  if (m.number_isotopic_atoms() > 0) {
    isosave = m.GetIsotopes();
  }

  m.set_isotopes(atype);

  output.write(m);

  m.transform_to_non_isotopic_form();

  if (isosave) {
    m.set_isotopes(isosave.get());
  }

  return 1;
}

int
write_empty_fingerprint(IWString_and_File_Descriptor& output) {
  output << tag << ">\n";

  return 1;
}

int
write_01_fingerprint(const IW_Bits_Base& bb, const IWString& tag,
                     IWString_and_File_Descriptor& output) {
  IWString tmp;
  bb.daylight_ascii_representation_including_nset_info(tmp);

  output << tag << tmp << ">\n";

  return output.length();
}

// Write `fp` as fixed form.
int
write_01_fingerprint(const Sparse_Fingerprint_Creator& fp, const IWString& tag,
                     IWString_and_File_Descriptor& output) {
  IW_Bits_Base bb(fixed_width_size);
  for (auto [bit, count] : fp.bits_found()) {
    bb.set(bit % fixed_width_size);
  }

  return write_01_fingerprint(bb, tag, output);
}

int
write_01_fingerprint(const Sparse_Fingerprint_Creator& fp, const IWString& tag,
                     const std::unordered_map<atype_t, unsigned int>& pairs_to_find,
                     IWString_and_File_Descriptor& output) {
  int nbits = pairs_to_find.size();
  if (0 == nbits) {
    cerr << "write_01_fingerprint:no pairs to find, impossible\n";
    output << tag << ">\n";
    return 1;
  }

  const Sparse_Fingerprint_Creator::FPHash& bits_found = fp.bits_found();

  if (0 == fixed_width_size) {
    fixed_width_size = nbits;
    if (0 != fixed_width_size % 8) {
      fixed_width_size = (nbits / 8 + 1) * 8;
      cerr << "Output fingerprints sized to " << fixed_width_size << " bits\n";
    }
  }

  IW_Bits_Base bb(fixed_width_size);

  for (Sparse_Fingerprint_Creator::FPHash::const_iterator f = bits_found.begin();
       f != bits_found.end(); ++f) {
    int zbit = (*f).first;

    std::unordered_map<atype_t, unsigned int>::const_iterator zcol =
        pairs_to_find.find(zbit);

    if (zcol == pairs_to_find.end())  // should not happen
    {
      cerr << "Huh, bit " << zbit << " missing from column hash\n";
      continue;
    }

    int b = (*zcol).second;

    assert(b >= 0 && b < nbits);

    bb.set(b);
  }

  return write_01_fingerprint(bb, tag, output);
}

int
do_write_descriptors(const Sparse_Fingerprint_Creator& fp,
                     const std::unordered_map<atype_t, unsigned int>& pairs_to_find,
                     IWString_and_File_Descriptor& output) {
  const int number_descriptors = pairs_to_find.size();
  std::unique_ptr<IWString[]> tokens = std::make_unique<IWString[]>(number_descriptors);

  const Sparse_Fingerprint_Creator::FPHash& bits_found = fp.bits_found();

  for (Sparse_Fingerprint_Creator::FPHash::const_iterator f = bits_found.begin();
       f != bits_found.end(); ++f) {
    int zbit = (*f).first;
    int zcnt = (*f).second;

    std::unordered_map<atype_t, unsigned int>::const_iterator zcol =
        pairs_to_find.find(zbit);

    if (zcol == pairs_to_find.end())  // should not happen
    {
      cerr << "Huh, bit " << zbit << " missing from column hash\n";
      continue;
    }

    int column = (*zcol).second;

    iwdigits.append_number(tokens[column], zcnt);
  }

  for (int i = 0; i < number_descriptors; i++) {
    if (tokens[i].length()) {
      output << tokens[i];
    } else {
      output << " 0";
    }
  }

  output << '\n';

  return 1;
}

aptype_t
form_atom_pair(atype_t atype1, const int distance, atype_t atype2) {
  if (atype1 > atype2) {
    std::swap(atype1, atype2);
  }

  if (stream_with_all_pairs.active()) {
    IWString& h = stream_with_all_pairs;
    iwdigits.append_number(h, atype1);
    iwdigits.append_number(h, distance);
    iwdigits.append_number(h, atype2);
    stream_with_all_pairs << '\n';
  }

  return 975535 * atype1 + distance * 7927 + atype2;
}

int
bond_constant(const Bond* b) {
  if (b->is_aromatic()) {
    return 4;
  }

  if (b->is_single_bond()) {
    return 1;
  }

  if (b->is_double_bond()) {
    return 2;
  }

  if (b->is_triple_bond()) {
    return 3;
  }

  return 5;
}

aptype_t
form_atom_pair(const Bond* b, const atype_t* atype) {
  atype_t atype1 = atype[b->a1()];
  atype_t atype2 = atype[b->a2()];

  if (atype1 > atype2) {
    std::swap(atype1, atype2);
  }

  auto bc = bond_constant(b);

  if ((bonded_atom_type | AP_BOND_PROPERTY_RING_MEMBERSHIP) && b->nrings()) {
    bc = bc * 11;
  }

  return 41202 * atype1 + bc * 2121 + atype2;
}

int
Fuzziness_Profile::set_bits(const int t1, const int t2,
                            Sparse_Fingerprint_Creator& sfp) const {
  for (int i = _minsep; i < _maxsep; i++) {
    const aptype_t map = form_atom_pair(t1, i, t2);
    assert(_value[i] > 0);
    sfp.hit_bit(map, _value[i]);
  }

  return 1;
}

/*
  Scaled back version that does not support all available features,
  but will allow multi-fragment molecules to be fingerprinted
*/

int
ExtendedAtomPairs_multi_fragment(Molecule& m, const IWString& tag, atype_t* atype,
                                 IWString_and_File_Descriptor& output) {
  multi_fragment_molecules_encountered++;

  Sparse_Fingerprint_Creator fp;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int fragi = m.fragment_membership(i);

    const auto atypei = atype[i];

    for (int j = i + 1; j < matoms; j++) {
      if (m.fragment_membership(j) != fragi) {
        continue;
      }

      int d = m.bonds_between(i, j);

      if (d < min_distance) {
        continue;
      }

      if (d <= max_distance) {
        if (all_distances_identical) {
          d = max_distance;
        }
      } else if (include_distances_beyond_max) {
        d = include_distances_beyond_max;
      } else {
        continue;
      }

      if (d > 1 || 0 == bonded_atom_type) {
        const aptype_t map = form_atom_pair(atypei, d, atype[j]);

        fp.hit_bit(map, default_number_bits_to_set);
      } else  // bonded atoms, and special processing in that case
      {
        const aptype_t map = form_atom_pair(m.bond_between_atoms(i, j), atype);

        fp.hit_bit(map, default_number_bits_to_set);
      }
    }
  }

  if (verbose) {
    nbits.extra(fp.nbits());
  }

  if (tag.starts_with("NC")) {
    IWString tmp;
    fp.daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << "\n";
    return 1;
  }

  return write_01_fingerprint(fp, tag, pairs_to_find, output);
}

// form a bit for the pair `a1` `a2` which are `d` bonds apart in `m`.
aptype_t
bit_number_same_or_different(const Molecule& m, const int d, const atom_number_t a1,
                             const atom_number_t a2, const atype_t* atype) {
  if (1 == d && 0 != bonded_atom_type) {
    const Bond* b = m.bond_between_atoms(a1, a2);

    const auto bc = bond_constant(b);

    if (atype[a1] == atype[a2]) {
      return form_atom_pair(734 + bc, 1, 734 + bc);
    } else {
      return form_atom_pair(58 + bc, 1, 158 + bc);
    }
  } else {
    if (atype[a1] == atype[a2]) {
      return form_atom_pair(101, d, 101);
    } else {
      return form_atom_pair(2, d, 50);
    }
  }

  cerr << "bit_number_same_or_different:gack!\n";
  return 0;  // should never happen
}

/*
  Fast version.
  Actually I did some timing and it really does not make much
  of a difference. Next time I need to work with
  this code, just consolidate things


  Scan for atom pairs in `m`.
*/
int
ExtendedAtomPairs_default(Molecule& m, const IWString& tag, atype_t* atype,
                          IWString_and_File_Descriptor& output) {
  assert(0 == explain_bits);

  Sparse_Fingerprint_Creator fp;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int atypei = atype[i];
    for (int j = i + 1; j < matoms; j++) {
      int d = m.bonds_between(i, j);
      //    cerr << d << " bonds between " << i << " and " << j << '\n';

      if (d < min_distance) {
        continue;
      }

      if (d <= max_distance) {
        if (all_distances_identical) {
          d = max_distance;
        }
      } else if (include_distances_beyond_max) {
        d = include_distances_beyond_max;
      } else {
        continue;
      }

      aptype_t atpair;

      if (same_or_different) {
        atpair = bit_number_same_or_different(m, d, i, j, atype);
      } else if (d > 1 || 0 == bonded_atom_type) {
        atpair = form_atom_pair(atypei, d, atype[j]);
      } else {  // bonded atoms, and special processing in that case
        atpair = form_atom_pair(m.bond_between_atoms(i, j), atype);
      }

      fp.hit_bit(atpair, default_number_bits_to_set);
    }
  }

  if (tag.starts_with("NC")) {
    IWString tmp;
    fp.daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << "\n";
  } else if (pairs_to_find.size() > 0) {
    return write_01_fingerprint(fp, tag, pairs_to_find, output);
  } else {
    return write_01_fingerprint(fp, tag, output);
  }

  return 1;
}

// Given a set of queries, set entries in `include_these_atoms` according
// to which atoms are matched by the queries.
int
identify_atoms_to_process(Molecule& m, int* include_these_atoms,
                          resizable_array_p<Substructure_Query>& only_process_queries) {
  Molecule_to_Match target(&m);

  const int nq = only_process_queries.number_elements();

  int rc = 0;

  for (int i = 0; i < nq; ++i) {
    Substructure_Results sresults;

    const int nhits = only_process_queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    rc++;

    sresults.each_embedding_set_vector(include_these_atoms, 1);
  }

  return rc;
}

// Atoms `a1` and `a2` are a possible atom pair. If they are in distance range,
// make a bit in `fp`.
void
form_atom_pair(Molecule& m, const atype_t* atype, const atom_number_t a1,
               const atom_number_t a2, Sparse_Fingerprint_Creator& fp,
               IWString_and_File_Descriptor& output) {
  int d = m.bonds_between(a1, a2);
  // cerr << d << " bonds between " << a1 << " and " << a2 << '\n';

  if (d < min_distance) {
    return;
  }

  if (d <= max_distance) {
    if (all_distances_identical) {
      d = max_distance;
    }
  } else if (include_distances_beyond_max) {
    d = include_distances_beyond_max;
  } else {
    return;
  }

  if (nullptr != fuzziness_profile && nullptr != fuzziness_profile[d]) {
    fuzziness_profile[d]->set_bits(atype[a1], atype[a2], fp);
  } else if (explain_bits) {
    const auto map = form_atom_pair(atype[a1], d, atype[a2]);
    do_explain_bits(pairs_to_find, m, map, d, a1, a2, atype[a1], atype[a2], output);
  } else if (same_or_different) {
    const auto amap = bit_number_same_or_different(m, d, a1, a2, atype);

    fp.hit_bit(amap, default_number_bits_to_set);
  } else if (d > 1 || 0 == bonded_atom_type) {
    const auto map = form_atom_pair(atype[a1], d, atype[a2]);

    if (pairs_to_find.size() && pairs_to_find.find(map) == pairs_to_find.end()) {
      return;
    }

    fp.hit_bit(map, default_number_bits_to_set);
  } else  // bonded atoms, and special processing in that case
  {
    const auto map = form_atom_pair(m.bond_between_atoms(a1, a2), atype);

    if (pairs_to_find.size() && pairs_to_find.find(map) == pairs_to_find.end()) {
      return;
    }

    fp.hit_bit(map, default_number_bits_to_set);
  }

  output.write_if_buffer_holds_more_than(8192);

  return;
}

int
ExtendedAtomPairs_matched_atoms(Molecule& m, const IWString& tag, atype_t* atype,
                                IWString_and_File_Descriptor& output) {
  assert(only_process_queries.number_elements() > 0);

  const int matoms = m.natoms();

  int* include_these_atoms = new_int(matoms);
  std::unique_ptr<int[]> free_include_these_atoms(include_these_atoms);

  if (!identify_atoms_to_process(m, include_these_atoms, only_process_queries)) {
    molecules_not_hit_by_any_queries++;
    if (verbose) {
      cerr << "Cannot identify any atoms to process in '" << m.name() << "'\n";
    }

    if (function_as_filter) {
      write_empty_fingerprint(output);
    }

    return 0;
  }

  const int atoms_matched =
      std::count(include_these_atoms, include_these_atoms + matoms, 1);

  if (1 == atoms_matched) {
    molecules_with_only_one_matched_atom++;
    if (verbose) {
      cerr << "Only one matched atom in '" << m.name() << "'\n";
    }

    if (2 == included_atoms_needed) {
      write_empty_fingerprint(output);
      return 0;
    }
  }

  if (stream_for_labelled_molecules.active()) {
    write_labelled_molecule(m, include_these_atoms, stream_for_labelled_molecules);
  }

  // cerr << "LIne " << __LINE__ << " D " << write_descriptors << " function_as_filter "
  // << function_as_filter << '\n';
  if (write_descriptors) {
    ;
  } else if (function_as_filter) {
    ;
  } else {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator fp;

  if (extra_bit_for_atoms_matched_by_queries) {
    fp.hit_bit(EXTRA_BIT_FOR_ATOMS_MATCHED_BY_QUERIES, atoms_matched);
  }

  if (2 == included_atoms_needed) {
    for (int i = 0; i < matoms; i++) {
      if (0 == include_these_atoms[i]) {
        continue;
      }

      for (int j = i + 1; j < matoms; j++) {
        if (0 == include_these_atoms[j]) {
          continue;
        }

        form_atom_pair(m, atype, i, j, fp, output);
      }
    }
  } else {
    for (int i = 0; i < matoms; i++) {
      for (int j = i + 1; j < matoms; ++j) {
        if (include_these_atoms[i] ||
            include_these_atoms[j]) {  // interesting to contemplate an XOR?
          form_atom_pair(m, atype, i, j, fp, output);
        }
      }
    }
  }

  if (explain_bits) {
    return 1;
  }

  if (verbose) {
    nbits.extra(fp.nbits());
  }

  if (write_descriptors) {
    return do_write_descriptors(fp, pairs_to_find, output);
  }

  if (stream_with_all_pairs.is_open()) {
    stream_with_all_pairs << "|\n";
    stream_with_all_pairs.write_if_buffer_holds_more_than(32768);
  }

  if (tag.starts_with("NC")) {
    IWString tmp;
    fp.daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << "\n";
  } else if (pairs_to_find.size() > 0) {
    return write_01_fingerprint(fp, tag, pairs_to_find, output);
  } else {
    return write_01_fingerprint(fp, tag, output);
  }

  if (!function_as_filter) {
    output << "|\n";
  }

  return 1;
}

// Fast version for forming atom pairs.
int
ExtendedAtomPairs(Molecule& m, const IWString& tag, const atype_t* atype,
                  IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  Sparse_Fingerprint_Creator fp;

  for (int i = 0; i < matoms; i++) {
    for (int j = i + 1; j < matoms; ++j) {
      form_atom_pair(m, atype, i, j, fp, output);
    }
  }

  if (explain_bits) {
    return 1;
  }

  if (verbose) {
    nbits.extra(fp.nbits());
  }

  if (write_descriptors) {
    return do_write_descriptors(fp, pairs_to_find, output);
  }

  if (stream_with_all_pairs.is_open()) {
    stream_with_all_pairs << "|\n";
    stream_with_all_pairs.write_if_buffer_holds_more_than(32768);
  }

#ifdef DEBUG_BITS_FORMED
  IWString tmp;
  fp.write_in_svml_form(tmp);
  cerr << tmp << '\n';
#endif

  if (tag.starts_with("NC")) {
    IWString tmp;
    fp.daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << "\n";
  } else if (pairs_to_find.size() > 0) {
    return write_01_fingerprint(fp, tag, pairs_to_find, output);
  } else {
    return write_01_fingerprint(fp, tag, output);
  }

  return 1;
}

int
ExtendedAtomPairs(Molecule& m, IWString_and_File_Descriptor& output) {
  // cerr << "Processing " << m.smiles() << '\n';
  if (write_descriptors) {
    append_first_token_of_name(m.name(), output);
  } else if (function_as_filter) {
  } else if (explain_bits) {
  } else if (only_process_queries.number_elements() > 0) {
  } else {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  m.recompute_distance_matrix();

  atype_t* atype = new atype_t[m.natoms()];
  std::unique_ptr<atype_t[]> free_atype(atype);

  if (!atom_typing_specification.assign_atom_types(m, atype)) {
    cerr << "Cannot assign atom types, '" << m.name() << "'\n";
    return 0;
  }
  // cerr << "After assigning atom types " << m.smiles() << '\n';

  if (stream_for_labelled_molecules.active()) {
    write_labelled_molecule(m, atype, stream_for_labelled_molecules);
  }

  if (stream_with_all_pairs.active()) {
    append_first_token_of_name(m.name(), stream_with_all_pairs);
    stream_with_all_pairs << '\n';
  }

  if (bonded_atom_type) {
    //  if (bonded_atom_type | AP_BOND_PROPERTY_RING_MEMBERSHIP)    not necessary given
    //  that aromaticity is done below
    //    m.ring_membership();
    m.compute_aromaticity_if_needed();
  }

  if (reduce_to_largest_fragment && 0 == only_process_queries.number_elements() &&
      0 == explain_bits && 0 == pairs_to_find.size() && 0 == verbose) {
    (void)ExtendedAtomPairs_default(m, tag, atype, output);
  } else if (only_process_queries.number_elements() > 0) {
    (void)ExtendedAtomPairs_matched_atoms(m, tag, atype, output);
  } else if (m.number_fragments() > 1 && 0 == explain_bits) {
    (void)ExtendedAtomPairs_multi_fragment(m, tag, atype, output);
  } else {
    (void)ExtendedAtomPairs(m, tag, atype, output);
  }

  if (write_descriptors) {
    ;
  } else if (function_as_filter) {
    ;
  } else if (explain_bits) {
    ;
  } else if (only_process_queries.number_elements() > 0) {
    ;
  } else {
    output << "|\n";
  }

  return 1;
}

int
parse_atom_type_specification(const const_IWSubstring& buffer) {
  assert(buffer.starts_with("$ATYPE"));

  int i = 0;
  IWString token;

  buffer.nextword(token, i);
  if (!buffer.nextword(token, i)) {
    cerr << "parse_atom_type_specification:directive must be followed by atom type\n";
    return 0;
  }

  return atom_typing_specification.build(token);
}

// Return true if the string representation of any of the tokens in `buffer`
// is a key in `pairs_to_find`.
int
hash_only_want_pair(const const_IWSubstring& buffer,
                    std::unordered_map<atype_t, unsigned int>& pairs_to_find)

{
  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i)) {
    atype_t p;
    if (!token.numeric_value(p)) {
      cerr << "hash_only_want_pair: invalid numeric '" << token << "'\n";
      return 0;
    }

    const auto s = pairs_to_find.size();

    pairs_to_find[p] = s;
    //  pairs_to_find.insert(std::make_pair<int,unsigned int>(p, s));

    //  I initially wrote this for accepting all tokens on the line, but
    //  that has not been needed. Should issue warning.

    return 1;
  }

  return 1;
}

// Populate `pairs_to_find` with the bit numbers read from `input`.
int
process_only_want_pairs_from_file(
    iwstring_data_source& input,
    std::unordered_map<atype_t, unsigned int>& pairs_to_find) {
  input.set_skip_blank_lines(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    buffer.strip_leading_blanks();

    if (buffer.starts_with("$ATYPE")) {
      if (!parse_atom_type_specification(buffer)) {
        cerr << "Invalid atom type specification '" << buffer << "'\n";
        return 0;
      }
    }

    if (!hash_only_want_pair(buffer, pairs_to_find)) {
      cerr << "Fatal error reading pair specification on line " << input.lines_read()
           << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
process_only_want_pairs_from_file(
    const const_IWSubstring& fname,
    std::unordered_map<atype_t, unsigned int>& pairs_to_find) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "process_only_want_pairs_from_file: cannot open input file '" << fname
         << "'\n";
    return 0;
  }

  return process_only_want_pairs_from_file(input, pairs_to_find);
}

int
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (etrans.active()) {
    etrans.process(m);
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return 1;
}

static void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(32768);
  }
}

int
ExtendedAtomPairs(data_source_and_type<Molecule>& input,
                  IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!ExtendedAtomPairs(*m, output)) {
      return 0;
    }

    MaybeFlush(output);
  }

  return 1;
}

// Process input from a fingerprint file.
int
ExtendedAtomPairs(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output) {
  const_IWSubstring smiles(buffer);

  assert(buffer.ends_with('>'));
  smiles.chop();
  smiles.remove_up_to_first('<');

  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    cerr << "Cannot interpret smiles '" << smiles << "'\n";
    return 0;
  }

  return ExtendedAtomPairs(m, output);
}

// Input from a fingerprint file.
int
ExtendedAtomPairs(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!ExtendedAtomPairs(buffer, output)) {
      cerr << "Fatal error '" << buffer << "' line " << input.lines_read() << '\n';
      return 0;
    }

    MaybeFlush(output);
  }

  return 1;
}

int
ExtendedAtomPairs(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ExtendedAtomPairs(input, output);
}

int
ExtendedAtomPairs(const char* fname, FileType input_type,
                  IWString_and_File_Descriptor& output) {
  if (function_as_filter) {
    return ExtendedAtomPairs(fname, output);
  }

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ExtendedAtomPairs(input, output);
}

// Given a set of bits to find in `pairs_to_find`, write a header
// to file descriptor `output_fd`.
int
write_header(const std::unordered_map<atype_t, unsigned int>& pairs_to_find,
             int output_fd) {
  IWString buffer;

  buffer << "name";

  const int number_descriptors = pairs_to_find.size();
  std::unique_ptr<IWString[]> tokens = std::make_unique<IWString[]>(number_descriptors);

  for (auto f : pairs_to_find) {
    IWString& s = tokens[f.second];
    s << descriptor_prefix << f.first;
  }

  for (int i = 0; i < number_descriptors; i++) {
    const IWString& s = tokens[i];

    if (0 == s.length()) {
      cerr << "Bad news, column " << i << " not specified\n";
      return 0;
    }

    buffer << ' ' << s;
  }

  buffer << '\n';

  return buffer.write(output_fd);
}

int
read_fuzziness_specifications_record(const const_IWSubstring& buffer) {
  Fuzziness_Profile* fp = new Fuzziness_Profile;

  if (!fp->build(buffer)) {
    delete fp;
    return 0;
  }

  int c = fp->central_separation();
  if (c > max_distance) {
    cerr << "Fuzziness specification out of range, max distance " << max_distance << '\n';
    return 0;
  }

  if (nullptr != fuzziness_profile[c]) {
    cerr << "Warning, duplicate specification of fuzziness for separation " << c << '\n';
  }

  fuzziness_profile[c] = fp;

  return 1;
}

int
read_fuzziness_specifications_from_file(iwstring_data_source& input) {
  input.set_dos(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (0 == buffer.length()) {
      continue;
    }

    if (!read_fuzziness_specifications_record(buffer)) {
      cerr << "Invalid fuzziness specification '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
read_fuzziness_specifications_from_file(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open fuzziness file specification '" << fname << "'\n";
    return 0;
  }

  return read_fuzziness_specifications_from_file(input);
}

void
DisplayDashGOptins(std::ostream& output) {
  output << " -G samedist     distances within range are identical\n";
  output << " -G help         this message\n";
}

int
ExtendedAtomPairs(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:lg:T:c:C:kJ:fP:O:dL:i:B:F:b:xm:e:2s:q:yuD:X:G:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('A')) {
    set_global_aromaticity_type(Daylight);
  } else if (!process_standard_aromaticity_options(cl, verbose > 1)) {
    cerr << "Cannot process -A option\n";
    usage(11);
  }

  if (!process_elements(cl, verbose > 1, 'E')) {
    cerr << "Cannot initialise elements\n";
    usage(8);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      usage(5);
    }
  }

  if (cl.option_present('2')) {
    set_use_version_2_augmented_atom_algorithm(1);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', min_distance) || min_distance < 1) {
      cerr << "The minimum distance option (-c) must be a whole positive number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will only fingerprint distances " << min_distance << " or shorter\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_distance) || max_distance < min_distance) {
      cerr << "The maximum distance option (-C) must be a whole positive number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will only fingerprint distances " << max_distance << " or longer\n";
    }
  }

  if (cl.option_present('D')) {
    if (!cl.option_present('C')) {
      cerr << "The fixed distance for out of range pairs option (-D) must be used with "
              "the -C option\n";
      usage(1);
    }

    include_distances_beyond_max = 5000;  // just an arbitrary number

    if (verbose) {
      cerr << "Pairs beyond " << max_distance
           << " fingerprinted with constant separation\n";
    }
  }

  if (cl.option_present('k')) {
    included_atoms_needed = 1;

    if (verbose) {
      cerr << "Only needs one matched atom in each pair\n";
    }
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', default_number_bits_to_set) || default_number_bits_to_set < 1) {
      cerr << "The default number of bits to set (-b) must be a +ve whole number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "By default will set one bit count per feature\n";
    }
  }

  if (cl.option_present('d')) {
    if (cl.option_present('f')) {
      cerr << "Cannot write descriptors when working as a filter\n";
      return 8;
    }

    if (!cl.option_present('O')) {
      cerr << "When writing descriptors, must specify the descriptors to produce via the "
              "-O option\n";
      usage(6);
    }

    write_descriptors = 1;

    if (verbose) {
      cerr << "Will write descriptors\n";
    }
  }

  if (cl.option_present('x') && cl.option_present('X')) {
    cerr << "Cannot use both -x and -X options\n";
    return 1;
  }

  if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring s;

    while (cl.value('s', s, i++)) {
      Substructure_Query* q = new Substructure_Query;
      if (!q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 4;
      }

      only_process_queries.add(q);
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, only_process_queries, verbose, 'q')) {
      cerr << "Cannot process queries (-q)\n";
      return 4;
    }
  }

  if (verbose) {
    cerr << "Have " << only_process_queries.size() << " queries for atoms to match\n";
  }

  if (cl.option_present('u')) {
    same_or_different = 1;

    if (verbose) {
      cerr << "Will produce same or different fingerprints\n";
    }
  }

  if (cl.option_present('G')) {
    const_IWSubstring g;
    for (int i = 0; cl.value('g', g, i); ++i) {
      if (g == "samedist") {
        all_distances_identical = 1;
        if (verbose) {
          cerr << "Add distances set to max\n";
        }
      } else if (g == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (g == "help") {
        DisplayDashGOptins(cerr);
      } else {
        cerr << "Unrecognised -G qualifier '" << g << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('e')) {
    const_IWSubstring e;

    for (auto i = 0; cl.value('e', e, i); ++i) {
      if ("type" == e) {
        bonded_atom_type |= AP_BOND_PROPERTY_BOND_TYPE;
      } else if (e.starts_with("ring")) {
        bonded_atom_type |= AP_BOND_PROPERTY_RING_MEMBERSHIP;
      } else if ("help" == e) {
        cerr << " -e type           include atom type in bit formation\n";
        cerr << " -e ring           include ring status in bit formation\n";
        return 1;
      } else {
        cerr << "Unrecognised -e qualifier '" << e << "'\n";
        usage(2);
      }
    }
  }

  if (cl.option_present('y')) {
    extra_bit_for_atoms_matched_by_queries = 1;

    if (verbose) {
      cerr << "When doing molecular subsets, will set an extra bit for the number of "
              "atoms hit\n";
    }
  }

  if (cl.option_present('O')) {
    if (cl.option_present('X')) {
      explain_bits = 1;

      const_IWSubstring x;
      for (int i = 0; cl.value('X', x, i); ++i) {
        if ("all" == x) {
          report_all_explain_bits = std::numeric_limits<int>::max();
        } else if ("amap" == x) {
          apply_isotopic_labels = 0;
        } else if ("1" == x) {
          ;
        } else {
          cerr << "Unrecognised -X directive '" << x << "'\n";
          return 1;
        }
      }

      if (verbose) {
        cerr << "Will only produce explanations for bits\n";
      }
    }

    //  pairs_to_find.resize(2000);

    IWString o;
    for (int i = 0; cl.value('O', o, i); ++i) {
      if (o.starts_with("F:")) {
        o.remove_leading_chars(2);
        if (!process_only_want_pairs_from_file(o, pairs_to_find)) {
          return i + 1;
        }
      } else {
        if (!hash_only_want_pair(o, pairs_to_find)) {
          return i + 1;
        }

        if (verbose) {
          cerr << "Will process pair '" << o << "'\n";
        }
      }
    }

    if (explain_bits) {
      ;
    } else if (write_descriptors) {
      write_header(pairs_to_find, 1);

      if (verbose) {
        for (auto i : pairs_to_find) {
          cerr << "Pair '" << i.first << "' is number " << i.second << '\n';
        }
      }
    }
  }

  if (cl.option_present('F')) {
    if (!cl.option_present('C')) {
      cerr << "In order to use fuzziness profiles (-F) must specify max separation via "
              "the -C option\n";
      usage(3);
    }

    if (!cl.option_present('b')) {
      cerr << "If you are using fuzziness, you must specify default bit set count via "
              "the -b option\n";
      usage(3);
    }

    fuzziness_profile = new Fuzziness_Profile*[max_distance + 1];
    std::unique_ptr<Fuzziness_Profile*[]> free_fuzziness_profile(Fuzziness_Profile);
    set_vector(fuzziness_profile, max_distance + 1,
               static_cast<Fuzziness_Profile*>(nullptr));

    int number_fuzziness_profiles = cl.option_count('F');

    if (1 == number_fuzziness_profiles && cl.string_value('F').starts_with("FILE=")) {
      const_IWSubstring f = cl.string_value('F');
      f.remove_leading_chars(5);
      if (!read_fuzziness_specifications_from_file(f)) {
        cerr << "Cannot read fuzziness specifications from '" << f << "'\n";
        return 4;
      }
    } else {
      const_IWSubstring f;
      for (int i = 0; i < number_fuzziness_profiles; i++) {
        cl.value('F', f, i);
        Fuzziness_Profile* fp = new Fuzziness_Profile;  // leaks if not used.

        if (!fp->build(f)) {
          cerr << "Malformed fuzziness specification '" << f << "'\n";
          return 3;
        }

        if (fp->central_separation() > max_distance) {
          cerr << "Max distance " << max_distance << " so cannot specify fuzziness at "
               << fp->central_separation() << '\n';
          return 3;
        }

        if (nullptr != fuzziness_profile[fp->central_separation()]) {
          cerr << "Warning, duplicate specification for separation "
               << fp->central_separation() << '\n';
        }

        fuzziness_profile[fp->central_separation()] = fp;
      }

      if (verbose) {
        cerr << "Defined " << number_fuzziness_profiles << " fuzziness profiles\n";
        for (int i = 0; i <= max_distance; i++) {
          if (nullptr != fuzziness_profile[i]) {
            fuzziness_profile[i]->debug_print(cerr);
          }
        }
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('f')) {
    function_as_filter = 1;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  int user_specified_tag = 0;

  if (cl.option_present('J')) {
    tag = cl.string_value('J');

    if (verbose) {
      cerr << "Fingerprint(s) produced with tag '" << tag << "'\n";
    }

    user_specified_tag = 1;

    if (!tag.ends_with('<')) {
      tag << '<';
    }
  }

  if (cl.option_present('m')) {
    if (!tag.starts_with("FP")) {
      cerr << "Fixed width fingerprint tag nbits specified (-m) but fingerprint tag not "
              "FP '"
           << tag << "'\n";
      return 2;
    }

    if (!cl.value('m', fixed_width_size) || fixed_width_size < 8 ||
        0 != fixed_width_size % 8) {
      cerr << "The fixed width fingerprint bit count (-m) must be a whole +ve number "
              "divisible by 8\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Fixed width fingerprints produced, " << fixed_width_size << " bits\n";
    }
  }

  if (!tag.starts_with("NC") && 0 == fixed_width_size) {
    cerr << "Fingerprint not sparse '" << tag << "' but no fixed width for fingerprint\n";
    usage(1);
  }

  int nfp = cl.option_count('P');

  if (0 == nfp) {
    atom_typing_specification.set_atom_type(IWATTYPE_COMPLEX);

    if (!user_specified_tag) {
      tag = "NCC0MAP";
      if (std::numeric_limits<int>::max() != max_distance) {
        tag << max_distance;
      }
      tag << '<';
    }

    if (verbose) {
      cerr << "Doing 'C' atom type by default, tag '" << tag << "'\n";
    }
  } else if (nfp > 1) {
    cerr << "Sorry, map no longer supports multiple atom types\n";
    exit(3);
  } else {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot determine atom typying to use '" << p << "'\n";
      return 3;
    }

    if (!user_specified_tag) {
      tag = "NC";
      atom_typing_specification.append_to_tag(tag);
      if (std::numeric_limits<int>::max() != max_distance) {
        tag << max_distance;
      }

      tag << '<';
    }

    if (verbose) {
      cerr << "Fingerprint tag set to '" << tag << "'\n";
    }
  }

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(100);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('L')) {
    const_IWSubstring l = cl.string_value('L');

    stream_for_labelled_molecules.add_output_type(FILE_TYPE_SMI);

    if (stream_for_labelled_molecules.would_overwrite_input_files(cl, l)) {
      cerr << "The -L file cannot overwrite its input\n";
      return 4;
    }

    if (!stream_for_labelled_molecules.new_stem(l)) {
      cerr << "Cannot initialise labelled molecule stream '" << l << "'\n";
      return 6;
    }

    if (verbose) {
      cerr << "Labelled molecules written to '" << l << "'\n";
    }
  }

  if (cl.option_present('B')) {
    const char* b = cl.option_value('B');

    if (!stream_with_all_pairs.open(b)) {
      cerr << "Cannot open stream for Howard data '" << b << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Raw type data written to '" << b << "'\n";
    }

    stream_with_all_pairs.resize(9000);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!ExtendedAtomPairs(cl[i], input_type, output)) {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      rc = i + 1;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";

    if (only_process_queries.number_elements()) {
      cerr << molecules_not_hit_by_any_queries << " molecules did not hit any of the "
           << only_process_queries.number_elements() << " atom selection queries\n";
      cerr << molecules_with_only_one_matched_atom
           << " molecules only matched one atom\n";
    }
  }

  if (explain_bits) {
    int pairs_found = 0;
    for (int i = 0; i < separation_of_bits_explained.number_elements(); i++) {
      if (separation_of_bits_explained[i]) {
        cerr << separation_of_bits_explained[i] << " bits were for atoms " << i
             << " bonds apart\n";
        pairs_found += separation_of_bits_explained[i];
      }
    }

    cerr << "Found " << pairs_found << " examples of instances of "
         << pairs_to_find.size() << " bits to find\n";

    pairs_found = 0;
    for (auto p : pairs_to_find) {
      const auto f = bit_explained.find(p.first);
      if (f == bit_explained.end()) {
        cerr << 0;
      } else {
        cerr << f->second;
        pairs_found++;
      }
      cerr << " instances of bit " << p.first << '\n';
    }

    cerr << "Found examples for " << pairs_found << " of " << pairs_to_find.size()
         << " pairs\n";
  } else if (write_descriptors) {
    ;
  } else if (verbose) {
    cerr << "Fingerprint '" << tag << "' between " << nbits.minval() << " and "
         << nbits.maxval();
    if (nbits.n() > 1) {
      cerr << " ave " << nbits.average();
    }
    cerr << " bits set\n";

    if (multi_fragment_molecules_encountered) {
      cerr << multi_fragment_molecules_encountered
           << " multi fragment molecules processed\n";
    }
  }

  return rc;
}

}  // namespace extended_atom_pairs

int
main(int argc, char** argv) {
  int rc = extended_atom_pairs::ExtendedAtomPairs(argc, argv);

  return rc;
}
