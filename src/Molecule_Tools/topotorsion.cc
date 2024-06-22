/*
  Computes topological torsions
*/

#include <assert.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

#include "absl/container/flat_hash_set.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

#include "topotorsion_support.h"
#include "torsion_hash.h"

using std::cerr;

static int verbose = 0;
static int molecules_read = 0;

static int function_as_filter = 0;

static IWDigits iwdigits;

/*
  If we are writing our results as a fingerprint, we need a tag
*/

static IWString tag;

static IWString topotorsion_prefix("tt_");

const char* prog_name = nullptr;

static Torsion_Hash torsion_hash;

static Sparse_Fingerprint_Creator sparse_fp_creator;

/*
  We can omit torsions which are not hit enough or hit always
*/

static int lower_torsion_count = 0;
static int lower_torsion_percent = 0;
static int upper_torsion_count = 0;
static int upper_torsion_percent = 0;

/*
  There are a number of kinds of output.
    We can write the actual torsions for each molecule.
    We can write torsion numbers with each molecule and then write a cross
    reference file.
    We can accumulate the torsions for each molecule and write everything
    at the end.
*/

static int write_torsion_hash = 0;
static int write_torsion_strings = 1;  // the default behaviour

static int produce_descriptor_file = 0;

/*
  When doing hit followup work, we may only want a small number of
  torsions. We keep Torsion_Hash of the desired torsions in
  that case
*/

static Torsion_Hash torsions_to_find;

/*
  One option is to directly produce a descriptor file.
  This is a two-phase operation.
  As each torsion is encountered, assign it a number.
  For each molecule, accumulate a list of the torsions associated with that
  molecule.
  Then, when we are done, go back to each molecule and fill out the descriptor
  matrix

  When I first thought of this, I thought I'd create a class which would hold
  the torsion number and the count for that molecule, but since this would need
  to a resizable_array_p<>, it would take a lot of storage. Therefore we
  go the dangerous route and have separate arrays
*/

static resizable_array<int>* moltor = nullptr;  // list of torsions in each molecule
static resizable_array<int>* moltorcount =
    nullptr;  // number of occurrences in each molecule
static IWString* stored_name = nullptr;

static int molecules_stored = 0;

/*
  We can generate two kinds of names. SAS requires descriptors to be eight
  characters or less, so the default names are compressed to 8 characters.
*/

static int eight_character_torsion_names = 1;

/*
  Dec 98. Further problems with SAS. It converts everything to uppercase
  We need case insensitive torsion names, so let's use numbers
*/

static int numeric_torsion_names = 0;

// Apr 2024. Allow writing a fixed width descriptor file.
static uint32_t fixed_ncols = 0;

// Keep track of the number of unique torsions encountered.
// Note that this is different from the number of columns filled
// since there will likely be collisions.
static absl::flat_hash_set<IWString> torsions_found;

/*
  As a quick check for presence of a torsion, we hash the atomic
  numbers in the target torsions
*/

static int* only_want_some = nullptr;

#define QUICK_HASH_SIZE (26 * 26 * 26 * 26)

static int
quick_torsion_hash(const const_IWSubstring& t)
{
  int h = 0;

  int istop;
  int delta;
  if (eight_character_torsion_names) {
    istop = 8;
    delta = 2;
  } else {
    istop = 12;
    delta = 3;
  }

  istop += topotorsion_prefix.length();

  // cerr << "Computing hash for '" << t << "'\n";

  for (int i = topotorsion_prefix.length(); i < istop; i += delta) {
    char c = t[i];

    //  cerr << " i = " << i << " char '" << c << "'\n";

    if ('_' == c) {
      c = 'Z';
    }

    h = 26 * h + (c - 'A');
  }

  // cerr << "Torsion '" << t << "' hash value " << h << '\n';

  assert(h >= 0 && h < QUICK_HASH_SIZE);

  return h;
}

/*
  Feb 99. Stupid SAS converts all the descriptor names to uppercase,
  which can lead to collisions.
  Use newdescriptornames to create new descriptor names and a
  cross reference file. We read that cross reference file
*/

static int
hash_only_want_torsion(const const_IWSubstring& t, IWString& header)
{
  if (0 == topotorsion_prefix.length()) {
    ;
  } else if (t.starts_with(topotorsion_prefix)) {
    ;
  } else {
    cerr << "Torsions must start with '" << topotorsion_prefix << "', '" << t
         << "' is invalid\n";
    return 0;
  }

  if (eight_character_torsion_names && (8 + topotorsion_prefix.length()) != t.length()) {
    cerr << "Wrong length for specified torsion, '" << t << "' is invalid\n";
    return 0;
  }

  if (nullptr == only_want_some) {
    only_want_some = new_int(QUICK_HASH_SIZE);
  }

  int h = quick_torsion_hash(t);

  only_want_some[h] = 1;

  torsions_to_find.add(t);

  // cerr << "Will be looking for torsion '" << t << "'\n";

  header << ' ' << t;

  return 1;
}

static int
read_cross_reference_file(iwstring_data_source& input, IWString& header)
{
  input.set_skip_blank_lines(1);
  input.set_strip_trailing_blanks(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (buffer.nwords() < 2) {
      cerr << "Bad record in cross ref file, line " << input.lines_read() << " '"
           << buffer << "'\n";
      return 0;
    }

    const_IWSubstring mangled_name;
    IWString torsion;  // needs to be an IWString for alignment reasons - see the hash
                       // function in torsion_hash

    buffer.word(0, mangled_name);
    buffer.word(1, torsion);

    if (verbose) {
      cerr << "Descriptor '" << mangled_name << "' is really torsion '" << torsion
           << "'\n";
    }

    if (!hash_only_want_torsion(torsion, header)) {
      return 0;
    }
  }

  return 1;
}

static int
read_cross_reference_file(const const_IWSubstring& fname, IWString& header)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "read_cross_reference_file: cannot open '" << fname << "'\n";
    return 0;
  }

  if (nullptr == only_want_some) {
    only_want_some = new_int(QUICK_HASH_SIZE);
  }

  return read_cross_reference_file(input, header);
}

static int
process_only_want_torsions_from_file(iwstring_data_source& input, IWString& header)
{
  input.set_skip_blank_lines(1);

  // Some logic to see that if we read multiple tokens off the first line, we don't
  // read any more lines

  int tokens_on_first_line = 0;

  IWString buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    buffer.strip_leading_blanks();

    if (buffer.matches_ignore_case("name ")) {
      buffer.remove_leading_words(1);
    }

    buffer.gsub(" TT_", " tt_");

    if (tokens_on_first_line > 0) {
      buffer.truncate_at_first(' ');
      if (!hash_only_want_torsion(buffer, header)) {
        return 0;
      }
    } else {
      int i = 0;
      const_IWSubstring t;
      while (buffer.nextword(t, i)) {
        if (!hash_only_want_torsion(t, header)) {
          return 0;
        }

        tokens_on_first_line++;
      }

      if (tokens_on_first_line > 1) {
        return 1;
      }
    }
  }

  return 1;
}

static int
process_only_want_torsions_from_file(const const_IWSubstring& fname, IWString& header)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "process_only_want_torsions_from_file: cannot open input file '" << fname
         << "'\n";
    return 0;
  }

  return process_only_want_torsions_from_file(input, header);
}

/*
  When outputting just a subset of torsions, we can get the output as an array
*/

static int* output_vector = nullptr;

static int number_torsions_wanted = 0;

static int
torsion_requested(const const_IWSubstring& t)
{
  assert(nullptr != only_want_some);

  int h = quick_torsion_hash(t);

  if (0 == only_want_some[h]) {
    return 0;
  }

  return torsions_to_find.contains(t);
}

static int
torsion_requested(const const_IWSubstring& t, int& id)
{
  assert(nullptr != only_want_some);

  int h = quick_torsion_hash(t);

  // cerr << "Do we want '" << t << "', hash " << only_want_some[h] << ", hash " <<
  // torsions_to_find.contains (t) << '\n';

  if (0 == only_want_some[h]) {
    return 0;
  }

  return torsions_to_find.fetch_unique_identifier(t, id);
}

static IWString_and_File_Descriptor stream_for_torsion_names;

/*
  We do a brute force checking of the hash function
*/

// #define CHECK_HASHING_FUNCTION
#ifdef CHECK_HASHING_FUNCTION

static IWString check_torsion_array[1000000];

int check_torsions_stored = 0;

static int
check_hash_function(int uid, const IWString& torsion)
{
  for (int i = 1; i <= check_torsions_stored; i++) {
    if (torsion == check_torsion_array[i]) {
      if (i == uid) {
        return 1;
      }

      cerr << "Hashed torsion identifier mismatch, uid = " << uid << " index " << i
           << '\n';
      abort();
    }
  }

  check_torsions_stored++;
  check_torsion_array[check_torsions_stored] = torsion;

  return 1;
}

#endif

/*
  The number of molecules for which each torsion is set.
*/

static extending_resizable_array<int> hits;

/*
  When working with hashing functions, it is sometimes useful to
  insert spaces between the components of the torsion when it is
  written out. This option is not for general use, and is not
  reported in usage();
*/

static int insert_spaces_in_torsions = 0;

typedef unsigned short tt_atom_t;

static int suppress_molecules_with_no_torsions = 1;

// The number of molecules that have a given number of torsions.
static extending_resizable_array<uint32_t> torsion_count;

// Given a torsion name, return a column number in [0,fixed_ncols)
static uint32_t
TorsionToIndex(const uint32_t fixed_ncols,
               const IWString& name) {
  auto iter = torsions_found.find(name);
  if (iter == torsions_found.end()) {
    torsions_found.emplace(name);
  }

  // Also works with an absl hash, but that is not guaranteed
  // to remain stable across releases. This seems to be as performant
  // and well behaved.
  static IWStringHash hasher;

  return hasher(name) % fixed_ncols;
}

class Topological_Torsion
{
 private:
  tt_atom_t _a1;
  tt_atom_t _a2;
  tt_atom_t _a3;
  tt_atom_t _a4;

 public:
};

// Within a tt_atom_t, the atomic number is the first 4 bits

#define tt_carbon 0x10
#define tt_nitrogen 0x20
#define tt_oxygen 0x30
#define tt_fluorine 0x40
#define tt_sulphur 0x50
#define tt_chlorine 0x60
#define tt_bromine 0x70
#define tt_iodine 0x80

// within a tt_atom_t the number of pi electrons occupies bits 4-5

#define tt_zero_pi 0x10
#define tt_one_pi 0x14
#define tt_three_pi 0x18
#define tt_four_pi 0x1c

// within a tt_atom_t the number of non-hydrogen atom branchings occupy bits 6-7

#define tt_zero_ncon 0x1c
#define tt_one_ncon 0x1d
#define tt_two_ncon 0x1e
#define tt_three_ncon 0x1f

static int
do_write_torsions(const IWString& mname, const resizable_array_p<TopoTorsion>& tt,
                  IWString_and_File_Descriptor& output)
{
  output << "Name " << mname << '\n';

  int nt = tt.number_elements();

  if (write_torsion_strings) {
    for (int i = 0; i < nt; i++) {
      output << *(tt[i]) << ' ' << tt[i]->count() << '\n';
    }
  } else if (stream_for_torsion_names.is_open()) {
    for (int i = 0; i < nt; i++) {
      const IWString& string_rep = *(tt[i]);
      int j = torsion_hash.unique_identifier(string_rep);

#ifdef CHECK_HASHING_FUNCTION
      check_hash_function(j, string_rep);
#endif

      output << j << ' ' << tt[i]->count() << '\n';

      hits[j]++;  // another molecule with this torsion
    }
  }

  output << '|' << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

static int
do_output_vector(const IWString& mname, const int* v, int n,
                 IWString_and_File_Descriptor& output)
{
  append_first_token_of_name(mname, output);

  for (int i = 0; i < n; i++) {
    output << ' ' << v[i];
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

/*
  if we suppress any torsions, we need a vector which marks those torsions to
  be ignored.
*/

static int* ignore_torsion = nullptr;

static void
filter_torsions()
{
  int low;
  if (lower_torsion_count) {
    low = lower_torsion_count;
  } else if (lower_torsion_percent) {
    low = lower_torsion_percent * molecules_read / 100;
    if (0 == low) {
      low = 1;
    }
  } else {
    low = 0;
  }

  int nt = hits.number_elements();

  int torsions_suppressed = 0;

  if (low) {
    for (int i = 0; i < nt; i++) {
      if (hits[i] <= low) {
        ignore_torsion[i] = 1;
        torsions_suppressed++;
        if (verbose) {
          cerr << "Torsion " << i << " too few hits " << hits[i] << " min " << low
               << '\n';
        }
      }
    }
  }

  int high;
  if (upper_torsion_count) {
    high = upper_torsion_count;
  } else if (upper_torsion_percent) {
    high = upper_torsion_percent * molecules_read / 100;
  } else {
    high = 0;
  }

  if (high) {
    for (int i = 0; i < nt; i++) {
      if (ignore_torsion[i]) {
        continue;
      }

      if (hits[i] >= high) {
        ignore_torsion[i] = 1;
        torsions_suppressed++;
        if (verbose) {
          cerr << "Torsion " << i << " too many hits " << hits[i] << " max " << high
               << '\n';
        }
      }
    }
  }

  return;
}

static int
WriteFixedWidthVector(Molecule& m,
                     int* output_vector,
                     uint32_t fixed_ncols,
                     const resizable_array_p<TopoTorsion>& torsions,
                     IWString_and_File_Descriptor& output) {
  append_first_token_of_name(m.name(), output);

  std::fill_n(output_vector, fixed_ncols, 0);

  for (const TopoTorsion* t : torsions) {
    uint32_t col = TorsionToIndex(fixed_ncols, *t);
    output_vector[col] += t->count();
  }

  static constexpr char kSep = ' ';

  for (uint32_t i = 0; i < fixed_ncols; ++i) {
    output << kSep << output_vector[i];
  }
  output << '\n';

  if (output.write_if_buffer_holds_more_than(4096)) {
    output.flush();
  }

  return 1;
}


static int
write_descriptor_file_header(IWString* torsion, IWString_and_File_Descriptor& output)
{
  if (verbose) {
    cerr << "Writing " << torsion_hash.size() << " torsions\n";
  }

  IW_STL_Hash_Map_int::const_iterator f;
  for (f = torsion_hash.begin(); f != torsion_hash.end(); ++f) {
    int column = (*f).second;
    if (ignore_torsion[column]) {
      continue;
    }

    //  cerr << "Torsion in column " << column << " is '" << (*f).first << "'\n";

    assert(0 == torsion[column].length());

    torsion[column] << (*f).first;
  }

  int n = torsion_hash.size();

  // Make sure every column got assigned

  int rc = 1;

  for (int i = 1; i < n; i++) {
    if (ignore_torsion[i]) {
      continue;
    }

    if (0 == torsion[i].length()) {
      cerr << "Yipes, no torsion name assigned to column " << i << '\n';
      rc = 0;
    }
  }

  if (0 == rc) {
    return 0;
  }

  // Now write them

  for (int i = 0; i < n; i++) {
    if (ignore_torsion[i]) {
      continue;
    }

    output << ' ' << torsion[i];
  }

  return output.good();
}

static int
write_descriptor_file_header(IWString_and_File_Descriptor& output)
{
  output << "Name";

  IWString* tmp = new IWString[torsion_hash.size()];
  std::unique_ptr<IWString[]> free_tmp(tmp);

  int rc = write_descriptor_file_header(tmp, output);

  output << '\n';

  return 0 == rc ? 0 : output.good();
}

static int
write_topotorsion_descriptors(const resizable_array<int>& mt,
                              const resizable_array<int>& mtc, int* tmp,
                              IWString_and_File_Descriptor& output)
{
  set_vector(tmp, torsion_hash.size(), 0);

  int n = mt.number_elements();
  assert(n == mtc.number_elements());

  for (int i = 0; i < n; i++) {
    int t = mt[i];  // which torsion

    if (ignore_torsion[t]) {
      continue;
    }

    if (mtc[i] <= 255) {
      tmp[t] = mtc[i];
    } else {
      tmp[t] = 255;
    }
  }

  int nt = torsion_hash.size();

  for (int i = 0; i < nt; i++) {
    if (ignore_torsion[i]) {
      continue;
    }

    if (0 == tmp[i]) {
      output.strncat(" 0", 2);
    } else {
      iwdigits.append_number(output, tmp[i]);
    }
  }

  output.add('\n');

  // output << output_buffer;

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

static int
write_topotorsion_descriptors(IWString_and_File_Descriptor& output)
{
  int* tmp = new int[torsion_hash.size()];
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < molecules_stored; i++) {
    const IWString& mname = stored_name[i];
    const resizable_array<int>& mt = moltor[i];
    const resizable_array<int>& mtc = moltorcount[i];

    output << mname;

    if (!write_topotorsion_descriptors(mt, mtc, tmp, output)) {
      cerr << "I/O error writing descriptors\n";
      return 0;
    }
  }

  return output.good();
}

static int
write_descriptor_file(IWString_and_File_Descriptor& os)
{
  assert(hits.number_elements() > 0);

  ignore_torsion = new_int(hits.number_elements());

  if (lower_torsion_count || lower_torsion_percent || upper_torsion_count ||
      upper_torsion_percent) {
    filter_torsions();
  }

  if (!write_descriptor_file_header(os)) {
    return 0;
  }

  return write_topotorsion_descriptors(os);
}

static int
store_torsions(const IWString& molecule_name, const resizable_array_p<TopoTorsion>& tt)
{
  int n = molecules_stored;

  assert(0 == moltor[n].number_elements());
  assert(0 == moltorcount[n].number_elements());
  assert(0 == stored_name[n].length());

  int ntt = tt.number_elements();

  moltor[n].resize(ntt);
  if (!moltorcount[n].resize(ntt)) {
    cerr << "store_torsions: memory failure, molecule " << n << " ntorsions = " << ntt
         << '\n';
    return 0;
  }

  for (int i = 0; i < ntt; i++) {
    const TopoTorsion* t = tt[i];

    int tid = torsion_hash.unique_identifier(*t);

    moltor[n].add(tid);
    moltorcount[n].add(t->count());

    hits[tid]++;  // another molecule with this torsion
  }

  stored_name[n] = molecule_name;
  stored_name[n].truncate_at_first(' ');

  molecules_stored++;

  return 1;
}

static void
append_atom(IWString& tt, Molecule& m, const IWString& asymbol, atom_number_t a, int ncon)
{
  int pi = 0;
  (void)m.pi_electrons(a, pi);

  if (numeric_torsion_names) {
    append_numeric_torsion_names(tt, m, a, pi, ncon);
  } else if (eight_character_torsion_names) {
    append_numbers_eight_chars(tt, m, asymbol, a, pi, ncon);
  } else {
    append_numbers_twelve_chars(tt, m, a, pi, ncon);
  }

  return;
}

static int
topotorsion(Molecule& m, const atomic_number_t* z, const IWString& asymbol,
            const int* ncon, atom_number_t a0, atom_number_t a1, atom_number_t a2,
            atom_number_t a3, resizable_array_p<TopoTorsion>& tt)
{
  form_canonical_order(m, z, ncon, a0, a1, a2, a3);

  TopoTorsion* t = new TopoTorsion();

  IWString& tmp = *t;

  tmp << topotorsion_prefix;

  append_atom(tmp, m, asymbol, a0, ncon[a0]);
  if (insert_spaces_in_torsions) {
    tmp << ' ';
  }
  append_atom(tmp, m, asymbol, a1, ncon[a1]);
  if (insert_spaces_in_torsions) {
    tmp << ' ';
  }
  append_atom(tmp, m, asymbol, a2, ncon[a2]);
  if (insert_spaces_in_torsions) {
    tmp << ' ';
  }
  append_atom(tmp, m, asymbol, a3, ncon[a3]);

  // we want all of them - the most common case
  if (nullptr == only_want_some) {
    tt.add(t);
    return 1;
  }

  if (fixed_ncols) {
    tt << t;
    return 1;
  }

  // We need to do things differently depending on whether we are doing
  // array output or regular output

  if (nullptr != output_vector)  // array output
  {
    int i;
    if (torsion_requested(tmp, i)) {
      output_vector[i]++;
    }

    delete t;
    return 1;
  }

  if (!torsion_requested(*t)) {
    delete t;
    return 1;
  }

  tt.add(t);

  return 1;
}

static int
topotorsion(Molecule& m, const atomic_number_t* z, const IWString& asymbol,
            const int* ncon, atom_number_t a0, atom_number_t a1, atom_number_t a2,
            resizable_array_p<TopoTorsion>& tt)
{
  int nc2 = ncon[a2];

  int rc = 0;

  const Atom* a = m.atomi(a2);

  assert(nc2 == a->ncon());

  for (int i = 0; i < nc2; i++) {
    atom_number_t a3 = a->other(a2, i);
    if (a3 == a1) {
      continue;
    }

    rc += topotorsion(m, z, asymbol, ncon, a0, a1, a2, a3, tt);
  }

  return rc;
}

static int
topotorsion(Molecule& m, const atomic_number_t* z, const IWString& asymbol,
            const int* ncon, const Bond* b, resizable_array_p<TopoTorsion>& tt)
{
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();
  int nc1 = ncon[a1];

  if (1 == nc1 || 1 == ncon[a2]) {  // A1 and A2 are in the middle of a torsion
    return 0;
  }

  const Atom* a = m.atomi(a1);

  int rc = 0;
  for (int i = 0; i < nc1; i++) {
    atom_number_t o = a->other(a1, i);
    if (o == a2) {
      continue;
    }

    rc += topotorsion(m, z, asymbol, ncon, o, a1, a2, tt);
  }

  return rc;
}

static int
topotorsion(Molecule& m, const atomic_number_t* z, const IWString& asymbol,
            const int* ncon, IWString_and_File_Descriptor& output)
{
  if (fixed_ncols) {
  } else if (nullptr != output_vector) {
    set_vector(output_vector, number_torsions_wanted, 0);
  }

  resizable_array_p<TopoTorsion> tt;

  tt.reserve(4 * m.nedges());

  for (const Bond* b : m.bond_list()) {
    topotorsion(m, z, asymbol, ncon, b, tt);
  }

  if (verbose) {
    ++torsion_count[tt.size()];
  }

  if (fixed_ncols) {
    return WriteFixedWidthVector(m, output_vector, fixed_ncols, tt, output);
  }

  if (nullptr != output_vector) {
    return do_output_vector(m.name(), output_vector, number_torsions_wanted, output);
  }

  if (verbose > 1) {
    cerr << "Found " << tt.number_elements() << " topological torsions\n";
  }

  if (0 == tt.number_elements() && suppress_molecules_with_no_torsions) {
    return 1;
  }

  if (verbose > 2) {
    for (int i = 0; i < tt.number_elements(); i++) {
      cerr << "i = " << i << " '" << *(tt[i]) << "'\n";
    }
  }

  if (tt.number_elements()) {
    sort_and_count_tts(tt);
  }

  if (produce_descriptor_file) {
    return store_torsions(m.name(), tt);
  }

  return do_write_torsions(m.name(), tt, output);
}

static int
topotorsion(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  assert(input.ok());

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    m->remove_all(1);

    const int matoms = m->natoms();

    if (0 == matoms) {
      cerr << "Skip empty molecule " << m->name() << '\n';
      continue;
    }

    int* ncon = new int[matoms];
    std::unique_ptr<int[]> free_ncon(ncon);

    (void)m->ncon(ncon);

    atomic_number_t* z = new atomic_number_t[matoms];
    std::unique_ptr<atomic_number_t[]> free_z(z);

    m->atomic_numbers(z);

    IWString asymbol;  // when producing 8 char torsions, pre-compute the element symbols

    if (eight_character_torsion_names) {
      initialise_single_character_atom_symbols(*m, asymbol);
    }

    if (!topotorsion(*m, z, asymbol, ncon, output)) {
      return 0;
    }
  }

  return output.good();
}

static int
initialise_cache(data_source_and_type<Molecule>& input)
{
  int m = input.molecules_remaining();

  if (0 == m) {
    return 0;
  }

  moltor = new resizable_array<int>[m];
  moltorcount = new resizable_array<int>[m];
  stored_name = new IWString[m];

  if (nullptr == moltorcount) {
    cerr << "initialise_cache: memory failure, cannot allocate array for " << m
         << " molecule\n";
    return 0;
  }

  if (verbose) {
    cerr << "Cache sized for " << m << " molecules\n";
  }

  return 1;
}

static int
topotorsion(const char* fname, FileType input_type, IWString_and_File_Descriptor& output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  if (produce_descriptor_file) {
    if (!initialise_cache(input)) {
      cerr << "Cannot initialise descriptor cache for '" << fname << "'\n";
      return 0;
    }
  }

  return topotorsion(input, output);
}

static void
usage(int rc)
{
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Usage : " << prog_name << " options file1 file2 file3 ....\n";
  cerr << "  -n 8           create 8 character torsion names (for SAS, the default)\n";
  cerr << "  -n 12          create 12 character torsion names (more comprehensible)\n";
  cerr << "  -n 36          a case insensitive representation using 36 chars\n";
  cerr << "  -P <prefix>    prefix for topotorsion names\n";
  cerr << "  -z             write data for molecules which don't have any torsions\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -Q <fname>     create torsion index output\n";
  cerr << "  -O <torsion>   only process <torsion>\n";
  cerr << "  -O F:<fname>   only process torsions in <fname>\n";
  cerr << "  -X <fname>     only process torsions in cross reference file <fname>\n";
  cerr << "  -a             write output as array (only valid with -O or -X option)\n";
  cerr << "  -Y             single pass, memory intensive, produces descriptor file\n";
  cerr << "  -H <ncol>      hash feature names to generate a fixed width descriptor file\n";
  cerr << "  -c <number>    suppress torsions with <number> or fewer hits\n";
  cerr << "  -c <number%>   suppress torsions with <percent>or fewer hits\n";
  cerr << "  -C <number>    suppress torsions with <number> or more hits\n";
  cerr << "  -C <number%>   suppress torsions with <percent>or more hits\n";
  cerr << "  -E <symbol>    create element with symbol\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  A valid input is '-X 10' or '-X 30%'
*/

static int
get_number_or_percent(Command_Line& cl, char flag, int& as_number, int& as_percent)
{
  const_IWSubstring x = cl.string_value(flag);

  if (0 == x.length()) {
    cerr << "get_number_or_percent: no '" << flag << "' option present\n";
    return 0;
  }

  int rc = 0;

  if (x.ends_with('%')) {
    x.chop(1);
    rc = x.numeric_value(as_percent);
  } else {
    rc = x.numeric_value(as_number);
  }

  return rc;
}

int
topotorsion(int argc, char** argv)
{
  Command_Line cl(argc, argv, "sn:A:DE:vi:o:zQ:O:X:aYc:C:fJ:P:qH:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  iwdigits.set_include_leading_space(1);

  iwdigits.initialise(255);

  process_elements(cl);

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(7);
  }

  if (cl.option_present('J')) {
    tag = cl.string_value('J');
    if (verbose) {
      cerr << "Torsions written as non-colliding sparse fingerprints, tag '" << tag
           << "'\n";
    }

    if (!tag.ends_with('<')) {
      tag += '<';
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('f')) {
    function_as_filter = 1;
  } else if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      return 4;
    }
  } else if (!process_input_type(cl, input_type)) {
    cerr << prog_name << ": cannot discern input type\n";
    usage(2);
  }

  if (cl.option_present('z')) {
    suppress_molecules_with_no_torsions = 0;
    if (verbose) {
      cerr << "Will write even molecules with no torsions\n";
    }
  }

  if (cl.option_present('P')) {
    cl.value('P', topotorsion_prefix);

    if (verbose) {
      cerr << "All torsions must start with '" << topotorsion_prefix << "'\n";
    }
  }

  if (cl.option_present('n')) {
    const_IWSubstring nchars;
    cl.value('n', nchars);
    if ('8' == nchars) {
      eight_character_torsion_names = 1;
      if (verbose) {
        cerr << "Eight character torsion names will be generated\n";
      }
    } else if ("12" == nchars) {
      eight_character_torsion_names = 0;
      if (verbose) {
        cerr << "Twelve character torsion names will be generated\n";
      }
    } else if ("36" == nchars) {
      numeric_torsion_names = 1;
      eight_character_torsion_names = 0;
      if (verbose) {
        cerr << "Will write torsions as numbers\n";
      }
    } else {
      cerr << "Unrecognised -n qualifier '" << nchars << "'\n";
      usage(11);
    }
  }

  if (cl.option_present('s')) {
    if (eight_character_torsion_names) {
      cerr << "Sorry, the -s option only works with long torsion names (-n 12)\n";
      return 8;
    }

    insert_spaces_in_torsions = 1;
    if (verbose) {
      cerr << "Torsions will have spaces included\n";
    }
  }

  if (cl.option_present('c')) {
    if (!get_number_or_percent(cl, 'c', lower_torsion_count, lower_torsion_percent)) {
      cerr << "Cannot parse lower torsion count option (-c)\n";
      usage(17);
    }

    if (verbose) {
      cerr << "Lower torsion counts " << lower_torsion_count << ", and "
           << lower_torsion_percent << " percent\n";
    }
  }

  if (cl.option_present('C')) {
    if (!get_number_or_percent(cl, 'C', upper_torsion_count, upper_torsion_percent)) {
      cerr << "Cannot parse upper torsion count option (-C)\n";
      usage(17);
    }

    if (verbose) {
      cerr << "Upper torsion counts " << upper_torsion_count << ", and "
           << upper_torsion_percent << " percent\n";
    }
  }

  if (cl.empty()) {
    cerr << prog_name << ": insufficient arguments\n";
    usage(3);
  }

  if (cl.option_present('Q')) {
    const char* q = cl.option_value('Q');

    if (!stream_for_torsion_names.open(q)) {
      cerr << "Cannot open torsion index file '" << q << "'\n";
      return 31;
    }

    if (verbose) {
      cerr << "Torsion index file '" << q << "' being created\n";
    }

    write_torsion_strings = 0;
    write_torsion_hash = 1;
  }

  if (cl.option_present('O') && cl.option_present('X')) {
    cerr << "The -O and -X options are mutually exclusive\n";
    usage(42);
  }

  if (cl.option_present('Y') && (cl.option_present('O') || cl.option_present('X'))) {
    cerr << "The single pass (-Y) option is incompatible with the -O or -X options\n";
    usage(33);
  }

  // If we are doing either -O or -X, we can write the results as an array

  if (cl.option_present('a')) {
    if (!cl.option_present('O') && !cl.option_present('X')) {
      cerr << "The -a option can only be used with the -O or -X options\n";
      usage(19);
    }
  }

  IWString_and_File_Descriptor output(1);
  output.reserve(8192);

  IWString header("Name");
  if (cl.option_present('O')) {
    IWString o;  // alignment problems on the Sun if we use const_IWSubstring!
    int i = 0;
    while (cl.value('O', o, i++)) {
      if (o.starts_with("F:")) {
        o.remove_leading_chars(2);
        if (!process_only_want_torsions_from_file(o, header)) {
          return i + 1;
        }
      } else {
        if (!hash_only_want_torsion(o, header)) {
          return i + 1;
        }

        if (verbose) {
          cerr << "Will process torsion '" << o << "'\n";
        }
      }
    }

    output << header << '\n';

    number_torsions_wanted = torsions_to_find.torsions_stored();

    if (verbose > 1) {
      for (IW_STL_Hash_Map_int::const_iterator i = torsions_to_find.begin();
           i != torsions_to_find.end(); ++i) {
        cerr << "Torsion '" << (*i).first << "' is number " << (*i).second << '\n';
      }
    }
  }

  if (cl.option_present('X')) {
    IWString header("Name");

    const_IWSubstring x = cl.string_value('X');

    if (!read_cross_reference_file(x, header)) {
      cerr << "Error processing cross reference file '" << x << "'\n";
      return 54;
    }

    output << header << '\n';

    number_torsions_wanted = torsions_to_find.torsions_stored();
  }

  if (cl.option_present('Y')) {
    if (cl.number_elements() > 1) {
      cerr << "Sorry, cannot produce a descriptor file with multiple files present\n";
      usage(41);
    }

    produce_descriptor_file = 1;

    if (verbose) {
      cerr << "Will produce a descriptor file\n";
    }
  }

  if (cl.option_present('H')) {
    if (cl.option_present('Y') || cl.option_present('c') || cl.option_present('C') ||
        cl.option_present('Q') || cl.option_present('O') || cl.option_present('X') ||
        cl.option_present('a')) {
      cerr << "The -Y -c -C -Q -O -X -a options are not compatible with the -H option\n";
      return 1;
    }

    if (! cl.value('H', fixed_ncols) || fixed_ncols < 10) {
      cerr << "The -H option must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will produce a fixed width descriptor file with " << fixed_ncols << " columns\n";
    }
    
    output_vector = new int[fixed_ncols];

    static constexpr char kSep = ' ';
    output << "Name";
    for (uint32_t i = 0; i < fixed_ncols; ++i) {
      output << kSep << "tt_TT" << i;
    }
    output << '\n';
  }

  if (cl.option_present('a')) {
    output_vector = new int[number_torsions_wanted];

    if (verbose) {
      cerr << number_torsions_wanted << " selected torsions output as an array\n";
    }
  }

  if (produce_descriptor_file && verbose) {
    cerr << "Computation started at " << time(NULL) << '\n';
  }

  for (const char * fname: cl) {
    if (verbose) {
      cerr << prog_name << " processing '" << fname << "'\n";
    }

    if (!topotorsion(fname, input_type, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    cerr << molecules_read << " molecules read";
    if (torsion_hash.size()) {
      cerr << ", found " << torsion_hash.size() << " torsions";
    }
    cerr << '\n';
    Accumulator_Int<uint32_t> acc;
    for (uint32_t i = 0; i < torsion_count.size(); ++i) {
      if (torsion_count[i]) {
        cerr << torsion_count[i] << " molecules had " << i << " torsions\n";
        acc.extra(i, torsion_count[i]);
      }
    }
    cerr << "Mean " << acc.average() << '\n';
    if (fixed_ncols) {
      cerr << "NUmber of different torsions " << torsions_found.size() << '\n';
    }
  }

  if (write_torsion_hash && molecules_read) {
    torsion_hash.do_write(stream_for_torsion_names);

    stream_for_torsion_names << "|\n";

    for (int i = 0; i < hits.number_elements(); i++) {
      stream_for_torsion_names << i << ' ' << hits[i] << '\n';
      stream_for_torsion_names.write_if_buffer_holds_more_than(32768);
    }

    stream_for_torsion_names.flush();
  }

  if (produce_descriptor_file && molecules_read) {
    if (verbose) {
      cerr << "Starting writing descriptor file at " << time(NULL) << '\n';
    }
    write_descriptor_file(output);
  }

  if (cl.option_present('q')) {
    _exit(0);
  }

  if (nullptr != output_vector) {
    delete[] output_vector;
  }

  if (nullptr != only_want_some) {
    delete[] only_want_some;
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = topotorsion(argc, argv);

  return rc;
}
