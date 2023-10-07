// Atom triple fingerprints. Augmented atom pairs.
// Has not been fully explored, but seemed to do OK in a couple
// of models.

#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/fingerprint_writer.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int function_as_filter = 0;

static int molecules_read = 0;

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int min_distance = 0;
static int max_distance = std::numeric_limits<int>::max();

static Atom_Typing_Specification atom_typing_specification;

static int form_triangles = 0;

static int perform_canonicalisation = 1;

static int ntest = 0;

static Report_Progress report_progress;

static int allow_linear_arrangements = 1;

static int additive = 0;

static int area_based_atom_triangles = 0;

static fingerprint_writer::FingerprintWriter fp_writer;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Atom triple fingerprints, 3 atoms and the distances between them\n";
  cerr << "  -P <type>     atom typing specification\n";
  cerr << "  -f            function as a filter\n";
  cerr << "  -J <tag>      tag for fingerprints, enter '-J help' for info\n";
  cerr << "  -a            create the fingerprints as triangles rather than triples\n";
  cerr << "  -h            set bits based on triangle areas - input must be 3d\n";
  cerr << "  -u            do NOT perform canonicalisation\n";
  cerr << "  -t <n>        perform <n> tests on each molecule\n";
  cerr << "  -r <n>        when running tests, report progress every <n> molecules processed\n";
  cerr << "  -z            exclude linear arrangements of atoms\n";
  cerr << "  -d            bits are formed additively - less precise\n";
  cerr << "  -c            minimum distance allowed\n";
  cerr << "  -C <dist>     maximum distance allowed\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
write_fingerprint(const Molecule& m,
                  const Sparse_Fingerprint_Creator& sfc,
                  IWString_and_File_Descriptor& output)
{
#ifdef CHECK_OVER_255
  const auto bits_found = sfc.bits_found();
  int over_255 = 0;
  for (auto i : bits_found) {
    if (i.second > 254) {
      cerr << " BIt " << i.first << " count " << i.second << endl;
      over_255++;
    }
  }
  if (over_255) {
    cerr << over_255 << " of " << bits_found.size() << " bits over 255\n";
  }
#endif

  fp_writer.WriteFingerprint(m.name(), sfc, output);

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

/*
  How to sort 3 items. They will come in as one of the following
  123
  132
  213
  232
  312
  321
  Swap first and last if needed
  123
  132
  213
  231   132
  312   213
  321   123
  We see that the 6 cases reduce to just 3. Those are then handled by two if statements.
*/

/*
  For the area, we canonicalise on the atom types. The two lowest atom types form the base
  and the highest atom type forms the peak.
*/

static void
do_area_based_atom_triangles(atom_number_t i, atom_number_t j, atom_number_t k,
                             const int* atype, int dij, int dik, int djk,
                             Sparse_Fingerprint_Creator& sfc)
{
  auto ati = atype[i];
  auto atj = atype[j];
  auto atk = atype[k];

  if (ati > atk) {
    std::swap(i, k);
    std::swap(ati, atk);
    std::swap(djk, dij);
  }

  if (atj > atk) {
    std::swap(j, k);
    std::swap(atj, atk);
    std::swap(dij, dik);
  } else if (ati > atj) {
    std::swap(i, j);
    std::swap(ati, atj);
    std::swap(dik, djk);
  }

  assert(ati <= atj && atj <= atk);

  // Use Heron's Formula (aka Hero's formula).
  // http://www.mathopenref.com/heronsformula.html

  const double p = static_cast<double>(dij + dik + djk) / static_cast<double>(2.0);

  const double area = sqrt(p * (p - dij) * (p - dik) * (p - djk));

  // linear path, not a triangle
  if (fabs(area) < 1.0e-05) {
    return;
  }

  const unsigned int b = ati + atj + atk;

  // cerr << "Area " << area << endl;
  // if (area < 1.0)
  //   cerr << dij << ' ' << dik << ' ' << djk << " p = " << p << endl;

  // We take the square root to try to avoid overflow on the count

  sfc.hit_bit(b, static_cast<int>(sqrt(area) + 1.0));

  return;
}

/*
  When forming triangle or atom triple fingerprints, we need to identify an edge
  For the triangle fingerprint, this will be the longest side.
  For the atom triple fingerprint, this edge will be omitted.

  We make sure that the edge to be removed is betwen atoms J and K

  First check to see if there is one edge that is longer than the other two.

  Then if there is an edge that is shorter than the other two
*/

static void
identify_canonical_edge(atom_number_t& i, atom_number_t& j, atom_number_t& k,
                        const int* atype, int& dij, int& dik, int& djk)
{
  if (djk > dij && djk > dik) {
    return;
  }

  if (dij > dik && dij > djk) {
    std::swap(i, k);
    std::swap(dij, djk);
    return;
  }

  if (dik > dij && dik > djk) {
    std::swap(i, j);
    std::swap(dik, djk);
    return;
  }

  // Is there one side that is shorter than the others.

  if (djk < dik && djk < dij) {
    return;
  }

  if (dij < dik && dij < djk) {
    std::swap(i, k);
    std::swap(dij, djk);
    return;
  }

  if (dik < dij && dik < djk) {
    std::swap(i, j);
    std::swap(dik, djk);
    return;
  }

  // At this stage there appears to be an equilateral triangle.

  if (atype[i] == atype[j] && atype[j] == atype[k]) {  // perfect equilateral triangle
    return;
  }

  // There are at least two different atom types

  // See if we can resolve by atom types. If two atom types are the same, that defines the
  // edge

  if (atype[j] == atype[k]) {
    return;
  }

  if (atype[i] == atype[j]) {
    std::swap(i, k);
    std::swap(dij, djk);
    return;
  }

  if (atype[i] == atype[k]) {
    std::swap(i, j);
    std::swap(dik, djk);
    return;
  }

  // Looks like we have three different atom types. These sums must be distinct

  const auto eij = atype[i] + atype[j];
  const auto eik = atype[i] + atype[k];
  const auto ejk = atype[j] + atype[k];

  // assert (eij != eik && eij != ejk && eik != ejk);
  if (eij != eik && eij != ejk && eik != ejk) {
    ;
  } else {
    cerr << "eij " << eij << " eik " << eik << " ejk " << ejk << ", types " << atype[i]
         << ' ' << atype[j] << " " << atype[k] << endl;
  }

  if (ejk < eij && ejk < eik) {
    return;
  }

  if (eij < eik && eij < ejk) {
    std::swap(i, k);
    return;
  }

  if (eik < eij && eik < ejk) {
    std::swap(i, j);
    return;
  }

  cerr << "Canonicalisation failed, atoms " << i << " type " << atype[i] << ' ' << j
       << " type " << atype[j] << ' ' << k << " type " << atype[k] << endl;
  return;
}

// #define DEBUG_SET_BIT_TRIPLE

static void
set_bit_triple(const atom_number_t i, const atom_number_t j, const atom_number_t k,
               const int* atype, const int dij, const int dik, const int djk,
               Sparse_Fingerprint_Creator& sfc)
{
#ifdef DEBUG_SET_BIT_TRIPLE
  cerr << "Using " << atype[i] << ' ' << atype[j] << ' ' << atype[k] << " dists: dij "
       << dij << " dik " << dik << " djk " << djk << endl;
#endif

  unsigned int b = (4 * atype[i] + atype[j] + atype[k]) + 5121 * (dij + dik);

  sfc.hit_bit(b, 1);

#ifdef DEBUG_SET_BIT_TRIPLE
  cerr << "Hit bit " << b << endl;
#endif

  return;
}

static void
set_bit_triangle(const atom_number_t i, const atom_number_t j, const atom_number_t k,
                 const int* atype, const int dij, const int dik, const int djk,
                 Sparse_Fingerprint_Creator& sfc)

{
  unsigned int b;
  if (additive) {
    b = (atype[i] + atype[j] + atype[k]) + 910 * (dij + dik + djk);
  } else {
    b = (1421021 * atype[i] + 71212 * atype[j] + 13111 * atype[k] + dij) * dik + 23 * djk;
  }

  sfc.hit_bit(b, 1);

  return;
}

static void
set_bit_triple_no_canonicalisation(const atom_number_t i, const atom_number_t j,
                                   const atom_number_t k, const int* atype, const int dij,
                                   const int dik, const int djk,
                                   Sparse_Fingerprint_Creator& sfc)
{
  unsigned int b;
  if (additive) {
    b = (atype[i] + atype[j] + atype[k]) + 5121 * (dij + dik);
  } else {
    b = (11 * atype[i] + 111 * atype[j] + 2012 * atype[k]) * (dij + 29 * dik);
  }

  sfc.hit_bit(b, 1);

  return;
}

static void
set_bit_triangle_no_canonicalisation(const atom_number_t i, const atom_number_t j,
                                     const atom_number_t k, const int* atype,
                                     const int dij, const int dik, const int djk,
                                     Sparse_Fingerprint_Creator& sfc)

{
  unsigned int b = (atype[i] + atype[j] + atype[k]) + 910 * (dij + dik + djk);

  sfc.hit_bit(b, 1);

  return;
}

static void
set_bit(atom_number_t i, atom_number_t j, atom_number_t k, const int* atype, int dij,
        int dik, int djk, Sparse_Fingerprint_Creator& sfc)
{
  if (area_based_atom_triangles) {
    do_area_based_atom_triangles(i, j, k, atype, dij, dik, djk, sfc);
    return;
  }

  if (!perform_canonicalisation) {
    if (form_triangles) {
      set_bit_triangle_no_canonicalisation(i, j, k, atype, dij, dik, djk, sfc);
    } else {
      identify_canonical_edge(i, j, k, atype, dij, dik, djk);
      set_bit_triple_no_canonicalisation(i, j, k, atype, dij, dik, djk, sfc);
    }

    return;
  }

#ifdef DEBUG_SET_BIT_TRIPLE
  cerr << "Start: " << i << ' ' << j << ' ' << k << " dij " << dij << " dik " << dik
       << " djk " << djk << endl;
#endif

  identify_canonical_edge(i, j, k, atype, dij, dik, djk);

#ifdef DEBUG_SET_BIT_TRIPLE
  cerr << "After: " << i << ' ' << j << ' ' << k << " dij " << dij << " dik " << dik
       << " djk " << djk << endl;
#endif

  if (form_triangles) {
    set_bit_triangle(i, j, k, atype, dij, dik, djk, sfc);
  } else {
    set_bit_triple(i, j, k, atype, dij, dik, djk, sfc);
  }

  return;
}

static int
ok_set_of_distances(const int dij, const int dik, const int djk)
{
  if (!allow_linear_arrangements) {
    if (dij + dik == djk || dij + djk == dik || dik + djk == dij) {
      return 0;
    }
  }

  if (dik > max_distance || djk > max_distance) {  // dij is checked elsewhere
    return 0;
  }

  if (0 == min_distance) {
    ;
  } else if (dij < min_distance || djk < min_distance) {
    return 0;
  }

  return 1;
}

static void
form_fingerprint(Molecule& m, const int* atype, Sparse_Fingerprint_Creator& sfc)
{
  const auto matoms = m.natoms();

  const auto dm = m.distance_matrix_warning_may_change();

  for (auto i = 0; i < matoms; i++) {
    assert(nullptr != dm);
    for (auto j = i + 1; j < matoms; ++j) {
      const auto dij = dm[i * matoms + j];

      if (dij > max_distance || dij < min_distance) {
        continue;
      }

      for (auto k = j + 1; k < matoms; ++k) {
        if (k == i || k == j) {
          continue;
        }
        const auto dik = dm[i * matoms + k];
        const auto djk = dm[j * matoms + k];

        if (!ok_set_of_distances(dij, dik, djk)) {
          continue;
        }

        set_bit(i, j, k, atype, dij, dik, djk, sfc);
      }
    }
  }

  return;
}

static int
do_tests(Molecule& m, int* atype)
{
  Sparse_Fingerprint_Creator sfc;

  // m.debug_print(cerr);
  form_fingerprint(m, atype, sfc);

  IWString initial_smiles = m.smiles();

  if (verbose > 1) {
    cerr << "Begin " << ntest << " tests\n";
  }

  for (auto i = 0; i < ntest; ++i) {
    IWString s = m.random_smiles();
    //  cerr << "Random smiles '" << s << "'\n";

    Molecule mtmp;

    if (!mtmp.build_from_smiles(s)) {
      cerr << "Yipes, cannot build from random smiles '"
           << "', molecule '" << m.name() << "'\n";
      return 0;
    }

    atom_typing_specification.assign_atom_types(mtmp, atype);

    Sparse_Fingerprint_Creator tmpsfc;

    form_fingerprint(mtmp, atype, tmpsfc);

    if (sfc == tmpsfc) {
      continue;
    }

    const auto& f1 = sfc.bits_found();
    const auto& f2 = tmpsfc.bits_found();

    cerr << "Fingerprint creation mismatch, reference " << f1.size() << " bits, copy "
         << f2.size() << ", molecule '" << m.name() << "'\n";
    cerr << initial_smiles << endl;
    cerr << s << endl;

    sfc.debug_print(cerr);
    tmpsfc.debug_print(cerr);
    mtmp.debug_print(cerr);

    return 0;
  }

  if (report_progress()) {
    cerr << "Completed testing " << molecules_read << " molecules\n";
  }

  return 1;
}

/*
  All routes come here
*/

static int
atom_triples(Molecule& m, IWString_and_File_Descriptor& output)
{
  m.recompute_distance_matrix();

  int* atype = new int[m.natoms()];
  std::unique_ptr<int[]> free_atype(atype);

  if (!atom_typing_specification.assign_atom_types(m, atype)) {
    cerr << "Cannot assign atom types, '" << m.name() << "'\n";
    return 0;
  }

  if (ntest) {
    if (!do_tests(m, atype)) {
      return 0;
    }

    return 1;
  }

  Sparse_Fingerprint_Creator sfc;

  form_fingerprint(m, atype, sfc);

  return write_fingerprint(m, sfc, output);
}

static int
atom_triples(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (0 == ntest) {
      output << smiles_tag << m->smiles() << ">\n";
      output << identifier_tag << m->name() << ">\n";
    }

    if (!atom_triples(*m, output)) {
      return 0;
    }

    if (0 == ntest) {
      output << "|\n";
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
atom_triples_filter(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output)
{
  const_IWSubstring smiles(buffer);

  assert(buffer.ends_with('>'));
  smiles.chop();
  smiles.remove_up_to_first('<');

  Molecule m;

  if (!m.build_from_smiles(smiles)) {
    cerr << "Cannot interpret smiles '" << smiles << "'\n";
    return 0;
  }

  return atom_triples(m, output);
}

static int
atom_triples_filter(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!atom_triples_filter(buffer, output)) {
      cerr << "Fatal error '" << buffer << "' line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
atom_triples_filter(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return atom_triples_filter(input, output);
}

static int
atom_triples(const char* fname, FileType input_type, IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return atom_triples(input, output);
}

static int
atom_triples(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lfJ:P:aut:e:r:zc:C:dh");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  set_use_version_2_augmented_atom_algorithm(1);

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot build atom typing specification (-P " << p << ")\n";
      return 2;
    }
  } else {
    atom_typing_specification.set_atom_type(IWATTYPE_Z);
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', ntest) || ntest < 1) {
      cerr << "The number of tests to perform (-t) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "For each molecule will perform " << ntest << " tests\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', min_distance) || min_distance < 1) {
      cerr << "The minimum distance value (-c) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will ignore arrangements where any atoms are closer than " << min_distance
           << " bonds\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_distance) || max_distance < min_distance) {
      cerr << "The maximum distance value (-C) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will ignore arrangements where any atoms are further than " << max_distance
           << " bonds\n";
    }
  }

  if (cl.option_present('a')) {
    form_triangles = 1;

    if (verbose) {
      cerr << "Fingerprints will be triangular\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting tool\n";
      usage(2);
    }
  }

  if (cl.option_present('h')) {
    area_based_atom_triangles = 1;
    allow_linear_arrangements = 0;

    if (verbose) {
      cerr << "Fingerprints will represent triangle areas\n";
    }
  }

  if (cl.option_present('e')) {
    random_number_seed_t seed;

    cl.value('e', seed);

    set_smiles_random_number_seed(seed);
  }

  if (cl.option_present('d')) {
    additive = 1;

    if (verbose) {
      cerr << "Bit numbers will be computed additively\n";
    }
  }

  if (cl.option_present('z')) {
    allow_linear_arrangements = 0;

    if (verbose) {
      cerr << "Will exclude linear arrangements of atoms\n";
    }
  }

  if (cl.option_present('u')) {
    perform_canonicalisation = 0;

    if (verbose) {
      cerr << "No canonicalisation performed\n";
    }
  }

  if (ntest) {
    ;
  } else if (!cl.option_present('J')) {
    cerr << "Must specify fingerprint tag via the -J option\n";
    usage(2);
  } else {
    if (! fp_writer.Initialise(cl, 'J', verbose)) {
      cerr << "Cannot initialise fingerprint (-J)\n";
      return 1;
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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (function_as_filter) {
    if (!atom_triples_filter(cl[0], output)) {
      rc = 1;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!atom_triples(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = atom_triples(argc, argv);

  return rc;
}
