/*
  Implementation of Ghose & Crippen clogp
  Is it from one or both of these papers:
    J. Comput. Chem. 7, 565, 1986
    J. Chem. Inf. Comput. Sci. 27, 21-35, 1987
*/

#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static Chemical_Standardisation chemical_standardisation;

static Charge_Assigner charge_assigner;

static IWString descriptor_name("ghose_ghcrlp");

static int reduce_to_largest_fragment = 0;

static int max_unclassified_atoms = 0;

static int write_to_stdout = 1;

/*
  Not sure we would ever want this. Note that the header generation never
  does this
*/

static int include_atoms_not_hit_and_confidence_level_in_descriptors = 0;

static IWDigits digits;

/*
  For automated processing we can insert a count of unclassified atoms
  in the output
*/

static int append_unclassified_count = 0;

static int use_query_name_as_descriptor_name = 0;

static extending_resizable_array<int> unclassified_atom_count;

static int write_as_array = 0;

/*
  The frequency vector makes a good descriptor
*/

static int include_frequency_vector_with_output = 0;

/*
  When testing, we often have the experimental value embedded in the
  name field.
*/

static int experimental_column = -1;

static Molecule_Output_Object stream_for_labeled_atoms;

/*
  There are several ways of getting atom type descriptors from this.
  The way implemented here is to populate a histogram of values
*/

static int do_histogram = 0;
static double histogram_start = 0.0;
static double histogram_stop = 0.0;
static double histogram_dx = 0.0;

/*
  There are 3 types of fingerprints we can write.
    1. Bit for each query
    2. Bit for each ISO= type (the label [] array)
    3. Bit for each bucket in a histogram
*/

static IWString bit_for_each_query_tag;
static IWString atom_label_tag;
static IWString histogram_tag;

/*
  If we are writing the frequency vector, we need an array for that.
*/

static int* frequency_vector = nullptr;

static std::ofstream stream_for_frequency_vector;

/*
  Sometimes it is helpful to have the labelled molecule in the fingerprint
  stream
*/

static IWString write_labelled_molecule_to_fingerprint;

/*
  By default, the isotopic labels are just the line number in the file.
  If the control file contains 'ISO=nnn' we set that for the isotopic label
*/

static resizable_array<int> label;

/*
  Sometimes it is convenient to know the highest numbered label present
*/

static int highest_atom_label = 0;

/*
  Each atom type can have a confidence level.
  These range from 0 (just a guess) to 100 (unambiguous)
*/

static resizable_array<int> confidence;

/*
  We allow a variable number of atoms to match per query.
*/

static resizable_array<int> atoms_to_match;

/*
  We can also insist on a minimum confidence level for each estimate
*/

static int minimum_confidence_level_needed = 0;

/*
  We keep track of (minimum) confidence levels encountered
*/

static extending_resizable_array<int> min_confidence_found;

/*
  We also offer the option of NOT marking the atoms as processed. this
  allows the possibility of multiple queries hitting the same atoms
*/

static resizable_array<int> mark_atoms;

// With the -H option, we make implicit hydrogens explicit

static int make_implicit_hydrogens_implicit = 0;

/*
  By default, we are an atom additivity method, and so just process
  the first atom in any embedding
*/

static int max_atoms_per_embedding_to_process = 1;

/*
  We can do a number of things with our output.
  We can read regular structure files, or we can be part of a pipeline.
*/

static IWString tag;

static int function_as_filter = 0;

static double intercept = 0.0;

static double multiplier = 1.0;

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
  cerr << "Usage: " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << " Implements atom based substituent models\n";
  cerr << "  -F <file>      specify control file of queries and values\n";
  cerr << "                 smarts value comments\n";
  cerr << "                 f:file.qry value comments\n";
  cerr << "  -u <number>    allow at most <number> unclassified atoms (default 0)\n";
  cerr << "  -c <number>    minimum confidence level needed (0, 100), default " << minimum_confidence_level_needed << '\n';
  cerr << "  -y             append unclassified atom count and confidence level to output\n";
  cerr << "  -T <tag>       function as filter, write results as dataitem <tag>\n";
  cerr << "  -a             output in array form - repeat to include frequency vector\n";
  cerr << "  -D <dataitem>  specify descriptor name (default '" << descriptor_name << "')\n";
  cerr << "  -G ...         various options for fingerprints, enter '-G help' for info\n";
  cerr << "  -f             function as a filter\n";
  cerr << "  -h             make implicit hydrogens explicit\n";
  cerr << "  -L <fname>     write molecules with labeled atoms to <fname>\n";
  cerr << "  -o <type>      specify output file type(s) for -L file\n";
  cerr << "  -m <number>    max atoms per query match to process (default 1)\n";
  cerr << "  -m all         process all atoms of each query match\n";
  cerr << "  -J ...         miscellaneous options, enter '-J help' for info\n";
  cerr << "  -R <file>      write atom type frequency vector to <file>\n";
  cerr << "  -S <fname>     specify output stream (default is stdout)\n";
  cerr << "  -I <intercept> intercept\n";
  cerr << "  -M <factor>    multiply all results by <factor>\n";
  cerr << "  -U             unique embeddings only\n";
  (void) display_standard_charge_assigner_options(cerr, 'N');
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static void
display_histogram_options(std::ostream& os)
{
  os << "                 Fingerprints are either bits according to queries hit\n";
  os << "                 or binned values of the atom contributions\n";
  os << "  -G fptag=<tag> fingerprint derived from one bit per query\n";
  os << "  -G altag=<tag> fingerprint with one bit per ISO= value\n";
  os << "  -G buckets=b   number of buckets to use in the fingerprint\n";
  os << "  -G hist=m,m,d  histogram for fingerprint, syntax 'min,max,dx'\n";

  return;
}

static void
display_dash_j_qualifiers(std::ostream& os)
{
  os << " -J blki         aromatic bonds lose their Kekule identity\n";
  os << " -J qndn         use query names as descriptor names\n";

  return;
}

static int
write_fingerprint(int nbits, const int* b, const IWString& tag, std::ostream& output)
{
  IW_Bits_Base fp;

  if (!fp.construct_from_array_of_ints(b, nbits)) {
    cerr << "Very strange, cannot create fingerprint from " << nbits << " bits\n";
    return 0;
  }

  return fp.write_daylight_ascii_representation(output, tag);
}

static int
write_atom_label_fingerprint(int matoms, const int* already_hit, const IWString& tag,
                             std::ostream& output)
{
  int* atom_type_vector = new_int(highest_atom_label + 1);
  std::unique_ptr<int[]> free_atom_type_vector(atom_type_vector);

  for (int i = 0; i < matoms; i++) {
    if (0 == already_hit[i]) {
      continue;
    }

    int j = label[already_hit[i]];

    atom_type_vector[j]++;
  }

  return write_fingerprint(highest_atom_label, atom_type_vector, tag, output);
}

static int
write_bit_for_each_query(int matoms, const int* already_hit, const IWString& tag,
                         std::ostream& output)
{
  int* atom_type_vector = new_int(queries.number_elements());
  std::unique_ptr<int[]> free_atom_type_vector(atom_type_vector);

  for (int i = 0; i < matoms; i++) {
    int j = already_hit[i];

    if (0 == j) {
      continue;
    }

    atom_type_vector[j]++;
  }

  return write_fingerprint(queries.number_elements(), atom_type_vector, tag, output);
}

/*
  There are various kinds of fingerprints we can write
*/

static int
write_fingerprint(Molecule& m, const IWString& smiles, const int* already_hit,
                  const IWHistogram& histogram, std::ostream& output)
{
  output << "$SMI<" << smiles << ">\n";
  if (m.name().length()) {
    output << "PCN<" << m.name() << ">\n";
  }

  if (write_labelled_molecule_to_fingerprint.length()) {
    output << write_labelled_molecule_to_fingerprint << m.smiles() << ">\n";
  }

  if (atom_label_tag.length()) {
    write_atom_label_fingerprint(m.natoms(), already_hit, atom_label_tag, output);
  }

  if (bit_for_each_query_tag.length()) {
    write_bit_for_each_query(m.natoms(), already_hit, bit_for_each_query_tag, output);
  }

  if (histogram_tag.length()) {
    histogram.write_as_daylight_fingerprint(histogram_tag, output);
  }

  if (!function_as_filter) {
    output << "|\n";
  }

  return output.good();
}

static int
write_histogram(const IWString& mname, std::ostream& output, double sum,
                const IWHistogram& histogram)
{
  write_space_suppressed_string(mname, output);

  int nb = histogram.nbuckets();

  IWString output_buffer;
  output_buffer.resize(nb * 4);

  output_buffer << ' ';
  output_buffer.append_number(sum, 4);

  const unsigned int* c = histogram.raw_counts();

  for (int i = 0; i < nb; i++) {
    digits.append_number(output_buffer, c[i]);
  }

  output << output_buffer << '\n';

  return output.good();
}

static int
write_histogram_header(std::ostream& output, const IWString& descriptor_name,
                       int nbuckets)
{
  output << "Name " << descriptor_name;

  for (int i = 0; i < nbuckets; i++) {
    output << " H" << descriptor_name[0] << i;
  }

  output << '\n';

  return output.good();
}

static int
write_the_identifier(const IWString& mname, std::ostream& output)
{
  write_first_token_of_string(mname, output);

  if (experimental_column < 0) {
    return output.good();
  }

  const_IWSubstring e;
  if (!mname.word(experimental_column, e)) {
    cerr << "Cannot extract token " << (experimental_column + 1) << " from '" << mname
         << "'\n";
    return 0;
  }

  output << ' ' << e;

  return output.good();
}

static void
append_frequency_vector_to_string(const int* frequency_vector, int nq, IWString& buffer)
{
  buffer.resize(buffer.length() + nq * 3);

  for (int i = 0; i < nq; i++) {
    digits.append_number(buffer, frequency_vector[i]);
  }

  return;
}

static int
write_frequency_vector(const int* frequency_vector, int nq, std::ostream& os)
{
  IWString output_buffer;

  append_frequency_vector_to_string(frequency_vector, nq, output_buffer);

  os << output_buffer;

  return os.good();
}

static int
apply_isotopic_labels(Molecule& m, const int* already_hit)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    if (0 == already_hit[i]) {
      continue;
    }

    m.set_isotope(i, label[already_hit[i] - 1]);
  }

  return 1;
}

static int
write_unclassified_atoms(Molecule& m, const int* already_hit, std::ostream& output)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (0 == already_hit[i]) {
      output << "Atom " << i << " " << m.smarts_equivalent_for_atom(i) << " not hit\n";
    }
  }

  return 1;
}

static int
too_many_unclassified_atoms_or_uncertainty(Molecule& m, const int* already_hit,
                                           int atoms_not_hit,
                                           int minimum_confidence_level)
{
  unclassified_atom_count[atoms_not_hit]++;

  if (atoms_not_hit)  // maybe problems
  {
    cerr << m.name() << ' ' << atoms_not_hit << " of " << m.natoms()
         << " atoms not classified\n";

    if (verbose) {
      write_unclassified_atoms(m, already_hit, cerr);
    }

    if (max_unclassified_atoms >= 0 && atoms_not_hit > max_unclassified_atoms) {
      return 1;
    }
  }

  if (minimum_confidence_level_needed &&
      minimum_confidence_level < minimum_confidence_level_needed) {
    cerr << m.name() << " confidence level too low " << minimum_confidence_level << '\n';
    return 1;
  }

  if ((atoms_not_hit || minimum_confidence_level < 100) && verbose) {
    cerr << atoms_not_hit << " of " << m.natoms() << " atoms not classified, confidence "
         << minimum_confidence_level << " (ok)\n";
  }

  return 0;  // the molecule is OK
}

/*
  Do any of the atoms in this embedding hit atoms already hit.
  We check the first ncheck atoms
*/

static int
all_atoms_unmatched(Set_of_Atoms& e, int ncheck, const int* already_hit)
{
  if (already_hit[e[0]]) {  // should knock out a lot of cases
    return 0;
  }

  if (e.number_elements() > ncheck) {
    e.resize(ncheck);
  } else {
    ncheck = e.number_elements();
  }

  if (1 == ncheck) {  // the most common case - we checked the first matched atom above
    return 1;
  }

  // check the rest of the atoms

  for (int i = 1; i < ncheck; i++) {
    if (already_hit[e[i]]) {
      return 0;
    }
  }

  return 1;
}

static int
ghose_crippin(Molecule& m, const IWString& unique_smiles, int* already_hit,
              std::ostream& output)
{
  Molecule_to_Match target(&m);

  int nq = queries.number_elements();

  double sum = 0.0;

  int minimum_confidence_level = 100;  // min confidence level encountered this molecule

  IWHistogram histogram;
  if (do_histogram) {
    histogram.initialise(histogram_start, histogram_stop, histogram_dx);
    histogram.set_put_out_of_range_values_in_first_or_last_bucket(1);
  }

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (verbose > 1 && nhits) {
      cerr << " query " << i << " '" << queries[i]->comment() << "' " << nhits
           << " hits to query " << i << " '" << queries[i]->comment() << "'\n";
    }

    if (0 == nhits) {
      continue;
    }

    double nvi;
    (void)queries[i]->numeric_value(nvi);

    for (int j = 0; j < nhits; j++) {
      Set_of_Atoms* e =
          const_cast<Set_of_Atoms*>(sresults.embedding(j));  // loss of const OK

      if (verbose > 1) {
        cerr << " qry " << i << " hit " << j << ' ' << *e << '\n';
      }

      if (all_atoms_unmatched(*e, atoms_to_match[i],
                              already_hit))  // may shorten E because of NMATCH value
      {
        sum += nvi;
        if (nullptr != frequency_vector) {
          frequency_vector[i]++;
        }
        if (mark_atoms[i]) {
          e->set_vector(already_hit, i + 1);
        }

#ifdef DEBUG_ATOMS_HIT
        for (auto x = 0; x < m.natoms(); ++x) {
          cerr << '(' << x << ") " << already_hit[x] << ' ';
        }
        cerr << '\n';
#endif

        if (confidence[i] < minimum_confidence_level) {
          minimum_confidence_level = confidence[i];
        }

        if (histogram.active()) {
          histogram.extra(nvi);
        }
      }
    }
  }

  min_confidence_found[minimum_confidence_level]++;

  sum = sum * multiplier + intercept;

  int atoms_not_hit = count_occurrences_of_item_in_array(0, target.natoms(), already_hit);

  // Write out labeled atoms now even if there is an error - useful for debugging

  if (stream_for_labeled_atoms.active() ||
      write_labelled_molecule_to_fingerprint.length()) {
    apply_isotopic_labels(m, already_hit);
    if (stream_for_labeled_atoms.active()) {
      stream_for_labeled_atoms.write(m);
    }
  }

  if (too_many_unclassified_atoms_or_uncertainty(m, already_hit, atoms_not_hit,
                                                 minimum_confidence_level)) {
    return 0;
  }

  if (bit_for_each_query_tag.length() || atom_label_tag.length() ||
      histogram_tag.length()) {
    write_fingerprint(m, unique_smiles, already_hit, histogram, output);
  }

  if (do_histogram) {
    write_histogram(m.name(), output, sum, histogram);
  } else if (write_as_array) {
    write_the_identifier(m.name(), output);

    output << ' ' << static_cast<float>(sum);

    if (include_atoms_not_hit_and_confidence_level_in_descriptors) {
      output << ' ' << atoms_not_hit << ' ' << minimum_confidence_level;
    }

    if (include_frequency_vector_with_output) {
      write_frequency_vector(frequency_vector, nq, output);
    }

    output << '\n';
  } else if (tag.length()) {
    output << tag << '<' << float(sum) << ">\n";
  } else {
    output << "Name " << m.name() << '\n';
    output << descriptor_name << ' ' << float(sum);
    if (append_unclassified_count) {
      output << ":UNCLASSIFIED=" << atoms_not_hit;
    }
    output << '\n';
    output << "|\n";
  }

  if (nullptr != frequency_vector && stream_for_frequency_vector.rdbuf()->is_open()) {
    write_space_suppressed_string(m.name(), stream_for_frequency_vector);
    write_frequency_vector(frequency_vector, nq, stream_for_frequency_vector);
    stream_for_frequency_vector << '\n';
  }

  return output.good();
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

  // cerr << m.smiles() << '\n';

  return;
}

static int
ghose_crippin(Molecule& m, std::ostream& output)
{
  preprocess(m);

  IWString unique_smiles;

  if (do_histogram) {  // we don't need the unique smiles
    unique_smiles = m.smiles();
  } else if (bit_for_each_query_tag.length() || atom_label_tag.length() ||
             histogram_tag.length()) {
    unique_smiles = m.smiles();
  }

  if (charge_assigner.active()) {
    charge_assigner.process(m);
  }

  if (make_implicit_hydrogens_implicit) {
    m.make_implicit_hydrogens_explicit();
  }

  int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "Ignoring empty molecule\n";
    return 1;
  }

  if (nullptr != frequency_vector) {
    set_vector(frequency_vector, queries.number_elements(), 0);
  }

  int* tmp = new_int(matoms);
  std::unique_ptr<int[]> free_tmp(tmp);

  return ghose_crippin(m, unique_smiles, tmp, output);
}

static int
ghose_crippin(data_source_and_type<Molecule>& input, std::ostream& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule()) && output.good()) {
    std::unique_ptr<Molecule> free_m(m);

    if (!ghose_crippin(*m, output)) {
      return 0;
    }
  }

  return output.good();
}

static int
ghose_crippin_tdt(const_IWSubstring& buffer, std::ostream& output)
{
  buffer.chop();  // remove trailing '>'

  buffer.remove_leading_chars(5);  // trim '$SMI<'

  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "Cannot parse smiles '" << buffer << "'\n";
    return 0;
  }

  return ghose_crippin(m, output);
}

static int
ghose_crippin_tdt(iwstring_data_source& input, std::ostream& output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer) && output.good()) {
    if (buffer.starts_with("$SMI<")) {
      ;
    } else {
      output << buffer << '\n';
      continue;
    }

    if (!ghose_crippin_tdt(buffer, output)) {  // fatal error
      return 0;
    }
  }

  return output.good();
}

static int
ghose_crippin_tdt(const char* fname, std::ostream& output)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ghose_crippin_tdt(input, output);
}

static int
ghose_crippin(const char* fname, FileType input_type, std::ostream& output)
{
  if (function_as_filter) {
    return ghose_crippin_tdt(fname, output);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ghose_crippin(input, output);
}

/*
  Extract a token like PREFIX=nnn
  from buffer
  We return 0 if an error is found
*/

static int
extract_qualifier(const const_IWSubstring& buffer, const char* prefix, int& result)
{
  int i = buffer.find(prefix);

  if (i < 0) {
    return 1;  // no value specified is OK
  }

  const_IWSubstring t = buffer.substr(i + strlen(prefix));

  if (t.starts_with('=')) {
    t++;
  }

  if (0 == t.length()) {
    cerr << "Cannot extract '" << prefix << "' value from '" << buffer << "'\n";
    return 0;
  }

  if (t.contains(' ')) {
    t.truncate_at_first(' ');
  }

  if (!t.numeric_value(result)) {
    cerr << "Invalid numeric for '" << prefix << "' value in '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
parse_control_file_record(const const_IWSubstring& buffer)
{
  if (buffer.nwords() < 3) {
    cerr << "Control file records must have at least 3 words\n";
    return 0;
  }

  const_IWSubstring fname_smarts, nv, comments;

  buffer.word(0, fname_smarts);

  buffer.word(1, nv);
  double nvv;
  if (!nv.numeric_value(nvv)) {
    cerr << "Bad numeric value '" << nv << "'\n";
    return 0;
  }

  comments = buffer;
  comments.remove_leading_words(2);

  Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;

  if (fname_smarts.starts_with("f:"))  // is a file name
  {
    if (!q->read(fname_smarts)) {
      cerr << "Cannot read query from '" << fname_smarts << "'\n";
      return 0;
    }
  } else {
    if (!q->create_from_smarts(fname_smarts)) {
      cerr << "Cannot read smarts '" << fname_smarts << "'\n";
      return 0;
    }
  }

  q->add_numeric_value(nvv);

  q->set_comment(comments);

  queries.add(q);

  if (verbose > 1) {
    cerr << "Created query '" << comments << "'\n";
  }

  if (comments.find("NOMARK") >= 0) {
    mark_atoms.add(0);  // don't mark the matched atoms (allows multiple matches);
  } else {
    mark_atoms.add(1);  // by default we mark all matched atoms
  }

  int iso = -1;

  if (!extract_qualifier(buffer, "ISO=", iso)) {  /// error
    return 0;
  }

  if (iso > 0) {
    label.add(iso);
  } else {
    label.add(queries.number_elements());
  }

  if (label.last_item() > highest_atom_label) {
    highest_atom_label = label.last_item();
  }

  int conf = -5;

  if (!extract_qualifier(buffer, "CONF=", conf)) {  /// error
    return 0;
  }

  if (conf >= 0) {
    confidence.add(conf);
  } else {
    confidence.add(100);
  }

  int nmatch = max_atoms_per_embedding_to_process;

  if (!extract_qualifier(buffer, "NMATCH=", nmatch)) {  /// error
    return 0;
  }

  atoms_to_match.add(nmatch);

  return 1;
}

static int
read_control_file(iwstring_data_source& input)
{
  int nr = input.records_remaining();

  queries.resize(queries.elements_allocated() + nr);
  label.resize(queries.elements_allocated() + nr);
  confidence.resize(queries.elements_allocated() + nr);
  mark_atoms.resize(queries.elements_allocated() + nr);
  atoms_to_match.resize(queries.elements_allocated() + nr);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    buffer.strip_trailing_blanks();
    buffer.strip_leading_blanks();
    if (buffer.starts_with('#')) {
      continue;
    }

    if (0 == buffer.length()) {
      continue;
    }

    if (verbose > 1) {
      cerr << "Examining " << input.lines_read() << " '" << buffer << "'\n";
    }

    if (!parse_control_file_record(buffer)) {
      cerr << "Cannot parse control file record '" << buffer << "'\n";
      return 0;
    }
  }

  assert(queries.number_elements() == mark_atoms.number_elements());
  assert(queries.number_elements() == label.number_elements());

  return queries.number_elements();
}

static int
read_control_file(IWString& fname)
{
  iwstring_data_source input(fname.c_str());

  if (!input.ok()) {
    cerr << "Cannot open control file '" << fname << "'\n";
    return 0;
  }

  return read_control_file(input);
}

static void
parse_dash_G_option(Command_Line& cl,
                    const Accumulator<double>& query_contribution_values,
                    std::ostream& output)
{
  const_IWSubstring g;
  int i = 0;

  int nbuckets = 0;

  while (cl.value('G', g, i++)) {
    if (g.starts_with("fptag=")) {
      bit_for_each_query_tag = g;
      bit_for_each_query_tag.remove_leading_chars(6);

      if (verbose) {
        cerr << "Atom type fingerprints written with tag '" << bit_for_each_query_tag
             << "'\n";
      }
    } else if (g.starts_with("altag=")) {
      atom_label_tag = g;
      g.remove_leading_chars(6);

      if (verbose) {
        cerr << "Atom label fingerprints written with tag '" << atom_label_tag << "'\n";
      }
    } else if (g.starts_with("fpmtag=")) {
      write_labelled_molecule_to_fingerprint = g;
      write_labelled_molecule_to_fingerprint.remove_leading_chars(7);
      write_labelled_molecule_to_fingerprint += '<';
    } else if (g.starts_with("buckets=")) {
      g.remove_leading_chars(8);
      if (!g.numeric_value(nbuckets) || nbuckets < 2) {
        cerr << "The buckets= histogram specifier must be a whole positive number\n";
        usage(4);
      }

      if (verbose) {
        cerr << "Histogram initialised with " << nbuckets << " buckets\n";
      }

      do_histogram = 1;

      histogram_start = query_contribution_values.minval();
      histogram_stop = query_contribution_values.maxval();
      histogram_dx =
          (histogram_stop - histogram_start) / static_cast<double>(nbuckets - 1);

      IWHistogram tmp;
      tmp.initialise(histogram_start, histogram_stop, histogram_dx);

      if (verbose) {
        cerr << "Histogram from " << tmp.minval() << " to " << tmp.maxval()
             << " dx = " << tmp.delta() << ", " << tmp.nbuckets() << " buckets\n";
      }
    } else if (g.starts_with("hist=")) {
      IWHistogram tmp;

      g.remove_leading_chars(5);

      if (!tmp.initialise(g)) {
        cerr << "Invalid histogram initialiser 'hist=" << g << "'\n";
        usage(11);
      }

      do_histogram = 1;

      histogram_start = tmp.minval();
      histogram_stop = tmp.maxval();
      histogram_dx = tmp.delta();

      if (verbose) {
        cerr << "Histogram from " << histogram_start << " to " << histogram_stop
             << " dx = " << histogram_dx << ", " << tmp.nbuckets() << " buckets\n";
      }

      nbuckets = tmp.nbuckets();
    } else if ("help" == g) {
      display_histogram_options(cerr);
      exit(0);
    } else {
      cerr << "Unrecognised histogram specifier '" << g << "'\n";
      usage(8);
    }
  }

  if (histogram_tag.length() && 0 == nbuckets) {
    cerr << "FPTAG specified, but no histogram or buckets\n";
    usage(5);
  }

  if (do_histogram) {
    digits.set_include_leading_space(1);
    digits.initialise(100);
    if (0 == histogram_tag.length()) {
      write_histogram_header(output, descriptor_name, nbuckets);
    }
  }

  if (write_to_stdout) {
    write_to_stdout = 0;
  }

  return;
}

static int
write_header_records(std::ostream& os)
{
  IWString buffer;

  if (experimental_column >= 0) {
    buffer << "expt ";
  }
  buffer << descriptor_name;

  if (0 == append_unclassified_count) {
    ;
  } else if (max_unclassified_atoms >= 0 || minimum_confidence_level_needed) {
    buffer << " unclassified confidence";
  }

  if (include_frequency_vector_with_output) {
    int nq = queries.number_elements();

    for (int i = 0; i < nq; i++) {
      if (!use_query_name_as_descriptor_name) {
        buffer << " ghose_GCAA" << i;
      } else if (queries[i]->comment().contains(' ')) {
        IWString tmp(queries[i]->comment());
        tmp.gsub(' ', '_');
        buffer.append_with_spacer(tmp);
      } else {
        buffer.append_with_spacer(queries[i]->comment());
      }
    }
  }

  if (write_to_stdout) {
    os << "Name " << buffer << '\n';
  }

  return os.good();
}

/*
 */

static int
ghose_crippin(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:o:g:F:R:D:u:L:lm:S:ayN:Hhc:T:d:nI:M:x:G:fJ:rU");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      return 17;
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot process element specification option(s), -E option\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      return 62;
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will reduce multi-fragment molecules to largest fragment\n";
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot initialise charge assigner\n";
      usage(3);
    }
  }

  if (cl.option_present('J')) {
    int i = 0;
    const_IWSubstring j;
    while (cl.value('J', j, i++)) {
      if ("blki" == j) {
        set_aromatic_bonds_lose_kekule_identity(1);
        if (verbose) {
          cerr << "Will use strict Daylight rules for aliphatic atoms\n";
        }
      } else if ("qndn" == j) {
        use_query_name_as_descriptor_name = 1;
        if (verbose) {
          cerr << "Will use query names as descriptor names\n";
        }
      } else if ("help" == j) {
        display_dash_j_qualifiers(cerr);
        return 0;
      } else {
        cerr << "Unrecognised -J qualifier '" << j << "'\n";
        return 5;
      }
    }
  }

  if (cl.option_present('r')) {
    set_aromatic_bonds_lose_kekule_identity(1);
    if (verbose) {
      cerr << "Will use strict Daylight rules for aliphatic atoms - use '-J blki' "
              "instead\n";
    }
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', max_unclassified_atoms)) {
      cerr << "The -u option must be followed by a whole number\n";
      usage(18);
    }

    if (!verbose) {
      ;
    } else if (max_unclassified_atoms > 0) {
      cerr << "At most " << max_unclassified_atoms << " unclassified atoms are allowed\n";
    } else {
      cerr << "Any number of unclassified atoms possible\n";
    }
  }

  if (cl.option_present('D')) {
    cl.value('D', descriptor_name);
    if (verbose) {
      cerr << "Values for '" << descriptor_name << "' will be written\n";
    }
  }

  if (cl.option_present('f')) {
    function_as_filter = 1;

    if (verbose) {
      cerr << "Will function as a filter\n";
    }
  }

  if (cl.option_present('h')) {
    make_implicit_hydrogens_implicit = 1;
    if (verbose) {
      cerr << "Implicit hydrogens will be made explicit\n";
    }
  }

  if (cl.option_present('I')) {
    if (!cl.value('I', intercept)) {
      cerr << "The intercept (-I option) must be a valid number\n";
      usage(6);
    }

    if (verbose) {
      cerr << "Intercept " << intercept << '\n';
    }
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', multiplier)) {
      cerr << "The multiplier (-M option) must be a valid number\n";
      usage(6);
    }

    if (verbose) {
      cerr << "Multiplier " << multiplier << '\n';
    }
  }

  // the -m option must be done before the -F option

  if (cl.option_present('m')) {
    const_IWSubstring m;
    cl.value('m', m);

    if ("all" == m) {
      max_atoms_per_embedding_to_process = std::numeric_limits<int>::max();
      if (verbose) {
        cerr << "Will process all matched atoms - fragment additivity method\n";
      }
    } else if (!m.numeric_value(max_atoms_per_embedding_to_process) ||
               max_atoms_per_embedding_to_process < 1) {
      cerr << "The -m option must be followed by a whole positive number\n";
      usage(38);
    } else if (verbose) {
      cerr << "A max of " << max_atoms_per_embedding_to_process
           << " atoms per query match will be processed\n";
    }
  }

  if (cl.option_present('x')) {
    if (!cl.value('x', experimental_column) || experimental_column < 1) {
      cerr << "The experimental column (-x option) must be a valid column number\n";
      usage(18);
    }

    if (verbose) {
      cerr << "The experimental values are the " << experimental_column
           << " token in the name field\n";
    }

    experimental_column--;
  }

  if (!cl.option_present('F')) {
    cerr << "Must specify one or more control files via the -F option\n";
    usage(19);
  }

  if (cl.option_present('F')) {
    IWString f;
    int i = 0;
    while (cl.value('F', f, i++)) {
      if (!read_control_file(f)) {
        cerr << "Cannot read control file '" << f << "'\n";
        return 73;
      }
    }

    if (verbose) {
      cerr << "Read " << queries.number_elements() << " queries\n";
    }
  }

  int unique_embeddings_only = 0;

  if (cl.option_present('U')) {
    unique_embeddings_only = 1;

    if (verbose) {
      cerr << "Will find only unique embeddings\n";
    }
  }

  int nq = queries.number_elements();
  assert(nq > 0);

  Accumulator<double> query_contribution_values;

  for (int i = 0; i < nq; i++) {
    queries[i]->set_find_one_embedding_per_atom(1);

    if (unique_embeddings_only) {
      queries[i]->set_find_unique_embeddings_only(
          1);  // do not want this all the time. Some queries are symmetric and since we
               // generally only process the first atom, we want all matches
    }

    double v;
    if (!queries[i]->numeric_value(v)) {
      cerr << "Yipes, no numeric value for query " << i << " '" << queries[i]->comment()
           << "'\n";
      return i + 1;
    }

    query_contribution_values.extra(v);
  }

  if (verbose) {
    cerr << "Atomic contributions between " << query_contribution_values.minval()
         << " and " << query_contribution_values.maxval() << '\n';
  }

  if (cl.option_present('G')) {
    parse_dash_G_option(cl, query_contribution_values, std::cout);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('R')) {
    IWString r;
    cl.value('R', r);

    stream_for_frequency_vector.open(r.null_terminated_chars(), std::ios::out);

    if (!stream_for_frequency_vector.good()) {
      cerr << "Cannot open frequency vector file '" << r << "'\n";
      return 23;
    }

    if (verbose) {
      cerr << "Frequency vector info written to '" << r << "'\n";
    }

    stream_for_frequency_vector << "Name";
    for (int i = 0; i < nq; i++) {
      if (use_query_name_as_descriptor_name) {
        write_space_suppressed_string(queries[i]->comment(), stream_for_frequency_vector);
      } else {
        stream_for_frequency_vector << ' ' << "GCAA" << i;
      }
    }
    stream_for_frequency_vector << '\n';

    frequency_vector = new int[queries.number_elements()];
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', minimum_confidence_level_needed) ||
        minimum_confidence_level_needed < 0 || minimum_confidence_level_needed > 100) {
      cerr << "The minimum confidence level option (-c) must be followed by a value "
              "confidence level (0, 100)\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Molecules with confidence levels below " << minimum_confidence_level_needed
           << " will be flagged as erroneous\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('L')) {
    if (!cl.option_present('o')) {
      stream_for_labeled_atoms.add_output_type(FILE_TYPE_SMI);
    } else if (!stream_for_labeled_atoms.determine_output_types(cl)) {
      cerr << "Cannot determine output types for labeled atom molecules\n";
      return 0;
    }

    const_IWSubstring l;
    cl.value('L', l);

    if (!stream_for_labeled_atoms.new_stem(l, 1)) {
      cerr << "Cannot open -L file(s) with stem '" << l << "'\n";
      return 53;
    }

    if (verbose) {
      cerr << "Molecules with labeled atoms written to '" << l << "'\n";
    }
  }

  // If the -S option is present, we respect that. Otherwise, use cout.
  // We need a pointer to the appropriate std::ostream

  std::ostream* destination = &std::cout;

  std::ofstream stream_for_output;
  if (cl.option_present('S')) {
    IWString stem_for_output;
    cl.value('S', stem_for_output);

    stream_for_output.open(stem_for_output.null_terminated_chars());
    if (!stream_for_output.good()) {
      cerr << "Cannot open '" << stem_for_output << "'\n";
      return 87;
    }

    if (verbose) {
      cerr << "Results written to '" << stem_for_output << "'\n";
    }

    destination = &stream_for_output;
  }

  if (cl.option_present('T') && cl.option_present('a')) {
    cerr << "The -a (write as array) and -T (work as TDT filter) options are mutually "
            "exclusive\n";
    usage(21);
  }

  if (cl.option_present('a')) {
    write_as_array = 1;

    if (verbose) {
      cerr << "Results output in array form\n";
    }

    if (cl.option_count('a') > 1) {
      include_frequency_vector_with_output = 1;
      if (verbose) {
        cerr << "Will include the frequency vector with the output\n";
      }

      if (nullptr == frequency_vector) {
        frequency_vector = new int[queries.number_elements()];
      }

      if (0 == digits.number_elements()) {
        digits.set_include_leading_space(1);
        digits.initialise(100);
      }
    }

    if (!write_header_records(*destination)) {
      return 8;
    }
  }

  if (cl.option_present('T')) {
    cl.value('T', tag);
    if (verbose) {
      cerr << "Results written to dataitem '" << tag << "'\n";
    }
  }

  if (cl.option_present('y')) {
    if (write_as_array) {
      cerr << "Sorry, the -a and -y options are mutually exclusive\n";
      usage(37);
    }

    append_unclassified_count = 1;
    if (verbose) {
      cerr << "Will append unclassified atom count to output\n";
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!ghose_crippin(cl[i], input_type, *destination)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Details on matches for " << nq << " queries\n";
    for (int i = 0; i < nq; i++) {
      queries[i]->report(cerr, verbose);
    }

    for (int i = 0; i < min_confidence_found.number_elements(); i++) {
      if (min_confidence_found[i]) {
        cerr << min_confidence_found[i] << " molecules had a confidence level of " << i
             << '\n';
      }
    }

    for (int i = 0; i < unclassified_atom_count.number_elements(); i++) {
      if (unclassified_atom_count[i]) {
        cerr << unclassified_atom_count[i] << " molecules had " << i
             << " unclassified atoms\n";
      }
    }
  }

  if (nullptr != frequency_vector) {
    delete[] frequency_vector;
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = ghose_crippin(argc, argv);

  return rc;
}
