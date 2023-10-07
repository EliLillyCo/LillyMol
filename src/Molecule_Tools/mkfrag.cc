/*
  Write all fragments in a molecule as individual molecules
*/

#include <iostream>
#include <limits>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/ematch.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

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
  cerr << "Writes a multi-fragment molecule as individual molecules\n";
  cerr << "  -f <natoms>    don't write fragments with <natoms> or fewer atoms\n";
  cerr << "  -F <natoms>    don't write fragments with <natoms> or more  atoms\n";
  cerr << "  -M <nfrag>     don't write molecules with <nfrag>  or more  fragments\n";
  cerr << "  -m <nfrag>     don't write molecules with <nfrag>  or fewer fragments\n";
  cerr << "  -O ''          \"organic\" fragments only\n";
  cerr << "  -O el          \"organic\" fragments only, but 'el' is OK (repeat for each OK ele)\n";
  cerr << "  -V             discard fragments containing invalid valences\n";
  cerr << "  -u             unique fragments only\n";
  cerr << "  -s             sort fragments by size\n";
  cerr << "  -a             append the fragment id to each component\n";
  cerr << "  -l             fragment id's are local to each molecule - will get duplicate numbers\n";
  cerr << "  -D isotope     ignore isotopic labels when determining duplicates\n";
  cerr << "  -D chiral      ignore chirality when determining duplicates\n";
  cerr << "  -I transform   convert all molecules to non isotopic form on input\n";
  cerr << "  -I discard     discard fragments containing isotopic atoms\n";   
  cerr << "  -w             write the whole molecule before fragmenting\n";
  cerr << "  -U <file>      write unique smiles of all fragments encountered to <file>\n";
  cerr << "  -R <sep>       reverse functionality - combine all input molecules into 1. <sep> is name field separator\n";
  cerr << "  -c             remove chirality\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -i <type>      input type\n";
  cerr << "  -o <type>      output type\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int verbose = 0;

static int molecules_read = 0;

static extending_resizable_array<int> frag_count;
static extending_resizable_array<int> frag_count_written;

static extending_resizable_array<int> atoms_in_fragment;

static extending_resizable_array<int> largest_fragment;
static extending_resizable_array<int> smaller_fragments;

static int min_frag_count = 0;

static int max_frag_count = 0;

static int min_frag_size_to_write = 0;
static int max_frag_size_to_write = std::numeric_limits<int>::max();

static int unique_fragments_only = 0;

static int append_frag_id = 0;

static int write_whole_molecule = 0;

static int sort_by_size = 0;

static int compare_with_chirality = 1;
static int compare_with_isotopes = 1;
static int transform_to_non_isotopic_form = 0;
static int discard_isotopic_fragments = 0;
static int discard_fragments_with_invalid_valences = 0;

/*
  We can optionally collect statistics on all fragments encountered
*/

static IW_STL_Hash_Map_int global_frag_count;

static Chemical_Standardisation chemical_standardisation;

static Set_of_Element_Matches ok_non_organics;

static int organic_fragments_only = 0;

static std::ofstream ufrag_stream;

/*
  Dec 2002. Qi noticed that numbers are within each fragment, and therefore repeated
  when multiple molecules are input.
*/

static int fragment_ids_are_global = 1;

static int molecules_written = 0;

static int reverse_mkfrag = 0;  // combine all input molecules into one

static IWString reverse_mkfrag_name_separator;

static int remove_chirality = 0;

static int
molecule_size_comparitor(Molecule* const* ppm1, Molecule* const* ppm2)
{
  Molecule* pm1 = const_cast<Molecule*>(*ppm1);
  Molecule* pm2 = const_cast<Molecule*>(*ppm2);

  int natoms1 = pm1->natoms();
  int natoms2 = pm2->natoms();

  if (natoms1 > natoms2) {
    return -1;
  }
  if (natoms1 < natoms2) {
    return 1;
  }

  molecular_weight_t amw1 = pm1->molecular_weight();
  molecular_weight_t amw2 = pm2->molecular_weight();

  if (amw1 > amw2) {
    return -1;
  }
  if (amw1 < amw2) {
    return 1;
  }

  return 0;
}

static void
remove_fragments_for_size(resizable_array_p<Molecule>& components)
{
  int nfrag = components.number_elements();
  for (int i = 0; i < nfrag; i++) {
    int n = components[i]->natoms();
    if (n <= min_frag_size_to_write || n >= max_frag_size_to_write) {
      components.remove_item(i);
      nfrag--;
      i--;
    }
  }

  return;
}

static int
contains_non_organic_atoms(const Molecule& m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Element* e = m.elementi(i);

    if (e->organic()) {
      continue;
    }

    if (!ok_non_organics.matches(e)) {
      return 1;
    }
  }

  return 0;
}

static void
remove_non_organic_fragments(resizable_array_p<Molecule>& components)
{
  int nc = components.number_elements();

  for (int i = 0; i < nc; i++) {
    const Molecule* m = components[i];

    if (contains_non_organic_atoms(*m)) {
      components.remove_item(i);
      nc--;
      i--;
    }
  }

  return;
}

static void
remove_isotopic_fragments(resizable_array_p<Molecule>& components)
{
  int nc = components.number_elements();

  for (int i = 0; i < nc; i++) {
    const Molecule* m = components[i];

    if (m->number_isotopic_atoms()) {
      components.remove_item(i);
      nc--;
      i--;
    }
  }

  return;
}

static void
do_discard_fragments_with_invalid_valences(resizable_array_p<Molecule>& components)
{
  int nc = components.number_elements();

  for (int i = 0; i < nc; i++) {
    Molecule* m = const_cast<Molecule*>(components[i]);

    if (!m->valence_ok()) {
      components.remove_item(i);
      nc--;
      i--;
    }
  }

  return;
}

static int
do_write_whole_molecule(Molecule& m, int nfrag, Molecule_Output_Object& output)
{
  IWString tmp;
  tmp.resize(m.name().length() + 10);

  tmp = m.name();

  tmp += " NFRAG=";
  tmp += nfrag;

  m.set_name(tmp);

  return output.write(m);
}

static int
form_unique_fragments(resizable_array_p<Molecule>& components, int* xref)
{
  int nfrag = components.number_elements();

  IW_STL_Hash_Map_int usmi;

  if (!compare_with_chirality) {
    set_include_chiral_info_in_smiles(0);
  }

  if (!compare_with_isotopes) {
    set_include_isotopes_in_smiles(0);
  }

  int duplicates_found = 0;
  for (int i = 0; i < nfrag; i++) {
    Molecule* m = components[i];

    const IWString& u = m->unique_smiles();  // careful, reference to string, don't
                                             // invalidate the smiles too soon

    int j;
    if (usmi.contains(u)) {
      j = usmi[u];
      duplicates_found++;
    } else {
      usmi[u] = i;
      j = i;
    }

    //  cerr << "xref for '" << u << "' is " << j << endl;

    m->invalidate_smiles();  // do this after all uses of 'U' are complete

    xref[i] = j;
  }

  set_include_chiral_info_in_smiles(1);

  set_include_isotopes_in_smiles(1);

  if (0 == duplicates_found) {
    return nfrag;
  }

  // Now count the number of occurrences of each fragment

  for (int i = 0; i < nfrag; i++) {
    if (xref[i] != i) {  // is a duplicate of something earlier in the list
      continue;
    }

    int count = 1;
    for (int j = i + 1; j < nfrag; j++) {
      if (xref[j] == i) {
        count++;
      }
    }

    IWString tmp = components[i]->name();
    tmp += "COUNT=";
    tmp += count;

    components[i]->set_name(tmp);
    components[i]->invalidate_smiles();
  }

  // Now get rid of the duplciates, anything whose fragment number does not match its
  // index

  for (int i = nfrag - 1; i > 0; i--) {
    if (xref[i] != i) {
      components.remove_item(i);
      nfrag--;
    }
  }

  return nfrag;
}

static int
form_unique_fragments(resizable_array_p<Molecule>& components)
{
  int nfrag = components.number_elements();

  int* tmp = new int[nfrag];
  std::unique_ptr<int[]> free_tmp(tmp);

  return form_unique_fragments(components, tmp);
}

/*
  We are writing fragment FRAG_NUMBER of a molecule.
*/

static int
do_output(Molecule& m, const IWString& mname, int frag_number,
          Molecule_Output_Object& output)
{
  if (append_frag_id) {
    if (fragment_ids_are_global) {
      frag_number = molecules_written;
    }

    IWString tmp = mname;
    tmp += " FRAG ";
    tmp += frag_number;
    m.set_name(tmp);
  } else {
    m.set_name(mname);
  }

  molecules_written++;

  return output.write(m);
}

static int
mkfrag(Molecule& m, int nfrag, Molecule_Output_Object& output)
{
  if (write_whole_molecule && !do_write_whole_molecule(m, nfrag, output)) {
    return 0;
  }

  if (1 == nfrag) {
    if (organic_fragments_only && contains_non_organic_atoms(m)) {
      return 1;
    }

    if (discard_fragments_with_invalid_valences && !m.valence_ok()) {
      return 1;
    }

    return do_output(m, m.name(), 0, output);
  }

  resizable_array_p<Molecule> components;
  nfrag = m.create_components(components);

  if (verbose) {
    int largest_fragment_this_molecule = 0;
    for (int i = 0; i < nfrag; i++) {
      int a = components[i]->natoms();

      atoms_in_fragment[a]++;

      if (a > largest_fragment_this_molecule) {
        largest_fragment_this_molecule = a;
      }
    }

    largest_fragment[largest_fragment_this_molecule]++;

    for (int i = 0; i < nfrag; i++) {
      int a = components[i]->natoms();

      if (a < largest_fragment_this_molecule) {
        smaller_fragments[a]++;
      }
    }
  }

  if (min_frag_size_to_write || max_frag_size_to_write) {
    remove_fragments_for_size(components);
    nfrag = components.number_elements();
  }

  if (organic_fragments_only) {
    remove_non_organic_fragments(components);
    nfrag = components.number_elements();
  }

  if (discard_isotopic_fragments && m.number_isotopic_atoms()) {
    remove_isotopic_fragments(components);
    nfrag = components.number_elements();
  }

  if (discard_fragments_with_invalid_valences) {
    do_discard_fragments_with_invalid_valences(components);
    nfrag = components.number_elements();
  }

  if (nfrag > 1 && unique_fragments_only) {
    nfrag = form_unique_fragments(components);
  }

  if (nfrag > 1 && sort_by_size) {
    components.sort(molecule_size_comparitor);
  }

  for (int i = 0; i < nfrag; i++) {
    Molecule* f = components[i];

    if (!do_output(*f, m.name(), i, output)) {
      return 0;
    }

    if (ufrag_stream.rdbuf()->is_open()) {
      global_frag_count[f->unique_smiles()]++;
    }
  }

  frag_count_written[nfrag]++;

  return output.good();
}

static void
preprocess(Molecule& m)
{
  if (transform_to_non_isotopic_form) {
    m.transform_to_non_isotopic_form();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (remove_chirality) {
    m.remove_all_chiral_centres();
  }

  return;
}

static int
mkfrag(data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    //  cerr << "Read " << m->smiles() << "'\n";

    int nfrag = m->number_fragments();

    if (verbose > 1) {
      cerr << molecules_read << " '" << m->name() << "' has " << nfrag << " fragments\n";
    }

    frag_count[nfrag]++;

    if (min_frag_count > 0 && nfrag <= min_frag_count) {
      ;
    } else if (max_frag_count > 0 && nfrag >= max_frag_count) {
      ;
    } else if (!mkfrag(*m, nfrag, output)) {
      return 0;
    }
  }

  return 1;
}

static int
mkfrag(const char* fname, FileType input_type, Molecule_Output_Object& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "mkfrag:cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return mkfrag(input, output);
}

static int
do_name_appending(Molecule& m, const IWString& s)
{
  IWString tmp(m.name());
  tmp.append_with_spacer(s, reverse_mkfrag_name_separator);
  m.set_name(tmp);

  return 1;
}

static int
do_reverse_mkfrag(data_source_and_type<Molecule>& input, IW_STL_Hash_Set& already_seen,
                  Molecule& tot)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    if (m->natoms() < min_frag_size_to_write) {
      continue;
    }

    if (m->natoms() > max_frag_size_to_write) {
      continue;
    }

    if (unique_fragments_only) {
      if (already_seen.contains(m->unique_smiles())) {
        continue;
      }

      already_seen.emplace(m->unique_smiles());
    }

    tot.add_molecule(m);
    do_name_appending(tot, m->name());
  }

  return 1;
}

static int
do_reverse_mkfrag(const char* fname, FileType input_type, IW_STL_Hash_Set& already_seen,
                  Molecule& tot)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "mkfrag:cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return do_reverse_mkfrag(input, already_seen, tot);
}

static int
do_reverse_mkfrag(Command_Line& cl, FileType input_type, Molecule_Output_Object& output)
{
  Molecule tot;
  IW_STL_Hash_Set already_seen;

  for (int i = 0; i < cl.number_elements(); ++i) {
    if (!do_reverse_mkfrag(cl[i], input_type, already_seen, tot)) {
      cerr << "Fatal erorr processing '" << cl[i] << "'\n";
      return 0;
    }
  }

  return output.write(tot);
}

static int
mkfrag(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:m:M:i:o:S:f:F:uasD:I:U:g:wO:VlR:c");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (!process_elements(cl, verbose)) {
    usage(1);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g) option\n";
      return 9;
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', min_frag_count) || min_frag_count <= 0) {
      cerr << "The minimum fragment count option (-m) must be followed by a whole "
              "positive number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Molecules with fewer than " << min_frag_count
           << " fragments will not be written\n";
    }
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', max_frag_count) || max_frag_count <= 0) {
      cerr << "The maximum fragment count option (-M) must be followed by a whole "
              "positive number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Molecules with more  than " << max_frag_count
           << " fragments will not be written\n";
    }
  }

  if (cl.option_present('f')) {
    if (!cl.value('f', min_frag_size_to_write) || min_frag_size_to_write < 1) {
      cerr << "The minimum fragment size to write option (-f) must be followed by a "
              "whole positive number\n";
      usage(8);
    }

    if (verbose) {
      cerr << "Fragments with fewer than " << min_frag_size_to_write
           << " will not be written\n";
    }
  }

  if (cl.option_present('F')) {
    if (!cl.value('F', max_frag_size_to_write) || max_frag_size_to_write < 1) {
      cerr << "The maximum fragment size to write option (-F) must be followed by a "
              "whole positive number\n";
      usage(8);
    }

    if (verbose) {
      cerr << "Fragments with more  than " << max_frag_size_to_write
           << " will not be written\n";
    }
  }

  if (cl.option_present('u')) {
    unique_fragments_only = 1;
    if (verbose) {
      cerr << "Will only write the unique fragments\n";
    }
  }

  if (cl.option_present('a')) {
    append_frag_id = 1;
    if (verbose) {
      cerr << "Will append the fragment ID to each fragment\n";
    }
  }

  if (cl.option_present('c')) {
    remove_chirality = 1;

    if (verbose) {
      cerr << "Will remove chirality\n";
    }
  }

  if (cl.option_present('s')) {
    sort_by_size = 1;
    if (verbose) {
      cerr << "Will sort fragments by size (natoms then amw)\n";
    }
  }

  if (cl.option_present('D')) {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('D', d, i++)) {
      if ("isotope" == d) {
        compare_with_isotopes = 0;

        if (verbose) {
          cerr << "Will compare fragments without considering isotopic labels\n";
        }
      } else if (d.starts_with("chiral")) {
        compare_with_chirality = 0;
        if (verbose) {
          cerr << "Will ignore chirality when determining duplicate fragments\n";
        }
      } else {
        cerr << "Unrecognised duplicate qualifier '" << d << "'\n";
        usage(14);
      }
    }
  }

  if (cl.option_present('I')) {
    const_IWSubstring i = cl.string_value('I');

    if ("transform" == i) {
      transform_to_non_isotopic_form = 1;

      if (verbose) {
        cerr << "All isotopes converted to natural form\n";
      }
    } else if ("discard" == i) {
      discard_isotopic_fragments = 1;

      if (verbose) {
        cerr << "Fragments with isotopes will be discarded\n";
      }
    } else {
      cerr << "Unrecognised isotope qualifier '" << i << "'\n";
      usage(6);
    }
  }

  if (cl.option_present('O')) {
    if (!ok_non_organics.construct_from_command_line(cl, verbose, 'O')) {
      cerr << "Could not ascertain list of OK non-organics from(s)\n";
      usage(9);
    }
    organic_fragments_only = 1;
    if (verbose) {
      cerr << "Only Organic fragments will be output\n";
    }
  }

  if (cl.option_present('V')) {
    discard_fragments_with_invalid_valences = 1;

    if (verbose) {
      cerr << "Will discard fragments with invalid valences\n";
    }
  }

  if (cl.option_present('w')) {
    write_whole_molecule = 1;

    if (verbose) {
      cerr << "Whole molecule written as well as fragments\n";
    }
  }

  if (cl.option_present('l')) {
    fragment_ids_are_global = 0;

    if (verbose) {
      cerr << "Fragment numbers are local within each molecule\n";
    }
  }

  if (cl.option_present('R')) {
    reverse_mkfrag = 1;
    reverse_mkfrag_name_separator = cl.string_value('R');
    if ("space" == reverse_mkfrag_name_separator) {
      reverse_mkfrag_name_separator = ' ';
    } else if ("tab" == reverse_mkfrag_name_separator) {
      reverse_mkfrag_name_separator = '\t';
    } else if ("empty" == reverse_mkfrag_name_separator) {
      ;
    }

    if (verbose) {
      cerr << "Beware, not all option supported by reverse\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      usage(1);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot recognise input files by name\n";
    return 3;
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 6;
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);
    if (verbose) {
      cerr << "Output type defaults to smi\n";
    }
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    usage(2);
  }

  const_IWSubstring s;
  if (cl.option_present('S')) {
    s = cl.string_value('S');
  } else {
    s = '-';
  }

  if (output.would_overwrite_input_files(cl, s)) {
    cerr << "Cannot overwrite our input, output stem = '" << s << "'\n";
    return 5;
  }

  if (!output.new_stem(s)) {
    cerr << "Cannot open output stream with stem '" << s << "'\n";
    return 14;
  }

  if (cl.option_present('U')) {
    const char* f = cl.option_value('U');

    ufrag_stream.open(f, std::ios::out);
    if (!ufrag_stream.good()) {
      cerr << "Cannot open unique fragment stream file '" << f << "'\n";
      return 19;
    }

    if (verbose) {
      cerr << "Smiles for all fragments encountered written to '" << f << "'\n";
    }
  }

  frag_count.resize(1000);  // let's hope this is enough
  frag_count_written.resize(1000);

  int rc = 0;
  if (reverse_mkfrag) {
    if (!do_reverse_mkfrag(cl, input_type, output)) {
      rc = 1;
    }

    if (verbose) {
      cerr << "Read " << molecules_read << " molecules\n";
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!mkfrag(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }

    if (verbose) {
      cerr << "Read " << molecules_read << " molecules\n";
      for (int i = 0; i < frag_count.number_elements(); i++) {
        if (frag_count[i]) {
          cerr << frag_count[i] << " molecules had " << i << " fragments\n";
        }
      }

      for (int i = 0; i < frag_count_written.number_elements(); i++) {
        if (frag_count_written[i]) {
          cerr << frag_count_written[i] << " molecules had " << i
               << " fragments written\n";
        }
      }

      cerr << "Wrote " << molecules_written << " molecules\n";

      for (int i = 1; i < atoms_in_fragment.number_elements(); i++) {
        if (atoms_in_fragment[i]) {
          cerr << atoms_in_fragment[i] << " fragments had " << i << " atoms\n";
        }
      }

      for (int i = 1; i < largest_fragment.number_elements(); i++) {
        if (largest_fragment[i]) {
          cerr << largest_fragment[i] << " molecules had " << i
               << " atoms in largest fragment\n";
        }
      }
    }
  }

  if (ufrag_stream.rdbuf()->is_open()) {
    IW_STL_Hash_Map_int::const_iterator i;
    for (i = global_frag_count.begin(); i != global_frag_count.end(); ++i) {
      const IWString& smi = (*i).first;
      int count = (*i).second;

      ufrag_stream << smi << " COUNT " << count << endl;
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = mkfrag(argc, argv);

  return rc;
}
