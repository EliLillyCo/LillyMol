/*
  Sometimes we need to determine the unique molecules from a file
  of structures.
*/

#include <stdlib.h>
#include <sys/stat.h>

#include <iostream>
#include <memory>

using std::cerr;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

static int verbose = 0;

const char* prog_name = nullptr;

static Elements_to_Remove elements_to_remove;

static Element_Transformations element_transformations;

static Number_Assigner number_assigner;

/*
  How should element transformations and removals be used?

  By default, we apply the removal and transformations and the resulting molecule
  is what gets written. But, what we probably want is for the transformed
  molecule to be used for comparisons only, and we want the original
  molecule written
*/

static int discard_molecule_changes = 0;

static Chemical_Standardisation chemical_standardisation;

static int molecules_read = 0;

static int duplicates_found = 0;

static int unique_molecules = 0;

static int reduce_to_largest_fragment = 0;

static Molecule_Output_Object unique_molecule_stream;

static Molecule_Output_Object duplicate_molecule_stream;

static int compare_as_graph = 0;

static int ignore_isotopes = 0;

// If there is an isotope, transform to a single value.
static int all_isotopes_become_identical = 0;

static constexpr isotope_t kFixedIsotope = 1;

static int default_primary_hash_size = 1000;

static int exclude_chiral_info = 0;

static int exclude_cis_trans_bonding_info = 0;

static int function_as_filter = 0;

static Report_Progress report_progress;

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

/*
  Jul 2005. Jibo wants to be able to enable the case where same
  structure but different name will be considered different
*/

static int only_same_if_structure_and_name_the_same = 0;

static int perform_formula_check = 0;  // slows things down!

static int use_atom_hash = 0;

/*
  Here's where we keep track of the previously encountered smiles
  Hashes get less and less efficient as they get larger, so have
  a different hash for every atom count
*/

static resizable_array_p<IW_STL_Hash_Set> formula_hash;
static resizable_array_p<IW_STL_Hash_Map_int> smiles_hash;
static resizable_array_p<std::unordered_set<uint64_t>> atom_hash;

static int store_dicer_fragments = 0;
static resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>> smiles_to_dicer_fragment;

template class resizable_array_p<IW_STL_Hash_Map_int>;
template class resizable_array_base<IW_STL_Hash_Map_int*>;
template class resizable_array_p<IW_STL_Hash_Set>;
template class resizable_array_base<IW_STL_Hash_Set*>;

static resizable_array_p<IWReaction> reaction;
static int number_reactions = 0;
static int molecules_changed_by_reactions = 0;

static int
UpdateDicerFragmentHash(Molecule& m, const IWString& usmi) {
  const int matoms = m.natoms();

  absl::flat_hash_map<IWString, dicer_data::DicerFragment>* smiles_hash = smiles_to_dicer_fragment[matoms];

  auto iter = smiles_hash->find(usmi);
  if (iter == smiles_hash->end()) {
    dicer_data::DicerFragment proto;
    proto.set_smi(m.smiles().data(), m.smiles().length());
    proto.set_par(m.name().data(), m.name().length());
    proto.set_n(1);
    proto.set_nat(m.natoms());
    smiles_hash->emplace(usmi, std::move(proto));
    ++unique_molecules;
    return 1;
  } else {
    uint32_t n = iter->second.n();
    iter->second.set_n(n + 1);
    ++duplicates_found;
    return 0;
  }
}

static int
UpdateSmilesHash(const Molecule& m, const IWString& usmi) {
  const int matoms = m.natoms();

  IW_STL_Hash_Map_int* h = smiles_hash[matoms];

  IW_STL_Hash_Map<IWString, int>::iterator f = h->find(usmi);

  if (f == h->end()) {
    (*h)[usmi] = 1;
    unique_molecules++;
    return 1;
  } else {
    duplicates_found++;
    (*f).second++;
    return 0;
  }
}

static int
UpdateHash(Molecule& m, const IWString& usmi) {
  if (store_dicer_fragments) {
    return UpdateDicerFragmentHash(m, usmi);
  } else {
    return UpdateSmilesHash(m, usmi);
  }
}

static int
is_unique_molecule(Molecule& m) {
  if (exclude_chiral_info) {
    set_include_chiral_info_in_smiles(0);
  }

  if (exclude_cis_trans_bonding_info) {
    set_include_cis_trans_in_smiles(0);
  }

  IWString usmi;

  if (compare_as_graph) {
    IWString formula = m.molecular_formula();

    Molecule g;
    g.add_molecule(&m);
    g.change_to_graph_form();

    usmi = g.unique_smiles();

    usmi << ':' << formula;

    //  cerr << "Graph " << usmi << '\n';
  } else {
    usmi = m.unique_smiles();
  }

  // cerr << "Testing unique smiles '" << usmi << "'\n";

  if (exclude_chiral_info) {  // reset back to useful value
    set_include_chiral_info_in_smiles(1);
  }

  if (exclude_cis_trans_bonding_info) {
    set_include_cis_trans_in_smiles(1);
  }

  if (only_same_if_structure_and_name_the_same) {
    usmi << m.name();
  }

  const int matoms = m.natoms();

  if (store_dicer_fragments) {
    while (matoms >= smiles_to_dicer_fragment.number_elements()) {
      smiles_to_dicer_fragment << new absl::flat_hash_map<IWString, dicer_data::DicerFragment>;
    }
  } else {
    while (matoms >= smiles_hash.number_elements()) {
      smiles_hash.add(new IW_STL_Hash_Map_int);
      if (use_atom_hash) {
        atom_hash.add(new std::unordered_set<uint64_t>());
      }
    }
  }

  // cerr << matoms << " smiles_hash " << smiles_hash.number_elements() << " formula_hash
  // " << formula_hash.number_elements() << '\n';

  if (use_atom_hash) {
    const uint64_t h = m.quick_bond_hash();

    std::unordered_set<uint64_t>* ha = atom_hash[matoms];

    const auto f = ha->find(h);

    if (f == ha->end()) {  // never seen this before, molecule must be unique
      ha->insert(h);

      return UpdateHash(m, usmi);
    }
  }

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
  if (perform_formula_check) {
    IWString mformula;
    m.formula_distinguishing_aromatic(mformula);

    // cerr << m.name() << " usmi " << usmi << "'\n";

    IW_STL_Hash_Set* s = formula_hash[matoms];

    if (s->end() == s->find(mformula)) {
      s->insert(mformula);

      // very important to update the smiles hash too
      IW_STL_Hash_Map_int* h = smiles_hash[matoms]; 
      (*h)[usmi] = 1;

      unique_molecules++;
      return 1;
    }
  }
#endif

  return UpdateHash(m, usmi);
}

/*
  No checking, no NOTHING!

  Beware if reaction queries overlap!!
*/

static int
perform_reactions(Molecule& m) {
  int rc = 0;

  for (int i = 0; i < number_reactions; i++) {
    Substructure_Results sresults;

    int nhits = reaction[i]->determine_matched_atoms(m, sresults);
    if (0 == nhits) {
      continue;
    }

    if (verbose > 2) {
      cerr << nhits << " hits for reaction " << i << '\n';
    }

    reaction[i]->in_place_transformations(m, sresults);

    rc++;
  }

  if (rc && verbose > 1) {
    cerr << m.name() << ", " << rc << " reactions changed\n";
  }
  if (rc) {
    molecules_changed_by_reactions++;
  }

  return rc;
}

// any non-zero isotope becomes `fixed`.
static int
TransformToIdenticalIsotope(Molecule& m, isotope_t fixed) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i)) {
      m.set_isotope(i, fixed);
      ++rc;
    }
  }

  return rc;
}

static void
apply_preprocessing(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (elements_to_remove.active()) {
    elements_to_remove.process(m);
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  if (number_reactions) {
    perform_reactions(m);
  }

  if (ignore_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  if (all_isotopes_become_identical) {
    TransformToIdenticalIsotope(m, kFixedIsotope);
  }

  return;
}

static int
is_unique(Molecule& m) {
  if (discard_molecule_changes) {
    Molecule tmp(m);  // make a copy
    apply_preprocessing(tmp);

    return is_unique_molecule(tmp);
  } else {
    apply_preprocessing(m);  // actually change it

    return is_unique_molecule(m);
  }
}

static int
unique_molecule(data_source_and_type<Molecule>& input) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (report_progress()) {
      cerr << " processed " << molecules_read << " molecules\n";
    }

    if (is_unique(*m)) {
      if (unique_molecule_stream.active()) {
        m->invalidate_smiles();

        if (number_assigner.active()) {
          number_assigner.process(*m);
        }

        unique_molecule_stream.write(m);
      }
    } else {
      if (duplicate_molecule_stream.active()) {
        m->invalidate_smiles();
        duplicate_molecule_stream.write(m);
      }

      if (verbose > 1) {
        cerr << "Is duplicate\n";
      }
    }
  }

  return 1;
}

// This does not seem to be necessary. TODO:ianwatson remove sometime.
static int
unique_molecule(iwstring_data_source& input, std::ostream& output) {
  IW_TDT tdt;
  while (tdt.next(input) && output.good()) {
    const_IWSubstring smi;
    if (!tdt.dataitem_value(smiles_tag, smi)) {
      cerr << "Yipes, cannot extract smiles from tdt\n";
      cerr << tdt;
      return 0;
    }

    Molecule m;
    if (!m.build_from_smiles(smi)) {
      cerr << "Very bad news, cannot parse smiles '" << smi << "'\n";
      cerr << tdt;
      return 0;
    }

    if (is_unique(m)) {
      output << tdt;
      continue;
    }

    //  If verbose we need to report the ID of the dup.

    if (verbose || duplicate_molecule_stream.active()) {
      IWString id;
      tdt.dataitem_value(identifier_tag, id);

      if (verbose > 1) {
        cerr << "Is duplicate '" << id << "'\n";
      }

      if (duplicate_molecule_stream.active()) {
        m.set_name(id);
        m.invalidate_smiles();
        duplicate_molecule_stream.write(m);
      }
    }
  }

  return output.good();
}

static int
unique_molecule(const char* fname, std::ostream& output) {
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return 0;
  }

  return unique_molecule(input, output);
}

static int
unique_molecule(const char* fname, FileType input_type) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(verbose);
  }

  return unique_molecule(input);
}

static int
build_previous_molecules(data_source_and_type<Molecule>& input, int& molecules) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules++;

    if (report_progress()) {
      cerr << " processed " << molecules << " previously known molecules\n";
    }

    (void)is_unique(*m);
  }

  return molecules;
}

static int
build_previous_molecules(const const_IWSubstring& fname, FileType input_type,
                         int& molecules) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return build_previous_molecules(input, molecules);
}

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Usage: " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << "Remove duplicate molecules from a set of molecules.\n";
  cerr << "The first instance of a structure is retained, the remaining instances are discarded,\n";
  cerr << "or perhaps written to the -D file\n";
  cerr << "Molecules may be transformed by various options\n";
  cerr << "  -l             strip to largest fragment\n";
  cerr << "  -a             compare as tautomers - skeleton and Hcount\n";
  cerr << "  -c             exclude chiral info - optical isomers will be duplicates\n";
  cerr << "  -z             exclude cis/trans bonding information\n";
  cerr << "  -I             ignore isotopic labels\n";
  cerr << "  -y             all non-zero isotopic values considered equivalent\n";
  cerr << "  -p <fname>     specify previously collected molecules\n";
  cerr << "  -s <size>      specify primary hash size (default " << default_primary_hash_size << ")\n";
  cerr << "  -S <name>      specify output file name stem\n";
  cerr << "  -D <name>      write duplicate structures to <name>\n";
  cerr << "  -R <rxn>       perform reaction(s) on molecules before comparing\n";
  cerr << "  -T             discard molecular changes after comparison\n";
  cerr << "  -r <number>    report progress every <number> molecules\n";
  cerr << "  -e             report all molecules together with counts\n";
  cerr << "  -d             gather data as DicerFragment protos, -U file is more informative\n";
  cerr << "  -U <fname>     write molecules and counts to <fname>, add '-U csv' for csv\n";
  cerr << "                 If using protos (-d) use '-U smiles' to write smiles + textproto\n";
  cerr << "  -j             items are the same only if both structure and name match\n";
  cerr << "  -t E1=E2       element transformations, enter '-t help' for details\n";
  cerr << "  -i <type>      specify input type\n";
  cerr << "  -o <type>      specify output type(s)\n";
  cerr << "  -A ...         standard aromaticity options, enter '-A help' for info\n";
  cerr << "  -K ...         standard smiles options, enter '-K help' for info\n";
  display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";

// These two options are deprecated - should never have been implemented.
//cerr << "  -f             function as filter (TDT input). Prefer -i tdt -i info -o tdt -o tdt -o info\n";
//cerr << "  -G <tag>       identifier tag when working as a filter\n";

  // clang-format on

  exit(rc);
}

static int
WriteTextproto(const resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>>& smiles_to_dicer_fragment,
               int write_smiles_and_textproto,
               IWString_and_File_Descriptor & output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  static constexpr char kSep = ' ';

  std::string buffer;

  for (const absl::flat_hash_map<IWString, dicer_data::DicerFragment>* f : smiles_to_dicer_fragment) {
    for (const auto& [usmi, proto] : *f) {
      if (! printer.PrintToString(proto, &buffer)) {
        cerr << "WriteTextproto:cannot write '" << proto.ShortDebugString() << "'\n";
        return 0;
      }
      if (write_smiles_and_textproto) {
        output << proto.smi() << kSep;
      }
      output << buffer << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  output.flush();

  return 1;
}

static int
WriteCounts(Command_Line& cl, int verbose) {
  int write_csv = 0;
  int write_smiles_and_textproto = 0;
  IWString fname;
  IWString u;
  for (int i = 0; cl.value('U', u, i); ++i) {
    if (u == "csv") {
      write_csv = 1;
    } else if (u == "smiles") {
      write_smiles_and_textproto = 1;
    } else if (fname.empty()) {
      fname = u;
    } else {
      cerr << "Unrecognised -U qualifier '" << u << "'\n";
      return 0;
    }
  }

  if (fname.empty()) {
    cerr << "WriteCounts:no file name specified\n";
    return 0;
  }

  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "WriteCounts:cannot open '" << fname << "'\n";
    return 0;
  }

  if (write_smiles_and_textproto && ! store_dicer_fragments) {
    cerr << "Writing smiles + textproto requested (-U smiles) but dicer fragments not accumulated (-d)\n";
    return 0;
  }

  if (store_dicer_fragments) {
    return WriteTextproto(smiles_to_dicer_fragment, write_smiles_and_textproto, output);
  }

  if (write_csv) {
    output << "Smiles,Count\n";
  }

  const char separator = write_csv ? ',' : ' ';

  for (const IW_STL_Hash_Map_int* by_atom : smiles_hash) {
    for (const auto& [smi, count] : *by_atom) {
      output << smi << separator << count << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  return 1;
}

static int
WriteDicerFragments(const resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>>& smiles_hash,
                    std::ostream& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const absl::flat_hash_map<IWString, dicer_data::DicerFragment>* f : smiles_hash) {
    for (const auto& [usmi, proto] : *f) {
      if (! printer.PrintToString(proto, &buffer)) {
        cerr << "DicerFragmentOutput::WriteTextcannot write '" << proto.ShortDebugString() << "'\n";
        return 0;
      }
      output << buffer << '\n';
    }
  }

  return output.good();
}

static int
WriteSmilesHash(const resizable_array_p<IW_STL_Hash_Map_int>& smiles_hash,
                std::ostream& output) {
  static constexpr char kSep = ' ';

  for (const IW_STL_Hash_Map_int* f : smiles_hash) {
    for (const auto& [usmi, count] : *f) {
      output << usmi << kSep << count << '\n';
    }
  }

  return output.good();
}

static int
WriteHash(std::ostream& output) {
  if (store_dicer_fragments) {
    return WriteDicerFragments(smiles_to_dicer_fragment, output);
  } else {
    return WriteSmilesHash(smiles_hash, output);
  }
}

static int
unique_molecule(int argc, char** argv) {
  Command_Line cl(argc, argv, "t:Tag:D:vS:A:E:X:i:o:ls:czfG:p:Ir:n:K:eR:jhU:yd");

  verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot discern elements, -E\n";
      usage(8);
    }
  }

  if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process standard aromaticity options\n";
    usage(2);
  }

  if (cl.option_present('K')) {
    if (!process_standard_smiles_options(cl, verbose, 'K')) {
      cerr << "Cannot initialise smiles options\n";
      return 5;
    }
  }

  if (cl.option_present('n')) {
    if (!number_assigner.initialise(cl, 'n', verbose > 1)) {
      cerr << "Cannot initialise number assigner (-n option)\n";
      usage(13);
    }
  }

  if (cl.option_present('j')) {
    only_same_if_structure_and_name_the_same = 1;
    if (verbose) {
      cerr << "Molecules will be considered identical only if both structure and name "
              "match\n";
    }
  }

  if (cl.option_present('h')) {
    use_atom_hash = 1;
    if (verbose) {
      cerr << "Will use atom hashes\n";
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will strip molecules to their largest fragment\n";
    }
  }

  if (cl.option_present('a')) {
    compare_as_graph = 1;
    if (verbose) {
      cerr << "molecules will be compared for tautomer equivalence\n";
    }
  }

  if (cl.option_present('c')) {
    exclude_chiral_info = 1;
    if (verbose) {
      cerr << "Optical isomers will be considered to be duplicates\n";
    }
  }

  if (cl.option_present('z')) {
    exclude_cis_trans_bonding_info = 1;
    if (verbose) {
      cerr << "Cis/trans isomers will be considered to be duplicates\n";
    }
  } else {
    exclude_cis_trans_bonding_info = 1;
  }

  if (cl.option_present('I')) {
    ignore_isotopes = 1;
    if (verbose) {
      cerr << "Isotopic variants are considered the same\n";
    }
  }

  if (cl.option_present('y')) {
    all_isotopes_become_identical = 1;
    if (verbose) {
      cerr << "Will compare isotopes as identical values\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', default_primary_hash_size) || default_primary_hash_size < 1) {
      cerr << "The -s option must be followed by a positive whole number\n";
      usage(8);
    }

    if (verbose) {
      cerr << "Primary hash size set to " << default_primary_hash_size << '\n';
    }
  }

  if (cl.option_present('X')) {
    if (!elements_to_remove.construct_from_command_line(cl, verbose, 'X')) {
      cerr << "Cannot discern elements to remove, -X\n";
      usage(3);
    }
  }

  if (cl.option_present('t')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 't')) {
      cerr << "Cannot process element transformations (-t option)\n";
      usage(32);
    }
  }

  if (cl.option_present('R')) {
    Sidechain_Match_Conditions sidechain_match_conditions;

    if (!read_reactions(cl, reaction, sidechain_match_conditions, 'R')) {
      cerr << "Cannot read reaction(s) (-R option)\n";
      return 4;
    }

    number_reactions = reaction.number_elements();

    if (verbose) {
      cerr << "Defined " << number_reactions << " reactions\n";
    }

    for (int i = 0; i < number_reactions; i++) {
      reaction[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisations (-g option)\n";
      usage(33);
    }
  }

  if (cl.option_present('T')) {
    if (!elements_to_remove.active() && !element_transformations.active() &&
        0 == number_reactions) {
      cerr << "The -T option only makes sense with the -X, -R or -t options\n";
      usage(41);
    }

    discard_molecule_changes = 1;
    if (verbose) {
      cerr << "Transformed molecule only used for comparisons\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr
          << "The report every option (-r) must be followed by a whole positive number\n";
      usage(8);
    }
  }

  // Deprecated, do not use. Prefer
  // -i tdt -i TDTID:tag -i info -o tdt -o info
  if (cl.option_present('f')) {
    if (cl.option_present('i')) {
      cerr << "The -i and -f options are incompatible\n";
      usage(18);
    }

    function_as_filter = 1;
    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }

    if (cl.option_present('G')) {
      cl.value('G', identifier_tag);
      if (!identifier_tag.ends_with('<')) {
        identifier_tag += '<';
      }

      if (verbose) {
        cerr << "Molecules identified by '" << identifier_tag << "' tag\n";
      }
    }
  }

  if (cl.option_present('d')) {
    store_dicer_fragments = 1;
    if (verbose) {
      cerr << "Data will be stored as DicerFragment protos\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_filter) {  // don't check anything about the input type
    ;
  } else if (!cl.option_present('i')) {
    if (1 == cl.number_elements() && 0 == strncmp(cl[0], "-", 1)) {
      input_type = FILE_TYPE_SMI;
    } else if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot automatically determine input type(s)\n";
      usage(8);
    }
  } else if (!process_input_type(cl, input_type)) {
    cerr << "Cannot determine input type\n";
    usage(1);
  }

  // Test this before opening any files

  if (cl.empty()) {
    cerr << "No input files specified\n";
    usage(1);
  }

  if (store_dicer_fragments) {
    for (int i = 0; i < 200; i++) {
      smiles_to_dicer_fragment << new absl::flat_hash_map<IWString, dicer_data::DicerFragment>;
    }
  } else {
    for (int i = 0; i < 200; i++) {
      smiles_hash.add(new IW_STL_Hash_Map_int);
      if (perform_formula_check) {
        formula_hash.add(new IW_STL_Hash_Set);
      }
    }
  }

  for (int i = 0; i < 200; i++) {
    //  smiles_hash[i]->resize(default_primary_hash_size);
    //  if (perform_formula_check)
    //   formula_hash.resize(default_primary_hash_size);
  }

  if (cl.option_present('p')) {
    int molecules = 0;

    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if (!build_previous_molecules(p, input_type, molecules)) {
        cerr << "Cannot process the -p option, '" << p << "'\n";
        return 72;
      }
    }

    if (verbose) {
      cerr << "read " << molecules << " molecules from previous set(s)\n";
    }

    if (duplicates_found) {
      if (verbose) {
        cerr << "Found " << duplicates_found << " duplicates in the -p file\n";
      }

      duplicates_found = 0;
    }
  }

  if (cl.option_present('S')) {
    if (function_as_filter) {
      cerr << "The -f and -S options are incompatible\n";
      usage(49);
    }

    if (!cl.option_present('o')) {
      unique_molecule_stream.add_output_type(FILE_TYPE_SMI);
    } else if (!unique_molecule_stream.determine_output_types(cl)) {
      cerr << "Cannot discern output types for duplicate stream\n";
      usage(12);
    }

    const_IWSubstring tmp;
    cl.value('S', tmp);

    if (unique_molecule_stream.would_overwrite_input_files(cl, tmp)) {
      cerr << "Cannot overwrite input file(s) '" << tmp << "'\n";
      return 7;
    }

    if (!unique_molecule_stream.new_stem(tmp, 1)) {  // causes files to be opened
      cerr << "Could not use stem '" << tmp << "' for duplicates\n";
      return 4;
    }

    if (verbose) {
      cerr << "Unique molecules written to stem '" << tmp << "'\n";
    }
  }

  if (cl.option_present('D')) {
    if (!cl.option_present('o')) {
      duplicate_molecule_stream.add_output_type(FILE_TYPE_SMI);
    } else if (!duplicate_molecule_stream.determine_output_types(cl)) {
      cerr << "Cannot discern output types for duplicate stream\n";
      usage(12);
    }

    IWString tmp;
    cl.value('D', tmp);

    if (!duplicate_molecule_stream.new_stem(tmp, 1)) {  // causes files to be opened
      cerr << "Could not use stem '" << tmp << "' for duplicates\n";
      return 4;
    }

    if (verbose) {
      cerr << "Duplicate molecules written to stem '" << tmp << "'\n";
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements() && 0 == rc; i++) {
    int rc;
    if (function_as_filter) {
      rc = unique_molecule(cl[i], std::cout);
    } else {
      rc = unique_molecule(cl[i], input_type);
    }

    if (0 == rc) {
      return i + 1;
    }
  }

  if (!unique_molecule_stream.active() && !duplicate_molecule_stream.active()) {
    verbose = 1;
  }

  if (verbose) {
    cerr << molecules_read << " molecules read, " << duplicates_found << " duplicates, "
         << (molecules_read - duplicates_found) << " unique structures\n";
    if (number_reactions) {
      cerr << molecules_changed_by_reactions << " molecules changed by reactions\n";
    }
  }

  if (cl.option_present('e')) {
    WriteHash(cerr);
  }

  if (cl.option_present('U')) {
    if (!WriteCounts(cl, verbose)) {
      cerr << "Cannot write hash (-U)\n";
      return 1;
    }
  }

  if (use_atom_hash) {
    unsigned int s = 0;
    for (int i = 0; i < atom_hash.number_elements(); ++i) {
      s += atom_hash[i]->size();
    }
    cerr << "Atom hash contains " << s << " discrete values\n";
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = unique_molecule(argc, argv);

  return rc;
}
