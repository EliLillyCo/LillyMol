/*
  Sometimes we need to determine the unique molecules from a file
  of structures.
*/

#include <stdlib.h>

#include <cstdint>
#include <iostream>
#include <memory>

#include <sys/stat.h>

using std::cerr;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
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

#include "absl/container/flat_hash_map.h"

#if defined(BUILD_INCHI)
#include "inchi_api.h"
#endif

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

static int verbose = 0;

const char* prog_name = nullptr;

static Element_Transformations element_transformations;

/*
  How should element transformations and removals be used?

  By default, we apply the removal and transformations and the resulting molecule
  is what gets written. But, what we probably want is for the transformed
  molecule to be used for comparisons only, and we want the original
  molecule written
*/

static int discard_molecule_changes = 0;

static Chemical_Standardisation chemical_standardisation;

static uint64_t molecules_read = 0;

static uint64_t duplicates_found = 0;

static uint64_t unique_molecules = 0;

static int reduce_to_largest_fragment = 0;

static Molecule_Output_Object unique_molecule_stream;

static Molecule_Output_Object duplicate_molecule_stream;

static int compare_as_graph = 0;

static int ignore_isotopes = 0;

// If there is an isotope, transform to a single value.
static int all_isotopes_become_identical = 0;

static constexpr isotope_t kFixedIsotope = 1;

static int exclude_chiral_info = 0;

static int exclude_cis_trans_bonding_info = 0;

static Report_Progress report_progress;

/*
  Jul 2005. Jibo wants to be able to enable the case where same
  structure but different name will be considered different
*/

static int only_same_if_structure_and_name_the_same = 0;

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
// Seemed like an interesting idea, but ultimately slows things down.
static int perform_formula_check = 0;
#endif

// Use quick_bond_hash for each molecule and store in a hash.
// Might be slightly quicker, but seems to not make much of a
// difference.
// Currently this is an undocumented option.
static int use_atom_hash = 0;

/*
  Here's where we keep track of the previously encountered smiles
  Hashes get less and less efficient as they get larger, so have
  a different hash for every atom count
*/

static resizable_array_p<IW_STL_Hash_Set> formula_hash;
static resizable_array_p<IW_STL_Hash_Map_int> smiles_hash;
static resizable_array_p<std::unordered_set<uint64_t>> atom_hash;

static int write_dicer_textproto = 0;
static resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>>
    smiles_to_dicer_fragment;

template class resizable_array_p<IW_STL_Hash_Map_int>;
template class resizable_array_base<IW_STL_Hash_Map_int*>;
template class resizable_array_p<IW_STL_Hash_Set>;
template class resizable_array_base<IW_STL_Hash_Set*>;

static resizable_array_p<IWReaction> reaction;
static int number_reactions = 0;
static int molecules_changed_by_reactions = 0;

#ifdef BUILD_INCHI
static int use_inchi = 0;
#endif

static int
UpdateDicerFragmentHash(Molecule& m, const IWString& usmi) {
  const int matoms = m.natoms();

  absl::flat_hash_map<IWString, dicer_data::DicerFragment>* smiles_hash =
      smiles_to_dicer_fragment[matoms];

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
  if (write_dicer_textproto) {
    return UpdateDicerFragmentHash(m, usmi);
  } else {
    return UpdateSmilesHash(m, usmi);
  }
}

#ifdef BUILD_INCHI
static int
IsUniqueInChI(Molecule& m) {
  IWString inchi;
  m.InChI(inchi);

  IWString inchi_key;
  InChIToInChIKey(inchi.data(), inchi_key);

  return UpdateSmilesHash(m, inchi_key);
}
#endif

static int
is_unique_molecule(Molecule& m) {
#ifdef BUILD_INCHI
  if (use_inchi) {
    return IsUniqueInChI(m);
  }
#endif
  if (exclude_chiral_info) {
    set_include_chiral_info_in_smiles(0);
  }

  if (exclude_cis_trans_bonding_info) {
    set_include_cis_trans_in_smiles(0);
  }

  IWString usmi;

  if (compare_as_graph) {
    IWString formula = m.molecular_formula();

    Molecule g(m);
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

  if (write_dicer_textproto) {
    while (matoms >= smiles_to_dicer_fragment.number_elements()) {
      smiles_to_dicer_fragment
          << new absl::flat_hash_map<IWString, dicer_data::DicerFragment>;
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
    assert(ha != nullptr);

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
    m.reduce_to_largest_fragment_carefully();
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
  cerr << R"(Remove duplicate molecules from a set of molecules.\n";
The first instance of a structure is retained, the remaining instances are discarded,
or perhaps written to the -D file.
Molecules may be transformed by various options - the first options below.
 -l             strip to largest fragment.
 -a             compare as graph/tautomers - concatenate molecular skeleton and formula.
 -c             exclude chiral info - optical isomers will be duplicates.
 -z             exclude cis/trans bonding information - true by default.
 -I             ignore isotopic labels.
 -y             all non-zero isotopic values considered equivalent.
 -R <rxn>       perform reaction(s) on molecules before comparing.
 -t E1=E2       element transformations, enter '-t help' for details.
 -j             items are the same only if both structure and name match.
)";
#ifdef BUILD_INCHI
  cerr << " -C             use InChi for structure comparisons.\n";
#endif
  cerr << R"(
 -p <fname>     specify previously collected molecules.
 -S <name>      specify output file name stem - the unique molecules.
 -D <name>      write duplicate structures to <name>.
 -T             discard molecular changes after comparison - writes initial molecule.
 -U <fname>     write molecules and counts to <fname>. for csv '-U csv' for dicer fragments '-U textproto'.
 -U csv         write the -U file a csv form.
 -U textproto   write the -U file as dicer_data::DicerFragment textproto
 -U smiles      if writing textprotos, write smiles + textproto
  
 -r <number>    report progress every <number> molecules.
 -i <type>      specify input type.
 -o <type>      specify output type(s).
 -E ...         standard element options, enter '-E help' for info.
 -A ...         standard aromaticity options, enter '-A help' for info.
 -K ...         standard smiles options, enter '-K help' for info.
 -g ...         chemical standardisation options, enter '-g help' for info.
 -v             verbose output.

With no arguments and no output, reports the number of duplicates.
unique_molecules file.smi

Specify structure comparison options and gather output.
unique_molecules -c -l -g all -I -r 10000 -S unique haystack.smi
)";

  // clang-format on

  exit(rc);
}

static int
WriteTextproto(
    const resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>>&
        smiles_to_dicer_fragment,
    int write_smiles_and_textproto, IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  static constexpr char kSep = ' ';

  std::string buffer;

  for (const absl::flat_hash_map<IWString, dicer_data::DicerFragment>* f :
       smiles_to_dicer_fragment) {
    for (const auto& [usmi, proto] : *f) {
      if (!printer.PrintToString(proto, &buffer)) {
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
      write_dicer_textproto = 1;
      if (verbose) {
        cerr << "Will write smiles + textproto to the -U file\n";
      }
    } else if (u == "textproto") {
      write_dicer_textproto = 1;
      if (verbose) {
        cerr << "Will write dicer_data::DicerFragment textprotos to the -U file\n";
      }
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

  if (write_dicer_textproto) {
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
WriteDicerFragments(
    const resizable_array_p<absl::flat_hash_map<IWString, dicer_data::DicerFragment>>&
        smiles_hash,
    std::ostream& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const absl::flat_hash_map<IWString, dicer_data::DicerFragment>* f : smiles_hash) {
    for (const auto& [usmi, proto] : *f) {
      if (!printer.PrintToString(proto, &buffer)) {
        cerr << "DicerFragmentOutput::WriteTextcannot write '" << proto.ShortDebugString()
             << "'\n";
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

#ifdef FOR_DEBUGGING
static int
WriteHash(std::ostream& output) {
  if (write_dicer_textproto) {
    return WriteDicerFragments(smiles_to_dicer_fragment, output);
  } else {
    return WriteSmilesHash(smiles_hash, output);
  }
}
#endif

static int
unique_molecule(int argc, char** argv) {
  // Note that the -C option for InChI is always present
  Command_Line cl(argc, argv, "t:Tag:D:vS:A:E:i:o:lczp:fIr:K:R:jhU:yC");

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

#ifdef BUILD_INCHI
  if (cl.option_present('C')) {
    use_inchi = 1;
    if (verbose) {
      cerr << "Structure comparisons done via InChI\n";
    }
  }
#endif

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
    if (!element_transformations.active() && 0 == number_reactions && reduce_to_largest_fragment == 0) {
      cerr << "The -T option only makes sense with the -R or -t options\n";
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

  FileType input_type = FILE_TYPE_INVALID;
  if (!cl.option_present('i')) {
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

  if (write_dicer_textproto) {
    for (int i = 0; i < 200; i++) {
      smiles_to_dicer_fragment
          << new absl::flat_hash_map<IWString, dicer_data::DicerFragment>;
    }
  } else {
    for (int i = 0; i < 200; i++) {
      smiles_hash.add(new IW_STL_Hash_Map_int);
#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
      if (perform_formula_check) {
        formula_hash.add(new IW_STL_Hash_Set);
      }
#endif
      if (use_atom_hash) {
        atom_hash.add(new std::unordered_set<uint64_t>());
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

  // We must determine before processing whether to store textprotos in the -U file
  // The -U option is fully parsed in WriteCounts
  if (cl.option_present('U')) {
    const_IWSubstring u;
    for (int i = 0; cl.value('U', u, i); ++i) {
      if (u == "textproto") {
        write_dicer_textproto = 1;
        if (verbose) {
          cerr << "Will write dicer_data::DicerFragment textprotos to the -U file\n";
        }
      } else if (u == "smiles") {
        write_dicer_textproto = 1;
        if (verbose) {
          cerr << "Will write smiles + textproto to the -U file\n";
        }
      }
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
    if (duplicate_molecule_stream.would_overwrite_input_files(cl, tmp)) {
      cerr << "Cannot overwrite input '" << tmp << "'\n";
      return 1;
    }

    if (!duplicate_molecule_stream.new_stem(tmp, 1)) {  // causes files to be opened
      cerr << "Could not use stem '" << tmp << "' for duplicates\n";
      return 4;
    }

    if (verbose) {
      cerr << "Duplicate molecules written to stem '" << tmp << "'\n";
    }
  }

  for (const char* fname : cl) {
    if (!unique_molecule(fname, input_type)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
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
