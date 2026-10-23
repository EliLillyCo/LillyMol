// CLI for the reusable rooted multi-arm fingerprint generator.
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "fingerprint_writer.h"
#include "multi_arm_fingerprint.h"

namespace multi_arm_fingerprint {
namespace {

void
Usage(int rc) {
  std::cerr << "Rooted three- and four-arm sparse fingerprints\n"
               " -S     synchronised arm expansion (default all length combinations)\n"
               " -R <n> maximum total bond length (default 6)\n"
               " -r <n> minimum total bond length (default 3)\n"
               " -L <n> maximum individual arm length (default total budget)\n"
               " -d <n> centre degree: 3 or 4 (default both)\n"
               " -c     omit ring-closing extensions\n"
               " -b     presence only, rather than occurrence counts\n"
               " -P ... atom typing (default UST:Y)\n"
               " -s ... SMARTS selecting roots via the first matched atom\n"
               " -q ... query files selecting roots via the first matched atom\n"
               " -z i   allow query nonmatches, writing empty fingerprints\n"
               " -J ... fingerprint output specification (default NCMARM)\n"
               " -f     TDT filter mode\n"
               " -l     largest fragment\n"
               " -g ... chemical standardisation\n"
               " -i ... input type (- is stdin; defaults to SMILES for stdin)\n"
               " -A ... aromaticity options\n"
               " -E ... element options\n"
               " -v     report enumeration statistics\n";
  exit(rc);
}

class Application {
 public:
  Generator _generator;
  Chemical_Standardisation _standardisation;
  resizable_array_p<Substructure_Query> _queries;
  fingerprint_writer::FingerprintWriter _writer;
  bool _largest_fragment = false;
  bool _filter = false;
  bool _ignore_nonmatches = false;
  uint64_t _molecules_not_matching_queries = 0;
  uint64_t _molecules = 0;
  Statistics _totals;
  Accumulator<double> _bits_per_molecule;

  // Queries select roots only; a matched root need not have eligible degree.
  int
  SelectRoots(Molecule& m, std::vector<int>& roots) {
    roots.assign(m.natoms(), 0);
    bool matched = false;
    Molecule_to_Match target(&m);
    for (Substructure_Query* query : _queries) {
      Substructure_Results matches;
      if (!query->substructure_search(target, matches)) {
        continue;
      }
      matched = true;
      for (const Set_of_Atoms* embedding : matches.embeddings()) {
        if (embedding->empty() || !m.ok_atom_number(embedding->front())) {
          std::cerr << "Query did not provide a valid first matched atom\n";
          return 0;
        }
        roots[embedding->front()] = 1;
      }
    }
    if (matched) {
      return 1;
    }
    ++_molecules_not_matching_queries;
    if (_ignore_nonmatches) {
      return 1;
    }
    std::cerr << m.smiles() << ' ' << m.name()
              << " no query matches (use -z i to allow)\n";
    return 0;
  }

  int
  Process(Molecule& m, IWString_and_File_Descriptor& output) {
    if (_largest_fragment) {
      m.reduce_to_largest_fragment();
    }
    if (_standardisation.active()) {
      _standardisation.process(m);
    }
    std::vector<int> roots;
    if (!_queries.empty() && !SelectRoots(m, roots)) {
      return 0;
    }
    Sparse_Fingerprint_Creator fingerprint;
    Statistics stats;
    // For an empty molecule a null mask is equivalent to an empty selection.
    if (!_generator.Fingerprint(m, fingerprint, roots.empty() ? nullptr : roots.data(),
                                &stats)) {
      return 0;
    }
    ++_molecules;
    _totals._centres += stats._centres;
    _totals._arm_paths += stats._arm_paths;
    _totals._combinations += stats._combinations;
    _totals._distinct_bits += stats._distinct_bits;
    _bits_per_molecule.extra(fingerprint.nbits());
    const bool tdt = !_writer.IsWritingDescriptors() && !_writer.IsWritingSvml();
    if (!_filter && tdt) {
      output << "$SMI<" << m.smiles() << ">\nPCN<" << m.name() << ">\n";
    }
    if (!_writer.WriteFingerprint(m.name(), fingerprint, output)) {
      return 0;
    }
    if (!_filter && tdt) {
      output << "|\n";
    }
    output.write_if_buffer_holds_more_than(32768);
    return 1;
  }

  int
  ReadFilter(const char* fname, IWString_and_File_Descriptor& output) {
    iwstring_data_source input(fname);
    if (!input.good()) {
      return 0;
    }
    const_IWSubstring line;
    while (input.next_record(line)) {
      output << line << '\n';
      if (!line.starts_with("$SMI<")) {
        continue;
      }
      if (!line.ends_with('>')) {
        return 0;
      }
      line.remove_leading_chars(5);
      line.chop();
      Molecule m;
      if (!m.build_from_smiles(line) || !Process(m, output)) {
        return 0;
      }
    }
    return 1;
  }

  int
  Read(const char* fname, FileType type, IWString_and_File_Descriptor& output) {
    if (type == FILE_TYPE_INVALID) {
      type = IWString(fname) == "-" ? FILE_TYPE_SMI : discern_file_type_from_name(fname);
    }
    if (type == FILE_TYPE_INVALID) {
      return 0;
    }
    data_source_and_type<Molecule> input(type, fname);
    if (!input.good()) {
      return 0;
    }
    while (Molecule* raw = input.next_molecule()) {
      std::unique_ptr<Molecule> m(raw);
      if (!Process(*m, output)) {
        return 0;
      }
    }
    return 1;
  }
};

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vSR:r:L:d:cbP:s:q:J:flg:i:A:E:z:");
  if (cl.unrecognised_options_encountered() || cl.empty()) {
    Usage(1);
  }
  const int verbose = cl.option_count('v');
  if (!process_elements(cl)) {
    return 1;
  }
  set_global_aromaticity_type(Daylight);
  if (!process_standard_aromaticity_options(cl, verbose)) {
    return 1;
  }
  Application app;
  const_IWSubstring z;
  for (int i = 0; cl.value('z', z, i); ++i) {
    if (z != "i") {
      std::cerr << "Unrecognised -z qualifier '" << z << "'\n";
      return 1;
    }
    app._ignore_nonmatches = true;
  }
  Options options;
  for (auto item : {std::pair<char, int*>{'R', &options._max_total_length},
                    {'r', &options._min_total_length},
                    {'L', &options._max_arm_length},
                    {'d', &options._centre_degree}}) {
    if (cl.option_present(item.first) &&
        (!cl.value(item.first, *item.second) || *item.second <= 0)) {
      std::cerr << "Invalid numeric option -" << item.first << '\n';
      return 1;
    }
  }
  options._ring_closures = !cl.option_present('c');
  options._presence_only = cl.option_present('b');
  options._synchronised = cl.option_present('S');
  if (!app._generator.SetOptions(options)) {
    std::cerr << "Incompatible length limits or centre degree\n";
    return 1;
  }
  if (cl.option_present('P') && !app._generator.SetAtomTyping(cl.string_value('P'))) {
    return 1;
  }
  if (cl.option_present('g') &&
      !app._standardisation.construct_from_command_line(cl, verbose, 'g')) {
    return 1;
  }
  app._largest_fragment = cl.option_present('l');
  app._filter = cl.option_present('f');
  if (cl.option_present('q') && !process_queries(cl, app._queries, verbose)) {
    return 1;
  }
  const_IWSubstring smarts;
  for (int i = 0; cl.value('s', smarts, i); ++i) {
    auto query = std::make_unique<Substructure_Query>();
    if (!query->create_from_smarts(smarts)) {
      return 1;
    }
    app._queries.add(query.release());
  }
  for (Substructure_Query* query : app._queries) {
    // Keep all root orientations, including embeddings on the same atom set.
    query->set_find_unique_embeddings_only(0);
    query->set_perceive_symmetry_equivalent_matches(1);
    query->set_find_one_embedding_per_atom(0);
    query->set_max_matches_to_find(std::numeric_limits<int>::max());
    query->set_save_matched_atoms(1);
  }
  if (cl.option_present('J')) {
    if (!app._writer.Initialise(cl, 'J', verbose)) {
      return 1;
    }
  } else {
    app._writer.SetSparseOutput("NCMARM");
  }
  if (app._filter &&
      (app._writer.IsWritingDescriptors() || app._writer.IsWritingSvml())) {
    std::cerr << "TDT filter mode requires fingerprint output\n";
    return 1;
  }
  FileType type = FILE_TYPE_INVALID;
  if (cl.option_present('i') && !process_input_type(cl, type)) {
    return 1;
  }
  IWString_and_File_Descriptor output(1);
  app._writer.WriteHeaderIfNeeded(output);
  for (const char* fname : cl) {
    const int rc =
        app._filter ? app.ReadFilter(fname, output) : app.Read(fname, type, output);
    if (!rc) {
      std::cerr << "Cannot fingerprint '" << fname << "'\n";
      return 1;
    }
  }
  output.flush();
  if (verbose) {
    std::cerr << "Processed " << app._molecules << " molecules, " << app._totals._centres
              << " centres, " << app._totals._arm_paths << " arm paths, "
              << app._totals._combinations << " combinations, "
              << app._totals._distinct_bits << " summed distinct bits\n";
    if (!app._queries.empty()) {
      std::cerr << app._molecules_not_matching_queries
                << " molecules did not match any queries\n";
    }
    if (!app._bits_per_molecule.empty()) {
      std::cerr << "Distinct bits per molecule: min " << app._bits_per_molecule.minval()
                << " max " << app._bits_per_molecule.maxval() << " average "
                << app._bits_per_molecule.average() << '\n';
    }
  }
  return 0;
}
}  // namespace
}  // namespace multi_arm_fingerprint

int
main(int argc, char** argv) {
  return multi_arm_fingerprint::Main(argc, argv);
}
