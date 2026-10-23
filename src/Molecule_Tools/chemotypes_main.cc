#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/md5.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule_preprocessing.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/util.h"

#include "Molecule_Tools/chemotypes.h"
#include "Molecule_Tools/iwecfp_lib.h"
#include "Molecule_Tools/dicer_fragments.pb.h"

#include "absl/container/flat_hash_map.h"
#include "google/protobuf/text_format.h"

namespace chemotypes_main {

using std::cerr;
using molecule_processing::MoleculePreprocessing;

constexpr int kDefaultHashWidth = 26;
constexpr int kMaxHashWidth = 26;

IWString
Base32Hash(const unsigned char* digest, int width) {
  static constexpr char kAlphabet[] = "0123456789ABCDEFGHJKMNPQRSTVWXYZ";

  IWString result;
  int bit = 0;
  for (int i = 0; i < width; ++i) {
    int value = 0;
    for (int j = 0; j < 5; ++j, ++bit) {
      value <<= 1;
      if (bit >= 128) {
        continue;
      }
      const int byte = bit / 8;
      const int offset = 7 - (bit % 8);
      value |= (digest[byte] >> offset) & 0x01;
    }
    result << kAlphabet[value];
  }

  return result;
}

IWString
ChemotypeHashString(const IWString& smiles, int width) {
  unsigned char digest[16];
  MD5_CTX ctx;
  MD5Init(&ctx);
  MD5Update(&ctx, reinterpret_cast<unsigned char*>(const_cast<char*>(smiles.data())),
            smiles.length());
  MD5Final(digest, &ctx);

  return Base32Hash(digest, width);
}

struct ChemotypeHashXref {
  IWString smiles;
  IWString parent;
  uint32_t n = 0;
};


enum class FingerprintKind {
  kNone,
  kFixed,
  kNonColliding,
};

struct ChemotypeFingerprintSpec {
  FingerprintKind kind = FingerprintKind::kNone;
  IWString tag;
  int ec_radius = 3;
  int expand = 0;
  int invert = 0;
};

int
EnsureFingerprintTag(IWString& tag) {
  if (! tag.ends_with('<')) {
    tag << '<';
  }

  return 1;
}

int
TrailingNumber(const IWString& tag, int default_value) {
  const int nchars = tag.length();
  int end = nchars - 1;
  if (end >= 0 && tag[end] == '<') {
    --end;
  }

  if (end < 0 || ! std::isdigit(static_cast<unsigned char>(tag[end]))) {
    return default_value;
  }

  int start = end;
  while (start > 0 && std::isdigit(static_cast<unsigned char>(tag[start - 1]))) {
    --start;
  }

  int result = 0;
  for (int i = start; i <= end; ++i) {
    result = 10 * result + tag[i] - '0';
  }

  return result;
}

int
ParseDashJQualifier(const const_IWSubstring& qualifier,
                    ChemotypeFingerprintSpec& spec) {
  if (qualifier.starts_with("EXPAND=")) {
    const_IWSubstring value(qualifier);
    value.remove_leading_chars(7);
    if (! value.numeric_value(spec.expand) || spec.expand < 0) {
      cerr << "Invalid -J EXPAND value '" << qualifier << "'\n";
      return 0;
    }
    return 1;
  }

  if (qualifier == "INVERT" || qualifier == "OUTSIDE") {
    spec.invert = 1;
    return 1;
  }

  if (qualifier.starts_with("FP")) {
    if (spec.kind != FingerprintKind::kNone) {
      cerr << "Only one -J fingerprint tag can be specified\n";
      return 0;
    }
    spec.kind = FingerprintKind::kFixed;
    spec.tag = qualifier;
    return EnsureFingerprintTag(spec.tag);
  }

  if (qualifier.starts_with("NC")) {
    if (spec.kind != FingerprintKind::kNone) {
      cerr << "Only one -J fingerprint tag can be specified\n";
      return 0;
    }
    spec.kind = FingerprintKind::kNonColliding;
    spec.tag = qualifier;
    spec.ec_radius = TrailingNumber(spec.tag, 3);
    return EnsureFingerprintTag(spec.tag);
  }

  cerr << "Unrecognised -J qualifier '" << qualifier << "'\n";
  return 0;
}


int
ParseDashYQualifier(const const_IWSubstring& qualifier,
                    chemotypes::ChemotypeOptions& options) {
  if (qualifier == "larger" || qualifier == "largest" || qualifier == "size") {
    options.prefer_larger_adjacent_ring_system = true;
    return 1;
  }

  cerr << "Unrecognised -Y qualifier '" << qualifier << "'\n";
  return 0;
}

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on

  cerr << R"(Extracts query-defined molecular chemotypes.
The first query that matches, and the first embedding of that query, define the seed ring system.

Options:
 -q <query>     LillyMol query file defining the chemotype seed.
 -s <smarts>    SMARTS defining the chemotype seed.
 -n <n>         include <n> directly adjacent ring systems, default 0.
 -r <n>         minimum number of rings required for processing, default 0.
 -D <dist>      maximum bond distance to an adjacent ring system.
 -u             include one-hop atoms attached to retained ring atoms.
 -x             ignore singly connected attached atoms with -u.
 -t             with -n, include all adjacent ring systems tied at the cutoff distance.
 -P <atype>     atom typing specification; non-terminal attachment atoms become labelled.
 -I <iso>       label both endpoints of chemotype boundary bonds with isotope <iso>; incompatible with -P.
 -C <stem>      write excluded atoms to a separate output stem; uses -o types.
 -J FP<tag>     generate fixed-width linear fingerprint of chemotype atoms.
 -J NC<tag>     generate non-colliding EC fingerprint of chemotype atoms.
                A trailing digit sets the EC radius, default 3.
 -J EXPAND=<n>  expand selected fingerprint atoms by <n> bonds.
 -J INVERT      fingerprint atoms outside the chemotype; OUTSIDE is a synonym.
 -Y larger      prefer larger same-distance adjacent ring-system candidates.
                Synonyms: largest, size.
 -H def         append a deterministic fixed-width alphanumeric hash of the chemotype.
 -H xref=<fname> also write hash, unique smiles, first parent, and count cross-reference.
 -H width=<n>   hash width, default 26, max 26.
 -f             work as a TDT filter; requires -J.
 -p <text>      write the parent before each chemotype; append <text> to parent name.
                Use -p . or -p def for no parent annotation.
 -F <fname>      write accumulated dicer_data::DicerFragment textproto summary.
 -U <prefix>     with -F, assign DicerFragment.id and append it to chemotype names.
                 Use -U . or -U def to append the bare id.
 -z i|ignore    ignore molecules not matching any query.
 -z f|first     use the first embedding if a query matches multiple ring systems.
 -S <stem>      output file stem; default stdout.
 -i <type>      input type.
 -o <type>      output type, default smiles.
 -g ...         chemical standardisation.
 -l             reduce to largest fragment.
 -c             remove chirality.
 -A ...         aromaticity options.
 -E ...         element options.
 -v             verbose output.
)";

  ::exit(rc);
}

class Options {
 private:
  int _verbose = 0;

  MoleculePreprocessing _preprocessing;

  int _min_rings = 0;
  int _ignore_molecules_not_matching_query = 0;
  int _write_parent = 0;
  Molecule_Output_Object _complement_output;
  bool _write_complement = false;
  uint64_t _complements_written = 0;
  IWString _parent_annotation;
  IWString _summary_fname;
  int _hash_chemotypes = 0;
  int _hash_width = kDefaultHashWidth;
  IWString _hash_xref_fname;
  int _write_unique_chemotype_id = 0;
  IWString _unique_chemotype_id_prefix;

  int _function_as_tdt_filter = 0;
  ChemotypeFingerprintSpec _fingerprint;

  resizable_array_p<Substructure_Query> _queries;
  chemotypes::ChemotypeOptions _chemotype_options;
  chemotypes::ChemotypeScratch _scratch;

  Atom_Typing_Specification _atom_typing;
  Atom_Typing_Specification* _atom_typing_ptr = nullptr;

  absl::flat_hash_map<IWString, dicer_data::DicerFragment> _chemotype;
  absl::flat_hash_map<IWString, IWString> _hash_from_usmi;
  absl::flat_hash_map<IWString, ChemotypeHashXref> _hash_xref;
  std::vector<IWString> _hash_order;

  uint64_t _molecules_read = 0;
  uint64_t _molecules_written = 0;
  uint64_t _parent_molecules_written = 0;
  uint64_t _unique_chemotypes = 0;
  uint64_t _duplicate_chemotypes = 0;
  uint64_t _molecules_not_matching_queries = 0;
  uint64_t _molecules_below_min_rings = 0;
  uint64_t _molecules_matched_query_without_ring_atom = 0;
  uint64_t _molecules_with_ambiguous_query_matches = 0;
  uint64_t _atom_typing_failures = 0;

 public:
  int Initialise(Command_Line& cl);
  int Preprocess(Molecule& m);
  int Process(Molecule& m, Molecule_Output_Object& output);
  int ProcessFingerprint(Molecule& m, IWString_and_File_Descriptor& output);
  int WriteFingerprint(Molecule& m, const int* include_atom,
                       const uint32_t* atom_type, IWString_and_File_Descriptor& output);
  int WriteEmptyFingerprint(IWString_and_File_Descriptor& output);
  int StartFingerprintRecord(Molecule& m, IWString_and_File_Descriptor& output);
  int FinishFingerprintRecord(Molecule& m, IWString_and_File_Descriptor& output);
  int AccumulateChemotype(Molecule& m, const IWString& parent_name, uint32_t& chemotype_id);
  int ChemotypeHash(Molecule& m, const IWString& parent_name, IWString& hash);
  int WriteSummary();
  int WriteHashXref();
  int Report(std::ostream& output) const;

  int verbose() const {
    return _verbose;
  }

  int fingerprinting_active() const {
    return _fingerprint.kind != FingerprintKind::kNone;
  }

  int function_as_tdt_filter() const {
    return _function_as_tdt_filter;
  }
};

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! _preprocessing.Initialise(cl)) {
    cerr << "Cannot initialise molecule preprocessing\n";
    return 0;
  }

  if (! cl.option_present('q') && ! cl.option_present('s')) {
    cerr << "Must specify one or more queries via -q or -s\n";
    return 0;
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> query = std::make_unique<Substructure_Query>();
      if (! query->create_from_smarts(smarts)) {
        cerr << "Invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << query.release();
    }
  }

  for (Substructure_Query* query : _queries) {
    query->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _chemotype_options.adjacent_ring_systems_to_include) ||
        _chemotype_options.adjacent_ring_systems_to_include < 0) {
      cerr << "The adjacent ring system count (-n) must be a non-negative integer\n";
      return 0;
    }
  }


  if (cl.option_present('r')) {
    if (! cl.value('r', _min_rings) || _min_rings < 0) {
      cerr << "The minimum ring count (-r) must be a non-negative integer\n";
      return 0;
    }
  }

  if (cl.option_present('D')) {
    if (! cl.value('D', _chemotype_options.max_distance) ||
        _chemotype_options.max_distance <= 0) {
      cerr << "The max distance (-D) must be a positive integer\n";
      return 0;
    }
  }

  if (cl.option_present('u')) {
    _chemotype_options.include_attached_atoms = 1;
  }

  if (cl.option_present('J')) {
    const_IWSubstring j;
    for (int i = 0; cl.value('J', j, i); ++i) {
      if (! ParseDashJQualifier(j, _fingerprint)) {
        return 0;
      }
    }
    if (_fingerprint.kind == FingerprintKind::kNone) {
      cerr << "The -J option must specify one fingerprint tag starting with FP or NC\n";
      return 0;
    }
  }


  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (! ParseDashYQualifier(y, _chemotype_options)) {
        return 0;
      }
    }
  }

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (! fingerprinting_active()) {
      cerr << "The TDT filter option (-f) requires fingerprint generation (-J)\n";
      return 0;
    }
  }

  if (fingerprinting_active() &&
      (cl.option_present('o') || cl.option_present('S') || cl.option_present('p') ||
       cl.option_present('F') || cl.option_present('U') || cl.option_present('I') ||
       cl.option_present('H') || cl.option_present('C'))) {
    cerr << "Fingerprint generation (-J) cannot be combined with molecule or summary output options\n";
    return 0;
  }

  if (cl.option_present('C')) {
    IWString stem = cl.string_value('C');
    if (stem == "-" || (cl.option_present('S') && stem == cl.string_value('S'))) {
      cerr << "Complement output must use a separate file stem\n";
      return 0;
    }
    if (! cl.option_present('o')) {
      _complement_output.add_output_type(FILE_TYPE_SMI);
    } else if (! _complement_output.determine_output_types(cl, 'o')) {
      return 0;
    }
    if (_complement_output.would_overwrite_input_files(cl, stem) ||
        ! _complement_output.new_stem(stem)) {
      cerr << "Cannot open complement output stem '" << stem << "'\n";
      return 0;
    }
    _write_complement = true;
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    const_IWSubstring p = cl.string_value('p');
    if (p != '.' && p != "def" && p != "default") {
      _parent_annotation << ' ' << p;
    }
  }

  if (cl.option_present('F')) {
    _summary_fname = cl.string_value('F');
  }


  if (cl.option_present('H')) {
    _hash_chemotypes = 1;
    IWString h;
    for (int i = 0; cl.value('H', h, i); ++i) {
      if (h == "def" || h == "default" || h == "hash") {
        continue;
      }
      if (h.starts_with("xref=")) {
        h.remove_leading_chars(5);
        if (h.empty()) {
          cerr << "The -H xref= qualifier requires a file name\n";
          return 0;
        }
        _hash_xref_fname = h;
      } else if (h.starts_with("width=")) {
        h.remove_leading_chars(6);
        if (! h.numeric_value(_hash_width) || _hash_width <= 0 ||
            _hash_width > kMaxHashWidth) {
          cerr << "The -H width= qualifier must be between 1 and " << kMaxHashWidth << '\n';
          return 0;
        }
      } else {
        cerr << "Unrecognised -H qualifier '" << h << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('U')) {
    _write_unique_chemotype_id = 1;
    const_IWSubstring u = cl.string_value('U');
    if (u != '.' && u != "def" && u != "default") {
      _unique_chemotype_id_prefix = u;
    }
  }

  if (_write_unique_chemotype_id && _summary_fname.empty()) {
    cerr << "The unique chemotype id option (-U) requires summary accumulation (-F)\n";
    return 0;
  }

  if (cl.option_present('x')) {
    _chemotype_options.ignore_singly_connected_attached_atoms = 1;
  }

  if (cl.option_present('t')) {
    _chemotype_options.include_tied_adjacent_ring_systems = 1;
  }


  if (cl.option_present('I')) {
    uint64_t iso;
    if (! cl.value('I', iso) || iso == 0) {
      cerr << "The exit point isotope (-I) must be a positive integer\n";
      return 0;
    }
    _chemotype_options.isotope_for_exit_points = static_cast<isotope_t>(iso);
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Invalid atom typing specification (-P) '" << p << "'\n";
      return 0;
    }
    _atom_typing_ptr = &_atom_typing;
  } else if (fingerprinting_active()) {
    const_IWSubstring p("UST:ARY");
    if (! _atom_typing.build(p)) {
      cerr << "Cannot initialise default fingerprint atom typing '" << p << "'\n";
      return 0;
    }
    _atom_typing_ptr = &_atom_typing;
  }


  if (_chemotype_options.isotope_for_exit_points != 0 && _atom_typing_ptr != nullptr &&
      ! fingerprinting_active()) {
    cerr << "The fixed exit point isotope option (-I) cannot be used with atom typing (-P)\n";
    return 0;
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i' || z == "ignore") {
        _ignore_molecules_not_matching_query = 1;
      } else if (z == 'f' || z == "first") {
        _chemotype_options.choose_first_embedding = true;
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (_verbose) {
    cerr << "Read " << _queries.number_elements() << " queries\n";
    cerr << "Will include " << _chemotype_options.adjacent_ring_systems_to_include
         << " adjacent ring systems\n";
    if (_min_rings > 0) {
      cerr << "Will only process molecules with at least " << _min_rings << " rings\n";
    }
    if (_chemotype_options.max_distance) {
      cerr << "Adjacent ring systems limited to distance "
           << _chemotype_options.max_distance << '\n';
    }
    if (_chemotype_options.include_attached_atoms) {
      cerr << "Will include one-hop atoms attached to retained ring atoms\n";
    }
    if (_chemotype_options.ignore_singly_connected_attached_atoms) {
      cerr << "Will ignore singly connected attached atoms\n";
    }
    if (_write_parent) {
      cerr << "Will write parent molecules before chemotypes";
      if (! _parent_annotation.empty()) {
        cerr << " with parent annotation '" << _parent_annotation << "'";
      }
      cerr << '\n';
    }
    if (! _summary_fname.empty()) {
      cerr << "Will write chemotype summary to '" << _summary_fname << "'\n";
    }
    if (_write_unique_chemotype_id) {
      cerr << "Will append unique chemotype ids to chemotype names";
      if (! _unique_chemotype_id_prefix.empty()) {
        cerr << " with prefix '" << _unique_chemotype_id_prefix << "'";
      }
      cerr << '\n';
    }
    if (_hash_chemotypes) {
      cerr << "Will append " << _hash_width << " character chemotype hashes\n";
      if (! _hash_xref_fname.empty()) {
        cerr << "Will write chemotype hash cross-reference to '" << _hash_xref_fname << "'\n";
      }
    }
    if (_chemotype_options.include_tied_adjacent_ring_systems) {
      cerr << "Will include adjacent ring systems tied at the -n cutoff distance\n";
    }
    if (_chemotype_options.prefer_larger_adjacent_ring_system) {
      cerr << "Will prefer larger adjacent ring systems among same-distance candidates\n";
    }
    if (_chemotype_options.choose_first_embedding) {
      cerr << "Will use the first query embedding if multiple ring systems match\n";
    }
    if (_atom_typing_ptr != nullptr) {
      cerr << "Will label non-terminal attachment atoms with atom types\n";
    }
    if (_chemotype_options.isotope_for_exit_points != 0) {
      cerr << "Will label chemotype boundary atoms with isotope "
           << _chemotype_options.isotope_for_exit_points << '\n';
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(m);
  }

  return ! m.empty();
}

int
Options::WriteEmptyFingerprint(IWString_and_File_Descriptor& output) {
  output << _fingerprint.tag << ">\n";
  output.write_if_buffer_holds_more_than(4096);
  return output.good();
}



int
Options::StartFingerprintRecord(Molecule& m, IWString_and_File_Descriptor& output) {
  if (! _function_as_tdt_filter) {
    output << "$SMI<" << m.smiles() << ">\n";
  }

  return output.good();
}

int
Options::FinishFingerprintRecord(Molecule& m, IWString_and_File_Descriptor& output) {
  if (! _function_as_tdt_filter) {
    output << "PCN<" << m.name() << ">\n";
    output << "|\n";
  }

  ++_molecules_written;
  output.write_if_buffer_holds_more_than(4096);
  return output.good();
}

int
Options::WriteFingerprint(Molecule& m, const int* include_atom,
                          const uint32_t* atom_type,
                          IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  int included_atoms = 0;
  for (int i = 0; i < matoms; ++i) {
    if (include_atom[i]) {
      ++included_atoms;
    }
  }

  if (included_atoms == 0) {
    return WriteEmptyFingerprint(output);
  }

  if (_fingerprint.kind == FingerprintKind::kFixed) {
    IWMFingerprint fp;
    if (! fp.construct_fingerprint(m, atom_type, include_atom)) {
      cerr << "Cannot generate fixed fingerprint for '" << m.name() << "'\n";
      return 0;
    }

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << _fingerprint.tag << tmp << ">\n";
    output.write_if_buffer_holds_more_than(4096);
    return output.good();
  }

  if (_fingerprint.kind == FingerprintKind::kNonColliding) {
    iwecfp::Iwecfp generator;
    generator.set_max_radius(_fingerprint.ec_radius);

    Sparse_Fingerprint_Creator sfc;
    const iwecfp::FingerprintResult result =
        generator.Fingerprint(m, atom_type, include_atom, &sfc);
    if (result == iwecfp::FingerprintResult::kFatal) {
      cerr << "Cannot generate EC fingerprint for '" << m.name() << "'\n";
      return 0;
    }
    if (result == iwecfp::FingerprintResult::kNoStartAtoms) {
      return WriteEmptyFingerprint(output);
    }

    IWString tmp;
    sfc.daylight_ascii_form_with_counts_encoded(_fingerprint.tag, tmp);
    output << tmp << '\n';
    output.write_if_buffer_holds_more_than(4096);
    return output.good();
  }

  cerr << "Options::WriteFingerprint:no fingerprint type active\n";
  return 0;
}

int
Options::ProcessFingerprint(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  m.compute_aromaticity_if_needed();

  if (_min_rings > 0 && m.nrings() < _min_rings) {
    ++_molecules_below_min_rings;
    if (! StartFingerprintRecord(m, output) || ! WriteEmptyFingerprint(output)) {
      return 0;
    }
    return FinishFingerprintRecord(m, output);
  }

  chemotypes::ChemotypeQueryMatch match;
  const chemotypes::ChemotypeQueryMatchStatus status = chemotypes::FirstChemotypeQueryMatch(
      m, _queries, match, _chemotype_options.choose_first_embedding);

  switch (status) {
    case chemotypes::ChemotypeQueryMatchStatus::kMatched:
      break;

    case chemotypes::ChemotypeQueryMatchStatus::kNoQueryMatch:
      ++_molecules_not_matching_queries;
      if (_ignore_molecules_not_matching_query) {
        if (! StartFingerprintRecord(m, output) || ! WriteEmptyFingerprint(output)) {
          return 0;
        }
        return FinishFingerprintRecord(m, output);
      }
      cerr << "No query matched '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom:
      ++_molecules_matched_query_without_ring_atom;
      cerr << "Query matched no ring atom in '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAmbiguousQueryMatches:
      ++_molecules_with_ambiguous_query_matches;
      cerr << "Query matched multiple ring systems in '" << m.name()
           << "', use -z first to use the first embedding\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAtomTypingFailed:
      ++_atom_typing_failures;
      cerr << "Atom typing failed for '" << m.name() << "'\n";
      return 0;
  }

  std::vector<int> chemotype_mask = chemotypes::ChemotypeAtomMask(
      m, match, _chemotype_options, _scratch);

  const int matoms = m.natoms();
  std::vector<int> include_atom(matoms, 0);
  Set_of_Atoms seeds;
  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    const int selected = _fingerprint.invert
        ? chemotype_mask[atom] == chemotypes::kChemotypeNotKept
        : chemotype_mask[atom] != chemotypes::kChemotypeNotKept;
    if (selected) {
      include_atom[atom] = 1;
      seeds << atom;
    }
  }

  if (_fingerprint.expand > 0) {
    lillymol::SetAtomsWithinRadius(m, seeds, _fingerprint.expand, 1,
                                   include_atom.data());
  }

  std::vector<uint32_t> atom_type(matoms, 0);
  if (! _atom_typing_ptr->assign_atom_types(m, atom_type.data())) {
    ++_atom_typing_failures;
    cerr << "Atom typing failed for '" << m.name() << "'\n";
    return 0;
  }

  if (! StartFingerprintRecord(m, output) ||
      ! WriteFingerprint(m, include_atom.data(), atom_type.data(), output)) {
    return 0;
  }

  return FinishFingerprintRecord(m, output);
}

int
Options::Process(Molecule& m, Molecule_Output_Object& output) {
  ++_molecules_read;
  m.compute_aromaticity_if_needed();

  if (_min_rings > 0 && m.nrings() < _min_rings) {
    ++_molecules_below_min_rings;
    return 1;
  }

  const IWString parent_name = m.name();

  std::unique_ptr<Molecule> parent;
  if (_write_parent) {
    parent = std::make_unique<Molecule>(m);
  }

  Molecule complement;
  chemotypes::ChemotypeQueryMatch match;
  const chemotypes::ChemotypeQueryMatchStatus status = chemotypes::ReduceToChemotype(
      m, _queries, _chemotype_options, _scratch, match, _atom_typing_ptr,
      _write_complement ? &complement : nullptr);

  switch (status) {
    case chemotypes::ChemotypeQueryMatchStatus::kMatched: {
      if (_write_complement && ! complement.empty()) {
        if (! _complement_output.write(complement)) {
          return 0;
        }
        ++_complements_written;
      }
      if (parent) {
        if (! _parent_annotation.empty()) {
          parent->append_to_name(_parent_annotation);
        }
        if (! output.write(*parent)) {
          return 0;
        }
        ++_parent_molecules_written;
      }
      uint32_t chemotype_id = 0;
      if (! AccumulateChemotype(m, parent_name, chemotype_id)) {
        return 0;
      }
      if (_write_unique_chemotype_id) {
        IWString suffix;
        suffix << ' ' << _unique_chemotype_id_prefix << chemotype_id;
        m.append_to_name(suffix);
      }
      if (_hash_chemotypes) {
        IWString hash;
        if (! ChemotypeHash(m, parent_name, hash)) {
          return 0;
        }
        IWString suffix;
        suffix << ' ' << hash;
        m.append_to_name(suffix);
      }
      ++_molecules_written;
      return output.write(m);
    }

    case chemotypes::ChemotypeQueryMatchStatus::kNoQueryMatch:
      ++_molecules_not_matching_queries;
      if (_ignore_molecules_not_matching_query) {
        return 1;
      }
      cerr << "No query matched '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom:
      ++_molecules_matched_query_without_ring_atom;
      cerr << "Query matched no ring atom in '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAmbiguousQueryMatches:
      ++_molecules_with_ambiguous_query_matches;
      cerr << "Query matched multiple ring systems in '" << m.name()
           << "', use -z first to use the first embedding\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAtomTypingFailed:
      ++_atom_typing_failures;
      cerr << "Atom typing failed for '" << m.name() << "'\n";
      return 0;
  }

  return 0;
}

dicer_data::Isotope
DicerIsotope(const chemotypes::ChemotypeOptions& options,
             const Atom_Typing_Specification* atom_typing) {
  if (atom_typing != nullptr) {
    return dicer_data::ATYPE;
  }

  if (options.isotope_for_exit_points != 0) {
    return dicer_data::ATT;
  }

  return dicer_data::NONE;
}

int
Options::AccumulateChemotype(Molecule& m, const IWString& parent_name, uint32_t& chemotype_id) {
  if (_summary_fname.empty()) {
    return 1;
  }

  const IWString& usmi = m.unique_smiles();
  auto iter = _chemotype.find(usmi);
  if (iter != _chemotype.end()) {
    iter->second.set_n(iter->second.n() + 1);
    chemotype_id = iter->second.id();
    ++_duplicate_chemotypes;
    return 1;
  }

  chemotype_id = _chemotype.size();

  dicer_data::DicerFragment proto;
  proto.set_smi(usmi.data(), usmi.length());
  proto.set_par(parent_name.data(), parent_name.length());
  proto.set_nat(m.natoms());
  proto.set_n(1);
  if (_write_unique_chemotype_id) {
    proto.set_id(chemotype_id);
  }

  const dicer_data::Isotope iso = DicerIsotope(_chemotype_options, _atom_typing_ptr);
  if (iso != dicer_data::NONE) {
    proto.set_iso(iso);
  }

  _chemotype.emplace(usmi, std::move(proto));
  ++_unique_chemotypes;
  return 1;
}


int
Options::ChemotypeHash(Molecule& m, const IWString& parent_name, IWString& hash) {
  const IWString& usmi = m.unique_smiles();

  auto existing = _hash_from_usmi.find(usmi);
  if (existing != _hash_from_usmi.end()) {
    hash = existing->second;
    auto xref = _hash_xref.find(hash);
    if (xref != _hash_xref.end()) {
      ++xref->second.n;
    }
    return 1;
  }

  hash = ChemotypeHashString(usmi, _hash_width);
  auto collision = _hash_xref.find(hash);
  if (collision != _hash_xref.end()) {
    cerr << "Chemotype hash collision for '" << hash << "'\n";
    cerr << " existing " << collision->second.smiles << '\n';
    cerr << " new      " << usmi << '\n';
    return 0;
  }

  ChemotypeHashXref xref;
  xref.smiles = usmi;
  xref.parent = parent_name;
  xref.n = 1;
  _hash_xref.emplace(hash, std::move(xref));
  _hash_from_usmi.emplace(usmi, hash);
  _hash_order.push_back(hash);
  return 1;
}

int
Options::WriteSummary() {
  if (_summary_fname.empty()) {
    return 1;
  }

  IWString_and_File_Descriptor output;
  if (! output.open(_summary_fname.null_terminated_chars())) {
    cerr << "Options::WriteSummary:cannot open '" << _summary_fname << "'\n";
    return 0;
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const auto& [usmi, proto] : _chemotype) {
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteSummary:cannot write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();
  return 1;
}


int
Options::WriteHashXref() {
  if (_hash_xref_fname.empty()) {
    return 1;
  }

  IWString_and_File_Descriptor output;
  if (! output.open(_hash_xref_fname.null_terminated_chars())) {
    cerr << "Options::WriteHashXref:cannot open '" << _hash_xref_fname << "'\n";
    return 0;
  }

  output << "chemotype_hash smiles n parent\n";
  for (const IWString& hash : _hash_order) {
    auto iter = _hash_xref.find(hash);
    if (iter == _hash_xref.end()) {
      cerr << "Options::WriteHashXref:internal error, missing hash '" << hash << "'\n";
      return 0;
    }
    const ChemotypeHashXref& xref = iter->second;
    output << hash << ' ' << xref.smiles << ' ' << xref.n << ' ' << xref.parent << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();
  return output.good();
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (fingerprinting_active()) {
    output << "Wrote " << _molecules_written << " fingerprints\n";
  } else {
    output << "Wrote " << _molecules_written << " chemotypes\n";
  }
  if (_write_complement) {
    output << "Wrote " << _complements_written << " complements\n";
  }
  if (_parent_molecules_written) {
    output << "Wrote " << _parent_molecules_written << " parent molecules\n";
  }
  if (! _summary_fname.empty()) {
    output << "Accumulated " << _unique_chemotypes << " unique chemotypes and "
           << _duplicate_chemotypes << " duplicate chemotypes\n";
  }
  if (_molecules_not_matching_queries) {
    output << _molecules_not_matching_queries << " molecules did not match any query\n";
  }
  if (_molecules_below_min_rings) {
    output << _molecules_below_min_rings << " molecules had fewer than " << _min_rings
           << " rings\n";
  }
  if (_molecules_matched_query_without_ring_atom) {
    output << _molecules_matched_query_without_ring_atom
           << " molecules matched a query with no ring atom\n";
  }
  if (_molecules_with_ambiguous_query_matches) {
    output << _molecules_with_ambiguous_query_matches
           << " molecules matched queries in multiple ring systems\n";
  }
  if (_atom_typing_failures) {
    output << _atom_typing_failures << " atom typing failures\n";
  }

  return 1;
}

int
ChemotypeFingerprints(Options& options, data_source_and_type<Molecule>& input,
                      IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! options.ProcessFingerprint(*m, output)) {
      return 0;
    }
  }

  return output.good();
}

int
ChemotypeFingerprints(Options& options, const char* fname, FileType input_type,
                      IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ChemotypeFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ChemotypeFingerprints(options, input, output);
}

int
ProcessFilterRecord(Options& options, const const_IWSubstring& buffer,
                    IWString_and_File_Descriptor& output) {
  const_IWSubstring smiles(buffer);
  smiles.remove_up_to_first('<');
  smiles.chop();

  Molecule m;
  if (! m.build_from_smiles(smiles)) {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  if (! options.Preprocess(m)) {
    return options.WriteEmptyFingerprint(output);
  }

  return options.ProcessFingerprint(m, output);
}

int
ChemotypeFingerprintFilter(Options& options, iwstring_data_source& input,
                           IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (! buffer.starts_with("$SMI<")) {
      output.write_if_buffer_holds_more_than(4096);
      continue;
    }

    if (! ProcessFilterRecord(options, buffer, output)) {
      cerr << "Fatal error processing TDT line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return output.good();
}

int
ChemotypeFingerprintFilter(Options& options, const char* fname,
                           IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ChemotypeFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return ChemotypeFingerprintFilter(options, input, output);
}

int
Chemotypes(Options& options, data_source_and_type<Molecule>& input,
           Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! options.Process(*m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Chemotypes(Options& options, const char* fname, FileType input_type,
           Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Chemotypes:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return Chemotypes(options, input, output);
}

int
Chemotypes(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:g:clfi:o:S:F:U:H:q:s:n:r:D:P:I:C:up:xtz:J:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    Usage(1);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process element options\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (options.function_as_tdt_filter()) {
    // TDT input may not have a molecule-file suffix.
  } else if (cl.number_elements() == 1 && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (options.fingerprinting_active()) {
    IWString_and_File_Descriptor output(1);
    for (const char* fname : cl) {
      int ok;
      if (options.function_as_tdt_filter()) {
        ok = ChemotypeFingerprintFilter(options, fname, output);
      } else {
        ok = ChemotypeFingerprints(options, fname, input_type, output);
      }
      if (! ok) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
    output.flush();
    if (verbose) {
      options.Report(cerr);
    }
    return 0;
  }

  Molecule_Output_Object output;
  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    return 1;
  }

  if (cl.option_present('S')) {
    IWString stem = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, stem)) {
      cerr << "Cannot overwrite input file(s) with stem '" << stem << "'\n";
      return 1;
    }
    if (! output.new_stem(stem)) {
      cerr << "Cannot open output stem '" << stem << "'\n";
      return 1;
    }
  } else {
    output.new_stem("-");
  }

  for (const char* fname : cl) {
    if (! Chemotypes(options, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (! options.WriteSummary()) {
    return 1;
  }

  if (! options.WriteHashXref()) {
    return 1;
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace chemotypes_main

int
main(int argc, char** argv) {
  return chemotypes_main::Chemotypes(argc, argv);
}
