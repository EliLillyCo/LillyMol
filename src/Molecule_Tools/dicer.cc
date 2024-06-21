/*
   Recursive breaking of a molecule.
*/

#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define IWQSORT_FO_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IW_IMPLEMENTATIONS_EXPOSED 1

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwmisc/iw_time.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/dicer_lib.h"

using std::cerr;
using std::endl;

using dicer_lib::USPVPTR;
using dicer_lib::InitialAtomNumber;
using dicer_lib::is_amide;
using dicer_lib::is_cf3_like;

//#define SPECIAL_THING_FOR_PARENT_AND_FRAGMENTS
#ifdef SPECIAL_THING_FOR_PARENT_AND_FRAGMENTS

static Molecule current_parent_molecule;
static Molecule_Output_Object stream_for_parent_joinged_to_fragment;

#endif

static const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 1;   // always

static Element_Transformations element_transformations;

static int discard_chirality = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int add_user_specified_queries_to_default_rules = 0;

static int fragments_discarded_for_not_containing_user_specified_queries = 0;

static int include_amide_like_bonds_in_hard_coded_queries = 1;

static int break_ring_chain_bonds = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries_for_bonds_to_not_break;

static int nq = 0;

static resizable_array_p<Substructure_Hit_Statistics> user_specified_queries_that_fragments_must_contain;

static int fragments_must_contain_heteroatom = 0;

static int append_atom_count = 0;

static int max_recursion_depth = 10;

static uint64_t max_fragments_per_molecule = 0;

static int ignore_molecules_not_matching_any_queries = 0;

static int molecules_not_hitting_any_queries = 0;

// Global mapping from fragment unique smiles to fragment id.
static IW_STL_Hash_Map_uint smiles_to_id;
// Given a fragment id, what is the smiles.
static resizable_array_p<IWString> id_to_smiles;

static int add_new_strings_to_hash = 1;

static IW_STL_Hash_Map_uint molecules_containing;

// If we wish to use atom map numbers to keep track of atom mappings
// from parents to fragments, we cannot use unique smiles, since that
// suppresses atom map numbers.
// Note that the resulting smiles are NOT unique, so any global
// accumulation will be invalid.
static int fragments_indexed_by_unique_smiles = 1;

/*
 During writing the fragstat file, we want to know the identity of the
 molecule that gives rise to each fragment.
 */

static IW_STL_Hash_Map_String starting_parent;

static IW_STL_Hash_Map_int atoms_in_parent;

static int accumulate_starting_parent_information = 0;

static int write_fragstat_as_proto = 0;
// It can be convenient to be able to read the proto file
// as a smiles. But the leading token will need to be stripped
// before parsing as a textproto.
static int prepend_smiles_to_fragstat_proto = 0;

// If we are writing binary protos - recommended for larger
// datasets.
static int write_fragstat_as_binary_proto = 0;
static iw_tf_data_record::TFDataWriter tfdata_writer;

static int work_like_recap = 0;

static int write_diced_molecules = 1;

static Accumulator_Base<clock_t, unsigned long> time_acc;

static int collect_time_statistics = 0;

static IWString_and_File_Descriptor stream_for_post_breakage;

/*
   We can avoid re-slicing fragments if we keep track of which
   fragments have already been encountered in the current molecule.
   So, if record_presence_and_absence_only is set, once a given
   fragment is encountered a 2nd time, it will not be subdivided.

   Note that when we do this, any counting information will be
   truncated
 */

static int record_presence_and_absence_only = 0;

static isotope_t isotope_for_join_points = 0;

static int increment_isotope_for_join_points = 0;

static int apply_atom_type_isotopic_labels = 0;

static int isotopic_label_is_recursion_depth = 0;

static int apply_isotopes_to_complementary_fragments = 0;

#define ISO_CF_ATYPE -2

static int number_by_initial_atom_number = 0;

static int apply_environment_isotopic_labels = 0;

static int add_environment_atoms = 0;

static int apply_atomic_number_isotopic_labels = 0;

static int write_only_smallest_fragment = 0;

static Accumulator<double> whole_set;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fingerprint_tag;

static int function_as_filter = 0;

static extending_resizable_array<int> global_counter_bonds_to_be_broken;

static int lower_atom_count_cutoff = 1;
static int upper_atom_count_cutoff = std::numeric_limits<int>::max();

// When limiting the size of fragments produced, we might want to make
// sure we get ring systems, but having a fixed upper_atom_count_cutoff can
// make that difficult. Limiting the number of non ring atoms allows
// flexibility.
// A value of -1 means not specified. Note that a value of 0 is OK, and
// would result in only rings being produced.
static int max_non_ring_atoms = -1;

static int max_atoms_lost_from_parent = -1;

static float min_atom_fraction = -1.0f;
static float max_atom_fraction = 1.0f;

static int lower_atom_count_cutoff_current_molecule = 0;
static int upper_atom_count_cutoff_current_molecule = 0;

static int reject_breakages_that_result_in_fragments_too_small = 0;

static int write_smiles = 1;
static int write_parent_smiles = 1;

// We can write fragments in several forms.
//   1. classical form, smiles one per line
//   2. textproto, one per line
//   3. serialized proto form.
// This class holds that information and the associated streams.

class DicerFragmentOutput {
  private:
    // pointer because if we need to use stdout, we need to use the
    // constructor.
    std::unique_ptr<IWString_and_File_Descriptor> _text_output;

    // If true, we flush _text_output whenever we feel like it.
    int _buffered_output;

    iw_tf_data_record::TFDataWriter _as_tfdata;

    // True if a proto needs to be formed. We will write it as either
    // text or tfdata.
    int _output_is_proto;

    // private functions

    int WriteText(const dicer_data::DicedMolecule& proto, IWString_and_File_Descriptor& output);
    int WriteSerialized(const dicer_data::DicedMolecule& proto);

  public:
    DicerFragmentOutput();

    void set_output_is_proto(int s) {
      _output_is_proto = s;
    }
    int output_is_proto() const {
      return _output_is_proto;
    }

    void set_buffered_output(int s) {
      _buffered_output = s;
    }
    int buffered_output() const {
      return _buffered_output;
    }

    void flush() {
      if (_text_output) {
        _text_output->flush();
      }
    }

    // After processing a molecule do whatever flushing is needed.
    int MaybeFlush();

    int write_if_buffer_holds_more_than(uint32_t s) {
      return _text_output->write_if_buffer_holds_more_than(s);
    }

    IWString_and_File_Descriptor& text_output() {
      return *_text_output;
    }

    int SetUseStdout();

    int SetOutputStream(IWString& fname);

    int SetTFDataFname(IWString& fname);

    int Write(const dicer_data::DicedMolecule& proto);

    template <typename T>
    DicerFragmentOutput& operator << (const T rhs) {
      *_text_output << rhs;
      return *this;
    }
};

DicerFragmentOutput::DicerFragmentOutput() {
  _output_is_proto = 0;
  _buffered_output = 1;
}

int
DicerFragmentOutput::SetUseStdout() {
  _text_output = std::make_unique<IWString_and_File_Descriptor>(1);

  return _text_output->good();
}

int
DicerFragmentOutput::SetOutputStream(IWString& fname) {
  _text_output = std::make_unique<IWString_and_File_Descriptor>();
  if (! _text_output->open(fname.null_terminated_chars())) {
    cerr << "DicerFragmentOutput::SetOutputStream:cannot open '" << fname << "'\n";
    return 0;
  }

  return _text_output->good();
}

int
DicerFragmentOutput::SetTFDataFname(IWString& fname) {
  if (! _as_tfdata.Open(fname)) {
    cerr << "DicerFragmentOutput::SetTFDataFname:cannot open '" << fname << "'\n";
    return 0;
  }

  _output_is_proto = 1;

  return 1;
}

int
DicerFragmentOutput::Write(const dicer_data::DicedMolecule& proto) {
  if (_text_output) {
    return WriteText(proto, *_text_output);
  } else {
    return WriteSerialized(proto);
  }
}

int
DicerFragmentOutput::WriteText(const dicer_data::DicedMolecule& proto,
                IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (! printer.PrintToString(proto, &buffer)) {
    cerr << "DicerFragmentOutput::WriteTextcannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  output << buffer;
  output << '\n';

  if (_buffered_output) {
    output.write_if_buffer_holds_more_than(8192);
  } else {
    output.flush();
  }

  return 1;
}

int
DicerFragmentOutput::WriteSerialized(const dicer_data::DicedMolecule& proto) {
  return _as_tfdata.WriteSerializedProto<dicer_data::DicedMolecule>(proto);
}

int
DicerFragmentOutput::MaybeFlush() {
  if (! _text_output) {
    return 0;
  }

  if (_buffered_output) {
    _text_output->write_if_buffer_holds_more_than(8192);
  } else {
    _text_output->flush();
  }

  return 1;
}

/*
   When viewing fragments, it is convenient to have the parent
   structure on each page.
 */

static int vf_structures_per_page = 0;

static int write_smiles_and_complementary_smiles = 0;
// We can impose limits on the size of complementary fragments generated.
static int upper_atom_count_cutoff_complementary_fragment = 0;

static IWString_and_File_Descriptor bit_smiles_cross_reference;

static int allow_cc_bonds_to_break = 0;

static int use_terminal_atom_type = 0;
static int perceive_terminal_groups_only_in_original_molecule = 0;    // misguided, the terminal status is about the complementary fragment...

static int change_element_for_heteroatoms_with_hydrogens = 0;

static int eliminate_fragment_subset_relations = 0;

static const Element * holmium = nullptr;
static const Element * neptunium = nullptr;
static const Element * strontium = nullptr;

static int break_fused_rings = 0;

static int break_ring_bonds = 0;

static int break_spiro_rings = 0;

static int dearomatise_aromatic_rings = 0;

static int top_level_kekule_forms_switched = 0;

static extending_resizable_array<int> fragments_produced;

/*
   One mode of functioning is to just write the bonds to be broken
 */


static Molecule_Output_Object bbrk_file;
static IWString bbrk_tag;
static int read_bonds_to_be_broken_from_bbrk_file = 0;

static int check_for_lost_chirality = 0;

static Atom_Typing_Specification atom_typing_specification;

static Molecule_Output_Object stream_for_fully_broken_parents;

// Feb 2022
// When writing the file of all fragments discovered, we can control
/// what tokens are written to that file.
class AllFragmentsStream {
  // Output token separator. Could be made configurable.
  static constexpr char kSep = ' ';

  private:
    // Which columns/features to write.
    IW_STL_Hash_Set _fragment_directives_selected;

    // Destination for output.
    IWString_and_File_Descriptor _output;

    // Names of valid columns/features. Used during setup to identify
    // valid names.
    IW_STL_Hash_Set _valid_all_fragments_directives;

  // private functions
    int OpenFile(const const_IWSubstring& fname);
    int GetFeatureToWrite(const const_IWSubstring& token);

  public:
    AllFragmentsStream();

    int ParseDirective(const_IWSubstring b);  // By value.

    // If no features have been selected, use a default set.
    void UseDefaultFeaturesUnlessSpecified();

    int active() const {
      return _output.is_open();
    }

    // If feature `feature_name` has been selected for writing, write `feature_value`.
    // If `need_separator` is true, prepend output with kSep, and set `need_separator`
    // to true.
    template <typename T> int MaybeWrite(const char * feature_name,
                                         const T& feature_value, bool& need_separator);
    int newline() {
      _output << '\n';
      _output.write_if_buffer_holds_more_than(8192);
      return _output.good();
    }
};

AllFragmentsStream::AllFragmentsStream() {
  _valid_all_fragments_directives.insert("smiles");
  _valid_all_fragments_directives.insert("name");
  _valid_all_fragments_directives.insert("rdepth");
  _valid_all_fragments_directives.insert("psmiles");
  _valid_all_fragments_directives.insert("frag");
  _valid_all_fragments_directives.insert("frag_id");
}

int
AllFragmentsStream::ParseDirective(const_IWSubstring b) {
  if (b.starts_with("WRITEALL=")) {
    b.remove_leading_chars(9);
    return OpenFile(b);
  } else if (b.starts_with("WRITEALL:")) {
    b.remove_leading_chars(9);
    return GetFeatureToWrite(b);
  } else {
    return 1;
  }
}

int
AllFragmentsStream::OpenFile(const const_IWSubstring& fname) {
  if (_output.is_open()) {
    cerr << "AllFragmentsStream::OpenFile:already open\n";
    return 0;
  }

  IWString tmp(fname);
  if (! _output.open(tmp.null_terminated_chars())) {
    cerr << "AllFragmentsStream:OpenFile:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
AllFragmentsStream::GetFeatureToWrite(const const_IWSubstring& token) {
  int i = 0;
  IWString feature;
  while (token.nextword(feature, i, ',')) {
    if (! _valid_all_fragments_directives.contains(feature)) {
      cerr << "AllFragmentsStream::GetFeatureToWrite:unrecognised feature '" << feature << "'\n";
      return 0;
    }
    _fragment_directives_selected.insert(feature);
  }
  return 1;
}

void
AllFragmentsStream::UseDefaultFeaturesUnlessSpecified() {
  if (! _fragment_directives_selected.empty()) {
    return;
  }

  _fragment_directives_selected.insert("name");
  _fragment_directives_selected.insert("frag");
  _fragment_directives_selected.insert("frag_id");
}

template <typename T>
int
AllFragmentsStream::MaybeWrite(const char * feature_name,
                               const T & feature_value,
                               bool& need_separator) {
  if (! _fragment_directives_selected.contains(feature_name)) {
    return 1;
  }

  if (need_separator) {
    _output << kSep;
  } else {
    need_separator = true;
  }
  _output << feature_value;
  return _output.good();
}

static AllFragmentsStream all_fragments;

static void
reset_variables()
{
     verbose = 0;
     molecules_read = 0;
     reduce_to_largest_fragment = 1;   // always
     discard_chirality = 0;
     add_user_specified_queries_to_default_rules = 0;
     fragments_discarded_for_not_containing_user_specified_queries = 0;
     include_amide_like_bonds_in_hard_coded_queries = 1;
     break_ring_bonds = 0;
     nq = 0;
     append_atom_count = 0;
     max_recursion_depth = 1;
     max_fragments_per_molecule = 0;
     ignore_molecules_not_matching_any_queries = 0;
     molecules_not_hitting_any_queries = 0;
     add_new_strings_to_hash = 1;
     accumulate_starting_parent_information = 0;
     work_like_recap = 0;
     collect_time_statistics = 0;
     record_presence_and_absence_only = 0;
     isotope_for_join_points = 0;
     increment_isotope_for_join_points = 0;
     apply_atom_type_isotopic_labels = 0;
     isotopic_label_is_recursion_depth = 0;
     apply_isotopes_to_complementary_fragments = 0;
     number_by_initial_atom_number = 0;
     apply_environment_isotopic_labels = 0;
     add_environment_atoms = 0;
     apply_atomic_number_isotopic_labels = 0;
     smiles_tag = "$SMI<";
     identifier_tag = "PCN<";
     function_as_filter = 0;
     lower_atom_count_cutoff = 1;
     upper_atom_count_cutoff = std::numeric_limits<int>::max();
     max_non_ring_atoms = 0;
     max_atoms_lost_from_parent = -1;
     min_atom_fraction = -1.0f;
     lower_atom_count_cutoff_current_molecule = 0;
     upper_atom_count_cutoff_current_molecule = 0;
     reject_breakages_that_result_in_fragments_too_small = 0;
     write_smiles = 1;
     write_parent_smiles = 1;
     vf_structures_per_page = 0;
     write_smiles_and_complementary_smiles = 0;
     allow_cc_bonds_to_break = 0;
     use_terminal_atom_type = 0;
     perceive_terminal_groups_only_in_original_molecule = 0;
     change_element_for_heteroatoms_with_hydrogens = 0;
     eliminate_fragment_subset_relations = 0;
     holmium = nullptr;
     neptunium = nullptr;
     strontium = nullptr;
     break_fused_rings = 0;
     break_ring_bonds = 0;
     break_spiro_rings = 0;
     dearomatise_aromatic_rings = 0;
     top_level_kekule_forms_switched = 0;
     read_bonds_to_be_broken_from_bbrk_file = 0;
	 bbrk_file.resize(0);
	 bbrk_tag.resize(0);
     check_for_lost_chirality = 0;
     return;
}

/*
   We can have the support level for reporting a fragment
   specified as either a fraction of molecules hit, or as
   a minimum number of molecules that must contain the
   fragment. This object just encapsulates that
 */


class Fragment_Statistics_Report_Support_Specification
{
  private:
    float _fraction;
    int _int;

//  We work just fine if not initialised, but sometimes we need
//  to know if the object has in fact been initialised

    int _initialised;

  public:
    Fragment_Statistics_Report_Support_Specification();

    int debug_print (std::ostream & os) const;

    int initialised () const
    {
      return _initialised;
    }

    int build (const const_IWSubstring &);

    int fractional_support_level_specified () const
    {
      return _fraction >= 0.0F;
    }

    int set_population_size (int);

    int meets_support_requirement (int c) const
    {
      return c >= _int;
    }
};

Fragment_Statistics_Report_Support_Specification::Fragment_Statistics_Report_Support_Specification()
{
  _fraction = -1.0F;
  _int = 0;

  _initialised = 0;

  return;
}

int
Fragment_Statistics_Report_Support_Specification::debug_print (std::ostream & os) const
{
  if (_int > 0)
    os << "Support level is " << _int << " items\n";
  else if (_fraction >= 0.0F)
    os << "Support level is fractional " << _fraction << '\n';
  else
    os << "Support level not specified\n";

  return 1;
}

int
Fragment_Statistics_Report_Support_Specification::build (const const_IWSubstring & s)
{
  _initialised = 1;

  if (s.numeric_value(_int))
  {
    if (_int <= 0)
    {
      cerr << "Fragment_Statistics_Report_Support_Specification::build:the support level must be a whole +ve number\n";
      return 0;
    }

    return 1;
  }

  if (s.numeric_value(_fraction))
  {
    if (_fraction < 0.0F || _fraction > 1.0F)
    {
      cerr << "Fragment_Statistics_Report_Support_Specification::build:the fractional support level must be a valid fraction\n";
      return 0;
    }

    return 1;
  }

  cerr << "Fragment_Statistics_Report_Support_Specification::build:the support level specification must be a fraction or a +ve integer\n";

  _initialised = 0;

  return 0;
}

int
Fragment_Statistics_Report_Support_Specification::set_population_size (int p)
{
  if (_int > 0)          // if specified as an int, the population size does not matter
    return _int;

  _int = static_cast<int>(static_cast<float>(p) * _fraction + 0.4999f);

  return _int;
}

static int
parse_bond_specification_token (const const_IWSubstring &token,
                                int matoms,
                                atom_number_t & a)
{
  int j;

  if (!token.numeric_value(j) || j < 1 || j > matoms)
  {
    cerr << "Invalid atom specifier '" << token << "'\n";
    return 0;
  }

  a = j - 1;

  return 1;
}
/*
   When we trim the number of bonds to be broken, or reduce the recursion
   depth in order to accommodate the max_fragments_per_molecule constraint,
   it is handy to be able to figure out ahead of time what the number
   of fragments generated will be.

   After much scribbling on paper, I determined that the pattern looks like:

   1:  1 + 1
   2:  1 + 2 +  1
   3:  1 + 3 +  3 +  1
   4:  1 + 4 +  6 +  4 + 1
   5:  1 + 5 + 10 + 10 + 5 +  1
   6:  1 + 6 + 15 + 20 + 15 + 6 +  1
   7:  1 + 7 + 21 + 35 + 35 + 21 + 7 + 1

   where each row corresponds to the given number of bonds b:, and each
   column is a recursion level.  First column is always 1 which is the
   case for zero bonds being broken.  If just one bond can be broken,
   then the number of fragments produced will be the number of bonds.
   To get the total number of fragments produced at a given recursion
   level, sum the columns out to that value to get the number of
   fragments produced

 */

class Fragments_At_Recursion_Depth
{
  private:
    int _max_recursion;           // from the command line
    int _max_bonds;               // just some convenient number to keep the array sizes down. Suggest 20

    uint64_t * _fragments;

  public:
    Fragments_At_Recursion_Depth();
    ~Fragments_At_Recursion_Depth();

    int initialise (int k, int b);

    int debug_print (std::ostream &) const;

    uint64_t fragments_at_recursion_depth (int k, int b) const;
};

Fragments_At_Recursion_Depth::Fragments_At_Recursion_Depth()
{
  _max_recursion = 0;
  _max_bonds = 0;

  _fragments = nullptr;
}

/*
   The number of fragments produced at any given combination of
   k and n is related to (n choose k)
 */

int
Fragments_At_Recursion_Depth::initialise (int k, int b)
{
  _max_recursion = k;
  _max_bonds = b;

  const uint64_t k1 = k + 1;
  const uint64_t b1 = b + 1;

  _fragments = new uint64_t[k1 * b1];

  set_vector(_fragments, k1 * b1, static_cast<uint64_t>(0));

  _fragments[0] = 1;         // for the recurrence relation, need to initialise something

  for (uint64_t i = 1; i < b1; i++)         // loop over number of bonds
  {
    const uint64_t * fprev = _fragments + (i - 1) * k1;

    uint64_t * f = _fragments + i * k1;

    f[0] = 1;
    for (uint64_t j = 1; j < k1; j++)
    {
      f[j] = fprev[j - 1] + fprev[j];
      if (0 == f[j])
        break;
    }

    assert (i == f[1]);                // of one bond broken, can produce i fragments
  }

  // Now sum them

  for (uint64_t i = 1; i < b1; i++)
  {
    uint64_t * f = _fragments + i * k1;
    for (uint64_t j = 1; j < k1; j++)
    {
      f[j] += f[j - 1];
    }
  }

  return 1;
}

int
Fragments_At_Recursion_Depth::debug_print (std::ostream & os) const
{
  os << "Fragments_At_Recursion_Depth:_max_recursion " << _max_recursion << " _max_bonds " << _max_bonds << '\n';
  for (int b = 1; b <= _max_bonds; b++)
  {
    for (int k = 0; k <= _max_recursion; k++)
    {
      int f = fragments_at_recursion_depth(k, b);

      os << "For " << b << " bonds and k = " << k << " will get " << f << '\n';
      if (0 == f)
        break;
    }
  }

  os << '\n';

  return 1;
}

Fragments_At_Recursion_Depth::~Fragments_At_Recursion_Depth()
{
  if (nullptr != _fragments)
    delete [] _fragments;

  return;
}

uint64_t
Fragments_At_Recursion_Depth::fragments_at_recursion_depth (int k, int b) const
{
  int multiplier;

  if (b > _max_bonds)          // too horrible to think about
    multiplier = 11 * (b - _max_bonds);                // just some random thing
  else
    multiplier = 1;

  const uint64_t * f = _fragments + b * (_max_recursion + 1);

  return multiplier * f[k];
}

static Fragments_At_Recursion_Depth fard;

/*
   We need a means of passing around arguments
 */

class Dicer_Arguments
{
  private:
    Molecule * _current_molecule;            // not const because we might ask for a smiles

    int _max_recursion_depth_this_molecule;

    int _recursion_depth;

//  As each fragment is discovered, we assign it a numeric identifier
//  and store the count in a hash. Found that if I tried a string
//  hash on the smiles, things ran 10% slower. This seems to be the
//  fastest way

    typedef IW_Hash_Map<unsigned int, unsigned int> ff_map;

    ff_map _fragments_found_this_molecule;

//  We need a fragment membership array for every recursion depth.
//  We just allocate a single array and hand out slices of it

    int _atoms_in_molecule;

//  We need to quickly be able to figure out if an atom in the
//  parent is aromatic

    int * _aromatic_in_parent;

//  We want a convenient nbits that holds all the atoms in the molecule

    int _nbits;

//  For each subset found, we generate a bit vector based on the atoms in
//  the subset. Quick way for determining uniqueness.

    resizable_array_p<IW_Bits_Base> _stored_bit_vectors;

//  At every step of the recursion, we need to keep track of
//  the correspondence between the atom numbers in the fragment
//  molecules, and the atom number of that atom in the initial molecule

    int * _atom_number_in_initial_molecule;

//  The xref array is used by create_subset()

    int * _xref;

//  When translating atom numbers from parent to subset, we need a temporary
//  array. It could be allocated each invocation, but we make it a class
//  member just for efficiency

    int * _local_xref;
    int * _local_xref2;

//  For each pair of atoms (bond) we keep track of the largest example of
//  that bond from the database

//  When doing ring splitting, we can get huge numbers of fragments,
//  so we keep track of the number of times we have attempt it

    int _ring_splitting_attempted;

//  When we have an atom typing active, we record the types assigned

    uint32_t * _atom_type;

    uint32_t * _isosave;           // for quickly saving and restoring isotopes - not allocated separately

    int _lower_atom_count_cutoff_current_molecule;
    int _upper_atom_count_cutoff_current_molecule;
    int _max_non_ring_atoms;

    // If the task is to just generate a cross reference file, we may
    // not need the usual output.
    int _write_diced_molecules;


// private functions

    int _write_fragment_and_complement (Molecule & m, const IW_Bits_Base & b,
                      int * subset_specification, IWString_and_File_Descriptor & output) const;
    int _write_fragment_and_hydrogen (const Molecule & m, int * subset_specification,
                      IWString_and_File_Descriptor & output) const;
    int _write_fragment_and_hydrogen (Molecule & m, int * subset_specification,
                      IWString_and_File_Descriptor & output) const;
    int _do_apply_isotopes_to_complementary_fragments(Molecule & parent,
                                                      Molecule & subset, const int xref[]) const;
    int _do_apply_isotopes_to_complementary_fragments(Molecule & parent,
                                                      Molecule & subset, Molecule & comp) const;

    int _do_apply_isotopes_to_complementary_fragments_canonical(Molecule & parent, int n,
                                                  Molecule & fragment, Molecule& comp) const;
    int _is_smallest_fragment(const IW_Bits_Base & b) const;
    void _add_new_smiles_to_global_hashes (const IWString & smiles);

    int WriteAsProto(Molecule& n, int breakable_bonds, DicerFragmentOutput& output) const;

    int ToProto(Molecule& m, int breakable_bonds,
                         dicer_data::DicedMolecule& proto) const;
    int AppendFragmentAndComplement(Molecule& m,
                                  const IW_Bits_Base & b,
                                  int * subset_specification,
                                  dicer_data::DicedMolecule& proto) const;
    int WriteFragmentsAndComplementsProto(Molecule& m,
           DicerFragmentOutput& output) const;

    int MaybeAddToHash(const IWString& smiles);

  public:
    Dicer_Arguments (int mrd);            // arg is max recursion depth
    ~Dicer_Arguments ();

    void set_current_molecule(Molecule * m)
    {
      _current_molecule = m;
    }

    int atoms_in_molecule () const
    {
      return _atoms_in_molecule;
    }

    int set_array_size (int);

    int initialise_parent_molecule_aromaticity (Molecule & m);

    int is_aromatic_in_parent (atom_number_t a) const
    {
      return _aromatic_in_parent[a];
    }
    int is_aromatic_in_parent (atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4) const
    {
      return _current_molecule->bond_between_atoms(a1, a2)->is_aromatic()
             && _current_molecule->bond_between_atoms(a2, a3)->is_aromatic()
             && _current_molecule->bond_between_atoms(a3, a4)->is_aromatic();
    }

    int  lower_atom_count_cutoff () const
    {
      return _lower_atom_count_cutoff_current_molecule;
    }
    void set_lower_atom_count_cutoff (int s)
    {
      _lower_atom_count_cutoff_current_molecule = s;
    }

    int  upper_atom_count_cutoff () const
    {
      return _upper_atom_count_cutoff_current_molecule;
    }
    void set_upper_atom_count_cutoff (int s)
    {
      _upper_atom_count_cutoff_current_molecule = s;
    }

    int ok_atom_count (int) const;
    int OkAtomCount(Molecule& m) const;

    void set_write_diced_molecules(int s) {
      _write_diced_molecules = s;
    }

    int max_recursion_depth_this_molecule () const
    {
      return _max_recursion_depth_this_molecule;
    }
//  void set_max_recursion_depth_this_molecule(int s) { _max_recursion_depth_this_molecule = s;}

    void increment_recursion_depth()
    {
      _recursion_depth++;
    }
    void decrement_recursion_depth()
    {
      _recursion_depth--;
    }
    int recursion_depth() const
    {
      return _recursion_depth;
    }

    uint64_t adjust_max_recursion_depth (int bonds_to_break_this_molecule, uint64_t max_fragments_per_molecule);

    int can_go_deeper() const
    {
      return _recursion_depth < _max_recursion_depth_this_molecule;
    }

    const int * parent_atom_number_in_initial_molecule() const
    {
      return _atom_number_in_initial_molecule + (_recursion_depth - 1) * _atoms_in_molecule;
    }

    int * atom_number_in_initial_molecule()
    {
      return _atom_number_in_initial_molecule + _recursion_depth * _atoms_in_molecule;
    }

    int * xref_array_this_level()
    {
      return _xref + _recursion_depth * _atoms_in_molecule;
    }

    int update_numbering_for_subset(const int * xref,
      const Set_of_Atoms & a0,
      const Set_of_Atoms & a1,
      Set_of_Atoms & b0,
      Set_of_Atoms & b1);

//  int number_fragments() const { return _fragments_found_this_molecule.size();}

    int is_unique(int atoms_in_parent, const int * xref, Molecule & c);
    int is_unique (Molecule & m);

    int contains_fragment(Molecule & smiles);

//  const resizable_array<int> & fragments_found_this_molecule() const { return _fragments_found_this_molecule;}

//  void sort_fragments_found_this_molecule() { _fragments_found_this_molecule.sort(int_comparitor_larger);}

    int try_more_ring_splitting ();

    int write_fragments_found_this_molecule(Molecule & m, int breakable_bonds, DicerFragmentOutput & output) const;

    int produce_fingerprint(Molecule & m, IWString_and_File_Descriptor & output) const;

    int write_fragments_and_complements (Molecule & m, DicerFragmentOutput & output) const;

    int store_atom_types (Molecule & m,  Atom_Typing_Specification & ats);

    const uint32_t * atom_types() const { return _atom_type;}

    int fragments_found () const
    {
      return _fragments_found_this_molecule.size();
    }

    // If `all_fragments_stream` is active, write information to it.
    int WriteAllFragments(Molecule& parent,
                          Molecule& fragment,
                          AllFragmentsStream& all_fragments_stream) const;
};

Dicer_Arguments::Dicer_Arguments(int mrd)
{
  _current_molecule = nullptr;

  _aromatic_in_parent = nullptr;

  _recursion_depth = 0;

  _max_recursion_depth_this_molecule = mrd;

  _atom_number_in_initial_molecule = nullptr;

  _xref = nullptr;

  _local_xref = nullptr;
  _local_xref2 = nullptr;

  _ring_splitting_attempted = 0;

  _atom_type = nullptr;

  _write_diced_molecules = 1;

  // Copy from file scope value.
  _max_non_ring_atoms = max_non_ring_atoms;

  return;
}

Dicer_Arguments::~Dicer_Arguments()
{
  if (nullptr != _aromatic_in_parent)
    delete [] _aromatic_in_parent;

  if (nullptr != _xref)
    delete [] _xref;

  if (nullptr != _local_xref)
    delete [] _local_xref;

  delete []_local_xref2;

  if (nullptr != _atom_number_in_initial_molecule)
    delete [] _atom_number_in_initial_molecule;

  if (nullptr != _atom_type)
    delete [] _atom_type;

  return;
}

int
Dicer_Arguments::set_array_size(int s)
{
  assert (s > 1);

  _atoms_in_molecule = s;

  if (0 == _atoms_in_molecule % IW_BITS_PER_WORD)
    _nbits = _atoms_in_molecule;
  else
    _nbits = (_atoms_in_molecule / IW_BITS_PER_WORD + 1) * IW_BITS_PER_WORD;

  // In many of the arrays, level 0 never gets used, so the arrays need to
  // be one larger than the max recursion depth

  int array_size;
  if (_max_recursion_depth_this_molecule > 0)
    array_size = (_max_recursion_depth_this_molecule + 1) * s;
  else if (max_recursion_depth > 0)
    array_size = (max_recursion_depth + 1) * s;
  else
    array_size = s + s;

  //cerr << "array_size " << array_size << '\n';

  _atom_number_in_initial_molecule = new int [array_size];

  // Set the first chunk of the atom number array

  for (int i = 0; i < s; i++)
  {
    _atom_number_in_initial_molecule[i] = i;
  }

  _xref = new_int(array_size, -1);

  _local_xref = new int[_atoms_in_molecule];
  _local_xref2 = new int[_atoms_in_molecule];

  return (nullptr != _xref && nullptr != _local_xref);
}

int
Dicer_Arguments::initialise_parent_molecule_aromaticity (Molecule & m)
{
  _aromatic_in_parent = new int[_atoms_in_molecule];

  for (int i = 0; i <  _atoms_in_molecule; i++)
  {
    if (m.is_aromatic(i))
      _aromatic_in_parent[i] = 1;
    else
      _aromatic_in_parent[i] = 0;
  }

  return 1;
}

// Fill in a cross reference array that converts from atom numbers n the parent
// to atom numbers in the current molecule

static void
initialise_xref (const Molecule & m,
                 const Dicer_Arguments & dicer_args,
                 int * xref)
{
  std::fill_n(xref, dicer_args.atoms_in_molecule(), -1);

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    std::optional<atom_number_t> x = InitialAtomNumber(m.atomi(i));

    xref[*x] = i;
  }

  return;
}

// Gets used by all transformations to convert atom numbers in the
// starting structure to atom numbers in the current structure

static int
convert_to_atom_numbers_in_child (const Set_of_Atoms & rbb,
                                  const int * xref,
                                  Set_of_Atoms & s)
{
  const int n = rbb.number_elements();

  s.resize_keep_storage(0);

  for (int i = 0; i < n; i++)
  {
    const atom_number_t jparent = rbb[i];

    const atom_number_t j = xref[jparent];
    if (j < 0)
      return 0;

    s.add(j);
  }

  return 1;
}

/*
   All dicer transformations are built on this class
 */

#define TRANSFORMATION_TYPE_BREAK_CHAIN_BOND 1
#define TRANSFORMATION_TYPE_BREAK_RING_BOND 2
#define TRANSFORMATION_TYPE_FUSED_RING_BREAKAGE 4
#define TRANSFORMATION_TYPE_FUSED_DEAROMATISE_RING 8
#define TRANSFORMATION_TYPE_SPIRO_RING_BREAKAGE 16

class Breakages;
class Breakages_Iterator;

class Dicer_Transformation
{
  protected:
    int _ttype;           // type of transformations, derived classes will fill in

    virtual int _write_bonds_to_be_broken (std::ostream &) const;

  public:
    Dicer_Transformation();

    virtual ~Dicer_Transformation();

    int process (Molecule & m, Dicer_Arguments & dicer_args, Breakages & b, Breakages_Iterator & bi) const
    {
      return _process(m, dicer_args, b, bi);
    }

    int break_bonds (Molecule & m) const
    {
      return _break_bonds(m);
    }

    virtual int _debug_print (std::ostream &) const = 0;
    virtual int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator & bi) const = 0;
    virtual int _break_bonds (Molecule &) const = 0;

    int write_bonds_to_be_broken (std::ostream & os) const
    {
      return _write_bonds_to_be_broken (os);
    }

    int debug_print (std::ostream & os) const
    {
      return _debug_print (os);
    }

    int ttype () const
    {
      return _ttype;
    }
};

Dicer_Transformation::Dicer_Transformation()
{
  _ttype = 0;
}

Dicer_Transformation::~Dicer_Transformation()
{
  return;
}

int
Dicer_Transformation::_write_bonds_to_be_broken (std::ostream & os) const
{
  return 1;
}

#ifdef __GNUG__
template class resizable_array_p<Dicer_Transformation>;
template class resizable_array_base<Dicer_Transformation *>;
#endif

class Bond_Breakage : public Dicer_Transformation
{
  protected:
    const atom_number_t _a1;
    const atom_number_t _a2;
    int _symmetry_equivalent;

  public:
    Bond_Breakage(atom_number_t z1, atom_number_t z2) : _a1(z1), _a2(z2)
    {
      _symmetry_equivalent = 0;
    }
    ~Bond_Breakage()  {
    };

    atom_number_t a1 () const
    {
      return _a1;
    }
    atom_number_t a2 () const
    {
      return _a2;
    }

    int involves (atom_number_t, atom_number_t) const;

    int symmetry_equivalent () const
    {
      return _symmetry_equivalent;
    }
    void set_symmetry_equivalent (int s)
    {
      _symmetry_equivalent = s;
    }
};

int
Bond_Breakage::involves (atom_number_t z1,
                         atom_number_t z2) const
{
  if (_a1 == z1 && _a2 == z2)
    return 1;

  if (_a1 == z2 && _a2 == z1)
    return 1;

  return 0;
}

class Chain_Bond_Breakage : public Bond_Breakage
{
  private:

    int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator &) const;
    int _debug_print (std::ostream &) const;
    int _write_bonds_to_be_broken (std::ostream &) const;
    int _process(Molecule& m0, resizable_array<Molecule *> & components, const IWString & mname,
                               Dicer_Arguments & dicer_args, Breakages & breakages,
                               Breakages_Iterator & bi) const;
    int _break_bonds(Molecule & m) const;

  public:
    Chain_Bond_Breakage(atom_number_t, atom_number_t);
};

Chain_Bond_Breakage::Chain_Bond_Breakage(atom_number_t j1, atom_number_t j2) : Bond_Breakage(j1, j2)
{
  _ttype = TRANSFORMATION_TYPE_BREAK_CHAIN_BOND;

  return;
}

int
Chain_Bond_Breakage::_debug_print (std::ostream & os) const
{
  os << "Chain_Bond_Breakage: atoms " << _a1 << " and " << _a2 << '\n';

  return 1;
}

int
Chain_Bond_Breakage::_write_bonds_to_be_broken (std::ostream & os) const
{
  os << (_a1 + 1) << ' ' << (_a2 + 1) << '\n';

  return 1;
}

int
Chain_Bond_Breakage::_break_bonds (Molecule & m) const
{
  m.remove_bond_between_atoms(_a1, _a2);

  return 1;
}

class Fused_Ring_Breakage : public Bond_Breakage
{
  private:
    const atom_number_t _a11;
    const atom_number_t _a12;
    const atom_number_t _a21;
    const atom_number_t _a22;

    int _aromatic_in_parent;

//  private functions

    int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator &) const;
    int _debug_print (std::ostream &) const;
    int _break_bonds(Molecule & m) const;

  public:
    Fused_Ring_Breakage (atom_number_t, atom_number_t, atom_number_t, atom_number_t, atom_number_t, atom_number_t);
    ~Fused_Ring_Breakage()
    {
      return;
    }

    atom_number_t a11 () const
    {
      return _a11;
    }
    atom_number_t a12 () const
    {
      return _a12;
    }
    atom_number_t a21 () const
    {
      return _a21;
    }
    atom_number_t a22 () const
    {
      return _a22;
    }

    int  aromatic_in_parent () const
    {
      return _aromatic_in_parent;
    }
    void set_aromatic_in_parent (int s)
    {
      _aromatic_in_parent = s;
    }
};

template class resizable_array_p<Fused_Ring_Breakage>;
template class resizable_array_base<Fused_Ring_Breakage *>;

Fused_Ring_Breakage::Fused_Ring_Breakage (atom_number_t z11, atom_number_t z1, atom_number_t z12,
                                          atom_number_t z21, atom_number_t z2, atom_number_t z22) :
                                                Bond_Breakage(z1, z2),
                                                _a11(z11), _a12(z12), _a21(z21), _a22(z22)
{
  _ttype = TRANSFORMATION_TYPE_FUSED_RING_BREAKAGE;

  _aromatic_in_parent = 0;

  return;
}

int
Fused_Ring_Breakage::_debug_print (std::ostream & os) const
{
  os << "Fused_Ring_Breakage: atoms a11 " << _a11 << " a1 " << _a1 << " a12 " << _a12 << '\n';
  os << "                           a21 " << _a21 << " a2 " << _a2 << " a22 " << _a22 << '\n';

  return 1;
}

int
Fused_Ring_Breakage::_break_bonds (Molecule & m) const
{
  return 1;
}

// Some transforms must act symmetrically. Otherwise the output may be different
// according to different input smiles to the same structure.
// Such process includes -B BRB -B BRF -B BRS

class PairTransformation : public Dicer_Transformation
{
  public:
    PairTransformation(Dicer_Transformation* t1, Dicer_Transformation* t2);
    virtual ~PairTransformation();
    virtual int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator & bi) const;
    virtual int _debug_print (std::ostream &) const;
    int _break_bonds(Molecule & m) const;

  private:
    Dicer_Transformation* _t1;
    Dicer_Transformation* _t2;
};

PairTransformation::PairTransformation (Dicer_Transformation* t1, Dicer_Transformation* t2)
{
  _t1 = t1;
  _t2 = t2;
}

PairTransformation::~PairTransformation()
{
  delete _t1;
  delete _t2;
}


int PairTransformation::_debug_print (std::ostream & stream) const
{
  if (_t1)
  {
    _t1->_debug_print(stream);
  }
  if (_t2)
  {
    _t2->_debug_print(stream);
  }
  return 1;
}


int
PairTransformation::_break_bonds (Molecule & m) const
{
  return 1;
}

///////////////////////////////////////////////////////////////////////
// Added by madlee 2011-12-04
class Spiro_Ring_Breakage : public Dicer_Transformation
{
  private:
    const atom_number_t _aCenter;
    const atom_number_t _a11;
    const atom_number_t _a12;
    const atom_number_t _a21;
    const atom_number_t _a22;


    int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator &) const;
    int _debug_print (std::ostream &) const;
    int _break_bonds(Molecule & m) const;

  public:
    Spiro_Ring_Breakage (atom_number_t, atom_number_t, atom_number_t, atom_number_t, atom_number_t);

    atom_number_t aCenter () const
    {
      return _aCenter;
    }
    atom_number_t a11 () const
    {
      return _a11;
    }
    atom_number_t a12 () const
    {
      return _a12;
    }
    atom_number_t a21 () const
    {
      return _a21;
    }
    atom_number_t a22 () const
    {
      return _a22;
    }

    bool involves (atom_number_t, atom_number_t) const;
};

template class resizable_array_p<Spiro_Ring_Breakage>;
template class resizable_array_base<Spiro_Ring_Breakage *>;

Spiro_Ring_Breakage::Spiro_Ring_Breakage (atom_number_t zCenter, atom_number_t z11, atom_number_t z12,
                                          atom_number_t z21, atom_number_t z22) :
                                                  _aCenter(zCenter),
                                                  _a11(z11), _a12(z12), _a21(z21), _a22(z22)
{
  _ttype = TRANSFORMATION_TYPE_SPIRO_RING_BREAKAGE;
  return;
}

int
Spiro_Ring_Breakage::_debug_print (std::ostream & os) const
{
  os << "Spiro_Ring_Breakage: Center:   " << _aCenter << '\n'
     << "                     Neighbors:" << _a11 << ", " << _a12 << ", " << _a21 << ", " << _a22 << '\n';

  return 1;
}

bool
Spiro_Ring_Breakage::involves (atom_number_t a1, atom_number_t a2) const
{
  assert (a1 != a2);
  if (a1 == aCenter())
  {
    return a11() == a2 || a12() == a2;
  }
  else if (a2 == aCenter())
  {
    return a11() == a1 || a12() == a1;
  }
  return false;
}

int
Spiro_Ring_Breakage::_break_bonds (Molecule & m) const
{
  return 1;
}

// Added by madlee
///////////////////////////////////////////////////////////////////////


class Dearomatise : public Dicer_Transformation
{
  private:
    Set_of_Atoms _ring;

//  private functions

    int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator &) const;
    int _debug_print (std::ostream &) const;
    int _break_bonds(Molecule & m) const;

  public:
    Dearomatise(const Ring &);
    ~Dearomatise()
    {
      return;
    }

    void set_ring_atoms(const Set_of_Atoms & r)
    {
      _ring = r;
    }
};

Dearomatise::Dearomatise (const Ring & r) : _ring(r)
{
  _ttype = TRANSFORMATION_TYPE_FUSED_DEAROMATISE_RING;

  return;
}

int
Dearomatise::_debug_print (std::ostream & os) const
{
  os << "Dearomatise ring";
  for (int i = 0; i < _ring.number_elements(); i++)
  {
    os << ' ' << _ring[i];
  }

  os << '\n';

  return os.good();
}

int
Dearomatise::_break_bonds (Molecule & m) const
{
  return 1;
}

/*
   A ring bond breakage is specified by the identities of the
   atoms to remove, and the identity of an atom to keep
 */

class Ring_Bond_Breakage : public Dicer_Transformation, public Set_of_Atoms
{
  private:
    atom_number_t _atom_that_must_be_a_ring_atom;
    bool _is_arom;

//  private functions

    int _do_ring_breaking (Molecule & m,
                           const Set_of_Atoms & rbb,
                           Dicer_Arguments & dicer_args,
                           Breakages & breakages,
                           Breakages_Iterator & bi) const;

    int _process (Molecule &, Dicer_Arguments &, Breakages &, Breakages_Iterator &) const;
    int _debug_print (std::ostream &) const;
    int _break_bonds(Molecule & m) const;

  public:
    Ring_Bond_Breakage(const Set_of_Atoms & s, bool is_arom);
    ~Ring_Bond_Breakage() { };

    void set_atom_that_must_be_a_ring_atom (atom_number_t a)
    {
      _atom_that_must_be_a_ring_atom = a;
    }
};

Ring_Bond_Breakage::Ring_Bond_Breakage(const Set_of_Atoms & s, bool is_arom) : ::Set_of_Atoms(s)
{
  _is_arom = is_arom;
  _ttype = TRANSFORMATION_TYPE_BREAK_RING_BOND;
  _atom_that_must_be_a_ring_atom = INVALID_ATOM_NUMBER;

  return;
}

int
Ring_Bond_Breakage::_debug_print (std::ostream & os) const
{
  os << "Ring_Bond_Breakage:: atoms ";

  for (int i = 0; i < _number_elements; i++)
  {
    os << ' ' << _things[i];
  }

  os << '\n';

  return 1;
}

int
Ring_Bond_Breakage::_break_bonds (Molecule & m) const
{
  return 1;
}

template class resizable_array_p<Ring_Bond_Breakage>;
template class resizable_array_base<Ring_Bond_Breakage *>;

class Breakages
{
  private:
    resizable_array_p<Dicer_Transformation> _transformations;

//  private functions

    int _parse_bond_specification_record (const IWString & t, int matoms);
    void _remove_fused_ring_breakage_if_present (atom_number_t a1, atom_number_t a2);
    void _remove_spiro_ring_breakage_if_present (atom_number_t a1, atom_number_t a2);
    void _remove_breakage_if_present (atom_number_t a1, atom_number_t a2);
    int _destroys_too_much_aromaticity (Molecule & m,
    const int aromatic_atom_count, const Fused_Ring_Breakage & f) const;
    int _destroys_too_much_aromaticity (Molecule & m, const int aromatic_atom_count,
                atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4) const;

  public:
    int number_transformations() const
    {
      return _transformations.number_elements();
    }

    const Dicer_Transformation * transformation (int i) const
    {
      return _transformations[i];
    }

    int debug_print (std::ostream &) const;

    int get_bonds_to_break_from_text_info(const Molecule & m, int istart);
    int get_bonds_to_break_from_text_info(const Molecule & m);

    int identify_bonds_to_break (Molecule & m);
    int identify_bonds_to_break_hard_coded_rules (Molecule & m);

    int identify_bonds_to_not_break(Molecule & m);

    int identify_cross_ring_bonds_to_be_broken (Molecule & m);

    int identify_spiro_ring_bonds_to_be_broken (Molecule & m);

    int identify_ring_bonds_to_be_broken (Molecule & m);

    int identify_aromatic_rings_to_de_aromatise (Molecule & m);

    int identify_symmetry_exclusions(Molecule & m);

    int assign_bond_priorities (Molecule & m, int * priority) const;

    int discard_breakages_that_result_in_fragments_too_small (Molecule & m,
                                      const Dicer_Arguments & dicer_args);

    void lower_numbers_if_needed (Molecule & m, int & bonds_to_break_this_molecule,
                          int & max_recursion_depth_this_molecule);

    int write_bonds_to_be_broken (Molecule & m, Molecule_Output_Object & bbrk_file) const;

    int do_recap (Molecule & m, Dicer_Arguments & dicer_args) const;
};

int
Breakages::debug_print (std::ostream & os) const
{
  int n = _transformations.number_elements();

  os << "Breakages contain " << n << " transformations\n";

  if (0 == n)
    return 1;

  IW_Hash_Map<int, int> ttypes;

  for (int i = 0; i < n; i++)
  {
    const Dicer_Transformation * ti = _transformations[i];

    int t = ti->ttype();

    ttypes[t]++;

    ti->debug_print(os);
  }

  for (IW_Hash_Map<int, int>::const_iterator i = ttypes.begin(); i != ttypes.end(); ++i)
  {
    os << (*i).second << " of type " << (*i).first << "\n";
  }

  return 1;
}

int
Breakages::_parse_bond_specification_record (const IWString & t,
  int matoms)
{
  if (2 != t.nwords())
  {
    cerr << "Bond breaking specifications must have at least two tokens\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  t.nextword(token, i);

  atom_number_t a1 = INVALID_ATOM_NUMBER;

  if (!parse_bond_specification_token(token, matoms, a1))
    return 0;

  t.nextword(token, i);

  atom_number_t a2 = INVALID_ATOM_NUMBER;

  if (!parse_bond_specification_token(token, matoms, a2))
    return 0;

  _transformations.add(new Chain_Bond_Breakage(a1, a2));

  return 1;
}

int
Breakages::get_bonds_to_break_from_text_info(const Molecule & m,
  int istart)
{
  int rc = 0;

  int n = m.number_records_text_info();

  for (int i = istart; i < n; i++)
  {
    const IWString & t = m.text_info(i);

    if (0 == t.length())
      break;

    if (!_parse_bond_specification_record (t, m.natoms()))
    {
      cerr << "Invalid bond breaking specification '" << t << "'\n";
      return 0;
    }

    rc++;
  }

  return rc;
}

int
Breakages::get_bonds_to_break_from_text_info(const Molecule & m)
{
  int n = m.number_records_text_info();

  if (n < 2)
  {
    cerr << "No text info with molecule '" << m.name() << "'\n";
    return 0;
  }

  IWString tag;
  tag << ">  <" << bbrk_tag << ">";

  for (int i = 0; i < n; i++)
  {
    const IWString & t = m.text_info(i);

    if (!t.starts_with(tag))
      continue;

    return get_bonds_to_break_from_text_info (m, i + 1);
  }

  cerr << "Did not find '" << tag << "' in text info for '" << m.name() << "'\n";
  return 0;
}


int
Breakages::write_bonds_to_be_broken (Molecule & m,
                                     Molecule_Output_Object & bbrk_file) const
{
  MDL_File_Supporting_Material * mdlfs = global_default_MDL_File_Supporting_Material();
  mdlfs->set_write_mdl_dollars (0);
  bbrk_file.write(m);
  mdlfs->set_write_mdl_dollars (1);

  std::ostream & os = bbrk_file.stream_for_type(FILE_TYPE_SDF);

  os << ">  <" << bbrk_tag << ">\n";
  for (int i = 0; i < _transformations.number_elements(); i++)
  {
    const Dicer_Transformation * ti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != ti->ttype())
      continue;

    ti->write_bonds_to_be_broken (os);
  }

  os << '\n';
  os << "$$$$\n";

  os.flush();

  if(!os.good())
     cerr<<"Breakages::write_bonds_to_be_broken:bad file stream!\n";

  return os.good();
}

/*
   We need an iterator into a Breakages object
 */

class Breakages_Iterator
{
  private:
    const int _n;
    int _state;

  public:
    Breakages_Iterator(const Breakages &);
    Breakages_Iterator(const Breakages_Iterator &);

    const Dicer_Transformation * current_transformation (const Breakages &);

    void increment()
    {
      _state++;
    }
    void decrement()
    {
      _state--;
    }

    int  state () const
    {
      return _state;                 // just used for debugging
    }

    int exhausted() const
    {
      return _state >= _n;
    }
};

Breakages_Iterator::Breakages_Iterator(const Breakages & b) :
                                          _n(b.number_transformations())
{
  _state = 0;

  return;
}

Breakages_Iterator::Breakages_Iterator (const Breakages_Iterator & rhs) :
                                          _n (rhs._n)
{
  _state = rhs._state;

  return;
}

const Dicer_Transformation *
Breakages_Iterator::current_transformation (const Breakages & b)
{
  if (_state >= _n)
    return nullptr;

  const Dicer_Transformation * rc = b.transformation(_state);

  _state++;

  return rc;
}

/*
   We can pre-compute the number of fragments expected for a given number
   of bonds in the molecule

   Some cases are easy.
   If just one bond, there are two fragments that can be formed.
 */

static int * fragments_expected = nullptr;

static int
initialise_fragments_per_number_of_bonds_array()
{
  assert (max_recursion_depth > 0);

  fragments_expected = new int[max_recursion_depth + 1];

  fragments_expected[0] = 0;
  fragments_expected[1] = 2;
  if (max_recursion_depth >= 2)
    fragments_expected[2] = 5;

  for (int i = 3; i <= max_recursion_depth; i++)
  {
    fragments_expected[i] = 2 * i;
    for (int j = 1; j < i; j++)
    {
      fragments_expected[i] += 2 * fragments_expected[j];
    }
  }

  if (verbose > 1)
  {
    for (int i = 0; i <= max_recursion_depth; i++)
    {
      cerr << "Max fragments for " << i << " bonds " << fragments_expected[i] << '\n';
    }
  }

  return 1;
}

int
Dicer_Arguments::try_more_ring_splitting ()
{
  _ring_splitting_attempted++;

  if (max_fragments_per_molecule > 0 && static_cast<uint64_t>(_ring_splitting_attempted) > max_fragments_per_molecule)
    return 0;

  return 1;
}

int
Dicer_Arguments::ok_atom_count (int matoms) const
{
  if (matoms < _lower_atom_count_cutoff_current_molecule)
    return 0;

  if (matoms > _upper_atom_count_cutoff_current_molecule)
    return 0;

  return 1;
}

int
Dicer_Arguments::OkAtomCount(Molecule& m) const {
  const int matoms = m.natoms();
  if (matoms < _lower_atom_count_cutoff_current_molecule) {
    return 0;
  }

  if (matoms > _upper_atom_count_cutoff_current_molecule) {
    return 0;
  }

  if (_max_non_ring_atoms < 0) {
    return 1;
  }

  int non_ring_atoms = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) > 0) {
      continue;
    }

    ++non_ring_atoms;
    if (non_ring_atoms > _max_non_ring_atoms) {
      return 0;
    }
  }

  return 1;
}

/*
   For a given number of bonds (B), and simultaneous breaks (K), generate
   the number of fragments produced at each -k value
*/

uint64_t
Dicer_Arguments::adjust_max_recursion_depth(int bonds_to_break_this_molecule,
                                            uint64_t max_fragments_per_molecule)
{
#ifdef NO_LONGER_USEDQWEUQYWE
  int m = max_recursion_depth;

  if (bonds_to_break_this_molecule < max_recursion_depth)
    m = bonds_to_break_this_molecule;
#endif

  auto fragments_that_will_be_produced = fard.fragments_at_recursion_depth(max_recursion_depth, bonds_to_break_this_molecule);

  //#define DEBUG_ADJUST_MAX_RECURSION_DEPTH
#ifdef DEBUG_ADJUST_MAX_RECURSION_DEPTH
  cerr << bonds_to_break_this_molecule << " bonds, max_recursion_depth " << max_recursion_depth << " yields " << fragments_that_will_be_produced << " fragments\n";
#endif

  // Commit out @ 20120713.
  // We should consist to the input restriction of recursion_depth.
  // if (fragments_that_will_be_produced <= max_fragments_per_molecule)
  // {
  //      _max_recursion_depth_this_molecule = bonds_to_break_this_molecule;
  //      return fragments_that_will_be_produced;
  // }

  if (_max_recursion_depth_this_molecule > bonds_to_break_this_molecule)
    _max_recursion_depth_this_molecule = bonds_to_break_this_molecule;

  while (fragments_that_will_be_produced > max_fragments_per_molecule && _max_recursion_depth_this_molecule > 1)
  {
    _max_recursion_depth_this_molecule--;
    fragments_that_will_be_produced = fragments_that_will_be_produced / _max_recursion_depth_this_molecule;

#ifdef DEBUG_ADJUST_MAX_RECURSION_DEPTH
    cerr << "at " << _max_recursion_depth_this_molecule << ", " << bonds_to_break_this_molecule << " will produce " << fragments_that_will_be_produced << '\n';
#endif
  }

#ifdef DEBUG_ADJUST_MAX_RECURSION_DEPTH
  cerr << "At " << _max_recursion_depth_this_molecule << ", " << bonds_to_break_this_molecule << " will produce " << fragments_that_will_be_produced << endl;
#endif

  return fragments_that_will_be_produced;
}

void
Dicer_Arguments::_add_new_smiles_to_global_hashes(const IWString & smiles)
{
  if (!add_new_strings_to_hash)
    return;

  unsigned int s = smiles_to_id.size();

  assert (s == static_cast<unsigned int>(id_to_smiles.number_elements()));

  smiles_to_id[smiles] = s;

  IWString * tmp = new IWString(smiles);
  id_to_smiles.add(tmp);

  _fragments_found_this_molecule[s] = 1;

  molecules_containing[smiles]++;

  return;
}

static void
add_name_of_parent_to_starting_parent_array(const IWString & smiles,
                                            const IWString & name_of_parent,
                                            IW_STL_Hash_Map_String & starting_parent)
{
  if (! name_of_parent.contains(' '))
  {
    starting_parent[smiles] = name_of_parent;
    return;
  }

  IWString tmp(name_of_parent);
  tmp.truncate_at_first(' ');

  starting_parent[smiles] = tmp;

  return;
}

static int
ok_heteromatom_count (const Molecule & m,
                      const int minh)
{
  const int matoms = m.natoms();

  int heteroatoms = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (6 == m.atomic_number(i))
      continue;

    heteroatoms++;

    if (heteroatoms >= minh)
      return 1;
  }

  return 0;
}

static int
contains_user_specified_queries (Molecule & m)
{
  int nq = user_specified_queries_that_fragments_must_contain.number_elements();

  if (0 == nq)
    return 1;

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++)
  {
    if (user_specified_queries_that_fragments_must_contain[i]->substructure_search(target))
      return 1;
  }

  fragments_discarded_for_not_containing_user_specified_queries++;

  return 0;
}

static int
do_check_for_lost_chirality (Molecule & m)
{
  int rc = 0;

  const int nc = m.chiral_centres();
  if (nc == 0) {
    return 0;
  }

  for (int i = nc - 1; i >= 0; --i)
  {
    const Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    if (! is_actually_chiral(m, c->a()))
    {
      m.remove_chiral_centre_at_atom(c->a());
      rc++;
    }
  }

  return rc;
}

/*
   We have encountered a new smiles. Add it to all the counters
 */

int
Dicer_Arguments::contains_fragment(Molecule & m)
{
  //cerr << "Found fragment '" << smiles << "'\n";

  if (! contains_user_specified_queries(m)) {
    return 0;
  }

  if (! OkAtomCount(m)) {
    return 0;
  }

  if (fragments_must_contain_heteroatom && ! ok_heteromatom_count(m, fragments_must_contain_heteroatom)) {
    return 0;
  }

  if (check_for_lost_chirality) {
    do_check_for_lost_chirality(m);
  }

  if (fragments_indexed_by_unique_smiles) {
    return MaybeAddToHash(m.unique_smiles());
  } else {
    return MaybeAddToHash(m.smiles());
  }
}

int
Dicer_Arguments::MaybeAddToHash(const IWString& smiles) {

  // std::cerr << smiles << " ####" << _recursion_depth << '\n';

  IW_STL_Hash_Map_uint::iterator f1 = smiles_to_id.find(smiles);

  if (f1 == smiles_to_id.end())         // fragment never been seen before
  {
    if (! add_new_strings_to_hash)
      return 0;

    _add_new_smiles_to_global_hashes(smiles);

    //  cerr << "New fragment generated from '" << _current_molecule->name() << "', acc " << accumulate_starting_parent_information <<endl;
    if (accumulate_starting_parent_information)
      add_name_of_parent_to_starting_parent_array(smiles, _current_molecule->name(), starting_parent);

    return 1;                // new smiles encountered this molecule
  }

  // Fragment had previously been encountered - either this molecule or elsewhere

  unsigned int s = (*f1).second;          // get the number for this smiles string

  if (0 == _fragments_found_this_molecule[s])
    molecules_containing[smiles]++;

  _fragments_found_this_molecule[s]++;

  return 0;          // smiles had been encountered before
}

//#define DEBUG_BIT_COMPARISONS

int
Dicer_Arguments::is_unique(Molecule & m)
{
#ifdef DEBUG_BIT_COMPARISONS
  cerr << "Examining '" << m.smiles() << "'\n";
#endif

  IW_Bits_Base * b = new IW_Bits_Base(_nbits);

  for (const Atom* a: m) {
    std::optional<atom_number_t> x = InitialAtomNumber(a);
    if (! x) {
      continue;
    }

    b->set(*x);
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    std::optional<atom_number_t> x = InitialAtomNumber(m.atomi(i));
    if (! x) {
      continue;
    }

    b->set(*x);
  }

//assert (b->nset() == matoms);    if we have added environment atoms, this will not be true

#ifdef DEBUG_BIT_COMPARISONS
  cerr << "Formed bit   ";
  b->printon(cerr);
  cerr << endl;
#endif

  int n = _stored_bit_vectors.number_elements();

  for (int i = 0; i < n; i++) {
#ifdef DEBUG_BIT_COMPARISONS
    cerr << "Compare with ";
    _stored_bit_vectors[i]->printon(cerr);
    cerr << " ? " << ((*b) == (*(_stored_bit_vectors[i]))) << endl;
#endif

    if ((*b) == *(_stored_bit_vectors[i])) {
      if (record_presence_and_absence_only) {                     // sept 2011, I don't think is useful, make both branches do the same thing
        delete b;
        return 1;
      }
      else
      {
        //break;
        delete b;
        return 1;
      }
    }
  }

#ifdef DEBUG_BIT_COMPARISONS
  cerr << "New fragment accepted\n";
#endif

  _stored_bit_vectors.add(b);

  // The following code was changed by madlee@20120801.
  // we set the -B atype=<atomtype> new meaning now.
  // Maybe we will create a new opiton for this function in future.
  contains_fragment(m);
  return 1;
  // end of change.

}

int
Dicer_Arguments::write_fragments_found_this_molecule(Molecule & m,
                                          int breakable_bonds,
                                          DicerFragmentOutput & output) const
{
  if (! _write_diced_molecules) {
    return 1;
  }

  if (! output.output_is_proto()) {
  } else if (write_smiles_and_complementary_smiles) {
    return 1;
  } else {
    return WriteAsProto(m, breakable_bonds, output);
  }

  int smiles_written = write_parent_smiles;

  for (ff_map::const_iterator f = _fragments_found_this_molecule.begin(); f != _fragments_found_this_molecule.end(); ++f) {
    if (vf_structures_per_page) {
      if (smiles_written > 0 && 0 == smiles_written % vf_structures_per_page)
        output << _current_molecule->smiles() << ' ' << _current_molecule->name() << " PARENT\n";
      smiles_written++;
    }

    int fragment_number = (*f).first;
    int count = (*f).second;

    const IWString * s = id_to_smiles[fragment_number];

    output << (*s) << ' ' << m.name() << " FRAGID=" << fragment_number << " N= " << count;
    if (append_atom_count) {
      const int n = count_atoms_in_smiles(*s);
      output << " NATOMS " << n;
      output << " FF " << (static_cast<float>(n) / static_cast<float>(m.natoms()));
    }
    output << '\n';
  }

  if (vf_structures_per_page && 0 != smiles_written % vf_structures_per_page)
  {
    do
    {
      output << "C\n";
      smiles_written++;
    }
    while (0 != smiles_written % vf_structures_per_page);
  }

  return 1;
}

// Load all the information about the breakages for this molecule into `proto`.
int
Dicer_Arguments::ToProto(Molecule& m,
                         int breakable_bonds,
                         dicer_data::DicedMolecule& proto) const {
  proto.set_name(m.name().data(), m.name().length());
  if (fragments_indexed_by_unique_smiles) {
    const IWString& usmi = m.unique_smiles();
    proto.set_smiles(usmi.data(), usmi.length());
  } else {
    const IWString& smi = m.smiles();
    proto.set_smiles(smi.data(), smi.length());
  }
  proto.set_natoms(m.natoms());
  proto.set_nbonds(breakable_bonds);

  for (const auto& [fragment_number, count] : _fragments_found_this_molecule) {
    const IWString * s = id_to_smiles[fragment_number];

    dicer_data::DicerFragment* frag = proto.add_fragment();
    frag->set_smi(s->data(), s->length());
    frag->set_nat(count_atoms_in_smiles(*s));
    frag->set_n(count);
    frag->set_id(fragment_number);
  }

  return 1;
}

int
Dicer_Arguments::WriteAsProto(Molecule& m,
                              int breakable_bonds,
                              DicerFragmentOutput& output) const {
  dicer_data::DicedMolecule proto;
  ToProto(m, breakable_bonds, proto);
  return output.Write(proto);
}

int
Dicer_Arguments::produce_fingerprint(Molecule & m,
                                     IWString_and_File_Descriptor & output) const
{
  Sparse_Fingerprint_Creator sfp;

  for (ff_map::const_iterator f = _fragments_found_this_molecule.begin(); f != _fragments_found_this_molecule.end(); ++f)
  {
    int fragment_number = (*f).first;
    int count = (*f).second;

    sfp.hit_bit(fragment_number, count);
  }

  if (! function_as_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  IWString tmp;
  sfp.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);
  output << tmp << '\n';

  //sfp.append_daylight_ascii_form_with_counts_encoded(fingerprint_tag, buffer);

  if (! function_as_filter)
    output << "|\n";

  return 1;
}

int
Dicer_Arguments::store_atom_types(Molecule & m,
                                  Atom_Typing_Specification & ats)
{
  if (nullptr == _atom_type)
  {
    const int matoms = m.natoms();

    _atom_type = new uint32_t[matoms+matoms];
    _isosave = _atom_type + matoms;
  }

  return ats.assign_atom_types(m, _atom_type);
}

int
Breakages::_destroys_too_much_aromaticity (Molecule & m,
                                          const int aromatic_atom_count,
                                          const Fused_Ring_Breakage & f) const
{
  if (_destroys_too_much_aromaticity(m, aromatic_atom_count, f.a11(), f.a1(), f.a2(), f.a21()))
    return 1;

  if (_destroys_too_much_aromaticity(m, aromatic_atom_count, f.a12(), f.a1(), f.a2(), f.a22()))
    return 1;

  return 0;
}

int
Breakages::_destroys_too_much_aromaticity (Molecule & m,
                                           const int aromatic_atom_count,
                                           atom_number_t a1,
                                           atom_number_t a2,
                                           atom_number_t a3,
                                           atom_number_t a4) const
{
  const Bond * b = m.bond_between_atoms(a1, a2);
  const bond_type_t b12 = b->btype();
  b = m.bond_between_atoms(a3, a4);
  const bond_type_t b34 = b->btype();

  assert (m.are_bonded(a1, a2));
  assert (m.are_bonded(a3, a4));
  m.remove_bond_between_atoms(a1, a2);
  m.remove_bond_between_atoms(a3, a4);

  //cerr << "After removing bonds " << m.smiles() << "\n";

  int a = m.aromatic_atom_count();         // after removing the bonds

  int rc = 0;

  if (0 == a)          // aromaticity destroyed, very bad
    rc = 1;
  else if (aromatic_atom_count - a > 6)         // more than 6 aromatic atoms lost, bad
    rc = 1;

  m.add_bond(a1, a2, b12);
  m.add_bond(a3, a4, b34);

  return rc;
}

/*
   We are trying to find two atoms adjacent to zatom, where the bonds are
   in just one ring.
 */

static int
identify_attached (const Atom & a,
                   const atom_number_t zatom,
                   const atom_number_t avoid,
                   atom_number_t & a1,
                   atom_number_t & a2)
{
  assert (3 == a.ncon());

  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = a[i];

    if (1 != b->nrings())
      continue;

    const atom_number_t j = b->other(zatom);

    if (j == avoid)               // not really necessary, this bond should be in > 1 ring
      continue;

    if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else
    {
      a2 = j;
      return 1;
    }
  }

  return INVALID_ATOM_NUMBER != a2;
}

static bool
two_ring_of_spiro (Molecule& m, const Atom & center, int icenter,
                   atom_number_t atoms[4])
{
  if (4 == center.ncon())
  {
    for (int i = 0; i < 4; ++i)
    {
      if (center[i]->nrings() != 1)
      {
        return 0;
      }
      if (center[i]->btype() != SINGLE_BOND)
      {
        return 0;
      }
    }
  }
  else
  {
    return 0;
  }

  atoms[0] = center[0]->other(icenter);
  atoms[1] = center[1]->other(icenter);
  atoms[2] = center[2]->other(icenter);
  atoms[3] = center[3]->other(icenter);

  if (! m.in_same_ring(atoms[0], atoms[1]))
  {
    if (m.in_same_ring(atoms[0], atoms[2]))
      std::swap(atoms[1], atoms[2]);
    else
      std::swap(atoms[1], atoms[3]);
  }

  return m.in_same_ring(atoms[0], atoms[1]) && m.in_same_ring(atoms[2], atoms[3]);
}


static int
contains_heteroatoms(const Molecule & m,
                     const Ring & r)
{
  const int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    if (6 != m.atomic_number(j))
      return 1;
  }

  return 0;
}

int
Breakages::identify_aromatic_rings_to_de_aromatise (Molecule & m)
{
  const int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);

    if (!r->is_aromatic())
      continue;

    if (!contains_heteroatoms(m, *r))
      continue;

    Dearomatise * d = new Dearomatise(*r);

    _transformations.add(d);
  }

  return rc;
}

//#define DEBUG_IDENTIFY_CROSS_RING_BONDS_TO_BE_BROKEN

static Fused_Ring_Breakage *
create_fused_ring_breakage(Molecule& m,
                           const atom_number_t a11, const atom_number_t a1, const atom_number_t a12,
                           const atom_number_t a21, const atom_number_t a2, const atom_number_t a22)
{
  const Bond* bcenter = m.bond_between_atoms(a1, a2);
  const Bond* b11 = m.bond_between_atoms(a1, a11);
  const Bond* b21 = m.bond_between_atoms(a2, a21);

  bool can_creat = true;

  if (b11->btype() != b21->btype())
  {
    // Only cut if the two bonds are in same bond order;
    can_creat = false;
  }
  else if (b11->btype() != SINGLE_BOND || b21->btype() != SINGLE_BOND)
  {
    can_creat = bcenter->is_aromatic() && b11->is_aromatic() && b21->is_aromatic();
  }

  if (can_creat)
  {
    Fused_Ring_Breakage * f1 = new Fused_Ring_Breakage(a11, a1, a12, a21, a2, a22);
    const Bond* b12 = m.bond_between_atoms(a1, a12);
    const Bond* b22 = m.bond_between_atoms(a2, a22);
    if (bcenter->is_aromatic() && b12->is_aromatic() && b22->is_aromatic())
    {
      f1->set_aromatic_in_parent(1);
    }

    return f1;
  }
  else
  {
    return nullptr;
  }


}
int
Breakages::identify_cross_ring_bonds_to_be_broken (Molecule & m)
{
  const int nr = m.nrings();

  if (nr < 2)
    return 0;

  m.compute_aromaticity_if_needed();

  const int nb = m.nedges();

  int result = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (2 != b->nrings())
      continue;

    atom_number_t a1 = b->a1();

    const Atom * aa1 = m.atomi(a1);

    if (3 != aa1->ncon())
      continue;

    atom_number_t a2 = b->a2();

    const Atom * aa2 = m.atomi(a2);

    if (3 != aa2->ncon())
      continue;

    atom_number_t a11, a12;
    if (!identify_attached(*aa1, a1, a2, a11, a12))
      continue;

    atom_number_t a21, a22;
    if (!identify_attached(*aa2, a2, a1, a21, a22))
      continue;

    if (a11 == a22 || m.in_same_ring(a11, a22))
      iwswap(a11, a12);

    Fused_Ring_Breakage * f1 = create_fused_ring_breakage(m, a11, a1, a12, a21, a2, a22);
    Fused_Ring_Breakage * f2 = create_fused_ring_breakage(m, a12, a1, a11, a22, a2, a21);

    if (f1 || f2)
    {
      _transformations.add(new PairTransformation(f1, f2));
      ++result;
    }
  }

  if (verbose > 1)
    cerr << m.name() << " identified " << result << " fused ring bonds to break\n";

  return result;
}


///////////////////////////////////////////////////////////////////////
// int Breakages::identify_spiro_ring_bonds_to_be_broken (Molecule & m)
// by madlee

// #define DEBUG_IDENTIFY_SPIRO_RING_TO_BE_BROKEN
int
Breakages::identify_spiro_ring_bonds_to_be_broken (Molecule & m)
{
  const int nr = m.nrings();

  if (nr < 2)
    return 0;

  const int na = m.natoms();

  int result = 0;

  for (int i = 0; i < na; i++)
  {
    const Atom * a = m.atomi(i);

    if (a->ncon() != 4)
      continue;

    atom_number_t neighbors[4];
    if (!two_ring_of_spiro(m, *a, i, neighbors))
    {
      continue;
    }

    Spiro_Ring_Breakage* s1 = new Spiro_Ring_Breakage(i, neighbors[0], neighbors[1], neighbors[2], neighbors[3]);
    Spiro_Ring_Breakage* s2 = new Spiro_Ring_Breakage(i, neighbors[2], neighbors[3], neighbors[0], neighbors[1]);
    _transformations.add(new PairTransformation(s1, s2));
    ++result;
  }


  if (verbose > 1)
    cerr << m.name() << " identified " << result << " spiro ring bonds to break\n";

  return result;
}
// int Breakages::identify_spiro_ring_bonds_to_be_broken (Molecule & m)
///////////////////////////////////////////////////////////////////////

/*
   We want to identify those parts of the ring that are the edge.
   First atom in the result must be one of the fusion points, then
   the atoms around the edge, and finally the atom at the other
   side of the fusion.
 */
static int
identify_ring_edge (Molecule & m,
                    const Ring & r,
                    Set_of_Atoms & ring_edge)
{
  int ring_size = r.number_elements();

  int first_nrings2 = -1;
  int second_nrings2 = -1;

  for (int i = 0; i < ring_size; i++)
  {
    if (2 != m.nrings(r[i]))
      continue;

    if (first_nrings2 < 0)
      first_nrings2 = i;
    else if (second_nrings2 < 0)
      second_nrings2 = i;
    else
      return 0;                     // gack, more than two fusions
  }

  int istart;
  int direction;

  if (first_nrings2 + 1 == second_nrings2)
  {
    istart = first_nrings2;
    direction = 1;
  }
  else if (0 == first_nrings2 && (ring_size - 1) == second_nrings2)
  {
    istart = ring_size - 1;
    direction = -1;
  }
  else         // spiro fusion or something? Cannot imagine
  {
    return 0;
  }

  ring_edge.add(r[istart]);
  for (int i = 0; i < (ring_size - 2); i++)
  {
    ring_edge.add(r.next_after_wrap(istart, direction));
  }
  ring_edge.add(r.next_after_wrap(istart, direction));

  return 1;
}


/*
   We have identified a ring to break apart. but we do not want to totally
   destroy a ring system, so we identify an atom in an adjacent ring. Only
   if that atom is still a ring atom will we destroy our ring
 */


atom_number_t
identify_atom_in_adjacent_ring (Molecule & m,
                                const Set_of_Atoms & s)
{
  atom_number_t zatom = s[0];

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  assert (acon > 2);          // first atom in the set is supposed to be fused

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (s.contains(j))
      continue;

    if (m.is_non_ring_atom(j))
    {
      continue;
    }

    return j;
  }

  return INVALID_ATOM_NUMBER;
}

#ifdef FROMCEX

static Ring_Bond_Breakage*
create_ring_bond_breakage (Molecule & m, const Ring* r)
{
  if (r->strongly_fused_ring_neighbours())         // too complex
    return nullptr;

  if (1 != r->fused_ring_neighbours())         // not a terminal fing
    return nullptr;

  if (r->fused_ring_neighbours() > 1)          // too complex, but, potentially limiting
    return nullptr;

  if (r->is_aromatic() && !r->is_fused())          // don't handle these
    return nullptr;

  //  At this stage, we have a "terminal" fused ring. Identify the atoms that are NOT
  //  part or adjacent to the ring

  Set_of_Atoms not_part_of_adjacent_ring;
  if (!identify_ring_edge (m, *r, not_part_of_adjacent_ring))
    return nullptr;

  not_part_of_adjacent_ring.add(not_part_of_adjacent_ring[0]);        // close the ring, we will be processing bonds

  Ring_Bond_Breakage * tmp = new Ring_Bond_Breakage(not_part_of_adjacent_ring, r->is_aromatic());

  atom_number_t atom_in_adjacent_ring = identify_atom_in_adjacent_ring(m, not_part_of_adjacent_ring);

  assert (INVALID_ATOM_NUMBER != atom_in_adjacent_ring);

  tmp->set_atom_that_must_be_a_ring_atom(atom_in_adjacent_ring);
  return tmp;
}

static int
get_two_ring_on_bond(Molecule & m, const Bond* b, const Ring* result[])
{
  assert (b->nrings() == 2);
  const int nr = m.nrings();
  int n = 0;
  for (int i = 0; i < nr; ++i)
  {
    const Ring* r = m.ringi(i);
    if (r->contains_bond(b->a1(), b->a2()))
    {
      result[n++] = r;
    }
  }
  return n;
}
#endif

int
Breakages::identify_ring_bonds_to_be_broken (Molecule & m)
{
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = m.ringi(i);

    if (r->strongly_fused_ring_neighbours())               // too complex
      continue;

    if (1 != r->fused_ring_neighbours())               // not a terminal ring
      continue;

    if (r->fused_ring_neighbours() > 1)                // too complex, but, potentially limiting
      continue;

    if (r->is_aromatic() && !r->is_fused())                // don't handle these
      continue;

    //  At this stage, we have a "terminal" fused ring. Identify the atoms that are NOT
    //  part or adjacent to the ring

    Set_of_Atoms not_part_of_adjacent_ring;
    if (!identify_ring_edge (m, *r, not_part_of_adjacent_ring))
      continue;

    not_part_of_adjacent_ring.add(not_part_of_adjacent_ring[0]);              // close the ring, we will be processing bonds

    int n = not_part_of_adjacent_ring.number_elements();

    Ring_Bond_Breakage * tmp = new Ring_Bond_Breakage(not_part_of_adjacent_ring, r->is_aromatic());

    atom_number_t atom_in_adjacent_ring = identify_atom_in_adjacent_ring(m, not_part_of_adjacent_ring);

    assert (INVALID_ATOM_NUMBER != atom_in_adjacent_ring);

    tmp->set_atom_that_must_be_a_ring_atom(atom_in_adjacent_ring);

    _transformations.add(tmp);

    //  cerr << "Identified ring bond breakages " << not_part_of_adjacent_ring << ' ' << m.isotopically_labelled_smiles() << " ex ring " << atom_in_adjacent_ring << endl;

    //  We have identified something like 0 1 2 3 4 5 0 (break bond 0-5),
    //  we now also need to get 1 0 5 4 3 2 1 (break bond 2 1)

    Set_of_Atoms s;
    s.resize(n);

    s.add(not_part_of_adjacent_ring[1]);
    s.add(not_part_of_adjacent_ring[0]);
    for (int j = n-2; j > 1; j--)
    {
      s.add(not_part_of_adjacent_ring[j]);
    }
    s.add(not_part_of_adjacent_ring[1]);
    //  cerr << "From " << not_part_of_adjacent_ring << " generated " << s << endl;

    tmp = new Ring_Bond_Breakage(s, r->is_aromatic());
    tmp->set_atom_that_must_be_a_ring_atom(atom_in_adjacent_ring);
    _transformations.add(tmp);
  }

  return rc;
}

#ifdef FROMCEX

int
Breakages::identify_ring_bonds_to_be_broken(Molecule & m)
{
  const int nr = m.nrings();

  if (nr < 2)
    return 0;

  m.compute_aromaticity_if_needed();

  const int nb = m.nedges();

  int result = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (2 != b->nrings())
      continue;

    const Ring* ring[2];

    if (2 != get_two_ring_on_bond(m, b, ring))
      continue;

    Ring_Bond_Breakage* f1 = create_ring_bond_breakage(m, ring[0]);
    Ring_Bond_Breakage* f2 = create_ring_bond_breakage(m, ring[1]);

    if (f1 || f2)
    {
      _transformations.add(new PairTransformation(f1, f2));
      ++result;
    }
  }

  return result;
}
#endif

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "  -s <smarts>   smarts for cut points - breaks bond between 1st and 2nd matched atom\n";
  cerr << "  -q <...>      query specification (alternative to -s option)\n";
  cerr << "  -n <smarts>   smarts for bonds to NOT break\n";
  cerr << "  -N <...>      query specification(s) for bonds to NOT break\n";
  cerr << "  -k <number>   maximum number of bonds to simultaneously break (def 10)\n";
  cerr << "  -X <number>   maximum number of fragments per molecule to produce\n";
  cerr << "  -z i          ignore molecules not matching any of the queries\n";
  cerr << "  -I <iso>      specify isotope for join points, enter '-I help' for more info\n";
  cerr << "  -J <tag>      produce fingerprint with tag <tag>\n";
  cerr << "  -K <fname>    READ=<fname>, WRITE=<fname> manage cross reference files bit#->smiles\n";
  cerr << "  -m <natoms>   discard fragments with fewer than <natoms> atoms\n";
  cerr << "  -M <natoms>   discard fragments with more than <natoms> atoms\n";
  cerr << "  -M maxnr=<nat> discard fragments with more than <nat> non ring atoms\n";
  cerr << "  -c            discard chirality from input molecules - also discards cis-trans bonds\n";
  cerr << "  -C def        write smiles and complementary smiles\n";
  cerr << "  -C auto       write smiles and complementary smiles with auto label on break points.\n";
  cerr << "  -C iso=<n>    write smiles and complementary smiles with isotope <n> on break points.\n";
  cerr << "  -C atype      write smiles and complementary smiles with atom type labels on break points\n";
  cerr << "  -a            write only smallest fragments to -C ... output\n";
  cerr << "  -T <e1=e2>    one or more element transformations\n";
  cerr << "  -P ...        atom typing specification (used for isotopes and/or complementary smiles)\n";
  cerr << "  -B ...        various other options, enter '-B help' for info\n";
  cerr << "  -Z <fname>    write fully broken parent molecules to <fname>\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -G ...        standard smiles options\n";
  cerr << "  -S <fname>    output file (default to stdout)\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

/*
   We need to be careful about how we do this. My first implementation just
   created a molecule for each, a Molecule_to_Match, and a query. But this
   quickly becomes HUGE.

   So, scan through and process things by groups with common counts and
   being very careful about which molecules can actually be compared
 */

static int * atomic_number_to_array_position = nullptr;

static int
initialise_atomic_number_to_array_position()
{
  atomic_number_to_array_position = new_int(HIGHEST_ATOMIC_NUMBER + 1);

  atomic_number_to_array_position[1] = 1;
  atomic_number_to_array_position[6] = 2;
  atomic_number_to_array_position[7] = 3;
  atomic_number_to_array_position[8] = 4;
  atomic_number_to_array_position[9] = 5;
  atomic_number_to_array_position[15] = 6;
  atomic_number_to_array_position[16] = 7;
  atomic_number_to_array_position[17] = 8;
  atomic_number_to_array_position[35] = 9;
  atomic_number_to_array_position[53] = 10;

  return 1;
}

class Fragment_Info
{
  private:
    extending_resizable_array<int> _mf;

    int _natoms;            // atoms in fragment

    int _ndx;

    int _times_hit;

    int _is_subset_of_another;

  public:
    Fragment_Info();

    void set_ndx(int s)
    {
      _ndx = s;
    }
    int ndx() const
    {
      return _ndx;
    }

    void set_times_hit(int s)
    {
      _times_hit = s;
    }
    int times_hit() const
    {
      return _times_hit;
    }

    int natoms() const
    {
      return _natoms;
    }
    void set_natoms(int s)
    {
      _natoms = s;
    }

    int initialise_mf_array(Molecule & m);

    int is_subset_of(const Fragment_Info &) const;

    int is_subset_of_another() const
    {
      return _is_subset_of_another;
    }
    void set_is_subset_of_another(int s)
    {
      _is_subset_of_another = s;
    }
};

Fragment_Info::Fragment_Info()
{
  _ndx = -1;

  _times_hit = 0;

  _natoms = 0;

  _is_subset_of_another = 0;

  return;
}

int
Fragment_Info::initialise_mf_array(Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    _mf[atomic_number_to_array_position[z]]++;

    int h = m.hcount(i);

    _mf[atomic_number_to_array_position[1]] += h;
  }

  return matoms;
}

int
Fragment_Info::is_subset_of (const Fragment_Info & rhs) const
{
  if (_times_hit != rhs._times_hit)          // must be present in exactly the same number of molecules
    return 0;

  // These two fragments are present in the same number of molecule. Let's check formulae

  int n = _mf.number_elements();

  if (n > rhs._mf.number_elements())         // we have more elements than RHS
    return 0;

  for (int i = 0; i < n; i++)         // same number of elements, check counts
  {
    if (_mf[i] > rhs._mf[i])
      return 0;
  }

  // Our formula is a subset of RHS

  return 1;
}

static void
mark_duplicates (Fragment_Info * frag_info,
                 const resizable_array<int> & s)
{
  int n = s.number_elements();

#ifdef DEBUG_MARK_DUPLICATES
  cerr << "Looking for subset relationships in " << n << " fragments\n";
#endif

  // First initialise the arrays of molecules and targets

  resizable_array_p<Molecule> m;
  resizable_array_p<Molecule_to_Match> target;

  m.resize(n);
  target.resize(n);

  for (int i = 0; i < n; i++)
  {
    int j = s[i];

    const Fragment_Info & fj = frag_info[j];

    int ndx = fj.ndx();

    const IWString & smiles = *(id_to_smiles[ndx]);
    //  cerr <<  smiles << " ndx " << ndx << endl;

    Molecule * tm = new Molecule;
    tm->build_from_smiles(smiles);

    //  Get rid of isotopes if we have applied any.

    if (isotope_for_join_points || apply_environment_isotopic_labels || isotopic_label_is_recursion_depth)
      tm->transform_to_non_isotopic_form();

    m.add(tm);

    Molecule_to_Match * t = new Molecule_to_Match(tm);
    target.add(t);
  }

  Molecule_to_Query_Specifications mqs;

  for (int i = 0; i < n; i++)
  {
    Substructure_Query q;
    q.create_from_molecule(*(m[i]), mqs);
    q.set_max_matches_to_find(1);

//  cerr << "Begin molecule " << i << " contains " << m[i]->natoms() << " atoms\n";
    for (int j = i + 1; j < n; j++)
    {
      if (0 == q.substructure_search(*(target[j])))
        continue;

      frag_info[s[i]].set_is_subset_of_another(1);
      break;
    }
  }
}

class Frag_Info_Comparator
{
public:
int operator() (const Fragment_Info &, const Fragment_Info &) const;
};

int
Frag_Info_Comparator::operator() (const Fragment_Info & f1,
  const Fragment_Info & f2) const
{
  if (f1.times_hit() < f2.times_hit())
    return -1;

  if (f1.times_hit() > f2.times_hit())
    return 1;

  if (f1.natoms() < f2.natoms())
    return -1;

  if (f1.natoms() > f2.natoms())
    return 1;

  if (f1.is_subset_of(f2))
    return -1;

  if (f2.is_subset_of(f1))
    return 1;

  return 0;              // not resolved
}

static int
do_eliminate_fragment_subset_relations()
{
  int n = smiles_to_id.size();

  Fragment_Info * frag_info = new Fragment_Info[n];

  if (nullptr == frag_info)
  {
    cerr << "Cannot allocate " << n << " fragment info objects\n";
    return 0;
  }

  std::unique_ptr<Fragment_Info[]> free_frag_info(frag_info);

  set_input_aromatic_structures(1);

  int ndx = 0;

  for (int i = 0; i < n; i++)
  {
    const IWString & smiles = *(id_to_smiles[i]);

    Molecule tmp;
    if (!tmp.build_from_smiles(smiles))
    {
      cerr << "Cannot interpret smiles '" << smiles << "'\n";
      return 0;
    }

    frag_info[ndx].initialise_mf_array(tmp);
    frag_info[ndx].set_ndx(i);
    frag_info[ndx].set_natoms(tmp.natoms());
    ndx++;
  }

  if (verbose)
    cerr << "Will filter " << ndx << " fragments for possible subset relationships\n";

  // We make this easier if we sort things

  Frag_Info_Comparator frag_info_comparator;
  iwqsort(frag_info, ndx, frag_info_comparator);

  //#define ECHO_SORTED_LIST
#ifdef ECHO_SORTED_LIST
  for (int i = 0; i < ndx; i++)
  {
    const Fragment_Info & fi = frag_info[i];
    int j = fi.ndx();

    const IWString & smiles = *(id_to_smiles[j]);

    cerr << smiles << " ndx = " << fi.ndx() << " times_hit " << fi.times_hit() << " natoms " << fi.natoms() << endl;
  }
#endif

  // Now that things are sorted, we can do duplicate elimination within groups

  int rc = 0;          // the number of duplicates we find

  for (int i = 0; i < ndx; i++)
  {
    const Fragment_Info & fi = frag_info[i];

    if (fi.is_subset_of_another())
      continue;

    resizable_array<int> do_comparisons;
    do_comparisons.add(i);

    for (int j = i + 1; j < ndx; j++)
    {
      const Fragment_Info & fj = frag_info[j];

      if (fj.is_subset_of_another())
        continue;

      if (fi.is_subset_of(fj))
        do_comparisons.add(j);
    }

    if (do_comparisons.number_elements() > 1)
      mark_duplicates(frag_info, do_comparisons);
  }

  int duplicates_found = 0;
  for (int i = 0; i < ndx; i++)
  {
    const Fragment_Info & fi = frag_info[i];

    if (fi.is_subset_of_another())
    {
      //    int j = fi.ndx();
      //    cerr << " j = " << j << " is a subset\n";
      duplicates_found++;
    }
  }

  if (verbose)
    cerr << duplicates_found << " of the " << ndx << " fragments are full subsets\n";

  return rc;
}


enum class Env {
  kUndefined = 0,
  kTerminal = 1,
  kAromatic = 2,
  kSaturatedCarbon = 3,
  kSaturatedNitrogen = 4,
  kSaturatedOxygen = 5,
  kUnsaturatedCarbon = 6,
  kUnsaturatedNitrogen = 7,
  // This will likely never be used, since we usually
  // do not break =O bonds.
  kUnsaturatedOxygen = 8,
  kAmideCarbon = 9,
  kAmideNitrogen = 10,
  kHalogen = 11,
  kNitro = 12,
  kAliphaticRingCarbon = 13,
  kAliphaticRingNitrogen = 14,
  // Used to size the elemements array.
  kSize = 15,
};

// An array of the environment elements, index by the Env enum.
// set up in InitialiseEnvironmentElements
static const Element* env_element[static_cast<int>(Env::kSize)];

static int
do_apply_atomic_number_isotopic_labels(Molecule & m,
                                       const atom_number_t j1,
                                       const atom_number_t j2)
{
  const atomic_number_t z2 = m.atomic_number(j2);

  return m.set_isotope(j1, z2);
}

static void
do_apply_isotopic_labels (Molecule & m,
                          const atom_number_t a1,
                          const atom_number_t a2)
{
  if (increment_isotope_for_join_points)
  {
    m.set_isotope(a1, m.isotope(a1) + increment_isotope_for_join_points + isotope_for_join_points);
    m.set_isotope(a2, m.isotope(a2) + increment_isotope_for_join_points + isotope_for_join_points);
  }
  else
  {
    m.set_isotope(a1, isotope_for_join_points);
    m.set_isotope(a2, isotope_for_join_points);
  }

  return;
}

// Return a result if `zatom` is singly bonded to an
// atom with atomic number `target`.
std::optional<atom_number_t>
SinglyBondedTo(const Molecule& m,
               atom_number_t zatom,
               atomic_number_t target) {
  for (const Bond* b : m.atom(zatom)) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (m.atomic_number(j) == target) {
      return j;
    }
  }

  return std::nullopt;
}

// Return a result if `zatom` is doubly bonded to an
// atom with atomic number `target`.
std::optional<atom_number_t>
DoublyBondedTo(const Molecule& m,
               atom_number_t zatom,
               atomic_number_t target) {
  for (const Bond* b : m.atom(zatom)) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (m.atomic_number(j) == target) {
      return j;
    }
  }

  return std::nullopt;
}

int
IsAmideCarbon(const Molecule& m,
              atom_number_t zatom) {
  const Atom* a = m.atomi(zatom);
  if (a->ncon() != 3) {
    return 0;
  }

  atom_number_t doubly_bonded_oxygen = INVALID_ATOM_NUMBER;
  atom_number_t nitrogen = INVALID_ATOM_NUMBER;
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (b->is_double_bond() && m.atomic_number(j) == 8) {
      doubly_bonded_oxygen = j;
    } else if (b->is_single_bond() && m.atomic_number(j) == 7) {
      nitrogen = j;
    }
  }

  if (INVALID_ATOM_NUMBER == doubly_bonded_oxygen) {
    return 0;
  }

  if (INVALID_ATOM_NUMBER == nitrogen) {
    return 0;
  }

  return 1;
}

int
IsAmideNitrogen(const Molecule& m,
                atom_number_t zatom) {
  std::optional<atom_number_t> carbon = SinglyBondedTo(m, zatom, 6);
  if (! carbon) {
    return 0;
  }

  std::optional<atom_number_t> oxygen = DoublyBondedTo(m, *carbon, 8);
  if (! oxygen) {
    return 0;
  }

  return 1;
}

int
IsAmideOxygen(const Molecule& m,
              atom_number_t zatom) {
  std::optional<atom_number_t>carbon = DoublyBondedTo(m, zatom, 6);
  if (! carbon) {
    return 0;
  }

  std::optional<atom_number_t> nitrogen = SinglyBondedTo(m, *carbon, 7);
  if (! nitrogen) {
    return 0;
  }

  return 1;
}

int
IsAmide(const Molecule& m,
        atom_number_t zatom) {
  const Atom* a = m.atomi(zatom);
  if (a->atomic_number() == 6) {
    return IsAmideCarbon(m, zatom);
  } else if (a->atomic_number() == 7) {
    return IsAmideNitrogen(m, zatom);
  } else if (a->atomic_number() == 8) {
    return IsAmideOxygen(m, zatom);
  }

  return 0;
}

int
IsSulfonamideNitrogen(const Molecule& m,
                      atom_number_t zatom) {
  std::optional<atom_number_t> sulphur = SinglyBondedTo(m, zatom, 16);
  if (! sulphur) {
    return 0;
  }

  const Atom* a = m.atomi(*sulphur);
  if (a->ncon() != 4) {
    return 0;
  }

  int doubly_bonded_oxygens = 0;
  for (const Bond * b : *a) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(*sulphur);
    if (m.atomic_number(j) == 8) {
      ++doubly_bonded_oxygens;
    }
  }

  return doubly_bonded_oxygens == 2;
}

int
IsSulfonamideSulphur(const Molecule& m,
                     atom_number_t zatom) {
  int doubly_bonded_oxygens = 0;
  int nitrogens = 0;
  for (const Bond* b : m.atom(zatom)) {
    atom_number_t j = b->other(zatom);
    if (b->is_double_bond() && m.atomic_number(j) == 8) {
      ++doubly_bonded_oxygens;
    } else if (b->is_single_bond() && m.atomic_number(j) == 7) {
      ++nitrogens;
    }
  }

  // We also captures sulfamid
  if (doubly_bonded_oxygens == 2 && nitrogens > 0) {
    return 1;
  }

  return 0;
}

int
IsSulfonamideOxygen(const Molecule& m,
                    atom_number_t zatom) {
  std::optional<atom_number_t> sulphur = DoublyBondedTo(m, zatom, 16);
  if (! sulphur) {
    return 0;
  }

  return IsSulfonamideSulphur(m, *sulphur);
}

int IsSulfonamide(const Molecule& m,
                  atom_number_t zatom) {
  const atomic_number_t z = m.atomic_number(zatom);
  if (z == 7) {
    return IsSulfonamideNitrogen(m, zatom);
  }
  if (z == 8) {
    return IsSulfonamideOxygen(m, zatom);
  }
  if (z == 16) {
    return IsSulfonamideSulphur(m, zatom);
  }

  return 0;
}

int
IsNitroNitrogen(const Molecule& m,
                atom_number_t zatom) {
  const Atom& a = m.atom(zatom);
  if (a.ncon() != 3) {
    return 0;
  }

  int doubly_bonded_oxygen = 0;
  for (const Bond* b : a) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (m.atomic_number(j) != 8) {
      return 0;
    }

    ++doubly_bonded_oxygen;
  }

  return doubly_bonded_oxygen == 2;
}

int
IsNitroOxygen(const Molecule& m,
              atom_number_t zatom) {
  const Atom& a = m.atom(zatom);
  if (a.ncon() != 1) {
    return 0;
  }

  const Bond* b = a[0];
  if (! b->is_double_bond()) {
    return 0;
  }

  atom_number_t nitrogen = b->other(zatom);
  if (m.atomic_number(nitrogen) != 7) {
    return 0;
  }

  return IsNitroNitrogen(m, b->other(zatom));
}

int
IsNitro(const Molecule& m,
        atom_number_t zatom) {
  if (m.atomic_number(zatom) == 7) {
    return IsNitroNitrogen(m, zatom);
  }
  if (m.atomic_number(zatom) == 8) {
    return IsNitroOxygen(m, zatom);
  }

  return 0;
}

// Return true if `zatom` is atomic number `target` and
// it is in an aliphatic ring.
int
IsAliphaticRingAtom(Molecule& m,
                const atom_number_t zatom,
                atomic_number_t target) {
  if (m.atomic_number(zatom) != target) {
    return 0;
  }

  // Must be in a ring.
  if (m.ring_bond_count(zatom) == 0) {
    return 0;
  }

  // Must not be aromatic.
  return ! m.is_aromatic(zatom);
}

static Env
IdentifyEnvironment(Molecule& m,
                    const atom_number_t zatom) {
  if (use_terminal_atom_type && m.ncon(zatom) == 0) {
    return Env::kTerminal;
  }

  if (m.is_aromatic(zatom))
    return Env::kAromatic;

  if (IsAmideNitrogen(m, zatom)) {
    return Env::kAmideNitrogen;
  }
  if (IsAmideCarbon(m, zatom)) {
    return Env::kAmideCarbon;
  }
  if (IsSulfonamideNitrogen(m, zatom)) {
    return Env::kAmideNitrogen;
  }
  if (IsSulfonamideSulphur(m, zatom)) {
    return Env::kAmideCarbon;
  }

  if (IsAliphaticRingAtom(m, zatom, 6)) {
    return Env::kAliphaticRingCarbon;
  }

  if (IsAliphaticRingAtom(m, zatom, 7)) {
    return Env::kAliphaticRingNitrogen;
  }

  const Atom * atom = m.atomi(zatom);

  if (atom->is_halogen()) {
    return Env::kHalogen;
  }

  if (IsNitro(m, zatom)) {
    return Env::kNitro;
  }

  if (atom->ncon() == atom->nbonds()) {         // fully saturated
    if (atom->atomic_number() == 6) {
      return Env::kSaturatedCarbon;
    } else if (atom->atomic_number() == 7) {
      return Env::kSaturatedNitrogen;
    } else if (atom->atomic_number() == 8 || atom->atomic_number() == 16) {
      return Env::kSaturatedOxygen;
    }
  } else {
    if (atom->atomic_number() == 6) {
      return Env::kUnsaturatedCarbon;
    } else if (atom->atomic_number() == 7) {
      return Env::kUnsaturatedNitrogen;
    } else if (atom->atomic_number() == 8 || atom->atomic_number() == 16) {
      return Env::kUnsaturatedOxygen;
    }
  }

  return Env::kUndefined;
}

static int
do_add_environment_atoms(Molecule & m,
                         const atom_number_t j1,
                         const atom_number_t j2)
{
  const auto env = static_cast<int>(IdentifyEnvironment(m, j1));
 
  const int ndx = static_cast<int>(env);

  if (ndx == 0) {
    return 0;
  }

  m.add(env_element[ndx]);

  return m.add_bond(m.natoms() - 1, j2, SINGLE_BOND);
}

static int
do_apply_environment_isotopic_labels(Molecule & m,
                                     const atom_number_t j1,
                                     const atom_number_t j2)
{
  const int e = static_cast<int>(IdentifyEnvironment(m, j1));

  if (e > 0) {
    m.set_isotope(j2, e);
  }

  return 1;
}

static void
ApplyAtomTypeIsotopicLabel(Molecule& m,
                const atom_number_t j1,
                const atom_number_t j2,
                const uint32_t* atype) {
  m.set_isotope(j1, atype[j2]);
}

int
set_isotopes_if_needed(Molecule & m,
                       atom_number_t j1,
                       atom_number_t j2,
                       const Dicer_Arguments & dicer_args)
{
  if (isotope_for_join_points > 0 || increment_isotope_for_join_points) {
    do_apply_isotopic_labels(m, j1, j2);

    return 1;
  }

  if (apply_environment_isotopic_labels) {
    do_apply_environment_isotopic_labels(m, j1, j2);
    do_apply_environment_isotopic_labels(m, j2, j1);

    return 1;
  }

  if (apply_atom_type_isotopic_labels) {
    ApplyAtomTypeIsotopicLabel(m, j1, j2, dicer_args.atom_types());
    ApplyAtomTypeIsotopicLabel(m, j2, j1, dicer_args.atom_types());
    return 1;
  }

  if (add_environment_atoms) {
    do_add_environment_atoms(m, j1, j2);
    do_add_environment_atoms(m, j2, j1);

    return 1;
  }

  if (apply_atomic_number_isotopic_labels) {
    do_apply_atomic_number_isotopic_labels(m, j1, j2);
    do_apply_atomic_number_isotopic_labels(m, j2, j1);

    return 1;
  }

  if (isotopic_label_is_recursion_depth) {
    m.set_isotope(j1, dicer_args.recursion_depth());
    m.set_isotope(j2, dicer_args.recursion_depth());

    return 1;
  }

  return 0;         // nothing to do
}

int
Breakages::do_recap(Molecule & m,
                    Dicer_Arguments & dicer_args) const
{
  const int n = _transformations.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Dicer_Transformation * ti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != ti->ttype())
      continue;

    const Bond_Breakage * b = dynamic_cast<const Bond_Breakage *>(ti);

    atom_number_t a0 = b->a1();
    atom_number_t a1 = b->a2();

    if (!m.are_bonded(a0, a1))
      continue;

    m.remove_bond_between_atoms(a0, a1);

    set_isotopes_if_needed(m, a0, a1, dicer_args);
  }


  resizable_array_p<Molecule> c;

  const int nc = m.create_components(c);


  for (int i = 0; i < nc; i++)
  {
    Molecule & ci = *(c[i]);
//  cerr << "checking component " << ci.smiles() <<endl;

    if (!dicer_args.OkAtomCount(ci))
      continue;

    dicer_args.is_unique(ci);
  }

  // Now put things back together. Not sure this is necessary...

  for (int i = 0; i < n; i++)
  {
    const Dicer_Transformation * ti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != ti->ttype())
      continue;

    const Bond_Breakage * b = dynamic_cast<const Bond_Breakage *>(ti);

    atom_number_t a0 = b->a1();
    atom_number_t a1 = b->a2();

    if (m.are_bonded(a0, a1))
      continue;

    m.add_bond(a0, a1, SINGLE_BOND);

    m.set_isotope(a0, 0);
    m.set_isotope(a1, 0);
  }

  return 1;
}

int
Breakages::identify_symmetry_exclusions(Molecule & m)
{
  const int * s = m.symmetry_classes();

  int n = _transformations.number_elements();

  for (int i = 0; i < n; i++)
  {
    Dicer_Transformation * ti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != ti->ttype())
      continue;

    Chain_Bond_Breakage * b = dynamic_cast<Chain_Bond_Breakage *>(ti);

    if (nullptr == b)
    {
      cerr << "Skipping nullptr chain bond breakage ptr\n";
      continue;
    }

    if (b->symmetry_equivalent())
      continue;

    int s0 = s[b->a1()];
    int s1 = s[b->a2()];

    for (int j = i + 1; j < n; j++)
    {
      Dicer_Transformation * tj = _transformations[j];

      if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != tj->ttype())
        continue;

      Chain_Bond_Breakage * bj = dynamic_cast<Chain_Bond_Breakage *>(tj);
      if (nullptr == bj)
      {
        cerr << "Skipping nullptr ptr chain bond breakage\n";
        continue;
      }

      if (bj->symmetry_equivalent())
        continue;

      if (s[bj->a1()] == s0 && s[bj->a2()] == s1)
        bj->set_symmetry_equivalent(1);
      else if (s[bj->a2()] == s0 && s[bj->a1()] == s1)
        bj->set_symmetry_equivalent(1);
    }
  }

  return 1;
}


int
Dicer_Arguments::_is_smallest_fragment (const IW_Bits_Base & b) const
{
  int n = _stored_bit_vectors.number_elements();

  for (int i =0; i<n; i++)
  {
    if (b == *(_stored_bit_vectors[i]))
      continue;

    if (b.is_subset(*_stored_bit_vectors[i]))
      return 0;
  }

  return 1;
}

int
Dicer_Arguments::write_fragments_and_complements(Molecule & m,
                                                 DicerFragmentOutput& output) const
{
  if (output.output_is_proto()) {
    return WriteFragmentsAndComplementsProto(m, output);
  }

  int * tmp1 = new int[_nbits]; std::unique_ptr<int[]> free_tmp1(tmp1);

  int n = _stored_bit_vectors.number_elements();

  //cerr << "Writing " << n << " stored bit vectors\n";

  for (int i = 0; i < n; i++)
  {
    const IW_Bits_Base * b = _stored_bit_vectors[i];

    assert (b->nbits() > 0);

    if (write_only_smallest_fragment && (!_is_smallest_fragment(*b)))
      continue;

    _write_fragment_and_complement(m, *b, tmp1, output.text_output());
  }

  if (0 == lower_atom_count_cutoff())
    _write_fragment_and_hydrogen(m, tmp1, output.text_output());

  return 1;
}

int
Dicer_Arguments::WriteFragmentsAndComplementsProto(Molecule & m,
                                                  DicerFragmentOutput & output) const {

  std::unique_ptr<int[]> tmp1 = std::make_unique<int[]>(_nbits);

  dicer_data::DicedMolecule proto;

  proto.set_name(m.name().data(), m.name().length());
  const IWString& usmi = m.unique_smiles();
  proto.set_smiles(usmi.data(), usmi.length());
  proto.set_natoms(m.natoms());
  // We do not know the number of breakable bonds here.
  // proto.set_nbonds(breakable_bonds);

  for (const IW_Bits_Base* b : _stored_bit_vectors) {
    assert (b->nbits() > 0);

    if (write_only_smallest_fragment && (!_is_smallest_fragment(*b))) {
      continue;
    }

    AppendFragmentAndComplement(m, *b, tmp1.get(), proto);
  }

  return output.Write(proto);
}

static void
dump_iw_atom_type (Molecule& frag, const Molecule& comp,
                   const int parent_atoms, const int xref[], const uint32_t atypes[],
                   IWString_and_File_Descriptor & output)
{
//Molecule mcopy(comp);
//cerr << "dump_iw_atom_type  processing " << frag.smiles() << " and comp " << mcopy.smiles() << endl;
  IWString v[10];
  for (int i = 0; i < parent_atoms; ++i)
  {
    if (xref[i] < 0)
      continue;

    const uint32_t atype = atypes[i];

//  cerr << " i = " << i << " xref " << xref[i] << " iso " << comp.isotope(xref[i]) << " atype " << atypes[i] << endl;
    for (int iso = comp.isotope(xref[i]); iso > 0; iso /= 10)
    {
      int breakage = iso % 10;

      if (v[breakage].length())
        v[breakage] << ',';
      v[breakage] << atype;
    }
  }

  const int matoms = frag.natoms();
  output << " AT=[";

  frag.unique_smiles();            // force unique smiles

  const resizable_array<atom_number_t> & atom_order_in_smiles = frag.atom_order_in_smiles();

  assert (atom_order_in_smiles.number_elements() == matoms);

  int count = 0;
  for (int i = 0; i < matoms; ++i)
  {
    const atom_number_t j = atom_order_in_smiles[i];

    int iso = frag.isotope(j);

//  cerr << " i = " << i << " atom " << j << " iso " << iso << endl;

    for (; iso > 0; iso /= 10) {
      if (count > 0)
        output << "|";

      count++;

      const int breakage = iso % 10;
      if (0 == breakage || 0 == v[breakage].length())
        output << breakage << ":" << "INVALID";
      else
        output << breakage << ":" << v[breakage];
    }
  }

  output << "]";

  return;
}

static void
dump_atom_type (Molecule& frag, const Molecule& comp,
                const int parent_atoms, const int xref[], const uint32_t atypes[],
                IWString_and_File_Descriptor & output)
{
  if (atypes == nullptr) {
    return;
  }

  if (ISO_CF_ATYPE == apply_isotopes_to_complementary_fragments)
  {
    dump_iw_atom_type(frag, comp, parent_atoms, xref, atypes, output);
    return;
  }

//static const char* UNKNOWN_TYPE = "X";
  static const char* SYBYL_TYPES[] = {
    "H", "C3", "C2", "CAR", "C1", "N4", "N3", "N2", "N1", "NAR", "NPL3","NAM", "O3", "O2", "S3", "S2", "P3",
    "BR", "CL", "F", "I", "SO", "SO2", "OCO2", "CCAT", "X"
  };

  static int NUMBER_OF_SYBYL_TYPES = sizeof(SYBYL_TYPES) / sizeof(SYBYL_TYPES[0]) - 1;

  IWString v[10];
  for (int i = 0; i < parent_atoms; ++i) {
    if (xref[i] != -1) {
      const Atom& atom_i = comp.atom(xref[i]);
      uint32_t atype = atypes[i];
      if (atype > static_cast<uint32_t>(NUMBER_OF_SYBYL_TYPES)) {
        atype = NUMBER_OF_SYBYL_TYPES;
      }
      for (int iso = atom_i.isotope(); iso > 0; iso /= 10) {
        int breakage = iso % 10;
        v[breakage].append_with_spacer(SYBYL_TYPES[atype], ',');
      }
    }
  }

  const int matoms = frag.natoms();
  output << " AT=[";

  frag.unique_smiles();            // force unique smiles

  const resizable_array<atom_number_t> & atom_order_in_smiles = frag.atom_order_in_smiles();

  assert (atom_order_in_smiles.number_elements() == matoms);

  int count = 0;
  for (int i = 0; i < matoms; ++i)
  {
    const atom_number_t j = atom_order_in_smiles[i];

    const Atom& atom_j = frag.atom(j);

    for (int iso = atom_j.isotope(); iso > 0; iso /= 10) {
      if (count > 0)
        output << "|";

      count++;

      int breakage = iso % 10;
      if (breakage == 0 || v[breakage] == "" || v[breakage].c_str() == nullptr)
        output << breakage << ":" << "INVALID";
      else
        output << breakage << ":" << v[breakage];
    }
  }

  output << "]";
}

/*
  Create complementary fragment for every atom that has a Hydrogen attachment
*/

int
Dicer_Arguments::_write_fragment_and_hydrogen(Molecule & m,
                                              int * subset_specification,
                                              IWString_and_File_Descriptor & output) const
{
  const int matoms = m.natoms();
  std::fill_n(subset_specification, matoms, 1);
  for (int i = 0; i < matoms; i++)
  {
    if (0 == m.hcount(i))
      continue;

    output << '[';

    Molecule s1;
    m.create_subset(s1, subset_specification, 1, _local_xref);
    if (apply_isotopes_to_complementary_fragments > 0) {
      s1.set_isotope(_local_xref[i], apply_isotopes_to_complementary_fragments);
      output << apply_isotopes_to_complementary_fragments << "H]";
    }
    else if (apply_isotopes_to_complementary_fragments < 0)
    {
      s1.set_isotope(_local_xref[i], 1);
      output << "1H]";
    }

    output << ' ' << m.name();

    if (append_atom_count)
      output << " 1 " << 1.0f / static_cast<float>(m.natoms());

    if (apply_isotopes_to_complementary_fragments < 0)
      dump_atom_type(s1, s1, matoms, _local_xref, _atom_type, output);

    output << ' ';
    output << s1.unique_smiles() << " COMP " << m.name();

    if (nullptr == _atom_type)
      ;
    else if (apply_isotopes_to_complementary_fragments < 0)
      output << " AT=[1:H]";

    if (append_atom_count)
      output << ' ' << s1.natoms() << ' ' << (static_cast<float>(s1.natoms()) / static_cast<float>(m.natoms()));

    output << '\n';
  }

  return 1;
}


int
Dicer_Arguments::_write_fragment_and_complement (Molecule & m,
                                  const IW_Bits_Base & b,
                                  int * subset_specification,
                                  IWString_and_File_Descriptor & output) const
{
  const int matoms = m.natoms();

  assert (b.nbits() >= matoms);

  if (0 == b.nset() || matoms == b.nset())
    return 1;

  if (! ok_atom_count(b.nset()))
    return 0;

  std::fill_n(subset_specification, matoms, 0);

  b.set_vector(subset_specification);

  //cerr << "Nset " <<b.nset() << ", natoms " << matoms <<endl;
  //cerr << "First has " << count_non_zero_occurrences_in_array(subset_specification, matoms) << " atoms\n";

#ifdef THIS_CANNOT_HAPPEN
  if (0 == count_non_zero_occurrences_in_array(subset_specification, matoms))
  {
    cerr << "NO atoms in molecule\n";
    b.printon(cerr);
    cerr << endl;
    for (int i = 0; i < matoms; i++)
    {
      cerr << subset_specification[i];
    }
    cerr << endl;
  }
#endif

  Molecule s1;
  m.create_subset(s1, subset_specification, 1, _local_xref);

  Molecule s2;
  m.create_subset(s2, subset_specification, 0, _local_xref2);

//cerr << "Subsets " << s1.smiles() << " and " << s2.smiles() << " apply_isotopes_to_complementary_fragments " << apply_isotopes_to_complementary_fragments << endl;

  if (apply_isotopes_to_complementary_fragments > 0) {
    _do_apply_isotopes_to_complementary_fragments(m, s1, _local_xref);
    _do_apply_isotopes_to_complementary_fragments(m, s2, _local_xref2);
  }
  else if (apply_isotopes_to_complementary_fragments < 0)
    _do_apply_isotopes_to_complementary_fragments(m, s1, s2);

//cerr << "LINE " << __LINE__ << " have " << s1.smiles() << " and " << s2.smiles() << endl;

  output << s1.unique_smiles() << ' ' << m.name();

  if (append_atom_count)
    output << ' ' << s1.natoms() << ' ' << (static_cast<float>(s1.natoms())/static_cast<float>(matoms));

  if (apply_isotopes_to_complementary_fragments == -1)
    dump_atom_type(s1, s2, matoms, _local_xref2, _atom_type, output);

  output << ' ' << s2.unique_smiles() << " COMP " << m.name();

  if (apply_isotopes_to_complementary_fragments < 0)
    dump_atom_type(s2, s1, matoms, _local_xref, _atom_type, output);

  if (append_atom_count)
    output << ' ' << s2.natoms() << ' ' << (static_cast<float>(s2.natoms())/static_cast<float>(matoms));

  output << '\n';

  return 1;
}

// `b` describes a subset of `m`. Generate that subset and the complement
// fragment and load those into `proto`.
// subset_specification is just used as working storate - convert `b`
// to an array.
// Maybe should try to share code with _write_fragment_and_complement
int
Dicer_Arguments::AppendFragmentAndComplement(Molecule& m,
                                  const IW_Bits_Base & b,
                                  int * subset_specification,
                                  dicer_data::DicedMolecule& proto) const {
  const int matoms = m.natoms();

  assert (b.nbits() >= matoms);

  if (0 == b.nset() || matoms == b.nset()) {
    return 1;
  }

  if (! ok_atom_count(b.nset())) {
    return 0;
  }

  std::fill_n(subset_specification, matoms, 0);

  b.set_vector(subset_specification);

  //cerr << "Nset " <<b.nset() << ", natoms " << matoms <<endl;
  //cerr << "First has " << count_non_zero_occurrences_in_array(subset_specification, matoms) << " atoms\n";

  Molecule s1;
  m.create_subset(s1, subset_specification, 1, _local_xref);

  Molecule s2;
  m.create_subset(s2, subset_specification, 0, _local_xref2);

  if (upper_atom_count_cutoff_complementary_fragment > 0 &&
      s2.natoms() > upper_atom_count_cutoff_complementary_fragment) {
    return 1;
  }

  if (apply_isotopes_to_complementary_fragments > 0) {
    _do_apply_isotopes_to_complementary_fragments(m, s1, _local_xref);
    _do_apply_isotopes_to_complementary_fragments(m, s2, _local_xref2);
  }
  else if (apply_isotopes_to_complementary_fragments < 0)
    _do_apply_isotopes_to_complementary_fragments(m, s1, s2);

  dicer_data::DicerFragment* frag = proto.add_fragment();
  const IWString& usmi1 = s1.unique_smiles();
  frag->set_smi(usmi1.data(), usmi1.length());
  frag->set_nat(s1.natoms());
  const IWString& usmi2 = s2.unique_smiles();
  frag->set_comp(usmi2.data(), usmi2.length());

  // We do not know the ID of this fragment, and we don't care about the number of instances in the molecule.

  return 1;
}

/*
   We have generated two complementary fragments. The _local_xref array
   contains the cross reference from initial atom numbers
 */

#ifdef VERSION_USING_ATOMS_RATHER_THAN_BONDS
int
Dicer_Arguments::_do_apply_isotopes_to_complementary_fragments (Molecule & parent,
                                                Molecule & subset, const int xref[]) const
{

 const int matoms = parent.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (xref[i] < 0)               // atom I not part of subset
      continue;

    const Atom * a = parent.atomi(i);

    const int acon = a->ncon();

    for (int j = 0; j < acon; j++)               // are any atoms bonded to I not in the subset
    {
      atom_number_t k = a->other(i, j);

      if (xref[k] >= 0)     // is in the subset
        continue;

//    cerr << "Setting isotope " << xref[i] << " to " << apply_isotopes_to_complementary_fragments << endl;
      subset.set_isotope(xref[i], apply_isotopes_to_complementary_fragments);
      break;
    }
  }

  return 1;
}
#endif

int
Dicer_Arguments::_do_apply_isotopes_to_complementary_fragments (Molecule & parent,
                                                Molecule & subset, const int xref[]) const
{
  const int nb = parent.nedges();

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = parent.bondi(i);

    const atom_number_t xrefa1 = xref[b->a1()];    // atom number of atom in subset
    const atom_number_t xrefa2 = xref[b->a2()];    // atom number of atom in subset

    if (xrefa1 < 0 && xrefa2 < 0)    // neither side of bond is in subset
      continue;

    if (xrefa1 >= 0 && xrefa2 >= 0)   // both sides of bond are in subset
      continue;

    if (xrefa1 >= 0)    // a1 in subset
      subset.set_isotope(xrefa1, apply_isotopes_to_complementary_fragments);
    else            // if (xref[a2] >= 0) must be true
      subset.set_isotope(xrefa2, apply_isotopes_to_complementary_fragments);
  }

  return 1;
}

/*
  _local_xref generated FRAGMENT and _local_xref2 generated COMP
*/

int
Dicer_Arguments::_do_apply_isotopes_to_complementary_fragments(Molecule & parent,
                                                        Molecule & fragment,
                                                        Molecule & comp) const
{
  const int matoms = parent.natoms();

  int n = 0;

// First identify those atoms that are at a break point

  int * labels = new_int(matoms); std::unique_ptr<int[]> free_labels(labels);

  const int nbonds = parent.nedges();

  n = 0;
  for (int i = 0; i < nbonds; ++i)
  {
    const Bond * b = parent.bondi(i);

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    const atom_number_t xref2a1 = _local_xref2[a1];
    const atom_number_t xref2a2 = _local_xref2[a2];

    if (xref2a1 >= 0 && xref2a2 < 0)    // atom a1 is in COMP, atom a2 is not
    {
      comp.set_isotope(xref2a1, ++n);
      labels[a1] = n;
    }
    else if (xref2a1 < 0 && xref2a2 >= 0)
    {
      comp.set_isotope(xref2a2, ++n);
      labels[a2] = n;
    }
  }

//cerr << "Line " << __LINE__ << " n = " << n << endl;
  if (n > 1)
    return _do_apply_isotopes_to_complementary_fragments_canonical(parent, n, fragment, comp);

// The easier case of just one break detected
// Put isotopes on FRAGMENT

//#define OLD_VERSION_OQUWEQWE
#ifdef OLD_VERSION_OQUWEQWE
  for (int i = 0; i < matoms; i++)
  {
    if (_local_xref[i] < 0)               // atom I not part of frag
      continue;

    const Atom * a = parent.atomi(i);

    const int acon = a->ncon();
    int v0 = 0;

    for (int j = 0; j < acon; j++)               // are any atoms bonded to I not in the subset
    {
      atom_number_t k = a->other(i, j);

      if (labels[k] > 0) {
        v0 = v0*10+labels[k];
      }
    }
    fragment.set_isotope(_local_xref[i], v0);
  }
#endif


#define NEW_VERSION_QWEQWEQ
#ifdef NEW_VERSION_QWEQWEQ
  for (int i = 0; i < nbonds; ++i)
  {
    const Bond * b = parent.bondi(i);

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    const atom_number_t xref1a1 = _local_xref[a1];
    const atom_number_t xref1a2 = _local_xref[a2];

    if (xref1a1 >= 0 && xref1a2 < 0)    // atom a1 is in FRAGMENT, atom a2 is in COMP
      fragment.set_isotope(xref1a1, 10 * fragment.isotope(xref1a1) + labels[a2]);
    else if (xref1a1 < 0 && xref1a2 >= 0)
      fragment.set_isotope(xref1a2, 10 * fragment.isotope(xref1a2) + labels[a1]);
  }
#endif

  return 1;
}


/*
  There are multiple join points in the complementary fragment. We need to get canonical
  isotope labelling.
*/

class Atom_Rank_Atom
{
  private:
    atom_number_t _in_frag;
    int _canon;       // in parent
    atom_number_t _in_comp;

  public:
    Atom_Rank_Atom () {};

    void set (atom_number_t af, int c, atom_number_t ac)
    {
      _in_frag = af;
      _canon = c;
      _in_comp = ac;
    }

    int canon () const { return _canon;}
    atom_number_t in_frag () const { return _in_frag;}
    atom_number_t in_comp () const { return _in_comp;}
};

static void
reset_isotope_multiple_connections (Molecule & m,
                                    const atom_number_t zatom,
                                    const int * isotope_translation_table)
{
  int iso = m.isotope(zatom);

  resizable_array<int>attached;

  while (iso > 0)
  {
    int i = iso % 10;
    attached.add(isotope_translation_table[i]);
    iso = iso / 10;
  }

  Int_Comparator_Larger icl;
  attached.iwqsort(icl);

//attached.iwqsort_lambda([] (int i1, int i2) { if (i1 < i2) return -1; if (i1 > i2) return 1; return 0;});

  int new_iso = 0;

  for (unsigned int i = 0; i < attached.size(); ++i)
  {
    new_iso = 10 * new_iso + attached[i];
  }

  m.set_isotope_no_perturb_canonical_ordering(zatom, new_iso);

  return;
}

class Atom_Rank_Atom_Comparator
{
  private:
  public:
    int operator () (const Atom_Rank_Atom & ara1, const Atom_Rank_Atom & ara2) const
    {
      return ara1.canon() < ara2.canon();
    }
};

/*
  This is really strange. We need to ensure that the isotopic atoms in the
  complementary fragment are in sequential ascending order.
*/

int
Dicer_Arguments::_do_apply_isotopes_to_complementary_fragments_canonical(Molecule & parent,
                                                  int n,
                                                  Molecule & frag, Molecule& comp) const
{
  const int matoms = parent.natoms();

  Atom_Rank_Atom * ara = new Atom_Rank_Atom[n]; std::unique_ptr<Atom_Rank_Atom[]> free_ara(ara);

  int ndx = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (_local_xref2[i] < 0)               // atom I not part of comp
      continue;

    const Atom * a = parent.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)               // are any atoms bonded to I not in comp
    {
      atom_number_t k = a->other(i, j);

      if (_local_xref2[k] < 0) {
        ara[ndx].set(_local_xref[k], parent.canonical_rank(i), _local_xref2[i]);
        ndx++;
        break;
      }
    }
  }

#ifdef IW_HAVE_LAMBDA
  std::sort (ara, ara + n, [] (const Atom_Rank_Atom & ara1, const Atom_Rank_Atom & ara2) { return ara1.canon() < ara2.canon();});
#else
  Atom_Rank_Atom_Comparator arac;
  std::sort(ara, ara + n, arac);
#endif

  assert (n == ndx);

  for (int i = 0; i < n; ++i)
  {
    const Atom_Rank_Atom & arai = ara[i];

//  cerr << " i = " << i << " rand " << arai.canon() << endl;

    comp.set_isotope(arai.in_comp(), i + 1);

    const atom_number_t af = arai.in_frag();

    int iso = frag.isotope(af) * 10 + i + 1;

    frag.set_isotope(af, iso);
  }

  comp.unique_smiles();    // force smiles generation

#ifdef DEBUG_DO_APPLY_ISOTOPES_TO_COMPLEMENTARY_FRAGMENTS_CANONICAL
  cerr << "_do_apply_isotopes_to_complementary_fragments_canonical:processing " << comp.unique_smiles() << endl;
#endif

  const resizable_array<atom_number_t> & atom_order_in_smiles = comp.atom_order_in_smiles();

  const int cmatoms = comp.natoms();

  int prev_isotope = 0;

  int in_order = 1;

  for (int i = 0; i < cmatoms; ++i)
  {
    int j = atom_order_in_smiles[i];

    const int iso = comp.isotope(j);

    if (0 == iso)
      continue;

    if (iso <= prev_isotope)
    {
      in_order = 0;
      break;
    }

    prev_isotope = iso;
  }

#ifdef DEBUG_DO_APPLY_ISOTOPES_TO_COMPLEMENTARY_FRAGMENTS_CANONICAL
  cerr << "in_order? " << in_order << endl;
#endif

  if (in_order)     // great!
    return 1;

  int * isotope_translation_table = new int[n+1]; std::unique_ptr<int[]> free_isotope_translation_table(isotope_translation_table);

  ndx = 1;

  for (int i = 0; i < cmatoms; ++i)
  {
    atom_number_t j = atom_order_in_smiles[i];

    const int iso = comp.isotope(j);

    if (0 == iso)
      continue;

    isotope_translation_table[iso] = ndx;
#ifdef DEBUG_DO_APPLY_ISOTOPES_TO_COMPLEMENTARY_FRAGMENTS_CANONICAL
    cerr << "Will translate isotope " << iso << " to " << ndx << endl;
#endif
    ndx++;
  }

  for (int i = 0; i < cmatoms; ++i)
  {
    const int iso = comp.isotope(i);

    if (0 == iso)
      continue;

    if (iso != isotope_translation_table[iso])
      comp.set_isotope_no_perturb_canonical_ordering(i, isotope_translation_table[iso]);
  }

#ifdef DEBUG_DO_APPLY_ISOTOPES_TO_COMPLEMENTARY_FRAGMENTS_CANONICAL
  cerr << " comp unique smiles " << comp.unique_smiles() << endl;
#endif

  const int fmatoms = frag.natoms();

  for (int i = 0; i < fmatoms; ++i)
  {
    const int iso = frag.isotope(i);

    if (0 == iso)
      continue;

    if (iso > 10)   // special handling, complicated
      reset_isotope_multiple_connections(frag, i, isotope_translation_table);
    else if (iso != isotope_translation_table[iso])
      frag.set_isotope_no_perturb_canonical_ordering(i, isotope_translation_table[iso]);
  }

  return 1;
}


int
OpenStreaForAllFragments(IWString& fname,
                         IWString_and_File_Descriptor& stream_for_all_fragments) {
  return stream_for_all_fragments.open(fname.null_terminated_chars());
}

template <typename T>
void
AppendMaybeWithSeparator(bool& need_separator,
                        char sep,
                        const T& s,
                        IWString_and_File_Descriptor& output) {

  if (need_separator) {
    output << sep;
  } else {
    need_separator = true;
  }
  output << s;
}

int
Dicer_Arguments::WriteAllFragments(Molecule& parent,
                                   Molecule& fragment,
                                   AllFragmentsStream& all_fragments_stream) const {
  if (! all_fragments_stream.active()) {
    return 1;
  }

  auto iter_frag = smiles_to_id.find(fragment.unique_smiles());

  if (iter_frag == smiles_to_id.end()) {
    // Should complain here.
    return 0;
  }

  bool need_separator = false;
  all_fragments_stream.MaybeWrite("smiles", _current_molecule->unique_smiles(), need_separator);
  all_fragments_stream.MaybeWrite("name", _current_molecule->name(), need_separator);
  all_fragments_stream.MaybeWrite("rdepth", _recursion_depth, need_separator);
  all_fragments_stream.MaybeWrite("psmiles", parent.unique_smiles(), need_separator);
  all_fragments_stream.MaybeWrite("frag", fragment.unique_smiles(), need_separator);
  all_fragments_stream.MaybeWrite("frag_id", iter_frag->second, need_separator);

  return all_fragments_stream.newline();
}

int
WriteAllFragments(Molecule& parent,
                  Molecule& fragment,
                  IWString_and_File_Descriptor& stream_for_all_fragments) {
  if (! stream_for_all_fragments.active()) {
    return 1;
  }

  auto iter_frag = smiles_to_id.find(fragment.unique_smiles());

  if (iter_frag == smiles_to_id.end()) {
    // Should complain here.
    return 0;
  }

  constexpr char kSep = ' ';

  stream_for_all_fragments << parent.unique_smiles() << kSep << parent.name() << kSep <<
                               fragment.unique_smiles() << kSep << iter_frag->second << '\n';

  stream_for_all_fragments.write_if_buffer_holds_more_than(4096);

  return stream_for_all_fragments.good();
}

/*
   We need to group together bonds that are missing
 */


static int
write_hash_data(IWString_and_File_Descriptor & output,
                const IW_STL_Hash_Map_uint & smiles_to_bit)
{
  if (verbose)
    cerr << "Writing " << smiles_to_bit.size() << " smiles->bit relationships\n";

  for (IW_STL_Hash_Map_uint::const_iterator i = smiles_to_bit.begin(); i != smiles_to_bit.end(); ++i)
  {
    const IWString & smi = (*i).first;
    unsigned int b = (*i).second;

    output << smi << " FRAGID " << b << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

// Populate `proto` with information about a fragment `smiles` found `count`
// times in the current molecule.
void
ToProto(const IWString& smiles,
        uint32_t count,
        const IW_STL_Hash_Map_String& starting_parent,
        dicer_data::DicerFragment& proto) {
  proto.set_smi(smiles.AsString());
  proto.set_n(count);
  proto.set_nat(count_atoms_in_smiles(smiles));
  if (auto f = starting_parent.find(smiles);
      f != starting_parent.end()) {
    proto.set_par(f->second.AsString());
  }

  if (isotope_for_join_points) {
    proto.set_iso(dicer_data::ATT);
  } else if (apply_atomic_number_isotopic_labels) {
    proto.set_iso(dicer_data::Z);
  } else if (apply_environment_isotopic_labels) {
    proto.set_iso(dicer_data::ENV);
  }
}

static int
WriteFragmentStatisticsAsProto(const IWString& smiles,
                               uint32_t count,
                               IWString_and_File_Descriptor& output) {
  dicer_data::DicerFragment proto;
  ToProto(smiles, count, starting_parent, proto);

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (! printer.PrintToString(proto, &buffer)) {
    cerr << "WriteFragmentStatisticsAsProto:cannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  if (prepend_smiles_to_fragstat_proto) {
    output << smiles << ' ';
  }

  if (buffer.ends_with(' ')) {
    buffer.pop_back();
  }

  output << buffer;
  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
WriteFragmentStatisticsAsBinaryProto(const IWString& fragment_smiles,
                                     const uint32_t count,
                                     iw_tf_data_record::TFDataWriter& output) {
  dicer_data::DicerFragment proto;
  ToProto(fragment_smiles, count, starting_parent, proto);

  return output.WriteSerializedProto<dicer_data::DicerFragment>(proto);
}

static void
ReportSuppressed(const IW_STL_Hash_Map_uint& molecules_containing,
                 uint32_t suppressed_by_support_requirement,
                 std::ostream& output) {
  float fraction = iwmisc::Fraction<float>(suppressed_by_support_requirement,
                        static_cast<uint32_t>(molecules_containing.size()));

  output << "Suppressed " << suppressed_by_support_requirement << " of "
            << molecules_containing.size() <<
            " fragments (" << fraction << ") because of support requirement\n";
}

static int
WriteFragmentStatisticsAsBinaryProto(const IW_STL_Hash_Map_uint & molecules_containing,
                          Fragment_Statistics_Report_Support_Specification & fsrss,
                          iw_tf_data_record::TFDataWriter& output) {
  fsrss.set_population_size(molecules_read);

  uint32_t suppressed_by_support_requirement = 0;

  for (const auto& [s, count] : molecules_containing) {
    if (!fsrss.meets_support_requirement(count)) {
      suppressed_by_support_requirement++;
      continue;
    }

    WriteFragmentStatisticsAsBinaryProto(s, count, output);
  }

  if (verbose && fsrss.initialised()) {
    ReportSuppressed(molecules_containing, suppressed_by_support_requirement, cerr);
  }

  return 1;
}

static int
write_fragment_statistics(const IW_STL_Hash_Map_uint & molecules_containing,
                          Fragment_Statistics_Report_Support_Specification & fsrss,
                          IWString_and_File_Descriptor & output)
{
  fsrss.set_population_size(molecules_read);

  //float float_molecules_read = static_cast<float>(molecules_read);

  uint32_t suppressed_by_support_requirement = 0;

  for (const auto& [s, count] : molecules_containing) {
    if (!fsrss.meets_support_requirement(count)) {
      suppressed_by_support_requirement++;
      continue;
    }

    if (write_fragstat_as_binary_proto) {
      WriteFragmentStatisticsAsBinaryProto(s, count, tfdata_writer);
      continue;
    } else if (write_fragstat_as_proto) {
      WriteFragmentStatisticsAsProto(s, count, output);
      continue;
    }

    output << s << " FRAGID=" << smiles_to_id[s];

    if (append_atom_count) {
      output << " natoms " << count_atoms_in_smiles(s);
    }

    output << " in " << count << " molecules";

    //  Nice idea, but it makes the output too messy

    //  if (fsrss.fractional_support_level_specified())
    //  {
    //    set_default_iwstring_float_concatenation_precision(3);
    //    output << ' ' << static_cast<float>(count) / float_molecules_read;
    //    set_default_iwstring_float_concatenation_precision(7);
    //  }

    if (accumulate_starting_parent_information) {
      IW_STL_Hash_Map_String::const_iterator f = starting_parent.find(s);
      if (f == starting_parent.end())
        cerr << "Huh, no starting parent for '" << s << "'\n";
      else
        output << " FP " << (*f).second;

      IW_STL_Hash_Map_int::const_iterator f2 = atoms_in_parent.find((*f).second);
      if (f2 == atoms_in_parent.end())
        cerr << "Huh, no atoms for parent '" << (*f).second << "', frag " << s << "'\n";
      else
        output << " PNA " << (*f2).second;
    }

    output << "\n";

    if (count > static_cast<unsigned int>(molecules_read))
      cerr << "Warning, out of range, " << molecules_read << " molecules read\n";

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  if (verbose && fsrss.initialised()) {
    ReportSuppressed(molecules_containing, suppressed_by_support_requirement, cerr);
  }

  return 1;
}

// Printable output, either free text or text proto.
static int
write_fragment_statistics(const IW_STL_Hash_Map_uint & molecules_containing,
                          Fragment_Statistics_Report_Support_Specification & fsrss,
                          IWString & name_for_fragment_statistics_file)
{
  IWString_and_File_Descriptor output;

  if (!output.open(name_for_fragment_statistics_file.null_terminated_chars()))
  {
    cerr << "Cannot open fragment statistics file '" << name_for_fragment_statistics_file << "'\n";
    return 0;
  }

  return write_fragment_statistics(molecules_containing, fsrss, output);
}

static int
WriteFragmentStatisticsAsBinaryProto(const IW_STL_Hash_Map_uint& molecules_containing,
                Fragment_Statistics_Report_Support_Specification& fsrss,
                IWString& fname) {
  iw_tf_data_record::TFDataWriter output;
  if (! output.Open(fname)) {
    cerr << "WriteFragmentStatisticsAsBinaryProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteFragmentStatisticsAsBinaryProto(molecules_containing, fsrss, output);
}

static int
read_existing_hash_data_record (const const_IWSubstring & buffer,
                                IW_STL_Hash_Map_uint & smiles_to_bit)
{
  if (buffer.nwords() < 3)
  {
    cerr << "Existing dicer data must have at least 3 tokens on each line\n";
    return 0;
  }

  IWString * smi = new IWString;

  int i = 0;

  buffer.nextword(*smi, i);

  const_IWSubstring token;

  buffer.nextword(token, i);

  if ("FRAGID" != token)
    cerr << "Warning, 2nd token should be 'FRAGID', not '" << token << "'\n";

  buffer.nextword(token, i);

  unsigned int b;
  if (!token.numeric_value(b))
  {
    cerr << "Invalid bit number '" << token << "'\n";
    return 0;
  }

  smiles_to_bit[*smi] = b;

  id_to_smiles.add(smi);

  return 1;
}

static int
read_existing_hash_data(iwstring_data_source & input,
                        IW_STL_Hash_Map_uint & smiles_to_bit)
{
  const int nr = input.records_remaining();
  if (0 == nr)
    return 0;

  id_to_smiles.resize(nr);

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (0 == buffer.length() || buffer.starts_with('#'))
      continue;

    if (!read_existing_hash_data_record(buffer, smiles_to_bit))
    {
      cerr << "Invalid hash data '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }
  }

  return smiles_to_bit.size();
}

static int
read_existing_hash_data (const IWString & fname,
                        IW_STL_Hash_Map_uint & smiles_to_bit)
{
  iwstring_data_source input(fname);

  if (!input.good())
  {
    cerr << "Cannot open existing smiles->bit cross reference file '" << fname << "'\n";
    return 0;
  }

  return read_existing_hash_data(input, smiles_to_bit);
}

static int
do_write_to_stream_for_post_breakage(Molecule & m,
                                     IWString_and_File_Descriptor & stream_for_post_breakage)
{
  stream_for_post_breakage << m.unique_smiles() << ' ' << m.name() << '\n';
  stream_for_post_breakage.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
convert_to_atom_numbers_in_current_molecule(const Molecule & m,
                                            const atom_number_t a1,
                                            atom_number_t & t1,
                                            const atom_number_t a2,
                                            atom_number_t & t2)
{
  t1 = INVALID_ATOM_NUMBER;
  t2 = INVALID_ATOM_NUMBER;

  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Atom * a = m.atomi(i);

    const void * v = a->user_specified_void_ptr();
    if (nullptr == v)
      continue;

    atom_number_t x = *(reinterpret_cast<const atom_number_t *>(v));
    //  cerr << "Atom " << i << " was atom " << x << endl;

    if (x == a1)
    {
      t1 = i;
      if (INVALID_ATOM_NUMBER != t2)
        return 1;
    }
    else if (x == a2)
    {
      t2 = i;
      if (INVALID_ATOM_NUMBER != t1)
        return 1;
    }
  }

  //cerr << "Gack, no match for atom numbers " << a1 << " (" << t1 << ") and/or " << a2 << " (" << t2 << ")\n";

  return 0;
}

#ifdef NOT_NEEDED_HERE
atom_number_t
atom_number_in_current_molecule (const Molecule & m,
  const int * xref,
  atom_number_t a)
{
  if (xref[a] < 0)
    return INVALID_ATOM_NUMBER;

  return xref[a];
}
#endif

static int
dicer (Molecule & m, Dicer_Arguments & dicer_args, Breakages & breakages, Breakages_Iterator bi);



static int
call_dicer (Molecule & m, Dicer_Arguments & dicer_args, Breakages & breakages,
            Breakages_Iterator & bi)
{
#ifdef DEBUG_CALL_DICER
  cerr << "Calling dicer, depth " << dicer_args.recursion_depth() << ", molecule is '" << m.smiles() << "', can_go_deeper " << dicer_args.can_go_deeper() << endl;
  cerr << "max_recursion_depth " << max_recursion_depth << " max recursion_depth this m " << dicer_args.max_recursion_depth_this_molecule() << endl;
#endif

  if (!dicer_args.can_go_deeper())
    return 0;

#ifdef DEBUG_CALL_DICER
  cerr << "Compare " << m.natoms()  << " with " << dicer_args.lower_atom_count_cutoff() << endl;
#endif

  if (m.natoms() <= dicer_args.lower_atom_count_cutoff())         // cannot be split
    return 0;

  //dicer_args.increment_recursion_depth();
  //bi.increment();

  dicer(m, dicer_args, breakages, bi);

  //bi.decrement();
  //dicer_args.decrement_recursion_depth();

  return 1;
}

int PairTransformation::_process (Molecule & mol, Dicer_Arguments & arg, Breakages & breakages, Breakages_Iterator & bi) const
{
  Breakages_Iterator bi2(bi);
  if (_t1)
  {
    _t1->_process(mol, arg, breakages, bi2);
  }
  if (_t2)
  {
    _t2->_process(mol, arg, breakages, bi2);
  }
  return 1;
}


int
Chain_Bond_Breakage::_process(Molecule & m0,
                              Dicer_Arguments & dicer_args,
                              Breakages & breakages,
                              Breakages_Iterator & bi) const
{
  if (0 != dicer_args.recursion_depth())         // can only do symmetry suppression at first level
    ;
  else if (this->symmetry_equivalent())
    return 0;

  Molecule m(m0);

  atom_number_t j1, j2;
  if (! convert_to_atom_numbers_in_current_molecule(m, a1(), j1, a2(), j2)) {         // not in this fragment
    return 0;
  }

  if (! m.are_bonded(j1, j2)) {     // should not happen
    return 0;
  }

  set_isotopes_if_needed(m, j1, j2, dicer_args);

#ifdef DEBUG_BOND_BREAKING
  cerr << "Before breaking bond '" << m.smiles() << "', level " << dicer_args.recursion_depth() << endl;
#endif

  assert (m.are_bonded(j1, j2));
  m.remove_bond_between_atoms(j1, j2);

  resizable_array<Molecule *> components;

  if (stream_for_post_breakage.is_open()) {
    do_write_to_stream_for_post_breakage(m, stream_for_post_breakage);

    resizable_array_p<Molecule> c;
    if (2 != m.create_components(c)) {
      cerr << "Huh, created " << c.number_elements() << " fragments during chain bond breakage " << m.smiles() << endl;
      return 0;
    }

    components.add(c[0]);
    components.add(c[1]);
    _process(m0, components, m.name(), dicer_args, breakages, bi);
  }
  else
  {
    Molecule frags;
//  IWString s = m.smiles();

    m.split_off_fragments(0, frags);

//  cerr << "Split " << s << " into " << m.smiles() << " and " << frags.smiles() << endl;

    components.add(&m);
    components.add(&frags);

//  m.unique_smiles();
//  frags.unique_smiles();
//  cerr << "Created " << m.unique_smiles() << " and " << frags.unique_smiles() << endl;

    _process(m0, components, m.name(), dicer_args, breakages, bi);
  }

  return 1;
}

int
Chain_Bond_Breakage::_process(Molecule& m0,
                              resizable_array<Molecule *> & components,
                              const IWString & mname,
                              Dicer_Arguments & dicer_args,
                              Breakages & breakages,
                              Breakages_Iterator & bi) const
{
  assert (2 == components.number_elements());

  for (int i = 0; i < 2; i++)
  {
    Molecule & c = *(components[i]);

    if (c.natoms() < dicer_args.lower_atom_count_cutoff()) {               // cannot be subdivided any more
      continue;
    }

    c.set_name(mname);

#ifdef SPECIAL_THING_FOR_PARENT_AND_FRAGMENTS
    if (ci.natoms() > upper_atom_count_cutoff)
      ;
    else if (!smiles_to_id.contains(c.unique_smiles()))
    {
      Molecule tmp(current_parent_molecule);
      tmp.add_molecule(&c);
      //    cerr << " c contains " << c.number_fragments() << " fragments, tmp " << tmp.number_fragments() << endl;
      stream_for_parent_joinged_to_fragment.write(tmp);
    }
#endif

    if (c.natoms() > upper_atom_count_cutoff)
      ;
    else if (!contains_user_specified_queries(c))
      ;
    else if (!dicer_args.is_unique(c))
      continue;

    dicer_args.WriteAllFragments(m0, c, all_fragments);

#ifdef DEBUG_BOND_BREAKING
    cerr << "Produced '" << c.smiles() << "'\n";
#endif

    call_dicer(c, dicer_args, breakages, bi);
  }

  return 1;
}

static int
remove_fragments_not_containing (Molecule & m,
                                 atom_number_t zatom)
{
  if (1 == m.number_fragments())
    return 0;

  int f = m.fragment_membership(zatom);

  m.delete_all_fragments_except(f);

  return 1;
}

/*
   Breaking a ring is hard.
   First and last items in RBB are the fused ring atoms that do not change
 */

//#define DEBUG_DO_RING_BREAKING

int
Ring_Bond_Breakage::_do_ring_breaking(Molecule & m,
                                      const Set_of_Atoms & rbb,
                                      Dicer_Arguments & dicer_args,
                                      Breakages & breakages,
                                      Breakages_Iterator & bi) const
{
  int * xref = dicer_args.xref_array_this_level();

  initialise_xref(m, dicer_args, xref);

  Set_of_Atoms r;          // atom numbers in the current molecule

  if (! convert_to_atom_numbers_in_child(rbb, xref, r))
    return 0;

  Molecule mcopy(m);

#ifdef DEBUG_DO_RING_BREAKING
  cerr << "Starting _do_ring_breaking with " << r << " mc = " << &mcopy << endl;
  cerr << mcopy.isotopically_labelled_smiles() << endl;
  cerr << "Atoms " << r << endl;
#endif

  int n = r.number_elements();

  if (6 == m.atomic_number(r[n-2]) && 6 == m.atomic_number(r[n-1]))        // don't break CC ring bonds, just to keep numbers down
    return 0;

  if (! mcopy.are_bonded(r[n-2], r[n-1]))       // should not happen
    return 0;

  if (SINGLE_BOND != mcopy.btype_between_atoms(r[n-2], r[n-1]))
    return 0;

  if (! mcopy.is_ring_atom(r[n-2]))
    return 0;

  assert (mcopy.are_bonded(r[n-2], r[n-1]));
  mcopy.remove_bond_between_atoms(r[n-2], r[n-1]);

#ifdef DEBUG_DO_RING_BREAKING
  cerr << "After breaking final bond in ring " << mcopy.isotopically_labelled_smiles() << endl;
#endif

  if (!dicer_args.is_unique(mcopy))
    return 1;

  call_dicer(mcopy, dicer_args, breakages, bi);

#ifdef DEBUG_DO_RING_BREAKING
  cerr << "Begin removal loop " << mcopy.isotopically_labelled_smiles() << endl;
#endif

  for (int i = n-2; i >= 1; i--)
  {
    if (!dicer_args.ok_atom_count(mcopy.natoms() - 1))
      break;

#ifdef  DEBUG_DO_RING_BREAKING
    cerr << "i = " << i << " Processing atom " << r[i] << " in '" << mcopy.isotopically_labelled_smiles() << endl;
#endif

    int need_to_look_for_fragments = (mcopy.ncon(r[i]) > 1);

    if (mcopy.is_aromatic(r[i]))
      continue;

    //  Feb 2012, for now, just pass over these problems, but they should be prevented
    //  {
    //    cerr << "GACK, removing aromatic atom, i = " << i << " atom " << r[i] << endl;
    //    cerr << mcopy.isotopically_labelled_smiles() << endl;
    //      iwabort();
    //   }

    mcopy.remove_atom(r[i]);
    r.adjust_for_loss_of_atom(r[i], 1);               // 1 means keep item in list, do not bother removing it

    if (need_to_look_for_fragments)
      remove_fragments_not_containing(mcopy, r[0]);

    if (!dicer_args.OkAtomCount(mcopy))
      break;

    if (!dicer_args.is_unique(mcopy))
      break;

    call_dicer(mcopy, dicer_args, breakages, bi);

#ifdef DEBUG_DO_RING_BREAKING
    cerr << "Generated " << mcopy.smiles() << endl;
#endif

    if (i > 1)
    {
      initialise_xref(mcopy, dicer_args, xref);
      (void) convert_to_atom_numbers_in_child(rbb, xref, r);                     // expected that not all items present
    }
  }

  return 1;
}

/*
   We need to identify all extra-ring atoms that are attached to atom zatom
 */

static void
identify_extra_ring_attachment (const Molecule & m,
                                const Set_of_Atoms & r,
                                atom_number_t zatom,
                                Set_of_Atoms & exr)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  if (2 == acon)
    return;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->nrings())
      continue;

    exr.add(b->other(zatom));
  }

  return;
}

static int
atoms_in_two_aromatic_rings(Molecule & m,
                            atom_number_t a1,
                            atom_number_t a2)
{
  int nr = m.nrings();

  int aromatic_rings_containing_both = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (!ri->is_aromatic())
      continue;

    if (!ri->contains_bond(a1, a2))
      continue;

    if (1 == aromatic_rings_containing_both)               // we have now found the second one
      return 1;

    aromatic_rings_containing_both++;
  }

  return 0;
}

int
Ring_Bond_Breakage::_process(Molecule & m,
                             Dicer_Arguments & dicer_args,
                             Breakages & breakages,
                             Breakages_Iterator & bi) const
{
  int n = this->number_elements();

  int * xref = dicer_args.xref_array_this_level();

  initialise_xref(m, dicer_args, xref);

  Set_of_Atoms r;           // we make this variable do a bunch of different things...

  if (! convert_to_atom_numbers_in_child(*this, xref, r) || r.empty())
  {
    call_dicer(m, dicer_args, breakages, bi);
    return 0;
  }

  //cerr << "RBB:Processing atoms " << r << " in " << m.isotopically_labelled_smiles() << " other ring " << xref[_atom_that_must_be_a_ring_atom] << endl;

  if (xref[_atom_that_must_be_a_ring_atom] < 0)
    return 0;

  if (m.is_non_ring_atom(xref[_atom_that_must_be_a_ring_atom]))
    return 0;

  if (!m.is_ring_atom(r[0]))
  {
    call_dicer(m, dicer_args, breakages, bi);
    return 0;
  }

  m.compute_aromaticity_if_needed();

  n = r.number_elements();

  int * process_this_bond = new_int(n, 1); std::unique_ptr<int[]> free_process_this_bond(process_this_bond);

  atom_number_t prev = r[0];
  for (int i = 1; i < n; prev = r[i], i++)
  {
    if (!m.are_bonded(prev, r[i]))
    {
      process_this_bond[i] = 0;
      continue;
    }

    if (SINGLE_BOND != m.btype_between_atoms(prev, r[i]))
    {
      process_this_bond[i] = 0;
      continue;
    }

    if (1 == m.nrings(prev) && 1 == m.nrings(r[i]))               // no fusions here
      continue;

    if (!m.is_aromatic(prev) || !m.is_aromatic(r[i]))                 // nothing to worry about
      continue;

    if (atoms_in_two_aromatic_rings(m, prev, r[i]))
    {
      process_this_bond[i] = 0;
      continue;
    }

    if (!_is_arom && m.bond_between_atoms(prev, r[i])->is_aromatic())
      process_this_bond[i] = 0;
  }

  Molecule mcopy(m);

  if (this->_is_arom)           // I don't think we want to break aromatic bonds, not sure this ever happens, besides, the code is wrong, prev does not get updated
  {
    prev = r[0];
    for (int i = 1; i < n; prev = r[i], i++)                // first de-aromatise
    {
      if (!process_this_bond[i])
        continue;

      mcopy.set_bond_type_between_atoms(prev, r[i], SINGLE_BOND);
    }

    if (!dicer_args.is_unique(mcopy))
      return 1;

    call_dicer(mcopy, dicer_args, breakages, bi);
  }

  _do_ring_breaking(mcopy, *this, dicer_args, breakages, bi);

  // Now we want to create the fragment with the end atoms only still attached to the ring
  // Must have at least 5 atoms in the ring for that to work

  if (r.number_elements() < 5)
    return 1;

  mcopy.ring_membership();          // force sssr

  if (!convert_to_atom_numbers_in_child(*this, xref, r))
    return 0;

  Set_of_Atoms exr1, exr2;
  identify_extra_ring_attachment(mcopy, r, r[0], exr1);
  identify_extra_ring_attachment(mcopy, r, r[r.number_elements() - 2], exr2);

  for (int i = 0; i < exr1.number_elements(); i++)
  {
    if (mcopy.are_bonded(r[0], exr1[i]))
      mcopy.remove_bond_between_atoms(r[0], exr1[i]);
  }

  for (int i = 0; i < exr2.number_elements(); i++)
  {
    if (mcopy.are_bonded(r[r.number_elements() - 2], exr2[i]))
      mcopy.remove_bond_between_atoms(r[r.number_elements() - 2], exr2[i]);
  }

  const Atom * retain = mcopy.atomi(r[0]);         // an atom that must be retained

  r.remove_item(0);
  //r.remove_item(0); // removed by madlee so there is NO antenna for a breaked ring.
  r.pop();
  r.pop();

  //cerr << "Residual ring contains " << r.number_elements() << " " << r << endl;

  if (isotope_for_join_points > 0)
  {
    if (increment_isotope_for_join_points)              // not implemented here
    {
      mcopy.set_isotope(r[0], isotope_for_join_points);
      mcopy.set_isotope(r.last_item(), isotope_for_join_points);
    }
    else
    {
      mcopy.set_isotope(r[0], isotope_for_join_points);
      mcopy.set_isotope(r.last_item(), isotope_for_join_points);
    }
  }

  mcopy.remove_atoms(r);

  //cerr << "After removing residual ring " << mcopy.smiles() << "'\n";

  atom_number_t retained = mcopy.which_atom(retain);

  remove_fragments_not_containing(mcopy, retained);

  if (!dicer_args.is_unique(mcopy))
    return 1;

  return call_dicer(mcopy, dicer_args, breakages, bi);
}

int
Fused_Ring_Breakage::_process(Molecule & m,
                             Dicer_Arguments & dicer_args,
                             Breakages & breakages,
                             Breakages_Iterator & bi) const
{
  atom_number_t a11_in_parent = this->a11();
  atom_number_t a12_in_parent = this->a12();
  atom_number_t a1_in_parent  = this->a1();
  atom_number_t a2_in_parent  = this->a2();
  atom_number_t a21_in_parent = this->a21();
  atom_number_t a22_in_parent = this->a22();

  int * xref = dicer_args.xref_array_this_level();

  initialise_xref(m, dicer_args, xref);

  atom_number_t a11, a12, a1, a2, a21, a22;

  a11 = xref[a11_in_parent];
  if (a11 < 0)
    return 0;
  a12 = xref[a12_in_parent];
  if (a12 < 0)
    return 0;
  a1 = xref[a1_in_parent];
  if (a1 < 0)
    return 0;
  a2 = xref[a2_in_parent];
  if (a2 < 0)
    return 0;
  a21 = xref[a21_in_parent];
  if (a21 < 0)
    return 0;
  a22 = xref[a22_in_parent];
  if (a22 < 0)
    return 0;

  if (!m.are_bonded(a11, a1) || !m.are_bonded(a21, a2) || !m.are_bonded(a1, a2))
    return 0;

  Molecule mcopy(m);
  assert (mcopy.are_bonded(a11, a1));
  assert (mcopy.are_bonded(a21, a2));

  mcopy.remove_bond_between_atoms(a11, a1);
  mcopy.remove_bond_between_atoms(a21, a2);

  //#define DEBUG_DO_RING_SPLITTING
#ifdef DEBUG_DO_RING_SPLITTING
  cerr << "From '" << m.smiles() << "' generate '" << mcopy.smiles() << "'\n";
#endif

  if (isotope_for_join_points > 0)
    do_apply_isotopic_labels(mcopy, a1, a2);

  if (! this->aromatic_in_parent())          // no need to check anything
    ;
  else if (mcopy.bond_between_atoms(a1, a2)->is_single_bond() &&
    2 == mcopy.nbonds(a1) && 2 == mcopy.nbonds(a2))
    mcopy.set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);

  if (dicer_args.is_aromatic_in_parent(a12_in_parent, a1_in_parent, a2_in_parent, a22_in_parent))
  {
    if (!mcopy.is_aromatic(a12) || !mcopy.is_aromatic(a1) || !mcopy.is_aromatic(a2) || !mcopy.is_aromatic(a22))
    {
      return 0;
    }
  }

  if (1 == mcopy.number_fragments())
  {
    cerr << "Fused_Ring_Breakage::molecule has just one fragment " << mcopy.isotopically_labelled_smiles() << " " << mcopy.name() << endl;
    cerr << "Removed bonds between atoms " << a11 << " and " << a1 << " AND " << a21 << " and " << a2 << endl;
    cerr << "Starting molecule " << m.isotopically_labelled_smiles() << endl;
  }

  mcopy.remove_fragment_containing_atom(a11);

#ifdef DEBUG_RING_SPLITTING
  cerr << "From '" << m.smiles() << "' generate '" << mcopy.smiles() << "'\n";
#endif

  if (!dicer_args.OkAtomCount(mcopy))
    return 1;

  if (!dicer_args.is_unique(mcopy))
    return 1;

#ifdef DEBUG_RING_SPLITTING
  if (mcopy.aromatic_atom_count() < initial_aromatic_atom_count - 4)
    cerr << "Loses too much aromaticity " << mcopy.aromatic_atom_count() << " vs " << initial_aromatic_atom_count << endl;
#endif

  if (mcopy.natoms() > dicer_args.lower_atom_count_cutoff())
    return call_dicer(mcopy, dicer_args, breakages, bi);
  else
    return 0;
}

///////////////////////////////////////////////////////////////////////
// Spiro_Ring_Breakage::_process
// by madlee @ 2011.12.01
int
Spiro_Ring_Breakage::_process(Molecule & m,
                             Dicer_Arguments & dicer_args,
                             Breakages & breakages,
                             Breakages_Iterator & bi) const
{
  atom_number_t aCenter_in_parent = this->aCenter();

  atom_number_t a11_in_parent = this->a11();
  atom_number_t a12_in_parent = this->a12();

  atom_number_t a21_in_parent = this->a21();
  atom_number_t a22_in_parent = this->a22();

  int * xref = dicer_args.xref_array_this_level();

  initialise_xref(m, dicer_args, xref);

  atom_number_t aCenter, a11, a12, a21, a22;

  aCenter = xref[aCenter_in_parent];
  if (aCenter < 0)
    return 0;

  a11 = xref[a11_in_parent];
  if (a11 < 0)
    return 0;
  a12 = xref[a12_in_parent];
  if (a12 < 0)
    return 0;
  a21 = xref[a21_in_parent];
  if (a21 < 0)
    return 0;
  a22 = xref[a22_in_parent];
  if (a22 < 0)
    return 0;

  if (! m.are_bonded(a11, aCenter) || !m.are_bonded(a12, aCenter) ||  !m.are_bonded(a21, aCenter) || !m.are_bonded(a22, aCenter))
    return 0;

  Molecule mcopy(m);
  mcopy.remove_bond_between_atoms(a11, aCenter);
  mcopy.remove_bond_between_atoms(a12, aCenter);

//#define DEBUG_DO_RING_SPLITTING
#ifdef DEBUG_DO_RING_SPLITTING
  cerr << "From '" << m.smiles() << "' generate '" << mcopy.smiles() << "'\n";
#endif

  if (isotope_for_join_points > 0)
  {
    mcopy.set_isotope(aCenter, isotope_for_join_points);
  }

  if (1 == mcopy.number_fragments())
  {
    cerr << "Spiro_Ring_Breakage::molecule has just one fragment " << mcopy.isotopically_labelled_smiles() << " " << mcopy.name() << endl;
    cerr << "Removed bonds between atoms " << a11 << " and " << aCenter << " AND " << a12 << " and " << aCenter << endl;
    cerr << "Starting molecule " << m.isotopically_labelled_smiles() << endl;
  }

  mcopy.remove_fragment_containing_atom(a11);
  // mcopy.remove_fragment_containing_atom(a12);

#ifdef DEBUG_RING_SPLITTING
  cerr << "From '" << m.smiles() << "' generate '" << mcopy.smiles() << "'\n";
#endif

  if (dicer_args.OkAtomCount(mcopy)) {
    dicer_args.is_unique(mcopy);
  }

  if (mcopy.natoms() > dicer_args.lower_atom_count_cutoff()) {
    return call_dicer(mcopy, dicer_args, breakages, bi);
  } else {
    return 0;
  }
}
//
///////////////////////////////////////////////////////////////////////


int
Dearomatise::_process(Molecule & m,
                      Dicer_Arguments & dicer_args,
                      Breakages & breakages,
                      Breakages_Iterator & bi) const
{
  int * xref = dicer_args.xref_array_this_level();

  initialise_xref(m, dicer_args, xref);

  Set_of_Atoms r;
  r.resize(_ring.number_elements() + 1);

  if (!convert_to_atom_numbers_in_child(_ring, xref, r))
    return 0;

  r.add(r[0]);          // makes processing easier below

  int n = r.number_elements();

  int * process_this_bond = new_int(n, 1);
  std::unique_ptr<int[]> free_process_this_bond(process_this_bond);

  m.compute_aromaticity_if_needed();

  // First make sure we don't destroy an adjacent fused aromatic ring

  atom_number_t prev = r[0];
  for (int i = 1; i < n; prev = r[i], i++)
  {
    if (!m.are_bonded(prev, r[i]))
      return call_dicer(m, dicer_args, breakages, bi);

    if (1 == m.nrings(prev) && 1 == m.nrings(r[i]))               // no fusions here
      continue;

    if (!m.is_aromatic(prev) || !m.is_aromatic(r[i]))                 // nothing to worry about
      continue;

    if (atoms_in_two_aromatic_rings(m, prev, r[i]))
      process_this_bond[i] = 0;
  }

  Molecule mcopy(m);

  prev = r[0];
  for (int i = 1; i < n; prev = r[i], i++)
  {
    if (process_this_bond[i])
      mcopy.set_bond_type_between_atoms(prev, r[i], SINGLE_BOND);
  }

  if (!dicer_args.is_unique(mcopy))
    return 1;

  return call_dicer(mcopy, dicer_args, breakages, bi);
}

static int
dicer(Molecule & m, Dicer_Arguments & dicer_args,
      Breakages & breakages, Breakages_Iterator bi)
{
//#define DEBUG_DICER
#ifdef DEBUG_DICER
  cerr << "Entering dicer " << dicer_args.recursion_depth() <<  ' ' << m.smiles() << " iterator " << bi.state() << endl;
#endif

  const Dicer_Transformation * dt = bi.current_transformation(breakages);
  while (nullptr != dt)
  {
#ifdef DEBUG_DICER
    std::cerr << "####" << dicer_args.recursion_depth() << "  --  " << bi.state() << std::endl;
#endif
    dicer_args.increment_recursion_depth();
    dt->process(m, dicer_args, breakages, bi);
    dicer_args.decrement_recursion_depth();
    dt = bi.current_transformation(breakages);
  }

#ifdef DEBUG_DICER
  cerr << "Returning from dicer " << dicer_args.recursion_depth() << endl;
#endif

  return 1;
}

static int
do_change_element_for_heteroatoms_with_hydrogens (Molecule & m)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = m.atomic_number(i);

    if (6 == zi)
      continue;

    if (0 == m.hcount(i))
      continue;

    if (m.is_ring_atom(i))
      continue;

    if (7 == zi)
      m.set_element(i, neptunium);
    else if (8 == zi)
      m.set_element(i, holmium);
    else if (16 == zi)
      m.set_element(i, strontium);

    rc++;
  }

  return rc;
}

static int
singly_connected_atoms_only (const Molecule & m,
                             atom_number_t zatom,
                             atom_number_t avoid)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == avoid)
      continue;

    if (1 != m.ncon(j))
      return 0;
  }

  return 1;
}

/*
   Detect cases like
 *-zatom-C#N
 *-zatom-C(=O)[NH2]
 */

static int
adjacent_to_terminal_atom (const Molecule & m,
                           atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  if (2 != a->ncon())
    return 0;

  atom_number_t a1 = a->other(zatom, 0);
  atom_number_t a2 = a->other(zatom, 1);

  if (singly_connected_atoms_only(m, a1, zatom))
    return 1;

  if (singly_connected_atoms_only(m, a2, zatom))
    return 1;

  return 0;
}

static int
is_terminal_heteroatom (const Molecule & m,
                        atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  if (1 != a->ncon())
    return 0;

  if (6 == a->atomic_number())
    return 0;

  return 1;
}

/*
   We lower the priority of internal bonds within things like
   amines, sulfonamide, ester, amidine, guanidine

   Basically anything that is [!C]-*=*.
   The caller guarantees that neither atom is in a ring
 */

static int
assign_bond_priority (Molecule & m,
                      atom_number_t a1,
                      atom_number_t a2)
{
  int rc = 25;         // just some random number

  int r1 = 0;          // is a1 in a ring
  if (m.is_ring_atom(a1))
  {
    rc += 10;
    r1 = 1;
  }

  int r2 = 0;          // is a2 in a ring
  if (m.is_ring_atom(a2))
  {
    rc += 10;
    r2 = 1;
  }

  if (r1 && r2)          // bridging bond, highest priority
    return rc;

  if (r1 || r2)         // cannot be amide
    ;
  else if (is_amide(m, a1, a2))         // lower priority of breaking these
    rc -= 10;

  // If it is a terminal heteroatom, on a chain, lower the priority

  if (r1 || r2)
    ;
  else if (is_terminal_heteroatom(m, a1))
    return rc - 5;
  else if (is_terminal_heteroatom(m, a2))
    return rc - 5;

  if (r1)
    ;
  else if (adjacent_to_terminal_atom(m, a1))
    return rc - 4;

  if (r2)
    ;
  else if (adjacent_to_terminal_atom(m, a2))
    return rc - 4;

  // Small lowering for anything not associated with a ring

  if (r1 || r2)
    ;
  else
    rc--;

  return rc;
}

int
Breakages::assign_bond_priorities (Molecule & m,
                                   int * priority) const
{
  int n = _transformations.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    Dicer_Transformation * dti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != dti->ttype())
    {
      priority[i] = 0;
      continue;
    }

    Chain_Bond_Breakage * cbb = dynamic_cast<Chain_Bond_Breakage *>(dti);

    atom_number_t b0i = cbb->a1();
    atom_number_t b1i = cbb->a2();

    priority[i] = assign_bond_priority(m, b0i, b1i);
    if (priority[i] > rc)
      rc = priority[i];
  }

  return rc;
}

/*
   Check to see whether or not this number of bonds to break is consistent with
   the maximum recursion depth and the maximum number of fragments that can be
   produced.
 */

void
Breakages::lower_numbers_if_needed (Molecule & m,
                                    int & bonds_to_break_this_molecule,
                                    int & max_recursion_depth_this_molecule)
{
  const auto f = fard.fragments_at_recursion_depth(max_recursion_depth_this_molecule, bonds_to_break_this_molecule);

  //cerr << "At recursion depth " << max_recursion_depth_this_molecule << " will create " << f << " fragments, cf max " << max_fragments_per_molecule << endl;

  if (f <= max_fragments_per_molecule)        // great
    return;

  // First strategy is to omit undesirable bond breakages.  Only after
  // that will we lower the recursion depth

  int n = _transformations.number_elements();

  int * priority = new int[n]; std::unique_ptr<int[]> free_priority(priority);

  int * to_remove = new_int(n); std::unique_ptr<int[]> free_to_remove(to_remove);

  const int highest_priority = assign_bond_priorities(m, priority);

  int bonds_removed = 0;
  while (1)
  {
    int lowest_priority = std::numeric_limits<int>::max();
    int number_with_min = 0;
    int first_with_min = -1;

    for (int i = 0; i < n; i++)
    {
      if (to_remove[i])
        continue;

      if (0 == priority[i])                   //   probably not a chain bond breakage
        continue;

      if (highest_priority == priority[i])                     // cannot remove those
        continue;

      if (priority[i] > lowest_priority)                       // uninteresting
        continue;
      else if (priority[i] == lowest_priority)
        number_with_min++;
      else
      {
        lowest_priority = priority[i];
        number_with_min = 1;
        first_with_min = i;
      }
    }

    if (first_with_min < 0)                          // nothing found
      break;

    for (int j = first_with_min; j < n; j++)
    {
      if (priority[j] == lowest_priority)
      {
        to_remove[j] = 1;
        bonds_removed++;
      }
    }

    if (fard.fragments_at_recursion_depth(max_recursion_depth_this_molecule, bonds_to_break_this_molecule - bonds_removed) <= max_fragments_per_molecule)
      break;
  }

  for (int i = n - 1; i >= 0; i--)
  {
    if (to_remove[i])
      _transformations.remove_item(i);
  }

  bonds_to_break_this_molecule -= bonds_removed;

  if (fard.fragments_at_recursion_depth(max_recursion_depth_this_molecule, bonds_to_break_this_molecule) <= max_fragments_per_molecule)
    return;

  if (fard.fragments_at_recursion_depth(max_recursion_depth_this_molecule - 1, bonds_to_break_this_molecule) <= max_fragments_per_molecule)
  {
    max_recursion_depth_this_molecule--;
    return;
  }

  // Wow, what to do at this stage!! We have removed low priority breakages, lowering the recursion level will
  // not solve our problem.

  max_recursion_depth_this_molecule--;
  cerr << "Lowered bonds to be broken to " << bonds_to_break_this_molecule << " lowering recursion level to " << max_recursion_depth_this_molecule << " frags " << fard.fragments_at_recursion_depth(max_recursion_depth_this_molecule, bonds_to_break_this_molecule) << endl;

  return;
}

static int
count_side_of_bond (const Molecule & m,
                    atom_number_t astart,
                    int & encountered_ring,
                    int * tmp)
{
  int rc = 1;

  const Atom * a = m.atomi(astart);

  int acon = a->ncon();

  //cerr << "Continue exploring " << astart << endl;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(astart, i);

    if (0 == tmp[j])               // great
      ;
    else if (1 == tmp[j])               // already visited
      continue;
    else if (2 == tmp[j])               // back to the starting atom
      continue;
    else if (2 == tmp[astart])               // starting atom will see the other bond, skip
      continue;
    else if (3 == tmp[j])               // very bad, encountered the starting bond
    {
      cerr << "Atom " << astart << " encountered ring\n";
      encountered_ring = 1;
      return 0;
    }
    else
      continue;

    tmp[j] = 1;
    rc += count_side_of_bond(m, j, encountered_ring, tmp);

    if (encountered_ring)
      return 0;
  }

  //cerr << "From atom " << astart <<" returning " << rc <<endl;
  return rc;
}

static int
atoms_on_side_of_bond_satisfy_min_frag_size(const Molecule & m,
                                  const Dicer_Arguments & dicer_args,
                                  atom_number_t b1,
                                  atom_number_t b2,
                                  int * tmp)
{
  if (1 == m.ncon(b1)) {
    return 1 >= dicer_args.lower_atom_count_cutoff();
  }

  int matoms = m.natoms();

  set_vector(tmp, matoms, 0);

  tmp[b1] = 2;
  tmp[b2] = 3;          // special flag for this atom

  int encountered_ring = 0;

  int c = count_side_of_bond(m, b1, encountered_ring, tmp);

  //cerr << "From " << b2 << " to " << b1 << " count " << c << ", ring " << encountered_ring << endl;

  // If there is a ring, no violation of lower_atom_count_cutoff possible
  if (encountered_ring) {
    return 1;
  }

  if (c < dicer_args.lower_atom_count_cutoff()) {  // does not satisfy constraint
    return 0;
  }

  //cerr << "Constraint is OK!\n";
  return 1;
}

int
Breakages::discard_breakages_that_result_in_fragments_too_small(Molecule & m,
                                          const Dicer_Arguments & dicer_args)
{
  int matoms = m.natoms();

  int * tmp = new int[matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  int rc = 0;

  for (int i = _transformations.number_elements() - 1; i >= 0; i--)
  {
    const Dicer_Transformation * dti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != dti->ttype()) {
      continue;
    }

    const Chain_Bond_Breakage * b = dynamic_cast<const Chain_Bond_Breakage *>(dti);

    atom_number_t b0 = b->a1();
    atom_number_t b1 = b->a2();

    //  If either is singly bonded, we know the number of atoms removed by
    //  breaking the bond

    if (lower_atom_count_cutoff > 1 && (1 == m.ncon(b0) || 1 == m.ncon(b1)))
    {
      _transformations.remove_item(i);
      rc++;
      continue;
    }

    if (atoms_on_side_of_bond_satisfy_min_frag_size(m, dicer_args, b0, b1, tmp) &&
      atoms_on_side_of_bond_satisfy_min_frag_size(m, dicer_args, b1, b0, tmp))
      ;
    else
    {
      _transformations.remove_item(i);
      rc++;
      //    cerr << "Removing breakage " << b0 << " to " << b1 << endl;
    }
  }

  return rc;
}

/*
   The hard coded rules are
   aromatic - *
   heteroatom - *
   unsaturated - *
   heteroatom -C - *
 */

int
Breakages::identify_bonds_to_break_hard_coded_rules(Molecule & m)
{
  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (const Bond* b : m.bond_list()) {
    if (b->nrings()) {
      continue;
    }

    if (!b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    const int rbc1 = m.ring_bond_count(a1);
    const int rbc2 = m.ring_bond_count(a2);

    if (include_amide_like_bonds_in_hard_coded_queries)
      ;
    else if (rbc1 || rbc2)
      ;
    else if (is_amide(m, a1, a2))
      continue;

    if (m.is_aromatic(a1) || m.is_aromatic(a2))
    {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    if (!m.saturated(a1) || !m.saturated(a2)) {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    if (! break_ring_chain_bonds) {
    } else if ((rbc1 == 0 && rbc2 > 0) || (rbc1 > 0 && rbc2 == 0)) {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    atomic_number_t z1 = m.atomic_number(a1);
    atomic_number_t z2 = m.atomic_number(a2);

    // heteroatoms at both ends of bond
    if (6 != z1 && 6 != z2) {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    //  At least one heteroatom. Don't split up CF3 groups

    if (m.is_halogen(a1) || m.is_halogen(a2))
      continue;

    //  Do allow the bond to a CF3 or t-Butyl to break

    if (is_cf3_like(m, *b))
    {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    //  Breaking any bond to a heteroatom can generate a lot of molecules. Add
    //  Some heuristics about needed number of atoms in the fragment

    if (6 != z1 || 6 != z2)
    {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }

    if (allow_cc_bonds_to_break && 6 == z1 && 6 == z2)
    {
      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(a1, a2);
      _transformations.add(b);
      rc++;
      continue;
    }
  }

//#define ECHO_BONDS_BROKEN
#ifdef ECHO_BONDS_BROKEN
  cerr << m.isotopically_labelled_smiles() << ' ' << m.name() << endl;
  debug_print(cerr);
#endif

  return rc;
}

int
Breakages::identify_bonds_to_break (Molecule & m)
{
  if (0 == nq) {
    return identify_bonds_to_break_hard_coded_rules(m);
  }

  int rc = 0;

  if (add_user_specified_queries_to_default_rules) {
    rc += identify_bonds_to_break_hard_coded_rules(m);
  }

  Molecule_to_Match target(&m);

  int queries_hitting = 0;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    sresults.set_remove_invalid_atom_numbers_from_new_embeddings(1);

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    queries_hitting++;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      if (e->number_elements() < 2)
      {
        cerr << "dicer:gack, query " << i << " only has " << e->number_elements() << " matched atoms\n";
        return 0;
      }

      if (m.in_same_ring(e->item(0), e->item(1)))
      {
        cerr << "Cannot process ring bonds, atoms " << e->item(0) << " and " << e->item(1) << ", query " << i << " against '" << m.name() << "'\n";
        return 0;
      }

      Chain_Bond_Breakage * b = new Chain_Bond_Breakage(e->item(0), e->item(1));
      _transformations.add(b);

      //    if (m.is_ring_atom(e->item(0)) && m.is_ring_atom(e->item(1)))
      //      _ring_bonds.add(b);
      //    else
      //      _chain_bonds.add(b);

      rc++;
    }
  }

  if (verbose > 1)
    cerr << queries_hitting << " of " << nq << " queries produced " << rc << " matches to '" << m.name() << "'\n";

  if (0 == rc)
  {
    if (ignore_molecules_not_matching_any_queries)
    {
      molecules_not_hitting_any_queries++;
      return 1;
    }

    if (verbose)
      cerr << "None of " << nq << " queries matched '" << m.name() << "'\n";

    return 0;
  }

  return rc;
}



void
Breakages::_remove_fused_ring_breakage_if_present (atom_number_t a1,
  atom_number_t a2)
{
  for (int i = _transformations.number_elements() - 1; i >= 0; i--)
  {
    const Dicer_Transformation * dti = _transformations[i];

    if (TRANSFORMATION_TYPE_FUSED_RING_BREAKAGE != dti->ttype())
      continue;

    const Fused_Ring_Breakage * frb = dynamic_cast<const Fused_Ring_Breakage *>(dti);

    if (!frb->involves(a1, a2))                 // check central atom only
      continue;

    _transformations.remove_item(i);
    return;
  }

  return;
}

void
Breakages::_remove_spiro_ring_breakage_if_present (atom_number_t a1,
  atom_number_t a2)
{
  for (int i = _transformations.number_elements() - 1; i >= 0; i--)
  {
    const Dicer_Transformation * dti = _transformations[i];

    if (TRANSFORMATION_TYPE_SPIRO_RING_BREAKAGE != dti->ttype())
      continue;

    const Spiro_Ring_Breakage * frb = dynamic_cast<const Spiro_Ring_Breakage *>(dti);

    if (!frb->involves(a1, a2))                 // check central atom only
      continue;

    _transformations.remove_item(i);
    return;
  }

  return;
}

void
Breakages::_remove_breakage_if_present (atom_number_t a1,
                                        atom_number_t a2)
{
  for (int i = _transformations.number_elements() - 1; i >= 0; i--)
  {
    const Dicer_Transformation * dti = _transformations[i];

    if (TRANSFORMATION_TYPE_BREAK_CHAIN_BOND != dti->ttype())
      continue;

    const Chain_Bond_Breakage * cbi = dynamic_cast<const Chain_Bond_Breakage *>(dti);

    if (!cbi->involves(a1, a2))
      continue;

    _transformations.remove_item(i);
    return;
  }

  return;
}

int
Breakages::identify_bonds_to_not_break(Molecule & m)
{
  Molecule_to_Match target(&m);

  int queries_hitting = 0;

  for (int i = 0; i < queries_for_bonds_to_not_break.number_elements(); i++)
  {
    Substructure_Results sresults;

    int nhits = queries_for_bonds_to_not_break[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

//  cerr << nhits << " hits to query " << i << " in '" << m.name() << "'\n";

    queries_hitting++;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      if (e->number_elements() < 2)
      {
        cerr << "dicer:gack, query " << i << " only has " << e->number_elements() << " matched atoms\n";
        continue;
      }

//    cerr << "Processing query match " << *e << endl;

      _remove_breakage_if_present(e->item(0), e->item(1));
      _remove_fused_ring_breakage_if_present(e->item(0), e->item(1));
      _remove_spiro_ring_breakage_if_present(e->item(0), e->item(1));
    }
  }

  return _transformations.number_elements();
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  m.revert_all_directional_bonds_to_non_directional();         // by default

  if (discard_chirality) {
    m.remove_all_chiral_centres();
  }

  if (element_transformations.active())
    element_transformations.process(m);

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (perceive_terminal_groups_only_in_original_molecule) {
    int matoms = m.natoms();

    for (int i = 0; i < matoms; i++)
    {
      if (1 != m.ncon(i))
        continue;

      if (m.nbonds(i) == m.ncon(i))
        m.set_isotope(i, static_cast<int>(Env::kTerminal));
    }
  }

  return;
}


static int
do_number_by_initial_atom_number(Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 1; i < matoms; i++)
  {
    m.set_isotope(i, i);
  }

  return 1;
}

static int
write_fully_broken_parent(const Molecule & m,
                          const Breakages &breakages,
                          Molecule_Output_Object & stream_for_fully_broken_parents)
{
  Molecule mcopy(m);

  for (int i = 0; i < breakages.number_transformations(); ++i)
  {
    const auto t = breakages.transformation(i);

    t->break_bonds(mcopy);
  }

  return stream_for_fully_broken_parents.write(mcopy);
}

static void
add_parent_atom_count_to_global_array(const IWString & mname,
                                      int matoms,
                                      IW_STL_Hash_Map_int & atoms_in_parent)
{
  int ispace = mname.index(' ');

  if (ispace < 0)
  {
    atoms_in_parent[mname] = matoms;
    return;
  }

  IWString tmp(mname);
  tmp.truncate_at_first(' ');

  atoms_in_parent[tmp] = matoms;

  return;
}

static int
dicer(Molecule & m,
      DicerFragmentOutput & output)
{
  static bool first_call = true;
  if (first_call) {
    if (m.highest_coordinate_dimensionality() > 1) {
      set_append_coordinates_after_each_atom(1);
    }
  }

  preprocess(m);
  // std::cerr << "####" << m.name() << std::endl;

  if (number_by_initial_atom_number)
    do_number_by_initial_atom_number(m);

  if (accumulate_starting_parent_information)
    add_parent_atom_count_to_global_array(m.name(), m.natoms(), atoms_in_parent);

  Breakages breakages;

  int bonds_to_break_this_molecule;

  if (read_bonds_to_be_broken_from_bbrk_file)
    bonds_to_break_this_molecule = breakages.get_bonds_to_break_from_text_info(m);
  else
    bonds_to_break_this_molecule = breakages.identify_bonds_to_break(m);

  if (queries_for_bonds_to_not_break.number_elements())
    bonds_to_break_this_molecule = breakages.identify_bonds_to_not_break(m);

  if (break_fused_rings || break_ring_bonds)
  {
    // top_level_kekule_forms_switched += arrange_kekule_forms_in_fused_rings(m);
    // *****************     disabled by madlee @2012.03.13     **********************
    // re arrange the kekule is dangers and we try to change the bond order in fragment
    // if possible.
    //  cerr << "After toggling bonds for preferred Kekule form '" << m.smiles() << "'\n";
  }

  if (break_ring_bonds)
    bonds_to_break_this_molecule += breakages.identify_ring_bonds_to_be_broken(m);

  if (break_fused_rings)
    bonds_to_break_this_molecule += breakages.identify_cross_ring_bonds_to_be_broken(m);

  if (break_spiro_rings)
    bonds_to_break_this_molecule += breakages.identify_spiro_ring_bonds_to_be_broken(m);

  if (dearomatise_aromatic_rings)
    bonds_to_break_this_molecule += breakages.identify_aromatic_rings_to_de_aromatise(m);

  Dicer_Arguments dicer_args(max_recursion_depth);

  if (max_atoms_lost_from_parent > 0)
    dicer_args.set_lower_atom_count_cutoff(m.natoms() - max_atoms_lost_from_parent);
  else if (min_atom_fraction > 0.0f)
    dicer_args.set_lower_atom_count_cutoff(static_cast<int>(static_cast<float>(m.natoms()) * min_atom_fraction));
  else
    dicer_args.set_lower_atom_count_cutoff(lower_atom_count_cutoff);

  dicer_args.set_upper_atom_count_cutoff(std::min(upper_atom_count_cutoff,
                                                  int(max_atom_fraction*m.natoms())));

  if (0 == lower_atom_count_cutoff)         // no need to check
    ;
  else if (reject_breakages_that_result_in_fragments_too_small)
    bonds_to_break_this_molecule -= breakages.discard_breakages_that_result_in_fragments_too_small(m, dicer_args);

  int max_recursion_depth_this_molecule = max_recursion_depth;

  if (max_fragments_per_molecule > 0)
    breakages.lower_numbers_if_needed(m, bonds_to_break_this_molecule, max_recursion_depth_this_molecule);

  if (bbrk_file.active())
    return breakages.write_bonds_to_be_broken(m, bbrk_file);

  global_counter_bonds_to_be_broken[bonds_to_break_this_molecule]++;

  if (verbose > 1)
    cerr << "In '" << m.name() << "', identified " << bonds_to_break_this_molecule << " bonds to break\n";

  /*if (bonds_to_break_this_molecule > max_bonds_allowed_in_input)
     {
     return 1;
     }*/

  if (fingerprint_tag.length()) {
  } else if (! write_diced_molecules) {
  } else if (write_parent_smiles) {
    output << m.smiles() << ' ' << m.name() << " B=" << bonds_to_break_this_molecule << '\n';
  }

  if (0 == bonds_to_break_this_molecule)
  {
    if (fingerprint_tag.length())
      dicer_args.produce_fingerprint(m, output.text_output());

    return 1;
  }

  if (stream_for_fully_broken_parents.active())
    write_fully_broken_parent(m, breakages, stream_for_fully_broken_parents);

  dicer_args.set_current_molecule(&m);

  if (change_element_for_heteroatoms_with_hydrogens)
    do_change_element_for_heteroatoms_with_hydrogens(m);

  if (atom_typing_specification.active()) {
    (void) dicer_args.store_atom_types(m, atom_typing_specification);
  }

  if (max_fragments_per_molecule > 0)
    dicer_args.adjust_max_recursion_depth(bonds_to_break_this_molecule, max_fragments_per_molecule);

  dicer_args.set_write_diced_molecules(write_diced_molecules);

  breakages.identify_symmetry_exclusions(m);

  dicer_args.set_array_size(m.natoms());

  if (break_fused_rings)
    dicer_args.initialise_parent_molecule_aromaticity(m);

  USPVPTR * uspvptr = new USPVPTR[m.natoms()]; std::unique_ptr<USPVPTR[]> free_uspvptr(uspvptr);

  dicer_lib::initialise_atom_pointers(m, dicer_args.atom_types(), uspvptr);

  // int * atom_numbers = new int[m.natoms()];
  // std::unique_ptr<int[]> free_atom_numbers(atom_numbers);
  // initialise_atom_pointers(m, atom_numbers);

  if (work_like_recap) {
    breakages.do_recap(m, dicer_args);

    // The number of breakable bonds is not used with recap output, set to -1.
    if (write_smiles) {
      dicer_args.write_fragments_found_this_molecule(m, -1, output);
    }
    if (fingerprint_tag.length()) {
      dicer_args.produce_fingerprint(m, output.text_output());
    }

    return 1;
  }

  if (m.number_fragments() > 1) {
    cerr << m.smiles() << ' ' << m.name() << '\n';
  }

  assert (1 == m.number_fragments());

  Breakages_Iterator bi(breakages);

  if (! dicer(m, dicer_args, breakages, bi))
    return 0;

  // verbose is checked here because we only produce the count summary if verbose

  //if (write_smiles || || verbose)
  //  dicer_args.sort_fragments_found_this_molecule();

  if (write_smiles) {
    dicer_args.write_fragments_found_this_molecule(m, bonds_to_break_this_molecule, output);
  }

  if (fingerprint_tag.length())
    dicer_args.produce_fingerprint(m, output.text_output());

  if (write_smiles_and_complementary_smiles) {
    dicer_args.write_fragments_and_complements(m, output);
  }

  if (verbose)
    fragments_produced[dicer_args.fragments_found()]++;

  return 1;
}

static int
dicer(data_source_and_type<Molecule> & input,
      DicerFragmentOutput & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    clock_t t0;
    if (collect_time_statistics)
      t0 = clock();
    else
      t0 = static_cast<clock_t>(0);                     // keep the compiler quiet

#ifdef SPECIAL_THING_FOR_PARENT_AND_FRAGMENTS
    current_parent_molecule = *m;
#endif

    if (!dicer(*m, output)) {
      cerr << "dicer:error processing " << m->name() << '\n';
      return 0;
    }

    if (collect_time_statistics) {
      clock_t t = (clock() - t0) / 1000;

      // have observed this going negative at times. Presumably something wrapping
      if (t >= 0) {
        cerr << m->name() << " took " << t << " ms\n";
        time_acc.extra(t);
      }
    }

    output.MaybeFlush();
  }

  return 1;
}

static int
dicer_tdt_filter_(const_IWSubstring buffer,    // local copy
                  DicerFragmentOutput & output)
{
  assert (buffer.ends_with('>'));

  buffer.chop();

  buffer.remove_leading_chars(smiles_tag.length());

  Molecule m;

  if (!m.build_from_smiles(buffer))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return dicer(m, output);
}

static int
dicer_tdt_filter(iwstring_data_source & input,
                 DicerFragmentOutput & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag))
      continue;

    if (!dicer_tdt_filter_(buffer, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
dicer_tdt_filter(const char * fname,
                DicerFragmentOutput & output)
{
  iwstring_data_source input(fname);

  if (!input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return dicer_tdt_filter(input, output);
}

static int
InitialiseEnvironmentElements() {
  env_element[0] = nullptr;

  env_element[static_cast<int>(Env::kTerminal)] = get_element_from_symbol_no_case_conversion("Te");
  env_element[static_cast<int>(Env::kAromatic)] = get_element_from_symbol_no_case_conversion("Ar");
  env_element[static_cast<int>(Env::kSaturatedCarbon)] = get_element_from_symbol_no_case_conversion("Cs");
  env_element[static_cast<int>(Env::kSaturatedNitrogen)] = get_element_from_symbol_no_case_conversion("Nh");
  env_element[static_cast<int>(Env::kSaturatedOxygen)] = get_element_from_symbol_no_case_conversion("Os");

  env_element[static_cast<int>(Env::kUnsaturatedCarbon)] = get_element_from_symbol_no_case_conversion("Cu");
  env_element[static_cast<int>(Env::kUnsaturatedNitrogen)] = get_element_from_symbol_no_case_conversion("Ne");
  env_element[static_cast<int>(Env::kUnsaturatedOxygen)] = get_element_from_symbol_no_case_conversion("Og");

  env_element[static_cast<int>(Env::kAmideCarbon)] = get_element_from_symbol_no_case_conversion("Ac");
  env_element[static_cast<int>(Env::kAmideNitrogen)] = get_element_from_symbol_no_case_conversion("Am");

  env_element[static_cast<int>(Env::kHalogen)] = get_element_from_symbol_no_case_conversion("Hg");
  env_element[static_cast<int>(Env::kNitro)] = get_element_from_symbol_no_case_conversion("No");

  env_element[static_cast<int>(Env::kAliphaticRingCarbon)] = get_element_from_symbol_no_case_conversion("Rb");
  env_element[static_cast<int>(Env::kAliphaticRingNitrogen)] = get_element_from_symbol_no_case_conversion("Rn");

  return 1;
}

static int
display_misc_B_options (std::ostream & os)
{
  os << " -B atype=<tag>    Calculate atom type and dump it in frag/comp pairs\n";
  os << " -B term           perceive terminal groups\n";
  os << " -B bscb           allow Carbon-Carbon single bonds to break\n";
  os << " -B eH             change O,N and S with Hydrogens to other elements\n";
  os << " -B xsub           do not report fragments that are exact subsets of others\n";
  os << " -B nosmi          suppress output of fragment smiles\n";
  os << " -B noparent       suppress output of parent smiles\n";
  os << " -B proto          main output is a dicer_data.DicedMolecule proto\n";
  os << " -B serialized_proto main output is tfdatarecord of serialized dicer_data::DicedMolecule\n";
  os << " -B time           run timing\n";
  os << " -B addq           run the -q queries in addition to the default rules\n";
  os << " -B WB=fname       write smiles of just broken molecules to <fname>\n";
  os << " -B BRB            break some kinds of cross ring bond pairs\n";
  os << " -B BRF            break fused rings\n";
  os << " -B BRS            break spiro rings\n";
  os << " -B BBRK_FILE=<fname> write bonds to break to <fname> - no computation\n";
  os << " -B BBRK_TAG=<tag> read bonds to be broken from sdf tag <tag> (def DICER_BBRK)\n";
  os << " -B nousmi         do NOT uuse unique smiles to keep track of fragments - dangerous\n";
  os << " -B nbamide        do NOT break amide and amide-like bonds\n";
  os << " -B brcb           break all bonds from a chain to a ring\n";
  // the nbfts option is not working, TODO ianwatson fix.
  //os << " -B nbfts          do NOT break bonds that would yield fragments too small\n";
  os << " -B appnatoms      append the atom count to each fragment produced\n";
  os << " -B MAXAL=<n>      max atoms lost when creating any fragment\n";
  os << " -B MINFF=<f>      discard fragments that comprise less than <f> fraction of atoms in parent\n";
  os << " -B MAXFF=<f>      discard fragments that comprise more than <f> fraction of atoms in parent\n";
  os << " -B fragstat=<fname> write fragment statistics to <fname>\n";
  os << " -B fstsp=<float>  fractional support level for inclusion in fragstat= file\n";
  os << " -B fragstatproto  write the fragstat file as a DicerFragments text proto\n";
  os << " -B fragstatbinproto write the fragstat file as TFDatarecord DicerFragments serialized proto's\n";
  os << " -B smiles_proto   write the fragstat file as 'smiles DicerFragments text proto'\n";
  os << " -B WRITEALL=<fname> write all fragments found to <fname>. Can be HUGE!\n";
  os << " -B FMC:<qry>      queries that fragments reported must contain\n";
  os << " -B FMC:heteroatom fragments reported must contain a heteroatom, :n for count\n";
  os << " -B VF=<n>         output will be viewed in vf with <n> structures per page\n";
  os << " -B flush          flush output after each molecule\n";
  os << " -B lostchiral     check each fragment for lost chirality\n";
  os << " -B recap          work like Recap - break all breakable bonds at once\n";
  os << " -B nooutput       suppress writing diced fragments, useful with the '-B fragstat' option\n";
  os << " -B help           this message\n";

  exit(3);
}

static int
display_dash_K_options(char flag, std::ostream & os)
{
  os << " -" << flag << " READ=<fname>    read existing data from <fname>\n";
  os << " -" << flag << " WRITE=<fname>   when complete, write data to <fname>\n";
  os << " -" << flag << " ADD             add new strings to the hash (default)\n";
  os << " -" << flag << " NOADD           do NOT add new strings to the hash\n";
  os << "                       useful for creating fixed dictionary fingerprints\n";
  os << " -" << flag << " help            this message\n";

  exit(3);
}

static int
dicer(const char * fname, FileType input_type,
      DicerFragmentOutput& output)
{
  assert (nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return dicer(input, output);
}

static void
display_dash_i_options(std::ostream & os)
{
  os << " -I env           atoms get isotopic labels according to their environment\n";
  os << " -I enva          specific atoms are added to indicate the environment\n";
  os << " When the 'env' or 'enva' directives are used, these isotopes or elements are added\n";
  os << " Descripton               iso element\n";
  os << " Terminal                 1   Te\n";
  os << " Aromatic                 2   Ar\n";
  os << " SaturatedCarbon          3   Cs\n";
  os << " SaturatedNitrogen        4   Nh\n";
  os << " SaturatedOxygen          5   Os\n";
  os << " UnsaturatedCarbon        6   Cu\n";
  os << " UnsaturatedNitrogen      7   Cu\n";
  os << " UnsaturatedOxygen        8   Cu\n";
  os << " AmideCarbon              9   Ac\n";
  os << " AmideNitrogen           10   Am\n";
  os << " Halogen                 11   Hg\n";
  os << " Nitro                   12   No\n";
  os << " RingCarbonAliphatic     13   Cr\n";
  os << " RingNitrogenAliphatic   14   Na\n";
  os << " -I z             atoms labelled by atomic number of neighbour\n";
  os << " -I ini           atoms labelled by initial atom number (debugging uses)\n";
  os << " -I <number>      constant isotopic label applied to all join points\n";
  os << " -I inc=<n>       increment existing isotopic label before labelling join points\n";

  exit(0);
}

static int
DisplayComplementaryFragmentOptions(std::ostream& output) {
  output << " -C def            complementary fragments with all default conditions\n";
  output << " -C auto           complementary fragments have same isotope as fragment\n";
  output << " -C atype          complementary fragments have isotopes based on atom type\n";
  output << " -C iso=<iso>      apply isotope <n> to the join point on the complementary fragment\n";
  output << " -C maxat=<n>      do not write complementary fragments with more than n> atoms\n";

  ::exit(0);
}

static int
ParseComplementaryFragmentOptions(Command_Line& cl) {
  if (! cl.option_present('C')) {
    return 1;
  }

  const_IWSubstring c;
  for (int i = 0; cl.value('C', c, i); ++i) {
    if (c.starts_with("maxat=") || c.starts_with("maxat:")) {
      c.remove_leading_chars(6);
      if (! c.numeric_value(upper_atom_count_cutoff_complementary_fragment) ||
          upper_atom_count_cutoff_complementary_fragment < 1) {
        cerr << "The upper atom count limit for complementary fragments must be a +ve whole number\n";
        return 0;
      }
      if (verbose) {
        cerr << "Will not write complementary fragmes with more than " << 
               upper_atom_count_cutoff_complementary_fragment << " atoms\n";
      }
    } else if (c == "def") {
    } else if (c == "auto") {
      apply_isotopes_to_complementary_fragments = -1;
      if (verbose)
        cerr << "Same isotope will be set to the fragments and complementary fragments\n";
    } else if (c == "atype") {
      apply_isotopes_to_complementary_fragments = ISO_CF_ATYPE;
      if (verbose)
        cerr << "Isotopes at complementary atoms will be initial molecule atom type\n";
    } else if (c.starts_with("iso=") || c.starts_with("iso:")) {
      c.remove_leading_chars(4);

      if (!c.numeric_value(apply_isotopes_to_complementary_fragments) ||
           apply_isotopes_to_complementary_fragments < 0) {
        cerr << "Invalid complementary set isotope specification '" << c << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Isotope " << apply_isotopes_to_complementary_fragments << " applied to join points of complementary fragments\n";
      }
    } else if (c == "help") {
      DisplayComplementaryFragmentOptions(cerr);
    } else {
      cerr << "Unrecognised -C qualifier '" << c << "'\n";
      DisplayComplementaryFragmentOptions(cerr);
    }
  }

  write_smiles_and_complementary_smiles = 1;
  if (verbose) {
    cerr << "Will write smiles and complementary smiles\n";
  }

  return 1;
}

static int
dicer (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:ls:q:z:ek:X:I:J:cm:M:T:B:K:n:N:d:C:afG:S:Z:P:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E'))
  {
    if (!process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('G'))
  {
    if (! process_standard_smiles_options(cl, verbose, 'G'))
    {
      cerr << "Cannot initialise smiles specific options (-G)\n";
      return 2;
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  set_reset_implicit_hydrogens_known_on_bond_removal(1);

  if (cl.option_present('T'))
  {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T'))
    {
      cerr << "Cannot initialise element transformations (-T)\n";
      usage(3);
    }
  }

  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);
  set_copy_atom_based_user_specified_void_pointers_during_add_molecule(1);

  if (cl.option_present('k'))
  {
    if (!cl.value('k', max_recursion_depth) || max_recursion_depth < 1)
    {
      cerr << "The maximum recursion depth (-k) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will break a max of " << max_recursion_depth << " bonds simultaneously\n";

    initialise_fragments_per_number_of_bonds_array();
  }

  fard.initialise(max_recursion_depth, 64);

  if (verbose > 1)
    fard.debug_print(cerr);

  if (cl.option_present('X'))
  {
    if (!cl.value('X', max_fragments_per_molecule) || max_fragments_per_molecule < 1)
    {
      cerr << "The maximum fragments per molecule (-X) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will produce a max of " << max_fragments_per_molecule << " fragments per molecule\n";
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring i;

    for (int j = 0; cl.value('I', i, j); ++j)
    {
      if ("env" == i)
      {
        apply_environment_isotopic_labels = 1;
        if (verbose)
          cerr << "Broken bonds marked with isotopes to indicate environment\n";
      }
      else if ("ini" == i)
      {
        number_by_initial_atom_number = 1;
        if (verbose)
          cerr << "Will apply isotopic labels according to initial atom numbers\n";
      }
      else if ('z' == i)
      {
        apply_atomic_number_isotopic_labels = 1;
        if (verbose)
          cerr << "Will label cut points by atomic number of removed atom\n";
      }
      else if ("enva" == i)
      {
        add_environment_atoms = 1;
      }
      else if ("depth" == i)
      {
        isotopic_label_is_recursion_depth = 1;
        if (verbose)
          cerr << "Isotopic labels will be according to bonds cut\n";
      }
      else if (i.starts_with("inc="))
      {
        i.remove_leading_chars(4);
        if (! i.numeric_value(increment_isotope_for_join_points) || increment_isotope_for_join_points < 1)
        {
          cerr << "The increment isotope values for attachment points directive (-I inc=) must be a whole +ve number\n";
          return 2;
        }

        if (verbose)
          cerr << "Will add " << increment_isotope_for_join_points << " to existing isotopic labels at join points\n";
      }
      else if ("atype" == i)
      {
        apply_atom_type_isotopic_labels = 1;
        if (verbose)
          cerr << "Isotopic labels will be atom types\n";
      }
      else if ("help" == i)
      {
        display_dash_i_options(cerr);
      }
      else
      {
        if (! i.numeric_value(isotope_for_join_points) || isotope_for_join_points < 0) {
          cerr << "The isotope for join points value (-I) must be a whole +ve number '" << i << "' invalid\n";
          usage(4);
        }

        if (verbose)
          cerr << "Attachment points indicated with isotope " << isotope_for_join_points << endl;
      }
    }

    if (add_environment_atoms && apply_environment_isotopic_labels) {
      cerr << "Adding environment atoms and applying isotopic labels doesn't make sense\n";
      usage(4);
    }

    if (add_environment_atoms) {
      InitialiseEnvironmentElements();
    }
  }

  if (cl.option_present('r')) {
    work_like_recap = 1;

    if (verbose)
      cerr << "Will work like recap - single level recursion\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', fingerprint_tag);

    if (! fingerprint_tag.starts_with("NC")) {
      cerr << "Dicer produces sparse fingerprints only, tag must start with 'NC', '" << fingerprint_tag << " invalid\n";
      return 1;
    }

    if (verbose)
      cerr << "Will create fingerprints with tag '" << fingerprint_tag << "'\n";

    if (!fingerprint_tag.ends_with('<'))
      fingerprint_tag.add('<');

    write_smiles = 0;

    if (cl.option_present('f'))
    {
      function_as_filter = 1;
      if (verbose)
        cerr << "Will work as a TDT filter\n";
    }
  }

  if (cl.option_present('m'))
  {
    if (!cl.value('m', lower_atom_count_cutoff) || lower_atom_count_cutoff < 0)
    {
      cerr << "The minimum number of atoms in a fragment (-m) option must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard fragments with fewer than " << lower_atom_count_cutoff << " atoms\n";

    lower_atom_count_cutoff_current_molecule = lower_atom_count_cutoff;
  }

  if (cl.option_present('M'))
  {
    const_IWSubstring m;
    for (int i = 0; cl.value('M', m, i); ++i) {
      if (m.starts_with("maxnr=")) {
        m.remove_leading_chars(6);
        if (! m.numeric_value(max_non_ring_atoms) || max_non_ring_atoms < 0) {
          cerr << "Invalid maxnr= qualifier '" << m << "'\n";
          return 1;
        }
      } else if (! m.numeric_value(upper_atom_count_cutoff) ||
                upper_atom_count_cutoff < 1 || upper_atom_count_cutoff < lower_atom_count_cutoff) {
        cerr << "The maximum number of atoms in a fragment (-M) option must be a whole +ve number\n";
        return 1;
      } else if (verbose) {
        if (verbose) {
          cerr << "Will discard fragments with more than " << upper_atom_count_cutoff << " atoms\n";
        }
      }

      upper_atom_count_cutoff_current_molecule = upper_atom_count_cutoff;
    }
  }

  if (cl.option_present('c'))
  {
    discard_chirality = 1;

    if (verbose)
      cerr << "Chirality will be discarded\n";
  }

  IWString name_for_fragment_statistics_file;

  Fragment_Statistics_Report_Support_Specification fsrss;

  DicerFragmentOutput output;

  // We can write the fragments generated as binary protos.
  int write_fragments_as_binary_protos = 0;

  // Too hard to check support levels since some may be numeric, some may be percentages....

  if (cl.option_present('B'))
  {
    IWString bond_break_file;

    int i = 0;
    const_IWSubstring b;
    while (cl.value('B', b, i++))
    {
      if ("term" == b)
      {
        use_terminal_atom_type = 1;
        if (verbose)
          cerr << "Will use the terminal atom type\n";
      }
      else if ("bscb" == b)
      {
        allow_cc_bonds_to_break = 1;
        if (verbose)
          cerr << "Will allow Carbon-Carbon single bonds to break\n";
      }
      else if ("eH" == b)
      {
        change_element_for_heteroatoms_with_hydrogens = 1;
        if (verbose)
          cerr << "Heteroatoms with Hydrogens will be transformed\n";

        holmium = get_element_from_symbol_no_case_conversion("Ho");
        neptunium = get_element_from_symbol_no_case_conversion("Np");
        strontium = get_element_from_symbol_no_case_conversion("Sr");
      }
      else if ("xsub" == b)
      {
        eliminate_fragment_subset_relations = 1;
        if (verbose)
          cerr << "Will suppress output of fragments that are exact subsets of other fragments\n";

        initialise_atomic_number_to_array_position();
      }
      else if ("nosmi" == b)
      {
        write_smiles = 0;
        write_parent_smiles = 0;

        if (verbose)
          cerr << "No smiles written\n";
      }
      else if ("noparent" == b)
      {
        write_parent_smiles = 0;

        if (verbose)
          cerr << "Will not write the parent smiles\n";
      }
      else if ("time" == b)
      {
        collect_time_statistics = 1;
        if (verbose)
          cerr << "Will collect time statistics\n";
      }
      else if (b.starts_with("fragstat="))
      {
        b.remove_leading_chars(9);
        name_for_fragment_statistics_file = b;
        if (verbose)
          cerr << "Fragment statistics will be written to '" << name_for_fragment_statistics_file << "'\n";

        accumulate_starting_parent_information = 1;
      }
      else if (b.starts_with("fstsp="))
      {
        b.remove_leading_chars(6);

        if (!fsrss.build(b))
        {
          cerr << "Cannot initialise fragment report support specification '" << b << "'\n";
          return 2;
        }

        if (verbose)
          fsrss.debug_print(cerr);
      }
      else if (b == "fragstatproto") {
        write_fragstat_as_proto = 1;
        if (verbose) {
          cerr << "Fragstat file written as proto\n";
        }
      }
      else if (b == "smiles_proto") {
        write_fragstat_as_proto = 1;
        prepend_smiles_to_fragstat_proto = 1;
        if (verbose) {
          cerr << "Will prepend the smiles to the fragstat proto\n";
        }
      }
      else if (b == "fragstatbinproto") {
        write_fragstat_as_binary_proto = 1;
        if (verbose) {
          cerr << "Will write fragstat file as a TFDatarecord proto\n";
        }
      }
      else if (b.starts_with("WRITEALL")) {
        if (! all_fragments.ParseDirective(b)) {
          cerr << "Unrecognosed '" << b << "'\n";
          return 0;
        }
      }
      else if (b == "addq")
      {
        add_user_specified_queries_to_default_rules = 1;
        if (verbose)
          cerr << "User specified bond breaking queries will be run in addition to the default set\n";
      }
      else if ("lostchiral" == b)
      {
        check_for_lost_chirality = 1;
        if (verbose)
          cerr << "Will check for lost chirality on output\n";
      }
      else if (b.starts_with("WB="))
      {
        IWString fname = b;
        fname.remove_leading_chars(3);
        if (!fname.ends_with(".smi"))
          fname << ".smi";

        if (!stream_for_post_breakage.open(fname.null_terminated_chars()))
        {
          cerr << "Cannot open stream for post breakage '" << fname << "'\n";
          return 5;
        }

        if (verbose)
          cerr << "Immediately broken molecules written to '" << fname << "'\n";
      }
      else if ("BRF" == b)
      {
        break_fused_rings = 1;
        set_display_abnormal_valence_messages(0);
        set_display_messages_about_unable_to_compute_implicit_hydgogens(0);

        if (verbose)
          cerr << "Will dice across some rings\n";
      }
      else if ("BRS" == b)
      {
        break_spiro_rings = 1;
        if (verbose)
          cerr << "Will dice across spiro rings\n";
      }
      else if ("BRB" == b)
      {
        break_ring_bonds = 1;
        if (verbose)
          cerr << "Will allow ring bonds to break\n";
      }
      else if ("DEAROM" == b || b.starts_with("DEAR") || b.starts_with("dear"))
      {
        dearomatise_aromatic_rings = 1;
        if (verbose)
          cerr << "Will de-aromatise aromatic rings\n";
      }
      else if (b.starts_with("BBRK_FILE="))
      {
        b.remove_leading_chars(10);
        bond_break_file = b;
      }
      else if (b.starts_with("BBRK_TAG="))
      {
        b.remove_leading_chars(9);
        bbrk_tag = b;
        read_bonds_to_be_broken_from_bbrk_file = 1;
      }
      else if ("BBRK_TAG" == b)
      {
        bbrk_tag = "DICER_BBRK";
        read_bonds_to_be_broken_from_bbrk_file = 1;
      }
      else if ("nbamide" == b)
      {
        include_amide_like_bonds_in_hard_coded_queries = 0;

        if (verbose)
          cerr << "Amide and amide-like bonds will not be broken\n";
      }
      else if ("nbfts" == b)
      {
        reject_breakages_that_result_in_fragments_too_small = 1;
        cerr << "Warning the nbfts option is not working\n";

        if (verbose)
          cerr << "Breakages that would yield fragments too small not allowed\n";
      }
      else if ("appnatoms" == b)
      {
        append_atom_count = 1;
        if (verbose)
          cerr << "Will append the atom count to the output\n";
      }
      else if (b.starts_with("atype="))
      {
        b.remove_leading_chars(6);
        if (! atom_typing_specification.build(b))
        {
          cerr << "Invalid atom typing specification '" << b << "'\n";
          return 2;
        }
      }
      else if (b.starts_with("MAXAL="))
      {
        b.remove_leading_chars(6);
        if (!b.numeric_value(max_atoms_lost_from_parent) || max_atoms_lost_from_parent <= 0)
        {
          cerr << "The max atoms lost directive (MAXAL=) must be a whole +ve number\n";
          return 3;
        }

        if (verbose)
          cerr << "Fragments must represent a loss of " << max_atoms_lost_from_parent << " or fewer atoms from parent\n";
      }
      else if (b.starts_with("MINFF="))
      {
        b.remove_leading_chars(6);
        if (!b.numeric_value(min_atom_fraction) || min_atom_fraction < 0.0f || min_atom_fraction > 1.0f)
        {
          cerr << "The min fraction atoms (MINFF=) directive must specify a fraction\n";
          return 2;
        }

        if (verbose)
          cerr << "Will only produce fragments that cover at least " << min_atom_fraction << " of the atoms in the parent\n";
      }
      else if (b.starts_with("MAXFF="))
      {
        b.remove_leading_chars(6);
        if (!b.numeric_value(max_atom_fraction) || max_atom_fraction < 0.0f || max_atom_fraction > 1.0f)
        {
          cerr << "The max fraction atoms (MAXFF=) directive must specify a fraction\n";
          return 2;
        }

        if (verbose)
          cerr << "Will only produce fragments that cover at least " << max_atom_fraction << " of the atoms in the parent\n";
      }
      else if (b.starts_with("FMC:heteroatom"))
      {
        b.remove_leading_chars(14);
        if (0 == b.length())
          fragments_must_contain_heteroatom = 1;
        else
        {
          if (':' != b[0] || 1 == b.length())
          {
            cerr << "To specify the number of heteroatoms, must be of the form 'FMC:heteroatoms:n' '" << b << "' invalid\n";
            usage(1);
          }

          b.remove_leading_chars(1);

          if (! b.numeric_value(fragments_must_contain_heteroatom) || fragments_must_contain_heteroatom < 0)
          {
            cerr << "The FMC:heteroatom: directive must be followed by a whole +ve number\n";
            usage(1);
          }
        }
        if (verbose)
          cerr << "Fragments reported must contain at least " << fragments_must_contain_heteroatom << " heteroatoms\n";
      }
      else if (b.starts_with("FMC:"))
      {
        b.remove_leading_chars(4);
        if (! process_cmdline_token(' ', b, user_specified_queries_that_fragments_must_contain, verbose))
        {
          cerr << "Cannot process Fragment Must Contain specification '" << b << "'\n";
          return 3;
        }
      }
      else if (b == "brcb") {
        break_ring_chain_bonds = 1;
        if (verbose) {
          cerr << "Will break ring/chain bonds\n";
        }
      }
      else if (b.starts_with("VF="))
      {
        b.remove_leading_chars(3);

        if (!b.numeric_value(vf_structures_per_page) || vf_structures_per_page < 1)
        {
          cerr << "The insert parent every n written value must be a whole +ve integer, '" << b << "' is invalid\n";
          return 3;
        }

        if (verbose)
          cerr << "Will insert parent structure into output every " << vf_structures_per_page << " fragments written\n";
      }
      else if ("flush" == b)
      {
        output.set_buffered_output(0);
        if (verbose)
          cerr << "Output flushed after each molecule\n";
      }
      else if ("recap" == b)
      {
        work_like_recap = 1;
        if (verbose)
          cerr << "Will work like Recap\n";
      }
      else if (b == "nooutput") {
        write_diced_molecules = 0;
        if (verbose) {
          cerr << "Normal output suppressed\n";
        }
      }
      else if (b == "proto") {
        write_parent_smiles = 0;
        if (verbose) {
          cerr << "Output is dicer_data::DicedMolecule textproto\n";
        }
        output.set_output_is_proto(1);
      }
      else if (b == "serialized_proto") {
        write_fragments_as_binary_protos = 1;
        write_parent_smiles = 0;
        if (verbose) {
          cerr << "Will write fragments as serialised dicer_data::DicedMolecule protos\n";
        }
        output.set_output_is_proto(1);
      }
      else if (b == "nousmi") {
        fragments_indexed_by_unique_smiles = 0;
        if (verbose) {
          cerr << "Fragments will NOT be indexed by unique smiles - do not use fragstat!\n";
        }
      }
      else if("help" == b)
      {
        display_misc_B_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        return 4;
      }

      if (user_specified_queries_that_fragments_must_contain.number_elements())
        cerr << "Definied " << user_specified_queries_that_fragments_must_contain.number_elements() << " queries for fragments\n";
    }

    if (bond_break_file.length() && bbrk_tag.length())
    {
      cerr << "Cannot specify both BBRK_FILE and BBRK_TAG\n";
      return 3;
    }

    if (bond_break_file.length())
    {
      if (!bond_break_file.ends_with(".sdf"))
        bond_break_file << ".sdf";

      bbrk_file.add_output_type(FILE_TYPE_SDF);
      if (!bbrk_file.new_stem(bond_break_file))
      {
        cerr << "Cannot open stream for bond breaks '" << bond_break_file << "'\n";
        return 8;
      }

      bbrk_tag = "DICER_BBRK";

      if (verbose)
        cerr << "Bonds to be broken written in .sdf form to '" << bond_break_file << "', tag '" << bbrk_tag << "'\n";
    }
    else
    {
      moleculeio::set_read_extra_text_info(1);
      MDL_File_Supporting_Material * mdlfs = global_default_MDL_File_Supporting_Material();
      mdlfs->set_report_unrecognised_records(0);
    }

    all_fragments.UseDefaultFeaturesUnlessSpecified();
  }

  if (write_fragments_as_binary_protos && ! cl.option_present('S')) {
    cerr << "When writing serialized protos, must specify output file via the -S option\n";
    return 1;
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p)) {
      cerr << "Cannot initialise atom typing '" << p << "'\n";
      return 1;
    }
  }

  if (cl.option_present('a'))
  {
    write_only_smallest_fragment = 1;
  }

  if (cl.option_present('C')) {
    // will exit if it fails.
    ParseComplementaryFragmentOptions(cl);
  }

#ifdef NOW_IN_FUNCTION
  if (cl.option_present('C')) {
    const_IWSubstring c;
    for (int i = 0; cl.value('C', c, i); ++i) {
      if (c.starts_with("maxat=") || c.starts_with("maxat:")) {
        c.remove_leading_chars(6);
        if (! c.numeric_value(upper_atom_count_cutoff_complementary_fragment) ||
            upper_atom_count_cutoff_complementary_fragment < 1) {
          cerr << "The upper atom count limit for complementary fragments must be a +ve whole number\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will not write complementary fragmes with more than " << 
                 upper_atom_count_cutoff_complementary_fragment << " atoms\n";
        }
      } else if (c == "def") {
      } else if (c == "auto") {
        apply_isotopes_to_complementary_fragments = -1;
        if (verbose)
          cerr << "Same isotope will be set to the fragments and complementary fragments\n";
      } else if (c == "atype") {
        apply_isotopes_to_complementary_fragments = ISO_CF_ATYPE;
        if (verbose)
          cerr << "Isotopes at complementary atoms will be initial molecule atom type\n";
      } else if (c.starts_with("iso=") || c.starts_with("iso:")) {
        c.remove_leading_chars(4);

        if (!c.numeric_value(apply_isotopes_to_complementary_fragments) ||
             apply_isotopes_to_complementary_fragments < 0) {
          cerr << "Invalid complementary set isotope specification '" << c << "'\n";
          return 0;
        }
        if (verbose) {
          cerr << "Isotope " << apply_isotopes_to_complementary_fragments << " applied to join points of complementary fragments\n";
        }
      } else if (c == "help") {
        DisplayComplementaryFragmentOptions(cerr);
      } else {
        cerr << "Unrecognised -C qualifier '" << c << "'\n";
        DisplayComplementaryFragmentOptions(cerr);
      }
    }

    write_smiles_and_complementary_smiles = 1;
    if (verbose) {
      cerr << "Will write smiles and complementary smiles\n";
    }
  }
#endif  // NOW_IN_FUNCTION

  if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(smarts))
      {
        delete q;
        return 60 + i;
      }

      if (verbose)
        cerr << "Created query from smarts '" << smarts << "'\n";

      queries.add(q);
    }
  }

  if (cl.option_present('q'))
  {
    if (!process_queries(cl, queries, verbose, 'q'))
    {
      cerr << "Cannot read queries, -q option\n";
      usage(3);
    }
  }

  nq = queries.number_elements();

  if (nq)
  {
    if (verbose)
      cerr << "Read " << nq << " queries\n";

    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_find_unique_embeddings_only(1);
    }

  }
  else
  {
    if (add_user_specified_queries_to_default_rules)
    {
      cerr << "Have been asked to run queries in addition to default rules, but no extra queries (-q)\n";
      return 3;
    }

    if (verbose)
      cerr << "No bond breaking queries specified, will use default rules\n";
  }

  if (cl.option_present('n'))
  {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('n', smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(smarts))
      {
        delete q;
        return 60 + i;
      }

      if (verbose)
        cerr << "Created query for non-breaking motif from smarts '" << smarts << "'\n";

      queries_for_bonds_to_not_break.add(q);
    }
  }

  if (cl.option_present('N'))
  {
    if (!process_queries(cl, queries_for_bonds_to_not_break, verbose, 'N'))
    {
      cerr << "Cannot read queries, -q option\n";
      usage(3);
    }
  }

  for (int i = 0; i < queries_for_bonds_to_not_break.number_elements(); i++)
  {
    queries_for_bonds_to_not_break[i]->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++))
    {
      if ('i' == z)
      {
        ignore_molecules_not_matching_any_queries = 1;
        if (verbose)
          cerr << "Will ignore molecules not hitting any of the queries\n";
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(8);
      }
    }
  }

  if (write_smiles && fingerprint_tag.length())
  {
    cerr << "Sorry, doesn't make sense to write both a smiles and fingerprints\n";
    return 3;
  }

  if (write_smiles) {
  } else if (fingerprint_tag.length()) {
  } else if (write_smiles_and_complementary_smiles) {
  } else if (! name_for_fragment_statistics_file.empty()) {
  } else {
    cerr << "NO output specified\n";
    usage(4);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (!process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (function_as_filter)
    ;
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (all_files_recognised_by_suffix(cl))
    ;
  else
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString write_fname;

  if (cl.option_present('K'))
  {
    IWString read_fname;

    int i = 0;
    const_IWSubstring k;
    while(cl.value('K', k, i++))
    {
      if (k.starts_with("READ="))
      {
        read_fname = k;
        read_fname.remove_leading_chars(5);
      }
      else if (k.starts_with("WRITE="))
      {
        write_fname = k;
        write_fname.remove_leading_chars(6);
      }
      else if ("ADD" == k)
      {
        add_new_strings_to_hash = 1;
        if (verbose)
          cerr << "New features added to hash\n";
      }
      else if ("NOADD" == k)
      {
        add_new_strings_to_hash = 0;
        if (verbose)
          cerr << "New features not added to hash\n";
      }
      else if ("help" == k)
      {
        display_dash_K_options('K', cerr);
      }
      else
      {
        cerr << "Unrecognised -K qualifier '" << k << "'\n";
        usage(4);
      }
    }

    if (read_fname.length())
    {
      if (!read_existing_hash_data(read_fname, smiles_to_id))
      {
        cerr << "Cannot read existing feature->bit number data from '" << read_fname << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Read " << smiles_to_id.size() << " smiles to bit number relationships from '" << read_fname << "'\n";
    }

    if (write_fname.length())
    {
      if (!bit_smiles_cross_reference.open(write_fname.null_terminated_chars()))
      {
        cerr << "Cannot open bit->fragment smiles cross reference file '" << write_fname << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Bit ->fragment smiles cross reference written to '" << write_fname << "'\n";
    }
  }

#ifdef SPECIAL_THING_FOR_PARENT_AND_FRAGMENTS
  stream_for_parent_joinged_to_fragment.add_output_type(FILE_TYPE_SDF);
  if (!stream_for_parent_joinged_to_fragment.new_stem("dcpjoined"))
  {
    cerr << "Cannot open 'dcpjoined' file\n";
    return 4;
  }
  cerr << "Warning, running special version for parent/fragment joins\n";
#endif

  if (cl.option_present('Z'))
  {
    IWString z;
    cl.value('Z', z);

    if (z.ends_with(".smi"))
      stream_for_fully_broken_parents.add_output_type(FILE_TYPE_SMI);
    else if (z.ends_with(".sdf"))
      stream_for_fully_broken_parents.add_output_type(FILE_TYPE_SDF);
    else
    {
      z << ".sdf";
      stream_for_fully_broken_parents.add_output_type(FILE_TYPE_SDF);
    }

    if (stream_for_fully_broken_parents.would_overwrite_input_files(cl, z))
    {
      cerr << "Stream for fully broken parent molecules '" << z << "' would overwrite input file(s)\n";
      return 1;
    }

    if (! stream_for_fully_broken_parents.new_stem(z))
    {
      cerr << "Cannot open stream for fully broken parent molecules '" << z << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Will write fully broken parent molecules to '" << z << "'\n";
  }

  IW_Time tstart;

  if (write_fragments_as_binary_protos) {
    IWString fname = cl.option_value('S');
    if (! output.SetTFDataFname(fname)) {
      cerr << "Cannot open stream for serialized protos '" << fname << "'\n";
      return 1;
    }
  } else if (cl.option_present('S')) {
    IWString fname = cl.option_value('S');

    if (! output.SetOutputStream(fname)) {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 1;
    }
  } else {
    output.SetUseStdout();
  }

  for (const char * fname: cl) {
    int rc;
    if (function_as_filter)
      rc = dicer_tdt_filter(fname, output);
    else
      rc = dicer(fname, input_type, output);

    if (0 == rc) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (eliminate_fragment_subset_relations)
    do_eliminate_fragment_subset_relations();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    for (int i = 0; i < global_counter_bonds_to_be_broken.number_elements(); i++)
    {
      if (global_counter_bonds_to_be_broken[i])
        cerr << global_counter_bonds_to_be_broken[i] << " molecules had " << i << " bonds to be broken\n";
    }

    if (smiles_to_id.size()) {
      cerr << smiles_to_id.size() << ' ';
      if (fragments_indexed_by_unique_smiles) {
        cerr << "unique";
      } 
      cerr << " fragments encountered\n";
    }

    if (break_fused_rings)
      cerr << top_level_kekule_forms_switched << " top level Kekule forms switched\n";
  }

  if (collect_time_statistics)
  {
    cerr << "Processing took between " << time_acc.minval() << " and " << time_acc.maxval() << ", ave " << time_acc.average() << endl;
  }

  cerr.flush();

  if (bit_smiles_cross_reference.is_open())
    write_hash_data(bit_smiles_cross_reference, smiles_to_id);

  IW_Time tend;

  if (verbose)
    cerr << "Computation took " << (tend - tstart) << " seconds\n";

  if (0 == molecules_containing.size())
    ;
  else if (name_for_fragment_statistics_file.empty()) {
  } else if (write_fragstat_as_binary_proto) {
    WriteFragmentStatisticsAsBinaryProto(molecules_containing, fsrss, name_for_fragment_statistics_file);
  }
  else if (name_for_fragment_statistics_file.length()) {
    write_fragment_statistics(molecules_containing, fsrss, name_for_fragment_statistics_file);
  }

  if (stream_for_post_breakage.is_open())
    stream_for_post_breakage.close();

  if (verbose)
  {
    for (int i = 0; i < fragments_produced.number_elements(); i++)
    {
      if (fragments_produced[i])
        cerr << fragments_produced[i] << " molecules produced " << i << " fragments\n";
    }
  }

  return 0;
}

#ifndef _WIN32_

int
main(int argc, char ** argv) {
  prog_name = argv[0];

  int rc = dicer(argc, argv);

  return rc;
}

#else

extern "C"
{
  int dicer_csharp(int argc, char **argv)
  {
    return dicer(argc,argv);
  }

}

int
Extern_Dicer::extern_dicer(int argc, char ** argv)
{
  reset_variables();

  verbose=1;
  int rc = dicer(argc, argv);
        
  return rc;
}
#endif
