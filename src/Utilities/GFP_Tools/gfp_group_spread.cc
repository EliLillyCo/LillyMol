// Spread variant where molecules must be selected in groups

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"


#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp_standard.h"

namespace gfp_group_spread {

using std::cerr;

static Fraction_as_String fraction_as_string;

// clang-format off
static void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Spread variant that selects molecules in groups.
For each molecule, a group membership must be specified. That can come from one of two sources:
  Via the -G option specifying a descriptor file containing a mapping of ID->group:
Header Header
ID1 group1
ID2 group1
ID2 group2

Or the group name can be a part of the name in the fingerprint
PCN<ID1 group1>
with '-G col=2' which means that the group membership is column 2 of the name.
Every fingerprint must have a group assignment.

The -D option specifies a mapping from id to desirability. This is a descriptor file.
Each group forms the average desirability of all fingerprints in it. Not all ids need to 
have a desirability assigned, but those without will have a default zero value assigned.

During the spread selection each group keeps track of distances to previously selected groups.
With the -P option, you can specify a file of distances associated with each identifier. If the
-P file name ends with '.gfp' it is treated as a fingerprint file and computations are performed.

A typical usage might be

gfp_group_spread -G file.grp file.gfp

where 'file.grp' contains group memberships.

The following options are recognised.

-G <fname>      Descriptor file containing group assignements.
-G col=<col>    The group membership of each molecule is the <col> column of the name.
-P <fname>      Descriptor file containing previous distances.
-P <fname.gfp>  GFP file containing previously selected items and distances will be computed.
-D <fname>      Descriptor file containing relative desirability values.
-T <dist>       Truncate long distances at <dist>. Long gfp distances may not be meaningful.
-n <nsel>       Number of groups to select.
-w ...          How to weight the various factors in the score. Enter '-w help' for info.
-v              verbose output. 
)";

  ::exit(rc);
}

class Fingerprint : public GFP_Standard {
  private:
    IWString _name;

    IWString _group;

  public:
    Fingerprint();

    void set_group(const IWString& s) {
      _group = s;
    }
    const IWString& group() const {
      return _group;
    }

    void set_name(const IWString& s) {
      _name = s;
    }
    const IWString& name() const {
      return _name;
    }

    // Name is currently 'name .. grp'.
    // Set _group and truncate the name to the first token.
    int GroupNameIsColumn(int col);
};

Fingerprint::Fingerprint() {
}

int
Fingerprint::GroupNameIsColumn(int col) {
  if (_name.nwords() <= col) {
    cerr << "Fingerprint::GroupNameIsColumn:not enough columns '" << _name << "'\n";
    return 0;
  }
  _group = _name.word(col);
  _name.truncate_at_first(' ');

  return 1;
}

class Group {
  private:
    // Will be set once this group has been selected
    // It will be the order in which our group was found.
    uint32_t _selected;

    IWString _group_name;

    resizable_array<Fingerprint*> _fp;

    // Each fingerprint can have an associated desirability based on the
    // desirability of each member of the group.
    Accumulator<double> _desirability;

    // The distances to previously selected for all our fingerprints.
    Accumulator<double> _distance;

    // We can optionally compute a measure of internal diversity.
    double _intra_group_ave_distance;

    // Maybe it would be interesting to keep track of all distances.
    // Convert to a number in the range [0,100)
    // We allocate 101 values just in case a distance of 1.0 ever happens
    // which is hard to imagine.
    uint32_t* _distance_histogram;

    // We don't want to give too much weight to long distances
    // that are likely to be non meaningful.
    static inline float _truncate_long_distances = std::numeric_limits<float>::max();
    // If _truncate_long_distances is set, we can do a one-time computation of the
    // index into _distance_histogram.
    static inline uint32_t _truncate_index = 0;

    // While the closest distance is of interest, as might be the average,
    // we might also be interested in the prevalence of short distances.
    static inline float _short_distance_cutoff = 0.0f;
    uint64_t _short_distances;

    static inline char _output_separator = ' ';

    // We may have a penalty on exact matches.
    uint32_t _zero_distances;

    // private functions
  public:
    Group(const IWString& id);
    ~Group();

    const IWString& name() const {
      return _group_name;
    }

    uint32_t size() const {
      return _fp.size();
    }

    void set_selected(uint32_t s) {
      _selected = s;
    }
    uint32_t selected() const {
      return _selected;
    }

    void SetShortDistanceCutoff(float s);
    void SetTruncateLongDistances(float s);

    void AnotherDistance(float d);
    void AnotherDesirability(float d) {
      _desirability.extra(d);
    }

    int EstablishInitialDistances(Group& rhs);

    int AnotherMember(Fingerprint* fp);

    void GroupHasBeenSelected(const Group* sel);

    // Will compute this if it is not already computed.
    double IntraGroupAverageDistance();

    // Note the empty case returns 1 (very distant) rather than 0 (very close).
    // The previously measured distance option does not work without this.
    double AverageDistance() const {
      if (_distance.empty()) {
        return 1.0;
      } else {
        return _distance.average();
      }
    }

    // If we are keeping track of short distances, we can return the
    // fraction of distances that are outside that threshold.
    // The reason for such a strange measure is that we want to return something
    // where larger values are better.
    double FractionNotShort() const;

    const Accumulator<double> DistanceAccumulator() const {
      return _distance;
    }

    double AverageDesirability() const {
      if (_desirability.empty()) {
        return 0.0;
      } else {
        return _desirability.average();
      }
    }

    int Write(IWString_and_File_Descriptor& output) const;
};

Group::Group(const IWString& id) : _group_name(id) {
  _selected = 0;
  _zero_distances = 0;
  _distance_histogram = new uint32_t[101];
  std::fill_n(_distance_histogram, 101, 0);
  _short_distances = 0;
  _intra_group_ave_distance = -1.0;  // Initialised negative so we know it has not been computed.

  _fp.reserve(96);   // 96 well plates is a common use case.

  // Ran into problems with the -P option, previously measured distances.
  // If some fingerprints do not have a value, they will appear as being zero distance,
  // which is not really correct. Therefore initialise all groups with an infinite distance.
  AnotherDistance(1.0f);
}

Group::~Group() {
  if (_distance_histogram != nullptr) {
    delete [] _distance_histogram;
  }
}

void
Group::SetTruncateLongDistances(float s) {
  _truncate_long_distances = s;
  _truncate_index = static_cast<uint32_t>(s * 100.0f);
}

void
Group::SetShortDistanceCutoff(float s) {
  _short_distance_cutoff = s;
}

int
Group::AnotherMember(Fingerprint* fp) {
  _fp << fp;

  return 1;
}

// Note that this does not respect _truncate_long_distances.
// Not sure if that is a bug or a feature.
// It is also unclear what should be returned if this is a singleton group?
double
Group::IntraGroupAverageDistance() {
  if (_intra_group_ave_distance > 0.0) {
    return _intra_group_ave_distance;
  }

  _intra_group_ave_distance = 0.0;
  uint32_t n = _fp.size();
  for (uint32_t i = 0; i < n; ++i) {
    for (uint32_t j = i + 1; j < n; ++j) {
      float d = _fp[i]->tanimoto_distance(*_fp[j]);
      _intra_group_ave_distance += d;
    }
  }

  _intra_group_ave_distance /= static_cast<double>(n);

  return _intra_group_ave_distance;
}

double
Group::FractionNotShort() const {
  if (_short_distances == 0) {
    return 0.0;
  }

  if (_distance.empty()) {  // not sure what we should do here.
    return 0.0;
  }

  return iwmisc::Fraction<double>(_short_distances, _distance.n());
}

int
Group::EstablishInitialDistances(Group& rhs) {
  const uint32_t n1 = _fp.size();
  const uint32_t n2 = rhs._fp.size();
  for (uint32_t i = 0; i < n1; ++i) {
    const Fingerprint* fpi = _fp[i];
    for (uint32_t j = 0; j < n2; ++j) {
      float d = fpi->tanimoto_distance(*rhs._fp[j]);
      AnotherDistance(d);
      rhs.AnotherDistance(d);
    }
  }

  return 1;
}


void
Group::AnotherDistance(float d) {
  if (d > _truncate_long_distances) {
    _distance.extra(_truncate_long_distances);
    if (_distance_histogram != nullptr) {
      ++_distance_histogram[_truncate_index];
    }
  } else {
    _distance.extra(d);
    if (d == 0.0f) {
      ++_zero_distances;
    }

    if (_distance_histogram != nullptr) {
      uint32_t ndx = static_cast<uint32_t>(d * 100.0f);
      // cerr << "Updating histogram, ndx " << ndx << "\n'";
      _distance_histogram[ndx] += 1;
      if (d < _short_distance_cutoff) {
        ++_short_distances;
      }
    }
  }
}

void
Group::GroupHasBeenSelected(const Group* sel) {
  for (const Fingerprint* rhs : sel->_fp) {
    for (Fingerprint* lhs : _fp) {
      float d = rhs->tanimoto_distance(*lhs);
      AnotherDistance(d);
    }
  }
}

int
Group::Write(IWString_and_File_Descriptor& output) const {
  output << _group_name;
  fraction_as_string.append_number(output, AverageDistance());
  fraction_as_string.append_number(output, AverageDesirability());

  for (const Fingerprint* f : _fp) {
    output << _output_separator << f->name();
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

class Options {
  private:
    // The number of groups to select - the -n option.
    uint32_t _nsel;

    uint32_t _number_fingerprints;
    Fingerprint* _fp;

    // A mapping from names to index in the _fp array.
    absl::flat_hash_map<IWString, uint32_t> _id_to_ndx;

    IWString _identifier_tag;

    // Will be set if -G COL= is specified.
    int _group_is_column;

    // Each unique group name will populate one of these.
    resizable_array_p<Group> _group;

    // A mapping from group names to an index in the _group array.
    absl::flat_hash_map<IWString, uint32_t> _group_name_to_ndx;

    // The overall desirability of a group of molecules is a composite
    // of how desirable the molecules are and the distance to previously
    // selected items.
    // Note that the desirability should be in the range [0,1] in order
    // to ensure comparability.
    double _ave_distance_weight;
    // The larger the minimum distance, the more desirable the group.
    double _min_distance_weight;
    // For each group compute the average intra-member distance. If this is 
    // large it is desirable.
    double _intra_group_weight;
    // Weight given to the average desirability of members of the group.
    double _desirability_weight;

    // If we have used the -t option, there will be a concept of "short distances".
    // That is the FractionNotShort() value returned from a group.
    double _non_short_distance_weight;

    int _verbose;

    // private functions
    uint32_t ReadFingerprints(iwstring_data_source& input);
    int ReadDesirabilty(IWString& fname);
    int ReadDesirabilty(iwstring_data_source& input);
    int ReadDesirabiltyRecord(const const_IWSubstring& buffer);
    int ReadGroupMembership(iwstring_data_source& input);
    int ReadPreviousDistance(iwstring_data_source& input);
    int ReadGroupMembershipRecord(const const_IWSubstring& buffer);
    int ReadPreviousDistanceRecord(const const_IWSubstring& buffer);
    int GroupsFromNameToken();
    int ComputePreviousDistances(IWString& fname);
    int ComputePreviousDistances(iwstring_data_source& input);

    int DetermineInitialDistances();

    Group* GetNextGroup() const;
    double Score(Group& g) const;

  public:
    Options();
    ~Options();

    int Initialise(const Command_Line& cl);

    uint32_t ReadFingerprints(const char* fname);
    int ReadGroupMembership(IWString& fname);
    int ReadPreviousDistance(IWString& fname);

    int Spread(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _number_fingerprints = 0;
  _fp = nullptr;
  _verbose = 0;
  _nsel = 0;

  _group_is_column = 0;

  _ave_distance_weight = 0.5;
  _min_distance_weight = 0.0;
  _desirability_weight = 0.5;
  _non_short_distance_weight = 0.0;
  _intra_group_weight = 0.0;

  _identifier_tag = "PCN<";
}

Options::~Options() {
  if (_fp != nullptr) {
    delete [] _fp;
  }
}

void
DisplayDashwOptions(int rc) {
  cerr << R"(The -w option controls the relative weight given to average desirability and distance to previously selected.
These two properties need to be combined in order to come up with a unified score for each group of molecules.
The -w option specifies the relative weights assigned to each measure of fitness.
All weights are relative weights and should sum to 1.0, although this is not enforced.

 -w zero                set all weights to zero - should be the first -w option - otherwise will clear previous -w values.
 -w avdist=<w>          the relative weight assigned to average distance from previously selected.
 -w mindist=<w>         the relative weight assigned to the shortest distance from previously selected.
 -w nonshort=<w>        the relative weight assigned to the fraction of non-short distances (-t).
 -w desr=<w>            the relative weight assigned to the score of the molecules in the group.
 -w intra=<w>           the relative weight of the average intra-plate distance.
)";

  ::exit(rc);
}

int
Options::Initialise(const Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! ReadFingerprints(cl[0])) {
    cerr << "Options::Initialise:cannot read fingerprints '" << cl[0] << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Read " << _number_fingerprints << " fingerprints\n";
  }

  // The -w nonshort= directive is only meaningful if this is set.
  bool short_distance_cutoff_specified = false;

  if (cl.option_present('t')) {
    float tshort;
    if (! cl.value('t', tshort) || tshort < 0.0f || tshort > 1.0f) {
      cerr << "Options::Initialise:invalid short distance cutoff (-t)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will count distances less than " << tshort << '\n';
    }

    Group notused("notused");
    notused.SetShortDistanceCutoff(tshort);
    short_distance_cutoff_specified = true;
  }

  if (cl.option_present('T')) {
    float tld;
    if (! cl.value('T', tld) || tld < 0.0f || tld > 1.0f) {
      cerr << "Options::Initialise:invalid distance cutoff (-T)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will truncate distances at " << tld << '\n';
    }

    Group notused("notused");
    notused.SetTruncateLongDistances(tld);
  }

  // This is not needed.
  // DetermineInitialDistances();

  if (! cl.option_present('G')) {
    cerr << "Options::Initialise:Must specify Group memberships with the -G option\n";
    Usage(1);
  }

  if (cl.option_present('G')) {
    IWString fname = cl.string_value('G');
    if (fname.starts_with("COL=")) {
      fname.remove_leading_chars(4);
      if (! fname.numeric_value(_group_is_column) || _group_is_column < 2) {
        cerr << "The Group column directive must be a valid column number > 1\n";
        return 0;
      }
      if (_verbose) {
        cerr << "Group membership is column " << _group_is_column << " of the name\n";
      }
      --_group_is_column;
      if (! GroupsFromNameToken()) {
        cerr << "Options::Initialise:cannot initialise group membership from column " <<
                 (_group_is_column + 1) << " of the name\n";
      }
    } else if (! ReadGroupMembership(fname)) {
      cerr << "Options::Initialise:cannot read group membership '" << fname << "'\n";
      return 1;
    }
  }

  // THis will be an index into the _group array.
  for (uint32_t i = 0; i < _number_fingerprints; ++i) {
    const IWString& g = _fp[i].group();
    if (g.empty()) {
      cerr << "Options::Initialise:no group " << _fp[i].name() << '\n';
      return 0;
    }

    auto iter = _group_name_to_ndx.find(g);
    if (iter == _group_name_to_ndx.end()) {
      uint32_t n = _group_name_to_ndx.size();
      _group_name_to_ndx[g] = n;
      // cerr << "Group '" << g << " index " << n << '\n';
      _group << new Group(g);
    }
  }

  const uint32_t ngroups = _group_name_to_ndx.size();

  if (ngroups < 2) {
    cerr << "Options::Initialise:only " << _group_name_to_ndx.size() << " groups\n";
    return 0;
  }

  _group.resize(ngroups);

  for (uint32_t i = 0; i < _number_fingerprints; ++i) {
    const IWString& g = _fp[i].group();
    auto iter = _group_name_to_ndx.find(g);
    assert(iter != _group_name_to_ndx.end());
    // cerr << "group '" << g << "' found at " << iter->second << " index\n";
    _group[iter->second]->AnotherMember(_fp + i);
  }

  if (cl.option_present('P')) {
    IWString fname = cl.string_value('P');
    if (fname.ends_with(".gfp")) {
      if (! ComputePreviousDistances(fname)) {
        cerr << "Options::Initialise:cannot read '" << fname << "'\n";
        return 0;
      }
    } else if (! ReadPreviousDistance(fname)) {
      cerr << "Options::Initialise:cannot read previously computed distances file '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('D')) {
    IWString fname = cl.string_value('D');
    if (! ReadDesirabilty(fname)) {
      cerr << "Options::Initialise:cannot read desirability data '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _nsel) || _nsel > _group.size()) {
      cerr << "Invalid number to select (-n) must be less than ngroups " <<
              _group.size() << '\n';
      return 0;
    }
  } else {
    _nsel = _group.size();
  }

  if (cl.option_present('w')) {
    IWString w;
    for (int i = 0; cl.value('w', w, i); ++i) {
      if (w == "clear" || w == "zero") {
        _ave_distance_weight = 0.0;
        _desirability_weight = 0.0;
        _min_distance_weight = 0.0;
        _non_short_distance_weight = 0.0;
        _intra_group_weight = 0.0;
      } else if (w.starts_with("avedist=")) {
        w.remove_leading_chars(5);
        if (! w.numeric_value(_ave_distance_weight) || _ave_distance_weight < 0.0) {
          cerr << "Options::Initialise:invalid distance weight '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Average distance weight " << _ave_distance_weight << '\n';
        }
      } else if (w.starts_with("desr")) {
        w.remove_leading_chars(5);
        if (! w.numeric_value(_desirability_weight) || _desirability_weight < 0.0) {
          cerr << "Options::Initialise:invalid desirability weight '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Desirability weight " << _desirability_weight << '\n';
        }
      } else if (w.starts_with("mindist=")) {
        w.remove_leading_chars(8);
        if (! w.numeric_value(_min_distance_weight) || _min_distance_weight < 0.0) {
          cerr << "Options::Initialise:invalid mindistance weight\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Minimum distance weight " << _min_distance_weight << '\n';
        }
      } else if (w.starts_with("nonshort=")) {
        if (! short_distance_cutoff_specified) {
          cerr << "The nonshort= weight directive only makes sense if a short distance cutoff has been specified with the -t option\n";
          return 0;
        }
        w.remove_leading_chars(9);
        if (! w.numeric_value(_non_short_distance_weight) || _non_short_distance_weight < 0.0) {
          cerr << "Options::Initialise:invalid nonshort= qualifier '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Weight give to the non-sort fraction " << _non_short_distance_weight << '\n';
        }
      } else if (w.starts_with("intra=")) {
        w.remove_leading_chars(6);
        if (! w.numeric_value(_intra_group_weight) || _intra_group_weight < 0.0) {
          cerr << "Options::Initialise:invalid intra group ave distance weight '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Average intra group average distance weight " << _intra_group_weight << '\n';
        }
      } else if (w == "help") {
        DisplayDashwOptions(0);
      } else {
        cerr << "Options::Initialise:unrecognised -w qualifier '" << w << "'\n";
        DisplayDashwOptions(1);
      }
    }
  }

  fraction_as_string.set_leading_string(' ');
  if (! fraction_as_string.initialise(0.0f, 1.0f, 4)) {
    cerr << "Cannot initialise fractional output cache\n";
    return 0;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  return 1;
}

int
Options::DetermineInitialDistances() {
  const uint32_t n = _group.size();
  for (uint32_t i = 0; i < n; ++i) {
    Group* gi = _group[i];
    for (uint32_t j = i + 1; j < n; ++j) {
      _group[j]->EstablishInitialDistances(*gi);
    }
  }

  return 1;
}

uint32_t
Options::ReadFingerprints(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  _number_fingerprints = input.count_records_starting_with(_identifier_tag);

  if (0 == _number_fingerprints) {
    cerr << "Options::ReadFingerprints:no occurrences of " << _identifier_tag << "' in input\n";
    return 0;
  }

  _fp = new Fingerprint[_number_fingerprints];

  return ReadFingerprints(input);
}

uint32_t
Options::ReadFingerprints(iwstring_data_source& input) {
  IW_TDT tdt;

  uint32_t ndx = 0;

  for (; tdt.next(input) && ndx < _number_fingerprints; ndx++) {
//  _fp[ndx].set_initial_ndx(ndx);

//  tdt.dataitem_value(smiles_tag, smiles[ndx]);

//  tdt.dataitem_value(identifier_tag, pcn[ndx]);

    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    if (0 == ndx) {
      if (!standard_fingerprints_present()) {
        return 0;
      }
    }

    _fp[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    _fp[ndx].build_iw(gfp[0]);
    _fp[ndx].build_mk(gfp[1]);
    _fp[ndx].build_mk2(gfp[2]);
    _fp[ndx].set_name(gfp.id());

    _id_to_ndx[gfp.id()] = ndx;
  }

  return _number_fingerprints;
}

int
Options::ReadGroupMembership(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadGroupMembership:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadGroupMembership(input);
}

int
Options::ReadGroupMembership(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadGroupMembership:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadGroupMembershipRecord(buffer)) {
      cerr << "Options::ReadGroupMembership:error " << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
Options::ReadGroupMembershipRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, group;
  if (! buffer.nextword(id, i) || ! buffer.nextword(group, i) ||
      id.empty() || group.empty()) {
    return 0;
  }

  const auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadGroupMembershipRecord:no fingerprint for '" << id << "'\n";
    for (const auto& [k, v] : _id_to_ndx) {
      cerr << "key '" << k << "' ndx " << v << '\n';
    }
    return 0;
  }

  _fp[iter->second].set_group(group);

  return 1;
}

int
Options::GroupsFromNameToken() {
  for (uint32_t i = 0; i < _number_fingerprints; ++i) {
    if (! _fp[i].GroupNameIsColumn(_group_is_column)) {
      cerr << "Options::GroupsFromNameToken:error processing '" << _fp[i].name() << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadPreviousDistance(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadPreviousDistance:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadPreviousDistance(input);
}

int
Options::ReadPreviousDistance(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadPreviousDistance:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadPreviousDistanceRecord(buffer)) {
      cerr << "Options::ReadPreviousDistance:error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadPreviousDistanceRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, dist;
  if (! buffer.nextword(id, i) || ! buffer.nextword(dist, i) ||
      id.empty() || dist.empty()) {
    return 0;
  }

  auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadPreviousDistanceRecord:no fingerprint for '" << id << "'\n";
    return 0;
  }

  float d;
  if (! dist.numeric_value(d) || d < 0.0f) {
    cerr << "Options::ReadPreviousDistanceRecord:invalid numeric '" << dist << "'\n";
    return 0;
  }

//_fp[iter->second].set_distance(d);

  const IWString& group = _fp[iter->second].group();
  auto iter2 = _group_name_to_ndx.find(group);
  if (iter2 == _group_name_to_ndx.end()) {
    cerr << "Options::ReadPreviousDistanceRecord:no group for '" << group << "'\n";
    return 0;
  }

  // cerr << "Group " << iter2->second << " gets distance " << d << '\n';
  _group[iter2->second]->AnotherDistance(d);

  return 1;
}

int
Options::ComputePreviousDistances(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ComputePreviousDistances:cannot open '" << fname << "'\n";
    return 0;
  }

  return ComputePreviousDistances(input);
}

int
Options::ComputePreviousDistances(iwstring_data_source& input) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Options::ComputePreviousDistances:Cannot read fingerprint\n";
      return 0;
    }

    GFP_Standard sfp;
    sfp.build_molecular_properties(gfp.molecular_properties_integer());
    sfp.build_iw(gfp[0]);
    sfp.build_mk(gfp[1]);
    sfp.build_mk2(gfp[2]);

    for (uint32_t i = 0; i < _number_fingerprints; ++i) {
      float d = sfp.tanimoto_distance(_fp[i]);
      auto iter = _group_name_to_ndx.find(_fp[i].group());
      if (iter == _group_name_to_ndx.end()) [[ unlikely ]] {  // Impossible.
        cerr << "Options::ComputePreviousDistances:no group for " << _fp[i].group() << '\n';
        return 0;
      }

      _group[iter->second]->AnotherDistance(d);
    }
  }

  return 1;
}

int
Options::ReadDesirabilty(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadDesirabilty:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDesirabilty(input);
}

int
Options::ReadDesirabilty(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadDesirabilty:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadDesirabiltyRecord(buffer)) {
      cerr << "Options::ReadDesirabilty:error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadDesirabiltyRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, desirability;
  if (! buffer.nextword(id, i) || ! buffer.nextword(desirability, i) ||
      id.empty() || desirability.empty()) {
    return 0;
  }

  auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadDesirabiltyRecord:no fingerprint for '" << id << "'\n";
    return 0;
  }

  float d;
  if (! desirability.numeric_value(d) || d < 0.0f) {
    cerr << "Options::ReadDesirabiltyRecord:invalid desirability '" << desirability << "\n";
    return 0;
  }

  const IWString& group = _fp[iter->second].group();
  auto iter2 = _group_name_to_ndx.find(group);

  _group[iter2->second]->AnotherDesirability(d);

  return 1;
}

int
Options::Spread(IWString_and_File_Descriptor & output) {
  uint32_t nsel = 0;
  if (_verbose) {
    cerr << "Will select " << _nsel << " groups\n";
  }

  for ( ; nsel < _nsel; ++nsel) {
    Group* sel = GetNextGroup();
    // cerr << "To select " << sel << '\n';
    if (sel == nullptr) {
      return 1;
    }

    if (_verbose) {
      cerr << "Selected group " << sel->name() << '\n';
    }

    for (Group* g : _group) {
      if (g == sel) {
        continue;
      }
      g->GroupHasBeenSelected(sel);
    }

    sel->set_selected(nsel + 1);
    sel->Write(output);
  }

  // This is not necessary, groups were written as they formed.
  _group.iwqsort_lambda([](const Group* g1, const Group* g2) {
      if (g1->selected() < g2->selected()) {
        return -1;
      }
      if (g1->selected() > g2->selected()) {
        return 1;
      }
      return 0;
    });

  return 1;
}

double
Options::Score(Group& g) const {
  double rc = 0.0;
  const Accumulator<double>& acc_dist = g.DistanceAccumulator();

  if (_ave_distance_weight > 0 && acc_dist.n() > 0) {
    rc += _ave_distance_weight * acc_dist.average();
  }

  if (_desirability_weight > 0) {
    rc += _desirability_weight * g.AverageDesirability();
  }

  if (_min_distance_weight > 0.0 && acc_dist.n() > 0) {
    rc += _min_distance_weight + acc_dist.minval();
  }

  if (_non_short_distance_weight > 0.0) {
    rc += _non_short_distance_weight + g.FractionNotShort();
  }

  if (_intra_group_weight > 0) {
    rc += _intra_group_weight * g.IntraGroupAverageDistance();
  }

  return rc;
}

Group*
Options::GetNextGroup() const {
  double highest_score = -1.0;
  Group* highest_score_group = nullptr;

  // cerr << "GetNextGroup\n";
  for (Group* g : _group) {
    // cerr << " group? " << g->name() << " sel " << g->selected() << '\n';
    if (g->selected()) {
      continue;
    }

    double score = Score(*g);
    // cerr << "Score for group " << score << '\n';
    if (score > highest_score) {
      highest_score = score;
      highest_score_group = g;
    }
  }

  if (highest_score_group == nullptr) {
    return nullptr;
  }

  return highest_score_group;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vP:G:w:n:D:T:t:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    Usage(2);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! options.Spread(output)) {
    cerr << "gfp_group_spread failure\n";
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_group_spread

int
main(int argc, char** argv) {
  return gfp_group_spread::Main(argc, argv);
}
