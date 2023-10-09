// Converts a nearest neighbours file to csv

#include <cctype>
#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"

namespace nn2csv {
using std::cerr;

IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";
IWString distance_tag = "DIST<";
IWString nbr_tag = "NBR<";

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Converts a .nn file to csv form\n";
  cerr << " -o <sep>      set token separator (default ,)\n";
  cerr << " -n <nbrs>     only write <nbrs> neighbours\n";
  cerr << " -s ...        sort the targets, enter '-s help' for info\n";
  cerr << " -z            do not write molecules with no neighbours\n";
  cerr << " -l            strip trailing neighbour count from needle names\n";
  cerr << " -v            verbose output\n";
// clang-format on

  exit(rc);
}

class Neighbour {
 private:
  IWString _smiles;
  IWString _id;
  float _distance;
  int _nbr_index;

 public:
  Neighbour();

  // If `reading_nbr_indices` is set, look for nbr_tag in the input
  // and store in `_nbr_index`.
  int Build(iwstring_data_source& input, int reading_nbr_indices);

  float distance() const {
    return _distance;
  }

  int nbr_index() const {
    return _nbr_index;
  }

  int Write(char sep, IWString_and_File_Descriptor& output) const;
};

Neighbour::Neighbour() {
  _distance = -1.0f;
  _nbr_index = -1;
}

int
Neighbour::Write(char sep, IWString_and_File_Descriptor& output) const {
  output << _smiles;
  output << sep << _id;
  output << sep << _distance;
  return 1;
}

// We assume that `input` is pointing at the start of a neighbour
// and we expect to read smiles,id,distance, but we do not assume
// that order.
int
Neighbour::Build(iwstring_data_source& input,
                 int reading_nbr_indices) {
  const_IWSubstring buffer;

  // We need 3 components initialised, and then we are done.
  int components_initialised = 0;
  while (input.next_record(buffer)) {
    // cerr << "Neighbour::Build:reading '" << buffer << "'\n";
    if (_smiles.empty() && buffer.starts_with(smiles_tag)) {
      buffer.remove_leading_chars(smiles_tag.length());
      buffer.chop();
      _smiles = buffer;
      ++components_initialised;
    } else if (_id.empty() && buffer.starts_with(identifier_tag)) {
      buffer.remove_leading_chars(identifier_tag.length());
      buffer.chop();
      _id = buffer;
      ++components_initialised;
    } else if (_distance < 0.0f && buffer.starts_with(distance_tag)) {
      buffer.remove_leading_chars(distance_tag.length());
      buffer.chop();
      if (!buffer.numeric_value(_distance) || _distance < 0.0f || _distance > 1.0f) {
        cerr << "Neighbour::Build:invalid distance '" << buffer << "'\n";
        return 0;
      }
      ++components_initialised;
    } else if (reading_nbr_indices && _nbr_index < 0 && buffer.starts_with(nbr_tag)) {
      buffer.remove_leading_chars(nbr_tag.length());
      buffer.chop();
      if (! buffer.numeric_value(_nbr_index) || _nbr_index < 0) {
        cerr << "Neighbour::Build:invalid neighbour index '" << buffer << "'\n";
        return 0;
      }
    } else {
      cerr << "Neighbour::Build:invalid input, line " << input.lines_read() << '\n';
      cerr << reading_nbr_indices << ' ' << _nbr_index << '\n';
      cerr << buffer << '\n';
      return 0;
    }

    // cerr << reading_nbr_indices << ' ' << _smiles << ' ' << _id << " " << _distance << ' ' << components_initialised << '\n';

    // Reading indices, may have smiles, but must have ndx and distance.
    if (reading_nbr_indices) {
      if (_nbr_index >= 0 && _distance >= 0.0) {
        return 1;
      }
      continue;
    }

    // Reading normal output, must have smiles, id, distance
    if (components_initialised < 3) {
      continue;
    }
    
    if (_smiles.length() > 0 && _id.length() > 0 && _distance >= 0.0) {
      return 1;
    }

    continue;
  }

  cerr << "Neighbour::Build:read beyond nndata\n";
  return 0;
}

// An individual entry in the NN file.
// smiles, id and a list of neighbours.
class NNData {
 private:
  IWString _smiles;
  IWString _id;
  resizable_array_p<Neighbour> _neighbour;

  // private functions.

 public:
  // If `reading_nbr_indices` is set, look for nbr_tag and store with nbrs.
  int Build(iwstring_data_source& input, int reading_nbr_indices);

  int number_neighbours() const {
    return _neighbour.number_elements();
  }

  const Neighbour* operator[](int ndx) const {
    return _neighbour[ndx];
  }

  // Many tools add the number of neighbours to the name of the needle.
  // Remove that if present.
  int RemoveTrailingNumberFromName();

  int Write(int max_nbrs, char sep, IWString_and_File_Descriptor& output) const;
  int WriteIndices(int max_nbrs, char sep, IWString_and_File_Descriptor& output) const;

};

int
FetchRecord(iwstring_data_source& input, const IWString& tag, IWString& destination) {
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "FetchRecord:cannot fetch\n";
    return 0;
  }
  if (!buffer.starts_with(tag)) {
    cerr << "FetchRecord:tag mismatch, expected '" << tag << "' got " << buffer << '\n';
    return 0;
  }

  buffer.remove_leading_chars(tag.length());
  buffer.chop();
  destination = buffer;

  return 1;
}

int
NNData::Build(iwstring_data_source& input,
              int reading_nbr_indices) {
  if (!FetchRecord(input, smiles_tag, _smiles) ||
      !FetchRecord(input, identifier_tag, _id)) {
    cerr << "NNData::Build:cannot read first two records\n";
    return 0;
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer == '|') {
      return 1;
    }
    input.push_record();
    std::unique_ptr<Neighbour> nbr = std::make_unique<Neighbour>();
    if (!nbr->Build(input, reading_nbr_indices)) {
      cerr << "NNData::Build:fatal error at line " << input.lines_read() << '\n';
      return 0;
    }
    _neighbour << nbr.release();
  }

  return 1;
}

// Write the neighbour data. Tokens separated by `sep`
// and we must produce `max_nbrs` columns.
int
NNData::Write(int max_nbrs, char sep, IWString_and_File_Descriptor& output) const {
  static constexpr char kMissing = '*';

  output << _smiles;
  output << sep << _id;
  output << sep << max_nbrs;

  Accumulator<float> acc;
  for (const Neighbour* nbr : _neighbour) {
    acc.extra(nbr->distance());
  }
  if (acc.n() > 0) {
    output << sep << static_cast<float>(acc.average());
    output << sep << acc.maxval();
  } else {
    output << sep << kMissing;
    output << sep << kMissing;
  }

  int nwrite = max_nbrs;
  if (nwrite > _neighbour.number_elements()) {
    nwrite = _neighbour.number_elements();
  }

  for (int i = 0; i < nwrite; ++i) {
    output << sep;
    _neighbour[i]->Write(sep, output);
  }

  for (int i = nwrite; i < max_nbrs; ++i) {
    output << sep << kMissing;  // smiles
    output << sep << kMissing;  // id
    output << sep << kMissing;  // distance
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
NNData::WriteIndices(int max_nbrs, char sep, IWString_and_File_Descriptor& output) const {
  static constexpr char kMissing = '*';

  output << _id;

  int nwrite = max_nbrs;
  if (nwrite > _neighbour.number_elements()) {
    nwrite = _neighbour.number_elements();
  }

  for (int i = 0; i < nwrite; ++i) {
    output << sep << _neighbour[i]->nbr_index();
  }

  // Fill missing columns (if any).
  for (int i = nwrite; i < max_nbrs; ++i) {
    output << sep << kMissing;
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
NNData::RemoveTrailingNumberFromName() {
  // Min case would be 'a n'.
  if (_id.length() < 3) {
    return 0;
  }

  for (int i = _id.length() - 1; i > 0; --i) {
    char c = _id[i];
    if (std::isspace(c)) {
      _id.iwtruncate(i);
      return 1;
    }
    if (!std::isdigit(c)) {
      return 0;
    }
  }

  return 0;
}

enum class SortType {
  kNone = 0,
  kShortestDistance = 1,
  kLongestDistance = 2,
  kNumberNeighbours = 3,
};

class NN2Csv {
 private:
  int _verbose;
  // The output separator.
  char _sep;

  // We can limit the number of neighbours output
  int _neighbours_to_write;

  // All the nearest neighbour data read from the input.
  resizable_array_p<NNData> _nn_data;

  SortType _sort_type;

  // If the -z option is present, do not write a record
  // if there are zero neighbours.
  int _write_molecules_with_zero_neighbours;

  // Many utilities append the number of neighbours to the id of the needle.
  int _strip_trailing_neighbour_count_in_name;

  int _reading_nbr_indices;

  // private functions.

  int Sort(SortType& s);
  int SortByNumberNeighbours();
  int SortByClosestDist();

  int WriteIndices(int max_nbrs, IWString_and_File_Descriptor& output) const;

 public:
  NN2Csv();

  int Initialise(Command_Line& cl);

  int Accumulate(const char* fname);
  int Accumulate(iwstring_data_source& input);

  int number_targets() const {
    return _nn_data.number_elements();
  }

  int SortIfRequested() {
    return Sort(_sort_type);
  }

  int Write(IWString_and_File_Descriptor& output) const;
};

NN2Csv::NN2Csv() {
  _sep = ',';
  _sort_type = SortType::kNone;
  _neighbours_to_write = std::numeric_limits<int>::max();
  _write_molecules_with_zero_neighbours = 1;
  _strip_trailing_neighbour_count_in_name = 0;
  _reading_nbr_indices = 0;
}

void
DisplaySortOptions(std::ostream& output) {
  output << " -s nbrs      sort by the number of neighbours\n";
  output << " -s dist      sort by the closest distance\n";
  exit(0);
}

int
NN2Csv::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    char_name_to_char(o);
    _sep = o[0];
    if (_verbose) {
      cerr << "Token separator '" << _sep << "'\n";
    }
  }

  if (cl.option_present('s')) {
    const IWString s = cl.string_value('s');
    if (s == "none") {
      _sort_type = SortType::kNone;
    } else if (s == "nbrs") {
      _sort_type = SortType::kNumberNeighbours;
    } else if (s == "dist") {
      _sort_type = SortType::kShortestDistance;
    } else if (s == "maxdist") {
      _sort_type = SortType::kLongestDistance;
    } else if (s == "help") {
      DisplaySortOptions(cerr);
    } else {
      cerr << "NN2Csv:Initialise:unrecognised sort type '" << s << "'\n";
      DisplaySortOptions(cerr);
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', _neighbours_to_write) || _neighbours_to_write < 0) {
      cerr << "NN2Csv::Initialise:the number of neighbours to write (-n) must be a whole "
              "+ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will write a max of " << _neighbours_to_write << " neighbours\n";
    }
  }

  if (cl.option_present('z')) {
    _write_molecules_with_zero_neighbours = 0;
    if (_verbose) {
      cerr << "MOlecules with zero neighbours will not be written\n";
    }
  }

  if (cl.option_present('x')) {
    const_IWSubstring x = cl.string_value('x');
    if (x == '.') {
    } else {
      nbr_tag = x;
      if (! nbr_tag.ends_with('<')) {
        nbr_tag << '<';
      }
    }

    _reading_nbr_indices = 1;
    if (_verbose) {
      cerr << "Processing input that contains nbr indices in '" << nbr_tag << "'\n";
    }
  }

  if (cl.option_present('l')) {
    _strip_trailing_neighbour_count_in_name = 1;
  }

  return 1;
}

int
NN2Csv::Accumulate(const char* fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "NN2Csv::Accumulate:cannot open '" << fname << "'\n";
    return 0;
  }

  return Accumulate(input);
}

int
NN2Csv::Accumulate(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    std::unique_ptr<NNData> nn = std::make_unique<NNData>();
    input.push_record();

    if (!nn->Build(input, _reading_nbr_indices)) {
      cerr << "Cannot fetch start of nbr list\n";
      return 0;
    }

    if (_strip_trailing_neighbour_count_in_name) {
      nn->RemoveTrailingNumberFromName();
    }

    _nn_data << nn.release();
  }

  if (_nn_data.empty()) {
    cerr << "NN2Csv:Accumulate:no data\n";
    return 0;
  }

  return 1;
}

int
NN2Csv::Sort(SortType& sort_type) {
  switch (sort_type) {
    case SortType::kNone:
      return 1;
    case SortType::kNumberNeighbours:
      return SortByNumberNeighbours();
    case SortType::kShortestDistance:
      return SortByClosestDist();
    default:
      cerr << "Unimplemented sort\n";
      return 0;
  }
}

int
NN2Csv::SortByNumberNeighbours() {
  _nn_data.iwqsort_lambda([](const NNData* n1, const NNData* n2) {
    if (n1->number_neighbours() < n2->number_neighbours()) {
      return -1;
    } else if (n1->number_neighbours() > n2->number_neighbours()) {
      return 1;
    } else {
      return 0;
    }
  });

  return 1;
}

// Sorting by closest neighbour is complicate by the fact
// that there may be no neighbours.
int
NN2Csv::SortByClosestDist() {
  _nn_data.iwqsort_lambda([](const NNData* n1, const NNData* n2) {
    float d1 = std::numeric_limits<float>::max();
    if (n1->number_neighbours() > 0) {
      d1 = (*n1)[0]->distance();
    }
    float d2 = std::numeric_limits<float>::max();
    if (n2->number_neighbours() > 0) {
      d2 = (*n2)[0]->distance();
    }
    if (d1 < d2) {
      return -1;
    } else if (d1 > d2) {
      return 1;
    } else {
      return 0;
    }
  });

  return 1;
}

int
NN2Csv::Write(IWString_and_File_Descriptor& output) const {
  int max_nbrs = 0;
  for (const NNData* nn_data : _nn_data) {
    if (nn_data->number_neighbours() > max_nbrs) {
      max_nbrs = nn_data->number_neighbours();
    }
  }

  if (_verbose) {
    cerr << "Across " << _nn_data.size() << " targets max nbr count " << max_nbrs << '\n';
  }

  if (max_nbrs > _neighbours_to_write) {
    max_nbrs = _neighbours_to_write;
    if (_verbose) {
      cerr << "Output trimmed to " << max_nbrs << " via command line\n";
    }
  }

  if (_reading_nbr_indices) {
    return WriteIndices(max_nbrs, output);
  }

  output << "smiles";
  output << _sep << "id";
  output << _sep << "nbrs";
  output << _sep << "avedist";
  output << _sep << "maxdist";
  constexpr char kJoin = '_';
  for (int i = 0; i < max_nbrs; ++i) {
    output << _sep << "smiles" << kJoin << (i + 1);
    output << _sep << "id" << kJoin << (i + 1);
    output << _sep << "dist" << kJoin << (i + 1);
  }

  output << '\n';

  int targets_with_no_neighbours = 0;

  for (const NNData* nn_data : _nn_data) {
    if (!_write_molecules_with_zero_neighbours && nn_data->number_neighbours() == 0) {
      ++targets_with_no_neighbours;
      continue;
    }
    nn_data->Write(max_nbrs, _sep, output);
  }

  if (_verbose && !_write_molecules_with_zero_neighbours) {
    cerr << "Skipped " << targets_with_no_neighbours << " targets with no neighbours\n";
  }

  return 1;
}

int
NN2Csv::WriteIndices(int max_nbrs, IWString_and_File_Descriptor& output) const {
  output << "Id";
  for (int i = 0; i < max_nbrs; ++i) {
    output << _sep << "nbr" << i;
  }
  output << '\n';

  int targets_with_no_neighbours = 0;

  for (const NNData* nn_data: _nn_data) {
    if (!_write_molecules_with_zero_neighbours && nn_data->number_neighbours() == 0) {
      ++targets_with_no_neighbours;
      continue;
    }
    nn_data->WriteIndices(max_nbrs, _sep, output);
  }

  if (_verbose && !_write_molecules_with_zero_neighbours) {
    cerr << "Skipped " << targets_with_no_neighbours << " targets with no neighbours\n";
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vo:s:n:zx:l");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    Usage(1);
  }

  NN2Csv nn_to_csv;
  if (!nn_to_csv.Initialise(cl)) {
    cerr << "Cannot initialise NN2Csv conditions\n";
    return 1;
  }

  for (const char* fname : cl) {
    if (!nn_to_csv.Accumulate(fname)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read data on " << nn_to_csv.number_targets() << " targets\n";
  }

  IWString_and_File_Descriptor output(1);

  nn_to_csv.SortIfRequested();

  nn_to_csv.Write(output);

  return 0;
}

}  // namespace nn2csv

int
main(int argc, char** argv) {
  int rc = nn2csv::Main(argc, argv);
  return rc;
}
