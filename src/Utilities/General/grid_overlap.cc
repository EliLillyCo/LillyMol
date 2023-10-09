// Computes the overlap between two grids.
// Grids are read in as x y z value' tuples.
// They are not aligned. Form an overlap
// between them.

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/compile_time.h"
#include "Foundational/iwmisc/misc.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"

namespace grid_overlap {

using std::cerr;

// The grid is read in as (probably) unordered xyz tuples.
// We discern the spacing and re-order the entries in _value
// so they can be indexed in a regular form. Note that if
// the grid comes in sorted, this is probably not necessary.
class Grid {
  private:
    // the data as read in.
    uint32_t _npoints;
    float* _x;
    float* _y;
    float* _z;
    float* _value;

    // When looking at inter-grid spacing, we need a tolerance for
    // when values look the same. this is implemented as an absolute
    // tolerance.
    float _tol;

    // When finding the location of coordinates in the grid, we
    // need to know the dimensions of the grid.
    float _xmin, _xmax;
    float _ymin, _ymax;
    float _zmin, _zmax;

    // When we convert to indexed form, we need these values
    double _dx, _dy, _dz;
    // the integer dimensions of the grid.
    uint32_t _nx, _ny, _nz;

    // We store the file name.
    IWString _name;

  // private functions
    void DetermineExtremities();
    int DetermineDyXyz();
    int ToIndexedForm();

    // Given coordinates, where is that point in this grid.
    // Returns nullopt if the point is out of range.
    std::optional<uint32_t> Index(float x, float y, float z) const;

  public:
    Grid();
    ~Grid();

    int Build(const char * fname);
    int Build(iwstring_data_source& input);

    uint32_t npoints() const {
      return _npoints;
    }

    void set_tolerance(float s) {
      _tol = s;
    }

    std::array<float, 6> Extremities() const;

    std::array<float, 3> DxDyDz() const;

    // Return the value of the point closest to [x,y,z], If out of
    // range, return nullopt.
    std::optional<float> ValueClosest(float x, float y, float z) const;

    // Identify the 6 points that enclose `x,y,z` and return the mean
    // value.
    std::optional<float> MeanNbrValue(float x, float y, float z) const;

    // For each of our points, fetch the corresponding value from `rhs`
    // and pass that to `fn`. Note that std::nullopt will be passed
    // for those points of ours that are not in range of `rhs`.
    template <typename T> void CompareViaClosest(const Grid& rhs, T& fn) const;
    // Same idea, but this time the comparison is via the average of
    // the closest point and the neighbours.
    template <typename T> void CompareViaNNs(const Grid& rhs, T& fn) const;

    void set_name(const char* s) {
      _name = s;
    }
    void set_name(const const_IWSubstring& s) {
      _name = s;
    }
    const IWString& name() const {
      return _name;
    }
};

Grid::Grid() {
  _x = nullptr;
  _y = nullptr;
  _z = nullptr;
  _value = nullptr;

  _tol = std::numeric_limits<float>::epsilon();

  _xmin = std::numeric_limits<float>::max();
  _xmax = -std::numeric_limits<float>::max();
  _ymin = std::numeric_limits<float>::max();
  _ymax = -std::numeric_limits<float>::max();
  _zmin = std::numeric_limits<float>::max();
  _zmax = -std::numeric_limits<float>::max();

  _dx = 0.0;
  _dy = 0.0;
  _dz = 0.0;
  _nx = 0;
  _ny = 0;
  _nz = 0;
}

Grid::~Grid() {
  delete [] _x;
  delete [] _y;
  delete [] _z;
  delete [] _value;
}

int
Grid::Build(const char * fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Grid::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
Parse(const const_IWSubstring& buffer,
      float& x,
      float& y,
      float& z,
      float& value) {

  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); ++col) {
    float f;
    if (! token.numeric_value(f)) {
      cerr << "Parse:invalid float '" << token << "'\n";
      return 0;
    }
    if (col == 0) {
      x = f;
    } else if (col == 1) {
      y = f;
    } else if (col == 2) {
      z = f;
    } else if (col == 3) {
      value = f;
    }
  }

  return 1;
}

int
Grid::Build(iwstring_data_source& input) {
  _npoints = input.records_remaining();
  if (_npoints == 0) {
    cerr << "Grid::Build:empty file\n";
    return 0;
  }

  _x = new float[_npoints];
  _y = new float[_npoints];
  _z = new float[_npoints];
  _value = new float[_npoints];

  const_IWSubstring buffer;
  for (int line = 0; input.next_record(buffer); ++line) {
    if (! Parse(buffer, _x[line], _y[line], _z[line], _value[line])) {
      cerr << "Grid::Build:cannot parse '" << buffer << "'\n";
      return 0;
    }
  }

  uint32_t nonzero_values = 0;
  for (uint32_t i = 0; i < _npoints; ++i) {
    if (_value[i] > 0.0) {
      ++nonzero_values;
    }
  }

  cerr << nonzero_values << " of " << _npoints << " non zero "
       << iwmisc::Fraction<float>(nonzero_values, _npoints) << '\n';

  DetermineExtremities();
  if (! DetermineDyXyz()) {
    cerr << "Grid::Build:cannot determine dx, dy, dz\n";
    return 0;
  }

  if (! ToIndexedForm()) {
    cerr << "Grid::Build:canot make indexed form\n";
    return 0;
  }

  return 1;
}

void
Grid::DetermineExtremities() {
  _xmin = _x[0];
  _xmax = _x[0];
  _ymin = _y[0];
  _ymax = _y[0];
  _zmin = _z[0];
  _zmax = _z[0];
  for (uint32_t i = 1; i < _npoints; ++i) {
    if (_x[i] < _xmin) {
      _xmin = _x[i];
    } else if (_x[i] > _xmax) {
      _xmax = _x[i];
    }
    if (_y[i] < _ymin) {
      _ymin = _y[i];
    } else if (_y[i] > _ymax) {
      _ymax = _y[i];
    }
    if (_z[i] < _zmin) {
      _zmin = _z[i];
    } else if (_z[i] > _zmax) {
      _zmax = _z[i];
    }
  }
}

std::array<float, 6>
Grid::Extremities() const {
  std::array<float, 6> result;
  result[0] = _xmin;
  result[1] = _xmax;
  result[2] = _ymin;
  result[3] = _ymax;
  result[4] = _zmin;
  result[5] = _zmax;

  return result;
}

// Determining the inter layer separation is tricky.
// Sort `values` and then examine the differences between
// adjacent values. If they are within `near_zero`, we
// assume in the same layer. Look at all the diffs that
// are not within the same layer, and hopefully we find
// a single value.
std::optional<float>
DetermineSeparation(const float* values,
                    uint32_t nvalues,
                    float near_zero,
                    float* storage) {
  std::copy_n(values, nvalues, storage);
  std::sort(storage, storage + nvalues, [] (float v1, float v2) {
    return v1 < v2;
  });

  Accumulator<double> diffs;
  // Load `storage` with the diffs.
  uint32_t ndx = 0;
  for (uint32_t i = 0; i < (nvalues - 1); ++i) {
    double d = storage[i + 1] - storage[i];
    if (d < near_zero) {
      continue;
    }
    storage[ndx] = d;
    diffs.extra(d);
    ++ndx;
  }

  // If the diffs are all near each other, we are done.
  if (diffs.range() < near_zero) {
    return static_cast<float>(diffs.average());
  }

  cerr << "DetermineSeparation:possibly irregular grid spacing\n";
  cerr << "Range of diffs " << diffs.range() << " tolerance allowed " << near_zero << '\n';

  // Sort the diffs
  std::sort(storage, storage + ndx, [] (float v1, float v2) {
    return v1 < v2;
  });

  cerr << "Find " << ndx << " differences\n";
  for (uint32_t i = 0; i < ndx; ++i) {
    cerr << storage[i];
    if (i > 0 && (storage[i] - storage[i-1] > near_zero)) {
      cerr << " *";
    }
    cerr << '\n';
  }

  return std::nullopt;
}

// For each axis, deterine the separation.
int
Grid::DetermineDyXyz() {
  std::unique_ptr<float[]> storage = std::make_unique<float[]>(_npoints);
  std::optional<float> sx = DetermineSeparation(_x, _npoints, _tol, storage.get());
  std::optional<float> sy = DetermineSeparation(_y, _npoints, _tol, storage.get());
  std::optional<float> sz = DetermineSeparation(_z, _npoints, _tol, storage.get());

  if (! sx || ! sy || !sz) {
    cerr << "Grid::DxDyDz:cannot determine inter grid spacing\n";
    return 0;
  }

  _dx = *sx;
  _dy = *sy;
  _dz = *sz;

  return 1;
}

std::array<float, 3>
Grid::DxDyDz() const {
  std::array<float, 3> result;
  result[0] = _dx;
  result[1] = _dy;
  result[2] = _dz;

  return result;
}

std::optional<uint32_t>
Grid::Index(float x, float y, float z) const {
  if (x < _xmin || x > _xmax) {
    return std::nullopt;
  }
  if (y < _ymin || y > _ymax) {
    return std::nullopt;
  }
  if (z < _zmin || z > _zmax) {
    return std::nullopt;
  }

  uint32_t ix = static_cast<uint32_t>((x - _xmin) / _dx + 0.49999);
  uint32_t iy = static_cast<uint32_t>((y - _ymin) / _dy + 0.49999);
  uint32_t iz = static_cast<uint32_t>((z - _zmin) / _dz + 0.49999);

  uint32_t ndx = ix * (_ny * _nz) + iy * _nz + iz;

  return ndx;
}

std::optional<float>
Grid::ValueClosest(float x, float y, float z) const {
  std::optional<uint32_t> ndx = Index(x, y, z);
  if (! ndx) {
    return std::nullopt;
  }

  return _value[*ndx];
}

std::optional<float>
Grid::MeanNbrValue(float x, float y, float z) const {
  Accumulator<double> result;
  std::optional<uint32_t> ndx = Index(x, y, z);
  if (!  ndx) {
    return std::nullopt;
  }
  std::optional<uint32_t> north = Index(x, y + _dy, z);
  if (north) {
    result.extra(_value[*north]);
  }
  std::optional<uint32_t> east = Index(x + _dx, y, z);
  if (east) {
    result.extra(_value[*east]);
  }
  std::optional<uint32_t> south = Index(x, y - _dy, z);
  if (south) {
    result.extra(_value[*south]);
  }
  std::optional<uint32_t> west = Index(x - _dx, y, z);
  if (west) {
    result.extra(_value[*west]);
  }
  std::optional<uint32_t> top = Index(x, y, z + _dz);
  if (top) {
    result.extra(_value[*top]);
  }
  std::optional<uint32_t> bot = Index(x, y, z - _dz);
  if (bot) {
    result.extra(_value[*bot]);
  }

  return result.average();
}

// Assuming that _dx, _dy, _dz and the range has been determined,
// place the points at their proper index in _values.
int
Grid::ToIndexedForm() {
  // first check that the number of points is OK.
  _nx = static_cast<uint32_t> ((_xmax - _xmin) / _dx + 0.4999) + 1;
  _ny = static_cast<uint32_t> ((_ymax - _ymin) / _dy + 0.4999) + 1;
  _nz = static_cast<uint32_t> ((_zmax - _zmin) / _dz + 0.4999) + 1;
  uint32_t calculated_npoints = _nx * _ny * _nz;
  if (calculated_npoints != _npoints) {
    cerr << "Grid::ToIndexedForm:npoints mismatch. Calculalte "
         << calculated_npoints << " but have " << _npoints << '\n';
    return 0;
  }

  // Keep track of whether or not each index is assigned.
  std::unique_ptr<int[]> assigned(new_int(_npoints));
  std::unique_ptr<float[]> new_values(new_float(_npoints));
  for (uint32_t i = 0; i < _npoints; ++i) {
    std::optional<uint32_t> ndx = Index(_x[i], _y[i], _z[i]);
    // this should not happen.
    if (! ndx) {
      cerr << "Grid::ToIndexedForm:out of range?\n";
      return 0;
    }
    if (assigned[*ndx]) {
      cerr << "Grid::ToIndexedForm:index " << *ndx << " already assigned, i = " << i << '\n';
      return 0;
    }
    assigned[*ndx] = 1;
    new_values[*ndx] = _value[i];
  }

  int not_assigned = 0;
  for (uint32_t i = 0; i < _npoints; ++i) {
    if (! assigned[i]) {
      cerr << "No value assigned " << i << '\n';
      ++not_assigned;
    }
  }

  if (not_assigned) {
    cerr << not_assigned << " values in " << _npoints << " grid not assigned\n";
    return 0;
  }

  std::copy_n(new_values.get(), _npoints, _value);

  return 1;
}

// Comparison between `this` and `rhs` where `fn` accumulates
// the degree of overlap based on the value of the closest
// point in `rhs`.
template <typename T>
void
Grid::CompareViaClosest(const Grid& rhs, T& fn) const {
  for (uint32_t i = 0; i < _npoints; ++i) {
    fn(_value[i], rhs.ValueClosest(_x[i], _y[i], _z[i]));
  }
}

// Comparison between `this` and `rhs` where `fn` accumulates
// the degree of overlap based on mean neighbour values in `rhs`.
template <typename T>
void
Grid::CompareViaNNs(const Grid& rhs, T& fn) const {
  for (uint32_t i = 0; i < _npoints; ++i) {
    fn(_value[i], rhs.MeanNbrValue(_x[i], _y[i], _z[i]));
  }
}

// When comparing against multiple grids, we need to keep track
// of the grid number and the associated value.
// When comparing against multiple grids, we need to keep track
// of the grid number and the associated valuee
struct GridNumberAndValue {
  int grid_number;
  float value;
};

// A query `query_name` has been compared against all grids in
// `comparison_grids`. The results are stored in `grid_and_value`.
// Sort the results and write to `outupt`.
int
WriteNeighbours(const IWString& query_name,
                const resizable_array_p<Grid>& comparison_grids,
                GridNumberAndValue* grid_and_value,
                IWString_and_File_Descriptor& output) {
  const int ncomparison = comparison_grids.size();

  // Sort by similarity, higher values -> more similar.
  std::sort(grid_and_value, grid_and_value + ncomparison,
      [](const GridNumberAndValue& gv1, const GridNumberAndValue& gv2) {
        return gv1.value > gv2.value;
  });

  nnbr::NearNeighbours proto;
  proto.set_name(query_name.data(), query_name.length());
  for (int i = 0; i < ncomparison; ++i) {
    auto* nbr = proto.add_nbr();
    const int grid_num = grid_and_value[i].grid_number;
    const float tanimoto = grid_and_value[i].value;
    const IWString& grid_name = comparison_grids[grid_num]->name();
    nbr->set_id(grid_name.data(), grid_name.length());
    nbr->set_dist(tanimoto);
  }

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (! printer.PrintToString(proto, &buffer)) {
    cerr << "WriteNeighbours:cannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }
  output << buffer;
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

// Generalized comparison of one or more query grids against
// one or more comparison grids.
// Loop over `query_grids` and treat each individually.
// Comparison function is CompareViaNNs
template <typename T>
int
CompareViaNNs(const resizable_array_p<Grid>& comparison_grids,
              const resizable_array_p<Grid>& query_grids,
              T& fn,
              IWString_and_File_Descriptor& output) {
  const int ncomparison = comparison_grids.size();
  std::unique_ptr<GridNumberAndValue[]> grid_value =
     std::make_unique<GridNumberAndValue[]>(ncomparison);

  // Compare each query against all of the comparison grids.
  for (const Grid* q : query_grids) {
    for (int i = 0; i < ncomparison; ++i) {
      fn.Reset();
      q->CompareViaNNs(*comparison_grids[i], fn);
      grid_value[i].grid_number = i;
      grid_value[i].value = fn.Tanimoto();
    }
    WriteNeighbours(q->name(), comparison_grids, grid_value.get(), output);
  }

  return 1;
}

// Generalized comparison of one or more query grids against
// one or more comparison grids.
// Loop over `query_grids` and treat each individually.
// Comparison function is CompareViaClosest
template <typename T>
int
CompareViaClosest(const resizable_array_p<Grid>& comparison_grids,
              const resizable_array_p<Grid>& query_grids,
              T& fn,
              IWString_and_File_Descriptor& output) {
  const int ncomparison = comparison_grids.size();
  std::unique_ptr<GridNumberAndValue[]> grid_value =
     std::make_unique<GridNumberAndValue[]>(ncomparison);

  // Compare each query against all of the comparison grids.
  for (const Grid* q : query_grids) {
    for (int i = 0; i < ncomparison; ++i) {
      fn.Reset();
      q->CompareViaClosest(*comparison_grids[i], fn);
      grid_value[i].grid_number = i;
      grid_value[i].value = fn.Tanimoto();
    }
    WriteNeighbours(q->name(), comparison_grids, grid_value.get(), output);
  }

  return 1;
}

// Functor class for Grid::Compare.
// Accumulates comparison results between two grids.
class CompareGrids {
  private:
    // the number of comparisons done.
    uint32_t _comparisons;

    // the number where the rhs value is in range.
    uint32_t _in_range;

    // Accumulators for lhs, rhs and the signed differences between them.
    Accumulator<double> _acc_lhs;
    Accumulator<double> _acc_rhs;
    Accumulator<double> _acc_diffs;
    // For computing a Tanimoto coefficient.
    double _in_common;
    // the number of times rhs is greater than lhs.
    uint32_t _rhs_larger;

  // private functions.
    void ResetScalars();

  public:
    CompareGrids();

    void operator() (const float lhs, std::optional<float> rhs);

    int Report(std::ostream& output) const;

    void Reset();

    double Tanimoto() const;
};

CompareGrids::CompareGrids() {
  ResetScalars();
}

void
CompareGrids::ResetScalars() {
  _comparisons = 0;
  _in_range = 0;
  _rhs_larger = 0;
  _in_common = 0.0;
}

void
CompareGrids::operator() (const float lhs, std::optional<float> rhs) {
  ++_comparisons;
  if (! rhs) {
    return;
  }
  ++_in_range;
  _acc_lhs.extra(lhs);
  _acc_rhs.extra(*rhs);
  _acc_diffs.extra(*rhs - lhs);

  if (*rhs > lhs) {
    ++_rhs_larger;
    _in_common += lhs;
  } else {
    _in_common += *rhs;
  }
}

int
CompareGrids::Report(std::ostream& output) const {
  output << "made " << _comparisons << " comparisons, "
         << _in_range << " in range "
         << iwmisc::Fraction<float>(_in_range, _comparisons) << '\n';
  output << "Lhs btw " << _acc_lhs.minval() << " and "
         << _acc_lhs.maxval() << " mean "
         << static_cast<float>(_acc_lhs.average()) << '\n';
  output << "Rhs btw " << _acc_rhs.minval() << " and "
         << _acc_rhs.maxval() << " mean "
         << static_cast<float>(_acc_rhs.average()) << '\n';
  output << "diffs btw " << _acc_diffs.minval() << " and "
         << _acc_diffs.maxval() << " mean "
         << static_cast<float>(_acc_diffs.average()) << '\n';
  output << _rhs_larger << " rhs larger "
         << iwmisc::Fraction<float>(_rhs_larger, _comparisons) << '\n';

  const double tanimoto = Tanimoto();

  output << "Tanimoto " << static_cast<float>(tanimoto) << '\n';

  return 1;
}

double
CompareGrids::Tanimoto() const {
  return _in_common / (_acc_lhs.sum() + _acc_rhs.sum() - _in_common);
}

void
CompareGrids::Reset() {
  ResetScalars();

  _acc_lhs.reset();
  _acc_rhs.reset();
  _acc_diffs.reset();
}

class Options {
  private:
    int _verbose;

  public:
    Options();

    int Initialise(Command_Line& cl);

};

void
Usage(int rc) {
#ifdef GIT_HASH
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  cerr << "Compares values on two rectangular grids\n";
  cerr << "Takes two arguments, each must be 'x y z value' grids\n";
  cerr << "and must define a full grid - must be full rectangular box\n";
  cerr << " -t <tol>    tolerance for determining equivalent values.\n";
  cerr << "             default is epsilon, suggest something like 0.001\n";
  cerr << " -n          compare via nearest neighbour average. Default is to\n";
  cerr << "             compare via the closest point in the other grid\n";
  cerr << " -F <fname>  file containing file names of multiple grids.\n";
  cerr << "             The last argument is compared against each grid in <fname>\n";
  cerr << " -v          verbose output\n";

  ::exit(rc);
}

// Across an array of grids, return the extremeties.
std::array<float, 6>
FindExtremeties(const resizable_array_p<Grid>& grids) {
  std::array<float, 6> result;
  result[0] = std::numeric_limits<float>::max();
  result[1] = -std::numeric_limits<float>::max();
  result[2] = std::numeric_limits<float>::max();
  result[3] = -std::numeric_limits<float>::max();
  result[4] = std::numeric_limits<float>::max();
  result[5] = -std::numeric_limits<float>::max();
  for (const Grid* grid : grids) {
    std::array<float, 6> ext = grid->Extremities();
    result[0] = std::min(result[0], ext[0]);
    result[1] = std::max(result[1], ext[1]);
    result[2] = std::min(result[2], ext[2]);
    result[3] = std::max(result[3], ext[3]);
    result[4] = std::min(result[4], ext[4]);
    result[5] = std::max(result[5], ext[5]);
  }

  return result;
}

int
GetGrid(const char* fname, float tol,
        resizable_array_p<Grid>& destination) {
  std::unique_ptr<Grid> grid = std::make_unique<Grid>();
  grid->set_tolerance(tol);

  if (! grid->Build(fname)) {
    cerr << "GetGrid:cannot build grid from '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring string_fname(fname);
  if (string_fname.contains('/')) {
    const_IWSubstring basename;
    string_fname.iwbasename(basename);
    grid->set_name(basename);
  } else {
    grid->set_name(fname);
  }


  destination << grid.release();

  return destination.size();
}

// Read all but the last token on `cl` as grids into `destination`.
int
GridsFromCommandLine(const Command_Line& cl,
                     float tol,
                     resizable_array_p<Grid>& destination) {
  const int n = cl.size() - 1;
  for (int i = 0; i < n; ++i) {
    if (! GetGrid(cl[i], tol, destination)) {
      cerr << "GridsFromCommandLine:cannot build '" << cl[i] << "'\n";
      return 0;
    }
  }

  return destination.size();
}

int
GridsFromFile(iwstring_data_source& input, float tol,
        resizable_array_p<Grid>& destination) {

  IWString buffer;
  while (input.next_record(buffer)) {
    if (! GetGrid(buffer.null_terminated_chars(), tol, destination)) {
      cerr << "GridsFromFile:cannot read grid from '" << buffer << "'\n";
      return 0;
    }
  }

  return destination.size();
}

int
GridsFromFile(IWString& fname, float tol,
              resizable_array_p<Grid>& destination) {
  iwstring_data_source input;
  if (! input.open(fname.null_terminated_chars())) {
    cerr << "GridsFromFile:cannot open '" << fname << "'\n";
    return 0;
  }

  return GridsFromFile(input, tol, destination);
}

// If the -F option is present, read all tokens in `cl`
// as grids. Otherwise read just the last one.
int
GetQueryGrids(Command_Line& cl,
              float tol,
              resizable_array_p<Grid>& destination) {
  int istart;
  if (cl.option_present('F')) {
    istart = 0;
  } else {
    istart = cl.size() - 2;
  }

  int istop = cl.size();
  for (int i = istart; i < istop; ++i) {
    if (! GetGrid(cl[i], tol, destination)) {
      cerr << "GetQueryGrids:cannot read '" << cl[i] << "'\n";
      return 0;
    }
  }
  return destination.size();
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vt:nF:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  float tol = std::numeric_limits<float>::epsilon();

  if (cl.option_present('t')) {
    if (! cl.value('t', tol) || tol < 0.0f) {
      cerr << "The tolerance must be > 0\n";
      return 1;
    }
    if (verbose) {
      cerr << "Grid spacings determined to within " <<tol << '\n';
    }
  }

  // Grids against which we compare come either from the command line,
  // or from the -F option.
  resizable_array_p<Grid> comparison_grids;

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! GridsFromFile(fname, tol, comparison_grids)) {
      cerr << "Could not read multiple grids from '" << fname << "'\n";
      return 1;
    }
  } else {
    if (! GridsFromCommandLine(cl, tol, comparison_grids)) {
      cerr << "Could not read grids from command line\n";
      return 1;
    }
  }

  if (comparison_grids.empty()) {
    cerr << "No comparison grids\n";
    Usage(1);
  }

  if (verbose) {
    cerr << "Read " << comparison_grids.size() << " comparison grids\n";
    for (const Grid* grid : comparison_grids) {
      cerr << "Grid contains " << grid->npoints() << " points\n";
      std::array<float, 6> exti = grid->Extremities();
      cerr << " x " << exti[0] << ' ' << exti[1] << '\n';
      cerr << " y " << exti[2] << ' ' << exti[3] << '\n';
      cerr << " z " << exti[4] << ' ' << exti[5] << '\n';
    }
  }

  if (verbose) {
    cerr << " resolution\n";
    for (const Grid* grid : comparison_grids) {
      auto xyz = grid->DxDyDz();
      cerr << "Grid " << ' ' << grid->name();
      for (int j = 0; j < 3; ++j) {
        cerr << ' ' << xyz[j];
      }
      cerr << '\n';
    }
  }

  if (verbose) {
    std::array<float, 6> extremeties = FindExtremeties(comparison_grids);
    cerr << "Overal extremeties\n";
    cerr << " x " << extremeties[0] << ' ' << extremeties[1] << '\n';
    cerr << " y " << extremeties[2] << ' ' << extremeties[3] << '\n';
    cerr << " z " << extremeties[4] << ' ' << extremeties[5] << '\n';
  }

  resizable_array_p<Grid> query_grids;
  if (! GetQueryGrids(cl, tol, query_grids)) {
    cerr << "NO query grids\n";
    return 1;
  }

  if (verbose) {
    cerr << "Read " << query_grids.size() << " query grids\n";
    cerr << "Read " << comparison_grids.size() << " comparison grids\n";
  }

  IWString_and_File_Descriptor output(1);

  // Preserve original behaviour, two grid files on the command line, 
  // and no -F option.

  CompareGrids compare_grids;

  if (cl.size() == 2 && ! cl.option_present('F')) {
    if (cl.option_present('n')) {
      comparison_grids[0]->CompareViaNNs(*comparison_grids[1], compare_grids);
    } else {
      comparison_grids[0]->CompareViaClosest(*comparison_grids[1], compare_grids);
    }

    compare_grids.Report(std::cout);

    return 0;
  }

  // Branch according to how many grids are involved.
  if (cl.option_present('n')) {
    CompareViaNNs(comparison_grids, query_grids, compare_grids, output);
  } else {
    CompareViaClosest(comparison_grids, query_grids, compare_grids, output);
  }

#ifdef OLDDDDD
#endif

  return 0;
}

}  // namespace grid_overlap

int
main(int argc, char** argv) {
  int rc = grid_overlap::Main(argc, argv);

  return rc;
}
