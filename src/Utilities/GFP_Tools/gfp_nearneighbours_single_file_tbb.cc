/*
  Finds near neighbours within a single fingerprint file
*/

#include <omp.h>
#include <pthread.h>
#include <stdint.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>

#include "tbb/scalable_allocator.h"
#include "tbb/task_group.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/sparse_collection.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int neighbours_to_find = 0;

static int write_neighbours_as_index_numbers = 0;

static Accumulator<float> distance_stats;

static Fraction_as_String fraction_as_string;

/*
  Sometimes it is useful to know the statistics of the nearest neighbour of each molecule
*/

static Accumulator<similarity_type_t> nearest_neighbour_distance_stats;

static IWHistogram histogram_nearnest_neighbour_distances;

static int create_histogram = 0;

static int write_minimal_histogram = 0;

static extending_resizable_array<int> neighbour_count;

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

static IWString number_neighbours_tag;

static float float_max_32_but_unsigned_int =
    static_cast<float>(std::numeric_limits<uint32_t>::max() - 130);

/*
  If we are writing neighbours as index numbers
*/

static IWString neighbour_tag("NBR<");

/*
  May 99. For each input TDT, I need to know the average distance
  of the neighbours within the pool
*/

static IWString tag_for_average_distance;

static const_IWSubstring distance_tag("DIST<");

/*
  We ignore distances longer than DISTANCE_THRESHOLD
*/

static similarity_type_t upper_distance_threshold = std::numeric_limits<float>::max();
static similarity_type_t lower_distance_threshold = 0.0f;

static int allow_arbitrary_distances = 0;

/*
  When we have thresholds, we may choose to not write molecules with no neighbours
*/

static int write_molecules_with_no_neighbours = 1;

// static int molecules_with_no_neighbours = 0;
static std::atomic<int> molecules_with_no_neighbours(0);

static int write_smiles = 1;

static int include_neighbour_count_with_target_identifier = 0;

static int nworkers = 2;

/*
  We can encode an integer ID and a floating point distance into a single 64 bit int
  We need methods to compare these
*/

class Encoded_Dist_ID_Sorter
{
 private:
 public:
  int
  operator()(const uint64_t i1, const uint64_t i2) const;
};

static uint64_t
encode_distance_id(const float d, const int id)
{
  uint64_t rc = d * float_max_32_but_unsigned_int;
  // cerr << "With distance encoded " << std::hex << rc << ' ' << std::dec << rc << endl;

  rc = rc << 32;
  // cerr << "Shifted " << std::hex << rc << ' ' << std::dec << rc << endl;

  uint64_t idt = id;

  return rc | idt;
}

static float
decode_distance(const uint64_t did)
{
  uint32_t d = (did >> 32);

  return static_cast<float>(d) / float_max_32_but_unsigned_int;
}

static uint32_t
decode_id(const uint64_t did)
{
  const uint32_t* i = reinterpret_cast<const uint32_t*>(&did);
  return i[0];
}

int
Encoded_Dist_ID_Sorter::operator()(const uint64_t i1, const uint64_t i2) const
{
  if (i1 < i2) {
    return -1;
  }
  if (i1 > i2) {
    return 1;
  }

  return 0;
}

template <typename F>
int
write_neighbour_list(const F& target, const F* pool, IWString_and_File_Descriptor& output)
{
  const resizable_array<uint64_t>& nbrs = target.nbrs();

  const int n = nbrs.number_elements();

  if (0 == n && !write_molecules_with_no_neighbours) {
    return 1;
  }

  target.write_smiles_id(include_neighbour_count_with_target_identifier, output);

  for (int i = 0; i < n; ++i) {
    const auto j = nbrs[i];

    const uint32_t id = decode_id(j);

    if (write_neighbours_as_index_numbers) {
      if (write_smiles) {
        output << smiles_tag << pool[id].smiles() << ">\n";
      }
      output << neighbour_tag << id << ">\n";
    } else {
      pool[id].write_smiles_id(0, output);
    }

    const float d = decode_distance(j);
    //  cerr << "From " << j << " get id " << id << " and distance " << d << endl;

    fraction_as_string.append_number(output, d);

    if (verbose) {
      distance_stats.extra(d);
    }

    if (0 == i && (create_histogram || verbose)) {
      nearest_neighbour_distance_stats.extra(d);
      if (create_histogram) {
        histogram_nearnest_neighbour_distances.extra(d);
      }
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  output << "|\n";

  return 1;
}

/*
  We work in one of two modes.
  Fixed number of neighbours
  ALl molecules within a distance
*/

class FP_Fixed_Distance : public IW_General_Fingerprint
{
 private:
  resizable_array<uint64_t> _dist_id;

  IWString _smiles;

 public:
  FP_Fixed_Distance();
  ~FP_Fixed_Distance();

  void
  set_smiles(const IWString& s)
  {
    _smiles = s;
  }

  const IWString&
  smiles() const
  {
    return _smiles;
  }

  int
  write_smiles_id(const int write_nbr_count, IWString_and_File_Descriptor&) const;

  float
  tanimoto(FP_Fixed_Distance& rhs);

  void
  extra(const float d, const int id);

  void
  sort_nbrs();

  const resizable_array<uint64_t>&
  nbrs() const
  {
    return _dist_id;
  };

  int
  do_write(const FP_Fixed_Distance* pool, IWString_and_File_Descriptor& output) const
  {
    return write_neighbour_list(*this, pool, output);
  };
};

class FP_Number_Nbrs : public IW_General_Fingerprint
{
 private:
  //  std::priority_queue<int64_t, std::vector<int64_t>, std::greater<int64_t>> _nbrs;
  std::priority_queue<uint64_t> _nbrs;

  resizable_array<uint64_t> _sorted;

  float _dmax;

  IWString _smiles;

 public:
  FP_Number_Nbrs();
  ~FP_Number_Nbrs();

  void
  set_smiles(const IWString& s)
  {
    _smiles = s;
  }

  const IWString&
  smiles() const
  {
    return _smiles;
  }

  int
  write_smiles_id(const int write_nbr_count, IWString_and_File_Descriptor&) const;

  float
  tanimoto(FP_Number_Nbrs& rhs);

  void
  extra(const float d, const int id);
  void
  _extra(const float d, const int id);

  void
  sort_nbrs();

  const resizable_array<uint64_t>&
  nbrs() const
  {
    return _sorted;
  };

  int
  do_write(const FP_Number_Nbrs* pool, IWString_and_File_Descriptor& output) const
  {
    return write_neighbour_list(*this, pool, output);
  };
};

FP_Fixed_Distance::FP_Fixed_Distance()
{
  return;
}

FP_Fixed_Distance::~FP_Fixed_Distance()
{
}

void
FP_Fixed_Distance::extra(const float d, const int id)
{
  // cerr << "cmp " << d << " with " << upper_distance_threshold << " size " <<
  // _dist_id.size() << endl;
  if (d > upper_distance_threshold || d < lower_distance_threshold) {
    return;
  }

  // if (4664 == id && fabs(d - 0.2416) < 0.0001)
  //   cerr << "Got " << id << " at dist " << d << endl;

  uint64_t e = encode_distance_id(d, id);

  _dist_id.add(e);

  return;
}

float
FP_Fixed_Distance::tanimoto(FP_Fixed_Distance& rhs)
{
  IW_General_Fingerprint& gfp_lhs = *this;
  IW_General_Fingerprint& gfp_rhs = rhs;

  return gfp_lhs.tanimoto(gfp_rhs);
}

void
FP_Fixed_Distance::sort_nbrs()
{
  Encoded_Dist_ID_Sorter edids;

  _dist_id.iwqsort(edids);

  return;
}

FP_Number_Nbrs::FP_Number_Nbrs()
{
  _sorted.resize(neighbours_to_find);
  _dmax = 1.0f;
}

FP_Number_Nbrs::~FP_Number_Nbrs()
{
}

#ifdef VARIANT_NEEDED_FOR_TESTING
void
FP_Number_Nbrs::extra(const float d, const int id)
{
  _extra(d, id);

#define CHECK_SORTED_FP
#ifdef CHECK_SORTED_FP
  uint64_t e = encode_distance_id(d, id);
  for (int i = 1; i < _sorted.number_elements(); ++i) {
    if (_sorted[i] <= _sorted[i - 1]) {
      cerr << "List out of order, i = " << i << " prev " << _sorted[i - 1] << " now "
           << _sorted[i] << " e = " << e << endl;
    }
  }
  const float x = decode_distance(_sorted.last_item());
  if (x != _dmax) {
    cerr << "FP_Number_Nbrs::extra:dmax mismatch, found " << x << " dmax " << _dmax
         << " n = " << _sorted.size() << endl;
  }
#endif

  return;
}
#endif

void
FP_Number_Nbrs::extra(const float d, const int id)
{
  // cerr << "FP_Number_Nbrs:extra distance " << d << " have " <<
  // _sorted.number_elements() << " max " << _dmax << " id " << id << endl;

  if (_sorted.number_elements() == neighbours_to_find &&
      d >= _dmax) {  // most common case
    return;
  }

  if (d > upper_distance_threshold || d < lower_distance_threshold) {
    return;
  }

  const uint64_t e = encode_distance_id(d, id);

  // cerr << "encoded " << e << ' ' << (1.0f - d) << endl;

  if (0 == _sorted.number_elements()) {
    _sorted.add(e);
    _dmax = d;
    return;
  }

  if (_sorted.number_elements() == neighbours_to_find)  // make room
  {
    if (1 == neighbours_to_find)  // must be treated specially
    {
      _sorted[0] = e;
      _dmax = d;
      return;
    }

    _sorted.pop();
    _dmax = decode_distance(_sorted.last_item());
  }

  if (e < _sorted.first()) {
    _sorted.insert_at_beginning(e);
    return;
  }

  if (e > _sorted.last_item()) {
    _sorted.add(e);
    _dmax = d;
    return;
  }

  int left = 0;
  int right = _sorted.number_elements() - 1;
  int middle = (left + right) / 2;

  const auto s = _sorted.rawdata();

  // cerr << "INserting " << e << " into list of size " << _sorted.size() << endl;

  while (middle > left) {
    const auto m = s[middle];

    if (e < m) {
      right = middle;
    } else if (e > m) {
      left = middle;
    }

    middle = (left + right) / 2;

    //  cerr << "Update " << left << " " << middle << " " << right << endl;
  }

  _sorted.insert_before(middle + 1, e);

  return;
}

float
FP_Number_Nbrs::tanimoto(FP_Number_Nbrs& rhs)
{
  IW_General_Fingerprint& gfp_lhs = *this;
  IW_General_Fingerprint& gfp_rhs = rhs;

  return gfp_lhs.tanimoto(gfp_rhs);
}

void
FP_Number_Nbrs::sort_nbrs()
{
  return;
  assert(0 == _sorted.size());

  const int s = _nbrs.size();

  _sorted.extend(s);

  for (int i = 0; i < s; ++i) {
    _sorted[s - i - 1] = _nbrs.top();
    _nbrs.pop();
  }

  return;
}

int
FP_Number_Nbrs::write_smiles_id(const int write_nbr_count,
                                IWString_and_File_Descriptor& output) const
{
  output << smiles_tag << _smiles << ">\n";
  output << identifier_tag << id();
  if (write_nbr_count) {
    output << ' ' << _sorted.number_elements();  // make sure sorted!
  }
  output << ">\n";

  return 1;
}

int
FP_Fixed_Distance::write_smiles_id(const int write_nbr_count,
                                   IWString_and_File_Descriptor& output) const
{
  output << smiles_tag << _smiles << ">\n";
  output << identifier_tag << id();
  if (write_nbr_count) {
    output << ' ' << _dist_id.number_elements();
  }
  output << ">\n";

  return 1;
}

template <typename F>
void
gfp_nearneighbours(F* pool, const int istart, const int istop, const int jstart,
                   const int jstop)
{
  // cerr << istart << ' ' << istop << ' ' << jstart << ' ' << jstop << endl;
  if (property_based_windows_present()) {
    for (int i = istart; i < istop; ++i) {
      auto& pi = pool[i];

      for (int j = jstart; j < jstop; ++j) {
        auto& pj = pool[j];
        if (!can_be_compared(pi, pj)) {
          continue;
        }

        const auto d = 1.0f - pi.tanimoto(pj);
        //      cerr << " i = " << i << " j = " << j << " dist " << d << endl;

        pi.extra(d, j);
        pj.extra(d, i);
      }
    }
  } else {
    for (int i = istart; i < istop; ++i) {
      auto& pi = pool[i];

      for (int j = jstart; j < jstop; ++j) {
        auto& pj = pool[j];

        const auto d = 1.0f - pi.tanimoto(pj);
        //      cerr << " i = " << i << " j = " << j << " dist " << d << endl;

        pi.extra(d, j);
        pj.extra(d, i);
      }
    }
  }

  return;
}

template <typename F>
void
gfp_nearneighbours_diagonal(F* pool, const int istart, const int istop)
{
  if (property_based_windows_present()) {
    for (int i = istart; i < istop; ++i) {
      auto& pi = pool[i];

      for (int j = i + 1; j < istop; ++j) {
        auto& pj = pool[j];

        if (!can_be_compared(pi, pj)) {
          continue;
        }

        const auto d = 1.0f - pi.tanimoto(pj);
        //      cerr << " i = " << i << " j = " << j << " dist " << d << endl;

        pi.extra(d, j);
        pj.extra(d, i);
      }
    }
  } else {
    for (int i = istart; i < istop; ++i) {
      auto& pi = pool[i];

      for (int j = i + 1; j < istop; ++j) {
        auto& pj = pool[j];

        const auto d = 1.0f - pi.tanimoto(pj);
        //      cerr << " i = " << i << " j = " << j << " dist " << d << endl;

        pi.extra(d, j);
        pj.extra(d, i);
      }
    }
  }

  return;
}

/*
  Not really parallel at all
*/

template <typename F>
void
gfp_nearneighbours_parallel1(F* pool, const int pool_size)
{
  gfp_nearneighbours_diagonal(pool, 0, pool_size);

  return;
}

template <typename F>
void
gfp_nearneighbours_parallel2(F* pool, const int pool_size)
{
  const int half = pool_size / 2;

  tbb::task_group g1;
  g1.run([&] { gfp_nearneighbours_diagonal(pool, 0, half); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, half, pool_size); });
  g1.wait();
  if (verbose) {
    cerr << "End 1\n";
  }

  gfp_nearneighbours(pool, 0, half, half, pool_size);
  if (verbose) {
    cerr << "End 2\n";
  }

  return;
}

template <typename F>
void
gfp_nearneighbours_parallel3(F* pool, const int pool_size)
{
  const int third = pool_size / 3;
  const int two_thirds = (2 * pool_size) / 3;

  tbb::task_group g1;
  g1.run([&] { gfp_nearneighbours_diagonal(pool, 0, third); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, third, two_thirds); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, two_thirds, pool_size); });
  g1.wait();

  gfp_nearneighbours(pool, 0, third, third, two_thirds);

  gfp_nearneighbours(pool, third, two_thirds, two_thirds, pool_size);

  gfp_nearneighbours(pool, 0, third, two_thirds, pool_size);

  return;
}

/*
 Greco Latin square design
 http://www.unh.edu/halelab/BIOL933/supp_mats/Useful_LSs.pdf

  a b c d e f g h
  b a d c f e h g
  c d a b g h e f
  d c b a h g f e
  e f g h a b c d
  f e h g b a d c
  g h e f c d a b
  h g f e d c b a

  Idea from Steve Ruberg
*/

template <typename F>
void
gfp_nearneighbours_parallel8_latin_square(F* pool, const int pool_size)
{
  int p[9];
  for (int i = 0; i < 8; ++i) {
    p[i] = (i * pool_size) / 8;
  }

  p[8] = pool_size;

  tbb::task_group g1;
  g1.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[7], p[8]); });
  g1.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[6], p[7]); });
  g1.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[5], p[6]); });
  gfp_nearneighbours(pool, p[3], p[4], p[4], p[5]);
  g1.wait();
  if (verbose) {
    cerr << "End 1\n";
  }

  tbb::task_group g2;
  g2.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[6], p[7]); });
  g2.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[7], p[8]); });
  g2.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[4], p[5]); });
  gfp_nearneighbours(pool, p[3], p[4], p[5], p[6]);
  g2.wait();
  if (verbose) {
    cerr << "End 2\n";
  }

  tbb::task_group g3;
  g3.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[5], p[6]); });
  g3.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[4], p[5]); });
  g3.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[7], p[8]); });
  gfp_nearneighbours(pool, p[3], p[4], p[6], p[7]);
  g3.wait();
  if (verbose) {
    cerr << "End 3\n";
  }

  tbb::task_group g4;  // D
  g4.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[4], p[5]); });
  g4.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[5], p[6]); });
  g4.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[6], p[7]); });
  gfp_nearneighbours(pool, p[3], p[4], p[7], p[8]);
  g4.wait();
  if (verbose) {
    cerr << "End 4\n";
  }

  tbb::task_group g5;  // E
  g5.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[3], p[4]); });
  g5.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[2], p[3]); });
  g5.run([&] { gfp_nearneighbours(pool, p[4], p[5], p[7], p[8]); });
  gfp_nearneighbours(pool, p[5], p[6], p[6], p[7]);
  g5.wait();
  if (verbose) {
    cerr << "End 5\n";
  }

  tbb::task_group g6;  // F
  g6.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[2], p[3]); });
  g6.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[3], p[4]); });
  g6.run([&] { gfp_nearneighbours(pool, p[4], p[5], p[6], p[7]); });
  gfp_nearneighbours(pool, p[5], p[6], p[7], p[8]);
  g6.wait();
  if (verbose) {
    cerr << "End 6\n";
  }

  tbb::task_group g7;  // G
  g7.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[1], p[2]); });
  g7.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[3], p[4]); });
  g7.run([&] { gfp_nearneighbours(pool, p[4], p[5], p[5], p[6]); });
  gfp_nearneighbours(pool, p[6], p[7], p[7], p[8]);
  g7.wait();
  if (verbose) {
    cerr << "End 7\n";
  }

  tbb::task_group g8;
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[0], p[1]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[1], p[2]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[2], p[3]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[3], p[4]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[4], p[5]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[5], p[6]); });
  g8.run([&] { gfp_nearneighbours_diagonal(pool, p[6], p[7]); });
  gfp_nearneighbours_diagonal(pool, p[7], p[8]);
  g8.wait();
  if (verbose) {
    cerr << "End 8\n";
  }

  return;
}

/*
  Latin square size 6
  1 2 3 4 5 6
  2 4 1 6 3 5
  3 1 5 2 6 4
  4 6 2 5 1 3
  5 3 6 1 4 2
  6 5 4 3 2 1

  From
  http://www.stat.purdue.edu/docs/research/tech-reports/1980/tr80-26.pdf
*/

template <typename F>
void
gfp_nearneighbours_parallel6(F* pool, const int pool_size)
{
  int p[7];
  for (int i = 0; i < 6; ++i) {
    p[i] = (i * pool_size) / 6;
  }

  p[6] = pool_size;

  tbb::task_group g1;
  g1.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[5], p[6]); });
  g1.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[3], p[4]); });
  g1.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[4], p[5]); });
  g1.wait();
  if (verbose) {
    cerr << "End 1\n";
  }

  tbb::task_group g2;
  g2.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[4], p[5]); });
  g2.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[5], p[6]); });
  g2.run([&] { gfp_nearneighbours_diagonal(pool, p[2], p[3]); });
  g2.run([&] { gfp_nearneighbours_diagonal(pool, p[3], p[4]); });
  g2.wait();
  if (verbose) {
    cerr << "End 2\n";
  }

  tbb::task_group g3;
  g3.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[3], p[4]); });
  g3.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[5], p[6]); });
  g3.run([&] { gfp_nearneighbours_diagonal(pool, p[1], p[2]); });
  g3.run([&] { gfp_nearneighbours_diagonal(pool, p[4], p[5]); });
  g3.wait();
  if (verbose) {
    cerr << "End 3\n";
  }

  tbb::task_group g4;
  g4.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[2], p[3]); });
  g4.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[4], p[5]); });
  g4.run([&] { gfp_nearneighbours(pool, p[3], p[4], p[5], p[6]); });
  g4.wait();
  if (verbose) {
    cerr << "End 4\n";
  }

  tbb::task_group g5;
  g5.run([&] { gfp_nearneighbours(pool, p[0], p[1], p[1], p[2]); });
  g5.run([&] { gfp_nearneighbours(pool, p[2], p[3], p[3], p[4]); });
  g5.run([&] { gfp_nearneighbours(pool, p[4], p[5], p[5], p[6]); });
  g5.wait();
  if (verbose) {
    cerr << "End 5\n";
  }

  tbb::task_group g6;
  g6.run([&] { gfp_nearneighbours_diagonal(pool, p[0], p[1]); });
  g6.run([&] { gfp_nearneighbours(pool, p[1], p[2], p[2], p[3]); });
  g6.run([&] { gfp_nearneighbours(pool, p[3], p[4], p[4], p[5]); });
  g6.run([&] { gfp_nearneighbours_diagonal(pool, p[5], p[6]); });
  g6.wait();
  if (verbose) {
    cerr << "End 6\n";
  }

  return;
}

template <typename F>
void
gfp_nearneighbours_parallel4(F* pool, const int pool_size)
{
  int p[5];
  for (int i = 0; i < 4; ++i) {
    p[i] = (i * pool_size) / 4;
  }

  p[4] = pool_size;

  const int quarter = pool_size / 4;
  const int half = pool_size / 2;
  const int three_quarters = (3 * pool_size) / 4;

  tbb::task_group g1;
  g1.run([&] { gfp_nearneighbours_diagonal(pool, 0, p[1]); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, p[1], p[2]); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, p[2], p[3]); });
  g1.run([&] { gfp_nearneighbours_diagonal(pool, p[3], p[4]); });
  g1.wait();
  if (verbose) {
    cerr << "End 1\n";
  }

  tbb::task_group g2;
  g2.run([&] { gfp_nearneighbours(pool, 0, quarter, quarter, half); });
  g2.run(
      [&] { gfp_nearneighbours(pool, half, three_quarters, three_quarters, pool_size); });
  g2.wait();
  if (verbose) {
    cerr << "End 2\n";
  }

  tbb::task_group g3;
  g3.run([&] { gfp_nearneighbours(pool, 0, quarter, three_quarters, pool_size); });
  g3.run([&] { gfp_nearneighbours(pool, quarter, half, half, three_quarters); });
  // g3.run([&] {gfp_nearneighbours(pool, half, three_quarters, quarter, half);});
  // g3.run([&] {gfp_nearneighbours(pool, three_quarters, pool_size, 0, quarter);});
  g3.wait();
  if (verbose) {
    cerr << "End 3\n";
  }

  tbb::task_group g4;
  g4.run([&] { gfp_nearneighbours(pool, 0, quarter, half, three_quarters); });
  g4.run([&] { gfp_nearneighbours(pool, quarter, half, three_quarters, pool_size); });
  g4.wait();
  if (verbose) {
    cerr << "End 4\n";
  }

  return;
}

template <typename F>
void
gfp_nearneighbours(F* pool, const int pool_size)
{
  for (int i = 0; i < pool_size; ++i) {
    auto& pi = pool[i];

    for (int j = i + 1; j < pool_size; ++j) {
      auto& pj = pool[j];

      const auto d = 1.0f - pi.tanimoto(pj);
      //    cerr << " i = " << i << " j = " << j << " dist " << d << endl;

      pi.extra(d, j);
      pj.extra(d, i);
    }
  }

  return;
}

template <typename F>
int
gfp_nearneighbours(F* pool, const int pool_size, IWString_and_File_Descriptor& output)
{
  // cerr << "nworkers " << nworkers << endl;
  if (2 == nworkers) {
    gfp_nearneighbours_parallel2(pool, pool_size);
  } else if (3 == nworkers) {
    gfp_nearneighbours_parallel3(pool, pool_size);
  } else if (4 == nworkers) {
    gfp_nearneighbours_parallel4(pool, pool_size);
  } else if (6 == nworkers) {
    gfp_nearneighbours_parallel6(pool, pool_size);
  } else if (8 == nworkers) {
    gfp_nearneighbours_parallel8_latin_square(pool, pool_size);
  } else if (1 == nworkers) {
    gfp_nearneighbours_parallel1(pool, pool_size);
  } else {
    cerr << "Do not know how to process " << nworkers << " threads, see Ian\n";
    return 0;
  }

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < pool_size; ++i) {
      pool[i].sort_nbrs();
    }
  }  // pragma omp parallel

  for (int i = 0; i < pool_size; ++i) {
    pool[i].do_write(pool, output);
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

template <typename F>
int
gfp_nearneighbours(const char* fname, int& pool_size,
                   IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "gfp_nearneighbours:cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }
  }

  F* pool = new F[pool_size];
  std::unique_ptr<F[]> free_pool(pool);

  if (NULL == pool) {
    cerr << "gfp_nearneighbours:cannot allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  IW_TDT tdt;
  int fatal;
  for (int ndx = 0; tdt.next(input); ++ndx) {
    fatal = 0;
    if (!pool[ndx].construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot build pool, ndx " << ndx << endl;
      return 0;
    }

    IWString smiles;
    if (!tdt.dataitem_value(smiles_tag, smiles)) {
      cerr << "Cannot extract smiles, ndx " << ndx << endl;
      return 0;
    }

    pool[ndx].set_smiles(smiles);
  }

  return gfp_nearneighbours(pool, pool_size, output);
}

/*
  When doing near neighbour determinations within a single set of
  molecules, we gain great efficiencies by writing the index of the
  item rather than its name.  This helps programmes that read the nn
  file
*/

class IW_GFP_D_ID : public IW_GFP_D
{
 private:
  int _ndx;
  IWString _smiles;

 public:
  IW_GFP_D_ID();

  void
  set_index(int n)
  {
    _ndx = n;
  }

  int
  index_in_pool() const
  {
    return _ndx;
  }

  int
  write_smiles_and_id(IWString&) const;

  void
  set_smiles(const const_IWSubstring& s)
  {
    _smiles = s;
  }

  const IWString&
  smiles() const
  {
    return _smiles;
  }
};

IW_GFP_D_ID::IW_GFP_D_ID()
{
  _ndx = -1;
}

int
IW_GFP_D_ID::write_smiles_and_id(IWString& output) const
{
  if (_smiles.length()) {
    output << smiles_tag << _smiles << ">\n";
  }

  output << identifier_tag << IW_General_Fingerprint::id() << ">\n";

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static IW_GFP_D_ID* pool = nullptr;

static int pool_size = 0;

static int
build_pool(iwstring_data_source& input)
{
  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;

    int fatal;
    if (!pool[items_in_pool].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    pool[items_in_pool].set_index(items_in_pool);

    if (write_smiles) {
      const_IWSubstring smi;
      if (!tdt.dataitem_value(smiles_tag, smi)) {
        cerr << "Cannot extract smiles\n";
        cerr << tdt;
        return 0;
      }

      pool[items_in_pool].set_smiles(smi);
    }

    items_in_pool++;

    if (items_in_pool == pool_size) {
      if (verbose) {
        cerr << "Pool is full, max " << pool_size << endl;
      }
      break;
    }
  }

  pool_size = items_in_pool;

  if (verbose) {
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size
         << " fingerprints\n";
  }

  return 1;
}

static int
build_pool(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    pool = new IW_GFP_D_ID[pool_size];
    if (NULL == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose) {
      cerr << "Pool automatically sized to " << pool_size << endl;
    }
  }

  return build_pool(input);
}

#ifdef NOT_IN_USE
static int
write_average_neighbour_distance(IW_GFP_D_ID** neighbours, int number_neighbours,
                                 IWString& output)
{
  Accumulator<similarity_type_t> acc;

  for (int i = 0; i < number_neighbours; i++) {
    IW_GFP_D_ID* n = neighbours[i];
    assert(NULL != neighbours[i]);

    acc.extra(n->distance());
  }

  output << tag_for_average_distance << acc.average_if_available_minval_if_not() << ">\n";

  return output.good();
}
#endif

/*
  Common function for doing the distance computation
*/

similarity_type_t
compute_the_distance(IW_GFP_D_ID& fp1, IW_GFP_D_ID& fp2)
{
  return static_cast<similarity_type_t>(1.0) - fp1.tanimoto(fp2);
}

static int
three_column_output_all_pairs(IWString_and_File_Descriptor& output)
{
  for (int i = 0; i < pool_size; i++) {
    IW_GFP_D_ID& fpi = pool[i];

    for (int j = i + 1; j < pool_size; j++) {
      similarity_type_t d = compute_the_distance(fpi, pool[j]);

      if (d > upper_distance_threshold || d < lower_distance_threshold) {
        continue;
      }

      output << fpi.id() << ' ' << pool[j].id() << ' ' << d << '\n';
      output.write_if_buffer_holds_more_than(32768);
    }
  }

  output.flush();

  return 1;
}

#ifdef TIS_IS_NOT_USED
static int
nearneighbours(IW_GFP_D_ID& fp, IW_GFP_D_ID** neighbours,
               int neighbours_to_find_this_fingerprint)
{
  int neighbours_found = 0;
  assert(neighbours_to_find_this_fingerprint > 0);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < pool_size; ++i) {
      pool[i].set_distance(2.0f);
      neighbours[i] = pool + i;
    }
  }  // pragma omp parallel

  neighbours_found = 0;

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < pool_size; i++) {
      if (!can_be_compared(fp, pool[i])) {
        continue;
      }

      similarity_type_t t;
      t = compute_the_distance(fp, pool[i]);

// #define DEBUG_NN
#ifdef DEBUG_NN
      cerr << "Distance between '" << fp.id() << " and pool " << i << " '" << pool[i].id()
           << "' is " << t << endl;
#endif

      pool[i].set_distance(t);
    }
  }  // pragma omp parallel

  std::partial_sort(neighbours, neighbours + pool_size,
                    neighbours + neighbours_to_find_this_fingerprint,
                    [](const IW_GFP_D_ID* fp1, const IW_GFP_D_ID* fp2) {
                      return fp1->distance() < fp2->distance();
                    });

  if (std::numeric_limits<float>::max() == upper_distance_threshold) {
    return neighbours_to_find_this_fingerprint;
  }

  for (int i = 0; i < neighbours_to_find_this_fingerprint; ++i) {
    const auto d = neighbours[i]->distance();

    if (d > upper_distance_threshold || d < lower_distance_threshold) {
      return i - 1;
    }
  }

  return neighbours_to_find_this_fingerprint;
}
#endif

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

  cerr << "Finds near neighbours of a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -T <dis>         discard distances longer than <dis>\n";
  cerr << " -z               don't write molecules with no neighbours\n";
  cerr << " -I <tag>         specify identifier dataitem (default '" << identifier_tag << ")\n";
  cerr << " -A <TAG>         write average neighbour distance to <TAG>\n";
//cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -o               cross referencing a single file. Write neighbours as index numbers\n";
  cerr << " -H <fname>       write histogram of closest distances to <fname>\n";
  cerr << " -b               write minimal histogram data - two columns\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -N <tag>         write number neighbours as <tag>\n";
  cerr << " -p               write all pair-wise distances in 3 column form\n";
  cerr << " -j <precision>   output precision for distances\n";
  cerr << " -x               exclude smiles from the output\n";
  cerr << " -y               allow arbitrary distances\n";
  cerr << " -u               include number of neighbours with target identifier\n";
  // cerr << " -C <nthreads>    number of TBB threads to use\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_write_histogram(std::ostream& os)
{
  if (write_minimal_histogram) {
    return histogram_nearnest_neighbour_distances.write_terse(os);
  } else {
    return histogram_nearnest_neighbour_distances.write(os);
  }
}

static int
nearneighbours(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:n:I:T:t:P:F:W:Q:A:zoxH:N:bypj:C:u");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('I')) {
    (void)cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose) {
      cerr << "Identifiers tagged as '" << identifier_tag << "'\n";
    }
  }

  if (!iw_little_endian()) {
    cerr << "Sorry, this programme only works on little endian machines, contact Ian\n";
    return 1;
  }

// #define DEBUG_ENCODING_STUFF
#ifdef DEBUG_ENCODING_STUFF
  for (int i = 0; i <= 100; ++i) {
    const uint64_t e1 = encode_distance_id(0.01 * i, 0);
    const uint64_t e2 = encode_distance_id(0.01 * i, 100);

    cerr << (0.01 * i) << ' ' << e1 << ' ' << e2 << endl;

    const float d1 = decode_distance(e1);
    const float d2 = decode_distance(e2);

    const int id1 = decode_id(e1);
    const int id2 = decode_id(e2);

    cerr << (0.01 * i) << " id1 " << id1 << " d1 " << d1 << " id2 " << id2 << " d2 " << d2
         << endl;
  }
#endif

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  Set_of_Sparse_Fingerprint_Collection_Profile sfcp;

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  } else if (!initialise_fingerprints(cl[0], verbose)) {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  // std::unique_ptr<tbb::task_scheduler_init[]> free_init(init);   seems to expose
  // problems

  if (cl.option_present('u')) {
    include_neighbour_count_with_target_identifier = 1;

    if (verbose) {
      cerr << "Will include neighbour count with target identifier\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new IW_GFP_D_ID[pool_size];
    if (NULL == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose) {
      cerr << "system sized to " << pool_size << endl;
    }
  }

  if (!build_pool(cl[0])) {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 5;
  }

  if (0 == pool_size) {
    cerr << "No fingerprints\n";
    return 5;
  }

  std::unique_ptr<IW_GFP_D_ID[]> free_pool(pool);

  if (cl.option_present('y')) {
    allow_arbitrary_distances = 1;

    if (verbose) {
      cerr << "Distances not constrained to [0-1]\n";
    }
  }

  std::ofstream stream_for_nearest_neighbour_histogram;

  if (cl.option_present('H')) {
    if (allow_arbitrary_distances) {
      cerr << "Sorry the histogram option is inconsistent with the -y option\n";
      usage(4);
    }

    create_histogram = 1;

    if (cl.option_present('b')) {
      write_minimal_histogram = 1;
      if (verbose) {
        cerr << "Only a minimal histogram will be written\n";
      }
    }

    const char* h = cl.option_value('H');
    stream_for_nearest_neighbour_histogram.open(h, std::ios::out);
    if (!stream_for_nearest_neighbour_histogram.good()) {
      cerr << "Sorry, cannot open histogram stream '" << h << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Histogram data written to '" << h << "'\n";
    }
  }

  if (cl.option_present('N')) {
    cl.value('N', number_neighbours_tag);

    if (verbose) {
      cerr << "The number of neighbours written as '" << number_neighbours_tag << "'\n";
    }

    if (!number_neighbours_tag.ends_with('<')) {
      number_neighbours_tag << '<';
    }
  }

  if (cl.option_present('K')) {
    for (int i = 0; i < pool_size; i++) {
      sfcp.build_profile(pool[i]);
    }

    sfcp.finished_profiling(verbose);

    if (verbose) {
      sfcp.report(cerr);
    }

    for (int i = 0; i < pool_size; i++) {
      pool[i].convert_to_non_sparse_forms(sfcp);
    }
  }

  // process the -o option here, after the -e but before the -E option

  if (cl.option_present('o')) {
    write_neighbours_as_index_numbers = 1;
    write_smiles = 0;

    if (verbose) {
      cerr << "Neighbours written as index numbers\n";
    }
  }

  if (cl.option_present('x')) {
    write_smiles = 0;

    if (verbose) {
      cerr << "Smiles excluded from output\n";
    }
  }

  if (cl.option_present('A')) {
    tag_for_average_distance = cl.string_value('A');

    if (verbose) {
      cerr << "The average neighbour distance will be written to '"
           << tag_for_average_distance << "'\n";
    }

    if (!tag_for_average_distance.ends_with('<')) {
      tag_for_average_distance += '<';
    }
  }

  if (cl.option_present('T')) {
    if (!cl.value('T', upper_distance_threshold) || upper_distance_threshold < 0.0 ||
        upper_distance_threshold > 1.0) {
      cerr << "The -T option must be followed by a valid distance\n";
      usage(12);
    }

    if (verbose) {
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', lower_distance_threshold) || lower_distance_threshold < 0.0 ||
        lower_distance_threshold > 1.0) {
      cerr << "The -t option must be followed by a valid distance\n";
      usage(12);
    }

    if (verbose) {
      cerr << "Lower distance threshold set to " << lower_distance_threshold << endl;
    }
  }

  if (lower_distance_threshold >= upper_distance_threshold) {
    cerr << "Inconsistent lower " << lower_distance_threshold << " and upper "
         << upper_distance_threshold << " distance thresholds\n";
    return 1;
  }

  if (cl.option_present('z') && !cl.option_present('t') && !cl.option_present('T') &&
      !cl.option_present('W')) {
    cerr << "The don't write molecules with no neighbours option (-z) only makes sense "
            "with thresholds and/or windows\n";
    usage(13);
  }

  if (cl.option_present('z')) {
    write_molecules_with_no_neighbours = 0;
    if (verbose) {
      cerr << "Will not write molecules with no neighbours\n";
    }
  }

  if (cl.option_present('n')) {
    const_IWSubstring nvalue = cl.string_value('n');
    if ("all" == nvalue) {
      neighbours_to_find = pool_size;

      if (verbose) {
        cerr << "May get as many as " << pool_size << " neighbours\n";
      }
    } else {
      if (!nvalue.numeric_value(neighbours_to_find) || neighbours_to_find < 1) {
        cerr << "Invalid neighbours to find specifier '" << nvalue << "'\n";
        usage(19);
      }

      if (neighbours_to_find > pool_size) {
        cerr << "You asked for " << neighbours_to_find
             << " neighbours, but pool only contains " << pool_size << ". Shortened\n";
        neighbours_to_find = pool_size;
      }

      if (verbose) {
        cerr << "A maximum of " << neighbours_to_find
             << " neighbours of each molecule will be found\n";
      }
    }

    neighbour_count.resize(neighbours_to_find + 1);
  }

  if (!cl.option_present('n') && !cl.option_present('T')) {
    neighbours_to_find = 1;
  }

  // If verbose and a threshold specified, they still need the neighbour characteristics

  if (0 == neighbour_count.elements_allocated()) {
    neighbour_count.resize(pool_size + 1);
  }

  if (allow_arbitrary_distances) {
    ;
  } else if (create_histogram) {
    histogram_nearnest_neighbour_distances.initialise(0.0, 1.0, 0.01);
  }

  fraction_as_string.set_leading_string(distance_tag);

  if (cl.option_present('j')) {
    int j;
    if (!cl.value('j', j) || j < 2) {
      cerr << "The output precision option (-j) must be a whole +ve number\n";
      usage(3);
    }

    fraction_as_string.initialise(0.0f, 1.0f, j);

    set_default_iwstring_float_concatenation_precision(j);

    if (verbose) {
      cerr << "Default float concatenation precision " << j << endl;
    }
  } else {
    fraction_as_string.initialise(0.0f, 1.0f, 4);
  }

  fraction_as_string.append_to_each_stored_string(">\n");

  int rc = 0;

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('p')) {
    if (!three_column_output_all_pairs(output)) {
      cerr << "Cannot create three column output form\n";
      rc = 4;
    }
  } else if (neighbours_to_find > 0) {
    rc = gfp_nearneighbours<FP_Number_Nbrs>(cl[0], pool_size, output);
  } else {
    rc = gfp_nearneighbours<FP_Fixed_Distance>(cl[0], pool_size, output);
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << pool_size << " fingerprints\n";
    cerr << "Neighbour distances for " << distance_stats.n() << " neighbours between "
         << distance_stats.minval() << " and " << distance_stats.maxval() << endl;
    if (distance_stats.n() > 1) {
      cerr << "Average " << distance_stats.average() << " variance "
           << distance_stats.variance();
    }
    cerr << endl;

    cerr << nearest_neighbour_distance_stats.n()
         << " nearest neighbour distances between "
         << nearest_neighbour_distance_stats.minval() << " and "
         << nearest_neighbour_distance_stats.maxval();
    if (nearest_neighbour_distance_stats.n() > 1) {
      cerr << " ave " << nearest_neighbour_distance_stats.average() << " std dev "
           << static_cast<float>(sqrt(nearest_neighbour_distance_stats.variance()));
    }
    cerr << endl;

    for (int i = 0; i < neighbour_count.number_elements(); i++) {
      if (neighbour_count[i]) {
        cerr << neighbour_count[i] << " molecules had " << i << " neighbours\n";
      }
    }
  }

  if (allow_arbitrary_distances) {  // don't have a histogram
    ;
  } else if (stream_for_nearest_neighbour_histogram.rdbuf()->is_open()) {
    do_write_histogram(stream_for_nearest_neighbour_histogram);
  }

  delete_gfp_file_scope_static_objects();

  if (rc) {
    return 0;
  }

  return 1;
}

int
main(int argc, char** argv)
{
  int rc = nearneighbours(argc, argv);

  return rc;
}
