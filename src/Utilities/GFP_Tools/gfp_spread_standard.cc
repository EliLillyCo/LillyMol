/*
  Spread implementation with fixed fingerprints
*/

#include <algorithm>
#include <iostream>
#include <optional>
#include <random>

#define USE_OMP
#ifdef USE_OMP
#include <omp.h>
#endif

#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "gfp_standard.h"
#include "smiles_id_dist.h"
#include "sparse_collection.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int choose_first_point_randomly = 0;

static int nsel = 0;

static float stop_once_distance_drops_below = 0.0f;

static int previously_selected_file_specified = 0;

static int previous_computed_distance_column = -1;

static IWString previously_computed_nn_distance_tag;
static IWString nn_smiles_tag("NNSMI<");
static IWString nn_id_tag("NNID<");

static int brief_output = 0;

static Numeric_Data_From_File<float> previously_computed_distances;

static Report_Progress report_progress;

static int squeeze = 0;
static int next_squeeze = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString mk_tag("FPMK<");
static IWString mk2_tag("FPMK2<");
static IWString iw_tag("FPIW<");

static Fraction_as_String fraction_as_string;

/*
  When we squeeze the pool, we must update the sdistance and nsn arrays, but
  we do not want to update the smiles and pcn arrays
  So each item must keep track of their initial index into those arrays
*/

class SSpread_Item : public GFP_Standard {
 private:
  int _selected;

  float _distance;

  int _nearest_previously_selected;

  int _initial_ndx;

 public:
  SSpread_Item();

  int selected() const {
    return _selected;
  }

  void set_selected(int s) {
    _selected = s;
  }

  float distance() const {
    return _distance;
  }

  int nearest_previously_selected() const {
    return _nearest_previously_selected;
  }

  void IsFirstSelected();

  void SetDistance(SSpread_Item& sel) {
    _distance = GFP_Standard::tanimoto_distance(sel);
    _nearest_previously_selected = sel.initial_ndx();
  }

  SSpread_Item& operator=(const SSpread_Item& rhs);
  // In retrospect, there seems to be no advantage to having
  // a move assignment operator.
  SSpread_Item& operator=(SSpread_Item&& rhs);

  void set_initial_ndx(int s) {
    _initial_ndx = s;
  }

  int initial_ndx() const {
    return _initial_ndx;
  }

  int SelectedOrZeroDistance() const;

  int NofitySelected(const SSpread_Item& rhs) {
    float d = GFP_Standard::tanimoto_distance(rhs);
    if (d >= _distance) {
      return 0;
    }
    _distance = d;
    _nearest_previously_selected = rhs.initial_ndx();
    return 1;
  }

  int NotifySelectedUseThreshold(const SSpread_Item& rhs) {
    std::optional<float> d = GFP_Standard::tanimoto_distance_if_less(rhs, _distance);
    if (!d) {
      return 0;
    }
    _distance = *d;
    _nearest_previously_selected = rhs.initial_ndx();
    return 1;
  }

  const IWString& smiles(const IWString* smiles) const {
    return smiles[_initial_ndx];
  }

  const IWString& pcn(const IWString* pcn) const {
    return pcn[_initial_ndx];
  }

  int DebugPrint(std::ostream& output) const;
};

SSpread_Item::SSpread_Item() {
  _selected = 0;
  _distance = 1.0;
  _initial_ndx = -1;
  _nearest_previously_selected = -1;
}

SSpread_Item&
SSpread_Item::operator=(const SSpread_Item& rhs) {
  _selected = rhs._selected;
  _distance = rhs._distance;
  _nearest_previously_selected = rhs._nearest_previously_selected;
  _initial_ndx = rhs._initial_ndx;
  GFP_Standard* me = this;
  const GFP_Standard* other = &rhs;
  *me = *other;
  return *this;
}

SSpread_Item&
SSpread_Item::operator=(SSpread_Item&& rhs) {
  _initial_ndx = rhs._initial_ndx;
  _selected = rhs._selected;
  _nearest_previously_selected = rhs._nearest_previously_selected;
  _distance = rhs._distance;
  GFP_Standard* me = this;
  const GFP_Standard* other = &rhs;
  *me = std::move(*other);
  return *this;
}

void
SSpread_Item::IsFirstSelected() {
  _selected = 1;
  _distance = 1.0f;
  _nearest_previously_selected = -1;
}

int
SSpread_Item::SelectedOrZeroDistance() const {
  if (_selected || 0.0f == _distance) {
    return 1;
  }

  return 0;
}

int
SSpread_Item::DebugPrint(std::ostream& output) const {
  output << "SSpread_Item:item " << _initial_ndx << " sel? " << _selected << " dist "
         << _distance << " nsn " << _nearest_previously_selected << '\n';
  return 1;
}

// A struct to describe the results of the computation.
struct Selected_Item {
  int _sel;
  int _nsn;
  float _dist;
};

static int pool_size = 0;

// Each of these will be sized for the job.
static SSpread_Item* fingerprints = nullptr;
static IWString* smiles = nullptr;
static IWString* pcn = nullptr;
static int* selected = nullptr;
static float* distances = nullptr;
static int* nearest_previously_selected = nullptr;
static float* weight = nullptr;
static Selected_Item* selected_item = nullptr;

// clang-format off
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
  cerr << "Spread implementation requiring MPR IW MK MK2\n";
  cerr << " -n <nsel>       number of items to select\n";
  cerr << " -A <fname>      file of previously selected fingerprints\n";
  cerr << " -p FILE=<fname> read weight values from <fname>\n";
  cerr << " -p COL=<col>    weight values are column <col> in the name\n";
  cerr << " -N <tag>        gfp_nearneighbours has been run and initial distances are in <tag> "
          "(also FILE=, COL=)\n";
  cerr << " -S rand         choose first item randomly\n";
  cerr << " -s <size>       size of fingerprint file\n";
  cerr << " -t <dist>       stop selection once distance drops below <dist>\n";
  cerr << " -h <n>          maximum number of OMP threads to use\n";
  cerr << " -r <n>          report progress every <n> items selected\n";
  cerr << " -q <n>          squeeze out already selected items evern <n> selections (def 1000)\n";
  cerr << " -b              brief output (smiles id dist)\n";
  cerr << " -v              verbose output\n";
// clang-format on

  exit(rc);
}

// clang-format on

static int
get_previously_computed_nearest_neighbour(
    const IWString& id, const IW_STL_Hash_Map_float& previously_computed_distances,
    float& d) {
  IW_STL_Hash_Map_float::const_iterator f = previously_computed_distances.find(id);

  if (f ==
      previously_computed_distances.end()) {  // OK if no previously computed distance
    return 1;
  }

  d = f->second;
  return 1;
}

static int
get_previously_computed_nearest_neighbour(
    const IW_TDT& tdt, const IWString& previously_computed_nn_distance_tag, float& d) {
  if (!tdt.dataitem_value(previously_computed_nn_distance_tag, d)) {
    return 0;
  }

  return 1;
}

static int
get_previously_computed_from_name_column(const IWString& id,
                                         int previous_computed_distance_column,
                                         float& d) {
  const_IWSubstring c;
  if (!id.word(previous_computed_distance_column, c)) {
    cerr << "get_previously_computed_from_name_column:cannot extract column "
         << (previous_computed_distance_column + 1) << " from '" << id << "'\n";
    return 0;
  }

  if (!c.numeric_value(d) || d < 0.0f) {
    cerr << "get_previously_computed_from_name_column:invalid numeric '" << id << "'\n";
    return 0;
  }

  return 1;
}

template <typename F>
void
do_output(const F* pool, const int isel, const int nsn, const float d,
          IWString_and_File_Descriptor& output) {
  const int ini = pool[isel].initial_ndx();

  if (brief_output) {
    output << smiles[ini] << ' ' << pcn[ini] << ' ' << d << "\n";
  } else {
    output << smiles_tag << smiles[isel] << ">\n";
    output << identifier_tag << pcn[isel] << ">\n";

    if (nsn >= 0) {
      const int n = pool[nsn].initial_ndx();

      output << smiles_tag << smiles[n] << ">\n";
      output << identifier_tag << pcn[n] << ">\n";
      output << distance_tag << d << ">\n";
    } else {
      output << smiles_tag << "*>\n";
      output << identifier_tag << "*>\n";
      output << distance_tag << "1>\n";
    }
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  return;
}

static void
set_selected_item(Selected_Item& s, const int isel, const int nsn, const float d) {
  s._sel = isel;
  s._nsn = nsn;
  s._dist = d;

  return;
}

#ifdef NO_LONGER_NEEDED
// Item pool[selected] has been selected. Update the info into 's`.
static void
set_selected_item(SSpread_Item* pool, int selected,
                  const int* nearest_previously_selected, const float* distances,
                  Selected_Item& s) {
  s._sel = pool[selected].initial_ndx();
  if (nearest_previously_selected[selected] >= 0) {
    s._nsn = pool[nearest_previously_selected[selected]].initial_ndx();
  }
  s._dist = distances[selected];
}
#endif

template <typename F>
int
mark_all_remaining_items_selected(F* pool, const int pool_size,
                                  Selected_Item* selected_item, int ndx,
                                  const int* selected,
                                  const int* nearest_previously_selected,
                                  const float* distances) {
  int rc = 0;

  if (verbose) {
    cerr << "Processing " << pool_size << " zero distance items\n";
  }
  for (int i = 0; i < pool_size; ++i) {
    if (selected[i]) {
      continue;
    }

    rc++;

    const int n = pool[nearest_previously_selected[i]].initial_ndx();

    // cerr << "Zero distance item i = " << i << " ndx " << pool[i].initial_ndx() << "
    // prev " << n
    // << " dist " << distances[i] << '\n';

    set_selected_item(selected_item[ndx], pool[i].initial_ndx(), n, distances[i]);
    ndx++;
  }

  return rc;
}

int
RandomPoolMember(int pool_size) {
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> u(0, pool_size - 1);
  return u(generator);
}

template <typename F>
int
spread_with_weights(F* fingerprints, const int pool_size, const float* weight,
                    Selected_Item* selected_item) {
  int s;

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = 1;
    float dmax = distances[0];
    for (int i = 1; i < pool_size; i++) {
      if (weight[i] * distances[i] > dmax) {
        dmax = distances[i];
        s = i;
      }
    }
  } else {
    s = 0;
  }

  if (!previously_selected_file_specified) {
    // #pragma omp parallel for schedule(dynamic,256)
    for (int i = 0; i < pool_size; i++) {
      distances[i] = weight[i] * fingerprints[s].tanimoto_distance(fingerprints[i]);
      assert(distances[i] >= 0.0 && distances[i] <= 1.0);
      nearest_previously_selected[i] = s;
    }

    distances[s] = 1.0f;
    nearest_previously_selected[s] = -1;

    //  iw_write_array(distances, pool_size, "initial distances", cerr);
  }

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    //  cerr << " s = " << s << '\n';
    //  iw_write_array(spread_order, pool_size, "selected", cerr);

    set_selected_item(selected_item[items_selected], fingerprints[s].initial_ndx(),
                      fingerprints[nearest_previously_selected[s]].initial_ndx(),
                      distances[s]);

    selected[s] = 1;

    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

#pragma omp for schedule(dynamic, 256) nowait
      for (int i = 0; i < pool_size; i++) {
        if (selected[i] || 0.0f == distances[i]) {
          continue;
        }

        float d = weight[i] * fingerprints[s].tanimoto_distance(fingerprints[i]);
        if (d < distances[i]) {
          distances[i] = d;
          nearest_previously_selected[i] = s;
        }

        if (distances[i] > mydmax) {
          mydmax = distances[i];
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

    //  cerr << "Selected " << s << " dist " << distances[s] << '\n';

    if (nexts_shared < 0) {
      break;
    }

    s = nexts_shared;

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints\n";
    }
  }

  return 1;
}

template <typename T>
int
do_squeeze(T* pool, int& pool_size, int next_selected) {
  if (pool_size < 1000) {
    return next_selected;
  }

  int ndx = 0;
  int rc = -1;
  for (int i = 0; i < pool_size; ++i) {
    if (pool[i].selected()) {
      continue;
    }

    if (i == next_selected) {
      rc = ndx;
    }

    if (i == ndx) {
      ++ndx;
      continue;
    }

    pool[ndx] = pool[i];
    ++ndx;
  }

  pool_size = ndx;

  return rc;
}

template <typename F>
int
do_squeeze(F* pool, int& pool_size, int next_selected, int* selected,
           int* nearest_previously_selected, float* distances) {
  if (pool_size < 1000) {
    return next_selected;
  }

  cerr << "Performing squeeze pool size " << pool_size << '\n';
  int rc = -1;

  int ndx = 0;
  for (int i = 0; i < pool_size; ++i) {
    cerr << i << " selected[i] " << selected[i] << '\n';
    if (selected[i]) {
      cerr << " wtihin squeeze item " << i << " is selected\n";
    }
    if (selected[i]) {
      continue;
    }

    if (i == next_selected) {
      rc = ndx;
    }

    if (i == ndx) {  // No movement
      ndx++;
      continue;
    }

    pool[ndx] = pool[i];
    assert(pool[ndx].initial_ndx() == pool[i].initial_ndx());
    selected[ndx] = selected[i];  // Should always be 0.
    nearest_previously_selected[ndx] = nearest_previously_selected[i];
    distances[ndx] = distances[i];
    ndx++;
  }

  pool_size = ndx;
  cerr << " pool_size updated to " << pool_size << '\n';

  return rc;
}

template <typename T>
int
WriteSelectedBrief(const T* pool, int selected, IWString_and_File_Descriptor& output) {
  const T& sel = pool[selected];
  output << sel.smiles(smiles) << ' ' << sel.pcn(pcn) << ' ' << sel.distance() << '\n';
  output.write_if_buffer_holds_more_than(4096);
  return 1;
}

template <typename T>
int
WriteSelected(T* pool, int selected, IWString_and_File_Descriptor& output) {
  if (brief_output) {
    return WriteSelectedBrief(pool, selected, output);
  }

  const T& sel = pool[selected];

  output << smiles_tag << sel.smiles(smiles) << ">\n";
  output << identifier_tag << sel.pcn(pcn) << ">\n";

  if (sel.nearest_previously_selected() >= 0) {
    int ndx = sel.nearest_previously_selected();
    output << smiles_tag << smiles[ndx] << ">\n";
    output << identifier_tag << pcn[ndx] << ">\n";
    output << distance_tag << sel.distance() << ">\n";
  } else {
    output << smiles_tag << "*>\n";
    output << identifier_tag << "*>\n";
    output << distance_tag << "1>\n";
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

template <typename F>
int
spread4(F* pool, int& pool_size, IWString_and_File_Descriptor& output) {
  int s;

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = std::max_element(
            pool, pool + pool_size,
            [](const F& fp1, const F& fp2) { return fp1.distance() > fp2.distance(); }) -
        pool;
  } else {
    s = 0;
  }

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    WriteSelected(pool, s, output);
    pool[s].set_selected(1);

    // Variables shared between threads.
    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

// #pragma omp for schedule(dynamic, 1024) nowait
#pragma omp for schedule(static) nowait
      for (int i = 0; i < pool_size; i++) {
        if (pool[i].SelectedOrZeroDistance()) {
          continue;
        }

#define USE_THRESHOLD
#ifdef USE_THRESHOLD
        pool[i].NotifySelectedUseThreshold(pool[s]);
#else
        pool[i].NofitySelected(pool[s]);
#endif
        if (pool[i].distance() > mydmax) {
          mydmax = pool[i].distance();
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

#ifdef DEBUG_SPREAD4
    cerr << "Selected " << s << " dist " << pool[s].distance() << " nexts_shared "
         << nexts_shared << '\n';
#endif

    if (nexts_shared < 0) {  // we are ignoring nsel, unlikely to be a problem...
      for (int i = 0; i < pool_size; ++i) {
        if (!pool[i].selected()) {
          WriteSelected(pool, i, output);
        }
      }
      return items_selected;
    }

    s = nexts_shared;

    if (dmax_shared < stop_once_distance_drops_below) {
      return items_selected;
    }

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints, dist " << dmax_shared
           << '\n';
    }

    if (next_squeeze == items_selected) {
      // cerr << "Before squeeze s " << s << " pool_size " << pool_size << '\n';
      s = do_squeeze(pool, pool_size, s);
      // cerr << "After  squeeze s " << s << " pool_size " << pool_size << '\n';
      next_squeeze += squeeze;
    }
  }

  return nsel;
}

template <typename F>
int
spread3(F* pool, int& pool_size, Selected_Item* selected_item) {
  int s;

  // iw_write_array(distances, pool_size, "distances", cerr);

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = index_of_largest(pool_size, distances);
  } else {
    s = 0;
  }

  if (!previously_selected_file_specified) {
#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < pool_size; i++) {
      distances[i] = fingerprints[s].tanimoto_distance(fingerprints[i]);
      nearest_previously_selected[i] = s;
    }

    distances[s] = 1.0f;
    nearest_previously_selected[s] = -1;
  }

  //  iw_write_array(distances, pool_size, "initial distances", cerr);

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    //  cerr << " s = " << s << '\n';
    //  iw_write_array(spread_order, pool_size, "selected", cerr);

    set_selected_item(pool, s, nearest_previously_selected, distances,
                      selected_item[items_selected]);

    selected[s] = 1;

    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

#pragma omp for schedule(dynamic, 1024) nowait
      for (int i = 0; i < pool_size; i++) {
        if (selected[i] || 0.0f == distances[i]) {
          continue;
        }

#ifdef USE_THRESHOLD
        std::optional<float> d = pool[s].tanimoto_distance_if_less(pool[i], distances[i]);
        if (!d) {
          continue;
        }
        if (*d < distances[i]) {
          distances[i] = *d;
          nearest_previously_selected[i] = s;
        }
#else
        float d = pool[s].tanimoto_distance(pool[i]);
        if (d < distances[i]) {
          distances[i] = d;
          nearest_previously_selected[i] = s;
        }
#endif

        if (distances[i] > mydmax) {
          mydmax = distances[i];
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

    cerr << "Selected " << s << " dist " << distances[s] << " nexts_shared "
         << nexts_shared << '\n';

    if (nexts_shared < 0) {  // we are ignoring nsel, unlikely to be a problem...
      // for (int i = 0; i < pool_size; ++i) {
      //   cerr << i << " Selected " << selected[i] << " distance " << distances[i] << "
      //   nsn " << nearest_previously_selected[i] << '\n';
      // }
      return items_selected;
    }

    s = nexts_shared;

    if (dmax_shared < stop_once_distance_drops_below) {
      return items_selected;
    }

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints, dist " << dmax_shared
           << '\n';
    }

    if (next_squeeze == items_selected) {
      // cerr << "Before squeeze s " << s << " pool_size " << pool_size << '\n';
      s = do_squeeze(pool, pool_size, s, selected, nearest_previously_selected,
                     distances);
      // cerr << "After  squeeze s " << s << " pool_size " << pool_size << '\n';
      next_squeeze += squeeze;
    }
  }

  return nsel;
}

template <typename F>
int
spread2(F* pool, int& pool_size, Selected_Item* selected_item) {
  int s;

  // iw_write_array(distances, pool_size, "distances", cerr);

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = index_of_largest(pool_size, distances);
  } else {
    s = 0;
  }

  if (!previously_selected_file_specified) {
#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < pool_size; i++) {
      distances[i] = fingerprints[s].tanimoto_distance(fingerprints[i]);
      nearest_previously_selected[i] = s;
    }

    distances[s] = 1.0f;
    nearest_previously_selected[s] = -1;
  }

  //  iw_write_array(distances, pool_size, "initial distances", cerr);

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    //  cerr << " s = " << s << '\n';
    //  iw_write_array(spread_order, pool_size, "selected", cerr);

    set_selected_item(pool, s, nearest_previously_selected, distances,
                      selected_item[items_selected]);

    selected[s] = 1;

    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

#pragma omp for schedule(dynamic, 1024) nowait
      for (int i = 0; i < pool_size; i++) {
        if (selected[i] || 0.0f == distances[i]) {
          continue;
        }

        std::optional<float> d = pool[s].tanimoto_distance_if_less(pool[i], distances[i]);
        if (!d) {
          continue;
        }
        if (*d < distances[i]) {
          distances[i] = *d;
          nearest_previously_selected[i] = s;
        }

        if (distances[i] > mydmax) {
          mydmax = distances[i];
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

    //  cerr << "Selected " << s << " dist " << distances[s] << '\n';

    if (nexts_shared < 0)  // we are ignoring nsel, unlikely to be a problem...
    {
      mark_all_remaining_items_selected(pool, pool_size, selected_item,
                                        items_selected + 1, selected,
                                        nearest_previously_selected, distances);
      return pool_size;
    }

    s = nexts_shared;

    if (dmax_shared < stop_once_distance_drops_below) {
      return items_selected;
    }

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints, dist " << dmax_shared
           << '\n';
    }

    if (next_squeeze == items_selected) {
      s = do_squeeze(pool, pool_size, s, selected, nearest_previously_selected,
                     distances);
      next_squeeze += squeeze;
    }
  }

  return nsel;
}

/*
  Not quite as parallel as the version above
*/

#ifdef SLIGHTLY_SLOWER
static int
spread() {
  int s;

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = index_of_largest(pool_size, distances);
  } else {
    s = 0;
  }

  // cerr << "Initial point is " << s << " distance " << distances[s] << '\n';

  if (!previously_selected_file_specified) {
#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < pool_size; i++) {
      distances[i] = fingerprints[s].tanimoto_distance(fingerprints[i]);
      nearest_previously_selected[i] = s;
    }

    distances[s] = 1.0f;
    nearest_previously_selected[s] = -1;

    //  iw_write_array(distances, pool_size, "initial distances", cerr);
  }

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    //  cerr << " s = " << s << " distance " << distances[s] << '\n';
    //  iw_write_array(spread_order, pool_size, "selected", cerr);

    spread_order[items_selected] = s;
    selected[s] = 1;

#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < pool_size; i++) {
      if (selected[i]) {
        continue;
      }

      std::optional<float> d =
          fingerprints[s].tanimoto_distance_if_less(fingerprints[i], distances[i]);
      if (!d) {
        continue;
      }
      nearest_previously_selected[i] = s;
      distances[i] = d;
    }

    int nexts = -1;
    float dmax = -1.0f;

    for (int i = 0; i < pool_size; i++) {
      if (selected[i]) {
        continue;
      }

      if (distances[i] > dmax) {
        nexts = i;
        dmax = distances[i];
      }
    }

    //  cerr << "Selected " << s << " dist " << distances[s] << '\n';

    if (nexts < 0) {
      break;
    }

    s = nexts;

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints\n";
    }
  }

  return 1;
}
#endif

static int
assign_scaling_factors_from_token(const IWString* pcn, const int pool_size,
                                  int scaling_factor_column, float* weight) {
  scaling_factor_column--;
  cerr << "From column " << scaling_factor_column << '\n';

  for (int i = 0; i < pool_size; i++) {
    const_IWSubstring s;

    if (!pcn[i].word(scaling_factor_column, s)) {
      cerr << "Cannot extract column " << (scaling_factor_column + 1) << " from '"
           << pcn[i] << "'\n";
      return 0;
    }

    float v;

    if (!s.numeric_value(v) || v < 0.0f) {
      cerr << "Invalid scaling factor '" << pcn[i] << "'\n";
      return 0;
    }

    weight[i] = v;
  }

  return 1;
}

static int
assign_scaling_factors_from_file(iwstring_data_source& input, int col,
                                 const IWString* pcn, int verbose) {
  IW_STL_Hash_Map_int id_to_ndx;

  for (int i = 0; i < pool_size; i++) {
    IWString id = pcn[i];

    id.truncate_at_first(' ');

    id_to_ndx[id] = i;
  }

  const_IWSubstring buffer;

  int rc = 0;

  while (input.next_record(buffer)) {
    IWString id, s;

    if (!buffer.split(id, ' ', s)) {
      cerr << "Cannot tokenise fingerprint weight record '" << buffer << "'\n";
      return 0;
    }

    if (1 == col) {
      s.truncate_at_first(' ');
    } else {
      const_IWSubstring tmp(s);
      tmp.word(col - 1, s);
    }

    float w;

    if (!s.numeric_value(w) || w < 0.0f) {
      if (1 == input.lines_read()) {  // header record
        continue;
      }

      cerr << "Invalid fingerprint weight '" << buffer << "'\n";
      return 0;
    }

    IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find(id);

    if (f == id_to_ndx.end()) {
      if (verbose > 2) {
        cerr << "Weight specified for '" << id << "' but no fingerprint\n";
      }
      continue;
    }

    weight[(*f).second] = w;

    rc++;
  }

  if (rc != pool_size) {
    cerr << "Only got " << rc << " weight values for " << pool_size << " fingerprints\n";
    return 0;
  }

  return 1;
}

static int
assign_scaling_factors_from_file(const char* fname, int col, const IWString* pcn,
                                 int verbose) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open scaling data file '" << fname << "'\n";
    return 0;
  }

  return assign_scaling_factors_from_file(input, col, pcn, verbose);
}

static void
do_previously_selected(const IW_General_Fingerprint& gfp) {
  GFP_Standard sgfp;

  sgfp.build_molecular_properties(gfp.molecular_properties_integer());
  sgfp.build_iw(gfp[0]);
  sgfp.build_mk(gfp[1]);
  sgfp.build_mk2(gfp[2]);

#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < pool_size; i++) {
    float d = sgfp.tanimoto_distance(fingerprints[i]);

    assert(d >= 0.0f);

    if (d < distances[i]) {
      distances[i] = d;
    }
  }

  // iw_write_array(distances, pool_size, "After previous selections", cerr);

  return;
}

static int
do_previously_selected(iwstring_data_source& input) {
  IW_TDT tdt;

  while (tdt.next(input)) {
    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    do_previously_selected(gfp);
  }

  return 1;
}

static int
do_previously_selected(const char* fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected(input);
}

template <typename F>
int
read_pool(iwstring_data_source& input, F* pool, int& pool_size) {
  IW_TDT tdt;

  int ndx = 0;

  for (; tdt.next(input) && ndx < pool_size; ndx++) {
    pool[ndx].set_initial_ndx(ndx);

    tdt.dataitem_value(smiles_tag, smiles[ndx]);

    tdt.dataitem_value(identifier_tag, pcn[ndx]);

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

    fingerprints[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    fingerprints[ndx].build_iw(gfp[0]);
    fingerprints[ndx].build_mk(gfp[1]);
    fingerprints[ndx].build_mk2(gfp[2]);

    if (previously_computed_nn_distance_tag.length()) {
      if (!get_previously_computed_nearest_neighbour(
              tdt, previously_computed_nn_distance_tag, distances[ndx])) {
        return 0;
      }
    } else if (previously_computed_distances.size()) {
      if (!get_previously_computed_nearest_neighbour(
              pcn[ndx], previously_computed_distances, distances[ndx])) {
        return 0;
      }
    } else if (previous_computed_distance_column > 0) {
      if (!get_previously_computed_from_name_column(
              pcn[ndx], previous_computed_distance_column, distances[ndx])) {
        return 0;
      }
    }
  }

  for (int i = 0; i < ndx; ++i) {
    if (pool[i].initial_ndx() != i) {
      cerr << "Initial index wrong\n";
    }
  }

  if (ndx < 2) {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  pool_size = ndx;

  return pool_size;
}

static int
allocate_pool(int s, int scaling_factor_specified) {
  selected = new_int(s);
  fingerprints = new SSpread_Item[s];
  distances = new float[s];
  std::fill_n(distances, s, 1.0f);
  nearest_previously_selected = new_int(s, -1);
  smiles = new IWString[s];
  pcn = new IWString[s];
  selected_item = new Selected_Item[s];

  if (nullptr == selected || nullptr == fingerprints || nullptr == distances ||
      nullptr == nearest_previously_selected || nullptr == smiles) {
    cerr << "Yipes, could not allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  if (scaling_factor_specified) {
    weight = new_float(s, 1.0f);
    if (nullptr == weight) {
      cerr << "Cannot allocate weight array\n";
      return 0;
    }
  }

  pool_size = s;

  return 1;
}

static int
gfp_spread_standard(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:s:r:n:bh:p:S:N:t:q:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  set_report_fingerprint_status(0);

  set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    usage(2);
  }

  int scaling_factor_specified = 0;
  IWString scaling_factor_file_name;
  int scaling_factor_column_in_file = 1;
  int scaling_factor_column = -1;

  if (cl.option_present('p')) {
    IWString col;

    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++)) {
      if (p.starts_with("FILE=") && 0 == scaling_factor_file_name.length()) {
        p.remove_leading_chars(5);
        scaling_factor_file_name = p;
      } else if (p.starts_with("FCOL=")) {
        p.remove_leading_chars(5);
        if (!p.numeric_value(scaling_factor_column_in_file) ||
            scaling_factor_column_in_file < 1) {
          cerr << "Invalid scaling factor column number 'FCOL=" << p << "'\n";
          return 2;
        }
        scaling_factor_column_in_file--;
      } else if (p.starts_with("COL=") && 0 == col.length()) {
        p.remove_leading_chars(4);
        col = p;
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(4);
      }
    }

    if (scaling_factor_file_name.length() && col.length()) {
      cerr << "Cannot specify both scaling factor file name (FILE=) and scaling factor "
              "column "
              "(COL=)\n";
      usage(1);
    }

    if (col.length()) {
      if (!col.numeric_value(scaling_factor_column) || scaling_factor_column < 2) {
        cerr << "Invalid scaling factor column '" << col << "'\n";
        return 4;
      }

      scaling_factor_column--;
    }

    scaling_factor_specified = 1;
  }

  if (cl.option_present('N')) {
    int nset = 0;

    IWString fname;

    int i = 0;
    const_IWSubstring n;

    while (cl.value('N', n, i++)) {
      if (n.starts_with("FILE=")) {
        n.remove_leading_chars(5);
        fname = n;
        nset++;
      } else if (n.starts_with("COL=")) {
        n.remove_leading_chars(4);
        if (!n.numeric_value(previous_computed_distance_column) ||
            previous_computed_distance_column < 1) {
          cerr << "The previously computed distance column must be a whole +ve number\n";
          return 3;
        }
        if (verbose) {
          cerr << "Previously computed near neighbour distances in column "
               << previous_computed_distance_column << '\n';
        }

        previous_computed_distance_column--;
        nset++;
      } else if (0 == previously_computed_nn_distance_tag.length()) {
        previously_computed_nn_distance_tag = n;
        if (verbose) {
          cerr << "Previously computed near neighbour distances in the '"
               << previously_computed_nn_distance_tag << "' tag\n";
        }

        if (!previously_computed_nn_distance_tag.ends_with('<')) {
          previously_computed_nn_distance_tag += '<';
        }

        nset++;
      } else {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(3);
      }
    }

    if (nset > 1) {
      cerr << "Can specify just one of FILE=, COL= or tag for previously computed "
              "distances\n";
      usage(3);
    }

    if (fname.length()) {
      if (!previously_computed_distances.read_data(fname)) {
        cerr << "Cannot read previously computed distances from '" << fname << "'\n";
        return 4;
      }

      if (verbose) {
        cerr << "Read " << previously_computed_distances.size()
             << " previously computed nn distances from '" << fname << "'\n";
      }
    }

    previously_selected_file_specified = 1;
  }

  const char* fname = cl[0];

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 2;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    if (!allocate_pool(pool_size, scaling_factor_specified)) {
      return 3;
    }

    if (verbose) {
      cerr << "Job automatically sized to " << pool_size << " fingerprints\n";
    }
  } else {
    if (!allocate_pool(pool_size, scaling_factor_specified)) {
      return 3;
    }
  }

  if (!read_pool(input, fingerprints, pool_size)) {
    cerr << "Cannot read fingerprints\n";
    return 0;
  }

  if (scaling_factor_column > 1) {
    if (!assign_scaling_factors_from_token(pcn, pool_size, scaling_factor_column,
                                           weight)) {
      cerr << "Cannot extract scaling factors from column " << scaling_factor_column
           << '\n';
      return 3;
    }
  } else if (scaling_factor_file_name) {
    if (!assign_scaling_factors_from_file(scaling_factor_file_name,
                                          scaling_factor_column_in_file, pcn, verbose)) {
      cerr << "Cannot assign scaling factors from '" << scaling_factor_file_name << "'\n";
      return 3;
    }
  }

#ifdef USE_OMP
  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h < 0) {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve "
              "number\n";
      usage(2);
    }

    omp_set_num_threads(h);
  }
#endif

  if (cl.option_present('A')) {
    const char* fname = cl.option_value('A');

    set_vector(distances, pool_size, 1.0f);

    if (!do_previously_selected(fname)) {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Finished processing already selected file\n";
    }

    previously_selected_file_specified = 1;
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if ("rand" == s) {
      choose_first_point_randomly = 1;

      if (verbose) {
        cerr << "Will choose first item randomly\n";
      }
    } else {
      cerr << "Unrecognised -S qualifier '" << s << "'\n";
      return 2;
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', nsel) || nsel < 2) {
      cerr << "The number if fingerprints to select must be a whole +ve number\n";
      usage(2);
    } else if (nsel > pool_size) {
      cerr << "ONly " << pool_size << " fingerprints read, nsel " << nsel
           << " impossible\n";
      return 3;
    }
  } else {
    nsel = pool_size;
  }

  if (cl.option_present('q')) {
    if (!cl.value('q', squeeze) || squeeze < 0) {
      cerr << "The squeeze every option (-q) must be a whole non negative number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will squeeze selected items every " << squeeze << " items selected\n";
    }

    if (squeeze == 0) {
      next_squeeze = -1;
    } else {
      next_squeeze = squeeze;
    }
  } else {
    squeeze = 1000;
    next_squeeze = 1000;
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', stop_once_distance_drops_below) ||
        stop_once_distance_drops_below < 0.0f || stop_once_distance_drops_below >= 1.0f) {
      cerr << "The stop selection distance option (-t) must be a valid distance\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will stop selection once distance drops below "
           << stop_once_distance_drops_below << '\n';
    }
  }

  fraction_as_string.initialise(0.0, 1.0, 4);
  fraction_as_string.set_leading_string("DIST<");
  fraction_as_string.append_to_each_stored_string(">\n");

  IWString_and_File_Descriptor output(1);

  cerr.flush();
  if (scaling_factor_specified) {
    spread_with_weights(fingerprints, pool_size, weight, selected_item);
  } else {
    spread4(fingerprints, pool_size, output);
  }

// #define DEBUG_SPREAD
#ifdef DEBUG_SPREAD
  for (int i = 0; i < nsel; i++) {
    cerr << " i = " << i << " selected " << spread_order[i] << " dist "
         << distances[spread_order[i]] << '\n';
  }
#endif

  if (verbose) {
    cerr << "Writing results\n";
  }

  if (cl.option_present('b')) {
    brief_output = 1;

    if (verbose) {
      cerr << "Brief output\n";
    }
  }

  return 0;

  if (brief_output) {
    for (int i = 0; i < nsel; ++i) {
      const Selected_Item& s = selected_item[i];

      output << smiles[s._sel] << ' ' << pcn[s._sel] << ' ' << s._dist << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  } else {
    for (int i = 0; i < nsel; ++i) {
      const Selected_Item& s = selected_item[i];

      const int id = s._sel;
      const int nsn = s._nsn;

      output << smiles_tag << smiles[id] << ">\n";
      output << identifier_tag << pcn[id] << ">\n";
      if (nsn >= 0) {
        output << smiles_tag << smiles[nsn] << ">\n";
        output << identifier_tag << pcn[nsn] << ">\n";
        output << distance_tag << s._dist << ">\n";
      } else {
        output << smiles_tag << "*>\n";
        output << identifier_tag << "*>\n";
        if (s._dist > 0.0f) {
          output << distance_tag << s._dist << ">\n";
        } else {
          output << distance_tag << "1>\n";
        }
      }
      output << "|\n";

      output.write_if_buffer_holds_more_than(4096);
    }
  }

  output.flush();

  delete[] fingerprints;
  delete[] smiles;
  delete[] pcn;
  delete[] selected;
  delete[] distances;
  delete[] nearest_previously_selected;
  delete[] weight;
  delete[] selected_item;
  delete_gfp_file_scope_static_objects();

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = gfp_spread_standard(argc, argv);

  return rc;
}
