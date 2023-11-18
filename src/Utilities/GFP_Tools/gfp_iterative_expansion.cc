/*
  Howard Broughton idea of iterative expansion around a set of needles
*/

#include <stdlib.h>

#include <algorithm>
#include <iomanip>
#include <memory>
#include <vector>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::cout;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static Tversky tversky;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

// static Accumulator_Int<int>* neighbour_count = nullptr;

static Report_Progress report_progress;

static similarity_type_t forever_discard = 2.0f;
static int nkeep = -1;

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
  cerr << "Howard Broughton's idea of iterative near neighbour expansion around a set of seeds\n";
  cerr << " -T <dist>      search radius for each iteration (use multiple -t options)\n";
  cerr << " -n <niter>     or just use one -T and <niter> expansions\n";
  cerr << " -Y <dist>      on first run through, discard all haystack molecules more than <dist> from any needle\n";
  cerr << " -k <n>         on first run through, keep only the <n> closest haystack molecules\n";
  cerr << " -r <n>         report progress every <n> computations done\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

struct NBR {
  int _generation;
  int _times_found;
  IWString _provenance;
  int _ancestor;
  const IW_General_Fingerprint* _parent;
  similarity_type_t _dist_from_parent;
};

class Fingerprint_Provenance
{
 private:
  const IW_General_Fingerprint* _fp;
  Fingerprint_Provenance* _next;

 public:
  Fingerprint_Provenance(const IW_General_Fingerprint*);
  ~Fingerprint_Provenance();

  void
  set_next(const IW_General_Fingerprint*);
};

Fingerprint_Provenance::Fingerprint_Provenance(const IW_General_Fingerprint* f)
{
  _fp = f;
  _next = nullptr;

  return;
}

Fingerprint_Provenance::~Fingerprint_Provenance()
{
  if (nullptr != _next) {
    delete _next;
  }

  return;
}

void
Fingerprint_Provenance::set_next(const IW_General_Fingerprint* f)
{
  if (nullptr != _next) {
    _next->set_next(f);
  }

  _next = new Fingerprint_Provenance(f);

  return;
}

template <typename T>
class Set_of_Fingerprints
{
 private:
  T* _pool;
  IWString* _smiles;
  int _n;

 public:
  Set_of_Fingerprints();
  ~Set_of_Fingerprints();

  int
  set_size(int);

  size_t
  size() const
  {
    return _n;
  }

  int
  build(const char*);
  int
  build(iwstring_data_source&);

  const T&
  operator[](const int ndx) const
  {
    return _pool[ndx];
  }

  T&
  operator[](const int ndx)
  {
    return _pool[ndx];
  }

  T*
  pointer_to_fp(const int ndx)
  {
    return _pool + ndx;
  }

  int
  copy_fingerprint_pointers(resizable_array<T*>&) const;

  void
  add_ptr(const int ndx, resizable_array<T*>& fp)
  {
    fp.add(_pool + ndx);
  }

  const IWString&
  smiles(int ndx) const
  {
    return _smiles[ndx];
  }

  const IWString&
  id(int ndx) const
  {
    return _pool[ndx].id();
  }
};

template <typename T>
Set_of_Fingerprints<T>::Set_of_Fingerprints()
{
  _pool = nullptr;
  _smiles = nullptr;
  _n = 0;

  return;
}

template <typename T>
Set_of_Fingerprints<T>::~Set_of_Fingerprints()
{
  if (nullptr != _pool) {
    delete[] _pool;
  }

  if (nullptr != _smiles) {
    delete[] _smiles;
  }

  return;
}

template <typename T>
int
Set_of_Fingerprints<T>::set_size(int s)
{
  assert(s > 0);

  if (nullptr != _pool) {
    delete[] _pool;
  }

  if (nullptr != _smiles) {
    delete[] _smiles;
  }

  _pool = new (std::nothrow) T[s];
  _smiles = new (std::nothrow) IWString[s];

  if (nullptr == _pool || nullptr == _smiles) {
    cerr << "Set_of_Fingerprints::set_size:cannot allocate " << s << " fingerprints\n";
    return 0;
  }

  _n = s;

  return 1;
}

template <typename T>
int
Set_of_Fingerprints<T>::build(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Set_of_Fingerprints::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

template <typename T>
int
Set_of_Fingerprints<T>::build(iwstring_data_source& input)
{
  if (0 == _n) {
    IWString pcn;
    pcn << '^' << identifier_tag;
    std::unique_ptr<re2::RE2> pcn_rx;
    iwre2::RE2Reset(pcn_rx, pcn);
    auto x = input.grep(*pcn_rx);
    if (!set_size(x)) {
      return 0;
    }
  }

  IW_TDT tdt;
  int ndx = 0;

  while (tdt.next(input)) {
    int fatal;
    if (_pool[ndx].construct_from_tdt(tdt, fatal)) {
      tdt.dataitem_value(smiles_tag, _smiles[ndx]);
      ndx++;
    } else if (!fatal) {
      continue;
    } else {
      cerr << "Fatal error building fingerprint " << ndx << endl;
      return 0;
    }
  }

  _n = ndx;

  return _n;
}

template <typename T>
int
Set_of_Fingerprints<T>::copy_fingerprint_pointers(resizable_array<T*>& fps) const
{
  if (!fps.resize(_n)) {
    cerr << "Set_of_Fingerprints::copy_fingerprint_pointers:cannot allocate " << _n
         << " indices\n";
    return 0;
  }

  for (auto i = 0; i < _n; ++i) {
    fps.add(_pool + i);
  }

  return 1;
}

static similarity_type_t
compute_the_distance(IW_General_Fingerprint& fp1, IW_General_Fingerprint& fp2,
                     const Tversky& tversky)
{
  if (tversky.active()) {
    return static_cast<similarity_type_t>(1.0) -
           fp1.IW_General_Fingerprint::tversky(fp2, tversky);
  }

  return static_cast<similarity_type_t>(1.0) - fp1.tanimoto(fp2);
}

static int distance_histogram[100];

static int
gfp_iterative_expansion_generation_0(
    Set_of_Fingerprints<IW_General_Fingerprint>& haystack, NBR* nbr_info,
    const resizable_array<IW_General_Fingerprint*>& needles,
    const similarity_type_t mythreshold, resizable_array<int>* neighbours)
{
  std::fill(distance_histogram, distance_histogram + 100, 0);

  const int pool_size = static_cast<int>(haystack.size());
  const auto nneedles = static_cast<int>(needles.size());

  similarity_type_t* mindist = new similarity_type_t[pool_size + pool_size];
  std::unique_ptr<similarity_type_t[]> free_mindist(mindist);
  std::fill(mindist, mindist + pool_size, static_cast<similarity_type_t>(2.0f));

  int rc = 0;

  for (auto i = 0; i < pool_size; ++i) {
    for (auto j = 0; j < nneedles; ++j) {
      const auto t = compute_the_distance(haystack[i], *(needles[j]), tversky);

      const auto ndx = static_cast<int>(t * 100.0f + 0.4999f);
      distance_histogram[ndx]++;

      if (t < mindist[i]) {
        mindist[i] = t;
      }

      if (t > mythreshold) {
        continue;
      }

      neighbours[j].add(i);
      rc++;
    }

    if (report_progress()) {
      const auto f =
          static_cast<float>(i) / static_cast<float>(pool_size);  // fraction processed
      cerr << "Generation 0 processing " << nneedles << " needles, completed " << i
           << " of " << pool_size << " fingerprints (" << std::setprecision(3) << f
           << "), found " << rc << " neighbours\n";
    }
  }

  if (0 == rc) {
    return 0;
  }

  int molecules_permanently_excluded = 0;

  if (forever_discard < 1.0f) {
    for (auto i = 0; i < pool_size; ++i) {
      if (mindist[i] > forever_discard) {
        nbr_info[i]._times_found = -1;
        molecules_permanently_excluded++;
      }
    }
  } else if (nkeep > 0) {
    std::copy(mindist, mindist + pool_size, mindist + pool_size);
    std::sort(mindist + pool_size, mindist + pool_size + pool_size);

    similarity_type_t distance_cutoff = mindist[pool_size + nkeep];

    if (verbose) {
      cerr << "In order to retain " << nkeep << " distance cutoff set to "
           << distance_cutoff << endl;
    }

    for (auto i = 0; i < pool_size; ++i) {
      if (mindist[i] > distance_cutoff) {
        nbr_info[i]._times_found = -1;
        molecules_permanently_excluded++;
      }
    }
  }

  if (verbose) {
    cerr << "Generation zero distance histogram\n";
    cerr << "Distance Molecules\n";
    for (auto i = 0; i < 100; ++i) {
      cerr << static_cast<float>(i) / 100.0f << ' ' << distance_histogram[i] << endl;
    }
    if (molecules_permanently_excluded) {
      cerr << "Permanently exclude " << molecules_permanently_excluded
           << " from all subsequent iterations\n";
    }
  }

  return rc;
}

static int
gfp_iterative_expansion(Set_of_Fingerprints<IW_General_Fingerprint>& haystack,
                        const NBR* nbr_info,
                        const resizable_array<IW_General_Fingerprint*>& needles,
                        const similarity_type_t mythreshold,
                        resizable_array<int>* neighbours)
{
  const int pool_size = static_cast<int>(haystack.size());
  const int nneedles = static_cast<int>(needles.size());

  int rc = 0;

  for (auto i = 0; i < pool_size; ++i) {
    if (nbr_info[i]._times_found) {  // already captured, or turned off
      continue;
    }

    IW_General_Fingerprint& fpi = haystack[i];

    for (auto j = 0; j < nneedles; ++j) {
      const auto t = compute_the_distance(fpi, *(needles[j]), tversky);

      if (t > mythreshold) {
        continue;
      }

      neighbours[j].add(i);
      rc++;
    }

    if (report_progress()) {
      const auto f =
          static_cast<float>(i) / static_cast<float>(pool_size);  // fraction processed
      cerr << "Processing " << nneedles << " needles, completed " << i << " of "
           << pool_size << " fingerprints (" << f << "), found " << rc << " neighbours\n";
    }
  }

  return rc;
}

static int
mask_off_neighbours(Set_of_Fingerprints<IW_General_Fingerprint>& haystack,
                    resizable_array<int>* neighbours, int generation, NBR* nbr_info,
                    resizable_array<IW_General_Fingerprint*>& fp,
                    resizable_array<int>& fp_ancestor)
{
  const int nneedles = static_cast<int>(fp.size());

  // We want to clean out the previous needles array (fp), but we need the ids and the
  // ancestor information, so just save them

  int* old_fp_ancestor = new (std::nothrow) int[nneedles];
  std::unique_ptr<int[]> free_old_fp_ancestor(old_fp_ancestor);

  if (nullptr == old_fp_ancestor) {
    cerr << "Cannot allocate " << nneedles << " identifier strings\n";
    return 0;
  }

  for (auto i = 0; i < nneedles; ++i) {
    old_fp_ancestor[i] = fp_ancestor[i];
  }

  fp_ancestor.resize_keep_storage(0);

  resizable_array<IW_General_Fingerprint*> newfp;
  newfp.resize(nneedles * 10);

  int rc = 0;

  for (auto i = 0; i < nneedles; ++i) {
    resizable_array<int>& nbrs = neighbours[i];

    const int n = static_cast<int>(nbrs.size());

    for (auto j = 0; j < n; ++j) {
      int x = nbrs[j];

      if (nbr_info[x]._times_found < 0) {  // permanently discarded
        continue;
      }

      nbr_info[x]._times_found++;

      if (1 == nbr_info[x]._times_found) {
        newfp.add(haystack.pointer_to_fp(x));

        fp_ancestor.add(old_fp_ancestor[i]);

        nbr_info[x]._generation = generation;
        nbr_info[x]._parent = fp[i];
        nbr_info[x]._ancestor = old_fp_ancestor[i];
        rc++;
      }
    }
    nbrs.resize_keep_storage(0);
  }

  fp.resize_keep_storage(0);

  for (int i = 0; i < static_cast<int>(newfp.size()); ++i) {
    fp.add(newfp[i]);
  }

  return rc;
}

static int
gfp_iterative_expansion(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vT:n:F:V:W:p:s:r:Y:k:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress report option (-r)\n";
      usage(1);
    }
  }

  if (cl.option_present('Y') && cl.option_present('k')) {
    cerr << "Can use only one of the -Y and -k options\n";
    usage(2);
  }

  if (cl.option_present('Y')) {
    if (!cl.value('Y', forever_discard) || forever_discard < 0.0 ||
        forever_discard > 1.0) {
      cerr << "The forever discard option (-Y) must be a valid distance\n";
      usage(2);
    }

    if (verbose) {
      cerr << "All molecules beyond " << forever_discard
           << " on first iteration permanently discarded\n";
    }
  }

  if (cl.option_present('k')) {
    if (!cl.value('k', nkeep) || nkeep < 1) {
      cerr << "The number of haystack molecules to keep at generation zero (-k) must be "
              "a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will keep only the " << nkeep << " most close molecules at generation 0\n";
    }
  }

  std::vector<similarity_type_t> threshold;

  if (!cl.option_present('T')) {
    cerr << "Must specify search radius via the -t option\n";
    usage(1);
  }

  if (cl.option_present('T')) {
    similarity_type_t t;
    for (auto i = 0; cl.value('T', t, i); ++i) {
      threshold.push_back(t);
    }
  }

  if (cl.option_present('n')) {
    if (1 != threshold.size()) {
      cerr << "If you specify a number of expansions (-n) you can have only one radius "
              "(-T)\n";
      usage(2);
    }

    int n;
    if (!cl.value('n', n) || n < 1) {
      cerr << "The number of iterations (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will perform " << n << " expansions\n";
    }

    while (static_cast<int>(threshold.size()) < n) {
      threshold.push_back(threshold[0]);
    }
  }

  if (!cl.option_present('p')) {
    cerr << "Must specify 'needles' via the -p option\n";
    usage(2);
  }

  Set_of_Fingerprints<IW_General_Fingerprint> needles;

  if (cl.option_present('p')) {
    const char* p = cl.option_value('p');

    if (!needles.build(p)) {
      cerr << "Cannot read needles fingerprints '" << p << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Read " << needles.size() << " fingerprints from '" << p << "'\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  Set_of_Fingerprints<IW_General_Fingerprint> haystack;

  if (cl.option_present('s')) {
    int s;

    if (!cl.value('s', s) || s < 2) {
      cerr << "The haystack size option (-s) must be a whole +ve number\n";
      usage(1);
    }

    if (!haystack.set_size(s)) {
      cerr << "Cannot allocate " << s << " fingerprints\n";
      return 1;
    }
  }

  if (!haystack.build(cl[0])) {
    cerr << "Cannot build haystack '" << cl[0] << "'\n";
    return 0;
  }

  const int haystack_size = static_cast<int>(haystack.size());

  auto nbrs = new (std::nothrow) NBR[haystack_size];
  std::unique_ptr<NBR[]> free_nbrs(nbrs);
  if (nullptr == nbrs) {
    cerr << "Cannot allocate " << haystack_size << " nbr structs\n";
    return 2;
  }

  for (auto i = 0; i < haystack_size; ++i) {
    nbrs[i]._times_found = 0;
  }

  // For each molecule ultimately found, we can keep track of the ancestor

  int size_of_nbr_array = 10 * needles.size();

  resizable_array<int>* neighbours =
      new (std::nothrow) resizable_array<int>[size_of_nbr_array];
  if (nullptr == neighbours) {
    cerr << "Cannot allocate " << size_of_nbr_array << " neighbour list arrays\n";
    return 2;
  }

  for (auto i = 0; i < size_of_nbr_array; ++i) {
    neighbours[i].resize(1000);
  }

  // First round is a little different because we use the fingerprints in the needles
  // object

  resizable_array<IW_General_Fingerprint*> fp;
  if (!needles.copy_fingerprint_pointers(fp)) {
    return 2;
  }

  auto nextra =
      gfp_iterative_expansion_generation_0(haystack, nbrs, fp, threshold[0], neighbours);
  if (0 == nextra) {
    cerr << "No neighbours found at generation 0, distance " << threshold[0] << endl;
    return 1;
  }

  resizable_array<int> fp_ancestor;

  for (auto i = 0; i < static_cast<int>(needles.size()); ++i) {
    fp_ancestor.add(i);
  }

  mask_off_neighbours(haystack, neighbours, 1, nbrs, fp, fp_ancestor);

  if (verbose) {
    cerr << "At generation 0, radius " << threshold[0] << " found " << nextra
         << " neighbours\n";
  }

  for (auto i = 1; i < static_cast<int>(threshold.size()); ++i) {
    if (nextra > size_of_nbr_array) {
      delete[] neighbours;
      size_of_nbr_array = 10 * nextra;

      neighbours = new (std::nothrow) resizable_array<int>[size_of_nbr_array];
      if (nullptr == neighbours) {
        cerr << "Cannot allocate " << size_of_nbr_array << " neighbour arrays\n";
        return 2;
      }
    }

    gfp_iterative_expansion(haystack, nbrs, fp, threshold[i], neighbours);

    nextra = mask_off_neighbours(haystack, neighbours, i + 1, nbrs, fp, fp_ancestor);

    if (verbose) {
      cerr << "At generation " << i << ", radius " << threshold[i] << " found " << nextra
           << " neighbours\n";
    }

    if (0 == nextra) {
      if (0 == verbose) {
        cerr << "No neighbours found at generation " << i << " dist " << threshold[i]
             << endl;
      }
      break;
    }
  }

  Accumulator<double> acc;

  int written = 0;
  for (auto i = 0; i < haystack_size; ++i) {
    if (nbrs[i]._times_found <= 0) {
      continue;
    }

    const auto ancestor = nbrs[i]._ancestor;

    const auto t = compute_the_distance(haystack[i], needles[ancestor], tversky);

    acc.extra(t);

    std::cout << haystack.smiles(i) << ' ' << haystack.id(i) << ' ' << nbrs[i]._generation
              << ' ' << nbrs[i]._times_found << ' ' << nbrs[i]._parent->id() << ' '
              << needles.id(ancestor) << endl;

    written++;
  }

  if (verbose && written > 0) {
    cerr << "From " << needles.size() << " needle fingerprints, found " << written
         << " neighbours\n";

    //  for (auto i = 0; i < threshold.size(); ++i)
    //  {
    //    cerr << neighbour_count[i] << " neighbours found at generation " << i << " dist
    //    " << threshold[i] << endl;
    //  }

    cerr << "Distances from ancestor between " << acc.minval() << " and " << acc.maxval()
         << " ave " << static_cast<float>(acc.average()) << endl;
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_iterative_expansion(argc, argv);

  return rc;
}
