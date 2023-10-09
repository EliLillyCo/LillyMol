/*
  Performs gfp_standard calculations based on smiles as input.
*/

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#define IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/maccskeys_fn5.h"
#include "Molecule_Tools/mpr.h"

#include "Utilities/GFP_Tools/gfp_standard.h"

using std::cerr;

static const char* prog_name = nullptr;

/*
  We need a lightweight structure for holding distance/id information.
*/

struct Distance_ID {
  float _d;
  int _ndx;  // index into pool
};

class Compare_DID
{
 private:
 public:
  bool
  operator()(const Distance_ID& lhs, const Distance_ID& rhs) const
  {
    return lhs._d > rhs._d;
  }
};

class Set_of_Fingerprints
{
 private:
  int _pool_size;
  GFP_Standard* _pool;
  IWString* _smiles_id;  // may hold just the ID's

  int _tmp[2048];

  Molecular_Properties_Generator _mpr;
  MACCSKeys _mk;

  Chemical_Standardisation _chemical_standardisation;

  //  private functions

  void
  _preprocess(Molecule& m);

  int
  _build_fingerprint(Molecule& m, GFP_Standard& gfp);

 public:
  Set_of_Fingerprints();
  ~Set_of_Fingerprints();

  int
  pool_size() const
  {
    return _pool_size;
  }

  int
  build(const char*);
  int
  build(data_source_and_type<Molecule>& input);

  const GFP_Standard&
  operator[](int s) const
  {
    return _pool[s];
  }

  const GFP_Standard*
  pool() const
  {
    return _pool;
  }

  const IWString&
  id(int s) const
  {
    return _smiles_id[s];
  }

  int
  find_id(const IWString& s) const;
};

Set_of_Fingerprints::Set_of_Fingerprints()
{
  _pool = nullptr;
  _smiles_id = nullptr;

  return;
}

Set_of_Fingerprints::~Set_of_Fingerprints()
{
  if (nullptr != _pool) {
    delete[] _pool;
  }

  if (nullptr != _smiles_id) {
    delete[] _smiles_id;
  }

  return;
}

void
Set_of_Fingerprints::_preprocess(Molecule& m)
{
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  _chemical_standardisation.process(m);

  return;
}

int
Set_of_Fingerprints::find_id(const IWString& s) const
{
  for (int i = 0; i < _pool_size; ++i) {
    if (s == _smiles_id[i]) {
      return i;
    }
  }

  return -1;
}

int
Set_of_Fingerprints::_build_fingerprint(Molecule& m, GFP_Standard& gfp)
{
  _preprocess(m);

  _mpr(m, gfp.molecular_properties());

  for (int i = 0; i < 8; ++i) {
    if (gfp.molecular_properties()[i] < 0 || gfp.molecular_properties()[i] > 255) {
      cerr << "out of range molecular property " << i << " value "
           << gfp.molecular_properties()[i] << '\n';
    }
  }

  IWMFingerprint iwfp;
  iwfp.construct_fingerprint(m);  // _iwfp is a bit vector
  gfp.build_iwfp(iwfp.bits(), iwfp.nset());

  _mk(m, _tmp);
  gfp.build_mk(_tmp, _mk.nbits());
  _mk.set_level_2_fingerprint(_tmp);
  gfp.build_mk2(_tmp, _mk.nbits());

  return 1;
}

int
Set_of_Fingerprints::build(const char* fname)
{
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);

  if (!input.good()) {
    cerr << "Set_of_Fingerprints::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

int
Set_of_Fingerprints::build(data_source_and_type<Molecule>& input)
{
  _pool_size = input.molecules_remaining();

  if (0 == _pool_size) {
    cerr << "Set_of_Fingerprints::build:no molecules in input\n";
    return 0;
  }

  _pool = new GFP_Standard[_pool_size];
  _smiles_id = new IWString[_pool_size];

  Molecule* m;

  int ndx = 0;

  for (; NULL != (m = input.next_molecule()); ndx++) {
    std::unique_ptr<Molecule> free_m(m);

    if (!_build_fingerprint(*m, _pool[ndx])) {
      cerr << "Set_of_Fingerprints::fatal error processing '" << m->name() << "'\n";
      return 0;
    }

    //  _smiles_id[ndx] << m->smiles() << ' ' << m->name();
    _smiles_id[ndx] << m->name();
  }

  _pool_size = ndx;

  return _pool_size;
}

class Search_Conditions
{
 private:
  int _nbrs;
  float _upper_distance_cutoff;
  float _lower_distance_cutoff;
  int _distance_constraints_present;
  int _avoid_self_neighbours;
  int _distance_matrix;

  resizable_array_p<IWString> _needle;

  // private functions

  void
  _usage(int rc);

 public:
  Search_Conditions();

  int
  build(Command_Line& cl);

  int
  nbrs() const
  {
    return _nbrs;
  }

  float
  upper_distance_cutoff() const
  {
    return _upper_distance_cutoff;
  }

  float
  lower_distance_cutoff() const
  {
    return _lower_distance_cutoff;
  }

  int
  distance_constraints_present() const
  {
    return _distance_constraints_present;
  }

  int
  avoid_self_neighbours() const
  {
    return _avoid_self_neighbours;
  }

  int
  distance_matrix() const
  {
    return _distance_matrix;
  }

  int
  number_needles() const
  {
    return _needle.number_elements();
  }

  const IWString&
  needle(const int i) const
  {
    return *_needle[i];
  }

  void
  set_nbrs(int s)
  {
    _nbrs = s;
  }

  void
  set_upper_distance_cutoff(float s)
  {
    _upper_distance_cutoff = s;
    _distance_constraints_present = 1;
  }

  void
  set_lower_distance_cutoff(float s)
  {
    _lower_distance_cutoff = s;
    _distance_constraints_present = 1;
  }

  void
  set_avoid_self_neighbours(int s)
  {
    _avoid_self_neighbours = s;
  }

  void
  set_distance_matrix(int s)
  {
    _distance_matrix = s;
  }

  int
  between_lower_and_upper(float s) const
  {
    if (s < _lower_distance_cutoff) {
      return 0;
    }
    if (s > _upper_distance_cutoff) {
      return 0;
    }

    return 1;
  }
};

Search_Conditions::Search_Conditions()
{
  _nbrs = 0;
  _upper_distance_cutoff = std::numeric_limits<float>::max();
  _lower_distance_cutoff = 0.0f;
  _distance_constraints_present = 0;
  _avoid_self_neighbours = 1;
  _distance_matrix = 0;

  return;
}

void
Search_Conditions::_usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  // clang-format on

  exit(rc);
}

int
Search_Conditions::build(Command_Line& cl)
{
  const auto verbose = cl.option_present('v');

  if (cl.option_present('d') || cl.option_present('D')) {
    if (cl.option_present('h') || cl.option_present('n') || cl.option_present('t') ||
        cl.option_present('T')) {
      cerr << "When computing a distance matrix, other options cannot be used\n";
      return 0;
    }

    if (cl.option_present('d') && cl.option_present('D')) {
      cerr << "Sorry, can have only one of -d or -D options\n";
      return 0;
    }

    if (cl.option_present('d')) {
      if (2 == cl.number_elements()) {
        cerr << "Search_Conditions::build:cannot do triangular matrix (-d) between two "
                "files\n";
        return 0;
      }
      _distance_matrix = 1;
    } else if (cl.option_present('D')) {
      _distance_matrix = 2;
    }

    if (verbose) {
      cerr << "Will write distance matrix\n";
    }

    return 1;
  }

  if (cl.option_present('h')) {
    _avoid_self_neighbours = 1;

    if (verbose) {
      cerr << "Will suppress self neighbours\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', _nbrs) || _nbrs < 1) {
      cerr << "The number of neighbours (-n) must be a whole +ve number\n";
      _usage(1);
    }
  }

  if (cl.option_present('T')) {
    if (!cl.value('T', _upper_distance_cutoff) || _upper_distance_cutoff < 0.0f ||
        _upper_distance_cutoff > 1.0f) {
      cerr << "The upper distance cutoff (-T) must be a valid distance\n";
      _usage(1);
    }

    _distance_constraints_present = 1;
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', _lower_distance_cutoff) || _lower_distance_cutoff < 0.0f ||
        _lower_distance_cutoff > 1.0f) {
      cerr << "The lower distance cutoff (-t) must be a valid distance\n";
      _usage(1);
    }

    _distance_constraints_present = 1;
  }

  if (_lower_distance_cutoff > _upper_distance_cutoff) {
    cerr << "Search_Conditions::build:inconsistent distance constraints "
         << _lower_distance_cutoff << " and " << _upper_distance_cutoff << '\n';
    return 0;
  }

  if (cl.option_present('e')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('e', s, i); ++i) {
      _needle.add(new IWString(s));
    }
  }

  if (!_distance_constraints_present && 0 == _nbrs) {
    _nbrs = 1;
  }

  return 1;
}

/*
  Compares sets of fingerprints. There will always be a haystack and one or more needles
*/

class GFP_Comparisons
{
 private:
  Set_of_Fingerprints _needles;
  Set_of_Fingerprints _haystack;

  //  All searches of the haystack go through this pointer. That way
  //  it can point either to _haystack or to _needles

  Set_of_Fingerprints* _hp;

  //  We need a structure for accumulating results

  Distance_ID* _did;

  //  What kind of output is required

  int _write_smiles;

  Fraction_as_String _fraction_as_string;

  int _output_precision;

  int _write_similarities;

  char _output_separator;

  int _verbose;

  Accumulator_Int<int> _acc_nbrs;
  Accumulator<float> _acc_shortest_dist;
  Accumulator<float> _acc_longest_dist;
  Accumulator<float> _acc_dist;

  //  private functions

  void
  _usage(int rc) const;

  int
  _tanimoto_nearest_nbr(const Set_of_Fingerprints& haystack, const Search_Conditions& sc,
                        IWString_and_File_Descriptor& output);
  int
  _tanimoto_nearest_nbr(const Set_of_Fingerprints& haystack, const GFP_Standard& sfp,
                        const IWString& needle_id, IWString_and_File_Descriptor& output);
  int
  _tanimoto_nearest_nbr_avoid_self_neighbours(const Set_of_Fingerprints& haystack,
                                              const GFP_Standard& sfp,
                                              const IWString& needle_id,
                                              IWString_and_File_Descriptor& output);

  int
  _tanimoto_no_distance_constraints(const Set_of_Fingerprints& haystack,
                                    const Search_Conditions& sc,
                                    IWString_and_File_Descriptor& output);
  int
  _tanimoto_no_distance_constraints(const Set_of_Fingerprints& haystack,
                                    const Search_Conditions& sc, const GFP_Standard& sfp,
                                    const IWString& needle_id,
                                    IWString_and_File_Descriptor& output);

  int
  _tanimoto_within_distance(const Set_of_Fingerprints& haystack,
                            const Search_Conditions& sc,
                            IWString_and_File_Descriptor& output);
  int
  _tanimoto_within_distance(const Set_of_Fingerprints& haystack, const GFP_Standard& gfp,
                            const IWString& needle_id, const float cutoff,
                            IWString_and_File_Descriptor& output);

  int
  _tanimoto_distance_constraints_present(const Set_of_Fingerprints& haystack,
                                         const Search_Conditions& sc,
                                         IWString_and_File_Descriptor& output);
  int
  _tanimoto_distance_constraints_present(const Set_of_Fingerprints& haystack,
                                         const Search_Conditions& sc,
                                         const GFP_Standard& sfp,
                                         const IWString& needle_id,
                                         IWString_and_File_Descriptor& output);

  int
  _do_search_vs_internal_needles(const Search_Conditions& sc,
                                 IWString_and_File_Descriptor& output);

  int
  _do_distance_matrix(const Search_Conditions& sc, IWString_and_File_Descriptor& output);
  int
  _row_of_distance_matrix(const GFP_Standard& sfp, const int istart, const int istop,
                          IWString_and_File_Descriptor& output);

  int
  _write_distance_matrix_header(const Set_of_Fingerprints& fp,
                                IWString_and_File_Descriptor& output) const;
  int
  _do_inter_file_distance_matrix(const Search_Conditions& sc,
                                 IWString_and_File_Descriptor& output) const;

  int
  _do_output(const Set_of_Fingerprints& haystack, const IWString& needle_id,
             const int nbrs, IWString_and_File_Descriptor& output);

  int
  _tanimoto_within_distance(const GFP_Standard&, float);
  int
  _tanimoto(const GFP_Standard& sfp, int nbrs);

 public:
  GFP_Comparisons();
  ~GFP_Comparisons();

  int
  build(const Command_Line& cl);

  int
  doit(const Search_Conditions& sc, IWString_and_File_Descriptor& output);

  int
  report(std::ostream& os) const;
};

GFP_Comparisons::GFP_Comparisons()
{
  _hp = &_haystack;

  _did = nullptr;

  _output_separator = ' ';

  _output_precision = 3;

  _write_similarities = 0;

  _verbose = 0;

  return;
}

GFP_Comparisons::~GFP_Comparisons()
{
  if (nullptr != _did) {
    delete[] _did;
  }

  return;
}

void
GFP_Comparisons::_usage(int rc) const
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Performs similarity calculations on smiles input(s)\n";
  cerr << " -n <nbrs>      number of neighbours to find\n";
  cerr << " -h             avoid self neighbours\n";
  cerr << " -t <dist>      lower distance cutoff\n";
  cerr << " -T <dist>      upper distance cutoff\n";
  cerr << " -d             write a lower triangular distance matrix\n";
  cerr << " -D             write a square distance matrix\n";
  // cerr << " -z             do not write molecules that have no neighbours\n";
  cerr << " -y             write similarities rather than distances\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

int
GFP_Comparisons::build(const Command_Line& cl)
{
  _verbose = cl.option_present('v');

  if (cl.option_present('y')) {
    _write_similarities = 1;
    if (_verbose) {
      cerr << "Will write similarities\n";
    }
  }

  set_global_aromaticity_type(Daylight);

  if (0 == cl.number_elements()) {
    cerr << "GFP_Comparisons::build:must specify one or two smiles files\n";
    _usage(1);
  }

  if (!_needles.build(cl[0])) {
    cerr << "GFP_Comparisons::build:cannot process '" << cl[0] << "'\n";
    return 0;
  }

  if (1 == cl.number_elements()) {
    _did = new Distance_ID[_needles.pool_size()];
    _hp = &_needles;
  } else {
    if (!_haystack.build(cl[1])) {
      cerr << "GFP_Comparisons::build:cannot process '" << cl[1] << "'\n";
      return 0;
    }
    _did = new Distance_ID[_haystack.pool_size()];
  }

  set_default_iwstring_float_concatenation_precision(3);

  return 1;
}

int
GFP_Comparisons::doit(const Search_Conditions& sc, IWString_and_File_Descriptor& output)
{
  if (sc.distance_matrix()) {
    return _do_distance_matrix(sc, output);
  }

  if (sc.number_needles()) {
    return _do_search_vs_internal_needles(sc, output);
  }

  if (!sc.distance_constraints_present() && 1 == sc.nbrs()) {
    return _tanimoto_nearest_nbr(*_hp, sc, output);
  }

  if (!sc.distance_constraints_present()) {
    return _tanimoto_no_distance_constraints(*_hp, sc, output);
  }

  return _tanimoto_distance_constraints_present(*_hp, sc, output);
}

int
GFP_Comparisons::_write_distance_matrix_header(const Set_of_Fingerprints& fp,
                                               IWString_and_File_Descriptor& output) const
{
  const auto pool_size = fp.pool_size();

  output << fp.id(0);

  for (auto i = 0; i < pool_size; ++i) {
    output << _output_separator << fp.id(i);
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
GFP_Comparisons::_do_distance_matrix(const Search_Conditions& sc,
                                     IWString_and_File_Descriptor& output)
{
  IWString tmp(_output_separator);

  _fraction_as_string.set_leading_string(tmp);
  _fraction_as_string.initialise(0.0f, 1.0f, _output_precision);

  if (_haystack.pool_size() > 0) {
    return _do_inter_file_distance_matrix(sc, output);
  }

  const GFP_Standard* pool = _needles.pool();

  const auto pool_size = _needles.pool_size();

  _write_distance_matrix_header(_needles, output);

  if (1 == sc.distance_matrix()) {
    for (auto i = 0; i < pool_size; ++i) {
      output << _needles.id(i);

      _row_of_distance_matrix(pool[i], 0, i, output);

      if (_write_similarities) {
        output << " 1\n";
      } else {
        output << " 0\n";
      }

      output.write_if_buffer_holds_more_than(4096);
    }
  } else {
    for (auto i = 0; i < pool_size; ++i) {
      output << _needles.id(i);

      _row_of_distance_matrix(pool[i], 0, i, output);

      if (_write_similarities) {
        output << " 1";
      } else {
        output << " 0";
      }

      _row_of_distance_matrix(pool[i], i + 1, pool_size, output);

      output << '\n';

      output.write_if_buffer_holds_more_than(4096);
    }
  }

  return 1;
}

int
GFP_Comparisons::_do_search_vs_internal_needles(const Search_Conditions& sc,
                                                IWString_and_File_Descriptor& output)
{
  const int needles = sc.number_needles();

  for (int i = 0; i < needles; ++i) {
    const IWString& ni = sc.needle(i);

    const int j = _haystack.find_id(ni);

    if (j < 0) {
      cerr << "GFP_Comparisons::_do_search_vs_internal_needles:cannot find '" << ni
           << "'\n";
      return 0;
    }

    //  _do_search_vs_internal_needle(sc, j, output);
  }

  return 1;
}

#ifdef IMPLEMENT_SOMETIME_IF_NEEDED
int
GFP_Comparisons::_do_search_vs_internal_needle(const Search_Conditions& sc,
                                               const GFP_Standard& sfp,
                                               IWString_and_File_Descriptor& output)
{
  return 1;
}
#endif

int
GFP_Comparisons::_row_of_distance_matrix(const GFP_Standard& sfp, const int istart,
                                         const int istop,
                                         IWString_and_File_Descriptor& output)
{
  const auto pool = _needles.pool();

  for (auto i = istart; i < istop; ++i) {
    auto d = sfp.tanimoto(pool[i]);

    if (!_write_similarities) {
      d = 1.0f - d;
    }

    _fraction_as_string.append_number(output, d);
  }

  return 1;
}

int
GFP_Comparisons::_do_inter_file_distance_matrix(
    const Search_Conditions& sc, IWString_and_File_Descriptor& output) const
{
  _write_distance_matrix_header(_haystack, output);

  const auto pool_size = _needles.pool_size();

  const auto nhaystack = _haystack.pool_size();

  for (auto i = 0; i < pool_size; ++i) {
    output << _needles.id(i);

    const GFP_Standard& sfp = _needles[i];

    for (auto j = 0; j < nhaystack; ++j) {
      auto d = sfp.tanimoto(_haystack[j]);

      if (!_write_similarities) {
        d = 1.0f - d;
      }

      _fraction_as_string.append_number(output, d);
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// static Compare_DID compare_did;

int
GFP_Comparisons::_tanimoto_no_distance_constraints(const Set_of_Fingerprints& haystack,
                                                   const Search_Conditions& sc,
                                                   IWString_and_File_Descriptor& output)
{
  for (auto i = 0; i < _needles.pool_size(); ++i) {
    const GFP_Standard& nfpi = _needles[i];

    _tanimoto_no_distance_constraints(haystack, sc, nfpi, _needles.id(i), output);
  }

  return 1;
}

int
GFP_Comparisons::_tanimoto_no_distance_constraints(const Set_of_Fingerprints& haystack,
                                                   const Search_Conditions& sc,
                                                   const GFP_Standard& sfp,
                                                   const IWString& needle_id,
                                                   IWString_and_File_Descriptor& output)
{
  const auto pool_size = haystack.pool_size();
  const GFP_Standard* pool = haystack.pool();

  auto nbrs = sc.nbrs();
  if (nbrs > pool_size) {
    nbrs = pool_size;
  }

  // First load up nbrs neighbours.

  int istart = 0;

  if (sc.avoid_self_neighbours()) {
    for (auto i = 0; i < pool_size; ++i) {
      const auto d = sfp.tanimoto(pool[i]);
      if (1.0f == d && needle_id == haystack.id(i)) {
        continue;
      }

      _did[istart]._d = d;
      _did[istart]._ndx = i;
      istart++;

      if (istart >= nbrs) {
        break;
      }
    }
  } else {
    for (auto i = 0; i < nbrs; ++i) {
      _did[i]._d = sfp.tanimoto(pool[i]);
      _did[i]._ndx = i;
    }
    istart = nbrs;
  }

  std::priority_queue<Distance_ID, std::vector<Distance_ID>, Compare_DID> pq(_did,
                                                                             _did + nbrs);

  Distance_ID tmp;

  for (auto i = istart; i < pool_size; ++i) {
    float s = sfp.tanimoto(pool[i]);

    if (1.0f == s && sc.avoid_self_neighbours() && needle_id == haystack.id(i)) {
      continue;
    }

    const Distance_ID& t = pq.top();

    //  cerr << "Top of queue " << t._d << ", compare " << s << '\n';

    if (s < t._d) {
      continue;
    }

    pq.pop();
    tmp._d = s;
    tmp._ndx = i;
    pq.push(tmp);
  }

  for (auto ndx = nbrs - 1; pq.size() > 0; ndx--) {
    const Distance_ID& d = pq.top();
    _did[ndx]._d = d._d;
    _did[ndx]._ndx = d._ndx;
    pq.pop();
  }

  return _do_output(haystack, needle_id, nbrs, output);
}

int
GFP_Comparisons::_tanimoto_nearest_nbr(const Set_of_Fingerprints& haystack,
                                       const Search_Conditions& sc,
                                       IWString_and_File_Descriptor& output)
{
  const auto nneedles = _needles.pool_size();

  const GFP_Standard* needle = _needles.pool();

  if (sc.avoid_self_neighbours()) {
    for (auto i = 0; i < nneedles; ++i) {
      _tanimoto_nearest_nbr_avoid_self_neighbours(haystack, needle[i], _needles.id(i),
                                                  output);
    }
  } else {
    for (auto i = 0; i < nneedles; ++i) {
      _tanimoto_nearest_nbr(haystack, needle[i], _needles.id(i), output);
    }
  }

  return 1;
}

int
GFP_Comparisons::_tanimoto_nearest_nbr(const Set_of_Fingerprints& haystack,
                                       const GFP_Standard& sfp, const IWString& needle_id,
                                       IWString_and_File_Descriptor& output)
{
  const auto pool_size = haystack.pool_size();

  const GFP_Standard* pool = haystack.pool();

  float max_sim = -1.0f;
  int idmax = -1;

  for (auto i = 0; i < pool_size; ++i) {
    const auto d = sfp.tanimoto(pool[i]);

    if (d > max_sim) {
      max_sim = d;
      idmax = i;
    }
  }

  assert(idmax >= 0);

  if (idmax < 0) {  // should not happen
    return 0;
  }

  _did[0]._ndx = idmax;
  _did[0]._d = max_sim;

  return _do_output(haystack, needle_id, 1, output);
}

int
GFP_Comparisons::_tanimoto_nearest_nbr_avoid_self_neighbours(
    const Set_of_Fingerprints& haystack, const GFP_Standard& sfp,
    const IWString& needle_id, IWString_and_File_Descriptor& output)
{
  const auto pool_size = haystack.pool_size();

  const GFP_Standard* pool = haystack.pool();

  float max_sim = -1.0f;
  int idmax = -1;

  for (auto i = 0; i < pool_size; ++i) {
    const auto d = sfp.tanimoto(pool[i]);

    if (1.0f == d && haystack.id(i) == needle_id) {
      continue;
    }

    if (d > max_sim) {
      max_sim = d;
      idmax = i;
    }
  }

  assert(idmax >= 0);

  if (idmax < 0) {  // should not happen
    return 0;
  }

  _did[0]._ndx = idmax;
  _did[0]._d = max_sim;

  return _do_output(haystack, needle_id, 1, output);
}

/*
  We assume that the _did array has been filled
*/

int
GFP_Comparisons::_do_output(const Set_of_Fingerprints& haystack,
                            const IWString& needle_id, const int nbrs,
                            IWString_and_File_Descriptor& output)
{
  if (_verbose) {
    _acc_nbrs.extra(nbrs);
  }

  if (0 == nbrs) {
    return 1;
  }

  if (_verbose) {
    _acc_shortest_dist.extra(_did[0]._d);
  }

  if (1 == nbrs) {
    const auto n = _did[0]._ndx;

    output << needle_id << _output_separator << haystack.id(n) << _output_separator;
    if (_write_similarities) {
      output << _did[0]._d << '\n';
    } else {
      output << (1.0f - _did[0]._d) << '\n';
    }

    output.write_if_buffer_holds_more_than(4096);

    if (_verbose) {
      _acc_longest_dist.extra(_did[0]._d);
      _acc_dist.extra(_did[0]._d);
    }

    return 1;
  }

  IWString tmp;
  tmp << needle_id << _output_separator;

  for (auto i = 0; i < nbrs; ++i) {
    output << tmp;
    const auto n = _did[i]._ndx;
    output << haystack.id(n) << _output_separator;
    if (_write_similarities) {
      output << _did[i]._d << '\n';
    } else {
      output << (1.0f - _did[i]._d) << '\n';
    }

    if (_verbose) {
      const auto d = _did[i]._d;

      if (0 == i) {
        _acc_shortest_dist.extra(d);
      } else if (i == nbrs - 1) {
        _acc_longest_dist.extra(d);
      }

      _acc_dist.extra(d);
    }
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
GFP_Comparisons::report(std::ostream& os) const
{
  os << "GFP_Comparisons::report:processed " << _acc_nbrs.n() << " needles\n";
  if (0 == _acc_nbrs.n()) {
    return 1;
  }

  os << " needles had between " << _acc_nbrs.minval() << " and " << _acc_nbrs.maxval()
     << " nbrs, ave " << static_cast<float>(_acc_nbrs.average()) << " nbrs\n";
  os << " shortest distances btw " << _acc_shortest_dist.minval() << " and "
     << _acc_shortest_dist.maxval() << " ave "
     << static_cast<float>(_acc_shortest_dist.average()) << '\n';
  os << " longtest distances btw " << _acc_longest_dist.minval() << " and "
     << _acc_longest_dist.maxval() << " ave "
     << static_cast<float>(_acc_longest_dist.average()) << '\n';
  os << " all      distances btw " << _acc_dist.minval() << " and " << _acc_dist.maxval()
     << " ave " << static_cast<float>(_acc_dist.average()) << '\n';

  return 1;
}

class Did_Comparator
{
 private:
 public:
  int
  operator()(const Distance_ID&, const Distance_ID&) const;
};

int
Did_Comparator::operator()(const Distance_ID& did1, const Distance_ID& did2) const
{
  if (did1._d < did2._d) {
    return 1;
  }
  if (did1._d > did2._d) {
    return -1;
  }

  return 0;
}

static Did_Comparator did_comparator;

int
GFP_Comparisons::_tanimoto_within_distance(const Set_of_Fingerprints& haystack,
                                           const Search_Conditions& sc,
                                           IWString_and_File_Descriptor& output)
{
  const auto nnedles = _needles.pool_size();

  const auto needles = _needles.pool();

  const auto d = sc.upper_distance_cutoff();

  for (auto i = 0; i < nnedles; ++i) {
    output << _needles.id(i);

    if (!_tanimoto_within_distance(haystack, needles[i], _needles.id(i), d, output)) {
      return 0;
    }
  }

  return 1;
}

int
GFP_Comparisons::_tanimoto_within_distance(const Set_of_Fingerprints& haystack,
                                           const GFP_Standard& gfp,
                                           const IWString& needle_id, float cutoff,
                                           IWString_and_File_Descriptor& output)
{
  const auto pool_size = haystack.pool_size();

  const auto* pool = haystack.pool();

  int did_ndx = 0;

  cutoff = 1.0f - cutoff;  // convert to similarity

  for (auto i = 0; i < pool_size; ++i) {
    float s = pool[i].tanimoto(gfp);
    if (s >= cutoff) {
      _did[did_ndx]._d = pool[i].tanimoto(gfp);
      _did[did_ndx]._ndx = i;
      did_ndx++;
    }
  }

  if (0 == did_ndx) {
    output << _output_separator << ".\n";
    return 1;
  }

  if (1 == did_ndx) {
    return _do_output(haystack, needle_id, 1, output);
  }

  iwqsort(_did, did_ndx, did_comparator);

  return _do_output(haystack, needle_id, did_ndx, output);
}

int
GFP_Comparisons::_tanimoto_distance_constraints_present(
    const Set_of_Fingerprints& haystack, const Search_Conditions& sc,
    IWString_and_File_Descriptor& output)
{
  const auto nneedles = _needles.pool_size();

  for (auto i = 0; i < nneedles; ++i) {
    _tanimoto_distance_constraints_present(haystack, sc, _needles[i], _needles.id(i),
                                           output);
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
GFP_Comparisons::_tanimoto_distance_constraints_present(
    const Set_of_Fingerprints& haystack, const Search_Conditions& sc,
    const GFP_Standard& sfp, const IWString& needle_id,
    IWString_and_File_Descriptor& output)
{
  const auto pool = haystack.pool();

  const auto pool_size = haystack.pool_size();

  int did_ndx = 0;

  for (auto i = 0; i < pool_size; ++i) {
    cerr << "Processing pool member " << i << '\n';
    const auto d = sfp.tanimoto(pool[i]);
    if (!sc.between_lower_and_upper(1.0f - d)) {
      continue;
    }

    if (1.0f == d && sc.avoid_self_neighbours() && haystack.id(i) == needle_id) {
      continue;
    }

    _did[did_ndx]._ndx = i;
    _did[did_ndx]._d = d;
    did_ndx++;
  }

  if (did_ndx > 1) {
    iwqsort(_did, did_ndx, did_comparator);
  }

  if (sc.nbrs() > 0 && did_ndx > sc.nbrs()) {
    did_ndx = sc.nbrs();
  }

  return _do_output(haystack, needle_id, did_ndx, output);
}

static int
runJob(const Search_Conditions& sc, GFP_Comparisons& server,
       IWString_and_File_Descriptor& output)
{
  int rc = server.doit(sc, output);

  output.flush();

  return rc;
}

static int
gfp_standalone_main(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vydDt:T:n:hzS:E:e:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    return 1;
  }

  int verbose = cl.option_count('v');

  initialise_fingerprints(cl, 0);

  GFP_Comparisons server;

  if (!server.build(cl)) {
    cerr << "Cannot initialise command line arguments\n";
    return 1;
  }

  Search_Conditions sc;

  if (!sc.build(cl)) {
    cerr << "Cannot build search conditions object\n";
    return 1;
  }

  int rc = 0;

  if (cl.option_present('S')) {
    const char* fname = cl.option_value('S');
    IWString_and_File_Descriptor output;
    if (!output.open(fname)) {
      cerr << "Cannot open " << fname << '\n';
      return 1;
    }
    if (!runJob(sc, server, output)) {
      rc = 1;
    }
  } else {
    IWString_and_File_Descriptor output(1);
    if (!runJob(sc, server, output)) {
      rc = 1;
    }
  }

  if (verbose) {
    server.report(cerr);
  }

  return rc;
}

#ifdef _WIN32_
#include "gfp_standalone.h"
extern "C" {
int
gfp_standalone_csharp(int argc, char** argv)
{
  return gfp_standalone_main(argc, argv);
}
}

int
Extern_Gfp_Standalone::extern_gfp_standalone(int argc, char** argv)
{
  int rc = gfp_standalone_main(argc, argv);

  return rc;
}
#else

#endif
int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_standalone_main(argc, argv);

  return rc;
}
