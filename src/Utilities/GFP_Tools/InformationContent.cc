/*
  Compute the mutual information content of bits
*/

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <unordered_map>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"

#include "Utilities/GFP_Tools/sparsefp.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::unordered_map;

static int verbose = 0;

static IWString identifier_tag("PCN<");

static IWString fingerprint_tag;

static int show_zero_hits = 1;

#define SORT_BY_INFORMATION 1
#define SORT_BY_BIT_NUMBER 2

static int sort_bits = SORT_BY_INFORMATION;

/*
  Activity data comes from a separate file, or a token in the name
*/

static Numeric_Data_From_File<int> activity;

static IWString activity_tag;

static int activity_column = -1;

/*
  Fred sometimes doesn't want the data on multiple records.
*/

static int newline_in_output = 1;

/*
  When names have come from a descriptor computation, they may have
  underscores in them
*/

static int gsub_underscores_in_name = 0;

/*
  We can read sparse fingerprints in one of two forms.
  They can be ascii representation,   1,3,15-20;2048
  Or they can be non-colliding counted form
*/

static int sparse_ascii_representation = 0;

/*
  We may not want to print out all the bits
*/

static int nprint = 0;

/*
  The MACCS keys start at 0, but mk0 is actually maccs key 6
*/

static int print_offset = 0;

static int errors_encountered = 0;

static int errors_to_ignore = 0;

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

  cerr << " -F <tag>         tag for fingerprints\n";
  cerr << " -s <size>        number of molecules to be processed\n";
  cerr << " -A <fname>       activity file\n";
  cerr << " -A col=nn        activities in column <nn> in the activity file\n";
  cerr << " -u               translate underscores in the name when read in\n";
  cerr << " -o <number>      add <number> to bit numbers when doing output (suggest 6 for MACCS keys)\n";
//cerr << " -z               suppress display of bits with no hits\n";
  cerr << " -T NC            non-colliding counted fingerprints - Daylight representation\n";
  cerr << " -T sparse        sparse ASCII representation  \"1,4,55-66;99\"\n";
  cerr << " -p <number>      number of the best bits to print\n";
  cerr << " -y               suppress newline in output\n";
  cerr << " -B <number>      number of otherwise fatal errors to ignore\n";
  cerr << " -b               sort output by bit number - useful for sparse fingerprints\n";
  cerr << " -e               do NOT sort the output by significance\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

class FBFingerprint
{
 private:
  int _activity_class;
  IWString _id;

 public:
  FBFingerprint();

  const IWString&
  id() const
  {
    return _id;
  }

  int
  activity_class() const
  {
    return _activity_class;
  }

  int
  construct_from_tdt(const IW_TDT& tdt, int&);
};

FBFingerprint::FBFingerprint()
{
  _activity_class = -1;

  return;
}

int
FBFingerprint::construct_from_tdt(const IW_TDT& tdt, int& fatal)
{
  if (!tdt.dataitem_value(identifier_tag, _id)) {
    cerr << "FBFingerprint::construct_from_tdt: missing or invalid identifier '"
         << identifier_tag << "'\n";
    const_IWSubstring notused;
    if (tdt.dataitem("$FPG<", notused)) {
      cerr << "FBFingerprint::construct_from_tdt: ignored\n";
      fatal = 0;
      return 0;
    }

    fatal = 1;
    return 0;
  }

  if (gsub_underscores_in_name) {
    _id.gsub('_', ' ');
  }

  if (activity_tag.length()) {
    if (!tdt.dataitem_value(activity_tag, _activity_class)) {
      cerr << "FBFingerprint::construct_from_tdt: missing or invalid activity '"
           << activity_tag << "'\n";
      fatal = 1;
      return 0;
    }
  } else if (activity_column > 0) {
    const_IWSubstring token;
    if (!_id.word(activity_column, token)) {
      cerr << "FBFingerprint::construct_from_tdt: cannot extract token "
           << activity_column << " from '" << _id << "'\n";
      fatal = 1;
      return 0;
    }

    if (!token.numeric_value(_activity_class)) {
      cerr << "FBFingerprint::construct_from_tdt: invalid activity class '" << token
           << "'\n";
      fatal = 1;
      return 0;
    }
  } else {
    IW_STL_Hash_Map_int::const_iterator f = activity.find(_id);
    if (f == activity.end()) {
      cerr << "NO activity data for '" << _id << "', cannot continue\n";
      return 0;
    }

    _activity_class = (*f).second;
  }

  return 1;
}

class FBFingerprint_dense : public FBFingerprint, public IW_Bits_Base
{
 private:
 public:
  int
  construct_from_tdt(const IW_TDT& tdt, int&);

  int
  next_bit_set(int&, unsigned int&, int&) const;
};

class FBFingerprint_sparse : public FBFingerprint, public Sparse_Fingerprint
{
 private:
 public:
  int
  construct_from_tdt(const IW_TDT& tdt, int&);
};

int
FBFingerprint_dense::construct_from_tdt(const IW_TDT& tdt, int& fatal)
{
  if (!FBFingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  const_IWSubstring fp;
  if (!tdt.dataitem_value(fingerprint_tag, fp)) {
    cerr << "FBFingerprint_dense::construct_from_tdt: cannot extract '" << fingerprint_tag
         << "' from tdt\n";
    fatal = 1;
    return 0;
  }

  int rc;
  if (sparse_ascii_representation) {
    rc = IW_Bits_Base::construct_from_sparse_representation(fp);
  } else {
    fp.truncate_at_first(';');

    rc = IW_Bits_Base::construct_from_daylight_ascii_representation(fp);
  }

  if (0 == rc) {
    cerr << "FBFingerprint_dense::construct_from_daylight_ascii_representation: invalid "
            "fp\n";
    cerr << fp << endl;
    fatal = 1;
    return 0;
  }

  return 1;
}

int
FBFingerprint_dense::next_bit_set(int& istart, unsigned int& zbit, int& hits) const
{
  if (istart >= IW_Bits_Base::nbits()) {
    return 0;
  }

  int nb = IW_Bits_Base::nbits();

  while (istart < nb) {
    if (IW_Bits_Base::is_set(istart)) {
      zbit = istart;
      istart++;
      return 1;
    }

    istart++;
  }

  return 0;
}

int
FBFingerprint_sparse::construct_from_tdt(const IW_TDT& tdt, int& fatal)
{
  if (!FBFingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  const_IWSubstring fp;
  if (!tdt.dataitem_value(fingerprint_tag, fp)) {
    cerr << "FBFingerprint_sparse::construct_from_tdt: cannot extract '"
         << fingerprint_tag << "' from tdt\n";
    fatal = 1;
    return 0;
  }

  if (0 == fp.length()) {
    cerr << "FBFingerprint_sparse::construct_from_tdt: ignoring empty fingerprint\n";
    return 1;
  }

  if (!Sparse_Fingerprint::construct_from_daylight_ascii_representation(fp)) {
    cerr << "FBFingerprint_sparse::construct_from_daylight_ascii_representation: invalid "
            "fp\n";
    cerr << fp << endl;
    fatal = 1;
    return 0;
  }

  return 1;
}

static double factor = 1.0 / log10(0.5);

static double
compute_entropy(int n, double population)
{
  if (0 == n) {
    return 0.0;
  }

  double nd = static_cast<double>(n);

  return nd * log10(nd / population) * factor;
}

/*
  The two classes in our input
*/

static int c1 = -7;  // the numbers assigned to the classes
static int c2 = -7;

static int n1 = 0;
static int n2 = 0;

class EResults
{
 private:
  unsigned int _zbit;
  int _hit_1;
  int _hit_2;
  int _nothit_1;
  int _nothit_2;

  float _iab;

  float _explained_by_descriptor_split;

  //  private functions

  void
  _default_values();

 public:
  EResults();
  EResults(unsigned int);

  //  We rely on the copy operator being a bit-wise copy

  int
  set_bit(unsigned int b);

  int
  zbit() const
  {
    return _zbit;
  }

  int
  bit_hit_in_class(int c);

  int
  report(int, ostream&) const;

  float
  iab() const
  {
    return _iab;
  }

  int
  set_results(int, double, double);
};

EResults::EResults()
{
  _zbit = 0;
  _default_values();
}

EResults::EResults(unsigned int b)
{
  _zbit = b;
  _default_values();
}

void
EResults::_default_values()
{
  _hit_1 = 0;
  _hit_2 = 0;
  _nothit_1 = 0;
  _nothit_2 = 0;

  _iab = static_cast<float>(0.0);

  return;
}

int
EResults::bit_hit_in_class(int c)
{
  if (c1 == c) {
    _hit_1++;
  } else if (c2 == c) {
    _hit_2++;
  } else {
    cerr << "EResults::bit_hit_in_class: invalid class " << c << " must be " << c1
         << " or " << c2 << endl;
    abort();
  }

  return 1;
}

int
EResults::set_bit(unsigned int b)
{
  _zbit = b;

  _default_values();

  return 1;
}

/*
  We have read in all our molecules, and now know how many there are in the whole set
*/

int
EResults::set_results(int pool_size, double e1, double e2)
{
  // cerr << "Bit " << _zbit << ", " << pool_size << " molecules, " << _hit_1 << " hits in
  // class 1, " << _hit_2 << " hits in class 2\n";
  _nothit_1 = n1 - _hit_1;
  assert(_nothit_1 >= 0);

  _nothit_2 = n2 - _hit_2;
  assert(_nothit_2 >= 0);

  double T = static_cast<double>(pool_size);

  double eh = compute_entropy(_hit_1 + _hit_2, T);        // entropy of the hits
  double em = compute_entropy(_nothit_1 + _nothit_2, T);  // entropy of the misses

  double e = compute_entropy(_hit_1, T) + compute_entropy(_hit_2, T) +
             compute_entropy(_nothit_1, T) + compute_entropy(_nothit_2, T);

  // The output variables are all float

  _iab = (e1 + e2) + (eh + em) - e;

  if (_iab >= static_cast<float>(0.0)) {
    ;
  } else if (_iab >= static_cast<float>(-0.01)) {
    _iab = static_cast<float>(0.0);
  } else {
    cerr << "Fatal error, e = " << e << " _iab " << _iab << endl;
    assert(NULL == "this is very bad");
  }

  _explained_by_descriptor_split = _iab / (e1 + e2);

  return 1;
}

int
EResults::report(int pool_size, ostream& os) const
{
  double T = static_cast<double>(pool_size);

  double eh = compute_entropy(_hit_1 + _hit_2, T);        // entropy of the hits
  double em = compute_entropy(_nothit_1 + _nothit_2, T);  // entropy of the misses

  os << "Bit " << (_zbit + print_offset) << " hit " << (_hit_1 + _hit_2)
     << " times. hit1 " << _hit_1 << " hit2 " << _hit_2 << ", nothit1 " << _nothit_1
     << ", nothit2 " << _nothit_2 << ". EH " << static_cast<float>(eh) << " EM "
     << static_cast<float>(em) << " IAB " << _iab;
  if (newline_in_output) {
    os << endl;
  }

  float bits_per_compound = _iab / T;
  float descriptor_info_that_predicts_activity = _iab / (eh + em);

  os << " bits/cpd " << bits_per_compound << " explained by split "
     << _explained_by_descriptor_split << " predicts "
     << descriptor_info_that_predicts_activity << endl;

  return os.good();
}

/*static int
eresults_comparitor (const void * p1, const void * p2)
{
  const EResults * e1 = (const EResults *) p1;
  const EResults * e2 = (const EResults *) p2;

  float iab1 = e1->iab ();
  float iab2 = e2->iab ();

  if (iab1 < iab2)
    return 1;
  if (iab1 > iab2)
    return -1;

  return 0;
}*/

static int
increment_class_counters(int c)
{
  if (c == c1) {
    n1++;
    return 1;
  }

  if (c == c2) {
    n2++;
    return 1;
  }

  if (c1 < 0) {
    c1 = c;
    n1 = 1;
    return 1;
  }

  if (c2 < 0) {
    c2 = c;
    n2 = 1;
    return 1;
  }

  cerr << "Invalid class " << c << " classes are " << c1 << " and " << c2 << endl;

  return 0;
}

/*
  We need to keep track of all the bits encountered
*/

static resizable_array_p<EResults> bits_found;

#ifdef __GNUG__
template class resizable_array_p<EResults>;
template class resizable_array_base<EResults*>;
#endif

/*
  We need a quick means of getting the results for a particular bit
*/

// static hash_map<unsigned int, EResults *> bit_hash;
static unordered_map<unsigned int, EResults*> bit_hash;

static int
eresults_comparitor(EResults* const* v1, EResults* const* v2)
{
  const EResults* e1 = *v1;
  const EResults* e2 = *v2;

  float iab1 = e1->iab();
  float iab2 = e2->iab();

  if (iab1 < iab2) {
    return 1;
  }
  if (iab1 > iab2) {
    return -1;
  }

  return 0;
}

static int
bit_number_comparitor(EResults* const* v1, EResults* const* v2)
{
  const EResults* e1 = *v1;
  const EResults* e2 = *v2;

  int b1 = e1->zbit();
  int b2 = e2->zbit();

  if (b1 < b2) {
    return -1;
  }
  if (b1 > b2) {
    return 1;
  }

  return 0;
}

static int fingerprints_read = 0;

static int
InformationContent(ostream& output)
{
  int nb = bits_found.number_elements();

  if (verbose) {
    cerr << "Input contains " << nb << " different bits in " << fingerprints_read
         << " fingerprints\n";
    cerr << "N1 = " << n1 << ", n2 = " << n2 << endl;
  }

  if (0 == fingerprints_read) {
    cerr << "No fingerprints read\n";
    return 0;
  }

  double e1 = compute_entropy(n1, static_cast<double>(fingerprints_read));
  double e2 = compute_entropy(n2, static_cast<double>(fingerprints_read));

  if (verbose) {
    cerr << "E1 = " << static_cast<float>(e1) << " E2 = " << static_cast<float>(e2)
         << endl;
  }

  Accumulator<float> acc;
  extending_resizable_array<int> bucket;
  bucket.resize(int(bits_found[0]->iab()) + 1);

  for (int i = 0; i < nb; i++) {
    EResults* eri = bits_found[i];

    eri->set_results(fingerprints_read, e1, e2);
    if (verbose) {
      acc.extra(eri->iab());
    }

    int b = static_cast<int>(eri->iab());

    bucket[b]++;
  }

  if (sort_bits == SORT_BY_INFORMATION) {
    bits_found.sort(eresults_comparitor);
  } else if (SORT_BY_BIT_NUMBER == sort_bits) {
    bits_found.sort(bit_number_comparitor);
  }

  if (0 == nprint) {
    nprint = nb;
  }

  for (int i = 0; i < nprint; i++) {
    const EResults* eri = bits_found[i];

    eri->report(fingerprints_read, cout);
  }

  if (verbose) {
    cerr << "IAB values between " << acc.minval() << " and " << acc.maxval();
    if (acc.n() > 1) {
      cerr << " ave " << acc.average();
    }
    cerr << endl;

    for (int i = 0; i < bucket.number_elements(); i++) {
      if (0 == bucket[i]) {
        continue;
      }

      int larger = 0;
      for (int j = i + 1; j < bucket.number_elements(); j++) {
        if (bucket[j]) {
          larger++;
        }
      }

      cerr << bucket[i] << " bits with IAB values at " << i << ", " << larger
           << " higher\n";
    }
  }

  return output.good();
}

template <typename T>
int
InformationContent(T& fp)
{
  fingerprints_read++;

  int c = fp.activity_class();

  if (!increment_class_counters(c)) {
    return 0;
  }

  int i = 0;
  unsigned int zbit;
  int hits;
  while (fp.next_bit_set(i, zbit, hits)) {
    EResults* er;

    unordered_map<unsigned int, EResults*>::const_iterator f = bit_hash.find(zbit);
    if (f == bit_hash.end()) {
      er = new EResults(zbit);
      bits_found.add(er);
      bit_hash[zbit] = er;

      f = bit_hash.find(zbit);
    } else {
      er = (*f).second;
    }

    int c = fp.activity_class();

    er->bit_hit_in_class(c);
  }

  return 1;
}

static int
InformationContent_sparse(iwstring_data_source& input, ostream& output)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    FBFingerprint_sparse fp;
    int fatal;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (!fatal) {
        continue;
      }

      cerr << "Invalid TDT, line " << input.lines_read() << endl;
      cerr << tdt;

      errors_encountered++;
      if (errors_encountered > errors_to_ignore) {
        return 0;
      }

      cerr << "Ignored " << errors_encountered << endl;
      continue;
    }

    if (!InformationContent(fp)) {
      cerr << "Fatal error, " << input.lines_read() << " lines read\n";
      return 0;
    }
  }

  return InformationContent(output);
}

static int
InformationContent_dense(iwstring_data_source& input, ostream& output)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    FBFingerprint_dense fp;
    int fatal;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (!fatal) {
        continue;
      }

      cerr << "INvalid TDT, line " << input.lines_read() << endl;
      cerr << tdt;

      errors_encountered++;

      if (errors_encountered > errors_to_ignore) {
        return 0;
      }

      cerr << "Ignored " << errors_encountered << endl;
      continue;
    }

    if (!InformationContent(fp)) {
      cerr << "Fatal error, " << input.lines_read() << " lines read\n";
      return 0;
    }
  }

  return InformationContent(output);
}

static int
InformationContent_sparse(const char* fname, ostream& output)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return InformationContent_sparse(input, output);
}

static int
InformationContent_dense(const char* fname, ostream& output)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return InformationContent_dense(input, output);
}

static int
InformationContent(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vF:A:o:zT:p:yuB:eb");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('F')) {
    cerr << "Must specify fingerprint tag via the -F option\n";
    usage(5);
  }

  if (cl.option_present('F')) {
    fingerprint_tag = cl.string_value('F');

    if (verbose) {
      cerr << "Fingerprints in '" << fingerprint_tag << "' dataitem processed\n";
    }
  }

  if (cl.option_present('A')) {
    IWString activity_fname;

    int i = 0;
    const_IWSubstring a;
    while (cl.value('A', a, i++)) {
      if (a.starts_with("col=")) {
        a.remove_leading_chars(4);
        int activity_column;
        if (!a.numeric_value(activity_column) || activity_column < 1) {
          cerr << "INvalid activity column '" << a << "'\n";
          usage(3);
        }

        if (verbose) {
          cerr << "The activity is the " << activity_column << " token in the name\n";
        }

        activity.set_column(activity_column - 1);
      } else if (activity_fname.length()) {
        cerr << "There can be only one activity file, '" << a << "' not allowed\n";
        usage(3);
      } else {
        activity_fname = a;
      }
    }

    if (0 == activity_fname.length()) {
      cerr << "Must specify activity file name via the -A option\n";
      usage(3);
    }

    if (!activity.read_data(activity_fname)) {
      cerr << "Cannot read activity data from '" << activity_fname << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Read " << activity.size() << " activity values from '" << activity_fname
           << "'\n";
    }
  }

  if (cl.option_present('B')) {
    if (!cl.value('B', errors_to_ignore) || errors_to_ignore < 0) {
      cerr << "The number of errors to ignore must be a whole non-negative number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will ignore a maximum of " << errors_to_ignore
           << " otherwise fatal input errors\n";
    }
  }

  if (cl.option_present('o')) {
    if (!cl.value('o', print_offset)) {
      cerr << "Invalid print offset (-o option)\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Bit numbers offset by " << print_offset << endl;
    }
  }

  if (cl.option_present('y')) {
    newline_in_output = 0;

    if (verbose) {
      cerr << "Suppressing newline in output\n";
    }
  }

  if (cl.option_present('u')) {
    gsub_underscores_in_name = 1;

    if (verbose) {
      cerr << "Will translate underscores to spaces in names\n";
    }
  }

  if (cl.option_present('z')) {
    show_zero_hits = 0;

    if (verbose) {
      cerr << "Will suppress display of zero hits\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', nprint) || nprint < 0) {
      cerr << "Invalid number of bits to print (-p option)\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will print only the best " << nprint << " bits\n";
    }
  }

  if (cl.option_present('e') && cl.option_present('b')) {
    cerr << "The -e (no sorting) and -b (sort by bit number) options are mutually "
            "exclusive\n";
    usage(5);
  }

  if (cl.option_present('e')) {
    sort_bits = 0;

    if (verbose) {
      cerr << "Bits not sorted\n";
    }
  }

  if (cl.option_present('b')) {
    sort_bits = SORT_BY_BIT_NUMBER;

    if (verbose) {
      cerr << "Will sort output by bit number\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  int rc;

  if (cl.option_present('T')) {
    const_IWSubstring t = cl.string_value('T');
    if ("NC" == t) {
      rc = InformationContent_sparse(cl[0], cout);
    } else if ("sparse" == t) {
      sparse_ascii_representation = 1;
      rc = InformationContent_dense(cl[0], cout);
    } else {
      cerr << "Unrecognised fingerprint type '" << t << "'\n";
      usage(5);
    }
  } else if (fingerprint_tag.starts_with("NC")) {
    rc = InformationContent_sparse(cl[0], cout);
  } else if (fingerprint_tag.starts_with("FP")) {
    rc = InformationContent_dense(cl[0], cout);
  } else {
    cerr << "Cannot infer fingerprint type from '" << fingerprint_tag
         << "', use the -T option\n";
    usage(4);
  }

  if (verbose) {
    if (errors_encountered) {
      cerr << errors_encountered << " read errors encountered\n";
    }
  }

  return 0;
}

int
main(int argc, char** argv)
{
  int rc = InformationContent(argc, argv);

  return rc;
}
