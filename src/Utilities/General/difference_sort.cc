/*
  We have two columns in a file, presumably obs and pred. Sort the
  file by the difference between the two
*/

#include <stdlib.h>

#include <cmath>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int is_descriptor_file = 0;

static int absolute_difference = 0;

static int col1 = -1;
static int col2 = -1;

static IWString d1, d2;

static resizable_array_p<IWString> header_records;

static int write_difference = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Sorts a file by the difference between two columns\n";
  cerr << " -c1 <col>      first  column specification\n";
  cerr << " -c2 <col>      second column specification\n";
  cerr << " -d1 <dname>    first  column specification\n";
  cerr << " -d2 <dname>    second column specification\n";
  cerr << " -j             input is a descriptor file\n";
  cerr << " -a             absolute differences\n";
  cerr << " -w             write the difference as an extra column\n";
  cerr << " -i <sep>       input file separator\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

template <typename T>
class File_Record_Template {
 private:
  IWString _zdata;
  T _v1;
  T _v2;
  T _diff;

 public:
  File_Record_Template();

  int initialise(const const_IWSubstring&, int, int, const char sep);

  T
  diff() const {
    return _diff;
  }

  const IWString&
  zdata() const {
    return _zdata;
  }
};

typedef File_Record_Template<double> File_Record;

template <typename T>
File_Record_Template<T>::File_Record_Template() {
  _v1 = static_cast<T>(0);
  _v2 = static_cast<T>(0);
  _diff = static_cast<T>(0);

  return;
}

template <typename T>
int
File_Record_Template<T>::initialise(const const_IWSubstring& buffer, int c1, int c2,
                                    const char sep) {
  int got_c2 = 0;

  int i = 0;
  const_IWSubstring token;

  for (int col = 0; buffer.nextword_single_delimiter(token, i, sep); col++) {
    if (col == c1 || col == c2) {
      T v;
      if (!token.numeric_value(v)) {
        cerr << "File_Record_Template::initialise:invalid data in column " << (col + 1)
             << '\n';
        return 0;
      }

      if (col == c1) {
        _v1 = v;
      } else {
        _v2 = v;
        got_c2 = 1;
      }
    }
  }

  if (!got_c2) {
    cerr << "File_Record_Template::initialise:did not find column " << (c2 + 1) << '\n';
    return 0;
  }

  _zdata = buffer;

  if (absolute_difference) {
    _diff = fabs(_v1 - _v2);
  } else {
    _diff = _v2 - _v1;
  }

  return 1;
}

class File_Record_Sorter {
 private:
 public:
  int operator()(const File_Record*, const File_Record*) const;
};

int
File_Record_Sorter::operator()(const File_Record* f1, const File_Record* f2) const {
  double d1 = f1->diff();
  double d2 = f2->diff();

  if (d1 < d2) {
    return -1;
  }
  if (d1 > d2) {
    return 1;
  }

  return 0;
}

static int
determine_columns(const const_IWSubstring& buffer, const IWString& d1, int& col1,
                  const IWString& d2, int& col2, const char sep) {
  int i = 0;
  const_IWSubstring token;

  int got1 = 0;
  int got2 = 0;

  for (int col = 0; buffer.nextword_single_delimiter(token, i, sep); col++) {
    if (token == d1) {
      col1 = col;
      got1++;
    } else if (token == d2) {
      col2 = col;
      got2++;
    }
  }

  if (1 == got1 && 1 == got2) {
    if (verbose) {
      cerr << "Descriptor '" << d1 << "' in column " << (col1 + 1) << ", '" << d2
           << "' in column " << (col2 + 1) << '\n';
    }
    return 1;
  }

  cerr << "Did not find both '" << d1 << "' and/or '" << d2
       << "' in header, cannot continue\n";
  return 0;
}

#ifdef __GNUG__
template class File_Record_Template<double>;
#endif

static int
gather_input_records(iwstring_data_source& input, const char sep,
                     resizable_array_p<File_Record>& zdata) {
  const_IWSubstring buffer;

  if (is_descriptor_file && 0 == header_records.number_elements()) {
    if (!input.next_record(buffer)) {
      cerr << "Cannot read header record\n";
      return 0;
    }

    header_records.add(new IWString(buffer));

    if (0 == d1.length()) {
      ;
    } else if (!determine_columns(buffer, d1, col1, d2, col2, sep)) {
      cerr << "Cannot identify columns from descriptor names '" << d1 << "' and '" << d2
           << "'\n";
      return 0;
    }
  }

  while (input.next_record(buffer)) {
    File_Record* f = new File_Record;

    if (!f->initialise(buffer, col1, col2, sep)) {
      cerr << "Cannot initialise '" << buffer << "'\n";
      return 0;
    }

    zdata.add(f);
  }

  return zdata.number_elements();
}

static int
gather_input_records(const char* fname, const char sep,
                     resizable_array_p<File_Record>& zdata) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gather_input_records(input, sep, zdata);
}

static int
difference_sort(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-c1=ipos-c2=ipos-d1=s-d2=s-a-j-w");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a')) {
    absolute_difference = 1;

    if (verbose) {
      cerr << "Differences reported will be absolute differences\n";
    }
  }

  if (cl.option_present('j')) {
    is_descriptor_file = 1;

    if (verbose) {
      cerr << "Files will be processed as descriptor files\n";
    }
  }

  int column_specifications = 0;

  if (cl.option_present("c1") && cl.option_present("c2")) {
    column_specifications++;
  }

  // cerr << cl.option_present("c1") << cl.option_present("c2") << cl.option_present("d1")
  // << cl.option_present("d2") << '\n';
  if (cl.option_present("d1") && cl.option_present("d2")) {
    column_specifications++;
    is_descriptor_file = 1;
  }

  if (1 != column_specifications) {
    cerr << "To specify columns to compare must specify -c1 and -c2, OR -d1 and -d2\n";
    usage(2);
  }

  if (cl.option_present("c1")) {
    if (!cl.value("c1", col1) || col1 < 1) {
      cerr << "Invalid column specification -c1\n";
      usage(2);
    }
    if (!cl.value("c2", col2) || col2 < 1 || col2 == col1) {
      cerr << "Invalid column specification -c2\n";
      usage(2);
    }

    if (col2 < col1) {
      iwswap(col1, col2);
    }

    if (verbose) {
      cerr << "Data in columns " << col1 << " and " << col2 << '\n';
    }

    col1--;
    col2--;
  }

  if (cl.option_present("d1")) {
    cl.value("d1", d1);
    cl.value("d2", d2);

    if (verbose) {
      cerr << "Difference columns are descriptors " << d1 << " and " << d2 << '\n';
    }
  }

  char input_separator = ' ';

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (!char_name_to_char(i)) {
      cerr << "Unrecognised input file separator '" << i << "'\n";
      return 1;
    }

    input_separator = i[0];
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  resizable_array_p<File_Record> zdata;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!gather_input_records(cl[i], input_separator, zdata)) {
      rc = i + 1;
      break;
    }
  }

  int n = zdata.number_elements();

  if (verbose) {
    cerr << "Read " << n << " data records\n";
  }

  if (0 == n) {
    cerr << "No data, cannot continue\n";
    return 3;
  }

  if (cl.option_present('w')) {
    write_difference = 1;

    if (verbose) {
      cerr << "Will create an extra column with the difference\n";
    }
  }

  File_Record_Sorter frs;

  zdata.iwqsort(frs);

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < header_records.number_elements(); i++) {
    output << *(header_records[i]);

    if (write_difference) {  // only works if 1 header record
      output << " DIFF";
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(4196);
  }

  for (int i = 0; i < n; i++) {
    const File_Record* f = zdata[i];

    output << f->zdata();

    if (write_difference) {
      output << ' ' << f->diff();
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  if (verbose) {
    Accumulator<double> acc;
    int zero_diff_count = 0;
    for (const File_Record* f : zdata) {
      if (f->diff() == 0.0) {
        ++zero_diff_count;
      } else {
        acc.extra(abs(f->diff()));
      }
    }

    cerr << "Read " << n << " records, " << zero_diff_count << 
            " no difference " << iwmisc::Fraction<float>(zero_diff_count, n) << '\n';
    cerr << acc.n() << " abs diffs btw " << acc.minval() << ' ' << acc.maxval() <<
            " ave " << static_cast<float>(acc.average()) << '\n';
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = difference_sort(argc, argv);

  return rc;
}
