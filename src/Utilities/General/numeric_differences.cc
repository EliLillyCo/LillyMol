// Numeric differences between files.

#include <math.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;
using std::endl;

static int verbose = 0;

typedef double diff_type_t;
static diff_type_t tolerance = 1.0e-07;

/*
  We can do absolute or relative errors
*/

static int absolute_differences = 0;

static int append_difference = 0;

/*
  Sometimes we just want to know the overall behaviour of the files
*/

static int summary_only = 0;

/*
  If this is a descriptor file, we process it specially
*/

static int descriptor_file = 0;

static int break_upon_finding_a_difference = 0;

static Accumulator<diff_type_t> differences;

/*
  Sometimes it is interesting to keep track of the columns in which differences
  are found
*/

static extending_resizable_array<int> differences_in_column;

static int differences_found = 0;

static int records_differing = 0;

static int records_read = 0;

/*
  For descriptor files, we know the number of columns in each row
*/

static int columns_in_input = 0;

/*
  When processing a descriptor file, we need the column headings
  Note that we access this array from 1 to columns_in_input
*/

static IWString* column_heading = NULL;

/*
  We also keep track of the differences in each column
*/

static Accumulator<diff_type_t>* column_differences = NULL;

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
  cerr << "Compares files for numeric differences. Files must be line by line comparable\n";
  cerr << "Might be obsolate, see jfilecompare instead\n";
  cerr << " -t <tolerance>     specify tolerance (default " << tolerance << ")\n";
  cerr << " -a                 compare absolute rather than relative differences\n";
  cerr << " -d                 append the actual numeric difference to the output\n";
  cerr << " -j                 files are descriptor files\n";
  cerr << " -b                 stop processing a record once a difference is found\n";
  cerr << " -s                 summary only - no per-record info\n";
  cerr << " -v                 verbose output\n";
  // clang-format on

  exit(rc);
}

static int
establish_column_headings(const const_IWSubstring& h1, const const_IWSubstring& h2) {
  assert(NULL == column_heading);
  assert(0 == columns_in_input);

  if (h1 != h2) {
    IWString ch1(h1);
    IWString ch2(h2);
    ch1.remove_leading_words(1);
    ch2.remove_leading_words(1);
    if (ch1 != ch2) {
      cerr << "Descriptor files must have identical header records\n";
      cerr << h1 << endl;
      cerr << h2 << endl;
      differences_found++;
      return 0;
    }
  }

  columns_in_input = h1.nwords();

  if (verbose) {
    cerr << "Descriptor files contain " << columns_in_input << " columns\n";
  }

  if (0 == columns_in_input) {
    cerr << "Huh, empty header records '" << h1 << "'\n";
    return 0;
  }

  column_heading = new IWString[columns_in_input + 1];

  column_differences = new Accumulator<diff_type_t>[columns_in_input + 1];

  int i = 0;
  int zcol = 1;

  const_IWSubstring token;

  while (h1.nextword(token, i)) {
    column_heading[zcol] = token;
    if (column_heading[zcol].length() < 8) {
      column_heading[zcol].extend(8, ' ');
    }
    zcol++;
  }

  assert(zcol - 1 == columns_in_input);

  return 1;
}

static int
do_append_difference(const const_IWSubstring& token1, const const_IWSubstring& token2,
                     std::ostream& output) {
  double v1;
  double v2;

  token1.numeric_value(v1);
  token2.numeric_value(v2);

  output << ' ' << fabs(v1 - v2);

  return output.good();
}

static int
are_the_same(int zcol, const_IWSubstring& token1, const_IWSubstring& token2) {
  diff_type_t v1, v2;

  if (!token1.numeric_value(v1) ||
      !token2.numeric_value(v2)) {  // cannot compare as numeric
    return 0;                       // these tokens are different
  }

  if (v1 ==
      v2) {  // probably a leading 0, or maybe less than the precision of diff_type_t
    return 1;
  }

  // cerr << "Token " << token1 << "' becomes " << v1 << " and '" << token2 << "' becomes
  // " << v2 << endl; cerr << "Adding " << (v1 - v2) << endl;

  diff_type_t d = fabs(v1 - v2);

  if (absolute_differences) {
    if (d < tolerance) {
      return 1;  // these are the same
    }

    differences.extra(v1 - v2);
    if (descriptor_file) {
      column_differences[zcol].extra(v1 - v2);
    }

    return 0;  // not the same
  }

  diff_type_t tmp = fabs(v1 + v2) * 0.5;  // to be used as the denominator

  if (0.0 == tmp) {
    tmp = fabs(v1);
  }

  if ((fabs(v1 - v2) / tmp) < tolerance) {
    return 1;
  }

  differences.extra(d);
  if (descriptor_file) {
    column_differences[zcol].extra(d);
  }

  return 0;  // not the same
}

/*
  The two records are different at the text level
*/

static int
ndiff(const_IWSubstring& buffer1, const_IWSubstring& buffer2, int zline,
      std::ostream& output) {
  if (buffer1.nwords() != buffer2.nwords()) {
    output << "Token count mismatch on line " << zline << endl;
    records_differing++;

    if (descriptor_file) {  // this is fatal for a descriptor file
      return 0;
    }

    return 1;
  }

  int i1 = 0;
  int i2 = 0;

  int zcol = 0;

  int differences_found_this_record = 0;

  const_IWSubstring token1, token2;

  IWString id;

  while (buffer1.nextword(token1, i1)) {
    (void)buffer2.nextword(token2, i2);

    zcol++;

    if (1 == zcol && descriptor_file) {
      id = token1;
    }

    if (token1 == token2) {
      continue;
    }

    if (are_the_same(zcol, token1, token2)) {
      continue;
    }

    if (!summary_only) {
      output << "Difference on line " << zline << " token " << zcol;
      if (descriptor_file) {
        output << ' ' << id << ' ' << column_heading[zcol];
      }

      output << " '" << token1 << "' vs '" << token2 << "'";

      if (append_difference) {
        do_append_difference(token1, token2, output);
      }

      output << '\n';
    }

    differences_found++;

    if (0 == differences_found_this_record) {
      records_differing++;
    }

    differences_found_this_record++;

    if (descriptor_file)  // keep looking for differences in this row
    {
      if (break_upon_finding_a_difference) {
        break;
      }
    } else if (verbose || summary_only) {
      differences_in_column[zcol]++;  // column numbering starting with column 1
    } else {  // don't bother looking for more differences in this record
      break;
    }
  }

  return output.good();
}

static int
ndiff(iwstring_data_source& input1, iwstring_data_source& input2, std::ostream& output) {
  const_IWSubstring buffer1, buffer2;

  while (input1.next_record(buffer1)) {
    if (!input2.next_record(buffer2)) {
      cerr << "Premature EOF on file2\n";
      return 0;
    }

    records_read++;

    if (records_read > 1) {  // beyond the first record
      ;
    } else if (descriptor_file) {
      if (!establish_column_headings(buffer1, buffer2)) {
        return 0;
      }

      continue;
    } else  // not a descriptor file
    {
      int nw = buffer1.nwords();
      differences_in_column.resize((10 * nw < 100 ? 100 : 10 * nw));
    }

    if (buffer1 == buffer2) {
      continue;
    }

    if (!ndiff(buffer1, buffer2, input1.lines_read(), output)) {
      return 0;
    }
  }

  if (input2.next_record(buffer2)) {
    cerr << "Premature EOF on file1\n";
  }

  return output.good();
}

static int
ndiff(const char* fname1, const char* fname2, std::ostream& output) {
  iwstring_data_source input1(fname1);
  if (!input1.ok()) {
    cerr << "Cannot open '" << fname1 << "'\n";
    return 0;
  }

  input1.set_strip_trailing_blanks();

  iwstring_data_source input2(fname2);
  if (!input2.ok()) {
    cerr << "Cannot open '" << fname2 << "'\n";
    return 0;
  }

  input2.set_strip_trailing_blanks();

  return ndiff(input1, input2, output);
}

static int
ndiff(int argc, char** argv) {
  Command_Line cl(argc, argv, "vt:ajbds");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('t')) {
    if (!cl.value('t', tolerance) || tolerance < 0.0) {
    }

    if (verbose) {
      cerr << "Tolerance set to " << tolerance << endl;
    }
  }

  if (cl.option_present('a')) {
    absolute_differences = 1;
    if (verbose) {
      cerr << "Will compare absolute differences\n";
    }
  }

  if (cl.option_present('j')) {
    descriptor_file = 1;

    if (verbose) {
      cerr << "Will process as a descriptor file\n";
    }

    break_upon_finding_a_difference = 0;
  }

  if (cl.option_present('b')) {
    break_upon_finding_a_difference = 1;

    if (verbose) {
      cerr << "Will break on finding a difference\n";
    }
  }

  if (cl.option_present('d')) {
    append_difference = 1;

    if (verbose) {
      cerr << "Will append differences\n";
    }
  }

  if (cl.option_present('s')) {
    if (cl.option_present('d')) {
      cerr << "The -s and -d options are mutually incompatible\n";
      usage(2);
    }

    summary_only = 1;

    if (verbose) {
      cerr << "Will only write a summary of differences\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  if (2 != cl.number_elements()) {
    cerr << "Takes exactly two file name arguments\n";
    usage(4);
  }

  (void)ndiff(cl[0], cl[1], std::cout);

  std::cout.flush();

  if (verbose || summary_only) {
    cerr << differences_found << " differences on " << records_differing << " records\n";
    cerr << "Differences between " << differences.minval() << " and "
         << differences.maxval();
    if (differences.n() > 1) {
      cerr << " ave " << differences.average();
    }
    cerr << endl;

    if (descriptor_file) {
      for (int i = 1; i <= columns_in_input; i++) {
        Accumulator<diff_type_t>& diffs = column_differences[i];
        if (0 == diffs.n()) {
          continue;
        }

        cerr << " column " << i << ' ' << column_heading[i] << ' ' << diffs.n()
             << " differences between " << diffs.minval() << " and " << diffs.maxval()
             << endl;
      }
    } else {
      for (int i = 0; i < differences_in_column.number_elements(); i++) {
        if (differences_in_column[i]) {
          cerr << differences_in_column[i] << " differences in column " << i << endl;
        }
      }
    }
  }

  if (NULL != column_heading) {
    delete[] column_heading;
  }

  if (NULL != column_differences) {
    delete[] column_differences;
  }

  return differences_found;
}

int
main(int argc, char** argv) {
  int rc = ndiff(argc, argv);

  return rc;
}
