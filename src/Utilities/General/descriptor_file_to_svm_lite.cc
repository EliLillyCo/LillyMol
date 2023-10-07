/*
  Converts a descriptor file to the format needed by svm-lite
*/

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int strip_leading_zeros = 0;

static int activity_column = 1;

typedef IW_STL_Hash_Map_String ID_to_Activity;

static IWDigits iwdigits;

static int ignore_no_activity_data = 0;

static IWString missing_value('.');

static int reject_records_with_missing_values = 1;

static int records_with_missing_values = 0;

static double ignore_close_to_zero = 0.0;

static int reading_categorical_data = 0;

static IWString write_header_fname;

static IWString existing_header;

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
  cerr << "Converts a descriptor file to the form needed by svm-lite\n";
  cerr << " -A <fname>     activity file\n";
  cerr << " -H <fname>     write the descriptor file header to <fname>\n";
  cerr << " -e             processing a test set, no -A file needed\n";
  cerr << " -U <fname>     check descriptors are in <fname> - created with -H option\n";
  cerr << " -z             strip leading zero's when comparing identifiers\n";
  cerr << " -C             activity data is categorical - skip numeric conversion\n";
  cerr << " -t <tol>       ignore values within <tol> of zero\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
descriptor_file_record_to_svm_lite(const const_IWSubstring& buffer,
                                   const ID_to_Activity& id_to_activity,
                                   const int* column_to_feature_number,
                                   IWString_and_File_Descriptor& output)
{
  const int length_of_output_on_entry = output.length();

  IWString id;
  const_IWSubstring token;
  int i = 0;

  buffer.nextword(id, i);

  if (id_to_activity.size() > 0) {
    ID_to_Activity::const_iterator f = id_to_activity.find(id);
    if (f == id_to_activity.end()) {
      const_IWSubstring tmp(id);
      tmp.remove_leading_chars('0');
      f = id_to_activity.find(tmp);

      if (f == id_to_activity.end()) {
        cerr << "No activity data for '" << id << "'\n";
        if (ignore_no_activity_data) {
          return 1;
        }

        return 0;
      }

      id = tmp;
    }

    output << (*f).second;
  } else {
    output << '0';
  }

  int missing_values_this_record = 0;

  for (int col = 1; buffer.nextword(token, i); col++) {
    if (missing_value == token) {
      missing_values_this_record++;
      if (reject_records_with_missing_values) {
        records_with_missing_values++;
        output.iwtruncate(length_of_output_on_entry);
        return 1;
      }
      continue;
    }

    if ('0' == token) {
      continue;
    }

    if (NULL != column_to_feature_number && column_to_feature_number[col] < 0) {
      continue;
    }

    double a;
    if (!token.numeric_value(a)) {
      cerr << "Invalid numeric value '" << token << "'\n";
      return 0;
    }

    if (0.0 == a) {
      continue;
    }

    if (ignore_close_to_zero > 0.0 && fabs(a) <= ignore_close_to_zero) {
      continue;
    }

    output << ' ';

    if (NULL == column_to_feature_number) {
      output << col;
    } else {
      output << column_to_feature_number[col];
    }

    output << ':' << token;
  }

  output << " # " << id << '\n';

  if (missing_values_this_record) {
    records_with_missing_values++;
  }

  return 1;
}

static int
do_write_header(const const_IWSubstring& buffer, IWString& write_header_fname)
{
  IWString_and_File_Descriptor output;
  if (!output.open(write_header_fname.null_terminated_chars())) {
    cerr << "do_write_header:cannot write to '" << write_header_fname << "'\n";
    return 0;
  }

  output << buffer << '\n';

  return 1;
}

/*
 */

static int
check_test_set_header_against_training_set_header(
    const const_IWSubstring& test_set_header, const IWString& training_set_header,
    int* column_to_feature_number)
{
  if (training_set_header.nwords() > test_set_header.nwords()) {
    cerr << "Training set had " << training_set_header.nwords() << " tokens\n";
    cerr << "Test     set has " << test_set_header.nwords()
         << " tokens, cannot continue\n";
    return 0;
  }

  IW_STL_Hash_Map_int descriptor_to_feature_number;

  int i = 0;
  IWString token;

  training_set_header.nextword(token, i);  // first token is the identifier

  int feature_number = 1;  // first descriptor gets feature number 1

  while (training_set_header.nextword(token, i)) {
    descriptor_to_feature_number[token] = feature_number;
    feature_number++;
  }

  int col = 1;  // see the 'for (int col = 1;...' loop above to see why se start this at 1
  i = 0;

  int test_descriptors_assigned = 0;
  int previous_feature_number_assigned =
      -1;  // feature numbers must be in ascending order

  while (test_set_header.nextword(token, i)) {
    IW_STL_Hash_Map_int::iterator f = descriptor_to_feature_number.find(token);

    if (f != descriptor_to_feature_number
                 .end())  // descriptor has a feature number in the training set
    {
      int feature_number = (*f).second;

      if (feature_number < previous_feature_number_assigned) {
        cerr << "Feature numbers out of order, column " << col << endl;
        return 0;
      }

      column_to_feature_number[col] = feature_number;
      (*f).second =
          feature_number;  // indicates that this descriptor has been assigned to a column
      test_descriptors_assigned++;
      previous_feature_number_assigned = feature_number;
    }

    col++;
  }

  if (test_descriptors_assigned != training_set_header.nwords() - 1) {
    cerr << "Training set had " << (training_set_header.nwords() - 1) << " descriptors\n";
    cerr << "Only found " << test_descriptors_assigned << " in the test set input file\n";

    for (IW_STL_Hash_Map_int::const_iterator i = descriptor_to_feature_number.begin();
         i != descriptor_to_feature_number.end(); i++) {
      if ((*i).second < 0) {
        cerr << "Descriptor '" << (*i).first
             << "' not found in test set descriptor file\n";
      }
    }

    return 0;
  }

  return 1;
}

static int
descriptor_file_to_svm_lite(iwstring_data_source& input,
                            const ID_to_Activity& id_to_activity,
                            IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record from descriptor file\n";
    return 0;
  }

  int* column_to_feature_number = NULL;

  if (write_header_fname.length()) {
    if (!do_write_header(buffer, write_header_fname)) {
      return 0;
    }
  } else if (existing_header.length()) {
    column_to_feature_number = new_int(buffer.nwords(), -1);

    if (!check_test_set_header_against_training_set_header(buffer, existing_header,
                                                           column_to_feature_number)) {
      return 0;
    }
  }

  std::unique_ptr<int> free_column_to_feature_number(
      column_to_feature_number);  // it IS safe to send a possibly NULL pointer to an
                                  // auto_ptr

  while (input.next_record(buffer)) {
    if (!descriptor_file_record_to_svm_lite(buffer, id_to_activity,
                                            column_to_feature_number, output)) {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    if (output.size() > 32768) {
      output.write_whole_blocks_shift_unwritten();
    }
  }

  return 1;
}

static int
descriptor_file_to_svm_lite(const char* fname, const ID_to_Activity& id_to_activity,
                            IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  return descriptor_file_to_svm_lite(input, id_to_activity, output);
}

static int
process_activity_record(const const_IWSubstring& buffer, ID_to_Activity& id_to_activity)
{
  if (buffer.nwords() < 2) {
    cerr << "Activity data must have at least two tokens\n";
    return 0;
  }

  int i = 0;
  IWString id;

  buffer.nextword(id, i);

  if (strip_leading_zeros) {
    id.remove_leading_chars('0');
  }

  assert(1 == activity_column);  // implement this sometime

  const_IWSubstring token;

  buffer.nextword(token, i);

  id_to_activity[id] = token;

  if (reading_categorical_data) {
    return 1;
  }

  float a;

  if (token.numeric_value(a)) {
    ;
  } else if (1 == id_to_activity.size())  // reading first record of file, probably header
  {
    //  cerr << "Skipping non-numeric header in activity file '" << buffer << "'\n";
    return 1;
  } else {
    cerr << "Non numeric activity value '" << buffer << "'\n";
    cerr << id_to_activity.size() << endl;
    return 0;
  }

  return 1;
}

static int
read_activity_data(iwstring_data_source& input, ID_to_Activity& id_to_activity)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (!process_activity_record(buffer, id_to_activity)) {
      cerr << "Fatal error reading activity data '" << buffer << "'\n";
      return 0;
    }
  }

  return id_to_activity.size();
}

static int
read_activity_data(const const_IWSubstring& fname, ID_to_Activity& id_to_activity)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_activity_data(input, id_to_activity);
}

static int
read_existing_header(const char* u, IWString& existing_header)
{
  iwstring_data_source input(u);
  if (!input.good()) {
    cerr << "Cannot open existing header file '" << u << "'\n";
    return 0;
  }

  if (!input.next_record(existing_header)) {
    cerr << "Cannot read header from existing header file\n";
    return 0;
  }

  return existing_header.length();
}

static int
descriptor_file_to_svm_lite(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:zgt:CeH:U:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  int n = 0;
  if (cl.option_present('e')) {
    n++;
  }
  if (cl.option_present('A')) {
    n++;
  }

  if (1 != n) {
    cerr << "Must specify exactly one of the -A or -e options\n";
    usage(4);
  }

  if (cl.option_present('z')) {
    strip_leading_zeros = 1;
    if (verbose) {
      cerr << "Will strip leading zero's from identifiers\n";
    }
  }

  if (cl.option_present('C')) {
    reading_categorical_data = 1;
    if (verbose) {
      cerr << "Reading categorical data\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', activity_column) || activity_column < 2) {
      cerr << "The activity column must be > 1\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Activity data in column " << activity_column << " of the activity file\n";
    }
  }

  if (cl.option_present('H') && cl.option_present('U')) {
    cerr << "The -H (create header file) and -U (use an existing header file) options "
            "cannot be used together\n";
    usage(4);
  }

  if (cl.option_present('H')) {
    write_header_fname = cl.string_value('H');
    if (verbose) {
      cerr << "Header record written to '" << write_header_fname << "'\n";
    }
  }

  if (cl.option_present('U')) {
    const char* u = cl.option_value('U');

    if (!read_existing_header(u, existing_header)) {
      cerr << "Cannot read existing header info from '" << u << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Header record written to '" << write_header_fname << "'\n";
    }
  }

  ID_to_Activity id_to_activity;

  if (cl.option_present('A')) {
    const_IWSubstring a = cl.string_value('A');

    if (!read_activity_data(a, id_to_activity)) {
      cerr << "Cannot read activity data in '" << a << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Read " << id_to_activity.size() << " activity records from '" << a
           << "'\n";
    }
  }

  if (cl.option_present('g')) {
    ignore_no_activity_data = 1;
    if (verbose) {
      cerr << "Will skip records for which there is no activity data\n";
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', ignore_close_to_zero) || ignore_close_to_zero <= 0.0) {
      cerr << "The tolerance for close to zero must be a +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will skip values within " << ignore_close_to_zero << " of zero\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!descriptor_file_to_svm_lite(cl[i], id_to_activity, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_to_svm_lite(argc, argv);

  return rc;
}
