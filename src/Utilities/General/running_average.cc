/*
  Computes running averages in a descriptor file
*/

#include <math.h>

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int chunk_size = 1000;

static IWString missing_value('.');

static int number_descriptors = 0;
static IWString* descriptor = NULL;
static Accumulator<double>* acc;
static Accumulator<double>* global_acc;
static extending_resizable_array<double>* average;
static int* missing_value_count = NULL;
static int* column_xref = NULL;

static int report_cumulative_average = 0;

static int records_read = 0;

static int output_counts = 0;

static int report_min_and_max = 0;

static IWString* ids_this_chunk = NULL;

static int write_header_record = 0;

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
  cerr << " -c <records>     chunk size (default " << chunk_size << ")\n";
  cerr << " -d <desc>        descriptor to process\n";
  cerr << " -u               report cumulative average\n";
  cerr << " -o               write the X data as counts rather than bucket number\n";
  cerr << " -x               write the min and max values within each band\n";
  cerr << " -I               retain the ID's of the chunk-centre items\n";
  cerr << " -H               include a header record\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_write_header_record(std::ostream& os)
{
  return 1;
}

static int
write_centre_id(const IWString* ids_this_chunk, int records_read_this_chunk,
                std::ostream& output)
{
  int nmid = records_read_this_chunk / 2;

  output << ids_this_chunk[nmid];

  return output.good();
}

static int
process_chunk(int chunks_processed, int records_read_this_chunk, std::ostream& output)
{
  if (NULL != ids_this_chunk) {
    write_centre_id(ids_this_chunk, records_read_this_chunk, output);
  } else if (output_counts) {
    output << (chunks_processed * chunk_size);
  } else {
    output << chunks_processed;
  }

  for (int i = 0; i < number_descriptors; i++) {
    Accumulator<double>& ai = acc[i];
    if (0 == ai.n()) {
      output << ' ' << missing_value;
      if (report_min_and_max) {
        output << ' ' << missing_value << ' ' << missing_value;
      }
    } else {
      output << ' ' << ai.average();

      if (report_min_and_max) {
        output << ' ' << ai.minval() << ' ' << ai.maxval();
      }
    }

    average[i].add(ai.average());
    global_acc[i].extra(ai);

    ai.reset();

    if (!report_cumulative_average) {
      ;
    } else if (0 == global_acc[i].n()) {
      output << ' ' << missing_value;
    } else {
      output << ' ' << global_acc[i].average();
    }
  }

  output << '\n';

  return output.good();
}

static int
establish_column_headers(const const_IWSubstring& buffer)
{
  assert(number_descriptors > 0);
  assert(NULL == column_xref);

  int nw = buffer.nwords();
  column_xref = new_int(nw, -1);

  const_IWSubstring token;
  int i = 0;
  int col = 0;

  int descriptors_found = 0;

  while (buffer.nextword(token, i)) {
    if (0 == col)  // first column is the identifiers
    {
      col++;
      continue;
    }

    for (int j = 0; j < number_descriptors; j++) {
      if (token == descriptor[j]) {
        column_xref[col] = j;
        descriptors_found++;
        break;
      }
    }

    col++;
  }

  if (descriptors_found != number_descriptors) {
    cerr << "Yipes, could not find all the descriptors in the header\n";
    return 0;
  }

  if (verbose) {
    for (int i = 0; i < nw; i++) {
      int descriptor_number = column_xref[i];

      if (descriptor_number >= 0) {
        cerr << "Column " << (i + 1) << " is descriptor " << descriptor[descriptor_number]
             << '\n';
      }
    }
  }

  return 1;
}

static int
running_average(const const_IWSubstring& buffer, int records_read_this_chunk)
{
  int i = 0;
  const_IWSubstring token;
  int col = 0;

  if (NULL != ids_this_chunk) {
    const_IWSubstring tmp(buffer);
    tmp.truncate_at_first(' ');
    ids_this_chunk[records_read_this_chunk] = tmp;
  }

  while (buffer.nextword(token, i)) {
    int descriptor_number = column_xref[col];
    col++;

    if (descriptor_number < 0) {
      continue;
    }

    if (missing_value == token) {
      missing_value_count[descriptor_number]++;
      continue;
    }

    double d;
    if (!token.numeric_value(d)) {
      cerr << "INvalid numeric '" << token << "'\n";
      return 0;
    }

    acc[descriptor_number].extra(d);
  }

  return 1;
}

static int
running_average(iwstring_data_source& input, std::ostream& output)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (!establish_column_headers(buffer)) {
    return 0;
  }

  int chunks_processed = 0;
  int records_read_this_chunk = 0;
  while (input.next_record(buffer)) {
    records_read++;

    if (!running_average(buffer, records_read_this_chunk)) {
      cerr << "Fatal error on line " << input.lines_read() << '\n';
      return 0;
    }

    records_read_this_chunk++;

    if (records_read_this_chunk == chunk_size) {
      process_chunk(chunks_processed, records_read_this_chunk, output);
      chunks_processed++;
      records_read_this_chunk = 0;
    }
  }

  if (records_read_this_chunk) {
    process_chunk(chunks_processed, records_read_this_chunk, output);
  }

  return output.good();
}

static int
running_average(const char* fname, std::ostream& output)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return running_average(input, output);
}

static int
report_trend(const IWString& dname, const extending_resizable_array<double>& average,
             std::ostream& output)
{
  int chunks_read = average.number_elements();

  Accumulator<double> accy, accx;
  double r = 0.0;
  for (int i = 0; i < chunks_read; i++) {
    double a = average[i];

    accy.extra(a);
    accx.extra(static_cast<double>(i));
    r += i * a;
  }

  double xbar = accx.average();
  double ybar = accy.average();

  double nx1bx2b = static_cast<double>(chunks_read) * xbar * ybar;

  double v1 = accx.sum_of_squares() - chunks_read * xbar * xbar;
  double v2 = accy.sum_of_squares() - chunks_read * ybar * ybar;

  if (0.0 == v1 * v2) {
    output << "Trend in '" << dname << "' is *\n";
  } else {
    double rho = (r - nx1bx2b) / sqrt(v1 * v2);

    output << "Trend in '" << dname << "' is " << rho << '\n';
  }

  return output.good();
}

static int
running_average(int argc, char** argv)
{
  Command_Line cl(argc, argv, "c:d:vuoxIH");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(2);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    if (!cl.value('c', chunk_size) || chunk_size < 2) {
      cerr << "The chunk size must be a whole positive number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will report running averages for groups of " << chunk_size << " records\n";
    }
  }

  if (cl.option_present('I')) {
    ids_this_chunk = new IWString[chunk_size];
  }

  if (cl.option_present('u')) {
    report_cumulative_average = 1;

    if (verbose) {
      cerr << "WIll report cumulative averages\n";
    }
  }

  if (cl.option_present('o')) {
    output_counts = 1;

    if (verbose) {
      cerr << "Will write number of items processed as the bucket number\n";
    }
  }

  if (cl.option_present('x')) {
    report_min_and_max = 1;
    if (verbose) {
      cerr << "Will report min and max values\n";
    }
  }

  number_descriptors = cl.option_count('d');

  if (0 == number_descriptors) {
    cerr << "Must specify descriptors via the -d option\n";
    usage(3);
  }

  descriptor = new IWString[number_descriptors];
  acc = new Accumulator<double>[number_descriptors];
  global_acc = new Accumulator<double>[number_descriptors];
  average = new extending_resizable_array<double>[number_descriptors];
  missing_value_count = new_int(number_descriptors);

  for (int i = 0; i < number_descriptors; i++) {
    descriptor[i] = cl.string_value('d', i);
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!running_average(cl[i], std::cout)) {
      rc = i + 1;
      break;
    }
  }

  if (0 == verbose) {
    ;
  } else if (records_read < 3) {  // nothing to report
    cerr << "Insufficient data\n";
  } else {
    cerr << "Read " << records_read << " records\n";

    for (int i = 0; i < number_descriptors; i++) {
      const Accumulator<double>& acci = global_acc[i];
      cerr << "Descriptor '" << descriptor[i] << "' between " << acci.minval() << " and "
           << acci.maxval();
      if (acci.n() > 1) {
        cerr << ", ave " << acci.average();
      }
      cerr << '\n';
      report_trend(descriptor[i], average[i], cerr);
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = running_average(argc, argv);

  return rc;
}
