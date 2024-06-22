/*
  Computes the distance matrix for a set of fingerprints
*/

#include <stdlib.h>

#include <iomanip>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "gfp.h"
#include "sparse_collection.h"
#include "tversky.h"

using std::cerr;

static int verbose = 0;

static IW_General_Fingerprint *pool = nullptr;

static int pool_size = 0;

static Accumulator<similarity_type_t> stats;

static IWHistogram histogram;

static Report_Progress report_progress;

static int write_multi_token_names_as_separate_columns = 0;

static int equal_weight_tanimoto = 0;

/*
  When doing the full distance matrix, we don't want to compute everything
  twice. But if the dataset is large, we can suppress use of this
  and just recompute everything
*/

static similarity_type_t *stored_result = nullptr;

static int small_memory_behaviour = 0;

static similarity_type_t zero_distance_value = 0.0;

static int lower_triangular_form = 0;

static int output_distances = 1;

static Tversky tversky;

/*
  Do they want fast I/O or not
*/

static int fast_io = 0;

/*
  When doing C++ I/O, we can have any character as the inter-token separator
*/

static char token_separator = ' ';

/*
  By default the fast I/O option has the numbers nicely lined up. That may not be
  appropriate in all cases
*/

static int squeeze_multiple_blanks_from_fast_io = 0;

// not yet implemented. TODO:ianwatson
static Fraction_as_String fraction_as_string;

/*
  C++ I/O is awful.
*/

static const char *string_result[] = {
    "   0", " .01", " .02", " .03", " .04", " .05", " .06", " .07", " .08", " .09",
    " .10", " .11", " .12", " .13", " .14", " .15", " .16", " .17", " .18", " .19",
    " .20", " .21", " .22", " .23", " .24", " .25", " .26", " .27", " .28", " .29",
    " .30", " .31", " .32", " .33", " .34", " .35", " .36", " .37", " .38", " .39",
    " .40", " .41", " .42", " .43", " .44", " .45", " .46", " .47", " .48", " .49",
    " .50", " .51", " .52", " .53", " .54", " .55", " .56", " .57", " .58", " .59",
    " .60", " .61", " .62", " .63", " .64", " .65", " .66", " .67", " .68", " .69",
    " .70", " .71", " .72", " .73", " .74", " .75", " .76", " .77", " .78", " .79",
    " .80", " .81", " .82", " .83", " .84", " .85", " .86", " .87", " .88", " .89",
    " .90", " .91", " .92", " .93", " .94", " .95", " .96", " .97", " .98", " .99",
    " 1.0"};

/*
  This is pretty awful, but it's fairly fast.
  This array of int's holds the string values above, eg " 0.12", as int's.
*/

static int int_sr[101];

static void
initialise_int_sr() {
  for (int i = 0; i < 101; i++) {
    char *s = (char *)&(int_sr[i]);

    strncpy(s, string_result[i], 4);
  }

  return;
}

/*
  Our output buffer is an int too. When we get a result, just copy an element from
  int_sr[] into the appropriate place in output_buffer;
*/

static int *output_buffer = nullptr;

static int write_header_record = 0;

// clang-format off
static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Computes the distance matrix for a pool of fingerprints\n";
 -h               write a header record
 -l               lower triangular matrix only
 -m               output as similarities rather than distances
 -s <size>        specify the number of molecules (usually not needed)
 -f               fast I/O, repeat to squeeze multiple blanks
 -g <char>        inter-token output separator when doing C++ I/O
 -C <fname>       comparator file: default with one file is to use that as comparitor.
 -j               output multi-token names as separate columns (-C only)
 -e               small memory option, requires recomputation instead of storage
 -F -P -W -Q      standard gfp options, enter '-F help' for info
 -q               use equal weight tanimoto function on composite fingerprints
 -r <rpt>         report progress every <rpt> fingerprints processed
 -v               verbose output
)";
// clang-format on

  exit (rc);
}

// clang-format on

static int
allocate_pool(IW_General_Fingerprint *&pool, int pool_size) {
  assert(pool_size > 0 && nullptr == pool);

  pool = new IW_General_Fingerprint[pool_size];

  if (verbose) {
    cerr << "Pool sized for " << pool_size << " molecules\n";
  }

  assert(nullptr != pool);

  return 1;
}

static int
do_squeeze_multiple_blanks_from_fast_io(const int *output_buffer, int jstop,
                                        IWString_and_File_Descriptor &output) {
  const char *s = reinterpret_cast<const char *>(output_buffer);

  IWString tmp(s, 4 * jstop);

  tmp.compress_blanks();

  output << tmp;

  return output.good();
}

static void
write_identifier(const IWString &id, IWString_and_File_Descriptor &output) {
  if (write_multi_token_names_as_separate_columns) {
    output << id;
  } else {
    append_first_token_of_name(id, output);
  }

  return;
}

/*
  jstop might be controlling triangular output
*/

static void
write_line_of_distance_matrix(const int *output_buffer, int jstop,
                              IWString_and_File_Descriptor &output) {
  if (fast_io) {
    if (squeeze_multiple_blanks_from_fast_io) {
      do_squeeze_multiple_blanks_from_fast_io(output_buffer, jstop, output);
    } else {
      output.write((char *)output_buffer, 4 * jstop);
    }
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return;
}

static similarity_type_t
compute_the_distance(IW_General_Fingerprint &fpi, IW_General_Fingerprint &fpj) {
  similarity_type_t result;

  if (tversky.active()) {
    result = fpi.tversky(fpj, tversky);
  } else if (equal_weight_tanimoto) {
    result = fpi.equal_weight_tanimoto(fpj);
  } else {
    result = fpi.tanimoto(&fpj);
  }

  if (output_distances) {
    result = static_cast<similarity_type_t>(1.0) - result;
  }

  return result;
}

/*
  We have N1 rows each of N2 columns wide
*/

static int
distance_matrix(IW_General_Fingerprint *pool1, int n1, IW_General_Fingerprint *pool2,
                int n2, IWString_and_File_Descriptor &output) {
  // output.resize(40000);

  if (pool1 == pool2 && n1 != n2) {
    cerr << "Size mismatch for identical fingerprint sets, " << n1 << " and " << n2
         << '\n';
    return 0;
  }

  if (small_memory_behaviour) {
    ;
  } else if (pool1 != pool2) {  // no symmetry to exploit
    ;
  } else if (0 == lower_triangular_form) {
    stored_result = new similarity_type_t[n1 * n2];
    assert(nullptr != stored_result);
    set_vector(stored_result, n1 * n2, static_cast<similarity_type_t>(8.0));
  }

  if (fast_io) {
    output_buffer = new int[n2];
    assert(nullptr != output_buffer);
  }

  if (verbose) {
    cerr << "N1 " << n1 << " and N2 " << n2 << '\n';
  }

  similarity_type_t result;  // scope here for efficiency
  int jstop;                 // scope here for efficiency

  for (int i = 0; i < n1; i++) {
    IW_General_Fingerprint &fpi = pool1[i];

    write_identifier(fpi.id(), output);

    if (lower_triangular_form) {
      jstop = i;
    } else {
      jstop = n2;
    }

    assert(jstop <= n2);

    for (int j = 0; j < jstop; j++) {
      if (j > i || lower_triangular_form) {
        result = compute_the_distance(fpi, pool2[j]);

        if (nullptr != stored_result) {
          stored_result[j * n1 + i] = result;
        }

        if (verbose) {
          stats.extra(result);
          histogram.extra(result);
        }
      } else if (j < i) {
        if (nullptr != stored_result) {
          result = stored_result[i * n1 + j];
        } else {
          result = compute_the_distance(fpi, pool2[j]);
        }
      } else if (j == i && fpi.id() != pool2[j].id()) {  // must be using the -C option
        result = compute_the_distance(fpi, pool2[j]);
      } else {
        result = zero_distance_value;
      }

      assert(result >= 0.0 && result <= 1.0);
      assert(int(100.0 * result + 0.499999) >= 0 &&
             int(100.0 * result + 0.499999) <= 100);
      if (fast_io) {
        output_buffer[j] = int_sr[int(100.0 * result + 0.499999)];
      } else {
        output << token_separator << result;
      }
      //    cerr << "Completed i = " << i << " j = " << j << '\n';
    }

    //  output.flush();
    write_line_of_distance_matrix(output_buffer, jstop, output);

    if (report_progress()) {
      cerr << "Processed " << i << " fingerprints\n";
    }
  }

  return output.good();
}

static int
build_pool(const char *fname, IW_General_Fingerprint *&pool, int &fingerprints_in_pool,
           int &pool_size) {
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with("PCN<");

    if (0 == pool_size) {
      cerr << "zero occurrences of '"
           << "PCN<"
           << "' in input\n";
      return 0;
    }

    if (!allocate_pool(pool, pool_size)) {
      return 0;
    }
  }

  return build_pool(input, pool, pool_size, fingerprints_in_pool);
}

static int
all_ids_have_same_token_count(const IW_General_Fingerprint *pool, int n) {
  assert(n > 0);

  int nwords = pool[0].id().nwords();

  for (int i = 1; i < n; i++) {
    if (nwords == pool[i].id().nwords()) {
      continue;
    }

    cerr << "Token count mismatch in compartor set\n";
    cerr << "Found '" << pool[i].id() << "' at i = " << i << " expected " << nwords
         << " tokens\n";
    return 0;
  }

  return nwords;
}

static int
WriteHeaderRecord(const IW_General_Fingerprint* pool, int pool_size,
                  IWString_and_File_Descriptor& output) {
  output << "Name";
  for (int i = 0; i < pool_size; ++i) {
    output << token_separator << pool[i].id();

    output.write_if_buffer_holds_more_than(8192);
  }

  output << '\n';

  return 1;
}


static int
distance_matrix(int argc, char **argv) {
  Command_Line cl(argc, argv, "vs:lmV:fF:P:W:Q:C:ejr:g:qh");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (verbose) {
    histogram.initialise(0.0, 1.0, 0.01);
  }

  if (!initialise_fingerprints(cl, verbose)) {
    cerr << "Cannot initialise fingerprint specifications\n";
    usage(23);
  }

  if (cl.option_present('l')) {
    lower_triangular_form = 1;
    if (verbose) {
      cerr << "Output will be lower triangular form\n";
    }
  }

  if (cl.option_present('m')) {
    output_distances = 0;
    zero_distance_value = 1.0;

    if (verbose) {
      cerr << "Will output similarity values rather than distances\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting object\n";
      return 3;
    }
  }

  if (cl.option_present('g')) {
    IWString g;
    cl.value('g', g);

    if (! char_name_to_char(g)) {
      cerr << "Unrecognised token separator '" << token_separator << "'\n";
      return 1;
    }

    token_separator = g[0];

    if (verbose) {
      cerr << "Inter token output separator set to '" << token_separator << "'\n";
    }
  }

  if (cl.option_present('e')) {
    small_memory_behaviour = 1;

    if (verbose) {
      cerr << "Will NOT allocate large arrays, recompute instead\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The -s option must be followed by a whole number > 2\n";
      usage(5);
    }

    if (!allocate_pool(pool, pool_size)) {
      return 12;
    }
  }

  if (cl.option_present('q')) {
    equal_weight_tanimoto = 1;

    if (verbose) {
      cerr << "Will use equal weight Tanimoto measure\n";
    }
  }

  if (cl.option_present('h')) {
    write_header_record = 1;
    if (verbose) {
      cerr << "Will write a header record\n";
    }
  }

  if (cl.option_present('f')) {
    fast_io = 1;

    initialise_int_sr();

    if (verbose) {
      cerr << "Fast I/O\n";
    }

    if (cl.option_count('f') > 1) {
      squeeze_multiple_blanks_from_fast_io = 1;

      if (verbose) {
        cerr << "Multiple blanks squeezed from fast I/O output\n";
      }
    }
  }

  if (cl.option_present('V')) {
    if (!tversky.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot get Tversky specification (-V option)\n";
      usage(18);
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements()) {
    cerr << "Takes only one command line argument\n";
    usage(3);
  }

  int nfingerprints = 0;

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!build_pool(cl[i], pool, nfingerprints, pool_size)) {
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << nfingerprints << " fingerprints\n";
  }

  if (0 == nfingerprints) {
    cerr << "No fingerprints read\n";
    return 32;
  }

  set_default_iwstring_float_concatenation_precision(4);

  IWString_and_File_Descriptor output(1);
  
  if (write_header_record) {
    WriteHeaderRecord(pool, pool_size, output);
  }

  IW_General_Fingerprint *pool2 = nullptr;
  int ncomparators = 0;

  if (cl.option_present('C')) {
    if (lower_triangular_form) {
      cerr << "When using a comparator set (-C), cannot have triangular output (-l)\n";
      return 3;
    }

    const char *fname = cl.option_value('C');

    int pool2_size = 0;

    if (!build_pool(fname, pool2, ncomparators, pool2_size)) {
      cerr << "Cannot build comparator pool '" << fname << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Read " << ncomparators << " comparator fingerprints from '" << fname
           << "'\n";
    }

    if (0 == ncomparators) {
      cerr << "No comparator fingerprints read\n";
      return 3;
    }

    if (cl.option_present('j')) {
      int nwords = all_ids_have_same_token_count(pool2, ncomparators);

      if (0 == nwords) {
        cerr << "Request multi column output but names do not all have same token count\n";
        return 2;
      }

      write_multi_token_names_as_separate_columns = 1;

      if (verbose) {
        cerr << "Will write multi-column names as separate columns\n";
      }

      for (int i = 0; i < nwords; i++) {
        if (i > 0) {
          output << ' ';
        }
        output << "ID" << i;
      }
    } else {
      output << "ID";
    }

    for (int i = 0; i < pool_size; i++) {
      output << ' ';
      append_first_token_of_name(pool[i].id(), output);
      output.write_if_buffer_holds_more_than(32768);
    }

    output << '\n';
  }

  fraction_as_string.initialise(0.0f, 1.0f, 100);
  fraction_as_string.set_include_leading_space(token_separator);

  if (0 == ncomparators) {
    (void)distance_matrix(pool, nfingerprints, pool, pool_size, output);
  } else {
    (void)distance_matrix(pool2, ncomparators, pool, pool_size, output);
  }

  if (verbose) {
    cerr << "Computed " << stats.n() << " values, between " << stats.minval() << " and "
         << stats.maxval() << '\n';
    if (stats.n() > 1) {
      cerr << "Average " << stats.average() << " variance " << stats.variance() << '\n';
    }
    histogram.debug_print(cerr);
  }

  delete[] pool;
  if (nullptr != pool2) {
    delete[] pool2;
  }

  if (nullptr != stored_result) {
    delete[] stored_result;
  }

  return 0;
}

int
main(int argc, char **argv) {
  int rc = distance_matrix(argc, argv);

  return rc;
}
