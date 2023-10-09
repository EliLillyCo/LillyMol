#include <stdlib.h>

#include <fstream>
#include <iostream>

#include "assert.h"

using std::cerr;
using std::cout;
using std::endl;

static int verbose = 0;
const char *prog_name = nullptr;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "tokens_in_quoted_string.h"

extending_resizable_array<int> token_counts;

static Accumulator_Int<int> record_length;

static int records_to_process = 0;

/*
  Sometimes it is the end of the file which is of interest
*/

static std::streampos offset_from_eof = 0;

static int abort_on_different_token_count = 0;

static int report_different_token_count = 0;

static int report_different_from_first_record = 0;

static char token_delimiter = ' ';

static char output_token_delimiter = ' ';

static int consecutive_delimiters_mean_empty_words_between = 0;

static char record_delimiter = '\n';

/*
  When reporting a record different, do we just report the record
  number of do we also write the record
*/

static int write_reported_records_to_cerr = 0;

static IWString new_fname_unreported;
static IWString new_fname_reported;

static std::ofstream new_fstream_unreported;
static std::ofstream new_fstream_reported;

/*
  When dealing with multiple files, it is sometimes convenient to get
  a count of the number of non-rectangular files
*/

static int report_non_rectangular_files = 0;

static int non_rectangular_files_found = 0;

static int data_is_quoted_tokens = 0;

static std::ofstream stream_for_count;

static int token_length_n = 0;
static Accumulator_Int<int> *acc_token_length = nullptr;
static int *zero_length = nullptr;

static int contains_backslashed_separators = 0;

/*
  When writing different column counts to different files, we need a
  means of specifying that
*/

class Stream_for_Columns : public std::ofstream {
 private:
  const int _tokens;

 public:
  Stream_for_Columns(int, const char *);

  int ok() const;
  int debug_print(std::ostream &) const;

  int tokens() const {
    return _tokens;
  }
};

Stream_for_Columns::Stream_for_Columns(int nc, const char *fname) : _tokens(nc) {
  open(fname, std::ios::out);

  if (!good()) {
    cerr << "Cannot open '" << fname << "' for " << _tokens << " stream\n";
  } else if (verbose) {
    cerr << "Stream '" << fname << "' initialised for " << _tokens << " tokens\n";
  }

  return;
}

static resizable_array_p<Stream_for_Columns> cstreams;

static int
write_translate_to_output_delimiter(const const_IWSubstring &buffer,
                                    Stream_for_Columns &s) {
  const int n = buffer.length();

  if (0 == n) {
    s << '\n';
    return 1;
  }

  IWString tmp(buffer);

  if (!data_is_quoted_tokens) {
    tmp.gsub(token_delimiter, output_token_delimiter);
    s << tmp << '\n';

    return 1;
  }

  // The more complex case of quoted text

  const char dquote = '"';

  int in_quote = tmp[0] == dquote;

  for (int i = 1; i < n; ++i) {
    if (dquote == tmp[i]) {
      in_quote = !in_quote;
    } else if (in_quote) {
      ;
    } else if (token_delimiter == tmp[i]) {
      tmp[i] = output_token_delimiter;
    }
  }

  s << tmp << '\n';

  return s.good();
}

static int
process_to_streams(const const_IWSubstring &buffer, int ntokens) {
  int ns = cstreams.number_elements();
  for (int i = 0; i < ns; i++) {
    Stream_for_Columns *s = cstreams[i];
    if (ntokens == s->tokens()) {
      if (token_delimiter == output_token_delimiter) {
        *s << buffer << endl;
      } else {
        write_translate_to_output_delimiter(buffer, *s);
      }
      return 1;
    }
  }

  return 0;
}

template class resizable_array_p<Stream_for_Columns>;
template class resizable_array_base<Stream_for_Columns *>;

static void
usage(int rc, int showDelimiterOptions = 0) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  // cerr << "sizeof (off_t) " << sizeof (off_t) << " sizeof (streampos) " << sizeof
  // (streampos) << endl;
  cerr << "Usage: " << prog_name << " <options> <fname1> <fname2> ...\n";
  cerr << "  -n <number>       process only the first <number> records from each file\n";
  cerr << "  -e <offset>       start processing <offset> bytes before EOF\n";
  cerr << "  -r                report different token count from previous record\n";
  cerr << "  -f                report different token count from first record\n";
  cerr << "  -R <name>         write reported   records to <name>\n";
  cerr << "  -U <name>         write unreported records to <name>\n";
  cerr << "  -t <count>        write records with <count> tokens to -T <fname>\n";
  cerr << "  -T <fname>        file name for corresponding -t\n";
  cerr << "  -o <char>         when writing -T files, use <char> as output delimiter\n";
  cerr << "  -C <fname>        write a descriptor file of 'line.number tcount'\n";
  cerr << "  -b                abort on different token count from previous record\n";
  cerr << "  -w                write reported records to stderr\n";
  if (!showDelimiterOptions) {
    cerr << "  -d <delim>        delimiter is any character, or a token word(e.g. "
            "'tab'). see '-d help' for a list (default ' ')\n";
  }
  cerr << "  -c <delim>        record delimiter (default \\n, use ^m for DOS)\n";
  cerr << "  -u                consecutive delimiters mean empty words in between\n";
  cerr << "  -g                report non rectangular files\n";
  cerr << "  -q                data is quoted tokens\n";
  cerr << "  -k                data may contain backslashed token separators\n";
  cerr << "  -v                verbose output\n";
  if (showDelimiterOptions) {
    cerr << "  -d <delim>        delimiter is any character, or one of these token "
            "words:\n";
    IWString optionVal = "-d";
    char_name_to_char_usage("-d");
  }
  exit(rc);
}

/*
  An R file might have one fewer token on the first record than all other records
*/

static int
header_record_has_only_different_token_count(
    const extending_resizable_array<int> &token_counts,
    const int tokens_on_first_record) {
  resizable_array<int> non_zero_token_counts;

  int nt = token_counts.number_elements();

  for (int i = 1; i < nt; i++) {
    if (0 == token_counts[i]) {
      continue;
    }

    non_zero_token_counts.add(i);

    if (1 == non_zero_token_counts.number_elements())  // first non zero token count
    {
      if (1 != token_counts[i]) {  // if an R file, there must be a single header record
                                   // with that token count
        return 0;
      }

      if (tokens_on_first_record != i) {
        return 0;
      }
    } else if (non_zero_token_counts.number_elements() > 2) {
      return 0;
    }
  }

  return 2 == non_zero_token_counts.number_elements();
}

/*
  No provision for quoted tokens, or consecutive delimiters meaning anything
*/

static int
tokens_in_presence_of_backslashed_separator(const const_IWSubstring &buffer,
                                            const char token_delimiter) {
  int rc = 0;

  const int n = buffer.length();

  int previous_was_backslash = 0;

  for (int i = 0; i < n; ++i) {
    const char c = buffer[i];

    if (token_delimiter == c) {
      if (!previous_was_backslash) {
        rc++;
      }
    } else if ('\\' == c && !previous_was_backslash) {
      previous_was_backslash = true;
    } else {
      previous_was_backslash = false;
    }
  }

  return rc;
}

static void
determine_per_token_lengths_def(const const_IWSubstring &buffer,
                                const char token_delimiter) {
  int col = 0;
  const_IWSubstring token;

  for (int i = 0; buffer.nextword(token, i, token_delimiter); ++col) {
    if (col >= token_length_n) {
      break;
    }

    acc_token_length[col].extra(token.length());

    if (0 == token.length()) {
      zero_length[col]++;
    }
  }

  return;
}

static void
determine_per_token_lengths_cdewb(const const_IWSubstring &buffer,
                                  const char token_delimiter) {
  const_IWSubstring token;
  int col = 0;

  for (int i = 0; buffer.nextword_single_delimiter(token, i, token_delimiter); ++col) {
    if (col >= token_length_n) {
      break;
    }

    acc_token_length[col].extra(token.length());

    if (0 == token.length()) {
      zero_length[col]++;
    }
  }

  return;
}

static void
determine_per_token_lengths(const const_IWSubstring &buffer, const char token_delimiter) {
  if (consecutive_delimiters_mean_empty_words_between) {
    determine_per_token_lengths_cdewb(buffer, token_delimiter);
  }
  // else if (data_is_quoted_tokens)
  //   determine_per_token_lengths_quoted(buffer, token_delimiter);
  else {
    determine_per_token_lengths_def(buffer, token_delimiter);
  }

  return;
}

static int
tcount(iwstring_data_source &input) {
  input.set_dos(1);
  token_counts.resize_keep_storage(0);
  const_IWSubstring buffer;

  int prev_token_count = -1;

  int tokens_on_first_record = 0;

  int rectangular = 1;

  while (input.next_record(buffer)) {
    int tokens;
    if (consecutive_delimiters_mean_empty_words_between) {
      tokens = buffer.nwords_single_delimiter(token_delimiter);
    } else if (data_is_quoted_tokens) {
      tokens = tokens_in_quoted_string(buffer, token_delimiter, input.lines_read());
    } else if (contains_backslashed_separators) {
      tokens = tokens_in_presence_of_backslashed_separator(buffer, token_delimiter);
    } else {
      tokens = buffer.nwords(token_delimiter);
    }

    if (verbose > 2) {
      cerr << "Record " << input.lines_read() << " had " << tokens << " words and "
           << buffer.length() << " characters\n";
    }

    process_to_streams(buffer, tokens);

    if (token_length_n > 0) {
      determine_per_token_lengths(buffer, token_delimiter);
    }

    if (stream_for_count.is_open()) {
      stream_for_count << input.lines_read() << ' ' << tokens << '\n';
    }

    token_counts[tokens]++;

    if (prev_token_count < 0) {
      tokens_on_first_record = tokens;
    }

    int reported = 0;

    if (prev_token_count >= 0 && prev_token_count != tokens) {
      if (abort_on_different_token_count) {
        cerr << "Line " << input.lines_read() << " has " << tokens << " tokens\n";
        cerr << "Previous line had " << prev_token_count << endl;
        non_rectangular_files_found++;
        return 0;
      }

      rectangular = 0;

      if (report_different_token_count) {
        cerr << "Token count mismatch line " << input.lines_read() << ". " << tokens
             << " vs " << prev_token_count << endl;

        reported = 1;
      }
    }

    if (report_different_from_first_record && tokens != tokens_on_first_record) {
      cerr << "Token count mismatch, line " << input.lines_read() << ". First record had "
           << tokens_on_first_record << endl;

      reported++;
    }

    prev_token_count = tokens;

    record_length.extra(buffer.length());

    if (reported) {
      if (write_reported_records_to_cerr) {
        cerr << buffer << endl;
      }

      if (new_fname_reported.length()) {
        new_fstream_reported << buffer << endl;
      }
    } else if (new_fname_unreported.length()) {
      new_fstream_unreported << buffer << endl;
    }

    if (records_to_process > 0 && input.lines_read() > records_to_process) {
      break;
    }
  }

  if (!rectangular) {
    non_rectangular_files_found++;
  }

  int rc = 0;

  int nt = token_counts.number_elements();

  // If the only record that has a different token count is the header record, note that.

  int looks_like_R =
      header_record_has_only_different_token_count(token_counts, tokens_on_first_record);

  for (int i = 0; i < nt; i++) {
    if (0 == token_counts[i]) {
      continue;
    }

    cerr << token_counts[i] << " records had " << i << " tokens";
    if (i == tokens_on_first_record && looks_like_R) {
      cerr << " HEADER";
    }
    cerr << endl;
    rc++;
  }

  return rc - 1;
}

static int
tcount(const char *fname) {
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose) {
    cerr << fname << endl;
  }

  input.set_record_delimiter(record_delimiter);

  if (offset_from_eof) {
    const auto fs = input.file_size();

    if (fs < offset_from_eof) {
      cerr << "Requested offset from EOF " << offset_from_eof << " larger than file "
           << fs << endl;
      return 0;
    }

    const auto tmp = fs - offset_from_eof;
    if (!input.seekg(tmp)) {
      cerr << "Yipes, cannot seek to offset " << tmp << "\n";
      return 0;
    }

    const_IWSubstring notused;
    if (!input.next_record(notused)) {
      cerr << "Yipes, could not read record from offset " << tmp << endl;
      return 0;
    }
  }

  int rc = tcount(input);

  cerr << input.lines_read() << " records read from '" << fname << "'";
  if (verbose) {
    cerr << " records between " << record_length.minval() << " and "
         << record_length.maxval();
    if (record_length.n() > 1) {
      cerr << " ave " << record_length.average();
    }
  }
  cerr << endl;

  return rc;
}

static int
tcount(int argc, char **argv) {
  Command_Line cl(argc, argv, "c:d:wt:T:U:R:fbrvn:e:guqC:o:L:k");

  if (cl.unrecognised_options_encountered()) {
    usage(2);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('d')) {
    IWString tmp;
    cl.value('d', tmp);
    tmp.to_lowercase();
    if (tmp == "help") {
      usage(1, 1);
    }
    if (!char_name_to_char(tmp)) {
      cerr << "Sorry, the delimiter character can be one character only, or one of a "
              "list of tokens.  see -d help for details\n";
      usage(13);
    }

    token_delimiter = tmp[0];
    //  if (' ' != token_delimiter)
    //    consecutive_delimiters_mean_empty_words_between = 1;
  }

  if (cl.option_present('c')) {
    const_IWSubstring c = cl.string_value('c');

    if (c.length() > 1) {
      cerr << "Sorry, the record delimiter can be 1 character at most\n";
      usage(15);
    }

    record_delimiter = c[0];
    if (verbose) {
      cerr << "Record delimiter is '" << record_delimiter << "'\n";
    }
  }

  if (cl.option_present('u')) {
    consecutive_delimiters_mean_empty_words_between = 1;
    if (verbose) {
      cerr << "Consecutive delimiters mean empty words in between adjacent delimiters\n";
    }
  }
  if (cl.option_present('q')) {
    data_is_quoted_tokens = 1;
    if (verbose) {
      cerr << "Data is quoted tokens\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', records_to_process) || records_to_process <= 0) {
      cerr << "The -n option requires a positive whole number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will read only the first " << records_to_process
           << " records from each file\n";
    }
  }

  if (cl.option_present('L')) {
    if (!cl.value('L', token_length_n) || token_length_n < 1) {
      cerr << "The number of columns for which you want length stats (-L) must be a "
              "whole +ve number\n";
      usage(1);
    }

    acc_token_length = new Accumulator_Int<int>[token_length_n];
    zero_length = new int[token_length_n];
    std::fill_n(zero_length, token_length_n, 0);
  }

  if (cl.option_present('e')) {
    unsigned long tmp;
    if (!cl.value('e', tmp) || tmp < 3) {
      cerr << "Must specify a value offset from eof\n";
      usage(8);
    }

    offset_from_eof = tmp;

    if (verbose) {
      cerr << "Will start processing " << offset_from_eof << " bytes from eof\n";
    }
  }

  if (cl.option_present('b')) {
    abort_on_different_token_count = 1;
    if (verbose) {
      cerr << "Will stop processing on a differing token count\n";
    }
  }

  if (cl.option_present('r')) {
    report_different_token_count = 1;
    if (verbose) {
      cerr << "Will report all instances where token count is different from previous "
              "record\n";
    }
  }

  if (cl.option_present('f')) {
    report_different_from_first_record = 1;

    if (verbose) {
      cerr << "Will report instances where token count is different from first record\n";
    }
  }

  if (cl.option_present('g')) {
    report_non_rectangular_files = 1;
    if (verbose) {
      cerr << "Will report non-rectangular files\n";
    }
  }

  if (cl.option_present('k')) {
    contains_backslashed_separators = 1;

    if (verbose) {
      cerr << "Input may contain backslashed token separators\n";
    }
  }

  if (cl.option_present('U')) {
    if (cl.option_present('b')) {
      cerr << "The -U and -b options are inconsistent\n";
      usage(18);
    }

    if (!cl.option_present('r') && !cl.option_present('f')) {
      cerr << "The -U option only makes sense with either the -r or -f options\n";
      usage(73);
    }

    cl.value('U', new_fname_unreported);

    new_fstream_unreported.open(new_fname_unreported.null_terminated_chars(),
                                std::ios::out);
    if (!new_fstream_unreported.good()) {
      cerr << "Cannot open '" << new_fname_unreported << "' for output\n";
      return 8;
    }

    if (verbose) {
      cerr << "Unreported records will be written to '" << new_fname_unreported << "'\n";
    }
  }

  if (cl.option_present('R')) {
    if (cl.option_present('b')) {
      cerr << "The -R and -b options are inconsistent\n";
      usage(13);
    }

    if (!cl.option_present('r') && !cl.option_present('f')) {
      cerr << "The -R option only makes sense with either the -r or -f options\n";
      usage(73);
    }

    cl.value('R', new_fname_reported);

    new_fstream_reported.open(new_fname_reported.null_terminated_chars(), std::ios::out);
    if (!new_fstream_reported.good()) {
      cerr << "Cannot open reported stream '" << new_fname_reported << "'\n";
      return 9;
    }

    if (verbose) {
      cerr << "Reported records will be written to '" << new_fname_reported << "'\n";
    }
  }

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    if (!char_name_to_char(o)) {
      cerr << "Unrecognised output token directive '" << o << "'\n";
      usage(1);
    }
    output_token_delimiter = o[0];
  }

  int nt = cl.option_count('t');
  if (nt != cl.option_count('T')) {
    cerr << "There must be the same number of -t and -T options\n";
    usage(34);
  }

  for (int i = 0; i < nt; i++) {
    int t;
    if (!cl.value('t', t, i) || t < 0) {
      cerr << "The -t option requires a whole positive number\n";
      usage(5);
    }

    IWString fname;
    (void)cl.value('T', fname, i);

    Stream_for_Columns *s = new Stream_for_Columns(t, fname.null_terminated_chars());
    if (!s->good()) {
      return 9;
    }

    cstreams.add(s);
  }

  if (cl.option_present('w')) {
    write_reported_records_to_cerr = 1;
    if (verbose) {
      cerr << "Will echo reported records to stderr\n";
    }
  }

  if (0 == cl.number_elements()) {
    usage(1);
  }

  if (cl.option_present('C')) {
    const char *c = cl.option_value('C');

    stream_for_count.open(c, std::ios::out);

    if (!stream_for_count.good()) {
      cerr << "Cannot open stream for tcount '" << c << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Data on column counts written to '" << c << "'\n";
    }

    stream_for_count << "Line Tcount\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    rc += tcount(cl[i]);
  }

  if (report_non_rectangular_files) {
    cerr << non_rectangular_files_found << " non rectangular files found\n";
  }

  return 0;
}

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = tcount(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
