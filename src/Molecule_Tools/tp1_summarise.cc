/*
  Examine the output of tp1_pipe.sh and create a formatted report
*/

#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString bad_file_stem ("bad");

static IWString bad_file_stem_array ("BQTP");

static int molecules_written = 0;

static int rejected_molecules = 0;

static int demerited_molecules = 0;

static int include_reason = 0;

static int include_zero_demerit_molecules = 0;

static int include_header_record = 0;

static IWString prepend_string;

static int prepend_d = 1;

static IW_STL_Hash_Map_int reason_bad;
static IW_STL_Hash_Map_int reason_survivor;

static int accumulate_reasons = 0;

static extending_resizable_array<int> demerits_per_molecule;

static char output_separator = ' ';

static int latex_table = 0;

static int asciidoc_table = 0;

static int gsub_spaces_in_reason_to_underscore = 0;

static int bad0_demerit = 200;
static int bad12_demerit = 100;

static int produce_demerit_scaling_file = 0;

static int suppress_normal_output = 0;

static int process_rejection_files = 1;

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
  cerr << " -r             include the rejection reason with the output\n";
  cerr << " -z             include zero demerit molecules in the output\n";
  cerr << " -h             include a header record\n";
  cerr << " -D             exclude the D prefix\n";
  cerr << " -B <stem>      bad file stem (default 'bad')\n";
  cerr << " -s <string>    separator between output fields\n";
  cerr << " -u             change spaces in reason fields to underscores\n";
  cerr << " -t             give report of which reasons hit\n";
  cerr << " -T <fname>     write report on rejection reasons to <fname>\n";
  cerr << " -b             in the -T file, produce two column format rejections | demerits\n";
  cerr << " -X             produce table in LaTex format\n";
  cerr << " -A             produce table in AsciiDoc format\n";
  cerr << " -m <n>         in the table file, discard any reason with <n> or fewer examples\n";
  cerr << " -c             produce a demerit based scale factor file\n";
  cerr << " -f <n>         numeric demerit value for rejections (default 100)\n";
  cerr << " -k             only process the survivors file (do not process bad0, bad1...)\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
write_demerit_value(const const_IWSubstring & id,
                    const int demerit,
                    IWString_and_File_Descriptor & output)
{
  molecules_written++;

  if (demerit >= 100)
    rejected_molecules++;
  else if (demerit > 0)
    demerited_molecules++;

  if (suppress_normal_output)
    return 1;

  output << id << output_separator;

  if (prepend_string.length())
    output << prepend_string;

  if (prepend_d)
    output << 'D';

  if (produce_demerit_scaling_file)
  {
    if (demerit >= 100)
      output << '0';
    else
      output << static_cast<float>(100 - demerit) * 0.01f;
  }
  else
    output << demerit;

  output.write_if_buffer_holds_more_than(32768);

  return output.good ();
}

re2::RE2 d_parentheses("^D\\([0-9]+\\)$");
//static IW_Regular_Expression d_parentheses ("^D\\([0-9]+\\)$");

static int
process_from_iwdemerit(const const_IWSubstring & buffer,
                       int must_have_demerit,
                       IWString_and_File_Descriptor & output)
{
  int i = 0;

  const_IWSubstring token;

  buffer.nextword (token, i);

  const_IWSubstring id;

  buffer.nextword (id, i);

//cerr << "Processing '" << id << "'\n";

  int previous_was_colon = 0;

  const_IWSubstring dmrt;

  while (buffer.nextword (token, i))
  {
    if (':' == token)
    {
      previous_was_colon = 1;
      continue;
    }

    if (! previous_was_colon)
      continue;

    if (! iwre2::RE2PartialMatch(token, d_parentheses)) 
//  if (! d_parentheses.matches (token))
    {
      previous_was_colon = 0;
      continue;
    }

    dmrt = token;
    dmrt.remove_leading_chars (2);
    dmrt.chop ();
    break;
  }

//cerr << " Demerit for '" << id << " is '" << dmrt << "'\n";

  if (dmrt.length ())
    ;
  else if (must_have_demerit)
  {
    cerr << "NO demerit value found for '" << id << "'\n";
    return 0;
  }
  else
  {
    if (include_zero_demerit_molecules)
    {
      write_demerit_value (id, 0, output);
      if (! suppress_normal_output)
        output << '\n';
    }

    return output.good ();
  }

  int d;
  (void) dmrt.numeric_value (d);

  if (d > 100)
    d = 100;

  write_demerit_value (id, d, output);

  if (include_reason)
  {
    int demerit_reasons_this_molecule = 0;

    static IWString myreason;

    while (buffer.nextword (myreason, i, ':'))
    {
      myreason.strip_leading_blanks();

      if (! suppress_normal_output)
        output << output_separator << myreason;

      if (! accumulate_reasons)
        ;
      else if (d >= 100)
        reason_bad[myreason]++;
      else
        reason_survivor[myreason]++;

      demerit_reasons_this_molecule++;
    }

    demerits_per_molecule[demerit_reasons_this_molecule]++;
  }

  if (! suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(23768);

  return output.good ();
}

static int
process_from_iwdemerit (iwstring_data_source & input,
                        int must_have_demerit,
                        IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! process_from_iwdemerit (buffer, must_have_demerit, output))
    {
      cerr << "Invalid from iwdemerit record, line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return output.good ();
}

//static IW_Regular_Expression open_paren ("^\\([0-9]+$");
static re2::RE2 open_paren("^\\([0-9]+$");

static int
process_bad12 (const const_IWSubstring & buffer,
               IW_STL_Hash_Map_int & reason,
               IWString_and_File_Descriptor & output)
{
  int i = 0;
  const_IWSubstring token;

  buffer.nextword (token, i);

  const_IWSubstring id;
  buffer.nextword (id, i);

  int got_open_paren = 0;

  while (buffer.nextword (token, i))
  {
    if (! iwre2::RE2PartialMatch(token, open_paren))
      continue;

    got_open_paren = 1;
    break;
  }

  if (! got_open_paren)
  {
    cerr << "Never found a '(\\d+' token\n";
    return 0;
  }

  buffer.nextword (token, i);

  if ("matches" != token)
  {
    cerr << "Expected 'matches' but got '" << token << "'\n";
    return 0;
  }

  buffer.nextword (token, i);

  if ("to" != token)
  {
    cerr << "Expected 'to' but got '" << token << "'\n";
    return 0;
  }

  write_demerit_value (id, bad12_demerit, output);

  if (include_reason)
  {
    static IWString myreason;

    myreason.resize_keep_storage (0);

    while (buffer.nextword (token, i))
    {
      if (token.starts_with ('\''))
        token.remove_leading_chars (1);

      if (token.ends_with ("')"))
        token.chop (2);

      myreason.append_with_spacer (token);
    }

    if (! suppress_normal_output)
      output << output_separator << myreason;

    if (accumulate_reasons)
      reason[myreason]++;
  }

  if (! suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good ();
}

static int
process_bad12 (iwstring_data_source & input,
               IW_STL_Hash_Map_int & reason,
               IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! process_bad12 (buffer, reason, output))
    {
      cerr << "Invalid bad12 record '" << buffer << "'\n";
      return 0;
    }
  }

  return output.good ();
}

static int
process_bad12 (const IWString & fname,
               IW_STL_Hash_Map_int & reason,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return process_bad12 (input, reason, output);
}

static int
process_bad1 (const IWString & bad_stem,
              IW_STL_Hash_Map_int & reason,
              IWString_and_File_Descriptor & output)
{
  IWString fname;

  fname = bad_stem;
  fname << "1.smi";

  return process_bad12 (fname, reason, output);
}

static int
process_bad2 (const IWString & bad_stem,
              IW_STL_Hash_Map_int & reason,
              IWString_and_File_Descriptor & output)
{
  IWString fname;

  fname = bad_stem;
  fname << "2.smi";

  return process_bad12 (fname, reason, output);
}

static int
process_from_iwdemerit (const char * fname,
                        int must_have_demerit,
                        IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return process_from_iwdemerit (input, must_have_demerit, output);
}

static int
process_from_iwdemerit (const IWString & bad_stem,
                        IWString_and_File_Descriptor & output,
                        int must_have_demerit)
{
  IWString fname;

  fname = bad_stem;
  fname << "3.smi";

  return process_from_iwdemerit (fname.null_terminated_chars (), must_have_demerit, output);
}

static int
process_bad0 (const const_IWSubstring & buffer,
              IW_STL_Hash_Map_int & reason,
              IWString_and_File_Descriptor & output)
{
  int i = 0;
  const_IWSubstring token;

  buffer.nextword (token, i);

  const_IWSubstring id;
  buffer.nextword (id, i);

  int got_tp1 = 0;

  while (buffer.nextword (token, i))
  {
    if ("TP1" != token)
      continue;

    got_tp1 = 1;
    break;
  }

  if (! got_tp1)
  {
    cerr << "No 'TP1' token\n";
    return 0;
  }

  write_demerit_value (id, bad0_demerit, output);

  if (include_reason)
  {
    static IWString myreason;

    myreason.resize_keep_storage (0);

    while (buffer.nextword (token, i))
    {
      myreason.append_with_spacer (token);
    }

    if (gsub_spaces_in_reason_to_underscore)
      myreason.gsub (' ', '_');

    if (! suppress_normal_output)
      output << output_separator << myreason;

    if (accumulate_reasons)
      reason[myreason]++;
  }

  if (! suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good ();
}

static int
process_bad0 (iwstring_data_source & input, 
              IW_STL_Hash_Map_int & reason,
              IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! process_bad0 (buffer, reason, output))
    {
      cerr << "Invalid bad0 record '" << buffer << "'\n";
      return 0;
    }
  }

  return output.good ();
}

static int
process_bad01 (const IWString & fname,
               IW_STL_Hash_Map_int & reason,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return process_bad0 (input, reason, output);
}

static int
process_bad0 (const IWString & bad_stem,
              IW_STL_Hash_Map_int & reason,
              IWString_and_File_Descriptor & output)
{
  IWString fname;

  fname = bad_stem;
  fname << "0.smi";

  return process_bad01 (fname, reason, output);
}

static int
all_files_present (const char * okfile,
                   const IWString & bad_stem)
{
  for (int i = 0; i < 3; i++)
  {
    IWString fname;
    fname << bad_stem << i << ".smi";

    if (verbose > 1)
      cerr << "Checking '" << fname << "', result " << dash_f (fname.null_terminated_chars ()) << endl;

    if (! dash_f (fname.null_terminated_chars ()))
    {
      if (verbose)
        cerr << "File '" << fname << "' not present\n";
      return 0;
    }
  }

  return dash_f (okfile);
}

static int
process_single_run (const char * okfile,
                    const IWString & bad_stem,
                    int ndx,
                    IWString_and_File_Descriptor & output)
{
  if (process_rejection_files)
  {
    process_bad0 (bad_stem, reason_bad, output);
    process_bad1 (bad_stem, reason_bad, output);
    process_bad2 (bad_stem, reason_bad, output);
  }
  process_from_iwdemerit (bad_stem, output, 1);
  process_from_iwdemerit (okfile, 0, output);

  return output.good ();
}

class Reason_and_Count
{
  private:
    const IWString _reason;
    const int _count;

  public:
    Reason_and_Count (const IWString & r, int c) : _reason(r), _count(c) {}

    const IWString & reason () const { return _reason;}
    int count () const { return _count;}

    template <typename T> int latex_table (T &) const;
    template <typename T> int asciidoc_table (T & os) const;
};

template <typename T>
T &
operator << (T & os, const Reason_and_Count & rc)
{
  os << rc.count() << ' ' << rc.reason();

  return os;
}

template <typename T>
int
Reason_and_Count::latex_table (T & os) const
{
  os << _count << " & ";

  if (_reason.contains('_'))
  {
    IWString tmp(_reason);
    tmp.gsub("_", "\\_");
    os << tmp;
  }
  else
    os << _reason;

  os << "\\\\\n";

  return 1;
}

template <typename T>
int
Reason_and_Count::asciidoc_table (T & os) const
{
  os << '|' << _count << " | " << _reason;

  return 1;
}

class Reason_and_Count_Comparator
{
  private:
  public:
    int operator () (const Reason_and_Count *, const Reason_and_Count *) const;
};

int
Reason_and_Count_Comparator::operator () (const Reason_and_Count * rc1, 
                                          const Reason_and_Count * rc2) const
{
  if (rc1->count() < rc2->count())
    return 1;
  if (rc1->count() > rc2->count())
    return -1;

  return 0;
}

static int
last_item_meeting_support_requirement (const resizable_array_p<Reason_and_Count> & r,
                                       const int min_reasons_needed_for_output)
{
  for (int i = r.number_elements() - 1; i >= 0; --i)
  {
    if (r[i]->count() >= min_reasons_needed_for_output)
      return i;
  }

  return -1;
}

static void
sorted_list_of_reasons (const IW_STL_Hash_Map_int & reason,
                        resizable_array_p<Reason_and_Count> & r)
{
  int n = reason.size();

  r.resize(n);

  for (IW_STL_Hash_Map_int::const_iterator i = reason.begin(); i != reason.end(); ++i)
  {
    Reason_and_Count * t = new Reason_and_Count( (*i).first, (*i).second);
    r.add(t);
  }

  Reason_and_Count_Comparator rcc;

  r.iwqsort(rcc);

  return;
}

template <typename T>
void
write_rejection_reason (const resizable_array_p<Reason_and_Count> & r,
                        const int ndx,
                        T & output)
{
  if (ndx < r.number_elements())
  {
    if (latex_table)
      r[ndx]->latex_table(output);
    else if (asciidoc_table)
      r[ndx]->asciidoc_table(output);
    else
      output << *(r[ndx]);
  }
  else
  {
    if (latex_table)
      output << " & &";
    else if (asciidoc_table)
      output << "|.|.";
  }

  return;
}

/*
  WRite a multi column file containing rejection reasons, and demerit reasons
*/

template <typename T>
int
write_reasons (const IW_STL_Hash_Map_int & reason_bad,
               const IW_STL_Hash_Map_int & reason_survivor,
               const int min_reasons_needed_for_output,
               T & output)
{
  resizable_array_p<Reason_and_Count> rej, dem;
  sorted_list_of_reasons(reason_bad, rej);
  sorted_list_of_reasons(reason_survivor, dem);

  int dstop = last_item_meeting_support_requirement(dem, min_reasons_needed_for_output);
  int rstop = last_item_meeting_support_requirement(rej, min_reasons_needed_for_output);

  if (dstop < 0 && rstop < 0)
    return 1;

  if (asciidoc_table)
  {
    output << "[width=\"30%\", options=\"header\"]\n";
    output << "|===============\n";
    output << "|N |Rejected |N |Demerits\n";
  }
  else if (latex_table)
  {
    output << "\\begin{center}\n";
    output << "\\begin{tabular}{r c r c}\n";
    output << "Molecules & Reason Rej & Molecules & Reason Demerit \\\\\n";
    output << "\\hline\n";
  }

  int istop = dstop;
  if (istop < rstop)
    istop = rstop;

  for (int i = 0; i < istop; ++i)
  {
    write_rejection_reason(rej, i, output);
    write_rejection_reason(dem, i, output);
    output << '\n';
  }

  if (asciidoc_table)
    output << "|===============\n";
  else if (latex_table)
  {
    output << "\\hline\n";
    output << "\\end{tabular}\n";
    output << "\\end{center}\n";
  }

  return 1;
}

template <typename T>
int
write_reasons (const IW_STL_Hash_Map_int & reason,
               const int min_reasons_needed_for_output,
               T & output)
{
  resizable_array_p<Reason_and_Count> r;

  sorted_list_of_reasons(reason, r);

  if (asciidoc_table)
  {
    output << "[width=\"30%\"]\n";
    output << "|===============\n";
  }
  else if (latex_table)
  {
    output << "\\begin{center}\n";
    output << "\\begin{tabular}{r l}\n";
    output << "Molecules & Reason \\\\\n";
    output << "\\hline\n";
  }

  for (int i = 0; i < r.number_elements(); i++)
  {
    const Reason_and_Count * ri = r[i];

    if (ri->count() < min_reasons_needed_for_output)
      break;

    if (latex_table)
      ri->latex_table(output);
    else if (asciidoc_table)
      ri->asciidoc_table(output);
    else
      output << (*ri);

    output << '\n';
  }

  if (asciidoc_table)
  {
    output << "|===============\n";
  }
  else if (latex_table)
  {
    output << "\\hline\n";
    output << "\\end{tabular}\n";
    output << "\\end{center}\n";
  }

  return 1;
}

static int
tp1_summarise (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vB:S:rhzDs:tT:uj:f:XAP:cnkbm:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('B'))
  {
    cl.value ('B', bad_file_stem);

    if (verbose)
      cerr << "Bad files have stem '" << bad_file_stem << "'\n";

    bad_file_stem_array = bad_file_stem;
  }

  if (cl.option_present ('r'))
  {
    include_reason = 1;

    if (verbose)
      cerr << "The reason for the rejection will be included\n";

    if (cl.option_present ('u'))
    {
      gsub_spaces_in_reason_to_underscore = 1;

      if (verbose)
        cerr << "Spaces in reason converted to underscore\n";
    }
  }

  if (cl.option_present ('h'))
  {
    include_header_record = 1;

    if (verbose)
      cerr << "Will write a header record\n";
  }

  if (cl.option_present ('z'))
  {
    include_zero_demerit_molecules = 1;

    if (verbose)
      cerr << "Will include zero demerit molecules\n";
  }

  if (cl.option_present ('D'))
  {
    prepend_d = 0;

    if (verbose)
      cerr << "No D prefix\n";
  }

  if (cl.option_present('P'))
  {
    cl.value('P', prepend_string);

    if (verbose)
      cerr << "Will prepend '" << prepend_string << "' to all output\n";
  }

  if (cl.option_present ('t'))
  {
    accumulate_reasons = 1;
    include_reason = 1;

    if (verbose)
      cerr << "Will accumulate reasons for rejections\n";
  }

  const char * fname_for_reason_summary = nullptr;

  if (cl.option_present('T'))
  {
    fname_for_reason_summary = cl.option_value('T');

    if (verbose)
      cerr << "Summary of rejection reasons written to '" << fname_for_reason_summary << "'\n";

    accumulate_reasons = 1;
    include_reason = 1;
  }

  if (cl.option_present('c'))
  {
    produce_demerit_scaling_file = 1;

    if (verbose)
      cerr << "Will produce a file where demerit values have been converted to scaling factors\n";
  }

  if (cl.option_present('n'))
  {
    suppress_normal_output = 1;

    if (verbose)
      cerr << "Will suppress normal output\n";
  }

  if (cl.option_present('k'))
  { 
    process_rejection_files = 0;

    if (verbose)
      cerr << "Will only process survivor files (not bad0, bad1...)\n";
  }

// Initial implementation used -j, but it should be -f to be compatible with
// what is used in iwdemerit

  const_IWSubstring fj;

  if (cl.option_present('j'))
    cl.value('j', fj);
  else if (cl.option_present('f'))
    cl.value('f', fj);

  if (fj.length())
  {
    int j;
    if (! fj.numeric_value(j) || j < 1)
    {
      cerr << "The rejection threshold value (-f) must be a whole +ve number\n";
      usage (4);
    }

    bad0_demerit = j;
    bad12_demerit = j;

    if (verbose)
      cerr << "Demerit value assigned to hard rejections " << j << endl;
  }

  if (cl.option_present ('s'))
  {
    IWString s;
    cl.value ('s', s);

    if (! char_name_to_char(s))
    {
      cerr << "Unrecognised output separator '" << s << "'\n";
      return 1;
    }

    output_separator = s[0];
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  if (include_header_record)
  {
    output << "Name" << output_separator;
    if (produce_demerit_scaling_file)
      output << "scale\n";
    else
      output << "demerit\n";
  }

  int rc = 0;

  cerr << "Processing " << cl.number_elements() << " files\n";

  if (1 == cl.number_elements ())
    process_single_run (cl[0], bad_file_stem, -1, output);
  else
  {
    for (int i = 0; i < cl.number_elements (); i++)
    {
      IWString bstem = bad_file_stem_array;
      bstem << (i + 1) << '_';

      if (! all_files_present (cl[i], bstem))
        break;

      if (verbose)
        cerr << "Processing okfile '" << cl[i] << "', bad stem '" << bstem << "'\n";

      if (! process_single_run (cl[i], bstem, i, output))
      {
        cerr << "Fatal error processing files with stem '" << bstem << "'\n";
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (cl.option_present('X'))
  {
    latex_table = 1;
  }

  if (cl.option_present('A'))
  {
    asciidoc_table = 1;
  }

  if (accumulate_reasons)
  {
    cerr << "Encountered " << (reason_bad.size() + reason_survivor.size()) << " different reasons\n";

    int min_reasons_needed_for_output = 0;

    if (cl.option_present('m'))
    {
      if (! cl.value('m', min_reasons_needed_for_output) || min_reasons_needed_for_output < 0)
      {
        cerr << "The minimum reasons for output (-m) option must be a whole +ve number\n";
        usage(1);
      }

      if (verbose)
        cerr << "Will not write any rejection/demerit reason to the -T file with fewer than " << min_reasons_needed_for_output << endl;
    }

    if (nullptr == fname_for_reason_summary)
    {
      for (IW_STL_Hash_Map_int::const_iterator i = reason_bad.begin (); i != reason_bad.end (); ++i)
      {
        if (i->second > min_reasons_needed_for_output)
          cerr << (*i).second << " occurrences of '" << (*i).first << "'\n";
      }
      for (IW_STL_Hash_Map_int::const_iterator i = reason_survivor.begin (); i != reason_survivor.end (); ++i)
      {
        if (i->second > min_reasons_needed_for_output)
          cerr << (*i).second << " occurrences of '" << (*i).first << "'\n";
      }
    }
    else
    {
      IWString_and_File_Descriptor tfile;
      if (! tfile.open(fname_for_reason_summary))
      {
        cerr << "Cannot open reason summary file '" << fname_for_reason_summary<< "'\n";
        return 5;
      }

      if (cl.option_present('b'))
        write_reasons(reason_bad, reason_survivor, min_reasons_needed_for_output, tfile);
      else
      {
        if (verbose)
          cerr << reason_bad.size() << " reasons associated with rejected molecules, " << reason_survivor.size() << " with survivors\n";
        write_reasons (reason_bad, min_reasons_needed_for_output, tfile);
        write_reasons (reason_survivor, min_reasons_needed_for_output, tfile);
      }
    }
  }

  if (verbose)
  {
    cerr << "Wrote " << molecules_written << " molecules, ";
    cerr << rejected_molecules << " rejected, " << demerited_molecules << " demerited\n";

    for (int i = 0; i < demerits_per_molecule.number_elements (); i++)
    {
      if (demerits_per_molecule[i])
        cerr << demerits_per_molecule[i] << " molecules had " << i << " demerits\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tp1_summarise (argc, argv);

  return rc;
}

// arch-tag: da1723f4-6a29-4479-8a6e-ba7ee05882a6
