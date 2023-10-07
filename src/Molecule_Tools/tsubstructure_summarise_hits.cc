/*
  Scan an output file from tsubstructure and summarise the hits
*/

#include <iostream>
#include <fstream>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwstring/iw_stl_hash_multimap.h"
#include "Foundational/data_source/iwstring_data_source.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int ignore_unrecognised_records = 0;

static int unrecognised_records_encountered = 0;

static int append_reason = 0;

static int first_token_of_multi_token_query_names = 0;

static int create_matches_files = 1;

static IWString suffix(".smi");

static IWString prefix("TSUB_");

static int write_names_of_molecules_matched = 0;

static IWString_and_File_Descriptor stream_for_molecule_and_hits;

static IWString queries_matching_separator(':');

static int include_smiles_in_M_file = 1;

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
  cerr << " -P <prefix>    create files with common prefix (default '" << prefix << "')\n";
  cerr << " -z             ignore unrecognised records\n";
  cerr << " -f             take the first token of multi-token query names\n";
  cerr << " -T <fname>     create table of hit statistics\n";
  cerr << " -t             write names of molecules hit in the -T file\n";
  cerr << " -M <fname>     for each molecule table of queries that hit\n";
  cerr << " -M sep=.       separator for multiple queries matching in -M file (def :)\n";
  cerr << " -M nosmi       exclude smiles from -M file\n";
  cerr << " -M hdr         write a header record in -M file (suppresses smiles)\n";
  cerr << " -n             suppress writing individual files when using -T option\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
strip_to_first_token (IWString & qry_name)
{
  if (qry_name.nwords() > 0)
  {
    qry_name.truncate_at_first(' ');
    return 1;
  }

  int ndx = qry_name.index('	');      // look for tab
  if (ndx > 0)
  {
    qry_name.iwtruncate(ndx);
    return 1;
  }

  return 0;
}

static int
tp_first_pass_summarise_hits (const const_IWSubstring & buffer,
                              IW_STL_Hash_Multimap_IWString & count)
{
  IWString to_store;

  const_IWSubstring token;

  int i = 0;

  while (buffer.nextword(token, i))
  {
    if ("TP1" == token)
      break;

    to_store.append_with_spacer(token);
  }

  IWString zreason;

  while (buffer.nextword(token, i))
  {
    zreason.append_with_spacer(token);
  }

  if (append_reason)
    to_store.append_with_spacer(zreason);

  count.insert(IW_STL_Hash_Multimap_IWString::value_type(zreason, to_store));

//cerr << "storing '" << to_store << "' for '" << zreason << "'\n";

  return count.size();
}

static int
get_reason (const const_IWSubstring & buffer,
            const_IWSubstring token,
            int & i,
            IWString & zreason)
{
  zreason.resize_keep_storage(0);

  while (buffer.nextword(token, i))
  {
    zreason.append_with_spacer(token);

    if (zreason.ends_with(':'))
    {
      zreason.chop();
      return 1;
    }
  }

  return zreason.length();
}

static re2::RE2 d_parentheses("^D\\([0-9]+\\)$");

static int
iwdemerit_summarise_hits (const const_IWSubstring & buffer,
                          IW_STL_Hash_Multimap_IWString & count)
{
  IWString to_store;

  int i = 0;

  const_IWSubstring token;

  while (buffer.nextword(token, i))
  {
    to_store.append_with_spacer(token);

    if (':' == token)
      break;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "End of buffer before D(nn)\n";
    return 0;
  }

  if (! iwre2::RE2PartialMatch(token, d_parentheses))
  {
    cerr << "Not D(nn) form '" << token << "'\n";
    return 0;
  }

  IWString reason;

  while (get_reason(buffer, token, i, reason))
  { 
    if (append_reason)
    {
      IWString tmp(to_store);
      tmp.append_with_spacer(reason);
      count.insert(IW_STL_Hash_Multimap_IWString::value_type(tmp, to_store));
    }
    else
      count.insert(IW_STL_Hash_Multimap_IWString::value_type(reason, to_store));
  }

  return count.size();
}

static re2::RE2 open_paren("^\\([0-9]+$");

static int
number_matches_to(const const_IWSubstring & buffer,
                  int & i,
                  const const_IWSubstring & s)   // the token just read
{
  if (! iwre2::RE2PartialMatch(s, open_paren))
  {
    cerr << "Should be '" << open_paren.pattern() << "', got '" << s << "'\n";
    return 0;
  }

  const_IWSubstring token;

  if (! buffer.nextword(token, i))
  {
    cerr << "Premature termination\n";
    return 0;
  }

  if ("matches" != token)
  {
    cerr << "Should be 'matches', got '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "Premature termination\n";
    return 0;
  }

  if ("to" != token)
  {
    cerr << "Should be 'to', got '" << token << "'\n";
    return 0;
  }

  return 1;
}

static int
get_query_name (const const_IWSubstring & buffer,
                int & i,
                IWString & qry_name)
{
  qry_name.resize_keep_storage(0);

  const_IWSubstring token;

  if (! buffer.nextword(token, i))
  {
    cerr << "Premature termination\n";
    return 0;
  }

  if (! token.starts_with('\''))
  {
    cerr << "Must start with a quote, got '" << token << "'\n";
    return 0;
  }

  token.remove_leading_chars(1);

  while (1)
  {
    qry_name.append_with_spacer(token);

    if (qry_name.ends_with("')"))
    {
      qry_name.chop(2);
      return 1;
    }

    if (! buffer.nextword(token, i))
    {
      cerr << "No closure to query name\n";
      return 0;
    }
  }
}

static int
tsubstructure_summarise_hits2 (const const_IWSubstring & buffer,
                               IW_STL_Hash_Multimap_IWString & count)
{
  if (buffer.nwords() < 6)
  {
    cerr << "Invalid input '" << buffer << "'\n";
    return 0;
  }

  int i = 0;

  const_IWSubstring token;

  IWString to_store;

  while (buffer.nextword(token, i))
  {
    if (iwre2::RE2PartialMatch(token, open_paren))
      break;

    to_store.append_with_spacer(token);
  }

  int matches_found = 0;

  IWString queries_matching;

  while (1)
  {
    if (! number_matches_to(buffer, i, token))
    {
      cerr << "Not '(nn matches to'\n";
      return 0;
    }

    IWString qry_name;

    if (! get_query_name(buffer, i, qry_name))
    {
      cerr << "Cannot extract query name\n";
      return 0;
    }

    if (first_token_of_multi_token_query_names)
      strip_to_first_token(qry_name);

    if (stream_for_molecule_and_hits.is_open())
      queries_matching.append_with_spacer(qry_name, queries_matching_separator);

    matches_found++;

//  count[qry_name] = to_store;


    if (append_reason)
    {
      IWString tmp(to_store);
      tmp.append_with_spacer(qry_name);
      count.insert(IW_STL_Hash_Multimap_IWString::value_type(qry_name, tmp));
    }
    else
      count.insert(IW_STL_Hash_Multimap_IWString::value_type(qry_name, to_store));

    if (! buffer.nextword(token, i))
      break;
  }

  if (stream_for_molecule_and_hits.is_open())
  {
    if (! include_smiles_in_M_file)
      to_store.remove_leading_words(1);

    stream_for_molecule_and_hits << to_store << ' ' << queries_matching << '\n';
    stream_for_molecule_and_hits.write_if_buffer_holds_more_than(4096);
  }

  if (0 == matches_found)
  {
    cerr << "Huh, no query match details found\n";
    return 0;
  }

  return 1;
}

static re2::RE2 from_tsubstructure("[0-9] matches to '");
static re2::RE2 from_iwdemerit(" : D\\([0-9]+\\)");
static re2::RE2 from_tp_first_pass(" TP1 ");

static int
tsubstructure_summarise_hits(const const_IWSubstring & buffer,
                              IW_STL_Hash_Multimap_IWString & count)
{
  if (iwre2::RE2PartialMatch(buffer, from_tsubstructure))
    return tsubstructure_summarise_hits2(buffer, count);

  if (iwre2::RE2PartialMatch(buffer, from_iwdemerit))
    return iwdemerit_summarise_hits(buffer, count);

  if (iwre2::RE2PartialMatch(buffer, from_tp_first_pass))
    return tp_first_pass_summarise_hits(buffer, count);

  unrecognised_records_encountered++;

  if (ignore_unrecognised_records)
    return 1;

  cerr << "Not sure what to do with '" << buffer << "'\n";
  return 0;
}

static int
tsubstructure_summarise_hits (iwstring_data_source & input,
                              IW_STL_Hash_Multimap_IWString & count)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (tsubstructure_summarise_hits(buffer, count))
      continue;

    if (ignore_unrecognised_records && input.eof())
      return 1;

    return 0;
  }

  return count.size();
}

static int
tsubstructure_summarise_hits (const char * fname,
                              IW_STL_Hash_Multimap_IWString & count)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tsubstructure_summarise_hits(input, count);
}

static int
translate_special_characters (IWString & fname)
{
  fname.gsub(' ', '_');
  fname.gsub('>', ".gt.");
  fname.gsub('<', ".lt.");
  fname.gsub('&', ".and.");
  fname.gsub('|', ".or.");
  fname.gsub('/', ".fs.");
  fname.gsub('\\', ".bs.");

  if (fname.starts_with('.'))
    fname.insert("DOT", 0);

  return 1;
}

static int
report_items_per_query (const IW_STL_Hash_Multimap_IWString & count,
                        std::ostream & output)
{
  IWString previous_name;

  int items_current_query = 0;

  extending_resizable_array<int> items_per_query;

  for (IW_STL_Hash_Multimap_IWString::const_iterator i = count.begin(); i != count.end(); ++i)
  {
    const IWString & qry_name = (*i).first;

    if (previous_name == qry_name)
      items_current_query++;
    else
    {
      if (previous_name.length())
        items_per_query[items_current_query]++;
      previous_name = qry_name;
      items_current_query = 1;
    }
  }

  if (items_current_query)
    items_per_query[items_current_query]++;

  for (int i = 0; i < items_per_query.number_elements(); i++)
  {
    if (items_per_query[i])
      output << items_per_query[i] << " queries had " << i << " hits\n";
  }

  return output.good();

}

static int
append_name_field_from_smiles (const IWString & buffer,
                               IWString & output)
{
  const_IWSubstring id;
  buffer.word(1, id);

  output.append_with_spacer(id);

  return 1;
}

static int
do_create_table_output (const IW_STL_Hash_Multimap_IWString & count,
                        IWString_and_File_Descriptor & output)
{
  IWString previous_name;

  int different_queries_found = 0;

  int items_current_query = 0;

  IWString names_of_molecules_matching;

  for (IW_STL_Hash_Multimap_IWString::const_iterator i = count.begin(); i != count.end(); ++i)
  {
    const IWString & qry_name = (*i).first;

//  cerr << "Examining '" << qry_name << "'\n";

    if (previous_name == qry_name)
    {
      items_current_query++;
    }
    else
    {
      different_queries_found++;
      if (items_current_query)
      {
        output << previous_name << ' ' << items_current_query;
        if (write_names_of_molecules_matched)
        {
          output << ' ' << names_of_molecules_matching;
          names_of_molecules_matching.resize_keep_storage(0);
        }
        output << '\n';
      }

      if (verbose && previous_name.length())
        cerr << items_current_query << " hits to '" << previous_name << "'\n";

      items_current_query = 1;
      previous_name = qry_name;
    }

    if (write_names_of_molecules_matched)
      append_name_field_from_smiles((*i).second, names_of_molecules_matching);

    output.write_if_buffer_holds_more_than(4096);
  }

  if (verbose && previous_name.length())
    cerr << items_current_query << " hits to '" << previous_name << "'\n";

  if (items_current_query)
  {
    output << previous_name << ' ' << items_current_query;
    if (write_names_of_molecules_matched)
      output << ' ' << names_of_molecules_matching;
    output << '\n';
  }

  return output.good();
}

static int
do_create_matches_files (const IW_STL_Hash_Multimap_IWString & count)
{
  IWString previous_name;

  IWString_and_File_Descriptor current_output_stream;

  int different_queries_found = 0;

  int items_written_current_file = 0;

  IWString fname;

  for (IW_STL_Hash_Multimap_IWString::const_iterator i = count.begin(); i != count.end(); ++i)
  {
    const IWString & qry_name = (*i).first;
    if (previous_name == qry_name)
    {
      current_output_stream << (*i).second << '\n';
      items_written_current_file++;
      continue;
    }

    different_queries_found++;

    if (current_output_stream.is_open())
    {
      current_output_stream.close();
      if (verbose > 1)
        cerr << "Wrote " << items_written_current_file << " items to '" << fname << "'\n";
    }

    if (prefix.length())
      fname = prefix;

    IWString tmp(qry_name);
    translate_special_characters(tmp);

    fname << tmp;

    if (suffix.length())
      fname << suffix;

    if (! current_output_stream.open(fname.null_terminated_chars()))
    {
      cerr << "Cannot open '" << fname << "'\n";
      return 4;
    }

    if (verbose > 2)
      cerr << "Creating '" << fname << "'\n";

    previous_name = qry_name;
    items_written_current_file = 1;

    current_output_stream << (*i).second << '\n';
  }

  current_output_stream.close();

  if (items_written_current_file)
  {
    if (verbose > 1)
      cerr << "Wrote " << items_written_current_file << " items to '" << fname << "'\n";
  }

  if (verbose)
    cerr << "Created " << different_queries_found << " files\n";

  return 1;
}


static int
tsubstructure_summarise_hits (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vP:zpfT:tnM:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('P'))
  {
    cl.value('P', prefix);

    if (verbose)
      cerr << "Files created with prefix '" << prefix << "'\n";
  }

  if (cl.option_present('z'))
  {
    ignore_unrecognised_records = 1;

    if (verbose)
      cerr << "Will ignore unrecognised records\n";
  }

  if (cl.option_present('p'))
  {
    append_reason = 1;

    if (verbose)
      cerr << "Will append the reason for the match\n";
  }

  if (cl.option_present('f'))
  {
    first_token_of_multi_token_query_names = 1;

    if (verbose)
      cerr << "Will take the first token of multi-token query names\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor stream_for_table_of_matches;

  if (cl.option_present('T'))
  {
    const char * t = cl.option_value('T');

    if (! stream_for_table_of_matches.open(t))
    {
      cerr << "Cannot open table file '" << t << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Table data on matches written to '" << t << "'\n";

    if (cl.option_present('n'))
    {
      create_matches_files = 0;
      if (verbose)
        cerr << "Individual file matches not written\n";
    }

    if (cl.option_present('t'))
    {
      write_names_of_molecules_matched = 1;

      if (verbose)
        cerr << "Table will include names of molecules matched\n";
    }
  }

  if (cl.option_present('M'))
  {
    IWString fname;
    const_IWSubstring m;
    int add_header = 0;
    for (int i = 0; cl.value('M', m, i); ++i)
    {
      if (m.starts_with("nosmi"))
      {
        include_smiles_in_M_file = 0;
        if (verbose)
          cerr << "Smiles suppressed from -M output file\n";
      }
      else if ("hdr" == m)
      {
        add_header = 1;
        include_smiles_in_M_file = 0;
      }
      else if (m.starts_with("sep="))
      {
        m.remove_leading_chars(4);
        queries_matching_separator = m;
      }
      else
      {
        fname = m;
      }
    }

    if (0 == fname.length())
    {
      cerr << "NO file name specified for per molecule match file\n";
      usage(1);
    }

    if (! stream_for_molecule_and_hits.open(fname.null_terminated_chars()))
    {
      cerr << "Cannot open stream for molecules and matches '" << fname << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "MOlecules and matches written to '" << fname << "'\n";

    if (add_header)
      stream_for_molecule_and_hits << "ID Matched_By\n";

    if (cl.option_present('n'))
    {
      create_matches_files = 0;
      if (verbose)
        cerr << "Individual file matches not written\n";
    }
  }

  IW_STL_Hash_Multimap_IWString count;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! tsubstructure_summarise_hits(cl[i], count))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
    cerr << "Read data on " << count.size() << " molecules\n";

  if (create_matches_files)
    do_create_matches_files(count);

  if (cl.option_present('T'))
    do_create_table_output(count, stream_for_table_of_matches);

  if (verbose)
    report_items_per_query(count, cerr);

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tsubstructure_summarise_hits(argc, argv);

  return rc;
}

// arch-tag: de26da19-9b6d-4ee9-936d-2d2a8ec4b00f
