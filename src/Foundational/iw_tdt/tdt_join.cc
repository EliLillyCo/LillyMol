/*
  We have two TDT files and need a merged TDT
*/

#include <iostream>
#include <unordered_map>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwhash.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iw_tdt/iw_tdt.h"

const char * prog_name = nullptr;

static int verbose = 0;

static IWString identifier_tag1("PCN<");
static int identifier_column1 = -1;
static IWString identifier_tag2("PCN<");
static int identifier_column2 = -1;

static int remove_duplicate_data = 0;

static int ignore_duplicate_identifiers = 0;

static int tdts_read = 0;
static int tdts_written = 0;
static int trim_leading_zeros_from_identifiers = 0;

static int only_write_records_when_all_files_have_data = 0;
static int records_skipped_because_of_incomplete_data = 0;

static Report_Progress report_progress;

class TDT_File : public std::unordered_map<IWString, off_t, IWStringHash>, public IW_TDT
{
  private:
    iwstring_data_source _input;

//  private functions

    int _determine_offsets();

  public:
    int determine_offsets(const char * fname);

    int fetch_data_for(const IWString & id);
};

static int
extract_identifier(const IW_TDT & tdt,
                    const IWString & tag,
                    int col,
                    IWString & id)
{
  if (! tdt.dataitem_value(tag, id))
  {
    cerr << "Cannot extract data for '" << tag << "'\n";
    return 0;
  }

  if (col < 0)
    return 1;

  int nw = id.nwords();

  if (col >= nw)
  {
    cerr << "Cannot extract column " << (col + 1) << " from '" << id << "'\n";
    return 0;
  }

  if (col > 0)
    id.remove_leading_words(col);

  id.truncate_at_first(' ');

  if (trim_leading_zeros_from_identifiers)
    id.remove_leading_chars('0');

  return id.length();
}

int
TDT_File::_determine_offsets()
{
  off_t o = _input.tellg();

  IW_TDT tdt;

  while (tdt.next(_input))
  {
    IWString id1;

    if (! extract_identifier(tdt, identifier_tag2, identifier_column2, id1))
    {
      cerr << "Bad news, cannot extract identifier data '" << identifier_tag2 << "' from TDT\n";
      cerr << tdt;
      return 0;
    }

//  std::unordered_map<IWString, off_t, IWStringHash>::const_iterator f = find(id1);
    const auto f = find(id1);
    if (f != end())
    {
      cerr << "TDT_File::_determine_offsets: duplicate identifier '" << id1 << "'\n";
      if (! ignore_duplicate_identifiers)
        return 0;
    }
    else
    {
      operator[] (id1) = o;
      if (report_progress())
        cerr << "Accumulated offset data on " << size() << " identifiers\n";
    }

    o = _input.tellg();
  }

  return size();
}

int
TDT_File::determine_offsets(const char * fname)
{
  if (! _input.open(fname))
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return _determine_offsets();
}

int
TDT_File::fetch_data_for(const IWString & id)
{
//std::unordered_map<IWString, off_t, IWStringHash>::const_iterator f = find(id);
  const auto f = find(id);

  if (f == end())
    return 0;

  const off_t o = (*f).second;

  if (! _input.seekg(o))
  {
    cerr << "TDT_File::fetch_data_for:cannot seek to " << o << " for '" << id << "'\n";
    return 0;
  }

  return IW_TDT::next(_input);
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Merges TDT files\n";
  cerr << " -t <tag>       identifier tag in first file\n";
  cerr << " -c <col>       identifier is column <col> in first file\n";
  cerr << " -T <tag>       identifier tag in second file if different\n";
  cerr << " -C <col>       identifier is column <col> in second file\n";
  cerr << " -I             only write TDT's for which identifier is present in every file\n";
  cerr << " -g             ignore duplicate identifiers in files\n";
  cerr << " -d             eliminate duplicate tags - info from 2nd+ file discarded\n";
  cerr << " -r <number>    report progress during offset scans\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
write_without_vbar(const IW_TDT & tdt,
                   IWString_and_File_Descriptor & output)
{
  const IWString & zdata = tdt.rawdata();

  if (zdata.length() <= 2)
    return output.good();

  output.write(zdata.rawchars(), zdata.length() - 2);

  return output.good();
}

static int
do_remove_duplicate_data(const IW_TDT & tdt,
                         TDT_File * files,
                         const int nfiles)
{
  const int n = tdt.number_elements();

  for (int i = 0; i < n; i++)
  {
    const_IWSubstring tag;

    tdt.item(i, tag);

    const int oab = tag.index('<');

    tag.iwtruncate(oab + 1);

//  cerr << "Removing tag '" << tag << "'\n";

    for (int j = 0; j < nfiles; j++)
    {
      files[j].remove_items_with_tag(tag);
    }
  }

  return 1;
}

static int
tdt_join(IW_TDT & tdt,
         const IWString & id,
         TDT_File * files,
         const int nfiles,
         IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < nfiles; i++)
  {
    if (files[i].fetch_data_for(id))
      continue;

    if (only_write_records_when_all_files_have_data)
    {
      records_skipped_because_of_incomplete_data++;
      if (verbose > 1)
        cerr << "No data for '" << id << "' in file " << i << endl;
      return output.good();
    }
  }

// At this stage, we have the data for each TDT

  if (remove_duplicate_data)
    do_remove_duplicate_data(tdt, files, nfiles);

  write_without_vbar(tdt, output);

  for (int i = 0; i < nfiles; i++)
  {
    write_without_vbar(files[i], output);
  }

  output << "|\n";

  tdts_written++;

  return output.good();
}

static int
tdt_join(IW_TDT & tdt,
         TDT_File * files,
         const int nfiles,
         IWString_and_File_Descriptor & output)
{
  IWString id;
  if (! extract_identifier(tdt, identifier_tag1, identifier_column1, id))
  {
    cerr << "Cannot determine identifier\n";
    cerr << tdt;
    return 0;
  }

  return tdt_join(tdt, id, files, nfiles, output);
}

static int
tdt_join(iwstring_data_source & input,
         TDT_File * files,
         const int nfiles,
         IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    tdts_read++;

    if (! tdt_join(tdt, files, nfiles, output))
    {
      cerr << "Fatal error processing TDT at line " << input.lines_read() << endl;
      cerr << tdt;
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

static int
tdt_join(const char * fname,
         TDT_File * files,
         const int nfiles,
         IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return tdt_join(input, files, nfiles, output);
}

static int
tdt_join(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vt:c:T:C:Idgzr:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('t'))
  {
    cl.value('t', identifier_tag1);
    if (verbose)
      cerr << "First file identifiers in the '" << identifier_tag1 << "' tag\n";

    if (! identifier_tag1.ends_with('<'))
      identifier_tag1 += '<';
  }

  if (cl.option_present('T'))
  {
    cl.value('T', identifier_tag2);
    if (verbose)
      cerr << "Second file identifiers in the '" << identifier_tag2 << "' tag\n";

    if (! identifier_tag2.ends_with('<'))
      identifier_tag2 += '<';
  }
  else
  {
    identifier_tag2 = identifier_tag1;
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "The report progress option (-r) must be a whole +ve number\n";
      usage(4);
    }
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', identifier_column1) || identifier_column1 < 1)
    {
      cerr << "INvalid identifier column (-c option) specification\n";
      usage(4);
    }

    if (verbose)
      cerr << "Identifier token for file 1 in column " << identifier_column1 << endl;

    identifier_column1--;
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', identifier_column2) || identifier_column2 < 1)
    {
      cerr << "INvalid identifier column (-C option) specification\n";
      usage(4);
    }

    if (verbose)
      cerr << "Identifier token for file 2 in column " << identifier_column2 << endl;

    identifier_column2--;
  }

  if (cl.option_present('d'))
  {
    remove_duplicate_data = 1;

    if (verbose)
      cerr << "Tags already present in first file will be removed\n";
  }

  if (cl.option_present('I'))
  {
    only_write_records_when_all_files_have_data = 1;

    if (verbose)
      cerr << "Will only write records when each input file has data for the identifier\n";
  }

  if (cl.option_present('z'))
  {
    trim_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Leading zero's will be trimmed from identifiers\n";
  }

  if (cl.option_present('g'))
  {
    ignore_duplicate_identifiers = cl.option_count('g');

    if (verbose)
      cerr << "Will ignore duplicate identifiers\n";
  }

  int nfiles = cl.number_elements();

  if (nfiles < 2)
  {
    cerr << "Joins TDT files, must have at least two arguments\n";
    usage(2);
  }

  TDT_File * files = new TDT_File[nfiles - 1]; std::unique_ptr<TDT_File[]> free_files(files);

  assert (nullptr != files);

  for (int i = 1; i < cl.number_elements(); i++)
  {
    if (! files[i - 1].determine_offsets(cl[i]))
    {
      cerr << "Cannot determine offsets in '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc;

  if (! tdt_join(cl[0], files, nfiles - 1, output))
  {
    cerr << "Failure processing '" << cl[0] << "'\n";
    rc = 3;
  }
  else
  {
    rc = 0;
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << tdts_read << " tdts, wrote " << tdts_written << endl;
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tdt_join(argc, argv);

  return rc;
}
