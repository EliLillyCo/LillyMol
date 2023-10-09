/*
  Check and report on characteristics of a tdt
*/

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/data_source/iwstring_data_source.h"
using std::cerr;
using std::cout;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int tdts_read = 0;

static extending_resizable_array<int> items_in_tdt;

static IW_STL_Hash_Map_int tdts_with_tag;

typedef extending_resizable_array<int> extdintarray;

typedef IW_STL_Hash_Map<IWString, extdintarray> IW_STL_Hash_Map_extdingarray;

static IW_STL_Hash_Map_extdingarray tag_count;

//template _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > ** fill_n<_Hashtable_node<pair<IWString const, extending_resizable_array<int> > > **, unsigned int, _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > *>(_Hashtable_node<pair<IWString const, extending_resizable_array<int> > > **, unsigned int, _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > * const &);
//template void fill<_Hashtable_node<pair<IWString const, extending_resizable_array<int> > > **, _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > *>(_Hashtable_node<pair<IWString const, extending_resizable_array<int> > > **, _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > **, _Hashtable_node<pair<IWString const, extending_resizable_array<int> > > * const &);
//template class hashtable<pair<IWString const, extending_resizable_array<int> >, IWString, IWhash<IWString>, _Select1st<pair<IWString const, extending_resizable_array<int> > >, equal_to<IWString>, allocator<extending_resizable_array<int> > >;
//template class vector<_Hashtable_node<pair<IWString const, extending_resizable_array<int> > > *, allocator<extending_resizable_array<int> > >;

//template class _Hashtable_const_iterator<pair<IWString const, extending_resizable_array<int> >, IWString, IWhash<IWString>, _Select1st<pair<IWString const, extending_resizable_array<int> > >, equal_to<IWString>, allocator<extending_resizable_array<int> > >;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Scans a TDT file and does quality checking and reporting\n";
  cerr << " -v           verbose output\n";

  exit (rc);
}

static int
check_tdt (iwstring_data_source & input)
{
  IW_STL_Hash_Map_int found_this_tdt;

  const_IWSubstring buffer;

  int items_current_tdt = 0;
  while (input.next_record (buffer))
  {
    if ('|' == buffer)
    {
      tdts_read++;
      for (IW_STL_Hash_Map_int::const_iterator i = found_this_tdt.begin (); i != found_this_tdt.end (); i++)
      {
        const IWString & tag = (*i).first;
        tdts_with_tag[tag]++;
        int          count   = (*i).second;
        extdintarray & e = tag_count[tag];
        e[count]++;
      }
      found_this_tdt.clear ();

      items_in_tdt[items_current_tdt]++;
      items_current_tdt = 0;
    }
    else if (! buffer.ends_with ('>'))
    {
      cerr << "Invalid TDT record '" << buffer << "'\n";
      cerr << "line " << input.lines_read () << endl;
      return 0;
    }
    else
    {
      buffer.truncate_at_first ('<');
      found_this_tdt[buffer]++;
      items_current_tdt++;
    }
  }

  if (items_current_tdt)
    cerr << "Last TDT not terminated - ignored\n";

  return 1;
}

static int
check_tdt (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return check_tdt (input);
}

static int
check_tdt (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! check_tdt (cl[i]))
    {
      rc = i + 1;
      break;
    }
  }

  cout << "Read " << tdts_read << " TDT's\n";
  for (int i = 0; i < items_in_tdt.number_elements (); i++)
  {
    if (items_in_tdt[i])
    {
      cout << items_in_tdt[i];
      if (tdts_read == items_in_tdt[i])
        cout << " (all)";
      cout << " TDT's had " << i << " dataitems\n";
    }
  }

  for (IW_STL_Hash_Map_int::const_iterator i = tdts_with_tag.begin (); i != tdts_with_tag.end (); i++)
  {
    const IWString & tag = (*i).first;
    int count = (*i).second;

    cout << "Tag '" << tag << "' appeared in " << count << " TDTS";
    if (count == tdts_read)
      cout << " - all";
    cout << endl;
  }

  for (IW_STL_Hash_Map_extdingarray::const_iterator i = tag_count.begin (); i != tag_count.end (); i++)
  {
    const IWString & tag = (*i).first;
    const extdintarray & e = (*i).second;

    for (int j = 0; j < e.number_elements (); j++)
    {
      if (e.item (j))
      {
        cout << e[j];
        if (tdts_read == e[j])
          cout << " (all)";
        else
          cout << "      ";
        cout << " tdts had " << j << " ocurrences of '" << tag << "'\n";
      }
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = check_tdt (argc, argv);

  return rc;
}
