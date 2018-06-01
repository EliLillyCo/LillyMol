#include <stdlib.h>

/*
  Tester for STL stuff
*/

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "iw_stl_hash_map.h"

using std::cerr;
using std::cout;
using std::ostream;

static void
usage (int rc)
{
  cerr << "Tester for IWString stl hash maps\n";
  cerr << "Input file needs two columns: key value\n";
  cerr << " -e              echo the hash\n";
  cerr << " -n              test numeric hash\n";
  cerr << " -F <string>     check for <string> as a key in the hash\n";
  cerr << " -v              verbose output\n";

  exit (rc);
}

static int verbose = 0;

static int echo_hash = 0;

static int test_numeric = 0;

/*
  Once the hash-map is built, we can test for the presence of an item
*/

static IWString string_to_find;

typedef IWString mystringtype;

static int
tstl_string (iwstring_data_source & input,
             ostream & output)
{
  input.set_strip_trailing_blanks ();
  input.set_skip_blank_lines ();
  input.set_ignore_pattern ("^#");

  IW_STL_Hash_Map_String foo;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.nwords () < 2)
    {
      cerr << "Hmmm, need at least 2 words on each line\n";
      return 0;
    }

    IWString key = buffer.word (0);
    mystringtype kvalue = buffer.word (1);

    foo[key] = kvalue;

    output << "Key '" << key << "' value '" << kvalue << "' becomes '" << foo[key] << "'\n";
    foo[key] += " changed";
    output << " Now " << foo[key] << endl;
  }

  if (echo_hash)
  {
    cerr << "Hash contains " << foo.size () << " items\n";

    int item_count = 0;

    IW_STL_Hash_Map_String::const_iterator f;
    for (f = foo.begin (); f != foo.end (); f++)
    {
      cerr << (*f).first << " value " << (*f).second << endl;

      item_count++;
    }

    if (item_count != foo.size ())
    {
      cerr << "Huh, hash contains " << foo.size () << " items, but we counted " << item_count << endl;
    }
  }

  if (0 == string_to_find.length ())
    return output.good ();

  IW_STL_Hash_Map_String::const_iterator f = foo.find (string_to_find);
  if (f == foo.end ())
    cerr << "Item '" << string_to_find << "' not in hash\n";
  else
    cerr << "Found '" << (*f).second << "'\n";

  int tmp = 0; //foo.contains (string_to_find);
  cerr << "The contains member function reports " << tmp << endl;
//cerr << "The contains member function reports " << foo.contains (string_to_find) << endl;

  return output.good ();
}

typedef long mytype;

static int
tstl_numeric (iwstring_data_source & input,
              ostream & output)
{
  input.set_strip_trailing_blanks ();
  input.set_skip_blank_lines ();
  input.set_ignore_pattern ("^#");

  IW_STL_Hash_Map<IWString, mytype> foo;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.nwords () < 2)
    {
      cerr << "Hmmm, need at least 2 words on each line\n";
      return 0;
    }

    IWString key = buffer.word (0);
    const_IWSubstring kvalue = buffer.word (1);

    mytype x;
    if (! kvalue.numeric_value (x))
    {
      cerr << "Invalid numeric '" << kvalue << "'\n";
      return 0;
    }

    foo[key] = x;

    cerr << "Set array value for '" << key << "' to " << foo[key] << " from '" << x << "'\n";
    foo[key] += 1;
    cerr << " mow " << foo[key] << endl;

  }

  return output.good ();
}

static int
tstl (const char * fname,
      ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (test_numeric)
    return tstl_numeric (input, output);
  else 
    return tstl_string (input, output);
}

static int
tstl (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vnF:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options present\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('n'))
  {
    test_numeric = 1;

    cerr << "Will test numeric functionality ";
    mytype tmp = 1;
    if (tmp + 0.1 == tmp)
      cerr << "integer\n";
    else
      cerr << "float\n";
  }

  if (cl.option_present ('F'))
  {
    string_to_find = cl.string_value ('F');

    if (verbose)
      cerr << "Will look for key '" << string_to_find << "'\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = tstl (cl[0], cout);

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = tstl (argc, argv);

  return rc;
}
