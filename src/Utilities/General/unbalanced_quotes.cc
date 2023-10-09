#include <stdlib.h>
#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/cmdline/cmdline.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int changes = 0;

static void
usage(int rc)
{
  cerr << "Identifies unbalanced double quotes in a file\n";
  cerr << " -t <from>-<to>     convert occurrences of <from> to <to> when quoted\n";
  cerr << " -p <ud>            test for balanced pairs of Up and Down chars\n";
  cerr << "                    the following ar recognised 'paren', 'sqb', 'brace'\n";
  cerr << " -v                 verbose output\n";

  exit(rc);
}

static int
unbalanced_whatevers(const const_IWSubstring & buffer,
                     const char up,
                     const char down)
{
  int inside_dquote = 0;
  int inside_squote = 0;
  int xlevel = 0;

  for (int i = 0; i < buffer.length(); ++i)
  {
    const char c = buffer[i];

    if ('"' == c)
    {
      if (inside_dquote)
        inside_dquote = 0;
      else
        inside_dquote = 1;
    }
    else if ('\'' == c)
    {
      if (inside_squote)
        inside_squote = 0;
      else
        inside_squote = 1;
    }
    else if (inside_dquote || inside_squote)
      ;
    else if (up == c)
      xlevel++;
    else if (down == c)
      xlevel--;
  }

  if (0 == xlevel)
    return 0;

  return 1;   // something wrong
}

// Sometime I need to allow for "" being an escaped double quote

static int
unbalanced_quotes(const const_IWSubstring & buffer)
{
  int inside_quote = 0;

  for (int i = 0; i < buffer.length(); i++)
  {
    if ('"' != buffer[i])
      continue;

    if (i > 0 && '\\' == buffer[i - 1])
    {
      if (verbose > 1)
        cerr << "Escaped double quote found\n";
      continue;
    }

    if (i > 0 && '\'' == buffer[i - 1] && i < buffer.length() - 1 && '\'' == buffer[i + 1])
    {
      if (verbose > 1)
        cerr << "Singly quoted double quote found\n";

      continue;
    }

    if (inside_quote)
      inside_quote = 0;
    else
      inside_quote = 1;
  }

  if (inside_quote)
    return 1;

  return 0;
}

static int
unbalanced_quotes(IWString & buffer,
                  const char * translation_table,
                  IWString_and_File_Descriptor & output)
{
  int inside_quote = 0;

  for (int i = 0; i < buffer.length(); i++)
  {
    const auto c = static_cast<int>(buffer[i]);

//  cerr << " c = '" << c << "' tr '" << translation_table[c] << "'\n";
    if (inside_quote && '\0' != translation_table[c])
    {
//    cerr << "Change '" << c << "' to " << translation_table[c] << "'\n";
      buffer[i] = translation_table[c];
      changes++;
    }

    if ('"' != buffer[i])
      continue;

    if (i > 0 && '\\' == buffer[i - 1])
    {
      if (verbose > 1)
        cerr << "Escaped double quote found\n";
      continue;
    }

    if (i > 0 && '\'' == buffer[i - 1] && i < buffer.length() - 1 && '\'' == buffer[i + 1])
    {
      if (verbose > 1)
        cerr << "Singly quoted double quote found\n";

      continue;
    }

    if (inside_quote)
      inside_quote = 0;
    else
      inside_quote = 1;
  }

  if (inside_quote)
    return 1;

  return 0;
}

static int
unbalanced_whatevers(iwstring_data_source & input,
                    const char up,
                    const char down)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (unbalanced_whatevers(buffer, up, down))
    {
      cerr << "Unbalanced " << up << down << " pair(s) found on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
unbalanced_quotes (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! buffer.contains ('"'))
      continue;

    if (unbalanced_quotes (buffer))
    {
      cerr << "Unbalanced quotes found on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
unbalanced_quotes (iwstring_data_source & input,
                   const char * translation_table,
                   IWString_and_File_Descriptor & output)
{
  IWString buffer;

  while (input.next_record (buffer))
  {
    if (! buffer.contains('"'))
      output << buffer << '\n';
    else if (! unbalanced_quotes(buffer, translation_table, output))
      output << buffer << '\n';
    else
    {
      cerr << "Unbalanced quotes found\n";
      cerr << buffer << endl;
      cerr << "Line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
unbalanced_quotes (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return unbalanced_quotes (input);
}

static int
unbalanced_quotes (const char * fname,
                   const char * translation_table,
                   IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return unbalanced_quotes (input, translation_table, output);
}

static int
unbalanced_whatevers (const char * fname,
                      const char up,
                      const char down)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return unbalanced_whatevers (input, up, down);
}

static int
convert_to_char (IWString & s)
{
  if (1 == s.length())
    return 1;

  if ("space" == s)
  {
    s = ' ';
    return 1;
  }

  if ("tab" == s)
  {
    s = '\t';
    return 1;
  }

  if ("vbar" == s)
  {
    s = '|';
    return 1;
  }

  if ("comma" == s)
  {
    s = ',';
    return 1;
  }

  if ("squote" == s)
  {
    s = '\'';
    return 1;
  }

  if ("dquote" == s)
  {
    s = '"';
    return 1;
  }

  return 0;
}

static int
unbalanced_quotes (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vt:p:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  char translation_table[256];

  for (int i = 0; i < 256; ++i)
  {
    translation_table[i] = '\0';
  }

  if (cl.option_present('t'))
  {
    const_IWSubstring s;
    for (int i = 0; cl.value('t', s, i); ++i)
    {
      IWString f,t;
      if (! s.split(f, '-', t))
      {
        cerr << "Translate directives must be of the form 'f-t', '" << s << "' invalid\n";
        usage(1);
      }
      if (1 == f.length() && 1 == t.length())
        translation_table[static_cast<int>(f[0])] = t[0];
      else
      {
        if (! convert_to_char(f) || ! convert_to_char(t))
        {
          cerr << "Unrecognised translation directive '" << f << "' and/or '" << t << "'\n";
          usage(1);
        }

        translation_table[static_cast<int>(f[0])] = t[0];
        if (verbose)
          cerr << "Will translate quoted '" << f[0] << "' to '" << translation_table[static_cast<int>(f[0])] << "'\n";
      }
    }
  }

  char up = ' ';
  char down = ' ';

  if (cl.option_present('p'))
  {
    const_IWSubstring p = cl.string_value('p');

    if ("paren" == p)
    {
      up = '(';
      down = ')';
    }
    else if ("sqb" == p)
    {
      up = '[';
      down = ']';
    }
    else if ("brace" == p)
    {
      up = '{';
      down = '}';
    }
    else if (2 != p.length())
    {
      cerr << "The -p option must be one of 'paren', 'sqb', 'brace' or a two character combination, '" << p << "' invalid\n";
      usage(1);
    }
    else
    {
      up = p[0];
      down = p[1];
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  if (cl.option_present('t'))
  {
    IWString_and_File_Descriptor output(1);
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! unbalanced_quotes(cl[i], translation_table, output))
      {
        rc = i + 1;
        break;
      }
    }
    output.flush();
  }
  else if (' ' != up)
  {
    for (int i = 0; i < cl.number_elements(); ++i)
    {
      if (! unbalanced_whatevers(cl[i], up, down))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! unbalanced_quotes(cl[i]))
      {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose && cl.option_present('t'))
    cerr << "Made " << changes << " translation changes\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = unbalanced_quotes(argc, argv);

  return rc;
}
