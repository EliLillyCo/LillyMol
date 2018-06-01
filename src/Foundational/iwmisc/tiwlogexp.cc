/*
  Tester for logical expressions
*/

#include <iostream>
using std::cerr;
using std::endl;

#include "iwstring_data_source.h"
#include "cmdline.h"
#include "logical_expression.h"

static int verbose = 0;

static void
usage (int rc)
{
  cerr << "Tester for logical expressions\n";
  cerr << "Usage: tiwlogexp <file1> <file2> ...\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
test_result (IW_Logical_Expression & l,
             int i,
             const const_IWSubstring & expected_result,
             std::ostream & output,
             int & completed)
{
  int result;
  if (! l.evaluate (result))     // that's OK, must be incomplete
  {
    completed = 0;
    return 1;
  }

  completed = 1;

  if (1 == result && '1' == expected_result)
    return 1;

  if (0 == result && '0' == expected_result)
    return 1;

  output << "Result mismatch, evaluated after " << i << " results\n";
  output << " result is " << result << " expected result is " << expected_result << endl;
  l.debug_print (output);

  return 0;
}

static int
tiwlogexp (const_IWSubstring & buffer,
           int line_number,
           std::ostream & output)
{
  if (buffer.nwords () < 2)
  {
    cerr << "Input must have expression followed by expected result\n";
    return 0;
  }

  const_IWSubstring logexp, expected_result;
  buffer.word (0, logexp);
  buffer.word (1, expected_result);

  if (verbose)
    cerr << "Processing '" << logexp << "', line " << line_number << " expected result " << expected_result << endl;

  if ('0' == expected_result || '1' == expected_result)
    ;
  else
  {
    cerr << "The expected result must be '0' or '1'\n";
    return 0;
  }

  IW_Logical_Expression l;

// scan the string to extract the operators and the values

  IWString zvalues;
  zvalues.resize (logexp.nchars ());

  for (int i = 0; i < logexp.nchars (); i++)
  {
    if ('1' == logexp[i] || '0' == logexp[i])
    {
      zvalues += logexp[i];
      continue;
    }

    if ('!' == logexp[i])
    {
      l.set_unary_operator (l.number_operators (), 0);
      continue;
    }

    if (! l.add_operator (logexp[i]))
    {
      cerr << "Bad operator '" << logexp[i] << "' in '" << logexp << "'\n";
      return 0;
    }
  }

// Now feed in the results

  for (int i = 0; i < zvalues.nchars (); i++)
  {
    char v = zvalues[i];

    if ('1' == v)
      l.set_result (i, 1);
    else if ('0' == v)
      l.set_result (i, 0);
    else
    {
      cerr << "Results in test strings must be '0' or '1'. Invalid '" << logexp << "'\n";
      return 0;
    }

    int completed;
    if (! test_result (l, i + 1, expected_result, output, completed))
      return 0;

    if (completed)
      return 1;
  }

  int completed;
  if (! test_result (l, zvalues.nchars (), expected_result, output, completed))
  {
    l.debug_print (cerr);
    return 0;
  }

  if (! completed)
  {
    cerr << "Yipes, cannot evaluate completed expression\n";
    l.debug_print (cerr);
    return 0;
  }

  return 1;
}

static int
tiwlogexp (iwstring_data_source & input, std::ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer) && output.good ())
  {
    buffer.strip_leading_blanks ();
    buffer.strip_trailing_blanks ();

    if (0 == buffer.length ())
      continue;

    if (buffer.starts_with ('#'))
      continue;

    if (! tiwlogexp (buffer, input.lines_read (), output))
      return 0;
  }

  return output.good ();
}


static int
tiwlogexp (const char * fname, std::ostream & output)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tiwlogexp (input, output);
}

static int
tiwlogexp (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (2);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    return 1;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! tiwlogexp (cl[i], std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = tiwlogexp (argc, argv);

  return rc;
}
