#include "molecule.h"
#include "qry_wstats.h"
#include "accumulator.h"
#include "istream_and_type.h"

/*
 * By default, the isotopic labels are just the line number in the file.
 * If the control file contains 'ISO=nnn' we set that for the isotopic label
 *
 * Each atom type can have a confidence level.
 * These range from 0 (just a guess) to 100 (unambiguous)
 *
 * We also offer the option of NOT marking the atoms as processed. this
 * allows the possibility of multiple queries hitting the same atoms
 *
 * We allow a variable number of atoms to match per query. 
 *
 * By default, we are an atom additivity method, and so just process
 * the first atom in any embedding
 */

extern int max_atoms_per_embedding_to_process;
extern resizable_array<int> atoms_to_match;
extern resizable_array<int> mark_atoms;
extern resizable_array<int> label;
extern resizable_array<int> confidence;
extern resizable_array_p<Substructure_Hit_Statistics> queries;

extern int verbose;

/*
 *  Extract a token like PREFIX=nnn
 *  from buffer
 *  We return 0 if an error is found
 */

static int
extract_qualifier (const const_IWSubstring & buffer,
                   const char * prefix,
                   int & result)
{
  int i = buffer.find (prefix);

  if (i < 0)
    return 1;     // no value specified is OK

  const_IWSubstring t = buffer.substr (i + strlen (prefix));

  if (t.starts_with ('='))
    t++;

  if (0 == t.length ())
  {
    cerr << "Cannot extract '" << prefix << "' value from '" << buffer << "'\n";
    return 0;
  }

  if (t.contains (' '))
    t.truncate_at_first (' ');

  if (! t.numeric_value (result))
  {
    cerr << "Invalid numeric for '" << prefix << "' value in '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
parse_control_file_record (const const_IWSubstring & buffer)
{
  if (buffer.nwords () < 3)
  {
    cerr << "Control file records must have at least 3 words\n";
    return 0;
  }

  const_IWSubstring fname_smarts, nv, comments;

  buffer.word (0, fname_smarts);

  buffer.word (1, nv);
  double nvv;
  if (! nv.numeric_value (nvv))
  {
    cerr << "Bad numeric value '" << nv << "'\n";
    return 0;
  }

  comments = buffer;
  comments.remove_leading_words (2);

  Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;

  if (fname_smarts.starts_with ("f:"))    // is a file name
  {
    if (! q->read (fname_smarts))
    {
      cerr << "Cannot read query from '" << fname_smarts << "'\n";
      return 0;
    }
  }
  else
  {
    if (! q->create_from_smarts (fname_smarts))
    {
      cerr << "Cannot read smarts '" << fname_smarts << "'\n";
      return 0;
    }
  }

  q->add_numeric_value (nvv);

  q->set_comment (comments);

  queries.add (q);

  if (verbose > 1)
    cerr << "Created query '" << comments << "'\n";

  if (comments.find ("NOMARK") >= 0)
    mark_atoms.add (0);     // don't mark the matched atoms (allows multiple matches);
  else
    mark_atoms.add (1);     // by default we mark all matched atoms

  int iso = -1;

  if (! extract_qualifier (buffer, "ISO=", iso))    /// error
    return 0;

  if (iso > 0)
    label.add (iso);
  else
    label.add (queries.number_elements ());

  int conf = -5;

  if (! extract_qualifier (buffer, "CONF=", conf))    /// error
    return 0;

  if (conf >= 0)
    confidence.add (conf);
  else
    confidence.add (100);

  int nmatch = max_atoms_per_embedding_to_process;

  if (! extract_qualifier (buffer, "NMATCH=", nmatch))  /// error
    return 0;

  atoms_to_match.add (nmatch);

  return 1;
}

static int
read_control_file (iwstring_data_source & input)
{
  int nr = input.records_remaining ();

  queries.resize (queries.elements_allocated () + nr);
  label.resize (queries.elements_allocated () + nr);
  confidence.resize (queries.elements_allocated () + nr);
  mark_atoms.resize (queries.elements_allocated () + nr);
  atoms_to_match.resize (queries.elements_allocated () + nr);

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    buffer.strip_trailing_blanks ();
    buffer.strip_leading_blanks ();
    if (buffer.starts_with ('#'))
      continue;

    if (0 == buffer.length ())
      continue;

    if (verbose > 1)
      cerr << "Examining " << input.lines_read () << " '" << buffer << "'\n";

    if (! parse_control_file_record (buffer))
    {
      cerr << "Cannot parse control file record '" << buffer << "'\n";
      return 0;
    }
  }

  assert (queries.number_elements () == mark_atoms.number_elements ());
  assert (queries.number_elements () == label.number_elements ());

  return queries.number_elements ();
}

int
read_control_file (const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open control file '" << fname << "'\n";
    return 0;
  }

  return read_control_file (input);
}

