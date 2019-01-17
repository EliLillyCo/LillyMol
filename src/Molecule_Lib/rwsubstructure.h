
#ifndef RW_SUBSTRUCTURE_H
#define RW_SUBSTRUCTURE_H

#include <memory>

using std::cerr;
using std::endl;

/*
  Various functions for getting Substructure_Queries (and their derived types)
  from the outside world.

  This .h file is too long, need to figure out some way of getting this into 
  a separate file.
*/

#include "iwstring.h"
#include "cmdline.h"
#include "msi_object.h"

#include "istream_and_type.h"
#include "molecule_to_query.h"
#include "aromatic.h"
#include "mdl_molecule.h"

template <typename T>
int
build_query_from_smiles (const const_IWSubstring & smiles,
                         resizable_array_p<T> & queries,
                         int verbose)
{
  Molecule m;

  set_input_aromatic_structures(1);   // don't bother saving and resetting

  if (! m.build_from_smiles(smiles))
  {
    cerr << "build_query_from_smiles:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);

  T * q = new T;

  if (! q->create_from_molecule(m, mqs))
  {
    cerr << "build_query_from_smiles:invalid molecule?? '" << smiles << "'\n";
    delete q;
    return 0;
  }

  queries.add(q);

  return 1;
}

template <typename T>
int
queries_from_file_of_molecules (MDL_Molecule & m,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<T> & queries,
                                int verbose)
{
  T * q = new T;

  if (! q->create_from_molecule(m, mqs))
  {
    cerr << "queries_from_file_of_molecules:cannot create query from '" << m.name() << "'\n";
    delete q;
    return 0;
  }

  queries.add(q);

  return 1;
}

template <typename T>
int
queries_from_file_of_molecules (data_source_and_type<MDL_Molecule> & input,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<T> & queries,
                                int verbose)
{
  set_input_aromatic_structures(1);

  MDL_Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! queries_from_file_of_molecules(*m, mqs, queries, verbose))
    {
      cerr << "Cannot create query from '" << m->name() << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Created query from '" << m->name() << "'\n";
  }

  if (verbose)
    cerr << "Read " << queries.number_elements() << " queries\n";

  return queries.number_elements();
}


/*template <typename T>
int
queries_from_file_of_isis_queries(const const_IWSubstring & fname,
                                   Molecule_to_Query_Specifications & mqs,
                                   resizable_array_p<T> & queries, 
                                   int verbose)
{
  int input_type = discern_file_type_from_name(fname);

  if (0 == input_type)
    input_type = SMI;

// Don't follow any seeking or such from the command line

  off_t o = seek_to_from_command_line();

  set_seek_to (static_cast<off_t>(0));

  data_source_and_type<MDL_Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  int rc = queries_from_file_of_molecules(input, mqs, queries, verbose);

  set_seek_to(o);

  return rc;
}*/

template <typename T>
int
queries_from_file_of_molecules (const const_IWSubstring & fname,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<T> & queries, 
                                int verbose)
{
  int input_type = discern_file_type_from_name(fname);

  if (0 == input_type)
    input_type = SMI;

// Don't follow any seeking or such from the command line

  const off_t o = seek_to_from_command_line();

  set_seek_to(static_cast<off_t>(0));

  data_source_and_type<MDL_Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  int rc = queries_from_file_of_molecules(input, mqs, queries, verbose);

  set_seek_to(o);

  return rc;
}

/*
*/

/*template <typename T>
int
queries_from_file_of_isis_queries (const const_IWSubstring & fname,
                                   resizable_array_p<T> & queries,
                                   int verbose)
{
  Molecule_to_Query_Specifications mqs;

  if (fname.contains(DIRECTIVE_SEPARATOR_TOKEN))
  {
    const_IWSubstring fname2, directives;

    fname.split (fname2, DIRECTIVE_SEPARATOR_TOKEN, directives);

    cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      cerr << "INvalid molecule to query directives '" << directives << "'\n";
      return 0;
    }

    return queries_from_file_of_molecules(fname2, mqs, queries, verbose);
  }

  return queries_from_file_of_isis_queries(fname, mqs, queries, verbose);
}*/

template <typename T>
int
query_from_ISIS_query_file(MDL_Molecule & m,
                            Molecule_to_Query_Specifications & mqs,
                            resizable_array_p<T> & queries,
                            int verbose)
{
  T * q = new T;

  if (! q->create_from_molecule(m, mqs))
  {
    delete q;
    return 0;
  }

  queries.add(q);

  if (verbose > 1 && m.name().length())
    cerr << "Created query from '" << m.name() << "'\n";

  return 1;
}

template <typename T>
int
queries_from_ISIS_query_file(data_source_and_type<MDL_Molecule> & input,
                              Molecule_to_Query_Specifications & mqs,
                              resizable_array_p<T> & queries,
                              int verbose)
{
  MDL_Molecule *m;

  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! query_from_ISIS_query_file(*m, mqs, queries, verbose))
    {
      cerr << "queries_from_ISIS_query_file:cannot process '" << m->name() << "'\n";
      return 0;
    }
  }

  return queries.number_elements();
}

template <typename T>
int
queries_from_ISIS_query_file(const const_IWSubstring & fname,
                              Molecule_to_Query_Specifications & mqs,
                              resizable_array_p<T> & queries,
                              int verbose)
{
  data_source_and_type<MDL_Molecule> input(MDL, fname);

  if (! input.good())
  {
    cerr << "queries_from_ISIS_query_file:cannot open '" << fname << "'\n";
    return 0;
  }

  input.seekg(0);    // do not allow any seeking from the command line

  int rc = queries_from_ISIS_query_file(input, mqs, queries, verbose);

  if (0 == rc)
    return 0;

  if (verbose)
    cerr << "Created " << queries.number_elements() << " queries from '" << fname << "'\n";

  return rc;
}

template <typename T>
int
queries_from_ISIS_query_file(const const_IWSubstring & fname,
                              resizable_array_p<T> & queries,
                              int verbose)
{
  Molecule_to_Query_Specifications mqs;

  if (fname.contains(DIRECTIVE_SEPARATOR_TOKEN))
  {
    const_IWSubstring fname2, directives;

    fname.split (fname2, DIRECTIVE_SEPARATOR_TOKEN, directives);

    cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      cerr << "INvalid molecule to query directives '" << directives << "'\n";
      return 0;
    }

    return queries_from_ISIS_query_file(fname2, mqs, queries, verbose);
  }

  return queries_from_ISIS_query_file(fname, mqs, queries, verbose);
}

template <typename T>
int
queries_from_file_of_molecules(const const_IWSubstring & fname,
                                resizable_array_p<T> & queries,
                                int verbose)
{
  Molecule_to_Query_Specifications mqs;

  if (fname.contains(DIRECTIVE_SEPARATOR_TOKEN))
  {
    const_IWSubstring fname2, directives;

    fname.split(fname2, DIRECTIVE_SEPARATOR_TOKEN, directives);

    cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      cerr << "INvalid molecule to query directives '" << directives << "'\n";
      return 0;
    }

    return queries_from_file_of_molecules(fname2, mqs, queries, verbose);
  }

  return queries_from_file_of_molecules(fname, mqs, queries, verbose);
}

template <typename T>
int
smarts_from_file(iwstring_data_source & input,
                  resizable_array_p<T> & queries,
                  int verbose)
{
  return smarts_or_smiles_from_file(input,
                  queries, 
                  0);  // parse as smarts
}

template <typename T>
int
smiles_from_file(iwstring_data_source & input,
                  resizable_array_p<T> & queries,
                  int verbose)
{
  return smarts_or_smiles_from_file(input,
                  queries, 
                  1);  // parse as smiles
}
template <typename T>
int
smarts_or_smiles_from_file(iwstring_data_source & input,
                  resizable_array_p<T> & queries, 
                  int smilesNotSmarts)
{
  int rc = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    buffer.strip_leading_blanks();
    if (buffer.starts_with('#'))
      continue;

    buffer.strip_trailing_blanks();

    if (0 == buffer.length())
      continue;


    T * q = new T;
    if (smilesNotSmarts)
    {
      //cerr << "Creating query from smiles '" << buffer << "'\n";

      if (! q->create_from_smiles(buffer))
      {
        cerr << "smarts_or_smiles_from_file: cannot parse '" << buffer << "'\n";
        delete q;
        return 0;
      }      
    }
    else
    {
      //cerr << "Creating query from smarts '" << buffer << "'\n";
      if (! q->create_from_smarts(buffer))
      {
        cerr << "smarts_or_smiles_from_file: cannot parse '" << buffer << "'\n";
        delete q;
        return 0;
      }
    }

    queries.add(q);
    //cerr << "Created query from '" << buffer << "'\n";
    rc++;
  }

  return rc;
}

template <typename T>
int
smarts_from_file(const const_IWSubstring & fname, resizable_array_p<T> & queries, int verbose)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "smarts_from_file: cannot open '" << fname << "'\n";
    return 0;
  }

  return smarts_from_file(input, queries, verbose);
}

template <typename T>
int
smiles_from_file(const const_IWSubstring & fname, resizable_array_p<T> & queries, int verbose)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "smiles_from_file: cannot open '" << fname << "'\n";
    return 0;
  }

  return smiles_from_file(input, queries, verbose);
}

template <typename T>
int
file_record_is_smarts(resizable_array_p<T> & queries,
                       IWString & buffer,
                       int verbose)
{
  T * tmp = new T;
  buffer.remove_leading_chars(7);

  if (! tmp->create_from_smarts(buffer))
  {
    cerr << "Invalid smarts 'SMARTS:" << buffer << "'\n";
    delete tmp;
    return 0;
  }

  if (verbose)
    cerr << "Created query '" << tmp->comment() << "' from SMARTS:" << buffer << endl;

  queries.add(tmp);

  return 1;
}

template <typename T>
int 
read_one_or_more_queries_from_file(resizable_array_p<T> & queries,
                                    iwstring_data_source & input,
                                    int verbose)
{
  off_t file_size = input.file_size();

  msi_object msi;

  int rc = 0;

  input.set_ignore_pattern("^#");
  input.set_skip_blank_lines(1);


  while (msi.read(input))
  {
    T * tmp = new T();
    if (! tmp->construct_from_msi_object(msi))
    {
      cerr << "process_queries: cannot build query from '" << msi << "'\n";
      return 0;
    }

    assert (tmp->ok());
    queries.add(tmp);

    if (verbose)
      cerr << "Created query " << (queries.number_elements() - 1) << " '" << tmp->comment() << "'\n";

    rc++;

    if (input.tellg() == file_size)   // avoids error messages in the msi.read call
      break;
  }

  return rc++;
}

template <typename T>
int 
read_one_or_more_queries_from_file(resizable_array_p<T> & queries,
                                    const const_IWSubstring & fname,
                                    int verbose)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_one_or_more_queries_from_file::cannot open '" << fname << "'\n";
    return 0;
  }

  int rc = read_one_or_more_queries_from_file(queries, input, verbose);

  if (verbose)
    cerr << "Read " << rc << " queries from '" << fname << "'\n";

  return rc;
}

template <typename T>
int
file_record_is_file(resizable_array_p<T> & queries,
                     const IWString & directory_path,
                     IWString & buffer,
                     int verbose)
{
  IWString fname;
  if (! buffer.word(0, fname))
  {
    cerr << "file_record_is_file: cannot get first word from '" << buffer << "'\n";
    return 0;
  }

  IWString pathname;
  if (fname.starts_with('/'))
    pathname=fname;
  else if (directory_path.length())
    pathname = directory_path + fname;
  else
    pathname = fname;

  T * tmp = new T;

  if (! tmp->read(pathname))
  {
    cerr << "Queries_from_file: cannot read file '" << fname << "'\n";
    delete tmp;
    return 0;
  }

  if (verbose)
    cerr << "Created query '" << tmp->comment() << "' from '" << pathname << "'\n";

  queries.add(tmp);

  return 1;
}

template <typename T>
int
queries_from_file(iwstring_data_source & input, resizable_array_p<T> & queries,
                   const IWString & directory_path,
                   int verbose)
{
  input.set_strip_leading_blanks(1);

  int rc = 0;

  IWString buffer;
  while (input.next_record(buffer))
  {
    if (0 == buffer.length())
      continue;

    if ('#' == buffer[0])
      continue;

    int rc_this_record;

    if (buffer.starts_with("SMARTS:"))
      rc_this_record = file_record_is_smarts(queries, buffer, verbose);
    else
      rc_this_record = file_record_is_file(queries, directory_path, buffer, verbose);

    if (0 == rc_this_record)
    {
      cerr << "Queries_from_file: fatal error on line " << input.lines_read() << endl;
      return 0;
    }

    rc++;
  }

  return rc;
}

/*
  Read a series of substructure queries from a file.
  We return the number read.
*/

/*template <typename T>
int
queries_from_file (const char * fname, resizable_array_p<T> & queries,
                   int inherit_directory_path,
                   int verbose)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open file '" << fname << "'\n";
    return 0;
  }

  return queries_from_file(input, queries, directory_path, verbose);
}*/

template <typename T>
int
queries_from_file (const const_IWSubstring & fname, resizable_array_p<T> & queries,
                   int inherit_directory_path,
                   int verbose)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open file '" << fname << "'\n";
    return 0;
  }

  IWString directory_path;
  if (inherit_directory_path)
  {
    int i = fname.rindex('/');
    if (i < 0)
      directory_path = "./";
    else
    {
      directory_path = fname;
      directory_path.iwtruncate(i + 1);
    }
  }

  return queries_from_file(input, queries, directory_path, verbose);
}

template <typename T>
int
queries_from_file (const IWString & fname, resizable_array_p<T> & queries,
                   int inherit_directory_path,
                   int verbose)
{
  const const_IWSubstring s = fname;

  return queries_from_file(s, queries, inherit_directory_path, verbose);
}

/*
  A Command_Line has an option which specifies one or more files containing
  lists of queries.

  Initially, one needed to have the full pathname of each query in the file.
  This was inconvenient when it came to moving things from one machine to
  the next.

  Therefore we have an option to inherit the directory path from the name
  of the containing file
*/

template <typename T>
int
process_files_of_queries (Command_Line & cl, resizable_array_p<T> & queries,
                 int inherit_directory_path,
                 int verbose, char option)
{
  int i = 0;
  int rc = 0;
  const_IWSubstring fname;
  while (cl.value(option, fname, i++))
  {
    int tmp = queries_from_file(fname, queries, inherit_directory_path, verbose);

    if (0 == tmp)
    {
      cerr << "process_files_of_queries: could not read queries from file '" <<
              cl.option_value(option, i-1) << "'\n";
      return rc;
    }

    rc += tmp;
  }

  return rc;
}

/*
  If the token starts with 'F:' it is a file of queries.
  Something starting with 'S:' is a file of smarts.
  Otherwise we assume that the query can create itself...
*/

template <typename T>
int
process_cmdline_token (char option,
                       const const_IWSubstring & token,
                       resizable_array_p<T> & queries,
                       int verbose)
{
  const_IWSubstring mytoken(token);

//cerr << "Examining token '" << mytoken << "'\n";
  if (mytoken.starts_with("F:") || mytoken.starts_with("Q:"))
  {
    mytoken.remove_leading_chars(2);

    if (0 == mytoken.length())
    {
      cerr << "Must follow S: specification with file name of queries\n";
      return 0;
    }

    if (! queries_from_file(mytoken, queries, 1, verbose))   // queries always in same directory as controlling file
    {
      cerr << "process_queries: cannot read queries from file specifier 'F:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("S:"))
  {
    mytoken.remove_leading_chars(2);

    if (0 == mytoken.length())
    {
      cerr << "Must follow S: specification with file name of queries\n";
      return 0;
    }

    if (! smarts_from_file(mytoken, queries, verbose))
    {
      cerr << "process_queries::cannot read smarts from file of smarts specifier 'S:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("M:"))
  {
    mytoken.remove_leading_chars(2);

    if (! mytoken.length())
    {
      cerr << "Must follow M: specification with file name of molecules\n";
      return 0;
    }

    if (! queries_from_file_of_molecules(mytoken, queries, verbose))
    {
      cerr << "process_queries::cannot read queries from file of molecules specifier 'M:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("smiles:"))
  {
    mytoken.remove_leading_chars(7);

    if (0 == mytoken.length())
    {
      cerr << "Must follow smiles: specification with smiles\n";
      return 0;
    }

    if (! build_query_from_smiles(mytoken, queries, verbose))
    {
      cerr << "process_queries:cannot build query from 'smiles:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("I:"))
  {
    mytoken.remove_leading_chars(2);
    if (0 == mytoken.length())
    {
      cerr << "Must follow I: specification with query file\n";
      return 0;
    }

    if (! queries_from_ISIS_query_file(mytoken, queries, verbose))
    {
      cerr << "process_queries::cannot read queries from isis query 'I:" << mytoken << "'\n";
      return 0;
    }
  }
/*else if (mytoken.starts_with("ISIS:"))
  {
    mytoken.remove_leading_chars(5);
    if (0 == mytoken.length())
    {
      cerr << "Must follow ISIS: specification with file name of queries\n";
      return 0;
    }

    if (! queries_from_file_of_isis_queries(mytoken, queries, verbose))
    {
      cerr << "process_queries::cannot read queries from file of isis queries 'ISIS:" << mytoken << "'\n";
      return 0;
    }
  }*/
  else if ("help" == mytoken)
  {
    cerr << "The following query specifications are recognised\n";
    cerr << " -" << option <<" SMARTS:smarts          smarts (use quotes to hide special characters)\n";
    cerr << " -" << option <<" S:file                 file of smarts queries\n";
    cerr << " -" << option <<" Q:file                 file of query object queries (also F: recognised)\n";
    cerr << " -" << option <<" M:file                 file of molecules that will be converted to query objects\n";
    cerr << " -" << option <<" I:file                 an ISIS query file\n";
    cerr << " -" << option <<" file                   single query file\n";
    ::exit (0);
  }
  else if (mytoken.starts_with("SMARTS:"))
  {
    mytoken.remove_leading_chars(7);

    T * q = new T;
    if (! q->create_from_smarts(mytoken))
    {
      cerr << "process_queries::invalid smarts '" << mytoken << "'\n";
      delete q;
      return 0;
    }

    queries.add(q);
  }
  else
  {
    if (! read_one_or_more_queries_from_file(queries, mytoken, verbose))
    {
      cerr << "process_queries::cannot read query/queries from '" << mytoken << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  For each occurrence of -q in a command_line object, read the accompanying
  query, and add it to the resizable array.
*/

template <typename T>
int
process_queries (Command_Line & cl, resizable_array_p<T> & queries,
                 int verbose, const char option)
{
  int nqueries = cl.option_count(option);

  if (queries.elements_allocated() < nqueries)
    queries.resize(nqueries);

  int i = 0;
  const_IWSubstring c;
  while (cl.value(option, c, i++))
  {
    if (! process_cmdline_token(option, c, queries, verbose))
    {
      cerr << "Cannot process -" << option << " option '" << c << "'\n";
      return 0;
    }
  }

  return queries.number_elements();
}

#endif
