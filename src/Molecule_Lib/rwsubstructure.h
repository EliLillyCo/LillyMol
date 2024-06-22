#ifndef MOLECULE_LIB_RWSUBSTRUCTURE_H_
#define MOLECULE_LIB_RWSUBSTRUCTURE_H_

#include <cassert>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

/*
  Various functions for getting Substructure_Queries (and their derived types)
  from the outside world.

  This .h file is too long, need to figure out some way of getting this into 
  a separate file.
*/

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/msi_object.h"
#include "Foundational/iwmisc/proto_support.h"

#include "istream_and_type.h"
#include "molecule_to_query.h"
#include "aromatic.h"
#include "mdl_molecule.h"

template <typename T>
T *
ReadProtoQueryFile(const const_IWSubstring fname)
{
  IWString s(fname);
  std::optional<SubstructureSearch::SubstructureQuery> maybe_proto =
      iwmisc::ReadTextProtoCommentsOK<SubstructureSearch::SubstructureQuery>(s);
  if (! maybe_proto) {
    return 0;
  }

  T * to_be_returned = new T;

  if (! to_be_returned->ConstructFromProto(*maybe_proto)) 
  {
    std::cerr << "ReadProtoQueryFile:cannot build query from proto\n";
    std::cerr << maybe_proto->ShortDebugString() << '\n';
    delete to_be_returned;
    return nullptr;
  }

  return to_be_returned;
}

// Normally, if there are multiple "queries {" directives in a file, that
// means multiple components of a single composite query. But it is also
// convenient to have multiple textproto queries in a file.
// Read the contents of `fname` interpreting individual "query" entries
// as separate queries.
template <typename T>
int
MultipleQueriesFromTextProto(const const_IWSubstring fname,
                             int verbose,
                             resizable_array_p<T>& queries) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    std::cerr << "MultipleQueriesFromTextProto:cannot open '" << fname;
    return 0;
  }

  return MultipleQueriesFromTextProto(input, verbose, queries);
}

template <typename T>
int
MultipleQueriesFromTextProto(iwstring_data_source& input,
                             int verbose,
                             resizable_array_p<T>& queries) {
  while (true) {
    std::optional<std::string> string_proto = iwsubstructure::GetNextQueryTextProto(input);
    if (! string_proto) {
      break;
    }

    SubstructureSearch::SubstructureQuery proto;
    if (! google::protobuf::TextFormat::ParseFromString(*string_proto, &proto)) {
      std::cerr << "MultipleQueriesFromTextProto:cannot parse proto\n";
      std::cerr << *string_proto << '\n';
      return 0;
    }

    std::unique_ptr<T> query = std::make_unique<T>();
    if (! query->ConstructFromProto(proto)) {
      std::cerr << "MultipleQueriesFromTextProto:cannot build query\n";
      std::cerr << *string_proto << '\n';
      return 0;
    }
    if (query->comment().empty()) {
      query->set_comment(query->item(0)->comment());
    }
    queries << query.release();
  }

  if (verbose) {
    std::cerr << "MultipleQueriesFromTextProto::read " << queries.size() << " queries\n";
  }

  return queries.number_elements();
}

template <typename T>
int
ReadFileOfProtoQueries(const const_IWSubstring fname, resizable_array_p<T>& queries)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    std::cerr << "ReadFileOfProtoQueries:cannot open '" << fname << "'\n";
    return 0;
  }

  const IWString dirname = iwmisc::IWDirname(fname);

  IWString line;
  while (input.next_record(line))
  {
    std::optional<IWString> file_that_exists = iwmisc::FileOrPath(fname, line);
    if (! file_that_exists) {
      std::cerr << "ReadFileOfProtoQueries:no file '" << line << "'\n";
      return 0;
    }
//  IWString fname;
//  fname << dirname << std::filesystem::path::preferred_separator << line;
    T* q = ReadProtoQueryFile<T>(*file_that_exists);

    if (nullptr == q) {
      std::cerr << "ReadFileOfProtoQueries:cannot process '" << line << "'\n";
      return 0;
    }

    queries.add(q);
  }

  return 1;
}

template <typename T>
int
ReadQueriesFromTFDataRecordProtoFile(iw_tf_data_record::TFDataReader& input,
                        resizable_array_p<T>& queries) {
  while (1) {
    std::optional<SubstructureSearch::SubstructureQuery> proto =
         input.ReadProto<SubstructureSearch::SubstructureQuery>();
    if (! proto) {
      return queries.number_elements();
    }
    std::unique_ptr<T> query = std::make_unique<T>();
    if (! query->ConstructFromProto(*proto)) {
      std::cerr << "ReadQueriesFromTFDataRecordProtoFile:cannot build query\n";
      return 0;
    }
    queries.add(query.release());
  }
  // Should never come here.
}

template <typename T>
int
ReadQueriesFromTFDataRecordProtoFile(const const_IWSubstring& fname,
                        resizable_array_p<T>& queries) {
  iw_tf_data_record::TFDataReader input(fname);
  if (! input.good()) {
    std::cerr << "ReadQueriesFromTFDataRecordProtoFile:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadQueriesFromTFDataRecordProtoFile(input, queries);
}
template <typename T>
int
build_query_from_smiles(const const_IWSubstring & smiles,
                        resizable_array_p<T> & queries,
                        int verbose)
{
  Molecule m;

  set_input_aromatic_structures(1);   // don't bother saving and resetting

  if (! m.build_from_smiles(smiles))
  {
    std::cerr << "build_query_from_smiles:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);

  T * q = new T;

  if (! q->create_from_molecule(m, mqs))
  {
    std::cerr << "build_query_from_smiles:invalid molecule?? '" << smiles << "'\n";
    delete q;
    return 0;
  }

  queries.add(q);

  return 1;
}

template <typename T>
int
queries_from_file_of_molecules(MDL_Molecule & m,
                               Molecule_to_Query_Specifications & mqs,
                               resizable_array_p<T> & queries,
                               int verbose)
{
  T * q = new T;

  if (! q->create_from_molecule(m, mqs))
  {
    std::cerr << "queries_from_file_of_molecules:cannot create query from '" << m.name() << "'\n";
    delete q;
    return 0;
  }

  queries.add(q);

  return 1;
}

template <typename T>
int
queries_from_file_of_molecules(data_source_and_type<MDL_Molecule> & input,
                               Molecule_to_Query_Specifications & mqs,
                               resizable_array_p<T> & queries,
                               int verbose)
{
  set_input_aromatic_structures(1);

  MDL_Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! queries_from_file_of_molecules(*m, mqs, queries, verbose))
    {
      std::cerr << "Cannot create query from '" << m->name() << "'\n";
      return 0;
    }

    if (verbose > 1) {
      std::cerr << "Created query from '" << m->name() << "'\n";
    }
  }

  if (verbose)
    std::cerr << "Read " << queries.number_elements() << " queries\n";

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
    std::cerr << "cannot open '" << fname << "'\n";
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
  FileType input_type = discern_file_type_from_name(fname);

  if (input_type == FILE_TYPE_INVALID)
    input_type = FILE_TYPE_SMI;

// Don't follow any seeking or such from the command line

  const off_t o = moleculeio::seek_to_from_command_line();

  moleculeio::set_seek_to(static_cast<off_t>(0));

  data_source_and_type<MDL_Molecule> input(input_type, fname);
  if (! input.good())
  {
    std::cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  int rc = queries_from_file_of_molecules(input, mqs, queries, verbose);

  moleculeio::set_seek_to(o);

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

    std::cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      std::cerr << "INvalid molecule to query directives '" << directives << "'\n";
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
    std::cerr << "Created query from '" << m.name() << "'\n";

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

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! query_from_ISIS_query_file(*m, mqs, queries, verbose))
    {
      std::cerr << "queries_from_ISIS_query_file:cannot process '" << m->name() << "'\n";
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
  data_source_and_type<MDL_Molecule> input(FILE_TYPE_MDL, fname);

  if (! input.good())
  {
    std::cerr << "queries_from_ISIS_query_file:cannot open '" << fname << "'\n";
    return 0;
  }

  input.seekg(0);    // do not allow any seeking from the command line

  int rc = queries_from_ISIS_query_file(input, mqs, queries, verbose);

  if (0 == rc)
    return 0;

  if (verbose) {
    std::cerr << "Created " << queries.number_elements() << " queries from '" << fname << "'\n";
  }

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

    std::cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      std::cerr << "INvalid molecule to query directives '" << directives << "'\n";
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

    std::cerr << "Split into '" << fname2 << "' and '" << directives << "'\n";

    if (! mqs.parse_directives(directives))
    {
      std::cerr << "INvalid molecule to query directives '" << directives << "'\n";
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
      //std::cerr << "Creating query from smiles '" << buffer << "'\n";

      if (! q->create_from_smiles(buffer))
      {
        std::cerr << "smarts_or_smiles_from_file: cannot parse '" << buffer << "'\n";
        delete q;
        return 0;
      }      
    }
    else
    {
      //std::cerr << "Creating query from smarts '" << buffer << "'\n";
      if (! q->create_from_smarts(buffer))
      {
        std::cerr << "smarts_or_smiles_from_file: cannot parse '" << buffer << "'\n";
        delete q;
        return 0;
      }
    }

    queries.add(q);
    //std::cerr << "Created query from '" << buffer << "'\n";
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
    std::cerr << "smarts_from_file: cannot open '" << fname << "'\n";
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
    std::cerr << "smiles_from_file: cannot open '" << fname << "'\n";
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
    std::cerr << "Invalid smarts 'SMARTS:" << buffer << "'\n";
    delete tmp;
    return 0;
  }

  if (verbose)
    std::cerr << "Created query '" << tmp->comment() << "' from SMARTS:" << buffer << '\n';

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
      std::cerr << "process_queries: cannot build query from '" << msi << "'\n";
      return 0;
    }

    assert (tmp->ok());
    queries.add(tmp);

    if (verbose)
      std::cerr << "Created query " << (queries.number_elements() - 1) << " '" << tmp->comment() << "'\n";

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
    std::cerr << "read_one_or_more_queries_from_file::cannot open '" << fname << "'\n";
    return 0;
  }

  int rc = read_one_or_more_queries_from_file(queries, input, verbose);

  if (verbose)
    std::cerr << "Read " << rc << " queries from '" << fname << "'\n";

  return rc;
}

// A record in a file of queries starts with "PROTO:". Build the
// file name, read the query and append to `queries`.
template <typename T>
int
ReadProtoQuery(const IWString& fname,
               const IWString& directory_path,
               resizable_array_p<T>& queries, 
               int verbose) {
  assert(fname.starts_with("PROTO:"));

  const_IWSubstring myfname(fname);
  myfname.remove_leading_chars(6);

  IWString pathname;
  if (myfname.starts_with('/'))
    pathname=myfname;
  else if (directory_path.length())
    pathname = directory_path + myfname;
  else
    pathname = myfname;

  T * q = ReadProtoQueryFile<T>(pathname);
  if (nullptr == q) {
    std::cerr << "ReadProtoQuery::cannot read proto file '" << pathname << "'\n";
    return 0;
  }

  queries.add(q);

  return 1;
}

// Read a single record from a file containing queries.
template <typename T>
int
file_record_is_file(resizable_array_p<T> & queries,
                     const IWString & directory_path,
                     IWString & buffer,
                     int verbose)
{
  if (buffer.starts_with("PROTO:")) {
    return ReadProtoQuery(buffer, directory_path, queries, verbose);
  }

  IWString fname;
  if (! buffer.word(0, fname)) {
    std::cerr << "file_record_is_file: cannot get first word from '" << buffer << "'\n";
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
    std::cerr << "Queries_from_file: cannot read file '" << fname << "'\n";
    delete tmp;
    return 0;
  }

  if (verbose)
    std::cerr << "Created query '" << tmp->comment() << "' from '" << pathname << "'\n";

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
    if (buffer.empty())
      continue;

    if ('#' == buffer[0])
      continue;

    int rc_this_record;

    if (buffer.starts_with("SMARTS:"))
      rc_this_record = file_record_is_smarts(queries, buffer, verbose);
    else
      rc_this_record = file_record_is_file(queries, directory_path, buffer, verbose);

    if (0 == rc_this_record) {
      std::cerr << "Queries_from_file: fatal error on line " << input.lines_read() << '\n';
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
    std::cerr << "Cannot open file '" << fname << "'\n";
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
    std::cerr << "Cannot open file '" << fname << "'\n";
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

  const int rc = queries_from_file(input, queries, directory_path, verbose);
  if (rc == 0) {
    return 0;
  }
  if (verbose) {
    std::cerr << "After processing '" << fname << "' have " << queries.size() << " queries\n";
  }
  return rc;
}

template <typename T>
int
ReadTextProtoQuery(const const_IWSubstring& token, 
                   resizable_array_p<T>& queries) {
  SubstructureSearch::SubstructureQuery proto;

  const std::string s(token.data(), token.length());
  if (!google::protobuf::TextFormat::ParseFromString(s, &proto))
  {
    std::cerr << "ReadTextProtoQuery::cannot parse proto string '" << token << "'\n";
    return 0;
  }

  T * query = new T();
  if (! query->ConstructFromProto(proto)) {
    std::cerr << "ReadTextProtoQuery::cannot parse proto data '" << token << "'\n";
    delete query;
    return 0;
  }

  return queries.add(query);
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
      std::cerr << "process_files_of_queries: could not read queries from file '" <<
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
process_cmdline_token(char option,
                      const const_IWSubstring & token,
                      resizable_array_p<T> & queries,
                      int verbose)
{
  const_IWSubstring mytoken(token);

//std::cerr << "Examining token '" << mytoken << "'\n";
  if (mytoken.starts_with("F:") || mytoken.starts_with("Q:"))
  {
    mytoken.remove_leading_chars(2);

    if (0 == mytoken.length())
    {
      std::cerr << "Must follow S: specification with file name of queries\n";
      return 0;
    }

    if (! queries_from_file(mytoken, queries, 1, verbose))   // queries always in same directory as controlling file
    {
      std::cerr << "process_queries: cannot read queries from file specifier 'F:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("S:"))
  {
    mytoken.remove_leading_chars(2);

    if (0 == mytoken.length())
    {
      std::cerr << "Must follow S: specification with file name of queries\n";
      return 0;
    }

    if (! smarts_from_file(mytoken, queries, verbose))
    {
      std::cerr << "process_queries::cannot read smarts from file of smarts specifier 'S:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("M:"))
  {
    mytoken.remove_leading_chars(2);

    if (! mytoken.length())
    {
      std::cerr << "Must follow M: specification with file name of molecules\n";
      return 0;
    }

    if (! queries_from_file_of_molecules(mytoken, queries, verbose))
    {
      std::cerr << "process_queries::cannot read queries from file of molecules specifier 'M:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("smiles:"))
  {
    mytoken.remove_leading_chars(7);

    if (0 == mytoken.length())
    {
      std::cerr << "Must follow smiles: specification with smiles\n";
      return 0;
    }

    if (! build_query_from_smiles(mytoken, queries, verbose))
    {
      std::cerr << "process_queries:cannot build query from 'smiles:" << mytoken << "'\n";
      return 0;
    }
  }
  else if (mytoken.starts_with("I:"))
  {
    mytoken.remove_leading_chars(2);
    if (0 == mytoken.length())
    {
      std::cerr << "Must follow I: specification with query file\n";
      return 0;
    }

    if (! queries_from_ISIS_query_file(mytoken, queries, verbose))
    {
      std::cerr << "process_queries::cannot read queries from isis query 'I:" << mytoken << "'\n";
      return 0;
    }
  }
/*else if (mytoken.starts_with("ISIS:"))
  {
    mytoken.remove_leading_chars(5);
    if (0 == mytoken.length())
    {
      std::cerr << "Must follow ISIS: specification with file name of queries\n";
      return 0;
    }

    if (! queries_from_file_of_isis_queries(mytoken, queries, verbose))
    {
      std::cerr << "process_queries::cannot read queries from file of isis queries 'ISIS:" << mytoken << "'\n";
      return 0;
    }
  }*/
  else if ("help" == mytoken)
  {
    std::cerr << "The following query specifications are recognised\n";
    std::cerr << " -" << option <<" SMARTS:smarts          smarts (use quotes to hide special characters)\n";
    std::cerr << " -" << option <<" S:file                 file of smarts queries\n";
    std::cerr << " -" << option <<" Q:file                 file of query object queries (also F: recognised)\n";
    std::cerr << " -" << option <<" M:file                 file of molecules that will be converted to query objects\n";
    std::cerr << " -" << option <<" I:file                 an ISIS query file\n";
    std::cerr << " -" << option <<" proto:proto            inline proto 'proto:query{nrings:3}'\n";
    std::cerr << " -" << option <<" PROTO:file             a query in text proto form\n";
    std::cerr << " -" << option <<" PROTOFILE:file         a file containing names of proto query files\n";
    std::cerr << " -" << option <<" TFPROTO:file           TFDataRecord file containing serialized protos\n";
    std::cerr << " -" << option <<" file                   single query file\n";
    ::exit (0);
  }
  else if (mytoken.starts_with("SMARTS:"))
  {
    mytoken.remove_leading_chars(7);

    T * q = new T;
    if (! q->create_from_smarts(mytoken))
    {
      std::cerr << "process_queries::invalid smarts '" << mytoken << "'\n";
      delete q;
      return 0;
    }

    queries.add(q);
  }
  else if (mytoken.starts_with("PROTO:"))
  {
    mytoken.remove_leading_chars(6);
    T * q = ReadProtoQueryFile<T>(mytoken);
    if (nullptr == q) {
      std::cerr << "process_queries::cannot read proto file '" << mytoken << "'\n";
      return 0;
    }
    queries.add(q);
  } else if (mytoken.starts_with("PROTOFILE:")) {
    mytoken.remove_leading_chars(10);
    if (!ReadFileOfProtoQueries(mytoken, queries)) 
    {
      std::cerr << "process_queries:cannot read file of proto files '" << mytoken << "'\n";
      return 0;
    }
  } else if (mytoken.starts_with("TFPROTO:")) {
    mytoken.remove_leading_chars(8);
    if (! ReadQueriesFromTFDataRecordProtoFile(mytoken, queries)) {
      std::cerr << "process_queries::cannot read TFdataREcord proto file '" << mytoken << "'\n";
      return 0;
    }
  } else if (mytoken.starts_with("proto:")) {
    // for example -q 'proto:query{min_natoms: 4}'
    mytoken.remove_leading_chars(6);
    if (! ReadTextProtoQuery(mytoken, queries)) {
      std::cerr << "process_queries:cannot read text proto option '" << mytoken << "'\n";
      return 0;
    }
  } else if (mytoken.starts_with("MPROTO:")) {
    mytoken.remove_leading_chars(7);
    if (! MultipleQueriesFromTextProto(mytoken, verbose, queries)) {
      std::cerr << "process_queries:cannot read multiple text proto option '" << mytoken << "'\n";
      return 0;
    }
  } else {
    if (! read_one_or_more_queries_from_file(queries, mytoken, verbose))
    {
      std::cerr << "process_queries::cannot read query/queries from '" << mytoken << "'\n";
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
      std::cerr << "Cannot process -" << option << " option '" << c << "'\n";
      return 0;
    }
  }

  return queries.number_elements();
}

#endif  // MOLECULE_LIB_RWSUBSTRUCTURE_H_
