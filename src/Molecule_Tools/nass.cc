#include <stdlib.h>
#include <memory>
using namespace std;

#include "cmdline.h"
#include "msi_object.h"
#include "misc.h"

#include "nass.h"
#include "iw_stl_hash_map.h"
#include "iwstring_data_source.h"
#include "target.h"

NA_Substructure_Query::NA_Substructure_Query()
{
}

int
NA_Substructure_Query::read (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "NA_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  int rc = Substructure_Query::read(input);

  if (0 == rc)
    cerr << "NA_Substructure_Query::read: cannot read from '" << fname << "'\n";

  return rc;
}

#define NASS_NAME_OF_AVAILABLE_ATTRIBUTE "makes_available"
#define NASS_NAME_OF_NEEDS_ATTRIBUTE "needs"

int
NA_Substructure_Query::read (iwstring_data_source & input)
{
  msi_object msi;

  if (! msi.read(input))
  {
    cerr << "NA_Substructure_Query::read: cannot fetch msi object\n";
    return 0;
  }

  if (! Substructure_Query::construct_from_msi_object(msi))
  {
    cerr << "NA_Substructure_Query::read: cannot parse msi object\n";
    cerr << msi;
    return 0;
  }

  (void) msi.string_value_for_attribute(NASS_NAME_OF_AVAILABLE_ATTRIBUTE, _makes_available);

  int nn = msi.attribute_count(NASS_NAME_OF_NEEDS_ATTRIBUTE);
  if (0 == nn)
    return 1;

  _string_needs.resize(nn);

  const msi_attribute * att;
  int i = 0;
  while (NULL != (att = msi.attribute(NASS_NAME_OF_NEEDS_ATTRIBUTE, i++)))
  {
    IWString * n = new IWString;
    att->value(*n);

    _string_needs.add(n);
  }

  return _string_needs.number_elements();
}

int
NA_Substructure_Query::also_needs (const const_IWSubstring & n)
{
  IWString * tmp = new IWString(n);

  _string_needs.add(tmp);

  return 1;
}

#define NASS_UNDETERMINED -1

/*
  Mostly during testing, we need to distinguish between those queries
  which were actually evaluated and those which were implied as zero
  hits because preconditions were not met
*/

#define NASS_NOT_EVALUATED -2
#define NASS_PRECONDITIONS_NOT_MET -1

Set_of_NA_Substructure_Query::Set_of_NA_Substructure_Query()
{
  _queries_not_needing_anything = NASS_UNDETERMINED;

  _sresults = NULL;

  _break_at_first_match = 0;

  _break_at_first_non_match = 0;

  _ignore_precondition_matches_for_break = 0;

  _searches_performed = 0;

  _optimisation_level = 0;
}

Set_of_NA_Substructure_Query::~Set_of_NA_Substructure_Query()
{
}

int
Set_of_NA_Substructure_Query::debug_print (ostream & os) const
{
  os << "Set of N/A queries with " << _number_elements << " queries\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const NA_Substructure_Query & q = *(_things[i]);

    os << ' ' << i << " '" << q.comment() << "' ";

    if (q.makes_available().length())
      os << "makes available '" << q.makes_available() << "' ";

    const resizable_array_p<IWString> & needs = q.string_needs();

    if (needs.number_elements())
    {
      os << "NEEDS:";
      for (int j = 0; j < needs.number_elements(); j++)
      {
        const IWString & n = *(needs[j]);
  
        os << ' ' << n;
      }
    }

    os << endl;
  }

  return os.good();
}

void
Set_of_NA_Substructure_Query::_count_queries_not_needing_anything()
{
  assert (_number_elements > 0);

  _queries_not_needing_anything = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == _things[i]->string_needs().number_elements())
      _queries_not_needing_anything++;
  }

  return;
}

int
Set_of_NA_Substructure_Query::build_from_command_line (Command_Line & cl,
                                                       char cflag,
                                                       int verbose)
{
  int i = 0;
  const_IWSubstring c;
  while (cl.value(cflag, c, i++))
  {
    if (! build_from_option(c, verbose))
    {
      cerr << "Set_of_NA_Substructure_Query::build_from_command_line: cannot parse '" << c << "'\n";
      return 0;
    }
  }

  return _number_elements;
}

int
Set_of_NA_Substructure_Query::build_from_option (const const_IWSubstring & c, int verbose)
{
  if (c.starts_with("F:"))
    return _build_from_query_file(c, verbose);

  if (c.starts_with("S:"))
    return _build_from_smarts_file(c, verbose);

  if (c.starts_with("SMARTS:"))
    return _build_from_smarts(c, verbose);

// Must be a query file

  if (_elements_allocated == _number_elements)
    resize(_number_elements + 100);

  NA_Substructure_Query * q = new NA_Substructure_Query;

// A query file record may look like
// fname   NEED=foo AVAIL=bar
// so extract the first token as the file name

  const_IWSubstring fname = c;
  fname.truncate_at_first(' ');
    
  if (! q->read(fname))
  {
    cerr << "Set_of_NA_Substructure_Query::build_from_option: cannot read query file '" << fname << "'\n";
    delete q;
    return 0;
  }

  if (! _get_need_avail(q, c))
  {
    cerr << "Set_of_NA_Substructure_Query::build_from_option: cannot interpret '" << c << "'\n";
    delete q;
    return 0;
  }

  add(q);

  return _number_elements;
}

int
Set_of_NA_Substructure_Query::_build_from_query_file (const const_IWSubstring & c, int verbose)
{
  IWString fname(c);
  assert (fname.starts_with("F:"));

  fname.remove_leading_chars(2);

  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_query_file: cannot open '" << fname << "'\n";
    return 0;
  }

// we need to form the directory name

  int i = fname.rindex('/');

  if (i < 0)
    fname = "./";
  else
    fname.iwtruncate(i + 1);

  return _build_from_query_file(input, fname, verbose);
}

int
Set_of_NA_Substructure_Query::_build_from_query_file (iwstring_data_source & input,
                                                      const IWString & directory,
                                                      int verbose)
{
  int nr = input.records_remaining();

  if (_number_elements + nr > _elements_allocated)
    resize(_number_elements + nr);

  input.set_strip_trailing_blanks(1);
  input.set_strip_leading_blanks(1);
  input.set_skip_blank_lines(1);

  int rc = 0;     // the number of queries we read

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    IWString tmp;    // we may need to pre-pend the directory name

    if (buffer.starts_with('/'))     // looks like a full path name
      tmp = buffer;
    else
    {
      tmp.resize(directory.length() + buffer.length());

      tmp = directory;
      tmp += buffer;
    }

    NA_Substructure_Query * q = new NA_Substructure_Query;

    const_IWSubstring fname = tmp;
    fname.truncate_at_first(' ');    // first token is the file name

    if (! q->read(fname))
    {
      cerr << "Set_of_NA_Substructure_Query::_build_from_query_file: cannot read query from '" << fname << "'\n";
      delete q;
      return 0;
    }

    if (! _get_need_avail(q, tmp))
    {
      cerr << "Set_of_NA_Substructure_Query::build_from_option: cannot interpret '" << tmp << "'\n";
      delete q;
      return 0;
    }

    add(q);
    rc++;
  }

  return rc;
}

int
Set_of_NA_Substructure_Query::_build_from_smarts_file (const const_IWSubstring & c, int verbose)
{
  IWString fname(c);

  assert (fname.starts_with("S:"));

  fname.remove_leading_chars(2);

  if (0 == fname.length())
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts_file:no file name\n";
    return 0;
  }

  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts_file: cannot open '" << fname << "'\n";
    return 0;
  }

  return _build_from_smarts_file(input, verbose);
}

int
Set_of_NA_Substructure_Query::_build_from_smarts_file (iwstring_data_source & input, int verbose)
{
  int nr = input.records_remaining();

  if (_number_elements + nr < _elements_allocated)
    resize(_number_elements + nr);

  input.set_skip_blank_lines(1);
  input.set_strip_leading_blanks(1);
  input.set_strip_trailing_blanks(1);

  int rc = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if(buffer.starts_with('#'))
      continue;

    NA_Substructure_Query * q = new NA_Substructure_Query;

    if (! q->create_from_smarts(buffer))
    {
      cerr << "Set_of_NA_Substructure_Query::_build_from_smarts_file: cannot parse '" << buffer << "'\n";
      delete q;
      return 0;
    }

    if (! _get_need_avail(q, buffer))
    {
      cerr << "Set_of_NA_Substructure_Query::_build_from_smarts_file: cannot parse need/avail '" << buffer << "'\n";
      delete q;
      return 0;
    }

    add(q);

    if (verbose > 1)
      cerr << "Added smarts query '" << q->comment() << "'\n";

    rc++;
  }

  return rc;
}


int
Set_of_NA_Substructure_Query::_build_from_smarts (const const_IWSubstring & c, int verbose)
{
  assert (c.starts_with("SMARTS:"));

  const_IWSubstring smarts(c);

  smarts.remove_leading_chars(7);

  if (_number_elements == _elements_allocated)
    resize(_number_elements + 100);

  NA_Substructure_Query * q = new NA_Substructure_Query;

  if (! q->create_from_smarts(smarts))
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts: cannot parse smarts '" << c << "'\n";
    delete q;
    return 0;
  }

  if (! _get_need_avail(q, c))
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts: cannot process need/avail '" << c << "'\n";
    delete q;
    return 0;
  }

  add(q);

  return 1;
}

int
Set_of_NA_Substructure_Query::_get_need_avail (NA_Substructure_Query * q,
                                               const const_IWSubstring & buffer)
{
  if (! _parse_need_token(q, buffer))
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts: cannot parse need token '" << buffer << "'\n";
    return 0;
  }

  if (! _parse_avail_token(q, buffer))
  {
    cerr << "Set_of_NA_Substructure_Query::_build_from_smarts: cannot parse avail token '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
get_token (const const_IWSubstring & buffer, 
           const char * zprefix,
           IWString & zresult)
{
  zresult = "";

  int i = buffer.index(zprefix);

  if (i < 0)
    return 1;

  buffer.from_to(i + ::strlen(zprefix), zresult);

  if (zresult.starts_with('='))
    zresult.remove_leading_chars(1);

  zresult.truncate_at_first(' ');    // we extract the first token only

  return zresult.length();      // we must have some kind of valid token
}

int
Set_of_NA_Substructure_Query::_parse_avail_token (NA_Substructure_Query * q,
                                                  const const_IWSubstring & buffer)
{
  IWString avail;

  if (! get_token(buffer, NASS_AVAIL_TOKEN, avail))
  {
    cerr << "Set_of_NA_Substructure_Query::_parse_avail_token: cannot parse avail token '" << buffer << "'\n";
    return 0;
  }

  if (0 != q->makes_available().length())
    return 1;

  q->set_makes_available(avail);

  return 1;
}
               

int
Set_of_NA_Substructure_Query::_parse_need_token (NA_Substructure_Query * q,
                                                 const const_IWSubstring & buffer)
{
  const_IWSubstring token;
  int i = 0;
  while (buffer.nextword(token, i))
  {
    if (! token.starts_with(NASS_NEED_TOKEN))
      continue;

    token.remove_leading_chars(int(::strlen(NASS_NEED_TOKEN)));

    if (0 == token.length())
    {
      cerr << "Set_of_NA_Substructure_Query::_parse_need_token: empty need value\n";
      return 0;
    }

    if (! q->also_needs(token))
    {
      cerr << "Set_of_NA_Substructure_Query::_parse_need_token: yipes also_needs failed\n";
      return 0;
    }
  }

  return 1;
}

int
Set_of_NA_Substructure_Query::establish_cross_references()
{
  IW_STL_Hash_Map_int avail_hash;

//avail_hash.resize(_number_elements * 2);

  for (int i = 0; i < _number_elements; i++)
  {
    NA_Substructure_Query & q = *(_things[i]);
    const IWString & avail = q.makes_available();

    if (0 == avail.length())
      continue;

    if (avail_hash.contains(avail))
    {
      cerr << "Set_of_NA_Substructure_Query::establish_cross_references: '" << avail << "' multiply defined availablity\n";
      cerr << "Queries " << avail_hash[avail] << " and " << i << " both make this available\n";
      return 0;
    }

    avail_hash[avail] = i;
  }

  if (avail_hash.empty())
  {
    cerr << "Set_of_NA_Substructure_Query::establish_cross_references: no available queries!\n";
  }

//#define ECHO_AVAIL_HASH
#ifdef ECHO_AVAIL_HASH
  cerr << avail_hash.size() << " items in hash\n";
  IW_STL_Hash_Map_int::const_iterator m;
  for (m = avail_hash.begin(); m != avail_hash.end(); ++m)
  {
    cerr << "Hash item " << (*m).first << " value " << (*m).second << endl;
  }
#endif

// Now make sure that everything which is needed is in fact available.
// Convert the string "need" values to int's

  for (int i = 0; i < _number_elements; i++)
  {
    const resizable_array_p<IWString> & string_need = _things[i]->string_needs();
    resizable_array<int> & int_need = _things[i]->needs();    // we will fill this array with the corresponding INT

    int nn = string_need.number_elements();

    int_need.resize(nn);

    for (int j = 0; j < nn; j++)
    {
      const IWString & n = *(string_need[j]);

      IW_STL_Hash_Map_int::const_iterator f = avail_hash.find(n);

      if (f == avail_hash.end())
      {
        cerr << "Set_of_NA_Substructure_Query::establish_cross_references: query " << i << " needs '" << n << "' which has not been made available\n";
        return 0;
      }
      
      int_need.add((*f).second);
    }
  }

// The final check is to look for circular dependencies

  int * tmp = new int[_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);

  return _find_circular_dependencies(tmp);
}

/*

*/

int
Set_of_NA_Substructure_Query::_find_circular_dependency (NA_Substructure_Query * q,
                                                         int * already_hit)
{
  const resizable_array<int> & needs = q->needs();

  for (int i = 0; i < needs.number_elements(); i++)
  {
    int j = needs[i];

    if (already_hit[j])
    {
      cerr << "Set_of_NA_Substructure_Query::_find_circular_dependency: query '" << q->comment() << " needs " << j << " but already set\n";
      return 1;
    }

    already_hit[j] = 1;

    if (_find_circular_dependency(_things[j], already_hit))
      return 1;

    already_hit[j] = 0;
  }

  return 0;    // no circular dependency found
}

int
Set_of_NA_Substructure_Query::_find_circular_dependencies (int * already_hit)
{
  for (int i = 0; i < _number_elements; i++)
  {
    const NA_Substructure_Query & q = *(_things[i]);

    const resizable_array<int> & needs = q.needs();

    if (0 == needs.number_elements())
      continue;
      
    set_vector(already_hit, _number_elements, 0);

    already_hit[i] = 1;

    if (_find_circular_dependency(_things[i], already_hit))
    {
      cerr << "Set_of_NA_Substructure_Query::_find_circular_dependencies: circular dependency involving query " << i << ", '" << q.comment() << "'\n";
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_SUBSTRUCTURE_SEARCH

/*
  Record the number of hits for each query
*/

int
Set_of_NA_Substructure_Query::substructure_search (Molecule & m,
                                                   int * result)
{
  assert (NULL != result);

  if (NASS_UNDETERMINED == _queries_not_needing_anything)
  {
    _count_queries_not_needing_anything();
    if (! establish_cross_references())
      return 0;
  }

  if (NULL == _sresults)
  {
    _sresults = new Substructure_Results[_number_elements];
    for (int i = 0; i < _number_elements; i++)
    {
      _sresults[i].set_save_query_atoms_matched(1);    // must save matched atoms in case any queries need unique matches
    }
  }

  set_vector(result, _number_elements, NASS_NOT_EVALUATED);

  for (int i = 0; i < _number_elements; i++)
  {
    _sresults[i].reset();
  }

  int rc = _substructure_search(m, result);

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " Query " << i;
    if (NASS_PRECONDITIONS_NOT_MET == result[i])
      cerr << " not evaluated\n";
    else if (result[i] >= 0)
      cerr << ' ' << result[i] << " hits\n";
    else
      cerr << " result = " << result[i] << endl;
  }
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if (NASS_PRECONDITIONS_NOT_MET == result[i])
    {
      result[i] = 0;
      _optimisation_level++;
    }
    else if (NASS_NOT_EVALUATED == result[i] && (_break_at_first_match || _break_at_first_non_match))    // hmmm, what should we really do here?
      ;
    else if (result[i] < 0)
    {
      cerr << "Set_of_NA_Substructure_Query::substructure_search: huh, query " << i << " not set!! result " << result[i] << endl;
      rc = 0;
    }
  }

  return rc;
}

int
Set_of_NA_Substructure_Query::_substructure_search (Molecule & m, int * result)
{
  _searches_performed++;

  Molecule_to_Match target(&m);

// First do all the queries which don't need anything


#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << _queries_not_needing_anything << " queries with no preconditions\n";
#endif

  if (_queries_not_needing_anything)     // these will always be evaluated
  {
    for (int i = 0; i < _number_elements; i++)
    {
      NA_Substructure_Query & q = *(_things[i]);

      if (q.makes_available().length())    // is a precondition, may never be needed
        continue;

      if (q.needs().number_elements())
        continue;
  
      int nhits = q.substructure_search(target, _sresults[i]);

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
      cerr << nhits << " hits to query " << i << endl;
#endif

      result[i] = nhits;

      if (_ignore_precondition_matches_for_break && q.makes_available().length())
        continue;

      if (_break_at_first_match && nhits)
        return 1;

      if (_break_at_first_non_match && 0 == nhits)
        return 1;
    }
  }

// Now we need to work on the queries which need other queries

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Checking dependent queries\n";
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if (result[i] >= 0)     // already evaluated
      continue;

    NA_Substructure_Query * q = _things[i];

    if (! _establish_preconditions(q, target, result))
      result[i] = NASS_PRECONDITIONS_NOT_MET;
    else
      result[i] = q->substructure_search(target, _sresults[i]);

//  We need to decide whether to return based on one of the _break conditions

    if (_ignore_precondition_matches_for_break && q->makes_available().length())
      continue;

    if (result[i] > 0)    // got a match
    {
      if (_break_at_first_match)
        return 1;
    }
    else    // no match
    {
      if (_break_at_first_non_match)
        return 1;
    }
  }

  return 1;
}

int
Set_of_NA_Substructure_Query::_establish_preconditions (NA_Substructure_Query * q,
                                                        Molecule_to_Match & target,
                                                        int * result)
{
#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Establishing preconditions for '" << q->comment() << "'\n";
#endif

  const resizable_array<int> & needs = q->needs();

  Substructure_Results sresults;

  for (int i = 0; i < needs.number_elements(); i++)
  {
    int j = needs[i];

    if (result[j] > 0)     // already computed and hit
      continue;

    if (0 == result[j] || NASS_PRECONDITIONS_NOT_MET == result[j])     // already computed and no hits
      return 0;

//  Query J has not yet been evaluated.

    if (_things[j]->has_preconditions())
    {
      if (! _establish_preconditions(_things[j], target, result))
      {
        result[j] = NASS_PRECONDITIONS_NOT_MET;
        return 0;
      }
    }

    result[j] = _things[j]->substructure_search(target, _sresults[j]);
    if (0 == result[j])
      return 0;
  }

  return 1;
}
