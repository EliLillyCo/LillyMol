#include <stdlib.h>

#include "Foundational/iwmisc/iwre2.h"
#include "iw_tdt.h"

using std::cerr;
using std::endl;

IW_TDT::IW_TDT()
{
};

/*
  We can turn off some checks
*/

static int iwtdt_careful_mode = 0;

void
set_iwtdt_careful_mode (int c)
{
  iwtdt_careful_mode = c;
}

/*
  It is optional whether or not we store newlines in the TDT.
  But, if they are NOT stored, writing will come out in Daylight dump format!

  The value must be 0 or 1, 
*/

static int _include_newlines_in_tdt = 1;

void
set_include_newlines_in_tdt (int i)
{
  if (i)
    _include_newlines_in_tdt = 1;
  else
    _include_newlines_in_tdt = 0;
}

int 
include_newlines_in_tdt()
{
  return _include_newlines_in_tdt;
}

static int
valid_tdt_form (const const_IWSubstring & buffer)
{
  if (! buffer.ends_with('>'))
  {
    cerr << "valid_tdt_form: TDT's must end with a > character\n";
    return 0;
  }

  int openangle = buffer.index('<');
  if (openangle < 1)
  {
    cerr << "valid_tdt_form: TDT's must start with a valid tag\n";
    return 0;
  }

  if (openangle == buffer.length() - 2)
    cerr << "valid_tdt_form: warning, TDT with tag but no data\n";

  return 1;
}

int
IW_TDT::ok() const
{
  if (0 == _zdata.length() && 0 == _end.number_elements())
    return 1;

  if (_zdata.length() && _end.number_elements())
  {
    if (_include_newlines_in_tdt)
    {
      if (_zdata.ends_with(">\n|\n"))
        return 1;

      cerr << "TDT data does not end in newline VBAR\n";
      cerr << "***'" << _zdata << "'***\n";
      return 0;
    }
    else if (_zdata.ends_with(">|"))
      return 1;
    else
      return 0;
  }

  return 0;
}

int
IW_TDT::next (iwstring_data_source & input)
{
  input.set_skip_blank_lines(1);

  _zdata.resize_keep_storage(0);
  _end.resize_keep_storage(0);

//_offset.set (input.tellg());

  int got_vbar = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
//  cerr << "Read ---'" << buffer << "'---\n";
    if (_zdata.length())
      _end.add(_zdata.length());
      
    _zdata += buffer;

//  cerr << "Data now '" << _zdata << "'\n";

    if (_include_newlines_in_tdt)
      _zdata += '\n';

    if ('|' == buffer)
    {
      got_vbar = 1;
      break;
    }

    if (iwtdt_careful_mode && ! valid_tdt_form(buffer))
    {
      cerr << "Invalid tdt form '" << buffer << "'\n";
      return 0;
    }
  }

  if (_zdata.length() > 0 && ! got_vbar)
  {
    cerr << "IW_TDT::next:improperly terminated TDT, line " << input.lines_read() << endl;
    return 0;
  }

  return _end.number_elements();
}

int
IW_TDT::next_dataitem (const_IWSubstring & zresult,
                       int & ptr) const
{
  assert (ok());

  if (ptr < 0)
    return 0;

  if (ptr >= _end.number_elements())
    return 0;

  if (0 == ptr)
    zresult.set(_zdata.rawchars(), _end[0]);
  else
    zresult.set(_zdata.rawchars() + _end[ptr - 1], _end[ptr] - _end[ptr - 1]);

  ptr++;

  return 1;
}

int
IW_TDT::next_dataitem_value (const_IWSubstring & ztag,
                             const_IWSubstring & zdata,
                             int & ptr) const
{
  const_IWSubstring tmp;
  if (! next_dataitem(tmp, ptr))
    return 0;

  int open_angle_bracket = tmp.index('<');

  assert (open_angle_bracket > 0);
  assert (open_angle_bracket < tmp.length());

  tmp.from_to(0, open_angle_bracket - 1, ztag);
  tmp.from_to(open_angle_bracket + 1, tmp.length() - 3, zdata);

//cerr << "next_dataitem_value: tag '" << ztag << "' result '" << tmp << "'\n";

  return 1;
}

int
IW_TDT::do_write (std::ostream & output) const
{
  assert (ok());

  if (0 == _zdata.length())
  {
    cerr << "IW_TDT::do_write: write on empty TDT\n";
    return 0;
  }

  output << _zdata;

  return output.good();
}

template <typename S>
int
IW_TDT::_common_write_all_except_vbar(S & output) const
{
  assert (ok());

  if (0 == _zdata.length())
  {
    cerr << "IW_TDT::write_all_except_vbar:empty tdt\n";
    return 0;
  }

  int chars_to_write = _zdata.length() - 1;

  if (_include_newlines_in_tdt)
    chars_to_write--;

  output.write(_zdata.rawchars(), chars_to_write);

  return output.good();
}


template int IW_TDT::_common_write_all_except_vbar(std::ostream & output) const;
template int IW_TDT::_common_write_all_except_vbar(IWString_and_File_Descriptor & output) const;

int
IW_TDT::write_all_except_vbar (std::ostream & output) const
{
  return _common_write_all_except_vbar (output);
}

int
IW_TDT::write_all_except_vbar (IWString_and_File_Descriptor & output) const
{
  return _common_write_all_except_vbar (output);
}

int
IW_TDT::_find_index_in_end_array (const char * dkey,
                                   int klen,
                                   int which_one) const
{
  if (0 == _zdata.length())
    return -1;

//cerr << "IW_TDT::_find_index_in_end_array\n";

  const char * s = _zdata.rawchars();

  int istart = 0;

  int instances_found = 0;

  int last_one_found = -1;

  for (int i = 0; i < _end.number_elements(); i++)
  {
    int istop = _end[i];

//  cerr << "istart " << istart << " istop " << istop << " klen " << klen << endl;

    if (istop - istart - 2 < klen)   // must have at least <>
      ;
    else if (0 == strncmp(s + istart, dkey, klen))
    {
      if (which_one == instances_found)
        return i;

      last_one_found = i;

      instances_found++;
    }

    istart = istop;     // get ready for next one
  }

  if (-1 == which_one && last_one_found >= 0)
    return last_one_found;

  return -1;
}

std::ostream &
operator << (std::ostream & os, const IW_TDT & tdt)
{
  assert (tdt.ok());

  tdt.do_write(os);

  return os;
}

IWString &
operator << (IWString & output, const IW_TDT & tdt)
{
  output << tdt.rawdata();

  return output;
}

int
IW_TDT::remove_all(RE2 & rx)
{
  assert (ok());

  if (0 == _zdata.length())
    return 0;

  IWString new_data;
  resizable_array<int> new_end;
  
  int i = 0;
  const_IWSubstring token;

  while (next_dataitem(token, i))
  {
//  cerr << "Does '" << token << "' match '" << rx.source() << "'\n";

    if (! iwre2::RE2PartialMatch(token, rx))  // we remove all the ones that match rx
    {
      if (new_data.length())
        new_end.add(new_data.length());

      new_data += token;
      if (_include_newlines_in_tdt && ! new_data.ends_with('\n'))
        new_data += '\n';
    }
  }

  int rc = _end.number_elements() - new_end.number_elements();

  if (0 == rc)    // nothing removed
    return 0;

  _zdata = new_data;

  _zdata += '|';
  if (_include_newlines_in_tdt)
    _zdata += '\n';

  new_end.add(_zdata.length());

//cerr << "End array contains " << new_end.number_elements() << endl;

//_end.resize_keep_storage (0);

  _end = new_end;
//for (int i = 0; i < new_end.number_elements(); i++)
//{
//  _end.add (new_end[i]);
//}

  return rc;
}

int
IW_TDT::dataitem_value (RE2 & rx,
                        const_IWSubstring & dataitem,
                        const_IWSubstring & zresult,
                        int which_to_return) const
{
  if (! ok())
  {
    cerr << "Not OK\n";
    cerr << _zdata;
  }

  assert (ok());

  if (0 == _zdata.length())
    return 0;

  int nfound = 0;

//#define DEBUG_DATAITEM_VALUE
#ifdef DEBUG_DATAITEM_VALUE
  cerr << "Checking " << _end.number_elements() << " items\n";
#endif

  int i = 0;
  const_IWSubstring zitem;

  while (next_dataitem(zitem, i))
  {
#ifdef DEBUG_DATAITEM_VALUE
    cerr << "What about '" << zitem << "'\n";
#endif

    const_IWSubstring tag(zitem);

    int openangle = tag.index('<');
    assert (openangle > 0);

    tag.iwtruncate (openangle);

#ifdef DEBUG_DATAITEM_VALUE
    cerr << "Does '" << rx.source() << "' match '" << tag << "'\n";
#endif

    if (iwre2::RE2PartialMatch(tag, rx))
      continue;

    if (nfound == which_to_return)
    {
      zitem.from_to(0, openangle - 1, dataitem);
      zitem.from_to(openangle + 1, zitem.length() - 2 - _include_newlines_in_tdt, zresult);

#ifdef DEBUG_DATAITEM_VALUE
      cerr << "Results '" << dataitem << "' and '" << zresult << "'\n";
#endif

      return 1;
    }

    nfound++;
  }

  return 0;
}

int
IW_TDT::index_of_dataitem (const const_IWSubstring & tag) const
{
  int i = 0;
  const_IWSubstring token;

  int rc = 0;    // we could use the value of the variable i, but it should remain an opaque iterator-type-thing

  while (next_dataitem(token, i))
  {
    if (token.starts_with(tag))
      return rc;

    rc++;
  }

  return -1;
}

int
IW_TDT::remove_item (int to_remove)
{
  assert (ok());

  if (0 == _zdata.length())
    return 0;

  if (! _end.ok_index(to_remove))
  {
    cerr << "IW_TDT::remove_item: invalid item to remove " << to_remove << " has " << _end.number_elements() << " components\n";
    return 0;
  }

  int characters_removed;

  if (0 == to_remove)
  {
    _zdata.remove_leading_chars(_end[0]);
    characters_removed = _end[0];
  }
  else
  {
    _zdata.remove_from_to(_end[to_remove - 1], _end[to_remove] - 1);
    characters_removed = _end[to_remove] - _end[to_remove - 1];
  }

  int ne = _end.number_elements() - 1;

  for (int i = to_remove; i < ne; i++)
  {
    _end[i] = _end[i + 1] - characters_removed;
  }

  _end.pop();

  return characters_removed;
}

int
IW_TDT::insert_before (int b4,
                       const const_IWSubstring & to_insert)
{
  assert (ok());
  assert (_end.ok_index(b4));

  cerr << "IW_TDT::insert_before:not implemented\n";

  return 0;
}

int
IW_TDT::remove_items_with_tag (const IWString & tag)
{
  if (tag.ends_with('<'))
    return remove_items_with_tag(tag.rawchars(), tag.length());

  IWString tmp(tag);
  tmp += '<';

  return remove_items_with_tag(tmp.rawchars(), tmp.length());
}

int
IW_TDT::remove_items_with_tag (const const_IWSubstring & tag)
{
  if (tag.ends_with('<'))
    return remove_items_with_tag(tag.rawchars(), tag.length());

  IWString tmp(tag);
  tmp += '<';

  return remove_items_with_tag(tmp.rawchars(), tmp.length());
}

int
IW_TDT::remove_items_with_tag (const char * s, int len)
{
  int rc = 0;

  for (int i = _end.number_elements() - 1; i >= 0; i--)
  {
    const char * d = _zdata.rawchars();

    if (i > 0)
      d += _end[i - 1];

    if (0 != strncmp(d, s, len))
      continue;

    remove_item(i);
    rc++;
  }

  return rc;
}
 
int
IW_TDT::item (int ndx,
              const_IWSubstring & v) const
{
//#define DEBUG_ITEM
#ifdef DEBUG_ITEM
  cerr << "Looking for item " << ndx << " data contains " << _zdata.length() << " bytes\n";
  for (int i = 0; i < _end.number_elements(); i++)
  {
    cerr << " item " << i << " ends at " << _end[i] << endl;
  }
#endif

  if (! _end.ok_index(ndx))
  {
    cerr << "IW_TDT::item: item " << ndx << " is invalid\n";
    return 0;
  }

#ifdef JAN2008_THIS_LOOKS_WRONG
  int istart = _end[ndx];

  int istop;
  if (ndx == _end.number_elements() - 1)
    istop = _zdata.length() - 1;
  else
    istop = _end[ndx + 1];
#else
  int istart;
  if (0 == ndx)
    istart = 0;
  else
    istart = _end[ndx - 1];

#endif

  _zdata.from_to(istart, _end[ndx] - 1, v);

#ifdef DEBUG_ITEM
  cerr << "Item " << ndx << " from " << istart << " to " << _end[ndx] << " got " << v << "'\n";
#endif

  return 1;
}

template <typename S>
int
IW_TDT::echo_dataitem (const char * tag,
                       int len_tag,
                       int which_one,
                       S & output) const
{
  int i = _find_index_in_end_array(tag, len_tag, which_one);

  if (i < 0)
  {
    cerr << "IW_TDT::echo_dataitem:no " << which_one << " instance of '";
    cerr.write(tag, len_tag);
    cerr << "' in TDT\n";
    return 0;
  }

  int istart, istop;
  if (0 == i)
  {
    istart = 0;
    istop = _end[0];
  }
  else
  {
    istart = _end[i - 1];
    istop = _end[i];
  }

  output.write(_zdata.rawchars() + istart, istop - istart);

  return output.good();
}

template int IW_TDT::echo_dataitem(const char * tag, int len_tag, int which_one, IWString_and_File_Descriptor & output) const;
template int IW_TDT::echo_dataitem(const char * tag, int len_tag, int which_one, std::ostream & output) const;

int
IW_TDT::count_dataitems (const char * tag,
                         int len_tag) const
{
  int i = 0;
  const_IWSubstring token;

  int rc = 0;

  while (next_dataitem(token, i))
  {
    if (0 == token.strncmp(tag, len_tag))
      rc++;
  }

  return rc;
}

int
IW_TDT::count_dataitems(RE2 & rx) const
{
  int i = 0;
  const_IWSubstring token;

  int rc = 0;
  while (next_dataitem(token, i))
  {
    if (iwre2::RE2PartialMatch(token, rx))
      rc++;
  }

  return rc;
}

bool
IW_TDT::operator == (const IW_TDT & rhs) const
{
  if (_zdata != rhs._zdata)
    return false;

#ifdef DEBUG_TDT_EQUALS
  cerr << "Operator == content matches " << _end.number_elements() << " vs " << rhs._end.number_elements() << endl;

  for (int i = 0; i < _end.number_elements(); i++)
  {
    cerr << _end[i] << " vs " << rhs._end[i] << endl;
  }
#endif

  return _end == rhs._end;
}

int
IW_TDT::Build(const const_IWSubstring& input) {

  _zdata.resize_keep_storage(0);
  _end.resize_keep_storage(0);

  int got_vbar = 0;

  const_IWSubstring line;
  int i = 0;
  while (input.nextword(line, i, '\n'))
  {
//  cerr << "Read ---'" << buffer << "'---\n";
    if (_zdata.length())
      _end.add(_zdata.length());
      
    _zdata += line;

//  cerr << "Data now '" << _zdata << "'\n";

    if (_include_newlines_in_tdt)
      _zdata += '\n';

    if ('|' == line)
    {
      got_vbar = 1;
      break;
    }

    if (iwtdt_careful_mode && ! valid_tdt_form(line))
    {
      cerr << "Invalid tdt form '" << line << "'\n";
      return 0;
    }
  }

  if (_zdata.length() > 0 && ! got_vbar)
  {
    cerr << "IW_TDT::Build:improperly terminated TDT, " << input << "\n";
    return 0;
  }

  return _end.number_elements();
}
