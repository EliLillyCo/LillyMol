#ifndef IWCONST_TDT_H
#define IWCONST_TDT_H

#include "re2/re2.h"

#include "Foundational/iwstring/iwstring.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/set_or_unset.h"

/*
  Variant on the TDT class, hopefully with better efficiency

  We retain all info in the TDT, including newlines, and the vertical bar
*/

class IW_TDT
{
  private:
    IWString _zdata;

//  Where in _zdata does each dataitem end

    resizable_array<int> _end;

    Set_or_Unset<off_t> _offset;

//  private functions

    int _find_index_in_end_array (const char *, int n, int) const;

    template <typename S> int _dataitem_value_string (const char * tag, int len_tag, S & s, int which_one) const;
    template <typename S> int _dataitem_value (const char * tag, int len_tag, S & s, int which_one) const;
    template <typename S> int _common_write_all_except_vbar (S &) const;
    template <typename S> int _add_dataitem (const char * tag, int len_tag, const S & s, int where_to_put = -1);

  public:
    IW_TDT ();

    int ok () const;

    int number_elements () const { return _end.number_elements ();}

    int index_of_dataitem (const const_IWSubstring &) const;

    int remove_item (int);
    int remove_items_with_tag (const IWString &);
    int remove_items_with_tag (const const_IWSubstring &);

    int remove_items_with_tag (const char *, int);    // does not check that tag ends with '<', so potentially dangerous

    int insert_before (int, const const_IWSubstring &);

    int item (int, const_IWSubstring &) const;

    int next (iwstring_data_source &);

    int Build(const const_IWSubstring & s);    // build from string with embedded newlines
    int Build(const char *, int);    // build from string with embedded newlines

    int next_dataitem (const_IWSubstring &, int &) const;

    int next_dataitem_value (const_IWSubstring & ztag, const_IWSubstring & zdata, int &) const;

    int do_write (std::ostream &) const;

    template <typename S>
    int echo_dataitem (const char * tag, int len_tag, int which_one, S & output) const;

    int echo_dataitem (const IWString & tag, int which_one, std::ostream & output) const { return echo_dataitem (tag.rawchars (), tag.length (), which_one, output);}
    int echo_dataitem (const const_IWSubstring & tag, int which_one, std::ostream & output) const { return echo_dataitem (tag.rawchars (), tag.length (), which_one, output);}

    int echo_dataitem (const const_IWSubstring & tag, int which_one, IWString_and_File_Descriptor & output) const { return echo_dataitem (tag.rawchars (), tag.length (), which_one, output);}

    int write_all_except_vbar (std::ostream &) const;
    int write_all_except_vbar (IWString_and_File_Descriptor &) const;

    int dataitem_value (const char * tag, int len_tag, IWString & s, int which_one = 0) const { return _dataitem_value_string (tag, len_tag,                         s, which_one);}
    int dataitem_value (const char * tag,              IWString & s, int which_one = 0) const { return _dataitem_value_string (tag, static_cast<int>(::strlen(tag)), s, which_one);}
    int dataitem_value (const const_IWSubstring & tag, IWString & s, int which_one = 0) const { return _dataitem_value_string (tag.rawchars (), tag.length (),       s, which_one);}
    int dataitem_value (const IWString & tag,          IWString & s, int which_one = 0) const { return _dataitem_value_string (tag.rawchars (), tag.length (),       s, which_one);}

    int dataitem_value (const char * tag, int len_tag, const_IWSubstring & s, int which_one = 0) const { return _dataitem_value_string (tag,             len_tag,                         s, which_one);}
    int dataitem_value (const char * tag,              const_IWSubstring & s, int which_one = 0) const { return _dataitem_value_string (tag,             static_cast<int>(::strlen(tag)), s, which_one);}
    int dataitem_value (const const_IWSubstring & tag, const_IWSubstring & s, int which_one = 0) const { return _dataitem_value_string (tag.rawchars (), tag.length (),                   s, which_one);}
    int dataitem_value (const IWString & tag,          const_IWSubstring & s, int which_one = 0) const { return _dataitem_value_string (tag.rawchars (), tag.length (),                   s, which_one);}

    int dataitem_value (const char * tag, int len_tag, int & s, int which_one = 0) const { return _dataitem_value (tag,             len_tag,                         s, which_one);}
    int dataitem_value (const char * tag,              int & s, int which_one = 0) const { return _dataitem_value (tag,             static_cast<int>(::strlen(tag)), s, which_one);}
    int dataitem_value (const const_IWSubstring & tag, int & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}
    int dataitem_value (const IWString & tag,          int & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}

    int dataitem_value (const char * tag, int len_tag, float & s, int which_one = 0) const { return _dataitem_value (tag,             len_tag,                         s, which_one);}
    int dataitem_value (const char * tag,              float & s, int which_one = 0) const { return _dataitem_value (tag,             static_cast<int>(::strlen(tag)), s, which_one);}
    int dataitem_value (const const_IWSubstring & tag, float & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}
    int dataitem_value (const IWString & tag,          float & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}

    int dataitem_value (const char * tag, int len_tag, double & s, int which_one = 0) const { return _dataitem_value (tag,             len_tag,                         s, which_one);}
    int dataitem_value (const char * tag,              double & s, int which_one = 0) const { return _dataitem_value (tag,             static_cast<int>(::strlen(tag)), s, which_one);}
    int dataitem_value (const const_IWSubstring & tag, double & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}
    int dataitem_value (const IWString & tag,          double & s, int which_one = 0) const { return _dataitem_value (tag.rawchars (), tag.length (),                   s, which_one);}

    int remove_all (RE2 &);

    int dataitem_value (RE2 & rx,
                        const_IWSubstring & dataitem,
                        const_IWSubstring & zresult,
                        int which_to_return) const;

    template <typename S> int dataitem (const char * tag, int len_tag, S & s, int w = 0) const;
    template <typename S> int dataitem (const char * tag,              S & s, int w = 0) const { return dataitem (tag, static_cast<int>(::strlen(tag)), s, w);}
    template <typename S> int dataitem (const const_IWSubstring & tag, S & s, int w = 0) const { return dataitem (tag.rawchars (), tag.length (),       s, w);}
    template <typename S> int dataitem (const IWString & tag,          S & s, int w = 0) const { return dataitem (tag.rawchars (), tag.length (),       s, w);}

    const IWString & rawdata () const { return _zdata;}

    template <typename S> int set_dataitem_value (const char * tag, int len_tag, const S & s, int which_one = 0);
    template <typename S> int set_dataitem_value (const char * tag, const S & s, int which_one = 0)              { return set_dataitem_value (tag, static_cast<int>(::strlen(tag)), s, which_one);}
    template <typename S> int set_dataitem_value (const IWString & tag, const S & s, int which_one = 0)          { return set_dataitem_value (tag.rawchars (), tag.length (), s, which_one);}
    template <typename S> int set_dataitem_value (const const_IWSubstring & tag, const S & s, int which_one = 0) { return set_dataitem_value (tag.rawchars (), tag.length (), s, which_one);}

    template <typename S> int add_dataitem (const char * tag, int len_tag, const S & s, int where_to_put = -1) { return _add_dataitem(tag, len_tag, s, where_to_put);}

    template <typename S> int add_dataitem (const char * tag, const S & s, int where_to_put = -1) { return _add_dataitem (tag, static_cast<int>(::strlen(tag)), s, where_to_put);}
    template <typename S> int add_dataitem (const IWString & tag, const S & s, int where_to_put = -1) { return _add_dataitem (tag.rawchars (), tag.length (), s, where_to_put);}
    template <typename S> int add_dataitem (const const_IWSubstring & tag, const S & s, int where_to_put = -1) { return _add_dataitem (tag.rawchars (), tag.length (), s, where_to_put);}

    int count_dataitems (const char * tag, int len_tag) const;
    int count_dataitems (const char * tag) const { return count_dataitems (tag, static_cast<int>(strlen (tag)));}
    int count_dataitems (const const_IWSubstring & tag) const { return count_dataitems (tag.rawchars (), tag.length ());}
    int count_dataitems (const IWString & tag) const { return count_dataitems (tag.rawchars (), tag.length ());}

    int count_dataitems (RE2 &) const;

    bool operator == (const IW_TDT &) const;
};

extern std::ostream & operator << (std::ostream &, const IW_TDT &);
extern IWString & operator << (IWString &, const IW_TDT &);

extern void set_iwtdt_careful_mode (int);

extern void set_include_newlines_in_tdt (int);

extern int include_newlines_in_tdt ();


#ifdef DATAITEM_VALUE_TEMPLATE_IMPLEMENTATION
/*
  Fetch via string
*/

//#define DEBUG_DATAITEM_VALUE
#ifdef DEBUG_DATAITEM_VALUE
using std::cerr;
using std::endl;
#endif

template <typename S>
int
IW_TDT::_dataitem_value_string (const char * dkey,
                               int klen,
                               S & zvalue,
                               int which_one) const
{
  int nend = _find_index_in_end_array (dkey, klen, which_one);

#ifdef DEBUG_DATAITEM_VALUE
  cerr << "Looking for '";
  for (int i = 0; i < klen; i++)
  {
    cerr << dkey[i];
  }
  cerr << "' in chunk " << nend << endl;
#endif

  if (nend < 0)
    return 0;

  const char * s = _zdata.rawchars ();

  if (0 == nend)
    zvalue.set (s + klen, _end[0] - 1 - klen);
  else
    zvalue.set (s + _end[nend - 1] + klen, _end[nend] - _end[nend - 1] - 1 - klen);

// If they come in with a tag without a < character, 

  if (zvalue.starts_with ('<') && '<' != dkey[klen - 1])
    zvalue.remove_leading_chars (1);

  if (zvalue.ends_with ('\n'))
    zvalue.chop ();

  if (zvalue.ends_with ('>'))
    zvalue.chop ();

#ifdef DEBUG_DATAITEM_VALUE
  cerr << "_dataitem_value_string: tag '";
  cerr.write (dkey, klen);
  cerr << "' result '" << zvalue << "'\n";
#endif

  return 1;
}

template <typename S>
int
IW_TDT::dataitem (const char * tag,
                  int len_tag,
                  S & s,
                  int which_to_find) const
{
  const const_IWSubstring mytag(tag, len_tag);

  int i = 0;
  const_IWSubstring token;

  int nfound = 0;

  const_IWSubstring last_one_found;

  while (next_dataitem (token, i))
  {
    if (! token.starts_with(mytag))
      continue;

    if ('<' == tag[len_tag - 1])     // TAG included the opening < character
      ;
    else if (token.length () == len_tag)    // hard to imagine this happening
      continue;
    else if ('<' != token[len_tag])
      continue;

    if (nfound == which_to_find)
    {
      s = token;
      return 1;
    }

    nfound++;
    last_one_found = s;
  }

  if (-1 == which_to_find && last_one_found.length())
  {
    s = last_one_found;
    return 1;
  }

  return 0;
}

template <typename T>
int
IW_TDT::_dataitem_value (const char * tag,
                         int len_tag,
                         T & v,
                         int which_one) const
{
  int nend = _find_index_in_end_array (tag, len_tag, which_one);

  if (nend < 0)
    return 0;

  const char * s = _zdata.rawchars ();

  const_IWSubstring tmp;

  if (0 == nend)
    tmp.set (s + len_tag, _end[0] - 2 - len_tag);
  else
    tmp.set (s + _end[nend - 1] + len_tag, _end[nend] - _end[nend - 1] - 2 - len_tag);

  return tmp.numeric_value (v);
}

#endif

#define IWTDT_ADD_DATAITEM_IMPLEMENTATION
#ifdef IWTDT_ADD_DATAITEM_IMPLEMENTATION

template <typename S>
int
IW_TDT::_add_dataitem (const char * tag,
                       int len_tag,
                       const S & s,
                       int where_to_put)
{
//cerr << "Inserting " << s << "\n";

  IWString to_insert;
  to_insert.strncat (tag, len_tag);
  to_insert << s;
  to_insert << '>';
  if (include_newlines_in_tdt ())
    to_insert.add ('\n');

  if (where_to_put < 0)
    where_to_put = _end.number_elements () + where_to_put;
   
  assert (where_to_put >= 0);

//cerr << "Will put '" << to_insert << "' at position " << where_to_put << "size " << _zdata.length() << endl;

  if (where_to_put > 0)
    where_to_put = _end[where_to_put];

  return _zdata.insert (to_insert, where_to_put);
}

#endif

#ifdef SET_DATAITEM_VALUE_IMPLEMENTATION

template <typename S>
int
IW_TDT::set_dataitem_value (const char * tag,
                            int len_tag,
                            const S & new_data,
                            int which_one)
{
  int i = _find_index_in_end_array (tag, len_tag, which_one);
  if (i < 0)
  {
    std::cerr << "IW_TDT::set_dataitem_value:no " << which_one << " instance of ";
    std::cerr.write (tag, len_tag) << "' found\n";
    return 0;
  }

  int e = _end[i];   // end of the current data

  int s;             // start of current data
  if (0 == i)
    s = 0;
  else
    s = _end[i - 1];

  return _zdata.change (s, e, new_data);
}

#endif

#endif
