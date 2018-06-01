#ifndef IW_TDT_IMPLEMENTATION_H
#define IW_TDT_IMPLEMENTATION_H

#include "iw_tdt.h"

/*
  Fetch via string
*/

//#define DEBUG_DATAITEM_VALUE

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
  int i = 0;
  const_IWSubstring token;

  int nfound = 0;

  while (next_dataitem (token, i))
  {
    if (token.length () < len_tag)
      continue;

    if (0 != token.strncmp (tag, len_tag))
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
