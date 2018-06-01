#include "iwtokeniser.h"

static const char dquote('"');

IWTokeniser::IWTokeniser(const const_IWSubstring & s) : _buffer(s)
{
  _i = 0;
  _empty_fields_valid = 0;
  _quoted_tokens = 0;
  _in_quote = 0;

  _sep = ' ';
}

template <typename S>
int
IWTokeniser::next_token(S & token)
{
  if (_i >= _buffer.length())
    return 0;

  if (_quoted_tokens)
    return _next_token_with_quotes(token);

  if (_empty_fields_valid)
    return _buffer.nextword_single_delimiter(token, _i, _sep);

  return _buffer.nextword(token, _i, _sep);
}

template <typename S>
int
IWTokeniser::_next_token_with_quotes(S & token)
{
  const int iprev = _i;

  for (; _i < _buffer.length(); ++_i)
  {
    if (dquote == _buffer[_i])
      _in_quote = ! _in_quote;
    else if (_in_quote)
      ;
    else if (_sep == _buffer[_i])
    {
      if (dquote == _buffer[_i-1] && dquote == _buffer[iprev])
        _buffer.from_to(iprev + 1, _i-2, token);
      else
        _buffer.from_to(iprev, _i-1, token);
      _i++;

      return 1;
    }
  }

  if (0 == _buffer.length())
    return 0;

//cerr << "_next_token_with_quotes i = " << _i << endl;

  if (dquote == _buffer[iprev] && dquote == _buffer[_i-1])
    _buffer.from_to(iprev + 1, _i-2, token);
  else
    _buffer.from_to(iprev, _i-1, token);

  return 1;
}

template int IWTokeniser::_next_token_with_quotes(IWString&);
template int IWTokeniser::_next_token_with_quotes(const_IWSubstring&);
template int IWTokeniser::next_token(IWString&);
template int IWTokeniser::next_token(const_IWSubstring&);
