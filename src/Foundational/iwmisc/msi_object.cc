#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#define MSI_OBJECT_READ_IMPLEMENTATION

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/iwmmap.h"

#include "msi_object.h"

static int convert_tags_to_lowercase = 1;

void
set_convert_msi_tags_to_lowercase (int s)
{
  convert_tags_to_lowercase = s;

  return;
}

int
convert_msi_tags_to_lowercase()
{
  return convert_tags_to_lowercase;
}

/*
  Parse a string as a numeric and add result to an array
*/

template <class T>
int
_fetch_numeric_value (const IWString & buffer, resizable_array<T> & vals)
{
  T tmp;

  if (! buffer.numeric_value(tmp))
  {
    cerr << "Fetch_numeric_value: '" << buffer << "' cannot be parsed as numeric " << sizeof(T) << endl;
    return 0;
  }

  vals.add(tmp);

  return 1;
}

template int _fetch_numeric_value(const IWString &, resizable_array<int> &);
template int _fetch_numeric_value(const IWString &, resizable_array<double> &);

//#define DEBUG_MSI_ATTRIBUTE_CREATION

msi_attribute::msi_attribute ()
{
  _default_values();
}

void
msi_attribute::_default_values ()
{
  _valid = 0;
  _type = 0;

  return;
}

int
msi_attribute::build (const IWString & buf)
{
#ifdef DEBUG_MSI_ATTRIBUTE_CREATION
  cerr << "Creating msi attribute from '" << buffer << "'\n";
#endif

  _default_values();

  int nw = buf.nwords();

  if (nw < 4)
  {
    cerr << "msi_attribute::build:must be at least 4 tokens '" << buf << "'\n";
    return 0;
  }

  const_IWSubstring token;
  int i = 0;

  (void) buf.nextword(token, i);

  int needs_closing_paren = 0;

  if ("(A" == token)
  {
    if (! buf.ends_with(')'))
    {
      cerr << "msi_attribute::build:unbalanced parentheses '" << buf << "'\n";
      return 0;
    }

    needs_closing_paren = 1;
  }
  else if ('A' == token)
    ;
  else
  {
    cerr << "msi_attribute::buid:unrecognised '" << buf << "'\n";
    return 0;
  }

  (void) buf.nextword(token, i);

  if ('I' == token)
    _type = MSI_ATTRIBUTE_TYPE_INT;
  else if ('F' == token)
    _type = MSI_ATTRIBUTE_TYPE_FLOAT;
  else if ('D' == token)
    _type = MSI_ATTRIBUTE_TYPE_DOUBLE;
  else if ('C' == token)
    _type = MSI_ATTRIBUTE_TYPE_STRING;
  else if ('O' == token)
    _type = MSI_ATTRIBUTE_TYPE_OBJECT;
  else
  {
    cerr << "msi_attribute: unrecognised type '" << token << "' from " << buf << "\n";
    return 0;
  }
    
  buf.nextword(_name, i);

  if (convert_tags_to_lowercase)
    _name.to_lowercase();

// Grab the string representation - all the rest of the tokens on the line

  _string.resize_keep_storage(0);

  while (buf.nextword(token, i))
  {
    _string.append_with_spacer(token);
  }

  if (needs_closing_paren)
  {
    if (_string.ends_with(')'))
      _string.chop();
    else
    {
      cerr << "msi_attribute::buid:unbalanced paren '" << buf << "'\n";
      return 0;
    }
  }
  else if (_string.ends_with(')'))    // does not need a closing paren, but one present. Too dangerous
  {
    cerr << "msi_attribute::build:possibly mismatched paren '" << buf << "'\n";
    return 0;
  }

  if (MSI_ATTRIBUTE_TYPE_STRING == _type)
  {
    if (_string.starts_with('"') && _string.ends_with('"'))
    {
      _string.remove_leading_chars(1);
      _string.chop();
    }
    else if (_string.starts_with('\'') && _string.ends_with('\''))
    {
      _string.remove_leading_chars(1);
      _string.chop();
    }
    else if (_string.starts_with('"') || _string.ends_with('"'))
    {
      cerr << "msi_attribute::build:unbalanced quotes '" << buf << "'\n";
      return 0;
    }

    return 1;
  }

  if (_string.starts_with('(') && _string.ends_with(')'))
  {
    _string.remove_leading_chars(1);
    _string.chop();
  }
  else if (_string.starts_with('(') || _string.ends_with(')'))
  {
    cerr << "msi_attribute::build:unbalanced parentheses '" << buf << "'\n";
    return 0;
  }

  int n = _string.split(_string_values);

  for (int i = 0; i < n; i++)
  {
    const const_IWSubstring & token = *(_string_values[i]);

    int rc;

    if (MSI_ATTRIBUTE_TYPE_INT == _type)
      rc = _fetch_numeric_value(token, _int_values);
    else if (MSI_ATTRIBUTE_TYPE_DOUBLE == _type || MSI_ATTRIBUTE_TYPE_FLOAT == _type)
      rc = _fetch_numeric_value(token, _double_values);
    else
    {
      cerr << "msi_object::build:what kind is this " << _type << "\n";
      return 0;
    }

    if (0 == rc)
    {
      cerr << "msi_attribute::build:invalid numeric '" << buf << "'\n";
      return 0;
    }

    if (MSI_ATTRIBUTE_TYPE_INT == _type)
      _double_values.add(static_cast<double>(_int_values.last_item()));
  }

  _valid = 1;

  return 1;
}

msi_attribute::~msi_attribute ()
{
  _valid = 0;
  _type = 0;
}

/*
  number_string_values is complicated by the fact that if we only
  get one token on input, we do not bother putting a copy in _string_values
*/

int
msi_attribute::number_string_values () const
{
  if (_string_values.number_elements())
    return _string_values.number_elements();

  return 1;
}

const IWString *
msi_attribute::string_multi_value (int i) const
{
  if (0 == _string_values.number_elements())
  {
    if (0 == i)
      return &_string;
    else
      return nullptr;
  }

  if (_string_values.ok_index(i))
    return _string_values[i];
  else
    return nullptr;
}

int
msi_attribute::value (int & i) const
{
  if (1 != _int_values.number_elements())
  {
    cerr << "msi_attribute::cannot get int 'value' for attribute with " << _int_values.number_elements() << " values\n";
    return 0;
  }

  i = _int_values[0];
  return 1;
}

/*
  Note that treating unsigned int's this was is not really very
  safe....
*/

int 
msi_attribute::value (unsigned int & i) const
{
  if (1 != _int_values.number_elements())
  {
    cerr << "msi_attribute::cannot get unsigned int 'value' for attribute with " << _int_values.number_elements() << " values\n";
    return 0;
  }

  if (_int_values[0] < 0)
    return 0;

  i = (unsigned int) _int_values[0];

  return 1;
}

int
msi_attribute::value (float & x) const
{
  if (1 != _double_values.number_elements())
  {
    cerr << "msi_attribute::cannot get float 'value' for attribute with " << _double_values.number_elements() << " values\n";
    return 0;
  }

  x = float(_double_values[0]);
  return 1;
}


int
msi_attribute::value (double & x) const
{
  if (1 != _double_values.number_elements())
  {
    cerr << "msi_attribute::cannot get double 'value' for attribute with " << _double_values.number_elements() << " values\n";
    return 0;
  }

  x = _double_values[0];
  return 1;
}

int
msi_attribute::value (IWString & x) const
{
  x = _string;

  return x.length();
}

int
msi_attribute::value (const_IWSubstring & x) const
{
  x = _string;

  return x.length();
}

int
msi_attribute::fetch_double_multiple_values (int n, double * result) const
{
  assert (n <= _double_values.number_elements());

  for (int i = 0; i < n; i++)
  {
    *result++ = _double_values[i];
  }

  return n;
}

template <typename T>
int
fetch_next_value (const resizable_array<T> & values,
                  T & v,
                  int & ndx)
{
  if (ndx >= values.number_elements())
    return 0;

  assert (ndx >= 0);

  v = values[ndx];

  ndx++;

  return 1;
}

template int fetch_next_value (const resizable_array<int> &, int &, int &);
template int fetch_next_value (const resizable_array<double> &, double &, int &);

int
msi_attribute::next_value (int & v, int & ndx) const
{
  return fetch_next_value(_int_values, v, ndx);
}

int
msi_attribute::next_value (double & v, int & ndx) const
{
  return fetch_next_value(_double_values, v, ndx);
}

int
msi_attribute::next_value (IWString & v, int & ndx) const
{
  if (ndx >= _string_values.number_elements())
    return 0;

  assert (ndx >= 0);

  v = *(_string_values[ndx]);

  ndx++;

  return 1;
}

std::ostream &
operator << (std::ostream & os, const msi_attribute & msi)
{
  assert (os.good());

  int add_quotes = 0;
  int add_closing_paren = 0;

  os << "  (A ";      // this is an attribute

  if (1 == msi.number_int_values())
    os << "I ";
  else if (msi.number_int_values() > 1)
  {
    os << "I ";
    add_closing_paren = 1;
  }
  else if (1 == msi.number_double_values())
    os << "D ";
  else if (msi.number_double_values() > 1)
  {
    os << "D ";
    add_closing_paren = 1;
  }
  else if (msi._name == "ACL")   // for some reason, MSI think this is an Int
  {
    add_quotes = 1;
    os << "I ";
  }
  else
  {
    os << "C ";
    add_quotes = 1;
  }

  os << msi.name() << " ";

  if (add_closing_paren)
    os << '(';
  else if (add_quotes)
    os << "\"";

  os << msi.stringval();

  if (add_closing_paren)
    os << ")";
  else if (add_quotes)
    os << "\"";

  return os << ")\n";
}

/*
  Several of these do not need to be member functions, but template instantiation goes
  better if they are known
*/

int
msi_object::_looks_like_new_object (const IWString & buffer) const
{
  if (! buffer.starts_with('('))
    return 0;

 if (2 != buffer.nwords())
    return 0;

  return 1;
}

int
msi_object::_looks_like_attribute (const IWString & buffer) const
{
  if (buffer.starts_with("(A "))
    return 1;

  if (buffer.starts_with("A "))
    return 1;

  return 0;
}

//template class resizable_array_p<msi_object>;
//template class resizable_array_base<msi_object *>;
//template class resizable_array_p<msi_attribute>;
//template class resizable_array_base<msi_attribute *>;

#define MSI_OBJECT_ID_INVALID -1

/*
  For efficiency, we resize the attribute array to size 6, as that seems to
  be what Cerius puts out for atoms.
*/

msi_object::msi_object ()
{
  _attributes.resize(6);

  _object_id = MSI_OBJECT_ID_INVALID;
  _valid = 0;

  _display_no_data_error_message = 1;
}

msi_object::~msi_object ()
{
  if (! ok())
  {
    cerr << "Destroying non-ok msi object\n" << (*this);
    abort();
  }

  _object_id = MSI_OBJECT_ID_INVALID;
  _valid = 0;
}

int
msi_object::ok () const
{
  if (0 == _number_elements && 0 == _attributes.number_elements())
    return 1;

  return _valid;
}

int
msi_object::ignore (const char * c)
{
  assert (nullptr != c);

  IWString * tmp = new IWString(c);
  return _attributes_to_ignore.add(tmp);
}

int
msi_object::_matches_a_rejection (const IWString & buffer) const
{
  for (int i = 0; i < _attributes_to_ignore.number_elements(); i++)
  {
    if (buffer.contains(*(_attributes_to_ignore[i])))
      return 1;
  }

  return 0;
}

int
msi_object::_trim_stuff_beyond_last_paren(IWString & buffer)
{
  int paren_level = 0;

  int s = buffer.length();

  for (int i = 0; i < s; i++)
  {
    if ('(' == buffer[i])
      paren_level++;
    else if (')' == buffer[i])
    {
      paren_level--;
      if (0 == paren_level)
      {
        buffer.resize_keep_storage(i + 1);
        return 1;
      }
    }
  }

  return 1;
}
/*
  When I did the lowercase conversion stuff I didn't change this
  function. Beware of potential problems...
*/

int
msi_object::attribute_count (const char * attribute_name) const
{
  int number_attributes = _attributes.number_elements();

  int nfound = 0;

  for (int i = 0; i < number_attributes; i++)
  {
    const msi_attribute * att = _attributes[i];
    if (att->name() == attribute_name)
      nfound++;
  }

  return nfound;
}

int
msi_object::object_count (int j) const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
    if (j == _things[i]->_object_id)
      rc++;

  return rc;
}

/*
  No allowance made for lowercase conversion
*/

int
msi_object::object_count (const char * s) const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (s == _things[i]->_name)
      rc++;
  }

  return rc;
}

/*
  No allowance for lowercase
*/

static const msi_attribute *
_locate_attribute (const const_IWSubstring & attribute_name,
                   int which_one_to_find,
                   const resizable_array_p<msi_attribute> & attributes)
{
  int number_attributes = attributes.number_elements();
  int nfound = 0;
  for (int i = 0; i < number_attributes; i++)
  {
    const msi_attribute * att = attributes[i];
    if (att->name() == attribute_name)
    {
      if (nfound == which_one_to_find)
        return att;
      nfound++;
    }
  }

  return nullptr;
}

const msi_attribute *
msi_object::attribute (const char * attribute_name, int which_one_to_find) const
{
  return _locate_attribute(attribute_name, which_one_to_find, _attributes);
}

const msi_attribute *
msi_object::attribute (const IWString & attribute_name, int which_one_to_find) const
{
  return _locate_attribute(attribute_name, which_one_to_find, _attributes);
}

const msi_attribute *
msi_object::attribute (int i) const
{
  assert (i >= 0);

  int number_attributes = _attributes.number_elements();
  if (i >= number_attributes)
    return nullptr;

  return _attributes[i];
}

std::ostream &
operator << (std::ostream & os, const msi_object & msi)
{
  assert (os.good());
  os << "(" << msi.object_id() << " " << msi.name() << endl;

  int n = msi._attributes.number_elements();
  for (int i = 0; i < n; i++)
  {
    const msi_attribute * att = msi._attributes[i];
    os << " ";
    os << *att;
  }

  n = msi.number_elements();
  for (int i = 0; i < n; i++)
  {
    const msi_object * child = msi[i];
    child->print(os, 2);
//  os << " " << *child;
  }

  return os << ")\n";
}

int
msi_object::print (std::ostream & os, int indentation) const
{
  assert (ok());
  assert (os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  os << ind << "(" << _object_id << " " << _name << endl;
  int na = _attributes.number_elements();
  for (int i = 0; i < na; i++)
  {
    const msi_attribute * a = _attributes[i];
    os << ind << (*a);
  }

  for (int i = 0; i < _number_elements; i++)
  {
    const msi_object * m = _things[i];
    m->print(os, indentation + 2);
  }
  os << ind << ")\n";

  return os.good();
}

int
msi_object::string_value_for_attribute (const char * attribute_name, 
               IWString & zstring) const
{
  assert (zstring.ok());
  int na = _attributes.number_elements();

  for (int i = 0; i < na; i++)
  {
    if (_attributes[i]->name() == attribute_name)
    {
      zstring = _attributes[i]->stringval();
      return 1;
    }
  }

  return 0;
}

int
msi_object::object_count () const
{
  assert (ok());

  int rc = _number_elements;

  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->object_count();
  }

  return rc;
}

#include "iwminmax.h"

int
msi_object::highest_object_id () const
{
  assert (ok());

  iwmax<int> rc (_object_id);

  for (int i = 0; i < _number_elements; i++)
    rc.extra(_things[i]->highest_object_id());

  return rc.maxval();
}

const msi_object *
msi_object::component (const const_IWSubstring & id, int which_one_to_find) const
{
  int nfound = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    const msi_object * m = _things[i];
    if (id == m->name())
    {
      if (nfound == which_one_to_find)
        return m;

      nfound++;
    }
  }

  return nullptr;
}

void
msi_object::names_to_lowercase ()
{
  assert (ok());

  _name.to_lowercase();

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->names_to_lowercase();
  }

  return;
}

template int msi_object::read(iwstring_data_source&);
