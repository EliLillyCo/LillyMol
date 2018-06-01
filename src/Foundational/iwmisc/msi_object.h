#ifndef IW_MSI_OBJECT_H
#define IW_MSI_OBJECT_H 1

#include <iostream>

using std::cerr;
using std::endl;

#include "iwstring.h"
#include "iwstring_data_source.h"

#define MSI_ATTRIBUTE_TYPE_INT     3
#define MSI_ATTRIBUTE_TYPE_FLOAT   4
#define MSI_ATTRIBUTE_TYPE_DOUBLE  5
#define MSI_ATTRIBUTE_TYPE_STRING  6
#define MSI_ATTRIBUTE_TYPE_OBJECT  8

class msi_attribute
{
  friend
    std::ostream &
      operator << (std::ostream &, const msi_attribute &);
  private:
    IWString _name;
    int    _type;    // bad design, yes....

    IWString _string;

    resizable_array_p<IWString> _string_values;
    resizable_array<int>        _int_values;
    resizable_array<double>     _double_values;

    int    _valid;

    void _default_values ();

  public:
    msi_attribute (); //char *);
    ~msi_attribute ();

    int ok () const;

    int build (const IWString &);

    int valid_as_int    () const { return _int_values.number_elements ();}
    int valid_as_double () const { return _double_values.number_elements ();}

    int number_double_values   () const { return _double_values.number_elements ();}
    int number_int_values      () const { return _int_values.number_elements ();}
    int number_string_values   () const;

    const IWString & name () const { return _name;}
    const IWString & stringval () const { return _string;}

    int value (int &) const;
    int value (unsigned int &) const;
    int value (float &) const;
    int value (double &) const;
    int value (IWString &) const;
    int value (const_IWSubstring &) const;

    int  fetch_int_multiple_values (int, int *) const;
    int  int_multi_value (int i) const { return _int_values[i];}

    int  fetch_double_multiple_values (int, double *) const;
    int  fetch_double_multiple_values (resizable_array<double> & result) const
           { return result.copy (_double_values);}
    double double_multi_value (int i) const { return _double_values[i];}

    const IWString * string_multi_value (int i) const;

//  First arg is the value to be returned, the 2nd arg is just an index

    int next_value (int &, int &) const;
    int next_value (double &, int &) const;
    int next_value (IWString &, int &) const;
};

class msi_object : public resizable_array_p<msi_object>
{
  friend
    std::ostream & operator << (std::ostream &, const msi_object &);

  private:
    resizable_array_p<msi_attribute> _attributes;
    int      _object_id;
    IWString _name;
    int      _valid;

//  For efficiency, we can ignore various strings when reading attributes
//  ideally, these would be regular expressions, but I don't trust the rx
//  class yet.

    resizable_array_p<IWString> _attributes_to_ignore;

    int _display_no_data_error_message;

//  Private functions

    int _matches_a_rejection (const IWString &) const;
    int _looks_like_new_object (const IWString & buffer) const;
    int _looks_like_attribute (const IWString & buffer) const;
    int _trim_stuff_beyond_last_paren(IWString & buffer);
    
  public:
    msi_object ();
    ~msi_object ();

    int ok () const;
    int print (std::ostream &, int = 0) const;

    void set_display_no_data_error_message(const int s) { _display_no_data_error_message = s;}

    int object_id () const { return _object_id;}

    int active () const { return object_count () + attribute_count ();}

    int object_count () const;
    int highest_object_id () const;

    const IWString & name () const { return _name;}

    int ignore (const char *);

    template <typename T> int read (T &);

    void names_to_lowercase ();

    const msi_attribute * attribute (const char *, int = 0) const;
    const msi_attribute * attribute (const IWString &, int = 0) const;

    const msi_attribute * attribute (int) const;

    const msi_object * component (const const_IWSubstring &, int = 0) const;

    int  add_attribute (msi_attribute * extra) {return _attributes.add (extra);}

    int  attribute_count () const { return _attributes.number_elements ();};
    int  number_attributes () const { return _attributes.number_elements ();}
    int  attribute_count (const char *) const;

    int  object_count    (int) const;
    int  object_count    (const char *) const;

    int  string_value_for_attribute (const char *, IWString &) const;
};

extern void set_convert_msi_tags_to_lowercase (int);
extern int convert_msi_tags_to_lowercase();

#if defined(MSI_OBJECT_READ_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

template <typename T>
int
msi_object::read (T & input)
{
  assert (input.ok());

  input.set_strip_leading_blanks(1);
  input.set_skip_blank_lines(1);

  resize_keep_storage(0);
  _attributes.resize_keep_storage(0);
  _object_id = -7;
  _name.resize_keep_storage(0);

  _valid = 0;

  IWString buffer;

  if (! input.next_record(buffer))
  {
    if (_display_no_data_error_message)
      cerr << "msi_object::read:no data\n";
    return 0;
  }

  assert (buffer.length());

  if ('(' != buffer[0])
    return 0;

  _name = buffer.word(1);
  if (convert_msi_tags_to_lowercase())
    _name.to_lowercase();

  if (! is_int(buffer.chars() + 1, &_object_id))
  {
    cerr << "msi_object: cannot discern object id '" << buffer << "'\n";
    return 0;
  }

  assert (_object_id >= 0);

  while (input.next_record(buffer))
  {
    assert (buffer.length() > 0);

    if (')' == buffer[0])    // end of this object
    {
      _valid = 1;
      return 1;
    }

    buffer.strip_trailing_blanks();

    if (_looks_like_new_object(buffer))
    {
      input.push_record();
      msi_object * msi = new msi_object();
      if (! msi->read(input))
      {
        delete msi;
        return 0;
      }

      add(msi);
    }
    else if (_matches_a_rejection(buffer))
      ;             // ignore it
    else if (_looks_like_attribute(buffer))
    {
      _trim_stuff_beyond_last_paren(buffer);

      msi_attribute * att = new msi_attribute();
      if (! att->build(buffer))
      {
        cerr << "msi_object::read: cannot build msi attribute\n";
        cerr << buffer << endl;
        return 0;
      }
      add_attribute(att);
    }
    else if ('#' == buffer[0])    // comment line
      ;
    else
    {
      cerr << "msi_object::read:unrecognised input '" << buffer << "'\n";
      cerr << "IGNORED\n";
    }
  }

  cerr << "msi_object::read  unterminated object\n";
  return 0;
}

#endif


#endif
