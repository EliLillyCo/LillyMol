#ifndef IW_TABULAR_DATA_H
#define IW_TABULAR_DATA_H

#include "iwstring_data_source.h"
#include "iwbits.h"

template <typename T>
class IW_Tabular_Data
{
  private:
    int _nrows;
    int _ncols;
    T * _zdata;

    int _has_header;
    IWString * _header;

    IWString _missing_value;

    IW_Bits_Base _missing;

    int _first_column_is_identifier;

    IWString * _id;

//  private functions

    int _determine_number_columns(const const_IWSubstring & buffer, const char sep);
    int _parse_row(const const_IWSubstring & buffer, const char sep, const int r);
    int _determine_rows_and_columns(iwstring_data_source & input, const char sep, const_IWSubstring & buffer);
    int _parse_header(const const_IWSubstring & buffer, const char sep);

  public:
    IW_Tabular_Data();
    ~IW_Tabular_Data();

    int nrows() const { return _nrows;}
    int ncols() const { return _ncols;}

    void set_missing_value_string(const const_IWSubstring & s) { _missing_value = s;}

    void set_has_header(const int s) { _has_header = s;}
    void set_first_column_is_identifier(const int s) { _first_column_is_identifier = s;}

    int resize (const int nr, const int nc);

    void set(const int r, const int c, T v) { _zdata[r * _ncols + c] = v;}

    const T * data () const { return _zdata;}
    T * data () { return _zdata;}

    const IWString * ids () const { return _id;}

    int build(const char * fname, const char sep);
    int build(iwstring_data_source &, const char sep);

    int remove_column(const int c);

    template <typename OP> void const_each_row(OP &) const;
    template <typename OP> void change_each_row(OP &);

    template <typename O> int do_write(O &, const char sep) const;
};

#if defined(IW_TABULAR_DATA_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

template <typename T> template <typename OP>
void
IW_Tabular_Data<T>::const_each_row(OP & op) const
{
  for (int i = 0; i < _nrows; ++i)
  {
    op(_zdata + i * _ncols);
  }

  return;
}

template <typename T> template <typename OP>
void
IW_Tabular_Data<T>::change_each_row(OP & op)
{
  for (int i = 0; i < _nrows; ++i)
  {
    op(_zdata + i * _ncols);
  }

  return;
}

#endif

#endif
