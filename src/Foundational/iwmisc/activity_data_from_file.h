#ifndef FOUNDATIONAL_IWMISC_ACTIVITY_DATA_FROM_FILE_H_
#define FOUNDATIONAL_IWMISC_ACTIVITY_DATA_FROM_FILE_H_

/*
  I have several programmes that read activity data from a file
*/

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

template <typename T>
class Activity_Data_From_File : public IW_STL_Hash_Map<IWString, T> {
 private:
  int _activity_column;
  int _strip_leading_zeros;
  IWString _activity_column_header;

  int _input_separator_from_file_name;
  char _input_separator;

  // private functions

  int _read_activity_data_record(const const_IWSubstring& buffer);
  int _read_activity_data(iwstring_data_source& input);

 public:
  Activity_Data_From_File();

  void
  set_activity_column(int c) {
    _activity_column = c;
  }

  void
  set_strip_leading_zeros(int s) {
    _strip_leading_zeros = s;
  }

  const IWString&
  activity_column_header() const {
    return _activity_column_header;
  }

  int read_activity_data(const const_IWSubstring& fname);

  int get_activity(const const_IWSubstring&, T&) const;

  template <typename C>
  int construct_from_command_line(const C& cl, char, int verbose);
};

#ifdef ACTIVITY_DATA_IMPLEMENATION_H

template <typename T>
Activity_Data_From_File<T>::Activity_Data_From_File() {
  _activity_column = 1;

  _strip_leading_zeros = 0;

  _input_separator_from_file_name = 0;
  _input_separator = ' ';

  return;
}

template <typename T>
template <typename C>
int
Activity_Data_From_File<T>::construct_from_command_line(const C& cl, char flag,
                                                        int verbose) {
  if (!cl.option_present(flag)) {
    std::cerr << "Activity_Data_From_File::construct_from_command_line: no option '" << flag
         << "' specified\n";
    return 0;
  }

  IWString fname;

  int i = 0;
  IWString token;
  while (cl.value(flag, token, i++)) {
    if (token.starts_with("col=")) {
      token.remove_leading_chars(4);
      if (!token.numeric_value(_activity_column) || _activity_column < 1) {
        std::cerr << "Activity_Data_From_File::construct_from_command_line:invalid column "
                "directive '"
             << token << "'\n";
        return 0;
      }
      _activity_column--;
      continue;
    }

    if (token.starts_with("sep=")) {
      token.remove_leading_chars(4);
      if (token == "auto") {
        _input_separator_from_file_name = 1;
      } else if (!char_name_to_char(token)) {
        std::cerr << "Activity_Data_From_File::construct_from_command_line:unrecognised "
                "separator '"
             << token << "'\n";
        return 0;
      } else {
        _input_separator = token[0];
      }
      continue;
    }

    fname = token;
    if (!read_activity_data(fname)) {
      std::cerr << "Activity_Data_From_File::construct_from_command_line:could not read "
              "activity data from '"
           << fname << "'\n";
      return 0;
    }
  }

  int rc = static_cast<int>((*this).size());
  if (verbose) {
    std::cerr << "Read " << rc << " activity values from " << fname
         << '\n';  // Wrong if multiple files read.
  }

  return rc;
}

template <typename T>
int
Activity_Data_From_File<T>::read_activity_data(const const_IWSubstring& fname) {
  if (_input_separator_from_file_name) {
    _input_separator_from_file_name = iwstring::SeparatorFromFileName(fname);
  }

  iwstring_data_source input(fname);

  if (!input.good()) {
    std::cerr << "Activity_Data_From_File::_read_activity_data:cannot open '" << fname
         << "'\n";
    return 0;
  }

  return _read_activity_data(input);
}

template <typename T>
int
Activity_Data_From_File<T>::_read_activity_data(iwstring_data_source& input) {
  input.set_dos(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (!_read_activity_data_record(buffer)) {
      std::cerr << "Activity_Data_From_File::_read_activity_data_record:cannot read activity "
              "data, line "
           << input.lines_read() << '\n';
      std::cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return static_cast<int>((*this).size());
}

template <typename T>
int
Activity_Data_From_File<T>::_read_activity_data_record(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id;

  if (!buffer.NextWord(id, i, _input_separator)) {
    std::cerr << "Cannot extract identifier\n";
    return 0;
  }

  IWString token;

  if (!buffer.nextword(token, i)) {
    std::cerr << "Not enough tokens on experimental data record\n";
    return 0;
  }

  if (1 != _activity_column) {
    if (!buffer.word(_activity_column, token)) {
      std::cerr
          << "Activity_Data_From_File::_read_activity_data_record:cannot extract column '"
          << (_activity_column + 1) << " from record\n";
      return 0;
    }
  }

  float a;
  if (token.numeric_value(a))
    ;
  else if (0 == (*this).size()) {
    _activity_column_header = token;
  } else {
    std::cerr << "Activity_Data_From_File::_read_activity_data_record:non numeric activity "
            "value '"
         << token << "', id '" << id << "'\n";
    return 0;
  }

  if (_strip_leading_zeros) {
    id.remove_leading_chars('0');
  }

  (*this)[id] = a;

  // std::cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token <<
  // "'\n";

  return 1;
}

template <typename T>
int
Activity_Data_From_File<T>::get_activity(const const_IWSubstring& id, T& a) const {
  typename Activity_Data_From_File<T>::const_iterator f = (*this).find(id);

  if (f != (*this).end()) {
    a = (*f).second;
    return 1;
  }

  if (!_strip_leading_zeros) {
    return 0;
  }

  IWString tmp(id);
  tmp.remove_leading_chars('0');

  f = (*this).find(tmp);

  if (f == (*this).end()) {
    return 0;
  }

  a = (*f).second;

  return 1;
}

#endif /* implementation */

#endif  // FOUNDATIONAL_IWMISC_ACTIVITY_DATA_FROM_FILE_H_
