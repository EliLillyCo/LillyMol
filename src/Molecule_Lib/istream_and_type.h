#ifndef IW_ISTREAM_AND_TYPE
#define IW_ISTREAM_AND_TYPE

#include "Foundational/data_source/iwstring_data_source.h"

#include "iwmtypes.h"

template <typename T>
class data_source_and_type : public iwstring_data_source
{
  private:
    FileType _input_type;
    IWString _fname;
    int _valid;
    int _molecules_read;
    int _connection_table_errors_encountered;
    int _connection_table_errors_allowed;
    int _dos;
    int _verbose;

//  Jun 2000. Often seem to want the ability to skip over molecules and/or
//  only process a certain number of them.

    int _skip_first;
    int _do_only;

    IWString_and_File_Descriptor _stream_for_connection_table_errors;

//  If we are logging connection table errors, we need to know where the
//  most recently read molecule started

    off_t _offset_for_most_recent_molecule;
  
//  private functions

    int _default_values (FileType);
    int _file_opened_ok (const IWString &);

    void _do_skip_first();
    int  _do_merck_skip_records();
    int  _skip_csv_header_record();
    int  _skip_merck_connection_table();

    RE2  _rx_for_input_type () const;

    size_t _average_size (off_t offset, RE2 & rx, int n);

  public:
    data_source_and_type (const IWString &);
    data_source_and_type (FileType, const char *);
    data_source_and_type (FileType, const IWString &);
    data_source_and_type ();
    ~data_source_and_type();

    int ok() const;
//  int debug_print (ostream &) const;

    void set_skip_first(int s) {
      _skip_first = s;
    }
    void set_do_only(int s) {
      _do_only = s;
    }

    int do_open (const char * fname, FileType input_type);

    void set_verbose (int verbose) {_verbose = verbose;}

    int molecules_read() const { return _molecules_read;}

    void set_connection_table_errors_allowed (int i);
    int  set_connection_table_error_file (const IWString &);

//  int next_molecule (Molecule &);
    T * next_molecule();

    int molecules_remaining();

    int connection_table_errors_encountered() const { return _connection_table_errors_encountered;}

    int stopped_because_of_error () const { return _connection_table_errors_encountered > _connection_table_errors_allowed;}

    int estimate_molecules_in_file ();
};

#if defined(ISTREAM_AND_TYPE_IMPLEMENTATION) || defined (IW_IMPLEMENTATIONS_EXPOSED)

#include <iostream>
#include "molecule.h"

//#define INPUT_TYPE_FROM_NO_ARG_CONSTRUCTOR -6273

template <typename T>
int
data_source_and_type<T>::_default_values(FileType input_type)
{
  _valid = 0;

  _molecules_read = 0;
  _verbose = 0;

  _dos = 0;

  _connection_table_errors_encountered = 0;
  _connection_table_errors_allowed = ::number_connection_table_errors_to_skip();

  _skip_first = skip_first_molecules();
  _do_only = do_only_n_molecules();

  if (FILE_TYPE_NO_ARG_CONSTRUCTOR == input_type)
    ;
  else if (! valid_file_type(input_type))
  {
    std::cerr << "data_source_and_type::_default_values: unknown file type " << input_type << '\n';
    return 0;
  }

  _input_type = input_type;
      
  if (FILE_TYPE_SMI == _input_type)
    set_skip_blank_lines(1);

  set_strip_trailing_blanks(1);

  _offset_for_most_recent_molecule = static_cast<off_t>(0);

  return 1;
}

/*
  Called by the constructors to determine whether or not the file
  was opened successfully.
*/

template <typename T>
int
data_source_and_type<T>::_file_opened_ok (const IWString & fname)
{
  if (! good())
  {
    std::cerr << "data_source_and_type::_file_opened_ok: cannot open '" << fname << "' for input\n";
    return 0;
  }

  _fname = fname;

  _valid = 1;

  iwstring_data_source::set_record_delimiter (input_file_delimiter());

  if (input_is_dos_mode())
    iwstring_data_source::set_dos (1);

  if (seek_to_from_command_line() > 0)
  {
    if (! iwstring_data_source::seekg (seek_to_from_command_line()))
    {
      std::cerr << "data_source_and_type::_file_opened_ok: cannot seek to " << seek_to_from_command_line() << '\n';
    }
  }

  return 1;
}

template <typename T>
data_source_and_type<T>::data_source_and_type(FileType input_type,
                                              const char * fname) :
                                      iwstring_data_source (fname)
{
  if (0 == input_type) {
    input_type = discern_file_type_from_name(fname);
    if (0 == input_type)
    {
      std::cerr << "Data_Source_and_Type:: cannot discern input type '" << fname << "'\n";
      return;
    }
  }

  if (!_default_values(input_type)) {
    return;
  }

  if (! _file_opened_ok (fname)) {
    return;
  }

  _do_skip_first();

  return;
}

template <typename T>
data_source_and_type<T>::data_source_and_type(FileType input_type,
                                              const IWString & fname) :
                           iwstring_data_source (fname)
{
  if (! _default_values(input_type))
    return;

  if (! _file_opened_ok(fname))
    return;

  _do_skip_first();

  return;
}

template <typename T>
data_source_and_type<T>::data_source_and_type (const IWString & fname) 
                           : iwstring_data_source (fname)
{
  FileType tmp;
  if (FILE_TYPE_INVALID == (tmp = discern_file_type_from_name (fname)))
  {
    std::cerr << "Data_source_and_type::constructor: cannot determine file type '" << fname << "'\n";
    return;
  }

  if (! _default_values(tmp))
    return;

  if (! _file_opened_ok (fname))
    return;

  _do_skip_first();

  return;
}

template <typename T>
data_source_and_type<T>::data_source_and_type()
{
  _input_type = FILE_TYPE_INVALID;

  _default_values(FILE_TYPE_NO_ARG_CONSTRUCTOR);

  return;
}

template <typename T>
data_source_and_type<T>::~data_source_and_type()
{

  if (_verbose)
    std::cerr << "Read " << _molecules_read << " molecules from '" << _fname << "'\n";

  if (_connection_table_errors_encountered)
    std::cerr << "data_source_and_type:: " << _connection_table_errors_encountered << " connection table errors encountered\n";

  return;
}

template <typename T>
int
data_source_and_type<T>::ok() const
{
  return iwstring_data_source::ok();
}

template <typename T>
int
data_source_and_type<T>::do_open(const char * fname,
                                 FileType input_type)
{
  if (valid_file_type(input_type))
    _input_type = input_type;
  else if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    if (input_type <= 0)
    {
      std::cerr << "data_source_and_type::do_open:cannot discern input type from '" << fname << "'\n";
      return 0;
    }
    _input_type = input_type;
  }
  else
  {
    std::cerr << "data_source_and_type::do_open:invalid file type " << input_type << '\n';
    return 0;
  }

  if (! iwstring_data_source::open(fname))
  {
    std::cerr << "data_source_and_type::do_open:cannot open '" << fname << "'\n";
    return 0;
  }

  if (! _file_opened_ok(fname)) {
    return 0;
  }

  _do_skip_first();

  return 1;
}

template <typename T>
T *
data_source_and_type<T>::next_molecule()
{
//std::cerr << "Reading next molecule\n";
//debug_print(std::cerr);

  if (! _valid)
    return nullptr;

  if (_do_only > 0 && _molecules_read >= _do_only)
  {
    std::cerr << "data_source_and_type::next_molecule: already read " << _molecules_read << " molecules\n";
    _valid = 0;
    return nullptr;
  }

  if (! good())
  {
    _valid = 0;
    return nullptr;
  }

  while (_connection_table_errors_encountered <= _connection_table_errors_allowed)
  {
    if (! is_pipe())   // tellg does not work on stdin
    {
      _offset_for_most_recent_molecule = tellg();

      if (_offset_for_most_recent_molecule >= max_offset_from_command_line())
        return nullptr;
    }

    T * m = new T;
    if (m->read_molecule_ds (*this, _input_type))
    {
      _molecules_read++;
      if (_verbose)
      {
        std::cerr << _molecules_read;
        if (m->name().length())
          std::cerr << " read '" << m->name() << "'\n";
        else
          std::cerr << " no name\n";
      }

      return m;
    }

    if (at_eof())     // normal termination
    {
      delete m;
      return 0;
    }

    _connection_table_errors_encountered++;

    std::cerr << "data_source_and_type::next_molecule: Skipping connection table error " << _connection_table_errors_encountered << 
            " record " << lines_read() << '\n';

    if (m->name().length()) {
      std::cerr << "Molecule name '" << m->name() << "'\n";
    }

    delete m;

    if (_stream_for_connection_table_errors.is_open()) {
      off_t here = tellg();

      size_t bytes_to_echo = here - _offset_for_most_recent_molecule;

//    std::cerr << "From " << here << " go back to " << _offset_for_most_recent_molecule << " echo " << bytes_to_echo << " bytes\n";

      if (! seekg (_offset_for_most_recent_molecule))
        std::cerr << "data_source_and_type::next_molecule:sorry, cannot seek back to " << _offset_for_most_recent_molecule << '\n';
      else
      {
        echo (_stream_for_connection_table_errors, static_cast<off_t>(bytes_to_echo));
        seekg (here);

        std::cerr << "from " << _offset_for_most_recent_molecule << " to " << here << " now " << tellg() << '\n';
      }
    }

    if (_connection_table_errors_encountered > _connection_table_errors_allowed)
    {
      std::cerr << "data_source_and_type::next_molecule:too many connection table errors " << _connection_table_errors_allowed << '\n';
      return 0;
    }
  }

  return nullptr;
}

template <typename T>
void
data_source_and_type<T>::set_connection_table_errors_allowed(int i)
{
  assert (i >= 0);

  _connection_table_errors_allowed = i;

  return;
}

template <typename T>
int
data_source_and_type<T>::set_connection_table_error_file(const IWString & fname)
{
  if (_stream_for_connection_table_errors.is_open())
    _stream_for_connection_table_errors.close();

  IWString tmp(fname);

  if (! _stream_for_connection_table_errors.open(tmp.null_terminated_chars()))
  {
    std::cerr << "data_source_and_type:set_connection_table_error_file: cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

template <typename T>
RE2
data_source_and_type<T>::_rx_for_input_type() const
{
  switch (_input_type)
  {
    case FILE_TYPE_MDL:
    case FILE_TYPE_SDF:
      return RE2("^\\$\\$");

    case FILE_TYPE_SMI:
      return RE2(".");

    case FILE_TYPE_TDT:
      return RE2("^|");
      break;
      
    case FILE_TYPE_MRK:
      return RE2("^[ 0-9]{5} [ 0-9]{5} *$");
      break;

    case FILE_TYPE_PDB:
      return RE2("^END$");
      break;

    default:
      std::cerr << "data_source_and_type::molecules_remaining: no rx for type " << _input_type << '\n';
      return 0;
  }

  return RE2(".");
}

template <typename T>
int
data_source_and_type<T>::molecules_remaining()
{
  if (FILE_TYPE_INVALID == _input_type)
  {
    std::cerr << "data_source_and_type::molecules_remaining: no type specified\n";
    return 0;
  }

  if (FILE_TYPE_SMI == _input_type)
    return records_remaining();

  RE2 rx = _rx_for_input_type();

  return iwstring_data_source::grep(rx);
}

// If _skip_first is specified, skip that many structures.
// Reading csv is a special case.
// TODO: ianwatson read the csv header record and preserve to enable Info collection.
template <typename T>
void
data_source_and_type<T>::_do_skip_first()
{
  assert (_valid);
  if (FILE_TYPE_CSV == _input_type) {
    if (! _skip_csv_header_record()) {
      return;
    }
  }

  if (_skip_first == 0) {
    return;
  }

  if (FILE_TYPE_SMI == _input_type || _input_type == FILE_TYPE_CSV)
  {
    if (! skip_records(_skip_first))
      _valid = 0;

    return;
  }

  if (FILE_TYPE_MRK == _input_type)
  {
    _do_merck_skip_records();
    return;
  }

  RE2 rx = _rx_for_input_type();

  if (! skip_records(rx, _skip_first))
  {
    _valid = 0;

    return;
  }

  return;
}

template <typename T>
int
data_source_and_type<T>::_skip_csv_header_record() {
  const_IWSubstring buffer;
  if (! next_record(buffer)) {
    return 0;
  }

  if (! lillymol_csv::ok_header_record(buffer)) {
    std::cerr << "data_source_and_type::_skip_csv_header_record:invalid csv header '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
data_source_and_type<T>::_do_merck_skip_records()
{
  for (int i = 0; i < _skip_first; i++)
  {
    if (! _skip_merck_connection_table())
    {
      std::cerr << "data_source_and_type::_do_merck_skip_records: fatal error trying to skip " << _skip_first << " records, got " << i << '\n';
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
data_source_and_type<T>::_skip_merck_connection_table()
{
  const_IWSubstring buffer;

  if (! iwstring_data_source::next_record (buffer))
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: premature eof\n";
    return 0;
  }

  if (! iwstring_data_source::next_record (buffer))
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: cannot read 2nd header record\n";
    return 0;
  }

  if (! buffer.starts_with ("MOL "))
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: 2nd record looks invalid, line " << lines_read() << '\n';
    std::cerr << buffer << '\n';
    return 0;
  }

  if (! iwstring_data_source::next_record (buffer))
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: cannot read atom/bond count record\n";
    return 0;
  }

  if (buffer.nwords() < 2)
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: invalid na,nb record, line " << lines_read() << '\n';
    std::cerr << buffer << '\n';
    return 0;
  }

  const_IWSubstring token;
  int i = 0;

  (void) buffer.nextword (token, i);
  int na;

  if (! token.numeric_value (na) || na < 1)
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: invalid number of atoms\n";
    std::cerr << buffer << '\n';
    return 0;
  }

  (void) buffer.nextword (token, i);
  int nb;

  if (! token.numeric_value (nb) || nb < 1)
  {
    std::cerr << "data_source_and_type::_skip_merck_connection_table: invalid number of bonds\n";
    std::cerr << buffer << '\n';
    return 0;
  }

  int bond_records = nb / 5;
  if (0 != nb % 5)
    bond_records++;

//std::cerr << "Skipping " << na << " atoms and " << nb << " bonds or " << (na + bond_records) << " records\n";

  for (int i = 0; i < (na + bond_records); i++)
  {
    if (! data_source_and_type::next_record (buffer))
    {
      std::cerr << "data_source_and_type::_skip_merck_connection_table: eof im middle of molecule\n";
      return 0;
    }
  }

  return 1;
}

/*
  Starting at Offset, read N items and return the average size of those
  items in bytes
*/

template <typename T>
size_t
data_source_and_type<T>::_average_size (off_t offset,
                                     RE2 & rx,
                                     int n)
{
  if (! iwstring_data_source::seekg(offset))
  {
    std::cerr << "data_source_and_type::_average_size:cannot seek to " << offset << '\n';
    return 0;
  }

  assert (n > 0);

  bool found_matching_record = false;
  const_IWSubstring buffer;
  while (next_record(buffer))
  {
    re2::StringPiece tmp(buffer.data(), buffer.length());
    if (RE2::PartialMatch(tmp, rx))
      continue;

    found_matching_record = true;
    break;
  }

  if (! found_matching_record)
    return 0;

  off_t start_pos = iwstring_data_source::tellg();

  int nfound = 0;
  while (next_record(buffer))
  {
    re2::StringPiece tmp(buffer.data(), buffer.length());
    if (RE2::PartialMatch(tmp, rx))
      continue;

    nfound++;
    if (nfound >= n)
      break;
  }

  if (nfound < n)
    return 0;

  off_t end_pos = iwstring_data_source::tellg();

  size_t rc = (end_pos - start_pos) / n;

  return rc;
}

template <typename T>
int
data_source_and_type<T>::estimate_molecules_in_file()
{
  RE2 rx = _rx_for_input_type();

  size_t fsize = iwstring_data_source::file_size();

  if (fsize < 2000000)
    return iwstring_data_source::grep(rx);

// Get estimates from start, middle and end of file

  size_t start_of_file = _average_size(0, rx, 1000);
  size_t middle_of_file = _average_size(fsize / 2, rx, 1000);
  size_t end_of_file = _average_size(fsize - 2000 * middle_of_file, rx, 1000);

  std::cerr << "Sizes: start " << start_of_file << " middle " << middle_of_file << " and " << end_of_file << '\n';

  size_t ave = (start_of_file + middle_of_file + end_of_file) / 3;

  int items_in_file = fsize / ave;

  return items_in_file;
}

#endif

#endif


/* arch-tag: 2f164a46-03e8-45fa-b4a0-5b4e4724503d

*/
