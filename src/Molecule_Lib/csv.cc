// I/O with smiles CSV form.

#include <iostream>

#include "istream_and_type.h"
#include "molecule.h"
#include "rwmolecule.h"

using std::cerr;
using std::endl;

namespace lillymol_csv {

char csv_separator = ',';

void
set_csv_separator(char s) {
  csv_separator = s;
}

int smiles_column = 0;
int id_column = 1;

// Columns can also be identified by name

IWString smiles_column_name;
IWString id_column_name;

int
set_smiles_column(int smi_col) {
  if (smi_col < 0) {
    cerr << "set_smiles_col:invalid smiles column " << smi_col << "\n";
    return 0;
  }

  smiles_column = smi_col;
  return 1;
}

int
set_id_column(int id_col) {
  if (id_col < 0) {
    cerr << "set_id_col:invalid id column " << id_col << "\n";
    return 0;
  }

  id_column = id_col;

  return 1;
}

void
set_smiles_column_name(const const_IWSubstring& s) {
  smiles_column_name = s;
  smiles_column = -1;
}

void
set_id_column_name(const const_IWSubstring& s) {
  id_column_name = s;
  id_column = -1;
}

int
ok_header_record(const const_IWSubstring& header) {
  const_IWSubstring token;
  int i = 0;
  for (int col = 0; header.nextword(token, i, csv_separator); ++col) {
    if (token == smiles_column_name) {
      smiles_column = col;
    } else if (token == id_column_name) {
      id_column = col;
    }
  }
  if (smiles_column < 0) {
    cerr << "lillymol_csv::ok_header_record:no match for " << smiles_column_name << " in '" << header << "'\n";
    return 0;
  }

  return 1;
}

// It can be useful to concatenate all the other tokens into the name.
int all_tokens_to_name = 0;
void set_all_tokens_to_name(int s) {
   all_tokens_to_name = s;

   if (all_tokens_to_name) {
     id_column = -1;
   } else {
     id_column = 1;
   }
}

}  // namespace lillymol_csv

int
Molecule::write_molecule_csv(std::ostream& output) {
  output << smiles() << lillymol_csv::csv_separator << name() << '\n';

  return output.good();
}

int
Molecule::read_molecule_csv_ds(iwstring_data_source & input) {
  if (input.eof())
    return 0;

  const_IWSubstring buffer;
  EXTRA_STRING_RECORD(input, buffer, "read mol csv");

  const_IWSubstring token;
  const_IWSubstring smi;
  int i = 0;
  for (int col = 0; buffer.nextword(token, i, lillymol_csv::csv_separator); ++col) {
    if (col == lillymol_csv::smiles_column) {
      smi = token;
      if (_molecule_name.length() > 0 && ! lillymol_csv::all_tokens_to_name) {
        break;
      }
    } else if (col == lillymol_csv::id_column) {
      _molecule_name = token;
      if (smi.length() > 0) {
        break;
      }
    } else if (lillymol_csv::all_tokens_to_name) {
      _molecule_name.append_with_spacer(token, lillymol_csv::csv_separator);
    }
  }

  if (smi.empty()) {
    cerr << "Molecule::read_molecule_csv_ds:no smiles in column " << lillymol_csv::smiles_column << " '" << buffer << "'\n";
    return 0;
  }

  if (! build_from_smiles(smi)) {
    cerr << "Cannot parse smiles '" << smi << "'\n";
    return 0;
  }

  return 1;
}

// included in this file because the file is not big.
namespace lillymol_textproto {

IWString smiles_tag("smi:");

void
set_smiles_tag(const_IWSubstring s) {
  smiles_tag = s;
  if (! smiles_tag.ends_with(':')) {
    smiles_tag << ':';
  }
}

}  // namespace lillymol_textproto

int
Molecule::write_molecule_textproto(std::ostream& output) {
  output << lillymol_textproto::smiles_tag << " \"" << smiles() <<
         "\" id: \"" << _molecule_name << "\"\n";
  return output.good();
}

int
Molecule::read_molecule_textproto_ds(iwstring_data_source & input) {
  if (input.eof()) {
    return 0;
  }

  const_IWSubstring buffer;
  EXTRA_STRING_RECORD(input, buffer, "read mol textproto");

  static constexpr char kSpace = ' ';

  const_IWSubstring smi;

  const_IWSubstring token;
  bool next_is_smiles = false;
  for (int i = 0; buffer.nextword(token, i, kSpace); ) {
    if (next_is_smiles) {
      smi = token;
      next_is_smiles = false;
    } else if (token == lillymol_textproto::smiles_tag) {
      next_is_smiles = true;
    } else {
      _molecule_name.append_with_spacer(token);
    }
  }
  
  if (smi.empty()) {
    cerr << "Molecule::read_molecule_textproto_ds:no smiles '" << buffer << "'\n";
    return 0;
  }

  // Very likely written as "C"
  if (smi.starts_with('"')) {
    smi++;
  }
  if (smi.ends_with('"')) {
    smi.chop();
  }

  if (! build_from_smiles(smi)) {
    cerr << "Cannot parse smiles '" << smi << "'\n";
    return 0;
  }

  return 1;
}
