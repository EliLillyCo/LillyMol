// Reads a TDT file and when certain values are encountered, looks up
// an identifier in a Berkeley DB database and inserts the value.
// This was primarily used for inserting inventory data into a .gfp
// file. That use case is no longer of interest, and this will
// be deprecated soon.

#include <cstdlib>
#include <iostream>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "db_cxx.h"

const char* prog_name = NULL;

static IWString input_tag("PCN<");

static IWString output_tag;

static int verbose = 0;

static Db** database = NULL;
static int number_databases = 0;

static int tdts_read = 0;

static int identifiers_processed = 0;

static int items_inserted = 0;

static int items_without_data = 0;

static int strip_leading_zeros_from_identifiers = 0;

static int pad_identifier_with_leading_zeros = 0;

static IWString missing_value_string;

static int ignore_missing_data = 0;

static int tdts_without_identifiers = 0;

static int special_processing_for_lilly_inventory = 0;

static int special_processing_for_demerit_data = 0;

static int remove_existing_data = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << DB_VERSION_STRING << endl;

  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Inserts values from a Berkeley database into a TDT\n";
  cerr << DB_VERSION_STRING << endl;
  cerr << " -d <dbname>    name of database\n";
  cerr << " -I <tag>       input tag\n";
  cerr << " -O <tag>       write results with tag <tag>\n";
  cerr << " -z             strip leading zero's from identifiers before doing lookup\n";
  cerr << " -p <width>     pad identifiers to <width> characters with leading 0's\n";
  cerr << " -g             ignore identifiers with no data in the database\n";
  cerr << " -M <string>    what to insert if no value in database\n";
  cerr << " -P inv         special processing for Lilly inventory\n";
  cerr << " -P dmrt        special processing for demerit data\n";
  cerr << " -r             remove any existing data in the input\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_the_insertion_lilly_inventory(const IWString& to_insert, IWString& output) {
  const_IWSubstring tmp(to_insert);
  tmp.truncate_at_first(' ');

  double d;
  if (!tmp.numeric_value(d)) {
    cerr << "Invalid numeric, cannot do conversion '" << to_insert << "'\n";
    return 0;
  }

  if (d < 0.0) {
    return 1;
  }

  output << output_tag << static_cast<int>(d * 1000.0 + 0.4999) << ">\n";

  return 1;
}

static int
do_the_insertion_demerit_data(const IWString& to_insert, IWString& output) {
  assert(to_insert.starts_with('D'));

  const_IWSubstring tmp(to_insert);

  tmp.remove_leading_chars(1);

  output << output_tag << tmp << ">\n";

  return 1;
}

static int
do_the_insertion(const IWString& to_insert, IWString& output) {
  if (special_processing_for_lilly_inventory) {
    return do_the_insertion_lilly_inventory(to_insert, output);
  }

  if (special_processing_for_demerit_data) {
    return do_the_insertion_demerit_data(to_insert, output);
  }

  output << output_tag << to_insert << ">\n";

  return 1;
}

static int
iwbdb_merge_into_tdt_item_1(const const_IWSubstring& id, IWString& output) {
  Dbt dbkey(const_cast<char*>(id.rawchars()), id.length());

  for (int i = 0; i < number_databases; i++) {
    Dbt fromdb;
    int rc = database[i]->get(NULL, &dbkey, &fromdb, 0);

    if (0 != rc) {
      continue;
    }

    IWString to_insert;
    to_insert.set(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
    items_inserted++;
    return do_the_insertion(to_insert, output);
  }

  items_without_data++;

  if (missing_value_string.length()) {
    output << missing_value_string;
  }

  if (ignore_missing_data) {
    return 1;
  }

  cerr << "No data for '" << id << "'\n";

  return 0;
}

static int
iwbdb_merge_into_tdt_item(const_IWSubstring& buffer, IWString& output) {
  const_IWSubstring id(buffer);

  id.remove_leading_chars(input_tag.length());
  id.chop();

  if (strip_leading_zeros_from_identifiers && id.starts_with('0')) {
    id.remove_leading_chars('0');
    return iwbdb_merge_into_tdt_item_1(id, output);
  }

  if (pad_identifier_with_leading_zeros &&
      id.length() < pad_identifier_with_leading_zeros) {
    IWString tmp;

    return 1;
  }

  return iwbdb_merge_into_tdt_item_1(id, output);
}

static int
iwbdb_merge_into_tdt(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring input_buffer;

  int got_identifier_this_tdt = 0;

  while (input.next_record(input_buffer)) {
    if (remove_existing_data && input_buffer.starts_with(output_tag)) {
      continue;
    }

    output << input_buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    //  cerr << "Read record of length " << input_buffer.length() << " buffer " <<
    //  output.length() << ", " << tdts_read << " tdts\n";

    if ('|' == input_buffer) {
      tdts_read++;
      if (!got_identifier_this_tdt) {
        tdts_without_identifiers++;
      }
      got_identifier_this_tdt = 0;
      continue;
    }

    if (!input_buffer.starts_with(input_tag)) {
      continue;
    }

    identifiers_processed++;

    got_identifier_this_tdt = 1;
    if (!iwbdb_merge_into_tdt_item(input_buffer, output)) {
      cerr << "Failure doing insert, line " << input.lines_read() << endl;
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_merge_into_tdt(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return iwbdb_merge_into_tdt(input, output);
}

static int
iwbdb_merge_into_tdt(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:I:O:zM:p:gP:r");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('d')) {
    cerr << "Must specify database with the -d option\n";
    usage(2);
  }

  if (!cl.option_present('O')) {
    cerr << "Must specify the output tag via the -O option\n";
    usage(4);
  }

  number_databases = cl.option_count('d');

  if (number_databases) {
    database = new Db*[number_databases];
    int env_flags = DB_CXX_NO_EXCEPTIONS;

    int i = 0;
    IWString dbname;
    while (cl.value('d', dbname, i)) {
      database[i] = new Db(NULL, env_flags);

      int rc = database[i]->open(NULL, dbname.null_terminated_chars(), NULL, DB_UNKNOWN,
                                 DB_RDONLY, 0);

      if (0 != rc) {
        cerr << "Cannot open database '" << dbname << "' '";
        database[i]->err(rc, "");
        return 0;
      }

      if (verbose) {
        cerr << "database " << i << " is '" << dbname << "'\n";
      }

      i++;
    }
  }

  if (cl.option_present('I')) {
    cl.value('I', input_tag);
    if (verbose) {
      cerr << "Input tag '" << input_tag << "'\n";
    }

    if (!input_tag.ends_with('<')) {
      input_tag << '<';
    }
  }

  if (cl.option_present('O')) {
    cl.value('O', output_tag);

    if (verbose) {
      cerr << "Output tag '" << output_tag << "'\n";
    }

    if (!output_tag.ends_with('<')) {
      output_tag << '<';
    }
  }

  if (cl.option_present('r')) {
    remove_existing_data = 1;
    if (verbose) {
      cerr << "Existing data removed\n";
    }
  }

  if (cl.option_present('g')) {
    ignore_missing_data = 1;

    if (verbose) {
      cerr << "Will ignore identifiers with no data\n";
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m;
    cl.value('M', m);

    if (verbose) {
      cerr << "Missing data written s '" << m << "'\n";
    }

    missing_value_string << output_tag << m << ">\n";
    ignore_missing_data = 1;
  }

  if (cl.option_present('z')) {
    strip_leading_zeros_from_identifiers = 1;
    if (verbose) {
      cerr << "Leading 0's will be stripped from identifiers\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', pad_identifier_with_leading_zeros) ||
        pad_identifier_with_leading_zeros < 2) {
      cerr << "The pad with leading zero's value(-p) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will pad identifiers with leading 0's to width "
           << pad_identifier_with_leading_zeros << endl;
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (p == "inv") {
      special_processing_for_lilly_inventory = 1;
      if (verbose) {
        cerr << "Special processing for Lilly inventory\n";
      }
    } else if ("dmrt" == p) {
      special_processing_for_demerit_data = 1;
      if (verbose) {
        cerr << "Special processing for demerit data\n";
      }
    } else {
      cerr << "Unrecognised special processing directive '" << p << "'\n";
      usage(4);
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!iwbdb_merge_into_tdt(cl[i], output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  for (int i = 0; i < number_databases; i++) {
    database[i]->close(0);
  }

  delete[] database;

  if (verbose) {
    cerr << "Read " << tdts_read << " tdts, processed " << identifiers_processed
         << " identifiers\n";
    cerr << items_without_data << " items had no data, " << items_inserted
         << " items inserted\n";
    if (tdts_without_identifiers) {
      cerr << "Warning, " << tdts_without_identifiers << " tdts had no identifier\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = iwbdb_merge_into_tdt(argc, argv);

  return rc;
}
