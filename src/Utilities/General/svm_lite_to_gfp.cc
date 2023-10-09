/*
  Converts svm lite format to .gfp form
*/

#include <stdlib.h>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static IWString tag("NCFC<");

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int lines_read = 0;

static Sparse_Fingerprint_Creator sfpc;

static char input_separator = ' ';

static int response_in_column_one = 1;

static Accumulator_Int<unsigned int> acc_nbits;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Converts svm-lite form to .gfp form\n";
  cerr << "    feature:number feature:number feature:number .... # ID1\n";
  cerr << "    feature:number feature:number feature:number .... # ID1\n";
  cerr << " -J <tag>       tag to create\n";
  cerr << " -S <fname>     merge in smiles information from <fname>\n";
  cerr << " -i <sep>       input file separator\n";
  cerr << " -G             input file does NOT have a response as the first column\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_write(const IWString& id, const IWString& smiles, Sparse_Fingerprint_Creator& sfpc,
         IWString_and_File_Descriptor& output)
{
  output << identifier_tag << id << ">\n";
  if (smiles.length()) {
    output << smiles_tag << smiles << ">\n";
  }

  IWString tmp;

  sfpc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << "\n";

  output << "|\n";

  if (verbose) {
    acc_nbits.extra(sfpc.nbits());
  }

  return 1;
}

static int
process_bit(const const_IWSubstring& token, Sparse_Fingerprint_Creator& sfpc)
{
  const_IWSubstring bs, cs;

  if (!token.split(bs, ':', cs) || 0 == bs.length() || 0 == cs.length()) {
    cerr << "process_bit:cannot split into bit:count '" << token << "'\n";
    return 0;
  }

  unsigned int b, c;
  if (!bs.numeric_value(b) || !cs.numeric_value(c)) {
    cerr << "process_bit::invalid numeric '" << token << "'\n";
    return 0;
  }

  // cerr << "Setting bit " << b << " " << c << " times\n";
  sfpc.hit_bit(b, c);

  return 1;
}

static int
get_smiles(const IW_STL_Hash_Map_String& smiles, const IWString& id, IWString& mysmiles)
{
  const auto f = smiles.find(id);

  if (f == smiles.end()) {
    return 0;
  }

  mysmiles = f->second;

  return 1;
}

static int
svm_lite_to_gfp_record(const const_IWSubstring& buffer,
                       const IW_STL_Hash_Map_String& smiles,
                       IWString_and_File_Descriptor& output)
{
  int i = 0;

  Sparse_Fingerprint_Creator sfpc;
  const_IWSubstring token;

  if (response_in_column_one) {
    if (!buffer.nextword(token, i, input_separator)) {
      return 0;
    }
  }

  while (buffer.nextword(token, i, input_separator)) {
    if ('#' == token) {
      break;
    }

    if (!process_bit(token, sfpc)) {
      cerr << "Cannot process '" << token << "'\n";
      return 0;
    }
  }

  IWString id;
  IWString mysmiles;

  if (buffer.nextword(id, i, input_separator)) {
    if (smiles.size()) {
      if (!get_smiles(smiles, id, mysmiles)) {
        cerr << "svm_lite_to_gfp_record:cannot get smiles for '" << id << "'\n";
        return 0;
      }
    }
  }

  return do_write(id, mysmiles, sfpc, output);
}

static int
svm_lite_to_gfp(iwstring_data_source& input, const IW_STL_Hash_Map_String& smiles,
                IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    lines_read++;

    if (!svm_lite_to_gfp_record(buffer, smiles, output)) {
      cerr << "Fatal error processing\n";
      cerr << buffer << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
svm_lite_to_gfp(const char* fname, const IW_STL_Hash_Map_String& smiles,
                IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return svm_lite_to_gfp(input, smiles, output);
}

static int
read_smiles_record(const const_IWSubstring& buffer, IW_STL_Hash_Map_String& smiles_hash)
{
  IWString smiles, pcn;

  if (!buffer.split(smiles, ' ', pcn)) {
    cerr << "Cannot split line into smiles and id\n";
    return 0;
  }

  pcn.truncate_at_first(' ');

  if (verbose > 2) {
    cerr << "Reading smiles identifier '" << pcn << "'\n";
  }

  if (smiles_hash.contains(pcn)) {
    cerr << "Ignoring duplicate identifier '" << pcn << "'\n";
    return 1;
  }

  smiles_hash[pcn] = smiles;

  return smiles_hash.size();
}

static int
read_smiles(iwstring_data_source& input, IW_STL_Hash_Map_String& smiles_hash)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!read_smiles_record(buffer, smiles_hash)) {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return smiles_hash.size();
}

static int
read_smiles(const const_IWSubstring& fname, IW_STL_Hash_Map_String& smiles)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles(input, smiles);
}

static int
svm_lite_to_gfp(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vJ:S:i:G");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('G')) {
    response_in_column_one = 0;

    if (verbose) {
      cerr << "First column is NOT the response\n";
    }
  }

  if (cl.option_present('i')) {
    IWString tmp = cl.string_value('i');
    if (!char_name_to_char(tmp)) {
      cerr << "Unrecognised input separator specification '" << tmp << "'\n";
      return 1;
    }

    input_separator = tmp[0];

    if (verbose) {
      cerr << "Input file separator set to " << input_separator << endl;
    }
  }

  IW_STL_Hash_Map_String smiles;

  if (cl.option_present('S')) {
    const char* s = cl.option_value('S');

    if (!read_smiles(s, smiles)) {
      cerr << "Cannot read smiles from '" << s << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Read " << smiles.size() << " smiles from '" << s << "'\n";
    };
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!svm_lite_to_gfp(cl[i], smiles, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Processed " << lines_read << " lines\n";
    if (acc_nbits.n() > 0) {
      cerr << "Fingerprints had btw " << acc_nbits.minval() << " and "
           << acc_nbits.maxval() << " ave " << static_cast<float>(acc_nbits.average())
           << " bits\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = svm_lite_to_gfp(argc, argv);

  return rc;
}
