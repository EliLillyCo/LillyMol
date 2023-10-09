/*
  Converts the output of parallel_nn_search to the input for
  gfp_spread
*/

#include <stdlib.h>

#include <iostream>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

typedef IW_STL_Hash_Map<IWString, IW_TDT*> ID_to_TDT;

static IWString identifier_tag("PCN<");
static IWString smiles_tag("$SMI<");
static IWString distance_tag("DIST<");

static IWString nnsmi_tag("NNSMI<");
static IWString nnid_tag("NNID<");
static IWString nndist_tag("NNDIST<");

static int items_read = 0;
static int items_written = 0;

static Accumulator<double> first_nbr_dist;

static float lower_distance_threshold = static_cast<float>(0.0);

static int items_closer_than_lower_distance_threshold = 0;

static int need_to_convert_string_distance_to_float = 0;

static int molecules_with_no_nbr_data = 0;

static int take_first_token_of_name = 0;

static int remove_leading_zeros = 0;

static int ignore_identifiers_not_in_fingerprint_file = 0;

static int missing_identifiers_ignored = 0;

/*
  Note we include MPR in the fingerprint rx....
*/

static std::unique_ptr<re2::RE2> fingerprint_rx;

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
  cerr << "Converts a .nn file to proto form\n";
  cerr << " -o <sep>      set token separator (default ,)\n";
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts the output of parallel_nn_search to the -A input file for gfp_spread\n";
  cerr << " -F <fname>     fingerprint file from which splits were made\n";
  cerr << " -t <dist>      discard molecules with distances less than <dist>\n";
  cerr << " -j             only consider first token of all name fields\n";
  cerr << " -z             remove leading 0's from identifiers\n";
  cerr << " -X             just skip molecules for which there is no identifier in the -F file\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
add_fingerprint_to_hash(IW_TDT* tdt, ID_to_TDT& fp)
{
  IWString id;

  if (!tdt->dataitem_value(identifier_tag, id)) {
    cerr << "Cannot extract identifier '" << identifier_tag << "'\n";
    return 0;
  }

  if (take_first_token_of_name) {
    id.truncate_at_first(' ');
  }

  if (remove_leading_zeros) {
    id.remove_leading_chars('0');
  }

  fp[id] = tdt;

  return 1;
}

static int
read_fingerprints(iwstring_data_source& input, ID_to_TDT& fp)
{
  IW_TDT tdt;

  while (1) {
    IW_TDT* tdt = new IW_TDT;

    if (!tdt->next(input))  // no next tdt, might be normal eof
    {
      delete tdt;
      if (input.eof()) {
        return fp.size();
      } else {
        return 0;
      }
    }

    if (!add_fingerprint_to_hash(tdt, fp)) {
      cerr << "Fatal error processing fingerprint\n";
      cerr << *tdt;
      return 0;
    }
  }

  return fp.size();
}

static int
read_fingerprints(const char* fname, ID_to_TDT& fp)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  return read_fingerprints(input, fp);
}

static int
dataitem_to_value(const const_IWSubstring& buffer, int len_tag, IWString& s)
{
  assert(buffer.ends_with('>'));

  s = buffer;
  s.chop();

  s.remove_leading_chars(len_tag);

  return 1;
}

static int
insert_fingerprints(const IW_TDT& tdt, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  int rc = 0;
  int n = tdt.number_elements();

  for (int i = 0; i < n; i++) {
    tdt.item(i, buffer);

    //  cerr << "Checking '" << buffer << "'\n";

    if (!iwre2::RE2PartialMatch(buffer, *fingerprint_rx)) {
      continue;
    }

    if (buffer.ends_with('\n')) {
      buffer.chop();
    }

    output << buffer << '\n';
    rc++;
  }

  return rc;
}

/*
Note that we can read both the target and the nbr this way
*/

static int
read_next_molecule(iwstring_data_source& input, IWString& smiles, IWString& pcn)
{
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    return 0;
  }

  if ('|' == buffer) {
    return 0;
  }

  if (!buffer.starts_with(smiles_tag)) {
    cerr << "Input TDT's must begin with smiles, not '" << buffer << "'\n";
    return 0;
  }

  dataitem_to_value(buffer, smiles_tag.length(), smiles);

  if (!input.next_record(buffer)) {
    cerr << "Premature EOF\n";
    return 0;
  }

  if (!buffer.starts_with(identifier_tag)) {
    cerr << "Second line of input tdts must be pcn, not '" << buffer << "'\n";
    return 0;
  }

  return dataitem_to_value(buffer, identifier_tag.length(), pcn);
}

static int
parallel_nn_search_to_gfp_spread(iwstring_data_source& input,
                                 const IWString& target_smiles,
                                 const IWString& target_pcn, ID_to_TDT& id_to_tdt,
                                 IWString_and_File_Descriptor& output)
{
  IWString nbr_smiles, nbr_pcn;

  if (!read_next_molecule(input, nbr_smiles, nbr_pcn)) {
    //  cerr << "Could not read nbr data for '" << target_pcn << "'\n";
    molecules_with_no_nbr_data++;
    return 1;
  }

  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Premature eof\n";
    return 0;
  }

  if (!buffer.starts_with(distance_tag)) {
    cerr << "Following nbr must find distance tag, not '" << buffer << "'\n";
    return 0;
  }

  IWString dist;

  if (!dataitem_to_value(buffer, distance_tag.length(), dist)) {
    return 0;
  }

  input.skip_past("|");

  if (need_to_convert_string_distance_to_float) {
    float d;
    if (!dist.numeric_value(d) || d < static_cast<float>(0.0) ||
        d > static_cast<float>(1.0)) {
      cerr << "Invalid distance '" << buffer << "'\n";
      return 0;
    }

    if (d < lower_distance_threshold) {
      items_closer_than_lower_distance_threshold++;
      return 1;
    }

    if (verbose) {
      first_nbr_dist.extra(d);
    }
  }

  IW_TDT* fp = nullptr;

  if (id_to_tdt.size()) {
    ID_to_TDT::const_iterator f = id_to_tdt.find(target_pcn);

    if (f != id_to_tdt.end()) {
      ;
    } else if (ignore_identifiers_not_in_fingerprint_file) {
      missing_identifiers_ignored++;
      return 1;
    } else {
      cerr << "No fingerprint for '" << target_pcn << "'\n";
      return 0;
    }

    fp = (*f).second;
  }

  output << smiles_tag << target_smiles << ">\n";
  output << identifier_tag << target_pcn << ">\n";

  output << nnsmi_tag << nbr_smiles << ">\n";
  output << nnid_tag << nbr_pcn << ">\n";
  output << nndist_tag << dist << ">\n";

  if (NULL != fp) {
    insert_fingerprints(*fp, output);
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(16384);

  items_written++;

  return 1;
}

static int
parallel_nn_search_to_gfp_spread(iwstring_data_source& input, ID_to_TDT& id_to_tdt,
                                 IWString_and_File_Descriptor& output)
{
  IWString target_smiles, target_pcn;

  while (read_next_molecule(input, target_smiles, target_pcn)) {
    items_read++;

    if (take_first_token_of_name) {
      target_pcn.truncate_at_first(' ');
    } else if (remove_leading_zeros) {
      target_pcn.remove_leading_chars('0');
    } else {
      target_pcn.truncate_at_last(' ');
    }

    if (!parallel_nn_search_to_gfp_spread(input, target_smiles, target_pcn, id_to_tdt,
                                          output)) {
      cerr << "Fatal error processing data for '" << target_pcn << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
parallel_nn_search_to_gfp_spread(const char* fname, ID_to_TDT& id_to_tdt,
                                 IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return parallel_nn_search_to_gfp_spread(input, id_to_tdt, output);
}

static int
parallel_nn_search_to_gfp_spread(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vF:t:jzX");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (verbose) {
    need_to_convert_string_distance_to_float = 1;
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('z')) {
    remove_leading_zeros = 1;

    if (verbose) {
      cerr << "Will strip leading 0's from identifiers\n";
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', lower_distance_threshold) ||
        lower_distance_threshold < static_cast<float>(0.0)) {
      cerr << "The lower distance threshold (-t) option must be a valid distance\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will discard items having nearest neighbour distances less than "
           << lower_distance_threshold << endl;
    }

    need_to_convert_string_distance_to_float = 1;
  }

  IWString rx("^(FP[A-Z,0-9]+|HX[A-Z,0-9]+|NC[A-Z,0-9]+|FC[A-Z,0-9]+|MC[A-Z,0-9]+|MPR)<");
  iwre2::RE2Reset(fingerprint_rx, rx);

  ID_to_TDT id_to_tdt;

  if (cl.option_present('F')) {
    const char* f;
    for (int i = 0; (f = cl.option_value('F', i)); i++) {
      if (!read_fingerprints(f, id_to_tdt)) {
        cerr << "Cannot read fingerprints '" << f << "'\n";
        return 3;
      }
    }

    if (verbose) {
      cerr << "Read " << id_to_tdt.size() << " fingerprints\n";
    }
  }

  if (cl.option_present('j')) {
    take_first_token_of_name = 1;

    if (verbose) {
      cerr << "Will only use the first token of all name fields (in -F file and in "
              "input)\n";
    }
  }

  if (cl.option_present('X')) {
    ignore_identifiers_not_in_fingerprint_file = 1;

    if (verbose) {
      cerr << "Will skip molecules where there is no fingerprint in the -F file\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!parallel_nn_search_to_gfp_spread(cl[i], id_to_tdt, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << items_read << " items";
    if (cl.option_present('t')) {
      cerr << ", wrote " << items_written;
    }
    cerr << endl;

    if (items_written > 0) {
      cerr << "Distances between " << first_nbr_dist.minval() << " and "
           << first_nbr_dist.maxval() << " ave "
           << static_cast<float>(first_nbr_dist.average_if_available_minval_if_not())
           << endl;
    }

    if (ignore_identifiers_not_in_fingerprint_file && missing_identifiers_ignored) {
      cerr << "Skipped " << missing_identifiers_ignored
           << " molecules with missing fingerprints in the -F file\n";
    }
  }

  if (molecules_with_no_nbr_data) {
    cerr << "Warning, " << molecules_with_no_nbr_data << " had no nbr data\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = parallel_nn_search_to_gfp_spread(argc, argv);

  return rc;
}
