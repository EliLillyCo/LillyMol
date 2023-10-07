/*
  Tester for the NEED/AVAILABLE substructure queries
*/

#include <iostream>
#include <memory>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"

#include "nass.h"

using std::cerr;
using std::endl;

static const char * prog_name = nullptr;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Does substructure searches with NEED/AVAILABLE flags\n";

  cerr << prog_name << " <options> <input_file1> <input_file2> ... > cout\n";
  cerr << "  -m <fname>     write molecules matching at least one query to <fname>\n";
  cerr << "  -D             append matching details to the names\n";
  cerr << "  -n <fname>     write non matching molecules to <fname>\n";
  cerr << "  -r             ignore queries with AVAILABLE features when considering matches\n";
  cerr << "  -b             for each molecule, break after finding a query which matches\n";
  cerr << "  -a             write as array (0,1). Repeat to write counts\n";
  cerr << "  -Y <name>      descriptor names start with <name>\n";
  cerr << "  -J <tag>       specify fingerprint dataitem\n";
  cerr << "  -y <nbits>     creation size of fingerprint in bits\n";
  cerr << "  -u             find unique matches only\n";
  cerr << "  -G <fname>     file for lists of matched atoms\n";
  cerr << "  -F             function as filter\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -E autocreate  automatically create new elements when encountered\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_etrans_options(cerr, 't');
  (void) display_standard_charge_assigner_options(cerr, 'N');
  cerr << "  -H <...>       donor acceptor options, enter '-H help' for details\n";
  cerr << "  -W <bytes>     store as bytes\n";
  cerr << "  -W <sparse_bytes> store as sparse bytes, non zero results only\n";
  cerr << "  -k             special processing for CNK descriptor\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -X ...         miscellaneous options, enter '-Y help' for info\n";
  cerr << "  -x <nhits>     when producing descriptors, limit to <nhits> occurrences\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int verbose = 0;

static int molecules_read = 0;

static int molecules_matching = 0;

static int * hits_to_query;

static int append_query_details_to_name = 0;

static int reduce_to_largest_fragment = 0;

static int write_as_array = 0;

static int is_filter = 0;

static IWString fingerprint_tag;

static IWString non_colliding_fingerprint_tag;

static int default_fingerprint_nbits;

static Chemical_Standardisation chemical_standardisation;

static Element_Transformations element_transformations;

static Charge_Assigner charge_assigner;

static Molecule_Output_Object stream_for_non_matching_structures;
static Molecule_Output_Object stream_for_matching_structures;

static Set_of_NA_Substructure_Query queries;

static int * results = nullptr;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static int * hits_per_molecule = nullptr;

static int max_hits = std::numeric_limits<int>::max();

static IWDigits iwdigits;

static int special_processing_for_cnk = 0;

/*
  Sometimes the AVAILABLE queries are just there for efficiency, so we ignore
  matches to them when considering results
*/

static int ignore_available_queries = 0;

static IWString_and_File_Descriptor bob_coner_stream;

static int flush_after_each_molecule = 0;

static void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(8192);
  }
}

static int
write_bob_coner_special (const IWString & mname)
{
  append_first_token_of_name(mname, bob_coner_stream);

  bob_coner_stream << ' ';

  const Substructure_Results * sresults = queries.sresults();

  for (int i = 0; i < queries.number_elements(); i++)
  {
    if (i > 0)
      bob_coner_stream << ',';

    int nhits = results[i];

    const Substructure_Results & s = sresults[i];

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = s.embedding(j);

      bob_coner_stream << '(';
      for (int k = 0; k < e->number_elements(); k++)
      {
        if (k > 0)
          bob_coner_stream << ' ';

        bob_coner_stream << e->item(k);
      }
      bob_coner_stream << ')';
    }
  }

  bob_coner_stream << '\n';

  bob_coner_stream.write_if_buffer_holds_more_than(32768);

  return bob_coner_stream.good();
}

static int
tnass (Molecule & m)
{
  if (! queries.substructure_search(m, results))    // catastrophic error
    return 0;

  int queries_matching = 0;

  IWString new_name;
  if (append_query_details_to_name)
  {
    new_name.resize(m.name().length() + 1000);
    new_name = m.name();
  }

  int nq = queries.number_elements();
  for (int i = 0; i < nq; i++)
  {
    if (verbose > 1)
      cerr << "Results from query " << i << " are " << results[i] << endl;

    if (results[i] <= 0)
    {
      results[i] = 0;
      continue;
    }

    const NA_Substructure_Query & q = *(queries[i]);

    if (ignore_available_queries && q.makes_available().length())    // this query is just for efficiency
    {
      results[i] = 0;
      continue;
    }

    if (append_query_details_to_name)
    {
      new_name += " (";
      new_name << results[i] << " matches to '" << queries[i]->comment();
      new_name += "')";
    }

    queries_matching++;

    if (verbose > 2)
      cerr << m.name() << " " << results[i] << " matches to query " << i << " '" << queries[i]->comment() << "'\n";
  }

  if (verbose > 1)
    cerr << m.name() << " matched " << queries_matching << " queries\n";

  if (verbose)
    hits_per_molecule[queries_matching]++;

  if (append_query_details_to_name)
    m.set_name(new_name);

  if (queries_matching)
    molecules_matching++;

  if (verbose)
  {
    for (int i = 0; i < nq; i++)
    {
      if (results[i])
        hits_to_query[i]++;
    }
  }

  if (queries_matching && stream_for_matching_structures.active())
    stream_for_matching_structures.write(m);
  else if (0 == queries_matching && stream_for_non_matching_structures.active())
    stream_for_non_matching_structures.write(m);

  if (bob_coner_stream.is_open())
    write_bob_coner_special(m.name());

  return 1;
}

/*
  In a very arbitrary move, we capture the unique smiles before the charge
  assigner and donor-acceptor assigner.
*/

static int
preprocess_molecule (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (element_transformations.active())
    (void) element_transformations.process(m);

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (charge_assigner.active())
    (void) charge_assigner.process(m);

  if (donor_acceptor_assigner.active())
    donor_acceptor_assigner.process(m);

  return 1;
}

static int
write_fingerprint (IWString & output_buffer)
{
  if (fingerprint_tag.length())
  {
    IW_Bits_Base fp(default_fingerprint_nbits);

    fp.construct_from_array_of_ints(results, default_fingerprint_nbits);

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output_buffer << fingerprint_tag << tmp << ">\n";
  }
  else if (non_colliding_fingerprint_tag.length())
  {
    Sparse_Fingerprint_Creator fp;

    fp.create_from_array_of_ints(results, default_fingerprint_nbits);

    IWString tmp;
    fp.daylight_ascii_form_with_counts_encoded(tmp);
    output_buffer << non_colliding_fingerprint_tag << tmp << ">\n";
  }
  else
  {
    cerr << "tnass::write_fingerprint: what am I supposed to be doing\n";
    abort();
  }

  return 0;
}

static int
do_output (Molecule & m,
           IWString & output)
{
  int nq = queries.number_elements();

  if (fingerprint_tag.length() || non_colliding_fingerprint_tag.length())
  {
    output << "$SMI<" << m.smiles() << ">\n";
    output << "PCN<" << m.name() << ">\n";

    write_fingerprint(output);

    output << "|\n";

    return 1;
  }

  if (0 == write_as_array)
    return 1;

  append_first_token_of_name(m.name(), output);

  if (1 == write_as_array)
  {
    for (int i = 0; i < nq; i++)
    {
      if (0 == results[i])
        output << " 0";
      else
        output << " 1";
    }
  }
  else
  {
    for (int i = 0; i < nq; i++)
    {
      if (results[i] > max_hits)
        iwdigits.append_number(output, max_hits);
      else
        iwdigits.append_number(output, results[i]);
    }
  }

  output.add('\n');

  return 1;
}

static int
tnass (data_source_and_type<Molecule> & input,
       IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess_molecule(*m);

    if (! tnass(*m))
      return 0;

    if (write_as_array || fingerprint_tag.length() || non_colliding_fingerprint_tag.length())
    {
      do_output(*m, output);
      MaybeFlush(output);
    }
  }

  output.flush();

  return 1;
}

static int
tnass_filter (iwstring_data_source & input, 
              IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with("$SMI<"))
      continue;

    buffer.remove_leading_chars(5);    // get rid if $SMI<
    buffer.chop();

    Molecule m;
    if (! m.build_from_smiles(buffer))
    {
      cerr << "tnass_filter: yipes cannot parse smiles '" << buffer << "'\n";
      cerr << " line " << input.lines_read() << endl;
      return 0;
    }

    (void) preprocess_molecule(m);

    (void) tnass(m);     // compute the results

    (void) write_fingerprint(output);
  }

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
tnass_filter (const char * fname, 
              IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << " tnass_filter::cannot open '" << fname << "'\n";
    return 0;
  }

  return tnass_filter(input, output);
}

static int
tnass (const char * fname, FileType input_type,
       IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tnass(input, output);
}

static int
handle_file_opening(Molecule_Output_Object & zstream,
                     Command_Line & cl,
                     char cflag)
{
  IWString fname;
  cl.value(cflag, fname);

  if (! cl.option_present('o'))
    zstream.add_output_type(FILE_TYPE_SMI);
  else if (! zstream.determine_output_types(cl))
  {
    cerr << "Cannot determine output types\n";
    return 0;
  }

  if (zstream.would_overwrite_input_files(cl, fname))
  {
    cerr << "Input and output (-" << cflag << ") must be distinct\n";
    return 19;
  }

  if (! zstream.new_stem(fname))
  {
    cerr << "Cannot open stream for -" << cflag << " molecules '" << fname << "'\n";
    return 21;
  }


  if (verbose)
    cerr << "Molecules for the -" << cflag << " option written to '" << fname << "'\n";

  return 1;
}

static void
create_header (IWString & header,
               const IWString & descriptor_stem = "")
{
  int nq = queries.number_elements();

  header.resize(nq * 20);

  for (int i = 0; i < nq; i++)
  {
    if (i > 0)
      header += ' ';

    if (descriptor_stem.length())
      header << descriptor_stem;

    IWString c = queries[i]->comment();

    if (special_processing_for_cnk)
      header << i;
    else if (0 == c.length())
    {
      if (0 == descriptor_stem.length())
        header << "qry";
      header << i;
    }
    else
    {
      c.gsub(' ', '_');
      header << c;
    }
  }

  return;
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush        flush output after each molecule\n";
  ::exit(0);
}

static int
tnass (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "N:vA:E:i:o:q:lm:n:g:t:W:D:ruG:bH:aJ:y:Y:Fx:kY:X:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(257);

  if (! process_elements(cl))
  {
    usage(2);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('t'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 't'))
    {
      cerr << "Cannot parse -t option\n";
      usage(11);
    }
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', max_hits) || max_hits < 1)
    {
      cerr << "The maximum number of hits per bit option (-x) must be a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Max hits per bit set to " << max_hits << endl;
  }

  if (cl.option_present('k'))
  {
    special_processing_for_cnk = 1;

    if (verbose)
      cerr << "Special processing for cnk descriptor\n";
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(33);
    }
  }

  if (cl.option_present('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      usage(33);
    }
  }

  if (0 == cl.option_count('q'))
  {
    cerr << "Must specify one or more queries via the -q option\n";
    usage(4);
  }

  if (cl.option_present('q'))
  {
    if (! queries.build_from_command_line(cl, 'q', verbose))
    {
      cerr << "Cannot parse -q option\n";
      return 32;
    }

    if (verbose > 1)
      queries.debug_print(cerr);
  }

  int nq = queries.number_elements();

  int nr = nq;
  if (default_fingerprint_nbits > nr)
    nr = default_fingerprint_nbits;

  results = new_int(nr);     // important to be zero initialised in case default_fingerprint_nbits > nq

  if (verbose)
  {
    hits_to_query = new_int(nq);
    hits_per_molecule = new_int(nq);
  }

  if (cl.option_present('b'))
  {
    queries.set_break_at_first_match(1);
    if (verbose)
      cerr << "Will break after first match\n";
  }

  if (cl.option_present('F') && ! cl.option_present('J'))
  {
    cerr << "The -F (work as filter) option requires the -J option\n";
    usage(32);
  }

  if (cl.option_present('J') && (cl.option_present('a') || cl.option_present('Y')))
  {
    cerr << "The -J (write as fingerprint) and -a (write as array) options are mutually exclusive\n";
    usage(27);
  }

  if (cl.option_present('y') && ! cl.option_present('J'))
  {
    cerr << "The -y option (number of bits to write) can only be used with the -J option (fingerprint tag)\n";
    usage(17);
  }

  IWString header;

  if (cl.option_present('J'))     // will write fingerprints
  {
    const_IWSubstring j = cl.string_value('J');

    if (j.starts_with("NC"))
    {
      non_colliding_fingerprint_tag = j;

      if (verbose)
        cerr << "Fingerprints written as non-colliding form, tag '" << non_colliding_fingerprint_tag << "'\n";

      if (! non_colliding_fingerprint_tag.ends_with('<'))
        non_colliding_fingerprint_tag += '<';
    }
    else
    {
      fingerprint_tag = j;

      if (! fingerprint_tag.ends_with('<'))
        fingerprint_tag << '<';

      if (verbose)
        cerr << "Fingerprints produced with the '" << fingerprint_tag << " dataitem\n";
    }

    if (cl.option_present('F'))
    {
      is_filter = 1;
      if (verbose)
        cerr << "Will work as a TDT filter\n";
    }

    if (cl.option_present('y'))
    {
      if (! cl.value('y', default_fingerprint_nbits) ||
            default_fingerprint_nbits < 1 ||
            default_fingerprint_nbits < queries.number_elements())
      {
        cerr << "The -y option must be followed by a whole number between 2 and " << queries.number_elements() << endl;
        usage(15);
      }
    }
    else
      default_fingerprint_nbits = queries.number_elements();

    if (0 != default_fingerprint_nbits % 8)
      default_fingerprint_nbits = (default_fingerprint_nbits / 8 + 1) * 8;

    if (verbose)
      cerr << "Query matches written as fingerprints of size " << default_fingerprint_nbits << endl;
  }

  if (cl.option_present('Y'))
  {
    if (! cl.option_present('a'))
    {
      cerr << "Header requested (-Y) but no array (-a) output\n";
      usage(5);
    }

    create_header(header, cl.string_value('Y'));
  }

  if (cl.option_present('a'))
  {
    write_as_array = cl.option_count('a');

    if (0 == header.length())
      create_header(header);
  }

  if (cl.option_present('n'))
  {
    if (! handle_file_opening(stream_for_non_matching_structures, cl, 'n'))
    {
      cerr << "Cannot process the -n option\n";
      return 32;
    }
  }

  if (cl.option_present('D') && ! cl.option_present('m'))
  {
    cerr << "The -D option (append query match details) only makes sense with the -m option\n";
    usage(32);
  }

  if (cl.option_present('m'))
  {
    if (! handle_file_opening(stream_for_matching_structures, cl, 'm'))
    {
      cerr << "Cannot process the -m option\n";
      return 33;
    }

    if (cl.option_present('D'))
    {
      append_query_details_to_name = 1;
      if (verbose)
        cerr << "Query match details appended to molecule names\n";
    }
  }

  if (cl.option_present('r'))
  {
    ignore_available_queries = 1;
    queries.set_ignore_precondition_matches_for_break(1);
    if (verbose)
      cerr << "Queries making features available ignored when determining matches\n";
  }

  if (cl.option_present('u'))
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_find_unique_embeddings_only(1);
    }

    if (verbose)
      cerr << "Will only find unique embeddings\n";
  }

  if (cl.option_present('G'))
  {
    const char * g = cl.option_value('G');

    if (! bob_coner_stream.open(g))
    {
      cerr << "Cannot open matched atom stream '" << g << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Matched atom lists written to '" << g << "'\n";
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (is_filter)     // working as a filter, no need to check input type
    ;
  else if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (header.length() && write_as_array)
    output << "Name " << header << '\n';

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (is_filter)
    {
      if (! tnass_filter(cl[i], output))
      {
        return i + 1;
      }
    }
    else if (! tnass(cl[i], input_type, output))
    {
      return i + 1;
    }
  }

  if (0 == verbose && ! cl.option_present('m') && ! cl.option_present('n') && ! is_filter)
    cerr << "Read " << molecules_read << " molecules, " << molecules_matching << " matched one or more queries\n";

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, " << molecules_matching << " matched one or more queries\n";
    for (int i = 0; i < nq; i++)
    {
      cerr << hits_to_query[i] << " molecules hit query " << i << " (" << queries[i]->comment() << ")\n";
    }

    cerr << "Over " << queries.searches_performed() << " optimisation " << queries.optimisation_level() << endl;

    Accumulator_Int<int> acc;
    for (int i = 0; i < nq; i++)
    {
      if (hits_per_molecule[i])
      {
        cerr << hits_per_molecule[i] << " molecules hit " << i << " queries\n";
        acc.extra(i, hits_per_molecule[i]);
      }
    }

    if (acc.n() > 1)
      cerr << "Average hits per molecule " << acc.average() << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tnass(argc, argv);

  return rc;
}
