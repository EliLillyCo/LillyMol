#include <assert.h>
#include <math.h>
#include <errno.h>

#include <iostream>
#include <memory>

#if (__GNUC_MINOR__ == 95)
#include <strstream>
#else
#include <sstream>
#endif

#include <google/protobuf/message.h>
#include <google/protobuf/text_format.h>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/demerit.pb.h"
#include "substructure_demerits.h"
#include "demerit.h"

using std::cerr;
using std::endl;

static int verbose = 0;     // visible in substructure_demerits

static int keep_going_after_rejection = 0;

static int molecules_read = 0;

const char * prog_name = nullptr;

static int molecules_receiving_demerits = 0;
static int molecules_rejected = 0;

static Elements_to_Remove elements_to_remove;

static Molecule_Output_Object stream_for_non_rejected_molecules;
static Molecule_Output_Object stream_for_rejected_molecules;

// Write all molecules as protos to a single output stream.
static IWString_and_File_Descriptor proto_stream;

static IWString stem_for_demerit_file;

static int do_hard_coded_substructure_queries = 1;

static int append_demerit_text_to_name = 0;

static int write_rejection_reason_like_tsubstructure = 0;

static int write_rejection_reasons_separately = 0;

static int write_sorted_vbar_list = 0;

static int skip_molecules_with_abnormal_valences;
static int molecules_with_abnormal_valences = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*
  For molecules that lie between the hard and low atom count cutoffs, we apply demerits
*/

static int soft_lower_atom_count_cutoff = 0;
static int hard_lower_atom_count_cutoff = 0;
static int lower_atom_count_demerit = 100;

static int soft_upper_atom_count_cutoff = 0;
static int hard_upper_atom_count_cutoff = 0;
static int upper_atom_count_demerit = 100;

/*
  Apr 2017. For PAINS processing
*/

static int demerit_numeric_value_index = 0;

static int make_implicit_hydrogens_explicit = 0;

static void
do_atom_count_demerits(Molecule & m,
                       Demerit & demerit)
{
  int nf = m.number_fragments();

  int matoms;
  if (1 == nf)
    matoms = m.natoms();
  else
  {
    matoms = 0;
    for (int i = 0; i < nf; i++)
    {
      int f = m.atoms_in_fragment(i);
      if (f > matoms)
        matoms = f;
    }
  }

  if (hard_lower_atom_count_cutoff > 0 && matoms <= hard_lower_atom_count_cutoff)
  {
    demerit.extra(lower_atom_count_demerit, "too_few_atoms");

    return;
  }

// Feb 2005. Heuristic for adding demerits beyond the hard atom count cutoff

  if (hard_upper_atom_count_cutoff > 0 && matoms >= hard_upper_atom_count_cutoff)
  {
    demerit.extra(rejection_threshold() + 6 * (matoms - hard_upper_atom_count_cutoff), "too_many_atoms");

//  cerr << hard_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff, matoms = " << matoms << endl;
    return;
  }

  if (matoms > hard_lower_atom_count_cutoff && matoms < soft_lower_atom_count_cutoff)
  {
    float r = static_cast<float>(soft_lower_atom_count_cutoff - matoms) / static_cast<float>(soft_lower_atom_count_cutoff - hard_lower_atom_count_cutoff);
    int d = static_cast<int>(lower_atom_count_demerit * r);

    if (0 == d)
      d = 1;

    demerit.extra(d, "too_few_atoms");

    return;
  }

  if (matoms > soft_upper_atom_count_cutoff && matoms < hard_upper_atom_count_cutoff)
  {
    float r = static_cast<float>(matoms - soft_upper_atom_count_cutoff) / static_cast<float>(hard_upper_atom_count_cutoff - soft_upper_atom_count_cutoff);
    int d = static_cast<int>(upper_atom_count_demerit * r);
//  cerr << "too_many_atoms: " << soft_upper_atom_count_cutoff << " soft_upper_atom_count_cutoff, matoms = " << matoms << ", d = " << d << endl;
//  cerr << "upper_atom_count_demerit " << upper_atom_count_demerit << endl;
//  cerr << "too_many_atoms: soft_upper_atom_count_cutoff " << soft_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff " << hard_upper_atom_count_cutoff << " matoms " << matoms << " D = " << d << endl;
    if (0 == d)
      d = 1;

    demerit.extra(d, "too_many_atoms");

    return;
  }

  return;
}


static void
run_a_set_of_queries(Molecule_to_Match & target,
                     Demerit & demerit,
                     resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  int nqueries = queries.number_elements();

  for (int i = 0; i < nqueries; i++)
  {
    Substructure_Hit_Statistics * q = queries[i];

    int nhits = q->substructure_search(target);
    if (0 == nhits)
      continue;

    if (verbose > 1)
      cerr << nhits << " matches to query " << i << ' ' << q->comment() << endl;

    double d;
    (void) q->numeric_value(d, demerit_numeric_value_index);
    
    d = d * nhits;

    int intd = static_cast<int>(d * 1.0001);

    if (intd >= rejection_threshold())
    {
      demerit.reject(q->comment());
      if (0 == keep_going_after_rejection)
        return;
    }
    else
      demerit.extra(intd, q->comment());

    if (demerit.rejected() && 0 == keep_going_after_rejection)
      return;
  }

  return;
}

static void
iwdemerit(Molecule & m,
          resizable_array_p<Substructure_Hit_Statistics> & q1,
          resizable_array_p<Substructure_Hit_Statistics> & q2,
          Demerit & demerit)
{
  if (soft_lower_atom_count_cutoff > 0 || soft_upper_atom_count_cutoff > 0)
  {
    do_atom_count_demerits(m, demerit);
    if (demerit.rejected() && 0 == keep_going_after_rejection)
      return;
  }

  if (skip_molecules_with_abnormal_valences && ! m.valence_ok())
  {
    demerit.reject("valence");
    molecules_with_abnormal_valences++;
    if (0 == keep_going_after_rejection)
      return;
  }

  if (do_hard_coded_substructure_queries)
  {
    substructure_demerits::hard_coded_queries(m, demerit);
//  cerr << "After hard coded queries, score is " << demerit.score() << endl;

    if (demerit.rejected() && 0 == keep_going_after_rejection)
      return;
  }

  if (q1.number_elements())
  {
    Molecule_to_Match target(&m);

    run_a_set_of_queries(target, demerit, q1);

//  cerr << "After command line queries, score is " << demerit.score() << " rej? " << demerit.rejected() << endl;
    if (demerit.rejected())
      return;
  }

  if (q2.number_elements())
  {
    m.reduce_to_largest_fragment_carefully();
    Molecule_to_Match target(&m);

    run_a_set_of_queries(target, demerit, q2);

    if (demerit.rejected())
      return;
  }

  return;
}

/*
  In jul 97, people wanted to know how many molecules received multiple
  demerits, and how many were rejected as a result of multiple demerits
*/

static extending_resizable_array<int> demerits_per_molecule;
static extending_resizable_array<int> demerits_per_rejected_molecule;

/*
  Later in Jul 97 they wanted a list of molecules for which there were
  rejected as a result of multiple different types of demerit
*/

static Molecule_Output_Object stream_for_multiple_demerits;

#ifdef NOT_USED_ANY_MOPRE_ASDASDASD
static int
do_append_demerit_text_to_name(Molecule & m,
                               const Demerit & demerit,
                               const int score,
                               const IWString & reason)
{
//cerr << "do_append_demerit_text_to_name " << score << ' ' << reason << endl;

  IWString new_name(m.molecule_name());

  if (write_rejection_reasons_separately && write_rejection_reason_like_tsubstructure)
    new_name << ' ' << reason;
  else if (write_rejection_reason_like_tsubstructure)
    new_name << " (" << demerit.number_different_demerits_applied() << " matches to 'D" << score << ' ' << reason << "')";
  else
    new_name << " : D(" << score << ") " << reason;

  m.set_name(new_name);

//cerr << "new_name " << new_name << endl;

  return 1;
}
#endif

static int
do_append_demerit_text_to_name(Molecule & m,
                               Demerit & demerit)
{
  demerit.do_sort();

  IWString new_name(m.name());

  const auto & demerits = demerit.demerits();

  const int n = demerits.number_elements();

  const int score = demerit.score();

  if (write_sorted_vbar_list)
  {
    new_name << " : D(" << score << ") " << score;
    for (int i = 0; i < n; ++i)
    {
      new_name << ":" << demerits[i]->demerit() << '|' << demerits[i]->reason();
    }
  }
  else if (write_rejection_reasons_separately)
  {
    new_name << " : D(" << score << ")";

    for (int i = 0; i < n; ++i)
    {
      new_name << ' ' << demerits[i]->reason();
    }
  }
  else if (write_rejection_reason_like_tsubstructure)
  {
    new_name << " : D(" << score << ")";
    for (int i = 0; i < n; ++i)
    {
      new_name << " (1 matches to 'D" << demerits[i]->demerit() << ' ' << demerits[i]->reason() << "')";
    }
  }
  else
  {
    new_name << " : D(" << score << ") ";
    for (int i = 0; i < n; ++i)
    {
      if (i > 0)
        new_name << ':';
      new_name << demerits[i]->reason();
    }
  }

  m.set_name(new_name);

  return 1;
}

// Write results in `demerit` to `output`.
// `mname` will be set if m.name() has been over-written
static int
WriteAsProto(Molecule& m,
             const IWString& mname,
             const Demerit& demerit,
             IWString_and_File_Descriptor& output) {
  MedchemRules::Molecule proto;
  proto.set_smiles(m.smiles().AsString());

  // If the molecule name was NOT saved, then use the regular name.
  if (mname.empty()) {
    proto.set_name(m.name().AsString());
  } else {
    proto.set_name(mname.AsString());
  }

  if (demerit.rejected()) {
    proto.set_rejected(true);
  }

  const resizable_array_p<Demerit_and_Reason> & demerits = demerit.demerits();
  for (const Demerit_and_Reason* dar : demerits) {
    MedchemRules::QueryMatch * q = proto.add_query_match();
    q->set_name(dar->reason().AsString());
    q->set_demerit(dar->demerit());
  }

  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);
  std::string as_string;
  printer.PrintToString(proto, &as_string);
  output.write(as_string.data(), as_string.size());
  output << '\n';
  output.write_if_buffer_holds_more_than(8192);
  return 1;
}

static int
iwdemerit(Molecule & m,
          resizable_array_p<Substructure_Hit_Statistics> & q1,
          resizable_array_p<Substructure_Hit_Statistics> & q2,
          std::ofstream & output)
{
  elements_to_remove.process(m);

  IWString save_name;
  if (proto_stream.active() && append_demerit_text_to_name) {
    save_name = m.name();
  }

  Demerit demerit;

  iwdemerit(m, q1, q2, demerit);

  if (demerit.score())
  {
    molecules_receiving_demerits++;
    demerits_per_molecule[demerit.number_different_demerits_applied()]++;

    if (demerit.rejected())
    {
      molecules_rejected++;
      demerits_per_rejected_molecule[demerit.number_different_demerits_applied()]++;
    }

    if (append_demerit_text_to_name)
      do_append_demerit_text_to_name(m, demerit);
  }

//cerr << m.name() << " rejected? " << demerit.rejected() << " stream is " << stream_for_non_rejected_molecules.active() << endl;

  if (proto_stream.active()) {
    WriteAsProto(m, save_name, demerit, proto_stream);
    // Should we return now?
  }

  if (demerit.rejected())
  {
    if (stream_for_rejected_molecules.active())
      stream_for_rejected_molecules.write(&m);
    if (stream_for_multiple_demerits.active() && ! demerit.rejected_by_single_rule())
      stream_for_multiple_demerits.write(&m);
  }
  else if (stream_for_non_rejected_molecules.active())
    stream_for_non_rejected_molecules.write(&m);

  if (verbose > 1 && demerit.score())
    demerit.debug_print(cerr);

  if (! output.rdbuf()->is_open())
    return 1;

  output << "$SMI<" << m.smiles() << ">\n";
  output << "PCN<" << m.name() << ">\n";

  if (demerit.score())
    demerit.write_in_tdt_form(output);

  output << "|\n";
  output.flush();

  return output.good();
}

static void
preprocess(Molecule & m)
{
  m.remove_all(1);

  if (make_implicit_hydrogens_explicit)
    m.make_implicit_hydrogens_explicit();

  return;
}

static int
iwdemerit(data_source_and_type<Molecule> & input,
          resizable_array_p<Substructure_Hit_Statistics> & q1,
          resizable_array_p<Substructure_Hit_Statistics> & q2,
          std::ofstream & output)
{
  assert (input.good());

  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    if (verbose > 1)
      cerr << input.molecules_read() <<  " processing '" << m->name() << "'\n";

    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! iwdemerit(*m, q1, q2, output))
      return 0;
  }

  if (verbose)
  {
    cerr << input.molecules_read() << " molecules read\n";
    elements_to_remove.report(cerr);
    cerr << molecules_receiving_demerits << " molecules were assigned demerits\n";
    cerr << molecules_rejected << " molecules were rejected\n";

    cerr << "Details on queries\n";

    if (do_hard_coded_substructure_queries)
      substructure_demerits::hard_coded_queries_statistics(cerr);

    for (int i = 0; i < queries.number_elements(); i++)
    {
      cerr << "Results of query " << i << endl;
      queries[i]->report(cerr, verbose);
    }
  }

  return 1;
}

static int
iwdemerit(const char * fname, FileType input_type,
          resizable_array_p<Substructure_Hit_Statistics> & q1,
          resizable_array_p<Substructure_Hit_Statistics> & q2)
{
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  std::ofstream output;
  
  if (stem_for_demerit_file.length())
  {
    IWString output_file_name = stem_for_demerit_file;

    change_suffix(output_file_name, "demerit");

    output.open(output_file_name.null_terminated_chars());
    if (! output.good())
    {
      cerr << prog_name << ": cannot open '" << output_file_name << "' for output\n";
      return 2;
    }
  }

  return iwdemerit(input, q1, q2, output);
}

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
  cerr << "Usage : " << prog_name << " options file1 file2 file3 ....\n";
  cerr << "  -q <file>      specify substructure query file\n";
  cerr << "  -R <stem>      write rejected structures to 'stem.otype'\n";
  cerr << "  -G <stem>      write non rejected (good) structures to 'stem.otype'\n";
  cerr << "  -S <stem>      specify stem for .demerit file\n";
  cerr << "  -M <stem>      write molecules rejected by cumulative demerits to <stem>\n";
  cerr << "  -P <fname>     write results in proto form to <fname>\n";
  cerr << "  -k             check all criteria even if a molecule is already rejected\n";
  cerr << "  -c <number>    atom count cutoffs - enter -c help for info\n";
  cerr << "  -O hard        skip all the hard coded substructure queries\n";
  cerr << "  -V             reject molecules containing unusual valences\n";
  cerr << "  -t             append demerit text to molecule names\n";
  cerr << "  -f <n>         molecules are rejected when they have <n> demerits\n";
  cerr << "  -x             atom count demerits remain scaled to 100\n";
//cerr << "  -I <nrings>    set threshold for too many rings rejection\n";  deprecated
  cerr << "  -Z <rsize>     set threshold for C7 ring size (default 7)\n";
  cerr << "  -z <length>    set threshold for long carbon chain (default 7)\n";
  cerr << "  -C <file>      control file for what demerits to apply\n";
  cerr << "  -r             only do rejection rules - no demerits\n";
  cerr << "  -u             write rejection reasons like tsubstructure\n";
  cerr << "  -y             collate individual demerits into name\n";
  cerr << "  -d <number>    set all demerit numeric values to <number>\n";
  cerr << "  -N ...         charge assigner specifications, enter '-N help' for info\n";
  cerr << "  -W ...         miscellaneous options, enter -W for help\n";
  cerr << "  -E <symbol>    create element with symbol\n";
  cerr << "  -X <symbol>    before processing, delete all <symbol> atoms\n";
  cerr << "  -o <type>      file type for structures written\n";
  cerr << "  -i <type>      specify input file type\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static void
separate_depending_on_fragment_match(resizable_array_p<Substructure_Hit_Statistics> & queries,
                                     resizable_array_p<Substructure_Hit_Statistics> & q1,
                                     resizable_array_p<Substructure_Hit_Statistics> & q2)
{
  int n = queries.number_elements();

  for (int i = n - 1; i >= 0; i--)
  {
    Substructure_Hit_Statistics * q = queries.remove_no_delete(i);

    if (q->only_keep_matches_in_largest_fragment())
      q2.add(q);
    else
      q1.add(q);
  }

  if (verbose > 1)
    cerr << "Separated queries into " << q1.number_elements() << " and " << q2.number_elements() << " queries\n";

  return;
}

static void
display_dash_W_options(std::ostream & os)
{
  os << " -W dnv=<ndx>         which numeric value in the query file to use as the demerit score\n";
  os << " -W maxe=<n>          maximum number of substructure queries to identify\n";
  os << " -W imp2exp           make implicit Hydrogens explicit\n";
  os << " -W nokekule          aromatic bonds no longer remember their Kekule forms\n";
  os << " -W slist             write a sorted list of demerit values and reasons\n";

  exit(1);
}

static void
display_atom_cutoff_options(std::ostream & output)
{
  output << " -c hmin=<n>           hard min cutoff - molecules with fewer atoms rejected\n";
  output << " -c smin=<h>           soft min cutoff - molecules with fewer atoms demerited\n";
  output << " -c smax=<h>           soft max cutoff - molecules with more  atoms demerited\n";
  output << " -c hmax=<h>           hard max cutoff - molecules with more  atoms rejected\n";
}

int
SetupProtoStream(IWString& fname,
                IWString_and_File_Descriptor& output)  {
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "SetupProtoStream:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
SetupOutputFile(Command_Line& cl,
                const char flag,
                Molecule_Output_Object& output ) {
  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl)) {
    cerr << "Cannot discern file types for -R output\n";
    return 0;
  }

  IWString fname = cl.string_value(flag);
  if (! output.new_stem(fname)) {
    cerr << "Cannot set -R output name to '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Rejected structures written to '" << fname << "'\n";
  }

  return 1;
}

int
iwdemerit(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "M:VX:tA:S:R:G:O:kd:Dq:E:vi:o:c:C:N:uyf:xrlI:Z:z:W:P:");

  if (cl.unrecognised_options_encountered())
    usage(1);

  verbose = cl.option_count('v');

  if (verbose)
    substructure_demerits::set_verbose(verbose);

  if (cl.option_present('C'))
  {
    IWString fname = cl.string_value('C');

    if (! substructure_demerits::initialise_hard_coded_queries_to_do(fname))
    {
      cerr << "Cannot initialise hard coded queries subset (-C option) from '" << fname << "'\n";
      return 8;
    }
  }

  if (cl.option_present('r'))
  {
    substructure_demerits::set_only_apply_rejection_rules();
    if (verbose)
      cerr << "Only rejection rules will be applied\n";
  }

  if (cl.option_present('N'))
  {
    Charge_Assigner & charge_assigner = substructure_demerits::charge_assigner();

    if (! charge_assigner.construct_from_command_line(cl, verbose > 0, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      return 3;
    }
  }

  if (cl.option_present('f'))
  {
    int f;
    if (! cl.value('f', f) || f < 1)
    {
      cerr << "The rejection threshold (-f) must be a whole +ve number\n";
      usage(4);
    }

    set_rejection_threshold(f);

    if (verbose)
      cerr << "Rejection threshold set to " << f << endl;

    if (! cl.option_present('x'))
    {
      lower_atom_count_demerit = f;
      upper_atom_count_demerit = f;
    }
  }

  if (cl.option_present('I'))
  {
    int r;
    if (! cl.value('I', r) || r < 1)
    {
      cerr << "The too many rings rejection threshold (-I) must be a whole +ve number\n";
      usage(1);
    }

    substructure_demerits::set_substructure_demerits_too_many_rings(r);

    if (verbose)
      cerr << "Will reject molecules having more than " << r << " rings\n";
  }

  if (cl.option_present('Z'))
  {
    int r;
    if (! cl.value('Z', r) || r < 3)
    {
      cerr << "The size of C7 ring rejection threshold (-Z) must be a valid ring size\n";
      usage(1);
    }

    substructure_demerits::set_substructure_demerits_ring_size_too_large(r);

    if (verbose)
      cerr << "Will reject molecules having a mostly saturated Carbon ring containing " << r << " or more atoms\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (! cl.option_present('i'))
  {
    if (! all_files_recognised_by_suffix(cl))
    {
      cerr << "Cannot discern input type(s)\n";
      return 3;
    }
  }
  else if (! process_input_type(cl, input_type))
  {
    cerr << prog_name << ": cannot discern input type\n";
    usage(2);
  }
  
  if (cl.option_present('R'))
  {
    if (! SetupOutputFile(cl, 'R', stream_for_rejected_molecules)) {
      cerr << "Cannot setup file for rejected molecules (-R)\n";
      return 1;
    }
  }

  if (cl.option_present('G'))
  {
    if (! SetupOutputFile(cl, 'G', stream_for_non_rejected_molecules)) {
      cerr << "Cannot setup file for ok molecules (-G)\n";
      return 1;
    }
  }

  if (cl.option_present('P')) {
    IWString fname = cl.option_value('P');
    if (! SetupProtoStream(fname, proto_stream)) {
      cerr << "Cannot initialise proto output '" << fname << "\n";
      return 1;
    }
    if (verbose) {
      cerr << "Output written in proto form to '" << fname << "'\n";
    }
  }

  if (cl.option_present('S'))
  {
    cl.value('S', stem_for_demerit_file);

    if (verbose)
      cerr << "Demerit file will have name stem '" << stem_for_demerit_file << "'\n";
  }

  (void) process_elements(cl);

  if (cl.option_present('M'))
  {
    const_IWSubstring fname;
    cl.value('M', fname);

    if (! stream_for_multiple_demerits.determine_output_types(cl))
    {
      cerr << "Cannot discern output types for -M output\n";
      return 76;
    }

    if (stream_for_multiple_demerits.would_overwrite_input_files(cl, fname))
    {
      cerr << "Cannot overwrite input file(s) '" << fname << "' with -M file\n";
      return 1;
    }

    if (! stream_for_multiple_demerits.new_stem(fname))
    {
      cerr << "Cannot set -M output stem to '" << fname << "'\n";
      return 18;
    }

    if (verbose)
      cerr << "Molecules rejected for multiple purposes written to '" << fname << "'\n";
  }

  if (cl.option_present('X'))
  {
    if (! elements_to_remove.construct_from_command_line(cl, verbose, 'X'))
    {
      cerr << "Cannot discern elements to be removed\n";
      usage(21);
    }
  }

  if (cl.option_present('t'))
  {
    append_demerit_text_to_name = 1;
    if (verbose)
      cerr << "Demerit types will be appended to molecule names\n";
  }

  if (cl.option_present('u'))
  {
    append_demerit_text_to_name = 1;
    write_rejection_reason_like_tsubstructure = 1;
    set_store_demerit_reasons_like_tsubstructure(1);
    if (verbose)
      cerr << "Rejection reasons written like tsubstructure\n";
  }

  if (cl.option_present('y'))
  {
    set_demerit_reason_contains_individual_demerits(1);
    write_rejection_reasons_separately = 1;
    if (verbose)
      cerr << "Demerit reasons separated\n";
  }

  if (cl.option_present('c'))
  {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('c', c, i++))
    {
      if ("help" == c)
      {
        display_atom_cutoff_options(cerr);
        return 0;
      }

      const_IWSubstring directive;
      int dvalue;
      if (! c.split_into_directive_and_value(directive, '=', dvalue))
      {
        cerr << "Invalid -c directive '" << c << "'\n";
        usage(5);
      }

      if (dvalue < 0)
      {
        cerr << "INvalid numeric directive, cannot be negative '" << c << "'\n";
        usage(6);
      }

      if ("hmin" == directive)
      {
        hard_lower_atom_count_cutoff = dvalue;
      }
      else if ("smin" == directive)
      {
        soft_lower_atom_count_cutoff = dvalue;
      }
      else if ("smax" == directive)
      {
        soft_upper_atom_count_cutoff = dvalue;
      }
      else if ("hmax" == directive)
      {
        hard_upper_atom_count_cutoff = dvalue;
      }
      else if ("mindmrt" == directive)
      {
      }
      else if ("maxdmrt" == directive)
      {
      }
      else
      {
        cerr << "Unrecognised -c qualifier '" << c << "'\n";
        usage(5);
      }
    }

    if (0 == soft_lower_atom_count_cutoff && 0 == hard_lower_atom_count_cutoff)
      ;
    else if (hard_lower_atom_count_cutoff <= soft_lower_atom_count_cutoff)
      ;
    else
    {
      cerr << "Invalid lower atom count cutoff values soft:" << soft_lower_atom_count_cutoff << " hard:" << hard_lower_atom_count_cutoff << endl;
      return 5;
    }

    if (0 == soft_upper_atom_count_cutoff && 0 == hard_upper_atom_count_cutoff)
      ;
    else if (soft_upper_atom_count_cutoff < hard_upper_atom_count_cutoff)
      ;
    else
    {
      cerr << "Invalid upper atom count cutoff values soft:" << soft_upper_atom_count_cutoff << " hard:" << hard_upper_atom_count_cutoff << endl;
      return 5;
    }

    if (0 == soft_lower_atom_count_cutoff && 0 == soft_upper_atom_count_cutoff)
      ;
    else if (soft_lower_atom_count_cutoff < soft_upper_atom_count_cutoff)
      ;
    else
    {
      cerr << "Soft cutoffs invalid, lower:" << soft_lower_atom_count_cutoff << " upper:" << soft_upper_atom_count_cutoff << endl;
      return 8;
    }

//  Should do more checks here...
  }

  if (! process_standard_smiles_options(cl, verbose))
  {
    usage(7);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    usage(8);
  }

  if (cl.option_present('k'))
  {
    keep_going_after_rejection = 1;
    if (verbose)
      cerr << "Will scan all criteria, ignoring previous rejections\n";

    substructure_demerits::set_keep_going_after_rejection(1);
  }

  if (cl.option_present('O'))
  {
    IWString tmp;
    int i = 0;
    while (cl.value('O', tmp, i++))
    {
      if ("hard" == tmp)
      {
        do_hard_coded_substructure_queries = 0;
        if (verbose)
          cerr << "Hard coded substructure queries will be skipped\n";
      }
      else
      {
        cerr << "Unrecognised skip specifier '" << tmp << "'\n";
        usage(19);
      }
    }
  }

  int max_embeddings = 0;

  if (cl.option_present('W'))
  {
    const_IWSubstring w;
    for (int i = 0; cl.value('W', w, i); ++i)
    {
      if (w.starts_with("dnv="))
      {
        w.remove_leading_chars(4);
        if (! w.numeric_value(demerit_numeric_value_index) || demerit_numeric_value_index < 0)
        {
          cerr << "The demerit numeric value index must be a valid index '" << w << "' invalid\n";
          return 1;
        }

        if (verbose)
          cerr << "Demerit numeric value index set to " << demerit_numeric_value_index << endl;
      }
      else if (w.starts_with("maxe="))
      {

        w.remove_leading_chars(5);
        if (! w.numeric_value(max_embeddings) || max_embeddings < 1)
        {
          cerr << "The maximum number of query embeddings 'maxe=' must be a whole +ve number\n";
          return 1;
        }

        if (verbose)
          cerr << "Queries limited to " << max_embeddings << " embeddings\n";
      }
      else if ("imp2exp" == w)
      {
        make_implicit_hydrogens_explicit = 1;

        if (verbose)
          cerr << "Implicit hydrogens will be made explicit\n";
      }
      else if ("nokekule" == w)
      {
        set_aromatic_bonds_lose_kekule_identity(1);

        if (verbose)
          cerr << "Aromatic bonds will no longer match Kekule forms\n";
      }
      else if ("slist" == w)
      {
        write_sorted_vbar_list = 1;
        append_demerit_text_to_name = 1;
        if (verbose)
          cerr << "Will write a sorted list of demerits and reasons\n";
      }
      else if ("help" == w)
      {
        display_dash_W_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -W qualifier '" << w << "'\n";
        display_dash_W_options(cerr);
      }
    }
  }

  int output_formats = 0;
  if (write_sorted_vbar_list)
    output_formats++;
  if (write_rejection_reasons_separately)
    output_formats++;
  if (write_rejection_reason_like_tsubstructure)
    output_formats++;

  if (output_formats > 1)
  {
    cerr << "Multiple output formats specified, cannot continue\n";
    usage(1);
  }

  if (cl.option_present('V'))
  {
    skip_molecules_with_abnormal_valences = 1;
    if (verbose)
      cerr << "Molecules containing unusual valences will be rejected\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << prog_name << ": no files specified\n";
    usage(2);
  }

  int all_demerits_same_numeric_value = -1;

  if (cl.option_present('d'))
  {
    if (! cl.value('d', all_demerits_same_numeric_value) || all_demerits_same_numeric_value < 0 || all_demerits_same_numeric_value > 100)
    {
      cerr << "The set all numeric demerits option (-d) must be a whole number between 0 and 100\n";
      usage(1);
    }

    if (verbose)
      cerr << "All demerits assigned numeric value " << all_demerits_same_numeric_value << endl;

    substructure_demerits::set_all_numeric_demerit_values(all_demerits_same_numeric_value);
  }

  if (cl.option_present('z'))
  {
    int t;
    if (! cl.value('z', t) || t < 2)
    {
      cerr << "THe length of Carbon chain to be rejected (-z) must be a whole +ve number > 1\n";
      usage(1);
    }

    substructure_demerits::set_cx_chain_rejection_length(t);
  }

// As we read in queries, we make sure that each component has a numeric
// value - saves checking them later...

  if (cl.option_present('r'))
    ;
  else if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose, 'q'))
    {
      cerr << prog_name << ": processing of -q option failed\n";
      return 8;
    }

    for (int i = 0; i < queries.number_elements(); i++)
    {
      Substructure_Hit_Statistics * q = queries[i];

      q->set_find_unique_embeddings_only(1);

      for (int j = 0; j < q->number_elements(); j++)
      {
        Single_Substructure_Query * sq = q->item(j);

        double d;
        if (! sq->numeric_value(d, demerit_numeric_value_index))
        {
          cerr << "Yipes, query '" << q->comment() << " has no demerit value\n";
          return 13;
        }
  
        if (d < 0.0)
        {
          cerr << "Hmmm, query '" << q->comment() << "' has a negative demerit " << d << endl;
          return 14;
        }

        if (all_demerits_same_numeric_value >= 0)
          sq->set_numeric_value(static_cast<double>(all_demerits_same_numeric_value), 0);
      }
    }
  }

  set_remove_hits_not_in_largest_fragment_behaviour(1);    // to get reproducible behaviour with multiple instances of largest fragment

  if (max_embeddings > 0)
  {
    for (int i = 0; i < queries.number_elements(); ++i)
    {
      queries[i]->set_max_matches_to_find(max_embeddings);
    }
  }

  int rc = 0;

  if (cl.option_present('l'))
  {
    resizable_array_p<Substructure_Hit_Statistics> q1, q2;
    separate_depending_on_fragment_match(queries, q1, q2);

    for (int i = 0; i < cl.number_elements(); i++)
    {
      const char *fname = cl[i];

      if (verbose)
        cerr << prog_name << " processing '" << fname << "'\n";

      if (! iwdemerit(fname, input_type, q1, q2))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    resizable_array_p<Substructure_Hit_Statistics> notused;

    for (int i = 0; i < cl.number_elements(); i++)     // each argument is a file
    {
      const char *fname = cl[i];

      if (verbose)
        cerr << prog_name << " processing '" << fname << "'\n";

      if (! iwdemerit(fname, input_type, queries, notused))
      {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose)
  {
    if (cl.number_elements() > 1)
      cerr << molecules_read << " molecules read\n";

    if (molecules_with_abnormal_valences)
      cerr << "Skipped " << molecules_with_abnormal_valences << " molecules with abnormal valences\n";

    for (int i = 0; i < demerits_per_molecule.number_elements(); i++)
    {
      int j = demerits_per_molecule[i];
      if (j)
        cerr << j << " molecules had " << i << " sources of demerits\n";
    }

    for (int i = 0; i < demerits_per_rejected_molecule.number_elements(); i++)
    {
      int j = demerits_per_rejected_molecule[i];
      if (j)
        cerr << j << " rejected molecules had " << i << " sources of demerits\n";
    }
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = iwdemerit(argc, argv);

  return rc;
}
