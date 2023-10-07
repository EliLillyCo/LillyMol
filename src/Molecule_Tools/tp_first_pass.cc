/*
  First phase for 3rd party processing.

  Can do the following:   (although driven by command line options)

    Remove any molecule in which no fragment contains more than 11 atoms.
    Remove any molecule with covalently bound non-organic atoms
    Remove any molecule containing atoms types we don't want
    Remove any molecule which does not contain C and (O or N)
    Remove hydrogens.
    Apply all chemical standardisations
*/

#include <cassert>
#include <iostream>
#include <limits>
#include <memory>

#include <google/protobuf/message.h>
#include <google/protobuf/text_format.h>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "tp_first_pass_lib.h"
#include "Molecule_Tools/demerit.pb.h"

const char *prog_name;

namespace tp_first_pass {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "The following options are recognised\n";
  cerr << "  -c <number>    exclude molecules with atom count below <number>\n";
  cerr << "  -C <number>    exclude molecules with atom count above <number>\n";
  cerr << "  -r <number>    omit molecules with fewer than <number> rings\n";
  cerr << "  -R <number>    omit molecules with more than <number> rings\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -Z <ringsize>  upper ring size limit\n";
  cerr << "  -I <0,1>       <exclude,include> molecules containing isotopes\n";
  cerr << "  -I change      change all isotopic atoms to their natural form\n";
  cerr << "  -s ...         chirality related options, enter '-s help' for info\n";
//cerr << "  -n <number>    assign sequential numbers R(%d) starting with <number>\n";  deprecated
  cerr << "  -b <ratio>     skip molecules with ring bond ratio's >= than <ratio>\n";
  cerr << "  -w             run all checks - normally discards molecules once problem found\n";
  cerr << "  -V             skip any molecule with abnormal valences\n";
  cerr << "  -k             allow molecules having no \"interesting\" atoms to pass\n";
  cerr << "  -f <fraction>  minimum fraction of interesting atoms required\n";
  cerr << "  -y             allow non periodic table elements if they are not connected\n";
  cerr << "  -F zz:nn       discard if more than nn instances of atomic number zz\n";
  cerr << "  -F zz:fraction discard if more than fraction instances of atomic number zz\n";
  cerr << "  -H <size>      discard if a ring system with > size rings\n";
  cerr << "  -L <fname>     write rejected molecules to <fname>\n";
  cerr << "  -L appreason   append rejection reason to name in reject file\n";
  cerr << "  -L tsub        write the rejection reason like tsubstructure\n";
  cerr << "  -L proto       write rejected molecules as Demerit protos\n";
  cerr << "  -K <text>      prepend 'text' before rejection reason\n";
//cerr << "  -a             append rejection reason to name in reject log\n";   deprecated.
//cerr << "  -u             write the rejection reason like tsubstructure\n";   deprecated.
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -P <text>      append <text> to the name of every molecule\n";
  cerr << "  -p remove      remove non printing characters in molecule names\n";
  cerr << "  -p <char>      convert non printing characters in molecule names to '<char>'\n";
  cerr << "  -x <natoms>    reject if molecule contains > 1 fragment with at least <natoms> atoms\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol> (use -E '*' for auto create)\n";
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  (void) display_standard_etrans_options (cerr, 't');
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -v             verbose output\n";

  exit(rc);
}

struct Parameters {
  int verbose = 0;

  // We may not want strange characters in molecule names.
  char translate_non_printing_chars = '\0';
  int remove_non_printing_chars = 0;

  int molecules_read = 0;
  int molecules_written = 0;

  Chemical_Standardisation chemical_standardisation;

  int run_all_checks = 0;

  Accumulator_Int<int> natoms_accumulator;

  int append_rejection_reason_to_name = 0;

  int write_rejection_reason_like_tsubstructure = 0;

  // We can make the output easier to parse if we have a fixed string before 
  // the rejection reason

  IWString prepend_before_reason;

  // With the -B option, one specify the number of connection table errors
  // llowed before programme exit.

  int connection_table_errors_allowed = 0;

  IWString output_file_stem;

  // We can append any arbitrary text to the name of each molecule

  IWString text_to_append;

  Molecule_Output_Object rejections_output_object;

  IWString_and_File_Descriptor stream_for_demerit_protos;

  // Functions.
  private:
    int SetupProtoRejectionFile(IWString& fname);
    int SetupRejectionFile(Command_Line& cl);
    int WriteProtoRejection(Molecule& m, const IWString& reason, IWString_and_File_Descriptor& output);

  public:
    int Initialise(Command_Line& cl);

    int GetOutputFileStem(Command_Line& cl, char flag, Molecule_Output_Object& output);

    // Returns 1 if opening `output` does not collide with `fname`.
    // This is not a great idea, should just restrict things to a single
    // output stream for all input files.
    int OkFileName(Molecule_Output_Object& output, const char * fname);

    int NonPrintingCharactersInName(Molecule & m);

    int AppendRejectionReason(Molecule & m,
                      const IWString & rejection_reason) const;
    int HandleRejectedMolecule(Molecule& m, const IWString& rejection_reason);
};

int
Parameters::SetupProtoRejectionFile(IWString& fname) {
  IWString with_suffix;
  with_suffix << fname << ".txtproto";
  if (! stream_for_demerit_protos.open(with_suffix.null_terminated_chars())) {
    cerr << "Parameters::SetupProtoRejectionFile:cannot open proto rejection file '" << with_suffix << "'\n";
    return 0;
  }
  return 1;
}

int
Parameters::SetupRejectionFile(Command_Line& cl) {
  IWString l;
  bool write_proto_output = false;
  IWString reject_log_file_name;
  for (int i = 0; cl.value('L', l, i); ++i) {
    if (l == "tsub") {
      append_rejection_reason_to_name = 1;
      write_rejection_reason_like_tsubstructure = 1;
      if (verbose) {
        cerr << "Rejection reasons written like tsubstructure\n";
      }
    } else if (l == "appreason") {
      append_rejection_reason_to_name = 1;
      if (verbose) {
        cerr << "The reason for rejection will be appended to the molecule name\n";
      }
    } else if (l == "proto") {
      write_proto_output = true;
    } else {
      reject_log_file_name = l;
    }

    if (reject_log_file_name.empty()) {
      cerr << "No rejection stem specified (-L)\n";
      return 0;
    }

    if (write_proto_output) {
      SetupProtoRejectionFile(reject_log_file_name);
    }

    if (! cl.option_present('o')) {
      rejections_output_object.add_output_type(FILE_TYPE_SMI);
    }
    else if (! rejections_output_object.determine_output_types(cl)) {
      cerr << "Cannot discern output types for rejections file\n";
      return 0;
    }

    if (! rejections_output_object.new_stem(reject_log_file_name)) {
      cerr << "Rejections file cannot use stem '" << reject_log_file_name << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Rejected structures will be written to '" << reject_log_file_name << "'\n";

    // Support for deprecated options.
    if (cl.option_present('a')) {
      append_rejection_reason_to_name = 1;
      if (verbose)
        cerr << "The reason for rejection will be appended to the molecule name\n";
    }

    if (cl.option_present('u')) {
      append_rejection_reason_to_name = 1;
      write_rejection_reason_like_tsubstructure = 1;
      if (verbose)
        cerr << "Rejection reasons written like tsubstructure\n";
    }
  }

  return 1;
}

int
Parameters::Initialise(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      return 0;
    }
  }

  if (cl.option_present('L')) {
    if (! SetupRejectionFile(cl)) {
      cerr << "Parameters::Initialise:cannot initialise -L conditions\n";
      return 0;
    }
  }

  if (cl.option_present('K')) {
    prepend_before_reason = cl.string_value('K');

    append_rejection_reason_to_name = 1;

    if (verbose)
      cerr << "Will prepend '" << prepend_before_reason << "' before rejection reasons\n";

    if (! prepend_before_reason.ends_with(' '))
      prepend_before_reason << ' ';
  }

  if (cl.option_present('w')) {
    run_all_checks = 1;

    if (verbose)
      cerr << "Will run all checks\n";
  }

  if (cl.option_present('B')) {
    if (1 != cl.option_count('B')) {
      cerr << "Only one -B option is allowed\n";
      Usage(57);
    }

    if (! cl.value('B', connection_table_errors_allowed) ||
          connection_table_errors_allowed < 0) {
      cerr << "The -B option requires a non negative integer argument\n";
      Usage(34);
    }

    if (verbose)
      cerr << connection_table_errors_allowed << " connection table errors allowed\n";
  }

  if (cl.option_present('P')) {
    cl.value('P', text_to_append);
    if (verbose)
      cerr << "Will append '" << text_to_append << "' to the name of each molecule written\n";

    if (! text_to_append.starts_with(' '))
      text_to_append.insert_at_beginning(' ');
  }

  if (cl.option_present('q')) {
  }

  return 1;
}

int
Parameters::AppendRejectionReason(Molecule & m,
                      const IWString & rejection_reason) const {
  IWString tmp = m.name();
  tmp += ' ';

  if (write_rejection_reason_like_tsubstructure) {
    tmp << "(1 matches to '" << rejection_reason << "')";
  } else {
    if (prepend_before_reason.length())
      tmp << prepend_before_reason;
    tmp += rejection_reason;
  }

  m.set_name(tmp);

  return 1;
}

#ifdef NOT_NEEDED
int
apply_all_filters(Molecule & m, int molecule_number)
{
  IWString rejection_reason;
  int rc = ApplyAllFiltersInner(m, rejection_reason);

  if (run_all_checks && rejection_reason.length() > 0)
    ;
  else if (rc)      // molecule is OK.
    return rc;

  if (append_rejection_reason_to_name && rejection_reason.length() > 0)
    AppendRejectionReason(m, rejection_reason);

  if (rejections_output_object.good())
    rejections_output_object.write(m);

  return 1;
}
#endif

/*
  We return the number of non-printing characters in the name
*/

int
Parameters::NonPrintingCharactersInName(Molecule & m) {
  IWString mname = m.name();

  int rc = 0;
  for (int i = mname.length() - 1; i >= 0; i--)
  {
    char c = mname[i];

    if (isprint(c))
      continue;

    if ('\0' != translate_non_printing_chars) {
      mname[i] = translate_non_printing_chars;
    }
    else if (remove_non_printing_chars) {
      mname.remove_item(i);
      i++;
    }

    rc++;
  }

  if (rc == 0) {
    return 0;
  }

  if ('\0' != translate_non_printing_chars || remove_non_printing_chars)
    m.set_name(mname);

  if (verbose > 1)
    cerr << "Removed/changed " << rc << " non printing chars in '" << mname << "'\n";

  return rc;
}

int
Parameters::GetOutputFileStem(Command_Line& cl, char flag,
                              Molecule_Output_Object& output) {
  if (! cl.option_present(flag)) {
    return 1;
  }
  cl.value (flag, output_file_stem);

  if (output.would_overwrite_input_files(cl, output_file_stem)) {
    cerr << "Cannot overwrite input file(s)\n";
    return 0;
  }

  if (verbose) {
    cerr << "New files will be created with stem '" << output_file_stem << "'\n";
  }

  return 1;
}

int
Parameters::OkFileName(Molecule_Output_Object& output,
                       const char * fname) {
  if (output_file_stem.length()) {
    if (output.would_use_name(output_file_stem.null_terminated_chars(), fname)) {
      cerr << "TpFirstPass: input '" << fname << "' and output stem '" << 
               output_file_stem << "' must not collide\n";
      return 0;
    }
  } else if (output.would_use_name(fname)) {
    cerr << "TpFirstPass: input '" << fname << "' and output must be distinct\n";
    return 0;
  }

  return 1;
}

int
Parameters::HandleRejectedMolecule(Molecule& m,
                const IWString& rejection_reason) {
  if (stream_for_demerit_protos.active()) {
    WriteProtoRejection(m, rejection_reason, stream_for_demerit_protos);
    // Or should we return here.
  }

  if (append_rejection_reason_to_name && rejection_reason.length() > 0)
    AppendRejectionReason(m, rejection_reason);

  if (rejections_output_object.good())
    rejections_output_object.write(m);

  return 1;
}

int
Parameters::WriteProtoRejection(Molecule& m,
     const IWString& reason,
     IWString_and_File_Descriptor& output) {
  MedchemRules::Molecule proto;
  proto.set_smiles(m.smiles().AsString());
  proto.set_name(m.name().AsString());
  MedchemRules::QueryMatch * q = proto.add_query_match();
  q->set_name(reason.AsString());
  q->set_rejected(true);

  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);
  std::string as_string;
  printer.PrintToString(proto, &as_string);
  output.write(as_string.data(), as_string.size());
  output << '\n';
  output.write_if_buffer_holds_more_than(8192);
  return 1;
}

int
TpFirstPass(Molecule & m,
            Parameters& params,
            lilly_medchem_rules::MCFirstPass& rules,
            lilly_medchem_rules::MCFirstPassCounter& counters,
            Molecule_Output_Object & output)
{
  if (params.chemical_standardisation.active()) {
    params.chemical_standardisation.process(m);
  }

  params.NonPrintingCharactersInName(m);

  params.natoms_accumulator.extra(m.natoms());   // Ignore atoms removed.

  IWString rejection_reason;
  if (rules.Rejected(m, counters, rejection_reason)) {
    return params.HandleRejectedMolecule(m, rejection_reason);
  }

  if (params.text_to_append.length())
    m.append_to_name(params.text_to_append);

  if (! output.write(m)) {
    return 0;
  }

  return 1;
}

int
TpFirstPass(data_source_and_type<Molecule> & input,
            Parameters& params,
            lilly_medchem_rules::MCFirstPass& rules,
            lilly_medchem_rules::MCFirstPassCounter& counters,
            Molecule_Output_Object & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    if (! TpFirstPass(*m, params, rules, counters, output)) {
      return 0;
    }
  }
  return 1;
}

int
TpFirstPass(const char *fname,
            FileType input_type,
            Parameters& params,
            lilly_medchem_rules::MCFirstPass& rules,
            lilly_medchem_rules::MCFirstPassCounter& counters,
            Molecule_Output_Object & output)
{
  assert (nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (params.verbose > 1)
    input.set_verbose(1);

  if (params.connection_table_errors_allowed)
    input.set_connection_table_errors_allowed(params.connection_table_errors_allowed);

// Avoid name collisions before they occur

  if (! params.OkFileName(output, fname)) {
    return 0;
  }

// Set up the output object for this file stem
// Make sure we deal properly with the case of multiple input files
// and a single output file (via the -S option)

  static int first_call = 1;

  int rc = 0;
  if (params.output_file_stem.length()) {
    if (first_call)
      rc = output.new_stem(params.output_file_stem);
    else
      rc = output.ok();

    first_call = 0;
  } else {
    rc = output.new_stem(fname);
  }

  if (rc == 0) {
    cerr << "Output object could not open file\n";
    return 0;
  }

  return TpFirstPass(input, params, rules, counters, output);
}

int
TpFirstPass(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:E:aI:g:t:L:S:K:X:c:C:Vi:o:r:R:B:P:p:b:kyue:x:Z:wf:F:H:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  if (! process_elements(cl)) {
    Usage(2);
  }

  if (! process_standard_aromaticity_options(cl)) {
    Usage(3);
  }

  Parameters params;
  if (! params.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  lilly_medchem_rules::MCFirstPass rules;

  if (! rules.Build(cl)) {
    cerr << "Cannot initialise first pass\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  }

  if (0 == input_type && ! all_files_recognised_by_suffix(cl))
    return 4;

  Molecule_Output_Object output;

  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl)) {
    cerr << "Cannot determine output types\n";
    Usage(8);
  }

  if (! params.GetOutputFileStem(cl, 'S', output)) {
    cerr << "Cannot establish output (-S)\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << prog_name << ": insufficient arguments " << argc << "\n";
    Usage(56);
  }

  lilly_medchem_rules::MCFirstPassCounter counters;

  for (const char * fname : cl) {
    if (params.verbose)
      cerr << prog_name << " processing '" << fname << "'\n";

    if (!TpFirstPass(fname, input_type, params, rules, counters, output)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  if (params.verbose) {
    cerr << params.natoms_accumulator.n() << " molecules had between " <<
          params.natoms_accumulator.minval() << " and " <<
          params.natoms_accumulator.maxval() << " atoms, mean " << 
          params.natoms_accumulator.average() << '\n';

    rules.Report(counters, cerr);
  }

  return 0;
}

}  // namespace tp_first_pass

int
main(int argc, char **argv)
{
  prog_name = argv[0];

  int rc = tp_first_pass::TpFirstPass(argc, argv);

  return rc;
}
