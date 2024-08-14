/*
  File conversion utility.
*/

#include <iostream>
#include <iostream>
#include <memory>
#include <limits>
#include <assert.h>


using std::cerr;

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"

#include "fileconv_opts.h"

static const char *prog_name;

namespace fileconv {

// When fileconv is processing, there can be a number of output streams.
struct LocalOptions {
  int verbose = 0;
  uint64_t molecules_read = 0;
  uint64_t molecules_changed = 0;
  uint64_t molecules_written = 0;
  int audit_input = 0;
  int debug_print_each_molecule = 0;

  // If a -S option is specified.
  IWString output_file_stem;

  // The main output stream - the -S option, or a per file stream.
  Molecule_Output_Object output;

  // Log messages about rejections.
  IWString_and_File_Descriptor reject_log;

  // Rejected molecules
  Molecule_Output_Object rejections_output;

  // Smiles before any filtering.
  // As part of the smiles infrastructure, I wanted to simultaneously
  // produce llyg.smi and llygO.smi. This allows that...
  IWString_and_File_Descriptor stream_for_smiles_before_filters;

  int connection_table_errors_allowed = 0;

  // Optionally we can write connection table errors to a destination.
  IWString connection_table_error_file;

  Report_Progress report_progress;

  // functions
  private:
    // If requested, write `smiles` to stream_for_smiles_before_filters
    int MaybeWritePostTransformSmiles(const IWString& smiles);

    // If there is a rejection stream open, write `m`.
    int MaybeWriteRejection(Molecule& m,
                             const IWString& rejection_reason);

    // Parse the -L option.
    int SetupRejectionsFile(Command_Line& cl, char flag);

    // Parse the -Y B4F= option.
    int OpenPreFilteredSmilesFile(Command_Line& cl, char flag);

    void SetDebugPrintEachMolecule(Command_Line& cl, char flag);

    int ParseBadHandlingoptions(Command_Line& cl, char flag);

    int SetupReporting(Command_Line& cl, char flag);

  public:
    // Initialise settings from command line options.
    int Build(Command_Line& cl);

    // Set verbosity, connection table errors... on a newly opened
    // input file.
    void SetupInputFile(data_source_and_type<Molecule>& input) const;

    // Does it seem likely that we will be able to open an output stream
    // for processing input from `fname`.
    int OutputFileSeemsOk(const char * fname) const;

    // If needed, open the output stream needed for processing `fname`.
    int MaybeOpenNewOutputFIle(const char * fname, int first_call);

    // A successful processing. If active, write it.
    int MaybeWrite(Molecule& m);

    // The most commonly used method.
    int Fileconv(Molecule& m, fileconv::FileconvConfig& config);

    int ReportResults(const Command_Line& cl, std::ostream& output) const;
};

int
LocalOptions::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (! SetupRejectionsFile(cl, 'L')) {
    cerr << "Cannot initialise rejection files (-L)\n";
    return 0;
  }

  if (! OpenPreFilteredSmilesFile(cl, 'Y')) {
    return 0;
  }

  SetDebugPrintEachMolecule(cl, 'Y');

  SetupReporting(cl, 'Y');

  audit_input = cl.option_present('a');

  if (! ParseBadHandlingoptions(cl, 'B')) {
    return 0;
  }

  // If just auditing the input, no output to be set up.
  if (audit_input) {
    return 1;
  }

  // Setup output options.
  if (cl.option_present('S')) {
    if (1 != cl.option_count('S')) {
      cerr << "Sorry, only one -S option possible\n";
      return 0;
    }

    cl.value('S', output_file_stem);
    if (verbose)
      cerr << "New files will be created with stem '" << output_file_stem << "'\n";
  }

  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);

    if (verbose)
      cerr << "Default output type smiles\n";
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output types\n";
    return 0;
  }

  if (cl.option_present('G')) {
    output.set_molecules_per_file(1);
    if (verbose)
      cerr << "Each molecule written to its own file (-G option)\n";
  }

  if ("-" == output_file_stem) {
  } else if (output.would_overwrite_input_files(cl, output_file_stem)) {
    cerr << "fileconv:cannot overwrite input file(s)\n";
    return 4;
  }

  return 1;
}

void
DisplayDashBOptions(std::ostream& os) {
  os << " -B <nn>        ignore as many as <nn> otherwise fatal input errors\n";
  os << " -B log=<fname> echo (unchanged) rejected input records to <fname>\n";

  exit(1);
}

int
LocalOptions::ParseBadHandlingoptions(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring b;
  for (int i = 0; cl.value(flag, b, i); ++i) {
    if (b.starts_with("log=")) {
      b.remove_leading_chars(4);
      connection_table_error_file = b;

      if (verbose)
        cerr << "Will log connection table errors to '" << connection_table_error_file << "'\n";
    } else if ("help" == b) {
      DisplayDashBOptions(cerr);
    } else if (b.numeric_value(connection_table_errors_allowed) &&
               connection_table_errors_allowed >= 0) {
      if (verbose)
        cerr << connection_table_errors_allowed << " connection table errors allowed\n";
    } else {
      cerr << "Invalid -B qualifier '" << b << "'\n";
      return 6;
    }
  }

  return 1;
}

void
LocalOptions::SetDebugPrintEachMolecule(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return;
  }
  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (y == "dbg") {
      debug_print_each_molecule = 1;
      return;
    }
  }
}

int
LocalOptions::SetupRejectionsFile(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  if (!cl.option_present('o'))
    rejections_output.add_output_type(FILE_TYPE_SMI);
  else if (!rejections_output.determine_output_types(cl)) {
    cerr << "Cannot discern output types for rejections file\n";
    return 0;
  }

  IWString reject_log_file_name = cl.option_value(flag);

  if (rejections_output.would_overwrite_input_files(cl, reject_log_file_name)) {
    cerr << "Reject file '" << reject_log_file_name << "' cannot overwrite input file(s)\n";
    return 0;
  }

  if (!rejections_output.new_stem(reject_log_file_name, 1)) {
    cerr << "Rejections file cannot use stem '" << reject_log_file_name << "'\n";
    return 0;
  }

  if (!reject_log.open(reject_log_file_name.null_terminated_chars())) {
    cerr << "Cannot open reject file '" << reject_log_file_name << "'\n";
    return 0;
  }

  if (cl.option_present('v')) { // verbose
    cerr << "Rejected structures will be written to '" << reject_log_file_name << "'\n";
  }

  return 1;
}

// Parse '-Y B4F=<fname>' and open stream_for_smiles_before_filters.
int
LocalOptions::OpenPreFilteredSmilesFile(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return 1;
  }
  IWString y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (! y.starts_with("B4F=")) {
      continue;
    }

    y.remove_leading_chars(4);
    IWString fname(y);
    if ('-' == y)
      ;
    else if (!fname.ends_with(".smi"))
      fname << ".smi";

    if (!stream_for_smiles_before_filters.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for pre-filtered smiles '" << fname << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Pre-filtered smiles written to '" << fname << "'\n";
    return 1;
  }

  return 1;
}

int
LocalOptions::SetupReporting(Command_Line& cl, char flag) {
  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (! y.starts_with("rpt=")) {
      continue;
    }
    y.remove_leading_chars(4);
    uint32_t rpt;
    if (! y.numeric_value(rpt)) {
      cerr << "Invalid -Y rpt=... directive '" << y << "'\n";
      return 0;
    }
    report_progress.set_report_every(rpt);
    return 1;
  }

  return 1;
}

void
LocalOptions::SetupInputFile(data_source_and_type<Molecule>& input) const {
  if (verbose > 1) {
    input.set_verbose(1);
  }

  if (connection_table_errors_allowed) {
    input.set_connection_table_errors_allowed(connection_table_errors_allowed);
    if (connection_table_error_file.length())
      input.set_connection_table_error_file(connection_table_error_file);
  }
}

int
LocalOptions::OutputFileSeemsOk(const char * fname) const {
  if (audit_input) {
    return 1;
  }

  // output file names generated by the output object
  if (output.name_token_for_file_name() >= 0)  {
    return 1;
  }

  // Single output file for each input.
  if (output_file_stem.length()) {
    return 1;
  }

  if (output.would_use_name(fname)) {
    cerr << "LocalOptions::OutputFileSeemsOk: input '" << fname << "' and output must be distinct\n";
    return 0;
  }

  return 1;
}


int
LocalOptions::MaybeOpenNewOutputFIle(const char * fname, int first_call) {
  if (audit_input) {  // Nothing to do.
    return 1;
  }

  // One file for all input files, open once.
  if (output_file_stem.length()) {
    if (first_call)
      return output.new_stem(output_file_stem);
    else
      return output.ok();
  }

  // New output file for each molecule.
  if (output.name_token_for_file_name() >= 0) {
    return 1;
  }

  // New output file for this input.
  if (!  output.new_stem(fname)) {
    return 0;
  }

  return 1;
}

int
LocalOptions::MaybeWritePostTransformSmiles(const IWString& smiles) {
  if (smiles.empty()) {
    return 1;
  }
  if (! stream_for_smiles_before_filters.is_open()) {
    return 1;
  }

  // organic subset
  stream_for_smiles_before_filters << smiles << '\n';
  stream_for_smiles_before_filters.write_if_buffer_holds_more_than(32768);

  return stream_for_smiles_before_filters.good();
}

int
LocalOptions::MaybeWriteRejection(Molecule& m,
                             const IWString& rejection_reason) {
  if (reject_log.is_open()) {
    reject_log << molecules_read << " '" << m.name() << "' REASON " << rejection_reason << '\n';
    reject_log.write_if_buffer_holds_more_than(4096);
  }

  if (rejections_output.active()) {
    m << " REJ " << rejection_reason;
    rejections_output.write(m);
  }

  return 1;
}

int
LocalOptions::ReportResults(const Command_Line& cl, std::ostream& output) const {
  if (cl.number_elements() > 1)
    output << molecules_read << " molecules read, ";
  if (!audit_input)
    output << molecules_written << " molecules written";
  if (molecules_changed)
    output << ' ' << molecules_changed << " molecules changed";
  output << '\n';

  return output.good();
}

int
LocalOptions::Fileconv(Molecule& m,
                       fileconv::FileconvConfig& config) {

  if (report_progress()) {
    cerr << "Read " << molecules_read << " molecules";
    if (! audit_input) {
      cerr << ", wrote " << molecules_written;
      if (molecules_changed) {
        cerr << molecules_changed << " molecules changed";
      }
    }
    cerr << '\n';
  }

  if (debug_print_each_molecule) {
    m.compute_aromaticity_if_needed();
    m.debug_print(cerr);
  }
    
  FileconvResult result = config.Process(m);
  if (result.error) {
    cerr << "Fatal error processing " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  MaybeWritePostTransformSmiles(result.smiles_before_filters);

  if (result.rejected) {
    return MaybeWriteRejection(m, result.rejection_reason);
  }

  if (result.molecule_changed) {
    molecules_changed++;
  }

  if (audit_input) {  // Nothing more to do.
    return 1;
  }

  if (! output.write(m)) {
    return 0;
  }

  molecules_written++;

  return 1;
}

int
Fileconv(data_source_and_type<Molecule> & input,
         fileconv::FileconvConfig& config,
         LocalOptions& local_options) {
  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    local_options.molecules_read++;

    if (! local_options.Fileconv(*m, config)) {
      return 0;
    }
  }

  if (input.stopped_because_of_error())
    return 0;

  return 1;
}

int
Fileconv(const char *fname,
         FileType input_type,
         fileconv::FileconvConfig& config,
         LocalOptions & local_options) {
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  local_options.SetupInputFile(input);

  if (! local_options.OutputFileSeemsOk(fname)) {
    return 0;
  }
// Set up the output object for this file stem
// If there is a new stem for output files, make sure we get it.
// Make sure we deal properly with the case of multiple input files
// and a single output file (via the -S option)

  static bool first_call = true;

  if (! local_options.MaybeOpenNewOutputFIle(fname, first_call)) {
    cerr << "Cannot establish output for '" << fname << "'\n";
    first_call = false;
    return 0;
  }

  first_call = false;

  return Fileconv(input, config, local_options);
}

int
FileconvListOfFiles(iwstring_data_source & input,
                    FileType input_type,
                    fileconv::FileconvConfig& config,
                    LocalOptions& local_options) {
  input.set_strip_trailing_blanks(1);
  input.set_skip_blank_lines(1);
  input.set_strip_leading_blanks(1);

  IWString buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#'))
      continue;

    if (local_options.verbose > 1)
      cerr << "Processing '" << buffer << "'\n";

    if (! Fileconv(buffer.null_terminated_chars(), input_type, config, local_options)) {
      cerr << "Fatal error processing file '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
FileconvListOfFiles(const char * fname,
                    FileType input_type,
                    fileconv::FileconvConfig& config,
                    LocalOptions& local_options) {
  const char * tmp = fname + 2;

  iwstring_data_source input(tmp);

  if (! input.good()) {
    cerr << "Cannot open list of files to process '" << fname << "'\n";
    return 0;
  }

  return FileconvListOfFiles(input, input_type, config, local_options);
}

int
AllFilesRecognisedByType(const Command_Line& cl) {
  for (const char * fname: cl) {
    const_IWSubstring tmp(fname);
    if (tmp.starts_with("F:")) {
      continue;
    }
    if (!discern_file_type_from_name(fname)) {
      cerr << "Cannot determine file type from '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
Fileconv(int argc, char ** argv) {
  Command_Line cl(argc, argv, "PN:T:eB:I:s:g:h:H:t:n:L:S:GA:K:w:W:m:X:c:C:Q:O:E:f:F:vVi:ao:r:R:p:J:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    fileconv::Usage(1);
  }

  fileconv::FileconvConfig config;
  if (! config.Build(cl)) {
    cerr << "Cannot initialise fileconv\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;
  LocalOptions local_options;
  if (! local_options.Build(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      fileconv::Usage(6);
    }
  }

  if (FILE_TYPE_INVALID != input_type) {  // great, explicitly specified
  } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  // reading from a pipe, assume smiles input
    if (local_options.verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  } else if (AllFilesRecognisedByType(cl)) {
  } else {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    fileconv::Usage(1);
  }

  int rc = 0;
  for (const char * fname: cl) {
    const const_IWSubstring as_string(fname);

    if (local_options.verbose)
      cerr << "Processing '" << fname << "'\n";

    int myrc;
    if (as_string.length() > 2 && as_string.starts_with("F:"))
      myrc = FileconvListOfFiles(fname, input_type, config, local_options);
    else
      myrc = Fileconv(fname, input_type, config, local_options);

    if (!myrc) {
      rc = 1;
      break;
    }
  }

  if (local_options.verbose) {
    local_options.ReportResults(cl, cerr);
    config.ReportResults(cl, std::cerr);
  }

  return rc;
}

}  // namespace fileconv
//#define DLL_FLAG=1
#ifndef DLL_FLAG

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = fileconv::Fileconv(argc, argv);

  return rc;
}

#endif

#ifdef HZ_CSHARP

extern "C" int
fileconv_csharp(int argc, char **argv) {
  return fileconv::Fileconv(argc, argv);
}

#endif

