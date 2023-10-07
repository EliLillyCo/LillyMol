// Convert reaction smiles to OpenReaction form
// https://github.com/Open-Reaction-Database/ord-schema/blob/master/proto/reaction.proto

#include <fstream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"

#include "Molecule_Tools/rxn_to_openrxn.h"

#include "external/ord_schema/proto/reaction.pb.h"
#include "external/ord_schema/proto/dataset.pb.h"

namespace rxn_to_openrxn {

struct JobResults
{
  int bad_reactions_ignored = 0;
  int reactions_read = 0;
  int next_report = 0;
};

int verbose = 0;

const char* prog_name = nullptr;

void
usage(int rc)
{
  cerr << "Converts reaction smiles to Open-Reaction proto form\n";
  cerr << "  -r <number>   report progress every <number> reactions\n";
  cerr << "  -p            reaction grouping is via + rather than .\n";
  cerr << "  -D            write DebugString of each proto (debugging)\n";
  cerr << "  -S <fname>    Output file to produce\n";
  cerr << "  -d <name>     Write a Dataset with `name` as id\n";
  cerr << "  -z            ignore erroneous reaction data\n";
  cerr << "  -v            verbose output\n";
  exit(rc);
}

void
ReportProgress(const JobOptions& job_options, JobResults* job_results, std::ostream& output)
{
  if (job_results->reactions_read < job_results->next_report)
  {
    return;
  }

  output << "Processed " << job_results->reactions_read << " reactions\n";
  job_results->next_report += job_options.report;

  return;
}

int
rxn_to_openrxn(const JobOptions& job_options, JobResults& job_results, RXN_File& rxn,
               const const_IWSubstring& buffer,
               std::unique_ptr<ord::Dataset>& dataset, std::ostream& output)
{
  ord::Reaction ord_proto = BuildOrdReaction(job_options, buffer, rxn);

  if (job_options.debug_string)
  {
    cerr << ord_proto.DebugString() << '\n';
  }

  if (dataset) {
    ord::Reaction * r = dataset->add_reactions();
    *r = std::move(ord_proto);
    return 1;
  }

  if (! ord_proto.SerializeToOstream(&output))
  {
    cerr << "rxn_to_openrxn:cannot write\n";
    return 0;
  }

  ReportProgress(job_options, &job_results, cerr);

  return 1;
}

int
rxn_to_openrxn_record(const JobOptions& job_options, JobResults& job_results,
                      const const_IWSubstring& buffer,
                      std::unique_ptr<ord::Dataset>& dataset,
                      std::ostream& output)
{
  RXN_File rxn;
  if (! rxn.build_from_reaction_smiles(buffer, job_options.component_grouping_is_plus))
  {
    cerr << "rxn_to_openrxn_record:cannot parse reaction\n";
    return 0;
  }

  job_results.reactions_read++;

  return rxn_to_openrxn(job_options, job_results, rxn, buffer, dataset, output);
}

int
rxn_to_openrxn(const JobOptions& job_options, JobResults& job_results, iwstring_data_source& input,
               std::unique_ptr<ord::Dataset>& dataset,
               std::ostream& output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! rxn_to_openrxn_record(job_options, job_results, buffer, dataset, output))
    {
      cerr << "Cannot process '" << buffer << "'\n";
      if (! job_options.ignore_bad_reactions)
      {
        job_results.bad_reactions_ignored++;
        return 0;
      }
    }
  }

  return 1;
}

int
rxn_to_openrxn(const JobOptions& job_options, JobResults& job_results, const char* fname,
               std::unique_ptr<ord::Dataset>& dataset,
               std::ostream& output)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "rxn_to_openrxn:cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return rxn_to_openrxn(job_options, job_results, input, dataset, output);
}

int
rxn_to_openrxn(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:S:zr:Dpd:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }
  else
    set_auto_create_new_elements(1);

  JobOptions job_options;
  JobResults job_results;

  if (cl.option_present('z'))
  {
    job_options.ignore_bad_reactions = true;

    if (verbose)
      cerr << "Will skip errors\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', job_options.report) || job_options.report <= 0)
    {
      cerr << "INvalid report option (-r)\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will report progress every " << job_options.report << " reactions processed\n";

    job_results.next_report = job_options.report;
  }

  if (cl.option_present('D'))
  {
    job_options.debug_string = true;
    if (verbose)
      cerr << "Will write DebugString for each proto\n";
  }

  if (cl.option_present('p'))
  {
    job_options.component_grouping_is_plus = true;
    if (verbose)
      cerr << "Reagent component grouping is + rather than .\n";
  }

  if (cl.number_elements() == 0)
  {
    cerr << "Must specify one or more reaction smiles files\n";
    usage(1);
  }

  set_reasonable_formal_charge_range(-7, 7);

  std::ofstream output;

  if (cl.option_present('S'))
  {
    IWString fname;
    cl.value('S', fname);
    output.open(fname.null_terminated_chars());
    if (! output.good())
    {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 1;
    }
    if (verbose)
      cerr << "Output written to '" << fname << "'\n";
  }

  std::unique_ptr<ord::Dataset> dataset;

  if (cl.option_present('d')) {
    dataset = std::make_unique<ord::Dataset>();

    const char * d = cl.option_value('d');
    dataset->set_name(d);
  }

  for (int i = 0; i < cl.number_elements(); ++i)
  {
    if (! rxn_to_openrxn(job_options, job_results, cl[i], dataset, output))
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (dataset) {
    output << dataset->DebugString();
  }

  if (verbose)
  {
    cerr << "Read " << job_results.reactions_read << " reactions";
    if (job_options.ignore_bad_reactions)
      cerr << " skipped " << job_results.bad_reactions_ignored << " bad reactions";
    cerr << "\n";
  }

  return 0;
}

}    // namespace rxn_to_openrxn

int
main(int argc, char** argv)
{
  rxn_to_openrxn::prog_name = argv[0];

  int rc = rxn_to_openrxn::rxn_to_openrxn(argc, argv);

  return rc;
}
