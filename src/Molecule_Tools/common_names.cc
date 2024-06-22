/*
  We have a file of structures and want all structures which
  are the same to have the same name
*/

#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/common_names.pb.h"

namespace common_names {

using std::cerr;

int verbose = 0;

int reduce_to_largest_fragment = 0;

int compare_as_graph = 0;

int remove_isotopes = 0;

int exclude_chiral_info = 0;

int exclude_cis_trans_info = 0;

IWString separator(':');

int molecules_read = 0;

Report_Progress report_progress;

Element_Transformations element_transformations;

/*
  In order to avoid recomputing the unique smiles, we store the
  unique smiles for each molecule
*/

int max_molecules = 0;

Chemical_Standardisation chemical_standardisation;

IWString* usmi = nullptr;

int write_duplicate_molecules = 0;

int first_name_and_count = 0;

// Optionally append a column with the number of examples.
int append_count = 0;

IWString_and_File_Descriptor stream_for_protos;

// For each unique smiles, data about the duplicate molecules
// found.
struct Instances {
 public:
  // As duplicates are found, we may build up a concatenated name.
  IWString concatenated_names;
  // Or we might just accumulate the number of instances.
  int number_instances;

 public:
  Instances(const IWString& name);

  void AddName(const IWString& name);

  void IncrementCount() {
    ++number_instances;
  }

  const IWString& names() const {
    return concatenated_names;
  }

  int FirstNameAndCount(Molecule& m) const;
};

Instances::Instances(const IWString& name) {
  concatenated_names = name;
  number_instances = 1;
}

void
Instances::AddName(const IWString& name) {
  concatenated_names.append_with_spacer(name, separator);
  ++number_instances;
}

int
Instances::FirstNameAndCount(Molecule& m) const {
  IWString new_name;
  new_name << concatenated_names << ' ' << number_instances;
  m.set_name(new_name);

  return 1;
}

/*
  The key will be the unique smiles, the value will be the molecule name
*/

// static IW_STL_Hash_Map_String usmi_hash;
IW_STL_Hash_Map<IWString, Instances> usmi_hash;

time_t tzero = static_cast<time_t>(0);

void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << "  -a             compare graph forms - add 2nd -a option to include H count\n";
  cerr << "  -c             exclude chirality information\n";
  cerr << "  -x             exclude directional bonds\n";
  cerr << "  -l             strip to largest fragment\n";
  cerr << "  -I             remove isotopes before storing\n";
  cerr << "  -T ...         standard element transformation options, enter '-T help'\n";
  cerr << "  -D <separator> separator for when storing duplicate entries\n";
  cerr << "  -f             single pass operation, smiles output only\n";
  cerr << "  -y             write first name and count of smiles only\n";
  cerr << "  -s <size>      maximum number of molecules to process\n";
  cerr << "  -r <number>    report progress every <number> molecules processed\n";
  cerr << "  -X ...         miscellaneous and obscure options\n";
  cerr << "  -S <name>      specify name for output\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type\n";
  cerr << "  -E ...         standard element options\n";
  (void)display_standard_aromaticity_options(cerr);
  (void)display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";

  exit(rc);
}

/*
  All the various transformations get applied very carefully.
  We may be comparing based on the isotope removed unique smiles,
  but we want to preserve the isotopes in the molecules written.
  But even that seems wrong, because molecules that were compared
  as being the same may have different isotops, but it will be the
  first molecule that gets written.... Oh well...
*/

void
do_conversions_needed_for_unique_smiles_generation(Molecule& m, int& hcount) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (compare_as_graph) {
    hcount = m.implicit_hydrogens();
    m.change_to_graph_form();
  }

  if (remove_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  return;
}

class Molecule_to_be_Written {
 private:
  IWString _smiles;
  IWString _mname;
  IWString _unique_smiles;

 public:
  int initialise(Molecule& m);
  int do_write(IWString_and_File_Descriptor&) const;
};

int
generate_unique_smiles(Molecule& m, int hcount, IWString& usmi) {
  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  if (compare_as_graph > 1) {
    IWString tmp = m.unique_smiles();
    tmp << ':' << hcount;
    usmi = tmp;
  } else {
    usmi = m.unique_smiles();
  }

  return 1;
}

int
update_global_unique_smiles_to_name_hash(const IWString& usmi, const IWString& mname) {
  auto f = usmi_hash.find(usmi);

  if (f == usmi_hash.end()) {
    Instances instances(mname);

    if (verbose > 2) {
      cerr << mname << " is unique, '" << usmi << "'\n";
    }
    usmi_hash.emplace(usmi, instances);

    return 1;
  }

  if (verbose > 1) {
    cerr << mname << " duplicate with '" << f->second.names() << "'\n";
  }

  if (first_name_and_count) {
    f->second.IncrementCount();
  } else {
    f->second.AddName(mname);
  }

  return 1;
}

int
Molecule_to_be_Written::initialise(Molecule& m) {
  _smiles = m.smiles();
  _mname = m.name();

  int hcount;

  do_conversions_needed_for_unique_smiles_generation(m, hcount);

  generate_unique_smiles(m, hcount, _unique_smiles);

  return update_global_unique_smiles_to_name_hash(_unique_smiles, _mname);
}

int
Molecule_to_be_Written::do_write(IWString_and_File_Descriptor& output) const {
  auto f = usmi_hash.find(_unique_smiles);

  if (f == usmi_hash.end()) {  // already written
    return 1;
  }

  output << _smiles << ' ';

  output << (*f).second.names();

  usmi_hash.erase(_unique_smiles);

  if (append_count) {
    output << ' ' << f->second.number_instances;
  }

  output << '\n';

  return 1;
}

void
preprocess(Molecule& m) {
  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

int
establish_names(Molecule& m, resizable_array_p<Molecule_to_be_Written>& mtbw) {
  Molecule_to_be_Written* t = new Molecule_to_be_Written;
  if (!t->initialise(m)) {
    delete t;
    return 0;
  }

  mtbw.add(t);

  return 1;
}

int
common_names_single_pass(data_source_and_type<Molecule>& input,
                         resizable_array_p<Molecule_to_be_Written>& mtbw) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!establish_names(*m, mtbw)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  if (verbose) {
    cerr << "Read " << mtbw.number_elements() << " molecules, " << usmi_hash.size()
         << " unique molecules\n";
  }

  return mtbw.number_elements();
}

int
write_molecules_to_be_written(const resizable_array_p<Molecule_to_be_Written>& mtbw,
                              IWString_and_File_Descriptor& output) {
  int n = mtbw.number_elements();

  if (verbose) {
    cerr << "Writing " << n << " molecules\n";
  }

  for (int i = 0; i < n; i++) {
    mtbw[i]->do_write(output);

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
write_molecules_to_be_written(Command_Line& cl,
                              const resizable_array_p<Molecule_to_be_Written>& mtbw) {
  if (cl.option_present('S')) {
    IWString fname = cl.option_value('S');
    if (!fname.ends_with(".smi")) {
      fname << ".smi";
    }

    IWString_and_File_Descriptor output;

    if (!output.open(fname.null_terminated_chars())) {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }

    return write_molecules_to_be_written(mtbw, output);
  } else {
    IWString_and_File_Descriptor output(1);

    return write_molecules_to_be_written(mtbw, output);
  }
}

int
common_names_single_pass(const char* fname, FileType input_type,
                         resizable_array_p<Molecule_to_be_Written>& mtbw) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return common_names_single_pass(input, mtbw);
}

int
WriteAsProto(Molecule& m, const IWString& usmi, const Instances& instance,
             IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(1);

  common_names::CommonNames proto;
  proto.set_smiles(m.smiles().AsString());
  proto.set_key(usmi.AsString());

  // Unfortunately we must unconcatenate the names.
  int i = 0;
  IWString id;
  while (instance.concatenated_names.nextword(id, i, separator[0])) {
    proto.add_id(id.AsString());
  }

  std::string as_string;
  printer.PrintToString(proto, &as_string);
  output << as_string << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
CommonNames(Molecule& m, Molecule_Output_Object& output) {
  IWString& u = usmi[molecules_read];

  auto f = usmi_hash.find(u);

  if (f == usmi_hash.end()) {
    if (!write_duplicate_molecules) {
      return 1;
    }

    cerr << "Huh: unique smiles '" << u << "' not in hash\n";
    return 0;
  }

  if (stream_for_protos.active()) {
    WriteAsProto(m, u, f->second, stream_for_protos);
    if (!write_duplicate_molecules) {
      usmi_hash.erase(u);
    }
    return 1;
  }

  if (first_name_and_count) {
    f->second.FirstNameAndCount(m);
  } else {
    m.set_name((*f).second.names());
    if (append_count) {
      m << ' ' << f->second.number_instances;
    }
  }

  if (!write_duplicate_molecules) {
    usmi_hash.erase(u);
  }

  return output.write(m);
}

int
read_dash_p_file(const IWString& usmi, const IWString& mname) {
  auto f = usmi_hash.find(usmi);

  if (f != usmi_hash.end()) {
    cerr << "Ignoring duplicate molecule '" << mname << "' in -p file\n";
    return 1;
  }

  Instances instances(mname);

  usmi_hash.emplace(usmi, std::move(instances));

  return 1;
}

int
read_dash_p_file(Molecule& m, const IWString& append_to_dash_p) {
  if (0 == append_to_dash_p.length()) {
    return read_dash_p_file(m.unique_smiles(), m.name());
  }

  IWString tmp(m.name());
  tmp << append_to_dash_p;

  return read_dash_p_file(m.unique_smiles(), tmp);
}

int
read_dash_p_file(data_source_and_type<Molecule>& input,
                 const IWString& append_to_dash_p) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    preprocess(*m);

    int hcount;

    do_conversions_needed_for_unique_smiles_generation(*m, hcount);

    if (!read_dash_p_file(*m, append_to_dash_p)) {
      return 0;
    }
  }

  return usmi_hash.size();
}

int
read_dash_p_file(const IWString& fname, FileType input_type,
                 const IWString& append_to_dash_p) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open -p file '" << fname << "'\n";
    return 0;
  }

  return read_dash_p_file(input, append_to_dash_p);
}

int
CommonNames(data_source_and_type<Molecule>& input, Molecule_Output_Object& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    (void)CommonNames(*m, output);
  }

  return 1;
}

int
CommonNames(FileType input_type, const char* fname, Molecule_Output_Object& output) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return CommonNames(input, output);
}

int
size_problem() {
  if (0 == max_molecules) {
    cerr << "No molecules\n";
    return 0;
  }

  if (max_molecules < 2) {
    cerr << "Strange to use common_names with just one molecule\n";
  }

  // usmi_hash.resize (max_molecules);     // how to check for failure?

  assert(nullptr == usmi);

  usmi = new IWString[max_molecules + 1];  // we don't use the 0 element

  if (nullptr == usmi) {
    cerr << "Cannot allocate " << max_molecules << " strings\n";
    return 0;
  }

  if (verbose) {
    cerr << "Problem sized for " << max_molecules << " molecules\n";
  }

  return 1;
}

int
establish_names(Molecule& m, int hcount) {
  if (molecules_read > max_molecules) {
    cerr << "Problem sized too small, max molecules = " << max_molecules << '\n';
    return 0;
  }

  IWString& usmi_ref = usmi[molecules_read];

  generate_unique_smiles(m, hcount, usmi_ref);

  return update_global_unique_smiles_to_name_hash(usmi_ref, m.name());
}

int
establish_names(data_source_and_type<Molecule>& input) {
  if (exclude_chiral_info) {
    set_include_chiral_info_in_smiles(0);
  }

  if (exclude_cis_trans_info) {
    set_include_cis_trans_in_smiles(0);
  }

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    int hcount;

    do_conversions_needed_for_unique_smiles_generation(*m, hcount);

    if (!establish_names(*m, hcount)) {
      return 0;
    }

    if (report_progress()) {
      cerr << "Processed " << molecules_read
           << " molecules. Time = " << (time(NULL) - tzero) << '\n';
    }
  }

  if (exclude_chiral_info) {
    set_include_chiral_info_in_smiles(1);
  }

  if (exclude_cis_trans_info) {
    set_include_cis_trans_in_smiles(1);
  }

  return 1;
}

int
establish_names(FileType input_type, const char* fname) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  if (0 == max_molecules) {
    max_molecules = input.molecules_remaining();

    if (0 == max_molecules) {
      cerr << "Yipes, no molecules in the input\n";
      return 0;
    }

    if (!size_problem()) {
      return 0;
    }
  }

  return establish_names(input);
}

int
OpenStreamForProtos(IWString& fname, IWString_and_File_Descriptor& output) {
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "OpenStreamForProtos:cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Proto output written to '" << fname << "'\n";
  }

  return 1;
}

void
DisplayDashXOptions(std::ostream& output) {
  output << " -X num         append an extra column to the output with the number of "
            "exemplars\n";
  output << " -X proto=<fname>  write results as common_names::CommonNames textprotos to "
            "<fname>\n";

  ::exit(0);
}

int
CommonNames(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:E:ag:Ilo:D:s:S:zcr:T:Zxp:fyK:X:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl)) {
    cerr << "Cannot parse -A options\n";
    usage(1);
  }

  if (cl.option_present('K')) {
    if (!process_standard_smiles_options(cl, verbose, 'K')) {
      cerr << "Cannot initialise smiles options\n";
      return 1;
    }
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot parse -E option\n";
    usage(3);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, (verbose > 1), 'g')) {
      cerr << "Cannot initialise chemical standardisations (-g option)\n";
      usage(18);
    }
  }

  if (cl.option_present('D')) {
    cl.value('D', separator);
    if (verbose) {
      cerr << "Concatenated items separated by '" << separator << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The report progress (-r) option must be followed by a whole positive "
              "number\n";
      usage(14);
    }

    tzero = time(NULL);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will strip to largest fragment before doing lookup\n";
    }
  }

  if (cl.option_present('I')) {
    remove_isotopes = 1;
    if (verbose) {
      cerr << "Isotopes stripped\n";
    }
  }

  if (cl.option_present('a')) {
    compare_as_graph = cl.option_count('a');
    if (0 == verbose) {
      ;
    } else if (1 == compare_as_graph) {
      cerr << "Comparisons will be based on the molecular graph\n";
    } else {
      cerr << "Comparisons will be based on the molecular graph\n";
    }
  }

  if (cl.option_present('c')) {
    exclude_chiral_info = 1;
    if (verbose) {
      cerr << "Optical isomers will be considered to be duplicates\n";
    }
  }

  if (cl.option_present('Z') || cl.option_present('x')) {
    exclude_cis_trans_info = 1;
    if (verbose) {
      cerr << "Will exclude cis/trans bonding information\n";
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Cannot initialise element transformations (-T option)\n";
      usage(5);
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "num") {
        append_count = 1;
        if (verbose) {
          cerr << "Output will include count of exemplars\n";
        }
      } else if (x.starts_with("proto=")) {
        x.remove_leading_chars(6);
        IWString proto_fname = x;
        if (!OpenStreamForProtos(proto_fname, stream_for_protos)) {
          cerr << "Cannot initialise stream for protos '" << proto_fname << "'\n";
          return 1;
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
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (FILE_TYPE_INVALID == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 1;
  }

  if (cl.option_present('p')) {
    IWString append_to_dash_p;
    IWString fname;

    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if (p.starts_with("suffix=")) {
        append_to_dash_p = p;
        append_to_dash_p.remove_leading_chars(7);
      } else if (fname.length()) {
        cerr << "ONly one previously selected molecule file allowed";
        usage(4);
      } else {
        fname = p;
      }
    }

    if (0 == fname.length()) {
      cerr << "Must specify file name of previous molecules\n";
      usage(5);
    }

    if (!read_dash_p_file(fname, input_type, append_to_dash_p)) {
      cerr << "Cannot read previous molecules from '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('z')) {
    cerr << "The -z option is no longer needed\n";
    write_duplicate_molecules = 0;
    if (verbose) {
      cerr << "Duplicate molecules will not be written\n";
    }
  }

  if (cl.option_present('y')) {
    first_name_and_count = 1;

    if (verbose) {
      cerr << "Will write the first name and count of instances\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('f')) {
    if (cl.option_present('o')) {
      cerr << "Output specification(s) ignored with -f, smiles output only\n";
    }

    if (cl.option_present('s')) {
      cerr << "Problem size, -s option, ignored with -f\n";
    }

    resizable_array_p<Molecule_to_be_Written> mtbw;
    for (const char* fname : cl) {
      if (!common_names_single_pass(fname, input_type, mtbw)) {
        cerr << "Cannot read molecules from '" << fname << "'\n";
        return 1;
      }
    }

    if (!write_molecules_to_be_written(cl, mtbw)) {
      return 1;
    }

    return 0;
  }

  // Normal two phase operation

  if (stream_for_protos.active() && cl.option_present('S')) {
    cerr << "Proto output and regular -S output cannot both be used\n";
    usage(1);
  }

  if (stream_for_protos.active()) {
  } else if (!cl.option_present('S')) {
    cerr << "Must use the -S option to indicate output file name\n";
    usage(1);
  }

  Molecule_Output_Object output;

  if (!cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
    if (verbose > 1) {
      cerr << "Output defaults to .smi\n";
    }
  } else if (!output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    usage(11);
  }

  if (!cl.option_present('s') && cl.number_elements() > 1) {
    cerr << "Must specify the number of molecules in the files via the -s option\n";
    usage(4);
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', max_molecules) || max_molecules < 2) {
      cerr << "The -s option (max molecules to process) option must be followed by a "
              "whole number > 1\n";
      usage(13);
    }

    if (!size_problem()) {
      return 1;
    }
  }

  for (const char* fname : cl) {
    if (!establish_names(input_type, fname)) {
      cerr << "Cannot read '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules from " << cl.number_elements()
         << " files. Found " << usmi_hash.size() << " distinct structures\n";
  }

  if (cl.option_present('S'))  // always
  {
    const_IWSubstring s = cl.string_value('S');
    if (!output.new_stem(s, 1)) {
      cerr << "Cannot open output stream '" << s << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "output to '" << s << "'\n";
    }
  }

  molecules_read = 0;  // reset so we can access the stored unique smiles

  for (const char* fname : cl) {
    if (!CommonNames(input_type, fname, output)) {
      cerr << "Huh, error writing '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}

}  // namespace common_names

int
main(int argc, char** argv) {
  int rc = common_names::CommonNames(argc, argv);

  return rc;
}
