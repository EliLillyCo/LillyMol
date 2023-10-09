// Sequentially filter a collection of molecules, discarding any
// molecule that is a minor variant of an existing molecule.
// The variants are specified by reactions
// Undecided if this would work better if we read in the entire set of
// molecules and then operated pairwise on them. For now, this seems to
// work OK

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/standardise.h"

namespace molecular_variants {
using std::cerr;

void
Usage(int rc) {
  cerr << "Filter a set of sorted molecules based on externally specified reactions\n";
  cerr << "If the product has been seen before, the starting molecule is discarded\n";
  cerr << " -a              read all the input into a vector of Molecules (recommended)\n";
  cerr << " -r <fname>      reaction file (proto only)\n";
  cerr << " -R <fname>      file containing reaction files\n";
  cerr << " -X <fname>      write discarded molecules to <fname>\n";
  cerr << " -p              in the -X file, write the predecessor on a separate line\n";
  cerr << " -g ...          standard chemical standardisation options\n";
  cerr << " -l              reduce to largest fragment\n";
  cerr << " -c              discard chirality\n";
  cerr << " -o <sep>        output separator\n";
  cerr << " -v              verbose output\n";
  exit(rc);
}

// When the reaction product is found, we record the Predecessor - the previously
// encountered molecule that is the reaction product.
struct Predecessor {
  IWString smiles;
  IWString name;

  int times_seen = 0;

  Predecessor() {
  }
  Predecessor(const Predecessor& rhs) {
    smiles = rhs.smiles;
    name = rhs.name;
  }
};

struct JobParameters {
  int verbose = 0;

  bool load_all_molecules = false;

  int reduce_to_largest_fragment = 0;

  int remove_chirality = 0;

  // Optionally we can perform the reactions on the reaction products.
  // Only tested this to value 1.
  int recursion = 0;

  Chemical_Standardisation chemical_standardisation;

  int molecules_read = 0;

  // If we encounter duplicates within the input set, record here.
  int starting_molecule_duplicates_discarded = 0;

  // Data on what has been encountered before. Indexed by unique smiles.
  IW_STL_Hash_Map<IWString, Predecessor> seen;

  // The reactions applied to each input molecule.
  resizable_array_p<IWReaction> reactions;

  // Data on how many times each reaction matches an input.
  extending_resizable_array<int> rxn_hits;

  IWString_and_File_Descriptor stream_for_rejected;

  bool predecessor_on_separate_line = false;

  char sep = ' ';
};

void
Preprocess(Molecule& m,
           JobParameters& job_parameters) {
  if (job_parameters.reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
  if (job_parameters.remove_chirality) {
    m.remove_all_chiral_centres();
  }
  if (job_parameters.chemical_standardisation.active()) {
    job_parameters.chemical_standardisation.process(m);
  }
}

IWString
MaybeQuoted(const IWString& s, char sep) {
  if (! s.contains(sep)) {
    return s;
  }

  IWString result;
  result << '"' << s << '"';
  return result;
}

// Initially intended this to hold more pieces of information,
// but for now, only one item is held.
struct Args {
  // If the molecule has been changed by a previous reaction, that
  // will be recorded here.
  resizable_array<IWString> changes;
};

// Molecule `m` has been discarded after being processed by `rxn`.
// `predecessor` is the previously encountered molecule that is the
// same as the reaction product.
// Writes `m` to job_parameters.stream_for_rejected if it is active.
int
HandleDiscarded(Molecule& m,
                const IWReaction& rxn,
                Predecessor& predecessor,
                const Args& args,
                JobParameters& job_parameters) {
  predecessor.times_seen++;

  char sep = job_parameters.sep;

  IWString_and_File_Descriptor& output = job_parameters.stream_for_rejected;

  if (output.is_open()) {
    output << m.smiles() << sep << m.name() << sep << MaybeQuoted(rxn.comment(), sep);
    for (const IWString& changed_by :  args.changes) {
      output << sep << changed_by;
    }
    if (job_parameters.predecessor_on_separate_line) {
      output << '\n';
    } else {
      output << sep;
    }
    output << predecessor.smiles << sep << predecessor.name << '\n';
  }
  return 1;
}

// Returns 1 if applying job_parameters.reactions[ndx] to `m` results in
// something that has already been encountered.
// Depending on the value of job_parameters.recursion and `depth` this
// may recurse.
int
ProductEncountered(Molecule& m,
                   int ndx,
                   int depth,
                   Args& args,
                   JobParameters& job_parameters) {
  IWReaction& rxn = *job_parameters.reactions[ndx];

  Substructure_Results sresults;
  const int nhits = rxn.substructure_search(m, sresults);
  if (job_parameters.verbose > 2) {
    cerr << nhits << " hits to reaction " << rxn.comment() << '\n';
  }

  if (nhits == 0) {
    return 0;
  }

  job_parameters.rxn_hits[ndx]++;

  for (int i = 0; i < nhits; ++i) {
    Molecule mcopy(m);
    if (! rxn.in_place_transformation(mcopy, sresults.embedding(i))) {
      continue;
    }
    const auto iter = job_parameters.seen.find(mcopy.unique_smiles());
    if (iter != job_parameters.seen.end()) {
      return HandleDiscarded(m, rxn, iter->second, args, job_parameters);
    }

    if (depth == job_parameters.recursion) {  // No recursion.
      continue;
    }

    if (depth == 0) {
      args.changes.resize_keep_storage(0);
    }
    args.changes << rxn.name();
    for (int j = 0; j < job_parameters.reactions.number_elements(); ++j) {
      if (ProductEncountered(mcopy, j, depth + 1, args, job_parameters)) {
        return 1;
      }
    }
  }

  return 0;
}


// Applies all the reactions on job_parameters to `m`.
// If the reaction products are all new, then write to `output`.
int
MolecularVariants(Molecule& m,
                  JobParameters& job_parameters,
                  std::ostream& output) {
  if (job_parameters.seen.contains(m.unique_smiles())) {
    job_parameters.starting_molecule_duplicates_discarded++;
    return 1;
  }

  for (int i = 0; i < job_parameters.reactions.number_elements(); ++i) {
    Args args;
    if (ProductEncountered(m, i, 0, args, job_parameters)) {
      return 1;
    }
  }

  // Not rejected by any predecessor, add to map and write.

  Predecessor p;
  p.smiles = m.smiles();
  p.name = m.name();
  p.times_seen = 1;
  job_parameters.seen.emplace(m.unique_smiles(), std::move(p));

  output << m.smiles() << job_parameters.sep << m.name() << '\n';
  return 1;
}

int
RemoveDuplicates(resizable_array_p<Molecule>& mols,
                 JobParameters& job_parameters) {
  int initial_mols = mols.number_elements();
  for (int i = initial_mols - 1; i >= 0; --i) {
    Molecule& m = *mols[i];
    const IWString& usmi = m.unique_smiles();
    if (job_parameters.seen.contains(usmi)) {
      mols.remove_item(i);
    } else {
      Predecessor p;
      p.smiles = m.smiles();
      p.name = m.name();
      p.times_seen = 1;
      job_parameters.seen.emplace(std::move(usmi), std::move(p));
    }
  }

  if (job_parameters.verbose) {
    cerr << "Removed " << (initial_mols - mols.size()) << " exact duplicates\n";
  }
  return initial_mols - mols.size();
}

int
MolecularVariantsAll(resizable_array_p<Molecule>& mols,
                      JobParameters& job_parameters,
                      std::ostream& output) {
  RemoveDuplicates(mols, job_parameters);

  for (int i = mols.size() - 1; i >= 0; --i) {
    Molecule& m = *mols[i];

    for (int j = 0; j < job_parameters.reactions.number_elements(); ++j) {
      Args args;
      if (ProductEncountered(m, j, 0, args, job_parameters)) {
        mols.remove_item(i);
        break;
      }
    }
  }

  for (Molecule* m : mols) {
    output << m->smiles() << job_parameters.sep << m->name() << '\n';
  }

  if (job_parameters.verbose) {
    cerr << "Wrote " << mols.size() << " molecules\n";
  }

  return 1;
}

int
MolecularVariantsAll(data_source_and_type<Molecule>& input,
                      JobParameters& job_parameters,
                      std::ostream& output) {
  resizable_array_p<Molecule> mols;
  mols.resize(10000);

  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    Preprocess(*m, job_parameters);
    mols.add(m);
  }

  job_parameters.molecules_read = mols.size();
  if (job_parameters.verbose) {
    cerr << "Read " << mols.size() << " molecules\n";
  }

  return MolecularVariantsAll(mols, job_parameters, output);
}

int
MolecularVariants(data_source_and_type<Molecule>& input,
                  JobParameters& job_parameters,
                  std::ostream& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    job_parameters.molecules_read++;
    Preprocess(*m, job_parameters);
    if (! MolecularVariants(*m, job_parameters, output)) {
      cerr << "Fatal error processing " << m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }
  }
  return 1;
}

int
MolecularVariants(const char * fname,
                  FileType input_type,
                  JobParameters& job_parameters,
                  std::ostream& output) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  if (job_parameters.load_all_molecules) {
    return MolecularVariantsAll(input, job_parameters, output);
  } else {
    return MolecularVariants(input, job_parameters, output);
  }
}

// Reads the ReactionProto from `fname`.
std::unique_ptr<IWReaction>
ReadReaction(IWString& fname) {
  std::optional<ReactionProto::Reaction> proto = iwmisc::ReadTextProto<ReactionProto::Reaction>(fname);
  if (! proto) {
    cerr << "Cannot read reaction " << fname << '\n';
    return nullptr;
  }

  std::unique_ptr<IWReaction> result = std::make_unique<IWReaction>();
  if (! result->ConstructFromProto(proto.value(), fname)) {
    cerr << "Cannot parse proto\n";
    cerr << proto.value().ShortDebugString() << '\n';
    return nullptr;
  }

  return result;
}

// Reads the reactions in the file `input` which is in `dirname` and
// places the result in `reactions`.
int
ReactionsFromFile(const IWString& dirname,
                  iwstring_data_source& input,
                  resizable_array_p<IWReaction>& reactions) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    IWString fname;
    fname << dirname << buffer;
    std::unique_ptr<IWReaction> rxn = ReadReaction(fname);
    if (! rxn) {
      cerr << "ReactionsFromFile:cannot read " << fname << '\n';
      return 0;
    }

    reactions << rxn.release();
  }

  return 1;
}

// Reads the reactions specified in `fname` and places the results
// into `reactions`.
int
ReactionsFromFile(IWString& fname,
                  resizable_array_p<IWReaction>& reactions) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ReactionsFromFile:cannot open " << fname << '\n';
    return 0;
  }

  input.set_skip_blank_lines(1);
  input.set_ignore_pattern("^#");

  IWString dirname;
  if (fname.contains('/')) {
    dirname = fname;
    int length = dirname.rindex('/');
    dirname.iwtruncate(length + 1);
  } else {
    dirname = "./";
  }

  return ReactionsFromFile(dirname, input, reactions);
}

int
MolecularVariants(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vi:A:g:lcr:R:X:b:apo:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_SMI;
  if (!  process_standard_aromaticity_options(cl)) {
    cerr << "Cannot process aromaticity specifications\n";
    return 1;
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      return 1;
    }
  }

  JobParameters job_parameters;

  job_parameters.verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    job_parameters.reduce_to_largest_fragment = 1;
    if (job_parameters.verbose) {
      cerr << "Will strip to the largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    job_parameters.remove_chirality = 1;
    if (job_parameters.verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o)) {
      cerr << "Unrecognised character specification '" << o << "'\n";
      return 1;
    }
    job_parameters.sep = o[0];
    if (job_parameters.verbose) {
      cerr << "output separator '" << job_parameters.sep << "'\n";
    }
  }

  if (cl.option_present('g')) {
    if (! job_parameters.chemical_standardisation.construct_from_command_line(cl, false, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 1;
    }
  }

  if (cl.option_present('r')) {
    IWString fname;
    for (int i = 0; cl.value('r', fname, i); ++i) {
      std::unique_ptr<IWReaction> rxn = ReadReaction(fname);
      if (! rxn) {
        cerr << "Cannot read reaction from " << fname << '\n';
        return 1;
      }
      job_parameters.reactions << rxn.release();
    }
  }

  if (cl.option_present('R')) {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (! ReactionsFromFile(fname, job_parameters.reactions)) {
        cerr << "Cannot read reactions from " << fname << '\n';
        return 1;
      }
    }
  }

  if (cl.option_present('a')) {
    job_parameters.load_all_molecules = true;
    if (job_parameters.verbose) {
      cerr << "All molecules loaded into memory\n";
    }
  }

  if (cl.option_present('p')) {
    job_parameters.predecessor_on_separate_line = true;
    if (job_parameters.verbose) {
      cerr << "In the rejected molecules file (-X) write predecessor on separate line\n";
    }
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', job_parameters.recursion) || job_parameters.recursion < 0) {
      cerr << "Invalid recursion specification (-b)\n";
      Usage(1);
    }
    if (job_parameters.verbose)
      cerr << "Recursion set to " << job_parameters.recursion << '\n';
  }

  if (job_parameters.verbose) {
    cerr << "Read " << job_parameters.reactions.size() << " reactions\n";
  }

  if (job_parameters.reactions.size() == 0) {
    cerr << "Must specify reaction(s) via the -r and -R options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.option_present('X')) {
    const char * fname = cl.option_value('X');
    if (! job_parameters.stream_for_rejected.open(fname)) {
      cerr << "Cannot open stream for rejected molecules " << fname << '\n';
      return 1;
    }
    if (job_parameters.verbose) {
      cerr << "Rejected molecules written to " << fname << '\n';
    }
    const char sep = job_parameters.sep;

    job_parameters.stream_for_rejected << "smiles" << sep << "id" << sep <<
                "reaction" << sep << "product" << sep << "product_id\n";
  }

  for (const char * fname : cl) {
    if (! MolecularVariants(fname, input_type, job_parameters, std::cout)) {
      cerr << "Error processing " << fname << '\n';
      return 1;
    }
  }

  if (job_parameters.verbose) {
    if (! job_parameters.load_all_molecules) {  // Info has already been written.
      cerr << "Read " << job_parameters.molecules_read << " molecules\n";
    }
    cerr << job_parameters.starting_molecule_duplicates_discarded << " starting molecule duplicates discarded\n";
    cerr << "Generated " << job_parameters.seen.size() << " molecules\n";
    for (int i = 0; i < job_parameters.reactions.number_elements(); ++i) {
      cerr << job_parameters.rxn_hits[i] << " hits to reaction " << job_parameters.reactions[i]->comment() << '\n';
    }
  }

  return 0;
}

}  // namespace molecular_variants

int
main(int argc, char ** argv) {
  return molecular_variants::MolecularVariants(argc, argv);
}
