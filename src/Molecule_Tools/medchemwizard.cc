/*
  Implementation of Medchem Wizard tool
*/

#include <cstring>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/medchemwizard_lib.h"

using std::cerr;

const char* prog_name = nullptr;

static void
usage(int rc) {
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(  -R <fname>    file with list of reactions
  -D <maxdepth> recursively generated transformed molecules to <maxdepth> (def 0 - no recursion)
  -x <max_hits> allow no more than <max_hits> matches to any reaction
  -x quiet      do NOT report when multiple hits are truncated
  -U ...        uniqueness. Either 'all' or 'each'
  -m <natoms>   discard products having fewer than <natoms> atoms
  -M <natoms>   discard products having more  than <natoms> atoms
  -c            remove chirality from all input molecules
  -W <sep>      append reaction names to products
  -q <query>    atoms matching query must not be changed
  -s <smarts>   atoms matching SMARTS must not be changed
  -z i          ignore molecules not matching any -q/-s query
  -V .          discard molecules with bad valences
  -V <fname>    discard molecules with bad valences, write to <fname>
  -y            standardise molecules produced
  -l            reduce to largest fragment
  -i <type>     input specification
  -g ...        chemical standardisation options
  -E ...        standard element specifications
  -A ...        standard aromaticity specifications
  -v            verbose output
)";

  exit(rc);
}

static int
MedchemWizard(medchemwizard::MedchemWizard& wizard, data_source_and_type<Molecule>& input,
              IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (! wizard.ProcessToStream(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
MedchemWizard(medchemwizard::MedchemWizard& wizard, const char* fname,
              FileType input_type, IWString_and_File_Descriptor& output, int verbose) {
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return MedchemWizard(wizard, input, output);
}

static int
MedchemWizard(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lR:D:x:U:yV:cm:M:W:q:s:z:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  medchemwizard::MedchemWizard wizard;
  medchemwizard::Options& options = wizard.options();
  options.verbose = cl.option_count('v');
  const int verbose = options.verbose;

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!options.chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    options.reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    options.remove_all_chiral_centres = 1;
    if (verbose) {
      cerr << "Will remove all chirality from input molecules\n";
    }
  }

  if (cl.option_present('W')) {
    options.append_names = 1;

    cl.value('W', options.sep);

    char_name_to_char(options.sep, false /* no error messages */ );

    if (verbose) {
      cerr << " Will append reaction names, separtor '" << options.sep << "'\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', options.min_atoms) || options.min_atoms < 1) {
      cerr << "The minimum number of atoms option (-m) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will discard products having fewer than " << options.min_atoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', options.max_atoms) || options.max_atoms < options.min_atoms) {
      cerr << "The maximum number of atoms option (-M) must be a whole +ve number "
              "greater than "
           << options.min_atoms << '\n';
      usage(1);
    }

    if (verbose) {
      cerr << "Will discard products having more than " << options.max_atoms << " atoms\n";
    }
  }

  if (!cl.option_present('R')) {
    cerr << "Must specify file of reactions via the -R option\n";
    usage(1);
  }

  if (cl.option_present('D')) {
    if (!cl.value('D', options.max_depth) || options.max_depth < 1) {
      cerr << "The maximum depth option (-D) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "max depth " << options.max_depth << '\n';
    }
  }

  if (cl.option_present('x')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('x', x, i); ++i) {
      if ("quiet" == x) {
        options.report_multiple_hits_truncated = 0;
        if (verbose) {
          cerr << "Will suppress warnings about too many substructure matches\n";
        }
      } else if (!x.numeric_value(options.max_hits) || options.max_hits < 1) {
        cerr << "The maximum number of substructure matches (-x) must be a whole +ve "
                "number\n";
        usage(1);
      } else if (verbose) {
        cerr << "A maximum of " << options.max_hits << " substructure matches will be used\n";
      }
    }
  }

  if (cl.option_present('U')) {
    const_IWSubstring u = cl.string_value('U');

    if ("all" == u) {
      options.unique_across_all_molecules = 1;
    } else if ("each" == u) {
      options.unique_within_molecule = 1;
    } else {
      cerr << "Unrecognised -U qualifier '" << u << "'\n";
      usage(1);
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, wizard.do_not_change_queries(), verbose, 'q')) {
      cerr << "Cannot read do-not-change queries (-q)\n";
      return 1;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> query = std::make_unique<Substructure_Query>();
      if (! query->create_from_smarts(smarts)) {
        cerr << "Invalid do-not-change SMARTS (-s) '" << smarts << "'\n";
        return 1;
      }
      wizard.do_not_change_queries() << query.release();
    }
  }

  if (verbose && ! wizard.do_not_change_queries().empty()) {
    cerr << "Read " << wizard.do_not_change_queries().number_elements()
         << " do-not-change queries\n";
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i' || z == "ignore") {
        options.ignore_do_not_change_queries_not_matching = 1;
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(1);
      }
    }
  }

  if (cl.option_present('y')) {
    options.postprocess_molecules_produced = 1;

    if (verbose) {
      cerr << "Will post process molecules produced just like input molecules\n";
    }
  }

  IWString_and_File_Descriptor stream_for_bad_valence;
  if (cl.option_present('V')) {
    IWString v = cl.string_value('V');

    options.discard_bad_valences = 1;

    if ('.' != v) {
      if (!v.ends_with(".smi")) {
        v << ".smi";
      }

      if (!stream_for_bad_valence.open(v.null_terminated_chars())) {
        cerr << "Cannot open stream for bad valence '" << v << "'\n";
        return 1;
      }

      if (verbose) {
        cerr << "Molecules with bad valences discarded, written to '" << v << "'\n";
      }
    }
  }

  wizard.set_bad_valence_stream(&stream_for_bad_valence);

  const char* r = cl.option_value('R');
  if (! wizard.ReadReactions(r)) {
    cerr << "Cannot read reactions from '" << r << "'\n";
    return 1;
  }

  if (verbose) {
    cerr << "Read " << wizard.number_reactions() << " reactions from '" << r << "'\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!MedchemWizard(wizard, cl[i], input_type, output, verbose)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    const medchemwizard::Stats& stats = wizard.stats();
    cerr << "Read " << stats.molecules_read << " molecules, produced " << stats.molecules_produced
         << '\n';
    if (stats.truncated_to_max_hits > 0) {
      cerr << stats.truncated_to_max_hits << " substructure searches truncated to "
           << options.max_hits << '\n';
    }
    if (stats.duplicate_molecules_suppressed) {
      cerr << stats.duplicate_molecules_suppressed << " duplicate molecules suppressed\n";
    }

    const int* molecules_hitting_reaction = wizard.molecules_hitting_reaction();
    for (int i = 0; i < wizard.number_reactions(); ++i) {
      float f = static_cast<float>(molecules_hitting_reaction[i]) /
                static_cast<float>(stats.molecules_read);

      cerr << molecules_hitting_reaction[i] << " molecules hit " << wizard.reaction(i)->comment()
           << ' ' << f << '\n';
    }
    cerr << stats.bad_valences_discarded << " products with bad valences discarded\n";
    if (stats.discarded_for_too_many_atoms) {
      cerr << stats.discarded_for_too_many_atoms << " discarded_for_too_many_atoms\n";
    }
    if (stats.discarded_for_too_few_atoms) {
      cerr << stats.discarded_for_too_few_atoms << " discarded_for_too_few_atoms\n";
    }
    if (stats.molecules_not_matching_do_not_change_queries) {
      cerr << stats.molecules_not_matching_do_not_change_queries
           << " molecules did not match do-not-change queries\n";
    }
    if (stats.embeddings_rejected_for_changing_protected_atoms) {
      cerr << stats.embeddings_rejected_for_changing_protected_atoms
           << " reaction embeddings rejected for changing protected atoms\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = MedchemWizard(argc, argv);

  return rc;
}
