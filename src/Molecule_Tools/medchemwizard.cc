/*
  Implementation of Medchem Wizard tool
*/

#include <iostream>
#include <memory>
#include <limits>
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int max_depth = 0;

static int molecules_produced = 0;

static int max_hits = std::numeric_limits<int>::max();
static int truncated_to_max_hits = 0;

static int report_multiple_hits_truncated = 1;

static int unique_across_all_molecules = 0;
static int unique_within_molecule = 0;

static int duplicate_molecules_suppressed = 0;

static int postprocess_molecules_produced = 0;

static int discard_bad_balences = 0;

static int bad_valences_discarded = 0;

static IWString_and_File_Descriptor stream_for_bad_valence;

static int max_atoms = std::numeric_limits<int>::max();
static int min_atoms = 1;

static int discarded_for_too_many_atoms = 0;
static int discarded_for_too_few_atoms = 0;

static int remove_all_chiral_centres = 0;

static int append_names = 0;

static IWString sep;

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
  cerr << "  -R <fname>    file with list of reactions\n";
  cerr << "  -D <maxdepth> recursively generated transformed molecules to <maxdepth> (def 0 - no recursion)\n";
  cerr << "  -x <max_hits> allow no more than <max_hits> matches to any reaction\n";
  cerr << "  -x quiet      do NOT report when multiple hits are truncated\n";
  cerr << "  -U ...        uniqueness. Either 'all' or 'each'\n";
  cerr << "  -m <natoms>   discard products having fewer than <natoms> atoms\n";
  cerr << "  -M <natoms>   discard products having more  than <natoms> atoms\n";
  cerr << "  -c            remove chirality from all input molecules\n";
  cerr << "  -W <sep>      append reaction names to products\n";
  cerr << "  -V .          discard molecules with bad valences\n";
  cerr << "  -V <fname>    discard molecules with bad valences, write to <fname>\n";
  cerr << "  -y            standardise molecules produced\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
read_reaction(const const_IWSubstring & buffer,
              const IWString & dir,
              const Sidechain_Match_Conditions & smc,
              resizable_array_p<IWReaction> & rxn)
{
  IWString fname;
  fname << dir << '/' << buffer;

  iwstring_data_source input(fname.null_terminated_chars());

  if (! input.good())
  {
    cerr << "Cannot open reaction '" << fname << "'\n";
    return 0;
  }

  if (verbose > 2)
    cerr << "Reading '" << fname << "'\n";

  msi_object msi;
  msi.set_display_no_data_error_message(0);

  while(msi.read(input))
  {
    IWReaction * r = new IWReaction();
    if (! r->construct_from_msi_object(msi, smc))
    {
      cerr << "Cannot build reaction\n";
      cerr << msi;
      delete r;

      return 0;
    }

    rxn.add(r);
  }


  return 1;
}

static int
read_reactions(iwstring_data_source & input,
               const IWString & dir,
               const Sidechain_Match_Conditions & smc,
               resizable_array_p<IWReaction> & rxn)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (0 == buffer.length() || buffer.starts_with('#'))
      continue;

    if (! read_reaction(buffer, dir, smc, rxn))
    {
      cerr << "Cannot read reaction '" << buffer << "'\n";
      return 0;
    }
  }

  return rxn.number_elements();
}

static int
read_reactions(const char * fname,
               const Sidechain_Match_Conditions & smc,
               resizable_array_p<IWReaction> & rxn)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_reactions:cannot open '" << fname << "'\n";
    return 0;
  }

  IWString dir;

  const_IWSubstring tmp(fname);

  const int slash = tmp.rindex('/');

  if (slash >= 0)
  {
    dir = fname;
    dir.iwtruncate(slash);
  }
  else
    dir = "./";

  return read_reactions(input, dir, smc, rxn);
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
check_uniqueness(Molecule & m,
                 const resizable_array_p<IWReaction> & rxn,
                 const int query_just_run,
                 IW_STL_Hash_Set & smiles_hash)
{
  if (postprocess_molecules_produced)
    preprocess(m);

  if (discard_bad_balences && ! m.valence_ok())
  {
    if (verbose > 2)
      cerr << "bad valence " << m.smiles() << ' ' << m.name() << ", last rxn " << query_just_run << ' ' << rxn[query_just_run]->comment() << endl;
    bad_valences_discarded++;

    if (stream_for_bad_valence.is_open())
    {
      stream_for_bad_valence << m.smiles() << ' ' << m.name() << '\n';
      stream_for_bad_valence.write_if_buffer_holds_more_than(8192);
    }

    return 0;
  }

  const int matoms = m.natoms();

  if (matoms > max_atoms)
  {
    discarded_for_too_many_atoms++;

    return 0;
  }

  if (matoms < min_atoms)
  {
    discarded_for_too_few_atoms++;

    return 0;
  }

  if (0 == unique_across_all_molecules && 0 == unique_within_molecule)
    return 1;

  const auto & usmi = m.unique_smiles();

  if (smiles_hash.contains(usmi))
  {
    duplicate_molecules_suppressed++;
    return 0;
  }

  smiles_hash.insert(usmi);

  return 1;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
medchem_wizard(Molecule & m,
               const IWString & mname,
               resizable_array_p<IWReaction> & rxn,
               int * molecules_hitting_reaction,
               const int depth,
               IW_STL_Hash_Set & smiles_hash,
               IWString_and_File_Descriptor & output)
{
  Molecule_to_Match target(&m);

  const int n = rxn.number_elements();

  IWString product_molecule_name(mname);
  int initial_name_length = mname.length();

  for (int i = 0; i < n; ++i)
  {
    const auto r = rxn[i];

    Substructure_Results sresults;

    int nhits = r->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    molecules_hitting_reaction[i]++;

    if (nhits > max_hits)
    {
      if (report_multiple_hits_truncated)
        cerr << nhits << " matches to " << r->comment() << " in " << m.name() << " in " << m.name() << " truncated\n";
      nhits = max_hits;
      truncated_to_max_hits++;
    }

    for (int j = 0; j < nhits; ++j)
    {
      Molecule product;
      if (! r->perform_reaction(&m, sresults.embedding(j), product))
        continue;

      molecules_produced++;

      if (! check_uniqueness(product, rxn, i, smiles_hash))
        continue;

      if (append_names)
        product_molecule_name << sep << r->comment();

      output << product.smiles() << ' ' << product_molecule_name << '\n';

      output.write_if_buffer_holds_more_than(8192);

      if (depth < max_depth)
        medchem_wizard(product, product_molecule_name, rxn, molecules_hitting_reaction, depth+1, smiles_hash, output);

      product_molecule_name.resize(initial_name_length);
    }
  }

  return output.good();
}

static int
medchem_wizard(data_source_and_type<Molecule> & input,
               resizable_array_p<IWReaction> & rxn,
               int * molecules_hitting_reaction,
               IWString_and_File_Descriptor & output)
{
  IW_STL_Hash_Set smiles_hash;

  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    IWString mname(m->name());

    if (! medchem_wizard(*m, mname, rxn, molecules_hitting_reaction, 0, smiles_hash, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);

    if (unique_within_molecule)
      smiles_hash.clear();
  }

  return 1;
}

static int
medchem_wizard(const char * fname, FileType input_type, 
               resizable_array_p<IWReaction> & rxn,
               int * molecules_hitting_reaction,
               IWString_and_File_Descriptor & output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return medchem_wizard(input, rxn, molecules_hitting_reaction, output);
}

static int
medchem_wizard(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lR:D:x:U:yV:cm:M:W:");

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

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c'))
  {
    remove_all_chiral_centres = 1;
    if (verbose)
      cerr << "Will remove all chirality from input molecules\n";
  }

  if (cl.option_present('W'))
  {
    append_names = 1;

    cl.value('W', sep);

    if (verbose)
      cerr << " Will append reaction names, separtor '" << sep << "'\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_atoms) || min_atoms < 1)
    {
      cerr << "The minimum number of atoms option (-m) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will discard products having fewer than " << min_atoms << " atoms\n";
  }

  if (cl.option_present('M'))
  {
    if (! cl.value('M', max_atoms) || max_atoms < min_atoms)
    {
      cerr << "The maximum number of atoms option (-M) must be a whole +ve number greater than " << min_atoms << endl;
      usage(1);
    }

    if (verbose)
      cerr << "Will discard products having more than " << max_atoms << " atoms\n";
  }

  if (! cl.option_present('R'))
  {
    cerr << "Must specify file of reactions via the -R option\n";
    usage(1);
  }

  Sidechain_Match_Conditions smc;

  resizable_array_p<IWReaction> rxn;

  if (cl.option_present('R'))
  {
    const char * r = cl.option_value('R');

    if (! read_reactions(r, smc, rxn))
    {
      cerr << "Cannot read reactions from '" << r << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << rxn.size() << " reactions from '" << r << "'\n";
  }

  if (cl.option_present('D'))
  {
    if (! cl.value('D', max_depth) || max_depth < 1)
    {
      cerr << "The maximum depth option (-D) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "max depth " << max_depth << endl;
  }

  if (cl.option_present('x'))
  {
    const_IWSubstring x;
    for (int i = 0; cl.value('x', x, i); ++i)
    {
      if ("quiet" == x)
      {
        report_multiple_hits_truncated = 0;
        if (verbose)
          cerr << "Will suppress warnings about too many substructure matches\n";
      }
      else if (! x.numeric_value(max_hits) || max_hits < 1)
      {
        cerr << "The maximum number of substructure matches (-x) must be a whole +ve number\n";
        usage(1);
      }
      else if (verbose)
        cerr << "A maximum of " << max_hits << " substructure matches will be used\n";
    }
  }

  if (cl.option_present('U'))
  {
    const_IWSubstring u = cl.string_value('U');

    if ("all" == u)
    {
      unique_across_all_molecules = 1;
    }
    else if ("each" == u)
    {
      unique_within_molecule = 1;
    }
    else
    {
      cerr << "Unrecognised -U qualifier '" << u << "'\n";
      usage(1);
    }
  }

  if (cl.option_present('y'))
  { 
    postprocess_molecules_produced = 1;

    if (verbose)
      cerr << "Will post process molecules produced just like input molecules\n";
  }

  if (cl.option_present('V'))
  {
    IWString v = cl.string_value('V');

    discard_bad_balences = 1;

    if ('.' != v)
    {
      if (! v.ends_with(".smi"))
        v << ".smi";

      if (! stream_for_bad_valence.open(v.null_terminated_chars()))
      {
        cerr << "Cannot open stream for bad valence '" << v << "'\n";
        return 1;
      }

      if (verbose)
        cerr << "MOlecules with bad valences discarded, written to '" << v << "'\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  const int nr = rxn.number_elements();

  int * molecules_hitting_reaction = new_int(nr); std::unique_ptr<int[]> free_molecules_hitting_reaction(molecules_hitting_reaction);

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! medchem_wizard(cl[i], input_type, rxn, molecules_hitting_reaction, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, produced " << molecules_produced << endl;
    if (truncated_to_max_hits > 0)
      cerr << truncated_to_max_hits << " substructure searches truncated to " << max_hits << endl;
    if (duplicate_molecules_suppressed)
      cerr << duplicate_molecules_suppressed << " duplicate molecules suppressed\n";

    for (int i = 0; i < nr; ++i)
    {
      float f = static_cast<float>(molecules_hitting_reaction[i]) / static_cast<float>(molecules_read);

      cerr << molecules_hitting_reaction[i] << " molecules hit " << rxn[i]->comment() << ' ' << f << endl;
    }
    cerr << bad_valences_discarded << " produces with bad valences discarded\n";
    if (discarded_for_too_many_atoms)
      cerr << discarded_for_too_many_atoms << " discarded_for_too_many_atoms\n";
    if (discarded_for_too_few_atoms)
      cerr << discarded_for_too_few_atoms << " discarded_for_too_few_atoms\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = medchem_wizard(argc, argv);

  return rc;
}
