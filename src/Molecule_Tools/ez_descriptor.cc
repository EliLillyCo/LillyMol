/*
  Produces descriptor(s) based on E/Z configurations
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::cout;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int ignore_molecules_not_matching = 0;

static int take_first_of_multiple_hits = 0;

static int skip_molecules_having_multiple_hits = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static angle_t tolerance = static_cast<angle_t> (10.0 * DEG2RAD);

static int use_coordinates = 0;

static int work_as_filter = 0;

int nq = 0;

static Molecule_Output_Object stream_for_matches[2];

static int write_descriptors = 1;

static IWString fingerprint_tag;

static int nbits = 0;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");

/*
  For each query, we need to keep track of how many E, Z and undetermined
  matches there are
*/

class EZ_Matches
{
  private:
    int _n;

    int _e;
    int _z;
    int _0;
  public:
    EZ_Matches();

    void extra (int);

    int report (std::ostream &) const;
};

EZ_Matches::EZ_Matches()
{
  _n = 0;

  _e = 0;
  _z = 0;
  _0 = 0;

  return;
}

void
EZ_Matches::extra (int s)
{
  if (s > 0)
    _z++;
  else if (s < 0)
    _e++;
  else
    _0++;

  _n++;

  return;
}

int 
EZ_Matches::report (std::ostream & output) const
{
  output << "EZ_Matches::report: data for " << _n << " molecules\n";
  output << _e << " E, " << _z << " Z, and " << _0 << " 0 matches\n";

  return output.good();
}

static EZ_Matches * ez_matches = nullptr;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id: ez_descriptor.cc,v 1.1.1.1 2005/05/10 17:38:28 rx87851 Exp $\n";
  cerr << "Produces descriptor(s) based on E/Z bonding configuration\n";
  cerr << "  -q <qry>      query file specification\n";
  cerr << "  -s <smarts>   give query as smarts\n";
  cerr << "  -z i          ignore molecules not matching a query\n";
  cerr << "  -z f          take the first of multiple query matches\n";
  cerr << "  -z nom        skip molecules having more than one hit to the query\n";
  cerr << "  -J <tag>      write results as a fingerprint with tag <tag>\n";
  cerr << "  -f            work as a TDT filter\n";
  cerr << "  -c            use the molecule's coordinates\n";
  cerr << "  -i ...        input specification(s)\n";
  cerr << "  -S <stem>     write matches to two files starting with <stem>\n";
  cerr << "  -o ...        output specifications\n";
  cerr << "  -A ...        aromaticity specifications, enter '-A help' for details\n";
  cerr << "  -E ...        element specifications, enter '-E help' for details\n";
  cerr << "  -v            verbose output\n";

  exit (rc);
}

static int
write_matches (Molecule & m,
               int ez)
{
  if (ez < 0)
    return stream_for_matches[0].write (m);
  else if (ez > 0)
    return stream_for_matches[1].write (m);

  return 0;
}

static int
do_write_descriptors (const IWString & mname,
                      const int * result_vector,
                      std::ostream & output)
{
  write_space_suppressed_string (mname, output);

  for (int i = 0; i < nq; i++)
  {
    output << ' ' << result_vector[i];
  }

  output << endl;

  return output.good();
}

static int
do_write_fingerprints (Molecule & m,
                       const int * result_vector,
                       std::ostream & output)
{
  if (! work_as_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  IW_Bits_Base fp (nbits + nbits);

  for (int i = 0; i < nq; i++)
  {
    if (result_vector[i] > 0)
      fp.set (2 * i);
    else if (result_vector[i] < 0)
      fp.set (2 * i + 1);
  }

  fp.write_daylight_ascii_representation (output, fingerprint_tag);

  if (! work_as_filter)
    output << "|\n";

  return output.good();
}

static int
ez_descriptor (const Molecule & m,
               atom_number_t a1,
               atom_number_t a2,
               atom_number_t a3,
               atom_number_t a4)
{
  const Bond * b23 = m.bond_between_atoms (a2, a3);
  if (! b23->is_double_bond())
  {
    cerr << "ez_descriptor:not a double bond, atoms " << a2 << " and " << a3 << " in '" << m.name() << "'\n";
    return 0;
  }

  if (! b23->part_of_cis_trans_grouping())
  {
    if (verbose)
      cerr << "ez_descriptor:not part of directional bond, atoms " << a2 << " and " << a3 << " in '" << m.name() << "'\n";
    return 0;
  }

  const Bond * b21 = m.bond_between_atoms (a1, a2);
  const Bond * b34 = m.bond_between_atoms (a3, a4);

  if (! b21->is_directional() || ! b34->is_directional())
  {
    cerr << "ez_descriptor:one or more bonds not directional '" << m.name() << "'\n";
    return 0;
  }

  int direction21;

  if (a2 == b21->a1())
  {
    if (b21->is_directional_up())
      direction21 = 1;
    else 
      direction21 = -1;
  }
  else
  {
    if (b21->is_directional_up())
      direction21 = -1;
    else 
      direction21 = 1;
  }

  int direction34;

  if (a3 == b34->a1())
  {
    if (b34->is_directional_up())
      direction34 = 1;
    else 
      direction34 = -1;
  }
  else
  {
    if (b34->is_directional_up())
      direction34 = -1;
    else 
      direction34 = 1;
  }

  if (direction21 == direction34)
    return 1;
  else 
    return -1;
}

static int
ez_descriptor (Molecule & m,
               const Set_of_Atoms & s,
               int & zresult)
{
  int n = s.number_elements();

  if (n < 4)
  {
    cerr << "Only " << n << " matched atoms in query match, '" << m.name() << "'\n";
    return 0;
  }

  atom_number_t s1 = s[1];
  atom_number_t s2 = s[2];

  const Bond * b = m.bond_between_atoms (s1, s2);
  if (! b->is_double_bond())
  {
    cerr << "Atoms " << s1 << " and " << s2 << " in '" << m.name() << "' not double bond\n";
    return 0;
  }

  if (b->is_cis_trans_either_double_bond())
  {
    if (verbose > 1)
      cerr << "Molecule '" << m.name() << " bond between atoms " << s1 << " and " << s2 << " is indeterminate\n";

    zresult = 0;

    return 1;
  }

  if (use_coordinates)
    zresult = m.ez_by_geometry (s[0], s1, s2, s[3], tolerance);
  else
    zresult = ez_descriptor (m, s[0], s1, s2, s[3]);

  return 1;
}

static int
ez_descriptor (Molecule & m,
               int * result_vector)
{
  Molecule_to_Match target (&m);

  int queries_matching = 0;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search (target, sresults);

    if (1 == nhits)     // great!
      ;
    else if (nhits > 1 && skip_molecules_having_multiple_hits)
      continue;
    else if (nhits > 1 && take_first_of_multiple_hits)
    {
      if (verbose)
        cerr << "Taking first of " << nhits << " query hits to '" << m.name() << "'\n";
    }
    else if (0 == nhits && ignore_molecules_not_matching)
      continue;
    else
    {
      cerr << "Fatal error, " << nhits << " hits to query " << i << " molecule '" << m.name() << "'\n";
      return 0;
    }

    queries_matching++;

    if (! ez_descriptor (m, *(sresults.embedding(0)), result_vector[i]))
      return 0;

    ez_matches[i].extra (result_vector[i]);
  }

  return 1;
}

/*
  Very crude test. If the first two atoms have identical coordinates,
  we assume no coordinates
*/

static int
determine_has_coordinates (Molecule & m)
{
  const Atom * a0 = m.atomi (0);
  const Atom * a1 = m.atomi (1);

  if (a0->x() != a1->x())
    return 1;

  if (a0->y() != a1->y())
    return 1;

  if (a0->z() != a1->z())
    return 1;

  return 0;
}

static int
ez_descriptor (Molecule & m,
               int * result_vector,
               std::ostream & output)
{
  const int has_coordinates = determine_has_coordinates(m);

  if (! has_coordinates && use_coordinates)
  {
    cerr << "Molecule '" << m.name() << "' has no coordinates, but request to use coordinates\n";
    return 0;
  }

  set_vector (result_vector, nq, 0);

  if (m.natoms() < 4)    // cannot have a cis-trans bond
    ;
  else if (0 == ez_descriptor (m, result_vector))
    return 0;

  if (stream_for_matches[0].active() && 0 != result_vector[0])
    write_matches (m, result_vector[0]);

  if (write_descriptors)
    do_write_descriptors (m.name(), result_vector, output);
  else if (nbits)
    do_write_fingerprints (m, result_vector, output);

  return output.good();
}

static int
ez_descriptor (data_source_and_type<Molecule> & input,
               int * result_vector,
               std::ostream & output)
{
  Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m (m);

    molecules_read++;

    if (! ez_descriptor (*m, result_vector, output))
      return 0;
  }

  return output.good();
}

static int
ez_descriptor (const const_IWSubstring & buffer,
               int * result_vector, 
               std::ostream & output)
{
  Molecule m;
  if (! m.build_from_smiles (buffer))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return ez_descriptor (m, result_vector, output);
}

static int
ez_descriptor (iwstring_data_source & input,
               int * result_vector,
               std::ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    output << buffer << endl;

    if (! buffer.starts_with (smiles_tag))
      continue;

    molecules_read++;

    buffer.remove_leading_chars (smiles_tag.length());
    assert (buffer.ends_with ('>'));
    buffer.chop();

    if (! ez_descriptor (buffer, result_vector, output))
    {
      cerr << "Fatal error processing line " << input.lines_read() << endl;
      return 0;
    }
  }

  return output.good();
}

static int
ez_descriptor (const char * fname,
               int * result_vector,
               std::ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ez_descriptor (input, result_vector, output);
}

static int
ez_descriptor (const char * fname,
               FileType input_type,
               int * result_vector,
               std::ostream & output)
{
  if (work_as_filter)
    return ez_descriptor (fname, result_vector, output);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name (fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ez_descriptor (input, result_vector, output);
}

static int
ez_descriptor (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vi:E:A:q:s:J:z:cfo:S:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! process_standard_aromaticity_options (cl, verbose > 1))
  {
    cerr << "Cannot process -A option\n";
    usage (11);
  }

  if (! process_elements (cl, verbose > 1, 'E'))
  {
    cerr << "Cannot initialise elements\n";
    usage (8);
  }

  if (cl.option_present ('c'))
  {
    use_coordinates = 1;

    if (verbose)
      cerr << "Will use coordinates to determine configuration\n";
  }

  if (cl.option_present ('f'))
  {
    work_as_filter = 1;

    if (! cl.option_present ('J'))
    {
      cerr << "To work as a filter, only fingerprints can be produced\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will work as a filter\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (work_as_filter)
    ;
  else if (cl.option_present ('i'))
  {
    if (! process_input_type (cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix (cl))
    return 4;

  if (cl.option_present ('s'))
  {
    int i = 0;
    const_IWSubstring s;
    while (cl.value ('s', s, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts (s))
      {
        cerr << "Invalid smarts '" << s << "'\n";
        return 6;
      }

      queries.add (q);
    }
  }

  if (cl.option_present ('q'))
  {
    if (! process_queries (cl, queries, verbose > 1, 'q'))
    {
      cerr << "Cannot process queries (-q option)\n";
      return 3;
    }
  }

  nq = queries.number_elements();

  if (verbose)
    cerr << "Defined " << nq << " queries\n";

  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_unique_embeddings_only (1);
  }

  ez_matches = new EZ_Matches[nq];

  if (cl.option_present ('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while (cl.value ('z', z, i++))
    {
      if ('i' == z)
      {
        ignore_molecules_not_matching = 1;
      }
      else if ('f' == z)
      {
        take_first_of_multiple_hits = 1;
      }
      else if ("nom" == z)
      {
        skip_molecules_having_multiple_hits = 1;
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage (5);
      }
    }
  }

  int size_of_result_vector = nq;

  if (cl.option_present ('J'))
  {
    cl.value ('J', fingerprint_tag);
    
    if (verbose)
      cerr << "Output written as fingerprints, tag '" << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag.add ('<');

    nbits = nq;
    if (nbits != nbits / 8 * 8)
      nbits = (nbits / 8 + 1) * 8;

    if (verbose)
      cerr << "Fingerprints will contain " << nbits << " bits\n";

    size_of_result_vector = nbits + nbits;

    write_descriptors = 0;
  }

  if (0 == size_of_result_vector)
  {
    cerr << "No queries defined, cannot continue\n";
    usage (9);
  }

  int * result_vector = new_int (size_of_result_vector); std::unique_ptr<int[]> free_result_vector (result_vector);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present ('S'))
  {
    if (nq > 1)
    {
      cerr << "Sorry, can only write match files with one query present\n";
      return 6;
    }

    const_IWSubstring s = cl.string_value ('S');
    for (int i = 0; i < 2; i++)
    {
      if (! cl.option_present ('o'))
        stream_for_matches[i].add_output_type (FILE_TYPE_SMI);
      else if (! stream_for_matches[i].determine_output_types (cl, 'o'))
      {
        cerr << "Cannot determine output type(s)\n";
        return 8;
      }

      IWString fname;
      fname << s << i;
      if (stream_for_matches[i].would_overwrite_input_files (cl, fname))
      {
        cerr << "The output file specification (-S) cannot overwrite input file(s) '" << fname << "'\n";
        return 6;
      }

      if (! stream_for_matches[i].new_stem (fname))
      {
        cerr << "Cannot initialise output file '" << fname << "'\n";
        return 8;
      }

      if (verbose)
        cerr << "Matches written to '" << fname << "'\n";
    }
  }

  if (0 == fingerprint_tag.length())
  {
    cout << "Name";
    for (int i = 0; i < nq; i++)
    {
      const Substructure_Hit_Statistics * q = queries[i];
      if (0 == q->comment().length())
        cout << " EZ" << i;
      else
        cout << ' ' << q->comment();
    }

    cout << endl;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ez_descriptor (cl[i], input_type, result_vector, cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < nq; i++)
    {
      ez_matches[i].report (cerr);
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ez_descriptor (argc, argv);

  return rc;
}
