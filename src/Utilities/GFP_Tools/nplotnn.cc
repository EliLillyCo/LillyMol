/*
  Compiled version of plotnn
*/

#include <stdlib.h>
#include <limits>

#include <fstream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iw_stl_hash_set.h"
#include "iw_stl_hash_map.h"
#include "cmdline.h"
#include "iwrandom.h"
#include "iwstring_data_source.h"
#include "accumulator.h"
#include "iwhistogram.h"
#include "iwdigits.h"

#include "smiles_id_dist.h"
#include "distance_scaling.h"

#include "mol/molecule.h"
#include "mol/smiles.h"

static IWString identifier_tag("PCN<");

static IWString smiles_tag("$SMI<");

static IWString distance_tag("DIST<");

static IWString scale_tag;

static int output_precision = 2;

static int from_gfp_leader = 0;
static int clusters_found = 0;

static int tabular_leader_output = 0;

static int molecules_written = 0;

static IWString append_to_target_record;

static int append_neighbour_number_to_each_neighbour = 1;

// Sometimes people want things printed all on one line

static IWString neighbour_separator('\n');

static IWString space_or_tab(' ');

static IWString cluster(" CLUSTER ");

static IWDigits iwdigits;   // used for neighbour numbers
static Fraction_as_String fraction_as_string;

/*
  We may want to only write (say) the first token in the names
*/

static int target_column_to_write = -1;
static resizable_array<int> neighbour_columns_to_write;

static int mandatory_neighbour_count = -1;

#define IW_FLUSH_BUFFER 8192

static IWString_and_File_Descriptor stream_for_targets_not_written;

/*
  Sometimes we just want a table of near neighbour distances
*/

static int tabular_output = 0;
static int tabular_output_nearest_neighbour_only = 0;

/*
  We may get jobs from programmes that don't have smiles. We can
  fill in missing smiles
*/

static IW_STL_Hash_Map_String missing_smiles;

/*
  Sept 2005. In doing Tversky searches on fragments, I want to be able to
  specify a minimum and maximum number of atoms difference between the
  needle and the members of the haystack
*/

static int min_extra_atoms = 0;
static int max_extra_atoms = 0;

static int initial_nbr_list_size = 100;

static int take_first_token_of_name_field = 0;

static Distance_Scaling distance_scaling;

/*
  Jarvis Patrick does not have any distances
*/

static int distances_present_in_input = 1;

static char three_column_output_separator = ' ';
static int three_column_output = 0;

/*
  In addition to the smiles, the ID and the distance, we need a unique identifier for each neighbour.
  It may be just the ID, or it may be the unique smiles of the molecule
*/

class Smiles_ID_Dist_UID : public Smiles_ID_Dist
{
  private:
    IWString _uid;

  public:
    Smiles_ID_Dist_UID (const IWString & s, const IWString & i, similarity_type_t d) : Smiles_ID_Dist (s, i, d) {};

    void set_unique_identifier (const IWString & s) { _uid = s;}

    const IWString & unique_identifier() const { return _uid;}
};

#ifdef __GNUG__
template class resizable_array_p<Smiles_ID_Dist_UID>;
template class resizable_array_base<Smiles_ID_Dist_UID *>;
#endif

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << endl;

  cerr << "Processes the output from gfp_nearneighbours and produces a smiles file\n";
  cerr << " -n <number>      number of neighbours per structure to display\n";
  cerr << " -t <dist>        discard neighbours with distance <dist> or less\n";
  cerr << " -T <dist>        discard neighbours with distance <dist> or greater\n";
  cerr << " -f <number>      number of structures per page in vf\n";
//cerr << " -r               choose neighbours at random (requires -n option)\n";
//cerr << "                  favours close neighbours\n";
  cerr << " -x               don't print target molecules\n";
  cerr << " -z               don't write molecules with no neighbours\n";
  cerr << " -W N,D           don't write a molecule unless it has N nbrs within distance D\n";
  cerr << " -p               if necessary, pad each neighbour list to -n length\n";
  cerr << " -D <fname>       distance scaling data file\n";
  cerr << " -u ...           options for uniqueness among neighbours, enter '-u help'\n";
  cerr << " -Y <string>      append <string> to each target record\n";
  cerr << " -O <number>      write distances as '<number> - distance' - useful for writing\n";
  cerr << "                  similarities rather than distances\n";
  cerr << " -s               collect statistics on the neighbour distances\n";
  cerr << " -c <precision>   output precision for distances (default " << output_precision << " decimal places)\n";
  cerr << " -L def           output is from gfp_leader_v2, use '-L tbl' for tabular output\n";
  cerr << " -j <string>      separator between records (default newline) - use ' ' for tabular output. 'space' and 'tab' OK\n";
  cerr << " -w               remove the constraint that distances must be in [0,1]\n";
  cerr << " -X ...           more obscure options, enter '-X help' for info\n";
  cerr << " -3               produce three column output 'id1 id2 dist'\n";
  cerr << " -H <fname>       write the nearest neighbour histogram to <fname>\n";
  cerr << " -h               discard neighbours with zero distance and same ID as target\n";
  cerr << " -M <fname>       smiles file to be used for filling missing smiles\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static void
display_dash_x_options (ostream & os)
{
  os << " -X nonnum          suppress printing of neighbour number with each neighbour\n";
  os << " -X tcol=<col>      just write column <col> from the target molecules\n";
  os << " -X ncol=<col>      just write column <col> from the neighbour molecules\n";
  os << " -X ftn             take the first token of both target and neighbour molecules\n";
  os << " -X TABS            create tab separated output rather than space separated\n";
  os << " -X minextra=n      neighbours must have at least     <n> atoms more than target\n";
  os << " -X maxextra=n      neighbours must have no more than <n> atoms more than target\n";
  os << " -X nosmi           input does not contain smiles\n";
  os << " -X random          random subset of neighbours (requires -n option)\n";
  os << " -X brs             biased random subset of neighbours (requires -n option)\n";
  os << " -X SMITAG=<tag>    Set smiles     tag to <tag>, default " << smiles_tag << '\n';
  os << " -X IDTAG=<tag>     Set identifier tag to <tag>, default " << identifier_tag << '\n';
  os << " -X DISTAG=<tag>    Set distance   tag to <tag>, default " << distance_tag << '\n';
  os << " -X NOTW=<fname>    write targets that don't get written to <fname>\n";
  os << " -X SCALE=<tag>     scaled spread output, unscale the distances\n";
  os << " -X table           produce tabular output of near neighbour distances\n";
  os << " -X table1          produce tabular output of nearest neighbour distances\n";
  os << " -X normh           normalise the -H file to largest count\n";
  os << " -X allh            include all distances in the -H file (not just nearest)\n";
  os << " -X HR=<fname>      create an R file for plotting the -H file\n";
  os << " -X cumh            the data in the HR file is cumulative\n";
  os << " -X ckds            check that all incoming distances are sorted\n";
  os << " -X 0               ignore leading 0's when comparing identifiers\n";
  os << " -X sdcensor=<d>    in tabular output, censor distances shorter than <d>\n";
  os << " -X mnc=<n>         only write targets that have exactly<n> nbrs after distance filtering\n";
  os << " -X sdo             with the -s option, only write the shortest distance\n";

  return;
}

static void
display_dash_u_options (ostream & os)
{
  cerr << " -u id            unique neighbours only - note that only identifiers are\n";
  cerr << "                  checked, NOT the unique smiles\n";
  cerr << " -u smi           unique neighbours only - uniqueness by unique smiles\n";
  cerr << " -u nochiral      unique neighbours only - uniqueness by unique smiles without chirality\n";
  cerr << " -u APP=XXX       write duplicate neighbours but append XXX to each\n";

  return;
}

class Fatal_Error
{
  private:
    IWString _fname;
    int _line_number;

  public:
    Fatal_Error (const char *, int);

    int line_number() const { return _line_number;}
    const IWString & fname() const { return _fname;}
};

Fatal_Error::Fatal_Error (const char * file_name, int n) : _fname(file_name), _line_number(n)
{
}

ostream & 
operator << (ostream & os, const Fatal_Error & f)
{
  os << "Fatal error at line " << f.line_number() << " in '" << f.fname() << "'\n";

  return os;
}

static int verbose = 0;

static int molecules_processed = 0;

static int write_molecules_with_no_neighbours = 1;

static int print_target_molecules = 1;

static int neighbours_per_structure = 0;

static int suppress_neighbours = 0;

static int min_neighbours_per_structure = 0;

static int vf_per_page = 0;

static int pad_neighbour_list = 0;

static int choose_neighbours_at_random = 0;

static int biased_subset = 0;

static Accumulator_Int<int> neighbour_statistics;

static Accumulator<float> distance_stats;

static Accumulator<float> nearest_neighbour_stats;
static Accumulator<float> furthest_neighbour_stats;

static IWHistogram nearest_neighbour_histogram;

static int place_all_distances_in_nn_histogram = 0;

static int discard_self_neighbours = 0;

static int zero_distance_neighbours = 0;

static extending_resizable_array<int> neighbour_count;

static float lower_distance_threshold = -1.0;

static float upper_distance_threshold = -1.0;

/*
  Sept 04. I want to be able to discard molecules that don't have at least N neighbours within
  a given distance
*/

static int number_needed_within_distance = -1;
static float distance_for_number_needed_within_distance = static_cast<float>(0.0);

/*
  When we get files from tools other than gfp, we may get distances beyond the range [0,1]
*/

static float max_possible_distance = 1.0;

#define UNIQUE_BY_ID 1
#define UNIQUE_BY_SMILES 2

static int unique_neighbours_only = 0;

static IWString append_to_non_unique_neighbours;

static int collect_statistics = 0;

static int write_shortest_distance_with_stats = 0;

/*
  The identifiers encountered may be either text identifiers or unique smiles
*/

static IWString_STL_Hash_Set identifiers_encountered;

/*
  We keep a count of the number of molecules suppressed because they were duplicates
*/

static int duplicate_neighbours_suppressed = 0;

static int ignore_leading_zeros_in_identifiers = 0;

/*
  Some people want to see distances written as similarities
*/

static float offset = -1.0;

/*
  Nov 2010. We are returning distances to a collaborator and want
  to censor all distances below a cutoff
*/

static float censor_distances_shorter_than = 0.0f;

static int check_distance_ordering = 0;

/*
  Because the neighbour list is sorted, this test is pretty easy
*/

static int
passes_number_needed_within_distance (const resizable_array_p<Smiles_ID_Dist_UID> & neighbours)
{
  int nn = neighbours.number_elements();

  if (0 == number_needed_within_distance)    // we want zero neighbours within the given distance
  {
    if (0 == nn)    // can this happen?
      return 1;

    if (neighbours[0]->distance() < distance_for_number_needed_within_distance)   // if the first nbr is too close, then we are done
      return 0;

    return 1;
  }

  if (nn < number_needed_within_distance)
    return 0;

  return neighbours[number_needed_within_distance - 1]->distance() <= distance_for_number_needed_within_distance;
}

static int
write_neighbour_list (const resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                      const IWString & target_smiles,
                      IWString_and_File_Descriptor & output)
{
  if (suppress_neighbours)
    return 1;

  int n = neighbours.number_elements();

  int lines_written_this_molecule;
  if (print_target_molecules)
    lines_written_this_molecule = 1;
  else
    lines_written_this_molecule = 0;

  for (int i = 0; i < n; i++)
  {
    if (i > 0 && vf_per_page > 0 && 0 == i % (vf_per_page - 1))
    {
      lines_written_this_molecule++;
      output << target_smiles << neighbour_separator;
    }

    const Smiles_ID_Dist_UID * sid = neighbours[i];

    if (smiles_tag.length())
      output << sid->smiles() << space_or_tab;

    output << sid->id() << space_or_tab;
    if (append_neighbour_number_to_each_neighbour)
      output << (i + 1) << space_or_tab;

    if (distances_present_in_input)
    {
      float d = sid->distance();

      if (static_cast<float>(0.0) == d)
        zero_distance_neighbours++;

      if (offset > static_cast<float>(0.0))
        d = (offset - sid->distance());

      if (fraction_as_string.active())
        fraction_as_string.append_number(output, d);
      else
      {
        output.append_number(d, output_precision);

        output += neighbour_separator;
      }
    }
    else
      output += neighbour_separator;

    molecules_written++;

    lines_written_this_molecule++;

    if (unique_neighbours_only)
      identifiers_encountered.insert(sid->unique_identifier());

    output.write_if_buffer_holds_more_than(32768);
  }

//cerr << "Printed " << n << " neighbours, pad = " << pad_neighbour_list << " per = " << neighbours_per_structure << endl;

  if (pad_neighbour_list && n < neighbours_per_structure)
  {
    for (int i = n; i < neighbours_per_structure; i++)
    {
      output << '*' << neighbour_separator;
    }
  }

  if (vf_per_page > 0 && 0 != lines_written_this_molecule % vf_per_page)
  {
    while (0 != lines_written_this_molecule % vf_per_page)
    {
      output << '*' << neighbour_separator;
      lines_written_this_molecule++;
    }
  }

  return 1;
}

/*
  Append statistics about its neighbours to the ID of a target molecule
*/

template <typename T>
void
write_statistics_for_neighbour_list (const resizable_array_p<T> & neighbours,
                                     IWString & output)
{
  int n = neighbours.number_elements();

  if (tabular_output_nearest_neighbour_only)
  {
    if (0 == n)
    {
      output << " 1";
      return;
    }

    float d = neighbours[0]->distance();
    if (d < censor_distances_shorter_than)
      output << " 0";
    else
    {
      output << ' ';
      output.append_number(d, output_precision);
    }

    return;
  }

  if (write_shortest_distance_with_stats)
  {
    if (0 == n)
    {
      output << " 1";
      return;
    }

    if (n > 1)
      output << " N=" << n;

    output << ' ';

    output.append_number(neighbours[0]->distance(), output_precision);

    return;
  }

  Accumulator<float> d;

  for (int i = 0; i < n; i++)
  {
    const T * sid = neighbours[i];

    d.extra(sid->distance());
  }

  if (tabular_output)
  {
    if (from_gfp_leader)
      output << ' ' << (n + 1) << ' ';
    else
      output << ' ' << n << ' ';
  }
  else
    output << " N=" << n << ' ';

  output.append_number(d.minval(), output_precision);
  output << ' ';
  output.append_number(d.maxval(), output_precision);
  output << ' ';
  if (d.n() > 0)
    output.append_number(static_cast<float>(d.average()), output_precision);
  else
    output << '.';

  return;
}

#ifdef __GNUG__
template void write_statistics_for_neighbour_list (const resizable_array_p<Smiles_ID_Dist_UID> &, IWString &);
#endif

static void
do_biased_subset (resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                  int nkeep)
{
  while (neighbours.number_elements() > nkeep)
  {
    random_number_t r = iwrandom();
    r = r * r;

    int i = int(r * (neighbours.number_elements() - 1));
    if (0 == i || i == neighbours.number_elements() - 1)
      continue;

    neighbours.remove_item(i);
  }
}

template <typename T>
void
random_subset (resizable_array_p<T> & neighbours,
               int nkeep)
 
{
  assert (nkeep > 0);

  while (neighbours.number_elements() > nkeep)
  {
    int i = intbtwij(0, neighbours.number_elements() - 1);

    neighbours.remove_item(i);
  }

  return;
}

static int
write_targets_not_written_if_requested (const IWString & smiles,
                                        const IWString & id,
                                        IWString_and_File_Descriptor & output)
{
  if (! output.is_open())
    return 1;

  output << smiles << ' ' << id << '\n';

  output.write_if_buffer_holds_more_than(IW_FLUSH_BUFFER);

  return 1;
}

static int
do_check_distance_ordering (const resizable_array_p<Smiles_ID_Dist_UID> & neighbours)
{
  const int n = neighbours.number_elements();

  if (0 == n)
    return 1;

  float dprev = neighbours[0]->distance();

  int rc = 1;

  for (int i = 1; i < n; ++i)
  {
    const float d = neighbours[i]->distance();

    if (d < dprev)
    {
      cerr << "do_check_distance_ordering:out of order neighbours\n";
      cerr << "parent " << neighbours[0]->id() << " nbr " << i << " dprev " << dprev << " d " << d << endl;
      rc = 0;
    }

    dprev = d;
  }

  return rc;
}

#ifdef __GNUG__
template void random_subset (resizable_array_p<Smiles_ID_Dist_UID> &, int);
#endif

static int
do_three_column_output (const IWString & id,
                        const resizable_array_p<Smiles_ID_Dist_UID> & neighbours, 
                        IWString_and_File_Descriptor & output)
{
  const int nn = neighbours.number_elements();

  if (0 == nn)
    return 1;

  nearest_neighbour_stats.extra(neighbours[0]->distance());
  furthest_neighbour_stats.extra(neighbours.last_item()->distance());

  const_IWSubstring first_token_id(id);
  if (first_token_id.nwords() > 1)
    first_token_id.truncate_at_first(' ');

  for (int i = 0; i < nn; ++i)
  {
    output << first_token_id << three_column_output_separator << neighbours[i]->id();

    if (fraction_as_string.active())
      fraction_as_string.append_number(output, neighbours[i]->distance());
    else
      output << three_column_output_separator << neighbours[i]->distance();

    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
    if (0.0f == neighbours[i]->distance())
      zero_distance_neighbours++;
  }

  molecules_written += nn;

  return output.good();
}

static int
process_molecule (const IWString & smiles,
                  const IWString & id,
                  const resizable_array_p<Smiles_ID_Dist_UID> & neighbours, 
                  IWString_and_File_Descriptor & output)
{
  molecules_processed++;

  int nn = neighbours.number_elements();

  neighbour_statistics.extra(nn);
  neighbour_count[nn]++;

  if (verbose > 1)
    cerr << id << " has " << nn << " neighbours\n";

  if (three_column_output)
    return do_three_column_output(id, neighbours, output);

  if (number_needed_within_distance >= 0 && ! passes_number_needed_within_distance(neighbours))
    return write_targets_not_written_if_requested(smiles, id, stream_for_targets_not_written);

  if (mandatory_neighbour_count >= 0 && nn != mandatory_neighbour_count)
    return 1;

  if (0 == nn && ! write_molecules_with_no_neighbours)
    return write_targets_not_written_if_requested(smiles, id, stream_for_targets_not_written);

  if (tabular_output || tabular_output_nearest_neighbour_only)
  {
    append_first_token_of_name(id, output);

    if (collect_statistics)
      write_statistics_for_neighbour_list(neighbours, output);

    output << '\n';

    molecules_written++;

    return 1;
  }

  if (print_target_molecules)
  {
    if (smiles_tag.length())
      output << smiles << space_or_tab;

    output << id;

    if (collect_statistics)
      write_statistics_for_neighbour_list(neighbours, output);

    if (append_to_target_record.length())
      output << space_or_tab << append_to_target_record;

    output << neighbour_separator;

    molecules_written++;
  }

  if (nn > 0)
  {
    nearest_neighbour_stats.extra(neighbours[0]->distance());

    if (nearest_neighbour_histogram.active())
    {
      if (place_all_distances_in_nn_histogram)
      {
        for (int i = 0; i < nn; i++)
        {
          nearest_neighbour_histogram.extra(neighbours[i]->distance());
        }
      }
      else
        nearest_neighbour_histogram.extra(neighbours[0]->distance());
    }

    furthest_neighbour_stats.extra(neighbours.last_item()->distance());
  }

  int rc = write_neighbour_list(neighbours, smiles, output);

//cerr << "After write_neighbour_list '" << output << "'\n";

  if (! neighbour_separator.contains('\n'))
    output << '\n';

  return rc;
}

/*
  Probably need to re-do this whole thing sometime.
  The identifiers_encountered set was set up to deal with
  neighbours that show up as neighbours of different
  targets. But, unique smiles may show up multiple times
  as neighbours of one molecule.
*/

static int
same_structure_this_target (const resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                            const IWString & usmi)
{
  int n = neighbours.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (usmi == neighbours[i]->unique_identifier())
      return 1;
  }

  return 0;
}

/*
  Check whether or not an identifier is a duplicate. Update the hash and
  global counters as needed
*/

static int
is_duplicate (const IWString_STL_Hash_Set & identifiers_encountered,
              const IWString & id)
{
//cerr << "Checking duplicate? '" << id << "'\n";
  if (identifiers_encountered.contains(id))
  {
    duplicate_neighbours_suppressed++;
    return 1;
  }

//identifiers_encountered.insert(id);    not correct, perhaps this neighbour won't be printed

  return 0;
}

static int
identifiers_the_same (const IWString & id1,
                      const IWString & id2)
{
  if (id1 == id2)
    return 1;

  const_IWSubstring s1(id1);
  const_IWSubstring s2(id2);

  if (ignore_leading_zeros_in_identifiers)
  {
    s1.remove_leading_chars('0');
    s2.remove_leading_chars('0');
  }

  if (s1 == s2)
    return 1;

  s1.truncate_at_first(' ');    // this is kind of dangerous, but will be OK almost all the time
  s2.truncate_at_first(' ');

//cerr << "Comparing truncated forms '" << s1 << "' and '" << s2 << "'\n";

  return s1 == s2;
}

static int
neighbour_is_suppressed (const IWString & id_of_target,
                         const IWString & id_of_nbr,
                         resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                         float distance)
{
//cerr << "neighbour_is_suppressed " << discard_self_neighbours << " cmp " << id_of_target << "' and '" << id_of_nbr << "'\n";

  if (0 == discard_self_neighbours)
    ;
  else if (1 == discard_self_neighbours)
  {
    if (static_cast<float>(0.0) == distance && identifiers_the_same(id_of_target, id_of_nbr))
      return 1;
  }
  else if (2 == discard_self_neighbours && identifiers_the_same(id_of_target, id_of_nbr))
    return 1;

  if (neighbours.number_elements() < min_neighbours_per_structure)  // must write the molecule
    return 0;

  if (lower_distance_threshold >= 0.0 && distance <= lower_distance_threshold)
    return 1;

  if (upper_distance_threshold > 0.0 && distance >= upper_distance_threshold)
    return 1;

  return 0;
}

static void
create_neighbour_item (const IWString & id_of_target,
                       resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                       const IWString & smiles,
                       const IWString & id,
                       float distance,
                       float scale)
{
  if (neighbour_is_suppressed(id_of_target, id, neighbours, distance))
    return;

  IWString usmi;      // may not get set, but needs to be scoped here

  int neighbour_is_duplicate = 0;

  if (0 == unique_neighbours_only)
    ;
  else if (UNIQUE_BY_ID == unique_neighbours_only)
  {
    neighbour_is_duplicate = is_duplicate(identifiers_encountered, id);
  }
  else if (UNIQUE_BY_SMILES == unique_neighbours_only)
  {
    Molecule m;
    if (! m.build_from_smiles(smiles))
    {
      cerr << "Yipes, cannot parse smiles '" << smiles << "'\n";
      throw Fatal_Error (__FILE__, __LINE__);
    }

    usmi = m.unique_smiles();

    neighbour_is_duplicate = is_duplicate(identifiers_encountered, usmi);
    if (! neighbour_is_duplicate)
      neighbour_is_duplicate = same_structure_this_target(neighbours, usmi);
  }

  IWString myid;
  if (! neighbour_is_duplicate)   // we will write it
    myid = id;
  else if (append_to_non_unique_neighbours.length())
  {
    myid = id;
    myid.append_with_spacer(append_to_non_unique_neighbours);
  }
  else
    return;

  distance = distance / scale;

  distance_stats.extra(distance);

  Smiles_ID_Dist_UID * sid = new Smiles_ID_Dist_UID(smiles, myid, distance);

// Should always use the ID first even if the ultimate comparison is via unique smiles.
// fix sometime...

  if (unique_neighbours_only)
  {
  }

  if (0 == unique_neighbours_only)
    ;
  else if (UNIQUE_BY_ID == unique_neighbours_only)
    sid->set_unique_identifier(id);
  else if (UNIQUE_BY_SMILES == unique_neighbours_only)
    sid->set_unique_identifier(usmi);

  neighbours.add(sid);

  return;
}

static int
extract_tdt_value (const const_IWSubstring & buffer,
                   int tag_length,
                   IWString & result)
{
  if (! buffer.ends_with('>'))
  {
    cerr << "TDT items must end in > '" << buffer << "'\n";
    throw Fatal_Error(__FILE__, __LINE__);
  }

  int bstop = buffer.length() - 2;
  if (bstop < tag_length)    // happens with 'PCN<>'
  {
    result = "";
    return 1;
  }

  buffer.from_to(tag_length, buffer.length() - 2, result);

  return 1;
}

static int
process_neighbour_list_item (const const_IWSubstring & buffer,
                             const IWString & tag,
                             IWString & result)
{
  if (result.length())
  {
    cerr << "Consecutive '" << tag << "' dataitems in neighbour list. Impossible\n";
    return 0;
  }

  if (! extract_tdt_value(buffer, tag.length(), result))
  {
    cerr << "Invalid neighbour list record '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
fetch_missing_smiles (const IWString & id,
                      IWString & smiles)
{
  IW_STL_Hash_Map_String::const_iterator f = missing_smiles.find(id);

//cerr << "Fetching smiles for '" << id << "' " << (f == missing_smiles.end()) << endl;

  if (f != missing_smiles.end())
  {
    smiles = (*f).second;
    return 1;
  }

  if (1 == id.nwords())
    return 0;

  const_IWSubstring tmp(id);
  tmp.truncate_at_first(' ');

  f = missing_smiles.find(tmp);

  if (f != missing_smiles.end())
  {
    smiles = (*f).second;
    return 1;
  }

  return 0;
}

static int
reduce_to_token (IWString & s, 
                 int w)
{
  if (1 == s.nwords())
    return 0;

  if (0 == w)
  {
    s.truncate_at_first(' ');
    return 1;
  }

  s.remove_leading_words(w);
  s.truncate_at_first(' ');

  return 1;
}

static int
reduce_to_token (IWString & s,
                 const resizable_array<int> & c)
{
  IWString rc;

  int i = 0;
  const_IWSubstring token;

  for (int col = 0; s.nextword(token, i); ++col)
  {
    if (! c.contains(col))
      continue;

    rc.append_with_spacer(token);
  }

  s = rc;

  return 1;
}

/*
  If smiles are missing from the input file, we just get the next ID
*/

static int
get_next_id (iwstring_data_source & input,
             IWString & id,
             int & fatal)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    if (input.eof())
      fatal = 0;
    else
      fatal = 1;
    return 0;
  }

  if (! buffer.starts_with(identifier_tag))
  {
    cerr << "First token in non-smiles not '" << identifier_tag << "'\n";
    fatal = 1;
    return 0;
  }

  if (! extract_tdt_value(buffer, identifier_tag.length(), id))
  {
    cerr << "Cannot extract identifier from '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (0 == id.length())
  {
    cerr << "Blank ID on line " << input.lines_read() << endl;
    id << "line " << input.lines_read();
    fatal = 1;
    return 0;
  }

  return 1;
}

/*
  Each TDT grouping starts with the smiles and the ID of the target molecule
*/

static int
get_smiles_and_id (iwstring_data_source & input,
                   IWString & smiles,
                   IWString & id,
                   int & fatal)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with(smiles_tag))
    {
      if (! extract_tdt_value(buffer, smiles_tag.length(), smiles))
      {
        cerr << "Cannot extract smiles from '" << buffer << "'\n";
        return 0;
      }

      if (smiles.length() && id.length())
        return 1;
    }
    else if (buffer.starts_with(identifier_tag))
    {
      if (id.length())
      {
        cerr << "Duplicate '" << identifier_tag << "' items, line " << input.lines_read() << endl;
        return 0;
      }

      if (! extract_tdt_value(buffer, identifier_tag.length(), id))
      {
        cerr << "Cannot extract identifier from '" << buffer << "'\n";
        return 0;
      }

      if (0 == id.length())
      {
        cerr << "Blank ID on line " << input.lines_read() << endl;
        id << "line " << input.lines_read();
      }

      if (smiles.length() && id.length())
        return 1;

      if (missing_smiles.size())
        return fetch_missing_smiles(id, smiles);
    }
    else if ('|' == buffer)
      return 0;
  }

// Should not have just one of SMILES or ID

  if (smiles.length() || id.length())
  {
    fatal = 1;
    return 0;
  }

  if (input.eof())
  {
    fatal = 0;
    return 0;
  }

  fatal = 1;

  return 0;
}

// Handles both SCALE and DIST tags.

static int
process_neighbour_list_item (const const_IWSubstring & buffer,
                             const IWString & tag,
                             float & dist)
{
  IWString tmp;

  if (! extract_tdt_value(buffer, tag.length(), tmp))
  {
    cerr << "Invalid neighbour list record '" << buffer << "'\n";
    return 0;
  }

  if (! tmp.numeric_value(dist) || dist < 0.0 || dist > max_possible_distance)
  {
    cerr << "Invalid distance '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
process_neighbour_list_record (const IWString & id_of_target,
                               resizable_array_p<Smiles_ID_Dist_UID> & neighbours,
                               IWString & smiles,
                               IWString & id,
                               float & distance,
                               const const_IWSubstring & buffer,
                               float & scale,
                               int record_number)
{
//#define DEBUG_NLR
#ifdef DEBUG_NLR
  cerr << "Processing '" << buffer << "'\n";
  cerr << "smiles '" << smiles << "' id '" << id << "' dist " << distance << endl;
#endif

  if (smiles_tag.length() && buffer.starts_with(smiles_tag))
  {
    if (! process_neighbour_list_item(buffer, smiles_tag, smiles))
      return 0;
  }
  else if (buffer.starts_with(identifier_tag))
  {
    if (! process_neighbour_list_item(buffer, identifier_tag, id))
      return 0;

    if (0 == id.length())
    {
      cerr << "Blank ID on line " << record_number << endl;
      id << "HUH, line " << record_number;
    }

    if (take_first_token_of_name_field)
      id.truncate_at_first(' ');

    if (0 == smiles_tag.length())
      ;
    else if (0 == smiles.length() && missing_smiles.size() > 0)
      fetch_missing_smiles(id, smiles);
  }
  else if (buffer.starts_with(distance_tag))
  {
    if (! process_neighbour_list_item(buffer, distance_tag, distance))
      return 0;
  }
  else if (scale_tag.length() && buffer.starts_with(scale_tag))
  {
    if (! process_neighbour_list_item(buffer, scale_tag, scale))
      return 0;
  }
  else
    return 1;

  if (0 == id.length())
    return 1;
  if (distances_present_in_input && distance < static_cast<float>(0.0))
    return 1;
  if (smiles_tag.length() && 0 == smiles.length())
    return 1;

  if (! distances_present_in_input)
    distance = 0.0f;

// We have all 3 items. Create a new neighbour item if the distance constraints are satisfied

  if (neighbour_columns_to_write.size() > 0)
    reduce_to_token(id, neighbour_columns_to_write);

// Some processing needs a Molecule to be built from the smiles

//Molecule * nbr_molecule = NULL;

  if (distance_scaling.active())
    distance = distance_scaling.convert(distance);

  create_neighbour_item(id_of_target, neighbours, smiles, id, distance, scale);

  smiles.resize_keep_storage(0);
  id.resize_keep_storage(0);
  distance = static_cast<float>(-1.0);
  scale = static_cast<float>(1.0);

  return 1;
}

static int
process_molecule (iwstring_data_source & input,
                  int & fatal,
                  IWString_and_File_Descriptor & output)
{
  IWString smiles, id;

  if (smiles_tag.length())
  {
     if (! get_smiles_and_id(input, smiles, id, fatal))
      return 0;
  }
  else if (! get_next_id(input, id, fatal))
    return 0;

  if (target_column_to_write >= 0)
    reduce_to_token(id, target_column_to_write);
  else if (take_first_token_of_name_field)
    id.truncate_at_first(' ');

  resizable_array_p<Smiles_ID_Dist_UID> neighbours;

  neighbours.resize(initial_nbr_list_size);

  const_IWSubstring buffer;

  IWString nsmiles, nid;     // smiles and ID of neighbours
  float distance = -1.0;
  float scale = 1.0F;

  while (input.next_record(buffer))
  {
    if ('|' == buffer)
    {
      if (0 == nsmiles.length() && 0 == nid.length() && distance >= 0.0f)   // case of just a distance, but no nbrs - like from gfp_spread_standard
        create_neighbour_item(id, neighbours, smiles, id, distance, scale);
        
      break;
    }

    if (! process_neighbour_list_record(id, neighbours, nsmiles, nid, distance, buffer, scale, input.lines_read()))
    {
      cerr << "Invalid neighbour list record, line '" << input.lines_read() << "'\n";
      cerr << buffer << endl;
      return 0;
    }
  }

// Wow, all kinds of interesting things about which of these operations should be done
// first. Do thresholding first and then resize, or should it be done the other way round?

//if (lower_distance_threshold >= 0.0 || upper_distance_threshold >= 0.0)
//  filter_to_thresholds (neighbours);

  if (neighbours_per_structure > 0 && neighbours_per_structure < neighbours.number_elements())
  {
    if (choose_neighbours_at_random)
      random_subset(neighbours, neighbours_per_structure);
    else if (biased_subset)
      do_biased_subset(neighbours, neighbours_per_structure);
    else
      neighbours.resize(neighbours_per_structure);
  }

  if (check_distance_ordering)
    do_check_distance_ordering(neighbours);

  return process_molecule(smiles, id, neighbours, output);
}

static int
process_molecule_gfp_leader (iwstring_data_source & input,
                             int & fatal,
                             IWString_and_File_Descriptor & output)
{
  IWString smiles, id;

  if (! get_smiles_and_id(input, smiles, id, fatal))
    return 0;

  if (target_column_to_write >= 0)
    reduce_to_token(id, target_column_to_write);
  else if (take_first_token_of_name_field)
    id.truncate_at_first(' ');

  resizable_array_p<Smiles_ID_Dist_UID> neighbours;

  const_IWSubstring buffer;

  IWString nsmiles, nid;     // smiles and ID of neighbours
  float distance = -1.0F;
  float scale = 1.0F;

  while (input.next_record(buffer))
  {
    if ('|' == buffer)
      break;

    if (! process_neighbour_list_record(id, neighbours, nsmiles, nid, distance, buffer, scale, input.lines_read()))
    {
      cerr << "Invalid neighbour list record, line '" << input.lines_read() << "'\n";
      cerr << buffer << endl;
      return 0;
    }
  }

// Wow, all kinds of interesting things about which of these operations should be done
// first. Do thresholding first and then resize, or should it be done the other way round?

//if (lower_distance_threshold >= 0.0 || upper_distance_threshold >= 0.0)
//  filter_to_thresholds (neighbours);

  int neighbours_initially_found = neighbours.number_elements();

  if (neighbours_per_structure > 0 && neighbours_per_structure < neighbours.number_elements())
  {
    if (choose_neighbours_at_random)
      random_subset(neighbours, neighbours_per_structure);
    else if (biased_subset)
      do_biased_subset(neighbours, neighbours_per_structure);
    else
      neighbours.resize(neighbours_per_structure);
  }

  if (tabular_output || tabular_output_nearest_neighbour_only)
    ;
  else if (! distances_present_in_input && tabular_leader_output)
    id << space_or_tab << clusters_found << space_or_tab << 'P';
  else if (tabular_leader_output)
    id << space_or_tab << clusters_found << space_or_tab << 'P' << space_or_tab << '0';
  else if (' ' == space_or_tab)
    id << cluster << clusters_found << " (" << (neighbours_initially_found + 1) << " members)";
  else
    id << cluster << clusters_found << "\t(" << (neighbours_initially_found + 1) << "\tmembers)";

  int n = neighbours.number_elements();

  for (int i = 0; i < n; i++)
  {
    Smiles_ID_Dist_UID * ni = neighbours[i];

    IWString tmp(ni->id());
    if (tabular_leader_output)
      tmp << space_or_tab << clusters_found << space_or_tab << 'C';
    else
      tmp << cluster << clusters_found << '.' << (i + 1);
    ni->set_id(tmp);
  }

  clusters_found++;

  return process_molecule(smiles, id, neighbours, output);
}

static int
plotnn_gfp_leader (iwstring_data_source & input,
                   IWString_and_File_Descriptor & output)
{
  int fatal;
  while (process_molecule_gfp_leader(input, fatal, output))
  {
    output.write_if_buffer_holds_more_than(IW_FLUSH_BUFFER);
  }

  if (fatal)
    return 0;

  return 1;
}

static int
plotnn (iwstring_data_source & input,
        IWString_and_File_Descriptor & output)
{
  int fatal;

  while (process_molecule(input, fatal, output))
  {
    output.write_if_buffer_holds_more_than(IW_FLUSH_BUFFER);
  }

  if (fatal)
    return 0;

  return 1;
}

static int
plotnn (const char * fname,
        IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  if (from_gfp_leader)
    return plotnn_gfp_leader(input, output);

  return plotnn(input, output);
}

static int
read_possibly_missing_smiles_record (const const_IWSubstring & buffer,
                                     IW_STL_Hash_Map_String & missing_smiles)
{
  const_IWSubstring smiles, id;

  if (! buffer.split(smiles, ' ', id))
  {
    cerr << "Smiles record must have at least two tokens\n";
    return 0;
  }

  id.truncate_at_first(' ');

  missing_smiles[id] = smiles;

//cerr << "Smiles for '" << id << "' is '" << smiles << "'\n";

  return 1;
}

static int
read_possibly_missing_smiles (iwstring_data_source & input,
                              IW_STL_Hash_Map_String & missing_smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_possibly_missing_smiles_record(buffer, missing_smiles))
    {
      cerr << "Erroneous smiles record '" << buffer << "'\n";
      return 0;
    }
  }

  return missing_smiles.size();
}

static int
read_possibly_missing_smiles (const const_IWSubstring & fname,
                              IW_STL_Hash_Map_String & missing_smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open missing smiles file '" << fname << "'\n";
    return 0;
  }

  return read_possibly_missing_smiles(input, missing_smiles);
}

static int
write_normalised_histogram (const IWHistogram & nearest_neighbour_histogram,
                            ostream & stream_for_nearest_neighbour_histogram)
{
  const int b = nearest_neighbour_histogram.nbuckets();

  const unsigned int * raw_counts = nearest_neighbour_histogram.raw_counts();

  unsigned int max_count = raw_counts[0];

  for (int i = 1; i < b; i++)
  {
    if (raw_counts[i] > max_count)
      max_count = raw_counts[i];
  }

  float float_max_count = static_cast<float>(max_count);

  for (int i = 0; i < b; i++)
  {
    float d = static_cast<float>(i) * 0.01;

    float y = static_cast<float>(raw_counts[i]) / float_max_count;

    stream_for_nearest_neighbour_histogram << d << ' ' << y << endl;
  }

  return 1;
}

// Remove?
static int
do_create_julia_file_for_histogram_plot(const IWHistogram & nearest_neighbour_histogram,
                                  const int normalise_h_file,
                                  const int cumulative,
                                  const IWString & fname,
                                  IWString_and_File_Descriptor & output)
{
  const int b = nearest_neighbour_histogram.nbuckets();

  const unsigned int * raw_counts = nearest_neighbour_histogram.raw_counts();

  unsigned int tot = 0;

  int last_non_zero = 0;

  for (int i = 0; i < b; i++)
  {
    if (0 == raw_counts[i])
      continue;

    tot += raw_counts[i];
    last_non_zero = i;
  }

  if (verbose)
    cerr << "Last non zero distance " << (last_non_zero * 0.01) << endl;

  IWString stem = fname;
  assert (stem.ends_with(".r"));
  stem.chop(2);

  IWString png = stem;
  png << ".png";

  output << "x=[0";
  for (int i = 1; i < b; ++i)
  {
    output << ',' << (static_cast<float>(i) * 0.01);
  }
  output << "]\n";

  output << "y=[" << raw_counts[0];
  for (int i = 1; i < b; i++)
  {
    if (normalise_h_file)
      output << ',' << (static_cast<float>(raw_counts[i]) / static_cast<float>(tot));
    else
      output << ',' << raw_counts[i];
  }
  output << "]\n";

  int sum_for_cumulative = raw_counts[0];

  output << "ycum=[" << raw_counts[0];

  for (int i = 1; i < b; i++)
  {
    sum_for_cumulative += raw_counts[i];

    if (normalise_h_file)
      output << ',' << (static_cast<float>(sum_for_cumulative) / static_cast<float>(tot));
    else
      output << ',' << sum_for_cumulative;
  }
  output << "]\n";

  output << "dmax=" << (last_non_zero * 0.01) << '\n';

  output << "using Plots\n";
  output << "p = plot(x,";
  if (cumulative)
    output << "ycum,";
  else
    output << "y,";

  output << "size=(800,600),\n";
  output << "background_color = :ivory,\n";
  output << "fillrange = 0,\n";
  output << "fillalpha = 0.25,\n";
  output << "fillcolor = :lightgoldenrod,\n";
  output << "label = \"\",\n";
  output << "xlims=[0,dmax+0.02],\n";
  output << "xticks=0:0.05:dmax,\n";


  output << "lw=2,color=:red,xlabel=\"Distance\",title=\"" << stem << "\",";
  if (normalise_h_file)
    output << "ylabel=\"Fraction\")\n";
  else
    output << "ylabel=\"Number\")\n";

  output << "png(p,\"" << png << "\")\n";

  return 1;
}

static int
do_create_rfle_for_histogram_plot(const IWHistogram & nearest_neighbour_histogram,
                                  const int normalise_h_file,
                                  const int cumulative,
                                  const IWString & fname,
                                  IWString_and_File_Descriptor & output)
{
  const int b = nearest_neighbour_histogram.nbuckets();

  const unsigned int * raw_counts = nearest_neighbour_histogram.raw_counts();

  unsigned int tot = 0;

  int last_non_zero = 0;

  for (int i = 0; i < b; i++)
  {
    if (0 == raw_counts[i])
      continue;

    tot += raw_counts[i];
    last_non_zero = i;
  }

  if (verbose)
    cerr << "Last non zero distance " << (last_non_zero * 0.01) << endl;

  IWString stem = fname;
  assert (stem.ends_with(".r"));
  stem.chop(2);

  IWString png = stem;
  png << ".png";

  output << "png('" << png << "',width=600,height=500)\n";

  output << "x=c(0";
  for (int i = 1; i < b; ++i)
  {
    output << ',' << (static_cast<float>(i) * 0.01);
  }
  output << ")\n";

  output << "y=c(" << raw_counts[0];
  for (int i = 1; i < b; i++)
  {
    if (normalise_h_file)
      output << ',' << (static_cast<float>(raw_counts[i]) / static_cast<float>(tot));
    else
      output << ',' << raw_counts[i];
  }
  output << ")\n";

  int sum_for_cumulative = raw_counts[0];

  output << "ycum=c(" << raw_counts[0];

  for (int i = 1; i < b; i++)
  {
    sum_for_cumulative += raw_counts[i];

    if (normalise_h_file)
      output << ',' << (static_cast<float>(sum_for_cumulative) / static_cast<float>(tot));
    else
      output << ',' << sum_for_cumulative;
  }
  output << ")\n";

  output << "dmax=" << (last_non_zero * 0.01) << '\n';

  output << "plot(x,";
  if (cumulative)
    output << "ycum";
  else
    output << 'y';

  output << ",lwd=2,type='l',col='red',xlab='Distance',las=1,xlim=c(0,1),main='" << stem << "',";
  if (normalise_h_file)
    output << "ylab='Fraction')\n";
  else
    output << "ylab='Number')\n";

  output << "dev.off()\n";

  return 1;
}

static int
do_create_rfle_for_histogram_plot(const IWHistogram & nearest_neighbour_histogram,
                                  const int normalise_h_file,
                                  const int cumulative,
                                  IWString & fname)
{
  IWString_and_File_Descriptor output;
  
  if (! output.open(fname.null_terminated_chars()))
  {
    cerr << "do_create_rfle_for_histogram_plot:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_create_rfle_for_histogram_plot(nearest_neighbour_histogram, normalise_h_file, cumulative, fname, output);
//return do_create_julia_file_for_histogram_plot(nearest_neighbour_histogram, normalise_h_file, cumulative, fname, output);
}

static int
plotnn (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vn:t:T:xzD:pO:ru:sf:c:L:m:M:wY:j:X:H:hW:3");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', neighbours_per_structure) || neighbours_per_structure < 0)
    {
      cerr << "The neighbours per structure (-n) option must be a whole positive number\n";
      usage(3);
    }

    if (0 == neighbours_per_structure)
    {
      cerr << "Warning, you have chosen to plot NO neighbours of each target molecule\n";
      suppress_neighbours = 1;
    }
    else if (verbose)
      cerr << "Will plot " << neighbours_per_structure << " neighbours of each target molecule\n";
  }

  if (cl.option_present('w'))
  {
    max_possible_distance = numeric_limits<float>::max();

    if (verbose)
      cerr << "Distances can be outside the range [0,1]\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_neighbours_per_structure) || min_neighbours_per_structure < 1)
    {
      cerr << "The minimum number of molecules per target (-m option) must be a whole positive number\n";
      usage(17);
    }

    if (verbose)
      cerr << "At least " << min_neighbours_per_structure << " neighbours per molecule processed\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.option_present('n'))
    {
      cerr << "The pad option (-p) only makes sense with the -n option\n";
      usage(11);
    }

    pad_neighbour_list = 1;

    if (verbose)
      cerr << "All neighbour lists will be padded to " << neighbours_per_structure << " structures\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', lower_distance_threshold) || lower_distance_threshold < 0.0 || lower_distance_threshold >= max_possible_distance)
    {
      cerr << "The lower distance threshold (-t) option must be a valid distance\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will discard neighbours " << lower_distance_threshold << " or closer\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', upper_distance_threshold) || upper_distance_threshold < lower_distance_threshold || upper_distance_threshold >= max_possible_distance)
    {
      cerr << "The upper distance threshold (-T) option must be a valid distance\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will discard neighbours " << upper_distance_threshold << " or further\n";
  }

  if (cl.option_present('W'))
  {
    const_IWSubstring w = cl.string_value('W');

    const_IWSubstring n, d;
    if (! w.split(n, ',', d) || 0 == n.length() || 0 == d.length())
    {
      cerr << "Invalid -W specification '" << w << "'\n";
      usage(4);
    }

    if (! n.numeric_value(number_needed_within_distance) || number_needed_within_distance < 0)
    {
      cerr << "The number specification to the -W option must be a whole non-negativeve number\n";
      usage(5);
    }

    if (! d.numeric_value(distance_for_number_needed_within_distance) || distance_for_number_needed_within_distance < 0.0)
    {
      cerr << "The distance specification to the -W option must be +ve\n";
      usage(6);
    }

    if (verbose)
      cerr << "Must have at least " << number_needed_within_distance << " neighbours within " << distance_for_number_needed_within_distance << " for output\n";
  }

  if (cl.option_present('O'))
  {
    if (! cl.value('O', offset) || offset < 0.0)
    {
      cerr << "Offset values (-O) must be >= 0.0\n";
      usage(7);
    }

    if (verbose)
      cerr << "Distances will be written as " << offset << " - distance\n";
  }

  if (cl.option_present('z'))
  {
    write_molecules_with_no_neighbours = 0;

    if (verbose)
      cerr << "Molecules with no neighbours not written\n";
  }

  if (cl.option_present('x'))     // must be handled after the -n option
  {
    print_target_molecules = 0;

    if (verbose)
      cerr << "Printing of target molecules suppressed\n";

    if (cl.option_present('n') && 0 == neighbours_per_structure)
    {
      cerr << "HUH, you have selected zero neighbours (-n 0) and to NOT output the target molecule, no output is possible\n";
      return 2;
    }
  }

  if (cl.option_present('3'))
  {
    three_column_output = 1;

    if (verbose)
      cerr << "Will produce three column output\n";
  }

  if (cl.option_present('D'))
  {
    const char * fname = cl.option_value('D');

    if (! distance_scaling.build(fname))
    {
      cerr << "Cannot build distance scaling data from '" << fname << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Distance scaling initialised '" << fname << "'\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.option_present('n'))
    {
      cerr << "The choose neighbours at random (-r) option doesn't make sense without the -n option\n";
      usage(3);
    }

    choose_neighbours_at_random = 1;

    if (verbose)
      cerr << "A random sampling of near neighbours will be shown\n";

    iw_random_seed();
  }

  if (cl.option_present('u'))
  {
    int i = 0;
    const_IWSubstring u;
    while (cl.value('u', u, i++))
    {
      if ("id" == u)
      {
        unique_neighbours_only = UNIQUE_BY_ID;
        if (verbose)
          cerr << "Only neighbours with unique ID's will be shown\n";
      }
      else if ("smi" == u)
      {
        unique_neighbours_only = UNIQUE_BY_SMILES;
        if (verbose)
          cerr << "Only neighbours with unique smiles will be shown\n";
      }
      else if ("nochiral" == u)
      {
        set_include_chiral_info_in_smiles(0);
  
        unique_neighbours_only = UNIQUE_BY_SMILES;
        if (verbose)
          cerr << "Only neighbours with unique smiles will be shown\n";
      }
      else if (u.starts_with("APP="))
      {
        u.remove_leading_chars(4);
        append_to_non_unique_neighbours = u;

        if (verbose)
          cerr << "Will append '" << append_to_non_unique_neighbours << "' to non-unique neighbours\n";
      }
      else if ("help" == u)
      {
        display_dash_u_options(cerr);
        return 0;
      }
      else
      {
        cerr << "Unrecognised -u qualifier '" << u << "'\n";
        usage(14);
      }
    }

    if (append_to_non_unique_neighbours.length() && ! unique_neighbours_only)
    {
      cerr << "Append specified, but no means for determining uniqueness\n";
      usage(14);
    }
  }

  if (cl.option_present('s'))
  {
    collect_statistics = 1;

    if (verbose)
      cerr << "Will print neigbhbour statistics with each target\n";
  }

  if (cl.option_present('f'))
  {
    if (! cl.value('f', vf_per_page) || vf_per_page < 1)
    {
      cerr << "The vf plot per page (-f) option must be followed by a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will display the target molecule every " << vf_per_page << " structures\n";
  }

  if (cl.option_present('j'))
  {
    IWString j = cl.string_value('j');
    if (! char_name_to_char(j))
    {
      cerr << "Invalid character specifier (-j) '" << j << "'\n";
      return 1;
    }

     neighbour_separator = j[0];

    if (verbose)
      cerr << "Will separate neighbours with '" << neighbour_separator << "'\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', output_precision) || output_precision < 1)
    {
      cerr << "The output precision must be a whole positive number\n";
      usage(6);
    }

    if (verbose)
      cerr << "Output precision " << output_precision << endl;
  }

  if (! cl.option_present('w'))
  {
    if (three_column_output)
    {
      IWString tmp(three_column_output_separator);
      fraction_as_string.set_leading_string(tmp);
    }

    fraction_as_string.initialise(0.0, 1.0, output_precision);

    if (! three_column_output)
      fraction_as_string.append_to_each_stored_string(neighbour_separator);
  }

  if (cl.option_present('L'))
  {
    from_gfp_leader = 1;

    int i = 0;
    const_IWSubstring l;
    while (cl.value('L', l, i++))
    {
      if ("def" == l)
      {
      }
      else if ("tbl" == l)
      {
        tabular_leader_output = 1;
        append_neighbour_number_to_each_neighbour = 0;

        if (verbose)
          cerr << "Will produce tabular leader output\n";
      }
      else
      {
        cerr << "Unrecognised -L qualifier '" << l << "'\n";
        usage(11);
      }
    }
  }

  if (cl.option_present('Y'))
  {
    append_to_target_record = cl.string_value('Y');

    if (verbose)
      cerr << "Will append '" << append_to_target_record << "' to each target molecule\n";
  }

  discard_self_neighbours = cl.option_count('h');

  if (0 == discard_self_neighbours)
    ;
  else if (1 == discard_self_neighbours)
  {
    if (verbose)
      cerr << "Will discard neighbours with zero distance and the same id as the target\n";
  }
  else
  {
    if (verbose)
      cerr << "Will discard neighbours with same id, regardless of distance\n";
  }

  int normalise_h_file = 0;
  int cumulative_rfile = 0;

  IWString rfile_for_histogram_plot;

  if (cl.option_present('X'))
  {
    int i = 0;
    const_IWSubstring x;
    while (cl.value('X', x, i++))
    {
      if ("nonnum" == x)
      {
        append_neighbour_number_to_each_neighbour = 0;

        if (verbose)
          cerr << "No neighbour numbers with each neighbour\n";
      }
      else if (x.starts_with("tcol="))
      {
        x.remove_leading_chars(5);
        if (! x.numeric_value(target_column_to_write) || target_column_to_write < 1)
        {
          cerr << "Invalid target column to write '" << x << "'\n";
          usage(5);
        }
        if (verbose)
          cerr << "Target identifier in column " << target_column_to_write << endl;

        target_column_to_write--;
      }
      else if (x.starts_with("ncol="))
      {
        x.remove_leading_chars(5);
        const_IWSubstring token;
        for (int j = 0; x.nextword(token, j, ','); )
        {
          int c;
          if (! token.numeric_value(c) || c < 1)
          {
            cerr << "Invalid neighbour column to write '" << token << "'\n";
            usage(1);
          }

          neighbour_columns_to_write.add_if_not_already_present(c - 1);
          if (verbose)
            cerr << "Will write neighbour name token " << c << endl;
        }
      }
      else if (x.starts_with("minextra="))
      {
        x.remove_leading_chars(9);
        if (! x.numeric_value(min_extra_atoms) || min_extra_atoms < 0)
        {
          cerr << "The mininum atom count difference must be a whole +ve number\n";
          usage(5);
        }

        if (verbose)
          cerr << "Neighbours must have at least " << min_extra_atoms << " extra atoms over the target\n";
      }
      else if (x.starts_with("maxextra="))
      {
        x.remove_leading_chars(9);
        if (! x.numeric_value(max_extra_atoms) || max_extra_atoms < min_extra_atoms)
        {
          cerr << "The mininum atom count difference must be a whole +ve number > " << min_extra_atoms << "\n";
          usage(5);
        }

        if (verbose)
          cerr << "Will discard neighbours that have " << max_extra_atoms << " more atoms than the target\n";
      }
      else if ("TABS" == x)
      {
        space_or_tab = '\t';
        cluster = "\tCLUSTER\t";
      }
      else if ("nosmi" == x)
      {
        smiles_tag.resize(0);
         
        if (verbose)
          cerr << "No smiles in input\n";
      }
      else if ("brs" == x)
      {
        if (! cl.option_present('n'))
        {
          cerr << "The biased random neighbours option doesn't make sense without the -n option\n";
          usage(3);
        }

        biased_subset = 1;

        if (verbose)
          cerr << "Will choose a random set of neighbours, biased toward close neighbours\n";

        iw_random_seed();
      }
      else if ("random" == x)
      {
        choose_neighbours_at_random = 1;

        if (verbose)
          cerr << "A random sampling of near neighbours will be shown\n";

        iw_random_seed();
      }
      else if (x.starts_with("SMITAG="))
      {
        smiles_tag = x;
        smiles_tag.remove_leading_chars(7);

        if (verbose)
          cerr << "Smiles tag set to '" << smiles_tag << "'\n";

        if (! smiles_tag.ends_with('<'))
          smiles_tag.add('<');
      }
      else if (x.starts_with("IDTAG="))
      {
        identifier_tag = x;
        identifier_tag.remove_leading_chars(6);

        if (verbose)
          cerr << "Identifier tag set to '" << identifier_tag << "'\n";

        if (! identifier_tag.ends_with('<'))
          identifier_tag.add('<');
      }
      else if (x.starts_with("DISTAG="))
      {
        distance_tag = x;
        distance_tag.remove_leading_chars(7);

        if (verbose)
          cerr << "Distance tag set to '" << distance_tag << "'\n";

        if (! distance_tag.ends_with('<'))
          distance_tag.add('<');
      }
      else if (x.starts_with("SCALE="))
      {
        scale_tag = x;
        scale_tag.remove_leading_chars(6);

        if (verbose)
          cerr << "Scale tag set to '" << scale_tag << "'\n";

        if (! scale_tag.ends_with('<'))
          scale_tag << '<';
      }
      else if (x.starts_with("NOTW="))
      {
        IWString tmp(x);
        tmp.remove_leading_chars(5);
        if (! tmp.ends_with(".smi"))
          tmp << ".smi";
        if (! stream_for_targets_not_written.open(tmp.null_terminated_chars()))
        {
          cerr << "Cannot open stream for targets not written '" << tmp << "'\n";
          return 3;
        }

        if (verbose)
          cerr << "Targets not otherwise output, written to '" << tmp << "'\n";
      }
      else if (x.starts_with("INLS="))
      {
        x.remove_leading_chars(5);
        if (! x.numeric_value(initial_nbr_list_size) || initial_nbr_list_size < 1)
        {
          cerr << "The initial neighbour list size (INLS) must be a whole +ve number\n";
          return 3;
        }

        if (verbose)
          cerr << "Neighbour list initially sized to " << initial_nbr_list_size << endl;
      }
      else if ("table" == x)
      {
        tabular_output = 1;

        if (verbose)
          cerr << "Will produce tabular output of near neighbour distances\n";

        neighbours_per_structure = 0;
        suppress_neighbours = 1;
        collect_statistics = 1;

        if (offset > 0.0)
        {
          cerr << "Offset does not work with tabular output\n";
          return 3;
        }
      }
      else if ("table1" == x)
      {
        tabular_output_nearest_neighbour_only = 1;

        if (verbose)
          cerr << "Will produce tabular output of near neighbour distances\n";

        neighbours_per_structure = 0;
        suppress_neighbours = 1;
        collect_statistics = 1;

        if (offset > 0.0)
        {
          cerr << "Offset does not work with tabular output\n";
          return 3;
        }
      }
      else if ("normh" == x)
      {
        normalise_h_file = 1;

        if (verbose)
          cerr << "The -H file will be normalised\n";
      }
      else if ("cumh" == x)
      {
        cumulative_rfile = 1;

        if (verbose)
          cerr << "The -HR= file will be cumulative\n";
      }
      else if ("allh" == x)
      {
        place_all_distances_in_nn_histogram = 1;

        if (verbose)
          cerr << "All distances placed into histogram\n";
      }
      else if (x.starts_with("HR="))
      {
        x.remove_leading_chars(3);
        rfile_for_histogram_plot = x;
        if (! rfile_for_histogram_plot.ends_with(".r"))
          rfile_for_histogram_plot << ".r";
      }
      else if ("ckds" == x)
      {
        check_distance_ordering = 1;
        if (verbose)
          cerr << "WIll check distance ordering\n";
      }
      else if ('0' == x)
      {
        ignore_leading_zeros_in_identifiers = 1;

        if (verbose)
          cerr << "Will remove leading zeros from identifiers when doing id comparisons\n";

        if (0 == discard_self_neighbours)
          discard_self_neighbours = 1;
      }
      else if (x.starts_with("sdcensor="))
      {
        x.remove_leading_chars(9);
        if (! x.numeric_value(censor_distances_shorter_than) || censor_distances_shorter_than <= 0.0f || censor_distances_shorter_than >= 1.0f)
        {
          cerr << "The short distance censor value (sdcensor=) must be a valid distance\n";
          return 3;
        }

        if (verbose)
          cerr << "Will censor distances shorter than " << censor_distances_shorter_than << endl;
      }
      else if (x.starts_with("mnc="))
      {
        x.remove_leading_chars(4);
        if (! x.numeric_value(mandatory_neighbour_count) || mandatory_neighbour_count < 0)
        {
          cerr << "The mandatory neighbour count flag must be a whole non negative number\n";
          usage(1);
        }

        if (verbose)
          cerr << "Will only write a target if it has precisely " << mandatory_neighbour_count << " neighbours\n";
      }
      else if ("nodist" == x)
      {
        distances_present_in_input = 0;

        if (verbose)
          cerr << "No distance values present in input\n";
      }
      else if ("ftn" == x)
      {
        take_first_token_of_name_field = 1;
        if (verbose)
          cerr << "Will take only the first token of name fields\n";
      }
      else if ("sdo" == x)
      {
        write_shortest_distance_with_stats = 1;

        if (verbose)
          cerr << "Will write the shortest distance only with the -s output\n";
      }
      else if ("help" == x)
      {
        display_dash_x_options(cerr);
        return 0;
      }
      else
      {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        display_dash_x_options(cerr);
        return 4;
      }
    }
  }

  if (append_neighbour_number_to_each_neighbour)
  {
    iwdigits.initialise(200);
    iwdigits.append_to_each_stored_string(space_or_tab);
  }

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (! read_possibly_missing_smiles(m, missing_smiles))
      {
        cerr << "Cannot process smiles file '" << m << "'\n";
        return 4;
      }
    }

    if (verbose)
      cerr << "Read " << missing_smiles.size() << " smiles\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "INsufficient arguments\n";
    usage(2);
  }

  ofstream stream_for_nearest_neighbour_histogram;

  if (cl.option_present('H'))
  {
    const char * h = cl.option_value('H');

    stream_for_nearest_neighbour_histogram.open(h, ios::out);
    if (! stream_for_nearest_neighbour_histogram.good())
    {
      cerr << "Cannot open histogram stream file '" << h << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Histogram written to '" << h << "'\n";

    nearest_neighbour_histogram.initialise(0.0, 1.0, 0.01);
  }

  IWString_and_File_Descriptor output(1);

  if (tabular_output)
  {
    output << "ID ";
    if (from_gfp_leader)
    {
      output << "CSIZE";
    }
    else
    {
      output << "NBRS";
    }

    output << " Min Max\n";
  }
  else if (tabular_output_nearest_neighbour_only)
  {
    output << "ID Min\n";
  }

  int rc = 0;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! plotnn(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (nearest_neighbour_histogram.active())
  {
    if (normalise_h_file)
      write_normalised_histogram(nearest_neighbour_histogram, stream_for_nearest_neighbour_histogram);
    else
      nearest_neighbour_histogram.write_terse(stream_for_nearest_neighbour_histogram, 0);

    if (rfile_for_histogram_plot.length())
      do_create_rfle_for_histogram_plot(nearest_neighbour_histogram, normalise_h_file, cumulative_rfile, rfile_for_histogram_plot);
  }

  if (verbose)
  {
    cerr << "Read " << molecules_processed << " molecules\n";
    cerr << "Molecules had between " << neighbour_statistics.minval() << " and " << neighbour_statistics.maxval() << " neighbours";
    if (neighbour_statistics.n() > 1 && neighbour_statistics.minval() < neighbour_statistics.maxval())
      cerr << " ave " << neighbour_statistics.average();
    cerr << endl;

    cerr << "distances between " << distance_stats.minval() << " and " << distance_stats.maxval();
    if (distance_stats.n() > 1)
      cerr << " ave " << distance_stats.average();
    cerr << endl;

    cerr << "Nearest distances between " << nearest_neighbour_stats.minval() << " and " << nearest_neighbour_stats.maxval();
    if (nearest_neighbour_stats.n() > 1)
      cerr << " ave " << nearest_neighbour_stats.average();
    cerr << endl;

    cerr << zero_distance_neighbours << " exact matches\n";

    cerr << "Furthest distances between " << furthest_neighbour_stats.minval() << " and " << furthest_neighbour_stats.maxval();
    if (furthest_neighbour_stats.n() > 1)
      cerr << " ave " << furthest_neighbour_stats.average();
    cerr << endl;

    for (int i = 0; i < neighbour_count.number_elements(); i++)
    {
      if (neighbour_count[i])
        cerr << neighbour_count[i] << " molecules had " << i << " neighbours\n";
    }

    if (duplicate_neighbours_suppressed)
      cerr << duplicate_neighbours_suppressed << " duplicate neighbours suppressed by -u option\n";

    if (three_column_output)
      cerr << "Wrote " << molecules_written << " distance values\n";
    else
      cerr << "Wrote " << molecules_written << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = 0; 

  try
  {
    rc = plotnn(argc, argv);
  }
  catch ( const Fatal_Error & f)
  {
    cerr << "Caught " << f << endl;
  }

  return rc;
}
