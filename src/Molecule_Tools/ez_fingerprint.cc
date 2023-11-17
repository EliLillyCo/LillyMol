/*
  Produces sparse fingerprints based on discovered cis-trans bonds
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwbits/iwbits.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/output.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_fingerprinted = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static angle_t tolerance = static_cast<angle_t> (10.0 * DEG2RAD);

static int use_coordinates = 0;

static int work_as_filter = 0;

static int reduce_to_largest_fragment = 0;

static int discover_new_queries = 1;

/*
  We can control just how different different configurations across a
  bond are perceived to be.
*/

static int positive_direction_value = 2;
static int negative_direction_value = 1;

/*
  We typically generate just a couple of fingerprints. These can get overwhelmed
  by larger numbers of bits from other fingerprints. We can create multiple bits
*/

static int replicate_bits = 1;

int nq = 0;

static Molecule_Output_Object stream_for_matches[2];

static IWString fingerprint_tag("NCEZF<");

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

static resizable_array_p<EZ_Matches> ez_matches;

template class resizable_array_p<EZ_Matches>;
template class resizable_array_base<EZ_Matches *>;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces descriptor(s) based on E/Z bonding configuration\n";
  cerr << "  -q <qry>      query file specification\n";
  cerr << "  -s <smarts>   give query as smarts\n";
  cerr << "  -J <tag>      write results as a fingerprint with tag <tag>\n";
  cerr << "  -f            work as a TDT filter\n";
  cerr << "  -S <fname>    write final queries to <fname>\n";
  cerr << "  -c            use the molecule's coordinates\n";
  cerr << "  -d <number>   difference between configurations (between 1 and 254)\n";
  cerr << "  -l            strip to largest fragment\n";
  cerr << "  -n            do NOT discover new queries\n";
  cerr << "  -i ...        input specification(s)\n";
  cerr << "  -o ...        output specifications\n";
  cerr << "  -A ...        aromaticity specifications, enter '-A help' for details\n";
  cerr << "  -E ...        element specifications, enter '-E help' for details\n";
  cerr << "  -v            verbose output\n";

  exit (rc);
}

static int
write_queries (const resizable_array_p<Substructure_Hit_Statistics> & queries,
               std::ostream & output)
{
  for (int i = 0; i  < queries.number_elements(); i++)
  {
    queries[i]->write_msi(output);
  }

  return 1;
}

static int
write_queries (const resizable_array_p<Substructure_Hit_Statistics> & queries,
               const char * fname)
{
  std::ofstream output(fname, std::ios::out);

  if (! output.good())
  {
    cerr << "write_queries:cannot open '" << fname << "'\n";
    return 0;
  }

  return write_queries (queries, output);
}

static int
compute_invariant (Molecule & m,
                   atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 100000 * acon;

  if (m.is_aromatic(zatom))
    rc += 20000;
  else if (m.is_ring_atom(zatom))
    rc += 10000;

  int ahc = 0;
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (6 != m.atomic_number(j))
      ahc++;
  }

  rc += 500 * ahc;

  rc += a->atomic_number();

  return rc;
}

static int
copy_bond_attributes (const Bond & b,
                      Substructure_Bond & s)
{
  if (b.is_aromatic())
    s.set_bond_type(AROMATIC_BOND);
  else if (b.is_single_bond())
    s.set_bond_type(SINGLE_BOND);
  else if (b.is_double_bond())
    s.set_bond_type(DOUBLE_BOND);
  else
  {
    cerr << "Huh, strange bond between " << b.a1() << " and " << b.a2() << endl;
    s.set_bond_type(SINGLE_BOND|DOUBLE_BOND);   // makes no sense, should abort...
    return 0;
  }

  return 1;
}

static int
fill_query_atom_properties(Molecule & m,
                           atom_number_t zatom,
                           Substructure_Atom & s)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  s.set_element(a->element());

  s.set_ncon(acon);

  s.set_attached_heteroatom_count(m.attached_heteroatom_count(zatom));

  if (m.is_aromatic(zatom))
    s.set_aromaticity(AROMATIC);
  else if (m.is_ring_atom(zatom))
    s.set_min_nrings(1);
  else
    s.set_nrings(0);

  return 1;
}

/*
  a1        a4
    \      /
     a2==a3
*/

static int
setup_query_for_matched_atoms (Molecule & m,
                               atom_number_t a1,
                               atom_number_t a2,
                               atom_number_t a3,
                               atom_number_t a4)
{
  Substructure_Hit_Statistics * c = new Substructure_Hit_Statistics;
  c->set_comment(m.name());

  Single_Substructure_Query * q = new Single_Substructure_Query;
  q->set_comment(m.name());

  Substructure_Atom * q1 = new Substructure_Atom;
  Substructure_Atom * q2 = new Substructure_Atom;
  Substructure_Atom * q3 = new Substructure_Atom;
  Substructure_Atom * q4 = new Substructure_Atom;

  q1->set_initial_atom_number(0);
  q2->set_initial_atom_number(1);
  q3->set_initial_atom_number(2);
  q4->set_initial_atom_number(3);

  fill_query_atom_properties (m, a1, *q1);
  fill_query_atom_properties (m, a2, *q2);
  fill_query_atom_properties (m, a3, *q3);
  fill_query_atom_properties (m, a4, *q4);

  Substructure_Bond * b12 = new Substructure_Bond;
  Substructure_Bond * b23 = new Substructure_Bond;
  Substructure_Bond * b34 = new Substructure_Bond;

  const Bond * b = m.bond_between_atoms(a1, a2);
  copy_bond_attributes(*b, *b12);

  b = m.bond_between_atoms(a2, a3);
  copy_bond_attributes(*b, *b23);

  b = m.bond_between_atoms(a3, a4);
  copy_bond_attributes(*b, *b34);

  q2->set_parent(q1, b12);
  q1->add_child(q2);
  q3->set_parent(q2, b23);
  q2->add_child(q3);
  q4->set_parent(q3, b34);
  q3->add_child(q4);

  q->add_root_atom(q1);
  c->add(q);

  if (verbose > 1)
    cerr << "Found new cis-trans bond in " << m.name() << "'\n";

//std::ofstream newly_found("newly_found.qry", std::ios::out);
//c->write_msi(newly_found);

  queries.add(c);
  ez_matches.add(new EZ_Matches);
  nq++;

  assert (nq == queries.number_elements());

  return nq;
}

/*
  a1        a5
    \      /
     a3==a4
    /      \
  a2        a6
*/

static int
setup_query_for_matched_atoms (Molecule & m,
                               atom_number_t a3, 
                               atom_number_t a4)
{
  const Atom * a = m.atomi(a3);

  int acon = a->ncon();

  if (acon > 3)   // how could that happen
    return 0;

  atom_number_t a1 = INVALID_ATOM_NUMBER;
  atom_number_t a2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a3, i);

    if (j == a4)
      continue;

    if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else
      a2 = j;
  }

  atom_number_t a5 = INVALID_ATOM_NUMBER;
  atom_number_t a6 = INVALID_ATOM_NUMBER;

  a = m.atomi(a4);

  acon = a->ncon();

  if (acon > 3)
    return 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a4, i);

    if (j == a3)
      continue;

    if (INVALID_ATOM_NUMBER == a5)
      a5 = j;
    else
      a6 = j;
  }

// We need to identify 4 atoms to form a query. That's hard...

  if (INVALID_ATOM_NUMBER == a2 && INVALID_ATOM_NUMBER == a6)
    return setup_query_for_matched_atoms (m, a1, a3, a4, a5);

// We can make things easier for ourselves if we ensure that the LHS is the 3 connected form

  if (INVALID_ATOM_NUMBER == a2)
  {
    std::swap(a3, a4);
    std::swap(a1, a5);
    std::swap(a2, a6);
  }

  int inv1 = compute_invariant(m, a1);
  int inv2 = compute_invariant(m, a2);

  if (inv1 == inv2)    // cannot be resolved
    return 0;

  return setup_query_for_matched_atoms (m, a1, a3, a4, a5);
}

static int
do_write_fingerprints (Molecule & m,
                       const Sparse_Fingerprint_Creator & sfc,
                       IWString_and_File_Descriptor & output)
{
  if (! work_as_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);

  output << tmp << '\n';

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
               Sparse_Fingerprint_Creator & sfp,
               int bit_number,
               const Set_of_Atoms & s)
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

    return 0;
  }

  int zresult;
  if (use_coordinates)
    zresult = m.ez_by_geometry (s[0], s1, s2, s[3], tolerance);
  else
    zresult = ez_descriptor (m, s[0], s1, s2, s[3]);

  if (0 == zresult)
    return 0;

  int c;
  if (zresult > 0)
  {
    c = positive_direction_value;
    ez_matches[bit_number]->extra(positive_direction_value);
  }
  else
  {
    c = negative_direction_value;
    ez_matches[bit_number]->extra(negative_direction_value);
  }

  for (int i = 0; i < replicate_bits; i++)
  {
    sfp.hit_bit(bit_number + i * 100000, c);
  }

  return 1;
}

static int
ez_descriptor (Molecule & m,
               Sparse_Fingerprint_Creator & sfp)
{
  if (! m.cis_trans_bonds_present())
    return 0;

  int queries_matching = 0;

  Molecule_to_Match target(&m);

  int matoms = m.natoms();

  resizable_array<int> hit_by_queries;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    if (0 == queries[i]->substructure_search(target, sresults))
      continue;

    const Set_of_Atoms * s = sresults.embedding(0);

    atom_number_t a2 = s->item(1);
    atom_number_t a3 = s->item(2);

    if (a2 > a3)
      std::swap(a2, a3);

    if (! hit_by_queries.add_if_not_already_present(a2 * matoms + a3))
      continue;

//  cerr << m.name() << " matches atoms " << (*s) << endl;
    ez_descriptor (m, sfp, i, *s);

    queries_matching++;
  }

  if (! discover_new_queries)
    return queries_matching;

// We keep track of whether or not any new queries are defined

  int initial_number_queries = nq;

  int nbonds = m.nedges();

  for (int i = 0; i < nbonds; i++)
  {
    const Bond * b = m.bondi(i);

    if (! b->is_double_bond())
      continue;

    if (! b->part_of_cis_trans_grouping())
      continue;

    if (b->is_cis_trans_either_double_bond())
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (a1 > a2)
      std::swap(a1, a2);

    if (hit_by_queries.contains(a1 * matoms + a2))
      continue;

    setup_query_for_matched_atoms(m, a1, a2);
  }

  if (initial_number_queries == queries.number_elements())
    return queries_matching;

  for (int i = initial_number_queries; i < nq; i++)
  {
    Substructure_Results sresults;

    if (0 == queries[i]->substructure_search (target, sresults))
      continue;

    const Set_of_Atoms * s = sresults.embedding(0);

    atom_number_t a1 = s->item(1);
    atom_number_t a2 = s->item(2);

    if (a1 > a2)
      std::swap(a1, a2);

    if (! hit_by_queries.add_if_not_already_present(a1 * matoms + a2))
      continue;

    queries_matching++;

//  cerr << m.name() << " matches atoms " << (*s) << endl;
    ez_descriptor (m, sfp, i, *s);
  }

  return queries_matching;
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

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  return;
}

static int
ez_descriptor (Molecule & m,
               IWString_and_File_Descriptor & output)
{
  preprocess(m);

  const int has_coordinates = determine_has_coordinates(m);

  if (! has_coordinates && use_coordinates)
  {
    cerr << "Molecule '" << m.name() << "' has no coordinates, but request to use coordinates\n";
    return 0;
  }

  Sparse_Fingerprint_Creator sfp;

  if (m.natoms() >= 4)    // must have at least 4 atoms for a cis-trans bond
  {
    if (ez_descriptor (m, sfp))
      molecules_fingerprinted++;
  }

  return do_write_fingerprints (m, sfp, output);
}

static int
ez_descriptor (data_source_and_type<Molecule> & input,
               IWString_and_File_Descriptor & output)
{
  Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m (m);

    molecules_read++;

    if (! ez_descriptor (*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static int
ez_descriptor (const const_IWSubstring & buffer,
               IWString_and_File_Descriptor & output)
{
  Molecule m;
  if (! m.build_from_smiles (buffer))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return ez_descriptor (m, output);
}

static int
ez_descriptor (iwstring_data_source & input,
               IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with (smiles_tag))
      continue;

    molecules_read++;

    buffer.remove_leading_chars (smiles_tag.length());
    assert (buffer.ends_with ('>'));
    buffer.chop();

    if (! ez_descriptor (buffer, output))
    {
      cerr << "Fatal error processing line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static int
ez_descriptor (const char * fname,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ez_descriptor (input, output);
}

static int
ez_descriptor (const char * fname,
               FileType input_type,
               IWString_and_File_Descriptor & output)
{
  if (work_as_filter)
    return ez_descriptor (fname, output);

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

  return ez_descriptor (input, output);
}

static int
ez_descriptor (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vi:E:A:q:s:J:cflS:d:nr:");

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

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present ('f'))
  {
    work_as_filter = 1;

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

  if (cl.option_present('d'))
  {
    int d;
    if (! cl.value('d', d) || d <= 0)
    {
      cerr << "The difference between configurations option (-d) must be a whole +ve number\n";
      usage(4);
    }

    if (d > 254)
    {
      cerr << "Sorry, the difference between configurations option (-d) must be less than 255\n";
      return 4;
    }

    negative_direction_value = 1;
    positive_direction_value = 1 + d;

    if (verbose)
      cerr << "Negative direction " << negative_direction_value << " positive " << positive_direction_value << '\n';
  }

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
    cerr << "Defined " << nq << " starting queries\n";

  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_unique_embeddings_only (1);
    ez_matches.add(new EZ_Matches);
  }

  if (cl.option_present('n'))
  {
    discover_new_queries = 0;

    if (0 == nq)
    {
      cerr << "NO queries defined, but -n specified, cannot continue\n";
      usage(3);
    }

    if (verbose)
      cerr << "Identification of new motifs suppressed\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', replicate_bits) || replicate_bits < 1)
    {
      cerr << "The replicate bits option (-r) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will create " << replicate_bits << " replicates of each bit\n";
  }

  if (cl.option_present ('J'))
  {
    cl.value ('J', fingerprint_tag);
    
    if (verbose)
      cerr << "Output written as fingerprints, tag '" << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag.add ('<');

    if (verbose)
      cerr << "Fingerprints created with tag '" << fingerprint_tag << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ez_descriptor (cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (cl.option_present('S'))
  {
    IWString s = cl.string_value('S');

    if (! s.ends_with(".qry"))
      s << ".qry";

    if (! write_queries(queries, s.null_terminated_chars()))
    {
      cerr << "Cannot write queries to '" << s << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Wrote queries to '" << s << "'\n";
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, finished with " << queries.number_elements() << " queries\n";
    cerr << molecules_fingerprinted << " molecules fingerprinted\n";
    for (int i = 0; i < nq; i++)
    {
      queries[i]->report(cerr, 0);
    }

    assert (nq == queries.number_elements());
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
