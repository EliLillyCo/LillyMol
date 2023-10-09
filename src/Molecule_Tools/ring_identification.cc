/*
  Identify rings in molecules
*/

#include <stdlib.h>
#include <memory>

#include "cmdline.h"
#include "iw_auto_array.h"
#include "iw_stl_hash_set.h"
#include "misc.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "path.h"
#include "target.h"
#include "substructure.h"
#include "output.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;
static int molecules_matching = 0;
static int molecules_written = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int remove_chirality = 0;

static int remove_cis_trans_bonds = 0;

static resizable_array_p<Substructure_Query> must_be_in_ring;

static int number_queries_that_must_match = 0;

static resizable_array_p<Substructure_Query> must_not_be_in_ring;

static Min_Max_Specifier<int> fused_system_size;

static Min_Max_Specifier<int> number_connections;

static Min_Max_Specifier<int> number_spinach_attachments;
static Min_Max_Specifier<int> number_non_spinach_attachments;

static Min_Max_Specifier<int> spinach_size_specification;

static extending_resizable_array<int> nrings;
static int nq = 0;

static int cut_out_matching_rings = 0;

static int cut_out_ring_and_spinach = 0;

static int isotope_label = 0;

static int isotope_for_spinach_connections = 0;

static int isotope_for_embedding_atoms = 0;

static Molecule_Output_Object stream_for_molecules_not_matching;

static int unique_fragments_only = 0;

static IW_STL_Hash_Set already_seen;

static int duplicate_rings_discarded = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Identifies rings in molecules\n";
  cerr << "  -q <query>    query specification for atom(s) that must be in the ring\n";
  cerr << "  -s <smarts>   smarts for atom(s) that must be in the ring\n";
  cerr << "  -n <smarts>   smarts for atom(s) that must NOT be in the ring\n";
//cerr << "  -a <n>        how many of the queries must match\n";
  cerr << "  -f <size>     fused system size for the ring\n";
  cerr << "  -c <n>        total number of connections from ring\n";
  cerr << "  -h <n>        total number of     spinach connections from ring. e.g. -h 0,1,2\n";
  cerr << "  -j <n>        total number of non spinach connections from ring. e.g. -j 1,2\n";
  cerr << "  -H ...        specification of total atoms in spinach\n";
  cerr << "  -M ...        what to do with rings that match - default is write whole molecule\n";
  cerr << "  -M <n>        apply isotopic label <n> to all ring system atoms\n";
  cerr << "  -M spch=<n>   apply isotopic label <n> to points where spinach grows\n";
  cerr << "  -M embd=<n>   apply isotopic label <n> to points where the ring is embedded\n";
  cerr << "  -M cut        cut the ring (system) atoms out of the molecule\n";
  cerr << "  -M cutrs      cut ring (system) and associated spinach\n";
  cerr << "  -u            when excising pieces, only write unique fragments\n";
  cerr << "  -B <fname>    write molecules not matching criteria to <fname>\n";
  cerr << "  -z            remove all chirality and cis/trans bonding information\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

class Ring_Identification
{
  private:
    int * _rs;
    int * _spinach;
    Substructure_Results * _results_for_must_match_queries;
    Substructure_Results * _results_for_must_not_match_queries;

//  private functions

    int  _do_cut_ring_and_spinach (const Molecule & m, Molecule & r) const;
    void _include_all_spinach (const Molecule & m, atom_number_t zatom, int * to_transfer) const;

    int _process (Molecule & m, int system_size, IWString_and_File_Descriptor & output);
    int _handle_match (Molecule & m, int system_size, IWString_and_File_Descriptor & output);
    int _apply_isotopic_labels (Molecule & m) const;
    int _matches_spinach_specification (Molecule & m);

  public:
    Ring_Identification();
    ~Ring_Identification();

    int process (Molecule &, IWString_and_File_Descriptor & output);
};

Ring_Identification::Ring_Identification()
{
  _rs = nullptr;
  _spinach = nullptr;
  _results_for_must_match_queries = nullptr;
  _results_for_must_not_match_queries = nullptr;

  return;
}

Ring_Identification::~Ring_Identification()
{
  if (nullptr != _rs)
    delete [] _rs;

  if (nullptr != _spinach)
    delete [] _spinach;

  if (nullptr != _results_for_must_match_queries)
    delete [] _results_for_must_match_queries;

  if (nullptr != _results_for_must_not_match_queries)
    delete [] _results_for_must_not_match_queries;

  return;
}

/*
  We always return 1
*/

int
Ring_Identification::process (Molecule & m,
                              IWString_and_File_Descriptor & output)
{
  const int nr = m.nrings();

  nrings[nr]++;

  if (0 == nr)
    return 1;

  Molecule_to_Match target(&m);

  _results_for_must_match_queries = new Substructure_Results[nq];

  int hits_this_molecule = 0;

  for (int i = 0; i < nq; i++)
  {
    int nhits = must_be_in_ring[i]->substructure_search(target, _results_for_must_match_queries[i]);

    if (nhits)
      hits_this_molecule++;
    else if (number_queries_that_must_match == nq)
    {
      if (verbose)
        cerr << "In '" << m.name() << "' no hits to query " << i << " '" << must_be_in_ring[i]->comment() << "'\n";
      return 1;
    }
  }

  if (hits_this_molecule < number_queries_that_must_match)
  {
    if (verbose)
      cerr << "In '" << m.name() << "' not enough hits " << hits_this_molecule << " to " << nq << " queries\n";
    return 1;
  }

  int mnb = must_not_be_in_ring.number_elements();
  if (mnb)
  {
    _results_for_must_not_match_queries = new Substructure_Results[mnb];

    for (int i = 0; i < mnb; i++)
    {
      must_not_be_in_ring[i]->substructure_search(target, _results_for_must_not_match_queries[i]);
    }
  }

  const int matoms = m.natoms();

  _spinach = new int[matoms];

  m.identify_spinach(_spinach);

// Now allowance for spiro fusions, too hard, change if needed

  int * ring_already_done = new_int(nr); iw_auto_array<int> free_ring_already_done(ring_already_done);

  _rs = new_int(matoms);

  int matches_this_molecule = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    set_vector(_rs, matoms, 0);

    const Ring * ri = m.ringi(i);

    ri->set_vector(_rs, 1);

    if (! ri->is_fused())
    {
      matches_this_molecule += _process(m, 1, output);
      continue;
    }

    int system_size = 1;

    for (int j = i + 1; j < nr; j++)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        continue;

      rj->set_vector(_rs, 1);

      ring_already_done[j] = 1;

      system_size++;
    }

    matches_this_molecule += _process(m, system_size, output);
  }

  if (matches_this_molecule)
    molecules_matching++;
  else if (stream_for_molecules_not_matching.is_open())
    stream_for_molecules_not_matching.write(m);

  return 1;
}

int
Ring_Identification::_process (Molecule & m,
                               int system_size,
                               IWString_and_File_Descriptor & output)
{
  if (! fused_system_size.matches(system_size))
    return 0;

  for (int i = 0; i < nq; i++)
  {
    const Substructure_Results & sri = _results_for_must_match_queries[i];

    int nhits = sri.number_embeddings();

    int got_match_this_system = 0;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sri.embedding(j);

      atom_number_t k = e->item(0);

      if (_rs[k])
      {
        got_match_this_system = 1;
        break;
      }
    }

    if (! got_match_this_system)
      return 0;
  }

  for (int i = 0; i < must_not_be_in_ring.number_elements(); i++)
  {
    const Substructure_Results & sri = _results_for_must_not_match_queries[i];

    int nhits = sri.number_embeddings();

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sri.embedding(j);

      atom_number_t k = e->item(0);

      if (_rs[k])
        return 0;
    }
  }

  if (number_spinach_attachments.is_set() ||
      number_non_spinach_attachments.is_set() ||
      spinach_size_specification.is_set() ||
      number_connections.is_set())
  {
    if (! _matches_spinach_specification(m))
      return 0;
  }

  return _handle_match (m, system_size, output);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (remove_chirality)
    m.remove_all_chiral_centres();

  if (remove_cis_trans_bonds)
    m.revert_all_directional_bonds_to_non_directional();

  return;
}

static int
is_unique (Molecule & m)
{
  const IWString & usmi = m.unique_smiles();

  if (already_seen.contains(usmi))
  {
    duplicate_rings_discarded++;
    return 0;
  }

  already_seen.insert(usmi);

  m.invalidate_smiles();

  return 1;
}

/*
  This is quite simple because we know we are working with atoms in the molecular spinach
*/

static int
atoms_in_spinach (const Molecule & m,
                  atom_number_t avoid,
                  atom_number_t zatom)
{
  int rc = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  if (1 == acon)
    return 1;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == avoid)
      continue;

    rc += atoms_in_spinach(m, zatom, j);
  }

  return rc;
}

void
Ring_Identification::_include_all_spinach (const Molecule & m,
                                           atom_number_t zatom,
                                           int * to_transfer) const
{
  to_transfer[zatom] = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (to_transfer[j])
      continue;

    _include_all_spinach(m, j, to_transfer);
  }

  return;
}

int
Ring_Identification::_do_cut_ring_and_spinach (const Molecule & m,
                                               Molecule & r) const
{
  int matoms = m.natoms();

  int * to_transfer = new_int(matoms); iw_auto_array<int> free_to_transfer(to_transfer);

  for (int i = 0; i < matoms; i++)
  {
    if (0 == _rs[i])
      continue;

    to_transfer[i] = 1;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (_rs[k])
        continue;

      if (! _spinach[k])
        continue;

      _include_all_spinach(m, k, to_transfer);
    }
  }

  m.create_subset(r, to_transfer);

  return 1;
}

int
Ring_Identification::_apply_isotopic_labels (Molecule & m) const
{
  if (0 == isotope_label && 0 == isotope_for_spinach_connections && 0 == isotope_for_embedding_atoms)
    return 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == _rs[i])
      continue;

    if (isotope_label)
      m.set_isotope(i, isotope_label);   // the default action

    if (0 == isotope_for_spinach_connections && 0 == isotope_for_embedding_atoms)
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (_rs[k])
        continue;

      if (isotope_for_spinach_connections && _spinach[k])
        m.set_isotope(i, isotope_for_spinach_connections);
      else if (isotope_for_embedding_atoms && ! _spinach[k])
        m.set_isotope(i, isotope_for_embedding_atoms);

      break;
    }
  }

  return 1;
}

static int
expand_subset_to_include_doubly_bonded_to_ring (const Molecule & m,
                                                int * subset)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == subset[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    if (2 == acon)
      continue;

    if (a->nbonds() == acon)
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(i);

      if (subset[k])
        continue;

      subset[k] = 1;
      rc++;
    }
  }

  return rc;
}

int
Ring_Identification::_handle_match (Molecule & m,
                                    int system_size,
                                    IWString_and_File_Descriptor & output)
{
  molecules_written++;

  int need_to_transform_to_non_isotopic_form = _apply_isotopic_labels(m);

  if (cut_out_matching_rings)
  {
    expand_subset_to_include_doubly_bonded_to_ring(m, _rs);
    Molecule r;
    m.create_subset(r, _rs);

    if (unique_fragments_only && ! is_unique(r))
      return 1;

    output << r.smiles() << ' ' << m.name() << '\n';
  }
  else if (cut_out_ring_and_spinach)
  {
    Molecule r;
    _do_cut_ring_and_spinach(m, r);

    if (unique_fragments_only && ! is_unique(r))
      return 1;

    output << r.smiles() << ' ' << m.name() << '\n';
  }
  else
    output << m.smiles() << ' ' << m.name() << '\n';

  if (need_to_transform_to_non_isotopic_form)
    m.transform_to_non_isotopic_form();

  return 1;
}

int
Ring_Identification::_matches_spinach_specification (Molecule & m)
{
  int spinach_attachments = 0;
  int non_spinach_attachments = 0;
  int total_atoms_in_spinach = 0;
  int ncon = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == _rs[i])
      continue;

    const Atom * ai = m.atomi(i);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (_rs[k])    // in the ring system
        continue;

      ncon++;
      if (_spinach[k])
      {
        spinach_attachments++;
        int spinach_size = atoms_in_spinach(m, i, k);
        total_atoms_in_spinach += spinach_size;
      }
      else
      {
        non_spinach_attachments++;
      }
    }
  }

  if (! number_connections.is_set())
    ;
  else if (! number_connections.matches(ncon))
    return 0;

  if (! number_spinach_attachments.is_set())
    ;
  else if (! number_spinach_attachments.matches(spinach_attachments))
    return 0;

  if (! number_non_spinach_attachments.is_set())
    ;
  else if (! number_non_spinach_attachments.matches(non_spinach_attachments))
    return 0;

  if (! spinach_size_specification.is_set())
    ;
  else if (! spinach_size_specification.matches(total_atoms_in_spinach))
    return 0;

  return 1;
}

static int
ring_identification (Molecule & m,
                     IWString_and_File_Descriptor & output)
{
  Ring_Identification ri;

  return ri.process(m, output);
}

static int
ring_identification (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ring_identification(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
ring_identification (const char * fname, int input_type, 
                IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return ring_identification(input, output);
}

static int
ring_identification (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:ls:q:n:f:h:j:c:H:M:B:uz");

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

  if (cl.option_present('z'))
  {
    remove_chirality = 1;
    remove_cis_trans_bonds = 1;

    if (verbose)
      cerr << "All chirality and cis trans bonding information removed\n";
  }

  if (cl.option_present('u'))
  {
    if (! cl.option_present('M'))
    {
      cerr << "The unique fragments only option only makes sense with the -M option\n";
      usage(2);
    }

    unique_fragments_only = 1;

    if (verbose)
      cerr << "Will only write unique fragments\n";
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, must_be_in_ring, verbose, 'q'))
    {
      cerr << "Cannot process queries for atoms that must be in the ring\n";
      usage(2);
    }
  }

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring s;
    while (cl.value ('s', s, i++))
    {
      Substructure_Query * tmp = new Substructure_Query;

      if (! tmp->create_from_smarts (s))
      {
        cerr << "Invalid smarts '" << s << "'\n";
        return 6;
      }

      must_be_in_ring.add (tmp);
    }
  }

  if (cl.option_present('n'))
  {
    int i = 0;
    const_IWSubstring n;
    while (cl.value ('n', n, i++))
    {
      Substructure_Query * tmp = new Substructure_Query;

      if (! tmp->create_from_smarts (n))
      {
        cerr << "Invalid smarts '" << n << "'\n";
        return 6;
      }

      must_not_be_in_ring.add (tmp);
    }

    for (int i = 0; i < must_not_be_in_ring.number_elements(); i++)
    {
      must_not_be_in_ring[i]->set_find_unique_embeddings_only(1);
    }
  }

  nq = must_be_in_ring.number_elements();

  number_queries_that_must_match = nq;

  for (int i = 0; i < nq; i++)
  {
    must_be_in_ring[i]->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('f'))
  {
    const_IWSubstring f = cl.string_value('f');

    if (! fused_system_size.initialise(f))
    {
      cerr << "Invalid fused system size (-f) specification '" << f << "'\n";
      usage(2);
    }
  }

  if (cl.option_present('h'))
  {
    const_IWSubstring h = cl.string_value('h');

    if (! number_spinach_attachments.initialise(h))
    {
      cerr << "Invalid number spinach attachments specification (-h) '" << h << "'\n";
      usage(2);
    }
  }

  if (cl.option_present('j'))
  {
    const_IWSubstring j = cl.string_value('j');

    if (! number_non_spinach_attachments.initialise(j))
    {
      cerr << "Invalid number non spinach attachments specification (-j) '" << j << "'\n";
      usage(2);
    }
  }

  if (cl.option_present('H'))
  {
    const_IWSubstring h = cl.string_value('H');

    if (! spinach_size_specification.initialise(h))
    {
      cerr << "Invalid total spinach size specification (-H) '" << h << "'\n";
      usage(2);
    }
  }

  if (cl.option_present('c'))
  {
    const_IWSubstring c = cl.string_value('c');

    if (! number_connections.initialise(c))
    {
      cerr << "Invalid total ring substitution specification (-c) '" << c << "'\n";
      usage(2);
    }
  }

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if ("cut" == m)
      {
        cut_out_matching_rings = 1;

        if (verbose)
          cerr << "Matching rings will be cut out\n";
      }
      else if ("cutrs" == m)
      {
        cut_out_ring_and_spinach = 1;

        if (verbose)
          cerr << "Will cut out ring and spinach\n";
      }
      else if (m.starts_with("spch="))
      {
        m.remove_leading_chars(5);
        if (! m.numeric_value(isotope_for_spinach_connections) || isotope_for_spinach_connections <= 0)
        {
          cerr << "The isotope for spinach connections (spch=) must be a whole +ve number\n";
          usage(2);
        }

        if (verbose)
          cerr << "Will label spinach connections with " << isotope_for_spinach_connections << endl;
      }
      else if (m.starts_with("embd="))
      {
        m.remove_leading_chars(5);
        if (! m.numeric_value(isotope_for_embedding_atoms) || isotope_for_embedding_atoms <= 0)
        {
          cerr << "The isotope for embedding connections (embd=) must be a whole +ve number\n";
          usage(2);
        }

        if (verbose)
          cerr << "Will label embedding connections with " << isotope_for_embedding_atoms << endl;
      }
      else if (m.numeric_value(isotope_label) && isotope_label > 0)
      {
        if (verbose)
          cerr << "Will isotopically label matched rings with " << isotope_label << endl;
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        return 2;
      }
    }
  }
  else
    isotope_label = 1;

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('B'))
  {
    const const_IWSubstring b = cl.option_value('B');

    stream_for_molecules_not_matching.add_output_type(SMI);
    if (stream_for_molecules_not_matching.would_overwrite_input_files(cl, b))
    {
      cerr << "The -B file cannot overwrite any input file(s)\n";
      return 2;
    }

    if (! stream_for_molecules_not_matching.new_stem(b))
    {
      cerr << "Cannot initialise stream for non matching molecules '" << b << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Molecules not matching constraints written to '" << b << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ring_identification(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_matching << " molecules matched, wrote " << molecules_written << " molecules\n";

    if (unique_fragments_only)
      cerr << "Discarded " << duplicate_rings_discarded << " duplicate fragments, wrote " << already_seen.size() << " fragments\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ring_identification (argc, argv);

  return rc;
}
