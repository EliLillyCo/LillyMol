/*
  For Ying Yang.
  has fragments with values computed on certain atoms - mostly heteroatoms
*/

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static char output_separator = ' ';

static int numeric_attributes = 1;

static int omit_unmatched_heteroatoms = 0;

class Three_Accumulators
{
  private:
    Accumulator<double> _acc[3];

  public:
    Three_Accumulators();

    void extra (const Atom * a);

    int n() const { return _acc[0].n();}

    template <typename T> int report (const int i, const char output_separator, T & output) const;
};

Three_Accumulators::Three_Accumulators()
{
  return;
}

void
Three_Accumulators::extra (const Atom * a) 
{
  _acc[0].extra(a->x());
  _acc[1].extra(a->y());
  _acc[2].extra(a->z());
}

template <typename T>
int
Three_Accumulators::report (const int i, const char output_separator, T & output) const
{
  output << _acc[i].n() << output_separator << _acc[i].minval() << output_separator << _acc[i].maxval();
  if (_acc[i].n() > 0)
    output << output_separator << static_cast<float>(_acc[i].average()) << '\n';
  else
    output << '\n';

  return 1;
}

/*
  Certain atoms within the query can have numeric values attached
*/

class Query_and_Values : public Substructure_Query
{
  private:
    Molecule _m;

    resizable_array<double> * _v;

    int * _bonds_to_attachment_point;

    int _number_attachment_points;

//  private functions

    int _extract_numeric_values();
    int _identify_bonds_to_attachment_point();

  public:
    Query_and_Values();
    ~Query_and_Values();

    int build (const const_IWSubstring & buffer);

    template <typename T> int process (const Molecule & m, Molecule_to_Match & target, int * already_matched, const int query_number, Three_Accumulators * acc, T & output);
    template <typename T> int _process_match (const Molecule & m, int * already_matched, const int query_number, Three_Accumulators * acc, const Set_of_Atoms & e, T & output) const;
};

Query_and_Values::Query_and_Values()
{
  _v = nullptr;

  _bonds_to_attachment_point = nullptr;

  return;
}

Query_and_Values::~Query_and_Values()
{
  if (nullptr != _v)
    delete [] _v;

  if (nullptr != _bonds_to_attachment_point)
    delete [] _bonds_to_attachment_point;
}

int
Query_and_Values::build (const const_IWSubstring & buffer)
{
  if (! _m.build_from_smiles(buffer))
  {
    cerr << "Query_and_Values::build:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  if (! _extract_numeric_values())
  {
    cerr << "Query_and_Values::build:cannot identify atom with numeric values\n";
    return 0;
  }

  if (! _identify_bonds_to_attachment_point())
  {
    cerr << "Query_and_Values::build:cannot identify bonds to attachment point\n";
    return 0;
  }

  Molecule_to_Query_Specifications mqs;
  mqs.set_substituents_only_at_isotopic_atoms(1);

  if (! Substructure_Query::create_from_molecule(_m, mqs))
  {
    cerr << "Query_and_Values::build:cannot build query from " << _m.smiles() << endl;
    return 0;
  }

  return 1;
}

int
Query_and_Values::_extract_numeric_values ()
{
  const int matoms = _m.natoms();

  assert (nullptr == _v);

  _v = new resizable_array<double>[matoms];
  _bonds_to_attachment_point = new_int(matoms, matoms + 1);

  int atoms_with_coordinates = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = _m.atomi(i);

    if (0.0f == a->x() && 0.0f == a->y() && 0.0f == a->z())
      continue;

    _v[i].add(a->x());
    _v[i].add(a->y());
    _v[i].add(a->z());

    atoms_with_coordinates++;
  }

  if (0 == atoms_with_coordinates)
  {
    cerr << "Query_and_Values::_extract_numeric_values:no atom with non-zero coordinates found\n";
    return 0;
  }

  return atoms_with_coordinates;
}

int
Query_and_Values::_identify_bonds_to_attachment_point()
{
  if (1 != _m.number_fragments())
  {
    cerr << "Query_and_Values::_identify_bonds_to_attachment_point:multi fragment molecule " << _m.smiles() << endl;
    return 0;
  }

  _number_attachment_points = 0;

  const int matoms = _m.natoms();

  Set_of_Atoms special_atoms;

  for (int i = 0; i < matoms; ++i)
  {
    if (_v[i].number_elements())
      special_atoms.add(i);
  }

  for (int i = 0; i < matoms; ++i)
  {
    if (_v[i].number_elements() > 0)    // values associated with this atom
      continue;

    if (0 == _m.isotope(i))
      continue;

    _number_attachment_points++;

    for (int j = 0; j < special_atoms.number_elements(); ++j)
    {
      const atom_number_t s = special_atoms[j];

      const  int b = _m.bonds_between(s, i);

      if (b < _bonds_to_attachment_point[s])
        _bonds_to_attachment_point[s] = b;
    }
  }

  if (0 == _number_attachment_points)
  {
    cerr << "Query_and_Values::_identify_bonds_to_attachment_point:no isotopic atoms found\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
Query_and_Values::process (const Molecule & m,
                           Molecule_to_Match & target,
                           int * already_matched,
                           const int query_number,
                           Three_Accumulators * acc,
                           T & output)
{
  Substructure_Results sresults;

  const int nhits = Substructure_Query::substructure_search(target, sresults);

  if (0 == nhits)
    return 0;

  for (int j = 0; j < nhits; ++j)
  {
    const Set_of_Atoms * e = sresults.embedding(j);
    _process_match(m, already_matched, query_number, acc, *e, output);
  }

  return nhits;
}

template <typename T>
int
Query_and_Values::_process_match (const Molecule & m,
                                  int * already_matched,
                                  const int query_number,
                                  Three_Accumulators * acc,
                                  const Set_of_Atoms & e,
                                  T & output) const
{
  const int n = _m.natoms();
  for (int i = 0; i < n; ++i)
  {
    if (0 == _v[i].number_elements())   // we want to scan the special atoms
      continue;

    const atom_number_t x = e[i];

    output << query_number << output_separator << i << output_separator << x << output_separator << m.smarts_equivalent_for_atom(x) << output_separator << _m.natoms() << output_separator << _bonds_to_attachment_point[i];

    acc[x].extra(_m.atomi(i));

    for (int j = 0; j < _v[i].number_elements(); ++j)
    {
      output << output_separator << _v[i][j];
    }
    output << '\n';

    already_matched[x]++;
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Overlapping fragment model - for Ying Yang\n";
  cerr << "  -F ...        isotopically labelled fragments - values are X and Y coordinates\n";
  cerr << "  -n <n>        number of numeric attributes to process (btw 1 and 3, default 1)\n";
  cerr << "  -b            omit unmatched heteroatoms from reports\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
read_queries (iwstring_data_source & input,
              resizable_array_p<Query_and_Values> & queries)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    Query_and_Values * q = new Query_and_Values;
    if (! q->build(buffer))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      delete q;
      return 0;
    }

    q->set_find_unique_embeddings_only(1);

    queries.add(q);
  }

  return queries.number_elements();
}

static int
read_queries (IWString & fname,
              resizable_array_p<Query_and_Values> & queries)
{
  iwstring_data_source input(fname.null_terminated_chars());

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_queries(input, queries);
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
overlapping_fragment_model (Molecule & m,
                            resizable_array_p<Query_and_Values> & queries,
                            IWString_and_File_Descriptor & output)
{
  const int matoms = m.natoms();
  int * already_matched = new_int(matoms); std::unique_ptr<int[]> free_already_matched(already_matched);
  Three_Accumulators * acc = new Three_Accumulators[matoms]; std::unique_ptr<Three_Accumulators[]> free_acc(acc);

  output << m.smiles() << output_separator << m.name() << output_separator << "PARENT\n";
  write_isotopically_labelled_smiles(m, false, output);
  output << output_separator << "AtomNo\n";

  Molecule_to_Match target(&m);

  const int nq = queries.number_elements();

  for (int i = 0; i < nq; ++i)
  {
    queries[i]->process(m, target, already_matched, i, acc, output);
  }

  m.set_isotopes(already_matched);

  int heteroatoms = 0;
  int heteroatoms_not_matched = 0;

  IWString coord("XYZ");

  for (int i = 0; i < matoms; ++i)
  {
    if (6 == m.atomic_number(i))
      continue;

    heteroatoms++;

    if (0 == already_matched[i])
      heteroatoms_not_matched++;

    for (int j = 0; j < 3; ++j)
    {
      if (0 == acc[i].n() && omit_unmatched_heteroatoms)
        continue;

      output << "AVE" << output_separator << i << output_separator << m.smarts_equivalent_for_atom(i) << output_separator << coord[j] << output_separator;
      acc[i].report(j, output_separator, output);
    }
  }

  output << m.smiles() << output_separator << m.name() << output_separator << "MATCHES" << output_separator << heteroatoms << output_separator << heteroatoms_not_matched << '\n';

//for (int i = 0; i < matoms; ++i)
//{
//}

  return output.good();
}

static int
overlapping_fragment_model (data_source_and_type<Molecule> & input,
                            resizable_array_p<Query_and_Values> & queries,
                            IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! overlapping_fragment_model(*m, queries, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
overlapping_fragment_model (const char * fname, FileType input_type, 
                            resizable_array_p<Query_and_Values> & queries,
                            IWString_and_File_Descriptor & output)
{
  assert(NULL != fname);

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

  return overlapping_fragment_model(input, queries, output);
}

static int
overlapping_fragment_model (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lF:o:n:s:b");

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

  if (cl.option_present('n'))
  {
    if (! cl.value('n', numeric_attributes) || numeric_attributes < 1 || numeric_attributes > 3)
    {
      cerr << "The number of numeric attributes must be a +ve whole number between 1 and 3\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will process " << numeric_attributes << " numeric attributes\n";
  }

  if (cl.option_present('b'))
  {
    omit_unmatched_heteroatoms = 1;

    if (verbose)
      cerr << "Will omit unmatched heteroatoms from output\n";
  }

  if (! cl.option_present('F'))
  {
    cerr << "Must specify fragments via the -F option\n";
    usage(1);
  }

  resizable_array_p<Query_and_Values> queries;

  if (cl.option_present('F'))
  {
    IWString f = cl.string_value('F');

    if (! read_queries(f, queries))
    {
      cerr << "Cannot read query fragments from '" << f << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << queries.number_elements() << " query fragments from '" << f << "'\n";
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

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! overlapping_fragment_model(cl[i], input_type, queries, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = overlapping_fragment_model(argc, argv);

  return rc;
}
