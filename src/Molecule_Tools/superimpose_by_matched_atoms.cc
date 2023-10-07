/*
  superimpose molecules by best fit to matched atoms
*/

#include <iostream>
#include <memory>
#include <limits>

using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define ISTREAM_AND_TYPE_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/target.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int report_rms_b4_and_after = 0;

static resizable_array_p<Substructure_Query> queries;

static Molecule_Output_Object stream_for_combined_molecules;

// Unimplemented feature
#ifdef MATCH_ATOMS_BY_ISOTOPES
static int match_atoms_by_isotopes = 0;
#endif

static resizable_array<int> matched_query_atoms;

static IWString_and_File_Descriptor stream_for_matched_atoms;

/*
  Apply weight different from 1.0 to matched atoms 1 and 2
  0-1.2-3
*/

static double special_processing_for_jibo = 0.0;

/*
  The -m option generates a list of matched query atoms. These will be picked up by
  each of the template molecules
*/

static Set_of_Atoms default_matched_query_atoms;

/*
*/

class Template_Molecule : public Molecule
{
  private:
    double * _a1;
    double * _a2;
    double * _weight;

    ::resizable_array_p<Set_of_Atoms>  _matched_query_atoms;

//  private functions

//  int _fill_a_array (const Set_of_Atoms & e);
    double _report_rms (int template_embedding_number, const Molecule & m, const Set_of_Atoms & m_embedding, std::ostream & output) const;
    double _compute_rms (int ndx, const Molecule & m, const Set_of_Atoms & m_embedding, double & maxd) const;

    int _write_matched_atoms (const int template_embedding, Molecule & m, const Set_of_Atoms & e, IWString_and_File_Descriptor & output) const;

    int _superimpose_by_matched_atoms (int template_embedding_number, Molecule & m, const Set_of_Atoms & e);
    int _superimpose_by_matched_atoms (Molecule & m, const Set_of_Atoms & e, Molecule_Output_Object & output);
    int _superimpose_by_matched_atoms (int template_embedding_number, Molecule & m,
                                                  const Set_of_Atoms & e, Molecule_Output_Object & output);

  public:
    Template_Molecule ();
    ~Template_Molecule ();

    int allocate_arrays(int size);

    int identify_query_matches (::resizable_array_p<Substructure_Query> &);

    int superimpose_by_matched_atoms (Molecule & m, const Set_of_Atoms & e, Molecule_Output_Object & output);

//  int report_rms (const Molecule & m, const Set_of_Atoms & m_embedding, ostream & output) const;
};

Template_Molecule::Template_Molecule ()
{
  _a1 = nullptr;
  _a2 = nullptr;
  _weight = nullptr;

  return;
}

Template_Molecule::~Template_Molecule ()
{
  if (nullptr != _a1)
  {
    delete [] _a1;
    delete [] _a2;
    delete [] _weight;
  }

  return;
}

int
Template_Molecule::allocate_arrays(int s)
{
  assert (nullptr == _a1);

  const auto matoms = this->natoms();

  _a1 = new double[matoms];
  _a2 = new double[matoms];
  _weight = new double[matoms];

  return 1;
}

static resizable_array_p<Template_Molecule> template_molecule;

static int ignore_molecules_not_matching_template_query = 0;

static int take_first_of_multiple_template_hits = 0;

static int ignore_queries_hitting_multiple_times = 0;

static int process_all_query_matches = 0;

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
  cerr << "  -T <fname>    file containing one template molecule\n";
  cerr << "  -s <smarts>   specify smarts for matched atoms\n";
  cerr << "  -M <smiles>   read the query as smiles - will be converted internally\n";
  cerr << "  -m n,n,n,n    specify which matched atoms to use\n";
  cerr << "  -u            only find unique substructure search embeddings\n";
#ifdef MATCH_ATOMS_BY_ISOTOPES
  cerr << "  -I            input molecules already have isotopic labels\n";
#endif
  cerr << "  -z i          ignore molecules not matching the query\n";
  cerr << "  -z f          take the first of multiple query matches\n";
  cerr << "  -z nom        ignore all matches if multiple matches\n";
  cerr << "  -z all        process all hits\n";
  cerr << "  -j <weight>   weight applied to matched atoms 1and2 0-1.2-3\n";
  cerr << "  -S <fname>    output file name\n";
  cerr << "  -o <type>     specify output type\n";
  cerr << "  -C <fname>    file for combined molecules\n";
  cerr << "  -Z ...        miscellaneous options\n";
  cerr << "  -r            produce RMS both before and after fitting\n";
  cerr << "  -D <fname>    write the matched atoms to <fname>\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

int
Template_Molecule::identify_query_matches (::resizable_array_p<Substructure_Query> & queries)
{
  Molecule_to_Match target(this);

  for (int i = 0; i < queries.number_elements(); i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
    {
      if (verbose)
        cerr << "Query " << i << " only matched " << sresults.max_query_atoms_matched_in_search() << " query atoms\n";
      continue;
    }

    if (verbose > 2)
      cerr << "Query " << i << " matched " << this->name() << " " << nhits << " times\n";

    if (1 == nhits)
      ;
    else if (take_first_of_multiple_template_hits)
      nhits = 1;
    else if (ignore_queries_hitting_multiple_times)
      continue;

    for (auto j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      Set_of_Atoms * a = new Set_of_Atoms;

      if (0 == default_matched_query_atoms.size())    // just take all the matched atoms
        *a = *e;
      else
      {
        for (unsigned int k = 0; k < default_matched_query_atoms.size(); k++)   // we only want some of the matched atoms
        {
          const auto l = default_matched_query_atoms[k];
          if (! e->ok_index(l))
          {
            cerr << "Template_Molecule::identify_query_matches:matched atom " << l << " invalid index into embedding " << (*e) << endl;
            return 0;
          }
          a->add(e->item(l));
        }
      }

      _matched_query_atoms.add(a);
    }

    int array_size = _matched_query_atoms[0]->size();

    _a1 = new double[array_size * 3];
    _a2 = new double[array_size * 3];
    _weight = new double[array_size];

    set_vector(_weight, array_size, 1.0);
    if (0.0 != special_processing_for_jibo)
    {
      _weight[1] = special_processing_for_jibo;
      _weight[2] = special_processing_for_jibo;
    }

    if (verbose)
      cerr << "Template molecule '" << Molecule::name() << "' " << nhits << " hits to substructure query\n";

    return 1;
  }

  cerr << "No queries matched '" << Molecule::name() << "'\n";

  return 0;
}

double
Template_Molecule::_report_rms (int template_embedding_number,
                                const Molecule & m,
                                const Set_of_Atoms & m_embedding,
                                std::ostream & output) const
{
  double maxd = 0.0;

  const double rms = _compute_rms(template_embedding_number, m, m_embedding, maxd);

  output << Molecule::name() << " to " << m.name() << ' ' << m_embedding.size() << " matched atoms, RMS " << rms << ", max " << maxd << "\n";

  return rms;
}

/*
  Compute RMS vs one of our embeddings
*/

double
Template_Molecule::_compute_rms (int template_embedding_number,
                                 const Molecule & m,
                                 const Set_of_Atoms & m_embedding,
                                 double & maxd) const
{
  const Set_of_Atoms & template_embedding = *(_matched_query_atoms[template_embedding_number]);

  const int n = m_embedding.number_elements();

  assert (n == static_cast<int>(template_embedding.size()));

  double sum = 0.0;

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t a1 = template_embedding[i];
    const atom_number_t a2 = m_embedding[i];

    const Atom * at1 = Molecule::atomi(a1);
    const Atom * at2 = m.atomi(a2);

    double d = at1->distance(*at2);

    if (d > maxd)
      maxd = d;

    sum += d * d;
  }

  return sqrt(sum / static_cast<double>(n));
}

static int
fill_a_array (const Molecule & m,
              const Set_of_Atoms & e,
              double * a)
{
  int n = e.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

    const Atom * ai = m.atomi(j);

    a[i * 3] = ai->x();
    a[i * 3 + 1] = ai->y();
    a[i * 3 + 2] = ai->z();
  }

  return 1;
}

extern "C" void u3b_(const double * w, double * c1, double * c2, const long int* n, const long int* mode, double *rms, double * u, double * t, long int* ier);

int
Template_Molecule::_superimpose_by_matched_atoms (int template_embedding_number,
                                                  Molecule & m,
                                                  const Set_of_Atoms & e)
{
  const Set_of_Atoms &  template_embedding = *(_matched_query_atoms[template_embedding_number]);

  fill_a_array(*this, template_embedding, _a1);
  fill_a_array(m, e, _a2);

  long int n = e.number_elements();

  long int mode = 1;
  double u[9];
  double rms;
  double t[3];
  long int ier = 0;

  u3b_(_weight, _a1, _a2, &n, &mode, &rms, u, t, &ier);

//if (verbose > 1)
//  cerr << "RMS fit for '" << m.name() << "' " << rms << endl;

  if (0 == ier)   // good
    ;
  else if (-1 == ier)
    cerr << "superposition not unique, but optimal\n";
  else
  {
    cerr << "u3b failed\n";
    return 0;
  }

  const int matoms = m.natoms();

  double rotmat11 = u[0];
  double rotmat12 = u[1];
  double rotmat13 = u[2];
  double rotmat21 = u[3];
  double rotmat22 = u[4];
  double rotmat23 = u[5];
  double rotmat31 = u[6];
  double rotmat32 = u[7];
  double rotmat33 = u[8];

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    double x0 = a->x() - t[0];
    double y0 = a->y() - t[1];
    double z0 = a->z() - t[2];

    double xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    double yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    double zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    m.setxyz(i, xx , yy , zz );
  }

  return 1;
}



/*int
Template_Molecule::_fill_a_array (const Set_of_Atoms & e)
{
  const int n = e.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

    const Atom * ai = this->atomi(j);

    _a1[i * 3] = ai->x();
    _a1[i * 3 + 1] = ai->y();
    _a1[i * 3 + 2] = ai->z();
  }

  return 1;
}*/

static void
write_atom (const char sep,
            const Molecule & m,
            const atom_number_t a,
            IWString_and_File_Descriptor & output)
{
  output << m.name() << sep << a << sep << m.smarts_equivalent_for_atom(a) << sep << m.x(a) << sep << m.y(a) << sep << m.z(a);

  return;
}

int
Template_Molecule::_write_matched_atoms (const int template_embedding_number,
                                         Molecule & m,
                                         const Set_of_Atoms & e,
                                         IWString_and_File_Descriptor & output) const
{ 
  const auto sep = ' ';

  output << "Template:" << sep << this->name() << sep << " embedding " << sep << template_embedding_number << sep << "molecule:" << sep << m.name() << '\n';

  const Set_of_Atoms * template_embedding = _matched_query_atoms[template_embedding_number];

  const auto n = e.size();

  double maxd = 0.0;        // we redo these computations
  double sum = 0.0;

  for (unsigned int i = 0; i < n; ++i)
  {
    auto j = template_embedding->item(i);

    write_atom(sep, *this, j, output);

    const auto at = Molecule::atomi(j);

    j = e[i];

    output << sep;

    write_atom(sep, m, j, output);

    const auto am = m.atomi(j);

    const auto d = at->distance(*am);

    if (d > maxd)
      maxd = d;

    sum += (d * d);

    output << '\n';
  }

  output << "RMS" << sep << static_cast<float>(sqrt(sum / static_cast<double>(n))) << sep << "maxd" << sep << maxd << '\n';

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Template_Molecule::_superimpose_by_matched_atoms (int template_embedding_number,
                                                  Molecule & m,
                                                  const Set_of_Atoms & e,
                                                  Molecule_Output_Object & output)
{
  if (stream_for_matched_atoms.is_open())
    _write_matched_atoms(template_embedding_number, m, e, stream_for_matched_atoms);

  if (report_rms_b4_and_after)
    _report_rms(template_embedding_number, m, e, cerr);

  if (! _superimpose_by_matched_atoms(template_embedding_number, m, e))
    return 0;

  const auto rms = _report_rms(template_embedding_number, m, e, cerr);

  if (stream_for_combined_molecules.active())
  {
    Molecule tmp(*this);
    tmp.add_molecule(&m);
    IWString mname;
    mname << Molecule::name() << ':' << m.name() << " RMS " << static_cast<float>(rms);
    tmp.set_name(mname);
    stream_for_combined_molecules.write(tmp);
  }

  IWString save_name(m.name());
  IWString tmp;
  tmp << m.name() << " RMS " << static_cast<float>(rms);
  const auto rc = output.write(m);
  m.set_name(save_name);
  return rc;
}

int
Template_Molecule::_superimpose_by_matched_atoms (Molecule & m,
                                                  const Set_of_Atoms & e,
                                                  Molecule_Output_Object & output)
{
  for (unsigned int i = 0; i < _matched_query_atoms.size(); ++i)
  {
    _superimpose_by_matched_atoms(i, m, e, output);
  }

  return 1;
}

/*
  The main function entry point.
  First task is to figure out whether we are using all the atoms in the embedding
  or just a subset.
*/

int
Template_Molecule::superimpose_by_matched_atoms (Molecule & m,
                                                 const Set_of_Atoms & m_embedding,
                                                 Molecule_Output_Object & output)

{
  const auto dmqa = default_matched_query_atoms.size();

  if (0 == dmqa)
    return _superimpose_by_matched_atoms(m, m_embedding, output);

  Set_of_Atoms s;
  s.resize(dmqa);

  for (unsigned int i = 0; i < dmqa; ++i)
  {
    const auto j = default_matched_query_atoms[i];
    s.add(m_embedding[j]);
  }

  return _superimpose_by_matched_atoms(m, s, output);
}

static int
superimpose_by_matched_atoms (Molecule & m,
                              Molecule_Output_Object & output)
{
  int rc = 0;

  for (int i = 0; i < queries.number_elements(); i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(m, sresults);

    if (0 == nhits)
      continue;

    if (1 == nhits)
      ;
    else if (take_first_of_multiple_template_hits)
      nhits = 1;
    else if (ignore_queries_hitting_multiple_times)
      continue;
    else if (process_all_query_matches)
      ;
    else
    {
      cerr << nhits << " hits to query " << i << " for molecule '" << m.name() << "', not sure what to do, skipping\n";
      continue;
    }

    if (verbose > 1)
      cerr << m.name() << " " << nhits << " hits to query " << i << endl;

    rc++;

    for (auto j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      for (unsigned int k = 0; k < template_molecule.size(); ++k)
      {
        template_molecule[k]->superimpose_by_matched_atoms(m, *e, output);
      }
    }
  }

  if (rc > 0)
    return rc;

  cerr << "None of " << queries.number_elements() << " queries matched '" << m.name() << "'\n";

  if (ignore_molecules_not_matching_template_query)
    return 1;

  return 0;
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
ok_dimensionality (const Molecule & m)
{
  if (3 == m.highest_coordinate_dimensionality())
    return 1;

  if (2 == m.highest_coordinate_dimensionality())
  {
    cerr << "Warning, only 2D coordinates in '" << m.name() << "'\n";
    return 1;
  }

  cerr << "NO coordinate info in '" << m.name() << "', cannot process\n";
  return 0;
}

static int
superimpose_by_matched_atoms(data_source_and_type<Molecule> & input,
                             Molecule_Output_Object & output)
{
  Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ok_dimensionality(*m))
      return 0;

    if (! superimpose_by_matched_atoms(*m, output))
    {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
superimpose_by_matched_atoms(const char * fname,
                             FileType input_type,
                             Molecule_Output_Object & output)
{

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if(! input.good())
  {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  return superimpose_by_matched_atoms(input, output);
}

static int
read_template_molecule(data_source_and_type<Template_Molecule> & input,
                       resizable_array_p<Template_Molecule> & template_molecule)
{
  Template_Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    preprocess(*m);

    if (3 != m->highest_coordinate_dimensionality())
    {
      cerr << "Template molecule '" << m->name() << "' is not 3D\n";
      delete m;
      return 0;
    }

    if (! m->identify_query_matches(queries))
    {
      cerr << "Cannot identify query matches in '" << m->name() << "'\n";
      delete m;
      return 0;
    }

    template_molecule.add(m);
  }

  return template_molecule.size();
}

static int
read_template_molecule(const const_IWSubstring & fname, 
                       resizable_array_p<Template_Molecule> & template_molecule,
                       FileType input_type)
{
  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }
  
  data_source_and_type<Template_Molecule> input(input_type, fname);

  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return read_template_molecule(input, template_molecule);
}

static void
display_misc_options (char flag, std::ostream & output)
{
  output << " -" << flag << " vH           the smarts v directive includes implicit Hydrogens\n";
  output << " -" << flag << " help         this message\n";

  exit(0);
}

static int
superimpose_molecules (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:o:g:lT:z:S:Iq:s:m:C:j:ryM:D:uUZ:");

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

    chemical_standardisation.deactivate_lactim_lactam();    // too troublesome
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('y'))
  {
    set_ignore_chirality_in_smarts_input(1);

    if (verbose)
      cerr << "Will ignore chirality in input queries\n";
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('z', m, i++))
    {
      if ("first" == m || 'f' == m)
      {
        take_first_of_multiple_template_hits = 1;
        if (verbose)
          cerr << "Will take the first of multiple hits for user query\n";
      }
      else if ("ignore" == m || "i" == m)
      {
        ignore_molecules_not_matching_template_query = 1;
        if (verbose)
          cerr << "Will ignore queries not matching\n";
      }
      else if ("nom" == m || "ignmmatch" == m)
      {
        ignore_queries_hitting_multiple_times = 1;
        if (verbose)
          cerr << "Will ignore multiple query matches in the user query\n";
      }
      else if ("all" == m)
      {
        process_all_query_matches = 1;

        if (verbose)
          cerr << "Will process all query matches\n";
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << m << "'\n";
        usage(17);
      }
    }
  }

  if (cl.option_present('j'))
  {
    if (! cl.value('j', special_processing_for_jibo) || special_processing_for_jibo <= 0.0)
    {
      cerr << "The weight for middle matched atoms (-j) must be a valid weight\n";
      usage(3);
    }

    if (verbose)
      cerr << "Weight for matched atoms 1 and 2 set to " << special_processing_for_jibo << endl;
  }

  if (cl.option_present('r'))
  {
    report_rms_b4_and_after = 1;

    if (verbose)
      cerr << "Will report RMS values before and after superimposition\n";
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
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('s'))
  {
    const_IWSubstring smarts;
    int i = 0;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Query * q = new Substructure_Query;
      if (! q->create_from_smarts(smarts))
      {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 62;
      }

      if (cl.option_present('u'))
        q->set_find_unique_embeddings_only(1);

      queries.add(q);
    }
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose, 'q'))
    {
      cerr << "Cannot read queries (-q)\n";
      return 2;
    }
  }

  if (verbose)
    cerr << "Defined " << queries.number_elements() << " queries\n";

  if (cl.option_present('U'))
  {
    for (int i = 0; i < queries.number_elements(); ++i)
    {
      queries[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('Z'))
  {
    const_IWSubstring z;
    for (int i = 0; cl.value('Z', z, i); ++i)
    {
      if ("vH" == z)
      {
        set_global_setting_nbonds_includes_implicit_hydrogens(1);
      }
      else if ("help" == z)
      {
        display_misc_options('Z', cerr);
      }
      else
      {
        cerr << "Unrecognised -Z qualifier '" << z << "'\n";
        display_misc_options('Z', cerr);
      }
    }
  }

  if (cl.option_present('M'))
  {
    const_IWSubstring smiles;
    int i = 0;
    while (cl.value('M', smiles, i++))
    {
      Molecule m;
      if (! m.build_from_smiles(smiles))
      {
        cerr << "Cannot interpret smiles '" << smiles << "'\n";
        return 2;
      }

      if (cl.option_present('y'))
        m.remove_all_chiral_centres();

      if (verbose > 1)
        cerr << "Building query from '" << m.smiles() << "'\n";

      if (chemical_standardisation.active())
        chemical_standardisation.process(m);

      Molecule_to_Query_Specifications mqs;
      mqs.set_ignore_molecular_hydrogen_information(1);
      mqs.set_make_embedding(1);

      Substructure_Query * q = new Substructure_Query;

      q->create_from_molecule(m, mqs);

      if (cl.option_present('u'))
        q->set_find_unique_embeddings_only(1);

      queries.add(q);
    }
  }

  if (cl.option_present('m'))
  {
    const_IWSubstring m = cl.string_value('m');

    int i = 0;
    const_IWSubstring token;
    while (m.nextword(token, i, ','))
    {
      int j;
      if (! token.numeric_value(j) || j < 0)
      {
        cerr << "Matched atom numbers must be whole non-negative numbers\n";
        usage(4);
      }

      default_matched_query_atoms.add(j);
    }
  }

#ifdef MATCH_ATOMS_BY_ISOTOPES
  if (cl.option_present('I') && queries.number_elements())
  {
    cerr << "Cannot do alignment by isotopes (-I) and alignment by queries (-s)\n";
    usage(3);
  }

  if (cl.option_present('I'))
  {
    match_atoms_by_isotopes = 1;

    if (verbose)
      cerr << "Will match atoms by isotopic labels\n";
  }
#endif

  if (cl.option_present('T'))
  {
    const_IWSubstring t = cl.string_value('T');

    if (! read_template_molecule(t, template_molecule, input_type))
    {
      cerr << "Cannot read template molecule from '" << t << "'\n";
      return 6;
    }

    const auto ntemplate = template_molecule.size();

    if (verbose)
      cerr << ntemplate <<  " template molecule(s) and queries built from '" << t << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('C'))
  {
    const_IWSubstring c = cl.string_value('C');
    stream_for_combined_molecules.add_output_type(FILE_TYPE_SDF);

    if (! stream_for_combined_molecules.new_stem(c))
    {
      cerr << "Cannot open stream for combined files '" << c << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Combined molecules written to '" << c << ".sdf\n";
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');

    if (! stream_for_matched_atoms.open(d))
    {
      cerr << "Cannot open stream for matched atoms '" << d << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Matched atoms written to '" << d << "'\n";
  }

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SDF);
  else if (! output.determine_output_types(cl))
  {
    cerr << "Cannot determine output type(s)\n";
    return 3;
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) with the -S option\n";
      return 8;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot initialise -S output stem '" << s << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Output written to '" << s << "'\n";
  }
  else 
    output.new_stem("-");

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! superimpose_by_matched_atoms(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }


  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = superimpose_molecules(argc, argv);

  return rc;
}
