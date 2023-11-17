/*
  Produce subsets of fused rings
*/

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <limits>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1
//#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define IWQSORT_IMPLEMENTATION 1

#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
using std::cerr;
using std::endl;
using std::ostream;
using std::unique_ptr;
using std::numeric_limits;

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_written = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int min_nrings = 1;

static int max_nrings = numeric_limits<int>::max();

static int write_parent_molecule = 0;

static int unique_rings_only = 0;

static int remove_chirality = 0;

static int preserve_aromaticity = 0;

static extending_resizable_array<int> ring_subsets;

static IW_STL_Hash_Map_int global_counter;

static int accumulate_global_counter = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces fused subsets of fused ring systems\n";
  cerr << "  -p            write parent molecules\n";
  cerr << "  -r <n>        min number of rings in the output\n";
  cerr << "  -R <n>        max number of rings in the output\n";
  cerr << "  -a            subsets must preserve aromaticity from the parent\n";
  cerr << "  -u            use unique smiles to discard duplicate structures\n";
  cerr << "  -G <fname>    accumulate count of all ring subsets found\n";
  cerr << "  -c            discard chirality\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -o <type>     specify output type\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

class IWBitsHash
{
  private:
  public:

#if defined (IW_INTEL_COMPILER)

    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;
    bool  operator () (const IW_Bits_Base &, const IW_Bits_Base &) const;

#endif

    size_t operator () (const IW_Bits_Base &) const;
};

typedef std::unordered_set<IW_Bits_Base, IWBitsHash> Bits_Hash;

size_t
IWBitsHash::operator () (const IW_Bits_Base & b) const
{
  const unsigned char * bits = b.bits();

  const unsigned int * u = reinterpret_cast<const unsigned int *>(bits);

  return static_cast<size_t>(u[0]);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_chirality)
  {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

class Atoms_Hit
{
  private:
    static int _natoms;

    int * _hit;

    int _nrings;

  public:
    Atoms_Hit(int n);
    Atoms_Hit(const Atoms_Hit &);
    ~Atoms_Hit();

    Atoms_Hit & operator = (const Atoms_Hit &);

    void operator += (const Ring &);
    void operator -= (const Ring &);
    int operator == (const Atoms_Hit &) const;

    int nrings() const { return _nrings;}

    void bits (IW_Bits_Base & b) const;
};

Atoms_Hit::Atoms_Hit(int n)
{
  _natoms = n;
  _nrings = 0;
  _hit = new_int(_natoms);
}

Atoms_Hit::~Atoms_Hit()
{
  delete [] _hit;
}

Atoms_Hit::Atoms_Hit (const Atoms_Hit & rhs)
{
  _hit = new int[_natoms];

  _nrings = rhs._nrings;

  copy_vector(_hit, rhs._hit, _natoms);
}

Atoms_Hit &
Atoms_Hit::operator = (const Atoms_Hit & rhs)
{
  copy_vector(_hit, rhs._hit, _natoms);

  _nrings = rhs._nrings;

  return *this;
}

void
Atoms_Hit::operator += (const Ring & r)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    int j = r[i];

    _hit[j]++;
  }

  _nrings++;

  return;
}

void
Atoms_Hit::operator -= (const Ring & r)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    int j = r[i];

    _hit[j]--;
  }

  _nrings--;

  return;
}

int
Atoms_Hit::operator == (const Atoms_Hit & rhs) const
{
  if (_nrings != rhs._nrings)
    return 0;

  for (int i = 0; i < _natoms; i++)
  {
    if (_hit[i] && rhs._hit[i])
      ;
    else if (0 == _hit[i] && 0 == rhs._hit[i])
      ;
    else
      return 0;
  }

  return 1;
}

void
Atoms_Hit::bits (IW_Bits_Base & b) const
{
  b.clear();

  b.construct_from_array_of_ints(_hit, _natoms);

  return;
}

int Atoms_Hit::_natoms;

static int
aromaticity_lost (Molecule & parent,
                  const int * xref,
                  Molecule & subset)
{
  int matoms = parent.natoms();

  for (int i = 0; i < matoms; i++)
  {
    int j = xref[i];

    if (j < 0)    // not in subset
      continue;

    if (! parent.is_aromatic(i))
      continue;

    if (! subset.is_aromatic(j))    // bad, aromaticity lost
      return 1;
  }

  return 0;    // no aromaticity lost
}

class Atom_Count_Comparator
{
  private:
  public:
    int operator () (Molecule *, Molecule *) const;
};

int
Atom_Count_Comparator::operator() (Molecule * m1,
                                   Molecule * m2) const
{
  int matoms1 = m1->natoms();
  int matoms2 = m2->natoms();

  if (matoms1 < matoms2)
    return -1;

  if (matoms1 > matoms2)
    return 1;

  int a1 = m1->aromatic_atom_count();
  int a2 = m2->aromatic_atom_count();

  if (a1 < a2)
    return -1;

  if (a1 > a2)
    return 1;

  int h1 = m1->natoms(6);
  int h2 = m2->natoms(6);

  if (h1 < h2)
    return -1;

  if (h1 > h2)
    return 1;

  return 0;
}

static Atom_Count_Comparator acc;

template void resizable_array_base<Molecule*>::iwqsort<Atom_Count_Comparator>(Atom_Count_Comparator&);

static int
form_subset (Molecule & m,
             int ndx,
             const IW_Bits_Base & b,
             Molecule & subset,
             int * tmp,
             int * xref,
             IW_STL_Hash_Set & usmiles)
{
  set_vector(tmp, m.natoms(), 0);

#ifdef DEBUG_DO_OUTPUT
  cerr << "Bit vector creating subset " << b.nset() << " atoms, molecule contains " << m.natoms() << " atoms\n";
  b.printon(cerr);
  cerr << endl;
#endif

  b.set_vector(tmp);

#ifdef DEBUG_DO_OUTPUT
  for (int i = 0; i < m.natoms(); i++)
  {
    if (tmp[i])
      cerr << "Includes atom " << i << endl;
  }
#endif

  m.create_subset(subset, tmp, 1, xref);

  if (subset.nrings() > max_nrings)
    return 0;

  if (subset.nrings() < min_nrings)
    return 0;

  if (preserve_aromaticity && aromaticity_lost (m, xref, subset))
    return 0;

  if (accumulate_global_counter)
    global_counter[subset.unique_smiles()]++;

  if (unique_rings_only)
  {
    if (usmiles.contains(subset.unique_smiles()))
      return 0;
    else
      usmiles.insert(subset.unique_smiles());

    subset.invalidate_smiles();
  }

  return 1;
}

static int
do_output (Molecule & m,
           Bits_Hash & seen,
           Molecule_Output_Object & output)
{
  if (write_parent_molecule)
  {
    output.write(m);
    molecules_written++;
  }

  IW_STL_Hash_Set usmiles;

  int * tmp = new int[m.natoms()]; std::unique_ptr<int[]> free_tmp(tmp);
  int * xref = new int[m.natoms()]; std::unique_ptr<int[]> free_xref(xref);

  resizable_array_p<Molecule> formed;

  int ndx = 1;
  for (Bits_Hash::const_iterator i = seen.begin(); i != seen.end(); ++i)
  {
    Molecule * subset = new Molecule;

    if (! form_subset (m, ndx, *i, *subset, tmp, xref, usmiles))
    {
      delete subset;
      continue;
    }

    formed.add(subset);

    ndx++;
  }

  formed.iwqsort(acc);

  int n = formed.number_elements();

  if (verbose)
    ring_subsets[n]++;

  for (int i = 0; i < n; i++)
  {
    Molecule * s = formed[i];

    IWString mname(m.name());
    mname << '_' << i;
    s->set_name(mname);

    output.write(*s);
  }

  molecules_written += n;

  return 1;
}

#ifdef NOT_BEING_USED
static void
set_bits (const Ring & r,
          IW_Bits_Base & b)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    b.set(j, 1);
  }

  return;
}
#endif

static int
already_formed (const Atoms_Hit & atoms_hit,
                Bits_Hash & seen)
{
  IW_Bits_Base b;
  atoms_hit.bits(b);

#ifdef DEBUG_ALREADY_FORMED
  cerr << "New ring(s), " << b.nset() << " atoms set\n"; b.printon(cerr);
  cerr << endl;
#endif

  if (seen.find(b) != seen.end())
    return 1;
  
  seen.insert(b);
  return 0;
}

static int
expand_ring (Molecule & m,
             const Ring & r,
             int nrings,
             int * ring_already_done,
             Bits_Hash & seen,
             Atoms_Hit & atoms_hit)
{
  int n = r.fused_ring_neighbours();

  for (int i = 0; i < n; i++)
  {
    const Ring * ri = r.fused_neighbour(i);

    int rn = ri->ring_number();

    if (ring_already_done[rn])
      continue;

    atoms_hit += *ri;

    if (already_formed (atoms_hit, seen))
      continue;

    ring_already_done[rn] = 1;

    if (nrings < max_nrings)
      expand_ring(m, *ri, nrings + 1, ring_already_done, seen, atoms_hit);

    ring_already_done[rn] = 0;
    atoms_hit -= *ri;
  }

  return 1;
}

static int
ring_combinations (Molecule & m,
                   Bits_Hash & seen,
                   int nbits,
                   int * ring_already_done)
{
  int nr = m.nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    int rn = ri->ring_number();

    if (ring_already_done[rn])
      continue;

    ring_already_done[rn] = 1;

    Atoms_Hit atoms_hit(m.natoms());
    atoms_hit += *ri;

    (void) already_formed(atoms_hit, seen);

    if (max_nrings > 1)
      expand_ring (m, *ri, 2, ring_already_done, seen, atoms_hit);

    ring_already_done[rn] = 0;
  }

  return seen.size();
}

static int
ring_combinations (Molecule & m,
                Molecule_Output_Object & output)
{
  int nr = m.nrings();

  if (nr < 2)
    return 1;

  int matoms = m.natoms();

  int nbits;
  if (0 == matoms % 32)
    nbits = matoms;
  else
    nbits = (matoms / 32 + 1) * 32;

  Bits_Hash seen;

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_rad(ring_already_done);

  if (0 == ring_combinations (m, seen, nbits, ring_already_done))
    return 0;

  if (0 == seen.size())   // all isolated rings maybe
    return 1;

  return do_output (m, seen, output);
}

static int
ring_combinations (data_source_and_type<Molecule> & input,
                Molecule_Output_Object & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ring_combinations(*m, output))
      return 0;
  }

  return 1;
}

static int
ring_combinations (const char * fname, FileType input_type, 
                Molecule_Output_Object & output)
{
  assert (NULL != fname);

  if (input_type == FILE_TYPE_INVALID)
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

  return ring_combinations(input, output);
}
static int
ring_combinations (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lr:R:o:pcuaG:");

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
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_nrings) || min_nrings < 1)
    {
      cerr << "The min rings in output (-r) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Output must contain at least " << min_nrings << " rings\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_nrings) || max_nrings < min_nrings)
    {
      cerr << "The max rings option (-R) must be a whole +ve number >= " << min_nrings << endl;
      usage(4);
    }

    if (verbose)
      cerr << "Will produce rings with at most " << max_nrings << " rings\n";
  }

  if (cl.option_present('p'))
  {
    write_parent_molecule = 1;

    if (verbose)
      cerr << "Will write parent molecules\n";
  }

  if (cl.option_present('u'))
  {
    unique_rings_only = 1;

    if (verbose)
      cerr << "Will suppress ring/ring combinations with duplicate smiles\n";
  }

  if (cl.option_present('a'))
  {
    preserve_aromaticity = 1;

    if (verbose)
      cerr << "Subsets must preserve all parent aromaticity\n";
  }

  if (cl.option_present('c'))
  {
    remove_chirality = 1;

    if (verbose)
      cerr << "Will remove chirality from input molecules\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor gfile;

  if (cl.option_present('G'))
  {
    accumulate_global_counter = 1;

    const char * fname = cl.option_value('G');

    if (! gfile.open(fname))
    {
      cerr << "Cannot open file for global counter '" << fname << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Global counters written to '" << fname << "'\n";
  }

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SMI);
  else if (! output.determine_output_types(cl, 'o'))
  {
    cerr << "Cannot determine output types (-o)\n";
    usage(3);
  }

  if (! output.new_stem("-"))
  {
    cerr << "Yipes, cannot open stdout\n";
    return 3;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ring_combinations(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << "Wrote " << molecules_written << " molecules\n";

    for (int i = 0; i < ring_subsets.number_elements(); i++)
    {
      if (ring_subsets[i])
        cerr << ring_subsets[i] << " molecules had " << i << " subsets\n";
    }
  }

  if (gfile.is_open())
  {
    int highest_count = 1;

    for (IW_STL_Hash_Map_int::const_iterator i = global_counter.begin(); i != global_counter.end(); ++i)
    {
      const IWString & s = (*i).first;
      int c = (*i).second;

      gfile << s << ' ' << c << '\n';
      gfile.write_if_buffer_holds_more_than(32768);

      if (c > highest_count)
        highest_count = c;
    }

    gfile.close();

    cerr << global_counter.size() << " unique ring combinations, as many as " << highest_count << " instances\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ring_combinations (argc, argv);

  return rc;
}
