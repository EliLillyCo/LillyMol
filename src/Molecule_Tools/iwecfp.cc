/*
 Implementation of extended connectivity
*/

#include <limits>

#include <assert.h>
#include <memory>
#include <unordered_map>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "iwaray.h"

#include "cmdline.h"
#include "sparse_fp_creator.h"
#include "accumulator.h"
#include "misc.h"
#include "iw_auto_array.h"
#include "iwdigits.h"

#include "molecule.h"
#include "iwstandard.h"
#include "aromatic.h"
#include "istream_and_type.h"
#include "atom_typing.h"

static Chemical_Standardisation chemical_standardisation;

static int verbose=0;

static int molecules_read = 0;

static int add_tails = 0;

static int reduce_to_largest_fragment = 0;    // Feb 2010, allow multi fragment molecules

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString tag;

static Atom_Typing_Specification atom_type;

static int function_as_filter = 0;

static int min_shell_radius = 0;
static int max_shell_radius = 0;

static int additive = 1;

static int only_set_bits_for_max_radius_shell = 0;

static Accumulator_Int<int> nbits_acc;

static int each_shell_gets_own_fingerprint = 0;

static Sparse_Fingerprint_Creator * global_sparse_fingerprint = NULL;

static int all_bonds_same_type = 0;

static IWDigits iwdigits_center;    // first on line, no leading blank
static IWDigits iwdigits;

/*
  There are two ways we can update the global fingerprint.
  We can count all the bits present in each molecule, or we can
  just record presence of the bit in yet another molecule
*/

static int update_global_fingerprint_presence_only = 0;

/*
  Sometimes we want to determine the structural features that give rise
  to certain bits.
*/

static std::unordered_map<unsigned int, unsigned int> bits_to_investigate;
static int looking_for_bit_meanings = 0;
static int bits_found = 0;

static IWString_and_File_Descriptor stream_for_bit_meanings;

static IWString_and_File_Descriptor stream_for_all_bits;

/*
  Jan 2017. Play with the idea of generating extra pairs to equalise times an atom is used
*/

static int equalise_atom_coverage = 0;

static int label_by_visited = 0;

typedef unsigned int atype_t;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Compute the Extended Connectivity fingerprints for molecules\n";
  cerr << "  -r <len>       min shell width for writing a fingerprint\n";
  cerr << "  -R <length>    set the maximum step for the connected shell\n";
  cerr << "  -m             multiplicative formation of bits\n";
  cerr << "  -s             each radius gets its own fingerprint\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -J <tag>       set the tag for the name tag of the fingerprints\n";
  cerr << "  -P z           atom types are atomic numbers (default)\n";
  cerr << "  -P complex     atom types are a convolution of type and environment\n";
  cerr << "  -P tt          atom types are topological torsion types\n";
  cerr << "  -P bs          atom types are basic types\n";
  cerr << "  -f             filter existing TDT/fingerprint file\n";
  cerr << "  -G <fname>     write the global fingerprint\n";
  cerr << "  -M READ=<fname>  identify structural features making bits in <fname>\n";
  cerr << "  -M WRITE=<fname> write identified features to <fname>\n";
  cerr << "  -B <fname>     write info on all bits produced to <fname> (large!)\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -x             only set bits for max radius shell\n";
  cerr << "  -b             all bond type info lost\n";
  cerr << "  -Q ...         options controlling equalisation of times atoms are fingerprinted\n";
  //  (void) display_standard_aromaticity_options (cerr);
  //  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  //  (void) display_standard_sparse_fingerprint_options (cerr, 'F');
  cerr << "  -E ...         standard element options\n";
  cerr << "  -v             verbose output\n";
  
  exit(rc);
}

static int
read_bits_to_investigate(iwstring_data_source & input,
                         std::unordered_map<unsigned int, unsigned int> & bits_to_investigate)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;
    buffer.truncate_at_first(' ');

    unsigned int b;
    if (! buffer.numeric_value(b))
    {
      cerr << "Invalid bit number '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    bits_to_investigate[b] = 0;
  }

  return static_cast<int>(bits_to_investigate.size());
}

static int
read_bits_to_investigate(const char * fname,
                          std::unordered_map<unsigned int, unsigned int> & bits_to_investigate)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open bits to be examined file name '" << fname << "'\n";
    return 0;
  }

  return read_bits_to_investigate(input, bits_to_investigate);
}

static int
check_against_list(Molecule & m,
                   const IWString & smarts,
                   const atom_number_t centre_of_shell,
                   unsigned int sum_so_far,
                   const int radius)
{
  std::unordered_map<unsigned int, unsigned int>::iterator f = bits_to_investigate.find(sum_so_far);

  if (f == bits_to_investigate.end())
    return 0;

  bits_found++;

  (*f).second++;

  if (radius > 0)
  m.set_atom_map_number(centre_of_shell, radius);   //  does not work for zero radius

  stream_for_bit_meanings << m.smiles() << ' ' << m.name() << " bit " << sum_so_far << " atom " << centre_of_shell << " " << smarts << " radius " << radius << '\n';

  if (radius > 0)
    m.set_atom_map_number(centre_of_shell, 0);

  stream_for_bit_meanings.write_if_buffer_holds_more_than(32768);

  return 1;
}

#define USE_IWDIGITS
#ifdef USE_IWDIGITS
static void
write_bit(const int centre_of_shell, 
           const IWString & smarts_for_centre_of_shell,
           const int radius,
           unsigned int b, 
           IWString_and_File_Descriptor & output)
{
  iwdigits_center.append_number(output, centre_of_shell);
  iwdigits.append_number(output, radius);
  iwdigits.append_number(output, b);    // caching will hardly ever be useful

  output << ' ' << smarts_for_centre_of_shell << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return;
}
#else
static void
write_bit(const int centre_of_shell, 
           const IWString & smarts_for_centre_of_shell,
           const int radius,
           unsigned int b, 
           IWString_and_File_Descriptor & output)
{
  char sep = ' ';

  output << centre_of_shell << sep << radius << sep << b << sep << smarts_for_centre_of_shell << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return;
}
#endif

static int
write_global_fingerprint(const Sparse_Fingerprint_Creator::FPHash & bits_found,
                          int radius,
                          int molecules_read,
                          IWString_and_File_Descriptor & output)
{
  IWString string_radius;

  string_radius << radius << ' ';

  for (Sparse_Fingerprint_Creator::FPHash::const_iterator i = bits_found.begin(); i != bits_found.end(); ++i)
  {
    unsigned int b = (*i).first;
    int c = (*i).second;

    float ave = static_cast<float>(static_cast<double>(c) / static_cast<double>(molecules_read));

    output << string_radius << b << ' ' << (*i).second << ' ' << ave << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static int
preprocess_molecule(Molecule & m)
{
  m.remove_all(1);
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();  // always reduce to largest fragment
  
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return 1;
}

static int 
bond_constant(const Bond * bondi)
{
  if (all_bonds_same_type)
    return 1;

  if (bondi->is_aromatic())
    return 11;
  if (bondi->is_triple_bond())
    return 7;
  if (bondi->is_double_bond())
    return 5;
  
  return 3;
}

static void
increment(unsigned int & sum_so_far,
          const int bc,
          const atype_t atom_constant)
{
  if (additive)
    sum_so_far += bc * atom_constant;
  else
    sum_so_far *= bc * atom_constant;

  return;
}

/*
  For the mode in which we look for bits, we need to know the 
  current molecule. Ideally, this would be passed as
  an argument, but the argument list is already long enough, so
  we do a very bad thing and just make it a file scope static
  variable
*/

static Molecule * current_molecule;
static atom_number_t centre_of_shell = INVALID_ATOM_NUMBER;
static IWString smarts_for_centre_of_shell;

#define PROCESSING_FINISHED 1
#define READY_TO_PROCESS 2
#define NEXT_TIME 3

static int
generate_shells(const int matoms,
                int radius,
                const int max_radius,
                const Atom * const * a,
                const atype_t * atom_constant,
                int * processing_status,
                unsigned int sum_so_far,
                Sparse_Fingerprint_Creator * sfc)
{
  radius++;

  if (additive)
    sum_so_far *= 7879;   // an arbitrary prime number

//#define DEBUG_ECFP_BITS 1
#ifdef DEBUG_ECFP_BITS
  cerr << "On entry sum_so_far " << sum_so_far << " radius " << radius << endl;
#endif

// Check for tail addition outside the loop

  int add_tails_here;
  if (add_tails > 0 && radius <= add_tails)
    add_tails_here = 1;
  else
    add_tails_here = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (READY_TO_PROCESS != processing_status[i])
      continue;

    const Atom * ai = a[i];

    int acon = ai->ncon();

    const Bond * const * bonds = ai->rawdata();    // for efficiency

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = bonds[j];

      atom_number_t k = b->other(i);

      if (PROCESSING_FINISHED == processing_status[k])   // we are extending the shell
      {
        int bc = bond_constant(b);
//      cerr << "     BC " << bc << " atom " << i << " atype " << atom_constant[i] << " begin " << sum_so_far << " extra " << (bc * atom_constant[i]) << " expect " << (sum_so_far + bc*atom_constant[i]) << endl;
        increment(sum_so_far, bc, atom_constant[i]);
//      cerr << "at level " << radius << " added atom " << i << " atype " << atom_constant[i] << " BC " << bc << " sum " << sum_so_far << endl;
        if (add_tails_here)
          sfc->hit_bit(sum_so_far);
      }
      else if (READY_TO_PROCESS == processing_status[k])
        ;
      else
        processing_status[k] = NEXT_TIME;
    }
  }

  if (radius > min_shell_radius)
  {
#ifdef DEBUG_ECFP_BITS
    cerr << "Hit bit " << sum_so_far << " at radius " << radius << endl;
#endif
    if (! only_set_bits_for_max_radius_shell || radius == max_radius)
      sfc->hit_bit(sum_so_far);

    if (looking_for_bit_meanings)
      check_against_list(*current_molecule, smarts_for_centre_of_shell, centre_of_shell, sum_so_far, radius);   // horrible hack with file scope variable
    else if (stream_for_all_bits.is_open())
      write_bit(centre_of_shell, smarts_for_centre_of_shell, radius, sum_so_far, stream_for_all_bits);
  }

  if (max_radius > 0 && radius >= max_radius)
    return 1;

  int continue_processing = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (READY_TO_PROCESS == processing_status[i])
      processing_status[i] = PROCESSING_FINISHED;
    else if (NEXT_TIME == processing_status[i])
    {
      processing_status[i] = READY_TO_PROCESS;
      continue_processing = 1;
    }
  }

  if (! continue_processing)
    return 1;

  if (each_shell_gets_own_fingerprint)
    return generate_shells(matoms, radius, max_radius, a, atom_constant, processing_status, sum_so_far, sfc + 1);
  else
    return generate_shells(matoms, radius, max_radius, a, atom_constant, processing_status, sum_so_far, sfc);
}

static int
add_fingerprint(const Sparse_Fingerprint_Creator & sfrom,
                 Sparse_Fingerprint_Creator & sto)
{
  const Sparse_Fingerprint_Creator::FPHash & bits_found = sfrom.bits_found();

  if (update_global_fingerprint_presence_only)
  {
    for (Sparse_Fingerprint_Creator::FPHash::const_iterator i = bits_found.begin(); i != bits_found.end(); ++i)
    {
      unsigned int b = (*i).first;

      sto.hit_bit(b);
    }
  }
  else
  {
    for (Sparse_Fingerprint_Creator::FPHash::const_iterator i = bits_found.begin(); i != bits_found.end(); ++i)
    {
      unsigned int b = (*i).first;
      int c = (*i).second;

      sto.hit_bit(b, c);
    }
  }

  return 1;
}

static int
update_global_sparse_fingerprint(const Sparse_Fingerprint_Creator * sfc)
{
  if (each_shell_gets_own_fingerprint)
  {
    for (int i = 0; i < max_shell_radius; i++)
    {
      add_fingerprint(sfc[i], global_sparse_fingerprint[i]);
    }
  }
  else
    add_fingerprint(sfc[0], global_sparse_fingerprint[0]);

  return 1;
}

static int
write_array_of_fingerprints(Sparse_Fingerprint_Creator * sfc,
                            IWString_and_File_Descriptor & output)
{
  IWString tmp;

  if (each_shell_gets_own_fingerprint)
  {
//  cerr << "Writing " << max_shell_radius << " fingerprints\n";
    for (int i = 0; i <= max_shell_radius; i++)
    {
      IWString tag_this_fp(tag);
      tag_this_fp << i;

      tmp.resize_keep_storage(0);
//    cerr << "Fingerprint " << i << " gets tag '" << tag_this_fp << "'\n";
      sfc[i].daylight_ascii_form_with_counts_encoded(tag_this_fp, tmp);
      output << tmp << '\n';
    }
  }
  else
  {
    sfc[0].daylight_ascii_form_with_counts_encoded(tag, tmp);
    output << tmp << '\n';
  }

  return output.size();
}

static int
do_output (Molecule &m,
           Sparse_Fingerprint_Creator * sfc, 
           IWString_and_File_Descriptor & output)
{
  if (verbose)
    nbits_acc.extra(sfc[0].nbits());

  if (function_as_filter)
    write_array_of_fingerprints(sfc, output);
  else
  {
    output << smiles_tag << m.smiles() << ">\n";
    write_array_of_fingerprints(sfc, output);
    output << identifier_tag << m.name() << ">\n";
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(32768);

//sfc->debug_print(cerr);

  if (NULL != global_sparse_fingerprint)
    update_global_sparse_fingerprint(sfc);

  return output.good();
}

static void
form_bit(Molecule & m,
         const atype_t * atom_constant,
         const Atom * const * atoms,
         const atom_number_t zatom,
         int max_r,
         int * processing_status,
         Sparse_Fingerprint_Creator * sfc)
{
  const int matoms = m.natoms();

  std::fill_n(processing_status, matoms, 0);

  processing_status[zatom] = PROCESSING_FINISHED;

  const auto e = atom_constant[zatom];

  sfc[0].hit_bit(e);

  const Atom * a = atoms[zatom];

  const int acon = a->ncon();

  if (0 == max_r)
    return;

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const auto j = b->other(zatom);

    processing_status[j] = READY_TO_PROCESS;
  }

  if (each_shell_gets_own_fingerprint)
    generate_shells(matoms, 0, max_shell_radius, atoms, atom_constant, processing_status, e, sfc + 1);
  else
    generate_shells(matoms, 0, max_shell_radius, atoms, atom_constant, processing_status, e, sfc);

  return;
}

static void
identify_atoms_within_range(Molecule & m,
                            Set_of_Atoms * atoms_within_range)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)    // separate loop to add the central atom first
  {
    Set_of_Atoms * ai = atoms_within_range + i * (max_shell_radius + 1);

    for (int r = 0; r <= max_shell_radius; ++r)    // add central atom to each shell
    {
      ai[r].add(i);
    }
  }

  for (int i = 0; i < matoms; ++i)
  {
    Set_of_Atoms * ai = atoms_within_range + i * (max_shell_radius + 1);

    for (int j = i+1; j < matoms; ++j)
    {
      const int d = m.bonds_between(i, j);

      if (d > max_shell_radius)
        continue;

      Set_of_Atoms * aj = atoms_within_range + j * (max_shell_radius+1);

      for (int r = 1; r <= max_shell_radius; ++r)
      {
        if (d <= r)
        {
          ai[r].add(j);
          aj[r].add(i);
        }
      }
    }
  }

#ifdef DEBUG_IDENTIFY_ATOMS_WITHIN_RANGE
  for (int i = 0; i < matoms; ++i)
  {
    for (int r = 0; r <= max_shell_radius; ++r)
    {
      cerr << " atom " << i << " rad " << r << " atms " << atoms_within_range[i * (max_shell_radius+1)+r] << endl;
    }
  }
#endif
}

class Atom_Radius_Variance
{
  private:
    atom_number_t _a;
    int _radius;
    int _var;

  public:
    Atom_Radius_Variance() { _var = 0;}

    void set_atom(const atom_number_t s) { _a = s;}
    void set_radius (const int s) { _radius = s;}
    void set_var (const double s) { _var = s;}

    int var() const { return _var;}
};

/*
  We do not really compute the variance, since that would be a floating point number. We want an int.
  So, what we compute is

  matoms * sumv - (sumv * sumv)

  Divide by matoms to get a more standard form
*/

static int
our_own_custom_variance(const int matoms,
                        const int sumv,
                        const int sumv2)
{
  const int rc = matoms * sumv2 - (sumv * sumv);

  if (rc >= 0)
    return rc;

  cerr << "Negative variance: matoms " << matoms << " sumv " << sumv << " sumv2 " << sumv2 << endl;

  return std::numeric_limits<int>::max();
}

static int
compute_variance(const Set_of_Atoms & s, 
                 const int * visited,
                 const int matoms,
                 int sumv,
                 int sumv2)
{
  const int n = s.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = s[i];

    const int v = visited[j];

    sumv += 1;
    sumv2 = sumv2 - (v*v) + (v+1)*(v+1);    // take away the previous value, and add the value for visited one extra time - could simplify.....
  }

  return our_own_custom_variance(matoms, sumv, sumv2);
}

static void
recompute_sums(const int * visited,
               const int matoms,
               int & sumv,
               int & sumv2)
{
  sumv = 0;
  sumv2 = 0;

  for (int i = 0; i < matoms; ++i)
  {
    sumv += visited[i];
    sumv2 += visited[i] * visited[i];
  }

  return;
}

static int
do_equalise_atom_coverage(Molecule & m,
                          const Atom * const * atoms,
                          const atype_t * atom_constant,
                          int * processing_status,
                          Sparse_Fingerprint_Creator * sfc)
{
  const int matoms = m.natoms();

  Set_of_Atoms * atoms_within_range = new Set_of_Atoms[matoms * (max_shell_radius+1)]; std::unique_ptr<Set_of_Atoms[]> free_atoms_within_range(atoms_within_range);

  identify_atoms_within_range(m, atoms_within_range);

// Work out initial visited values

  int * visited = new_int(matoms, 1); std::unique_ptr<int[]> free_visited(visited);   // note initialised with 1 (zero radius shell)

  for (int i = 0; i < matoms; ++i)
  {
    for (int j = 0; j < matoms; ++j)
    {
      if (j == i)
        continue;

      const int d = m.bonds_between(i, j);

      if (d > max_shell_radius)
        continue;

      visited[j] += max_shell_radius - d + 1;
    }
  }

  int max_visited = visited[0];
  int min_visited = visited[0];

  for (int i = 0; i < matoms; ++i)
  {
    if (visited[i] > max_visited)
      max_visited = visited[i];
    else if (visited[i] < min_visited)
      min_visited = visited[i];
  }

  if ((static_cast<float>(max_visited - min_visited) / static_cast<float>(max_visited)) > 0.8f)   // close enough
    return 1;

// partial sums needed for variance computation

  int sumv = 0;
  int sumv2 = 0;

  const int narv = matoms * (max_shell_radius+1);

  int * var = new int[narv]; std::unique_ptr<int[]> free_var(var);

  recompute_sums(visited, matoms, sumv, sumv2);

  if (label_by_visited && verbose > 1)
  {
    for (int i = 0; i < matoms; ++i)
    {
      m.set_atom_map_number(i, visited[i]);
    }
    cerr << m.smiles() << ' ' << m.name() << " before equalisation, var " << our_own_custom_variance(matoms, sumv, sumv2) << endl;
  }

  for (int i = 0; i < equalise_atom_coverage; ++i)
  {
    for (int j = 0; j < matoms; ++j)
    {
      for (int r = 0; r <= max_shell_radius; ++r)
      {
        const Set_of_Atoms & s = atoms_within_range[j * (max_shell_radius+1) + r];
        var[j * (max_shell_radius+1) + r] = compute_variance(s, visited, matoms, sumv, sumv2);
      }
    }

    int min_variance = var[0];
    resizable_array<int> mins;

    for (int j = 1; j < narv; ++j)
    {
      if (0 == j % (max_shell_radius+1))   // skip radius 0 things
        continue;

      if (var[j] > min_variance)    // definitely not of interest
        continue;

      if (var[j] < min_variance)    // a new lowest value
      {
        mins.resize_keep_storage(0);
        mins.add(j);
        min_variance = var[j];
      }
      else
        mins.add(j);
    }

//  cerr << "Iteration " << i << " nmin " << mins.number_elements() << " var " << min_variance << endl;

    for (int j = 0; j < mins.number_elements(); ++j)
    {
      const Set_of_Atoms & s = atoms_within_range[mins[j]];
//    cerr << "    atoms " << s << endl;
      s.increment_vector(visited, 1);
      form_bit(m, atom_constant, atoms, s[0], s.number_elements() - 1, processing_status, sfc);
    }

    recompute_sums(visited, matoms, sumv, sumv2);
  }

  if (label_by_visited)
  {
    for (int i = 0; i < matoms; ++i)
    {
      m.set_atom_map_number(i, visited[i]);
    }

    if (verbose > 1)
      cerr << m.smiles() << ' ' << m.name() << " after equalisation, var " << our_own_custom_variance(matoms, sumv, sumv2) << endl;
  }

  return 1;
}

/* 
*/

static int
iwecfp (Molecule &m, 
        IWString_and_File_Descriptor & output)

{
  m.compute_aromaticity_if_needed();

//cerr << "Processing " << m.unique_smiles() << "'\n";
  int matoms = m.natoms();

  Atom ** atoms = new Atom * [matoms]; iw_auto_array<Atom *> free_atoms(atoms);

  m.atoms ((const Atom **) atoms);   // disregard of const OK

  atype_t * atom_constant = new atype_t[matoms]; std::unique_ptr<atype_t[]> free_atom_constant(atom_constant);

//#define CHECK_ATOM_TYPES_STUFF
#ifdef CHECK_ATOM_TYPES_STUFF
  Atom_Typing_Specification a1, a2;
  a1.build("Z");
  a2.build("UST:z");
  a2.swap_atomic_number_atom_type_to_atomic_number_prime();

  a1.assign_atom_types(m, atom_constant, NULL);
  iw_write_array(atom_constant, matoms, "C", cerr);
  a2.assign_atom_types(m, atom_constant, NULL);
  iw_write_array(atom_constant, matoms, "HPAC", cerr);
#endif

  if (! atom_type.assign_atom_types(m, atom_constant, NULL))
  {
    cerr << "Cannot assign atom types '" << m.name() << "'\n";
    return 0;
  }

  int set_centre_atom_global_variable = 0;
  if (looking_for_bit_meanings || stream_for_all_bits.is_open())
  {
    set_centre_atom_global_variable = 1;
    current_molecule = &m;
    if (stream_for_all_bits.is_open())
      stream_for_all_bits << m.name() << '\n';
  }

  int * processing_status = new int[matoms]; iw_auto_array<int> free_processing_status(processing_status);

  Sparse_Fingerprint_Creator * sfc;
  if (each_shell_gets_own_fingerprint)
    sfc = new Sparse_Fingerprint_Creator[max_shell_radius + 1];
  else
    sfc = new Sparse_Fingerprint_Creator;

  for (int i = 0; i < matoms; i++)
  {
    const auto e = atom_constant[i];

#ifdef DEBUG_ECFP_BITS
    if (0 == min_shell_radius)
      cerr << "Starting with atom " << i << " bit " << e << endl;
#endif
    if (set_centre_atom_global_variable)
    {
      smarts_for_centre_of_shell = m.smarts_equivalent_for_atom(i);
      centre_of_shell = i;
    }

    if (0 == min_shell_radius)
    {
      if (! only_set_bits_for_max_radius_shell && 0 != max_shell_radius)
        sfc[0].hit_bit(e);

      if (looking_for_bit_meanings)
        check_against_list(m, smarts_for_centre_of_shell, i, e, 0);
      else if (stream_for_all_bits.is_open())
        write_bit(i, smarts_for_centre_of_shell, 0, e, stream_for_all_bits);
    }

    std::fill_n(processing_status, matoms, 0);

    processing_status[i] = PROCESSING_FINISHED;

    const Atom * ai = atoms[i];

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = ai->other(i, j);

      processing_status[k] = READY_TO_PROCESS;
    }

    if (each_shell_gets_own_fingerprint)
      generate_shells(matoms, 0, max_shell_radius, atoms, atom_constant, processing_status, e, sfc + 1);
    else
      generate_shells(matoms, 0, max_shell_radius, atoms, atom_constant, processing_status, e, sfc);
  }

#ifdef DEBUG_ECFP_BITS
  cerr << "Produced " << sfc[0].nbits() << " bits\n";
#endif

  if (equalise_atom_coverage)
    do_equalise_atom_coverage(m, atoms, atom_constant, processing_status, sfc);

  int rc = do_output (m, sfc, output);

  if (each_shell_gets_own_fingerprint)
    delete [] sfc;
  else
    delete sfc;

  if (stream_for_all_bits.is_open())
    stream_for_all_bits << "|\n";

  return rc;
}

static int
iwecfp(data_source_and_type<Molecule> & input, 
       IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (verbose > 1)
      cerr << molecules_read << " processing '" << m->name() << "'\n";

    if (! preprocess_molecule(*m))
    {
      cerr << "Skipping non organic or too large '" << m->name() << "'\n";
      continue;
    }

    if (! iwecfp(*m, output))
    {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
iwecfp_filter(const const_IWSubstring & buffer, 
           IWString_and_File_Descriptor & output)
{
  assert (buffer.ends_with('>'));
  const_IWSubstring smiles(buffer);

  smiles.remove_up_to_first('<');
  smiles.chop();
  
  Molecule m;
  
  if (! m.build_from_smiles(smiles))
  {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  preprocess_molecule(m);
  
  molecules_read++;
  
  return iwecfp(m, output);
}

static int
iwecfp_filter(iwstring_data_source & input, 
           IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  
  while (input.next_record(buffer))
  {
    output << buffer << '\n';
    
    if (! buffer.starts_with(smiles_tag))
      continue;
    
    if (! iwecfp_filter(buffer, output))
      {
        cerr << "Fatal error on line " << input.lines_read() << endl;
        cerr << buffer << endl;
        return 0;
      }
  }

  return output.good();
}

static int
iwecfp_filter(const char * fname, 
           IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.good())
    {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }

  return iwecfp_filter(input, output);
}

static int
iwecfp(const char * fname,
         int input_type,
         IWString_and_File_Descriptor & output)
{
  if (function_as_filter)
    return iwecfp_filter(fname, output);

  data_source_and_type<Molecule> input(input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return iwecfp(input, output);

}

static void
display_dash_G_options(std::ostream & os)
{
  os << "  -G presence   record presence of bits only, not count\n";
  os << "  -G WRITE=fname write the global fingerprint to <fname>\n";

  exit(0);
}

static int
iwecfp(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vE:A:g:J:i:fr:R:P:mt:lxsG:M:bB:Q:");
  
  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');
  
  if (! process_elements(cl))
  {
    usage(2);
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_shell_radius) || min_shell_radius < 0)
    {
      cerr << "The min shell radius (-r) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will only fingerprint paths larger than " << min_shell_radius << " bonds\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_shell_radius) || max_shell_radius < 0)
    {
      cerr << "The max shell radius (-R) must be a whole +ve number\n";
      usage(2);
    }

    if (max_shell_radius < min_shell_radius)
    {
      cerr << "Inconsistent min " << min_shell_radius << " and max " << max_shell_radius << " shell radius\n";
      usage(6);
    }

    if (verbose)
      cerr << "Max radius " << max_shell_radius << endl;
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', add_tails) || add_tails < 1)
    {
      cerr << "The add tails radius (-t) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
     cerr << "Will fingerprint tails out to radius " << add_tails << endl;
  }

  if (cl.option_present('s'))
  {
    if (! cl.option_present('R'))
    {
      cerr << "Sorry, the -s option only works when the -R option is specified\n";
      usage(4);
    }

    each_shell_gets_own_fingerprint = 1;

    if (verbose)
      cerr << "Each radius will get its own fingerprint\n";
  }

  if (cl.option_present('x'))
  {
    only_set_bits_for_max_radius_shell = 1;
    if (verbose)
      cerr << "Only bits for the largest radius will be set\n";
  }

  if (cl.option_present('b'))
  {
    all_bonds_same_type = 1;

    if (verbose)
      cerr << "All bonds considered identical\n";
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      usage(14);
    }
  }

  if (cl.option_present('m'))
  {
    additive = 0;

    if (verbose)
      cerr << "Fingerprints formed with multiplication operations\n";
  }

// If neither atom type nor tag specified, take a default

  if (! cl.option_present('J') && ! cl.option_present('P'))
  {
    atom_type.set_atom_type(IWATTYPE_COMPLEX);
    tag = "NCECC<";
  }
  else if (! cl.option_present('P') && cl.option_present('J'))
  {
    cerr << "Must specify atom typing to use with the -P option\n";
    usage(3);
  }
  else      // -P present, and maybe also -J
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_type.build(p))
    {
      cerr << "Cannot discern atom type '" << p << "'\n";
      return 3;
    }

    atom_type.swap_atomic_number_atom_type_to_atomic_number_prime();

    if (cl.option_present('J'))
    {
      tag = cl.string_value('J');
    }
    else
    {
      tag = "NCEC";
      atom_type.append_to_tag(tag);
    }

    if (! tag.ends_with('<'))
      tag += '<';

    if (verbose)
      cerr<< "Extended connectivity index written as non-colliding sparse fingerprints, tag '"<<tag <<"'\n";
  }
  
  if (each_shell_gets_own_fingerprint)
    tag.chop();

  int input_type = 0;

  if (cl.option_present('f'))
  {
    function_as_filter = 1;
  }
  else if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
      {
        cerr << "Cannot determine input type\n";
        usage(6);
      }
  }
  else if (! all_files_recognised_by_suffix(cl))
  {
    cerr << "Cannot determine input type(s)\n";
    return 7;
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }
  
  set_global_aromaticity_type(Daylight);
  set_input_aromatic_structures(1);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor stream_for_global_fingerprint;

  if (cl.option_present('G'))
  {
    const char * fname = NULL;
    const char * dbname = NULL;

    int i = 0;
    IWString g;
    while (cl.value('G', g, i++))
    {
      if ("presence" == g)
      {
        update_global_fingerprint_presence_only = 1;
        if (verbose)
          cerr << "Global fingerprint updated for presence only\n";
      }
      else if (g.starts_with("WRITE="))
      {
        g.remove_leading_chars(6);
        fname = g.c_str();
      }
      else if ("help" == g)
        display_dash_G_options(cerr);
      else
      {
        cerr << "Unrecognised -G qualifier '" << g << "'\n";
        usage(3);
      }
    }

    if (NULL == fname)
    {
      cerr << "No output specified for -G\n";
      return 4;
    }

    if (each_shell_gets_own_fingerprint)
      global_sparse_fingerprint = new Sparse_Fingerprint_Creator[max_shell_radius];
    else
      global_sparse_fingerprint = new Sparse_Fingerprint_Creator[1];

    if (NULL != fname)
    {
      if (! stream_for_global_fingerprint.open(fname))
      {
        cerr << "Cannot open stream for global fingerprint '" << fname << "'\n";
        return 3;
      }

      if (verbose)
        cerr << "Global fingerprint written to '" << fname << "'\n";
    }
    else if (NULL != dbname)
    {
    }
  }

  if (cl.option_present('M'))
  {
    IWString read_fname, write_fname;

    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (m.starts_with("READ="))
      {
        m.remove_leading_chars(5);
        read_fname = m;
      }
      else if (m.starts_with("WRITE="))
      {
        m.remove_leading_chars(6);
        write_fname = m;
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        usage(2);
      }
    }

    if (0 == read_fname.length() || 0 == write_fname.length())
    {
      cerr << "For identifying features, must have both read and write file names (-M)\n";
      usage(3);
    }

    if (read_fname == write_fname)
    {
      cerr << "Cannot use same file for reading and writing (-M)\n";
      usage(2);
    }

    if (! read_bits_to_investigate(read_fname.null_terminated_chars(), bits_to_investigate))
    {
      cerr << "Cannot read bits to investigate from '" << read_fname << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Read " << bits_to_investigate.size() << " bits to investigate from '" << read_fname << "'\n";

    if (! stream_for_bit_meanings.open(write_fname.null_terminated_chars()))
    {
      cerr << "Cannot open stream for bit interpretations '" << write_fname << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Bit interpretations written to '" << write_fname << "'\n";

    looking_for_bit_meanings = 1;
  }
  
  if (cl.option_present('B'))
  {
    const char * b = cl.option_value('B');

    if (! stream_for_all_bits.open(b))
    {
      cerr << "Cannot open stream for all bits '" << b << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Info on all bits written to '" << b << "'\n";

    iwdigits_center.initialise(100);

    iwdigits.set_include_leading_space(1);
    iwdigits.initialise(100);
  }

  if (cl.option_present('Q'))
  {
    const_IWSubstring q;
    for (int i = 0; cl.value('Q', q, i); ++i)
    {
      if ("label" == q)
      {
        label_by_visited = 1;
      }
      else if ("help" == q)
      {
      }
      else if (! q.numeric_value(equalise_atom_coverage) || equalise_atom_coverage < 1)
      {
        cerr << "The number of visited equalisation steps (-Q) must be a whole +ve number\n";
        usage(1);
      }
    }
  }

  int rc = 0;

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! iwecfp(cl[i], input_type, output))
        {
          rc = i + 1;
          break;
        }
    }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    if (nbits_acc.n() > 0)
      cerr << "Fingerprints had between " << nbits_acc.minval() << " and " << nbits_acc.maxval() << " ave " << static_cast<float>(nbits_acc.average_if_available_minval_if_not()) << " bits set\n";

    if (looking_for_bit_meanings)
    {
      cerr << "Found " << bits_found << " bits in lookup file\n";
      for (std::unordered_map<unsigned int, unsigned int>::const_iterator i = bits_to_investigate.begin(); i != bits_to_investigate.end(); ++i)
      {
        cerr << "Found " << (*i).second << " instances of " << (*i).first << '\n';
      }
    }
  }

  if (NULL != global_sparse_fingerprint)
  {
    int istop;
    if (each_shell_gets_own_fingerprint)
      istop = max_shell_radius;
    else
      istop = 1;

    for (int i = 0; i < istop; i++)
    {
      const Sparse_Fingerprint_Creator::FPHash & bits_found = global_sparse_fingerprint[i].bits_found();
      write_global_fingerprint(bits_found, i, molecules_read, stream_for_global_fingerprint);
    }

    stream_for_global_fingerprint.close();
  }
  
  return rc;
}

int
main(int argc, char ** argv)
{
  int rc = iwecfp(argc, argv);
  return rc;
}
