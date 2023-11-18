/*
  Compute E/Z fingerprints
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Utilities/GFP_Tools/sparsefp.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_tested = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static  Atom_Typing_Specification atom_typing;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString tag("NCCT<");

static int function_as_filter = 0;

static extending_resizable_array<int> nbits;

static int ntest = 0;

static int keep_going_after_test_failure = 0;

static int test_failures = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes Cis/Trans fingerprints\n";
  cerr << "  -P <type>     atom typing to use for computing tie breaks\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -J <tag>      tag to use (default '" << tag << "'\n";
  cerr << "  -t <ntest>    perform tests\n";
  cerr << "  -k            keep going after test failure\n";
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
count_bonds (const Bond * b)
{
  if (b->is_aromatic())
    return 31;

  if (b->is_single_bond())
    return 1312;

  if (b->is_double_bond())
    return 290322;

  if (b->is_triple_bond())
    return 6092;

  return 9;    // HUH?
}

#define EXPAND_FROM1 1
#define NEXT_SHELL1 2
#define COMPLETED1 4
#define EXPAND_FROM2 8
#define NEXT_SHELL2 16
#define COMPLETED2 32

/*
  To avoid passing large numbers of arguments, we
  set up a small class
*/

class Resolve
{
  private:
    const atom_number_t _start;
    const atom_number_t _a1;
    const atom_number_t _a2;
    unsigned int _s1;
    unsigned int _s2;

//  private functions

    void _reset_for_next_expansion (int n, int * p, unsigned int, unsigned int, unsigned int) const;
    void _reset_for_next_expansion (int n, int * p) const;
    int _expand(const Molecule & m, int * p, const int * atype, unsigned int, unsigned int, unsigned int) const;

  public:
    Resolve (atom_number_t start, atom_number_t a1, atom_number_t a2);

    unsigned int s1 () const { return _s1;}
    unsigned int s2 () const { return _s2;}

    int compute_resolution (const Molecule & m, int * processed, const int *);
};

Resolve::Resolve (atom_number_t start,
                  atom_number_t a1,
                  atom_number_t a2) : _start(start), _a1(a1), _a2(a2)
{
  assert (INVALID_ATOM_NUMBER != _a1);

  return;
}

int
Resolve::_expand (const Molecule & m,
                  int * processed,
                  const int * atype,
                  unsigned int expand_from,
                  unsigned int next_shell,
                  unsigned int completed) const
{
  int ne = m.nedges();

  int rc = 0;

#ifdef DEBUG_SCAN
  cerr << "Scan " << ne << " edges, from " << _start << endl;
#endif

  for (int i = 0; i < ne; i++)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

#ifdef DEBUG_SCAN
    cerr << "Atoms " << a1 << " value " << processed[a1] << " and " << a2 << " value " << processed[a2] << endl;
#endif

    if (completed & processed[a1])
      ;
    else if (expand_from & processed[a1])   // great, expand from here
    {
      if (0 == (completed & processed[a2]) || (next_shell & processed[a2]))
      {
        rc += atype[a2] * count_bonds(b);
        processed[a2] |= next_shell;
      }
    }

    if (completed & processed[a2])
      ;
    else if (expand_from & processed[a2])    // great, expand from here
    {
      if (0 == (completed & processed[a1]) || (next_shell & processed[a1]))
      {
        rc += atype[a1] * count_bonds(b);
        processed[a1] |= next_shell;
      }
    }
  }

#ifdef DEBUG_SCAN
  for (int i = 0; i < m.natoms(); i++)
  {
    if (next_shell & processed[i])
      cerr << "Will expand from " << i << endl;
  }
#endif

  return rc;
}

void
Resolve::_reset_for_next_expansion (int n,
                                    int * processed,
                                    unsigned int expand_from,
                                    unsigned int next_shell,
                                    unsigned int completed) const
{
  for (int i = 0; i < n; i++)
  {
    if (0 == processed[i])
      continue;

    if (expand_from & processed[i])
    {
      processed[i] ^= expand_from;
      processed[i] |= completed;
    }
    else if (next_shell & processed[i])
    {
      processed[i] ^= next_shell;
      processed[i] |= expand_from;
    }
  }

  return;
}

void
Resolve::_reset_for_next_expansion (int n, int * processed) const
{
  _reset_for_next_expansion(n, processed, EXPAND_FROM1, NEXT_SHELL1, COMPLETED1);
  _reset_for_next_expansion(n, processed, EXPAND_FROM2, NEXT_SHELL2, COMPLETED2);

//#define DEBUG_RESOLUTION
#ifdef DEBUG_RESOLUTION
  for (int i = 0; i < n; i++)
  {
    if (0 != processed[i])
      cerr << "Atom " << i << " value " <<processed[i] <<endl;
  }
#endif

  return;
}

int
Resolve::compute_resolution (const Molecule & m,
                             int * processed,
                             const int * atype)
{
  set_vector(processed, m.natoms(), 0);

  processed[_start] = (COMPLETED1 | COMPLETED2);

  if (INVALID_ATOM_NUMBER == _a1)
    _s1 = 0;
  else 
  {
    _s1 = atype[_a1];
    processed[_a1] = EXPAND_FROM1;
  }

  if (INVALID_ATOM_NUMBER == _a2)
    _s2 = 0;
  else 
  {
    _s2 = atype[_a2];
    processed[_a2] = EXPAND_FROM2;
  }

#ifdef DEBUG_RESOLUTION
    cerr << "Initial score for " << _a1 << " is " << _s1 << endl;
    cerr << "Initial score for " << _a2 << " is " << _s2 << endl;
#endif

  while (_s1 == _s2)
  {
    int d1 = _expand(m, processed, atype, EXPAND_FROM1, NEXT_SHELL1, COMPLETED1);
    int d2 = _expand(m, processed, atype, EXPAND_FROM2, NEXT_SHELL2, COMPLETED2);

#ifdef DEBUG_RESOLUTION
    cerr << "After expansion  s1 " << _s1 << " d1 " << d1 << " s2 " << _s2 << " d2 " << d2 << endl;
#endif

    if (0 == d1 && 0 == d2)    // cannot expand, done, not resolved
      return 0;

    _s1 = _s1 * 1011 + d1;    // just some arbitrary number
    _s2 = _s2 * 1011 + d2;

    if (d1 != d2)    // resolved, done
      break;

    _reset_for_next_expansion(m.natoms(), processed);
  }

  if (_s1 < _s2)
    return -1;
  else if (_s1 > _s2)
    return 1;
  else
    return 0;
}

static int
identify_attached_atoms (const Molecule & m,
                         atom_number_t avoid,
                         atom_number_t start,
                         atom_number_t & a1,
                         atom_number_t & a2)
{
  const Atom * as = m.atomi(start);

  int scon = as->ncon();

  if (scon > 3)
    return 0;

  for (int i = 0; i < scon; i++)
  {
    const Bond * b = as->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(start);

    if (avoid == j)   // should never happen due to double bond check above
      continue;

    if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else if (INVALID_ATOM_NUMBER == a2)
      a2 = j;
    else         // should not happen
      return 0;
  }

  return (INVALID_ATOM_NUMBER != a1);
}

/*
  Looks like

  a3            a1
    \          /  
     \        /
      a1 == a2         <- the bond's a1 and a2
     /        \
    /          \
  a4            a2
*/

//#define DEBUG_EZ 1

static int
ez_fingerprint_v2 (Molecule & m,
                   const Bond & b,
                   int * processed,
                   Sparse_Fingerprint_Creator & sfp)
{
  set_vector(processed, m.natoms(), 0);

  int * atype = new int[m.natoms()]; std::unique_ptr<int[]> free_atype(atype);

  (void) atom_typing.assign_atom_types(m, atype);

  atom_number_t a1 = INVALID_ATOM_NUMBER;
  atom_number_t a2 = INVALID_ATOM_NUMBER;

  if (! identify_attached_atoms (m, b.a1(), b.a2(), a1, a2))    // first arg is avoid
    return 0;

#ifdef DEBUG_RESOLUTION
  cerr << "Line " << __LINE__ << " atoms " << a1 << " and " << a2 << endl;
#endif

  Resolve rhs (b.a2(), a1, a2);
  if (!  rhs.compute_resolution(m, processed, atype))
    return 0;

#ifdef DEBUG_RESOLUTION
  cerr << "Line " << __LINE__ << endl;
#endif

  atom_number_t a3 = INVALID_ATOM_NUMBER;
  atom_number_t a4 = INVALID_ATOM_NUMBER;

  if (! identify_attached_atoms (m, b.a2(), b.a1(), a3, a4)) 
    return 0;

#ifdef DEBUG_RESOLUTION
  cerr << "Line " << __LINE__ << " atoms " << a3 << " and " << a4 << endl;
#endif

  Resolve lhs (b.a1(), a3, a4);

  if (! lhs.compute_resolution(m, processed, atype))
    return 0;

#ifdef DEBUG_RESOLUTION
  cerr << "Line " << __LINE__ << endl;
#endif

// Great, we have each of the two ends resolved. Now we need
// to overlay that on the cis-trans arrangement. Try to find
// a V-like shape. A1 and A3 will NOT be INVALID_ATOM_NUMBER

  const Bond * b13 = m.bond_between_atoms (b.a1(), a3);
  const Bond * b21 = m.bond_between_atoms (b.a2(), a1);

// We need to figure out which atoms are at the 4 "corners"

  int lhsdirection;

  if (b13->is_directional_up())
  {
    if (b13->a1() == b.a1())
      lhsdirection =  1;
    else 
      lhsdirection = -1;
  }
  else if (b13->is_directional_down())
  {
    if (b13->a1() == b.a1())
      lhsdirection = -1;
    else
      lhsdirection =  1;
  }
  else
    return 0;

  int rhsdirection;

  if (b21->is_directional_up())
  {
    if (b21->a1() == b.a2())
      rhsdirection =  1;
    else
      rhsdirection =  -1;
  }
  else if (b21->is_directional_down())
  {
    if (b21->a1() == b.a2())
      rhsdirection = -1;
    else
      rhsdirection =  1;
  }
  else
    return 0;

  unsigned int nw, ne, se, sw;

  if (lhsdirection > 0)
  {
    nw = lhs.s1();
    sw = lhs.s2();
  }
  else
  {
    nw = lhs.s2();
    sw = lhs.s1();
  }

  if (rhsdirection > 0)
  {
    ne = rhs.s1();
    se = rhs.s2();
  }
  else
  {
    ne = rhs.s2();
    se = rhs.s1();
  }

#ifdef DEBUG_EZ
  cerr << "Initial assignment - random\n";
  cerr <<  nw << "\t\t" << ne << endl;
  cerr <<  sw << "\t\t" << se << endl;
#endif

// Now we need to impose a canonical arrangement of the vertices

  atom_number_t aa1 = b.a1();
  atom_number_t aa2 = b.a2();

  unsigned int t;

  if (nw + sw + (nw * sw) + 5232 * atype[b.a1()] < ne + se + (ne * se) + 5232 * atype[b.a2()])   // swap across N-S
  {
    t = nw;
    nw = ne;
    ne = t;
    t = sw;
    sw = se;
    se = t;
    std::swap(aa1, aa2);
  }
  
  if (nw + ne + (nw * ne) < sw + se + (sw * se))    // swap across W-E
  {
    t = nw;
    nw = sw;
    sw = t;
    t = ne;
    ne = se;
    se = t;
  }

  if (nw > sw)
  {
    t = nw;
    nw = sw;
    sw = t;
    t = ne;
    ne = se;
    se = t;
  }

#ifdef DEBUG_EZ
  cerr <<  nw << "\t\t" << ne << endl;
  cerr <<  sw << "\t\t" << se << endl;
  cerr << "atomic numbers " << m.atomic_number(aa1) << " and " << m.atomic_number(aa2) << endl;
#endif

  unsigned int bb = nw + 23 * atype[aa1] + 79 * sw +
                    367 * atype[aa2] + 4243 * ne + 7927 * se;

#ifdef DEBUG_EZ
  cerr << "Bit " << bb << endl;
#endif

  sfp.hit_bit(bb);

  return 1;
}

static int
ez_fingerprint_v2 (Molecule & m,
                   Sparse_Fingerprint_Creator & sfp)
{
#ifdef DEBUG_EZ
  cerr << "Processing '" << m.name() << "'\n";
#endif

  int * processed = new int[m.natoms()]; std::unique_ptr<int[]> free_processed(processed);

  int ne = m.nedges();

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < ne; i++)
  {
    const Bond * b = m.bondi(i);

    if (! b->part_of_cis_trans_grouping())
      continue;

    if (m.is_ring_atom(b->a1()) && m.is_ring_atom(b->a2()))  // too hard
      continue;

#ifdef DEBUG_EZ
    cerr << "Bond between " << b->a1() << " and " << b->a2() << endl;
#endif
    if (! ez_fingerprint_v2(m, *b, processed, sfp))
      return 0;
  }

#ifdef DEBUG_EZ
  cerr << "Final form has " << sfp.nbits() << " bits\n";
#endif

  return 1;
}

static int
ez_fingerprint_v2 (Molecule & m,
                   IWString_and_File_Descriptor & output)
{
  if (! function_as_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfp;

  (void) ez_fingerprint_v2(m, sfp);

  if (verbose)
    nbits[sfp.nbits()]++;

  IWString tmp;

  sfp.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << '\n';

  if (! function_as_filter)
    output << "|\n";

  output.write_if_buffer_holds_more_than(16388);

  return output.good ();
}

static int
ez_fingerprint_v2 (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ez_fingerprint_v2(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
ez_fingerprint_v2_tester (Molecule & m)
{
  Sparse_Fingerprint_Creator sfp;
  if (! ez_fingerprint_v2(m, sfp))
  {
    if (verbose > 1)
      cerr << "Nothing to test in '" << m.name() << "'\n";
    return 1;
  }

  if (0 == sfp.nbits())   // nothing to test
  {
    if (verbose > 1)
      cerr << "No bits set in '" << m.name() << "'\n";
    return 1;
  }

  IWString initial_smiles = m.smiles();

  molecules_tested++;

  for (int i = 0; i < ntest; i++)
  {
    IWString s = m.random_smiles();

    Molecule mcopy;
    mcopy.build_from_smiles(s);

    Sparse_Fingerprint_Creator sfpc;
    ez_fingerprint_v2(mcopy, sfpc);

    if (sfpc != sfp)
    {
      cerr << "Test failure '" << m.name() << "'\n";
      cerr << initial_smiles << endl;
      cerr << s << endl;
      sfp.debug_print(cerr);
      sfpc.debug_print(cerr);
      return 0;
    }
  }

  return 1;
}

static int
ez_fingerprint_v2_tester (data_source_and_type<Molecule> & input)
{
  Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ez_fingerprint_v2_tester(*m))
    {
      cerr << "Test failure processing '" << m->name() << "'\n";
      test_failures++;
      if (! keep_going_after_test_failure)
        return 0;
    }
  }

  return 1;
}

static int
ez_fingerprint_v2 (const char * fname, FileType input_type, 
                IWString_and_File_Descriptor & output)
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

  if (ntest)
    return ez_fingerprint_v2_tester(input);

  return ez_fingerprint_v2(input, output);
}

static int
ez_fingerprint_v2_filter_record (const_IWSubstring buffer,  // local copy
                          IWString_and_File_Descriptor & output)
{
  buffer.remove_leading_chars(smiles_tag.length());
  buffer.chop();

  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Cannot interpret smiles '" << buffer << "'\n";
    return 0;
  }

  return ez_fingerprint_v2(m, output);
}

static int
ez_fingerprint_v2_filter (iwstring_data_source & input,
                          IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with(smiles_tag))
      continue;

    if (! ez_fingerprint_v2_filter_record (buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
ez_fingerprint_v2_filter (IWString_and_File_Descriptor & output)
{
  iwstring_data_source input("-");
  if (! input.good())
  {
    cerr << "Huh, ez_fingerprint_v2_filter cannot open stdin\n";
    return 0;
  }

  return ez_fingerprint_v2_filter(input, output);
}

static int
ez_fingerprint_v2 (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lfJ:t:kP:");

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

  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing.build(p))
    {
      cerr << "Cannot discern atom typing specification '" << p << "'\n";
      return 4;
    }
  }
  else
    atom_typing.set_atom_type(IWATTYPE_TT);

  if (cl.option_present('f'))
  {
    function_as_filter = 1;

    if (verbose)
      cerr << "Will work as a TDT filter\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', tag);

    if (! tag.ends_with ('<'))
      tag << '<';

    if (verbose)
      cerr << "Will write fingerprints with tag '" << tag << "'\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', ntest) || ntest < 1)
    {
      cerr << "The number of tests to perform (-t) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will perform " << ntest << " tests per molecule\n";

    if (cl.option_present('k'))
    {
      keep_going_after_test_failure = 1;

      if (verbose)
        cerr << "Will keep going after test failures\n";
    }
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
  else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))
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

  if (function_as_filter)
  {
    if (! ez_fingerprint_v2_filter(output))
      rc = 6;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! ez_fingerprint_v2(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < nbits.number_elements(); i++)
    {
      if (nbits[i])
        cerr << nbits[i] << " molecules had " << i << " cis/trans bonds\n";
    }

    if (ntest)
      cerr << "Tested " << molecules_tested << " molecules, " <<  test_failures << " test failures\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ez_fingerprint_v2 (argc, argv);

  return rc;
}
