/*
  We want to fingerprint just a substructure - useful for
  grouping molecules
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "circular_fingerprint.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static Chemical_Standardisation chemical_standardisation;

static Element_Transformations element_transformations;

static int reduce_to_largest_fragment = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static IWString fingerprint_tag;

static IWString empty_fingerprint_with_closing_angle_bracket_and_newline;

static int molecules_processed = 0;

static int break_at_first_match = 0;

static int fingerprint_each_substructure_match = 0;

static Molecule_Output_Object stream_for_labelled_molecules;

static int extend_shell = 0;

static IWString input_smiles_tag("$SMI<");
static IWString output_smiles_tag("SBSMI<");
static IWString tag_for_isotopically_labelled_parent;

static IWString identifier_tag("PCN<");

static IWString natoms_tag;

static int produce_atom_pair_fingerprint = 0;

static int ignore_molecules_not_matching_any_queries = 0;

static int discard_molecules_not_matching_any_queries = 0;

static int molecules_not_matching_queries = 0;

static int join_disconected_sections = 0;

static Atom_Typing_Specification atom_typing_specification;

static extending_resizable_array<int> disconnected_groups;

static int truncate_counted_fingerprints_to_one = 0;

static Circular_Fingerprint_Generator circular_fingerprint_generator;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Fingerprints just a subset of the atoms in a molecule\n";
  cerr << "  -q ...        query specifications for identifying the subset\n";
  cerr << "  -s <smarts>   smarts to identify the subset\n";
  cerr << "  -x <number>   include atoms within <number> atoms of the subset\n";
  cerr << "  -z i          ignore molecules not hitting any queries\n";
  cerr << "  -z d          \n";
  cerr << "  -z e          fingerprint each substructure match\n";
  cerr << "  -z f          if multiple matches, take the first\n";
  cerr << "  -S in=TAG     tag for reading smiles\n";
  cerr << "  -S out=TAG    tag for writing smiles\n";
  cerr << "  -S iso=TAG    tag for isotopically labelled molecules in output\n";
  cerr << "  -N <tag>      tag for the number of atoms in the subset\n";
  cerr << "  -I <stem>     file name stem for isotopically labelled subset molecules\n";
  cerr << "  -Y ...        standard fingerprint options\n";
  cerr << "  -j <n>        join disconnected sections <n> or closer bonds\n";
  cerr << "  -J <tag>      tag to use for fingerprints\n";
  cerr << "  -P ...        atom typing specification, enter '-P help' for info\n";
  cerr << "  -M            produce atom pair fingerprints\n";
  cerr << "  -f            work as a filter\n";
  cerr << "  -T ...         standard element transformations -T I=Cl -T Br=Cl ...\n";
  cerr << "  -e            truncate any counted fingerprint to 0/1\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
write_sparse_fingerprint (Sparse_Fingerprint_Creator & sfc,
                          IWString_and_File_Descriptor & output)
{
  if (truncate_counted_fingerprints_to_one)
    sfc.flatten_to_01();

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);

  output << tmp << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_natoms_if_needed (int natoms,
                        IWString_and_File_Descriptor & output)
{
  if (natoms_tag.length())
    output << natoms_tag << natoms << ">\n";

  return output.good();
}

static int
write_empty_molecule_result (IWString_and_File_Descriptor & output)
{
  if (output_smiles_tag.length())
    output << output_smiles_tag << ".>\n";

  output << empty_fingerprint_with_closing_angle_bracket_and_newline;

  write_natoms_if_needed(0, output);

  return output.good();
}

/*
  Same code from map.cc, although I don't think there is any 
  reason it needs to be the same
*/

static int
form_atom_pair (int atype1,
                int distance,
                int atype2)
{
  if (atype1 > atype2)
  {
    int tmp = atype1;
    atype1 = atype2;
    atype2 = tmp;
  }

  return 975535 * atype1 + distance * 7927 + atype2;
}

static int
do_produce_atom_pair_fingerprint (Molecule & m,
                                  const int * include_these_atoms,
                                  IWString_and_File_Descriptor & output)
{
  int matoms = m.natoms();

  int * atype = new int[matoms]; std::unique_ptr<int[]> free_atype(atype);

  (void) atom_typing_specification.assign_atom_types(m, atype);

  Sparse_Fingerprint_Creator sfp;

  for (int i = 0; i < matoms; i++)
  {
    if (! include_these_atoms[i])
      continue;

    for (int j = i + 1; j < matoms; j++)
    {
      if (! include_these_atoms[j])
        continue;

      int d = m.bonds_between(i, j);

      int b = form_atom_pair(atype[i], atype[j], d);

      sfp.hit_bit(b);
    }
  }

  return write_sparse_fingerprint(sfp, output);
}

static int
fingerprint_substructure2 (Molecule & m,
                           const int * atype,
                           IWString_and_File_Descriptor & output)
{
//cerr << "SM<ILES is " << m.smiles() << " contains " << m.natoms() << " atoms " << output_smiles_tag << m.smiles() << ">\n";
  if (output_smiles_tag.length())
    output << output_smiles_tag << m.smiles() << ">\n";

  IWMFingerprint fp;

  if (nullptr != atype)
    fp.construct_fingerprint(m, atype, nullptr);
  else
    fp.construct_fingerprint(m);

  IWString tmp;

  fp.daylight_ascii_representation_including_nset_info(tmp);

  output << fingerprint_tag << tmp << ">\n";

  write_natoms_if_needed(m.natoms(), output);

  return output.good();
}

/*
  We mark any atoms we extend with -1
*/

static int
do_extend_shell (Molecule & m,
                 int * include_these_atoms,
                 int depth)
{
  int matoms = m.natoms();

  int atoms_added = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (include_these_atoms[i] > 0)
      ;
    else if (-99 == include_these_atoms[i])
      continue;
    else if (0 == include_these_atoms[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);
      if (0 != include_these_atoms[k])
        continue;

      include_these_atoms[k] = -99;

      atoms_added++;
    }
  }

  if (0 == atoms_added)
    return 1;

  for (int i = 0; i < matoms; i++)
  {
    if (-99 == include_these_atoms[i])
      include_these_atoms[i] = -(depth + 1);
  }

  if (0 == depth)
    return 1;

  return do_extend_shell(m, include_these_atoms, depth - 1);
}

static int
join_disconnected_fragments_2(Molecule & m,
                            atom_number_t a1,
                            atom_number_t astop,
                            int * include_these_atoms,
                            int flag)
{
  int d = m.bonds_between(a1, astop);

#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
  cerr << "join_disconnected_fragments, atom " << a1 << " to " << astop << " " << d << " bonds\n";
#endif

  const Atom * a = m.atomi(a1);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a1, i);

    if (j == astop)   // got to the end
      return 1;

    if (d - 1 != m.bonds_between(j, astop))
      continue;

    if (0 == include_these_atoms[j])
      include_these_atoms[j] = flag;

    join_disconnected_fragments_2(m, j, astop, include_these_atoms, flag);
  }

  return 1;
}

static int
join_disconnected_fragments(Molecule & m,
                            const Set_of_Atoms & a1,
                            const Set_of_Atoms & a2,
                            int * include_these_atoms,
                            int flag)
{
  int n = a1.number_elements();

  assert (n > 0);

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = a1[i];

    for (int k = 0; k < n; k++)
    {
      atom_number_t l = a2[k];

      join_disconnected_fragments_2(m, j, l, include_these_atoms, flag);
    }
  }

  return 1;
}

static int 
grow_group(const Molecule & m,
           atom_number_t zatom,
           int * include_these_atoms, 
           int flag)
{
  include_these_atoms[zatom] = flag;

  int rc = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (1 != include_these_atoms[j])
      continue;

    rc += grow_group(m, j, include_these_atoms, flag);
  }

  return rc;
}

/*
  We have possibly disconnected sections. Join them if they are
  close enough. We use the flag -2 to denote atoms flagged this way
*/

static int
do_join_disconected_sections (Molecule & m,
                              int * include_these_atoms)
{
  int f = 2;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (include_these_atoms[i] <= 0)
      continue;

    if (include_these_atoms[i] > 1)
      continue;

    grow_group(m, i, include_these_atoms, f);
    f++;
  }

//cerr << "Found " << (f - 2) << " disconnected groups\n";
  disconnected_groups[f - 2]++;

  if (3 == f)    // just one grouping
    return 1;

// We have possibly separated groups.  Join them by finding shortest
// paths between the disconnected sections

//iw_write_array (include_these_atoms, matoms, "include_these_atoms", cerr);

  for (int f1 = 2; f1 < f; f1++)
  {
    for (int f2 = f1 + 1; f2 < f; f2++)
    {
      int closest_separation = matoms + 2;
      Set_of_Atoms a1;
      Set_of_Atoms a2;

      for (int i = 0; i < matoms; i++)
      {
        if (f1 != include_these_atoms[i])
        continue;

        for (int j = i + 1; j < matoms; j++)
        {
          if (f2 != include_these_atoms[j])
            continue;

          int dij = m.bonds_between(i, j);
#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
          cerr << dij << " bonds between " << i << " and " << j << endl;
#endif

          if (dij < closest_separation)
          {
            a1.resize_keep_storage(0);
            a2.resize_keep_storage(0);
            a1.add(i);
            a2.add(j);
            closest_separation = dij;
          }
          else if (dij == closest_separation)
          {
            a1.add(i);
            a2.add(j);
          }
        }
      }

      assert(closest_separation > 0);

#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
      cerr << f1 << " to " << f2 << ", closest separation " << closest_separation << endl;
#endif

      if (closest_separation <= 1)    // should never be 0
        continue;

      if (closest_separation > join_disconected_sections)
        continue;
      
      assert (a1.number_elements() == a2.number_elements());
      assert (a1.number_elements() > 0);

      join_disconnected_fragments(m, a1, a2, include_these_atoms, -2);
    }
  }

  return 1;
}

static int
do_fingerprint_each_substructure_match (Molecule & m,
                                        IWString_and_File_Descriptor & output)
{
  const int matoms = m.natoms();

  m.recompute_distance_matrix();

  int * tmp = new int[matoms + matoms]; std::unique_ptr<int[]> free_tmp(tmp);
  int * atype = tmp + matoms;

  atom_typing_specification.assign_atom_types(m, atype);

  Sparse_Fingerprint_Creator sfc;

  Molecule_to_Match target(&m);

  const int nq = queries.number_elements();

  int queries_matching = 0;

  for (int i = 0; i < nq; ++i)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    if (verbose > 2)
      cerr << m.name() << " " << nhits << " hits to query " << i << endl;

    queries_matching++;

//  mark all atoms within extend_shell of the matched atoms with the number 2. Matched atoms get 1

    set_vector(tmp, matoms, 0);

    for (int j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms & e = (*sresults.embedding(j));

      e.set_vector(tmp, 1);    // note value of 1

      const int n = e.number_elements();

      for (int k = 0; k < n; ++k)
      {
        const atom_number_t l = e[k];     // atom l is a matched atom

        for (int a = 0; a < matoms; ++a)
        {
          if (0 != tmp[a])     // already marked
            continue;

          if (m.bonds_between(l, a) <= extend_shell)
            tmp[a] = 2;     // note value 2
        }
      }
    }

    circular_fingerprint_generator.generate_fingerprint(m, atype, tmp, sfc);
  }

  if (0 == queries_matching)    // should do more about this, but it really does not work...
  {
    molecules_not_matching_queries++;
    if (discard_molecules_not_matching_any_queries)
      return 1;
  }

  return write_sparse_fingerprint(sfc, output);
}

static int
fingerprint_substructure (Molecule & m,
                          int * include_these_atoms,
                          int embeddings_processed,
                          IWString_and_File_Descriptor & output)
{
  if (extend_shell)
    do_extend_shell(m, include_these_atoms, extend_shell - 1);

  if (join_disconected_sections && embeddings_processed > 1)
    do_join_disconected_sections(m, include_these_atoms);

  const int matoms = m.natoms();

  int * isotopes = new_int(matoms); std::unique_ptr<int[]> free_isotopes(isotopes);

  if (join_disconected_sections || extend_shell)
  {
    for (int i = 0; i < matoms; i++)
    {
      if (include_these_atoms[i] < 0)
      {
        isotopes[i] = - include_these_atoms[i];
        include_these_atoms[i] = 1;
      }
      else if (include_these_atoms[i] > 0)
      {
        isotopes[i] = include_these_atoms[i];
        include_these_atoms[i] = 1;
      }
    }
  }

  if (stream_for_labelled_molecules.active())
  {
    std::unique_ptr<isotope_t[]> isosave = m.GetIsotopes();

    m.transform_to_non_isotopic_form();
    m.set_isotopes(include_these_atoms);
    stream_for_labelled_molecules.write(m);

    m.set_isotopes(isosave.get());
  }

  if (produce_atom_pair_fingerprint)
    return do_produce_atom_pair_fingerprint(m, include_these_atoms, output);

  int * xref = new int[matoms]; std::unique_ptr<int[]> free_xref(xref);

  int * atype = nullptr;

  if (atom_typing_specification.active())
  {
    atype = new int[matoms];
    atom_typing_specification.assign_atom_types(m, atype);
  }

  std::unique_ptr<int[]> free_atype(atype);           // harmless if null

  Molecule subset;

  if (! m.create_subset(subset, include_these_atoms, 1, xref))
  {
    cerr << "fingerprint_substructure::cannot create subset\n";
    return 0;
  }
//cerr << "Created subset " << subset.smiles() << endl;

  for (int i = 0; i < matoms; i++)
  {
    const atom_number_t j = xref[i];

    if (j < 0)
      continue;

    if (! m.is_aromatic(i))
      continue;
      
    subset.set_permanent_aromatic(j, 1);
  }

// Now bonds

  for (int i = 0; i < matoms; i++)
  {
    const atom_number_t j = xref[i];

    if (j < 0)
      continue;

    if (! m.is_aromatic(i))
      continue;

    const Atom * a = m.atomi(i);

    const int acon = a->ncon();

    for (int k = 0; k < acon; k++)
    {
      const atom_number_t l = a->other(i, k);

      const atom_number_t n = xref[l];

      if (n < 0)
        continue;

      if (! m.in_same_aromatic_ring(i, l))
        continue;

      Bond * b = const_cast<Bond *>(subset.bond_between_atoms(j, n));

      b->set_permanent_aromatic(1);       // just one of setting the bond to perm arom, or setting to single should be enough

      subset.set_bond_type_between_atoms(j, n, SINGLE_BOND);
    }
  }

  if (nullptr != atype)
  {
    for (int i = 0; i < matoms; ++i)
    {
      if (xref[i] < 0)
        continue;

      assert(xref[i] <= i);    // this only works because we have removed atoms and the new atom number should be less than where we started

      atype[xref[i]] = atype[i];
    }
  }

  if (tag_for_isotopically_labelled_parent.length())
  {
    m.set_isotopes(isotopes);
    output << tag_for_isotopically_labelled_parent << m.smiles() << ">\n";
  }
//cerr << "Processing " << subset.smiles() << endl;

  return fingerprint_substructure2(subset, atype, output);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  return;
}

static int
fingerprint_substructure (Molecule & m,
                          IWString_and_File_Descriptor & output)
{
  output.write_if_buffer_holds_more_than(8192);

  molecules_processed++;

  m.revert_all_directional_bonds_to_non_directional();   // until I get that working properly

  if (fingerprint_each_substructure_match)
    return do_fingerprint_each_substructure_match(m, output);

  int nq = queries.number_elements();

  Molecule_to_Match target(&m);

  int * include_these_atoms = new_int(m.natoms()); std::unique_ptr<int[]> free_include_these_atoms(include_these_atoms);

  int queries_matching = 0;

  int embeddings_processed = 0;

  int inc = 1;
  if (extend_shell && tag_for_isotopically_labelled_parent.length())
    inc = extend_shell + 1;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    if (verbose > 2)
      cerr << m.name() << " " << nhits << " hits to query " << i << endl;

    queries_matching++;

    int jstop = nhits;
    if (break_at_first_match)
      jstop = 1;

    for (int j = 0; j < jstop; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);
      e->set_vector(include_these_atoms, inc);
      // cerr << "Doing embedding " << (*e) << endl;

      embeddings_processed++;
    }

    if (break_at_first_match)
      return fingerprint_substructure(m, include_these_atoms, embeddings_processed, output);
  }

  if (0 == queries_matching)
  {
    molecules_not_matching_queries++;
    if (ignore_molecules_not_matching_any_queries)
      return write_empty_molecule_result(output);

    cerr << "None of " << nq << " queries matched '" << m.name() << "'\n";
    return 0;
  }

  return fingerprint_substructure(m, include_these_atoms, embeddings_processed, output);
}

static int
fingerprint_substructure_filter (const const_IWSubstring & smiles,
                                 IWString_and_File_Descriptor & output)
{
  Molecule m;

  if (! m.build_from_smiles(smiles))
  {
    cerr << "fingerprint_substructure_filter:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return fingerprint_substructure(m, output);
}

static int
fingerprint_substructure_filter (iwstring_data_source & input,
                                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with(input_smiles_tag))
      continue;

    buffer.remove_leading_chars(input_smiles_tag.length());
    buffer.chop();

    if (! fingerprint_substructure_filter(buffer, output))
    {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      return 0;
    }

  }

  return output.good();
}

static int
fingerprint_substructure_filter (IWString_and_File_Descriptor & output)
{
  iwstring_data_source input("-");

  if (! input.good())
  {
    cerr << "fingerprint_substructure_filter:cannot open stdin?\n";
    return 0;
  }

  return fingerprint_substructure_filter(input, output);
}

static int
fingerprint_substructure (data_source_and_type<Molecule> & input,
                          IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    output << input_smiles_tag << m->smiles() << ">\n";
    output << identifier_tag << m->name() << ">\n";

    preprocess(*m);

    if (! fingerprint_substructure(*m, output))
      return 0;

    output << "|\n";
  }

  return output.good();
}

static int
fingerprint_substructure (const char * fname, FileType input_type, 
                          IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

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

  return fingerprint_substructure(input, output);
}

static int
fingerprint_substructure (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lq:s:Y:fx:S:z:T:N:i:I:o:j:P:MJ:e");

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

  if (cl.option_present('N'))
  {
    cl.value('N', natoms_tag);
    if (verbose)
      cerr << "Will write the number of atoms in the subset to the '" << natoms_tag << "' dataitem\n";

    if (! natoms_tag.ends_with('<'))
      natoms_tag.add('<');
  }

  if (cl.option_present('S'))
  {
    int i = 0;
    IWString s;
    while (cl.value('S', s, i++))
    {
      s.to_uppercase();

      if (s.starts_with("IN="))
      {
        s.remove_leading_chars(3);
        input_smiles_tag = s;
        if (! input_smiles_tag.ends_with('<'))
          input_smiles_tag.add('<');
      }
      else if (s.starts_with("OUT="))
      {
        s.remove_leading_chars(4);
        output_smiles_tag = s;
        if (0 == output_smiles_tag.length())
          ;
        else if (! output_smiles_tag.ends_with('<'))
          output_smiles_tag.add('<');
      }
      else if (s.starts_with("ISO="))
      {
        s.remove_leading_chars(4);
        tag_for_isotopically_labelled_parent = s;
        if (! tag_for_isotopically_labelled_parent.ends_with('<'))
          tag_for_isotopically_labelled_parent << '<';
      }
      else
      {
        cerr << "Unrecognised -S qualifier '" << s << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('T')) {
    if (! element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Cannot initialise element transformations (-T)\n";
      return 0;
    }
  }

  if (cl.option_present('e'))
  {
    truncate_counted_fingerprints_to_one = 1;

    if (verbose)
      cerr << "Will truncate counted fingerprints to 0/1\n";
  }

  if (cl.option_present('T')) {
    cl.value('T', fingerprint_tag);

    if (verbose)
      cerr << "Fingerprints produced with tag '" << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with('<'))
      fingerprint_tag.add('<');
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', extend_shell) || extend_shell < 1)
    {
      cerr << "The extend shell option (-x) must be a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will extend the shell " << extend_shell << " atoms\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', fingerprint_tag);

    if (! fingerprint_tag.ends_with('<'))
      fingerprint_tag << '<';
  }

  if (cl.option_present('Y'))
  {
    set_include_bits_for_rings(3);   // set here so it can be changed on the command line

    if (! parse_misc_fingerprint_options(cl, 'Y', verbose))
    {
      cerr << "Cannot parse fingerprint options (-Y)\n";
      return 4;
    }

    if (0 == fingerprint_tag.length())
      fingerprint_tag = "FPSUB<";
  }
  else if (cl.option_present('M'))
  {
    produce_atom_pair_fingerprint = 1;

    if (verbose)
      cerr << "Will produce atom pair fingerprints\n";

    if (0 == fingerprint_tag.length())
      fingerprint_tag = "NCAPSUB<";

    output_smiles_tag = "";   // the subset molecule is never formed

    if (! cl.option_present('P'))
    {
      atom_typing_specification.set_user_specified_type(IWATTYPE_USP_Y);
    }
  }
  else
  {
    set_iwmfingerprint_nbits(1024);

    if (0 == fingerprint_tag.length())
      fingerprint_tag = "FPSUB<";
  }

  if (produce_atom_pair_fingerprint)
  {
    empty_fingerprint_with_closing_angle_bracket_and_newline << fingerprint_tag << ">\n";
  }
  else
  {
    int produce_fingerprints = iwmfingerprint_nbits();
    assert (produce_fingerprints > 0);

    IW_Bits_Base fp(produce_fingerprints);
    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    empty_fingerprint_with_closing_angle_bracket_and_newline  << fingerprint_tag << tmp << ">\n";
  }

  if (cl.option_present('P'))
  {
    const const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "Cannot initialise atom typing specification '" << p << "'\n";
      return 2;
    }
  }

  if (! cl.option_present('q') && ! cl.option_present('s'))
  {
    cerr << "Must specify a substructure query via the -q or -s option\n";
    usage(3);
  }

  queries.resize (cl.option_count('q') + cl.option_count('s') + 100);

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose))
    {
      cerr << prog_name << ": cannot process queries from -q option(s)\n";
      return 6;
    }
  }

  if (cl.option_present('s'))
  {
    const_IWSubstring smarts;
    int i = 0;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts(smarts))
      {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 62;
      }

      queries.add(q);
    }
  }

  if (0 == queries.number_elements())
  {
    cerr << "No queries read, use the -q and/or -s options to specify queries\n";
    usage(4);
  }

  for (int i = 0; i < queries.number_elements(); i++)
  {
    queries[i]->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++))
    {
      if ("i" == z)
      {
        ignore_molecules_not_matching_any_queries = 1;
        if (verbose)
          cerr << "Will ignore molecules not hitting any queries\n";
      }
      else if ('d' == z)
      {
        discard_molecules_not_matching_any_queries = 1;
        if (verbose)
          cerr << "Will discard molecules not hitting any queries\n";
      }
      else if ('f' == z)
      {
        break_at_first_match = 1;

        if (verbose)
          cerr << "Will only process the first of multiple matches\n";
      }
      else if ('e' == z)
      {
        fingerprint_each_substructure_match = 1;

        if (verbose)
          cerr << "Will fingerprint all substructure matches\n";
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(5);
      }
    }

    if (fingerprint_each_substructure_match &&(0 == fingerprint_tag.length() || fingerprint_tag.starts_with("FP")))
      fingerprint_tag = "NCSUB<";

    if (fingerprint_each_substructure_match && ! atom_typing_specification.active())
      atom_typing_specification.set_user_specified_type(IWATTYPE_USP_Y);
  }

  if (cl.option_present('j'))
  {
    if (! cl.value('j', join_disconected_sections) || join_disconected_sections < 1)
    {
      cerr << "The join disconnected sections option (-j) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will join disconnected sections " << join_disconected_sections << " or closer\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring i = cl.string_value('I');

    if (! cl.option_present('o'))
    {
      stream_for_labelled_molecules.add_output_type(FILE_TYPE_SMI);
    }
    else if (! stream_for_labelled_molecules.determine_output_types(cl, 'o'))
    {
      cerr << "Cannot determine output type(s)\n";
      return 4;
    }

    if (stream_for_labelled_molecules.would_overwrite_input_files(cl, i))
    {
      cerr << "-I option '" << i << "' cannot overwrite input file(s)\n";
      return 4;
    }

    if (! stream_for_labelled_molecules.new_stem(i))
    {
      cerr << "Cannot initialise stream for labelled molecules '" << i << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Labelled molecules written to '" << i << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f'))
  {
    tag_for_isotopically_labelled_parent.resize(0);

    rc = fingerprint_substructure_filter(output);
  }
  else
  {
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

    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! fingerprint_substructure(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "processed " << molecules_processed << " molecules\n";
    if (molecules_not_matching_queries)
      cerr << molecules_not_matching_queries << " did not match any queries\n";

    if (join_disconected_sections)
    {
      for (int i = 0; i < disconnected_groups.number_elements(); i++)
      {
        if (disconnected_groups[i])
          cerr << disconnected_groups[i] << " molecules had " << i << " disconnected groups\n";
      }
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fingerprint_substructure(argc, argv);

  return rc;
}
