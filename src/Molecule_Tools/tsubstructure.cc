#include <assert.h>
#include <iostream>
#include <memory>
#include <limits>
#include <time.h>

#include "iw_tdt.h"
#include "iwbits.h"
#include "cmdline.h"
#include "report_progress.h"

#include "molecule.h"
#include "path.h"
#include "qry_wstats.h"
#include "molecule_to_query.h"
#include "target.h"
#include "misc.h"
#include "istream_and_type.h"
#include "aromatic.h"
#include "output.h"
#include "iwstandard.h"
#include "charge_assigner.h"
#include "donor_acceptor.h"
#include "etrans.h"
#include "smiles.h"
#include "numass.h"
#include "rmele.h"
#include "atom_typing.h"

#include "tsubstructure_fp.h"

typedef unsigned int atype_t;

static const char * prog_name;

namespace TSubstructure {
static Number_Assigner matched_structures_number_assigner;
static Number_Assigner non_matched_structures_number_assigner;

static Chemical_Standardisation chemical_standardisation;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static Charge_Assigner charge_assigner;

static Element_Transformations element_transformations;

static Elements_to_Remove elements_to_remove;

static Atom_Typing_Specification atom_typing_specification;
static int useUserAtomTypes = 0;


/*
  If we write the results as a fingerprint, we need a dataitem
*/

static IWString fingerprint_tag("FP");

/*
  Apr 2003. I want to be able to put the number of hits into a TDT dataitem
*/

static IWString tag_for_nhits;

/*
  Jan 2004. I want to write each matching query to its own tag
*/

static IWString tag_for_hits;

/*
  Feb 2004. Stop everything once you have a given number of molecules
  that have matched any query
*/

static int stop_processing_after_this_many_molecules_matching = 0;

/*
  Jun 2009. Want to be able to easily do a self search of a file.
*/

static int perform_search_even_if_names_the_same = 1;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage: " << prog_name << " <options> <input_file1> <input_file2>...\n";
//cerr << "  -m REPOrt      report (to the screen) matches\n";
  cerr << "  -m SEPArate    write matches to each query to a separate stream\n";
  cerr << "  -m QDT         append query match details to molecule name\n";
//cerr << "  -m NONM        write non-matching molecules to the matched atom stream\n";
//cerr << "  -m NONMX       write non-matching molecules to the matched atom stream, but\n";
  cerr << "                 do not append non-match query details\n";
  cerr << "  -m alist       write matched atoms as MDL V30 atom lists in the output file\n";
  cerr << "  -m blist       write matched atoms as MDL V30 bond lists in the output file\n";
//cerr << "  -m <number>    assign sequential numbers R(%d) to matches written\n";
  cerr << "  -m <file>      write matched structures to <file>\n";
  cerr << endl;
//cerr << "  -n REPOrt      report (to the screen) non matches\n";
  cerr << "  -n SEPArate    write matches to each query to a separate stream\n";
  cerr << "  -n QDT         append query non-match details to molecule name\n";
//cerr << "  -n <number>    assign sequential numbers R(%d) to non matches written\n";
  cerr << "  -n <file>      write non matching structures to <file>\n";
  cerr << endl;
  cerr << "  -q <query>     specify query file\n";
  cerr << "  -q F:file      specify file of queries\n";
  cerr << "  -q S:file      specify file of smarts\n";
  cerr << "  -q M:file      specify file of molecules\n";
  cerr << "  -P ...        atom typing specification - determine changing atoms and searching match conditions (default=UST:AZUCORS)\n";
//cerr << "  -Q <file>      specify a file of queries (same as -q F:<file>)\n";
//cerr << "  -h             queries are in same directory as -Q <file>\n";
  cerr << "  -s <smarts>    specify smarts for search\n";
//cerr << "  -S <smarts>    specify smarts for search (use -s)\n";
  cerr << endl;
  cerr << "  -j <...>       options for labelling matched atoms. Enter '-j help' for details\n";
  cerr << "  -J <tag>       output Daylight fingerprints with tag <TAG>\n";
  cerr << "  -y <nbits>     specify the number of bits in the fingerprint file\n";
  cerr << "  -a             output hits as array\n";
  cerr << "  -Y <stem>      use <stem> as descriptor name stem\n";
  cerr << "  -c             remove chirality before doing any matching\n";
  cerr << "  -l             reduce to largest fragment before doing search\n";
  cerr << "  -f             only find one embedding of the query (default is to find all)\n";
  cerr << "  -u             find unique matches only\n";
  cerr << "  -k             don't perceive symmetrically equivalent matches\n";
  cerr << "  -r             find one embedding per root atom match\n";
//cerr << "  -p             respect aromaticity of smarts queries\n";
//cerr << "  -K             aromatic bonds lose their Kekule forms - no longer match single or double\n";
  cerr << "  -M <...>       miscellaneous query conditions, enter '-M help' for details\n";
  cerr << "  -R <fname>     create report file, how many times each query matches <fname>\n";
  cerr << "  -b             for each molecule, break after finding a query which matches\n";
  cerr << "  -B             for each molecule, break after finding a query which DOESN't match\n";
  cerr << "  -G <fname>     file for lists of matched atoms\n";
  cerr << "  -i <type>      specify input type, enter '-i help' for info\n";
  cerr << "  -o <type>      output type for all molecules written\n";
  display_standard_aromaticity_options(cerr);
  display_standard_chemical_standardisation_options(cerr, 'g');
//display_standard_charge_assigner_options(cerr, 'N');
//cerr << "  -H <...>       donor acceptor options. Enter '-H help' for details\n";
  cerr << "  -X <element>   remove elements of type <element> before doing matches\n";
  cerr << "  -T ...         element transformations, enter '-T help' for info\n";
  cerr << "  -E <symbol>    create element with symbol <symbol>\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static void
display_dash_M_options()
{
  cerr << "  -M impring     matches must include an implicit ring\n";
  cerr << "  -M noimpring   no implicit rings allowed between matched atoms\n";
  cerr << "  -M nrnr        non ring atoms only match non ring targets, M: only\n";
  cerr << "  -M print       print all embeddings found\n";
  cerr << "  -M debug       debug print every query read in (yipes!)\n";
  cerr << "  -M maxe=nn     specify max embeddings for each query to find\n";
  cerr << "  -M nosave      don't save matched atoms. Repeat for don't save embeddings\n";
  cerr << "  -M kekule      aromatic bonds will still match their Kekule form\n";
  cerr << "  -M nokekule    aromatic bonds will no longer match their Kekule form\n";
  cerr << "  -M fkekule     aromatic bonds in fully aromatic rings no longer match Kekule forms\n";
  cerr << "  -M organic     discard hits that are in a non-organic fragment\n";
  cerr << "  -M TDT         input is a TDT, will add fingerprint. Must use -J to specify tag\n";
  cerr << "  -M bitrep=<n>  number of bit replicats when creating fingerprints\n";
  cerr << "  -M nhtag=<tag> output TDT, insert number of hits in the <tag> dataitem\n";
  cerr << "  -M htag=<tag>  output TDT, insert id of each query that matches\n";
  cerr << "  -M time        report timing\n";
  cerr << "  -M report=nn   report progress every <nn> molecules processed\n";
  cerr << "  -M meach       match each query of a multi-component query - ignores operators\n";
#ifdef DO_GETPSF
  cerr << "  -M getpsf      special processing for Dan Robertson's getpsf script\n";
#endif
  cerr << "  -M ncon=xxx    number of connections to matched atoms\n";
  cerr << "  -M dbe=xxx     distance between embeddings\n";
//cerr << "  -M aeuao       Atom Environments match Unmatched Atoms Only - unmatched during matching\n";
  cerr << "  -M aecmm       Atom Environments Can Match Matched Atoms- matched during matching\n";
  cerr << "  -M imp2exp     make implicit Hydrogens explicit\n";
  cerr << "  -M Hmmin       H in a smarts means minimum of 1 Hydrogen\n";
  cerr << "  -M vH          the v smarts directive includes implicit Hydrogens\n";
  cerr << "  -M stopm=nn    stop processing after <nn> molecules have matched\n";
  cerr << "  -M omlf        only make matches in the largest fragment\n";
  cerr << "  -M onlysubiso  isotopic atoms in query molecules indicate substitution points\n";
  cerr << "  -M edno        embedings do not overlap\n";
  cerr << "  -M nosm        do not record self matches - ID's the same\n";
  cerr << "  -M xrianum     turn off the respect initial atom numbering flag (obscure)\n";
#ifdef INTERNAL_DISTRIBUTION
  cerr << "  -M usefp       use fingerprints for queries created from molecules\n";
#endif // INTERNAL_DISTRIBUTION
  cerr << "  -M ucez        interpret uppercase smarts letters as element specifiers only\n";
  cerr << "  -M ecount      initialise element counts in targets - may help '-q M:...'\n";
  cerr << "  -M mmaq        must match all queries in order for a match to be perceived\n";
  cerr << "  -M nsssr       within target objects, nrings includes non-sssr rings\n";
  cerr << "  -M DCF=fname   create a csib directcolorfile <fname>\n";
  cerr << "  -M CEH         Condense Explicit Hydrogens to anchor atom(s)\n";
  cerr << "  -M owdmm       only write descriptors for molecules that match\n";
  cerr << "  -M anmatch     array output will be the number of hits (rather than column per query)\n";

  return;
}

static void
display_dash_j_options()
{
  cerr << "  -j <number>    label matched atoms as isotopes. Offset hits by <number>\n";
  cerr << "  -j <symbol>    changed matched atoms to type <symbol>\n";
  cerr << "  -j n=<number>  label only the first <number> of each set of matched atoms\n";
  cerr << "  -j same        all matched atoms get the same isotopic label\n";
  cerr << "  -j -n          existing isotopic labels on matched atoms are decremented by n\n";
  cerr << "  -j +n          existing isotopic labels on matched atoms are incremented by n\n";
  cerr << "  -j query       label by unique id of the query atoms\n";
  cerr << "  -j iquery      label by initial query atom number of the query atoms\n";
  cerr << "  -j mqoffset=<n> same as iquery, but offset each numbering by <n> * query number\n";
  cerr << "  -j qnum        label by query number\n";
  cerr << "  -j writeach=fname write individually labelled matches to <fname>\n";
  cerr << "                 each match written as a separate molecule\n";
  cerr << "  -j qmatch      isotopic label will be number of times any query hits an atom. No other -j possible\n";
  cerr << "  -j amap        use atom map numbers rather than isotopes\n";

  return;
}

int verbose = 0;
static int molecules_read = 0;
static int molecules_which_match = 0;
static int print_embeddings = 0;

static Report_Progress report_progress;

static int work_as_filter  = 0;

static TSubstructure_FP tsubfp;

/*
  When the report option is used, this is its stream
*/

static std::ofstream stream_for_report;

static int default_fingerprint_nbits = 0;

static int write_as_array = 0;

static int array_output_is_number_of_matches = 0;

static int only_write_array_for_molecules_with_hits = 0;

static int debug_print_queries = 0;

/*
  When processing multiple queries we have the option of stopping once
  we get a match, or matching all queries with the molecule.
  The default is to run all queries
*/

static int break_at_first_match = 0;
static int break_at_first_non_match = 0;

static int must_match_all_queries = 0;

static int reduce_to_largest_fragment = 0;

static int remove_all_chiral_centres = 0;

static int discard_hits_in_non_organic_fragments = 0;

static int hits_in_non_organic_fragments_discarded = 0;

/*
  There are two pre-defined output stream, controlled by the 'file'
  argument to the -M and -N switches.
*/

static Molecule_Output_Object stream_for_non_matching_structures;
static Molecule_Output_Object stream_for_matching_structures;
static Molecule_Output_Object stream_for_multiple_matches;

/*
  Sometimes people are interested in having the match count appended
  to the name regardless of the number of matches
*/

static int write_non_matches_to_matched_stream = 0;

/*
  If you are also using QDT, you may or may not want the non-match
  details appended to the molecules
*/

static int NONM_append_non_match_query_details = 0;

static int write_matched_atoms_as_mdl_v30_atom_lists = 0;
static int write_matched_atoms_as_mdl_v30_bond_lists = 0;

static int report_matches = 0;
static int report_non_matches = 0;
static int report_multiple_matches = -1;

/*
  When reporting multiple matches, the default is to consider two
  matches to one query the same as one match to each of two different
  queries.
*/

static int multiple_matches_each_query_counts_one = 0;

/*
  We can also specify a minimum number of hits and a maximum number of
  hits needed for a "match"
*/

static int min_hits_needed = 0;
static int molecules_with_too_few_hits = 0;
static int max_hits_needed = std::numeric_limits<int>::max();
static int molecules_with_too_many_hits = 0;

/*
  Sometimes it is interesting to label the matched atoms
*/

static int label_matched_atoms = 0;

/*
  But we may only want to label the first N of the matched atoms
*/

static int label_matched_atoms_stop = -1;

/*
  Aug 2001, Dan Robertson wanted each labelled match written separately
  to a file
*/

static Molecule_Output_Object stream_for_individually_labelled_matches;

static int report_timing = 0;

/*
  The matched atoms can be labelled as isotopes, or transformed into
  a new atom type
*/

static const Element * matched_atoms_element = NULL;

/*
  By default, the atoms of each hit are labelled with different isotopes.
  As an option, we can label them all the same
*/

static int all_matched_atoms_get_same_isotope = 0;

/*
  Nov 99. The software can apply an isotopic offset to matched atoms
*/

static int positive_j_offset = 0;
static int negative_j_offset = 0;

/*
  Aug 2017. Allow use of atom map numbers for matched atoms
*/

static int label_matched_atoms_via_atom_map_numbers = 0;

/*
  Note that the directcolorfile implementation is flawed because we don't check
  on there being the same number of entries in the dcf file as whatever file it
  is being compared with
*/

static IWString default_directcolorfile_colour("red");

static IWString_and_File_Descriptor stream_for_directcolorfile;

/*
  Jan 2000. Want a means of determining Substructure_Atom identifiers

  Set to 1 to use inique_id(), set to 2 to use initial_atom_numbers
*/

static int label_by_query_atom_number = 0;

/*
  Jul 2002. Want to label by which query has hit
*/

static int label_by_query_number = 0;

/*
  Dec 2013. For each molecule, we want to know how many of the queries hit a given atom
*/

static int increment_isotopic_labels_to_indicate_queries_matching = 0;

/*
  Feb 2003. Special thing for Dan Robertson
*/

#ifdef DO_GETPSF
static IWString getpsf_directory;
#endif

static int make_implicit_hydrogens_explicit = 0;

static int
write_mdl_v30_bond_list (Molecule & m,
                         const const_IWSubstring & query_name,
                         const Set_of_Atoms & e,
                         int * btmp,
                         std::ostream & output)
{
  int nb = m.nedges();

  set_vector(btmp, nb, 0);

  int n = e.number_elements();

  cerr << "Embedding contains " << n << " atoms\n";

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

    for (int k = i + 1; k < n; k++)
    {
      atom_number_t l = e[k];

      int ndx = m.which_bond(j, l);

      cerr << "Atoms " << j << " to " << l << " bond number " << ndx << endl;

      if (ndx >= 0)
        btmp[ndx] = 1;
    }
  }

  return m.write_set_of_bonds_as_mdl_v30_collection(btmp, query_name, "", output);
}

static int
write_matched_atoms_as_mdl_v30_lists (Molecule & m,
                                resizable_array_p<Substructure_Hit_Statistics> & queries,
                                Substructure_Results * sresults,
                                int * btmp,
                                std::ostream & output)
{
  int nq = queries.number_elements();

  m.write_molecule_mdl_v30(output, "", 0);

  for (int i = 0; i < nq; i++)
  {
    const Substructure_Results & sri = sresults[i];

    int nhits = sri.number_embeddings();

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sri.embedding(j);

      if (write_matched_atoms_as_mdl_v30_atom_lists)
        e->write_as_mdl_v30_collection_block(queries[i]->comment(), "", output);
      else
        write_mdl_v30_bond_list(m, queries[i]->comment(), *e, btmp, output);
    }
  }

  output << "M  V30 END CTAB\n";
  output << "M  END\n";
  output << "$$$$\n";

  return output.good();
}

static int
write_matched_atoms_as_mdl_v30_lists (Molecule & m,
                                resizable_array_p<Substructure_Hit_Statistics> & queries,
                                Substructure_Results * sresults,
                                std::ostream & output)
{
  int * btmp;

  if (write_matched_atoms_as_mdl_v30_bond_lists)
    btmp = new int[m.nedges()];
  else
    btmp = NULL;

  int rc = write_matched_atoms_as_mdl_v30_lists(m, queries, sresults, btmp, output);

  if (NULL != btmp)
    delete [] btmp;

  return rc;
}

/*
  Somewhat dangerous here. We ask for fragment operations on the Molecule
  underlying the Molecule_to_Match object. But since we are done searching
  we are probably OK
*/

static int
do_discard_hits_in_non_organic_fragments (Molecule_to_Match & target,
                                          Substructure_Results & sresults)
{
  const Molecule * m = target.molecule();

  int nhits = sresults.number_embeddings();

  if (m->organic_only())
    return nhits;

  resizable_array<int> non_organic_fragments;

  int matoms = target.natoms();
  for (int i = 0; i < matoms; i++)
  {
    Target_Atom & ai = target[i];

    if (! ai.element()->organic())
      non_organic_fragments.add(ai.fragment_membership());
  }

  for (int i = 0; i < sresults.number_embeddings(); i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    atom_number_t a = e->item(0);

    int f = target[a].fragment_membership();

    if (non_organic_fragments.contains(f))
    {
      sresults.remove_embedding(i);
      i--;

      hits_in_non_organic_fragments_discarded++;
    }
  }

  return sresults.number_embeddings();
}

static int
do_hits_tag_output (Molecule & m,
                    const resizable_array_p<Substructure_Hit_Statistics> & queries,
                    const int * hits,
                    std::ostream & output)
{
  if (! work_as_filter)
    write_smiles_and_pcn(m, output);

  int nq = queries.number_elements();

  for (int i = 0; i < nq; i++)
  {
    if (0 == hits[i])
      continue;

    output << tag_for_hits << queries[i]->comment() << " hits=" << hits[i] << ">\n";
  }

  if (! work_as_filter)
    output << "|\n";

  return output.good();
}

static int
do_nhits_tag_output (Molecule & m,
                     int nq,
                     const int * hits,
                     std::ostream & output)
{
  if (! work_as_filter)
    write_smiles_and_pcn(m, output);

  output << tag_for_nhits << sum_vector(hits, nq) << ">\n";

  if (! work_as_filter)
    output << "|\n";

  return output.good();
}

static int
do_array_output (const IWString & mname,
                 int nhits,
                 int nq,
                 const int * hits,
                 std::ostream & output)
{
  if (only_write_array_for_molecules_with_hits && 0 == nhits)
    return 1;

  write_space_suppressed_string(mname, output);

  if (array_output_is_number_of_matches)
  {
    int rc = 0;
    for (int i = 0; i < nq; ++i)
    {
      rc += hits[i];
    }

    output << ' ' << rc << endl;

    return 1;
  }

  for (int i = 0; i < nq; i++)
  {
    output << ' ' << hits[i];
  }

  output << endl;

  return output.good();
}

#ifdef DO_GETPSF
static int
do_getpsf (Molecule & m,
           const Substructure_Results & sresults)
{
  int matoms = m.natoms();

  IWString fname(getpsf_directory);
  fname << m.name();
  fname.truncate_at_first(' ');

  std::ofstream output(fname.null_terminated_chars(), ios::out);
  if (! output.good())
  {
    cerr << "Cannot open getpsf file '" << fname << "'\n";
    return 0;
  }

  int * already_written = new_int(matoms); std::unique_ptr<int[]> free_already_written(already_written);

  int ne = sresults.number_embeddings();

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    int n = e->number_elements();

    for (int j = 0; j < n; j++)
    {
      atom_number_t a = e->item(j);

      if (already_written[a])
        continue;

      output << m.atomic_symbol(a) << (a + 1) << endl;

      already_written[a] = 1;
    }
  }

  return output.good();
}

#endif

/*
  Apr 99. The ability to get the matched atoms out in a 
  parsable form is requested. A typical line might look like

  ID (1 2 3)(4 5 6),,(2),,,,,,,(3 2 1 4 5 6)

  There were two matched to query 0: atoms (1 2 3) and (4 5 6). No 
  matches to query 1, one match to query 2, which matched just a 
  single atom, etc.., etc...
*/

static std::ofstream bob_coner_stream;

static int
write_bob_coner_special (const Substructure_Results * sresults,
                         int nq, 
                         const IWString & mname)
{
  if (mname.contains(' '))
  {
    IWString tmp(mname);
    tmp.gsub(' ', '_');
    bob_coner_stream << tmp;
  }
  else
    bob_coner_stream << mname;

  bob_coner_stream << ' ';

  for (int i = 0; i < nq; i++)
  {
    const Substructure_Results & s = sresults[i];
    if (i > 0)
      bob_coner_stream << ',';

    int nhits = s.number_embeddings();
    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = s.embedding(j);

      bob_coner_stream << '(';
      for (int k = 0; k < e->number_elements(); k++)
      {
        if (k > 0)
          bob_coner_stream << ' ';

        bob_coner_stream << e->item(k);
      }
      bob_coner_stream << ')';
    }
  }

  bob_coner_stream << '\n';

  return bob_coner_stream.good();
}

static int
write_directcolorfile (const Molecule & m,
                       const Substructure_Results & sresults,
                       IWString_and_File_Descriptor & output)
{
  int nhits = sresults.number_embeddings();

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    int n = e->number_elements();

    for (int j = 0; j < n; j++)
    {
      atom_number_t k = e->item(j);

      output << k << ' ' << default_directcolorfile_colour << '\n';

      const Atom * ak = m.atomi(k);

      int kcon = ak->ncon();

      for (int l = 0; l < kcon; l++)
      {
        atom_number_t a2 = ak->other(k, l);

        if (a2 < k)
          continue;

        if (e->contains(a2))
          output << k << ' ' << a2 << ' ' << default_directcolorfile_colour << '\n';
      }
    }
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
do_increment_isotopic_labels_to_indicate_queries_matching (const Molecule & m,
                                const Set_of_Atoms & embedding,
                                const int na,
                                int * atom_isotopic_label)
{
  const int matoms = m.natoms();

  for (auto i = 0; i < na; ++i)
  {
    const atom_number_t j = embedding[i];

    if (j < 0 || j >= matoms)
      continue;

    atom_isotopic_label[j]++;
  }

  return 1;
}


/*
  apply the labels for a single embedding
*/

static int
apply_isotopic_labels_embedding (const Molecule & m,
                                 int query_number,
                                 const Substructure_Results & sresults,
                                 int embedding_number,
                                 int * atom_isotopic_label)
{
  const Set_of_Atoms * embedding = sresults.embedding(embedding_number);

//cerr << "Processing embedding " << (*embedding) << endl;

  int na = embedding->number_elements();    // how many of the matched atoms do we process
  if (label_matched_atoms_stop > 0 && na > label_matched_atoms_stop)
    na = label_matched_atoms_stop;

  if (increment_isotopic_labels_to_indicate_queries_matching)
    return do_increment_isotopic_labels_to_indicate_queries_matching(m, *embedding, na, atom_isotopic_label);

  const Query_Atoms_Matched * qam;
  if (label_by_query_atom_number)
    qam = sresults.query_atoms_matching(embedding_number);
  else
    qam = NULL;

  int matoms = m.natoms();

// If we are numbering by initial atom number, there may be gaps in
// the embedding, for example, someone specifies initial atom numbers
// of 0 1 and 9.  The embedding will have 10 entries, but most of them
// will be INVALID_ATOM_NUMBER

// May 2005.  Note that there is a bizzare and complex interplay
// between our isotopic labelling and the value of
// _respect_initial_atom_numbering within the query object.  The
// 'iquery' and 'query' may not work as expected with
// _respect_initial_atom_numbering is set.  Too complex to figure out
// for something used so seldom - just beware

  int ndx_for_initial_atom_number = 0;

  for (int j = 0; j < na; j++)
  {
    const atom_number_t k = embedding->item(j);
    if (k >= 0 && k < matoms)
      ;
    else if (INVALID_ATOM_NUMBER == k)
      continue;
    else
    {
      cerr << "apply_isotopic_labels_embedding:invalid atom number in embedding " << k << ", matoms " << matoms << endl;
      continue;
    }

    if (0 == atom_isotopic_label[k])    // definitely not already set
      ;
    else if (positive_j_offset || negative_j_offset)   // we cannot distinguish an isotope in the initial molecule from something added elsewhere here.
      ;
    else
      continue;   // already set by a different query

    const Atom * ak = m.atomi(k);

//  Make sure we don't generate the 0 isotope

    int iso;
    if (positive_j_offset || negative_j_offset)
    {
      if (ak->is_isotope())
        iso = ak->isotope() + positive_j_offset + negative_j_offset;
      else
        iso = label_matched_atoms + positive_j_offset + negative_j_offset;

      if (iso < 0)
        iso = 0;
    }
    else if (all_matched_atoms_get_same_isotope)
      iso = label_matched_atoms;
    else if (1 == label_by_query_atom_number)
      iso = qam->item(j)->unique_id() + 1;
    else if (label_by_query_number)
      iso = query_number + 1;
    else if (2 == label_by_query_atom_number)
    {
      iso = qam->item(ndx_for_initial_atom_number)->initial_atom_number();
      ndx_for_initial_atom_number++;
    }
    else if (label_by_query_atom_number > 2)
    {
      iso = (label_by_query_atom_number * (query_number + 1)) + label_matched_atoms * embedding_number + j + 1;
      ndx_for_initial_atom_number++;
    }
    else
      iso = label_matched_atoms * embedding_number + j + 1;

    atom_isotopic_label[k] = iso;
  }

  return 1;
}

static void
do_label_matched_atoms_via_atom_map_numbers(Molecule & m,
                                const int * atom_isotopic_label)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    m.set_atom_map_number(i, atom_isotopic_label[i]);
  }

  return;
}

static int
label_match_and_write_to_stream_for_individually_labelled_matches (Molecule & m,
                        int query_number,
                        const Substructure_Results & sresults,
                        int embedding_number,
                        int * atom_isotopic_label)
{
  set_vector(atom_isotopic_label, m.natoms(), 0);

  apply_isotopic_labels_embedding(m, query_number, sresults, embedding_number, atom_isotopic_label);

  if (label_matched_atoms_via_atom_map_numbers)
    do_label_matched_atoms_via_atom_map_numbers(m, atom_isotopic_label);
  else
    m.set_isotopes(atom_isotopic_label);

  return stream_for_individually_labelled_matches.write(m);
}

static int
label_matches_and_write_to_stream_for_individually_labelled_matches (Molecule & m,
                              int query_number,
                              const Substructure_Results & sresults,
                              int * atom_isotopic_label)
{
  for (int i = 0; i < sresults.number_embeddings(); i++)
  {
    if (i > 0)
      m.transform_to_non_isotopic_form();

    if (! label_match_and_write_to_stream_for_individually_labelled_matches(m, query_number, sresults, i, atom_isotopic_label))
      return 0;
  }

  return 1;
}

static int
label_matches_and_write_to_stream_for_individually_labelled_matches (Molecule * m,
                              int query_number,
                              const Substructure_Results & sresults)
{
  Molecule mcopy(*m);
  mcopy.set_name(m->name());

  int * tmp = new int[mcopy.natoms()]; std::unique_ptr<int[]> free_tmp(tmp);

  return label_matches_and_write_to_stream_for_individually_labelled_matches(mcopy, query_number, sresults, tmp);
}

static int
label_atoms_hit (Molecule * m,
                 int query_number,
                 const Substructure_Results & sresults,
                 int * atom_isotopic_label)
{
  int ne = sresults.number_embeddings();

  for (int i = 0; i < ne; i++)
  {
    apply_isotopic_labels_embedding(*m, query_number, sresults, i, atom_isotopic_label);
  }

  return 1;
}

static int
transform_atoms_hit (const Substructure_Results & sresults,
                     const Element ** element_labels)
{
  int ne = sresults.number_embeddings();

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * s = sresults.embedding(i);
    int na = s->number_elements();

    if (label_matched_atoms_stop > 0 && na > label_matched_atoms_stop)
      na = label_matched_atoms_stop;

    for (int j = 0; j < na; j++)
    {
      atom_number_t k = s->item(j);

      element_labels[k] = matched_atoms_element;
    }
  }

  return 1;
}

static int
do_single_query (Molecule_to_Match & target,
                 int query_number,
                 Substructure_Hit_Statistics * query,
                 Substructure_Results & sresults,
                 const Element ** element_labels,
                 int * atom_isotopic_label)
{
  int nmatches = query->substructure_search(target, sresults);

  if (0 == nmatches)
    return 0;

  if (discard_hits_in_non_organic_fragments)
    nmatches = do_discard_hits_in_non_organic_fragments(target, sresults);

  if ((verbose > 1 || report_matches) ||
      (0 == multiple_matches_each_query_counts_one && report_multiple_matches > 0 && nmatches > report_multiple_matches))
    cerr << molecules_read << ": '" << target.molecule()->molecule_name() << "' " << nmatches <<
            " matches to query '" << query->comment() << "' (total matches = " << query->molecules_which_match()
            << ")\n";

  if (0 == multiple_matches_each_query_counts_one && nmatches >= report_multiple_matches && stream_for_multiple_matches.good())
    stream_for_multiple_matches.write(target.molecule());

  if (stream_for_individually_labelled_matches.active())
    label_matches_and_write_to_stream_for_individually_labelled_matches(target.molecule(), query_number, sresults);
    
  if (label_matched_atoms || label_by_query_atom_number || label_by_query_number || increment_isotopic_labels_to_indicate_queries_matching || positive_j_offset || negative_j_offset)
    (void) label_atoms_hit(target.molecule(), query_number, sresults, atom_isotopic_label);
  else if (matched_atoms_element)
    (void) transform_atoms_hit(sresults, element_labels);

  return nmatches;
}

static int
do_all_queries (Molecule & m,
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * hits,
               const Element ** element_labels,
               int * atom_isotopic_label)
{
  int rc = 0;

  if (useUserAtomTypes)
  {
		atype_t *atype = new atype_t[m.natoms()]; std::unique_ptr<atype_t[]> free_atype(atype);

		atom_typing_specification.assign_atom_types(m, atype);

		for (int atomIndex = 0 ; atomIndex != m.natoms() ; ++atomIndex)
		{
		  m.set_userAtomType(atomIndex, atype[atomIndex]);  
		}
	}

  Molecule_to_Match target(&m);
  
	
  IWString mname;
  if (! perform_search_even_if_names_the_same)
    mname = m.name();

  if (element_transformations.active())
    element_transformations.process(target);

  int total_hits_across_all_queries = 0;

  int nq = queries.number_elements();
  for (int i = 0; i < nq; i++)
  {
    if (perform_search_even_if_names_the_same)
      ;
    else if (0 == mname.length())
      cerr << "Molecule with no name detected, cannot determine self match, searching...\n";
    else if (mname == queries[i]->comment())
      continue;

    int nhits;
    if ((nhits = do_single_query(target, i, queries[i], sresults[i], element_labels, atom_isotopic_label)))
    {
      rc++;
      hits[i] = nhits;

      if (print_embeddings)
        sresults[i].print_embeddings(cerr, verbose > 1);

      total_hits_across_all_queries += nhits;

      if (break_at_first_match)
        break;
    }
    else if (break_at_first_non_match)
      break;
  }

  if (bob_coner_stream.rdbuf()->is_open())
    write_bob_coner_special(sresults, nq, m.name());

  if (stream_for_directcolorfile.is_open())
    write_directcolorfile(m, *sresults, stream_for_directcolorfile);

#ifdef DO_GETPSF
  if (getpsf_directory.length())
    do_getpsf(m, sresults[0]);
#endif

  if (multiple_matches_each_query_counts_one && rc > 1)
  {
    if (verbose)
      cerr << rc << " of " << nq << " queries matched\n";
    stream_for_multiple_matches.write(m);
  }

  if (total_hits_across_all_queries < min_hits_needed)
  {
    molecules_with_too_few_hits++;
    return 0;
  }

  if (total_hits_across_all_queries > max_hits_needed)
  {
    molecules_with_too_many_hits++;
    return 0;
  }

  if (tsubfp.active())
    tsubfp.do_fingerprint_output(m, nq, hits, std::cout);
  else if (tag_for_nhits.length())
    do_nhits_tag_output(m, nq, hits, std::cout);
  else if (tag_for_hits.length())
    do_hits_tag_output(m, queries, hits, std::cout);
  else if (write_as_array)
    do_array_output(m.name(), rc, nq, hits, std::cout);

  if (element_labels)
  {
    int matoms = m.natoms();
    for (int i = 0; i < matoms; i++)
    {
      if (NULL != element_labels[i])
        m.set_element(i, element_labels[i]);
    }
  }

  if (atom_isotopic_label)
  {
    const int matoms = m.natoms();
    for (int i = 0; i < matoms; i++)
    {
      if (label_matched_atoms_via_atom_map_numbers)
        m.set_atom_map_number(i, atom_isotopic_label[i]);
      else
        m.set_isotope(i, atom_isotopic_label[i]);
    }
  }

  return rc;
}

static int
tsubstructure (Molecule & m, 
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched)
{
  if (report_progress())
    cerr << "Processed " << molecules_read << " molecules, " << molecules_which_match << " molecules matched\n";

  int tsize = queries.number_elements();
  if (tsize < default_fingerprint_nbits)
    tsize = default_fingerprint_nbits;

  int * tmp = new_int(tsize); std::unique_ptr<int[]> free_tmp(tmp);

  int matoms = m.natoms();

  int * atom_isotopic_label;

  if (label_matched_atoms || label_by_query_atom_number || label_by_query_number || increment_isotopic_labels_to_indicate_queries_matching || negative_j_offset || positive_j_offset)
  {
    atom_isotopic_label = new_int(matoms);
    if (negative_j_offset|| positive_j_offset || increment_isotopic_labels_to_indicate_queries_matching)
      m.get_isotopes(atom_isotopic_label);
  }
  else
    atom_isotopic_label = NULL;
    
 
  const Element ** new_elements = NULL;

  if (matched_atoms_element)
  {
    new_elements = new const Element *[matoms];
    for (int i = 0; i < matoms; i++)
    {
      new_elements[i] = NULL;
    }
  }

  int nmatched = do_all_queries(m, queries, sresults, tmp, new_elements, atom_isotopic_label);

  if (NULL != atom_isotopic_label)
    delete [] atom_isotopic_label;

  if (new_elements)
    delete [] new_elements;

  queries_matched[nmatched]++;

  if (verbose > 1 && queries.number_elements() > 1)
    cerr << "Matched " << nmatched << " of " << queries.number_elements() << " queries\n";

  if (must_match_all_queries && nmatched != queries.number_elements())
    nmatched = 0;

  if (nmatched)
  {
    molecules_which_match++;
    if (matched_structures_number_assigner.active())
      matched_structures_number_assigner.process(m);
    if (write_matched_atoms_as_mdl_v30_atom_lists || write_matched_atoms_as_mdl_v30_bond_lists)
      write_matched_atoms_as_mdl_v30_lists(m, queries, sresults, stream_for_matching_structures.stream_for_type(SDF));
    else if (stream_for_matching_structures.good())
      stream_for_matching_structures.write(m); 
  }
  else
  {
    if (report_non_matches)
      cerr << molecules_read << ": " << nmatched << " matches: " << m.name() << endl;
    if (non_matched_structures_number_assigner.active())
      non_matched_structures_number_assigner.process(m);
    if (stream_for_non_matching_structures.good())
      stream_for_non_matching_structures.write(m);

    if (write_non_matches_to_matched_stream && stream_for_matching_structures.good())
      stream_for_matching_structures.write(m);
  }

  return 0;
}

static void
preprocess (Molecule & m)
{
  if (elements_to_remove.number_elements())
    (void) elements_to_remove.process(m);

  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (charge_assigner.active())
    charge_assigner.process(m);

  if (donor_acceptor_assigner.active())
    donor_acceptor_assigner.process(m);

  if (make_implicit_hydrogens_explicit)
    m.make_implicit_hydrogens_explicit();

  return;
}

static int
tsubstructure (data_source_and_type<Molecule> & input,
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched)
{
  assert (input.good());

  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;
    if (verbose > 1)
      cerr << "Processing " << molecules_read << " " << m->name() << endl;

    preprocess(*m);

    if (0 == m->natoms())
    {
      if (verbose)
        cerr << "Empty molecule skipped\n";
      continue;
    }

    (void) tsubstructure(*m, queries, sresults, queries_matched);

    if (stop_processing_after_this_many_molecules_matching > 0 && molecules_which_match >= stop_processing_after_this_many_molecules_matching)
    {
      if (verbose)
        cerr << "Processing halted because " << molecules_which_match << " molecules already matched\n";

      return 1;
    }
  }

  return 1;
}

static IWString smiles_tag("$SMI<");

static int
tsubstructure_filter (const const_IWSubstring & smi,
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched,
               std::ostream & output)
{
  Molecule m;
  if (! m.build_from_smiles(smi))
  {
    cerr << "Invalid smiles '" << smi << "'\n";
    return 0;
  }

  molecules_read++;

  preprocess(m);

  if (0 == m.natoms())
  {
    cerr << "Empty molecule skipped\n";
    return 0;
  }

  (void) tsubstructure(m, queries, sresults, queries_matched);

  return 1;
}

static int
tsubstructure_filter (iwstring_data_source & input,
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched,
               std::ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with(smiles_tag))
      continue;

    assert (buffer.ends_with('>'));

    molecules_read++;

    buffer.chop();
    buffer.remove_up_to_first('<');

    if (! tsubstructure_filter(buffer, queries, sresults, queries_matched, output))
    {
      cerr << "Fatal error at line " << input.lines_read() << endl;
      return 0;
    }
  }

  return output.good();
}

static int
tsubstructure_filter (const char * input_fname,
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched,
               std::ostream & output)
{
  iwstring_data_source input(input_fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << input_fname << "'\n";
    return 0;
  }

  return tsubstructure_filter(input, queries, sresults, queries_matched, output);
}

static int
tsubstructure (const char * input_fname,
               const int input_type, 
               resizable_array_p<Substructure_Hit_Statistics> & queries,
               Substructure_Results * sresults,
               int * queries_matched)
{
  assert (NULL != input_fname);

  if (work_as_filter)
    return tsubstructure_filter(input_fname, queries, sresults, queries_matched, std::cout);

  data_source_and_type<Molecule> input(input_type, input_fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << input_fname << "'\n";
    return 1;
  }

  if (verbose)
    cerr << "Processing '" << input_fname << "'\n";

  return tsubstructure(input, queries, sresults, queries_matched);
}

//template class resizable_array_p<Substructure_Hit_Statistics>;
//template class resizable_array_base<Substructure_Hit_Statistics *>;

static int
open_labelled_multiple_matches_file (Command_Line & cl,
                                     const const_IWSubstring & fname)
{
  if (cl.option_present('o'))
  {
    if (! stream_for_individually_labelled_matches.determine_output_types(cl, 'o'))
    {
      cerr << "Cannot discern output types for individually labelled matches\n";
      return 0;
    }
  }
  else
    stream_for_individually_labelled_matches.add_output_type(SMI);

  if (stream_for_individually_labelled_matches.would_overwrite_input_files(cl, fname))
  {
    cerr << "Cannot overwrite input file(s)\n";
    return 0;
  }

  if (! stream_for_individually_labelled_matches.new_stem(fname))
  {
    cerr << "Cannot open stream for indivually labelled matches '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Indivually labelled matches written to '" << fname << "'\n";

  return 1;
}

static int
open_match_or_non_match_stream (Command_Line & cl,
                                const char * match_or_non_match,
                                IWString & fname,
                                Molecule_Output_Object & stream_for)
{
  stream_for.set_verbose(verbose);

  if (! cl.option_present('o'))
  {
    if (write_matched_atoms_as_mdl_v30_atom_lists || write_matched_atoms_as_mdl_v30_bond_lists)
      stream_for.add_output_type(SDF);
    else
      stream_for.add_output_type(SMI);
  }
  else if (! stream_for.determine_output_types(cl))
  {
    cerr << "Cannot determine output types for " << match_or_non_match << endl;
    return 0;
  }

  int otype = stream_for.first_output_type();

  const char * s = suffix_for_file_type(otype);

  if (fname.ends_with(s))
    fname.chop(::strlen(s));

  if (verbose)
    cerr << "Setting stem for " << match_or_non_match << " to " << fname << endl;

  if (stream_for.would_overwrite_input_files(cl, fname))
  {
    cerr << "Cannot overwrite input file(s): '" << fname << "'\n";
    return 0;
  }

  if (! stream_for.new_stem(fname, 1))
  {
    cerr << "Cannot set new stem to '" << fname << "'\n";
    return 0;
  }

  return 1;
}

/*
  There are some common tasks associated with opening the various streams
  for molecule output
*/

static int
handle_file_opening (Command_Line & cl, 
                     IWString & fname,
                     resizable_array_p<Substructure_Hit_Statistics> & queries,
                     Molecule_Output_Object & stream_for,
                     const char * match_or_non_match,
                     int separate_file_for_each_query)
{
  if (0 == fname.length())
  {
    cerr << "No " << match_or_non_match << " structure file name specified\n";
    return 0;
  }

  if (0 == separate_file_for_each_query)
    return open_match_or_non_match_stream(cl, match_or_non_match, fname, stream_for);

// At this stage, molecules will be written to streams associated with
// the queries. These only process one output type

  int output_type;
  if (! cl.option_present('o'))
    output_type = SMI;
  else if (! process_output_type(cl, output_type))
  {
    cerr << "Must specify output type via -o option\n";
    usage(9);
  }

  if (0 == output_type)
    output_type = SMI;

  if (strstr(match_or_non_match, "non"))
  {
    for (int i = 0; i < queries.number_elements(); i++)
    {
      Substructure_Hit_Statistics * q = queries[i];
      IWString tmp = fname;
      append_digit(tmp, i);
      tmp += '.';
      append_appropriate_suffix(tmp, output_type);

      if (! q->set_stream_for_non_matches(output_type, tmp.chars()))
      {
        cerr << "Cannot set output stream for non matches for query " << i << endl;
        return 0;
      }

      if (verbose)
        cerr << "Output stream for " << match_or_non_match << " to query " << i << " is '" << tmp << "'\n";
    }
  }
  else     // matched structures
  {
    for (int i = 0; i < queries.number_elements(); i++)
    {
      Substructure_Hit_Statistics * q = queries[i];
      IWString tmp = fname;
      append_digit(tmp, i);
      tmp += '.';
      append_appropriate_suffix(tmp, output_type);

      if (! q->set_stream_for_matches(output_type, tmp.chars()))
      {
        cerr << "Cannot set output stream for non matches for query " << i << endl;
        return 0;
      }
      cerr << "Output stream for matches to query " << i << " '" <<
              q->comment() << "' is '" << tmp << "'\n";
    }
  }

  return 1;
}

/*
  The options associated with -n and -m are the same.
*/

static int
parse_set_of_options (Command_Line & cl,
                      char flag,
                      resizable_array_p<Substructure_Hit_Statistics> & queries,
                      const char * match_or_non_match,
                      int & report_matches,
                      Number_Assigner & number_assigner,
                      int & append_query_details,
                      IWString & fname,
                      Molecule_Output_Object & output_stream)
{
  if (! cl.option_present(flag))
    return 1;

  int separate_file_for_each_query = 0;

// Almost all of the options require a file name to be specified

  int need_to_open_file = 0;

  int use_vertical_bars_for_query_details = 0;

  int i = 0;
  IWString tmp;
  while (cl.value(flag, tmp, i))
  {
//  cerr << "processing -mn option '" << tmp << "' sw REPO = " << tmp.starts_with ("REPO") << endl;
    
    if (tmp.starts_with("REPO"))
    {
      report_matches = 1;
      if (verbose)
        cerr << match_or_non_match << " will be reported\n";
    }
    else if ("QDT" == tmp)
    {
      append_query_details = 1;
      need_to_open_file = 1;
      if (verbose)
        cerr << "Details of " << match_or_non_match << " will be appended to the molecule name\n";
    }
    else if ("QDTVB" == tmp)
    {
      append_query_details = 1;
      need_to_open_file = 1;
      if (verbose)
        cerr << "Details of " << match_or_non_match << " will be appended to the molecule name\n";

      for (int i = 0; i < queries.number_elements(); i++)  // bit of a kludge doing this here
      {
        queries[i]->set_use_vertical_bars_for_query_details(1);
      }
    }
    else if (tmp.starts_with("SEPA"))
    {
      separate_file_for_each_query = 1;
      need_to_open_file = 1;
      if (verbose)
        cerr << match_or_non_match << " written to separate file for each query\n";
    }
    else if ("NONM" == tmp)
    {
      if ('m' != flag)
      {
        cerr << "The 'NONM' qualifier is only recognised with the -m option\n";
        usage(7);
      }

      write_non_matches_to_matched_stream = 1;
      NONM_append_non_match_query_details = 1;

      if (verbose)
        cerr << "Molecules which don't match the query will be written to the match stream\n";
    }
    else if ("NONMX" == tmp)
    {
      if ('m' != flag)
      {
        cerr << "The 'NONMX' qualifier is only recognised with the -m option\n";
        usage(7);
      }

      write_non_matches_to_matched_stream = 1;
      NONM_append_non_match_query_details = 0;

      if (verbose)
        cerr << "Molecules which don't match the query will be written to the match stream (no QDT)\n";
    }
    else if ("alist" == tmp)
    {
      write_matched_atoms_as_mdl_v30_atom_lists = 1;

      if (verbose)
        cerr << "Matched atoms written as MDL V30 atom lists\n";
    }
    else if ("blist" == tmp)
    {
      write_matched_atoms_as_mdl_v30_bond_lists = 1;

      if (verbose)
        cerr << "Matched atoms written as MDL V30 bond lists\n";
    }
#ifdef NUMBER_ASSIGNER_NOT_USEFUL_HERE
    else if (cl.value(flag, j, i))
    {
      if (! number_assigner.initialise(j))
      {
        cerr << "Number assigner initialisation failed\n";
        return 0;
      }
      need_to_open_file = 1;
      if (verbose)
        cerr << match_or_non_match << "ed molecules assigned R(%d) numbers starting with " << j << endl;
    }
#endif
    else if (0 == fname.number_elements())    // not yet set
    {
      fname = tmp;
      if (verbose)
        cerr << match_or_non_match << " will be written to stem '" << fname << "'\n";
      need_to_open_file = 1;
    }
    else
    {
      cerr << "Unrecognised -" << flag << " option '" << tmp << "'\n";
      cerr << "File name already set to '" << fname << "'\n";
      return 0;
    }

    i++;
  }

  if (0 == need_to_open_file)
    return 1;

  if (0 == fname.number_elements())
  {
    cerr << "No file name specified for -" << flag << " option\n";
    return 0;
  }

  return handle_file_opening(cl, fname, queries, output_stream,
                     match_or_non_match, separate_file_for_each_query);
}

static int
assign_name_if_needed (Substructure_Hit_Statistics & q,
                       int query_number)
{
  if (q.comment().length() > 0)   // strip to first token
  {
    IWString zname = q.comment();
    zname.truncate_at_first(' ');
    q.set_comment(zname);
    return 0;
  }

  IWString zname;

  zname << "QRY" << query_number;

  q.set_comment(zname);

  return 1;
}

static int
tsubstructure (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "w:W:y:R:X:hj:J:E:rcls:S:t:A:K:bBfm:n:uvo:i:q:Q:ag:pkG:Y:N:H:M:x:T:P:");

  if (cl.unrecognised_options_encountered())
  {
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl, verbose))
  {
    cerr << "Cannot parse element specifications\n";
    usage(2);
  }

  int input_type = 0;

  if (cl.option_present('K'))
  {
    if (! process_standard_smiles_options(cl, verbose, 'K'))
      usage(5);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
    usage(6);

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner\n";
      return 31;
    }
  }

  if (cl.option_present('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor acceptor assigner (-H option)\n";
      usage(5);
    }
  }

  int match_first = 0;
  if (cl.option_present('f'))
  {
    match_first = 1;
    if (verbose)
      cerr << "Only the first match will be made\n";
  }

  int unique_matches_only = 0;
  if (cl.option_present('u'))
  {
    unique_matches_only = 1;
    if (verbose)
      cerr << "Only unique matches will be reported\n";
  }

  int do_not_perceive_symmetry_equivalent_matches = 0;
  if (cl.option_present('k'))
  {
    do_not_perceive_symmetry_equivalent_matches = 1;
    if (verbose)
      cerr << "Will not perceive symmetry equivalent matches\n";
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', min_hits_needed) || min_hits_needed < 0)
    {
      cerr << "The -w option must be followed by a whole number\n";
      usage(6);
    }

    if (verbose)
      cerr << "Fewer than " << min_hits_needed << " matches will be reported as non matches\n";
  }

  if (cl.option_present('W'))
  {
    if (! cl.value('W', max_hits_needed) || max_hits_needed < min_hits_needed)
    {
      cerr << "The -W option requires a positive whole number at least " << min_hits_needed << endl;
      usage(17);
    }

    if (verbose)
      cerr << "More than " << max_hits_needed << " matches will be reported as non matches\n";
  }

  if (cl.option_present('T'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 'T'))
    {
      cerr << "Cannot initialise element transformations (-T option)\n";
      usage(4);
    }
  }

  if (cl.option_present('g'))
  {
    if(! chemical_standardisation.construct_from_command_line(cl, verbose))
    {
      cerr << "Cannot parse -g option\n";
      return 61;
    }
  }

  int find_one_embedding_per_atom = 0;
  if (cl.option_present('r'))
  {
    find_one_embedding_per_atom = 1;
    if (verbose)
      cerr << "Only one embedding per root atom will be reported\n";
  }

#ifdef OLD_CODE
  if (cl.option_present('K'))
  {
    set_aromatic_bonds_lose_kekule_identity(1);
    if (verbose)
      cerr << "Aromatic bonds lose their Kekule identity\n";
  }
#endif

  if (cl.option_present('b'))
  {
    break_at_first_match = 1;
    if (verbose)
      cerr << "Will cease processing queries after first match\n";
  }

  if (cl.option_present('B'))
  {
    break_at_first_non_match = 1;
    if (verbose)
      cerr << "Will cease processing queries after first non-match\n";
  }
  
  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "INvalid atom typing specification '" << p << "'\n";
      return 1;
    }
    useUserAtomTypes = 1;
  }
  

  int report_environment_matches = 0;

// Still support the -x option

  int max_matches_to_find = 0;

  if (cl.option_present('x'))
  {
    if (! cl.value('x', max_matches_to_find) || max_matches_to_find < 1)
    {
      cerr << "The -x option requires a whole positive number\n";
      usage(23);
    }

    if (verbose)
      cerr << "A maximum of " << max_matches_to_find << " matches will be found\n";
  }

  int implicit_ring_condition = -1;

  int search_each_component = 0;

  int only_match_largest_fragment = 0;

  int embeddings_do_not_overlap = 0;

  static int respect_initial_atom_numbering = 1;

  time_t tzero = 0;

// once upon a time, we used -z to control whether the embeddings were saved
// or whether the matched atoms were saved. For historical reasons, the variable
// is still called zoption

  Min_Max_Specifier<int> ncon, dbe;

  static int max_atom_count = 0;
  static int min_atom_count = 0;

  int zoption = 0;

  static IWString echo_query_stem;

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if ("impring" == m)
      {
        implicit_ring_condition = 1;
        if (verbose)
          cerr << "An implicit ring will be required in every match\n";
      }
      else if ("noimpring" == m)
      {
        implicit_ring_condition = 0;
        if (verbose)
          cerr << "Matches must not include implicit rings\n";
      }
      else if ("nrnr" == m)
      {
        set_respect_ring_membership(1);
        if (verbose)
          cerr << "Non ring atoms in molecule queries will not match ring atoms\n";
      }
      else if ("env" == m)
      {
        report_environment_matches = 1;
        if (verbose)
          cerr << "Statistics on environment matches will be printed\n";
      }
      else if ("print" == m)
      {
        print_embeddings = 1;

        if (verbose)
          cerr << "Embeddings will be printed\n";
      }
      else if ("debug" == m)
      {
        debug_print_queries = 1;
      }
      else if (m.starts_with("maxe="))
      {
        m.remove_leading_chars(5);
        if (! m.numeric_value(max_matches_to_find) || max_matches_to_find < 1)
        {
          cerr << "The maxe= directive must be a whole positive number\n";
          usage(4);
        }

        if (verbose)
          cerr << "Each query will find a maximum of " << max_matches_to_find << " matches\n";
      }
      else if ("nosave" == m)
      {
        zoption++;
      }
      else if ("nokekule" == m)
      {
        set_aromatic_bonds_lose_kekule_identity(1);

        if (verbose)
          cerr << "Aromatic bonds will no longer match Kekule forms\n";
      }
      else if ("fkekule" == m)
      {
        set_aromatic_bonds_lose_kekule_identity(2);

        if (verbose)
          cerr << "Aromatic bonds in alternating Kekule rings will no longer match Kekule forms\n";
      }
      else if ("kekule" == m)
      {
        set_aromatic_bonds_lose_kekule_identity(0);

        if (verbose)
          cerr << "Aromatic bonds will no longer match Kekule forms\n";
      }
      else if ("time" == m)
      {
        report_timing = 1;

        if (verbose)
          cerr << "Will report timing\n";

        tzero = time(NULL);
      }
      else if ("organic" == m)
      {
        discard_hits_in_non_organic_fragments = 1;

        if (verbose)
          cerr << "Molecules containing non organics in the same fragment as the matched atoms don't match\n";
      }
      else if ("tdt" == m || "TDT" == m)
      {
        work_as_filter = 1;

        if (verbose)
          cerr << "Will add fingerprints to an existing TDT file\n";
      }
      else if (m.starts_with("bitrep="))
      {
        m.remove_leading_chars(7);

        int bit_replicates;
        if (! m.numeric_value(bit_replicates) || bit_replicates < 1)
        {
          cerr << "The bit replicates directive (bitrep=) must be a whole +ve number\n";
          usage(2);
        }

        if (verbose)
          cerr << "Will create " << bit_replicates << " bit replicates\n";

        tsubfp.set_bit_replicates(bit_replicates);
      }
      else if (m.starts_with("nhtag="))
      {
        m.remove_leading_chars(6);
        tag_for_nhits = m;
        if (verbose)
          cerr << "Number of hits inserted as '" << tag_for_nhits << "' dataitem\n";

        if (! tag_for_nhits.ends_with('<'))
          tag_for_nhits += '<';
      }
      else if (m.starts_with("htag="))
      {
        m.remove_leading_chars(5);
        tag_for_hits = m;
        if (verbose)
          cerr << "Queries matching written as tag '" << tag_for_hits << "'\n";

        if (! tag_for_hits.ends_with('<'))
          tag_for_hits.add('<');
      }
      else if ("aeuao" == m)
      {
        set_atom_environment_only_matches_unmatched_atoms(1);
        set_query_environment_must_match_unmatched_atoms(1);
      }
      else if ("aecmm" == m)
      {
        set_atom_environment_only_matches_unmatched_atoms(0);
        set_query_environment_must_match_unmatched_atoms(0);
      }
      else if ("imp2exp" == m)
      {
        make_implicit_hydrogens_explicit = 1;
      }
      else if ("Hmmin" == m)
      {
        set_h_means_exactly_one_hydrogen(0);
      }
      else if ("vH" == m)
      {
        set_global_setting_nbonds_includes_implicit_hydrogens(1);
      }
      else if (m.starts_with("report="))
      {
        m.remove_leading_chars(7);
        int r;
        if (! m.numeric_value(r) || r < 1)
        {
          cerr << "The report option must be followed by a valid whole number\n";
          usage(4);
        }

        if (verbose)
          cerr << "Will report progress every " << r << " molecules read\n";

        report_progress.set_report_every(r);
      }
      else if ("meach" == m)
      {
        search_each_component = 1;
      }
      else if (m.starts_with("echo="))
      {
        echo_query_stem = m;
        echo_query_stem.remove_leading_chars(5);
      }
#ifdef DO_GETPSF
      else if (m.starts_with("getpsf="))
      {
        getpsf_directory = m;
        getpsf_directory.remove_leading_chars(7);

        if (! getpsf_directory.ends_with('/'))
          getpsf_directory += '/';
      }
#endif
      else if (m.starts_with("ncon="))
      {
        m.remove_leading_chars(5);

        if (! ncon.initialise(m))
        {
          cerr << "Invalid ncon specification '" << m << "'\n";
          return 9;
        }
      }
      else if (m.starts_with("dbe="))
      {
        m.remove_leading_chars(4);

        if (! dbe.initialise(m))
        {
          cerr << "Invalid dbe specification '" << m << "'\n";
          return 11;
        }
      }
      else if (m.starts_with("stopm="))
      {
        m.remove_leading_chars(6);
        if (! m.numeric_value(stop_processing_after_this_many_molecules_matching) || stop_processing_after_this_many_molecules_matching < 1)
        {
          cerr << "Invalid stop match criterion '" << m << "'\n";
          return 13;
        }

        if (verbose)
          cerr << "Will stop processing after " << stop_processing_after_this_many_molecules_matching << " molecules that match are found\n";
      }
      else if ("omlf" == m)
      {
        only_match_largest_fragment = 1;
        if (verbose)
          cerr << "Will only make matches in the largest fragment of a molecule\n";
      }
      else if ("onlysubiso" == m)
      {
        set_substituents_only_at_isotopic_atoms(1);
        if (verbose)
          cerr << "When building queries from molecules, isotopic atoms indicate substitution points\n";
      }
      else if ("edno" == m)
      {
        embeddings_do_not_overlap = 1;
      }
      else if ("usefp" == m)
      {
        set_use_fingerprints_for_screening_substructure_searches(1);
      }
      else if ("ucez" == m)
      {
        set_respect_aliphatic_smarts(0);
        if (verbose)
          cerr << "Will interpret uppercase elements as atomic numbers only\n";
      }
      else if ("ecount" == m)
      {
        set_initialise_element_counts(1);
        if (verbose)
          cerr << "Will initialise element counts in targets\n";
      }
      else if ("mmaq" == m)
      {
        if (break_at_first_match || break_at_first_non_match)
        {
          cerr << "The -b and/or -B options cannot be used with 'mmaq'\n";
          usage(3);
        }

        must_match_all_queries = 1;
        if (verbose)
          cerr << "Molecules must match all queries\n";
      }
      else if ("rian" == m)
      {
      }
      else if ("nsssr" == m)
      {
        set_global_setting_nrings_includes_non_sssr_rings(1);
        if (verbose)
          cerr << "nrings includes non-sssr rings\n";
      }
      else if (m.starts_with("DCF="))
      {
        IWString dcf_fname = m;
        dcf_fname.remove_leading_chars(4);

        if (! stream_for_directcolorfile.open(dcf_fname.null_terminated_chars()))
        {
          cerr << "Cannot open directcolorfile '" << dcf_fname << "'\n";
          return 4;
        }

        if (verbose)
          cerr << "Will create directcolorfile '" << dcf_fname << "'\n";
      }
      else if (m.starts_with("DCFCOL="))
      {
        m.remove_leading_chars(7);

        default_directcolorfile_colour = m;

        if (verbose)
          cerr << "Directcolorfile colour set to '" << m << "'\n";
      }
      else if ("CEH" == m)
      {
        set_molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms (1);
        if (verbose)
          cerr << "Queries from molecules always condense explicit Hydrogens\n";
      }
      else if (m.starts_with("maxat="))
      {
        m.remove_leading_chars(6);
        if (! m.numeric_value(max_atom_count) || max_atom_count < 1)
        {
          cerr << "Invalid maxat= qualifier '" << m << "'\n";
          usage(3);
        }
      }
      else if (m.starts_with("minat="))
      {
        m.remove_leading_chars(6);
        if (! m.numeric_value(max_atom_count) || max_atom_count < 1)
        {
          cerr << "Invalid maxat= qualifier '" << m << "'\n";
          usage(3);
        }
      }
      else if ("owdmm" == m)
      {
        only_write_array_for_molecules_with_hits = 1;
        if (verbose)
          cerr << "Will only write descriptors for molecules that match queries\n";
      }
      else if ("nosm" == m)
      {
        perform_search_even_if_names_the_same = 0;
        if (verbose)
          cerr << "Performing self search, molecules do NOT check themselves\n";
      }
      else if ("xrianum" == m)
      {
        respect_initial_atom_numbering = 0;
        if (verbose)
          cerr << "Will relax strict adherence to initial atom numbering\n";
      }
      else if ("ama" == m)
      {
        set_only_aromatic_atoms_match_aromatic_atoms(1);
        if (verbose)
          cerr << "Molecule to query transformations such that only aromatic atoms match aromatic atoms\n";
      }
      else if (m.starts_with("rmm="))
      {
        m.remove_leading_chars(4);
        unsigned int r;
        if (! m.numeric_value(r) || r < 1)
        {
          cerr << "The report multiple matches value (rmm=) must be a whole +ve number '" << m << "' not recognised\n";
          return 2;
        }

        set_report_multiple_hits_threshold(r);
      }
      else if ("anmatch" == m)
      {
        array_output_is_number_of_matches = 1;
        if (verbose)
          cerr << "Array output will be number of matches\n";
      }
      else if ("help" == m)
      {
        display_dash_M_options();
        return 0;
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_M_options();
        return 0;
      }
    }

    if (work_as_filter)
    {
      if (! cl.option_present('J') && 0 == tag_for_nhits.length())
      {
        cerr << "Must specify fingerprint tag via the -J option\n";
        usage(5);
      }
    }

    if (0 == max_atom_count)
      ;
    else if (0 == min_atom_count)
      ;
    else if (min_atom_count >= max_atom_count)
    {
        cerr << "Inconsistent min (" << min_atom_count << ") and max (" << max_atom_count << ") atom count specifications\n";
        return 4;
      }
      else if (verbose)
        cerr << "Will only match atoms with between " << min_atom_count << " and " << max_atom_count << " atoms\n";
  }

  if (zoption > 1)
  {
    if (cl.option_present('j') || cl.option_present('u') || cl.option_present('k'))
    {
      cerr << "The '-M nosave' option cannot be used with the -k, -u or -j options\n";
      usage(17);
    }
  }

// the atom numbering options are complex
 
  if (cl.option_present('j'))
  {
    int i = 0;
    const_IWSubstring j;
    while (cl.value('j', j, i++))
    {
      if (j.starts_with("n="))
      {
        j.remove_leading_chars(2);
        if (! j.numeric_value(label_matched_atoms_stop) || label_matched_atoms_stop < 1)
        {
          cerr << "The -j n= specifier must be followed by a whole positive number\n";
          usage(3);
        }
        if (verbose)
          cerr << "The first " << label_matched_atoms_stop << " atoms of each match will be labelled\n";
      }
      else if ("same" == j)
      {
        all_matched_atoms_get_same_isotope = 1;
        if (verbose)
          cerr << "All atoms will be assigned the same isotopic label\n";
      }
      else if ("query" == j)
      {
        label_by_query_atom_number = 1;
        if (verbose)
          cerr << "Will label by matched substructure atom unique identifier\n";
      }
      else if ("iquery" == j)
      {
        label_by_query_atom_number = 2;
        if (verbose)
          cerr << "Will label by matched substructure atom initial atom number\n";
      }
      else if (j.starts_with("mqoffset="))
      {
        j.remove_leading_chars(9);
        if (! j.numeric_value(label_by_query_atom_number) || label_by_query_atom_number <= 2)
        {
          cerr << "Invalid query label offset '" << j << "'\n";
          return 2;
        }
        if (verbose)
          cerr << "Multiplicative offset for query number " << label_by_query_atom_number << endl;
      }
      else if (j.starts_with("qnum"))
      {
        label_by_query_number = 1;
        if (verbose)
          cerr << "Will label matched atoms by query number\n";
      }
      else if (j.starts_with('+'))
      {
        j++;     // '+' is not a valid part of a number
        if (! j.numeric_value(positive_j_offset) || positive_j_offset <= 0)
        {
          cerr << "With -j +nnn, the number, nnn must be non-zero\n";
          usage(44);
        }

        if (verbose)
          cerr << "Positive isotopic offset " << positive_j_offset << endl;
      }
      else if (j.starts_with('-'))
      {
        if (! j.numeric_value(negative_j_offset) || negative_j_offset <= 0)
        {
          cerr << "With -j -nnn, the number, nnn must be non-zero\n";
          usage(44);
        }

        if (verbose)
          cerr << "Negative isotopic offset " << negative_j_offset << endl;
      }
      else if (j.starts_with("writeach="))
      {
        j.remove_up_to_first ('=');

        if (! open_labelled_multiple_matches_file(cl, j))
        {
          cerr << "Cannot open file for labelled matches '" << j << "'\n";
          return 0;
        }
      }
      else if ("qmatch" == j)
      {
        increment_isotopic_labels_to_indicate_queries_matching = 1;

        if (verbose)
          cerr << "Will increment isotopic labels by the number of queries hitting that atom\n";
      }
      else if ("amap" == j)
      {
        label_matched_atoms_via_atom_map_numbers = 1;
        if (verbose)
          cerr << "Will label matched atoms via atom map number\n";
      }
      else if ("help" == j)
      {
        display_dash_j_options();
        return 0;
      }
      else if (j.numeric_value(label_matched_atoms))
      {
        if (label_matched_atoms < 0)
        {
          cerr << "The -j option requires a whole non-negative number\n";
          usage(37);
        }
        if (verbose)
          cerr << "Matched atoms will be isotopically labelled, multiple matches offset " << label_matched_atoms << endl;
      }
      else
      {
        matched_atoms_element = get_element_from_symbol_no_case_conversion(j);
        if (NULL == matched_atoms_element)
        {
          cerr << "Cannot retrieve element '" << j << "', -j option - should be a whole number\n";
          return 9;
        }

        if (verbose)
          cerr << "Matched atoms will be changed to type '" << matched_atoms_element->symbol() << "'\n";
      }
    }

//  make sure the various option combinations are consistent

    if (label_by_query_atom_number)
      ;
    else if (label_by_query_number)
      ;
    else if (positive_j_offset || negative_j_offset)
      ;
    else if (increment_isotopic_labels_to_indicate_queries_matching)
      ;
    else if (0 == label_matched_atoms && NULL == matched_atoms_element)
    {
      cerr << "Must specify '-j <number>' or '-j <element>' as a -j option\n";
      usage(13);
    }

    if (positive_j_offset && negative_j_offset)
    {
      cerr << "Cannot specify both a positive and negative isotopic offset\n";
      usage(22);
    }

    if ((positive_j_offset || negative_j_offset) && 0 == label_matched_atoms)
    {
      cerr << "Warning, positive (" << positive_j_offset << ") and/or negative (" << negative_j_offset << ") specified, but no label for non isotopic atoms. Any resulting negative isotopes ignored\n";
//    cerr << "In order to isotopically offset matched atoms, an initial value for non-isotopic atoms is needed\n";
//    usage(29);
    }

    if (label_by_query_atom_number && all_matched_atoms_get_same_isotope)
    {
      cerr << "Cannot label to query atom number and have all get the same isotope\n";
      usage(33);
    }

//  put more checks on label_by_query_atom_number here....
  }

  if (work_as_filter)    // no need to check input specifications
    ;
  else if (! cl.option_present('i'))
  {
    if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))   // reading a pipe, assume smiles
      input_type = SMI;
    else if (! all_files_recognised_by_suffix(cl))
    {
      cerr << "Cannot discern all file types, use the -i option\n";
      return 4;
    }
  }
  else if (! process_input_type(cl, input_type))
    usage(3);


  resizable_array_p<Substructure_Hit_Statistics> queries;

  queries.resize(256);    // hopefully large enough to avoid extra mallocs

  if (cl.option_present('q') && ! process_queries(cl, queries, verbose))
  {
    cerr << prog_name << ": cannot process queries from -q option(s)\n";
    return 6;
  }

  if (cl.option_present('Q'))
  {
    if (! process_files_of_queries(cl, queries, cl.option_present('h'), verbose))
    {
      cerr << prog_name << ": cannot process files of queries, -Q option\n";
      usage(7);
    }
  }

  if (cl.option_present('s') && cl.option_present('S'))
  {
    cerr << "Use one of -s or -S but not both\n";
    usage(32);
  }

  if (cl.option_present('S') || cl.option_present('s'))
  {
    char s;
    if (cl.option_present('S'))
      s = 'S';
    else
      s = 's';

    int i = 0;
    const_IWSubstring smarts;
    while (cl.value(s, smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts(smarts))
      {
        delete q;
        return 60 + i;
      }

      if (verbose > 1)
        cerr << "Created query from smarts '" << smarts << "'\n";

      queries.add(q);
    }
  }

  int nq = queries.number_elements();

  if (0 == nq)
  {
    cerr << prog_name << ": No queries specified, use -s, -q or -Q options\n";
    usage(12);
  }

  if (! perform_search_even_if_names_the_same)
  {
    for (int i = 0; i < nq; i++)
    {
      if (0 == queries[i]->comment().length())
      {
        cerr << "Cannot process query with no name when looking at name matches\n";
        return 3;
      }
    }
  }

//IWString fname("q0.msi");
//queries[0]->write_msi(fname);

  if (ncon.is_set())
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_ncon(ncon);
    }
  }

  if (dbe.is_set())
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_distance_between_hits(dbe);
    }
  }

  if (only_match_largest_fragment)
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_only_keep_matches_in_largest_fragment(1);
    }
  }

  if (embeddings_do_not_overlap)
  {
    for (int i = 0; i < queries.number_elements(); i++)
    {
      queries[i]->set_embeddings_do_not_overlap(1);
    }
  }

  if (min_atom_count > 0)
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_min_atoms_to_match(min_atom_count);
    }
  }

  if (max_atom_count > 0)
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_max_atoms_to_match(max_atom_count);
    }
  }

  if (0 == respect_initial_atom_numbering)
  {
    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_respect_initial_atom_numbering(0);
    }
  }

  if (echo_query_stem.length())
  {
    for (int i = 0; i < nq; i++)
    {
      IWString fname;
      fname << echo_query_stem << i << ".qry";
      queries[i]->write_msi(fname);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will perform searches on the largest fragment\n";
  }

// Process the -X option after the queries are constructed. Otherwise '-X R' will make
// R an element, so any smarts with R in it will be recognised as element R rather
// than a ring directive.

  if (cl.option_present('X'))
  {
    if (! elements_to_remove.construct_from_command_line(cl, verbose, 'X'))
    {
      usage(53);
    }
  }

  if (cl.option_present('c'))
  {
    remove_all_chiral_centres = 1;
    if (verbose)
      cerr << "Chiral centres will be removed\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', report_multiple_matches) || report_multiple_matches < 1)
    {
      cerr << "The -t option must be followed by a positive number\n";
      usage(19);
    }
    if (verbose)
      cerr << "Queries hitting more than " << report_multiple_matches << " will be reported\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(8);
  }

  if (cl.option_present('R'))
  {
    IWString fname;
    cl.value('R', fname);

    stream_for_report.open(fname.null_terminated_chars(), std::ios::out);
    if (! stream_for_report.good())
    {
      cerr << "Cannot open report file '" << fname << "'\n";
      return 63;
    }
  }

  if (cl.option_present('o'))
  {
    if (! cl.option_present('m') && ! cl.option_present('n') && ! cl.option_present('T'))
    {
      cerr << "Output type specified (-o option) but no output requested (-T, -m or -n)\n";
      usage(63);
    }
  }

  int append_match_details_to_molecule_name = 0;
  IWString fname_for_m;
  if (! parse_set_of_options(cl, 'm', queries, "matches",
                              report_matches,
                              matched_structures_number_assigner,
                              append_match_details_to_molecule_name,
                              fname_for_m,
                              stream_for_matching_structures))
  {
    cerr << "Cannot process -m options\n";
    usage(15);
  }

  int append_non_match_details_to_molecule_name = 0;
  IWString fname_for_n;
  if (! parse_set_of_options(cl, 'n', queries, "non-matches",
                              report_non_matches,
                              non_matched_structures_number_assigner,
                              append_non_match_details_to_molecule_name,
                              fname_for_n,
                              stream_for_non_matching_structures))
  {
    cerr << "Cannot process -n options\n";
    usage(16);
  }

  if (fname_for_m.length() && fname_for_n.length() && fname_for_m == fname_for_n)
  {
    cerr << "Cannot use the same file name for the -m and -n options\n";
    return 5;
  }

  if (cl.option_present('G'))
  {
    const char * g = cl.option_value('G');

    bob_coner_stream.open(g, std::ios::out);
    if (! bob_coner_stream.good())
    {
      cerr << "Cannot open matched atom stream '" << g << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Matched atom lists written to '" << g << "'\n";
  }

// Now that all queries are in, make any global changes to them.

  for (int i = 0; i < queries.number_elements(); i++)
  {
    Substructure_Hit_Statistics * qi = queries[i];

    if (max_matches_to_find)
      qi->set_max_matches_to_find(max_matches_to_find);

    if (search_each_component)
      qi->set_each_component_search(search_each_component);

//  If we are writing subsets, make sure everything has a name

    if (write_matched_atoms_as_mdl_v30_atom_lists || write_matched_atoms_as_mdl_v30_bond_lists)
      assign_name_if_needed(*qi, i);

//  cerr << "Process " << qi->number_elements() << " components of query " << i << endl;

    for (int j = 0; j < qi->number_elements(); j++)
    {
      Single_Substructure_Query * sj = qi->item(j);

      if (match_first)
      {
        sj->set_max_matches_to_find(1);
        sj->set_find_one_embedding_per_atom(1);
      }
      else if (find_one_embedding_per_atom)
        sj->set_find_one_embedding_per_atom(1);

      if (unique_matches_only)
        sj->set_find_unique_embeddings_only(1);

      if (do_not_perceive_symmetry_equivalent_matches)
        sj->set_do_not_perceive_symmetry_equivalent_matches(1);

      if (implicit_ring_condition >= 0)
        sj->set_implicit_ring_condition (implicit_ring_condition);
    }

    qi->set_append_match_details_to_molecule_name (append_match_details_to_molecule_name);
    qi->set_append_non_match_details_to_molecule_name (append_non_match_details_to_molecule_name);

//  Kludge alert. We want '-m QDT -m NONM' to work as expected

//  if (write_non_matches_to_matched_stream && append_match_details_to_molecule_name)
    if (NONM_append_non_match_query_details)
      qi->set_append_non_match_details_to_molecule_name(1);

    if (verbose > 1)
      qi->set_verbose(verbose - 1);   // change OCT 97

    if (debug_print_queries)
    {
      cerr << "\ntsubstructure: Info on query " << i << "\n\n";
      qi->debug_print(cerr);
    }

    assert (qi->ok());
  }

#ifdef OLD_DASH_T_OPTION
  if (cl.option_present('T'))
  {
    if (! report_multiple_matches)
    {
      cerr << "The -T option requires the -t option also\n";
      usage(9);
    }

    IWString fname;
    int i = 0;
    while (cl.value('T', fname, i++))
    {
      if ("QDT" == fname)
      {
        for (int j = 0; j < queries.number_elements(); j++)
        {
          queries[j]->set_append_match_details_to_molecule_name(report_multiple_matches);
        }
      }
      else if (fname.starts_with("SING"))
      {
        multiple_matches_each_query_counts_one = 1;
        if (verbose)
          cerr << "When looking for multiple matches each query counts once regardless of embeddings\n";
      }
      else
      {
        if (! handle_file_opening(cl, fname, queries, stream_for_multiple_matches,
                                   "multiple matches", 0))
        usage(7);
      }
    }
  }
#endif

  if (cl.option_present('J') && cl.option_present('a'))
  {
    cerr << "The -J (write as fingerprint) and -a (write as array) options are mutually exclusive\n";
    usage(27);
  }

  if (cl.option_present('y') && ! cl.option_present('J'))
  {
    cerr << "The -y option (number of bits to write) can only be used with the -J option (fingerprint tag)\n";
    usage(17);
  }

  if (cl.option_present('J'))
  {
    IWString fingerprint_tag;

    cl.value('J', fingerprint_tag);

    if (verbose)
      cerr << "Fingerprints produced with the '" << fingerprint_tag << " dataitem\n";

    if (! fingerprint_tag.ends_with('<'))
      fingerprint_tag += '<';

    tsubfp.set_fingerprint_tag(fingerprint_tag);

    if (cl.option_present('y'))
    {
      if (! cl.value('y', default_fingerprint_nbits) ||
            default_fingerprint_nbits < 1 ||
            default_fingerprint_nbits < queries.number_elements())
      {
        cerr << "The -y option must be followed by a whole number between 2 and " << queries.number_elements() << endl;
        usage(15);
      }
      tsubfp.set_default_fingerprint_nbits(default_fingerprint_nbits);
    }
    else
      default_fingerprint_nbits = queries.number_elements();

    if (0 != default_fingerprint_nbits % 8)
      default_fingerprint_nbits = (default_fingerprint_nbits / 8 + 1) * 8;

    if (verbose)
      cerr << "Query matches written as fingerprints to '" << cl.option_value('J') <<
              "' fingerprints size " << default_fingerprint_nbits << endl;

    tsubfp.set_default_fingerprint_nbits(default_fingerprint_nbits);
  }
  else if (cl.option_present('a') || cl.option_present('Y') || array_output_is_number_of_matches)     //  Will we write Daylight fingerprints or an array of hits
  {
    write_as_array = 1;
    if (verbose)
      cerr << "Will write hit counts as an array of int's\n";

    std::cout << "Name";      // write the descriptor header line

    IWString descriptor_stem;

    if (array_output_is_number_of_matches)
      std::cout << " Matches";
    else if (cl.option_present('Y'))
    {
      descriptor_stem = cl.string_value('Y');
      if (verbose)
        cerr << "Query descriptors written with stem '" << descriptor_stem << "'\n";

      for (int i = 0; i < queries.number_elements(); i++)
      {
        std::cout << ' ' << descriptor_stem << i;
      }
    }
    else
    {
      for (int i = 0; i < queries.number_elements(); i++)
      {
        IWString c = queries[i]->comment();
  
        if (0 == c.length())
          std::cout << " qry" << i;
        else
        {
//        c.truncate_at_first(' ');
          c.gsub(' ', '_');
          std::cout << ' ' << c;
        }
      }
    }

    std::cout << endl;
  }

// This keeps track of the number of molecules which match I queries.

  int nqueries = queries.number_elements();
  int * queries_matched = new_int(nqueries + 1); std::unique_ptr<int[]> free_queries_matched(queries_matched);

  Substructure_Results * sresults = new Substructure_Results[nqueries]; std::unique_ptr<Substructure_Results[]> free_sresults(sresults);

  if (zoption > 0)
  {
    for (int i = 0; i < nqueries; i++)
    {
      queries[i]->set_save_matched_atoms(0);

      sresults[i].set_save_matched_atoms(0);

      if (zoption > 1)
        sresults[i].set_save_query_atoms_matched(0);
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if(! tsubstructure(cl[i], input_type, queries, sresults, queries_matched))
    {
      rc = i + 1;
      break;
    }

    if (stop_processing_after_this_many_molecules_matching > 0 && molecules_which_match >= stop_processing_after_this_many_molecules_matching)
      break;
  }

  if (stream_for_matching_structures.active())
    stream_for_matching_structures.do_flush();

// report our results

  if (molecules_read > 0)
  {
    float fraction = static_cast<int>(molecules_which_match) / static_cast<float>(molecules_read);

    cerr << molecules_read << " molecules read, " << molecules_which_match << " molecules match fraction " << fraction << endl;
  }

  if (report_timing)
  {
    time_t execution_time = time(NULL) - tzero;

    cerr << "Took " << execution_time << " seconds, " << (molecules_read / execution_time) << " molecules per second\n";
  }

  if (verbose)
  {
    for (int i = 0; i < nqueries + 1; i++)
    {
      if (queries_matched[i])
        cerr << queries_matched[i] << " molecules matched " << i << " queries\n";
    }
  }

  if (verbose && molecules_with_too_few_hits)
    cerr << molecules_with_too_few_hits << " molecules had fewer than " << min_hits_needed << " hits (-w option)\n";

  if (verbose && molecules_with_too_many_hits)
    cerr << molecules_with_too_many_hits << " molecules had more than " << max_hits_needed << " hits (-W option)\n";

  if (hits_in_non_organic_fragments_discarded)
    cerr << "Discarded " << hits_in_non_organic_fragments_discarded << " hits in non-organic fragments\n";

  for (int i = 0; i < nqueries; i++)
  {
    Substructure_Hit_Statistics * q = queries[i];

    if (verbose)
    {
      cerr << i << ' ';
      q->report(cerr, report_environment_matches);
    }
    if (stream_for_report.rdbuf()->is_open())
    {
      stream_for_report << q->molecules_which_match() << ' ' << q->comment() << endl;
    }
  }

  return rc;
}
}     // end namespace

//#define DLL_FLAG=1
#ifndef DLL_FLAG

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = TSubstructure::tsubstructure(argc, argv);

  return rc;
}

#endif

#ifdef HZ_CSHARP

extern "C"
{

int tsubstructure_csharp(int argc, char **argv)
{
  return tsubstructure(argc,argv);
}

}

#endif
