/*
  Implementation of verloop descriptors
*/

#include <iomanip>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "sterimol.h"
#include "topological_sterimol.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int molecules_read = 0;

static int ignore_molecules_not_matching = 0;

static int multiple_hits_ok = 0;

static int just_take_first_of_multiple_hits = 0;

static Charge_Assigner charge_assigner;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static Elements_to_Remove elements_to_remove;

static Chemical_Standardisation chemical_standardisation;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int do_3d_computations = 0;

static int do_partial_charge_computation = 0;

#define PARTIAL_CHARGE_ABRAHAM 1
#define PARTIAL_CHARGE_GASTEIGER 2
#define PARTIAL_CHARGE_GH 3

/*
  We need a way of storing the 2d and 3d results
*/

static int nbonds = 0;

static Molecule_Output_Object stream_for_labelled_molecules;

static Topological_Sterimol* results_2d = nullptr;
static Sterimol* results_3d = nullptr;

/*
  We may also want to do a full sterimol computation on the whole molecule
*/

static int do_whole_molecule_sterimol = 0;

static int make_implicit_hydrogens_explicit = 0;

/*
  We need any number of matched atom pairs associated with each query
*/

class Matched_Atom_Pair
{
 private:
  int _a1, _a2;

  Matched_Atom_Pair* _next;

  //  private functions

  int _add(const const_IWSubstring&);

 public:
  Matched_Atom_Pair();
  ~Matched_Atom_Pair();

  int
  active() const
  {
    return INVALID_ATOM_NUMBER != _a1 && INVALID_ATOM_NUMBER != _a2;
  }

  int debug_print(std::ostream&) const;

  //  The default is to use the first two matched atoms in the query

  void
  use_defaults()
  {
    _a1 = 0;
    _a2 = 1;
  }

  Matched_Atom_Pair*
  next() const
  {
    return _next;
  }

  int number_nodes() const;

  int add(const const_IWSubstring&);

  int
  a1() const
  {
    return _a1;
  }

  int
  a2() const
  {
    return _a2;
  }
};

Matched_Atom_Pair::Matched_Atom_Pair()
{
  _a1 = INVALID_ATOM_NUMBER;
  _a2 = INVALID_ATOM_NUMBER;

  _next = nullptr;
}

Matched_Atom_Pair::~Matched_Atom_Pair()
{
  if (nullptr != _next) {
    delete _next;
  }

  return;
}

int
Matched_Atom_Pair::debug_print(std::ostream& os) const
{
  os << "  Matched_Atom_Pair ";

  if (!active()) {
    os << "inactive\n";
  } else {
    os << _a1 << " to " << _a2 << endl;
  }

  if (nullptr == _next) {
    return os.good();
  }

  return _next->debug_print(os);
}

int
Matched_Atom_Pair::add(const const_IWSubstring& s)
{
  if (!active()) {
    return _add(s);
  }

  Matched_Atom_Pair* p = this;

  while (nullptr != p->_next) {
    p = p->_next;
  }

  p->_next = new Matched_Atom_Pair();

  return p->_next->_add(s);
}

int
Matched_Atom_Pair::number_nodes() const
{
  if (nullptr == _next) {
    return 1;
  }

  return 1 + _next->number_nodes();
}

int
Matched_Atom_Pair::_add(const const_IWSubstring& s)
{
  assert(nullptr == _next);

  const_IWSubstring s1, s2;

  if (!s.split(s1, '-', s2) || 0 == s1.length() || 0 == s2.length()) {
    cerr << "Matched_Atom_Pair::_add:invalid matched atom pair specification '" << s
         << "\n";
    return 0;
  }

  if (!s1.numeric_value(_a1) || !s2.numeric_value(_a2) || _a1 == _a2 || _a1 < 0 ||
      _a2 < 0) {
    cerr << "matched::_add:invalid atom numbers '" << s << "'\n";
    return 0;
  }

  return 1;
}

/*
  We will have one matched atom pair linked list for each query
*/

static Matched_Atom_Pair* map = nullptr;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Computes verloop descriptors based on query matches\n";
  cerr << "  -B <bond>      specify bonds defining fragments\n";
  cerr << "                 <bond> is either 'aaa-bbb' meaning matched atoms aaa and bbb\n";
  cerr << "                 or 'qqq:aaa-bbb' which means query qqq\n";
  cerr << "                 NOTE THAT QUERIES AND MATCHED ATOMS START AT 0\n";
  cerr << "  -3             do 3D sterimol computations\n";
  cerr << "  -q <query>     specify query(s)\n";
  cerr << "  -s <smarts>    specify smarts\n";
  cerr << "  -w             do whole molecule sterimol computation\n";
  cerr << "  -z first       just take the first match when multiple matches found\n";
  cerr << "  -z i           ignore molecules not matching any of the queries\n";
  cerr << "  -Q ...         partial charge specification, 'abr', 'gast' or 'gh'\n";
  (void) display_standard_charge_assigner_options(cerr, 'N');
  cerr << "  -H <..>        donor acceptor assignment, enter '-H help' for info\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -h             implicit hydrogens become explicit\n";
  cerr << "  -X <element>   remove atoms of type <element>\n";
  cerr << "  -E ...         standard element options, enter '-E help' for details\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
verloop(Molecule& m, const int* don_acc, int* tmp, const atom_number_t a1,
        const atom_number_t a2, Topological_Sterimol& x2d, Sterimol& x3d)
{
  if (!x2d.compute_descriptors(m, don_acc, tmp, a1, a2)) {
    cerr << "verloop: cannot compute 2D verloop descriptors, atoms " << a1 << " and "
         << a2 << endl;
    return 0;
  }

  if (do_3d_computations) {
    if (!x3d.do_computation(m, don_acc, tmp, a1)) {
      cerr << "verloop: cannot compute 3D verloop descriptors\n";
      return 0;
    }
  }

  return 1;
}

static int
verloop(Molecule& m, const int* don_acc, int* tmp, const Set_of_Atoms& embedding,
        const Matched_Atom_Pair& mm, Topological_Sterimol& x2d, Sterimol& x3d)
{
  atom_number_t a1 = embedding[mm.a1()];
  atom_number_t a2 = embedding[mm.a2()];

  if (verbose > 2) {
    cerr << "Molecule '" << m.name() << "' with " << m.natoms() << " atoms\n";
    cerr << "Embedding " << embedding << endl;
    atom_number_t a1 = embedding[mm.a1()];
    cerr << "Atom 1 " << mm.a1() << " atom " << a1 << " which is " << m.atomic_symbol(a1)
         << " with " << m.ncon(a1) << " connections\n";
    atom_number_t a2 = embedding[mm.a2()];
    cerr << "Atom 2 " << mm.a2() << " atom " << a2 << " which is " << m.atomic_symbol(a2)
         << " with " << m.ncon(a2) << " connections\n";
  }

  assert(m.are_bonded(a1, a2));

  int rc = verloop(m, don_acc, tmp, a1, a2, x2d, x3d);

  if (rc) {
    return rc;
  }

  cerr << "Cannot compute 2D verloop descriptors\n";
  cerr << "Embedding " << embedding << " items " << mm.a1() << " and " << mm.a2() << endl;

  return 0;
}

static int
verloop(Molecule& m, const int* don_acc, int* tmp, std::ostream& output)
{
  Molecule_to_Match target(&m);

  int nq = queries.number_elements();

  // We need to keep track of which result we are currently working on

  int whichresult = 0;

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;
    int nhits = queries[i]->substructure_search(target, sresults);

    if (verbose > 1) {
      cerr << m.name() << " " << nhits << " hits to query " << i << endl;
    }

    if (0 == nhits) {
      if (ignore_molecules_not_matching) {
        continue;
      }

      cerr << "No hits to query " << i << ", only matched "
           << queries[i]->max_query_atoms_matched_in_search()
           << " query atoms, molecule '" << m.name() << "', fatal\n";
      return 0;
    }

    if (1 == nhits) {  // hopefully the most common case
      ;
    } else if (just_take_first_of_multiple_hits) {
      nhits = 1;
    } else if (multiple_hits_ok) {
      ;
    } else {
      cerr << nhits << " hits to query " << i << " '" << queries[i]->comment()
           << "' in molecule '" << m.name() << "'\n";
      return 0;
    }

    nhits = 1;  // too hard otherwise

    const Matched_Atom_Pair* mapj = &(map[i]);

    while (mapj) {
      //    cerr << "Processing result " << whichresult << endl;
      for (int k = 0; k < nhits; k++)  // doesn't do anything, nhits was set to 1
      {
        const Set_of_Atoms* s = sresults.embedding(k);
        if (!verloop(m, don_acc, tmp, *s, *mapj, results_2d[whichresult],
                     results_3d[whichresult])) {
          return 0;
        }
      }

      mapj = mapj->next();

      whichresult++;
    }
  }

  write_space_suppressed_string(m.name(), output);

  for (int i = 0; i < nbonds; i++) {
    results_2d[i].write_descriptors(output);
    if (do_3d_computations) {
      results_3d[i].write_descriptors(output);
    }
  }

  if (do_whole_molecule_sterimol) {
    results_3d[nbonds].do_computation(m);
    results_3d[nbonds].write_descriptors(output);
  }

  output << endl;

  return output.good();
}

static void
preprocess(Molecule& m)
{
  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (elements_to_remove.active()) {
    elements_to_remove.process(m);
  }

  if (stream_for_labelled_molecules.active()) {
    stream_for_labelled_molecules.write(m);
  }

  return;
}

static int
verloop(Molecule& m, std::ostream& output)
{
  (void)m.ring_membership();  // always needed

  preprocess(m);

  if (charge_assigner.active()) {
    charge_assigner.process(m);
  }

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  const int matoms = m.natoms();

  int* don_acc = new_int(matoms + matoms);
  std::unique_ptr<int[]> free_don_acc(don_acc);
  int* tmp = don_acc + matoms;

  if (donor_acceptor_assigner.active()) {
    donor_acceptor_assigner.process(m, don_acc);
  }

  return verloop(m, don_acc, tmp, output);
}

static int
verloop(data_source_and_type<Molecule>& input, std::ostream& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (1 == molecules_read && do_3d_computations) {
      if (m->highest_coordinate_dimensionality() < 3) {
        cerr << "3D computations requested, but molecule does not have 3D coordinates\n";
        return 0;
      }
    }

    if (!verloop(*m, output)) {
      return 0;
    }
  }

  return 1;
}

static int
verloop(const char* fname, FileType input_type, std::ostream& output)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  return verloop(input, output);
}

static int
process_query(const_IWSubstring& s)
{
  Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;

  if (!q->create_from_smarts(s)) {
    cerr << "Cannot parse smarts '" << s << "'\n";
    return 0;
  }

  queries.add(q);

  return 1;
}

/*
  A bond specification can be either
  qqq:aaa-bbb
  or just
  aaa-bbb

  'aaa' and 'bbb' are the matched atom numbers. If present, 'qqq' is the query
  number


*/

static int
parse_dash_b_option(const const_IWSubstring& b)
{
  if (!b.contains(':')) {
    return map[0].add(b);
  }

  int query_number;

  const_IWSubstring s1, s2;

  if (!b.split(s1, ':', s2) || 0 == s1.length() || 0 == s2.length()) {
    cerr << "The -B optin must be of the form nn:ii-jj, '" << b << "' invalid\n";
    return 0;
  }

  if (!s1.numeric_value(query_number) || query_number < 0 ||
      query_number >= queries.number_elements()) {
    cerr << "Invalid query number specification '" << b << "', have "
         << queries.number_elements() << " queries\n";
    return 0;
  }

  return map[query_number].add(s2);
}

static int
verloop(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vq:A:E:i:V:H:N:B:3hs:X:wz:Q:L:o:g:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose)) {
      cerr << "Cannot process aromaticity options (-A)\n";
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      return 8;
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose > 1, 'N')) {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(33);
    }
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose > 1)) {
      cerr << "Cannot process donor acceptor (-H) option(s)\n";
      usage(12);
    }
    donor_acceptor_assigner.set_apply_isotopic_labels(0);
  }

  if (cl.option_present('X')) {
    if (!elements_to_remove.construct_from_command_line(cl, verbose, 'X')) {
      cerr << "Cannot process element removal specification (-X)\n";
      usage(5);
    }
  }

  if (cl.option_present('3')) {
    do_3d_computations = 1;
    if (verbose) {
      cerr << "Will perform 3D computations\n";
    }
  }

  if (cl.option_present('w')) {
    do_whole_molecule_sterimol = 1;
    do_3d_computations = 1;

    if (verbose) {
      cerr << "Will do whole molecule sterimol compuation - implies 3D\n";
    }
  }

  if (cl.option_present('z')) {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++)) {
      if ("first" == z || 'f' == z) {
        just_take_first_of_multiple_hits = 1;
        if (verbose) {
          cerr << "Will process just the first match of multiple hits\n";
        }
      } else if ("i" == z) {
        ignore_molecules_not_matching = 1;
        if (verbose) {
          cerr << "Will ignore molecules not matching any queries\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(11);
      }
    }
  }

  if (cl.option_present('h')) {
    make_implicit_hydrogens_explicit = 1;
    if (verbose) {
      cerr << "Implicit Hydrogens will be converted to explicit Hydrogens\n";
    }

    if (cl.option_present('3')) {
      cerr << "WARNING:required placement of Hydrogens in 3D. Not reliable\n";
    }
  }

  if (cl.option_present('Q')) {
    const_IWSubstring q = cl.string_value('Q');

    if (q.starts_with("abr")) {
      do_partial_charge_computation = PARTIAL_CHARGE_ABRAHAM;
      if (verbose) {
        cerr << "Will use Abraham partial charges\n";
      }
    } else if ("gast" == q) {
      do_partial_charge_computation = PARTIAL_CHARGE_GASTEIGER;

      if (verbose) {
        cerr << "Will use Gasteiger partial charges\n";
      }
    } else if ("gh" == q) {
      do_partial_charge_computation = PARTIAL_CHARGE_GH;

      if (verbose) {
        cerr << "Will use Gasteiger-Huckel partial charges\n";
      }
    } else {
      cerr << "Unrecognised partial charge specification '" << q << "'\n";
    }
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify one or more queries via the -q or -s options\n";
    usage(4);
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose, 'q')) {
      cerr << "Cannot read query specification(s) (-q option)\n";
      usage(5);
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    int i = 0;
    while (cl.value('s', s, i++)) {
      if (!process_query(s)) {
        cerr << "Cannot process smarts '" << s << "'\n";
        return 15;
      }
    }
  }

  int nq = queries.number_elements();

  if (0 == nq) {
    cerr << "No queries. Use the -q option\n";
    usage(6);
  }

  map = new Matched_Atom_Pair[nq];
  std::unique_ptr<Matched_Atom_Pair[]> free_map(map);
  if (nullptr == map) {
    cerr << "Memory failure, cannot allocate " << nq << " matched atom pair objects\n";
    return 41;
  }

  for (int i = 0; i < nq; i++) {
    queries[i]->set_find_unique_embeddings_only(1);
    queries[i]->set_find_one_embedding_per_atom(1);  // should this be optional?
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    int i = 0;
    while (cl.value('B', b, i++)) {
      if (!parse_dash_b_option(b)) {
        cerr << "Cannot parse -B option '" << b << "'\n";
        return 14;
      }
    }
  }

  // Some queries may have had matched atoms specified via the -B option. For
  // others, we need to initialise them

  nbonds = 0;  // the total number of Verloop bonds defined

  for (int i = 0; i < nq; i++) {
    Matched_Atom_Pair& m = map[i];

    if (!m.active()) {  // must have been done via the -B option
      m.use_defaults();
    }

    nbonds += m.number_nodes();
  }

  if (0 == nbonds) {
    cerr << "Yipes, no bonds defined!!!!\n";
    return 17;
  }

  if (verbose) {
    cerr << nbonds << " in " << nq << " queries defined\n";
  }

  if (verbose > 1) {
    for (int i = 0; i < nq; i++) {
      cerr << "Query " << i << " '" << queries[i]->comment() << "'\n";
      map[i].debug_print(cerr);
    }
  }

  if (cl.option_present('L')) {
    const_IWSubstring l = cl.string_value('L');

    if (cl.option_present('o')) {
      if (!stream_for_labelled_molecules.determine_output_types(cl)) {
        cerr << "Cannot discern output type(s) for -L file\n";
        return 6;
      }
    } else if (do_3d_computations) {
      stream_for_labelled_molecules.add_output_type(FILE_TYPE_SDF);
    } else {
      stream_for_labelled_molecules.add_output_type(FILE_TYPE_SMI);
    }

    if (stream_for_labelled_molecules.would_overwrite_input_files(cl, l)) {
      cerr << "Cannot overwrite input file(s) via the -L option, '" << l
           << "' is invalid\n";
      return 4;
    }

    if (!stream_for_labelled_molecules.new_stem(l)) {
      cerr << "Cannot open stream for labelled molecules '" << l << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Labelled molecules written to '" << l << "'\n";
    }
  }

  //  Allocate 2D and 3D results

  results_2d = new Topological_Sterimol[nbonds];
  std::unique_ptr<Topological_Sterimol[]> free_2d(results_2d);
  results_3d = new Sterimol[nbonds + do_whole_molecule_sterimol];
  std::unique_ptr<Sterimol[]> free_3d(results_3d);

  std::cout << "name";

  for (int i = 0; i < nbonds; i++) {
    if (do_partial_charge_computation) {
      results_2d[i].set_do_partial_charge_descriptors(1);
    }

    IWString vi;
    vi = "v";
    vi.append_number(i);

    std::cout << ' ' << vi << "natoms";
    std::cout << ' ' << vi << "hetero";
    std::cout << ' ' << vi << "donor";
    std::cout << ' ' << vi << "accept";
    std::cout << ' ' << vi << "nrings";
    std::cout << ' ' << vi << "arom";
    std::cout << ' ' << vi << "mxdist";
    std::cout << ' ' << vi << "mxaccd";
    std::cout << ' ' << vi << "mxdond";
    std::cout << ' ' << vi << "npos";
    std::cout << ' ' << vi << "nneg";
    std::cout << ' ' << vi << "nunsat";
    std::cout << ' ' << vi << "cigar";
    std::cout << ' ' << vi << "avd10";
    std::cout << ' ' << vi << "atoms1";
    std::cout << ' ' << vi << "atoms2";
    std::cout << ' ' << vi << "atoms3";
    std::cout << ' ' << vi << "minarom";
    std::cout << ' ' << vi << "maxarom";
    std::cout << ' ' << vi << "minhtro";
    std::cout << ' ' << vi << "maxhtro";

    if (do_partial_charge_computation) {
      std::cout << ' ' << vi << "q0";
      std::cout << ' ' << vi << "q1";
      std::cout << ' ' << vi << "qmax";
      std::cout << ' ' << vi << "qmin";
      std::cout << ' ' << vi << "qext";
      std::cout << ' ' << vi << "tabs";
    }

    if (do_3d_computations) {
      write_sterimol_descriptor_headers(std::cout, vi, nq);
      if (do_partial_charge_computation) {  // should implement this sometime
        ;
      }
    }
  }

  if (do_whole_molecule_sterimol) {
    write_sterimol_descriptor_headers(std::cout, "vw");
  }

  std::cout << endl;

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!verloop(cl[i], input_type, std::cout)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    if (verbose > 1) {
      for (int i = 0; i < queries.number_elements(); i++) {
        queries[i]->report(cerr, verbose);
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = verloop(argc, argv);

  return rc;
}
