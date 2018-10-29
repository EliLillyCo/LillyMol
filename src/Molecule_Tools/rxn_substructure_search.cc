/*
  Substructure Searching over Reactions
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"

#include "istream_and_type.h"
#include "rwsubstructure.h"
#include "misc2.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "substructure.h"
#include "target.h"
#include "path.h"

const char * prog_name = NULL;

static int verbose = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int reactions_with_not_enough_components_for_queries = 0;

static int reactions_read = 0;

static int reactions_matching = 0;

static int single_query_search_anywhere = 0;

static IWString_and_File_Descriptor stream_for_matches;
static IWString_and_File_Descriptor stream_for_non_matches;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Substructure Searching over reactions\n";
  cerr << "  -q <q>        qry+qry+qry>qry>qry+qry\n";
  cerr << "                where qry can be\n";
  cerr << "                S:fname     file of smarts\n";
  cerr << "                Q:fname     file of query file names\n";
  cerr << "                q:fname     single query file\n";
  cerr << "                            otherwise interpreted as smarts\n";
  cerr << "  -b .          if single query, allow it to match anywhere in reagents/agents/products\n";
  cerr << "  -m <fname>    write reactions that        match to <fname>\n";
  cerr << "  -n <fname>    write reactions that do not match to <fname>\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

class Reaction_to_Search
{
  private:
    resizable_array_p<Molecule> _reagent;
    resizable_array_p<Molecule> _agent;
    resizable_array_p<Molecule> _product;

    IWString _name;

//  private functions

    int _build_set_of_molecules(const_IWSubstring s, resizable_array_p<Molecule> & m);

  public:
    Reaction_to_Search();
    ~Reaction_to_Search();

    int build(const const_IWSubstring & buffer);

    int nr() const { return _reagent.number_elements();}
    int na() const { return _agent.number_elements();}
    int np() const { return _product.number_elements();}

    resizable_array_p<Molecule> & reagent() { return _reagent;}
    resizable_array_p<Molecule> & agent() { return _agent;}
    resizable_array_p<Molecule> & product() { return _product;}
};

Reaction_to_Search::Reaction_to_Search()
{
  return;
}

Reaction_to_Search::~Reaction_to_Search()
{
  return;
}

int
Reaction_to_Search::build(const const_IWSubstring & buffer)
{
  const_IWSubstring s, n;
  if (! buffer.split(s, ' ', n) || 0 == s.length())
  {
    cerr << "Reaction_to_Search::build:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (2 != s.ccount('>'))
  {
    cerr << "Reaction_to_Search::build:reaction must have three components separated by >, " << buffer << " invalid\n";
    return 0;
  }

  _name = n;

  const_IWSubstring r, a, p;

  int i = 0;
  s.nextword_single_delimiter(r, i, '>');
  s.nextword_single_delimiter(a, i, '>');
  s.nextword_single_delimiter(p, i, '>');

  if (! _build_set_of_molecules(r, _reagent) || ! _build_set_of_molecules(a, _agent) || ! _build_set_of_molecules(p, _product))
  {
    cerr << "Reaction_to_Search::build:invalid input '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Reaction_to_Search::_build_set_of_molecules(const_IWSubstring s,
                                            resizable_array_p<Molecule> & m)
{
  if (0 == s.length())    // nothing to do
    return 1;

  resizable_array<int> pos;

  identify_plus_positions(s, pos);
  pos.add(s.length() - 1);

  const int n = pos.number_elements();

  for (int i = 0; i < n; ++i)
  {
    int cstart;          // bunch of code lifted from build_molecules in rxnfile2.cc
    if (0 == i)
      cstart = 0;
    else
      cstart = pos[i-1] + 1;

    int cstop;

    if (i == n - 1)
      cstop = s.length() - 1;
    else
      cstop = pos[i] - 1;

    const_IWSubstring smiles;

    s.from_to(cstart, cstop, smiles);

    Molecule * t = new Molecule();

    if (! t->build_from_smiles(smiles))
    {
      cerr << "build_set_of_molecules:invalid smiles " << s << "\n";
      delete t;
      return 0;
    }

    m.add(t);
  }

  return n;
}

class RXN_Substructure_Search
{
  private:
    resizable_array_p<Substructure_Query> * _reagent;
    resizable_array_p<Substructure_Query> * _agent;
    resizable_array_p<Substructure_Query> * _product;

    int _nr, _na, _np;

    IWString _name;

//  private functions

    int _build_set_of_queries(const const_IWSubstring & q, resizable_array_p<Substructure_Query> * & queries, int & n);
    int _substructure_search(resizable_array_p<Substructure_Query> * q, const int nq, resizable_array_p<Molecule> & m);

  public:
    RXN_Substructure_Search();
    ~RXN_Substructure_Search();

    int build(const const_IWSubstring & buffer);

    int substructure_search(Reaction_to_Search & to_search);
};

RXN_Substructure_Search::RXN_Substructure_Search()
{
  _reagent = nullptr;
  _nr = 0;
  _agent = nullptr;
  _na = 0;
  _product = nullptr;
  _np = 0;

  return;
}

RXN_Substructure_Search::~RXN_Substructure_Search()
{
  if (nullptr != _reagent)
    delete [] _reagent;
  if (nullptr != _agent)
    delete [] _agent;
  if (nullptr != _product)
    delete [] _product;

  return;
}

int
RXN_Substructure_Search::build(const const_IWSubstring & buffer)
{
  const_IWSubstring query;

  int i = 0;
  if (! buffer.nextword(query, i))
  {
    cerr << "RXN_Substructure_Search::build:empty input\n";
    return 0;
  }

  buffer.nextword(_name, i);

  if (2 != query.ccount('>'))
  {
    cerr << "RXN_Substructure_Search::build:query must have three components separated by > chars. '" << query << "' invalid\n";
    return 0;
  }

  i = 0;

  const_IWSubstring r, a, p;

  query.nextword_single_delimiter(r, i, '>');
  query.nextword_single_delimiter(a, i, '>');
  query.nextword_single_delimiter(p, i, '>');

  if (! _build_set_of_queries(r, _reagent, _nr) || ! _build_set_of_queries(a, _agent, _na) || ! _build_set_of_queries(p, _product, _np))
  {
    cerr << "RXN_Substructure_Search::build:invalid input '" << query << "'\n";
    return 0;
  }

  return 1;
}



int
RXN_Substructure_Search::_build_set_of_queries(const const_IWSubstring & q,
                                               resizable_array_p<Substructure_Query> * & queries,
                                               int & n)
{
  if (0 == q.length())    // nothing to do
    return 1;

  resizable_array<int> pos;

  identify_plus_positions(q, pos);
  pos.add(q.length() - 1);

  n = pos.number_elements();

  queries = new resizable_array_p<Substructure_Query>[n];

  for (int i = 0; i < n; ++i)
  {
    int cstart;          // bunch of code lifted from build_molecules in rxnfile2.cc
    if (0 == i)
      cstart = 0;
    else
      cstart = pos[i-1] + 1;

    int cstop;

    if (i == n - 1)
      cstop = q.length() - 1;
    else
      cstop = pos[i] - 1;

    const_IWSubstring s;

    q.from_to(cstart, cstop, s);

    int rc = 1;
    if (s.starts_with("S:"))
    {
      s.remove_leading_chars(2);
      rc = smarts_from_file(s, queries[i], verbose);
    }
    else if (s.starts_with("Q:") || s.starts_with("F:"))
    {
      s.remove_leading_chars(2);
      rc = queries_from_file(s, queries[i], 1, verbose);
    }
    else if (s.starts_with("q:"))
    {
      s.remove_leading_chars(2);
      rc = read_one_or_more_queries_from_file(queries[i], s, verbose);
    }
    else
    {
      Substructure_Query * q = new Substructure_Query();
      if (! q->create_from_smarts(s))
      {
        cerr << "RXN_Substructure_Search::build:invalid smarts '" << s << "'\n";
        delete q;
        return 0;
      }

      queries[i].add(q);
    }
  }

  return n;
}

#ifdef OLD_STUFFQWEQWE
int
RXN_Substructure_Search::substructure_search(const const_IWSubstring & buffer)
{
  const_IWSubstring q, n;
  if (! buffer.split(q, ' ', n) || 0 == q.length())
  {
    cerr << "RXN_Substructure_Search::substructure_search:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (2 != q.ccount('>'))
  {
    cerr <<
  }

  const_IWSubstring r, a, p;

  int i = 0;
  q.nextword_single_delimiter(r, i, '>');
  q.nextword_single_delimiter(a, i, '>');
  q.nextword_single_delimiter(p, i, '>');

  Reaction_to_Search to_search;

  if (! to_search.build(buffer))
  {
    cerr << "RXN_Substructure_Search::substructure_search:invalid input '" << buffer << "'\n";
    return 0;
  }

  return _substructure_search(to_search);
}
#endif

int
RXN_Substructure_Search::substructure_search(Reaction_to_Search & to_search)
{
  if (! _substructure_search(_reagent, _nr, to_search.reagent()))
    return 0;

  if (! _substructure_search(_agent, _na, to_search.agent()))
    return 0;

  if (! _substructure_search(_product, _np, to_search.product()))
    return 0;

  return 1;
}

//#define DEBUG_SS_RXN_MATCH

static int
matches(resizable_array_p<Substructure_Query> & q,
        Molecule & m)
{
#ifdef DEBUG_SS_RXN_MATCH
  cerr << "Searching " << m.smiles() << " with " << q.number_elements() << " queries\n";
#endif

  Molecule_to_Match target(&m);

  const int nq = q.number_elements();

  for (int i = 0; i < nq; ++i)
  {
    Substructure_Results sresults;

    if (q[i]->substructure_search(target, sresults))
      return 1;

#ifdef DEBUG_SS_RXN_MATCH
    cerr << " query " << i << " only matched " << sresults.max_query_atoms_matched_in_search() << " query atoms\n";
#endif
  }

  return 0;
}

int
RXN_Substructure_Search::_substructure_search(resizable_array_p<Substructure_Query> * q,
                                              const int nq,
                                              resizable_array_p<Molecule> & m)
{
  const int nm = m.number_elements();

// Special case of 1 query, but multiple reagents/products.

  if (single_query_search_anywhere && 1 == nq)
  {
    for (int i = 0; i < nm; ++i)
    {
      if (matches(q[0], *m[i]))
        return 1;
    }

    return 0;
  }

  static const char * role[] = {"reagent", "agent", "product"};

  if (nq <= nm)      // fewer queries than components in the reaction is easy
  {
    for (int i = 0; i < nq; ++i)
    {
      if (matches(q[i], *m[i]))
        continue;

      if (verbose > 1)
        cerr << role[i] << " no match, " << q[i].number_elements() << " queries\n";

      return 0;
    }

    return 1;
  }

// More queries than molecules to search, do not know what to do...

  if (verbose > 1)
    cerr << _name << " cannot search " << nq << " queries against " << nm << " molecules\n";

  reactions_with_not_enough_components_for_queries++;

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
rxn_substructure_search (iwstring_data_source & input,
                         RXN_Substructure_Search & query,
                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    reactions_read++;

    Reaction_to_Search to_search;

    if (! to_search.build(buffer))
    {
      cerr << "Cannot read reaction " << buffer << endl;
      return 0;
    }

    if (query.substructure_search(to_search))
    {
      reactions_matching++;

      if (stream_for_matches.is_open())
      {
        stream_for_matches << buffer << '\n';
        stream_for_matches.write_if_buffer_holds_more_than(4096);
      }
    }
    else
    {
      if (stream_for_non_matches.is_open())
      {
        stream_for_non_matches << buffer << '\n';
        stream_for_non_matches.write_if_buffer_holds_more_than(4096);
      }
    }
  }

  return 1;
}

static int
rxn_substructure_search (const char * fname,
                         RXN_Substructure_Search & query,
                         IWString_and_File_Descriptor & output)
{
  assert(NULL != fname);

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return rxn_substructure_search(input, query, output);
}

static int
rxn_substructure_search (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:g:lq:n:m:b:");

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

  if (! cl.option_present('q'))
  {
    cerr << "Must specify query via the -q option\n";
    usage(1);
  }

  RXN_Substructure_Search rxn_query;

  if (cl.option_present('q'))
  {
    const_IWSubstring q = cl.string_value('q');

    if (! rxn_query.build(q))
    {
      cerr << "Cannot build reaction query '" << q << "'\n";
      return 1;
    }
  }

  if (cl.option_present('b'))
  {
    single_query_search_anywhere = 1;

    if (verbose)
      cerr << "A single query can match anywhere in reagents/agents/products\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int output_produced = 0;

  if (cl.option_present('m'))
  {
    IWString tmp = cl.string_value('m');

    if (! tmp.ends_with(".rxnsmi"))
      tmp << ".rxnsmi";

    if (! stream_for_matches.open(tmp.null_terminated_chars()))
    {
      cerr << "Cannot open stream for matches '" << tmp << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Reactions matching written to '" << tmp << "'\n";

    output_produced++;
  }

  if (cl.option_present('n'))
  {
    IWString tmp = cl.string_value('n');

    if (! tmp.ends_with(".rxnsmi"))
      tmp << ".rxnsmi";

    if (! stream_for_non_matches.open(tmp.null_terminated_chars()))
    {
      cerr << "Cannot open stream for non matches '" << tmp << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Reactions not matching written to '" << tmp << "'\n";

    output_produced++;
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! rxn_substructure_search(cl[i], rxn_query, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose || ! output_produced)
  {
    cerr << "Read " << reactions_read << " reactions, " << reactions_matching << " matched\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rxn_substructure_search(argc, argv);

  return rc;
}
