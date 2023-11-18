/*
  A server for doing gfp searches
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <memory>
#include <fstream>
#include <queue>

//#include "boost/heap/heap_concepts.hpp"

#include "zmq.hpp"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#define IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwbits/iwbits.h"

#include "gfp_standard.h"

#include "mol/molecule.h"
#include "mol/mpr.h"
#include "mol/iwmfingerprint.h"
#include "mol/maccskeys_fn5.h"
#include "mol/iwstandard.h"
#include "mol/aromatic.h"

const char * prog_name = nullptr;

static int verbose = 0;

//static std::unique_ptr<re2::RE2> fingerprint_rx("^\\|");

static int searches_done = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

#define RELOAD_POOL "RELOAD"
#define TERMINATE_PROCESSING "printf"    // so it looks like a regular string in the binary

static Chemical_Standardisation chemical_standardisation;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Server for molecular similarity\n";
  cerr << " -p <port>      port for both input (smiles or gfp) and results\n";
  cerr << " -i             also read the smiles\n";
  cerr << " -s <int>       size of the input fingerprint file\n";
  cerr << " -f             fault tolerance mode\n";
  cerr << " -L <fname>     file name for logging\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
global_logging (const IWString & s)
{
  return 1;
}

static int
global_error (const IWString & s)
{
  global_logging(s);

  return 0;
}

/*
  We need a lightweight structure for holding distance/id information.
*/

struct Distance_ID
{
  float _d;
  int _ndx;    // index into pool
};

class Compare_DID
{
  private:
  public:
    bool operator () (const Distance_ID & lhs, const Distance_ID & rhs) const
    {
      return lhs._d > rhs._d;
    }
};

static void
preprocess (Molecule & m)
{
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  chemical_standardisation.process(m);

  return;
}

class GFP_Server
{
  private:
    int _pool_size;
    GFP_Standard * _pool;

//  We need a structure for accumulating results

    Distance_ID * _did;

    IWString _output_buffer;

    int _fault_tolerant;
    std::ofstream _stream_for_logging;
    bool _stream_for_logging_open;

    int _input_port, _output_port;

    Molecular_Properties_Generator _mpr;
    MACCSKeys _mk;

    int _tmp[2048];    // used by maccskeys and iwfp

//  the names of the fingerprints in the pool

    IWString * _id;

    int _read_smiles;

    IWString * _smiles;

    unsigned int _searches_done;

//  private functions

    void _free_arrays ();

    void _log_message(const IWString &);
    void _error_exit (const IWString &);

    int _insert_smiles_into_output_buffer (const IW_TDT & tdt);
    int _insert_smiles_into_output_buffer (const IWString & s);
    int _fill_output_buffer_failure ();
    int _fill_output_buffer (int);
    int _tanimoto_single_nbr(const GFP_Standard &);
    int _tanimoto_within_distance(const GFP_Standard &, float);
    int _tanimoto (const GFP_Standard & sfp, int nbrs);
    int _extract_number_nbrs_tdt (IWString &, int &, float &);
    int _extract_number_nbrs_smiles (IWString &, int &, float &);
    int _rebuild_pool (const IWString &);
    int _do_shutdown ();

  public:
    GFP_Server();
    ~GFP_Server();

    int parse_command_line (Command_Line & cl);

    int build (const char * fname);
    int build (iwstring_data_source &);

    int tanimoto_from_tdt (IWString &);
    int tanimoto_from_smiles (IWString &);

    void doit ();
};

GFP_Server::GFP_Server()
{
  _pool_size = 0;
  _pool = nullptr;
  _fault_tolerant = 0;

  _did = new Distance_ID[1];    // the most common case

  _input_port = 0;
  _output_port = 0;

  _stream_for_logging_open = false;

  _id = nullptr;
  _read_smiles = 0;
  _smiles = nullptr;

  _searches_done = 0;

  return;
}

void
GFP_Server::_free_arrays ()
{
  if (NULL != _pool)
  {
    delete [] _pool;
    _pool = nullptr;
  }

  _pool_size = 0;

  if (NULL != _did)
  {
    delete [] _did;
    _did = nullptr;
  }

  if (NULL != _id)
  {
    delete [] _id;
    _id = nullptr;
  }

  if (NULL != _smiles)
  {
    delete [] _smiles;
    _smiles = nullptr;
  }

  return;
}

GFP_Server::~GFP_Server()
{
  _free_arrays();

  return;
}

void
GFP_Server::_log_message (const IWString & msg)
{
  if (! _stream_for_logging_open)
    return;

  _stream_for_logging << msg;

  if (! msg.ends_with('\n'))
    _stream_for_logging << '\n';

  _stream_for_logging.flush();

  return;
}

void
GFP_Server::_error_exit (const IWString & msg)
{
  _log_message(msg);

  _do_shutdown();

  exit(1);
}

int
GFP_Server::parse_command_line (Command_Line & cl)
{
  if (! cl.option_present('p'))
  {
    cerr << "Must specify port to use via the -p option\n";
    return 0;
  }
  else if (! cl.value('p', _input_port) || _input_port < 30)
  {
    cerr << "Invalid input port specification (-i)\n";
    return 0;
  }
  else if (verbose)
  {
    IWString msg;
    msg << "Input from port " << _input_port;
    cerr << msg << "\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', _pool_size) || _pool_size < 1)
    {
      cerr << "The pool size must be a whole +ve number\n";
      return 0;
    }

    if (verbose)
      cerr << "Fingerprint pool explicitly sized to " << _pool_size << " fingerprints\n";
  }

  if (cl.option_present('f'))
  {
    _fault_tolerant = 1;

    if (verbose)
      cerr << "Will run in fault tolerant mode\n";
  }

  if (cl.option_present('L'))
  {
    const char * l = cl.option_value('L');

    _stream_for_logging.open(l, std::fstream::out);

    if (! _stream_for_logging.good())
    {
      cerr << "GFP_Server::parse_command_line:cannot initialise logging stream '" << l << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Logging messages written to '" << l << "'\n";

    _stream_for_logging_open = true;
  }

  if (cl.option_present('i'))
  {
    _read_smiles = 1;

    if (verbose)
      cerr << "Will include smiles in memory\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "GFP_Server::parse_command_line:no fingerprint file specified\n";
    return 0;
  }

  return build(cl[0]);
}

int
GFP_Server::build (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    IWString msg;
    msg << "GFP_Server::build:cannot open '" << fname << "'";
    _error_exit(msg);
  }

  return build (input);
}

static int
build_gfp_standard (GFP_Standard & sfp,
                    IW_General_Fingerprint & gfp)
{
  sfp.build_molecular_properties(gfp.molecular_properties_integer());
  sfp.build_iw(gfp[0]);
  sfp.build_mk(gfp[1]);
  sfp.build_mk2(gfp[2]);

  return 1;
}

static int
build_gfp_standard (GFP_Standard & sfp,
                    IWString & id,
                    IW_TDT & tdt)
{
  int fatal;
  IW_General_Fingerprint gfp;

  if (! gfp.construct_from_tdt(tdt, fatal))
    return 0;

  id = gfp.id();

  return build_gfp_standard(sfp, gfp);
}

static int
build_gfp_standard (GFP_Standard & sfp,
                    IW_TDT & tdt)
{
  int fatal;
  IW_General_Fingerprint gfp;

  if (! gfp.construct_from_tdt(tdt, fatal))
    return 0;

  return build_gfp_standard (sfp, gfp);
}

int
GFP_Server::build (iwstring_data_source & input)
{
  if (0 == _pool_size)
  {
    _pool_size = input.count_records_starting_with("|");
//  cerr << "Pool size " << _pool_size << endl;
  }

  if (0 == _pool_size)
  {
    IWString msg;
    msg << "GFP_Server::build:no fingerprints in file\n";
    _error_exit(msg);
  }

  if (verbose)
  {
    IWString msg;
    msg << "GFP_Server::build:fingerprint file contains " << _pool_size << " fingerprints";
    _log_message(msg);
  }

  _pool = new GFP_Standard[_pool_size];
  _did = new Distance_ID[_pool_size];
  _id = new IWString[_pool_size];

  if (_read_smiles)
    _smiles = new IWString[_pool_size];

  if (NULL == _pool || NULL == _did || NULL == _id)
  {
    IWString msg;
    msg << "GFP_Server::build:cannot allocate " << _pool_size << " fingerprints";
    _error_exit(msg);
  }

  IWString msg;
  msg << "Begin building pool size " << _pool_size;
  _log_message(msg);

  for (auto i = 0; i < _pool_size; ++i)
  {
    IW_TDT tdt;
    if (! tdt.next(input))
    {
      IWString msg;
      msg << "GFP_Server::build:premature, expected " << _pool_size << " got " << i;
      _log_message(msg);
      _pool_size = i - 1;
      break;
    }

    if (! build_gfp_standard (_pool[i], _id[i], tdt))
    {
      IWString msg;
      msg << "GFP_Server::build:fatal error reading fingerprint " << i;
      _error_exit(msg);
    }

    if (NULL != _smiles)
      tdt.dataitem_value(smiles_tag, _smiles[i]);
  }

  return _pool_size;
}

/*
  We ignore errors where we can just leave the existing pool in place
*/

int
GFP_Server::_rebuild_pool (const IWString & s)
{
  IWString notused, fname;
  if (! s.split(notused, ' ', fname))
  {
    IWString msg;
    msg << "GFP_Server::_rebuild_pool:invalid directive '" << s << "'\n";
    _error_exit(msg);
    return 1;
  }

  iwstring_data_source input(fname.null_terminated_chars());

  if (! input.good())
  {
    IWString msg;
    msg << "GFP_Server::_rebuild_pool:cannot open '" << fname << "'\n";
    _error_exit(msg);
    return 1;
  }

  int nps = input.count_records_starting_with("|");
  if (nps < (9 * _pool_size / 10))     // pool cannot shrink by more than 10%
  {
    IWString msg;
    msg << "GFP_Server::_rebuild_pool:new fingerprint file too small '" << fname << "' count " << nps;
    _error_exit(msg);
    return 1;
  }

  _free_arrays();

  _pool_size = nps;

  if (! build (input))
  {
    IWString msg("Cannot rebuilt build");
    _error_exit(msg);
  }

  _output_buffer = fname << " reloaded.";

  return 1;
}

int
GFP_Server::_do_shutdown ()
{
  _output_buffer << "ACK";

  IWString msg;
  msg << "Shutdown, performed " << _searches_done << " searches\n";
  _log_message(msg);

  return 1;
}

int
GFP_Server::_fill_output_buffer (int nbrs)
{
#ifdef DEBUG_SERVER
  cerr << "Server filling output buffer for " << nbrs << " nbrs\n";
#endif

  for (auto i = 0; i < nbrs; ++i)
  {
    const Distance_ID & didi = _did[i];

    if (NULL != _smiles)
      _output_buffer << _smiles[didi._ndx] << ' ';

    _output_buffer << _id[didi._ndx] << ' ' << didi._d << "\n";
  }

  return 1;
}

int
GFP_Server::_fill_output_buffer_failure ()
{
  _output_buffer.resize_keep_storage(0);

  _output_buffer << "-1.0 -1\n";

  return 0;     // always
}

/*
  Our convention is that after the vertical bar there will be either
  1. nothing, return single nearest nbr
  2. an integer - return this many nbrs
  3. a distance - return all nbrs within that distance
*/

int
GFP_Server::_extract_number_nbrs_tdt (IWString & s,
                     int & nbrs,
                     float & dist)
{
  dist = -1.0f;

  if (s.ends_with("|\n"))
  {
    nbrs = 1;
    return 1;
  }

// Either a digit or a fraction has been added

  nbrs = -1;

  int last_newline = s.rindex('\n');

  if (last_newline > 0)
    ;
  else if (_fault_tolerant)
  {
    IWString msg;
    msg << "GFP_Server::_extract_number_nbrs:no vbar in TDT input '" << s << "'";
    _log_message(msg);
    return 0;
  }
  else
  {
    IWString msg;
    msg << "GFP_Server::_extract_number_nbrs:no vbar in TDT input '" << s << "'";
    _error_exit(msg);
  }

  const_IWSubstring token;
  s.from_to(last_newline + 1, s.length() - 1, token);
  s.iwtruncate(last_newline);
  s << '\n';

#ifdef DEBUG_SERVER
  cerr << "Nbr information contained in '" << token << "'\n";
#endif

  if (token.contains('.'))
  {
    if (token.numeric_value(dist) || dist < 0.0f)
      ;
    else
    {
      IWString msg;
      msg << "GFP_Server::_extract_number_nbrs:invalid distance '" << s << "'";
      if (_fault_tolerant)
      {
        _log_message(msg);
        return 0;
      }
      else
        _error_exit(msg);
    }
  }
  else if (token.numeric_value(nbrs) && nbrs > 0)
  {
  }
  else
  {
    IWString msg;
    msg << "GFP_Server::_extract_number_nbrs:invalid neighbour specification '" << s << "'";
    if (_fault_tolerant)
    {
      _log_message(msg);
      return 0;
    }
    else
      _error_exit(msg);
  }

  return 1;
}

/*
  Smiles consists of smiles id <token>
*/

int
GFP_Server::_extract_number_nbrs_smiles (IWString & s, int & nbrs, float & dist)
{
  dist = -1.0f;

  const_IWSubstring token;
  int i = 0;

  if (s.nextword(token, i) && s.nextword(token, i))
    ;
  else if (_fault_tolerant)
    return 0;
  else
  {
    IWString msg;
    msg << "GFP_Server::_extract_number_nbrs_smiles:invalid input '" << s << "'\n";
    _error_exit(msg);
  }

  if (! s.nextword(token, i))
  {
    nbrs = 1;
    return 1;
  }

  if (token.ends_with('\n'))
    token.chop();

  int rc = 1;

  if (token.contains('.'))
  {
    nbrs = 0;
    if (token.numeric_value(dist) && dist > 0.f && dist < 1.0f)
      return 1;
  }
  else if (token.numeric_value(nbrs) && nbrs > 0)
    return 1;

// Something not right
  if (_fault_tolerant)
  {
    nbrs = 1;
    return 1;
  }

  IWString msg;
  msg << "GFP_Server::_extract_number_nbrs_smiles:invalid nbrs '" << s << "'\n";
  _error_exit(msg);
}

int
GFP_Server::tanimoto_from_tdt (IWString & s)
{
  searches_done++;

  int nbrs = -1;
  float dist = -1.0f;

  if (_extract_number_nbrs_tdt (s, nbrs, dist))
    ;
  else if (_fault_tolerant)
    return _fill_output_buffer_failure ();
  else
  {
    IWString msg;
    msg << "GFP_Server::tanimoto:cannot discern number of nbrs in '" << s << "'\n";
    _error_exit(msg);
  }

  IW_TDT tdt;

  if (! tdt.build(s))
  {
    IWString msg;
    msg << "GFP_Server::tanimoto:cannot build tdt from '" << s << "'";
    _error_exit(msg);
  }

  GFP_Standard gfp;

  if (build_gfp_standard (gfp, tdt))
    ;
  else if (_fault_tolerant)
    return _fill_output_buffer_failure();
  else
  {
    IWString msg;
    msg << "GFP_Server::tanimoto:cannot build fingerprint from '" << s << "'\n";
    _error_exit(msg);
  }

  if (NULL != _smiles)
    _insert_smiles_into_output_buffer(tdt);

  if (1 == nbrs)
    return _tanimoto_single_nbr(gfp);

  if (dist > 0.0f)
    return _tanimoto_within_distance(gfp, dist);

  if (nbrs > _pool_size)
    nbrs = _pool_size;

  return _tanimoto(gfp, nbrs);
}

static inline void
array_to_bits(const int * bfrom, 
              int n,
              unsigned char * destination)
{
  int nb = n / IW_BITS_PER_BYTE;

  for (auto i = 0; i < nb; ++i)
  {
    for (auto j = 0; j < IW_BITS_PER_BYTE; ++j)
    {
      if (*bfrom)
        destination[i] |= one_bit_8[j];

      bfrom++;
    }
  }

  return;
}

int
GFP_Server:: _insert_smiles_into_output_buffer (const IW_TDT & tdt)
{
  if (! tdt.dataitem_value(smiles_tag, _output_buffer))
    return 0;

  IWString pcn;
  if (tdt.dataitem_value(identifier_tag, pcn))
    _output_buffer << ' ' << pcn << '\n';

  return 1;
}

int
GFP_Server::tanimoto_from_smiles (IWString & s)
{
  searches_done++;

  int nbrs = -1;
  float dist = -1.0f;

  if (_extract_number_nbrs_smiles (s, nbrs, dist))
    ;
  else if (_fault_tolerant)
    return _fill_output_buffer_failure ();
  else
  {
    IWString msg;
    msg << "GFP_Server::tanimoto:cannot discern number of nbrs in '" << s << "'\n";
    _error_exit(msg);
  }

  _insert_smiles_into_output_buffer(s);

  Molecule m;
  if (! m.build_from_smiles(s))
  {
    IWString msg;
    msg << "GFP_Server::tanimoto:cannot build molecule from '" << s << "'";
    _error_exit(msg);
  }

  preprocess(m);

  GFP_Standard gfp;

  _mpr(m, gfp.molecular_properties());

  IWMFingerprint iwfp;
  iwfp.construct_fingerprint(m);     // _iwfp is a bit vector
  gfp.build_iwfp(iwfp.bits(), iwfp.nset());

  _mk(m, _tmp);
  gfp.build_mk(_tmp);
  _mk.set_level_2_fingerprint(_tmp);
  gfp.build_mk2(_tmp);

  if (1 == nbrs)
    return _tanimoto_single_nbr(gfp);

  if (dist > 0.0f)
    return _tanimoto_within_distance(gfp, dist);

  if (nbrs > _pool_size)
    nbrs = _pool_size;

  return _tanimoto(gfp, nbrs);
}

int
GFP_Server::_insert_smiles_into_output_buffer (const IWString & s)
{
  if (2 == s.nwords())
  {
    _output_buffer << s << '\n';
    return 1;
  }

  const_IWSubstring token;
  int i = 0;
  s.nextword(token, i);
  _output_buffer << token;
  s.nextword(token, i);
  _output_buffer << ' ' << token << '\n';

  return 1;
}

static Compare_DID compare_did;

int
GFP_Server::_tanimoto (const GFP_Standard & sfp, int nbrs)
{
  Distance_ID * did = new Distance_ID[nbrs]; std::unique_ptr<Distance_ID[]> free_did(did);
  
  for (auto i = 0; i < nbrs; ++i)
  {
    did[i]._d = sfp.tanimoto(_pool[i]);
    did[i]._ndx = i;

//  cerr << "Initial seeding i = " << i << " sim " << did[i]._d << endl;
  }

  std::priority_queue<Distance_ID, std::vector<Distance_ID>, Compare_DID> pq(did, did + nbrs);

  Distance_ID tmp;

  int n = 0;
  for (auto i = nbrs; i < _pool_size; ++i)
  {
    float s = sfp.tanimoto(_pool[i]);

    const Distance_ID & t = pq.top();

//  cerr << "Top of queue " << t._d << ", compare " << s << endl;

    if (s < t._d)
      continue;

    pq.pop();
    tmp._d = s;
    tmp._ndx = i;
    pq.push(tmp);
  }

  for (auto ndx = nbrs - 1; pq.size() > 0; ndx--)
  {
    const Distance_ID & d = pq.top();
    _did[ndx]._d = 1.0f - d._d;
    _did[ndx]._ndx = d._ndx;
    pq.pop();
  }

  return _fill_output_buffer(nbrs);
}

static int
did_comparatorq (const Distance_ID * did1,
                const Distance_ID * did2)
{
  if (did1->_d < did2->_d)
    return -1;
  if (did1->_d > did2->_d)
    return 1;

  return 0;
}

class Did_Comparator
{
  private:
  public:
    int operator ()(const Distance_ID &, const Distance_ID &) const;
};

int
Did_Comparator::operator () (const Distance_ID & did1,
                const Distance_ID & did2) const
{
  if (did1._d < did2._d)
    return -1;
  if (did1._d > did2._d)
    return 1;

  return 0;
}

static Did_Comparator did_comparator;

int 
GFP_Server::_tanimoto_single_nbr (const GFP_Standard & gfp)
{
  float rc = _pool[0].tanimoto(gfp);
  int idmax = 0;

  for (auto i = 1; i < _pool_size; ++i)
  {
    float s = _pool[i].tanimoto(gfp);
    if (s > rc)
    {
      rc = s;
      idmax = i;
    }
  }

  _did[0]._d = 1.0f - rc;
  _did[0]._ndx = idmax;

  return _fill_output_buffer(1);
}

int
GFP_Server::_tanimoto_within_distance (const GFP_Standard & gfp,
                                       float cutoff)
{
  int did_ndx = 0;

  cutoff = 1.0f - cutoff;     // convert to similarity

  for (auto i = 0; i < _pool_size; ++i)
  {
    float s = _pool[i].tanimoto(gfp);
    if (s >= cutoff)
    {
      _did[did_ndx]._d = 1.0f - _pool[i].tanimoto(gfp);
      _did[did_ndx]._ndx = i;
      did_ndx++;
    }
  }

  if (0 == did_ndx)
    return _fill_output_buffer(0);

  if (1 == did_ndx)
    return _fill_output_buffer(1);

  iwqsort(_did, did_ndx, did_comparator);

  return _fill_output_buffer(did_ndx);
}

/*
  zeromq thinks it owns the buffer, but needs a free function. In fact we own it, so the free function does nothing
*/

static void
free_fn (void * v,
         void * hint)
{
  return;
}

void
GFP_Server::doit ()
{
  zmq::context_t context(1);

  IWString tmp;
  tmp << "tcp://*:" << _input_port;

// Socket to receive messages on
  zmq::socket_t receiver(context, ZMQ_REP);
  receiver.bind(tmp.null_terminated_chars());

  IWString msg;
  msg << "Listening on '" << tmp << "'\n";
  _log_message(msg);

  msg.resize_keep_storage(0);
  msg << "Begin processing look " << _pool_size;
  _log_message(msg);

// Process tasks forever
  while (1)
  {
    _output_buffer.resize_keep_storage(0);

    zmq::message_t message;

    bool rc = receiver.recv(&message);

#ifdef DEBUG_SERVER
    cerr << "Received message of size " << message.size() << ", rc " << rc << endl;
#endif

    IWString iss(static_cast<char*>(message.data()), message.size());

//  Do the work

    int doing_shutdown = 0;

    if (iss.starts_with(smiles_tag))
      tanimoto_from_tdt (iss);
    else if (iss.starts_with(RELOAD_POOL))
    {
      if (! _rebuild_pool(iss))    // because the first thing we do is delete the arrays
        break;
    }
    else if (TERMINATE_PROCESSING == iss)
    {
      _do_shutdown();
      doing_shutdown = 1;
    }
    else
      tanimoto_from_smiles (iss);

//  Send results to sink
    message.rebuild(_output_buffer.rawdata(), static_cast<size_t>(_output_buffer.length()), free_fn, NULL);
    receiver.send(message);

    if (doing_shutdown)
      break;

//  Simple progress indicator for the viewer
//  std::cout << "." << std::flush;

    _searches_done++;
  }

  return;
}

static int
gfp_server_main (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vfL:p:s:i");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  set_global_aromaticity_type (Daylight);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  GFP_Server server;

  if (! server.parse_command_line(cl))
  {
    cerr << "Cannot initialise command line arguments\n";
    usage(2);
  }

  server.doit();

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_server_main(argc, argv);

  return rc;
}
