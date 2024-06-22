/*
  A server for doing gfp searches
  Comminucates with clients via ZeroMQ with 
  data transferred via serialized protocol buffers.
*/

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <queue>
#include <string>

#include "absl/strings/string_view.h"

#include "zmq.hpp"

using std::cerr;

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#define IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/maccskeys_fn5.h"
#include "Molecule_Tools/mpr.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#include "Utilities/GFP_Tools/nn_request.pb.h"

const char* prog_name = nullptr;

static int verbose = 0;

// static std::unique_ptr<re2::RE2> fingerprint_rx("^\\|");

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static void
usage(int rc) {
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Server for molecular similarity\n";
  cerr << " -p <port>      port for both input (smiles or gfp) and results\n";
  cerr << " -i             also read the smiles\n";
  cerr << " -s <int>       size of the input fingerprint file\n";
  cerr << " -f             fault tolerance mode\n";
  cerr << " -L <fname>     file name for logging\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*
  We need a lightweight structure for holding distance/id information.
*/

struct Distance_ID {
  float _d;
  int _ndx;  // index into pool
};

class Compare_DID {
 private:
 public:
  bool
  operator()(const Distance_ID& lhs, const Distance_ID& rhs) const {
    return lhs._d > rhs._d;
  }
};

class GFP_Server {
 private:
  // The gfp file from which we read fingerprints.
  IWString _fname;

  int _pool_size;
  GFP_Standard* _pool;

  //  We need a structure for accumulating results

  Distance_ID* _did;

  // Some measure of control over how errors are handled.
  int _fault_tolerant;
  std::ofstream _stream_for_logging;
  bool _stream_for_logging_open;

  // The port number over which we communicate
  uint32_t _input_port;

  Molecular_Properties_Generator _mpr;
  MACCSKeys _mk;

  int _tmp[2048];  // used by maccskeys and iwfp

  //  the names of the fingerprints in the pool

  std::string* _id;

  int _read_smiles;

  std::string* _smiles;

  unsigned int _searches_done;

  Chemical_Standardisation _chemical_standardisation;

  //  private functions

  void _free_arrays();

  void _log_message(const IWString&);
  void _error_exit(const IWString&);

  int FillResultProto(int nbrs, nnbr::NearNeighbours& result);
  int _tanimoto_single_nbr(const GFP_Standard& gfp, nnbr::NearNeighbours& result);
  int _tanimoto_within_distance(const GFP_Standard& gfp, float cutoff, nnbr::NearNeighbours& result);
  int _tanimoto(const GFP_Standard& sfp, int nbrs);
  int _tanimoto(const GFP_Standard& sfp, int nbrs, nnbr::NearNeighbours& result);
  int ReloadPool();
  int _do_shutdown();
  void Preprocess(Molecule& m);
  int DoNearNeighbours(const gfp_server::NnRequest& req,
                             gfp_server::Reply &reply);

  int DoServerRequest(const gfp_server::ServerRequest& req, int& doing_shutdown,
                             gfp_server::Reply &reply);
  int SmilesToGfp(const std::string& smiles, GFP_Standard& gfp);

 public:
  GFP_Server();
  ~GFP_Server();

  int parse_command_line(Command_Line& cl);

  int build(const char* fname);
  int build(iwstring_data_source&);

  int tanimoto_from_tdt(IWString&);
  int tanimoto_from_smiles(IWString&);

  void doit();
};

GFP_Server::GFP_Server() {
  _pool_size = 0;
  _pool = nullptr;
  _fault_tolerant = 0;

  _did = new Distance_ID[1];  // the most common case

  _input_port = 0;

  _stream_for_logging_open = false;

  _id = nullptr;
  _read_smiles = 0;
  _smiles = nullptr;

  _searches_done = 0;

  _chemical_standardisation.activate_all();

  return;
}

void
GFP_Server::_free_arrays() {
  if (NULL != _pool) {
    delete[] _pool;
    _pool = nullptr;
  }

  _pool_size = 0;

  if (NULL != _did) {
    delete[] _did;
    _did = nullptr;
  }

  if (NULL != _id) {
    delete[] _id;
    _id = nullptr;
  }

  if (NULL != _smiles) {
    delete[] _smiles;
    _smiles = nullptr;
  }

  return;
}

GFP_Server::~GFP_Server() {
  _free_arrays();

  return;
}

void
GFP_Server::Preprocess(Molecule& m) {
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  _chemical_standardisation.process(m);

  return;
}

void
GFP_Server::_log_message(const IWString& msg) {
  if (!_stream_for_logging_open) {
    return;
  }

  _stream_for_logging << msg;

  if (!msg.ends_with('\n')) {
    _stream_for_logging << '\n';
  }

  _stream_for_logging.flush();

  return;
}

void
GFP_Server::_error_exit(const IWString& msg) {
  _log_message(msg);

  _do_shutdown();

  exit(1);
}

int
GFP_Server::parse_command_line(Command_Line& cl) {
  if (!cl.option_present('p')) {
    cerr << "Must specify port to use via the -p option\n";
    return 0;
  } else if (!cl.value('p', _input_port) || _input_port < 30) {
    cerr << "Invalid input port specification (-i)\n";
    return 0;
  } else if (verbose) {
    IWString msg;
    msg << "Input from port " << _input_port;
    cerr << msg << "\n";
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', _pool_size) || _pool_size < 1) {
      cerr << "The pool size must be a whole +ve number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Fingerprint pool explicitly sized to " << _pool_size << " fingerprints\n";
    }
  }

  if (cl.option_present('f')) {
    _fault_tolerant = 1;

    if (verbose) {
      cerr << "Will run in fault tolerant mode\n";
    }
  }

  if (cl.option_present('L')) {
    const char* l = cl.option_value('L');

    _stream_for_logging.open(l, std::fstream::out);

    if (!_stream_for_logging.good()) {
      cerr << "GFP_Server::parse_command_line:cannot initialise logging stream '" << l
           << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Logging messages written to '" << l << "'\n";
    }

    _stream_for_logging_open = true;
  }

  if (cl.option_present('i')) {
    _read_smiles = 1;

    if (verbose) {
      cerr << "Will include smiles in memory\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "GFP_Server::parse_command_line:no fingerprint file specified\n";
    return 0;
  }

  return build(cl[0]);
}

int
GFP_Server::build(const char* fname) {
  _fname = fname;

  iwstring_data_source input(fname);

  if (!input.good()) {
    IWString msg;
    msg << "GFP_Server::build:cannot open '" << fname << "'";
    _error_exit(msg);
  }

  return build(input);
}

static int
build_gfp_standard(GFP_Standard& sfp, IW_General_Fingerprint& gfp) {
  sfp.build_molecular_properties(gfp.molecular_properties_integer());
  sfp.build_iw(gfp[0]);
  sfp.build_mk(gfp[1]);
  sfp.build_mk2(gfp[2]);

  return 1;
}

static int
build_gfp_standard(GFP_Standard& sfp, IWString& id, IW_TDT& tdt) {
  int fatal;
  IW_General_Fingerprint gfp;

  if (!gfp.construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  id = gfp.id();

  return build_gfp_standard(sfp, gfp);
}

#ifdef NOP_LONGER_USED_Q
static int
build_gfp_standard(GFP_Standard& sfp, IW_TDT& tdt) {
  int fatal;
  IW_General_Fingerprint gfp;

  if (!gfp.construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  return build_gfp_standard(sfp, gfp);
}
#endif

int
GFP_Server::build(iwstring_data_source& input) {
  if (0 == _pool_size) {
    _pool_size = input.count_records_starting_with("|");
    //  cerr << "Pool size " << _pool_size << '\n';
  }

  if (0 == _pool_size) {
    IWString msg;
    msg << "GFP_Server::build:no fingerprints in file\n";
    _error_exit(msg);
  }

  if (verbose) {
    IWString msg;
    msg << "GFP_Server::build:fingerprint file contains " << _pool_size
        << " fingerprints";
    _log_message(msg);
  }

  _pool = new GFP_Standard[_pool_size];
  _did = new Distance_ID[_pool_size];
  _id = new std::string[_pool_size];

  if (_read_smiles) {
    _smiles = new std::string[_pool_size];
  }

  if (NULL == _pool || NULL == _did || NULL == _id) {
    IWString msg;
    msg << "GFP_Server::build:cannot allocate " << _pool_size << " fingerprints";
    _error_exit(msg);
  }

  IWString msg;
  msg << "Begin building pool size " << _pool_size;
  _log_message(msg);

  for (auto i = 0; i < _pool_size; ++i) {
    IW_TDT tdt;
    if (!tdt.next(input)) {
      IWString msg;
      msg << "GFP_Server::build:premature, expected " << _pool_size << " got " << i;
      _log_message(msg);
      _pool_size = i - 1;
      break;
    }

    IWString tmp;
    if (!build_gfp_standard(_pool[i], tmp, tdt)) {
      IWString msg;
      msg << "GFP_Server::build:fatal error reading fingerprint " << i;
      _error_exit(msg);
    }
    _id[i] = tmp.AsString();

    if (NULL != _smiles) {
      IWString tmp;
      tdt.dataitem_value(smiles_tag, tmp);
      _smiles[i] = tmp.AsString();
    }
  }

  return _pool_size;
}

/*
  We ignore errors where we can just leave the existing pool in place
*/

int
GFP_Server::ReloadPool() {

  iwstring_data_source input(_fname.null_terminated_chars());

  if (!input.good()) {
    IWString msg;
    msg << "GFP_Server::_rebuild_pool:cannot open '" << _fname << "'\n";
    _error_exit(msg);
    return 0;
  }

  int nps = input.count_records_starting_with("|");
  if (nps < (9 * _pool_size / 10))  // pool cannot shrink by more than 10%
  {
    IWString msg;
    msg << "GFP_Server::_rebuild_pool:new fingerprint file too small '" << _fname
        << "' count " << nps;
    _error_exit(msg);
    return 0;
  }

  _free_arrays();

  _pool_size = nps;

  if (!build(input)) {
    IWString msg("Cannot rebuilt build");
    _error_exit(msg);
  }

  _log_message("Pool reloaded");
  return 1;
}

int
GFP_Server::_do_shutdown() {

  IWString msg;
  msg << "Shutdown, performed " << _searches_done << " searches\n";
  _log_message(msg);

  return 1;
}

int
GFP_Server::FillResultProto(int nbrs, nnbr::NearNeighbours& result) {
#ifdef DEBUG_SERVER
  cerr << "Server filling output buffer for " << nbrs << " nbrs\n";
#endif

  result.mutable_nbr()->Reserve(nbrs);

  for (auto i = 0; i < nbrs; ++i) {
    const Distance_ID& didi = _did[i];

    nnbr::Nbr* n = result.add_nbr();
    n->set_id(_id[didi._ndx]);
    n->set_dist(didi._d);
    if (NULL != _smiles) {
      n->set_smi(_smiles[didi._ndx]);
    }
  }

  return 1;
}

static inline void
array_to_bits(const int* bfrom, int n, unsigned char* destination) {
  int nb = n / IW_BITS_PER_BYTE;

  for (auto i = 0; i < nb; ++i) {
    for (auto j = 0; j < IW_BITS_PER_BYTE; ++j) {
      if (*bfrom) {
        destination[i] |= one_bit_8[j];
      }

      bfrom++;
    }
  }

  return;
}

int
GFP_Server::_tanimoto(const GFP_Standard& sfp, int nbrs, nnbr::NearNeighbours& result) {
  Distance_ID* did = new Distance_ID[nbrs];
  std::unique_ptr<Distance_ID[]> free_did(did);

  for (auto i = 0; i < nbrs; ++i) {
    did[i]._d = sfp.tanimoto(_pool[i]);
    did[i]._ndx = i;

    //  cerr << "Initial seeding i = " << i << " sim " << did[i]._d << '\n';
  }

  std::priority_queue<Distance_ID, std::vector<Distance_ID>, Compare_DID> pq(did,
                                                                             did + nbrs);

  Distance_ID tmp;

  for (auto i = nbrs; i < _pool_size; ++i) {
    float s = sfp.tanimoto(_pool[i]);

    const Distance_ID& t = pq.top();

    //  cerr << "Top of queue " << t._d << ", compare " << s << '\n';

    if (s < t._d) {
      continue;
    }

    pq.pop();
    tmp._d = s;
    tmp._ndx = i;
    pq.push(tmp);
  }

  for (auto ndx = nbrs - 1; pq.size() > 0; ndx--) {
    const Distance_ID& d = pq.top();
    _did[ndx]._d = 1.0f - d._d;
    _did[ndx]._ndx = d._ndx;
    pq.pop();
  }

  result.mutable_nbr()->Reserve(nbrs);
  for (int i = 0; i < nbrs; ++i) {
  }

  return FillResultProto(nbrs, result);
}

#ifdef NOP_LONGER_USED_Q
static int
did_comparatorq(const Distance_ID* did1, const Distance_ID* did2) {
  if (did1->_d < did2->_d) {
    return -1;
  }
  if (did1->_d > did2->_d) {
    return 1;
  }

  return 0;
}
#endif

class Did_Comparator {
 private:
 public:
  int operator()(const Distance_ID&, const Distance_ID&) const;
};

int
Did_Comparator::operator()(const Distance_ID& did1, const Distance_ID& did2) const {
  if (did1._d < did2._d) {
    return -1;
  }
  if (did1._d > did2._d) {
    return 1;
  }

  return 0;
}

static Did_Comparator did_comparator;

int
GFP_Server::_tanimoto_single_nbr(const GFP_Standard& gfp, nnbr::NearNeighbours& result) {
  float rc = _pool[0].tanimoto(gfp);
  int idmax = 0;

  for (auto i = 1; i < _pool_size; ++i) {
    float s = _pool[i].tanimoto(gfp);
    if (s > rc) {
      rc = s;
      idmax = i;
    }
  }

  _did[0]._d = 1.0f - rc;
  _did[0]._ndx = idmax;

  return FillResultProto(1, result);
}

int
GFP_Server::_tanimoto_within_distance(const GFP_Standard& gfp, float cutoff, nnbr::NearNeighbours& result) {
  int did_ndx = 0;

  cutoff = 1.0f - cutoff;  // convert to similarity

  for (auto i = 0; i < _pool_size; ++i) {
    float s = _pool[i].tanimoto(gfp);
    if (s >= cutoff) {
      _did[did_ndx]._d = 1.0f - _pool[i].tanimoto(gfp);
      _did[did_ndx]._ndx = i;
      did_ndx++;
    }
  }

  if (0 == did_ndx) {
    return FillResultProto(0, result);
  }

  if (1 == did_ndx) {
    return FillResultProto(1, result);
  }

  iwqsort(_did, did_ndx, did_comparator);

  return FillResultProto(did_ndx, result);
}

/*
  zeromq thinks it owns the buffer, but needs a free function. In fact we own it, so the
  free function does nothing
*/

static void
free_fn(void* v, void* hint) {
  return;
}

void
GFP_Server::doit() {
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
  msg << "Begin processing loop " << _pool_size << " fingerprints\n";
  _log_message(msg);

  std::string reply_buffer;

  // Process tasks forever
  while (1) {
    zmq::message_t message;

    std::optional<long unsigned int> rc = receiver.recv(message, zmq::recv_flags::none);

    if (! rc) {
      _log_message("recv failed");
      continue;
    }
#ifdef DEBUG_SERVER
    cerr << "Received message of size " << message.size() << ", rc " << rc << '\n';
#endif

    absl::string_view sview((const char*) message.data(), message.size());
    gfp_server::Request req;
    if (! req.ParseFromString(sview)) {
      msg.resize_keep_storage(0);
      msg << "Cannot decode " << message.size() << " bytes";
      _log_message(msg);
      continue;
    }

    int doing_shutdown = 0;

    gfp_server::Reply reply;
    if (req.has_server_request()) {
      DoServerRequest(req.server_request(), doing_shutdown, reply);
    } else {
      ++_searches_done;
      DoNearNeighbours(req.nn_request(), reply);
    }

    //  Return results to client
    reply.SerializeToString(&reply_buffer);
    message.rebuild(reply_buffer.data(),
                    reply_buffer.size(), free_fn, NULL);
    receiver.send(message, zmq::send_flags::none);

    if (doing_shutdown) {
      break;
    }
  }

  return;
}

int
GFP_Server::DoServerRequest(const gfp_server::ServerRequest& req, int& doing_shutdown,
                             gfp_server::Reply &reply) {
  switch (req.request()) {
    case gfp_server::ServerRequest::NONE:
      reply.set_status(gfp_server::Status::OK);
      return 1;
    case gfp_server::ServerRequest::SHUTDOWN:
      doing_shutdown = 1;
      reply.set_status(gfp_server::Status::OK);
      return 1;
    case gfp_server::ServerRequest::RELOAD:
      if (ReloadPool()) {
        reply.set_status(gfp_server::Status::OK);
        return 1;
      }
      reply.set_status(gfp_server::Status::RELOAD_FAILED);
      return 0;
    default:
      reply.set_status(gfp_server::Status::NO_DIRECTIVE);
      return 0;
  }

  return 1;
}

int
GFP_Server::DoNearNeighbours(const gfp_server::NnRequest& req,
                             gfp_server::Reply &reply) {
  if (req.smiles().empty()) {
    reply.set_status(gfp_server::Status::NO_SMILES);
    return 0;
  }

  GFP_Standard gfp;
  if (! SmilesToGfp(req.smiles(), gfp)) {
    reply.set_status(gfp_server::Status::BAD_SMILES);
    return 0;
  }

  if (req.has_nbrs()) {
  } else if (req.has_distance()) {
  } else {
    reply.set_status(gfp_server::Status::NO_DIRECTIVE);
    return 0;
  }

  if (req.has_id()) {
    reply.mutable_result()->set_name(req.id());
  }

#ifdef DEBUG_SERVER
  cerr << "looking for " << req.ShortDebugString() << '\n';
#endif

  if (req.nbrs() == 1) {
    _tanimoto_single_nbr(gfp, *reply.mutable_result());
    reply.set_status(gfp_server::Status::OK);
    return 1;
  }

  if (req.has_distance() && req.distance() > 0.0f) {
    _tanimoto_within_distance(gfp, req.distance(), *reply.mutable_result());
    reply.set_status(gfp_server::Status::OK);
    return 1;
  }

  int nbrs = req.nbrs();

  if (nbrs > _pool_size) {
    nbrs = _pool_size;
  }

  _tanimoto(gfp, nbrs, *reply.mutable_result());

  reply.set_status(gfp_server::Status::OK);
  return 1;
}

int
GFP_Server::SmilesToGfp(const std::string& smiles, GFP_Standard& gfp) {
  Molecule m;
  if (! m.build_from_smiles(smiles)) {
    // TODO
    _log_message("");
    return 0;
  }

  Preprocess(m);

  IWMFingerprint iwfp;
  iwfp.construct_fingerprint(m);
  gfp.build_iwfp(iwfp.bits(), iwfp.nset());

  std::fill_n(_tmp, 2048, 0);

  _mk(m, _tmp);
  gfp.build_mk(_tmp, _mk.nbits());

  _mk.set_level_2_fingerprint(_tmp);
  gfp.build_mk2(_tmp, _mk.nbits());

  _mpr(m, gfp.molecular_properties());

  return 1;
}

static int
gfp_server_main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vfL:p:s:i");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  set_global_aromaticity_type(Daylight);

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  GFP_Server server;

  if (!server.parse_command_line(cl)) {
    cerr << "Cannot initialise command line arguments\n";
    usage(2);
  }

  server.doit();

  if (verbose) {
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  GOOGLE_PROTOBUF_VERIFY_VERSION;

  int rc = gfp_server_main(argc, argv);

  return rc;
}
