# Client to talk to a host running gfp_server.cc

from absl import app
from absl import flags
from absl import logging

import zmq

import Utilities.GFP_Tools.nn_request_pb2 as nn_request

from google.protobuf.json_format import MessageToJson

FLAGS = flags.FLAGS

flags.DEFINE_string("port", "", "Port tcp://host:port")
flags.DEFINE_integer("niter", 1, "Number of repeat lookups (for benchmarking) def 1")
flags.DEFINE_string("smiles", "C1=NC(=CC2=C1C=CC=C2)NC(=O)NCC", "Smiles to use")
flags.DEFINE_string("id", "CHEMBL4076111", "Name of the smiles")
flags.DEFINE_integer("nbrs", 1, "Number of neighbours to find")
flags.DEFINE_boolean("json", False, "Write output in JSON form")
flags.DEFINE_boolean("shutdown", False, "Shut down the server")

def gfp_client(argv):
  context = zmq.Context()
  socket = context.socket(zmq.REQ)
  if len(FLAGS.port) == 0:
    logging.error("Must specify port via the -p option")
    return 1
  socket.connect(FLAGS.port)

  if FLAGS.shutdown:
    req = nn_request.Request()
    req.server_request.request = nn_request.ServerRequest.SHUTDOWN
    serialised = req.SerializeToString()
    socket.send(serialised)
    return

  for i in range(FLAGS.niter):
    req = nn_request.Request()
    req.nn_request.smiles = FLAGS.smiles
    req.nn_request.id = FLAGS.id
    req.nn_request.nbrs = FLAGS.nbrs
    serialised = req.SerializeToString()
    socket.send(serialised)

    message = socket.recv()
    proto = nn_request.Reply()
    proto.ParseFromString(message);
    if FLAGS.json:
      print(MessageToJson(proto))
    else:
      print(proto)


if __name__ == "__main__":
  app.run(gfp_client)
