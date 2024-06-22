# gfp_server
A common need is to find the nearest neighbour(s) for a molecule with respect to one or
more collections.

At the command line this can be done via `gfp_lnearneighbours`, but this is not necessarily
fast. For example, to find the five nearest neighbours of a single molecule vs a recent Chembl,
2.249M molecules
```
gfp_make needle.smi > needle.gfp
gfp_lnearneighbours -p needle.gfp  -n 5 /path/to/chembl/chembl.gfp > needle.nn
```
Takes about 5 seconds on fairly old hardware, but with the database residing on an SSD -
the fingerprint file is 1.3GB in size.

But most of that time was spent reading the file of fingerprints, not doing similarity
calculations.

## gfp_server
This tool takes a simplified approach to the problem, by supporting only one GFP fingerprint,
the so-called "standard" (default) fingerprint, a composite fingerprint that consists of

* Molecular Properties
* iwfp linear fingerprint (2048 bits)
* mk fingerprint (192 bits)
* mk2 fingerprint (192 bits)

These can be computed within C++. There is no allowance for using other kinds of fingerprints.

That said, this set of fingerprints has been in use for many years and is widely accepted
as a reasonable approximation to human perception of molecular similarity.

The system works by launching a C++ server that slurps an existing fingerprint file
into RAM, and then waits for similarity requests on a port. The haystack is read
once, and then re-used as requested.
Communication with the
server is via a serialised Protocol Buffer, which means it is accessible from
many languages.

Performance seems good enough. Making 100 calls for any number of neighbours
up to 10 takes less than 10 seconds - 100 msec per call, 10 per second.

### Architecture
Socket communication is via ZeroMQ, a widely used inter process communcation
system, that is supported from a wide variety of languages. The other choice
might have been grpc. This post
[stack overflow](https://stackoverflow.com/questions/39350681/grpc-and-zeromq-comparsion)
provides a fairly detailed comparison. For this purpose ZermMQ seems like a reasonable
choice.

The client builds a `gfp_server.Request` proto
[nn_request.proto](/src/Utilities/GFP_Tools/nn_request.proto). The client
must then serialise that proto and send it to the server.

The server decodes the proto and does what is requested.

The server responds with a `gfp_server.Reply` serialised proto. This
will contain a status value, and if successful, a `nnbr::NearNeighbours`
message with the results.

Lanuching a server instance might look like
```
gfp_server -v -p 5555 -i -L /tmp/logfile /path/to/chembl/chembl.gfp > /tmp/out 2>&1
```
Any client will need to point to the port number chosen.
This can of course be done across hosts thanks to ZermMQ. At this
stage I have not really figured out logging, with some messgaes coming to the
`-L` file, and others appearing on stderr. Maybe it would be cleaner if it just
used stdout/stderr.

A python client for this server this might look like
```
  # Some boilerplate zmq startup
  context = zmq.Context()
  socket = context.socket(zmq.REQ)
  socket.connect("tcp://localhost:5555")

  # Enter some kind of processing loop...
    # Populate the proto, we want the 10 nearest neighbours.
    req = nn_request_pb2.Request()
    req.nn_request.smiles = "C1=NC(=CC2=C1C=CC=C2)NC(=O)NCC"
    req.nn_request.id = "CHEMBL4076111"
    req.nn_request.nbrs = 10

    # Serialise and send to the host.
    serialised = req.SerializeToString()
    socket.send(serialised)

    # Wait for a response
    message = socket.recv()

    # Decode the response and process it
    proto = nn_request_pb2.Reply()
    proto.ParseFromString(message);
    print(f"Received {proto}")
```

The example client is in the repo as
[gfp_client.py](src/Utilities/GFP_Tools/gfp_client_app.py)

ZermMQ for python can be installed via 
```
pip install pyzmq
```

### Server Control
The Request message sent to the server can be one of two forms

1. Perform a similarity search on the smiles provided.
2. Perform a server related action.

Currently two server actions are supported

. Shutdown the server
. Reload the server - from the same fingerprint file from which it was launched.

The reload functionality is designed for the case where a collection gets
rebuilt periodically. This has not been tested extensively, and it might
be easier to restart the server.

There is no security or authentication associated with the server, it is
assumed to be running in a secured network.

### Failure Modes
This has not been tested under stressful situations. I do not know what
would happen if a client sent a request, and then exited before
accepting the reply. Presumably this is a well known failure
mode for ZeroMQ.

A server could be overwhelmed by requests. ZeroMQ will queue requests
up to a point. Again, we are relying on default ZeroMQ behaviour.
