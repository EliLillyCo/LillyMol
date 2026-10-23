#include <stdexcept>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "Utilities/GFP_Tools/gfp_server_lib.h"

namespace nb = nanobind;

NB_MODULE(lillymol_gfp_server, m) {
  nb::class_<gfp_server::GfpNearNeighbourServer>(m, "GfpNearNeighbourServer")
      .def(nb::init<>())
      .def("__init__",
           [](gfp_server::GfpNearNeighbourServer* server, const std::string& fname) {
             new (server) gfp_server::GfpNearNeighbourServer();
             if (!server->Build(fname.c_str(), 0, true)) {
               throw std::runtime_error(
                   "Cannot initialise GfpNearNeighbourServer from '" + fname + "'");
             }
           },
           nb::arg("fname"),
           "Load a GFP pool into memory")
      .def("__init__",
           [](gfp_server::GfpNearNeighbourServer* server, const std::string& fname,
              int pool_size_hint, bool store_smiles) {
             new (server) gfp_server::GfpNearNeighbourServer();
             if (!server->Build(fname.c_str(), pool_size_hint, store_smiles)) {
               throw std::runtime_error(
                   "Cannot initialise GfpNearNeighbourServer from '" + fname + "'");
             }
           },
           nb::arg("fname"), nb::arg("pool_size_hint"), nb::arg("store_smiles"),
           "Load a GFP pool into memory")
      .def("build",
           [](gfp_server::GfpNearNeighbourServer& server, const std::string& fname,
              int pool_size_hint, bool store_smiles) {
             return server.Build(fname.c_str(), pool_size_hint, store_smiles);
           },
           nb::arg("fname"), nb::arg("pool_size_hint") = 0,
           nb::arg("store_smiles") = true,
           "Load or replace the in-memory GFP pool")
      .def("pool_size", &gfp_server::GfpNearNeighbourServer::pool_size,
           "Return the number of fingerprints loaded in the search pool")
      .def("search_proto",
           [](gfp_server::GfpNearNeighbourServer& server, nb::bytes serialized_request) {
             const std::string request(serialized_request.c_str(), serialized_request.size());
             std::string reply;
             {
               nb::gil_scoped_release release;
               reply = server.SearchSerialized(request);
             }
             return nb::bytes(reply.data(), reply.size());
           },
           nb::arg("request"),
           "Search one serialized gfp_server.NnRequest and return serialized gfp_server.NnReply")
      .def("search_batch_proto",
           [](gfp_server::GfpNearNeighbourServer& server, nb::bytes serialized_request) {
             const std::string request(serialized_request.c_str(), serialized_request.size());
             std::string reply;
             {
               nb::gil_scoped_release release;
               reply = server.SearchBatchSerialized(request);
             }
             return nb::bytes(reply.data(), reply.size());
           },
           nb::arg("request"),
           "Search one serialized gfp_server.NnBatchRequest and return serialized gfp_server.NnBatchReply");
}
