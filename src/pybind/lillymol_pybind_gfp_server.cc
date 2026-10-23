#include <memory>
#include <stdexcept>
#include <string>

#include "pybind11/pybind11.h"

#include "Utilities/GFP_Tools/gfp_server_lib.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_gfp_server, m) {
  py::class_<gfp_server::GfpNearNeighbourServer>(m, "GfpNearNeighbourServer")
    .def(py::init<>())
    .def(py::init([](const std::string& fname, int pool_size_hint, bool store_smiles) {
        auto result = std::make_unique<gfp_server::GfpNearNeighbourServer>();
        if (! result->Build(fname.c_str(), pool_size_hint, store_smiles)) {
          throw std::runtime_error("Cannot initialise GfpNearNeighbourServer from '" + fname + "'");
        }
        return result;
      }),
      py::arg("fname"),
      py::arg("pool_size_hint") = 0,
      py::arg("store_smiles") = true,
      "Load a GFP pool into memory"
    )
    .def("build",
      [](gfp_server::GfpNearNeighbourServer& server,
         const std::string& fname,
         int pool_size_hint,
         bool store_smiles)->bool {
        return server.Build(fname.c_str(), pool_size_hint, store_smiles);
      },
      py::arg("fname"),
      py::arg("pool_size_hint") = 0,
      py::arg("store_smiles") = true,
      "Load or replace the in-memory GFP pool"
    )
    .def("pool_size", &gfp_server::GfpNearNeighbourServer::pool_size,
      "Return the number of fingerprints loaded in the search pool"
    )
    .def("search_proto",
      [](gfp_server::GfpNearNeighbourServer& server, py::bytes serialized_request)->py::bytes {
        const std::string request = serialized_request;
        std::string reply;
        {
          py::gil_scoped_release release;
          reply = server.SearchSerialized(request);
        }
        return py::bytes(reply);
      },
      py::arg("request"),
      "Search one serialized gfp_server.NnRequest and return serialized gfp_server.NnReply"
    )
    .def("search_batch_proto",
      [](gfp_server::GfpNearNeighbourServer& server, py::bytes serialized_request)->py::bytes {
        const std::string request = serialized_request;
        std::string reply;
        {
          py::gil_scoped_release release;
          reply = server.SearchBatchSerialized(request);
        }
        return py::bytes(reply);
      },
      py::arg("request"),
      "Search one serialized gfp_server.NnBatchRequest and return serialized gfp_server.NnBatchReply"
    )
  ;
}
