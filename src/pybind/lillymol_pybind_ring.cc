#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "Molecule_Lib/path.h"

namespace py = pybind11;

static std::string
RingToString(const Ring& r) {
  IWString result;
  result << "<Ring size " << r.size();
  if (r.is_aromatic()) {
    result << " arom";
  } else if (r.is_non_aromatic()) {
    result << " aliph";
  }
  result << " frag " << r.fragment_membership();
  result << " fsid " << r.fused_system_identifier();
  result << " fused " << r.is_fused();
  result << " atoms";
  for (int i = 0; i < r.number_elements(); ++i) {
    result << ' ' << r[i];
  }

  result << '>';

  return std::string(result.data(), result.size());
}

PYBIND11_MODULE(lillymol_ring, r)
{
  py::class_<Ring, Set_of_Atoms>(r, "Ring")
    .def(py::init<>())
    .def("size", &Ring::size, "size")
    .def("ring_number", &Ring::ring_number, "Unique ring number")
    .def("fragment_membership", &Ring::fragment_membership, "Fragment membership")
    .def("fused_system_identifier", &Ring::fused_system_identifier, "Fused system identifier")
    .def("is_fused",
      [](const Ring& r)->bool{
        return r.is_fused();
      },
      "True if the ring is fused"
    )
    .def("is_fused_to",
      [](const Ring& r1, const Ring* r2)->bool{
        return r1.is_fused_to(r2);
      },
      "True if rings are fused"
    )
    .def("fused_ring_neighbours", &Ring::fused_ring_neighbours, "Number of fused ring neighbours")
    .def("largest_number_of_bonds_shared_with_another_ring", &Ring::largest_number_of_bonds_shared_with_another_ring, "largest_number_of_bonds_shared_with_another_ring")
    .def("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours, "strongly_fused_ring_neighbours")
    .def("contains_bond",
      [](const Ring& r, atom_number_t a1, atom_number_t a2)->bool{
        return r.contains_bond(a1, a2);
      },
      "True if a1 and a2 are adjacent"
    )
    .def("contains_both",
      [](const Ring& r, atom_number_t a1, atom_number_t a2)->bool{
        return r.contains_both(a1, a2);
      },
      "True if ring contains both a1 and a2"
    )
    .def("is_aromatic",
      [](const Ring& r)->bool{
        return r.is_aromatic();
      },
      "True if the ring is aromatic"
    )
    .def("__getitem__",
      [](const Ring& me, int ndx)->atom_number_t{
        return me.item(ndx);
      })
    .def("__repr__",
      [](const Ring &r) {
        return RingToString(r);
      }
    )
    .def("__str__",
      [](const Ring& r) {
        return RingToString(r);
      }
    )
    .def("__iter__",
      [](const Ring&r) {
        return py::make_iterator(r.begin(), r.end());
      },
      py::keep_alive<0, 1>()
    )
  ;
}
