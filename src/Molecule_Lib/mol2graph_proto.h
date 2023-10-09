#ifndef MOLECULE_TOOLS_MOL2GRAPH_PROTO_H_
#define MOLECULE_TOOLS_MOL2GRAPH_PROTO_H_


#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Lib/mol2graph.pb.h"

// interconversions between the two different means of controlling graph conversions.
// Build a Mol2Graph class from the proto representation.
extern Mol2Graph Mol2GraphFromProto(const LLYMol::Mol2Graph& mol2graph);
// Build a Mol2Graph proto from the class representation.
extern LLYMol::Mol2Graph Mol2GraphToProto(const Mol2Graph& mol2graph);

#endif // MOLECULE_TOOLS_MOL2GRAPH_PROTO_H_
