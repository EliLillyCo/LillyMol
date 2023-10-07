#ifndef MOLECULE_TOOLS_RXN_TO_OPENRXN_H
#define MOLECULE_TOOLS_RXN_TO_OPENRXN_H

#include <memory>

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/rxn_file.h"

//#include "external/ord_schema/proto/reaction.pb.h"
//#include "external/ord_schema/proto/dataset.pb.h"
#include "proto/dataset.pb.h"
#include "proto/reaction.pb.h"

namespace rxn_to_openrxn {

struct JobOptions
{
  // By default, we die upon failing to interpret a reaction.
  bool ignore_bad_reactions = false;
  // If the reactions have already passed through out software.
  bool component_grouping_is_plus = false;
  // The key values assigned for the various reaction components.
  IWString reagent_stem = "REAGENT_";
  IWString agent_stem = "AGENT_";
  // Report processing periodically.
  int report = 0;
  // Write out the DebugString of each proto formed.
  bool debug_string = false;
};


// Given `job_options` as possible behaviour modifiers, convert `rxn` to
// an Open Reaction Database proto.
// `buffer` is a record containing the reaction smiles as the first token.
// Mostly job_options provides the name stems for the different
// kinds of reaction inputs.
ord::Reaction BuildOrdReaction(const JobOptions& job_options,
                   const const_IWSubstring& buffer,
                   RXN_File& rxn);

}  // namespace rxn_to_openrxn

#endif  // MOLECULE_TOOLS_RXN_TO_OPENRXN_H
