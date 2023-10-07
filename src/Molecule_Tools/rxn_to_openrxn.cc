#include "Molecule_Tools/rxn_to_openrxn.h"

namespace rxn_to_openrxn {

// Return 'stem << i' as std::string.
std::string
FormName(const IWString& stem, const int i)
{
  IWString tmp(stem);
  tmp << i;

  return std::string(tmp.data(), tmp.length());
}

// Just a cast.
template <typename T>
std::string
ToStdString(const T& s)
{
  return std::string(s.data(), s.length());
}

// Add a SMILES CompoundIdentifier to `component`.
void
AddSmilesToComponent(const IWString& smiles, ord::Compound* component)
{
  ord::CompoundIdentifier* identifier = component->add_identifiers();
  identifier->set_type(ord::CompoundIdentifier::SMILES);
  identifier->set_value(ToStdString(smiles));
}

// Return a new ReactionInput which has a component consisting
// of the structure `smiles` in the role of `role`.
ord::ReactionInput
MakeReactionInput(const IWString& smiles, const ord::Compound::ReactionRole::ReactionRoleType role)
{
  ord::ReactionInput reaction_input;    // To be returned.

  ord::Compound* component = reaction_input.add_components();
  component->set_reaction_role(role);

  AddSmilesToComponent(smiles, component);

  return reaction_input;
}

ord::Reaction
BuildOrdReaction(const JobOptions& job_options,
  const const_IWSubstring& buffer, RXN_File& rxn)
{
  ord::Reaction ord_proto;    // To be returned.

  ord::ReactionIdentifier* id = ord_proto.add_identifiers();
  id->set_type(ord::ReactionIdentifier::ATOM_MAPPED_SMILES);
  const_IWSubstring rxn_smiles;
  int i = 0;
  buffer.nextword(rxn_smiles, i);
  const_IWSubstring maybe_f;
  buffer.nextword(maybe_f, i);
  const_IWSubstring name;
  if (maybe_f.starts_with("|f:")) {
    IWString tmp;
    tmp << rxn_smiles << ' ' << maybe_f;
    id->set_value(ToStdString(tmp));
    buffer.nextword(name, i);
  } else {
    id->set_value(ToStdString(rxn_smiles));
    name = maybe_f;   // Was not f.
  }

  ord::ReactionProvenance* provenence = ord_proto.mutable_provenance();
  provenence->set_patent(ToStdString(name));

  resizable_array<float> yields;

  const_IWSubstring year;
  if (buffer.nextword(year, i)) {
    const_IWSubstring maybe_yield;
    while (buffer.nextword(maybe_yield, i)) {
      if (! maybe_yield.ends_with('%'))
        continue;

      maybe_yield.chop();
      float yield;
      if (! maybe_yield.numeric_value(yield) || yield < 0.0f || yield > 100.0f) {
        continue;
      }
      yields.add(yield);
    }
  }

  for (int i = 0; i < rxn.number_reagents(); ++i)
  {
    ISIS_RXN_FILE_Molecule& reagent = rxn.reagent(i);
    if (reagent.natoms() == 0)
    {
      continue;
    }

    const std::string key = FormName(job_options.reagent_stem, i);

    ord::ReactionInput reaction_input =
        MakeReactionInput(reagent.smiles(), ord::Compound::ReactionRole::REACTANT);

    (*ord_proto.mutable_inputs())[key] = std::move(reaction_input);
  }

  for (int i = 0; i < rxn.number_agents(); ++i)
  {
    ISIS_RXN_FILE_Molecule& agent = rxn.agent(i);
    if (agent.natoms() == 0)
    {
      continue;
    }

    const std::string key = FormName(job_options.agent_stem, i);

    // We really do not know if an Agent is a catalyst, solvent, or what?
    ord::ReactionInput reaction_input =
        MakeReactionInput(agent.smiles(), ord::Compound::ReactionRole::CATALYST);

    (*ord_proto.mutable_inputs())[key] = std::move(reaction_input);
  }

  ord::ReactionOutcome* reaction_outcome = ord_proto.add_outcomes();

  for (int i = 0; i < rxn.number_products(); ++i)
  {
    ISIS_RXN_FILE_Molecule& product = rxn.product(i);
    if (product.natoms() == 0)
    {
      continue;
    }

    ord::ReactionProduct* ord_product = reaction_outcome->add_products();

    ord::CompoundIdentifier* identifier = ord_product->mutable_compound()->add_identifiers();
    identifier->set_type(ord::CompoundIdentifier::SMILES);
    identifier->set_value(ToStdString(product.smiles()));
    for (float yield : yields) {
      ord_product->mutable_compound_yield()->set_value(yield);
    }
  }

  return ord_proto;
}

}  // namespace rxn_to_openrxn
