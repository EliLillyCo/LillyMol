// Builds a Mol2Graph proto

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Lib/mol2graph.pb.h"

using std::cerr;
using std::endl;

int BuildMol2Graph(Command_Line & cl, const char flag, const int verbose,
        LLYMol::Mol2Graph& destination) {

  // Until we get cis-trans canonicalization worked out.
  destination.set_revert_all_directional_bonds_to_non_directional(true);

  const_IWSubstring h;
  for (int i = 0; cl.value(flag, h, i); ++i)
  {
    if ('s' == h)
    {
      destination.set_preserve_cc_double_bonds_saturated(true);
      if (verbose)
        cerr
            << "During graph reduction, will preserve double bonds adjacent to fully saturated atoms\n";
    }
    else if ('c' == h)
    {
      destination.set_preserve_cc_double_bonds_no_heteroatoms(true);
      if (verbose)
        cerr
            << "During graph reduction, will preserve double bonds adjacent to only carbon atoms\n";
    }
    else if ("kptriple" == h) {
      destination.set_exclude_triple_bonds_from_graph_reduction(true);
      if (verbose)
        cerr << "Triple bonds not changed during graph formation\n";
    }
    else if ("chiral" == h)
    {
      destination.set_remove_chiral_centres(false);
      if (verbose)
        cerr << "Will NOT remove chiral centres during graph reduction\n";
    }
    else if ("rmchiral" == h)
    {
      destination.set_remove_chiral_centres(true);
      if (verbose)
        cerr << "Will remove chiral centres during graph reduction\n";
    }
    else if ("help" == h)
    {
      // clang-format off
      cerr << "Mol2Graph::construct:the following options are recognised\n";
      cerr << " -" << flag << " s          preserve C=C bonds that are adjacent to fully saturated atoms\n";
      cerr << " -" << flag << " c          preserve C=C bonds that are adjacent to only carbon atoms\n";
      cerr << " -" << flag << " kptriple   preserve triple bonds during graph formation\n";
      cerr << " -" << flag << " chiral     preserve chiral centres during graph reduction\n";
      cerr << " -" << flag << " rmchiral   remove chiral centres during graph reduction\n";
      // clang-format on
      exit(0);
    }
    else
    {
      cerr << "Mol2Graph::construct:unrecognised -" << flag << " qualifier '" << h << "'\n";
      return 0;
    }
  }

  destination.set_active(true);

  return 1;
}

Mol2Graph Mol2GraphFromProto(const LLYMol::Mol2Graph& mol2graph) {
  Mol2Graph to_be_returned;
  if (mol2graph.exclude_triple_bonds_from_graph_reduction()) {
    to_be_returned.set_exclude_triple_bonds_from_graph_reduction(true);
  }

  if (mol2graph.remove_chiral_centres()) {
    to_be_returned.set_remove_chiral_centres(true);
  }

  if (mol2graph.revert_all_directional_bonds_to_non_directional()) {
    to_be_returned.set_revert_all_directional_bonds_to_non_directional(true);
  }

  if (mol2graph.preserve_cc_double_bonds_saturated()) {
    to_be_returned.set_preserve_cc_double_bonds_saturated(true);
  }

  if (mol2graph.preserve_cc_double_bonds_no_heteroatoms()) {
    to_be_returned.set_preserve_cc_double_bonds_no_heteroatoms(true);
  }

  if (mol2graph.append_molecular_formula()) {
    to_be_returned.set_append_molecular_formula(true);
  }

  if (mol2graph.aromatic_distinguishing_formula()) {
    to_be_returned.set_aromatic_distinguishing_formula(true);
  }

  to_be_returned.set_active(true);

  return to_be_returned;
}

LLYMol::Mol2Graph Mol2GraphToProto(const Mol2Graph& mol2graph) {
  LLYMol::Mol2Graph to_be_returned;

  if (mol2graph.exclude_triple_bonds_from_graph_reduction()) {
    to_be_returned.set_exclude_triple_bonds_from_graph_reduction(true);
  }

  if (mol2graph.remove_chiral_centres()) {
    to_be_returned.set_remove_chiral_centres(true);
  }

  if (mol2graph.revert_all_directional_bonds_to_non_directional()) {
    to_be_returned.set_revert_all_directional_bonds_to_non_directional(true);
  }

  if (mol2graph.preserve_cc_double_bonds_saturated()) {
    to_be_returned.set_preserve_cc_double_bonds_saturated(true);
  }

  if (mol2graph.preserve_cc_double_bonds_no_heteroatoms()) {
    to_be_returned.set_preserve_cc_double_bonds_no_heteroatoms(true);
  }

  if (mol2graph.append_molecular_formula()) {
    to_be_returned.set_append_molecular_formula(true);
  }
  
  if (mol2graph.aromatic_distinguishing_formula()) {
    to_be_returned.set_aromatic_distinguishing_formula(true);
  }

  to_be_returned.set_active(true);

  return to_be_returned;
}
