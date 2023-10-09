#define COMPILING_RXN_FILE

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "atom_typing.h"
#include "misc2.h"
#include "path.h"
#include "rxn_file.h"
#include "smiles.h"
#include "toggle_kekule_form.h"

using std::cerr;
using std::endl;

template <typename T>
int write_rxn_smiles(const Reaction_Smiles_Options &, T &output);

Reaction_Smiles_Options::Reaction_Smiles_Options() {
  _reagent_product_plus_rather_than_dot = true;
  //_orphan_plus_rather_than_dot = 1;
  _write_reaction_name = 1;
  _write_agent = 1;
  _output_separator = ' ';
}

template <typename T>
void
write_set_of_reactions(ISIS_RXN_FILE_Molecule *m, const int n,
                       const bool plus_rather_than_dot, T &output) {
  for (int i = 0; i < n; ++i) {
    if (0 == i) {
      ;
    } else if (plus_rather_than_dot) {
      output << '+';
    } else {
      output << '.';
    }

    m[i].invalidate_smiles();  // make sure we do NOT get unique smiles
    output << m[i].smiles();
  }

  return;
}

template <typename T>
int
RXN_File::write_rxn_smiles(const Reaction_Smiles_Options &opts, T &output) {
  // cerr << "RXN_File::write_rxn_smiles:reaction has " << _nr << " reagents and " << _np
  // << " products. plus_rather_than_dot " << plus_rather_than_dot << endl;
  write_set_of_reactions(_reagent, _nr, opts.reagent_product_plus_rather_than_dot(),
                         output);

  if (_orphan_atoms.natoms()) {
    if (opts.reagent_product_plus_rather_than_dot()) {
      output << '+' << _orphan_atoms.smiles();
    } else {
      output << '.' << _orphan_atoms.smiles();
    }
  }

  output << '>';

  if (opts.write_agent()) {
    write_set_of_reactions(_agent, _na, opts.reagent_product_plus_rather_than_dot(),
                           output);
  }

  output << '>';

  write_set_of_reactions(_product, _np, opts.reagent_product_plus_rather_than_dot(),
                         output);

  if (opts.write_reaction_name()) {
    output << opts.output_separator() << _comment << '\n';
  }

  return 1;
}

template int RXN_File::write_rxn_smiles(const Reaction_Smiles_Options &,
                                        IWString_and_File_Descriptor &);
template int RXN_File::write_rxn_smiles(const Reaction_Smiles_Options &, std::ostream &);

// Given a reaction, determine whether or not + is present as a
// component separator. As long as it is there outside a smiles []
// then it is counted as a component separator.
//    for (const auto& c : rxn_smiles)  does not work, not sure why
static int
PlusPresentAsComponentSeparator(const const_IWSubstring &rxn_smiles) {
  int in_square_brackets = 0;
  const int n = rxn_smiles.length();
  for (int i = 0; i < n; ++i) {
    const char c = rxn_smiles[i];

    if (c == '[') {
      in_square_brackets = 1;
    } else if (in_square_brackets) {
      if (c == ']') {
        in_square_brackets = 0;
      }
    } else if ('+' == c) {
      return 1;
    }
  }

  return 0;
}

/*
  When building a reaction smiles, [N:1] is valid. But under strict smiles interpretation,
  this is a Nitrogen atom with zero Hydrogens. We relax that.
*/

static void
unset_implicit_hydrogen_information(Molecule &m) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (m.implicit_hydrogens_known(i)) {
      m.unset_all_implicit_hydrogen_information(i);
    }
  }

  return;
}

#ifdef NOLONGERNEEDED_ASDASD
template <typename T>
int
build_molecules(ISIS_RXN_FILE_Molecule *m, const resizable_array_p<T> &component) {
  const int n = component.number_elements();

  for (int i = 0; i < n; ++i) {
#ifdef DEBUG_BUILD_MOLECULES_q
    cerr << " i = " << i << " cstart " << cstart << " cstop " << cstop << " '" << s
         << "'\n";
#endif

    if (!m[i].build_from_smiles(*component[i])) {
      cerr << "build_molecules::cannot interpret smiles '" << *component[i] << "'\n";
      return 0;
    }

    m[i].transfer_atom_map_data_to_global_array();
    m[i].allocate_mdl_atom_array();

    unset_implicit_hydrogen_information(m[i]);
  }

  return 1;
}
#endif

// Given smiles in 'component' set 'n' and 'molecule' and then build
// the molecules from the smiles in 'component'.
static int
build_molecules(const resizable_array_p<const_IWSubstring> &component, int &n,
                ISIS_RXN_FILE_Molecule *&molecule) {
  n = component.number_elements();

  if (n == 0) {
    molecule = nullptr;
    return 1;
  }

  molecule = new ISIS_RXN_FILE_Molecule[n];

  for (int i = 0; i < n; ++i) {
    if (!molecule[i].build_from_smiles(*component[i])) {
      cerr << "build_molecules::cannot interpret smiles '" << *component[i] << "'\n";
      return 0;
    }

    molecule[i].transfer_atom_map_data_to_global_array();
    molecule[i].allocate_mdl_atom_array();

    unset_implicit_hydrogen_information(molecule[i]);
  }

  return 1;
}

#ifdef NOLONGERNEEDED_ASDASD
static int
build_molecules(ISIS_RXN_FILE_Molecule *m, const int n, const const_IWSubstring &buffer,
                const resizable_array<int> &pos) {
#ifdef DEBUG_BUILD_MOLECULES_q
  cerr << "Processing " << buffer << endl;
  cerr << "           ";
  for (int i = 0; i < buffer.length(); ++i) {
    cerr << (i % 10);
  }
  cerr << endl;
#endif

  for (int i = 0; i < n; ++i) {
    int cstart;
    if (0 == i) {
      cstart = 0;
    } else {
      cstart = pos[i - 1] + 1;
    }

    int cstop;

    if (i == n - 1) {
      cstop = buffer.length() - 1;
    } else {
      cstop = pos[i] - 1;
    }

    const_IWSubstring s;

    buffer.from_to(cstart, cstop, s);

#ifdef DEBUG_BUILD_MOLECULES_q
    cerr << " i = " << i << " cstart " << cstart << " cstop " << cstop << " '" << s
         << "'\n";
#endif

    if (!m[i].build_from_smiles(s)) {
      cerr << "build_molecules::cannot interpret smiles '" << s << "'\n";
      return 0;
    }

    m[i].transfer_atom_map_data_to_global_array();
    m[i].allocate_mdl_atom_array();

    unset_implicit_hydrogen_information(m[i]);
  }

  return 1;
}
#endif

// #define DEBUG_BUILD_FROM_REACTION_SMILES

int
RXN_File::build_from_reaction_smiles(const const_IWSubstring &buffer,
                                     int component_grouping_is_plus) {
  _nr = 0;
  _na = 0;
  _np = 0;

  int i = 0;
  IWString smiles;

  buffer.nextword(smiles, i);  // skip over smiles, but store it

  // If there is a plus present, override the argument.
  if (component_grouping_is_plus) {
    ;
  } else if (PlusPresentAsComponentSeparator(smiles)) {
    component_grouping_is_plus = 1;
  }

#ifdef DEBUG_BUILD_FROM_REACTION_SMILES
  cerr << "RXN_File::build_from_reaction_smiles:smiles '" << smiles << " plus? "
       << component_grouping_is_plus << endl;
#endif

  ChemAxonFragmentData chemaxon_fragment_data;

  const_IWSubstring token;
  if (!buffer.nextword(token, i)) {  // strange, reaction with no name
    ;
  } else if (token.starts_with('|') && token.ends_with('|')) {
    if (component_grouping_is_plus) {
      cerr << "RXN_File::build_from_reaction_smiles:plus present as component grouping, "
              "and Chemaxon extension, impossible\n";
      return 0;
    }

    if (!_parse_ChemAxon_Extensions(token, chemaxon_fragment_data)) {
      cerr << "RXN_File::build_from_reaction_smiles:invalid ChemAxon extensions '"
           << token << "'\n";
      return 0;
    }
    buffer.nextword(_comment, i);
  } else {
    _comment = token;
  }

  while (buffer.nextword(token, i)) {
    _comment.append_with_spacer(token);
  }

#ifdef DEBUG_BUILD_FROM_REACTION_SMILES
  cerr << "COmment set to '" << _comment << "', chemaxon ? "
       << chemaxon_fragment_data.active() << endl;
#endif

  // If Chemaxon fragment grouping is present, create a new smiles with + signs.
  if (chemaxon_fragment_data.active()) {
    IWString tmp;
    chemaxon_fragment_data.ChangeToPlusForm(smiles, '+', tmp);
    smiles = tmp;
    component_grouping_is_plus = 1;
    //  cerr << "New smiles '" << smiles << "'\n";
  }

  // What character is the component separator?
  char component_separator = '.';
  if (component_grouping_is_plus) {
    component_separator = '+';
  }

  // ReactionStringRep needs a const_IWSubstring, so avoid any problems with
  // temporary variables by making our own here.
  const_IWSubstring _smiles(smiles);

  ReactionStringRep reaction_string_rep;
  if (!reaction_string_rep.Build(_smiles, component_separator)) {
    cerr << "RXN_File::build_from_reaction_smiles:cannot make string form '" << smiles
         << "'\n";
    return 0;
  }

  if (!build_molecules(reaction_string_rep.reagents(), _nr, _reagent)) {
    cerr << "RXN_File::build_from_reaction_smiles:invalid reagents '" << smiles << "'\n";
    return 0;
  }

  if (!build_molecules(reaction_string_rep.agents(), _na, _agent)) {
    cerr << "RXN_File::build_from_reaction_smiles:invalid agents '" << smiles << "'\n";
    return 0;
  }

  if (!build_molecules(reaction_string_rep.products(), _np, _product)) {
    cerr << "RXN_File::build_from_reaction_smiles:invalid products '" << smiles << "'\n";
    return 0;
  }

  _do_automatic_atom_mapping = 0;

  // If something has a chiral centre, we may be able to infer something about the bonding
  // Turns out this is wrong. The unspecified connection in the product may be an
  // unspecified

  return 1;
}

ChemAxonFragmentData::ChemAxonFragmentData() {
  _only_fragents_that_are_grouped_become_components = 1;
}

ChemAxonFragmentData::~ChemAxonFragmentData() {
}

int
ChemAxonFragmentData::_add_group(const const_IWSubstring &buffer) {
  std::unique_ptr<ConnectedGroup> group = std::make_unique<ConnectedGroup>();

  // cerr << "ChemAxonFragmentData::_add_group: processing '" << buffer << "'\n";
  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i, '.')) {
    //  cerr << "Token is '" << token << "'\n";
    if (token.empty()) {
      cerr << "ChemAxonFragmentData::_add_group:empty connected group '" << buffer
           << "'\n";
      return 0;
    }

    int k;
    if (!token.numeric_value(k) || k < 0) {
      cerr << "ChemAxonFragmentData::_add_group:invalid fragment '" << buffer << "'\n";
      return 0;
    }

    group->add(k);
  }

  if (group->number_elements() < 2) {
    cerr << "ChemAxonFragmentData::_add_group:invalid group '" << buffer << "'\n";
    return 0;
  }

  _connected_group.add(group.release());

  return 1;
}

int
ChemAxonFragmentData::Build(const const_IWSubstring &buffer) {
  if (!_connected_group.empty()) {
    _connected_group.resize_keep_storage(0);
  }

  if (buffer.empty()) {
    return 1;
  }

  int f_seen = 0;

  int i = 0;
  const_IWSubstring token;
  int previous_was_f = 0;
  while (buffer.nextword(token, i, ',')) {
    //  cerr << "ChemAxonFragmentData::Build token '" << token << "'\n";
    if (token.starts_with("f:")) {
      token.remove_leading_chars(2);
      previous_was_f = 1;
      if (!_add_group(token)) {
        cerr << "ChemAxonFragmentData::Build:invalid group '" << buffer << "'\n";
        return 0;
      }
      f_seen++;
    } else if (previous_was_f && isdigit(token[0])) {  // && token.contains('.') &&
                                                       // isdigit(token.last_item())) {
      if (!_add_group(token)) {
        cerr << "ChemAxonFragmentData::Build:invalid continuation group '" << buffer
             << "'\n";
        return 0;
      }
    } else {
      previous_was_f = 0;
    }
  }

  if (!f_seen) {
    return 1;
  }
  //{
  //  return 0;
  //}

  // Check for dupes

  resizable_array<int> seen;
  for (const auto *group : _connected_group) {
    for (const int i : *group) {
      if (seen.contains(i)) {
        cerr << "ChemAxonFragmentData::Build:duplicate fragment '" << buffer << "'\n";
        return 0;
      }
      seen.add(i);
    }
  }

  // Cannot return number of groups since zero fragment groups is OK.
  return 1;
}

void
ChemAxonFragmentData::debug_print(std::ostream &output) const {
  output << "ChemAxonFragmentData:with " << _connected_group.size()
         << " connected groups\n";
  for (int i = 0; i < _connected_group.number_elements(); ++i) {
    output << "  group " << i << " :";
    for (auto j : *_connected_group[i]) {
      output << ' ' << j;
    }
    output << "\n";
  }
}

int
ChemAxonFragmentData::ChangeToPlusForm(const const_IWSubstring &buffer, char sep,
                                       IWString &new_rxnsmiles) {
  const_IWSubstring reagents, agents, products;
  int i = 0;
  buffer.nextword_single_delimiter(reagents, i, '>');
  buffer.nextword_single_delimiter(agents, i, '>');
  buffer.nextword_single_delimiter(products, i, '>');
#ifdef DEBUG_CHANGE_TO_PLUS_FORM
  cerr << "reagents " << reagents << endl;
  cerr << "agents " << agents << endl;
  cerr << "products " << products << endl;
  debug_print(cerr);
#endif

  if (reagents.empty() || products.empty()) {
    cerr << "ChemAxonFragmentData::ChangeToPlusForm:empty reagents " << reagents
         << " and/or products " << products << endl;
    return 0;
  }

  resizable_array_p<const_IWSubstring> tokens;
  resizable_array<RAP> rap;

  i = 0;
  const_IWSubstring token;
  while (reagents.nextword(token, i, '.')) {
    tokens.add(new const_IWSubstring(token));
    rap.add(Reagent);
  }
  i = 0;
  while (agents.nextword(token, i, '.')) {
    tokens.add(new const_IWSubstring(token));
    rap.add(Agent);
  }
  i = 0;
  while (products.nextword(token, i, '.')) {
    tokens.add(new const_IWSubstring(token));
    rap.add(Product);
  }

  if (!_only_fragents_that_are_grouped_become_components) {
    _create_groups_for_ungrouped(rap);
  }

  const int n = tokens.number_elements();

#ifdef DEBUG_CHANGE_TO_PLUS_FORM
  for (int i = 0; i < n; ++i) {
    cerr << "token " << i << " is " << *tokens[i] << " rap " << rap[i] << endl;
  }
#endif

  new_rxnsmiles.resize_keep_storage(0);
  new_rxnsmiles.resize(buffer.length());

  // At this stage, there are two ways of processing the input.
  if (_only_fragents_that_are_grouped_become_components) {
    _create_groups_for_ungrouped(rap);
  }

  // Keep track of which tokens have been written already.
  int *done = new_int(n);
  std::unique_ptr<int[]> free_done(done);

  bool first_in_component = true;
  for (int i = 0; i < n; ++i) {
    //  cerr << " i = " << i << " result now " << new_rxnsmiles << endl;
    if (i > 0 && rap[i] != rap[i - 1]) {
      new_rxnsmiles << '>';
      // Handle the case of no agents.
      if (rap[i - 1] == Reagent && rap[i] == Product) {
        new_rxnsmiles << '>';
      }
      first_in_component = true;
    }

    if (done[i]) {
      continue;
    }

    done[i] = 1;

    if (first_in_component) {
      first_in_component = false;
    } else {
      new_rxnsmiles << '+';
    }

    //  cerr << "Looking for group containing token " << i << endl;

    new_rxnsmiles << *tokens[i];
    for (const auto &group : _connected_group) {
      if (!group->contains(i)) {
        continue;
      }

      for (const auto f : *group) {
        //      cerr << "  in group with " << i << " is " << f << endl;
        if (f == i) {
          continue;
        }

        if (!tokens.ok_index(f)) {
          cerr << "ChemAxonFragmentData::ChangeToPlusForm:group " << f << " but have "
               << tokens.number_elements() << " tokens\n";
          return 0;
        }
        new_rxnsmiles << '.' << *tokens[f];
        done[f] = 1;
      }
    }
  }

  return 1;
}

int
ChemAxonFragmentData::_is_in_group(const int f) const {
  for (const auto &group : _connected_group) {
    if (group->contains(f)) {
      return 1;
    }
  }

  return 0;
}

// If we are operating in the mode of ungrouped fragments become one component,
// then identify those tokens that are not part of a grouping, and for each
// component, possibly add a new ConnectedGroup that contains all the items
// not in a group.
// Note that this changes _connected_group.
int
ChemAxonFragmentData::_create_groups_for_ungrouped(const resizable_array<RAP> &rap) {
  const int n = rap.number_elements();
  // For each type, those tokens that are not part of a connected group.
  std::unique_ptr<ConnectedGroup> reagents = std::make_unique<ConnectedGroup>();
  std::unique_ptr<ConnectedGroup> agents = std::make_unique<ConnectedGroup>();
  std::unique_ptr<ConnectedGroup> products = std::make_unique<ConnectedGroup>();

  for (int i = 0; i < n; ++i) {
    if (_is_in_group(i)) {
      continue;
    }

    if (rap[i] == Reagent) {
      reagents->add(i);
    } else if (rap[i] == Agent) {
      agents->add(i);
    } else if (rap[i] == Product) {
      products->add(i);
    }
  }

  if (!reagents->empty() && reagents->number_elements() > 1) {
    _connected_group.add(reagents.release());
  }
  if (!agents->empty() && agents->number_elements() > 1) {
    _connected_group.add(agents.release());
  }
  if (!products->empty() && products->number_elements() > 1) {
    _connected_group.add(products.release());
  }

  // cerr << "After adding extra ";
  // debug_print(cerr);

  return 1;
}

ReactionStringRep::ReactionStringRep() {
}

void
ReactionStringRep::_reset() {
  _reagent.resize(0);
  _agent.resize(0);
  _product.resize(0);
}

ReactionStringRep::~ReactionStringRep() {
  _reset();
}

int
ReactionStringRep::Build(const const_IWSubstring &rxn_smiles,
                         const char component_separator) {
  // Reagents, agents and products.
  const_IWSubstring r, a, p;
  int i = 0;
  if (!rxn_smiles.nextword_single_delimiter(r, i, '>') ||
      !rxn_smiles.nextword_single_delimiter(a, i, '>') ||
      !rxn_smiles.nextword_single_delimiter(p, i, '>')) {
    cerr << "ReactionStringRep:Build:invalid reaction '" << rxn_smiles << "'\n";
    return 0;
  }

  if (!SplitOnPlusses(r, _reagent, component_separator)) {
    cerr << "ReactionStringRep::Build:invalid reagents '" << rxn_smiles << "'\n";
    return 0;
  }

  if (!a.empty() && !SplitOnPlusses(a, _agent, component_separator)) {
    cerr << "ReactionStringRep::Build:invalid agents '" << rxn_smiles << "'\n";
    return 0;
  }

  if (!SplitOnPlusses(p, _product, component_separator)) {
    cerr << "ReactionStringRep::Build:invalid product '" << rxn_smiles << "'\n";
    return 0;
  }

  return 1;
}

// For now, just looking at fragment grouping, "f:0.1,5.6"
// Note that if there is some other directive and the fragment grouping is
// not the first, this will fail.
int
RXN_File::_parse_ChemAxon_Extensions(const_IWSubstring buffer,  // note local copy
                                     ChemAxonFragmentData &cxfd) {
  assert(buffer.starts_with('|') && buffer.ends_with('|'));

  buffer.remove_leading_chars(1);
  buffer.chop();

  return cxfd.Build(buffer);
}

/*
  The ISIS_RXN_FILE_Molecule has an array of atom maps. If the object has been built from
  a smiles, we need to transfer that atom map information to the central array
*/

int
ISIS_RXN_FILE_Molecule::transfer_atom_map_data_to_global_array() {
  const int matoms = natoms();

  if (0 == matoms) {
    return 0;
  }

  if (nullptr == _atom_map) {
    _atom_map = new int[matoms];
  }

  for (int i = 0; i < matoms; ++i) {
    const Atom *a = atomi(i);

    _atom_map[i] = a->atom_map();

    //  cerr << "  atom " << i << " transfered " << _atom_map[i] << endl;
  }

  return 1;
}

/*
  We have built the object from a smiles, but none of the underlying MDL_Molecule things
  are present We also allocate the bond array...
*/

int
ISIS_RXN_FILE_Molecule::allocate_mdl_atom_array() {
  const int matoms = natoms();

  _mdl_atom.resize(matoms);
  _nbonds = new_int(matoms);

  for (int i = 0; i < matoms; ++i) {
    MDL_Atom_Data *x = new MDL_Atom_Data();

    x->set_atom_number(i);
    x->set_atom_map(atom_map(i));

    _mdl_atom.add(x);

    _nbonds[i] = nbonds(i);
  }

  const int nbonds = nedges();

  _mdl_bond.resize(nbonds);

  for (int i = 0; i < nbonds; ++i) {
    MDL_Bond_Data *x = new MDL_Bond_Data();

    const Bond *b = bondi(i);

    x->set_btype(b->btype());

    _mdl_bond.add(x);
  }

  _compute_implicit_hydrogens();

  return 1;
}

int
Atom_Locations::initialise(ISIS_RXN_FILE_Molecule *reagents, const int nr,
                           ISIS_RXN_FILE_Molecule *products, const int np,
                           const int *product_locator) {
  for (int i = 0; i < nr; ++i) {
    const int matoms = reagents[i].natoms();

    const int *amap = reagents[i].atom_map();

    for (int j = 0; j < matoms; ++j) {
      const int x = amap[j];

      const int p = product_locator[x];

      if (p < 0) {
        continue;
      }

      const auto f = _reagent_to_product.find(x);

      if (f != _reagent_to_product.end()) {
        cerr << "Atom_Locations::initialise:atom map " << x
             << " found in multiple reagents, impossible\n";
        return 0;
      }

      unsigned int u;
      unsigned short *s = reinterpret_cast<unsigned short *>(&u);

      s[0] = p;
      s[1] = products[p].which_is_mapped_atom(x);

      _reagent_to_product[x] = u;
    }
  }

  return 1;
}

int
Atom_Locations::get_product_and_atom_number(const int amap, int &product,
                                            atom_number_t &zatom) const {
  const auto f = _reagent_to_product.find(amap);

  if (f == _reagent_to_product.end()) {
    return 0;
  }

  const auto u = f->second;

  const unsigned short *s = reinterpret_cast<const unsigned short *>(&u);

  product = s[0];
  zatom = s[1];

  return 1;
}

int
RXN_File::contains_isotopic_reagent_atoms() const {
  for (int i = 0; i < _nr; ++i) {
    if (_reagent[i].number_isotopic_atoms() > 0) {
      return 1;
    }
  }

  return 0;
}

int
RXN_File::contains_isotopic_product_atoms() const {
  for (int i = 0; i < _np; ++i) {
    if (_product[i].number_isotopic_atoms() > 0) {
      return 1;
    }
  }

  return 0;
}

/*
  I encountered reactions where none of the mapped atoms in the LHS showed up in the RHS
*/

int
RXN_File::at_least_some_mapped_atoms_common_btw_reagents_and_products() {
  if (nullptr == _product_locator) {
    _reestablish_reagent_locator_array();
  }

  for (int i = 0; i < _nr; ++i) {
    const int *amap = _reagent[i].atom_map();

    const int matoms = _reagent[i].natoms();

    for (int j = 0; j < matoms; ++j) {
      const int jmap = amap[j];
      //    cerr << " atom " << j << " map " << jmap << " product " <<
      //    _product_locator[jmap] << endl;

      if (_product_locator[jmap] >= 0) {
        return 1;
      }
    }
  }

  cerr << "RXN_File::at_least_some_mapped_atoms_common_btw_reagents_and_products "
          "returning false, nr "
       << _nr << endl;
  return 0;
}

int
RXN_File::_all_atoms_still_aromatic(ISIS_RXN_FILE_Molecule &m, const Ring &r) {
  const int ring_size = r.number_elements();

  const int *amap = m.atom_map();

  for (int i = 0; i < ring_size; ++i) {
    const auto j = r[i];

    const auto amapj = amap[j];
    if (amapj <= 0) {
      continue;
    }

    const auto p = _product_locator[amapj];

    if (p < 0) {
      return 0;
    }

    const auto pj = _product[p].which_is_mapped_atom(amapj);

    if (pj < 0)  //
    {
      cerr << "RXN_File::_all_atoms_still_aromatic:did not find mapped atom " << amapj
           << " in product " << _comment << endl;
      return 0;
    }

    if (!_product[p].is_aromatic(pj)) {
      return 0;
    }
  }

  return 1;  // all atoms are still aromatic
}

/*
  Try to put reactions back into consistent Kekule forms
*/

int
RXN_File::fix_kekule_differences() {
  prepare_for_reaction_construction();

  int rc = 0;

  for (int i = 0; i < _nr; ++i) {
    rc += _fix_kekule_differences_reagent(i);
  }

  return rc;
}

static int
is_duplicate(const atom_number_t a1, const atom_number_t a2,
             const Set_of_Atoms &already_found) {
  for (int i = 0; i < already_found.number_elements(); i += 2) {
    const atom_number_t f1 = already_found[i];
    const atom_number_t f2 = already_found[i + 1];

    if (a1 == f1 && a2 == f2) {
      return 1;
    }

    if (a1 == f2 && a2 == f1) {
      return 1;
    }
  }

  return 0;
}

/*
  This is a minor variant of _all_bonds_unchanged

  Changed to a template function so it can be used for checking both the case
  of two bonds being the same or two bonds being different
*/

template <typename CMP>
int
RXN_File::_identify_atoms_with_different_bonding(
    const ISIS_RXN_FILE_Molecule &m, const Ring &r, atom_number_t &a1, atom_number_t &a2,
    CMP compare, const Set_of_Atoms &already_found) const {
  // cerr << "_identify_atoms_with_different_bonding\n";
  // Molecule mcopy(m);
  // write_isotopically_labelled_smiles(mcopy, 0, cerr);

  const int *amap = m.atom_map();

  for (Ring_Bond_Iterator i(r); i != r.zend(); i++) {
    const atom_number_t ra1 = i.a1();
    const atom_number_t ra2 = i.a2();

    const int amap1 = amap[ra1];
    if (amap1 < 1) {
      continue;
    }

    const int amap2 = amap[ra2];
    if (amap2 < 1) {
      continue;
    }

    const int p1 = _product_locator[amap1];
    const int p2 = _product_locator[amap2];

    //  cerr << "_identify_atoms_with_different_bonding:considering atom " << ra1 << "
    //  (map " << amap[ra1] << ") and " << ra2 << " (map " << amap[ra2] << "), products "
    //  << p1 << " and " << p2 << endl;

    if (p1 < 0 || p1 != p2) {
      return 0;
    }

    const atom_number_t pa1 = _product[p1].which_is_mapped_atom(amap1);
    const atom_number_t pa2 = _product[p1].which_is_mapped_atom(amap2);

    const Bond *b1 = m.bond_between_atoms_if_present(ra1, ra2);
    if (nullptr == b1) {
      return 0;
    }

    const Bond *b2 = _product[p1].bond_between_atoms_if_present(pa1, pa2);
    if (nullptr == b2) {
      return 0;
    }

    if (compare(b1->btype(), b2->btype())) {
      continue;
    }
    //  if (b1->btype() == b2->btype())
    //    continue;

    if (is_duplicate(ra1, ra2, already_found)) {
      continue;
    }

    a1 = ra1;
    a2 = ra2;

    //  cerr << "Bond changes, bonds " << *b1 << " and " << *b2 << endl;

    return 1;
  }

  return 0;  // did not detect any ring bonds with different bonds
}

#ifdef NO_LONGER_USED_____________
static void
copy_bonding(const Molecule &mfrom, Molecule &mto) {
  const auto nb = mfrom.nedges();

  assert(nb == mto.nedges());

  for (int i = 0; i < nb; ++i) {
    const Bond *b1 = mfrom.bondi(i);
    const Bond *b2 = mto.bondi(i);

    if (BOND_TYPE_ONLY(b1->btype()) == BOND_TYPE_ONLY(b2->btype())) {
      continue;
    }

    if (b1->is_single_bond()) {
      mto.set_bond_type_between_atoms(b2->a1(), b2->a2(), DOUBLE_BOND);
    } else {
      mto.set_bond_type_between_atoms(b2->a1(), b2->a2(), SINGLE_BOND);
    }
  }

  return;
}
#endif

static void
grow_aromatic_system(Molecule &m, const int rng, int *fused_aromatic_system_id,
                     const int id) {
  const Ring *r = m.ringi(rng);

  const int nfused = r->fused_ring_neighbours();

  for (int i = 0; i < nfused; ++i) {
    const Ring *ri = r->fused_neighbour(i);

    const int ring_number = ri->ring_number();

    if (fused_aromatic_system_id[ring_number] > 0) {
      continue;
    }

    if (!ri->is_aromatic()) {
      continue;
    }

    fused_aromatic_system_id[ring_number] = id;

    grow_aromatic_system(m, ring_number, fused_aromatic_system_id, id);
  }
}

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
static int
identify_fused_aromatic_systems(Molecule &m, int *fused_aromatic_system_id) {
  const int nr = m.nrings();

  int id = 0;

  for (int i = 0; i < nr; ++i) {
    if (fused_aromatic_system_id[i] > 0) {
      continue;
    }

    const Ring *ri = m.ringi(i);

    if (!ri->is_aromatic()) {
      continue;
    }

    id++;
    fused_aromatic_system_id[i] = id;

    if (ri->is_fused()) {
      grow_aromatic_system(m, i, fused_aromatic_system_id, id);
    }
  }

  return id;  // the number of groups we identify
}
#endif

/*
  for each aromatic ring, how many fused aromatic neighbours does it have
*/

#ifdef QWEQWE_NOT_USED
static void
identify_fused_aromatic_neighbours(Molecule &m, int *fused_aromatic_neighbours) {
  const int nr = m.nrings();

  for (int i = 0; i < nr; ++i) {
    const Ring *ri = m.ringi(i);

    if (!ri->is_aromatic()) {
      continue;
    }

    if (!ri->is_fused()) {
      continue;
    }

    const int nfused = ri->fused_ring_neighbours();

    for (int j = 0; j < nfused; ++j) {
      const Ring *nj = ri->fused_neighbour(j);
      if (nj->is_aromatic()) {
        fused_aromatic_neighbours[i]++;
      }
    }
  }

  return;
}
#endif

static int
alternating_single_and_double_bonds(Molecule &m, const Ring &r) {
  bond_type_t previous_bond = INVALID_BOND_TYPE;

  for (Ring_Bond_Iterator i(r); i != r.zend(); i++) {
    const atom_number_t ra1 = i.a1();
    const atom_number_t ra2 = i.a2();

    const Bond *b = m.bond_between_atoms(ra1, ra2);

    if (INVALID_BOND_TYPE == previous_bond) {
      previous_bond = b->btype();
    } else if (b->btype() == previous_bond) {
      return 0;
    } else {
      previous_bond = b->btype();
    }
  }

  return 1;
}

/*
  things with =O and =N bonds outside the ring cannot have Kekule forms
*/

static int
exocyclic_non_aromatic_multiple_bond(Molecule &m, const atom_number_t zatom) {
  const Atom *a = m.atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i) {
    const Bond *b = a->item(i);

    if (b->is_aromatic()) {
      continue;
    }

    if (b->is_single_bond()) {
      continue;
    }

    return 1;  // found exocyclic, non aromatic double bond
  }

  return 0;  // none found on this atom
}

/*
  Consider
  [C:1]1(=[CH:2][CH:3]=[C:6]2[C:5](=[CH:4]1)[N:10]=[C:9]([Cl:25])[C:8](=[N:7]2)[CH3:11])[C:12]([F:13])([F:14])[F:15]>>[C:1]1(=[CH:4][C:5]2=[N:10][C:9](=[C:8]([CH3:11])[N:7]=[C:6]2[CH:3]=[CH:2]1)[OH:16])[C:12]([F:13])([F:14])[F:15]
  82J-E03701-086

  the ring is not alternating single and double, but can take Kekule forms
*/

static int
could_have_kekule_forms(Molecule &m, const Ring &r) {
  if (alternating_single_and_double_bonds(m, r)) {
    return 1;
  }

  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t j = r[i];

    const Atom *aj = m.atomi(j);

    const int jcon = aj->ncon();

    if (7 == aj->atomic_number()) {
      if (3 == jcon || m.hcount(j)) {
        return 0;
      }
    } else if (8 == aj->atomic_number()) {
      return 0;
    } else if (6 == aj->atomic_number() && 3 == jcon) {
      if (exocyclic_non_aromatic_multiple_bond(m, j)) {
        return 0;
      }
    }
  }

  return 1;
}

/*
  This is complicated by the need to do it ring-system at a time
*/

// #define DEBUG_FIX_KEKULE_DIFFERENCES

int
RXN_File::_fix_kekule_differences_reagent(const int rgnt) {
  ISIS_RXN_FILE_Molecule &r = _reagent[rgnt];

  if (!r.contains_aromatic_atoms()) {
    return 0;
  }

  int rc = 0;

  const int nr = r.nrings();

  int *ring_already_done = new_int(nr + nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);
  int *fused_aromatic_system = ring_already_done + nr;  // re-use previous allocation

  r.compute_aromaticity_if_needed();

  int aromatic_rings_to_be_processed = 0;

  for (int i = 0; i < nr; ++i) {
    const Ring *ri = r.ringi(i);

    if (!ri->is_aromatic()) {
      ring_already_done[i] = 1;
      continue;
    }

    if (6 != ri->number_elements()) {  // only 6 membered rings can toggle
      continue;
    }

    if (!_all_atoms_still_aromatic(r, *ri))  // aromaticity lost
    {
      ring_already_done[i] = 1;
      continue;
    }

#ifdef _MAYBE_WE_DO_NOT_NEED_TO_BE_SO_CAREFUL
    if (_involves_exocyclic_bonding_changes(rgnt,
                                            *ri))  // cannot handle Kekule changes here
    {
      ring_already_done[i] = 1;
      continue;
    }
#endif

    if (!could_have_kekule_forms(r, *ri)) {
      ring_already_done[i] = 1;
      continue;
    }

    aromatic_rings_to_be_processed++;
  }

  if (0 == aromatic_rings_to_be_processed) {
    return 0;
  }

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
  const int number_aromatic_groups = identify_fused_aromatic_systems(
      r, fused_aromatic_system);  // so we can process them together

  cerr << "Contains " << number_aromatic_groups << " aromatic groups\n";
  for (int i = 0; i < nr; ++i) {
    cerr << "fused_aromatic_system " << fused_aromatic_system[i] << " ring "
         << *r.ringi(i) << " ring_already_done " << ring_already_done[i] << endl;
  }
#endif

  for (int i = 0; i < nr; ++i) {
    if (0 == fused_aromatic_system[i] || ring_already_done[i]) {
      continue;
    }

    Toggle_Kekule_Form tkf;
    Set_of_Atoms s;
    int ndx = 0;

    int differences_found = 0;

    for (int j = i; j < nr; ++j)  // note we start with ring I
    {
      if (fused_aromatic_system[i] != fused_aromatic_system[j]) {
        continue;
      }

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
      cerr << "System with " << i << " and " << j << " ring_already_done "
           << ring_already_done[i] << endl;
#endif
      const Ring *rj = r.ringi(j);

      if (6 != rj->number_elements()) {
        continue;
      }

      if (ring_already_done[j]) {
        continue;
      }

      ring_already_done[j] = 1;

      //    cerr << "Looking for changing atoms j = " << j << endl;

      int diff = 0;
      atom_number_t a1, a2;
      if (_identify_atoms_with_different_bonding(r, *rj, a1, a2, std::equal_to<int>(),
                                                 s)) {  // bonding has changed
        diff = 1;
      } else if (_identify_atoms_with_different_bonding(r, *rj, a1, a2,
                                                        std::not_equal_to<int>(),
                                                        s)) {  // bonding the same
        diff = 0;
      } else {
        continue;
      }

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
      cerr << "Changing atoms " << a1 << " and " << a2 << " found, diff? " << diff
           << endl;
#endif

      const Bond *b = r.bond_between_atoms(a1, a2);

      bond_type_t bt;
      if (b->is_single_bond()) {
        if (diff) {
          bt = DOUBLE_BOND;
        } else {
          bt = SINGLE_BOND;
        }
      } else  // bond must be double
      {
        if (diff) {
          bt = SINGLE_BOND;
        } else {
          bt = DOUBLE_BOND;
        }
      }

      tkf.add_bond(ndx, ndx + 1, bt);
      s.add(a1);
      s.add(a2);
      ndx += 2;

      if (diff) {
        differences_found++;
      }

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
      cerr << "Changing bond btw " << a1 << " (map " << r.atom_map()[a1] << ") and " << a2
           << "  (map " << r.atom_map()[a2] << ") to " << bt << " fsid "
           << rj->fused_system_identifier() << endl;
      set_include_atom_map_with_smiles(0);
      write_isotopically_labelled_smiles(r, false, cerr);
      set_include_atom_map_with_smiles(1);
      cerr << endl;
#endif
    }

    if (!tkf.active() || 0 == differences_found) {
      continue;
    }

    int changed = 0;

    tkf.set_display_error_messages(2);
    tkf.set_unset_unnecessary_implicit_hydrogens_known_values(1);

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
    cerr << "Toggling with " << s << endl;
    write_isotopically_labelled_smiles(r, false, cerr);
#endif

    if (!tkf.process(r, s, changed) || !changed) {
      continue;
    }

    r.compute_aromaticity_if_needed();

#ifdef DEBUG_FIX_KEKULE_DIFFERENCES
    cerr << "After toggling " << r.smiles() << endl;
#endif

    rc++;
  }

  return rc;
}

int
RXN_File::_bond_the_same_in_product(const ISIS_RXN_FILE_Molecule &r,
                                    const Bond &b) const {
  const auto ma1 = r.atom_map()[b.a1()];
  const auto ma2 = r.atom_map()[b.a2()];

  const auto p = _product_locator[ma1];
  if (p < 0) {
    return 0;
  }

  if (p != _product_locator[ma2]) {
    return 0;
  }

  const atom_number_t pa1 = _product[p].which_is_mapped_atom(ma1);
  const atom_number_t pa2 = _product[p].which_is_mapped_atom(ma2);

  if (pa1 < 0 || pa2 < 0) {
    return 0;
  }

  const Bond *b2 = _product[p].bond_between_atoms_if_present(pa1, pa2);

  if (nullptr == b2) {
    return 0;
  }

  return (BOND_TYPE_ONLY(b.btype()) == BOND_TYPE_ONLY(b2->btype()));
}

int
RXN_File::_involves_exocyclic_bonding_changes(const int rgnt, const Ring &zring) const {
  ISIS_RXN_FILE_Molecule &r = _reagent[rgnt];

  const int ring_size = zring.number_elements();

  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t j = zring[i];

    //  cerr << "Checking connectivity atom " << j << ' ' <<
    //  r.smarts_equivalent_for_atom(j) << endl;

    if (_same_connectivity_and_bonds_in_product(r, j)) {  // no change
      ;
    } else if (_loss_or_gain_of_singly_bonded_atom(
                   rgnt, r.atom_map()[j])) {  // gain or loss of heavy atom, but Hydrogen
                                              // atom present
      ;
    } else {
      return 1;
    }

#ifdef STUFF_ABOVE_WORKS_BETTER_DISCARD_HERE
    for (int k = 0; k < jcon; ++k) {
      const Bond *b = a->item(k);

      const atom_number_t x = b->other(j);

      if (zring.contains(x)) {  // we are looking for exocyclic changes
        continue;
      }

      if (_bond_the_same_in_product(r, *b)) {
        continue;
      }

      return 1;  // exocyclic bond changed
    }
#endif
  }

  return 0;
}

int
RXN_File::_same_connectivity_and_bonds_in_product(const ISIS_RXN_FILE_Molecule &r,
                                                  const atom_number_t x1) const {
  const int mx1 = r.atom_map()[x1];

  const auto p = _product_locator[mx1];

  if (p < 0) {
    return 0;
  }

  const auto x2 = _product[p].which_is_mapped_atom(mx1);

  assert(x2 >= 0);

  // cerr << "In product " << p << " atom number " << x << " con " << _product[p].ncon(x)
  // << endl;

  if (r.ncon(x1) != _product[p].ncon(x2)) {
    return 0;
  }

  return r.nbonds(x1) == _product[p].nbonds(x2);
}

int
RXN_File::_loss_or_gain_of_singly_bonded_atom(const int rgnt, const int am) const {
  const auto p = _product_locator[am];

  if (p < 0) {
    return 0;
  }

  const auto x1 = _reagent[rgnt].which_is_mapped_atom(am);
  const auto x2 = _product[p].which_is_mapped_atom(am);

  assert(x1 >= 0 && x2 >= 0);

  const Atom *a1 = _reagent[rgnt].atomi(x1);
  const Atom *a2 = _product[p].atomi(x2);

  const int acon1 = a1->ncon();
  const int acon2 = a2->ncon();

  // Complicated by the fact that a reaction might involve changing a charged N atom

  if (acon1 > acon2)  // losing a connection, product must have a Hydrogen atm
  {
    if (_product[p].hcount(x2)) {
      return 1;
    }

    if (7 == a1->atomic_number() && 1 == a1->formal_charge() && 3 == acon1 &&
        2 == acon2 && 0 == a2->formal_charge()) {
      return 1;
    }

    return 0;
  } else if (acon1 < acon2)  // gaining a connection, reagent must have a Hydrogen atom
  {
    if (_reagent[rgnt].hcount(x1)) {
      return 1;
    }

    if (7 == a1->atomic_number() && 1 == a2->formal_charge() && 3 == acon2 &&
        2 == acon1 && 0 == a1->formal_charge()) {
      return 1;
    }

    return 0;
  }

  // Do we have a bonding change

  return a1->nbonds() == a2->nbonds();
}

static int
max_atom_in_any(const int n, ISIS_RXN_FILE_Molecule *m) {
  int rc = 0;

  for (int i = 0; i < n; ++i) {
    if (m[i].natoms() > rc) {
      rc = m[i].natoms();
    }
  }

  return rc;
}

int
RXN_File::max_atom_in_any_reagent() const {
  if (0 == _nr) {
    return 0;
  }

  return max_atom_in_any(_nr, _reagent);
}

int
RXN_File::max_atom_in_any_product() const {
  if (0 == _np) {
    return 0;
  }

  return max_atom_in_any(_np, _product);
}

RXN_File_Create_Reaction_Options::RXN_File_Create_Reaction_Options() {
  _only_create_query_from_first_reagent = 0;
}

static int
_remove_duplicates_ignore_atom_map(ISIS_RXN_FILE_Molecule *reagent_or_product, int &n) {
  set_include_chiral_info_in_smiles(0);
  for (int i = 0; i < n; ++i) {
    reagent_or_product[i].recompute_unique_smiles();
  }
  set_include_chiral_info_in_smiles(1);

  int rc = 0;

  for (int i = 0; i < n; ++i) {
    // cerr << "all parts: ";

    //      for (int k = 0; k < n; ++k)
    //    {
    //      const IWString & tempSmi = reagent_or_product[k].unique_smiles();
    //      cerr << " " << tempSmi;
    //    }
    //    cerr << endl << endl;

    const IWString &iusmi = reagent_or_product[i].unique_smiles();
    // cerr << "Reagent " << i << " usmi " << iusmi << endl;

    for (int j = (i + 1); j < n; ++j) {
      //      cerr << "  reagent " << j << " usmi " <<
      //      reagent_or_product[j].unique_smiles() << endl;

      if (reagent_or_product[j].unique_smiles() != iusmi) {  // definitely not the same
        continue;
      }

      //    cerr << "in " << name() << " duplicate reagents i = " << i << " j = " << j <<
      //    endl;

      for (int k = j; k < (n - 1);
           ++k)  // reagent J is a duplicate of I. Shift everything down
      {
        reagent_or_product[k] = std::move(reagent_or_product[k + 1]);
      }
      n--;
      j--;
      rc++;
    }
  }

  return rc;
}

int
RXN_File::remove_duplicate_reagents_ignore_atom_map() {
  const int rc = _remove_duplicates_ignore_atom_map(_reagent, _nr);

  if (0 == rc) {
    return 0;
  }

  _reestablish_reagent_locator_array();

  return rc;
}

int
RXN_File::remove_duplicate_products_ignore_atom_map() {
  const int rc = _remove_duplicates_ignore_atom_map(_product, _np);

  if (0 == rc) {
    return 0;
  }

  _reestablish_reagent_locator_array();

  return rc;
}

int
RXN_File::remove_duplicate_agents() {
  const int rc = _remove_duplicates_ignore_atom_map(_agent, _na);

  if (0 == rc) {
    return 0;
  }

  _reestablish_reagent_locator_array();

  return rc;
}

int
RXN_File::remove_all_agents() {
  _na = 0;

  _reestablish_reagent_locator_array();

  return 1;
}

static int
_contains_duplicate_atom_map_numbers(const ISIS_RXN_FILE_Molecule &m, int *mapped) {
  const int matoms = m.natoms();

  const int *amap = m.atom_map();

  if (nullptr == amap) {
    return 0;
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    const int mi = amap[i];
    if (mi <= 0) {
      continue;
    }

    if (0 != mapped[mi]) {
      //    cerr << "Atom " << i << " map " << mi << " is a dup\n";
      rc++;
    }

    mapped[mi] += 1;
  }

  return rc;
}

static int
_contains_duplicate_atom_map_numbers(const int n, const ISIS_RXN_FILE_Molecule *m,
                                     int *mapped) {
  int rc = 0;

  for (int i = 0; i < n; ++i) {
    if (_contains_duplicate_atom_map_numbers(m[i], mapped)) {
      return rc = 1;
    }
  }

  return rc;
}

int
RXN_File::contains_duplicate_atom_map_numbers() const {
  const int h = highest_atom_map_number();

  if (h <= 0) {
    return 0;
  }

  int *mapped = new_int(h + 1);
  std::unique_ptr<int[]> free_mapped(mapped);

  if (_contains_duplicate_atom_map_numbers(_nr, _reagent, mapped)) {
    return 1;
  }

  std::fill_n(mapped, h + 1, 0);

  return _contains_duplicate_atom_map_numbers(_np, _product, mapped);
}

static int
_unmap_duplicate_atom_map_numbers(ISIS_RXN_FILE_Molecule &m, int *mapped) {
  const int matoms = m.natoms();

  const int *amap = m.atom_map();

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    const int mi = amap[i];

    if (mapped[mi] <= 1) {
      continue;
    }

    m.set_atom_map(i, 0);
    m.set_atom_map_number(i, 0);

    // mapped[mi]--;   //tad:  this will cause the last of the duplicates to be NOT
    // removed because the mapped count reduced to 1
    rc++;
  }
  if (rc > 0) {
    m.invalidate_smiles();
  }

  return rc;
}

static int
_unmap_duplicate_atom_map_numbers(const int n, ISIS_RXN_FILE_Molecule *m, int *mapped) {
  int contains_dups = 0;

  for (int i = 0; i < n; ++i) {
    if (_contains_duplicate_atom_map_numbers(m[i], mapped)) {
      contains_dups++;
    }
  }

  if (0 == contains_dups) {
    return 0;
  }

  for (int i = 0; i < n; ++i) {
    _unmap_duplicate_atom_map_numbers(m[i], mapped);
  }

  return 1;
}

/*
  Note that for each duplicated atom map we leave one atom mapped. We do not do this
  in any kind of careful way, could be improved substantially
*/

int
RXN_File::unmap_duplicate_atom_map_numbers() {
  const int h = highest_atom_map_number();

  if (h <= 0) {
    return 1;
  }
  int rc = 0;

  int *mapped = new_int(h + 1);
  std::unique_ptr<int[]> free_mapped(mapped);

  if (_unmap_duplicate_atom_map_numbers(_nr, _reagent, mapped)) {
    rc = 1;
  }

  std::fill_n(mapped, h + 1, 0);

  if (_unmap_duplicate_atom_map_numbers(_np, _product, mapped)) {
    rc = 1;
  }

  return rc;
}

static void
_assign_unmapped_atoms(ISIS_RXN_FILE_Molecule &m, int &highestAtomMapping) {
  const int matoms = m.natoms();

  const int *amap = m.atom_map();

  for (int i = 0; i < matoms; ++i) {
    const int mi = amap[i];

    if (mi > 0) {
      continue;
    }

    highestAtomMapping++;
    m.set_atom_map(i, highestAtomMapping);
    m.set_atom_map_number(i, highestAtomMapping);
  }

  return;
}

static void
_assign_unmapped_atoms(const int n, ISIS_RXN_FILE_Molecule *m, int &highestAtomMapping) {
  for (int i = 0; i < n; ++i) {
    _assign_unmapped_atoms(m[i], highestAtomMapping);
    m->invalidate_smiles();
  }

  return;
}

/*
  Note that for each duplicated atom map we leave one atom mapped. We do not do this
  in any kind of careful way, could be improved substantially
*/

void
RXN_File::assign_unmapped_atoms() {
  int h = highest_atom_map_number();

  if (h <= 0) {
    return;
  }

  _assign_unmapped_atoms(_nr, _reagent, h);

  _assign_unmapped_atoms(_np, _product, h);
}

static int
identify_largest_fragment(const int n, const ISIS_RXN_FILE_Molecule *m) {
  int zbest = 0;
  int most_atoms = m[0].natoms();

  for (int i = 1; i < n; ++i) {
    if (m[i].natoms() < most_atoms) {
      continue;
    }

    zbest = i;
    most_atoms = m[i].natoms();
  }

  return zbest;
}

int
RXN_File::largest_fragment_is_changed() const {
  if (0 == _nr || 0 == _np) {
    return 0;
  }

  int largest_reagent = identify_largest_fragment(_nr, _reagent);
  int largest_product = identify_largest_fragment(_np, _product);

  if (_reagent[largest_reagent].natoms() !=
      _product[largest_product].natoms()) {  // definitely changes
    return 1;
  }

  if (_reagent[largest_reagent].quick_atom_hash() !=
      _product[largest_product].quick_atom_hash()) {
    return 1;
  }

  _reagent[largest_reagent].invalidate_smiles();
  _product[largest_product].invalidate_smiles();

  set_include_atom_map_with_smiles(0);

  _reagent[largest_reagent].recompute_unique_smiles();
  _product[largest_product].recompute_unique_smiles();
  auto ruSmi = _reagent[largest_reagent].unique_smiles();
  auto puSmi = _product[largest_product].unique_smiles();

  const auto rc = (ruSmi != puSmi);

  set_include_atom_map_with_smiles(1);

  _reagent[largest_reagent].invalidate_smiles();
  _product[largest_product].invalidate_smiles();

  return rc;
}

// #define DEBUG_REMOVE_UNCHANGING_FRAGMENTS
/*
  Consider:

  [IH:1].[NH:15]([NH2:16])[C:4]1=[N:10][CH2:9][CH2:8][CH2:7][CH2:6][NH:5]1>>[IH:1].CS[C:4]1=[N:10][CH2:9][CH2:8][CH2:7][CH2:6][NH:5]1.[NH2:15][NH2:16]

  It has a fragment that does not change. It should be removed.

  Note that if a whole multi-fragment reagent/product is the same, then what gets left is
  arbitrary - since we do not remove the last fragment

  Note that this is not handled. The reason is that on one side, the common fragment is a
  fragment, but on the other side it is a reagent.
  [CH3:16][S:17](=[O:18])(=[O:19])[OH:31].[Cl:1][C:2]1=[CH:3][CH:4]=[C:5]2[N:6]([CH2:7][CH2:8][CH2:9][OH:10])[C:11](=[O:12])[NH:13][C:14]2=[CH:15]1>[CH3:20][CH2:21][N:22]([CH2:23][CH3:24])[CH2:25][CH3:26]+[Cl:28][CH2:29][Cl:30]>[Cl:1][C:2]1=[CH:3][CH:4]=[C:5]2[N:6]([CH2:7][CH2:8][CH2:9][OH:10])[C:11](=[O:12])[NH:13][C:14]2=[CH:15]1+[CH3:16][S:17](=[O:18])(=[O:19])[Cl:27]
  04160836__dup_42167 Maybe I need to write a different function to handle those cases....
*/

int
RXN_File::remove_unchanging_fragments() {
  if (nullptr == _product_locator) {
    prepare_for_reaction_construction();
  }

  int rc = 0;

  for (int i = 0; i < _nr; ++i) {
    ISIS_RXN_FILE_Molecule &r = _reagent[i];

    if (1 == r.number_fragments()) {
      continue;
    }

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
    cerr << "RXN_File::remove_unchanging_fragments:reagent " << i << " had "
         << r.number_fragments() << " fragments\n";
    cerr << r.smiles() << endl;
#endif

    for (int j = 0; j < _np; ++j) {
      ISIS_RXN_FILE_Molecule &p = _product[j];

      if (1 == p.number_fragments()) {
        continue;
      }

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
      cerr << "RXN_File::remove_unchanging_fragments:product " << j << " has "
           << p.number_fragments() << " fragments\n";
#endif

      rc += _remove_common_fragments(r, p);
    }
  }

  return rc;
}

int
RXN_File::remove_unchanging_components() {
  if (nullptr == _product_locator) {
    prepare_for_reaction_construction();
  }

  int rc = 0;

  set_include_chiral_info_in_smiles(0);
  for (int i = 0; i < _nr; ++i) {
    _reagent[i].invalidate_smiles();
    _reagent[i].recompute_unique_smiles();
  }
  for (int i = 0; i < _np; ++i) {
    _product[i].invalidate_smiles();
    _product[i].recompute_unique_smiles();
  }
  set_include_chiral_info_in_smiles(1);

  for (int i = 0; i < _nr; ++i) {
    const IWString &iusmi = _reagent[i].unique_smiles();

    //  cerr << "Reagent " << i << " usmi " << iusmi << endl;

    for (int j = 0; j < _np; ++j) {
      if (_product[j].unique_smiles() != iusmi) {  // definitely not the same
        continue;
      }

      //    cerr << "in " << name() << " duplicate reagents i = " << i << " j = " << j <<
      //    endl;

      for (int k = i; k < (_nr - 1);
           ++k)  // reagent I and product j are the same.  Remove them
      {
        _reagent[k] = std::move(_reagent[k + 1]);
      }
      for (int k = j; k < (_np - 1);
           ++k)  // reagent I and product j are the same.  Remove them
      {
        _product[k] = std::move(_product[k + 1]);
      }

      i--;
      _nr--;
      j--;
      _np--;
      rc++;
      break;
    }
  }

  return rc;
}

// Remove fragments where the mapped atoms are the same in
// reagents and products. The idea of 'same' is defined in
// the functions called below.
int
RXN_File::_remove_common_fragments(ISIS_RXN_FILE_Molecule &r, ISIS_RXN_FILE_Molecule &p) {
  const int rf = r.number_fragments();
  const int pf = p.number_fragments();
  r.compute_aromaticity_if_needed();
  p.compute_aromaticity_if_needed();

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
  cerr << "RXN_File::_remove_common_fragments:remving common frags btw " << r.smiles()
       << " and " << p.smiles() << endl;
#endif

  resizable_array<int> reagent_fragment_to_remove, product_fragment_to_remove;

  for (int i = 0; i < rf; ++i) {
    for (int j = 0; j < pf; ++j) {
      if (!_same_atoms_and_bonds(r, i, p, j)) {
        continue;
      }

      reagent_fragment_to_remove.add(i);
      product_fragment_to_remove.add(j);
    }
  }

  if (reagent_fragment_to_remove.empty()) {
    return 0;
  }

  r.delete_fragments(reagent_fragment_to_remove);
  p.delete_fragments(product_fragment_to_remove);

  return reagent_fragment_to_remove.number_elements();
}

// Do the mapped atoms in fragment 'rf' in 'r' look very much like
// the mapped atoms in fragment 'pf' in 'p'.
int
RXN_File::_same_atoms_and_bonds(ISIS_RXN_FILE_Molecule &r, const int rf,
                                ISIS_RXN_FILE_Molecule &p, const int pf) const {
  const int matoms = r.natoms();

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
  cerr << "_same_atoms_and_bonds:reagent fragment " << rf << " product fragment " << pf
       << endl;
  cerr << "atoms_in_fragment " << r.atoms_in_fragment(rf) << " and "
       << p.atoms_in_fragment(pf) << endl;
#endif

  for (int i = 0; i < matoms; ++i) {
    if (r.fragment_membership(i) != rf) {
      continue;
    }

    // Must have matching atom maps.
    const int amap = r.atom_map(i);

    if (amap <= 0) {
      return 0;
    }

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
    cerr << " atom " << i << " in reagent fragment " << r.smarts_equivalent_for_atom(i)
         << " map " << amap << endl;
#endif

    const int x = p.which_is_mapped_atom(amap);
    if (x < 0) {
      return 0;
    }

#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
    cerr << " atom " << x << " in product fragment " << p.smarts_equivalent_for_atom(x)
         << endl;
#endif

    if (p.fragment_membership(x) != pf) {
      return 0;
    }

    if (!_atoms_the_same(r, i, p, x)) {
      return 0;
    }
  }

  return 1;
}

static int
bond_numeric_constant(const Bond *b) {
  if (b->is_aromatic()) {
    return 4;
  }

  if (b->is_single_bond()) {
    return 1;
  }

  if (b->is_double_bond()) {
    return 2;
  }

  if (b->is_triple_bond()) {
    return 3;
  }

  return 9;
}

int
RXN_File::_atoms_the_same(ISIS_RXN_FILE_Molecule &m1, const atom_number_t a1,
                          ISIS_RXN_FILE_Molecule &m2, const atom_number_t a2) const {
  const Atom *aa1 = m1.atomi(a1);
  const Atom *aa2 = m2.atomi(a2);
#ifdef DEBUG_REMOVE_UNCHANGING_FRAGMENTS
  cerr << "_atoms_the_same comparing atom " << a1 << " "
       << m1.smarts_equivalent_for_atom(a1) << " with " << a2 << " "
       << m2.smarts_equivalent_for_atom(a2) << endl;
#endif

  if (aa1->element() != aa2->element()) {
    return 0;
  }

  if (aa1->ncon() != aa2->ncon()) {
    return 0;
  }

  if (m1.hcount(a1) != m2.hcount(a2)) {
    return 0;
  }

  if (aa1->nbonds() != aa2->nbonds()) {
    return 0;
  }

  if (aa1->formal_charge() != aa2->formal_charge()) {
    return 0;
  }

  if (m1.ring_bond_count(a1) != m2.ring_bond_count(a2)) {
    return 0;
  }

  if (m1.is_aromatic(a1) != m2.is_aromatic(a2)) {
    return 0;
  }

  resizable_array<int> c1;

  const int map1 = m1.atom_map(a1);
  const int map2 = m2.atom_map(a2);

  for (const Bond *b : *aa1) {
    const atom_number_t j = b->other(a1);

    const int x = m1.atom_map(j);

    if (x <= 0) {
      continue;
    }

    c1.add(map1 * bond_numeric_constant(b) * x);
  }

  for (const Bond *b : *aa2) {
    const atom_number_t j = b->other(a2);

    const int x = m2.atom_map(j);
    if (x <= 0) {
      continue;
    }

    const auto y = map2 * bond_numeric_constant(b) * x;

    if (!c1.contains(y)) {
      return 0;
    }
  }

  return 1;
}

int
RXN_File::remove_cis_trans_bonding() {
  int rc = 0;

  for (int i = 0; i < _nr; ++i) {
    rc += _reagent[i].revert_all_directional_bonds_to_non_directional();
  }
  for (int i = 0; i < _np; ++i) {
    rc += _product[i].revert_all_directional_bonds_to_non_directional();
  }

  return rc;
}

Product_Atom_Types::Product_Atom_Types() {
  _pstart = nullptr;
  _atom_type = nullptr;
}

Product_Atom_Types::~Product_Atom_Types() {
  if (nullptr != _atom_type) {
    delete[] _atom_type;
  }

  if (nullptr != _pstart) {
    delete[] _pstart;
  }
}

int
Product_Atom_Types::initialise(RXN_File &rxn, Atom_Typing_Specification &ats) {
  const int np = rxn.number_products();
  if (0 == np) {  // hard to imagine
    return 0;
  }

  if (nullptr != _atom_type) {
    delete[] _atom_type;
  }

  _pstart = new int[np];

  _pstart[0] = 0;

  int atoms_in_products = 0;

  for (int i = 0; i < np; ++i) {
    const auto x = rxn.product(i).natoms();

    if (i != (np - 1)) {
      _pstart[i + 1] = _pstart[i] + x;
    }

    atoms_in_products += x;
  }

  if (nullptr != _atom_type) {
    delete[] _atom_type;
  }

  _atom_type = new int[atoms_in_products];

  for (int i = 0; i < np; ++i) {
    ats.assign_atom_types(rxn.product(i), _atom_type + _pstart[i]);
  }

  return 1;
}

int
Product_Atom_Types::atom_type(const int f, const atom_number_t a) const {
  return _atom_type[_pstart[f] + a];  // no checking here!
}
