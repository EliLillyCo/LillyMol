/*
  Make specific molecules from a combinatorially derived set
*/

#include <stdlib.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/make_these_molecules.pb.h"

using std::cerr;

const char *prog_name = nullptr;

/*
  A set of reagents is just a set of molecules that can be retrieved by name
  Note that the Molecule_and_Embedding objects are never deleted
*/

typedef IW_STL_Hash_Map<IWString, Molecule_and_Embedding *> ID_to_Molecule;

class Set_of_Reagents : public ID_to_Molecule {
 private:
  IWString _fname;

 public:
  ~Set_of_Reagents();

  const IWString &fname() const {
    return _fname;
  }

  int read_reagents(const char *, FileType);
  int read_reagents(data_source_and_type<Molecule_and_Embedding> &);

  int number_reagents() const {
    return size();
  }

  Molecule_and_Embedding *operator[](const IWString &) const;

  int ensure_sidechains_match_reaction_queries(Reaction_Site &q, const IWString &rname,
                                               int take_first_of_multiple_matches);
};

Set_of_Reagents::~Set_of_Reagents() {
  for (auto i : *this) {
    delete i.second;
  }

  return;
}

int
Set_of_Reagents::read_reagents(const char *fname, FileType input_type) {
  data_source_and_type<Molecule_and_Embedding> input(input_type, fname);

  if (!input.good()) {
    cerr << "Set_of_Reagents::read_reagents:cannot open '" << fname << "'\n";
    return 0;
  }

  _fname = fname;

  return read_reagents(input);
}

int
Set_of_Reagents::read_reagents(data_source_and_type<Molecule_and_Embedding> &input) {
  Molecule_and_Embedding *m;

  while (nullptr != (m = input.next_molecule())) {
    //  _preprocess(*m);
    const IWString &mname = m->name();
    if (!mname.contains(' ')) {
      ID_to_Molecule::insert(
          value_type(m->name(), m));  // no checking for duplicate names!
    } else {
      IWString tmp(mname);
      tmp.truncate_at_first(' ');
      ID_to_Molecule::insert(value_type(tmp, m));  // no checking for duplicate names!
    }
  }

  return size();
}

Molecule_and_Embedding *
Set_of_Reagents::operator[](const IWString &s) const {
  const auto f = ID_to_Molecule::find(s);

  if (f == ID_to_Molecule::end()) {
    return nullptr;
  }

  return (*f).second;
}

/*
  We want to do the substructure search once and store the result
*/

int
Set_of_Reagents::ensure_sidechains_match_reaction_queries(
    Reaction_Site &q, const IWString &rname, int take_first_of_multiple_matches) {
  for (auto i : *this) {
    Molecule_and_Embedding *m = i.second;

    Substructure_Results sresults;

    const int nhits = q.substructure_search(m, sresults);

    if (1 == nhits) {  // great
      m->set_embedding(*(sresults.embedding(0)));
    } else if (0 == nhits) {
      cerr << "No hits to reaction '" << rname << "' in reagent " << m->smiles() << ' '
           << m->name() << "', only matched "
           << sresults.max_query_atoms_matched_in_search() << " atoms\n";
      cerr << m->smiles() << '\n';
      return 0;
    } else if (take_first_of_multiple_matches) {
      m->set_embedding(*(sresults.embedding(0)));
    } else {
      cerr << nhits << " hits to '" << rname << "' in reagent '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

class Make_These_Molecules {
 private:
  int _verbose;

  Chemical_Standardisation _chemical_standardisation;
  int _reduce_reagents_to_largest_fragment;
  int _reduce_products_to_largest_fragment;

  Sidechain_Match_Conditions _smc;

  int _take_first_of_multiple_matches;

  int _number_sidechains;

  /*!
   * \var int _total_sidechain
   * \brief Total sidechains listed in all the reactions
   */
  int _total_sidechain;

  /*!
   * \var int _current_sidechain_id
   * \brief Index for the current sidechain id from all reactions
   */
  int _current_sidechain_id;

  /*!
   * \var Sidechain_Reaction_Site * _current_backup_sidechain
   * \brief Pointer for the original sidechains info from reaction file.
   *        This pointer is based on the _current_sidechain_id.
   *        It always points to the original sidechain coresponding to
   * _current_sidechain_id.
   */
  Sidechain_Reaction_Site *_current_backup_sidechain;

  /*!
   * \var IWReaction * _reaction
   * \brief Store the reaction data.
   *        This reaction data maybe modified for perform reaction
   *        This reaction data can be recovered using _reaction_backup data
   */
  IWReaction *_reaction;

  /*!
   * \var IWReaction * _reaction_backup
   * \brief Store the original reaction data
   */
  IWReaction *_reaction_backup;

  int _nr;

  Set_of_Reagents *_reagent;

  int _convert_isotopes_in_product_molecules;

  int _molecules_written;

  int _return_if_no_match;

  int _single_reaction;

  IWString _component_separator;

  /*!
   * \var _input_is_proto
   * \brief Input is a MakeTheseMoleculesData::Product proto.
   *
   */
  int _input_is_proto;

  // private functions

  void _default_values();
  void _usage(int rc);
  void _preprocess(Molecule &m);
  /*!
   * \fn int _ensure_sidechains_match_reaction_queries()
   * \brief Make sure all molecules in all sidechains match their query conditions.
   *        It also captures the embeddings that will be used later
   */
  int _ensure_sidechains_match_reaction_queries();

  /*!
   * \fn int _make_these_molecules (const std::vector<Molecule_and_Embedding *> &
   * reagents, Molecule_Output_Object & output)
   * \brief Perform reaction on each sidechain.
   *        Generate name for each new molecule
   *        Reaction site shall be prepared for the reaction before calling this function
   * \param reagents  The reagent list used for generate new molecule name
   * \param name_stem  Prepended to the molecule name.
   * \param output The new molecule
   */
  int _make_these_molecules(const std::vector<Molecule_and_Embedding *> &reagents,
                            const IWString& name_stem,
                            Molecule_Output_Object &output);

  /*!
   * \fn int _make_these_molecules (const const_IWSubstring & buffer,
   *                                Molecule_Output_Object & output)
   * \brief Prepare each sidechain before reaction
   * \var buffer List of the molecule ID for reaction
   * \var output The new molecule
   */
  int _make_these_molecules(const const_IWSubstring &buffer,
                            Molecule_Output_Object &output);

  int _make_these_molecules(iwstring_data_source &input, Molecule_Output_Object &output);

  int _make_these_molecules(const char *fname, Molecule_Output_Object &output);
  /*!
   * \fn Sidechain_Reaction_Site * _get_current_sidechain(void)
   * \brief Return the pointer for current sidechain from all reactions.
   *        The pointer for the backup sidechain is updated
   *        It will return NULL pointer at the end of all sidechain
   */
  Sidechain_Reaction_Site* _get_current_sidechain();

 public:
  Make_These_Molecules();
  ~Make_These_Molecules();

  int operator()(int argc, char** argv);
};

Make_These_Molecules::Make_These_Molecules() {
  _default_values();

  return;
}

Make_These_Molecules::~Make_These_Molecules() {
  if (nullptr != _reaction) {
    delete[] _reaction;
  }
  if (nullptr != _reagent) {
    delete[] _reagent;
  }
  if (nullptr != _reaction_backup) {
    delete[] _reaction_backup;
  }

  return;
}

void
Make_These_Molecules::_default_values() {
  _verbose = 0;
  _reduce_reagents_to_largest_fragment = 0;
  _reduce_products_to_largest_fragment = 0;
  _number_sidechains = 0;
  _take_first_of_multiple_matches = 0;
  _convert_isotopes_in_product_molecules = 0;
  _molecules_written = 0;
  _return_if_no_match = 0;
  _single_reaction = 0;
  _current_sidechain_id = 0;
  _reaction_backup = nullptr;
  _component_separator = " + ";
  _input_is_proto = 0;

  _reaction = nullptr;

  _nr = 0;
  _reagent = nullptr;

  return;
}

void
Make_These_Molecules::_usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Takes any number of sidechains (N) and (N-1) reactions\n";
  cerr << "  -R <rxn>      one or more reactions descibing how to assemble molecules\n";
  cerr << "  -p            reactions are proto form\n";
  cerr << "  -M <fname>    file containing products to be made\n";
  cerr << "  -S <stem>     output stream\n";
  cerr << "  -z f          take first of any multiple matches in sidechains\n";
  cerr << "  -z i          ignore no match errors in sidechains\n";
  cerr << "  -I            change isotopes to natural form in product molecules\n";
  cerr << "  -s            single reaction, many sidechains\n";
  cerr << "  -W <string>   token put between names of products (default \" + \")\n";
  cerr << "  -Y ...        miscellaneous options, enter '-M help' for info\n";
  cerr << "  -l            reduce reagents to largest fragment\n";
  cerr << "  -L            reduce products to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

void
Make_These_Molecules::_preprocess(Molecule &m) {
  if (_reduce_reagents_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return;
}

/*
  Here is the guts of things.

  We have changed the input record into a set of Molecule_and_Embedding objects.

  First task is to assign the reagents identified to the appropriate reaction object
*/

int
Make_These_Molecules::_make_these_molecules(
    const std::vector<Molecule_and_Embedding *> &reagents,
    const IWString& name_stem,
    Molecule_Output_Object &output) {
  // Product of all reactions
  Molecule product;
  // Product name for new molecule
  IWString mname(name_stem);

  // Perform the first reaction. We should have at least one reaction
  int current_reagent_id = 0;

  // Start with scaffold id for the new molecule name
  mname = reagents[current_reagent_id]->name();
  // Run the reaction for the first sidechain with the scaffold
  if (!_reaction[0].perform_reaction(reagents[current_reagent_id],
                                     reagents[current_reagent_id]->embedding(),
                                     product)) {
    cerr << "initial reaction failed!" << reagents[0]->name() << '\n';
    return _return_if_no_match;
  }

  current_reagent_id++;
  if (1 == _single_reaction) {
    // Single reaction
    // All the reaction is done and for single reaction
    // Each reagent id will be appended to the new molecule name
    for (unsigned int ui = 1; ui < reagents.size(); ++ui) {
      mname << _component_separator << reagents[ui]->name();
    }
  } else {
    // Multiple reactions
    // Append molecule name for the first reaction
    mname << _component_separator << reagents[current_reagent_id]->name();
    current_reagent_id++;
    // Run rest of reactions 1 ... _total_sidechain -1
    for (int i = 1; i < _total_sidechain; i++) {
      Substructure_Results sresults;
      int nhits = _reaction[i].substructure_search(product, sresults);

      if (0 == nhits) {
        cerr << "no hits to reaction " << i << ", " << _reaction[i].comment()
             << "for intermediate product " << mname << ' ' << product.smiles() << '\n';
        return _return_if_no_match;
      }

      // Make the reaction
      Molecule new_Product;
      _reaction[i].perform_reaction(&product, sresults.embedding(0), new_Product);
      product = new_Product;
      mname << _component_separator << reagents[current_reagent_id]->name();
      current_reagent_id++;
    }
  }

  product.set_name(mname);

  if (_convert_isotopes_in_product_molecules) {
    product.transform_to_non_isotopic_form();
  }
  if (_reduce_products_to_largest_fragment) {
    product.reduce_to_largest_fragment_carefully();
  }

  _molecules_written++;

  return output.write(product);
}

Sidechain_Reaction_Site *
Make_These_Molecules::_get_current_sidechain(void) {
  Sidechain_Reaction_Site *ptr = nullptr;
  if (_current_sidechain_id < _total_sidechain) {
    if (1 == _single_reaction) {
      ptr = _reaction[0].sidechain(_current_sidechain_id);
      _current_backup_sidechain = _reaction_backup[0].sidechain(_current_sidechain_id);
    } else {
      ptr = _reaction[_current_sidechain_id].sidechain(0);
      _current_backup_sidechain = _reaction_backup[_current_sidechain_id].sidechain(0);
    }
    _current_sidechain_id++;
  } else {
    // Reset _current_sidechain_id to 0
    _current_sidechain_id = 0;
  }
  return ptr;
}

// populate `reagent_names` and `name_stem` from values in `proto`.
int
ReadReagents(const const_IWSubstring& buffer,
             resizable_array_p<IWString>& reagent_names,
             IWString& name_stem) {
  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());
  MakeTheseMoleculesData::Product proto;
  if (! google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "ReadReagents:cannot parse '" << buffer << "'\n";
    return 0;
  }

  name_stem = proto.name_stem();

  for (const std::string& r : proto.reagent()) {
    reagent_names << new IWString(r);
  }

  return reagent_names.size();
}

int
Make_These_Molecules::_make_these_molecules(const const_IWSubstring &buffer,
                                            Molecule_Output_Object &output) {
  resizable_array_p<IWString> reagent_names;
  IWString name_stem;

  if (_input_is_proto) {
    if (! ReadReagents(buffer, reagent_names, name_stem)) {
      cerr << "Cannot read proto '" << buffer << "'\n";
      return 0;
    }
  } else {
    const_IWSubstring token;
    int i = 0;
    while (buffer.nextword(token, i)) {
      reagent_names << new IWString(token);
    }
  }

  // Check number for id in the target ID file with number of provided reagent files and
  // embedded smile in reaction file The id count shall equal to _number_sidechains
  if (reagent_names.number_elements() != _number_sidechains) {
    cerr << "Incorrect number reagents specified, found " << reagent_names.number_elements()
         << " expected " << _number_sidechains << '\n';
    return 0;
  }

  // Check the scaffold id
  int current_reagent_id = 0;  // An index into the reagent_names vector.
  const IWString& token = *reagent_names[current_reagent_id];

  Molecule_and_Embedding *m = _reagent[current_reagent_id][token];
  if (m == nullptr) {
    cerr << "Cannot find molecule for '" << token << "', in scaffold\n";
    return 0;
  }

  // Populate the reagent list for reaction. Convert names to Molecule_and_Embedding
  std::vector<Molecule_and_Embedding *> reagents;

  reagents.push_back(m);
  ++current_reagent_id;

  // Check all reagent id
  _current_sidechain_id = 0;
  Sidechain_Reaction_Site *site_ptr = _get_current_sidechain();
  while (site_ptr != nullptr) {
    const IWString& token = *reagent_names[current_reagent_id];

    // Use backup sidechain data as footprint
    int reagent_count = _current_backup_sidechain->number_reagents();
    if (0 == reagent_count) {
      // Smart in sidechain
      const auto iter = _reagent[current_reagent_id].find(token);
      if (iter == _reagent[current_reagent_id].end()) {
        cerr << "Cannot find molecule for '" << token << "', reagent "
             << current_reagent_id << "\n";
        cerr << "Hash contains " << _reagent[current_reagent_id].size() << " items\n";
        for (const auto& [k, v] : _reagent[current_reagent_id]) {
          cerr << " key " << k << " value " << v->name() << '\n';
        }
        return 0;
      }
      Molecule_and_Embedding* m = iter->second;
      // Populate the reagent list for reaction
      reagents.push_back(m);
      if (site_ptr->number_reagents() > 0) {
        site_ptr->remove_first_reagent_no_delete();
      }
      // Regenerate reagents with match id. TODO: maybe it is here that we get
      // the double memory free problem if a reagent is not found. 
      site_ptr->add_reagent_embedding_identified(m, _smc);
      // Move to next reagent file
      current_reagent_id++;
    } else if (1 == reagent_count) {
      // Smile in the sidechain
      Molecule_and_Embedding *m = _current_backup_sidechain->reagent(0);
      reagents.push_back(m);
    } else {
      // This code shall never be reached. This is for failsafe only
      cerr << "Unexpected smiles in the the sidechain\n";
      return 0;
    }
    site_ptr = _get_current_sidechain();
  }

  return _make_these_molecules(reagents, name_stem, output);
}

//  The input stream is the file containing the records with individual
//  reagent combinations
int
Make_These_Molecules::_make_these_molecules(iwstring_data_source &input,
                                            Molecule_Output_Object &output) {
  input.set_translate_tabs(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {  // skip comments
      continue;
    }

    if (!_make_these_molecules(buffer, output)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Make_These_Molecules::_ensure_sidechains_match_reaction_queries() {
  // Add this to count how many reagent file is required based on user input
  int required_reagent_file = 0;
  // If there is multiple reaction with single sidechain
  _total_sidechain = _nr;
  // Ensure the number of reagent matches total sidechain in all reactions
  for (auto i = 0; i < _nr; i++) {
    // Read in number of side chains in the reaction file
    int number_of_sidechain = _reaction[i].number_sidechains();

    // Check the side chain with the embedded smile
    int sidechain_with_reagent = _reaction[i].number_sidechains_with_reagents();

    if (0 == _single_reaction) {
      // Assumption: one side chain for each reaction file
      // Only allow one side chain for each reaction file if multiple reaction files are
      // provided
      if (1 != number_of_sidechain) {
        // Require at least for side chain for any reaction file
        cerr << "Warning, can only process reactions with a single sidechain, look for "
                "unpredictable results\n";
        // TODO check the standard return error code
        return 0;
      } else {
        // We do not need the reagent file if the side chain has embedded smile
        required_reagent_file += (number_of_sidechain - sidechain_with_reagent);
      }
    } else {
      // One reaction file with multiple side chains
      required_reagent_file = number_of_sidechain - sidechain_with_reagent;
      _total_sidechain = number_of_sidechain;
      break;
      // This code shall only be executed once because the _nr is 1 for this codition
    }
  }
  // Check if there is enough reagent file from input
  // _number_sidechain includes one scaffold file and all reagent file
  if (required_reagent_file + 1 != _number_sidechains) {
    cerr << required_reagent_file << " reagent files are required , only  "
         << (_number_sidechains - 1) << " reagent files are found in the input.\n";
    return 0;
  }

  if (_total_sidechain > (_number_sidechains - 1)) {
    cerr << "Generated molecules will be complemented with smiles provided in the "
            "reaction\n";
  }

  Reaction_Site &s = _reaction[0];

  if (!_reagent[0].ensure_sidechains_match_reaction_queries(
          s, _reaction[0].comment(), _take_first_of_multiple_matches)) {
    cerr << "Non matching reagent, reaction 0 in sidechain 0, file '"
         << _reagent[0].fname() << "'\n";
    return 0;
  }

  if (1 == _single_reaction) {
    // Only one reaction file for this case
    // The first reagent file is scaffold. Start with the second file here
    int reagent_index = 1;
    for (int i = 0; i < _reaction[0].number_sidechains(); ++i) {
      auto &s = *(_reaction[0].sidechain(i));
      if (0 == s.number_reagents()) {
        // Only check reaction match if now reagent in the side chain
        if (!_reagent[reagent_index].ensure_sidechains_match_reaction_queries(
                s, _reaction[0].comment(), _take_first_of_multiple_matches)) {
          cerr << "Non matching reagent, reaction 0 in sidechain " << i << ", file '"
               << _reagent[0].fname() << "'\n";
          return 0;
        } else {
          reagent_index++;
        }
      }
    }
  } else {
    // The first reagent file is scaffold. Start with the second file here
    int reagent_index = 1;
    for (auto i = 0; i < _nr; i++)  // loop over reactions
    {
      Sidechain_Reaction_Site *s = _reaction[i].sidechain(0);

      if (0 == s->number_reagents()) {
        if (!_reagent[reagent_index].ensure_sidechains_match_reaction_queries(
                *s, _reaction[i].comment(), _take_first_of_multiple_matches)) {
          cerr << "Non matching reagent, reaction " << i << " in sidechain " << (i + 1)
               << ", file '" << _reagent[i + 1].fname() << "'\n";
          return 0;
        } else {
          reagent_index++;
        }
      }
    }
  }

  return 1;
}

int
Make_These_Molecules::operator()(int argc, char **argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lLR:pS:z:M:IsW:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    _usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, _verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      _usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, _verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l')) {
    _reduce_reagents_to_largest_fragment = 1;

    if (_verbose) {
      cerr << "Will reduce reagents to largest fragment\n";
    }
  }

  if (cl.option_present('l')) {
    _reduce_products_to_largest_fragment = 1;

    if (_verbose) {
      cerr << "Will reduce to products largest fragment\n";
    }
  }

  if (!cl.option_present('R')) {
    cerr << "Must specify one or more reactions via the -R option\n";
    _usage(2);
  }

  // Reading in number of provided reaction file from arguments
  _nr = cl.option_count('R');
  if (cl.option_present('s')) {
    // This is the case for only one reaction file with multiple sidechain
    _single_reaction = 1;
    if (1 != _nr) {
      cerr << "Only one reaction file is allowed for option -s\n";
      _usage(2);
    }
    if (_verbose) {
      cerr << "Input interpreted as single reaction with multiple sidechains\n";
    }
  }

  int reactions_are_proto_form = cl.option_present('p');

  // Number of provided side chain files
  _number_sidechains = cl.number_elements();
  _reaction = new IWReaction[_nr];
  _reaction_backup = new IWReaction[_nr];
  Sidechain_Match_Conditions smc;
  // Populate the reaction data
  // TODO: Is it possible to just copy _reaction_backup from _reaction?
  // Unclear, so for now, read both from the file.
  IWString fname;
  for (int i = 0; cl.value('R', fname, i); i++) {
    int ok;
    if (reactions_are_proto_form) {
      ok = _reaction[i].Read(fname);
      _reaction_backup[i].Read(fname);
    } else {
      ok = _reaction[i].do_read(fname, smc);
      _reaction_backup[i].do_read(fname, smc);
    }

    if (!ok) {
      cerr << "Cannot build reaction from '" << fname << "'\n";
      return 2;
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (auto i = 0; cl.value('z', z, i); i++) {
      if ('f' == z) {
        _take_first_of_multiple_matches = 1;
        if (_verbose) {
          cerr << "Will take the first of any multiple matches in sidechains\n";
        }
      } else if ('i' == z) {
        _return_if_no_match = 1;
        if (_verbose) {
          cerr << "Will ignore if no matches in sidechains\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        _usage(2);
      }
    }
  }

  if (cl.option_present('W')) {
    _component_separator = cl.string_value('W');

    set_component_separator(_component_separator);

    if (_verbose) {
      cerr << "Component separator set to '" << _component_separator << "'\n";
    }
  }

  if (!cl.option_present('M')) {
    cerr << "Must specify molecules to make via the -M option\n";
    _usage(2);
  }

  iwstring_data_source input;

  if (cl.option_present('M')) {
    const char *m = cl.option_value('M');

    if (!input.open(m)) {
      cerr << "Cannot open file of molecules to be made '" << m << "'\n";
      return 2;
    }
  }

  if (cl.option_present('I')) {
    _convert_isotopes_in_product_molecules = 1;

    if (_verbose) {
      cerr << "Isotopes in product molecules will be converted\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      _usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  _reagent = new Set_of_Reagents[_number_sidechains];

  for (auto i = 0; i < _number_sidechains; i++) {
    if (!_reagent[i].read_reagents(cl[i], input_type)) {
      cerr << "Cannot assemble reagents '" << cl[i] << "'\n";
      return i + 1;
    }

    if (_verbose) {
      cerr << "Read " << _reagent[i].number_reagents() << " reagents from '" << cl[i]
           << "'\n";
    }
  }

  // The reactions must react with the sidechains we have read
  if (0 == _ensure_sidechains_match_reaction_queries()) {
    cerr << "One or more sidechains do match reaction query conditions\n";
    // return 3;
    // Check the error code
    _usage(3);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  Molecule_Output_Object output;

  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      _usage(2);
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s)) {
      cerr << "Cannot overwrite input file(s), stem '" << s << "'\n";
      return 2;
    }

    if (!output.new_stem(s)) {
      cerr << "Cannot initialise output stream '" << s << "'\n";
      return 2;
    }
  } else {
    if (!output.new_stem("-")) {
      cerr << "huh, cannot initialise stdout\n";
      return 1;
    }
  }

  if (!_make_these_molecules(input, output)) {
    cerr << "Processing failed\n";
    return 2;
  }

  if (_verbose) {
    cerr << "Wrote " << _molecules_written << " molecules\n";
  }

  // Before the descructor gets called, we need to make sure none of the reactions think
  // they own any of the reagents

  if (_single_reaction) {
    for (int i = 1; i <= _reaction[0].number_sidechains(); ++i) {
      Sidechain_Reaction_Site *s = _reaction[0].sidechain(i - 1);

      if (s->number_reagents()) {  // if already got a reagent, remove it, but do not
                                   // delete it
        s->remove_first_reagent_no_delete();
      }
    }
  } else {
    for (auto i = 0; i < (_number_sidechains - 1); i++)  // loop over reactions
    {
      Sidechain_Reaction_Site *s = _reaction[i].sidechain(0);

      if (s->number_reagents()) {  // if already got a reagent, remove it, but do not
                                   // delete it
        s->remove_first_reagent_no_delete();
      }
    }
  }

  return 0;
}

int
main(int argc, char **argv) {
  prog_name = argv[0];

  Make_These_Molecules Make_These_Molecules;

  return Make_These_Molecules(argc, argv);
}
