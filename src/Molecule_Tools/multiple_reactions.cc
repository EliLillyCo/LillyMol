//  Apply a set of reactions to input molecules.

#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/multiple_reactions.pb.h"

namespace apply_multiple_reactions {

using std::cerr;
namespace fs = std::filesystem;

void
Usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Apply a set of reactions. Reactions are applied separately.
Only in place transformations are applied, so there can be no externally specified sidechains.
 -C <fname>        textproto containing multiple_reactions.Options proto
 -R <rxn>          specify one or more reactions
 -R F:<fname>      file containing multiple reactions
 -R TFDATA:<fname> file containing serialised ReactionProto::Reaction protos
 -G <rxn>          one or more old format reactions
 -r <recursion>    number of recursive invocations (def 0)
 -z i              ignore molecules not matching any queries
 -z w              write molecules not reaction to the output
 -S <fname>        file name stem for output
 -V <fname>        write products with bad valence to <fname>
 -V drop           discard products with bad valence - no writing
 -p                write starting molecule in addition to products
 -b                perform the substructure search whilever there are matches
 -x                apply all reactions to the starting molecule, generating 1 product
 -J ...            other options, enter '-J help' for info
 -c                remove chirality
 -l                strip to largest fragment
 -o <type>         specify output type(s)
 -v                verbose output
  )";

  ::exit(rc);
}

class MultipleReactions {
 private:

  int _verbose;

  int _molecules_read;

  int _reduce_to_largest_fragment;

  int _remove_chirality;

  resizable_array_p<IWReaction> _rxn;
  //  Keep track of how many times each reaction matches
  extending_resizable_array<int> _rxn_matched;

  //  We keep track of how many of the reactions match each molecule.
  extending_resizable_array<int> _reactions_matching;

  // A molecule may generate products if none of the reactions match,
  // or it it only produces duplicate products.
  int _molecules_not_generating_products;

  // By default, reactions are applied separately, generating a new
  // product for each reaction.
  // Alternatively, we can apply all reactions to the molecule, generating
  // just one product.
  int _apply_all_reactions;

  // By default, we use in_place_transformation for all reactions.
  // But if there are reactions that have multiple sidechains, we need
  // to do those reactions via an iterator
  int _perform_reactions_via_iterator;

  //  Used to keep identify duplicates.
  IW_STL_Hash_Set _seen;

  //  The number of duplicates we identify;
  int _duplicates_discarded;

  // If reading sidechains from a proto config file.
  Sidechain_Match_Conditions _smc;

  int _append_reaction_name_to_product;

  // The first match to a reaction might mess up the molecule to the point
  // where matches found are no longer valid.
  int _re_search_each_scaffold_hit;

  Chemical_Standardisation _chemical_standardisation;

  // We can optionally discard molecules with bad valences.
  int _discard_bad_valence;
  IWString_and_File_Descriptor _stream_for_bad_valence;
  int _bad_valence_detected;

  //  private functions.

  // Reading old format reactions. Minimally supported.
  int ReadReactionOldFormat(IWString& fname);
  int ReadReactionOldFormat(iwstring_data_source& input, const IWString& fname);

  // Read proto forms.
  int ReadReaction(IWString& fname);
  int ReadFileOfReactions(IWString& fname);
  int ReadFileOfReactions(const IWString& fname, iwstring_data_source& input);

  int ReadTfDataReactions(IWString& fname);
  int ReadTfDataReactions(iw_tf_data_record::TFDataReader& input, const IWString& dirname);

  int ReadProtoConfig(IWString& fname);
  int ReadReaction(const multiple_reactions::ReactionAndSidechains& proto);

  int Process(Molecule& m, Molecule_to_Match& target, IWReaction& rxn,
              resizable_array_p<Molecule>& result);

  int Process(Molecule& m, const Set_of_Atoms& embedding, IWReaction& rxn,
              resizable_array_p<Molecule>& result);
  int ProcessReSearch(Molecule& m, IWReaction& rxn,
                      resizable_array_p<Molecule>& result);
  int ProcessAllReactionsApplied(Molecule& m, resizable_array_p<Molecule>& result);

  int ViaIterator(Molecule& m, resizable_array_p<Molecule>& results);
  int ViaIterator(Molecule& m,
                        int ndx,
                        Reaction_Iterator* iter,
                        resizable_array_p<Molecule>& results);

  int IsDuplicate(Molecule& m);
  void MaybeAppendRxnName(const IWString& rxnname, Molecule& product);

 public:

  MultipleReactions();

  int Initialise(Command_Line& cl);

  int Preprocess(Molecule& m);

  int Process(Molecule& m, resizable_array_p<Molecule>& result);

  int Report(std::ostream& output) const;
};

MultipleReactions::MultipleReactions()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _duplicates_discarded = 0;
  _molecules_not_generating_products = 0;
  _append_reaction_name_to_product = 0;
  _perform_reactions_via_iterator = 0;
  _apply_all_reactions = 0;
  _discard_bad_valence = 0;
  _bad_valence_detected = 0;
  _re_search_each_scaffold_hit = 0;
}

void
DisplayDashJQualifiers(std::ostream& output) {
  output << R"(The following options are recognisde
 -J noschmsg    do NOT write warning messages about no sidechain substructure matches
  )";

  ::exit(0);
}

int
MultipleReactions::Initialise(Command_Line& cl)
{
  set_copy_name_in_molecule_copy_constructor(1);

  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('R')) {
  } else if (cl.option_present('G')) {
  } else if (cl.option_present('C')) {
  } else {
    cerr << "Must specify one or more reactions via the -R option\n";
    return 0;
  }

  if (cl.option_present('R')) {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (fname.starts_with("F:")) {
        fname.remove_leading_chars(2);
        if (!ReadFileOfReactions(fname)) {
          cerr << "Cannot read reactions from '" << fname << "'\n";
          return 0;
        }
      }
      else if (fname.starts_with("TFDATA:")) {
        fname.remove_leading_chars(7);
        if (! ReadTfDataReactions(fname)) {
          cerr << "Cannot read TFDataRecord file '" << fname << "'\n";
          return 0;
        }
      }
      else {
        if (!ReadReaction(fname)) {
          cerr << "Cannot read reaction '" << fname << "'\n";
          return 0;
        }
      }
    }
  }

  if (cl.option_present('G')) {
    IWString fname;
    for (int i = 0; cl.value('G', fname, i); ++i) {
      if (! ReadReactionOldFormat(fname)) {
        cerr << "Cannot read old style reaction '" << fname << "'\n";
        return 0;
      }
    }
  }

  for (const IWReaction* r : _rxn) {
    if (r->number_sidechains() == 0) {
      continue;
    }

    for(const Sidechain_Reaction_Site* s : r->sidechains()) {
      if (s->number_reagents() > 1) {
        ++_perform_reactions_via_iterator;
      }
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    if (! ReadProtoConfig(fname)) {
      cerr << "Cannot read proto config file '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('a')) {
    _append_reaction_name_to_product = 1;
    if (_verbose) {
      cerr << "Will append the reaction name to the product\n";
    }
  }

  if (cl.option_present('V')) {
    _discard_bad_valence = 1;
    IWString fname = cl.string_value('V');
    if (fname == "drop") {
    } else {
      if (! fname.ends_with(".smi")) {
        fname << ".smi";
      }
      if (! _stream_for_bad_valence.open(fname)) {
        cerr << "Cannot open stream for bad valence '" << fname << "'\n";
        return 0;
      }
      if (_verbose) {
        cerr << "Will write bad valence molecules to '" << fname << "'\n";
      }
    }
  }

  if (cl.option_present('b')) {
    _re_search_each_scaffold_hit = 1;
    if (_verbose) {
      cerr << "Will sequentially perform substructure searches\n";
    }
  }

  if (cl.option_present('x')) {
    _apply_all_reactions = 1;
    if (_verbose) {
      cerr << "Will apply all reactions to the starting molecule, generating 1 product\n";
    }
  }

  if (cl.option_present('J')) {
    const_IWSubstring j;
    for (int i = 0; cl.value('J', j, i); ++i) {
      if (j == "noschmsg") {
        _smc.set_issue_sidechain_no_match_warnings(0);
        if (_verbose) {
          cerr << "Will NOT warn about no sidechain query matches\n";
        }
      } else if (j == "help") {
        DisplayDashJQualifiers(cerr);
      } else {
        cerr << "Unrecognised -J qualifier '" << j << "'\n";
        DisplayDashJQualifiers(cerr);
        return 1;
      }
    }
  }

  return 1;
}

int
MultipleReactions::ReadFileOfReactions(IWString& fname)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "MultipleReactions::ReadFileOfReactions:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFileOfReactions(fname, input);
}

int
MultipleReactions::ReadFileOfReactions(const IWString& fname, iwstring_data_source& input)
{
  IWString dirname(fname);
  dirname.truncate_at_last('/');
  dirname << '/';
  IWString buffer;
  while (input.next_record(buffer)) {
    IWString fname;
    fname << dirname << buffer;
    if (!ReadReaction(fname)) {
      cerr << "MultipleReactions::ReadFileOfReactions:Cannot read '" << buffer << "'\n";
      return 0;
    }
  }

  return _rxn.size();
}

int
MultipleReactions::ReadReaction(IWString& fname)
{
  std::optional<ReactionProto::Reaction> maybe_rxn =
      iwmisc::ReadTextProto<ReactionProto::Reaction>(fname);
  if (!maybe_rxn) {
    cerr << "MultipleReactions::ReadTextProto:cannot open '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<IWReaction> r = std::make_unique<IWReaction>();
  IWString dirname;  // not making allowance for embedding query file names here.
  if (!r->ConstructFromProto(*maybe_rxn, dirname)) {
    cerr << "MultipleReactions::ReadReaction:cannot parse '" << maybe_rxn->ShortDebugString()
         << '\n';
    return 0;
  }

  _rxn << r.release();

  return 1;
}

int
MultipleReactions::ReadReactionOldFormat(IWString& fname)
{
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "MultipleReactions::ReadReactionOldFormat:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReactionOldFormat(input, fname);
}

int
MultipleReactions::ReadReactionOldFormat(iwstring_data_source& input,
                        const IWString& fname) {

  // Currently dirname not being used...
  IWString dirname(fname);
  dirname.truncate_at_last('/');

  std::unique_ptr<IWReaction> r = std::make_unique<IWReaction>();
  if (r->do_read(input, _smc)) {
    _rxn << r.release();
    return 1;;
  }

  cerr << "MultipleReactions::ReadReactionOldFormat:cannot read reaction\n";
  return 0;
}

int
MultipleReactions::ReadTfDataReactions(IWString& fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (! input.good()) {
    cerr << "MultipleReactions::ReadTfDataReactions:cannot open '" << fname << "'\n";
    return 0;
  }

  fs::path dir = fs::path(fname.null_terminated_chars());
  const IWString dirname(dir.parent_path());

  return ReadTfDataReactions(input, dirname);
}

int
MultipleReactions::ReadTfDataReactions(iw_tf_data_record::TFDataReader& input,
                        const IWString& dirname) {
  while (true) {
    std::optional<const_IWSubstring> data = input.Next();
    if (! data) {
      return _rxn.size();
    }

    // Thought I would need this, but not so. Keep as an example.
    // using google::protobuf::io::ZeroCopyInputStream;
    // using google::protobuf::io::ArrayInputStream;
    // ArrayInputStream zero_copy_input(data->data(), data->length());

    ReactionProto::Reaction proto;
    if (! proto.ParseFromArray(data->data(), data->length())) {
      cerr << "MultipleReactions::ReadTfDataReactions:cannot decode\n";
      return 0;
    }

    std::unique_ptr<IWReaction> rxn = std::make_unique<IWReaction>();
    if (! rxn->ConstructFromProto(proto, dirname)) {
      cerr << "MultipleReactions::ReadTfDataReactions:invalid proto " << 
              proto.ShortDebugString();
      return 0;
    }

    _rxn << rxn.release();
  }

  return _rxn.size();
}

int
MultipleReactions::ReadProtoConfig(IWString& fname) {
  std::optional<multiple_reactions::Options> proto = 
     iwmisc::ReadTextProto<multiple_reactions::Options>(fname);
  if (! proto) {
    cerr << "MultipleReactions::ReadProtoConfig:cannot read '" << fname << "'\n";
    return 0;
  }

  for (const auto& rxn : proto->rxn()) {
    if (! ReadReaction(rxn)) {
      cerr << "MultipleReactions::ReadProtoConfig:cannot parse " << rxn.ShortDebugString() << '\n';
      return 0;
    }
  }

  _apply_all_reactions = proto->apply_all_reactions_to_reagent();

  return 1;
}

int
MultipleReactions::ReadReaction(const multiple_reactions::ReactionAndSidechains& proto) {
  if (proto.has_proto_reaction_file()) {
  } else if (proto.has_legacy_reaction_file()) {
  } else {
    cerr << "MultipleReactions::ReadReaction:no reaction specified\n";
    return 0;
  }

  std::unique_ptr<IWReaction> rxn = std::make_unique<IWReaction>();
  if (proto.has_proto_reaction_file()) {
    IWString fname(proto.proto_reaction_file());
    std::unique_ptr<ReactionProto::Reaction> proto = 
                iwmisc::ReadTextProtoPtr<ReactionProto::Reaction>(fname);
    if (! proto) {
      cerr << "MultipleReactions::ReadReaction:cannot read textproto '" << fname << "'\n";
      return 0;
    }
    if (! rxn->ConstructFromProto(*proto, fname)) {
      cerr << "MultipleReactions::ReadReaction:cannot build reaction from '" << fname << "'\n";
      return 0;
    }
  } else if (proto.has_legacy_reaction_file()) {
    IWString fname(proto.legacy_reaction_file());
    if (! rxn->do_read(fname, _smc)) {
      cerr << "MultipleReactions::ReadReaction:cannot read msi file '" << fname << "'\n";
      return 0;
    }
  }

  assert(rxn);

  if (proto.sidechain_file_size() == 0) {
    _rxn << rxn.release();
    return 1;
  }

  if (rxn->number_sidechains() == 0) {
    cerr << "MultipleReactions::ReadReaction:no sidechains in reaction, but sidechain file specified\n";
    return 0;
  }

  // In index into proto.sidechain_file()
  int ndx = 0;
  for (int i = 0; i < rxn->number_sidechains(); i++) {
    Sidechain_Reaction_Site* s = rxn->sidechain(i);

    if (s->single_reagent_only()) {
      continue;
    }

    IWString fname = proto.sidechain_file(ndx);
    ++ndx;
    if (! s->add_reagents(fname, FILE_TYPE_SMI, _smc)) {
      cerr << "MultipleReactions::ReadReaction:cannot read sidechain '" << fname << "'\n";
      return 0;
    }
  }

  if (ndx != proto.sidechain_file_size()) {
    cerr << "MultipleReactions::ReadReaction:size mismatch between sidechains " << ndx << '\n';
    return 0;
  }

  _rxn << rxn.release();

  return 1;
}

int
MultipleReactions::Preprocess(Molecule& m)
{
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
MultipleReactions::Process(Molecule& m, resizable_array_p<Molecule>& results)
{
  ++_molecules_read;

  _seen.insert(m.unique_smiles());

  if (_perform_reactions_via_iterator) {
    return ViaIterator(m, results);
  }

  if (_apply_all_reactions) {
    return ProcessAllReactionsApplied(m, results);
  }

  Molecule_to_Match target(&m);

  int matches = 0;
  const int nrxn = _rxn.number_elements();
  for (int i = 0; i < nrxn; ++i) {
    if (Process(m, target, *_rxn[i], results)) {
      ++_rxn_matched[i];
      ++matches;
    }
  }

  ++_reactions_matching[matches];

  // Initially I had the check for zero matches being a failure, but something
  // can fail if it only generates products that have already been seen.
  // We could detect that, but too messy.
  if (matches == 0) {
    if (_verbose > 1) {
      cerr << "MultipleReactions::Process:Duplicates or none of " << _rxn.size() << " reactions reacted with '"
           << m.smiles() << ' ' << m.name() << "'\n";
      cerr << "Ignored\n";
      ++_molecules_not_generating_products;
    }

    return 1;
  }

  return 1;
}

int
MultipleReactions::Process(Molecule& m, Molecule_to_Match& target, IWReaction& rxn,
                           resizable_array_p<Molecule>& result)
{
  if (_re_search_each_scaffold_hit) {
    return ProcessReSearch(m, rxn, result);
  }

  Substructure_Results sresults;
  const int nhits = rxn.substructure_search(target, sresults);
  if (_verbose > 1) {
    cerr << nhits << " hits to " << rxn.name() << " in " << m.name() << '\n';
  }

  if (nhits == 0) {
    return 0;
  }

  int rc = 0;
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    if (Process(m, *e, rxn, result)) {
      ++rc;
    }
  }

#ifdef DEBUG_PROCESS
  cerr << "From " << nhits << " rc " << rc << '\n';
#endif
  return rc;
}

// Apply all reactions to `m`, generating a single product.
int
MultipleReactions::ProcessAllReactionsApplied(Molecule& m,
                resizable_array_p<Molecule>& result) {
  int got_match = 0;

  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
  Substructure_Results sresults;

  for (int i = 0; i < _rxn.number_elements(); ++i) {
    const int nhits = _rxn[i]->substructure_search(*mcopy, sresults);
    if (nhits == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      _rxn[i]->in_place_transformation(*mcopy, e);
    }
    ++got_match;
    ++_rxn_matched[i];
  }

  if (! got_match) {
    cerr << "MultipleReactions::ProcessAllReactionsApplied:no matches " << m.name() << '\n';
    return 0;
  }

  result << mcopy.release();

  return got_match;
}

int
MultipleReactions::ProcessReSearch(Molecule& m,
                                   IWReaction& rxn,
                                   resizable_array_p<Molecule>& result) {
  // In order to avoid an infinite loop, only look this many times for a motif.
  constexpr int kMaxTries = 100;

  int rc = 0;

  Substructure_Results sresults;
  for (int i = 0; i < kMaxTries; ++i) {
    const int nhits = rxn.substructure_search(&m, sresults);
    if (nhits == 0) {
      return 0;
    }

    if (Process(m, *sresults.embeddings().first(), rxn, result)) {
      ++rc;
    }
  }

  return rc;
}

int
MultipleReactions::Process(Molecule& m, const Set_of_Atoms& embedding, IWReaction& rxn,
                           resizable_array_p<Molecule>& result)
{
  std::unique_ptr<Molecule> prod = std::make_unique<Molecule>(m);

  if (!rxn.in_place_transformation(*prod.get(), &embedding)) {
    if (_verbose > 2) {
      cerr << "in_place_transformation failed\n";
    }
    return 0;
  }

  if (IsDuplicate(*prod.get())) {
    if (_verbose > 2) {
      cerr << "Duplicate found\n";
    }
    return 0;
  }

  MaybeAppendRxnName(rxn.name(), *prod.get());

  if (! _discard_bad_valence) {
  } else if (prod->valence_ok()) {
  } else if (_stream_for_bad_valence.active()) {
    _stream_for_bad_valence << prod->smiles() << ' ' << prod->name() << '\n';
    _stream_for_bad_valence.write_if_buffer_holds_more_than(32768);
    ++_bad_valence_detected;
    return 1;
  } else {
    ++_bad_valence_detected;
    return 1;
  }

  result << prod.release();
  return 1;
}

void
MultipleReactions::MaybeAppendRxnName(const IWString& rxnname, Molecule& product)
{
  if (!_append_reaction_name_to_product) {
    return;
  }

  IWString new_name(product.name());
  new_name << " %% " << rxnname;
  product.set_name(new_name);
}

int
MultipleReactions::IsDuplicate(Molecule& m)
{
  const IWString& usmi = m.unique_smiles();
  const auto iter = _seen.find(usmi);
  //  if we have seen it before, return now.
  if (iter != _seen.end()) {
    ++_duplicates_discarded;
    return 1;
  }

  _seen.insert(usmi);
  //  Have not seen it before.
  return 0;
}

int
MultipleReactions::ViaIterator(Molecule& m, resizable_array_p<Molecule>& results) {
  const int nrxn = _rxn.number_elements();

  std::unique_ptr<Reaction_Iterator[]> iter = std::make_unique<Reaction_Iterator[]>(nrxn);

  for (int i = 0; i < nrxn; ++i) {
    iter[i].initialise(*_rxn[i]);
  }

  Molecule_to_Match target(&m);

  Substructure_Results sresults;

  for (int i = 0; i < nrxn; ++i) {
    if (_rxn[i]->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      for (iter[i].reset(); iter[i].active(); iter[i]++) {
        std::unique_ptr<Molecule> result = std::make_unique<Molecule>();
        if (! _rxn[i]->perform_reaction(&m, embedding, iter[i], *result)) {
          continue;
        }

        results << result.release();
      }
    }

    ViaIterator(m, i, iter.get(), results);
  }

  return 1;
}

int
MultipleReactions::ViaIterator(Molecule& m,
                        int ndx,
                        Reaction_Iterator* iter,
                        resizable_array_p<Molecule>& results) {
  const int nrxn = _rxn.number_elements();

  Molecule_to_Match target(&m);

  Substructure_Results sresults;

  for (; ndx < nrxn; ++ndx) {
    if (_rxn[ndx]->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      for (iter[ndx].reset(); iter[ndx].active(); iter[ndx]++) {
        std::unique_ptr<Molecule> result = std::make_unique<Molecule>();
        if (! _rxn[ndx]->perform_reaction(&m, embedding, iter[ndx], *result)) {
          continue;
        }

        if (ndx < (_rxn.number_elements() - 1)) {
          ViaIterator(*result, ndx + 1, iter, results);
        }
        results << result.release();
      }
    }
  }

  return 1;
}

int
MultipleReactions::Report(std::ostream& output) const
{
  output << "MultipleReactions:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }
  output << "Discarded " << _duplicates_discarded << " duplicates\n";
  output << "Discarded " << _bad_valence_detected << " products with bad valence\n";
  output << _molecules_not_generating_products << " molecules did not generate any products\n";

  for (int i = 0; i < _reactions_matching.number_elements(); ++i) {
    if (_reactions_matching[i]) {
      output << _reactions_matching[i] << " molecules matched " << i << " reactions\n";
    }
  }

  int istop = _rxn_matched.number_elements();
  if (istop > _rxn.number_elements()) {
    istop = _rxn.number_elements();
  }
  for (int i = 0; i < istop; ++i) {
    output << _rxn_matched[i] << " molecules matched " << _rxn[i]->name() << '\n';
  }

  return 1;
}

//  Options for the main program.
struct Options {
  int verbose;

  int write_starting_molecule;

  int ignore_molecules_not_matching_reactions;
  int write_molecules_not_reacting;

  int molecules_read;
  int molecules_not_reacting;

  int recursion_depth;

  extending_resizable_array<int> products_generated;

  FileType _input_type;

 public:

  Options();

  int Initialise(Command_Line& cl);

  int Report(std::ostream& output) const;

  FileType input_type() const {
    return _input_type;
  }
};

Options::Options()
{
  verbose = 0;
  write_starting_molecule = 0;
  ignore_molecules_not_matching_reactions = 0;
  write_molecules_not_reacting = 0;

  recursion_depth = 0;

  molecules_read = 0;
  molecules_not_reacting = 0;

  _input_type = FILE_TYPE_INVALID;
}

int
Options::Initialise(Command_Line& cl)
{
  verbose = cl.option_present('v');

  if (cl.option_present('p')) {
    write_starting_molecule = 1;
    if (verbose) {
      cerr << "Will write the starting molecule\n";
    }
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        ignore_molecules_not_matching_reactions = 1;
      }
      else if (z == 'w') {
        write_molecules_not_reacting = 1;
      }
      else {
        cerr << "Options::Initialise:unrecognized -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', recursion_depth) || recursion_depth < 0) {
      cerr << "option_present::Initialise:the recursion depth (-r) option must be a whole number\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will recurse " << recursion_depth << " times\n";
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  //  reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  }
  else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  }
  else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const
{
  output << "Options:read " << molecules_read << " molecules\n";
  output << molecules_not_reacting << " molecules did not produce results\n";
  for (int i = 0; i < products_generated.number_elements(); ++i) {
    if (products_generated[i]) {
      output << products_generated[i] << " molecules generated " << i << " products\n";
    }
  }

  return 1;
}

void
HandleMoleculesNotMatching(Options& options, Molecule& m,
                           const resizable_array_p<Molecule>& results)
{
  ++options.molecules_not_reacting;
  ++options.products_generated[results.size()];  //  should be zero...
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options, Molecule& m,
                       int recursion_depth, resizable_array_p<Molecule>& results)
{
  if (many_reactions.Process(m, results) == 0) {
    HandleMoleculesNotMatching(options, m, results);
    if (recursion_depth > 0) {
      return 1;
    }

    if (!options.ignore_molecules_not_matching_reactions) {
      cerr << "Error processing " << m.name() << '\n';
      return 0;
    }

    return 1;
  }

  ++options.products_generated[results.size()];

  if (recursion_depth < options.recursion_depth) {
    const int existing = results.number_elements();
    for (int i = 0; i < existing; ++i) {
      ApplyMultipleReactions(many_reactions, options, *results[i], recursion_depth + 1, results);
    }
  }

  return 1;
}

//  Write the results of transformations.
int
DoOutput(const Options& options, Molecule& m, resizable_array_p<Molecule>& results,
         Molecule_Output_Object& output)
{
  if (options.write_starting_molecule) {
    output.write(m);
  }

  for (Molecule* r : results) {
    output.write(*r);
  }

  return 1;
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options,
                       data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!many_reactions.Preprocess(*m)) {
      return 0;
    }

    ++options.molecules_read;

    resizable_array_p<Molecule> results;
    if (!ApplyMultipleReactions(many_reactions, options, *m, 0, results)) {
      return 0;
    }

    DoOutput(options, *m, results, output);
  }

  return 1;
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options, const char* fname,
                       Molecule_Output_Object& output)
{
  FileType input_type = options.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ApplyMultipleReactions(many_reactions, options, input, output);
}

int
Main(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:A:i:g:lcR:S:az:pr:V:bG:xJ:C:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  MultipleReactions many_reactions;
  if (!many_reactions.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    Usage(1);
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determing output type(s) (-o)\n";
      return 1;
    }
  }
  else {
    output.add_output_type(FILE_TYPE_SMI);
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    return 1;
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 1;
    }

    if (!output.new_stem(fname)) {
      cerr << "Cannot open output file (-S) '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written to '" << fname << "'\n";
    }
  }

  for (const char* fname : cl) {
    if (!ApplyMultipleReactions(many_reactions, options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
    many_reactions.Report(cerr);
    cerr << "Write " << output.molecules_written() << " molecules\n";
  }

  return 0;
}

}  //  namespace apply_multiple_reactions

int
main(int argc, char** argv)
{
  int rc = apply_multiple_reactions::Main(argc, argv);

  return rc;
}
