#include "Molecule_Tools/medchemwizard_lib.h"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/target.h"

namespace medchemwizard {

using std::cerr;

namespace {

int
ReadTextProtoReaction(const IWString& fname,
                      resizable_array_p<IWReaction>& destination) {
  IWString proto_fname(fname);
  std::optional<ReactionProto::Reaction> maybe_proto =
      iwmisc::ReadTextProto<ReactionProto::Reaction>(proto_fname);
  if (! maybe_proto) {
    cerr << "ReadTextProtoReaction:cannot read '" << fname << "'\n";
    return 0;
  }

  Sidechain_Match_Conditions smc;

  std::unique_ptr<IWReaction> rxn = std::make_unique<IWReaction>();
  if (!rxn->ConstructFromProto(*maybe_proto, fname, smc)) {
    cerr << "ReadTextProtoReaction:cannot build reaction from proto\n";
    cerr << maybe_proto->ShortDebugString() << '\n';
    return 0;
  }

  destination << rxn.release();

  return 1;
}

int
ReadReaction(const_IWSubstring& buffer, const IWString& dir,
             const Sidechain_Match_Conditions& smc, resizable_array_p<IWReaction>& rxn,
             int verbose) {
  bool reaction_is_textproto;
  if (buffer.starts_with("PROTO:")) {
    buffer.remove_leading_chars(6);
    reaction_is_textproto = true;
  } else {
    reaction_is_textproto = false;
  }

  IWString fname;
  fname << dir << '/' << buffer;

  if (verbose > 2) {
    cerr << "Reading '" << fname << "'\n";
  }

  if (reaction_is_textproto) {
    return ReadTextProtoReaction(fname, rxn);
  }

  iwstring_data_source input(fname.null_terminated_chars());

  if (!input.good()) {
    cerr << "Cannot open reaction '" << fname << "'\n";
    return 0;
  }

  msi_object msi;
  msi.set_display_no_data_error_message(0);

  while (msi.read(input)) {
    IWReaction* r = new IWReaction();
    if (!r->construct_from_msi_object(msi, smc)) {
      cerr << "Cannot build reaction\n";
      cerr << msi;
      delete r;

      return 0;
    }

    rxn.add(r);
  }

  return 1;
}

int
ReadReactionsFromList(iwstring_data_source& input, const IWString& dir,
                      const Sidechain_Match_Conditions& smc,
                      resizable_array_p<IWReaction>& rxn, int verbose) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (0 == buffer.length() || buffer.starts_with('#')) {
      continue;
    }

    if (!ReadReaction(buffer, dir, smc, rxn, verbose)) {
      cerr << "Cannot read reaction '" << buffer << "'\n";
      return 0;
    }
  }

  return rxn.number_elements();
}

}  // namespace

MedchemWizard::MedchemWizard() {
}

int
MedchemWizard::ReadReactions(const char* fname) {
  set_warn_no_query_in_reaction_object(0);

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "MedchemWizard::ReadReactions:cannot open '" << fname << "'\n";
    return 0;
  }

  IWString dir;

  const_IWSubstring tmp(fname);

  const int slash = tmp.rindex('/');

  if (slash >= 0) {
    dir = fname;
    dir.iwtruncate(slash);
  } else {
    dir = "./";
  }

  Sidechain_Match_Conditions smc;
  if (! ReadReactionsFromList(input, dir, smc, _rxn, _options.verbose)) {
    return 0;
  }

  _molecules_hitting_reaction = std::make_unique<int[]>(_rxn.number_elements());

  return _rxn.number_elements();
}

int
MedchemWizard::ReadReactionFile(const char* fname, bool textproto) {
  set_warn_no_query_in_reaction_object(0);
  const_IWSubstring path(fname);
  const int slash = path.rindex('/');
  IWString dir(".");
  IWString record;
  if (textproto) {
    record << "PROTO:";
  }
  if (slash >= 0) {
    dir = fname;
    dir.iwtruncate(slash == 0 ? 1 : slash);
    record << (fname + slash + 1);
  } else {
    record << fname;
  }

  // Legacy files may contain multiple reactions. Do not keep a partially read file.
  resizable_array_p<IWReaction> pending;
  Sidechain_Match_Conditions smc;
  const_IWSubstring buffer(record);
  if (!ReadReaction(buffer, dir, smc, pending, _options.verbose) || pending.empty()) {
    return 0;
  }
  const int added = pending.number_elements();
  const int previous = _rxn.number_elements();
  _rxn.transfer_in(pending);
  auto hits = std::make_unique<int[]>(_rxn.number_elements());
  for (int i = 0; i < previous; ++i) {
    hits[i] = _molecules_hitting_reaction[i];
  }
  _molecules_hitting_reaction = std::move(hits);
  return added;
}

int
MedchemWizard::InitialiseFromEnvironment() {
  const char* lillymol_home = getenv("LILLYMOL_HOME");
  if (lillymol_home == nullptr) {
    cerr << "MedchemWizard::InitialiseFromEnvironment:LILLYMOL_HOME not defined\n";
    return 0;
  }

  IWString fname(lillymol_home);
  fname << "/data/MedchemWizard/REACTIONS";

  return ReadReactions(fname.null_terminated_chars());
}

void
MedchemWizard::Preprocess(Molecule& m) {
  if (_options.reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_options.remove_all_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  if (_options.chemical_standardisation.active()) {
    _options.chemical_standardisation.process(m);
  }
}

int
MedchemWizard::IdentifyAtomsThatMustNotChange(Molecule& m,
                                              std::unique_ptr<int[]>& do_not_change) {
  if (_do_not_change_query.empty()) {
    return 1;
  }

  Molecule_to_Match target(&m);
  int matched = 0;

  for (Substructure_Query* query : _do_not_change_query) {
    Substructure_Results sresults;
    const int nhits = query->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    if (! do_not_change) {
      do_not_change = std::make_unique<int[]>(m.natoms());
    }
    sresults.each_embedding_set_vector(do_not_change.get(), 1);
    matched += nhits;
  }

  if (matched > 0) {
    return 1;
  }

  ++_stats.molecules_not_matching_do_not_change_queries;
  if (_options.ignore_do_not_change_queries_not_matching) {
    return 1;
  }

  cerr << "MedchemWizard::IdentifyAtomsThatMustNotChange:no do-not-change query matched '"
       << m.name() << "'\n";
  return 0;
}

int
MedchemWizard::CheckUniqueness(Molecule& m, int query_just_run,
                               IW_STL_Hash_Set& smiles_hash) {
  if (_options.postprocess_molecules_produced) {
    Preprocess(m);
  }

  if (_options.discard_bad_valences && !m.valence_ok()) {
    if (_options.verbose > 2) {
      cerr << "bad valence " << m.smiles() << ' ' << m.name() << ", last rxn "
           << query_just_run << ' ' << _rxn[query_just_run]->comment() << '\n';
    }
    _stats.bad_valences_discarded++;

    if (_stream_for_bad_valence != nullptr && _stream_for_bad_valence->is_open()) {
      (*_stream_for_bad_valence) << m.smiles() << ' ' << m.name() << '\n';
      _stream_for_bad_valence->write_if_buffer_holds_more_than(8192);
    }

    return 0;
  }

  const int matoms = m.natoms();

  if (matoms > _options.max_atoms) {
    _stats.discarded_for_too_many_atoms++;
    return 0;
  }

  if (matoms < _options.min_atoms) {
    _stats.discarded_for_too_few_atoms++;
    return 0;
  }

  if (0 == _options.unique_across_all_molecules && 0 == _options.unique_within_molecule) {
    return 1;
  }

  const auto& usmi = m.unique_smiles();

  if (smiles_hash.contains(usmi)) {
    _stats.duplicate_molecules_suppressed++;
    return 0;
  }

  smiles_hash.insert(usmi);

  return 1;
}

int
MedchemWizard::MaybeTruncateNhits(int nhits, IWReaction& reaction, Molecule& m) {
  if (nhits <= _options.max_hits) {
    return nhits;
  }

  if (_options.report_multiple_hits_truncated) {
    cerr << nhits << " matches to " << reaction.comment() << " in " << m.name()
         << " in " << m.name() << " truncated\n";
  }
  ++_stats.truncated_to_max_hits;

  return _options.max_hits;
}

int
MedchemWizard::GenerateProducts(Molecule& m, const IWString& mname, int depth,
                       const int* do_not_change, IW_STL_Hash_Set& smiles_hash,
                       ProductCallback& callback) {
  Molecule_to_Match target(&m);

  const int n = _rxn.number_elements();

  IWString product_molecule_name(mname);
  const int initial_name_length = mname.length();

  for (int i = 0; i < n; ++i) {
    const auto r = _rxn[i];

    Substructure_Results sresults;

    int nhits = r->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    _molecules_hitting_reaction[i]++;

    nhits = MaybeTruncateNhits(nhits, *r, m);

    for (int j = 0; j < nhits; ++j) {
      const Set_of_Atoms* embedding = sresults.embedding(j);
      if (do_not_change != nullptr && embedding->any_members_set_in_array(do_not_change)) {
        ++_stats.embeddings_rejected_for_changing_protected_atoms;
        continue;
      }

      std::unique_ptr<Molecule> product = std::make_unique<Molecule>();
      if (!r->perform_reaction(&m, embedding, *product)) {
        continue;
      }

      _stats.molecules_produced++;

      if (!CheckUniqueness(*product, i, smiles_hash)) {
        continue;
      }

      if (_options.append_names) {
        product_molecule_name << _options.sep << r->comment();
      }

      product->set_name(product_molecule_name);
      Molecule* product_ptr = product.get();

      if (! callback(product)) {
        return 0;
      }

      if (depth < _options.max_depth) {
        if (! GenerateProducts(*product_ptr, product_molecule_name, depth + 1, do_not_change,
                              smiles_hash, callback)) {
          return 0;
        }
      }

      product_molecule_name.resize_keep_storage(initial_name_length);
    }
  }

  return 1;
}

int
MedchemWizard::Process(Molecule& m, ProductCallback callback) {
  if (_molecules_hitting_reaction == nullptr) {
    cerr << "MedchemWizard::Process:no reactions loaded\n";
    return 0;
  }
  if (! callback) {
    cerr << "MedchemWizard::Process:no product callback specified\n";
    return 0;
  }

  _stats.molecules_read++;
  Preprocess(m);

  std::unique_ptr<int[]> do_not_change;
  if (! IdentifyAtomsThatMustNotChange(m, do_not_change)) {
    return 0;
  }

  IW_STL_Hash_Set within_molecule_hash;
  IW_STL_Hash_Set& smiles_hash =
      _options.unique_within_molecule ? within_molecule_hash : _smiles_hash;

  IWString mname(m.name());
  return GenerateProducts(m, mname, 0, do_not_change.get(), smiles_hash, callback);
}

int
MedchemWizard::ProcessToArray(Molecule& m, resizable_array_p<Molecule>& output) {
  return Process(m, [&output](std::unique_ptr<Molecule>& product) {
    output << product.release();
    return 1;
  });
}

int
MedchemWizard::ProcessToStream(Molecule& m, IWString_and_File_Descriptor& output) {
  return Process(m, [&output](std::unique_ptr<Molecule>& product) {
    output << product->smiles() << ' ' << product->name() << '\n';
    output.write_if_buffer_holds_more_than(8192);
    return output.good();
  });
}

}  // namespace medchemwizard
