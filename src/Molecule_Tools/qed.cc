#include <array>
#include <iostream>
#include <optional>

#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/qed.h"

namespace qed {

using std::cerr;

void
QEDProperties::Reset() {
  amw = 0.0f;
  alogp = 0.0f;
  hba = 0;
  hbd = 0;
  psa = 0.0f;
  rotb = 0;
  arom = 0;
  alerts = 0;
}

std::array<float, 8> qed_max {
 0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.0
};

std::array<float, 8> qed_mean = {
  0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95
};

std::array<float, 8> qed_none = {
  1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
};


ADSparameter amw_ads = {2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677,
               65.37051707, 104.9805561};
ADSparameter alogp_ads {3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154,
               0.576295591, 131.3186604};
ADSparameter hba_ads = {2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953,
               1.300669958, 148.7763046};
ADSparameter hbd_ads = {1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843,
               0.920922555, 258.1632616};
ADSparameter psa_ads = {1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824,
               28.51324732, 104.5686167};
ADSparameter rotb_ads = {0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166,
               2.758063707, 105.4420403};
ADSparameter arom_ads = {3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384,
               0.375760881, 312.3372610};
ADSparameter alerts_ads = {0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477,
               0.875193782, 417.7253140};

// ADS function
float
ads(float x, const ADSparameter& p) {
  double exp1 = 1.0 + exp(-1.0 * (x - p.c + p.d / 2.0) / p.e);
  double exp2 = 1.0 + exp(-1.0 * (x - p.c - p.d / 2.0) / p.f);
  double dx = p.a + p.b / exp1 * (1 - 1 / exp2);
  return dx / p.dmax;
}

void
Usage(std::ostream& output) {
}

Qed::Qed() {
  _alogp.set_rdkit_charged_nitrogen(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_rdkit_phoshoric_acid_hydrogen(1);

  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);

  // For compatibility with RDKit.
  nvrtspsa::set_zero_for_all_sulphur_atoms(1);
  nvrtspsa::set_display_psa_unclassified_atom_mesages(0);
}

int
Qed::Initialise(Command_Line& cl, char flag) {
  IWString s;
  for (int i = 0; cl.value(flag, s, i); ++i) {
    if (s.starts_with("QUERY=")) {
      s.remove_leading_chars(6);
      static constexpr int kVerbose = 0;
      if (! process_cmdline_token(flag, s, _queries, kVerbose)) {
        cerr << "Qed::Initialise:cannot read queries '" << s << "'\n";
        return 0;
      }
    } else if (s.starts_with("ACC=")) {
      s.remove_leading_chars(4);
      static constexpr int kVerbose = 0;
      if (! process_cmdline_token(flag, s, _acceptor_queries, kVerbose)) {
        cerr << "Qed::Initialise:cannot read acceptor queries '" << s << "'\n";
        return 0;
      }
    } else if (s.starts_with("DON=")) {
      s.remove_leading_chars(4);
      static constexpr int kVerbose = 0;
      if (! process_cmdline_token(flag, s, _donor_queries, kVerbose)) {
        cerr << "Qed::Initialise:cannot read donor queries '" << s << "'\n";
        return 0;
      }
    } else if (s.starts_with("ROTB=")) {
      s.remove_leading_chars(5);
      static constexpr int kVerbose = 0;
      if (! process_cmdline_token(flag, s, _rdkit_rotb_queries, kVerbose)) {
        cerr << "Qed::Initialise:cannot read rdkit rotb queries '" << s << "'\n";
        return 0;
      }
    } else if (s == "help") {
      Usage(cerr);
      return 0;
    } else {
      cerr << "Qed::Initialise:unrecogised directive '" << s << "'\n";
      Usage(cerr);
      return 0;
    }
  }

  if (_queries.empty()) {
    cerr << "Qed::Initialise:no queries specified\n";
    return 0;
  }

  for (Substructure_Query* q : _queries) {
    q->set_max_matches_to_find(1);
    q->set_save_matched_atoms(0);
  }

  for (Substructure_Query* q : _rdkit_rotb_queries) {
    q->set_find_unique_embeddings_only(1);
  }

  return 1;
}

int
Qed::Alerts(Molecule& m) {
  int rc = 0;

  Molecule_to_Match target(&m);

  Substructure_Results sresults;

  for (Substructure_Query* q : _queries) {
    if (q->substructure_search(target, sresults)) {
      ++rc;
      // rc += sresults.number_embeddings();
    }
  }

  return rc;
}

int
Qed::ExternalQueryCount(Molecule& m,
                resizable_array_p<Substructure_Query>& queries) {
  Molecule_to_Match target(&m);

  Substructure_Results sresults;

  int rc = 0;
  for (Substructure_Query* q : queries) {
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    rc += sresults.number_embeddings();
  }

  return rc;
}

void
Lipinski(Molecule& m, uint32_t& hba, uint32_t& hbd) {
  const int matoms = m.natoms();

  hba = 0;
  hbd = 0;

  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      continue;
    }

    int hcount = m.hcount(i);
    if (hcount == 0) {
      ++hba;
    } else {
      hbd += hcount;
    }
  }

  return;
}

bool
Qed::CalculateProperties(Molecule& m, QEDProperties& result) {

  result.Reset();

  m.reduce_to_largest_organic_fragment();

  if (! m.organic_only()) {
    result.alerts += 1;
  }
  if (m.ContainsIsotopicAtoms()) {
    result.alerts += 1;
  }

  result.amw = m.molecular_weight_ignore_isotopes();

  bool rc = true;

  std::optional<float> alogp = _alogp.LogP(m);
  if (alogp) {
    result.alogp = *alogp;
  } else {
    rc = false;
  }

  std::optional<double> psa = _nvrtspsa.PolarSurfaceArea(m);
  if (psa) {
    result.psa = *psa;
  } else {
    rc = false;
  }

  if (_acceptor_queries.empty() && _donor_queries.empty()) {
    Lipinski(m, result.hba, result.hbd);
  }

  if (! _acceptor_queries.empty()) {
    result.hba = ExternalQueryCount(m, _acceptor_queries);
  }

  if (! _donor_queries.empty()) {
    result.hbd = ExternalQueryCount(m, _donor_queries);
  }

  if (_rdkit_rotb_queries.empty()) {
    result.rotb = _rotbond.Process(m);
  } else {
    result.rotb = ExternalQueryCount(m, _rdkit_rotb_queries);
  }

  m.compute_aromaticity_if_needed();

  result.arom = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      ++result.arom;
    }
  }

  result.alerts += Alerts(m);

  return rc;
}

// Does not really need to be a member function.
float
Qed::ComputeQed(const QEDProperties& properties) const {

  std::array<float, 8> individual;
  individual[0] = ads(properties.amw, amw_ads);
  individual[1] = ads(properties.alogp, alogp_ads);
  individual[2] = ads(properties.hba, hba_ads);
  individual[3] = ads(properties.hbd, hbd_ads);
  individual[4] = ads(properties.psa, psa_ads);
  individual[5] = ads(properties.rotb, rotb_ads);
  individual[6] = ads(properties.arom, arom_ads);
  individual[7] = ads(properties.alerts, alerts_ads);

  // same variable name as in the python
  float t = 0.0f;
  float sumw = 0.0f;
  for (int i = 0; i < 8; ++i) {
//  cerr << "Individual " << i << " " << individual[i] << '\n';
    t += qed_mean[i] * log(individual[i]);
    sumw += qed_mean[i];
  }

  return exp(t / sumw);
}

std::optional<float>
Qed::qed(Molecule& m) {

  m.remove_all(1);

  QEDProperties properties;

  bool rc = CalculateProperties(m, properties);

  if (! rc) {
    return std::nullopt;
  }

  return ComputeQed(properties);
}

}  // namespace qed
