#include <iostream>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "spread_weights.h"

namespace spread_weights {

using std::cerr;

ScalingFactor::ScalingFactor() {
  _scaling_factor_column = -1;
  _column_in_weight_file = 1;
  _active = false;
}

void
DisplayHelpMessage(char flag, std::ostream& output) {
  output << "Specifies per item weighting values for spread\n";
  output << " -" << flag << " COL=<col>         weight values are a column in the name\n";
  output << " -" << flag << " NOHDR             the weight file does NOT have a header record\n";
  output << " -" << flag << " FNAME=<fname>     weight values are in an external file\n";
  output << "                   assumed to be space separated with weights in column 2\n";
  output << "                   Id Weight    (header)\n";
  output << "                   id1 w1\n";
  output << "                   id2 w2\n";

}

int
ScalingFactor::Initialise(Command_Line& cl, char flag, bool verbose) {
  if (! cl.option_present(flag)) {
    _active = false;
    return 1;
  }

  int weight_file_has_header = 1;
  IWString fname;
  const_IWSubstring s;
  for (int i = 0; cl.value(flag, s, i); ++i) {
    if (s.starts_with("FILE=")) {
      s.remove_leading_chars(5);
      fname = s;
    } else if (s == "NOHDR") {
      weight_file_has_header = 0;
      if (verbose) {
        cerr << "Weight file does NOT have a header record\n";
      }
    } else if (s.starts_with("COL=")) {
      s.remove_leading_chars(4);
      if (! s.numeric_value(_scaling_factor_column) || _scaling_factor_column < 1) {
        cerr << "ScalingFactor::Initialise::invalid scalling factor column '" << s << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Weights are column " << _scaling_factor_column << " in the name\n";
      }
    } else if (s == "help") {
      DisplayHelpMessage(flag, cerr);
    } else {
      cerr << "Unrecognised -" << flag << " qualifier - scaling factor\n";
      DisplayHelpMessage(flag, cerr);
    }
  }

  if (! fname.empty()) {
    if (! ReadWeights(fname, weight_file_has_header)) {
      cerr << "ScalingFactor::Initialise:cannot read weight data from '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Read " << _id_to_weight.size() << " id to weight values from '" << fname << "'\n";
    }
  }

  return 1;
}

int
ScalingFactor::ReadWeights(const IWString& fname, int weight_file_has_header) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ScalingFactor::ReadWeights:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadWeights(input, weight_file_has_header);
}

int
ScalingFactor::ReadWeights(iwstring_data_source& input, int weight_file_has_header) {
  const_IWSubstring buffer;

  // No checking here.
  if (weight_file_has_header) {
    input.next_record(buffer);
  }

  while (input.next_record(buffer)) {
    if (! ReadWeightRecord(buffer)) {
      cerr << "ScalingFactor::ReadWeightRecord:Invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  if (_id_to_weight.empty()) {
    cerr << "ScalingFactor::ReadWeights:no data\n";
    return 0;
  }

  return _id_to_weight.size();
}

int
ScalingFactor::ReadWeightRecord(const const_IWSubstring& buffer) {
  IWString id, weight_as_string;

  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i, ' '); ++col) {
    if (col == 0) {
      id = token;
    } else if (col == _column_in_weight_file) {
      weight_as_string = token;
    }
  }

  if (id.empty() || weight_as_string.empty()) {
    cerr << "ScalingFactor::ReadWeightRecord:empty id '" << id << " or weight '" <<
             weight_as_string << "'\n";
    return 0;
  }

  float w;
  if (! weight_as_string.numeric_value(w)) {
    cerr << "Invalid weight '" << buffer << "'\n";
    return 0;
  }

  _id_to_weight[id] = w;

  return 1;
}

std::optional<float>
ScalingFactor::Weight(const IWString& name) const {
  if (! _active) {
    return std::nullopt;
  }

  if (_scaling_factor_column >= 0) {
    return WeightFromTokenInName(name);
  }

  auto iter = _id_to_weight.find(name);
  if (iter != _id_to_weight.end()) {
    return iter->second;
  }

  int ndx = name.index(' ');
  if (ndx < 0) {
    return std::nullopt;
  }

  IWString tmp(name.data(), ndx - 1);
  iter = _id_to_weight.find(tmp);
  if (iter == _id_to_weight.end()) {
    return std::nullopt;
  }

  return iter->second;
}

std::optional<float>
ScalingFactor::WeightFromTokenInName(const IWString& name) const {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; name.nextword(token, i, ' '); ++col) {
    if (col != _scaling_factor_column) {
      continue;
    }

    float w;
    if (! token.numeric_value(w) || w <= 0.0f) {
      cerr << "ScalingFactor::WeightFromTokenInName:invalid weight in name '" << name << "'\n";
      return std::nullopt;
    }

    return w;
  }

  return std::nullopt;
}

}  // namespace spread_weights
