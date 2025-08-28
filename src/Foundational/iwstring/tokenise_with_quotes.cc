#include <iostream>

#include "Foundational/iwstring/iwstring.h"

namespace iwstring {

using std::cerr;

constexpr int kInvalid = -1;

int
TokeniseWithQuotes(const const_IWSubstring& buffer,
                   char sep,
                   resizable_array<int>& tstart,
                   resizable_array<int>& tstop) {
  tstart.resize_keep_storage(0);
  tstop.resize_keep_storage(0);

  static constexpr char kDQuote = '"';

  const int nchars = buffer.length();
  if (nchars == 0) {
    return 0;
  }

  // Maybe allow for empty token at start?
  if (buffer[0] == sep) {
    return kInvalid;
  }

  bool inside_quoted_string = false;
  int ntokens = 1;

  if (buffer[0] == kDQuote) {
    tstart << 1;
    inside_quoted_string = true;
  } else {
    tstart << 0;
  }

  for (int i = 1; i < nchars; ++i) {
    const char c = buffer[i];
    char next_char;
    if (i == nchars - 1) {
      next_char = '\0';
    } else {
      next_char = buffer[i + 1];
    }

    if (inside_quoted_string) {
      if (c == kDQuote && (next_char == sep || next_char == '\0')) {
        inside_quoted_string = false;
      }
    } else if (c == sep) {
      if (buffer[i-1] == kDQuote) {
        tstop << (i - 1);
      } else {
        tstop << (i - 0);
      }
      if (next_char == kDQuote) {
        tstart << (i + 2);
      } else {
        tstart << (i + 1);
      }
      ++ntokens;
    } else if (c == kDQuote && buffer[i-1] == sep) {
      inside_quoted_string = true;
    }
  }

  if (inside_quoted_string) {
    cerr << "TokeniseWithQuotes:unclosed quote '" << buffer << "'\n";
    return kInvalid;
  }

  if (buffer.ends_with(kDQuote)) {
    tstop << (buffer.length() - 1);
  } else {
    tstop << (buffer.length() - 0);
  }

  if (tstart.size() != tstop.size()) {
    cerr << "TokeniseWithQuotes::Mismatch between opening and closing tokens\n";
    cerr << tstart.size() << " vs " << tstop.size() << '\n';
    return kInvalid;
  }

  if (tstart.number_elements() != ntokens) {
    cerr << "TokeniseWithQuotes:Mismatch btw tokens " << ntokens <<
            " and array size " << tstart.size() << '\n';
    return kInvalid;
  }

// #define DEBUG_TOKENISE_WITH_QUOTES
#ifdef DEBUG_TOKENISE_WITH_QUOTES
  cerr << "Found " << ntokens << " tokens\n";
  for (int i = 0; i < tstart.number_elements(); ++i) {
    cerr << ' ' << i << " start " << tstart[i] << ' ' << buffer[tstart[i]] << 
            " stop " << tstop[i] << ' ' << buffer[tstop[i]] << '\n';
  }
#endif

  return ntokens;
}
}  // namespace iwstring
