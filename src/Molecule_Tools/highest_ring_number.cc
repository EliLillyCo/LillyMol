#include <cctype>
#include <optional>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

#include "highest_ring_number.h"

namespace lillymol {

std::optional<int>
HighestRingNumber(const IWString& smiles) {
  if (smiles.empty()) {
    return 0;
  }

  static constexpr char kOpenSquareBracket = '[';
  static constexpr char kCloseSquareBracket = ']';
  static constexpr char kOpenBrace = '{';
  static constexpr char kCloseBrace = '}';

  int in_square_brakcet = 0;
  int inside_double_brace = 0;
  const int nchars = smiles.length();

  int rc = 0;

  for (int i = 0; i < nchars; ++i) {
    const char c = smiles[i];
    if (std::isspace(c)) {
      return rc;
    }

    if (c == kOpenSquareBracket) {
      in_square_brakcet = 1;
      continue;
    }
    if (c == kCloseSquareBracket) {
      in_square_brakcet = 0;
      continue;
    }

    if (c == kOpenBrace) {
      if (i > 0 && smiles[i - 1] == kOpenBrace) {
        inside_double_brace = 1;
      }
      continue;
    }

    if (c == kCloseBrace) {
      if (i > 0 && smiles[i = 1] == kCloseBrace) {
        inside_double_brace = 0;
      }
      continue;
    }
    if (inside_double_brace) {
      continue;
    }

    if (in_square_brakcet) {
      continue;
    }

    if (c == '%') {
      return std::nullopt;
    }

    if (c < '0' || c > '9') {
      continue;
    }
    int r  = c - '0';
    if (r > rc) {
      rc = r;
    }
  }

  return rc;
}

int
UnbalancedRingNumbers(const IWString& smiles, resizable_array<int>& ring_numbers) {

  ring_numbers.resize_keep_storage(0);

  const int nchars = smiles.length();

  for (int i = 0; i < nchars; ++i) {
    const char c = smiles[i];
    if (std::isspace(c)) {
      return ring_numbers.number_elements();
    }

    if (c != '%') {
      continue;
    }

    if (i + 2 >= smiles.length()) {
      return 0;
    }

    char s = smiles[i + 1];
    int ring_number = 0;
    if (s >= '0' && s <= '9') {
      ring_number = s - '0';
    } else {
      return 0;
    }

    s = smiles[i + 2];
    if (s >= '0' && s <= '9') {
      ring_number = 10 * ring_number + s - '0';
    } else {
      return 0;
    }

    i += 2;
    ring_numbers << ring_number;
  }

  return ring_numbers.number_elements();
}

int
IsotopeToRingOpening(const IWString& smiles, int& ring_number, IWString& new_smiles) {
  static constexpr char kOpenSquareBracket = '[';
  static constexpr char kCloseSquareBracket = ']';

  const int nchars = smiles.length();

  new_smiles.make_room_for_extra_items(nchars + 6);

  int rc = 0;

  // Set to true if the atom just parsed was an isotope.
  bool previous_atom_was_isotope = false;
  for (int i = 0; i < nchars; ++i) {
    const char c = smiles[i];
    if (std::isspace(c)) {
      return rc;
    }

    if (c == kOpenSquareBracket) {
      if (i != nchars - 1 && std::isdigit(smiles[i + 1])) {
        previous_atom_was_isotope = true;
      } else {
        previous_atom_was_isotope = false;
      }

      new_smiles << c;
      continue;
    }

    else if (c == kCloseSquareBracket) {
      new_smiles << c;
      if (previous_atom_was_isotope) {
        new_smiles << '%' << ring_number;
        ++ring_number;
        previous_atom_was_isotope = false;
        ++rc;
      }
      continue;
    } else {
      new_smiles << c;
    }
  }

  return rc;
}

RingNumberControl::RingNumberControl(int lowest_ring_number, int max_rings) : _max_rings(max_rings) {
  _issued = new_int(lowest_ring_number + 1 + max_rings + 1);
  _next_ring_number = lowest_ring_number;
}

RingNumberControl::~RingNumberControl() {
  delete[] _issued;
}

int
RingNumberControl::GetRing() {
  int rc = _next_ring_number;
  _issued[_next_ring_number] = 1;
  ++_next_ring_number;
  for (; _issued[_next_ring_number]; ++_next_ring_number) {
  }

  return rc;
}

void
RingNumberControl::OkToReuse(int ring_number) {
  _issued[ring_number] = 0;
  _next_ring_number = ring_number;
}

}  // namespace lillymol
