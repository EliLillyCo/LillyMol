#ifndef MOLECULE_LIB_MISC2_H_
#define MOLECULE_LIB_MISC2_H_

#include <cstdint>
#include <memory>
#include <optional>

#include "google/protobuf/text_format.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"

#include "iwmtypes.h"

extern void iwabort();

extern int iw_rename(const char*, const char*);

extern int iw_getpid();

extern int int_comparitor_larger(const int*, const int*);

extern int uint64_comparitor_smaller(const uint64_t*, const uint64_t*);

extern void iwxor(const int*, int*, int);

class const_IWSubstring;

extern int fetch_numeric(const const_IWSubstring& string, int& value, int max_chars = 0);
extern int fetch_numeric_char(const char* string, int& value, int max_chars);

// Starting at s[ndx], see if there is a number.
// If a number is detected, it will be returned.
// 'ndx' will be incremented, and on success will be one past the last digit.
// If there is no digit at s[ndx], ndx will be unchanged.
extern std::optional<int> FetchNumeric(const const_IWSubstring& s, int& ndx);

/*
  Sometimes we need to compute combinatorial permutations and we
  may be dealing with numbers larger than can be held in an int
*/

extern uint64_t iw_combinatorial_combinations(int n, int k);

template <typename T>
int skip_to_string(T& input, const char* target, int report_discard);

// identify the 'split' characters in a reaction smiles.
// Complicated by the presence of + signs inside square brackets
extern int identify_plus_positions(const const_IWSubstring& buffer,
                                   resizable_array<int>& pos, char split);

// Return a vector of substrings delimited by 'split'.
// Does special handling of '+' in smiles, respecting square brackets.
// Beware scoping, const_IWSubstring are returned, pointint at 'buffer'.
extern int SplitOnPlusses(const const_IWSubstring& buffer,
                          resizable_array_p<const_IWSubstring>& parts, const char split);

namespace misc2 {

template <typename PROTO, typename T>
int
ReadProtoFile(const const_IWSubstring fname, T& destination) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    std::cerr << "ReadProtoFile:cannot open '" << fname;
    return 0;
  }

  IWString file_contents;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    file_contents << buffer;
  }

  const std::string string_proto(file_contents.rawdata(),
                                 file_contents.number_elements());

  PROTO proto;

  if (!google::protobuf::TextFormat::ParseFromString(string_proto, &proto)) {
    std::cerr << "ReadProtoFile:cannot parse proto\n";
    std::cerr << string_proto << '\n';
    return 0;
  }

  if (!destination.ConstructFromProto(proto)) {
    std::cerr << "ReadProtoFile:cannot build object from proto\n";
    std::cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

// Starting at buffer[i], which must be `open`, return the
// index of the matching `close` character.
// Returns a negative number if not found.
// Accounts for nested `open` and `close`.
int MatchingOpenCloseChar(const const_IWSubstring& buffer, int i, const char open,
                          const char close);
}  // namespace misc2
#endif  // MOLECULE_LIB_MISC2_H_
