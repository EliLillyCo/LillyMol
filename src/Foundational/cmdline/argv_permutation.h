#ifndef FOUNDATIONAL_CMDLINE_ARGV_PERMUTATION_H_
#define FOUNDATIONAL_CMDLINE_ARGV_PERMUTATION_H_

#include <cstring>
#include <vector>

namespace cmdline_internal {

// Move options and their values before operands, preserving both orders.
// Only argv pointers move; getopt and Command_Line retain the original strings.
inline int PermuteArguments(int argc, char** argv, const char* options) {
  // Explicit ordering modes should not acquire default GNU permutation.
  if (options[0] == '+' || options[0] == '-') {
    return argc;
  }
  std::vector<char*> flags;
  std::vector<char*> operands;
  bool missing_value = false;
  int end_of_options = argc;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--") == 0) {
      end_of_options = i;
      break;
    }
    if (argv[i][0] != '-' || argv[i][1] == '\0') {
      operands.push_back(argv[i]);
      continue;
    }
    flags.push_back(argv[i]);
    // A cluster ends when an option takes a value. That value may be attached
    // or be the next token, even if the next token starts with a dash.
    for (const char* flag = argv[i] + 1; *flag; ++flag) {
      const char* spec = std::strchr(options, *flag);
      if (*flag == ':' || spec == nullptr || spec[1] != ':') {
        continue;
      }
      // Double-colon optional values, if supported, must be attached.
      if (flag[1] == '\0' && spec[2] != ':' && i + 1 < argc) {
        flags.push_back(argv[++i]);
      } else if (flag[1] == '\0' && spec[2] != ':') {
        missing_value = true;
      }
      break;
    }
  }
  int dest = 1;
  for (char* flag : flags) {
    argv[dest++] = flag;
  }
  // Place the delimiter before all operands, including those originally before
  // it, so getopt cannot interpret any operand as an option.
  if (end_of_options < argc) {
    argv[dest++] = argv[end_of_options];
  }
  for (char* operand : operands) {
    argv[dest++] = operand;
  }
  // The suffix after -- was never moved and cannot be overwritten above.
  // A missing value must not consume an operand moved behind the options.
  return missing_value ? static_cast<int>(flags.size()) + 1 : argc;
}

}  // namespace cmdline_internal
#endif  // FOUNDATIONAL_CMDLINE_ARGV_PERMUTATION_H_
