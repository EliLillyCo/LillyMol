// Tool for seeing that two files contain the same structural information.
// Performs token-wise comparisons across lines.
// Tokens are the same if
//   They are the same as text
//   They can be interpreted as smiles and the unique_smiles are the same.
//   They can be interpreted as numeric, and the values are within a tolerance.
// NB: today this stops once a difference is found. TODO:ianwatson enable
// processing to continue;

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <unordered_set>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"

namespace same_structures {
using std::cerr;
namespace fs = std::filesystem;

void
Usage(int rc)
{
  cerr << "Performs chemically and numerically aware comparisons of files\n";
  cerr << " -D d1,d2     specify directories to be scanned for same files\n";
  cerr << " -r <tol>      relative tolerance for numeric values\n";
  cerr << " -h <nskip>    number of header records to skip\n";
  cerr << " -I            remove isotopes before comparing\n";
  cerr << " -s            slurp and sort the contents of the files\n";
  cerr << " -k            do NOT stop after encountering a difference, detect all diffs\n";
  cerr << " -S <col>      slupr and sort the file contents based on column <col>\n";
  cerr << " -i <sep>      name of input separator -i ,   -i tab\n";
  cerr << " -A ...        standard aromaticity options\n";
  cerr << " -E ...        standard elements options\n";
  cerr << " -K ...        standard smiles options\n";
  cerr << " -v            verbose output\n";
  ::exit(rc);
}

// Return 0 or 1 according to whether token1 is the same molecular
// graph as token2
// Return nullopt if the inputs do not look like graphs.
std::optional<int>
CompareAsGraphs(const IWString& token1, const IWString& token2) {
  const_IWSubstring smi1, h1;
  if (! token1.split(smi1, ':', h1)) {
    return std::nullopt;
  }
  const_IWSubstring smi2, h2;
  if (! token2.split(smi2, ':', h2)) {
    return std::nullopt;
  }

  if (h1 != h2) {
    return 0;
  }

  Molecule m1, m2;
  if (! m1.build_from_smiles(smi1) || ! m2.build_from_smiles(smi2)) {
    return std::nullopt;
  }

  return m1.unique_smiles() == m2.unique_smiles();
}

// If `token` looks like smiles:Hnn, then build `mol` from smiles
// and `hcount` from Hnn
int
BuildGraph(const IWString& token,
           IWString& usmi,
           IWString& hcount)
{
  IWString smi;
  if (! token.split(smi, ':', hcount) || smi.empty() || hcount.length() < 2) {
    return 0;
  }

  if (hcount[0] != 'H') {
    return 0;
  }

  hcount.remove_leading_chars(1);
  uint32_t notused;
  if (! hcount.numeric_value(notused)) {
    return 0;
  }

  Molecule m;
  if (! m.build_from_smiles(smi)) {
    return 0;
  }

  usmi = m.unique_smiles();

  return 1;
}

// Each file being compared is represented by this structure.
struct Afile {
  // File name
  std::string fname;
  // A stream for input from `fname`.
  iwstring_data_source input;

  // Each file will be asked to read a next_record, and will hold
  // the contents in `buffer`.
  IWString buffer;
  // When using nextword to scan through `buffer`, we need an int to
  // hold the state. Must be reset in `next_record`.
  int pos = 0;

  // Optionally, we can slurp the file - useful for sorting the contents
  std::vector<IWString> file_contents;
  // If we have slurped the file contents, we need to keep track of which
  // record has been copied to `buffer`.
  int buffer_is = -1;

  private:
    // There may be transformations performed on each line as it is read.
    void DoAnyRecordTransformations(IWString& line);
  public:

  // Opening and closing
  int Open(const char* file_name);
  int Close() {
    return input.do_close();
  }

  // Read the file contents into `file_contents`.
  int Slurp();

  // Sort `file_contents`. If `sort_column` is >= 0, sort on
  // the contents of that column.
  int SortFileContents(int sort_column);

  // Fetch the next record from `input` or `file_contents` into `buffer`.
  int NextRecord();

  int eof() const  {
    return input.eof();
  }

  const IWString& Buffer() const {
    return buffer;
  }

  // return into `token` the next word in `buffer` given that the
  // token separator is `sep`.
  int nextword(const_IWSubstring& token, char sep);
};

int
Afile::Open(const char * file_name) {
  if (! input.open(file_name)) {
    cerr << "Afile::Open:cannot open '" << file_name << "'\n";
    return 0;
  }

  fname = file_name;
  return 1;
}

void
Afile::DoAnyRecordTransformations(IWString& line) {
}

int
Afile::Slurp() {
  IWString line;
  while (input.next_record(line)) {
    DoAnyRecordTransformations(line);
    file_contents.emplace_back(line);
  }
  return file_contents.size();
}

int
Afile::SortFileContents(int sort_column) {
  if (sort_column < 0) {
    std::sort(file_contents.begin(), file_contents.end());
  } else {
    // Very inefficient.
    std::sort(file_contents.begin(), file_contents.end(), [sort_column](const IWString& line1, const IWString& line2) {
      const_IWSubstring token1 = line1.word(sort_column);
      const_IWSubstring token2 = line2.word(sort_column);
      return token1.strcmp(token2);
    });
  }
#define ECHO_SORT
#ifdef ECHO_SORT
  cerr << "File sorted\n";
  for (const IWString& line : file_contents) {
    cerr << line << '\n';
  }
#endif
  return file_contents.size();
}

int
Afile::NextRecord() {
  pos = 0;
  if (file_contents.empty()) {  // Reading line at a time.
    if (! input.next_record(buffer)) {
      return 0;
    }
    DoAnyRecordTransformations(buffer);
    return 1;
  }

  // From slurped data

  ++buffer_is;
  if (buffer_is >= static_cast<int>(file_contents.size())) {
    return 0;
  }
  buffer = file_contents[buffer_is];
  return 1;
}

int
Afile::nextword(const_IWSubstring& token, char sep) {
  return buffer.nextword(token, pos, sep);
}

// If comparisons are being done on a directory basis, a struct
// to hold that information.
struct Directory {
  std::string dirname;
  // Full path name of each item in the directory.
  std::vector<std::string> paths;

  // All the file names.
  std::unordered_set<std::string> file_names;

  int Initialise(IWString& dname, IWString& ignore_pattern);

  // Returns true if this directory has a file `fname`.
  bool HaveFile(const std::string& fname) const;

  // Given a file name `fname` that is in `dirname`, return the
  // full path name.
  std::string PathName(const std::string& fname);

  // All files in the directory.
  const std::unordered_set<std::string> FileNames() const {
    return file_names;
  }
};

// Scan `dname` for files, populating `paths` and `file_names`.
// Returns the number of files encountered.
// TODO: ianwatson enable exclusion conditions.
int
Directory::Initialise(IWString& dname, IWString& ignore_pattern)
{
  std::unique_ptr<re2::RE2> ignore_rx;
  if (ignore_pattern.length() > 0) {
    ignore_rx = std::make_unique<RE2>(ignore_pattern.null_terminated_chars());
  }

  dirname.assign(dname.data(), dname.length());
  paths.reserve(6);  // 6 is just a guess, adquate for gc3tk tests.
  for (const auto& file : fs::directory_iterator(dirname)) {
    std::string path_name = file.path();
    if (path_name[0] == '.') {
      continue;
    }

    const fs::file_status status = fs::status(path_name);
    // cerr << " file " << path_name << " dir ? " << (status.type() == fs::file_type::directory) << '\n';
    if (status.type() == fs::file_type::directory) {
      continue;
    }

    if (ignore_rx > 0 && RE2::PartialMatch(path_name, *ignore_rx)) {
      continue;
    }

    paths.emplace_back(std::move(path_name));

    std::filesystem::path tmp(file.path());
    file_names.insert(tmp.filename());
  }

  if (paths.size() == 0) {
    cerr << "Directory::Initialise:no files in " << dname << '\n';
  }

  return paths.size();
}

bool
Directory::HaveFile(const std::string& fname) const {
  const auto iter = file_names.find(fname);
  return iter != file_names.end();
}

std::string
Directory::PathName(const std::string& fname) {
  fs::path path_name;
  path_name.concat(dirname);
  path_name.append(fname);
  return path_name.string();
}

// Holds all information needed for a same_structures run.
struct Options {
  int verbose = 0;
  // Lines are read from each input file, this is the total in
  // each file.
  int lines_read = 0;

  int header_records_to_skip = 0;

  // The number of tokens interpreted as smiles.
  int molecules_read = 0;

  // If comparisons are done other than as text, how many such
  // comparisons are there.
  int same_via_numeric = 0;
  int same_via_smiles = 0;

  // Relative tolerance for numeric comparisons.
  double rtol = 0.0;

  // Either number of files or number of directories.
  int nfiles = 0;

  // We can either read each file line at a time, or slurp
  // the contents, sort them, and then serve line at a time.
  bool slurp_and_sort = false;

  // Sort based on the content of a specific column.
  int slurp_and_sort_by_column = -1;

  // If working via directory scanning, the number of directories
  // and Directory struct's for each.
  // From the -D option
  int ndir = 0;
  std::unique_ptr<Directory[]> directory;

  // If functioning by scanning files from the command line.
  std::vector<std::string> file_names;

  // The number of files absent in one or more directories.
  int missing_files = 0;

  int remove_isotopes = 0;

  // by default, we fail upon encountering an error.
  int detect_all_failures = 0;
  int records_with_differences = 0;

  // Input token separator.
  char sep = ' ';

  // Write all valid smiles here, one per file being handled.
  IWString_and_File_Descriptor stream_for_smiles[2];

  // As `nfiles` files are each processed, arrays to hold the values
  // from each file.
  std::unique_ptr<const_IWSubstring[]> token;
  std::unique_ptr<double[]> numeric;
  std::unique_ptr<Molecule[]> mol;
  std::unique_ptr<IWString[]> usmi;
  std::unique_ptr<IWString[]> hcount;

  private:

  int AllStructuresTheSame(Afile* files);
  int AllStructuresTheSameToken(const Afile* files);

  int CompareNumbers(const Afile* files);
  int AboutTheSame(double v1, double v2) const;

  int CompareFiles(const std::vector<std::string>& fnames);
  int CompareFilesInner(Afile* files);
  int CompareDirectories();

  public:

  ~Options();

  // Returns true if all files are EOF.
  int AllFilesEof(const Afile* files) const;

  int Initialise(Command_Line& cl);

  // Once Initialise'd, and files/directories specified, process
  // that information.
  int Process();

  void Report(std::ostream& output) const;
};

Options::~Options() {
}

int
Options::Initialise(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('E')) {
    if (! process_elements(cl, verbose, 'E')) {
      cerr << "Cannot discern elements, -E\n";
      Usage(1);
    }
  }

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process standard aromaticity options\n";
    return 0;
  }

  if (cl.option_present('K')) {
    if (! process_standard_smiles_options(cl, verbose, 'K')) {
      cerr << "Cannot initialise smiles options\n";
      return 0;
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', rtol) || rtol <= 0.0) {
      cerr << "The relative tolerance value (-r) must be a valid relative tolerance\n";
      return 0;
    }
    if (verbose) {
      cerr << "Numeric values the same if rtol within " << rtol << '\n';
    }
  }

  if (cl.option_present('h')) {
    if (! cl.value('h', header_records_to_skip) || header_records_to_skip <= 0) {
      cerr << "Header records to skip must be a whole +ve number\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will skip the first " <<header_records_to_skip << " records of each file\n";
    }
  }

  if (cl.option_present('i')) {
    IWString i = cl.option_value('i');
    if (! char_name_to_char(i)) {
      cerr << "Unrecognised input separator '" << i << "'\n";
      return 0;
    }
    sep = i[0];
  }

  if (cl.option_present('s') &&
      cl.option_present('S')) {
    cerr << "Cannot use both -s and -S options\n";
    return 0;
  }

  if (cl.option_present('s')) {
    slurp_and_sort = true;
  }

  if (cl.option_present('S')) {
    if (! cl.value('S', slurp_and_sort_by_column) || slurp_and_sort_by_column < 1) {
      cerr << "Invalid slurp and sort column (-S)\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will slurp file contents and sort on column " << slurp_and_sort_by_column << '\n';
    }
    --slurp_and_sort_by_column;
  }

  if (cl.option_present('I')) {
    remove_isotopes = 1;
    if (verbose) {
      cerr << "Will remove isotopes\n";
    }
  }

  if (cl.option_present('k')) {
    detect_all_failures = 1;
    if (verbose) {
      cerr << "Will detect all failures\n";
    }
  }

  if (cl.option_present('W')) {
    const IWString fname = cl.string_value('W');

    for (int i = 0; i < 2; ++i) {
      IWString my_fname(fname);
      my_fname << i << ".smi";
      if (! stream_for_smiles[i].open(my_fname.null_terminated_chars())) {
        cerr << "Cannot open smiles stream '" << my_fname << "'\n";
        return 0;
      }
    }
    if (verbose) {
      cerr << "All smiles writen to files starting with " << fname << "'\n";
    }
  }

  if (cl.option_present('D')) {
    IWString ignore_pattern;
    if (cl.option_present('X')) {
      ignore_pattern = cl.string_value('X');
      if (verbose) {
        cerr << "Will skip files matching " << ignore_pattern << '\n';
      }
    }

    IWString dir_names = cl.string_value('D');
    ndir = dir_names.nwords(',');
    if (ndir < 2) {
      cerr << "The number of directories to scan must be > 1\n";
      return 0;
    }
    directory = std::make_unique<Directory[]>(ndir);
    int i = 0;
    IWString dir;
    for (int dir_index = 0; dir_names.nextword(dir, i, ','); ++dir_index ) {
      if (! directory[dir_index].Initialise(dir, ignore_pattern)) {
        cerr << "Cannot initialise directory '" << dir << "'\n";
        return 0;
      }
    }
    nfiles = ndir;
  } else if (cl.empty()) {
    cerr << "No input files specified\n";
    Usage(1);
  } else {
    nfiles = cl.number_elements();
    if (nfiles < 2) {
      cerr << "Must specify multiple input files, only " << nfiles << " specified\n";
      Usage(1);
    }
    for (const char * fname : cl) {
      file_names.emplace_back(fname);
    }
  }

  token = std::make_unique<const_IWSubstring[]>(nfiles);
  numeric = std::make_unique<double[]>(nfiles);
  mol = std::make_unique<Molecule[]>(nfiles);
  usmi = std::make_unique<IWString[]>(nfiles);
  hcount = std::make_unique<IWString[]>(nfiles);

  return 1;
}

int
Options::CompareFiles(const std::vector<std::string>& fnames) {
  const int nfiles = fnames.size();
  Afile* files = new Afile[nfiles]; std::unique_ptr<Afile[]> free(files);
  int empty_files = 0;
  for (int i = 0; i < nfiles; ++i) {
    if (! files[i].Open(fnames[i].c_str())) {
      cerr << "Cannot open '" << fnames[i] << "'\n";
      return 1;
    }
    if (slurp_and_sort || slurp_and_sort_by_column >= 0) {
      if (files[i].Slurp() == 0) {
        ++empty_files;
      } else {
        files[i].SortFileContents(slurp_and_sort_by_column);
      }
    }
  }

  if (empty_files == nfiles) {
    return 1;
  }

  return CompareFilesInner(files);
}

int
Options::CompareFilesInner(Afile* files)
{
  if (header_records_to_skip > 0) {
    for (int i = 0; i < nfiles; ++i) {
      for (int j = 0; j < header_records_to_skip; ++j) {
        if (! files[i].NextRecord()) {
          cerr << "Premature eof on file " << i << '\n';
          return 0;
        }
      }
    }
  }

  if (verbose > 1) {
    cerr << "Opened " << nfiles << " files for comparison\n";
  }

  for (int current_line_number = 0; ; ++current_line_number, ++lines_read) {
    int got_data = 0;
    for (int i = 0; i < nfiles; ++i) {
      if (files[i].NextRecord()) {
        ++got_data;
      }
    }

    if (verbose > 1) {
      cerr << current_line_number << " lines read, " << got_data << " got data\n";
    }
    // All files are at eof, good.
    if (got_data == 0) {
      return 1;
    }

    if (got_data != nfiles) {
      cerr << "Premature EOF, got data on " << got_data << " of " << nfiles << " files\n";
      Report(cerr);
      return 0;
    }

    if (! AllStructuresTheSame(files)) {
      cerr << "Differences detected on line " << current_line_number << '\n';
      ++records_with_differences;
      if (detect_all_failures) {
        continue;
      }
      return 0;
    }
  }

  Report(cerr);

  return 1;
}

// A set of directories has been loaded into `directory`. Scan all files in
// those directories to look for equivalence.
// Note that file name scans are done relative to the first directory, so if
// there are new files in some of the subsequent directories, those are not
// examined. Bug or feature depending... 
int
Options::CompareDirectories()
{
  const std::unordered_set<std::string> file_names = directory[0].FileNames();
  for (const std::string& fname : file_names) {
    std::vector<std::string> path_names;
    if (verbose > 1) {
      cerr << "Looking for " << fname << " across " << ndir << " directories\n";
    }
    for (int i = 0; i < ndir; ++i) {
      if (directory[i].HaveFile(fname)) {
        path_names.emplace_back(directory[i].PathName(fname));
      } else {
        cerr << "Options::CompareDirectories:directory " << directory[i].dirname << " does not have " << fname << '\n';
      }
    }

    if (static_cast<int>(path_names.size()) != ndir) {
      if (verbose > 1) {
        cerr << "Options::CompareDirectories: " << path_names.size() << " of " << ndir << " directories had " << fname << '\n';
      }
      ++missing_files;
      return 0;
    }

    if (! CompareFiles(path_names)) {
      cerr << "Error processing file " << fname << " across " << ndir << " directories\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Process()
{
  if (ndir > 0) {
    return CompareDirectories();
  } else {
    return CompareFiles(file_names);
  }
}

// Each file in `files` has a buffer read from file. Determine
// whether or not those buffers are all the same.
// First try comparison as text, then other means.
int
Options::AllStructuresTheSame(Afile* files)
{
  bool all_records_identical = true;
  for (int i = 1; i < nfiles; ++i) {
    if (files[i - 1].Buffer() != files[i].Buffer()) {
      all_records_identical = false;
      break;
    }
  }

  if (all_records_identical) {
    return 1;
  }

  // Lines are not identical, compare token by token.

#ifdef EEBUG_SAME_STRUCTURES
  for (int i = 0; i < nfiles; ++i) {
    cerr << "Comparing as token " << files[i].Buffer() << '\n';
  }
#endif

  while (true) {
    int got_token = 0;
    for (int i = 0; i < nfiles; ++i) {
      if (files[i].nextword(token[i], sep)) {
        ++got_token;
      }
    }

    if (got_token == 0) {  // end of line.
        return 1;
    }

    if (got_token != nfiles) {
      cerr << "Options::AllStructuresTheSame:cannot extract token from all files\n";
      return 0;
    }

    if (! AllStructuresTheSameToken(files)) {
      cerr << "Options::AllStructuresTheSame: mismatch on line " << lines_read << '\n';
      for (int i = 0; i < nfiles; ++i) {
        cerr << files[i].Buffer() << "    :  " << token[i] << '\n';
      }
      return 0;
    }
  }

  return 1;
}

int
Options::AllFilesEof(const Afile* files) const {
  for (int i = 0; i < nfiles; ++i) {
    if (! files[i].eof()) {
      return 0;
    }
  }

  return 1;
}

// The `token` array has been populated from the files in `files`.
// Can those tokens be considered the same.
int
Options::AllStructuresTheSameToken(const Afile* files)
{
  int same_text = 1;
  for (int i = 1; i < nfiles; ++i) {
    if (token[i-1] == token[i]) {
      ++same_text;
    }
  }

  if (same_text == nfiles) {
    return 1;
  }

  int is_numeric = 0;
  for (int i = 0; i < nfiles; ++i) {
    if (token[i].numeric_value(numeric[i])) {
      ++is_numeric;
    }
  }

  if (is_numeric == nfiles) {
    return CompareNumbers(files);
  }

  int is_graph = 0;
  for (int i = 0; i < nfiles; ++i) {
    if (BuildGraph(token[i], usmi[i], hcount[i])) {
      ++is_graph;
    }
  }

  if (is_graph == nfiles) {
    bool all_graphs_the_same = true;
    for (int i = 1; i < nfiles; ++i) {
      if (hcount[i-1] != hcount[i] || usmi[i-1] != usmi[i]) {
        all_graphs_the_same = false;
        break;
      }
    }

    if (all_graphs_the_same) {
      return 1;
    }
  }

  int valid_molecule = 0;
  for (int i = 0; i < nfiles; ++i) {
    if (mol[i].build_from_smiles(token[i])) {
      ++molecules_read;
      ++valid_molecule;
      if (stream_for_smiles[i].active()) {
        stream_for_smiles[i] << token[i] << '\n';
      }
    }
  }

  if (valid_molecule == nfiles) {  // All valid smiles
  } else if (valid_molecule == 0) {  //  None valid smiles
    return 0;
  } else {
    cerr << "Some tokens valid as smiles\n";
    for (int i = 0; i < nfiles; ++i) {
      if (mol[i].build_from_smiles(token[i])) {
        cerr << token[i] << " OK\n";
      } else {
        cerr << token[i] << " INVALID\n";
      }
    }
    return 0;
  }

  if (remove_isotopes) {
    for (int i = 0; i < nfiles; ++i) {
      mol[i].transform_to_non_isotopic_form();
    }
  }

  bool all_smiles_the_same = true;
  for (int i = 0; i < nfiles; ++i) {
    usmi[i] = mol[i].unique_smiles();
    // cerr << "generated " << usmi[i] << '\n';
    if (i > 0 && usmi[i-1] != usmi[i]) {
      cerr << "Smiles mismatch " << usmi[i-1] << " vs " << usmi[i] << '\n';
      cerr << " aromatic atom count " << mol[i-1].aromatic_atom_count() << " and " << mol[i].aromatic_atom_count() << '\n';
      all_smiles_the_same = false;
    }
  }

  if (! all_smiles_the_same) {
    return 0;
  }

  ++same_via_smiles;
  return 1;
}

// The `numeric` array has been filled with a value from each file in `files`.
// Return true if those numbers are equivalent.
int
Options::CompareNumbers(const Afile* files)
{
  bool all_numbers_equal = true;
  for (int i = 1; i < nfiles; ++i) {
    if (numeric[i-1] != numeric[i]) {
      all_numbers_equal = false;
    }
  }

  if (all_numbers_equal) {
    return 1;
  }

  ++same_via_numeric;
  for (int i = 0; i < nfiles; ++i) {
    for (int j = i + 1; j < nfiles; ++j) {
      if (! AboutTheSame(numeric[i], numeric[j])) {
        cerr << "Options::CompareNumbers:not the same " << numeric[i] << " " << numeric[j] << " line " << lines_read << '\n';
        return 0;
      }
    }
  }

  return 1;
}

// Return true if `v1` and `v2` are within `rtol` of each other.
int
Options::AboutTheSame(double v1, double v2) const {
  if (v1 == 0.0) {
    return std::abs(v2) < rtol;
  }
  if (v2 == 0.0) {
    return std::abs(v1) < rtol;
  }

  const double diff = std::abs(v1 - v2);
  const double denominator = std::max(v1, v2);

  return (diff / denominator) < rtol;
}

void
Options::Report(std::ostream& output) const {
  output << "Read " <<  lines_read << " lines, across " << nfiles << " files " << molecules_read << " molecules\n";
  cerr << same_via_smiles << " same via unique smiles\n";
  cerr << same_via_numeric << " same via rtol " << rtol << '\n';
  if (verbose && detect_all_failures) {
    cerr << records_with_differences << " of " << lines_read << " lines contained differences\n";
  }
}

int
SameStructures(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:g:E:K:r:h:i:D:sS:X:IkW:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  set_display_smiles_interpretation_error_messages(0);

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  if (! options.Process()) {
    cerr << "Error during processing\n";
    return 1;
  }

  options.Report(cerr);

  return 0;
}

}  // namespace same_structures

int
main(int argc, char** argv)
{
  int rc = same_structures::SameStructures(argc, argv);

  return rc;
}
