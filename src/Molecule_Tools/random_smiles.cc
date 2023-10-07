#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_map>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define ISTREAM_AND_TYPE_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int molecules_read = 0;

static extending_resizable_array<int> permutations_formed;

static Chemical_Standardisation chemical_standardisation;

static int number_permutations = 10;

static int max_iterations = 100;

static int append_permutation_number_to_name = 0;

static int echo_initial_molecule = 0;

static std::random_device rd;

static Atom_Typing_Specification atom_typing_specification;

static int no_atom_with_type = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << "  -n <number>    number of random smiles to make for each input molecule\n";
  cerr << "  -I <number>    try a maximum of <number> times to find the permuations\n";
  cerr << "  -a             append permutation number to name\n";
  cerr << "  -e             echo initial molecule\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -E autocreate  automatically create new elements when encountered\n";
  (void)display_standard_aromaticity_options(cerr);
  (void)display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -K ...         standard smiles control options, enter '-K help' for info\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

class Molecule_and_Atom_Type : public Molecule
{
 private:
  uint32_t* _atype;

 public:
  Molecule_and_Atom_Type();
  ~Molecule_and_Atom_Type();

  int assign_atom_types(Atom_Typing_Specification&);

  const uint32_t*
  atype() const
  {
    return _atype;
  }
};

Molecule_and_Atom_Type::Molecule_and_Atom_Type()
{
  _atype = nullptr;
}

Molecule_and_Atom_Type::~Molecule_and_Atom_Type()
{
  if (nullptr != _atype) {
    delete[] _atype;
  }
}

int
Molecule_and_Atom_Type::assign_atom_types(Atom_Typing_Specification& ats)
{
  if (nullptr != _atype) {  // should not happen
    delete[] _atype;
  }

  const int matoms = this->natoms();
  assert(matoms > 0);

  _atype = new uint32_t[matoms];

  if (!ats.assign_atom_types(*this, _atype)) {
    cerr << "Molecule_and_Atom_Type::assign_atom_types:cannot assign atom types\n";
    return 0;
  }

  return 1;
}

class Smiles_and_Index
{
 private:
  IWString _smiles;

  std::unordered_map<int, int> _atype_to_smiles;  // atom type to index in _smiles

 public:
  void build(Molecule_and_Atom_Type& m, int* tmp, const std::unordered_set<int>& todo);

  int place_smiles(const int atype, const int centre, const int width,
                   IWString& output) const;
};

void
Smiles_and_Index::build(Molecule_and_Atom_Type& m, int* tmp,
                        const std::unordered_set<int>& todo)
{
  m.random_smiles();

  const resizable_array<atom_number_t>& sorder = m.atom_order_in_smiles();

  const int matoms = m.natoms();

  const uint32_t* atype = m.atype();

  for (int i = 0; i < matoms; ++i) {
    tmp[sorder[i]] = i;

    const auto f = todo.find(atype[i]);
    if (f == todo.end()) {
      continue;
    }

    _atype_to_smiles[atype[i]] = sorder[i];
  }

  for (int i = 0; i < matoms; ++i) {
    const atom_number_t j = tmp[i];  // which atom is the I'th atom in the smiles
    const Element& e = m.element(j);

    if (m.is_aromatic(j)) {
      _smiles << e.aromatic_symbol();
    } else {
      _smiles << e.symbol();
    }
  }

  return;
}

int
Smiles_and_Index::place_smiles(const int atype, const int centre, const int width,
                               IWString& output) const
{
  const auto f = _atype_to_smiles.find(atype);

  if (f == _atype_to_smiles.end())  // molecule does not contain the atom type
  {
    if (verbose > 2) {
      cerr << "Smiles_and_Index::place_smiles:no atom with atom type " << atype << endl;
    }
    no_atom_with_type++;
    return 0;
  }

  const int j = f->second;

  cerr << "Atype " << atype << " found at position " << j << " " << _smiles[j] << endl;

  if (j - width < 0 || j + width >= _smiles.length()) {
    return 0;
  }

  for (int i = -width; i < width; ++i) {
    output[centre + i] = _smiles[j + i];
  }

  return 1;
}

class Set_of_Random_Smiles
{
 private:
  std::mt19937_64 _rng;

  IWString _name;

  int _n;
  Smiles_and_Index* _m;

 public:
  Set_of_Random_Smiles();
  ~Set_of_Random_Smiles();

  int build(Molecule_and_Atom_Type& m, const int n, const std::unordered_set<int>& todo);

  int place_smiles(const int atype, const int centre, const int width, IWString& output);

  const IWString&
  name() const
  {
    return _name;
  }
};

Set_of_Random_Smiles::Set_of_Random_Smiles() : _rng(rd())
{
  _n = 0;
  _m = nullptr;
}

Set_of_Random_Smiles::~Set_of_Random_Smiles()
{
  if (nullptr != _m) {
    delete[] _m;
  }
}

int
Set_of_Random_Smiles::build(Molecule_and_Atom_Type& m, const int ncopy,
                            const std::unordered_set<int>& todo)
{
  if (nullptr != _m) {  // should not happen
    delete[] _m;
  }

  _m = new Smiles_and_Index[ncopy];
  _n = ncopy;

  int* tmp = new int[m.natoms()];
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < ncopy; ++i) {
    _m[i].build(m, tmp, todo);
  }

  _name = m.name();

  return 1;
}

int
Set_of_Random_Smiles::place_smiles(const int atype, const int centre, const int width,
                                   IWString& output)
{
  std::uniform_int_distribution<int> u(0, _n - 1);

  int s = u(_rng);

  for (int i = s; i < _n; ++i) {
    if (_m[i].place_smiles(atype, centre, width, output)) {
      return 1;
    }
  }

  for (int i = 0; i < s; ++i) {
    if (_m[i].place_smiles(atype, centre, width, output)) {
      return 1;
    }
  }

  return 0;
}

class For_Convolution
{
 private:
  int _n;
  Set_of_Random_Smiles* _s;

  resizable_array<int> _bits;  // which bits are we producing - centre of the window

  int _frame_width;  // how many atoms either side of the mid-point

  int _replicates;  // within each Set_of_Random_Smiles how many duplicates should they
                    // build

  std::unordered_set<int> _todo;  // the actual atom types we will process

  int _replicate_each_type;  // how many instances of a given type do we ask for

  //  private functions

  void _display_for_convolution_options(const char flag, std::ostream& output) const;
  int _build(Set_of_Random_Smiles& sors, Molecule& m,
             std::unordered_map<int, int>& atype_count);
  int _generate_row(const int r, IWString& buffer);

 public:
  For_Convolution();
  ~For_Convolution();

  int debug_print(std::ostream&) const;

  int build(Command_Line& cl, char flag);

  int run(std::ostream& output);
};

For_Convolution::For_Convolution()
{
  _n = 0;
  _s = nullptr;
  _frame_width = 4;
  _replicates = 10;
  _replicate_each_type = 1;
}

For_Convolution::~For_Convolution()
{
  if (nullptr != _s) {
    delete[] _s;
  }
}

int
For_Convolution::debug_print(std::ostream& output) const
{
  output << "For_Convolution::debug_print: " << _frame_width << " frame width, "
         << _replicate_each_type << " replicates each atom type\n";
  output << _todo.size() << " atom types defined for output\n";
  output << _replicates << " replicates of each of " << _n << " smiles\n";

  return 1;
}

static void
count_atom_types(const Molecule_and_Atom_Type& m,
                 std::unordered_map<int, int>& atype_count)
{
  const int matoms = m.natoms();

  const uint32_t* atype = m.atype();

  for (int i = 0; i < matoms; ++i) {
    atype_count[atype[i]]++;
  }

  return;
}

int
For_Convolution::build(Command_Line& cl, char flag)
{
  int nbits = 0;

  const_IWSubstring s;
  for (int i = 0; cl.value(flag, s, i); ++i) {
    if (s.starts_with("width=")) {
      s.remove_leading_chars(6);
      if (!s.numeric_value(_frame_width) || _frame_width < 1) {
        cerr << "The frame width must be positive '" << s << "' invalid\n";
        return 0;
      }
    } else if (s.starts_with("N=")) {
      s.remove_leading_chars(2);
      if (!s.numeric_value(_replicates) || _replicates < 1) {
        cerr << "For_Convolution::build:the number of replicates must ba positive '" << s
             << "' invalid\n";
        return 0;
      }
    } else if (s.starts_with("nbits=")) {
      s.remove_leading_chars(6);
      if (!s.numeric_value(nbits) || nbits < 2) {
        cerr << "For_Convolution::build:the number of bits must be a +ve whole number '"
             << s << "' invalid\n";
        return 0;
      }
    } else if (s.starts_with("each=")) {
      s.remove_leading_chars(5);
      if (!s.numeric_value(_replicate_each_type) || _replicate_each_type < 1) {
        cerr << "For_Convolution::build:the number of replicates each type must be "
                "positive\n";
        return 0;
      }
    } else if ("help" == s) {
      _display_for_convolution_options(flag, cerr);
      exit(0);
    } else {
      cerr << "For_Convolution::build:unrecognised directive\n";
      return 0;
    }
  }

  if (0 == _frame_width || 0 == _replicates || 0 == nbits) {
    cerr << "For_Convolution::build:not completely initialised\n";
    return 0;
  }

  if (verbose) {
    debug_print(cerr);
  }

  if (cl.empty()) {
    cerr << "For_Convolution::build:no input file\n";
    return 0;
  }

  const char* fname = cl[0];

  data_source_and_type<Molecule_and_Atom_Type> input(FILE_TYPE_SMI, fname);

  if (!input.good()) {
    cerr << "For_Convolution::build:cannot open '" << fname << "'\n";
    return 0;
  }

  _n = input.molecules_remaining();

  if (0 == _n) {
    cerr << "For_Convolution:;build:empty file\n";
    return 0;
  }

  _s = new Set_of_Random_Smiles[_n];

  Molecule_and_Atom_Type* m;

  resizable_array_p<Molecule_and_Atom_Type> molecules;

  for (int i = 0; nullptr != (m = input.next_molecule()); ++i) {
    if (!m->assign_atom_types(atom_typing_specification)) {
      cerr << "For_Convolution::build:fatal error processing " << m->smiles() << ' '
           << m->name() << endl;
      delete m;
      ;
      return 0;
    }

    molecules.add(m);
  }

  assert(molecules.number_elements() == _n);

  std::unordered_map<int, int>
      atype_count;  // how frequently do we find a given atom type

  for (int i = 0; i < _n; ++i) {
    count_atom_types(*molecules[i], atype_count);
  }

  if (verbose) {
    cerr << "Found " << atype_count.size() << " different atom types across " << _n
         << " molecules\n";
  }

  const int natypes = atype_count.size();

  std::pair<int, int>* pairs = new std::pair<int, int>[natypes];

  int ndx = 0;
  for (auto f : atype_count) {
    pairs[ndx].first = f.first;
    pairs[ndx].second = f.second;
    ndx++;
  }

  std::sort(pairs, pairs + natypes,
            [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
              return p1.second > p2.second;
            });

  if (nbits > natypes) {
    cerr << "For_Convolution::build:only found " << natypes
         << " different atom types, cannot do " << nbits << " bits, truncating\n";
    nbits = natypes;
  }

  _bits.resize(nbits);

  assert(nbits <= natypes);

  for (int i = 0; i < nbits; ++i) {
    const auto x = pairs[i].first;  // one of the selected atom types
    cerr << " atype " << x << " found " << pairs[i].second << " times\n";

    _bits.add(x);
    _todo.insert(x);
  }

  for (int i = 0; i < _n; ++i) {
    if (!_s[i].build(*molecules[i], _replicates, _todo)) {
      return 0;
    }
  }

  return 1;
}

void
For_Convolution::_display_for_convolution_options(const char flag,
                                                  std::ostream& output) const
{
  output << " -" << flag << " width=w    frame width for each smiles fragment\n";
  output << " -" << flag
         << " N=n        number of random smiles replicates (default 10)\n";
  output << " -" << flag << " nbits=n    number of atom types to process\n";
  return;
}

#ifdef NO_LONGER_USED_MMMMMMMMM2
int
For_Convolution::_build(Set_of_Random_Smiles& sors, Molecule& m,
                        std::unordered_map<int, int>& atype_count)
{
  const int matoms = m.natoms();

  int* atype = new int[matoms];
  std::unique_ptr<int[]> free_atype(atype);

  if (!atom_typing_specification.assign_atom_types(m, atype)) {
    cerr << "For_Convolution::_build:cannot assign atom types to " << m.smiles() << ' '
         << m.name() << endl;
    return 0;
  }

  m.set_isotopes(atype);

  for (int i = 0; i < matoms; ++i) {
    atype_count[atype[i]]++;
  }

  return sors.build(m, _replicates);
}
#endif

int
For_Convolution::run(std::ostream& output)
{
  IWString buffer;

  const int row_size =
      _replicate_each_type * _bits.number_elements() * (_frame_width + _frame_width + 1);

  buffer.resize(row_size);

  for (int i = 0; i < _n; ++i)  // the number of rows we will generate
  {
    buffer.resize_keep_storage(0);
    buffer.extend(row_size, ' ');

    _generate_row(i, buffer);

    //  output << _s[i].name() << ' ' <<  buffer << '\n';
    output << buffer << '\n';

    buffer.resize_keep_storage(0);
  }

  return output.good();
}

int
For_Convolution::_generate_row(const int r, IWString& buffer)
{
  const int two_width = _frame_width + _frame_width;

  cerr << "Will generate " << _replicate_each_type << " replicates of each of "
       << _todo.size() << " atom types\n";

  int centre = _frame_width + 1;

  for (auto t : _bits)  // each of the atom types we are going to do
  {
    for (int j = 0; j < _replicate_each_type; ++j) {
      cerr << "placing smiles at " << centre << ", buf " << buffer.size() << endl;

      _s[r].place_smiles(t, centre, _frame_width, buffer);
      centre += two_width;
    }
  }

  return 1;
}

static int
write_the_permutation(const IWString& mname, const IWString& rsmi, int permutation_number,
                      std::ostream& output)
{
  output << rsmi << ' ' << mname;

  if (append_permutation_number_to_name) {
    output << " PERMUTATION " << permutation_number;
  }

  output << endl;

  return output.good();
}

static void
preprocess(Molecule& m)
{
  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
make_random_smiles(Molecule& m, std::ostream& output)
{
  preprocess(m);

  IWString_STL_Hash_Set generated;

  if (echo_initial_molecule) {
    output << m.smiles() << ' ' << m.name() << endl;
  }

  int permutations_generated = 0;

  int iterations = 0;

  while (permutations_generated < number_permutations) {
    IWString rsmi = m.random_smiles();
    //  cerr << "rsmi " << rsmi << endl;

    if (!generated.contains(rsmi)) {
      generated.insert(rsmi);
      permutations_generated++;

      write_the_permutation(m.name(), rsmi, permutations_generated, output);

      if (permutations_generated >= number_permutations) {
        break;
      }
    }

    iterations++;

    if (iterations >= max_iterations) {
      if (verbose) {
        cerr << "After " << iterations << " iterations, only have "
             << permutations_generated << " permutations '" << m.name() << "'\n";
      }
      break;
    }
  }

  permutations_formed[permutations_generated]++;  // update global statistic

  return output.good();
}

static int
make_random_smiles(data_source_and_type<Molecule>& input, std::ostream& output)
{
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (!make_random_smiles(*m, output)) {
      return 0;
    }
  }

  return 1;
}

static int
make_random_smiles(const char* fname, FileType input_type, std::ostream& output)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return make_random_smiles(input, output);
}

static int
do_convolution(Command_Line& cl, std::ostream& output)
{
  For_Convolution fc;

  if (!fc.build(cl, 'C')) {
    cerr << "Cannot initialise for convolution object\n";
    return 0;
  }

  return fc.run(output);
}

static int
make_random_smiles(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:g:i:n:I:eaK:P:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(3);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose > 1, 'A')) {
    cerr << "Cannot process standard aromaticity options (-A)\n";
    usage(4);
  }

  if (!process_elements(cl, verbose > 1, 'E')) {
    cerr << "Cannot initialise elements (-E option)\n";
    usage(2);
  }

  if (!process_standard_smiles_options(cl, verbose)) {
    usage(4);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      usage(6);
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', number_permutations) || number_permutations < 1) {
      cerr << "The number of permutations (-n) option must be followed by a whole "
              "positive number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will make " << number_permutations << " random smiles of each molecule\n";
    }
  }

  if (cl.option_present('I')) {
    if (!cl.value('I', max_iterations) || max_iterations < number_permutations) {
      cerr << "The maximum number of iterations (-I) option must be a whole number >= "
           << number_permutations << endl;
      usage(6);
    }

    if (verbose) {
      cerr << "Will try a maximum of " << max_iterations
           << " iterations for each molecule\n";
    }
  } else {
    max_iterations = 10 * number_permutations;
  }

  if (cl.option_present('a')) {
    append_permutation_number_to_name = 1;

    if (verbose) {
      cerr << "permutation number appended to name\n";
    }
  }

  if (cl.option_present('e')) {
    echo_initial_molecule = 1;

    if (verbose) {
      cerr << "Will echo initial smiles\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Invalid atom type '" << p << "'\n";
      return 1;
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('C')) {
    if (!atom_typing_specification.active()) {
      atom_typing_specification.build("UST:AY");
      if (verbose) {
        cerr << "Default atom typing 'UST:AY'\n";
      }
    }

    if (!do_convolution(cl, std::cout)) {
      return 1;
    }

    return 0;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!make_random_smiles(cl[i], input_type, std::cout)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";

    for (int i = 0; i < permutations_formed.number_elements(); i++) {
      if (permutations_formed[i]) {
        cerr << permutations_formed[i] << " molecules yielded " << i
             << " random smiles\n";
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = make_random_smiles(argc, argv);

  return rc;
}
