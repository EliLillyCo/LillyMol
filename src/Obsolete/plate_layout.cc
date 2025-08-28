/*
  We are trying to select a set of plates which promises the best
  compromise between diversity and score.
  We need:
    a file with 'plate score'
    an intra plate diversity file 'plate diversity'
    an inter plate diversity file 'plate diversity'
*/

#include <stdlib.h>

#include <iostream>
#include <memory>
#include <random>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

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
  // clang-format off
  cerr << "Diversity vs score plate layout optimisation\n";
  cerr << " -R <fname>     file with intRa plate diversity scores\n";
  cerr << " -E <fname>     file with intEr plate diversity scores\n";
  cerr << " -S <fname>     file with plate and score values\n";
  cerr << " -h <number>    size of diversity histograms (usually 10)\n";
  cerr << " -H <number>    size of score histograms (usually 10)\n";
  cerr << " -p <number>    the number of plates to select\n";
  cerr << " -n <number>    number of random configurations to create\n";
  cerr << " -V <fname>     evaluate existing layout in <fname>\n";
  cerr << " -P <options>   spread options\n";
  cerr << " -P random      choose first plate at random\n";
  cerr << " -P best        choose the plate with the highest score first\n";
  cerr << " -b <float>     bias towards diversity\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int verbose = 0;

typedef float score_t;

/*
  We need to know the histogram sizes
*/

static int diversity_histogram_size = 0;
static int score_histogram_size = 0;

/*
  When comparing scores and diversities, we may not want to consider
  the entire distribution
*/

static int score_histogram_consider = 6;
static int diversity_histogram_consider = 3;

/*
  When doing close comparisons between diversity and distances, we can
  bias the trade-off towards the model
*/

static float model_bias = 0.0;

/*
  When doing a spread run, how do we choose the first point
*/

static int choose_first_plate_randomly = 0;

/*
  Or we may choose the plate with the best score
*/

static int choose_first_plate_best = 0;

class Set_of_Values : public iwaray<int>
{
 private:
 public:
  Set_of_Values& operator=(const Set_of_Values&);
  Set_of_Values& operator+=(const Set_of_Values&);
  Set_of_Values& operator-=(const Set_of_Values&);

  int operator==(const Set_of_Values&) const;
  int operator!=(const Set_of_Values&) const;
};

class Score : public Set_of_Values
{
 private:
 public:
  Score();

  int read_score(const const_IWSubstring&);

  Score& operator=(const Score&);

  int operator!=(const Score&) const;

  int operator>(const Score&) const;
};

Score::Score()
{
  assert(score_histogram_size > 0);

  resize(score_histogram_size);

  for (int i = 0; i < _number_elements; i++) {
    _things[i] = 0;
  }

  return;
}

Score&
Score::operator=(const Score& rhs)
{
  Set_of_Values::operator=(rhs);

  return *this;
}

int
Score::operator!=(const Score& rhs) const
{
  return Set_of_Values::operator!=(rhs);
}

/*
  A typical score record will be
    plateid score1 score2 score3 ...
*/

int
Score::read_score(const const_IWSubstring& buffer)
{
  int i = 0;  // pointer into buffer
  const_IWSubstring token;
  (void)buffer.nextword(token, i);

  int j = 0;  // index into the _things array

  while (buffer.nextword(token, i)) {
    if (!token.numeric_value(_things[j]) || _things[j] < 0) {
      cerr << "Score::read_score: invalid token '" << token << "'\n";
      return 0;
    }

    j++;
    if (j >= score_histogram_size) {
      return 1;
    }
  }

  if (j < score_histogram_size) {
    cerr << "Score::read_score: too few values, j = " << j
         << " score hisogram = " << score_histogram_size << endl;
    return 0;
  }

  return 1;
}

int
Score::operator>(const Score& rhs) const
{
  for (int i = 0; i < score_histogram_size; i++) {
    if (_things[i] > rhs._things[i]) {
      return 1;
    }
  }

  return 0;
}

/*
  We need a class which can hold a histogram of diversity values
*/

typedef float similarity_type_t;

class Diversity : public Set_of_Values
{
 private:
 public:
  Diversity();

  int read_histogram(iwstring_data_source&);
};

Set_of_Values&
Set_of_Values::operator=(const Set_of_Values& rhs)
{
  for (int i = 0; i < _number_elements; i++) {
    _things[i] = rhs._things[i];
  }

  return *this;
}

Set_of_Values&
Set_of_Values::operator+=(const Set_of_Values& rhs)
{
  for (int i = 0; i < _number_elements; i++) {
    _things[i] += rhs._things[i];
  }

  return *this;
}

Set_of_Values&
Set_of_Values::operator-=(const Set_of_Values& rhs)
{
  for (int i = 0; i < _number_elements; i++) {
    _things[i] -= rhs._things[i];
  }

  return *this;
}

Diversity::Diversity()
{
  assert(diversity_histogram_size > 0);

  resize(diversity_histogram_size);

  for (int i = 0; i < _number_elements; i++) {
    _things[i] = 0;
  }

  return;
}

int
Set_of_Values::operator==(const Set_of_Values& rhs) const
{
  for (int i = 0; i < _number_elements; i++) {
    if (_things[i] != rhs._things[i]) {
      return 0;
    }
  }

  return 1;
}

int
Set_of_Values::operator!=(const Set_of_Values& rhs) const
{
  for (int i = 0; i < _number_elements; i++) {
    if (_things[i] != rhs._things[i]) {
      return 1;
    }
  }

  return 0;
}

std::ostream&
operator<<(std::ostream& os, const Set_of_Values& rhs)
{
  for (int i = 0; i < rhs.number_elements(); i++) {
    if (i) {
      os << ' ';
    }
    os << rhs[i];
  }

  return os;
}

class Plate
{
 private:
  IWString _id;
  //  score_t _score;

  Score _score;

  Diversity _intra_plate_diversity;

  //  During a spread run, each plate keeps track of its distances to
  //  the current configuration as it grows

  Diversity _diversity_to_current_config;

  /*
    We store the inter plate diversity with each plate. Note that each plate
    knows its plate number and it stores the diversity relative to each plate
    with a higher plate number
  */

  int _my_plate_number;

  Diversity* _inter_plate_diversity;

  int _selected;

  //  private functions

  Diversity&
  private_inter_plate_diversity(int) const;

 public:
  Plate();
  ~Plate();

  int plate_number() const {
    return _my_plate_number;
  }

  void set_plate_number(int p) {
    _my_plate_number = p;
  }

  const IWString& id() const {
    return _id;
  }

  void set_id(const IWString& i) {
    _id = i;
  }

  const Score& score() const {
    return _score;
  }

  int read_score(const const_IWSubstring&);

  int selected() const {
    return _selected;
  }

  void set_selected(int s) {
    _selected = s;
  }

  int allocate_inter_plate_diversity_array(int);

  int read_intra_plate_histogram(iwstring_data_source&);
  int read_inter_plate_diversity(iwstring_data_source&, int);

  const Diversity& intra_plate_diversity() const {
    return _intra_plate_diversity;
  }

  const Diversity& inter_plate_diversity(int) const;

  Diversity& diversity_to_current_config() {
    return _diversity_to_current_config;
  }
};

static Plate* plate = nullptr;

static int nplates = 0;

static int plates_to_choose = 0;

int
Diversity::read_histogram(iwstring_data_source& input)
{
  int rc = 0;  // how many components have we read

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    const_IWSubstring token = buffer.word(1);
    if ("samples" != token) {
      cerr << "Invalid histogram record '" << buffer << "'\n";
      return 0;
    }

    token = buffer.word(0);
    if (!token.numeric_value(_things[rc]) || _things[rc] < 0) {
      cerr << "Invalid count in histogram '" << buffer << "'\n";
      return 0;
    }

    rc++;
    if (rc == diversity_histogram_size) {
      return rc;
    }
  }

  cerr << "Diversity::read_histogram: premature eof\n";
  return 0;
}

Plate::Plate()
{
  _my_plate_number = -1;

  _inter_plate_diversity = nullptr;
}

Plate::~Plate()
{
  if (NULL != _inter_plate_diversity) {
    delete[] _inter_plate_diversity;
  }

  return;
}

int
Plate::read_intra_plate_histogram(iwstring_data_source& input)
{
  return _intra_plate_diversity.read_histogram(input);
}

/*
 */

int
Plate::allocate_inter_plate_diversity_array(int np)
{
  assert(np == nplates);
  assert(NULL == _inter_plate_diversity);

  if (_my_plate_number == np - 1) {
    return 1;
  }

  _inter_plate_diversity = new Diversity[np - _my_plate_number];

  if (NULL == _inter_plate_diversity) {
    return 0;
  }

  return 1;
}

Diversity&
Plate::private_inter_plate_diversity(int other_plate) const
{
  assert(other_plate > _my_plate_number);
  assert(NULL != _inter_plate_diversity);

  return _inter_plate_diversity[other_plate - _my_plate_number - 1];
}

const Diversity&
Plate::inter_plate_diversity(int other_plate) const
{
  return private_inter_plate_diversity(other_plate);
}

int
Plate::read_inter_plate_diversity(iwstring_data_source& input, int other_plate)
{
  Diversity& d = private_inter_plate_diversity(other_plate);

  return d.read_histogram(input);
}

int
Plate::read_score(const const_IWSubstring& buffer)
{
  return _score.read_score(buffer);
}

/*
  A configuration is characterised by a model score and a diversity
  score.
  The resizable array holds the indices of the plates which are
  in this configuration
*/

class Configuration : public resizable_array<int>
{
 private:
  Score _model_score;

  Diversity _diversity;

 public:
  Configuration();

  Score& model_score()
  {
    return _model_score;
  }

  const Score& model_score() const {
    return _model_score;
  }

  Diversity& diversity() {
    return _diversity;
  }

  const Diversity& diversity() const {
    return _diversity;
  }

  int write(std::ostream&) const;

  Configuration& operator=(const Configuration&);

  int operator<(const Configuration& rhs) const;
  int operator<=(const Configuration& rhs) const;
  int operator==(const Configuration& rhs) const;
  int operator>(const Configuration& rhs) const;
  int operator>=(const Configuration& rhs) const;
};

Configuration::Configuration()
{
  resize(nplates);

  return;
}

int
Configuration::write(std::ostream& os) const
{
  for (int i = 0; i < _model_score.number_elements(); i++) {
    os << ' ' << _model_score[i];
  }

  os << " :";

  for (int i = 0; i < _diversity.number_elements(); i++) {
    os << ' ' << _diversity[i];
  }

  os << " Selected";

  for (int i = 0; i < _number_elements; i++) {
    const Plate& p = plate[_things[i]];
    os << ' ' << p.id();
  }

  return os.good();
}

Configuration&
Configuration::operator=(const Configuration& rhs)
{
  resizable_array<int>::operator=(rhs);

  _model_score = rhs._model_score;

  _diversity = rhs._diversity;

  return *this;
}

int
Configuration::operator==(const Configuration& rhs) const
{
  if (!resizable_array<int>::operator==(rhs)) {
    return 0;
  }

  if (_model_score != rhs._model_score) {
    return 0;
  }

  if (_diversity != rhs._diversity) {
    return 0;
  }

  return 1;
}

/*
  During a spread run, we need some means of keeping track the desirability
  of plates to be added.
  Because searching through the configuration for the point at which to
  insert a plate number is expensive, the Added object keeps track of
  the insertion point
*/

class Added
{
 private:
  int _plate_number;

  Score _score;
  Diversity _diversity;

 public:
  Added();

  int print_considered(std::ostream&) const;

  int plate_number() const {
    return _plate_number;
  }

  void set_plate_number(int p) {
    _plate_number = p;
  }

  Diversity& diversity() {
    return _diversity;
  }

  const Diversity& diversity() const {
    return _diversity;
  }

  Score& score() {
    return _score;
  }

  const Score& score() const {
    return _score;
  }

  Added& operator=(const Added&);

  int operator>(const Added&) const;
  int model_bias_compare(const Added&) const;
};

Added::Added()
{
  _plate_number = -1;

  return;
}

Added&
Added::operator=(const Added& rhs)
{
  _plate_number = rhs._plate_number;

  _score = rhs._score;

  _diversity = rhs._diversity;

  return *this;
}

/*
  This is the guts of all optimisations. How do we compare two
  possible additions.
  We want to add the plate which adds the largest number of
  high scores (low score buckets) with the smallest number of
  close distances (low diversity buckets)
*/

static int score_multiplier[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};
static int diversity_multiplier[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};

int
Added::operator>(const Added& rhs) const
{
  // If either configuration has no data, then it cannot be greater

  if (_plate_number < 0) {
    return 0;
  }

  if (rhs._plate_number < 0) {
    return 1;
  }

  // We want high scores

  int s = 0;
  for (int i = 0; i < score_histogram_consider; i++) {
    s += score_multiplier[i] * (_score[i] - rhs._score[i]);
  }

  int d = 0;
  for (int i = 0; i < diversity_histogram_consider; i++) {
    d += diversity_multiplier[i] * (_diversity[i] - rhs._diversity[i]);
  }

  if (s > 0 && d < 0) {  // My score is better and I have fewer close distances
    return 1;
  }

  if (s < 0 && d > 0) {  // his score is better and he has fewer close distances
    return 0;
  }

  if (0 == s) {
    return (d < 0);
  }

  if (0 == d) {
    return (s > 0);
  }

  return model_bias_compare(rhs);
}

/*
  This is designed to be invoked from operator> when score says one
  thing and diversity say the opposite.
  Produce something in the range 0-1
*/

int
Added::model_bias_compare(const Added& rhs) const
{
  float fs = 0.0;
  for (int i = 0; i < score_histogram_consider; i++) {
    if (_score[i] == rhs._score[i]) {
      continue;
    }

    if (0 == rhs._score[i]) {
      fs += 1.0 / float(i + 1);
    } else if (0 == _score[i]) {
      fs -= 1.0 / float(i + 1);
    } else if (_score[i] > rhs._score[i]) {
      fs += float(rhs._score[i]) / float(_score[i]) / float(i + 1);
    } else {
      fs -= float(_score[i]) / float(rhs._score[i]) / float(i + 1);
    }
  }

  float ds = 0.0;
  for (int i = 0; i < diversity_histogram_consider; i++) {
    if (_diversity[i] == rhs._diversity[i]) {
      continue;
    }

    if (0 == rhs._diversity[i]) {
      ds -= 1.0 / float(i + 1);
    } else if (0 == _diversity[i]) {
      ds += 1.0 / float(i + 1);
    } else if (_diversity[i] > rhs._diversity[i]) {
      ds -= float(rhs._diversity[i]) / float(_diversity[i]) / float(i + 1);
    } else {
      ds += float(_diversity[i]) / float(rhs._diversity[i]) / float(i + 1);
    }
  }

#ifdef DEBUG_ADDED_COMPARISON
  for (int i = 0; i < score_histogram_consider; i++) {
    cerr << " (" << _score[i] << ' ' << rhs._score[i] << ')';
  }
  cerr << " :";
  for (int i = 0; i < diversity_histogram_consider; i++) {
    cerr << " (" << _diversity[i] << ' ' << rhs._diversity[i] << ')';
  }
  cerr << " Result: s = " << s << " d = " << d << endl;
  cerr << "fs = " << fs << " ds = " << ds << endl;
#endif

  if (fs > 0.0 && ds < 0.0) {
    return 1;
  }
  if (fs < 0.0 && ds > 0.0) {
    return 0;
  }

  // At this stage fs and ds have the same sign (both +ve or both -ve)

  float tmp = fs - ds + model_bias;

  // cerr << "fs = " << fs << " ds = " << ds << " tmp = " << tmp << endl;

  if (tmp > 0.0) {
    return 1;
  } else {
    return 0;
  }
}

int
Added::print_considered(std::ostream& os) const
{
  os << "plate " << _plate_number << " S:";

  for (int i = 0; i < score_histogram_consider; i++) {
    os << ' ' << _score[i];
  }
  os << " D:";
  for (int i = 0; i < diversity_histogram_consider; i++) {
    os << ' ' << _diversity[i];
  }

  return os.good();
}

/*
  We need a cross reference between plate id's and their index ino the plate array
*/

static IW_STL_Hash_Map_int Plate_ID_index_hash;

static int
read_inter_plate_diversities(iwstring_data_source& input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.nwords() < 2) {
      cerr << "Invalid inter plate diveristy record '" << buffer << "'\n";
      return 0;
    }

    IWString token = buffer.word(0);
    if (!Plate_ID_index_hash.contains(token)) {
      cerr << "Inter plate diversity: plate '" << token << "' not defined\n";
      return 0;
    }

    int i = Plate_ID_index_hash[token];
    assert(i >= 0 && i < nplates);

    token = buffer.word(1);
    if (!Plate_ID_index_hash.contains(token)) {
      cerr << "Inter plate diversity: plate '" << token << "' not defined\n";
      return 0;
    }

    int j = Plate_ID_index_hash[token];
    assert(j >= 0 && j < nplates);

    if (verbose > 1) {
      cerr << "Reading inter plate diversity between '" << plate[i].id() << "' and '"
           << plate[j].id() << "'\n";
    }

    if (i > j) {
      int tmp = j;
      j = i;
      i = tmp;
    }

    if (!plate[i].read_inter_plate_diversity(input, j)) {
      cerr << "Cannot read inter plate diversity '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_inter_plate_diversities(const const_IWSubstring& fname)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open inter plate diverisity file '" << fname << "'\n";
    return 0;
  }

  return read_inter_plate_diversities(input);
}

static int
allocate_inter_plate_diversity_arrays()
{
  assert(nplates > 0);

  for (int i = 0; i < nplates - 1; i++) {
    Plate& p = plate[i];

    assert(i == p.plate_number());

    if (!p.allocate_inter_plate_diversity_array(nplates)) {
      return 0;
    }
  }

  return 1;
}

static int
read_intra_plate_diversities(iwstring_data_source& input)
{
  int rc = 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    IWString plate_id = buffer.word(0);
    if (Plate_ID_index_hash.end() == Plate_ID_index_hash.find(plate_id)) {
      cerr << "Yipes, plate '" << plate_id << "' not defined\n";
      return 0;
    }

    Plate& p = plate[Plate_ID_index_hash[plate_id]];

    if (!p.read_intra_plate_histogram(input)) {
      cerr << "Cannot read intra plate diversity score for plate '" << plate_id << "'\n";
      return 0;
    }

    if (verbose > 1) {
      cerr << "Read intra-plate diversity for plate '" << p.id() << "'\n";
    }

    rc++;
  }

  if (rc != nplates) {
    cerr << "huh, read intra plate diversity scores for " << rc
         << " plates, but there are " << nplates << " plates\n";
    return 0;
  }

  return rc;
}

static int
read_intra_plate_diversities(const const_IWSubstring& fname)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open intra plate diversity score file '" << fname << "'\n";
    return 0;
  }

  return read_intra_plate_diversities(input);
}

static int
read_scores(iwstring_data_source& input)
{
  input.set_skip_blank_lines();
  input.set_strip_trailing_blanks();
  input.set_ignore_pattern("^#");

  nplates = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.nwords() < 2) {
      cerr << "Invalid score record '" << buffer << "'\n";
      return 0;
    }

    IWString p = buffer.word(0);

    if (Plate_ID_index_hash.end() != Plate_ID_index_hash.find(p)) {
      cerr << "Duplicate score for plate '" << p << "'\n";
      return 0;
    }

    Plate_ID_index_hash[p] = nplates;

    if (!plate[nplates].read_score(buffer)) {
      cerr << "Cannot read score values from '" << buffer << "'\n";
      return 0;
    }

    plate[nplates].set_plate_number(nplates);
    plate[nplates].set_id(p);

    nplates++;
  }

  if (verbose) {
    cerr << "Read scores for " << nplates << " plates\n";
  }

  return nplates;
}

static int
read_scores(const const_IWSubstring& fname)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Yipes, cannot open scores file '" << fname << "'\n";
    return 0;
  }

  assert(NULL == plate);

  int nr = input.records_remaining();

  plate = new Plate[nr];

  if (verbose) {
    cerr << "Scores file contains " << nr << " records\n";
  }

  return read_scores(input);
}

static void
do_choose_first_plate_best(int& first_plate)
{
  Score best = plate[0].score();
  first_plate = 0;

  for (int i = 1; i < nplates; i++) {
    Plate& pi = plate[i];

    if (pi.score() > best) {
      first_plate = i;
      best = pi.score();
    }
  }

  return;
}

static void
unselect_all()
{
  for (int i = 0; i < nplates; i++) {
    plate[i].set_selected(0);
  }

  return;
}

/*
  We are thinking of selecting plate ZPLATE and need to know
  what it will contribute to the configuration
*/

/*static void
compute_added_score_and_diversity (const Configuration & config,
                                    Added & added)
{
  int zplate = added.plate_number ();

  Plate & pi = plate[zplate];

  added.score () = pi.score ();
  added.diversity () == pi.intra_plate_diversity ();

  int csize = config.number_elements ();

  int insert_point = -1;

  for (int i = 0; i < csize; i++)
  {
    int j = config[i];

    assert (zplate != j);     // cannot have the same plate included twice

    if (zplate < j)
    {
      insert_point = i;
      added.diversity () += pi.inter_plate_diversity (j);
    }
    else
    {
      added.diversity () += plate[j].inter_plate_diversity (zplate);
    }
  }

  added.set_insert_point (insert_point);

  return;

}*/

/*
  A plate has been selected. We need to adjust the diversity_to_current_config
  value for every unselected plate.
*/

static void
select_plate(Configuration& config, int zplate)
{
  Plate& p = plate[zplate];

  assert(!p.selected());

  p.set_selected(1);

  config.model_score() += p.score();
  config.diversity() += p.intra_plate_diversity();

  int csize = config.number_elements();

  if (0 == csize)  // must be just starting
  {
    config.add(zplate);
    return;
  }

  int insert_point = -1;
  for (int i = 0; i < csize; i++) {
    int j = config[i];

    Plate& pj = plate[j];

    Diversity& d = pj.diversity_to_current_config();

    if (j < zplate) {
      insert_point = i;
      d += pj.inter_plate_diversity(zplate);
      config.diversity() += pj.inter_plate_diversity(zplate);
    } else {
      d += p.inter_plate_diversity(j);
      config.diversity() += p.inter_plate_diversity(j);
    }
  }

  if (insert_point < 0) {  // new plate greater than all existing plate numbers
    config.add(zplate);
  } else {
    config.insert_before(insert_point, zplate);
  }

  return;
}

static int
spread(std::ostream& output)
{
  unselect_all();

  Configuration config;

  int first_plate;
  if (choose_first_plate_randomly) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> u(0, nplates - 1);
    first_plate = u(rng);
  } else if (choose_first_plate_best) {
    do_choose_first_plate_best(first_plate);
  } else {
    first_plate = 0;
  }

  select_plate(config, first_plate);

  output << "Selected 1 " << plate[first_plate].id() << " plate " << first_plate << endl;

  assert(plates_to_choose > 0);

  while (config.number_elements() < plates_to_choose && output.good()) {
    Added best;

    for (int i = 0; i < nplates; i++) {
      Plate& p = plate[i];

      if (p.selected()) {
        continue;
      }

      Added d;
      d.set_plate_number(i);

      d.score() += p.score();
      d.diversity() += p.intra_plate_diversity();
      d.diversity() += p.diversity_to_current_config();

      //    cerr << "Checking plate " << i << endl;
      if (d > best) {
        best = d;
        //      cerr << "Plate " << i << " was better\n";
        //      best.print_considered (cerr);
        //      cerr << endl;
      }
    }

    if (best.plate_number() < 0) {
      cerr << "Huh, no better plate found!!!\n";
      return 0;
    }

    int p = best.plate_number();

    select_plate(config, p);

    output << "Selected " << config.number_elements() << ' ' << plate[p].id() << ' ';
    best.print_considered(output);
    output << endl;
  }

  output << endl;

  config.write(output);
  output << endl;

  return output.good();
}

static void
compute_score(Configuration& config)
{
  int csize = config.number_elements();

  for (int i = 0; i < csize; i++) {
    const Plate& pi = plate[config[i]];

    config.model_score() += pi.score();

    config.diversity() += pi.intra_plate_diversity();

    for (int j = i + 1; j < csize; j++) {
      config.diversity() += pi.inter_plate_diversity(config[j]);
    }
  }

  return;
}

static void
select_initial_plates(Configuration& config)
{
  unselect_all();

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> u(0, nplates - 1);

  while (config.number_elements() < plates_to_choose) {
    const int i = u(rng);

    Plate& p = plate[i];
    if (p.selected()) {
      continue;
    }

    p.set_selected(1);
    config.add(i);
  }

  compute_score(config);

  return;
}

/*static int
single_optimisation_step (std::ostream & output)
{
  return output.good ();
}*/

/*static int
optimise_layout (Configuration & config,
                 std::ostream & output)
{
  Score score;
  Diversity d;
  compute_score (score, d);

  int nopt = 10;
  for (int i = 0; i < nopt; i++)
  {
    single_optimisation_step (output);
  }

  return 1;
}*/

static int
plate_layout(int nsteps, std::ostream& output)
{
  for (int i = 0; i < nsteps; i++) {
    Configuration config;
    select_initial_plates(config);
    config.write(output);
  }

  if (verbose) {
  }

  return output.good();
}

static int
int_comparitor_larger(const int* p1, const int* p2)
{
  if (*p1 < *p2) {
    return -1;
  } else if (*p1 > *p2) {
    return 1;
  }

  return 0;
}

static int
evaluate_existing_plate(iwstring_data_source& input, std::ostream& output)
{
  Configuration config;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    const_IWSubstring p = buffer.word(0);

    if (Plate_ID_index_hash.end() == Plate_ID_index_hash.find(p)) {
      cerr << "Plate '" << p << "' has not been defined\n";
      return 0;
    }

    int i = Plate_ID_index_hash[p];

    config.add(i);

    Plate& pl = plate[i];
    pl.set_selected(1);
  }

  config.sort(int_comparitor_larger);

  compute_score(config);

  return config.write(output);
}

static int
evaluate_existing_plate(const const_IWSubstring& fname, std::ostream& output)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open existing plate layout '" << fname << "'\n";
    return 0;
  }

  return evaluate_existing_plate(input, output);
}

static int
plate_layout(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vR:E:S:h:H:n:p:V:P:b:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('S')) {
    cerr << "Must specify plate scores via the -S option\n";
    usage(2);
  }

  if (!cl.option_present('E')) {
    cerr << "Must specify inter plate diversities via the -E option\n";
    usage(2);
  }

  if (!cl.option_present('R')) {
    cerr << "Must specify intra plate diversities via the -R option\n";
    usage(4);
  }

  if (!cl.option_present('h')) {
    cerr << "You must specify the diversity histogram size via the -h option\n";
    usage(5);
  }

  if (!cl.option_present('H')) {
    cerr << "You must specify the score histogram size via the -H option\n";
    usage(5);
  }

  // If we are checking a pre-computed configuration, we don't need a number
  // of plates to select

  if (cl.option_present('V')) {
    ;
  } else if (cl.option_present('P')) {  // running spread
    ;
  } else if (!cl.option_present('p')) {
    cerr << "Must specify the number of plates to choose via the -p option\n";
    usage(11);
  }

  if (cl.option_present('h')) {
    if (!cl.value('h', diversity_histogram_size) || diversity_histogram_size < 1) {
      cerr << "The diversity histogram size option (-h) must be followed by a positive "
              "whole number\n";
      usage(7);
    }
    if (verbose) {
      cerr << "Diversity histogram sizes to " << diversity_histogram_size << endl;
    }
  }

  if (cl.option_present('H')) {
    if (!cl.value('H', score_histogram_size) || score_histogram_size < 1) {
      cerr << "The score histogram size option (-H) must be followed by a positive whole "
              "number\n";
      usage(7);
    }
    if (verbose) {
      cerr << "Score histogram sized to " << score_histogram_size << endl;
    }
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');
    if (!read_scores(s)) {
      cerr << "Cannot read score file\n";
      return 31;
    }

    if (0 == nplates) {
      cerr << "Huh, no plates after reading scores\n";
      return 11;
    }
  }

  if (cl.option_present('R')) {
    const_IWSubstring r = cl.string_value('R');
    if (!read_intra_plate_diversities(r)) {
      cerr << "Cannot read intra plate diversity scores from '" << r << "'\n";
      return 21;
    }
  }

  if (cl.option_present('E')) {
    if (!allocate_inter_plate_diversity_arrays()) {
      cerr << "Cannot allocate inter plate diversity arrays\n";
      return 0;
    }

    const_IWSubstring e = cl.string_value('E');
    if (!read_inter_plate_diversities(e)) {
      cerr << "Cannot read inter plate diversity scores from '" << e << "'\n";
      return 21;
    }
  }

  if (cl.option_present('V')) {
    int i = 0;
    const_IWSubstring v;
    while (cl.value('V', v, i++)) {
      if (!evaluate_existing_plate(v, std::cout)) {
        cerr << "Cannot process existing plate layout in file '" << v << "'\n";
        return 12;
      }
    }

    return 0;
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', model_bias)) {
      cerr << "The model bias option (-b) must be followed by a number\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Model bias " << model_bias << endl;
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', plates_to_choose) || plates_to_choose < 1 ||
        plates_to_choose >= nplates) {
      cerr << "The plates to choose option (-p) must be a value between 1 and "
           << (nplates - 1) << endl;
      usage(17);
    }

    if (verbose) {
      cerr << "Will select " << plates_to_choose << " of " << nplates << " plates\n";
    }
  }

  if (cl.option_present('P'))  // doing a spread computation
  {
    int i = 0;
    const_IWSubstring p;
    while (cl.value('P', p, i++)) {
      if ("random" == p) {
        choose_first_plate_randomly = 1;
        if (verbose) {
          cerr << "First spread plate chosen randomly\n";
        }
      } else if ("best" == p) {
        choose_first_plate_best = 1;
        if (verbose) {
          cerr << "The highest scoring plate will be selected first\n";
        }
      } else if ("" == p) {
        ;
      } else {
        cerr << "Unrecognised spread qualifier '" << p << "'\n";
        usage(4);
      }
    }

    if (0 == plates_to_choose) {
      plates_to_choose = nplates;
    }

    spread(std::cout);

    return 0;
  }

  int nsteps;
  if (cl.option_present('n')) {
    if (!cl.value('n', nsteps) || nsteps < 1) {
      cerr << "The number of steps option (-n) must be followed by a positive number\n";
      usage(16);
    }

    if (verbose) {
      cerr << "Will perform " << nsteps << " of optimisation\n";
    }
  } else {
    nsteps = 1000;
    if (verbose) {
      cerr << "Will perform 1000 steps by default\n";
    }
  }

  if (!plate_layout(nsteps, std::cout)) {
    return 18;
  }

  return 0;
}

int
main(int argc, char** argv)
{
  int rc = plate_layout(argc, argv);

  return rc;
}
