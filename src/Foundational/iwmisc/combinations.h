#ifndef FOUNDATIONAL_IWMISC_COMBINATIONS_H
#define FOUNDATIONAL_IWMISC_COMBINATIONS_H

#include <iostream>
#include <vector>

namespace combinations {

// Class to facilitate exploring combinations of things.
template <typename T>
class Combinations {
  private:
    // The number if items in each slot.
    std::vector<int> _count;

    bool _first_call;

  public:
    Combinations(const std::vector<T>& count);

    int Next(std::vector<T>& state);

    void Reset();
};

template <typename T>
Combinations<T>::Combinations(const std::vector<T>& count) : _count(count) {
  _first_call = true;
}

template <typename T>
int
Combinations<T>::Next(std::vector<T>& state) {
  if (_first_call) {
    std::fill_n(state.begin(), _count.size(), 0);
    _first_call = false;
    return 1;
  }

  //std::cerr << " sizeof _count " << _count.size() << " state " << state.size() << '\n';
  for (int i = state.size() - 1; i >= 0; --i) {
    //std::cerr << " i = " << i << '\n';
    if (state[i] < _count[i] - 1) {
      ++state[i];
      return 1;
    }
    state[i] = 0;
  }

  return 0;
}

template <typename T>
void
Combinations<T>::Reset() {
  _first_call = true;
}

};  // namespace combinations

#endif  // FOUNDATIONAL_IWMISC_COMBINATIONS_H
