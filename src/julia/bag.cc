// Box interface

#include <iostream>
#include <vector>

#include "jlcxx/jlcxx.hpp"

namespace bag {

struct Bag {
  std::vector<int> contents;

  Bag(size_t mysize) : contents(mysize) {
  }

  int& operator[](int ndx) {
    if (ndx < 0 || ndx >= contents.size()) { std::cerr << "Out of range " << ndx << '\n';}
    return contents[ndx];
  }
  int& item(int ndx) {
    return contents[ndx];
  }
};

Bag
BagMaker(size_t s) {
  return Bag(s);
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.add_type<Bag>("Bag")
    .constructor<size_t>()
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](Bag& bag, size_t ndx) {
      std::cerr << "getindex ndx " << ndx << '\n';
      return bag[ndx];
    }
  );
  mod.method("length",
    [](const Bag& bag) {
      std::cerr << "length " << bag.contents.size() << '\n';
      return bag.contents.size();
    }
  );
  mod.method("setindex!",
    [](Bag& bag, int64_t value, int64_t ndx)->void {
      std::cerr << "setindex ndx " << ndx << " value " << value << '\n';
      bag[ndx] = value;
    }
  );
#ifdef SETINDECX_
  mod.method("setindex!",
    [](jlcxx::BoxedValue<Bag>& boxed_bag, int64_t value, int64_t ndx)->void{
      std::cerr << "setindex Boxed<box> ndx " << ndx << " value " << value << '\n';
      Bag& bag = jlcxx::unbox<Bag&>(boxed_bag);
      bag[ndx] = value;
    }
  );
#endif
  mod.unset_override_module();
  mod.method("new_bag",
    [](size_t s) {
      return BagMaker(s);
    }
  );
}

}  // namespace bag
