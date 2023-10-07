#include "iwstring.h"

// Small test thing to see if
// const const_IWSubstring& = IWString
// creates a temporary. Indeed, it seems to do just that.

namespace test_pass_by_ref {

class Foo {
 public:
  Foo() {};
  void set_foo(const char * s) {
    foo_ = s;
  }
  const IWString& foo() const { return foo_;}
 private:
  IWString foo_;
};

void RunTestPassByRef() {
  Foo foo;
  foo.set_foo("bar");
  const const_IWSubstring & ref = foo.foo();
  if (ref == "bar")
    exit(0);
  exit(1);
}

}

int main() {
  test_pass_by_ref::RunTestPassByRef();
}
