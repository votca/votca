
#if !defined(EXAMPLE_H_)
#define EXAMPLE_H_

#include "parent.hpp"

namespace pyxtp {

class Example : public Parent {
 public:
  Example() = default;
  ~Example() = default;
};
}  // namespace pyxtp
#endif  // EXAMPLE_H_
