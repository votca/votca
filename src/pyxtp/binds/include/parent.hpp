
#if !defined(PARENT_H_)
#define PARENT_H_

namespace pyxtp {

class Parent {
 public:
  Parent() = default;
  ~Parent() = default;

  void set_value(int n) { _value = n; }

 protected:
  int _value = 0;
};

}  // namespace pyxtp

#endif  // PARENT_H_
