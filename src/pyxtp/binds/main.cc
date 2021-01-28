
#include "pyxtp.hpp"

int main() {
  pyxtp::PyXTP pyxtp;
  std::cout << "starting:\n";
  pyxtp.Initialize("eanalyze", 1);

  pyxtp::call_calculator("eanalyze", 1);
}
