
#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <votca/tools/elements.h>

using namespace std;
using namespace votca::tools;

int main(void) {

  cout << "Testing: constructor" << endl;
  { Elements ele; }

  cout << "Testing: getVdWChelpG" << endl;
  {
    Elements ele;
    assert(ele.getVdWChelpG("H") == 1.45);
    bool test = false;
    try {
      ele.getVdWChelpG("Blah");
    } catch (...) {
      test = true;
    }
    assert(test);
  }

  cout << "Testing: getNucCrgECP" << endl;
  {
    Elements ele;
    bool     test = false;
    try {
      ele.getNucCrgECP("He");
    } catch (...) {
      test = true;
    }
    assert(test);
  }

  cout << "Testing: getMass" << endl;
  {
    Elements ele;
    assert(ele.getMass("K") = 39.098);
  }

  cout << "Testing: getEleNum" << endl;
  {
    Elements ele;
    assert(ele.getEleNum("Li") == 3);
  }

  cout << "Testing: getEleName" << endl;
  {
    Elements ele;
    assert(ele.getEleName(17) == "Cl");
  }

  cout << "Testing: getEleShort" << endl;
  {
    Elements ele;
    assert(ele.getEleShort("MAGNESIUM") == "Mg");
  }

  cout << "Testing: getEleFull" << endl;
  {
    Elements ele;
    assert(ele.getEleFull("Ge") == "GERMANIUM");
  }

  cout << "Testing: getVdWMK" << endl;
  {
    Elements ele;
    assert(ele.getVdWMK("F") == 1.35);
    bool test = false;
    try {
      ele.getVdWMK("Pb");
    } catch (...) {
      test = true;
    }
    assert(test);
  }
  return 0;
}
