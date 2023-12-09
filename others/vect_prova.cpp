#include <chrono>
#include <iostream>
#include <vector>

#include "../include/matplotlibcpp.h"
#include "../include/vector.hpp"

using namespace std;
namespace plt = matplotlibcpp;

int main() {
  vector<double> y = {1, 3, 2, 4, 3, 5, 4, 6, 5, 7};
  plt::plot(y);
  // plt::show();
  plt::save("minimal.pdf");
  plt::detail::_interpreter::kill();  // To avoid a final segmentation fault
}
