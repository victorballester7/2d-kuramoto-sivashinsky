#include <time.h>

#include <iostream>
#include <vector>
using namespace std;

double time_elapsed(timespec& start, timespec& end) {
  return ((1e9 * end.tv_sec + end.tv_nsec) -
          (1e9 * start.tv_sec + start.tv_nsec)) /
         1.0e9;
}

int main() {
  vector<int> v(1e9, 42);

  timespec start, end;

  //  Old-fashioned loop.
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  size_t size = v.size();
  for (size_t i = 0; i < size; i++) {
    v[i] *= v[i];
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  cout << "Old-fashioned loop: " << time_elapsed(start, end) << " seconds\n";

  //  Range-based loop.
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  for (int& val : v) {
    val *= val;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  cout << "Range-based loop: " << time_elapsed(start, end) << " seconds.\n";
}