#include <chrono>
#include <iostream>
#include <vector>

#include "../include/vector.hpp"

#define UNIT_SECONDS 1000.  // in microseconds
#define LABEL_SECONDS "ms"

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

using namespace std;

int main() {
  vector<double> x, y;
  double a = 0.5;
  uint length = 100000;
  uint rep = 1000, count = 0;
  // random double entries
  for (uint i = 0; i < length; i++)
    x.push_back(rand() % 100);

  chrono::steady_clock::time_point begin, end;
  uint64_t total = 0;

  START_TIMER();
  while (count < rep) {
    y = x % a;
    count++;
  }
  END_TIMER();
  ADD_TIME_TO(total);
  cout << "Time for " << rep << " multiplications (with for): " << total / UNIT_SECONDS << LABEL_SECONDS << endl;

  count = 0;
  total = 0;
  START_TIMER();
  while (count < rep) {
    y = x * a;
    count++;
  }
  END_TIMER();
  ADD_TIME_TO(total);
  cout << "Time for " << rep << " multiplications (naive): " << total / UNIT_SECONDS << LABEL_SECONDS << endl;
}
