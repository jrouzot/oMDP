#include "../header/Algorithm.hpp"
#include <limits>
#include <sys/resource.h>

// template <typename SequenceOfRankings>
dataflow::Solution::Solution(const vector<vector<vector<int>>> &prio) {
  copy(prio);
}

void dataflow::Solution::copy(const vector<vector<vector<int>>> &prio) {
  priority.resize(prio.size());
  for (auto i{0}; i < prio.size(); ++i) {
    priority[i].resize(prio[i].size());
    for (auto r{0}; r < prio[i].size(); ++r) {
      priority[i][r] = prio[i][r];
    }
  }
}

dataflow::Solution::Solution() {}

dataflow::Solution::~Solution() {}

double dataflow::cpu_time(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static unsigned long x = 123456789, y = 362436069, z = 521288629;

void dataflow::seed(const unsigned long s) {
  x = s;
  y = 362436069;
  z = 521288629;
}

unsigned long dataflow::random(void) { // period 2^96-1
  unsigned long t;
  x ^= x << 16;
  x ^= x >> 5;
  x ^= x << 1;

  t = x;
  x = y;
  y = z;
  z = t ^ x ^ y;

  return z;
}

double dataflow::random(const double lb, const double ub) {
  // unsigned long r{};
  return static_cast<double>(random()) /
             static_cast<double>(numeric_limits<unsigned long>::max()) *
             (ub - lb) +
         lb;
}
