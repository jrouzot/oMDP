
#ifndef __Algorithm_HPP
#define __Algorithm_HPP

#include <iostream>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>

#include "MonotoneVector.hpp"
#include "Options.hpp"
#include "Simulator.hpp"
#include "Event.hpp"
#include "State.hpp"

using namespace std;

namespace dataflow {

/**********************************************
* Algorithm
**********************************************/

double cpu_time(void);

void seed(const unsigned long s);

unsigned long random(void);

double random(const double lb, const double ub);

class Solution {

public:
  Solution();
  ~Solution();

  // template <typename SequenceOfRankings>
  Solution(const vector<vector<vector<int>>> &prio);

  void copy(const vector<vector<vector<int>>> &prio);

  vector<vector<vector<int>>> priority;
};

// // template <typename SequenceOfRankings>
// Solution::Solution(const vector<vector<vector<int>>> priority &prio) {
//   priority.resize(prio.size());
//   for (auto i{0}; i < prio.size(); ++i) {
//     priority[i].resize(prio[i].size());
//     for (auto r{0}; r < prio[i].size(); ++r) {
//       priority[i][r] = prio[i][r];
//     }
//   }
// }
//
// Solution::Solution() {}
//
// Solution::~Solution() {}

template <typename TimeT, typename DataT> class Algorithm {
public:
  Algorithm(const Instance<TimeT, DataT> &i, Options &opt);

  void setTargetMargin(const double margin);

  Solution *readSolution(const string &filename);
  double playSolution(vector<vector<vector<int>>> &prio,
                      const bool verbose = false);
  double playDownlink(const int i, const int j, vector<vector<int>> &prio,
                      const bool verbose = false);
  double playDownlinkWithoutTransfer(const int i, const int j,
                                     const bool verbose = false);
  double playDownlinkWithParallelTransfer(const int i, const int j,
                                          const bool verbose = false);

  void initialise();

  double multipleWindows(const double m);
  double randomizedMultipleWindows(const double m, const int nruns);
  // double multipleWindowsPlus(const double m);

  double downlinkCountHeuristic(const double target_margin = 0,
                                const int flag = 0);

  Solution *lnsDescent();
  Solution *iterativeLeveling();
  Solution *randomWalk();
  Solution *maxMultipleWindows();

  double relaxBandwidth();
  double greedyOptimize();
  bool singleWindow(const int i);
  double singleWindowMax(const int i);

  ostream &display(ostream &os) const;

  // int option.verbosity{0};
  // int option.runs{100};
  // double option.factor{1};

  double real_margin(const double margin) const;

private:
  const Instance<TimeT, DataT> &data;

  Options &option;

  // simulator for the whole horizon
  Simulator<TimeT, DataT> simulator;

  // target handover usage
  vector<vector<DataT>> handover;

  vector<vector<double>> hdelta;
  vector<double> marginat;
  vector<vector<DataT>> maxusage;
  // vector<vector<DataT>> fill;
  // the usage at the end of the previous downlink (or at time 0 for the
  // first
  // downlink)
  vector<vector<DataT>> usage;
  // the pointer in data to the end of the previous downlink (or at time 0
  // for
  // the first downlink)
  vector<vector<int>> pointer;

  // all events, ordered  backward
  vector<Event<TimeT>> events;

  // all events, ordered  backward
  vector<DataT> target;
  // in order to ignore some capacities
  vector<DataT> nolimit;

  // pointers to first, last and all end-of-downlink event
  vector<typename vector<Event<TimeT>>::const_iterator> window;

  //
  vector<vector<vector<int>>> priority;

  vector<vector<int>> best_priority;

  bool debug_flag{false};

  SparseSet<> thebuffers;
  // buffers.reserve(data.numBuffers());

  bool fail_flag{false};

  void postHandoverConstraints(const double margin);
  double getHandoverDeltas();
  double repair(const double base, Solution *repaired, const double epsilon = .00001);

  // double firstLossHeuristic(const double target_margin = 0);
  double randomizedDownlinkCountHeuristic(const double best_margin,
                                          const double target_margin,
                                          int &nruns);

  double computeMinimumMargin(const int i) const;
  void getMaxUsage(const int i);

  void clearHandover();

  void computeDCPriority(const int i);
  void computeFLPriority(const int i);

  DataT capacity(const int b) const;

  double marginOf(const int i, const int b) const;
  double percentFull(const int i, const int b) const;
  DataT delta(const int i, const int b) const;

  TimeT startOf(const int i) const { return (i ? (window[i] - 1)->time : 0); }
};

template <typename TimeT, typename DataT>
DataT Algorithm<TimeT, DataT>::capacity(const int b) const {
  return option.factor * data.capacity(b);
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::real_margin(const double margin) const {
  return option.factor * margin - option.factor + 1;
}

template <typename TimeT, typename DataT>
Algorithm<TimeT, DataT>::Algorithm(const Instance<TimeT, DataT> &i,
                                   Options &opt)
    : data(i), option(opt), simulator(data) {}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::setTargetMargin(const double margin) {
  for (auto b{0}; b < data.numBuffers(); ++b) {
    assert(nolimit[b] > capacity(b));
    target[b] = capacity(b) * (1.0 - margin);
  }
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::marginOf(const int i, const int b) const {
  return static_cast<double>(capacity(b) - maxusage[i][b]) /
         static_cast<double>(capacity(b));
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::percentFull(const int i, const int b) const {
  return static_cast<double>(maxusage[i][b]) / static_cast<double>(capacity(b));
}

template <typename TimeT, typename DataT>
DataT Algorithm<TimeT, DataT>::delta(const int i, const int b) const {
  return maxusage[i][b] - usage[i][b];
}

template <typename TimeT, typename DataT>
Solution *Algorithm<TimeT, DataT>::readSolution(const string &filename) {
  ifstream ifs(filename);

  Solution *sol = new Solution();

  int i;
  // int nb{data.numBuffers()};

  sol->priority.resize(data.numDownlinks() + 1);

  for (auto w{0}; w < data.numDownlinks(); ++w) {

    ifs >> i;
    assert(w == i);

    sol->priority[i].clear();

    auto nb{0};

    while (nb < data.numBuffers()) {
      int n;
      ifs >> n;

      nb += n;

      sol->priority[i].resize(sol->priority[i].size() + 1);
      for (auto k{0}; k < n; ++k) {
        int b;
        ifs >> b;

        assert(b >= 0 and b < data.numBuffers());

        sol->priority[i].back().push_back(b);
      }
    }
  }

  return sol;
}

template <typename TimeT, typename DataT>
Solution* Algorithm<TimeT, DataT>::lnsDescent() {

  // compute, for every buffer i and downlink j, the minimum delta in usage for buffer i in window j
  getHandoverDeltas();

  // compute a trivial optimistic bound
  double ubmargin{relaxBandwidth()};
  // if (option.verbosity >= 0)
  //   cout << "ub " << real_margin(ubmargin) << " " << cpu_time()*1000 << endl;

  // compute a pessimistic bound using downlink count
  double lbmargin{downlinkCountHeuristic(0, 0)};
  // if (option.verbosity >= 0)
  //   cout << "sol " << real_margin(lbmargin) << " " << cpu_time()*1000 << endl;

  cout << "obj=" << (1 - real_margin(lbmargin)) * 100 << endl;
  cout << "time=" << cpu_time()*1000 << endl;
  // cout << "sol " << real_margin(lbmargin) << " " << cpu_time()*1000 << endl;

  Solution *best_solution = new Solution(priority);

  // stop when the gap is small enough or when the iteration limit is reached
  // for (auto i{0}; lbmargin + .000001 < ubmargin and i < option.runs; ++i) {
  
  // stop when the iteration limit is reached or time limit is reached
  // for (auto i{0}; i < option.runs;++i) {
  while(cpu_time()*1000 < option.htime_limit) {

    // if (option.verbosity > 0)
    //   cout << "run " << i << ": " << real_margin(lbmargin) << ".."
    //        << real_margin(ubmargin) << " " << cpu_time()*1000 << endl;

		int choice{option.heuristic};
		if(choice > 3 or choice < 0)
			choice = (random() % 4);


    // compute a new solution using downlink count (default for heuristic is 1, so w.r.t. the current pessimistic bound)
    double current;
    if (choice == 0)
      current = downlinkCountHeuristic(0, -1);
    else if (choice == 1)
      current = downlinkCountHeuristic(lbmargin, -1);
    else if (choice == 2)
      current = downlinkCountHeuristic(ubmargin, -1);
    else
      current = downlinkCountHeuristic((lbmargin + ubmargin) / 2.0, -1);

    if (current > lbmargin) {
      lbmargin = current;
      best_solution->copy(priority);
      cout << "obj=" << (1 - real_margin(lbmargin)) * 100 << endl;
      cout << "time=" << cpu_time()*1000 << endl;
      // if (option.verbosity >= 0)
      // cout << "sol " << real_margin(lbmargin) << " " << cpu_time()*1000 << endl;
    }

    // repair the solution so that the margin is strictly greater than the current pessimistic bound
    lbmargin = repair(lbmargin, best_solution, option.epsilon);

  }

  // playSolution(priority, true);

  auto real_margin = playSolution(best_solution->priority, false);


  cout << "c real margin = " << real_margin << endl;

  return best_solution;
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::repair(const double base, Solution *best_solution, const double epsilon) {

  double current{base};

  auto fail{false};
  while (not fail && cpu_time()*1000 < option.htime_limit) {

    double minimum_margin{1};

    if (option.verbosity > 1)
      cout << "restart: " << minimum_margin << endl;

    double target_margin = current + epsilon;
    setTargetMargin(target_margin);
    postHandoverConstraints(target_margin);

    if (option.verbosity > 1)
      cout << "set new target: " << target_margin << endl;

    int i{0};
    while (i <= data.numDownlinks()) {

      auto margin{playDownlink(i, i + 1, priority[i])};

      auto success{simulator.excess_buffers.empty() and
                   margin >= target_margin};

      if (not success and option.verbosity > 2)
        cout << "window " << i << " failed (" << margin << "/" << target_margin
             << "); " << simulator.excess_buffers.count() << "\n";

      if (not success and singleWindow(i)) {
        priority[i].clear();
        priority[i] = best_priority;
        margin = playDownlink(i, i + 1, priority[i]);

        if (option.verbosity > 2)
          cout << "window " << i << " repaired (" << margin << ")\n";

        success = true;
      }

      if (success) {
        minimum_margin = min(margin, minimum_margin);
        if (option.verbosity > 2)
          cout << "window " << i << " ok (" << margin << "/" << minimum_margin
               << ")\n";
        ++i;
      } else {
        if (option.verbosity > 2) {
          cout << "window " << i << " infeasible at " << target_margin << endl;
          cout << setw(31) << " ";
          for (auto b{0}; b < data.numBuffers(); ++b) {
            cout << setw(6) << setprecision(3)
                 << (handover[i][b] / capacity(b));
          }
          cout << endl;
        }

        if (i == 0) {
          fail = true;
          break;
        }

        priority[i].clear();
        priority[i] = best_priority;
        margin = playDownlink(i, i + 1, priority[i]);

        int b = simulator
                    .excess_buffers[random() % simulator.excess_buffers.size()];

        auto delta(min(target[b] - simulator.getMaxUsage(b),
                       handover[i][b] - simulator.usage(b)));

        if (option.verbosity > 2) {
          cout << "data loss on buffers " << simulator.excess_buffers
               << " put a constraint on " << b << " @" << (i - 1) << " ("
               << handover[i - 1][b] << "->" << (usage[i][b] + delta)
               << ")\n and backtrack to " << (i - 1) << endl;
        }

        handover[i - 1][b] = (usage[i][b] + delta);

        --i;
      }
    }

    if (not fail) {
      current = minimum_margin;
      // delete best_solution;
      best_solution->copy(priority);
      // if (option.verbosity >= 0)
      //   cout << "sol " << real_margin(current) << " " << cpu_time()*1000 << endl;
      cout << "obj=" << (1 - playSolution(best_solution->priority, false)) * 100 << endl;
      cout << "time=" << cpu_time()*1000 << endl;
      // cout << "sol " << real_margin(current) << " " << cpu_time()*1000 << endl;
    }
  }


  return current;
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::playDownlink(const int i, const int j,
                                             vector<vector<int>> &prio,
                                             const bool verbose) {

  simulator.initialise(pointer[i], usage[i]);

  assert(simulator.excess_buffers.empty());



  if (verbose)
    simulator.print_flag = true;
  simulator.run(window[i], window[j], startOf(i), prio, target, handover[i]);
  if (verbose)
    simulator.print_flag = false;

  for (auto b{0}; b <= data.numBuffers(); ++b) {
    pointer[j][b] = simulator.pointer(b);
    usage[j][b] = simulator.usage(b);
  }

  getMaxUsage(j);
  marginat[j] = computeMinimumMargin(j);

  fail_flag = not simulator.excess_buffers.empty();

  return marginat[j];
}

template <typename TimeT, typename DataT>
double
Algorithm<TimeT, DataT>::playDownlinkWithoutTransfer(const int i, const int j,
                                                     const bool verbose) {

  simulator.initialise(pointer[i], usage[i]);

  if (verbose)
    simulator.print_flag = true;
  simulator.runWithoutTransfer(window[i], window[j], startOf(i), target);
  if (verbose)
    simulator.print_flag = false;

  for (auto b{0}; b <= data.numBuffers(); ++b) {
    pointer[j][b] = simulator.pointer(b);
    usage[j][b] = simulator.usage(b);
  }

  getMaxUsage(j);
  marginat[j] = computeMinimumMargin(j);

  return marginat[j];
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::playDownlinkWithParallelTransfer(
    const int i, const int j, const bool verbose) {

  simulator.initialise(pointer[i], usage[i]);

  if (verbose)
    simulator.print_flag = true;
  simulator.runWithParallelTransfer(window[i], window[j], startOf(i));
  if (verbose)
    simulator.print_flag = false;

  for (auto b{0}; b <= data.numBuffers(); ++b) {
    pointer[j][b] = simulator.pointer(b);
    usage[j][b] = simulator.usage(b);
  }

  getMaxUsage(j);
  marginat[j] = computeMinimumMargin(j);

  return marginat[j];
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::playSolution(vector<vector<vector<int>>> &prio,
                                             const bool verbose) {

  // cout << "PLAY SOLUTION\n" ;

  setTargetMargin(0);
  clearHandover();

  // simulator.print_flag = true;

  double minimum_margin{1};

  auto n{data.numDownlinks()};
  for (auto i{0}; i < n; ++i) {

    double margin =
        playDownlink(i, i + 1, prio[i], verbose); // sol->priority[i]);
    if (margin < minimum_margin)
      minimum_margin = margin;

    if (not simulator.excess_buffers.empty()) {

      assert(fail_flag);

      // cout << simulator.excess_buffers << endl;

      // cout << "FAIL AT WINDOW " << i << endl;
      // exit(1);
    } else {
      assert(not fail_flag);
    }
  }

  double margin{playDownlinkWithoutTransfer(n, n + 1, verbose)};

  // simulator.initialise(pointer[n], usage[n]);
  // simulator.runWithoutTransfer(window[n], window[n + 1], startOf(n),
  // target);
  // getMaxUsage(n + 1);
  //
  // double margin = computeMinimumMargin(n + 1);
  if (margin < minimum_margin)
    minimum_margin = margin;

  return minimum_margin;
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::initialise() {
  // fill functions
  for (auto b{0}; b < data.numBuffers(); ++b) {
    target.push_back(capacity(b));
    for (auto e{data.event_begin(b)}; e != data.event_end(b); ++e) {
      events.push_back({e->time, b});
    }
  }

  // downlink windows
  for (auto w{data.downlink_begin()}; w != data.downlink_end(); ++w) {
    events.push_back({w->time, static_cast<int>(data.numBuffers())});
  }

  sort(events.begin(), events.end());

  auto e{events.begin()};
  auto w{data.downlink_begin() + 2};

  while (e != events.end() and e->time == 0)
    ++e;
  window.push_back(e);

  // cout << "look for first event of window " << window.size() << ": " <<
  // w->time << endl;
  for (; e != events.end(); ++e) {
    if (w < data.downlink_end() and w->time < e->time) {
      window.push_back(e);
      ++w;
      ++w;

      // // cout << "-> " << e->time << endl;
      // if(e!=events.end() and w < data.downlink_end())
      // 	cout << "look for first event of window " << window.size() << ":
      // " << w->time << endl;
    }
    // else if(e!=events.end())
  }

  window.push_back(events.end());

  // cout << window.size() << " / " << (data.numDownlinks()) << endl;

  // assert(window.size() == data.numDownlinks() + 2);
  if (window.size() < data.numDownlinks() + 2) {
    window.push_back(events.end());
    // there's no event after the last downlink
  }

  marginat.resize(data.numDownlinks() + 2, 0);

  pointer.resize(data.numDownlinks() + 2);
  pointer[0].resize(data.numBuffers() + 1, 0);
  usage.resize(data.numDownlinks() + 2);
  maxusage.resize(data.numDownlinks() + 2);
  handover.resize(data.numDownlinks() + 1);
  // handover[0].resize(data.numBuffers(), capacity()
  //                    Simulator<TimeT, DataT>::infinite);

  for (auto b{data.buffer_begin()}; b != data.buffer_end(); ++b) {
    usage[0].push_back(b->initial_usage);
    maxusage[0].push_back(b->initial_usage);
  }
  usage[0].resize(data.numBuffers() + 1, 0);
  maxusage[0].resize(data.numBuffers() + 1, 0);

  for (auto w{1}; w <= data.numDownlinks() + 1; ++w) {
    maxusage[w].resize(data.numBuffers() + 1, 0);
    usage[w].resize(data.numBuffers() + 1, 0);
    pointer[w].resize(data.numBuffers() + 1, 0);
    // handover[w].resize(data.numBuffers(),
    //                    Simulator<TimeT, DataT>::infinite);
  }

  priority.resize(data.numDownlinks() + 1);

  nolimit.resize(data.numBuffers(), Simulator<TimeT, DataT>::infinite);

  thebuffers.reserve(data.numBuffers());

  clearHandover();

  simulator.initialise(pointer[0], usage[0]);
  double minimum_margin{1};

  for (auto b{0}; b < data.numBuffers(); ++b) {
    double margin = static_cast<double>(capacity(b) - usage[0][b]) /
                    static_cast<double>(capacity(b));
    if (margin < minimum_margin)
      minimum_margin = margin;
  }
  marginat[0] = minimum_margin;

  // fill.resize(data.numDownlinks() + 1);
  // for (auto w{0}; w <= data.numDownlinks(); ++w) {
  //   fill[w].resize(data.numBuffers(), 0);
  // }
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::clearHandover() {
  for (auto w{0}; w <= data.numDownlinks(); ++w) {
    handover[w].clear();
    for (auto b{0}; b < data.numBuffers(); ++b)
      handover[w].push_back(capacity(b));
  }
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::computeMinimumMargin(const int i) const {
  double minimum_margin{1};
  for (auto b{0}; b < data.numBuffers(); ++b) {
    assert(maxusage[i][b] == simulator.getMaxUsage(b));

    double margin{marginOf(i, b)};
    //     //
    //     //
    //     //
    //     double margin =
    //         static_cast<double>(capacity(b) - simulator.getMaxUsage(b)) /
    //         static_cast<double>(capacity(b));
    //
    // assert(margin == )

    if (margin < minimum_margin)
      minimum_margin = margin;
  }
  return minimum_margin;
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::getMaxUsage(const int i) {

  // cout << maxusage.size() << endl;

  for (auto b{0}; b < data.numBuffers(); ++b) {

    // assert(i>=0);
    // assert(i<maxusage.size());
    // assert(b>=0);
    // assert(b<maxusage[i].size());
    maxusage[i][b] = simulator.getMaxUsage(b);
  }
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::computeDCPriority(const int i) {
  best_priority.clear();
  best_priority.resize(1);
  int p{0};
  thebuffers.clear();

  // run from the current point without downlinks to compute the loss ordering
  for (auto j{i}; j < data.numDownlinks(); ++j) {

    // cout << "probe " << j ;

    // double margin{
    playDownlinkWithoutTransfer(j, j + 1);
    // };

    // simulator.initialise(pointer[j], usage[j]);
    // simulator.runWithoutTransfer(window[j], window[j + 1], startOf(j),
    // target);
    //
    // for (auto b{0}; b <= data.numBuffers(); ++b) {
    //   pointer[j + 1][b] = simulator.pointer(b);
    //   usage[j + 1][b] = simulator.usage(b);
    // }

    auto save{thebuffers.count()};
    for (auto b : simulator.excess_buffers) {
      if (not thebuffers.has(b)) {
        thebuffers.add(b);
        best_priority[p].push_back(b);
      }
    }

    // cout << " " << p << " / " << buffers << endl;

    if (thebuffers.count() == data.numBuffers())
      break;

    if (save < thebuffers.count()) {
      best_priority.resize(++p + 1);
    }
  }
  if (thebuffers.count() < data.numBuffers()) {
    for (auto b{thebuffers.bbegin()}; b != thebuffers.bend(); ++b) {
      best_priority[p].push_back(*b);
    }
  }

  priority[i].clear();
  for (auto batch{best_priority.rbegin()}; batch != best_priority.rend();
       ++batch) {
    priority[i].push_back(*batch);
  }
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::computeFLPriority(const int i) {
  best_priority.clear();
  // best_priority.resize(1);

  // run from the current point without downlinks to compute the loss ordering
  playDownlinkWithoutTransfer(i, data.numDownlinks() + 1);
  // simulator.initialise(pointer[i], usage[i]);
  // simulator.runWithoutTransfer(window[i], window.back(), startOf(i),
  // target);

  // extract the priority ordering
  // for (auto &batch : best_priority) {
  //   batch.clear();
  // }
  best_priority.resize(simulator.excess_buffers.count() + 1);
  // (simulator.excess_buffers.count() < data.numBuffers() ? 1 : 0));

  // if (simulator.excess_buffers.count() < data.numBuffers()) {
  for (auto b{simulator.excess_buffers.end()};
       b != simulator.excess_buffers.bend(); ++b) {
    best_priority[0].push_back(*b);
  }
  int k{0};
  for (auto b{simulator.excess_buffers.rbegin()};
       b != simulator.excess_buffers.rend(); ++b) {
    best_priority[++k].push_back(*b);
  }
  // }

  priority[i] = best_priority;
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::randomizedDownlinkCountHeuristic(
    const double best_margin, const double target_margin, int &nruns) {
  // for (auto i{0}; i < nruns; ++i) {

  // double best{best_margin};
  while (nruns > 0) {
    auto margin{downlinkCountHeuristic(target_margin, -1)};

		if (option.verbosity >= 1)
			cout << "r " << margin << " / " << target_margin << endl;

    if (margin > best_margin)
      return margin;
    --nruns;
  }
  return best_margin - .00001;
}

template <typename TimeT, typename DataT>
Solution *Algorithm<TimeT, DataT>::randomWalk() {
	
	// cout << "HELLO (" << option.verbosity << ")\n";

  Solution *best{NULL};

  double ubmargin{relaxBandwidth()};

  double lbmargin{downlinkCountHeuristic(0, 0)};
  cout << "obj=" << (1 - real_margin(lbmargin)) * 100 << endl;
  cout << "time=" << cpu_time()*1000 << endl;
  best = new Solution(priority);

  while (cpu_time()*1000 < option.htime_limit) {
    auto maxiter = option.runs;
    auto margin{randomizedDownlinkCountHeuristic(lbmargin, 0, maxiter)};

    if (margin > lbmargin) {
      lbmargin = margin;
      // if (option.verbosity >= 0)
      cout << "obj=" << (1 - real_margin(lbmargin)) * 100 << endl;
      cout << "time=" << cpu_time()*1000 << endl;
      // delete best;
      // best = new Solution(priority);
      best->copy(priority);
    }
  }

  return best;
}

template <typename TimeT, typename DataT>
double
Algorithm<TimeT, DataT>::downlinkCountHeuristic(const double target_margin,
                                                const int flag) {

  if (flag >= 0)
    clearHandover();
  priority.clear();
  priority.resize(data.numDownlinks() + 1);

  setTargetMargin(target_margin);

  double minimum_margin = 1.0;

  for (auto i = 0; i < data.numDownlinks(); ++i) {

    int choice{flag};
    if (choice != 0 and choice != 1)
      choice = random() % 2;

    if (choice == 1) {
			// cout << "fl\n";
      computeFLPriority(i);
    } else if (choice == 0) {
			// cout << "dc\n";
      computeDCPriority(i);
    }

    auto margin{playDownlink(i, i + 1, priority[i])};

    if (margin < minimum_margin) {
      minimum_margin = margin;
    }

    // if (option.verbosity > 0) {
    //   cout << setw(10) << minimum_margin << " " <<
    //   real_margin(minimum_margin)
    //        << "\n";
    // }

    if (minimum_margin < 0)
      break;
  }

  // run residual fills without transfer
  auto n{data.numDownlinks()};

  // if (minimum_margin >= 0) {
  // simulator.runWithoutTransfer(window[n], window[n + 1], startOf(n),
  // target);
  //
  // getMaxUsage(n + 1);
  // auto margin{computeMinimumMargin(n + 1)};
  auto margin{playDownlinkWithoutTransfer(n, n + 1)};

  if (margin < minimum_margin)
    minimum_margin = margin;

  // if (option.verbosity > ) {
  //   cout << setw(10) << minimum_margin << " " <<
  //   real_margin(minimum_margin)
  //        << "\n";
  // }
  // }

  return minimum_margin;
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::relaxBandwidth() {

  double minimum_margin{1};
  for (auto b{0}; b < data.numBuffers(); ++b) {

    if (option.verbosity > 1)
      cout << "run " << b << " at maximum priority\n";

    best_priority.clear();
    best_priority.resize(2);
    for (auto a{0}; a < data.numBuffers(); ++a) {
      if (a != b)
        best_priority[0].push_back(a);
    }
    best_priority[1].push_back(b);

    // cout << "hi " << b << endl;
    // assert(handover.back().size() == data.numBuffers());

    // if (option.verbosity > 3 or b == 12)
    //   simulator.print_flag = true;
    simulator.initialise(pointer[0], usage[0]);
    simulator.run(window[0], window.back(), 0, best_priority, target,
                  handover.back());

    // if (option.verbosity < 5 or b == 12)
    //   simulator.print_flag = false;

    auto margin{static_cast<double>(capacity(b) - simulator.getMaxUsage(b)) /
                static_cast<double>(capacity(b))};

    if (option.verbosity > 1)
      cout << "margin=" << margin << "\n";

    // cout << b << ": " << margin << endl;

    if (margin < minimum_margin)
      minimum_margin = margin;
  }

  // cout << " ==> " << minimum_margin << endl;

  return minimum_margin;
}

template <typename TimeT, typename DataT>
void Algorithm<TimeT, DataT>::postHandoverConstraints(const double margin) {

  auto n{data.numDownlinks()};

  for (auto i{0}; i < n; ++i) {
    for (auto b{0}; b < data.numBuffers(); ++b) {
      // handover[i][b] =
      //     min(handover[i][b], capacity(b) * (1.0 - margin - hdelta[i +
      //     1][b]));
      handover[i][b] = capacity(b) * (1.0 - margin - hdelta[i + 1][b]);
      // min(handover[i][b], );
    }
  }

  // for (auto i{0}; i < n; ++i) {
  //   cout << setw(3) << i << " @" << setw(9)
  //        << static_cast<long int>(data.getDownlinkStart(i)) << " : (" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkRate(i)) << "*" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkDuration(i)) << ")";
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12)
  // 				   << (handover[i][b] / capacity(b))
  //          // << static_cast<long int>(handover[i][b]);
  //   }
  //   cout << endl;
  // }

  // cout << "handover of 3 @57 = " << handover[57][3] << " ("
  //      << (handover[57][3] / capacity(3)) << ")\n";
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::getHandoverDeltas() {
  double minimum_margin{1};

  auto n{data.numDownlinks()};

  // // hdelta[i][b] is the minimum increase of the usage (relative to
  // capacity)
  // // from the start to the end of window i
  // hdelta.resize(n + 1);
  // // hdelta[0].resize()
  // for (auto i{0}; i <= n + 1; ++i)
  //   hdelta[i].resize(data.numBuffers(), 0);
  //
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //
  //   if (option.verbosity > 0)
  //     cout << "run " << b << " at maximum priority\n";
  //
  //   best_priority.clear();
  //   best_priority.resize(2);
  //   for (auto a{0}; a < data.numBuffers(); ++a) {
  //     if (a != b)
  //       best_priority[0].push_back(a);
  //   }
  //   best_priority[1].push_back(b);
  //
  //   for (auto i{0}; i < n; ++i) {
  //
  //     if (option.verbosity > 1)
  //       cout << "window " << i << endl;
  //
  //     auto margin{playDownlink(i, i + 1, best_priority)};
  //     if (margin < minimum_margin)
  //       minimum_margin = margin;
  //
  //     auto increase{(maxusage[i + 1][b] - usage[i][b])};
  //
  //     // if (increase > 0)
  //     hdelta[i][b] = increase / capacity(b);
  //     // else
  //     // hdelta[i][b] = 0;
  //
  //     // optimistic_usage[i + 1][b] = usage[i + 1][b];
  //   }
  //
  //   auto margin{playDownlinkWithoutTransfer(n, n + 1)};
  //   if (margin < minimum_margin)
  //     minimum_margin = margin;
  //
  //   if (maxusage[n + 1][b] > 0)
  //     hdelta[n][b] = maxusage[n + 1][b] / capacity(b);
  //   else
  //     hdelta[n][b] = 0;
  // }
  //
  // cout << "cap:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>(capacity(b));
  // }
  // cout << endl;
  // cout << "deb:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>(usage[0][b]);
  // }
  // cout << endl;
  // for (auto i{0}; i < n; ++i) {
  //   cout << setw(3) << i << " @" << setw(9)
  //        << static_cast<int>(data.getDownlinkStart(i)) << " : (" << setw(6)
  //        << static_cast<int>(data.getDownlinkRate(i)) << "*" << setw(6)
  //        << static_cast<int>(data.getDownlinkDuration(i)) << ")";
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << hdelta[i][b];
  //   }
  //   cout << endl;
  // }
  // cout << "fin:                            ";
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << hdelta[n][b];
  // }
  // cout << endl;
  //
  // cout << "deb:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>((1 - hdelta[0][b]) *
  //   capacity(b));
  // }
  // cout << endl;
  // for (auto i{0}; i < n; ++i) {
  //   cout << setw(3) << i << " @" << setw(9)
  //        << static_cast<long int>(data.getDownlinkStart(i)) << " : (" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkRate(i)) << "*" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkDuration(i)) << ")";
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12)
  //          << static_cast<long int>((1 - hdelta[i + 1][b]) * capacity(b));
  //   }
  //   cout << endl;
  // }
  // // cout << "fin:                            ";
  // // for (auto b{0}; b < data.numBuffers(); ++b) {
  // //   cout << setw(12)
  // //        << static_cast<long int>((1 - hdelta[n + 1][b]) * capacity(b));
  // // }
  // // cout << endl;
  //
  // // exit(1);

  vector<DataT> dummy(data.numBuffers() + 1, 0);
  // for (auto b{0}; b < data.numBuffers(); ++b)
  // 	usage[i][b] =

  // newdelta[i][b] is the minimum increase of the usage (relative to
  // capacity)
  // from the start to the end of window i

  // vector<vector<double>> newdelta;
  hdelta.resize(n + 1);
  // hdelta[0].resize()
  for (auto i{0}; i <= n; ++i)
    hdelta[i].resize(data.numBuffers(), 0);

  simulator.bottom = false;
  for (auto i{0}; i < n; ++i) {

    swap(dummy, usage[i]);

    for (auto b{0}; b < data.numBuffers(); ++b)
      usage[i][b] = 0;

    playDownlinkWithParallelTransfer(i, i + 1);

    for (auto b{0}; b < data.numBuffers(); ++b) {
      hdelta[i][b] = maxusage[i + 1][b] / capacity(b);
    }

    swap(dummy, usage[i]);
  }

  swap(dummy, usage[n]);
  playDownlinkWithParallelTransfer(n, n + 1);
  for (auto b{0}; b < data.numBuffers(); ++b) {
    hdelta[n][b] = maxusage[n + 1][b] / capacity(b);
  }
  swap(dummy, usage[n]);
  simulator.bottom = true;

  // cout << "cap:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>(capacity(b));
  // }
  // cout << endl;
  // cout << "deb:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>(usage[0][b]);
  // }
  // cout << endl;
  // for (auto i{0}; i < n; ++i) {
  //   cout << setw(3) << i << " @" << setw(9)
  //        << static_cast<int>(data.getDownlinkStart(i)) << " : (" << setw(6)
  //        << static_cast<int>(data.getDownlinkRate(i)) << "*" << setw(6)
  //        << static_cast<int>(data.getDownlinkDuration(i)) << ")";
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << hdelta[i][b];
  //   }
  //   cout << endl;
  // }
  // cout << "fin:                            ";
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << hdelta[n][b];
  // }
  // cout << endl;
  //
  // cout << "deb:                            " << setprecision(3) << fixed;
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   cout << setw(12) << static_cast<long int>((1 - hdelta[0][b]) *
  //   capacity(b));
  // }
  // cout << endl;
  // for (auto i{0}; i < n; ++i) {
  //   cout << setw(3) << i << " @" << setw(9)
  //        << static_cast<long int>(data.getDownlinkStart(i)) << " : (" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkRate(i)) << "*" <<
  //        setw(6)
  //        << static_cast<long int>(data.getDownlinkDuration(i)) << ")";
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12)
  //          << static_cast<long int>((1 - hdelta[i + 1][b]) * capacity(b));
  //   }
  //   cout << endl;
  // }
  // // cout << "fin:                            ";
  // // for (auto b{0}; b < data.numBuffers(); ++b) {
  // //   cout << setw(12)
  // //        << static_cast<long int>((1 - hdelta[n + 1][b]) * capacity(b));
  // // }
  // // cout << endl;
  //
  // // exit(1);

  return minimum_margin;
}

template <typename TimeT, typename DataT>
double
Algorithm<TimeT, DataT>::randomizedMultipleWindows(const double allowed_margin,
                                                   const int nruns) {
  for (auto i{0}; i < nruns; ++i) {
    auto margin{multipleWindows(allowed_margin)};
    if (margin > allowed_margin)
      return margin;
  }
  return -static_cast<double>(data.numBuffers());
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::multipleWindows(const double allowed_margin) {

  setTargetMargin(allowed_margin);
  clearHandover();

  int i{0};

  DataT tcapa{0};
  for (auto b{0}; b < data.numBuffers(); ++b)
    tcapa += capacity(b);

  // auto iter{100};

  double minimum_margin{1};

  while (i <= data.numDownlinks()) {

    if (singleWindow(i)) {

      debug_flag = false;

      // cout << " " << i;

      if (option.verbosity > 0)
        cout << "window " << i << " ok at " << allowed_margin << endl;

      priority[i].clear();
      priority[i] = best_priority;

      double margin{playDownlink(i, i + 1, priority[i])};
      if (margin < minimum_margin) {
        // cout << i << ": " << minimum_margin << " -> " << margin << endl;
        minimum_margin = margin;
      }

      ++i;

    } else {

      debug_flag = false;

      if (option.verbosity > 0) {
        cout << "window " << i << " infeasible at " << allowed_margin << endl;
        cout << setw(31) << " ";
        for (auto b{0}; b < data.numBuffers(); ++b) {
          cout << setw(6)
               << static_cast<int>(1000.0 * handover[i][b] / capacity(b));
        }
        cout << endl;
      }

      if (i == 0) {
        // cout << " fail\n";
        return -static_cast<double>(data.numBuffers());
      }

      if (option.verbosity > 0)
        cout << simulator.excess_buffers << endl;

      int b =
          simulator.excess_buffers[random() % simulator.excess_buffers.size()];

      if (option.verbosity > 0)
        cout << "data loss on buffer " << b << endl;

      auto delta(min(target[b] - simulator.getMaxUsage(b),
                     handover[i][b] - simulator.usage(b)));
      handover[i - 1][b] = (usage[i][b] + delta);

      if (option.verbosity > 0)
        cout << "prev window " << (i - 1) << endl;

      // if(i == 25 and b == 27) {
      // 	cout << "DEBUG!\n";
      // 		debug_flag = true;
      // }

      --i;
    }
  }

  // cout << endl;
  return minimum_margin;
}

template <typename TimeT, typename DataT>
Solution *Algorithm<TimeT, DataT>::maxMultipleWindows() {

  auto proved_ub{relaxBandwidth()};
  // double proved_ub{1};

  if (option.verbosity >= 0)
    cout << "ub " << real_margin(proved_ub) << " " << cpu_time()*1000 << " \n";

  auto lb{downlinkCountHeuristic()};
  // if (lb < 0)
  //   lb = 0;

  Solution *best_solution = new Solution(priority);

  // if (option.verbosity >= 0)
  //   cout << "sol " << real_margin(lb) << " " << cpu_time()*1000 << endl;
  cout << "obj=" << (1 - real_margin(lb)) * 100 << endl;
  cout << "time=" << cpu_time()*1000 << endl;

  // double lb{0};
  double ub{proved_ub};

  if (option.verbosity > 0)
    cout << lb << " .. " << ub << "\n";

  double epsilon{0.00001};

  while (lb + epsilon < ub) {

    double target_margin{(lb + ub) / 2.0};

    // double actual_margin{multipleWindows(target_margin)};

    // cout << target_margin << endl;

    double actual_margin{randomizedMultipleWindows(target_margin, option.runs)};

    // cout << "actual_margin = " << actual_margin << endl;

    if (actual_margin >= target_margin) {
      lb = actual_margin;
      // delete best_solution;
      // best_solution = new Solution(priority);
      best_solution->copy(priority);
      if (option.verbosity >= 0)
        cout << "sol " << real_margin(lb) << " " << cpu_time()*1000 << endl;
      cout << "obj=" << (1 - real_margin(lb)) * 100 << endl;
      cout << "time=" << cpu_time()*1000 << endl;
    } else
      ub = target_margin;

    if (option.verbosity > 0)
      cout << lb << " .. " << ub << "\n";
  }

  // for(auto margin{lb}; margin<proved_ub; margin+=.0001) {
  //
  // 	double actual_margin{multipleWindows(margin)};
  // 	if(actual_margin > lb) {
  // 		cout << "IMPROVE LB!!!: " << actual_margin << "\n";
  // 		lb = actual_margin;
  // 	}
  // 	// else {
  // 	// 	cout << margin << " -> " << actual_margin << endl;
  // 	//  		}
  // }
  //
  // cout << "HERE\n";

  // for(auto margin{lb}; margin<proved_ub; margin+=.0001) {
  //
  // 	double actual_margin{downlinkCountHeuristic(margin)};
  // 	if(actual_margin > lb) {
  // 		cout << "IMPROVE LB!!!: " << actual_margin << "\n";
  // 		lb = actual_margin;
  // 	}
  // 	// else {
  // 	// 	cout << margin << " -> " << actual_margin << endl;
  // 	//  		}
  // }

  return best_solution;
}

template <typename TimeT, typename DataT>
Solution *Algorithm<TimeT, DataT>::iterativeLeveling() {

  double margin{0};
  double proved_ub{relaxBandwidth()};

  // if (option.verbosity >= 0)
  //   cout << "ub " << real_margin(proved_ub) << " " << cpu_time()*1000 << " \n";

  double delta{proved_ub};
  double prev{delta};

  double best{downlinkCountHeuristic()};
  // if (option.verbosity >= 0)
  //   cout << "sol " << real_margin(best) << " " << cpu_time()*1000 << " \n";

  cout << "obj=" << (1 - real_margin(best)) * 100 << endl;
  cout << "time=" << cpu_time()*1000 << endl;

  Solution *best_solution = new Solution(priority);

  // Solution *best_solution{NULL};

  while (cpu_time()*1000 < option.htime_limit) {

    // cout << "call heuristic(" << margin << ")\n";

    // if (option.verbosity > 0)
    //   cout << "margin=" << margin << " delta=" << delta << endl;

    // auto result{downlinkCountHeuristic(margin)};
    int r{option.runs};
    auto result{randomizedDownlinkCountHeuristic(best, margin, r)};
    // auto result{firstLossHeuristic(margin)};
    // best = max(best, result);
    if (result > best) {
      // if (best_solution)
      // delete best_solution;
      // best_solution = new Solution(priority);
      best_solution->copy(priority);
      best = result;
      cout << "obj=" << (1 - real_margin(best)) * 100 << endl;
      cout << "time=" << cpu_time()*1000 << endl;
    }

    delta = abs(margin - prev) / 2;
    prev = margin;

    // cout << "result=" << result << "; best=" << best << "; delta=" << delta
    // << "; prev=" << prev << endl;

    if (result > margin) {
      margin += delta;
    } else {
      margin -= delta;
    }
  }

  // for (auto margin{best}; margin < proved_ub; margin += .0001) {
  //
  //   double actual_margin{downlinkCountHeuristic(margin)};
  //   if (actual_margin > best) {
  //     cout << "IMPROVE LB!!!: " << actual_margin << "\n";
  //     best = actual_margin;
  //     if (best_solution)
  //       delete best_solution;
  //     best_solution = new Solution(priority);
  //   }
  //   // else {
  //   // 	cout << margin << " -> " << actual_margin << endl;
  //   //  		}
  // }

  // if(best >= 0) {
  // 	auto verification = playSolution(best_solution);
  // 	assert(verification == best);
  // }

  return best_solution;
}

template <typename TimeT, typename DataT>
double Algorithm<TimeT, DataT>::singleWindowMax(const int i) {

  if (debug_flag)
    cout << "\nOPTIMIZE WINDOW " << i << endl;

  setTargetMargin(0);

  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   assert(nolimit[b] > capacity(b));
  //   target[b] = capacity(b);
  // }

  if (debug_flag)
    cout << "\nfeasibility test at full capacity\n";
  if (not singleWindow(i)) {
    if (debug_flag)
      cout << "INFEASIBLE!\n";
    return -static_cast<double>(data.numBuffers());
  }

  double ub{1};
  double lb{0};

  while (lb + 0.001 < ub) {
    double m{(lb + ub) / 2.0};

    // for (auto b{0}; b < data.numBuffers(); ++b) {
    //   assert(nolimit[b] > capacity(b));
    //   target[b] = capacity(b) * (1.0 - m);
    // }
    setTargetMargin(m);

    if (debug_flag)
      cout << "\n [" << lb << ".." << ub << "] dichotomic call with target "
           << m << "\n";
    if (singleWindow(i)) {
      lb = m;

      priority[i].clear();
      priority[i] = best_priority;

    } else
      ub = m;
  }

  if (debug_flag) {
    cout << "\n run final prio on window " << i << ":\n";
    auto batch{priority[i].rbegin()};
    while (batch != priority[i].rend()) {
      cout << "p" << (batch - priority[i].rbegin()) << ": ";
      for (auto b : *batch)
        cout << " " << b;
      cout << endl;
      ++batch;
    }
  }

  setTargetMargin(lb);
  // for (auto b{0}; b < data.numBuffers(); ++b) {
  //   assert(nolimit[b] > capacity(b));
  //   target[b] = capacity(b) * (1.0 - lb);
  // }

  // simulator.setAllActive();

  simulator.initialise(pointer[i], usage[i]);
  simulator.run(window[i], window[i + 1], startOf(i), priority[i], target,
                handover[i]);

  for (auto b{0}; b <= data.numBuffers(); ++b) {
    pointer[i + 1][b] = simulator.pointer(b);
    usage[i + 1][b] = simulator.usage(b);
  }

  if (debug_flag)
    for (auto b{0}; b < data.numBuffers(); ++b) {
      cout << "usage: " << usage[i][b] << " -> " << usage[i + 1][b]
           << " (cap. = " << capacity(b) << ")"
           << " pointer: " << pointer[i][b] << " -> " << pointer[i + 1][b]
           << " dl pointer " << pointer[i][data.numBuffers()] << " -> "
           << pointer[i + 1][data.numBuffers()] << endl;
    }

  for (auto b{0}; b < data.numBuffers(); ++b) {
    if (simulator.getMaxUsage(b) > target[b]) {
      cout << "target exceeded at window " << i << " for " << b << ": "
           << simulator.getMaxUsage(b) << "/" << target[b] << endl;
      exit(1);
    }
  }

  return lb;
}

// // tries to minimize the handover of buffer b while ensuring that the
// template <typename TimeT, typename DataT>
// double Algorithm<TimeT, DataT>::singleWindowMax(const int i) {
//
//   if (debug_flag)
//     cout << "\nOPTIMIZE WINDOW " << i << endl;
//
//   setTargetMargin(0);
//
//   // for (auto b{0}; b < data.numBuffers(); ++b) {
//   //   assert(nolimit[b] > capacity(b));
//   //   target[b] = capacity(b);
//   // }
//
//   if (debug_flag)
//     cout << "\nfeasibility test at full capacity\n";
//   if (not singleWindow(i)) {
//     if (debug_flag)
//       cout << "INFEASIBLE!\n";
//     return -1;
//   }
//
//   double ub{1};
//   double lb{0};
//
//   while (lb + 0.001 < ub) {
//     double m{(lb + ub) / 2.0};
//
//     // for (auto b{0}; b < data.numBuffers(); ++b) {
//     //   assert(nolimit[b] > capacity(b));
//     //   target[b] = capacity(b) * (1.0 - m);
//     // }
//     setTargetMargin(m);
//
//     if (debug_flag)
//       cout << "\n [" << lb << ".." << ub << "] dichotomic call with target
//       "
//            << m << "\n";
//     if (singleWindow(i)) {
//       lb = m;
//
//       priority[i].clear();
//       priority[i] = best_priority;
//
//     } else
//       ub = m;
//   }
//
//   if (debug_flag) {
//     cout << "\n run final prio on window " << i << ":\n";
//     auto batch{priority[i].rbegin()};
//     while (batch != priority[i].rend()) {
//       cout << "p" << (batch - priority[i].rbegin()) << ": ";
//       for (auto b : *batch)
//         cout << " " << b;
//       cout << endl;
//       ++batch;
//     }
//   }
//
//   setTargetMargin(lb);
//   // for (auto b{0}; b < data.numBuffers(); ++b) {
//   //   assert(nolimit[b] > capacity(b));
//   //   target[b] = capacity(b) * (1.0 - lb);
//   // }
//
//   // simulator.setAllActive();
//
//   simulator.initialise(pointer[i], usage[i]);
//   simulator.run(window[i], window[i + 1], startOf(i), priority[i], target,
//                 handover[i]);
//
//   for (auto b{0}; b <= data.numBuffers(); ++b) {
//     pointer[i + 1][b] = simulator.pointer(b);
//     usage[i + 1][b] = simulator.usage(b);
//   }
//
//   if (debug_flag)
//     for (auto b{0}; b < data.numBuffers(); ++b) {
//       cout << "usage: " << usage[i][b] << " -> " << usage[i + 1][b]
//            << " (cap. = " << capacity(b) << ")"
//            << " pointer: " << pointer[i][b] << " -> " << pointer[i + 1][b]
//            << " dl pointer " << pointer[i][data.numBuffers()] << " -> "
//            << pointer[i + 1][data.numBuffers()] << endl;
//     }
//
//   for (auto b{0}; b < data.numBuffers(); ++b) {
//     if (simulator.getMaxUsage(b) > target[b]) {
//       cout << "target exceeded at window " << i << " for " << b << ": "
//            << simulator.getMaxUsage(b) << "/" << target[b] << endl;
//       exit(1);
//     }
//   }
//
//   return lb;
// }

template <typename TimeT, typename DataT>
bool Algorithm<TimeT, DataT>::singleWindow(const int i) {

  auto inf{Simulator<TimeT, DataT>::infinite};

  for (auto b{0}; b < data.numBuffers(); ++b) {
    assert(target[b] < nolimit[b]);
    assert(nolimit[b] == inf);
  }

  // debug_flag = false;

  if (debug_flag)
    cout << "\n at window " << i << endl;

  best_priority.clear();
  best_priority.reserve(data.numBuffers());
  best_priority.resize(1);
  for (auto b{0}; b < data.numBuffers(); ++b)
    best_priority[0].push_back(b);

  int p{0};
  while (p < best_priority.size() and not best_priority[p].empty()) {

    for (auto b : best_priority[p]) {
      if (p) {
        // assert(target[b] > capacity(b));
        // assert(nolimit[b] <= capacity(b));
        assert(target[b] == inf);
        assert(nolimit[b] != inf);
        swap(target[b], nolimit[b]);
      } else {
        assert(target[b] != inf);
      }
      // else {
      //         if (target[b] > capacity(b)) {
      //
      //           cout << i << " " << b << " " << target[b] << " " <<
      //           capacity(b)
      //                << endl;
      //           exit(1);
      //         }
      //       }
    }

    if (debug_flag)
      cout << "solve at prio " << p << endl;

    do {

      simulator.initialise(pointer[i], usage[i]);
      if (debug_flag)
        cout << "downlink pointer = " << pointer[i][data.numBuffers()] << endl;

      if (debug_flag) {
        cout << "run for buffers of priority " << p << endl;
        auto batch{best_priority.rbegin()};
        while (batch != best_priority.rend()) {
          cout << "p" << (batch - best_priority.rbegin()) << ": ";
          for (auto b : *batch)
            cout << " " << b;
          cout << endl;
          ++batch;
        }
      }

      // if (debug_flag)
      //   simulator.print_flag = true;
      simulator.run(window[i], window[i + 1], startOf(i), best_priority, target,
                    handover[i]);

      // simulator.print_flag = false;

      bool no_overflow{true};
      for (auto b{best_priority[p].begin()}; b < best_priority[p].end();) {
        if (simulator.excess_buffers.has(*b)) {

          if (debug_flag) {
            cout << "raise prio of " << *b << " tgt=" << target[*b]
                 << " nlt=" << nolimit[*b] << endl;
          }

          auto buf{*b};

          swap(target[*b], nolimit[*b]);

          best_priority.resize(p + 2);
          best_priority[p + 1].push_back(*b);

          swap(*b, *best_priority[p].rbegin());
          assert(best_priority[p].back() == buf);
          best_priority[p].pop_back();
          no_overflow = false;

          assert(target[buf] > capacity(buf));

        } else
          ++b;
      }

      if (debug_flag) {
        for (auto xp{0}; xp < best_priority.size(); ++xp) {
          if (xp == p)
            cout << "**";
          else
            cout << "  ";
          for (auto b : best_priority[xp])
            cout << " " << b;
          cout << endl;
        }
      }

      if (no_overflow) {
        if (debug_flag) {
          cout << "no overflow\n";
        }
        break;
      }

      if (best_priority[p].empty()) {
        if (debug_flag) {
          cout << "UNSAT!\n";
        }

        for (auto b{best_priority[p + 1].begin()};
             b < best_priority[p + 1].end(); ++b) {
          // assert(target[*b] > capacity(*b));
          assert(target[*b] == inf);
          swap(target[*b], nolimit[*b]);
        }

        for (auto b{0}; b < data.numBuffers(); ++b) {
          assert(target[b] != inf);
          assert(nolimit[b] == inf);
          // assert(target[b] <= capacity(b));
          // assert(nolimit[b] > capacity(b));
        }

        debug_flag = false;
        return false;
      }

    } while (true);

    ++p;
  }

  debug_flag = false;
  return true;
}

template <typename TimeT, typename DataT>
ostream &Algorithm<TimeT, DataT>::display(ostream &os) const {
  os << "hello there";
  return os;
}

template <typename TimeT, typename DataT>
ostream &operator<<(ostream &os, const Algorithm<TimeT, DataT> &i) {
  return i.display(os);
}
}

#endif // __Algorithm_HPP
