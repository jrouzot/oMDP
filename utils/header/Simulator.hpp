
#ifndef __Simulator_HPP
#define __Simulator_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <math.h>

#include "Instance.hpp"
#include "SparseSet.hpp"
#include "Event.hpp"
#include "State.hpp"

using namespace std;

#define PRINTDEBUG printState(false);

namespace dataflow {

/**********************************************
* Simulator
**********************************************/

template <typename U, typename T> T round_and_cast(const U x) {
  if (not is_integral<U>::value and is_integral<T>::value) {
    return static_cast<T>(ceil(x));
  }
  return static_cast<T>(x);
}

template<typename TimeT, typename DataT>
class Simulator
{
public:
  // template <typename Container, typename Timepoint, typename Ranking>
  // DataT run(Container &buffers, Timepoint be, Timepoint ee, TimeT start_time,
  //           Ranking &priority, vector<DataT> &target);
  // template <typename Container, typename Timepoint>
  // void runWithoutTransfer(Container &buffers, Timepoint begin_event,
  //                        Timepoint end_event, TimeT start_time,
  //                        vector<DataT> &target);

  // template <typename Timepoint, typename Ranking>
  // bool runRelaxed(Timepoint be, Timepoint ee, TimeT start_time, Ranking &priority, vector<int> subsetBuffer);
  
  template <typename Timepoint, typename Ranking>
  void run(Timepoint be, Timepoint ee, TimeT start_time, Ranking &priority,
           vector<DataT> &target, vector<DataT> &htarget);
  template <typename Timepoint>
  void runWithoutTransfer(Timepoint begin_event, Timepoint end_event,
                          TimeT start_time, vector<DataT> &target);

  // used to find upper bounds
  template <typename Timepoint>
  void runWithParallelTransfer(Timepoint begin_event, Timepoint end_event,
                               TimeT start_time);

  // constexpr static const DataT infinite = numeric_limits<DataT>::max();
  // constexpr static const TimeT never = numeric_limits<TimeT>::max();
  // constexpr static const DataT epsilon = static_cast<DataT>(0.001);
  // constexpr static const DataT infinite; //= numeric_limits<DataT>::max();
  // constexpr static const TimeT never; //= numeric_limits<TimeT>::max();
  // constexpr static const DataT epsilon; // = static_cast<DataT>(0.001);
  static const DataT infinite; //= numeric_limits<DataT>::max();
  static const TimeT never;    //= numeric_limits<TimeT>::max();
  static const DataT epsilon;  // = static_cast<DataT>(0.001);

  Simulator(const Instance<TimeT, DataT> &i);

  // template <typename Container>
  void initialise(vector<int> &initial_pointers, vector<DataT> &initial_usage);

  // void computeMargin();
  TimeT currentTime() const;
  template <typename Ranking> void computeRates(Ranking &priority);
  void advanceTo(const TimeT t);

  SparseSet<> excess_buffers;

  ostream &display(ostream &os) const;

  DataT usage(const int b) const;
  int pointer(const int b) const;

  void printState(const bool all = false) const;

  DataT getMaxUsage(const int b) const;

  // void setAllActive();
  //
  // template <typename Container> void setActive(Container &buffers);

  // void setAllChecked() const;
  //
  // template <typename Container>
  // void setChecked(Container &buffers);

  bool print_flag{false};
  bool debug_flag{false};
  int debug_buffer{-1};
  bool bottom{true};

  void verification() const;

private:
  const Instance<TimeT, DataT> &data;

  vector<DataT> capacity;
  vector<State<DataT>> current;
  TimeT now{0};

  vector<State<DataT>> parallel;

  // vector<vector<State<DataT>>> saved_state;
  // vector<TimeT> saved_time;

  vector<DataT> handover;

  vector<DataT> max_usage;
  // DataT overall_max_usage;

  // buffers that are run in the simulation (must be a subset of
  // {0,...,data.numBuffers()-1})
  // SparseSet<> active_buffers;

  // // buffers that are checked for overflow
  // SparseSet<> checked_buffers;

  // whether we stop at the first overlfow
  // bool stop_early{false};

  TimeT getEmptyTime(const int b) const;
  DataT fillRate(const int b) const;
  DataT transferRate(const int b) const;

  DataT dumpRate() const;

  void getOverflow();
};

template <typename TimeT, typename DataT>
const DataT Simulator<TimeT, DataT>::infinite = numeric_limits<DataT>::max();

template <typename TimeT, typename DataT>
const TimeT Simulator<TimeT, DataT>::never = numeric_limits<TimeT>::max();

template <typename TimeT, typename DataT>
const DataT Simulator<TimeT, DataT>::epsilon = static_cast<DataT>(0.001);

template <typename TimeT, typename DataT>
DataT Simulator<TimeT, DataT>::getMaxUsage(const int b) const {
  return max_usage[b];
}

// template <typename TimeT, typename DataT>
// void Simulator<TimeT, DataT>::setAllActive() {
//   active_buffers.fill();
// }
//
// template <typename TimeT, typename DataT>
// template <typename Container>
// void Simulator<TimeT, DataT>::setActive(Container &buffers) {
//   active_buffers.clear();
//   for (auto b : buffers)
//     active_buffers.add(b);
// }

template <typename TimeT, typename DataT>
TimeT Simulator<TimeT, DataT>::currentTime() const {
  return now;
}

template <typename TimeT, typename DataT>
DataT Simulator<TimeT, DataT>::dumpRate() const {
  return data.getDumpRate(current.back().pointer);
}

template <typename TimeT, typename DataT>
DataT Simulator<TimeT, DataT>::fillRate(const int b) const {
  return data.getFillRate(b, current[b].pointer);
}

template <typename TimeT, typename DataT>
DataT Simulator<TimeT, DataT>::transferRate(const int b) const {
  return current[b].transferrate;
}

template <typename TimeT, typename DataT>
DataT Simulator<TimeT, DataT>::usage(const int b) const {
  return current[b].usage;
}

template <typename TimeT, typename DataT>
int Simulator<TimeT, DataT>::pointer(const int b) const {
  return current[b].pointer;
}

template <typename TimeT, typename DataT>
TimeT Simulator<TimeT, DataT>::getEmptyTime(const int b) const {
  if (fillRate(b) >= transferRate(b))
    return Simulator<TimeT, DataT>::never;

  // assert(usage(b) > 0);

  if (usage(b) < epsilon) {
    cout << b << endl;
    printState(true);
    exit(1);
  }

  TimeT t = now + round_and_cast<DataT, TimeT>(usage(b) /
                                               (transferRate(b) - fillRate(b)));

  if (t < now) {
    cout << "overflow!!\n";
    return Simulator<TimeT, DataT>::never;
  }
  return t;
}

template <typename TimeT, typename DataT>
template <typename Ranking>
void Simulator<TimeT, DataT>::computeRates(Ranking &priority) {
  DataT bandwidth = dumpRate();
  DataT bandwidthAlloc;

  auto batch{priority.rbegin()};
  while (batch != priority.rend()) {

    // if(print_flag)
    // 	cout << "batch " << (batch - priority.rbegin()) << " (" << batch->size()
    // << ") bandwidth=" << bandwidth << endl;

    // cout << "alloc for";
    // for (auto b : *batch)
    //   cout << " " << b;
    // cout << endl;

    if (bandwidth >= epsilon) {

      // 	// if(print_flag) {
      // 		for(auto b : *batch) {
      // 			cout << " " << b;
      // 			assert(current[b].pointer <
      // data.numBufferEvents(b));
      // 			cout << " (" << usage(b) << "/" << fillRate(b)
      // <<
      // ")";
      // 		}
      // 	cout << endl;
      // // }

      sort(batch->begin(), batch->end(), [&](const int a, const int b) {

        // // cout << "compare " << a << " and " << b << endl;
        //
        // assert(a >= 0);
        // assert(a < data.numBuffers());
        // assert(b >= 0);
        // assert(b < data.numBuffers());
        // // assert(a != b);

        auto ua{usage(a)};
        auto ub{usage(b)};
        auto fa(fillRate(a));
        auto fb(fillRate(b));

        if (ua < epsilon and ub < epsilon)
          return fa < fb;
        return ua < ub;

        // return (usage(a) < usage(b) or (usage(a) >= epsilon and usage(b) >=
        // epsilon and fillRate(a) < fillRate(b)));
      });

      // if(print_flag)
      // 	cout << "end of sort\n";

      bandwidthAlloc = (bandwidth / static_cast<DataT>(batch->size()));
    } else {
      bandwidth = 0;
      bandwidthAlloc = 0;
    }

    // if(print_flag)
    // 	cout << "bandwidth=" << bandwidth << endl;

    // cout << "standard bandwidth alloc: " << bandwidthAlloc << endl;

    for (auto b{batch->begin()}; b != batch->end(); ++b) {

      // if(print_flag)
      // 	cout << " --buffer " << *b << endl;

      if (bandwidth < epsilon) {
        current[*b].transferrate = 0;
      } else if (usage(*b) > epsilon or bandwidthAlloc < fillRate(*b)) {
        current[*b].transferrate = bandwidthAlloc;
        bandwidth -= transferRate(*b);
      } else {
        current[*b].transferrate = fillRate(*b);
        bandwidth -= transferRate(*b);
        auto d{batch->end() - b - 1};
        if (d > 0)
          bandwidthAlloc = bandwidth / static_cast<DataT>(d);
      }

      // cout << *b << ": U=" << usage(*b) << " F=" << fillRate(*b)
      //      << " T=" << transferRate(*b) << " residual bandwidth = " <<
      //      bandwidth
      //      << " fair alloc = " << bandwidthAlloc << endl;
    }

    ++batch;
  }

  // if(print_flag)
  // 	cout << "the end" << endl;
}

template <typename TimeT, typename DataT>
void Simulator<TimeT, DataT>::advanceTo(const TimeT t) {
  // for (auto b{0}; b < data.numBuffers(); ++b) {

  // cout << "avdance from " << now << " to " << t << endl;

  for (auto b{0}; b < data.numBuffers(); ++b) {
    // for (auto b : active_buffers) {

    // assert(current.size() > b);
    //
    // assert(current[b].pointer >= 0 and current[b].pointer <
    // data.numBufferEvents(b));
    //
    //     cout << b << ": " << usage(b) << " += " << "(" << fillRate(b) << " -
    //     " <<
    //     transferRate(b) <<" ) * (" << t << " - " << now << ")\n" ;

    current[b].usage +=
        (fillRate(b) - transferRate(b)) * static_cast<DataT>(t - now);

    // cout << current[b].usage << endl;

    if (bottom and current[b].usage < epsilon)
      current[b].usage = 0;

    // cout << "here\n";

    max_usage[b] = max(max_usage[b], usage(b));
    if (not excess_buffers.has(b) and max_usage[b] > capacity[b]) {
      excess_buffers.add(b);
      // cout << excess_buffers << endl;
      if (debug_flag)
        cout << "data loss for " << b << "!\n";
    }
  }
  now = t;
}

// template <typename TimeT, typename DataT>
// void Simulator<TimeT, DataT>::computeMargin() {
//   for (auto b : active_buffers) {
//     max_usage[b] = max(max_usage[b], usage(b));
//     if (not excess_buffers.has(b) and max_usage[b] > capacity[b]) {
//       excess_buffers.add(b);
// 			// cout << excess_buffers << endl;
// 			cout << "data loss for " << b << "!\n";
//     }
//   }
// }

// template <typename TimeT, typename DataT>
// DataT Simulator<TimeT, DataT>::getOverflow() {
//
//   assert(capacity.size() == data.numBuffers());
//
//   for (auto b{excess_buffers.end()}; b != excess_buffers.bend(); ++b) {
//     if (usage(*b) > capacity[*b] + epsilon) {
//       // cout << "ADD " << *b << " to " << excess_buffers << endl;
//       excess_buffers.add(*b);
//     }
//   }
//   // cout << excess_buffers << endl;
// }

template <typename TimeT, typename DataT>
template <typename Timepoint, typename Ranking>
void Simulator<TimeT, DataT>::run(Timepoint begin_event, Timepoint end_event,
                                  TimeT start_time, Ranking &priority,
                                  vector<DataT> &target,
                                  vector<DataT> &htarget) {

  if(debug_flag)
    cout << "run simu\n";

  // saved_state.clear();
  // saved_time.clear();

  // assert(active_buffers.count() == data.numBuffers());

  swap(capacity, target);
  swap(handover, htarget);

  now = start_time;

  // saved_time.push_back(now);
  // saved_state.resize(saved_state.size() + 1);
  // saved_state.back() = current;

  PRINTDEBUG

  auto cur{begin_event};

  auto dump{data.numBuffers()};

  TimeT next{cur->time};

  while (cur != end_event) {

    // cout << "advance to " << next << endl;
    assert(next > now);

    // update the usage to be at Timepoint "next" and set "now = next"

    advanceTo(next);

    // if (cur->time != now) {
    //   cout << "EMPTY EVENT! @" << now << endl;
    // } else if (cur->buffer == dump) {
    //   cout << "DOWNLINK EVENT @" << now << endl;
    // } else {
    //   cout << "BUFFER " << cur->buffer << " EVENT @" << now << endl;
    // }

    // apply all events up to now
    while (cur->time == now) {
      auto e{*cur};
      ++cur;
      ++current[e.buffer].pointer;

      if (e.buffer == dump) {
        if (dumpRate() < epsilon) {
          for (auto b{0}; b < data.numBuffers(); ++b) {
            // for (auto b : active_buffers) {
            current[b].transferrate = 0;
          }
        }
      }

      if (cur == end_event)
        break;
      else
        assert(cur->time >= now);
    }

    next = cur->time;

    // // computeMargin();
    // if (stop_early and not excess_buffers.empty())
    //   break;

    // compute the transfer rates
    if (dumpRate() >= epsilon) {
      computeRates(priority);

      // add the "empty memory" events
      for (auto b{0}; b < data.numBuffers(); ++b) {
        // for (auto b : active_buffers) {
        auto t{getEmptyTime(b)};
        if (t < next) {
          next = t;

          // cout << "next empty event for " << b << " at " << t << endl;
        }
      }
    }

    // if (debug_buffer >= 0) {
    //   cout << "buffer " << debug_buffer << ": + " << fillRate(debug_buffer)
    //        << " - " << transferRate(debug_buffer) << endl;
    // }

    // saved_time.push_back(now);
    // saved_state.resize(saved_state.size() + 1);
    // saved_state.back() = current;

    PRINTDEBUG

    // if(not excess_buffers.empty()) {
    // 	for(auto b : excess_buffers) {
    // 		cout << b << ": " << static_cast<long int>(usage(b)) << " (" <<
    // static_cast<long int>(max_usage[b]) << "/" << static_cast<long
    // int>(capacity[b]) << ")\n";
    // 	}
    // }
  }

  // cout << "hello " << handover.size() << endl;
  for (auto b{0}; b < data.numBuffers(); ++b) {
    if (usage(b) > handover[b] and not excess_buffers.has(b)) {

      if (debug_flag)
        cout << "HANDOVER FAIL ON " << b << " "
             << handover[b] / data.capacity(b) << " "
             << usage(b) / data.capacity(b) << endl;
      excess_buffers.add(b);
    }
  }

  // cout << "return " << minimum_margin << endl;

  swap(capacity, target);
  swap(handover, htarget);
  //
  // return minimum_margin;
}

template <typename TimeT, typename DataT>
template <typename Timepoint>
void Simulator<TimeT, DataT>::runWithoutTransfer(Timepoint begin_event,
                                                 Timepoint end_event,
                                                 TimeT start_time,
                                                 vector<DataT> &target) {

  // cout << "no dump simu\n";

  swap(capacity, target);
  // excess_buffers.clear();

  PRINTDEBUG

  now = start_time;

  auto cur{begin_event};

  // cout << (end_event - begin_event) << endl;

  TimeT next{cur->time};

  while (cur != end_event) {

    // cout << "advance to " << next << endl;
    assert(now == 0 or next > now);

    // update the usage to be at Timepoint "next" and set "now = next"
    advanceTo(next);

    // apply all events up to now
    while (cur->time == now) {
      auto e{*cur};
      ++cur;
      if (e.buffer < data.numBuffers())
        ++current[e.buffer].pointer;

      if (cur == end_event)
        break;
      else
        assert(cur->time >= now);
    }

    next = cur->time;

    // computeMargin();

    PRINTDEBUG
  }

  swap(capacity, target);
}

template <typename TimeT, typename DataT>
template <typename Timepoint>
void Simulator<TimeT, DataT>::runWithParallelTransfer(Timepoint begin_event,
                                                      Timepoint end_event,
                                                      TimeT start_time) {

  // cout << "no dump simu\n";

  // for (auto b{0}; b <= data.numBuffers(); ++b)
  // parallel[b].pointer = current[b].pointer;

  // swap(current, parallel);

  for (auto b{0}; b <= data.numBuffers(); ++b)
    assert(current[b].usage == 0);

  // excess_buffers.clear();

  PRINTDEBUG

  now = start_time;

  auto cur{begin_event};

  auto dump{data.numBuffers()};

  // cout << (end_event - begin_event) << endl;

  TimeT next{cur->time};

  while (cur != end_event) {

    // cout << "advance to " << next << endl;
    assert(now == 0 or next > now);

    // cout << "advance\n";

    // update the usage to be at Timepoint "next" and set "now = next"
    advanceTo(next);

    // cout << "events\n";

    // apply all events up to now
    while (cur->time == now) {
      auto e{*cur};
      ++cur;
      ++current[e.buffer].pointer;

      if (e.buffer == dump) {
        for (auto b{0}; b < data.numBuffers(); ++b)
          current[b].transferrate = dumpRate();
      }

      if (cur == end_event)
        break;
      else
        assert(cur->time >= now);
    }

    next = cur->time;

    // computeMargin();

    PRINTDEBUG
  }

  // swap(current, parallel);
}

template <typename TimeT, typename DataT>
Simulator<TimeT, DataT>::Simulator(const Instance<TimeT, DataT> &i) : data(i) {

  max_usage.resize(data.numBuffers(), 0);

  current.resize(data.numBuffers() + 1, {0, 0, 0});

  parallel.resize(data.numBuffers() + 1, {0, 0, 0});

  for (auto b{0}; b < data.numBuffers(); ++b) {
    capacity.push_back(data.capacity(b));
  }

  excess_buffers.reserve(data.numBuffers());

  // active_buffers.reserve(data.numBuffers());
  // active_buffers.fill();
}

// template <typename TimeT, typename DataT>
// void Simulator<TimeT, DataT>::initialise(vector<Event<TimeT>> &horizon,
//                                          vector<int> &initial_pointer,
//                                          vector<DataT> &initial_usage) {
template <typename TimeT, typename DataT>
void Simulator<TimeT, DataT>::initialise( // Container &buffers,
    vector<int> &initial_pointer, vector<DataT> &initial_usage) {

  // cout << "init simu\n";

  max_usage.clear();
  for (auto b{0}; b < data.numBuffers(); ++b)
    max_usage.push_back(initial_usage[b]);

  // max_usage.resize(data.numBuffers(), 0);

  // current.resize(initial_pointer.size(), {0, 0, 0});

  // assert(initial_pointer.size() == data.numBuffers() + 1);

  // one state for each buffer
  // for (auto b{0}; b < initial_pointer.size(); ++b) {
  //   current[b].pointer = initial_pointer[b];
  //   current[b].usage = initial_usage[b];
  //   current[b].transferrate = 0;
  // }
  for (auto b{0}; b < data.numBuffers(); ++b) {
    // for (auto b : active_buffers) {
    current[b].pointer = initial_pointer[b];
    current[b].usage = initial_usage[b];
    current[b].transferrate = 0;
  }

  auto dump{data.numBuffers()};
  current[dump].pointer = initial_pointer[dump];

  // cout << "1: " << excess_buffers << endl;

  // excess_buffers.reserve(data.numBuffers());
  excess_buffers.clear();

  // active_buffers.reserve(data.numBuffers());
  // active_buffers.fill();

  // cout << "INIT: " << active_buffers << endl;

  // for(auto b : buffers)
  // 	excess_buffers

  //   // excess_buffers.fill();
  //
  // cout << "2: " << excess_buffers << " - " << buffers << endl;
  //
  //   for (auto b{buffers.bbegin()}; b != buffers.bend(); ++b)
  //     excess_buffers.remove_front(*b);
  //
  // cout << "3: " << excess_buffers << endl;
  //
  //   excess_buffers.clear_front();
  //
  // cout << "4: " << excess_buffers << endl;
  //
  //   // cout << buffers << endl << excess_buffers << endl;
}

template <typename TimeT, typename DataT>
void Simulator<TimeT, DataT>::verification() const {
  DataT trate{0};
  for (auto b{0}; b < data.numBuffers(); ++b) {
    // for (auto b : active_buffers) {
    trate += transferRate(b);
    assert(current[b].pointer < data.numBufferEvents(b));
  }
  assert(trate <= dumpRate() + epsilon);
  // if(trate < dumpRate() - epsilon) {
  // 	cout << "dumprate (" << dumpRate() << ") not fully used by " <<
  // active_buffers << ":\n";
  // 	for (auto b : active_buffers) {
  // 		cout << b << " -> " << transferRate(b) << " (" << usage(b) <<
  // "|"
  // <<
  // fillRate(b) << ")"<< endl;
  // 	}
  // 	exit(1);
  // }
  // // cout << trate << " < " << (dumpRate() - epsilon) << endl;
  //
  // assert(trate >= dumpRate() - epsilon);
  assert(current[data.numBuffers()].pointer < data.numDownlinkEvents());
}

template <typename TimeT, typename DataT>
void Simulator<TimeT, DataT>::printState(const bool all) const {

  // verification();

  cout << fixed << setprecision(3);

  if (print_flag) {
    // cout << left << "@" << setw(12) << static_cast<long int>(now)
    //      << " rate=" << setw(8) << static_cast<long int>(dumpRate()) <<
    //      right;
    // for (auto b{0}; b < data.numBuffers(); ++b) {
    //   // cout << setw(12) << static_cast<int>(1000.0 * usage(b) /
    //   // data.capacity(b));
    //   cout << setw(12) << static_cast<long int>(usage(b));
    // }
    // cout << endl;

    cout << left << "@" << setw(12) << static_cast<long int>(now)
         << " rate=" << setw(8) << static_cast<long int>(dumpRate()) << right;
    for (auto b{0}; b < data.numBuffers(); ++b) {
      // cout << setw(12) << static_cast<int>(1000.0 * usage(b) /
      // data.capacity(b));
      cout << setw(12) << (usage(b) / data.capacity(b));
    }
    cout << endl;

    if (all) {
      cout << left << setw(27) << "margin:" << right << setprecision(3);
      for (auto b{0}; b < data.numBuffers(); ++b) {
        cout << setw(12)
             << (data.capacity(b) - getMaxUsage(b)) / data.capacity(b);
        // cout << setw(12) << static_cast<long int>(usage(b));
      }
      cout << endl;

      cout << left << setw(27) << "rates:" << right << setprecision(3);
      for (auto b{0}; b < data.numBuffers(); ++b) {
        cout << setw(12) << transferRate(b);
        // cout << setw(12) << static_cast<long int>(usage(b));
      }
      cout << endl;
    }

    // for (auto b{0}; b < data.numBuffers(); ++b) {
    //   // cout << setw(12) << static_cast<int>(1000.0 * usage(b) /
    //   // data.capacity(b));
    //   cout << setw(12) << static_cast<long int>(transferRate(b));
    // }
    // cout << endl;
    // cout << left << setw(31) << "fill:" << right;
    // for (auto b{0}; b < data.numBuffers(); ++b) {
    //   // cout << setw(12) << static_cast<int>(1000.0 * usage(b) /
    //   // data.capacity(b));
    //   cout << setw(12) << static_cast<long int>(fillRate(b));
    // }
    // cout << endl;

    // cout << left << setw(13) << " "
    //      << " rate=" << setw(12) << dumpRate() << right;
    // for (auto b{0}; b < data.numBuffers(); ++b) {
    //   cout << setw(6)
    //        << static_cast<int>(1000.0 * getMaxUsage(b) / data.capacity(b));
    // }
    // cout << endl;
  }

  // if (all) {
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << (data.capacity(b));
  //   }
  //   cout << endl;
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << (usage(b));
  //   }
  //   cout << endl;
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << (fillRate(b));
  //   }
  //   cout << endl;
  //   for (auto b{0}; b < data.numBuffers(); ++b) {
  //     cout << setw(12) << (transferRate(b));
  //   }
  //   cout << endl;
  // }
}

template<typename TimeT, typename DataT>
ostream& Simulator<TimeT, DataT>::display(ostream& os) const {
  os << "hello there";
  return os;
}

template<typename TimeT, typename DataT>
ostream& operator<<(ostream& os, const Simulator<TimeT, DataT>& i) {
  return i.display(os);
}

}

#endif // __Simulator_HPP
