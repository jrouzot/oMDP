#ifndef __Instance_HPP
#define __Instance_HPP

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <assert.h>

using namespace std;

namespace dataflow
{

  /*
   *  An event is caracterized by its time of occurence and the value of the fill rate
   */
  template <typename TimeT, typename DataT>
  struct EventData
  {
    TimeT time;
    DataT value;
  };

  /*
   *  Every buffer have a max and initial capacity
   */
  template <typename DataT>
  struct BufferData
  {
    DataT initial_usage;
    DataT capacity;
  };

  /**********************************************
   * Instance header
   **********************************************/

  template <typename TimeT, typename DataT>
  class Instance
  {
  public:
    Instance() {}

    /*
     *  The instance is parsed from the mtp files
     */
    Instance(const string &filename);
    ~Instance();

    /*
     *  Create a sub-instance of a given instance from the windows indexes 
     *  and the memory usage at the begin of the sub-instance
     */
    Instance(Instance ins, int w_start, int w_end, vector<DataT> memory_state);

    /*
     *  Read the file and update the class variables
     */
    void read(const string &filename);
    void readin(const string &filename, vector<BufferData<DataT>> &bf,
                vector<vector<EventData<TimeT, DataT>>> &evt,
                vector<EventData<TimeT, DataT>> &dl, const TimeT time_offset);

    /*
     *  Concatenate the current instance with another one in a file
     */
    void concatenate(const string &filename);

    size_t numBuffers() const;                 // N buffer
    size_t numDownlinks() const;               // D downlinks
    size_t numDownlinkEvents() const;          // Number of downlinks events
    size_t numBufferEvents(const int b) const; // Number of buffer events for a particular buffer

    typename vector<BufferData<DataT>>::const_iterator buffer_begin() const; // index of first buffer
    typename vector<BufferData<DataT>>::const_iterator buffer_end() const;   // end of last buffer

    typename vector<EventData<TimeT, DataT>>::const_iterator
    event_begin(const int b) const; // first event index for buffer b
    typename vector<EventData<TimeT, DataT>>::const_iterator
    event_end(const int b) const; // last event index for buffer b

    typename vector<EventData<TimeT, DataT>>::const_iterator
    downlink_begin() const; // first downlink event
    typename vector<EventData<TimeT, DataT>>::const_iterator
    downlink_end() const; // last downlink event

    typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
    event_rbegin(const int b) const; // ???
    typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
    event_rend(const int b) const; // ???

    typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
    downlink_rbegin() const; // ???
    typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
    downlink_rend() const; // ???

    DataT capacity(const int b) const; // capacity of current buffer
    DataT initial_usage(const int b) const; // initial memory of current buffer

    DataT getFillRate(const int b, const int p) const; // fill rate of current buffer for event p
    DataT getDumpRate(const int j) const;              // dump rate for ???

    TimeT getDownlinkStart(const int j) const;    // downlink start time for ???
    TimeT getDownlinkDuration(const int j) const; // downlink duration for ???
    DataT getDownlinkRate(const int j) const;     // downlink rate for ???

    void cloneBuffer(const int b, const double capacity_coeff, // clone buffer and apply a time shift
                     const double fill_coeff, const TimeT time_shift);

    TimeT horizon() const; // Total duration for the instance

    ostream &write(ostream &os) const; // Write instance in a file
    template <typename Container>
    ostream &writeSubset(ostream &os, Container &subset, const TimeT from,
                         const TimeT to) const;

    void print();                                         // Print the instance (buffers, downlinks, events)
    void setIntitialValue(int b, DataT usage);            // Update the initial memory usage of the instance 
    vector<BufferData<DataT>> getBuffer();                // Return the buffer list
    vector<EventData<TimeT, DataT>> getDownlink();        // Return the downlink events list
    vector<vector<EventData<TimeT, DataT>>> getEvents();  // Return the buffer events list

  private:
    vector<BufferData<DataT>> buffer;              // List of buffer state (capacity, start)
    vector<vector<EventData<TimeT, DataT>>> event; // List of events for each buffer
    vector<EventData<TimeT, DataT>> downlink;      // List of downlinks events
  };

  template <typename TimeT, typename DataT>
  vector<BufferData<DataT>> Instance<TimeT, DataT>::getBuffer() {
    return buffer;
  } 

  template <typename TimeT, typename DataT>
  vector<EventData<TimeT, DataT>> Instance<TimeT, DataT>::getDownlink() {
    return downlink;
  } 

  template <typename TimeT, typename DataT>
  vector<vector<EventData<TimeT, DataT>>> Instance<TimeT, DataT>::getEvents() {
    return event;
  }

  template <typename TimeT, typename DataT>
  Instance<TimeT, DataT>::Instance(const string &filename)
  {
    read(filename);
  }

  template <typename TimeT, typename DataT>
  Instance<TimeT, DataT>::~Instance()
  {
    event.clear();
  }


  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::read(const string &filename)
  {
    readin(filename, buffer, event, downlink, 0);
  }

  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::concatenate(const string &filename)
  {
    vector<BufferData<DataT>> dummy;

    readin(filename, dummy, event, downlink, horizon());
  }

  /*
   *  Read the instance file (mtp format) and update the class variables
   *
   *  The mtp file is divided into differents parts : 
   *  buffers : gives the number of buffers, their capacities and the starting memory usage
   *  downlinks : gives the downlinks informations (index, time_start, time_end, bandwidth available)
   *  events : for each buffer a list of fill rates breakpoints (start_time, fill_rate)
   */
  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::readin(
      const string &filename, vector<BufferData<DataT>> &bf,
      vector<vector<EventData<TimeT, DataT>>> &evt,
      vector<EventData<TimeT, DataT>> &dl, const TimeT time_offset)
  {
    ifstream ifs(filename); // File stream

    string flag;  // Indicate the current part

    size_t n;     // Number of buffers 
    ifs >> n;     // Read the number of buffers from the file
    ifs >> flag;  // Read the flag from the file (should be "buffer" at this point)

    auto buffer_offset{bf.size()};  // ???

    assert(flag == "buffers");      // Ensure its the buffer description part

    for (auto k{0}; k < n; ++k)     // For each buffer update buffer infos
    {
      DataT c;                      // capacity
      DataT i;                      // start memory
      ifs >> c;                     // read capacity
      ifs >> i;                     // read start memory
      bf.push_back({i, c});         // push back to the buffer list
    }

    size_t m;                       // Number of donwlinks
    ifs >> m;                       // read the number of downlinks
    ifs >> flag;                    // read the new flag (should be "downlinks")

    assert(flag == "downlinks");    // Ensure it's the downlinks part

    if (time_offset == 0)           // Insert a dummy start element at the begining
      dl.push_back({0, 0});

    for (auto k{0}; k < m; ++k)           // Read the downlinks windows
    {
      int x;                              // Downlink index
      TimeT s;                            // Start time for the downlink
      TimeT e;                            // End time for the downlink
      DataT r;                            // Bandwidth available per time unit for this window
      ifs >> x;                           // Read index
      assert(x == k);                     // Index equals downlink index in file
      ifs >> s;                           // Read start time 
      ifs >> e;                           // Read end time 
      ifs >> r;                           // Read bandwidth
      dl.push_back({s + time_offset, r}); // Push back start event
      dl.push_back({e + time_offset, 0}); // Push back end event
    }

    evt.resize(bf.size());                // Set the size of buffer events to the number of buffers
    for (auto k{0}; k < n; ++k)           // For every buffers read buffer events
    {
      size_t m;                           // Number of events for current buffer
      ifs >> m;                           // Read number of events 
      ifs >> flag;                        // Read flag (should be "events")
      assert(flag == "events");           
      ifs >> flag;                        // Read flag (should be "for")
      assert(flag == "for");
      int b;                              // buffer id
      ifs >> b;
      // assert(b == k);

      // from time 0 to the first inflexion point
      if (time_offset == 0)
        evt[buffer_offset + k].push_back({0, 0}); // Insert an event with 0 fill rate if no data at 0

      DataT v;                          // Fill_rate
      TimeT t;                          // Start time 
      for (auto j{0}; j < m; ++j)       // For each events
      {
        ifs >> t;                       // Read start time
        ifs >> v;                       // Read fill_rate
        evt[buffer_offset + k].push_back({t + time_offset, v}); // Insert buffer event
      }
    }
  }

  /*
   *  Return the number of buffer for this instance
   */
  template <typename TimeT, typename DataT>
  size_t Instance<TimeT, DataT>::numBuffers() const
  {
    return buffer.size();
  }

  /*
   *  Return the number of downlinks for this instance
   */
  template <typename TimeT, typename DataT>
  size_t Instance<TimeT, DataT>::numDownlinks() const
  {
    return (downlink.size() - 1) / 2;
  }

  /*
   *  Return the number of downlinks events for this instance
   */
  template <typename TimeT, typename DataT>
  size_t Instance<TimeT, DataT>::numDownlinkEvents() const
  {
    return downlink.size();
  }

  /*
   *  Return the number of buffer events for a particular buffer
   */
  template <typename TimeT, typename DataT>
  size_t Instance<TimeT, DataT>::numBufferEvents(const int b) const
  {
    return event[b].size();
  }

  /*
   *  Return the capacity for the given buffer
   */
  template <typename TimeT, typename DataT>
  DataT Instance<TimeT, DataT>::capacity(const int b) const
  {
    return buffer[b].capacity;
  }

  /*
   *  Return the initial memory for the given buffer
   */
  template <typename TimeT, typename DataT>
  DataT Instance<TimeT, DataT>::initial_usage(const int b) const
  {
    return buffer[b].initial_usage;
  }

  /*
   *  Return the horizon which is the time for the latest event
   */
  template <typename TimeT, typename DataT>
  TimeT Instance<TimeT, DataT>::horizon() const
  {
    TimeT maximum{0};
    for (auto i{0}; i < numBuffers(); ++i)
      if (numBufferEvents(i) > 0)
        maximum = max(maximum, event_rbegin(i)->time);

    if (numDownlinks() > 0)
      maximum = max(maximum, downlink_rbegin()->time);

    return maximum;
  }

  /*
   *  Return the dumprate for a particular downlink window
   */
  template <typename TimeT, typename DataT>
  DataT Instance<TimeT, DataT>::getDumpRate(const int j) const
  {
    return downlink[j].value;
  }

  /*
   *  Return the dump rate for a particular window 
   */
  template <typename TimeT, typename DataT>
  DataT Instance<TimeT, DataT>::getDownlinkRate(const int j) const
  {
    return downlink[2 * j + 1].value;
  }

  /*
   *  Return the duration of a given downlink window
   */
  template <typename TimeT, typename DataT>
  TimeT Instance<TimeT, DataT>::getDownlinkDuration(const int j) const
  {
    return (downlink[2 * j + 2].time - downlink[2 * j + 1].time);
  }

  /*
   *  Get the start time of the given downlink
   */
  template <typename TimeT, typename DataT>
  TimeT Instance<TimeT, DataT>::getDownlinkStart(const int j) const
  {
    return downlink[2 * j + 1].time;
  }

  /*
   *  Get the fill rate for a given buffer and a given time index 
   */
  template <typename TimeT, typename DataT>
  DataT Instance<TimeT, DataT>::getFillRate(const int b, const int p) const
  {
    return event[b][p].value;
  }

  template <typename TimeT, typename DataT>
  typename vector<BufferData<DataT>>::const_iterator Instance<TimeT, DataT>::buffer_begin() const
  {
    return buffer.begin();
  }

  template <typename TimeT, typename DataT>
  typename vector<BufferData<DataT>>::const_iterator Instance<TimeT, DataT>::buffer_end() const
  {
    return buffer.end();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_iterator
  Instance<TimeT, DataT>::downlink_begin() const
  {
    return downlink.begin();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_iterator
  Instance<TimeT, DataT>::downlink_end() const
  {
    return downlink.end();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_iterator Instance<TimeT, DataT>::event_begin(const int b) const
  {
    return event[b].begin();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_iterator Instance<TimeT, DataT>::event_end(const int b) const
  {
    return event[b].end();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
  Instance<TimeT, DataT>::downlink_rbegin() const
  {
    return downlink.rbegin();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
  Instance<TimeT, DataT>::downlink_rend() const
  {
    return downlink.rend();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
  Instance<TimeT, DataT>::event_rbegin(const int b) const
  {
    return event[b].rbegin();
  }

  template <typename TimeT, typename DataT>
  typename vector<EventData<TimeT, DataT>>::const_reverse_iterator
  Instance<TimeT, DataT>::event_rend(const int b) const
  {
    return event[b].rend();
  }

  /*
   *  Print the instance (buffers, downlinks, events)
   */
  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::print() {
    cout << numBuffers() << " buffers" << endl;
    cout << numDownlinks() << " downlinks" << endl;
    cout << "horizon : " << horizon() << endl;

    for(size_t i(0); i < numBuffers(); ++i) {
      cout << "Buffer #" << i << " : " << buffer[i].initial_usage << " " << buffer[i].capacity << endl;
      cout << "Events : ";
      for(size_t e(0); e < numBufferEvents(i); ++e) {
        cout << "(" << event[i][e].time << ", " << event[i][e].value << ") ";
      }
      cout << endl;
    }

    cout << "Downlinks events : ";
    for(size_t d(0); d < numDownlinkEvents(); ++d) {
      cout << "(" << downlink[d].time << ", " << downlink[d].value << ") ";
    }
    cout << endl;
  }

  /*
   *  Update the buffer intial memory usage
   */
  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::setIntitialValue(int b, DataT usage) {
    buffer[b].initial_usage = usage;
  }

  /*
   *  Clone a buffer with a coefficient for de capacity and fill rate and a time shift
   */
  template <typename TimeT, typename DataT>
  void Instance<TimeT, DataT>::cloneBuffer(const int b,
                                           const double capacity_coeff,
                                           const double fill_coeff,
                                           const TimeT time_shift)
  {
    buffer.push_back(buffer[b]);
    buffer.back().capacity *= capacity_coeff;
    buffer.back().initial_usage *= fill_coeff;

    event.push_back(event[b]);
    for (auto &e : event.back())
    {
      if (e.time > 0)
        e.time += time_shift;
    }
  }


  /*
   *  Write an instance to a file
   */
  template <typename TimeT, typename DataT>
  ostream &Instance<TimeT, DataT>::write(ostream &os) const
  {
    os << numBuffers() << " buffers\n";
    for (auto &b : buffer)
    {
      os << boost::lexical_cast<string>(b.capacity) << " "
         << boost::lexical_cast<string>(b.initial_usage) << endl;
    }
    os << numDownlinks() << " downlinks\n";
    int n{0};
    DataT rate;
    for (auto dl{downlink_begin() + 1}; dl != downlink_end(); ++dl)
    {
      if ((n % 2) == 0)
      {
        assert(dl->value > 0);
        os << n / 2 << " " << boost::lexical_cast<string>(dl->time);
        rate = dl->value;
      }
      else
      {
        assert(dl->value == 0);
        os << " " << boost::lexical_cast<string>(dl->time) << " "
           << boost::lexical_cast<string>(rate) << endl;
      }
      ++n;
    }
    assert(n == 2 * numDownlinks());
    for (auto b{0}; b < numBuffers(); ++b)
    {
      os << (numBufferEvents(b) - 1) << " events for " << b << endl;
      for (auto evt{event_begin(b) + 1}; evt != event_end(b); ++evt)
        os << boost::lexical_cast<string>(evt->time) << " "
           << boost::lexical_cast<string>(evt->value) << endl;
    }

    // for (auto w{downlink_begin() + 1}; w != downlink_end(); w += 2)
    //   os << " [" << w->time << "," << (w + 1)->time << "]: " << w->value;
    // os << endl;
    // for (auto b{buffer_begin()}; b != buffer_end(); ++b) {
    //   os << b->capacity << " (" << b->initial_usage << "):";
    //   for (auto e{event_begin(b - buffer_begin())};
    //        e != event_end(b - buffer_begin()); ++e) {
    //     os << " " << e->value << "@" << e->time;
    //   }
    //   os << endl;
    // }
    return os;
  }


  /*
   *  Write a sub intance with start time @from and end time @to
   */
  template <typename TimeT, typename DataT>
  template <typename Container>
  ostream &Instance<TimeT, DataT>::writeSubset(ostream &os, Container &subset,
                                               const TimeT from,
                                               const TimeT to) const
  {
    os << subset.size() << " buffers\n";
    for (auto &b : subset)
    {
      os << boost::lexical_cast<string>(buffer[b].capacity) << " "
         << boost::lexical_cast<string>(buffer[b].initial_usage) << endl;
    }

    auto beg_dl{downlink_begin() + 1};

    // cout << "search for start " << from << " (" << beg_dl->time << " " <<
    // beg_dl->value << ")";
    while (beg_dl->time < from)
    {
      assert(beg_dl->value > 0);
      beg_dl += 2;

      // cout << " " << " (" << beg_dl->time << " " << beg_dl->value << ")";
    }
    // cout << endl;

    auto end_dl{downlink_end()};

    // cout << "search for end " << to << " (" << (end_dl - 1)->time << " " <<
    // (end_dl - 1)->value << ")";
    while ((end_dl - 1)->time > to)
    {
      assert((end_dl - 1)->value == 0);
      end_dl -= 2;

      // cout << " (" << (end_dl - 1)->time << " " << (end_dl - 1)->value << ")";
    }
    // cout << endl;

    assert(beg_dl < end_dl);

    os << (end_dl - beg_dl) / 2 << " downlinks\n";
    int n{0};
    DataT rate;
    for (auto dl{beg_dl}; dl != end_dl; ++dl)
    {
      if ((n % 2) == 0)
      {
        assert(dl->value > 0);
        os << n / 2 << " " << boost::lexical_cast<string>(dl->time - from);
        rate = dl->value;
      }
      else
      {
        assert(dl->value == 0);
        os << " " << boost::lexical_cast<string>(dl->time - from) << " "
           << boost::lexical_cast<string>(rate) << endl;
      }
      ++n;
    }
    assert(n == (end_dl - beg_dl));
    auto k{0};
    for (auto &b : subset)
    {

      // cout << "buffer " << b << endl;

      auto beg_evt{event_begin(b) + 1};
      auto end_evt{event_end(b)};

      while (beg_evt < end_evt and beg_evt->time < from)
      {
        ++beg_evt;
      }

      while (end_evt > beg_evt and (end_evt - 1)->time > to)
      {
        --end_evt;
      }

      if (beg_evt < end_evt)
      {
        os << (end_evt - beg_evt) << " events for " << k << endl;
        for (auto evt{beg_evt}; evt != end_evt; ++evt)
          os << boost::lexical_cast<string>(evt->time - from) << " "
             << boost::lexical_cast<string>(evt->value) << endl;
      }
      else
      {
        os << "0 events for " << k << endl;
      }

      ++k;
    }

    assert(k == subset.size());

    // for (auto w{downlink_begin() + 1}; w != downlink_end(); w += 2)
    //   os << " [" << w->time << "," << (w + 1)->time << "]: " << w->value;
    // os << endl;
    // for (auto b{buffer_begin()}; b != buffer_end(); ++b) {
    //   os << b->capacity << " (" << b->initial_usage << "):";
    //   for (auto e{event_begin(b - buffer_begin())};
    //        e != event_end(b - buffer_begin()); ++e) {
    //     os << " " << e->value << "@" << e->time;
    //   }
    //   os << endl;
    // }
    return os;
  }

  template <typename TimeT, typename DataT>
  ostream &operator<<(ostream &os, const Instance<TimeT, DataT> &i)
  {
    return i.write(os);
  }

}

#endif // __Instance_HPP
