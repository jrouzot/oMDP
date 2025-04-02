#ifndef SIMULATION_SINGLE_WINDOW
#define SIMULATION_SINGLE_WINDOW

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

#include <cmath>

#include "Instance.hpp"
#include "InstanceSingleWindow.hpp"
#include "RoundAndCast.hpp"
#include "SparseSet.hpp"
#include "Event.hpp"
#include "State.hpp"

using namespace std;

#define TABBED_TABBED_COUT std::cout << "\t\t"

#define PRINTDEBUG printState(false);

namespace dataflow
{

    /**********************************************
     * SimulatorSingleWindow
     **********************************************/

    /*
     *  Handle the simulation algorithm
     */
    template <typename TimeT, typename DataT>
    class SimulatorSingleWindow
    {
    public:      

        /*
         *  Run simulation and track buffers that have been allocated bandwidth
         */
        template <typename Timepoint, typename Ranking>
        void runAndTrack(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, std::vector<bool> &hasBandwidth);
        
        /*
         *  Check if best priority pool relax the bandwidth
         */
        template <typename Timepoint, typename Ranking>
        bool runAndExit(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, vector<int> bufferset);

        /*
         *  Run simulation algorithm with a particular solution with buffer capacities as targets.
         *  Keeps track of the peaks
         */
        template <typename Timepoint, typename Ranking>
        void run(Timepoint begin_event, Timepoint end_event,
            TimeT start_time, Ranking &priority, vector<DataT> &target, vector<DataT> &peaks);

        /*
         *  Run simulation algorithm with a particular solution with buffer maxpeak as targets.
         *  return false if the simulation exceed the max peak
         */
        template <typename Timepoint, typename Ranking>
        void runSingleWindow(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, 
                                            vector<DataT> &target, vector<DataT> &handover);

        // /*
        //  *  Run simulation algorithm with no data transfer to find upper bounds
        //  */
        // template <typename Timepoint>
        // void runWithoutTransfer(Timepoint begin_event, Timepoint end_event,
        //                         TimeT start_time, vector<DataT> &target);

        // /*
        //  *  Run simulation algorithm with parralel data transfer to find lower bounds
        //  */
        // template <typename Timepoint>
        // void runWithParallelTransfer(Timepoint begin_event, Timepoint end_event,
        //                             TimeT start_time);

        static const DataT infinite;
        static const TimeT never;    
        static const DataT epsilon; 

        SimulatorSingleWindow(InstanceSingleWindow &i);

        void initialise(vector<int> &initial_pointers, vector<DataT> &initial_usage);

        TimeT currentTime() const;

        template <typename Ranking>
        void computeRates(Ranking &priority);

        void advanceTo(const TimeT t);

        SparseSet<> excess_buffers;

        ostream &display(ostream &os) const;

        DataT usage(const int b) const;
        int pointer(const int b) const;

        void printState(const bool all = false) const;

        DataT getMaxUsage(const int b) const;

        bool print_flag{false};
        bool debug_flag{false};
        int debug_buffer{-1};
        bool bottom{true};

        void verification() const;

    private:
        InstanceSingleWindow &data;

        vector<DataT> capacity;
        vector<DataT> peaks;
        vector<int> peaks_percent;
        vector<State<DataT>> current;
        vector<State<DataT>> baseCurrent;       // Save the initial values for current for fast initialization
        TimeT now{0};

        double dump_rate;

        vector<State<DataT>> parallel;
        vector<DataT> handover;
        vector<DataT> max_usage;

        TimeT getEmptyTime(const int b) const;
        DataT fillRate(const int b) const;
        DataT transferRate(const int b) const;

        void getOverflow();
    };

    template <typename TimeT, typename DataT>
    const DataT SimulatorSingleWindow<TimeT, DataT>::infinite = numeric_limits<DataT>::max();

    template <typename TimeT, typename DataT>
    const TimeT SimulatorSingleWindow<TimeT, DataT>::never = numeric_limits<TimeT>::max();

    template <typename TimeT, typename DataT>
    const DataT SimulatorSingleWindow<TimeT, DataT>::epsilon = static_cast<DataT>(0.001);

    template <typename TimeT, typename DataT>
    DataT SimulatorSingleWindow<TimeT, DataT>::getMaxUsage(const int b) const
    {
        return max_usage[b];
    }

    template <typename TimeT, typename DataT>
    TimeT SimulatorSingleWindow<TimeT, DataT>::currentTime() const
    {
        return now;
    }

    template <typename TimeT, typename DataT>
    DataT SimulatorSingleWindow<TimeT, DataT>::fillRate(const int b) const
    {
        return data.getFillRate(b, current[b].pointer);
    }

    template <typename TimeT, typename DataT>
    DataT SimulatorSingleWindow<TimeT, DataT>::transferRate(const int b) const
    {
        return current[b].transferrate;
    }

    template <typename TimeT, typename DataT>
    DataT SimulatorSingleWindow<TimeT, DataT>::usage(const int b) const
    {
        return current[b].usage;
    }

    template <typename TimeT, typename DataT>
    int SimulatorSingleWindow<TimeT, DataT>::pointer(const int b) const
    {
        return current[b].pointer;
    }

    /*
     *  Get the time when the buffer b is empty
     */
    template <typename TimeT, typename DataT>
    TimeT SimulatorSingleWindow<TimeT, DataT>::getEmptyTime(const int b) const
    {
        if (fillRate(b) >= transferRate(b))
            return SimulatorSingleWindow<TimeT, DataT>::never;

        if (usage(b) < epsilon)
        {
            printState(true);
            exit(1);
        }

        TimeT t = now + usage(b) / (transferRate(b) - fillRate(b));

        if (t < now)
        {
            cout << "overflow!!\n";
            return SimulatorSingleWindow<TimeT, DataT>::never;
        }
        return t;
    }

    template <typename TimeT, typename DataT>
    template <typename Ranking>
    void SimulatorSingleWindow<TimeT, DataT>::computeRates(Ranking &priority)
    {
        DataT bandwidth = dump_rate;
        DataT bandwidthAlloc;

        auto batch{priority.rbegin()};
        while (batch != priority.rend())
        {
            if (bandwidth >= epsilon)
            {
                sort(batch->begin(), batch->end(), [&](const int a, const int b)
                {
                    auto ua{usage(a)};
                    auto ub{usage(b)};
                    auto fa(fillRate(a));
                    auto fb(fillRate(b));

                    if (ua < epsilon and ub < epsilon)
                        return fa < fb;
                    return ua < ub;
                });
                bandwidthAlloc = (bandwidth / static_cast<DataT>(batch->size()));
            }
            else
            {
                bandwidth = 0;
                bandwidthAlloc = 0;
            }

            for (auto b{batch->begin()}; b != batch->end(); ++b)
            {
                if (bandwidth < epsilon)
                {
                    current[*b].transferrate = 0;
                }
                else if (usage(*b) > epsilon or bandwidthAlloc < fillRate(*b))
                {
                    current[*b].transferrate = bandwidthAlloc;
                    bandwidth -= transferRate(*b);
                }
                else
                {
                    current[*b].transferrate = fillRate(*b);
                    bandwidth -= transferRate(*b);
                    auto d{batch->end() - b - 1};
                    if (d > 0)
                        bandwidthAlloc = bandwidth / static_cast<DataT>(d);
                }
            }
            ++batch;
        }
    }

    template <typename TimeT, typename DataT>
    void SimulatorSingleWindow<TimeT, DataT>::advanceTo(const TimeT t)
    {
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            current[b].usage +=
                (fillRate(b) - transferRate(b)) * static_cast<DataT>(t - now); // Compute the usage from last time to now
            if(current[b].usage > peaks[b])
                peaks[b] = current[b].usage;

            if (bottom and current[b].usage < epsilon)
                current[b].usage = 0; // 0.00000... = 0

            max_usage[b] = max(max_usage[b], usage(b)); // Compute peak
            if (not excess_buffers.has(b) and max_usage[b] > capacity[b])
            { // Compute buffers that overflow
                // TABBED_TABBED_COUT << "buffer(" << b << ") overflows for target " << capacity[b] << " with mx usage " << max_usage[b] << " at time " << t << std::endl;
                excess_buffers.add(b);
            }
        }
        now = t;
    }


    /*
     *  Run Simulation and keep track of buffers that have been allocated bandwidth
     */
    template <typename TimeT, typename DataT>
    template <typename Timepoint, typename Ranking>
    void SimulatorSingleWindow<TimeT, DataT>::runAndTrack(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, std::vector<bool> &hasBandwidth) {
        
        now = start_time;
        auto cur{begin_event};

        auto dump{data.getNumBuffers()};

        computeRates(priority);

        for(int i{0}; i < data.getNumBuffers(); ++i) {
            if(!hasBandwidth[i] and current[i].transferrate > epsilon) {
                hasBandwidth[i] = true;
            }
        }

        advanceTo(now);

        while (cur->time == start_time)
        {
            ++cur;
        }

        TimeT next{cur->time};

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            auto t{getEmptyTime(b)};
            if (t < next)
            {
                next = t;
            }
        }

        while (cur != end_event)
        {
            advanceTo(next);

            while (cur->time == now)
            {
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;

                if (cur == end_event)
                    break;
                else
                    assert(cur->time >= now);
            }

            next = cur->time;

            // compute the transfer rates
            if (dump_rate >= epsilon and now < end_event->time)
            {
                computeRates(priority);

                for(int i{0}; i < data.getNumBuffers(); ++i) {
                    if(!hasBandwidth[i] and current[i].transferrate > epsilon) {
                        hasBandwidth[i] = true;
                    }
                }

                for (auto b{0}; b < data.getNumBuffers(); ++b)
                {
                    auto t{getEmptyTime(b)};
                    if (t < next)
                    {
                        next = t;
                    }
                }
            }
        }
    }

    /*
     *  Check if best priority pool relax the bandwidth
     */
    template <typename TimeT, typename DataT>
    template <typename Timepoint, typename Ranking>
    bool SimulatorSingleWindow<TimeT, DataT>::runAndExit(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, vector<int> bufferset) {
        
        now = start_time;
        auto cur{begin_event};

        auto dump{data.getNumBuffers()};

        computeRates(priority);
        advanceTo(now);

        while (cur->time == start_time)
        {
            ++cur;
        }

        TimeT next{cur->time};

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            auto t{getEmptyTime(b)};
            if (t < next)
            {
                next = t;
            }
        }

        while (cur != end_event)
        {
            advanceTo(next);

            while (cur->time == now)
            {
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;

                if (cur == end_event)
                    break;
                else
                    assert(cur->time >= now);
            }

            next = cur->time;

            // compute the transfer rates
            if (dump_rate >= epsilon and now < end_event->time)
            {
                computeRates(priority);

                for(auto b : bufferset) {
                    if(current[b].transferrate > epsilon) {
                        return true;
                    }
                }


                for (auto b{0}; b < data.getNumBuffers(); ++b)
                {
                    auto t{getEmptyTime(b)};
                    if (t < next)
                    {
                        next = t;
                    }
                }
            }
        }
        return false;
    }


    /*
     *  Main simulation algorithm
     */
    template <typename TimeT, typename DataT>
    template <typename Timepoint, typename Ranking>
    void SimulatorSingleWindow<TimeT, DataT>::run(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, 
                                      vector<DataT> &target, vector<DataT> &peaks_out)
    {
        swap(capacity, target);
        now = start_time;

        // Init peaks 
        for(int b(0); b < data.getNumBuffers(); ++b)
        {
            peaks[b] = current[b].usage;
        }

        PRINTDEBUG
        
        auto cur{begin_event};

        auto dump{data.getNumBuffers()};

        computeRates(priority);
        advanceTo(now);

        while (cur->time == start_time)
        {
            ++cur;
        }

        TimeT next{cur->time};

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            auto t{getEmptyTime(b)};
            if (t < next)
            {
                next = t;
            }
        }

        while (cur != end_event)
        {
            advanceTo(next);

            while (cur->time == now)
            {
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;

                if (cur == end_event)
                    break;
                else
                    assert(cur->time >= now);
            }

            next = cur->time;

            // compute the transfer rates
            if (dump_rate >= epsilon)
            {
                // // TABBED_TABBED_COUT << "compute rates" << std::endl;
                computeRates(priority);

                for (auto b{0}; b < data.getNumBuffers(); ++b)
                {
                    auto t{getEmptyTime(b)};
                    if (t < next)
                    {
                        next = t;
                    }
                }
            }
            PRINTDEBUG
        }

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            if (usage(b) > capacity[b] and not excess_buffers.has(b))
            {
                excess_buffers.add(b);
            }
            // Update peaks_out 
            peaks_out[b] = peaks[b];
        }
        swap(capacity, target);
    }


    /*
     *  Main simulation algorithm
     */
    template <typename TimeT, typename DataT>
    template <typename Timepoint, typename Ranking>
    void SimulatorSingleWindow<TimeT, DataT>::runSingleWindow(Timepoint begin_event, Timepoint end_event, TimeT start_time, Ranking &priority, 
                                      vector<DataT> &target, vector<DataT> &handover)
    {
        swap(capacity, target);
        now = start_time;

        auto cur{begin_event};

        computeRates(priority);
        advanceTo(now);

        while (cur->time == start_time)
        {
            ++cur;
        }

        TimeT next{cur->time};

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            auto t{getEmptyTime(b)};
            if (t < next)
            {
                next = t;
            }
        }

        while (cur != end_event)
        {
            advanceTo(next);

            while (cur->time == now)
            {
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;

                if (cur == end_event)
                {
                    break;
                }
                else
                    assert(cur->time >= now);
            }

            if(cur == end_event)
            {   
                break;
            }
                
            next = cur->time;

            // compute the transfer rates
            if (dump_rate >= epsilon)
            {
                computeRates(priority);

                for (auto b{0}; b < data.getNumBuffers(); ++b)
                {
                    auto t{getEmptyTime(b)};
                    if (t < next)
                    {
                        next = t;
                    }
                }
            }
        }

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            if (usage(b) > handover[b] and not excess_buffers.has(b))
            { 
                excess_buffers.add(b);
            }
        }
        swap(capacity, target);
    }

    /*
     *  Run with no transfer to find upper bounds
     */
    // template <typename TimeT, typename DataT>
    // template <typename Timepoint>
    // void SimulatorSingleWindow<TimeT, DataT>::runWithoutTransfer(Timepoint begin_event,
    //                                                  Timepoint end_event,
    //                                                  TimeT start_time,
    //                                                  vector<DataT> &target)
    // {
    //     swap(capacity, target);

    //     PRINTDEBUG

    //     now = start_time;

    //     auto cur{begin_event};

    //     while (cur->time == start_time)
    //     {
    //         ++cur;
    //     }

    //     TimeT next{cur->time};

    //     while (cur != end_event)
    //     {
    //         assert(now == 0 or next > now);
    //         advanceTo(next);

    //         while (cur->time == now)
    //         {
    //             auto e{*cur};
    //             ++cur;
    //             if (e.buffer < data.getNumBuffers())
    //                 ++current[e.buffer].pointer;

    //             if (cur == end_event)
    //                 break;
    //             else
    //                 assert(cur->time >= now);
    //         }

    //         next = cur->time;

    //         PRINTDEBUG
    //     }

    //     swap(capacity, target);
    // }


    /*
     *  Run with parralel transfer to find lower bounds
     */
    // template <typename TimeT, typename DataT>
    // template <typename Timepoint>
    // void SimulatorSingleWindow<TimeT, DataT>::runWithParallelTransfer(Timepoint begin_event,
    //                                                       Timepoint end_event,
    //                                                       TimeT start_time)
    // {
    //     for (auto b{0}; b <= data.getNumBuffers(); ++b)
    //         assert(current[b].usage == 0);

    //     PRINTDEBUG

    //     now = start_time;

    //     auto cur{begin_event};

    //     auto dump{data.getNumBuffers()};

    //     while (cur->time == start_time)
    //     {
    //         ++cur;
    //     }

    //     TimeT next{cur->time};
    //     for (auto b{0}; b < data.getNumBuffers(); ++b)
    //         current[b].transferrate = dump_rate;
            
    //     while (cur != end_event)
    //     {
    //         assert(now == 0 or next > now);
    //         advanceTo(next);

    //         while (cur->time == now)
    //         {
    //             auto e{*cur};
    //             ++cur;
    //             ++current[e.buffer].pointer;

    //             if (e.buffer == dump)
    //             {
    //                 for (auto b{0}; b < data.getNumBuffers(); ++b)
    //                     current[b].transferrate = dump_rate;
    //             }

    //             if (cur == end_event)
    //                 break;
    //             else
    //                 assert(cur->time >= now);
    //         }

    //         next = cur->time;

    //         PRINTDEBUG
    //     }
    // }

    /*
     *  Constructor
     */
    template <typename TimeT, typename DataT>
    SimulatorSingleWindow<TimeT, DataT>::SimulatorSingleWindow(InstanceSingleWindow &i) : data(i) {
        current = {};
        baseCurrent = {};
        for(short b(0); b < data.getNumBuffers() + 1; ++b)
        {
            baseCurrent.push_back({0, 0, 0});
        }

        current.resize(baseCurrent.size());

        capacity = {};
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            capacity.push_back(data.getBufferCapacity(b));
        }

        max_usage.resize(data.getNumBuffers());

        dump_rate = data.getDumpRate();
    }



    template <typename TimeT, typename DataT>
    void SimulatorSingleWindow<TimeT, DataT>::initialise(vector<int> &initial_pointer, vector<DataT> &initial_usage)
    {
        std::copy(baseCurrent.begin(), baseCurrent.end(), current.begin());

        excess_buffers.clear_back();
        excess_buffers.reserve(data.getNumBuffers());

        peaks.clear();
        peaks.reserve(data.getNumBuffers());

        peaks_percent.clear();
        peaks_percent.reserve(data.getNumBuffers());
    
        std::copy(initial_usage.begin(), initial_usage.end(), max_usage.begin()); // init max_usage to initial usage 

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            current[b].pointer = initial_pointer[b];
            current[b].usage = initial_usage[b];
            current[b].transferrate = 0;
        }

        auto dump{data.getNumBuffers()};
        current[dump].pointer = initial_pointer[dump];
    }

    // template <typename TimeT, typename DataT>
    // void SimulatorSingleWindow<TimeT, DataT>::verification() const
    // {
    //     DataT trate{0};
    //     for (auto b{0}; b < data.getNumBuffers(); ++b)
    //     {
    //         trate += transferRate(b);
    //         assert(current[b].pointer < data.numBufferEvents(b));
    //     }
    //     assert(trate <= dump_rate + epsilon);
    //     assert(current[data.getNumBuffers()].pointer < data.numDownlinkEvents());
    // }


    /*
     *  Debug tool for simulation algorithm
     */
    template <typename TimeT, typename DataT>
    void SimulatorSingleWindow<TimeT, DataT>::printState(const bool all) const
    {
        // cout << fixed << setprecision(3);

        if (print_flag)
        {
            string filename = data.getName();
            string filepath = "logs/" + filename;
            FILE * log_file = fopen(filepath.c_str(), "ae");

            // cout << left << "@" << setw(12) << now
            //     << " rate=" << setw(8) << round_and_cast_up<double, int>(dump_rate) << right;
            fprintf(log_file, "%f ", now);
            for (auto b{0}; b < data.getNumBuffers(); ++b)
            {
                fprintf(log_file, "%i ", round_and_cast_up<double, int>(100.0 * usage(b) / capacity[b]));
                // cout << setw(12) << round_and_cast_up<double, int>(100.0 * usage(b) / capacity[b]) << "%";
            }
            // cout << endl;
            fprintf(log_file, "\n");
            fclose(log_file);

            // if (all)
            // {
            //     cout << left << setw(27) << "margin:" << right << setprecision(3);
            //     for (auto b{0}; b < data.getNumBuffers(); ++b)
            //     {
            //         cout << setw(12)
            //              << (capacity[b] - getMaxUsage(b)) / capacity[b];
            //     }
            //     cout << endl;

            //     cout << left << setw(27) << "rates:" << right << setprecision(3);
            //     for (auto b{0}; b < data.getNumBuffers(); ++b)
            //     {
            //         cout << setw(12) << transferRate(b);
            //     }
            //     cout << endl;
            // }
        }
    }
}

#endif