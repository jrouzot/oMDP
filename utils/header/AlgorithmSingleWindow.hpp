#ifndef ALGO_SINGLE_WINDOW
#define ALGO_SINGLE_WINDOW

#include <iostream>
#include <limits>
#include <vector>
#include <tuple>

#include "ortools/constraint_solver/constraint_solver.h"

#include "SimulationSingleWindow.hpp"
#include "InstanceSingleWindow.hpp"
#include "PrintVector.hpp"
#include "Event.hpp"
#include "State.hpp"

#define TABBED_COUT std::cout << "\t"

using namespace std;
using namespace operations_research;

namespace dataflow
{
    template <typename TimeT, typename DataT>
    class AlgorithmSingleWindow
    {
    public:
        AlgorithmSingleWindow(InstanceSingleWindow &i);

        void setTargetPeak(const double peak);
        double getCurrentUsage(int buffer);

        void playDownlinkAndTrack(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<bool> &bufferset);
        bool playDownlinkAndExit(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<int> bufferset);
        void playDownlink(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<DataT> &peaks, const bool verbose);

        void initialise();

        bool singleWindowWithVar(vector<DataT> &memory_state, vector<DataT> &handover, vector<IntVar *> &priorities, DataT &peak, const bool verbose);

        void printAlgorithmState();

    private:
        InstanceSingleWindow &data;

        SimulatorSingleWindow<TimeT, DataT> simulator;

        vector<double> hdelta;
        vector<DataT> maxusage;
        vector<DataT> usage;

        vector<int> pointer;
        vector<Event<TimeT>> events;

        vector<DataT> target;
        vector<DataT> nolimit;

        // We store the initial values for fast initialization 
        vector<DataT> baseUsage;
        vector<DataT> baseMaxusage;
        vector<DataT> basePointer;
        vector<DataT> baseNolimit;

        vector<vector<int>> priority;
        vector<vector<int>> best_priority;
        vector<vector<int>> baseBestPriority;

        bool debug_flag{false};
        bool fail_flag{false};

        double computeMinimumMargin() const;
        double computeMaxPeak() const;
        void getMaxUsage();

        double marginOf(const int b) const;
        double peakOf(const int b) const;
        double percentFull(const int b) const;
        DataT delta(const int b) const;
    };

    /*
     *  Constructor
     */
    template <typename TimeT, typename DataT>
    AlgorithmSingleWindow<TimeT, DataT>::AlgorithmSingleWindow(InstanceSingleWindow &i)
        : data(i), simulator(data) {
            // Capacity target
            target = {};
            events = {};

            // Get all buffer events
            for (auto b{0}; b < data.getNumBuffers(); ++b)
            {
                target.push_back(static_cast<double>(data.getBufferCapacity(b)));
                for (auto e{data.getBufferBegin(b)}; e != data.getBufferEnd(b); ++e)
                {
                    events.push_back({e->time, b});
                }
            }

            // Sort all events by time 
            sort(events.begin(), events.end());
            
            events.insert(events.begin(), {static_cast<double>(data.getTimeOffset()), data.getNumBuffers()});       // Add start downlink event
            events.push_back({static_cast<double>(data.getHorizon()), data.getNumBuffers()});                       // Add end downlink event

            // Actual variables
            pointer = {};
            usage = {};
            maxusage = {};
            nolimit = {};

            // Initialization variables
            basePointer = {};
            baseUsage = {};
            baseMaxusage = {};
            baseNolimit = {};

            for(short b(0); b < data.getNumBuffers(); ++b)
            {
                basePointer.push_back(0);
                baseUsage.push_back(data.getMemoryStartBuffer(b));
                baseMaxusage.push_back(data.getMemoryStartBuffer(b));
                baseNolimit.push_back(SimulatorSingleWindow<TimeT, DataT>::infinite);
            }

            basePointer.push_back(0);

            pointer.resize(basePointer.size());
            usage.resize(baseUsage.size());
            maxusage.resize(baseMaxusage.size());
            nolimit.resize(baseNolimit.size());

            best_priority = {};
            best_priority.resize(data.getNumBuffers());
            baseBestPriority = {};
            for(int i(0); i < data.getNumBuffers(); ++i)
            {
                baseBestPriority.push_back({});
                baseBestPriority.reserve(data.getNumBuffers());
            }
        }

    /*
     *  Target margin
     */
    template <typename TimeT, typename DataT>
    void AlgorithmSingleWindow<TimeT, DataT>::setTargetPeak(const double peak)
    {
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            assert(nolimit[b] > data.getBufferCapacity(b));
            target[b] = data.getBufferCapacity(b) * peak;
        }
    }

    /*
     *  Margin of buffer b
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::marginOf(const int b) const
    {
        return static_cast<double>(data.getBufferCapacity(b) - maxusage[b]) /
               static_cast<double>(data.getBufferCapacity(b));
    }

    /*
     *  Peak of buffer b
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::peakOf(const int b) const
    {
        return static_cast<double>(maxusage[b]) /
               static_cast<double>(data.getBufferCapacity(b));
    }


    /*
     *  Percentage of memory usage
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::percentFull(const int b) const
    {
        return static_cast<double>(maxusage[b]) / static_cast<double>(data.getBufferCapacity(b));
    }

    /*
     *  Delta margin
     */
    template <typename TimeT, typename DataT>
    DataT AlgorithmSingleWindow<TimeT, DataT>::delta(const int b) const
    {
        return maxusage[b] - usage[b];
    }

    /*
     *  Run simulation and keep track of buffer that have been allocated bandwidth
     */
    template <typename TimeT, typename DataT>
    void AlgorithmSingleWindow<TimeT, DataT>::playDownlinkAndTrack(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<bool> &hasBandwidth) {
        // Update usage from the memory at the begin of the window 
        for(short b(0); b < data.getNumBuffers(); ++b)
        {
            usage[b] += memory_state[b];
        }

        simulator.initialise(pointer, usage);
        simulator.runAndTrack(events.begin(), events.end(), data.getTimeOffset(), prio, hasBandwidth);
    }


    /*
     *  Run simulation and exit if priority pool 1 get bandwidth
     */
    template <typename TimeT, typename DataT>
    bool AlgorithmSingleWindow<TimeT, DataT>::playDownlinkAndExit(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<int> bufferset) {
        // Update usage from the memory at the begin of the window 
        for(short b(0); b < data.getNumBuffers(); ++b)
        {
            usage[b] += memory_state[b];
        }

        simulator.initialise(pointer, usage);

        bool exit = simulator.runAndExit(events.begin(), events.end(), data.getTimeOffset(), prio, bufferset);
        return exit;
    }

    /*
     *  Run simulation algorithm for the downlink
     *  for a given set of priorities.
     */
    template <typename TimeT, typename DataT>
    void AlgorithmSingleWindow<TimeT, DataT>::playDownlink(vector<vector<int>> &prio, vector<DataT> &memory_state, vector<DataT> &peaks,
                                                 const bool verbose)
    {
        // Update usage from the memory at the begin of the window 
        for(short b(0); b < data.getNumBuffers(); ++b)
        {
            usage[b] += memory_state[b];
        }

        simulator.initialise(pointer, usage);

        if (verbose)
            simulator.print_flag = true;

        simulator.run(events.begin(), events.end(), data.getTimeOffset(), prio, target, peaks);
        
        if (verbose)
            simulator.print_flag = false;

        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            memory_state[b] = simulator.usage(b);
        }

        getMaxUsage();
        fail_flag = not simulator.excess_buffers.empty();
    }


    /*
     *  Get the current memory usage of buffer b
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::getCurrentUsage(int b)
    {
        return usage[b];
    }

    /*
     *  Initialisation
     */
    template <typename TimeT, typename DataT>
    void AlgorithmSingleWindow<TimeT, DataT>::initialise()
    {
        std::copy(basePointer.begin(), basePointer.end(), pointer.begin());
        std::copy(baseUsage.begin(), baseUsage.end(), usage.begin());
        std::copy(baseMaxusage.begin(), baseMaxusage.end(), maxusage.begin());
        std::copy(baseNolimit.begin(), baseNolimit.end(), nolimit.begin());

        // for (auto b{0}; b < data.getNumBuffers(); ++b)
        // {
        //     double p = static_cast<double>(usage[b]) /
        //                     static_cast<double>(data.getBufferCapacity(b));
        // }
    }

    /*
     *  Minimum margin
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::computeMinimumMargin() const
    {
        double minimum_margin{1};
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            assert(maxusage[b] == simulator.getMaxUsage(b));

            double margin{marginOf(b)};

            if (margin < minimum_margin)
                minimum_margin = margin;
        }
        return minimum_margin;
    }

    /*
     *  Max peak
     */
    template <typename TimeT, typename DataT>
    double AlgorithmSingleWindow<TimeT, DataT>::computeMaxPeak() const
    {
        double max_peak{0};
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            assert(maxusage[b] == simulator.getMaxUsage(b));

            double peak{peakOf(b)};

            if (peak > max_peak)
                max_peak = peak;
        }
        return max_peak;
    }

    /*
     *  Max usage
     */
    template <typename TimeT, typename DataT>
    void AlgorithmSingleWindow<TimeT, DataT>::getMaxUsage()
    {
        for (auto b{0}; b < data.getNumBuffers(); ++b)
        {
            maxusage[b] = simulator.getMaxUsage(b);
        }
    }


    /*
     *  Single window with domains
     *  return <bool, bool> fail or maj of domains 
     */
    template <typename TimeT, typename DataT>
    bool AlgorithmSingleWindow<TimeT, DataT>::singleWindowWithVar(  vector<DataT> &memory_state, 
                                                                    vector<DataT> &handover,
                                                                    vector<IntVar *> &priorities, 
                                                                    DataT &peak,                                                                   
                                                                    const bool verbose)
    {
        // Update usage from the memory at the begin of the window 
        for(short b(0); b < data.getNumBuffers(); ++b)
        {
            usage[b] += memory_state[b];
        }

        bool proc = false;          // Can we filter some priorities domains 

        setTargetPeak(peak);

        // Find the highest upper bound among the IntVar domains
        int worst_prio = 0;

        best_priority = baseBestPriority;

        for (int b(0); b < data.getNumBuffers(); ++b) {
            best_priority[data.getNumBuffers() - 1 - priorities[b]->Max()].push_back(b);
        }  

        bool overflow = true;

        while (overflow)
        {
            overflow = false;

            // std::cout << "Running Simulation with prio:" << std::endl;
            // print_2Dvec(best_priority);

            simulator.initialise(pointer, usage);
            simulator.runSingleWindow(events.begin(), events.end(), data.getTimeOffset(), best_priority, target, handover);

            for (auto b{0}; b < data.getNumBuffers(); ++b)
            {
                if (simulator.excess_buffers.has(b))
                {
                    // std::cout << "Buffer " << b << " overflows" << std::endl;

                    int current_index = data.getNumBuffers() - 1 - priorities[b]->Max();

                    // Check if current buffer can have a better priority
                    if(priorities[b]->Bound())
                    {
                        return true;    // Overflow
                    }

                    // Otherwise we raise the priority of current buffer
                    int offset = 1;
                    while(!(priorities[b]->Contains(priorities[b]->Max() - offset)))
                        ++offset;

                    // Find the position of the value in best_priorities
                    auto it = std::find(best_priority[current_index].begin(), best_priority[current_index].end(), b);

                    // Check if the value was found
                    if (it != best_priority[current_index].end()) {
                        best_priority[current_index + offset].insert(best_priority[current_index + offset].end(), b);
                        best_priority[current_index].erase(it);
                    }

                    priorities[b]->SetMax(priorities[b]->Max() - offset);
                    
                    // Indicates that a buffer overflows
                    overflow = true;
                }
            }
        }
        return false;   // No overflow
    }
}

#endif