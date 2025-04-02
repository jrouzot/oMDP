#ifndef __INSTANCE_SINGLE_WINDOW
#define __INSTANCE_SINGLE_WINDOW

#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "Instance.hpp"

using namespace std;

namespace dataflow
{

class InstanceSingleWindow
{
    public:
        InstanceSingleWindow();                                 // Default constructor
        InstanceSingleWindow(Instance<double, double> ins, string filename, int window_index, int factor);                
                                                                // We build the InstanceSingleWindow from the Instance class 

        void buildInstanceEnd(Instance<double, double> ins, int factor);    // Get the memory usage for each buffer after the last visibility window
        void print();                                           // Print the instance
        void roundThis(double coef);                            // Round the instance values with @coef

        vector<vector<EventData<double, double>>> getEvents();
        vector<EventData<double, double>>::const_iterator getBufferBegin(int b);
        vector<EventData<double, double>>::const_iterator getBufferEnd(int b); 
        vector<double> getMemoryStart();
        double getMemoryStartBuffer(int b);
        int64_t getMemoryStartIntBuffer(int b);
        double getInitialStart(int b);
        void setInitalStart(int b, int value);
        vector<double> getCapacity();
        double getBufferCapacity(int b);
        int64_t getBufferCapacityInt(int b);
        double getFillRate(const int b, const int p) const; 
        double getDumpRate();
        int64_t getDumpRateInt();
        int getNumBuffers();
        int getId();
        double getTotalFill(int i);
        int64_t getTotalFillInt(int i);
        double getTotalDump();
        int64_t getTotalDumpInt();
        double getHorizon();
        void setMemoryStart(int index, double value);
        string getName();
        double getTimeOffset();
        double getDeltaMin(int index);

    private:
        /*
         *  List of events for each buffers.
         *  For each buffer, the first event must have a non-zero value and the last event must have a zero value. 
         *  The number of events is pair and the events alternate between zero values and non-zero values.
         *  The time for each event must be strictly ordered (no duplicates).
         */ 
        vector<vector<EventData<double, double>>> buffer_event; 
        vector<double> delta_min;

        int num_buffers;                        // Number of buffers
        vector<double> memory_start;            // Memory state for each buffer at the start of the downlink window 
        vector<int64_t> memory_start_int;       // Memory start rounded to int
        vector<double> initial_start;           // Memory state at time 0 on the initial instance 
        vector<double> capacity;                // Capacity for each buffer
        vector<int64_t> capacity_int;           // Capacity rounded to int64_t for each buffer
        double horizon;                         // End time for the downlink window
        double dump_rate;                       // Dump rate during the downlink window
        int64_t dump_rate_int;                  // Dump rate rounded up to integer
        int id;                                 // Window id
        vector<double> total_fill;              // Total memory during the window
        vector<int64_t> total_fill_int;         // Total memory int during the window
        string name;                            // Name of the instance
        double time_offset;                     // Start time of the instance

        void computeDeltaMin();                 // Compute minimum delta for the memory usage for each buffer at each window 
};

}

#endif