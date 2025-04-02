#include <stdlib.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <sys/stat.h>

#include "../header/InstanceSingleWindow.hpp"
#include "../header/Instance.hpp"
#include "../header/SimulationSingleWindow.hpp"
#include "../header/PrintVector.hpp"
#include "../header/RoundAndCast.hpp"

using namespace std;
using namespace dataflow;

namespace dataflow
{

/*
 *  Default constructor
 */
InstanceSingleWindow::InstanceSingleWindow() {
    buffer_event = {};
    memory_start = {};
    dump_rate = 0.0;
    num_buffers = 0;
    name = "";
}

/*
 *  Convert an Instance for a given window into an InstanceSingleWindow 
 */
InstanceSingleWindow::InstanceSingleWindow( Instance<double, double> ins, string filename, int window_index, int factor) : 
                                            num_buffers(ins.numBuffers()) {

    // Find the last occurrence of '/' in the string
    size_t pos = filename.find_last_of('/');

    // Extract the filename from the string
    if (pos != string::npos) {
        name = filename.substr(pos + 1);
    } else {
        // If '/' is not found, the string contains only the filename
        name = "";
    }

    vector<EventData<double, double>> downlink = ins.getDownlink();         // Get the downlink events
    vector<vector<EventData<double, double>>> event = ins.getEvents();      // Get the buffer events

    initial_start = {};
    for(size_t b(0); b < num_buffers; ++b) {
        initial_start.push_back(ins.initial_usage(b));
    }

    if(2*window_index+1 == downlink.size()) {
        buildInstanceEnd(ins, factor);
        return;
    }

    buffer_event = {};
    memory_start = {};
    
    num_buffers = ins.numBuffers();                                         // Get the number of buffers

    double end_last_downlink = downlink[2*window_index].time;               // Get time points 
    double start_downlink = downlink[2*window_index+1].time;
    double end_downlink = downlink[2*window_index+2].time;
    time_offset = start_downlink;
    horizon = end_downlink;

    dump_rate = downlink[2*window_index+1].value;                   // Get the dump rate for the current window   
    dump_rate_int = round_and_cast_up<double, int64_t>(dump_rate); 

    // Fix the instances generated
    for(size_t b(0); b < num_buffers; ++b)
    {
        if(!(event[b].back().value == 0))
        {
            event[b].push_back({ins.horizon(), 0});
        }
    }          

    // Init memory start and capacity for each buffers
    for(size_t b(0); b < num_buffers; ++b) {
        memory_start.push_back(0.0);
        memory_start_int.push_back(0);
        capacity.push_back(ins.getBuffer()[b].capacity * factor);
        capacity_int.push_back(round_and_cast_up<double, int64_t>(ins.getBuffer()[b].capacity * factor));
    } 

    // Init buffer events
    for(size_t b(0); b < num_buffers; ++b) {
        buffer_event.push_back({});
    } 

    // Get memory start for each buffer and get the buffer events
    for(size_t b(0); b < num_buffers; ++b) 
    {
        // If no event add dummy start event and return
        if(!event[b].size())
        {
            buffer_event[b].push_back({start_downlink, 0});
            break;
        }   

        auto e = event[b].begin(); 

        // Get memory start
        while(e->time <= end_last_downlink)
        {
            ++e;
            if(e == event[b].end())
                break;
        }    

        if(e->time > start_downlink)
        {
            memory_start[b] += (e-1)->value*(start_downlink - end_last_downlink);
        } 
        else 
        {
            while(e->time < start_downlink && e != event[b].end())
            {
                memory_start[b] += (e-1)->value*(e->time - max(end_last_downlink, (e-1)->time));
                ++e;
            }   
            memory_start[b] += (e-1)->value*(start_downlink - max(end_last_downlink, (e-1)->time));
        }
        
        // Get events
        if(e->time > start_downlink && (e-1)->time < start_downlink)
        {
            buffer_event[b].push_back({start_downlink, (e-1)->value});
        } 
        else if(e->time >= end_downlink)
        {
            buffer_event[b].push_back({start_downlink, (e-1)->value});
        } 

        while(e != event[b].end() && e->time < end_downlink)
        {
            buffer_event[b].push_back({e->time, e->value});
            ++e;
        }   
        if(!buffer_event[b].size())
            buffer_event[b].push_back({start_downlink, 0});
        if(buffer_event[b].back().value != 0)
            buffer_event[b].push_back({end_downlink, 0});
    }

    for(size_t b(0); b < num_buffers; ++b) {
        memory_start_int[b] = round_and_cast_up<double, int64_t>(memory_start[b]);
    }

    time_offset = downlink[2*window_index+1].time;

    // // Re-scale the time horizon
    // for(size_t b(0); b < num_buffers; ++b) {
    //     for(size_t e(0); e < buffer_event[b].size(); ++e) {
    //         buffer_event[b][e] = {buffer_event[b][e].time - time_offset, buffer_event[b][e].value};
    //     }
    // }

    // Add start events
    for(size_t b(0); b < num_buffers; ++b) {
        if(!buffer_event[b].size() || (*buffer_event[b].begin()).time != time_offset) {
            buffer_event[b].insert(buffer_event[b].begin(), {time_offset, 0});
        } else {

        }
    }


    // Set the total fill
    for(short int i(0); i < num_buffers; ++i)
    {
        total_fill.push_back(memory_start[i]);
        for(short int e(0); e < buffer_event[i].size(); ++e)
        {
            total_fill[i] += buffer_event[i][e].value*(buffer_event[i][e+1].time - buffer_event[i][e].time);
        }
    }   

    // Set the total fill int
    for(short int i(0); i < num_buffers; ++i)
    {
        total_fill_int.push_back(round_and_cast_up<double, int64_t>(memory_start[i]));
        for(short int e(0); e < buffer_event[i].size(); ++e)
        {
            total_fill_int[i] += round_and_cast_up<double, int64_t>(buffer_event[i][e].value*(buffer_event[i][e+1].time - buffer_event[i][e].time));
        }
    }   

    computeDeltaMin();
}


/*
 *  Print the instance. 
 */
void InstanceSingleWindow::print() {
    cout << num_buffers << " buffers" << endl;

    for(size_t b(0); b < num_buffers; ++b) {
      cout << "Buffer #" << b << " (start : " << memory_start[b] << ", capacity : " << capacity[b] << ")" << endl;
      if(buffer_event.size()) {
        cout << "Events : { ";
        for(size_t e(0); e < buffer_event[b].size(); ++e) {
            cout << "(" << buffer_event[b][e].time << ", " << buffer_event[b][e].value << ") ";
        }
        cout << "}" << endl;
      }
    }
    cout << "Downlinks dump_rate : " << dump_rate << endl;
    cout << "Horizon : " << horizon << endl;
}


/*
 *  Round the instance. 
 */
void InstanceSingleWindow::roundThis(double coef) {

    for(size_t b(0); b < num_buffers; ++b) {
        memory_start[b] = round(memory_start[b]*coef);
        capacity[b] = round(capacity[b]*coef);
        if(buffer_event.size()) {
            for(size_t e(0); e < buffer_event[b].size(); ++e) {
                buffer_event[b][e].value = round(buffer_event[b][e].value*coef);
            }
        }
    }
    dump_rate = round(dump_rate*coef);
}


/*
 *  Get the memory usage after the last non visibility window
 */
void InstanceSingleWindow::buildInstanceEnd(Instance<double, double> ins, int factor) {

    buffer_event = {};
    memory_start = {};
    dump_rate = 0;
    
    num_buffers = ins.numBuffers(); 

    vector<EventData<double, double>> downlink = ins.getDownlink();         // Get the downlink events
    vector<vector<EventData<double, double>>> event = ins.getEvents();      // Get the buffer events

    double end_last_window = downlink[downlink.size()-1].time;              // Get the last window time

    horizon = ins.horizon();

    dump_rate = 0.0;
    dump_rate_int = 0;

    // Init capacity for each buffers
    for(size_t b(0); b < num_buffers; ++b) {
        capacity.push_back(ins.getBuffer()[b].capacity * factor);
        capacity_int.push_back(round_and_cast_up<double, int64_t>(ins.getBuffer()[b].capacity * factor));
    } 

    // Init buffer events
    for(size_t b(0); b < num_buffers; ++b) {
        buffer_event.push_back({});
    } 


    // Fix the instances generated
    for(size_t b(0); b < num_buffers; ++b)
    {
        if(!(event[b].back().value == 0))
        {
            event[b].push_back({ins.horizon(), 0});
        }
    }   

    for(size_t b(0); b < num_buffers; ++b) {
        memory_start.push_back(0);                                  // Init memory variables
        memory_start_int.push_back(0);
        for(int e(0); e < event[b].size()-1; ++e) {

            // Skip the events that are not in the current window
            while(event[b][e+1].time < end_last_window) {
                ++e;
                if(e >= event[b].size()-1)
                    break;
            }  

            if(e >= event[b].size()-1)
                break;

            // The event that start before the end of the last window and that ends after
            if(event[b][e].value > 0 && event[b][e].time < end_last_window && event[b][e+1].time > end_last_window) {
                memory_start[b] += event[b][e].value * (event[b][e+1].time - end_last_window);
            }

            // The events that starts end ends after the last window 
            if(event[b][e].value > 0 && event[b][e].time >= end_last_window && event[b][e+1].time > end_last_window) {
                memory_start[b] += event[b][e].value * (event[b][e+1].time - event[b][e].time);
            }
        }   
    }

    for(size_t b(0); b < num_buffers; ++b) {
        memory_start_int[b] = round_and_cast_up<double, int64_t>(memory_start[b]);
    }


    // Add a dummy event for each buffer 
    for(int b(0); b < num_buffers; ++b)
    {
        buffer_event[b].push_back({horizon, 0});
    }

    // Set the total fill and print
    for(short int i(0); i < num_buffers; ++i)
    {
        total_fill.push_back(memory_start[i]);
    }  

    // Set the total fill and print
    for(short int i(0); i < num_buffers; ++i)
    {
        total_fill_int.push_back(round_and_cast_up<double, int64_t>(memory_start[i]));
    }  

    computeDeltaMin();
}

/*
 *  Compute minimum delta for the memory usage for each buffer at each window
 *
 *  As the maximum dump rate and the fill rate are know we can compute the minimum memory raise during a given window (worst case scenario) 
 */
void InstanceSingleWindow::computeDeltaMin()
{
    delta_min = {};
    for (size_t w(0); w < num_buffers; ++w)
    {
        double delta = total_fill[w] - getTotalDump();
        delta_min.emplace_back(delta);
    }
}

/*
 *  GETTERS
 */

vector<vector<EventData<double, double>>> InstanceSingleWindow::getEvents() {
    return buffer_event;
} 

vector<EventData<double, double>>::const_iterator InstanceSingleWindow::getBufferBegin(int b) {
    return buffer_event[b].begin();
} 

vector<EventData<double, double>>::const_iterator InstanceSingleWindow::getBufferEnd(int b) {
    return buffer_event[b].end();
} 

double InstanceSingleWindow::getHorizon() {
    return horizon;
}

vector<double> InstanceSingleWindow::getMemoryStart() {
    return memory_start;
}

double InstanceSingleWindow::getInitialStart(int b) {
    return initial_start[b];
}

void InstanceSingleWindow::setInitalStart(int b, int value) {
    initial_start[b] = value;
}

double InstanceSingleWindow::getMemoryStartBuffer(int b) {
    return memory_start[b];
}

int64_t InstanceSingleWindow::getMemoryStartIntBuffer(int b) {
    return memory_start_int[b];
}


vector<double> InstanceSingleWindow::getCapacity() {
    return capacity;
}

double InstanceSingleWindow::getBufferCapacity(int b) {
    return capacity[b];
}

int64_t InstanceSingleWindow::getBufferCapacityInt(int b) {
    return capacity_int[b];
}

double InstanceSingleWindow::getFillRate(const int b, const int p) const {
    if(b < num_buffers) {
        return buffer_event[b][p].value;
    } else {
        return 0;
    }
    
}

double InstanceSingleWindow::getDumpRate() {
    return dump_rate;
}

int64_t InstanceSingleWindow::getDumpRateInt() {
    return dump_rate_int;
}

int InstanceSingleWindow::getNumBuffers() {
    return num_buffers;
}
        
int InstanceSingleWindow::getId() {
    return id;
}

void InstanceSingleWindow::setMemoryStart(int index, double value) {
    memory_start[index] = value;
}

/*
 *  Returns the total fill during the window for a given buffer
 */
double InstanceSingleWindow::getTotalFill(int i)
{
    return total_fill[i];
}

/*
 *  Returns the total fill int during the window for a given buffer
 */
int64_t InstanceSingleWindow::getTotalFillInt(int i)
{
    return total_fill_int[i];
}

/*
 *  Returns the total dump during the window
 */
double InstanceSingleWindow::getTotalDump()
{
    return dump_rate*(horizon-time_offset);
}

int64_t InstanceSingleWindow::getTotalDumpInt()
{
    return dump_rate_int*(horizon-time_offset);
}

/*
 *  Returns the name of the instance
 */
string InstanceSingleWindow::getName()
{
    return name;
}

/*
 *  Returns the time offset of the instance
 */
double InstanceSingleWindow::getTimeOffset()
{
    return time_offset;
}

double InstanceSingleWindow::getDeltaMin(int index)
{
    return delta_min[index];
}

}