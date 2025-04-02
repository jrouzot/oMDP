#ifndef BREAK_INSTANCE
#define BREAK_INSTANCE

#include "AlgorithmSingleWindow.hpp"
#include "InstanceSingleWindow.hpp"
#include "RoundAndCast.hpp"

#include <cmath>

using namespace dataflow;


/*
 *  Find the breakpoints for the initial instances.
 *  A breakpoint can be found at the end of a window if all buffer are empty even in the worst case scenario.
 *  
 *  Return a vector of window indexes. 
 */
auto GetBreakpoints(std::vector<InstanceSingleWindow *> &instances) -> std::vector<int>
{
    int N = instances[0]->getNumBuffers();
    int W = instances.size() - 1;               // The last window is not relevant

    // Store all breakpoints for the instance
    std::vector<int> breakpoints;

    // Keep track of the total memory accumulated
    double mem_state = 0;  

    // Init mem state with initial memory usage
    for(int b(0); b < N; ++b)
    {
        mem_state += instances[0]->getInitialStart(b);
    }

    for(int w(0); w < W; ++w)
    {
        for(int b(0); b < N; ++b)
        {
            mem_state += instances[w]->getTotalFill(b);
        }

        double dump = instances[w]->getTotalDump();
        mem_state = max(mem_state - dump, 0.0);

        // Add breakpoint iff all memory are empty at the end of the window even in the worst case scenario
        if(mem_state == 0.0)
            breakpoints.emplace_back(w);
    }
    breakpoints.emplace_back(W);
    return breakpoints;
}

#endif