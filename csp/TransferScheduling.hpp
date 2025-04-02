#ifndef TRANSFER_SCHEDULING
#define TRANSFER_SCHEDULING

#include "../utils/header/InstanceSingleWindow.hpp"

namespace operations_research {

class TransferScheduling
{
    public:

        int resolution_time;                            // The resolution time of the last solve()
        string status;                                  // The status of the last solve() : UNKNOWN, SAT, UNSAT, TIMEOUT 
        vector<vector<int>> solution;                   // Solution of last solve()
        int obj_value;                                  // Objective value of last solve()
        tuple<int, int, int64_t> rmax_index;            // Store the window index, winodw index, and peak upper bound 


        TransferScheduling();


    /*
    *  Build the CSP and solve it.
    *  @mtp : the instance to solve
    *  @precision : Ratio to round the instances
    *  @time_limit : time limit to solve the instance  
    *  @obj_value_ub : The objective value treshold (we must find a max_peak value lower than @obj_value)
    *  @dc_heuristic : Use of downlink count heuristic to branch
    *  @firstfail_heuristic : Use the classical firstfail heuristic to branch
    *  @random_heuristic : Use random heuristic to branch
    *  @dense_ranking : Use of dense ranking global constraint
    *  @global_constraint : Enable priority transfer global constraint
    *  @symmetry_constraint : Enable symmetry breaking for priority transfer
    *  @luby_factor : Scale factor for the Luby sequence to perform the restarts
    *  @lds_max : Set the number of differences from the initial solution the solver can allow (if 0, lds is disabled)
    *  @delta : Enable/Disable delta constraints
    *  @initial_solution : Initial solution given by the heuristics to use as a starting point for lds 
    *  @lb : lower bound on the max peak
    *  @lns : Enables local neighborhood search. We add : peak(r_window, r_buffer) < rmax and relax the global rmax constraint. 
    *  @debug : print debug and info messages
    */
    void solve(  std::vector<dataflow::InstanceSingleWindow*> &mtp, 
                                            double precision, 
                                            int obj_value_ub,
                                            int time_limit,
                                            bool dc_heuristic,
                                            bool firstfail_heuristic,
                                            bool random_heuristic,
                                            bool dense_ranking,
                                            bool global_constraint,
                                            bool singlewindow_constraint,
                                            bool symmetry_constraint,
                                            int luby_factor,
                                            int lds_max,
                                            bool delta,
                                            vector<vector<int>> &initial_solution,
                                            int lb,
                                            bool lns,
                                            bool debug);

};

}

#endif