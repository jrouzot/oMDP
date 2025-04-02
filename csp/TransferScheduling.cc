#include <ostream>
#include <stdlib.h>
#include <chrono>
#include <tuple>
#include <absl/time/time.h> 

#include "ortools/constraint_solver/constraint_solver.h"
#include "utils/header/Options.hpp"
#include "../utils/header/Instance.hpp"
#include "../utils/header/SimulationSingleWindow.hpp"
#include "../utils/header/RoundAndCast.hpp"

#include "TransferScheduling.hpp"
#include "DenseRankingConstraint.hpp"
#include "PeakConstraint.hpp"
#include "DumbPeakConstraint.hpp"
#include "DownlinkCountDecisionBuilder.hpp"
#include "SymmetryConstraint.hpp"

using namespace dataflow;
using namespace std::chrono;


/**********************************************************************
 *  CSP                                                               *
 **********************************************************************/

namespace operations_research
{
    /*
     *  Flatten the decisions variable into a 1D vector (format needed by OR-Tools solver)
     */
    std::vector<IntVar *> flatten(std::vector<std::vector<IntVar *>> matrix)
    {
        std::vector<IntVar *> flat({});
        for (auto line : matrix)
        {
            for (auto item : line)
            {
                flat.push_back(item);
            }
        }
        return flat;
    }


    /*
    *  Constructor. Init the class attributes
    */
    TransferScheduling::TransferScheduling()
    {
        resolution_time = 0;                                // The resolution time of the last solve()
        status = "UNKNOWN";                                 // The status of the last solve() : UNKNOWN, SAT, UNSAT, TIMEOUT 
        vector<vector<int>> solution = {};                  // Solution of last solve()
        int obj_value = -1;                                 // Objective value of last solve()
        tuple<int, int, int64_t> rmax_index = {};           // Store the window index, winodw index, and peak upper bound 
    }


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
    void TransferScheduling::solve(  std::vector<dataflow::InstanceSingleWindow*> &mtp, 
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
                                            bool debug)

    {
        Solver solver("TransferScheduling");                            // Instanciate a new solver
        Solver *const solver_p = &solver;                               // Pointer on the solver

        int W = mtp.size();                                             // Number of downlink windows
        int N = mtp[0]->getNumBuffers();                                // Number of buffers

        std::cout << "buffers=" << N << std::endl;
        std::cout << "windows=" << W << std::endl; 

        // std::cout << "Building the model for " << W << " windows," << N << " buffers and rmax " << obj_value_ub << std::endl;


        /*********************************
         *  PRIORITIES                   *
         *********************************/

        std::vector<std::vector<IntVar *>> priorities({});      // The priority vars for each window for each buffer

        for (size_t w(0); w < W-1; ++w)                           // For each downlink windows (W), creates the priority variables
        { 
            std::vector<IntVar *> window_priorities({});
            for (size_t b(0); b < N; ++b)                       // For each buffer, create the priority variable
            {            
                // Create the priority variable with a given ub                           
                IntVar *p = solver.MakeIntVar(0, N-1, "priority_" + to_string(w) + "_" + to_string(b)); 
                window_priorities.push_back(p);                 // Push the priority for the current window
            }
            priorities.push_back(window_priorities);            // Push the prirorities for each dwonlink windows
        }
        std::vector<IntVar *> window_priorities({});
        for (size_t b(0); b < N; ++b)                       // For each buffer, force the last priority variable to 0
        {            
            // Create the priority variable with a given ub                           
            IntVar *p = solver.MakeIntVar(0, 0, "priority_" + to_string(W-1) + "_" + to_string(b)); 
            window_priorities.push_back(p);                 // Push the priority for the current window
        }
        priorities.push_back(window_priorities);            // Push the prirorities for each dwonlink windows


        /*********************************
         *  MEMORY                       *
         *********************************/

        std::vector<std::vector<IntVar *>> memory({});                  // The memory vars for each window for each buffer (at each start of the new window)
        std::vector<double> capacities_d = mtp[0]->getCapacity();       // The original capacity for each buffer 
        std::vector<int64_t> capacities = {};                           // The rounded capacity for each buffer 

        for(short b(0); b < capacities_d.size(); ++b)
        {
            capacities.push_back(round_and_cast_down<double, int64_t>(capacities_d[b]));
        }

        std::vector<IntVar *> window_memory({});            // The memory for the first window is known
        for (size_t b(0); b < N; ++b)                       // For each buffer, create the memory variable
        { 
            IntVar *m = solver.MakeIntVar(round_and_cast_up<double, int64_t>(mtp[0]->getInitialStart(b)), 
                                          round_and_cast_up<double, int64_t>(mtp[0]->getInitialStart(b)),
                                          "memory_" + to_string(0) + "_" + to_string(b));   // Create the first memory variables
            window_memory.push_back(m);                                                     // Push the memory for the first window
        }
        memory.push_back(window_memory);                                    // Push the memory variables for the first window

        for (size_t w(1); w < W+1; ++w)
        {                                                          // For each downlinks (but the first), creates the memory variables
            std::vector<IntVar *> window_memory({});
            for (size_t b(0); b < N; ++b)
            {                                                               // For each buffer, create the memory variable
                IntVar *m = solver.MakeIntVar(0, round_and_cast_up<double, int64_t>(capacities[b] * (obj_value_ub-1) * precision),
                    "memory_" + to_string(w) + "_" + to_string(b));         // Create the memory variable
                window_memory.push_back(m);                                 // Push the memory for the current window

                if(delta)
                {
                    // Add the delta min constraints
                    solver.AddConstraint(solver.MakeGreaterOrEqual( solver.MakeDifference(m, memory[w-1][b]), 
                                                                    round_and_cast_down<double, int64_t>(mtp[w-1]->getDeltaMin(b))));
                }
            }
            memory.push_back(window_memory);                                        // Push the memory variables for all windows
        }


        /*********************************
         *  PEAKS                        *
         *********************************/

        std::vector<std::vector<IntVar *>> peaks({});                               // The peaks vars for each window
        
        IntVar * maxpeak = solver.MakeIntVar(0, obj_value_ub-1, "maxpeak");        // Constraint on the maxpeak upper bound

        solver.AddConstraint(solver.MakeGreaterOrEqual(maxpeak, lb));
        
        for (size_t w(0); w < W; ++w)                                               // For each downlinks (w), creates the peaks variables
        {       
            peaks.push_back({});
            for (size_t b(0); b < N; ++b)
            {                                                               
                IntVar *peak = solver.MakeIntVar(0, obj_value_ub-1, "peak_" + to_string(w) + "_" + to_string(b));    // Create the peak variable
                peaks[w].push_back(peak);                                           // Push the peaks variables for the current window
                // if(std::get<0>(rmax_index) == w && std::get<1>(rmax_index) == b)                 
                // {
                //     int64_t peak_ub = std::get<2>(rmax_index);
                //     solver.AddConstraint(solver.MakeLess(peak, std::get<2>(rmax_index)));
                //     // std::cout << "[CSP] Added peak(" << w << "," << b << ") < " << peak_ub << std::endl;
                // }  
                solver.AddConstraint(solver.MakeGreaterOrEqual(maxpeak, peak));     // The maxpeak is greater or equal than the other peaks
            }                                  
        }



        /***************************************
         *  CONSTRAINTS                        *
         ***************************************/

        if(lns)
        {
            for(int w(0); w < W-1; ++w)
            {
                for(int b(0); b < N; ++b)
                {
                    if(initial_solution[w][b] != -1)
                    {
                        solver.AddConstraint(solver.MakeEquality(priorities[w][b], initial_solution[w][b]));
                    }
                }
            }
        }

        if(symmetry_constraint)
        {
            for(int w(0); w < W; ++w)
            {
                solver.AddConstraint(solver.RevAlloc( 
                    new SymmetryConstraint(solver_p, priorities[w], memory[w], mtp[w])));
            }
        }

        if(dense_ranking) {
            for(int w(0); w < W; ++w)
            {
                solver.AddConstraint(solver.RevAlloc( 
                    new DenseRankingConstraint(solver_p, priorities[w], memory[w], mtp[w])));
            }
        }

        // if(global_constraint)
        //     std::cout << "Priority Transfer enabled" << std::endl;

        // if(singlewindow_constraint)
        //     std::cout << "Single Window enabled" << std::endl;

        if(global_constraint || singlewindow_constraint)
        {
            for (size_t w(0); w < W; ++w)               // Add the peak constraints to compute the peaks and next memory state from the priorities domain
            {
                solver.AddConstraint(solver.RevAlloc(
                    new PeakConstraint(solver_p, memory[w], memory[w + 1], priorities[w], peaks[w], maxpeak, mtp[w], global_constraint, singlewindow_constraint, precision)));
            }
        }
        else
        {
            for (size_t w(0); w < W; ++w)               // Add the peak dummy constraints to compute the peaks and next memory state from the priorities values
            {
                solver.AddConstraint(solver.RevAlloc(
                        new DumbPeakConstraint(solver_p, memory[w], memory[w + 1], priorities[w], peaks[w], mtp[w], precision)));
            }
        }


        /***************************************
         *  LDS                                *
         ***************************************/

        if (lds_max >= 0)
        {

            vector<IntVar *> differences = {}; // This vector represent the the differences from the initial solution for each prio var

            for (size_t w(0); w < W; ++w)
            {
                for (size_t b(0); b < N; ++b)
                {
                    IntVar *diff = solver.MakeIntVar(0, 1, "diff_" + to_string(w) + "_" + to_string(b));
                    solver.AddConstraint(solver.MakeIsDifferentCstCt(priorities[w][b], initial_solution[w][b], diff));
                    differences.push_back(diff);
                }
            }

            // Add the number of max differences constraint
            solver.AddConstraint(solver.MakeSumEquality(differences, lds_max));
        }

        /***************************************
         *  Monitors                           *
         ***************************************/

        std::vector<SearchMonitor *> monitors;

        // Search log
        SearchMonitor * log = nullptr;

        // if(debug)
        // {
        //     log = solver.MakeSearchLog(1000);
        //     monitors.push_back(log);
        // } 

        // log = solver.MakeSearchLog(1000);
        // monitors.push_back(log);
        
        // Restarts with Luby sequence
        SearchMonitor* restart = solver.MakeLubyRestart(luby_factor);
        monitors.push_back(restart);

        // Search limit 
        const absl::Duration absl_time_limit = absl::Milliseconds(time_limit);
        SearchLimit *const limit = solver.MakeTimeLimit(absl_time_limit); 
        monitors.push_back(limit);


        /*********************************
         *  SOLVE                        *
         *********************************/

        // Flatten the decision variable to branch on
        std::vector<IntVar *> decision_variables = flatten(priorities);

        DecisionBuilder * db = nullptr; 
        
        if(dc_heuristic) {
            // std::cout << "Decision builder : Downlink Count" << std::endl;
            db = new DownlinkCountDecisionBuilder(mtp, priorities, memory, maxpeak);     // Downlink Count decision builder
        }
        else if(firstfail_heuristic) {
            // std::cout << "Decision builder : Min dom" << std::endl;
            db = solver.MakePhase(
                decision_variables, Solver::CHOOSE_MIN_SIZE, Solver::ASSIGN_MIN_VALUE);
        } else if(random_heuristic) {
            // std::cout << "Decision builder : Random" << std::endl;
            db = solver.MakePhase( // Dummy decision builder
                decision_variables, Solver::CHOOSE_RANDOM, Solver::ASSIGN_RANDOM_VALUE);
        }

        auto start = high_resolution_clock::now();              // Save time before solve starts
        solver.NewSearch(db, monitors);                         // Solve the model


        /*********************************
         *  SOLUTION DISPLAY             *
         *********************************/

        resolution_time = 0;
        status = "UNSAT";
        solution = {};

        if (solver.NextSolution())
        {            
            status = "SAT";    

            // Update objective value
            obj_value = -1;
            for (size_t w(0); w < peaks.size(); ++w)
            {
                for (size_t b(0); b < peaks[w].size(); ++b)
                {
                    if(peaks[w][b]->Value() > obj_value)
                        obj_value = peaks[w][b]->Value();
                }
            }

            // Update solution
            for (size_t w(0); w < priorities.size(); ++w)
            {
                solution.push_back({});
                for (size_t b(0); b < priorities[w].size(); ++b)
                {
                    solution[w].push_back(priorities[w][b]->Value());
                }
            } 

            for (size_t w(0); w < priorities.size(); ++w)
            {
                for (size_t b(0); b < N; ++b)
                {   
                    if(peaks[w][b]->Value() == maxpeak->Min())
                    {
                        rmax_index = std::make_tuple(w, b, obj_value);
                    }
                }
            }
        }

        solver.EndSearch();                                         // We reached the time limit or the optimal solution

        std::cout << "branches=" << solver.branches() << std::endl;
        std::cout << "failures=" << solver.failures() << std::endl;

        auto stop = high_resolution_clock::now();                   // Stop the chrono
        auto duration = duration_cast<microseconds>(stop - start);  // Compute duration
        resolution_time = round(duration.count() / 1000);           // Update resolution time

        if (resolution_time >= time_limit)
            status = "TIMEOUT";

        if(debug)
            std::cout << status << std::endl;
    } 
}