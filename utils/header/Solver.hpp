#ifndef SOVLER_H
#define SOLVER_H

#include <vector>
#include <tuple>
#include <filesystem>

#include "Instance.hpp"
#include "InstanceSingleWindow.hpp"
#include "Options.hpp"
#include "../../csp/TransferScheduling.hpp"
#include "BreakInstance.hpp"
#include "SplitVector.hpp"
#include "Algorithm.hpp"
#include "PlaySolution.hpp"
#include "BreakInstance.hpp"

using namespace std;
using namespace dataflow;
using namespace operations_research;

const double PRECISION_SOLVE = 0.001;       // How much do we approximate the peaks
const double PRECISION_INSTANCE = 1.0;      // How do we approximate the instance

// Compute current time
int now()
{
    return chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
}

// Convert the solution format to match a vector<vector<int>> respecting a dense ranking
vector<vector<int>> convert(vector<vector<vector<int>>> sol, int N)
{
    vector<vector<int>> solution({});

    // Remove empty vectors (to have a Dense Ranking soution)
    for (size_t i(0); i < sol.size(); ++i)
    {
        sol[i].erase(
            remove_if(
                sol[i].begin(),
                sol[i].end(),
                [](const vector<int> &vec)
                { return vec.empty(); }),
            sol[i].end());
    }

    // Convert solution format
    for (size_t i(0); i < sol.size(); ++i)
    {
        reverse(sol[i].begin(), sol[i].end());
        vector<int> line(N, 0);
        solution.push_back(line);
        for (size_t j(0); j < sol[i].size(); ++j)
        {
            for (size_t k(0); k < sol[i][j].size(); ++k)
            {
                solution[i][sol[i][j][k]] = j;
            }
        }
    }

    return solution;
}


// Print solution Algorithm 
void printSolutionAlgo(vector<vector<vector<int>>> solution_algo)
{
    for(auto window : solution_algo)
    {
        cout << "{" << endl;
        for(auto rank : window)
        {
            cout << "\t{ ";
            for(auto b : rank)
            {
                cout << b << " ";
            }
            cout << "}" << endl;
        }
        cout << "}" << endl;
    }
}

namespace mysolver
{

class Solver
{
private:
    vector<vector<int>> solution_;                          // Heuristic algorithm to find an initial solution
    Options opt_;                                           // Options for the heuristic
    Instance<double, double> mtp_;                          // Instance
    vector<InstanceSingleWindow*> instance_;                // The instance, divided into single windows       
    int peak_lb_;                                           // Rounded peak found by relaxing the bandwidth share constraints
    double real_lb_;                                        // Peak found by relaxing the bandwidth share constraints
    int N_;                                                 // Number of buffers
    int W_;                                                 // Number of downlinks

public:
    Solver( Instance<double, double> mtp, 
            vector<InstanceSingleWindow*> instance, 
            Options opt);

    vector<vector<vector<int>>> dummySolution3D();                                                                  // Build a dummy solution with equals priorities (3D format)
    vector<vector<int>> dummySolution2D();                                                                          // Build a dummy solution with equals priorities (2D format)
    double playSolution(vector<vector<int>> sol);                                                                   // Play a solution and get objective   
    double getRealLb();
    tuple<vector<vector<int>>, double, int> solveRandomWalk();                                                      // Solve with RandomWalk
    tuple<vector<vector<int>>, double, int> solveIterativeLeveling();                                               // Solve with iterative  leveling heuritic and get the solution, max peak and resolution time
    tuple<vector<vector<int>>, double, int> solveRepairDescent();                                                   // Solve with repair descent heuritic and get the solution, max peak and resolution time
    tuple<vector<vector<int>>, double, int> solveCSP(int time_limit);                                               // Solve with our constraint model and get the solution, max peak and resolution time
    tuple<vector<vector<int>>, double, int> solveLNS(int time_limit, int time_restart);                             // Solve with our LNS method and get the solution, max peak and resolution time                                                              // Try to prove unfaisability 
};

Solver::Solver( Instance<double, double> mtp, 
                vector<InstanceSingleWindow*> instance, 
                Options opt)
                : mtp_(mtp), instance_(instance), opt_(opt), N_(mtp_.numBuffers()), W_(mtp_.numDownlinks()+1)
{
    // Get the peak lower_bound
    Algorithm<double, double> algostart(mtp_, opt_);  // Algorithm class to run optimistic simulation
    algostart.initialise();
    real_lb_ = 1.0 - algostart.relaxBandwidth();
    peak_lb_ = round_and_cast_up<double, int>((real_lb_) / PRECISION_SOLVE);
}

// Build dumb solution (all buffer at best priorities for all windows)
vector<vector<vector<int>>> Solver::dummySolution3D()
{
    vector<vector<vector<int>>> dummy_solution({});
    for (size_t j(0); j < W_; ++j)
    {
        dummy_solution.push_back({});
        vector<int> line({});
        for (size_t i(0); i < N_; ++i)
        {
            line.push_back(i);
        }
        dummy_solution[j].push_back({line});
    }
    return dummy_solution;
}

vector<vector<int>> Solver::dummySolution2D()
{
    vector<vector<int>> dummy_solution({});
    for (size_t w(0); w < W_; ++w)
    {
        dummy_solution.push_back({});
        for (size_t b(0); b < N_; ++b)
        {
            dummy_solution[w].push_back(0);
        }
    }
    return dummy_solution;
}


double Solver::playSolution(vector<vector<int>> sol)
{
    vector<vector<vector<int>>> solution_algo({});

    for(int w(0); w < W_; ++w)
    {
        solution_algo.push_back({});
        for(int b(0); b < N_; ++b)
        {
            solution_algo[w].push_back({});
        }
        for(int b(0); b < N_; ++b)
        {
            solution_algo[w][sol[w][b]].push_back(b);
        }
    }

    // Remove empty vectors (to have a Dense Ranking soution)
    for (size_t i(0); i < solution_algo.size(); ++i)
    {
        reverse(solution_algo[i].begin(), solution_algo[i].end()); // solution ranks must be reversed for simulation algo
        solution_algo[i].erase(
            remove_if(
                solution_algo[i].begin(),
                solution_algo[i].end(),
                [](const vector<int> &vec)
                { return vec.empty(); }),
            solution_algo[i].end());
    }

    Algorithm<double, double> algo(mtp_, opt_);
    algo.initialise();
    return (1 - algo.playSolution(solution_algo, false)) * 100;
}

double Solver::getRealLb()
{
    return real_lb_;
}

tuple<vector<vector<int>>, double, int> Solver::solveRandomWalk()
{
    seed(1);       
    Algorithm<double, double> algo(mtp_, opt_);
    algo.initialise();
    int start = now();
    vector<vector<vector<int>>> dummy_solution = dummySolution3D();
    int peak_ub = round_and_cast_up<double, int>((1 - algo.playSolution(dummy_solution, false)) / PRECISION_SOLVE + 1);
    std::cout << "obj=" << peak_ub << std::endl;
    std::cout << "time=" << now() - start << std::endl;
    vector<vector<vector<int>>> solution = algo.randomWalk()->priority;
    double peak = (1 - algo.playSolution(solution, false)) * 100;
    int resolution_time = now() - start;
    return make_tuple(convert(solution, N_), peak, resolution_time);
}

tuple<vector<vector<int>>, double, int> Solver::solveIterativeLeveling()
{
    seed(1);       
    Algorithm<double, double> algo(mtp_, opt_);
    algo.initialise();
    int start = now();
    vector<vector<vector<int>>> dummy_solution = dummySolution3D();
    int peak_ub = round_and_cast_up<double, int>((1 - algo.playSolution(dummy_solution, false)) / PRECISION_SOLVE + 1);
    std::cout << "obj=" << peak_ub << std::endl;
    std::cout << "time=" << now() - start << std::endl;
    vector<vector<vector<int>>> solution = algo.iterativeLeveling()->priority;
    double peak = (1 - algo.playSolution(solution, false)) * 100;
    int resolution_time = now() - start;
    return make_tuple(convert(solution, N_), peak, resolution_time);
}


tuple<vector<vector<int>>, double, int> Solver::solveRepairDescent()
{
    seed(1);                
    Algorithm<double, double> algo(mtp_, opt_);
    algo.initialise();
    int start = now();
    vector<vector<vector<int>>> dummy_solution = dummySolution3D();
    int peak_ub = round_and_cast_up<double, int>((1 - algo.playSolution(dummy_solution, false)) / PRECISION_SOLVE + 1);
    std::cout << "obj=" << peak_ub << std::endl;
    std::cout << "time=" << now() - start << std::endl;
    vector<vector<vector<int>>> solution = algo.lnsDescent()->priority;
    // std::cout << "W_=" << W_ << std::endl;
    // std::cout << "solution size=" << solution.size() << std::endl;
    // std::cout << algo.playSolution(solution, false) << std::endl;
    double peak = (1 - algo.playSolution(solution, false)) * 100;
    // std::cout << "peak=" << peak << std::endl;
    int resolution_time = now() - start;

    return make_tuple(convert(solution, N_), peak, resolution_time);
}


tuple<vector<vector<int>>, double, int> Solver::solveCSP(int time_limit)
{
    TransferScheduling csp;                             // Create a new solver
    string status = "UNKOWN";                           // Keep track of the solve status    
    int resolution_time = 0;                            // Keep track of the resolution time
    vector<vector<int>> solution = dummySolution2D();   // Keep track of the solution

    // Run dummy solution to get an upper bound on the peak
    Algorithm<double, double> algo(mtp_, opt_);
    algo.initialise();
    vector<vector<vector<int>>> solution_warmstart;
    if(opt_.warm_start) {
        solution_warmstart = algo.lnsDescent()->priority;
        solution = convert(solution_warmstart, N_);
    } else {
        solution_warmstart = dummySolution3D();
    }
    int peak_ub = round_and_cast_up<double, int>((1 - algo.playSolution(solution_warmstart, false)) / PRECISION_SOLVE + 1);

    std::cout << "peak-ub=" << peak_ub << std::endl;

    while (status != "UNSAT" && status != "TIMEOUT")
    {
        csp.solve(  instance_,
                    PRECISION_SOLVE, 
                    peak_ub,
                    opt_.time_limit - resolution_time,
                    opt_.dc_heuristic,
                    opt_.firstfail_heuristic,
                    opt_.random_heuristic,
                    opt_.dense_ranking,
                    opt_.global,
                    opt_.singlewindow,
                    opt_.symmetry,
                    opt_.luby_factor,
                    -1,
                    opt_.delta,
                    solution,
                    peak_lb_,
                    false,
                    opt_.debug);

        status = csp.status;                                // Update status
        resolution_time += csp.resolution_time;             // Update resolution time

        if(status == "SAT")
        {
            solution = csp.solution;                        // Update best solution
            peak_ub = csp.obj_value;                        // Update best peak
            std::cout << "obj=" << peak_ub << std::endl;
            std::cout << "time=" << resolution_time << std::endl;
        }                      
    }

    std::cout << "original_status=" << status << std::endl;
    std::cout << "original_lb=" << peak_lb_ << std::endl;
    std::cout << "obj=" << peak_ub << std::endl;
    std::cout << "time=" << resolution_time << std::endl;

    if(peak_ub <= peak_lb_+1 || status == "UNSAT") {        // We found the optimal solution
        status = "OPTIMAL";
    }
    std::cout << "status=" << status << std::endl;

    double real_peak = playSolution(solution);
    return make_tuple(solution, real_peak, resolution_time);
}


/*
 *  Modify the solution according to which window/buffer overflows.
 *  @solution : the last solution found by the solver.
 *  @depth : how much we relax the solution (how many window do we backtrack)
 */
void relaxSolution(vector<vector<int>> &solution, int depth, tuple<int, int, int64_t> r_index, int N, int W)
{
    if(get<0>(r_index) - depth < 0)
        return;
    for(int w(get<0>(r_index) - depth); w <= get<0>(r_index) + depth / 2; ++w)
    {
        if(w < W) {
            int prio = solution[w][get<1>(r_index)];
            for(int b(0); b < N; ++b)
            {
                if(solution[w][b] <= prio || b == get<1>(r_index))
                {
                    solution[w][b] = -1;
                }
            }
        }
    }
}

tuple<vector<vector<int>>, double, int> Solver::solveLNS(int time_limit, int time_restart)
{
    std::cout << "Starting LNS..." << std::endl; 

    TransferScheduling csp;                                     // Create a new solver
    string status = "SAT";                                      // Keep track of the solve status 
    bool lns = true;                                            // Enable large neighborhood search
    int depth = 0;                                              // Each time we find UNSAT, we relax more variables (ie. @depth windows)
    int resolution_time = 0;                                    // Init solve time

    // Run repair descent to get an initial solution
    auto res = solveRepairDescent();
    vector<vector<int>> solution;
    double obj;
    int time;
    std::tie(solution, obj, time) = res;
    vector<vector<int>> tmp = solution;
    int peak_ub = round_and_cast_up<double, int>(obj*10 + 1);

    resolution_time += time;

    // first solve to get the max peak index
    csp.solve(  instance_,
                PRECISION_SOLVE, 
                peak_ub,         // Make a valid objective with the initial solution
                time_limit - resolution_time,
                opt_.dc_heuristic,
                opt_.firstfail_heuristic,
                opt_.random_heuristic,
                opt_.dense_ranking,
                opt_.global,
                opt_.singlewindow,
                opt_.symmetry,
                opt_.luby_factor,
                -1,
                opt_.delta,
                tmp,               
                peak_lb_,
                lns,
                opt_.debug);

    while(peak_ub > peak_lb_ && resolution_time < time_limit)  // Solve until we reached the time limit
    {
        // print_2Dvec(tmp);
        csp.solve(  instance_,
                    PRECISION_SOLVE, 
                    peak_ub,
                    time_restart,
                    opt_.dc_heuristic,
                    opt_.firstfail_heuristic,
                    opt_.random_heuristic,
                    opt_.dense_ranking,
                    opt_.global,
                    opt_.singlewindow,
                    opt_.symmetry,
                    opt_.luby_factor,
                    -1,
                    opt_.delta,
                    tmp,               
                    peak_lb_,
                    lns,
                    opt_.debug);

        status = csp.status;                                // Update status
        resolution_time += csp.resolution_time;             // Update resolution time
        if(status == "SAT")
        {
            depth = 1;
            solution = csp.solution;                                    // Update global solution
            tmp = solution;

            peak_ub = csp.obj_value;                                    // Update best peak only if its better
            std::cout << "obj=" << peak_ub  << std::endl;               // Print current obj
            std::cout << "time=" << resolution_time << std::endl;       // Print current time
            relaxSolution(tmp, depth, csp.rmax_index, N_, W_);
        }
        if(status == "TIMEOUT" || status == "UNSAT")
        {
            // Relax some variables of the solution (set to -1, then relaxed during the model creation in csp)
            relaxSolution(tmp, depth, csp.rmax_index, N_, W_);  
            if(depth > get<0>(csp.rmax_index))
            {
                double real_peak = playSolution(solution);
                return make_tuple(solution, real_peak, resolution_time);
            } 
            ++depth;                                        // Increase depth
        }
    }

    std::cout << "obj=" << peak_ub << std::endl;         // Print current obj
    std::cout << "time=" << resolution_time << std::endl;       // Print current time
    
    if(peak_ub <= peak_lb_+1 || status == "UNSAT") {        // We found the optimal solution
        status = "OPTIMAL";
        std::cout << "lb=" << peak_ub << std::endl;
    } else {
        std::cout << "lb=" << peak_lb_ << std::endl;
    }
    std::cout << "status=" << status << std::endl;

    double real_peak = playSolution(solution);
    return make_tuple(solution, real_peak, resolution_time);
}
}

#endif