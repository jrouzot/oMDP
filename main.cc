#include <iostream>
#include <limits>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <ctime>

#include "utils/header/Instance.hpp"
#include "utils/header/InstanceSingleWindow.hpp"
#include "utils/header/AlgorithmSingleWindow.hpp"
#include "utils/header/Options.hpp"
#include "utils/header/PlaySolution.hpp"
#include "utils/header/FilterInstances.hpp"


using namespace std;
using namespace dataflow;
using namespace mysolver;

void printStat(string method, tuple<vector<vector<int>>, double, int> res)
{
    vector<vector<int>> sol;
    double obj;
    int time;
    std::tie(sol, obj, time) = res;

    std::cout << "real_obj=" << obj << std::endl;
    std::cout << "real_time=" << time << std::endl;
}

// void printBounds(string method, tuple<vector<vector<int>>, double, double, int> res)
// {
//     vector<vector<int>> sol;
//     double lb;
//     double ub;
//     int time;
//     std::tie(sol, lb, ub, time) = res;

//     std::cout << method << ":" << std::endl;
//     std::cout << "---lb=" << lb << std::endl;
//     std::cout << "---ub=" << ub << std::endl;
//     std::cout << "---time=" << time << std::endl;
// }


int main(int argc, char *argv[])
{
    Options opt = parse(argc, argv);                        // Get the instances from the command line
    auto mtp = Instance<double, double>(opt.instance_file);
    vector<InstanceSingleWindow*> instances;                // List of single window instances

    // std::cout << "Running filter for instance file: " << opt.instance_file << std::endl;
    // std::cout << "Instance N=" << mtp.numBuffers() << ", W=" << mtp.numDownlinks() << std::endl;

    // FilterInstances filter = FilterInstances(opt, mtp);
    // if(filter.FilterLb()) {
    //     return 1;           // Filter the file
    // }
    // return 0;               // Do not filter the file

    // Create the single window instances
    for (short w(0); w <= mtp.numDownlinks(); ++w)
    {
        instances.push_back(new InstanceSingleWindow(mtp, opt.instance_file, w, opt.factor));
    }

    // Instanciate the solver
    mysolver::Solver solver = mysolver::Solver(mtp, instances, opt);

    std::cout << "global=" << opt.global << std::endl;
    std::cout << "denseranking=" << opt.dense_ranking << std::endl;
    std::cout << "symmetry=" << opt.symmetry << std::endl;
    std::cout << "single_window=" << opt.singlewindow << std::endl;
    
    // // Solve with IteratedLeveling heuristic
    // std::cout << "method=IteratedLeveling" << std::endl;
    // printStat("IteratedLeveling", solver.solveIteratedLeveling());

    // // Solve with RepairDescent heuristic
    // std::cout << "method=RepairDescent" << std::endl;
    // printStat("RepairDescent", solver.solveRepairDescent());

    // // Find bounds and instance pre-processing
    // printBounds("CSPBounds", solver.breakAndSolveCSPBounds(obj, opt.time_limit));

    // Solve with our CSP
    if(opt.randomwalk) {
        printStat("RandomWalk", solver.solveRandomWalk());
    }
    if(opt.iterativeleveling) {
        printStat("IterativeLeveling", solver.solveIterativeLeveling());
    }
    if(opt.repairdescent) {
        printStat("RepairDescent", solver.solveRepairDescent());
    }
    if(opt.lns) {
        printStat("LNS", solver.solveLNS(opt.time_limit, opt.time_limit*0.1));
    }
    if(!opt.randomwalk && !opt.iterativeleveling && !opt.repairdescent && !opt.lns) {
        printStat("CSP", solver.solveCSP(opt.time_limit));
    }



    // if(opt.lns) {
    //     printStat("LNS", solver.solveLNS(opt.time_limit, opt.time_limit/10));
    // } else {
    //     printStat("CSP", solver.solveCSP(opt.time_limit));
    // }

    // // Solve with LDS
    // printStat("LDS", solver.solveLDS(opt.lds_max, opt.time_limit));

    // // Solve with CSP and instance pre-processing
    // std::cout << "method=CSP+Preprocessing" << std::endl;
    // printStat("CSP+Preprocessing", solver.breakAndSolveCSP(opt.time_limit));

    // // Solve with CSP and instance pre-processing
    // std::cout << "method=LDS+Preprocessing" << std::endl;
    // printStat("LDS+Preprocessing", solver.breakAndSolveLDS(opt.lds_max, opt.time_limit));

    // // Solve with LNS
    // std::cout << "method=LNS" << std::endl;
    // printStat("LNS", solver.solveLNS(opt.time_limit, opt.time_limit/10));

    // // Solve with LNS and instance pre-processing
    // std::cout << "method=LNS+Preprocessing" << std::endl;
    // printStat("LNS+Preprocessing", solver.solveLNS(opt.time_limit, opt.time_limit/1000));
}