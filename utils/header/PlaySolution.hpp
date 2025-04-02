#ifndef PLAY_SOLUTION
#define PLAY_SOLUTION

#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "AlgorithmSingleWindow.hpp"
#include "InstanceSingleWindow.hpp"
#include "RoundAndCast.hpp"

using namespace dataflow;

/*
 *  Return the solution in a batch format
 */
std::vector<std::vector<int>> GetSolutionBatch(int N, std::vector<int> sol)
{
    std::vector<std::vector<int>> priorities({});

    for (int i(0); i < N; ++i)
    {
        priorities.emplace_back(std::vector<int>{});
    }

    for (int i(0); i < N; ++i)
    {
        priorities[sol[i]].emplace_back(i);
    }
    reverse(priorities.begin(), priorities.end());
    return priorities;
}

void playSolutionSingleWindow(std::vector<std::vector<int>> solution, std::vector<InstanceSingleWindow*> instances)
{
    const int N = solution[0].size();
    vector<double> mem_state(N, 0);

    for(int b(0); b < N; ++b)        // Get initial memory state
    {
        mem_state[b] += instances[0]->getInitialStart(b);
    }

    std::vector<double> peaks(N, 0);
    
    // Write solution
    const string filename = instances[0]->getName();
    const string sol_filepath = "solutions/" + filename;
    FILE * sol_log_file = fopen(sol_filepath.c_str(), "we");

    for(auto item : solution)
    {
        for(auto i : item)
        {
            fprintf(sol_log_file, "%i ", i);
        }
        fprintf(sol_log_file, "\n");
    }
    fclose(sol_log_file);

    // Write logs details
    const string filepath = "logs/" + filename;
    FILE * log_file_reset = fopen(filepath.c_str(), "we");

    // Print initial state at t = 0
    fprintf(log_file_reset, "%f ", 0.0);
    for(short b(0); b < instances[0]->getNumBuffers(); ++b)
    {
        int64_t state = round_and_cast_up<double, int64_t>(mem_state[b]*100/instances[0]->getBufferCapacity(b));
        fprintf(log_file_reset, "%ld ", state);
    }
    fprintf(log_file_reset, "\n");
    fclose(log_file_reset);

    // Print the memory state for every window using simulation algorithm 
    for(int i(0); i < instances.size() - 1; ++i)
    {
        unique_ptr<AlgorithmSingleWindow<double, double>> algo = make_unique<AlgorithmSingleWindow<double, double>>(*instances[i]);
        algo->initialise();
        vector<vector<int>> p = GetSolutionBatch(N, solution[i]);
        algo->playDownlink(p, mem_state, peaks, true);
        algo.reset();
    }

    FILE * log_file = fopen(filepath.c_str(), "ae");
    fprintf(log_file, "%f ", instances[instances.size() - 1]->getHorizon());
    for(int b(0); b < N; ++b)
    {
        const int64_t state = round_and_cast_up<double, int64_t>(mem_state[b] + instances[instances.size() - 1]->getMemoryStart()[b])*100/instances[instances.size() - 1]->getBufferCapacity(b);
        fprintf(log_file, "%ld ", state);
    }
    fclose(log_file);   
} 

#endif